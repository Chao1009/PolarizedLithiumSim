"""Sampler-level checks: distributions, rates, reproducibility, weights."""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen import bookkeeping as bk  # noqa: E402
from polligen.sample import InclusiveSampler  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.polarized import toy_b1, toy_delta_gluon  # noqa: E402


UNPOL = bk.SpinCategory("unpol", 1.0, (1 / 3, 1 / 3, 1 / 3))


@pytest.fixture(scope="module")
def sampler():
    kern = InclusiveKernel(
        beams.LI6, b1_func=toy_b1,
        delta_func=lambda x, q2, f1: toy_delta_gluon(x, q2, f1, scale=1e-2))
    config = beams.default_configs("6Li")[1]
    return InclusiveSampler(kern, config, fom.Scenario(), nx=24, nq2=18,
                            x_range=(1e-3, 0.5), q2_range=(1.5, 100.0))


def test_unpolarized_cell_distribution(sampler):
    rng = np.random.default_rng(1)
    n = 150_000
    ev = sampler.sample_category(UNPOL, n=n, rng=rng)
    counts = np.bincount(ev["cell"], minlength=sampler.xsec_flat.size)
    # population-averaged rates = unpolarized cell cross sections
    expected = n * sampler.xsec_flat / sampler.xsec_flat.sum()
    use = expected > 50
    z = (counts[use] - expected[use]) / np.sqrt(expected[use])
    assert np.abs(z).max() < 5.0
    assert z.std() == pytest.approx(1.0, abs=0.25)


def test_event_kinematics_consistent(sampler):
    ev = sampler.sample_category(UNPOL, n=20_000,
                                 rng=np.random.default_rng(2))
    np.testing.assert_allclose(ev["y"], ev["q2"] / (sampler.s * ev["x"]),
                               rtol=1e-12)
    assert np.all((ev["phi"] >= 0) & (ev["phi"] < 2 * np.pi))
    # events stay inside their cells
    assert np.all(ev["x"] >= np.exp(sampler.logx_lo[ev["cell"]]))
    assert np.all(ev["x"] <= np.exp(sampler.logx_hi[ev["cell"]]))


def test_tensor_rate_asymmetry_matches_analytic(sampler):
    # pure-state accepted cross sections reproduce the sigma-weighted Azz
    from polli_fastsim.asymmetries import azz
    t = sampler.tables
    y = sampler.q2_cells / (sampler.s * sampler.x_cells)
    truth = np.average(azz(t["b1"], t["f1"], t["f2"], sampler.x_cells, y),
                       weights=sampler.xsec_flat)
    pure = {m: sampler.sigma_tot_pb(UNPOL, m=m) for m in (1.0, 0.0, -1.0)}
    measured = ((pure[1.0] + pure[-1.0] - 2 * pure[0.0])
                / (pure[1.0] + pure[-1.0] + pure[0.0]))
    assert measured == pytest.approx(truth, rel=1e-9)


def test_poisson_rate_normalization(sampler):
    lumi_pb = 1e-4 * sampler.scenario.lumi_fb_per_nucleon * 1e3
    mu = sampler.expected_events(UNPOL, lumi_pb)
    rng = np.random.default_rng(3)
    draws = [sampler.sample_category(UNPOL, lumi_pb=lumi_pb,
                                     rng=rng)["x"].size
             for _ in range(60)]
    assert np.mean(draws) == pytest.approx(mu, abs=4 * np.sqrt(mu / 60))


def test_bunch_rng_reproducible(sampler):
    ev1 = sampler.sample_category(UNPOL, n=5000, rng=bk.bunch_rng(7, 1, 42))
    ev2 = sampler.sample_category(UNPOL, n=5000, rng=bk.bunch_rng(7, 1, 42))
    ev3 = sampler.sample_category(UNPOL, n=5000, rng=bk.bunch_rng(7, 1, 43))
    np.testing.assert_array_equal(ev1["x"], ev2["x"])
    np.testing.assert_array_equal(ev1["phi"], ev2["phi"])
    assert not np.array_equal(ev1["x"], ev3["x"])


def test_vector_rate_sign_and_labels(sampler):
    """Deterministic guard on the population <-> m pairing: category rates
    must equal the kernel evaluated at EXPLICITLY constructed spin states,
    including the sign of the vector-sector rate shift.  (A mutation that
    reverses the m-ordering convention passes every statistical closure
    test -- the A_par truths are smaller than the MC tolerances -- but
    fails here exactly.)"""
    from polli_fastsim.asymmetries import a_parallel
    from polligen.xsec import EventSpinState
    plus = bk.SpinCategory("v+", 1.0, (1.0, 0.0, 0.0), lam_e=+1, pe=1.0)
    minus = bk.SpinCategory("v-", 1.0, (0.0, 0.0, 1.0), lam_e=+1, pe=1.0)
    expected = {}
    for name, m in (("v+", 1.0), ("v-", -1.0)):
        w, _, _ = sampler.kernel.amplitudes(
            sampler.tables, sampler.x_cells, sampler.q2_cells, sampler.s,
            EventSpinState(+1, 1.0, 1.0, m))
        expected[name] = float((sampler.xsec_flat * (1.0 + w)).sum())
    assert sampler.sigma_tot_pb(plus) == pytest.approx(expected["v+"],
                                                       rel=1e-12)
    assert sampler.sigma_tot_pb(minus) == pytest.approx(expected["v-"],
                                                        rel=1e-12)
    assert expected["v+"] != expected["v-"]
    t = sampler.tables
    y = sampler.q2_cells / (sampler.s * sampler.x_cells)
    apar_w = float((sampler.xsec_flat * a_parallel(
        t["g1"], t["f1"], y, sampler.x_cells, sampler.q2_cells)).sum())
    assert np.sign(expected["v+"] - expected["v-"]) == np.sign(apar_w)
    # event spin labels carry the m = +J ... -J convention
    ev = sampler.sample_category(plus, n=2000, rng=np.random.default_rng(9))
    assert np.all(ev["m"] == 1.0)


def test_mode_w_weights_reproduce_rate_ratio(sampler):
    """Reweighting an unpolarized sample must reproduce polarized rates."""
    rng = np.random.default_rng(4)
    ev = sampler.sample_category(UNPOL, n=200_000, rng=rng)
    tensor_cat = bk.SpinCategory("t0", 1.0,
                                 (0.15, 0.7, 0.15))  # m0-enriched
    w = sampler.weights_for(ev, [tensor_cat])[:, 0]
    expected = (sampler.sigma_tot_pb(tensor_cat)
                / sampler.sigma_tot_pb(UNPOL))
    assert w.mean() == pytest.approx(expected,
                                     abs=4 * w.std() / np.sqrt(w.size))
