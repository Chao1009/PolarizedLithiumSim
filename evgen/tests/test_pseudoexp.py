"""Gate 3 (plans/05 §5.4): pseudo-experiment estimator closure.

Full-chain closure of the Step-5.A deliverable: spin-labeled events ->
counting/fit estimators -> (i) pulls unbiased against the sigma-weighted
analytic asymmetries, (ii) spreads matching the fastsim error formulas,
(iii) phi-fit amplitude unbiased under a holey phi acceptance, and
(iv) the relative-luminosity systematic: the naive estimators are biased
by exactly the analytic first-order formula, the luminosity-corrected
estimators are not.
"""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen import bookkeeping as bk  # noqa: E402
from polligen import estimators as est  # noqa: E402
from polligen.sample import InclusiveSampler, run_pseudo_experiment  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.asymmetries import (  # noqa: E402
    a_parallel, azz, a_cos2phi, err_a_parallel, err_azz,
    err_cos2phi_amplitude)
from polli_fastsim.polarized import toy_b1, toy_delta_gluon  # noqa: E402


CONFIG = beams.default_configs("6Li")[1]


def _window_sampler(kern, x_range, q2_range):
    return InclusiveSampler(kern, CONFIG, fom.Scenario(), nx=10, nq2=8,
                            x_range=x_range, q2_range=q2_range)


def _sigma_weighted(sampler, values):
    return float(np.average(values, weights=sampler.xsec_flat))


@pytest.fixture(scope="module")
def vector_setup():
    kern = InclusiveKernel(beams.LI6)  # vector sector only
    sampler = _window_sampler(kern, (0.03, 0.3), (4.0, 60.0))
    t = sampler.tables
    y = sampler.q2_cells / (sampler.s * sampler.x_cells)
    truth = _sigma_weighted(sampler, a_parallel(t["g1"], t["f1"], y,
                                                sampler.x_cells,
                                                sampler.q2_cells))
    return sampler, truth


@pytest.fixture(scope="module")
def tensor_setup():
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    sampler = _window_sampler(kern, (0.03, 0.3), (4.0, 60.0))
    t = sampler.tables
    y = sampler.q2_cells / (sampler.s * sampler.x_cells)
    truth = _sigma_weighted(sampler, azz(t["b1"], t["f1"], t["f2"],
                                         sampler.x_cells, y))
    return sampler, truth


def _lumi_for(sampler, plan, n_target):
    sigma = sum(sampler.sigma_tot_pb(c) * c.lumi_fraction
                for c in plan.categories)
    return n_target / sigma


def test_apar_closure(vector_setup):
    sampler, truth = vector_setup
    pe, pz = 0.7, 0.6
    plan = bk.helicity_flip_plan(1.0, pz, pe)
    lumi = _lumi_for(sampler, plan, 60_000)
    rng = np.random.default_rng(11)
    vals, ns = [], []
    for _ in range(240):
        ev, _ = run_pseudo_experiment(sampler, plan, lumi, rng=rng)
        np_, nm = ev["apar+"]["x"].size, ev["apar-"]["x"].size
        vals.append(est.apar_flip(np_, nm, pe, pz))
        ns.append(np_ + nm)
    vals = np.array(vals)
    expected_err = err_a_parallel(np.mean(ns), pe, pz)
    # unbiased pulls
    assert vals.mean() == pytest.approx(
        truth, abs=4 * vals.std() / np.sqrt(vals.size))
    # spread matches the analytic error formula
    assert 0.85 < vals.std() / expected_err < 1.15


def test_azz_closure(tensor_setup):
    sampler, truth = tensor_setup
    pz, pzz = 0.7, 0.6
    plan = bk.tensor_thirds_plan(pz, pzz)
    lumi = _lumi_for(sampler, plan, 50_000)
    rng = np.random.default_rng(12)
    vals, ns = [], []
    for _ in range(280):
        ev, _ = run_pseudo_experiment(sampler, plan, lumi, rng=rng)
        n = {k: ev["azz" + k]["x"].size for k in "+-0"}
        vals.append(est.azz_thirds(n["+"], n["-"], n["0"], pzz))
        ns.append(sum(n.values()))
    vals = np.array(vals)
    expected_err = err_azz(np.mean(ns), pzz)
    assert vals.mean() == pytest.approx(
        truth, abs=4 * vals.std() / np.sqrt(vals.size))
    assert 0.85 < vals.std() / expected_err < 1.15


def test_azz_relative_lumi_bias(tensor_setup):
    """The plans/05 10^-4 rel-lumi systematic, validated at delta = 2%:
    naive thirds estimator biased by -(2/3) delta / Pzz; corrected clean."""
    sampler, truth = tensor_setup
    pz, pzz, delta = 0.7, 0.6, 0.02
    plan = bk.tensor_thirds_plan(pz, pzz, rel_lumi_offset=delta)
    lumi = _lumi_for(sampler, plan, 90_000)
    rng = np.random.default_rng(13)
    naive, corrected = [], []
    for _ in range(140):
        ev, shares = run_pseudo_experiment(sampler, plan, lumi, rng=rng)
        n = {k: ev["azz" + k]["x"].size for k in "+-0"}
        naive.append(est.azz_thirds(n["+"], n["-"], n["0"], pzz))
        corrected.append(est.azz_thirds(
            n["+"], n["-"], n["0"], pzz,
            lumis=(shares["azz+"], shares["azz-"], shares["azz0"])))
    naive, corrected = np.array(naive), np.array(corrected)
    se = naive.std() / np.sqrt(naive.size)
    bias = bk.azz_rel_lumi_bias(delta, pzz)
    assert abs(bias) > 10 * se  # the test can actually resolve the bias
    assert naive.mean() - truth == pytest.approx(bias, abs=5 * se)
    assert corrected.mean() == pytest.approx(truth, abs=5 * se)


def test_apar_relative_lumi_bias(vector_setup):
    sampler, truth = vector_setup
    pe, pz, delta = 0.7, 0.6, 0.02
    plan = bk.helicity_flip_plan(1.0, pz, pe, rel_lumi_offset=delta)
    lumi = _lumi_for(sampler, plan, 80_000)
    rng = np.random.default_rng(14)
    naive, corrected = [], []
    for _ in range(140):
        ev, shares = run_pseudo_experiment(sampler, plan, lumi, rng=rng)
        np_, nm = ev["apar+"]["x"].size, ev["apar-"]["x"].size
        naive.append(est.apar_flip(np_, nm, pe, pz))
        corrected.append(est.apar_flip(np_, nm, pe, pz,
                                       l_plus=shares["apar+"],
                                       l_minus=shares["apar-"]))
    naive, corrected = np.array(naive), np.array(corrected)
    se = naive.std() / np.sqrt(naive.size)
    bias = bk.apar_rel_lumi_bias(delta, pe, pz)
    assert abs(bias) > 10 * se
    assert naive.mean() - truth == pytest.approx(bias, abs=5 * se)
    assert corrected.mean() == pytest.approx(truth, abs=5 * se)


@pytest.fixture(scope="module")
def cos2phi_setup():
    kern = InclusiveKernel(
        beams.LI6,
        delta_func=lambda x, q2, f1: toy_delta_gluon(x, q2, f1, scale=0.05))
    sampler = _window_sampler(kern, (0.02, 0.12), (8.0, 80.0))
    t = sampler.tables
    y = sampler.q2_cells / (sampler.s * sampler.x_cells)
    truth = _sigma_weighted(sampler, a_cos2phi(t["delta"], t["f1"], t["f2"],
                                               sampler.x_cells, y))
    return sampler, truth


def test_cos2phi_moment_closure(cos2phi_setup):
    sampler, truth = cos2phi_setup
    pzz = 0.8
    plan = bk.transverse_tensor_plan(pzz, phi_s=0.6)
    cat = plan.categories[0]
    rng = np.random.default_rng(15)
    n = 40_000
    vals = []
    for _ in range(260):
        ev = sampler.sample_category(cat, n=n, rng=rng)
        vals.append(est.cos2phi_moment(ev["phi"] - cat.phi_s, pzz))
    vals = np.array(vals)
    expected_err = err_cos2phi_amplitude(n, pzz)
    assert vals.mean() == pytest.approx(
        truth, abs=4 * vals.std() / np.sqrt(vals.size))
    assert 0.85 < vals.std() / expected_err < 1.15


def test_cos2phi_fit_unbiased_with_holey_acceptance(cos2phi_setup):
    """Amplitude recovery with two dead phi sectors (gate: phi-modulation
    recovery with uniform AND holey acceptance)."""
    sampler, truth = cos2phi_setup
    pzz = 0.8
    plan = bk.transverse_tensor_plan(pzz)
    cat = plan.categories[0]
    holes = ((0.4, 1.0), (3.3, 4.1))  # asymmetric in 2*phi -> moment breaks

    def acceptance(phi):
        ok = np.ones_like(phi, dtype=bool)
        for lo, hi in holes:
            ok &= ~((phi >= lo) & (phi <= hi))
        return ok

    rng = np.random.default_rng(16)
    n = 60_000
    fit_vals, moment_vals = [], []
    for _ in range(120):
        ev = sampler.sample_category(cat, n=n, rng=rng)
        phi = ev["phi"] - cat.phi_s
        keep = acceptance(np.mod(phi, 2 * np.pi))
        fit_vals.append(est.cos2phi_fit(phi[keep], pzz,
                                        acceptance=acceptance))
        moment_vals.append(est.cos2phi_moment(phi[keep], pzz))
    fit_vals = np.array(fit_vals)
    moment_vals = np.array(moment_vals)
    se = fit_vals.std() / np.sqrt(fit_vals.size)
    assert fit_vals.mean() == pytest.approx(truth, abs=5 * se)
    # the naive moment IS biased under these holes (sanity of the gate)
    assert abs(moment_vals.mean() - truth) > 10 * se


def test_polarimetry_bookkeeping():
    plan = bk.tensor_thirds_plan(0.7, 0.6)
    assert plan.measured["pzz"] == 0.6  # no smearing by default
    smeared = bk.RunPlan(plan.categories, pz_true=0.7, pzz_true=0.6,
                         delta_p_over_p=0.03, polarimetry_seed=5)
    assert smeared.measured["pzz"] != 0.6
    again = bk.RunPlan(plan.categories, pz_true=0.7, pzz_true=0.6,
                       delta_p_over_p=0.03, polarimetry_seed=5)
    assert smeared.measured["pzz"] == again.measured["pzz"]
