"""Coherent (intact-6Li) channel: scenario model, tagging, phi pseudo-exp."""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen import bookkeeping as bk  # noqa: E402
from polligen import coherent as coh  # noqa: E402
from polligen.estimators import cos2phi_fit_binned  # noqa: E402
from polligen.sample import InclusiveSampler, phi_histogram_pseudo  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.asymmetries import a_cos2phi  # noqa: E402
from polli_fastsim.farforward import HIGH_ACCEPTANCE, HIGH_DIVERGENCE  # noqa: E402
from polli_fastsim.polarized import toy_delta_gluon  # noqa: E402


def test_gaussian_slope_li6_magnitude():
    b = coh.gaussian_slope(2.4)
    assert 45.0 < b < 55.0  # R^2/3 with R ~ 2.4 fm


def test_coherent_fraction_shape():
    sc = coh.CoherentScenario()
    f = sc.coherent_fraction
    assert f(1e-4) == pytest.approx(sc.f0, rel=1e-3)
    assert f(sc.x_coh) == pytest.approx(sc.f0 / 2.0)
    x = np.logspace(-4, 0, 40)
    assert np.all(np.diff(f(x)) < 0)
    assert f(0.1) < sc.f0 / 50.0


def test_sample_t_matches_slope_and_acceptance():
    sc = coh.CoherentScenario(slope_b=50.0)
    rng = np.random.default_rng(1)
    t = sc.sample_t(400000, rng)
    assert np.mean(t) == pytest.approx(1.0 / 50.0, rel=0.02)
    cut = HIGH_ACCEPTANCE.pt_cut_near_beam
    frac = np.mean(t > cut * cut)
    assert frac == pytest.approx(sc.tag_acceptance(cut), rel=0.02)


def test_tag_acceptance_sampled_through_router():
    """Since 2026-08-25 the router's near-beam cut is ANGULAR, so the
    sampled acceptance must reproduce `tag_acceptance_angular` at the
    configuration's beam momentum -- not the constant 0.20 GeV proton
    number, which for 6Li at 50 GeV/u is 0.218 GeV and 1.45x smaller
    acceptance (code review S8)."""
    sc = coh.CoherentScenario(slope_b=50.0)
    p_u = beams.default_configs("6Li")[1].ion_momentum_per_nucleon
    acc = coh.tag_acceptance_sampled(sc, HIGH_ACCEPTANCE, p_u, n=200000,
                                     rng=np.random.default_rng(2))
    assert acc == pytest.approx(
        float(sc.tag_acceptance_angular(HIGH_ACCEPTANCE.sigma_theta, p_u)),
        rel=0.03)
    assert acc < sc.tag_acceptance(0.20)
    acc_hd = coh.tag_acceptance_sampled(sc, HIGH_DIVERGENCE, p_u, n=200000,
                                        rng=np.random.default_rng(3))
    assert acc_hd < 1e-4


def test_angular_cut_scales_with_the_beam_momentum():
    """The whole point of S8: the same optics gives a very different tag
    acceptance at the three beam energies, because the envelope is an
    angle and the recoil momentum is A p_u."""
    sc = coh.CoherentScenario(slope_b=50.0)
    sig = HIGH_ACCEPTANCE.sigma_theta
    accs = [float(sc.tag_acceptance_angular(sig, p_u))
            for p_u in (20.5, 50.0, 137.5)]
    assert accs[0] == pytest.approx(0.67, abs=0.02)
    assert accs[1] == pytest.approx(0.093, abs=0.005)
    assert accs[2] < 1e-7
    assert accs[0] > accs[1] > accs[2]


def test_a2_deformation_scaling():
    sc = coh.CoherentScenario(slope_b=50.0, eps_b0=-0.08)
    # linear in |t| and in pzz; sign: eps_b0 < 0 (6Li) with pzz > 0 -> +
    a = sc.a2_deformation(0.1, 0.6)
    assert a == pytest.approx(0.6 / 4.0 * 0.08 * 50.0 * 0.1)
    assert a > 0
    assert sc.a2_deformation(0.2, 0.6) == pytest.approx(2 * a)
    assert sc.a2_deformation(0.1, 0.0) == 0.0
    # tagged mean uses <|t|> = cut^2 + 1/B exactly (linear regime)
    cut = HIGH_ACCEPTANCE.pt_cut_near_beam
    assert sc.a2_tagged(cut, 0.6) == pytest.approx(
        float(sc.a2_deformation(sc.mean_t_tagged(cut), 0.6)))


def test_mantysaari_table_m_state_relation():
    # wave-function symmetry: a_2(0) ~ -2 a_2(+-1) in the linear regime
    for t, (a0, a1) in coh.MANTYSAARI_A2_DEUTERON.items():
        if t <= 0.2:
            assert a0 == pytest.approx(-2.0 * a1, rel=0.25)
        assert a0 * a1 < 0  # opposite signs at all tabulated |t|


def test_recoil_stays_in_near_beam_band():
    lab = coh.recoil_lab(np.array([0.01, 0.05, 0.2]), 0.0, 100.0,
                         x_pom=0.01)
    assert np.all(np.abs(lab["R"] - 1.0) < 0.05)
    assert np.all(lab["theta"] < 5e-3)


def test_project_coherent_rates_nested():
    config = beams.default_configs("6Li")[1]
    sc = coh.CoherentScenario()
    proj, n_coh, tagged = coh.project_coherent(
        config, fom.Scenario(lumi_fb_per_nucleon=100.0), sc,
        optics_list=(HIGH_ACCEPTANCE, HIGH_DIVERGENCE))
    assert np.all(n_coh <= proj.n_events + 1e-9)
    n_tag = tagged[HIGH_ACCEPTANCE.name]
    assert np.all(n_tag <= n_coh + 1e-9)
    assert np.all(tagged[HIGH_DIVERGENCE.name] <= n_tag + 1e-9)
    # at x << x_coh the suppression is just f0 * exp(-B cut^2)
    lowx = proj.accepted & (proj.x < 1e-3)
    ratio = n_tag[lowx] / np.maximum(proj.n_events[lowx], 1e-30)
    assert ratio == pytest.approx(
        sc.f0 * sc.tag_acceptance(0.20), rel=0.05)


def test_veto_table_routing():
    table = coh.veto_table()
    frag = {name: {f: (r, dest) for f, r, dest in rows}
            for name, rows in table.items()}
    r_a, dest_a = frag["alpha+d"]["alpha"]
    r_d, dest_d = frag["alpha+d"]["d"]
    assert r_a == pytest.approx(1.0) and "beam-blind" in dest_a
    assert r_d == pytest.approx(1.0) and "beam-blind" in dest_d
    r_he3, dest_he3 = frag["3He+t"]["3He"]
    assert r_he3 == pytest.approx(0.75) and dest_he3 == "RomanPots"
    r_t, dest_t = frag["3He+t"]["t"]
    assert r_t == pytest.approx(1.5) and "over-rigid" in dest_t
    assert frag["alpha+p+n"]["n"][1] == "ZDC"
    assert frag["alpha+p+n"]["p"][0] == pytest.approx(0.5)
    assert frag["alpha+p+n"]["p"][1] == "OMD"


# --- binned phi pseudo-experiments (full-luminosity projections) ----------


def test_phi_histogram_exact_amplitude_recovery():
    counts, edges = phi_histogram_pseudo(1e8, a2=0.03, nbins=36,
                                         poisson=False)
    assert counts.sum() == pytest.approx(1e8, rel=1e-12)
    amp = cos2phi_fit_binned(counts, edges, pzz=0.6)
    assert amp == pytest.approx(0.03 / 0.6, rel=1e-9)


def test_phi_histogram_poisson_spread():
    n, a2, pzz = 4e6, 0.01, 0.6
    rng = np.random.default_rng(7)
    amps = []
    for _ in range(80):
        counts, edges = phi_histogram_pseudo(n, a2=a2, rng=rng)
        amps.append(cos2phi_fit_binned(counts, edges, pzz))
    amps = np.asarray(amps)
    expect = np.sqrt(2.0 / n) / pzz
    assert np.mean(amps) == pytest.approx(a2 / pzz, abs=3 * expect / 9)
    assert np.std(amps) == pytest.approx(expect, rel=0.35)


def test_effective_modulation_matches_analytic():
    kern = InclusiveKernel(
        beams.LI6,
        delta_func=lambda x, q2, f1: toy_delta_gluon(x, q2, f1, scale=1e-2))
    config = beams.default_configs("6Li")[1]
    sampler = InclusiveSampler(kern, config, fom.Scenario(), nx=16, nq2=12,
                               x_range=(3e-3, 0.3), q2_range=(1.5, 50.0))
    pzz = 0.6
    plan = bk.transverse_tensor_plan(pzz)
    cat = plan.categories[0]
    sigma_pb, a1_eff, a2_eff = sampler.effective_modulation(cat)
    assert sigma_pb == pytest.approx(sampler.sigma_tot_pb(cat), rel=1e-12)
    assert a1_eff == 0.0
    # b1 = 0 here, so the rate weight is just sigma_cell and a2_eff is the
    # rate-weighted pzz * a_cos2phi over the accepted cells
    t = sampler.tables
    x, q2 = sampler.x_cells, sampler.q2_cells
    y = q2 / (sampler.s * x)
    amp = a_cos2phi(t["delta"], t["f1"], t["f2"], x, y)
    expect = pzz * float((sampler.xsec_flat * amp).sum()
                         / sampler.xsec_flat.sum())
    assert a2_eff == pytest.approx(expect, rel=1e-9)


def test_cos2phi_fit_binned_raises_on_dead_acceptance():
    counts, edges = phi_histogram_pseudo(1e5, a2=0.02, poisson=False)
    with pytest.raises(ValueError):
        cos2phi_fit_binned(counts, edges, 0.6,
                           acceptance=lambda phi: np.zeros_like(
                               np.asarray(phi), dtype=bool))
