"""Reconstructed-level pseudo-experiments (polligen.recopseudo): response
normalization, migration, exact expected counts, ratio-fit closure,
Delta bin-centering, and the coherent two-azimuth closure."""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen import bookkeeping as bk  # noqa: E402
from polligen import coherent as coh  # noqa: E402
from polligen import reco, recopseudo as rp  # noqa: E402
from polligen.sample import InclusiveSampler  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.polarized import toy_b1, toy_delta_gluon  # noqa: E402

CONFIG = beams.default_configs("6Li")[1]


@pytest.fixture(scope="module")
def response():
    kern = InclusiveKernel(
        beams.LI6, b1_func=toy_b1,
        delta_func=lambda x, q2, f1: toy_delta_gluon(x, q2, f1, scale=3e-2))
    analysis = fom.Scenario(lumi_fb_per_nucleon=10.0, pol_ion_tensor=0.6)
    gen = rp.generator_scenario(analysis)
    sampler = InclusiveSampler(kern, CONFIG, gen, nx=20, nq2=15,
                               x_range=(3e-3, 0.5), q2_range=(0.7, 100.0))
    resp = rp.RecoResponse(sampler, rp.RecoModel(), n_mc_per_cell=300,
                           rng=np.random.default_rng(1))
    return sampler, resp


def test_response_conserves_generator_rate(response):
    sampler, resp = response
    assert resp.w.sum() == pytest.approx(sampler.xsec_flat.sum(), rel=1e-9)
    assert np.all(resp.eff >= 0) and np.all(resp.eff <= 1)
    assert 0.3 < (resp.w * resp.eff).sum() / resp.w.sum() < 1.0


def test_low_y_bin_needs_hadronic_y(response):
    sampler, resp = response
    plan = bk.tensor_flip_plan(0.6)
    # a bin at y ~ 0.04-0.1 (Q2 ~ 1.4, x ~ 0.014 at s = 2000): away from
    # the y = 0.01 edge, where the reco cut truncates the bin (efficiency
    # 0.45 and a 9% amplitude shift for spot 1 of money plot 5R)
    xlo, xhi = 0.01, 0.02
    q2lo, q2hi = 1.0, 2.0
    mixed = resp.bin_summary(xlo, xhi, q2lo, q2hi, plan.categories[0])
    assert 0.6 < mixed["purity"] < 0.95
    assert 0.9 < mixed["a_reco_bin"] / mixed["a_true_bin"] < 1.05
    # 3 mrad direction resolution at theta' ~ 0.1 rad -> sigma_phi ~ 25 mrad
    # -> <cos 2 dphi> = 0.9987: negligible, but not 1
    assert 0.995 < mixed["dilution_phi"] < 0.9995
    # the electron method alone fails at the y ~ 0.01-0.03 edge bin
    ele = rp.RecoResponse(sampler, rp.RecoModel(y_method="electron"),
                          n_mc_per_cell=150, rng=np.random.default_rng(2))
    low = (0.03, 0.06, 1.0, 2.0)
    e_sum = ele.bin_summary(*low, plan.categories[0])
    m_sum = resp.bin_summary(*low, plan.categories[0])
    assert e_sum["purity"] < m_sum["purity"] - 0.15
    assert e_sum["efficiency"] < m_sum["efficiency"]


def test_expected_counts_match_generator_level_when_unsmeared(response):
    """With truth binning and no efficiency the exact expected counts
    reproduce the sampler's effective_modulation for the fill."""
    sampler, resp = response
    plan = bk.tensor_flip_plan(0.6)
    cat = plan.categories[0]
    mask_t = resp.mask_true(0.03, 0.1, 2.0, 10.0) & (resp.eff > 0)
    edges = np.linspace(0, 2 * np.pi, 25)
    mu = resp.expected_counts([cat], 1.0, mask_t, edges)
    # unsmeared reference: sum over cells of sigma_c(den + num cos2phi)
    cells = np.unique(resp.cell[mask_t])
    cmask = np.zeros(sampler.xsec_flat.size, dtype=bool)
    cmask[cells] = True
    sig, _a1, a2 = sampler.effective_modulation(cat, mask=cmask)
    # total rate: eff-weighted -> compare the modulation SHAPE only
    amp_mc = ((mu[0] * np.cos(2 * 0.5 * (edges[:-1] + edges[1:]))).sum()
              / mu[0].sum() * 2.0 / (np.sin(np.pi / 24) / (np.pi / 24)))
    assert amp_mc == pytest.approx(a2, rel=0.05)


def test_ratio_fit_closure_on_reco_bin(response):
    sampler, resp = response
    plan = bk.tensor_flip_plan(0.6)
    mask = resp.mask_reco(0.03, 0.1, 2.0, 10.0)
    summ = resp.bin_summary(0.03, 0.1, 2.0, 10.0, plan.categories[0])
    lumi_pb = 1.0e4
    eff = lambda ph: 1.0 + 0.04 * np.cos(2 * ph) + 0.02 * np.cos(ph)  # noqa: E731
    exact = rp.measure_inclusive(resp, plan, lumi_pb, mask, poisson=False,
                                 phi_eff=eff)
    assert exact["amp"] == pytest.approx(summ["a_reco_bin"], rel=2e-3)
    rng = np.random.default_rng(3)
    est = np.array([rp.measure_inclusive(resp, plan, lumi_pb, mask, rng=rng,
                                         phi_eff=eff)["amp"]
                    for _ in range(60)])
    n = exact["n"]
    analytic = reco.err_harmonic_ratio(n, [0.6, -1.2])
    assert est.mean() == pytest.approx(summ["a_reco_bin"],
                                       abs=4 * est.std() / np.sqrt(est.size))
    assert est.std() == pytest.approx(analytic, rel=0.35)


def test_delta_bin_centering_exact_without_noise(response):
    sampler, resp = response
    plan = bk.tensor_flip_plan(0.6)
    mask = resp.mask_reco(0.03, 0.1, 2.0, 10.0)
    summ = resp.bin_summary(0.03, 0.1, 2.0, 10.0, plan.categories[0])
    exact = rp.measure_inclusive(resp, plan, 1.0e4, mask, poisson=False)
    out = rp.delta_from_amplitude(exact, summ, delta_center=-0.05)
    assert out["delta"] == pytest.approx(-0.05, rel=2e-3)


def test_tensor_flip_plan_populations():
    plan = bk.tensor_flip_plan(0.6, share_plus=0.5, rel_lumi_offset=0.01)
    cats = plan.categories
    assert cats[0].moments()[1] == pytest.approx(0.6)
    assert cats[1].moments()[1] == pytest.approx(-1.2)
    assert all(abs(c.moments()[0]) < 1e-12 for c in cats)
    assert cats[0].lumi_fraction == pytest.approx(0.505)
    assert all(np.all(np.asarray(c.populations) >= 0) for c in cats)


# --- coherent -----------------------------------------------------------------

@pytest.fixture(scope="module")
def coherent_response():
    sc = coh.CoherentScenario()
    return sc, rp.CoherentResponse(sc, CONFIG, reco.SIGMA_THETA_HA,
                                   aspect=1.25, n_mc=200000,
                                   rng=np.random.default_rng(4))


def test_coherent_response_acceptance(coherent_response):
    sc, cr = coherent_response
    # rectangle with aspect 1.25 accepts less than the isotropic circle
    circ = sc.tag_acceptance_angular(reco.SIGMA_THETA_HA,
                                     CONFIG.ion_momentum_per_nucleon)
    assert 0.2 < cr.acceptance / circ < 1.0
    assert cr.t_reco.min() > 0.0
    # single-fill fake harmonic from the cutout is large (negative w.r.t.
    # the VERTICAL axis: the taller cutout accepts more along x)
    fake = np.mean(np.cos(2.0 * cr.beta_reco))
    assert fake < -0.3


def test_slot_cutout_flips_the_fake_harmonic():
    """The ePIC pots surround a horizontal slot (wide in x, tight in y):
    the tagged sample is dominated by vertical recoils, so the fake
    harmonic w.r.t. the vertical spin axis is large and POSITIVE, and the
    acceptance differs from the isotropic circular cut."""
    sc = coh.CoherentScenario()
    cr = rp.CoherentResponse(sc, CONFIG, reco.SIGMA_THETA_HA, aspect=1.0,
                             cut_scale_xy=(2.5, 1.0), n_mc=100000,
                             rng=np.random.default_rng(9))
    assert np.mean(np.cos(2.0 * cr.beta_reco)) > 0.3
    circ = sc.tag_acceptance_angular(reco.SIGMA_THETA_HA,
                                     CONFIG.ion_momentum_per_nucleon)
    assert 0.1 < cr.acceptance / circ < 1.0
    assert cr.cut_pt_xy[0] == pytest.approx(2.5 * cr.cut_pt_xy[1], rel=1e-6)


def test_two_azimuth_fit_closure(coherent_response):
    sc, cr = coherent_response
    plan = bk.tensor_flip_plan(0.6)
    a_t = lambda t: sc.a2_deformation(t, 1.0)  # noqa: E731  (per unit Pzz)
    fit = rp.measure_coherent(cr, 1e9, plan, 0.05, 0.12, a_e=0.01,
                              a_t_func=a_t, a_m=0.004, u1=0.05, u2=0.02,
                              kappa=0.002, poisson=False)
    tr = fit["truth"]
    # exact expected counts + acceptance-weighted basis: closure at the
    # level of the a_t(t) variation inside the t bin (~1e-3)
    assert fit["a_e"] == pytest.approx(0.01, rel=3e-3)
    assert fit["a_t"] == pytest.approx(tr["a_t"], rel=5e-3)
    assert fit["a_m"] == pytest.approx(tr["a_m"], rel=2e-2, abs=1e-4)
    assert fit["const"] == pytest.approx(0.002, abs=2e-4)
    # without the u correction the mixed term is contaminated at O(a_t u1)
    fit0 = reco.harmonic_ratio_fit_2d(fit["expected"], [0.5, 0.5],
                                      [0.6, -1.2], fit["alpha_edges"],
                                      fit["beta_edges"], u_coeffs=None,
                                      beta_means=fit["beta_means"])
    # (with the acceptance-shaped beta distribution the leakage reaches
    # a_e as well: ~9% here)
    assert abs(fit0["a_e"] / 0.01 - 1.0) > 2e-2
    assert abs(fit0["a_m"] - tr["a_m"]) > 5e-4
    # and with the ANALYTIC (uniform in-bin) basis the cutout's in-bin
    # acceptance shape biases a_t visibly
    fit1 = reco.harmonic_ratio_fit_2d(fit["expected"], [0.5, 0.5],
                                      [0.6, -1.2], fit["alpha_edges"],
                                      fit["beta_edges"], u_coeffs=(0.05, 0.02))
    assert abs(fit1["a_t"] / tr["a_t"] - 1.0) > 2e-3


def test_two_azimuth_fit_pulls(coherent_response):
    sc, cr = coherent_response
    plan = bk.tensor_flip_plan(0.6)
    a_t = lambda t: sc.a2_deformation(t, 1.0)  # noqa: E731
    rng = np.random.default_rng(5)
    vals = []
    for _ in range(40):
        fit = rp.measure_coherent(cr, 3e6, plan, 0.05, 0.12, a_e=0.01,
                                  a_t_func=a_t, u1=0.05, u2=0.02, rng=rng)
        vals.append((fit["a_t"], fit["err_t"], fit["a_e"], fit["err_e"]))
    vals = np.asarray(vals)
    tr = fit["truth"]
    assert vals[:, 0].mean() == pytest.approx(
        tr["a_t"], abs=4 * vals[:, 0].std() / np.sqrt(vals.shape[0]))
    assert vals[:, 0].std() == pytest.approx(vals[:, 1].mean(), rel=0.4)
    assert vals[:, 2].mean() == pytest.approx(
        0.01, abs=4 * vals[:, 2].std() / np.sqrt(vals.shape[0]))


def test_project_coherent_angular_cut():
    sc = coh.CoherentScenario()
    scen = fom.Scenario(lumi_fb_per_nucleon=10.0)
    proj, n_coh, tagged = coh.project_coherent(
        CONFIG, scen, sc, sigma_theta_list=(reco.SIGMA_THETA_HA,))
    key = "sigma_theta=73urad"
    ratio = tagged[key].sum() / n_coh.sum()
    cut = reco.tag_pt_cut(reco.SIGMA_THETA_HA, CONFIG.ion_momentum_per_nucleon)
    assert ratio == pytest.approx(np.exp(-sc.slope_b * cut * cut), rel=1e-9)
