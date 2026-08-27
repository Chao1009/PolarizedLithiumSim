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
# The near-beam envelope these machinery tests run at.  Corrected 2026-08-27
# (plans/10): with the machine's real per-configuration divergence and the
# gamma-matched ion energies, the coherent tag has an acceptance of ~5e-7 at
# best, so a test at the nominal optics would be fitting empty histograms.
# `sigma_theta_tagging` is the analytic tagging optimum -- the beta* at which
# L x acceptance peaks, giving acc = 1/e -- and is the only working point at
# which this measurement exists.  These tests check the FIT MACHINERY, so
# they run there.
SIG_TAG = reco.sigma_theta_tagging(CONFIG)


def _sigma_for_cut(cfg, pt_cut, n_sigma=10.0):
    """The sigma_theta whose n-sigma near-beam cut is `pt_cut` [GeV].

    Several tests below are about what a TIGHT CUTOUT DOES -- the fake
    cos 2beta harmonic, the fill-dependent envelope, the u2 leakage, the
    folded fit against a wrong prior.  Those pathologies scale with the
    cutout severity, not with the beam energy, so they are pinned to the
    0.218 GeV cut they were written at (the old 50 GeV/u working point)
    rather than to a sigma_theta that the 2026-08-27 corrections moved."""
    return pt_cut / (n_sigma * cfg.ion.A * cfg.ion_momentum_per_nucleon)


SIG_TIGHT = _sigma_for_cut(CONFIG, 0.218)


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
    return sc, rp.CoherentResponse(sc, CONFIG, SIG_TIGHT,
                                   aspect=1.25, n_mc=200000,
                                   rng=np.random.default_rng(4))


def test_coherent_response_acceptance(coherent_response):
    sc, cr = coherent_response
    # rectangle with aspect 1.25 accepts less than the isotropic circle
    circ = sc.tag_acceptance_angular(SIG_TIGHT,
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
    cr = rp.CoherentResponse(sc, CONFIG, SIG_TIGHT, aspect=1.0,
                             cut_scale_xy=(2.5, 1.0), n_mc=100000,
                             rng=np.random.default_rng(9))
    assert np.mean(np.cos(2.0 * cr.beta_reco)) > 0.3
    circ = sc.tag_acceptance_angular(SIG_TIGHT,
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
        CONFIG, scen, sc, sigma_theta_list=(SIG_TIGHT,))
    key = "sigma_theta=%.0furad" % (1e6 * SIG_TIGHT)
    ratio = tagged[key].sum() / n_coh.sum()
    cut = reco.tag_pt_cut(SIG_TIGHT, CONFIG.ion_momentum_per_nucleon)
    assert ratio == pytest.approx(np.exp(-sc.slope_b * cut * cut), rel=1e-9)


# --- what the coherent ratio does NOT cancel (code review F1, F6) -----------

@pytest.fixture(scope="module")
def slot_response():
    """The 6R default: slot-like cutout, isotropic divergence."""
    sc = coh.CoherentScenario(amp=0.01, eps_b0=-0.08)
    return sc, rp.CoherentResponse(sc, CONFIG, SIG_TIGHT,
                                   aspect=1.0, cut_scale_xy=(2.5, 1.0),
                                   n_mc=300000, rng=np.random.default_rng(11))


def _slot_fit(cr, sc, responses=None, u_assumed=None, lumi_assumed=None,
              n=8e7):
    plan = bk.tensor_flip_plan(0.6)
    a_t = lambda t: sc.cos2phi_coefficient_deformation(t, 1.0)  # noqa: E731
    fit = rp.measure_coherent(cr, n, plan, 0.05, 0.08, a_e=0.010,
                              a_t_func=a_t, u1=0.05, u2=0.02, poisson=False,
                              responses=responses, u_coeffs_assumed=u_assumed,
                              lumi_assumed=lumi_assumed)
    return fit, fit["truth"]


def test_with_cut_is_a_view_of_the_same_recoils(slot_response):
    _sc, cr = slot_response
    same = cr.with_cut(cr.cut_scale_xy)
    np.testing.assert_array_equal(same.t_reco, cr.t_reco)
    np.testing.assert_array_equal(same.beta_reco, cr.beta_reco)
    assert same.acceptance == cr.acceptance
    assert cr.cut_scale_xy == (2.5, 1.0)          # the original is untouched
    tighter = cr.with_cut((2.5, 1.02))
    assert tighter.acceptance < cr.acceptance
    assert tighter.cut_pt_xy[1] == pytest.approx(1.02 * cr.cut_pt_xy[1],
                                                 rel=1e-9)


def test_fill_dependent_envelope_biases_a_t(slot_response):
    """The coherent counterpart of F1.  The spin-state ratio cancels a
    COMMON cutout exactly; a difference between the fills does not, and
    the slot amplifies it far beyond the naive d<cos 2beta>/(P+ - P0):
    only half the beta bins are live and the t-shape template is ~99%
    anti-correlated with the constant there."""
    sc, cr = slot_response
    ref, truth = _slot_fit(cr, sc)
    assert ref["a_t"] == pytest.approx(truth["a_t"], rel=1e-5)   # closure
    cov = ref["cov"]
    d = np.sqrt(np.diag(cov))
    assert cov[0, 2] / (d[0] * d[2]) < -0.9      # const/a_t degeneracy

    bias = {}
    for delta in (1e-4, 1e-3, 1e-2):
        f, _ = _slot_fit(cr, sc, responses=[cr.with_cut((2.5, 1.0 + delta)),
                                            cr])
        bias[delta] = (f["a_t"] - truth["a_t"]) / truth["a_t"]
    assert 0.0 < bias[1e-4] < bias[1e-3] < bias[1e-2]
    # 1e-3 of the vertical envelope already moves a_t by more than 10%,
    # i.e. by more than the one-year statistical error -- the documented
    # "1% envelope -> 1.3% a_t" estimate is two orders of magnitude out
    assert bias[1e-3] > 0.10
    assert bias[1e-2] > 1.0


def test_assumed_u2_leaks_into_a_e_at_first_order(slot_response):
    """u2 is NOT a second-order nuisance: an error in it reaches a_e as
    a_t x du2 x <cos 2beta>, and the slot's fake <cos 2beta> = +0.77 with
    a_t/a_e ~ 12 makes a ZEUS-1sigma du2 a ~20% shift of a_e."""
    sc, cr = slot_response
    ref, truth = _slot_fit(cr, sc)
    fake_c2 = float(np.mean(np.cos(2.0 * cr.beta_reco)))
    assert fake_c2 > 0.5

    for du2 in (-0.024, +0.024):
        f, _ = _slot_fit(cr, sc, u_assumed=(0.05, 0.02 + du2))
        shift = f["a_e"] - 0.010
        assert np.sign(shift) == np.sign(du2)
        assert abs(shift / 0.010) == pytest.approx(0.20, abs=0.06)
        # first order: the estimate a_t * du2 * <cos 2beta> gets the size
        assert abs(shift) == pytest.approx(
            abs(truth["a_t"] * du2 * fake_c2), rel=0.35)
        # a_t itself is immune
        assert abs(f["a_t"] / truth["a_t"] - 1.0) < 5e-3

    # linearity: a ten-times-smaller error gives a ten-times-smaller shift
    small, _ = _slot_fit(cr, sc, u_assumed=(0.05, 0.02 + 0.0024))
    big, _ = _slot_fit(cr, sc, u_assumed=(0.05, 0.02 + 0.024))
    assert (big["a_e"] - 0.010) / (small["a_e"] - 0.010) == \
        pytest.approx(10.0, rel=0.05)


def test_assumed_u1_moves_the_mixed_term_not_a_e(slot_response):
    sc, cr = slot_response
    ref, truth = _slot_fit(cr, sc)
    f, _ = _slot_fit(cr, sc, u_assumed=(0.05 + 0.03, 0.02))
    assert abs(f["a_e"] / 0.010 - 1.0) < 0.03           # a_e barely moves
    assert abs(f["a_m"] - ref["a_m"]) > 1e-3            # a_m does


def test_assumed_luminosity_is_second_order_in_the_coherent_ratio(slot_response):
    """As on the inclusive side, a relative-luminosity error enters the
    ratio only through Pbar: a constant offset plus a second-order
    rescaling."""
    sc, cr = slot_response
    ref, truth = _slot_fit(cr, sc)
    f, _ = _slot_fit(cr, sc, lumi_assumed=[0.5 * 1.001, 0.5 * 0.999])
    assert abs(f["a_t"] / truth["a_t"] - 1.0) < 2e-3
    assert abs(f["a_e"] / 0.010 - 1.0) < 5e-3
    assert abs(f["const"] - ref["const"]) > 1e-4        # absorbed by kappa


def test_measure_coherent_rejects_a_wrong_number_of_responses(slot_response):
    sc, cr = slot_response
    plan = bk.tensor_flip_plan(0.6)
    with pytest.raises(ValueError):
        rp.measure_coherent(cr, 1e7, plan, 0.05, 0.08, 0.01,
                            lambda t: 0.1 * t, responses=[cr])


# --- the azimuthal-alignment null test (plans/08 A2) -----------------------

def _rotated_fit(cr, sc, delta, roll_only=False, drop_s2t=False, n=8e7):
    """Generate with the tensor modulation rotated by `delta` and fit with
    the sin columns.  A spin-axis error rotates BOTH tensor harmonics; a
    Roman-Pot azimuthal roll rotates only the recoil one."""
    plan = bk.tensor_flip_plan(0.6)
    base = lambda t: sc.cos2phi_coefficient_deformation(t, 1.0)  # noqa: E731
    c, s_ = np.cos(2 * delta), np.sin(2 * delta)
    kw = dict(u1=0.05, u2=0.02, poisson=False, with_sin=True,
              a_t_s_func=lambda t: base(t) * s_)
    if not roll_only:
        kw["a_e_s"] = 0.010 * s_
    if drop_s2t:
        be = np.linspace(0.0, 2.0 * np.pi, 25)
        ae = np.linspace(0.0, 2.0 * np.pi, 13)
        pzz = [cat.moments()[1] for cat in plan.categories]
        frac = [cat.lumi_fraction for cat in plan.categories]
        mu = cr.expected_counts_2d(n, pzz, frac, 0.05, 0.08, ae, be,
                                   0.010 if roll_only else 0.010 * c,
                                   lambda t: base(t) * c, u1=0.05, u2=0.02,
                                   a_t_s_func=lambda t: base(t) * s_)
        bm = dict(cr.basis_means(0.05, 0.08, be))
        bm.pop("s2t")
        return reco.harmonic_ratio_fit_2d(mu, frac, pzz, ae, be,
                                          u_coeffs=(0.05, 0.02),
                                          beta_means=bm, with_sin=True)
    return rp.measure_coherent(cr, n, plan, 0.05, 0.08,
                               0.010 if roll_only else 0.010 * c,
                               lambda t: base(t) * c, **kw)


def test_sin_nulls_vanish_and_leave_the_cos_sector_alone(slot_response):
    """Unpolarized leptons + a headless axis forbid sin 2a, sin 2b and
    sin(a+b) exactly, so all three must come out zero -- and turning the
    columns on must not perturb the coefficients that are measured."""
    sc, cr = slot_response
    off, _ = _slot_fit(cr, sc)
    on = rp.measure_coherent(cr, 8e7, bk.tensor_flip_plan(0.6), 0.05, 0.08,
                             a_e=0.010,
                             a_t_func=lambda t: sc.cos2phi_coefficient_deformation(t, 1.0),
                             u1=0.05, u2=0.02, poisson=False, with_sin=True)
    for key in ("a_e_s", "a_t_s", "a_m_s"):
        assert abs(on[key]) < 1e-8
    assert on["a_e"] == pytest.approx(off["a_e"], abs=1e-12)
    assert on["a_t"] == pytest.approx(off["a_t"], abs=1e-9)
    assert on["err_t"] == pytest.approx(off["err_t"], rel=1e-3)


def test_spin_axis_error_rotates_both_tensor_harmonics(slot_response):
    sc, cr = slot_response
    for delta in (0.02, 0.05):
        f = _rotated_fit(cr, sc, delta)
        assert f["a_e_s"] / f["a_e"] == pytest.approx(np.tan(2 * delta),
                                                      rel=1e-3)
        assert f["a_t_s"] / f["a_t"] == pytest.approx(np.tan(2 * delta),
                                                      rel=1e-3)


def test_pot_roll_shows_only_in_the_recoil_harmonic(slot_response):
    """The signature that separates the two misalignments: a Roman-Pot
    azimuthal roll leaves the lepton-plane harmonic exactly alone."""
    sc, cr = slot_response
    for delta in (0.02, 0.05):
        f = _rotated_fit(cr, sc, delta, roll_only=True)
        assert abs(f["a_e_s"] / f["a_e"]) < 1e-5
        assert f["a_t_s"] / f["a_t"] == pytest.approx(np.tan(2 * delta),
                                                      rel=1e-3)


def test_null_closure_needs_the_t_weighted_sine_template(slot_response):
    """`t_s` must use the same t-shape template as `t`: with the plain
    <sin 2beta> the null misses its own closure by ~5%."""
    sc, cr = slot_response
    good = _rotated_fit(cr, sc, 0.05, roll_only=True)
    bad = _rotated_fit(cr, sc, 0.05, roll_only=True, drop_s2t=True)
    tan2d = np.tan(0.10)
    assert good["a_t_s"] / good["a_t"] == pytest.approx(tan2d, rel=1e-3)
    assert abs(bad["a_t_s"] / bad["a_t"] / tan2d - 1.0) > 0.03


def test_basis_2d_sin_columns_are_orthogonal_to_the_cos_ones():
    """On a full uniform grid the six columns are mutually orthogonal, so
    with_sin cannot inflate the cos errors."""
    ae = np.linspace(0.0, 2.0 * np.pi, 13)
    be = np.linspace(0.0, 2.0 * np.pi, 25)
    b = reco.basis_2d(ae, be)
    keys = ("e", "t", "m", "e_s", "t_s", "m_s")
    for i, k1 in enumerate(keys):
        for k2 in keys[i + 1:]:
            assert abs(float(np.dot(b[k1], b[k2]))) < 1e-10


# --- importance sampling above the near-beam envelope (plans/08 A4) --------

def test_t_floor_reproduces_the_acceptance_and_reaches_tight_envelopes():
    """The recoil spectrum is exponential, so sampling from |t_floor|
    upwards has a CONSTANT likelihood ratio exp(-B t_floor): the estimate
    is unchanged and the accepted sample grows.  Without it a tight
    envelope is not simulable at all -- the 0.60 GeV cut of the
    top-energy configuration leaves zero accepted recoils."""
    sc = coh.CoherentScenario(slope_b=50.0)
    kw = dict(shape="ellipse", n_mc=200000, rng=np.random.default_rng(3))
    plain = rp.CoherentResponse(sc, CONFIG, SIG_TIGHT, **kw)
    lifted = rp.CoherentResponse(sc, CONFIG, SIG_TIGHT,
                                 t_floor=0.25 * plain.pt_cut ** 2, **kw)
    assert lifted.acceptance == pytest.approx(plain.acceptance, rel=0.05)
    assert lifted.t_true.size > 1.5 * plain.t_true.size
    assert lifted._t_weight < 1.0 and plain._t_weight == 1.0
    # the analytic value for the circular cut, within the smearing effect
    analytic = float(sc.tag_acceptance_angular(
        SIG_TIGHT, CONFIG.ion_momentum_per_nucleon))
    assert 0.8 < lifted.acceptance / analytic < 1.5

    # the top configuration at its REAL divergence (plans/10): 92 urad on a
    # 825 GeV 6Li is a 0.76 GeV envelope, and a plain sampler finds nothing
    top = beams.default_configs("6Li")[2]
    sig_top = reco.sigma_theta_for(top)[0]
    dead = rp.CoherentResponse(sc, top, sig_top, **kw)
    assert dead.t_true.size == 0
    alive = rp.CoherentResponse(sc, top, sig_top,
                                t_floor=0.30, **kw)
    assert alive.t_true.size > 1000
    assert 0.0 < alive.acceptance < 1e-6


def test_t_floor_weights_leave_the_fit_unbiased():
    """The per-event weight is uniform, so nothing downstream of the
    sampler needs to know about it."""
    sc = coh.CoherentScenario(amp=0.01, eps_b0=-0.08)
    plan = bk.tensor_flip_plan(0.6)
    a_t = lambda t: sc.cos2phi_coefficient_deformation(t, 1.0)  # noqa: E731
    cr = rp.CoherentResponse(sc, CONFIG, SIG_TAG,
                             cut_scale_xy=(2.5, 1.0), n_mc=200000,
                             rng=np.random.default_rng(11), t_floor=0.012)
    fit = rp.measure_coherent(cr, 8e7, plan, 0.05, 0.08, a_e=0.010,
                              a_t_func=a_t, u1=0.05, u2=0.02, poisson=False)
    assert fit["a_t"] == pytest.approx(fit["truth"]["a_t"], rel=1e-4)
    assert fit["a_e"] == pytest.approx(0.010, rel=1e-4)


# --- the forward-folded response replaces the bin-by-bin K (plans/08 A6) ----

def _truth_delta(x, q2, f1):
    """The fixture's injected shape, interpretation A of the repository's
    ansatz: Delta = scale * F1(x, Q2) * x^0.3 (1-x)^4."""
    return toy_delta_gluon(x, q2, f1, scale=3e-2)


def _prior_delta(x, q2, f1):
    """Interpretation B of the same ansatz -- the same x shape WITHOUT
    the F1 factor (delta_models module docstring).  A and B is the
    repository's own Delta ambiguity, so it is the wrong prior the
    closure tests below correct with."""
    x = np.asarray(x, dtype=float)
    return 1e-2 * x ** 0.3 * (1.0 - x) ** 4


def _slice_bins(resp, cat, q2lo, q2hi, n_min=50):
    """Every x bin of one Q2 slice with enough response MC, as the bin
    dicts fold_shape_fit consumes (amplitude = the exact reco-bin one)."""
    xe = resp.sampler.proj.x_edges
    out = []
    for i in range(xe.size - 1):
        mask = resp.mask_reco(xe[i], xe[i + 1], q2lo, q2hi)
        if mask.sum() < n_min:
            continue
        summ = resp.bin_summary(xe[i], xe[i + 1], q2lo, q2hi, cat)
        if not np.isfinite(summ["a_reco_bin"]):
            continue
        out.append({"mask": mask, "x": np.sqrt(xe[i] * xe[i + 1]),
                    "q2": np.sqrt(q2lo * q2hi), "amp": summ["a_reco_bin"],
                    "err": 0.02 * abs(summ["a_reco_bin"]), "summ": summ})
    return out


def test_amplitude_is_exactly_linear_in_delta(response):
    """dA/dDelta against the published master formula, written out here
    from the paper rather than taken from the code: for the transverse
    alignment axis of the flip plan,

        A = a_2/P_zz = -[(1-y)/y^2] <c_eff> sin^2(theta_S) Delta / D_phi,
        D_phi = F1 + (1-y)/(x y^2) F2,  c_eff(m) = 3 m^2 - 2 for J = 1

    (Hoodbhoy-Jaffe-Manohar NPB 312:571 (1989) Eq. (30); the same
    normalization as polli_fastsim.asymmetries.a_cos2phi).  Delta appears
    linearly and nowhere else, so the derivative is Delta-independent and
    the forward fold is exact."""
    sampler, resp = response
    cat = bk.tensor_flip_plan(0.6).categories[0]
    x, q2 = sampler.x_cells, sampler.q2_cells
    y = q2 / (sampler.s * x)
    kern = sampler.kernel
    f1 = kern.nf2.f1a(x, q2) / kern.ion.A
    f2 = kern.nf2.f2a(x, q2) / kern.ion.A
    d_phi = f1 + (1.0 - y) / (x * y * y) * f2
    c_eff = sum(p * (3.0 * m * m - 2.0)
                for p, m in zip(cat.populations, (-1.0, 0.0, 1.0)))
    st2 = np.sin(cat.theta_s) ** 2
    hand = (-(1.0 - y) / (y * y) * c_eff * st2 / d_phi
            / float(cat.moments()[1]))
    np.testing.assert_allclose(resp.delta_response(cat), hand, rtol=1e-12)
    # and the per-event amplitude the response actually uses is that
    # derivative times the model, to float roundoff
    model = _truth_delta
    np.testing.assert_allclose(
        (resp.delta_response(cat) * resp.delta_cells(model))[resp.cell],
        resp.amplitude_per_event(cat), rtol=1e-12)


def test_fold_reproduces_the_reco_bin_amplitude_of_any_model(response):
    """The identity that makes `fold` a drop-in for the bin-by-bin K:
    folding the model the response was built from returns that bin's
    a_reco_bin, and folding a DIFFERENT model returns what an independent
    response built with THAT model measures -- the same events, the same
    smearing, only the structure function changed."""
    sampler, resp = response
    cat = bk.tensor_flip_plan(0.6).categories[0]
    model = _truth_delta
    for edges in ((0.01, 0.03, 1.0, 3.0), (0.03, 0.1, 2.0, 10.0),
                  (0.1, 0.3, 5.0, 50.0)):
        summ = resp.bin_summary(*edges, cat)
        mask = resp.mask_reco(*edges)
        assert resp.fold(model, mask, cat) == pytest.approx(
            summ["a_reco_bin"], rel=1e-12)
    # an independent response whose kernel carries 3.7 x the same Delta:
    # identical events (same seed, Delta does not enter the sampling of
    # the response), so its reco-bin amplitude must be exactly 3.7 x
    other = InclusiveKernel(
        beams.LI6, b1_func=toy_b1,
        delta_func=lambda x, q2, f1: 3.7 * toy_delta_gluon(x, q2, f1,
                                                           scale=3e-2))
    smp2 = InclusiveSampler(other, CONFIG, rp.generator_scenario(
        fom.Scenario(lumi_fb_per_nucleon=10.0, pol_ion_tensor=0.6)),
        nx=20, nq2=15, x_range=(3e-3, 0.5), q2_range=(0.7, 100.0))
    resp2 = rp.RecoResponse(smp2, rp.RecoModel(), n_mc_per_cell=300,
                            rng=np.random.default_rng(1))
    scaled = lambda x, q2, f1: 3.7 * model(x, q2, f1)  # noqa: E731
    for edges in ((0.01, 0.03, 1.0, 3.0), (0.03, 0.1, 2.0, 10.0)):
        assert resp.fold(scaled, resp.mask_reco(*edges), cat) == \
            pytest.approx(resp2.bin_summary(*edges, cat)["a_reco_bin"],
                          rel=1e-12)


def _stratified_bootstrap_spread(resp, model, mask, cat, nrep=600, seed=99):
    """sigma of `fold` from resampling the response IN EVERY SAMPLER CELL
    -- the same variance `fold_mc_error` computes by the delta method,
    derived instead by resampling, and written here from the definition
    of the response's sampling (`RecoResponse.__init__`: n_c independent
    draws in cell c, of which only some reconstruct into the bin) rather
    than from the estimator under test.

    The n_c draws of a cell are resampled with equal probabilities, the
    out-of-bin ones included -- that is the whole content of the
    stratification, and an estimator that resampled only the IN-BIN
    events would miss the dominant term wherever a bin takes a small
    fraction of a cell.
    """
    rng = np.random.default_rng(seed)
    ncell = resp.sampler.xsec_flat.size
    n_c = np.bincount(resp.cell, minlength=ncell)
    u = (resp.w * resp.eff)[mask]
    a = resp.delta_response(cat) * resp.delta_cells(model)
    z = u * resp.dil[mask] * a[resp.cell[mask]]
    cells = resp.cell[mask]
    num = np.zeros(nrep)
    den = np.zeros(nrep)
    for c in np.unique(cells):
        j = np.flatnonzero(cells == c)
        p = np.full(j.size + 1, 1.0 / n_c[c])
        p[-1] = 1.0 - j.size / n_c[c]           # the draws outside the bin
        m = rng.multinomial(n_c[c], p, size=nrep)[:, :j.size]
        den += m.dot(u[j])
        num += m.dot(z[j])
    return float((num / den).std(ddof=1))


def test_fold_mc_error_is_the_stratified_bootstrap_of_the_same_response(
        response):
    """The response MC error is the headline new term in the 7R bars
    (plans/08 A6), so pin its ALGEBRA, not just its order of magnitude:
    the delta-method variance must equal a cell-by-cell resampling of the
    same response to the bootstrap's own precision (600 replicas know a
    sigma to 2.9%).

    Two bins, one of which takes only ~5% of each cell it touches: there
    the stratification index matters, and an estimator that used the
    in-bin count instead of the cell's draw count n_c comes out 11x too
    small (0.087 of the bootstrap against 0.98 for the shipped one).
    """
    sampler, resp = response
    cat = bk.tensor_flip_plan(0.6).categories[0]
    model = _truth_delta
    for edges in ((0.03, 0.1, 2.0, 10.0), (0.03, 0.04, 2.0, 3.0)):
        mask = resp.mask_reco(*edges)
        boot = _stratified_bootstrap_spread(resp, model, mask, cat)
        est = resp.fold_mc_error(model, mask, cat)
        assert 0.90 < est / boot < 1.12, (edges, est, boot)


def test_fold_mc_error_matches_the_spread_over_response_seeds(response):
    """The same error against an OUTSIDE measurement: the scatter of
    `fold` over independent response samples.  16 seeds determine that
    scatter itself only to 18%, so this cannot be tight -- it is the
    check that the bootstrap above is bootstrapping the right thing, and
    the window still excludes the factor-two mis-normalizations
    (est/spread measured 1.07 and 1.12 on the two bins)."""
    sampler, resp = response
    cat = bk.tensor_flip_plan(0.6).categories[0]
    model = _truth_delta
    reps = [rp.RecoResponse(sampler, rp.RecoModel(), n_mc_per_cell=300,
                            rng=np.random.default_rng(700 + k))
            for k in range(16)]
    for edges in ((0.03, 0.1, 2.0, 10.0), (0.03, 0.04, 2.0, 3.0)):
        vals = np.array([r.fold(model, r.mask_reco(*edges), cat)
                         for r in reps])
        est = np.mean([r.fold_mc_error(model, r.mask_reco(*edges), cat)
                       for r in reps])
        spread = vals.std(ddof=1)
        assert 0.7 < est / spread < 1.5, (edges, est, spread)
        assert 1e-4 < est / abs(vals.mean()) < 0.05


def test_folded_fit_is_the_bin_by_bin_K_when_the_prior_is_the_truth(response):
    """Opt-in must not move anything: with the injected shape as prior the
    tilt sits at zero and K is the published bin-by-bin factor."""
    sampler, resp = response
    cat = bk.tensor_flip_plan(0.6).categories[0]
    model = _truth_delta
    bins = _slice_bins(resp, cat, 2.0, 10.0)
    assert len(bins) >= 5
    fit = rp.fold_shape_fit(resp, cat, bins, model)
    assert abs(fit["params"][1]) < 1e-3 and abs(fit["params"][2]) < 1e-3
    assert fit["chi2"] < 1e-12
    for b, k in zip(bins, fit["k"]):
        f1c = resp.f1_center(b["x"], b["q2"])
        assert k == pytest.approx(float(model(b["x"], b["q2"], f1c))
                                  / b["summ"]["a_reco_bin"], rel=2e-3)


def test_folded_fit_recovers_the_truth_from_a_wrong_prior(response):
    """The closure that matters (code review F5): noise-free pseudo-data
    generated with Delta = scale F1 x^0.3 (1-x)^4 (interpretation A),
    corrected with a prior that has no F1 factor (interpretation B) --
    the repository's own ansatz ambiguity.  The bin-by-bin K inherits the
    prior's shape error; the folded fit measures the shape from the data
    and must beat it by a large factor."""
    sampler, resp = response
    cat = bk.tensor_flip_plan(0.6).categories[0]
    model = _truth_delta
    prior = _prior_delta
    bins = _slice_bins(resp, cat, 2.0, 10.0)
    xs = np.array([b["x"] for b in bins])
    q2s = np.array([b["q2"] for b in bins])
    truth = np.asarray(model(xs, q2s, resp.f1_center(xs, q2s)), dtype=float)
    amps = np.array([b["amp"] for b in bins])

    k_bin = np.array([float(prior(b["x"], b["q2"], 0.0))
                      / resp.fold(prior, b["mask"], cat) for b in bins])
    bias_bin = amps * k_bin / truth - 1.0
    fit = rp.fold_shape_fit(resp, cat, bins, prior)
    bias_fit = amps * fit["k"] / truth - 1.0
    # Re-measured 2026-08-27 at the corrected beam energy (99.5 GeV/u,
    # gamma-matched, plans/10 -- the fixture ran at 50 GeV/u before, which
    # is not a machine configuration).  The prior's K is wrong by -22% at
    # the low-x edge and +16% at the high-x edge; the fitted shape by
    # -6.3% / +3.5%, and by 2.7% at the median.  The CLAIM is unchanged and
    # is what these assertions encode: the folded fit beats bin-by-bin by
    # ~3.5x on the worst bin and ~2.2x at the median.
    assert np.max(np.abs(bias_bin)) > 0.15          # the prior does hurt
    assert np.max(np.abs(bias_fit)) < 0.35 * np.max(np.abs(bias_bin))
    assert np.max(np.abs(bias_fit)) < 0.10
    assert np.median(np.abs(bias_fit)) < 0.55 * np.median(np.abs(bias_bin))
    # K is a ratio of the shape to its own fold, so the fitted
    # normalization cannot enter it
    big = rp.fold_shape_fit(resp, cat, bins,
                            lambda x, q2, f1: 137.0 * prior(x, q2, f1))
    np.testing.assert_allclose(big["k"] / fit["k"], 1.0, rtol=1e-6)


def test_prior_spread_covers_the_residual_the_tilt_cannot_absorb(response):
    """The folded fit removes most of the prior dependence but not all of
    it -- the tilt cannot reproduce the F1(x, Q2) factor that separates
    interpretations A and B -- so the 7R bar has to carry what is left.

    The claim under test is that refitting from the other shapes of the
    family and taking the spread of K is a bound on that residual: the
    bias is measured here against the INJECTED truth, which the fit never
    sees, and the spread is built from K alone."""
    sampler, resp = response
    cat = bk.tensor_flip_plan(0.6).categories[0]
    model, prior = _truth_delta, _prior_delta
    bins = _slice_bins(resp, cat, 2.0, 10.0)
    xs = np.array([b["x"] for b in bins])
    q2s = np.array([b["q2"] for b in bins])
    truth = np.asarray(model(xs, q2s, resp.f1_center(xs, q2s)), dtype=float)
    amps = np.array([b["amp"] for b in bins])

    plain = rp.fold_shape_fit(resp, cat, bins, prior)
    assert np.all(plain["k_prior_rel"] == 0.0)       # inert by default
    assert plain["alt_chi2"] == []
    np.testing.assert_allclose(plain["k_rel_all"], plain["k_rel"], rtol=1e-12)

    fit = rp.fold_shape_fit(resp, cat, bins, prior, alt_bases=(model,))
    np.testing.assert_array_equal(fit["k"], plain["k"])   # K itself unmoved
    bias = amps * fit["k"] / truth - 1.0
    assert np.max(np.abs(bias)) > 0.05               # there IS a residual
    assert np.all(fit["k_prior_rel"]
                  >= np.abs(bias) / (1.0 + np.abs(bias)) - 1e-9)
    assert np.all(fit["k_rel_all"] >= fit["k_rel"])
    # the alternative here IS the shape that generated the amplitudes, so
    # its own folded fit is the closure fit and the data prefer it
    assert len(fit["alt_chi2"]) == 1
    assert fit["alt_chi2"][0] < 1e-6 * fit["chi2"]


def test_folded_fit_does_not_oscillate(response):
    """The conditioning guard (plans/08 A6 (ii)): a matrix unfolding of
    this response oscillates to +-30% on noise-free pseudo-data because
    the migration matrix is ill-conditioned at purity ~0.6.  A bounded
    two-parameter tilt cannot: the residual must stay smooth in x, with
    no sign-alternating structure and no excursion beyond the bounds."""
    sampler, resp = response
    cat = bk.tensor_flip_plan(0.6).categories[0]
    model = _truth_delta
    prior = _prior_delta
    bins = _slice_bins(resp, cat, 2.0, 10.0)
    xs = np.array([b["x"] for b in bins])
    q2s = np.array([b["q2"] for b in bins])
    truth = np.asarray(model(xs, q2s, resp.f1_center(xs, q2s)), dtype=float)
    amps = np.array([b["amp"] for b in bins])
    fit = rp.fold_shape_fit(resp, cat, bins, prior)
    assert not fit["at_bound"]
    bias = amps * fit["k"] / truth - 1.0
    k_bin = np.array([float(prior(b["x"], b["q2"], 0.0))
                      / resp.fold(prior, b["mask"], cat) for b in bins])
    rough = np.max(np.abs(np.diff(amps * k_bin / truth - 1.0)))
    # measured 2026-08-27 at the corrected energy: 0.014 bin-to-bin against
    # a bin-by-bin roughness of 0.088, and 2 sign changes over 15 steps
    assert np.max(np.abs(np.diff(bias))) < 0.05     # smooth, bin to bin
    assert np.max(np.abs(np.diff(bias))) < 0.25 * rough
    signs = np.sign(np.diff(bias))
    assert np.sum(signs[1:] != signs[:-1]) <= 2     # monotone, not wobbling
    assert np.all(np.isfinite(fit["k_err"])) and np.all(fit["k_err"] >= 0.0)


def test_delta_from_amplitude_keeps_the_published_path(response):
    """The new keywords must be inert when unused, and additive when
    used: same K, same error, plus the K error in quadrature."""
    sampler, resp = response
    plan = bk.tensor_flip_plan(0.6)
    mask = resp.mask_reco(0.03, 0.1, 2.0, 10.0)
    summ = resp.bin_summary(0.03, 0.1, 2.0, 10.0, plan.categories[0])
    exact = rp.measure_inclusive(resp, plan, 1.0e4, mask, poisson=False)
    old = rp.delta_from_amplitude(exact, summ, delta_center=-0.05)
    assert old["delta"] == exact["amp"] * (-0.05 / summ["a_reco_bin"])
    assert old["err"] == exact["err"] * abs(-0.05 / summ["a_reco_bin"])
    same = rp.delta_from_amplitude(exact, summ, -0.05,
                                   k_conv=old["k"], k_rel_err=0.0)
    assert same["delta"] == old["delta"] and same["err"] == old["err"]
    with_k = rp.delta_from_amplitude(exact, summ, -0.05, k_conv=old["k"],
                                     k_rel_err=0.03)
    assert with_k["delta"] == old["delta"]
    assert with_k["err"] == pytest.approx(
        np.hypot(old["err"], 0.03 * abs(with_k["delta"])), rel=1e-12)
    assert with_k["err"] > old["err"]


def test_shape_fit_error_is_the_delta_chi2_band(response):
    """The K error the folded fit propagates comes from the curvature of
    the profiled chi2; check it against the textbook definition it stands
    for -- the spread of K over every shape within Delta chi2 <= 1,
    scanned by brute force."""
    sampler, resp = response
    cat = bk.tensor_flip_plan(0.6).categories[0]
    prior = _prior_delta
    bins = _slice_bins(resp, cat, 2.0, 10.0)
    fit = rp.fold_shape_fit(resp, cat, bins, prior, mc_error=False)
    n0, c1, c2 = fit["params"]
    amp = np.array([b["amp"] for b in bins])
    wgt = 1.0 / np.array([b["err"] for b in bins]) ** 2
    kers = [resp.fold_kernel(b["mask"], cat) for b in bins]
    xs = np.array([b["x"] for b in bins])
    q2s = np.array([b["q2"] for b in bins])
    base_c = np.asarray(prior(xs, q2s, 0.0), dtype=float)
    lo, hi = [], []
    sig = np.sqrt(np.diag(fit["cov"]))       # scan window only; the band
    for d1 in np.linspace(-2.5 * sig[0], 2.5 * sig[0], 25):   # itself comes
        for d2 in np.linspace(-2.5 * sig[1], 2.5 * sig[1], 25):  # from chi2
            shp = rp.tilted_shape(prior, fit["x_ref"], (1.0, c1 + d1, c2 + d2))
            fold = np.array([resp.fold(shp, b["mask"], cat, kernel=k)
                             for b, k in zip(bins, kers)])
            norm = float((amp * fold * wgt).sum() / (fold * fold * wgt).sum())
            chi2 = float((wgt * (amp - norm * fold) ** 2).sum())
            if chi2 - fit["chi2"] <= 1.0:
                k = (base_c * np.exp((c1 + d1) * np.log(xs / fit["x_ref"])
                                     + (c2 + d2) * np.log(1 - xs)) / fold)
                lo.append(k)
                hi.append(k)
    band = 0.5 * (np.max(hi, axis=0) - np.min(lo, axis=0)) / np.abs(fit["k"])
    assert len(lo) > 20
    np.testing.assert_allclose(fit["k_fit_rel"], band, rtol=0.25, atol=1e-4)


# --- the two systematics this programme specifies for itself (report 3 §7) ---

def test_pzz_scale_error_propagates_one_to_one_not_quadratically():
    """CORRECTED 2026-08-27.  An earlier version of this test asserted only
    that sigma_P^2 scales as (1+d)^2 and concluded the amplitude did too.
    That is true of sigma_P^2 and IRRELEVANT to the amplitude, because it
    never ran the estimator on rescaled data.

    The estimator's weights w_f = P_f - Pbar are built from the ASSUMED
    polarizations, so the ratio R carries one power of the assumed scale
    while sigma_P^2 carries two.  One power cancels:

        A_hat / A = P_zz(true) / P_zz(assumed),  EXACTLY.

    So delta P_zz / P_zz costs the SAME fraction on the amplitude, not
    twice it.  The quadratic is real but belongs to REACH rather than bias:
    delta_A ~ 1/sigma_P, so the luminosity needed at fixed delta_A goes as
    1/P_zz^2."""
    sc = coh.CoherentScenario(amp=0.01, eps_b0=-0.08)
    cfg = beams.default_configs("6Li")[0]
    cr = rp.CoherentResponse(sc, cfg, reco.sigma_theta_tagging(cfg),
                             shape="rectangle", n_mc=120000,
                             cut_scale_xy=(1.0, 1.0),
                             rng=np.random.default_rng(11))
    plan = bk.tensor_flip_plan(0.6)
    pzz = [float(c.moments()[1]) for c in plan.categories]
    frac = [c.lumi_fraction for c in plan.categories]
    a_t = lambda t: sc.cos2phi_coefficient_deformation(t, 1.0)   # noqa: E731
    ae = np.linspace(0.0, 2.0 * np.pi, 13)
    be = np.linspace(0.0, 2.0 * np.pi, 25)

    def fit(pzz_true):
        mu = np.concatenate(
            [cr.expected_counts_2d(1e9, [pf], [lf], 0.05, 0.12, ae, be,
                                   0.01, a_t, a_m=0.0, u1=0.05, u2=0.02)
             for pf, lf in zip(pzz_true, frac)], axis=0)
        f = reco.harmonic_ratio_fit_2d(mu, list(frac), list(pzz), ae, be,
                                       u_coeffs=(0.05, 0.02),
                                       beta_means=cr.basis_means(0.05, 0.12, be),
                                       with_sin=False)
        return f["a_e"], f["a_t"]

    base_e, base_t = fit(pzz)
    assert base_e == pytest.approx(0.01, rel=1e-6)      # closes with no error
    for d in (0.02, 0.05, 0.10, 0.30):
        e, t = fit([p * (1.0 + d) for p in pzz])
        assert e / base_e == pytest.approx(1.0 + d, rel=1e-6)
        assert t / base_t == pytest.approx(1.0 + d, rel=1e-6)


def test_the_flip_plan_luminosity_split_is_optimal():
    """sigma_P^2 = f(1-f)(P1-P2)^2 for two states, so an equal split maximises
    it and the plan is already there.  Worth pinning: an unequal split would
    cost statistical precision for nothing."""
    plan = bk.tensor_flip_plan(0.6)
    pzz = [float(c.moments()[1]) for c in plan.categories]
    assert plan.categories[0].lumi_fraction == pytest.approx(0.5)

    def sigma_p2(f):
        lum = np.array([f, 1.0 - f])
        pbar = float((lum * np.asarray(pzz)).sum())
        return float((lum * (np.asarray(pzz) - pbar) ** 2).sum())

    best = max(np.linspace(0.05, 0.95, 91), key=sigma_p2)
    assert best == pytest.approx(0.5, abs=0.02)
    # and the closed form
    assert sigma_p2(0.5) == pytest.approx(0.25 * (pzz[0] - pzz[1]) ** 2)


def test_relative_luminosity_bias_is_one_third_per_unit_ratio_error():
    """CORRECTED 2026-08-27.  An earlier version quoted "~1.4 x delta" and
    asserted only a loose 0.01-0.20 band, which was too wide to notice that
    the number was wrong and which conflated two conventions.

    To first order the bias is

        dA/A = -[(P1 + P2) / (P1 - P2)] x delta_ratio

    where delta_ratio is the fractional error on L1/L2.  For the flip plan's
    (+0.6, -1.2) that coefficient is exactly 1/3.

    MIND THE CONVENTION.  The scripts' --rel-lumi-offset d sets the ASSUMED
    shares to [0.5(1+d), 0.5(1-d)] against equal truth, which is a RATIO
    error of 2d/(1-d^2) ~ 2d -- so the apparent coefficient in that
    convention is ~2/3, not 1/3.  Quote which one you mean."""
    sc = coh.CoherentScenario(amp=0.01, eps_b0=-0.08)
    cfg = beams.default_configs("6Li")[0]
    cr = rp.CoherentResponse(sc, cfg, reco.sigma_theta_tagging(cfg),
                             shape="rectangle", n_mc=120000,
                             cut_scale_xy=(1.0, 1.0),
                             rng=np.random.default_rng(11))
    plan = bk.tensor_flip_plan(0.6)
    pzz = [float(c.moments()[1]) for c in plan.categories]
    a_t = lambda t: sc.cos2phi_coefficient_deformation(t, 1.0)   # noqa: E731
    ae = np.linspace(0.0, 2.0 * np.pi, 13)
    be = np.linspace(0.0, 2.0 * np.pi, 25)

    def fit(lum_assumed):
        mu = np.concatenate(
            [cr.expected_counts_2d(1e9, [pf], [0.5], 0.05, 0.12, ae, be,
                                   0.01, a_t, a_m=0.0, u1=0.05, u2=0.02)
             for pf in pzz], axis=0)
        f = reco.harmonic_ratio_fit_2d(mu, list(lum_assumed), list(pzz),
                                       ae, be, u_coeffs=(0.05, 0.02),
                                       beta_means=cr.basis_means(0.05, 0.12, be),
                                       with_sin=False)
        return f["a_t"]

    analytic = -(pzz[0] + pzz[1]) / (pzz[0] - pzz[1])
    assert analytic == pytest.approx(1.0 / 3.0, rel=1e-9)

    base = fit([0.5, 0.5])
    for d in (1e-3, 1e-2):
        t = fit([0.5 * (1.0 + d), 0.5 * (1.0 - d)])
        ratio_err = (1.0 + d) / (1.0 - d) - 1.0
        assert (t / base - 1.0) / ratio_err == pytest.approx(analytic, rel=0.05)
    # and in absolute terms it is small: a per-cent ratio error costs ~0.3%
    t = fit([0.5 * 1.005, 0.5 * 0.995])
    assert abs(t / base - 1.0) < 0.005
