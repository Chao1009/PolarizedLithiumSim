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


# --- what the coherent ratio does NOT cancel (code review F1, F6) -----------

@pytest.fixture(scope="module")
def slot_response():
    """The 6R default: slot-like cutout, isotropic divergence."""
    sc = coh.CoherentScenario(amp=0.01, eps_b0=-0.08)
    return sc, rp.CoherentResponse(sc, CONFIG, reco.SIGMA_THETA_HA,
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
    plain = rp.CoherentResponse(sc, CONFIG, reco.SIGMA_THETA_HA, **kw)
    lifted = rp.CoherentResponse(sc, CONFIG, reco.SIGMA_THETA_HA,
                                 t_floor=0.25 * plain.pt_cut ** 2, **kw)
    assert lifted.acceptance == pytest.approx(plain.acceptance, rel=0.05)
    assert lifted.t_true.size > 1.5 * plain.t_true.size
    assert lifted._t_weight < 1.0 and plain._t_weight == 1.0
    # the analytic value for the circular cut, within the smearing effect
    analytic = float(sc.tag_acceptance_angular(
        reco.SIGMA_THETA_HA, CONFIG.ion_momentum_per_nucleon))
    assert 0.8 < lifted.acceptance / analytic < 1.5

    top = beams.default_configs("6Li")[2]          # 0.60 GeV envelope
    dead = rp.CoherentResponse(sc, top, reco.SIGMA_THETA_HA, **kw)
    assert dead.t_true.size == 0
    alive = rp.CoherentResponse(sc, top, reco.SIGMA_THETA_HA,
                                t_floor=0.30, **kw)
    assert alive.t_true.size > 1000
    assert 0.0 < alive.acceptance < 1e-6


def test_t_floor_weights_leave_the_fit_unbiased():
    """The per-event weight is uniform, so nothing downstream of the
    sampler needs to know about it."""
    sc = coh.CoherentScenario(amp=0.01, eps_b0=-0.08)
    plan = bk.tensor_flip_plan(0.6)
    a_t = lambda t: sc.cos2phi_coefficient_deformation(t, 1.0)  # noqa: E731
    cr = rp.CoherentResponse(sc, CONFIG, reco.SIGMA_THETA_HA,
                             cut_scale_xy=(2.5, 1.0), n_mc=200000,
                             rng=np.random.default_rng(11), t_floor=0.012)
    fit = rp.measure_coherent(cr, 8e7, plan, 0.05, 0.08, a_e=0.010,
                              a_t_func=a_t, u1=0.05, u2=0.02, poisson=False)
    assert fit["a_t"] == pytest.approx(fit["truth"]["a_t"], rel=1e-4)
    assert fit["a_e"] == pytest.approx(0.010, rel=1e-4)
