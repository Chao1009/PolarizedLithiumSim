"""Measured quantities and reconstruction (polligen.reco): frames,
covariant azimuths, electron-method resolution, spin-state-sorted
harmonic estimator, coherent recoil + Roman-Pot emulation."""

import pathlib
import sys

import math

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen import reco  # noqa: E402
from polligen.estimators import cos2phi_fit_binned  # noqa: E402

from polli_fastsim import beams  # noqa: E402
from polli_fastsim.kinematics import scattered_electron  # noqa: E402

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
S_NN = CONFIG.sqrt_s_per_nucleon ** 2


# --- frames ------------------------------------------------------------------

def test_lab_frame_beam_directions_and_roundtrip():
    k, p = reco.beam_fourvectors(CONFIG)
    kl, pl = reco.head_on_to_lab(k), reco.head_on_to_lab(p)
    e_dir = kl[1:] / np.linalg.norm(kl[1:])
    i_dir = pl[1:] / np.linalg.norm(pl[1:])
    np.testing.assert_allclose(e_dir, [0.0, 0.0, -1.0], atol=1e-12)
    # massive ion: direction agrees with the massless expectation to O(M^2/E^2)
    np.testing.assert_allclose(i_dir, [-np.sin(reco.XING_IP6), 0.0,
                                       np.cos(reco.XING_IP6)], atol=1e-5)
    np.testing.assert_allclose(reco.lab_to_head_on(kl), k, atol=1e-9)
    np.testing.assert_allclose(reco.lab_to_head_on(pl), p, atol=1e-9)


def test_lab_azimuth_about_detector_axis_only_odd_harmonics():
    """Uniform phi about the ion axis (head-on frame) is not exactly
    uniform about the detector axis, but because the ePIC axis coincides
    with the electron beam the boost and rotation of the head-on
    transformation cancel at first order for the scattered electron: the
    residual is a ~1e-3 cos phi distortion (odd harmonics only), while
    the cos 2phi harmonic -- the observable -- is untouched.  The
    transformation is applied anyway (standard and exact)."""
    phi = np.linspace(0.0, 2.0 * np.pi, 4000, endpoint=False)
    for x, y, a1_min in ((0.056, 0.0101, 3e-4), (0.1413, 0.0506, 1e-3)):
        kp = reco.electron_fourvector(x, y, S_NN, CONFIG.electron_energy, phi)
        phi_lab = reco.azimuth_about_z(reco.head_on_to_lab(kp))
        a1 = 2.0 * np.mean(np.cos(phi_lab))
        a2 = 2.0 * np.mean(np.cos(2.0 * phi_lab))
        resp = 2.0 * np.mean(np.cos(2.0 * phi_lab) * np.cos(2.0 * phi))
        assert a1_min < abs(a1) < 1e-2
        assert abs(a2) < 1e-12
        assert abs(resp - 1.0) < 1e-4
        back = reco.azimuth_about_z(reco.lab_to_head_on(reco.head_on_to_lab(kp)))
        assert np.abs(np.angle(np.exp(1j * (back - phi)))).max() < 1e-10


# --- covariant azimuths ------------------------------------------------------

def test_phi_s_equals_lab_angle_for_massless_target():
    k, p = reco.beam_fourvectors(CONFIG, ion_mass=0.0)
    rng = np.random.default_rng(1)
    n = 500
    x = 10 ** rng.uniform(-3, -0.7, n)
    y = rng.uniform(0.01, 0.95, n)
    phi_e = rng.uniform(0, 2 * np.pi, n)
    phi_s = rng.uniform(0, 2 * np.pi)
    kp = reco.electron_fourvector(x, y, S_NN, CONFIG.electron_energy, phi_e)
    phi = reco.azimuth_wrt_lepton_plane(k, kp, p, reco.spin_fourvector(phi_s))
    dev = np.angle(np.exp(1j * (phi - (phi_e - phi_s))))
    assert np.abs(dev).max() < 1e-12


def test_phi_s_massive_target_deviation_is_order_gamma2():
    k, p = reco.beam_fourvectors(CONFIG)
    rng = np.random.default_rng(2)
    n = 2000
    x = 10 ** rng.uniform(-3, -0.7, n)
    y = rng.uniform(0.01, 0.95, n)
    phi_e = rng.uniform(0, 2 * np.pi, n)
    kp = reco.electron_fourvector(x, y, S_NN, CONFIG.electron_energy, phi_e)
    phi = reco.azimuth_wrt_lepton_plane(k, kp, p,
                                        reco.spin_fourvector(np.pi / 2))
    dev = np.abs(np.angle(np.exp(1j * (phi - (phi_e - np.pi / 2)))))
    g2 = 4.0 * 0.9383 ** 2 * x * x / (S_NN * x * y)
    assert dev.max() < 0.5 * g2.max()          # bounded by gamma^2 / 2
    assert dev[g2 < 1e-3].max() < 1e-3


def _boost_to_rest(v, w):
    beta = w[1:] / w[0]
    b2 = beta @ beta
    g = 1.0 / np.sqrt(1.0 - b2)
    e, pv = v[0], v[1:]
    bp = pv @ beta
    return np.concatenate([[g * (e - bp)],
                           pv + ((g - 1.0) * bp / b2 - g * e) * beta])


def _rot_to_z(vec):
    a = vec / np.linalg.norm(vec)
    z = np.array([0.0, 0.0, 1.0])
    c = a @ z
    ax = np.cross(a, z)
    s = np.linalg.norm(ax)
    if s < 1e-15:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    kx = ax / s
    kmat = np.array([[0, -kx[2], kx[1]], [kx[2], 0, -kx[0]],
                     [-kx[1], kx[0], 0]])
    return np.eye(3) + s * kmat + (1.0 - c) * kmat @ kmat


def _explicit_collinear_azimuth(k, kp, p, v):
    """Independent construction: boost to the gamma*-A c.m. frame, rotate
    the photon onto +z (Trento orientation), read azimuths about z."""
    q = k - kp
    w = p + q
    kk, pp, qq, vv = (_boost_to_rest(u, w) for u in (k, p, q, v))
    rot = _rot_to_z(qq[1:])
    kk, pp, vv = (rot @ u[1:] for u in (kk, pp, vv))
    assert np.linalg.norm(pp[:2]) < 1e-9 * np.linalg.norm(pp)
    return np.angle(np.exp(1j * (np.arctan2(vv[1], vv[0])
                                 - np.arctan2(kk[1], kk[0]))))


def test_covariant_azimuth_matches_explicit_collinear_frame():
    k, p = reco.beam_fourvectors(CONFIG)
    rng = np.random.default_rng(5)
    for _ in range(60):
        x = 10 ** rng.uniform(-3, -0.7)
        y = rng.uniform(0.02, 0.9)
        kp = reco.electron_fourvector(x, y, S_NN, CONFIG.electron_energy,
                                      rng.uniform(0, 2 * np.pi))
        s_vec = reco.spin_fourvector(rng.uniform(0, 2 * np.pi),
                                     rng.uniform(1.2, 1.9))
        a = reco.azimuth_wrt_lepton_plane(k, kp, p, s_vec)
        b = _explicit_collinear_azimuth(k, kp, p, s_vec)
        assert abs(np.angle(np.exp(1j * (a - b)))) < 1e-10
        if x < 0.02:
            pp = reco.recoil_fourvector(-rng.uniform(0.04, 0.2),
                                        rng.uniform(0, 2 * np.pi),
                                        rng.uniform(0.001, 0.01), p)
            a = reco.azimuth_wrt_lepton_plane(k, kp, p, pp)
            b = _explicit_collinear_azimuth(k, kp, p, pp)
            assert abs(np.angle(np.exp(1j * (a - b)))) < 1e-9


def test_recoil_lab_shortcut_is_sub_mrad():
    k, p = reco.beam_fourvectors(CONFIG)
    kp = reco.electron_fourvector(0.002, 0.9, S_NN, CONFIG.electron_energy,
                                  np.linspace(0, 2 * np.pi, 13))
    pp = reco.recoil_fourvector(-0.04, 0.3, 0.01, p)
    err = reco.lab_azimuth_shortcut_error(k, kp, p, pp, 0.3)
    assert np.abs(err).max() < 1e-3


# --- electron reconstruction -------------------------------------------------

def test_electron_method_roundtrip():
    rng = np.random.default_rng(3)
    x = 10 ** rng.uniform(-3, -0.5, 300)
    y = rng.uniform(0.01, 0.95, 300)
    e_p, th, _ = scattered_electron(x, y, S_NN, CONFIG.electron_energy)
    q2, yr, xr = reco.electron_method(e_p, th, CONFIG.electron_energy, S_NN)
    np.testing.assert_allclose(q2, S_NN * x * y, rtol=1e-10)
    np.testing.assert_allclose(yr, y, rtol=1e-10)
    np.testing.assert_allclose(xr, x, rtol=1e-9)


def test_electron_method_resolution_formula_vs_mc():
    """At the y = 0.025 sweet spot the electron method loses y: the
    linear formula reproduces the MC spread and dy/y ~ 50% for a 1.2%
    energy resolution."""
    x, y = 0.0224, 0.0254
    e_p, th, _ = scattered_electron(x, y, S_NN, CONFIG.electron_energy)
    rng = np.random.default_rng(4)
    n = 200000
    e_s, th_s, _ = reco.smear_electron(np.full(n, e_p), np.full(n, th),
                                       np.zeros(n), 0.012, 1e-3, 1e-3, rng)
    _, y_r, x_r = reco.electron_method(e_s, th_s, CONFIG.electron_energy, S_NN)
    dq2, dy, dx = reco.electron_method_resolution(y, np.pi - th, 0.012, 1e-3)
    assert 0.4 < dy < 0.6
    assert np.std(y_r) / y == pytest.approx(dy, rel=0.1)
    ok = np.isfinite(x_r) & (y_r > 0.3 * y)   # linear regime only
    assert np.std(x_r[ok]) / x > 0.3          # x unusable from e' alone


def test_mixed_method_x_resolution_tracks_hadronic_y():
    rng = np.random.default_rng(6)
    y = np.full(100000, 0.0254)
    y_h = reco.hadronic_y(y, 0.2, rng)
    x_r = reco.mixed_method(S_NN * 0.0224 * 0.0254, y_h, S_NN)
    assert np.std(np.log(x_r)) == pytest.approx(0.2, rel=0.1)


# --- spin-state-sorted estimator ---------------------------------------------

def _acceptance(phi):
    return 1.0 + 0.03 * np.cos(2.0 * phi) + 0.02 * np.cos(phi)


def test_ratio_estimator_unbiased_under_phi_dependent_efficiency():
    edges = np.linspace(0.0, 2.0 * np.pi, 25)
    amp, pz, n_tot = 0.012, [0.6, -1.2], 4e6
    mu, frac = reco.expected_counts_by_fill(n_tot, pz, amp, edges,
                                            acceptance=_acceptance,
                                            const=0.01)
    rng = np.random.default_rng(7)
    est, err = [], []
    for _ in range(200):
        out = reco.harmonic_ratio_fit(rng.poisson(mu), frac, pz, edges)
        est.append(out["amp"])
        err.append(out["err"])
    est = np.asarray(est)
    analytic = reco.err_harmonic_ratio(n_tot, pz)
    assert est.mean() == pytest.approx(amp, abs=4 * est.std() / np.sqrt(est.size))
    assert est.std() == pytest.approx(analytic, rel=0.25)
    assert np.mean(err) == pytest.approx(analytic, rel=0.1)
    # two fills (0.6, -1.2) beat the single-fill sqrt(2/N)/0.6 by 1.5x
    w = np.pi / 24
    assert analytic * 0.6 / np.sqrt(2.0 / n_tot) == pytest.approx(
        2.0 * 0.6 / 1.8 / (np.sin(2 * w) / (2 * w)), rel=1e-6)
    # the single-fill fit under the same efficiency is biased by ~0.03/0.6
    mu1, _ = reco.expected_counts_by_fill(n_tot, [0.6], amp, edges,
                                          acceptance=_acceptance)
    single = np.mean([cos2phi_fit_binned(rng.poisson(mu1[0]), edges, 0.6)
                      for _ in range(50)])
    assert single - amp == pytest.approx(0.03 / 0.6, rel=0.1)


def test_ratio_estimator_immune_to_relative_luminosity_offset():
    edges = np.linspace(0.0, 2.0 * np.pi, 25)
    amp, pz = 0.012, [0.6, -1.2]
    mu, frac = reco.expected_counts_by_fill(1e9, pz, amp, edges,
                                            acceptance=_acceptance)
    # a 2% luminosity error on fill 0 that the analysis does not know about
    wrong = np.asarray(frac) * np.array([1.0, 1.02])
    out = reco.harmonic_ratio_fit(mu, wrong, pz, edges)
    # a 2% relative-luminosity error (200x the plans/05 reference 1e-4)
    # rescales the amplitude only at second order (~0.7%) and shifts the
    # constant; at 1e-4 the amplitude effect is 3e-5 relative
    assert out["amp"] == pytest.approx(amp, rel=1e-2)
    assert abs(out["const"]) > 1e-3       # the offset lands in the constant


# --- coherent recoil + Roman Pots --------------------------------------------

def test_recoil_fourvector_exact_t_and_mass():
    _, p = reco.beam_fourvectors(CONFIG)
    t = -np.array([0.02, 0.05, 0.1, 0.2])
    pp = reco.recoil_fourvector(t, 0.7, 0.005, p)
    np.testing.assert_allclose(reco.t_from_fourvectors(p, pp), t, rtol=1e-9)
    np.testing.assert_allclose(reco.mdot(pp, pp), reco.mdot(p, p), rtol=1e-9)
    with pytest.raises(ValueError):
        reco.recoil_fourvector(-1e-4, 0.0, 0.01, p)   # below t_min


def test_tag_pt_cut_scaling():
    """The LEGACY proton-derived constant, kept because every number
    published before 2026-08-27 used it (plans/10)."""
    assert reco.tag_pt_cut(reco.SIGMA_THETA_HA, 275.0, a_beam=1) == pytest.approx(0.20)
    assert reco.tag_pt_cut(reco.SIGMA_THETA_HA, 137.5, a_beam=6) == pytest.approx(0.60)
    assert reco.tag_pt_cut(reco.SIGMA_THETA_HA, 50.0, a_beam=6) == pytest.approx(0.218,
                                                                                 rel=1e-2)


def test_sigma_theta_for_reproduces_the_yellow_report_for_protons():
    """plans/10: the per-configuration divergence.  A PROTON must come back
    with Table 10.1's own numbers unchanged -- that is the check that the
    species step is not being applied where it should not be."""
    from polli_fastsim import beams as _b
    lo, mid, top = _b.default_configs("p")
    assert [round(1e6 * x) for x in reco.sigma_theta_for(top)] == [65, 65]
    assert [round(1e6 * x) for x in reco.sigma_theta_for(mid)] == [180, 180]
    assert [round(1e6 * x) for x in reco.sigma_theta_for(lo)] == [220, 380]
    assert [round(1e6 * x) for x in
            reco.sigma_theta_for(top, "high-divergence")] == [119, 119]


def test_sigma_theta_species_step_only_applies_where_rigidity_binds():
    """6Li is GAMMA-matched at the two lower configurations -- same speed as
    the proton, so the same divergence -- and only rigidity-capped at the
    top, where it picks up sqrt(2)."""
    from polli_fastsim import beams as _b
    for i, expect in ((0, 1.00), (1, 1.00), (2, 2.0 ** 0.5)):
        cfg = _b.default_configs("6Li")[i]
        p_cfg = _b.default_configs("p")[i]
        f = reco.sigma_theta_for(cfg)[0] / reco.sigma_theta_for(p_cfg)[0]
        assert f == pytest.approx(expect, rel=0.02)


def test_tagging_optics_is_the_analytic_optimum():
    """sigma_theta_tagging puts the n-sigma cut at t = 1/B, which is where
    L x acceptance peaks -- so the acceptance comes out at 1/e regardless of
    configuration.  That invariance IS the derivation."""
    from polli_fastsim import beams as _b
    from polligen import coherent as _c
    sc = _c.CoherentScenario()
    for cfg in _b.default_configs("6Li"):
        st = reco.sigma_theta_tagging(cfg, slope_b=sc.slope_b)
        assert reco.tag_pt_cut(st, cfg.ion_momentum_per_nucleon,
                               a_beam=cfg.ion.A) == pytest.approx(
                                   sc.slope_b ** -0.5, rel=1e-9)
        acc = sc.tag_acceptance_angular(st, cfg.ion_momentum_per_nucleon)
        assert acc == pytest.approx(math.exp(-1.0), rel=0.02)


def test_rp_hole_acceptance_harmonics():
    sq = reco.rp_hole_acceptance(50.0, 0.2, 0.2)
    assert abs(sq["a2"]) < 1e-9 and sq["a4"] > 0.2
    el = reco.rp_hole_acceptance(50.0, 0.2, 0.2, shape="ellipse")
    assert el["acc"] == pytest.approx(np.exp(-50.0 * 0.04), rel=1e-6)
    assert abs(el["a2"]) < 1e-9 and abs(el["a4"]) < 1e-9
    tall = reco.rp_hole_acceptance(50.0, 0.2, 0.3)
    wide = reco.rp_hole_acceptance(50.0, 0.3, 0.2)
    assert tall["a2"] > 0.3 and wide["a2"] == pytest.approx(-tall["a2"], rel=1e-9)


def test_rp_measure_matches_analytic_acceptance_and_smears_phi():
    _, p = reco.beam_fourvectors(CONFIG)
    rng = np.random.default_rng(8)
    n = 200000
    t = -rng.exponential(1.0 / 50.0, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    pp = reco.recoil_fourvector(t, phi, 0.0, p)
    sig = SIG_TAG
    m = reco.rp_measure(pp, p, (sig, sig), shape="ellipse", rng=rng)
    cut = reco.tag_pt_cut(sig, CONFIG.ion_momentum_per_nucleon)
    # the cut applies to the divergence-smeared pT: Gaussian in each
    # component with variance 1/(2B) + (p_A sigma)^2  ->  B_eff
    p_a = 6.0 * CONFIG.ion_momentum_per_nucleon
    b_eff = 1.0 / (1.0 / 50.0 + 2.0 * (p_a * sig) ** 2)
    assert m["accepted"].mean() == pytest.approx(np.exp(-b_eff * cut * cut),
                                                 rel=0.03)
    acc = m["accepted"]
    dphi = np.angle(np.exp(1j * (m["phi_t"][acc] - phi[acc])))
    assert np.std(dphi) == pytest.approx(p_a * sig / np.mean(m["pT"][acc]),
                                         rel=0.3)
    np.testing.assert_allclose(m["t_reco"][acc], -m["pT"][acc] ** 2)


# --- fill-dependent phi' acceptance (code review F1) ----------------------

PZZ_FLIP = (0.6, -1.2)          # m = +-1-rich / m = 0-rich pattern


def _harmonic_eff(e2, e1=0.0):
    return lambda phi: 1.0 + e2 * np.cos(2.0 * phi) + e1 * np.cos(phi)


def _fit_two_fills(amp, eps2_per_fill, n_total=1e9, nbins=24):
    """Noise-free expected counts under a per-fill efficiency, refitted."""
    edges = np.linspace(0.0, 2.0 * np.pi, nbins + 1)
    counts, frac = reco.expected_counts_by_fill(
        n_total, PZZ_FLIP, amp, edges,
        acceptance=[_harmonic_eff(e) for e in eps2_per_fill])
    return reco.harmonic_ratio_fit(counts, frac, PZZ_FLIP, edges,
                                   with_sin=True)


def test_common_phi_efficiency_cancels_exactly():
    """The estimator's whole point: a COMMON eps(phi') drops out."""
    fit = _fit_two_fills(1e-3, (0.03, 0.03))
    assert fit["amp"] == pytest.approx(1e-3, abs=1e-8)
    assert fit["amp_sin"] == pytest.approx(0.0, abs=1e-8)


def test_fill_dependent_efficiency_bias_matches_analytic():
    """dA = sum_f l_f (P_f - Pbar) e_f / sigma_P^2; for two equal-luminosity
    fills that is (e_+ - e_0)/(P_+ - P_0) -- the code review's F1 table."""
    for eplus, ezero in ((0.03, 0.02), (0.03, 0.029), (0.03, 0.0299)):
        predicted = reco.fill_acceptance_bias((eplus, ezero), PZZ_FLIP)
        assert predicted == pytest.approx(
            (eplus - ezero) / (PZZ_FLIP[0] - PZZ_FLIP[1]), rel=1e-12)
        fit = _fit_two_fills(1e-3, (eplus, ezero))
        assert fit["amp"] - 1e-3 == pytest.approx(predicted, rel=2e-3)
    # the tabulated case: 0.03 vs 0.02 fakes 5.6x a Delta/F1 ~ 1e-3 signal
    assert reco.fill_acceptance_bias((0.03, 0.02), PZZ_FLIP) == \
        pytest.approx(5.5556e-3, rel=1e-4)


def test_fill_acceptance_bias_is_zero_for_a_common_harmonic():
    assert reco.fill_acceptance_bias((0.03, 0.03), PZZ_FLIP) == 0.0
    assert reco.fill_acceptance_bias((0.03, 0.03, 0.03),
                                     (0.6, -1.2, 0.0)) == pytest.approx(0.0)


def test_fill_acceptance_bias_weights_by_luminosity():
    """For TWO fills the bias is (e_+ - e_0)/(P_+ - P_0) at any luminosity
    split -- the l_1 l_2 (P_1 - P_2) factors cancel between numerator and
    sigma_P^2 -- so a lopsided run plan buys no protection.  With three
    fills the luminosity weighting is real."""
    lf = (0.7, 0.3)
    edges = np.linspace(0.0, 2.0 * np.pi, 25)
    counts, frac = reco.expected_counts_by_fill(
        1e9, PZZ_FLIP, 1e-3, edges, lumi_fractions=lf,
        acceptance=[_harmonic_eff(0.03), _harmonic_eff(0.02)])
    fit = reco.harmonic_ratio_fit(counts, frac, PZZ_FLIP, edges)
    predicted = reco.fill_acceptance_bias((0.03, 0.02), PZZ_FLIP,
                                          lumi_fractions=lf)
    assert predicted == pytest.approx(0.01 / 1.8, rel=1e-12)
    assert fit["amp"] - 1e-3 == pytest.approx(predicted, rel=3e-3)

    three = (0.6, -1.2, 0.0)
    e3 = (0.03, 0.02, 0.03)
    assert reco.fill_acceptance_bias(e3, three, (1, 1, 1)) != \
        pytest.approx(reco.fill_acceptance_bias(e3, three, (1, 1, 4)),
                      rel=1e-3)


def test_per_fill_acceptance_normalizes_and_validates():
    assert reco.per_fill_acceptance(None, 3) == [None, None, None]
    f = _harmonic_eff(0.01)
    assert reco.per_fill_acceptance(f, 2) == [f, f]
    with pytest.raises(ValueError):
        reco.per_fill_acceptance([f, f, f], 2)


# --- eta-dependent calorimeter resolution (code review F4) ------------------

def test_emcal_eta_table_matches_the_yellow_report_regions():
    """Without `eta` the backward-endcap (PbWO4) specification applies
    everywhere -- the behaviour the chain had before 2026-08-25 and the
    default, so nothing published moves."""
    for e in (1.0, 10.0):
        assert reco.emcal_resolution(e) == reco.emcal_resolution(e, eta=-3.0)
    # backward endcap unchanged, barrel and transition coarser
    assert reco.emcal_resolution(10.0, eta=-2.5) == pytest.approx(0.0118,
                                                                 abs=1e-4)
    assert reco.emcal_resolution(10.0, eta=-1.6) / \
        reco.emcal_resolution(10.0, eta=-3.0) == pytest.approx(2.06, rel=0.05)
    assert reco.emcal_resolution(10.0, eta=0.0) / \
        reco.emcal_resolution(10.0, eta=-3.0) == pytest.approx(2.81, rel=0.05)
    # the ratio is energy-dependent: quoting "x3" without an energy is wrong
    assert reco.emcal_resolution(1.0, eta=0.0) / \
        reco.emcal_resolution(1.0, eta=-3.0) > 4.0
    # vectorized over eta
    out = reco.emcal_resolution(np.full(4, 10.0),
                                eta=np.array([-3.0, -1.6, 0.0, 2.0]))
    assert out.shape == (4,) and np.all(np.diff(out[:3]) > 0)


def test_best_energy_switches_to_the_tracker_in_the_barrel():
    """With the Yellow Report table on, min(cal, trk) IS the eta/E switch
    recommendation 4 endorses: crystals backward, tracker in the barrel."""
    for eta, expect_cal in ((-3.0, True), (-1.6, False), (0.0, False)):
        cal = reco.emcal_resolution(10.0, eta=eta)
        trk = reco.tracking_resolution(10.0, eta)
        assert bool(cal < trk) is expect_cal


def test_sweet_spot_electrons_are_insensitive_to_the_table():
    """The four sweet spots sit at eta = -2.9 to -1.6 with E' ~ 10 GeV;
    three of them are in the backward endcap, so the headline numbers do
    not move (code review F4)."""
    for eta in (-2.93, -2.92, -2.42):
        assert reco.emcal_resolution(9.9, eta=eta) == pytest.approx(
            reco.emcal_resolution(9.9), rel=1e-12)
