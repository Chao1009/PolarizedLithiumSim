"""Tests pinning the fixes of the 2026-08-28 code review of the
reconstruction chain (Report 2), each against an INDEPENDENT derivation or
an explicit measurability requirement rather than against the code itself.
"""

import numpy as np
import pytest

from polligen import hfs, reco, recopseudo as rp
from polligen import bookkeeping as bk
from polligen.sample import InclusiveSampler
from polligen.xsec import InclusiveKernel

from polli_fastsim import beams, fom
from polli_fastsim.polarized import toy_b1, toy_delta_gluon

CONFIG = beams.default_configs("6Li")[1]


# --- statistics ------------------------------------------------------------

def test_single_fill_bins_get_a_finite_positive_variance():
    """A bin populated by one fill only used to have var = 0 (R = w_f
    exactly, so the observed-count delta method collapsed) and hence an
    infinite weight in the fit.  The expected-count form is finite and
    equals sigma_P^2/N to leading order."""
    counts = np.array([[100.0, 50.0, 0.0, 3.0], [90.0, 60.0, 5.0, 0.0]])
    r, var, sig2, pbar = reco.spin_state_ratio(counts, [0.5, 0.5], [0.6, -1.2])
    assert np.all(np.isfinite(var)) and np.all(var > 0)
    # leading order: var ~ sigma_P^2 / N for a well-populated bin
    n0 = counts[:, 0].sum()
    assert abs(var[0] / (sig2 / n0) - 1.0) < 0.05
    # and the fit no longer raises on such a histogram
    edges = np.linspace(0.0, 2.0 * np.pi, 5)
    fit = reco.harmonic_ratio_fit(counts, [0.5, 0.5], [0.6, -1.2], edges)
    assert np.isfinite(fit["amp"]) and np.isfinite(fit["err"])


def test_variance_reproduces_the_poisson_spread_of_the_ratio():
    """The delta-method variance with expected counts must reproduce the
    spread of R over Poisson draws at the few-percent level."""
    rng = np.random.default_rng(3)
    mu = np.array([[4000.0, 3000.0], [3800.0, 3300.0]])
    rs = np.array([reco.spin_state_ratio(rng.poisson(mu), [0.5, 0.5],
                                         [0.6, -1.2])[0] for _ in range(4000)])
    _, var, _, _ = reco.spin_state_ratio(mu, [0.5, 0.5], [0.6, -1.2])
    ratio = rs.std(axis=0) / np.sqrt(var)
    assert np.all(np.abs(ratio - 1.0) < 0.06)


def test_ratio_inversion_jacobian_is_exact():
    """T = R (1 + u)/(sigma_P^2 - Pbar R): the propagated variance must
    equal (dT/dR)^2 var(R) with the exact derivative."""
    sig2, pbar, u = 0.81, -0.3, 0.02
    r = np.array([0.05, -0.02, 0.1])
    var = np.array([1e-4, 2e-4, 3e-4])
    t, var_t = reco._ratio_to_modulation(r, var, sig2, pbar, u=u, n_iter=60)
    t_exact = r * (1.0 + u) / (sig2 - pbar * r)
    dtdr = (1.0 + u) * sig2 / (sig2 - pbar * r) ** 2
    assert np.allclose(t, t_exact, rtol=1e-12)
    assert np.allclose(var_t, dtdr ** 2 * var, rtol=1e-9)


# --- measurability ---------------------------------------------------------

@pytest.fixture(scope="module")
def toy_library():
    rng = np.random.default_rng(11)
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1,
                           delta_func=toy_delta_gluon)
    analysis = fom.Scenario(lumi_fb_per_nucleon=10.0, pol_ion_tensor=0.6)
    sampler = InclusiveSampler(kern, CONFIG, rp.generator_scenario(analysis),
                               nx=24, nq2=18, q2_range=(0.7, 2e3))
    smp = hfs.toy_library_sample(sampler, 40000, rng)
    # a deliberately lossy response (calorimeters to |eta| = 2.5) so that
    # the capture bias the calibration must remove is visible on the toy
    resp = hfs.HadronResponse(noise_sigma=0.05, eta_cal=2.5, eta_track=2.5)
    lib = hfs.HFSLibrary(smp, resp, nx=16, nq2=12, rng=rng)
    return sampler, lib


def test_calibration_never_reads_the_true_cell(toy_library, monkeypatch):
    """HFSResponse(calibrate=True) must run without the truth-keyed map:
    the calibration is looked up at the reconstructed (x_mixed, Q2_e)."""
    sampler, lib = toy_library

    def forbidden(*a, **k):
        raise AssertionError("calibration keyed on the true (x, Q2) cell")
    monkeypatch.setattr(lib, "cell_means", forbidden)
    hr = hfs.HFSResponse(lib, calibrate=True)
    rng = np.random.default_rng(1)
    x = np.array([0.02, 0.05]); q2 = np.array([1.5, 3.0])
    y = q2 / (sampler.s * x)
    e_e = CONFIG.electron_energy
    kp = reco.electron_fourvector(x, y, sampler.s, e_e, np.zeros(2))
    th = np.arccos(kp[:, 3] / np.sqrt((kp[:, 1:] ** 2).sum(axis=1)))
    out = hr.hadronic(x, q2, y, kp[:, 0], th, e_e, sampler.s, rng)
    assert np.all(np.isfinite(out["y_sigma"]))


def test_reco_keyed_calibration_removes_the_scale_bias(toy_library):
    """On the library's own events the calibrated y_Sigma must be unbiased
    on average (median y_reco/y_true within a few % of 1) while the
    uncalibrated one carries the capture bias."""
    sampler, lib = toy_library
    rng = np.random.default_rng(5)
    n = 20000
    cell = rng.integers(0, sampler.xsec_flat.size, size=n)
    x = np.exp(sampler.logx_lo[cell] + rng.uniform(size=n)
               * (sampler.logx_hi[cell] - sampler.logx_lo[cell]))
    q2 = np.exp(sampler.logq2_lo[cell] + rng.uniform(size=n)
                * (sampler.logq2_hi[cell] - sampler.logq2_lo[cell]))
    ok = sampler._in_acceptance(x, q2) & (q2 / (sampler.s * x) > 0.05)
    x, q2 = x[ok], q2[ok]
    y = q2 / (sampler.s * x)
    e_e = CONFIG.electron_energy
    kp = reco.electron_fourvector(x, y, sampler.s, e_e, np.zeros(x.size))
    th = np.arccos(kp[:, 3] / np.sqrt((kp[:, 1:] ** 2).sum(axis=1)))
    med = {}
    for cal in (False, True):
        hr = hfs.HFSResponse(lib, calibrate=cal)
        yr, _ = hr.y_hadronic(x, q2, y, kp[:, 0], th, e_e, sampler.s,
                              np.random.default_rng(7))
        med[cal] = np.median(yr / y)
    assert med[False] < 0.995         # the toy response loses some Sigma
    assert abs(med[True] - 1.0) < 0.03


def test_target_mass_term_enters_the_transferred_sigma(toy_library):
    """The pseudo-event's true Sigma must carry the same m_N^2/(E + p) term
    the library's Sigma_true does, so that a captured FRACTION of 1 returns
    Sigma_true itself."""
    sampler, lib = toy_library
    hr = hfs.HFSResponse(lib, calibrate=False)
    lib.transfer = lambda x, q2, rng: (np.ones(np.size(x)), np.ones(np.size(x)),
                                       np.zeros(np.size(x)))
    hr.noise_sigma = 0.0
    x = np.array([0.03]); q2 = np.array([1.2]); y = q2 / (sampler.s * x)
    e_e = CONFIG.electron_energy
    kp = reco.electron_fourvector(x, y, sampler.s, e_e, np.zeros(1))
    th = np.arccos(kp[:, 3] / np.sqrt((kp[:, 1:] ** 2).sum(axis=1)))
    out = hr.hadronic(x, q2, y, kp[:, 0], th, e_e, sampler.s,
                      np.random.default_rng(0))
    p_u = CONFIG.ion_momentum_per_nucleon
    m_n = 0.938272
    expect = 2.0 * e_e * y + m_n ** 2 / (np.hypot(p_u, m_n) + p_u)
    assert np.allclose(out["sigma"], expect, rtol=1e-12)


def test_selected_bins_of_the_r_script_use_measurable_criteria_only():
    """The 5R/7R scripts must not select bins on the truth amplitude."""
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1]
           / "scripts" / "money_cos2phi_reco.py").read_text()
    assert 'abs(summ["a_reco_bin"]) < 1e-5' not in src
    assert 'abs(summ["a_reco_bin"]) >= 1e-5' not in src


# --- detector frame ---------------------------------------------------------

def test_hadron_acceptance_is_applied_in_the_lab_frame():
    """A particle at the forward calorimeter edge in the HEAD-ON frame is
    accepted or not according to its LAB pseudorapidity, and the sums are
    formed back in the head-on frame."""
    resp = hfs.HadronResponse(noise_sigma=0.0, eff_track=1.0, eff_nhad=1.0)
    # a photon at head-on eta = 3.75 (outside 3.7) along -x: the crossing
    # angle rotates the ion side toward -x, so in the lab it sits INSIDE
    eta = 3.75
    th = 2.0 * np.arctan(np.exp(-eta))
    e = 5.0
    p4 = np.array([[e, -e * np.sin(th), 0.0, e * np.cos(th)],
                   [e, +e * np.sin(th), 0.0, e * np.cos(th)]])
    lab = resp.lab_eta(p4)
    assert lab[0] < 3.7 < lab[1]
    out, w = resp.reconstruct_particles(p4, np.array([22, 22]),
                                        np.zeros(2), np.random.default_rng(0))
    assert w[0] == 1.0 and w[1] == 0.0
    # the accepted photon comes back in the head-on frame (smeared energy,
    # same direction to the smearing precision)
    assert abs(out[0, 3] / np.sqrt((out[0, 1:] ** 2).sum()) - np.cos(th)) < 1e-6
    # xing = 0 reproduces the head-on acceptance
    resp0 = hfs.HadronResponse(noise_sigma=0.0, xing=0.0)
    _, w0 = resp0.reconstruct_particles(p4, np.array([22, 22]), np.zeros(2),
                                        np.random.default_rng(0))
    assert w0.sum() == 0.0


def test_beam_energy_spread_has_the_analysis_sign(monkeypatch):
    """An event made at E_e (1 + d) and reconstructed with the nominal E_e
    has Q2_e LOW by (1 + d) and 1 - y_e HIGH by it, so with every detector
    smearing switched off Q2_e (1 - y_e) is exactly the true product.
    (Both carried the same sign until 2026-08-28.)"""
    monkeypatch.setattr(reco, "tracking_angular_resolution",
                        lambda eta: np.zeros_like(np.asarray(eta, float)))
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1, delta_func=toy_delta_gluon)
    analysis = fom.Scenario(lumi_fb_per_nucleon=10.0, pol_ion_tensor=0.6)
    sampler = InclusiveSampler(kern, CONFIG, rp.generator_scenario(analysis),
                               nx=12, nq2=9, q2_range=(0.7, 2e3))
    m = rp.RecoModel(beam_e_spread=0.02, y_method="electron",
                     emcal_stoch=0.0, emcal_const=0.0, eid=False)
    r = rp.RecoResponse(sampler, m, n_mc_per_cell=30, rng=np.random.default_rng(1))
    ok = np.isfinite(r.x_reco) & (r.y > 0.05) & (r.y < 0.9)
    q2_ratio = r.q2_reco[ok] / r.q2[ok]
    omy_ratio = (1.0 - r.y_reco[ok]) / (1.0 - r.y[ok])
    assert np.allclose(q2_ratio * omy_ratio, 1.0, rtol=1e-6)
    assert q2_ratio.std() > 0.01          # the spread is on


def test_beam_fourvectors_use_the_physical_nuclear_mass():
    _, p = reco.beam_fourvectors(CONFIG)
    m2 = reco.mdot(p, p)
    assert abs(np.sqrt(m2) - beams.NUCLEUS_MASS[("6Li", 6, 3)]) < 1e-9


def test_legacy_high_divergence_constants_agree():
    from polli_fastsim.farforward import HIGH_DIVERGENCE
    assert reco.SIGMA_THETA_HD == HIGH_DIVERGENCE.sigma_theta


# --- self-calibration of the spin axis ------------------------------------

def test_inclusive_axis_tilt_shows_up_as_the_sin_term():
    """An azimuthal misalignment delta of the assumed spin axis turns
    A cos 2phi' into A cos 2d cos 2phi' + A sin 2d sin 2phi': the fitted
    sin/cos ratio must return tan 2d at the reconstructed level."""
    delta = 0.10
    amp = 8e-3
    edges = np.linspace(0.0, 2.0 * np.pi, 25)
    n_tot = 4e9
    # generate with the axis at phi_S + delta: the modulation in the
    # ANALYSIS azimuth phi' = phi_e - phi_S is cos 2(phi' - delta)
    acc = lambda ph: 1.0 + 0.03 * np.cos(2 * ph)  # noqa: E731
    sub = 32
    lo, hi = edges[:-1], edges[1:]
    s = lo[:, None] + (np.arange(sub)[None, :] + 0.5) / sub * (hi - lo)[:, None]
    mu = []
    for pf, ff in ((0.6, 0.5), (-1.2, 0.5)):
        dens = acc(s) * (1.0 + pf * amp * np.cos(2.0 * (s - delta)))
        mu.append(n_tot * ff * (hi - lo) / (2 * np.pi) * dens.mean(axis=1))
    mu = np.array(mu)
    fit = reco.harmonic_ratio_fit(mu, [0.5, 0.5], [0.6, -1.2], edges, with_sin=True)
    assert abs(fit["amp_sin"] / fit["amp"] - np.tan(2 * delta)) < 1e-3
    assert abs(fit["phase"] - delta) < 1e-4


# --- hadronic kinematics against an independent derivation ------------------

def test_hadronic_methods_reproduce_the_exact_kinematics_by_hand():
    """Massless e + N at (x, y): compute Sigma and pT_h from the electron
    alone by four-momentum conservation and check every method returns
    the input (Sigma, JB, DA), using the textbook formulas written here
    independently of hfs.hadronic_kinematics."""
    e_e, p_u = 10.0, 99.5
    s = 4 * e_e * p_u
    x, y = 0.03, 0.2
    q2 = s * x * y
    kp = reco.electron_fourvector(np.array([x]), np.array([y]), s, e_e, np.array([0.3]))
    e_p = kp[0, 0]
    th = np.arccos(kp[0, 3] / np.sqrt((kp[0, 1:] ** 2).sum()))
    # hadronic system X = k + P - k' (massless target)
    k = np.array([e_e, 0, 0, -e_e]); P = np.array([p_u, 0, 0, p_u])
    X = k + P - kp[0]
    sigma = X[0] - X[3]
    ptx, pty = X[1], X[2]
    kin = hfs.hadronic_kinematics(np.array([sigma]), np.array([ptx]), np.array([pty]),
                                  np.array([e_p]), np.array([th]), e_e, s)
    # by hand
    y_jb = sigma / (2 * e_e)
    pt_h = np.hypot(ptx, pty)
    q2_jb = pt_h ** 2 / (1 - y_jb)
    y_sig = sigma / (sigma + e_p * (1 - np.cos(th)))
    tan_g2 = sigma / pt_h                      # tan(gamma_h/2)
    tan_e2 = np.tan(th / 2.0)                  # theta from the ion direction
    y_da = tan_g2 / (tan_g2 + tan_e2)
    q2_da = 4 * e_e ** 2 / tan_e2 / (tan_g2 + tan_e2)
    for key, val in (("y_jb", y_jb), ("q2_jb", q2_jb), ("y_sigma", y_sig),
                     ("y_da", y_da), ("q2_da", q2_da)):
        assert abs(kin[key][0] / val - 1.0) < 1e-9, key
    for key in ("y_jb", "y_sigma", "y_da"):
        assert abs(kin[key][0] / y - 1.0) < 1e-9
    for key in ("q2_jb", "q2_da", "q2_e"):
        assert abs(kin[key][0] / q2 - 1.0) < 1e-9


# --- the tagging optics helper against the pricing script -------------------

def test_tagging_optics_point_reproduces_the_priced_optimum():
    """reco.tagging_optics_point must give the optimum that
    scripts/tagging_optics.py prints (r_h = 46.5 / 164.1 / 89.3,
    eps = 0.374 / 0.251 / 0.332, L/L_HA = 1/6.8 / 1/12.8 / 1/9.5)."""
    # Since 2026-08-29 the dispersive envelope term is priced on each
    # configuration's own (R12, D) rather than on the 18 x 275 pair alone.
    # It enters as the RATIO D / R12 -- 1.62e-2 / 1.35e-2 / 9.74e-3 m/m --
    # so the top row is unmoved (the old default WAS the measured 18 x 275
    # ratio, within 0.6% of the September-2024 0.30 / 30.6 it inherited)
    # and the two lower ones lose acceptance: 0.423 -> 0.374 and
    # 0.323 -> 0.251.  The retired single-lever optimum, which is what
    # every tagging number published before that date carries, is
    # tagging_optics_point(cfg, levers="18x275").
    expect = ((46.5, 0.374, 6.8, 0.363, 3.80), (164.1, 0.251, 12.8, 0.192, 1.80),
              (89.3, 0.332, 9.5, 0.117, 0.92))
    for cfg, (r_h, eps, one_over_l, env_x, env_y) in zip(beams.default_configs("6Li"), expect):
        t = reco.tagging_optics_point(cfg)
        assert abs(t["r_h"] / r_h - 1.0) < 0.02
        assert abs(t["acceptance"] - eps) < 0.005
        assert abs(1.0 / t["lumi_fraction"] - one_over_l) < 0.15
        assert abs(1e3 * t["env_x"] - env_x) < 0.005
        assert abs(1e3 * t["env_y"] - env_y) < 0.01
    # the single-lever behaviour is still reachable and still exact
    old = reco.tagging_optics_point(beams.default_configs("6Li")[1],
                                    levers="18x275")
    assert abs(old["r_h"] / 175.6 - 1.0) < 0.02
    assert abs(old["acceptance"] - 0.323) < 0.005


def test_coherent_response_accepts_an_anisotropic_divergence():
    sc = __import__("polligen.coherent", fromlist=["CoherentScenario"]).CoherentScenario()
    cfg = beams.default_configs("6Li")[0]
    top = reco.tagging_optics_point(cfg)
    cr = rp.CoherentResponse(sc, cfg, (top["sigma_x_eff"], top["sigma_y"]),
                             n_mc=20000, rng=np.random.default_rng(2),
                             cut_theta_xy=(top["env_x"], top["env_y"]))
    assert abs(cr.aspect - top["sigma_y"] / top["sigma_x_eff"]) < 1e-12
    assert abs(cr.acceptance - top["acceptance"]) < 0.03
    # the per-fill perturbation acts on the binding cut -- the horizontal
    # envelope here, through which the recoils escape
    v = cr.with_cut(eff_scale_xy=(1.02, 1.0))
    assert v.cut_theta_eff[0] > cr.cut_theta_eff[0]
    assert v.acceptance < cr.acceptance
    # ... and is a view of the same recoils (no re-draw)
    assert v._m is cr._m


def test_sigma_weighted_merge_of_p_and_n():
    """A p + n merge weights events by sigma_gen/n_events per sample."""
    rng = np.random.default_rng(0)
    e_e, p_u = 10.0, 99.5
    toy = hfs.ToyHFS(e_e, p_u)
    x = np.full(50, 0.02); q2 = np.full(50, 2.0); ph = rng.uniform(0, 2 * np.pi, 50)
    a = toy.generate(x, q2, ph, rng); a.meta["sigma_gen_mb"] = 2.0
    b = toy.generate(x, q2, ph, rng); b.meta["sigma_gen_mb"] = 1.0
    m = hfs.HFSSample.concatenate([a, b])
    assert abs(m.weight[:50].mean() / m.weight[50:].mean() - 2.0) < 1e-12
    assert abs(m.weight.mean() - 1.0) < 1e-12
