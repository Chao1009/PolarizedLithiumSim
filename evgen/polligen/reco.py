"""Measured quantities and reconstruction for the cos 2phi analyses.

This module is the seed of plans/07 WP3 ("reconstructed-level closure"):
it turns the generator-level truth of polligen (x, Q2, y, the physics
azimuth phi', the spin state) into the quantities an experiment actually
records, and implements the reconstruction an experiment actually runs
-- with SIMPLE detector effects (Gaussian resolutions, geometric
acceptance) but REALISTIC measured quantities and estimators.  The
companion report `reports/reconstruction_chain_report` describes the
full chain and audits what exists elsewhere in the code.

Contents
--------
Frames and four-vectors
  beam_fourvectors, electron_fourvector, spin_fourvector,
  boost_x, rotate_y, head_on_to_lab, lab_to_head_on, azimuth_about_z
Covariant azimuths (Bacchetta et al. JHEP 02 (2007) 093 conventions)
  azimuth_wrt_lepton_plane -- phi_S, phi_t, phi_h for any four-vector
  lab_azimuth_shortcut_error -- the O(gamma^2) / recoil-mass shifts
Scattered-electron reconstruction
  electron_method, electron_method_resolution, smear_electron,
  hadronic_y (parametrized Sigma/Jacquet-Blondel y), mixed_method
Spin-state-sorted estimator
  spin_state_ratio, harmonic_ratio_fit, err_harmonic_ratio
Coherent recoil and Roman-Pot emulation
  recoil_fourvector, t_from_fourvectors, tag_pt_cut, rp_measure,
  rp_hole_acceptance

Conventions.  Four-vectors are arrays (..., 4) ordered (E, px, py, pz)
with metric (+,-,-,-).  The HEAD-ON frame has the ion along +z and the
electron along -z (the frame of polli_fastsim.kinematics); the LAB
(detector) frame has the electron beam exactly along -z and the ion beam
at -xing in x (ePIC: xing = 25 mrad at IP6, horizontal; the far-forward
gun scan of tools/fullsim uses the same sign).  Polar angles theta are
measured from +z (ion direction); theta' = pi - theta is the angle from
the electron beam.  Per-nucleon (x, y, s) throughout, as in the rest of
the package; whole-nucleus four-momenta where a mass matters.
"""

import numpy as np

from polli_fastsim.kinematics import scattered_electron
from polli_fastsim.spectator import M_U

XING_IP6 = 25.0e-3   # rad, horizontal crossing angle at IP6
XING_IP8 = 35.0e-3   # rad, IP8

# 10-sigma near-beam cut expressed as an ANGULAR cut: the documented
# Roman-Pot pT cuts for 275 GeV protons (0.20 GeV high-acceptance,
# ~0.41-0.45 GeV high-divergence; plans/06 SS6.5) divided by 10 x 275 GeV.
# The Li beam divergence is undocumented (plans/04 #11): these are the
# proton-derived placeholders that `tag_pt_cut` scales to any beam.
SIGMA_THETA_HA = 0.20 / (10.0 * 275.0)   # 73 microrad
SIGMA_THETA_HD = 0.41 / (10.0 * 275.0)   # 149 microrad


# --- Minkowski helpers -----------------------------------------------------

def mdot(a, b):
    """Minkowski product a.b = E_a E_b - vec a . vec b (broadcasting)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return a[..., 0] * b[..., 0] - (a[..., 1:] * b[..., 1:]).sum(axis=-1)


def eps4(a, b, c, d):
    """epsilon_{mu nu rho sigma} a^mu b^nu c^rho d^sigma with
    epsilon_{0123} = +1, i.e. det of the 4x4 matrix with rows (a,b,c,d)."""
    a, b, c, d = np.broadcast_arrays(*(np.asarray(v, dtype=float)
                                       for v in (a, b, c, d)))
    return np.linalg.det(np.stack([a, b, c, d], axis=-2))


def fourvector(e, px, py, pz):
    return np.stack(np.broadcast_arrays(*(np.asarray(v, dtype=float)
                                          for v in (e, px, py, pz))),
                    axis=-1)


# --- frames ----------------------------------------------------------------

def beam_fourvectors(config, ion_mass=None):
    """(k, P) in the head-on frame: massless electron along -z with
    energy E_e, whole-nucleus ion along +z with momentum A * p_u and
    mass ion_mass (default A * M_U; pass 0.0 for the massless-target
    limit of the master formula)."""
    ion = config.ion
    p_a = ion.A * config.ion_momentum_per_nucleon
    m_a = ion.A * M_U if ion_mass is None else float(ion_mass)
    k = fourvector(config.electron_energy, 0.0, 0.0, -config.electron_energy)
    p = fourvector(np.hypot(p_a, m_a), 0.0, 0.0, p_a)
    return k, p


def electron_fourvector(x, y, s, e_energy, phi_e):
    """Scattered-electron four-vector in the head-on frame from the
    per-nucleon (x, y) and the lab azimuth phi_e about the ion axis
    (kinematics.scattered_electron gives E' and theta; theta ~ pi)."""
    e_p, theta, _eta = scattered_electron(np.asarray(x, dtype=float),
                                          np.asarray(y, dtype=float), s,
                                          e_energy)
    st, ct = np.sin(theta), np.cos(theta)
    phi_e = np.asarray(phi_e, dtype=float)
    return fourvector(e_p, e_p * st * np.cos(phi_e), e_p * st * np.sin(phi_e),
                      e_p * ct)


def spin_fourvector(phi_s, theta_s=np.pi / 2.0):
    """Spin (alignment-axis) four-vector S = (0, n) for the axis
    n(theta_s, phi_s); transverse to the ion beam when theta_s = pi/2, so
    S.P = 0 exactly for a beam along z.  The tensor axis is headless:
    only 2 phi_s matters downstream."""
    st = np.sin(theta_s)
    return fourvector(0.0, st * np.cos(phi_s), st * np.sin(phi_s),
                      np.cos(theta_s))


def boost_x(v, beta):
    """Boost along +x with velocity beta (E' = g(E - b px), px' = g(px - b E))."""
    v = np.asarray(v, dtype=float)
    g = 1.0 / np.sqrt(1.0 - beta * beta)
    out = v.copy()
    out[..., 0] = g * (v[..., 0] - beta * v[..., 1])
    out[..., 1] = g * (v[..., 1] - beta * v[..., 0])
    return out


def rotate_y(v, angle):
    """Rotate spatial components about y: x' = x cos a + z sin a,
    z' = -x sin a + z cos a."""
    v = np.asarray(v, dtype=float)
    c, s = np.cos(angle), np.sin(angle)
    out = v.copy()
    out[..., 1] = c * v[..., 1] + s * v[..., 3]
    out[..., 3] = -s * v[..., 1] + c * v[..., 3]
    return out


def head_on_to_lab(v, xing=XING_IP6):
    """Head-on -> detector frame: boost along x by sin(xing/2), then rotate
    about y by -xing/2.  The electron beam ends exactly along -z and the
    ion beam along (-sin xing, 0, cos xing) -- "-25 mrad in x", the ePIC
    convention.  Energies change by 1/cos(xing/2) - 1 = 8e-5."""
    return rotate_y(boost_x(v, np.sin(0.5 * xing)), -0.5 * xing)


def lab_to_head_on(v, xing=XING_IP6):
    """Inverse of head_on_to_lab: the "head-on frame" transformation every
    ePIC inclusive analysis applies before computing kinematics."""
    return boost_x(rotate_y(v, 0.5 * xing), -np.sin(0.5 * xing))


def azimuth_about_z(v):
    """Naive lab azimuth atan2(py, px) about the frame's z axis."""
    v = np.asarray(v, dtype=float)
    return np.arctan2(v[..., 2], v[..., 1])


def polar_from_electron_beam(v):
    """theta' = angle from the electron-beam direction (-z) in the frame
    of v (pi - theta)."""
    v = np.asarray(v, dtype=float)
    p = np.sqrt((v[..., 1:] ** 2).sum(axis=-1))
    return np.arccos(np.clip(-v[..., 3] / np.maximum(p, 1e-300), -1.0, 1.0))


# --- covariant azimuths ----------------------------------------------------

def transverse_part(v, P, q):
    """g_perp^{mu nu} v_nu with the Bacchetta et al. projector
    g_perp = g - (q P + P q)/(P.q (1+g2)) + g2/(1+g2) (q q/Q2 - P P/M2),
    g2 = M2 Q2/(P.q)^2 (= 4 M^2 x^2/Q^2).  The P P/M2 term is dropped
    when M2 = 0 (massless target), where g2 = 0 anyway."""
    v = np.asarray(v, dtype=float)
    pq = mdot(P, q)
    q2 = -mdot(q, q)
    m2 = mdot(P, P)
    g2 = m2 * q2 / (pq * pq)
    pv = mdot(P, v)
    qv = mdot(q, v)
    out = (v - (q * pv[..., None] + P * qv[..., None])
           / (pq * (1.0 + g2))[..., None]
           + (g2 / (1.0 + g2))[..., None] * q * (qv / q2)[..., None])
    massive = np.abs(m2) > 1e-12
    if np.any(massive):
        corr = ((g2 / (1.0 + g2)) * pv / np.where(massive, m2, 1.0))[..., None] * P
        out = out - np.where(massive[..., None], corr, 0.0)
    return out


def azimuth_wrt_lepton_plane(k, kp, P, v):
    """Covariant azimuthal angle of the four-vector v about the virtual
    photon, measured from the lepton plane (Bacchetta et al. 2007 Eqs.
    2.4-2.5 applied to v; v = S gives phi_S, v = P' of a recoil gives
    phi_t, v = P_h gives phi_h):

        cos phi = -l_perp . v_perp / (|l_perp| |v_perp|),
        sin phi = -eps_{mu nu rho sigma} l^mu v^nu P^rho q^sigma
                  / (P.q sqrt(1+g2)) / (|l_perp| |v_perp|),  l = k.

    Sign convention fixed so that in the head-on frame with a transverse
    spin axis at lab azimuth phi_s and the electron at lab azimuth phi_e,
    phi_S = phi_e - phi_s exactly for a massless target (tested).  Only
    2 phi enters the tensor observables, so the orientation convention
    cannot change any result.  Returns phi in (-pi, pi]."""
    k, kp, P, v = (np.asarray(a, dtype=float) for a in (k, kp, P, v))
    q = k - kp
    lp = transverse_part(k, P, q)
    vp = transverse_part(v, P, q)
    norm = np.sqrt(np.maximum(-mdot(lp, lp), 1e-300)
                   * np.maximum(-mdot(vp, vp), 1e-300))
    pq = mdot(P, q)
    g2 = mdot(P, P) * (-mdot(q, q)) / (pq * pq)
    cos_phi = -mdot(lp, vp) / norm
    sin_phi = eps4(k, v, P, q) / (pq * np.sqrt(1.0 + g2)) / norm
    return np.arctan2(sin_phi, cos_phi)


def lab_azimuth_shortcut_error(k, kp, P, v, v_lab_phi):
    """phi(covariant) - (phi_e - phi_v_lab): the error of the lab-angle
    shortcut for the four-vector v whose head-on-frame azimuth is
    v_lab_phi.  Zero for a massless target and transverse S; O(gamma^2)
    for a massive target; O(M_A^2/(P.q)) for a heavy recoil."""
    phi_e = azimuth_about_z(kp)
    phi_cov = azimuth_wrt_lepton_plane(k, kp, P, v)
    return np.angle(np.exp(1j * (phi_cov - (phi_e - v_lab_phi))))


# --- scattered-electron reconstruction ------------------------------------

def electron_method(e_prime, theta, e_energy, s):
    """(Q2, y, x) from the scattered electron alone (theta from the ion
    direction): Q2 = 2 E_e E'(1 + cos theta), y = 1 - E'(1 - cos theta)/(2E_e)."""
    e_prime = np.asarray(e_prime, dtype=float)
    ct = np.cos(np.asarray(theta, dtype=float))
    q2 = 2.0 * e_energy * e_prime * (1.0 + ct)
    y = 1.0 - e_prime * (1.0 - ct) / (2.0 * e_energy)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.where(y > 0, q2 / (s * np.where(y > 0, y, 1.0)), np.nan)
    return q2, y, x


def electron_method_resolution(y, theta_prime, de_over_e, dtheta,
                               de_beam_over_e=0.0):
    """Relative resolutions of the electron method (linear propagation):

        dQ2/Q2 = dE'/E'  (+)  cot(theta'/2) dtheta'
        dy/y   = ((1-y)/y) [ dE'/E' (+) tan(theta'/2) dtheta' (+) dE_e/E_e ]
        dx/x   = dQ2/Q2 (+) dy/y

    theta' = angle from the electron beam; (+) = quadrature.  The (1-y)/y
    factor is why x is unmeasurable from the electron alone at y ~ 0.01:
    the electron carries 99% of the beam energy and y is a 1% difference."""
    y = np.asarray(y, dtype=float)
    th = np.asarray(theta_prime, dtype=float)
    dq2 = np.hypot(de_over_e, dtheta / np.tan(0.5 * th))
    dy = (1.0 - y) / y * np.sqrt(de_over_e ** 2
                                 + (np.tan(0.5 * th) * dtheta) ** 2
                                 + de_beam_over_e ** 2)
    dx = np.hypot(dq2, dy)
    return dq2, dy, dx


def emcal_resolution(e_prime, stoch=0.02, const=0.01, noise=0.0):
    """dE/E of an EM calorimeter: stoch/sqrt(E) (+) const (+) noise/E.
    Defaults are PbWO4-class (backward ePIC EMCal specification scale,
    eta < -2).  The barrel and forward calorimeters are coarser (Yellow
    Report requirements 7-12%/sqrt(E) (+) 1-3%): RecoModel applies these
    defaults at every eta, which is optimistic outside the backward endcap
    (code review 2026-08-25, docs/code_review_2026-08-25.md F4)."""
    e = np.maximum(np.asarray(e_prime, dtype=float), 1e-9)
    return np.sqrt((stoch / np.sqrt(e)) ** 2 + const ** 2 + (noise / e) ** 2)


def tracking_resolution(e_prime, eta, a=None, b=None):
    """dp/p = sqrt((a p)^2 + b^2) with the eta-piecewise ePIC-like table of
    fastsim/scripts/money_delta_20260729.py (tracking only).  PLACEHOLDER:
    that table has no published source ("provided in parent communication"
    per its docstring); replace by the ePIC full-simulation values when
    available (code review 2026-08-25)."""
    eta = np.asarray(eta, dtype=float)
    p = np.asarray(e_prime, dtype=float)
    if a is None or b is None:
        a = np.full_like(eta, 0.0005)
        b = np.full_like(eta, 0.005)
        for lo, hi, av, bv in ((-np.inf, -3.5, 0.0020, 0.030),
                               (-3.5, -2.5, 0.0010, 0.030),
                               (-2.5, -1.0, 0.0005, 0.010),
                               (1.0, 2.5, 0.0005, 0.010),
                               (2.5, np.inf, 0.0010, 0.025)):
            m = (eta >= lo) & (eta < hi)
            a = np.where(m, av, a)
            b = np.where(m, bv, b)
    return np.sqrt((a * p) ** 2 + b ** 2)


def tracking_angular_resolution(eta):
    """Track direction resolution sigma_theta [rad] (also the transverse
    direction resolution used for phi), eta-piecewise table of
    fastsim/scripts/money_delta_20260729.py: 5/3/2/1/2/3/5 mrad.
    PLACEHOLDER without a published source; it sets the electron-method
    Q2 resolution at the low-y sweet spots (cot(theta'/2) dtheta' = 5%
    at theta' = 0.1 rad for 3 mrad, 2% for 1 mrad -- code review
    2026-08-25, F3)."""
    eta = np.asarray(eta, dtype=float)
    out = np.full_like(eta, 1.0e-3)
    for lo, hi, val in ((-np.inf, -3.5, 5.0e-3), (-3.5, -2.5, 3.0e-3),
                        (-2.5, -1.0, 2.0e-3), (1.0, 2.5, 2.0e-3),
                        (2.5, 3.5, 3.0e-3), (3.5, np.inf, 5.0e-3)):
        out = np.where((eta >= lo) & (eta < hi), val, out)
    return out


def eps_eid(eta):
    """Electron-ID efficiency eps_eID(eta), linearly interpolated between
    the ATHENA (JINST 17 (2022) P10019, Table 5) / ECCE (NIM A 1055 (2023)
    168464, Sec. 3.5.2) anchors used by fastsim/scripts/
    money_delta_20260729.py; zero outside |eta| > 3.5.  No official ePIC
    curve exists (pCDR v1)."""
    eta = np.asarray(eta, dtype=float)
    eta_pts = np.array([-3.5, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 3.5])
    eps_pts = np.array([0.85, 0.92, 0.95, 0.93, 0.90, 0.90, 0.85, 0.80, 0.70])
    return np.where((eta < -3.5) | (eta > 3.5), 0.0,
                    np.interp(eta, eta_pts, eps_pts))


def smear_electron(e_prime, theta, phi, de_over_e, dtheta, dphi, rng):
    """Gaussian smearing of the measured (E', theta, phi) of the scattered
    electron; resolutions are absolute (rad) for angles, relative for E'."""
    e_prime = np.asarray(e_prime, dtype=float)
    n = np.broadcast(e_prime, theta, phi).shape
    e_s = e_prime * (1.0 + de_over_e * rng.standard_normal(n))
    th_s = np.asarray(theta, dtype=float) + dtheta * rng.standard_normal(n)
    ph_s = np.asarray(phi, dtype=float) + dphi * rng.standard_normal(n)
    return e_s, np.clip(th_s, 1e-9, np.pi - 1e-9), np.mod(ph_s, 2.0 * np.pi)


def hadronic_y(y_true, rel_res, rng, floor=1e-4):
    """Parametrized hadronic-method y (Jacquet-Blondel / Sigma method):
    y_had = y_true (1 + rel_res N(0,1)).  Stand-in for the hadronic final
    state polligen does not generate.  Sources for rel_res (refs/README.md,
    2026-08-25): the ePIC kinematic-fit study smears the hadronic
    E - p_z sum with sigma(delta_h)/delta_h = 25% (S. Maple, ePIC seminar
    Dec 2024, slide 44; Djangoh 18x275, Q2 > 1 GeV2); the ATHENA proposal
    (JINST 17 (2022) P10019, Sec. 3.1 / Fig. 22) uses e-Sigma or DA for
    y <~ 0.1 with ~25% y resolution at y ~ 0.01 improving to ~10% at
    y ~ 0.1, and quotes 20-30% for JB throughout; the ATHENA fast
    simulation of Arratia et al. (NIM A 1025 (2022) 166164, Fig. 5,
    Q2 > 200 GeV2) gives RMS(y)/y ~ 0.13-0.17 for ISigma/DA/JB at
    y ~ 0.05-0.2.  rel_res ~ 0.15-0.30 therefore brackets the documented
    expectation; 0.20 is the default, 0.25 the ePIC study's own value."""
    y = np.asarray(y_true, dtype=float)
    return np.maximum(y * (1.0 + rel_res * rng.standard_normal(y.shape)),
                      floor)


def mixed_method(q2_electron, y_hadronic, s):
    """x from the electron-method Q2 and the hadronic-method y (the eSigma
    / 'mixed' method used at HERA at low y): x = Q2_e / (s y_had)."""
    return np.asarray(q2_electron, dtype=float) / (s * np.asarray(y_hadronic,
                                                                    dtype=float))


# --- spin-state-sorted harmonic estimator ---------------------------------

def spin_state_ratio(counts, lumis, pzz):
    """Acceptance-free tensor ratio per bin.

    counts: array (F, K) of counts per fill type f and bin i; lumis: (F,)
    integrated luminosities (only their ratios matter); pzz: (F,) fill
    tensor polarizations.  With Pbar = sum_f L_f P_f / sum_f L_f and
    w_f = P_f - Pbar, the bin-wise luminosity-weighted ratio

        R_i = sum_f w_f N_fi / sum_f N_fi

    is, for yields N_fi = L_f eps_i sigma_i [1 + P_f T_i],

        R_i = sigma_P^2 T_i / (1 + Pbar T_i)   EXACTLY,
        sigma_P^2 = sum_f L_f (P_f - Pbar)^2 / sum_f L_f,

    because sum_f L_f w_f = 0: any phi-dependent acceptance/efficiency
    eps_i cancels bin by bin, and a relative-luminosity error enters
    only through Pbar, i.e. as a bin-INDEPENDENT offset of R plus a
    second-order rescaling (delta x Pbar/sigma_P^2).  Returns
    (R, var_R, sigma_P2, Pbar) with var_R from linear error propagation
    of Poisson counts: var(R) = sum_f ((w_f - R)/sum N)^2 N_f.
    """
    n = np.asarray(counts, dtype=float)
    lum = np.asarray(lumis, dtype=float)
    p = np.asarray(pzz, dtype=float)
    lw = lum / lum.sum()
    pbar = float((lw * p).sum())
    w = (p - pbar)[:, None]
    sig2 = float((lw * (p - pbar) ** 2).sum())
    num = (w * n).sum(axis=0)
    den = n.sum(axis=0)
    live = den > 0                      # empty bins (e.g. inside a
    den_safe = np.where(live, den, 1.0)  # Roman-Pot cutout) get weight 0
    r = np.where(live, num / den_safe, 0.0)
    var = np.where(live, (((w - r[None, :]) / den_safe[None, :]) ** 2
                          * n).sum(axis=0), np.inf)
    return r, var, sig2, pbar


def _ratio_to_modulation(r, var, sig2, pbar, u=0.0, n_iter=4):
    """Invert R = sigma_P^2 T / (1 + u + Pbar T) for T bin by bin
    (fixed-point iteration; u = known spin-independent modulation of the
    denominator, 0 for inclusive DIS), propagating the variance."""
    t = r / sig2
    for _ in range(n_iter):
        t = r * (1.0 + u + pbar * t) / sig2
    scale = (1.0 + u + pbar * t) / sig2
    return t, var * scale * scale


def harmonic_ratio_fit(counts, lumis, pzz, edges, with_sin=False,
                       denominator_correction=True):
    """Fit T_i = kappa + A cos 2phi_i (+ A_s sin 2phi_i) to the spin-state
    ratio inverted for the modulation T (R = sigma_P^2 T/(1 + Pbar T),
    exact; `denominator_correction=False` uses T = R/sigma_P^2), by
    weighted LSQ.  Returns {"amp": A_hat, "err": dA, "const": kappa_hat,
    "sigma_p2", "pbar"} (+ "amp_sin", "phase" with with_sin), with the
    finite-bin dilution sin(w)/w (w = half-width in 2phi) divided out as
    in estimators.cos2phi_fit_binned.  Two fills (P+, P0) with equal
    luminosity give dA = 2 sqrt(2/N) / (P+ - P0) / dilution: the
    m=0-enriched fill (P0 = -2 P+ for the same purity) makes this 1.5x
    BETTER than the single-fill fit sqrt(2/N)/P+.  A residual sin 2phi
    term measures the azimuthal misalignment of the assumed spin axis."""
    r, var, sig2, pbar = spin_state_ratio(counts, lumis, pzz)
    if denominator_correction:
        t, var_t = _ratio_to_modulation(r, var, sig2, pbar)
    else:
        t, var_t = r / sig2, var / sig2 ** 2
    edges = np.asarray(edges, dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    cols = [np.ones_like(centers), np.cos(2.0 * centers)]
    if with_sin:
        cols.append(np.sin(2.0 * centers))
    design = np.vstack(cols).T
    wgt = 1.0 / np.sqrt(np.maximum(var_t, 1e-300))
    coef, *_ = np.linalg.lstsq(design * wgt[:, None], t * wgt, rcond=None)
    cov = np.linalg.inv((design * wgt[:, None]).T @ (design * wgt[:, None]))
    w = 0.5 * (edges[1] - edges[0])
    dil = np.sin(2.0 * w) / (2.0 * w)
    out = {"amp": coef[1] / dil, "err": np.sqrt(cov[1, 1]) / dil,
           "const": coef[0], "sigma_p2": sig2, "pbar": pbar}
    if with_sin:
        out["amp_sin"] = coef[2] / dil
        out["phase"] = 0.5 * np.arctan2(coef[2], coef[1])
    return out


def _bin_dilutions(alpha_edges, beta_edges):
    ae = np.asarray(alpha_edges, dtype=float)
    be = np.asarray(beta_edges, dtype=float)
    wa, wb = 0.5 * (ae[1] - ae[0]), 0.5 * (be[1] - be[0])
    return (np.sin(wa) / wa, np.sin(wb) / wb,
            np.sin(2 * wa) / (2 * wa), np.sin(2 * wb) / (2 * wb))


def basis_2d(alpha_edges, beta_edges, beta_means=None):
    """Bin-averaged harmonic basis on the (Ka x Kb) grid, flattened in C
    order: {"e": <cos 2a>, "t": <cos 2b>, "m": <cos(a+b)>,
    "u1": <cos(a-b)>, "u2": <cos 2(a-b)>}.  The alpha average is analytic
    (uniform in-bin distribution).  The beta averages are analytic too
    unless `beta_means` = {"c1","s1","c2","s2"} (Kb,) supplies the
    ACCEPTANCE-WEIGHTED in-bin means <cos b>, <sin b>, <cos 2b>, <sin 2b>
    of the true beta of the events reconstructed into each beta bin --
    required when the acceptance varies strongly across a bin (the
    Roman-Pot cutout: x25 across beta), and the way the smearing of beta
    enters the response (the fit then estimates the UNSMEARED
    coefficients, as an MC-corrected analysis does)."""
    ae = np.asarray(alpha_edges, dtype=float)
    be = np.asarray(beta_edges, dtype=float)
    ac = 0.5 * (ae[:-1] + ae[1:])
    bc = 0.5 * (be[:-1] + be[1:])
    d1a, d1b, d2a, d2b = _bin_dilutions(ae, be)
    if beta_means is None:
        c1, s1 = np.cos(bc) * d1b, np.sin(bc) * d1b
        c2, s2 = np.cos(2 * bc) * d2b, np.sin(2 * bc) * d2b
        c2t = c2
    else:
        c1, s1, c2, s2 = (np.asarray(beta_means[k], dtype=float)
                          for k in ("c1", "s1", "c2", "s2"))
        # template basis of a t-dependent a_t (recopseudo.basis_means)
        c2t = np.asarray(beta_means.get("c2t", c2), dtype=float)
    ca1, sa1 = (np.cos(ac) * d1a)[:, None], (np.sin(ac) * d1a)[:, None]
    ca2, sa2 = (np.cos(2 * ac) * d2a)[:, None], (np.sin(2 * ac) * d2a)[:, None]
    ones_a = np.ones((ac.size, 1))
    return {"e": (ca2 * np.ones((1, bc.size))).ravel(),
            "t": (ones_a * c2t[None, :]).ravel(),
            "m": (ca1 * c1[None, :] - sa1 * s1[None, :]).ravel(),
            "u1": (ca1 * c1[None, :] + sa1 * s1[None, :]).ravel(),
            "u2": (ca2 * c2[None, :] + sa2 * s2[None, :]).ravel()}


def unpolarized_modulation_2d(alpha_edges, beta_edges, u1, u2,
                              beta_means=None):
    """Bin-averaged spin-independent modulation u(alpha - beta) =
    u1 cos(alpha-beta) + u2 cos 2(alpha-beta), flattened in C order."""
    b = basis_2d(alpha_edges, beta_edges, beta_means)
    return u1 * b["u1"] + u2 * b["u2"]


def harmonic_ratio_fit_2d(counts, lumis, pzz, alpha_edges, beta_edges,
                          u_coeffs=None, beta_means=None):
    """Two-azimuth version for the coherent channel: counts (F, Ka, Kb)
    per fill in bins of alpha = phi_e - phi_S and beta = phi_t - phi_S.
    The spin-sorted yields are modelled as
        N_f ~ eps(alpha,beta) [1 + u(alpha-beta) + P_f T(alpha,beta)],
        T = kappa + a_e cos 2alpha + a_t cos 2beta + a_m cos(alpha+beta),
    with u = u1 cos(alpha-beta) + u2 cos 2(alpha-beta) the spin-
    independent lepton-plane/recoil-plane modulation (u_coeffs = (u1, u2),
    known from the unpolarized analysis; None = 0).  The bin-wise ratio
    R = sigma_P^2 T/(1 + u + Pbar T) is inverted for the bin-averaged T
    and fitted linearly on the bin-averaged basis (basis_2d; pass
    `beta_means` from the MC response for an acceptance-shaped beta).
    Returns {"const", "a_e", "a_t", "a_m", "err_e", "err_t", "err_m",
    "cov", "sigma_p2", "pbar"}."""
    n = np.asarray(counts, dtype=float)
    nf, ka, kb = n.shape
    r, var, sig2, pbar = spin_state_ratio(n.reshape(nf, ka * kb), lumis, pzz)
    basis = basis_2d(alpha_edges, beta_edges, beta_means)
    u = (u_coeffs[0] * basis["u1"] + u_coeffs[1] * basis["u2"]
         if u_coeffs is not None else 0.0)
    t, var_t = _ratio_to_modulation(r, var, sig2, pbar, u=u)
    design = np.vstack([np.ones(ka * kb), basis["e"], basis["t"],
                        basis["m"]]).T
    wgt = 1.0 / np.sqrt(np.maximum(var_t, 1e-300))
    coef, *_ = np.linalg.lstsq(design * wgt[:, None], t * wgt, rcond=None)
    cov = np.linalg.inv((design * wgt[:, None]).T @ (design * wgt[:, None]))
    err = np.sqrt(np.diag(cov))
    return {"const": coef[0], "a_e": coef[1], "a_t": coef[2], "a_m": coef[3],
            "err_e": err[1], "err_t": err[2], "err_m": err[3], "cov": cov,
            "sigma_p2": sig2, "pbar": pbar}


def err_harmonic_ratio(n_total, pzz_list, lumi_fractions=None, nbins=24):
    """Analytic statistical error of harmonic_ratio_fit for n_total events
    shared between fill types with tensor polarizations pzz_list and
    luminosity fractions lumi_fractions (equal by default):
    dA = sqrt(2/N) / sigma_P / dilution  (sigma_P = rms of P_f about the
    luminosity-weighted mean)."""
    p = np.asarray(pzz_list, dtype=float)
    f = (np.full(p.size, 1.0 / p.size) if lumi_fractions is None
         else np.asarray(lumi_fractions, dtype=float) / np.sum(lumi_fractions))
    pbar = (f * p).sum()
    sig = np.sqrt((f * (p - pbar) ** 2).sum())
    w = np.pi / nbins
    dil = np.sin(2.0 * w) / (2.0 * w)
    return np.sqrt(2.0 / np.asarray(n_total, dtype=float)) / sig / dil


def expected_counts_by_fill(n_total, pzz_list, amp, edges, acceptance=None,
                            lumi_fractions=None, const=0.0):
    """Exact expected phi'-bin counts per fill type for pseudo-experiments:
    N_fi = L_f eps_i (1 + P_f (const + amp cos 2phi'))  integrated over
    the bin; `acceptance(phi)` is an arbitrary smooth efficiency
    (default 1) -- the thing the ratio estimator must cancel."""
    p = np.asarray(pzz_list, dtype=float)
    f = (np.full(p.size, 1.0 / p.size) if lumi_fractions is None
         else np.asarray(lumi_fractions, dtype=float) / np.sum(lumi_fractions))
    edges = np.asarray(edges, dtype=float)
    nsub = 32
    sub = edges[:-1, None] + (np.arange(nsub)[None, :] + 0.5) / nsub \
        * (edges[1:] - edges[:-1])[:, None]
    eps = np.ones_like(sub) if acceptance is None else acceptance(sub)
    base = eps.mean(axis=1)
    mod = (eps * np.cos(2.0 * sub)).mean(axis=1)
    dphi = (edges[1:] - edges[:-1]) / (2.0 * np.pi)
    out = np.empty((p.size, edges.size - 1))
    for k, (pf, ff) in enumerate(zip(p, f)):
        out[k] = n_total * ff * dphi * ((1.0 + pf * const) * base + pf * amp * mod)
    return out, f


# --- coherent recoil and Roman-Pot emulation ------------------------------

def recoil_fourvector(t, phi_t, x_pom, P):
    """Intact-recoil four-vector P' for momentum transfer t (< 0), recoil
    azimuth phi_t about the ion axis, and plus-momentum fraction
    x_L = 1 - x_pom, from the exact light-cone relation
        t = -[pT^2 + M^2 (1 - x_L)^2] / x_L,
    so pT^2 = -t x_L - M^2 (1 - x_L)^2 (requires |t| >= t_min =
    M^2 x_pom^2/(1 - x_pom))."""
    P = np.asarray(P, dtype=float)
    t = np.asarray(t, dtype=float)
    xl = 1.0 - np.asarray(x_pom, dtype=float)
    m2 = mdot(P, P)
    pt2 = -t * xl - m2 * (1.0 - xl) ** 2
    if np.any(pt2 < 0):
        raise ValueError("|t| below t_min for the requested x_pom")
    pt = np.sqrt(pt2)
    pplus = xl * (P[..., 0] + P[..., 3])
    pminus = (pt2 + m2) / pplus
    return fourvector(0.5 * (pplus + pminus), pt * np.cos(phi_t),
                      pt * np.sin(phi_t), 0.5 * (pplus - pminus))


def t_from_fourvectors(P, Pprime):
    """t = (P - P')^2 exactly."""
    d = np.asarray(P, dtype=float) - np.asarray(Pprime, dtype=float)
    return mdot(d, d)


def tag_pt_cut(sigma_theta, p_per_nucleon, a_beam=6, n_sigma=10.0):
    """Near-beam cut on the NUCLEUS transverse momentum implied by an
    angular beam envelope: pT_cut = n_sigma * sigma_theta * A * p_u.
    The documented 0.20 GeV (HA) is a proton number at 275 GeV; the same
    envelope on a 6Li nucleus at p_u = 137.5 GeV/u is 0.60 GeV."""
    return n_sigma * sigma_theta * a_beam * np.asarray(p_per_nucleon,
                                                        dtype=float)


def rp_measure(Pprime, P, sigma_theta_xy, n_sigma=10.0, rng=None,
               shape="rectangle", sigma_pos_over_l=0.0,
               cut_scale_xy=(1.0, 1.0)):
    """Emulated Roman-Pot measurement of an intact recoil at R ~ 1.

    The pots see the transverse displacement of the track from the beam
    orbit; the momentum (x_L ~ 1) is unresolved.  What the reconstruction
    returns is the angle pair at the IP relative to the beam axis,
    smeared by the bunch's own angular deviation (beam divergence,
    sigma_theta_xy = (sx, sy) rad -- the dominant resolution: the ePIC
    far-forward simulation with all beam effects gives dpT ~ 40 MeV for
    275 GeV protons with the detector alone at the 1% level, Jentsch
    DIS 2023 [refs/]; ZEUS's LPS was likewise beam-spread dominated) and
    an optional position term sigma_pos/L_eff.  Acceptance: outside the
    envelope cutout with half-widths (c_x, c_y) = (n sx, n sy) scaled by
    `cut_scale_xy` -- a RECTANGLE |theta_x| > c_x OR |theta_y| > c_y, or
    an ELLIPSE (theta_x/c_x)^2 + (theta_y/c_y)^2 > 1 (the "circular pT
    cut" of the routing code when c_x = c_y).  The ePIC pots are planes
    around a horizontal SLOT (wide in x for the beam's momentum spread
    and dispersion, tight in y at the beam size; refs/README.md), i.e. a
    rectangle with c_x >> c_y: cut_scale_xy = (2.5, 1.0) emulates it.
    Returns dict with theta_x/y, pT, phi_t, t_reco = -pT^2 (x_L set to 1),
    the acceptance mask and the cut half-widths in GeV.
    """
    rng = rng or np.random.default_rng(20260824)
    Pp = np.asarray(Pprime, dtype=float)
    P = np.asarray(P, dtype=float)
    p_beam = np.sqrt((P[..., 1:] ** 2).sum(axis=-1))
    thx = Pp[..., 1] / Pp[..., 3]
    thy = Pp[..., 2] / Pp[..., 3]
    sx, sy = sigma_theta_xy
    cx, cy = n_sigma * sx * cut_scale_xy[0], n_sigma * sy * cut_scale_xy[1]
    n = thx.shape
    s_eff_x = np.hypot(sx, sigma_pos_over_l)
    s_eff_y = np.hypot(sy, sigma_pos_over_l)
    thx_m = thx + s_eff_x * rng.standard_normal(n)
    thy_m = thy + s_eff_y * rng.standard_normal(n)
    if shape == "rectangle":
        acc = (np.abs(thx_m) > cx) | (np.abs(thy_m) > cy)
    elif shape == "ellipse":
        acc = (thx_m / cx) ** 2 + (thy_m / cy) ** 2 > 1.0
    else:
        raise ValueError("shape must be 'rectangle' or 'ellipse'")
    pt = p_beam * np.hypot(thx_m, thy_m)
    return {"theta_x": thx_m, "theta_y": thy_m, "pT": pt,
            "phi_t": np.arctan2(thy_m, thx_m), "t_reco": -pt * pt,
            "accepted": acc, "pT_true": p_beam * np.hypot(thx, thy),
            "phi_t_true": np.arctan2(thy, thx),
            "cut_pt_xy": (float(np.mean(p_beam) * cx), float(np.mean(p_beam) * cy))}


def rp_hole_acceptance(slope_b, cut_x, cut_y, shape="rectangle", nphi=3600):
    """Azimuthal acceptance of an exponential coherent recoil spectrum
    dN/d2pT ~ exp(-B pT^2) outside a near-beam cutout with half-widths
    (cut_x, cut_y) in pT [GeV]: eps(phi) = exp(-B rho(phi)^2) with
    rho the cutout boundary along phi.  Returns
    {"phi", "eps", "acc", "a2", "a4"}: the total acceptance and the
    cos 2phi / cos 4phi Fourier coefficients <cos n phi> of the TAGGED
    sample -- the fake modulation a single-fill fit would attribute to
    physics (a2 = 0 only for cut_x = cut_y; a4 != 0 even then)."""
    phi = (np.arange(nphi) + 0.5) * 2.0 * np.pi / nphi
    c, s = np.abs(np.cos(phi)), np.abs(np.sin(phi))
    if shape == "rectangle":
        rho = np.minimum(cut_x / np.maximum(c, 1e-300),
                         cut_y / np.maximum(s, 1e-300))
    elif shape == "ellipse":
        rho = 1.0 / np.sqrt((c / cut_x) ** 2 + (s / cut_y) ** 2)
    else:
        raise ValueError("shape must be 'rectangle' or 'ellipse'")
    eps = np.exp(-slope_b * rho * rho)
    return {"phi": phi, "eps": eps, "acc": float(eps.mean()),
            "a2": float((eps * np.cos(2.0 * phi)).mean() / eps.mean()),
            "a4": float((eps * np.cos(4.0 * phi)).mean() / eps.mean())}
