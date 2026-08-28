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
  spin_state_ratio, harmonic_ratio_fit, err_harmonic_ratio,
  fill_acceptance_bias -- the one systematic the ratio cannot cancel
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

from polli_fastsim.farforward import HIGH_ACCEPTANCE, HIGH_DIVERGENCE
from polli_fastsim.kinematics import scattered_electron
from polli_fastsim.spectator import M_U

XING_IP6 = 25.0e-3   # rad, horizontal crossing angle at IP6
XING_IP8 = 35.0e-3   # rad, IP8

# --- the near-beam angular envelope ---------------------------------------
#
# CORRECTED 2026-08-27 (plans/10).  The two constants below are the LEGACY
# proton-derived placeholders: a single energy-independent, isotropic number
# back-derived from the documented Roman-Pot pT cuts for a 275 GeV proton
# (0.20 GeV high-acceptance, 0.41 GeV high-divergence) divided by
# 10 x 275 GeV.  Every published number before 2026-08-27 used them.
#
# The Yellow Report's own beam tables are neither energy-independent nor
# isotropic.  Table 10.1 (e+p) gives, per configuration and per optics, the
# RMS divergence h/v that `sigma_theta_for` now returns.  Use that.
SIGMA_THETA_HA = HIGH_ACCEPTANCE.sigma_theta      # LEGACY  73 microrad
# The high-divergence legacy value is the SAME constant as farforward's
# (0.45 GeV / (10 x 275 GeV) = 164 microrad) since 2026-08-28: this module
# used to derive it from the 0.41 GeV end of the same band (149 microrad),
# so the two packages disagreed by 10% for one optics (plans/10 A1b).  Both
# are placeholders that only the pre-2026-08-27 figures used; the published
# per-configuration values are `sigma_theta_for`.
SIGMA_THETA_HD = HIGH_DIVERGENCE.sigma_theta      # LEGACY 164 microrad


def sigma_theta_tagging(config, slope_b=50.0, n_sigma=10.0):
    """The divergence a TAGGING-OPTIMISED optics would have [rad, isotropic].

    Acceptance goes as exp(-C/beta*) while luminosity goes as 1/beta*, so
    the figure of merit L x acceptance is maximised where the n-sigma cut
    sits at t = 1/B, i.e. pT = 1/sqrt(B) (plans/10).  This returns the
    sigma_theta that puts it there, which is the working point a dedicated
    high-beta* lithium store would run at -- and the only working point at
    which the coherent channel exists at all once the corrected energies and
    divergences are used.  The beta* factor it implies relative to
    `sigma_theta_for` is the square of the ratio."""
    p_ion = config.ion.A * config.ion_momentum_per_nucleon
    return 1.0 / (slope_b ** 0.5 * n_sigma * p_ion)


def sigma_theta_for(config, optics="high-acceptance"):
    """(sigma_theta_h, sigma_theta_v) [rad] for a beam configuration: the
    Yellow Report Table 10.1 proton divergence of the configuration with
    the gamma-matched species step (plans/10).  Lives in
    polli_fastsim.farforward.sigma_theta_for since 2026-08-28 so that the
    fast simulation's spectator routing uses the same numbers; kept here
    by name for every existing caller."""
    from polli_fastsim import farforward as _ff
    return _ff.sigma_theta_for(config, optics)


# Far-forward transport constants (tools/fullsim): the horizontal lever
# from an IP angle to the pot-plane position and the dispersion at the
# pots.  Defined in polli_fastsim.farforward; re-exported here.
from polli_fastsim.farforward import POT_R12, POT_DISPERSION  # noqa: E402


def tagging_optics_point(config, slope_b=50.0, n_sigma=10.0, r_max=2000.0,
                         n_grid=400, dispersion=True, optics="high-acceptance"):
    """The lithium TAGGING OPTICS of Report 1 Section 6.1 -- see
    polli_fastsim.farforward.tagging_optics_point, of which this is the
    re-export used by the coherent reconstruction chain."""
    from polli_fastsim import farforward as _ff
    return _ff.tagging_optics_point(config, slope_b=slope_b, n_sigma=n_sigma,
                                    r_max=r_max, n_grid=n_grid,
                                    dispersion=dispersion, optics=optics)


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
    mass ion_mass (default: the physical nuclear mass of beams.NUCLEUS_MASS
    -- the same mass `beams` and `spectator` use, A * M_U being 12.6 MeV
    high for 6Li; pass 0.0 for the massless-target limit of the master
    formula).  Only O(gamma^2) terms and t_min feel the difference."""
    ion = config.ion
    p_a = ion.A * config.ion_momentum_per_nucleon
    if ion_mass is None:
        try:
            m_a = ion.A * ion.mass_per_nucleon
        except (KeyError, AttributeError):
            m_a = ion.A * M_U
    else:
        m_a = float(ion_mass)
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


# Yellow Report electromagnetic-calorimeter RESOLUTION REQUIREMENTS by
# pseudorapidity region (stochastic term, constant term); the constant
# term is quoted as a 1-3% band and the optimistic end is taken, as for
# the backward crystals.  The backward endcap (PbWO4) is the ePIC
# specification the whole chain used to apply everywhere; the barrel
# imaging calorimeter and the forward endcap are 3-6x coarser in the
# stochastic term (code review F4).
#: EM calorimeter (stochastic, constant) per eta region.  Fact-checked
#: 2026-08-27 against Yellow Report Table 3.1, Fig. 8.3 and Table 8.20:
#:   eta < -2.0   2%/sqrt(E) (+) 1%   -- correct, and matches the PbWO4
#:                                       crystal expectation
#:   -2.0..-1.0   7%/sqrt(E)          -- stochastic exact (Fig. 8.3);
#:                                       CONSTANT 2x OPTIMISTIC vs T3.1 (2%)
#:   -1.0..1.0    10%/sqrt(E)         -- the OPTIMISTIC CORNER of the YR's
#:                                       (10-12)%/sqrt(E) (+) (1-3)% band,
#:                                       and equal to the ePIC BIC
#:                                       requirement's stochastic term;
#:                                       CONSTANT 2-3x OPTIMISTIC vs both
#:   eta > 1.0    7%/sqrt(E)          -- MIS-SOURCED.  7%/sqrt(E) belongs to
#:                                       -2.0 < eta < -1.0 only; neither YR
#:                                       table gives it for the forward
#:                                       endcap.  Optimistic by ~1.5x
#:                                       stochastic and 2x constant.
#: Also absent: the YR's 1-2%/E NOISE term, which at 1 GeV is 2% -- as large
#: as the whole backward resolution quoted here.
#:
#: WHY IT SURVIVES ANYWAY.  Every sweet spot this programme publishes puts
#: the scattered electron BACKWARD (eta = -2.9 to -1.6), where the table is
#: right to within 20%, and `RecoModel(energy="best")` hands those events to
#: the tracker in any case.  The forward row is latent, not active.
EMCAL_YR_TABLE = ((-np.inf, -2.0, 0.02, 0.01),     # backward endcap, PbWO4
                  (-2.0, -1.0, 0.07, 0.01),        # backward transition
                  (-1.0, 1.0, 0.10, 0.01),         # barrel, optimistic corner
                  (1.0, np.inf, 0.07, 0.01))       # forward endcap, mis-sourced

#: Beyond this the ePIC tracker and calorimeters do not reach (-4.0 < eta <
#: 3.7 in the current design; the Yellow Report calls |eta| > 4.0 "not
#: accessible" and 3.5-4.0 "reduced performance").  The tables above extend
#: to +-inf, which is latent rather than active here -- nothing in this
#: programme scatters an electron past |eta| = 3.5 -- but a consumer that
#: does should clip on this.
ETA_ACCEPTANCE_MAX = 3.5


def emcal_resolution(e_prime, stoch=0.02, const=0.01, noise=0.0, eta=None):
    """dE/E of an EM calorimeter: stoch/sqrt(E) (+) const (+) noise/E.

    With `eta` the stochastic and constant terms come from the Yellow
    Report requirement table per region (EMCAL_YR_TABLE) and `stoch`,
    `const` are ignored; without it the PbWO4-class defaults (the ePIC
    backward-endcap specification, eta < -2) apply at every eta, which is
    what the chain did before 2026-08-25 and is optimistic by 3-5x in the
    stochastic term for barrel electrons (code review F4).  The
    sweet-spot electrons are backward (eta = -2.9 to -1.6), so the
    headline numbers are insensitive; the amplitude-vs-x panels reach the
    barrel, where the tracker is the better measurement and
    RecoModel(energy="best") selects it."""
    e = np.maximum(np.asarray(e_prime, dtype=float), 1e-9)
    if eta is not None:
        eta = np.asarray(eta, dtype=float)
        stoch = np.full(eta.shape, EMCAL_YR_TABLE[0][2])
        const = np.full(eta.shape, EMCAL_YR_TABLE[0][3])
        for lo, hi, sv, cv in EMCAL_YR_TABLE:
            m = (eta >= lo) & (eta < hi)
            stoch = np.where(m, sv, stoch)
            const = np.where(m, cv, const)
    return np.sqrt((stoch / np.sqrt(e)) ** 2 + const ** 2 + (noise / e) ** 2)


def tracking_resolution(e_prime, eta, a=None, b=None):
    """dp/p = sqrt((a p)^2 + b^2), eta-piecewise.

    CORRECTED 2026-08-27: this table is NOT a placeholder.  It is the Yellow
    Report tracking REQUIREMENT -- Fig. 8.3 / Table 11.2 (arXiv:2103.05419
    pp. 261, 437) -- which ATHENA Table 4 also quotes as "Requirements".  The
    barrel row (a, b) = (0.0005, 0.005) and the 1.0 < eta < 2.5 row
    (0.0005, 0.010) match it exactly; the endcap constant terms are padded
    (0.030 against 0.020 backward, 0.025 against 0.020 forward), i.e.
    conservative.  An earlier docstring called the whole table unsourced.

    Three cautions, none of which change the headline numbers:
      * it is a REQUIREMENT, not an achievement, so agreement with it is not
        validation -- comparing this table to the requirement compares a
        number to itself.  ePIC full simulation reaches 0.45-0.6% at
        p = 1 GeV/c depending on eta slice, so the barrel smearing is
        realistic to about +-20%.
      * the parameterisation is in TOTAL momentum p, not pT.  The consuming
        code passes p (money_delta_20260729.py sets p_true = e_prime_true),
        which is right; quoting the same numbers "at pT = 1 GeV/c" is not.
      * the table extends past |eta| = 3.5, where the tracker ends
        (-4.0 < eta < 3.7 in the current design and "not accessible" above
        4.0 in the YR).  Nothing in this programme scatters an electron
        there, so it is latent rather than active."""
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
    direction resolution used for phi), eta-piecewise: 5/3/2/1/2/3/5 mrad.

    A genuine PLACEHOLDER -- no angular-resolution requirement exists in
    either Yellow Report table (Fig. 8.3's only angular entry is the low-Q2
    tagger's dtheta/theta < 1.5%), and no ePIC full-simulation angular
    resolution for backward electrons has been published.  It sets the
    electron-method Q2 resolution at the low-y sweet spots
    (cot(theta'/2) dtheta' = 5% at theta' = 0.1 rad for 3 mrad).

    Two things bound it.  The 0.1 mrad quoted in some ePIC talks as a
    smearing input is not a measured resolution and should not replace
    this table.  And the chain carries NO separate beam-divergence term on
    theta' (`beam_e_spread` enters only dy/y): the electron beam's
    angular divergence at the IP, which the Yellow Report puts at
    81-211 microrad h/v depending on energy and optics and states "cannot
    be corrected on an event-by-event basis", is an irreducible floor that
    this table is the only stand-in for.  The defensible bracket is
    therefore between that floor (0.08-0.21 mrad) and the 3 mrad used
    here -- a factor <~ 3 on dQ2/Q2 at the sweet spots (5.2% -> ~2%,
    purity 0.64 -> 0.70 at spot 1 with the table divided by three, code
    review 2026-08-25 F3), not a factor 30."""
    eta = np.asarray(eta, dtype=float)
    out = np.full_like(eta, 1.0e-3)
    for lo, hi, val in ((-np.inf, -3.5, 5.0e-3), (-3.5, -2.5, 3.0e-3),
                        (-2.5, -1.0, 2.0e-3), (1.0, 2.5, 2.0e-3),
                        (2.5, 3.5, 3.0e-3), (3.5, np.inf, 5.0e-3)):
        out = np.where((eta >= lo) & (eta < hi), val, out)
    return out


def eps_eid(eta):
    """Electron-ID efficiency eps_eID(eta); zero outside |eta| > 3.5.

    Re-documented 2026-08-27.  This is a CONSTRUCTED eta profile, not a
    published curve, and it is anchored at exactly one point.  No ePIC
    electron-ID efficiency curve exists in any ePIC document.

    The two proposal-era sources disagree with each other in the barrel by
    about 25 points:
      * ATHENA, JINST 17 (2022) P10019, Table 5 (printed p. 20) gives ONE
        number -- 95% electron efficiency at >99.8% pion rejection for the
        BARREL ECal, at all p >= 0.1 GeV/c, with NO eta dependence, and from
        a standalone calorimeter simulation with no material in front and no
        magnetic field.
      * ECCE, NIM A 1055 (2023) 168464 Sec. 3.5.2 IS eta-resolved: EEMC
        (-3.4 < eta < -1.5) ~95%, BEMC (-1.72 < eta < 1.31) ~70%, FEMC
        (1.3 < eta < 3.5) ~90-95%.
    This curve's 0.95 at eta = -2.0 is ECCE's EEMC value; its barrel 0.90
    matches neither source, and its forward tail (0.85/0.80/0.70) runs
    OPPOSITE to ECCE's FEMC.

    IT ALMOST CERTAINLY DOES NOT MATTER.  The tensor observable is a
    spin-state RATIO (reco.spin_state_ratio), so an eta-dependent but
    fill-independent efficiency cancels exactly, bin by bin; eps_eID should
    enter the statistical error and nothing else.  Check that before
    spending effort on the curve."""
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
    expectation; 0.25 (the ePIC study's own value) is the default of
    recopseudo.RecoModel and of the R-scripts."""
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
    eps_i COMMON TO ALL FILLS cancels bin by bin (a fill-dependent
    harmonic enters as (eps2_f - eps2_f')/(P_f - P_f'): bunch-by-bunch
    alternation is required, code review 2026-08-25), and a
    relative-luminosity error enters
    only through Pbar, i.e. as a bin-INDEPENDENT offset of R plus a
    second-order rescaling (delta x Pbar/sigma_P^2).  Returns
    (R, var_R, sigma_P2, Pbar) with var_R from linear error propagation
    of Poisson counts,

        var(R) = sum_f ((w_f - R)/sum N)^2 E[N_f],

    with the per-fill counts entering through their EXPECTED values
    E[N_f] = l_f sum N (1 + P_f T)/(1 + Pbar T), T = R/sigma_P^2 to first
    order, rather than the observed N_f.  The two agree in expectation, but
    the observed form is degenerate for a bin populated by ONE fill only
    (there R = w_f exactly, the (w_f - R)^2 N_f term vanishes and the other
    fill contributes nothing, so var = 0 and the bin gets an infinite
    weight in the fit -- which is what happened in low-count bins at the
    edge of a Roman-Pot cutout, code review 2026-08-28).
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
    # first-order T, clipped so that every expected count stays positive
    # (a single-fill bin has R = w_f, far outside the physical |T| << 1)
    t_max = 0.5 / max(float(np.max(np.abs(p))), 1e-12)
    t_lin = np.clip(r / sig2, -t_max, t_max)
    n_exp = (lw[:, None] * den_safe[None, :] * (1.0 + p[:, None] * t_lin[None, :])
             / (1.0 + pbar * t_lin)[None, :])
    var = np.where(live, (((w - r[None, :]) / den_safe[None, :]) ** 2
                          * n_exp).sum(axis=0), np.inf)
    return r, var, sig2, pbar


def _ratio_to_modulation(r, var, sig2, pbar, u=0.0, n_iter=4):
    """Invert R = sigma_P^2 T / (1 + u + Pbar T) for T bin by bin
    (fixed-point iteration; u = known spin-independent modulation of the
    denominator, 0 for inclusive DIS), propagating the variance."""
    t = r / sig2
    for _ in range(n_iter):
        t = r * (1.0 + u + pbar * t) / sig2
    # exact Jacobian of the inversion: T = R (1 + u)/(sigma_P^2 - Pbar R),
    # so dT/dR = (1 + u + Pbar T)/(sigma_P^2 - Pbar R); the earlier
    # (1 + u + Pbar T)/sigma_P^2 dropped the second-order denominator
    # (0.3% on the error at |Pbar T| ~ 1e-3, code review 2026-08-28)
    scale = (1.0 + u + pbar * t) / (sig2 - pbar * r)
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
    "u1": <cos(a-b)>, "u2": <cos 2(a-b)>} and the three PARITY-FORBIDDEN
    partners {"e_s": <sin 2a>, "t_s": <sin 2b>, "m_s": <sin(a+b)>}.  The
    alpha average is analytic (uniform in-bin distribution).  The beta
    averages are analytic too unless `beta_means` =
    {"c1","s1","c2","s2"} (Kb,) supplies the ACCEPTANCE-WEIGHTED in-bin
    means <cos b>, <sin b>, <cos 2b>, <sin 2b> of the true beta of the
    events reconstructed into each beta bin -- required when the
    acceptance varies strongly across a bin (the Roman-Pot cutout: x25
    across beta), and the way the smearing of beta enters the response
    (the fit then estimates the UNSMEARED coefficients, as an
    MC-corrected analysis does).

    Why the sin columns are a null test: with unpolarized leptons and a
    HEADLESS alignment axis the cross section is even under
    (alpha, beta) -> (-alpha, -beta), so sin 2alpha, sin 2beta and
    sin(alpha+beta) are exactly forbidden.  A spin-axis azimuth error
    delta makes all three ratios (tan 2d, tan 2d, tan d); a Roman-Pot
    azimuthal ROLL d_t gives (0, tan 2d_t, tan d_t), which separates
    them.  `t_s` uses the same t-weighted template as `t` (beta_means
    "s2t"): with the plain <sin 2b> the null test misses its own closure
    by ~7%, exactly as `c2t` is needed for `t`."""
    ae = np.asarray(alpha_edges, dtype=float)
    be = np.asarray(beta_edges, dtype=float)
    ac = 0.5 * (ae[:-1] + ae[1:])
    bc = 0.5 * (be[:-1] + be[1:])
    d1a, d1b, d2a, d2b = _bin_dilutions(ae, be)
    if beta_means is None:
        c1, s1 = np.cos(bc) * d1b, np.sin(bc) * d1b
        c2, s2 = np.cos(2 * bc) * d2b, np.sin(2 * bc) * d2b
        c2t, s2t = c2, s2
    else:
        c1, s1, c2, s2 = (np.asarray(beta_means[k], dtype=float)
                          for k in ("c1", "s1", "c2", "s2"))
        # template basis of a t-dependent a_t (recopseudo.basis_means)
        c2t = np.asarray(beta_means.get("c2t", c2), dtype=float)
        s2t = np.asarray(beta_means.get("s2t", s2), dtype=float)
    ca1, sa1 = (np.cos(ac) * d1a)[:, None], (np.sin(ac) * d1a)[:, None]
    ca2, sa2 = (np.cos(2 * ac) * d2a)[:, None], (np.sin(2 * ac) * d2a)[:, None]
    ones_a = np.ones((ac.size, 1))
    return {"e": (ca2 * np.ones((1, bc.size))).ravel(),
            "t": (ones_a * c2t[None, :]).ravel(),
            "m": (ca1 * c1[None, :] - sa1 * s1[None, :]).ravel(),
            "u1": (ca1 * c1[None, :] + sa1 * s1[None, :]).ravel(),
            "u2": (ca2 * c2[None, :] + sa2 * s2[None, :]).ravel(),
            "e_s": (sa2 * np.ones((1, bc.size))).ravel(),
            "t_s": (ones_a * s2t[None, :]).ravel(),
            "m_s": (sa1 * c1[None, :] + ca1 * s1[None, :]).ravel()}


def unpolarized_modulation_2d(alpha_edges, beta_edges, u1, u2,
                              beta_means=None):
    """Bin-averaged spin-independent modulation u(alpha - beta) =
    u1 cos(alpha-beta) + u2 cos 2(alpha-beta), flattened in C order."""
    b = basis_2d(alpha_edges, beta_edges, beta_means)
    return u1 * b["u1"] + u2 * b["u2"]


def _harmonic_rank_guard(a, cols, live, nbins, label="(alpha, beta)"):
    """Raise a NAMED LinAlgError when the weighted harmonic design has
    lost rank.  A cutout tight enough to empty whole regions of
    (alpha, beta) leaves the template columns linearly dependent, and
    `np.linalg.inv`/`solve` then raises a bare "Singular matrix" that says
    nothing about why.  The cause is a statement about the acceptance, not
    a numerical accident.  (Hit for real when the measured Roman-Pot
    aperture was put under the mid-energy configuration, where it leaves
    2.7e-3 of the recoils.)  `a` is the design weighted by the square root
    of the per-bin weight -- 1/var for the ratio fit, the diagonal of the
    expected information for the likelihood -- so the two estimators fail
    on the same designs and report the same cause (they must: the Fisher
    information of the profiled likelihood is the same Gram matrix as the
    weighted design, up to weights)."""
    if np.linalg.matrix_rank(a) < a.shape[1]:
        raise np.linalg.LinAlgError(
            "the %d-column harmonic design is rank-deficient: %d of %d "
            "%s bins carry weight, and the cutout has emptied "
            "enough of the circle that %s can no longer be separated. "
            "Loosen the cutout, widen the bins, or drop columns "
            "(with_sin=False)."
            % (a.shape[1], live, nbins, label, "/".join(cols)))


def harmonic_ratio_fit_2d(counts, lumis, pzz, alpha_edges, beta_edges,
                          u_coeffs=None, beta_means=None, with_sin=False):
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
    With `with_sin` the three parity-forbidden partners sin 2alpha,
    sin 2beta and sin(alpha+beta) are fitted alongside as a null test
    (basis_2d): all three vanish for an aligned axis and an unrolled pot,
    and their pattern identifies which of the two is misaligned.

    Returns {"const", "a_e", "a_t", "a_m", "err_e", "err_t", "err_m",
    "cov", "sigma_p2", "pbar"} (+ "a_e_s", "a_t_s", "a_m_s" and their
    errors with with_sin)."""
    n = np.asarray(counts, dtype=float)
    nf, ka, kb = n.shape
    r, var, sig2, pbar = spin_state_ratio(n.reshape(nf, ka * kb), lumis, pzz)
    basis = basis_2d(alpha_edges, beta_edges, beta_means)
    u = (u_coeffs[0] * basis["u1"] + u_coeffs[1] * basis["u2"]
         if u_coeffs is not None else 0.0)
    t, var_t = _ratio_to_modulation(r, var, sig2, pbar, u=u)
    cols = ["e", "t", "m"] + (["e_s", "t_s", "m_s"] if with_sin else [])
    design = np.vstack([np.ones(ka * kb)] + [basis[c] for c in cols]).T
    wgt = 1.0 / np.sqrt(np.maximum(var_t, 1e-300))
    a = design * wgt[:, None]
    live = int(np.count_nonzero(np.isfinite(var_t) & (var_t > 0.0)))
    _harmonic_rank_guard(a, cols, live, ka * kb)
    coef, *_ = np.linalg.lstsq(a, t * wgt, rcond=None)
    cov = np.linalg.inv(a.T @ a)
    err = np.sqrt(np.diag(cov))
    # err[0] belongs with coef[0]: the two estimators advertise the same
    # return keys, and until 2026-08-28 only the likelihood twin carried
    # "err_const", so code written to that contract raised KeyError on the
    # default (published) estimator
    out = {"const": coef[0], "err_const": err[0], "cov": cov,
           "sigma_p2": sig2, "pbar": pbar}
    for i, c in enumerate(cols, start=1):
        out["a_" + c] = coef[i]
        out["err_" + c] = err[i]
    return out


def _profile_likelihood_newton(n, lw, p, design, u, max_iter=50, tol=1e-10):
    """Newton iteration of the acceptance-profiled Poisson log-likelihood
    of `harmonic_likelihood_fit_2d`, from A = 0, with a step-halving guard
    that keeps every density q_{f,b} positive in the bins that carry
    counts.  The step uses the OBSERVED information (the exact Hessian of
    the profile), falling back on the expected information wherever a
    fluctuation makes that Hessian indefinite or singular -- both have the
    same zero, so the fixed point is unchanged.  `n` is (F, B), `design`
    (B, npar), `u` (B,).  Returns (A, n_iter, grad, q, D): the
    coefficients, the iterations used, the gradient AT the returned A (a
    convergence diagnostic, and it is recomputed after the loop for that
    reason), and the per-fill and spin-averaged densities at the solution.

    The loop exits on the step, so it cannot by itself tell convergence
    from a step-halving STALL: where the likelihood has no interior
    maximum the guard drives lambda to 2^-60 and lambda*step underflows
    below `tol` at a point whose gradient is huge.  `_check_profile_
    convergence` is what separates the two, on the gradient this returns;
    every caller runs it (2026-08-28 review, finding 5)."""
    pbar = float((lw * p).sum())
    theta = np.zeros(design.shape[1])
    live = n.sum(axis=0) > 0.0
    design_live, u_live = design[live], u[live]   # hoisted: the positivity
    grad = np.zeros(design.shape[1])              # guard runs every step
    it = 0
    for it in range(1, int(max_iter) + 1):
        t = design @ theta
        q = 1.0 + u[None, :] + p[:, None] * t[None, :]
        d = 1.0 + u + pbar * t
        score = (n * (p[:, None] / q - pbar / d[None, :])).sum(axis=0)
        grad = design.T @ score
        curv = (n * (p[:, None] ** 2 / q ** 2
                     - pbar ** 2 / d[None, :] ** 2)).sum(axis=0)
        info = design.T @ (curv[:, None] * design)
        try:
            np.linalg.cholesky(info)             # positive definite?
        except np.linalg.LinAlgError:
            w = _expected_information_weights(n, lw, p, pbar, q, d)
            info = design.T @ (w[:, None] * design)
        step = np.linalg.solve(info, grad)
        lam = 1.0
        for _ in range(60):
            qc = (1.0 + u_live
                  + p[:, None] * (design_live @ (theta + lam * step))[None, :])
            if np.all(qc > 1e-6):
                break
            lam *= 0.5
        theta = theta + lam * step
        if np.max(np.abs(lam * step)) < tol:
            break
    t = design @ theta
    q = 1.0 + u[None, :] + p[:, None] * t[None, :]
    d = 1.0 + u + pbar * t
    grad = design.T @ (n * (p[:, None] / q - pbar / d[None, :])).sum(axis=0)
    return theta, it, grad, q, d


def _check_profile_convergence(grad, covm, label="(alpha, beta)", tol=1e-3):
    """Raise a NAMED LinAlgError when the profiled likelihood did not
    reach its stationary point.  The Newton loop stops on the step, and a
    step the positivity guard has halved sixty times is below any
    tolerance whatever the gradient is; the fit would then return a point
    on the boundary 1 + u_b + P_f T_b = 0 with a spuriously tiny
    covariance and nothing to say so.  The test is on the gradient at the
    returned solution measured in units of the error it would move --
    |dlnL/dA_c| sigma_c, dimensionless -- so it does not depend on the
    luminosity or on the number of bins.  The margin is enormous where the
    estimator is used -- 6e-15 over 200 Poisson draws of the sparsest
    published |t| bin at its 48 counts per (alpha, beta) cell, and 2e-15
    at one count per cell, against a tolerance of 1e-3 -- and the guard
    first fires around a tenth of a count per cell, an occupancy no
    repository command reaches."""
    # |diag| so that this stays a test of the GRADIENT: with cov=
    # "observed" a fluctuation can make a variance negative, which is a
    # different complaint and surfaces as a nan error bar on its own
    scaled = (np.abs(np.asarray(grad, dtype=float))
              * np.sqrt(np.abs(np.diag(covm))))
    worst = float(np.max(scaled)) if scaled.size else 0.0
    if not np.isfinite(worst) or worst > tol:
        raise np.linalg.LinAlgError(
            "the profiled likelihood stalled on the positivity boundary "
            "of the %s design: the step-halving guard drove the Newton "
            "step to zero at a point whose gradient is still %.3g error "
            "bars (tolerance %g), so the returned coefficients and their "
            "covariance are meaningless.  That happens when a spin state "
            "has no counts where the model needs them -- an empty fill, a "
            "cutout that empties whole regions of the circle -- not at "
            "any occupancy this chain reaches." % (label, worst, tol))


def _expected_information_weights(n, lw, p, pbar, q, d):
    """Diagonal of the CONDITIONAL Fisher information at the observed bin
    totals: n_b sum_f p_{f,b} (P_f/q_{f,b} - Pbar/D_b)^2 with
    p_{f,b} = l_f q_{f,b}/D_b the conditional multinomial shares."""
    s = p[:, None] / q - pbar / d[None, :]
    return n.sum(axis=0) * (lw[:, None] * q / d[None, :] * s * s).sum(axis=0)


def harmonic_likelihood_fit_2d(counts, lumis, pzz, alpha_edges, beta_edges,
                               u_coeffs=None, beta_means=None, with_sin=False,
                               cov="expected", max_iter=50, tol=1e-10):
    """Maximum-likelihood twin of `harmonic_ratio_fit_2d`: the same model,
    the same basis (`basis_2d` with the same arguments, so the c2t/s2t
    templates and the bin dilutions are bit-for-bit shared), the same
    return keys -- and no low-count bias.

    THE MODEL.  Spin-sorted counts in the (Ka x Kb) bins of
    alpha = phi_e - phi_S and beta = phi_t - phi_S are Poisson with

        E[N_{f,b}] = l_f eps_b [1 + u_b + P_f (H A)_b],                (1)
        (H A)_b = kappa + a_e <cos 2a> + a_t <g(t) cos 2b>
                  + a_m <cos(a+b)>  [+ the three sin partners],

    which is EXACTLY what `recopseudo.CoherentResponse.expected_counts_2d`
    integrates: eps_b (acceptance x flux) is common to the fills, only the
    luminosity share l_f and the tensor polarization P_f are
    fill-dependent.  That common eps_b is what the spin-state ratio
    cancels bin by bin; here it is carried as a free nuisance parameter
    per bin and PROFILED OUT.

    WHY THE PROFILE IS EXACT.  eps_b enters (1) as a linear scale, so

        d lnL/d eps_b = 0  =>  eps_b(A) = n_b / D_b(A),
        n_b = sum_f N_{f,b},  D_b = sum_f l_f q_{f,b} = 1 + u_b + Pbar (H A)_b,

    and substituting it back leaves everything that depends on n_b alone
    outside the fit: the profile likelihood is exactly the CONDITIONAL
    multinomial given the bin totals,

        -ln L_prof(A) = - sum_{b,f} N_{f,b} ln p_{f,b},
        p_{f,b} = l_f q_{f,b} / D_b,  q_{f,b} = 1 + u_b + P_f (H A)_b.  (6)

    Profiling B nuisance parameters that grow with the data would normally
    raise the Neyman-Scott incidental-parameter worry; it does not here,
    because the conditional score has zero mean bin by bin at ANY count:

        E[d lnL/dA_c | n_b] = n_b sum_f p_{f,b} h_c(b)
                              [P_f/q_{f,b} - Pbar/D_b] = 0              (7)

    since sum_f p_f (P_f/q_f - Pbar/D) = Pbar/D - Pbar/D.  There is no
    1/nu_b term anywhere.  That is the whole point: the ratio estimator
    inverts R_b = sigma_P^2 T_b/(1 + u_b + Pbar T_b) bin by bin, and the
    inversion T = R(1+u)/(sigma_P^2 - Pbar R) is strongly CURVED over the
    range |R| <= max_f |P_f - Pbar| that R actually explores at low counts
    (a bin with one count has R = w_f exactly): d2T/dR2 = 2(1 + u)
    sigma_P^2 Pbar/(sigma_P^2 - Pbar R)^3 carries the sign of Pbar, so the
    flip plan's P0 = -2 P+ (Pbar = -0.3) makes it strictly CONCAVE and the
    Jensen offset (1 + u_b) Pbar/(sigma_P^2 nu_b) per bin negative, while
    the data-driven 1/var weights add an opposite-sign term of the same
    order.  The offset is flat in alpha (eps_b is: the Roman-Pot cutout
    modulates beta, not alpha) so the LSQ lands it on the constant and on
    a_t and leaves a_e alone.  Measured on the real chain at the tagging
    optics over two hundred one-year pseudo-experiments (5 x 40.8, 12 x 24
    bins, money_cos2phi_coherent_reco.py --config 0 --optics tagging
    --n-mc 6000000 --ensemble 200): the ratio's a_t is biased by
    +0.1 / -0.2 / -4.0 / -34.3 % in the four published |t| bins (1451,
    620, 192, 48 counts per (alpha, beta) bin), this fit by
    +0.4 / +0.3 / +0.3 / +0.7 %, with pulls of the mean
    +2.0 / +1.0 / +0.5 / +0.9 against the ratio's +0.4 / -0.7 / -9.0 /
    -42.1 -- and the first bin's residual is the response Monte Carlo's
    own floor, not the estimator.  (Table 5 of Report 2 quotes the same
    comparison on the TWENTY draws its published columns use, where the
    ratio reads 0.0 / -0.9 / -5.4 / -37 %; the two ensembles are not
    mixed anywhere.)  Push the sparsest bin to ONE count per
    (alpha, beta) cell and the ratio changes sign, mean -0.085 against
    +0.182 injected, while this fit still returns +0.181.  Its quoted
    errors are the ensemble spread at both occupancies, where the ratio's
    are 62% larger than its own spread at one count per bin (the
    compressed estimator fluctuates less than it says); and it is the
    smaller of the two even in variance, 0.0200 against 0.0211 in the
    sparsest published bin.

    Empty bins contribute exactly nothing to (6) -- no infinite weight, no
    ad-hoc masking, no selection on n_b > 0 -- and a bin populated by one
    fill contributes a finite, correct n ln p term, which is the
    single-fill pathology of `spin_state_ratio` removed structurally
    rather than patched.  Coarser (alpha, beta) bins only ATTENUATE the
    ratio's bias -- in the sparsest chain bin, over the same two hundred
    experiments, -34.3% at 12 x 24, -13.9% at 8 x 16 and -7.8% at 6 x 12
    (--n-alpha/--n-beta) -- cost 17% on err(a_e) through the wider-bin
    dilutions (0.0139 -> 0.0148 -> 0.0162), and run into a hard floor at
    Ka = 4, where <cos 2alpha> vanishes identically and the design is
    rank-deficient at 430 counts per bin: adaptive binning is a
    cross-check, not the fix.

    WHAT IT SHARES WITH THE RATIO, AND THEREFORE DOES NOT FIX.  Model (1)
    assumes eps_b COMMON to the fills, so this estimator is blind to a
    fill-dependent acceptance in exactly the same way, and a wrong
    (u1, u2) or a wrong luminosity share moves it the same way too: on
    exact counts the shifts agree to better than 1% of themselves and to
    2% of the statistical error, the residual being second order in a
    deliberately large 12% perturbation.  On exact (Asimov) counts the
    score is identically zero at the truth, so it returns the injected
    coefficients to machine precision (relative 4e-16 against the ratio's
    2e-6) and its errors agree with the ratio fit's to 0.24% on err_t and
    0.01% on err_e, INDEPENDENT of the luminosity -- no published error
    bar moves.  The residual is the difference between the delta-method
    variance of the inverted ratio and the exact Fisher information, and
    it goes the efficient way: this fit's error is the smaller.

    ARGUMENTS.  `counts` (F, Ka, Kb); `lumis`, `pzz` per fill;
    `u_coeffs` = (u1, u2) the spin-independent modulation the analysis
    subtracts (None = 0); `beta_means` the response's acceptance-weighted
    in-bin means (`CoherentResponse.basis_means`); `with_sin` adds the
    three parity-forbidden null columns.  `cov` selects the covariance:
    "expected" (default) is the conditional Fisher information at the
    observed bin totals, "observed" the Hessian at the solution; over
    Poisson draws of the sparsest published bin the two agree to 0.01% on
    average at 48 counts per (alpha, beta) cell (0.5% worst case) and to
    0.5% at one count per cell, where the observed form can fluctuate by
    14% -- which is why the expected one is the default.

    Returns {"const", "a_e", "a_t", "a_m", "err_const", "err_e", "err_t",
    "err_m", "cov", "sigma_p2", "pbar", "n_iter", "nll", "grad_max"} (+
    "a_e_s", "a_t_s", "a_m_s" and their errors with `with_sin`), with
    `cov` ordered (const, e, t, m [, e_s, t_s, m_s]) as in the ratio fit.
    `nll` is the conditional -ln L of (6) at the solution and `grad_max`
    the largest |dlnL/dA_c| there, which `_check_profile_convergence` has
    already tested against the error bars: a fit that stalled on the
    positivity boundary raises a named LinAlgError rather than returning
    a boundary point with a meaningless covariance.  Cost 2.3 ms at
    12 x 24 with seven columns against 0.4 ms for the ratio fit -- five
    times as much, and negligible either way against building the
    response."""
    n = np.asarray(counts, dtype=float)
    nf, ka, kb = n.shape
    n = n.reshape(nf, ka * kb)
    p = np.asarray(pzz, dtype=float)
    lum = np.asarray(lumis, dtype=float)
    lw = lum / lum.sum()
    pbar = float((lw * p).sum())
    sig2 = float((lw * (p - pbar) ** 2).sum())
    basis = basis_2d(alpha_edges, beta_edges, beta_means)
    u = (u_coeffs[0] * basis["u1"] + u_coeffs[1] * basis["u2"]
         if u_coeffs is not None else np.zeros(ka * kb))
    u = np.broadcast_to(np.asarray(u, dtype=float), (ka * kb,)).copy()
    cols = ["e", "t", "m"] + (["e_s", "t_s", "m_s"] if with_sin else [])
    design = np.vstack([np.ones(ka * kb)] + [basis[c] for c in cols]).T
    # the rank guard, on the same footing as the ratio fit's: the weight a
    # bin carries in the Fisher information at A = 0 is
    # n_b sigma_P^2/(1 + u_b)^2, so the live bins are exactly those with
    # counts and the two estimators fail on the same designs
    nb = n.sum(axis=0)
    w0 = nb * sig2 / (1.0 + u) ** 2
    _harmonic_rank_guard(design * np.sqrt(np.maximum(w0, 0.0))[:, None],
                         cols, int(np.count_nonzero(nb > 0.0)), ka * kb)
    theta, n_iter, grad, q, d = _profile_likelihood_newton(
        n, lw, p, design, u, max_iter=max_iter, tol=tol)
    if cov == "observed":
        w = (n * (p[:, None] ** 2 / q ** 2
                  - pbar ** 2 / d[None, :] ** 2)).sum(axis=0)
    else:
        w = _expected_information_weights(n, lw, p, pbar, q, d)
    covm = np.linalg.inv(design.T @ (w[:, None] * design))
    _check_profile_convergence(grad, covm)
    err = np.sqrt(np.diag(covm))
    nll = -float((n * np.log(lw[:, None] * q / d[None, :],
                             where=n > 0.0,
                             out=np.zeros_like(q))).sum())
    out = {"const": theta[0], "err_const": err[0], "cov": covm,
           "sigma_p2": sig2, "pbar": pbar, "n_iter": n_iter, "nll": nll,
           "grad_max": float(np.max(np.abs(grad)))}
    for i, c in enumerate(cols, start=1):
        out["a_" + c] = theta[i]
        out["err_" + c] = err[i]
    return out


def harmonic_likelihood_fit(counts, lumis, pzz, edges, with_sin=False):
    """One-azimuth twin of `harmonic_ratio_fit`, on the same Newton core
    as `harmonic_likelihood_fit_2d`: the per-bin efficiency is profiled
    out, which makes the estimator exactly the conditional multinomial
    given the bin totals and therefore unbiased at any count.  The
    inclusive phi' bins carry >= 10^4 counts, where it agrees with the
    ratio fit to a small fraction of the error, so this is the general
    statement rather than a change of the published inclusive numbers.

    Returns {"amp", "err", "const", "err_const", "sigma_p2", "pbar",
    "n_iter", "nll"} (+ "amp_sin", "phase" with `with_sin`), with the
    finite-bin dilution sin(2w)/2w divided out of the amplitude and its
    error exactly as in `harmonic_ratio_fit`."""
    n = np.asarray(counts, dtype=float)
    p = np.asarray(pzz, dtype=float)
    lw = np.asarray(lumis, dtype=float)
    lw = lw / lw.sum()
    pbar = float((lw * p).sum())
    sig2 = float((lw * (p - pbar) ** 2).sum())
    edges = np.asarray(edges, dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    cols = [np.ones_like(centers), np.cos(2.0 * centers)]
    if with_sin:
        cols.append(np.sin(2.0 * centers))
    design = np.vstack(cols).T
    u = np.zeros(centers.size)
    nb = n.sum(axis=0)
    _harmonic_rank_guard(design * np.sqrt(nb * sig2)[:, None],
                         ["cos2phi"] + (["sin2phi"] if with_sin else []),
                         int(np.count_nonzero(nb > 0.0)), centers.size,
                         label="phi'")
    theta, n_iter, grad, q, d = _profile_likelihood_newton(n, lw, p,
                                                           design, u)
    w = _expected_information_weights(n, lw, p, pbar, q, d)
    covm = np.linalg.inv(design.T @ (w[:, None] * design))
    _check_profile_convergence(grad, covm, label="phi'")
    err = np.sqrt(np.diag(covm))
    wdt = 0.5 * (edges[1] - edges[0])
    dil = np.sin(2.0 * wdt) / (2.0 * wdt)
    nll = -float((n * np.log(lw[:, None] * q / d[None, :], where=n > 0.0,
                             out=np.zeros_like(q))).sum())
    out = {"amp": theta[1] / dil, "err": err[1] / dil, "const": theta[0],
           "err_const": err[0], "sigma_p2": sig2, "pbar": pbar,
           "n_iter": n_iter, "nll": nll}
    if with_sin:
        out["amp_sin"] = theta[2] / dil
        out["phase"] = 0.5 * np.arctan2(theta[2], theta[1])
    return out


def unpolarized_insitu_fit_2d(counts, eps_mc, lumis, pzz, alpha_edges,
                              beta_edges, harmonics=None, beta_means=None,
                              with_sin=False, max_iter=50, tol=1e-12):
    """Fit the spin-independent modulation (u1, u2) IN SITU, from the
    spin-averaged counts of the same data.

    THE ASSUMPTION, STATED.  Once the per-bin acceptance eps_b is a free
    parameter -- which is what makes both the ratio and
    `harmonic_likelihood_fit_2d` acceptance-free -- u is NOT identifiable:
    u_b multiplies eps_b in exactly the same way, and the spin-sorted
    counts constrain only their product.  An in-situ (u1, u2) therefore
    NECESSARILY leans on the acceptance model, and this fit says so
    explicitly: it compares the spin-averaged counts with the response's
    own acceptance prediction,

        E[n_b] = N eps_b^MC [1 + u1 <cos(a-b)> + u2 <cos 2(a-b)>
                             + Pbar (H A)_b],                          (8)

    N free (so only the SHAPE of eps^MC is used, not its normalization),
    eps^MC taken from the acceptance Monte Carlo at u = 0, P_zz = 0 -- i.e.
    what an analysis takes from its own acceptance simulation -- and the
    small Pbar (H A)_b term held at the harmonic coefficients `harmonics`.
    The result is only as good as that acceptance shape; what it buys is
    that the ZEUS LPS measurement becomes a PRIOR on u rather than the
    input, and that the leakage a_t du2 <cos 2beta> into a_e is bounded by
    the data's own statistics.

    N profiles out in closed form, N(u) = n_tot/S with
    S = sum_b eps_b^MC m_b, leaving a two-parameter Poisson likelihood in
    (u1, u2) solved by Newton with the same positivity guard as the
    harmonic fit.  Bins with eps_b^MC <= 0 are dropped (the acceptance MC
    predicts nothing there).

    `counts` (F, Ka, Kb) or (Ka, Kb) -- summed over fills either way;
    `eps_mc` (Ka, Kb) or flat, up to any normalization; `harmonics` the
    tensor coefficients (kappa, a_e, a_t, a_m [, a_e_s, a_t_s, a_m_s]) or
    a fit dict, None for zero.  Returns {"u1", "u2", "cov", "err_u1",
    "err_u2", "norm", "n_iter", "nll", "n_live"}."""
    n = np.asarray(counts, dtype=float)
    if n.ndim == 3:
        n = n.sum(axis=0)
    n = n.ravel()
    eps = np.asarray(eps_mc, dtype=float).ravel()
    if eps.size != n.size:
        raise ValueError("eps_mc must have one entry per (alpha, beta) bin")
    p = np.asarray(pzz, dtype=float)
    lw = np.asarray(lumis, dtype=float)
    pbar = float((lw / lw.sum() * p).sum())
    basis = basis_2d(alpha_edges, beta_edges, beta_means)
    cols = ["e", "t", "m"] + (["e_s", "t_s", "m_s"] if with_sin else [])
    if harmonics is None:
        tens = np.zeros(n.size)
    else:
        if isinstance(harmonics, dict):
            a = [harmonics.get("const", 0.0)] + [harmonics.get("a_" + c, 0.0)
                                                 for c in cols]
        else:
            a = list(harmonics)
        tens = a[0] + sum(ai * basis[c] for ai, c in zip(a[1:], cols))
    live = eps > 0.0
    h = np.vstack([basis["u1"], basis["u2"]]).T[live]
    eps, n, base = eps[live], n[live], (1.0 + pbar * tens)[live]
    n_tot = float(n.sum())
    theta = np.zeros(2)
    it = 0
    for it in range(1, int(max_iter) + 1):
        m = base + h @ theta
        s = float((eps * m).sum())
        ge = eps @ h
        grad = (n / m) @ h - n_tot * ge / s
        info = (h.T @ ((n / m ** 2)[:, None] * h)
                - n_tot * np.outer(ge, ge) / s ** 2)
        step = np.linalg.solve(info, grad)
        lam = 1.0
        for _ in range(60):
            if np.all(base + h @ (theta + lam * step) > 1e-6):
                break
            lam *= 0.5
        theta = theta + lam * step
        if np.max(np.abs(lam * step)) < tol:
            break
    m = base + h @ theta
    s = float((eps * m).sum())
    norm = n_tot / s
    ge = eps @ h
    # expected information at the fitted expectation: a weighted
    # covariance of h/m, hence positive semi-definite by construction
    info = norm * (h.T @ ((eps / m)[:, None] * h) - np.outer(ge, ge) / s)
    covm = np.linalg.inv(info)
    err = np.sqrt(np.diag(covm))
    nu = norm * eps * m
    nll = -float((n * np.log(nu, where=n > 0.0, out=np.zeros_like(nu))
                  - nu).sum())
    return {"u1": float(theta[0]), "u2": float(theta[1]), "cov": covm,
            "err_u1": float(err[0]), "err_u2": float(err[1]),
            "norm": norm, "n_iter": it, "nll": nll,
            "n_live": int(live.sum())}

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


def per_fill_acceptance(acceptance, n_fills):
    """Normalize an `acceptance` argument to a list of n_fills callables
    (or Nones).  Accepts None, one callable common to every fill, or a
    sequence of one callable (or None) per fill -- the fill-DEPENDENT
    case, which is the only phi' efficiency the spin-state ratio cannot
    cancel (docs/code_review_2026-08-25.md F1)."""
    if acceptance is None or callable(acceptance):
        return [acceptance] * n_fills
    acc = list(acceptance)
    if len(acc) != n_fills:
        raise ValueError("acceptance sequence has %d entries for %d fills"
                         % (len(acc), n_fills))
    return acc


def fill_acceptance_bias(eps2, pzz_list, lumi_fractions=None):
    """Fake cos 2phi' amplitude from a FILL-DEPENDENT efficiency harmonic.

    For yields N_fi = L_f eps_f(phi_i) sigma_i [1 + P_f T_i] with
    eps_f = eps_0(phi)(1 + e_f cos 2phi'), the spin-state ratio inverts to
    T_hat = T + (sum_f l_f (P_f - Pbar) e_f) / sigma_P^2 cos 2phi', i.e.

        dA_hat = sum_f l_f (P_f - Pbar) e_f / sigma_P^2,

    to first order in e_f (l_f = luminosity fractions).  Two equal-
    luminosity fills give the closed form (e_+ - e_0)/(P_+ - P_0) of the
    code review: with the m = +-1-rich / m = 0-rich pattern
    (P_+ - P_0 = 1.8) a 1e-3 difference of the harmonic between the two
    samples fakes 5.6e-4 -- 56% of a Delta/F1 ~ 1e-3 signal, and ~5% of
    the sweet-spot amplitudes.  A COMMON harmonic (all e_f equal) gives
    exactly zero: that is the cancellation the estimator rests on, and
    the requirement it turns into is bunch-by-bunch alternation, or
    stability of the harmonic to ~1e-4 between fills."""
    e = np.asarray(eps2, dtype=float)
    p = np.asarray(pzz_list, dtype=float)
    lf = (np.full(p.size, 1.0 / p.size) if lumi_fractions is None
          else np.asarray(lumi_fractions, dtype=float)
          / np.sum(lumi_fractions))
    pbar = float((lf * p).sum())
    sig2 = float((lf * (p - pbar) ** 2).sum())
    return float((lf * (p - pbar) * e).sum() / sig2)


def expected_counts_by_fill(n_total, pzz_list, amp, edges, acceptance=None,
                            lumi_fractions=None, const=0.0):
    """Exact expected phi'-bin counts per fill type for pseudo-experiments:
    N_fi = L_f eps_f,i (1 + P_f (const + amp cos 2phi'))  integrated over
    the bin.  `acceptance` is None, one smooth eps(phi) common to every
    fill -- the thing the ratio estimator cancels -- or a SEQUENCE of one
    per fill, the fill-dependent case it cannot (fill_acceptance_bias)."""
    p = np.asarray(pzz_list, dtype=float)
    f = (np.full(p.size, 1.0 / p.size) if lumi_fractions is None
         else np.asarray(lumi_fractions, dtype=float) / np.sum(lumi_fractions))
    edges = np.asarray(edges, dtype=float)
    nsub = 32
    sub = edges[:-1, None] + (np.arange(nsub)[None, :] + 0.5) / nsub \
        * (edges[1:] - edges[:-1])[:, None]
    accs = per_fill_acceptance(acceptance, p.size)
    dphi = (edges[1:] - edges[:-1]) / (2.0 * np.pi)
    out = np.empty((p.size, edges.size - 1))
    for k, (pf, ff, acc) in enumerate(zip(p, f, accs)):
        eps = np.ones_like(sub) if acc is None else acc(sub)
        base = eps.mean(axis=1)
        mod = (eps * np.cos(2.0 * sub)).mean(axis=1)
        out[k] = n_total * ff * dphi * ((1.0 + pf * const) * base
                                        + pf * amp * mod)
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


# Geometric aperture of the ePIC Roman Pots for a track at beam rigidity,
# measured by shooting an intact 6Li through the epic-main geometry and
# reading where the hits start (tools/fullsim/README.md, 2026-08-26).
# Half-widths in ANGLE at the IP, per optics configuration; the pots are
# positioned per ring, so the aperture is not one number.
#
# It is not the beam envelope and does not scale with it: the beamline
# images an IP angle onto the pot plane with a horizontal lever of
# R12 = 30.6 m against a vertical few metres, so what clears the pots'
# horizontal slot is theta_x, and the boundary is a property of the
# optics and the mechanics.  The envelope and the aperture are separate
# constraints and a track must clear BOTH; `rp_measure` takes the larger
# per axis.
#
# Caveats, unchanged from the measurement: one event per scan point, 30
# degree azimuthal steps, no beam envelope in the scan, and a
# September-2024 epic-main whose pot stations have already moved once.
RP_APERTURE_MEASURED = {
    "5x41": (2.0e-3, 3.0e-3),
    "10x100": (1.35e-3, 3.0e-3),
    "18x275": (1.03e-3, 2.3e-3),
}


def rp_aperture_for(config, table=None):
    """The measured aperture of a machine configuration.

    `config` is a `beams.BeamConfig` OF ANY SPECIES, resolved through
    `farforward.yr_config_key` -- the aperture is a property of the ring
    and the pot mechanics, not of the beam in it, so every isotope at a
    given machine configuration sees the same slot.

    A bare per-nucleon momentum [GeV/u] is still accepted, for the callers
    written before 2026-08-28, and is matched against the 6Li
    configurations alone; it returns None off those three points rather
    than interpolating, since the pot positions are set per ring and not by
    a formula.  That path CANNOT serve 7Li: 7Li's top configuration is
    117.9 GeV/u against 6Li's 137.5, so `rp_aperture_for(117.9)` returned
    None and the near-beam scans could not be run for it at all (plans/09
    B3).  Pass the configuration.

    The three 6Li momenta are derived from `beams`, not hard-coded -- they
    moved on 2026-08-27 when the two lower configurations were corrected
    from rigidity-scaled (20.5, 50) to gamma-matched (40.8, 99.5) GeV/u
    (plans/10)."""
    from polli_fastsim import beams as _beams
    from polli_fastsim import farforward as _ff
    table = table or RP_APERTURE_MEASURED
    if hasattr(config, "ion_momentum_per_nucleon"):
        return table.get(_ff.yr_config_key(config))
    keys = ("5x41", "10x100", "18x275")
    for cfg, key in zip(_beams.default_configs("6Li"), keys):
        if abs(float(config) - cfg.ion_momentum_per_nucleon) < 1e-3:
            return table[key]
    return None


def rp_measure(Pprime, P, sigma_theta_xy, n_sigma=10.0, rng=None,
               shape="rectangle", sigma_pos_over_l=0.0,
               cut_scale_xy=(1.0, 1.0), cut_theta_xy=None):
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

    `cut_theta_xy` adds the pots' own GEOMETRIC aperture as absolute
    half-widths in angle (rad), and the cutout becomes the larger of the
    two per axis: a track has to clear the beam envelope AND the
    mechanics.  `RP_APERTURE_MEASURED` is that aperture, measured in the
    ePIC geometry, and it dominates the envelope at every configuration
    -- and inverts its aspect, because it is theta_x that clears the
    slot.  None keeps the envelope alone, which is what every number
    published before 2026-08-26 used.

    Returns dict with theta_x/y, pT, phi_t, t_reco = -pT^2 (x_L set to 1),
    the acceptance mask and the cut half-widths in GeV.

    Two idealisations are stated rather than modelled.  The pots measure a
    POSITION at the pot plane, x = R_12 theta_x + D delta, and the recoil
    is off-momentum by delta = -x_P (x_L = 1 - x_P): at the ePIC dispersion
    (D ~ 0.3 m against R_12 = 30.6 m, tools/fullsim) an x_P = 0.01 recoil
    is displaced by the same amount as a 0.1 mrad angle, i.e. below the
    divergence of every configuration; that degeneracy, and the beam's own
    momentum spread that enters the same way, are what `sigma_pos_over_l`
    and the dispersive term of `tagging_optics_point` stand in for.  And
    the pots see the beam-referenced angle, so the crossing angle and the
    reference orbit drop out by construction.
    """
    rng = rng or np.random.default_rng(20260824)
    Pp = np.asarray(Pprime, dtype=float)
    P = np.asarray(P, dtype=float)
    p_beam = np.sqrt((P[..., 1:] ** 2).sum(axis=-1))
    thx = Pp[..., 1] / Pp[..., 3]
    thy = Pp[..., 2] / Pp[..., 3]
    sx, sy = sigma_theta_xy
    cx, cy = n_sigma * sx * cut_scale_xy[0], n_sigma * sy * cut_scale_xy[1]
    if cut_theta_xy is not None:
        cx = max(cx, float(cut_theta_xy[0]))
        cy = max(cy, float(cut_theta_xy[1]))
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
    outside a near-beam cutout -- polli_fastsim.farforward.hole_acceptance,
    re-exported (2026-08-28) so that the fast simulation and the
    reconstruction chain share one implementation."""
    from polli_fastsim import farforward as _ff
    return _ff.hole_acceptance(slope_b, cut_x, cut_y, shape=shape, nphi=nphi)
