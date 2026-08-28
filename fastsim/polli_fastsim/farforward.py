"""Far-forward detector windows and fragment routing (Phase-1 stand-in).

Windows fetch-verified 2026-06-12 (YR detector matrix; arXiv:2108.08314
Table I; arXiv:2409.02811; details in plans/03 step 2.2):

  Roman Pots   z = 32.55/34.25 m  theta < 5 mrad  R in [0.60, 0.95]
               near-beam (|R-1| < 0.05): only theta > n_sigma sigma_theta
  OMD          z = 25.50/27.00 m  theta < 5 mrad  R in [0.45, 0.65], no cut
  B0           5.5 < theta < 20 mrad, any charged
  ZDC          theta < 4 mrad, neutrals
  no coverage  R > ~1.05 (bends less than beam), or 5-5.5 mrad gap

The z positions above are read from the CURRENT eic/epic main branch
(2026-08-27): roman_pots_eRD24_design.xml gives 32547.3 / 34245.5 mm and
offM_tracker.xml positions its four layers at 25500 / 25520 / 27000 /
27020 mm, i.e. two stations of two layers.  The 26/28 and 22.5/24.5 m this
docstring carried until then are the 2021 Yellow Report values, which the
YR, the ePIC wiki and several papers still quote -- the geometry moved and
the documents did not.  Note also that offM_tracker.xml DEFINES a constant
ForwardOffMTracker_zpos = B1APF_CenterPosition + B1APF_Length/2 + 10 cm =
22.16 m which its own <position> elements do not use; do not read it as the
station position.  Only the ANGULAR windows below enter the routing, so the
positions are documentation -- but they were wrong.

The near-beam cut is ANGULAR, not a momentum cut.  A fragment at beam
rigidity follows the beam's own optics, so what separates it from the
beam at the pots is its angle in units of the beam's angular divergence:
the documented "10 sigma" cut is theta > 10 sigma_theta.  The quoted
0.20 GeV (high-acceptance) and 0.41 GeV (high-divergence, derived from
the YR divergence tables and elsewhere rounded to 0.45) are that cut
expressed for a 275 GeV PROTON, i.e. sigma_theta = 73 and 149 microrad;
applying those numbers as a pT threshold to a nucleus of momentum A p_u
understates the envelope by A p_u / 275 GeV (code review S8).  For the
6Li alpha fragment at 137.5 GeV/u, 0.20 GeV is 5 sigma, not 10.
`Optics.pt_cut_near_beam` keeps the proton-referenced number and
`pt_cut_for(p)` gives the threshold at any momentum.

ASSUMPTION (flagged in plans/04 #11): off-rigidity tracks inside the RP
window are dispersion-separated from the beam, so no near-beam cut
applies to them; and the lithium divergence is taken equal to the
proton's -- only Phase-2 optics can refine either.
"""

from dataclasses import dataclass

import numpy as np

THETA_RP_MAX = 5.0e-3
THETA_B0_MIN, THETA_B0_MAX = 5.5e-3, 20.0e-3
THETA_ZDC_MAX = 4.0e-3
RP_R_WINDOW = (0.60, 0.95)
# Jentsch, Tu, Weiss, PRC 104 (2021) 065205 (arXiv:2108.08314) Table I:
# the off-momentum detectors cover zeta = p/p_beam in 0.45-0.65 (the
# paper's symbol is zeta, Eq. 58).  The 0.40-0.60 used before 2026-08-25
# came from an earlier reading (code review S7).
OMD_R_WINDOW = (0.45, 0.65)
NEAR_BEAM_BAND = 0.05  # |R-1| below this: inside the beam envelope

PROTON_REFERENCE_MOMENTUM = 275.0  # GeV, the momentum the cuts are quoted at


@dataclass(frozen=True)
class Optics:
    """Near-beam envelope of an optics setting, as the ANGLE it subtends.

    `sigma_theta` is the HORIZONTAL beam angular divergence at the IP and
    `sigma_theta_v` the vertical one (None = isotropic, the legacy form).
    The Roman Pots are planar, so the envelope a fragment at beam rigidity
    has to clear is a RECTANGLE of half-widths (n_sigma sigma_h,
    n_sigma sigma_v) in angle: it is tagged when |theta_x| > n sigma_h OR
    |theta_y| > n sigma_v.  Since 2026-08-28 the per-configuration Yellow
    Report values are available through `yr_optics(config)` and the
    lithium tagging optics through `tagging_optics(config)`; the two
    module-level constants below are the proton-derived legacy.
    """
    name: str
    sigma_theta: float          # rad, horizontal
    n_sigma: float = 10.0
    sigma_theta_v: float = None  # rad, vertical; None = isotropic
    lumi_fraction: float = 1.0   # luminosity relative to high acceptance

    @property
    def sigma_v(self):
        return (self.sigma_theta if self.sigma_theta_v is None
                else self.sigma_theta_v)

    @property
    def isotropic(self):
        return (self.sigma_theta_v is None
                or abs(self.sigma_theta_v - self.sigma_theta) < 1e-15)

    @property
    def envelope(self):
        """(n sigma_h, n sigma_v) [rad]: the cutout half-widths."""
        return (self.n_sigma * self.sigma_theta, self.n_sigma * self.sigma_v)

    def clears(self, theta, phi=None):
        """Boolean: does a fragment at polar angle `theta` (from the ion
        axis) and azimuth `phi` clear the envelope?  With `phi` None the
        cut is circular in theta at n sigma_h -- exact for an isotropic
        optics, and the only option when the azimuth is not known."""
        theta = np.asarray(theta, dtype=float)
        cx, cy = self.envelope
        if phi is None or self.isotropic:
            return theta > cx
        phi = np.asarray(phi, dtype=float)
        return (np.abs(theta * np.cos(phi)) > cx) | (np.abs(theta * np.sin(phi)) > cy)

    @property
    def pt_cut_near_beam(self):
        """The same cut expressed for a 275 GeV proton -- the number the
        Yellow Report quotes (0.20 / 0.41 GeV).  Use `pt_cut_for` for a
        fragment of any other momentum."""
        return self.pt_cut_for(PROTON_REFERENCE_MOMENTUM)

    def pt_cut_for(self, momentum):
        """Transverse-momentum threshold (horizontal) for a fragment of
        total momentum `momentum` [GeV] -- e.g. A * p_per_nucleon."""
        return self.n_sigma * self.sigma_theta * np.asarray(momentum,
                                                            dtype=float)


HIGH_ACCEPTANCE = Optics("high-acceptance", 0.20 / (10.0 * 275.0))    #  73 urad
# The high-divergence envelope is not a documented spec: plans/06 SS6.5
# derives ~0.41 GeV at 275 GeV from the YR divergence tables and the
# fast-sim has always quoted the rounded-up 0.45 (every published
# fast-sim number, plans/06 and the money-plot report).  polligen.reco
# quotes the 0.41 end of the same band.  Unifying the two is a
# documentation decision, not a code one (plans/04 #11, plans/08).
HIGH_DIVERGENCE = Optics("high-divergence", 0.45 / (10.0 * 275.0))    # 164 urad

# --- what the published tables actually say (plans/10) --------------------
#
# The two Optics above are ENERGY-INDEPENDENT and ISOTROPIC, and back-derived
# from a 275 GeV proton.  The Yellow Report's own beam tables are neither.
# Table 10.1 (e+p) and 10.2 (e+Au) give RMS divergence h/v and dp/p per
# configuration; these are those tables, verbatim, for the HADRON beam.
#
#   config: (HD_h, HD_v, HA_h, HA_v) [urad], dp_over_p [1e-4]
#
# Note that at 41 GeV the Yellow Report lists no separate high-acceptance
# option -- the two rows are identical -- and that the divergence is
# isotropic at 275 and 100 GeV but 220/380 at 41.
YR_PROTON_DIVERGENCE = {
    "18x275": ((150.0, 150.0, 65.0, 65.0), 6.8),    # e 18 x p 275, 290 bunches
    "10x100": ((220.0, 220.0, 180.0, 180.0), 9.7),  # e 10 x p 100 (e 5: 206/206)
    "5x41":   ((220.0, 380.0, 220.0, 380.0), 10.3),
}
YR_GOLD_DIVERGENCE = {                              # strong hadron cooling
    "110GeV/u": ((218.0, 379.0), 6.2),
    "41GeV/u":  ((275.0, 377.0), 10.0),
}

#: sigma_theta = sqrt(eps_N / (beta*gamma * beta*)).  At a given machine
#: configuration the lattice is set by RIGIDITY, so beta* is common to every
#: species and an A/Z = 2 ion sits at HALF the proton's beta*gamma -- hence
#: sqrt(2) more divergence at equal normalised emittance.  Calibrating that
#: assumption against the published gold rows gives eps_N(Au)/eps_N(p) = 0.85
#: horizontally and 2.6 vertically, and gold's IBS (~ N Z^4/A^2) is 450x
#: lithium's -- so equal eps_N is well supported for a light ion and the
#: lithium divergence is set by KINEMATICS, not by intrabeam scattering.
LIGHT_ION_DIVERGENCE_FACTOR = 1.409     # sqrt(beta*gamma_p / beta*gamma_Li)


def yr_divergence_for(config, a_over_z=2.0, optics="high-acceptance"):
    """Estimated (sigma_theta_h, sigma_theta_v) [rad] and dp/p for a fully
    stripped ion of the given A/Z in an EIC configuration, from the Yellow
    Report proton tables scaled by sqrt(A/Z) (plans/10).

    This is an ESTIMATE, not a machine specification: it assumes the ion's
    normalised emittance equals the proton's, which the gold calibration
    supports for a light ion but which no published light-ion optics
    confirms.  plans/10 D1 is the question that would replace it."""
    (hd_h, hd_v, ha_h, ha_v), dpp = YR_PROTON_DIVERGENCE[config]
    h, v = (ha_h, ha_v) if optics == "high-acceptance" else (hd_h, hd_v)
    f = float(np.sqrt(a_over_z))
    return (1e-6 * h * f, 1e-6 * v * f), 1e-4 * dpp


def route_charged(R, theta, pT, optics=HIGH_ACCEPTANCE, phi=None):
    """Classify charged fragments into far-forward systems.

    Returns an integer array: 0 lost, 1 RP, 2 OMD, 3 B0, 4 RP-near-beam
    (R ~ 1, accepted only outside the angular envelope -- a rectangle of
    half-widths n_sigma (sigma_h, sigma_v) when the fragment azimuth `phi`
    is given, a circle at n_sigma sigma_h otherwise; Optics.clears).  `pT`
    is not used for the near-beam decision: the envelope is angular.
    """
    R = np.asarray(R, dtype=float)
    theta = np.asarray(theta, dtype=float)
    pT = np.asarray(pT, dtype=float)
    out = np.zeros(R.shape, dtype=int)

    in_b0 = (theta >= THETA_B0_MIN) & (theta <= THETA_B0_MAX)
    out[in_b0] = 3

    small = theta < THETA_RP_MAX
    near = np.abs(R - 1.0) < NEAR_BEAM_BAND
    in_rp = small & ~near & (R >= RP_R_WINDOW[0]) & (R <= RP_R_WINDOW[1])
    in_omd = small & (R >= OMD_R_WINDOW[0]) & (R < OMD_R_WINDOW[1])
    rp_tail = small & near & optics.clears(theta, phi)
    out[in_omd] = 2
    out[in_rp] = 1
    out[rp_tail] = 4
    return out


ROUTE_LABELS = {0: "lost", 1: "RomanPots", 2: "OMD", 3: "B0",
                4: "RP (pT tail, R~1)"}


def route_neutral(theta):
    """Neutral fragments: ZDC inside THETA_ZDC_MAX, else lost (B0 EMCal
    photons not modeled here). Returns 5 for ZDC, 0 for lost."""
    theta = np.asarray(theta, dtype=float)
    return np.where(theta <= THETA_ZDC_MAX, 5, 0)


def neutral_summary(theta):
    route = route_neutral(theta)
    n = float(len(route))
    return {"ZDC": float(np.sum(route == 5)) / n,
            "lost": float(np.sum(route == 0)) / n}


def acceptance_summary(R, theta, pT, optics=HIGH_ACCEPTANCE, phi=None):
    """Fraction of spectators in each far-forward system."""
    route = route_charged(R, theta, pT, optics, phi=phi)
    n = float(len(route))
    return {label: float(np.sum(route == code)) / n
            for code, label in ROUTE_LABELS.items()}


# --- per-configuration optics (2026-08-28) ---------------------------------

#: Far-forward transport constants read off the ePIC geometry scan
#: (tools/fullsim, 18 x 275): the horizontal lever from an IP angle to the
#: pot-plane position and the dispersion at the pots.
POT_R12 = 30.6        # m
POT_DISPERSION = 0.30  # m


def yr_config_key(config):
    """Which Yellow Report configuration ("5x41", "10x100", "18x275") a
    BeamConfig belongs to: the proton energy whose gamma-matched (or
    rigidity-capped) per-nucleon momentum reproduces the ion's."""
    from . import beams as _beams
    key = {41.0: "5x41", 100.0: "10x100", 275.0: "18x275"}
    p_e = min(_beams.PROTON_CONFIG_ENERGIES,
              key=lambda e: abs(config.ion.momentum_per_nucleon_at(e)
                                - config.ion_momentum_per_nucleon))
    return key[p_e]


def sigma_theta_for(config, optics="high-acceptance"):
    """(sigma_theta_h, sigma_theta_v) [rad] for a beam configuration.

    Two steps, both from plans/10:

    1.  The PROTON divergence of that machine configuration, from Yellow
        Report Table 10.1 (YR_PROTON_DIVERGENCE).
    2.  The species step.  An ion is GAMMA-MATCHED to its proton
        configuration unless the ring rigidity caps it first
        (beams.Ion.momentum_per_nucleon_at).  A gamma-matched ion has the
        SAME beta*gamma as the proton, hence the same geometric emittance at
        equal eps_N, hence the SAME divergence -- no penalty.  A
        rigidity-capped one sits at lower beta*gamma and picks up
        sqrt(beta*gamma_p / beta*gamma_ion), which for 6Li at the top
        configuration is sqrt(2).

    So the correction is not a blanket factor: it is 1.00 at the two lower
    configurations and 1.41 at the top, and it rides on top of a proton
    divergence that itself varies 65 -> 180 -> 220 microrad.  (Moved here
    from polligen.reco on 2026-08-28 so that the fast simulation's own
    spectator routing can use it; reco.sigma_theta_for delegates.)
    """
    from . import beams as _beams
    key = yr_config_key(config)
    p_e = {"5x41": 41.0, "10x100": 100.0, "18x275": 275.0}[key]
    (hd_h, hd_v, ha_h, ha_v), _ = YR_PROTON_DIVERGENCE[key]
    h, v = (ha_h, ha_v) if optics == "high-acceptance" else (hd_h, hd_v)
    g_p = ((p_e ** 2 + _beams.PROTON_MASS ** 2) ** 0.5) / _beams.PROTON_MASS
    bg_p = (g_p ** 2 - 1.0) ** 0.5
    bg_i = config.ion_momentum_per_nucleon / config.ion.mass_per_nucleon
    f = (bg_p / bg_i) ** 0.5
    return 1e-6 * h * f, 1e-6 * v * f


def yr_optics(config, optics="high-acceptance", n_sigma=10.0):
    """The Yellow Report optics of a beam configuration as an `Optics`:
    the per-configuration, anisotropic divergence of `sigma_theta_for`
    with the 10 sigma envelope, for the species of `config`."""
    h, v = sigma_theta_for(config, optics)
    return Optics("%s %s" % (yr_config_key(config), optics), h, n_sigma, v)


def hole_acceptance(slope_b, cut_x, cut_y, shape="rectangle", nphi=3600):
    """Azimuthal acceptance of an exponential coherent recoil spectrum
    dN/d2pT ~ exp(-B pT^2) outside a near-beam cutout with half-widths
    (cut_x, cut_y) in pT [GeV]: eps(phi) = exp(-B rho(phi)^2) with rho the
    cutout boundary along phi.  Returns {"phi", "eps", "acc", "a2", "a4"}:
    the total acceptance and the cos 2phi / cos 4phi Fourier coefficients
    <cos n phi> of the TAGGED sample about the horizontal axis -- the
    fake modulation a single-fill fit would attribute to physics (a2 = 0
    only for cut_x = cut_y; a4 != 0 even then).  polligen.reco re-exports
    this as rp_hole_acceptance."""
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


def tagging_optics_point(config, slope_b=50.0, n_sigma=10.0, r_max=2000.0,
                         n_grid=400, dispersion=True, optics="high-acceptance"):
    """The lithium TAGGING OPTICS of Report 1 Section 6.1: the working
    point that maximises (tagged fraction) x (luminosity) when the
    HORIZONTAL beta* alone is raised by r over the high-acceptance value,
    the vertical plane is held, and the Roman Pots follow the n_sigma
    envelope in both planes.

    Returns a dict with r_h (beta*_x / beta*_x,HA), sigma_x (the horizontal
    RMS angle at the IP the de-squeeze leaves), sigma_x_eff (the same with
    the beam's momentum spread through the pot dispersion added in
    quadrature -- the angular smearing a pot-plane measurement sees, and
    the quantity the envelope is n_sigma of), sigma_y, env_x / env_y (the
    envelope half-widths in angle), acceptance, lumi_fraction
    (= 1/sqrt(r_h)), p_ion, and `optics`: the working point as an
    `Optics` (sigma_x_eff, sigma_y, the luminosity fraction), so that the
    spectator routing can be evaluated at it.  Identical to the scan of
    evgen/scripts/tagging_optics.py (same grid, same acceptance), pinned
    by a test.
    """
    sh, sv = sigma_theta_for(config, optics)
    p_ion = config.ion.A * config.ion_momentum_per_nucleon
    dpp = 1e-4 * YR_PROTON_DIVERGENCE[yr_config_key(config)][1]
    disp = (POT_DISPERSION * dpp / POT_R12) if dispersion else 0.0
    r = np.logspace(np.log10(0.25), np.log10(r_max), n_grid)
    best = None
    for rr in r:
        sx = sh / rr ** 0.5
        sx_eff = np.hypot(sx, disp)
        acc = hole_acceptance(slope_b, n_sigma * sx_eff * p_ion,
                              n_sigma * sv * p_ion)["acc"]
        prod = acc / rr ** 0.5
        if best is None or prod > best["product"]:
            best = {"r_h": float(rr), "sigma_x": float(sx),
                    "sigma_x_eff": float(sx_eff), "sigma_y": float(sv),
                    "env_x": float(n_sigma * sx_eff),
                    "env_y": float(n_sigma * sv), "acceptance": float(acc),
                    "lumi_fraction": float(1.0 / rr ** 0.5),
                    "product": float(prod), "p_ion": float(p_ion),
                    "n_sigma": float(n_sigma)}
    best["optics"] = Optics("%s tagging (beta*_x x %.0f)"
                            % (yr_config_key(config), best["r_h"]),
                            best["sigma_x_eff"], n_sigma, best["sigma_y"],
                            best["lumi_fraction"])
    return best


def tagging_optics(config, **kw):
    """The tagging optics of `tagging_optics_point` as an `Optics`."""
    return tagging_optics_point(config, **kw)["optics"]
