"""Far-forward detector windows and fragment routing (Phase-1 stand-in).

Windows fetch-verified 2026-06-12 (YR detector matrix; arXiv:2108.08314
Table I; arXiv:2409.02811; details in plans/03 step 2.2):

  Roman Pots   z = 32.55/34.25 m  theta < 5 mrad  R in [0.60, 0.95]
               near-beam (|R-1| < 0.05): only theta > n_sigma sigma_theta
  OMD          z = 25.50/27.00 m  theta < 5 mrad  R in [0.45, 0.65], no cut
  B0           5.5 < theta < 20 mrad, any charged
  ZDC          theta < 4 mrad, neutrals
  RP-inner     R > ~1.05 (bends less than beam): the OVER-RIGID branch,
               route code 6, whenever the dispersive displacement clears
               the per-configuration blind block (48 / 32 / 16 mm) and
               stays inside the last module at 144 mm -- see
               `over_rigid_route` below, measured 2026-08-28
  no coverage  the 5-5.5 mrad gap, and a dispersive displacement past the
               last module (e.g. R = 1.504)

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
0.20 GeV (high-acceptance) and 0.45 GeV (high-divergence) are that cut
expressed for a 275 GeV PROTON, i.e. sigma_theta = 73 and 164 microrad;
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
#: Outer angular edge of the Roman-Pot sensor package, as a polar angle at
#: the IP.  THETA_RP_MAX is the acceptance bound the Yellow Report detector
#: matrix quotes and it is what the routing has always used, so the default
#: keeps every published number where it is; the constant exists because a
#: real pot is a tiling of 16 x 16 mm modules a few tens of millimetres
#: across, not a 150 mm disc (5 mrad x R12 = 29.97 m at 18 x 275; 96 and
#: 106 mm at 5 x 41 and 10 x 100), and any topology
#: question involving a SECOND fragment is dominated by where that edge is
#: (plans/09 B1, D3).  Pass `theta_outer` to `route_charged` to price it.
THETA_RP_OUTER = THETA_RP_MAX
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
        axis) and azimuth `phi` clear the envelope?

        The envelope is the RECTANGLE of half-widths n_sigma (sigma_h,
        sigma_v) -- the same cutout `hole_acceptance` integrates and
        `recopseudo.CoherentResponse` applies -- so with `phi` given the
        cut is rectangular at every aspect ratio, an isotropic optics
        included: there the rectangle is a SQUARE of half-width n sigma_h,
        which is not its inscribed circle.  With `phi` None the cut falls
        back to that inscribed circle, the only option when the azimuth is
        not known; it is the more generous of the two (by 1.7x at the
        tagging optics, plans/09 B2) and no published number uses it.
        Until 2026-08-28 an `or self.isotropic` short-circuit sent every
        square envelope down the circular branch, so the azimuth
        `polligen.tagged` and the spectator routing take care to pass was
        discarded at 10 x 100 and 18 x 275, where the Yellow Report
        divergence is 180/180 and 92/92 microrad.
        """
        theta = np.asarray(theta, dtype=float)
        cx, cy = self.envelope
        if phi is None:
            return theta > cx
        phi = np.asarray(phi, dtype=float)
        return (np.abs(theta * np.cos(phi)) > cx) | (np.abs(theta * np.sin(phi)) > cy)

    @property
    def pt_cut_near_beam(self):
        """The same cut expressed for a 275 GeV proton -- the number the
        Yellow Report quotes (0.20 / 0.45 GeV).  Use `pt_cut_for` for a
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
# fast-sim rounds it up to 0.45 (every published fast-sim number,
# plans/06 and the money-plot report).  The two ends of that band were
# carried by two modules until 2026-08-28 (R10), when polligen.reco's
# SIGMA_THETA_HD became an alias of this constant: 164 urad is now the
# one high-divergence number in the repository (plans/04 #11, plans/08).
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

# --- the species step, and what calibrates it (plans/10) ------------------
#
# sigma_theta = sqrt(eps_N / (beta*gamma * beta*)).  At a given machine
# configuration the lattice is set by RIGIDITY, so beta* is common to every
# species and an ion at HALF the proton's beta*gamma picks up sqrt(2) more
# divergence at equal normalised emittance -- which is a statement about
# beta*gamma, not about A/Z, and applies only where the ring rigidity caps
# the ion below its gamma-matched momentum (sigma_theta_for; plans/10
# SS10.3).  A blanket sqrt(A/Z) at every configuration is the
# rigidity-scaling rule that plans/10 rejects, and the function that
# applied it was retired on 2026-08-28.
#
# Calibrating the equal-emittance assumption against the published gold
# rows gives eps_N(Au)/eps_N(p) = 0.85 horizontally and 2.6 vertically,
# and gold is far more IBS-prone than lithium under either normalisation
# of the same law: the per-particle growth rate goes as Z^4/A^2, which
# makes gold 446x lithium, and at fixed BEAM CURRENT (the comparison
# plans/10 SS10.3 states, one factor Z removed because a fixed current is
# fewer ions when each carries more charge) as Z^3/A^2, which makes it
# 17x.  Either way lithium is proton-class -- 2.25x a proton per particle,
# 0.75x at fixed current -- so equal eps_N is well supported for it and
# the lithium divergence is set by KINEMATICS, not by intrabeam
# scattering.


def route_charged(R, theta, pT, optics=HIGH_ACCEPTANCE, phi=None,
                  theta_outer=None, pot_config="18x275"):
    """Classify charged fragments into far-forward systems.

    Returns an integer array: 0 lost, 1 RP, 2 OMD, 3 B0, 4 RP-near-beam
    (R ~ 1, accepted only outside the angular envelope -- a rectangle of
    half-widths n_sigma (sigma_h, sigma_v) when the fragment azimuth `phi`
    is given, a circle at n_sigma sigma_h otherwise; Optics.clears), 6
    RP-inner (over-rigid).  `pT` is not used for the near-beam decision:
    the envelope is angular.

    `theta_outer` replaces THETA_RP_OUTER as the outer edge of the pot and
    off-momentum acceptance; the default leaves it at THETA_RP_MAX, where
    every published number was computed (`theta_rp_outer_for` gives the
    measured edge, 2.9-4.0 mrad, for callers that want to price it).  The
    B0 window is untouched by it -- a fragment between the pot edge and
    5.5 mrad is lost, as it is now.

    THE OVER-RIGID BRANCH (2026-08-28, plans/09 B1).  Until then every
    R > 1 + NEAR_BEAM_BAND fragment fell through to "lost" by
    construction, and Report 3 Table 6 said so in as many words: "the
    triton (R = 1.29) is over-rigid, with no window in route_charged,
    which carries no R > 1 branch ... its 'no coverage' is a routing
    assumption rather than a measurement".  It is now measured, and it was
    wrong: the pot dispersion carries an over-rigid fragment to the INNER
    side of the bend and an R = 1.286 triton is on the silicon in 60 of 60
    events at every configuration (`over_rigid_route`, which holds the
    numbers).  `pot_config` selects the configuration whose blind block
    the displacement has to clear; it defaults to 18 x 275, the most
    permissive of the three, and the fragments this changes are the 7Li
    triton at R = 1.290 (now RP-inner at all three) and NOT the 6Li 3He+t
    triton at R = 1.504, whose 151 mm displacement is past the last
    module and which stays lost.  `phi` enters the branch when it is
    given -- an over-rigid fragment can be pushed back into the central
    block by a horizontal angle of the opposite sign, which is a real
    measured hole: theta_x in -1.7 to -2.7 mrad in this model at
    18 x 275, and -1.55 to -2.53 mrad in the scan it was fitted to
    (`lad_triton_18x275`; RESULTS section 5).
    """
    R = np.asarray(R, dtype=float)
    theta = np.asarray(theta, dtype=float)
    pT = np.asarray(pT, dtype=float)
    out = np.zeros(R.shape, dtype=int)

    in_b0 = (theta >= THETA_B0_MIN) & (theta <= THETA_B0_MAX)
    out[in_b0] = 3

    small = theta < (THETA_RP_OUTER if theta_outer is None else theta_outer)
    near = np.abs(R - 1.0) < NEAR_BEAM_BAND
    in_rp = small & ~near & (R >= RP_R_WINDOW[0]) & (R <= RP_R_WINDOW[1])
    in_omd = small & (R >= OMD_R_WINDOW[0]) & (R < OMD_R_WINDOW[1])
    rp_tail = small & near & optics.clears(theta, phi)
    theta_x = theta * np.cos(phi) if phi is not None else 0.0
    over = (small & (R > 1.0 + NEAR_BEAM_BAND)
            & over_rigid_route(R, theta_x, pot_config))
    out[in_omd] = 2
    out[in_rp] = 1
    out[rp_tail] = 4
    out[over] = 6
    return out


ROUTE_LABELS = {0: "lost", 1: "RomanPots", 2: "OMD", 3: "B0",
                4: "RP (pT tail, R~1)", 6: "RP-inner (over-rigid)"}


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


def acceptance_summary(R, theta, pT, optics=HIGH_ACCEPTANCE, phi=None,
                       theta_outer=None, pot_config="18x275"):
    """Fraction of spectators in each far-forward system.

    `pot_config` is passed straight to `route_charged` and MUST be set by
    any caller that works one machine configuration at a time: the
    over-rigid branch tests the pot-plane displacement against that
    configuration's blind block (48 / 32 / 16 mm), and the default
    18 x 275 is the most permissive of the three.  Leaving it at the
    default while sweeping configurations inflates the R > 1.05 share at
    5 x 41 and 10 x 100 by a factor 8 (0.0097 -> 0.0011 and 0.0015 of the
    6Li alpha-tag sample), which is the single-number-everywhere error
    this measurement exists to remove.
    """
    route = route_charged(R, theta, pT, optics, phi=phi,
                          theta_outer=theta_outer, pot_config=pot_config)
    n = float(len(route))
    return {label: float(np.sum(route == code)) / n
            for code, label in ROUTE_LABELS.items()}


# --- per-configuration optics (2026-08-28) ---------------------------------

# --- the pot-plane transport, measured per configuration (2026-08-28) -----
#
# `tools/fullsim`, plans/09 B1, in eic_xl-nightly / epic-main
# 9aaa296976d3ad9de404f775ae89fc17a068c07c.  An intact 6Li was walked
# along theta = 0.20-6.00 mrad in 0.05 steps at phi = 0/90/180/270 at each
# ring's reference rigidity, and the first-crossing hit position regressed
# against the IP angle plane by plane:
#
#   (R12, R34, D)  =  d x / d theta_x ,  d y / d theta_y ,  d x / d R
#
# R12 residual rms 0.67-0.82 mm over 42-141 rows per plane, and the fitted
# intercept lands on the pot centre to 0.6 mm at every configuration and
# in every plane -- the check that the regression measures transport and
# not a mis-set reference orbit.  Hits are assigned to planes in the
# ROTATED station frame (the stations are tilted -0.04545 rad about y);
# assigning on global z instead, as the reader did until the same day,
# double-counts every hit beyond |dx| ~ 110 mm into both layers of its
# station and moves these levers by 0.04-0.1%.
#
# Three things this replaces.  R12 = 30.6 m was ONE number, measured at
# 18 x 275 in the September-2024 geometry and applied everywhere: the
# lower configurations are 19.2 and 21.3 m, a third of it, so every
# millimetre quoted at 5 x 41 and 10 x 100 before today was 45-60% high.
# R34 had never been measured at all and defaulted to R12: it is 2.6-3.4 m,
# an order of magnitude smaller, and the far-forward line is ~9x stiffer
# in x than in y -- which is why the pot aperture is a horizontal slot
# even though the pots insert vertically.  And D = 0.30 m, which WAS a
# single 18 x 275 number, survives: 0.31 / 0.29 / 0.29 m measured.
#
# R34 is None at 5 x 41 and that is a measurement, not a gap.  The
# per-energy insertion there holds the central silicon off to |y| >= 29.6
# mm and the intermediate bands to 27.5 and 18.0, so a vertical kick needs
# theta_y ~ 10 mrad to reach any sensor: the whole vertical ladder
# returned one accepted row in 0.2-6.0 mrad.  The vertical plane at 5 x 41
# is shut, and `separation_at_pots` falls back on R12 there with a warning
# in its docstring rather than inventing a lever.
POT_LEVERS = {
    "5x41":   (19.24, None, 0.311),
    "10x100": (21.25, 3.35, 0.287),
    "18x275": (29.97, 2.93, 0.292),
}

#: THE LIGHT-ION-LATTICE ALTERNATIVE at 5 x 41, and the size of the
#: systematic the choice of compact file carries.
#:
#: `epic-main` ships two different lattices for the 5 x 41 ring setting.
#: `compact/fields/beamline_5x41.xml` is the 41 GeV PROTON one and is what
#: `epic_craterlake_5x41.xml` includes, so it is what every scan in this
#: repository and every number in POT_LEVERS stands on.
#: `compact/fields/beamline_5x41_He4.xml` is the Z/A = 0.5 one -- 82 GV,
#: FieldScaleFactor 82/275 on the Yellow Report 275 GeV magnet set -- and
#: is the lattice a 6Li fill at the gamma-matched 40.8 GeV/u (81.6 GV)
#: might actually run in.  They are NOT a field scale of one another
#: (Q1APF_GradientMax -15.38 against -72.61 T/m), and the transport
#: differs by a factor 1.55: R12 29.81 m against 19.24, and the
#: horizontal silicon edge 1.60-1.70 mrad against 2.50 (RESULTS section
#: 6, `lad_he4_5x41`, 6Li at 40.8 GeV/u, theta 0.2-6.0 mrad step 0.1).
#: D = 0.23 m from the ladder intercept x0 = -1.13 mm at delta = -0.0049,
#: a one-point number good to ~20%.  R34 is None for the same reason as
#: the baseline: the vertical plane is shut at 5 x 41.
#:
#: 10 x 100 is None because it CANNOT be checked: the only Z/A = 0.5 file
#: near that setting is `beamline_10x110_H2.xml` at 220 GV against the
#: 6Li fill's 199 GV, so a scan there would confound energy with lattice.
#: 18 x 275 carries the baseline triple unchanged, and that is not a
#: placeholder: the 18 x 275 compact file IS the Yellow Report 275 GeV
#: magnet set that the He4 file scales, which is why the He4 5 x 41 run
#: reproduces the 18 x 275 levers to 1%.
#:
#: Which file is right is the open question for the far-forward working
#: group (plans/09 B1); the baseline is what the published numbers use,
#: and the ratio of the two tables is the systematic to carry.
POT_LEVERS_LIGHT_ION_LATTICE = {
    "5x41":   (29.81, None, 0.23),
    "10x100": None,
    "18x275": POT_LEVERS["18x275"],
}

#: Second-order pot dispersion, x = D delta + D2 delta^2 [m].  Irrelevant
#: inside the near-beam band (0.6 mm at |delta| = 0.05) and NOT irrelevant
#: for an over-rigid fragment: at the triton's delta = 0.286 it is -18 mm
#: against a linear +84, and dropping it misses the measured +66 mm by a
#: quarter.  `over_rigid_route` uses it; nothing else does.
#:
#: D and D2 are one FIT and come from one station: `mkdisp.py` puts the
#: quadratic through the three rigidities each configuration provides at
#: theta = 0 -- alpha R = 0.857 and triton R = 1.286 from the frag guns,
#: intact 6Li R = 1.000 from the ladder's fitted intercept.  Station 1,
#: the plane R12 is regressed on, except at 5 x 41 where the R = 0.857
#: alpha reaches no station-1 silicon (the 32-48 mm band is held off to
#: |y| >= 18 mm there) and only station 2 has three points.  The two
#: stations disagree by 8% in D and up to 14% in D2 -- 0.292 / -0.215 at
#: station 1 against 0.315 / -0.249 at station 2 at 18 x 275 -- which is
#: real (station 2 has the longer lever) and is why the pair must not be
#: mixed across stations, as it was until 2026-08-28.
POT_DISPERSION_2 = {"5x41": -0.190, "10x100": -0.206, "18x275": -0.215}

#: |R - 1| over which the quadratic above was fitted and may be used.
MEASURED_DELTA_MAX = 0.30

#: Horizontal half-width of the silicon band a track at y ~ 0 must clear,
#: per configuration [m]: the central block is 16 mm wide and the bands
#: outside it are held off vertically by the per-energy insertion, so the
#: first band with a zero y offset starts at 48 / 32 / 16 mm
#: (compact/fields/beamline_*.xml).  The outer edge of the last module is
#: 144 mm at every configuration.
POT_BLIND_HALF_WIDTH = {"5x41": 0.048, "10x100": 0.032, "18x275": 0.016}
POT_OUTER_HALF_WIDTH = 0.144

#: Outer angular edge of the pot acceptance, MEASURED as the last theta
#: of a DEBRIS-FREE contiguous run: the run is broken by a gap of more
#: than 0.30 mrad and each row must carry at most 3 hits in the plane and
#: sit within 2 mm of the fitted transport line.  Taken as the smaller of
#: the +x and -x sides.  It is NOT 144 mm / R12 -- that arithmetic gives
#: 7.5 / 6.8 / 4.8 mrad and the ion strikes the pipe or the magnet
#: aperture first, at |dx| = 54-127 mm.  On the under-rigid (-x) side at
#: 5 x 41 the primary survives only 0.35 mrad past the inner edge.
#:
#: IT IS NOT THE LAST ROW WITH A PRIMARY.  Isolated on-line rows survive
#: 0.05-0.35 mrad beyond it, and the cleanliness cut rather than the
#: transport is what ends the run: at 5 x 41, phi = 180, theta = 2.90 mrad
#: the S1L1 hit is at dx = -55.0 mm against a fitted -55.6 (a primary) but
#: carries 13 hits in the plane; at 18 x 275, phi = 180, theta = 4.35 mrad
#: all four planes give a consistent on-line row (-130.8 / -131.2 /
#: -132.6 / -132.8 mm) past the 4.00 mrad quoted here.  The constant is
#: therefore the last angle at which a SINGLE clean track is reconstructed
#: through debris, which is the quantity a partner-fragment veto needs,
#: and it is conservative by up to one 0.05 mrad step plus one gap.
THETA_RP_OUTER_MEASURED = {"5x41": 2.85e-3, "10x100": 3.85e-3,
                           "18x275": 4.00e-3}

#: The 18 x 275 values under the old scalar names, so every caller and
#: every published number written before 2026-08-28 keeps working.  Note
#: that POT_R12 moves 30.6 -> 29.97 m (2%) because it is now regressed on
#: 138 rows of the current geometry rather than read off four hit
#: positions of the September-2024 one.
POT_R12, POT_R34, POT_DISPERSION = POT_LEVERS["18x275"]


def pot_levers_for(config):
    """(R12, R34, D) at a machine configuration.

    `config` may be a `beams.BeamConfig` of any species -- the transport
    is a property of the ring, not of the beam in it -- a configuration
    key ("5x41", "10x100", "18x275"), or a bare 6Li per-nucleon momentum
    in GeV/u, which resolves at the gamma-matched 40.8 / 99.5 / 137.5 and
    raises off them rather than interpolating.  R34 is None at 5 x 41; see
    POT_LEVERS.
    """
    if isinstance(config, str):
        return POT_LEVERS[config]
    if hasattr(config, "ion_momentum_per_nucleon"):
        return POT_LEVERS[yr_config_key(config)]
    from . import beams as _beams
    for cfg, key in zip(_beams.default_configs("6Li"),
                        ("5x41", "10x100", "18x275")):
        if abs(float(config) - cfg.ion_momentum_per_nucleon) < 1e-3:
            return POT_LEVERS[key]
    raise KeyError("no measured pot levers at %s GeV/u; the scan covers "
                   "the three 6Li configurations only "
                   "(tools/fullsim, plans/09 B1)" % config)


def theta_rp_outer_for(config):
    """The MEASURED outer angular edge of the pot acceptance [rad].

    Not the default of `route_charged`: THETA_RP_OUTER stays at
    THETA_RP_MAX so that no published acceptance moves under this
    measurement (`test_two_hit.py`).  Pass this value as `theta_outer` to
    price the real edge -- it is what the partner-fragment veto of
    plans/09 B4 depends on most.
    """
    if isinstance(config, str):
        return THETA_RP_OUTER_MEASURED[config]
    if hasattr(config, "ion_momentum_per_nucleon"):
        return THETA_RP_OUTER_MEASURED[yr_config_key(config)]
    from . import beams as _beams
    for cfg, key in zip(_beams.default_configs("6Li"),
                        ("5x41", "10x100", "18x275")):
        if abs(float(config) - cfg.ion_momentum_per_nucleon) < 1e-3:
            return THETA_RP_OUTER_MEASURED[key]
    raise KeyError("no measured outer edge at %s GeV/u" % config)


def over_rigid_route(R, theta_x=0.0, config="18x275"):
    """Does an over-rigid fragment (R > 1) land on Roman-pot silicon?

    Measured 2026-08-28 (tools/fullsim, plans/09 B1) and it is the answer
    to Report 3 Table 6's open assumption.  An over-rigid fragment bends
    LESS than the beam, so the pot dispersion carries it to +x, the inner
    side of the bend, and an R = 1.286 triton is on the pots in 60 of 60
    events at EVERY configuration -- dx = +66 mm at station 1, +70 to
    +72 mm at station 2, inside the 48-144 mm outer band, which carries no
    vertical insertion anywhere, so the hit needs no pT at all.  It
    deposits in the ZDC in 80 / 83 / 98% of the same events.  The
    repository routed every R > 1.05 fragment as lost.

    The test is on the pot-plane displacement, x = D delta + D2 delta^2 +
    R12 theta_x with delta = R - 1: silicon between the configuration's
    blind half-width (48 / 32 / 16 mm) and the last module at 144 mm.
    The second-order term is what makes the arithmetic reproduce the
    measurement -- at 18 x 275, D = 0.292 m and D2 = -0.215 m give
    0.292 x 286 - 0.215 x 286^2/1000 = 83.5 - 17.6 = 65.9 mm against
    66.3 measured.  `theta_x` is the signed horizontal IP angle; leaving
    it 0 asks the theta << dispersion question, which is the one the
    rigidity-window routing asks.  It also exposes the measured HOLE: at
    18 x 275 the triton's +66 mm is cancelled by theta_x in -1.7 to
    -2.7 mrad (this model, scanned in 0.01 mrad steps) and it disappears
    into the central block there.  The scan the model is fitted to loses
    the triton after -1.50 mrad and recovers it at -2.50, i.e. -1.55 to
    -2.53 mrad; the 0.15-0.2 mrad difference is the triton's own R12
    being 8% stiffer than the beam's, which the model does not carry.
    """
    r12, _r34, disp = POT_LEVERS[config]
    d2 = POT_DISPERSION_2[config]
    delta = np.asarray(R, dtype=float) - 1.0
    # The quadratic is fitted on delta = -0.143, 0, +0.286 and MUST NOT be
    # extrapolated: it turns over at delta = -D/(2 D2) ~ 0.6, which is a
    # three-point artefact and not a dispersion.  Outside the measured
    # range fall back on the linear term, which is monotone in delta as a
    # dispersion has to be.  The one fragment this decides is the 6Li
    # 3He + t triton at R = 1.5044: linear puts it at 152 mm, past the
    # last module at 144, so it stays "lost (over-rigid)" -- unmeasured
    # and routed as it always was.
    in_range = np.abs(delta) <= MEASURED_DELTA_MAX
    x = np.where(in_range, disp * delta + d2 * delta * delta, disp * delta)
    x = x + r12 * np.asarray(theta_x, dtype=float)
    ax = np.abs(x)
    return (ax >= POT_BLIND_HALF_WIDTH[config]) & (ax <= POT_OUTER_HALF_WIDTH)


def separation_at_pots(frag_a, frag_b, r12=POT_R12, r34=None,
                       dispersion=POT_DISPERSION, config=None):
    """Transverse distance [m] at the pot plane between two fragments that
    leave the IP at (theta, phi) carrying rigidity ratio R.

    Linear transport from a point source: horizontally the angular lever
    and the dispersion, vertically the lever alone,

        x = R12 theta cos phi + D (R - 1),      y = R34 theta sin phi,

    so each fragment dict needs `theta`, `phi` and `R` -- what
    `spectator._boost_fragment` returns.  `dispersion=0.0` recovers the
    angular lever alone.

    THE DISPERSIVE TERM DOES NOT CANCEL BETWEEN THE TWO FRAGMENTS, which
    is the correction of 2026-08-28 (code review of plans/09 B4).  This
    function was first written without it, on the argument that both
    fragments of 6Li -> alpha + d sit within 0.7% of beam rigidity so
    their dispersive displacements are common.  The 0.7% is true and the
    conclusion is false: the fragments take OPPOSITE rest-frame momenta,
    and the longitudinal component enters the two rigidities with opposite
    sign and unequal weight.  Measured over the sampled k and at every
    configuration, R_alpha = 0.9988 + 0.268 k_z and R_d = 1.0071 -
    0.536 k_z (k_z in GeV/c), so R_alpha - R_d = -0.0083 + 0.804 k_z and
    D (R_alpha - R_d) has 16/50/84% = 2.6 / 9.2 / 22.7 mm -- comparable to
    the angular separation itself, which it moves by 5 / 23 / 39% at
    5 x 41 / 10 x 100 / 18 x 275, and never below the 1.9 mm the k = 0
    rigidity difference alone gives.  It is also what makes a merge
    possible at all: the angular term alone puts a RECORDED pair tens of
    millimetres apart, while the dispersive term can cancel it.

    The SIGN of D relative to R12 is not measured, and does not matter
    here: the joint distribution of (angular dx, R_a - R_b) is even in the
    first (the breakup azimuth is uniform), so the separation distribution
    is invariant under D -> -D.

    Pass `config` (a BeamConfig, a key, or a 6Li GeV/u) to use the
    per-configuration levers `POT_LEVERS` measured on 2026-08-28, which is
    what plans/09 D3 asked for and what every caller should now do; the
    bare defaults keep the 18 x 275 numbers, which is where R12 = 30.6 m
    and D = 0.30 m came from and where they were applied at EVERY
    configuration until then.  The correction is large and it runs the
    unintuitive way: R12 is 19.2 and 21.2 m at 5 x 41 and 10 x 100, a
    third smaller than 30.0 m at the top, so the millimetres the lower
    configurations carried were 45-60% HIGH.  R34 is now measured at two
    of the three -- 3.35 and 2.93 m, an order of magnitude under R12, the
    far-forward line being ~9x stiffer in x than in y -- and is None at
    5 x 41, where the 29.6 mm insertion shuts the vertical plane
    altogether and there is nothing to regress; there the fallback r34 =
    r12 remains, and it is the one place a millimetre in this function is
    still an assumption.
    Nothing in the ROUTING depends on any of it -- the acceptance is
    decided in angle (`route_charged`), which is the same assumption from
    the other side (plans/04 #11: a near-beam fragment is taken to be
    dispersion-blind, though its own D (R - 1) is 0.6 mm for the alpha and
    1.4 mm for the deuteron) -- so this enters figures and tables only.
    """
    if config is not None:
        r12, r34_m, dispersion = pot_levers_for(config)
        r34 = r12 if r34_m is None else r34_m
    r34 = r12 if r34 is None else r34
    ta, tb = np.asarray(frag_a["theta"]), np.asarray(frag_b["theta"])
    pa, pb = np.asarray(frag_a["phi"]), np.asarray(frag_b["phi"])
    dx = r12 * (ta * np.cos(pa) - tb * np.cos(pb))
    if dispersion:
        dx = dx + dispersion * (np.asarray(frag_a["R"])
                                - np.asarray(frag_b["R"]))
    dy = r34 * (ta * np.sin(pa) - tb * np.sin(pb))
    return np.sqrt(dx * dx + dy * dy)


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
                         per_config_levers=False,
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
    # The pot dispersion turns the beam's momentum spread into an apparent
    # ANGLE at the IP, D dp/p / R12.
    #
    # MEASURED PER CONFIGURATION SINCE 2026-08-28 (POT_LEVERS, plans/09
    # B1) AND NOT YET APPLIED HERE BY DEFAULT.  D/R12 is 1.62e-2 /
    # 1.35e-2 / 9.74e-3 m/m per configuration.  READ THE DEFAULT BRANCH
    # CAREFULLY: POT_R12 and POT_DISPERSION are now POT_LEVERS["18x275"],
    # so the default is 0.292 / 29.97 = 9.74e-3, i.e. the MEASURED
    # 18 x 275 pair -- the third of those three numbers, not a legacy one.
    # It lands 0.6% BELOW the 0.30 / 30.6 = 9.80e-3 that Report 1 SS6.1,
    # the reproduction manual and every published tagging number were
    # priced with, so no published tagging number moves under the
    # re-measurement: the 18 x 275 optimum is r_h = 89.3, eps = 0.332,
    # L/L_HA = 1/9.5 before and after.  That agreement is a result, not a
    # coincidence to be relied on -- the two levers each moved by 2% and
    # their ratio did not.
    #
    # What `per_config_levers=True` buys is therefore the two LOWER
    # configurations, where the smearing is 66% and 39% larger than the
    # 18 x 275 default.  It moves the optimum measurably -- r_h 49.7 /
    # 175.6 / 89.3 -> 46.5 / 164.1 / 89.3, eps 0.423 / 0.323 / 0.332 ->
    # 0.374 / 0.251 / 0.332, L/L_HA 1/7.1 / 1/13.3 / 1/9.5 -> 1/6.8 /
    # 1/12.8 / 1/9.5, identical at 18 x 275 by construction -- and the
    # 22% drop in eps at 10 x 100 propagates through
    # `money_cos2phi_coherent_reco.py` into Reports 1, 3 and 4.  It is
    # therefore OPT-IN until the coherent chain is re-run with it, which
    # is the one open item plans/09 B1 hands on rather than closes.
    if per_config_levers:
        r12_c, _r34_c, disp_c = pot_levers_for(config)
    else:
        r12_c, disp_c = POT_R12, POT_DISPERSION
    disp = (disp_c * dpp / r12_c) if dispersion else 0.0
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
