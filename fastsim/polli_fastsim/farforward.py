"""Far-forward detector windows and fragment routing (Phase-1 stand-in).

Windows fetch-verified 2026-06-12 (YR detector matrix; arXiv:2108.08314
Table I; arXiv:2409.02811; details in plans/03 step 2.2):

  Roman Pots   z = 26/28 m   theta < 5 mrad   R in [0.60, 0.95]
               near-beam (|R-1| < 0.05): only theta > n_sigma sigma_theta
  OMD          z = 22.5/24.5 m  theta < 5 mrad  R in [0.45, 0.65], no cut
  B0           5.5 < theta < 20 mrad, any charged
  ZDC          theta < 4 mrad, neutrals
  no coverage  R > ~1.05 (bends less than beam), or 5-5.5 mrad gap

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

    `sigma_theta` is the beam angular divergence at the IP; a fragment at
    beam rigidity clears the envelope when theta > n_sigma sigma_theta.
    """
    name: str
    sigma_theta: float          # rad
    n_sigma: float = 10.0

    @property
    def pt_cut_near_beam(self):
        """The same cut expressed for a 275 GeV proton -- the number the
        Yellow Report quotes (0.20 / 0.41 GeV).  Use `pt_cut_for` for a
        fragment of any other momentum."""
        return self.pt_cut_for(PROTON_REFERENCE_MOMENTUM)

    def pt_cut_for(self, momentum):
        """Transverse-momentum threshold for a fragment of total momentum
        `momentum` [GeV] -- e.g. A * p_per_nucleon for a nucleus."""
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


def route_charged(R, theta, pT, optics=HIGH_ACCEPTANCE):
    """Classify charged fragments into far-forward systems.

    Returns an integer array: 0 lost, 1 RP, 2 OMD, 3 B0, 4 RP-near-beam
    (R ~ 1, accepted only outside the angular envelope
    theta > n_sigma sigma_theta).  `pT` is no longer used for the
    near-beam decision -- the envelope is angular, and theta is exact per
    fragment where a pT threshold needs the fragment momentum.
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
    rp_tail = small & near & (theta > optics.n_sigma * optics.sigma_theta)
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


def acceptance_summary(R, theta, pT, optics=HIGH_ACCEPTANCE):
    """Fraction of spectators in each far-forward system."""
    route = route_charged(R, theta, pT, optics)
    n = float(len(route))
    return {label: float(np.sum(route == code)) / n
            for code, label in ROUTE_LABELS.items()}
