"""Coherent (intact-ground-state) e+6Li channel: scenario rates, recoil
tagging, and cos 2phi tensor-modulation projections.

Process: e + 6Li -> e' + X + 6Li(g.s.) -- coherent diffractive DIS, the
nucleus stays in its 1+ ground state and recoils with |t| ~= pT^2.  The
intact 6Li has A/Z = 2, exactly the beam rigidity (R = 1.000 to O(x_P),
and the diffractive momentum loss x_P < ~0.02 keeps it deep inside the
|R-1| < 0.05 near-beam band), so at IP6 the ONLY far-forward handle is
the Roman-Pot near-beam pT tail: the same beam-blindness that limits the
6Li alpha-tag (plans/00).  Under the exponential coherent t-slope the
tagging acceptance is analytic, acc = exp(-B pT_cut^2).

Everything here is a SCENARIO model in the sense of polarized.py: the
coherent fraction, the slope, and the tensor modulation amplitude are
placeholders with explicit bands, to be replaced by a real diffractive
model (nuclear shadowing / dipole) when one is adopted -- plans/06.

Conventions:
* Rates reuse fom.project_rates (per-nucleon DIS rate with the scattered-
  electron acceptance cuts already applied), scaled by the coherent
  fraction f_coh(x).  f_coh is written against Bjorken x while coherence
  is really controlled by x_P >= x, so f_coh(x) is optimistic for
  high-mass diffraction; the f0 band is meant to cover this.
* The cos 2phi' modulation of the coherent yield (phi' = azimuth relative
  to the transverse alignment axis) is quoted at P_zz = 1; a fill with
  tensor polarization P_zz measures amp * P_zz, exactly like the
  inclusive gluonometry convention (xsec.py).
"""

from dataclasses import dataclass

import numpy as np

from polli_fastsim import fom
from polli_fastsim.farforward import (NEAR_BEAM_BAND, OMD_R_WINDOW,
                                      RP_R_WINDOW, route_charged)

M_LI6 = 5.6015  # GeV (6.0151228 u atomic minus 3 m_e)
GEV_PER_FM_INV = 0.19733  # hbar c


def gaussian_slope(r_rms_fm):
    """Coherent |F(t)|^2 t-slope B [GeV^-2] for a Gaussian density:
    F(t) = exp(-R_rms^2 |t| / 6) -> |F|^2 = exp(-B |t|), B = R_rms^2/3."""
    return (r_rms_fm / GEV_PER_FM_INV) ** 2 / 3.0


@dataclass(frozen=True)
class CoherentScenario:
    """Scenario parameters for coherent diffractive e+6Li (SCENARIO).

    f0:      coherent fraction of the DIS rate at x -> 0.  HERA ep sees
             ~10% diffraction at low x of which the elastic-proton part
             dominates; for a loosely bound A = 6 nucleus the coherent
             share is smaller -- default 0.04 with a x2 band {0.02, 0.08}.
    x_coh:   coherence falloff scale: l_c ~ (hbar c)/(2 m_N x) = 0.105/x fm
             drops below the 6Li diameter (~5 fm) at x ~ 0.02.
    slope_b: |F(t)|^2 slope [GeV^-2], B = <r^2>/3 (gaussian_slope):
             matter radius 2.32-2.45 fm -> 45-51; charge radius
             2.589(39) fm -> ~57 (Sick 2015); the band {40, 60} covers
             both, and STAR d+Au (PRL 128:122303) suggests gluonic
             distributions can be wider than charge.  Exponential form
             valid below the first form-factor minimum at |t| ~ 0.31
             GeV^2 (q^2 = 8 fm^-2, Li-Sick-Whitney-Yearian 1971).
    amp:     cos 2phi' modulation of the coherent yield at P_zz = 1
             (gluonic-quadrupole / double-helicity-flip scenario).
    """
    f0: float = 0.04
    x_coh: float = 0.01
    slope_b: float = 50.0
    amp: float = 0.02

    def coherent_fraction(self, x):
        """f_coh(x) = f0 / (1 + (x/x_coh)^2): flat well below x_coh, dies
        once the coherence length drops below the nuclear size."""
        x = np.asarray(x, dtype=float)
        return self.f0 / (1.0 + (x / self.x_coh) ** 2)

    def sample_t(self, n, rng, t_max=0.5):
        """|t| ~ B exp(-B |t|), truncated at t_max [GeV^2]."""
        u = rng.uniform(size=n)
        c = 1.0 - np.exp(-self.slope_b * t_max)
        return -np.log(1.0 - c * u) / self.slope_b

    def tag_acceptance(self, pt_cut):
        """Fraction of coherent recoils above the near-beam pT cut:
        P(pT > cut) = exp(-B cut^2) for dsigma/dt ~ exp(-B pT^2)."""
        return float(np.exp(-self.slope_b * pt_cut * pt_cut))

    def mean_t_tagged(self, pt_cut):
        """<|t|> of the TAGGED sample: pT_cut^2 + 1/B (exponential tail)."""
        return pt_cut * pt_cut + 1.0 / self.slope_b


def recoil_lab(t_abs, phi_t, p_per_nucleon, x_pom=0.0):
    """Lab kinematics of the intact 6Li recoil.

    pT = sqrt(|t|) (t_min ~ (x_P M_A)^2/(1-x_P) is < 1e-3 GeV^2 for
    x_P < 0.005 and is neglected); the longitudinal momentum keeps
    (1 - x_pom) of the beam value, so R = p/Z over beam rigidity ~ 1.
    Returns dict with pT, theta, R, xL (per-nucleus fraction).
    """
    t_abs = np.asarray(t_abs, dtype=float)
    pt = np.sqrt(t_abs)
    p_beam = 6.0 * p_per_nucleon
    pz = (1.0 - np.asarray(x_pom, dtype=float)) * p_beam
    p_lab = np.hypot(pt, pz)
    theta = np.arctan2(pt, pz)
    r = p_lab / p_beam  # same Z as beam: rigidity ratio = momentum ratio
    return {"pT": pt, "theta": theta, "R": r, "xL": pz / p_beam,
            "phi_t": np.asarray(phi_t, dtype=float)}


def tag_acceptance_sampled(scenario, optics, p_per_nucleon, n=200000,
                           rng=None, x_pom=0.0):
    """Monte-Carlo cross-check of tag_acceptance through the actual
    far-forward router (route 4 = RP near-beam pT tail)."""
    rng = rng or np.random.default_rng(20260713)
    t = scenario.sample_t(n, rng)
    phi_t = rng.uniform(0.0, 2.0 * np.pi, size=n)
    lab = recoil_lab(t, phi_t, p_per_nucleon, x_pom=x_pom)
    route = route_charged(lab["R"], lab["theta"], lab["pT"], optics)
    return float(np.mean(route == 4))


def project_coherent(config, scenario_fom, coh, optics_list=(), nx=40,
                     nq2=30, x_range=(1e-4, 1.0), q2_range=(1.0, 2e3),
                     nuclear_f2=None):
    """Coherent event counts per (x, Q2) bin.

    Returns (proj, n_coh, tagged) where proj is the underlying DIS
    BinnedProjection (e' acceptance applied), n_coh = DIS counts times
    f_coh(x), and tagged maps each Optics to n_coh * exp(-B pT_cut^2).
    """
    proj = fom.project_rates(config, scenario_fom, nx=nx, nq2=nq2,
                             x_range=x_range, q2_range=q2_range,
                             nuclear_f2=nuclear_f2)
    n_coh = proj.n_events * coh.coherent_fraction(proj.x)
    tagged = {opt.name: n_coh * coh.tag_acceptance(opt.pt_cut_near_beam)
              for opt in optics_list}
    return proj, n_coh, tagged


# --- breakup-channel bookkeeping for the incoherent-background note -------

# Lowest-lying particle decompositions of 6Li with separation energies
# [MeV] (TUNL A=6 evaluation, Tilley et al. NPA 708:3: S(alpha d)=1.4743,
# alpha+p+n = S(alpha d)+S(d) = 3.699, 3He+t = 15.795; the unbound
# 5He/5Li intermediates, thresholds 4.50/5.39, end in these channels).
LI6_BREAKUP = (
    ("alpha+d", 1.474, (("alpha", 4, 2), ("d", 2, 1))),
    ("alpha+p+n", 3.70, (("alpha", 4, 2), ("p", 1, 1), ("n", 1, 0))),
    ("3He+t", 15.79, (("3He", 3, 2), ("t", 3, 1))),
)


def fragment_rigidity(a, z, beam_a=6, beam_z=3):
    """Rigidity ratio R of a beam-velocity fragment: (A/Z)/(A/Z)_beam."""
    if z == 0:
        return np.nan
    return (a / z) / (beam_a / beam_z)


def fragment_route_label(a, z):
    """Far-forward destination of a beam-velocity fragment at IP6, from
    the rigidity windows alone (theta << 5 mrad for Fermi-motion pT)."""
    if z == 0:
        return "ZDC"
    r = fragment_rigidity(a, z)
    if abs(r - 1.0) < NEAR_BEAM_BAND:
        return "RP pT-tail only (beam-blind)"
    if RP_R_WINDOW[0] <= r <= RP_R_WINDOW[1]:
        return "RomanPots"
    if OMD_R_WINDOW[0] <= r < OMD_R_WINDOW[1]:
        return "OMD"
    if r > 1.0 + NEAR_BEAM_BAND:
        return "lost (over-rigid)"
    return "lost"


def veto_table():
    """{channel: [(fragment, R, destination), ...]} for the incoherent-
    background veto discussion (plans/06)."""
    out = {}
    for name, _thr, frags in LI6_BREAKUP:
        out[name] = [(fname, fragment_rigidity(a, z),
                      fragment_route_label(a, z))
                     for fname, a, z in frags]
    return out
