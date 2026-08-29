#!/usr/bin/env python3
"""Can a superconducting nanowire do the near-beam job (plans/09)?

nearbeam_aperture_scan.py says what a closer approach is worth.  This
script asks whether the Argonne-style NbN device can be the thing that
delivers it, and whether the same device can answer open question #19
(charge identification for A/Z = 2 fragments), on four axes:

  1. SIGNAL -- what a relativistic 6Li, alpha or deuteron deposits in a
     12 nm NbN film, against the one-photon threshold.
  2. Z-ID -- the hot-spot threshold model.  A current-biased nanowire is
     NOT proportional: the pulse amplitude is the diverted bias current
     and is the same for hadrons, muons, pions and showering electrons
     (measured, arXiv:2510.11725 and 2410.00251).  What the deposit sets
     is the BIAS CURRENT at which the wire fires,

         r_s = sqrt(Q / (e pi c rho (Tc - T0)))        (hot-spot radius)
         I_th / I_c = 1 - 2 r_s / w                    (firing threshold)

     [Lee et al. arXiv:2312.13405 Eqs. 1-2; the same relation appears in
     Renema et al. arXiv:1301.3337 as E = (w/C)^2 (1 - I_b/I_c)^2, and in
     the ion literature as I_th/I_c = 1 - (zeV)^(1/2) C/w, Cristiano et
     al., SUST 28 (2015) 124004.]  At fixed beta, Q goes as z^2, so r_s
     goes as z LINEARLY and the species have separated turn-on curves.
     The scheme is a granted patent (US 8,872,109, Ohkubo & Suzuki),
     demonstrated on singly- vs doubly-charged lysozyme by bias-point
     subtraction alone.
  3. AREA -- how much sensor the near-beam region needs, in devices of
     each available granularity.
  4. The alpha tag, per optics.

The anchor is ANL's own MEASURED hot-spot radius, r_s = 134 nm for a
120 GeV proton in this film.  A 6Li at 137.5 GeV/u has beta*gamma = 148
against 128 for that proton -- the same velocity -- so the step to z = 3
is pure z with no velocity correction.  ANL have also run a 5.5 MeV
alpha, which their own scaling puts at ~1 um: the 6Li lands BETWEEN
their two measured points.  This is an interpolation, not a leap.

Usage:  python3 scripts/nearbeam_sensor_budget.py
"""

import argparse
import math
import pathlib
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_SCRIPTS = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS.parent))
sys.path.insert(0, str(_SCRIPTS.parent.parent / "fastsim"))

from polligen import nearbeam as nb  # noqa: E402
from polligen import reco  # noqa: E402
from polli_fastsim import beams  # noqa: E402
from polli_fastsim import farforward as ff  # noqa: E402
from polli_fastsim import spectator as sp  # noqa: E402

# --- the devices, from the papers (the FILM physics is polligen.nearbeam) ---
R_S_PROTON_NM = nb.R_S_PROTON_NM   # MEASURED, 120 GeV p in 12 nm NbN
W_MIP_OPTIMAL_NM = nb.W_MIP_OPTIMAL_NM
DARK_COUNT_WALL = nb.DARK_COUNT_WALL
SNSPD_THRESHOLD_EV = nb.SNSPD_THRESHOLD_EV
DEVICE_UM = 30.0               # nanowire active area per device, square
SMSPD_MM = 1.0                 # arXiv:2510.11725, 8 pixels per 1 x 1 mm^2
SMSPD_W_NM = 1000.0            # its wire width
R12_M = 29.97                  # IP angle -> pot-plane x, MEASURED at 18x275
                               # (2026-08-28, farforward.POT_LEVERS; 19.24
                               # and 21.25 m at the two lower ones)

SPECIES = (("6Li", 3), ("alpha", 2), ("d", 1), ("p", 1))
WIRE_WIDTHS_NM = (250.0, 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0)
# Beam configurations are DERIVED from polli_fastsim.beams, never
# hard-coded: the two lower 6Li energies moved from rigidity-scaled
# (20.5, 50) to gamma-matched (40.8, 99.5) GeV/u on 2026-08-27 (plans/10).
_CFG = beams.default_configs("6Li")
CONFIGS = tuple((n, c.ion_momentum_per_nucleon)
                for n, c in zip(("5 x 41", "10 x 100", "18 x 275"), _CFG))
# Threshold bias read off Fig. 7 of arXiv:2312.13405 for 120 GeV protons.
# A DIGITISATION of a published figure, not a tabulated number.
ANL_FIG7_NM_RATIO = ((300.0, 0.215), (400.0, 0.49), (600.0, 0.62),
                     (800.0, 0.70))
C_LI, C_AL, C_D, C_GREY = "#1F4E79", "#C0392B", "#B8860B", "#8A8A8A"

threshold_ratio = nb.threshold_ratio


def gain_profile(theta_new, theta_old, p_ion, slope_b, theta_y, fracs,
                 n=400000, rng=None):
    """Of the coherent recoils a closer approach NEWLY accepts, how far
    past the new edge must the sensor reach to collect each fraction?

    The gain itself is analytic (reco.rp_hole_acceptance); only the shape
    needs sampling, and a plain exp(-B t) throw finds nothing where the
    acceptance is 1e-17, so t is sampled from the exponential shifted up
    to t_min = (theta_new p)^2 -- exactly the conditional distribution."""
    rng = rng or np.random.default_rng(11)
    gain = (reco.rp_hole_acceptance(slope_b, theta_new * p_ion,
                                    theta_y * p_ion)["acc"]
            - reco.rp_hole_acceptance(slope_b, theta_old * p_ion,
                                      theta_y * p_ion)["acc"])
    t = (theta_new * p_ion) ** 2 + rng.exponential(1.0 / slope_b, n)
    pt = np.sqrt(t)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    thx = np.abs(pt * np.cos(phi)) / p_ion
    gained = (thx > theta_new) & (thx < theta_old)
    if not gained.any():
        return gain, {f: float("nan") for f in fracs}
    d = np.sort(thx[gained] - theta_new)
    return gain, {f: float(d[min(int(f * d.size), d.size - 1)]) for f in fracs}


def alpha_edge_profile(p_per_nucleon, theta_new, fracs, n=200000, rng=None):
    """The same sizing question for the 6Li alpha spectator."""
    k = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG, p_per_nucleon, n,
                                    beta=0.30,
                                    rng=rng or np.random.default_rng(7))
    d = np.sort(k["theta"][k["theta"] > theta_new] - theta_new)
    if not d.size:
        return {f: float("nan") for f in fracs}
    return {f: float(d[min(int(f * d.size), d.size - 1)]) for f in fracs}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slope-b", type=float, default=50.0)
    ap.add_argument("--near-beam-mrad", type=float, default=None,
                    help="horizontal half-width the near-beam layer reaches "
                         "[mrad]; default: the 10 sigma envelope of the "
                         "tagging optics of each configuration "
                         "(farforward.tagging_optics_point; 2026-08-28), "
                         "which is what the layer exists to reach")
    ap.add_argument("--vertical-mm", type=float, default=20.0,
                    help="vertical extent of a near-beam insert")
    ap.add_argument("--n-spectator", type=int, default=200000)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    beta_gamma = 137.5 / 0.9315               # 6Li at top rigidity, beta ~ 1

    def near_beam(cfg):
        """The horizontal half-width the layer reaches [rad]."""
        if args.near_beam_mrad is not None:
            return 1e-3 * args.near_beam_mrad
        return ff.tagging_optics_point(cfg, slope_b=args.slope_b)["env_x"]

    print("== 1. signal in one %s film (beta*gamma = %.0f) =="
          % (nb.NBN.name, beta_gamma))
    print("%-6s %3s %12s %12s %10s"
          % ("", "z", "<dE> [eV]", "photon-equiv", "r_s [nm]"))
    dep = {}
    for name, z in SPECIES:
        q = nb.bethe_mean_ev(z, beta_gamma)
        dep[name] = (q, nb.hot_spot_nm(z))
        print("%-6s %3d %12.1f %12.0f %10.0f"
              % (name, z, q, q / SNSPD_THRESHOLD_EV, nb.hot_spot_nm(z)))
    print("Every species is far above the one-photon threshold, so "
          "DETECTION is not the question; the deposit ratios are z^2 "
          "exactly, because 6Li at 137.5 GeV/u and ANL's 120 GeV "
          "calibration proton have the same velocity to 15%.")
    print("The Landau MPV is deliberately NOT quoted: xi/I = %.1e in this "
          "film, far outside the xi >> I regime the formalism needs."
          % nb.NBN.xi_over_i())

    print("\n== 2. Z-ID: firing threshold I_th/I_c = 1 - 2 r_s / w ==")
    print("anchored on ANL's MEASURED r_s = %.0f nm for a 120 GeV proton "
          "in this film.  0.00 means the wire fires at ANY bias and "
          "carries no charge information." % R_S_PROTON_NM)
    print("%9s %8s %8s %8s   %s" % ("w [nm]", "p / d", "alpha", "6Li",
                                    "separated turn-ons"))
    for w in WIRE_WIDTHS_NM:
        th = {n: threshold_ratio(dep[n][1], w) for n, _ in SPECIES}
        live = [n for n in ("6Li", "alpha", "d") if th[n] > 0.0]
        print("%9.0f %8.2f %8.2f %8.2f   %s"
              % (w, th["p"], th["alpha"], th["6Li"],
                 "none" if len(live) < 2 else
                 "%s (%d)" % (" > ".join(live), len(live))))
    print("ANL's MIP optimum, w = %.0f nm, is 2 r_s for a proton: maximal "
          "efficiency and ZERO charge information.  Z-ID needs "
          "DELIBERATELY WIDE wires -- which is the microwire (SMSPD) "
          "geometry that already exists (w = %.0f nm, arXiv:2510.11725)."
          % (W_MIP_OPTIMAL_NM, SMSPD_W_NM))
    w = SMSPD_W_NM
    b_li = threshold_ratio(dep["6Li"][1], w)
    b_al = threshold_ratio(dep["alpha"][1], w)
    b_p = threshold_ratio(dep["p"][1], w)
    print("At w = %.0f nm the turn-ons are %.2f (6Li), %.2f (alpha), %.2f "
          "(p/d) I_c.  Two bias points at %.2f and %.2f I_c tag Z by the "
          "FIRING PATTERN alone, and both sit below the %.2f I_c "
          "dark-count wall."
          % (w, b_li, b_al, b_p, 0.5 * (b_li + b_al), 0.5 * (b_al + b_p),
             DARK_COUNT_WALL))

    print("\n== 3. sizing: where a closer approach's gain lives ==")
    print("(the layer reaches the tagging optics' 10 sigma envelope from the "
          "measured silicon aperture; vertical = the larger of the measured "
          "aperture and the high-acceptance 10 sigma_v)")
    print("%-9s %8s %10s %11s %11s %9s %9s"
          % ("optics", "p_ion", "to [mrad]", "gain", "d50 [urad]", "d90", "d99"))
    prof = {}
    for name, pu in CONFIGS:
        cfg = beams.default_configs("6Li")[[c[0] for c in CONFIGS].index(name)]
        p_ion = 6.0 * pu
        meas = reco.RP_APERTURE_MEASURED[name.replace(" ", "")]
        th_nb = near_beam(cfg)
        ty = max(meas[1], ff.yr_optics(cfg).envelope[1])
        g, d = gain_profile(th_nb, meas[0], p_ion, args.slope_b, ty,
                            (0.50, 0.90, 0.99))
        prof[name] = (g, d)
        print("%-9s %8.1f %10.2f %11.3e %11.0f %9.0f %9.0f"
              % (name, p_ion, 1e3 * th_nb, g, 1e6 * d[0.50], 1e6 * d[0.90],
                 1e6 * d[0.99]))
    th_nb = near_beam(beams.default_configs("6Li")[2])
    a_prof = alpha_edge_profile(137.5, th_nb, (0.50, 0.90, 0.99),
                                n=args.n_spectator)
    print("6Li alpha spectator at 137.5 GeV/u: d50 %.0f, d90 %.0f, d99 %.0f "
          "urad" % (1e6 * a_prof[0.50], 1e6 * a_prof[0.90],
                    1e6 * a_prof[0.99]))
    print("The exponential piles the gain against the inner edge: a "
          "near-beam insert is a NARROW STRIP, not a plane.")

    print("\n== 3b. what a strip costs in channels ==")
    for span_mm, label in (
            # 9.3 mm is (1.03 - 0.727) mrad x 30.6 m: the Sep-2024
            # silicon aperture down to the RETIRED single 72.7 urad
            # envelope.  Kept so the published channel count stays
            # reproducible; the row below is its current-geometry
            # replacement, against the tagging envelope.
            (9.3, "9.3 mm, the Sep-2024 geometry's gap at 18x275 "
                  "(R12 = %.1f m)" % R12_M),
            (12.4, "12.4 mm, the gap in the CURRENT geometry at 18x275 -- "
                   "the measured 0.53 mrad silicon aperture down to the "
                   "0.12 mrad tagging envelope at R12 = %.1f m" % R12_M),
            (3.0, "3.0 mm, a minimal near-beam strip")):
        area = 2.0 * span_mm * args.vertical_mm       # two sides of the beam
        print("  %s -> %.0f mm2 total" % (label, area))
        print("      %8.0f nanowire devices (%.0f x %.0f um2)"
              % (area / (1e-3 * DEVICE_UM) ** 2, DEVICE_UM, DEVICE_UM))
        print("      %8.0f microwire tiles = %.0f channels (%.0f x %.0f "
              "mm2, 8 px)" % (area / SMSPD_MM ** 2, 8.0 * area / SMSPD_MM ** 2,
                              SMSPD_MM, SMSPD_MM))
    print("The nanowire granularity cannot tile a far-forward aperture; "
          "the microwire can.  The area answer and the Z-ID answer are "
          "the SAME device.")

    print("\n== 4. the alpha tag, per optics (routed through the far-forward "
          "systems; rectangular envelope, vertical as above) ==")
    print("%-9s %13s %14s %14s %8s" % ("optics", "silicon", "YR envelope",
                                       "tagging optics", "gain"))
    for name, pu in CONFIGS:
        cfg = beams.default_configs("6Li")[[c[0] for c in CONFIGS].index(name)]
        meas = reco.RP_APERTURE_MEASURED[name.replace(" ", "")]
        env = ff.yr_optics(cfg).envelope
        ty = max(meas[1], env[1])
        th_nb = near_beam(cfg)
        k = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG, pu,
                                        args.n_spectator, beta=0.30,
                                        rng=np.random.default_rng(7))
        # the routed tag (any far-forward system), the same quantity as
        # tagging_acceptance.py and nearbeam_aperture_scan.py.  pot_config
        # is the configuration whose blind block the OVER-RIGID branch
        # tests against (48 / 32 / 16 mm); leaving it at the 18 x 275
        # default inflates the two lower rows by 50%.
        key = name.replace(" ", "")
        acc = lambda hx: 1.0 - ff.acceptance_summary(                # noqa: E731
            k["R"], k["theta"], k["pT"], ff.Optics("cut", hx / 10.0, 10.0, ty / 10.0),
            phi=k["phi"], pot_config=key)["lost"]
        a_old, a_env, a_new = acc(meas[0]), acc(env[0]), acc(th_nb)
        print("%-9s %13.4f %14.4f %14.4f %8.2f"
              % (name, a_old, a_env, a_new, a_new / max(a_old, 1e-12)))

    zid_figure(dep, args.outdir)


def zid_figure(dep, outdir):
    """I_th/I_c against wire width for the three species, with ANL's own
    measured proton points on top."""
    w = np.linspace(150.0, 2200.0, 400)
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    for name, label, col in (("6Li", "$^{6}$Li  (z = 3)", C_LI),
                             ("alpha", r"$\alpha$  (z = 2)", C_AL),
                             ("p", "p, d  (z = 1)", C_D)):
        ax.plot(w, [threshold_ratio(dep[name][1], x) for x in w], "-",
                color=col, lw=1.8, label=label)
        ax.axvline(2.0 * dep[name][1], color=col, lw=0.9, ls=":", alpha=0.7)
    ax.plot([p[0] for p in ANL_FIG7_NM_RATIO],
            [p[1] for p in ANL_FIG7_NM_RATIO], "o", color=C_D, ms=6,
            mfc="white", mew=1.6,
            label="measured, 120 GeV p\n(arXiv:2312.13405 Fig. 7)")
    ax.axhspan(DARK_COUNT_WALL, 1.0, color=C_GREY, alpha=0.16, lw=0)
    ax.text(2150, 0.5 * (1.0 + DARK_COUNT_WALL), "dark-count wall",
            fontsize=7.5, color=C_GREY, ha="right", va="center")
    ax.axvline(W_MIP_OPTIMAL_NM, color="k", lw=1.0, ls="--", alpha=0.6)
    ax.annotate("ANL MIP optimum\n(no charge information)",
                xy=(W_MIP_OPTIMAL_NM, 0.93), xytext=(6, 0),
                textcoords="offset points", fontsize=7.2, va="top")
    ax.axvline(SMSPD_W_NM, color="#1B7F5B", lw=1.0, ls="--", alpha=0.8)
    ax.annotate("existing microwire\n(three separated turn-ons)",
                xy=(SMSPD_W_NM, 0.06), xytext=(6, 0),
                textcoords="offset points", fontsize=7.2, va="bottom",
                color="#1B7F5B")
    ax.set_xlim(150, 2200)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("wire width w [nm]")
    ax.set_ylabel(r"firing threshold $I_{\rm th}/I_c$")
    ax.set_title("Charge identification by firing threshold: "
                 r"$I_{\rm th}/I_c = 1 - 2r_s/w$, $r_s \propto z$",
                 fontsize=10)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(fontsize=7.5, loc="lower right")
    fig.tight_layout()
    out = pathlib.Path(outdir) / "nearbeam_zid_threshold.png"
    fig.savefig(out, dpi=140)
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
