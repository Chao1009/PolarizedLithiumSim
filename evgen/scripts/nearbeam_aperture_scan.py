#!/usr/bin/env python3
"""What the far-forward near-beam aperture is worth (plans/09).

The coherent intact-6Li tag and the 6Li alpha tag are both beam-blind:
A/Z = 2 exactly for the intact nucleus and a rigidity ratio of 0.99813
for the alpha, so neither is separated from the beam by dispersion.  What
separates them is the transverse ANGLE at the interaction point, and the
acceptance is therefore exp(-B t) evaluated at t = (theta p)^2 -- an
exponential in the SQUARE of the angle times the square of the momentum.
Approach distance is the whole measurement.

Two angles matter and they are not the same number:

  * the 10 sigma beam envelope, 0.727 mrad for the high-acceptance
    optics -- an operational retraction, and a machine-protection
    question rather than a detector one;
  * the ePIC Roman-Pot silicon APERTURE, measured by shooting an intact
    6Li through the geometry (tools/fullsim): |theta_x| > 1.03 mrad at
    18 x 275, 1.35 at 10 x 100, 2.0 at 5 x 41.  It is the sensor package
    -- the 1 mm aluminium RF shield, the ASIC, the thermal strips, and
    the 32 mm module granularity -- not the beam.

    CAVEAT, 2026-08-26: that measurement was made in the September-2024
    epic-main geometry, and the pot layout has since changed -- 16 mm
    modules with per-energy 10 sigma insertions replace 32 mm modules
    with a single energy-independent one, and the 1 mm aluminium RF
    shields are commented out.  The marked points are therefore a
    SUPERSEDED anchor.  The CURVES are not: they price any aperture, and
    that is the point of plotting them.  plans/09 D2/B1.

The aperture is the binding constraint at every configuration, and this
script says what closing the gap between the two would buy.  It is the
quantitative case for a near-beam layer thin enough to reach the
envelope, which is what plans/09 examines.

Usage:  python3 scripts/nearbeam_aperture_scan.py [--outdir .]
"""

import argparse
import pathlib
import sys

import numpy as np

_SCRIPTS = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS.parent))
sys.path.insert(0, str(_SCRIPTS.parent.parent / "fastsim"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from polligen import reco  # noqa: E402
from polli_fastsim import spectator as sp  # noqa: E402

C_TRUTH, C_FIT, C_ALT, C_GREY = "#1F4E79", "#C0392B", "#B8860B", "#8A8A8A"

# 6Li at the reference rigidity of each ring, and the aperture measured
# there (reco.RP_APERTURE_MEASURED, tools/fullsim/README.md).
CONFIGS = tuple((name, pu) + reco.RP_APERTURE_MEASURED[name.replace(" ", "")]
                + (col,)
                for name, pu, col in (("5 x 41", 20.5, C_TRUTH),
                                      ("10 x 100", 50.0, C_FIT),
                                      ("18 x 275", 137.5, C_ALT)))

SIGMA_THETA_HA = reco.SIGMA_THETA_HA        # 72.7 urad, high-acceptance optics
ENVELOPE = 10.0 * SIGMA_THETA_HA            # the 10 sigma retraction
R12 = 30.6                                  # m, IP angle -> pot-plane x (measured)


def coherent_acceptance(theta_x, p_ion, slope_b, theta_y):
    """Tagged fraction of an exp(-B t) recoil outside a rectangular
    cutout of half-widths (theta_x, theta_y) in ANGLE."""
    return np.array([reco.rp_hole_acceptance(slope_b, t * p_ion,
                                             theta_y * p_ion)["acc"]
                     for t in np.atleast_1d(theta_x)])


def alpha_acceptance(theta, p_per_nucleon, n=200000, beta=0.30, seed=7):
    """Fraction of 6Li alpha spectators above a near-beam ANGLE cut."""
    k = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG, p_per_nucleon, n,
                                    beta=beta, rng=np.random.default_rng(seed))
    return np.array([float((k["theta"] > t).mean())
                     for t in np.atleast_1d(theta)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slope-b", type=float, default=50.0,
                    help="coherent t-slope [GeV^-2]")
    ap.add_argument("--n-spectator", type=int, default=200000)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    th = np.linspace(0.15e-3, 2.2e-3, 300)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.4))

    lines = ["# near-beam aperture scan (plans/09)",
             "# 10 sigma envelope = %.3f mrad; R12 = %.1f m, so 1 mrad = %.1f mm"
             % (1e3 * ENVELOPE, R12, R12),
             "# the vertical half-width is held at the measured value per "
             "optics (3.0 / 3.0 / 2.3 mrad)"]
    for name, pu, meas, meas_y, col in CONFIGS:
        p_ion = 6.0 * pu
        acc = coherent_acceptance(th, p_ion, args.slope_b, meas_y)
        ax1.plot(1e3 * th, acc, "-", color=col, lw=1.6, label="%s (6Li %g GeV/u)"
                 % (name, pu))
        a_meas = coherent_acceptance(meas, p_ion, args.slope_b, meas_y)[0]
        a_env = coherent_acceptance(ENVELOPE, p_ion, args.slope_b, meas_y)[0]
        ax1.plot([1e3 * meas], [a_meas], "o", color=col, ms=6, mfc="white", mew=1.6)
        ax1.plot([1e3 * ENVELOPE], [a_env], "s", color=col, ms=5)
        lines.append("coherent %-9s aperture %.2f mrad -> %.3e ; envelope %.3f mrad "
                     "-> %.3e ; gain x%.3g"
                     % (name, 1e3 * meas, a_meas, 1e3 * ENVELOPE, a_env,
                        a_env / max(a_meas, 1e-300)))

        acc_a = alpha_acceptance(th, pu, n=args.n_spectator)
        ax2.plot(1e3 * th, acc_a, "-", color=col, lw=1.6, label="%s" % name)
        b_meas = alpha_acceptance(meas, pu, n=args.n_spectator)[0]
        b_env = alpha_acceptance(ENVELOPE, pu, n=args.n_spectator)[0]
        ax2.plot([1e3 * meas], [b_meas], "o", color=col, ms=6, mfc="white", mew=1.6)
        ax2.plot([1e3 * ENVELOPE], [b_env], "s", color=col, ms=5)
        lines.append("alpha    %-9s aperture %.2f mrad -> %.4f ; envelope -> %.4f "
                     "; gain x%.3g" % (name, 1e3 * meas, b_meas, b_env,
                                       b_env / max(b_meas, 1e-12)))

    for ax, ylab, title in (
            (ax1, "tagged fraction of the coherent recoil",
             "(a) coherent intact ⁶Li, exp(−B|t|), B = %g GeV$^{-2}$" % args.slope_b),
            (ax2, "⁶Li α-tag acceptance",
             "(b) ⁶Li α spectator, cluster model (β = 0.30)")):
        ax.axvline(1e3 * ENVELOPE, color=C_GREY, lw=1.1, ls="--")
        ax.annotate("10σ envelope", xy=(1e3 * ENVELOPE, 1.0),
                    xycoords=("data", "axes fraction"), xytext=(3, -9),
                    textcoords="offset points", fontsize=7.5, color=C_GREY,
                    ha="left", va="top")
        ax.set_yscale("log")
        ax.set_xlabel("near-beam aperture, half-width in angle [mrad]")
        ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=9.5)
        ax.grid(alpha=0.25, lw=0.5)
        ax.legend(fontsize=7.5, loc="lower left")
    ax1.set_ylim(1e-18, 2.0)
    ax2.set_ylim(1e-4, 2.0)
    fig.suptitle("What the near-beam aperture is worth.  Circles: the ePIC silicon "
                 "aperture measured per optics in the Sep-2024 geometry\n(tools/fullsim; "
                 "the pot layout has since changed — see plans/09).  Squares: the 10σ beam "
                 "envelope, the machine limit.", fontsize=9)
    fig.tight_layout()
    out = outdir / "nearbeam_aperture_6Li.png"
    fig.savefig(out, dpi=140)
    lines.append("wrote %s" % out)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
