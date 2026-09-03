#!/usr/bin/env python3
"""The two figures of report 3 (the machine and the detector).

(a) THE ION ENERGY MENU.  The HSR and the ESR must have equal revolution
    period and the electrons are ultrarelativistic, so the hadron GAMMA is
    fixed by the ring circumference and the magnets then supply whatever
    rigidity that gamma demands -- up to the 275 GeV proton rigidity cap.
    Ions are therefore GAMMA-MATCHED to the proton configurations, and each
    species' menu is the isolated low-gamma point plus a continuous band
    running up to its own rigidity limit.  For 6Li: "41, and 99-138 GeV/u",
    with nothing between (plans/10).

(b) THE BEAM ANGULAR DIVERGENCE, which the coherent tagged fraction of this
    programme is an exponential in the square of.  Yellow Report Table 10.1
    per configuration and optics for protons, with the species step applied
    only where the rigidity cap binds -- so a gamma-matched 6Li carries the
    proton's divergence exactly, and only the top configuration pays the
    sqrt(2).  The single energy-independent 72.7 microrad this repository
    used before 2026-08-27 is drawn for comparison.

Usage:  python3 scripts/eic_beam_figures.py [--outdir .]
"""

import argparse
import pathlib
import sys

import numpy as np

_S = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_S.parent))
sys.path.insert(0, str(_S.parent.parent / "fastsim"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from polli_fastsim import beams  # noqa: E402
from polligen import reco  # noqa: E402

C = {"p": "#8A8A8A", "d": "#B8860B", "3He": "#1B7F5B",
     "6Li": "#1F4E79", "7Li": "#C0392B"}
SPECIES = ("p", "d", "3He", "6Li", "7Li")


def energy_menu(ax):
    """Per-nucleon momentum against species, showing the accessible bands."""
    gam_lo = 41.0 / beams.PROTON_MASS                 # the isolated point
    gam_hi_bottom = 100.0 / beams.PROTON_MASS         # bottom of the band
    for i, name in enumerate(SPECIES):
        ion = beams.IONS[name]
        top = ion.momentum_per_nucleon_max
        lo = ion.momentum_per_nucleon_at(41.0)
        band_lo = min(ion.momentum_per_nucleon_at(100.0), top)
        ax.plot([i, i], [band_lo, top], "-", color=C[name], lw=7,
                solid_capstyle="butt", alpha=0.85)
        ax.plot([i], [lo], "o", color=C[name], ms=9)
        ax.annotate("%.0f" % lo, (i, lo), xytext=(0, -14),
                    textcoords="offset points", fontsize=7.5, color=C[name],
                    ha="center")
        ax.annotate("%.3g–%.0f" % (band_lo, top), (i, top),
                    xytext=(0, 6), textcoords="offset points", fontsize=7.5,
                    color=C[name], ha="center")
        # what rigidity-scaling would have given, for the two lower points
        for pe in (41.0, 100.0):
            ax.plot([i], [pe * ion.Z / ion.A], "x", color="0.55", ms=6,
                    mew=1.3, zorder=1)
    ax.axhspan(43, 97, color="#F2E2E2", zorder=0)
    ax.annotate("no configuration in this band —\nthe circumference cannot be\n"
                "adjusted far enough",
                xy=(-0.55, 65), fontsize=7.5, color="#8A3A3A",
                ha="left", va="center")
    ax.set_xticks(range(len(SPECIES)))
    ax.set_xticklabels(SPECIES)
    ax.set_xlim(-0.75, len(SPECIES) - 0.15)
    ax.set_ylabel("momentum per nucleon [GeV]")
    ax.set_yscale("log")
    ax.set_ylim(12, 480)
    ax.set_title("the ion energy menu: γ-matched, then rigidity-capped",
                 fontsize=9.5)
    ax.grid(alpha=0.25, lw=0.5, axis="y")
    ax.plot([], [], "o", color="0.3", label="γ-matched to the 41 GeV proton")
    ax.plot([], [], "-", color="0.3", lw=7, alpha=0.85,
            label="continuous band, up to the rigidity cap")
    ax.plot([], [], "x", color="0.55", mew=1.3,
            label="what rigidity scaling would give")
    ax.legend(fontsize=7, loc="lower left", framealpha=0.95)


def divergence(ax):
    labels = ("5 × 41", "10 × 100", "18 × 275")
    x = np.arange(len(labels), dtype=float)
    w = 0.19
    series = (("p", "high-divergence", -1.5, "//", 0.40),
              ("p", "high-acceptance", -0.5, None, 0.90),
              ("6Li", "high-divergence", 0.5, "//", 0.40),
              ("6Li", "high-acceptance", 1.5, None, 0.90))
    for name, opt, off, hatch, alpha in series:
        vals = [1e6 * reco.sigma_theta_for(c, opt)[0]
                for c in beams.default_configs(name)]
        ax.bar(x + off * w, vals, w, color=C[name], alpha=alpha, hatch=hatch,
               edgecolor="white", lw=0.6,
               label="%s, %s" % (name, opt.split("-")[1]))

    # the vertical plane, marked where it differs from the horizontal
    for i, c in enumerate(beams.default_configs("6Li")):
        h, v = reco.sigma_theta_for(c)
        if abs(v - h) > 1e-9:
            ax.plot([x[i] + 1.5 * w], [1e6 * v], "v", color=C["6Li"], ms=9,
                    mfc="white", mew=1.6, zorder=5)
            ax.annotate("vertical\n%.0f μrad" % (1e6 * v),
                        (x[i] + 1.5 * w, 1e6 * v), xytext=(9, 0),
                        textcoords="offset points", fontsize=7.2,
                        color=C["6Li"], va="center")

    ax.axhline(72.7, color="k", ls="--", lw=1.2, zorder=4)
    ax.annotate("72.7 μrad — the single energy-independent value used\n"
                "throughout this programme before 2026-08-27",
                xy=(-0.42, 72.7), xytext=(0, 6), textcoords="offset points",
                fontsize=7.2, ha="left")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlim(-0.55, len(labels) - 0.25)
    ax.set_ylim(0, 470)
    ax.set_ylabel(r"beam angular divergence $\sigma_\theta$, horizontal [μrad]")
    ax.set_title("the divergence the coherent acceptance is "
                 "exp(−B(10σ·A·p)²) in", fontsize=9.5)
    ax.grid(alpha=0.25, lw=0.5, axis="y")
    ax.legend(fontsize=7, ncol=2, loc="upper center")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    out = pathlib.Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    for fn, draw, size in (("eic_energy_menu.png", energy_menu, (6.8, 4.4)),
                           ("eic_beam_divergence.png", divergence, (6.8, 4.4))):
        fig, ax = plt.subplots(figsize=size)
        draw(ax)
        fig.tight_layout()
        fig.savefig(out / fn, dpi=140)
        plt.close(fig)
        print("wrote %s" % (out / fn))

    print("\n6Li menu: %.1f, and %.1f-%.1f GeV/u"
          % (beams.LI6.momentum_per_nucleon_at(41.0),
             beams.LI6.momentum_per_nucleon_at(100.0),
             beams.LI6.momentum_per_nucleon_max))
    for name in SPECIES:
        c = beams.default_configs(name)
        print("  %-5s %s" % (name, " / ".join("%.1f" % x.ion_momentum_per_nucleon
                                              for x in c)))


if __name__ == "__main__":
    main()
