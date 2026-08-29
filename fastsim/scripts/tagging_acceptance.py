#!/usr/bin/env python3
"""Spectator-tagging acceptance for e + 6,7Li at IP6, per beam configuration
and optics (plans/02 step 1.5.4; re-derived 2026-08-28 on the Yellow Report
divergences of plans/10).

Cluster-model spectators (polli_fastsim.spectator) folded with the
far-forward windows (polli_fastsim.farforward).  The headline physics
question: does the soft alpha from 6Li survive the Roman-Pot near-beam
envelope, and how cleanly does the off-rigidity alpha from 7Li land in the
Roman Pots?

The near-beam envelope is an ANGLE, a rectangle of half-widths
10 (sigma_h, sigma_v) at the IP that a fragment at beam rigidity must clear
in one of the two planes.  Until 2026-08-28 the fast simulation applied one
energy-independent, isotropic, proton-derived divergence (73 / 164 microrad)
at every configuration; the Yellow Report's own tables give 220/380,
180/180 and 65/65 microrad for the proton at 5 x 41 / 10 x 100 / 18 x 275,
which a gamma-matched ion inherits and a rigidity-capped one scales by
sqrt(beta*gamma_p / beta*gamma_ion) (farforward.sigma_theta_for).  This
script evaluates the tag at

  * the Yellow Report high-acceptance and high-divergence optics of each
    configuration (farforward.yr_optics);
  * the lithium TAGGING optics of Report 1 Section 6.1 -- the horizontal
    beta* raised to the optimum of acceptance x luminosity for the coherent
    recoil, pots following the envelope (farforward.tagging_optics), whose
    luminosity is 1/7 - 1/13 of the high-acceptance value;
  * the legacy proton-derived optics, for reproduction of the numbers
    published before 2026-08-28.

Outputs: acceptance tables (text) and a beta-scan of the wave-function
tail, plus pT / R / theta distributions per channel.
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams  # noqa: E402
from polli_fastsim import farforward as ff  # noqa: E402
from polli_fastsim import spectator as sp  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def channel_plots(channel, kin, outdir, tag):
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.6))
    axes[0].hist(kin["pT"], bins=120, range=(0, 0.6), histtype="step",
                 density=True)
    axes[0].set_xlabel(r"$p_T$ [GeV]")
    axes[1].hist(kin["R"], bins=120, histtype="step", density=True)
    for lo, hi, c in ((*ff.RP_R_WINDOW, "seagreen"),
                      (*ff.OMD_R_WINDOW, "orange")):
        axes[1].axvspan(lo, hi, alpha=0.15, color=c)
    axes[1].set_xlabel(r"rigidity ratio $R$")
    axes[2].hist(kin["theta"] * 1e3, bins=120, histtype="step", density=True)
    axes[2].axvline(5.0, color="gray", ls="--", lw=1)
    axes[2].set_xlabel(r"$\theta$ [mrad]")
    fig.suptitle(channel.name, fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f"spectator_{tag}.png", dpi=140)
    plt.close(fig)


def optics_menu(config):
    """The optics at which the tag is evaluated for one configuration."""
    return (("YR high-acceptance", ff.yr_optics(config, "high-acceptance")),
            ("YR high-divergence", ff.yr_optics(config, "high-divergence")),
            ("tagging optics", ff.tagging_optics(config)),
            ("legacy 73 urad", ff.HIGH_ACCEPTANCE))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nevents", type=int, default=400_000)
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    lines = ["# Spectator-tagging acceptance (cluster model, far-forward "
             "windows of plans/03 SS2.2), per configuration and optics "
             "(2026-08-28: Yellow Report divergences, plans/10)",
             "# the near-beam envelope is a rectangle 10 (sigma_h, sigma_v) "
             "in angle; 'tagged' = any far-forward system; "
             "'near-beam' = the R ~ 1 pT tail alone",
             "# beta = wave-function short-range scale [GeV] (tail "
             "uncertainty scan; 0.30 central)", ""]
    rng = np.random.default_rng(7)

    for channel in sp.CHANNELS:
        isotope = "6Li" if channel.beam_A == 6 else "7Li"
        for config in beams.default_configs(isotope):
            p_u = config.ion_momentum_per_nucleon
            lines.append(f"== {channel.name}  ({config.label()}, "
                         f"kappa = {channel.kappa*1e3:.1f} MeV, "
                         f"L = {channel.l_wave}) ==")
            for beta in (0.20, 0.30, 0.40):
                kin = sp.spectator_lab_kinematics(channel, p_u, args.nevents,
                                                  beta=beta, rng=rng)
                for name, optics in optics_menu(config):
                    acc = ff.acceptance_summary(kin["R"], kin["theta"],
                                                kin["pT"], optics,
                                                phi=kin["phi"],
                                                pot_config=ff.yr_config_key(config))
                    tagged = 1.0 - acc["lost"]
                    env = optics.envelope
                    parts = "  ".join(f"{k}:{v:6.3f}" for k, v in acc.items()
                                      if v > 5e-4 and k != "lost")
                    lines.append(
                        f"  beta={beta:.2f}  {name:20s} "
                        f"[{1e3*env[0]:.2f} x {1e3*env[1]:.2f} mrad, "
                        f"L/L_HA = {optics.lumi_fraction:.3f}]  "
                        f"tagged={tagged:6.4f}   {parts}")
                if beta == 0.30 and config is beams.default_configs(isotope)[2]:
                    tag = (channel.spectator + "_" + isotope).replace("/", "")
                    channel_plots(channel, kin, outdir, tag)
            lines.append("")

    text = "\n".join(lines)
    (outdir / "tagging_acceptance.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
