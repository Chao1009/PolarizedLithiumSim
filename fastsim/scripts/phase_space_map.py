#!/usr/bin/env python3
"""Produce x-Q2 phase-space / rate / FOM maps for e+6Li and e+7Li.

Usage:  python3 scripts/phase_space_map.py [--ion 7Li] [--lumi 10]
                [--run-share 1.0] [--outdir out]

`--lumi` is the programme luminosity and `--run-share` this observable's
share of it (plans/07 WP2); rates scale as their product and every
published map is at share 1.

Outputs per beam configuration:
  - x-Q2 coverage + event-rate heat map (log axes)
  - per-bin statistical precision maps for A_par (g1), Azz (b1), A_cos2phi (Delta)
  - a text summary table of integrated rates and example-bin precisions
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.polarized import ToyG1, toy_b1, toy_delta_gluon

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def plot_map(proj, values, title, label, path, vmin=None, vmax=None, log=True):
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    vals = np.ma.masked_where(~proj.accepted | (values <= 0), values)
    norm = LogNorm(vmin=vmin, vmax=vmax) if log else None
    pcm = ax.pcolormesh(proj.x_edges, proj.q2_edges, vals.T, norm=norm,
                        cmap="viridis", shading="auto")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]")
    ax.set_title(title, fontsize=10)
    fig.colorbar(pcm, ax=ax, label=label)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="7Li", choices=list(beams.IONS))
    ap.add_argument("--lumi", type=float, default=10.0,
                    help="programme integrated luminosity per nucleon "
                         "[fb^-1] (default: %(default)s, one EIC year)")
    ap.add_argument("--run-share", type=float, default=1.0, dest="run_share",
                    help="this observable's share of that programme "
                         "luminosity (plans/07 WP2; default 1.0, which the "
                         "published maps assume)")
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    if not args.run_share > 0:
        ap.error("--run-share must be positive")

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi,
                            run_share=args.run_share)
    lumi_eff = scenario.lumi_effective_fb_per_nucleon
    # a non-default share writes its own files: the published maps are the
    # share-1 ones
    share_key = fom.run_share_tag(args.run_share)
    suffix = ("_" + share_key) if share_key else ""
    g1_model = ToyG1()

    summary = []
    for cfg in beams.default_configs(args.ion):
        proj = fom.project_rates(cfg, scenario)
        obs = fom.project_observables(cfg, scenario, proj, g1_model,
                                      toy_b1, toy_delta_gluon)
        tag = cfg.label().replace(" ", "").replace("(", "").replace(")", "") \
                         .replace("/u", "u").replace("x", "_x_")
        plot_map(proj, proj.n_events,
                 f"{cfg.label()},  L={lumi_eff:g} fb$^{{-1}}$/u  "
                 f"$\\sqrt{{s_{{eN}}}}$={cfg.sqrt_s_per_nucleon:.0f} GeV",
                 "events / bin", outdir / f"rate_{tag}{suffix}.png",
                 vmin=1, vmax=None)
        plot_map(proj, obs["err_g1_over_f1"],
                 f"{cfg.label()}: stat. precision on $g_1^A/F_1^A$",
                 r"$\delta(g_1/F_1)$", outdir / f"err_g1f1_{tag}{suffix}.png",
                 vmin=1e-4, vmax=1)
        plot_map(proj, obs["err_azz"],
                 f"{cfg.label()}: stat. precision on $A_{{zz}}$",
                 r"$\delta A_{zz}$", outdir / f"err_azz_{tag}{suffix}.png",
                 vmin=1e-4, vmax=1)
        plot_map(proj, obs["err_a_cos2phi"],
                 f"{cfg.label()}: stat. precision on $A_{{\\cos 2\\phi}}$",
                 r"$\delta A_{\cos 2\phi}$",
                 outdir / f"err_cos2phi_{tag}{suffix}.png",
                 vmin=1e-4, vmax=1)

        n_tot = proj.n_events.sum()
        # example precision at x ~ 0.3, in the highest-statistics accepted
        # Q2 bin at that x (Q2=10 GeV2 falls below y_min at top energy)
        ix = np.argmin(np.abs(proj.x[:, 0] - 0.3))
        row_n = np.where(proj.accepted[ix], proj.n_events[ix], 0.0)
        if row_n.max() > 0:
            iq = int(np.argmax(row_n))
            q2_at = proj.q2[ix, iq]
            example = (
                f"d(g1/F1)|x=0.3,Q2={q2_at:5.1f}={obs['err_g1_over_f1'][ix, iq]:8.2e}  "
                f"dAzz={obs['err_azz'][ix, iq]:8.2e}  "
                f"dA2phi={obs['err_a_cos2phi'][ix, iq]:8.2e}"
            )
        else:
            example = "x=0.3 outside acceptance"
        summary.append(
            f"{cfg.label():26s} sqrt(s_eN)={cfg.sqrt_s_per_nucleon:5.1f} GeV  "
            f"N_DIS={n_tot:9.3e}  " + example
        )

    text = "\n".join(summary)
    (outdir / f"summary_{args.ion}{suffix}.txt").write_text(text + "\n")
    print(f"# ion={args.ion}  L={lumi_eff:g} fb^-1/nucleon  "
          f"Pe={scenario.pol_electron} Pz={scenario.pol_ion_vector} "
          f"Pzz={scenario.pol_ion_tensor}")
    print("# " + fom.run_share_header(args.lumi, args.run_share))
    print(text)
    print(f"\nplots written to {outdir}/")


if __name__ == "__main__":
    main()
