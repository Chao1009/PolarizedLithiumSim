#!/usr/bin/env python3
"""Money plot 1 (plan step 1.3.1): polarized EMC effect on 7Li.

Projected statistical precision on the ratio
  DR(x) = (g1A/F1A) / A1_naive,   A1_naive = [P_p g1p + P_n g1n]/F1A
combined over Q2 bins and the three energy settings, overlaid on the two
camps of the polarized-EMC literature. Since 2026-08-28 those are the
DIGITIZED published curves -- Cloet-Bentz-Thomas PLB 642:210 Fig. 6 (7Li at
Q2 = 5 GeV2, their Eq. 23) and Tronchin-Matevosyan-Thomas PLB 783:247 Fig. 4
(nuclear matter at Q2 = 10 GeV2, rescaled to 7Li strength) -- not the
constants 2 and 1 applied to a hand-written EMC table, which `--emc-mode
constant` restores. The discrimination significance
|DR_CBT - DR_TMT| / dDR is the FOM.

The two luminosities drawn (10 and 100 fb^-1/u) are PROGRAMME luminosities;
`--run-share` is the fraction of them this observable is given (plans/07
WP2), so dDR scales as 1/sqrt(share) and the published figure is at
share 1.
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.inputs import get_backends
from polli_fastsim.polarized import (CBT_TABLE, POLEMC_VALENCE_WINDOW,
                                     cbt_polarized_emc_ratio, curve_x_range,
                                     tmt_polarized_emc_ratio,
                                     tmt_published_emc_ratio, toy_b1,
                                     toy_delta_gluon)
from polli_fastsim.structure import NuclearF2

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def delta_dr_per_x(ion_name, lumi, backends, min_events=100, run_share=1.0):
    """Combined-over-(Q2, energies) statistical error on DR vs x, at
    `lumi` of programme luminosity of which this observable gets
    `run_share` (so the error scales as 1/sqrt(lumi * run_share))."""
    g1m = backends["g1"]
    inv2_tot = None
    x_centers = None
    for cfg in beams.default_configs(ion_name):
        sc = fom.Scenario(lumi_fb_per_nucleon=lumi, run_share=run_share)
        proj = fom.project_rates(
            cfg, sc, nuclear_f2=NuclearF2(cfg.ion, base=backends["base"]))
        obs = fom.project_observables(cfg, sc, proj, g1m, toy_b1,
                                      toy_delta_gluon)
        x_centers = proj.x[:, 0]
        # A1_naive per bin (no medium modification in the denominator)
        f1 = proj.extras["nf2"].f1a(proj.x, proj.q2) / cfg.ion.A
        g1_naive = g1m.g1_nucleus(cfg.ion, proj.x, proj.q2) / cfg.ion.A
        a1_naive = np.abs(g1_naive / np.maximum(f1, 1e-30))
        err_dr = obs["err_g1_over_f1"] / np.maximum(a1_naive, 1e-12)
        use = proj.accepted & (proj.n_events >= min_events)
        inv2 = np.where(use & (err_dr > 0), 1.0 / err_dr**2, 0.0)
        inv2_tot = inv2.sum(axis=1) if inv2_tot is None \
            else inv2_tot + inv2.sum(axis=1)
    err = np.full(x_centers.shape, np.inf)
    np.divide(1.0, np.sqrt(inv2_tot), out=err, where=inv2_tot > 0)
    return x_centers, err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="7Li", choices=["6Li", "7Li"])
    ap.add_argument("--pdf", default="toy", choices=["toy", "grid"])
    ap.add_argument("--emc-mode", default="digitized", dest="emc_mode",
                    choices=["digitized", "constant"],
                    help="polarized-EMC curves: 'digitized' (default) reads "
                         "the published figures, 'constant' the "
                         "pre-2026-08-28 2x / 1x stand-ins")
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--run-share", type=float, default=1.0, dest="run_share",
                    help="this observable's share of the programme "
                         "luminosity (plans/07 WP2; default 1.0, which the "
                         "published figure assumes)")
    args = ap.parse_args()
    if not args.run_share > 0:
        ap.error("--run-share must be positive")
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    backends = get_backends(args.pdf)
    def cbt(z):
        return cbt_polarized_emc_ratio(z, mode=args.emc_mode)

    def tmt(z):
        return tmt_polarized_emc_ratio(z, mode=args.emc_mode)

    def tmt_raw(z):
        """TMT's published nuclear-matter curve, before the 7Li transfer.

        The transfer factor is a single constant fitted in the valence
        window, so applying it two decades below that window manufactures
        a separation the two papers do not have: the published curves
        agree to better than 0.007 for x < 0.3.  Both separations are
        printed so the reader can tell which is which."""
        if args.emc_mode == "constant":
            return tmt(z)
        return tmt_published_emc_ratio(z)

    print("# money_polemc.py  ion=%s  pdf=%s  emc=%s  %s"
          % (args.ion, args.pdf, args.emc_mode,
             fom.run_share_header(10.0, args.run_share)))

    xs = np.logspace(np.log10(0.005), np.log10(0.9), 300)
    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(7, 7), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1]})
    tag = ("digitized" if args.emc_mode == "digitized" else "scenario")
    ax.plot(xs, cbt(xs), "crimson",
            label=f"CBT (PLB 642:210 Fig. 6, Eq. 23) [{tag}]")
    ax.plot(xs, tmt(xs), "navy", ls="--",
            label=f"TMT (PLB 783:247 Fig. 4 $\\to$ $^7$Li) [{tag}]")
    if args.emc_mode == "digitized":
        # The same curve BEFORE the valence-window transfer.  It is drawn
        # because the transferred pair is otherwise read as a low-x
        # discrimination the two papers do not have: they agree to better
        # than 0.008 for x < 0.3 and part only inside the shaded window.
        ax.plot(xs, tmt_raw(xs), "navy", ls=":", lw=1,
                label="TMT as published (nuclear matter, no transfer)")
    ax.axhline(1.0, color="gray", lw=0.5)
    vlo, vhi = POLEMC_VALENCE_WINDOW
    for a in (ax, ax2):
        a.axvspan(vlo, vhi, color="0.85", alpha=0.5, lw=0, zorder=0)
    ax.text(np.sqrt(vlo * vhi), 1.26, "valence window\n(transfer defined)",
            ha="center", va="top", fontsize=7, color="0.35")

    for lumi, color, dx in ((10.0, "black", 1.0), (100.0, "seagreen", 1.04)):
        x_c, err = delta_dr_per_x(args.ion, lumi, backends,
                                  run_share=args.run_share)
        ok = np.isfinite(err) & (x_c > 0.004) & (x_c < 0.9) & (err < 0.5)
        mid = 0.5 * (cbt(x_c[ok])
                     + tmt(x_c[ok]))
        ax.errorbar(x_c[ok] * dx, mid, yerr=err[ok], fmt=".", ms=4,
                    color=color, lw=1, capsize=0,
                    label=f"proj. $\\delta\\Delta R$, "
                          f"{lumi * args.run_share:g} fb$^{{-1}}$/u")
        sep_c = np.abs(cbt(x_c[ok])
                       - tmt(x_c[ok]))
        ax2.plot(x_c[ok], sep_c / err[ok], ".-", color=color, lw=1, ms=4,
                 label=f"{lumi * args.run_share:g} fb$^{{-1}}$/u")

    ax.set_xscale("log")
    ax.set_ylim(0.5, 1.3)
    ax.set_ylabel(r"$\Delta R(x) = g_1^A/(P_p g_1^p + P_n g_1^n)$")
    ax.set_title(
        f"Polarized EMC effect, {args.ion}: CBT vs TMT "
        f"({args.emc_mode})\n"
        f"(3 energy settings combined; $P_e$=0.7, $P_z$=0.7; "
        f"stat. only; {backends['tag'].upper()} inputs"
        + ("" if args.run_share == 1.0
           else f"; run share {args.run_share:g}") + ")", fontsize=10)
    ax.legend(fontsize=8, loc="lower left")
    ax2.axhline(5, color="gray", ls=":", lw=1)
    ax2.set_yscale("log")
    ax2.set_xlabel(r"$x$")
    ax2.set_ylabel(r"$|\Delta R_{CBT}-\Delta R_{TMT}|/\delta$")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    share_key = fom.run_share_tag(args.run_share)
    stem = f"money_polemc_{args.ion}_{backends['tag']}"
    if args.emc_mode != "digitized":    # never overwrite the published PNG
        stem = f"{stem}_{args.emc_mode}"
    if share_key:
        stem = f"{stem}_{share_key}"
    path = outdir / f"{stem}.png"
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")
    x_c, err10 = delta_dr_per_x(args.ion, 10.0, backends,
                                run_share=args.run_share)
    for xq in (0.1, 0.3, 0.5, 0.7):
        i = np.argmin(np.abs(x_c - xq))
        sep_i = abs(cbt(x_c[i]) - tmt(x_c[i]))
        raw_i = abs(cbt(x_c[i]) - tmt_raw(x_c[i]))
        print(f"  x={x_c[i]:.2f}: dDR({10.0 * args.run_share:g}/fb)="
              f"{err10[i]:.4f}  CBT-TMT sep={sep_i:.3f}  "
              f"({sep_i / err10[i]:.2f} sigma here, "
              f"{sep_i / err10[i] * np.sqrt(10.0):.2f} at 100 fb^-1/u)"
              f"  [published curves differ by {raw_i:.3f} here]")
    good = np.isfinite(err10) & (x_c > 0.004) & (x_c < 0.9) & (err10 < 0.5)
    signif = np.abs(cbt(x_c[good]) - tmt(x_c[good])) / err10[good]
    j = int(np.argmax(signif))
    print(f"  best bin x={x_c[good][j]:.3f}: {signif[j]:.2f} sigma at "
          f"{10.0 * args.run_share:g} fb^-1/u, "
          f"{signif[j] * np.sqrt(10.0):.2f} at 100 fb^-1/u; "
          f"{(signif > 1.0).sum()} of {good.sum()} bins above 1 sigma")
    # The transfer factor is defined in the valence window and constant
    # outside it, so the honest headline is the reach where the two
    # PUBLISHED curves actually differ.  Both are printed; Report 0
    # section 5.3 and the manual quote both.
    vlo, vhi = POLEMC_VALENCE_WINDOW
    win = good & (x_c > vlo) & (x_c < vhi)
    if win.any():
        s_win = np.abs(cbt(x_c[win]) - tmt(x_c[win])) / err10[win]
        k = int(np.argmax(s_win))
        print(f"  valence window {vlo:g}<x<{vhi:g}: best bin "
              f"x={x_c[win][k]:.3f}, {s_win[k]:.2f} sigma at "
              f"{10.0 * args.run_share:g} fb^-1/u, "
              f"{s_win[k] * np.sqrt(10.0):.2f} at 100 fb^-1/u; "
              f"{(s_win > 1.0).sum()} of {win.sum()} bins above 1 sigma")
    # Restricted to where BOTH tables have data: below CBT's x_min its
    # curve is frozen at its endpoint, so a difference there measures the
    # extrapolation and not the papers.
    x_lo_ok = curve_x_range(CBT_TABLE)[0]
    lo = good & (x_c > x_lo_ok) & (x_c < 0.3)
    hi = good & (x_c > vlo)
    print(f"  published-curve separation: max "
          f"{np.abs(cbt(x_c[lo]) - tmt_raw(x_c[lo])).max():.4f} over "
          f"{x_lo_ok:.3f} < x < 0.3, where the transferred pair shows "
          f"{np.abs(cbt(x_c[lo]) - tmt(x_c[lo])).min():.4f}-"
          f"{np.abs(cbt(x_c[lo]) - tmt(x_c[lo])).max():.4f}; max "
          f"{np.abs(cbt(x_c[hi]) - tmt_raw(x_c[hi])).max():.4f} above "
          f"x = {vlo:g}")


if __name__ == "__main__":
    main()
