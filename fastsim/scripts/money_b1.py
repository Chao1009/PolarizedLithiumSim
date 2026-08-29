#!/usr/bin/env python3
"""Money plot 2 (plan step 1.3.2): tensor asymmetry A_zz / b1 on 6Li.

Projected per-x statistical precision on A_zz (longitudinal tensor,
equal-thirds spin states) against the two camps of the b1 literature,
both now digitized: Miller's pion + hidden-colour total (PRC 89:045203
Fig. 5), which reproduces HERMES, and the standard convolution of
Cosyn-Dong-Kumano-Sargsian (PRD 95:074036 Fig. 4), which is an order of
magnitude below it at x >~ 0.2. Shown for two tensor polarizations --
this plot sets the P_zz requirement for the source.

`--transfer rank2` (the default) carries the deuteron b1 into 6Li with the
rank-2 tensor transfer of the tagged two-cluster model, 0.921947, times the
2-of-6 per-nucleon factor; `--transfer legacy` restores the 0.87 x 1 of the
figures published before 2026-08-28, which was the VECTOR dilution and no
per-nucleon factor at all -- signal and error on different normalisations
(plans/08 D9).

`--lumi` is the programme luminosity and `--run-share` the fraction of it
this observable is given (plans/07 WP2); the errors scale as
1/sqrt(lumi x share) and every published curve is at share 1.
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.asymmetries import azz as azz_formula
from polli_fastsim.inputs import get_backends
from polli_fastsim.polarized import (LI6_B1_LEGACY_TRANSFER,
                                     LI6_B1_PER_NUCLEON,
                                     LI6_B1_RANK2_TRANSFER, b1_convolution,
                                     b1_li6_from_deuteron, toy_b1,
                                     toy_delta_gluon)
from polli_fastsim.structure import NuclearF2

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lumi", type=float, default=10.0,
                    help="programme luminosity [fb^-1/nucleon] (default: "
                         "%(default)s, one EIC year)")
    ap.add_argument("--run-share", type=float, default=1.0, dest="run_share",
                    help="this observable's share of that programme "
                         "luminosity (plans/07 WP2; default 1.0, which the "
                         "published figure assumes)")
    ap.add_argument("--pdf", default="toy", choices=["toy", "grid"])
    ap.add_argument("--transfer", default="rank2",
                    choices=["rank2", "legacy"],
                    help="deuteron -> 6Li b1 transfer: 'rank2' (default) is "
                         "0.921947 x 2/6, 'legacy' the pre-2026-08-28 "
                         "0.87 x 1 (plans/08 D9)")
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    if not args.run_share > 0:
        ap.error("--run-share must be positive")
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    lumi_eff = args.lumi * args.run_share
    transfer, per_nucleon = ((LI6_B1_RANK2_TRANSFER, LI6_B1_PER_NUCLEON)
                             if args.transfer == "rank2"
                             else (LI6_B1_LEGACY_TRANSFER, 1.0))
    factor = transfer * per_nucleon
    print("# money_b1.py  pdf=%s  %s"
          % (args.pdf, fom.run_share_header(args.lumi, args.run_share)))
    print("# b1 transfer: %s, %.6f x %.6f = %.6f"
          % (args.transfer, transfer, per_nucleon, factor))

    backends = get_backends(args.pdf)
    g1m = backends["g1"]
    fig, ax = plt.subplots(figsize=(7, 5))

    # scenario Azz(x) curves at a representative Q2 slice (combined plot
    # uses the bin-by-bin values; curves drawn at Q2 = 4 GeV2)
    xs = np.logspace(np.log10(0.003), np.log10(0.7), 200)
    q2s = np.full_like(xs, 4.0)
    cfg0 = beams.default_configs("6Li")[0]
    s0 = cfg0.sqrt_s_per_nucleon**2
    y0 = np.clip(q2s / (s0 * xs), 1e-4, 1.0)
    nf2 = NuclearF2(cfg0.ion, base=backends["base"])
    f1 = nf2.f1a(xs, q2s) / cfg0.ion.A
    f2 = nf2.f2a(xs, q2s) / cfg0.ion.A
    signals = {}
    for b1f, color, label in (
            (toy_b1, "crimson", "Miller $b_1$ (HERMES-like) [digitized]"),
            (b1_convolution, "navy", "CDKS convolution $b_1$ [digitized]")):
        b1 = b1_li6_from_deuteron(b1f(xs, q2s, f1), transfer, per_nucleon)
        signals[label] = np.abs(azz_formula(b1, f1, f2, xs, y0))
        ax.plot(xs, signals[label], color=color, label=label)

    # projected per-x errors, combined over Q2 and the 3 energy settings
    errors = {}
    for pzz, fmt in ((0.60, "ko"), (0.80, "g^")):
        inv2_tot, x_c = None, None
        for cfg in beams.default_configs("6Li"):
            sc = fom.Scenario(lumi_fb_per_nucleon=args.lumi,
                              run_share=args.run_share,
                              pol_ion_tensor=pzz)
            proj = fom.project_rates(
                cfg, sc, nuclear_f2=NuclearF2(cfg.ion, base=backends["base"]))
            obs = fom.project_observables(cfg, sc, proj, g1m, toy_b1,
                                          toy_delta_gluon)
            x_c = proj.x[:, 0]
            use = proj.accepted & (proj.n_events >= 100)
            inv2 = np.where(use, 1.0 / obs["err_azz"]**2, 0.0).sum(axis=1)
            inv2_tot = inv2 if inv2_tot is None else inv2_tot + inv2
        err = np.full(x_c.shape, np.inf)
        np.divide(1.0, np.sqrt(inv2_tot), out=err, where=inv2_tot > 0)
        errors[pzz] = (x_c, err)
        ok = np.isfinite(err) & (x_c > 0.003) & (x_c < 0.7)
        ax.plot(x_c[ok], err[ok], fmt, ms=4, ls="-", lw=0.8,
                label=f"$\\delta A_{{zz}}$/x-bin, $P_{{zz}}$={pzz:g}, "
                      f"{lumi_eff:g} fb$^{{-1}}$/u")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$|A_{zz}|$  and  $\delta A_{zz}$")
    ax.set_ylim(1e-5, 0.3)
    share_note = ("" if args.run_share == 1.0
                  else f"; run share {args.run_share:g} of the programme year")
    ax.set_title("Tensor asymmetry on $^6$Li (embedded deuteron, "
                 f"transfer {args.transfer} = {factor:.4f})\n"
                 "3 energies combined; stat. only; "
                 f"{backends['tag'].upper()} inputs{share_note}", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    share_key = fom.run_share_tag(args.run_share)
    stem = f"money_b1_6Li_{backends['tag']}"
    if args.transfer != "rank2":        # never overwrite the published PNG
        stem = f"{stem}_{args.transfer}"
    if share_key:
        stem = f"{stem}_{share_key}"
    path = outdir / f"{stem}.png"
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")
    for label, sig in signals.items():
        row = "  " + label.split(" $")[0]
        for xq in (0.005, 0.01, 0.03, 0.07, 0.2, 0.5):
            row += f"  |Azz|({xq:g})={np.interp(xq, xs, sig):.3g}"
        print(row)
    for pzz, (x_c, err) in errors.items():
        row = f"  dAzz Pzz={pzz:g}"
        for xq in (0.0035, 0.01, 0.03, 0.28, 0.56):
            i = np.argmin(np.abs(x_c - xq))
            row += f"  {x_c[i]:.4f}:{err[i]:.3g}"
        print(row)
    x_c, err = errors[0.60]
    for label, sig in signals.items():
        s_on_bins = np.interp(x_c, xs, sig)
        row = "  significance (Pzz=0.6) " + label.split(" $")[0]
        for xq in (0.0035, 0.01, 0.03, 0.07, 0.2, 0.5):
            i = np.argmin(np.abs(x_c - xq))
            row += f"  {x_c[i]:.4f}:{s_on_bins[i] / err[i]:.1f}sig"
        print(row)


if __name__ == "__main__":
    main()
