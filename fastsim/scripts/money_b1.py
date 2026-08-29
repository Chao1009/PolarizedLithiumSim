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

`--signal-q2 binned` (the default since 2026-08-29) evaluates the signal
BIN BY BIN, at each x bin's rate-weighted Q2 and y over the accepted
(Q2, energy) cells -- the same cells and the same weights the error is
combined over -- with Miller's b1 carried to that Q2 by his own Fig.-6 Q2
set (`polarized.miller_b1_q2_scale`, interpolated in ln Q2). `--signal-q2
fixed` restores the single Q2 = 4 GeV2 slice at the leading
configuration's y, which is what the figures published before that date
drew.

What moves between the two is NOT that Q2 lever. The accepted bins sit at
<Q2> = 3.23 to 52 GeV2, at or above the top node of Miller's Fig.-6 set
(3.25), so `miller_b1_q2_scale` returns 1.0000 in every bin this
projection quotes and 0.9993 in the lowest -- no accepted bin reaches the
region where his b1 is smaller than its Fig.-5 value. What moves is the
bin's rate-weighted y, and its <Q2> inside F1 and F2: the fixed slice
takes y = Q2/(s x) clipped at 1, which at x = 0.0035 IS the clip, against
the bin's own weighted y = 0.234. The binned signal is smaller by 9.5% at
the lowest bin, by 3-5% through the middle and by 0.3% at x = 0.45 --
largest where the measurement lives and negligible at the top. The lever
is kept as machinery, for a Q2 set that does not top out; it is not what
the numbers are.

The quoted |A_zz| row is evaluated on the dense x grid at the bins'
interpolated (Q2, y), NOT by interpolating the ~39 coarse bin values:
|A_zz| passes through Miller's zero crossing at x = 0.577, so the
bin-to-bin interpolation moved the headline at x = 0.5 by 24% (2.55e-3
against the 3.33e-3 the same kinematics give). The plotted markers and
the per-bin significance row are on the bin centres, where they are exact.

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
                                     MILLER_B1_Q2_REF,
                                     LI6_B1_RANK2_TRANSFER, b1_convolution,
                                     b1_li6_from_deuteron, toy_b1,
                                     toy_delta_gluon)
from polli_fastsim.structure import NuclearF2

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def azz_signal(b1f, x, q2, y, nf2, ion, transfer, per_nucleon,
               q2_evolve=False):
    """|A_zz| of one b1 camp on a grid of (x, Q2, y).

    The one place the signal is formed, so that the binned and the
    fixed-slice paths of `main` cannot drift apart.  At Q2 = 4 the two
    agree bit for bit at equal y: the only Q2 dependence `q2_evolve` adds
    is Miller's Fig.-6 lever, and above his top node (Q2 = 3.25) that is
    exactly 1, while `b1_convolution` carries no lever at all
    (`tests/test_digitized_curves.py`)."""
    f1 = nf2.f1a(x, q2) / ion.A
    f2 = nf2.f2a(x, q2) / ion.A
    kw = {"q2_evolve": True} if (q2_evolve and b1f is toy_b1) else {}
    b1 = b1_li6_from_deuteron(b1f(x, q2, f1, **kw), transfer, per_nucleon)
    return np.abs(azz_formula(b1, f1, f2, x, y))


def bin_kinematics(lumi, run_share, base, min_events=100):
    """Rate-weighted <Q2> and <y> per x bin over the accepted cells of the
    three 6Li configurations, and the mask of bins that have any.

    `asymmetries.err_azz` is sqrt(2/N)/Pzz, so the inverse-variance weight
    of a cell is its event count: this is the same weighting the combined
    error uses, and it is independent of P_zz, so one set of bin
    kinematics serves both tensor polarizations."""
    w_tot = wq = wy = x_c = None
    for cfg in beams.default_configs("6Li"):
        sc = fom.Scenario(lumi_fb_per_nucleon=lumi, run_share=run_share)
        proj = fom.project_rates(cfg, sc,
                                 nuclear_f2=NuclearF2(cfg.ion, base=base))
        x_c = proj.x[:, 0]
        w = np.where(proj.accepted & (proj.n_events >= min_events),
                     proj.n_events, 0.0)
        terms = (w.sum(axis=1), (w * proj.q2).sum(axis=1),
                 (w * proj.extras["y"]).sum(axis=1))
        w_tot, wq, wy = terms if w_tot is None else \
            tuple(a + b for a, b in zip((w_tot, wq, wy), terms))
    ok = w_tot > 0
    q2b = np.where(ok, wq / np.where(ok, w_tot, 1.0), 4.0)
    yb = np.where(ok, wy / np.where(ok, w_tot, 1.0), 1e-3)
    return x_c, q2b, np.clip(yb, 1e-4, 1.0), ok


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
    ap.add_argument("--signal-q2", default="binned", dest="signal_q2",
                    choices=["binned", "fixed"],
                    help="where the signal curves are evaluated: 'binned' "
                         "(default) at each x bin's rate-weighted Q2 and y, "
                         "with Miller's Fig.-6 Q2 set carrying b1 there; "
                         "'fixed' at the pre-2026-08-29 Q2 = 4 slice")
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
    print("# b1 transfer: %s, %.6f x %.6f = %.6f;  signal Q2: %s"
          % (args.transfer, transfer, per_nucleon, factor, args.signal_q2))

    backends = get_backends(args.pdf)
    g1m = backends["g1"]
    fig, ax = plt.subplots(figsize=(7, 5))

    cfg0 = beams.default_configs("6Li")[0]
    nf2 = NuclearF2(cfg0.ion, base=backends["base"])

    # scenario Azz(x): bin by bin at each bin's rate-weighted (Q2, y), or
    # on the pre-2026-08-29 Q2 = 4 slice at the leading configuration's y
    xs = np.logspace(np.log10(0.003), np.log10(0.7), 200)
    q2s = np.full_like(xs, 4.0)
    s0 = cfg0.sqrt_s_per_nucleon**2
    y0 = np.clip(q2s / (s0 * xs), 1e-4, 1.0)
    binned = args.signal_q2 == "binned"
    if binned:
        xb, q2b, yb, okb = bin_kinematics(args.lumi, args.run_share,
                                          backends["base"])
        # The bin centres are ~39 coarse log-spaced points and |A_zz|
        # oscillates through Miller's zero crossing at x = 0.577, so
        # interpolating the BIN VALUES to a quoted x is an artefact of the
        # binning and not a value of the observable (it moved the headline
        # at x = 0.5 by 24%).  The quoted numbers are therefore evaluated
        # on the dense x grid at the bins' own kinematics, carried there
        # by the same interpolation the plot's dashed slice uses; the
        # plotted markers and the per-bin significance row stay on the
        # bin centres, where they are exact.
        q2d = np.interp(xs, xb[okb], q2b[okb])
        yd = np.interp(xs, xb[okb], yb[okb])
    camps = ((toy_b1, "crimson", "Miller $b_1$ (HERMES-like) [digitized]"),
             (b1_convolution, "navy", "CDKS convolution $b_1$ [digitized]"))
    signals, signal_x = {}, (xb[okb] if binned else xs)
    quoted, quote_x = {}, xs
    for b1f, color, label in camps:
        args_sig = (nf2, cfg0.ion, transfer, per_nucleon)
        slice_sig = azz_signal(b1f, xs, q2s, y0, *args_sig)
        if binned:
            signals[label] = azz_signal(b1f, xb[okb], q2b[okb], yb[okb],
                                        *args_sig, q2_evolve=True)
            quoted[label] = azz_signal(b1f, xs, q2d, yd, *args_sig,
                                       q2_evolve=True)
            ax.plot(xs, slice_sig, color=color, lw=0.8, alpha=0.35, ls="--")
            ax.plot(signal_x, signals[label], color=color, marker=".", ms=4,
                    lw=1.2, label=label)
        else:
            signals[label] = quoted[label] = slice_sig
            ax.plot(xs, slice_sig, color=color, label=label)
    if binned:
        ax.plot([], [], color="0.5", lw=0.8, alpha=0.6, ls="--",
                label=f"same at the fixed $Q^2$ = 4 slice")

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
                 f"transfer {args.transfer} = {factor:.4f}; signal "
                 f"{args.signal_q2} in $Q^2$)\n"
                 "3 energies combined; stat. only; "
                 f"{backends['tag'].upper()} inputs{share_note}", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    share_key = fom.run_share_tag(args.run_share)
    stem = f"money_b1_6Li_{backends['tag']}"
    if args.transfer != "rank2":        # never overwrite the published PNG
        stem = f"{stem}_{args.transfer}"
    if args.signal_q2 != "binned":      # nor the other way round
        stem = f"{stem}_q2{args.signal_q2}"
    if share_key:
        stem = f"{stem}_{share_key}"
    path = outdir / f"{stem}.png"
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")
    if binned:
        row = "  bin <Q2> [GeV^2]"
        for xq in (0.0035, 0.01, 0.03, 0.07, 0.2, 0.5):
            i = np.argmin(np.abs(xb[okb] - xq))
            row += f"  {xb[okb][i]:.4f}:{q2b[okb][i]:.2f}"
        print(row + f"   (Miller's Q2 set tops out at {MILLER_B1_Q2_REF:g})")
    for label, sig in quoted.items():
        row = "  " + label.split(" $")[0]
        for xq in (0.005, 0.01, 0.03, 0.07, 0.2, 0.5):
            row += f"  |Azz|({xq:g})={np.interp(xq, quote_x, sig):.3g}"
        print(row)
    for pzz, (x_c, err) in errors.items():
        row = f"  dAzz Pzz={pzz:g}"
        for xq in (0.0035, 0.01, 0.03, 0.28, 0.56):
            i = np.argmin(np.abs(x_c - xq))
            row += f"  {x_c[i]:.4f}:{err[i]:.3g}"
        print(row)
    x_c, err = errors[0.60]
    for label, sig in signals.items():
        s_on_bins = np.interp(x_c, signal_x, sig)
        row = "  significance (Pzz=0.6) " + label.split(" $")[0]
        for xq in (0.0035, 0.01, 0.03, 0.07, 0.2, 0.5):
            i = np.argmin(np.abs(x_c - xq))
            row += f"  {x_c[i]:.4f}:{s_on_bins[i] / err[i]:.1f}sig"
        print(row)


if __name__ == "__main__":
    main()
