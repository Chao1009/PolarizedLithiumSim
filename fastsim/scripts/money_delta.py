#!/usr/bin/env python3
"""Money plot 3 (plan step 1.3.3): gluonometry discovery reach.

5-sigma luminosity for the cos(2phi) double-helicity-flip amplitude as a
function of the Delta/F1 scale (scenario range 1e-3 .. 1e-2 spans the
Sather-Schmidt bag estimate). Assumes transversely tensor-polarized 6Li
running with unpolarized electrons; significance combines all accepted
(x, Q2) bins:  sig^2 = sum_bins A_bin^2 * P_zz^2 * N_bin / 2.

`--cl-band` (default off) adds the 95% CL exclusion contour asked for by
plans/02 Step 1.3 item 3.  sig^2 is linear in the delivered luminosity, so
the exclusion curve is the discovery curve scaled by (1.645/5)^2 = 0.108241
wherever the min-events floor does not bind; it is nevertheless re-solved
from the same per-bin terms at target = 1.645^2 rather than multiplied in,
so that a floor which DOES bind shows up instead of being papered over.
The two-contour figure is written on its own stem (`_cl95`) -- the published
single-contour PNG stays bit-for-bit.

`--run-share` is this observable's share of the programme year (plans/07
WP2), kept separate from the 10 fb^-1/u the year is worth.  L_5sigma is
quoted as the luminosity THIS OBSERVABLE must accumulate, so it is exactly
invariant under the share -- the share buys wall-clock time, not physics --
while the statistical error at a fixed programme luminosity grows as
1/sqrt(share).  Both are printed.
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.asymmetries import a_cos2phi
from polli_fastsim.inputs import get_backends
from polli_fastsim.polarized import toy_delta_gluon
from polli_fastsim.structure import NuclearF2

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# "one EIC year" in the Yellow Report accounting, the programme luminosity
# every published projection is quoted at (Report 0 section 6)
L_PROGRAMME_FB = 10.0

# significance targets as sig^2: 5 sigma for discovery, and the ONE-sided
# Gaussian 95% point z = 1.645 for exclusion (5% in one tail; two-sided,
# 1.645 is the 90% point and 95% would be 1.96).  The one-sided
# convention is the one plans/02 Step 1.3 item 3 quotes, "the same curve
# shifted by (1.645/5)^2"
Z_95 = 1.645
TARGET_5SIG = 25.0
TARGET_95 = Z_95 ** 2               # 2.706025


def cl_band_tag(cl_band):
    """Filename key for the two-contour figure ('' at the published
    default).

    The published PNG carries the 5 sigma contour alone; a run that adds
    the 95% CL contour appends this key so it cannot overwrite it -- the
    same guard as `fom.run_share_tag`, `money_cos2phi.tensor_leakage_tag`
    and `money_tagged_azz.output_stem`."""
    return "cl95" if cl_band else ""


def bin_terms(cfg, scale, pzz, base=None, run_share=1.0):
    """Per-bin (sig^2 contribution, event count), both at 1 fb^-1/nucleon
    of PROGRAMME luminosity (so both carry the factor `run_share`), over
    the accepted bins.

    Bins whose structure functions come back non-finite are dropped here,
    explicitly: a PDF grid outside its fit range returns NaN, and until
    2026-08-28 those bins were removed only as a side effect of the
    min-events comparison (NaN >= 10 is False), which stopped being true
    the moment the threshold was allowed to fall below one event."""
    sc = fom.Scenario(lumi_fb_per_nucleon=1.0, run_share=run_share,
                      pol_ion_tensor=pzz)
    nf2_in = NuclearF2(cfg.ion, base=base) if base is not None else None
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_in)
    nf2 = proj.extras["nf2"]
    f1 = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2 = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y = proj.extras["y"]
    delta = toy_delta_gluon(proj.x, proj.q2, f1, scale=scale)
    amp = a_cos2phi(delta, f1, f2, proj.x, y)
    terms = amp**2 * pzz**2 * proj.n_events / 2.0
    ok = (proj.accepted & np.isfinite(terms) & np.isfinite(proj.n_events))
    return terms[ok].ravel(), proj.n_events[ok].ravel()


def reach_from_terms(terms, n_events, min_events=10, target=TARGET_5SIG):
    """L_5sigma [fb^-1/nucleon] from the per-bin sig^2 at 1 fb^-1/u, with
    the MIN-EVENTS floor applied at the luminosity the reach is quoted at.

    The floor is there to keep a Gaussian counting significance honest, so
    it belongs at the luminosity of the measurement, not at the 1 fb^-1/u
    the sig^2 is normalised to.  Until 2026-08-28 it was applied at
    1 fb^-1/u and the reach was then scaled up from the truncated sum: a
    bin holding two events per fb^-1 -- forty at a 20 fb^-1 reach, and a
    perfectly good bin there -- was discarded.

    L sig2(L) is non-decreasing in L (raising L both scales the
    significance and admits bins), so `L sig2(L) = target` has one
    solution and it is found exactly rather than iterated: bins enter in
    order of decreasing rate, bin k at L = min_events / n_k, so a walk
    over the cumulative sum lands on the interval containing its own
    solution.  Where the target is already met at the left edge of that
    interval the answer is the edge itself -- the first luminosity at
    which the bins carrying it have min_events each."""
    n_events = np.asarray(n_events, dtype=float)
    terms = np.asarray(terms, dtype=float)
    if n_events.size == 0:
        return np.inf
    order = np.argsort(-n_events)
    n_sorted, cum = n_events[order], np.cumsum(terms[order])
    if not cum[-1] > 0:            # no bin carries any significance
        return np.inf
    enter = min_events / np.maximum(n_sorted, 1e-300)
    for k in range(n_sorted.size):
        lumi = target / cum[k] if cum[k] > 0 else np.inf
        hi = enter[k + 1] if k + 1 < n_sorted.size else np.inf
        if lumi < hi:
            return float(max(lumi, enter[k]))
    return np.inf


def sig2_per_fb_at(cfg, scale, pzz, base=None, min_events=10, lumi_fb=None,
                   run_share=1.0):
    """Significance^2 per fb^-1/nucleon of PROGRAMME luminosity at a
    reference Delta/F1 scale (linear in `run_share`).
    The toy Delta shape is linear in `scale`, so reach curves follow
    L_5sig(s) = 25 / (sig2 * (s/scale)^2) analytically -- but the
    min-events floor does not scale with it, so it is applied at
    `lumi_fb`, the luminosity the reach is evaluated at (None = the
    self-consistent one, `reach_fb`; see `reach_from_terms`)."""
    terms, n_events = bin_terms(cfg, scale, pzz, base=base,
                                run_share=run_share)
    if lumi_fb is None:
        lumi_fb = reach_from_terms(terms, n_events, min_events=min_events)
    return float(terms[n_events * lumi_fb >= min_events].sum())


def reach_fb(cfg, scale, pzz, base=None, min_events=10, target=TARGET_5SIG,
             run_share=1.0):
    """L_5sigma [fb^-1/nucleon DELIVERED to this observable] at this
    Delta/F1 scale (target = 25 is 5 sigma), with the min-events floor at
    the luminosity it is quoted at.

    Exactly invariant under `run_share`: the share scales `terms` and
    `n_events` together, `reach_from_terms` is homogeneous of degree -1 in
    both (the target ratio, the entry luminosities and their ordering all
    scale by 1/share), and multiplying the programme reach back by the
    share returns the delivered one unchanged.  The share costs programme
    luminosity -- `reach_from_terms(...)` alone, which is 1/share times as
    large -- not physics."""
    terms, n_events = bin_terms(cfg, scale, pzz, base=base,
                                run_share=run_share)
    return run_share * reach_from_terms(terms, n_events,
                                        min_events=min_events, target=target)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="6Li", choices=["6Li", "7Li"])
    ap.add_argument("--pdf", default="toy", choices=["toy", "grid"])
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--run-share", type=float, default=1.0, dest="run_share",
                    help="this observable's share of the programme year "
                         "(plans/07 WP2; default 1.0 = the whole of the "
                         "10 fb^-1/u year, which is what every published "
                         "number assumes)")
    ap.add_argument("--cl-band", action="store_true", dest="cl_band",
                    help="also draw the 95%% CL exclusion contour "
                         "(z = 1.645) under each 5-sigma curve; the "
                         "figure then goes to its own '_cl95' stem so the "
                         "published single-contour PNG is untouched")
    args = ap.parse_args()
    if not args.run_share > 0:
        ap.error("--run-share must be positive")
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("# money_delta.py  ion=%s  pdf=%s  %s"
          % (args.ion, args.pdf,
             fom.run_share_header(L_PROGRAMME_FB, args.run_share)))

    backends = get_backends(args.pdf)
    scales = np.logspace(-3.3, -1.7, 15)
    s0 = 1e-3
    fig, ax = plt.subplots(figsize=(7, 5))
    reach_ref, reach_ref95 = {}, {}
    for cfg, color in zip(beams.default_configs(args.ion),
                          ("crimson", "seagreen", "navy")):
        for pzz, ls in ((0.60, "-"), (0.80, "--")):
            # the amplitude is linear in the scale, so the per-bin sig^2
            # terms scale as (s/s0)^2 and only the min-events floor has to
            # be re-solved per point
            terms, n_events = bin_terms(cfg, s0, pzz, base=backends["base"],
                                        run_share=args.run_share)
            # delivered luminosity: share-invariant, so the plotted curve
            # is the published one at any share (see `reach_fb`)
            reach = args.run_share * np.array(
                [reach_from_terms(terms * (s / s0) ** 2, n_events)
                 for s in scales])
            reach_ref[(cfg.label(), pzz)] = args.run_share * reach_from_terms(
                terms, n_events)
            ax.plot(scales, reach, ls, color=color, lw=1.5,
                    label=f"{cfg.label()}, $P_{{zz}}$={pzz:g}")
            if args.cl_band:
                # re-solved at target = z^2, not scaled from `reach`, so a
                # binding min-events floor would be visible as a departure
                # from the flat (1.645/5)^2 offset
                reach95 = args.run_share * np.array(
                    [reach_from_terms(terms * (s / s0) ** 2, n_events,
                                      target=TARGET_95)
                     for s in scales])
                reach_ref95[(cfg.label(), pzz)] = (
                    args.run_share * reach_from_terms(terms, n_events,
                                                      target=TARGET_95))
                ax.plot(scales, reach95, ls, color=color, lw=0.9, alpha=0.75)
                if pzz == 0.80:     # shade the band at the spec P_zz only:
                    # six overlapping fills are a grey mass, and 0.80 is
                    # the source requirement the reach is quoted at
                    ax.fill_between(scales, reach95, reach, color=color,
                                    alpha=0.10, lw=0)
    if args.cl_band:
        # one proxy entry for the six thin curves, so the legend keeps its
        # published six rows plus one
        ax.plot([], [], "-", color="0.35", lw=0.9, alpha=0.75,
                label=r"95% CL exclusion ($z=1.645$)")
    ax.axhspan(1, 100, color="gold", alpha=0.12,
               label="1-100 fb$^{-1}$/u (plausible program)")
    ax.axvline(1e-3, color="gray", ls=":", lw=1)
    ax.text(1.05e-3, 0.93, "Sather-Schmidt\n$O(10^{-3})$", fontsize=7,
            transform=ax.get_xaxis_transform(), va="top")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta/F_1$ scale (peak of scenario shape)")
    ax.set_ylabel(r"$L_{5\sigma}$ (and $L_{95\%}$) [fb$^{-1}$/nucleon]"
                  if args.cl_band else
                  r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]")
    share_note = ("" if args.run_share == 1.0
                  else f"; run share {args.run_share:g} of the programme year")
    # kept short: the second title line is already at the figure width
    cl_note = ("" if not args.cl_band
               else "; thin: 95% CL, $z$ = 1.645")
    ax.set_title(f"Nuclear gluonometry reach, transversely polarized "
                 f"{args.ion}\n(cos 2$\\phi$ amplitude, all bins combined; "
                 f"{backends['tag'].upper()} inputs{share_note}"
                 f"{cl_note})", fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    stem = f"money_delta_{args.ion}_{backends['tag']}"
    for key in (fom.run_share_tag(args.run_share), cl_band_tag(args.cl_band)):
        if key:                         # never overwrite the published PNG
            stem = f"{stem}_{key}"
    path = outdir / f"{stem}.png"
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")
    for cfg in beams.default_configs(args.ion):
        l5 = reach_ref[(cfg.label(), 0.8)]
        l5_hi = reach_fb(cfg, 1e-2, 0.8, base=backends["base"],
                         run_share=args.run_share)
        # error on the Delta/F1 scale after one programme year at this
        # share.  Derived from the reach itself rather than from a second
        # sum over the bins, so the two printed numbers cannot disagree:
        # sig^2 is linear in the delivered luminosity and equals 25 at
        # L_5sigma, and the toy amplitude is linear in the scale, so
        # sig = 5 sqrt(L_delivered / L_5sigma) and delta = 1e-3 / sig.
        lumi_eff = L_PROGRAMME_FB * args.run_share
        d_scale = (2e-4 * np.sqrt(l5 / lumi_eff) if np.isfinite(l5)
                   else np.inf)
        print(f"  {cfg.label():26s} L_5sig(Delta/F1=1e-3, Pzz=0.8) = "
              f"{l5:9.1f} fb^-1/u ; (1e-2) = {l5_hi:7.3f} ; "
              f"delta(Delta/F1) after {L_PROGRAMME_FB:g} fb^-1/u x share = "
              f"{d_scale:9.3e} ; programme years to 5 sigma = "
              f"{l5 / lumi_eff:8.2f}")
        if args.cl_band:
            l95 = reach_ref95[(cfg.label(), 0.8)]
            ratio = l95 / l5 if np.isfinite(l5) and l5 > 0 else float("nan")
            print(f"  {'':26s} L_95%(Delta/F1=1e-3, Pzz=0.8)  = "
                  f"{l95:9.3f} fb^-1/u ; L_95%/L_5sig = {ratio:.5f} "
                  f"(expected (1.645/5)^2 = {TARGET_95 / TARGET_5SIG:.5f})")


if __name__ == "__main__":
    main()
