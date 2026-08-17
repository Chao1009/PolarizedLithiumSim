#!/usr/bin/env python3
"""Phase-space companion figure: where the events are and how the
detailed money plots bin them.

Two (x, Q2) event-rate maps on the 40x30 log-log analysis grid at the
1-year EIC program (10 fb^-1/u), with every bin used by the detailed
cos 2phi / Delta plots drawn on top:

* LEFT -- inclusive DIS (scattered-electron acceptance applied): the
  rate map behind money plots 5/7.  Overlaid: the four sweet-spot
  super-bins (3x3 analysis cells, picked per Q2 band from the
  significance map -- the phi' panels of money plot 5), and the three
  sweet-spot Q2 slices (Q2 x [1/1.6, 1.6]) with their merged x-bin
  pairs -- one comb cell = one data point of the amplitude-vs-x panel
  (money plot 5, right) and of the Delta extraction (money plot 7).
  Comb cells are kept by the exact money-plot selection (N >= 1e3,
  dA <= 8e-3, |A_truth| >= 1e-5, evaluated on the same 60x45 sampler
  grid the money plots use).
* RIGHT -- RP-tagged coherent channel e 6Li -> e' X 6Li(g.s.): the same
  DIS map times the scenario coherent fraction f_coh(x) (f0 = 0.04,
  dying at x_coh = 0.01) times the analytic tag acceptance
  exp(-B pT_cut^2) = 13.5% (high-acceptance optics) -- the rate map
  behind money plot 6.  Overlaid: the tagged-count-maximal 3x3
  super-bin whose phi' modulation is money plot 6(d); money plot 6(c)
  is this map summed over Q2, column by column.

Color = expected events per bin (log scale), drawn where >= 1 event is
expected.  Conventions as everywhere else: per-nucleon luminosity, whole
luminosity in the transverse-tensor fill, TOY structure functions.

Usage:  python3 scripts/phase_space_bins.py
"""

import argparse
import pathlib
import sys

import numpy as np

_SCRIPTS = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS.parent))
sys.path.insert(0, str(_SCRIPTS))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from money_cos2phi import (build_delta_model, measure,  # noqa: E402
                           pick_sweet_spots_banded, superbin_edges,
                           superbin_mask)
from money_cos2phi_coherent import best_superbin  # noqa: E402

from polligen import bookkeeping as bk  # noqa: E402
from polligen import coherent as coh  # noqa: E402
from polligen.sample import InclusiveSampler  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, delta_models as dm, fom  # noqa: E402
from polli_fastsim.farforward import HIGH_ACCEPTANCE  # noqa: E402
from polli_fastsim.kinematics import kinematic_mask  # noqa: E402
from polli_fastsim.polarized import toy_b1  # noqa: E402

C_BIN = "#D55E00"   # vermillion: analysis bins (cased in white for
C_COMB = "0.15"     # contrast on any viridis level); combs dark gray
HALO = [pe.withStroke(linewidth=2.2, foreground="white")]


def sci(v, sig=2):
    """Mathtext scientific notation: 4.6e9 -> 4.6\\times10^{9}."""
    exp = int(np.floor(np.log10(abs(v))))
    mant = v / 10.0 ** exp
    return r"%.*g\times10^{%d}" % (sig, mant, exp)


def cased_rect(ax, xlo, xhi, q2lo, q2hi, color=C_BIN, lw=1.5, zorder=6):
    """Rectangle with a white casing so it reads on any map color."""
    for c, w, dz in (("white", lw + 1.4, 0), (color, lw, 1)):
        ax.add_patch(Rectangle((xlo, q2lo), xhi - xlo, q2hi - q2lo,
                               fill=False, edgecolor=c, linewidth=w,
                               zorder=zorder + dz, joinstyle="miter"))


def draw_map(ax, proj, values, cbar_label, fig):
    vals = np.ma.masked_where(~proj.accepted | (values < 1.0), values)
    pcm = ax.pcolormesh(proj.x_edges, proj.q2_edges, vals.T,
                        norm=LogNorm(vmin=1.0), cmap="viridis",
                        shading="auto", rasterized=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(proj.x_edges[0], proj.x_edges[-1])
    ax.set_ylim(proj.q2_edges[0], proj.q2_edges[-1])
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$Q^2$  [GeV$^2$]")
    ax.tick_params(labelsize=8)
    cb = fig.colorbar(pcm, ax=ax, pad=0.012)
    cb.set_label(cbar_label, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    return pcm


def slope_label(ax, x_at, func, text, above=True, color="0.35"):
    """Small label rotated along the curve q2 = func(x) at x = x_at."""
    x2 = x_at * 1.3
    p = ax.transData.transform([(x_at, func(x_at)), (x2, func(x2))])
    ang = np.degrees(np.arctan2(p[1, 1] - p[0, 1], p[1, 0] - p[0, 0]))
    off = 1.25 if above else 1.0 / 1.25
    ax.annotate(text, xy=(x_at, func(x_at) * off), fontsize=6.5,
                color=color, rotation=ang, rotation_mode="anchor",
                ha="left", va="center", zorder=4, path_effects=HALO)


def draw_guides(ax, s, scenario):
    """Dashed guide lines for the y and W2 cuts bounding the acceptance."""
    xg = np.logspace(-4.0, np.log10(0.999), 300)
    m_n2 = 0.9383 ** 2
    lines = [
        (lambda x: scenario.y_max * s * x,
         r"$y=%.2f$" % scenario.y_max, 0.05, True),
        (lambda x: scenario.y_min * s * x,
         r"$y=%.2f$" % scenario.y_min, 0.28, False),
        (lambda x: (scenario.w2_min - m_n2) * x / (1.0 - x),
         r"$W^2=%.0f$ GeV$^2$" % scenario.w2_min, 0.42, True),
    ]
    for func, text, x_at, above in lines:
        q2 = func(xg)
        ax.plot(xg, q2, "--", color="0.55", lw=0.8, zorder=3)
        slope_label(ax, x_at, func, text, above=above)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2),
                    help="beam-energy point (default mid, e10 x 6Li50)")
    ap.add_argument("--delta-model", default="moment_A",
                    choices=dm.available())
    ap.add_argument("--variant", default="mid_x",
                    choices=sorted(dm.VARIANTS))
    ap.add_argument("--dilution", type=float, default=1.0 / 3.0)
    ap.add_argument("--scale", type=float, default=1e-2,
                    help="peak Delta/F1 (toy model only)")
    ap.add_argument("--lumi-1yr", type=float, default=10.0,
                    help="1-year EIC program [fb^-1/nucleon]")
    ap.add_argument("--pzz", type=float, default=0.60)
    ap.add_argument("--seed", type=int, default=20260817)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs("6Li")[args.config]
    lumi1_pb = args.lumi_1yr * 1e3
    rng = np.random.default_rng(args.seed)

    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    model, _q2_ref = build_delta_model(args, config, scenario)
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1, delta_func=model)

    # --- inclusive: rate map + the money-plot 5/7 bins -------------------
    proj = fom.project_rates(config, scenario)
    obs = fom.project_observables(config, scenario, proj,
                                  kern.g1_model, toy_b1, model)
    spots = pick_sweet_spots_banded(proj, obs["sig_a_cos2phi"])[:4]
    q2_slices = sorted({round(q2, 3) for _x, q2, _i, _j in spots})

    sampler = InclusiveSampler(kern, config, scenario, nx=60, nq2=45)
    plan = bk.transverse_tensor_plan(args.pzz)
    cat = plan.categories[0]

    # the x-bin combs of the Delta-extraction points, by the exact
    # money-plot selection (money_cos2phi right panel / money plot 7)
    xe = proj.x_edges
    slice_bins = {}
    for q2s in q2_slices:
        kept = []
        for i0 in range(0, xe.size - 2, 2):  # merge pairs of x bins
            xc = np.sqrt(xe[i0] * xe[i0 + 2])
            if not kinematic_mask(xc, q2s, sampler.s):
                continue
            mask = superbin_mask(sampler, xe[i0], xe[i0 + 2],
                                 q2s / 1.6, q2s * 1.6)
            if not mask.any():
                continue
            m1 = measure(sampler, cat, mask, lumi1_pb, plan.pzz_true, rng)
            if m1["n"] < 1e3 or m1["err"] > 8e-3:
                continue
            if abs(m1["truth"]) < 1e-5:
                continue
            kept.append((xe[i0], xe[i0 + 2]))
        slice_bins[q2s] = kept

    # --- tagged coherent: rate map + the money-plot 6(d) bin -------------
    sc = coh.CoherentScenario()
    proj_c, n_coh, tagged = coh.project_coherent(
        config, scenario, sc, optics_list=(HIGH_ACCEPTANCE,))
    n_tag = tagged[HIGH_ACCEPTANCE.name]
    sel, txlo, txhi, tq2lo, tq2hi = best_superbin(proj_c, n_tag)

    # --- figure ----------------------------------------------------------
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(12.6, 5.1))
    draw_map(axa, proj, proj.n_events,
             "inclusive DIS events / bin (1 yr)", fig)
    draw_map(axb, proj_c, n_tag,
             "RP-tagged coherent events / bin (1 yr)", fig)
    for ax in (axa, axb):
        draw_guides(ax, sampler.s, scenario)

    # (a) Q2-slice combs: the Delta-extraction x bins
    for q2s, kept in slice_bins.items():
        for xlo, xhi in kept:
            axa.add_patch(Rectangle(
                (xlo, q2s / 1.6), xhi - xlo, q2s * 1.6 - q2s / 1.6,
                fill=False, edgecolor=C_COMB, linewidth=0.55,
                zorder=5))
        if kept:
            axa.annotate(r"$Q^2\!\approx\!%.3g$" % q2s,
                         xy=(kept[0][0] * 0.88, q2s), fontsize=6.5,
                         color=C_COMB, ha="right", va="center",
                         zorder=5, path_effects=HALO)

    # (a) sweet-spot super-bins, numbered like the money-plot 5 panels
    for k, (_xs, _qs, i, j) in enumerate(spots):
        xlo, xhi, q2lo, q2hi = superbin_edges(proj, i, j)
        cased_rect(axa, xlo, xhi, q2lo, q2hi)
        axa.annotate("%d" % (k + 1), xy=(np.sqrt(xlo * xhi), q2hi * 1.25),
                     fontsize=8, fontweight="bold", color=C_BIN,
                     ha="center", va="bottom", zorder=7,
                     path_effects=HALO)

    axa.annotate(r"$N_{\rm DIS}=%s$ (1 yr)" % sci(proj.n_events.sum()),
                 xy=(0.045, 0.615), xycoords="axes fraction", fontsize=8,
                 zorder=7, path_effects=HALO)
    axa.set_title("inclusive DIS  (money plots 5/7)", fontsize=9.5)
    axa.legend(handles=[
        Line2D([], [], color=C_BIN, lw=1.5,
               label=r"$\phi'$ super-bins 1–4 (plot 5)"),
        Line2D([], [], color=C_COMB, lw=0.7,
               label=r"$\Delta$-extraction $x$ bins (plots 5R/7)"),
        Line2D([], [], color="0.55", lw=0.8, ls="--",
               label=r"$y$ / $W^2$ cut boundaries"),
    ], loc="upper left", fontsize=6.5, framealpha=0.92,
        borderpad=0.6, handlelength=1.6)

    # (b) tagged super-bin + coherence scale
    cased_rect(axb, txlo, txhi, tq2lo, tq2hi)
    axb.annotate(r"$\phi'$ bin (plot 6d)",
                 xy=(np.sqrt(txlo * txhi), tq2hi * 1.25), fontsize=7,
                 color=C_BIN, ha="center", va="bottom", zorder=7,
                 path_effects=HALO)
    axb.axvline(sc.x_coh, color="0.35", ls=":", lw=0.9, zorder=4)
    axb.annotate(r"$x_{\rm coh}$: coherence dies",
                 xy=(sc.x_coh * 1.15, 250), fontsize=6.5, color="0.35",
                 rotation=90, va="top", ha="left", zorder=4,
                 path_effects=HALO)
    axb.annotate(
        (r"$N_{\rm tag}=%s$ (1 yr)" "\n"
         r"acc $= %.1f$%% (HA optics), $f_0=%.2f$")
        % (sci(n_tag.sum()), 100 * sc.tag_acceptance(
            HIGH_ACCEPTANCE.pt_cut_near_beam), sc.f0),
        xy=(0.985, 0.035), xycoords="axes fraction", fontsize=8,
        ha="right", zorder=7, path_effects=HALO)
    axb.set_title(r"RP-tagged coherent, $e\,^6$Li$\to e'X\,^6$Li(g.s.)"
                  "  (money plot 6)", fontsize=9.5)

    fig.suptitle(
        r"$(x, Q^2)$ phase space and analysis binning — %s, 1-year EIC "
        r"program (%g fb$^{-1}$/u), $40\times30$ log-log grid"
        "\n"
        r"color: expected events per bin where $\geq 1$; $e'$ cuts "
        r"$|\eta|\leq%.1f$, $E'\geq%.1f$ GeV, $%.2f\leq y\leq%.2f$, "
        r"$W^2\geq%.0f$ GeV$^2$"
        % (config.label(), args.lumi_1yr, scenario.eta_max,
           scenario.e_prime_min, scenario.y_min, scenario.y_max,
           scenario.w2_min), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "phase_space_bins_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)

    print("N_DIS (1 yr) = %.3e   N_coh = %.3e   N_tag (HA) = %.3e"
          % (proj.n_events.sum(), n_coh.sum(), n_tag.sum()))
    for k, (xs, qs, i, j) in enumerate(spots):
        xlo, xhi, q2lo, q2hi = superbin_edges(proj, i, j)
        print("super-bin %d: x [%.3g, %.3g]  Q2 [%.3g, %.3g]  "
              "(center %.3g, %.3g)" % (k + 1, xlo, xhi, q2lo, q2hi,
                                       xs, qs))
    for q2s in q2_slices:
        kept = slice_bins[q2s]
        print("slice Q2=%.3g [%.3g, %.3g]: %d x bins, x [%.3g, %.3g]"
              % (q2s, q2s / 1.6, q2s * 1.6, len(kept),
                 kept[0][0] if kept else np.nan,
                 kept[-1][1] if kept else np.nan))
    print("tagged super-bin: x [%.3g, %.3g]  Q2 [%.3g, %.3g]  N=%.3g"
          % (txlo, txhi, tq2lo, tq2hi, n_tag[sel].sum()))


if __name__ == "__main__":
    main()
