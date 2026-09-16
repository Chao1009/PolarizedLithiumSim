#!/usr/bin/env python3
"""Letter Figure 3 -- the extracted x*Delta(x, Q2), condensed from money plot 7
(`evgen/scripts/money_delta_extraction.py`,
`evgen/money_delta_extracted_6Li.png`).

plans/07 WP6 asks for "money plot 7 -> Fig. 3 (two Q2 slices)".  The published
figure carries the three sweet-spot Q2 slices; this one carries the two that
make the letter's point -- the lowest, Q2 = 1.14 GeV^2, where the projection is
most precise, and the highest, Q2 = 14.3 GeV^2, where the moment-constrained
model and the no-F1 variant cross and the data separate them.  The middle slice
is still measured, in its place in the loop, because it is a published row of
Report 1 Section 5.1 and because the pseudo-random stream is
`money_delta_extraction`'s own and must not be shortened.

All machinery is imported: the sweet-spot selector, the super-bin masks and
`measure` come from `money_cos2phi` exactly as the producing script imports
them, the Delta model and the bin-centering factor K are that script's, and the
published settings -- seed 20260811 among them -- are read off its own parser.

Usage:  python3 paper/fig3_extraction.py          (writes figs/fig3.*)
"""

import numpy as np

import _sources as src
import figstyle as fs

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402

import money_cos2phi as m5                                       # noqa: E402
import money_delta_extraction as m7                              # noqa: E402

from polligen import bookkeeping as bk                           # noqa: E402
from polligen.sample import InclusiveSampler                     # noqa: E402
from polligen.xsec import InclusiveKernel                        # noqa: E402

from polli_fastsim import beams, delta_models as dm, fom         # noqa: E402
from polli_fastsim.kinematics import kinematic_mask              # noqa: E402
from polli_fastsim.polarized import toy_b1                       # noqa: E402

# The two Q2 slices the letter panel shows, in GeV^2 as the producing script
# rounds them.  All three are measured.
SHOWN = (1.14, 14.3)


def measure_money7():
    """Re-run money plot 7's extraction, in its order, at its defaults."""
    args = src.published_args(m7.main)

    config = beams.default_configs("6Li")[args.config]
    lumi1_pb = args.lumi_1yr * 1e3
    lumi10_pb = args.lumi_10yr * 1e3
    rng = np.random.default_rng(args.seed)

    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    backends = m5.pdf_backends(args, config.ion)
    model, _q2_ref = m5.build_delta_model(args, config, scenario,
                                          backends=backends)
    b3_func, b4_func = m5.b34_funcs(args)
    kern = InclusiveKernel(config.ion, b1_func=toy_b1, delta_func=model,
                           b3_func=b3_func, b4_func=b4_func,
                           tensor_gamma=args.tensor_gamma,
                           f2_source=backends["base"],
                           g1_model=backends["g1"],
                           nuclear_f2=backends["nuclear"],
                           r_func=backends["nuclear"].r_func)

    proj = fom.project_rates(config, scenario,
                             nuclear_f2=backends["nuclear"])
    obs = fom.project_observables(config, scenario, proj,
                                  kern.g1_model, toy_b1, model)
    spots = m5.pick_sweet_spots_banded(proj, obs["sig_a_cos2phi"])
    q2_slices = sorted({round(q2, 3) for _x, q2, _i, _j in spots})

    sampler = InclusiveSampler(kern, config, scenario, nx=60, nq2=45)
    plan = bk.transverse_tensor_plan(args.pzz)
    cat = plan.categories[0]
    alt = dm.make("moment_B", variant=args.variant, dilution=args.dilution)

    xe = proj.x_edges
    slices = []
    for q2s in q2_slices:
        pts = []
        for i0 in range(0, xe.size - 2, 2):        # merge pairs of x bins
            xc = np.sqrt(xe[i0] * xe[i0 + 2])
            if not kinematic_mask(xc, q2s, sampler.s):
                continue
            mask = m5.superbin_mask(sampler, xe[i0], xe[i0 + 2],
                                    q2s / 1.6, q2s * 1.6)
            if not mask.any():
                continue
            m1 = m5.measure(sampler, cat, mask, lumi1_pb, plan.pzz_true, rng)
            if m1["n"] < 1e3 or m1["err"] > 8e-3:
                continue
            if abs(m1["truth"]) < 1e-5:
                continue
            m10 = m5.measure(sampler, cat, mask, lumi10_pb, plan.pzz_true,
                             rng)
            f1c = kern.nf2.f1a(xc, q2s) / kern.ion.A
            delta_c = float(model(xc, q2s, f1c))
            k_conv = delta_c / m1["truth"]
            pts.append({"x": float(xc), "delta_c": delta_c,
                        "d1": float(m1["amp"] * k_conv),
                        "e1": float(m1["err"] * abs(k_conv)),
                        "d10": float(m10["amp"] * k_conv),
                        "e10": float(m10["err"] * abs(k_conv))})
        # the model curves of this slice, on the producing script's grid
        xg = np.logspace(np.log10(2e-4), np.log10(0.8), 250)
        q2g = np.full_like(xg, q2s)
        ok = kinematic_mask(xg, q2g, sampler.s)
        f1g = kern.nf2.f1a(xg, q2g) / kern.ion.A
        slices.append({
            "q2": float(q2s), "points": pts,
            "curve_x": xg[ok],
            "curve_a": (xg * np.asarray(model(xg, q2g, f1g)))[ok],
            "curve_b": (xg * np.asarray(alt(xg, q2g, f1g)))[ok],
        })
    info = {"model": model.info(), "pzz": float(plan.pzz_true),
            "backend": m5.describe_backends(args, backends),
            "label": config.label()}
    return args, info, slices


def draw(args, info, slices):
    shown = [s for s in slices if round(s["q2"], 2) in
             tuple(round(v, 2) for v in SHOWN)]
    fig, axes = fs.subplots(len(shown), 1, width=fs.COLUMN_W, height=4.15,
                            sharey=True)
    axes = np.atleast_1d(axes)
    xoff = 1.045
    handles = None
    for k, (ax, s) in enumerate(zip(axes, shown)):
        pts = s["points"]
        ax.errorbar([p["x"] for p in pts],
                    [1e3 * p["x"] * p["d1"] for p in pts],
                    yerr=[1e3 * p["x"] * p["e1"] for p in pts], fmt="s",
                    mfc="none", color=fs.FIT, ms=3.0, capsize=1.4, lw=0.8,
                    zorder=3, label=r"1 yr, %g fb$^{-1}$/u" % args.lumi_1yr)
        ax.errorbar([p["x"] * xoff for p in pts],
                    [1e3 * p["x"] * p["d10"] for p in pts],
                    yerr=[1e3 * p["x"] * p["e10"] for p in pts], fmt="o",
                    color=fs.DATA, ms=3.0, capsize=1.4, lw=0.8, zorder=4,
                    label=r"10 yr, %g fb$^{-1}$/u" % args.lumi_10yr)
        ax.plot(s["curve_x"], 1e3 * s["curve_a"], "-", color=fs.TRUTH,
                lw=1.3, label="moment-constrained (A)")
        ax.plot(s["curve_x"], 1e3 * s["curve_b"], "-", color=fs.ALT,
                lw=1.1, label="no-$F_1$ (B)")
        ax.set_xscale("log")
        ax.axhline(0.0, color="0.85", lw=0.5, zorder=0)
        ax.tick_params(pad=2)
        if k == 0:
            handles, labels = ax.get_legend_handles_labels()
        fs.panel_label(ax, "(%s)  $Q^2 = %.3g$ GeV$^2$"
                       % ("ab"[k], s["q2"]), x=0.035, y=0.93,
                       fontweight="normal")
    # one y range for both panels, with room above zero for the panel tags
    lo = min(ax.get_ylim()[0] for ax in axes)
    axes[0].set_ylim(lo, 0.85)
    axes[-1].set_xlabel(r"$x$", labelpad=1)
    fig.supylabel(r"$x\,\Delta(x,Q^2)$ per nucleon  $[10^{-3}]$",
                  fontsize=fs.FONT_LARGE, x=0.035)
    fig.legend(handles, labels, loc="upper center", ncol=2,
               bbox_to_anchor=(0.57, 1.005), fontsize=fs.FONT_SMALL,
               columnspacing=0.9)
    fig.subplots_adjust(left=0.185, right=0.985, bottom=0.085, top=0.895,
                        hspace=0.20)
    return fig


def caption(args, info, slices):
    """Figure 3's caption, from the same extraction the panels draw."""
    shown = [s for s in slices if round(s["q2"], 2) in
             tuple(round(v, 2) for v in SHOWN)]
    best = []
    for s in shown:
        b = min(s["points"], key=lambda p: p["e10"])
        best.append((s["q2"], b["x"], b["delta_c"], b["e1"],
                     100 * abs(b["e1"] / b["delta_c"])))
    return (
        r"The double-helicity-flip structure function extracted from the "
        r"amplitudes of Fig.~\ref{fig:inclusive}, at %s with "
        r"$P_{zz} = %.2f$: $x\Delta(x, Q^{2})$ per nucleon in the "
        r"$Q^{2} = %.3g$\,GeV$^{2}$ (a) and $Q^{2} = %.3g$\,GeV$^{2}$ (b) "
        r"slices, with independent one-year (open, %g\,fb$^{-1}$ per "
        r"nucleon) and ten-year (filled, %g\,fb$^{-1}$ per nucleon) "
        r"pseudo-measurements against the injected moment-constrained "
        r"$\Delta$ and the no-$F_{1}$ variant. Each point is the fitted "
        r"amplitude of a pair of merged $x$ bins converted to a point value "
        r"by the model bin-centring factor "
        r"$K = \Delta_{\rm model}(x_{c}, Q^{2}_{c})/A_{\rm model}$, so "
        r"that $\delta\Delta = \delta\hat A\,|K|$. The best bin of each "
        r"slice reaches $\Delta = %s \pm %s$ at $x = %.3g$ and "
        r"$%s \pm %s$ at $x = %.3g$ in one year, %s\%% and %s\%% "
        r"relative. The two curves are separated bin by bin at low "
        r"$Q^{2}$ and cross at large $x$ in the high-$Q^{2}$ slice, which "
        r"is what makes the $x$ dependence, rather than the size of the "
        r"amplitude, the discriminating measurement. Statistical errors "
        r"only; the third sweet-spot slice, $Q^{2} = %.3g$\,GeV$^{2}$, is "
        r"measured with the same draw and is not drawn here."
        % (src.beam_math(info["label"]), info["pzz"],
           best[0][0], best[1][0], args.lumi_1yr, args.lumi_10yr,
           src.dec(best[0][2], 3), src.dec(best[0][3], 3), best[0][1],
           src.dec(best[1][2], 4), src.dec(best[1][3], 4), best[1][1],
           src.dec(best[0][4], 1), src.dec(best[1][4], 1),
           [s["q2"] for s in slices
            if s["q2"] not in [b[0] for b in best]][0]))


def main():
    args, info, slices = measure_money7()
    with fs.rc():
        fig = draw(args, info, slices)
        fs.save(fig, "fig3")
        plt.close(fig)
    src.write_caption("fig3", "figthreecaption", caption(args, info, slices))

    payload = {
        "figure": "fig3",
        "source_script": "evgen/scripts/money_delta_extraction.py",
        "published_figure": "evgen/money_delta_extracted_6Li.png",
        "settings": src.describe(args),
        "beam": info["label"], "pzz": info["pzz"],
        "delta_model": info["model"], "backend": info["backend"],
        "panels_shown": list(SHOWN),
        "slices": [],
    }
    for s in slices:
        pts = s["points"]
        best = min(pts, key=lambda p: p["e10"]) if pts else None
        payload["slices"].append({
            "q2": s["q2"], "n_points": len(pts), "points": pts,
            "best": None if best is None else {
                "x": best["x"], "delta_c": best["delta_c"],
                "err_1yr": best["e1"], "err_10yr": best["e10"],
                "relative_1yr": abs(best["e1"] / best["delta_c"])},
        })
    src.dump("fig3", payload)

    print("beam: %s, P_zz = %.2f, %s"
          % (info["label"], info["pzz"], info["backend"]))
    for s in payload["slices"]:
        b = s["best"]
        print("Q2=%.3g: %d points; best x=%.3g Delta=%+.4f +- %.4f (1yr) "
              "+- %.4f (10yr), %.1f%% relative"
              % (s["q2"], s["n_points"], b["x"], b["delta_c"],
                 b["err_1yr"], b["err_10yr"], 100 * b["relative_1yr"]))


if __name__ == "__main__":
    main()
