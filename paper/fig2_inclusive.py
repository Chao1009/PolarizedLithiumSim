#!/usr/bin/env python3
"""Letter Figure 2 -- the inclusive cos 2phi' measurement, condensed from
money plot 5 (`evgen/scripts/money_cos2phi.py`, `evgen/money_cos2phi_6Li.png`).

plans/07 WP6 asks for "money plot 5 -> Fig. 2 (two phi' panels + amplitude vs
x)".  The published figure carries four phi' panels; this one carries the two
that span the accepted phase space -- super-bin 1, the highest-significance
bin, and super-bin 4, the highest-Q2 one -- and the amplitude-versus-x panel
beside them.  The other two super-bins are still MEASURED here, in the order
money plot 5 measures them, because they are rows of Table 1 and because the
pseudo-random stream must be the one that script draws: skipping a
measurement would shift every later draw and the letter would show a
different pseudo-experiment from the report.

Everything numerical comes from `money_cos2phi`: its sweet-spot selector, its
super-bin masks, its `measure`, its error formula, its Delta-model builder,
and its published defaults read off its own parser (`_sources.published_args`),
seed 20260810 included.  Nothing is reimplemented and no setting is restated.

The four rows of Table 1 are written to `figs/fig2.json`; `table1.py` reads
them from there.

Usage:  python3 paper/fig2_inclusive.py            (writes figs/fig2.*)
"""

import numpy as np

import _sources as src
import figstyle as fs

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
from matplotlib.gridspec import GridSpec                         # noqa: E402

import money_cos2phi as m5                                       # noqa: E402

from polligen import bookkeeping as bk                           # noqa: E402
from polligen.estimators import cos2phi_fit_err                  # noqa: E402
from polligen.sample import InclusiveSampler                     # noqa: E402
from polligen.xsec import InclusiveKernel                        # noqa: E402

from polli_fastsim import beams, delta_models as dm, fom         # noqa: E402
from polli_fastsim.asymmetries import a_cos2phi                  # noqa: E402
from polli_fastsim.kinematics import kinematic_mask, y_from_xq2  # noqa: E402
from polli_fastsim.polarized import toy_b1                       # noqa: E402

# The two super-bins the letter panel shows, 1-based as money plot 5 numbers
# them: the best bin and the highest-Q2 bin.  All four are measured.
SHOWN = (1, 4)


def measure_money5():
    """Re-run money plot 5's measurement, in its order, at its defaults.

    Returns (args, config, spots, per-spot results, amplitude-vs-x points,
    model curves).  The body follows `money_cos2phi.main()` line for line
    down to the calls it makes; every quantity is that module's.
    """
    args = src.published_args(m5.main)

    config = beams.default_configs("6Li")[args.config]
    lumi1_pb = args.lumi_1yr * 1e3
    lumi10_pb = args.lumi_10yr * 1e3
    rng = np.random.default_rng(args.seed)

    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    backends = m5.pdf_backends(args, config.ion)
    model, q2_ref = m5.build_delta_model(args, config, scenario,
                                         backends=backends)

    proj = fom.project_rates(config, scenario,
                             nuclear_f2=backends["nuclear"])
    b3_func, b4_func = m5.b34_funcs(args)
    kern = InclusiveKernel(config.ion, b1_func=toy_b1, delta_func=model,
                           b3_func=b3_func, b4_func=b4_func,
                           tensor_gamma=args.tensor_gamma,
                           f2_source=backends["base"],
                           g1_model=backends["g1"],
                           nuclear_f2=backends["nuclear"],
                           r_func=backends["nuclear"].r_func)
    obs = fom.project_observables(config, scenario, proj,
                                  kern.g1_model, toy_b1, model)
    spots = m5.pick_sweet_spots_banded(proj, obs["sig_a_cos2phi"])[:args.nspots]

    sampler = InclusiveSampler(kern, config, scenario, nx=60, nq2=45)
    plan = bk.transverse_tensor_plan(args.pzz)
    cat = plan.categories[0]

    # --- the four super-bins, in money plot 5's order --------------------
    results = []
    for k, (xs, qs, i, j) in enumerate(spots):
        xlo, xhi, q2lo, q2hi = m5.superbin_edges(proj, i, j)
        mask = m5.superbin_mask(sampler, xlo, xhi, q2lo, q2hi)
        m = m5.measure(sampler, cat, mask, lumi1_pb, plan.pzz_true, rng)
        m["err10"] = cos2phi_fit_err(lumi10_pb * m["sigma_pb"], plan.pzz_true)
        m.update(spot=k + 1, x=float(xs), q2=float(qs),
                 x_lo=float(xlo), x_hi=float(xhi),
                 q2_lo=float(q2lo), q2_hi=float(q2hi))
        results.append(m)

    # --- amplitude versus x along the sweet-spot Q2 slice ----------------
    q2_spot = spots[0][1]
    q2lo, q2hi = q2_spot / 1.6, q2_spot * 1.6
    xe = proj.x_edges
    pts = []
    for i0 in range(0, xe.size - 2, 2):            # merge pairs of x bins
        xc = np.sqrt(xe[i0] * xe[i0 + 2])
        if not kinematic_mask(xc, q2_spot, sampler.s):
            continue
        mask = m5.superbin_mask(sampler, xe[i0], xe[i0 + 2], q2lo, q2hi)
        if not mask.any():
            continue
        m = m5.measure(sampler, cat, mask, lumi1_pb, plan.pzz_true, rng)
        if m["n"] < 1e3 or m["err"] > 8e-3:
            continue
        m10 = m5.measure(sampler, cat, mask, lumi10_pb, plan.pzz_true, rng)
        pts.append({"x": float(xc), "amp": float(m["amp"]),
                    "err": float(m["err"]), "amp10": float(m10["amp"]),
                    "err10": float(m10["err"])})

    # --- the model curves of the same panel ------------------------------
    xg = np.logspace(np.log10(2e-4), np.log10(0.8), 250)
    q2g = np.full_like(xg, q2_spot)
    ok = kinematic_mask(xg, q2g, sampler.s)
    f2 = kern.nf2.f2a(xg, q2g) / kern.ion.A
    f1 = kern.nf2.f1a(xg, q2g) / kern.ion.A
    y = y_from_xq2(xg, q2g, sampler.s)
    curves = []
    for dfunc, label in (
            (model, "moment-constrained (A)"),
            (dm.make("moment_B", variant=args.variant,
                     dilution=args.dilution), "no-$F_1$ (B)"),
            (dm.make("toy", scale=1e-3), r"flat $\Delta/F_1=10^{-3}$")):
        amp = a_cos2phi(dfunc(xg, q2g, f1), f1, f2, xg, y)
        curves.append((xg[ok], np.asarray(amp)[ok], label))

    info = {"model": model.info(), "q2_ref": q2_ref,
            "backend": m5.describe_backends(args, backends),
            "q2_slice": float(q2_spot), "pzz": float(plan.pzz_true),
            "label": config.label()}
    return args, info, results, pts, curves


def draw(args, info, results, pts, curves):
    fig = fs.figure(fs.DOUBLE_W, 2.85)
    gs = GridSpec(1, 3, figure=fig, width_ratios=(1.0, 1.0, 1.32),
                  wspace=0.34, left=0.055, right=0.995, bottom=0.165,
                  top=0.80)

    shown = [r for r in results if r["spot"] in SHOWN]
    for col, m in enumerate(shown):
        ax = fig.add_subplot(gs[0, col])
        centers = 0.5 * (m["edges"][:-1] + m["edges"][1:])
        nbar = m["counts"].mean()
        mod = 1e3 * (m["counts"] / nbar - 1.0)
        mod_err = 1e3 * np.sqrt(np.maximum(m["counts"], 1.0)) / nbar
        phi = np.linspace(0, 2 * np.pi, 200)
        ax.errorbar(centers, mod, yerr=mod_err, fmt="o", color=fs.DATA,
                    ms=2.6, capsize=1.4, lw=0.8, zorder=3)
        ax.plot(phi, 1e3 * m["a2_eff"] * np.cos(2 * phi), "-",
                color=fs.TRUTH, lw=1.3, zorder=4)
        ax.plot(phi, 1e3 * m["amp"] * info["pzz"] * np.cos(2 * phi), "--",
                color=fs.FIT, lw=1.1, zorder=5)
        ax.set_xlim(0, 2 * np.pi)
        ax.set_xticks([0, np.pi, 2 * np.pi])
        ax.set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
        ax.set_xlabel(r"$\phi'=\phi-\phi_S$", labelpad=1)
        ax.axhline(0.0, color="0.85", lw=0.5, zorder=0)
        ax.tick_params(pad=2)
        if col == 0:
            ax.set_ylabel(r"$N(\phi')/\langle N\rangle-1$  $[10^{-3}]$",
                          labelpad=1)
        ax.set_title((r"$x=%.3g$, $Q^2=%.3g$ GeV$^2$" "\n"
                      r"$\hat A=(%.2f\pm%.2f)\times10^{-3}$")
                     % (m["x"], m["q2"], 1e3 * m["amp"], 1e3 * m["err"]),
                     pad=3)
        ax.margins(y=0.13)
        fs.panel_label(ax, "(%s)" % "ab"[col])

    # --- amplitude versus x ----------------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    xoff = 1.045
    ax.errorbar([p["x"] for p in pts], [1e3 * p["amp"] for p in pts],
                yerr=[1e3 * p["err"] for p in pts], fmt="s", mfc="none",
                color=fs.FIT, ms=3.0, capsize=1.4, lw=0.8, zorder=3,
                label=r"1 yr, %g fb$^{-1}$/u" % args.lumi_1yr)
    ax.errorbar([p["x"] * xoff for p in pts],
                [1e3 * p["amp10"] for p in pts],
                yerr=[1e3 * p["err10"] for p in pts], fmt="o",
                color=fs.DATA, ms=3.0, capsize=1.4, lw=0.8, zorder=4,
                label=r"10 yr, %g fb$^{-1}$/u" % args.lumi_10yr)
    for (xc, amp, label), color, ls in zip(
            curves, (fs.TRUTH, fs.ALT, fs.GREY), ("-", "-", "--")):
        ax.plot(xc, 1e3 * amp, ls, color=color, lw=1.2, label=label)
    ax.set_xscale("log")
    ax.set_xlabel(r"$x$", labelpad=1)
    ax.set_ylabel(r"$A^{\cos2\phi}$  $[10^{-3}]$", labelpad=1)
    ax.axhline(0.0, color="0.85", lw=0.5, zorder=0)
    ax.set_title((r"amplitude vs $x$" "\n"
                  r"$Q^2=%.3g$ GeV$^2$ slice") % info["q2_slice"], pad=3)
    ax.tick_params(pad=2)
    ax.margins(y=0.16)
    ax.legend(loc="upper left", fontsize=fs.FONT_SMALL,
              bbox_to_anchor=(0.10, 1.0))
    fs.panel_label(ax, "(c)")
    return fig


def caption(args, info, results, pts):
    """Figure 2's caption, from the same measurements the panels draw."""
    shown = [r for r in results if r["spot"] in SHOWN]
    a, b = shown[0], shown[1]
    return (
        r"Inclusive $\cos 2\phi'$ projections at %s with $P_{zz} = %.2f$. "
        r"(a), (b): the $\phi'$ distribution of the two super-bins that "
        r"span the accepted phase space, $x = %.3g$, $Q^{2} = %.3g$\,GeV$^2$ "
        r"and $x = %.3g$, $Q^{2} = %.3g$\,GeV$^2$, as one-year pseudo-data "
        r"(%g\,fb$^{-1}$ per nucleon) with statistical errors, against the "
        r"injected modulation (solid) and the binned fit (dashed); the fits "
        r"return $\hat A = (%s \pm %s)\times10^{-3}$ and "
        r"$(%s \pm %s)\times10^{-3}$ per unit $P_{zz}$. "
        r"(c): the amplitude extracted in %d merged $x$ bins along the "
        r"$Q^{2} = %.3g$\,GeV$^{2}$ slice, with independent one-year (open) "
        r"and ten-year (filled, %g\,fb$^{-1}$ per nucleon) "
        r"pseudo-measurements, against the moment-constrained $\Delta$ that "
        r"was injected, the no-$F_{1}$ variant and a flat "
        r"$\Delta/F_{1} = 10^{-3}$ reference. All four super-bins, "
        r"including the two not drawn here, are tabulated in "
        r"Table~\ref{tab:projections}. Statistical errors only; the "
        r"$\Delta$ model carries the $^{6}$Li per-nucleon dilution "
        r"$1/%.0f$."
        % (src.beam_math(info["label"]), info["pzz"],
           a["x"], a["q2"], b["x"], b["q2"], args.lumi_1yr,
           src.dec(1e3 * a["amp"], 2), src.dec(1e3 * a["err"], 2),
           src.dec(1e3 * b["amp"], 2), src.dec(1e3 * b["err"], 2),
           len(pts), info["q2_slice"], args.lumi_10yr,
           1.0 / args.dilution))


def main():
    args, info, results, pts, curves = measure_money5()
    with fs.rc():
        fig = draw(args, info, results, pts, curves)
        fs.save(fig, "fig2")
        plt.close(fig)
    src.write_caption("fig2", "figtwocaption", caption(args, info, results, pts))

    payload = {
        "figure": "fig2",
        "source_script": "evgen/scripts/money_cos2phi.py",
        "published_figure": "evgen/money_cos2phi_6Li.png",
        "settings": src.describe(args),
        "beam": info["label"], "pzz": info["pzz"],
        "delta_model": info["model"], "q2_ref": info["q2_ref"],
        "backend": info["backend"], "q2_slice": info["q2_slice"],
        "panels_shown": list(SHOWN),
        "spots": [{"spot": m["spot"], "x": m["x"], "q2": m["q2"],
                   "x_lo": m["x_lo"], "x_hi": m["x_hi"],
                   "q2_lo": m["q2_lo"], "q2_hi": m["q2_hi"],
                   "n_1yr": float(m["n"]),
                   "a_truth": float(m["truth"]),
                   "a_hat": float(m["amp"]),
                   "err_1yr": float(m["err"]),
                   "err_10yr": float(m["err10"]),
                   "significance_1yr": float(abs(m["truth"]) / m["err"]),
                   "pull": float((m["amp"] - m["truth"]) / m["err"])}
                  for m in results],
        "amplitude_vs_x": pts,
    }
    src.dump("fig2", payload)

    print("beam: %s, P_zz = %.2f, %s" % (info["label"], info["pzz"],
                                         info["backend"]))
    print("delta model: %s (<Q2> = %.3g GeV^2)"
          % (info["model"], info["q2_ref"]))
    for m in results:
        print("spot %d: x=%.3g Q2=%.3g  N_1yr=%.2e  A_truth=%+.2e  "
              "A_hat=%+.2e +- %.1e (1yr) +- %.1e (10yr)  sig=%.0f"
              % (m["spot"], m["x"], m["q2"], m["n"], m["truth"], m["amp"],
                 m["err"], m["err10"], abs(m["truth"]) / m["err"]))
    print("amplitude vs x: %d points on the Q2 = %.3g GeV^2 slice"
          % (len(pts), info["q2_slice"]))


if __name__ == "__main__":
    main()
