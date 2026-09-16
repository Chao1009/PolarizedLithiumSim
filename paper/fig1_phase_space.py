#!/usr/bin/env python3
"""Letter Figure 1 -- the (x, Q2) phase space and the analysis binning,
condensed from the phase-space companion
(`evgen/scripts/phase_space_bins.py`, `evgen/phase_space_bins_6Li.png`).

plans/07 WP6 asks for "phase space -> Fig. 1 (single inclusive panel with
bins, coherent support as contour or inset)".  The published figure is two
maps side by side; this one is the inclusive map alone, with the tagged
coherent channel folded in as contours of its own rate map over the same
grid -- the right-hand panel of the published figure reduced to the two lines
that say where the coherent sample lives.

Everything is `phase_space_bins`': its rate projection, its sweet-spot
selector, its Delta-extraction comb (kept by the exact money-plot selection,
including the pseudo-measurements that selection makes, which is why the
driver carries the same rng), its coherent projection and its published
settings, seed 20260817 among them.

Usage:  python3 paper/fig1_phase_space.py        (writes figs/fig1.*)
"""

import numpy as np

import _sources as src
import figstyle as fs

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
from matplotlib.colors import LogNorm                            # noqa: E402
from matplotlib.lines import Line2D                              # noqa: E402
from matplotlib.patches import Rectangle                         # noqa: E402

import money_cos2phi as m5                                       # noqa: E402
import phase_space_bins as psb                                   # noqa: E402
from money_cos2phi_coherent import best_superbin                 # noqa: E402

from polligen import bookkeeping as bk                           # noqa: E402
from polligen import coherent as coh                             # noqa: E402
from polligen.sample import InclusiveSampler                     # noqa: E402
from polligen.xsec import InclusiveKernel                        # noqa: E402

from polli_fastsim import beams, fom                             # noqa: E402
from polli_fastsim.farforward import HIGH_ACCEPTANCE             # noqa: E402
from polli_fastsim.kinematics import kinematic_mask              # noqa: E402
from polli_fastsim.polarized import toy_b1                       # noqa: E402

# Contours of the tagged coherent rate map, in events per bin per year.
COH_LEVELS = (1.0e4, 1.0e5)


def project():
    """Re-run the phase-space companion at its published defaults."""
    args = src.published_args(psb.main)

    config = beams.default_configs("6Li")[args.config]
    lumi1_pb = args.lumi_1yr * 1e3
    rng = np.random.default_rng(args.seed)

    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    grid = beams.analysis_grid(args.binning)
    backends = m5.pdf_backends(args, config.ion)
    model, _q2_ref = m5.build_delta_model(args, config, scenario,
                                          backends=backends)
    kern = InclusiveKernel(config.ion, b1_func=toy_b1, delta_func=model,
                           f2_source=backends["base"],
                           g1_model=backends["g1"],
                           nuclear_f2=backends["nuclear"],
                           r_func=backends["nuclear"].r_func)

    proj = fom.project_rates(config, scenario,
                             nuclear_f2=backends["nuclear"], **grid)
    obs = fom.project_observables(config, scenario, proj,
                                  kern.g1_model, toy_b1, model)
    spots = m5.pick_sweet_spots_banded(proj, obs["sig_a_cos2phi"])[:4]
    q2_slices = sorted({round(q2, 3) for _x, q2, _i, _j in spots})

    sampler = InclusiveSampler(kern, config, scenario, nx=60, nq2=45)
    plan = bk.transverse_tensor_plan(args.pzz)
    cat = plan.categories[0]

    xe = proj.x_edges
    slice_bins = {}
    for q2s in q2_slices:
        kept = []
        for i0 in range(0, xe.size - 2, 2):
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
            kept.append((float(xe[i0]), float(xe[i0 + 2])))
        slice_bins[q2s] = kept

    sc = coh.CoherentScenario()
    proj_c, n_coh, tagged = coh.project_coherent(
        config, scenario, sc, optics_list=(HIGH_ACCEPTANCE,),
        nuclear_f2=backends["nuclear"], **grid)
    n_tag = tagged[HIGH_ACCEPTANCE.name]
    sel, txlo, txhi, tq2lo, tq2hi = best_superbin(proj_c, n_tag)

    info = {
        "label": config.label(), "pzz": float(plan.pzz_true),
        "backend": m5.describe_backends(args, backends),
        "binning": args.binning, "nx": grid["nx"], "nq2": grid["nq2"],
        "n_dis": float(proj.n_events.sum()),
        "n_coh": float(n_coh.sum()), "n_tag": float(n_tag.sum()),
        "tag_acceptance": float(sc.tag_acceptance(
            HIGH_ACCEPTANCE.pt_cut_near_beam)),
        "pt_cut": float(HIGH_ACCEPTANCE.pt_cut_near_beam),
        "f0": float(sc.f0), "x_coh": float(sc.x_coh),
        "tagged_bin": {"x_lo": float(txlo), "x_hi": float(txhi),
                       "q2_lo": float(tq2lo), "q2_hi": float(tq2hi),
                       "n_1yr": float(n_tag[sel].sum())},
        "scenario": scenario,
        "cuts": {"eta_max": float(scenario.eta_max),
                 "e_prime_min": float(scenario.e_prime_min),
                 "y_min": float(scenario.y_min),
                 "y_max": float(scenario.y_max),
                 "w2_min": float(scenario.w2_min),
                 "lumi_1yr": float(scenario.lumi_fb_per_nucleon)},
    }
    spot_boxes = []
    for k, (xs, qs, i, j) in enumerate(spots):
        xlo, xhi, q2lo, q2hi = m5.superbin_edges(proj, i, j)
        spot_boxes.append({"spot": k + 1, "x": float(xs), "q2": float(qs),
                           "x_lo": float(xlo), "x_hi": float(xhi),
                           "q2_lo": float(q2lo), "q2_hi": float(q2hi)})
    return args, info, proj, proj_c, n_tag, spot_boxes, slice_bins, sampler


def draw(info, proj, proj_c, n_tag, spot_boxes, slice_bins, sampler):
    fig = fs.figure(fs.COLUMN_W, 3.0)
    ax = fig.add_axes((0.155, 0.135, 0.665, 0.845))
    cax = fig.add_axes((0.845, 0.135, 0.035, 0.845))

    vals = np.ma.masked_where(~proj.accepted | (proj.n_events < 1.0),
                              proj.n_events)
    pcm = ax.pcolormesh(proj.x_edges, proj.q2_edges, vals.T,
                        norm=LogNorm(vmin=1.0), cmap="viridis",
                        shading="auto", rasterized=True)
    cb = fig.colorbar(pcm, cax=cax)
    cb.set_label("DIS events / bin (1 yr)", fontsize=fs.FONT_LARGE,
                 labelpad=2)
    cb.ax.tick_params(labelsize=fs.FONT_SMALL, pad=1.5)

    # the Delta-extraction x-bin combs, one row per Q2 slice
    for q2s, kept in slice_bins.items():
        for xlo, xhi in kept:
            ax.add_patch(Rectangle((xlo, q2s / 1.6), xhi - xlo,
                                   q2s * 1.6 - q2s / 1.6, fill=False,
                                   edgecolor="0.12", linewidth=0.35,
                                   zorder=5))

    # the four phi' super-bins, cased in white so they read on any level
    for box in spot_boxes:
        for color, lw, dz in (("white", 1.9, 0), (fs.FIT, 0.9, 1)):
            ax.add_patch(Rectangle(
                (box["x_lo"], box["q2_lo"]),
                box["x_hi"] - box["x_lo"], box["q2_hi"] - box["q2_lo"],
                fill=False, edgecolor=color, linewidth=lw, zorder=6 + dz,
                joinstyle="miter"))

    # the coherent support, as contours of the tagged rate map
    cs = ax.contour(proj_c.x, proj_c.q2, n_tag, levels=list(COH_LEVELS),
                    colors=[fs.ALT], linewidths=(0.7, 1.2), zorder=8)
    tb = info["tagged_bin"]
    for color, lw, dz in (("white", 1.9, 0), (fs.ALT, 0.9, 1)):
        ax.add_patch(Rectangle(
            (tb["x_lo"], tb["q2_lo"]), tb["x_hi"] - tb["x_lo"],
            tb["q2_hi"] - tb["q2_lo"], fill=False, edgecolor=color,
            linewidth=lw, zorder=8 + dz, joinstyle="miter"))

    # The y and W^2 cut boundaries, drawn by the producing script's own
    # helper.  Its slope labels are placed and sized (6.5 pt) for a 13 in
    # canvas and collide with each other at one column, so they are hidden and
    # the caption names the cuts instead; the lines themselves are the
    # helper's.  Nothing is redrawn by hand.
    before = set(id(o) for o in ax.findobj(plt.Text))
    psb.draw_guides(ax, sampler.s, info["scenario"])
    n_hidden = 0
    for obj in ax.findobj(plt.Text):
        if id(obj) not in before and obj.get_text().strip():
            obj.set_visible(False)
            n_hidden += 1

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(proj.x_edges[0], proj.x_edges[-1])
    ax.set_ylim(proj.q2_edges[0], proj.q2_edges[-1])
    ax.set_xlabel(r"$x$", labelpad=1)
    ax.set_ylabel(r"$Q^2$  [GeV$^2$]", labelpad=1)
    ax.tick_params(pad=2)
    ax.legend(handles=[
        Line2D([], [], color=fs.FIT, lw=1.1,
               label=r"$\phi'$ super-bins"),
        Line2D([], [], color="0.12", lw=0.7,
               label=r"$\Delta$-extraction bins"),
        Line2D([], [], color=fs.ALT, lw=1.1,
               label=r"coherent tag, $10^{4}$, $10^{5}$/bin"),
        Line2D([], [], color="0.55", lw=0.8, ls="--",
               label=r"$y$, $W^2$ cut boundaries"),
    ], loc="upper left", fontsize=fs.FONT_SMALL, borderpad=0.3,
        handlelength=1.2, labelspacing=0.22, framealpha=0.9)
    del cs, n_hidden
    return fig


def caption(info):
    """Figure 1's caption, from the same arrays the panel draws."""
    c = info["cuts"]
    return (
        r"The $(x, Q^{2})$ phase space of the measurement at %s for a "
        r"one-year programme (%g\,fb$^{-1}$ per nucleon). Colour: expected "
        r"inclusive DIS events per bin of the $%d \times %d$ log--log "
        r"analysis grid, drawn where at least one event is expected, after "
        r"the scattered-electron cuts $|\eta| \le %.1f$, "
        r"$E' \ge %.1f$\,GeV, $%.2f \le y \le %.2f$ and "
        r"$W^{2} \ge %.0f$\,GeV$^{2}$ (dashed: the $y$ and $W^{2}$ "
        r"boundaries). Vermillion: the four $\cos 2\phi'$ super-bins of "
        r"Table~\ref{tab:projections}; thin dark lines: the merged $x$ bins "
        r"of the $\Delta$ extraction along the sweet-spot $Q^{2}$ slices. "
        r"Green: the coherent channel over the same grid, as contours of the "
        r"tagged rate at $10^{%d}$ and $10^{%d}$ recoils per bin per year "
        r"for a scenario coherent fraction $f_{0} = %.2f$ dying at "
        r"$x_{\rm coh} = %.2f$ and a near-beam $p_{T} > %.2f$\,GeV tag of "
        r"acceptance $%.1f\%%$, with the box of its highest-yield "
        r"super-bin. One year gives $%s$ inclusive events, $%s$ coherent and "
        r"$%s$ tagged."
        % (src.beam_math(info["label"]), c["lumi_1yr"],
           info["nx"], info["nq2"], c["eta_max"], c["e_prime_min"],
           c["y_min"], c["y_max"], c["w2_min"],
           int(round(np.log10(COH_LEVELS[0]))),
           int(round(np.log10(COH_LEVELS[1]))),
           info["f0"], info["x_coh"], info["pt_cut"],
           100 * info["tag_acceptance"],
           src.sci(info["n_dis"], 2), src.sci(info["n_coh"], 2),
           src.sci(info["n_tag"], 2)))


def main():
    (args, info, proj, proj_c, n_tag, spot_boxes, slice_bins,
     sampler) = project()
    with fs.rc():
        fig = draw(info, proj, proj_c, n_tag, spot_boxes, slice_bins, sampler)
        fs.save(fig, "fig1")
        plt.close(fig)
    src.write_caption("fig1", "figonecaption", caption(info))

    payload = {
        "figure": "fig1",
        "source_script": "evgen/scripts/phase_space_bins.py",
        "published_figure": "evgen/phase_space_bins_6Li.png",
        "settings": src.describe(args),
        "beam": info["label"], "pzz": info["pzz"],
        "backend": info["backend"],
        "binning": "%s, %d x %d" % (info["binning"], info["nx"], info["nq2"]),
        "n_dis_1yr": info["n_dis"], "n_coh_1yr": info["n_coh"],
        "n_tag_1yr": info["n_tag"],
        "tag_acceptance": info["tag_acceptance"], "pt_cut": info["pt_cut"],
        "f0": info["f0"], "x_coh": info["x_coh"],
        "coherent_contours": list(COH_LEVELS),
        "superbins": spot_boxes,
        "tagged_bin": info["tagged_bin"],
        "delta_bins": {("%.3g" % q2s): kept
                       for q2s, kept in slice_bins.items()},
    }
    src.dump("fig1", payload)

    print("beam: %s, %s" % (info["label"], info["backend"]))
    print("binning: %s, %d x %d" % (info["binning"], info["nx"], info["nq2"]))
    print("N_DIS = %.3e   N_coh = %.3e   N_tag = %.3e (acc %.1f%%, "
          "p_T > %.2f GeV, f0 = %.2f)"
          % (info["n_dis"], info["n_coh"], info["n_tag"],
             100 * info["tag_acceptance"], info["pt_cut"], info["f0"]))
    for b in spot_boxes:
        print("super-bin %d: x [%.3g, %.3g]  Q2 [%.3g, %.3g]"
              % (b["spot"], b["x_lo"], b["x_hi"], b["q2_lo"], b["q2_hi"]))
    for q2s, kept in slice_bins.items():
        print("slice Q2=%.3g: %d x bins, x [%.3g, %.3g]"
              % (q2s, len(kept), kept[0][0], kept[-1][1]))
    tb = info["tagged_bin"]
    print("tagged super-bin: x [%.3g, %.3g]  Q2 [%.3g, %.3g]  N=%.3g"
          % (tb["x_lo"], tb["x_hi"], tb["q2_lo"], tb["q2_hi"], tb["n_1yr"]))


if __name__ == "__main__":
    main()
