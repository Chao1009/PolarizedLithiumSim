#!/usr/bin/env python3
"""Letter Figure 4 -- the coherent intact-6Li channel, condensed from money
plot 6 (`evgen/scripts/money_cos2phi_coherent.py`,
`evgen/money_cos2phi_coherent_6Li.png`) with the WP5 optics curve folded in.

plans/07 WP6 asks for "money plot 6 -> Fig. 4 (a2 anchor/band + tagged phi',
acceptance inset)", and the run-19 ledger's T24 asks for the WP5
coherent-optics curve to appear in this figure with its assumption named in
the caption.  Both are here:

  (a) the deformation coefficient a_2(|t|) -- the digitized polarized-deuteron
      anchor per m state, the scaled 6Li ensemble band over
      eps_B0 in -(0.04-0.13), the analysis window and the form-factor zero --
      with an inset carrying the tag acceptance versus the near-beam
      envelope, the curve of `evgen/scripts/coherent_optics_scan.py` panel
      (a).  The inset's slot curve is an ASSUMED cutout: a rectangle 2.5
      times as wide as it is high, the ePIC horizontal slot, which is that
      script's `--cut-scale-x` default; the circular curve beside it is the
      exp(-B p_T^2) the rate model itself uses, and the marker is the
      published 0.20 GeV working point.
  (b) the phi' modulation of the tagged sample at the deformation-anchored
      <a_2>, with the flat gluon-transversity scenario for scale.

Every number is `money_cos2phi_coherent`'s: its scenario object, its
best-super-bin selector, its pseudo-histogram, its binned fit and its error
formula, at its published defaults (seed 20260810).  The inset is
`coherent_optics_scan.acceptance_curve` at that script's published defaults.

Usage:  python3 paper/fig4_coherent.py           (writes figs/fig4.*)
"""

import numpy as np

import _sources as src
import figstyle as fs

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402

import coherent_optics_scan as cos_scan                          # noqa: E402
import money_cos2phi as m5                                       # noqa: E402
import money_cos2phi_coherent as m6                              # noqa: E402

from polligen import coherent as coh                             # noqa: E402
from polligen.estimators import (cos2phi_fit_binned,             # noqa: E402
                                 cos2phi_fit_err)
from polligen.sample import phi_histogram_pseudo                 # noqa: E402

from polli_fastsim import beams, fom                             # noqa: E402
from polli_fastsim.farforward import HIGH_ACCEPTANCE             # noqa: E402

# The envelope grid of coherent_optics_scan panel (a), unchanged.
CUTS = np.linspace(0.05, 0.70, 140)


def measure_money6():
    """Re-run money plot 6's measurement, in its order, at its defaults."""
    args = src.published_args(m6.main)

    config = beams.default_configs("6Li")[args.config]
    rng = np.random.default_rng(args.seed)
    sc = coh.CoherentScenario(amp=args.amp, eps_b0=args.eps_b0)
    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    lumi_ratio = args.lumi_10yr / args.lumi_1yr
    backends = m5.pdf_backends(args, config.ion)
    proj, n_coh, tagged = coh.project_coherent(
        config, scenario, sc, optics_list=(HIGH_ACCEPTANCE,),
        nuclear_f2=backends["nuclear"])
    n_tag = tagged[HIGH_ACCEPTANCE.name]
    cut_ha = HIGH_ACCEPTANCE.pt_cut_near_beam

    sel, xlo, xhi, q2lo, q2hi = m6.best_superbin(proj, n_tag)
    n_bin = float(n_tag[sel].sum())
    a2_def = sc.a2_tagged(cut_ha, args.pzz)
    a2_band = [coh.CoherentScenario(eps_b0=e).a2_tagged(cut_ha, args.pzz)
               for e in m6.EPS_BAND]

    counts, edges = phi_histogram_pseudo(n_bin, a2_def, nbins=24, rng=rng)
    amp_hat = cos2phi_fit_binned(counts, edges, args.pzz)
    err = cos2phi_fit_err(n_bin, args.pzz)
    err10 = cos2phi_fit_err(lumi_ratio * n_bin, args.pzz)

    # panel (a)'s curves, all from polligen.coherent
    tt = np.linspace(0.0, 0.25, 100)
    td = np.array(sorted(coh.MANTYSAARI_A2_DEUTERON))
    result = {
        "args": args, "label": config.label(),
        "n_coh": float(n_coh.sum()), "n_tag": float(n_tag.sum()),
        "bin": {"x_lo": float(xlo), "x_hi": float(xhi),
                "q2_lo": float(q2lo), "q2_hi": float(q2hi),
                "n_1yr": n_bin},
        "a2_def": float(a2_def), "a2_band": [float(v) for v in a2_band],
        "a2_per_pzz": float(a2_def / args.pzz),
        "amp_hat": float(amp_hat), "err": float(err), "err10": float(err10),
        "counts": counts, "edges": edges,
        "gluonic": float(args.amp * args.pzz),
        "slope_b": float(sc.slope_b), "f0": float(sc.f0),
        "pt_cut": float(cut_ha),
        "tag_acceptance": float(sc.tag_acceptance(cut_ha)),
        "t_curve": tt,
        "a2_mid": np.asarray(sc.a2_deformation(tt, args.pzz), dtype=float),
        "a2_lo": np.asarray(coh.CoherentScenario(
            eps_b0=m6.EPS_BAND[0]).a2_deformation(tt, args.pzz), dtype=float),
        "a2_hi": np.asarray(coh.CoherentScenario(
            eps_b0=m6.EPS_BAND[1]).a2_deformation(tt, args.pzz), dtype=float),
        "anchor_t": td,
        "anchor_m0": np.array([coh.MANTYSAARI_A2_DEUTERON[k][0] for k in td]),
        "anchor_m1": np.array([coh.MANTYSAARI_A2_DEUTERON[k][1] for k in td]),
        "t_window": m6.T_WINDOW, "t_dip": m6.T_DIP,
        "eps_band": m6.EPS_BAND,
    }
    return result


def optics_curve():
    """The WP5 acceptance curve (T24), at coherent_optics_scan's defaults."""
    a = src.published_args(cos_scan.main)
    slot = {}
    for b in (40.0, 50.0, 60.0):
        acc, _a2 = cos_scan.acceptance_curve(CUTS, b, a.cut_scale_x, 1.0)
        slot[b] = acc
    circular = np.exp(-50.0 * CUTS ** 2)
    return a, slot, circular


def draw(r, scan_args, slot, circular):
    fig, (ax1, ax2) = fs.subplots(1, 2, width=fs.DOUBLE_W, height=2.75)
    args = r["args"]

    # --- (a) a2(|t|): anchor, band, window -------------------------------
    ax1.fill_between(r["t_curve"], r["a2_lo"], r["a2_hi"], color=fs.TRUTH,
                     alpha=0.20, lw=0,
                     label=r"$^6$Li, $\varepsilon_{B0}\in-(0.04$–$0.13)$")
    ax1.plot(r["t_curve"], r["a2_mid"], "-", color=fs.TRUTH, lw=1.4)
    ax1.plot(r["anchor_t"], r["anchor_m0"], "s", color="0.35", ms=3.4,
             label=r"$d$, $m=0$ (digitized)")
    ax1.plot(r["anchor_t"], r["anchor_m1"], "^", color="0.62", ms=3.4,
             label=r"$d$, $m=\pm1$")
    ax1.axvspan(r["t_window"][0], r["t_window"][1], color=fs.ALT, alpha=0.08,
                lw=0, zorder=0)
    ax1.axvline(r["t_dip"], color="0.3", ls=":", lw=0.9, zorder=1)
    ax1.axhline(0.0, color="0.85", lw=0.5, zorder=0)
    ax1.set_xlim(0.0, 0.33)
    ax1.set_ylim(-0.34, 0.72)
    ax1.set_xlabel(r"$|t|$  [GeV$^2$]", labelpad=1)
    ax1.set_ylabel(r"$a_2(|t|)$", labelpad=1)
    ax1.tick_params(pad=2)
    ax1.legend(loc="lower left", fontsize=fs.FONT_SMALL, borderpad=0.3,
               handlelength=1.2, labelspacing=0.22)
    fs.panel_label(ax1, "(a)")

    # --- the T24 inset: tag acceptance vs the near-beam envelope ---------
    # placed over the empty upper-left of the panel, clear of the digitized
    # anchor points at |t| = 0.20 and 0.30, which the letter must not hide
    ins = ax1.inset_axes((0.150, 0.575, 0.430, 0.385))
    for b, ls, lw in ((40.0, ":", 0.7), (50.0, "-", 1.2), (60.0, "--", 0.7)):
        ins.plot(CUTS, slot[b], ls, color=fs.TRUTH, lw=lw)
    ins.plot(CUTS, circular, "-", color=fs.GREY, lw=0.9)
    ins.plot([r["pt_cut"]], [r["tag_acceptance"]], "o", color=fs.FIT,
             ms=3.4, zorder=5)
    ins.set_yscale("log")
    ins.set_ylim(1e-3, 1.6)
    ins.set_xlim(CUTS[0], CUTS[-1])
    ins.set_xticks([0.2, 0.6])
    ins.set_yticks([1e-3, 1e-1])
    ins.tick_params(labelsize=fs.FONT_SMALL, pad=1.2, length=2.0)
    ins.set_xlabel("envelope [GeV]", fontsize=fs.FONT_SMALL, labelpad=0.5)
    ins.annotate("tagged fraction", xy=(0.96, 0.90),
                 xycoords="axes fraction", ha="right", va="top",
                 fontsize=fs.FONT_SMALL)

    # --- (b) the tagged phi' modulation ----------------------------------
    centers = 0.5 * (r["edges"][:-1] + r["edges"][1:])
    nbar = r["counts"].mean()
    phi = np.linspace(0, 2 * np.pi, 200)
    ax2.errorbar(centers, r["counts"] / nbar,
                 yerr=np.sqrt(np.maximum(r["counts"], 1.0)) / nbar, fmt="o",
                 color=fs.DATA, ms=2.6, capsize=1.4, lw=0.8, zorder=3)
    ax2.plot(phi, 1 + r["a2_def"] * np.cos(2 * phi), "-", color=fs.TRUTH,
             lw=1.3, zorder=4,
             label=r"$\langle a_2\rangle=%.3f$ (deformation)" % r["a2_def"])
    ax2.plot(phi, 1 + r["amp_hat"] * args.pzz * np.cos(2 * phi), "--",
             color=fs.FIT, lw=1.1, zorder=5, label="binned fit")
    ax2.plot(phi, 1 + r["gluonic"] * np.cos(2 * phi), "-", color=fs.GREY,
             lw=0.9, zorder=2,
             label=r"gluon transversity (%.3f)" % r["gluonic"])
    ax2.axhline(1.0, color="0.85", lw=0.5, zorder=0)
    ax2.set_xlim(0, 2 * np.pi)
    ax2.set_xticks([0, np.pi, 2 * np.pi])
    ax2.set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
    ax2.set_xlabel(r"$\phi'=\phi-\phi_S$", labelpad=1)
    ax2.set_ylabel(r"$N(\phi')/\langle N\rangle$", labelpad=1)
    ax2.tick_params(pad=2)
    ax2.set_ylim(0.955, 1.062)
    ax2.legend(loc="upper center", fontsize=fs.FONT_SMALL, borderpad=0.3,
               handlelength=1.2, labelspacing=0.22, ncol=1)
    ax2.set_title(r"tagged sample, 1 yr: $\hat A=(%.1f\pm%.1f)\times10^{-3}$"
                  r" per unit $P_{zz}$"
                  % (1e3 * r["amp_hat"], 1e3 * r["err"]), pad=3)
    fs.panel_label(ax2, "(b)")

    fig.subplots_adjust(left=0.075, right=0.995, bottom=0.165, top=0.905,
                        wspace=0.26)
    return fig


def caption(r, scan_args):
    """Figure 4's caption, from the same scenario the panels draw.

    T24 of the run-19 ledger asks for the WP5 optics curve to appear here
    "with the assumption named in its caption"; the inset sentence below is
    that naming, and it states the cutout aspect ratio, the slope band and
    the circular reference the rate model itself uses.
    """
    args = r["args"]
    b = r["bin"]
    return (
        r"The coherent channel $e\,^{6}\mathrm{Li} \to e' X\,"
        r"^{6}\mathrm{Li}(\mathrm{g.s.})$ at %s with $P_{zz} = %.2f$. "
        r"(a) The deformation coefficient $a_{2}(|t|)$: the digitized "
        r"polarized-deuteron anchor per magnetic substate, against the "
        r"$^{6}$Li ensemble band obtained by scaling it with "
        r"$\varepsilon_{B0} \in -(%.2f\text{--}%.2f)$ (central value "
        r"$%.2f$), whose sign is opposite to the deuteron's because "
        r"$Q(^{6}\mathrm{Li}) < 0$. The shaded strip is the analysis window "
        r"$%.3f$--$%.2f$\,GeV$^{2}$ and the dotted line the form-factor "
        r"zero at $|t| = %.2f$\,GeV$^{2}$, above which the linear model "
        r"does not hold. Inset: the tagged fraction of coherent recoils "
        r"against the near-beam envelope. The assumption it rests on is the "
        r"cutout geometry: the three blue curves assume a rectangular cutout "
        r"$%.1f$ times as wide as it is high, an illustrative aspect ratio "
        r"rather than the measured one, at $B = 40$, $50$ and $60$\,GeV$^{-2}$ (dotted, solid, "
        r"dashed), the grey curve is the circular $e^{-B p_{T}^{2}}$ the "
        r"rate model itself uses at $B = %g$\,GeV$^{-2}$, and the marker is "
        r"the $p_{T} > %.2f$\,GeV working point at which that model gives "
        r"$%.1f\%%$. No $^{6}$Li optics delivering this envelope is "
        r"published, which is why the acceptance is shown as a curve and "
        r"not as a point. "
        r"(b) The $\phi'$ distribution of the tagged sample in the "
        r"highest-yield super-bin, $x \in [%s, %s]\times10^{-3}$, "
        r"$Q^{2} \in [%.2f, %.2f]$\,GeV$^{2}$, holding $%s$ recoils in one "
        r"year (%g\,fb$^{-1}$ per nucleon) for a coherent fraction "
        r"$f_{0} = %.2f$: pseudo-data with statistical errors against the "
        r"deformation-anchored $\langle a_{2}\rangle_{\rm tag} = %.3f$ at "
        r"$P_{zz} = %.2f$ (solid) and the binned fit (dashed), with the flat "
        r"gluon-transversity scenario (grey, $%.3f$) for scale. The fit "
        r"returns $\hat A = (%s \pm %s)\times10^{-3}$ per unit $P_{zz}$ "
        r"against $%.3f$ injected, i.e.\ $5\sigma$ floors of "
        r"$%s\times10^{-3}$ in one year and $%s\times10^{-3}$ in ten "
        r"(Table~\ref{tab:projections}). Statistical errors only; "
        r"backgrounds, purity and tensor radiative corrections are not "
        r"included."
        % (src.beam_math(r["label"]), args.pzz,
           abs(r["eps_band"][0]), abs(r["eps_band"][1]), args.eps_b0,
           r["t_window"][0], r["t_window"][1], r["t_dip"],
           scan_args.cut_scale_x, r["slope_b"], r["pt_cut"],
           100 * r["tag_acceptance"],
           src.dec(1e3 * b["x_lo"], 2), src.dec(1e3 * b["x_hi"], 2),
           b["q2_lo"], b["q2_hi"], src.sci(b["n_1yr"], 3),
           args.lumi_1yr, r["f0"], r["a2_def"], args.pzz, r["gluonic"],
           src.dec(1e3 * r["amp_hat"], 1), src.dec(1e3 * r["err"], 1),
           r["a2_per_pzz"],
           src.dec(1e3 * 5 * r["err"], 1),
           src.dec(1e3 * 5 * r["err10"], 1)))


def main():
    r = measure_money6()
    scan_args, slot, circular = optics_curve()
    with fs.rc():
        fig = draw(r, scan_args, slot, circular)
        fs.save(fig, "fig4")
        plt.close(fig)
    src.write_caption("fig4", "figfourcaption", caption(r, scan_args))

    args = r["args"]
    payload = {
        "figure": "fig4",
        "source_script": "evgen/scripts/money_cos2phi_coherent.py",
        "published_figure": "evgen/money_cos2phi_coherent_6Li.png",
        "inset_source_script": "evgen/scripts/coherent_optics_scan.py",
        "settings": src.describe(args),
        "inset_settings": src.describe(scan_args),
        "beam": r["label"], "pzz": float(args.pzz),
        "n_coh_1yr": r["n_coh"], "n_tag_1yr": r["n_tag"],
        "slope_b": r["slope_b"], "f0": r["f0"],
        "pt_cut": r["pt_cut"], "tag_acceptance": r["tag_acceptance"],
        "cutout_aspect_x": float(scan_args.cut_scale_x),
        "t_window": [float(v) for v in r["t_window"]],
        "t_dip": float(r["t_dip"]),
        "eps_b0": float(args.eps_b0),
        "eps_band": [float(v) for v in r["eps_band"]],
        "best_bin": r["bin"],
        "a2_tagged": r["a2_def"], "a2_tagged_band": r["a2_band"],
        "a2_per_pzz": r["a2_per_pzz"],
        "amp_hat": r["amp_hat"],
        "err_1yr": r["err"], "err_10yr": r["err10"],
        "floor5_1yr": 5.0 * r["err"], "floor5_10yr": 5.0 * r["err10"],
        "gluonic_scenario": r["gluonic"],
        "acceptance_curve": {
            "cuts": [float(c) for c in CUTS],
            "slot_b50": [float(v) for v in slot[50.0]],
            "slot_b40": [float(v) for v in slot[40.0]],
            "slot_b60": [float(v) for v in slot[60.0]],
            "circular_b50": [float(v) for v in circular],
        },
    }
    src.dump("fig4", payload)

    b = r["bin"]
    print("beam: %s, P_zz = %.2f" % (r["label"], args.pzz))
    print("coherent produced: %.3g   RP-tagged: %.3g   (acc %.1f%% at "
          "p_T > %.2f GeV)" % (r["n_coh"], r["n_tag"],
                               100 * r["tag_acceptance"], r["pt_cut"]))
    print("best super-bin: x [%.3g, %.3g], Q2 [%.3g, %.3g], N_1yr=%.3g"
          % (b["x_lo"], b["x_hi"], b["q2_lo"], b["q2_hi"], b["n_1yr"]))
    print("dA: %.4f (1 yr) / %.4f (10 yr); 5-sigma floors %.4f / %.4f"
          % (r["err"], r["err10"], 5 * r["err"], 5 * r["err10"]))
    print("deformation <a2>_tag = %.4f (band %.4f..%.4f) at P_zz = %.2f, "
          "i.e. %.4f per unit P_zz; gluonic scenario %.4f"
          % (r["a2_def"], r["a2_band"][0], r["a2_band"][1], args.pzz,
             r["a2_per_pzz"], r["gluonic"]))
    print("A_hat = %.4f +- %.4f (truth %.4f)"
          % (r["amp_hat"], r["err"], r["a2_per_pzz"]))
    print("inset: cutout %.1f x 1 slot, B = 40/50/60, circular reference; "
          "acceptance at %.2f GeV = %.3f"
          % (scan_args.cut_scale_x, r["pt_cut"], r["tag_acceptance"]))


if __name__ == "__main__":
    main()
