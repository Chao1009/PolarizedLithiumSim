#!/usr/bin/env python3
"""Money plot 6: cos 2phi modulation in COHERENT e+6Li with the intact
recoil tagged at IP6.

Process: e + 6Li -> e' + X + 6Li(g.s.) -- the nucleus survives in its 1+
ground state and is detected in the Roman-Pot near-beam pT tail (it keeps
the exact beam rigidity, R = 1.000, so the pT cut is the only handle;
plans/06 discusses the backgrounds and the Z-discrimination this tag
needs).

The tensor modulation carries two anchored mechanisms (plans/06 SS6.4b):
a t-linear DEFORMATION term scaled from the polarized-deuteron CGC
calculation of Mantysaari et al. (PLB 858:139053, arXiv:2408.13213 --
digitized in coherent.MANTYSAARI_A2_DEUTERON), predicted sign-flipped
for 6Li (Q < 0), and a flat GLUON-TRANSVERSITY scenario bounded by
lattice + Drell-Yan estimates.  Rates come from the DIS maps scaled by
the scenario coherent fraction.

Panels: (a) recoil |t| spectrum vs the RP pT cuts; (b) a_2(|t|) --
deuteron anchor points vs the scaled 6Li band, with the tagged window
and the form-factor-minimum exclusion; (c) tagged coherent yield vs x;
(d) phi'-modulation pseudo-data at the deformation-anchored amplitude
with the gluonic scenario overlaid.

Error bars assume the whole luminosity in the transverse-tensor fill
(same convention as the inclusive plot and the analytic FOM maps; a
run-plan split scales N down proportionally).

Usage:  python3 scripts/money_cos2phi_coherent.py
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from polligen import coherent as coh  # noqa: E402
from polligen.estimators import (cos2phi_fit_binned,  # noqa: E402
                                 cos2phi_fit_err)
from polligen.sample import phi_histogram_pseudo  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.farforward import (HIGH_ACCEPTANCE,  # noqa: E402
                                      HIGH_DIVERGENCE)

C_TRUTH, C_FIT, C_ALT = "#0072B2", "#D55E00", "#009E73"
EPS_BAND = (-0.04, -0.13)  # eps_b0 band for 6Li (plans/06 SS6.4b)
T_DIP = 0.31               # 6Li form-factor minimum [GeV^2]


def best_superbin(proj, n_tag, pad=1):
    """(mask over the map, x_lo/hi, q2_lo/hi) of the tagged-count-maximal
    FOM bin padded by `pad` bins on each side."""
    i, j = np.unravel_index(np.argmax(n_tag), n_tag.shape)
    xe, qe = proj.x_edges, proj.q2_edges
    ilo, ihi = max(i - pad, 0), min(i + pad + 1, xe.size - 1)
    jlo, jhi = max(j - pad, 0), min(j + pad + 1, qe.size - 1)
    sel = np.zeros_like(n_tag, dtype=bool)
    sel[ilo:ihi, jlo:jhi] = True
    return sel, xe[ilo], xe[ihi], qe[jlo], qe[jhi]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2))
    ap.add_argument("--lumi-1yr", type=float, default=10.0,
                    help="1-year EIC program [fb^-1/nucleon]")
    ap.add_argument("--lumi-10yr", type=float, default=100.0,
                    help="10-year EIC program [fb^-1/nucleon]")
    ap.add_argument("--pzz", type=float, default=0.60)
    ap.add_argument("--eps-b0", type=float, default=-0.08,
                    help="m=0 relative slope modulation of 6Li "
                         "(deformation mechanism; band -0.04..-0.13)")
    ap.add_argument("--amp", type=float, default=0.01,
                    help="flat gluon-transversity scenario amplitude "
                         "at Pzz=1 (band 3e-3..1e-2)")
    ap.add_argument("--seed", type=int, default=20260810)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs("6Li")[args.config]
    rng = np.random.default_rng(args.seed)
    sc = coh.CoherentScenario(amp=args.amp, eps_b0=args.eps_b0)
    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    lumi_ratio = args.lumi_10yr / args.lumi_1yr
    proj, n_coh, tagged = coh.project_coherent(
        config, scenario, sc, optics_list=(HIGH_ACCEPTANCE, HIGH_DIVERGENCE))
    n_tag = tagged[HIGH_ACCEPTANCE.name]
    cut_ha = HIGH_ACCEPTANCE.pt_cut_near_beam

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11.6, 8.4))

    # --- (a) recoil t spectrum vs the RP pT cuts --------------------------
    t = np.linspace(0.0, 0.25, 300)
    ax1.plot(t, sc.slope_b * np.exp(-sc.slope_b * t), "-", color=C_TRUTH,
             lw=1.8, label=r"$e^{-B|t|}$, $B=%g$ GeV$^{-2}$" % sc.slope_b)
    for b in (40.0, 60.0):
        ax1.plot(t, b * np.exp(-b * t), "-", color=C_TRUTH, lw=0.8,
                 alpha=0.45)
    ax1.fill_between(t, 1e-3, sc.slope_b * np.exp(-sc.slope_b * t),
                     where=t >= cut_ha ** 2, color=C_TRUTH, alpha=0.12)
    for optics, color in ((HIGH_ACCEPTANCE, C_FIT),
                          (HIGH_DIVERGENCE, "0.35")):
        cut2 = optics.pt_cut_near_beam ** 2
        acc = sc.tag_acceptance(optics.pt_cut_near_beam)
        band = [coh.CoherentScenario(slope_b=b).tag_acceptance(
            optics.pt_cut_near_beam) for b in (60.0, 40.0)]
        ax1.axvline(cut2, color=color, ls="--", lw=1.2)
        ax1.annotate(
            "%s\n$p_T>%.2f$ GeV\nacc = %.3g\n[%.2g, %.2g]"
            % ({"high-acceptance": "near-beam envelope", "high-divergence": "wide envelope"}.get(optics.name, optics.name),
               optics.pt_cut_near_beam, acc, *band),
            xy=(cut2, 2e-3), xytext=(4, 6), textcoords="offset points",
            fontsize=7, color=color, va="bottom")
    ax1.set_yscale("log")
    ax1.set_ylim(1e-3, 80)
    ax1.set_xlabel(r"$|t| \simeq p_T^2$  [GeV$^2$]")
    ax1.set_ylabel(r"$(1/\sigma)\, d\sigma_{\rm coh}/d|t|$  [GeV$^{-2}$]")
    ax1.set_title("intact-recoil $t$ spectrum vs near-beam $p_T$ cut", fontsize=9)
    ax1.legend(fontsize=7, loc="upper right")
    ax1.tick_params(labelsize=8)

    # --- (b) a_2(|t|): deuteron anchor vs scaled 6Li band -----------------
    td = np.array(sorted(coh.MANTYSAARI_A2_DEUTERON))
    a0 = np.array([coh.MANTYSAARI_A2_DEUTERON[k][0] for k in td])
    a1 = np.array([coh.MANTYSAARI_A2_DEUTERON[k][1] for k in td])
    ax2.plot(td, a0, "s", color="0.45", ms=5,
             label=r"$d$, $m=0$ (arXiv:2408.13213, digitized)")
    ax2.plot(td, a1, "^", color="0.65", ms=5, label=r"$d$, $m=\pm1$")
    tt = np.linspace(0.0, 0.25, 100)
    lo = coh.CoherentScenario(eps_b0=EPS_BAND[0]).a2_deformation(tt,
                                                                 args.pzz)
    hi = coh.CoherentScenario(eps_b0=EPS_BAND[1]).a2_deformation(tt,
                                                                 args.pzz)
    ax2.fill_between(tt, lo, hi, color=C_TRUTH, alpha=0.20,
                     label=r"$^6$Li ensemble, $\epsilon_{B0}\in"
                           r"[-0.04,-0.13]$, $P_{zz}=%.1f$" % args.pzz)
    ax2.plot(tt, sc.a2_deformation(tt, args.pzz), "-", color=C_TRUTH,
             lw=1.6)
    ax2.axvspan(cut_ha ** 2, 0.25, color=C_ALT, alpha=0.07)
    ax2.annotate("tagged window", xy=(0.145, 0.58),
                 xycoords=("data", "axes fraction"), fontsize=7,
                 color=C_ALT, ha="center")
    ax2.axvline(T_DIP, color="0.3", ls=":", lw=1)
    ax2.annotate("$C_0$ zero $\\approx0.31$:\nlinear model invalid",
                 xy=(0.245, 0.44), xycoords=("data", "axes fraction"),
                 fontsize=7, color="0.3", ha="right")
    ax2.annotate(r"$^6$Li sign opposite to $d$ ($Q<0$):"
                 "\na testable prediction", xy=(0.02, 0.06),
                 xycoords="axes fraction", fontsize=7, color=C_TRUTH)
    ax2.axhline(0, color="0.85", lw=0.6, zorder=0)
    ax2.set_xlim(0, 0.33)
    ax2.set_xlabel(r"$|t|$  [GeV$^2$]")
    ax2.set_ylabel(r"$a_2(|t|)$")
    ax2.set_title(r"deformation-driven $\cos2\phi'$: $d$ anchor "
                  r"$\to$ scaled $^6$Li band", fontsize=9)
    ax2.legend(fontsize=6.5, loc="upper left")
    ax2.tick_params(labelsize=8)

    # --- (c) tagged coherent yield vs x ----------------------------------
    xc = proj.x[:, 0]
    y_coh = n_coh.sum(axis=1)
    y_tag = n_tag.sum(axis=1)
    ok = y_coh > 0
    ax3.plot(xc[ok], y_coh[ok], "--", color=C_ALT, lw=1.4,
             label="coherent, produced (1 yr)")
    ax3.plot(xc[ok], np.maximum(y_tag[ok], 1e-1), "-", color=C_TRUTH,
             lw=1.8, label="tagged, 1 yr")
    ax3.plot(xc[ok], np.maximum(lumi_ratio * y_tag[ok], 1e-1), "-.",
             color="black", lw=1.2, label="tagged, 10 yr")
    for f0 in (0.02, 0.08):
        sc_b = coh.CoherentScenario(f0=f0)
        _, nb, tb = coh.project_coherent(
            config, scenario, sc_b, optics_list=(HIGH_ACCEPTANCE,))
        ax3.plot(xc[ok], np.maximum(
            tb[HIGH_ACCEPTANCE.name].sum(axis=1)[ok], 1e-1),
            "-", color=C_TRUTH, lw=0.8, alpha=0.45)
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_ylim(1e2, None)
    ax3.set_xlabel(r"$x$")
    ax3.set_ylabel(r"events per $x$ bin ($Q^2$ summed)")
    ax3.set_title(r"coherent yield: 1-year (%g) / 10-year (%g fb$^{-1}$/u)"
                  % (args.lumi_1yr, args.lumi_10yr), fontsize=9)
    ax3.annotate(r"$N_{\rm tag}^{\rm tot}$ = %.2g (1 yr), %.2g (10 yr)"
                 % (y_tag.sum(), lumi_ratio * y_tag.sum()),
                 xy=(0.05, 0.08), xycoords="axes fraction", fontsize=8)
    ax3.legend(fontsize=7, loc="upper right")
    ax3.tick_params(labelsize=8)

    # --- (d) phi modulation of the tagged coherent sample ----------------
    sel, xlo, xhi, q2lo, q2hi = best_superbin(proj, n_tag)
    n_bin = float(n_tag[sel].sum())
    a2_def = sc.a2_tagged(cut_ha, args.pzz)
    a2_def_band = [coh.CoherentScenario(eps_b0=e).a2_tagged(cut_ha,
                                                            args.pzz)
                   for e in EPS_BAND]
    counts, edges = phi_histogram_pseudo(n_bin, a2_def, nbins=24, rng=rng)
    amp_hat = cos2phi_fit_binned(counts, edges, args.pzz)
    err = cos2phi_fit_err(n_bin, args.pzz)
    err10 = cos2phi_fit_err(lumi_ratio * n_bin, args.pzz)
    centers = 0.5 * (edges[:-1] + edges[1:])
    nbar = counts.mean()
    ax4.errorbar(centers, counts / nbar,
                 yerr=np.sqrt(np.maximum(counts, 1.0)) / nbar, fmt="o",
                 color="black", ms=3.5, capsize=2, lw=1, zorder=3)
    phi = np.linspace(0, 2 * np.pi, 200)
    ax4.plot(phi, 1 + a2_def * np.cos(2 * phi), "-", color=C_TRUTH,
             lw=1.6)
    ax4.plot(phi, 1 + amp_hat * args.pzz * np.cos(2 * phi), "--",
             color=C_FIT, lw=1.4)
    ax4.plot(phi, 1 + args.amp * args.pzz * np.cos(2 * phi), "-",
             color="0.55", lw=0.9)
    ax4.annotate(r"deformation-anchored $\langle a_2\rangle=%.3f$"
                 % a2_def, xy=(0.03, 0.94), xycoords="axes fraction",
                 color=C_TRUTH, fontsize=7.5)
    ax4.annotate("binned fit", xy=(0.03, 0.87), xycoords="axes fraction",
                 color=C_FIT, fontsize=7.5)
    ax4.annotate(r"gluon-transversity scenario (%.3f)"
                 % (args.amp * args.pzz), xy=(0.03, 0.80),
                 xycoords="axes fraction", color="0.4", fontsize=7.5)
    ax4.set_xlim(0, 2 * np.pi)
    ax4.set_xticks([0, np.pi, 2 * np.pi])
    ax4.set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
    ax4.set_xlabel(r"$\phi' = \phi - \phi_S$", fontsize=9)
    ax4.set_ylabel(r"$N(\phi')/\langle N\rangle$", fontsize=9)
    ax4.set_title(
        (r"tagged sample (1 yr), $x\in[%.2g, %.2g]$, "
         r"$Q^2\in[%.2g, %.2g]$, $N=%.2g$" "\n"
         r"$\hat A=(%.1f\pm%.1f)\times10^{-3}$ [1 yr]; $\pm%.1f$ [10 yr];"
         r"  $5\sigma$ floors $%.1f/%.1f\times10^{-3}$")
        % (xlo, xhi, q2lo, q2hi, n_bin, 1e3 * amp_hat, 1e3 * err,
           1e3 * err10, 5e3 * err, 5e3 * err10), fontsize=8)
    ax4.tick_params(labelsize=8)
    ax4.axhline(1.0, color="0.85", lw=0.6, zorder=0)

    fig.suptitle(
        r"Coherent $e\,^6$Li$\,\to\,e'X\,^6$Li(g.s.), transversely "
        r"tensor-polarized, %s, $P_{zz}=%.2f$"
        "\n(scenario rates; modulation anchored on the polarized-$d$ "
        "CGC calculation, PLB 858:139053; statistical errors only, backgrounds, "
        "purity and tensor radiative corrections not included)"
        % (config.label(), args.pzz), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "money_cos2phi_coherent_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print("coherent produced: %.3g   RP-tagged (HA): %.3g   (HD): %.3g"
          % (n_coh.sum(), n_tag.sum(),
             tagged[HIGH_DIVERGENCE.name].sum()))
    print("best super-bin: x [%.3g, %.3g], Q2 [%.3g, %.3g], N_1yr=%.3g" % (
        xlo, xhi, q2lo, q2hi, n_bin))
    print("dA: %.4f (1 yr) / %.4f (10 yr)" % (err, err10))
    print("deformation <a2>_tag = %.4f (band %.4f..%.4f); "
          "gluonic scenario %.4f" % (
              a2_def, a2_def_band[0], a2_def_band[1],
              args.amp * args.pzz))
    print("A_hat = %.4f +- %.4f (truth %.4f); 5-sigma floor %.4f"
          % (amp_hat, err, a2_def / args.pzz, 5 * err))


if __name__ == "__main__":
    main()
