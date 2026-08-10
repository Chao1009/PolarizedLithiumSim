#!/usr/bin/env python3
"""Money plot 6: cos 2phi modulation in COHERENT e+6Li with the intact
recoil tagged at IP6.

Process: e + 6Li -> e' + X + 6Li(g.s.) -- the nucleus survives in its 1+
ground state and is detected in the Roman-Pot near-beam pT tail (it keeps
the exact beam rigidity, R = 1.000, so the pT cut is the only handle;
plans/06 discusses the backgrounds and the Z-discrimination this tag
needs).  The transverse tensor modulation of the coherent yield is a
SCENARIO amplitude (gluonic-quadrupole / double-helicity-flip analog);
rates come from the DIS maps scaled by the scenario coherent fraction.

Panels: (a) recoil |t| spectrum vs the RP pT cuts (both optics) with the
slope band; (b) tagged coherent yield vs x at the reference luminosity;
(c) phi'-modulation pseudo-data with projected statistical error bars in
the best low-x super-bin.

Usage:  python3 scripts/money_cos2phi_coherent.py --lumi 100 --amp 0.02
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
from polligen.estimators import cos2phi_fit_binned  # noqa: E402
from polligen.sample import phi_histogram_pseudo  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.asymmetries import err_cos2phi_amplitude  # noqa: E402
from polli_fastsim.farforward import (HIGH_ACCEPTANCE,  # noqa: E402
                                      HIGH_DIVERGENCE)

C_TRUTH, C_FIT, C_ALT = "#0072B2", "#D55E00", "#009E73"


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
    ap.add_argument("--lumi", type=float, default=100.0,
                    help="integrated luminosity [fb^-1/nucleon]")
    ap.add_argument("--pzz", type=float, default=0.60)
    ap.add_argument("--amp", type=float, default=0.02,
                    help="scenario cos2phi amplitude of the coherent "
                         "yield at Pzz=1")
    ap.add_argument("--seed", type=int, default=20260810)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs("6Li")[args.config]
    rng = np.random.default_rng(args.seed)
    sc = coh.CoherentScenario(amp=args.amp)
    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi,
                            pol_ion_tensor=args.pzz)
    proj, n_coh, tagged = coh.project_coherent(
        config, scenario, sc, optics_list=(HIGH_ACCEPTANCE, HIGH_DIVERGENCE))
    n_tag = tagged[HIGH_ACCEPTANCE.name]

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12.8, 4.4))

    # --- (a) recoil t spectrum vs the RP pT cuts --------------------------
    t = np.linspace(0.0, 0.25, 300)
    ax1.plot(t, sc.slope_b * np.exp(-sc.slope_b * t), "-", color=C_TRUTH,
             lw=1.8, label=r"$e^{-B|t|}$, $B=%g$ GeV$^{-2}$" % sc.slope_b)
    for b in (40.0, 60.0):
        ax1.plot(t, b * np.exp(-b * t), "-", color=C_TRUTH, lw=0.8,
                 alpha=0.45)
    ax1.fill_between(t, 1e-3, sc.slope_b * np.exp(-sc.slope_b * t),
                     where=t >= HIGH_ACCEPTANCE.pt_cut_near_beam ** 2,
                     color=C_TRUTH, alpha=0.12)
    for optics, color in ((HIGH_ACCEPTANCE, C_FIT),
                          (HIGH_DIVERGENCE, "0.35")):
        cut2 = optics.pt_cut_near_beam ** 2
        acc = sc.tag_acceptance(optics.pt_cut_near_beam)
        band = [coh.CoherentScenario(slope_b=b).tag_acceptance(
            optics.pt_cut_near_beam) for b in (60.0, 40.0)]
        ax1.axvline(cut2, color=color, ls="--", lw=1.2)
        ax1.annotate(
            "%s\n$p_T>%.2f$ GeV\nacc = %.3g\n[%.2g, %.2g]"
            % (optics.name, optics.pt_cut_near_beam, acc, *band),
            xy=(cut2, 2e-3), xytext=(4, 6), textcoords="offset points",
            fontsize=7, color=color, va="bottom")
    ax1.set_yscale("log")
    ax1.set_ylim(1e-3, 80)
    ax1.set_xlabel(r"$|t| \simeq p_T^2$  [GeV$^2$]")
    ax1.set_ylabel(r"$(1/\sigma)\, d\sigma_{\rm coh}/d|t|$  [GeV$^{-2}$]")
    ax1.set_title("intact-recoil $t$ spectrum vs RP $p_T$ cut", fontsize=9)
    ax1.legend(fontsize=7, loc="upper right")
    ax1.tick_params(labelsize=8)

    # --- (b) tagged coherent yield vs x ----------------------------------
    xc = proj.x[:, 0]
    y_coh = n_coh.sum(axis=1)
    y_tag = n_tag.sum(axis=1)
    ok = y_coh > 0
    ax2.plot(xc[ok], y_coh[ok], "--", color=C_ALT, lw=1.4,
             label="coherent, produced")
    ax2.plot(xc[ok], np.maximum(y_tag[ok], 1e-1), "-", color=C_TRUTH,
             lw=1.8, label="RP-tagged (high-acc.)")
    for f0 in (0.02, 0.08):
        sc_b = coh.CoherentScenario(f0=f0)
        _, nb, tb = coh.project_coherent(
            config, scenario, sc_b, optics_list=(HIGH_ACCEPTANCE,))
        ax2.plot(xc[ok], np.maximum(
            tb[HIGH_ACCEPTANCE.name].sum(axis=1)[ok], 1e-1),
            "-", color=C_TRUTH, lw=0.8, alpha=0.45)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_ylim(1e2, None)
    ax2.set_xlabel(r"$x$")
    ax2.set_ylabel(r"events per $x$ bin ($Q^2$ summed)")
    ax2.set_title(r"coherent yield at $L=%g$ fb$^{-1}$/u"
                  % args.lumi, fontsize=9)
    ax2.annotate(r"$N_{\rm tag}^{\rm tot} = %.2g$" % y_tag.sum(),
                 xy=(0.05, 0.08), xycoords="axes fraction", fontsize=8)
    ax2.legend(fontsize=7, loc="upper right")
    ax2.tick_params(labelsize=8)

    # --- (c) phi modulation of the tagged coherent sample ----------------
    sel, xlo, xhi, q2lo, q2hi = best_superbin(proj, n_tag)
    n_bin = float(n_tag[sel].sum())
    a2_eff = args.amp * args.pzz
    counts, edges = phi_histogram_pseudo(n_bin, a2_eff, nbins=24, rng=rng)
    amp_hat = cos2phi_fit_binned(counts, edges, args.pzz)
    err = err_cos2phi_amplitude(n_bin, args.pzz)
    centers = 0.5 * (edges[:-1] + edges[1:])
    nbar = counts.mean()
    ax3.errorbar(centers, counts / nbar,
                 yerr=np.sqrt(np.maximum(counts, 1.0)) / nbar, fmt="o",
                 color="black", ms=3.5, capsize=2, lw=1, zorder=3)
    phi = np.linspace(0, 2 * np.pi, 200)
    ax3.plot(phi, 1 + a2_eff * np.cos(2 * phi), "-", color=C_TRUTH, lw=1.6)
    ax3.plot(phi, 1 + amp_hat * args.pzz * np.cos(2 * phi), "--",
             color=C_FIT, lw=1.4)
    ax3.annotate(r"scenario $A^{\rm coh}=%g$" % args.amp, xy=(0.03, 0.94),
                 xycoords="axes fraction", color=C_TRUTH, fontsize=7.5)
    ax3.annotate("binned fit", xy=(0.03, 0.87), xycoords="axes fraction",
                 color=C_FIT, fontsize=7.5)
    ax3.set_xlim(0, 2 * np.pi)
    ax3.set_xticks([0, np.pi, 2 * np.pi])
    ax3.set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
    ax3.set_xlabel(r"$\phi' = \phi - \phi_S$", fontsize=9)
    ax3.set_ylabel(r"$N(\phi')/\langle N\rangle$", fontsize=9)
    ax3.set_title(
        (r"tagged sample, $x\in[%.2g, %.2g]$, $Q^2\in[%.2g, %.2g]$" "\n"
         r"$N=%.2g$;  $\hat A=(%.1f\pm%.1f)\times10^{-3}$;  "
         r"$5\sigma$ floor $%.1f\times10^{-3}$")
        % (xlo, xhi, q2lo, q2hi, n_bin, 1e3 * amp_hat, 1e3 * err,
           5e3 * err), fontsize=8)
    ax3.tick_params(labelsize=8)
    ax3.axhline(1.0, color="0.85", lw=0.6, zorder=0)

    fig.suptitle(
        r"Coherent $e\,^6$Li$\,\to\,e'X\,^6$Li(g.s.), transversely "
        r"tensor-polarized, %s, $P_{zz}=%.2f$ "
        "(SCENARIO coherent model; backgrounds in plans/06)"
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
    print("best super-bin: x [%.3g, %.3g], Q2 [%.3g, %.3g], N=%.3g" % (
        xlo, xhi, q2lo, q2hi, n_bin))
    print("A_hat = %.4f +- %.4f (truth %.4f); 5-sigma floor %.4f"
          % (amp_hat, err, args.amp, 5 * err))


if __name__ == "__main__":
    main()
