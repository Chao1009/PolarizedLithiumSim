#!/usr/bin/env python3
"""Money plot 5: the cos 2phi gluonometry measurement as projected data.

Left 2x2: phi'-modulation pseudo-data (points with statistical error bars
at the reference luminosity) in the four sweet-spot (x, Q2) super-bins --
the bins with the highest projected cos 2phi significance, spread in
(x, Q2).  Right: the extracted amplitude vs x along the sweet-spot Q2
slice, with Delta/F1 scenario curves.

Detection assumption: scattered electron only (the inclusive gluonometry
needs no other tag; the Scenario e' acceptance cuts are applied).  The
companion coherent-channel plot (money_cos2phi_coherent.py) adds the
intact-6Li far-forward tag.

Statistics are generated as per-phi-bin Poisson draws from the exact
expected yields (sample.phi_histogram_pseudo) -- identical to event-level
sampling for the binned estimator, at any luminosity.  The whole
luminosity is assigned to the transverse-tensor fill (same convention as
the analytic FOM maps).

Usage:  python3 scripts/money_cos2phi.py --lumi 100 --scale 1e-2
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

from polligen import bookkeeping as bk  # noqa: E402
from polligen.estimators import cos2phi_fit_binned  # noqa: E402
from polligen.sample import InclusiveSampler, phi_histogram_pseudo  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.asymmetries import (a_cos2phi,  # noqa: E402
                                       err_cos2phi_amplitude)
from polli_fastsim.kinematics import kinematic_mask, y_from_xq2  # noqa: E402
from polli_fastsim.polarized import toy_b1, toy_delta_gluon  # noqa: E402

# Okabe-Ito (colorblind-safe): blue = truth, vermillion = fit, green = alt
C_TRUTH, C_FIT, C_ALT = "#0072B2", "#D55E00", "#009E73"


def pick_sweet_spots(proj, sig, n=4, min_dx=0.35):
    """Greedy top-significance bins, separated by min_dx in log10(x)."""
    order = np.argsort(sig, axis=None)[::-1]
    picked = []
    for flat in order:
        i, j = np.unravel_index(flat, sig.shape)
        if sig[i, j] <= 0:
            break
        x, q2 = proj.x[i, j], proj.q2[i, j]
        if not any(abs(np.log10(x / xs)) < min_dx
                   for xs, *_ in picked):
            picked.append((x, q2, i, j))
        if len(picked) == n:
            break
    return picked


def pick_sweet_spots_banded(proj, sig,
                            bands=((1.0, 3.0, 2), (3.0, 12.0, 1),
                                   (12.0, 200.0, 1))):
    """Best bins per Q2 band: the raw significance map clusters at the
    lowest Q2 (rate-dominated), but the money plot must also show the
    Q2 lever arm, so one spot is forced into each higher band."""
    picked = []
    for qlo, qhi, npick in bands:
        in_band = np.where((proj.q2 >= qlo) & (proj.q2 < qhi), sig, 0.0)
        picked += pick_sweet_spots(proj, in_band, n=npick)
    return picked


def superbin_edges(proj, i, j, pad=1):
    """(x_lo, x_hi, q2_lo, q2_hi) of the FOM bin padded by `pad` bins."""
    xe, qe = proj.x_edges, proj.q2_edges
    return (xe[max(i - pad, 0)], xe[min(i + pad + 1, xe.size - 1)],
            qe[max(j - pad, 0)], qe[min(j + pad + 1, qe.size - 1)])


def superbin_mask(sampler, xlo, xhi, q2lo, q2hi):
    return ((sampler.x_cells >= xlo) & (sampler.x_cells < xhi)
            & (sampler.q2_cells >= q2lo) & (sampler.q2_cells < q2hi))


def measure(sampler, cat, mask, lumi_pb, pzz, rng, nbins=24):
    """One full-luminosity binned pseudo-measurement on a super-bin.

    Returns dict with N, truth amplitude (physics convention, i.e. the
    P_zz-divided a_cos2phi), extracted amplitude, its error, and the
    phi histogram (counts, edges)."""
    sigma_pb, _a1, a2_eff = sampler.effective_modulation(cat, mask=mask)
    n_exp = lumi_pb * sigma_pb
    counts, edges = phi_histogram_pseudo(n_exp, a2_eff, nbins=nbins, rng=rng)
    amp_hat = cos2phi_fit_binned(counts, edges, pzz)
    return {"n": n_exp, "truth": a2_eff / pzz, "amp": amp_hat,
            "err": err_cos2phi_amplitude(n_exp, pzz),
            "counts": counts, "edges": edges, "a2_eff": a2_eff}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2),
                    help="beam-energy point (default mid, e10 x 6Li50)")
    ap.add_argument("--lumi", type=float, default=100.0,
                    help="integrated luminosity [fb^-1/nucleon]")
    ap.add_argument("--pzz", type=float, default=0.60)
    ap.add_argument("--scale", type=float, default=1e-2,
                    help="Delta/F1 scenario scale generated as truth")
    ap.add_argument("--nspots", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20260810)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs("6Li")[args.config]
    lumi_pb = args.lumi * 1e3
    rng = np.random.default_rng(args.seed)

    def delta_func(x, q2, f1):
        return toy_delta_gluon(x, q2, f1, scale=args.scale)

    # --- sweet spots from the analytic significance map -------------------
    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi,
                            pol_ion_tensor=args.pzz)
    proj = fom.project_rates(config, scenario)
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1, delta_func=delta_func)
    obs = fom.project_observables(config, scenario, proj,
                                  kern.g1_model, toy_b1, delta_func)
    spots = pick_sweet_spots_banded(proj, obs["sig_a_cos2phi"])[:args.nspots]

    # --- one sampler over the full range; super-bins are cell masks ------
    sampler = InclusiveSampler(kern, config, scenario, nx=60, nq2=45)
    plan = bk.transverse_tensor_plan(args.pzz)
    cat = plan.categories[0]

    fig = plt.figure(figsize=(12.5, 6.8))
    gs = GridSpec(2, 3, figure=fig, width_ratios=(1, 1, 1.35),
                  hspace=0.42, wspace=0.30)
    summary = []

    for k, (xs, qs, i, j) in enumerate(spots):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        xlo, xhi, q2lo, q2hi = superbin_edges(proj, i, j)
        mask = superbin_mask(sampler, xlo, xhi, q2lo, q2hi)
        m = measure(sampler, cat, mask, lumi_pb, plan.pzz_true, rng)
        centers = 0.5 * (m["edges"][:-1] + m["edges"][1:])
        nbar = m["counts"].mean()
        mod = 1e3 * (m["counts"] / nbar - 1.0)
        mod_err = 1e3 * np.sqrt(np.maximum(m["counts"], 1.0)) / nbar
        ax.errorbar(centers, mod, yerr=mod_err, fmt="o", color="black",
                    ms=3.5, capsize=2, lw=1, zorder=3)
        phi = np.linspace(0, 2 * np.pi, 200)
        ax.plot(phi, 1e3 * m["a2_eff"] * np.cos(2 * phi), "-",
                color=C_TRUTH, lw=1.6)
        fit_amp = m["amp"] * plan.pzz_true
        ax.plot(phi, 1e3 * fit_amp * np.cos(2 * phi), "--", color=C_FIT,
                lw=1.4)
        ax.plot(phi, 1e3 * 0.1 * m["a2_eff"] * np.cos(2 * phi), "-",
                color="0.55", lw=0.9)
        ax.set_xlim(0, 2 * np.pi)
        ax.set_xticks([0, np.pi, 2 * np.pi])
        ax.set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
        ax.set_xlabel(r"$\phi' = \phi - \phi_S$", fontsize=9, labelpad=1)
        if k % 2 == 0:
            ax.set_ylabel(r"$N(\phi')/\langle N\rangle - 1$"
                          r"  $[\times 10^{-3}]$", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.axhline(0.0, color="0.85", lw=0.6, zorder=0)
        ax.annotate(
            (r"$x\approx%.3g$, $Q^2\approx%.3g$ GeV$^2$" "\n"
             r"$N=%.1f\times10^{%d}$;  "
             r"$\hat A=(%.2f\pm%.2f)\times10^{-3}$")
            % (xs, qs, m["n"] / 10 ** int(np.log10(m["n"])),
               int(np.log10(m["n"])), 1e3 * m["amp"], 1e3 * m["err"]),
            xy=(0.5, 1.02), xycoords="axes fraction", ha="center",
            fontsize=7.5)
        if k == 0:
            # direct labels (fixed spots; the modulation shape is common)
            ax.annotate(r"injected $\Delta/F_1{=}%g$" % args.scale,
                        xy=(0.03, 0.93), xycoords="axes fraction",
                        color=C_TRUTH, fontsize=7)
            ax.annotate("binned fit", xy=(0.03, 0.84),
                        xycoords="axes fraction", color=C_FIT, fontsize=7)
            ax.annotate(r"$\Delta/F_1{=}%g$" % (0.1 * args.scale),
                        xy=(0.60, 0.57), xycoords="axes fraction",
                        color="0.4", fontsize=7)
        summary.append(
            "spot %d: x=%.3g Q2=%.3g  N=%.2e  A_truth=%+.2e  "
            "A_hat=%+.2e +- %.1e  (5sig min amp %.1e)"
            % (k + 1, xs, qs, m["n"], m["truth"], m["amp"], m["err"],
               5 * m["err"]))

    # --- right: amplitude vs x along the sweet-spot Q2 slice --------------
    ax = fig.add_subplot(gs[:, 2])
    q2_spot = spots[0][1]
    q2lo, q2hi = q2_spot / 1.6, q2_spot * 1.6
    xe = proj.x_edges
    pts = []
    for i0 in range(0, xe.size - 2, 2):  # merge pairs of x bins
        xc = np.sqrt(xe[i0] * xe[i0 + 2])
        if not kinematic_mask(xc, q2_spot, sampler.s):
            continue  # keep points on the fixed-Q2 curve's kinematic range
        mask = superbin_mask(sampler, xe[i0], xe[i0 + 2], q2lo, q2hi)
        if not mask.any():
            continue
        m = measure(sampler, cat, mask, lumi_pb, plan.pzz_true, rng)
        if m["n"] < 1e4 or m["err"] > 8e-3:
            continue
        pts.append((xc, m))
    ax.errorbar([p[0] for p in pts], [1e3 * p[1]["amp"] for p in pts],
                yerr=[1e3 * p[1]["err"] for p in pts], fmt="o",
                color="black", ms=4, capsize=2, lw=1, zorder=3,
                label="projected data")

    xg = np.logspace(np.log10(2e-4), np.log10(0.8), 250)
    q2g = np.full_like(xg, q2_spot)
    ok = kinematic_mask(xg, q2g, sampler.s)
    f2 = kern.nf2.f2a(xg, q2g) / kern.ion.A
    f1 = kern.nf2.f1a(xg, q2g) / kern.ion.A
    y = y_from_xq2(xg, q2g, sampler.s)
    for sc, color, ls, lab in (
            (args.scale, C_TRUTH, "-",
             r"$\Delta/F_1=%g$ (injected)" % args.scale),
            (args.scale * 0.3, C_ALT, "-",
             r"$\Delta/F_1=%g$" % (0.3 * args.scale)),
            (args.scale * 0.1, "0.45", "--",
             r"$\Delta/F_1=%g$" % (0.1 * args.scale))):
        amp = a_cos2phi(toy_delta_gluon(xg, q2g, f1, scale=sc),
                        f1, f2, xg, y)
        ax.plot(xg[ok], 1e3 * amp[ok], ls, color=color, lw=1.5, label=lab)
    ax.set_xscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$A^{\cos 2\phi}$  $[\times 10^{-3}]$")
    ax.axhline(0, color="0.85", lw=0.6, zorder=0)
    ax.tick_params(labelsize=8)
    ax.set_title(r"amplitude vs $x$,  $Q^2\approx%.3g$ GeV$^2$" % q2_spot,
                 fontsize=9)
    ax.legend(fontsize=7, loc="lower left")

    fig.suptitle(
        r"Nuclear gluonometry, transversely tensor-polarized $^6$Li, "
        "%s\n"
        r"pseudo-data at $L=%g$ fb$^{-1}$/u, $P_{zz}=%.2f$ "
        "(all luminosity in the transverse fill; TOY/scenario inputs)"
        % (config.label(), args.lumi, plan.pzz_true), fontsize=10)
    fig.subplots_adjust(top=0.86, bottom=0.09, left=0.06, right=0.985)
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "money_cos2phi_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    for line in summary:
        print(line)


if __name__ == "__main__":
    main()
