#!/usr/bin/env python3
"""Money plot 5: the cos 2phi gluonometry measurement as projected data.

Left 2x2: phi'-modulation pseudo-data (points with statistical error bars
at the 1-YEAR EIC luminosity) in the four sweet-spot (x, Q2) super-bins --
the bins with the highest projected cos 2phi significance, spread in
(x, Q2).  Right: the extracted amplitude vs x along the sweet-spot Q2
slice with 1-year AND 10-year error bars and the Delta-model curves.

The Delta(x, Q2) input comes from the unified model registry
`polli_fastsim.delta_models` (single home for all Delta models and
constraints): default `moment_A` -- the sum-rule-constrained ansatz of
the money_delta suite, Delta = A alpha_s F1 x^a(1-x)^b with A solved
from int x Delta dx = -0.012 alpha_s (Sather-Schmidt bag moment) --
with the 6Li per-nucleon dilution 1/3 (two polarized nucleons of six;
the whole-nucleus P_zz stays a separate factor).  It is a counting
fraction, not the vector polarization of `beams.LI6` -- see the
convention note in `delta_models`.  `moment_B` (no F1 factor,
the conservative reading) and the legacy `toy` shape are one flag away.

Detection assumption: scattered electron only (Scenario e' cuts).
Statistics via per-phi-bin Poisson draws from exact expected yields
(sample.phi_histogram_pseudo) -- identical to event-level sampling for
the binned estimator at any luminosity.  The whole luminosity is
assigned to the transverse-tensor fill (FOM-map convention).

Usage:  python3 scripts/money_cos2phi.py --delta-model moment_A
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
from polligen.estimators import (cos2phi_fit_binned,  # noqa: E402
                                 cos2phi_fit_err)
from polligen.sample import InclusiveSampler, phi_histogram_pseudo  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, delta_models as dm, fom  # noqa: E402
PRETTY = {"moment_A": "moment-constrained $\\Delta$", "moment_B": "no-$F_1$ variant", "toy": "flat toy"}  # figure labels
from polli_fastsim.asymmetries import a_cos2phi  # noqa: E402
from polli_fastsim.kinematics import kinematic_mask, y_from_xq2  # noqa: E402
from polli_fastsim.polarized import toy_b1  # noqa: E402
from polli_fastsim.structure import NuclearF2  # noqa: E402

# Okabe-Ito (colorblind-safe): blue = injected, vermillion = fit/1-yr,
# green = alternative interpretation
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
    Q2 lever arm, so one spot is forced into each higher band.

    Audit note (2026-08-10): the ranking runs on the 40x30 analysis map
    with cell-center acceptance; counts of bins touching the y_min edge
    carry ~20% grid bias there.  All PLOTTED numbers instead come from
    the finer 60x45 sampler grid (within ~3% of a 240x180 reference),
    so the bias can only influence which near-edge bin gets picked."""
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
    """One full-luminosity binned pseudo-measurement on a super-bin."""
    sigma_pb, _a1, a2_eff = sampler.effective_modulation(cat, mask=mask)
    n_exp = lumi_pb * sigma_pb
    counts, edges = phi_histogram_pseudo(n_exp, a2_eff, nbins=nbins, rng=rng)
    amp_hat = cos2phi_fit_binned(counts, edges, pzz)
    return {"n": n_exp, "truth": a2_eff / pzz, "amp": amp_hat,
            "err": cos2phi_fit_err(n_exp, pzz, nbins),
            "counts": counts, "edges": edges, "a2_eff": a2_eff,
            "sigma_pb": sigma_pb}


def build_delta_model(args, config, scenario):
    """Delta model from the unified registry; moment_A needs the
    rate-weighted <Q2> of the accepted phase space and a per-nucleon
    F1 handle."""
    if args.delta_model == "toy":
        return dm.make("toy", scale=args.scale), None
    nf2 = NuclearF2(beams.LI6)
    if args.delta_model == "moment_B":
        return dm.make("moment_B", variant=args.variant,
                       dilution=args.dilution), None
    proj0 = fom.project_rates(config, scenario)
    acc = proj0.accepted
    q2_ref = float((proj0.n_events[acc] * proj0.q2[acc]).sum()
                   / proj0.n_events[acc].sum())
    model = dm.make(
        "moment_A", f1_func=lambda x, q2: nf2.f1a(x, q2) / beams.LI6.A,
        q2_ref=q2_ref, variant=args.variant, dilution=args.dilution)
    return model, q2_ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2),
                    help="beam-energy point (default mid, e10 x 6Li 99.5 GeV/u)")
    ap.add_argument("--delta-model", default="moment_A",
                    choices=dm.available())
    ap.add_argument("--variant", default="mid_x",
                    choices=sorted(dm.VARIANTS))
    ap.add_argument("--dilution", type=float, default=1.0 / 3.0,
                    help="6Li per-nucleon dilution folded into Delta "
                         "(2 of 6 nucleons); P_zz stays whole-nucleus")
    ap.add_argument("--scale", type=float, default=1e-2,
                    help="peak Delta/F1 (toy model only)")
    ap.add_argument("--lumi-1yr", type=float, default=10.0,
                    help="1-year EIC program [fb^-1/nucleon]")
    ap.add_argument("--lumi-10yr", type=float, default=100.0,
                    help="10-year EIC program [fb^-1/nucleon]")
    ap.add_argument("--pzz", type=float, default=0.60)
    ap.add_argument("--nspots", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20260810)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs("6Li")[args.config]
    lumi1_pb = args.lumi_1yr * 1e3
    lumi10_pb = args.lumi_10yr * 1e3
    rng = np.random.default_rng(args.seed)

    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    model, q2_ref = build_delta_model(args, config, scenario)
    delta_func = model  # DeltaModel is (x, q2, f1)-callable

    # --- sweet spots from the analytic significance map -------------------
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
        m = measure(sampler, cat, mask, lumi1_pb, plan.pzz_true, rng)
        err10 = cos2phi_fit_err(lumi10_pb * m["sigma_pb"],
                                plan.pzz_true)
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
             r"$\hat A=(%.2f\pm%.2f)\times10^{-3}$ [1 yr]; "
             r"$\pm%.2f$ [10 yr]")
            % (xs, qs, 1e3 * m["amp"], 1e3 * m["err"], 1e3 * err10),
            xy=(0.5, 1.02), xycoords="axes fraction", ha="center",
            fontsize=7.5)
        if k == 0:
            ax.annotate("injected (%s)" % PRETTY.get(args.delta_model, args.delta_model),
                        xy=(0.03, 0.93), xycoords="axes fraction",
                        color=C_TRUTH, fontsize=7)
            ax.annotate("binned fit", xy=(0.03, 0.84),
                        xycoords="axes fraction", color=C_FIT, fontsize=7)
        summary.append(
            "spot %d: x=%.3g Q2=%.3g  N_1yr=%.2e  A_truth=%+.2e  "
            "A_hat=%+.2e +- %.1e (1yr) +- %.1e (10yr)"
            % (k + 1, xs, qs, m["n"], m["truth"], m["amp"], m["err"],
               err10))

    # --- right: amplitude vs x along the sweet-spot Q2 slice --------------
    ax = fig.add_subplot(gs[:, 2])
    q2_spot = spots[0][1]
    q2lo, q2hi = q2_spot / 1.6, q2_spot * 1.6
    xe = proj.x_edges
    pts = []
    for i0 in range(0, xe.size - 2, 2):  # merge pairs of x bins
        xc = np.sqrt(xe[i0] * xe[i0 + 2])
        if not kinematic_mask(xc, q2_spot, sampler.s):
            continue
        mask = superbin_mask(sampler, xe[i0], xe[i0 + 2], q2lo, q2hi)
        if not mask.any():
            continue
        m = measure(sampler, cat, mask, lumi1_pb, plan.pzz_true, rng)
        if m["n"] < 1e3 or m["err"] > 8e-3:
            continue
        # independent 10-yr pseudo-measurement (2026-08-11 audit: reusing
        # the 1-yr draw with 10-yr bars scattered points ~sqrt(10) sigma)
        m10 = measure(sampler, cat, mask, lumi10_pb, plan.pzz_true, rng)
        m["amp10"] = m10["amp"]
        m["err10"] = m10["err"]
        pts.append((xc, m))
    xoff = 1.045  # slight offset so the two series stay legible
    ax.errorbar([p[0] for p in pts], [1e3 * p[1]["amp"] for p in pts],
                yerr=[1e3 * p[1]["err"] for p in pts], fmt="s",
                mfc="none", color=C_FIT, ms=4, capsize=2, lw=1, zorder=3,
                label=r"1-year EIC (%g fb$^{-1}$/u)" % args.lumi_1yr)
    ax.errorbar([p[0] * xoff for p in pts],
                [1e3 * p[1]["amp10"] for p in pts],
                yerr=[1e3 * p[1]["err10"] for p in pts], fmt="o",
                color="black", ms=4, capsize=2, lw=1, zorder=4,
                label=r"10-year EIC (%g fb$^{-1}$/u)" % args.lumi_10yr)

    xg = np.logspace(np.log10(2e-4), np.log10(0.8), 250)
    q2g = np.full_like(xg, q2_spot)
    ok = kinematic_mask(xg, q2g, sampler.s)
    f2 = kern.nf2.f2a(xg, q2g) / kern.ion.A
    f1 = kern.nf2.f1a(xg, q2g) / kern.ion.A
    y = y_from_xq2(xg, q2g, sampler.s)
    curves = [(delta_func, C_TRUTH, "-",
               "%s (injected)" % PRETTY.get(args.delta_model, args.delta_model))]
    if args.delta_model != "moment_B":
        alt = dm.make("moment_B", variant=args.variant,
                      dilution=args.dilution)
        curves.append((alt, C_ALT, "-", "no-$F_1$ variant (conservative)"))
    if args.delta_model != "toy":
        curves.append((dm.make("toy", scale=1e-3), "0.45", "--",
                       r"flat toy, $\Delta/F_1=10^{-3}$"))
    for dfunc, color, ls, lab in curves:
        amp = a_cos2phi(dfunc(xg, q2g, f1), f1, f2, xg, y)
        ax.plot(xg[ok], 1e3 * amp[ok], ls, color=color, lw=1.5, label=lab)
    ax.set_xscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$A^{\cos 2\phi}$  $[\times 10^{-3}]$")
    ax.axhline(0, color="0.85", lw=0.6, zorder=0)
    ax.tick_params(labelsize=8)
    ax.set_title(r"amplitude vs $x$,  $Q^2\approx%.3g$ GeV$^2$" % q2_spot,
                 fontsize=9)
    ax.legend(fontsize=7, loc="upper left")

    fig.suptitle(
        r"Nuclear gluonometry, transversely tensor-polarized $^6$Li, "
        r"%s, $P_{zz}=%.2f$""\n"
        r"$\Delta$: %s;  $\phi'$ pseudo-data at 1 yr (%g fb$^{-1}$/u); "
        "statistical errors only, no backgrounds or tensor radiative corrections"
        % (config.label(), plan.pzz_true,
           PRETTY.get(args.delta_model, args.delta_model)
           + (" (bag moment, dilution 1/3)" if args.delta_model == "moment_A" else ""),
           args.lumi_1yr),
        fontsize=10)
    fig.subplots_adjust(top=0.86, bottom=0.09, left=0.06, right=0.985)
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "money_cos2phi_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print("delta model:", model.info(),
          "" if q2_ref is None else "(<Q2> = %.3g GeV^2)" % q2_ref)
    for line in summary:
        print(line)


if __name__ == "__main__":
    main()
