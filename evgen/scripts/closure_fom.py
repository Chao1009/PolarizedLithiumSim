"""Step-5.A deliverable: fastsim <-> generator closure of the FOM maps.

For each Phase-1 observable (A_par, A_zz, A_cos2phi) run many pseudo-
experiments at a reduced luminosity, form the per-x-bin estimator (q2
pooled inside the bin), and compare the trial-to-trial spread with the
analytic error formulas of `polli_fastsim.asymmetries` evaluated at the
same event counts -- the closure that certifies both the analytic FOM
maps (plans/02 step 1.3) and the generator chain (plans/05 step 5.A).

Usage:
    python3 scripts/closure_fom.py --ion 6Li --events 200000 --trials 80
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from polligen import bookkeeping as bk  # noqa: E402
from polligen import estimators as est  # noqa: E402
from polligen.sample import InclusiveSampler, run_pseudo_experiment  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.asymmetries import (  # noqa: E402
    err_a_parallel, err_azz, err_cos2phi_amplitude)
from polli_fastsim.polarized import toy_b1, toy_delta_gluon  # noqa: E402


def x_binned_counts(events, x_edges):
    return np.histogram(events["x"], bins=x_edges)[0]


def run_observable(name, sampler, plan, trials, lumi_pb, rng):
    """Per-x-bin estimator values and counts across pseudo-experiments."""
    x_edges = sampler.proj.x_edges
    nxb = x_edges.size - 1
    vals = np.full((trials, nxb), np.nan)
    ntot = np.zeros((trials, nxb))
    for it in range(trials):
        ev, _ = run_pseudo_experiment(sampler, plan, lumi_pb, rng=rng)
        if name == "apar":
            np_ = x_binned_counts(ev["apar+"], x_edges)
            nm = x_binned_counts(ev["apar-"], x_edges)
            n = np_ + nm
            with np.errstate(divide="ignore", invalid="ignore"):
                vals[it] = np.where(
                    n > 0, est.apar_flip(np_, nm, PE, PZ), np.nan)
        elif name == "azz":
            n3 = {k: x_binned_counts(ev["azz" + k], x_edges) for k in "+-0"}
            n = n3["+"] + n3["-"] + n3["0"]
            with np.errstate(divide="ignore", invalid="ignore"):
                vals[it] = np.where(
                    n > 0, est.azz_thirds(n3["+"], n3["-"], n3["0"], PZZ),
                    np.nan)
        else:  # cos2phi: moment per x bin
            cat = plan.categories[0]
            phi = ev[cat.name]["phi"] - cat.phi_s
            xb = np.digitize(ev[cat.name]["x"], x_edges) - 1
            n = x_binned_counts(ev[cat.name], x_edges)
            for b in range(nxb):
                sel = xb == b
                if sel.sum() > 1:
                    vals[it, b] = est.cos2phi_moment(phi[sel], PZZ)
        ntot[it] = n
    return vals, ntot


def analytic_err(name, n):
    if name == "apar":
        return err_a_parallel(n, PE, PZ)
    if name == "azz":
        return err_azz(n, PZZ)
    return err_cos2phi_amplitude(n, PZZ)


PE, PZ, PZZ = 0.7, 0.6, 0.6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="6Li")
    ap.add_argument("--config-index", type=int, default=1)
    ap.add_argument("--events", type=float, default=2e5,
                    help="target events per pseudo-experiment")
    ap.add_argument("--trials", type=int, default=80)
    ap.add_argument("--min-events", type=float, default=400,
                    help="per-x-bin count floor for the comparison")
    ap.add_argument("--seed", type=int, default=20260713)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs(args.ion)[args.config_index]
    print("closure on %s, %d trials x %g events" %
          (config.label(), args.trials, args.events))

    plans = {
        "apar": (InclusiveKernel(beams.IONS[args.ion]),
                 bk.helicity_flip_plan(beams.IONS[args.ion].spin, PZ, PE)),
        "azz": (InclusiveKernel(beams.LI6, b1_func=toy_b1),
                bk.tensor_thirds_plan(PZ, PZZ)),
        "cos2phi": (InclusiveKernel(
            beams.LI6,
            delta_func=lambda x, q2, f1: toy_delta_gluon(x, q2, f1,
                                                         scale=1e-2)),
            bk.transverse_tensor_plan(PZZ)),
    }
    # note: the tensor observables (azz, cos2phi) always run the spin-1
    # 6Li kernel; --ion only switches the vector (apar) beam and config.
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(13, 6), sharex=True,
                             gridspec_kw={"height_ratios": (2, 1)})
    labels = {"apar": r"$\delta A_\parallel$", "azz": r"$\delta A_{zz}$",
              "cos2phi": r"$\delta A_{\cos 2\phi}$"}
    summary = []
    for col, name in enumerate(("apar", "azz", "cos2phi")):
        kern, plan = plans[name]
        sampler = InclusiveSampler(kern, config, fom.Scenario())
        sigma = sum(sampler.sigma_tot_pb(c) * c.lumi_fraction
                    for c in plan.categories)
        lumi_pb = args.events / sigma
        rng = np.random.default_rng(args.seed + col)
        vals, ntot = run_observable(name, sampler, plan, args.trials,
                                    lumi_pb, rng)
        n_mean = np.nanmean(ntot, axis=0)
        use = n_mean >= args.min_events
        spread = np.nanstd(vals, axis=0)
        analytic = analytic_err(name, np.maximum(n_mean, 1))
        ratio = np.where(use, spread / analytic, np.nan)
        xc = np.sqrt(sampler.proj.x_edges[:-1] * sampler.proj.x_edges[1:])
        band = 1.96 / np.sqrt(2 * args.trials)

        ax = axes[0, col]
        ax.loglog(xc[use], analytic[use], "k-", label="analytic (fastsim)")
        ax.loglog(xc[use], spread[use], "o", ms=4, mfc="none",
                  label="pseudo-exp spread")
        ax.set_title("%s  (%s)" % (labels[name], name))
        ax.legend(fontsize=8)
        axr = axes[1, col]
        axr.axhspan(1 - band, 1 + band, color="0.9")
        axr.axhline(1.0, color="k", lw=0.8)
        axr.semilogx(xc[use], ratio[use], "o", ms=4)
        axr.set_ylim(0.6, 1.4)
        axr.set_xlabel("x")
        if col == 0:
            axes[0, col].set_ylabel("per-x-bin uncertainty")
            axr.set_ylabel("spread / analytic")
        worst_i = np.nanargmax(np.abs(ratio - 1.0))
        worst = ratio[worst_i] - 1.0
        nbin = int(use.sum())
        summary.append(
            "%-8s: %2d x-bins, max |spread/analytic - 1| = %+.2f "
            "(band %.2f) at x=%.3g, <N>=%.0f"
            % (name, nbin, worst, band, xc[worst_i], n_mean[worst_i]))

    fig.suptitle("Generator <-> fastsim FOM closure, %s (toy SFs)"
                 % config.label())
    fig.tight_layout()
    out = outdir / ("closure_fom_%s.png" % args.ion)
    fig.savefig(out, dpi=140)
    print("wrote", out)
    for line in summary:
        print(line)


if __name__ == "__main__":
    main()
