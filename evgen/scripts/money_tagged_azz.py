"""Money plot 4: tagged tensor asymmetry A_zz^tag for the alpha-tagged
embedded deuteron in 6Li -- the first tagged spin observable projected for
any A > 2 (plans/05 step 5.B; extends the four-gap list of plans/00).

Left: the analytic wave-function asymmetry A_zz^wf(k, theta_k = 90 deg)
from the alpha-d S/D interference, with the model band (beta and P_D
scans -- the same tail systematics that dominate the 6Li RP acceptance).
Right: acceptance-folded pseudo-experiment -- tensor-thirds fills, events
routed through the far-forward windows, Roman-Pot-accepted (main window +
pT tail) thirds estimator vs the spectator momentum k, for both optics.
The plan's question -- does the Cosyn-Weiss-style O(1) asymmetry at
p_s ~ 300 MeV/c survive the RP pT-tail acceptance that dominates 6Li? --
is answered by the right panel.

Usage:  python3 scripts/money_tagged_azz.py --events 400000
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
from polligen import tagged  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.asymmetries import err_azz  # noqa: E402
from polli_fastsim.farforward import (  # noqa: E402
    HIGH_ACCEPTANCE, HIGH_DIVERGENCE, route_charged)
from polli_fastsim.polarized import toy_b1  # noqa: E402

PZ, PZZ = 0.7, 0.6


def azz_wf_curve(model, ic):
    n = {m: model.n_of_kc(m)[:, ic] for m in (1.0, 0.0, -1.0)}
    return ((n[1.0] + n[-1.0] - 2 * n[0.0])
            / (n[1.0] + n[-1.0] + n[0.0]))


def folded_asymmetry(sampler, plan, n_events, k_edges, optics, rng):
    """RP-accepted thirds estimator vs k, per category counts."""
    counts = {}
    for cat in plan.categories:
        share = cat.lumi_fraction
        ev = sampler.sample_category(cat, n=int(n_events * share), rng=rng)
        route = route_charged(ev["R"], ev["theta"], ev["pT"], optics)
        acc = (route == 1) | (route == 4)
        counts[cat.name] = np.histogram(ev["k"][acc], bins=k_edges)[0]
    n_tot = sum(counts.values())
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.where(n_tot > 0, est.azz_thirds(
            counts["azz+"], counts["azz-"], counts["azz0"], PZZ), np.nan)
    return a, n_tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", type=float, default=4e5,
                    help="generated events per pseudo-experiment")
    ap.add_argument("--lumi-ref", type=float, default=10.0,
                    help="reference lumi [fb^-1/u] for projected errors")
    ap.add_argument("--seed", type=int, default=20260713)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs("6Li")[1]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # --- left: analytic band ------------------------------------------------
    central = tagged.TaggedModel(tagged.li6_alpha_channel())
    ic = np.argmin(np.abs(central.c))
    ax1.plot(central.k, azz_wf_curve(central, ic), "k-", lw=2,
             label=r"$\beta=0.3$, $P_D=%.3f$" % tagged.P_D_LI6)
    for beta in (0.20, 0.40):
        m = tagged.TaggedModel(tagged.li6_alpha_channel(beta=beta))
        ax1.plot(m.k, azz_wf_curve(m, np.argmin(np.abs(m.c))), "--", lw=1,
                 label=r"$\beta=%.2f$" % beta)
    for p_d in (0.04, 0.13):
        m = tagged.TaggedModel(tagged.li6_alpha_channel(p_d=p_d))
        ax1.plot(m.k, azz_wf_curve(m, np.argmin(np.abs(m.c))), ":", lw=1,
                 label=r"$P_D=%.2f$" % p_d)
    ax1.set_xlabel(r"$k$ [GeV/$c$]")
    ax1.set_ylabel(r"$A_{zz}^{\rm tag}(k,\ \theta_k=90^\circ)$")
    ax1.set_title(r"analytic: $\alpha$-$d$ S/D interference")
    ax1.set_xlim(0, 0.6)
    ax1.legend(fontsize=7)
    ax1.axhline(0, color="0.6", lw=0.5)

    # --- right: acceptance-folded pseudo-experiment -------------------------
    kern = InclusiveKernel(beams.DEUTERON, b1_func=toy_b1)
    sampler = tagged.TaggedSampler(central, kern, config, fom.Scenario())
    plan = bk.tensor_thirds_plan(PZ, PZZ)
    k_edges = np.linspace(0.0, 0.6, 13)
    kc = 0.5 * (k_edges[:-1] + k_edges[1:])
    rng = np.random.default_rng(args.seed)
    summary = []
    for optics, fmt in ((HIGH_ACCEPTANCE, "o"), (HIGH_DIVERGENCE, "s")):
        a, n = folded_asymmetry(sampler, plan, args.events, k_edges,
                                optics, rng)
        errs = np.where(n > 1, err_azz(np.maximum(n, 1), PZZ), np.nan)
        use = n > 50
        ax2.errorbar(kc[use], a[use], yerr=errs[use], fmt=fmt, ms=4,
                     capsize=2, label=optics.name)
        # projected errors at the reference luminosity
        sigma_pb = sum(sampler.sigma_tot_pb(c) * c.lumi_fraction
                      for c in plan.categories)
        scale = (args.lumi_ref * 1e3 * sigma_pb) / args.events
        n_ref = n * scale
        i300 = np.argmin(np.abs(kc - 0.3))
        if n[i300] > 0:
            summary.append(
                "%s: accepted frac %.3f; at k~0.3 GeV A_zz = %+.2f, "
                "delta(10 fb^-1/u) = %.3f"
                % (optics.name, n.sum() / args.events, a[i300],
                   err_azz(max(n_ref[i300], 1), PZZ)))
    ax2.plot(central.k, azz_wf_curve(central, ic), "k-", lw=1, alpha=0.5,
             label=r"unfolded $90^\circ$ curve")
    ax2.set_xlabel(r"$k_{\rm rec}$ [GeV/$c$]")
    ax2.set_ylabel(r"$A_{zz}^{\rm tag}$ (RP-accepted)")
    ax2.set_title("acceptance-folded, tensor-thirds fills")
    ax2.set_xlim(0, 0.6)
    ax2.legend(fontsize=7)
    ax2.axhline(0, color="0.6", lw=0.5)

    fig.suptitle(r"$^6$Li($e,e^\prime\alpha$)X: tagged tensor asymmetry of "
                 r"the embedded deuteron, %s (TOY/scenario inputs)"
                 % config.label())
    fig.tight_layout()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "money_tagged_azz_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    for line in summary:
        print(line)


if __name__ == "__main__":
    main()
