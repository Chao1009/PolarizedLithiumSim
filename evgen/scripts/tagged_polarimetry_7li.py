"""7Li alpha-tag bonus observables (plans/05 step 5.B deliverables):

Left: the in-situ alignment polarimeter -- the tagged-alpha polar
anisotropy <P2(cos theta_k)> = -T/5 vs the fill tensor moment T, with
Roman-Pot-folded pseudo-experiment markers (the 7Li alpha-tag has 96-99%
RP acceptance, so the folding barely dilutes the analyzing power).

Right: the tagged polarized-EMC companion -- A_par^tag(x) on the
quasi-free triton (helicity-flip fills, RP-accepted alpha tag), analytic
D(y) g1_t/F1_t overlay.  P_p(t) = 0.86 makes the tagged triton an
effective polarized proton at 86% polarization with an alpha in the RP
as the event tag.

Usage:  python3 scripts/tagged_polarimetry_7li.py --events 300000
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
from polligen.spin import spin32_populations  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.asymmetries import a_parallel, err_a_parallel  # noqa: E402
from polli_fastsim.farforward import (  # noqa: E402
    HIGH_ACCEPTANCE, HIGH_DIVERGENCE, route_charged)

PE, PZ = 0.7, 0.7


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", type=float, default=3e5)
    ap.add_argument("--seed", type=int, default=20260713)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs("7Li")[1]
    model = tagged.TaggedModel(tagged.li7_alpha_channel())
    kern = InclusiveKernel(tagged.TRITON)
    sampler = tagged.TaggedSampler(model, kern, config, fom.Scenario())
    rng = np.random.default_rng(args.seed)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # --- left: alignment polarimetry ----------------------------------------
    t_scan = np.linspace(-1.0, 1.0, 9)
    ax1.plot(t_scan, -t_scan / 5.0, "k-", lw=1.5,
             label=r"analytic $\langle P_2\rangle = -T/5$")
    n_per_point = int(args.events / t_scan.size)
    for optics, fmt in ((HIGH_ACCEPTANCE, "o"), (HIGH_DIVERGENCE, "s")):
        vals = []
        for t in t_scan:
            pops = spin32_populations(0.0, t, 0.0)
            cat = bk.SpinCategory("pol T=%.2f" % t, 1.5, pops)
            ev = sampler.sample_category(cat, n=n_per_point, rng=rng)
            route = route_charged(ev["R"], ev["theta"], ev["pT"], optics)
            acc = (route == 1) | (route == 4)
            c = ev["cos_theta_k"][acc]
            vals.append(np.mean(0.5 * (3 * c * c - 1.0)))
        ax1.plot(t_scan, vals, fmt, ms=4, label="RP-folded, " + optics.name)
    ax1.set_xlabel(r"fill tensor moment $T$")
    ax1.set_ylabel(r"$\langle P_2(\cos\theta_k)\rangle$ of the tagged "
                   r"$\alpha$")
    ax1.set_title("in-situ alignment polarimeter")
    ax1.legend(fontsize=8)
    ax1.axhline(0, color="0.6", lw=0.5)

    # --- right: tagged A_par on the quasi-free triton -----------------------
    plan = bk.helicity_flip_plan(1.5, PZ, PE)
    x_edges = np.logspace(-3, np.log10(0.7), 11)
    xc = np.sqrt(x_edges[:-1] * x_edges[1:])
    counts = {}
    for cat in plan.categories:
        ev = sampler.sample_category(
            cat, n=int(args.events * cat.lumi_fraction), rng=rng)
        route = route_charged(ev["R"], ev["theta"], ev["pT"],
                              HIGH_ACCEPTANCE)
        acc = (route == 1) | (route == 4)
        counts[cat.name] = np.histogram(ev["x"][acc], bins=x_edges)[0]
    n_tot = counts["apar+"] + counts["apar-"]
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.where(n_tot > 0, est.apar_flip(
            counts["apar+"], counts["apar-"], PE, PZ), np.nan)
    errs = np.where(n_tot > 1,
                    err_a_parallel(np.maximum(n_tot, 1), PE, PZ), np.nan)
    use = n_tot > 100
    ax2.errorbar(xc[use], a[use], yerr=errs[use], fmt="o", ms=4, capsize=2,
                 label="tagged pseudo-exp (RP-accepted)")
    # analytic overlay: sigma-weighted D g1t/F1t per x-bin from the inner grid
    inner = sampler.inner
    t = inner.tables
    y = inner.q2_cells / (inner.s * inner.x_cells)
    apar_cells = a_parallel(t["g1"], t["f1"], y, inner.x_cells,
                            inner.q2_cells)
    xb = np.digitize(inner.x_cells, x_edges) - 1
    analytic = np.array([
        np.average(apar_cells[xb == b], weights=inner.xsec_flat[xb == b])
        if np.any(xb == b) else np.nan for b in range(xc.size)])
    ax2.plot(xc, analytic, "k-", lw=1.5, label=r"$D(y)\,g_1^t/F_1^t$"
             r" ($\sigma$-weighted)")
    ax2.set_xscale("log")
    ax2.set_xlabel("x")
    ax2.set_ylabel(r"$A_\parallel^{\rm tag}$ (quasi-free $t$)")
    ax2.set_title(r"tagged polarized-EMC companion, $P_p(t)=0.86$")
    ax2.legend(fontsize=8)
    ax2.axhline(0, color="0.6", lw=0.5)

    fig.suptitle(r"$^7$Li($e,e^\prime\alpha$)X, %s (TOY inputs): "
                 r"$\alpha$-tag acceptance makes both observables "
                 "essentially free" % config.label())
    fig.tight_layout()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "tagged_polarimetry_7Li.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)

    # headline numbers
    acc_frac = {}
    cat = bk.SpinCategory("acc", 1.5, (0.25, 0.25, 0.25, 0.25))
    ev = sampler.sample_category(cat, n=100_000, rng=rng)
    for optics in (HIGH_ACCEPTANCE, HIGH_DIVERGENCE):
        route = route_charged(ev["R"], ev["theta"], ev["pT"], optics)
        acc_frac[optics.name] = float(np.mean((route == 1) | (route == 4)))
    print("7Li alpha-tag RP acceptance:", acc_frac)


if __name__ == "__main__":
    main()
