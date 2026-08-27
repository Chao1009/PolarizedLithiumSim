#!/usr/bin/env python3
"""How much charge information the 6Li/alpha separation actually needs.

reports/nanowire_far_forward (2026-08-26) argued that a superconducting
nanowire loses open question #19 to the incumbent because the nanowire
delivers ONE BIT per plane while the ePIC AC-LGAD digitises an 8-bit
charge over a 30 um active layer, quoting "4.8 sigma per plane".  Both
halves of that were wrong, and this script is the correction.

  * A SIGMA IS THE WRONG FIGURE OF MERIT.  "Gap over the quadrature sum
    of two Landau core widths" is neither a PID separation power nor a
    fake rate.  What puts an alpha inside a 6Li's window is the Landau
    UPPER TAIL, which falls as 1/lambda -- nothing like a Gaussian.  The
    honest statement is a fake rate at a stated efficiency.

  * ONE BIT IS ALMOST ENOUGH.  Against the Neyman-Pearson optimum (the
    per-plane Landau log-likelihood ratio), a one-bit-per-plane threshold
    with a majority-of-k decision loses only a small factor, because the
    two species are far apart (MPV 31.7 vs 75.2 keV) and the power comes
    from requiring COINCIDENCE across planes rather than from precision
    within one.  A truncated mean -- the standard analogue dE/dx
    estimator -- is two orders of magnitude WORSE than one bit, and a
    plain sum is worse still, because a Landau has no mean and one delta
    ray drags it.

  * WHERE THE NANOWIRE ACTUALLY LOSES IS GEOMETRIC FILL FACTOR.  The
    coincidence that makes one bit sufficient needs every plane to record
    the track.  A silicon pixel plane does so ~99% of the time; a
    superconducting wire comb does so only over its fill factor -- 25%
    (arXiv:2510.11725), 40% (arXiv:2410.00251), 50% (the ANL EIC-targeted
    device, arXiv:2312.13405).  At 50% fill a four-plane array cannot
    reach 95% 6Li efficiency at any working point.

That is a fabrication number, not an information-theoretic one, which
makes it a fairer thing to put to the group than "you only have one bit".

Usage:  python3 scripts/nearbeam_zid_power.py [--outdir .]
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from polligen import nearbeam as nb  # noqa: E402

C_OPT, C_BIT, C_TRUNC, C_SUM = "#1B7F5B", "#1F4E79", "#B8860B", "#C0392B"

DEVICES = (
    ("AC-LGAD, 8-bit charge (LLR)",   "llr",       0.99, C_OPT),
    ("one bit/plane, 99% efficient",  "threshold", 0.99, C_BIT),
    ("AC-LGAD, truncated mean",       "trunc",     0.99, C_TRUNC),
    ("AC-LGAD, plain sum",            "sum",       0.99, C_SUM),
)
FILLS = ((0.50, "ANL, 50% fill"), (0.40, "40% fill"), (0.25, "SMSPD, 25% fill"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-mc", type=int, default=1500000)
    ap.add_argument("--efficiency", type=float, default=0.95)
    ap.add_argument("--n-planes", type=int, default=4)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    floor = 1.0 / args.n_mc

    print("== alpha faking 6Li, %d planes of %s, matched %.0f%% efficiency =="
          % (args.n_planes, nb.SI_ACLGAD.name, 100 * args.efficiency))
    print("MC = %.1e, so %.1e is the statistical floor\n" % (args.n_mc, floor))
    print("%-32s %-14s %s" % ("device / readout", "alpha fake", "eff reached"))
    res = {}
    for lbl, ro, pe, _ in DEVICES:
        f, e = nb.zid_fake_rate(3, 2, efficiency=args.efficiency, readout=ro,
                                n_planes=args.n_planes, n_mc=args.n_mc,
                                plane_efficiency=pe,
                                rng=np.random.default_rng(args.seed))
        res[lbl] = (f, e)
        print("%-32s %-14s %.3f"
              % (lbl, ("%.2e" % f) if f == f else "cannot reach", e))
    print()
    for fill, lbl in FILLS:
        f, e = nb.zid_fake_rate(3, 2, efficiency=args.efficiency,
                                readout="threshold", n_planes=args.n_planes,
                                n_mc=args.n_mc, plane_efficiency=fill,
                                rng=np.random.default_rng(args.seed))
        res[lbl] = (f, e)
        print("%-32s %-14s %.3f"
              % ("nanowire, " + lbl, ("%.2e" % f) if f == f else "cannot reach", e))

    opt = res["AC-LGAD, 8-bit charge (LLR)"][0]
    bit = res["one bit/plane, 99% efficient"][0]
    print("\nONE BIT costs a factor %.1f against the Neyman-Pearson optimum."
          % (bit / opt))
    print("A TRUNCATED MEAN -- the standard analogue estimator -- costs a "
          "factor %.0f AGAINST ONE BIT."
          % (res["AC-LGAD, truncated mean"][0] / bit))
    print("The nanowire's loss is not bits: at 50%% fill a %d-plane array "
          "tops out at %.1f%% 6Li efficiency."
          % (args.n_planes, 100 * res["ANL, 50% fill"][1]))

    # --- efficiency reach versus fill factor, which is the real story ----
    fills = np.linspace(0.20, 1.0, 17)
    reach = []
    for fl in fills:
        _, e = nb.zid_fake_rate(3, 2, efficiency=0.999, readout="threshold",
                                n_planes=args.n_planes, n_mc=200000,
                                plane_efficiency=fl,
                                rng=np.random.default_rng(args.seed))
        reach.append(e)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(100 * fills, 100 * np.array(reach), "-", color=C_BIT, lw=1.8)
    for fill, lbl in FILLS:
        _, e = nb.zid_fake_rate(3, 2, efficiency=0.999, readout="threshold",
                                n_planes=args.n_planes, n_mc=200000,
                                plane_efficiency=fill,
                                rng=np.random.default_rng(args.seed))
        ax.plot([100 * fill], [100 * e], "o", color=C_BIT, ms=6, mfc="white",
                mew=1.6)
        ax.annotate(lbl, xy=(100 * fill, 100 * e), xytext=(6, -10),
                    textcoords="offset points", fontsize=7.5)
    ax.axhline(100 * args.efficiency, color="0.45", ls="--", lw=1.0)
    ax.annotate("%.0f%% required" % (100 * args.efficiency),
                xy=(100, 100 * args.efficiency), xytext=(-4, 4),
                textcoords="offset points", fontsize=7.5, ha="right",
                color="0.35")
    ax.plot([99], [100 * (1 - (1 - 0.99) ** 1)], "s", color=C_OPT, ms=6)
    ax.annotate("silicon pixel plane", xy=(99, 99), xytext=(-6, -12),
                textcoords="offset points", fontsize=7.5, ha="right",
                color=C_OPT)
    ax.set_xlabel("per-plane efficiency (geometric fill factor) [%]")
    ax.set_ylabel(r"best reachable $^{6}$Li efficiency, %d planes [%%]"
                  % args.n_planes)
    ax.set_title("Why the nanowire loses open question #19: fill factor, "
                 "not bits", fontsize=10)
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    out = pathlib.Path(args.outdir) / "nearbeam_zid_power.png"
    fig.savefig(out, dpi=140)
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
