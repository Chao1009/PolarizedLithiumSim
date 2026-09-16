#!/usr/bin/env python3
"""Money plot 8R: the RATE tensor and vector asymmetries at the
RECONSTRUCTED level (plans/03 2.4).

The cos 2phi' chain of `money_cos2phi_reco.py` measures an AZIMUTHAL
harmonic; this one measures the two RATE asymmetries of the same run
plan through the same response:

  A_zz   tensor thirds (`bookkeeping.tensor_thirds_plan`,
         `estimators.azz_thirds`): three fills at (P_z, +P_zz),
         (-P_z, +P_zz) and the m0-enriched (0, -2 P_zz), and
         (y+ + y- - 2 y0)/(y+ + y- + y0)/P_zz;
  A_par  helicity flips (`bookkeeping.helicity_flip_plan`,
         `estimators.apar_flip`): one longitudinally vector-polarized
         fill read with both electron helicities, and
         (y+ - y-)/(y+ + y-)/(P_e P_z).

Both run on `recopseudo.measure_azz` / `measure_apar`, i.e. on the exact
expected COUNTS of each fill in a RECONSTRUCTED (x, Q2) bin -- the same
migration, the same reconstructed-level selection and the same
eps_eID(eta) as 5R and 7R -- Poisson-fluctuated at the programme
luminosity.  Each point is plotted against two references: the
reconstructed-bin truth (what the measurement converges to) and the
true-bin truth (what a perfect detector would give), whose ratio is the
migration the bin carries.

THE POLARIMETRY BAND.  `--delta-p-over-p` threads `RunPlan.delta_p_over_p`
-- the polarimetry scale of plans/04 #5, 3% in the ring -- into the
figure.  Both estimators divide by a MEASURED polarization, so an error
delta on it scales every point by 1/(1 + delta) and nothing else: the
band is a pure scale, fully correlated across x.  Each edge is drawn by
RE-READING the estimator on the same counts with the polarization it
divides by moved to P (1 +- delta), not by multiplying the central value
by delta, so the printed half-width is a measurement of that claim and
not a restatement of it; it comes out at delta/(1 - delta^2) of the
central value, which is delta to 0.1% at the 3% ring value.  For A_par
the divisor is the product P_e P_z and delta is the scale on it.  OFF by
default, and a run with it on writes its own PNG stem, so the published
figure never carries it.

Outputs:
  money_azz_reco_6Li.png     (A_zz in reconstructed x bins, one Q2 slice)
  money_apar_reco_6Li.png    (--observable apar)

Usage:  python3 scripts/money_azz_reco.py [--observable {azz,apar}]
                                          [--delta-p-over-p 0.03]
"""

import argparse
import pathlib
import sys

import numpy as np

_SCRIPTS = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS.parent))
sys.path.insert(0, str(_SCRIPTS))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from polligen import bookkeeping as bk  # noqa: E402
from polligen import estimators as est  # noqa: E402
from polligen import recopseudo as rp  # noqa: E402
from polligen.sample import InclusiveSampler  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.polarized import toy_b1  # noqa: E402

C_TRUTH, C_FIT, C_ALT = "#0072B2", "#D55E00", "#009E73"
C_BAND = "#9E9E9E"


def output_stem(base, args):
    """File stem of one run.  The published artefact is the DEFAULT
    combination: the middle configuration at the full programme share
    with no polarimetry band.  Every other setting appends its key rather
    than overwriting it -- the guard convention of
    `money_tagged_azz.output_stem` and of
    `money_cos2phi_reco.reco_setting_keys`."""
    keys = []
    if args.config != 1:
        keys.append("cfg%d" % args.config)
    share_key = fom.run_share_tag(args.lumi_fraction)
    if share_key:
        keys.append(share_key)
    if args.delta_p_over_p:
        keys.append(("pol%g" % args.delta_p_over_p).replace(".", "p"))
    return "_".join([base] + keys)


def x_slice_bins(resp, q2lo, q2hi, x_edges, plan, lumi_pb, n_min=1.0e4):
    """The reconstructed x bins of one Q2 slice that hold enough events."""
    out = []
    for xlo, xhi in zip(x_edges[:-1], x_edges[1:]):
        mask = resp.mask_reco(xlo, xhi, q2lo, q2hi)
        if not mask.any():
            continue
        n = resp.expected_rates(plan.categories, lumi_pb, mask).sum()
        if n < n_min:
            continue
        out.append((xlo, xhi, mask, float(n)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--observable", default="azz", choices=("azz", "apar"),
                    help="azz: tensor thirds; apar: electron-helicity flips")
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2))
    ap.add_argument("--lumi-1yr", type=float, default=10.0)
    ap.add_argument("--lumi-10yr", type=float, default=100.0)
    ap.add_argument("--lumi-fraction", type=float, default=1.0,
                    dest="lumi_fraction",
                    help="this observable's share of the PROGRAMME "
                         "luminosity (plans/07 WP2; default 1.0, which "
                         "every published number assumes)")
    ap.add_argument("--pzz", type=float, default=0.60)
    ap.add_argument("--pz", type=float, default=0.70)
    ap.add_argument("--pe", type=float, default=0.70)
    ap.add_argument("--delta-p-over-p", type=float, default=0.0,
                    dest="delta_p_over_p",
                    help="polarimetry scale error (plans/04 #5; 0.03 in "
                         "the ring).  Both estimators divide by a MEASURED "
                         "polarization (the product P_e P_z for A_par), "
                         "so this is a pure scale on every point.  The band "
                         "edges are the estimator re-read at P(1 +- delta), "
                         "and the half-width they give is measured against "
                         "delta afterwards.  OFF by default; a run with "
                         "it on writes its own PNG stem")
    ap.add_argument("--q2-slice", type=float, nargs=2, default=(3.0, 10.0),
                    dest="q2_slice", metavar=("Q2LO", "Q2HI"))
    ap.add_argument("--nx", type=int, default=8,
                    help="reconstructed x bins across the slice")
    ap.add_argument("--x-range", type=float, nargs=2, default=(0.01, 0.4),
                    dest="x_range", metavar=("XLO", "XHI"))
    ap.add_argument("--n-mc-per-cell", type=int, default=400)
    ap.add_argument("--rel-lumi-offset", type=float, default=0.0,
                    help="relative-luminosity error on the m0-enriched "
                         "(A_zz) or lam_e = +1 (A_par) share, unknown to "
                         "the analysis")
    ap.add_argument("--n-min", type=float, default=1.0e4,
                    help="minimum expected counts for a bin to be drawn")
    ap.add_argument("--seed", type=int, default=20260915)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    if not args.lumi_fraction > 0:
        ap.error("--lumi-fraction must be positive")
    config = beams.default_configs("6Li")[args.config]
    lumi1_pb = args.lumi_1yr * args.lumi_fraction * 1e3
    lumi10_pb = args.lumi_10yr * args.lumi_fraction * 1e3
    rng = np.random.default_rng(args.seed)
    analysis = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            run_share=args.lumi_fraction,
                            pol_ion_tensor=args.pzz,
                            pol_ion_vector=args.pz, pol_electron=args.pe)
    print("run plan: programme %g / %g fb^-1/u (1 yr / 10 yr) x share %g "
          "-> %g / %g fb^-1/u delivered to this observable"
          % (args.lumi_1yr, args.lumi_10yr, args.lumi_fraction,
             args.lumi_1yr * args.lumi_fraction,
             args.lumi_10yr * args.lumi_fraction))

    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    gen = rp.generator_scenario(analysis)
    sampler = InclusiveSampler(kern, config, gen, nx=60, nq2=45,
                               q2_range=(0.7, 2e3))
    resp = rp.RecoResponse(sampler, rp.RecoModel(
        q2_min=analysis.q2_min, y_min=analysis.y_min, y_max=analysis.y_max,
        w2_min=analysis.w2_min, eta_min=analysis.eta_min,
        eta_max=analysis.eta_max, e_prime_min=analysis.e_prime_min),
        n_mc_per_cell=args.n_mc_per_cell, rng=rng)
    print("response: %d MC events, generator sigma = %.4g pb, "
          "selected+eID = %.4g pb"
          % (resp.x.size, resp.w.sum(), (resp.w * resp.eff).sum()))

    if args.observable == "azz":
        plan = bk.tensor_thirds_plan(args.pz, args.pzz,
                                     rel_lumi_offset=args.rel_lumi_offset)
        measure = rp.measure_azz
        nominal = [1.0 / 3.0] * 3

        def rescaled(m, f):
            """The estimate the SAME counts give when the analysis divides
            by a polarization that is wrong by the factor f."""
            c = np.asarray(m["counts"], dtype=float)
            return float(est.azz_thirds(c[0], c[1], c[2], m["pzz"] * f,
                                        lumis=m["lumi_fractions"]))

        label = r"$A_{zz}$"
        pol_label = r"$P_{zz} = %.2f$" % args.pzz
        base_stem = "money_azz_reco"
    else:
        plan = bk.helicity_flip_plan(1.0, args.pz, args.pe,
                                     rel_lumi_offset=args.rel_lumi_offset)
        measure = rp.measure_apar
        nominal = [0.5, 0.5]

        def rescaled(m, f):
            """The same, the divisor here being the product P_e P_z."""
            c = np.asarray(m["counts"], dtype=float)
            lum = m["lumi_fractions"]
            return float(est.apar_flip(c[0], c[1], m["pe"] * f, m["pz"],
                                       l_plus=lum[0], l_minus=lum[1]))

        label = r"$A_{\parallel}$"
        pol_label = r"$P_e P_z = %.2f$" % (args.pe * args.pz)
        base_stem = "money_apar_reco"

    q2lo, q2hi = args.q2_slice
    x_edges = np.geomspace(args.x_range[0], args.x_range[1], args.nx + 1)
    bins = x_slice_bins(resp, q2lo, q2hi, x_edges, plan, lumi1_pb,
                        n_min=args.n_min)
    if not bins:
        raise SystemExit("no reconstructed x bin holds %g events at this "
                         "luminosity" % args.n_min)

    xc, val, err1, err10, ref_reco, ref_true, nn = [], [], [], [], [], [], []
    meas = []
    for xlo, xhi, mask, _n in bins:
        true_mask = resp.mask_true(xlo, xhi, q2lo, q2hi)
        one = measure(resp, plan, lumi1_pb, mask, rng=rng,
                      lumi_assumed=nominal, true_mask=true_mask)
        ten = measure(resp, plan, lumi10_pb, mask, rng=rng,
                      lumi_assumed=nominal)
        meas.append(one)
        xc.append(np.sqrt(xlo * xhi))
        val.append(one["value"])
        err1.append(one["err"])
        err10.append(ten["err"])
        ref_reco.append(one["truth_reco_bin"])
        ref_true.append(one["truth_true_bin"])
        nn.append(one["n"])
    xc = np.array(xc)
    val, err1, err10 = np.array(val), np.array(err1), np.array(err10)
    ref_reco, ref_true, nn = (np.array(ref_reco), np.array(ref_true),
                              np.array(nn))

    print("%s at the reconstructed level, Q2 = %g-%g GeV^2, %s"
          % (args.observable.upper(), q2lo, q2hi, pol_label))
    for i, x in enumerate(xc):
        print("  x = %-7.4g N(1yr) = %.3e  hat = %+.5e +- %.2e (1 yr) "
              "+- %.2e (10 yr)  reco-bin truth %+.5e  true-bin %+.5e "
              "(migration %+.3f sigma_1yr)  pull %+.2f"
              % (x, nn[i], val[i], err1[i], err10[i], ref_reco[i],
                 ref_true[i], (ref_reco[i] - ref_true[i]) / err1[i],
                 (val[i] - ref_reco[i]) / err1[i]))
    pull = (val - ref_reco) / err1
    print("  %d bins: pull mean %+.3f, sd %.3f" % (pull.size, pull.mean(),
                                                   pull.std(ddof=1)))

    # --- the polarimetry scale band --------------------------------------
    band = None
    if args.delta_p_over_p:
        # The band is NOT drawn as delta x the central value -- that would
        # assume the very thing the work package asks to be checked.  Each
        # edge is the SAME drawn counts re-read by the SAME estimator with
        # the polarization it divides by moved to P (1 +- delta), which is
        # what a polarimetry error does to a rate asymmetry; the half-width
        # the figure then carries is measured against delta afterwards.
        # RunPlan draws the MEASURED polarization once per fill from its own
        # fixed-seed stream: that draw is one realization inside this
        # envelope, and re-reading the estimator at it must scale every bin
        # by the same 1/(1 + drawn) if the band is to be a pure scale at all.
        banded = bk.RunPlan(plan.categories, pe_true=plan.pe_true,
                            pz_true=plan.pz_true, pzz_true=plan.pzz_true,
                            delta_p_over_p=args.delta_p_over_p,
                            polarimetry_seed=args.seed)
        key = "pzz" if args.observable == "azz" else "pz"
        drawn = banded.measured[key] / (plan.pzz_true if key == "pzz"
                                        else plan.pz_true) - 1.0
        d = args.delta_p_over_p
        hi = np.array([rescaled(m, 1.0 + d) for m in meas])
        lo = np.array([rescaled(m, 1.0 - d) for m in meas])
        band = (np.minimum(hi, lo), np.maximum(hi, lo))
        half = 0.5 * (band[1] - band[0])
        frac = half / np.abs(np.where(val == 0.0, np.nan, val))
        worst = float(np.nanmax(np.abs(frac / d - 1.0)))
        # the plan's own draw, read through the same estimator
        moved = np.array([rescaled(m, 1.0 + drawn) for m in meas])
        scale = moved / np.where(val == 0.0, np.nan, val)
        pure = float(np.nanmax(np.abs(scale * (1.0 + drawn) - 1.0)))
        print("polarimetry scale: delta P/P = %.4g; the +-1 sigma envelope "
              "is the estimator re-read at P(1 +- delta), and its half-width "
              "/ |central| is %.6g at every one of the %d bins -- delta/"
              "(1 - delta^2), i.e. delta P/P to %.2g relative (gate 1%%)"
              % (d, float(np.nanmax(frac)), int(np.isfinite(frac).sum()),
                 worst))
        print("  this plan's own draw is %+.4g (%.2f sigma) and scales every "
              "bin by %+.4g%%, which is 1/(1 + draw) to %.2g relative"
              % (drawn, drawn / d, 100.0 * (np.nanmean(scale) - 1.0), pure))
        if worst > 0.01 or pure > 1e-9:
            raise SystemExit("polarimetry band is not a pure scale")

    # --- the figure --------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    if band is not None:
        ax.fill_between(xc, band[0], band[1], color=C_BAND, alpha=0.35,
                        lw=0, zorder=1,
                        label=r"polarimetry scale $\pm%.0f\%%$"
                              % (100.0 * args.delta_p_over_p))
    ax.plot(xc, ref_true, color=C_ALT, ls="--", lw=1.4, zorder=2,
            label="true-bin truth")
    ax.plot(xc, ref_reco, color=C_TRUTH, lw=1.8, zorder=3,
            label="reconstructed-bin truth")
    ax.errorbar(xc, val, yerr=err1, fmt="o", ms=5.0, color=C_FIT,
                capsize=2.5, zorder=4, label="1 year")
    ax.errorbar(xc, val, yerr=err10, fmt="none", ecolor="#7A2E00",
                capsize=4.5, elinewidth=2.4, zorder=5, label="10 years")
    ax.set_xscale("log")
    # a log axis over less than a decade and a half labels one decade and
    # nothing else by default, which leaves the points unreadable
    ticks = [v for v in (0.01, 0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3,
                         0.4) if xc.min() / 1.3 <= v <= xc.max() * 1.3]
    ax.set_xticks(ticks, minor=False)
    ax.set_xticklabels(["%g" % v for v in ticks])
    ax.set_xticks([], minor=True)
    ax.set_xlabel("reconstructed $x$")
    ax.set_ylabel(label)
    ax.axhline(0.0, color="0.6", lw=0.8, zorder=0)
    ax.set_title("%s at the reconstructed level — %s, $Q^2$ = %g–%g GeV$^2$, "
                 "%s" % (label, config.label(), q2lo, q2hi, pol_label),
                 fontsize=10.5)
    ax.legend(fontsize=8.5, loc="best", framealpha=0.9)
    ax.grid(alpha=0.25, lw=0.6)
    fig.tight_layout()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / ("%s_6Li.png" % output_stem(base_stem, args))
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
