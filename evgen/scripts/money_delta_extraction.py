#!/usr/bin/env python3
"""Money plot 7: the extracted Delta(x, Q2) structure function as
projected data points.

Money plot 5 shows the cos 2phi modulation and its amplitude; this
companion presents what that amplitude DELIVERS -- the double-helicity-
flip structure function itself, unfolded per (x, Q2) bin:

    A^cos2phi = -(1-y)/y^2 * Delta / D_phi   (D_phi = F1 + (1-y)/(xy^2) F2)
    =>  Delta_hat = -A_hat * y^2 D_phi / (1 - y)

implemented as model bin-centering: Delta_hat = A_hat * K with
K = Delta_model(x_c, Q2_c) / A_model(bin) -- identical to the pure
kinematic inversion at the bin center, plus the (few-%) in-bin-shape
correction an experiment would take from the same model it fits.
Statistical errors propagate as dDelta = dA * |K|; each luminosity gets
an INDEPENDENT pseudo-measurement (1-year and 10-year EIC programs).

Panels: Delta(x) at the three sweet-spot Q2 slices (the Q2-band picks of
money_cos2phi.py), with the injected moment_A curve and the conservative
moment_B alternative -- the extraction resolves the A-vs-B ansatz
spread bin by bin.  Delta is per nucleon with the 6Li counting dilution
1/3 included (2 of 6 nucleons); statistical only, backgrounds and
tensor RC unquantified (plans/06).

Usage:  python3 scripts/money_delta_extraction.py
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

from money_cos2phi import (add_tensor_leakage_args, b34_funcs,  # noqa: E402
                           build_delta_model,
                           check_tensor_leakage_args, leakage_amplitude,
                           measure, pick_sweet_spots_banded, superbin_mask,
                           tensor_leakage_tag, truth_leakage_route)

from polligen import bookkeeping as bk  # noqa: E402
from polligen.sample import InclusiveSampler  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, delta_models as dm, fom  # noqa: E402
PRETTY = {"moment_A": "moment-constrained $\\Delta$", "moment_B": "no-$F_1$ variant", "toy": "flat toy"}  # figure labels
from polli_fastsim.kinematics import kinematic_mask  # noqa: E402
from polli_fastsim.polarized import toy_b1  # noqa: E402

C_TRUTH, C_FIT, C_ALT = "#0072B2", "#D55E00", "#009E73"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2))
    ap.add_argument("--delta-model", default="moment_A",
                    choices=dm.available())
    ap.add_argument("--variant", default="mid_x",
                    choices=sorted(dm.VARIANTS))
    ap.add_argument("--dilution", type=float, default=1.0 / 3.0)
    ap.add_argument("--scale", type=float, default=1e-2,
                    help="peak Delta/F1 (toy model only)")
    ap.add_argument("--lumi-1yr", type=float, default=10.0)
    ap.add_argument("--lumi-10yr", type=float, default=100.0)
    ap.add_argument("--pzz", type=float, default=0.60)
    add_tensor_leakage_args(ap)
    ap.add_argument("--seed", type=int, default=20260811)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    check_tensor_leakage_args(ap, args)

    config = beams.default_configs("6Li")[args.config]
    lumi1_pb = args.lumi_1yr * 1e3
    lumi10_pb = args.lumi_10yr * 1e3
    rng = np.random.default_rng(args.seed)

    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    model, _q2_ref = build_delta_model(args, config, scenario)
    b3_func, b4_func = b34_funcs(args)
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1, delta_func=model,
                           b3_func=b3_func, b4_func=b4_func,
                           tensor_gamma=args.tensor_gamma)

    proj = fom.project_rates(config, scenario)
    obs = fom.project_observables(config, scenario, proj,
                                  kern.g1_model, toy_b1, model)
    spots = pick_sweet_spots_banded(proj, obs["sig_a_cos2phi"])
    q2_slices = sorted({round(q2, 3) for _x, q2, _i, _j in spots})

    sampler = InclusiveSampler(kern, config, scenario, nx=60, nq2=45)
    plan = bk.transverse_tensor_plan(args.pzz)
    cat = plan.categories[0]
    subtract = args.subtract_tensor_leakage != "none"

    alt = (dm.make("moment_B", variant=args.variant,
                   dilution=args.dilution)
           if args.delta_model != "moment_B" else
           dm.make("toy", scale=args.scale))
    alt_label = ("no-$F_1$ variant (conservative)"
                 if args.delta_model != "moment_B" else "flat toy")

    fig, axes = plt.subplots(1, len(q2_slices), figsize=(12.6, 4.5),
                             sharey=True)
    xe = proj.x_edges
    summary = []
    for ax, q2s in zip(np.atleast_1d(axes), q2_slices):
        pts = []
        for i0 in range(0, xe.size - 2, 2):  # merge pairs of x bins
            xc = np.sqrt(xe[i0] * xe[i0 + 2])
            if not kinematic_mask(xc, q2s, sampler.s):
                continue
            mask = superbin_mask(sampler, xe[i0], xe[i0 + 2],
                                 q2s / 1.6, q2s * 1.6)
            if not mask.any():
                continue
            m1 = measure(sampler, cat, mask, lumi1_pb, plan.pzz_true, rng)
            if m1["n"] < 1e3 or m1["err"] > 8e-3:
                continue
            if subtract:
                # additive on the amplitude, and on the truth K divides
                # by, so that K stays the pure Delta bin-centering and the
                # correction is applied once (plans/08 D2 risks R3, R8)
                leak = leakage_amplitude(sampler, cat, mask, m1["sigma_pb"])
                m1["amp"] -= leak
                m1["truth"] -= leak
            if abs(m1["truth"]) < 1e-5:
                continue
            m10 = measure(sampler, cat, mask, lumi10_pb, plan.pzz_true,
                          rng)
            if subtract:
                m10["amp"] -= leak
            # model bin-centering: K converts the bin-averaged amplitude
            # to Delta at the bin center (pure kinematic inversion +
            # in-bin-shape correction from the fitted model itself)
            f1c = kern.nf2.f1a(xc, q2s) / kern.ion.A
            delta_c = float(model(xc, q2s, f1c))
            k_conv = delta_c / m1["truth"]
            pts.append({
                "x": xc, "delta_c": delta_c,
                "d1": m1["amp"] * k_conv, "e1": m1["err"] * abs(k_conv),
                "d10": m10["amp"] * k_conv,
                "e10": m10["err"] * abs(k_conv)})
        xoff = 1.045
        ax.errorbar([p["x"] for p in pts],
                    [1e3 * p["x"] * p["d1"] for p in pts],
                    yerr=[1e3 * p["x"] * p["e1"] for p in pts], fmt="s",
                    mfc="none", color=C_FIT, ms=4, capsize=2, lw=1,
                    zorder=3,
                    label=r"1-yr EIC (%g fb$^{-1}$/u)" % args.lumi_1yr)
        ax.errorbar([p["x"] * xoff for p in pts],
                    [1e3 * p["x"] * p["d10"] for p in pts],
                    yerr=[1e3 * p["x"] * p["e10"] for p in pts], fmt="o",
                    color="black", ms=4, capsize=2, lw=1, zorder=4,
                    label=r"10-yr EIC (%g fb$^{-1}$/u)" % args.lumi_10yr)
        xg = np.logspace(np.log10(2e-4), np.log10(0.8), 250)
        q2g = np.full_like(xg, q2s)
        ok = kinematic_mask(xg, q2g, sampler.s)
        f1g = kern.nf2.f1a(xg, q2g) / kern.ion.A
        ax.plot(xg[ok], 1e3 * (xg * np.asarray(model(xg, q2g, f1g)))[ok],
                "-", color=C_TRUTH, lw=1.6,
                label="%s (injected)" % PRETTY.get(args.delta_model, args.delta_model))
        ax.plot(xg[ok], 1e3 * (xg * np.asarray(alt(xg, q2g, f1g)))[ok],
                "-", color=C_ALT, lw=1.4, label=alt_label)
        ax.set_xscale("log")
        ax.set_xlabel(r"$x$")
        ax.axhline(0, color="0.85", lw=0.6, zorder=0)
        ax.set_title(r"$Q^2 \approx %.3g$ GeV$^2$" % q2s, fontsize=9)
        ax.tick_params(labelsize=8)
        if pts:
            best = min(pts, key=lambda p: p["e10"])
            summary.append(
                "Q2=%.3g: %d points; best x=%.3g Delta=%+.4f "
                "+- %.4f (1yr) +- %.4f (10yr)"
                % (q2s, len(pts), best["x"], best["delta_c"],
                   best["e1"], best["e10"]))
    np.atleast_1d(axes)[0].set_ylabel(
        r"$x\,\Delta(x, Q^2)$ per nucleon  $[\times 10^{-3}]$",
        fontsize=9)
    np.atleast_1d(axes)[0].legend(fontsize=7, loc="lower left")

    fig.suptitle(
        r"Extracted double-helicity-flip $\Delta(x,Q^2)$ of $^6$Li, %s, "
        r"$P_{zz}=%.2f$""\n"
        r"$\hat\Delta = -\hat A\, y^2 D_\phi/(1-y)$ per bin "
        r"(model bin-centering; dilution 1/3 included); area = bag moment $-0.012\,\alpha_s/3$""\n"
        "statistical errors only; backgrounds and tensor radiative "
        "corrections not included"
        % (config.label(), plan.pzz_true), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    # a non-default tensor-leakage setting writes its own stem
    tag = tensor_leakage_tag(args)
    out = outdir / ("money_delta_extracted_6Li%s.png"
                    % ("_" + tag if tag else ""))
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print("delta model:", model.info())
    if args.tensor_gamma:
        print("tensor sector: EXACT finite-gamma kernel, b3 = %.3g b2, "
              "b4 = %.3g b2; leakage subtraction: %s"
              % (args.b3_frac, args.b4_frac, truth_leakage_route(args)))
    for line in summary:
        print(line)


if __name__ == "__main__":
    main()
