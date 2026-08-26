#!/usr/bin/env python3
"""Money plots 5R and 7R: the inclusive gluonometry measurement at the
RECONSTRUCTED level (plans/07 WP3), re-deriving money plots 5 and 7 with

  * reconstructed (x, Q2) bins: Q2 from the scattered electron (EMCal
    2%/sqrt(E) (+) 1% + track angles), y from the hadronic final state
    (parametrized Sigma/JB resolution, default 25% = the ePIC kinematic-
    fit study's smearing and the ATHENA Fig. 22 value at y ~ 0.01) -> x by
    the mixed method; reco-level analysis cuts; eps_eID(eta); the 25 mrad crossing
    angle undone by the head-on transformation; the physics azimuth from
    the covariant formula on the smeared four-vectors;
  * the spin-state-sorted ratio estimator (reco.harmonic_ratio_fit) on
    m = +-1-rich (P_zz = +0.6) and m = 0-rich (-1.2) fills instead of the
    single-fill fit -- with a smooth phi'-dependent efficiency
    (1 + 0.03 cos 2phi' + 0.02 cos phi', aligned with the spin axis) and
    an unknown relative-luminosity offset switched ON to show that both
    cancel;
  * Delta from the amplitude with an MC bin-centering factor that
    includes the migration (recopseudo.delta_from_amplitude).

Statistics are exact at any luminosity (expected counts per phi' bin per
spin state, Poisson-fluctuated) on an importance-sampled response with
n_mc pseudo-events per sampler cell.  Outputs:

  money_cos2phi_reco_6Li.png        (5R: phi' pseudo-data in the four
                                     sweet-spot super-bins + amplitude vs x)
  money_delta_extracted_reco_6Li.png (7R: x Delta(x, Q2) points at the
                                     three sweet-spot Q2 slices)

Usage:  python3 scripts/money_cos2phi_reco.py [--y-method electron]
"""

import argparse
from dataclasses import replace
import pathlib
import sys

import numpy as np

_SCRIPTS = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS.parent))
sys.path.insert(0, str(_SCRIPTS))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

from money_cos2phi import (build_delta_model, pick_sweet_spots_banded,  # noqa: E402
                           superbin_edges)

from polligen import bookkeeping as bk  # noqa: E402
from polligen import hfs, reco, recopseudo as rp  # noqa: E402
from polligen.sample import InclusiveSampler  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, delta_models as dm, fom  # noqa: E402
from polli_fastsim.asymmetries import a_cos2phi  # noqa: E402
from polli_fastsim.kinematics import kinematic_mask, y_from_xq2  # noqa: E402
from polli_fastsim.polarized import toy_b1  # noqa: E402

C_TRUTH, C_FIT, C_ALT = "#0072B2", "#D55E00", "#009E73"


def modulation_points(fit):
    """Acceptance-free modulation T(phi') - kappa from the fitted
    pseudo-data, with errors (for display)."""
    r, var, sig2, pbar = reco.spin_state_ratio(fit["counts"],
                                               fit["lumi_fractions"],
                                               fit["pzz"])
    t, var_t = reco._ratio_to_modulation(r, var, sig2, pbar)
    return t - fit["const"], np.sqrt(var_t)


def measure_bin(resp, plan, mask, lumi_pb, rng, phi_eff, lumi_assumed):
    return rp.measure_inclusive(resp, plan, lumi_pb, mask, rng=rng,
                                phi_eff=phi_eff, lumi_assumed=lumi_assumed)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2))
    ap.add_argument("--delta-model", default="moment_A",
                    choices=dm.available())
    ap.add_argument("--variant", default="mid_x",
                    choices=sorted(dm.VARIANTS))
    ap.add_argument("--dilution", type=float, default=1.0 / 3.0)
    ap.add_argument("--scale", type=float, default=1e-2)
    ap.add_argument("--lumi-1yr", type=float, default=10.0)
    ap.add_argument("--lumi-10yr", type=float, default=100.0)
    ap.add_argument("--pzz", type=float, default=0.60)
    ap.add_argument("--y-method", default="mixed", choices=("mixed", "electron"))
    ap.add_argument("--y-had-res", type=float, default=0.25,
                    help="hadronic-method dy/y: 0.25 = the ePIC kinematic-"
                         "fit study's smearing of delta_h and the ATHENA "
                         "Fig. 22 value at y ~ 0.01 (band 0.15-0.30; "
                         "refs/README.md)")
    ap.add_argument("--energy", default="emcal",
                    choices=("emcal", "tracking", "best"))
    ap.add_argument("--y-source", default="param", choices=("param", "hfs"),
                    help="hadronic y from the Gaussian stand-in (param) or "
                         "from a hadronic-final-state library through the "
                         "hadron-side detector response (hfs)")
    ap.add_argument("--hfs-sample", nargs="*", default=None,
                    help="HFS sample file(s) (.npz from tools/pythia8/"
                         "gen_dis_hfs.py); default: toy string fragmentation")
    ap.add_argument("--hfs-method", default="sigma", choices=("sigma", "jb", "da"))
    ap.add_argument("--hfs-noise", type=float, default=0.05,
                    help="calorimeter noise on Sigma and each pT component [GeV]")
    ap.add_argument("--hfs-library-events", type=int, default=300000,
                    help="toy library size when no sample is given")
    ap.add_argument("--tag", default=None,
                    help="output-name suffix (default: '_hfs' for --y-source hfs)")
    ap.add_argument("--n-mc-per-cell", type=int, default=400)
    ap.add_argument("--eff-cos2", type=float, default=0.03,
                    help="cos 2phi' harmonic of the phi' efficiency")
    ap.add_argument("--eff-cos1", type=float, default=0.02)
    ap.add_argument("--eff-cos2-split", type=float, default=0.0,
                    help="DIFFERENCE of the cos 2phi' efficiency harmonic "
                         "between the m=+-1-rich and m=0-rich fills.  A "
                         "common harmonic cancels in the spin-state ratio; "
                         "a difference does not, and fakes an amplitude "
                         "(e_+ - e_0)/(P_+ - P_0) (code review F1).  "
                         "0 = common (the default of money plot 5R)")
    ap.add_argument("--rel-lumi-offset", type=float, default=1e-3,
                    help="relative-luminosity error on the +Pzz fill, "
                         "unknown to the analysis")
    ap.add_argument("--e-scale", type=float, default=1.0,
                    help="electron energy-scale calibration error unknown "
                         "to the analysis.  The bigger of the two energy "
                         "levers: d ln x / d ln E' = 2 - y with the Sigma "
                         "method (only 1 with the Gaussian y stand-in, "
                         "which never sees E')")
    ap.add_argument("--hfs-scale", type=float, default=1.0,
                    help="hadronic energy-scale calibration error "
                         "(--y-source hfs only): d ln x / d ln scale "
                         "= -(1 - y)")
    ap.add_argument("--eid-tilt", type=float, default=0.0,
                    help="linear eta slope on eps_eID.  A FLAT eps_eID "
                         "error is identically null in the ratio; only a "
                         "shape moves a number")
    ap.add_argument("--emcal-eta-table", action="store_true",
                    help="Yellow Report EMCal resolution per eta region "
                         "instead of the backward-endcap specification "
                         "everywhere (code review F4).  The sweet-spot "
                         "electrons are backward, so the headline numbers "
                         "do not move; the amplitude-vs-x panels reach the "
                         "barrel, where --energy best picks the tracker")
    ap.add_argument("--syst-scan", action="store_true",
                    help="after the money plots, rebuild the response with "
                         "each nuisance varied and print the Delta shift "
                         "per sweet spot (common random numbers, so the "
                         "shifts are not seed noise)")
    ap.add_argument("--seed", type=int, default=20260824)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs("6Li")[args.config]
    lumi1_pb, lumi10_pb = args.lumi_1yr * 1e3, args.lumi_10yr * 1e3
    rng = np.random.default_rng(args.seed)
    analysis = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    model, q2_ref = build_delta_model(args, config, analysis)
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1, delta_func=model)

    # analysis grid, sweet spots (identical selection to money plot 5)
    proj = fom.project_rates(config, analysis)
    obs = fom.project_observables(config, analysis, proj, kern.g1_model,
                                  toy_b1, model)
    spots = pick_sweet_spots_banded(proj, obs["sig_a_cos2phi"])[:4]
    s = config.sqrt_s_per_nucleon ** 2

    # generator sample (loose cuts) + reconstructed-level response
    gen = rp.generator_scenario(analysis)
    sampler = InclusiveSampler(kern, config, gen, nx=60, nq2=45,
                               q2_range=(0.7, 2e3))
    rmodel = rp.RecoModel(energy=args.energy, y_method=args.y_method,
                          y_source=args.y_source, hadronic_method=args.hfs_method,
                          y_had_res=args.y_had_res, e_scale=args.e_scale,
                          eid_tilt=args.eid_tilt,
                          emcal_eta_table=args.emcal_eta_table,
                          q2_min=analysis.q2_min, y_min=analysis.y_min,
                          y_max=analysis.y_max, w2_min=analysis.w2_min,
                          eta_min=analysis.eta_min, eta_max=analysis.eta_max,
                          e_prime_min=analysis.e_prime_min)
    hfs_resp = None
    if args.y_source == "hfs":
        hresp = hfs.HadronResponse(noise_sigma=args.hfs_noise)
        if args.hfs_sample:
            smp = hfs.HFSSample.concatenate([hfs.HFSSample.load(f)
                                             for f in args.hfs_sample])
            src = "%s (%d events)" % (smp.meta.get("generator", "sample"),
                                      smp.n_events)
        else:
            smp = hfs.toy_library_sample(sampler, args.hfs_library_events, rng)
            src = "TOY string fragmentation (%d events)" % smp.n_events
        lib = hfs.HFSLibrary(smp, hresp, nx=48, nq2=36, rng=rng)
        hfs_resp = hfs.HFSResponse(lib, method=args.hfs_method,
                                   scale=args.hfs_scale)
        print("HFS library: %s; %s method; %s | coverage %s"
              % (src, args.hfs_method.upper(), hresp.describe(), lib.coverage()))
    resp = rp.RecoResponse(sampler, rmodel, n_mc_per_cell=args.n_mc_per_cell,
                           rng=rng, hfs=hfs_resp)
    suffix = (args.tag if args.tag is not None
              else ("_hfs" if args.y_source == "hfs" else ""))
    plan = bk.tensor_flip_plan(args.pzz, rel_lumi_offset=args.rel_lumi_offset)
    lumi_assumed = [0.5, 0.5]
    cat_plus = plan.categories[0]
    pzz_list = [0.6 * args.pzz / 0.6, -2.0 * args.pzz]

    def _eff(e2):
        return lambda ph: (1.0 + e2 * np.cos(2 * ph)
                           + args.eff_cos1 * np.cos(ph))

    if args.eff_cos2_split:
        eff2 = [args.eff_cos2 + 0.5 * args.eff_cos2_split,
                args.eff_cos2 - 0.5 * args.eff_cos2_split]
        phi_eff = [_eff(e) for e in eff2]
        bias = reco.fill_acceptance_bias(eff2, pzz_list, lumi_assumed)
        print("fill-dependent phi' efficiency: cos2 harmonic %.4f / %.4f "
              "-> predicted fake amplitude %+.3e (the spin-state ratio "
              "cancels only a COMMON acceptance)" % (eff2[0], eff2[1], bias))
    else:
        phi_eff = _eff(args.eff_cos2)

    # --- money plot 5R -----------------------------------------------------
    fig = plt.figure(figsize=(12.5, 6.8))
    gs = GridSpec(2, 3, figure=fig, width_ratios=(1, 1, 1.35), hspace=0.42,
                  wspace=0.30)
    summary = []
    for k, (xs, qs, i, j) in enumerate(spots):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        xlo, xhi, q2lo, q2hi = superbin_edges(proj, i, j)
        mask = resp.mask_reco(xlo, xhi, q2lo, q2hi)
        summ = resp.bin_summary(xlo, xhi, q2lo, q2hi, cat_plus)
        fit = measure_bin(resp, plan, mask, lumi1_pb, rng, phi_eff, lumi_assumed)
        err10 = reco.err_harmonic_ratio(fit["n"] * lumi10_pb / lumi1_pb,
                                        pzz_list)
        pts, perr = modulation_points(fit)
        centers = 0.5 * (fit["edges"][:-1] + fit["edges"][1:])
        ax.errorbar(centers, 1e3 * pts, yerr=1e3 * perr, fmt="o",
                    color="black", ms=3.5, capsize=2, lw=1, zorder=3)
        phi = np.linspace(0, 2 * np.pi, 200)
        ax.plot(phi, 1e3 * summ["a_reco_bin"] * np.cos(2 * phi), "-",
                color=C_TRUTH, lw=1.6)
        ax.plot(phi, 1e3 * fit["amp"] * np.cos(2 * phi), "--", color=C_FIT,
                lw=1.4)
        ax.set_xlim(0, 2 * np.pi)
        ax.set_xticks([0, np.pi, 2 * np.pi])
        ax.set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
        ax.set_xlabel(r"$\phi' = \phi_e - \phi_S$ (reconstructed)", fontsize=9,
                      labelpad=1)
        if k % 2 == 0:
            ax.set_ylabel(r"$T(\phi')-\hat\kappa$  $[\times10^{-3}]$",
                          fontsize=9)
        ax.tick_params(labelsize=8)
        ax.axhline(0.0, color="0.85", lw=0.6, zorder=0)
        ax.annotate(
            (r"reco bin $x\approx%.3g$, $Q^2\approx%.3g$: purity %.2f, "
             r"$D=%.3f$" "\n"
             r"$\hat A=(%.2f\pm%.2f)\times10^{-3}$ [1 yr]; $\pm%.2f$ [10 yr]")
            % (xs, qs, summ["purity"],
               summ["a_reco_bin"] / summ["a_true_bin"],
               1e3 * fit["amp"], 1e3 * fit["err"], 1e3 * err10),
            xy=(0.5, 1.02), xycoords="axes fraction", ha="center",
            fontsize=7.2)
        if k == 0:
            ax.annotate("reco-bin truth (%s)" % args.delta_model,
                        xy=(0.03, 0.93), xycoords="axes fraction",
                        color=C_TRUTH, fontsize=7)
            ax.annotate("ratio fit", xy=(0.03, 0.84),
                        xycoords="axes fraction", color=C_FIT, fontsize=7)
        summary.append(
            "spot %d: x=%.3g Q2=%.3g N_1yr=%.3g purity=%.2f eff=%.2f "
            "D=%.3f  A_true(bin)=%+.4e A_true(reco)=%+.4e  "
            "A_hat=%+.4e +- %.2e (1yr) +- %.2e (10yr)  sin-term=%+.1e"
            % (k + 1, xs, qs, fit["n"], summ["purity"], summ["efficiency"],
               summ["a_reco_bin"] / summ["a_true_bin"], summ["a_true_bin"],
               summ["a_reco_bin"], fit["amp"], fit["err"], err10,
               fit.get("amp_sin", 0.0)))

    # right: amplitude vs x in reco bins along the spot-1 Q2 slice
    ax = fig.add_subplot(gs[:, 2])
    q2_spot = spots[0][1]
    q2lo, q2hi = q2_spot / 1.6, q2_spot * 1.6
    xe = proj.x_edges
    pts = []
    for i0 in range(0, xe.size - 2, 2):
        xc = np.sqrt(xe[i0] * xe[i0 + 2])
        if not kinematic_mask(xc, q2_spot, s):
            continue
        mask = resp.mask_reco(xe[i0], xe[i0 + 2], q2lo, q2hi)
        if mask.sum() < 50:
            continue
        summ = resp.bin_summary(xe[i0], xe[i0 + 2], q2lo, q2hi, cat_plus)
        m1 = measure_bin(resp, plan, mask, lumi1_pb, rng, phi_eff, lumi_assumed)
        if m1["n"] < 1e3 or m1["err"] > 8e-3 or abs(summ["a_reco_bin"]) < 1e-5:
            continue
        m10 = measure_bin(resp, plan, mask, lumi10_pb, rng, phi_eff,
                          lumi_assumed)
        pts.append((xc, m1, m10, summ))
    xoff = 1.045
    ax.errorbar([p[0] for p in pts], [1e3 * p[1]["amp"] for p in pts],
                yerr=[1e3 * p[1]["err"] for p in pts], fmt="s", mfc="none",
                color=C_FIT, ms=4, capsize=2, lw=1, zorder=3,
                label=r"1-year EIC (%g fb$^{-1}$/u), reco bins" % args.lumi_1yr)
    ax.errorbar([p[0] * xoff for p in pts], [1e3 * p[2]["amp"] for p in pts],
                yerr=[1e3 * p[2]["err"] for p in pts], fmt="o", color="black",
                ms=4, capsize=2, lw=1, zorder=4,
                label=r"10-year EIC (%g fb$^{-1}$/u)" % args.lumi_10yr)
    ax.plot([p[0] for p in pts], [1e3 * p[3]["a_reco_bin"] for p in pts],
            "_", color=C_TRUTH, ms=9, mew=1.4, zorder=5,
            label="reco-bin truth (migration incl.)")
    xg = np.logspace(np.log10(2e-4), np.log10(0.8), 250)
    q2g = np.full_like(xg, q2_spot)
    ok = kinematic_mask(xg, q2g, s)
    f2 = kern.nf2.f2a(xg, q2g) / kern.ion.A
    f1 = kern.nf2.f1a(xg, q2g) / kern.ion.A
    y = y_from_xq2(xg, q2g, s)
    curves = [(model, C_TRUTH, "-", "%s (injected)" % args.delta_model)]
    if args.delta_model != "moment_B":
        curves.append((dm.make("moment_B", variant=args.variant,
                               dilution=args.dilution), C_ALT, "-",
                       "moment_B (conservative)"))
    if args.delta_model != "toy":
        curves.append((dm.make("toy", scale=1e-3), "0.45", "--",
                       r"toy, $\Delta/F_1=10^{-3}$"))
    for dfunc, color, ls, lab in curves:
        amp = a_cos2phi(dfunc(xg, q2g, f1), f1, f2, xg, y)
        ax.plot(xg[ok], 1e3 * amp[ok], ls, color=color, lw=1.5, label=lab)
    ax.set_xscale("log")
    ax.set_xlabel(r"$x$ (reconstructed-bin centre)")
    ax.set_ylabel(r"$\hat A^{\cos 2\phi}$  $[\times 10^{-3}]$")
    ax.axhline(0, color="0.85", lw=0.6, zorder=0)
    ax.tick_params(labelsize=8)
    ax.set_title(r"amplitude vs $x$, $Q^2\approx%.3g$ GeV$^2$, "
                 r"reco bins" % q2_spot, fontsize=9)
    ax.legend(fontsize=7, loc="upper left")

    if args.y_method != "mixed":
        ymeth = "electron method only"
    elif args.y_source == "hfs":
        ymeth = ("mixed: $Q^2_e$ + hadronic $y$ from the %s final state, "
                 "%s method, noise %.0f MeV"
                 % ("TOY" if not args.hfs_sample else "PYTHIA8",
                    args.hfs_method.upper(), 1e3 * args.hfs_noise))
    else:
        ymeth = "mixed: $Q^2_e$ + hadronic $y$ ($%.0f\\%%$)" % (100 * args.y_had_res)
    fig.suptitle(
        r"Reconstructed-level gluonometry (5R), %s, $P_{zz}=%.2f$, "
        r"$\Delta$: %s" "\n"
        r"%s; EMCal $E'$, track angles, $\varepsilon_{\rm eID}(\eta)$; "
        r"spin-state ratio of $P_{zz}=+%.1f/-%.1f$ fills" "\n"
        r"with $\varepsilon(\phi')=1+%.2f\cos2\phi'+%.2f\cos\phi'$ and a "
        r"$%.0e$ relative-luminosity offset switched on (both cancel); "
        r"statistical errors only"
        % (config.label(), plan.pzz_true, model.info(), ymeth, args.pzz,
           2 * args.pzz, args.eff_cos2, args.eff_cos1, args.rel_lumi_offset),
        fontsize=9)
    fig.subplots_adjust(top=0.84, bottom=0.09, left=0.06, right=0.985)
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / ("money_cos2phi_reco_6Li%s.png" % suffix)
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print("delta model:", model.info(),
          "" if q2_ref is None else "(<Q2> = %.3g GeV^2)" % q2_ref)
    print("response: %d MC events, generator sigma = %.4g pb, "
          "selected+eID = %.4g pb" % (resp.w.size, resp.w.sum(),
                                      (resp.w * resp.eff).sum()))
    if resp.sigma_capture is not None:
        fc = resp.sigma_capture[resp.eff > 0]
        print("HFS captured Sigma fraction (selected events): median %.3f, "
              "16/84%% %.3f/%.3f" % (np.median(fc), *np.percentile(fc, [16, 84])))
    for line in summary:
        print(line)

    # --- money plot 7R -----------------------------------------------------
    q2_slices = sorted({round(q2, 3) for _x, q2, _i, _j in spots})
    alt = (dm.make("moment_B", variant=args.variant, dilution=args.dilution)
           if args.delta_model != "moment_B" else dm.make("toy", scale=args.scale))
    alt_label = ("moment_B (conservative)" if args.delta_model != "moment_B"
                 else "toy")
    fig, axes = plt.subplots(1, len(q2_slices), figsize=(12.6, 4.5),
                             sharey=True)
    summary7 = []
    for ax, q2s in zip(np.atleast_1d(axes), q2_slices):
        pts = []
        for i0 in range(0, xe.size - 2, 2):
            xc = np.sqrt(xe[i0] * xe[i0 + 2])
            if not kinematic_mask(xc, q2s, s):
                continue
            mask = resp.mask_reco(xe[i0], xe[i0 + 2], q2s / 1.6, q2s * 1.6)
            if mask.sum() < 50:
                continue
            summ = resp.bin_summary(xe[i0], xe[i0 + 2], q2s / 1.6, q2s * 1.6,
                                    cat_plus)
            m1 = measure_bin(resp, plan, mask, lumi1_pb, rng, phi_eff,
                             lumi_assumed)
            if (m1["n"] < 1e3 or m1["err"] > 8e-3
                    or abs(summ["a_reco_bin"]) < 1e-5):
                continue
            m10 = measure_bin(resp, plan, mask, lumi10_pb, rng, phi_eff,
                              lumi_assumed)
            f1c = kern.nf2.f1a(xc, q2s) / kern.ion.A
            delta_c = float(model(xc, q2s, f1c))
            d1 = rp.delta_from_amplitude(m1, summ, delta_c)
            d10 = rp.delta_from_amplitude(m10, summ, delta_c)
            pts.append({"x": xc, "delta_c": delta_c, "d1": d1["delta"],
                        "e1": d1["err"], "d10": d10["delta"], "e10": d10["err"],
                        "purity": summ["purity"]})
        xoff = 1.045
        ax.errorbar([p["x"] for p in pts], [1e3 * p["x"] * p["d1"] for p in pts],
                    yerr=[1e3 * p["x"] * p["e1"] for p in pts], fmt="s",
                    mfc="none", color=C_FIT, ms=4, capsize=2, lw=1, zorder=3,
                    label=r"1-yr EIC (%g fb$^{-1}$/u), reco bins" % args.lumi_1yr)
        ax.errorbar([p["x"] * xoff for p in pts],
                    [1e3 * p["x"] * p["d10"] for p in pts],
                    yerr=[1e3 * p["x"] * p["e10"] for p in pts], fmt="o",
                    color="black", ms=4, capsize=2, lw=1, zorder=4,
                    label=r"10-yr EIC (%g fb$^{-1}$/u)" % args.lumi_10yr)
        xg = np.logspace(np.log10(2e-4), np.log10(0.8), 250)
        q2g = np.full_like(xg, q2s)
        ok = kinematic_mask(xg, q2g, s)
        f1g = kern.nf2.f1a(xg, q2g) / kern.ion.A
        ax.plot(xg[ok], 1e3 * (xg * np.asarray(model(xg, q2g, f1g)))[ok], "-",
                color=C_TRUTH, lw=1.6, label="%s (injected)" % args.delta_model)
        ax.plot(xg[ok], 1e3 * (xg * np.asarray(alt(xg, q2g, f1g)))[ok], "-",
                color=C_ALT, lw=1.4, label=alt_label)
        ax.set_xscale("log")
        ax.set_xlabel(r"$x$ (reconstructed-bin centre)")
        ax.axhline(0, color="0.85", lw=0.6, zorder=0)
        ax.set_title(r"$Q^2 \approx %.3g$ GeV$^2$" % q2s, fontsize=9)
        ax.tick_params(labelsize=8)
        if pts:
            best = min(pts, key=lambda p: p["e10"])
            summary7.append(
                "Q2=%.3g: %d points; best x=%.3g Delta=%+.4f +- %.4f (1yr) "
                "+- %.4f (10yr), purity %.2f"
                % (q2s, len(pts), best["x"], best["delta_c"], best["e1"],
                   best["e10"], best["purity"]))
    np.atleast_1d(axes)[0].set_ylabel(
        r"$x\,\Delta(x, Q^2)$ per nucleon  $[\times 10^{-3}]$", fontsize=9)
    np.atleast_1d(axes)[0].legend(fontsize=7, loc="lower left")
    fig.suptitle(
        r"Reconstructed-level $\Delta(x,Q^2)$ extraction (7R), %s, "
        r"$P_{zz}=%.2f$: $\hat\Delta=\hat A\,K$, $K$ = model bin-centering "
        r"with migration" "\n%s; spin-state ratio estimator; dilution 1/3 incl.; "
        "stat. only; area = S–S moment $-0.012\\,\\alpha_s/3$"
        % (config.label(), plan.pzz_true, ymeth), fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    out = outdir / ("money_delta_extracted_reco_6Li%s.png" % suffix)
    fig.savefig(out, dpi=140)
    print("wrote", out)
    for line in summary7:
        print(line)

    # --- systematics scan --------------------------------------------------
    if args.syst_scan:
        print("\nDetector systematics: Delta shift per sweet spot "
              "(A_hat is fitted on data generated with the nuisance ON and "
              "converted with the NOMINAL K, so the shift is "
              "a_reco(alt)/a_reco(nominal) - 1).  Common random numbers.")
        edges = [superbin_edges(proj, i, j) for _, _, i, j in spots]
        nominal = [resp.bin_summary(*e, cat_plus)["a_reco_bin"] for e in edges]
        variations = (
            ("electron scale +1%", dict(e_scale=1.01)),
            ("electron scale -1%", dict(e_scale=0.99)),
            ("eps_eID eta tilt +0.05", dict(eid_tilt=0.05)),
            ("EMCal: YR eta table", dict(emcal_eta_table=True)),
            ("hadronic y res 0.30 (correct with %.2f)" % args.y_had_res,
             dict(y_had_res=0.30)),
        )
        # MC noise floor from three independent response seeds
        floor = []
        for seed in (args.seed + 101, args.seed + 202, args.seed + 303):
            alt = rp.RecoResponse(sampler, rmodel,
                                  n_mc_per_cell=args.n_mc_per_cell,
                                  rng=np.random.default_rng(seed),
                                  hfs=hfs_resp)
            floor.append([alt.bin_summary(*e, cat_plus)["a_reco_bin"]
                          for e in edges])
        floor = np.std(np.array(floor) / np.array(nominal), axis=0)
        print("  %-42s %s" % ("MC noise floor (3 seeds)",
                              "  ".join("%+7.3f%%" % (100 * f)
                                        for f in floor)))
        for label, kw in variations:
            alt_model = replace(rmodel, **kw)
            alt = rp.RecoResponse(sampler, alt_model,
                                  n_mc_per_cell=args.n_mc_per_cell,
                                  rng=np.random.default_rng(args.seed),
                                  hfs=hfs_resp)
            shifts = [alt.bin_summary(*e, cat_plus)["a_reco_bin"] / a0 - 1.0
                      for e, a0 in zip(edges, nominal)]
            print("  %-42s %s" % (label,
                                  "  ".join("%+7.3f%%" % (100 * v)
                                            for v in shifts)))


if __name__ == "__main__":
    main()
