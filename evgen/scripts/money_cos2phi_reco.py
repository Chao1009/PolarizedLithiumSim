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
n_mc pseudo-events per sampler cell.

`--lumi-1yr` / `--lumi-10yr` are the PROGRAMME luminosities (one and ten
EIC years) and `--lumi-fraction` is the share of them this observable is
given in the run plan (plans/07 WP2).  It is a THIRD factor, distinct
from the 0.5 / 0.5 spin-state share inside `bookkeeping.tensor_flip_plan`
(which divides this measurement's own luminosity between its fills) and
from the optics luminosity fraction of the coherent channel (which prices
a de-squeezed beta*_x).  Every published number is at --lumi-fraction 1,
and a non-default share writes its own PNG.  Outputs:

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
from polligen import hfs, radiative as rad, reco, recopseudo as rp  # noqa: E402
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
    ap.add_argument("--lumi-fraction", type=float, default=1.0,
                    dest="lumi_fraction",
                    help="this observable's share of the PROGRAMME "
                         "luminosity (plans/07 WP2; default 1.0, which "
                         "every published number assumes).  Not the "
                         "0.5/0.5 spin-state share of the flip plan, "
                         "which divides this measurement's own luminosity")
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
    ap.add_argument("--hfs-calibrate", action="store_true",
                    help="divide the transferred hadronic sums by the library's per-cell "
                         "mean captured fraction (the analysis's own hadronic-scale "
                         "calibration); off reproduces the 2026-08-26 published numbers")
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
    ap.add_argument("--unfold", default="model", choices=("model", "folded"),
                    help="7R bin-centering.  'model' (default, the "
                         "published one): K = Delta_prior(x_c, Q2_c) / "
                         "A_reco-bin(prior), evaluated with the model that "
                         "generated the pseudo-data, so the points sit on "
                         "the injected curve by construction (code review "
                         "F5).  'folded': a 3-parameter Delta(x) shape "
                         "FITTED THROUGH THE RESPONSE per Q2 slice "
                         "(recopseudo.fold_shape_fit), so the shape comes "
                         "from the data and the shape-fit, response-MC "
                         "and prior-spread errors of K all enter the "
                         "error bar (plans/08 A6)")
    ap.add_argument("--unfold-prior", default=None, choices=dm.available(),
                    help="model the bin-centering STARTS from (default: "
                         "the injected one).  With --unfold model this is "
                         "the F5 model dependence itself; with --unfold "
                         "folded it is only the starting shape of the fit")
    ap.add_argument("--unfold-scan", action="store_true",
                    help="closure scan: correct the SAME pseudo-data with "
                         "each of the repository's Delta shapes, bin-by-bin "
                         "K against folded fit, and print the residual bias "
                         "against the injected truth (plans/08 A6)")
    ap.add_argument("--syst-scan", action="store_true",
                    help="after the money plots, rebuild the response with "
                         "each nuisance varied and print the Delta shift "
                         "per sweet spot (common random numbers, so the "
                         "shifts are not seed noise)")
    ap.add_argument("--isr", action="store_true",
                    help="after the money plots, print the RADIATIVE-"
                         "CORRECTION migration bound at the four sweet "
                         "spots (plans/07 WP4): the response is rebuilt "
                         "with collinear leading-log initial-state "
                         "radiation switched on, with common random "
                         "numbers (polligen.radiative draws z from its own "
                         "stream), and the shift of Delta_hat that an "
                         "ISR-free bin-centering would leave is compared "
                         "with the 5%% gate")
    ap.add_argument("--isr-gen-q2min", type=float, default=None,
                    help="rebuild the ISR pair from a generator window "
                         "reaching down to this Q2 [GeV^2] instead of the "
                         "0.7 the money plots use.  Events below the "
                         "window cannot feed the analysis bins under "
                         "radiation, so the published window UNDERSTATES "
                         "the bound; 0.05 is where it converges")
    ap.add_argument("--isr-empz", type=float, default=None,
                    help="also apply a HERA-style E - p_z > this fraction "
                         "of 2E_e to both members of the pair (0.85 is the "
                         "cut H1/ZEUS used against radiative events); the "
                         "chain does not apply it, so this is the "
                         "MITIGATION column, not the bound.  The retention "
                         "is printed three ways -- globally, in bands of "
                         "nominal y, and per analysis bin -- because the "
                         "global number is dominated by the high-y bulk "
                         "and is NOT what a sweet spot pays")
    ap.add_argument("--isr-seeds", default=None,
                    help="comma-separated RESPONSE seeds over which the "
                         "bound is averaged (mean +- sem).  One draw of "
                         "the response scatters by 4-14%% of the bound "
                         "(seed-to-seed sd 0.05-0.09 percentage points), "
                         "so the PUBLISHED numbers are the average over "
                         "this list; without it the single --seed run is "
                         "printed")
    ap.add_argument("--isr-seed", type=int, default=20260828)
    ap.add_argument("--seed", type=int, default=20260824)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    if not args.lumi_fraction > 0:
        ap.error("--lumi-fraction must be positive")
    config = beams.default_configs("6Li")[args.config]
    # programme luminosity x this observable's share of the run plan; the
    # spin-state share of `plan` divides what is left, and is a different
    # object (see the module docstring)
    lumi1_pb = args.lumi_1yr * args.lumi_fraction * 1e3
    lumi10_pb = args.lumi_10yr * args.lumi_fraction * 1e3
    rng = np.random.default_rng(args.seed)
    analysis = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            run_share=args.lumi_fraction,
                            pol_ion_tensor=args.pzz)
    print("run plan: programme %g / %g fb^-1/u (1 yr / 10 yr) x share %g "
          "-> %g / %g fb^-1/u delivered to this observable; spin-state "
          "share 0.5 / 0.5 within it"
          % (args.lumi_1yr, args.lumi_10yr, args.lumi_fraction,
             args.lumi_1yr * args.lumi_fraction,
             args.lumi_10yr * args.lumi_fraction))
    model, q2_ref = build_delta_model(args, config, analysis)
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1, delta_func=model)
    # the shape the bin-centering starts from: the injected model unless
    # asked otherwise (a different prior is the F5 model dependence)
    prior_name = args.unfold_prior or args.delta_model

    def make_prior(name):
        """One of the repository's Delta shapes, built with this run's
        beam configuration and analysis scenario."""
        if name == args.delta_model:
            return model
        return build_delta_model(
            argparse.Namespace(**dict(vars(args), delta_model=name)),
            config, analysis)[0]

    prior = model if args.unfold_prior in (None, args.delta_model) else \
        make_prior(args.unfold_prior)

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
                                   scale=args.hfs_scale, calibrate=args.hfs_calibrate)
        print("HFS library: %s; %s method; %s | coverage %s"
              % (src, args.hfs_method.upper(), hresp.describe(), lib.coverage()))
    resp = rp.RecoResponse(sampler, rmodel, n_mc_per_cell=args.n_mc_per_cell,
                           rng=rng, hfs=hfs_resp)
    suffix = (args.tag if args.tag is not None
              else ("_hfs" if args.y_source == "hfs" else ""))
    # published PNGs are the --lumi-fraction 1 ones: a non-default share
    # appends its key rather than overwriting them (the same guard as
    # money_tagged_azz.output_stem)
    share_key = fom.run_share_tag(args.lumi_fraction)
    if share_key and args.tag is None:
        suffix = "%s_%s" % (suffix, share_key)
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
        # the bins shown are chosen on MEASURABLE quantities only: the
        # expected count and the statistical error (an earlier version also
        # required the reco-bin truth amplitude to be non-zero, which is
        # inert for every registry model but is not a criterion an
        # experiment has; code review 2026-08-28)
        if m1["n"] < 1e3 or m1["err"] > 8e-3:
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
    folded = args.unfold == "folded"
    for ax, q2s in zip(np.atleast_1d(axes), q2_slices):
        cand = []
        for i0 in range(0, xe.size - 2, 2):
            xc = np.sqrt(xe[i0] * xe[i0 + 2])
            # the folded fit also measures the FEED-IN bins on either side
            # of the plotted range (bins whose CENTRE is outside the
            # kinematic mask, or that the plot's quality cuts drop): the
            # edge bins of the plot are fed by true x from there, and a
            # shape constrained only inside would extrapolate into exactly
            # the region that feeds them (plans/08 A6)
            in_range = bool(kinematic_mask(xc, q2s, s))
            if not (in_range or folded):
                continue
            mask = resp.mask_reco(xe[i0], xe[i0 + 2], q2s / 1.6, q2s * 1.6)
            if mask.sum() < 50:
                continue
            summ = resp.bin_summary(xe[i0], xe[i0 + 2], q2s / 1.6, q2s * 1.6,
                                    cat_plus)
            m1 = measure_bin(resp, plan, mask, lumi1_pb, rng, phi_eff,
                             lumi_assumed)
            good = in_range and not (m1["n"] < 1e3 or m1["err"] > 8e-3)
            if not (good or folded):
                continue
            m10 = measure_bin(resp, plan, mask, lumi10_pb, rng, phi_eff,
                              lumi_assumed)
            cand.append({"x": xc, "mask": mask, "summ": summ, "m1": m1,
                         "m10": m10, "good": good})
        shape_fit = {}
        if folded:
            # the fit is repeated from every OTHER shape the repository
            # offers, on the same data: the spread of K over that family
            # is the residual prior dependence the tilt cannot absorb,
            # and it is the largest of the three K errors (plans/08 A6)
            alt_names = [nm for nm in dm.available() if nm != prior_name]
            alts = tuple(make_prior(nm) for nm in alt_names)
            for key in ("m1", "m10"):
                shape_fit[key] = rp.fold_shape_fit(
                    resp, cat_plus,
                    [{"mask": c["mask"], "x": c["x"], "q2": q2s,
                      "amp": c[key]["amp"], "err": c[key]["err"]}
                     for c in cand], prior, alt_bases=alts)
        pts = []
        for b, c in enumerate(cand):
            if not c["good"]:
                continue
            xc, summ, mask = c["x"], c["summ"], c["mask"]
            f1c = kern.nf2.f1a(xc, q2s) / kern.ion.A
            delta_c = float(prior(xc, q2s, f1c))
            kw1 = kw10 = {}
            if folded:
                kw1 = {"k_conv": shape_fit["m1"]["k"][b],
                       "k_rel_err": shape_fit["m1"]["k_rel_all"][b]}
                kw10 = {"k_conv": shape_fit["m10"]["k"][b],
                        "k_rel_err": shape_fit["m10"]["k_rel_all"][b]}
            elif prior is not model:
                kw1 = kw10 = {"k_conv": delta_c / resp.fold(prior, mask,
                                                            cat_plus)}
            d1 = rp.delta_from_amplitude(c["m1"], summ, delta_c, **kw1)
            d10 = rp.delta_from_amplitude(c["m10"], summ, delta_c, **kw10)
            pts.append({"x": xc, "delta_c": delta_c, "d1": d1["delta"],
                        "e1": d1["err"], "d10": d10["delta"],
                        "e10": d10["err"], "e1_stat": d1["err_stat"],
                        "e1_k": d1["err_k"], "e10_stat": d10["err_stat"],
                        "e10_k": d10["err_k"], "purity": summ["purity"],
                        "b": b})
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
        if folded:
            sh = shape_fit["m10"]["shape"](xg, q2g, f1g)
            ax.plot(xg[ok], 1e3 * (xg * np.asarray(sh))[ok], "--",
                    color="0.25", lw=1.2,
                    label=r"folded fit (10 yr), $\chi^2/{\rm ndf}=%.1f/%d$"
                          % (shape_fit["m10"]["chi2"],
                             shape_fit["m10"]["ndof"]))
        ax.set_xscale("log")
        ax.set_xlabel(r"$x$ (reconstructed-bin centre)")
        ax.axhline(0, color="0.85", lw=0.6, zorder=0)
        ax.set_title(r"$Q^2 \approx %.3g$ GeV$^2$" % q2s, fontsize=9)
        ax.tick_params(labelsize=8)
        if pts:
            best = min(pts, key=lambda p: p["e10"])
            extra = ""
            if folded:
                fit10 = shape_fit["m10"]
                b1, b10 = abs(best["d1"]), abs(best["d10"])
                bb = best["b"]
                extra = ("; Delta_hat(10 yr) = %+.4f; folded fit %d bins "
                         "(%d plotted) tilt (%+.3f, %+.3f) chi2/ndf = "
                         "%.3g/%d%s; error at that bin [%%]: 1 yr stat "
                         "%.2f (+) K %.2f = %.2f, 10 yr stat %.2f (+) K "
                         "%.2f = %.2f; K (10 yr) = shape-fit %.2f (+) "
                         "response-MC %.2f (+) prior-spread %.2f"
                         % (best["d10"], len(cand), len(pts),
                            fit10["params"][1],
                            fit10["params"][2], fit10["chi2"], fit10["ndof"],
                            " AT BOUND" if fit10["at_bound"] else "",
                            100 * best["e1_stat"] / b1,
                            100 * best["e1_k"] / b1, 100 * best["e1"] / b1,
                            100 * best["e10_stat"] / b10,
                            100 * best["e10_k"] / b10,
                            100 * best["e10"] / b10,
                            100 * fit10["k_fit_rel"][bb],
                            100 * fit10["k_mc_rel"][bb],
                            100 * fit10["k_prior_rel"][bb]))
                # no goodness-of-fit penalty enters the prior spread, so
                # say how well each alternative actually fits: where its
                # chi2 is far worse than the prior's, the spread is an
                # over-estimate of the residual (recopseudo.fold_shape_fit)
                extra += ("; prior fits (10 yr) chi2/ndf: %s = %.3g/%d, %s"
                          % (prior_name, fit10["chi2"], fit10["ndof"],
                             ", ".join("%s = %.3g/%d" % (nm, c2,
                                                         fit10["ndof"])
                                       for nm, c2 in zip(alt_names,
                                                         fit10["alt_chi2"]))))
            summary7.append(
                "Q2=%.3g: %d points; best x=%.3g Delta=%+.4f +- %.4f (1yr) "
                "+- %.4f (10yr), purity %.2f%s"
                % (q2s, len(pts), best["x"], best["delta_c"], best["e1"],
                   best["e10"], best["purity"], extra))
    np.atleast_1d(axes)[0].set_ylabel(
        r"$x\,\Delta(x, Q^2)$ per nucleon  $[\times 10^{-3}]$", fontsize=9)
    np.atleast_1d(axes)[0].legend(fontsize=7, loc="lower left")
    if folded:
        method = (r"$\hat\Delta=\hat A\,K$, $K$ from a 3-parameter "
                  r"$\Delta(x)$ shape fitted THROUGH the response "
                  r"(prior %s)" % prior_name)
        errs = "stat. (+) shape-fit (+) response-MC (+) prior spread"
    else:
        method = (r"$\hat\Delta=\hat A\,K$, $K$ = model bin-centering "
                  r"with migration")
        errs = "stat. only"
    fig.suptitle(
        r"Reconstructed-level $\Delta(x,Q^2)$ extraction (7R), %s, "
        r"$P_{zz}=%.2f$: %s"
        "\n%s; spin-state ratio estimator; dilution 1/3 incl.; "
        "%s; area = S–S moment $-0.012\\,\\alpha_s/3$"
        % (config.label(), plan.pzz_true, method, ymeth, errs), fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    out = outdir / ("money_delta_extracted_reco_6Li%s.png" % suffix)
    fig.savefig(out, dpi=140)
    print("wrote", out)
    for line in summary7:
        print(line)

    # --- bin-centering closure scan (plans/08 A6) --------------------------
    if args.unfold_scan:
        print("\nBin-centering closure: the SAME noise-free pseudo-data "
              "(generated with %s, exact expected counts through the ratio "
              "fit) corrected with each of the repository's Delta shapes, "
              "bin-by-bin K against the folded shape fit.  Rows are the "
              "residual bias of Delta_hat against the injected truth at the "
              "bin centre; the prior=%s rows are the closure floor of the "
              "chain itself.  Common random numbers throughout (one "
              "response, one data set)." % (args.delta_model,
                                            args.delta_model))
        priors = {nm: make_prior(nm) for nm in dm.available()}
        def row(lab, vals, fmt="%+7.2f"):
            return ("  %-26s" % lab) + "".join((fmt + " ") % v for v in vals)

        spot_k = {}
        for q2s in q2_slices:
            sbins, plotted = [], []
            for i0 in range(0, xe.size - 2, 2):
                xc = np.sqrt(xe[i0] * xe[i0 + 2])
                mask = resp.mask_reco(xe[i0], xe[i0 + 2], q2s / 1.6, q2s * 1.6)
                if mask.sum() < 50:
                    continue
                summ = resp.bin_summary(xe[i0], xe[i0 + 2], q2s / 1.6,
                                        q2s * 1.6, cat_plus)
                exact = rp.measure_inclusive(resp, plan, lumi1_pb, mask,
                                             poisson=False, phi_eff=phi_eff,
                                             lumi_assumed=lumi_assumed)
                sbins.append({"x": xc, "q2": q2s, "mask": mask,
                              "amp": exact["amp"], "err": exact["err"],
                              "summ": summ})
                plotted.append(bool(kinematic_mask(xc, q2s, s))
                               and exact["n"] >= 1e3 and exact["err"] <= 8e-3)
            xs = np.array([b["x"] for b in sbins])
            amps = np.array([b["amp"] for b in sbins])
            e1 = np.array([b["err"] for b in sbins])
            keep = np.array(plotted)
            f1s = resp.f1_center(xs, np.full_like(xs, q2s))
            truth = np.asarray(model(xs, np.full_like(xs, q2s), f1s), float)
            print("\n Q2 = %.3g GeV^2: %d bins fitted, %d plotted "
                  "(x = %.2e to %.2e)"
                  % (q2s, len(sbins), int(keep.sum()), xs.min(), xs.max()))
            print(row("x centre", xs, "%7.1e"))
            print(row("plotted", keep.astype(int), "%7d"))
            print(row("stat err 1 yr [%]", 100 * e1 / np.abs(amps)))
            print(row("stat err 10 yr [%]", 100 * e1 / np.abs(amps)
                      * np.sqrt(lumi1_pb / lumi10_pb)))
            k_folded = {}
            for nm, pri in sorted(priors.items()):
                kbin = np.array([np.asarray(pri(b["x"], q2s, f), float)
                                 / resp.fold(pri, b["mask"], cat_plus)
                                 for b, f in zip(sbins, f1s)])
                bias = amps * kbin / truth - 1.0
                print(row("%-9s K bin-by-bin" % nm, 100 * bias)
                      + " | max %.1f%% (plotted)"
                      % (100 * np.max(np.abs(bias[keep]))))
                fit = rp.fold_shape_fit(resp, cat_plus, sbins, pri)
                bias = amps * fit["k"] / truth - 1.0
                print(row("%-9s folded fit" % nm, 100 * bias)
                      + " | max %.1f%% (plotted), chi2/ndf = %.3g/%d"
                      % (100 * np.max(np.abs(bias[keep])), fit["chi2"],
                         fit["ndof"]))
                if nm == args.delta_model:
                    print(row("  K response-MC err [%]",
                              100 * fit["k_mc_rel"]))
                print(row("  K shape-fit err [%]", 100 * fit["k_fit_rel"]))
                k_folded[nm] = fit["k"]
                spot_k[(q2s, nm)] = fit["shape"]
            # what --unfold folded puts in the bar as the third K error,
            # printed here against the folded-fit rows it has to cover:
            # the spread of K over the family, seen from each prior
            for nm in sorted(priors):
                sprd = np.max([np.abs(k_folded[o] / k_folded[nm] - 1.0)
                               for o in priors if o != nm], axis=0)
                print(row("%-9s K prior-spread [%%]" % nm, 100 * sprd)
                      + " | max %.1f%% (plotted)"
                      % (100 * np.max(sprd[keep])))
        # the F5 comparison: the same K spread at the four sweet spots
        print("\n Sweet spots (code review F5, same convention): "
              "K(%s)/K(prior) - 1 [%%], the model dependence of the "
              "bin-centering itself" % args.delta_model)
        edges4 = [superbin_edges(proj, i, j) for _x, _q2, i, j in spots]
        q2c = [round(q2, 3) for _x, q2, _i, _j in spots]
        xc4 = [xsp for xsp, _q2, _i, _j in spots]
        f1_4 = resp.f1_center(np.array(xc4), np.array(q2c))
        for nm in sorted(priors):
            if nm == args.delta_model:
                continue
            for lab, shp in (("K bin-by-bin", None), ("folded fit", "fit")):
                out = []
                for e, xsp, q2sp, f1sp in zip(edges4, xc4, q2c, f1_4):
                    mask = resp.mask_reco(*e)
                    ref = priors[args.delta_model] if shp is None else \
                        spot_k[(q2sp, args.delta_model)]
                    alt_s = priors[nm] if shp is None else spot_k[(q2sp, nm)]
                    k_ref = (float(ref(xsp, q2sp, f1sp))
                             / resp.fold(ref, mask, cat_plus))
                    k_alt = (float(alt_s(xsp, q2sp, f1sp))
                             / resp.fold(alt_s, mask, cat_plus))
                    out.append(100 * (k_ref / k_alt - 1.0))
                print(row("%-9s %s" % (nm, lab), out))

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

    # --- radiative-correction migration bound (plans/07 WP4) ---------------
    if args.isr:
        isr_sampler = sampler
        if args.isr_gen_q2min is not None:
            isr_gen = rp.generator_scenario(
                analysis, q2_min=args.isr_gen_q2min,
                y_min=min(gen.y_min, 0.2 * analysis.y_min))
            isr_sampler = InclusiveSampler(kern, config, isr_gen, nx=60,
                                           nq2=45,
                                           q2_range=(args.isr_gen_q2min, 2e3))

        def build_isr(seed, isr_model):
            """Both members from the SAME response seed: the ISR model has
            its own stream, so the pair sits on identical pseudo-events
            and the difference is the radiation, not seed noise."""
            return rp.RecoResponse(isr_sampler, rmodel,
                                   n_mc_per_cell=args.n_mc_per_cell,
                                   rng=np.random.default_rng(seed),
                                   hfs=hfs_resp, isr=isr_model)

        isr_seeds = ([int(v) for v in args.isr_seeds.split(",")]
                     if args.isr_seeds else None)

        edges4 = [superbin_edges(proj, i, j) for _x, _q2, i, j in spots]
        labels4 = ["x=%.3g Q2=%.3g" % (xs, qs) for xs, qs, _i, _j in spots]
        print("\nRadiative corrections (plans/07 WP4): the collinear-ISR "
              "MIGRATION bound.\nCollinear ISR fakes no cos phi' or "
              "cos 2phi' (the covariant azimuth is invariant under "
              "k -> (1-z)k), and\nanything common to the fills cancels in "
              "the spin-state ratio; what is left is the Q2_e LABEL "
              "migrating\nby 1/(1-z) while x = Q2_e/(s y_Sigma) stays "
              "exact.  dDelta is A_reco(ISR)/A_reco(no ISR) - 1, i.e. what "
              "an\nISR-free bin-centering leaves on Delta_hat if the data "
              "radiate.  Common random numbers; %d MC/cell, response "
              "seed %s.\nNOT bounded here: the TENSOR-sector radiative "
              "correction, which has never been calculated (plans/05 5.5)."
              % (args.n_mc_per_cell,
                 args.isr_seeds if isr_seeds else args.seed))
        empz = (None if args.isr_empz is None
                else (args.isr_empz, 2.0 - args.isr_empz))
        if isr_seeds:
            bound = rad.migration_bound_seeds(
                build_isr, edges4, cat_plus, isr_seeds,
                isr_seed=args.isr_seed, labels=labels4, empz_cut=empz)
        else:
            bound = rad.migration_bound(
                lambda m: build_isr(args.seed, m), edges4, cat_plus,
                labels=labels4, isr=rad.ISRModel(seed=args.isr_seed),
                empz_cut=empz)
        print(rad.format_bound(bound))
        on = bound["on"]
        w_on = on.w * on.eff
        z = on.isr_z
        print("  <z> = %.4f over the selected rate; P(z > 1e-4) = %.3f; "
              "<z | z > 1e-4> = %.4f; t = %.4f at <Q2> = %.3g GeV^2"
              % (np.average(z, weights=w_on),
                 np.average(z > 0.0, weights=w_on),
                 np.average(z[z > 0], weights=w_on[z > 0]),
                 np.average(rad.beta_ll(on.q2), weights=w_on),
                 np.average(on.q2, weights=w_on)))
        # The method comparison is a STRONG function of the nominal y of
        # the bin (the electron-method ratios go as (y + z)/y), so it is
        # printed at a stated y -- both the rate-weighted mean of the whole
        # selected sample, which the high-y cells the money plots do not
        # use dominate, and the four sweet spots themselves, which is where
        # the letter's numbers live.  `on.y` is the HARD y with ISR on; the
        # table wants the NOMINAL one (radiative.method_bias_table).
        z_typ = float(np.average(z[z > 0], weights=w_on[z > 0]))
        y_bar = float(np.average(on.y_nominal, weights=w_on))
        y_spot = float(np.mean([qs / (sampler.s * xs)
                                for xs, qs, _i, _j in spots]))
        d_phi = on.isr_dphi
        print("  covariant azimuth under k -> (1-z)k (physical 6Li mass): "
              "max |dphi'| = %.2g rad,\n    1 - <cos 2 dphi'> = %.2g "
              "rate-weighted over %d events -- the fake cos 2phi' collinear "
              "ISR can make"
              % (np.abs(d_phi).max(),
                 1.0 - np.average(np.cos(2.0 * d_phi), weights=w_on),
                 d_phi.size))
        print("  method comparison at z = <z | z > 1e-4> = %.4f "
              "(ratio observed/hard; the chain uses the last row):" % z_typ)
        for y_typ, what in ((y_bar, "rate-weighted <y> of the whole "
                                    "selected sample"),
                            (y_spot, "mean y of the four sweet spots")):
            print("    at y = %.4f (%s)" % (y_typ, what))
            print("      %-24s %10s %10s %10s" % ("method", "Q2", "y", "x"))
            for lab, rq2, ry, rx in rad.method_bias_table(y_typ, z_typ):
                print("      %-24s %10.4f %10.4f %10.4f" % (lab, rq2, ry, rx))


if __name__ == "__main__":
    main()
