#!/usr/bin/env python3
"""Money plot 6R: the coherent intact-6Li channel at the RECONSTRUCTED
level (plans/07 WP3 -> WP5) -- money plot 6 re-derived with

  * the ANGULAR near-beam cut, pT_cut = 10 sigma_theta A p_u, and the
    Roman-Pot cutout.  The PUBLISHED configuration (2026-08-28) is
    `--optics tagging`: the lithium tagging optics of Report 1 Section
    6.1 (reco.tagging_optics_point) -- the horizontal beta* de-squeezed
    to the optimum of acceptance x luminosity, the vertical plane at the
    Yellow Report high-acceptance divergence, the pots following the
    10 sigma envelope in both planes (0.33 x 3.8 mrad at 5 x 40.8) --
    because with the Yellow Report optics and the measured pot aperture
    no recoil survives at any configuration (plans/10).  The legacy
    defaults (proton-derived 73 microrad, a 2.5 : 1 slot) reproduce the
    pre-2026-08-27 figure and are kept only for that;
  * the Roman-Pot emulation: divergence smearing of the recoil angle,
    |t| = pT^2 reconstructed with x_L = 1, reco t bins.  The PUBLISHED
    binning (2026-08-28) is the seven bins of
    `recopseudo.T_EDGES_PUBLISHED`, 0.017 to 0.25 GeV^2: the aperture
    floor at the tagging optics is |t|_min = 0.0064 / 0.0098 / 0.0094
    GeV^2 at the three configurations, so the window is an ANALYSIS
    choice above it, and it stops at 0.017 because a 0.006-0.017 bin
    empties a third to a half of its beta cells.  a_t is quoted at the
    bin's t_ref, the rate-weighted mean TRUE |t|, and in the four bins
    whose t_ref falls below 0.05 GeV^2 -- 0.017-0.028 through 0.05-0.08,
    at t_ref = 0.023, 0.029, 0.036 and 0.045 -- it is the linear-in-|t|
    deformation model extrapolated below the lowest digitized anchor
    point (coherent.MANTYSAARI_A2_DEUTERON starts at 0.05).  The retired
    four-bin window's own lowest bin was already one of those four, so
    the wider window deepens the extrapolation rather than introducing
    it.  The run-13 window is `--t-edges 0.05,0.08,0.12,0.17,0.25`;
  * the TWO azimuths: alpha = phi_e - phi_S (electron) and beta =
    phi_t - phi_S (recoil); the deformation term modulates cos 2beta and
    the gluon-transversity term cos 2alpha; the unpolarized lepton-plane
    /recoil-plane harmonics u_1 cos(alpha-beta), u_2 cos 2(alpha-beta)
    (the L-T and T-T' interference of diffractive DIS, Nikolaev-
    Pronyaev-Zakharov hep-ph/9812212) are set at the 1-sigma edge of the
    ZEUS LPS measurement (A_LT, A_TT consistent with zero within
    +-0.03-0.05; NPB 816 (2009) 1, Sec. 10.2);
  * the anchor convention (arXiv:2408.13213 Eq. 9: 1 + 2 sum a_n cos n Phi):
    the cos 2beta COEFFICIENT is 2 a_2 -- twice what money plot 6 injects;
  * the spin-state-sorted 2-D ratio fit (reco.harmonic_ratio_fit_2d) of
    P_zz = +0.6 / -1.2 fills, which cancels the cutout's fake harmonics.
    `--fit likelihood` swaps it for the acceptance-profiled Poisson
    likelihood (reco.harmonic_likelihood_fit_2d), which is unbiased at
    any count and is what the sparse |t| bins need;
  * `--u-in-situ`, which measures the unpolarized harmonics (u1, u2)
    from the spin-averaged counts of the same data against the response's
    acceptance shape and propagates their covariance into the harmonic
    errors, instead of assuming ZEUS's values.

Three luminosity factors multiply and are NOT the same object:
`--lumi-1yr` is the PROGRAMME luminosity of one EIC year; `--lumi-fraction`
is this observable's share of it in the run plan (plans/07 WP2), ours to
propose and 1.0 in every published number; and the optics fraction
L/L_HA = 1/7 - 1/13 (`reco.tagging_optics_point`, printed beside it) is
what the de-squeezed beta*_x costs at fixed wall time.  Inside all three,
the flip plan divides the luminosity 0.5 / 0.5 between the two spin
states.

Panels: (a) raw spin-sorted yields vs beta (acceptance-dominated);
(b) the acceptance-free modulation T projected on beta and on alpha with
the fits; (c) the expected 2-D modulation map; (d) the extracted a_t(t)
and a_e per reco t bin (1 yr / 10 yr) against the injected curves.

Usage:  python3 scripts/money_cos2phi_coherent_reco.py
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
from polligen import coherent as coh  # noqa: E402
from polligen import reco, recopseudo as rp  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.farforward import HIGH_ACCEPTANCE  # noqa: E402

C_TRUTH, C_FIT, C_ALT, C_GREY = "#0072B2", "#D55E00", "#009E73", "0.45"


def output_stem(args):
    """File stem for one run.  Money plot 6R -- Report 2's __RC5__ --
    is ONE combination: `--config 0 --optics tagging` with the default
    ratio fit, the assumed (u1, u2) and the published |t| edges.  Every
    other combination gets its keys appended, because the reproduction
    manual documents non-default runs (`--fit likelihood`, `--u-in-situ`,
    `--t-edges`) whose figure would otherwise silently overwrite the
    published one -- the bug `money_tagged_azz.output_stem` fixed for the
    tagged scripts on 2026-08-28."""
    base = "money_cos2phi_coherent_reco_6Li"
    share_key = fom.run_share_tag(getattr(args, "lumi_fraction", 1.0))
    if (args.config == 0 and args.optics == "tagging"
            and args.fit == "ratio" and not args.u_in_situ
            and args.t_edges is None and not share_key):
        return base
    keys = ["c%d" % args.config, args.optics]
    if args.fit != "ratio":
        keys.append(args.fit)
    if args.u_in_situ:
        keys.append("uinsitu")
    if args.t_edges is not None:
        keys.append("tedges")
    if share_key:                      # a run-plan share is a non-default run
        keys.append(share_key)
    return "%s_%s" % (base, "_".join(keys))


def t_edges_for(args):
    """The reconstructed |t| edges of one run.  `args.t_edges is None`
    is the PUBLISHED sentinel -- the one `output_stem` keys on -- and it
    resolves to the seven-bin window adopted on 2026-08-28,
    `recopseudo.T_EDGES_PUBLISHED` (0.017-0.25 GeV^2).  The run-13
    window survives as `--t-edges 0.05,0.08,0.12,0.17,0.25`
    (`recopseudo.T_EDGES_LEGACY`), which appends its key to the stem."""
    if args.t_edges is None:
        return list(rp.T_EDGES_PUBLISHED)
    edges = [float(v) for v in args.t_edges.split(",")]
    if len(edges) < 2 or any(hi <= lo for lo, hi in zip(edges[:-1],
                                                        edges[1:])):
        raise SystemExit("--t-edges wants at least two increasing edges, "
                         "not %r" % (args.t_edges,))
    return edges


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2))
    ap.add_argument("--lumi-1yr", type=float, default=10.0)
    ap.add_argument("--lumi-10yr", type=float, default=100.0)
    ap.add_argument("--lumi-fraction", type=float, default=1.0,
                    dest="lumi_fraction",
                    help="this observable's share of the PROGRAMME "
                         "luminosity (plans/07 WP2; default 1.0, which "
                         "every published number assumes).  Distinct from "
                         "the optics fraction L/L_HA of --optics tagging, "
                         "which is printed beside it, and from the 0.5/0.5 "
                         "spin-state share of the flip plan")
    ap.add_argument("--pzz", type=float, default=0.60)
    ap.add_argument("--eps-b0", type=float, default=-0.08)
    ap.add_argument("--amp", type=float, default=0.01,
                    help="flat gluon-transversity cos 2alpha coefficient "
                         "per unit P_zz (band 3e-3..1e-2)")
    ap.add_argument("--a-m", type=float, default=0.0,
                    help="mixed cos(alpha+beta) coefficient per unit P_zz")
    ap.add_argument("--u1", type=float, default=0.05,
                    help="unpolarized cos(alpha-beta) coefficient (L-T "
                         "interference).  ZEUS LPS, NPB 816 (2009) 1 "
                         "[refs/0812.2003]: A_LT = -0.036 +- 0.036 (x_P < "
                         "0.01), +0.051 +- 0.024 (0.01-0.1); default at the "
                         "1-sigma edge")
    ap.add_argument("--u2", type=float, default=0.02,
                    help="unpolarized cos 2(alpha-beta) coefficient (T-T' "
                         "interference).  ZEUS LPS: A_TT = -0.030 +- 0.037, "
                         "-0.010 +- 0.024; default inside the 1-sigma band")
    ap.add_argument("--optics", default="legacy",
                    choices=("legacy", "tagging", "high-acceptance"),
                    help="beam optics.  'tagging': the lithium tagging "
                         "optics of Report 1 Section 6.1 -- horizontal "
                         "beta* at the optimum of acceptance x luminosity, "
                         "vertical at high acceptance, pots following the "
                         "10 sigma envelope in both planes -- the "
                         "PUBLISHED 6R (sets --sigma-theta, --aspect, "
                         "--cut-scale-x/y = 1, --rp-aperture measured and "
                         "--near-beam-mrad from reco.tagging_optics_point); "
                         "'high-acceptance': the Yellow Report divergence "
                         "per configuration (reco.sigma_theta_for), under "
                         "which the tag is dead; 'legacy': the pre-2026-08-27 "
                         "proton-derived 73 microrad")
    ap.add_argument("--sigma-theta", type=float, default=reco.SIGMA_THETA_HA,
                    help="horizontal beam angular divergence [rad] "
                         "(legacy default: proton-derived 73 microrad)")
    ap.add_argument("--aspect", type=float, default=1.0,
                    help="beam-divergence anisotropy sigma_y/sigma_x "
                         "(HERA: 100/45 MeV vertical/horizontal pT spread, "
                         "ZEUS NPB 816:1); set by --optics tagging")
    ap.add_argument("--ensemble", type=int, default=0,
                    help="repeat the one-year pseudo-experiment this many "
                         "times with fresh seeds and print the mean and "
                         "spread of the fitted coefficients against the "
                         "injected values -- the bias test of the template "
                         "basis (2026-08-28)")
    ap.add_argument("--cut-scale-x", type=float, default=2.5,
                    help="cutout half-width in x in units of 10 sigma_x: "
                         "the ePIC pots surround a horizontal SLOT (beam "
                         "momentum spread + dispersion), Jentsch DIS 2023 "
                         "slide 15 [illustrative value]")
    ap.add_argument("--cut-scale-y", type=float, default=1.0,
                    help="cutout half-height in y in units of 10 sigma_y")
    ap.add_argument("--shape", default="rectangle",
                    choices=("rectangle", "ellipse"))
    ap.add_argument("--rp-aperture", default="none",
                    choices=("none", "measured"),
                    help="add the pots' GEOMETRIC aperture, measured in "
                         "the ePIC geometry (reco.RP_APERTURE_MEASURED, "
                         "tools/fullsim), as a floor under the envelope "
                         "cutout.  'none' is the envelope alone, which is "
                         "every number published before 2026-08-26")
    ap.add_argument("--near-beam-mrad", type=float, default=None,
                    help="replace the pots' measured HORIZONTAL aperture "
                         "with this half-width [mrad], keeping the measured "
                         "vertical: what a near-beam layer reaching closer "
                         "than the silicon package would buy (plans/09; "
                         "the Yellow Report high-acceptance 10 sigma_h "
                         "envelope is 2.2 / 1.8 / 0.92 mrad per "
                         "configuration, the tagging optics' 0.33 / 0.17 "
                         "/ 0.12 -- farforward.yr_optics, "
                         "tagging_optics).  Needs --rp-aperture measured")
    ap.add_argument("--envelope-split", type=float, default=0.0,
                    help="RELATIVE difference of the Roman-Pot cutout "
                         "half-width (along --split-axis) between the "
                         "m=+-1-rich and m=0-rich fills.  The spin-state "
                         "ratio cancels a COMMON cutout exactly; a "
                         "difference is the one systematic it cannot cancel "
                         "(code review F1).  0 = common")
    ap.add_argument("--split-axis", default="x", choices=("x", "y"),
                    help="axis of --envelope-split: x is the binding plane "
                         "under the tagging optics (the recoils escape "
                         "horizontally), y was the slot's (2026-08-28)")
    ap.add_argument("--u1-assumed", type=float, default=None,
                    help="u1 the ANALYSIS subtracts (default: the "
                         "generated value, i.e. exactly known)")
    ap.add_argument("--u2-assumed", type=float, default=None,
                    help="u2 the ANALYSIS subtracts.  An error here "
                         "reaches a_e at FIRST order as a_t du2 "
                         "<cos 2beta>, amplified by the cutout")
    ap.add_argument("--rel-lumi-offset", type=float, default=0.0,
                    help="relative-luminosity error unknown to the "
                         "analysis (second order in the ratio)")
    ap.add_argument("--u-in-situ", action="store_true",
                    help="measure (u1, u2) from the SPIN-AVERAGED counts "
                         "of the same data against the response's own "
                         "acceptance shape, instead of assuming them, and "
                         "propagate their covariance into the harmonic "
                         "errors (plans/08 A3).  With a free per-bin "
                         "acceptance u is not identifiable, so this "
                         "necessarily uses the acceptance MC -- see "
                         "reco.unpolarized_insitu_fit_2d")
    ap.add_argument("--fit", default="ratio", choices=("ratio", "likelihood"),
                    help="estimator of the two-azimuth harmonics.  'ratio' "
                         "(default, and every published number) inverts "
                         "the bin-wise spin-state ratio and fits it by "
                         "weighted LSQ; 'likelihood' profiles the per-bin "
                         "acceptance out of the Poisson likelihood, which "
                         "is exactly the conditional multinomial given the "
                         "bin totals and therefore unbiased at ANY count "
                         "-- the fix for the low-count bias of the sparse "
                         "|t| bins (plans/08 A12; run it with --ensemble "
                         "to see the difference)")
    ap.add_argument("--no-sin", action="store_true",
                    help="drop the sin 2alpha / sin 2beta / sin(alpha+beta) "
                         "null columns.  They are exactly forbidden by the "
                         "reflection symmetry of an unpolarized lepton on a "
                         "headless axis, are orthogonal to the cos columns "
                         "on the full grid (so they cost nothing), and "
                         "separate a spin-axis error from a Roman-Pot "
                         "azimuthal roll")
    ap.add_argument("--axis-tilt", type=float, default=0.0,
                    help="inject a spin-axis azimuth error [rad]: rotates "
                         "BOTH tensor harmonics, so both sin ratios return "
                         "tan 2delta")
    ap.add_argument("--pot-roll", type=float, default=0.0,
                    help="inject a Roman-Pot azimuthal roll [rad]: rotates "
                         "only the RECOIL harmonic, so sin 2alpha stays "
                         "null -- the signature that separates the two")
    ap.add_argument("--exact", action="store_true",
                    help="use the exact expected counts instead of a Poisson "
                         "draw: the systematic shifts (--envelope-split, "
                         "--u2-assumed, --rel-lumi-offset) then come out "
                         "free of statistical noise")
    ap.add_argument("--n-alpha", type=int, default=12)
    ap.add_argument("--n-beta", type=int, default=24,
                    help="(alpha, beta) binning of the two-azimuth fit; the "
                         "bin-wise ratio of Poisson counts is biased below "
                         "~30 counts per bin (--fit ratio, the default). "
                         "Coarser bins only attenuate that, at a cost in "
                         "resolution; --fit likelihood removes it outright "
                         "at any count and needs no re-binning "
                         "(plans/08 A12)")
    ap.add_argument("--t-edges", default=None,
                    help="comma-separated reconstructed |t| bin edges in "
                         "GeV^2, replacing the published seven-bin "
                         "0.017,0.028,0.039,0.05,0.08,0.12,0.17,0.25 "
                         "(recopseudo.T_EDGES_PUBLISHED, adopted "
                         "2026-08-28).  The run-13 window is "
                         "--t-edges 0.05,0.08,0.12,0.17,0.25 "
                         "(recopseudo.T_EDGES_LEGACY); a bin reaching "
                         "below 0.017 empties beta cells and is what "
                         "plans/08 8.4 rules out")
    ap.add_argument("--n-mc", type=int, default=600000,
                    help="response recoils; the published Table 3 uses "
                         "6e6 so that the template basis's own MC "
                         "statistics stay below the ten-year errors")
    ap.add_argument("--seed", type=int, default=20260824)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    if not args.lumi_fraction > 0:
        ap.error("--lumi-fraction must be positive")
    config = beams.default_configs("6Li")[args.config]
    rng = np.random.default_rng(args.seed)
    sc = coh.CoherentScenario(amp=args.amp, eps_b0=args.eps_b0)
    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            run_share=args.lumi_fraction,
                            pol_ion_tensor=args.pzz)
    lumi_ratio = args.lumi_10yr / args.lumi_1yr
    lumi_scale = 1.0                        # luminosity of the optics vs HA
    optics_note = "proton-derived 73 microrad (legacy)"
    if args.optics == "tagging":
        top = reco.tagging_optics_point(config, slope_b=sc.slope_b)
        args.sigma_theta = top["sigma_x_eff"]
        args.aspect = top["sigma_y"] / top["sigma_x_eff"]
        args.cut_scale_x = args.cut_scale_y = 1.0
        # the pots FOLLOW the envelope in both planes (Report 1 Section
        # 6.1): the geometric aperture is the envelope itself, so the
        # measured silicon aperture is not applied
        args.rp_aperture = "none"
        args.near_beam_mrad = None
        lumi_scale = top["lumi_fraction"]
        optics_note = ("tagging optics: horizontal beta* x %.0f, sigma_x "
                       "%.0f (eff. %.0f with dispersion) / sigma_y %.0f "
                       "microrad, L/L_HA = 1/%.1f, pots follow the 10 sigma "
                       "envelope %.2f x %.2f mrad"
                       % (top["r_h"], 1e6 * top["sigma_x"],
                          1e6 * top["sigma_x_eff"], 1e6 * top["sigma_y"],
                          1.0 / top["lumi_fraction"], 1e3 * top["env_x"],
                          1e3 * top["env_y"]))
    elif args.optics == "high-acceptance":
        sx, sy = reco.sigma_theta_for(config, "high-acceptance")
        args.sigma_theta, args.aspect = sx, sy / sx
        optics_note = ("Yellow Report high-acceptance divergence %.0f/%.0f "
                       "microrad" % (1e6 * sx, 1e6 * sy))
    print("optics:", optics_note)
    print("run plan: programme %g fb^-1/u/yr x share %g -> %g fb^-1/u "
          "delivered; optics L/L_HA = %.4f on top of it; spin-state share "
          "0.5 / 0.5 within it"
          % (args.lumi_1yr, args.lumi_fraction,
             args.lumi_1yr * args.lumi_fraction, lumi_scale))
    proj, n_coh, tagged = coh.project_coherent(
        config, scenario, sc, optics_list=(HIGH_ACCEPTANCE,),
        sigma_theta_list=(args.sigma_theta,))
    # coherent recoils produced in one year AT THIS OPTICS' luminosity
    n_produced = float(n_coh.sum()) * lumi_scale
    aperture = (reco.rp_aperture_for(config)
                if args.rp_aperture == "measured" else None)
    if args.optics == "tagging":
        aperture = (top["env_x"], top["env_y"])
    if args.rp_aperture == "measured" and aperture is None:
        raise SystemExit("no measured aperture for %s: it is tabulated per "
                         "machine configuration (reco.RP_APERTURE_MEASURED, "
                         "keyed by farforward.yr_config_key)"
                         % config.label())
    if args.near_beam_mrad is not None:
        if aperture is None:
            raise SystemExit("--near-beam-mrad replaces the horizontal half "
                             "of the measured aperture, so it needs "
                             "--rp-aperture measured")
        print("near-beam layer: horizontal aperture %.3f mrad (was %.3f), "
              "vertical unchanged at %.3f"
              % (args.near_beam_mrad, 1e3 * aperture[0], 1e3 * aperture[1]))
        aperture = (1e-3 * args.near_beam_mrad, aperture[1])
    cresp = rp.CoherentResponse(sc, config, args.sigma_theta,
                                aspect=args.aspect, shape=args.shape,
                                n_mc=args.n_mc, rng=rng,
                                cut_scale_xy=(args.cut_scale_x,
                                              args.cut_scale_y),
                                cut_theta_xy=aperture)
    n_tag = n_produced * cresp.acceptance
    plan = bk.tensor_flip_plan(args.pzz)
    pzz_list = [args.pzz, -2.0 * args.pzz]

    def a_t_func(t):   # cos 2beta coefficient per unit P_zz (Eq. 9 convention)
        return sc.cos2phi_coefficient_deformation(t, 1.0)

    # what the analysis ASSUMES, where that differs from the truth
    responses = None
    if args.envelope_split:
        # the perturbation acts on the BINDING vertical half-height
        # (envelope or measured aperture, whichever the cut is), so it is
        # never a silent no-op (2026-08-28)
        k = 0 if args.split_axis == "x" else 1
        scl = [1.0, 1.0]
        scl[k] = 1.0 + args.envelope_split
        responses = [cresp.with_cut(eff_scale_xy=tuple(scl)), cresp]
        print("fill-dependent Roman-Pot cutout: %s half-width %+.2e relative "
              "on the +Pzz fill (%.4f vs %.4f GeV), acceptance %.5f vs %.5f; "
              "the ratio cancels only a COMMON cutout"
              % ("horizontal" if k == 0 else "vertical", args.envelope_split,
                 responses[0].cut_pt_xy[k], cresp.cut_pt_xy[k],
                 responses[0].acceptance, cresp.acceptance))
    u_assumed = None
    if args.u_in_situ:
        if args.u1_assumed is not None or args.u2_assumed is not None:
            raise SystemExit("--u-in-situ measures (u1, u2) from the data; "
                             "it cannot be combined with --u1/--u2-assumed")
        u_assumed = "in-situ"
        print("unpolarized harmonics (u1, u2) measured IN SITU from the "
              "spin-averaged counts against the response's acceptance "
              "shape, and propagated into the harmonic errors")
    elif args.u1_assumed is not None or args.u2_assumed is not None:
        u_assumed = (args.u1 if args.u1_assumed is None else args.u1_assumed,
                     args.u2 if args.u2_assumed is None else args.u2_assumed)
        print("assumed unpolarized harmonics (u1, u2) = (%.4f, %.4f) "
              "against generated (%.4f, %.4f)"
              % (u_assumed[0], u_assumed[1], args.u1, args.u2))
    lumi_assumed = None
    if args.rel_lumi_offset:
        lumi_assumed = [0.5 * (1.0 + args.rel_lumi_offset),
                        0.5 * (1.0 - args.rel_lumi_offset)]
        print("assumed luminosity shares %s against equal truth"
              % (["%.6f" % v for v in lumi_assumed],))

    # azimuthal misalignment: a rotation of the tensor modulation, which
    # is what the sin columns are there to catch (rotating phi_S in the
    # response injects nothing -- it only relabels beta)
    d_e, d_t = args.axis_tilt, args.axis_tilt + args.pot_roll
    amp_c = args.amp * np.cos(2.0 * d_e)
    kw_sin = dict(with_sin=not args.no_sin,
                  a_e_s=args.amp * np.sin(2.0 * d_e),
                  a_t_s_func=(lambda t: a_t_func(t) * np.sin(2.0 * d_t)))
    if d_e or d_t:
        print("injected misalignment: spin axis %+.4f rad, pot roll %+.4f "
              "rad -> expected sin/cos ratios %.5f (alpha) and %.5f (beta)"
              % (args.axis_tilt, args.pot_roll, np.tan(2 * d_e),
                 np.tan(2 * d_t)))

    def a_t_rot(t):
        return a_t_func(t) * np.cos(2.0 * d_t)

    t_edges = t_edges_for(args)
    fits1, fits10, kept_edges = [], [], []
    for tlo, thi in zip(t_edges[:-1], t_edges[1:]):
        # A cutout tight enough to empty part of the circle leaves the
        # harmonic templates linearly dependent in that t bin.  That is a
        # statement about the acceptance, so report it and carry on with
        # the bins that survive rather than losing the whole figure --
        # which is what happens under --rp-aperture measured.
        try:
            f1 = rp.measure_coherent(
                cresp, n_produced, plan, tlo, thi, amp_c, a_t_rot,
                a_m=args.a_m, u1=args.u1, u2=args.u2, rng=rng,
                responses=responses, u_coeffs_assumed=u_assumed,
                lumi_assumed=lumi_assumed, poisson=not args.exact,
                n_alpha=args.n_alpha, n_beta=args.n_beta, fit=args.fit,
                **kw_sin)
            f10 = rp.measure_coherent(
                cresp, n_produced * lumi_ratio, plan, tlo, thi, amp_c,
                a_t_rot, a_m=args.a_m, u1=args.u1, u2=args.u2, rng=rng,
                responses=responses, u_coeffs_assumed=u_assumed,
                lumi_assumed=lumi_assumed, poisson=not args.exact,
                n_alpha=args.n_alpha, n_beta=args.n_beta, fit=args.fit,
                **kw_sin)
        except np.linalg.LinAlgError as exc:
            n_acc = int(((cresp.t_reco >= tlo)
                         & (cresp.t_reco < thi)).sum())
            print("t in [%g,%g): DROPPED -- %d accepted response "
                  "recoils, and %s" % (tlo, thi, n_acc, exc))
            continue
        fits1.append(f1)
        fits10.append(f10)
        kept_edges.append((tlo, thi))
    if not fits1:
        window = ("the cutout leaves NO accepted recoil at all"
                  if cresp.t_reco.size == 0 else
                  "the cutout leaves |t| = %.3f to %.3f GeV^2, outside the "
                  "%.3f-%.3f window binned here"
                  % (float(cresp.t_reco.min()), float(cresp.t_reco.max()),
                     t_edges[0], t_edges[-1]))
        raise SystemExit(
            "no t bin survived: with this acceptance the two-azimuth "
            "harmonic fit cannot be done as specified -- %s.  Loosen the "
            "cutout, re-bin |t| inside the window it leaves, or "
            "importance-sample above it with CoherentResponse(t_floor=...), "
            "which is what the WP5 optics scan does for exactly this "
            "reason." % window)

    # --- ensemble: is the template basis unbiased? ------------------------
    if args.ensemble > 0:
        print("ensemble of %d one-year pseudo-experiments per |t| bin "
              "(same response, fresh Poisson draws, %s fit): mean fit - "
              "injected, the spread, and the pull of the mean"
              % (args.ensemble, args.fit))
        for (tlo, thi) in kept_edges:
            at, ae_, dt, de = [], [], [], []
            for k in range(args.ensemble):
                rk = np.random.default_rng(args.seed + 1000 + k)
                f = rp.measure_coherent(
                    cresp, n_produced, plan, tlo, thi, amp_c, a_t_rot,
                    a_m=args.a_m, u1=args.u1, u2=args.u2, rng=rk,
                    responses=responses, u_coeffs_assumed=u_assumed,
                    lumi_assumed=lumi_assumed, n_alpha=args.n_alpha,
                    n_beta=args.n_beta, fit=args.fit, **kw_sin)
                at.append(f["a_t"]); ae_.append(f["a_e"])
                dt.append(f["err_t"]); de.append(f["err_e"])
            tr = f["truth"]
            at, ae_, dt, de = map(np.array, (at, ae_, dt, de))
            print("  t in [%g,%g): a_t mean %.4f (inj %.4f) spread %.4f "
                  "quoted err %.4f pull-of-mean %+.1f | a_e mean %.4f "
                  "(inj %.4f) spread %.4f quoted err %.4f pull-of-mean %+.1f"
                  % (tlo, thi, at.mean(), tr["a_t"], at.std(), dt.mean(),
                     (at.mean() - tr["a_t"]) / (at.std() / np.sqrt(len(at))),
                     ae_.mean(), tr["a_e"], ae_.std(), de.mean(),
                     (ae_.mean() - tr["a_e"]) / (ae_.std() / np.sqrt(len(ae_)))))

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11.8, 8.6))
    f0 = fits1[0]
    ae, be = f0["alpha_edges"], f0["beta_edges"]
    bc = 0.5 * (be[:-1] + be[1:])
    ac = 0.5 * (ae[:-1] + ae[1:])
    # the displayed modulation uses what the ANALYSIS assumes (u, shares)
    if args.u_in_situ:
        u1_shown, u2_shown = f0["u_insitu"]  # the summary loop prints all bins
    elif u_assumed is not None:
        u1_shown, u2_shown = u_assumed
    else:
        u1_shown, u2_shown = args.u1, args.u2
    lum_shown = lumi_assumed if lumi_assumed is not None else [0.5, 0.5]

    # --- (a) raw yields per fill vs beta --------------------------------
    for f, (col, lab) in enumerate(((C_TRUTH, r"$P_{zz}=+%.1f$ fill" % args.pzz),
                                    (C_FIT, r"$P_{zz}=-%.1f$ fill" % (2 * args.pzz)))):
        yb = f0["counts"][f].sum(axis=0)
        ax1.errorbar(bc, yb / yb.mean(), yerr=np.sqrt(yb) / yb.mean(), fmt="o-",
                     color=col, ms=3.5, lw=1, capsize=2, label=lab)
    fake = np.mean(np.cos(2.0 * cresp.beta_reco))
    ax1.annotate(r"cutout $|p_x|<%.2f$, $|p_y|<%.2f$ GeV: "
                 r"$\langle\cos2\beta\rangle=%.2f$ of the tagged sample" "\n"
                 r"(a single-fill fit would report $a_t\approx%.2f$ vs "
                 r"truth %.3f)"
                 % (cresp.cut_pt_xy[0], cresp.cut_pt_xy[1], fake,
                    fake / args.pzz, f0["truth"]["a_t"]),
                 xy=(0.03, 0.05), xycoords="axes fraction", fontsize=7.2)
    ax1.set_xlim(0, 2 * np.pi)
    ax1.set_xticks([0, np.pi, 2 * np.pi])
    ax1.set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
    ax1.set_xlabel(r"$\beta=\phi_t-\phi_S$ (reconstructed)")
    ax1.set_ylabel(r"$N_f(\beta)/\langle N_f\rangle$")
    # three decimals: the published window's lowest bin is 0.017-0.028,
    # which two would round to the same pair of labels
    ax1.set_title(r"(a) raw spin-sorted yields, $|t|\in[%.3f,%.3f]$: the "
                  r"cutout dominates" % kept_edges[0], fontsize=9)
    ax1.legend(fontsize=7, loc="upper right")
    ax1.tick_params(labelsize=8)

    # --- (b) acceptance-free modulation, projections -------------------
    r, var, sig2, pbar = reco.spin_state_ratio(
        f0["counts"].reshape(2, -1), lum_shown, pzz_list)
    aa, bb = np.meshgrid(ac, bc, indexing="ij")
    bm = f0["beta_means"]
    u = reco.unpolarized_modulation_2d(ae, be, u1_shown, u2_shown, beta_means=bm)
    t2d, var2d = reco._ratio_to_modulation(r, var, sig2, pbar, u=u)
    t2d = t2d.reshape(aa.shape)
    var2d = var2d.reshape(aa.shape)
    # projections over LIVE bins only (bins inside the cutout are empty):
    # remove the fitted a_t template per bin before averaging over beta,
    # so that the alpha projection isolates a_e cos 2alpha
    live = np.isfinite(var2d)
    resid = t2d - f0["const"] - f0["a_t"] * bm["c2t"][None, :]
    nb_live = np.maximum(live.sum(axis=0), 1)
    tb = np.where(live, t2d - f0["const"], 0.0).sum(axis=0) / nb_live
    eb = np.sqrt(np.where(live, var2d, 0.0).sum(axis=0)) / nb_live
    na_live = np.maximum(live.sum(axis=1), 1)
    ta = np.where(live, resid, 0.0).sum(axis=1) / na_live
    ea = np.sqrt(np.where(live, var2d, 0.0).sum(axis=1)) / na_live
    tb = np.where(live.sum(axis=0) > 0, tb, np.nan)
    ax2.errorbar(bc, tb, yerr=eb, fmt="o", color=C_TRUTH, ms=3.5, capsize=2,
                 lw=1, label=r"$T(\beta)$ ($\alpha$-averaged)")
    ax2.errorbar(ac, ta, yerr=ea, fmt="s", color=C_FIT, ms=3.5, capsize=2,
                 lw=1, label=r"$T(\alpha)$ ($\beta$-averaged, $a_t$ offset removed)")
    phi = np.linspace(0, 2 * np.pi, 200)
    ax2.plot(phi, f0["a_t"] * np.cos(2 * phi), "-", color=C_TRUTH, lw=1.0,
             alpha=0.5, label=r"$a_t\cos2\beta$ (smooth)")
    ax2.plot(bc, f0["a_t"] * bm["c2t"], "_", color=C_TRUTH, ms=9, mew=1.6,
             label=r"template fit ($a_t\propto|t|$, MC basis): "
                   r"$a_t(t_{\rm ref})=%.4f\pm%.4f$" % (f0["a_t"], f0["err_t"]))
    ax2.plot(phi, f0["a_e"] * np.cos(2 * phi), "-", color=C_FIT, lw=1.4,
             label=r"fit $a_e\cos2\alpha$: $a_e=%.4f\pm%.4f$"
                   % (f0["a_e"], f0["err_e"]))
    ax2.set_xlim(0, 2 * np.pi)
    ax2.set_xticks([0, np.pi, 2 * np.pi])
    ax2.set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
    ax2.set_xlabel(r"$\beta$ or $\alpha$")
    ax2.set_ylabel(r"$T-\hat\kappa$ (per unit $P_{zz}$)")
    ax2.set_title(r"(b) spin-state ratio inverted for the modulation "
                  r"(1 yr, $N=%.2g$)" % f0["n"], fontsize=9)
    ax2.legend(fontsize=6.8, loc="upper right")
    ax2.tick_params(labelsize=8)
    ax2.axhline(0, color="0.85", lw=0.6, zorder=0)

    # --- (c) expected 2-D modulation map -------------------------------
    mu = f0["expected"]
    r0, v0, s0, p0 = reco.spin_state_ratio(mu.reshape(2, -1), lum_shown,
                                           pzz_list)
    t0, _ = reco._ratio_to_modulation(r0, v0, s0, p0, u=u)
    im = ax3.imshow(t0.reshape(aa.shape).T, origin="lower",
                    extent=(0, 2 * np.pi, 0, 2 * np.pi), cmap="RdBu_r",
                    aspect="auto")
    fig.colorbar(im, ax=ax3, fraction=0.046, pad=0.03).set_label(
        r"$T(\alpha,\beta)$", fontsize=8)
    ax3.set_xlabel(r"$\alpha=\phi_e-\phi_S$")
    ax3.set_ylabel(r"$\beta=\phi_t-\phi_S$")
    ax3.set_xticks([0, np.pi, 2 * np.pi])
    ax3.set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
    ax3.set_yticks([0, np.pi, 2 * np.pi])
    ax3.set_yticklabels(["0", r"$\pi$", r"$2\pi$"])
    ax3.set_title(r"(c) expected modulation map: $a_t\cos2\beta$ (rows) + "
                  r"$a_e\cos2\alpha$ (columns)", fontsize=9)
    ax3.tick_params(labelsize=8)

    # --- (d) extracted coefficients vs t ---------------------------------
    tc = np.array([f["truth"]["t_ref"] for f in fits1])
    tt = np.linspace(0.0, 0.27, 100)
    ax4.plot(tt, a_t_func(tt), "-", color=C_TRUTH, lw=1.4,
             label=r"injected $a_t(t)=2a_2/P_{zz}=-\frac{\epsilon_{B0}}{2}B|t|$")
    ax4.plot(tc, [f["truth"]["a_t"] for f in fits1], "_", color=C_TRUTH, ms=10,
             mew=1.4, label=r"$a_t(t_{\rm ref})$, $t_{\rm ref}$ = bin mean true $|t|$")
    ax4.errorbar(tc, [f["a_t"] for f in fits1], yerr=[f["err_t"] for f in fits1],
                 fmt="s", mfc="none", color=C_TRUTH, ms=5, capsize=2, lw=1,
                 label=r"$a_t$, 1 yr")
    ax4.errorbar(tc * 1.02, [f["a_t"] for f in fits10],
                 yerr=[f["err_t"] for f in fits10], fmt="o", color=C_TRUTH,
                 ms=4, capsize=2, lw=1, label=r"$a_t$, 10 yr")
    ax4.axhline(args.amp, color=C_FIT, lw=1.2, label=r"injected $a_e$ (flat)")
    ax4.errorbar(tc, [f["a_e"] for f in fits1], yerr=[f["err_e"] for f in fits1],
                 fmt="s", mfc="none", color=C_FIT, ms=5, capsize=2, lw=1,
                 label=r"$a_e$, 1 yr")
    ax4.errorbar(tc * 1.02, [f["a_e"] for f in fits10],
                 yerr=[f["err_e"] for f in fits10], fmt="o", color=C_FIT, ms=4,
                 capsize=2, lw=1, label=r"$a_e$, 10 yr")
    ax4.axhline(0, color="0.85", lw=0.6, zorder=0)
    ax4.set_xlim(0, 0.27)
    ax4.set_xlabel(r"$|t|$ [GeV$^2$] ($t_{\rm ref}$ = reco-bin mean of true $|t|$)")
    ax4.set_ylabel(r"coefficient per unit $P_{zz}$")
    ax4.set_title("(d) two-component extraction per reco $t$ bin", fontsize=9)
    ax4.legend(fontsize=6.6, loc="upper left")
    ax4.tick_params(labelsize=8)

    fig.suptitle(
        r"Coherent $e\,^6$Li$\to e'X\,^6$Li(g.s.) at the reconstructed level "
        r"(6R), %s, $P_{zz}=%.2f$: $\sigma_\theta$ = %.0f / %.0f $\mu$rad (h / v), "
        r"%s cutout $|p_x|<%.2f$, $|p_y|<%.2f$ GeV [%s]"
        "\n" r"$N_{\rm tag}$ = %.2g (1 yr) vs %.2g with the constant 0.20 GeV cut at $L_{HA}$; "
        r"deformation $c_2=2a_2$ (Eq. 9 convention), $a_e=%.3f$"
        "\n" r"$u_1=%.2f$, $u_2=%.2f$ (ZEUS LPS $1\sigma$ bounds); "
        r"two-fill %s fit with the MC template basis; statistical errors only"
        % (config.label(), args.pzz, 1e6 * args.sigma_theta,
           1e6 * args.sigma_theta * args.aspect,
           args.shape, cresp.cut_pt_xy[0], cresp.cut_pt_xy[1],
           args.optics + " optics", n_tag,
           tagged[HIGH_ACCEPTANCE.name].sum(), args.amp, args.u1, args.u2,
           args.fit),
        fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.905))
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / ("%s.png" % output_stem(args))
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print("coherent produced (1 yr at this optics' luminosity, L/L_HA = %.3f): "
          "%.3g; tagged: %s cutout |px|<%.3f, |py|<%.3f GeV = %.2f x %.2f mrad "
          "(sigma_theta h/v %.0f/%.0f microrad) -> acc %.4f, N_tag %.3g; "
          "constant-cut reference at L_HA %.3g"
          % (lumi_scale, n_produced, args.shape, cresp.cut_pt_xy[0],
             cresp.cut_pt_xy[1], 1e3 * cresp.cut_theta_eff[0],
             1e3 * cresp.cut_theta_eff[1], 1e6 * args.sigma_theta,
             1e6 * args.sigma_theta * args.aspect, cresp.acceptance, n_tag,
             tagged[HIGH_ACCEPTANCE.name].sum()))
    print("cutout fake <cos 2beta> = %.3f" % fake)
    # kept_edges, not a rebuilt contiguous list: a bin dropped from
    # anywhere but the tail would shift every label after it
    for (tlo, thi), f1, f10 in zip(kept_edges, fits1, fits10):
        tr = f1["truth"]
        print("t in [%g,%g]: N=%.3g  a_t truth %.4f  fit %.4f +- %.4f (1yr) "
              "%.4f +- %.4f (10yr) | a_e truth %.4f fit %.4f +- %.4f (1yr) "
              "%.4f +- %.4f (10yr) | a_m %.4f +- %.4f"
              % (tlo, thi, f1["n"], tr["a_t"], f1["a_t"], f1["err_t"],
                 f10["a_t"], f10["err_t"], tr["a_e"], f1["a_e"], f1["err_e"],
                 f10["a_e"], f10["err_e"], f1["a_m"], f1["err_m"]))
        # the design diagnostic the |t| window was chosen on (plans/08
        # 8.4): how many beta bins the acceptance leaves populated in
        # this |t| bin, the counts per (alpha, beta) cell that set the
        # ratio's low-count bias, and the condition number of the
        # weighted design the fit solves -- sigma_max/sigma_min of that
        # design is sqrt(cond) of the returned covariance, since
        # cov = (A^T A)^-1 has eigenvalues 1/sigma_i^2.  `rank` is the
        # rank reco._harmonic_rank_guard MEASURED on that weighted
        # design (fit key "design_rank"), not the parameter count
        # restated twice -- it printed n/n unconditionally until the
        # 2026-08-28 review
        mu1 = np.asarray(f1["expected"])
        cov1 = np.asarray(f1.get("cov_stat", f1["cov"]))
        n_beta_live = int((mu1.sum(axis=(0, 1)) > 0.0).sum())
        n_cells = mu1.shape[1] * mu1.shape[2]
        print("    design: %.0f counts per (alpha, beta) cell, %d of %d "
              "beta bins populated, cond %.2f, rank %d/%d"
              % (mu1.sum() / n_cells, n_beta_live, mu1.shape[2],
                 np.sqrt(np.linalg.cond(cov1)),
                 f1["design_rank"], f1["n_par"]))
        if args.u_in_situ:
            print("    in-situ (u1, u2) = (%.4f +- %.4f, %.4f +- %.4f) "
                  "against generated (%.4f, %.4f); propagated into a_e it "
                  "adds %.5f in quadrature to the %.5f statistical error"
                  % (f1["u_insitu"][0], f1["u_err"][0], f1["u_insitu"][1],
                     f1["u_err"][1], args.u1, args.u2,
                     np.sqrt(max(f1["cov_u"][1, 1], 0.0)),
                     np.sqrt(f1["cov_stat"][1, 1])))
        if not args.no_sin:
            print("    null test: a_e_s %+.4f +- %.4f, a_t_s %+.4f +- %.4f, "
                  "a_m_s %+.4f +- %.4f  -> sin/cos %+.4f (alpha) %+.4f "
                  "(beta); axis resolution %.0f / %.0f mrad (1 yr)"
                  % (f1["a_e_s"], f1["err_e_s"], f1["a_t_s"], f1["err_t_s"],
                     f1["a_m_s"], f1["err_m_s"],
                     f1["a_e_s"] / f1["a_e"], f1["a_t_s"] / f1["a_t"],
                     5e2 * f1["err_e_s"] / abs(f1["a_e"]),
                     5e2 * f1["err_t_s"] / abs(f1["a_t"])))
    # a_e is ONE constant across the |t| bins, so its errors combine; a_t
    # is a different number in each bin and does not (plans/08 8.4)
    de1 = np.array([f["err_e"] for f in fits1])
    print("combined one-year delta(a_e) over the %d |t| bins: %.5f (from %s)"
          % (len(de1), 1.0 / np.sqrt((1.0 / de1 ** 2).sum()),
             " / ".join("%.5f" % v for v in de1)))


if __name__ == "__main__":
    main()
