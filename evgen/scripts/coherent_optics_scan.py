#!/usr/bin/env python3
"""WP5: the coherent tag versus the near-beam envelope (plans/07 WP5,
plans/08 A4).

The projection report quotes the coherent acceptance at two points
(0.20 GeV high-acceptance, 0.41-0.45 GeV high-divergence) and the
reconstruction note adds that the envelope is an ANGLE, so the same
optics gives an angular cut of 10 sigma_theta A p_u on the 6Li recoil at
the three machine configurations (40.8 / 99.5 / 137.5 GeV/u -- gamma-
matched at the two lower ones, plans/10).  A referee will ask for a curve, not two points; the Li
optics are undocumented (plans/04 #11, #20), so a curve is also the
honest form of the statement.  This script produces it:

  (a) tag acceptance versus the near-beam cut, analytically from
      `reco.rp_hole_acceptance` -- exact for any cutout shape and any
      cut, with no Monte Carlo, over the B = 40-60 GeV^-2 band, for the
      slot-like ePIC cutout, a square, and the circular cut the routing
      code uses;
  (b) tagged yield in one year, with exp(-B t_min) and the 0.73
      rate-weighting systematic folded into the central curve, against
      the IR-8 secondary-focus alternative;
  (c) the statistical error on the deformation coefficient a_t and on
      the gluon-transversity coefficient a_e from the full response
      (Roman-Pot emulation, divergence smearing, template fit), with the
      |t| bin moved with the cut.  This needs the importance sampling of
      `CoherentResponse(t_floor=...)`: above ~0.3 GeV the plain sampler
      leaves no accepted recoils at all;
  (d) acceptance versus beam momentum per nucleon at fixed optics -- the
      statement that the coherent program chooses the beam energy.

Usage:  python3 scripts/coherent_optics_scan.py
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
from polli_fastsim import farforward as _ff  # noqa: E402
from polli_fastsim.farforward import HIGH_ACCEPTANCE  # noqa: E402

C_TRUTH, C_FIT, C_ALT, C_GREY = "#0072B2", "#D55E00", "#009E73", "0.45"

# eSTARlight IR-8 secondary-focus intact-recoil efficiency x acceptance
# (plans/06 SS6.5): no 6Li entry exists; the interpolation is ours.
IR8_PUBLISHED = {"d": 0.47, "3He": 0.32, "4He": 0.29, "7Li": 0.178}

# Beam configurations are DERIVED, never hard-coded: the two lower 6Li
# energies moved from rigidity-scaled (20.5, 50) to gamma-matched
# (40.8, 99.5) GeV/u on 2026-08-27 (plans/10).
_PU = tuple(c.ion_momentum_per_nucleon for c in beams.default_configs("6Li"))
IR8_LI6_INTERPOLATED = 0.20


def acceptance_curve(cuts, slope_b, aspect_x, aspect_y, shape="rectangle"):
    """Analytic tagged fraction and fake <cos 2beta> versus the vertical
    half-height `cut`, for a cutout of half-widths (aspect_x, aspect_y)
    times it.  `reco.rp_hole_acceptance` integrates exp(-B rho(phi)^2)
    over the cutout boundary exactly."""
    acc, a2 = [], []
    for c in cuts:
        h = reco.rp_hole_acceptance(slope_b, aspect_x * c, aspect_y * c,
                                    shape=shape)
        acc.append(h["acc"])
        a2.append(h["a2"])
    return np.array(acc), np.array(a2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2))
    ap.add_argument("--lumi-1yr", type=float, default=10.0)
    ap.add_argument("--pzz", type=float, default=0.60)
    ap.add_argument("--amp", type=float, default=0.01)
    ap.add_argument("--eps-b0", type=float, default=-0.08)
    ap.add_argument("--u1", type=float, default=0.05)
    ap.add_argument("--u2", type=float, default=0.02)
    ap.add_argument("--cut-scale-x", type=float, default=2.5,
                    help="cutout half-width in x in units of the vertical "
                         "half-height (the ePIC horizontal slot)")
    ap.add_argument("--n-mc", type=int, default=300000)
    ap.add_argument("--err-cuts", type=float, nargs="*",
                    default=[0.10, 0.15, 0.22, 0.30, 0.45, 0.60],
                    help="near-beam cuts [GeV] at which the fitted errors "
                         "are evaluated with the full response")
    ap.add_argument("--seed", type=int, default=20260825)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs("6Li")[args.config]
    outdir = pathlib.Path(args.outdir)
    sc = coh.CoherentScenario(amp=args.amp, eps_b0=args.eps_b0)
    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    plan = bk.tensor_flip_plan(args.pzz)
    _proj, n_coh, _tagged = coh.project_coherent(
        config, scenario, sc, optics_list=(HIGH_ACCEPTANCE,))
    n_produced = float(n_coh.sum())
    p_ion = config.ion.A * config.ion_momentum_per_nucleon
    cut_here = reco.tag_pt_cut(HIGH_ACCEPTANCE.sigma_theta,
                               config.ion_momentum_per_nucleon,
                               a_beam=config.ion.A)

    print("configuration %s: coherent recoils produced in %g fb^-1/u = %.3g"
          % (config.label(), args.lumi_1yr, n_produced))
    cut_yr = _ff.yr_optics(config).envelope[0] * p_ion
    cut_tag = _ff.tagging_optics_point(config)["env_x"] * p_ion
    print("angular envelope at this energy (10 sigma_h x A p_u): %.3f GeV at "
          "the Yellow Report high-acceptance optics, %.3f GeV at the tagging "
          "optics of Report 1 Section 6.1, %.3f GeV at the legacy 73 urad "
          "(reproduction only)" % (cut_yr, cut_tag, cut_here))
    # The pots' own GEOMETRIC aperture, measured in the ePIC geometry
    # (reco.RP_APERTURE_MEASURED, tools/fullsim).  It is a second, larger
    # constraint at every configuration, and the curves below are read at
    # it rather than at the envelope wherever it is the binding one.
    ap = reco.rp_aperture_for(config.ion_momentum_per_nucleon)
    cut_meas = None if ap is None else float(ap[0]) * p_ion
    if cut_meas is not None:
        print("measured pot aperture at this energy: %.3f GeV in p_T "
              "(|theta_x| > %.2f mrad), i.e. %.2fx the Yellow Report envelope "
              "and %.1fx the tagging envelope"
              % (cut_meas, 1e3 * ap[0], cut_meas / cut_yr, cut_meas / cut_tag))

    cuts = np.linspace(0.05, 0.70, 140)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11.8, 8.4))

    # --- (a) acceptance vs the near-beam cut ------------------------------
    for b, style in ((40.0, ":"), (50.0, "-"), (60.0, "--")):
        acc, _ = acceptance_curve(cuts, b, args.cut_scale_x, 1.0)
        ax1.plot(cuts, acc, style, color=C_TRUTH, lw=1.4,
                 label=(r"slot %.1f$\times$1, $B$ = %g GeV$^{-2}$"
                        % (args.cut_scale_x, b)) if style == "-"
                 else r"$B$ = %g" % b)
    acc_sq, _ = acceptance_curve(cuts, 50.0, 1.0, 1.0)
    ax1.plot(cuts, acc_sq, "-", color=C_FIT, lw=1.4, label="square, $B$ = 50")
    ax1.plot(cuts, np.exp(-50.0 * cuts ** 2), "-", color=C_GREY, lw=1.2,
             label=r"circular $e^{-B p_T^2}$ (routing code)")
    # the horizontal 10 sigma envelope of each configuration at the Yellow
    # Report high-acceptance optics (plans/10; 0.54 / 1.07 / 0.76 GeV) and
    # at the tagging optics of Report 1 Section 6.1 (0.08 / 0.10 / 0.10 GeV)
    for cfg_i, lab in zip(beams.default_configs("6Li"),
                          ("5$\\times$41", "10$\\times$100", "18$\\times$275")):
        pi = cfg_i.ion.A * cfg_i.ion_momentum_per_nucleon
        c_yr = _ff.yr_optics(cfg_i).envelope[0] * pi
        c_tag = _ff.tagging_optics_point(cfg_i)["env_x"] * pi
        for c, ls, col in ((c_yr, "-", "0.75"), (c_tag, ":", C_ALT)):
            if cuts[0] <= c <= cuts[-1]:
                ax1.axvline(c, color=col, lw=0.9, ls=ls, zorder=0)
                ax1.text(c, 1.4e-3, " " + lab + (" YR" if ls == "-" else " tagging"),
                         fontsize=6.0, rotation=90, color="0.35", va="bottom")
    if cut_meas is not None:
        ax1.axvline(cut_meas, color=C_ALT, lw=1.3, ls="-.")
        ax1.text(cut_meas, 1.4e-3, " measured pot aperture", fontsize=6.6,
                 rotation=90, color=C_ALT, va="bottom")
    ax1.set_yscale("log")
    ax1.set_ylim(1e-3, 1.2)
    ax1.set_xlim(cuts[0], cuts[-1])
    ax1.set_xlabel(r"near-beam envelope (vertical half-height) [GeV]")
    ax1.set_ylabel("tagged fraction of coherent recoils")
    ax1.set_title("(a) analytic tag acceptance vs the envelope", fontsize=9)
    ax1.legend(fontsize=6.6, loc="lower left")
    ax1.tick_params(labelsize=8)

    # --- (b) tagged yield, with t_min and the rate weighting --------------
    acc50, _ = acceptance_curve(cuts, 50.0, args.cut_scale_x, 1.0)
    # t_min = M_A^2 x_P^2/(1-x_P) at the top of the coherent x_P window
    # (x_P = 0.01, M_A = 5.6 GeV): coherent.recoil_lab neglects it, and it
    # suppresses the tagged yield by exp(-B t_min) there (plans/07 WP5)
    t_min = 3.2e-3
    supp = float(np.exp(-50.0 * t_min))
    ax2.plot(cuts, n_produced * acc50, "-", color=C_TRUTH, lw=1.6,
             label=r"$N_{\rm tag}$, slot, $B$ = 50")
    ax2.plot(cuts, n_produced * acc50 * supp, "--", color=C_TRUTH, lw=1.2,
             label=r"$\times e^{-B t_{\min}}$ (%.0f%%)" % (100 * supp))
    ax2.plot(cuts, n_produced * acc50 * supp * coh.RATE_WEIGHT_SYST, ":",
             color=C_TRUTH, lw=1.2,
             label=r"$\times$ %.2f rate weighting" % coh.RATE_WEIGHT_SYST)
    ax2.axhline(n_produced * IR8_LI6_INTERPOLATED, color=C_ALT, lw=1.4,
                label="IR-8 secondary focus, %.0f%% (interpolated, ours)"
                      % (100 * IR8_LI6_INTERPOLATED))
    lo = n_produced * min(IR8_PUBLISHED.values())
    hi = n_produced * max(IR8_PUBLISHED.values())
    ax2.axhspan(lo, hi, color=C_ALT, alpha=0.12, zorder=0)
    ax2.text(0.28, hi, "  IR-8 published: "
             + ", ".join("%s %.0f%%" % (k, 100 * v)
                         for k, v in IR8_PUBLISHED.items()),
             fontsize=6.2, color=C_ALT, va="bottom")
    ax2.axvline(cut_tag, color=C_ALT, lw=1.0, ls=":")
    if cut_meas is not None:
        ax2.axvline(cut_meas, color=C_ALT, lw=1.3, ls="-.")
    ax2.set_yscale("log")
    ax2.set_xlim(cuts[0], cuts[-1])
    ax2.set_xlabel("near-beam envelope [GeV]")
    ax2.set_ylabel(r"tagged recoils in %g fb$^{-1}$/u" % args.lumi_1yr)
    ax2.set_title("(b) tagged yield, with $t_{\\min}$ and the rate weighting",
                  fontsize=9)
    ax2.legend(fontsize=6.4, loc="lower left")
    ax2.tick_params(labelsize=8)

    # --- (c) fitted errors from the full response -------------------------
    print("\n  cut [GeV]   acc      N_tag      t bin        "
          "a_t(t_ref)   d(a_t)     d(a_e)")
    rows = []
    for cut in args.err_cuts:
        sigma_theta = cut / (10.0 * p_ion)
        tlo, thi = 1.05 * cut ** 2, 1.05 * cut ** 2 + 0.03
        cresp = rp.CoherentResponse(sc, config, sigma_theta,
                                    cut_scale_xy=(args.cut_scale_x, 1.0),
                                    n_mc=args.n_mc,
                                    rng=np.random.default_rng(args.seed),
                                    t_floor=0.25 * cut ** 2)
        fit = rp.measure_coherent(
            cresp, n_produced, plan, tlo, thi, args.amp,
            lambda t: sc.cos2phi_coefficient_deformation(t, 1.0),
            u1=args.u1, u2=args.u2, poisson=False)
        tr = fit["truth"]
        rows.append((cut, cresp.acceptance, n_produced * cresp.acceptance,
                     tr["a_t"], fit["err_t"], fit["err_e"]))
        print("  %-11.2f %-8.2e %-10.3g [%.3f,%.3f]  %-12.4f %-10.4f %.4f"
              % (cut, cresp.acceptance, n_produced * cresp.acceptance,
                 tlo, thi, tr["a_t"], fit["err_t"], fit["err_e"]))
    rows = np.array(rows)
    ax3.errorbar(rows[:, 0], rows[:, 4] / rows[:, 3], fmt="o-", color=C_TRUTH,
                 ms=4, lw=1.3,
                 label=r"$\delta a_t / a_t$ (deformation)")
    ax3.plot(rows[:, 0], rows[:, 5] / args.amp, "s--", color=C_FIT, ms=4,
             lw=1.3, label=r"$\delta a_e / a_e$ (gluon transversity)")
    ax3.axhline(1.0, color="0.8", lw=0.8, zorder=0)
    ax3.axvline(cut_tag, color=C_ALT, lw=1.0, ls=":")
    if cut_meas is not None:
        ax3.axvline(cut_meas, color=C_ALT, lw=1.3, ls="-.")
    ax3.set_yscale("log")
    ax3.set_xlabel("near-beam envelope [GeV]")
    ax3.set_ylabel("relative statistical error, 1 yr")
    ax3.set_title("(c) reach of the two coefficients vs the envelope",
                  fontsize=9)
    ax3.legend(fontsize=7, loc="upper left")
    ax3.tick_params(labelsize=8)

    # --- (d) acceptance vs beam momentum ----------------------------------
    p_us = np.linspace(10.0, 145.0, 200)
    # the Yellow Report divergences are per configuration, not a function
    # of the momentum: each curve is one configuration's high-acceptance
    # (sigma_h, sigma_v) swept in p_u, with the configuration's own point
    # marked, and the tagging optics of Report 1 Section 6.1 alongside
    for cfg_i, col in zip(beams.default_configs("6Li"), (C_TRUTH, C_FIT, C_ALT)):
        sh, sv = _ff.sigma_theta_for(cfg_i)
        acc_p = np.array([reco.rp_hole_acceptance(50.0, 10 * sh * config.ion.A * p,
                                                  10 * sv * config.ion.A * p)["acc"]
                          for p in p_us])
        ax4.plot(p_us, acc_p, "-", color=col, lw=1.5,
                 label="YR HA %.0f/%.0f $\\mu$rad (%s)" % (1e6 * sh, 1e6 * sv, cfg_i.label()))
        pu_i = cfg_i.ion_momentum_per_nucleon
        ax4.plot([pu_i], [reco.rp_hole_acceptance(50.0, 10 * sh * config.ion.A * pu_i,
                                                  10 * sv * config.ion.A * pu_i)["acc"]],
                 "o", color=col, ms=6)
        t = _ff.tagging_optics_point(cfg_i)
        ax4.plot([pu_i], [t["acceptance"]], "^", color=col, ms=7, mfc="white", mew=1.6)
    ax4.plot([], [], "^", color="0.3", mfc="white", mew=1.6, label="tagging optics (pots follow)")
    for p_u in _PU:
        ax4.axvline(p_u, color="0.75", lw=0.9, zorder=0)
    ax4.set_yscale("log")
    ax4.set_ylim(1e-30, 1.5)
    ax4.set_xlabel(r"$^6$Li beam momentum per nucleon [GeV]")
    ax4.set_ylabel("tagged fraction")
    ax4.set_title("(d) the coherent program chooses the beam energy",
                  fontsize=9)
    ax4.legend(fontsize=7, loc="lower left")
    ax4.tick_params(labelsize=8)

    fig.suptitle(
        "Coherent intact-$^6$Li tag versus the near-beam envelope (WP5), "
        "%s, $B$ = 50 GeV$^{-2}$, $f_0$ = %.2f\n"
        "the envelope is ANGULAR ($10\\sigma_\\theta A p_u$): %.2f GeV at the "
        "legacy 73 $\\mu$rad, %.2f GeV at the Yellow Report high-acceptance optics, "
        "%.2f GeV at the tagging optics"
        % (config.label(), sc.f0, cut_here,
           _ff.yr_optics(config).envelope[0] * p_ion,
           _ff.tagging_optics_point(config)["env_x"] * p_ion), fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out = outdir / "coherent_optics_scan_6Li.png"
    fig.savefig(out, dpi=140)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
