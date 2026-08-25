#!/usr/bin/env python3
"""Money plot 6R: the coherent intact-6Li channel at the RECONSTRUCTED
level (plans/07 WP3 -> WP5) -- money plot 6 re-derived with

  * the ANGULAR near-beam cut, pT_cut = 10 sigma_theta A p_u (proton-
    derived sigma_theta = 73 microrad, high-acceptance optics), and a
    slot-like Roman-Pot cutout (the ePIC sensor planes surround a
    horizontal slot for the beam's momentum spread and dispersion,
    Jentsch DIS 2023: default half-widths 2.5 x 10 sigma in x, 1 x 10
    sigma in y) whose sides are parallel to the vertical spin axis;
  * the Roman-Pot emulation: divergence smearing of the recoil angle,
    |t| = pT^2 reconstructed with x_L = 1, reco t bins;
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2))
    ap.add_argument("--lumi-1yr", type=float, default=10.0)
    ap.add_argument("--lumi-10yr", type=float, default=100.0)
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
    ap.add_argument("--sigma-theta", type=float, default=reco.SIGMA_THETA_HA,
                    help="beam angular divergence [rad] (proton-derived)")
    ap.add_argument("--aspect", type=float, default=1.0,
                    help="beam-divergence anisotropy sigma_y/sigma_x "
                         "(HERA: 100/45 MeV vertical/horizontal pT spread, "
                         "ZEUS NPB 816:1; Li optics undocumented, #20)")
    ap.add_argument("--cut-scale-x", type=float, default=2.5,
                    help="cutout half-width in x in units of 10 sigma_x: "
                         "the ePIC pots surround a horizontal SLOT (beam "
                         "momentum spread + dispersion), Jentsch DIS 2023 "
                         "slide 15 [illustrative value]")
    ap.add_argument("--cut-scale-y", type=float, default=1.0,
                    help="cutout half-height in y in units of 10 sigma_y")
    ap.add_argument("--shape", default="rectangle",
                    choices=("rectangle", "ellipse"))
    ap.add_argument("--envelope-split", type=float, default=0.0,
                    help="RELATIVE difference of the Roman-Pot vertical "
                         "half-height between the m=+-1-rich and m=0-rich "
                         "fills.  The spin-state ratio cancels a COMMON "
                         "cutout exactly; a difference is the one "
                         "systematic it cannot cancel (code review F1), "
                         "and the slot amplifies it far beyond the naive "
                         "d<cos 2beta>/(P+ - P0).  0 = common")
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
    ap.add_argument("--n-mc", type=int, default=600000)
    ap.add_argument("--seed", type=int, default=20260824)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs("6Li")[args.config]
    rng = np.random.default_rng(args.seed)
    sc = coh.CoherentScenario(amp=args.amp, eps_b0=args.eps_b0)
    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    lumi_ratio = args.lumi_10yr / args.lumi_1yr
    proj, n_coh, tagged = coh.project_coherent(
        config, scenario, sc, optics_list=(HIGH_ACCEPTANCE,),
        sigma_theta_list=(args.sigma_theta,))
    n_produced = float(n_coh.sum())          # coherent recoils, 1 yr
    cresp = rp.CoherentResponse(sc, config, args.sigma_theta,
                                aspect=args.aspect, shape=args.shape,
                                n_mc=args.n_mc, rng=rng,
                                cut_scale_xy=(args.cut_scale_x,
                                              args.cut_scale_y))
    n_tag = n_produced * cresp.acceptance
    plan = bk.tensor_flip_plan(args.pzz)
    pzz_list = [args.pzz, -2.0 * args.pzz]

    def a_t_func(t):   # cos 2beta coefficient per unit P_zz (Eq. 9 convention)
        return sc.cos2phi_coefficient_deformation(t, 1.0)

    # what the analysis ASSUMES, where that differs from the truth
    responses = None
    if args.envelope_split:
        responses = [cresp.with_cut((args.cut_scale_x,
                                     args.cut_scale_y
                                     * (1.0 + args.envelope_split))),
                     cresp]
        print("fill-dependent Roman-Pot envelope: vertical half-height "
              "%+.2e relative on the +Pzz fill (%.4f vs %.4f GeV); the "
              "ratio cancels only a COMMON cutout"
              % (args.envelope_split, responses[0].cut_pt_xy[1],
                 cresp.cut_pt_xy[1]))
    u_assumed = None
    if args.u1_assumed is not None or args.u2_assumed is not None:
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

    t_edges = [0.05, 0.08, 0.12, 0.17, 0.25]
    fits1, fits10 = [], []
    for tlo, thi in zip(t_edges[:-1], t_edges[1:]):
        fits1.append(rp.measure_coherent(
            cresp, n_produced, plan, tlo, thi, amp_c, a_t_rot, a_m=args.a_m,
            u1=args.u1, u2=args.u2, rng=rng, responses=responses,
            u_coeffs_assumed=u_assumed, lumi_assumed=lumi_assumed, **kw_sin))
        fits10.append(rp.measure_coherent(
            cresp, n_produced * lumi_ratio, plan, tlo, thi, amp_c, a_t_rot,
            a_m=args.a_m, u1=args.u1, u2=args.u2, rng=rng,
            responses=responses, u_coeffs_assumed=u_assumed,
            lumi_assumed=lumi_assumed, **kw_sin))

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11.8, 8.6))
    f0 = fits1[0]
    ae, be = f0["alpha_edges"], f0["beta_edges"]
    bc = 0.5 * (be[:-1] + be[1:])
    ac = 0.5 * (ae[:-1] + ae[1:])

    # --- (a) raw yields per fill vs beta --------------------------------
    for f, (col, lab) in enumerate(((C_TRUTH, r"$P_{zz}=+%.1f$ fill" % args.pzz),
                                    (C_FIT, r"$P_{zz}=-%.1f$ fill" % (2 * args.pzz)))):
        yb = f0["counts"][f].sum(axis=0)
        ax1.errorbar(bc, yb / yb.mean(), yerr=np.sqrt(yb) / yb.mean(), fmt="o-",
                     color=col, ms=3.5, lw=1, capsize=2, label=lab)
    fake = np.mean(np.cos(2.0 * cresp.beta_reco))
    ax1.annotate(r"slot cutout $|p_x|<%.2f$, $|p_y|<%.2f$ GeV: "
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
    ax1.set_title(r"(a) raw spin-sorted yields, $|t|\in[%.2f,%.2f]$: the "
                  r"cutout dominates" % (t_edges[0], t_edges[1]), fontsize=9)
    ax1.legend(fontsize=7, loc="upper right")
    ax1.tick_params(labelsize=8)

    # --- (b) acceptance-free modulation, projections -------------------
    r, var, sig2, pbar = reco.spin_state_ratio(
        f0["counts"].reshape(2, -1), [0.5, 0.5], pzz_list)
    aa, bb = np.meshgrid(ac, bc, indexing="ij")
    bm = f0["beta_means"]
    u = reco.unpolarized_modulation_2d(ae, be, args.u1, args.u2, beta_means=bm)
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
    r0, v0, s0, p0 = reco.spin_state_ratio(mu.reshape(2, -1), [0.5, 0.5],
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
        r"(6R), %s, $P_{zz}=%.2f$: $\sigma_\theta=%.0f\,\mu$rad "
        r"($10\sigma_\theta A p_u=%.2f$ GeV), slot-like %s cutout "
        r"$|p_x|<%.2f$, $|p_y|<%.2f$ GeV"
        "\n" r"$N_{\rm tag}$ = %.2g (1 yr) vs %.2g with the constant 0.20 GeV cut; "
        r"deformation $c_2=2a_2$ (Eq. 9 convention), $a_e=%.3f$"
        "\n" r"$u_1=%.2f$, $u_2=%.2f$ (ZEUS LPS $1\sigma$ bounds); "
        r"two-fill ratio fit with the MC template basis; statistical errors only"
        % (config.label(), args.pzz, 1e6 * args.sigma_theta, cresp.pt_cut,
           args.shape, cresp.cut_pt_xy[0], cresp.cut_pt_xy[1], n_tag,
           tagged[HIGH_ACCEPTANCE.name].sum(), args.amp, args.u1, args.u2),
        fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.905))
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "money_cos2phi_coherent_reco_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print("coherent produced (1 yr): %.3g; tagged: %s cutout |px|<%.3f, "
          "|py|<%.3f GeV (divergence aspect %.2f) -> acc %.4f, N_tag %.3g; "
          "constant-cut reference %.3g"
          % (n_produced, args.shape, cresp.cut_pt_xy[0], cresp.cut_pt_xy[1],
             args.aspect, cresp.acceptance, n_tag,
             tagged[HIGH_ACCEPTANCE.name].sum()))
    print("cutout fake <cos 2beta> = %.3f" % fake)
    for (tlo, thi), f1, f10 in zip(zip(t_edges[:-1], t_edges[1:]), fits1, fits10):
        tr = f1["truth"]
        print("t in [%.2f,%.2f]: N=%.3g  a_t truth %.4f  fit %.4f +- %.4f (1yr) "
              "%.4f +- %.4f (10yr) | a_e truth %.4f fit %.4f +- %.4f (1yr) "
              "%.4f +- %.4f (10yr) | a_m %.4f +- %.4f"
              % (tlo, thi, f1["n"], tr["a_t"], f1["a_t"], f1["err_t"],
                 f10["a_t"], f10["err_t"], tr["a_e"], f1["a_e"], f1["err_e"],
                 f10["a_e"], f10["err_e"], f1["a_m"], f1["err_m"]))
        if not args.no_sin:
            print("    null test: a_e_s %+.4f +- %.4f, a_t_s %+.4f +- %.4f, "
                  "a_m_s %+.4f +- %.4f  -> sin/cos %+.4f (alpha) %+.4f "
                  "(beta); axis resolution %.0f / %.0f mrad (1 yr)"
                  % (f1["a_e_s"], f1["err_e_s"], f1["a_t_s"], f1["err_t_s"],
                     f1["a_m_s"], f1["err_m_s"],
                     f1["a_e_s"] / f1["a_e"], f1["a_t_s"] / f1["a_t"],
                     5e2 * f1["err_e_s"] / abs(f1["a_e"]),
                     5e2 * f1["err_t_s"] / abs(f1["a_t"])))


if __name__ == "__main__":
    main()
