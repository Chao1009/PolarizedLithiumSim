#!/usr/bin/env python3
"""Money plot 6R: the coherent intact-6Li channel at the RECONSTRUCTED
level (plans/07 WP3 -> WP5) -- money plot 6 re-derived with

  * the ANGULAR near-beam cut, pT_cut = 10 sigma_theta A p_u (proton-
    derived sigma_theta = 73 microrad, high-acceptance optics), and a
    rectangular Roman-Pot cutout of aspect ratio r (default 1.25) whose
    sides are parallel to the vertical spin axis;
  * the Roman-Pot emulation: divergence smearing of the recoil angle,
    |t| = pT^2 reconstructed with x_L = 1, reco t bins;
  * the TWO azimuths: alpha = phi_e - phi_S (electron) and beta =
    phi_t - phi_S (recoil); the deformation term modulates cos 2beta and
    the gluon-transversity term cos 2alpha; the unpolarized lepton-plane
    /recoil-plane harmonics u_1 cos(alpha-beta), u_2 cos 2(alpha-beta)
    are placeholders (assumption);
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
                    help="unpolarized cos(alpha-beta) [placeholder]")
    ap.add_argument("--u2", type=float, default=0.02,
                    help="unpolarized cos 2(alpha-beta) [placeholder]")
    ap.add_argument("--sigma-theta", type=float, default=reco.SIGMA_THETA_HA,
                    help="beam angular divergence [rad] (proton-derived)")
    ap.add_argument("--aspect", type=float, default=1.25,
                    help="cutout aspect ratio c_y/c_x [assumption]")
    ap.add_argument("--shape", default="rectangle",
                    choices=("rectangle", "ellipse"))
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
                                n_mc=args.n_mc, rng=rng)
    n_tag = n_produced * cresp.acceptance
    plan = bk.tensor_flip_plan(args.pzz)
    pzz_list = [args.pzz, -2.0 * args.pzz]

    def a_t_func(t):   # cos 2beta coefficient per unit P_zz (Eq. 9 convention)
        return sc.cos2phi_coefficient_deformation(t, 1.0)

    t_edges = [0.05, 0.08, 0.12, 0.17, 0.25]
    fits1, fits10 = [], []
    for tlo, thi in zip(t_edges[:-1], t_edges[1:]):
        fits1.append(rp.measure_coherent(
            cresp, n_produced, plan, tlo, thi, args.amp, a_t_func, a_m=args.a_m,
            u1=args.u1, u2=args.u2, rng=rng))
        fits10.append(rp.measure_coherent(
            cresp, n_produced * lumi_ratio, plan, tlo, thi, args.amp, a_t_func,
            a_m=args.a_m, u1=args.u1, u2=args.u2, rng=rng))

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
    ax1.annotate(r"cutout $r=%.2f$: $\langle\cos2\beta\rangle=%.2f$ of the "
                 r"tagged sample" "\n" r"(a single-fill fit would report "
                 r"$a_t\approx%.2f$ vs truth %.3f)"
                 % (args.aspect, fake, fake / args.pzz, f0["truth"]["a_t"]),
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
    c2_mean = float(np.mean(bm["c2"]))   # a_t offset of the alpha projection
    tb = t2d.mean(axis=0) - f0["const"]
    eb = np.sqrt(var2d.sum(axis=0)) / aa.shape[0]
    ta = t2d.mean(axis=1) - f0["const"] - f0["a_t"] * c2_mean
    ea = np.sqrt(var2d.sum(axis=1)) / aa.shape[1]
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
        r"(6R), %s, $P_{zz}=%.2f$: angular cut $10\sigma_\theta A p_u$ "
        r"($\sigma_\theta=%.0f\,\mu$rad $\to p_T>%.2f$ GeV), %s cutout $r=%.2f$"
        "\n" r"$N_{\rm tag}$ = %.2g (1 yr) vs %.2g with the constant 0.20 GeV cut; "
        r"deformation $c_2=2a_2$ (Eq. 9 convention), $a_e=%.3f$, "
        r"$u_1=%.2f$, $u_2=%.2f$ placeholders; two-fill ratio fit; stat. only"
        % (config.label(), args.pzz, 1e6 * args.sigma_theta, cresp.pt_cut,
           args.shape, args.aspect, n_tag, tagged[HIGH_ACCEPTANCE.name].sum(),
           args.amp, args.u1, args.u2), fontsize=9.3)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "money_cos2phi_coherent_reco_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print("coherent produced (1 yr): %.3g; tagged: angular %s cutout r=%.2f "
          "-> acc %.4f, N_tag %.3g; constant-cut reference %.3g"
          % (n_produced, args.shape, args.aspect, cresp.acceptance, n_tag,
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


if __name__ == "__main__":
    main()
