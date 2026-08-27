#!/usr/bin/env python3
"""Pricing a lithium tagging optics at IP6 (Report 1, Section 6.1).

The coherent intact-6Li recoil keeps the beam rigidity, so the only handle
is its angle at the IP against the 10 sigma beam envelope and the Roman-Pot
aperture, and the tagged fraction is exp(-B (10 sigma_theta A p_u)^2).
With the Yellow Report divergences no published optics tags it (plans/10).
Two levers exist, and they are different kinds of change:

  * beta* at IP6 -- a machine lattice setting, exactly what distinguishes
    the Yellow Report's high-acceptance from its high-divergence optics.
    sigma_theta ~ 1/sqrt(beta*), so a factor r in beta* buys sqrt(r) in
    the envelope; but luminosity ~ 1/beta*, so the same r costs a factor
    r in rate.  The figure of merit is L x acceptance, and this script
    maximises it.
  * the far-forward side -- the pot insertion (and, at a secondary focus,
    the transport).  The pots must FOLLOW the envelope: with the silicon
    fixed at the aperture measured in the ePIC geometry
    (reco.RP_APERTURE_MEASURED) the acceptance saturates and every further
    beta* only loses luminosity.

Per 6Li configuration the script scans r = beta*/beta*_HA, with the
envelope treated as the planar pots see it -- a rectangle of half-widths
10 sigma_h x 10 sigma_v in angle (reco.rp_measure's model, exact through
reco.rp_hole_acceptance) -- and reports the optimum of eps/r for pots that
follow and for pots fixed at the current aperture, the luminosity fraction
there, the tagged events per year at the programme's placeholder of
10 fb^-1/u per year in the high-acceptance optics, and the best-super-bin
5 sigma floor on the cos 2phi amplitude per unit P_zz.  The IR-8 secondary
focus (~20% for 6Li, interpolated from Chang et al. PRD 113:032018 at the
second IR's own luminosity) is the horizontal reference.

Because a recoil escapes through the horizontal gap OR the vertical one, only
the plane it escapes through needs the tagging beta*: the default scan
de-squeezes the HORIZONTAL plane alone and holds the vertical at its
high-acceptance value, so L/L_HA = 1/sqrt(r_h) rather than 1/r -- the
both-planes case is drawn for comparison.  Two assumptions carry the result:
(i) L ~ 1/beta* per plane, which the Yellow Report's own HD/HA pairs
reproduce to 5% (10.0/3.14 at 275 GeV for a divergence ratio 119/65), with
the electron beta* raised in step to keep the spots matched; (ii) a
parallel-to-point far-forward transport for the de-squeezed lattice
(R11 sigma* << R12 sigma_theta at the pot planes, as TOTEM's high-beta*
optics), so that the 10 sigma envelope at the pots is 10 R12 sigma_theta --
the present lattice satisfies this at the Yellow Report optics, and whether a
de-squeezed one does is the first question for C-AD.  Dispersion at the pots
is added in quadrature to the horizontal envelope as 10 D dp/p / R12 with
D = 0.3 m (read off the alpha/triton positions at 18 x 275) and the Yellow
Report dp/p.  Absolute lithium luminosities are not published; everything
here is relative to the high-acceptance optics.

Usage:  python3 scripts/tagging_optics.py [--outdir .] [--slope-b 50] [--pzz 0.6]
"""

import argparse
import pathlib
import sys

import numpy as np

_SCRIPTS = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS.parent))
sys.path.insert(0, str(_SCRIPTS.parent.parent / "fastsim"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from polligen import coherent as coh  # noqa: E402
from polligen import reco  # noqa: E402
from polli_fastsim import beams, fom  # noqa: E402
sys.path.insert(0, str(_SCRIPTS))
from money_cos2phi_coherent import best_superbin  # noqa: E402  (the paper's own bin)

C = ("#0072B2", "#D55E00", "#009E73")
C_GREY = "0.45"
KEYS = ("5x41", "10x100", "18x275")
# eSTARlight IR-8 secondary-focus intact-recoil efficiency x acceptance:
# d 47%, 3He 32%, 4He 29%, 7Li 17.8% (Chang et al. PRD 113:032018); no 6Li
# entry exists and 0.20 is this programme's interpolation
# (coherent_optics_scan.py).
IR8_LI6_INTERPOLATED = 0.20
R12 = 30.6          # m, IP angle -> pot-plane x, measured at 18 x 275 (tools/fullsim)
D_POT = 0.30        # m, dispersion at the pots from the alpha/triton positions (tools/fullsim)


def rect_acceptance(slope_b, cx_gev, cy_gev):
    """Tagged fraction of exp(-B t) outside the rectangle |px| < cx, |py| < cy."""
    return float(reco.rp_hole_acceptance(slope_b, cx_gev, cy_gev,
                                         shape="rectangle")["acc"])


def best_superbin_fraction(n_map, proj):
    """Fraction of the tagged sample in the best super-bin as the paper's
    Section 6.2 defines it -- money_cos2phi_coherent.best_superbin: the
    maximal cell padded by one on each side and clipped at the grid edge
    (3 x 2 at the Q2 = 1 GeV2 edge), the bin that holds 1.75e6 of the
    1.67e7 tagged events of the 0.20 GeV reference at the mid configuration."""
    a = np.asarray(n_map, dtype=float)
    sel, *_edges = best_superbin(proj, a)
    return float(a[sel].sum() / a.sum())


def mean_t_tagged(slope_b, cx_gev, cy_gev):
    """<|t|> of the recoils outside the rectangle, for the deformation term
    of Section 6.3, a_2 = -(P_zz/4) eps_B0 B |t|: with dN ~ exp(-B p^2) p dp dphi
    accepted above rho(phi), <p^2> = sum_phi e^{-B rho^2}(rho^2 + 1/B) /
    sum_phi e^{-B rho^2}."""
    h = reco.rp_hole_acceptance(slope_b, cx_gev, cy_gev, shape="rectangle")
    eps = np.asarray(h["eps"], dtype=float)
    rho2 = -np.log(np.clip(eps, 1e-300, None)) / slope_b
    return float((eps * (rho2 + 1.0 / slope_b)).sum() / eps.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slope-b", type=float, default=50.0)
    ap.add_argument("--pzz", type=float, default=0.60)
    ap.add_argument("--lumi-1yr", type=float, default=10.0,
                    help="placeholder fb^-1/u per year in the high-acceptance optics")
    ap.add_argument("--gluonic-amp", type=float, default=0.01,
                    help="flat exotic-glue amplitude per unit P_zz to time (5 sigma)")
    ap.add_argument("--r-max", type=float, default=2000.0)
    ap.add_argument("--no-dispersion", action="store_true")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sc = coh.CoherentScenario(slope_b=args.slope_b)
    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    r = np.logspace(np.log10(0.25), np.log10(args.r_max), 400)
    lines = ["# tagging optics at IP6: r = beta*/beta*_HA, luminosity ~ 1/r, "
             "envelope = rectangle 10 sigma_h x 10 sigma_v (planar pots), "
             "B = %g GeV^-2, P_zz = %.2f, %g fb^-1/u per year at r = 1"
             % (args.slope_b, args.pzz, args.lumi_1yr)]
    results = []
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.6))

    for i, (cfg, key, col) in enumerate(zip(beams.default_configs("6Li"), KEYS, C)):
        ax = axes[i]
        pu = cfg.ion_momentum_per_nucleon
        p_ion = cfg.ion.A * pu
        sh, sv = reco.sigma_theta_for(cfg, "high-acceptance")
        sh_hd, sv_hd = reco.sigma_theta_for(cfg, "high-divergence")
        apx, apy = reco.rp_aperture_for(pu)
        _proj, n_coh, _ = coh.project_coherent(cfg, scenario, sc)
        n_produced = float(n_coh.sum())
        f_bin = best_superbin_fraction(n_coh, _proj)
        n_ref = n_produced * sc.tag_acceptance(0.20)      # the 0.20 GeV reference tag

        from polli_fastsim import farforward as _ff
        dpp = 1e-4 * _ff.YR_PROTON_DIVERGENCE[key][1]
        disp = 0.0 if args.no_dispersion else D_POT * dpp / R12      # rad, in quadrature

        def env_h(rr):                       # horizontal 10 sigma envelope at r_h = rr
            return 10 * ((sh / rr ** 0.5) ** 2 + disp ** 2) ** 0.5

        # one-plane de-squeeze: horizontal beta* x r, vertical held at HA
        eps_env = np.array([rect_acceptance(args.slope_b, env_h(rr) * p_ion, 10 * sv * p_ion)
                            for rr in r])
        eps_fix = np.array([rect_acceptance(args.slope_b, max(env_h(rr), apx) * p_ion,
                                            max(10 * sv, apy) * p_ion) for rr in r])
        lum = 1.0 / r ** 0.5
        prod_env, prod_fix = eps_env * lum, eps_fix * lum
        # both planes de-squeezed by r (the naive case), for comparison
        eps_both = np.array([rect_acceptance(args.slope_b, env_h(rr) * p_ion,
                                             10 * sv / rr ** 0.5 * p_ion) for rr in r])
        prod_both = eps_both / r

        def at(idx, eps, prod):
            n_tag = n_produced * prod[idx]                     # per year
            floor = 5.0 * (2.0 / max(n_tag * f_bin, 1e-300)) ** 0.5 / args.pzz
            years = (floor / args.gluonic_amp) ** 2
            env_hh, env_v = env_h(r[idx]), 10 * sv
            t_mean = mean_t_tagged(args.slope_b, env_hh * p_ion, env_v * p_ion)
            a2 = float(sc.a2_deformation(t_mean, 1.0))      # per unit P_zz, this window
            return dict(r=float(r[idx]), eps=float(eps[idx]), lum=float(lum[idx]),
                        n_tag=n_tag, floor=floor, years_to_5sig=years,
                        env_h=env_hh, env_v=env_v, t_mean=t_mean, a2=abs(a2),
                        shape_sigma=5.0 * abs(a2) / floor, ref_ratio=n_ref / max(n_tag, 1e-300))

        k_env = int(np.argmax(prod_env))
        k_fix = int(np.argmax(prod_fix))
        k_one = int(np.argmin(np.abs(r - 1.0)))
        k_both = int(np.argmax(prod_both))
        opt_env, opt_fix, base = at(k_env, eps_env, prod_env), at(k_fix, eps_fix, prod_fix), at(k_one, eps_fix, prod_fix)
        opt_both = dict(r=float(r[k_both]), eps=float(eps_both[k_both]), prod=float(prod_both[k_both]),
                        n_tag=n_produced * prod_both[k_both])
        r_circ = (sh / reco.sigma_theta_tagging(cfg, slope_b=args.slope_b)) ** 2
        r_hd = (sh / sh_hd) ** 2
        results.append((cfg.label(), key, opt_env, opt_fix, base, n_produced, f_bin, r_circ))

        lines += [
            "%-24s  N_coh(%g fb^-1/u) = %.3g, best super-bin fraction %.3f, "
            "sigma_HA h/v = %.0f/%.0f urad, aperture %.2f/%.2f mrad, dispersive envelope term %.0f urad"
            % (cfg.label(), args.lumi_1yr, n_produced, f_bin, 1e6 * sh, 1e6 * sv,
               1e3 * apx, 1e3 * apy, 1e6 * 10 * disp),
            "   both planes de-squeezed by r (L = 1/r): optimum r = %.1f, eps = %.3f, eps*L = %.4f, N_tag/yr = %.2e"
            % (opt_both["r"], opt_both["eps"], opt_both["prod"], opt_both["n_tag"]),
            "   r = 1 (HA, pots at the measured aperture): eps = %.2e, N_tag/yr = %.2e"
            % (base["eps"], base["n_tag"]),
            "   HORIZONTAL de-squeeze only, pots FOLLOW: optimum r_h = %.1f, eps = %.3f, L/L_HA = 1/%.1f, "
            "eps*L = %.4f, envelope 10sigma h/v = %.2f/%.2f mrad (%.2f/%.2f GeV, vertical at HA), N_tag/yr = %.2e "
            "(0.20 GeV reference at 10 fb^-1/u: %.2e, i.e. %.0fx more), best-bin 5sigma floor/yr = %.2e "
            "(per unit P_zz), <|t|> tagged = %.4f GeV^2 -> shape term a2 = %.3f per unit P_zz = %.1f sigma/yr, "
            "years to 5sigma on %.3f = %.1f, IR-8 at L_HA: 0.20 / (eps*L) = %.0fx"
            % (opt_env["r"], opt_env["eps"], 1 / opt_env["lum"], opt_env["eps"] * opt_env["lum"],
               1e3 * opt_env["env_h"], 1e3 * opt_env["env_v"], opt_env["env_h"] * p_ion,
               opt_env["env_v"] * p_ion, opt_env["n_tag"], n_ref, opt_env["ref_ratio"],
               opt_env["floor"], opt_env["t_mean"], opt_env["a2"], opt_env["shape_sigma"],
               args.gluonic_amp, opt_env["years_to_5sig"],
               IR8_LI6_INTERPOLATED / (opt_env["eps"] * opt_env["lum"])),
            "   pots FIXED at the current aperture: optimum r_h = %.1f, eps = %.2e, L/L_HA = 1/%.1f, "
            "N_tag/yr = %.2e, floor/yr = %.2e"
            % (opt_fix["r"], opt_fix["eps"], 1 / opt_fix["lum"], opt_fix["n_tag"], opt_fix["floor"]),
            "   plans/10 circular-isotropic optimum for reference: r = %.1f;  HD optics at r = %.2f"
            % (r_circ, r_hd),
        ]

        ax.plot(r, np.clip(eps_env, 1e-30, None), "--", color=col, lw=1.3,
                label="tagged fraction ε, horizontal β* × r, pots follow the 10σ envelope")
        ax.plot(r, np.clip(prod_env, 1e-30, None), "-", color=col, lw=2.2,
                label="ε × L/L$_{HA}$ (the yield), horizontal plane only, L = 1/√r")
        ax.plot(r, np.clip(prod_both, 1e-30, None), "-", color=col, lw=0.9, alpha=0.6,
                label="ε × L/L$_{HA}$ with both planes × r, L = 1/r")
        ax.plot(r, lum, "-", color=C_GREY, lw=0.8, label="L/L$_{HA}$ = 1/√r")
        ax.plot([opt_env["r"]], [opt_env["eps"] * opt_env["lum"]], "o", color=col, ms=7,
                mfc="white", mew=1.8, label="optimum of the yield")
        # pots fixed at the measured aperture: below the axis at every r
        # (max %.1e), so it is stated rather than drawn
        ax.axhline(IR8_LI6_INTERPOLATED, color="0.25", lw=1.0, ls="-.")
        ax.annotate("IR-8 secondary focus: ≈ 20% (our interpolation) if L = L$_{HA}$",
                    xy=(0.3, IR8_LI6_INTERPOLATED), xytext=(0, 4), textcoords="offset points",
                    fontsize=7, color="0.25", ha="left", va="bottom")
        ax.axvline(1.0, color=C_GREY, lw=0.8, ls=":")
        ax.axvline(r_hd, color=C_GREY, lw=0.8, ls=":")
        if abs(r_hd - 1.0) < 0.3:
            ax.text(1.0, 1.25, "HA ≈ HD" if r_hd != 1.0 else "HA = HD", fontsize=7,
                    color=C_GREY, ha="center")
        else:
            ax.text(1.0, 1.25, "HA", fontsize=7, color=C_GREY, ha="center")
            ax.text(r_hd, 1.25, "HD", fontsize=7, color=C_GREY, ha="center")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(r[0], r[-1]); ax.set_ylim(1e-6, 2.0)
        ax.set_xlabel(r"$r = \beta^*_x/\beta^*_{x,\rm HA}$")
        ax.set_title("%s  (σ$_{HA}$ = %.0f/%.0f µrad, pots at %.2f/%.2f mrad)"
                     % (cfg.label(), 1e6 * sh, 1e6 * sv, 1e3 * apx, 1e3 * apy), fontsize=8.8)
        ax.text(0.03, 0.03,
                "optimum, horizontal β* × %.0f: ε = %.2f, L/L$_{HA}$ = 1/%.0f\n"
                "N$_{tag}$ = %.1f×10$^{%d}$ / yr, 5σ floor %.1f%% / yr (shape term %.1f%%: %.0fσ / yr)\n"
                "pots fixed at the measured aperture: ε < %.0e at every r"
                % (opt_env["r"], opt_env["eps"], 1 / opt_env["lum"],
                   opt_env["n_tag"] / 10 ** int(np.floor(np.log10(opt_env["n_tag"]))),
                   int(np.floor(np.log10(opt_env["n_tag"]))), 100 * opt_env["floor"],
                   100 * opt_env["a2"], opt_env["shape_sigma"],
                   max(float(eps_fix.max()), 1e-30)),
                transform=ax.transAxes, fontsize=7.2, va="bottom", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=col, lw=0.8))
        ax.grid(alpha=0.25, lw=0.5)
        if i == 0:
            ax.set_ylabel("fraction")
            handles, labels = ax.get_legend_handles_labels()

    fig.suptitle("A lithium tagging optics at IP6, priced: horizontal β* de-squeezed, vertical at high acceptance; "
                 "σ$_θ$ ∝ 1/√β*, L ∝ 1/β* per plane, parallel-to-point transport assumed\n"
                 "B = %g GeV$^{-2}$, P$_{zz}$ = %.2f; yields per year at %g fb$^{-1}$/u in the high-acceptance optics"
                 % (args.slope_b, args.pzz, args.lumi_1yr), fontsize=9.2)
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=7.0, frameon=False,
               bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.09, 1, 0.92))
    out = outdir / "tagging_optics_6Li.png"
    fig.savefig(out, dpi=140)
    lines.append("wrote %s" % out)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
