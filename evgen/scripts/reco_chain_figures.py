#!/usr/bin/env python3
"""Figures for the reconstruction-chain report (reports/
reconstruction_chain_report): what the inclusive and coherent cos 2phi
measurements actually record, and what the reconstruction must do with
it -- quantified with the polligen.reco module on the money-plot beam
configuration.

reco_chain_inclusive_6Li.png
  (a) the (x, Q2) plane with the analysis acceptance, iso-y lines, the
      four sweet-spot super-bins of money plot 5, and the region where
      the scattered electron alone cannot measure y (hadronic method
      required);
  (b) electron-method resolution dy/y vs y for calorimeter, tracking and
      beam-energy-spread inputs against the hadronic-method band, with
      the sweet-spot y values of the three beam configurations;
  (c) bin migration with the mixed (e-Sigma) method: super-bin purity
      and the amplitude seen in the reconstructed bin relative to the
      true bin, vs the assumed hadronic-y resolution;
  (d) the single-fill fit vs the spin-state-sorted ratio estimator under
      a smooth phi-dependent efficiency aligned with the spin axis.

reco_chain_coherent_6Li.png
  (a) the recoil transverse-momentum plane with the Roman-Pot cutout;
  (b) fake a_2 of the tagged sample vs cutout aspect ratio against the
      anchored deformation amplitude and the statistical floors;
  (c) tag acceptance vs beam momentum per nucleon for an ANGULAR
      envelope cut (the pT cut scales with A p_u);
  (d) Roman-Pot emulation: reconstructed vs true |t| and the phi_t
      resolution from the beam divergence.

Usage:  python3 scripts/reco_chain_figures.py
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
from matplotlib.patches import Ellipse, Rectangle  # noqa: E402

from money_cos2phi import (build_delta_model, pick_sweet_spots_banded,  # noqa: E402
                           superbin_edges)

from polligen import bookkeeping as bk  # noqa: E402
from polligen import reco  # noqa: E402
from polligen.estimators import cos2phi_fit_binned  # noqa: E402
from polligen.sample import InclusiveSampler  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.asymmetries import a_cos2phi  # noqa: E402
from polli_fastsim.kinematics import scattered_electron, y_from_xq2  # noqa: E402
from polli_fastsim.polarized import toy_b1  # noqa: E402

C_BLUE, C_VERM, C_GREEN, C_GREY = "#0072B2", "#D55E00", "#009E73", "0.45"
C_PURPLE, C_ORANGE = "#CC79A7", "#E69F00"
UNPOL = bk.SpinCategory("unpol", 1.0, (1 / 3, 1 / 3, 1 / 3))


def sweet_spots(config_index, pzz=0.6):
    config = beams.default_configs("6Li")[config_index]
    scenario = fom.Scenario(lumi_fb_per_nucleon=10.0, pol_ion_tensor=pzz)
    args = argparse.Namespace(delta_model="moment_A", variant="mid_x",
                              dilution=1.0 / 3.0, scale=1e-2)
    model, _ = build_delta_model(args, config, scenario)
    proj = fom.project_rates(config, scenario)
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1, delta_func=model)
    obs = fom.project_observables(config, scenario, proj, kern.g1_model,
                                  toy_b1, model)
    spots = pick_sweet_spots_banded(proj, obs["sig_a_cos2phi"])[:4]
    return config, scenario, proj, obs, spots, kern


def migration_study(config, scenario, proj, spots, kern, rel_res_list,
                    n_events=1500000, seed=20260824):
    """Purity and amplitude ratio per super-bin for the mixed method
    (Q2 from e' with EMCal resolution, y hadronic with rel_res) and for
    the electron method alone (EMCal resolution)."""
    s = config.sqrt_s_per_nucleon ** 2
    e_e = config.electron_energy
    rng = np.random.default_rng(seed)
    sampler = InclusiveSampler(kern, config, scenario, nx=60, nq2=45)
    ev = sampler.sample_category(UNPOL, n=n_events, rng=rng)
    x, q2, y = ev["x"], ev["q2"], ev["y"]
    t = kern.tables(x, q2)
    amp = a_cos2phi(t["delta"], t["f1"], t["f2"], x, y)
    e_p, th, _eta = scattered_electron(x, y, s, e_e)
    de = reco.emcal_resolution(e_p)
    e_s, th_s, _ = reco.smear_electron(e_p, th, np.zeros_like(e_p), de,
                                       1e-3, 1e-3, rng)
    q2_e, y_e, x_e = reco.electron_method(e_s, th_s, e_e, s)
    results = []
    for rr in rel_res_list:
        if rr is None:   # electron method alone
            xr, q2r = x_e, q2_e
        else:
            xr = reco.mixed_method(q2_e, reco.hadronic_y(y, rr, rng), s)
            q2r = q2_e
        row = []
        for xs, qs, i, j in spots:
            xlo, xhi, qlo, qhi = superbin_edges(proj, i, j)
            true_in = (x >= xlo) & (x < xhi) & (q2 >= qlo) & (q2 < qhi)
            reco_in = (np.isfinite(xr) & (xr >= xlo) & (xr < xhi)
                       & (q2r >= qlo) & (q2r < qhi))
            purity = (true_in & reco_in).sum() / max(reco_in.sum(), 1)
            eff = (true_in & reco_in).sum() / max(true_in.sum(), 1)
            ratio = (amp[reco_in].mean() / amp[true_in].mean()
                     if reco_in.any() else np.nan)
            row.append((purity, eff, ratio, reco_in.sum() / max(true_in.sum(), 1)))
        results.append(row)
    return results


def inclusive_figure(outdir):
    config, scenario, proj, obs, spots, kern = sweet_spots(1)
    s = config.sqrt_s_per_nucleon ** 2
    e_e = config.electron_energy
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12.4, 8.8))

    # --- (a) phase space, iso-y, sweet spots, method zones ---------------
    xg = np.logspace(-4, 0, 400)
    de_cal = float(reco.emcal_resolution(e_e * 0.99))     # ~1.2% at 10 GeV
    y_30 = de_cal / (0.30 + de_cal)    # electron method dy/y = 30 %
    y_50 = de_cal / (0.50 + de_cal)
    acc = proj.accepted
    ax1.pcolormesh(proj.x_edges, proj.q2_edges,
                   np.where(acc, 1.0, np.nan).T, cmap="Greys", vmin=0,
                   vmax=6, shading="flat", zorder=0)
    ax1.fill_between(xg, scenario.y_min * s * xg, y_30 * s * xg,
                     color=C_VERM, alpha=0.18, lw=0, zorder=1,
                     label=r"$e'$ alone: $\delta y/y>30\%%$ "
                           r"($\delta E'/E'=%.1f\%%$)" % (100 * de_cal))
    for yv, ls, lab in ((scenario.y_min, "-", r"$y=0.01$ (cut)"),
                        (y_50, ":", r"$\delta y/y=50\%$"),
                        (y_30, "--", r"$\delta y/y=30\%$"),
                        (0.1, "-.", r"$y=0.1$"),
                        (scenario.y_max, "-", r"$y=0.95$ (cut)")):
        ax1.plot(xg, yv * s * xg, ls, color="0.25", lw=0.9, zorder=2)
        # label in the right margin where the line leaves the axes
        ax1.annotate(lab.replace(" (cut)", ""), xy=(1.06, yv * s),
                     fontsize=6.3, color="0.2", ha="left", va="center",
                     annotation_clip=False, zorder=6)
    for n, (xs, qs, i, j) in enumerate(spots):
        xlo, xhi, qlo, qhi = superbin_edges(proj, i, j)
        ax1.add_patch(Rectangle((xlo, qlo), xhi - xlo, qhi - qlo, fill=False,
                                ec=C_BLUE, lw=1.6, zorder=4))
        yv = y_from_xq2(xs, qs, s)
        e_p, th, eta = scattered_electron(xs, yv, s, e_e)
        left = (n == 1)
        ax1.annotate("%d: $y$=%.3f, $\\theta'$=%.0f mrad, $\\eta$=%.1f"
                     % (n + 1, yv, 1e3 * (np.pi - th), eta),
                     xy=(xlo / 1.12 if left else xhi * 1.12, qs),
                     fontsize=6.5, color=C_BLUE, va="center",
                     ha="right" if left else "left", zorder=5,
                     bbox=dict(boxstyle="round,pad=0.1", fc="white",
                               ec="none", alpha=0.8))
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(1e-4, 1.0)
    ax1.set_ylim(1.0, 2e3)
    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$Q^2$ [GeV$^2$]")
    ax1.set_title("(a) where the sweet spots sit: %s" % config.label(),
                  fontsize=9)
    ax1.legend(fontsize=7, loc="upper left")
    ax1.tick_params(labelsize=8)

    # --- (b) dy/y vs y ------------------------------------------------------
    yy = np.logspace(-2.3, -0.05, 300)
    th_ref = np.pi - scattered_electron(0.03, yy, s, e_e)[1]
    curves = ((0.005, 1e-3, C_GREEN,
               r"$\delta E'/E'=0.5\%$, $\delta\theta'=1$ mrad"),
              (de_cal, 1e-3, C_BLUE,
               r"EMCal $2\%%/\sqrt{E}\oplus1\%%$ ($%.1f\%%$ at "
               r"$E'\approx E_e$)" % (100 * de_cal)),
              (0.03, 3e-3, C_VERM,
               r"tracking-only $\sigma_p/p=3\%$, 3 mrad ($\eta\approx-3$)"))
    for de, dth, col, lab in curves:
        _, dy, _ = reco.electron_method_resolution(yy, th_ref, de, dth)
        ax2.plot(yy, dy, "-", color=col, lw=1.6, label=lab)
    ax2.plot(yy, (1 - yy) / yy * 1e-3, "--", color="0.3", lw=1.0,
             label=r"beam energy spread $10^{-3}$ alone")
    ax2.axhspan(0.15, 0.30, color=C_ORANGE, alpha=0.25, lw=0,
                label=r"hadronic $y$ (JB / $\Sigma$), 15–30% assumed")
    for ci, mk in ((0, "s"), (1, "o"), (2, "^")):
        cfg_i, _, _, _, spots_i, _ = sweet_spots(ci)
        s_i = cfg_i.sqrt_s_per_nucleon ** 2
        ys = np.array([y_from_xq2(xs, qs, s_i) for xs, qs, _, _ in spots_i])
        ths = np.array([np.pi - scattered_electron(xs, y_from_xq2(xs, qs, s_i),
                                                   s_i, cfg_i.electron_energy)[1]
                        for xs, qs, _, _ in spots_i])
        _, dys, _ = reco.electron_method_resolution(ys, ths, de_cal, 1e-3)
        ax2.plot(ys, dys, mk, color="black", ms=6,
                 mfc="black" if ci == 1 else "white",
                 label="sweet spots, %s (EMCal curve)" % cfg_i.label())
    ax2.axhline(1.0, color="0.85", lw=0.6)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(5e-3, 1.0)
    ax2.set_ylim(3e-3, 30)
    ax2.set_xlabel(r"$y$")
    ax2.set_ylabel(r"$\delta y/y$ (electron method)")
    ax2.set_title(r"(b) $y$ from $e'$ alone: $\delta y/y=\frac{1-y}{y}"
                  r"\,[\delta E'/E'\oplus\tan\frac{\theta'}{2}\delta\theta'"
                  r"\oplus\delta E_e/E_e]$", fontsize=9)
    ax2.legend(fontsize=6.3, loc="lower left")
    ax2.tick_params(labelsize=8)

    # --- (c) bin migration with the mixed method ---------------------------
    res_list = (0.15, 0.20, 0.30, None)
    res = migration_study(config, scenario, proj, spots, kern, res_list)
    xpos = np.arange(len(spots))
    wbar = 0.19
    cols = (C_GREEN, C_BLUE, C_VERM, "0.5")
    labs = (r"mixed, $\delta y_{\rm had}/y=15\%$", r"mixed, $20\%$",
            r"mixed, $30\%$", r"electron method alone (EMCal)")
    for k, (rr, col, lab) in enumerate(zip(res_list, cols, labs)):
        pur = [r[0] for r in res[k]]
        rat = [r[2] for r in res[k]]
        ax3.bar(xpos + (k - 1.5) * wbar, pur, wbar, color=col, alpha=0.75,
                label=lab)
        for xp_, rt, pu in zip(xpos + (k - 1.5) * wbar, rat, pur):
            if np.isfinite(rt):
                ax3.annotate("%.2f" % rt, xy=(xp_, pu + 0.015), fontsize=6.2,
                             ha="center", rotation=90, va="bottom")
    ax3.set_xticks(xpos)
    ax3.set_xticklabels(["spot %d\n$y$=%.3f" % (n + 1, y_from_xq2(xs, qs, s))
                         for n, (xs, qs, _, _) in enumerate(spots)],
                        fontsize=8)
    ax3.set_ylim(0, 1.05)
    ax3.set_ylabel("super-bin purity (bars)")
    ax3.set_title(r"(c) $(x,Q^2)$ migration: purity of each super-bin; "
                  r"numbers = $\langle A\rangle_{\rm reco\,bin}/"
                  r"\langle A\rangle_{\rm true\,bin}$", fontsize=9)
    ax3.legend(fontsize=6.8, loc="upper left", ncol=2)
    ax3.tick_params(labelsize=8)

    # --- (d) estimator under a phi-dependent efficiency --------------------
    edges = np.linspace(0.0, 2 * np.pi, 25)
    i0, j0 = spots[0][2], spots[0][3]                      # 3x3 super-bin
    n_tot = float(proj.n_events[max(i0 - 1, 0):i0 + 2,
                                max(j0 - 1, 0):j0 + 2].sum())
    amp = float(obs["a_cos2phi"][spots[0][2], spots[0][3]])
    eff = lambda ph: 1.0 + 0.03 * np.cos(2 * ph) + 0.02 * np.cos(ph)  # noqa: E731
    rng = np.random.default_rng(20260824)
    mu1, _ = reco.expected_counts_by_fill(n_tot, [0.6], amp, edges,
                                          acceptance=eff)
    mu2, f2 = reco.expected_counts_by_fill(n_tot, [0.6, -1.2], amp, edges,
                                           acceptance=eff)
    mu0, _ = reco.expected_counts_by_fill(n_tot, [0.6], amp, edges)
    single = np.array([cos2phi_fit_binned(rng.poisson(mu1[0]), edges, 0.6)
                       for _ in range(400)])
    ideal = np.array([cos2phi_fit_binned(rng.poisson(mu0[0]), edges, 0.6)
                      for _ in range(400)])
    ratio = np.array([reco.harmonic_ratio_fit(rng.poisson(mu2), f2,
                                              [0.6, -1.2], edges)["amp"]
                      for _ in range(400)])
    bins = np.linspace(amp - 6e-3, amp + 6e-2, 160)
    ax4.hist(1e3 * ideal, bins=1e3 * bins, color=C_GREEN, alpha=0.55,
             label=r"single fill, uniform $\varepsilon$: $\hat A$ = "
                   r"%.2f $\pm$ %.2f" % (1e3 * ideal.mean(), 1e3 * ideal.std()))
    ax4.hist(1e3 * single, bins=1e3 * bins, color=C_VERM, alpha=0.6,
             label=r"single fill, $\varepsilon=1+0.03\cos2\phi'+0.02\cos\phi'$:"
                   "\n" r"$\hat A$ = %.2f $\pm$ %.2f (bias 0.03/$P_{zz}$)"
                   % (1e3 * single.mean(), 1e3 * single.std()))
    ax4.hist(1e3 * ratio, bins=1e3 * bins, color=C_BLUE, alpha=0.6,
             label=r"two fills $P_{zz}=+0.6/-1.2$, same $\varepsilon$, "
                   r"ratio fit:" "\n" r"$\hat A$ = %.2f $\pm$ %.2f"
                   % (1e3 * ratio.mean(), 1e3 * ratio.std()))
    ax4.axvline(1e3 * amp, color="black", lw=1.2)
    ax4.set_xlabel(r"$\hat A^{\cos2\phi}$ $[\times10^{-3}]$   "
                   r"(sweet spot 1, $N=%.1f\times10^{8}$, 1 yr)" % (n_tot / 1e8))
    ax4.set_ylabel("pseudo-experiments")
    ax4.set_title("(d) single-fill fit vs spin-state-sorted ratio under a "
                  "spin-axis-aligned efficiency", fontsize=9)
    ax4.legend(fontsize=6.6, loc="upper right")
    ax4.tick_params(labelsize=8)
    ins = ax4.inset_axes([0.26, 0.20, 0.36, 0.38])
    zb = np.linspace(amp - 7e-4, amp + 7e-4, 41)
    ins.hist(1e3 * ideal, bins=1e3 * zb, color=C_GREEN, alpha=0.55)
    ins.hist(1e3 * ratio, bins=1e3 * zb, color=C_BLUE, alpha=0.6)
    ins.axvline(1e3 * amp, color="black", lw=1.0)
    ins.set_title("zoom: uniform-$\\varepsilon$ single fill (green) vs "
                  "ratio (blue, 1.5x narrower)", fontsize=6.3)
    ins.tick_params(labelsize=6)

    fig.suptitle("Inclusive cos 2φ measurement: what the scattered electron "
                 "can and cannot deliver, and what the estimator must cancel",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = pathlib.Path(outdir) / "reco_chain_inclusive_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print("spot-1 super-bin N_1yr = %.3g, A = %.4g; single-fill sigma "
          "%.2e (ideal %.2e), ratio %.2e, analytic ratio %.2e"
          % (n_tot, amp, single.std(), ideal.std(), ratio.std(),
             reco.err_harmonic_ratio(n_tot, [0.6, -1.2])))
    for k, rr in enumerate(res_list):
        print("migration, res=%s: " % rr
              + "  ".join("spot%d pur=%.2f eff=%.2f Aratio=%.3f Nratio=%.2f"
                          % (n + 1, *res[k][n]) for n in range(len(spots))))


def coherent_figure(outdir):
    config = beams.default_configs("6Li")[1]
    p_u = config.ion_momentum_per_nucleon
    slope_b = 50.0
    sig = reco.SIGMA_THETA_HA
    cut = reco.tag_pt_cut(sig, p_u)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12.4, 8.8))

    # --- (a) recoil pT plane with the cutout -------------------------------
    g = np.linspace(-0.5, 0.5, 400)
    px, py = np.meshgrid(g, g)
    dens = np.exp(-slope_b * (px ** 2 + py ** 2))
    ax1.imshow(dens, extent=(-0.5, 0.5, -0.5, 0.5), origin="lower",
               cmap="Blues", vmin=0, vmax=1.2, zorder=0)
    r_asp = 1.5
    ax1.add_patch(Rectangle((-cut, -cut * r_asp), 2 * cut, 2 * cut * r_asp,
                            fill=True, fc="white", ec=C_VERM, lw=1.8,
                            zorder=2, alpha=0.85))
    ax1.add_patch(Ellipse((0, 0), 2 * cut, 2 * cut, fill=False, ec="0.3",
                          ls=":", lw=1.3, zorder=3))
    ax1.annotate("", xy=(0.0, 0.46), xytext=(0.0, -0.46),
                 arrowprops=dict(arrowstyle="<->", color=C_GREEN, lw=1.6),
                 zorder=4)
    ax1.annotate(r"alignment axis $\hat n$ (vertical, headless)",
                 xy=(0.015, 0.42), color=C_GREEN, fontsize=7)
    phi_ex = 0.62
    pt_ex = 0.33
    ax1.annotate("", xy=(pt_ex * np.cos(phi_ex), pt_ex * np.sin(phi_ex)),
                 xytext=(0, 0), arrowprops=dict(arrowstyle="->",
                                                color="black", lw=1.4),
                 zorder=4)
    ax1.annotate(r"$\vec p_T(^6$Li$)$, $|t|\simeq p_T^2$, azimuth $\phi_t$",
                 xy=(0.24, 0.25), fontsize=7)
    ax1.annotate("cutout: $10\\sigma$ envelope\n$|p_x|<%.2f$, $|p_y|<%.2f$ GeV "
                 "($r=%.1f$)\nbeam-blind" % (cut, cut * r_asp, r_asp),
                 xy=(0, -0.14), fontsize=6.8, ha="center", va="center",
                 color=C_VERM, zorder=5)
    ax1.annotate("circular $p_T>%.2f$ GeV\n(routing-code cut)" % cut,
                 xy=(-0.47, -0.45), fontsize=6.8, color="0.3")
    ax1.annotate(r"tagged: $e^{-B p_T^2}$ tail, $B=%g$ GeV$^{-2}$" % slope_b,
                 xy=(-0.47, 0.42), fontsize=7, color=C_BLUE)
    ax1.set_xlabel(r"$p_x$ [GeV]  (horizontal, crossing-angle plane)")
    ax1.set_ylabel(r"$p_y$ [GeV]  (vertical)")
    ax1.set_title(r"(a) recoil $\vec p_T$ plane at the Roman Pots, "
                  r"%s, $\sigma_\theta=%.0f\,\mu$rad" % (config.label(),
                                                         1e6 * sig),
                  fontsize=9)
    ax1.tick_params(labelsize=8)

    # --- (b) fake a2 vs aspect ratio ---------------------------------------
    rr = np.linspace(0.6, 2.0, 141)
    for shape, ls, lab in (("rectangle", "-", "rectangular cutout"),
                           ("ellipse", "--", "elliptical cutout")):
        a2 = [reco.rp_hole_acceptance(slope_b, cut, cut * r, shape)["a2"]
              for r in rr]
        ax2.plot(rr, a2, ls, color=C_VERM, lw=1.8, label=lab)
    a4 = [reco.rp_hole_acceptance(slope_b, cut, cut * r)["a4"] for r in rr]
    ax2.plot(rr, a4, "-", color=C_GREY, lw=1.2,
             label=r"rectangle, $\cos4\phi_t$ coefficient")
    ax2.axhspan(0.018, 0.059, color=C_BLUE, alpha=0.18, lw=0,
                label=r"anchored $\langle a_2\rangle$ band, $P_{zz}=0.6$")
    ax2.axhline(0.036, color=C_BLUE, lw=1.2)
    for lvl, lab in ((5 * 0.6 * 1.9e-3, r"$5\sigma$ floor, 1 yr"),
                     (5 * 0.6 * 6e-4, r"$5\sigma$ floor, 10 yr")):
        ax2.axhline(lvl, color="black", lw=0.7, ls=":")
        ax2.annotate(lab, xy=(1.92, lvl * 1.1), fontsize=6.5, ha="right")
    ax2.axhline(0, color="0.8", lw=0.6)
    ax2.set_yscale("symlog", linthresh=1e-3)
    ax2.set_xlabel(r"cutout aspect ratio $r=c_y/c_x$ ($c_x=%.2f$ GeV)" % cut)
    ax2.set_ylabel(r"$\langle\cos2\phi_t\rangle$ of the tagged sample")
    ax2.set_title(r"(b) fake $\cos2(\phi_t-\phi_S)$ from the cutout "
                  r"geometry vs the physics amplitude", fontsize=9)
    ax2.legend(fontsize=6.8, loc="lower right")
    ax2.tick_params(labelsize=8)

    # --- (c) acceptance vs beam momentum per nucleon -----------------------
    pu = np.linspace(15.0, 140.0, 300)
    for sg, col, lab in ((reco.SIGMA_THETA_HA, C_BLUE,
                          r"$\sigma_\theta=73\,\mu$rad (0.20 GeV at 275 GeV p)"),
                         (reco.SIGMA_THETA_HD, C_VERM,
                          r"$\sigma_\theta=149\,\mu$rad (0.41 GeV at 275 GeV p)")):
        c = reco.tag_pt_cut(sg, pu)
        ax3.plot(pu, np.exp(-slope_b * c * c), "-", color=col, lw=1.8,
                 label=lab)
        ax3.fill_between(pu, np.exp(-60.0 * c * c), np.exp(-40.0 * c * c),
                         color=col, alpha=0.15, lw=0)
    ax3.axhline(np.exp(-slope_b * 0.04), color="0.3", ls=":", lw=1.0,
                label="constant 0.20 GeV cut (13.5%, current code)")
    for cfg_i, mk in zip(beams.default_configs("6Li"), ("s", "o", "^")):
        pi = cfg_i.ion_momentum_per_nucleon
        ci = reco.tag_pt_cut(reco.SIGMA_THETA_HA, pi)
        ax3.plot(pi, np.exp(-slope_b * ci * ci), mk, color="black", ms=6,
                 mfc="white" if cfg_i is not config else "black")
        ax3.annotate(cfg_i.label(), xy=(pi, np.exp(-slope_b * ci * ci)),
                     xytext=(6, 4), textcoords="offset points", fontsize=6.8)
    ax3.set_yscale("log")
    ax3.set_ylim(1e-8, 1.5)
    ax3.set_xlabel(r"$^6$Li beam momentum per nucleon $p_u$ [GeV]")
    ax3.set_ylabel(r"tag acceptance $\exp[-B\,(10\sigma_\theta\,6p_u)^2]$")
    ax3.set_title(r"(c) the near-beam cut is ANGULAR: $p_T^{\rm cut}="
                  r"10\sigma_\theta\,A\,p_u$ (bands: $B=40$–$60$)", fontsize=9)
    ax3.legend(fontsize=6.8, loc="lower left")
    ax3.tick_params(labelsize=8)

    # --- (d) RP emulation: t_reco vs t_true, phi_t resolution --------------
    _, p4 = reco.beam_fourvectors(config)
    rng = np.random.default_rng(20260824)
    n = 400000
    t_true = -rng.exponential(1.0 / slope_b, n)
    x_pom = 10 ** rng.uniform(-3, -2, n)
    t_min = reco.mdot(p4, p4) * x_pom ** 2 / (1 - x_pom)
    ok = -t_true > t_min
    t_true, x_pom = t_true[ok], x_pom[ok]
    phi_t = rng.uniform(0, 2 * np.pi, t_true.size)
    pp = reco.recoil_fourvector(t_true, phi_t, x_pom, p4)
    m = reco.rp_measure(pp, p4, (sig, sig), rng=rng)
    acc = m["accepted"]
    tb = np.linspace(0.0, 0.3, 61)
    ax4.hist(-t_true, bins=tb, histtype="step", color="0.3", lw=1.2,
             label=r"true $|t|$, all coherent recoils")
    ax4.hist(-t_true[acc], bins=tb, histtype="step", color=C_BLUE, lw=1.6,
             label=r"true $|t|$, RP-accepted (%.1f%%)" % (100 * acc.mean()))
    ax4.hist(-m["t_reco"][acc], bins=tb, histtype="stepfilled",
             color=C_VERM, alpha=0.35, lw=1.2, ec=C_VERM,
             label=r"reconstructed $|t|=p_T^2$ (divergence-smeared)")
    ax4.set_yscale("log")
    ax4.set_xlabel(r"$|t|$ [GeV$^2$]")
    ax4.set_ylabel("recoils")
    ax4.set_title(r"(d) Roman-Pot emulation: $10\sigma$ square cutout, "
                  r"$\delta p_T=6p_u\sigma_\theta=%.0f$ MeV" % (
                      1e3 * 6 * p_u * sig), fontsize=9)
    ax4.legend(fontsize=7, loc="upper right")
    ax4.tick_params(labelsize=8)
    ins = ax4.inset_axes([0.47, 0.40, 0.50, 0.28])
    dphi = np.angle(np.exp(1j * (m["phi_t"][acc] - phi_t[acc])))
    ptb = np.linspace(cut, 0.5, 9)
    ptc = 0.5 * (ptb[:-1] + ptb[1:])
    sd = [np.std(dphi[(m["pT_true"][acc] >= lo) & (m["pT_true"][acc] < hi)])
          for lo, hi in zip(ptb[:-1], ptb[1:])]
    ins.plot(ptc, sd, "o-", color=C_BLUE, ms=3, lw=1.2)
    ins.plot(ptc, 6 * p_u * sig / ptc, "--", color="0.4", lw=1.0)
    ins.set_xlabel(r"$p_T$ [GeV]", fontsize=6.5, labelpad=1)
    ins.set_ylabel(r"$\sigma(\phi_t)$ [rad]", fontsize=6.5, labelpad=1)
    ins.tick_params(labelsize=6)
    ins.set_title(r"$\phi_t$ resolution ($\langle\cos2\delta\phi\rangle"
                  r"\approx e^{-2\sigma^2}$)", fontsize=6.5)

    fig.suptitle("Coherent intact-⁶Li channel: the recoil measurement, its "
                 "acceptance geometry, and the angular near-beam cut",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = pathlib.Path(outdir) / "reco_chain_coherent_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    for r in (1.0, 1.25, 1.5):
        h = reco.rp_hole_acceptance(slope_b, cut, cut * r)
        print("rectangle r=%.2f: acc=%.3f a2=%.3f a4=%.3f" % (r, h["acc"],
                                                              h["a2"], h["a4"]))
    for cfg_i in beams.default_configs("6Li"):
        c = reco.tag_pt_cut(reco.SIGMA_THETA_HA, cfg_i.ion_momentum_per_nucleon)
        print("%s: pT_cut(HA) = %.3f GeV, acc = %.3g" % (
            cfg_i.label(), c, np.exp(-slope_b * c * c)))
    print("RP emulation acceptance (square, isotropic): %.4f; phi_t sigma "
          "at pT=0.25: %.3f rad" % (acc.mean(), 6 * p_u * sig / 0.25))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    inclusive_figure(outdir)
    coherent_figure(outdir)


if __name__ == "__main__":
    main()
