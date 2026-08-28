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

reco_chain_coherent_6Li.png  (at the lithium tagging optics, 5 x 40.8)
  (a) the recoil transverse-momentum plane with the three cutouts: the
      Yellow Report high-acceptance 10 sigma envelope, the pot aperture
      measured in the ePIC geometry, and the tagging optics with the
      pots following the envelope (0.33 x 3.8 mrad);
  (b) the fake <cos 2beta> about the vertical spin axis and the tagged
      fraction versus the horizontal half-width of the cutout, against
      the anchored deformation band;
  (c) tagged fraction versus the horizontal half-angle per configuration,
      with the three cutouts marked -- the tag is an angle, so the pots
      must reach the envelope and the optics must shrink it;
  (d) Roman-Pot emulation at the tagging optics: reconstructed vs true
      |t| and the phi_t resolution from the anisotropic divergence.

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
    """The coherent measurement at the tagging optics of Report 1 Section
    6.1 (5 x 40.8: horizontal beta* x 50, pots following the 10 sigma
    envelope 0.33 x 3.8 mrad), against the Yellow Report high-acceptance
    optics and the pot aperture measured in the ePIC geometry
    (2026-08-28; the earlier version drew the legacy 73 microrad isotropic
    envelope and a 1.5-aspect slot at the mid configuration, where the tag
    has 8e-5 acceptance and panel (d) was empty)."""
    cfgs = beams.default_configs("6Li")
    config = cfgs[0]
    p_u = config.ion_momentum_per_nucleon
    p_ion = config.ion.A * p_u
    slope_b = 50.0
    top = reco.tagging_optics_point(config, slope_b=slope_b)
    sx, sy = top["sigma_x_eff"], top["sigma_y"]
    cx, cy = top["env_x"] * p_ion, top["env_y"] * p_ion          # GeV
    apx, apy = reco.rp_aperture_for(p_u)
    ha_x, ha_y = reco.sigma_theta_for(config, "high-acceptance")
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12.4, 8.8))

    # --- (a) recoil pT plane with the three cutouts --------------------------
    g = np.linspace(-0.7, 0.7, 500)
    px, py = np.meshgrid(g, g)
    dens = np.exp(-slope_b * (px ** 2 + py ** 2))
    ax1.imshow(dens, extent=(-0.7, 0.7, -0.7, 0.7), origin="lower",
               cmap="Blues", vmin=0, vmax=1.2, zorder=0)
    for (hx, hy), col, ls, lab in (
            ((10 * ha_x * p_ion, 10 * ha_y * p_ion), "0.3", ":",
             "10σ envelope, YR high-acceptance optics (%.0f/%.0f μrad): %.2f × %.2f GeV"
             % (1e6 * ha_x, 1e6 * ha_y, 10 * ha_x * p_ion, 10 * ha_y * p_ion)),
            ((apx * p_ion, apy * p_ion), C_PURPLE, "--",
             "pot aperture measured in the ePIC geometry: %.2f × %.2f GeV"
             % (apx * p_ion, apy * p_ion)),
            ((cx, cy), C_VERM, "-",
             "tagging optics, pots following the 10σ envelope: %.2f × %.2f GeV" % (cx, cy))):
        ax1.add_patch(Rectangle((-hx, -hy), 2 * hx, 2 * hy, fill=False, ec=col,
                                ls=ls, lw=1.8, zorder=3, label=lab))
    ax1.add_patch(Rectangle((-cx, -cy), 2 * cx, 2 * cy, fill=True, fc="white",
                            ec="none", zorder=2, alpha=0.85))
    ax1.annotate("", xy=(0.0, 0.66), xytext=(0.0, -0.66),
                 arrowprops=dict(arrowstyle="<->", color=C_GREEN, lw=1.6), zorder=4)
    ax1.annotate(r"alignment axis $\hat n$ (vertical, headless)", xy=(0.02, 0.60),
                 color=C_GREEN, fontsize=7)
    ax1.annotate(r"tagged: $e^{-B p_T^2}$ tail, $B=%g$ GeV$^{-2}$" % slope_b,
                 xy=(-0.68, 0.62), fontsize=7, color=C_BLUE)
    ax1.set_xlim(-0.7, 0.7)
    ax1.set_ylim(-0.7, 0.7)
    ax1.set_xlabel(r"$p_x$ [GeV]  (horizontal, crossing-angle plane)")
    ax1.set_ylabel(r"$p_y$ [GeV]  (vertical, along $\hat n$)")
    ax1.set_title(r"(a) recoil $\vec p_T$ plane at the Roman Pots, %s" % config.label(),
                  fontsize=9)
    ax1.legend(fontsize=6.4, loc="lower left")
    ax1.tick_params(labelsize=8)

    # --- (b) fake <cos 2beta> vs the horizontal half-width ------------------
    hx_gev = np.logspace(np.log10(0.03), np.log10(1.0), 160)
    acc_b, a2_b = [], []
    for h in hx_gev:
        r = reco.rp_hole_acceptance(slope_b, h, cy)
        acc_b.append(r["acc"])
        a2_b.append(-r["a2"])          # about the VERTICAL spin axis
    ax2.plot(hx_gev, a2_b, "-", color=C_VERM, lw=1.8,
             label=r"$\langle\cos 2\beta\rangle$ of the tagged sample "
                   r"($c_y$ = %.2f GeV, vertical at high acceptance)" % cy)
    ax2.plot(hx_gev, acc_b, "-", color=C_BLUE, lw=1.4, label="tagged fraction")
    ax2.axhspan(0.018, 0.059, color=C_GREEN, alpha=0.18, lw=0,
                label=r"anchored deformation term $a_t$ band at $P_{zz}=0.6$")
    ax2.axhline(0.0, color="0.8", lw=0.6)
    for h, col, lab in ((cx, C_VERM, "tagging optics"), (apx * p_ion, C_PURPLE, "pot aperture"),
                        (10 * ha_x * p_ion, "0.3", "YR HA envelope")):
        ax2.axvline(h, color=col, lw=0.9, ls="--")
        ax2.annotate(lab, xy=(h, 0.92), xytext=(3, 0), textcoords="offset points",
                     fontsize=6.6, color=col, rotation=90, va="top")
    ax2.set_xscale("log")
    ax2.set_ylim(-1.0, 1.0)
    ax2.set_xlabel(r"horizontal half-width of the cutout $c_x$ [GeV]")
    ax2.set_ylabel("coefficient / fraction")
    ax2.set_title(r"(b) what the cutout fakes in the recoil azimuth, and what it leaves",
                  fontsize=9)
    ax2.legend(fontsize=6.6, loc="lower right")
    ax2.tick_params(labelsize=8)

    # --- (c) acceptance vs the horizontal envelope, per configuration -------
    thx = np.logspace(np.log10(0.05), np.log10(3.0), 200)      # mrad
    for cfg_i, col, mk in zip(cfgs, (C_BLUE, C_VERM, C_GREEN), ("s", "o", "^")):
        pi = cfg_i.ion.A * cfg_i.ion_momentum_per_nucleon
        t_i = reco.tagging_optics_point(cfg_i, slope_b=slope_b)
        hx_i, hy_i = reco.sigma_theta_for(cfg_i, "high-acceptance")
        ax_i, ay_i = reco.rp_aperture_for(cfg_i.ion_momentum_per_nucleon)
        acc = [reco.rp_hole_acceptance(slope_b, 1e-3 * th * pi, t_i["env_y"] * pi)["acc"]
               for th in thx]
        ax3.plot(thx, acc, "-", color=col, lw=1.6,
                 label="%s, vertical envelope %.2f mrad" % (cfg_i.label(), 1e3 * t_i["env_y"]))
        ax3.plot([1e3 * t_i["env_x"]], [t_i["acceptance"]], mk, color=col, ms=7,
                 mfc="white", mew=1.6)
        ax3.plot([1e3 * ax_i], [reco.rp_hole_acceptance(slope_b, ax_i * pi, t_i["env_y"] * pi)["acc"]],
                 mk, color=col, ms=6)
        ax3.plot([1e4 * hx_i], [reco.rp_hole_acceptance(slope_b, 10 * hx_i * pi, t_i["env_y"] * pi)["acc"]],
                 mk, color=col, ms=6, alpha=0.45)
    ax3.plot([], [], "k", marker="o", mfc="white", ls="none", label="tagging optics (pots follow)")
    ax3.plot([], [], "k", marker="o", ls="none", label="pot aperture measured in the ePIC geometry")
    ax3.plot([], [], "k", marker="o", ls="none", alpha=0.45, label="10σ envelope of the YR high-acceptance optics")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_ylim(1e-9, 1.5)
    ax3.set_xlabel(r"horizontal half-width of the cutout at the IP, $\theta_x$ [mrad]")
    ax3.set_ylabel(r"tagged fraction (rectangle, $B=50$ GeV$^{-2}$)")
    ax3.set_title(r"(c) the tag is an angle: pots must reach the envelope, and the envelope must shrink",
                  fontsize=9)
    ax3.legend(fontsize=6.4, loc="lower left")
    ax3.tick_params(labelsize=8)

    # --- (d) RP emulation at the tagging optics ------------------------------
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
    m = reco.rp_measure(pp, p4, (sx, sy), rng=rng, cut_theta_xy=(top["env_x"], top["env_y"]))
    acc = m["accepted"]
    tb = np.linspace(0.0, 0.3, 61)
    ax4.hist(-t_true, bins=tb, histtype="step", color="0.3", lw=1.2,
             label=r"true $|t|$, all coherent recoils")
    ax4.hist(-t_true[acc], bins=tb, histtype="step", color=C_BLUE, lw=1.6,
             label=r"true $|t|$, tagged (%.1f%%)" % (100 * acc.mean()))
    ax4.hist(-m["t_reco"][acc], bins=tb, histtype="stepfilled",
             color=C_VERM, alpha=0.35, lw=1.2, ec=C_VERM,
             label=r"reconstructed $|t|=p_T^2$ (divergence-smeared)")
    ax4.set_yscale("log")
    ax4.set_xlabel(r"$|t|$ [GeV$^2$]")
    ax4.set_ylabel("recoils")
    ax4.set_title(r"(d) Roman-Pot emulation at the tagging optics: $\sigma_\theta$ = %.0f / %.0f μrad, "
                  r"$\delta p_T$ = %.0f / %.0f MeV (h / v)"
                  % (1e6 * sx, 1e6 * sy, 1e3 * p_ion * sx, 1e3 * p_ion * sy), fontsize=8.6)
    ax4.legend(fontsize=7, loc="upper right")
    ax4.tick_params(labelsize=8)
    ins = ax4.inset_axes([0.47, 0.40, 0.50, 0.28])
    dphi = np.angle(np.exp(1j * (m["phi_t"][acc] - phi_t[acc])))
    ptb = np.linspace(0.08, 0.6, 12)
    ptc = 0.5 * (ptb[:-1] + ptb[1:])
    sd = [np.std(dphi[(m["pT_true"][acc] >= lo) & (m["pT_true"][acc] < hi)])
          for lo, hi in zip(ptb[:-1], ptb[1:])]
    ins.plot(ptc, sd, "o-", color=C_BLUE, ms=3, lw=1.2, label="measured")
    ins.plot(ptc, p_ion * sy / ptc, "--", color="0.4", lw=1.0, label=r"$\delta p_{T,y}/p_T$")
    ins.set_xlabel(r"$p_T$ [GeV]", fontsize=6.5, labelpad=1)
    ins.set_ylabel(r"$\sigma(\phi_t)$ [rad]", fontsize=6.5, labelpad=1)
    ins.tick_params(labelsize=6)
    ins.legend(fontsize=5.8, loc="upper right")
    ins.set_title(r"$\phi_t$ resolution of the tagged recoils", fontsize=6.5)

    fig.suptitle("Coherent intact-⁶Li channel: the recoil measurement, the cutout geometry and "
                 "the angular near-beam envelope, at the lithium tagging optics", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = pathlib.Path(outdir) / "reco_chain_coherent_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print("tagging optics at %s: r_h = %.1f, sigma h/v = %.0f/%.0f microrad, envelope %.2f x %.2f mrad "
          "= %.3f x %.3f GeV, acceptance %.3f, L/L_HA = 1/%.1f"
          % (config.label(), top["r_h"], 1e6 * sx, 1e6 * sy, 1e3 * top["env_x"], 1e3 * top["env_y"],
             cx, cy, top["acceptance"], 1.0 / top["lumi_fraction"]))
    h = reco.rp_hole_acceptance(slope_b, cx, cy)
    print("fake <cos 2beta> about the vertical axis at the tagging cutout: %.3f (a4 %.3f); "
          "at the measured aperture: %.3f" % (-h["a2"], h["a4"],
                                             -reco.rp_hole_acceptance(slope_b, apx * p_ion, apy * p_ion)["a2"]))
    for cfg_i in cfgs:
        pi = cfg_i.ion.A * cfg_i.ion_momentum_per_nucleon
        t_i = reco.tagging_optics_point(cfg_i, slope_b=slope_b)
        hx_i, hy_i = reco.sigma_theta_for(cfg_i, "high-acceptance")
        ax_i, ay_i = reco.rp_aperture_for(cfg_i.ion_momentum_per_nucleon)
        print("%s: YR-HA 10 sigma envelope %.2f x %.2f mrad -> acc %.2e; measured aperture %.2f x %.2f mrad "
              "-> acc %.2e; tagging optics %.2f x %.2f mrad -> acc %.3f"
              % (cfg_i.label(), 1e4 * hx_i, 1e4 * hy_i,
                 reco.rp_hole_acceptance(slope_b, 10 * hx_i * pi, 10 * hy_i * pi)["acc"],
                 1e3 * ax_i, 1e3 * ay_i,
                 reco.rp_hole_acceptance(slope_b, max(ax_i, 10 * hx_i) * pi, max(ay_i, 10 * hy_i) * pi)["acc"],
                 1e3 * t_i["env_x"], 1e3 * t_i["env_y"], t_i["acceptance"]))
    print("RP emulation at the tagging optics: acceptance %.4f; phi_t sigma at pT = 0.15 / 0.30 GeV: "
          "%.3f / %.3f rad (measured), %.3f / %.3f (delta p_T,y / p_T)"
          % (acc.mean(), sd[1], sd[5], p_ion * sy / ptc[1], p_ion * sy / ptc[5]))


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
