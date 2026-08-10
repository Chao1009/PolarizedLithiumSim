#!/usr/bin/env python3
"""Money-plot variant (2026-07-15): gluonometry discovery reach vs beam energy.

Parent script: fastsim/scripts/money_delta.py (money plot 3, plan step 1.3.3).

Whereas money_delta.py shows required luminosity vs Delta/F1 scale at fixed
EIC beam energy points, this script visualises how L_5sigma depends on the
actual beam energies (electron or 6Li ion) for several Delta/F1 scale values,
and also shows L_5σ vs √s, kinematic (x, Q²) coverage, and yield vs azimuthal
angle φ for 6Li DIS.

  Figure 1  (money_delta_20260715_vs_Ee.png)
      x-axis: electron energy E_e [GeV], linear
      y-axis: L_5sigma [fb^-1/nucleon], log
      ion fixed at top rigidity-limited momentum per nucleon
      P_zz = 0.80, three Delta/F1 scale curves

  Figure 2  (money_delta_20260715_vs_Eion.png)
      x-axis: 6Li ion momentum per nucleon p_ion [GeV/u], linear
      y-axis: L_5sigma [fb^-1/nucleon], log
      E_e fixed at 18 GeV
      P_zz = 0.80, three Delta/F1 scale curves

  Figure 3  (money_delta_20260715_L5sig_vs_sqrts.png)
      Required 5σ luminosity L_5σ vs √s per nucleon, derived from
      sig2_per_fb_at via L_5σ = 25 / sig2.
      One curve per (E_e, Δ/F1 scale) combination: E_e ∈ {5, 10, 18} GeV ×
      scale ∈ {1e-3, 3e-3, 1e-2} = 9 curves total, sweeping p_ion from
      20 GeV/u to the rigidity limit and computing √s at each point.
      Gold band marks the 1–100 fb⁻¹/u plausible EIC program range.
      P_zz = 0.80 fixed.

  Figure 4  (money_delta_20260715_xQ2_coverage.png)
      Accepted (x, Q²) bin coverage shown as event-count pcolormesh in a 2×2
      subplot grid at four representative beam configurations covering low and
      high ends of the energy scan.

  Figure 5  (money_delta_20260715_yield_vs_phi.png)
      Two-panel figure: top panel shows expected 6Li DIS yield dN/dφ vs
      azimuthal angle φ; bottom panel shows the asymmetry A(φ) with a ±1σ
      statistical band.  Fixed parameters: L = 10 fb⁻¹/u, E_e = 18 GeV,
      p_ion at the rigidity limit, P_zz = 0.80.
      Physics setup: the spin-1 tensor-polarized DIS cross section (transverse
      tensor polarization, θ_m = 90°, λ_m = ±1) has a cos(2φ) term
      proportional to Δ(x, Q²), as given by the Hoodbhoy-Jaffe-Manohar master
      formula (NPB 312:571).  The φ-differential normalized bin yield is:
          dN/dφ (bin) = (N_bin / 2π) · [1 + P_zz · A_cos2φ(bin) · cos(2φ)]
      where A_cos2φ = −((1−y)/y²) · Δ / D_φ from asymmetries.a_cos2phi, and
      N_bin is the accepted event count from project_rates at L = 10 fb⁻¹/u.
      Three curves for Δ/F1 scale ∈ {1e-3, 3e-3, 1e-2}, plus a flat reference
      (top panel).  Bottom panel shows A(φ) = P_zz · ⟨A_cos2φ⟩ · cos(2φ)
      with a gray ±1σ band from err_cos2phi_amplitude(N_total, P_zz).

Significance formula (same as money_delta.py):
    sig^2 = sum_bins A_bin^2 * P_zz^2 * N_bin / 2
    L_5sig = 25 / sig^2   (per 1 fb^-1/u reference)
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.beams import LI6, BeamConfig
from polli_fastsim.asymmetries import a_cos2phi, phi_averaged_density, err_cos2phi_amplitude
from polli_fastsim.inputs import get_backends
from polli_fastsim.kinematics import sqrt_s
from polli_fastsim.polarized import toy_delta_gluon
from polli_fastsim.structure import NuclearF2

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ── fixed run parameters ──────────────────────────────────────────────────────
PZZ = 0.80
SCALES = [1e-3, 3e-3, 1e-2]
SCALE_COLORS = ("royalblue", "darkorange", "crimson")

EE_SWEEP = np.linspace(5.0, 18.0, 14)          # GeV, Figure 1
PION_TOP = LI6.momentum_per_nucleon_max         # ~137.5 GeV/u
PION_SWEEP = np.linspace(20.0, PION_TOP, 14)   # GeV/u, Figure 2
EE_FIXED = 18.0                                 # GeV, Figure 2

# Figure 3 parameters
EE_CURVES = [5.0, 10.0, 18.0]                  # GeV, one curve each
EE_CURVE_COLORS = ("steelblue", "seagreen", "darkorchid")
LUMI_REF = 10.0                                 # fb^-1/nucleon, Figures 3-5

# Figure 4: four representative (E_e, p_ion) configurations
FIG4_CONFIGS = [
    (5.0,    41.0 * 0.5,      "low E_e, low p_ion"),
    (5.0,    PION_TOP,        "low E_e, top p_ion"),
    (18.0,   41.0 * 0.5,     "top E_e, low p_ion"),
    (18.0,   PION_TOP,       "top E_e, top p_ion"),
]

# Figure 5 parameters
EE_FIG5 = 18.0
PION_FIG5 = PION_TOP
PZZ_FIG5 = 0.80
PHI_VALS = np.linspace(0.0, 2.0 * np.pi, 200)


# ── physics helper ────────────────────────────────────────────────────────────

def sig2_per_fb_at(cfg, scale, pzz, base=None, min_events=10):
    """Significance^2 per fb^-1/nucleon at a reference Delta/F1 scale.

    Copied from money_delta.py.  The toy Delta shape is linear in `scale`,
    so reach curves follow L_5sig(s) = 25 / (sig2 * (s/scale)^2) analytically.
    Here we call it once per (cfg, scale) pair without the analytic shortcut so
    the sweep is self-contained and readable.
    """
    sc = fom.Scenario(lumi_fb_per_nucleon=1.0, pol_ion_tensor=pzz)
    nf2_in = NuclearF2(cfg.ion, base=base) if base is not None else None
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_in)
    nf2 = proj.extras["nf2"]
    f1 = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2 = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y = proj.extras["y"]
    delta = toy_delta_gluon(proj.x, proj.q2, f1, scale=scale)
    amp = a_cos2phi(delta, f1, f2, proj.x, y)
    use = proj.accepted & (proj.n_events >= min_events)
    return np.where(use, amp**2 * pzz**2 * proj.n_events / 2.0, 0.0).sum()


# ── figure builders ───────────────────────────────────────────────────────────

def _add_plausible_band(ax):
    """Gold axhspan for the 1–100 fb^-1/u plausible EIC program."""
    ax.axhspan(1, 100, color="gold", alpha=0.12,
               label=r"1$-$100 fb$^{-1}$/u (plausible program)")


def build_fig1(backends):
    """Figure 1: L_5sig vs E_e at top ion energy, P_zz=0.80."""
    p_ion = round(PION_TOP, 1)
    tag = backends["tag"]
    base = backends["base"]

    fig, ax = plt.subplots(figsize=(7, 5))
    summary_rows = []

    for scale, color in zip(SCALES, SCALE_COLORS):
        l5sig_vals = []
        for ee in EE_SWEEP:
            cfg = BeamConfig(electron_energy=ee, ion=LI6,
                             ion_momentum_per_nucleon=p_ion)
            s2 = sig2_per_fb_at(cfg, scale, PZZ, base=base)
            l5 = 25.0 / np.maximum(s2, 1e-30)
            l5sig_vals.append(l5)
        l5sig_vals = np.array(l5sig_vals)
        label = fr"$\Delta/F_1 = {scale:.0e}$"
        ax.plot(EE_SWEEP, l5sig_vals, "-o", color=color, lw=1.5,
                ms=4, label=label)
        summary_rows.append((scale, EE_SWEEP[-1], l5sig_vals[-1]))

    _add_plausible_band(ax)
    ax.set_yscale("log")
    ax.set_xlabel(r"$E_e$ [GeV]", fontsize=11)
    ax.set_ylabel(r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]", fontsize=11)
    ax.set_title(
        rf"$^6$Li, top ion energy ({p_ion:.1f} GeV/u), $P_{{zz}}={PZZ:g}$"
        f"\n(cos 2$\\phi$ amplitude, all bins combined; {tag.upper()} inputs)",
        fontsize=10,
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig, summary_rows


def build_fig2(backends):
    """Figure 2: L_5sig vs p_ion at E_e=18 GeV, P_zz=0.80."""
    tag = backends["tag"]
    base = backends["base"]

    fig, ax = plt.subplots(figsize=(7, 5))
    summary_rows = []

    for scale, color in zip(SCALES, SCALE_COLORS):
        l5sig_vals = []
        for p_ion in PION_SWEEP:
            cfg = BeamConfig(electron_energy=EE_FIXED, ion=LI6,
                             ion_momentum_per_nucleon=p_ion)
            s2 = sig2_per_fb_at(cfg, scale, PZZ, base=base)
            l5 = 25.0 / np.maximum(s2, 1e-30)
            l5sig_vals.append(l5)
        l5sig_vals = np.array(l5sig_vals)
        label = fr"$\Delta/F_1 = {scale:.0e}$"
        ax.plot(PION_SWEEP, l5sig_vals, "-o", color=color, lw=1.5,
                ms=4, label=label)
        summary_rows.append((scale, PION_SWEEP[-1], l5sig_vals[-1]))

    _add_plausible_band(ax)
    ax.set_yscale("log")
    ax.set_xlabel(r"$p_\mathrm{ion}$ [GeV/u]", fontsize=11)
    ax.set_ylabel(r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]", fontsize=11)
    ax.set_title(
        rf"$^6$Li, $E_e = {EE_FIXED:g}$ GeV, $P_{{zz}}={PZZ:g}$"
        f"\n(cos 2$\\phi$ amplitude, all bins combined; {tag.upper()} inputs)",
        fontsize=10,
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig, summary_rows


def plot_l5sig_vs_sqrts(backends, outdir):
    """Figure 3: required 5σ luminosity L_5σ vs √s per nucleon, P_zz=0.80.

    Nine curves: E_e ∈ {5, 10, 18} GeV × Δ/F1 scale ∈ {1e-3, 3e-3, 1e-2}.
    Color encodes E_e; linestyle encodes scale.
    L_5σ = 25 / sig2_per_fb_at(...).
    """
    tag = backends["tag"]
    base = backends["base"]
    pzz = 0.80

    # linestyle encodes scale
    scale_ls = {1e-3: "-", 3e-3: "--", 1e-2: ":"}

    fig, ax = plt.subplots(figsize=(7, 5))

    pion_sweep = np.linspace(20.0, LI6.momentum_per_nucleon_max, 14)

    for ee, color in zip(EE_CURVES, EE_CURVE_COLORS):
        for scale in SCALES:
            ls = scale_ls[scale]
            sqrts_vals = []
            l5sig_vals = []
            for p_ion in pion_sweep:
                cfg = BeamConfig(electron_energy=ee, ion=LI6,
                                 ion_momentum_per_nucleon=p_ion)
                s2 = sig2_per_fb_at(cfg, scale, pzz, base=base)
                l5 = 25.0 / max(s2, 1e-30)
                sqrts_vals.append(sqrt_s(ee, p_ion))
                l5sig_vals.append(l5)

            sqrts_vals = np.array(sqrts_vals)
            l5sig_vals = np.array(l5sig_vals)
            label = rf"$E_e={ee:g}$ GeV, $\Delta/F_1={scale:.0e}$"
            ax.plot(sqrts_vals, l5sig_vals, ls, color=color, lw=1.5,
                    ms=4, label=label)

    ax.axhspan(1, 100, color="gold", alpha=0.12,
               label=r"1$-$100 fb$^{-1}$/u (plausible program)")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\sqrt{s}$ per nucleon [GeV]", fontsize=11)
    ax.set_ylabel(r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]", fontsize=11)
    ax.set_title(
        rf"$L_{{5\sigma}}$ vs $\sqrt{{s}}$, $^6$Li, $P_{{zz}}={pzz:g}$"
        f" ({tag.upper()} inputs)",
        fontsize=10,
    )
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()

    path = outdir / "money_delta_20260715_L5sig_vs_sqrts.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"wrote {path}")
    return path


def plot_xq2_coverage(backends, outdir):
    """Figure 4: accepted (x, Q²) coverage as pcolormesh in 2×2 panels."""
    tag = backends["tag"]
    base = backends["base"]
    sc = fom.Scenario(lumi_fb_per_nucleon=LUMI_REF)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes_flat = axes.flatten()

    for ax, (ee, p_ion, desc) in zip(axes_flat, FIG4_CONFIGS):
        cfg = BeamConfig(electron_energy=ee, ion=LI6,
                         ion_momentum_per_nucleon=p_ion)
        nf2 = NuclearF2(LI6, base=base)
        proj = fom.project_rates(cfg, sc, nuclear_f2=nf2)

        z = np.where(proj.accepted, proj.n_events, np.nan)
        # Guard all-zero / all-NaN panels
        z_finite = z[np.isfinite(z) & (z > 0)]
        if z_finite.size == 0:
            ax.set_title(f"$E_e$={ee:g}, $p_{{ion}}$={p_ion:g} (no accepted bins)",
                         fontsize=8)
            ax.set_xlabel("$x$", fontsize=9)
            ax.set_ylabel("$Q^2$ [GeV$^2$]", fontsize=9)
            continue

        vmin = z_finite.min()
        vmax = z_finite.max()
        if vmin <= 0:
            vmin = 1e-1

        sqrts_val = sqrt_s(ee, p_ion)
        pcm = ax.pcolormesh(
            proj.x_edges, proj.q2_edges, z.T,
            norm=LogNorm(vmin=vmin, vmax=vmax),
            cmap="viridis",
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("$x$", fontsize=9)
        ax.set_ylabel("$Q^2$ [GeV$^2$]", fontsize=9)
        ax.set_title(
            rf"$E_e={ee:g}$, $p_{{ion}}={p_ion:g}$, $\sqrt{{s}}={sqrts_val:.1f}$ GeV",
            fontsize=8,
        )
        fig.colorbar(pcm, ax=ax, label="Events / bin")

    fig.suptitle(
        r"$^6$Li accepted $(x,\,Q^2)$ coverage, $\mathcal{L}=10$ fb$^{-1}$/u"
        f"  ({tag.upper()} inputs)",
        fontsize=11,
    )
    fig.tight_layout()

    path = outdir / "money_delta_20260715_xQ2_coverage.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"wrote {path}")
    return path


def plot_yield_vs_phi(backends, outdir):
    """Figure 5: two-panel dN/dφ (top) and A(φ) with ±1σ band (bottom).

    Fixed setup: E_e = 18 GeV, p_ion = top, L = 10 fb⁻¹/u, P_zz = 0.80.
    Top panel: dN/dφ vs φ for three Δ/F1 scales plus flat reference.
    Bottom panel: A(φ) = P_zz · ⟨A_cos2φ⟩ · cos(2φ) for each scale, with a
    shared gray ±1σ band from err_cos2phi_amplitude(N_total, P_zz).
    """
    tag = backends["tag"]
    base = backends["base"]

    cfg = BeamConfig(electron_energy=EE_FIG5, ion=LI6,
                     ion_momentum_per_nucleon=PION_FIG5)
    sc = fom.Scenario(lumi_fb_per_nucleon=LUMI_REF, pol_ion_tensor=PZZ_FIG5)
    nf2_obj = NuclearF2(LI6, base=base)
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_obj)

    nf2 = proj.extras["nf2"]
    f1 = nf2.f1a(proj.x, proj.q2) / LI6.A
    f2 = nf2.f2a(proj.x, proj.q2) / LI6.A
    y = proj.extras["y"]
    n_events = proj.n_events
    accepted = proj.accepted

    # Total flat yield (φ-averaged reference): N_total / (2π)
    n_total = np.where(accepted, n_events, 0.0).sum()
    flat_density = n_total / (2.0 * np.pi)

    # Statistical uncertainty on the cos(2φ) amplitude
    sigma_a_amp = err_cos2phi_amplitude(n_total, PZZ_FIG5)

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True, figsize=(8, 8))

    cos2phi = np.cos(2.0 * PHI_VALS)

    for scale, color in zip(SCALES, SCALE_COLORS):
        delta = toy_delta_gluon(proj.x, proj.q2, f1, scale=scale)
        amp = a_cos2phi(delta, f1, f2, proj.x, y)

        # ── top panel: dN/dφ ────────────────────────────────────────────────
        phi_yield = np.array([
            np.where(
                accepted,
                n_events / (2.0 * np.pi) * (1.0 + PZZ_FIG5 * amp * np.cos(2.0 * phi)),
                0.0,
            ).sum()
            for phi in PHI_VALS
        ])
        label = fr"$\Delta/F_1 = {scale:.0e}$"
        ax_top.plot(PHI_VALS, phi_yield, color=color, lw=1.5, label=label)

        # ── bottom panel: A(φ) ──────────────────────────────────────────────
        # yield-weighted mean amplitude over accepted bins
        n_acc = np.where(accepted, n_events, 0.0)
        n_sum = n_acc.sum()
        mean_amp = (n_acc * amp).sum() / np.maximum(n_sum, 1e-30)
        a_phi = PZZ_FIG5 * mean_amp * cos2phi
        ax_bot.plot(PHI_VALS, a_phi, color=color, lw=1.5, label=label)

    # Flat reference (top panel only)
    ax_top.axhline(flat_density, color="gray", lw=1.2, ls="--",
                   label=r"flat $N_\mathrm{total}/(2\pi)$")

    # ±1σ band on bottom panel — same for all scales, draw once in gray
    band_hi = sigma_a_amp * np.abs(cos2phi)
    ax_bot.fill_between(PHI_VALS, -band_hi, band_hi,
                        alpha=0.20, color="gray",
                        label=r"$\pm1\sigma$ ($L=10$ fb$^{-1}$/u)")
    ax_bot.axhline(0.0, color="black", lw=0.8, ls="--")

    # x-tick labels at multiples of π/2
    xticks = [0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi]
    xticklabels = ["$0$", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"]
    ax_bot.set_xticks(xticks)
    ax_bot.set_xticklabels(xticklabels)
    ax_bot.set_xlabel(r"$\phi$ [rad]", fontsize=11)

    ax_top.set_ylabel(r"$dN/d\phi$ [events rad$^{-1}$]", fontsize=11)
    ax_bot.set_ylabel(r"$A(\phi) = P_{zz}\langle A_{\cos 2\phi}\rangle\cos(2\phi)$",
                      fontsize=10)

    ax_top.legend(fontsize=9)
    ax_bot.legend(fontsize=9)

    fig.suptitle(
        rf"Expected $^6$Li DIS yield and asymmetry vs $\phi$, "
        rf"$\mathcal{{L}}=10$ fb$^{{-1}}$/u, $E_e={EE_FIG5:g}$ GeV, "
        rf"$p_{{ion}}=$ top ({tag.upper()} inputs)",
        fontsize=10,
    )
    fig.tight_layout()

    path = outdir / "money_delta_20260715_yield_vs_phi.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"wrote {path}")
    return path


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Gluonometry reach vs beam energy for 6Li at EIC."
    )
    ap.add_argument("--pdf", default="toy", choices=["toy", "grid"])
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    backends = get_backends(args.pdf)
    tag = backends["tag"]

    # ── Figure 1: vs E_e ──────────────────────────────────────────────────────
    print(f"Building Figure 1 (vs E_e, {len(EE_SWEEP)} points × "
          f"{len(SCALES)} scales) …")
    fig1, rows1 = build_fig1(backends)
    path1 = outdir / "money_delta_20260715_vs_Ee.png"
    fig1.savefig(path1, dpi=150)
    plt.close(fig1)
    print(f"wrote {path1}")

    # ── Figure 2: vs p_ion ────────────────────────────────────────────────────
    print(f"Building Figure 2 (vs p_ion, {len(PION_SWEEP)} points × "
          f"{len(SCALES)} scales) …")
    fig2, rows2 = build_fig2(backends)
    path2 = outdir / "money_delta_20260715_vs_Eion.png"
    fig2.savefig(path2, dpi=150)
    plt.close(fig2)
    print(f"wrote {path2}")

    # ── Figure 3: L_5σ vs √s ─────────────────────────────────────────────────
    print(f"Building Figure 3 (L_5σ vs √s, {len(EE_CURVES)} E_e × "
          f"{len(SCALES)} scales × 14 p_ion points) …")
    path3 = plot_l5sig_vs_sqrts(backends, outdir)

    # ── Figure 4: (x, Q²) coverage ────────────────────────────────────────────
    print(f"Building Figure 4 ((x, Q²) coverage, {len(FIG4_CONFIGS)} panels) …")
    path4 = plot_xq2_coverage(backends, outdir)

    # ── Figure 5: yield vs φ ──────────────────────────────────────────────────
    print(f"Building Figure 5 (yield vs φ, {len(SCALES)} Delta/F1 scales, "
          f"{len(PHI_VALS)} phi points) …")
    path5 = plot_yield_vs_phi(backends, outdir)

    # ── summary table ─────────────────────────────────────────────────────────
    p_ion_top = round(PION_TOP, 1)
    print()
    print(f"Summary ({tag.upper()} inputs, P_zz={PZZ:g})")
    print("-" * 72)
    print(f"{'Figure':<10}  {'scale':>10}  {'sweep endpoint':>20}  "
          f"{'L_5sig [fb^-1/u]':>18}")
    print("-" * 72)
    for scale, ee_end, l5 in rows1:
        print(f"{'vs Ee':<10}  {scale:>10.0e}  "
              f"{'Ee=' + str(ee_end) + ' GeV':>20}  {l5:>18.2f}")
    for scale, pion_end, l5 in rows2:
        print(f"{'vs Eion':<10}  {scale:>10.0e}  "
              f"{'p=' + f'{pion_end:.1f}' + ' GeV/u':>20}  {l5:>18.2f}")
    print("-" * 72)
    print(f"figure 3 (L_5σ vs √s): {path3}")
    print(f"figure 4: {path4}")
    print(f"figure 5: {path5}")


if __name__ == "__main__":
    main()
