#!/usr/bin/env python3
"""Grid-backend replica of all six money-Δ figures.

This script reproduces every figure produced by the two parent scripts:
  * fastsim/scripts/money_delta.py        → Figure 1 (scale scan)
  * fastsim/scripts/money_delta_20260715.py → Figures 2–6 (energy sweeps,
    √s reach, (x,Q²) coverage, yield vs φ)

Results differ from the toy-backend versions *only* through the PDF inputs:
  Unpolarized: CT18NLO   (via the `parton` package)
  Polarized:   NNPDFpol11_100 (via the `parton` package)

The physics chain, kinematic cuts, significance formula, and all helper logic
are identical to the parent scripts.

PDF validity note
-----------------
All ``Scenario`` objects in this script are constructed with ``q2_min=2.0``
(instead of the shared default of 1.0 GeV²).  CT18NLO PDFs return NaN for
Q² below roughly 1.7 GeV², so raising the kinematic floor to 2.0 GeV² keeps
every grid evaluation inside the PDF fit range.  This is a physics-motivated
fiducial cut: the affected bins (1 < Q² < 2 GeV²) represent a tiny fraction
of the total DIS cross-section at EIC kinematics and do not change any
physics conclusion.  The toy-backend scripts (money_delta.py and
money_delta_20260715.py) are unaffected and retain the default q2_min=1.0.

Output filenames (all written to --outdir, default "out"):
  money_delta_pdfgrid_scale_scan.png       (Figure 1)
  money_delta_pdfgrid_vs_Ee.png            (Figure 2)
  money_delta_pdfgrid_vs_Eion.png          (Figure 3)
  money_delta_pdfgrid_L5sig_vs_sqrts.png   (Figure 4)
  money_delta_pdfgrid_xQ2_coverage.png     (Figure 5)
  money_delta_pdfgrid_yield_vs_phi.png     (Figure 6)

NOTE: Uses Scenario(q2_min=2.0) throughout to stay inside the CT18NLO PDF fit range (~1.7 GeV^2 lower bound).
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

EE_SWEEP = np.linspace(5.0, 18.0, 14)          # GeV, Figure 2
PION_TOP = LI6.momentum_per_nucleon_max         # ~137.5 GeV/u
PION_SWEEP = np.linspace(20.0, PION_TOP, 14)   # GeV/u, Figure 3
EE_FIXED = 18.0                                 # GeV, Figure 3

# Figure 4 parameters
EE_CURVES = [5.0, 10.0, 18.0]                  # GeV, one colour per E_e
EE_CURVE_COLORS = ("steelblue", "seagreen", "darkorchid")
LUMI_REF = 10.0                                 # fb^-1/nucleon, Figures 4–6

# Figure 5: four representative (E_e, p_ion) panels
FIG5_CONFIGS = [
    (5.0,  41.0 * 0.5,   "low E_e, low p_ion"),
    (5.0,  PION_TOP,     "low E_e, top p_ion"),
    (18.0, 41.0 * 0.5,   "top E_e, low p_ion"),
    (18.0, PION_TOP,     "top E_e, top p_ion"),
]

# Figure 6 parameters
EE_FIG6   = 18.0
PION_FIG6 = PION_TOP
PZZ_FIG6  = 0.80
PHI_VALS  = np.linspace(0.0, 2.0 * np.pi, 200)


# ── physics helper (verbatim from money_delta.py / money_delta_20260715.py) ───

def sig2_per_fb_at(cfg, scale, pzz, base=None, min_events=10):
    """Significance^2 per fb^-1/nucleon at a reference Delta/F1 scale.

    The toy Delta shape is linear in `scale`, so reach curves follow
    L_5sig(s) = 25 / (sig2 * (s/scale)^2) analytically.
    """
    sc = fom.Scenario(lumi_fb_per_nucleon=1.0, pol_ion_tensor=pzz, q2_min=2.0)
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


# ── shared axis helper ────────────────────────────────────────────────────────

def _add_plausible_band(ax):
    """Gold axhspan for the 1–100 fb^-1/u plausible EIC program."""
    ax.axhspan(1, 100, color="gold", alpha=0.12,
               label=r"1$-$100 fb$^{-1}$/u (plausible program)")


# ── figure builders ───────────────────────────────────────────────────────────

def build_fig1_scale_scan(backends):
    """Figure 1: L_5sig vs Delta/F1 scale for canonical 6Li configs (GRID).

    Reproduces money_delta.py exactly for 6Li with the grid backend.
    Three default_configs × two P_zz values (0.60, 0.80) = six curves.
    Uses analytic shortcut: compute sig2 once at s0=1e-3, then rescale.
    Returns (fig, sig2_ref) where sig2_ref[(label, pzz)] = sig2 at s0.
    """
    scales = np.logspace(-3.3, -1.7, 15)
    s0 = 1e-3
    base = backends["base"]

    fig, ax = plt.subplots(figsize=(7, 5))
    sig2_ref = {}

    for cfg, color in zip(beams.default_configs("6Li"),
                          ("crimson", "seagreen", "navy")):
        for pzz, ls in ((0.60, "-"), (0.80, "--")):
            sig2 = sig2_per_fb_at(cfg, s0, pzz, base=base)
            sig2_ref[(cfg.label(), pzz)] = sig2
            reach = 25.0 / np.maximum(sig2 * (scales / s0) ** 2, 1e-30)
            ax.plot(scales, reach, ls, color=color, lw=1.5,
                    label=f"{cfg.label()}, $P_{{zz}}$={pzz:g}")

    ax.axhspan(1, 100, color="gold", alpha=0.12,
               label=r"1$-$100 fb$^{-1}$/u (plausible program)")
    ax.axvline(1e-3, color="gray", ls=":", lw=1)
    ax.text(1.05e-3, 0.93, "Sather-Schmidt\n$O(10^{-3})$", fontsize=7,
            transform=ax.get_xaxis_transform(), va="top")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta/F_1$ scale (peak of scenario shape)")
    ax.set_ylabel(r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]")
    ax.set_title(
        "Nuclear gluonometry reach, transversely polarized 6Li\n"
        r"(cos 2$\phi$ amplitude, all bins combined; GRID inputs)",
        fontsize=10,
    )
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    return fig, sig2_ref


def build_fig2_vs_Ee(backends):
    """Figure 2: L_5sig vs E_e at top ion energy, P_zz=0.80 (GRID)."""
    p_ion = round(PION_TOP, 1)
    base = backends["base"]

    fig, ax = plt.subplots(figsize=(7, 5))

    for scale, color in zip(SCALES, SCALE_COLORS):
        l5sig_vals = []
        for ee in EE_SWEEP:
            cfg = BeamConfig(electron_energy=ee, ion=LI6,
                             ion_momentum_per_nucleon=p_ion)
            s2 = sig2_per_fb_at(cfg, scale, PZZ, base=base)
            l5sig_vals.append(25.0 / np.maximum(s2, 1e-30))
        l5sig_vals = np.array(l5sig_vals)
        label = fr"$\Delta/F_1 = {scale:.0e}$"
        ax.plot(EE_SWEEP, l5sig_vals, "-o", color=color, lw=1.5,
                ms=4, label=label)

    _add_plausible_band(ax)
    ax.set_yscale("log")
    ax.set_xlabel(r"$E_e$ [GeV]", fontsize=11)
    ax.set_ylabel(r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]", fontsize=11)
    ax.set_title(
        rf"$^6$Li, top ion energy ({p_ion:.1f} GeV/u), $P_{{zz}}={PZZ:g}$"
        r" (cos 2$\phi$ amplitude, all bins combined; GRID inputs)",
        fontsize=10,
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def build_fig3_vs_Eion(backends):
    """Figure 3: L_5sig vs p_ion at E_e=18 GeV, P_zz=0.80 (GRID)."""
    base = backends["base"]

    fig, ax = plt.subplots(figsize=(7, 5))

    for scale, color in zip(SCALES, SCALE_COLORS):
        l5sig_vals = []
        for p_ion in PION_SWEEP:
            cfg = BeamConfig(electron_energy=EE_FIXED, ion=LI6,
                             ion_momentum_per_nucleon=p_ion)
            s2 = sig2_per_fb_at(cfg, scale, PZZ, base=base)
            l5sig_vals.append(25.0 / np.maximum(s2, 1e-30))
        l5sig_vals = np.array(l5sig_vals)
        label = fr"$\Delta/F_1 = {scale:.0e}$"
        ax.plot(PION_SWEEP, l5sig_vals, "-o", color=color, lw=1.5,
                ms=4, label=label)

    _add_plausible_band(ax)
    ax.set_yscale("log")
    ax.set_xlabel(r"$p_\mathrm{ion}$ [GeV/u]", fontsize=11)
    ax.set_ylabel(r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]", fontsize=11)
    ax.set_title(
        rf"$^6$Li, $E_e = {EE_FIXED:g}$ GeV, $P_{{zz}}={PZZ:g}$"
        r" (cos 2$\phi$ amplitude, all bins combined; GRID inputs)",
        fontsize=10,
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def build_fig4_l5sig_vs_sqrts(backends):
    """Figure 4: L_5sig vs sqrt(s) per nucleon, nine curves (GRID).

    E_e ∈ {5, 10, 18} GeV × Δ/F1 scale ∈ {1e-3, 3e-3, 1e-2} = 9 curves.
    Color encodes E_e; linestyle encodes scale.
    """
    base = backends["base"]
    pzz = 0.80
    scale_ls = {1e-3: "-", 3e-3: "--", 1e-2: ":"}
    pion_sweep = np.linspace(20.0, LI6.momentum_per_nucleon_max, 14)

    fig, ax = plt.subplots(figsize=(7, 5))

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
            label = rf"$E_e={ee:g}$ GeV, $\Delta/F_1={scale:.0e}$"
            ax.plot(np.array(sqrts_vals), np.array(l5sig_vals),
                    ls, color=color, lw=1.5, ms=4, label=label)

    ax.axhspan(1, 100, color="gold", alpha=0.12,
               label=r"1$-$100 fb$^{-1}$/u (plausible program)")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\sqrt{s}$ per nucleon [GeV]", fontsize=11)
    ax.set_ylabel(r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]", fontsize=11)
    ax.set_title(
        rf"$L_{{5\sigma}}$ vs $\sqrt{{s}}$, $^6$Li, $P_{{zz}}={pzz:g}$"
        " (GRID inputs)",
        fontsize=10,
    )
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    return fig


def build_fig5_xq2_coverage(backends):
    """Figure 5: accepted (x, Q²) coverage as pcolormesh in 2×2 panels (GRID)."""
    base = backends["base"]
    sc = fom.Scenario(lumi_fb_per_nucleon=LUMI_REF, q2_min=2.0)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes_flat = axes.flatten()

    for ax, (ee, p_ion, _desc) in zip(axes_flat, FIG5_CONFIGS):
        cfg = BeamConfig(electron_energy=ee, ion=LI6,
                         ion_momentum_per_nucleon=p_ion)
        nf2 = NuclearF2(LI6, base=base)
        proj = fom.project_rates(cfg, sc, nuclear_f2=nf2)

        z = np.where(proj.accepted, proj.n_events, np.nan)
        z_finite = z[np.isfinite(z) & (z > 0)]
        if z_finite.size == 0:
            ax.set_title(
                f"$E_e$={ee:g}, $p_{{ion}}$={p_ion:g} (no accepted bins)",
                fontsize=8,
            )
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
        "  (GRID inputs)",
        fontsize=11,
    )
    fig.tight_layout()
    return fig


def build_fig6_yield_vs_phi(backends):
    """Figure 6: two-panel dN/dφ (top) and A(φ) with ±1σ band (bottom) (GRID).

    Fixed: E_e = 18 GeV, p_ion = top rigidity limit, L = 10 fb⁻¹/u, P_zz = 0.80.
    Three Δ/F1 scale curves plus flat reference (top panel).
    Yield-weighted mean amplitude drives the bottom-panel asymmetry curve.
    """
    base = backends["base"]

    cfg = BeamConfig(electron_energy=EE_FIG6, ion=LI6,
                     ion_momentum_per_nucleon=PION_FIG6)
    sc = fom.Scenario(lumi_fb_per_nucleon=LUMI_REF, pol_ion_tensor=PZZ_FIG6, q2_min=2.0)
    nf2_obj = NuclearF2(LI6, base=base)
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_obj)

    nf2 = proj.extras["nf2"]
    f1 = nf2.f1a(proj.x, proj.q2) / LI6.A
    f2 = nf2.f2a(proj.x, proj.q2) / LI6.A
    y = proj.extras["y"]
    n_events = proj.n_events
    accepted = proj.accepted

    n_total = np.where(accepted, n_events, 0.0).sum()
    flat_density = n_total / (2.0 * np.pi)
    sigma_a_amp = err_cos2phi_amplitude(n_total, PZZ_FIG6)

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True, figsize=(8, 8))
    cos2phi = np.cos(2.0 * PHI_VALS)

    for scale, color in zip(SCALES, SCALE_COLORS):
        delta = toy_delta_gluon(proj.x, proj.q2, f1, scale=scale)
        amp = a_cos2phi(delta, f1, f2, proj.x, y)

        # top panel: dN/dφ
        phi_yield = np.array([
            np.where(
                accepted,
                n_events / (2.0 * np.pi) * (1.0 + PZZ_FIG6 * amp * np.cos(2.0 * phi)),
                0.0,
            ).sum()
            for phi in PHI_VALS
        ])
        label = fr"$\Delta/F_1 = {scale:.0e}$"
        ax_top.plot(PHI_VALS, phi_yield, color=color, lw=1.5, label=label)

        # bottom panel: A(φ) = P_zz · <A_cos2φ> · cos(2φ)
        n_acc = np.where(accepted, n_events, 0.0)
        n_sum = n_acc.sum()
        mean_amp = (n_acc * amp).sum() / np.maximum(n_sum, 1e-30)
        a_phi = PZZ_FIG6 * mean_amp * cos2phi
        ax_bot.plot(PHI_VALS, a_phi, color=color, lw=1.5, label=label)

    # flat reference (top panel)
    ax_top.axhline(flat_density, color="gray", lw=1.2, ls="--",
                   label=r"flat $N_\mathrm{total}/(2\pi)$")

    # ±1σ band (bottom panel)
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
    ax_bot.set_ylabel(
        r"$A(\phi) = P_{zz}\langle A_{\cos 2\phi}\rangle\cos(2\phi)$",
        fontsize=10,
    )

    ax_top.legend(fontsize=9)
    ax_bot.legend(fontsize=9)

    fig.suptitle(
        rf"Expected $^6$Li DIS yield and asymmetry vs $\phi$, "
        rf"$\mathcal{{L}}=10$ fb$^{{-1}}$/u, $E_e={EE_FIG6:g}$ GeV, "
        rf"$p_{{ion}}=$ top (GRID inputs)",
        fontsize=10,
    )
    fig.tight_layout()
    return fig


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Reproduce all six money-Δ figures with the GRID PDF backend "
            "(CT18NLO + NNPDFpol11_100 via the `parton` package)."
        )
    )
    ap.add_argument("--outdir", default="out",
                    help="Directory for output PNGs (default: out)")
    args = ap.parse_args()

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Force grid backend; fail informatively if parton grids are absent.
    try:
        backends = get_backends("grid")
    except Exception as exc:
        print(
            "ERROR: --pdf grid requires the `parton` package with CT18NLO and "
            "NNPDFpol11_100 grids installed.\n"
            "Install with:\n"
            "  pip install parton\n"
            "  python3 -m parton install CT18NLO\n"
            "  python3 -m parton install NNPDFpol11_100\n"
            f"Underlying error: {exc}",
            file=sys.stderr,
        )
        sys.exit(2)

    # ── Figure 1: scale scan (money_delta.py replica, 6Li, grid) ─────────────
    print(
        f"Building Figure 1 (scale scan, {len(beams.default_configs('6Li'))} configs "
        "× 2 P_zz × 15 scale points) …"
    )
    fig1, sig2_ref = build_fig1_scale_scan(backends)
    path1 = outdir / "money_delta_pdfgrid_scale_scan.png"
    fig1.savefig(path1, dpi=150)
    plt.close(fig1)
    print(f"wrote {path1}")

    # ── Figure 2: L_5sig vs E_e ──────────────────────────────────────────────
    print(
        f"Building Figure 2 (vs E_e, {len(EE_SWEEP)} points "
        f"× {len(SCALES)} scales) …"
    )
    fig2 = build_fig2_vs_Ee(backends)
    path2 = outdir / "money_delta_pdfgrid_vs_Ee.png"
    fig2.savefig(path2, dpi=150)
    plt.close(fig2)
    print(f"wrote {path2}")

    # ── Figure 3: L_5sig vs p_ion ─────────────────────────────────────────────
    print(
        f"Building Figure 3 (vs p_ion, {len(PION_SWEEP)} points "
        f"× {len(SCALES)} scales) …"
    )
    fig3 = build_fig3_vs_Eion(backends)
    path3 = outdir / "money_delta_pdfgrid_vs_Eion.png"
    fig3.savefig(path3, dpi=150)
    plt.close(fig3)
    print(f"wrote {path3}")

    # ── Figure 4: L_5sig vs sqrt(s) ──────────────────────────────────────────
    print(
        f"Building Figure 4 (L_5σ vs √s, {len(EE_CURVES)} E_e "
        f"× {len(SCALES)} scales × 14 p_ion points) …"
    )
    fig4 = build_fig4_l5sig_vs_sqrts(backends)
    path4 = outdir / "money_delta_pdfgrid_L5sig_vs_sqrts.png"
    fig4.savefig(path4, dpi=150)
    plt.close(fig4)
    print(f"wrote {path4}")

    # ── Figure 5: (x, Q²) coverage ───────────────────────────────────────────
    print(f"Building Figure 5 ((x, Q²) coverage, {len(FIG5_CONFIGS)} panels) …")
    fig5 = build_fig5_xq2_coverage(backends)
    path5 = outdir / "money_delta_pdfgrid_xQ2_coverage.png"
    fig5.savefig(path5, dpi=150)
    plt.close(fig5)
    print(f"wrote {path5}")

    # ── Figure 6: yield vs φ ──────────────────────────────────────────────────
    print(
        f"Building Figure 6 (yield vs φ, {len(SCALES)} Δ/F1 scales, "
        f"{len(PHI_VALS)} φ points) …"
    )
    fig6 = build_fig6_yield_vs_phi(backends)
    path6 = outdir / "money_delta_pdfgrid_yield_vs_phi.png"
    fig6.savefig(path6, dpi=150)
    plt.close(fig6)
    print(f"wrote {path6}")

    # ── summary table (Figure 1 analog of money_delta.py lines 86-89) ────────
    s0 = 1e-3
    print()
    print("Summary (GRID inputs, P_zz=0.8, Figure 1 scale-scan configs)")
    print("-" * 72)
    for cfg in beams.default_configs("6Li"):
        sig2 = sig2_ref[(cfg.label(), 0.8)]
        print(
            f"  {cfg.label():26s} L_5sig(Δ/F1=1e-3, Pzz=0.8) = "
            f"{25.0 / sig2:9.1f} fb^-1/u ; (1e-2) = {25.0 / sig2 / 100:7.3f}"
        )
    print("-" * 72)


if __name__ == "__main__":
    main()
