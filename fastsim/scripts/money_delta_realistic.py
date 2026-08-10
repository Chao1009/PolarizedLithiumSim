#!/usr/bin/env python3
"""Realistic 6Li gluonometry discovery reach: L_5σ vs Δ/F₁ scale.

Realism ingredients (all active):
──────────────────────────────────
R1. Grid backend (CT18NLO via the `parton` package).
    Scenario(q2_min=2.0) throughout to stay inside the CT18NLO fit range.
    Code: fastsim/polli_fastsim/inputs.py, structure.py

R2. Three Δ(x,Q²) shape variants, normalized so peak Δ/F₁ equals the
    `scale` parameter: low_x (α=0.3, β=4), mid_x (α=0.7, β=3),
    high_x (α=1.5, β=2). The spread gives a shape-sensitivity band.
    The toy_delta_gluon shape (α=0.3, β=4) is the same as "low_x".

R3. EMC ratio hook: NuclearF2 constructed with emc_ratio=unpolarized_emc_ratio.
    Code: fastsim/polli_fastsim/polarized.py  (unpolarized_emc_ratio)
    Reference: EPPS21-era qualitative shape; replace before publication.

R4. Two R = σ_L/σ_T parameterizations:
    - R1998: simplified form of Abe et al., PLB 452:194 (1999).
    - Christy-Bosted: simplified DIS-region form from PRC 81:055213 (2010).
    Both override the module-level r_sigma_lt in structure.py and asymmetries.py
    via context-manager monkey-patching.

[R5 (η-dependent electron-ID efficiency) was removed at user request 2026-07-16;
 efficiency is now uniform inside the |η| ≤ 3.5 acceptance.]

R6. Two P_zz conventions for the same nominal ring P_zz = 0.8:
    - Cloet-conservative (per-nucleon-normalized): P_zz_eff = 0.267 = 0.8/3
    - Cluster-deuteron picture: P_zz_eff = 0.70
    PROVISIONAL — the factor ~2.4 ambiguity is unresolved; see beams.py
    docstring "6Li: UNRESOLVED convention (factor 2.4!)".

Output:
    fastsim/out/money_delta/money_delta_realistic_mid.png   (10 GeV × 50 GeV/u)
    fastsim/out/money_delta/money_delta_realistic_top.png   (18 GeV × 137.5 GeV/u)
"""

import argparse
import pathlib
import sys
from contextlib import contextmanager

import numpy as np

# ── ensure polli_fastsim is importable ────────────────────────────────────────
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.beams import LI6, BeamConfig
from polli_fastsim.asymmetries import a_cos2phi
from polli_fastsim.inputs import get_backends
from polli_fastsim.polarized import toy_delta_gluon, unpolarized_emc_ratio
from polli_fastsim.structure import NuclearF2
import polli_fastsim.structure as structure_mod
import polli_fastsim.asymmetries as asymmetries_mod

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ══════════════════════════════════════════════════════════════════════════════
# R2: Three Δ shape variants
# ══════════════════════════════════════════════════════════════════════════════

# Peak of x^α (1-x)^β is at x_peak = α/(α+β), peak value = α^α β^β / (α+β)^(α+β).
# We compute peak values analytically for each variant.

def _xab_peak_value(alpha, beta):
    """Peak value of x^alpha * (1-x)^beta at x = alpha/(alpha+beta)."""
    xp = alpha / (alpha + beta)
    return (xp ** alpha) * ((1.0 - xp) ** beta)


_VARIANTS = {
    "low_x":  (0.3, 4.0),   # peak x ≈ 0.070, peak val ≈ 0.4056
    "mid_x":  (0.7, 3.0),   # peak x ≈ 0.189, peak val ≈ 0.1867
    "high_x": (1.5, 2.0),   # peak x ≈ 0.429, peak val ≈ 0.0876
}

# Precompute normalization denominators
_PEAK_VALS = {name: _xab_peak_value(a, b) for name, (a, b) in _VARIANTS.items()}

# Sanity check: at x = x_peak, f1 = 1, Δ should equal `scale`.
def _sanity_check_shapes():
    for name, (alpha, beta) in _VARIANTS.items():
        xp = alpha / (alpha + beta)
        test_scale = 1.0
        f1_test = 1.0
        norm = test_scale / _PEAK_VALS[name]
        delta_at_peak = norm * f1_test * (xp ** alpha) * ((1.0 - xp) ** beta)
        assert abs(delta_at_peak - test_scale) < 1e-10, (
            f"Shape sanity FAILED for variant '{name}': "
            f"Δ(x_peak, f1=1) = {delta_at_peak}, expected {test_scale}"
        )

_sanity_check_shapes()


def delta_shape(x, q2, f1, scale, variant="mid_x"):
    """Δ(x,Q²) normalized so peak Δ/F₁ = scale.

    Parameters
    ----------
    x, q2 : array-like  kinematic variables
    f1    : array-like  unpolarized F₁ per nucleon (same shape)
    scale : float       peak Δ/F₁ value (the x-axis parameter of the money plot)
    variant : str       one of 'low_x', 'mid_x', 'high_x'
    """
    x = np.asarray(x, dtype=float)
    alpha, beta = _VARIANTS[variant]
    norm = scale / _PEAK_VALS[variant]
    return norm * f1 * np.power(np.maximum(x, 1e-12), alpha) * np.power(np.maximum(1.0 - x, 0.0), beta)


# ══════════════════════════════════════════════════════════════════════════════
# R4: Two R = σ_L/σ_T parameterizations
# ══════════════════════════════════════════════════════════════════════════════

def r1998(x, q2):
    """Simplified form of the SLAC R1998 fit (Abe et al., PLB 452:194, 1999).

    Valid for typical DIS kinematics (Q² > 1 GeV²); produces R ~ 0.05-0.25
    in that range.  Clipped to [0, 1] for safety.
    """
    q2 = np.asarray(q2, dtype=float)
    x  = np.asarray(x,  dtype=float)
    log_q2 = np.log(np.maximum(q2, 1.01) / 0.04)
    theta = 1.0 + 12.0 * q2 / (q2 + 1.0) * 0.125**2 / (0.125**2 + x**2)
    part1 = 0.0485 * theta / log_q2
    part2 = 0.5470 * theta / (q2**2 + 2.0**4) ** 0.5
    part3 = 0.2379 * theta / (q2 + 0.3**2) * (
        1.0 - 0.00815 * np.log(np.maximum(x, 1e-12)) / np.log(0.5)
    )
    return np.clip(part1 + part2 + part3, 0.0, 1.0)


def r_christy_bosted(x, q2):
    """Simplified DIS-region form of Christy & Bosted PRC 81:055213 (2010).

    Valid for Q² > 1 GeV²; clipped to [0, 0.5].
    """
    q2  = np.asarray(q2,  dtype=float)
    x   = np.asarray(x,   dtype=float)
    mn2 = 0.9383 ** 2
    xi  = 2.0 * x / (1.0 + np.sqrt(1.0 + 4.0 * x**2 * mn2 / np.maximum(q2, 0.01)))
    r   = 0.32 * (1.0 - xi) ** 2 / (1.0 + q2 / 5.0) + 0.05
    return np.clip(r, 0.0, 0.5)


# Store the original r_sigma_lt before any monkey-patching
_original_r = structure_mod.r_sigma_lt


@contextmanager
def r_override(r_func):
    """Context manager: temporarily override r_sigma_lt in structure and asymmetries."""
    structure_mod.r_sigma_lt = r_func
    asymmetries_mod.r_sigma_lt = r_func
    try:
        yield
    finally:
        structure_mod.r_sigma_lt = _original_r
        asymmetries_mod.r_sigma_lt = _original_r


# ══════════════════════════════════════════════════════════════════════════════
# R6: P_zz dilution conventions (PROVISIONAL)
# Note: simplified representation of the underlying dilution ambiguity.
# The factor ~2.4 uncertainty is unresolved; see beams.py docstring on
# "6Li: UNRESOLVED convention (factor 2.4!)".
# ══════════════════════════════════════════════════════════════════════════════

# Two P_zz dilution scenarios for 6Li at nominal ring P_zz = 0.8.
# The convention ambiguity is unresolved pending nuclear-structure expert input.
#   Cloet 1/3   = per-nucleon averaged: (2 polarized / 6 total) × 0.8 = 0.267.
#   Cluster d, folded = deuteron-cluster picture, occupancy-folded:
#       deuteron-cluster nucleons see P_d ≈ 0.87 × 0.8 = 0.70,
#       folded into uniform per-bin sum: sqrt(2/6) × 0.70 ≈ 0.404.
# See fastsim/notes/money_delta_note_2026-07-16.md §7 for full discussion.
PZZ_SCENARIOS = [
    ("Cloet 1/3",        0.267),   # per-nucleon-normalized: 0.8 / 3
    ("Cluster d, folded", 0.40),   # deuteron-cluster picture, occupancy-folded: sqrt(2/6) × 0.70 ≈ 0.404
]


# ══════════════════════════════════════════════════════════════════════════════
# Core significance helper (realistic)
# ══════════════════════════════════════════════════════════════════════════════

def sig2_per_fb_at_realistic(cfg, scale, pzz, base, variant="mid_x",
                             min_events=10):
    """Sig² per fb⁻¹/nucleon with realistic ingredients.

    Parameters
    ----------
    cfg       : BeamConfig
    scale     : float  peak Δ/F₁ value (reference scale)
    pzz       : float  effective P_zz per nucleon
    base      : PartonF2  CT18NLO base for NuclearF2
    variant   : str    Δ shape variant ('low_x', 'mid_x', 'high_x')
    min_events   : int  minimum events to include a bin

    The R-override must be active in the calling context (R4).
    EMC ratio is applied via NuclearF2 constructor (R3).
    Electron-ID efficiency is uniform inside the |η| ≤ 3.5 acceptance (R5
    placeholder removed 2026-07-16).
    """
    sc = fom.Scenario(lumi_fb_per_nucleon=1.0, pol_ion_tensor=pzz, q2_min=2.0)
    # R3: apply EMC ratio hook
    nf2_in = NuclearF2(cfg.ion, base=base, emc_ratio=unpolarized_emc_ratio)
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_in)

    nf2 = proj.extras["nf2"]
    f1 = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2 = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y  = proj.extras["y"]

    # R2: parametric Δ shape
    delta = delta_shape(proj.x, proj.q2, f1, scale=scale, variant=variant)
    amp   = a_cos2phi(delta, f1, f2, proj.x, y)

    n_events = proj.n_events
    use = proj.accepted & (n_events >= min_events)
    return np.where(use, amp**2 * pzz**2 * n_events / 2.0, 0.0).sum()


def sig2_per_fb_toy(cfg, scale, pzz, base, min_events=10):
    """Toy-backend sig² for comparison (no EMC, no η eff, toy R, toy Δ shape).

    Mirrors money_delta_pdfgrid.py sig2_per_fb_at but uses the original
    (toy) r_sigma_lt and toy_delta_gluon shape.  The `base` here may be
    CT18NLO; we only skip the EMC ratio and η efficiency.
    """
    sc = fom.Scenario(lumi_fb_per_nucleon=1.0, pol_ion_tensor=pzz, q2_min=2.0)
    nf2_in = NuclearF2(cfg.ion, base=base)  # no EMC ratio
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_in)

    nf2 = proj.extras["nf2"]
    f1 = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2 = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y  = proj.extras["y"]

    delta = toy_delta_gluon(proj.x, proj.q2, f1, scale=scale)
    amp   = a_cos2phi(delta, f1, f2, proj.x, y)

    use = proj.accepted & (proj.n_events >= min_events)
    return np.where(use, amp**2 * pzz**2 * proj.n_events / 2.0, 0.0).sum()


# ══════════════════════════════════════════════════════════════════════════════
# Plot builder
# ══════════════════════════════════════════════════════════════════════════════

# x-axis: 15 log-spaced Δ/F₁ scale points
SCALES = np.logspace(-3.3, -1.7, 15)
S0 = 1e-3          # reference scale for the analytic rescaling

R_FUNCS = [
    ("R1998",         r1998),
    ("Christy-Bosted", r_christy_bosted),
]

VARIANT_NAMES = ["low_x", "mid_x", "high_x"]

# "central" combination for the thick solid overlay
CENTRAL_VARIANT = "mid_x"
CENTRAL_R_NAME  = "R1998"
CENTRAL_PZZ     = 0.267   # Cloet 1/3

# toy baseline P_zz (matching money_delta.py convention)
TOY_PZZ = 0.80


def build_realistic_plot(cfg, base, outdir, tag, title_beam):
    """Build one realistic money plot for a given BeamConfig.

    Parameters
    ----------
    cfg        : BeamConfig
    base       : PartonF2  CT18NLO backend
    outdir     : pathlib.Path
    tag        : str  file-name suffix ('mid' or 'top')
    title_beam : str  human-readable beam description for the title

    Returns
    -------
    pathlib.Path : output PNG path
    """

    print(f"\n  Config: {cfg.label()}")

    # ── Step 1: pre-compute sig² at S0 for every (variant, R, Pzz) combo ──
    # Total: 3 variants × 2 R-funcs × 2 Pzz = 12 combinations.
    # Then L_5σ(s) = 25 / (sig2_at_s0 × (s/s0)²) for any scale s.

    combo_sig2 = {}   # key: (variant, r_name, pzz_label) → sig2 at S0

    n_combos = len(VARIANT_NAMES) * len(R_FUNCS) * len(PZZ_SCENARIOS)
    done = 0
    for variant in VARIANT_NAMES:
        for r_name, r_func in R_FUNCS:
            with r_override(r_func):
                for pzz_label, pzz in PZZ_SCENARIOS:
                    done += 1
                    print(f"    [{done}/{n_combos}] variant={variant}, R={r_name}, "
                          f"Pzz={pzz_label}({pzz:.3f}) …", end=" ", flush=True)
                    s2 = sig2_per_fb_at_realistic(
                        cfg, S0, pzz, base,
                        variant=variant,
                    )
                    combo_sig2[(variant, r_name, pzz_label)] = s2
                    l5 = 25.0 / max(s2, 1e-30)
                    print(f"L_5σ(1e-3) = {l5:.2f} fb⁻¹/u")

    # Central curve (mid_x, R1998, Cloet 1/3)
    central_key = (CENTRAL_VARIANT, CENTRAL_R_NAME, "Cloet 1/3")
    sig2_central = combo_sig2[central_key]

    # ── Step 2: toy-backend comparison (original toy R, no EMC, no η eff) ──
    print(f"    [toy] toy_delta_gluon, original R, P_zz={TOY_PZZ} …", end=" ", flush=True)
    # Use original r (not overridden)
    sig2_toy = sig2_per_fb_toy(cfg, S0, TOY_PZZ, base)
    print(f"L_5σ(1e-3) = {25.0 / max(sig2_toy, 1e-30):.2f} fb⁻¹/u")

    # ── Step 3: build reach curves for all 12 combos ──
    all_curves = []   # list of L_5σ arrays (one per combination)
    for key, s2 in combo_sig2.items():
        reach = 25.0 / np.maximum(s2 * (SCALES / S0) ** 2, 1e-30)
        all_curves.append(reach)

    all_curves = np.array(all_curves)  # shape (12, 15)
    band_min = all_curves.min(axis=0)
    band_max = all_curves.max(axis=0)

    reach_central = 25.0 / np.maximum(sig2_central * (SCALES / S0) ** 2, 1e-30)
    reach_toy     = 25.0 / np.maximum(sig2_toy     * (SCALES / S0) ** 2, 1e-30)

    # ── Step 4: draw the plot ──────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))

    # Uncertainty band (min–max over all 12 combinations)
    ax.fill_between(SCALES, band_min, band_max,
                    color="steelblue", alpha=0.22,
                    label="Min-max band (3 shapes × 2 R × 2 $P_{zz}$)")

    # Central realistic curve
    ax.plot(SCALES, reach_central, "-", color="black", lw=2.2,
            label=r"Central realistic: mid-$x$ $\Delta$, R1998, Cloet $P_{zz}=0.267$")

    # Toy-backend comparison
    ax.plot(SCALES, reach_toy, "--", color="gray", lw=1.6,
            label=r"Toy backend, $P_{zz}=0.8$ (previous)")

    # Gold band: plausible EIC program
    ax.axhspan(1, 100, color="gold", alpha=0.12,
               label=r"1$-$100 fb$^{-1}$/u (plausible program)")

    # Sather-Schmidt reference line
    ax.axvline(1e-3, color="dimgray", ls=":", lw=1.2)
    ax.text(1.08e-3, 0.90, "Sather-Schmidt\n$O(10^{-3})$", fontsize=7,
            transform=ax.get_xaxis_transform(), va="top", color="dimgray")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta/F_1$ scale (peak of shape)", fontsize=11)
    ax.set_ylabel(r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]", fontsize=11)
    ax.set_title(
        f"$^6$Li gluonometry realistic reach — {title_beam}\n"
        r"CT18NLO + EMC + [R1998/CB] + $P_{zz}$ dilution band",
        fontsize=10,
    )
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"money_delta_realistic_{tag}.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")

    # ── Step 5: summary numbers at Δ/F₁ = 1e-3 ──────────────────────────
    l5_central = 25.0 / max(sig2_central, 1e-30)
    l5_band_lo = 25.0 / max(combo_sig2[max(combo_sig2, key=lambda k: combo_sig2[k])], 1e-30)
    l5_band_hi = 25.0 / max(combo_sig2[min(combo_sig2, key=lambda k: combo_sig2[k])], 1e-30)
    l5_toy     = 25.0 / max(sig2_toy, 1e-30)

    summary = {
        "tag":        tag,
        "label":      cfg.label(),
        "l5_central": l5_central,
        "l5_band_lo": l5_band_lo,
        "l5_band_hi": l5_band_hi,
        "l5_toy":     l5_toy,
    }
    return outpath, summary


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description="Realistic 6Li gluonometry reach: L_5σ vs Δ/F₁ (two EIC configs)."
    )
    ap.add_argument(
        "--outdir",
        default="fastsim/out/money_delta",
        help="Output directory for PNGs (default: fastsim/out/money_delta)",
    )
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)

    # ── Load CT18NLO backend (fail gracefully if parton not installed) ────
    try:
        backends = get_backends("grid")
    except Exception as exc:
        print(
            "ERROR: grid backend requires `parton` with CT18NLO grid.\n"
            "Install with:\n"
            "  pip install parton\n"
            "  python3 -m parton install CT18NLO\n"
            f"Underlying error: {exc}",
            file=sys.stderr,
        )
        sys.exit(2)
    base = backends["base"]

    # ── Two EIC beam configurations ────────────────────────────────────────
    # Plot A: mid config  (E_e=10 GeV, p_ion=50 GeV/u)
    # Plot B: top config  (E_e=18 GeV, p_ion=137.5 GeV/u)
    configs = [
        (BeamConfig(electron_energy=10.0, ion=LI6, ion_momentum_per_nucleon=50.0),
         "mid", "10 GeV $e$ × 50 GeV/u $^6$Li"),
        (BeamConfig(electron_energy=18.0, ion=LI6, ion_momentum_per_nucleon=137.5),
         "top", "18 GeV $e$ × 137.5 GeV/u $^6$Li"),
    ]

    summaries = []
    for cfg, tag, title_beam in configs:
        print(f"\n{'='*66}")
        print(f"Building realistic plot: {tag.upper()} config — {cfg.label()}")
        print(f"{'='*66}")
        _, summary = build_realistic_plot(cfg, base, outdir, tag, title_beam)
        summaries.append(summary)

    # ── Summary table ─────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("SUMMARY TABLE  (Δ/F₁ = 1e-3, all values in fb⁻¹/nucleon)")
    print("=" * 72)
    hdr = f"{'Config':<30} {'Central':>10} {'Band lo':>10} {'Band hi':>10} {'Toy (P_zz=0.8)':>16}"
    print(hdr)
    print("-" * 72)
    for s in summaries:
        print(
            f"  {s['label']:<28} "
            f"{s['l5_central']:>10.2f} "
            f"{s['l5_band_lo']:>10.2f} "
            f"{s['l5_band_hi']:>10.2f} "
            f"{s['l5_toy']:>16.2f}"
        )
    print("=" * 72)
    print("Note: 'Band lo/hi' = best/worst of 12 combinations at Δ/F₁=1e-3.")
    print("      Central = mid_x shape, R1998, Cloet P_zz=0.267.")
    print("      Toy     = toy_delta_gluon, default R, no EMC, uniform eff, P_zz=0.8.")


if __name__ == "__main__":
    main()
