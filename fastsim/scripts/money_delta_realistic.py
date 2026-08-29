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

R3. EMC ratio hook: NuclearF2 constructed with emc_ratio=legacy_emc_ratio,
    i.e. unpolarized_emc_ratio frozen at mode='table' -- the 12-point shape
    this note was written against, kept when the default moved to the
    nuclear PDF grids on 2026-08-29.
    Code: fastsim/polli_fastsim/polarized.py  (legacy_emc_ratio)
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

FROZEN BEAM MENU.  The three configurations below are 5 × 27.5, 10 × 50 and
18 × 137.5 GeV/u — the rigidity-scaled menu this programme used before
2026-08-27.  Ions are γ-matched, not rigidity-scaled (plans/10 A0), so only
the TOP one is a machine configuration; the ⁶Li menu is 40.8 / 99.5 / 137.5.
They are kept because the L_5σ numbers of the 2026-07 notes are quoted
against them and this script is their reproduction — `money_delta.py` is the
same figure of merit on the corrected menu.  fastsim/tests/test_beams.py
exempts this file from the stale-energy sweep for that reason, by name.
"""

import argparse
import importlib.util
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
from polli_fastsim.polarized import legacy_emc_ratio, toy_delta_gluon
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


def r1998_theta_log(x, q2):
    """`r1998` above with Θ restricted to the log term — the minimal repair.

    Abe et al. (E143), PLB 452 (1999) 194 (arXiv:hep-ex/9808028): Θ(x, Q²)
    of their Eq. (3) multiplies the a₁/ln(Q²/0.04) term of Eq. (2) and
    nothing else.  `r1998` above multiplies it into all three terms, which
    saturates the sum and clips R to 1.000 for x ≲ 0.1, Q² ≲ 5 GeV² — the
    region that carries 62–72% of Σa²P²N/2 here, so the L_5σ this script
    publishes is too pessimistic by 1.96× (mid, 131.26 → 67.11 fb⁻¹/u) and
    1.75× (top, 274.64 → 156.91), measured 2026-08-26 with `--r-model
    theta-log` (code review 2026-08-25 item S1).
    Same coefficients and clip as `r1998`, so the difference between the
    two isolates the defect from any change of parameterisation.
    """
    q2 = np.asarray(q2, dtype=float)
    x  = np.asarray(x,  dtype=float)
    log_q2 = np.log(np.maximum(q2, 1.01) / 0.04)
    theta = 1.0 + 12.0 * q2 / (q2 + 1.0) * 0.125**2 / (0.125**2 + x**2)
    part1 = 0.0485 * theta / log_q2
    part2 = 0.5470 / (q2**2 + 2.0**4) ** 0.5
    part3 = 0.2379 / (q2 + 0.3**2) * (
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

# The reach solver is money_delta.py's, imported so that the two cannot
# drift: the min-events floor belongs at the luminosity the reach is quoted
# at, not at the 1 fb^-1/u the sig^2 is normalised to (2026-08-28).  It
# leaves every published number of this script unchanged -- the bins it
# recovers, holding 0.08 to 10 events per fb^-1/u, carry a part in 10^13 of
# sig^2 -- which is why the frozen July reproduction still reproduces.
def _reach_solver():
    path = pathlib.Path(__file__).resolve().parent / "money_delta.py"
    spec = importlib.util.spec_from_file_location("_money_delta", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.reach_from_terms


reach_from_terms = _reach_solver()


def bin_terms_realistic(cfg, scale, pzz, base, variant="mid_x",
                        run_share=1.0):
    """Per-bin (sig² contribution, event count) at 1 fb⁻¹/nucleon, with the
    realistic ingredients.

    Parameters
    ----------
    cfg       : BeamConfig
    scale     : float  peak Δ/F₁ value (reference scale)
    pzz       : float  effective P_zz per nucleon
    base      : PartonF2  CT18NLO base for NuclearF2
    variant   : str    Δ shape variant ('low_x', 'mid_x', 'high_x')

    The R-override must be active in the calling context (R4).
    EMC ratio is applied via NuclearF2 constructor (R3).
    Electron-ID efficiency is uniform inside the |η| ≤ 3.5 acceptance (R5
    placeholder removed 2026-07-16).  Bins whose CT18NLO evaluation is
    non-finite are dropped explicitly (they used to fall out of the
    min-events comparison by accident).
    """
    sc = fom.Scenario(lumi_fb_per_nucleon=1.0, run_share=run_share,
                      pol_ion_tensor=pzz, q2_min=2.0)
    # R3: apply EMC ratio hook
    nf2_in = NuclearF2(cfg.ion, base=base, emc_ratio=legacy_emc_ratio)
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_in)

    nf2 = proj.extras["nf2"]
    f1 = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2 = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y  = proj.extras["y"]

    # R2: parametric Δ shape
    delta = delta_shape(proj.x, proj.q2, f1, scale=scale, variant=variant)
    amp   = a_cos2phi(delta, f1, f2, proj.x, y)

    terms = amp**2 * pzz**2 * proj.n_events / 2.0
    ok = proj.accepted & np.isfinite(terms) & np.isfinite(proj.n_events)
    return terms[ok].ravel(), proj.n_events[ok].ravel()


def sig2_per_fb_at_realistic(cfg, scale, pzz, base, variant="mid_x",
                             min_events=10, lumi_fb=None, run_share=1.0):
    """Sig² per fb⁻¹/nucleon with realistic ingredients, with the
    min-events floor applied at `lumi_fb` (None = the self-consistent
    reach)."""
    terms, n_events = bin_terms_realistic(cfg, scale, pzz, base,
                                          variant=variant,
                                          run_share=run_share)
    if lumi_fb is None:
        lumi_fb = reach_from_terms(terms, n_events, min_events=min_events)
    return float(terms[n_events * lumi_fb >= min_events].sum())


def bin_terms_toy(cfg, scale, pzz, base, run_share=1.0):
    """Toy-backend per-bin terms for comparison (no EMC, no η eff, toy R,
    toy Δ shape).

    Mirrors money_delta.py bin_terms but uses the original (toy)
    r_sigma_lt and toy_delta_gluon shape.  The `base` here may be
    CT18NLO; we only skip the EMC ratio and η efficiency.
    """
    sc = fom.Scenario(lumi_fb_per_nucleon=1.0, run_share=run_share,
                      pol_ion_tensor=pzz, q2_min=2.0)
    nf2_in = NuclearF2(cfg.ion, base=base)  # no EMC ratio
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_in)

    nf2 = proj.extras["nf2"]
    f1 = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2 = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y  = proj.extras["y"]

    delta = toy_delta_gluon(proj.x, proj.q2, f1, scale=scale)
    amp   = a_cos2phi(delta, f1, f2, proj.x, y)

    terms = amp**2 * pzz**2 * proj.n_events / 2.0
    ok = proj.accepted & np.isfinite(terms) & np.isfinite(proj.n_events)
    return terms[ok].ravel(), proj.n_events[ok].ravel()


def sig2_per_fb_toy(cfg, scale, pzz, base, min_events=10, lumi_fb=None,
                    run_share=1.0):
    """Toy-backend sig² for comparison, floor at `lumi_fb`."""
    terms, n_events = bin_terms_toy(cfg, scale, pzz, base,
                                    run_share=run_share)
    if lumi_fb is None:
        lumi_fb = reach_from_terms(terms, n_events, min_events=min_events)
    return float(terms[n_events * lumi_fb >= min_events].sum())


# ══════════════════════════════════════════════════════════════════════════════
# Plot builder
# ══════════════════════════════════════════════════════════════════════════════

# x-axis: 15 log-spaced Δ/F₁ scale points
SCALES = np.logspace(-3.3, -1.7, 15)
S0 = 1e-3          # reference scale for the analytic rescaling
# "one EIC year" in the Yellow Report accounting; the programme luminosity
# a run-plan share divides (plans/07 WP2)
L_PROGRAMME_FB = 10.0

R_FUNCS = [
    ("R1998",         r1998),
    ("Christy-Bosted", r_christy_bosted),
]

# --r-model: which function stands behind the "R1998" entry of R_FUNCS — the
# one that also carries CENTRAL_R_NAME and therefore the published central
# L_5σ.  "simplified" is the frozen July form and the DEFAULT, so the numbers
# money_delta_note_2026-07-16/-20.md record (131.26 / 274.64 fb⁻¹/u) come back
# unchanged; the other two exist because code review 2026-08-25
# recommendation 0(b) asks for this reach to be re-derived with a corrected R.
# Only the function is swapped: the KEY stays "R1998" so the curve, the
# 12-combination band and the central-key lookup are untouched.  What does move
# is the second element, the caption the figure carries: the PNG leaves the
# output directory that names the model, so a non-default run has to say on the
# plot which R produced it, while the default stays the bare "R1998" of the
# published figure (a default run is byte-identical, plots included).
R_MODELS = {
    "simplified": (r1998,               "R1998"),
    "theta-log":  (r1998_theta_log,     r"R1998 [$\Theta$ on log term only]"),
    "r1998":      (structure_mod.r1998, "R1998 [published fit]"),
}

DEFAULT_R_MODEL = "simplified"

VARIANT_NAMES = ["low_x", "mid_x", "high_x"]

# "central" combination for the thick solid overlay
CENTRAL_VARIANT = "mid_x"
CENTRAL_R_NAME  = "R1998"
CENTRAL_PZZ     = 0.267   # Cloet 1/3

# toy baseline P_zz (matching money_delta.py convention)
TOY_PZZ = 0.80


def build_realistic_plot(cfg, base, outdir, tag, title_beam, r_funcs=None,
                         r_caption="R1998", run_share=1.0):
    """Build one realistic money plot for a given BeamConfig.

    Parameters
    ----------
    cfg        : BeamConfig
    base       : PartonF2  CT18NLO backend
    outdir     : pathlib.Path
    tag        : str  file-name suffix ('mid' or 'top')
    title_beam : str  human-readable beam description for the title
    r_funcs    : list  (name, callable) R models; None = the frozen R_FUNCS
    r_caption  : str   R = σ_L/σ_T caption for the central-curve legend and the
                       title; the default is the frozen July caption, so an
                       unpatched caller reproduces the published figure
    run_share  : float this observable's share of the programme year
                       (plans/07 WP2).  Every L_5σ below is the luminosity
                       the measurement must ACCUMULATE and is exactly
                       invariant under it; the share divides the programme
                       luminosity, so the years to that reach scale as
                       1/share.

    Returns
    -------
    pathlib.Path : output PNG path
    """

    print(f"\n  Config: {cfg.label()}")
    print(f"  {fom.run_share_header(L_PROGRAMME_FB, run_share)}")

    def _reach(terms, n_ev, **kw):
        """L_5σ in fb⁻¹/u DELIVERED to this observable: share-invariant
        (see money_delta.reach_fb)."""
        return run_share * reach_from_terms(terms, n_ev, **kw)

    # ── Step 1: pre-compute sig² at S0 for every (variant, R, Pzz) combo ──
    # Total: 3 variants × 2 R-funcs × 2 Pzz = 12 combinations.
    # Then L_5σ(s) = 25 / (sig2_at_s0 × (s/s0)²) for any scale s.

    combo_terms = {}  # key: (variant, r_name, pzz_label) → per-bin terms at S0
    combo_l5 = {}     # the same keys → L_5σ at S0

    r_funcs = R_FUNCS if r_funcs is None else r_funcs
    n_combos = len(VARIANT_NAMES) * len(r_funcs) * len(PZZ_SCENARIOS)
    done = 0
    for variant in VARIANT_NAMES:
        for r_name, r_func in r_funcs:
            with r_override(r_func):
                for pzz_label, pzz in PZZ_SCENARIOS:
                    done += 1
                    print(f"    [{done}/{n_combos}] variant={variant}, R={r_name}, "
                          f"Pzz={pzz_label}({pzz:.3f}) …", end=" ", flush=True)
                    terms, n_ev = bin_terms_realistic(
                        cfg, S0, pzz, base, variant=variant,
                        run_share=run_share)
                    combo_terms[(variant, r_name, pzz_label)] = (terms, n_ev)
                    l5 = _reach(terms, n_ev)
                    combo_l5[(variant, r_name, pzz_label)] = l5
                    print(f"L_5σ(1e-3) = {l5:.2f} fb⁻¹/u")

    # Central curve (mid_x, R1998, Cloet 1/3)
    central_key = (CENTRAL_VARIANT, CENTRAL_R_NAME, "Cloet 1/3")

    # ── Step 2: toy-backend comparison (original toy R, no EMC, no η eff) ──
    print(f"    [toy] toy_delta_gluon, original R, P_zz={TOY_PZZ} …", end=" ", flush=True)
    # Use original r (not overridden)
    terms_toy, n_ev_toy = bin_terms_toy(cfg, S0, TOY_PZZ, base,
                                        run_share=run_share)
    print(f"L_5σ(1e-3) = {_reach(terms_toy, n_ev_toy):.2f} fb⁻¹/u")

    # ── Step 3: build reach curves for all 12 combos ──
    # the amplitude is linear in the Δ/F₁ scale, so the per-bin terms scale
    # as (s/S0)² and only the min-events floor is re-solved per point
    def _curve(terms, n_ev):
        return np.array([_reach(terms * (s / S0) ** 2, n_ev)
                         for s in SCALES])

    all_curves = np.array([_curve(t, n) for t, n in combo_terms.values()])
    band_min = all_curves.min(axis=0)
    band_max = all_curves.max(axis=0)

    reach_central = _curve(*combo_terms[central_key])
    reach_toy     = _curve(terms_toy, n_ev_toy)

    # ── Step 4: draw the plot ──────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))

    # Uncertainty band (min–max over all 12 combinations)
    ax.fill_between(SCALES, band_min, band_max,
                    color="steelblue", alpha=0.22,
                    label="Min-max band (3 shapes × 2 R × 2 $P_{zz}$)")

    # Central realistic curve
    # r_caption, not a literal "R1998": under a non-default --r-model the curve
    # is no longer the published one and the PNG has to say so (C3 review).
    ax.plot(SCALES, reach_central, "-", color="black", lw=2.2,
            label=(r"Central realistic: mid-$x$ $\Delta$, " + r_caption
                   + r", Cloet $P_{zz}=0.267$"))

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
        + r"CT18NLO + EMC + [" + r_caption + r"/CB] + $P_{zz}$ dilution band",
        fontsize=10,
    )
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)
    share_key = fom.run_share_tag(run_share)
    outpath = outdir / ("money_delta_realistic_%s%s.png"
                        % (tag, "_" + share_key if share_key else ""))
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")

    # ── Step 5: summary numbers at Δ/F₁ = 1e-3 ──────────────────────────
    l5_central = combo_l5[central_key]
    l5_band_lo = min(combo_l5.values())
    l5_band_hi = max(combo_l5.values())
    l5_toy     = _reach(terms_toy, n_ev_toy)

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

# LOW is not part of the published money plot (it was never drawn); it is
# selectable so the third beam configuration of money_delta_20260729.py can be
# quoted on the same footing.  The default list is the frozen pair.
ALL_CONFIGS = {
    "low": (BeamConfig(electron_energy=5.0, ion=LI6,
                       ion_momentum_per_nucleon=27.5),
            "5 GeV $e$ × 27.5 GeV/u $^6$Li"),
    "mid": (BeamConfig(electron_energy=10.0, ion=LI6,
                       ion_momentum_per_nucleon=50.0),
            "10 GeV $e$ × 50 GeV/u $^6$Li"),
    "top": (BeamConfig(electron_energy=18.0, ion=LI6,
                       ion_momentum_per_nucleon=137.5),
            "18 GeV $e$ × 137.5 GeV/u $^6$Li"),
}

DEFAULT_CONFIGS = "mid,top"


def _make_stdio_portable():
    """Stop a non-UTF-8 console from killing the run on the summary table.

    Every line this script prints carries literal Unicode (L_5σ, Δ/F₁,
    fb⁻¹).  On a Windows cp936/GBK console print() raises
    UnicodeEncodeError (code review 2026-08-25 item S11, reported against
    the sibling dated script).  Changing the ERROR HANDLER and not the
    encoding leaves the byte stream bit-for-bit identical wherever the
    glyphs already encode.
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="backslashreplace")
        except (AttributeError, ValueError, OSError):
            pass   # not a reconfigurable TextIOWrapper; nothing to do


def main():
    _make_stdio_portable()
    ap = argparse.ArgumentParser(
        description="Realistic 6Li gluonometry reach: L_5σ vs Δ/F₁ (two EIC configs)."
    )
    ap.add_argument(
        "--outdir",
        default=None,
        help="Output directory for PNGs (default: <fastsim>/out/money_delta, "
             "or <fastsim>/out/money_delta_r_<model> for a non-default "
             "--r-model, so a re-run never overwrites the July figures)",
    )
    ap.add_argument(
        "--r-model",
        choices=sorted(R_MODELS),
        default=DEFAULT_R_MODEL,
        help="R = sigma_L/sigma_T behind the 'R1998' curve, which is also "
             "the central curve. 'simplified' (default) reproduces the "
             "July L_5sigma exactly; 'theta-log' is that form with Theta "
             "restricted to the log term (code review S1 repair); 'r1998' "
             "is the published SLAC/E143 fit from polli_fastsim.structure.",
    )
    ap.add_argument(
        "--run-share", type=float, default=1.0, dest="run_share",
        help="this observable's share of the programme year (plans/07 WP2; "
             "default 1.0, which every published number assumes). L_5sigma "
             "is the luminosity the measurement must accumulate and is "
             "invariant under it; a non-default share writes to its own "
             "file so the published PNGs cannot be overwritten.",
    )
    ap.add_argument(
        "--configs",
        default=DEFAULT_CONFIGS,
        help=f"Comma-separated beam configurations from "
             f"{sorted(ALL_CONFIGS)} (default: {DEFAULT_CONFIGS}, the pair "
             f"the published money plot draws)",
    )
    args = ap.parse_args()
    if not args.run_share > 0:
        ap.error("--run-share must be positive")
    r_func, r_caption = R_MODELS[args.r_model]
    r_funcs = [("R1998", r_func),
               ("Christy-Bosted", r_christy_bosted)]
    if args.outdir is not None:
        outdir = pathlib.Path(args.outdir)
    else:
        # Anchored on the file, not on the caller's cwd (code review S11).
        out_root = pathlib.Path(__file__).resolve().parents[1] / "out"
        outdir = out_root / (
            "money_delta" if args.r_model == DEFAULT_R_MODEL
            else f"money_delta_r_{args.r_model}"
        )
    try:
        config_tags = [t.strip() for t in args.configs.split(",") if t.strip()]
        configs = [(ALL_CONFIGS[t][0], t, ALL_CONFIGS[t][1]) for t in config_tags]
    except KeyError as exc:
        ap.error(f"unknown config {exc.args[0]!r}; choose from {sorted(ALL_CONFIGS)}")
    if not configs:
        # "" and "," parse to an empty list, which used to print an empty
        # summary table and exit 0 -- a silent no-op where the caller asked
        # for a reach (C3 review).
        ap.error("--configs selected no configuration; choose from "
                 f"{sorted(ALL_CONFIGS)}")

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

    summaries = []
    for cfg, tag, title_beam in configs:
        print(f"\n{'='*66}")
        print(f"Building realistic plot: {tag.upper()} config — {cfg.label()}")
        print(f"{'='*66}")
        _, summary = build_realistic_plot(cfg, base, outdir, tag, title_beam,
                                          r_funcs=r_funcs, r_caption=r_caption,
                                          run_share=args.run_share)
        summaries.append(summary)

    # ── Summary table ─────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("SUMMARY TABLE  (Δ/F₁ = 1e-3, all values in fb⁻¹/nucleon)")
    print(fom.run_share_header(L_PROGRAMME_FB, args.run_share))
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
    if args.r_model != DEFAULT_R_MODEL:
        # Loud, and only when the frozen July reproduction is NOT what ran.
        print(f"      R1998 here is --r-model {args.r_model!r}, NOT the frozen "
              f"July form: these are not the published numbers.")
    print("      Toy     = toy_delta_gluon, default R, no EMC, uniform eff, P_zz=0.8.")


if __name__ == "__main__":
    main()
