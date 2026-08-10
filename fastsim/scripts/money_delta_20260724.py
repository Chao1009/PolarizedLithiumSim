#!/usr/bin/env python3
"""6Li sum-rule-constrained A_cos2φ reach — bag-model & lattice-QCD references.

Differences from money_delta_20260721.py
──────────────────────────────────────────
1. **Δ formula changed**: instead of a dimensionless peak Δ/F₁ = scale, the
   polarized tensor structure function is defined with an explicit α_s(Q²)
   prefactor:

       Δ(x, Q²) = scale · α_s(Q²) · F₁(x, Q²) · x^α · (1-x)^β

   The `scale` parameter remains dimensionless and is still the x-axis of the
   reach plots.

2. **Sum-rule constraints**: instead of scanning scale freely, the bag-model
   and lattice-QCD sum-rule predictions are used to derive physically-motivated
   reference amplitudes A_bag and A_lat:

       ∫₀¹ x · Δ(x, Q²) dx = c · α_s(Q²)

   With c_bag = -0.012 and c_lat = -0.009. Since α_s(Q²) appears on both sides
   of the sum rule, it cancels:

       A = c / ∫₀¹ x · F₁(x, <Q²>) · x^α · (1-x)^β dx

   The residual Q² dependence enters only through F₁'s DGLAP evolution (small,
   ~few%). A rate-weighted mean Q² is used for F₁ evaluation.

3. **Four reach plots (Plots 1-4)** in place of the stat-uncertainty Plots A:
   - Plot 1: δA/A vs scale, 9 curves (3 energy pairs × 3 shapes), 2 subplots
             (10 fb⁻¹, 100 fb⁻¹). Vertical dashed lines at |A_bag|, |A_lat|.
   - Plot 2: L_5σ vs scale, 9 curves, single panel. Same verticals.
   - Plot 3: bag vs lattice, mid config only, δA/A vs scale, 2 subplots.
   - Plot 4: bag vs lattice, mid config only, L_5σ vs scale, single panel.

4. **Four relative per-bin heatmaps (Plots 5-8)**: |δA_bin / A_cos2φ,bin| on
   log x × log Q² axes using the bag sum-rule Δ at mid_x shape, magma colormap.
   For (MID, TOP) × (10, 100 fb⁻¹/u).

5. **R parameterization fixed at R1998 throughout** (no Christy-Bosted).

6. **α_s(Q²) source**: primary = interpolation from the `parton` PDF set's
   tabulated AlphaS_Qs / AlphaS_Vals in pdfset.info (from CT18NLO). Fallback to
   LO analytic formula if the primary path fails.

Physics note on the α_s cancellation
──────────────────────────────────────
The sum rule ∫ x Δ dx = c · α_s(Q²) with Δ = scale · α_s(Q²) · F₁ · shape
gives: scale · α_s(Q²) · ∫ x F₁ shape dx = c · α_s(Q²). The α_s(Q²) cancels
on both sides. The resulting A = c / ∫ x F₁(x, <Q²>) shape dx depends on Q²
only through F₁'s DGLAP evolution — a small (~few%) effect evaluated at the
rate-weighted mean Q² of the accepted phase space.

Sign convention: c_bag = -0.012, c_lat = -0.009 (both negative). This gives
negative A values. Since sig² ∝ A², the sign does not affect the observable,
but is preserved for physics honesty. Plots use |A| for the vertical-line x
position.

Three EIC beam configurations for ⁶Li
───────────────────────────────────────
LOW : E_e = 5 GeV,  p_ion = 27.5 GeV/u
MID : E_e = 10 GeV, p_ion = 50 GeV/u
TOP : E_e = 18 GeV, p_ion = 137.5 GeV/u

Output
──────
fastsim/out/money_delta/money_delta_20260724_p1_bag_10fb.png
fastsim/out/money_delta/money_delta_20260724_p1_bag_100fb.png
fastsim/out/money_delta/money_delta_20260724_p2_bag.png
fastsim/out/money_delta/money_delta_20260724_p3_baglat.png
fastsim/out/money_delta/money_delta_20260724_p4_baglat.png
fastsim/out/money_delta/money_delta_20260724_perbin_mid_10fb.png
fastsim/out/money_delta/money_delta_20260724_perbin_mid_100fb.png
fastsim/out/money_delta/money_delta_20260724_perbin_top_10fb.png
fastsim/out/money_delta/money_delta_20260724_perbin_top_100fb.png
"""

import argparse
import pathlib
import sys
from contextlib import contextmanager

import numpy as np

# ── ensure polli_fastsim is importable ────────────────────────────────────────
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import fom
from polli_fastsim.beams import LI6, BeamConfig
from polli_fastsim.asymmetries import a_cos2phi
from polli_fastsim.inputs import get_backends
import polli_fastsim.structure as structure_mod
import polli_fastsim.asymmetries as asymmetries_mod

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# ══════════════════════════════════════════════════════════════════════════════
# NuclearF2FromGrid — local class; do NOT modify polli_fastsim/structure.py
# (copied verbatim from money_delta_20260721.py)
# ══════════════════════════════════════════════════════════════════════════════

def _safe_xfx_local(pdf, pid, x, q2):
    """Coerce xfxQ2 return value to a plain Python float.

    Guards against NumPy ≥ 2.0 where float() on a 0-d ndarray raises
    TypeError instead of extracting the scalar.  Mirrors _safe_xfx in
    polli_fastsim/structure.py.
    """
    v = pdf.xfxQ2(pid, x, q2)
    try:
        return float(v)
    except TypeError:
        return float(np.asarray(v).item())


class NuclearF2FromGrid:
    """F₂ᴬ from a nuclear PDF set loaded via the `parton` package.

    Drop-in replacement for NuclearF2 (polli_fastsim.structure) when the
    nucleus-integrated structure function comes directly from a nuclear PDF
    rather than from Z*F2p + N*F2n with a free-proton PDF.

    API contract (same methods called by project_rates and downstream code):
    - f2a(x, q2)  → whole-nucleus F₂ᴬ  (callers divide by ion.A for /nucleon)
    - f1a(x, q2)  → whole-nucleus F₁ᴬ  via Callan-Gross + active r_sigma_lt

    Normalization
    ─────────────
    The `parton` package returns nuclear PDFs *per nucleon* (confirmed by
    numerical cross-check against CT18NLO at DIS kinematics for Li6).
    f2a() therefore returns  A × Σ_q e_q² (x·q_A + x·q̄_A)  where the sum
    over A-averaged parton densities already includes all A nucleons but is
    expressed per nucleon.  The ×A restores the whole-nucleus normalisation
    expected by the callers.

    Parameters
    ----------
    ion : beams.Ion
        The target ion (used for ion.A and for the f1a / f2a /A convention).
    pdf_set_name : str
        Name of the nuclear PDF set as used by parton.mkPDF, e.g.
        'EPPS21nlo_CT18Anlo_Li6'.
    member : int
        PDF replica/member index (default 0 = central).
    """

    # e_q²: quark electric charges squared (same as PartonF2 in structure.py)
    _E2 = {1: 1 / 9, 2: 4 / 9, 3: 1 / 9, 4: 4 / 9, 5: 1 / 9}

    def __init__(self, ion, pdf_set_name, member=0):
        from parton import mkPDF  # lazy: optional dependency
        self.ion = ion
        self._pdf = mkPDF(pdf_set_name, member)
        self._pdf_set_name = pdf_set_name
        # vectorize the scalar worker for numpy-array inputs
        self._f2a_vec = np.vectorize(self._f2a_scalar)

    def _f2a_scalar(self, x, q2):
        """Per-call scalar worker: returns whole-nucleus F₂ᴬ at (x, Q²)."""
        if not (0.0 < x < 1.0):
            return 0.0
        tot = 0.0
        for pid, e2 in self._E2.items():
            tot += e2 * (
                _safe_xfx_local(self._pdf, pid,  x, q2)
                + _safe_xfx_local(self._pdf, -pid, x, q2)
            )
        # `parton` returns per-nucleon xf_A; multiply by A for whole-nucleus
        return max(tot, 0.0) * self.ion.A

    def f2a(self, x, q2):
        """Whole-nucleus F₂ᴬ(x, Q²) — accepts numpy arrays."""
        return self._f2a_vec(
            np.asarray(x, dtype=float),
            np.asarray(q2, dtype=float),
        )

    def f1a(self, x, q2):
        """Whole-nucleus F₁ᴬ via Callan-Gross + currently-active r_sigma_lt.

        Uses the module-level structure_mod.r_sigma_lt, which is monkey-patched
        by the r_override() context manager so that R1998/Christy-Bosted
        swapping works identically to the parent script.
        """
        x  = np.asarray(x,  dtype=float)
        q2 = np.asarray(q2, dtype=float)
        r  = structure_mod.r_sigma_lt(x, q2)
        return self.f2a(x, q2) / (2.0 * x * (1.0 + r))


# ══════════════════════════════════════════════════════════════════════════════
# Δ shape variants (copied verbatim from money_delta_20260720.py)
# ══════════════════════════════════════════════════════════════════════════════

# Peak of x^α (1-x)^β is at x_peak = α/(α+β), peak value = α^α β^β / (α+β)^(α+β).

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


def peak_shape_value(alpha, beta):
    """Peak value of x^alpha * (1-x)^beta at x_peak = alpha/(alpha+beta)."""
    x_peak = alpha / (alpha + beta)
    return (x_peak ** alpha) * ((1.0 - x_peak) ** beta)


def scale_to_peak_delta_f1(scale_val, q2_ref, base, alpha, beta):
    """Convert code-parameter `scale` to peak Δ/F₁.

    peak Δ/F₁ = scale · α_s(Q²_ref) · max_x[x^α (1-x)^β]

    Uses the rate-weighted mean Q² and the specified shape variant.
    """
    return scale_val * alpha_s(q2_ref, base=base) * peak_shape_value(alpha, beta)


# Sanity check: at x = x_peak, f1 = 1, Δ (without α_s) should equal `scale`.
def _sanity_check_shapes():
    for name, (alpha, beta) in _VARIANTS.items():
        xp = alpha / (alpha + beta)
        test_scale = 1.0
        f1_test = 1.0
        norm = test_scale / _PEAK_VALS[name]
        delta_at_peak = norm * f1_test * (xp ** alpha) * ((1.0 - xp) ** beta)
        assert abs(delta_at_peak - test_scale) < 1e-10, (
            f"Shape sanity FAILED for variant '{name}': "
            f"Δ_unnorm(x_peak, f1=1) = {delta_at_peak}, expected {test_scale}"
        )


_sanity_check_shapes()


def delta_shape(x, q2, f1, scale, variant="mid_x"):
    """Δ(x,Q²) without α_s: normalized so peak Δ/F₁ = scale.

    This is the base shape used by delta_shape_with_alphas.
    Parameters
    ----------
    x, q2 : array-like  kinematic variables
    f1    : array-like  unpolarized F₁ per nucleon (same shape)
    scale : float       dimensionless amplitude parameter (x-axis of reach plots)
    variant : str       one of 'low_x', 'mid_x', 'high_x'
    """
    x = np.asarray(x, dtype=float)
    alpha, beta = _VARIANTS[variant]
    norm = scale / _PEAK_VALS[variant]
    return norm * f1 * np.power(np.maximum(x, 1e-12), alpha) * np.power(np.maximum(1.0 - x, 0.0), beta)


# ══════════════════════════════════════════════════════════════════════════════
# R = σ_L/σ_T parameterization (R1998 only in this script)
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
# Fixed physics parameters
# ══════════════════════════════════════════════════════════════════════════════

# P_zz: Cloet per-nucleon-normalized convention = 0.8 / 3
PZZ = 0.267

# Sum-rule coefficients (both negative — see docstring)
C_BAG = -0.012
C_LAT = -0.009

# Luminosities to evaluate [fb⁻¹/nucleon]
LUMINOSITIES_FB = [10.0, 100.0]

# Grid parameters (identical to parent scripts)
# Anchor SCALES range to bag sum-rule reference for MID+mid_x: |A_bag| = 0.310
# Range: 0.1 * |A_bag_ref| to 10 * |A_bag_ref| = [0.031, 3.10]
_S0_BAG_REF = 0.310  # |A_bag| for MID+mid_x, computed a priori for range anchoring
SCALES = np.logspace(np.log10(0.1 * _S0_BAG_REF), np.log10(10.0 * _S0_BAG_REF), 15)   # [0.031, 3.10]
S0 = 1e-3           # reference scale for analytic rescaling

# Rate-weighted mean Q² for MID config (used as reference for peak Δ/F₁ conversion)
# This is the value from solve_A_from_sum_rule for mid config / mid_x shape.
# Used only for x-axis labelling in Plots 1-4; does NOT affect any physics computation.
_Q2_MEAN_MID_REF = 7.38   # GeV²

# Reference shape parameters for x-axis conversion (mid_x: α=0.7, β=3)
_REF_ALPHA = 0.7
_REF_BETA  = 3.0

# Min events cut (same as parent)
MIN_EVENTS = 10

# Nuclear PDF set name (fixed throughout)
EPPS21_SET = "EPPS21nlo_CT18Anlo_Li6"

# Shape variants for reach plots
VARIANT_NAMES = ["low_x", "mid_x", "high_x"]

# Reference variant for sum-rule vertical lines in Plots 1 & 2
REF_VARIANT = "mid_x"

# Plot styling: color = energy config, linestyle = shape
CONFIG_COLORS = {
    "low":  "crimson",
    "mid":  "seagreen",
    "top":  "navy",
}
VARIANT_LINESTYLES = {
    "low_x":  "-",
    "mid_x":  "--",
    "high_x": ":",
}
BAG_COLOR = "black"
LAT_COLOR = "darkred"


# ══════════════════════════════════════════════════════════════════════════════
# α_s(Q²) helper — primary: parton PDF table; fallback: LO analytic
# ══════════════════════════════════════════════════════════════════════════════

_ALPHAS_SOURCE = [None]   # mutable sentinel: filled once on first call


def _build_alphas_interpolator(base):
    """Build a NumPy interpolation table from the parton PDF info.

    Parameters
    ----------
    base : PartonF2 or None
        If None, fall back to LO analytic.

    Returns
    -------
    str  : 'parton' or 'LO'
    callable : alpha_s(q2_array) → array
    """
    if base is not None:
        try:
            info = base._pdf.pdfset.info
            qs_arr   = np.array(info["AlphaS_Qs"],   dtype=float)
            vals_arr = np.array(info["AlphaS_Vals"],  dtype=float)
            q2s_arr  = qs_arr ** 2  # convert Q → Q²

            def _from_parton(q2):
                q2 = np.asarray(q2, dtype=float)
                q  = np.sqrt(np.maximum(q2, q2s_arr[0]))
                return np.interp(q, qs_arr, vals_arr)

            return "parton", _from_parton
        except Exception as exc:
            print(
                f"[alpha_s] WARNING: parton table extraction failed ({exc}); "
                "falling back to LO analytic formula."
            )

    # LO analytic fallback
    LAMBDA_QCD = 0.22   # GeV
    N_F        = 4      # fixed n_f (simplified; valid for Q² << m_b²)
    BETA0      = 33 - 2 * N_F   # = 25

    def _lo_analytic(q2):
        q2 = np.asarray(q2, dtype=float)
        lnq2lam2 = np.log(np.maximum(q2, LAMBDA_QCD**2 * 1.01) / LAMBDA_QCD**2)
        return 12.0 * np.pi / (BETA0 * lnq2lam2)

    return "LO", _lo_analytic


def alpha_s(q2, base=None):
    """Return α_s(Q²) for scalar or array Q² inputs.

    Primary source: tabulated AlphaS_Qs / AlphaS_Vals from parton PDF info.
    Fallback: LO analytic: α_s = 12π / [(33 - 2 n_f) ln(Q²/Λ²)],
              Λ = 0.22 GeV, n_f = 4.

    At first call, prints a one-line note indicating which source is used.
    Subsequent calls are silent and use the cached interpolator.

    Parameters
    ----------
    q2   : float or array-like   Q² in GeV²
    base : PartonF2 or None      CT18NLO backend (for primary path)

    Returns
    -------
    float or ndarray
    """
    if _ALPHAS_SOURCE[0] is None:
        source, func = _build_alphas_interpolator(base)
        _ALPHAS_SOURCE[0] = (source, func)
        print(f"[alpha_s] Using source: {source}")
    _, func = _ALPHAS_SOURCE[0]
    result = func(np.asarray(q2, dtype=float))
    # Return scalar if scalar input was passed
    if result.ndim == 0:
        return float(result)
    return result


# ══════════════════════════════════════════════════════════════════════════════
# Δ with α_s(Q²) prefactor
# ══════════════════════════════════════════════════════════════════════════════

def delta_shape_with_alphas(x, q2, f1, scale, variant, base):
    """Δ(x, Q²) = scale · α_s(Q²) · F₁(x, Q²) · x^α · (1-x)^β.

    The normalization of x^α (1-x)^β uses the same _PEAK_VALS as delta_shape,
    so that at x = x_peak, f1 = 1, α_s = 1: Δ = scale.
    At physical α_s ~ 0.2, the peak Δ/F₁ ≈ scale × 0.2.

    Parameters
    ----------
    x, q2    : array-like  kinematic variables
    f1       : array-like  unpolarized F₁ per nucleon (same shape as x, q2)
    scale    : float       dimensionless amplitude (x-axis of reach plots)
    variant  : str         one of 'low_x', 'mid_x', 'high_x'
    base     : PartonF2    CT18NLO backend (for α_s primary path)

    Returns
    -------
    ndarray  : Δ(x, Q²) per nucleon
    """
    x   = np.asarray(x,   dtype=float)
    q2  = np.asarray(q2,  dtype=float)
    f1  = np.asarray(f1,  dtype=float)
    alpha_beta = _VARIANTS[variant]
    alpha_v, beta_v = alpha_beta
    norm = scale / _PEAK_VALS[variant]
    as_val = alpha_s(q2, base=base)   # array of α_s values at each bin's Q²
    return (norm * as_val * f1
            * np.power(np.maximum(x, 1e-12), alpha_v)
            * np.power(np.maximum(1.0 - x, 0.0), beta_v))


# ══════════════════════════════════════════════════════════════════════════════
# Sum-rule solver: A = c / ∫₀¹ x F₁(x, <Q²>) shape(x) dx
# ══════════════════════════════════════════════════════════════════════════════

def solve_A_from_sum_rule(cfg, nf2_obj, base, variant, c_coef):
    """Solve for A = scale such that the sum rule ∫ x Δ dx = c · α_s(Q²) holds.

    Because Δ = scale · α_s(Q²) · F₁ · x^α · (1-x)^β, the α_s(Q²) cancels:
        scale · ∫₀¹ x · F₁(x, <Q²>) · x^α · (1-x)^β dx = c
    where <Q²> is the rate-weighted mean Q² of the accepted phase space.

    Parameters
    ----------
    cfg      : BeamConfig
    nf2_obj  : NuclearF2FromGrid  pre-constructed EPPS21 object
    base     : PartonF2           CT18NLO backend (for α_s; not used in integral)
    variant  : str                shape variant ('low_x', 'mid_x', 'high_x')
    c_coef   : float              sum-rule coefficient (c_bag or c_lat; negative)

    Returns
    -------
    float : A = c_coef / integral  (negative for both bag and lattice)
    float : <Q²>  (rate-weighted mean Q² of accepted bins, in GeV²)
    """
    # Run project_rates at unit luminosity to get the acceptance mask and rates
    sc   = fom.Scenario(lumi_fb_per_nucleon=1.0, pol_ion_tensor=PZZ, q2_min=2.0)
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_obj)

    # Rate-weighted mean Q² over accepted bins
    n_events = proj.n_events   # shape (nx, nq2)
    accepted = proj.accepted   # bool mask

    n_acc = n_events[accepted]
    q2_acc = proj.q2[accepted]

    total_events = n_acc.sum()
    if total_events <= 0:
        raise RuntimeError(
            f"No accepted events in solve_A_from_sum_rule for config {cfg.label()}"
        )
    q2_mean = float((n_acc * q2_acc).sum() / total_events)

    # Numerical integral: ∫₀¹ x · F₁(x, <Q²>) · x^α · (1-x)^β dx
    # Use the accepted x-grid plus a fine extension to x → 0
    # Grid centers from project_rates: proj.x is shape (nx, nq2), proj.x[:, 0] is x-centers
    x_grid_accepted = np.unique(proj.x[accepted])   # sorted unique x values

    # Extend to smaller x for the low-x tail (important for low_x variant)
    x_lo_ext = np.logspace(np.log10(1e-5), np.log10(x_grid_accepted.min() * 0.99), 50)
    x_hi_ext = np.linspace(x_grid_accepted.max() * 1.001, 0.9999, 30)
    x_full = np.concatenate([x_lo_ext, x_grid_accepted, x_hi_ext])
    x_full = np.sort(np.unique(x_full))
    x_full = x_full[(x_full > 0) & (x_full < 1.0)]

    # Evaluate F₁ per nucleon at <Q²>
    q2_arr = np.full_like(x_full, q2_mean)
    f1_arr = nf2_obj.f1a(x_full, q2_arr) / cfg.ion.A   # per nucleon

    alpha_v, beta_v = _VARIANTS[variant]
    integrand = (x_full
                 * f1_arr
                 * np.power(np.maximum(x_full, 1e-12), alpha_v)
                 * np.power(np.maximum(1.0 - x_full, 0.0), beta_v)
                 / _PEAK_VALS[variant])   # includes the 1/_PEAK normalization
    # _PEAK_VALS normalization: the shape function here is
    # x^α (1-x)^β / _PEAK_VALS[variant]  so that at x_peak the shape is 1
    # matching the convention in delta_shape / delta_shape_with_alphas

    integral = float(np.trapezoid(integrand, x_full))

    if abs(integral) < 1e-30:
        raise RuntimeError(
            f"Sum-rule integral nearly zero for variant={variant}, "
            f"config={cfg.label()}"
        )

    A = c_coef / integral
    return A, q2_mean


# ══════════════════════════════════════════════════════════════════════════════
# sig² per fb⁻¹ with α_s(Q²) prefactor in Δ
# ══════════════════════════════════════════════════════════════════════════════

def sig2_per_fb_at_sumrule(cfg, scale, pzz, nf2_obj, base, variant="mid_x",
                            min_events=10):
    """Sig² per fb⁻¹/nucleon using Δ = scale · α_s(Q²) · F₁ · x^α · (1-x)^β.

    Parameters
    ----------
    cfg      : BeamConfig
    scale    : float  dimensionless amplitude (reference scale for analytic rescaling)
    pzz      : float  effective P_zz per nucleon
    nf2_obj  : NuclearF2FromGrid  pre-constructed EPPS21 object
    base     : PartonF2           CT18NLO backend (for α_s)
    variant  : str    Δ shape variant ('low_x', 'mid_x', 'high_x')
    min_events : int  minimum events to include a bin

    The R-override must be active in the calling context.

    Returns
    -------
    float : sig² per fb⁻¹/nucleon at the given scale
    """
    sc   = fom.Scenario(lumi_fb_per_nucleon=1.0, pol_ion_tensor=pzz, q2_min=2.0)
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_obj)

    nf2 = proj.extras["nf2"]
    f1  = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2  = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y   = proj.extras["y"]

    delta = delta_shape_with_alphas(proj.x, proj.q2, f1, scale=scale,
                                    variant=variant, base=base)
    amp   = a_cos2phi(delta, f1, f2, proj.x, y)

    n_events = proj.n_events
    use = proj.accepted & (n_events >= min_events)
    return float(np.where(use, amp**2 * pzz**2 * n_events / 2.0, 0.0).sum())


# ══════════════════════════════════════════════════════════════════════════════
# Plot 1: δA/A vs scale, 9 curves, 2 subplots (10 fb⁻¹ and 100 fb⁻¹)
# ══════════════════════════════════════════════════════════════════════════════

def _build_plot1_one_lumi(sig2_dict, a_vals, luminosity, outdir):
    """Helper: produce one single-panel PNG for Plot 1 at a given luminosity.

    Parameters
    ----------
    sig2_dict  : dict  (config_tag, variant) → float  (sig² per fb⁻¹ at S0)
    a_vals     : dict  (config_tag, variant, model) → float
    luminosity : float  [fb⁻¹/nucleon]
    outdir     : pathlib.Path

    Returns
    -------
    pathlib.Path : output PNG path
    """
    from matplotlib.lines import Line2D

    L = luminosity
    a_bag_ref = abs(a_vals.get(("mid", "mid_x", "bag"), np.nan))
    a_lat_ref = abs(a_vals.get(("mid", "mid_x", "lat"), np.nan))

    configs_ordered = ["low", "mid", "top"]
    config_labels = {
        "low": r"5$\times$27.5 GeV/u",
        "mid": r"10$\times$50 GeV/u",
        "top": r"18$\times$137.5 GeV/u",
    }

    fig, ax = plt.subplots(figsize=(9, 6))

    for cfg_tag in configs_ordered:
        color = CONFIG_COLORS[cfg_tag]
        for variant in VARIANT_NAMES:
            ls = VARIANT_LINESTYLES[variant]
            key = (cfg_tag, variant)
            s2 = sig2_dict.get(key)
            if s2 is None or s2 <= 0:
                continue
            sig2_scaled = s2 * (SCALES / S0) ** 2 * L
            rel_unc = 1.0 / np.sqrt(np.maximum(sig2_scaled, 1e-60))
            label = (f"{cfg_tag.upper()} {config_labels[cfg_tag]}, {variant}")
            ax.plot(SCALES, rel_unc, ls=ls, color=color, lw=1.4,
                    label=label, alpha=0.85)

    # 5σ discovery threshold
    ax.axhline(0.20, color="firebrick", ls="--", lw=1.5, zorder=4,
               label="5σ discovery (δA/A=0.20)")

    # Vertical lines at raw |A| (scale) positions for bag and lattice sum-rule references
    if np.isfinite(a_bag_ref) and a_bag_ref > 0:
        ax.axvline(a_bag_ref, color=BAG_COLOR, ls="--", lw=1.3, zorder=5,
                   label=f"bag sum rule: $s_0^{{\\rm bag}}$ = {a_bag_ref:.3f}")
        ax.text(a_bag_ref * 1.07, 0.92,
                "bag", fontsize=5.5, color=BAG_COLOR,
                transform=ax.get_xaxis_transform(), va="top")
    if np.isfinite(a_lat_ref) and a_lat_ref > 0:
        ax.axvline(a_lat_ref, color=LAT_COLOR, ls="--", lw=1.3, zorder=5,
                   label=f"lattice sum rule: $s_0^{{\\rm lat}}$ = {a_lat_ref:.3f}")
        ax.text(a_lat_ref * 1.07, 0.78,
                "lat", fontsize=5.5, color=LAT_COLOR,
                transform=ax.get_xaxis_transform(), va="top")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"scale parameter $s$", fontsize=10)
    ax.set_ylabel(
        r"$\delta A / |A_{\cos 2\phi}|$ (combined)",
        fontsize=10
    )
    ax.set_xlim(SCALES[0] * 0.85, SCALES[-1] * 1.2)

    ax.text(0.02, 0.05,
            r"$\Delta(x,Q^2) = s \cdot \alpha_s(Q^2) \cdot F_1 \cdot x^\alpha (1-x)^\beta$",
            transform=ax.transAxes, fontsize=6, color="dimgray", va="bottom")

    ax.legend(fontsize=7, ncol=2, loc="upper right", framealpha=0.7)

    ax.set_title(
        f"$^6$Li δA/|A_cos2φ| vs scale — L = {luminosity:g} fb⁻¹/u\n"
        f"EPPS21, R1998, Cloet $P_{{zz}}$=0.267; sum-rule refs $s_0^{{\\rm bag}}$, $s_0^{{\\rm lat}}$",
        fontsize=10,
    )
    fig.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)
    lumi_tag = f"{int(L)}fb" if L == int(L) else f"{L:g}fb"
    outpath = outdir / f"money_delta_20260724_p1_bag_{lumi_tag}.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")
    return outpath


def build_plot1(sig2_dict, a_vals, outdir):
    """Plot 1: δA/A vs scale, 3 configs × 3 shapes = 9 curves.

    Produces two separate single-panel PNGs — one per luminosity (10 and
    100 fb⁻¹/u).

    Parameters
    ----------
    sig2_dict : dict
        Keys: (config_tag, variant) → float  (sig² per fb⁻¹ at S0 = 1e-3)
    a_vals : dict
        Keys: (config_tag, variant, model) → float  (A from sum rule)
        model: 'bag' or 'lat'
    outdir : pathlib.Path

    Returns
    -------
    list[pathlib.Path] : output PNG paths [10 fb⁻¹ path, 100 fb⁻¹ path]
    """
    paths = []
    for L in LUMINOSITIES_FB:
        p = _build_plot1_one_lumi(sig2_dict, a_vals, L, outdir)
        paths.append(p)
    return paths


# ══════════════════════════════════════════════════════════════════════════════
# Plot 2: L_5σ vs scale, 9 curves, single panel
# ══════════════════════════════════════════════════════════════════════════════

def build_plot2(sig2_dict, a_vals, outdir):
    """Plot 2: L_5σ vs scale, 3 configs × 3 shapes = 9 curves, single panel.

    Parameters
    ----------
    sig2_dict : dict  (same as build_plot1)
    a_vals    : dict  (same as build_plot1)
    outdir    : pathlib.Path

    Returns
    -------
    pathlib.Path : output PNG path
    """
    a_bag_ref = abs(a_vals.get(("mid", "mid_x", "bag"), np.nan))
    a_lat_ref = abs(a_vals.get(("mid", "mid_x", "lat"), np.nan))

    configs_ordered = ["low", "mid", "top"]
    config_labels = {
        "low": r"5$\times$27.5 GeV/u",
        "mid": r"10$\times$50 GeV/u",
        "top": r"18$\times$137.5 GeV/u",
    }

    fig, ax = plt.subplots(figsize=(9, 6))

    for cfg_tag in configs_ordered:
        color = CONFIG_COLORS[cfg_tag]
        for variant in VARIANT_NAMES:
            ls = VARIANT_LINESTYLES[variant]
            key = (cfg_tag, variant)
            s2 = sig2_dict.get(key)
            if s2 is None or s2 <= 0:
                continue
            l5sigma = 25.0 / np.maximum(s2 * (SCALES / S0) ** 2, 1e-30)
            label = (f"{cfg_tag.upper()} {config_labels[cfg_tag]}, {variant}")
            ax.plot(SCALES, l5sigma, ls=ls, color=color, lw=1.4,
                    label=label, alpha=0.85)

    # Gold luminosity band: 1–100 fb⁻¹/u "plausible EIC program"
    ax.axhspan(1, 100, color="gold", alpha=0.15,
               label=r"1–100 fb$^{-1}$/u (plausible EIC program)")

    # Vertical lines at raw |A| (scale) positions for bag and lattice sum-rule references
    if np.isfinite(a_bag_ref) and a_bag_ref > 0:
        ax.axvline(a_bag_ref, color=BAG_COLOR, ls="--", lw=1.3, zorder=5,
                   label=f"bag sum rule: $s_0^{{\\rm bag}}$ = {a_bag_ref:.3f}")
        ax.text(a_bag_ref * 1.07, 0.92,
                "bag", fontsize=5.5, color=BAG_COLOR,
                transform=ax.get_xaxis_transform(), va="top")
    if np.isfinite(a_lat_ref) and a_lat_ref > 0:
        ax.axvline(a_lat_ref, color=LAT_COLOR, ls="--", lw=1.3, zorder=5,
                   label=f"lattice sum rule: $s_0^{{\\rm lat}}$ = {a_lat_ref:.3f}")
        ax.text(a_lat_ref * 1.07, 0.75,
                "lat", fontsize=5.5, color=LAT_COLOR,
                transform=ax.get_xaxis_transform(), va="top")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"scale parameter $s$", fontsize=11)
    ax.set_ylabel(r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]", fontsize=11)
    ax.set_xlim(SCALES[0] * 0.85, SCALES[-1] * 1.2)
    ax.text(0.02, 0.05,
            r"$\Delta(x,Q^2) = s \cdot \alpha_s(Q^2) \cdot F_1 \cdot x^\alpha (1-x)^\beta$",
            transform=ax.transAxes, fontsize=6, color="dimgray", va="bottom")
    ax.set_title(
        r"$^6$Li $L_{5\sigma}$ vs scale $s$ — 3 energies × 3 shapes"
        "\n"
        r"EPPS21, R1998, Cloet $P_{zz}$=0.267; sum-rule refs $s_0^{\rm bag}$, $s_0^{\rm lat}$",
        fontsize=10,
    )

    # Legend: two-column
    from matplotlib.lines import Line2D
    shape_handles = [
        Line2D([0], [0], color="gray", ls="-",  lw=1.4, label="low_x shape"),
        Line2D([0], [0], color="gray", ls="--", lw=1.4, label="mid_x shape"),
        Line2D([0], [0], color="gray", ls=":",  lw=1.4, label="high_x shape"),
    ]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles + shape_handles, labels + [h.get_label() for h in shape_handles],
              fontsize=5.5, ncol=2, loc="upper right", framealpha=0.7)

    fig.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "money_delta_20260724_p2_bag.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")
    return outpath


# ══════════════════════════════════════════════════════════════════════════════
# Plot 3: bag vs lattice, mid config, δA/A vs scale, 2 subplots
# ══════════════════════════════════════════════════════════════════════════════

def build_plot3(sig2_mid_midx, a_vals, outdir):
    """Plot 3: δA/A vs scale — mid config, mid_x shape, both luminosities.

    Produces one single-panel PNG with two curves (L=10 solid black,
    L=100 dashed black) and two vertical lines at |A_bag| and |A_lat|.

    Parameters
    ----------
    sig2_mid_midx : float  sig² per fb⁻¹ at S0 for mid config, mid_x shape
    a_vals        : dict   (config_tag, variant, model) → float
    outdir        : pathlib.Path

    Returns
    -------
    pathlib.Path : output PNG path
    """
    a_bag = abs(a_vals.get(("mid", "mid_x", "bag"), np.nan))
    a_lat = abs(a_vals.get(("mid", "mid_x", "lat"), np.nan))

    s2 = sig2_mid_midx

    fig, ax = plt.subplots(figsize=(9, 6))

    # Two curves: one per luminosity
    lumi_styles = [(10.0, "-", "L = 10 fb$^{-1}$/u"),
                   (100.0, "--", "L = 100 fb$^{-1}$/u")]
    for L, ls, lumi_label in lumi_styles:
        sig2_scaled = s2 * (SCALES / S0) ** 2 * L
        rel_unc = 1.0 / np.sqrt(np.maximum(sig2_scaled, 1e-60))
        ax.plot(SCALES, rel_unc, ls, color="black", lw=2.0, label=lumi_label)

    # 5σ threshold
    ax.axhline(0.20, color="firebrick", ls="--", lw=1.4, zorder=4,
               label="5σ discovery")

    # Vertical lines at raw |A| (scale) positions for bag and lattice sum-rule references
    if np.isfinite(a_bag) and a_bag > 0:
        ax.axvline(a_bag, color=BAG_COLOR, ls="--", lw=1.3, zorder=5,
                   label=f"bag sum rule: $s_0^{{\\rm bag}}$ = {a_bag:.3f}")
        ax.text(a_bag * 1.07, 0.92, "bag", fontsize=6, color=BAG_COLOR,
                transform=ax.get_xaxis_transform(), va="top")
    if np.isfinite(a_lat) and a_lat > 0:
        ax.axvline(a_lat, color=LAT_COLOR, ls="--", lw=1.3, zorder=5,
                   label=f"lattice sum rule: $s_0^{{\\rm lat}}$ = {a_lat:.3f}")
        ax.text(a_lat * 1.07, 0.78, "lat", fontsize=6, color=LAT_COLOR,
                transform=ax.get_xaxis_transform(), va="top")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"scale parameter $s$", fontsize=10)
    ax.set_ylabel(
        r"$\delta A / |A_{\cos 2\phi}|$ (relative stat. uncertainty)",
        fontsize=10
    )
    ax.set_xlim(SCALES[0] * 0.85, SCALES[-1] * 1.2)
    ax.text(0.02, 0.05,
            r"$\Delta(x,Q^2) = s \cdot \alpha_s(Q^2) \cdot F_1 \cdot x^\alpha (1-x)^\beta$",
            transform=ax.transAxes, fontsize=6, color="dimgray", va="bottom")
    ax.set_title(
        "$^6$Li δA/|A_cos2φ| vs scale — mid config (10 × 50 GeV/u), mid_x shape\n"
        "Bag vs lattice sum-rule references, EPPS21, R1998, Cloet $P_{zz}$=0.267",
        fontsize=10,
    )
    ax.legend(fontsize=8, loc="upper right", framealpha=0.7)
    fig.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "money_delta_20260724_p3_baglat.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")
    return outpath


# ══════════════════════════════════════════════════════════════════════════════
# Plot 4: bag vs lattice, mid config, L_5σ vs scale, single panel
# ══════════════════════════════════════════════════════════════════════════════

def build_plot4(sig2_mid_midx, a_vals, outdir):
    """Plot 4: L_5σ vs scale for mid config, mid_x shape, with bag/lat verticals.

    Since L_5σ vs scale is identical for bag and lattice (only the vertical
    line position differs), this plot has ONE curve with two vertical lines.

    Parameters
    ----------
    sig2_mid_midx : float  sig² per fb⁻¹ at S0 for mid config, mid_x shape
    a_vals        : dict   (config_tag, variant, model) → float
    outdir        : pathlib.Path

    Returns
    -------
    pathlib.Path : output PNG path
    """
    a_bag = abs(a_vals.get(("mid", "mid_x", "bag"), np.nan))
    a_lat = abs(a_vals.get(("mid", "mid_x", "lat"), np.nan))

    s2     = sig2_mid_midx
    l5sigma = 25.0 / np.maximum(s2 * (SCALES / S0) ** 2, 1e-30)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(SCALES, l5sigma, "-", color=BAG_COLOR, lw=2.2,
            label=r"MID config, mid_x shape, R1998, EPPS21")

    # Gold luminosity band
    ax.axhspan(1, 100, color="gold", alpha=0.15,
               label=r"1–100 fb$^{-1}$/u (plausible EIC program)")

    # Vertical lines at raw |A| (scale) positions for bag and lattice sum-rule references
    bag_label = (f"bag sum rule: $s_0^{{\\rm bag}}$ = {a_bag:.3f}"
                 if np.isfinite(a_bag) else "Bag: N/A")
    lat_label = (f"lattice sum rule: $s_0^{{\\rm lat}}$ = {a_lat:.3f}"
                 if np.isfinite(a_lat) else "Lattice: N/A")

    if np.isfinite(a_bag) and a_bag > 0:
        ax.axvline(a_bag, color=BAG_COLOR, ls="--", lw=1.5, zorder=5,
                   label=bag_label)
        ax.text(a_bag * 1.07, 0.92, "bag", fontsize=7, color=BAG_COLOR,
                transform=ax.get_xaxis_transform(), va="top")
    if np.isfinite(a_lat) and a_lat > 0:
        ax.axvline(a_lat, color=LAT_COLOR, ls="--", lw=1.5, zorder=5,
                   label=lat_label)
        ax.text(a_lat * 1.07, 0.78, "lat", fontsize=7, color=LAT_COLOR,
                transform=ax.get_xaxis_transform(), va="top")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"scale parameter $s$", fontsize=11)
    ax.set_ylabel(r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]", fontsize=11)
    ax.set_xlim(SCALES[0] * 0.85, SCALES[-1] * 1.2)
    ax.text(0.02, 0.05,
            r"$\Delta(x,Q^2) = s \cdot \alpha_s(Q^2) \cdot F_1 \cdot x^\alpha (1-x)^\beta$",
            transform=ax.transAxes, fontsize=6, color="dimgray", va="bottom")
    ax.set_title(
        r"$^6$Li $L_{5\sigma}$ vs scale $s$ — mid config, mid_x shape"
        "\n"
        f"Bag ($s_0^{{\\rm bag}}$={a_bag:.3f}) and lattice ($s_0^{{\\rm lat}}$={a_lat:.3f}) "
        r"sum-rule references, EPPS21, R1998",
        fontsize=10,
    )
    ax.legend(fontsize=8, loc="upper right", framealpha=0.7)
    fig.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "money_delta_20260724_p4_baglat.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")
    return outpath


# ══════════════════════════════════════════════════════════════════════════════
# Plots 5-8: relative per-bin heatmaps  |δA_bin / A_cos2φ,bin|
# ══════════════════════════════════════════════════════════════════════════════

def build_perbin_plot(cfg, luminosity, config_tag, config_label,
                      a_bag_midx, nf2_obj, base, outdir):
    """Build per-bin 1×3 subplot heatmap: |A_bin|, δA_bin, |δA_bin/A_bin|.

    Computes per-bin A_cos2φ,bin using scale = A_bag (mid_x shape) and
    per-bin α_s(Q²). Three subpanels: signal amplitude, absolute uncertainty,
    and relative uncertainty.

    Parameters
    ----------
    cfg           : BeamConfig
    luminosity    : float   integrated luminosity [fb⁻¹/nucleon]
    config_tag    : str     'mid' or 'top'
    config_label  : str     human-readable beam description
    a_bag_midx    : float   A_bag for mid_x shape (from solve_A_from_sum_rule)
    nf2_obj       : NuclearF2FromGrid  pre-constructed EPPS21 object
    base          : PartonF2           CT18NLO backend (for α_s)
    outdir        : pathlib.Path

    Returns
    -------
    pathlib.Path : output PNG path
    """
    L = luminosity
    lumi_tag = f"{int(L)}fb" if L == int(L) else f"{L:g}fb"

    # Run project_rates at target luminosity
    sc   = fom.Scenario(lumi_fb_per_nucleon=L, pol_ion_tensor=PZZ, q2_min=2.0)
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_obj)

    n_events = proj.n_events   # shape (nx, nq2)
    accepted = proj.accepted   # bool mask

    nf2 = proj.extras["nf2"]
    f1  = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2  = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y   = proj.extras["y"]

    # Per-bin A_cos2φ using scale = a_bag_midx (negative, so take abs for signal)
    scale_bag = a_bag_midx   # negative value
    delta_bag = delta_shape_with_alphas(proj.x, proj.q2, f1,
                                        scale=scale_bag,
                                        variant="mid_x",
                                        base=base)
    amp_bin = a_cos2phi(delta_bag, f1, f2, proj.x, y)

    # Validity mask: acceptance + minimum events + non-zero signal
    use = accepted & (n_events >= MIN_EVENTS)
    mask_valid = use & (np.abs(amp_bin) > 1e-10)

    # Build per-bin arrays (NaN outside valid mask)
    amp_abs  = np.full_like(n_events, np.nan, dtype=float)   # |A_bin|
    delta_a  = np.full_like(n_events, np.nan, dtype=float)   # δA_bin
    rel_unc  = np.full_like(n_events, np.nan, dtype=float)   # |δA_bin / A_bin|

    amp_abs[mask_valid]  = np.abs(amp_bin[mask_valid])
    da_vals = np.sqrt(2.0 / np.maximum(n_events[mask_valid], 1e-12)) / PZZ
    delta_a[mask_valid]  = da_vals
    rel_unc[mask_valid]  = da_vals / np.abs(amp_bin[mask_valid])

    # Bin edges for pcolormesh (log scale)
    x_edges  = proj.x_edges    # shape (nx+1,)
    q2_edges = proj.q2_edges   # shape (nq2+1,)

    # pcolormesh arrays: shape (nq2, nx) — transpose of (nx, nq2)
    X_edge_2d, Q2_edge_2d = np.meshgrid(x_edges, q2_edges)  # (nq2+1, nx+1)
    C_amp   = amp_abs.T    # shape (nq2, nx)
    C_delta = delta_a.T
    C_rel   = rel_unc.T

    def _lognorm(arr):
        """Return LogNorm from positive-min to max of arr (ignoring NaN)."""
        valid = arr[~np.isnan(arr) & np.isfinite(arr) & (arr > 0)]
        if valid.size == 0:
            return mcolors.LogNorm(vmin=1e-6, vmax=1.0)
        vmin = float(valid.min())
        vmax = float(valid.max())
        if vmin >= vmax:
            vmin = vmax / 100.0
        return mcolors.LogNorm(vmin=vmin, vmax=vmax)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # ── Left: |A_bin| ─────────────────────────────────────────────────────────
    ax = axes[0]
    pcm = ax.pcolormesh(X_edge_2d, Q2_edge_2d, C_amp,
                        norm=_lognorm(C_amp), cmap="viridis", shading="flat")
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label(r"$|A_{\cos 2\phi, \mathrm{bin}}|$", fontsize=10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$", fontsize=10)
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]", fontsize=10)
    ax.set_title("|A_bin| (signal)", fontsize=10)

    # ── Middle: δA_bin ────────────────────────────────────────────────────────
    ax = axes[1]
    pcm = ax.pcolormesh(X_edge_2d, Q2_edge_2d, C_delta,
                        norm=_lognorm(C_delta), cmap="plasma", shading="flat")
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label(r"$\delta A_{\cos 2\phi, \mathrm{bin}}$", fontsize=10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$", fontsize=10)
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]", fontsize=10)
    ax.set_title("δA_bin (uncertainty)", fontsize=10)

    # ── Right: |δA_bin / A_bin| ───────────────────────────────────────────────
    ax = axes[2]
    pcm = ax.pcolormesh(X_edge_2d, Q2_edge_2d, C_rel,
                        norm=_lognorm(C_rel), cmap="magma", shading="flat")
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label(
        r"$|\delta A_{\cos 2\phi, \mathrm{bin}} / A_{\cos 2\phi, \mathrm{bin}}|$",
        fontsize=10
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$", fontsize=10)
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]", fontsize=10)
    ax.set_title("|δA_bin / A_bin| (relative)", fontsize=10)

    fig.suptitle(
        f"$^6$Li per-bin decomposition — {config_label}, L = {L:g} fb$^{{-1}}$/u\n"
        f"Bag sum rule A = {a_bag_midx:.3f}, mid_x shape (α=0.7, β=3), "
        r"EPPS21, R1998, Cloet $P_{zz}$=0.267",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"money_delta_20260724_perbin_{config_tag}_{lumi_tag}.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")

    # ── Stdout diagnostics ────────────────────────────────────────────────────
    print(f"── perbin diagnostics: {config_label}, L={L:g} fb⁻¹/u ──")
    print(f"  |A_bin|:     min={np.nanmin(amp_abs):.2e}  "
          f"median={np.nanmedian(amp_abs):.2e}  max={np.nanmax(amp_abs):.2e}")
    print(f"  δA_bin:      min={np.nanmin(delta_a):.2e}  "
          f"median={np.nanmedian(delta_a):.2e}  max={np.nanmax(delta_a):.2e}")
    print(f"  |δA/A|_bin:  min={np.nanmin(rel_unc):.2e}  "
          f"median={np.nanmedian(rel_unc):.2e}  max={np.nanmax(rel_unc):.2e}")
    # Location of max |δA/A|
    flat_idx = np.nanargmax(rel_unc)
    ix, iq = np.unravel_index(flat_idx, rel_unc.shape)
    x_mid  = 0.5 * (x_edges[ix]  + x_edges[ix  + 1])  if ix  < len(x_edges)  - 1 else x_edges[ix]
    q2_mid = 0.5 * (q2_edges[iq] + q2_edges[iq + 1])  if iq < len(q2_edges) - 1 else q2_edges[iq]
    print(f"  Location of max |δA/A|:  x={x_mid:.2e},  Q²={q2_mid:.2f} GeV²")

    return outpath


# ══════════════════════════════════════════════════════════════════════════════
# K conversion factor: <Δ/F₁>_A = K * s
# ══════════════════════════════════════════════════════════════════════════════

def compute_K_factor(cfg, nf2_obj, base, variant):
    """Compute the rate-weighted K factor: <Δ/F₁>_A = K * s.

    Since Δ = s · α_s(Q²) · F₁ · x^α (1-x)^β / _PEAK_VALS[variant], we have:

        <Δ/F₁>_A = s · [Σ_bins N_bin · α_s(Q²) · x^α (1-x)^β / _PEAK_VALS]
                       / [Σ_bins N_bin · F₁]  × F₁   (F₁ cancels per bin)

    More precisely:
        K = Σ_bins [N_bin · α_s(Q²_bin) · x_bin^α · (1-x_bin)^β / _PEAK_VALS]
            / Σ_bins [N_bin · F₁(x_bin, Q²_bin) / nucleon]
        × F₁_ref

    But the rate-weighted mean <Δ/F₁> at scale s=1 is simply:
        K = Σ_bins [N_bin · (Δ_bin / F₁_bin) / s]
          = Σ_bins [N_bin · α_s(Q²_bin) · x_bin^α (1-x_bin)^β / _PEAK_VALS]
            / Σ_bins N_bin

    (The N_bin weights are rates, so the denominator is the total accepted rate;
    F₁ in Δ/F₁ cancels per bin, leaving only the shape and α_s factors.)

    Parameters
    ----------
    cfg      : BeamConfig
    nf2_obj  : NuclearF2FromGrid  pre-constructed EPPS21 object
    base     : PartonF2           CT18NLO backend (for α_s)
    variant  : str                shape variant ('low_x', 'mid_x', 'high_x')

    Returns
    -------
    float : K such that rate-weighted <Δ/F₁> = K * s
    """
    sc   = fom.Scenario(lumi_fb_per_nucleon=1.0, pol_ion_tensor=PZZ, q2_min=2.0)
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_obj)

    n_events = proj.n_events   # shape (nx, nq2)
    accepted = proj.accepted & (proj.n_events >= MIN_EVENTS)

    n_acc  = n_events[accepted]
    x_acc  = proj.x[accepted]
    q2_acc = proj.q2[accepted]

    total_n = n_acc.sum()
    if total_n <= 0:
        return np.nan

    alpha_v, beta_v = _VARIANTS[variant]
    as_acc = alpha_s(q2_acc, base=base)
    shape_acc = (
        as_acc
        * np.power(np.maximum(x_acc, 1e-12), alpha_v)
        * np.power(np.maximum(1.0 - x_acc, 0.0), beta_v)
        / _PEAK_VALS[variant]
    )

    K = float((n_acc * shape_acc).sum() / total_n)
    return K


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description=(
            "6Li sum-rule-constrained A_cos2phi reach — bag-model & lattice-QCD "
            "references. Produces 8 PNGs:\n"
            "  Plots 1-4: δA/A vs scale and L_5σ vs scale with sum-rule vertical "
            "lines at |A_bag| and |A_lat|.\n"
            "  Plots 5-8: relative per-bin |δA/A_cos2φ| heatmaps for (MID, TOP) × "
            "(10, 100 fb⁻¹/u) using bag sum-rule Δ.\n\n"
            "Physics: Δ(x,Q²) = scale · α_s(Q²) · F₁(x,Q²) · x^α (1-x)^β.\n"
            "Sum rule: ∫x·Δ dx = c·α_s(Q²) → A = c / ∫ x F₁ shape dx "
            "(α_s cancels).\n"
            "  Bag:     c = -0.012\n"
            "  Lattice: c = -0.009\n\n"
            "Estimated runtime: ~15-25 min (9 sig² evaluations × 3 configs "
            "+ 3 sum-rule solves × 3 variants × 3 configs "
            "+ 4 per-bin heatmap evals)."
        )
    )
    ap.add_argument(
        "--outdir",
        default="fastsim/out/money_delta",
        help="Output directory for PNGs (default: fastsim/out/money_delta)",
    )
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)

    # ── Startup banner ────────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("money_delta_20260724.py — sum-rule-constrained 6Li tensor reach")
    print("=" * 72)
    print()
    print("Observable: Δ(x, Q²) = scale · α_s(Q²) · F₁(x, Q²) · x^α (1-x)^β")
    print("Sum rule:   ∫ x · Δ dx = c · α_s(Q²)")
    print("  Bag:      c = -0.012")
    print("  Lattice:  c = -0.009")
    print("Solving for A = scale such that sum rule is satisfied.")
    print("  → α_s cancels on both sides; A = c / ∫ x F₁(x, <Q²>) shape dx")
    print()

    # ── Load CT18NLO backend ──────────────────────────────────────────────────
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

    # Initialize α_s source (triggers the one-line print)
    _ = alpha_s(10.0, base=base)

    # ── Verify nuclear PDF set is installed ───────────────────────────────────
    try:
        from parton import mkPDF as _mkPDF
        _p = _mkPDF(EPPS21_SET, 0)
    except Exception as exc:
        print(
            f"ERROR: nuclear PDF set '{EPPS21_SET}' could not be loaded.\n"
            "Install with:\n"
            f"  python3 -m parton install {EPPS21_SET}\n"
            f"Underlying error: {exc}",
            file=sys.stderr,
        )
        sys.exit(2)

    # ── Three EIC beam configurations for ⁶Li ─────────────────────────────
    all_configs = [
        (
            BeamConfig(electron_energy=5.0, ion=LI6, ion_momentum_per_nucleon=27.5),
            "low",
            r"5 GeV $e$ × 27.5 GeV/u $^6$Li",
        ),
        (
            BeamConfig(electron_energy=10.0, ion=LI6, ion_momentum_per_nucleon=50.0),
            "mid",
            r"10 GeV $e$ × 50 GeV/u $^6$Li",
        ),
        (
            BeamConfig(electron_energy=18.0, ion=LI6, ion_momentum_per_nucleon=137.5),
            "top",
            r"18 GeV $e$ × 137.5 GeV/u $^6$Li",
        ),
    ]

    # ── Step 1: Solve sum-rule A values for all (config, variant, model) ──────
    print("=" * 72)
    print("Step 1: Solving sum-rule A values")
    print("=" * 72)

    a_vals = {}   # (config_tag, variant, model) → float

    with r_override(r1998):
        for cfg, cfg_tag, cfg_label in all_configs:
            print(f"\n  Config: {cfg_tag.upper()} — {cfg_label}")
            nf2_obj = NuclearF2FromGrid(cfg.ion, EPPS21_SET)

            for variant in VARIANT_NAMES:
                A_bag, q2_mean = solve_A_from_sum_rule(
                    cfg, nf2_obj, base, variant, C_BAG
                )
                A_lat, _ = solve_A_from_sum_rule(
                    cfg, nf2_obj, base, variant, C_LAT
                )
                a_vals[(cfg_tag, variant, "bag")] = A_bag
                a_vals[(cfg_tag, variant, "lat")] = A_lat
                print(
                    f"    variant={variant:8s}  <Q²>={q2_mean:.2f} GeV²  "
                    f"A_bag={A_bag:.4e}  A_lat={A_lat:.4e}"
                )

    # ── Summary table of A values ─────────────────────────────────────────────
    print()
    print("=" * 72)
    print("A-value summary table (scale satisfying sum rule)")
    print("=" * 72)
    header = f"{'Config':<8} {'Variant':<10} {'A_bag':>12} {'A_lat':>12} {'|A_bag|':>12} {'|A_lat|':>12}"
    print(header)
    print("-" * 68)
    for cfg, cfg_tag, _ in all_configs:
        for variant in VARIANT_NAMES:
            A_bag = a_vals[(cfg_tag, variant, "bag")]
            A_lat = a_vals[(cfg_tag, variant, "lat")]
            print(
                f"  {cfg_tag:<6} {variant:<10} "
                f"{A_bag:>12.4e} {A_lat:>12.4e} "
                f"{abs(A_bag):>12.4e} {abs(A_lat):>12.4e}"
            )
    print()

    # ── K conversion factors: <Δ/F₁>_A = K * s ───────────────────────────────
    print("=" * 72)
    print("K conversion factors: <Δ/F₁>_A  =  K * s   (rate-weighted, accepted bins)")
    print("=" * 72)
    header_k = (
        f"{'Config':<8} {'Variant':<12} {'K':>12} "
        f"{'<Δ/F₁> @ s0^bag':>18} {'<Δ/F₁> @ s0^lat':>18}"
    )
    print(header_k)
    print("-" * 72)

    k_vals = {}   # (config_tag, variant) → float

    with r_override(r1998):
        for cfg, cfg_tag, _ in all_configs:
            nf2_obj = NuclearF2FromGrid(cfg.ion, EPPS21_SET)
            for variant in VARIANT_NAMES:
                K = compute_K_factor(cfg, nf2_obj, base, variant)
                k_vals[(cfg_tag, variant)] = K
                a_bag_k = a_vals.get((cfg_tag, variant, "bag"), np.nan)
                a_lat_k = a_vals.get((cfg_tag, variant, "lat"), np.nan)
                df1_bag = K * abs(a_bag_k) if np.isfinite(K) and np.isfinite(a_bag_k) else np.nan
                df1_lat = K * abs(a_lat_k) if np.isfinite(K) and np.isfinite(a_lat_k) else np.nan
                print(
                    f"  {cfg_tag:<6} {variant:<12} {K:>12.3e} "
                    f"{df1_bag:>18.3e} {df1_lat:>18.3e}"
                )

    print()
    print("Note: on Plots 1-4, x-axis is raw scale s (universal across all curves).")
    print("      For physical interpretation, multiply s by K to get rate-weighted <Δ/F₁>.")
    print()

    # ── Step 2: Compute sig² per fb⁻¹ at S0 for each (config, variant) ───────
    print("=" * 72)
    print("Step 2: Computing sig² per fb⁻¹/u at S0=1e-3 (9 evaluations)")
    print("=" * 72)

    sig2_dict = {}   # (config_tag, variant) → float

    with r_override(r1998):
        n_total = len(all_configs) * len(VARIANT_NAMES)
        done = 0
        for cfg, cfg_tag, cfg_label in all_configs:
            nf2_obj = NuclearF2FromGrid(cfg.ion, EPPS21_SET)
            for variant in VARIANT_NAMES:
                done += 1
                print(
                    f"  [{done}/{n_total}] config={cfg_tag}, variant={variant} …",
                    end=" ", flush=True,
                )
                s2 = sig2_per_fb_at_sumrule(
                    cfg, S0, PZZ, nf2_obj, base, variant=variant,
                )
                sig2_dict[(cfg_tag, variant)] = s2
                l5 = 25.0 / max(s2, 1e-30)
                print(f"L_5σ(1e-3) = {l5:.2f} fb⁻¹/u")

    # ── Step 3: Plots 1-4 ─────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("Step 3: Building Plots 1-4")
    print("=" * 72)

    print("\n  → Plot 1: δA/A vs scale, 9 curves, two single-panel PNGs (10 and 100 fb⁻¹/u)")
    p1_paths = build_plot1(sig2_dict, a_vals, outdir)
    for p in p1_paths:
        print(f"    wrote {p}")

    print("  → Plot 2: L_5σ vs scale, 9 curves, single panel")
    build_plot2(sig2_dict, a_vals, outdir)

    print("  → Plot 3: bag vs lattice, mid config, δA/A, single panel (both luminosities)")
    sig2_mid_midx = sig2_dict.get(("mid", "mid_x"), 0.0)
    build_plot3(sig2_mid_midx, a_vals, outdir)

    print("  → Plot 4: bag vs lattice, mid config, L_5σ, single panel")
    build_plot4(sig2_mid_midx, a_vals, outdir)

    # ── Step 4: Plots 5-8 — relative per-bin heatmaps (MID and TOP) ──────────
    print()
    print("=" * 72)
    print("Step 4: Building Plots 5-8 (per-bin heatmaps, MID and TOP)")
    print("=" * 72)

    perbin_configs = [c for c in all_configs if c[1] in ("mid", "top")]

    with r_override(r1998):
        for cfg, cfg_tag, cfg_label in perbin_configs:
            nf2_obj = NuclearF2FromGrid(cfg.ion, EPPS21_SET)
            a_bag_midx = a_vals.get((cfg_tag, "mid_x", "bag"), None)
            if a_bag_midx is None:
                print(f"  WARNING: A_bag for {cfg_tag}/mid_x not found; skipping heatmaps.")
                continue
            for L in LUMINOSITIES_FB:
                lumi_tag = f"{int(L)}fb" if L == int(L) else f"{L:g}fb"
                print(
                    f"  → Plot perbin {cfg_tag} {lumi_tag}: "
                    f"bag A={a_bag_midx:.4e}, L={L:g} fb⁻¹/u …"
                )
                build_perbin_plot(
                    cfg, L, cfg_tag, cfg_label,
                    a_bag_midx, nf2_obj, base, outdir,
                )

    # ── Final summary ─────────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("FINAL SUMMARY — sig² and reach at bag/lattice A values")
    print("=" * 72)
    print()
    print(f"{'Config':<8} {'Variant':<10} {'A_bag':>12} {'sig²(A_bag)':>14} "
          f"{'δA/A @10fb':>12} {'δA/A @100fb':>12} {'L_5σ(A_bag)':>12}")
    print("-" * 82)

    with r_override(r1998):
        for cfg, cfg_tag, cfg_label in all_configs:
            nf2_obj = NuclearF2FromGrid(cfg.ion, EPPS21_SET)
            for variant in VARIANT_NAMES:
                A_bag = a_vals.get((cfg_tag, variant, "bag"))
                A_lat = a_vals.get((cfg_tag, variant, "lat"))
                s2_s0 = sig2_dict.get((cfg_tag, variant), 0.0)
                if A_bag is None or s2_s0 <= 0:
                    continue
                # Analytic rescaling: sig²(A) = s2_s0 × (A/S0)²
                s2_at_bag = s2_s0 * (abs(A_bag) / S0) ** 2
                s2_at_lat = s2_s0 * (abs(A_lat) / S0) ** 2
                # δA/A at fixed luminosity
                rel_10  = 1.0 / np.sqrt(max(s2_at_bag * 10.0,  1e-60))
                rel_100 = 1.0 / np.sqrt(max(s2_at_bag * 100.0, 1e-60))
                # L_5σ at A_bag
                l5_bag = 25.0 / max(s2_at_bag, 1e-30)
                print(
                    f"  {cfg_tag:<6} {variant:<10} {A_bag:>12.4e} "
                    f"{s2_at_bag:>14.4e} "
                    f"{rel_10:>12.4f} {rel_100:>12.4f} {l5_bag:>12.2f}"
                )

    print()
    print("Note: sign(A) negative for both bag and lattice; |A| used for")
    print("      vertical-line positions in plots. sig² ∝ A² so sign cancels.")
    print("      Central = mid_x shape, R1998, EPPS21, Cloet P_zz=0.267.")
    print("      δA/A = 1/√(sig²×L); L_5σ = 25/sig² [fb⁻¹/u].")
    print()
    print("Output directory:", outdir)
    print("Done.")


if __name__ == "__main__":
    main()
