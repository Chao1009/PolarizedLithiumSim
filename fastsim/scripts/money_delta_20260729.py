#!/usr/bin/env python3
"""6Li φ-distribution plots — detector-realistic version with MC smearing.

Detector-realistic extension of money_delta_20260728.py. Produces the same
16 plots (3 heatmaps + 8 φ plots + 5 TOP+peakbin bin-width scan), but with:

  1. Monte Carlo smearing of the scattered electron's (E', θ, φ) using
     ePIC tracking resolutions (piecewise in η).
  2. Bin migration in (x, Q²) from reconstructed kinematics.
  3. η-dependent electron-ID efficiency applied per event using
     ATHENA+ECCE-anchored numbers.

Reco-analysis mask
──────────────────
All reconstructed-level outputs (peak-bin find, Q²-slice find, integrated
averages, per-bin heatmaps, φ-modulation plots, smearing diagnostics) use a
single, explicitly-defined **reconstructed-analysis mask**:

    reco_mask = proj.accepted & (n_reco >= MIN_EVENTS)

`proj.accepted` encodes the FULL detector acceptance as built by
`fom.project_rates`: DIS kinematic cuts (q2_min, y_min, y_max, w2_min) PLUS
detector acceptance (eta_min ≤ η ≤ eta_max, E'_e ≥ e_prime_min).  Bin-center
kinematics are shared between the true and reco grids, so reapplying the
acceptance at bin centers reduces to a logical AND with `proj.accepted`.  The
`n_reco >= MIN_EVENTS` factor is a reco-side statistics floor with no
truth-side analogue.

This mask is computed once per config in `main()` via `reco_analysis_mask()`
and passed explicitly into every reco-side helper.  The true-level
`proj.accepted` is never used as a reco selection; it is used only where
truth-level acceptance is genuinely needed (smearing loop source bins,
true-level comparison A_true).

The invariant `reco_mask ⊆ proj.accepted` is asserted in `main()`.  See
`fastsim/notes/money_delta_note_2026-07-29_fix.md` for the full derivation and
numerical impact summary.

Physics setup (identical to money_delta_20260728.py):
  L   = 10 fb⁻¹/nucleon
  Shape: mid_x (α=0.7, β=3)
  R:    R1998
  PDF:  EPPS21nlo_CT18Anlo_Li6
  P_zz: 0.267  (Cloet convention)
  Interp. A:  Δ = s · α_s(Q²) · F₁ · x^α · (1-x)^β

A_bag values (hardcoded from money_delta_20260724.py mid_x output):
  LOW (5 × 27.5 GeV/u): |A_bag| = 0.318
  MID (10 × 50 GeV/u):  |A_bag| = 0.310
  TOP (18 × 137.5 GeV/u): |A_bag| = 0.297

Detector model:
  Tracking momentum resolution: σ_p/p = √(a²·p² + b²), piecewise in η.
  Tracking angular resolution: σ_θ = σ_φ (mrad), piecewise in η.
  Electron-ID efficiency: η-dependent, linearly interpolated from 9 anchor
    points sourced from:
    ATHENA (JINST 17 (2022) P10019, arXiv:2210.09048, Table 5) +
    ECCE (NIM A 1055 (2023) 168464, arXiv:2207.09437, Sec. 3.5.2).
    No official ε_eID(η) curve is published for ePIC as of pCDR v1 (Dec 2024).

ECal decision:
  ECal resolution parameters are not used in this script; the user chose
  tracking-only for the scattered electron across all η regions.
  See tracking_ab() and tracking_ang_mrad() for the implemented model.
  Alternative: ECal S/C resolution could be used in the forward (electron-
  going, η < −1) region where tracking degrades. That was provided in
  parent communication but not implemented here per user decision.

Output (16 PNGs, prefixed _20260729_):
  Per-bin (x, Q²) heatmaps — motivation for the φ plots (3 plots):
    money_delta_20260729_perbin_low.png
    money_delta_20260729_perbin_mid.png
    money_delta_20260729_perbin_top.png
  LOW+MID cases 1-3 and TOP cases 2-3 (8 plots, 5° bins):
    money_delta_20260729_phimodulation_peakbin_low.png
    money_delta_20260729_phimodulation_peakbin_mid.png
    money_delta_20260729_phimodulation_q2slice_low.png
    money_delta_20260729_phimodulation_q2slice_mid.png
    money_delta_20260729_phimodulation_q2slice_top.png
    money_delta_20260729_phimodulation_integrated_low.png
    money_delta_20260729_phimodulation_integrated_mid.png
    money_delta_20260729_phimodulation_integrated_top.png
  TOP case 1 φ-bin-width scan (5 plots):
    money_delta_20260729_phimodulation_peakbin_top_5deg.png
    money_delta_20260729_phimodulation_peakbin_top_10deg.png
    money_delta_20260729_phimodulation_peakbin_top_20deg.png
    money_delta_20260729_phimodulation_peakbin_top_30deg.png
    money_delta_20260729_phimodulation_peakbin_top_45deg.png

Reproducibility: rng = np.random.default_rng(seed=42) inside smear_config().
Runtime estimate: ~5–10 min on a laptop (≈10⁶ MC events total).
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
import polli_fastsim.kinematics as kinematics_mod

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# ══════════════════════════════════════════════════════════════════════════════
# NuclearF2FromGrid — local class; do NOT modify polli_fastsim/structure.py
# (copied verbatim from money_delta_20260728.py)
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
# Δ shape variants (mid_x only for this script)
# ══════════════════════════════════════════════════════════════════════════════

_VARIANTS = {
    "low_x":  (0.3, 4.0),
    "mid_x":  (0.7, 3.0),
    "high_x": (1.5, 2.0),
}


def _xab_peak_value(alpha, beta):
    """Peak value of x^alpha * (1-x)^beta at x = alpha/(alpha+beta)."""
    xp = alpha / (alpha + beta)
    return (xp ** alpha) * ((1.0 - xp) ** beta)


_PEAK_VALS = {name: _xab_peak_value(a, b) for name, (a, b) in _VARIANTS.items()}


# ══════════════════════════════════════════════════════════════════════════════
# R = σ_L/σ_T parameterization (R1998 only)
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

# Fixed integrated luminosity [fb⁻¹/nucleon]
LUMI_FB = 10.0

# Shape variant for all plots
VARIANT = "mid_x"

# Nuclear PDF set name (fixed throughout)
EPPS21_SET = "EPPS21nlo_CT18Anlo_Li6"

# Minimum events cut
MIN_EVENTS = 10


def reco_analysis_mask(proj, n_reco, min_events=MIN_EVENTS):
    """Reconstructed-level analysis mask (shape (nx, nq2)).

    Definition:
        reco_mask = proj.accepted & (n_reco >= min_events)

    This intentionally reuses the FULL acceptance that
    `fom.project_rates` bakes into `proj.accepted` — namely the DIS
    kinematic cuts (q2_min, y_min, y_max, w2_min) AND the detector
    acceptance (eta_min <= eta <= eta_max, E'_e >= e_prime_min).
    Bin-center kinematics are identical for the true and reco grids,
    so reapplying the acceptance test at bin centers reduces to a
    logical AND with `proj.accepted`. The extra factor of
    `n_reco >= min_events` is a reco-side statistics floor with no
    truth-side analogue.

    Note: `proj.accepted` alone is the TRUE-level acceptance; do NOT
    apply it to reco arrays without also enforcing the reco statistics
    floor. This function is the single source of truth for reco selection.

    Invariant: reco_mask ⊆ proj.accepted (asserted by callers if desired).
    """
    return proj.accepted & (np.asarray(n_reco) >= float(min_events))


# Hardcoded A_bag values (abs) from money_delta_20260724.py mid_x output
A_BAG = {
    "low": 0.318,  # 5 × 27.5 GeV/u
    "mid": 0.310,  # 10 × 50 GeV/u
    "top": 0.297,  # 18 × 137.5 GeV/u
}

# Number of φ bins (5° each, full 2π)
N_PHI = 72


# ══════════════════════════════════════════════════════════════════════════════
# α_s(Q²) helper — primary: parton PDF table; fallback: LO analytic
# (copied verbatim from money_delta_20260728.py)
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
    """
    if _ALPHAS_SOURCE[0] is None:
        source, func = _build_alphas_interpolator(base)
        _ALPHAS_SOURCE[0] = (source, func)
        print(f"[alpha_s] Using source: {source}")
    _, func = _ALPHAS_SOURCE[0]
    result = func(np.asarray(q2, dtype=float))
    if result.ndim == 0:
        return float(result)
    return result


# ══════════════════════════════════════════════════════════════════════════════
# Δ with α_s(Q²) prefactor (copied verbatim from money_delta_20260728.py)
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
    scale    : float       dimensionless amplitude
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
# sig² helper (kept for API compatibility with parent)
# ══════════════════════════════════════════════════════════════════════════════

def sig2_per_fb_at_sumrule(cfg, scale, pzz, nf2_obj, base, variant="mid_x",
                            min_events=10):
    """Sig² per fb⁻¹/nucleon using Δ = scale · α_s(Q²) · F₁ · x^α · (1-x)^β.

    The R-override must be active in the calling context.
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
# ── DETECTOR MODEL HELPERS (new in 20260729) ──────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════

def tracking_ab(eta):
    """Return (a, b) for tracking momentum resolution σ_p/p = √(a²·p² + b²).

    Piecewise in pseudorapidity η, vectorised over numpy arrays.

    Parameters
    ----------
    eta : array-like
        Pseudorapidity of the scattered electron track.

    Returns
    -------
    a : ndarray  (same shape as eta)
    b : ndarray  (same shape as eta)

    η region                  a         b
    η < -3.5                  0.0020    0.030
    -3.5 ≤ η < -2.5          0.0010    0.030
    -2.5 ≤ η < -1.0          0.0005    0.010
    -1.0 ≤ η <  1.0 (barrel) 0.0005    0.005
     1.0 ≤ η <  2.5          0.0005    0.010
     2.5 ≤ η <  3.5          0.0010    0.025
    η ≥ 3.5                  0.0010    0.025
    """
    eta = np.asarray(eta, dtype=float)
    a = np.full_like(eta, 0.0005)
    b = np.full_like(eta, 0.005)

    # η < -3.5
    m = eta < -3.5
    a[m] = 0.0020
    b[m] = 0.030

    # -3.5 ≤ η < -2.5
    m = (eta >= -3.5) & (eta < -2.5)
    a[m] = 0.0010
    b[m] = 0.030

    # -2.5 ≤ η < -1.0
    m = (eta >= -2.5) & (eta < -1.0)
    a[m] = 0.0005
    b[m] = 0.010

    # -1.0 ≤ η < 1.0  (barrel) — already set as default (a=0.0005, b=0.005)

    # 1.0 ≤ η < 2.5
    m = (eta >= 1.0) & (eta < 2.5)
    a[m] = 0.0005
    b[m] = 0.010

    # 2.5 ≤ η < 3.5
    m = (eta >= 2.5) & (eta < 3.5)
    a[m] = 0.0010
    b[m] = 0.025

    # η ≥ 3.5
    m = eta >= 3.5
    a[m] = 0.0010
    b[m] = 0.025

    return a, b


def tracking_ang_mrad(eta):
    """Return angular resolution in mrad (σ_θ = σ_φ = result × 1e-3 rad).

    Piecewise in pseudorapidity η, vectorised over numpy arrays.

    Parameters
    ----------
    eta : array-like
        Pseudorapidity of the scattered electron track.

    Returns
    -------
    mrad : ndarray  (same shape as eta)
        Angular resolution in milliradians.

    η region              mrad
    η < -3.5               5
    -3.5 ≤ η < -2.5       3
    -2.5 ≤ η < -1.0       2
    -1.0 ≤ η <  1.0        1 (barrel)
     1.0 ≤ η <  2.5       2
     2.5 ≤ η <  3.5       3
    η ≥ 3.5                5
    """
    eta = np.asarray(eta, dtype=float)
    mrad = np.full_like(eta, 1.0)  # barrel default

    # η < -3.5
    m = eta < -3.5
    mrad[m] = 5.0

    # -3.5 ≤ η < -2.5
    m = (eta >= -3.5) & (eta < -2.5)
    mrad[m] = 3.0

    # -2.5 ≤ η < -1.0
    m = (eta >= -2.5) & (eta < -1.0)
    mrad[m] = 2.0

    # -1.0 ≤ η < 1.0 (barrel) — already default (1.0 mrad)

    # 1.0 ≤ η < 2.5
    m = (eta >= 1.0) & (eta < 2.5)
    mrad[m] = 2.0

    # 2.5 ≤ η < 3.5
    m = (eta >= 2.5) & (eta < 3.5)
    mrad[m] = 3.0

    # η ≥ 3.5
    m = eta >= 3.5
    mrad[m] = 5.0

    return mrad


def eps_eid(eta):
    """Electron-ID efficiency ε_eID(η), linearly interpolated from 9 anchor points.

    Source: ATHENA (JINST 17 (2022) P10019, arXiv:2210.09048, Table 5) +
            ECCE (NIM A 1055 (2023) 168464, arXiv:2207.09437, Sec. 3.5.2).
    No official ε_eID(η) curve is published for ePIC as of pCDR v1 (Dec 2024).

    η anchor points:
    η       ε_eID
    -3.5    0.85
    -3.0    0.92
    -2.0    0.95
    -1.0    0.93
     0.0    0.90
     1.0    0.90
     2.0    0.85
     3.0    0.80
     3.5    0.70

    Outside |η| > 3.5: ε_eID = 0.

    Parameters
    ----------
    eta : array-like
        Pseudorapidity (scalar or array).

    Returns
    -------
    ndarray : ε_eID values in [0, 1].
    """
    eta = np.asarray(eta, dtype=float)
    eta_pts = np.array([-3.5, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 3.5])
    eps_pts = np.array([0.85,  0.92,  0.95,  0.93, 0.90, 0.90, 0.85, 0.80, 0.70])
    return np.where(
        (eta < -3.5) | (eta > 3.5),
        0.0,
        np.interp(eta, eta_pts, eps_pts),
    )


# ══════════════════════════════════════════════════════════════════════════════
# smear_config — MC smearing for one beam configuration
# ══════════════════════════════════════════════════════════════════════════════

def smear_config(cfg, proj, base, nf2, A_bag_config, rng, n_mc=1000):
    """Run MC smearing for one beam configuration.

    For each accepted (x_true, Q²_true) bin, draw n_mc MC pseudo-events,
    smear (E', θ, φ) using ePIC tracking resolutions, apply η-dependent
    electron-ID efficiency, recompute reconstructed kinematics (x_reco,
    Q²_reco), and accumulate event counts + rate-weighted A_cos2φ into the
    reconstructed grid.

    Parameters
    ----------
    cfg          : BeamConfig
    proj         : BinnedProjection   output of project_rates at target L
    base         : PartonF2           CT18NLO backend (for α_s + Δ)
    nf2          : NuclearF2FromGrid
    A_bag_config : float              |A_bag| for this config
    rng          : np.random.Generator
        Seeded generator (use np.random.default_rng(seed=42) for reproducibility).
    n_mc         : int               number of MC pseudo-events per true bin

    Returns
    -------
    n_reco : ndarray, shape (nx, nq2)
        Reconstructed event count per bin (weighted, with efficiency applied).
    A_reco : ndarray, shape (nx, nq2)
        Rate-weighted A_cos2φ per reconstructed bin.
    diag   : dict
        Diagnostic statistics:
          'n_true_total'   : total events before smearing (float)
          'n_reco_total'   : total events after smearing + efficiency (float)
          'survival_frac'  : n_reco_total / n_true_total
          'same_bin_frac'  : fraction of weight landing in the original bin
          'n_rejected'     : total MC events rejected (physicality cuts)
    """
    s          = proj.extras["s"]
    E_e        = cfg.electron_energy
    p_ion      = cfg.ion_momentum_per_nucleon  # GeV/c per nucleon
    x_edges    = proj.x_edges     # shape (nx+1,)
    q2_edges   = proj.q2_edges    # shape (nq2+1,)
    nx         = len(x_edges) - 1
    nq2        = len(q2_edges) - 1

    # ── Per-bin true A_cos2φ (needed for the rate-weighted A_reco accumulation)
    x_2d    = proj.x    # shape (nx, nq2)
    q2_2d   = proj.q2
    y_2d    = proj.extras["y"]
    n_true  = proj.n_events  # shape (nx, nq2)
    acc     = proj.accepted  # bool mask

    f1_2d = nf2.f1a(x_2d, q2_2d) / cfg.ion.A
    f2_2d = nf2.f2a(x_2d, q2_2d) / cfg.ion.A

    delta_2d = delta_shape_with_alphas(
        x_2d, q2_2d, f1_2d, scale=A_bag_config, variant=VARIANT, base=base
    )
    # Compute per-bin true A_cos2φ
    A_true_2d = np.zeros_like(n_true, dtype=float)
    for ix in range(nx):
        for iq2 in range(nq2):
            if acc[ix, iq2]:
                A_true_2d[ix, iq2] = float(a_cos2phi(
                    delta_2d[ix, iq2], f1_2d[ix, iq2], f2_2d[ix, iq2],
                    x_2d[ix, iq2], y_2d[ix, iq2]
                ))

    # ── Reco accumulators
    n_reco      = np.zeros((nx, nq2), dtype=float)
    A_reco_num  = np.zeros((nx, nq2), dtype=float)  # numerator: Σ(A_true * w)
    A_reco_den  = np.zeros((nx, nq2), dtype=float)  # denominator: Σ(w)

    # ── Diagnostic counters
    n_true_total      = 0.0
    n_reco_total      = 0.0
    same_bin_sum      = 0.0
    n_rejected        = 0
    raw_same_bin_count = 0     # number of valid MC events landing in true bin (unweighted)
    raw_valid_count    = 0     # number of valid MC events total (unweighted)

    # ── Loop over accepted true bins ──────────────────────────────────────────
    for ix in range(nx):
        for iq2 in range(nq2):
            if not acc[ix, iq2]:
                continue

            n_bin_true = float(n_true[ix, iq2])
            if n_bin_true <= 0.0:
                continue
            n_true_total += n_bin_true

            x_c   = float(x_2d[ix, iq2])
            q2_c  = float(q2_2d[ix, iq2])
            y_c   = float(y_2d[ix, iq2])

            # True kinematics of the scattered electron at bin center
            e_prime_true, theta_true, eta_true = kinematics_mod.scattered_electron(
                x_c, y_c, s, E_e
            )
            e_prime_true = float(e_prime_true)
            theta_true   = float(theta_true)
            eta_true     = float(eta_true)

            # Treat p_true ≈ E'_true (electron mass negligible at DIS energies)
            p_true = e_prime_true

            # Get tracking resolutions at true η (scalar → 0-d arrays → floats)
            a_val, b_val = tracking_ab(np.array([eta_true]))
            a_val = float(a_val[0])
            b_val = float(b_val[0])
            mrad_val = float(tracking_ang_mrad(np.array([eta_true]))[0])

            sigma_p = p_true * np.sqrt(a_val**2 * p_true**2 + b_val**2)
            sigma_ang = mrad_val * 1.0e-3   # convert mrad → rad

            # Per-MC-event weight before efficiency
            w_base = n_bin_true / n_mc

            # True A_cos2φ for this bin (used in A_reco accumulation)
            A_true_bin = A_true_2d[ix, iq2]

            # ── Draw n_mc pseudo-events ───────────────────────────────────────
            # φ_true: uniformly on [0, 2π] (integrated, so φ direction is arbitrary)
            phi_true_arr = rng.uniform(0.0, 2.0 * np.pi, size=n_mc)

            # Smear momentum and angles
            p_reco_arr     = rng.normal(p_true, sigma_p,   size=n_mc)
            theta_reco_arr = rng.normal(theta_true, sigma_ang, size=n_mc)
            phi_reco_arr   = rng.normal(phi_true_arr, sigma_ang)  # same σ per event

            # Recompute η_reco
            theta_reco_safe = np.clip(theta_reco_arr, 1.0e-9, np.pi - 1.0e-9)
            tan_half = np.tan(theta_reco_safe / 2.0)
            tan_half = np.maximum(tan_half, 1.0e-12)
            eta_reco_arr = -np.log(tan_half)

            # E'_reco ≈ p_reco (massless)
            e_prime_reco_arr = p_reco_arr

            # Reconstructed Q², y, x
            # Convention matches polli_fastsim.kinematics.scattered_electron:
            # ion beam along +z, electron beam along -z, θ measured from +z.
            # Scattered electron sits at θ ≈ π (negative η).
            # Formulas: Q² = 2 E_e E' (1 + cos θ), y = 1 - (E'/2E_e)(1 - cos θ)
            theta_phys = np.clip(theta_reco_arr, 1.0e-9, np.pi - 1.0e-9)

            q2_reco_arr = (
                2.0 * E_e * e_prime_reco_arr
                * (1.0 + np.cos(theta_phys))
            )

            y_reco_arr = (
                1.0
                - (e_prime_reco_arr / (2.0 * E_e))
                * (1.0 - np.cos(theta_phys))
            )
            # Guard division by (y_reco * s)
            denom_arr = y_reco_arr * s
            with np.errstate(invalid="ignore", divide="ignore"):
                x_reco_arr = np.where(
                    np.abs(denom_arr) > 0.0,
                    q2_reco_arr / denom_arr,
                    -1.0,   # sentinel → will be rejected
                )

            # ── Physicality cuts ──────────────────────────────────────────────
            valid = (
                (x_reco_arr > 0.0) & (x_reco_arr < 1.0)
                & (q2_reco_arr > 0.0)
                & (y_reco_arr > 0.0) & (y_reco_arr < 1.0)
            )
            n_rejected += int(np.sum(~valid))

            n_valid_this_bin = int(np.sum(valid))
            raw_valid_count += n_valid_this_bin

            if not np.any(valid):
                continue

            # Select valid events
            x_reco_v   = x_reco_arr[valid]
            q2_reco_v  = q2_reco_arr[valid]
            eta_reco_v = eta_reco_arr[valid]

            # Electron-ID efficiency per event
            eid_v = eps_eid(eta_reco_v)

            # Weight per valid MC event
            w_v = w_base * eid_v  # shape (n_valid,)

            # ── Find reconstructed (x, Q²) bin via digitize ──────────────────
            # np.digitize returns 1-based index; 0 = below first edge,
            # len(edges) = above last edge; valid bins are 1..len(edges)-1
            ix_reco_arr  = np.digitize(x_reco_v,  x_edges)  - 1
            iq2_reco_arr = np.digitize(q2_reco_v, q2_edges) - 1

            # Clip to valid grid range [0, n-1]
            in_grid = (
                (ix_reco_arr  >= 0) & (ix_reco_arr  < nx)
                & (iq2_reco_arr >= 0) & (iq2_reco_arr < nq2)
            )

            for k in range(len(x_reco_v)):
                if not in_grid[k]:
                    continue
                ir  = ix_reco_arr[k]
                jr  = iq2_reco_arr[k]
                wk  = w_v[k]

                n_reco[ir, jr]     += wk
                A_reco_num[ir, jr] += A_true_bin * wk
                A_reco_den[ir, jr] += wk
                n_reco_total       += wk

                if ir == ix and jr == iq2:
                    same_bin_sum += wk           # weighted (existing)
                    raw_same_bin_count += 1      # unweighted (new)

    # ── Compute rate-weighted A_reco (NaN where no events landed) ─────────────
    with np.errstate(invalid="ignore", divide="ignore"):
        A_reco = np.where(
            A_reco_den > 0.0,
            A_reco_num / A_reco_den,
            0.0,
        )

    # ── Diagnostics ───────────────────────────────────────────────────────────
    survival_frac    = n_reco_total / n_true_total if n_true_total > 0.0 else 0.0
    same_bin_frac    = same_bin_sum / n_reco_total if n_reco_total > 0.0 else 0.0
    same_bin_frac_raw = (
        raw_same_bin_count / raw_valid_count if raw_valid_count > 0 else 0.0
    )

    diag = {
        "n_true_total": n_true_total,
        "n_reco_total": n_reco_total,
        "survival_frac": survival_frac,
        "same_bin_frac": same_bin_frac,         # efficiency-weighted
        "same_bin_frac_raw": same_bin_frac_raw, # unweighted, pure migration
        "n_rejected": n_rejected,
    }

    return n_reco, A_reco, diag


# ══════════════════════════════════════════════════════════════════════════════
# Case-specific helper functions (identical to 20260728, using reco arrays)
# ══════════════════════════════════════════════════════════════════════════════

def find_peak_bin_reco(n_reco, reco_mask):
    """Return (ix, iq2) of the reco-analysis-mask bin with the largest N_reco.

    Parameters
    ----------
    n_reco    : ndarray, shape (nx, nq2)  reconstructed event counts
    reco_mask : bool ndarray              reco-analysis mask from
                                          `reco_analysis_mask(proj, n_reco)`.
                                          Must not be the raw true-level
                                          `proj.accepted`; use the function
                                          above to obtain the correct mask.

    Returns
    -------
    (int, int) : (ix, iq2) indices

    Raises
    ------
    RuntimeError
        If no cell in `reco_mask` is True (empty mask). The caller (main)
        guards against this earlier via the empty-mask hard-fail block.
    """
    n = n_reco.copy()
    n[~reco_mask] = -1.0
    idx_flat = int(np.argmax(n))
    ix, iq2 = np.unravel_index(idx_flat, n.shape)
    return int(ix), int(iq2)


def find_peak_bin(proj):
    """Return (ix, iq2) indices of the accepted bin with the largest N_bin.

    Parameters
    ----------
    proj : BinnedProjection  (output of project_rates at target luminosity)

    Returns
    -------
    (int, int) : (ix, iq2) indices into proj.n_events
    """
    n = proj.n_events.copy()
    n[~proj.accepted] = -1.0   # mask out non-accepted bins
    idx_flat = int(np.argmax(n))
    ix, iq2 = np.unravel_index(idx_flat, n.shape)
    return int(ix), int(iq2)


def compute_A_cos2phi_at_bin(proj, nf2, cfg, ix, iq2, A_bag_config, base):
    """Return A_cos2φ at a single (ix, iq2) bin.

    Parameters
    ----------
    proj           : BinnedProjection
    nf2            : NuclearF2FromGrid
    cfg            : BeamConfig
    ix, iq2        : int   bin indices
    A_bag_config   : float  |A_bag| for this config (used as scale in Δ)
    base           : PartonF2  CT18NLO backend (for α_s)

    Returns
    -------
    float : A_cos2φ at the specified bin
    """
    x_val  = float(proj.x[ix, iq2])
    q2_val = float(proj.q2[ix, iq2])
    y_val  = float(proj.extras["y"][ix, iq2])

    f1_val = float(nf2.f1a(x_val, q2_val)) / cfg.ion.A
    f2_val = float(nf2.f2a(x_val, q2_val)) / cfg.ion.A

    # Use A_bag as the scale (signed convention: negative A_bag → negative Δ;
    # we use the abs value per the plan since A_bag is pre-negated in the
    # hardcoded constants)
    scale = A_bag_config   # already absolute value from A_BAG dict
    delta_val = float(delta_shape_with_alphas(
        np.array([x_val]), np.array([q2_val]), np.array([f1_val]),
        scale=scale, variant=VARIANT, base=base
    )[0])
    return float(a_cos2phi(delta_val, f1_val, f2_val, x_val, y_val))


def find_peak_q2_slice_reco(n_reco, reco_mask):
    """Return iq2 index of Q² bin with the largest sum(N_reco) across x.

    Only bins within `reco_mask` contribute to the per-Q² sum.

    Parameters
    ----------
    n_reco    : ndarray, shape (nx, nq2)  reconstructed event counts
    reco_mask : bool ndarray              reco-analysis mask from
                                          `reco_analysis_mask(proj, n_reco)`.
                                          Must not be the raw true-level
                                          `proj.accepted`.

    Returns
    -------
    int : iq2 index
    """
    n_masked = np.where(reco_mask, n_reco, 0.0)
    q2_sums = n_masked.sum(axis=0)   # shape (nq2,)
    return int(np.argmax(q2_sums))


def find_peak_q2_slice(proj):
    """Return iq2 index of the Q² bin with the largest sum(N_bin) across x.

    Only accepted bins contribute to the sum.

    Parameters
    ----------
    proj : BinnedProjection

    Returns
    -------
    int : iq2 index
    """
    n_masked = np.where(proj.accepted, proj.n_events, 0.0)
    q2_sums = n_masked.sum(axis=0)   # shape (nq2,)
    return int(np.argmax(q2_sums))


def compute_A_cos2phi_q2slice(proj, nf2, cfg, iq2, A_bag_config, base):
    """Return rate-weighted ⟨A_cos2φ⟩ over accepted x bins at a given Q² slice.

    ⟨A_cos2φ⟩ = Σ(A_cos2φ,bin · N_bin) / Σ(N_bin)   over accepted x at iq2

    Parameters
    ----------
    proj           : BinnedProjection
    nf2            : NuclearF2FromGrid
    cfg            : BeamConfig
    iq2            : int   Q² slice index
    A_bag_config   : float  |A_bag| for this config
    base           : PartonF2

    Returns
    -------
    float : rate-weighted ⟨A_cos2φ⟩ over x at the given Q²
    float : total N in the slice (sum over accepted x)
    """
    # Select accepted x indices in this Q² slice
    accepted_slice = proj.accepted[:, iq2]
    ix_list = np.where(accepted_slice)[0]

    if len(ix_list) == 0:
        return 0.0, 0.0

    x_arr  = proj.x[ix_list, iq2]
    q2_arr = proj.q2[ix_list, iq2]
    y_arr  = proj.extras["y"][ix_list, iq2]
    n_arr  = proj.n_events[ix_list, iq2]

    f1_arr = nf2.f1a(x_arr, q2_arr) / cfg.ion.A
    f2_arr = nf2.f2a(x_arr, q2_arr) / cfg.ion.A

    delta_arr = delta_shape_with_alphas(
        x_arr, q2_arr, f1_arr,
        scale=A_bag_config, variant=VARIANT, base=base
    )
    amp_arr = np.array([
        float(a_cos2phi(delta_arr[i], f1_arr[i], f2_arr[i], x_arr[i], y_arr[i]))
        for i in range(len(x_arr))
    ])

    total_n = float(n_arr.sum())
    if total_n <= 0.0:
        return 0.0, 0.0

    a_weighted = float((amp_arr * n_arr).sum() / total_n)
    return a_weighted, total_n


def compute_A_cos2phi_q2slice_reco(n_reco, A_reco, reco_mask, iq2):
    """Rate-weighted ⟨A_cos2φ⟩ over x in a Q² slice using reco arrays.

    Only bins within `reco_mask[:, iq2]` are included.  This replaces the
    former `n_slice > 0` guard, which applied no analytic acceptance cuts.

    Parameters
    ----------
    n_reco    : ndarray (nx, nq2) reconstructed event counts
    A_reco    : ndarray (nx, nq2) rate-weighted A_cos2φ per reco bin
    reco_mask : bool ndarray      reco-analysis mask from
                                  `reco_analysis_mask(proj, n_reco)`.
                                  Must not be the raw true-level
                                  `proj.accepted`.
    iq2       : int               Q² slice index

    Returns
    -------
    float : rate-weighted ⟨A_cos2φ⟩
    float : total N_reco in this slice
    """
    n_slice = n_reco[:, iq2]
    A_slice = A_reco[:, iq2]
    valid = reco_mask[:, iq2]
    if not np.any(valid):
        return 0.0, 0.0
    total = float(n_slice[valid].sum())
    a_w   = float((A_slice[valid] * n_slice[valid]).sum() / total)
    return a_w, total


def compute_A_cos2phi_integrated(proj, nf2, cfg, A_bag_config, base):
    """Return rate-weighted ⟨A_cos2φ⟩ over ALL accepted (x, Q²) bins.

    ⟨A_cos2φ⟩ = Σ(A_cos2φ,bin · N_bin) / Σ(N_bin)   over all accepted bins

    Parameters
    ----------
    proj           : BinnedProjection
    nf2            : NuclearF2FromGrid
    cfg            : BeamConfig
    A_bag_config   : float  |A_bag| for this config
    base           : PartonF2

    Returns
    -------
    float : rate-weighted ⟨A_cos2φ⟩ over all accepted bins
    float : total N in all accepted bins
    """
    ix_flat, iq2_flat = np.where(proj.accepted)

    if len(ix_flat) == 0:
        return 0.0, 0.0

    x_arr  = proj.x[ix_flat, iq2_flat]
    q2_arr = proj.q2[ix_flat, iq2_flat]
    y_arr  = proj.extras["y"][ix_flat, iq2_flat]
    n_arr  = proj.n_events[ix_flat, iq2_flat]

    f1_arr = nf2.f1a(x_arr, q2_arr) / cfg.ion.A
    f2_arr = nf2.f2a(x_arr, q2_arr) / cfg.ion.A

    delta_arr = delta_shape_with_alphas(
        x_arr, q2_arr, f1_arr,
        scale=A_bag_config, variant=VARIANT, base=base
    )
    amp_arr = np.array([
        float(a_cos2phi(delta_arr[i], f1_arr[i], f2_arr[i], x_arr[i], y_arr[i]))
        for i in range(len(x_arr))
    ])

    total_n = float(n_arr.sum())
    if total_n <= 0.0:
        return 0.0, 0.0

    a_weighted = float((amp_arr * n_arr).sum() / total_n)
    return a_weighted, total_n


def compute_A_cos2phi_integrated_reco(n_reco, A_reco, reco_mask):
    """Rate-weighted ⟨A_cos2φ⟩ over all reco-analysis-mask bins.

    Only bins within `reco_mask` are included.  This replaces the former
    `n_reco > 0` guard, which applied no analytic acceptance cuts and thus
    could include bins outside the detector's kinematic acceptance.

    Parameters
    ----------
    n_reco    : ndarray (nx, nq2) reconstructed event counts
    A_reco    : ndarray (nx, nq2) rate-weighted A_cos2φ per reco bin
    reco_mask : bool ndarray      reco-analysis mask from
                                  `reco_analysis_mask(proj, n_reco)`.
                                  Must not be the raw true-level
                                  `proj.accepted`.

    Returns
    -------
    float : rate-weighted ⟨A_cos2φ⟩
    float : total N_reco
    """
    valid = reco_mask
    if not np.any(valid):
        return 0.0, 0.0
    total = float(n_reco[valid].sum())
    a_w   = float((A_reco[valid] * n_reco[valid]).sum() / total)
    return a_w, total


# ══════════════════════════════════════════════════════════════════════════════
# φ-distribution plot builder
# ══════════════════════════════════════════════════════════════════════════════

def build_phi_plot(A_ref, N_flat, A_bag_config, config_label, case_label,
                   case_tag, config_tag, outdir,
                   extra_info=None,
                   n_phi_bins=None, bin_width_deg=5,
                   filename_suffix=None):
    """Build one φ-modulation plot for a single (config, case) combination.

    Parameters
    ----------
    A_ref          : float   reference A_cos2φ for this case/config
    N_flat         : float   events per φ bin at flat (unmodulated) yield
                             (already divided by n_phi_bins by the caller)
    A_bag_config   : float   |A_bag| for this config (positive)
    config_label   : str     human-readable config label for plot title
    case_label     : str     human-readable case label for plot title
    case_tag       : str     one of 'peakbin', 'q2slice', 'integrated'
    config_tag     : str     one of 'low', 'mid', 'top'
    outdir         : pathlib.Path
    extra_info     : str or None   optional extra annotation text (appended)
    n_phi_bins     : int or None   number of φ bins; defaults to N_PHI (72 → 5°)
    bin_width_deg  : float         bin width in degrees (used in title/annotation)
    filename_suffix: str or None   if given, appended before '.png' instead of
                                   the default '_{config_tag}' suffix, e.g.
                                   '_top_5deg' → filename ends with that string

    Returns
    -------
    pathlib.Path : output PNG path
    """
    if n_phi_bins is None:
        n_phi_bins = N_PHI   # default: 72 bins of 5°

    # φ grid: n_phi_bins bins; centers at bin_width/2, 3·bin_width/2, …
    phi = np.linspace(0, 2 * np.pi, n_phi_bins, endpoint=False) + np.pi / n_phi_bins

    # Three signal scales
    scales = [0.1 * A_bag_config, A_bag_config, 10.0 * A_bag_config]
    scale_labels = [
        r"$s = 0.1\,A_\mathrm{bag}$",
        r"$s = A_\mathrm{bag}$ (bag reference)",
        r"$s = 10\,A_\mathrm{bag}$",
    ]
    scale_colors = ["lightsteelblue", "black", "red"]
    scale_lws    = [1.2, 1.8, 1.5]

    # Compute modulation curves and check positivity
    modulations = []
    for i_s, s in enumerate(scales):
        # Modulation: y_model(φ, s) = P_zz · (s/A_bag) · A_ref · cos(2φ)
        amplitude = PZZ * (s / A_bag_config) * A_ref
        y_model = amplitude * np.cos(2.0 * phi)

        # Positivity check
        max_abs_y = float(np.max(np.abs(y_model)))
        if max_abs_y >= 1.0:
            print(
                f"  WARNING: |y_model| = {max_abs_y:.4f} >= 1 for "
                f"config={config_tag}, case={case_tag}, "
                f"s={s:.4f} (={['0.1','1','10'][i_s]}×A_bag={A_bag_config:.3f})"
            )

        # N_phi and statistical error
        n_phi  = N_flat * (1.0 + y_model)
        # sigma_phi: Poisson error propagated to fractional units
        sigma_phi = np.sqrt(np.maximum(n_phi, 0.0)) / N_flat

        modulations.append({
            "s": s,
            "label": scale_labels[i_s],
            "color": scale_colors[i_s],
            "lw": scale_lws[i_s],
            "y_model": y_model,
            "n_phi": n_phi,
            "sigma_phi": sigma_phi,
        })

    # ── Build plot ────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5.5))

    # Reference flat line at y=0
    ax.axhline(0.0, color="gray", ls="--", lw=0.8, alpha=0.7, zorder=1)

    # Three signal-scale curves
    for m in modulations:
        ax.plot(phi, m["y_model"], color=m["color"], lw=m["lw"],
                label=m["label"], zorder=3)

    # Error bars at the middle scale (A_bag, index 1)
    m_mid = modulations[1]
    ax.errorbar(
        phi, m_mid["y_model"], yerr=m_mid["sigma_phi"],
        fmt="o", color="black", markersize=3, capsize=1,
        label=r"error bars from Poisson $\sqrt{N_\phi}$ (shown at $s=A_\mathrm{bag}$)",
        zorder=4, alpha=0.85
    )

    # x-axis ticks at 0, π/2, π, 3π/2, 2π
    ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
    ax.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
    ax.set_xlim(0, 2 * np.pi)

    ax.set_xlabel(r"$\phi$ [rad]", fontsize=11)
    ax.set_ylabel(r"fractional yield modulation $y(\phi)$", fontsize=10)

    # Include bin-width in title only when it differs from the 5° default
    bin_info = f" ({bin_width_deg}° bins)" if bin_width_deg != 5 else ""
    title_line1 = (
        rf"$^6$Li fractional yield modulation vs $\phi$ — {config_label} config,"
        rf" {case_label}{bin_info} [detector-realistic]"
    )
    title_line2 = (
        r"$L = 10\,\mathrm{fb}^{-1}/\mathrm{u}$, mid$\_x$ shape, EPPS21, "
        r"R1998, Cloet $P_{zz}$=0.267; ePIC tracking + e-ID eff."
    )
    ax.set_title(title_line1 + "\n" + title_line2, fontsize=9)

    # Legend
    ax.legend(fontsize=7.5, loc="upper right", framealpha=0.8)

    # Text annotation in upper-left corner
    ann_lines = [
        rf"$\langle A_{{\cos 2\phi}}\rangle$ = {A_ref:.4f}",
        rf"$N_\mathrm{{flat}}$ = {N_flat:.2e} events/bin ({bin_width_deg}°)",
        rf"$|A_\mathrm{{bag}}|$ = {A_bag_config:.3f}",
        r"[reco, MC-smeared]",
    ]
    if extra_info is not None:
        ann_lines.append(extra_info)
    annotation = "\n".join(ann_lines)
    ax.text(
        0.02, 0.98, annotation,
        transform=ax.transAxes,
        fontsize=8, va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7, ec="gray"),
    )

    fig.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)
    if filename_suffix is not None:
        filename = f"money_delta_20260729_phimodulation_{case_tag}{filename_suffix}.png"
    else:
        filename = f"money_delta_20260729_phimodulation_{case_tag}_{config_tag}.png"
    outpath = outdir / filename
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return outpath


# ══════════════════════════════════════════════════════════════════════════════
# Per-bin (x, Q²) heatmap — motivation plots (reco arrays)
# ══════════════════════════════════════════════════════════════════════════════

def build_perbin_heatmap(cfg, base, config_tag, config_label, A_bag_config,
                         luminosity, outdir, proj, n_reco, A_reco,
                         reco_mask):
    """Build per-bin 1×3 subplot heatmap using detector-smeared event counts.

    Uses reco n_reco and A_reco in place of proj.n_events and per-bin A_cos2φ.
    Interpretation A: Δ = s · α_s(Q²) · F₁ · x^α · (1-x)^β.

    Parameters
    ----------
    cfg            : BeamConfig
    base           : PartonF2   CT18NLO backend (for α_s)
    config_tag     : str        'low', 'mid', or 'top'
    config_label   : str        human-readable beam description
    A_bag_config   : float      |A_bag| for this config (positive)
    luminosity     : float      integrated luminosity [fb⁻¹/nucleon]
    outdir         : pathlib.Path
    proj           : BinnedProjection   (for grid edges only; acceptance is
                                         encoded in `reco_mask`)
    n_reco         : ndarray (nx, nq2)  reconstructed event counts
    A_reco         : ndarray (nx, nq2)  rate-weighted A_cos2φ per reco bin
    reco_mask      : bool ndarray       reco-analysis mask from
                                        `reco_analysis_mask(proj, n_reco)`.
                                        The inline `accepted & (n >= MIN_EVENTS)`
                                        expression previously used here is
                                        replaced by this precomputed mask to
                                        ensure consistency with all other
                                        reco-side helpers.

    Returns
    -------
    pathlib.Path : output PNG path
    """
    L = luminosity

    n_events = n_reco     # use reconstructed counts
    amp_bin  = A_reco     # use reconstructed A_cos2φ

    # ── Validity mask ─────────────────────────────────────────────────────────
    # reco_mask already encodes proj.accepted AND n_reco >= MIN_EVENTS.
    use = reco_mask
    mask_valid = use & (np.abs(amp_bin) > 1e-10)

    # ── Build per-bin arrays (NaN outside valid mask) ─────────────────────────
    amp_abs = np.full_like(n_events, np.nan, dtype=float)   # |A_bin|
    delta_a = np.full_like(n_events, np.nan, dtype=float)   # δA_bin
    rel_unc = np.full_like(n_events, np.nan, dtype=float)   # |δA_bin / A_bin|

    amp_abs[mask_valid] = np.abs(amp_bin[mask_valid])
    da_vals = np.sqrt(2.0 / np.maximum(n_events[mask_valid], 1e-12)) / PZZ
    delta_a[mask_valid] = da_vals
    rel_unc[mask_valid] = da_vals / np.abs(amp_bin[mask_valid])

    # ── Bin edges for pcolormesh (log scale) ─────────────────────────────────
    x_edges  = proj.x_edges    # shape (nx+1,)
    q2_edges = proj.q2_edges   # shape (nq2+1,)

    # pcolormesh arrays: shape (nq2, nx) — transpose of (nx, nq2)
    X_edge_2d, Q2_edge_2d = np.meshgrid(x_edges, q2_edges)   # (nq2+1, nx+1)
    C_amp   = amp_abs.T   # shape (nq2, nx)
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
    cbar.set_label(r"$|A_{\cos 2\phi, \mathrm{bin}}|$ [reco]", fontsize=10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$", fontsize=10)
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]", fontsize=10)
    ax.set_title("|A_bin| (signal, reco)", fontsize=10)

    # ── Middle: δA_bin ────────────────────────────────────────────────────────
    ax = axes[1]
    pcm = ax.pcolormesh(X_edge_2d, Q2_edge_2d, C_delta,
                        norm=_lognorm(C_delta), cmap="plasma", shading="flat")
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label(r"$\delta A_{\cos 2\phi, \mathrm{bin}}$ [reco]", fontsize=10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$", fontsize=10)
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]", fontsize=10)
    ax.set_title("δA_bin (uncertainty, reco)", fontsize=10)

    # ── Right: |δA_bin / A_bin| ───────────────────────────────────────────────
    ax = axes[2]
    pcm = ax.pcolormesh(X_edge_2d, Q2_edge_2d, C_rel,
                        norm=_lognorm(C_rel), cmap="magma", shading="flat")
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label(
        r"$|\delta A_{\cos 2\phi, \mathrm{bin}} / A_{\cos 2\phi, \mathrm{bin}}|$ [reco]",
        fontsize=10,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$", fontsize=10)
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]", fontsize=10)
    ax.set_title("|δA_bin / A_bin| (relative, reco)", fontsize=10)

    fig.suptitle(
        f"$^6$Li per-bin decomposition [detector-realistic] — {config_label} config, "
        f"L = {L:g} fb$^{{-1}}$/u\n"
        f"Interpretation A: "
        r"$\Delta = s \cdot \alpha_s(Q^2) \cdot F_1 \cdot x^\alpha (1-x)^\beta$"
        f", $s = |A_\\mathrm{{bag}}|$ = {A_bag_config:.3f}, "
        "mid_x shape, EPPS21, R1998, Cloet $P_{{zz}}$=0.267; ePIC tracking + e-ID",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"money_delta_20260729_perbin_{config_tag}.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath.name}")

    # ── Stdout diagnostics ────────────────────────────────────────────────────
    print(f"── perbin diagnostics (reco): {config_label}, L={L:g} fb⁻¹/u ──")
    print(f"  |A_bin|:     min={np.nanmin(amp_abs):.2e}  "
          f"median={np.nanmedian(amp_abs):.2e}  max={np.nanmax(amp_abs):.2e}")
    print(f"  δA_bin:      min={np.nanmin(delta_a):.2e}  "
          f"median={np.nanmedian(delta_a):.2e}  max={np.nanmax(delta_a):.2e}")
    print(f"  |δA/A|_bin:  min={np.nanmin(rel_unc):.2e}  "
          f"median={np.nanmedian(rel_unc):.2e}  max={np.nanmax(rel_unc):.2e}")

    # Location of max |δA/A|
    flat_idx = np.nanargmax(rel_unc)
    ix_r, iq_r = np.unravel_index(flat_idx, rel_unc.shape)
    x_mid_r  = (0.5 * (x_edges[ix_r]  + x_edges[ix_r  + 1])
                if ix_r  < len(x_edges)  - 1 else x_edges[ix_r])
    q2_mid_r = (0.5 * (q2_edges[iq_r] + q2_edges[iq_r + 1])
                if iq_r < len(q2_edges) - 1 else q2_edges[iq_r])
    print(f"  Location of max |δA/A|:  x={x_mid_r:.2e},  Q²={q2_mid_r:.2f} GeV²")

    # Peak-rate bin (for cross-ref with Case 1)
    n_bin_map = np.where(mask_valid, n_events, np.nan)
    if not np.all(np.isnan(n_bin_map)):
        flat_peak = np.nanargmax(n_bin_map)
        ix_p, iq_p = np.unravel_index(flat_peak, n_bin_map.shape)
        x_pk  = float(proj.x[ix_p, iq_p])
        q2_pk = float(proj.q2[ix_p, iq_p])
        a_pk  = float(amp_abs[ix_p, iq_p]) if not np.isnan(amp_abs[ix_p, iq_p]) else float("nan")
        print(f"  Peak-rate bin (reco, cross-ref Case 1): "
              f"x={x_pk:.2e}, Q²={q2_pk:.2f}, |A_bin|={a_pk:.2e}")

    return outpath


# ══════════════════════════════════════════════════════════════════════════════
# Smearing diagnostic printer
# ══════════════════════════════════════════════════════════════════════════════

def print_smearing_diagnostics(config_label, diag,
                                proj, n_reco, A_reco,
                                ix_peak_true, iq2_peak_true,
                                A_true_peakbin,
                                reco_mask, true_accepted_count: int):
    """Print per-config smearing effects summary to stdout.

    Parameters
    ----------
    config_label       : str    human-readable config label
    diag               : dict   from smear_config()
    proj               : BinnedProjection  (grid geometry only; this helper
                                            must NOT access proj.accepted —
                                            see S1 in _check_reco_mask_invariants.py)
    n_reco             : ndarray (nx, nq2)
    A_reco             : ndarray (nx, nq2)
    ix_peak_true       : int    true-level peak-rate bin ix
    iq2_peak_true      : int    true-level peak-rate bin iq2
    A_true_peakbin     : float  true (unsmeared) A_cos2φ at peak-rate bin
    reco_mask          : bool ndarray   reco-analysis mask from
                                        `reco_analysis_mask(proj, n_reco)`.
                                        Forwarded to `find_peak_bin_reco`.
    true_accepted_count : int   number of true-accepted bins (precomputed by
                                caller as int(proj.accepted.sum())); passed
                                explicitly so this helper never reads
                                proj.accepted directly (S1 invariant).
    """
    n_true_total  = diag["n_true_total"]
    n_reco_total  = diag["n_reco_total"]
    survival_frac = diag["survival_frac"]
    same_bin_frac = diag["same_bin_frac"]

    # Reco peak-rate bin — use reco_mask (not proj.accepted) as required by S1
    ix_r, iq2_r = find_peak_bin_reco(n_reco, reco_mask)
    x_reco_pk   = float(proj.x[ix_r, iq2_r])
    q2_reco_pk  = float(proj.q2[ix_r, iq2_r])
    n_reco_pk   = float(n_reco[ix_r, iq2_r])
    A_reco_pk   = float(A_reco[ix_r, iq2_r])

    # True peak-rate bin
    x_true_pk   = float(proj.x[ix_peak_true, iq2_peak_true])
    q2_true_pk  = float(proj.q2[ix_peak_true, iq2_peak_true])
    n_true_pk   = float(proj.n_events[ix_peak_true, iq2_peak_true])

    n_reco_acc_bins = int(reco_mask.sum())

    print(f"── smearing diagnostics: {config_label} ──")
    print(f"  N_events before smearing:          {n_true_total:.3e}")
    print(f"  N_events after smearing + eff:     {n_reco_total:.3e}"
          f"  ({100.0 * survival_frac:.1f}% survival)")
    print(f"  Reco-accepted bins: {n_reco_acc_bins} (vs true-accepted: {true_accepted_count})")
    print(f"  Bin migration (raw, efficiency-independent):   "
          f"{100.0 * diag['same_bin_frac_raw']:.0f}% of MC events land in original bin")
    print(f"  Bin migration (efficiency-weighted yield):     "
          f"{100.0 * same_bin_frac:.0f}% of reco yield originates from original bin")
    print(f"  MC events rejected (physicality):  {diag['n_rejected']:,}")
    print(f"  Peak-rate bin (reco):  x={x_reco_pk:.2e},"
          f" Q²={q2_reco_pk:.2f} GeV², N_reco={n_reco_pk:.3e}")
    print(f"  Peak-rate bin (true):  x={x_true_pk:.2e},"
          f" Q²={q2_true_pk:.2f} GeV², N_true={n_true_pk:.3e}"
          f"  ← unchanged from money_delta_20260728")
    print(f"  Diagnostic: reco A_cos2φ at peak-rate bin (reco) = {A_reco_pk:.4e},"
          f" vs no-smearing A = {A_true_peakbin:.4e}")
    print()


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description=(
            "6Li φ-distribution plots — detector-realistic version with MC smearing.\n\n"
            "Detector-realistic extension of money_delta_20260728.py.\n"
            "Produces 16 individual PNG plots showing the observable φ-modulation\n"
            "that an experimenter would see in the cos(2φ) tensor asymmetry\n"
            "measurement, at three signal scales (0.1·A_bag, A_bag, 10·A_bag),\n"
            "for three EIC beam configurations and three levels of kinematic\n"
            "integration, with ePIC detector smearing applied.\n\n"
            "New in 20260729 vs 20260728:\n"
            "  1. Monte Carlo smearing of scattered electron (E', θ, φ) using\n"
            "     ePIC tracking resolutions (piecewise in η).\n"
            "  2. Bin migration in (x, Q²) from reconstructed kinematics.\n"
            "  3. η-dependent electron-ID efficiency (ATHENA+ECCE anchor points).\n\n"
            "Cases:\n"
            "  peakbin:    Peak (x, Q²) bin (reconstructed)\n"
            "  q2slice:    Peak Q² slice (rate-weighted ⟨A_cos2φ⟩ over x, reco)\n"
            "  integrated: Fully integrated (rate-weighted ⟨A_cos2φ⟩, reco)\n\n"
            "Configs:\n"
            "  low: 5 GeV e × 27.5 GeV/u ⁶Li\n"
            "  mid: 10 GeV e × 50 GeV/u ⁶Li\n"
            "  top: 18 GeV e × 137.5 GeV/u ⁶Li\n\n"
            "Physics: Δ(x,Q²) = s · α_s(Q²) · F₁(x,Q²) · x^α (1-x)^β (mid_x).\n"
            "A_bag hardcoded from money_delta_20260724.py mid_x output:\n"
            "  LOW=0.318, MID=0.310, TOP=0.297\n\n"
            "ECal note: ECal resolution not implemented; tracking-only per user\n"
            "decision. See module docstring for details."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--outdir",
        default="fastsim/out/money_delta",
        help="Output directory for PNGs (default: fastsim/out/money_delta)",
    )
    ap.add_argument(
        "--n-mc",
        type=int,
        default=1000,
        help="MC pseudo-events per true bin (default: 1000)",
    )
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    n_mc = args.n_mc

    # ── Startup banner ────────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("money_delta_20260729.py — 6Li φ-modulation projections [DETECTOR-REALISTIC]")
    print("=" * 72)
    print()
    print("Observable: dN/dφ = N_flat · [1 + P_zz · (s/A_bag) · <A_cos2φ> · cos(2φ)]")
    print(f"Luminosity: {LUMI_FB:g} fb⁻¹/nucleon  (fixed)")
    print(f"Shape:      mid_x (α=0.7, β=3)")
    print(f"R:          R1998")
    print(f"PDF:        {EPPS21_SET}")
    print(f"P_zz:       {PZZ} (Cloet)")
    print(f"N_φ bins:   {N_PHI} (5° each)")
    print()
    print("Detector model (new in 20260729):")
    print("  Tracking momentum resolution: σ_p/p = √(a²·p² + b²) piecewise in η")
    print("  Tracking angular resolution: σ_θ = σ_φ (mrad piecewise in η)")
    print("  Electron-ID efficiency: linear interp from ATHENA+ECCE anchor points")
    print("  ECal: NOT implemented (tracking-only per user decision)")
    print(f"  MC events per true bin: {n_mc}")
    print()
    print("Hardcoded |A_bag| values:")
    for k, v in A_BAG.items():
        print(f"  {k.upper()}: {v}")
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
    print()

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

    # ── Three EIC beam configurations for ⁶Li ────────────────────────────────
    all_configs = [
        (
            BeamConfig(electron_energy=5.0, ion=LI6, ion_momentum_per_nucleon=27.5),
            "low",
            r"LOW (5$\times$27.5 GeV/u)",
        ),
        (
            BeamConfig(electron_energy=10.0, ion=LI6, ion_momentum_per_nucleon=50.0),
            "mid",
            r"MID (10$\times$50 GeV/u)",
        ),
        (
            BeamConfig(electron_energy=18.0, ion=LI6, ion_momentum_per_nucleon=137.5),
            "top",
            r"TOP (18$\times$137.5 GeV/u)",
        ),
    ]

    # Cases: (tag, label)
    cases = [
        ("peakbin",    "peak bin"),
        ("q2slice",    "peak Q² slice"),
        ("integrated", "fully integrated"),
    ]

    # ── Summary data for stdout table ─────────────────────────────────────────
    summary_rows = []

    # ── Main loop: 3 configs ──────────────────────────────────────────────────
    sc = fom.Scenario(
        lumi_fb_per_nucleon=LUMI_FB,
        pol_ion_tensor=PZZ,
        q2_min=2.0,
    )

    with r_override(r1998):
        for cfg, config_tag, config_label in all_configs:
            print(f"{'='*72}")
            print(f"Config: {config_tag.upper()} — {config_label}")
            print(f"{'='*72}")

            a_bag_config = A_BAG[config_tag]
            nf2_obj = NuclearF2FromGrid(cfg.ion, EPPS21_SET)

            print(f"  Running project_rates at L = {LUMI_FB:g} fb⁻¹/nucleon …",
                  end=" ", flush=True)
            proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_obj)
            print("done.")

            # Accepted events grid info (true-level, for reference)
            n_total_acc_true = float(proj.n_events[proj.accepted].sum())
            n_acc_bins  = int(proj.accepted.sum())
            print(f"  Accepted bins (true): {n_acc_bins}  |  "
                  f"Total events (true): {n_total_acc_true:.3e}")
            print()

            # ── Monte Carlo smearing ──────────────────────────────────────────
            print(f"  Running MC smearing ({n_mc} events/bin × {n_acc_bins} bins) …",
                  end=" ", flush=True)
            rng = np.random.default_rng(seed=42)
            n_reco, A_reco, diag = smear_config(
                cfg=cfg,
                proj=proj,
                base=base,
                nf2=nf2_obj,
                A_bag_config=a_bag_config,
                rng=rng,
                n_mc=n_mc,
            )
            print("done.")
            print()

            # ── Reconstructed-analysis mask: the SINGLE source of truth for
            #    reco selection.
            #    = proj.accepted (full DIS + eta + E' acceptance) AND
            #      n_reco >= MIN_EVENTS.
            #    See reco_analysis_mask docstring for the acceptance decision.
            reco_mask = reco_analysis_mask(proj, n_reco, min_events=MIN_EVENTS)
            n_reco_acc_bins = int(reco_mask.sum())
            print(f"  Reco-accepted bins: {n_reco_acc_bins}  "
                  f"(vs true-accepted: {n_acc_bins})")

            # ── Cross-consistency assertions ──────────────────────────────────
            assert reco_mask.shape == n_reco.shape == A_reco.shape, (
                "shape mismatch: reco_mask / n_reco / A_reco must be identical")
            assert not np.any(reco_mask & (n_reco < MIN_EVENTS)), (
                "reco_mask leaked a bin below MIN_EVENTS floor")
            assert bool(np.all(reco_mask <= proj.accepted)), (
                "reco_mask must be a subset of proj.accepted (acceptance invariant)")

            # ── Fail loudly if the mask is empty — see plan §"Failure mode".
            #    NOTE: the f-string is the DIRECT argument of RuntimeError(...)
            #    (not assigned to a `msg` variable first) so that the AST check
            #    S5 in _check_reco_mask_invariants.py can verify its presence.
            if n_reco_acc_bins == 0:
                print(
                    f"ERROR: Empty reco-analysis mask for config {config_tag!r}: "
                    f"true-accepted bins = {n_acc_bins}, "
                    f"n_reco_total = {diag['n_reco_total']:.3e}.",
                    file=sys.stderr,
                )
                raise RuntimeError(
                    f"Empty reco-analysis mask for config {config_tag!r}: "
                    f"true-accepted bins = {n_acc_bins}, "
                    f"n_reco_total = {diag['n_reco_total']:.3e}. "
                    f"Cannot produce the expected PNG set. "
                    f"Check --n-mc, Scenario cuts, or smearing kernel."
                )
            print()

            # ── Smearing diagnostics (need true peak bin A for comparison) ────
            ix_peak_true, iq2_peak_true = find_peak_bin(proj)
            A_true_peakbin = compute_A_cos2phi_at_bin(
                proj, nf2_obj, cfg, ix_peak_true, iq2_peak_true,
                a_bag_config, base
            )
            print_smearing_diagnostics(
                config_label=config_label,
                diag=diag,
                proj=proj,
                n_reco=n_reco,
                A_reco=A_reco,
                ix_peak_true=ix_peak_true,
                iq2_peak_true=iq2_peak_true,
                A_true_peakbin=A_true_peakbin,
                reco_mask=reco_mask,
                true_accepted_count=n_acc_bins,
            )

            # ── Heatmap: per-bin motivation plot (uses reco arrays) ───────────
            build_perbin_heatmap(
                cfg=cfg,
                base=base,
                config_tag=config_tag,
                config_label=config_label,
                A_bag_config=a_bag_config,
                luminosity=LUMI_FB,
                outdir=outdir,
                proj=proj,
                n_reco=n_reco,
                A_reco=A_reco,
                reco_mask=reco_mask,
            )
            print()

            # ── Case 1: Peak bin (RECO) ───────────────────────────────────────
            ix_peak_r, iq2_peak_r = find_peak_bin_reco(n_reco, reco_mask)
            n_peak_reco  = float(n_reco[ix_peak_r, iq2_peak_r])
            x_peak_reco  = float(proj.x[ix_peak_r, iq2_peak_r])
            q2_peak_reco = float(proj.q2[ix_peak_r, iq2_peak_r])
            A_peakbin    = float(A_reco[ix_peak_r, iq2_peak_r])

            N_flat_peakbin = n_peak_reco / N_PHI
            mod_at_abag    = abs(PZZ * A_peakbin)
            mod_at_10abag  = abs(PZZ * 10.0 * A_peakbin)

            print(f"  Case 1 (peak bin, reco):")
            print(f"    x = {x_peak_reco:.2e},  Q² = {q2_peak_reco:.2f} GeV²  "
                  f"(ix={ix_peak_r}, iq2={iq2_peak_r})")
            print(f"    N_reco = {n_peak_reco:.2e} events")
            print(f"    A_cos2φ (reco) = {A_peakbin:.4e}")
            print(f"    P_zz·A_cos2φ (at A_bag) = {mod_at_abag:.4f}")
            print(f"    Max modulation at 10×A_bag = {mod_at_10abag:.4f}")

            pos_ok_1 = mod_at_10abag < 1.0
            print(f"    Positivity check: {'OK' if pos_ok_1 else 'WARN'}")

            extra_peakbin = (f"peak bin (reco): x={x_peak_reco:.2e}, "
                             f"Q²={q2_peak_reco:.2f} GeV²")

            if config_tag == "top":
                # ── TOP + Case 1 only: scan over φ-bin widths ─────────────────
                print()
                print(f"  → TOP config, Case 1: scanning φ bin widths "
                      f"{{5°, 10°, 20°, 30°, 45°}}")
                for bin_deg in [5, 10, 20, 30, 45]:
                    n_bins = int(360 / bin_deg)
                    n_flat_scan = n_peak_reco / n_bins
                    err_scan = (1.0 / np.sqrt(n_flat_scan)
                                if n_flat_scan > 0 else float("nan"))
                    print(
                        f"    {bin_deg}° bins ({n_bins:2d} total): "
                        f"N_flat = {n_flat_scan:.2e} events/bin (reco), "
                        f"error bar ≈ {err_scan:.2e}"
                    )
                print()
                phi_bin_paths = []
                for bin_deg in [5, 10, 20, 30, 45]:
                    n_bins = int(360 / bin_deg)
                    n_flat_scan = n_peak_reco / n_bins
                    suffix = f"_top_{bin_deg}deg"
                    outpath = build_phi_plot(
                        A_ref=A_peakbin,
                        N_flat=n_flat_scan,
                        A_bag_config=a_bag_config,
                        config_label=config_label,
                        case_label="peak bin",
                        case_tag="peakbin",
                        config_tag=config_tag,
                        outdir=outdir,
                        extra_info=extra_peakbin,
                        n_phi_bins=n_bins,
                        bin_width_deg=bin_deg,
                        filename_suffix=suffix,
                    )
                    phi_bin_paths.append(outpath)
                    print(f"    wrote {outpath.name}")
            else:
                # All other configs: single 5° plot
                outpath = build_phi_plot(
                    A_ref=A_peakbin,
                    N_flat=N_flat_peakbin,
                    A_bag_config=a_bag_config,
                    config_label=config_label,
                    case_label="peak bin",
                    case_tag="peakbin",
                    config_tag=config_tag,
                    outdir=outdir,
                    extra_info=extra_peakbin,
                )
                print(f"  wrote {outpath.name}")

            summary_rows.append({
                "config": config_tag.upper(),
                "case": "Case 1 (peak bin, reco)",
                "A_ref": A_peakbin,
                "N_flat": N_flat_peakbin,
                "mod_at_abag": mod_at_abag,
                "mod_at_10abag": mod_at_10abag,
                "pos_ok": pos_ok_1,
                "extra": (f"x = {x_peak_reco:.2e}, Q² = {q2_peak_reco:.2f} GeV², "
                          f"N_reco = {n_peak_reco:.2e}"),
            })

            print()

            # ── Case 2: Peak Q² slice (RECO) ──────────────────────────────────
            iq2_peak_slice_r = find_peak_q2_slice_reco(n_reco, reco_mask)
            q2_slice_val = float(proj.q2[0, iq2_peak_slice_r])

            A_q2slice, N_q2slice = compute_A_cos2phi_q2slice_reco(
                n_reco, A_reco, reco_mask, iq2_peak_slice_r
            )
            N_flat_q2slice = N_q2slice / N_PHI
            mod_at_abag_q2    = abs(PZZ * A_q2slice)
            mod_at_10abag_q2  = abs(PZZ * 10.0 * A_q2slice)

            print(f"  Case 2 (peak Q² slice, reco):")
            print(f"    Peak Q² slice: iq2={iq2_peak_slice_r}, "
                  f"Q² ≈ {q2_slice_val:.2f} GeV²")
            print(f"    Sum N_reco over x = {N_q2slice:.2e} events")
            print(f"    ⟨A_cos2φ⟩ (reco) = {A_q2slice:.4e}")
            print(f"    P_zz·⟨A_cos2φ⟩ (at A_bag) = {mod_at_abag_q2:.4f}")
            print(f"    Max modulation at 10×A_bag = {mod_at_10abag_q2:.4f}")

            pos_ok_2 = mod_at_10abag_q2 < 1.0
            print(f"    Positivity check: {'OK' if pos_ok_2 else 'WARN'}")

            extra_q2slice = f"peak Q² slice (reco): Q²≈{q2_slice_val:.2f} GeV²"
            outpath = build_phi_plot(
                A_ref=A_q2slice,
                N_flat=N_flat_q2slice,
                A_bag_config=a_bag_config,
                config_label=config_label,
                case_label=r"peak $Q^2$ slice",
                case_tag="q2slice",
                config_tag=config_tag,
                outdir=outdir,
                extra_info=extra_q2slice,
            )
            print(f"  wrote {outpath.name}")

            summary_rows.append({
                "config": config_tag.upper(),
                "case": "Case 2 (Q² slice, reco)",
                "A_ref": A_q2slice,
                "N_flat": N_flat_q2slice,
                "mod_at_abag": mod_at_abag_q2,
                "mod_at_10abag": mod_at_10abag_q2,
                "pos_ok": pos_ok_2,
                "extra": (f"iq2={iq2_peak_slice_r}, Q² ≈ {q2_slice_val:.2f} GeV², "
                          f"N_slice = {N_q2slice:.2e}"),
            })

            print()

            # ── Case 3: Fully integrated (RECO) ───────────────────────────────
            A_integrated, N_integrated = compute_A_cos2phi_integrated_reco(
                n_reco, A_reco, reco_mask
            )
            N_flat_integrated = N_integrated / N_PHI
            mod_at_abag_int    = abs(PZZ * A_integrated)
            mod_at_10abag_int  = abs(PZZ * 10.0 * A_integrated)

            print(f"  Case 3 (fully integrated, reco):")
            print(f"    Total N_reco over all reco bins = {N_integrated:.2e} events")
            print(f"    ⟨A_cos2φ⟩ (reco) = {A_integrated:.4e}")
            print(f"    P_zz·⟨A_cos2φ⟩ (at A_bag) = {mod_at_abag_int:.4f}")
            print(f"    Max modulation at 10×A_bag = {mod_at_10abag_int:.4f}")

            pos_ok_3 = mod_at_10abag_int < 1.0
            print(f"    Positivity check: {'OK' if pos_ok_3 else 'WARN'}")

            outpath = build_phi_plot(
                A_ref=A_integrated,
                N_flat=N_flat_integrated,
                A_bag_config=a_bag_config,
                config_label=config_label,
                case_label="fully integrated",
                case_tag="integrated",
                config_tag=config_tag,
                outdir=outdir,
                extra_info=None,
            )
            print(f"  wrote {outpath.name}")

            summary_rows.append({
                "config": config_tag.upper(),
                "case": "Case 3 (integrated, reco)",
                "A_ref": A_integrated,
                "N_flat": N_flat_integrated,
                "mod_at_abag": mod_at_abag_int,
                "mod_at_10abag": mod_at_10abag_int,
                "pos_ok": pos_ok_3,
                "extra": f"N_total (reco) = {N_integrated:.2e}",
            })

            print()

    # ── Final summary table ───────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("SUMMARY — φ-modulation projections [detector-realistic] at L = 10 fb⁻¹/nucleon")
    print("=" * 72)
    print()

    current_config = None
    for row in summary_rows:
        if row["config"] != current_config:
            current_config = row["config"]
            a_bag_val = A_BAG[current_config.lower()]
            print(f"Config: {current_config}  (|A_bag| = {a_bag_val:.3f})")

        print(f"  {row['case']}:")
        print(f"    {row['extra']}")
        print(f"    ⟨A_cos2φ⟩ = {row['A_ref']:.4e}")
        print(f"    N_flat (events/5° bin) = {row['N_flat']:.2e}")
        print(f"    P_zz·⟨A_cos2φ⟩ at A_bag = {row['mod_at_abag']:.4f}  "
              f"(positivity at A_bag: {'OK' if row['mod_at_abag'] < 1.0 else 'WARN'})")
        print(f"    Max modulation at 10×A_bag = {row['mod_at_10abag']:.4f}  "
              f"({'OK' if row['pos_ok'] else 'WARN'})")
        print()

    print("=" * 72)
    print()
    print("Output directory:", outdir)
    print()
    print("Files written (16 total, _20260729_ prefix):")
    # 3 per-bin heatmap motivation plots
    for cfg_tag in ["low", "mid", "top"]:
        print(f"  money_delta_20260729_perbin_{cfg_tag}.png")
    # 8 φ plots: LOW+MID all 3 cases, TOP cases 2+3
    for case in ["peakbin", "q2slice", "integrated"]:
        for cfg_tag in ["low", "mid"]:
            print(f"  money_delta_20260729_phimodulation_{case}_{cfg_tag}.png")
    for case in ["q2slice", "integrated"]:
        print(f"  money_delta_20260729_phimodulation_{case}_top.png")
    # 5 TOP+peakbin bin-width scan plots
    for bin_deg in [5, 10, 20, 30, 45]:
        print(f"  money_delta_20260729_phimodulation_peakbin_top_{bin_deg}deg.png")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
