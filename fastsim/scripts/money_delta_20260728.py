#!/usr/bin/env python3
"""6Li φ-distribution plots — fractional yield modulation dN/dφ vs φ.

Produces 12 individual PNG plots showing the observable φ-modulation that
an experimenter would see in the cos(2φ) tensor asymmetry measurement, at
three signal scales (0.1·A_bag, A_bag, 10·A_bag), for three EIC beam
configurations, and three levels of kinematic integration.

Physics setup (fixed):
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

Three cases per config:
  Case 1: Peak (x, Q²) bin
  Case 2: Peak Q² slice (rate-weighted ⟨A_cos2φ⟩ over x)
  Case 3: Fully integrated (rate-weighted ⟨A_cos2φ⟩ over all accepted bins)

Output (16 PNGs):
  Per-bin (x, Q²) heatmaps — motivation for the φ plots (3 plots):
    money_delta_20260728_perbin_low.png
    money_delta_20260728_perbin_mid.png
    money_delta_20260728_perbin_top.png
  LOW+MID cases 1-3 and TOP cases 2-3 (8 plots, 5° bins):
    money_delta_20260728_phimodulation_{case}_{config}.png
  TOP case 1 φ-bin-width scan (5 plots):
    money_delta_20260728_phimodulation_peakbin_top_5deg.png
    money_delta_20260728_phimodulation_peakbin_top_10deg.png
    money_delta_20260728_phimodulation_peakbin_top_20deg.png
    money_delta_20260728_phimodulation_peakbin_top_30deg.png
    money_delta_20260728_phimodulation_peakbin_top_45deg.png
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
# (copied verbatim from money_delta_20260724.py)
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
# (copied verbatim from money_delta_20260724.py)
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
# Δ with α_s(Q²) prefactor (copied verbatim from money_delta_20260724.py)
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
# sig² helper (used for the sig²-per-fb build — kept for possible reference;
# not required by the φ-plot builder but matches parent's API)
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
# Case-specific helper functions
# ══════════════════════════════════════════════════════════════════════════════

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
    accepted_flat = proj.accepted.ravel()
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
        rf" {case_label}{bin_info}"
    )
    title_line2 = (
        r"$L = 10\,\mathrm{fb}^{-1}/\mathrm{u}$, mid$\_x$ shape, EPPS21, "
        r"R1998, Cloet $P_{zz}$=0.267"
    )
    ax.set_title(title_line1 + "\n" + title_line2, fontsize=9)

    # Legend
    ax.legend(fontsize=7.5, loc="upper right", framealpha=0.8)

    # Text annotation in upper-left corner
    ann_lines = [
        rf"$\langle A_{{\cos 2\phi}}\rangle$ = {A_ref:.4f}",
        rf"$N_\mathrm{{flat}}$ = {N_flat:.2e} events/bin ({bin_width_deg}°)",
        rf"$|A_\mathrm{{bag}}|$ = {A_bag_config:.3f}",
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
        filename = f"money_delta_20260728_phimodulation_{case_tag}{filename_suffix}.png"
    else:
        filename = f"money_delta_20260728_phimodulation_{case_tag}_{config_tag}.png"
    outpath = outdir / filename
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return outpath


# ══════════════════════════════════════════════════════════════════════════════
# Per-bin (x, Q²) heatmap — motivation plots
# ══════════════════════════════════════════════════════════════════════════════

def build_perbin_heatmap(cfg, base, config_tag, config_label, A_bag_config,
                         luminosity, outdir, proj=None):
    """Build per-bin 1×3 subplot heatmap: |A_bin|, δA_bin, |δA_bin/A_bin|.

    Interpration A: Δ = s · α_s(Q²) · F₁ · x^α · (1-x)^β.
    Mirrors money_delta_20260725.py's build_perbin_plot structurally.

    Parameters
    ----------
    cfg            : BeamConfig
    base           : PartonF2   CT18NLO backend (for α_s)
    config_tag     : str        'low', 'mid', or 'top'
    config_label   : str        human-readable beam description
    A_bag_config   : float      |A_bag| for this config (positive)
    luminosity     : float      integrated luminosity [fb⁻¹/nucleon]
    outdir         : pathlib.Path
    proj           : BinnedProjection or None
                     If provided, reused (avoids a redundant project_rates call).
                     If None, a fresh project_rates call is made.

    Returns
    -------
    pathlib.Path : output PNG path
    """
    L = luminosity

    # ── Run / reuse project_rates ─────────────────────────────────────────────
    if proj is None:
        nf2_obj = NuclearF2FromGrid(cfg.ion, EPPS21_SET)
        sc = fom.Scenario(lumi_fb_per_nucleon=L, pol_ion_tensor=PZZ, q2_min=2.0)
        proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_obj)

    n_events = proj.n_events   # shape (nx, nq2)
    accepted = proj.accepted   # bool mask

    nf2 = proj.extras["nf2"]
    # Per-nucleon structure functions
    f1 = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2 = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y  = proj.extras["y"]

    # ── Per-bin physics (Interpretation A) ───────────────────────────────────
    # Δ = A_bag_config · α_s(Q²) · F₁ · x^0.7 · (1-x)^3  (mid_x shape)
    delta_bin = delta_shape_with_alphas(
        proj.x, proj.q2, f1,
        scale=A_bag_config, variant=VARIANT, base=base
    )
    # A_cos2φ per bin (vectorised over 2-D grid)
    # a_cos2phi expects scalars or matching arrays; use vectorised form
    amp_bin = np.array([
        float(a_cos2phi(delta_bin[ix, iq], f1[ix, iq], f2[ix, iq],
                        proj.x[ix, iq], y[ix, iq]))
        for ix in range(proj.x.shape[0])
        for iq in range(proj.x.shape[1])
    ]).reshape(proj.x.shape)

    # ── Validity mask ─────────────────────────────────────────────────────────
    use = accepted & (n_events >= MIN_EVENTS)
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
        fontsize=10,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$", fontsize=10)
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]", fontsize=10)
    ax.set_title("|δA_bin / A_bin| (relative)", fontsize=10)

    fig.suptitle(
        f"$^6$Li per-bin decomposition — {config_label} config, "
        f"L = {L:g} fb$^{{-1}}$/u\n"
        f"Interpretation A: "
        r"$\Delta = s \cdot \alpha_s(Q^2) \cdot F_1 \cdot x^\alpha (1-x)^\beta$"
        f", $s = |A_\\mathrm{{bag}}|$ = {A_bag_config:.3f}, "
        "mid_x shape, EPPS21, R1998, Cloet $P_{{zz}}$=0.267",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"money_delta_20260728_perbin_{config_tag}.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath.name}")

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
    ix_r, iq_r = np.unravel_index(flat_idx, rel_unc.shape)
    x_mid_r  = (0.5 * (x_edges[ix_r]  + x_edges[ix_r  + 1])
                if ix_r  < len(x_edges)  - 1 else x_edges[ix_r])
    q2_mid_r = (0.5 * (q2_edges[iq_r] + q2_edges[iq_r + 1])
                if iq_r < len(q2_edges) - 1 else q2_edges[iq_r])
    print(f"  Location of max |δA/A|:  x={x_mid_r:.2e},  Q²={q2_mid_r:.2f} GeV²")

    # Peak-rate bin (for cross-ref with Case 1)
    n_bin_map = np.where(mask_valid, n_events, np.nan)
    flat_peak = np.nanargmax(n_bin_map)
    ix_p, iq_p = np.unravel_index(flat_peak, n_bin_map.shape)
    x_pk  = float(proj.x[ix_p, iq_p])
    q2_pk = float(proj.q2[ix_p, iq_p])
    a_pk  = float(amp_abs[ix_p, iq_p]) if not np.isnan(amp_abs[ix_p, iq_p]) else float("nan")
    print(f"  Peak-rate bin (for cross-ref with Case 1): "
          f"x={x_pk:.2e}, Q²={q2_pk:.2f}, |A_bin|={a_pk:.2e}")

    return outpath


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description=(
            "6Li φ-distribution plots — fractional yield modulation dN/dφ vs φ.\n\n"
            "Produces 9 individual PNG plots (3 cases × 3 configs) showing the\n"
            "observable φ-modulation at three signal scales:\n"
            "  0.1·A_bag, A_bag, 10·A_bag.\n\n"
            "Cases:\n"
            "  peakbin:    Peak (x, Q²) bin\n"
            "  q2slice:    Peak Q² slice (rate-weighted ⟨A_cos2φ⟩ over x)\n"
            "  integrated: Fully integrated (rate-weighted ⟨A_cos2φ⟩ over all bins)\n\n"
            "Configs:\n"
            "  low: 5 GeV e × 27.5 GeV/u ⁶Li\n"
            "  mid: 10 GeV e × 50 GeV/u ⁶Li\n"
            "  top: 18 GeV e × 137.5 GeV/u ⁶Li\n\n"
            "Physics: Δ(x,Q²) = s · α_s(Q²) · F₁(x,Q²) · x^α (1-x)^β (mid_x shape).\n"
            "A_bag hardcoded from money_delta_20260724.py mid_x output:\n"
            "  LOW=0.318, MID=0.310, TOP=0.297"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
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
    print("money_delta_20260728.py — 6Li φ-modulation projections")
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
    # summary_rows: list of dicts, one per (config, case)
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

            # Accepted events grid info
            n_total_acc = float(proj.n_events[proj.accepted].sum())
            n_acc_bins  = int(proj.accepted.sum())
            print(f"  Accepted bins: {n_acc_bins}  |  Total events: {n_total_acc:.3e}")
            print()

            # ── Heatmap: per-bin motivation plot (before Case 1/2/3 φ plots) ──
            build_perbin_heatmap(
                cfg=cfg,
                base=base,
                config_tag=config_tag,
                config_label=config_label,
                A_bag_config=a_bag_config,
                luminosity=LUMI_FB,
                outdir=outdir,
                proj=proj,   # reuse already-computed projection
            )
            print()

            # ── Case 1: Peak bin ──────────────────────────────────────────────
            ix_peak, iq2_peak = find_peak_bin(proj)
            n_peak = float(proj.n_events[ix_peak, iq2_peak])
            x_peak_val  = float(proj.x[ix_peak, iq2_peak])
            q2_peak_val = float(proj.q2[ix_peak, iq2_peak])

            A_peakbin = compute_A_cos2phi_at_bin(
                proj, nf2_obj, cfg, ix_peak, iq2_peak, a_bag_config, base
            )
            N_flat_peakbin = n_peak / N_PHI
            mod_at_abag = abs(PZZ * A_peakbin)
            mod_at_10abag = abs(PZZ * (10.0 * a_bag_config / a_bag_config) * A_peakbin)

            print(f"  Case 1 (peak bin):")
            print(f"    x = {x_peak_val:.2e},  Q² = {q2_peak_val:.2f} GeV²  (ix={ix_peak}, iq2={iq2_peak})")
            print(f"    N_bin = {n_peak:.2e} events")
            print(f"    A_cos2φ = {A_peakbin:.4e}")
            print(f"    P_zz·A_cos2φ (at A_bag) = {mod_at_abag:.4f}")
            print(f"    Max modulation at 10×A_bag = {mod_at_10abag:.4f}")

            pos_ok_1 = mod_at_10abag < 1.0
            print(f"    Positivity check: {'OK' if pos_ok_1 else 'WARN'}")

            extra_peakbin = f"peak bin: x={x_peak_val:.2e}, Q²={q2_peak_val:.2f} GeV²"

            if config_tag == "top":
                # ── TOP + Case 1 only: scan over φ-bin widths ─────────────────
                print()
                print(f"  → TOP config, Case 1: scanning φ bin widths {{5°, 10°, 20°, 30°, 45°}}")
                for bin_deg in [5, 10, 20, 30, 45]:
                    n_bins = int(360 / bin_deg)
                    n_flat_scan = n_peak / n_bins
                    err_scan = 1.0 / np.sqrt(n_flat_scan) if n_flat_scan > 0 else float("nan")
                    print(
                        f"    {bin_deg}° bins ({n_bins:2d} total): "
                        f"N_flat = {n_flat_scan:.2e} events/bin, "
                        f"error bar ≈ {err_scan:.2e}"
                    )
                print()
                phi_bin_paths = []
                for bin_deg in [5, 10, 20, 30, 45]:
                    n_bins = int(360 / bin_deg)
                    n_flat_scan = n_peak / n_bins
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
                # All other configs: single 5° plot with original filename
                build_phi_plot(
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
                print(f"  wrote money_delta_20260728_phimodulation_peakbin_{config_tag}.png")

            summary_rows.append({
                "config": config_tag.upper(),
                "case": "Case 1 (peak bin)",
                "A_ref": A_peakbin,
                "N_flat": N_flat_peakbin,
                "mod_at_abag": mod_at_abag,
                "mod_at_10abag": mod_at_10abag,
                "pos_ok": pos_ok_1,
                "extra": f"x = {x_peak_val:.2e}, Q² = {q2_peak_val:.2f} GeV², N_bin = {n_peak:.2e}",
            })

            print()

            # ── Case 2: Peak Q² slice ─────────────────────────────────────────
            iq2_peak_slice = find_peak_q2_slice(proj)
            q2_slice_val = float(proj.q2[0, iq2_peak_slice])  # Q² value for this slice

            A_q2slice, N_q2slice = compute_A_cos2phi_q2slice(
                proj, nf2_obj, cfg, iq2_peak_slice, a_bag_config, base
            )
            N_flat_q2slice = N_q2slice / N_PHI
            mod_at_abag_q2 = abs(PZZ * A_q2slice)
            mod_at_10abag_q2 = abs(PZZ * 10.0 * A_q2slice)

            print(f"  Case 2 (peak Q² slice):")
            print(f"    Peak Q² slice: iq2={iq2_peak_slice}, Q² ≈ {q2_slice_val:.2f} GeV²")
            print(f"    Sum N_bin over x = {N_q2slice:.2e} events")
            print(f"    ⟨A_cos2φ⟩ = {A_q2slice:.4e}")
            print(f"    P_zz·⟨A_cos2φ⟩ (at A_bag) = {mod_at_abag_q2:.4f}")
            print(f"    Max modulation at 10×A_bag = {mod_at_10abag_q2:.4f}")

            pos_ok_2 = mod_at_10abag_q2 < 1.0
            print(f"    Positivity check: {'OK' if pos_ok_2 else 'WARN'}")

            extra_q2slice = f"peak Q² slice: Q²≈{q2_slice_val:.2f} GeV²"
            build_phi_plot(
                A_ref=A_q2slice,
                N_flat=N_flat_q2slice,
                A_bag_config=a_bag_config,
                config_label=config_label,
                case_label="peak $Q^2$ slice",
                case_tag="q2slice",
                config_tag=config_tag,
                outdir=outdir,
                extra_info=extra_q2slice,
            )

            summary_rows.append({
                "config": config_tag.upper(),
                "case": "Case 2 (Q² slice)",
                "A_ref": A_q2slice,
                "N_flat": N_flat_q2slice,
                "mod_at_abag": mod_at_abag_q2,
                "mod_at_10abag": mod_at_10abag_q2,
                "pos_ok": pos_ok_2,
                "extra": f"iq2={iq2_peak_slice}, Q² ≈ {q2_slice_val:.2f} GeV², N_slice = {N_q2slice:.2e}",
            })

            print()

            # ── Case 3: Fully integrated ──────────────────────────────────────
            A_integrated, N_integrated = compute_A_cos2phi_integrated(
                proj, nf2_obj, cfg, a_bag_config, base
            )
            N_flat_integrated = N_integrated / N_PHI
            mod_at_abag_int = abs(PZZ * A_integrated)
            mod_at_10abag_int = abs(PZZ * 10.0 * A_integrated)

            print(f"  Case 3 (fully integrated):")
            print(f"    Total N_bin over all accepted bins = {N_integrated:.2e} events")
            print(f"    ⟨A_cos2φ⟩ = {A_integrated:.4e}")
            print(f"    P_zz·⟨A_cos2φ⟩ (at A_bag) = {mod_at_abag_int:.4f}")
            print(f"    Max modulation at 10×A_bag = {mod_at_10abag_int:.4f}")

            pos_ok_3 = mod_at_10abag_int < 1.0
            print(f"    Positivity check: {'OK' if pos_ok_3 else 'WARN'}")

            build_phi_plot(
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

            summary_rows.append({
                "config": config_tag.upper(),
                "case": "Case 3 (integrated)",
                "A_ref": A_integrated,
                "N_flat": N_flat_integrated,
                "mod_at_abag": mod_at_abag_int,
                "mod_at_10abag": mod_at_10abag_int,
                "pos_ok": pos_ok_3,
                "extra": f"N_total = {N_integrated:.2e}",
            })

            print()

    # ── Final summary table ───────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("SUMMARY — φ-modulation projections at L = 10 fb⁻¹/nucleon")
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
    print("Files written (16 total):")
    # 3 per-bin heatmap motivation plots
    for cfg_tag in ["low", "mid", "top"]:
        print(f"  money_delta_20260728_perbin_{cfg_tag}.png")
    # 8 unchanged plots: LOW+MID all 3 cases, TOP cases 2+3
    for case in ["peakbin", "q2slice", "integrated"]:
        for cfg_tag in ["low", "mid"]:
            print(f"  money_delta_20260728_phimodulation_{case}_{cfg_tag}.png")
    for case in ["q2slice", "integrated"]:
        print(f"  money_delta_20260728_phimodulation_{case}_top.png")
    # 5 TOP+peakbin bin-width scan plots
    for bin_deg in [5, 10, 20, 30, 45]:
        print(f"  money_delta_20260728_phimodulation_peakbin_top_{bin_deg}deg.png")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
