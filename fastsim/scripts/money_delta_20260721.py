#!/usr/bin/env python3
"""6Li statistical uncertainty on A_cos2φ vs Δ/F₁ — EPPS21 & nNNPDF3.0 comparison.

Differences from money_delta_20260720.py
──────────────────────────────────────────
1. **Observable changed**: instead of L_5σ (required luminosity for 5σ
   discovery), this script shows **δA/A** — the relative statistical
   uncertainty on the cos(2φ) amplitude — at **fixed luminosity**.

2. **Two beam configs**: both mid (E_e = 10 GeV × p_ion = 50 GeV/u ⁶Li) and
   top (E_e = 18 GeV × p_ion = 137.5 GeV/u ⁶Li) configs are produced.

3. **Two luminosities**: 10 fb⁻¹/nucleon and 100 fb⁻¹/nucleon.

4. **Eight output PNGs**: four "stat-uncertainty vs Δ/F₁" plots (Plot A,
   log-log), and four "per-(x,Q²)-bin stat uncertainty" heatmaps (Plot B).

5. **No gold luminosity band**: the luminosity is now a fixed input, so the
   "plausible EIC program" axhspan from the parent is not shown.

6. **5σ threshold line**: Plot A shows a horizontal line at δA/A = 0.20,
   labeled "5σ discovery". Crossing this threshold at a given Δ/F₁ scale is
   equivalent to L_5σ = L in the parent money plot — a self-consistency check
   that is reported in the printed summary.

All physics helpers (NuclearF2FromGrid, r_override, sig2_per_fb_at_nuc, etc.)
are copied verbatim from money_delta_20260720.py. Analytic rescaling for
Plot A reuses the sig² per fb⁻¹/u values (computed at reference scale S0 = 1e-3)
exactly as in the parent, but converts them to δA/A via

    δA/A = 1 / sqrt( sig2_per_fb × (s/S0)² × L )

which is equivalent to Nσ⁻¹ at luminosity L.

Instantaneous luminosity convention (2026-07-21 note)
─────────────────────────────────────────────────────
Rate-per-hour plots (Plot C) assume INSTANT_LUMI_FB_PER_S = 1e-6 fb⁻¹/s per
nucleon, corresponding to ~10 fb⁻¹/year at 10⁷ s of physics running.  This
matches the "10 fb⁻¹ = 1 year" convention used in the 2026-07-21 note.

Output
──────
fastsim/out/money_delta/money_delta_20260721_statvsdelta_mid_10fb.png
fastsim/out/money_delta/money_delta_20260721_statvsdelta_mid_100fb.png
fastsim/out/money_delta/money_delta_20260721_statvsdelta_top_10fb.png
fastsim/out/money_delta/money_delta_20260721_statvsdelta_top_100fb.png
fastsim/out/money_delta/money_delta_20260721_perbin_mid_10fb.png
fastsim/out/money_delta/money_delta_20260721_perbin_mid_100fb.png
fastsim/out/money_delta/money_delta_20260721_perbin_top_10fb.png
fastsim/out/money_delta/money_delta_20260721_perbin_top_100fb.png
fastsim/out/money_delta/money_delta_20260721_rate_perbin_mid.png   (new: EPPS21, plasma cmap)
fastsim/out/money_delta/money_delta_20260721_rate_perbin_top.png   (new: EPPS21, plasma cmap)
fastsim/out/money_delta/money_delta_20260721_rate_perbin_mid.csv   (new: full per-bin rate table)
fastsim/out/money_delta/money_delta_20260721_rate_perbin_top.csv   (new: full per-bin rate table)
"""

import argparse
import csv
import datetime
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
import matplotlib.colors as mcolors


# ══════════════════════════════════════════════════════════════════════════════
# NuclearF2FromGrid — local class; do NOT modify polli_fastsim/structure.py
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
# R2: Three Δ shape variants (copied verbatim from money_delta_20260720.py)
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
# R4: Two R = σ_L/σ_T parameterizations (copied verbatim from parent)
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
# Fixed P_zz: Cloet convention
# ══════════════════════════════════════════════════════════════════════════════

# P_zz is fixed at the Cloet per-nucleon-normalized value.
# See beams.py docstring on "6Li: UNRESOLVED convention (factor 2.4!)".
PZZ = 0.267   # = 0.8 / 3: per-nucleon-normalized, Cloet convention

# toy baseline P_zz (kept for continuity with money_delta_realistic.py)
TOY_PZZ = 0.80


# ── Instantaneous luminosity assumption (per nucleon) ─────────────────────────
# 1e-6 fb^-1/s per nucleon corresponds to ~10 fb^-1/year at 10^7 s of physics
# running, matching the "10 fb^-1 = 1 year" convention used in the 2026-07-21 note.
# Per-nucleon; for 6Li this is roughly 1.67e-7 fb^-1/s per ion, i.e. ~10^32
# cm^-2 s^-1 per nucleon, within a factor of a few of EIC design e-A luminosity.
INSTANT_LUMI_FB_PER_S = 1e-6                          # fb^-1/nucleon/s
LUMI_PER_HOUR_FB     = INSTANT_LUMI_FB_PER_S * 3600.0 # fb^-1/nucleon/hour = 3.6e-3


# ══════════════════════════════════════════════════════════════════════════════
# Core significance helpers (verbatim from money_delta_20260720.py)
# ══════════════════════════════════════════════════════════════════════════════

def sig2_per_fb_at_nuc(cfg, scale, pzz, nf2_obj, variant="mid_x",
                        min_events=10):
    """Sig² per fb⁻¹/nucleon using an arbitrary NuclearF2-compatible object.

    Parameters
    ----------
    cfg      : BeamConfig
    scale    : float  peak Δ/F₁ value (reference scale)
    pzz      : float  effective P_zz per nucleon
    nf2_obj  : NuclearF2 | NuclearF2FromGrid  pre-constructed F₂ backend
    variant  : str    Δ shape variant ('low_x', 'mid_x', 'high_x')
    min_events : int  minimum events to include a bin

    The R-override must be active in the calling context (R4).
    """
    sc = fom.Scenario(lumi_fb_per_nucleon=1.0, pol_ion_tensor=pzz, q2_min=2.0)
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_obj)

    nf2 = proj.extras["nf2"]
    f1 = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2 = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y  = proj.extras["y"]

    delta = delta_shape(proj.x, proj.q2, f1, scale=scale, variant=variant)
    amp   = a_cos2phi(delta, f1, f2, proj.x, y)

    n_events = proj.n_events
    use = proj.accepted & (n_events >= min_events)
    return np.where(use, amp**2 * pzz**2 * n_events / 2.0, 0.0).sum()


def sig2_per_fb_at_ct18emc(cfg, scale, pzz, base, variant="mid_x",
                             min_events=10):
    """Sig² for the CT18NLO + EMC reference (matching money_delta_realistic.py).

    This is identical to sig2_per_fb_at_realistic in the parent script.
    Used for the gray dashed reference line only.
    """
    nf2_obj = NuclearF2(cfg.ion, base=base, emc_ratio=unpolarized_emc_ratio)
    return sig2_per_fb_at_nuc(cfg, scale, pzz, nf2_obj, variant=variant,
                               min_events=min_events)


def sig2_per_fb_toy(cfg, scale, pzz, base, min_events=10):
    """Toy-backend sig² for continuity comparison (no EMC, toy R, toy Δ shape).

    Mirrors sig2_per_fb_toy in money_delta_realistic.py.
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
# Grid configuration
# ══════════════════════════════════════════════════════════════════════════════

# x-axis: 15 log-spaced Δ/F₁ scale points (same as parent)
SCALES = np.logspace(-3.3, -1.7, 15)
S0 = 1e-3          # reference scale for the analytic rescaling

R_FUNCS = [
    ("R1998",          r1998),
    ("Christy-Bosted", r_christy_bosted),
]

VARIANT_NAMES = ["low_x", "mid_x", "high_x"]

# Central combination: mid_x shape + R1998 + EPPS21
CENTRAL_VARIANT = "mid_x"
CENTRAL_R_NAME  = "R1998"

# Nuclear PDF sets for the two bands
NUCLEAR_SETS = [
    ("EPPS21",    "EPPS21nlo_CT18Anlo_Li6"),
    ("nNNPDF3.0", "nNNPDF30_nlo_as_0118_A6_Z3"),
]

# Luminosities to evaluate
LUMINOSITIES_FB = [10.0, 100.0]   # fb⁻¹/nucleon

# Min events cut (same as parent)
MIN_EVENTS = 10


# ══════════════════════════════════════════════════════════════════════════════
# Plot A — stat uncertainty vs Δ/F₁ (log-log)
# ══════════════════════════════════════════════════════════════════════════════

def build_stat_vs_delta_plot(cfg, sig2_dict, luminosity, config_tag,
                              config_label, outdir):
    """Build one Plot A: δA/A vs Δ/F₁ at fixed luminosity.

    Parameters
    ----------
    cfg          : BeamConfig
    sig2_dict    : dict with keys:
                     combo_sig2["EPPS21"][(variant, r_name)]  → sig² at S0
                     combo_sig2["nNNPDF3.0"][(variant, r_name)]
                     "sig2_ref"   → float (CT18+EMC reference at S0)
                     "sig2_toy"   → float (toy backend at S0)
    luminosity   : float  integrated luminosity [fb⁻¹/nucleon]
    config_tag   : str    e.g. 'mid' or 'top'
    config_label : str    human-readable beam description
    outdir       : pathlib.Path

    Returns
    -------
    pathlib.Path : output PNG path
    dict         : summary numbers
    """
    L = luminosity
    lumi_tag = f"{int(L)}fb" if L == int(L) else f"{L:g}fb"

    def sig2_to_rel_unc(sig2_at_s0_per_fb, scales, lumi):
        """Convert sig² per fb at S0 → δA/A(scales) at fixed lumi."""
        sig2_at_scale = sig2_at_s0_per_fb * (scales / S0) ** 2 * lumi
        # Avoid division by zero: if sig2 = 0 → δA/A = inf (undetectable)
        return 1.0 / np.sqrt(np.maximum(sig2_at_scale, 1e-60))

    # ── EPPS21 band ──────────────────────────────────────────────────────────
    epps_rel_curves = []
    for (variant, r_name), s2 in sig2_dict["EPPS21"].items():
        epps_rel_curves.append(sig2_to_rel_unc(s2, SCALES, L))
    epps_rel_curves = np.array(epps_rel_curves)
    epps_min = epps_rel_curves.min(axis=0)   # smallest δA/A = most sensitive
    epps_max = epps_rel_curves.max(axis=0)   # largest δA/A = least sensitive

    # ── nNNPDF3.0 band ───────────────────────────────────────────────────────
    nnpdf_rel_curves = []
    for (variant, r_name), s2 in sig2_dict["nNNPDF3.0"].items():
        nnpdf_rel_curves.append(sig2_to_rel_unc(s2, SCALES, L))
    nnpdf_rel_curves = np.array(nnpdf_rel_curves)
    nnpdf_min = nnpdf_rel_curves.min(axis=0)
    nnpdf_max = nnpdf_rel_curves.max(axis=0)

    # ── Central curve: mid_x + R1998 + EPPS21 ────────────────────────────────
    sig2_central = sig2_dict["EPPS21"][(CENTRAL_VARIANT, CENTRAL_R_NAME)]
    rel_central  = sig2_to_rel_unc(sig2_central, SCALES, L)

    # ── CT18+EMC reference (gray dashed) ─────────────────────────────────────
    rel_ref = sig2_to_rel_unc(sig2_dict["sig2_ref"], SCALES, L)

    # ── Toy backend (gray dotted) ─────────────────────────────────────────────
    rel_toy = sig2_to_rel_unc(sig2_dict["sig2_toy"], SCALES, L)

    # ── Draw ─────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 6))

    # EPPS21 band (steelblue)
    ax.fill_between(
        SCALES, epps_min, epps_max,
        color="steelblue", alpha=0.22,
        label="EPPS21 band (3 shapes × 2 R)",
    )

    # nNNPDF3.0 band (darkorange)
    ax.fill_between(
        SCALES, nnpdf_min, nnpdf_max,
        color="darkorange", alpha=0.22,
        label="nNNPDF3.0 band (3 shapes × 2 R)",
    )

    # Central curve (black solid)
    ax.plot(
        SCALES, rel_central, "-", color="black", lw=2.2,
        label=r"Central: mid-$x$ $\Delta$, R1998, EPPS21, Cloet $P_{zz}=0.267$",
    )

    # CT18+EMC reference (gray dashed)
    ax.plot(
        SCALES, rel_ref, "--", color="gray", lw=1.6,
        label=r"CT18NLO + EMC reference (mid-$x$, R1998)",
    )

    # Toy-backend comparison (gray dotted)
    ax.plot(
        SCALES, rel_toy, ":", color="gray", lw=1.3,
        label=r"Toy backend, $P_{zz}=0.8$ (continuity ref.)",
    )

    # 5σ discovery threshold: δA/A = 1/5 = 0.20
    ax.axhline(
        0.20, color="firebrick", ls="--", lw=1.4, zorder=3,
        label=r"5$\sigma$ discovery ($\delta A/A = 0.20$)",
    )

    # Sather-Schmidt reference line (copied from parent)
    ax.axvline(1e-3, color="dimgray", ls=":", lw=1.2)
    ax.text(
        1.08e-3, 0.90, "Sather-Schmidt\n$O(10^{-3})$",
        fontsize=7,
        transform=ax.get_xaxis_transform(),
        va="top", color="dimgray",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta/F_1$ scale (peak of shape)", fontsize=11)
    ax.set_ylabel(r"$\delta A / A_{\cos 2\phi}$ (relative stat. uncertainty)",
                  fontsize=11)
    ax.set_title(
        f"$^6$Li stat uncertainty on $A_{{\\cos 2\\phi}}$ vs "
        r"$\Delta/F_1$"
        f" — {config_label}, $L = {L:g}$ fb$^{{-1}}$/u\n"
        r"EPPS21 vs nNNPDF3.0 nuclear PDFs, Cloet $P_{zz}=0.267$",
        fontsize=10,
    )
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"money_delta_20260721_statvsdelta_{config_tag}_{lumi_tag}.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")

    # ── Summary numbers ───────────────────────────────────────────────────────
    # 5σ crossing: find scale where central curve crosses δA/A = 0.20
    # δA/A = 0.20 ↔ sig2_central × (s*/S0)² × L = 25
    # → s* = S0 × sqrt(25 / (sig2_central × L))
    s_star = S0 * np.sqrt(25.0 / max(sig2_central * L, 1e-30))

    # EPPS21 band at S0 (as δA/A values)
    epps_s2_vals = list(sig2_dict["EPPS21"].values())
    epps_rel_at_s0 = [1.0 / np.sqrt(max(s2 * L, 1e-60)) for s2 in epps_s2_vals]
    epps_rel_lo = min(epps_rel_at_s0)   # smallest δA/A at S0
    epps_rel_hi = max(epps_rel_at_s0)   # largest δA/A at S0

    # Derived L_5σ cross-check: L_5σ(1e-3) from parent = L × (s*/1e-3)²
    l5_derived = L * (s_star / S0) ** 2

    summary = {
        "config_tag":    config_tag,
        "config_label":  config_label,
        "luminosity":    L,
        "s_star":        s_star,
        "epps_rel_lo":   epps_rel_lo,
        "epps_rel_hi":   epps_rel_hi,
        "l5_derived":    l5_derived,
        "sig2_central":  sig2_central,
    }
    return outpath, summary


# ══════════════════════════════════════════════════════════════════════════════
# Plot B — per-(x, Q²)-bin stat uncertainty heatmap
# ══════════════════════════════════════════════════════════════════════════════

def build_per_bin_plot(cfg, luminosity, config_tag, config_label, outdir,
                       epps21_nf2_obj):
    """Build one Plot B: per-bin δA_bin heatmap on the (x, Q²) plane.

    Uses EPPS21 nuclear PDF (nuclear-PDF choice negligible; see §2.1 of the
    2026-07-20 note).  P_zz = 0.267 (Cloet convention).

    Parameters
    ----------
    cfg           : BeamConfig
    luminosity    : float  integrated luminosity [fb⁻¹/nucleon]
    config_tag    : str    e.g. 'mid' or 'top'
    config_label  : str    human-readable beam description
    outdir        : pathlib.Path
    epps21_nf2_obj : NuclearF2FromGrid  pre-constructed EPPS21 object (R1998 active)

    Returns
    -------
    pathlib.Path : output PNG path
    dict         : summary numbers (per-bin δA statistics)
    """
    L = luminosity
    lumi_tag = f"{int(L)}fb" if L == int(L) else f"{L:g}fb"

    # Run project_rates at the target luminosity, EPPS21, R1998 (active outside)
    sc = fom.Scenario(lumi_fb_per_nucleon=L, pol_ion_tensor=PZZ, q2_min=2.0)
    proj = fom.project_rates(cfg, sc, nuclear_f2=epps21_nf2_obj)

    n_events = proj.n_events          # shape (nx, nq2)
    accepted = proj.accepted          # bool mask

    # Per-bin δA_bin = sqrt(2 / N_bin) / P_zz  (from asymmetries.err_cos2phi_amplitude)
    # Only bins with N_bin >= MIN_EVENTS inside acceptance are shown; rest → NaN
    use = accepted & (n_events >= MIN_EVENTS)

    delta_a_bin = np.full_like(n_events, np.nan)
    delta_a_bin[use] = np.sqrt(2.0 / n_events[use]) / PZZ

    # Bin edges for pcolormesh (log scale)
    x_edges  = proj.x_edges    # shape (nx+1,)
    q2_edges = proj.q2_edges   # shape (nq2+1,)

    # pcolormesh needs the grid transposed: (nq2, nx) with edges as 2D meshes
    X_edge_2d, Q2_edge_2d = np.meshgrid(x_edges, q2_edges)  # (nq2+1, nx+1)

    # delta_a_bin has shape (nx, nq2); pcolormesh with (nq2+1, nx+1) edges
    # expects C to be (nq2, nx) → transpose
    C = delta_a_bin.T   # shape (nq2, nx)

    # Determine color range from valid (non-NaN) values
    valid = C[~np.isnan(C)]
    if valid.size == 0:
        vmin, vmax = 1e-3, 1.0
    else:
        vmin = max(valid.min(), 1e-6)
        vmax = valid.max()
        if vmin >= vmax:
            vmin = vmax / 100.0

    fig, ax = plt.subplots(figsize=(8, 6))
    pcm = ax.pcolormesh(
        X_edge_2d, Q2_edge_2d, C,
        norm=mcolors.LogNorm(vmin=vmin, vmax=vmax),
        cmap="viridis",
        shading="flat",
    )
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label(r"$\delta A_{\cos 2\phi}$ per bin", fontsize=11)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$", fontsize=11)
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]", fontsize=11)
    ax.set_title(
        f"$^6$Li per-bin stat uncertainty on $A_{{\\cos 2\\phi}}$"
        f" — {config_label}, $L = {L:g}$ fb$^{{-1}}$/u\n"
        rf"EPPS21, Cloet $P_{{zz}}=0.267$, $N_{{\rm bin}} \geq {MIN_EVENTS}$",
        fontsize=10,
    )
    fig.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"money_delta_20260721_perbin_{config_tag}_{lumi_tag}.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")

    # ── Per-bin δA summary statistics ─────────────────────────────────────────
    n_total_accepted = int(accepted.sum())
    n_passing_cut    = int(use.sum())
    frac_passing     = n_passing_cut / max(n_total_accepted, 1)

    valid_da = delta_a_bin[use]   # 1D array of accepted δA_bin values
    if valid_da.size > 0:
        da_min    = float(valid_da.min())
        da_max    = float(valid_da.max())
        da_median = float(np.median(valid_da))
    else:
        da_min = da_max = da_median = np.nan

    summary = {
        "config_tag":        config_tag,
        "config_label":      config_label,
        "luminosity":        L,
        "n_total_accepted":  n_total_accepted,
        "n_passing_cut":     n_passing_cut,
        "frac_passing":      frac_passing,
        "da_min":            da_min,
        "da_max":            da_max,
        "da_median":         da_median,
    }
    return outpath, summary


# ══════════════════════════════════════════════════════════════════════════════
# Plot C — rate per hour per (x, Q²) bin heatmap
# ══════════════════════════════════════════════════════════════════════════════

def build_rate_perbin_plot(cfg, base, config_tag, config_label, outdir):
    """Build Plot C: events-per-hour per (x, Q²) bin heatmap (EPPS21, plasma).

    Runs project_rates once with Scenario luminosity set to LUMI_PER_HOUR_FB so
    that proj.n_events is directly in units of events/hour/bin.  No N ≥ 10 cut
    is applied — low-rate bins are informative for DAQ planning.

    Parameters
    ----------
    cfg          : BeamConfig
    base         : PartonF2  CT18NLO backend (not used for rates, passed for
                   API consistency with other build_* functions)
    config_tag   : str    e.g. 'mid' or 'top'
    config_label : str    human-readable beam description
    outdir       : pathlib.Path

    Returns
    -------
    BinnedProjection : projection object (reused by build_rate_table_and_summary)
    """
    # Use EPPS21 nuclear PDF (mirroring build_per_bin_plot convention)
    nf2_obj = NuclearF2FromGrid(cfg.ion, "EPPS21nlo_CT18Anlo_Li6")

    # Scenario luminosity = LUMI_PER_HOUR_FB → n_events are events/hour/bin
    sc = fom.Scenario(
        lumi_fb_per_nucleon=LUMI_PER_HOUR_FB,
        pol_ion_tensor=PZZ,
        q2_min=2.0,
    )
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_obj)

    # Rate map: events/hour/bin; bins outside acceptance → NaN (blank on plot)
    rate = np.where(proj.accepted, proj.n_events, np.nan)

    # Bin edges for pcolormesh (log scale) — same pattern as build_per_bin_plot
    x_edges  = proj.x_edges    # shape (nx+1,)
    q2_edges = proj.q2_edges   # shape (nq2+1,)

    # pcolormesh with (nq2+1, nx+1) edge meshes; C must be (nq2, nx) → transpose
    X_edge_2d, Q2_edge_2d = np.meshgrid(x_edges, q2_edges)  # (nq2+1, nx+1)
    C = rate.T   # shape (nq2, nx)

    # Color range from positive accepted values
    valid = C[np.isfinite(C) & (C > 0)]
    if valid.size == 0:
        vmin, vmax = 1e-3, 1.0
    else:
        vmin = float(valid.min())
        vmax = float(valid.max())
        if vmin >= vmax:
            vmin = vmax / 100.0

    fig, ax = plt.subplots(figsize=(8, 6))
    pcm = ax.pcolormesh(
        X_edge_2d, Q2_edge_2d, C,
        norm=mcolors.LogNorm(vmin=vmin, vmax=vmax),
        cmap="plasma",
        shading="flat",
    )
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label("Rate [events / hour / bin]", fontsize=11)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$", fontsize=11)
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]", fontsize=11)
    ax.set_title(
        f"$^6$Li DIS rate per hour per (x, Q²) bin — {config_label}\n"
        f"EPPS21 PDF, L_inst = {INSTANT_LUMI_FB_PER_S:.0e} fb$^{{-1}}$/s per nucleon "
        f"(≈ {LUMI_PER_HOUR_FB*1e3:.2f} pb$^{{-1}}$/hour)",
        fontsize=10,
    )
    fig.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"money_delta_20260721_rate_perbin_{config_tag}.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")

    return proj


def build_rate_table_and_summary(proj, config_tag, config_label, outdir):
    """Write per-bin rate CSV and print a compact stdout summary.

    Takes the BinnedProjection returned by build_rate_perbin_plot (where
    proj.n_events is in units of events/hour/bin) and:

    1. Writes a CSV with one row per (x, Q²) bin center.
    2. Prints accepted-bin statistics and the top 20 bins by rate/hour.

    Parameters
    ----------
    proj         : BinnedProjection  from build_rate_perbin_plot
    config_tag   : str    e.g. 'mid' or 'top'
    config_label : str    human-readable beam description
    outdir       : pathlib.Path
    """
    # proj.x, proj.q2 are 2D arrays of shape (nx, nq2) — flatten to 1D rows
    x_flat   = proj.x.ravel()           # shape (nx*nq2,)
    q2_flat  = proj.q2.ravel()
    y_flat   = proj.extras["y"].ravel()
    rate_flat = proj.n_events.ravel()   # events/hour/bin
    acc_flat  = proj.accepted.ravel()   # bool

    # ── Write CSV ─────────────────────────────────────────────────────────────
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / f"money_delta_20260721_rate_perbin_{config_tag}.csv"

    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )
    comment_lines = [
        f"# money_delta_20260721 per-bin rate table",
        f"# Config: {config_label}",
        f"# L_inst = 1e-6 fb^-1/s per nucleon (≈ 3.6e-3 fb^-1/hour = 3.6 pb^-1/hour)",
        f"# Generated: {timestamp}",
    ]
    header_row = ["x_center", "q2_center", "y_center", "rate_per_hour", "accepted"]

    with open(csv_path, "w", newline="") as fh:
        for line in comment_lines:
            fh.write(line + "\n")
        writer = csv.writer(fh)
        writer.writerow(header_row)
        for x_c, q2_c, y_c, rate_c, acc_c in zip(
            x_flat, q2_flat, y_flat, rate_flat, acc_flat
        ):
            writer.writerow([
                f"{x_c:.6e}",
                f"{q2_c:.6e}",
                f"{y_c:.6e}",
                f"{rate_c:.6e}",
                acc_c,
            ])
    print(f"  wrote {csv_path}")

    # ── Stdout summary ────────────────────────────────────────────────────────
    # Work only with accepted bins that have positive rate
    acc_mask      = acc_flat.astype(bool)
    pos_mask      = acc_mask & (rate_flat > 0)
    n_pos         = int(pos_mask.sum())
    total_rate    = float(rate_flat[pos_mask].sum())

    print()
    print(f"── {config_label}: rate per hour summary ──")
    print(f"Accepted bins with rate > 0: {n_pos}")
    print(f"Total integrated rate:       {total_rate:.4e} events/hour")

    # Top 20 by rate/hour among accepted bins
    accepted_indices = np.where(pos_mask)[0]
    sorted_idx = accepted_indices[np.argsort(rate_flat[accepted_indices])[::-1]]
    top_n = min(20, len(sorted_idx))

    print()
    print("Top 20 (x, Q²) bins by rate:")
    print(f"{'rank':>5}  {'x':>12}  {'Q^2':>12}  {'y':>12}  {'rate [events/hour]':>20}")
    for rank, idx in enumerate(sorted_idx[:top_n], start=1):
        print(
            f"{rank:>5}  "
            f"{x_flat[idx]:>12.3e}  "
            f"{q2_flat[idx]:>12.3e}  "
            f"{y_flat[idx]:>12.3e}  "
            f"{rate_flat[idx]:>20.3e}"
        )


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description=(
            "6Li stat uncertainty on A_cos2phi vs Δ/F₁ at fixed luminosity, "
            "comparing EPPS21 and nNNPDF3.0 nuclear PDFs. "
            "Produces 8 PNGs: four 'δA/A vs Δ/F₁' log-log plots (Plot A) and "
            "four per-(x,Q²)-bin heatmaps (Plot B), for two EIC beam configs "
            "(mid and top) at two luminosities (10 and 100 fb⁻¹/nucleon). "
            "Also produces 2 rate-per-hour heatmaps (Plot C, plasma colormap, EPPS21) "
            "and 2 per-bin rate CSV tables, one per beam config, using an instantaneous "
            "luminosity of 1e-6 fb^-1/s per nucleon (~3.6e-3 fb^-1/hour). "
            "Estimated runtime: 8–15 minutes (28 heavy grid evaluations + 2 rate evals)."
        )
    )
    ap.add_argument(
        "--outdir",
        default="fastsim/out/money_delta",
        help="Output directory for PNGs (default: fastsim/out/money_delta)",
    )
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)

    # ── Load CT18NLO backend (for CT18+EMC reference and toy lines) ──────────
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

    # ── Verify nuclear PDF sets are installed ────────────────────────────────
    try:
        from parton import mkPDF as _mkPDF
        for _, pdf_name in NUCLEAR_SETS:
            _p = _mkPDF(pdf_name, 0)
    except Exception as exc:
        print(
            "ERROR: nuclear PDF set could not be loaded.\n"
            "Install with:\n"
            "  python3 -m parton install EPPS21nlo_CT18Anlo_Li6\n"
            "  python3 -m parton install nNNPDF30_nlo_as_0118_A6_Z3\n"
            f"Underlying error: {exc}",
            file=sys.stderr,
        )
        sys.exit(2)

    # ── Two EIC beam configurations ────────────────────────────────────────
    configs = [
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

    # Collect all summaries for the final printout
    all_stata_summaries  = []   # Plot A summaries keyed by (config_tag, lumi)
    all_perbin_summaries = []   # Plot B summaries keyed by (config_tag, lumi)

    for cfg, config_tag, config_label in configs:
        print(f"\n{'='*70}")
        print(
            f"Config: {config_tag.upper()} — "
            f"e({cfg.electron_energy:.0f}) × 6Li({cfg.ion_momentum_per_nucleon:.1f}/u)"
        )
        print(f"{'='*70}")

        # ── Step 1: pre-compute sig² at S0 for all (nuclear-set, variant, R) ─
        # Total: 2 nuclear sets × 3 shapes × 2 R = 12 evaluations
        # Plus 1 CT18+EMC reference (mid_x, R1998) + 1 toy = 14 total
        combo_sig2 = {label: {} for label, _ in NUCLEAR_SETS}
        n_nuc_combos = len(NUCLEAR_SETS) * len(VARIANT_NAMES) * len(R_FUNCS)
        done = 0

        for nuc_label, pdf_set_name in NUCLEAR_SETS:
            nf2_obj = NuclearF2FromGrid(cfg.ion, pdf_set_name)
            for variant in VARIANT_NAMES:
                for r_name, r_func in R_FUNCS:
                    done += 1
                    print(
                        f"    [{done}/{n_nuc_combos}] nuc={nuc_label}, "
                        f"variant={variant}, R={r_name}, Pzz={PZZ:.3f} …",
                        end=" ", flush=True,
                    )
                    with r_override(r_func):
                        s2 = sig2_per_fb_at_nuc(
                            cfg, S0, PZZ, nf2_obj, variant=variant,
                        )
                    combo_sig2[nuc_label][(variant, r_name)] = s2
                    l5 = 25.0 / max(s2, 1e-30)
                    print(f"L_5σ(1e-3) = {l5:.2f} fb⁻¹/u")

        # ── Step 2: CT18NLO + EMC reference (mid_x, R1998) ──────────────────
        print(
            f"    [CT18+EMC ref] mid_x, R1998, CT18NLO+EMC, Pzz={PZZ:.3f} …",
            end=" ", flush=True,
        )
        with r_override(r1998):
            sig2_ref = sig2_per_fb_at_ct18emc(
                cfg, S0, PZZ, base, variant=CENTRAL_VARIANT,
            )
        print(f"L_5σ(1e-3) = {25.0 / max(sig2_ref, 1e-30):.2f} fb⁻¹/u")

        # ── Step 3: toy-backend comparison ──────────────────────────────────
        print(
            f"    [toy] toy_delta_gluon, original R, P_zz={TOY_PZZ} …",
            end=" ", flush=True,
        )
        # Use original r (not overridden)
        sig2_toy = sig2_per_fb_toy(cfg, S0, TOY_PZZ, base)
        print(f"L_5σ(1e-3) = {25.0 / max(sig2_toy, 1e-30):.2f} fb⁻¹/u")

        # Pack the sig² dict for downstream plotting functions
        sig2_dict = dict(combo_sig2)   # shallow copy of the two nuclear sets
        sig2_dict["sig2_ref"] = sig2_ref
        sig2_dict["sig2_toy"] = sig2_toy

        # ── Step 4: Build Plot A (analytic rescaling, no extra project_rates) ─
        for L in LUMINOSITIES_FB:
            print(f"\n  → Plot A: stat uncertainty vs Δ/F₁  "
                  f"[{config_tag}, L={L:g} fb⁻¹/u]")
            _, stata_sum = build_stat_vs_delta_plot(
                cfg, sig2_dict, L, config_tag, config_label, outdir,
            )
            all_stata_summaries.append(stata_sum)

        # ── Step 5: Build Plot B (one project_rates call per luminosity) ─────
        # EPPS21 nuclear F₂ with R1998 active (central convention)
        epps21_nf2 = NuclearF2FromGrid(cfg.ion, "EPPS21nlo_CT18Anlo_Li6")
        with r_override(r1998):
            for L in LUMINOSITIES_FB:
                print(f"\n  → Plot B: per-bin heatmap  "
                      f"[{config_tag}, L={L:g} fb⁻¹/u]")
                _, perbin_sum = build_per_bin_plot(
                    cfg, L, config_tag, config_label, outdir, epps21_nf2,
                )
                all_perbin_summaries.append(perbin_sum)

        # ── Step 6: Build Plot C (rate per hour per bin) + CSV table ──────────
        print(f"\n  → Plot C: rate per hour per (x, Q²) bin  [{config_tag}]")
        proj_rate = build_rate_perbin_plot(cfg, base, config_tag, config_label, outdir)
        print(f"  → CSV + stdout table  [{config_tag}]")
        build_rate_table_and_summary(proj_rate, config_tag, config_label, outdir)

    # ══════════════════════════════════════════════════════════════════════════
    # Final summary printout
    # ══════════════════════════════════════════════════════════════════════════

    print()
    print("=" * 76)
    print("SUMMARY — Plot A: δA/A vs Δ/F₁")
    print("=" * 76)

    for s in all_stata_summaries:
        L      = s["luminosity"]
        s_star = s["s_star"]
        l5_d   = s["l5_derived"]
        print()
        print(f"  Config: {s['config_label']},  L = {L:g} fb⁻¹/u")
        print(
            f"    5σ crossing (central):  Δ/F₁ = {s_star:.3e}  "
            f"(δA/A crosses 0.20)"
        )
        print(
            f"    Derived L_5σ(1e-3):     {l5_d:.2f} fb⁻¹/u  "
            f"[= L × (s*/1e-3)²; cross-check vs money_delta_20260720.py]"
        )
        print(
            f"    EPPS21 band at Δ/F₁=1e-3:  "
            f"δA/A ∈ [{s['epps_rel_lo']:.4f}, {s['epps_rel_hi']:.4f}]"
        )
        print(
            f"    Sig² (central, per fb):  {s['sig2_central']:.6g}"
        )

    print()
    print("=" * 76)
    print("SUMMARY — Plot B: per-bin δA_bin statistics")
    print("=" * 76)

    for s in all_perbin_summaries:
        L = s["luminosity"]
        print()
        print(f"  Config: {s['config_label']},  L = {L:g} fb⁻¹/u")
        print(
            f"    Accepted bins (kinematic mask): {s['n_total_accepted']}"
        )
        print(
            f"    Bins with N_bin ≥ {MIN_EVENTS}:            "
            f"{s['n_passing_cut']}  "
            f"(fraction = {s['frac_passing']:.2%})"
        )
        print(
            f"    δA_bin range:  "
            f"min = {s['da_min']:.4g},  "
            f"max = {s['da_max']:.4g},  "
            f"median = {s['da_median']:.4g}"
        )

    print()
    print("=" * 76)
    print("Note: δA/A = 1/Nσ at the given luminosity.")
    print("      5σ crossing at scale s* is equivalent to L_5σ(1e-3) = L×(s*/1e-3)²")
    print("      from money_delta_20260720.py — self-consistency cross-check.")
    print("      Central = mid_x shape, R1998, EPPS21, Cloet Pzz=0.267.")
    print("      Plot B uses EPPS21 + R1998 (nuclear-PDF choice <1% effect).")


if __name__ == "__main__":
    main()
