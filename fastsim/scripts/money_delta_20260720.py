#!/usr/bin/env python3
"""6Li gluonometry discovery reach: L_5σ vs Δ/F₁ — EPPS21 & nNNPDF3.0 comparison.

Differences from money_delta_realistic.py
──────────────────────────────────────────
1. Only the **mid config** (E_e = 10 GeV × p_ion = 50 GeV/u ⁶Li) is plotted.
   The top config is dropped.

2. P_zz is fixed at the Cloet convention: PZZ = 0.267 (= 0.8 / 3, per-nucleon
   normalized).  No P_zz scan.

3. Three nuclear F₂ sources are compared:
   - CT18NLO + smooth EMC ratio (existing "realistic" baseline from
     money_delta_realistic.py; shown as a gray dashed reference line only).
   - EPPS21nlo_CT18Anlo_Li6: nuclear PDF via `parton`, F₂ᴬ directly.
   - nNNPDF30_nlo_as_0118_A6_Z3: nuclear PDF via `parton`, F₂ᴬ directly.

4. Δ shapes and R conventions unchanged: 3 shapes (low_x, mid_x, high_x) × 2 R
   (R1998, Christy-Bosted) = 6 combinations per nuclear-PDF choice.

NuclearF2FromGrid
─────────────────
Drop-in replacement for NuclearF2 (from polli_fastsim.structure) when F₂ᴬ
comes directly from a nuclear PDF set.  Lives entirely in this script — the
polli_fastsim package is not modified.

Nuclear PDF normalization note: the `parton` package returns xf_A(x,Q²) in
units of *per nucleon* (confirmed by cross-check: EPPS21 Li6 / CT18NLO ≈ 0.97
at typical DIS kinematics, not × A = 6).  The downstream callers (fom.py
project_rates, sig2_per_fb_at_nuc) divide f2a(x,Q²) by ion.A to get the
per-nucleon cross section.  Therefore NuclearF2FromGrid.f2a() must return
A × F₂(per-nucleon from grid) to satisfy the contract.

Output
──────
fastsim/out/money_delta/money_delta_20260720_mid.png
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
# R2: Three Δ shape variants (copied verbatim from money_delta_realistic.py)
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


# ══════════════════════════════════════════════════════════════════════════════
# Core significance helpers
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

# x-axis: 15 log-spaced Δ/F₁ scale points
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
    ("EPPS21",   "EPPS21nlo_CT18Anlo_Li6"),
    ("nNNPDF3.0", "nNNPDF30_nlo_as_0118_A6_Z3"),
]


# ══════════════════════════════════════════════════════════════════════════════
# Plot builder
# ══════════════════════════════════════════════════════════════════════════════

def build_nuclear_pdf_plot(cfg, base, outdir):
    """Build the EPPS21 vs nNNPDF3.0 money plot for the mid config.

    Parameters
    ----------
    cfg    : BeamConfig
    base   : PartonF2  CT18NLO backend (for the CT18+EMC reference and toy)
    outdir : pathlib.Path

    Returns
    -------
    pathlib.Path : output PNG path
    dict         : summary numbers at Δ/F₁ = 1e-3
    """

    print(f"\n  Config: {cfg.label()}")

    # ── Step 1: pre-compute sig² at S0 for all (nuclear-set, variant, R) combos ──
    # Total: 2 nuclear sets × 3 shapes × 2 R = 12 evaluations.
    # Plus 1 CT18+EMC reference (mid_x, R1998) + 1 toy = 14 total.

    # combo_sig2[nuc_label][(variant, r_name)] = sig2 at S0
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

    # ── Step 2: CT18NLO + EMC reference (mid_x, R1998) ──────────────────────
    print(
        f"    [CT18+EMC ref] mid_x, R1998, CT18NLO+EMC, Pzz={PZZ:.3f} …",
        end=" ", flush=True,
    )
    with r_override(r1998):
        sig2_ref = sig2_per_fb_at_ct18emc(
            cfg, S0, PZZ, base, variant=CENTRAL_VARIANT,
        )
    print(f"L_5σ(1e-3) = {25.0 / max(sig2_ref, 1e-30):.2f} fb⁻¹/u")

    # ── Step 3: toy-backend comparison ──────────────────────────────────────
    print(
        f"    [toy] toy_delta_gluon, original R, P_zz={TOY_PZZ} …",
        end=" ", flush=True,
    )
    # Use original r (not overridden)
    sig2_toy = sig2_per_fb_toy(cfg, S0, TOY_PZZ, base)
    print(f"L_5σ(1e-3) = {25.0 / max(sig2_toy, 1e-30):.2f} fb⁻¹/u")

    # ── Step 4: build reach curves ───────────────────────────────────────────
    # EPPS21 band
    epps_curves = []
    for (variant, r_name), s2 in combo_sig2["EPPS21"].items():
        epps_curves.append(25.0 / np.maximum(s2 * (SCALES / S0) ** 2, 1e-30))
    epps_curves = np.array(epps_curves)
    epps_min = epps_curves.min(axis=0)
    epps_max = epps_curves.max(axis=0)

    # nNNPDF3.0 band
    nnpdf_curves = []
    for (variant, r_name), s2 in combo_sig2["nNNPDF3.0"].items():
        nnpdf_curves.append(25.0 / np.maximum(s2 * (SCALES / S0) ** 2, 1e-30))
    nnpdf_curves = np.array(nnpdf_curves)
    nnpdf_min = nnpdf_curves.min(axis=0)
    nnpdf_max = nnpdf_curves.max(axis=0)

    # Central curve: mid_x + R1998 + EPPS21
    sig2_central = combo_sig2["EPPS21"][(CENTRAL_VARIANT, CENTRAL_R_NAME)]
    reach_central = 25.0 / np.maximum(
        sig2_central * (SCALES / S0) ** 2, 1e-30
    )

    # CT18+EMC reference and toy
    reach_ref = 25.0 / np.maximum(sig2_ref * (SCALES / S0) ** 2, 1e-30)
    reach_toy = 25.0 / np.maximum(sig2_toy * (SCALES / S0) ** 2, 1e-30)

    # ── Step 5: draw the plot ────────────────────────────────────────────────
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
        SCALES, reach_central, "-", color="black", lw=2.2,
        label=r"Central: mid-$x$ $\Delta$, R1998, EPPS21, Cloet $P_{zz}=0.267$",
    )

    # CT18+EMC reference (gray dashed)
    ax.plot(
        SCALES, reach_ref, "--", color="gray", lw=1.6,
        label=r"CT18NLO + EMC reference (mid-$x$, R1998)",
    )

    # Toy-backend comparison (gray dotted)
    ax.plot(
        SCALES, reach_toy, ":", color="gray", lw=1.3,
        label=r"Toy backend, $P_{zz}=0.8$ (continuity ref.)",
    )

    # Gold band: plausible EIC luminosity band
    ax.axhspan(
        1, 100, color="gold", alpha=0.12,
        label=r"1$-$100 fb$^{-1}$/u (plausible program)",
    )

    # Sather-Schmidt reference line
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
    ax.set_ylabel(r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]", fontsize=11)
    ax.set_title(
        r"$^6$Li gluonometry realistic reach — 10 GeV $e$ × 50 GeV/u $^6$Li"
        "\n"
        r"EPPS21 vs nNNPDF3.0 nuclear PDFs, Cloet $P_{zz}=0.267$",
        fontsize=10,
    )
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "money_delta_20260720_mid.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"  wrote {outpath}")

    # ── Step 6: collect summary numbers at Δ/F₁ = 1e-3 ─────────────────────
    l5_central = 25.0 / max(sig2_central, 1e-30)

    # EPPS21 band extremes
    epps_s2_vals = list(combo_sig2["EPPS21"].values())
    epps_s2_max  = max(epps_s2_vals)   # highest sig² → lowest L_5σ
    epps_s2_min  = min(epps_s2_vals)   # lowest sig²  → highest L_5σ
    l5_epps_lo   = 25.0 / max(epps_s2_max, 1e-30)
    l5_epps_hi   = 25.0 / max(epps_s2_min, 1e-30)

    # nNNPDF3.0 band extremes
    nnpdf_s2_vals = list(combo_sig2["nNNPDF3.0"].values())
    nnpdf_s2_max  = max(nnpdf_s2_vals)
    nnpdf_s2_min  = min(nnpdf_s2_vals)
    l5_nnpdf_lo   = 25.0 / max(nnpdf_s2_max, 1e-30)
    l5_nnpdf_hi   = 25.0 / max(nnpdf_s2_min, 1e-30)

    l5_ref = 25.0 / max(sig2_ref, 1e-30)
    l5_toy = 25.0 / max(sig2_toy, 1e-30)

    summary = {
        "label":        cfg.label(),
        "l5_central":   l5_central,
        "l5_epps_lo":   l5_epps_lo,
        "l5_epps_hi":   l5_epps_hi,
        "l5_nnpdf_lo":  l5_nnpdf_lo,
        "l5_nnpdf_hi":  l5_nnpdf_hi,
        "l5_ref":       l5_ref,
        "l5_toy":       l5_toy,
    }
    return outpath, summary


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description=(
            "6Li gluonometry reach: L_5σ vs Δ/F₁, comparing EPPS21 and "
            "nNNPDF3.0 nuclear PDFs at the mid EIC config (10 GeV × 50 GeV/u)."
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

    # ── Mid config only: E_e = 10 GeV, p_ion = 50 GeV/u ────────────────────
    cfg = BeamConfig(
        electron_energy=10.0,
        ion=LI6,
        ion_momentum_per_nucleon=50.0,
    )

    print(f"\n{'='*70}")
    print("Building nuclear PDF comparison plot: MID config — e(10) × 6Li(50/u)")
    print(f"{'='*70}")

    _, summary = build_nuclear_pdf_plot(cfg, base, outdir)

    # ── Summary table ─────────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("SUMMARY  (Δ/F₁ = 1e-3, all values in fb⁻¹/nucleon)")
    print(f"Config: {summary['label']}")
    print("=" * 72)
    print(
        f"  Central (mid-x, R1998, EPPS21, Cloet Pzz=0.267):  "
        f"{summary['l5_central']:>8.2f} fb⁻¹/u"
    )
    print(
        f"  EPPS21 band:  "
        f"lo = {summary['l5_epps_lo']:>8.2f}  "
        f"hi = {summary['l5_epps_hi']:>8.2f}  fb⁻¹/u  "
        f"(lo = most sensitive, hi = least)"
    )
    print(
        f"  nNNPDF3.0 band:  "
        f"lo = {summary['l5_nnpdf_lo']:>8.2f}  "
        f"hi = {summary['l5_nnpdf_hi']:>8.2f}  fb⁻¹/u"
    )
    print(
        f"  CT18NLO + EMC reference (mid-x, R1998):            "
        f"{summary['l5_ref']:>8.2f} fb⁻¹/u"
    )
    print(
        f"  Toy backend (toy_delta_gluon, original R, Pzz=0.8):"
        f"{summary['l5_toy']:>8.2f} fb⁻¹/u"
    )
    print("=" * 72)
    print("Note: Band lo/hi = lowest/highest L_5σ over 6 combinations")
    print("      (3 Δ shapes × 2 R parameterizations), Cloet Pzz=0.267.")
    print("      Central = mid_x shape, R1998, EPPS21.")


if __name__ == "__main__":
    main()
