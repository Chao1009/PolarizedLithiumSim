"""Unified Delta(x, Q2) model registry — double-helicity-flip structure
function of the spin-1 nucleus, one place for every model and constraint.

All Delta models used anywhere in the program live here behind ONE
callable interface (the (x, q2, f1) -> Delta signature that
`polligen.xsec.InclusiveKernel` and `fom.project_observables` consume),
so new constraints and models land in a single file and every consumer
switches by name:

    from polli_fastsim import delta_models
    model = delta_models.make("moment_A", f1_func=f1a_per_nucleon,
                              q2_ref=7.4, variant="mid_x")
    kernel = InclusiveKernel(beams.LI6, delta_func=model)
    print(model.info())          # plot-annotation string

Registered models
-----------------
toy       : legacy scale-based shape (polarized.toy_delta_gluon) --
            Delta = scale * F1 * x^0.3 (1-x)^4; `scale` = peak Delta/F1.
moment_A  : sum-rule-constrained ansatz, Interpretation A (ported from
            the money_delta suite, fastsim/scripts/money_delta_20260729.py):
                Delta = A * alpha_s(Q2) * F1(x,Q2) * shape(x),
                shape = x^a (1-x)^b / peak,
            with A solved so that the first moment obeys
                int_0^1 x Delta dx = c_moment * alpha_s  (alpha_s cancels;
            F1 evaluated at the rate-weighted <Q2> of the analysis).
            Bag-model coefficient c_moment = C_BAG = -0.012
            (Sather-Schmidt PRD 42:1424).
moment_B  : Interpretation B -- same, without the F1 factor:
                Delta = A * alpha_s(Q2) * shape(x),
            A from the analytic Beta-function integral.  B is the more
            conservative reading (their notes: factor 2-8 below A in
            reach); the A-vs-B spread IS the ansatz systematic.

Shape variants (money_delta convention): low_x (a,b) = (0.3, 4.0),
mid_x (0.7, 3.0) [production default], high_x (1.5, 2.0); shapes are
peak-normalized so `A * alpha_s` is the peak Delta/F1 (interp. A).

Convention note (plans/04 #6).  The money_delta suite folds the 2-of-6
per-nucleon dilution of 6Li into P_zz (their PZZ = 0.8/3 = 0.267).  Our
kernels keep P_zz as the whole-nucleus fill polarization, so the same
physics is expressed here with the `dilution` factor multiplying Delta
(default 1/3 for 6Li, the conservative Cloet convention already used in
beams.LI6): P_zz = 0.8 with dilution = 1/3 reproduces their numbers
exactly.  Do not set dilution != 1 AND use a per-nucleon-normalized
P_zz at the same time.  Known residual inconsistency (audited
2026-08-10): the scenario b1 curves (polarized.toy_b1) carry NO such
dilution -- acceptable while b1 is a shape scenario, but the b1 and
Delta conventions must be unified before any combined tensor fit.

Assumption made explicit (2026-08-10 audit): the Sather-Schmidt moment
was computed for a single bag (Delta++) target; applying it to the
PER-NUCLEON Delta of 6Li -- inherited verbatim from the money_delta
suite so the two lines stay comparable -- is itself a scenario choice,
on top of which `dilution` restricts the exotic glue to the polarized
pair.  <Q2>-convention note: the money_delta production solve used
Scenario(q2_min=2.0)+EPPS21 (<Q2> = 7.4, A = -0.310); our callers use
their own accepted phase space (q2_min=1.0, toy F1: <Q2> = 3.9,
A = -0.294) -- the F1 model, not q2_ref, dominates the small
difference.
"""

from dataclasses import dataclass, field

import numpy as np

from .polarized import toy_delta_gluon

# Sather-Schmidt bag-model sum-rule coefficient (PRD 42:1424):
# int_0^1 x Delta dx = C_BAG * alpha_s(Q2)
C_BAG = -0.012

# shape variants: {name: (alpha, beta)} -- money_delta convention
VARIANTS = {
    "low_x": (0.3, 4.0),
    "mid_x": (0.7, 3.0),
    "high_x": (1.5, 2.0),
}


def _xab_peak_value(alpha, beta):
    """Peak value of x^alpha (1-x)^beta at x = alpha/(alpha+beta)."""
    xp = alpha / (alpha + beta)
    return xp ** alpha * (1.0 - xp) ** beta


def shape_normalized(x, variant):
    """x^a (1-x)^b / peak for the named variant (peak value 1)."""
    a, b = VARIANTS[variant]
    x = np.asarray(x, dtype=float)
    return (np.power(np.maximum(x, 1e-12), a)
            * np.power(np.maximum(1.0 - x, 0.0), b)
            / _xab_peak_value(a, b))


def alpha_s_lo(q2, lambda_qcd=0.22, n_f=4):
    """LO analytic alpha_s = 12 pi / [(33 - 2 n_f) ln(Q2/Lambda2)] --
    the money_delta fallback (Lambda = 0.22 GeV, n_f = 4)."""
    q2 = np.asarray(q2, dtype=float)
    lam2 = lambda_qcd * lambda_qcd
    ln = np.log(np.maximum(q2, lam2 * 1.01) / lam2)
    return 12.0 * np.pi / ((33.0 - 2.0 * n_f) * ln)


_ALPHAS_CACHE = {}


def alpha_s_table(setname="CT18NLO"):
    """alpha_s(Q2) interpolated from the parton PDF table (AlphaS_Qs /
    AlphaS_Vals in the set's .info) -- the money_delta PRIMARY source.
    Returns a callable, or None if the parton package/grids are absent.
    """
    if setname in _ALPHAS_CACHE:
        return _ALPHAS_CACHE[setname]
    func = None
    try:
        from parton import mkPDF
        info = mkPDF(setname, 0).pdfset.info
        qs = np.asarray(info["AlphaS_Qs"], dtype=float)
        vals = np.asarray(info["AlphaS_Vals"], dtype=float)

        def func(q2, _qs=qs, _vals=vals):
            q = np.sqrt(np.maximum(np.asarray(q2, dtype=float),
                                   _qs[0] ** 2))
            return np.interp(q, _qs, _vals)
    except Exception:
        func = None
    _ALPHAS_CACHE[setname] = func
    return func


def default_alpha_s():
    """Preferred alpha_s: parton table (money_delta production convention,
    ~14% below the LO analytic at Q2 ~ 2-15 GeV2), LO fallback."""
    return alpha_s_table() or alpha_s_lo


@dataclass
class DeltaModel:
    """Callable Delta(x, q2, f1) with metadata for plot annotations."""
    name: str
    func: callable
    params: dict = field(default_factory=dict)

    def __call__(self, x, q2, f1):
        return self.func(x, q2, f1)

    def info(self):
        pieces = ", ".join("%s=%.4g" % (k, v) if isinstance(v, float)
                           else "%s=%s" % (k, v)
                           for k, v in self.params.items())
        return "%s(%s)" % (self.name, pieces)


def make_toy(scale=1e-3):
    """Legacy toy shape; `scale` is the peak Delta/F1."""
    return DeltaModel(
        "toy", lambda x, q2, f1: toy_delta_gluon(x, q2, f1, scale=scale),
        {"scale": scale})


def solve_A_interp_a(f1_func, q2_ref, variant="mid_x", c_moment=C_BAG,
                     x_grid=None):
    """A for Interpretation A: A * int x F1(x, q2_ref) shape(x) dx = c.

    `f1_func(x, q2)` must return the per-nucleon F1 (F1A/A).  The
    integral runs over the full x range (the low-x tail matters for the
    low_x variant) on a log+linear grid, trapezoid rule -- same
    numerics as money_delta_20260729.solve_A_from_sum_rule.
    """
    if x_grid is None:
        x_grid = np.unique(np.concatenate([
            np.logspace(-5, np.log10(0.5), 220),
            np.linspace(0.5, 0.9999, 120)]))
    f1 = np.asarray(f1_func(x_grid, np.full_like(x_grid, q2_ref)),
                    dtype=float)
    integrand = x_grid * f1 * shape_normalized(x_grid, variant)
    integral = float(np.trapz(integrand, x_grid))
    if abs(integral) < 1e-30:
        raise ValueError("sum-rule integral vanished (variant=%s)" % variant)
    return c_moment / integral


def solve_A_interp_b(variant="mid_x", c_moment=C_BAG):
    """A for Interpretation B (no F1): analytic Beta-function integral,
    int_0^1 x * x^a (1-x)^b / peak dx = B(a+2, b+1) / peak."""
    from math import gamma
    a, b = VARIANTS[variant]
    beta_int = gamma(a + 2.0) * gamma(b + 1.0) / gamma(a + b + 3.0)
    return c_moment * _xab_peak_value(a, b) / beta_int


def make_moment_a(f1_func, q2_ref, variant="mid_x", c_moment=C_BAG,
                  dilution=1.0, alphas=None):
    """Interpretation A: Delta = dilution * A * alpha_s(Q2) * F1 * shape.

    alphas=None resolves to default_alpha_s() (parton table preferred,
    matching the money_delta production convention; LO fallback)."""
    alphas = alphas or default_alpha_s()
    a_solved = solve_A_interp_a(f1_func, q2_ref, variant, c_moment)

    def func(x, q2, f1):
        return (dilution * a_solved * alphas(q2) * np.asarray(f1)
                * shape_normalized(x, variant))
    return DeltaModel("moment_A", func,
                      {"variant": variant, "c": c_moment, "A": a_solved,
                       "q2_ref": q2_ref, "dilution": dilution})


def make_moment_b(variant="mid_x", c_moment=C_BAG, dilution=1.0,
                  alphas=None):
    """Interpretation B: Delta = dilution * A * alpha_s(Q2) * shape.

    alphas=None resolves to default_alpha_s()."""
    alphas = alphas or default_alpha_s()
    a_solved = solve_A_interp_b(variant, c_moment)

    def func(x, q2, f1):
        return (dilution * a_solved * alphas(q2)
                * shape_normalized(x, variant))
    return DeltaModel("moment_B", func,
                      {"variant": variant, "c": c_moment, "A": a_solved,
                       "dilution": dilution})


_REGISTRY = {
    "toy": make_toy,
    "moment_A": make_moment_a,
    "moment_B": make_moment_b,
}


def available():
    return sorted(_REGISTRY)


def make(name, **params):
    """Build a registered Delta model by name.

    toy:      make("toy", scale=1e-3)
    moment_A: make("moment_A", f1_func=..., q2_ref=..., variant="mid_x",
                   dilution=1/3)
    moment_B: make("moment_B", variant="mid_x", dilution=1/3)
    """
    if name not in _REGISTRY:
        raise KeyError("unknown Delta model %r; available: %s"
                       % (name, available()))
    return _REGISTRY[name](**params)
