"""Asymmetry formulas for spin-1/2 and spin-1 (tensor) targets.

Spin-1 master formula (unpolarized e, target spin at angle theta_m, target
spin projection lambda_m in {+1, 0, -1}; Hoodbhoy-Jaffe-Manohar NPB 312:571,
as in J. Maxwell's slides, docs/Discussions.pptx p.5):

  dsigma/(dx dy dphi)(lambda_m) = (2 y alpha^2 / Q^2) *
      [ F1 + (2/3) a_m b1 + (1-y)/(x y^2) (F2 + (2/3) a_m b2)
        - (1-y)/y^2 * c_m * sin^2(theta_m) * Delta(x,Q2) * cos(2 phi) ]

  a_m = (1/4) c_m (3 cos^2 theta_m - 1),  c_m = 3|lambda_m| - 2 -> (1,-2,1)

The b1, b2 terms carry the overall sign TENSOR_LL_SIGN below, which
since 2026-08-29 is the LITERATURE one (Cosyn et al. Eq. 27, HERMES) and
not the transcription above: the master formula as written has the
opposite sign of b1 to the b1 every published number is quoted in.

Conventions chosen here (document in any plot):
  Azz  = (s+ + s- - 2 s0) / (s+ + s- + s0)        (longitudinal, theta_m=0)
       = -(2/3) [b1 + (1-y)/(x y^2) b2] / [F1 + (1-y)/(x y^2) F2]
       = -(2/3) (b1/F1) / (1 + eps(y) R)  exactly -> measures b1
  A2phi(lambda=+-1, theta_m=90deg)
       = -(1-y)/y^2 * Delta / D_phi  with D_phi = F1 + (1-y)/(x y^2) F2
"""

import numpy as np

from .structure import r_sigma_lt

# Per-nucleon target mass for gamma = 2 M x / Q.  It is the FREE-nucleon
# mass because x is per-nucleon; the bound-nucleon mass
# (`beams.Ion.mass_per_nucleon`, 0.9336 GeV for 6Li) would move gamma^2 by
# 1.0%, i.e. 1% of a <= 10% correction.
M_NUCLEON = 0.9383

# Sign of the tensor RATE (b1, b2) sector, and the single place it is
# set for the whole program (polligen.xsec imports this constant).
#
#   -1  THE CONVENTION OF THIS PROGRAM SINCE 2026-08-29 (author decision,
#       plans/08 D1): the literature one, as written by Cosyn, Roldan
#       Tomei, Sosa and Zec, EPJ A 61 (2025) 83 (arXiv:2410.12764)
#       Eq. (27) and used by HERMES,
#
#           Azz = -(2/3) b1/F1     (axis along q, Bjorken limit)
#
#       with P_zz = n+ + n- - 2 n0 the tensor polarization, i.e. b1 > 0
#       means the m = 0 state is the one with the LARGER cross section.
#   +1  the repository's own private convention until that date -- the
#       transcription of Hoodbhoy-Jaffe-Manohar in docs/Discussions.pptx
#       p.5, giving Azz = +(2/3) b1/F1.  Kept reachable by setting this
#       constant back, which is the whole of the change: nothing else in
#       the program knows the sign.
#
# The decision was taken on the literature and not on a new derivation:
# Jaffe-Manohar PLB 223 (1989) 218 is still not in refs/, so what is
# adopted is the convention every published b1 number is quoted in (the
# HERMES extraction and the Cosyn et al. re-analysis of it), which is the
# only way an extracted b1 of this programme can be compared with them.
# The two differ by the sign of b1 itself, so |Azz| and the whole Delta
# (cos 2phi) sector are unaffected; what flips is the sign of Azz at
# fixed b1, of the by-product kappa of the spin-state ratio, and of any
# b-sector subtraction built on kappa -- including the O(gamma^2) tensor
# leakage into cos 2phi, which is why plans/08 D2 was gated on this.
# tests/test_tensor_convention.py pins the identity the constant
# controls, in Cosyn's own form.
TENSOR_LL_SIGN = -1.0


def depolarization_d(y, x, q2, r_func=None):
    """Virtual-photon depolarization factor D for A_par ~= D * A1.

    `r_func(x, q2) -> R = sigma_L/sigma_T`; None (the default) uses THIS
    module's `r_sigma_lt` global, resolved at call time.  The lookup has
    to stay in this module: the dated money scripts rebind
    `asymmetries.r_sigma_lt` and `structure.r_sigma_lt` separately
    (fastsim/scripts/money_delta_*.py `r_override`), and they are frozen
    reproductions of dated notes.  With r_func=None the value is
    bit-for-bit what it was before the hook existed.
    """
    r = (r_sigma_lt if r_func is None else r_func)(x, q2)
    return y * (2.0 - y) / (y * y + 2.0 * (1.0 - y) * (1.0 + r))


# --- finite-gamma (target-mass) kinematics, E143 PRD 58:112003 -----------
#
# ONE implementation, used by both packages: `polligen.xsec` imports these
# three functions rather than carrying its own copies, so the generator
# kernel and the fast simulation cannot drift apart in the O(gamma^2) term
# (they did not, but only because the two transcriptions were written on
# the same afternoon).  They are the exact lab-frame factors written in
# (x, y); `evgen/tests/test_target_mass.py` pins them to double precision
# against E143's own (E, E', theta) definitions
#   eps = 1/[1 + 2(1 + nu^2/Q^2) tan^2(theta/2)]
#   D   = (1 - E' eps/E)/(1 + eps R)
#   eta = eps sqrt(Q^2)/(E - E' eps).
# At gamma -> 0 they collapse to the massless set: eps -> (1-y)/(1-y+y^2/2)
# and D -> `depolarization_d` above.


def gamma_squared(x, q2, m=M_NUCLEON):
    """Target-mass parameter gamma^2 = 4 M^2 x^2 / Q^2 (= Q^2/nu^2)."""
    return 4.0 * m * m * np.asarray(x, dtype=float) ** 2 / np.asarray(
        q2, dtype=float)


def epsilon_gamma(y, gamma2):
    """Virtual-photon transverse polarization at finite gamma."""
    y = np.asarray(y, dtype=float)
    return ((1.0 - y - 0.25 * gamma2 * y * y)
            / (1.0 - y + 0.5 * y * y + 0.25 * gamma2 * y * y))


def depolarization_d_gamma(y, gamma2, r):
    """D_gamma = [1 - (1-y) eps]/(1 + eps R) with the finite-gamma eps."""
    eps = epsilon_gamma(y, gamma2)
    return (1.0 - (1.0 - np.asarray(y, dtype=float)) * eps) / (1.0 + eps * r)


def eta_gamma(y, gamma2):
    """eta = eps gamma y/[1 - (1-y) eps], the A2 admixture in A_par."""
    y = np.asarray(y, dtype=float)
    eps = epsilon_gamma(y, gamma2)
    return eps * np.sqrt(gamma2) * y / (1.0 - (1.0 - y) * eps)


def a_parallel_exact(g1, g2, f1, y, x, q2, r_func=None):
    """Longitudinal double-spin asymmetry at finite gamma (E143).

        A_par = D_gamma (A1 + eta A2),
        A1 = (g1 - gamma^2 g2)/F1,   A2 = gamma (g1 + g2)/F1,

    with gamma^2 = 4 M^2 x^2/Q^2 and the eps, D_gamma, eta above.  This is
    the DEFAULT longitudinal kernel of the programme since 2026-08-29
    (author decision 1 of run 15): it is exact given g2, it costs one
    g2^WW table per grid and nothing per call, and switching it on removes
    the O(gamma^2) bias the massless form left on every extracted g1/F1.

    Both O(gamma^2) pieces are kept.  eta A2 alone is about half the
    correction, and above x ~ 0.5 it is the half that vanishes, because
    g2^WW -> -g1 there kills A2 while the -gamma^2 g2/F1 inside A1
    survives.  At the small y of this programme the whole correction
    collapses to a multiplicative (1 + gamma^2) on A_par, independent of
    g2, and the W^2 >= 10 GeV^2 cut of `fom.Scenario` caps gamma^2
    everywhere at M^2/(W2_min - M^2) = 0.0965.

    What is left over is the g2 model, not the kinematics: g2 is taken
    Wandzura-Wilczek (`polarized.g2_ww`), and the twist-3 departure from
    WW is the residual systematic that replaced the bias.  It is measured
    by `evgen/scripts/target_mass_bound.py`, which repeats the extraction
    with g2 = 0 and g2 = 1.5 g2^WW.
    """
    r = (r_sigma_lt if r_func is None else r_func)(x, q2)
    g2v = gamma_squared(x, q2)
    f1 = np.maximum(f1, 1e-30)
    a1 = (g1 - g2v * g2) / f1
    a2 = np.sqrt(g2v) * (g1 + g2) / f1
    return depolarization_d_gamma(y, g2v, r) * (a1 + eta_gamma(y, g2v) * a2)


def a_parallel(g1, f1, y, x, q2, r_func=None, g2=None):
    """Longitudinal double-spin asymmetry A_par(x, y).

    With `g2` supplied this is `a_parallel_exact` -- the finite-gamma form
    the programme now uses by default.  With `g2=None` it is the massless
    limit A_par = D(y) g1/F1, bit-for-bit what it was before the term
    existed, which is how every figure published before 2026-08-29 was
    made and how `InclusiveKernel(..., target_mass=False)` still computes.

    The limit is not small enough to ignore, which is why it is no longer
    the default: the exact A_par is (1 + gamma^2) times this one at these
    y, so inverting it with the massless D left the extracted g1/F1 -- and
    the polarized-EMC Delta-R built from it -- HIGH by 0.1-1% across the
    published x range.  `depolarization_effective` is the divisor that
    removes exactly that.
    """
    if g2 is not None:
        return a_parallel_exact(g1, g2, f1, y, x, q2, r_func=r_func)
    d = depolarization_d(y, x, q2, r_func=r_func)
    return d * g1 / np.maximum(f1, 1e-30)


def depolarization_effective(y, x, q2, g2_over_g1=None, r_func=None):
    """The factor D_eff with A_par = D_eff * (g1/F1), and so the divisor
    that turns delta(A_par) into delta(g1/F1).

    Writing rho = g2/g1 and gamma^2 = 4 M^2 x^2/Q^2,

        A_par = D_gamma [ (g1 - gamma^2 g2)/F1 + eta gamma (g1 + g2)/F1 ]
              = D_gamma [ 1 - gamma^2 rho + eta gamma (1 + rho) ] (g1/F1)
              = D_eff (g1/F1),

    so an extraction that divides a measured A_par by D_eff returns g1/F1
    with no O(gamma^2) bias, where dividing by the massless D returned it
    high by (1 + gamma^2) + O(gamma^2 y).  The step from the first line to
    the second is legitimate BECAUSE rho is a property of the shape of g1
    and not of its normalization: g2^WW is linear in g1
    (`polarized.g2_ww`), and the polarized-EMC observable Delta-R is a
    multiplicative modification of g1 at fixed shape, so rho is the same
    for the model and for the measurement.  What rho is NOT known to be is
    the truth -- that is the twist-3 residual, quantified by varying rho.

    `g2_over_g1=None` returns `depolarization_d` unchanged, bit-for-bit,
    for the massless path.
    """
    if g2_over_g1 is None:
        return depolarization_d(y, x, q2, r_func=r_func)
    r = (r_sigma_lt if r_func is None else r_func)(x, q2)
    g2v = gamma_squared(x, q2)
    rho = np.asarray(g2_over_g1, dtype=float)
    bracket = (1.0 - g2v * rho
               + eta_gamma(y, g2v) * np.sqrt(g2v) * (1.0 + rho))
    return depolarization_d_gamma(y, g2v, r) * bracket


def phi_averaged_density(f1, f2, x, y):
    """The phi-averaged bracket D_phi = F1 + (1-y)/(x y^2) F2."""
    return f1 + (1.0 - y) / (x * y * y) * f2


def azz(b1, f1, f2, x, y, b2=None, theta_m=0.0):
    """Tensor asymmetry from the master formula. b2 defaults to 2x*b1.

    Carries TENSOR_LL_SIGN, so with the convention adopted 2026-08-29
    A_zz(theta_m = 0) (1 + eps R) = -(2/3) b1/F1 exactly, at every y
    (Cosyn et al. Eq. 27; pinned in evgen/tests/test_tensor_convention.py).
    """
    if b2 is None:
        b2 = 2.0 * x * b1
    geom = 0.5 * (3.0 * np.cos(theta_m) ** 2 - 1.0)  # =1 at theta_m=0
    num = (TENSOR_LL_SIGN * 2.0 / 3.0
           * (b1 + (1.0 - y) / (x * y * y) * b2) * geom)
    return num / phi_averaged_density(f1, f2, x, y)


def a_cos2phi(delta, f1, f2, x, y):
    """cos(2phi) amplitude for lambda_m = +-1, theta_m = 90 deg."""
    return -(1.0 - y) / (y * y) * delta / phi_averaged_density(f1, f2, x, y)


# --- statistical uncertainties (per kinematic bin with N events total) ---


def err_a_parallel(n, pe, pz):
    """delta(A_par): two-state +/- flips, equal luminosity halves."""
    n = np.maximum(np.asarray(n, dtype=float), 1e-12)
    return 1.0 / (pe * pz * np.sqrt(n))


def err_azz(n, pzz):
    """delta(Azz) for the (n+ + n- - 2 n0)/(n+ + n- + n0) estimator,
    equal luminosity thirds: Var(num) ~= 2N -> delta = sqrt(2/N) / Pzz."""
    n = np.maximum(np.asarray(n, dtype=float), 1e-12)
    return np.sqrt(2.0 / n) / pzz


def err_cos2phi_amplitude(n, pzz):
    """delta(amplitude) of a cos(2phi) fit: sqrt(2/N), scaled by tensor
    polarization of the transverse spin states."""
    n = np.maximum(np.asarray(n, dtype=float), 1e-12)
    return np.sqrt(2.0 / n) / pzz
