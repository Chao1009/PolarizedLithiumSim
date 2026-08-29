"""Asymmetry formulas for spin-1/2 and spin-1 (tensor) targets.

Spin-1 master formula (unpolarized e, target spin at angle theta_m, target
spin projection lambda_m in {+1, 0, -1}; Hoodbhoy-Jaffe-Manohar NPB 312:571,
as in J. Maxwell's slides, docs/Discussions.pptx p.5):

  dsigma/(dx dy dphi)(lambda_m) = (2 y alpha^2 / Q^2) *
      [ F1 + (2/3) a_m b1 + (1-y)/(x y^2) (F2 + (2/3) a_m b2)
        - (1-y)/y^2 * c_m * sin^2(theta_m) * Delta(x,Q2) * cos(2 phi) ]

  a_m = (1/4) c_m (3 cos^2 theta_m - 1),  c_m = 3|lambda_m| - 2 -> (1,-2,1)

Conventions chosen here (document in any plot):
  Azz  = (s+ + s- - 2 s0) / (s+ + s- + s0)        (longitudinal, theta_m=0)
       = (2/3) [b1 + (1-y)/(x y^2) b2] / [F1 + (1-y)/(x y^2) F2] -> measures b1
  A2phi(lambda=+-1, theta_m=90deg)
       = -(1-y)/y^2 * Delta / D_phi  with D_phi = F1 + (1-y)/(x y^2) F2
"""

import numpy as np

from .structure import r_sigma_lt

# Sign of the tensor RATE (b1, b2) sector, and the single place it is
# set for the whole program (polligen.xsec imports this constant).
#
#   +1  the program's transcription of Hoodbhoy-Jaffe-Manohar
#       (docs/Discussions.pptx p.5), giving Azz = +(2/3) b1/F1
#   -1  the HJM/HERMES convention as written by Cosyn, Roldan Tomei,
#       Sosa and Zec, EPJ A 61 (2025) 83 (arXiv:2410.12764) Eq. (27),
#       Azz = -(2/3) b1/F1
#
# The two differ by the sign of b1 itself, so |Azz| and the whole Delta
# (cos 2phi) sector are unaffected; what flips is the sign of the
# by-product kappa of the spin-state ratio, and with it the sign of any
# O(gamma^2) b-sector subtraction built on kappa.  Settling it needs
# Jaffe-Manohar PLB 223 (1989) 218, which is not in refs/ (code review
# G1, recommendation 0a; plans/08 D1).  Changing this constant is the
# whole of that decision: tests/test_tensor_convention.py pins the
# identity it controls.
TENSOR_LL_SIGN = +1.0


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


def a_parallel(g1, f1, y, x, q2, r_func=None):
    """Longitudinal double-spin asymmetry, A1 ~= g1/F1 approximation.

    Massless kinematics: this is the gamma -> 0 limit of the exact
    A_par = D_gamma (A1 + eta A2) with gamma^2 = 4 M^2 x^2/Q^2,
    A1 = (g1 - gamma^2 g2)/F1, A2 = gamma (g1 + g2)/F1 and
    eta = eps gamma y/[1 - (1-y) eps] (E143 PRD 58:112003).  Everything
    published is this function; the exact form lives behind
    `polligen.xsec.InclusiveKernel(..., target_mass=True)`, default off.

    What the limit drops, measured with g2 = g2_WW by
    `evgen/scripts/target_mass_bound.py` (toy backends, the published
    ones).  At the small y of the analysis the whole correction
    collapses to a multiplicative (1 + gamma^2) on A_par, independent of
    g2, and the W^2 >= 10 GeV^2 cut of `fom.Scenario` caps it everywhere
    at gamma^2 <= M^2/(W2_min - M^2) = 0.0965, with grid maxima 0.0854 /
    0.0577 / 0.0258 for 6Li at 5 x 40.8, 10 x 99.5 and 18 x 137.5 GeV/u
    (0.0854 / 0.0577 / 0.0332 for 7Li).  On the two observables built
    from A_par it is much smaller than that, and it has a sign: the
    exact A_par is (1 + gamma^2) times this one, so inverting it with
    the massless D leaves g1/F1 -- and the polarized-EMC Delta-R built
    from it -- HIGH, by 0.12 / 0.44 / 0.71 / 1.06 % at
    x = 0.089 / 0.282 / 0.447 / 0.708
    (inverse-variance weighted over Q2 and the three energies, i.e. with
    the weights of the published error bars), and 2.1 % at most on the
    sigma-weighted tagged-triton A_par overlay of the published 7Li
    configuration (5.0 % at 5 x 40.8, both in the top x bin).  The eta A2
    piece alone is 0.08-0.17 % and 0.56 % there: it is HALF the
    correction, and above x ~ 0.5 the half that vanishes, because
    g2_WW -> -g1 kills A2 while the -gamma^2 g2/F1 inside A1 survives.
    M is the free-nucleon mass; the bound-nucleon one moves gamma^2 by
    1.0 %.
    """
    d = depolarization_d(y, x, q2, r_func=r_func)
    return d * g1 / np.maximum(f1, 1e-30)


def phi_averaged_density(f1, f2, x, y):
    """The phi-averaged bracket D_phi = F1 + (1-y)/(x y^2) F2."""
    return f1 + (1.0 - y) / (x * y * y) * f2


def azz(b1, f1, f2, x, y, b2=None, theta_m=0.0):
    """Tensor asymmetry from the master formula. b2 defaults to 2x*b1."""
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
