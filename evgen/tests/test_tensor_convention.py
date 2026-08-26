"""External-convention tests for the tensor sector (code review G1).

Every other tensor assertion in the suite compares code to code -- the
kernel against `asymmetries.azz`, which shares the same transcription --
so patching the sign in both files left the whole suite green.  These
tests compare against a LITERATURE relation instead:

    Cosyn, Roldan Tomei, Sosa, Zec, EPJ A 61 (2025) 83
    (arXiv:2410.12764, refs/2410.12764v1.pdf), Eq. (27):

        A_zz = -(2/3) b1 / F1

The program's transcription of Hoodbhoy-Jaffe-Manohar gives the opposite
sign, carried by the single constant `asymmetries.TENSOR_LL_SIGN`.  The
tests below are written against that constant, so they are green today
and become the switch when the author settles the convention
(plans/08 D1); they also pin the exact kinematic factor that relates the
code's low-y form to the literature's simple ratio.
"""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen.xsec import EventSpinState, InclusiveKernel  # noqa: E402

from polli_fastsim import asymmetries as asym  # noqa: E402
from polli_fastsim import beams  # noqa: E402
from polli_fastsim.polarized import toy_b1  # noqa: E402
from polli_fastsim.structure import r_sigma_lt  # noqa: E402

CONFIG = beams.default_configs("6Li")[1]
S = CONFIG.sqrt_s_per_nucleon ** 2

# (x, Q2) points across the analysis window
POINTS = [(0.0224, 1.14), (0.056, 1.14), (0.141, 3.13), (0.141, 14.3),
          (0.005, 2.0), (0.30, 30.0)]


def _eps(y):
    """Longitudinal photon polarization epsilon(y) = (1-y)/(1-y+y^2/2)."""
    return (1.0 - y) / (1.0 - y + 0.5 * y * y)


def test_azz_matches_cosyn_eq27_up_to_the_program_sign():
    """A_zz(theta_S = 0) (1 + eps(y) R) == TENSOR_LL_SIGN * (2/3) b1/F1,
    EXACTLY and at every y.

    Derivation: with b2 = 2x b1 and F2 = 2x(1+R) F1 (the repository's own
    Callan-Gross-with-R relation), the (1-y)/(x y^2) terms of numerator
    and denominator combine into 2(1-y+y^2/2) and 2[(1-y+y^2/2)+R(1-y)],
    whose ratio is exactly 1/(1 + eps R).  Note that the widely-quoted
    form which ALSO divides by [1 + 2(1-y)/y^2] and multiplies by
    [1 + 2(1-y)(1+R)/y^2] double-counts R -- those two brackets are
    identically 1 + eps R -- and misses by a factor 1.17.
    """
    for x, q2 in POINTS:
        f1 = 1.0                      # the identity is F1-free
        b1 = 0.037                    # any value; only the ratio enters
        r = float(r_sigma_lt(x, q2))
        f2 = 2.0 * x * (1.0 + r) * f1
        for y in (0.01, 0.05, 0.2, 0.6, 0.9):
            azz = float(asym.azz(b1, f1, f2, x, y))
            assert azz * (1.0 + _eps(y) * r) == pytest.approx(
                asym.TENSOR_LL_SIGN * (2.0 / 3.0) * b1 / f1, rel=1e-12)


def test_the_program_sign_is_opposite_to_the_literature():
    """Stated as a test so that flipping TENSOR_LL_SIGN is a deliberate
    act with a visible consequence, not a silent one."""
    x, q2, y, b1 = 0.056, 1.14, 0.05, 0.037
    r = float(r_sigma_lt(x, q2))
    f2 = 2.0 * x * (1.0 + r)
    cosyn = -(2.0 / 3.0) * b1                     # Cosyn Eq. (27), F1 = 1
    program = float(asym.azz(b1, 1.0, f2, x, y)) * (1.0 + _eps(y) * r)
    assert asym.TENSOR_LL_SIGN == +1.0, (
        "TENSOR_LL_SIGN was changed: this is plans/08 D1.  Update the "
        "sign of kappa in the reconstruction-chain report and any "
        "kappa-based subtraction at the same time.")
    assert program == pytest.approx(-cosyn, rel=1e-12)


def test_kernel_thirds_combination_reproduces_azz_with_its_sign():
    """The kernel's per-state rate shifts must build the SAME A_zz, sign
    included, so that the one constant governs both.  A_zz is the thirds
    combination (s+ + s- - 2 s0)/(s+ + s- + s0), not a single state: the
    m = +-1 shift alone is half of it."""
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    x, q2 = 0.056, 1.14
    y = q2 / (S * x)
    xa, qa = np.array([x]), np.array([q2])
    t = kern.tables(xa, qa)
    w = {}
    for m in (1.0, 0.0, -1.0):
        state = EventSpinState(lam_e=0, pe=0.0, j=1.0, m=m, theta_s=0.0)
        w[m], _, _ = kern.amplitudes(t, xa, qa, S, state)
    num = w[1.0] + w[-1.0] - 2.0 * w[0.0]
    den = 3.0 + w[1.0] + w[-1.0] + w[0.0]
    expected = asym.azz(t["b1"], t["f1"], t["f2"], x, y, b2=t["b2"])
    np.testing.assert_allclose(num / den, expected, rtol=1e-12)
    np.testing.assert_allclose(w[1.0], 0.5 * expected, rtol=1e-12)
    assert np.sign(float(w[1.0])) == np.sign(asym.TENSOR_LL_SIGN
                                             * float(t["b1"]))


def test_delta_sector_does_not_depend_on_the_convention():
    """|A_zz| and the whole cos 2phi sector are convention-independent:
    only the sign of the b1 rate term (and hence of kappa) flips."""
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1,
                           delta_func=lambda x, q2, f1: -1e-2 * f1)
    x, q2 = 0.056, 1.14
    y = q2 / (S * x)
    t = kern.tables(np.array([x]), np.array([q2]))
    state = EventSpinState(lam_e=0, pe=0.0, j=1.0, m=1.0,
                           theta_s=np.pi / 2.0)
    _, _, a2 = kern.amplitudes(t, np.array([x]), np.array([q2]), S, state)
    expected = asym.a_cos2phi(t["delta"], t["f1"], t["f2"], x, y)
    np.testing.assert_allclose(a2, expected, rtol=1e-12)
    assert float(a2) > 0.0            # Delta < 0 with c_m = +1 -> a2 > 0
    assert "TENSOR_LL_SIGN" not in asym.a_cos2phi.__doc__.upper()


# --- one rank-2 geometry for both spins (plans/08 B3) ----------------------

def _q_nn(j, m):
    return (3.0 * m * m - j * (j + 1.0)) / 3.0


def test_rank2_moments_are_one_formula_for_both_spins():
    """Cosyn Eq. (9): t_ij = (Q_NN/2)(3 n_i n_j - d_ij) for any J, so
    T_LL = Q_NN P_2 and T_TT = (3/2) Q_NN sin^2.  The kernel returns
    (Q_NN, 3 Q_NN) and builds both channels from it."""
    for ion, j, ms in ((beams.LI6, 1.0, (1.0, 0.0, -1.0)),
                       (beams.LI7, 1.5, (1.5, 0.5, -0.5, -1.5))):
        kern = InclusiveKernel(ion)
        for m in ms:
            q, c = kern._tensor_moments(m)
            assert q == pytest.approx(_q_nn(j, m), rel=1e-12)
            assert c == pytest.approx(3.0 * q, rel=1e-12)
    # spin 1/2 carries no rank-2 alignment
    assert InclusiveKernel(beams.PROTON)._tensor_moments(0.5) == (0.0, 0.0)


def test_spin1_geometry_is_unchanged_by_the_unification():
    """The J = 1 branch must reproduce the Hoodbhoy-Jaffe-Manohar
    transcription digit for digit: Q_NN = c_m/3 and 3 Q_NN = c_m."""
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    for m in (1.0, 0.0, -1.0):
        c_m = 3.0 * m * m - 2.0
        q, c = kern._tensor_moments(m)
        assert q == pytest.approx(c_m / 3.0, rel=1e-13)
        assert c == pytest.approx(c_m, rel=1e-13)
    x, q2 = 0.056, 1.14
    xa, qa = np.array([x]), np.array([q2])
    t = kern.tables(xa, qa)
    state = EventSpinState(lam_e=0, pe=0.0, j=1.0, m=1.0, theta_s=0.0)
    w_avg, _, _ = kern.amplitudes(t, xa, qa, S, state)
    expected = 0.5 * asym.azz(t["b1"], t["f1"], t["f2"], x,
                              q2 / (S * x), b2=t["b2"])
    np.testing.assert_allclose(w_avg, expected, rtol=1e-13)


def test_spin32_rate_and_cos2phi_channels_are_now_consistent():
    """CHARACTERIZATION, not a physics assertion: the spin-3/2 rank-2
    normalization is plans/04 #14 and this records the adopted one.  What
    the test does pin is INTERNAL consistency -- before 2026-08-25 the
    J = 3/2 branch returned (Q_NN, Q_NN), so its rate and cos 2phi
    channels disagreed with each other by 3 while J = 1 did not."""
    def ratio(ion, j):
        kern = InclusiveKernel(
            ion, b1_32_func=lambda x, q2, f1: 0.05 * f1,
            delta_32_func=lambda x, q2, f1: -1e-2 * f1,
            b1_func=lambda x, q2, f1: 0.05 * f1,
            delta_func=lambda x, q2, f1: -1e-2 * f1)
        x, q2 = 0.056, 1.14
        xa, qa = np.array([x]), np.array([q2])
        t = kern.tables(xa, qa)
        m = j                                     # the stretched state
        long_ = kern.amplitudes(t, xa, qa, S,
                                EventSpinState(0, 0.0, j, m, theta_s=0.0))[0]
        trans = kern.amplitudes(t, xa, qa, S,
                                EventSpinState(0, 0.0, j, m,
                                               theta_s=np.pi / 2))[2]
        return float(trans) / float(long_)

    r1, r32 = ratio(beams.LI6, 1.0), ratio(beams.LI7, 1.5)
    # the ratio is Q_NN-free by construction, so it must be spin-independent
    assert r32 == pytest.approx(r1, rel=1e-10)
