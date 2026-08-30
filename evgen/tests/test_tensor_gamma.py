"""The exact finite-gamma tensor kernel (plans/08 D2, Cosyn et al.).

`InclusiveKernel(tensor_gamma=True)` replaces the massless
Hoodbhoy-Jaffe-Manohar b-sector with Cosyn, Roldan Tomei, Sosa and Zec,
EPJ A 61 (2025) 83 (arXiv:2410.12764, refs/2410.12764v1.pdf) Eqs. (9),
(10), (14), (16), (17) and (24) -- the four tensor structure functions
with explicit geometry, the b3/b4 slots, and the rest-frame angle
theta_q that tilts the spin axis away from the photon direction and so
leaks the b-sector into cos 2phi.

Two kinds of test are kept apart here, because only one of them can
catch a transcription error.  `test_table1_row2_along_q` and
`test_table1_row3_along_the_beam` are anchored on a DERIVED result of the
paper -- the two finite-gamma rows of its Table 1, closed-form functions
of (gamma^2, eps, theta_q) that the module has never seen -- and it was
those two rows, and not the re-typing below, that caught the two
transcription errors this file's history records.  Everything else is a
second ROUTE through the same transcription: `_sfs_paper` writes Eqs.
(17) out again and `_harmonics_numeric` projects the cross-section weight
numerically instead of expanding it in closed form, which tests the
projection and the geometry but shares the module's reading of the
equations.
"""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen.xsec import (EventSpinState, InclusiveKernel,  # noqa: E402
                           cosyn_tensor_sfs, cosyn_unpolarized_sfs,
                           gamma_squared, theta_q_cos_sin)

from polli_fastsim import asymmetries as asym  # noqa: E402
from polli_fastsim import beams  # noqa: E402
from polli_fastsim.asymmetries import M_NUCLEON, epsilon_gamma  # noqa: E402
from polli_fastsim.polarized import toy_b1  # noqa: E402

CONFIG = beams.default_configs("6Li")[1]
S = CONFIG.sqrt_s_per_nucleon ** 2

# (x, Q2, y): a low-Q2 sweet spot of the low configuration, a mid one, and
# a deliberately large-gamma^2 point at the edge of the W^2 >= 10 cut
POINTS = [(0.08913, 1.135, 0.01561), (0.1413, 3.127, 0.02713),
          (0.2, 1.0, 0.05)]


def _tables(x, b1=0.02, b3=0.0, b4=0.0, r=0.3):
    """A hand-made SF table: F1 = 1, F2 = 2x(1+R)F1, b2 = 2x b1."""
    return {"f1": np.array([1.0]), "f2": np.array([2.0 * x * (1.0 + r)]),
            "b1": np.array([b1]), "b2": np.array([2.0 * x * b1]),
            "b3": np.array([b3]), "b4": np.array([b4]),
            "delta": np.array([0.0]), "g1": np.array([0.0]),
            "g2": np.array([0.0])}


# --- a second, literal transcription of the paper -------------------------

def _sfs_paper(b1, b2, b3, b4, x, g2):
    """Cosyn Eqs. (17a)-(17e), written out again from the PDF, with the
    repository's per-nucleon map x_d -> x.  The leading 2 of Eq. (17a)
    multiplies b1 alone: the large bracket opens before it."""
    g = np.sqrt(g2)
    a = 1.0 + g2
    f_tll_t = -(2.0 * a * b1 - (g2 / x) * ((1.0 / 6.0) * b2
                                           - (1.0 / 2.0) * b3))
    f_tll_l = (1.0 / x) * (2.0 * a * x * b1
                           - a ** 2 * ((1.0 / 3.0) * b2 + b3 + b4)
                           - a * ((1.0 / 3.0) * b2 - b4)
                           - ((1.0 / 3.0) * b2 - b3))
    f_tlt = -(g / (2.0 * x)) * (a * ((1.0 / 3.0) * b2 - b4)
                                + ((2.0 / 3.0) * b2 - 2.0 * b3))
    f_ttt = -(g2 / x) * ((1.0 / 6.0) * b2 - (1.0 / 2.0) * b3)
    return f_tll_t, f_tll_l, f_tlt, f_ttt


def _harmonics_numeric(t, x, q2, y, theta_s, m=1.0, j=1.0, n=8192):
    """(h0, h1, h2) by NUMERICAL projection of the tensor weight over the
    lepton azimuth -- the same physics as the module's closed form, by a
    different route.  Eqs. (9), (10), (14), (16), (17), (19), (24)."""
    g2 = 4.0 * M_NUCLEON ** 2 * x ** 2 / q2
    eps = ((1.0 - y - 0.25 * g2 * y * y)
           / (1.0 - y + 0.5 * y * y + 0.25 * g2 * y * y))
    root = np.sqrt(1.0 + g2)
    c = (1.0 + 0.5 * g2 * y) / root
    s = np.sqrt(g2 * (1.0 - y - 0.25 * g2 * y * y)) / root
    phi = 2.0 * np.pi * np.arange(n) / n
    ct, st = np.cos(theta_s), np.sin(theta_s)
    nx = c * st * np.cos(phi) + s * ct
    ny = st * np.sin(phi)
    nz = c * ct - s * st * np.cos(phi)
    q_nn = (3.0 * m * m - j * (j + 1.0)) / 3.0
    pref = 1.5 * q_nn                       # Eq. (9) for a pure state
    t_ll = pref * (nz * nz - 1.0 / 3.0)                       # Eq. (14a)
    t_lt = pref * nx * nz                                     # Eq. (14b)
    t_tt = pref * (nx * nx - ny * ny)                         # Eq. (14c)
    f_t, f_l, f_lt, f_tt = _sfs_paper(float(t["b1"]), float(t["b2"]),
                                      float(t["b3"]), float(t["b4"]), x, g2)
    f1, f2 = float(t["f1"]), float(t["f2"])
    den = 2.0 * f1 + eps * ((1.0 + g2) * f2 / x - 2.0 * f1)   # Eq. (16)
    num = (t_ll * (f_t + eps * f_l)
           + t_lt * np.sqrt(2.0 * eps * (1.0 + eps)) * f_lt
           + t_tt * eps * f_tt)                               # Eq. (19)
    w = -asym.TENSOR_LL_SIGN * num / den
    return (w.mean(), 2.0 * (w * np.cos(phi)).mean(),
            2.0 * (w * np.cos(2.0 * phi)).mean())


# --- the massless limit ---------------------------------------------------

def test_cosyn_sfs_reduce_to_the_hjm_b_sector_at_gamma_zero():
    """At gamma = 0, -(F_TLL_T + eps F_TLL_L)/(F_UU_T + eps F_UU_L) is
    EXACTLY the massless (b1 + (1-y)/(x y^2) b2)/D_phi of the master
    formula -- for any b2, and with b3, b4 cancelling identically."""
    for x, q2, y in POINTS:
        for b2s, b3, b4 in ((2.0, 0.0, 0.0), (1.3, 0.0, 0.0),
                            (2.0, 0.011, -0.007)):
            t = _tables(x, b3=b3, b4=b4)
            t["b2"] = np.array([b2s * x * float(t["b1"])])
            eps = epsilon_gamma(y, 0.0)
            f_t, f_l, _, f_tt = cosyn_tensor_sfs(
                t["b1"], t["b2"], t["b3"], t["b4"], x, 0.0)
            fu_t, fu_l = cosyn_unpolarized_sfs(t["f1"], t["f2"], x, 0.0)
            got = -(f_t + eps * f_l) / (fu_t + eps * fu_l)
            kern = (float(t["b1"]) + (1.0 - y) / (x * y * y) * float(t["b2"]))
            dphi = float(t["f1"]) + (1.0 - y) / (x * y * y) * float(t["f2"])
            assert float(got) == pytest.approx(kern / dphi, rel=1e-12)
            assert float(np.asarray(f_tt)) == 0.0


def test_theta_q_is_a_unit_vector_and_vanishes_masslessly():
    for x, q2, y in POINTS:
        g2 = gamma_squared(x, q2)
        c, s = theta_q_cos_sin(y, g2)
        assert float(c * c + s * s) == pytest.approx(1.0, rel=1e-14)
        assert float(s) > 0.0
    c0, s0 = theta_q_cos_sin(0.3, 0.0)
    assert (float(c0), float(s0)) == (1.0, 0.0)


def test_exact_path_reduces_to_the_massless_path_at_large_q2():
    """gamma^2 -> 0 by taking Q^2 -> infinity at fixed x: the two kernels
    must agree, and at a physical Q^2 they must differ by O(gamma^2)."""
    x = np.array([0.1])
    for q2v, tol in ((np.array([1.0e8]), 1e-9), (np.array([1.0e4]), 5e-5)):
        # s is chosen so that y stays at 1e-3 while gamma^2 -> 0; the
        # kernels see s only through y
        s_big = float(q2v) / (1.0e-3 * float(x))
        pair = [InclusiveKernel(beams.LI6, b1_func=toy_b1, tensor_gamma=tg)
                for tg in (False, True)]
        state = EventSpinState(0, 0.0, 1.0, 1.0, theta_s=0.0)
        w = [k.amplitudes(k.tables(x, q2v), x, q2v, s_big, state)[0]
             for k in pair]
        assert float(w[1]) == pytest.approx(float(w[0]), rel=tol)
        # and the residual really is O(gamma^2), coefficient of order 3
        g2v = float(gamma_squared(x, q2v))
        assert abs(float(w[1]) / float(w[0]) - 1.0) < 5.0 * g2v


def test_massless_path_is_untouched_by_the_new_slots():
    """The default kernel is bit-for-bit the pre-D2 one: b3/b4 tables are
    built but no massless amplitude reads them."""
    x, q2 = np.array([0.056]), np.array([1.14])
    y = q2 / (S * x)
    base = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    with_b34 = InclusiveKernel(beams.LI6, b1_func=toy_b1,
                               b3_func=lambda x, q2, f1: 0.3 * f1,
                               b4_func=lambda x, q2, f1: -0.2 * f1)
    state = EventSpinState(0, 0.0, 1.0, 1.0, theta_s=0.4)
    for k in (base, with_b34):
        t = k.tables(x, q2)
        got = k.amplitudes(t, x, q2, S, state)
        expected = (asym.TENSOR_LL_SIGN * (1.0 / 3.0)
                    * 0.5 * (3.0 * np.cos(0.4) ** 2 - 1.0)
                    * (t["b1"] + (1.0 - y) / (x * y * y) * t["b2"])
                    / (t["f1"] + (1.0 - y) / (x * y * y) * t["f2"]))
        assert float(got[0]) == float(expected)      # bit for bit
        assert float(got[1]) == 0.0


# --- the harmonics against an independent projection ----------------------

def test_harmonics_match_a_direct_numerical_projection():
    for x, q2, y in POINTS:
        t = _tables(x, b3=0.004, b4=-0.002)
        kern = InclusiveKernel(beams.LI6, tensor_gamma=True)
        for theta_s in (0.0, 0.5 * np.pi, 0.7, 2.3):
            for m in (1.0, 0.0, -1.0):
                state = EventSpinState(0, 0.0, 1.0, m, theta_s=theta_s)
                got = kern._tensor_harmonics_gamma(t, x, q2, y, state)
                want = _harmonics_numeric(t, x, q2, y, theta_s, m=m)
                for g, w in zip(got, want):
                    assert float(g) == pytest.approx(w, rel=1e-10,
                                                     abs=1e-15)


def test_the_transverse_fill_has_no_cos_phi_harmonic():
    """At theta_S = 90 deg every cos phi' coefficient carries sin theta_S
    cos theta_S and vanishes: the leakage is cos 2phi' only, which is why
    it lands on the observable."""
    t = _tables(0.1, b3=0.004, b4=-0.002)
    kern = InclusiveKernel(beams.LI6, tensor_gamma=True)
    state = EventSpinState(0, 0.0, 1.0, 1.0, theta_s=0.5 * np.pi)
    h0, h1, h2 = kern._tensor_harmonics_gamma(t, 0.1, 1.14, 0.03, state)
    assert abs(float(h1)) < 1e-18
    assert abs(float(h2)) > 0.0


# --- pinned numbers -------------------------------------------------------

# (h0, h1, h2) at theta_S = 90 deg, m = +1, for _tables(x) with b1 = 0.02,
# F1 = 1, R = 0.3, b3 = b4 = 0 -- the hand-derived combination of Eqs.
# (17d) and (17e) with the T_LL leakage, at the three POINTS.
PINNED = [(2.5543780132391e-03, 0.0, -3.0869622799178e-05),
          (2.5558792618530e-03, 0.0, -2.7677909400843e-05),
          (2.5180085726965e-03, 0.0, -1.6754330602946e-04)]


def test_pinned_values_at_three_points():
    kern = InclusiveKernel(beams.LI6, tensor_gamma=True)
    state = EventSpinState(0, 0.0, 1.0, 1.0, theta_s=0.5 * np.pi)
    for (x, q2, y), want in zip(POINTS, PINNED):
        got = kern._tensor_harmonics_gamma(_tables(x), x, q2, y, state)
        assert float(got[0]) == pytest.approx(want[0], rel=1e-11)
        assert float(got[2]) == pytest.approx(want[2], rel=1e-11)


def test_the_two_leading_channels_cancel_to_the_17e_term():
    """The three channels of the cos 2phi' harmonic stand
    T_LL : T_LT : T_TT = 3 : -3 : 1 as gamma^2 and y go to zero, so the
    leading-twist rate leakage and the twist-3 Eq. (17d) one cancel almost
    exactly and what survives is the twist-4 Eq. (17e) term alone.  That
    cancellation is the whole size of the effect: plans/08 D2 had guessed
    a factor 6.9 the other way, from the T_LT term ADDING."""
    def channels(x, q2, y):
        g2 = gamma_squared(x, q2)
        eps = epsilon_gamma(y, g2)
        c, s = theta_q_cos_sin(y, g2)
        t = _tables(x)
        f_t, f_l, f_lt, f_tt = cosyn_tensor_sfs(t["b1"], t["b2"], t["b3"],
                                                t["b4"], x, g2)
        ll = 0.5 * s * s * (f_t + eps * f_l)
        lt = -0.5 * c * s * np.sqrt(2.0 * eps * (1.0 + eps)) * f_lt
        tt = 0.5 * (c * c + 1.0) * eps * f_tt
        return float(ll / tt), float(lt / tt), float((ll + lt + tt) / tt)
    for x, q2, y in POINTS:
        r_ll, r_lt, r_tot = channels(x, q2, y)
        assert 2.9 < r_ll < 3.2 and -3.0 < r_lt < -2.8
        assert 0.9 < r_tot < 1.2
    r_ll, r_lt, r_tot = channels(0.05, 5.0, 1.0e-5)
    assert r_ll == pytest.approx(3.0, abs=5e-3)
    assert r_lt == pytest.approx(-3.0, abs=5e-3)
    assert r_tot == pytest.approx(1.0, abs=1e-2)


def test_table1_row2_along_q():
    """The paper's Table 1, second row: for b3 = b4 = 0, b2 = 2x b1 and a
    target polarized along q (Eq. 20a: T_LL = Q/3, T_LT = T_TT = 0),

        b1/(F1 A_T) = -9(1 + eps g2)/(6 + 5 g2 + 2 eps g2^2),  g2 = gamma^2.

    That row is a DERIVED result the module has never seen, and it is what
    fixes the bracketing of Eq. (17a): the leading 2 multiplies b1 alone.
    Applying it to the whole bracket gives 4 g2 for the 5 g2 and fails
    here by up to 2.5% over the range the W^2 cut allows."""
    for x, q2, y in POINTS:
        g2 = float(gamma_squared(x, q2))
        eps = float(epsilon_gamma(y, g2))
        t = _tables(x, r=0.0)                 # Table 1 uses F2 = 2x F1
        f_t, f_l, _, _ = cosyn_tensor_sfs(t["b1"], t["b2"], t["b3"],
                                          t["b4"], x, g2)
        fu_t, fu_l = cosyn_unpolarized_sfs(t["f1"], t["f2"], x, g2)
        a_t = (2.0 / 3.0) * float(f_t + eps * f_l) / float(fu_t + eps * fu_l)
        printed = -9.0 * (1.0 + eps * g2) / (6.0 + 5.0 * g2
                                             + 2.0 * eps * g2 * g2)
        assert float(t["b1"]) / (float(t["f1"]) * a_t) == pytest.approx(
            printed, rel=1e-10)


def test_table1_row3_along_the_beam():
    """The paper's Table 1, third row: the same quantity for a target
    polarized along the BEAM, which the module reaches at theta_S = 0.
    This row is the stronger anchor of the two -- it carries F[U T_LT] and
    F[U T_TT] and the whole photon-frame geometry, and it is what fixes
    the frame: the incoming lepton sits at +sin theta_q, so that
    T_LT cos phi_TL = +(Q/2) cos theta_q sin theta_q for N = N_e.  Eq.
    (22b) as printed has the opposite sign and misses this row by up to a
    factor of three at the largest gamma^2 the W^2 cut allows."""
    kern = InclusiveKernel(beams.LI6, tensor_gamma=True)
    state = EventSpinState(0, 0.0, 1.0, 1.0, theta_s=0.0)   # N = N_e, Q = 1
    for x, q2, y in POINTS:
        g2 = float(gamma_squared(x, q2))
        g = np.sqrt(g2)
        eps = float(epsilon_gamma(y, g2))
        c, sn = (float(v) for v in theta_q_cos_sin(y, g2))
        t = _tables(x, r=0.0)
        h0, h1, h2 = kern._tensor_harmonics_gamma(t, np.array([x]),
                                                  np.array([q2]),
                                                  np.array([y]), state)
        # w = -TENSOR_LL_SIGN * (numerator/denominator) and A_T is twice
        # that numerator over the same denominator, at Q = 1
        a_t = -2.0 * float(h0) / asym.TENSOR_LL_SIGN
        printed = 18.0 * (1.0 + eps * g2) / (
            (3.0 * c * c - 1.0) * (-6.0 - 5.0 * g2 - 2.0 * eps * g2 * g2)
            - c * sn * np.sqrt(2.0 * eps * (1.0 + eps)) * (9.0 * g
                                                           + 3.0 * g ** 3)
            - sn * sn * eps * 3.0 * g2)
        assert float(t["b1"]) / (float(t["f1"]) * a_t) == pytest.approx(
            printed, rel=1e-10)
        assert float(h1) == 0.0 and float(h2) == 0.0


def test_b3_b4_default_to_zero_and_reach_the_kernel():
    x, q2 = np.array([0.1]), np.array([1.14])
    t0 = InclusiveKernel(beams.LI6, b1_func=toy_b1).tables(x, q2)
    assert float(t0["b3"]) == 0.0 and float(t0["b4"]) == 0.0
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1, tensor_gamma=True,
                           b3_func=lambda x, q2, f1: 0.05 * f1)
    y = q2 / (S * x)
    state = EventSpinState(0, 0.0, 1.0, 1.0, theta_s=0.5 * np.pi)
    t = kern.tables(x, q2)
    assert float(t["b3"]) != 0.0
    base = InclusiveKernel(beams.LI6, b1_func=toy_b1, tensor_gamma=True)
    a2 = kern.amplitudes(t, x, q2, S, state)[2]
    a2_0 = base.amplitudes(base.tables(x, q2), x, q2, S, state)[2]
    assert float(a2) != float(a2_0)
