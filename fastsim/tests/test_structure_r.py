"""R = sigma_L/sigma_T: the published R1998 fit, and the r_func hook.

Two things are pinned here.

(1) `structure.r1998` against the E143 measurement it comes from --
K. Abe et al., Phys. Lett. B 452 (1999) 194 (arXiv:hep-ex/9808028,
refs/hep-ex_9808028.pdf).  The anchor is their Table I: nine measured
values of R with statistical and systematic errors at 0.0325 <= x <=
0.105, 1.32 <= Q2 <= 2.73 GeV2.  That is exactly the corner where the
factor Theta(x, Q2) of their Eq. (3) is large (6.15-7.64 over the nine
points, measured), so the data
themselves decide where Theta belongs: the published placement (log term
only) sits on top of the measurements (max |pull| 1.56), while the
fast-sim scripts' placement (all three terms) overshoots every one of
them by 5.2-8.0 sigma AFTER being clipped to R = 1.  That is code
review 2026-08-25 item S1 / plans/08 C2, and it is checked below against
the numbers, not against the code.

(2) The `r_func` hook of `NuclearF2`, `dsigma_dx_dq2`,
`depolarization_d` and `ToyG1`: r_func=None must be bit-for-bit the
status quo, and a constant R must reproduce F1 = F2/(2x(1+R)) and
D(y) = y(2-y)/(y^2 + 2(1-y)(1+R)) written out by hand.
"""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import asymmetries as asym  # noqa: E402
from polli_fastsim import beams, polarized  # noqa: E402
from polli_fastsim import structure as st  # noqa: E402

# Abe et al. (arXiv:hep-ex/9808028) Table I: x, Q2 [GeV2], R, stat, sys.
E143_TABLE_I = [
    (0.0325, 1.32, 0.45, 0.01, 0.07),
    (0.0375, 1.47, 0.51, 0.02, 0.09),
    (0.0450, 1.67, 0.40, 0.01, 0.10),
    (0.0550, 1.90, 0.28, 0.01, 0.09),
    (0.0650, 2.11, 0.29, 0.02, 0.10),
    (0.0750, 2.29, 0.18, 0.02, 0.11),
    (0.0850, 2.46, 0.26, 0.03, 0.10),
    (0.0950, 2.60, 0.25, 0.03, 0.14),
    (0.1050, 2.73, 0.17, 0.03, 0.13),
]

# delta-R of Abe et al. (the unnumbered equation below their Eq. 3)
# evaluated BY HAND at six (x, Q2); see the test for the arithmetic.
# Q2 = 0.5 is the fit's lower support edge, where delta-R is largest.
DELTA_R_BY_HAND = [
    (0.005, 0.5, 0.03867477),        # 0.007735 + 0.0680675/2.2
    (0.040, 0.5, 0.03251636),        # 0.007280 + 0.0555200/2.2
    (0.299, 0.5, 0.01117241),        # 0.003913 + 0.0159707/2.2 (the minimum)
    (0.860, 0.5, 0.11131091),        # -0.003380 + 0.2523200/2.2
    (0.050, 20.0, 0.00955783),       # 0.007150 + 0.0522500/21.7
    (0.0562341, 1.13506, 0.02480484),  # sweet spot 1 of the cos 2phi map
]

# The five (x, Q2) points of code review 2026-08-25 item S1, with the two
# numbers it quotes at each: the fast-sim scripts' function as written
# (Theta in all three terms, pre-clip) and the same coefficients with
# Theta restricted to the log term.
S1_POINTS = [
    (0.005, 2.00, 2.17, 0.34),
    (0.010, 2.50, 2.04, 0.31),
    (0.050, 1.14, 2.17, 0.41),
    (0.100, 5.00, 1.00, 0.20),
    (0.300, 7.40, 0.27, 0.12),
]


def _script_r1998(x, q2, theta_everywhere=True):
    """`money_delta_realistic.r1998` (line 121), transcribed verbatim.

    `theta_everywhere=False` is the code review's repair -- Theta on the
    log term only, the other two coefficients left where the script put
    them.  Neither is Abe et al. Eq. (2): the script's second term is
    0.5470/(Q^4 + 2^4)^(1/2) rather than a2 (Q^8 + a3^4)^(-1/4) times
    [1 + a4 x + a5 x^2] x^a6, and its third term
    0.2379/(Q2 + 0.3^2) [1 - 0.00815 ln x / ln 0.5] has no counterpart in
    the paper at all.
    """
    x = np.asarray(x, dtype=float)
    q2 = np.asarray(q2, dtype=float)
    log_q2 = np.log(np.maximum(q2, 1.01) / 0.04)
    theta = 1.0 + 12.0 * q2 / (q2 + 1.0) * 0.125**2 / (0.125**2 + x**2)
    other = theta if theta_everywhere else 1.0
    part1 = 0.0485 * theta / log_q2
    part2 = 0.5470 * other / (q2**2 + 2.0**4) ** 0.5
    part3 = 0.2379 * other / (q2 + 0.3**2) * (
        1.0 - 0.00815 * np.log(np.maximum(x, 1e-12)) / np.log(0.5))
    return part1 + part2 + part3


def _pulls(model_r):
    """(R_model - R_meas)/sigma over the nine E143 Table I rows."""
    out = []
    for x, q2, r, stat, sys_ in E143_TABLE_I:
        sig = np.hypot(stat, sys_)
        out.append((float(np.asarray(model_r(x, q2))) - r) / sig)
    return np.array(out)


# --- (a) the published measurement ------------------------------------------


def test_r1998_reproduces_the_e143_measurements():
    """R1998 sits on the nine measured points of Abe et al. Table I."""
    pulls = _pulls(st.r1998)
    assert np.max(np.abs(pulls)) < 2.0, (
        "R1998 must reproduce its own input data (Abe et al., PLB 452 "
        "(1999) 194, Table I); pulls = %s"
        % np.array2string(pulls, precision=2))
    # the fit was made TO these points, so the chi2 per point is < 1
    assert float(np.mean(pulls**2)) < 1.0, np.mean(pulls**2)


def test_r1998_parameters_are_table_ii():
    """Transcription guard on Abe et al. Table II.

    The code took the coefficients from the published PDF; this copy was
    read out of the paper's own LaTeX source (arXiv e-print of
    hep-ex/9808028, file RTABLE_COEF666.tex), so a slipped digit in
    either transcription shows up here rather than as a 5% shift in a
    money plot.  The chi2/df column of that table is 0.9 / 0.9 / 1.0 for
    231 degrees of freedom.
    """
    assert st.R1998_PARAMS == {
        "a": (0.0485, 0.5470, 2.0621, -0.3804, 0.5090, -0.0285),
        "b": (0.0481, 0.6114, -0.3509, -0.4611, 0.7172, -0.0317),
        "c": (0.0577, 0.4644, 1.8288, 12.3708, -43.1043, 41.7415),
    }
    assert st.R1998_X_RANGE == (0.005, 0.86)      # paper, above Eq. (2)
    assert st.R1998_Q2_RANGE == (0.5, 130.0)


def test_the_three_forms_agree_within_the_fit_error():
    """R_a, R_b and R_c were fitted to the SAME 237 world data points,
    with chi2/df = 0.9, 0.9 and 1.0 (Abe et al. Table II).  Three fits
    that each describe the same data that well cannot disagree with one
    another by much more than the fit error delta-R where the data are;
    a mistyped coefficient or a wrong power of Q2 in any ONE of them
    breaks that agreement even when the average still looks plausible.
    """
    x = np.geomspace(0.01, 0.5, 50)[:, None]
    q2 = np.geomspace(1.0, 50.0, 50)[None, :]
    ratio = st.r1998_spread(x, q2) / np.broadcast_to(
        st.r1998_fit_error(x, q2), (50, 50))
    assert ratio.mean() < 1.3 and ratio.max() < 3.0, (
        "spread of R_a/R_b/R_c over delta_R: mean %.3f, max %.2f "
        "(published transcription gives 1.011 and 2.72)"
        % (ratio.mean(), ratio.max()))


def test_r1998_in_the_world_data_range_over_the_low_x_dis_region():
    """0.2 <= R <= 0.35 is the world-data band quoted in the code review
    for the low-x DIS region; R1998 must bracket it, not sit at 1."""
    x = np.geomspace(0.005, 0.1, 40)[:, None]
    q2 = np.geomspace(1.0, 5.0, 40)[None, :]
    r = st.r1998(x, q2)
    assert 0.15 < r.min() and r.max() < 0.45, (
        "R1998 over 0.005<x<0.1, 1<Q2<5 spans [%.3f, %.3f]; the world data "
        "for that region are R = 0.2-0.35 (code review 2026-08-25 S1)"
        % (r.min(), r.max()))


def test_r1998_positive_and_never_saturates_in_the_acceptance():
    """The defect's signature was R = 1.000 exactly; the fit must not go
    anywhere near 1 over the whole 6Li acceptance."""
    x = np.geomspace(1e-4, 0.8, 90)[:, None]
    q2 = np.geomspace(1.0, 100.0, 90)[None, :]
    r = st.r1998(x, q2)
    assert r.min() > 0.0
    assert r.max() < 0.5, r.max()
    # the same box through the function this replaces
    bad = np.clip(_script_r1998(x, q2), 0.0, 1.0)
    assert np.mean(bad >= 1.0) > 0.3, (
        "the money-script form is expected to be clipped to R = 1.000 over "
        "a large part of the acceptance; got %.3f" % np.mean(bad >= 1.0))


# --- (b) the bug that motivated plans/08 C2 ---------------------------------


def test_theta_placement_is_decided_by_the_e143_data():
    """Theta on all three terms is excluded by the measurements.

    Independent of any transcription of Eq. (2): put Theta everywhere and
    the SAME coefficients miss the nine measured R values of Table I --
    where Theta = 6.15-7.64 -- by 5.2 to 8.0 sigma (measured),
    every one of them high, even after the clip to R <= 1 that hides the
    divergence.
    """
    bad = _pulls(lambda x, q2: np.clip(_script_r1998(x, q2), 0.0, 1.0))
    assert np.min(bad) > 5.0, (
        "Theta in all three terms should be excluded by Abe et al. Table I; "
        "pulls = %s" % np.array2string(bad, precision=1))
    assert np.max(np.abs(_pulls(st.r1998))) < np.min(bad)


def test_code_review_s1_points():
    """The five (x, Q2) of code review S1, with both sets of numbers.

    The review's arithmetic is reproduced exactly -- for the function it
    evaluated.  That function is not R1998: at x >= 0.1 its extra
    0.2379/(Q2+0.09) term is 23-26% of the total, which is why the full
    three-form fit lands lower than the review's 0.20 and 0.12 by more
    than the R_a/R_b/R_c spread.  What the code now uses is R1998.
    """
    for x, q2, script_quoted, review_quoted in S1_POINTS:
        assert abs(float(_script_r1998(x, q2)) - script_quoted) < 0.01
        repaired = float(_script_r1998(x, q2, theta_everywhere=False))
        assert abs(repaired - review_quoted) < 0.01, (x, q2, repaired)
        r = float(st.r1998(x, q2))
        spread = float(st.r1998_spread(x, q2))
        assert 0.05 < r < 0.45, (x, q2, r)
        # where the review's simplified form and R1998 disagree by more
        # than the functional-form spread, the culprit is the review
        # function's extra term, which grows relative to the rest with x
        if abs(r - review_quoted) > spread:
            assert x >= 0.1, (
                "unexpected R1998 vs code-review disagreement at low x: "
                "x=%.3f Q2=%.2f R1998=%.3f review=%.2f spread=%.3f"
                % (x, q2, r, review_quoted, spread))


def test_r1998_forms_average_and_spread():
    x = np.array([0.01, 0.1, 0.5])[:, None]
    q2 = np.array([1.0, 10.0, 100.0])[None, :]
    forms = st.r1998_forms(x, q2)
    np.testing.assert_allclose(st.r1998(x, q2),
                               sum(forms) / 3.0, rtol=0, atol=0)
    for name, want in zip("abc", forms):
        np.testing.assert_array_equal(st.r1998(x, q2, form=name), want)
    stacked = np.stack(np.broadcast_arrays(*forms))
    np.testing.assert_array_equal(st.r1998_spread(x, q2),
                                  stacked.max(0) - stacked.min(0))
    with pytest.raises(ValueError):
        st.r1998(0.1, 10.0, form="d")


def test_r1998_outside_the_fit_support_is_frozen_not_extrapolated():
    lo_x, hi_x = st.R1998_X_RANGE
    lo_q2, hi_q2 = st.R1998_Q2_RANGE
    assert st.r1998(1e-4, 2.0) == st.r1998(lo_x, 2.0)
    assert st.r1998(0.95, 2.0) == st.r1998(hi_x, 2.0)
    assert st.r1998(0.1, 0.2) == st.r1998(0.1, lo_q2)
    assert st.r1998(0.1, 500.0) == st.r1998(0.1, hi_q2)
    # clip=False follows the analytic form instead; below x = 0.005 Theta
    # has saturated, so the two differ by only a few per cent
    raw = float(st.r1998(1e-4, 2.0, clip=False))
    assert 1.0 < raw / float(st.r1998(1e-4, 2.0)) < 1.12


def test_r1998_fit_error_is_the_published_delta_r():
    """delta-R evaluated by hand from the published expression.

    Abe et al. (arXiv:hep-ex/9808028), unnumbered equation below their
    Eq. (3):

        delta-R(x, Q2) = 0.0078 - 0.013 x
                         + (0.070 - 0.39 x + 0.70 x^2)/(1.7 + Q2)

    The expectations in DELTA_R_BY_HAND were worked out on paper from
    that expression, digit by digit, e.g. at (x, Q2) = (0.040, 0.5):
    0.0078 - 0.00052 = 0.007280, and (0.070 - 0.01560 + 0.001120)/2.2 =
    0.055520/2.2 = 0.02523636, summing to 0.03251636.  Two-sided to 8
    decimals, so a slipped digit or a flipped sign in ANY of the five
    constants shows up here -- the previous version of this test asserted
    only 0 < delta-R < hypot(stat, sys) and monotonicity in Q2, and three
    corruptions of those constants (0.070 -> 0.140, -0.013 -> +0.013,
    -0.39 -> +0.39) passed the whole fastsim suite.

    delta-R is the fit error the Delta/F1 bands of plans/08 C2 are built
    from, so it has to be pinned as tightly as the Table II coefficients.
    """
    for x, q2, want in DELTA_R_BY_HAND:
        got = float(st.r1998_fit_error(x, q2))
        assert abs(got - want) < 5e-9, (x, q2, got, want)

    # The turning point in x pins the -0.013 x, -0.39 x and +0.70 x^2
    # terms JOINTLY and independently of the overall size: setting
    # d(delta-R)/dx = 0 gives 1.4 x = 0.39 + 0.013 (1.7 + Q2), i.e.
    # x = 0.299000 at Q2 = 0.5 and x = 0.340786 at Q2 = 5.  (It is a
    # MINIMUM, which is the disagreement with the paper's own sentence
    # "largest at low Q2, reaching a maximum value for x ~ 0.04" recorded
    # in the refs_dict entry: the formula has no maximum near x = 0.04.)
    grid = np.linspace(0.005, 0.86, 855001)
    for q2, want_x in ((0.5, 0.299000), (5.0, 0.340786)):
        got_x = float(grid[np.argmin(st.r1998_fit_error(grid, q2))])
        assert abs(got_x - want_x) < 2e-6, (q2, got_x, want_x)
        assert st.r1998_fit_error(0.04, q2) > st.r1998_fit_error(want_x, q2)

    # ...and the sanity the fit itself has to satisfy: delta-R is small
    # against the errors of the measurements it was fitted to, and falls
    # with Q2 (the 1/(1.7 + Q2) denominator).
    for x, q2, _r, stat, sys_ in E143_TABLE_I:
        d = float(st.r1998_fit_error(x, q2))
        assert 0.0 < d < np.hypot(stat, sys_)
    assert (st.r1998_fit_error(0.05, 2.0)
            > st.r1998_fit_error(0.05, 20.0)
            > st.r1998_fit_error(0.05, 100.0))


# --- (c) the r_func hook ----------------------------------------------------


def _dense():
    x = np.geomspace(1e-4, 0.8, 47)[:, None]
    q2 = np.geomspace(0.3, 300.0, 31)[None, :]
    return x, q2


def test_r_func_none_is_bit_for_bit_the_toy():
    """Every hooked consumer with r_func=r_sigma_lt == r_func=None."""
    x, q2 = _dense()
    s = 4.0 * 10.0 * 50.0
    y = q2 / (s * x)
    ion = beams.LI6
    pairs = [
        (st.NuclearF2(ion).f1a(x, q2),
         st.NuclearF2(ion, r_func=st.r_sigma_lt).f1a(x, q2)),
        (st.dsigma_dx_dq2(x, q2, s, st.NuclearF2(ion).f2a(x, q2) / ion.A),
         st.dsigma_dx_dq2(x, q2, s, st.NuclearF2(ion).f2a(x, q2) / ion.A,
                          r_func=st.r_sigma_lt)),
        (asym.depolarization_d(y, x, q2),
         asym.depolarization_d(y, x, q2, r_func=st.r_sigma_lt)),
        (polarized.ToyG1().g1p(x, q2),
         polarized.ToyG1(r_func=st.r_sigma_lt).g1p(x, q2)),
        (polarized.ToyG1().g1_nucleus(ion, x, q2),
         polarized.ToyG1(r_func=st.r_sigma_lt).g1_nucleus(ion, x, q2)),
    ]
    for default, hooked in pairs:
        np.testing.assert_array_equal(default, hooked)


def test_constant_r_reproduces_the_analytic_forms():
    """A constant R gives F1 = F2/(2x(1+R)) and the textbook D(y)."""
    x, q2 = _dense()
    s = 4.0 * 10.0 * 50.0
    y = q2 / (s * x)
    r0 = 0.25

    def const_r(xx, qq):
        return np.full(np.broadcast(np.asarray(xx), np.asarray(qq)).shape, r0)

    nf2 = st.NuclearF2(beams.LI6, r_func=const_r)
    np.testing.assert_allclose(nf2.f1a(x, q2),
                               nf2.f2a(x, q2) / (2.0 * x * (1.0 + r0)),
                               rtol=1e-14)
    np.testing.assert_allclose(
        asym.depolarization_d(y, x, q2, r_func=const_r),
        y * (2.0 - y) / (y * y + 2.0 * (1.0 - y) * (1.0 + r0)), rtol=1e-14)
    f2 = nf2.f2a(x, q2) / beams.LI6.A
    np.testing.assert_array_equal(
        st.dsigma_dx_dq2(x, q2, s, f2, r_func=const_r),
        st.dsigma_dx_dq2(x, q2, s, f2, fl=f2 * r0 / (1.0 + r0)))
    g1 = polarized.ToyG1(r_func=const_r)
    base = g1.base
    np.testing.assert_allclose(
        g1.g1p(x, q2),
        g1.a1p(x) * base.f2p(x, q2) / (2.0 * x * (1.0 + r0)), rtol=1e-14)


def test_r1998_hook_moves_f1_by_the_expected_factor():
    """Swapping the toy R for R1998 rescales F1 by (1+R_toy)/(1+R_1998)."""
    x, q2 = _dense()
    ion = beams.LI6
    toy = st.NuclearF2(ion).f1a(x, q2)
    fit = st.NuclearF2(ion, r_func=st.r1998).f1a(x, q2)
    expected = (1.0 + st.r_sigma_lt(x, q2)) / (1.0 + st.r1998(x, q2))
    np.testing.assert_allclose(fit / toy, np.broadcast_to(expected, fit.shape),
                               rtol=1e-14)
    assert np.any(np.abs(fit / toy - 1.0) > 0.1)  # the swap is not cosmetic


def test_dated_script_monkey_patch_still_reaches_three_of_four():
    """The frozen `r_override` contract of fastsim/scripts/money_delta_*.

    Those scripts rebind `structure.r_sigma_lt` and
    `asymmetries.r_sigma_lt`, so the module-global lookup has to stay
    inside each consumer.  `polarized.ToyG1` is the fourth consumer they
    never reach -- that is the defect plans/08 C2 records, kept here as a
    characterization so a future refactor cannot change the dated
    scripts' output without tripping a test.
    """
    def const_r(xx, qq):
        shape = np.broadcast(np.asarray(xx), np.asarray(qq)).shape
        return np.full(shape, 0.42)

    x, q2, s = 0.05, 2.0, 2000.0
    nf2, g1 = st.NuclearF2(beams.LI6), polarized.ToyG1()
    saved = (st.r_sigma_lt, asym.r_sigma_lt)
    before = (float(nf2.f1a(x, q2)), float(asym.depolarization_d(0.3, x, q2)),
              float(st.dsigma_dx_dq2(x, q2, s, nf2.f2a(x, q2) / 6.0)),
              float(g1.g1p(x, q2)))
    try:
        st.r_sigma_lt, asym.r_sigma_lt = const_r, const_r
        after = (float(nf2.f1a(x, q2)),
                 float(asym.depolarization_d(0.3, x, q2)),
                 float(st.dsigma_dx_dq2(x, q2, s, nf2.f2a(x, q2) / 6.0)),
                 float(g1.g1p(x, q2)))
    finally:
        st.r_sigma_lt, asym.r_sigma_lt = saved
    for b, a in zip(before[:3], after[:3]):
        assert b != a          # f1a, depolarization_d, dsigma_dx_dq2 follow
    assert before[3] == after[3]   # ToyG1 does not: the missed consumer
    # ...and the hook is the way to move it
    assert float(polarized.ToyG1(r_func=const_r).g1p(x, q2)) != before[3]
