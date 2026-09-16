"""The tensor sign conventions: the tensor kernel's, and the wave function's.

Motivated by ``docs/consistency_review_2026-09-02.md`` section "5.2 Mechanical
checks tools/consistency_check.py could add", item 8 ("Sign-convention
guard"), which covers F013 and F002.

Development run 16 (2026-08-29) moved the simulation to the literature
convention, ``polli_fastsim.asymmetries.TENSOR_LL_SIGN = -1``, the sign every
published b1 is quoted in.  The master formula the reports display is written
in this programme's own transcription of Hoodbhoy-Jaffe-Manohar, whose b1, b2
terms carry the other sign, so a reader who takes the displayed formula for
what the code evaluates gets A_zz backwards.  Report 2 displayed it without
the caveat until the consistency review; Reports 0 and 1 carried it.

Four assertions:

  * the constant is still -1 -- if it is flipped back, every "opposite"
    sentence in the reports becomes wrong at once;
  * every template that displays the master formula states the convention
    within twelve lines of the display, either by calling the simulation's
    sign the "opposite" one or by printing the literature relation
    A_zz = -(2/3) b1/F1;
  * the O(gamma^2) tensor-leakage correction run 18 added
    (``polligen.xsec.tensor_leakage_amplitude``, plans/08 D2 / D10)
    REVERSES with the same constant.  The constant is the generator's, in
    ``polligen.xsec``, not ``polli_fastsim.asymmetries``', and nothing in
    the first two assertions reaches it; without this one the two halves
    of the convention could drift apart and the reports' "it cancels part
    of a negative Delta rather than faking one" would silently invert.
    The cheapest form of the statement is the one
    ``evgen/tests/test_tensor_gamma.py::
    test_the_correction_reverses_with_the_tensor_sign_constant`` makes at
    one sweet spot: both harmonics carry -TENSOR_LL_SIGN, so flipping it
    inside a try/finally flips h2 and leaves the kinematic ratio
    L = h2/h0 -- what the subtraction multiplies the fitted kappa by --
    untouched;
  * the OTHER sign in the tensor sector, the one inside the wave function,
    is still there.  Development run 19 (2026-09-15) found that
    ``polligen.tagged.TaggedModel._amp2_table`` summed its partial waves
    without the i^L phase of phi_L = i^L psi_L, so the S-D interference
    that carries the whole tagged tensor asymmetry entered with the wrong
    relative sign: A_zz^wf(90 deg) came out negative, the axial node of
    the spectator density sat in M = 0 instead of M = +-1, and the folded
    money-plot asymmetries read +0.49 / -0.07 where the corrected model
    reads -0.84 / +0.22.  The defect is invisible to every angle-
    integrated quantity (``tensor_dilution`` moved by 4.3e-6 relative, the
    k-marginals not at all), so nothing else in this sweep would see it
    come back.  The statement that does is Cosyn-Weiss II
    (arXiv:2603.23700) Eq. (6.12),

        A_T|| = (2 f0 + f2/sqrt2)(f2/sqrt2)/(f0^2 + f2^2) (1 - 3 cos^2 th),

    which is spin algebra alone -- it holds for ANY pair of radial
    functions -- so the module's own A_zz^wf must equal it to double
    precision when f0, f2 are the channel's own stored tables, under the
    mapping A_T|| = +1 x A_zz^wf.  Dropping the phase is exactly
    f2 -> -f2 in that expression, which is a different number at every
    point where f2 is not small, so the identity is a one-sided gate on
    the phase.  ``evgen/tests/test_tagged.py`` makes the same statement
    over the whole (k, cos th_k) grid and adds CW's TABLE II on the
    tabulated AV18 deuteron; this check is the three-point form of it, so
    that the sweep catches the regression without the test suite.
"""

import re

ROOT = checker.ROOT

# the b1 sector of the master formula, in the TeX the templates embed
# (data-tex="... F_1 + \frac{2}{3}a_m b_1 ...") and in a plain-text or
# MathJax spelling of the same display.
_MASTER = re.compile(
    r"\\frac\{2\}\{3\}\s*a_m\s*b_1"
    r"|\(2\s*/\s*3\)\s*a_m\s*b[_₁1]"
    r"|\\tfrac\{2\}\{3\}\s*a_m\s*b_1")
# the caveat: the simulation carries the "opposite" sign, or the literature
# relation is printed outright.
_CAVEAT = re.compile(
    r"opposite"
    r"|A\s*zz\s*=\s*[−-]\s*\(?\s*2\s*/\s*3\s*\)?\s*b[₁1]\s*/\s*F[₁1]"
    r"|[−-]\s*\(2/3\)\s*b[₁1]/F[₁1]")
_WINDOW = 12


def _text(html):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html)).strip()


@check("physics: the tensor kernel carries the literature sign TENSOR_LL_SIGN = -1")
def _():
    from polli_fastsim import asymmetries
    sign = getattr(asymmetries, "TENSOR_LL_SIGN", None)
    if sign is None:
        return ["polli_fastsim.asymmetries has no TENSOR_LL_SIGN; the b1 sector "
                "lost the constant run 16 (2026-08-29) introduced"]
    if float(sign) != -1.0:
        return ["polli_fastsim.asymmetries.TENSOR_LL_SIGN = %r, not -1: the "
                "simulation has left the convention every published b1 is "
                "quoted in (Azz = -(2/3) b1/F1), which Reports 0, 1 and 2 "
                "state at the master formula" % (sign,)]
    return []


@check("physics: every template displaying the master formula states the sign convention")
def _():
    bad = []
    for tpl in sorted(ROOT.glob("reports/*.template.html")):
        lines = tpl.read_text().splitlines()
        for i, raw in enumerate(lines):
            if not _MASTER.search(raw):
                continue
            lo, hi = max(0, i - _WINDOW), min(len(lines), i + _WINDOW + 1)
            near = _text(" ".join(lines[lo:hi]))
            if not _CAVEAT.search(near):
                bad.append("%s:%d displays the b1 master formula with no sign "
                           "caveat within %d lines -- the simulation carries "
                           "TENSOR_LL_SIGN = -1, the opposite of the displayed "
                           "b1, b2 terms (Azz = -(2/3) b1/F1)"
                           % (tpl.relative_to(ROOT), i + 1, _WINDOW))
    return bad


# one of the twelve money-plot sweet spots of `tensor_gamma_leakage.py`
# (5 x 40.8, the worst of them): (x, Q2, y).  One point is enough -- the
# statement is an algebraic identity in the constant, not a sampling of it.
_LEAK_POINT = (0.08913, 1.135, 0.01561)


@check("physics: the O(gamma^2) tensor-leakage correction reverses with "
       "TENSOR_LL_SIGN")
def _():
    import numpy as np
    from polligen import xsec as xs

    x, q2, y = (np.array([v]) for v in _LEAK_POINT)
    one = np.ones(1)
    # b1 = 1 with the tensor Callan-Gross b2 = 2x b1 and the higher-twist
    # slots at zero: the same table `xsec.tensor_leakage_ratio` builds, and
    # the ratio L is independent of it (tripling b1 leaves L at six digits).
    t = {"b1": one, "b2": 2.0 * x * one, "b3": 0.0 * one, "b4": 0.0 * one,
         "f1": one, "f2": 2.0 * x * one}

    sign = getattr(xs, "TENSOR_LL_SIGN", None)
    if sign is None:
        return ["polligen.xsec has no TENSOR_LL_SIGN: the finite-gamma tensor "
                "kernel has lost the constant the leakage correction reverses "
                "with (plans/08 D1/D2)"]
    h0, h2 = xs.tensor_leakage_amplitude(t, x, q2, y)
    ell = xs.tensor_leakage_ratio(x, q2, y)
    try:
        xs.TENSOR_LL_SIGN = -sign
        h0_f, h2_f = xs.tensor_leakage_amplitude(t, x, q2, y)
        ell_f = xs.tensor_leakage_ratio(x, q2, y)
    finally:
        xs.TENSOR_LL_SIGN = sign

    bad = []
    a, b = float(h2.item()), float(h2_f.item())
    if not (a != 0.0 and abs(b + a) <= 1e-12 * abs(a)):
        bad.append("polligen.xsec.tensor_leakage_amplitude: h2 = %.6g at "
                   "(x, Q2, y) = %s does not reverse when TENSOR_LL_SIGN is "
                   "flipped (it gives %.6g, not %.6g) -- the leakage "
                   "correction no longer follows the tensor convention "
                   "plans/08 D1 adopted" % (a, _LEAK_POINT, b, -a))
    c, d = float(ell.item()), float(ell_f.item())
    if not (c != 0.0 and abs(d - c) <= 1e-12 * abs(c)):
        bad.append("polligen.xsec.tensor_leakage_ratio: L = %.6g moves to "
                   "%.6g when TENSOR_LL_SIGN is flipped -- L is the ratio of "
                   "two harmonics that both carry the constant, so it must "
                   "not depend on it (the subtraction A - L kappa would "
                   "otherwise change size, not just sign)" % (c, d))
    if a >= 0.0 or float(h0.item()) <= 0.0:
        bad.append("polligen.xsec: at the worst sweet spot the leakage h2 = "
                   "%.6g and the constant h0 = %.6g; with the literature "
                   "convention (TENSOR_LL_SIGN = -1) and a positive b1 the "
                   "reports state h2 < 0 < h0, i.e. that the leakage cancels "
                   "part of a negative Delta rather than faking one"
                   % (a, float(h0.item())))
    return bad


# three (k, cos theta_k) cells of the model's own grid, chosen away from the
# P2 zero and spanning the tagged window: the identity is exact everywhere,
# so three points are a spot check of an algebraic statement, not a sample.
_CW_POINTS = ((0.15, -0.9), (0.30, 0.0), (0.50, 0.9))
# below this the phase-less alternative (f2 -> -f2) is too close to the
# corrected value for the point to discriminate; measured separations are
# 1.06 / 1.32 / 2.07 (deuteron) and 1.24 / 1.39 / 2.03 (6Li).  The third
# cell was (0.50, 0.6) until 2026-09-15, one grid cell from the P2 zero at
# cos th_k = 1/sqrt(3) = 0.5774, where the separation is only 0.08 against
# this 0.05: a change of nc or k_max that moved that cell would have made
# the check report "move the point" rather than a regression
_CW_GAP = 0.05


@check("physics: the tagged S-D interference carries the i^L partial-wave "
       "phase (Cosyn-Weiss Eq. 6.12)")
def _():
    import numpy as np
    from polligen import tagged

    s2 = np.sqrt(2.0)
    bad = []
    for name, channel in (("deuteron", tagged.deuteron_channel),
                          ("6Li alpha-d", tagged.li6_alpha_channel)):
        m = tagged.TaggedModel(channel())
        if 0 not in m._rad or 2 not in m._rad:
            bad.append("polligen.tagged: the %s channel no longer carries an "
                       "S and a D wave, so the Cosyn-Weiss Eq. (6.12) "
                       "identity cannot be formed" % name)
            continue
        for k_t, c_t in _CW_POINTS:
            ik = int(np.argmin(np.abs(m.k - k_t)))
            ic = int(np.argmin(np.abs(m.c - c_t)))
            f0, f2, c = m._rad[0][ik], m._rad[2][ik], m.c[ic]
            quad = (2.0 * f0 + f2 / s2) * (f2 / s2) / (f0 ** 2 + f2 ** 2)
            cw = quad * (1.0 - 3.0 * c ** 2)
            # the same expression with the phase dropped, i.e. f2 -> -f2
            flip = ((2.0 * f0 - f2 / s2) * (-f2 / s2) / (f0 ** 2 + f2 ** 2)
                    * (1.0 - 3.0 * c ** 2))
            got = float(tagged.azz_tensor_curve(m, ic)[ik])
            if abs(cw - flip) < _CW_GAP:               # the check must bite
                bad.append("polligen.tagged: at %s k = %.4f, cos th_k = %.4f "
                           "the phase-less form differs by only %.3g, so this "
                           "cell no longer discriminates -- move the point"
                           % (name, m.k[ik], c, abs(cw - flip)))
                continue
            if abs(got - cw) > 1e-10:
                bad.append("polligen.tagged (%s): A_zz^wf = %.9f at k = %.4f, "
                           "cos th_k = %.4f, against Cosyn-Weiss Eq. (6.12)'s "
                           "%.9f built from the model's own f0, f2 (the "
                           "phase-less f2 -> -f2 form gives %.9f). The i^L "
                           "partial-wave phase of _amp2_table, restored "
                           "2026-09-15, is the thing that makes these equal"
                           % (name, got, m.k[ik], c, cw, flip))
    return bad
