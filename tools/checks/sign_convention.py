"""The tensor sign convention: the constant, and the caveat at every display.

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

Three assertions:

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
    untouched.
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
