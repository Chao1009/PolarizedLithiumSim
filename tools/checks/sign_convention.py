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

Two assertions:

  * the constant is still -1 -- if it is flipped back, every "opposite"
    sentence in the reports becomes wrong at once;
  * every template that displays the master formula states the convention
    within twelve lines of the display, either by calling the simulation's
    sign the "opposite" one or by printing the literature relation
    A_zz = -(2/3) b1/F1.
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
