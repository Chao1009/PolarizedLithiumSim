"""Canonical far-forward triples and range statements, re-derived from the code.

Review section 5.2 item 3 (docs/consistency_review_2026-09-02.md): "Canonical
triples from the code."  Every far-forward number the reports quote as a
per-configuration triple -- the tagging envelope, the pot insertion it asks
for, the measured silicon aperture, the transport levers, the blind block, the
alpha tag -- is DERIVED: it comes out of `farforward.tagging_optics_point`,
`farforward.POT_LEVERS`, `farforward.POT_BLIND_HALF_WIDTH`,
`farforward.THETA_RP_OUTER_MEASURED`, `reco.RP_APERTURE_MEASURED` or a
fixed-seed replay of `fastsim/scripts/tagging_acceptance.py`.  A
re-measurement therefore moves several documents at once, and the review found
five findings (F051, F078, F113, F215, F238) that are exactly one of those
triples left behind in one document while the others moved.

This module evaluates each quantity at import and then scans the documents for
the STATEMENT of it, comparing digit by digit.  A statement is recognised by

  * its SHAPE -- a triple (or span, or scalar) of the right precision with the
    right unit or bracket after it, and no fourth number continuing the list;
  * optionally a CONTEXT word within a short window either side of it; and
  * PROXIMITY -- every captured number must lie within `tol` of the derived
    one, and a span must also be as WIDE as the derived span.

The proximity band is what makes the scan specific.  These reports carry
hundreds of numeric triples, several sharing a shape with a canonical one (the
Table 3 acceptance column 0.45 / 0.25 / 0.10, the Sigma-method resolutions
0.32 / 0.22 / 0.29, the September-2024 aperture 2.00 / 1.35 / 1.03, the 5 mrad
convention's 96 / 106 / 150 mm), and a shape-only scan flags all of them.  The
band is set per statement so that the RETIRED form is inside it -- 0.33 / 0.17
/ 0.12 mrad for the envelope, 27-35% for the alpha tag, 0.12-0.33 mrad for its
span, 3.6 mm for the 18 x 275 insertion, "-" for the 5 x 41 vertical lever --
while the neighbouring quantities are outside.  A number that has moved by
more than the band is a different quantity, not a stale digit; catching that
is the retired-strings list's business, not this one's.

The width rule on spans is the same idea one dimension up: Report 1 Section 8
quotes "a 33-37% tag" for the two configurations it proposes, which is a
statement about a subset and not a stale form of the three-configuration
25-37%.

History is not flagged: the Appendix A revision rows of the templates and the
older `## Development run N (date)` sections of the plans are the programme's
record of what the numbers USED to be -- plans/00 still carries the run-13
envelope 0.33 / 0.17 / 0.12 mrad and the run-8 tagging optics beside it -- and
both are cut out of the scanned text.

Not re-derived here because a built-in check already owns it: the alpha + d
median separation 17.3 / 10.7 / 10.9 mm ("drift: the alpha+d separation table
matches the current energies" in tools/consistency_check.py).
"""

import glob
import pathlib
import re

import numpy as np

ROOT = checker.ROOT                                            # noqa: F821
CFG = ("5x41", "10x100", "18x275")

# --- the documents, and the history inside them ----------------------------

DOCS = ([pathlib.Path(p) for p in
         sorted(glob.glob(str(ROOT / "reports/*.template.html")))]
        + [ROOT / "reports/index.html", ROOT / "README.md",
           ROOT / "evgen/README.md", ROOT / "fastsim/README.md",
           ROOT / "refs/README.md", ROOT / "docs/reproduction_manual.md"]
        + [pathlib.Path(p) for p in
           sorted(glob.glob(str(ROOT / "plans/*.md")))])

_APPENDIX_A = re.compile(r"<h2>\s*Appendix A")
_RUN_SECTION = re.compile(r"^## Development run [^\n(]*\((\d{4}-\d\d-\d\d)",
                          re.M)


def _live_spans(txt):
    """The [start, end) ranges of `txt` that state the CURRENT numbers."""
    end = len(txt)
    m = _APPENDIX_A.search(txt)
    if m:
        end = m.start()
    runs = [(mm.start(), mm.group(1)) for mm in _RUN_SECTION.finditer(txt[:end])]
    if not runs:
        return [(0, end)]
    newest = max(date for _s, date in runs)
    bounds = [s for s, _d in runs] + [end]
    spans, cursor = [], 0
    for i, (start, date) in enumerate(runs):
        if date == newest:
            continue
        spans.append((cursor, start))
        cursor = bounds[i + 1]
    spans.append((cursor, end))
    return [(a, b) for a, b in spans if b > a]


# --- the derived quantities ------------------------------------------------

def _derive():
    from polli_fastsim import beams
    from polli_fastsim import farforward as ff
    from polli_fastsim import spectator as sp
    from polligen import reco

    cfgs = beams.default_configs("6Li")
    pt = [ff.tagging_optics_point(c) for c in cfgs]
    d = {
        "env_x": [1e3 * p["env_x"] for p in pt],
        "env_y": [1e3 * p["env_y"] for p in pt],
        "eps": [p["acceptance"] for p in pt],
        "r12": [ff.POT_LEVERS[k][0] for k in CFG],
        "r34": [ff.POT_LEVERS[k][1] for k in CFG],
        "disp": [ff.POT_LEVERS[k][2] for k in CFG],
        "blind": [1e3 * ff.POT_BLIND_HALF_WIDTH[k] for k in CFG],
        "ap_h": [1e3 * reco.RP_APERTURE_MEASURED[k][0] for k in CFG],
        "ap_v": [1e3 * reco.RP_APERTURE_MEASURED[k][1] for k in CFG],
        "outer": [1e3 * ff.THETA_RP_OUTER_MEASURED[k] for k in CFG],
    }
    d["env_mm"] = [e * r for e, r in zip(d["env_x"], d["r12"])]
    # the scale a recorded alpha + d pair's separation sits at: three times
    # the alpha's own pot-plane displacement at the envelope, on the lever of
    # the tighter axis (`nearbeam_two_hit._merge_scale`)
    d["merge_mm"] = [3.0 * min(r12 * ex, r34 * ey) for r12, r34, ex, ey
                     in zip(d["r12"], d["r34"], d["env_x"], d["env_y"])]
    d["closer"] = [a / e for a, e in zip(d["ap_h"], d["env_x"])]
    d["ratio"] = [a / b for a, b in zip(d["r12"], d["r34"])]

    # The 6Li alpha tag at the tagging optics, as Report 3 Table 6 publishes
    # it.  `tagging_acceptance.py` draws 400k spectators per (channel,
    # configuration, beta) from ONE default_rng(7) and LI6_ALPHA_TAG is its
    # first channel, so replaying that channel's nine draws in order
    # reproduces the published 0.3145 / 0.2229 / 0.2913 exactly, in 1.1 s.
    # Reproducing the draw and not merely the physics is what makes the check
    # safe: 0.3145 sits on the 31/32 rounding boundary that F051 turns on, and
    # a fresh sample of the same size lands either side of it (the 8e6-event
    # converged value is 0.3147).
    rng = np.random.default_rng(7)
    tag = {}
    for c in cfgs:
        for beta in (0.20, 0.30, 0.40):
            kin = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG,
                                              c.ion_momentum_per_nucleon,
                                              400_000, beta=beta, rng=rng)
            if beta == 0.30:
                acc = ff.acceptance_summary(
                    kin["R"], kin["theta"], kin["pT"], ff.tagging_optics(c),
                    phi=kin["phi"], pot_config=ff.yr_config_key(c))
                tag[ff.yr_config_key(c)] = 1.0 - acc["lost"]
    d["tag"] = [tag[k] for k in CFG]
    return d


def _veto():
    """The partner-deuteron veto given an alpha that fakes a coherent tag, at
    the tagging optics and at the measured outer pot edge:
    `nearbeam_two_hit.py --veto-events 12000000` replayed (one
    default_rng(7) per configuration, in the 1e6 chunks the script draws).
    9 s, hence --full only; 1e6 events is not enough, the 10 x 100 entry
    sitting at 82.5% where the published range turns on the second digit."""
    from polli_fastsim import beams
    from polli_fastsim import farforward as ff
    from polli_fastsim import spectator as sp
    at_5mrad, at_edge = [], []
    for c in beams.default_configs("6Li"):
        key = ff.yr_config_key(c)
        opt = ff.tagging_optics(c)
        rng = np.random.default_rng(7)
        n = {"f": 0, "v": 0, "fe": 0, "ve": 0}
        for _chunk in range(12):
            ev = sp.breakup_lab_kinematics(sp.LI6_ALPHA_TAG,
                                           c.ion_momentum_per_nucleon,
                                           1_000_000, beta=0.30, rng=rng)
            a, p = ev["spectator"], ev["partner"]
            for th, kf, kv in ((None, "f", "v"),
                               (ff.THETA_RP_OUTER_MEASURED[key], "fe", "ve")):
                ra = ff.route_charged(a["R"], a["theta"], a["pT"], opt,
                                      phi=a["phi"], theta_outer=th,
                                      pot_config=key)
                rp = ff.route_charged(p["R"], p["theta"], p["pT"], opt,
                                      phi=p["phi"], theta_outer=th,
                                      pot_config=key)
                fake, seen = (ra == 4), ((rp == 1) | (rp == 4))
                n[kf] += int(fake.sum())
                n[kv] += int((fake & seen).sum())
        at_5mrad.append(n["v"] / max(n["f"], 1))
        at_edge.append(n["ve"] / max(n["fe"], 1))
    return at_5mrad, at_edge


D = _derive()


# --- the statement table ---------------------------------------------------

#: how far either side of a match a statement's context word may sit
CONTEXT_WINDOW = 70

SEP = r"\s*(?:,|/|·|,?\s*and)\s*"
NO_FOURTH = r"(?!\s*(?:,|/|·|\s+and)\s*\d)"
D1 = r"(\d+\.\d)"
D2 = r"(\d+\.\d\d)"
D3 = r"(\d+\.\d\d\d)"
D12 = r"(\d+\.\d\d?)"
INT = r"(\d+)"
SHUT = r"(shut|\(shut\)|—|\d+\.\d\d)"


def _trip(num, tail=True):
    return num + SEP + num + SEP + num + (NO_FOURTH if tail else "")


def _f(vals, spec):
    return [spec % v for v in vals]


def _span(vals, spec):
    return [spec % min(vals), spec % max(vals)]


def _interleave(*columns):
    return [x for row in zip(*columns) for x in row]


# (name, regex, expected tokens, tolerance, context regex or None, is_span)
STATEMENTS = [
    ("the tagging envelope in angle (tagging_optics_point env_x)",
     _trip(D2) + r"(?=\s*(?:mrad|\)))", _f(D["env_x"], "%.2f"), 0.20,
     r"envelope|follow(?:ing|s|\s+the)", False),

    ("the tagging envelope quoted as a span",
     r"(\d+\.\d\d)\s*[–-]\s*(\d+\.\d\d)\s*mrad", _span(D["env_x"], "%.2f"),
     0.35, r"envelope", True),

    # one OR two decimals here: F078 is this triple's 18 x 275 entry written
    # "3.6 mm" off the rounded 0.12 mrad instead of 3.52 off the envelope
    ("the pot insertion the tagging optics asks for (env_x x R12)",
     _trip(D12) + r"\s*mm", _f(D["env_mm"], "%.2f"), 0.20, None, False),

    ("the 18 x 275 tagging envelope at the pot plane, as a scalar",
     r"tagging envelope of (\d+\.\d+) mm at 18",
     ["%.2f" % D["env_mm"][2]], 0.20, None, False),

    ("how much closer the tagging optics asks the pots to sit (ap_h/env_x)",
     r"factor\s+" + _trip(D1) + r"\s*(?:×\s*)?closer",
     _f(D["closer"], "%.1f"), 0.30, None, False),

    # the same ratio stated the other way round -- how far outside the tagging
    # envelope the measured silicon sits (Report 4 §2, plans/10); it was the
    # phrasing that was right while the "factor ... closer" one was stale (F001)
    ("how far outside the tagging envelope the silicon sits (ap_h/env_x)",
     _trip(D1) + r"\s*(?:times\s+outside|envelope\s+widths?\s+out)",
     _f(D["closer"], "%.1f"), 0.30, None, False),

    ("the horizontal-to-vertical lever ratio (R12/R34)",
     r"factor\s+" + _trip(D1), _f(D["ratio"], "%.1f"), 0.30, None, False),

    ("the vertical pot lever R34 (POT_LEVERS)",
     r"R₃₄\s*(?:of|=)\s*" + _trip(SHUT) + r"\s*m\b",
     _f(D["r34"], "%.2f"), 0.20, None, False),

    ("the measured pot transport (R12, R34, D) per configuration",
     SEP.join([D2, D2, D3, D2, D2, D3, D2, D2, D3]),
     None, 0.20, None, False),

    ("the alpha + d merge scale (3 x the envelope on the tighter lever)",
     INT + SEP + INT + SEP + D1 + NO_FOURTH + r"\s*mm",
     _f(D["merge_mm"][:2], "%.0f") + _f(D["merge_mm"][2:], "%.1f"), 0.25,
     r"scale", False),

    ("the per-configuration blind block (POT_BLIND_HALF_WIDTH)",
     _trip(INT) + r"\s*mm", _f(D["blind"], "%.0f"), 0.20, None, False),

    ("the measured horizontal pot aperture (RP_APERTURE_MEASURED)",
     _trip(D2) + r"(?=\s*(?:mrad|\)))", _f(D["ap_h"], "%.2f"), 0.20, None,
     False),

    ("the measured vertical pot aperture, 5 x 41 shut",
     SHUT + SEP + D2 + SEP + D2 + NO_FOURTH + r"\s*mrad",
     [("shut", "%.2f" % D["ap_v"][0])] + _f(D["ap_v"][1:], "%.2f"),
     0.20, None, False),

    ("the measured pot aperture written as horizontal/vertical pairs",
     D2 + r"/" + SHUT + r",\s*" + D2 + r"/" + D2 + r",\s*" + D2 + r"/" + D2
     + r"\s*mrad",
     ["%.2f" % D["ap_h"][0], ("shut", "%.2f" % D["ap_v"][0])]
     + _interleave(_f(D["ap_h"][1:], "%.2f"), _f(D["ap_v"][1:], "%.2f")),
     0.20, None, False),

    ("the measured outer pot edge (THETA_RP_OUTER_MEASURED)",
     _trip(D2) + r"\s*mrad", _f(D["outer"], "%.2f"), 0.20, None, False),

    ("the coherent tagged fraction at the tagging optics",
     _trip(D2), _f(D["eps"], "%.2f"), 0.12,
     r"ε\s*=|tagged fraction|tagging|acceptance", False),

    ("the coherent tag quoted as a span",
     r"(\d\d)\s*[–-]\s*(\d\d)\s*%\s*(?:coherent|tag\b)",
     _span([100 * v for v in D["eps"]], "%.0f"), 0.35, None, True),

    ("the ⁶Li α tag quoted as a span (tagging_acceptance.py)",
     r"(\d\d)\s*[–-]\s*(\d\d)\s*%\s*(?:⁶Li\s*)?α tag",
     _span([100 * v for v in D["tag"]], "%.0f"), 0.35, None, True),
]

# the (R12, R34, D) statement interleaves three configurations, so its
# expectation is filled in here rather than in the table above; find it by
# name -- an index would silently move the moment a statement is inserted
_TRANSPORT = next(i for i, s in enumerate(STATEMENTS)
                  if s[0].startswith("the measured pot transport"))
STATEMENTS[_TRANSPORT] = STATEMENTS[_TRANSPORT][:2] + (
    _interleave(_f(D["r12"], "%.2f"), _f(D["r34"], "%.2f"),
                _f(D["disp"], "%.3f")),) + STATEMENTS[_TRANSPORT][3:]

# no statement may reach the scan without an expectation (a None here means a
# table entry was added without its derived value, and the scan would raise)
assert all(s[2] is not None for s in STATEMENTS), \
    "a canonical statement has no expectation: %s" % [
        s[0] for s in STATEMENTS if s[2] is None]


def _numeric(token):
    m = re.fullmatch(r"\d+(?:\.\d+)?", token.strip("()"))
    return float(m.group()) if m else None


def _spellings(want):
    """An expectation is one accepted string, or a tuple of them: the
    5 x 41 vertical aperture is written "shut", "(shut)" or as its 6.49 mrad
    value, and all three are the same statement."""
    return (want,) if isinstance(want, str) else tuple(want)


def _near(token, want, tol):
    """Is `token` an instance of this statement -- the same QUANTITY as
    `want`, to within `tol`?

    A token that is not a number at all ("—", "shut") always counts as an
    instance, so that a placeholder left where a measured value belongs is
    REPORTED rather than skipped; that is F215, the 5 x 41 vertical lever
    still written "—" after it was measured at 4.56 m."""
    spell = _spellings(want)
    tok = token.strip("()")
    if tok in {s.strip("()") for s in spell}:
        return True
    got = _numeric(tok)
    if got is None:
        return True
    return any(abs(got - w) <= tol * max(abs(w), 1e-12)
               for w in (_numeric(s) for s in spell) if w is not None)


def _accepts(token, want):
    return token.strip("()") in {s.strip("()") for s in _spellings(want)}


def _same_width(got, expect, tol):
    """A quoted span must be as WIDE as the derived one, not merely have both
    ends near it: a narrower span is a statement about a named subset of the
    three configurations, not a stale three-configuration span."""
    lo, hi = (_numeric(g) for g in got)
    wlo, whi = (_numeric(_spellings(w)[0]) for w in expect)
    return abs((hi - lo) - (whi - wlo)) <= tol * max(abs(whi - wlo), 1e-12)


def _scan(name, pattern, expect, tol, context, is_span=False, window=None):
    bad, seen = [], 0
    rx = re.compile(pattern)
    for path in DOCS:
        if not path.exists():
            continue
        txt = path.read_text(errors="ignore")
        spans = _live_spans(txt)
        for m in rx.finditer(txt):
            if not any(a <= m.start() < b for a, b in spans):
                continue
            got = list(m.groups())
            if len(got) != len(expect):
                continue
            if not all(_near(g, w, tol) for g, w in zip(got, expect)):
                continue
            if is_span and not _same_width(got, expect, tol):
                continue
            w = CONTEXT_WINDOW if window is None else window
            if context and not re.search(
                    context, txt[max(0, m.start() - w):m.end() + w]):
                continue
            seen += 1
            clean = [g.strip("()") for g in got]
            if not all(_accepts(g, w) for g, w in zip(got, expect)):
                bad.append("%s:%d %s reads %s, the code gives %s (\"%s\")"
                           % (path.relative_to(ROOT),
                              txt[:m.start()].count("\n") + 1, name,
                              " / ".join(clean),
                              " / ".join(_spellings(w)[0] for w in expect),
                              txt[m.start():m.end()].replace("\n", " ")[:70]))
    if not seen:
        # a pattern that matches nothing is a check gone blind, which is how
        # a triple drifts unnoticed after a rewording; re-word the pattern
        # here, or restore the statement, but do not leave it silent
        bad.append("no live document states %s any more -- the pattern in "
                   "tools/checks/canonical_numbers.py matches nothing" % name)
    return bad


@check("derived: every canonical far-forward triple matches the code")  # noqa: F821
def _():
    bad = []
    for name, pattern, expect, tol, context, is_span in STATEMENTS:
        bad += _scan(name, pattern, expect, tol, context, is_span)
    return bad


@check("derived: the α + d partner veto matches nearbeam_two_hit.py")  # noqa: F821
def _():
    if not checker.FULL:                                       # noqa: F821
        return []
    at_5mrad, at_edge = _veto()
    return (_scan("the partner veto quoted as a span (5 mrad convention)",
                  r"(\d\d)\s*[–-]\s*(\d\d)\s*%",
                  _span([100 * v for v in at_5mrad], "%.0f"), 0.35,
                  r"veto|partner d", True)
            + _scan("the partner veto per configuration (tagging optics)",
                    _trip(D3), _f(at_5mrad, "%.3f"), 0.05,
                    r"veto|partner d is recorded", window=300)
            # the edge triple sits 0.03-0.05 from the 5 mrad one and from the
            # beta = 0.20 / 0.40 scan the manual prints beside it, so it is
            # recognised only where the sentence names the edge, and the
            # window reaches back past the 2.85 / 3.85 / 4.00 mrad clause
            # that always introduces it
            + _scan("the partner veto at the measured outer pot edge",
                    _trip(D2), _f(at_edge, "%.2f"), 0.12,
                    r"(?:MEASURED |measured )?outer edge|those edges",
                    window=200))
