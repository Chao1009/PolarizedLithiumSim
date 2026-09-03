"""The sweet-spot kinematics, recomputed from the selector and kinematics.py.

Review section 5.2 item 7 (docs/consistency_review_2026-09-02.md): "Table 3 of
Report 2 and its kin regenerated."  The programme's twelve sweet spots -- four
per beam configuration, chosen by `money_cos2phi.pick_sweet_spots_banded` on
the figure-of-merit map and reached through `reco_chain_figures.sweet_spots` --
fix a whole family of published numbers: Report 2 Table 3's (x, Q2), y, W and
scattered-electron columns, the delta-y/y the electron method alone gives
there, the eta the detector model is quoted at in Report 2 Table 2 and Report
3 Table 8, and the angular resolution the eta table hands each spot.  Six of
the review's findings are one of those cells left at a value the selector no
longer gives (F033, F062, F089, F201, F207, F208) and one is an eta span that
was the mid configuration's four quoted for all twelve (F134).

Everything here is recomputed:

  x, Q2   the spot itself, `sweet_spots(0/1/2)[:4]`
  y       `kinematics.y_from_xq2(x, Q2, s)`
  W       `sqrt(kinematics.w2(x, Q2))`
  E', th', eta   `kinematics.scattered_electron(x, y, s, E_e)`, with th'
          measured from the ION direction (pi - theta), which is the
          convention Table 3 and Figure 4a print
  dy/y    `reco.electron_method_resolution(y, th', 0.012, 0.0,
          de_beam_over_e=1e-3)` -- the 1.2% calorimeter term of Section 4.2
          with the 1e-3 beam-energy spread under it, which is what
          reproduces the table's 1.18 / 0.46 / 1.07 / 0.46 (the 1.2% alone
          gives 1.17 in the first cell)
  dtheta  `reco.tracking_angular_resolution(eta)`

A quoted number is compared at ITS OWN precision: "14.3", "1.14" and "107" are
all correct readings of the same three quantities, so the test is that the
derived value rounds to the quoted one at the number of decimals the document
chose.  That keeps the check honest about the last digit -- which is the digit
every one of these findings turns on -- without legislating a format.

The Appendix A revision rows are cut out of every template first: they record
what a cell read before a re-derivation, which is the programme's convention
and not a stale number.

Cost: `sweet_spots` for the three configurations is 0.9 s, all of it the first
figure-of-merit projection; nothing here executes a producing script.
"""

import math
import pathlib
import re

ROOT = checker.ROOT                                            # noqa: F821

R0 = ROOT / "reports/polarized_li_primer.template.html"
R1 = ROOT / "reports/cos2phi_money_plots_report.template.html"
R2 = ROOT / "reports/reconstruction_chain_report.template.html"
R3 = ROOT / "reports/eic_epic_reference.template.html"
EVGEN_README = ROOT / "evgen/README.md"

_APPENDIX_A = re.compile(r"<h2>\s*Appendix A")
_WORDS = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6,
          "seven": 7, "eight": 8, "nine": 9, "ten": 10, "eleven": 11,
          "twelve": 12}


def _body(path):
    """A document with its Appendix A revision history removed."""
    txt = path.read_text(errors="ignore")
    m = _APPENDIX_A.search(txt)
    return txt[:m.start()] if m else txt


def _line(txt, pos):
    return txt[:pos].count("\n") + 1


def _at(path, txt, pos):
    return "%s:%d" % (path.relative_to(ROOT), _line(txt, pos))


# --- the twelve spots ------------------------------------------------------

def _spots():
    import sys
    sys.path.insert(0, str(ROOT / "evgen" / "scripts"))
    import reco_chain_figures as rcf
    from polli_fastsim.kinematics import scattered_electron, y_from_xq2, w2
    from polligen import reco

    out = []
    for ci in (0, 1, 2):
        config, _scen, _proj, _obs, spots, _kern = rcf.sweet_spots(ci)
        s = config.sqrt_s_per_nucleon ** 2
        for k, (x, q2, _i, _j) in enumerate(spots):
            y = float(y_from_xq2(x, q2, s))
            e_p, theta, eta = scattered_electron(x, y, s,
                                                 config.electron_energy)
            theta_ion = math.pi - float(theta)
            _dq2, dy, _dx = reco.electron_method_resolution(
                y, theta_ion, 0.012, 0.0, de_beam_over_e=1e-3)
            out.append({
                "config": ci, "spot": k + 1, "label": config.label(),
                "x": float(x), "q2": float(q2), "y": y,
                "w": float(w2(x, q2)) ** 0.5, "e_prime": float(e_p),
                "theta_mrad": 1e3 * theta_ion, "eta": float(eta),
                "dy_over_y": float(dy),
                "dtheta_mrad": 1e3 * float(reco.tracking_angular_resolution(eta)),
            })
    return out


SPOTS = _spots()
MID = [s for s in SPOTS if s["config"] == 1]
ETA = [s["eta"] for s in SPOTS]
ETA_SPAN = (min(ETA), max(ETA))
MID_ETA_SPAN = (min(s["eta"] for s in MID), max(s["eta"] for s in MID))


# --- comparing a quoted number with a derived one --------------------------

def _places(token):
    """Decimal places the document chose to print."""
    token = token.replace("−", "-")
    return len(token.split(".")[1]) if "." in token else 0


def _value(token):
    return float(token.replace("−", "-").replace(",", ""))


def _agrees(token, derived, ulp=0.5):
    """Does `derived` round to `token` at the precision `token` is written
    to?  Half a unit in the last place, plus an epsilon for exact ties.

    `ulp=1.0` relaxes that to a whole unit, and is used for the ENDPOINTS OF
    A SPAN only.  The reports write the mid configuration's four as
    "η = −2.9 to −1.6" off Table 3's own −2.93 … −1.65, truncating the
    endpoint towards the interior rather than rounding it away; that is the
    house style for a span and is not what F134 was about.  The span check
    that matters -- which of the derived spans a sentence is quoting -- is
    unaffected, the twelve-spot and four-spot spans being half a unit apart
    in neither endpoint."""
    return abs(_value(token) - derived) <= ulp * 10.0 ** -_places(token) + 1e-9


def _cmp(where, what, token, derived, bad):
    if not _agrees(token, derived):
        bad.append("%s %s reads %s, the selector gives %.4f"
                   % (where, what, token, derived))


# --- Report 2 Table 3 ------------------------------------------------------

_TABLE3_CAP = re.compile(r"<p class=\"tabcap\"><b>Table 3 — the four sweet "
                         r"spots")
_ROW = re.compile(r"<tr>(.*?)</tr>", re.S)
_CELL = re.compile(r"<td[^>]*>(.*?)</td>", re.S)


def _table_before(txt, cap_re):
    m = cap_re.search(txt)
    if not m:
        return None, None
    start = txt.rfind("<table>", 0, m.start())
    end = txt.find("</table>", start)
    if start < 0 or end < 0:
        return None, None
    return txt[start:end], start


@check("derived: Report 2 Table 3's sweet-spot cells are what kinematics.py gives")  # noqa: F821
def _():
    txt = _body(R2)
    table, at = _table_before(txt, _TABLE3_CAP)
    if table is None:
        return ["reconstruction_chain_report.template.html no longer carries "
                "a 'Table 3 — the four sweet spots' caption above a table"]
    rows = [_CELL.findall(r) for r in _ROW.findall(table)]
    rows = [r for r in rows if len(r) >= 5]
    if len(rows) != 4:
        return ["reconstruction_chain_report.template.html:%d Table 3 has %d "
                "data rows, the selector gives 4 sweet spots"
                % (_line(txt, at), len(rows))]
    bad = []
    where = "reconstruction_chain_report.template.html:%d" % _line(txt, at)
    for row, spot in zip(rows, MID):
        tag = "%s Table 3 row %d" % (where, spot["spot"])
        m = re.match(r"\s*([\d.]+),\s*([\d.]+)\s*$", row[0])
        if not m:
            bad.append("%s: '%s' is not an (x, Q²) pair" % (tag, row[0]))
            continue
        _cmp(tag, "x", m.group(1), spot["x"], bad)
        _cmp(tag, "Q²", m.group(2), spot["q2"], bad)
        _cmp(tag, "y", row[1].strip(), spot["y"], bad)
        _cmp(tag, "W", row[2].strip(), spot["w"], bad)
        e = re.match(r"\s*([\d.]+)\s*GeV,\s*([\d.]+)\s*mrad,\s*(−[\d.]+)\s*$",
                     row[3])
        if not e:
            bad.append("%s: '%s' is not an (E′, θ′, η) cell" % (tag, row[3]))
        else:
            _cmp(tag, "E′", e.group(1), spot["e_prime"], bad)
            _cmp(tag, "θ′", e.group(2), spot["theta_mrad"], bad)
            _cmp(tag, "η", e.group(3), spot["eta"], bad)
        _cmp(tag, "δy/y (e′ alone)", row[4].strip(), spot["dy_over_y"], bad)
    return bad


# --- the y and delta-y/y statements ----------------------------------------

@check("derived: every sweet-spot y and δy/y statement matches the selector")  # noqa: F821
def _():
    bad = []
    ys = [s["y"] for s in MID]
    dyy = [s["dy_over_y"] for s in MID]

    # the four-value y list, at whatever precision the sentence uses
    quad = re.compile(r"y = (0\.\d+), (0\.\d+), (0\.\d+) and (0\.\d+)")
    seen_quad = 0
    for path in (R0, R1, R2, R3, EVGEN_README):
        txt = _body(path)
        for m in quad.finditer(txt):
            seen_quad += 1
            for tok, val, k in zip(m.groups(), ys, (1, 2, 3, 4)):
                _cmp(_at(path, txt, m.start()),
                     "the sweet-spot y of spot %d" % k, tok, val, bad)
    if not seen_quad:
        bad.append("no report lists the four sweet-spot y values any more -- "
                   "the pattern in tools/checks/sweet_spot_kinematics.py "
                   "matches nothing")

    # the bands: y over the four spots, and delta-y/y as a fraction and as a
    # percentage.  "sweet spot" has to be nearby: y and delta-y/y spans are
    # quoted for other things too (the Sigma method, the noise scan).
    bands = (
        (re.compile(r"y = (0\.\d+)\s*[–-]\s*(0\.\d+)"), r"sweet spot",
         "the sweet-spot y band", min(ys), max(ys), 1.0),
        (re.compile(r"δy/y = (\d\.\d+)\s*[–-]\s*(\d\.\d+)"), r"e′ alone|"
         r"electron method|energy resolution|sweet spot",
         "the δy/y band of the electron method", min(dyy), max(dyy), 1.0),
        (re.compile(r"δy/y = (\d+)\s*[–-]\s*(\d+)\s*%"), r"e′ alone|sweet spot",
         "the δy/y band of the electron method", min(dyy), max(dyy), 100.0),
    )
    for rx, ctx, what, lo, hi, scale in bands:
        seen = 0
        for path in (R1, R2, EVGEN_README):
            txt = _body(path)
            for m in rx.finditer(txt):
                near = txt[max(0, m.start() - 200):m.end() + 120]
                if not re.search(ctx, near):
                    continue
                seen += 1
                tag = _at(path, txt, m.start())
                _cmp(tag, what + " (low end)", m.group(1), scale * lo, bad)
                _cmp(tag, what + " (high end)", m.group(2), scale * hi, bad)
        if not seen:
            bad.append("no report states %s any more -- the pattern in "
                       "tools/checks/sweet_spot_kinematics.py matches nothing"
                       % what)

    # how many of the four the electron method cannot reconstruct
    n_over = sum(1 for v in dyy if v > 1.0)
    for path in (R1, R2):
        txt = _body(path)
        for m in re.finditer(r"(\w+) bins? (?:with|have) δy/y &gt; 1", txt):
            got = _WORDS.get(m.group(1).lower())
            if got is not None and got != n_over:
                bad.append("%s says %s of the four sweet spots have "
                           "δy/y > 1; the selector gives %d"
                           % (_at(path, txt, m.start()), m.group(1), n_over))

    # Report 2 Table 2's angular-resolution bracket, spot by spot (F201)
    txt = _body(R2)
    m = re.search(r"\((\d) mrad at spots 1[–-]2, (\d) mrad at spots 3[–-]4\)",
                  txt)
    if not m:
        bad.append("reconstruction_chain_report.template.html Table 2 no "
                   "longer names the tracking angular resolution at the "
                   "sweet spots -- the pattern in "
                   "tools/checks/sweet_spot_kinematics.py matches nothing")
    else:
        tag = "reconstruction_chain_report.template.html:%d" % _line(txt,
                                                                     m.start())
        want = [s["dtheta_mrad"] for s in MID]
        if len(set("%.0f" % v for v in want[:2])) != 1 or \
           len(set("%.0f" % v for v in want[2:])) != 1:
            bad.append("%s Table 2 splits the tracking angular resolution as "
                       "spots 1-2 against 3-4, but the η table gives %s mrad "
                       "at the four spots"
                       % (tag, " / ".join("%.0f" % v for v in want)))
        else:
            _cmp(tag, "the tracking angular resolution at spots 1-2",
                 m.group(1), want[0], bad)
            _cmp(tag, "the tracking angular resolution at spots 3-4",
                 m.group(2), want[2], bad)
    return bad


# --- the eta spans ---------------------------------------------------------

@check("derived: every sweet-spot η span in Reports 2 and 3 is one the selector gives")  # noqa: F821
def _():
    """Every "eta = A to B" a report writes about the sweet spots must be
    either the span of the twelve spots or the span of one configuration's
    four, and both reports must carry the twelve-spot one: F134 is Report 3
    Table 8 and Report 2 Table 2 quoting the mid configuration's
    -2.9 to -1.6 as "every sweet spot of the programme"."""
    allowed = {"the twelve spots": ETA_SPAN}
    for ci, name in ((0, "the low configuration's four"),
                     (1, "the mid configuration's four"),
                     (2, "the top configuration's four")):
        e = [s["eta"] for s in SPOTS if s["config"] == ci]
        allowed[name] = (min(e), max(e))

    bad = []
    span_rx = re.compile(r"(−\d\.\d+)\s+to\s+(−\d\.\d+)")
    for path in (R2, R3):
        txt = _body(path)
        carries_twelve = False
        for m in span_rx.finditer(txt):
            near = txt[max(0, m.start() - 220):m.end() + 120]
            if not re.search(r"η|spot", near):
                continue
            lo, hi = m.group(1), m.group(2)
            match = [k for k, (a, b) in allowed.items()
                     if _agrees(lo, a, 1.0) and _agrees(hi, b, 1.0)]
            if not match:
                bad.append("%s the η span %s to %s is neither the twelve "
                           "spots' (%.2f to %.2f) nor any one "
                           "configuration's four"
                           % (_at(path, txt, m.start()), lo, hi,
                              ETA_SPAN[0], ETA_SPAN[1]))
            elif "the twelve spots" in match:
                carries_twelve = True
        if not carries_twelve:
            bad.append("%s no longer states the η span of the twelve sweet "
                       "spots (%.1f to %.1f); F134 is exactly the mid "
                       "configuration's four quoted for all of them"
                       % (path.name, ETA_SPAN[0], ETA_SPAN[1]))

    # Report 3 Table 8's caption counts the twelve by detector region
    txt = _body(R3)
    regions = ((r"(\w+) of the twelve sit in the backward endcap",
                "the backward endcap", lambda e: e < -2.0),
               (r"(\w+)[^.]{0,40}sit in the backward transition",
                "the backward transition", lambda e: -2.0 <= e < -1.0))
    for rx, what, sel in regions:
        m = re.search(rx, txt)
        if not m:
            bad.append("eic_epic_reference.template.html Table 8's caption no "
                       "longer counts the sweet spots in %s -- the pattern in "
                       "tools/checks/sweet_spot_kinematics.py matches nothing"
                       % what)
            continue
        want = sum(1 for e in ETA if sel(e))
        got = _WORDS.get(m.group(1).lower())
        if got is None:
            continue
        if got != want:
            bad.append("eic_epic_reference.template.html:%d says %s of the "
                       "twelve sweet spots sit in %s; the selector gives %d"
                       % (_line(txt, m.start()), m.group(1), what, want))

    # and it names the two transition spots and the one in the barrel
    m = re.search(r"at η = (−\d\.\d+) and (−\d\.\d+),\s*sit in the backward "
                  r"transition", txt)
    if m:
        want = sorted((e for e in ETA if -2.0 <= e < -1.0), reverse=True)
        got = [m.group(1), m.group(2)]
        if len(want) != 2:
            bad.append("eic_epic_reference.template.html:%d names two sweet "
                       "spots in the backward transition; the selector gives "
                       "%d" % (_line(txt, m.start()), len(want)))
        else:
            for tok, val in zip(sorted(got, key=_value), sorted(want)):
                _cmp("eic_epic_reference.template.html:%d"
                     % _line(txt, m.start()),
                     "the η of a backward-transition sweet spot", tok, val,
                     bad)
    m = re.search(r"spot at 5 × 40\.8, η = (−\d\.\d+)", txt)
    if m:
        want = [e for e in ETA if -1.0 <= e <= 1.0]
        if len(want) != 1:
            bad.append("eic_epic_reference.template.html:%d names one sweet "
                       "spot in the barrel; the selector gives %d"
                       % (_line(txt, m.start()), len(want)))
        else:
            _cmp("eic_epic_reference.template.html:%d" % _line(txt, m.start()),
                 "the η of the barrel sweet spot", m.group(1), want[0], bad)
    return bad
