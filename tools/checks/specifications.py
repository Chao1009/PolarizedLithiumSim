"""The specifications the programme sets, across Report 2 Section 6 and
Report 3 Table 9.

Review section 5.2 item 10 (docs/consistency_review_2026-09-02.md):
"Specification cross-check."  Three of the numbers in these reports are not
measurements or assumptions but REQUIREMENTS the programme puts on the machine
and the detector, because no EIC document addresses a tensor-polarized beam:
the tensor-polarization scale, the relative luminosity between the tensor
states, and the stability of the detector's cos 2phi' acceptance harmonic
between them.  They are derived in Report 2 Section 6 and collected in Report
3 Table 9, and the review found the two documents disagreeing about which
tolerance belongs to which (F008, F130: Table 9 row 12 carried the 10^-4 of
the acceptance harmonic as though it were the relative-luminosity ask, and
Report 0 repeated it), and the front page quoting a stale row count (F214).

Three things are checked, none of which needs the code to run:

  1. Report 2 Section 6's closing paragraph is parsed for how many
     specifications it says are the programme's to set -- the "Two of the
     resulting specifications ... The third is ..." construction -- and that
     count must equal the number of specifications this module knows about.
     A fourth specification added to the paragraph therefore fails here until
     it is added to SPECS below and to Table 9.
  2. Each specification's TOLERANCE is read out of Report 2 Section 6 and out
     of the Report 3 Table 9 row that carries it, and the two must agree.
     The row must also be one that says "ours to specify" or "ours to
     propose": a specification recorded as somebody else's ask is the F130
     confusion in the other direction.
  3. Report 3 Table 9's row count must equal the count the front page and the
     repository documents quote in words, and the caption's own enumeration
     of the rows ("Items 1-10 and 14-15 ...; 11, 12 and 16 ...; 13 ...") must
     cover every row exactly once.

Report 0 and Report 1 also state these tolerances in passing; they are left to
the retired-strings list, which is where a stale 10^-4 in a sentence with no
table row belongs.
"""

import re

ROOT = checker.ROOT                                            # noqa: F821

R2 = ROOT / "reports/reconstruction_chain_report.template.html"
R3 = ROOT / "reports/eic_epic_reference.template.html"
COUNT_DOCS = [ROOT / "reports/index.html", ROOT / "README.md",
              ROOT / "plans/00_README.md"]

_APPENDIX_A = re.compile(r"<h2>\s*Appendix A")
_WORDS = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6,
          "seven": 7, "eight": 8, "nine": 9, "ten": 10, "eleven": 11,
          "twelve": 12, "thirteen": 13, "fourteen": 14, "fifteen": 15,
          "sixteen": 16, "seventeen": 17, "eighteen": 18, "nineteen": 19,
          "twenty": 20}
_ORDINALS = {"first": 1, "second": 2, "third": 3, "fourth": 4, "fifth": 5,
             "sixth": 6, "seventh": 7, "eighth": 8}


def _body(path):
    txt = path.read_text(errors="ignore")
    m = _APPENDIX_A.search(txt)
    return txt[:m.start()] if m else txt


def _line(txt, pos):
    return txt[:pos].count("\n") + 1


# The three specifications, each as (subject, Report 2 Section 6 pattern,
# Report 3 Table 9 pattern), both patterns capturing the tolerance.  A
# superscript-minus power of ten and a percentage are compared as the strings
# the documents write, which is the form the reader compares them in.
SPECS = (
    ("the tensor-polarization scale",
     r"tensor-polarization scale,\s*δP<sub>zz</sub>/P<sub>zz</sub>\s*"
     r"≤\s*(\d+\s*%)",
     r"δP<sub>zz</sub>/P<sub>zz</sub>\s*≤\s*(\d+\s*%)"),
    ("the relative luminosity between the tensor states",
     r"relative luminosity between the tensor states,\s*≤\s*"
     r"(10[⁻−-][⁰¹²³⁴-⁹]+)",
     r"ours to specify:\s*≤\s*"
     r"(10[⁻−-][⁰¹²³⁴-⁹]+) on the two-state "
     r"ratio"),
    ("the cos 2φ′ acceptance-harmonic stability between the states",
     r"acceptance-harmonic stability of\s*"
     r"(10[⁻−-][⁰¹²³⁴-⁹]+)",
     r"acceptance harmonic stable to\s*"
     r"(10[⁻−-][⁰¹²³⁴-⁹]+)"),
)

_SEC6 = re.compile(r"<h2>6 · (.*?)<h2>7 · ", re.S)
_SET_PARA = re.compile(r"<p>((?:(?!</p>).)*this programme's to set"
                       r"(?:(?!</p>).)*)</p>", re.S)
_TABLE9_CAP = re.compile(r"<p class=\"tabcap\"><b>Table 9 — ")
_ROW = re.compile(r"<tr>(.*?)</tr>", re.S)
_CELL = re.compile(r"<td[^>]*>(.*?)</td>", re.S)


def _table9():
    """{row number: row text} of Report 3 Table 9, and the caption."""
    txt = _body(R3)
    m = _TABLE9_CAP.search(txt)
    if not m:
        return None, None, txt
    start = txt.rfind("<table>", 0, m.start())
    end = txt.find("</table>", start)
    cap_end = txt.find("</p>", m.start())
    rows = {}
    for r in _ROW.findall(txt[start:end]):
        cells = _CELL.findall(r)
        if len(cells) >= 2 and cells[0].strip().isdigit():
            rows[int(cells[0].strip())] = " ".join(cells[1:])
    return rows, txt[m.end():cap_end], txt


@check("specifications: Report 2 §6 and Report 3 Table 9 set the same tolerances")  # noqa: F821
def _():
    bad = []
    r2 = _body(R2)
    sec6 = _SEC6.search(r2)
    if not sec6:
        return ["reconstruction_chain_report.template.html no longer has a "
                "section 6 followed by a section 7"]
    para = _SET_PARA.search(sec6.group(1))
    if not para:
        return ["reconstruction_chain_report.template.html §6 no longer "
                "carries the paragraph saying which specifications are "
                "\"this programme's to set\""]
    text = para.group(1)
    at = "reconstruction_chain_report.template.html:%d" % _line(
        r2, sec6.start(1) + para.start(1))

    # 1. how many the paragraph says there are
    lead = re.search(r"(\w+) of the resulting specifications", text)
    declared = _WORDS.get(lead.group(1).lower(), 0) if lead else 0
    for m in re.finditer(r"The (\w+) is\b", text):
        declared = max(declared, _ORDINALS.get(m.group(1).lower(), 0))
    if declared != len(SPECS):
        bad.append("%s §6 says %d specifications are the programme's to set; "
                   "tools/checks/specifications.py knows %d (%s) -- add the "
                   "new one to SPECS and to Report 3 Table 9"
                   % (at, declared, len(SPECS),
                      ", ".join(s[0] for s in SPECS)))

    rows, caption, r3 = _table9()
    if rows is None:
        return bad + ["eic_epic_reference.template.html no longer carries a "
                      "'Table 9 — ' caption above a table"]

    # 2. each specification, in both documents, with the same tolerance
    for subject, r2_pat, r3_pat in SPECS:
        m2 = re.search(r2_pat, text)
        if not m2:
            bad.append("%s §6 no longer states %s -- the pattern in "
                       "tools/checks/specifications.py matches nothing"
                       % (at, subject))
            continue
        hits = [(n, re.search(r3_pat, body)) for n, body in sorted(rows.items())]
        hits = [(n, mm) for n, mm in hits if mm]
        if not hits:
            bad.append("eic_epic_reference.template.html Table 9 has no row "
                       "for %s, which Report 2 §6 sets at %s"
                       % (subject, m2.group(1)))
            continue
        if len(hits) > 1:
            bad.append("eic_epic_reference.template.html Table 9 states %s in "
                       "rows %s; it belongs in one"
                       % (subject, " and ".join(str(n) for n, _ in hits)))
        for n, m3 in hits:
            want = re.sub(r"\s+", "", m2.group(1)).replace("−", "⁻")
            got = re.sub(r"\s+", "", m3.group(1)).replace("−", "⁻")
            if got != want:
                bad.append("eic_epic_reference.template.html Table 9 row %d "
                           "sets %s at %s; Report 2 §6 sets it at %s"
                           % (n, subject, m3.group(1), m2.group(1)))
            if not re.search(r"ours to (?:specify|propose)", rows[n]):
                bad.append("eic_epic_reference.template.html Table 9 row %d "
                           "carries %s but does not say \"ours to specify\"; "
                           "Report 2 §6 lists it among the specifications "
                           "this programme sets" % (n, subject))

    # 3. the caption's own enumeration must cover every row exactly once
    body = re.sub(r"^.*?</b>", "", caption or "", flags=re.S)
    listed = []
    for m in re.finditer(r"(\d+)\s*[–-]\s*(\d+)|(\d+)", body):
        if m.group(3):
            listed.append(int(m.group(3)))
        else:
            listed += list(range(int(m.group(1)), int(m.group(2)) + 1))
    if listed and sorted(listed) != sorted(rows):
        missing = sorted(set(rows) - set(listed))
        extra = sorted(n for n in set(listed) if n not in rows)
        dupes = sorted({n for n in listed if listed.count(n) > 1})
        bad.append("eic_epic_reference.template.html Table 9's caption "
                   "enumerates %s of its %d rows%s%s%s"
                   % (len(set(listed) & set(rows)), len(rows),
                      "; unlisted: %s" % missing if missing else "",
                      "; not a row: %s" % extra if extra else "",
                      "; listed twice: %s" % dupes if dupes else ""))

    # the rows that say "ours to specify" must be exactly the ones the
    # caption calls the programme's
    ours = {n for n, b in rows.items()
            if re.search(r"ours to (?:specify|propose)", b)}
    m = re.search(r"([\d,\s]+(?:and\s*\d+)?)\s*are[^;.]*this programme's to "
                  r"specify", body)
    if m:
        called = {int(x) for x in re.findall(r"\d+", m.group(1))}
        if called != ours:
            bad.append("eic_epic_reference.template.html Table 9's caption "
                       "calls rows %s the programme's to specify; the rows "
                       "that say so are %s"
                       % (sorted(called), sorted(ours)))
    return bad


@check("specifications: the quoted count of Table 9's assumptions is its row count")  # noqa: F821
def _():
    """F214: the front page said "fifteen assumptions" through the run that
    added item 16.  Only count words of four or more are read, so that the
    "Two assumptions carry these numbers" of Report 1 §6.1 and Report 3 §4.2
    -- a different sentence about a different thing -- is not mistaken for a
    count of the table."""
    rows, _caption, _r3 = _table9()
    if rows is None:
        return ["eic_epic_reference.template.html no longer carries a "
                "'Table 9 — ' caption above a table"]
    bad, seen = [], 0
    docs = COUNT_DOCS + sorted((ROOT / "reports").glob("*.template.html"))
    for path in docs:
        if not path.exists():
            continue
        txt = _body(path)
        for m in re.finditer(r"\b(\w+)\s+assumptions\b", txt):
            word = m.group(1).lower()
            n = _WORDS.get(word, int(word) if word.isdigit() else None)
            if n is None or n < 4:
                continue
            seen += 1
            if n != len(rows):
                bad.append("%s:%d says \"%s assumptions\"; Report 3 Table 9 "
                           "has %d rows"
                           % (path.relative_to(ROOT), _line(txt, m.start()),
                              m.group(1), len(rows)))
    if not seen:
        bad.append("no document states the number of Report 3 Table 9's "
                   "assumptions any more -- reports/index.html carried it "
                   "(F214) and the pattern in "
                   "tools/checks/specifications.py matches nothing")
    return bad
