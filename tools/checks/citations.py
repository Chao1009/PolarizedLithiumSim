"""Citation and cross-pointer integrity of the five reports.

Motivated by docs/consistency_review_2026-09-02.md section 5.2 item 5
("Citation integrity"), which asks for the mechanical half of the scan the
semantic review had to do by hand:

  * expand "[a-b]" and "[a,b]" markers and assert every one of them resolves
    to an entry of the report's own reference list, and that every entry is
    cited at least once (F247);
  * assert every "Report N SS x.y" / "Report N Table k" / "Report N Figure k"
    pointer resolves to a heading, table or figure that exists in the target
    report (F028, F053, F057, F067, F081, F090, F154, F227);
  * assert every "plans/NN #k" pointer resolves to an item numbered k in that
    plans file (F177).

Scoping, so that the programme's own conventions are not flagged:

  * The reference list of each report is the <p>[n] ...</p> block inside
    <div class="refs">; markers inside that block are entry labels, not
    citations, so the block is excluded from the marker scan.
  * Appendix A revision rows are history.  A row dated before the newest one
    records what a past run said, and the review's convention is that such a
    row is not corrected when a section is renumbered or a reference list
    re-ordered, so markers and pointers inside those older rows are not
    resolved.  The newest row is the current revision note and is checked
    like the body.
  * "Every entry is cited" is judged on the body proper -- everything before
    the Appendix A heading.  A reference that only a revision row names is
    an orphan of the paper, which is exactly what F247 was.
  * A bracketed token counts as a citation marker only when it holds nothing
    but digits, commas and dashes, so "[10^-4]" and "[Erratum: ...]" are not
    markers.  "Table 10.1" and its kin are never read as "Table 10" or
    "Table 1": a table number may not be followed by a digit or a decimal.
"""

import glob
import pathlib
import re

ROOT = checker.ROOT
REPORTS = sorted(pathlib.Path(p) for p in
                 glob.glob(str(ROOT / "reports/*.template.html")))

# --- markup of the templates ------------------------------------------------

REFS_OPEN = re.compile(r'<div\s+class="refs')
REFS_ENTRY = re.compile(r'<p>\s*\[(\d+)\]')
APPENDIX_A = re.compile(r'<h2>\s*Appendix\s+A\b', re.I)
REV_ROW = re.compile(r'<td\s+class="mono">\s*(\d{4}-\d\d-\d\d)\s*</td>')
HEADING = re.compile(r'<h[234]>\s*([0-9]+(?:\.[0-9]+)*)\s*(?:·|&middot;)')
CAPTION = re.compile(r'<b>\s*(Table|Figure)\s+(\d+)\b')

# a citation marker: only digits, commas and dashes between the brackets
MARKER = re.compile(r'\[\s*(\d{1,3}(?:\s*[,–—-]\s*\d{1,3})*)\s*\]')

# "Report 3 [26] Table 6", "Report 0, Table 3", "Report 2 SS 5.2".  A full
# stop is deliberately not allowed between the two halves: "... names
# Report 4. SS 6.1 and Figure 4's caption ..." is a sentence boundary
# followed by a pointer into the reporting document itself.
POINTER = re.compile(
    r'Report\s+([0-4])\b'
    r'(?:\s*\[[0-9,\s–—-]+\])?'
    r'[\s,]*'
    r'(?:§\s*([0-9]+(?:\.[0-9]+)*)'
    r'|(Table|Figure|Fig\.)\s*([0-9]+)\b)')
# "... Report 3 SS 5 and Table 6", "... Report 3 SS 7, Table 9": the second
# half belongs to the same report.  Only a bare ", " or " and " continues a
# pointer; anything else (a closing bracket, a full stop) ends it.
CONTINUATION = re.compile(r'(?:,\s+|\s+and\s+)(Table|Figure|Fig\.)\s*([0-9]+)\b')

PLANS_ITEM = re.compile(r'plans/(\d\d)\s*#(\d+)\b')


def _lines(path):
    return path.read_text(encoding="utf-8").splitlines()


def _refs_span(lines):
    """(first, last) line indices of the reference-list <div>, or None."""
    start = None
    for i, line in enumerate(lines):
        if start is None:
            if REFS_OPEN.search(line):
                start = i
        elif "</div>" in line:
            return start, i
    return None


def _reference_entries(lines):
    """{entry number: line number} of the report's reference list."""
    span = _refs_span(lines)
    if span is None:
        return {}, None
    lo, hi = span
    out = {}
    for i in range(lo, hi + 1):
        m = REFS_ENTRY.search(lines[i])
        if m:
            out.setdefault(int(m.group(1)), i + 1)
    return out, span


def _appendix_start(lines):
    for i, line in enumerate(lines):
        if APPENDIX_A.search(line):
            return i
    return len(lines)


def _stale_appendix_lines(lines):
    """Line indices of Appendix A revision rows older than the newest one.

    Those rows are the programme's history: they are not restated when a
    section is renumbered or a reference list re-ordered, so nothing in them
    is resolved.  The newest row is the current revision note and is left in.
    """
    apx = _appendix_start(lines)
    dated = [(REV_ROW.search(lines[i]).group(1), i)
             for i in range(apx, len(lines)) if REV_ROW.search(lines[i])]
    if not dated:
        return set()
    newest = max(d for d, _ in dated)
    return {i for d, i in dated if d != newest}


def _expand(token):
    """"35-37" -> [35, 36, 37]; "3,24" -> [3, 24]."""
    out = []
    for part in re.split(r"\s*,\s*", token):
        m = re.match(r"^(\d+)\s*[–—-]\s*(\d+)$", part)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            out.extend(range(a, b + 1) if a <= b else [a, b])
        else:
            out.append(int(part))
    return out


# --- report numbering -------------------------------------------------------

def _report_numbers():
    """{report number: template stem}, from build_report.py's own map."""
    txt = (ROOT / "reports/build_report.py").read_text(encoding="utf-8")
    return {int(m.group(1)): m.group(2) for m in re.finditer(
        r"^\s*(\d)\s+([A-Za-z0-9_]+)\.template\.html", txt, re.M)}


def _index(stem):
    """(section numbers, table numbers, figure numbers) of one report.

    Sections carry their number in the heading; tables and figures carry
    theirs in the caption, which is where the numbering has to be read from:
    Report 2 has eight <figure> elements but four numbered figures, the rest
    being "Schematic 1" and the three "Money plot NR" panels.
    """
    txt = (ROOT / "reports" / (stem + ".template.html")).read_text(encoding="utf-8")
    secs = {m.group(1) for m in HEADING.finditer(txt)}
    tables, figures = set(), set()
    for m in CAPTION.finditer(txt):
        (tables if m.group(1) == "Table" else figures).add(int(m.group(2)))
    return secs, tables, figures


NUMBERS = _report_numbers()
INDEX = {n: _index(s) for n, s in NUMBERS.items()}


# --- the checks -------------------------------------------------------------

@check("citations: every [n] marker resolves to a reference entry")
def _():
    bad = []
    for path in REPORTS:
        lines = _lines(path)
        entries, span = _reference_entries(lines)
        if not entries:
            bad.append("%s has no <div class=\"refs\"> reference list" % path.name)
            continue
        stale = _stale_appendix_lines(lines)
        for i, line in enumerate(lines):
            if span[0] <= i <= span[1] or i in stale:
                continue
            for m in MARKER.finditer(line):
                for n in _expand(m.group(1)):
                    if n not in entries:
                        bad.append("%s:%d: %s cites [%d], which the reference "
                                   "list does not define (it ends at [%d])"
                                   % (path.name, i + 1, m.group(0), n,
                                      max(entries)))
    return bad


@check("citations: every reference entry is cited at least once")
def _():
    bad = []
    for path in REPORTS:
        lines = _lines(path)
        entries, span = _reference_entries(lines)
        if not entries:
            continue
        apx = _appendix_start(lines)
        cited = set()
        for i, line in enumerate(lines):
            if span[0] <= i <= span[1] or i >= apx:
                continue          # entry labels, and Appendix A history
            for m in MARKER.finditer(line):
                cited.update(_expand(m.group(1)))
        for n in sorted(set(entries) - cited):
            bad.append("%s:%d: reference [%d] is defined but never cited in "
                       "the body (an Appendix A row does not count as a "
                       "citation)" % (path.name, entries[n], n))
    return bad


@check("citations: every Report N pointer resolves in the report it names")
def _():
    bad = []
    for path in REPORTS:
        lines = _lines(path)
        stale = _stale_appendix_lines(lines)
        for i, line in enumerate(lines):
            if i in stale:
                continue
            for m in POINTER.finditer(line):
                num = int(m.group(1))
                if num not in INDEX:
                    bad.append("%s:%d: %s names a report build_report.py does "
                               "not know" % (path.name, i + 1, m.group(0)))
                    continue
                secs, tables, figures = INDEX[num]
                target = NUMBERS[num] + ".template.html"
                if m.group(2):
                    if m.group(2) not in secs:
                        bad.append("%s:%d: \"%s\" points at a section %s does "
                                   "not have (its headings run %s)"
                                   % (path.name, i + 1, m.group(0), target,
                                      ", ".join(sorted(
                                          secs, key=lambda s: [int(x) for x
                                                               in s.split(".")]))))
                    tail = m.end()
                else:
                    kind, k = m.group(3), int(m.group(4))
                    have = tables if kind == "Table" else figures
                    if k not in have:
                        bad.append("%s:%d: \"%s\" points at a %s %s does not "
                                   "have (it is numbered to %d)"
                                   % (path.name, i + 1, m.group(0),
                                      "table" if kind == "Table" else "figure",
                                      target, max(have) if have else 0))
                    tail = m.end()
                # "Report 3 SS 5 and Table 6" -- the tail names the same report
                while True:
                    c = CONTINUATION.match(line, tail)
                    if not c:
                        break
                    kind, k = c.group(1), int(c.group(2))
                    have = tables if kind == "Table" else figures
                    if k not in have:
                        bad.append("%s:%d: \"Report %d ... %s %d\" points at a "
                                   "%s %s does not have (it is numbered to %d)"
                                   % (path.name, i + 1, num, kind, k,
                                      "table" if kind == "Table" else "figure",
                                      target, max(have) if have else 0))
                    tail = c.end()
    return bad


@check("citations: every plans/NN #k pointer resolves to that item")
def _():
    bad = []
    for path in REPORTS:
        lines = _lines(path)
        stale = _stale_appendix_lines(lines)
        for i, line in enumerate(lines):
            if i in stale:
                continue
            for m in PLANS_ITEM.finditer(line):
                nn, k = m.group(1), int(m.group(2))
                hits = sorted(glob.glob(str(ROOT / "plans" / ("%s_*.md" % nn))))
                if not hits:
                    bad.append("%s:%d: %s names plans/%s, which does not exist"
                               % (path.name, i + 1, m.group(0), nn))
                    continue
                body = pathlib.Path(hits[0]).read_text(encoding="utf-8")
                if not re.search(r"^\s{0,3}%d\.\s" % k, body, re.M):
                    bad.append("%s:%d: %s points at item %d, which %s does not "
                               "number" % (path.name, i + 1, m.group(0), k,
                                           pathlib.Path(hits[0]).name))
    return bad
