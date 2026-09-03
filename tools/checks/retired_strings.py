"""Retired numbers and phrasings must not come back.

Motivated by ``docs/consistency_review_2026-09-02.md`` section "5.2 Mechanical
checks tools/consistency_check.py could add", item 2 ("A retired-strings
list").  The review's root-cause section 5.1 found the same pattern four
times: a development run re-derived a number, restated the section that owned
it, and left the copies in the other reports, the READMEs, the manual's second
table, the plans' live summaries and the front page where they were.  A
grep-able list of what each run retired is what catches that on the next
sweep rather than on the next review.

``tools/retired_numbers.json`` carries one entry per retired form::

    {"retired": "27\\u201335%",          # a regex
     "current": "22-31%",                # what replaced it
     "since":   "2026-08-29",            # when it was replaced
     "note":    "F133: the 6Li alpha tag at the lithium tagging optics",
     "allow":   ["in angle"]}            # optional: legitimate contexts

and this module greps the live documents for every ``retired`` pattern.

What counts as live.  The programme keeps its history in three places, and
none of them is a live claim:

  * an Appendix A revision row of a report template (``<tr><td class="mono">
    2026-08-28</td>…``) -- the row exists to say what the report used to say;
  * a ``## Development run N (DATE)`` section of ``plans/00_README.md`` dated
    before the entry's ``since`` -- the run log is written once and not
    rewritten.  A run section dated on or after ``since`` may quote a retired
    form only when the same line withdraws it ("withdrawn", "retired",
    "superseded") or quotes it as what the value is *not*;
  * a paragraph carrying a dated closing stamp earlier than ``since``
    ("**Superseded 2026-08-27.**", "*Done 2026-08-25.*") in a plan.

Two document classes are outside the corpus for the same reason: the dated
audits ``docs/code_review_*.md`` and ``docs/consistency_review_*.md`` are the
record of what was retired, so every retired form appears in them by
construction.
"""

import ast
import json
import pathlib
import re

ROOT = checker.ROOT
DATA = ROOT / "tools/retired_numbers.json"

# the dated audits are the record of the retirement itself
_HISTORY_DOC = re.compile(r"(?:code_review|consistency_review)_\d{4}-\d{2}-\d{2}\.md$")

_APPENDIX_A_ROW = re.compile(r'<tr>\s*<td class="mono">\d{4}-\d{2}-\d{2}</td>')
_RUN_HEADING = re.compile(r"^##\s+Development run\b[^(\n]*\((\d{4}-\d{2}-\d{2})")
_WITHDRAWN = re.compile(
    r"(?i)\b(?:withdrawn|withdraws|withdrew|retired|retires|superseded|supersedes)\b")
_NEGATED_QUOTE = re.compile(
    r"(?i)\b(?:not|rather than|in place of|instead of|no longer)\s*[\"“‘']?\s*$")
# a document may quote a retired form in order to refute it, as
# nearbeam_aperture_scan.py's docstring does: '"... the machine binds
# everywhere" was true of the September-2024 table and is not true of this
# one.'  That is allowed when the form is inside quotation marks and the
# refutation is in the same sentence.
_REFUTATION = re.compile(
    r"(?i)\b(?:was true of|is not true|not true of|no longer true|withdrawn|"
    r"retired|superseded|in place of|rather than)\b")
# a dated closing stamp on a plans entry: "**Superseded 2026-08-27.**",
# "*Done 2026-08-25.*", "| closed 2026-08-28"
_DATED_STAMP = re.compile(
    r"(?im)(?:^|[*_(|]|☑)\s*(?:superseded|withdrawn|retired|done|closed)\b"
    r"[^.\n]{0,25}?(\d{4}-\d{2}-\d{2})")

_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


def _norm(s):
    return re.sub(r"\s+", " ", s).strip()


def _module_docstring(path):
    """(source lines of the module docstring, 0-based offset), or (None, 0)."""
    try:
        src = path.read_text()
        tree = ast.parse(src, filename=str(path))
    except (OSError, SyntaxError, UnicodeDecodeError):
        return None, 0
    if not tree.body or not isinstance(tree.body[0], ast.Expr):
        return None, 0
    node = tree.body[0].value
    if not isinstance(node, ast.Constant) or not isinstance(node.value, str):
        return None, 0
    end = getattr(node, "end_lineno", node.lineno)
    return src.splitlines()[node.lineno - 1:end], node.lineno - 1


def _corpus():
    """[(path, lines, line-number offset)] over the documents the review names:
    the five templates, README.md, reports/index.html, docs/*.md, plans/*.md,
    the three package READMEs, and every module docstring under evgen/ and
    fastsim/."""
    out = []

    def take(p, offset=0, lines=None):
        if lines is None:
            try:
                lines = p.read_text().splitlines()
            except (OSError, UnicodeDecodeError):
                return
        out.append((p, lines, offset))

    for p in sorted(ROOT.glob("reports/*.template.html")):
        take(p)
    for rel in ("README.md", "reports/index.html",
                "evgen/README.md", "fastsim/README.md", "refs/README.md"):
        if (ROOT / rel).is_file():
            take(ROOT / rel)
    for p in sorted(ROOT.glob("docs/*.md")) + sorted(ROOT.glob("plans/*.md")):
        if _HISTORY_DOC.search(p.name):
            continue
        take(p)
    for base in ("evgen", "fastsim"):
        for p in sorted((ROOT / base).rglob("*.py")):
            doc, off = _module_docstring(p)
            if doc:
                take(p, off, doc)
    return out


def _run_sections(path):
    """[(first line, last line, date)] for plans/00's development-run log."""
    try:
        lines = path.read_text().splitlines()
    except OSError:
        return []
    heads = [(i, m.group(1)) for i, ln in enumerate(lines)
             for m in [_RUN_HEADING.match(ln)] if m]
    out = []
    for k, (i, date) in enumerate(heads):
        end = heads[k + 1][0] - 1 if k + 1 < len(heads) else len(lines) - 1
        out.append((i, end, date))
    return out


def _paragraph(lines, i):
    """The contiguous block of non-blank lines containing line i."""
    a = b = i
    while a > 0 and lines[a - 1].strip():
        a -= 1
    while b + 1 < len(lines) and lines[b + 1].strip():
        b += 1
    return "\n".join(lines[a:b + 1])


def _load():
    entries = json.loads(DATA.read_text())
    bad = []
    for k, e in enumerate(entries):
        for field in ("retired", "current", "since", "note"):
            if not isinstance(e.get(field), str) or not e[field]:
                bad.append("%s entry %d has no %s"
                           % (DATA.relative_to(ROOT), k, field))
        if isinstance(e.get("since"), str) and not _DATE.match(e["since"]):
            bad.append("%s entry %d: since = %r is not YYYY-MM-DD"
                       % (DATA.relative_to(ROOT), k, e["since"]))
        try:
            # case-insensitive: a retired phrase recurs just as readily at
            # the head of a sentence ("Tritons need IR-8") as inside one, and
            # a number carries no case, so nothing is lost by folding it
            e["_rx"] = re.compile(e["retired"], re.IGNORECASE)
            e["_allow"] = [re.compile(a) for a in e.get("allow", [])]
        except re.error as exc:
            bad.append("%s entry %d: %r does not compile: %s"
                       % (DATA.relative_to(ROOT), k, e.get("retired"), exc))
            e["_rx"] = None
    return entries, bad


def _windows(lines):
    """[(index, raw line, window text, length of the line's own part)] -- a
    two-line window with the whitespace collapsed, so a form the prose wraps
    ("the machine binds / everywhere") is still one string.  A match that
    begins past the line's own part belongs to the next window."""
    norm = [_norm(x) for x in lines]
    out = []
    for i, head in enumerate(norm):
        if not head:
            continue
        tail = norm[i + 1] if i + 1 < len(norm) else ""
        out.append((i, lines[i], (head + " " + tail) if tail else head,
                    len(head)))
    return out


@check("drift: no document carries a retired number or phrasing "
       "(tools/retired_numbers.json)")
def _():
    entries, bad = _load()
    if bad:
        return bad
    live = [e for e in entries if e["_rx"] is not None]
    try:                          # one alternation as a prefilter over the text
        any_rx = re.compile("|".join("(?:%s)" % e["retired"] for e in live),
                            re.IGNORECASE)
    except re.error:
        any_rx = None
    runs = _run_sections(ROOT / "plans/00_README.md")
    for path, lines, off in _corpus():
        is_template = path.name.endswith(".template.html")
        is_run_log = path == ROOT / "plans/00_README.md"
        for i, raw, window, headlen in _windows(lines):
            if any_rx is not None and not any_rx.search(window):
                continue
            for e in live:
                since = e["since"]
                m = e["_rx"].search(window)
                if not m or m.start() > headlen:
                    continue
                if any(a.search(window) for a in e["_allow"]):
                    continue
                # an Appendix A revision row is history, not a claim
                if is_template and _APPENDIX_A_ROW.search(raw):
                    continue
                # the form quoted, in quotation marks, in order to refute it
                before = window[:m.start()]
                if (before.count('"') + before.count("“") - before.count("”")
                        ) % 2 == 1 and _REFUTATION.search(window):
                    continue
                # plans/00's run log: sections written before the change stand
                sect = next((d for a, b, d in runs if a <= i <= b), None) \
                    if is_run_log else None
                if sect is not None:
                    if sect < since:
                        continue
                    if _WITHDRAWN.search(window[:headlen]) or \
                            _NEGATED_QUOTE.search(
                                window[max(0, m.start() - 24):m.start()]):
                        continue
                # a plans entry closed and dated before the change
                if [s for s in _DATED_STAMP.findall(_paragraph(lines, i))
                        if s < since]:
                    continue
                bad.append('%s:%d carries the retired "%s" -- since %s it is '
                           '"%s" (%s)'
                           % (path.relative_to(ROOT), off + i + 1,
                              m.group(0), since, e["current"], e["note"]))
    return bad
