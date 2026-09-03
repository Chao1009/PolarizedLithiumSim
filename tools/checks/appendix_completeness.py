"""Appendix A revision rows against what git says actually changed.

Motivated by docs/consistency_review_2026-09-02.md section 5.2 item 9
("Appendix A completeness against git"): the newest revision row of each
report is supposed to say which tables the revision moved, and three of the
review's findings (F109, F111, F116) were tables whose cells had changed
under a row that never named them.

For every reports/*.template.html this reads the diff of the revision the
newest Appendix A row describes -- `git log -n 3 --format=%H -- <template>`
for the previous commit that touched the file, `git diff <hash> -- <template>`
against the working tree, or `git diff HEAD -- <template>` when the working
tree is dirty -- keeps the changed lines that lie inside a <table> element and
carry a <td> cell, maps each to the table it belongs to, and asserts the row
names every one of them as "Table k".

Scoping, so that the programme's conventions are not flagged:

  * Only the newest row is read, so the older rows -- history, which is not
    restated when a later run moves a table -- are untouched by construction.
  * A table is identified by its own caption, <p class="tabcap"><b>Table k,
    which follows the </table> it belongs to.  The Appendix A revision table
    itself has no such caption, so its cells -- which change at every
    revision -- are not counted, and neither is any other uncaptioned table.
  * Only <td> lines count.  A revision that restates a caption, a heading or
    a paragraph is not a change to the table's cells, which is the
    distinction the review draws ("the row need only name tables whose cells
    changed, not captions"): Report 3's Table 2 caption and its beta* source
    cell moved together on 2026-09-02, and it is the cell that obliges the
    row to name the table.
  * "Table 10.1", "YR Table 10.3" and "Report 2 Table 2" are not this
    report's tables: a number may not be followed by a digit or a decimal
    point, and a "YR" / "Yellow Report" / "Report N" prefix disqualifies it.
"""

import glob
import pathlib
import re
import subprocess

ROOT = checker.ROOT
REPORTS = sorted(pathlib.Path(p) for p in
                 glob.glob(str(ROOT / "reports/*.template.html")))

TABLE_OPEN = re.compile(r"<table\b")
TABLE_CLOSE = re.compile(r"</table>")
TABLE_CAPTION = re.compile(r'class="tabcap"[^>]*>\s*<b>\s*Table\s+(\d+)\b')
APPENDIX_A = re.compile(r"<h2>\s*Appendix\s+A\b", re.I)
REV_ROW = re.compile(r'<td\s+class="mono">\s*(\d{4}-\d\d-\d\d)\s*</td>\s*<td>(.*)$')
# this report's own tables only, and never the "10" of "Table 10.1"
OWN_TABLE = re.compile(r"(?:(YR|Yellow Report|Report\s+\d)\s+)?Table\s+(\d+)(?!\.?\d)")
HUNK = re.compile(r"^@@ -(\d+)(?:,\d+)? \+(\d+)(?:,\d+)? @@")


def _git(*args):
    """git stdout, or None when the command (or git itself) fails."""
    try:
        p = subprocess.run(("git",) + args, cwd=str(ROOT),
                           capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return p.stdout if p.returncode == 0 else None


def _table_spans(lines):
    """[(first, last, number or None)] for every <table> element, in order."""
    out, start = [], None
    for i, line in enumerate(lines):
        if TABLE_OPEN.search(line):
            start = i
        if TABLE_CLOSE.search(line) and start is not None:
            number = None
            for j in range(i + 1, min(i + 4, len(lines))):
                m = TABLE_CAPTION.search(lines[j])
                if m:
                    number = int(m.group(1))
                    break
                if lines[j].strip():
                    break
            out.append((start, i, number))
            start = None
    return out


def _table_of(spans, idx):
    for lo, hi, number in spans:
        if lo <= idx <= hi:
            return number
    return None


def _newest_row(lines):
    """(date, line number, text) of the latest-dated Appendix A row."""
    for i, line in enumerate(lines):
        if APPENDIX_A.search(line):
            break
    else:
        return None
    rows = []
    for j in range(i, len(lines)):
        m = REV_ROW.search(lines[j])
        if m:
            rows.append((m.group(1), j + 1, m.group(2)))
    if not rows:
        return None
    return max(rows, key=lambda r: r[0])


@check("appendix: the newest revision row names every table whose cells changed")
def _():
    if _git("rev-parse", "--git-dir") is None:
        return []
    status = _git("status", "--porcelain", "--", "reports") or ""
    dirty = {line[3:].strip() for line in status.splitlines()}
    bad = []
    for path in REPORTS:
        rel = "reports/" + path.name
        row = _newest_row(path.read_text(encoding="utf-8").splitlines())
        if row is None:
            bad.append("%s has no dated Appendix A revision row" % path.name)
            continue
        date, row_line, text = row
        log = _git("log", "-n", "3", "--format=%H", "--", rel)
        if not log:
            continue                       # not committed yet: nothing to diff
        hashes = log.split()
        if rel in dirty:
            base = "HEAD"
        elif len(hashes) > 1:
            base = hashes[1]
        else:
            continue                       # only ever one commit: no previous
        old = _git("show", "%s:%s" % (base, rel))
        diff = _git("diff", "-U0", base, "--", rel)
        if old is None or diff is None:
            continue
        new_spans = _table_spans(path.read_text(encoding="utf-8").splitlines())
        old_spans = _table_spans(old.splitlines())
        changed = {}                       # table number -> example line
        old_line = new_line = 0
        for line in diff.splitlines():
            m = HUNK.match(line)
            if m:
                old_line, new_line = int(m.group(1)), int(m.group(2))
                continue
            if line.startswith(("+++", "---")):
                continue
            if line.startswith("+"):
                if "<td" in line:
                    k = _table_of(new_spans, new_line - 1)
                    if k is not None and changed.get(k) is None:
                        changed[k] = new_line
                new_line += 1
            elif line.startswith("-"):
                if "<td" in line:
                    k = _table_of(old_spans, old_line - 1)
                    if k is not None:
                        changed.setdefault(k, None)
                old_line += 1
        named = {int(m.group(2)) for m in OWN_TABLE.finditer(text)
                 if not m.group(1)}
        for k in sorted(set(changed) - named):
            where = ("line %d" % changed[k]) if changed[k] else "a deleted row"
            since = ("in the working tree" if base == "HEAD"
                     else "since %s" % base[:8])
            bad.append("%s:%d: the %s revision row does not name Table %d, "
                       "whose cells changed %s (%s); the row names %s"
                       % (path.name, row_line, date, k, since, where,
                          ", ".join("Table %d" % n for n in sorted(named))
                          or "no table"))
    return bad
