"""The public front page and the repository map, against the reports.

Motivated by ``docs/consistency_review_2026-09-02.md`` section "5.2 Mechanical
checks tools/consistency_check.py could add", item 4 ("Front page and README
map").  ``reports/index.html`` and ``README.md`` are the two hand-maintained
documents of the programme, and section 5.1 records that they were the two
most often left behind: the front page still carried card dates that stopped
at 2026-08-27/28, a "27-35% alpha tag" the reports had restated to 22-31%,
"fifteen assumptions" against Report 3's sixteen, and the repository map
listed three unnumbered reports against five numbered ones (F119, F133, F138,
F144, F214, F238, F124, C004).

Four things are asserted here:

  * every date in a report's dateline appears in its card's date line, so a
    revision cannot be published without the front page saying so;
  * every number a card prints appears in that report's abstract -- the card
    is a precis of the abstract, so it may not carry a figure of its own;
  * every companion link in the footer resolves to a file in this repository;
  * ``README.md``'s repository-map row for ``reports/`` names every
    ``reports/*.template.html`` by its number, counts them correctly, and does
    not date the set earlier than the newest report dateline.
"""

import pathlib
import re

ROOT = checker.ROOT
INDEX = ROOT / "reports/index.html"
README = ROOT / "README.md"

_CARD = re.compile(r'<div class="card">(.*?)</div>\s*</div>', re.S)
_NUM_OF_CARD = re.compile(r'class="num">Report (\d+)<')
_CARD_DATE = re.compile(r'<div class="date">(.*?)</div>', re.S)
_CARD_BODY = re.compile(r"<p>(.*?)</p>", re.S)
_DATELINE = re.compile(r'<p class="dateline">(.*?)</p>', re.S)
_ABSTRACT = (re.compile(r'<div class="abstract">(.*?)</div>', re.S),
             re.compile(r'<p class="lede">(.*?)</p>', re.S))
_FOOTER = re.compile(r'<div class="foot">(.*?)</div>', re.S)
_BUILDER = re.compile(r'"stem":\s*"([a-z_0-9]+)",\s*"number":\s*(\d+)')
_DATE = re.compile(r"\d{4}-\d{2}-\d{2}")
# a bare quantity: not a "5R" money-plot label, not part of a longer number
_QUANTITY = re.compile(r"\d+(?:\.\d+)?(?![A-Za-z0-9.])")
_COUNT_WORD = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five",
               6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten"}


def _text(html):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html)).strip()


def _stems():
    """{report number: template path}, taken from the builder's own order."""
    src = (ROOT / "reports/build_report.py").read_text()
    return dict((int(n), ROOT / ("reports/%s.template.html" % s))
                for s, n in _BUILDER.findall(src))


def _cards():
    """[(number, date line text, body text, template path)] for each card."""
    out = []
    stems = _stems()
    for block in _CARD.findall(INDEX.read_text()):
        m = _NUM_OF_CARD.search(block)
        d = _CARD_DATE.search(block)
        b = _CARD_BODY.search(block)
        if not (m and d and b):
            out.append((None, block, None, None))
            continue
        out.append((int(m.group(1)), _text(d.group(1)), _text(b.group(1)),
                    stems.get(int(m.group(1)))))
    return out


def _dateline(tpl):
    m = _DATELINE.search(tpl.read_text())
    return _text(m.group(1)) if m else None


def _abstract(tpl):
    txt = tpl.read_text()
    for rx in _ABSTRACT:
        m = rx.search(txt)
        if m:
            return _text(m.group(1))
    return None


@check("front page: every index.html card names every date of its report's dateline")
def _():
    bad = []
    for num, date, _body, tpl in _cards():
        if num is None:
            bad.append("reports/index.html: a card carries no report number, "
                       "date line or summary paragraph")
            continue
        if tpl is None or not tpl.exists():
            bad.append("reports/index.html: card 'Report %d' has no template "
                       "in build_report.py" % num)
            continue
        dl = _dateline(tpl)
        if dl is None:
            bad.append("%s carries no <p class=\"dateline\">"
                       % tpl.relative_to(ROOT))
            continue
        for d in sorted(set(_DATE.findall(dl))):
            if d not in date:
                bad.append("reports/index.html: the Report %d card's date line "
                           "(%s) does not carry %s, which %s's dateline does"
                           % (num, date, d, tpl.name))
    return bad


@check("front page: every number on an index.html card appears in the report's abstract")
def _():
    bad = []
    for num, _date, body, tpl in _cards():
        if num is None or tpl is None or not tpl.exists() or body is None:
            continue                                # reported by the check above
        ab = _abstract(tpl)
        if ab is None:
            bad.append("%s carries no abstract (<div class=\"abstract\"> or "
                       "<p class=\"lede\">)" % tpl.relative_to(ROOT))
            continue
        have = set(_QUANTITY.findall(ab))
        for n in sorted(set(_QUANTITY.findall(body)) - have):
            bad.append("reports/index.html: the Report %d card prints %s, "
                       "which is nowhere in %s's abstract"
                       % (num, n, tpl.name))
    return bad


@check("front page: every companion link in the index.html footer resolves")
def _():
    txt = INDEX.read_text()
    m = _FOOTER.search(txt)
    if not m:
        return ["reports/index.html carries no <div class=\"foot\"> block"]
    line0 = txt[:m.start()].count("\n") + 1
    bad = []
    for href, label in re.findall(r'<a href="([^"]+)"[^>]*>(.*?)</a>',
                                  m.group(1), re.S):
        rel = re.sub(r"^https?://github\.com/[^/]+/[^/]+/blob/[^/]+/", "", href)
        if rel == href and "://" in href:
            continue                       # an external link, not ours to check
        if not (ROOT / rel).exists():
            bad.append("reports/index.html:%d+ the footer links %s, which does "
                       "not exist" % (line0, rel))
        elif _text(label) and "/" in _text(label) and _text(label) != rel:
            bad.append("reports/index.html:%d+ the footer link to %s is "
                       "labelled %s" % (line0, rel, _text(label)))
    return bad


@check("front page: README's repository map names every report template by number")
def _():
    lines = README.read_text().splitlines()
    rows = [(i + 1, ln) for i, ln in enumerate(lines)
            if re.match(r"\s*\|\s*`reports/`\s*\|", ln)]
    if not rows:
        return ["README.md has no repository-map row for `reports/`"]
    if len(rows) > 1:
        return ["README.md has %d repository-map rows for `reports/` (lines %s)"
                % (len(rows), ", ".join(str(n) for n, _ in rows))]
    ln, row = rows[0]
    stems = _stems()
    templates = sorted(ROOT.glob("reports/*.template.html"))
    bad = []
    for tpl in templates:
        nums = [n for n, p in stems.items() if p == tpl]
        if not nums:
            bad.append("README.md:%d — %s is not numbered in build_report.py, "
                       "so the map row cannot name it" % (ln, tpl.name))
            continue
        n = nums[0]
        if not re.search(r"\*\*%d\*\*|\bReport %d\b" % (n, n), row):
            bad.append("README.md:%d — the reports/ row does not name Report "
                       "%d (%s)" % (ln, n, tpl.name))
    word = re.search(r"\b(one|two|three|four|five|six|seven|eight|nine|ten)\b"
                     r"\s+circulate-able reports", row)
    if word and word.group(1) != _COUNT_WORD.get(len(templates)):
        bad.append("README.md:%d — the reports/ row says %s circulate-able "
                   "reports; reports/ holds %d templates (%s)"
                   % (ln, word.group(1), len(templates),
                      _COUNT_WORD.get(len(templates), len(templates))))
    # the row dates the whole set; that date may not predate a report's own
    newest = ""
    for tpl in templates:
        dl = _dateline(tpl) or ""
        newest = max([newest] + _DATE.findall(dl))
    claim = re.search(r"\ball\b[^.;|]{0,40}?(\d{4}-\d{2}-\d{2})", row)
    if claim and newest:
        # a row may date the set more than once -- "all restated <date> ...
        # and corrected <later date> ..." -- and it is the LAST of those
        # dates that has to reach the newest dateline; requiring it of the
        # first would fail a row that names the newest revision too.
        carried = _DATE.findall(row[claim.start():])
        if carried and max(carried) < newest:
            bad.append("README.md:%d — the reports/ row dates the set '%s', "
                       "but the newest report dateline is %s"
                       % (ln, _text(claim.group(0)), newest))
    return bad
