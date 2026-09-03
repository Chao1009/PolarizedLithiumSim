"""Figures against captions -- docs/consistency_review_2026-09-02.md section 5.2,
numbered item 6 ("Figures against captions").

Every embedded PNG is drawn by one evgen/scripts/*.py and described by one
<figcaption> per template that embeds it.  The two drift apart silently: a
caption keeps a panel letter the script stopped drawing (F186, F233), quotes a
label the panel never carried (F155, F180), or reads numbers off an axis whose
floor hides them (F131, F249).  Nothing in the build notices, because the
caption is prose and the panel is a raster.

The five checks below tie the two together through build_report.py's
__TAG__ -> PNG registry: the template's <figure> block gives the caption, the
registry gives the PNG, and checker.script_of_figure() gives the script.

  (b) every "(a)".."(d)" a caption uses is drawn by the script, and no script
      whose every figure holds a single axis draws panel letters at all;
  (c) every label a caption puts in quotes is a substring of a string the
      script carries;
  (d) a caption -- or prose citing one of its panels by letter -- that reads a
      value off a panel below that panel's set_ylim floor must record the
      floor, or say the markers are off the axis.  Only floors at or below
      1e-6 are treated as clipping (DEEP_FLOOR): a floor of 0.003 on a
      delta-y/y axis frames the data, it does not hide markers.

Review item (a), the title-extent check, is implemented under --full, but not
in the form the review sketches: re-rendering "the script's suptitle strings"
is impossible without running the script, because every long title in this
repository is a %-format whose runtime substitutions (the configuration label,
the fitted amplitudes) are most of its width -- rendering the literal would
measure a string nobody publishes and would have missed F183, F184 and F211
outright.  The published PNG is measured instead: matplotlib does not wrap or
shrink an overlong title, it draws it past the canvas and the raster cuts it,
so ink in the outermost two pixel rows or columns of the PNG is exactly the
clipping those three findings report.
"""

import ast
import html
import pathlib
import re

ROOT = checker.ROOT                                          # noqa: F821
REPORTS = sorted(ROOT.glob("reports/*.template.html"))

# --- the template side ------------------------------------------------------

FIGURE_RE = re.compile(r"<figure>(.*?)</figure>", re.S)
FIGCAP_RE = re.compile(r"<figcaption>(.*?)</figcaption>", re.S)
TAG_RE = re.compile(r"__[A-Z0-9_]+__")
SUP = {"\u2070": "0", "\u00b9": "1", "\u00b2": "2", "\u00b3": "3", "\u2074": "4",
       "\u2075": "5", "\u2076": "6", "\u2077": "7", "\u2078": "8", "\u2079": "9"}


def _plain(fragment):
    """A caption's HTML as the reader sees it: <sup>-9</sup> folded to the
    superscript digits the rest of the prose uses, tags dropped, entities
    resolved, whitespace collapsed."""
    def sup(m):
        return "".join({v: k for k, v in SUP.items()}.get(c, "\u207b" if c in "-\u2212" else c)
                       for c in m.group(1))
    fragment = re.sub(r"<sup>\s*([\-\u2212\d]+)\s*</sup>", sup, fragment)
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", fragment))).strip()


def _figures():
    """(template, line of the <figcaption>, tag, png path, caption text) for
    every PNG a template embeds.  A <figure> holding an inline SVG carries no
    tag and is skipped; a tag build_report.py does not register is the
    business of the artefacts check, not this one."""
    registry = dict(checker.registered_figures())             # noqa: F821
    out = []
    for tpl in REPORTS:
        text = tpl.read_text()
        for blk in FIGURE_RE.finditer(text):
            tags = TAG_RE.findall(blk.group(1))
            cap = FIGCAP_RE.search(blk.group(1))
            if not tags or not cap or tags[0] not in registry:
                continue
            line = text[:blk.start() + cap.start()].count("\n") + 1
            out.append((tpl, line, tags[0], registry[tags[0]], _plain(cap.group(1))))
    return out


_BODY = {}
APPENDIX_A = re.compile(r"<h2>\s*Appendix\s+A\b", re.I)


def _body_before_appendix(tpl):
    """The template as prose, cut at the Appendix A heading.  A revision row
    recounts what a figure used to show and quotes the numbers it used to
    carry; that history is the programme's convention and must not be read as
    a live claim about the panel."""
    if tpl not in _BODY:
        text = tpl.read_text()
        cut = APPENDIX_A.search(text)
        _BODY[tpl] = _plain(text[:cut.start()] if cut else text)
    return _BODY[tpl]


def _caption_number(caption):
    """The "Figure 3" a caption opens with, so prose citing panel 3c can be
    found; None when the caption does not number itself."""
    m = re.match(r"\s*Figure\s+(\d+)\b", caption)
    return m.group(1) if m else None


# --- the script side --------------------------------------------------------

TITLE_CALLS = {"set_title", "suptitle", "title", "annotate", "text"}
_CACHE = {}


def _tree(script):
    if script not in _CACHE:
        try:
            _CACHE[script] = ast.parse(pathlib.Path(script).read_text(), filename=str(script))
        except (OSError, SyntaxError):
            _CACHE[script] = None
    return _CACHE[script]


def _strings_in(node):
    return [n.value for n in ast.walk(node)
            if isinstance(n, ast.Constant) and isinstance(n.value, str)]


def _norm(s):
    r"""Mathtext and TeX spacing removed, so a caption's plain "legacy, 275 GeV
    p" can be found in the script's "(legacy, 275 GeV $p$)"."""
    s = s.replace("$", "")
    s = re.sub(r"\\[,;!:> ]", "", s)
    s = re.sub(r"\\(?:rm|mathrm|bf|it|langle|rangle)\b", "", s)
    return re.sub(r"\s+", " ", s).strip().lower()


def _drawn_strings(script):
    """Strings the script draws on the canvas: the arguments of set_title /
    suptitle / title / annotate / text, plus -- for the common idiom that
    carries the titles in a loop tuple,

        for ax, ylab, title in ((ax1, "...", "(a) ..."), (ax2, ...)):
            ax.set_title(title, fontsize=9.5)

    -- every string bound anywhere in the module to a name such a call passes.
    Without the second half nearbeam_aperture_scan.py's "(a)" and "(b)" are
    invisible, because its set_title sees only the loop variable."""
    tree = _tree(script)
    if tree is None:
        return []
    bound = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            targets = [t for t in node.targets]
        elif isinstance(node, ast.For):
            targets = [node.target]
        else:
            continue
        names = [n.id for t in targets for n in ast.walk(t) if isinstance(n, ast.Name)]
        value = node.value if isinstance(node, ast.Assign) else node.iter
        strings = _strings_in(value)
        for name in names:
            bound.setdefault(name, []).extend(strings)
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        name = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", "")
        if name not in TITLE_CALLS:
            continue
        out.extend(_strings_in(node))
        for arg in ast.walk(node):
            if isinstance(arg, ast.Name):
                out.extend(bound.get(arg.id, ()))
    return out


def _all_strings(script):
    tree = _tree(script)
    return _strings_in(tree) if tree is not None else []


def _axes_per_figure(script):
    """How many axes each figure the script creates holds, as a static list.
    A plt.subplots() with no grid is one axis; a literal grid is its product;
    anything the parser cannot evaluate, and any add_subplot / GridSpec /
    subplot_mosaic, is counted as a multi-panel figure (2), because guessing
    "single axis" is the answer that would raise a false alarm."""
    tree = _tree(script)
    if tree is None:
        return []
    counts = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        name = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", "")
        if name in ("add_subplot", "subplot", "subplot_mosaic", "GridSpec", "add_gridspec"):
            counts.append(2)
        elif name == "subplots":
            grid = list(node.args[:2]) + [k.value for k in node.keywords
                                          if k.arg in ("nrows", "ncols")]
            if not grid:
                counts.append(1)                     # plt.subplots(figsize=...)
            elif all(isinstance(a, ast.Constant) and isinstance(a.value, int)
                     and not isinstance(a.value, bool) for a in grid):
                n = 1
                for a in grid:
                    n *= a.value
                counts.append(max(n, 1))
            else:
                counts.append(2)                     # plt.subplots(1, len(slices))
    return counts


def _stem_candidates(rel):
    """The spellings of a figure's name a script may carry, in the order
    checker.script_of_figure() tries them."""
    stem = pathlib.Path(rel).stem
    bare = stem.replace("_6Li", "").replace("_7Li", "")
    cands = [stem, bare]
    trimmed = re.sub(r"_[a-z0-9]+$", "", bare)
    if trimmed != bare:
        cands.append(trimmed)
    return cands


def _ylim_floors(script, rel):
    """The literal, positive lower bounds of the set_ylim calls that belong to
    the figure at `rel`.  A script that writes more than one PNG writes each in
    its own function (reco_chain_figures.inclusive_figure / coherent_figure),
    so the function naming this PNG scopes the search; otherwise the whole
    module does."""
    tree = _tree(script)
    if tree is None:
        return []
    funcs = [n for n in ast.walk(tree)
             if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
    scope = tree
    for cand in _stem_candidates(rel):      # most specific spelling first, so
        named = [f for f in funcs           # "reco_chain_coherent_6Li" is not
                 if any(cand in s for s in _strings_in(f))]   # served by the
        if named:                           # "reco_chain" of its sibling
            if len(named) == 1:
                scope = named[0]
            break                           # ambiguous: the whole module
    floors = []
    for node in ast.walk(scope):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        name = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", "")
        if name not in ("set_ylim", "ylim"):
            continue
        lows = []
        if node.args and isinstance(node.args[0], ast.Constant):
            lows.append(node.args[0].value)
        elif node.args and isinstance(node.args[0], (ast.Tuple, ast.IfExp)):
            for t in ast.walk(node.args[0]):
                if isinstance(t, ast.Tuple) and t.elts and isinstance(t.elts[0], ast.Constant):
                    lows.append(t.elts[0].value)
        for k in node.keywords:
            if k.arg in ("bottom", "ymin") and isinstance(k.value, ast.Constant):
                lows.append(k.value.value)
        floors.extend(float(v) for v in lows
                      if isinstance(v, (int, float)) and not isinstance(v, bool) and v > 0)
    return floors


# --- numbers written in prose ----------------------------------------------

SCI_RE = re.compile(
    r"(?:(\d+(?:\.\d+)?)\s*[\u00d7x\u2715]\s*)?"
    r"10\s*([\u207b\u2212\-]?)\s*([\u2070\u00b9\u00b2\u00b3\u2074-\u2079]+)")


def _sci_values(text):
    """Every "2.0 x 10^-19" the prose writes with superscript digits."""
    out = []
    for mant, sign, digits in SCI_RE.findall(text):
        exp = int("".join(SUP[c] for c in digits))
        if sign:
            exp = -exp
        out.append(float(mant or 1.0) * 10.0 ** exp)
    return out


def _mentions(value, text):
    """Does the prose write `value` -- as a power of ten, or as a decimal?"""
    for v in _sci_values(text):
        if abs(v - value) <= 1e-9 * max(abs(v), abs(value)):
            return True
    for tok in re.findall(r"\d+(?:\.\d+)?(?:e[+\-]?\d+)?", text.replace(",", "")):
        try:
            if abs(float(tok) - value) <= 1e-9 * max(abs(float(tok)), abs(value)):
                return True
        except ValueError:
            pass
    return False


# --- (b) panel letters ------------------------------------------------------

# "(d)" straight after a figure number -- "Figure 5(d)", "money plot 6R(d)" --
# points at another figure's panel and says nothing about this one.
PANEL_RE = re.compile(r"(?<![A-Za-z0-9])\(([a-d])\)")
# a caption may carry the letters itself when the panels are titled rather than
# lettered, provided it says where each one sits (Report 1 Figure 5, F186)
DECLARED_RE = re.compile(
    r"\(([a-d])\)\s*(?:is\s+|the\s+)?(?:top|bottom|upper|lower|left|right)\b")


@check("figures: every panel letter a caption uses is drawn in the figure")
def _():
    bad = []
    for tpl, line, tag, rel, cap in _figures():
        script = checker.script_of_figure(rel)                        # noqa: F821
        if script is None:
            continue
        used = set(PANEL_RE.findall(cap))
        if not used:
            continue
        declared = set(DECLARED_RE.findall(cap))
        drawn = _norm(" | ".join(_drawn_strings(script)))
        missing = sorted(L for L in used - declared if "(%s)" % L not in drawn)
        if missing:
            bad.append("%s:%d (%s -> %s) uses panel %s, which %s never draws in a "
                       "title, annotation or text string, and the caption does not "
                       "say where the panel sits"
                       % (tpl.name, line, tag, rel,
                          ", ".join("(%s)" % L for L in missing),
                          pathlib.Path(script).name))
    return bad


@check("figures: no panel letters are drawn in a script that saves one PNG per axis")
def _():
    bad = []
    seen = set()
    for tpl, line, tag, rel, _cap in _figures():
        script = checker.script_of_figure(rel)                        # noqa: F821
        if script is None or script in seen:
            continue
        seen.add(script)
        counts = _axes_per_figure(script)
        if not counts or max(counts) > 1:
            continue
        drawn = _drawn_strings(script)
        letters = sorted({m for s in drawn for m in PANEL_RE.findall(_norm(s))})
        if letters:
            bad.append("%s draws panel %s although every figure it saves holds one "
                       "axis (%s is one of them, %s:%d) -- the letters read as panels "
                       "of one image and the reports print the PNGs as separate figures"
                       % (pathlib.Path(script).name,
                          ", ".join("(%s)" % L for L in letters), rel, tpl.name, line))
    return bad


# --- (c) quoted labels ------------------------------------------------------

QUOTE_RE = re.compile(r"\"([^\"\n]{2,90})\"|\u201c([^\u201d\n]{2,90})\u201d")
LABEL_WORDS = re.compile(r"label|legend|title|annotat|axis|panel|marked|reads?\b",
                         re.I)


@check("figures: every label a caption quotes appears in the producing script")
def _():
    bad = []
    for tpl, line, tag, rel, cap in _figures():
        script = checker.script_of_figure(rel)                        # noqa: F821
        if script is None:
            continue
        pool = [_norm(s) for s in _all_strings(script)]
        for m in QUOTE_RE.finditer(cap):
            quoted = m.group(1) or m.group(2)
            near = cap[max(0, m.start() - 60):m.end() + 60]
            if not LABEL_WORDS.search(near):
                continue          # a quotation, not a claim about the panel
            if not any(_norm(quoted) in s for s in pool):
                bad.append("%s:%d (%s -> %s) attributes the label \u201c%s\u201d to the "
                           "panel, but no string in %s contains it"
                           % (tpl.name, line, tag, rel, quoted,
                              pathlib.Path(script).name))
    return bad


# --- (d) axis floors --------------------------------------------------------

# Only a floor this deep clips markers rather than merely framing the data: a
# log axis stopped at 1e-6 or below is decades short of the values the prose
# reads off it (F131's 1e-9, F249's 1e-18), while reco_chain_figures' 0.003 on
# a delta-y/y axis and money_cos2phi_coherent's 0.001 on a dsigma/dt one are
# plot ranges, and a caption's unrelated 4e-5 acceptance is not a marker.
DEEP_FLOOR = 1e-6

# How far back of a "Figure 3c" citation its numbers may sit.  A sentence
# split is not enough: the clause F131 turns on carries a semicolon inside
# its own parenthesis, which would cut the three values off the citation.
CITE_WINDOW = 400

OFF_AXIS = re.compile(
    r"off the (?:axis|scale)|below the axis|below panel|under the axis|"
    r"not drawn|rather than drawing|states its own maximum|"
    r"outside the (?:axis|panel)", re.I)


@check("figures: a caption quoting a value below the axis floor records the floor")
def _():
    bad = []
    for tpl, line, tag, rel, cap in _figures():
        script = checker.script_of_figure(rel)                        # noqa: F821
        if script is None:
            continue
        floors = [f for f in _ylim_floors(script, rel) if f <= DEEP_FLOOR]
        if not floors:
            continue
        floor = min(floors)
        # the caption's own numbers, plus those of any sentence that cites one
        # of this figure's panels ("Figure 3c", "Figure 5(d)") -- the values
        # F131 reads off panel (c) sit in Report 2 section 4.5, not in the caption
        prose = [cap]
        num = _caption_number(cap)
        if num:
            body = _body_before_appendix(tpl)
            cite = re.compile(r"Figures?\s+%s\s*\(?[a-d]\)(?![\w])|"
                              r"Figures?\s+%s\s*[a-d](?![\w])" % (num, num))
            for m in cite.finditer(body):
                prose.append(body[max(0, m.start() - CITE_WINDOW):m.end()])
        below = sorted({v for p in prose for v in _sci_values(p) if 0 < v < floor})
        if not below:
            continue
        if _mentions(floor, cap) or OFF_AXIS.search(cap):
            continue
        bad.append("%s:%d (%s -> %s) reads %s off a panel whose %s floor is %g "
                   "(%s set_ylim), and the caption neither gives the floor nor says "
                   "the markers are off the axis"
                   % (tpl.name, line, tag, rel,
                      ", ".join("%.3g" % v for v in below[:3]),
                      "y-axis", floor, pathlib.Path(script).name))
    return bad


# --- (a) title extent, --full only -----------------------------------------

@check("figures: no title overruns the canvas of the figure it is drawn on")
def _():
    if not checker.FULL:                                              # noqa: F821
        return checker.SKIP
    try:
        import numpy as np
        from PIL import Image
    except ImportError:                       # no raster reader: nothing to say
        return []
    bad = []
    for tpl, line, tag, rel, _cap in _figures():
        png = ROOT / rel
        if not png.exists():
            continue
        ink = np.asarray(Image.open(png).convert("L")) < 200
        edges = [("top", ink[:2, :]), ("bottom", ink[-2:, :]),
                 ("left", ink[:, :2]), ("right", ink[:, -2:])]
        hit = [name for name, band in edges if band.any()]
        if hit:
            script = checker.script_of_figure(rel)                    # noqa: F821
            bad.append("%s touches the %s edge of its canvas (%s:%d, %s) -- a title "
                       "or label is cut off; widen the figure or shorten the string "
                       "in %s"
                       % (rel, " and ".join(hit), tpl.name, line, tag,
                          pathlib.Path(script).name if script else "the producing script"))
    return bad
