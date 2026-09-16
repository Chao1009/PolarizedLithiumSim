#!/usr/bin/env python3
"""The letter's figure style: one rcParams context, two column widths, one
palette.

plans/07 WP6 asks the figure pass for "a shared matplotlib style (column
widths, 8-9 pt fonts, consistent Okabe-Ito accents)", and the run-19 ledger
(T23) states the acceptance as: the four figures build at the journal column
width with every font >= 8 pt, asserted in the driver.  This module is that
style and that assertion; it draws nothing and knows no physics.

  * `rc()` is the one context manager every driver opens.  Nothing in it is
    set anywhere else: a driver that wants a different size says so through a
    keyword of `figure()`, not by touching rcParams.
  * `COLUMN_W` / `DOUBLE_W` are the elsarticle two-column text widths of
    Physics Letters B, in inches: 3.375 in for a `figure`, 7.0 in for a
    `figure*`.  A figure is drawn AT its final size and never scaled in
    LaTeX, which is the only way a declared 8 pt stays 8 pt on the page.
  * `PALETTE` is Okabe-Ito (Okabe & Ito, "Color Universal Design", 2002), the
    palette the producing scripts already use: `money_cos2phi.py` and its two
    companions carry C_TRUTH = #0072B2, C_FIT = #D55E00, C_ALT = #009E73, and
    `TRUTH` / `FIT` / `ALT` below are those three under the names the letter
    uses, so a letter figure and the published figure it condenses are the
    same colour for the same thing.
  * `save()` writes the PDF LaTeX includes and a PNG preview beside it, and
    asserts the font floor on the way out; `assert_min_fontsize()` is the
    assertion on its own for a driver that wants it earlier.

Every font size below is 8 or 9 pt.  Tick labels, legends and annotations sit
at the 8 pt floor; axis labels and panel titles at 9 pt.  There is no 6.5 pt
annotation anywhere in `paper/`, which is the one real difference from the
published figures: those are 11-13 in wide and are read on a screen.
"""

import contextlib
import logging
import pathlib

import matplotlib

# matplotlib subsets the embedded fonts through fontTools, which logs one INFO
# line per glyph table; a driver that writes four PDFs would bury its own
# output under a few hundred of them.
logging.getLogger("fontTools").setLevel(logging.WARNING)
logging.getLogger("fontTools.subset").setLevel(logging.WARNING)


# --- widths -----------------------------------------------------------------

COLUMN_W = 3.375   # in -- one elsarticle column (PLB two-column preprint)
DOUBLE_W = 7.0     # in -- the full text width, for a figure* float

# --- the Okabe-Ito palette --------------------------------------------------

PALETTE = {
    "black":     "#000000",
    "orange":    "#E69F00",
    "sky":       "#56B4E9",
    "green":     "#009E73",
    "yellow":    "#F0E442",
    "blue":      "#0072B2",
    "vermilion": "#D55E00",
    "purple":    "#CC79A7",
}

# The three roles the producing scripts fix, under the letter's names.
TRUTH = PALETTE["blue"]        # injected / model curve        (#0072B2)
FIT = PALETTE["vermilion"]     # fitted value, 1-year points   (#D55E00)
ALT = PALETTE["green"]         # alternative interpretation    (#009E73)
DATA = PALETTE["black"]        # 10-year points, pseudo-data
GREY = "0.45"                  # anchors, references, guides

# The cycle a driver gets if it plots without saying a colour.
CYCLE = [TRUTH, FIT, ALT, PALETTE["orange"], PALETTE["purple"],
         PALETTE["sky"], PALETTE["yellow"]]

# --- font sizes -------------------------------------------------------------

FONT_MIN = 8.0     # the floor the drivers assert
FONT_SMALL = 8.0   # ticks, legends, annotations
FONT_LARGE = 9.0   # axis labels, panel titles

RC = {
    # Fonts.  STIXGeneral ships with matplotlib, so the figures look the same
    # on any machine that can run the drivers, and it is the serif face that
    # sits closest to the journal's own.
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "text.usetex": False,
    "font.size": FONT_SMALL,
    "axes.labelsize": FONT_LARGE,
    "axes.titlesize": FONT_LARGE,
    "figure.titlesize": FONT_LARGE,
    "xtick.labelsize": FONT_SMALL,
    "ytick.labelsize": FONT_SMALL,
    "legend.fontsize": FONT_SMALL,
    "legend.title_fontsize": FONT_SMALL,

    # Lines and marks, scaled for a 3.4 in column rather than a 13 in screen.
    "axes.prop_cycle": matplotlib.cycler(color=CYCLE),
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.3,
    "lines.markersize": 3.2,
    "patch.linewidth": 0.7,
    "errorbar.capsize": 1.6,
    "grid.linewidth": 0.5,

    # Ticks: inward on all four sides, the journal convention.
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.minor.size": 1.6,
    "ytick.minor.size": 1.6,

    # Legends: no shadow, thin frame, tight.
    "legend.frameon": True,
    "legend.framealpha": 0.92,
    "legend.edgecolor": "0.8",
    "legend.fancybox": False,
    "legend.borderpad": 0.35,
    "legend.labelspacing": 0.3,
    "legend.handlelength": 1.5,
    "legend.handletextpad": 0.5,
    "legend.columnspacing": 1.0,

    # Output: embed TrueType (type 42) rather than type 3, which is what
    # publishers ask for, and keep the PDF vector.
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    # NOT "tight".  A tight bounding box crops the canvas to the ink, so the
    # saved PDF is no longer COLUMN_W wide, LaTeX scales it back up to the
    # column, and every 8 pt glyph lands on the page at some other size.  The
    # drivers place their axes with subplots_adjust instead, and what is saved
    # is exactly the figure that was declared.
    "savefig.bbox": "standard",
    "savefig.pad_inches": 0.0,
    "figure.dpi": 200,
    "savefig.dpi": 400,
    "axes.unicode_minus": False,
}


@contextlib.contextmanager
def rc(**overrides):
    """The letter's rcParams.  Every driver draws inside this context."""
    params = dict(RC)
    params.update(overrides)
    with matplotlib.rc_context(params):
        yield


def figure(width=COLUMN_W, height=None, **kwargs):
    """A figure at a journal width.  `height` defaults to the golden ratio."""
    import matplotlib.pyplot as plt
    if height is None:
        height = width / 1.618
    return plt.figure(figsize=(width, height), **kwargs)


def subplots(nrows=1, ncols=1, width=COLUMN_W, height=None, **kwargs):
    """`plt.subplots` at a journal width."""
    import matplotlib.pyplot as plt
    if height is None:
        height = width / 1.618
    return plt.subplots(nrows, ncols, figsize=(width, height), **kwargs)


def panel_label(ax, text, x=0.025, y=0.965, **kwargs):
    """The (a) / (b) tag in a panel's corner, at the 8 pt floor."""
    kwargs.setdefault("fontsize", FONT_SMALL)
    kwargs.setdefault("ha", "left")
    kwargs.setdefault("va", "top")
    kwargs.setdefault("fontweight", "bold")
    return ax.annotate(text, xy=(x, y), xycoords="axes fraction", **kwargs)


# --- the font-floor assertion ----------------------------------------------

def font_sizes(fig):
    """[(size, text)] for every non-empty visible Text in the figure."""
    import matplotlib.text as mtext
    out = []
    for obj in fig.findobj(mtext.Text):
        if not obj.get_visible():
            continue
        s = obj.get_text()
        if not s or not s.strip():
            continue
        out.append((float(obj.get_fontsize()), s))
    return out


def undersized(fig, minimum=FONT_MIN, tol=1e-6):
    """Every text in `fig` below the floor, smallest first."""
    bad = [(sz, s) for sz, s in font_sizes(fig) if sz < minimum - tol]
    return sorted(bad)


def raise_small_text(artist, minimum=FONT_MIN):
    """Bring every text in `artist` (a figure or an axes) up to the floor.

    The drivers re-use drawing helpers of the producing scripts -- the y and
    W^2 guide lines of `phase_space_bins.draw_guides`, for one -- and those
    scripts label at 6.5 pt because they draw at 13 in wide.  Their labels are
    raised here rather than redrawn, so that the letter figure keeps the
    helper as its single source and still satisfies the floor.  Returns the
    number of texts raised.
    """
    import matplotlib.text as mtext
    n = 0
    for obj in artist.findobj(mtext.Text):
        if obj.get_visible() and obj.get_text().strip() \
                and float(obj.get_fontsize()) < minimum:
            obj.set_fontsize(minimum)
            n += 1
    return n


def assert_min_fontsize(fig, minimum=FONT_MIN):
    """T23's acceptance, as an assertion: no glyph below `minimum` points.

    Called by `save()`, so a driver gets it for free; call it directly to
    fail before a slow save.  The figure must already be drawn far enough
    for its text to exist -- after tight_layout, before savefig.
    """
    bad = undersized(fig, minimum)
    if bad:
        raise AssertionError(
            "%d text object(s) below the %.3g pt floor: %s"
            % (len(bad), minimum,
               "; ".join("%.3g pt %r" % (sz, s[:40]) for sz, s in bad[:6])))
    return True


# --- output -----------------------------------------------------------------

FIGDIR = pathlib.Path(__file__).resolve().parent / "figs"


def assert_journal_width(fig, tol=1e-3):
    """The figure is exactly one column or exactly the text width.

    The font floor only means something if the figure reaches the page at the
    size it was drawn at, which it does only when its own width is the width
    LaTeX gives it.  A driver that drifts to 3.5 in and is then scaled to
    3.375 in prints 7.7 pt type while asserting 8.
    """
    w = float(fig.get_figwidth())
    if min(abs(w - COLUMN_W), abs(w - DOUBLE_W)) > tol:
        raise AssertionError(
            "figure is %.4f in wide, which is neither the %.3f in column nor "
            "the %.3f in text width" % (w, COLUMN_W, DOUBLE_W))
    return True


def save(fig, stem, outdir=None, minimum=FONT_MIN, png=True):
    """Write `<outdir>/<stem>.pdf` (what LaTeX includes) and a PNG preview.

    Returns the list of paths written.  Asserts the font floor and the journal
    width first, so a figure that breaks T23's acceptance never reaches the
    disk.
    """
    assert_journal_width(fig)
    assert_min_fontsize(fig, minimum)
    outdir = pathlib.Path(outdir) if outdir is not None else FIGDIR
    outdir.mkdir(parents=True, exist_ok=True)
    written = []
    pdf = outdir / ("%s.pdf" % stem)
    # `CreationDate: None` omits the timestamp the PDF backend would otherwise
    # stamp into every file.  Without it a rebuild that changes no number still
    # changes all four PDFs' bytes, so "the figures are unchanged" could not be
    # answered with a checksum -- and this repository's rule for a published
    # figure stem is that it is bit for bit.  With it, `build.sh` run twice on
    # one tree leaves the PDFs, the PNGs, the JSONs and the captions identical.
    fig.savefig(pdf, metadata={"CreationDate": None})
    written.append(pdf)
    if png:
        preview = outdir / ("%s.png" % stem)
        fig.savefig(preview)
        written.append(preview)
    for p in written:
        print("wrote", p)
    return written
