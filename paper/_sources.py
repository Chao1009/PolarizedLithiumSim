#!/usr/bin/env python3
"""Plumbing shared by the four letter-figure drivers: where the machinery
lives, what arguments the published commands run with, and where the numbers
a driver measures are written down.

There is no physics in this file and none in the drivers.  A letter figure is
a re-plot: the driver imports the functions of the producing script named in
`figs/README.md`, calls them in the order that script calls them, and draws
the arrays that come back at the journal's column width.  Nothing is
recomputed by a second implementation, which is the only way the letter and
the reports can be guaranteed to carry the same numbers.

`published_args` is the part worth reading.  A driver must run at the
producing script's published defaults -- the same seed, the same beam
configuration, the same Delta model, the same luminosities -- or its figure
is a different measurement wearing the same caption.  Restating those
defaults here would be a copy that drifts.  Instead the producing script's own
`main()` is entered and stopped at its `argparse` call, and the namespace it
was about to receive is handed to the driver.  The published command is
therefore executed, as far as the arguments; the driver cannot hold a default
the script no longer has, and a flag added to the script arrives here on its
own.
"""

import argparse
import decimal
import json
import pathlib
import re
import sys

PAPER = pathlib.Path(__file__).resolve().parent
ROOT = PAPER.parent
EVGEN = ROOT / "evgen"
SCRIPTS = EVGEN / "scripts"
FASTSIM = ROOT / "fastsim"
FIGDIR = PAPER / "figs"

# evgen/scripts first (the producing scripts import each other by bare name),
# then evgen (for `polligen`, whose __init__ puts fastsim/ on the path).
for _p in (str(SCRIPTS), str(EVGEN), str(FASTSIM)):
    if _p not in sys.path:
        sys.path.insert(0, _p)


class _StopAtArgs(Exception):
    """Raised inside the producing script's main() once its arguments are
    parsed, to stop it before it does any work."""


def published_args(main, argv=()):
    """The namespace the producing script's own `main()` would parse.

    `main` is that function; `argv` is the published command line, which for
    every figure of this letter is the empty one -- the four published PNGs
    are all produced by a bare `python3 scripts/<name>.py --outdir .`, so the
    published settings ARE the parser's defaults, and reading them off the
    parser is what keeps this letter pinned to them.
    """
    real = argparse.ArgumentParser.parse_args
    captured = {}

    def spy(self, args=None, namespace=None):
        captured["ns"] = real(self, list(argv), namespace)
        raise _StopAtArgs()

    argparse.ArgumentParser.parse_args = spy
    try:
        main()
    except _StopAtArgs:
        pass
    finally:
        argparse.ArgumentParser.parse_args = real
    if "ns" not in captured:
        raise RuntimeError(
            "%s.%s returned without parsing arguments: the producing script "
            "no longer starts with argparse, and the driver's settings can "
            "no longer be read off it" % (main.__module__, main.__name__))
    return captured["ns"]


def describe(args):
    """The published settings a caption or a log line needs to state."""
    return {k: (v if isinstance(v, (int, float, str, bool, type(None)))
                else str(v))
            for k, v in sorted(vars(args).items())}


# --- the numbers a driver measures -----------------------------------------

def dump(stem, payload):
    """Write `figs/<stem>.json`: everything the driver measured.

    Table 1 is assembled from these files and from nothing else, so a number
    in the table and a number in the figure beside it are the same draw of the
    same pseudo-experiment.
    """
    FIGDIR.mkdir(parents=True, exist_ok=True)
    path = FIGDIR / ("%s.json" % stem)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
    print("wrote", path)
    return path


# --- formatting shared by the captions and by Table 1 ----------------------

def dec(v, nd):
    """`v` to `nd` decimals, rounding half away from zero.

    Python's %-formatting rounds half to even through the binary value, so
    0.0585 prints as 0.058 while the documents print 0.059.  The value is
    first read at ten significant digits -- far beyond any precision this
    programme carries -- so that floating-point dust (0.0585 arrives as
    0.058499999999999996) cannot decide the last printed digit.
    """
    q = decimal.Decimal(1).scaleb(-nd)
    d = decimal.Decimal("%.10g" % float(v)).quantize(
        q, rounding=decimal.ROUND_HALF_UP)
    return "%s" % d


def sci(v, sig=3):
    """`1.9e8` -> `1.90\\times10^{8}` (math mode, no dollars)."""
    import math
    if v == 0:
        return "0"
    exp = int(math.floor(math.log10(abs(v))))
    mant = v / 10.0 ** exp
    return r"%.*f\times10^{%d}" % (max(sig - 1, 0), mant, exp)


def sig3(v):
    """Three significant figures, fixed point, half away from zero."""
    import math
    if v == 0:
        return "0"
    exp = int(math.floor(math.log10(abs(v))))
    return dec(v, max(2 - exp, 0))


_BEAM = re.compile(r"e\(([\d.]+)\)\s*x\s*(\d+)([A-Za-z]+)\(([\d.]+)/u\)")


def beam_math(label):
    """`e(10) x 6Li(99.5/u)` -> the same thing in math mode.

    The label is `beams.BeamConfig.label()`; it is parsed rather than
    restated, so a change of beam reaches every caption on its own.
    """
    m = _BEAM.match(label.strip())
    if m is None:
        return r"\texttt{%s}" % label
    e, a, sym, pu = m.groups()
    return (r"$e(%s\,\mathrm{GeV}) \times {}^{%s}\mathrm{%s}"
            r"\,(%s\,\mathrm{GeV/u})$" % (e, a, sym, pu))


def write_caption(stem, macro, text):
    """Write `figs/<stem>_caption.tex`, defining `\<macro>` as `text`.

    The letter's figure captions quote numbers -- yields, cuts, acceptances --
    and a number typed into `main.tex` by hand is a number that drifts from
    the figure beside it the next time a driver runs.  Each driver therefore
    emits its own caption from the same arrays it plots, `main.tex` inputs the
    four files in its preamble and writes `\caption{\figonecaption}`, and
    plans/07 SS 7.8's "regenerated by one build script" covers the captions as
    well as the figures.
    """
    FIGDIR.mkdir(parents=True, exist_ok=True)
    path = FIGDIR / ("%s_caption.tex" % stem)
    body = " ".join(text.split())
    path.write_text(
        "%% paper/figs/%s_caption.tex -- GENERATED by the %s driver.\n"
        "%% Do not edit: rebuild with paper/build.sh.\n"
        "\\newcommand{\\%s}{%%\n%s%%\n}\n" % (stem, stem, macro, body),
        encoding="utf-8")
    print("wrote", path)
    return path


def load(stem):
    path = FIGDIR / ("%s.json" % stem)
    if not path.exists():
        raise SystemExit(
            "%s does not exist: run the driver that writes it "
            "(paper/build.sh runs all four) before building Table 1." % path)
    return json.loads(path.read_text(encoding="utf-8"))
