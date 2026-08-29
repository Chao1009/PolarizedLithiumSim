#!/usr/bin/env python3
"""Extract curves from a VECTOR figure in a published PDF into a CSV table.

Every target figure of plans/02 step 1.2 (Cloet-Bentz-Thomas PLB 642:210
Fig. 6, Tronchin-Matevosyan-Thomas PLB 783:247 Fig. 4, Cosyn-Dong-Kumano-
Sargsian PRD 95:074036 Figs. 4/5, Miller PRC 89:045203 Figs. 5/6) is drawn
as PDF path operators, not as a raster.  So the published curve can be read
back exactly -- no pixel picking, no eyeballing -- by walking the drawing
operators of the page and mapping the plot-frame coordinates onto the axis
values printed on the tick labels.

    # 1. see what is inside the frame and how the legend labels it
    python3 tools/digitize_figure.py --pdf refs/nucl-th_0605061.pdf --page 7 \
        --frame 82.43 55.92 296.11 187.99 --inspect

    # 2. extract, one column per --curve, on a common x grid
    python3 tools/digitize_figure.py --pdf refs/nucl-th_0605061.pdf --page 7 \
        --frame 82.43 55.92 296.11 187.99 --xrange 0 1 --yrange 0.55 1.2 \
        --curve "R_unpol:0.182,0.19,0.573:[ 1.50769 ] 0" \
        --out fastsim/polli_fastsim/data/cbt_polemc_7Li_Q5.csv

`--frame` is the plot box in PDF points (`--inspect` prints the candidate
boxes it finds on the page) and `--xrange`/`--yrange` are the axis values at
its edges.  The calibration is CHECKED, not assumed: every numeric word in
the axis margins is mapped through it and printed next to its printed value,
so a mis-set range shows up at once as a column of non-round numbers.  The
y-axis words carry a known constant offset -- a text bounding box is centred
on the font's ascent-to-descent box, roughly 0.9 pt above the tick it labels,
while the x-axis words are centred horizontally and show no bias.

PyMuPDF is needed and is a DEV-TIME dependency only: what the repository
uses at runtime is the committed CSV, never this script.  Install it into a
scratch prefix (`pip install --target <dir> pymupdf`) and put that on
PYTHONPATH; nothing is added to the fast-sim or generator requirements.
Provenance for every committed table -- paper, page, figure, frame box, axis
ranges, curve keys and the exact command -- is recorded in
`fastsim/polli_fastsim/data/SOURCES.md`.
"""

import argparse
import pathlib
import re
import sys

import numpy as np

_NUM = re.compile(r"^[-−+]?\d+(?:\.\d+)?$")


def _open(pdf, page):
    try:
        import pymupdf
    except ImportError:  # pragma: no cover - dev-time tool
        sys.exit("PyMuPDF is required: pip install --target <dir> pymupdf")
    return pymupdf.open(pdf)[page - 1]


def _key(drawing):
    """(colour, dash, width) identity of a stroked path."""
    col = drawing.get("color")
    return (tuple(round(c, 3) for c in col) if col else None,
            drawing.get("dashes"), round(drawing.get("width") or 0.0, 2))


def _points(drawing):
    """Vertices of a path IN DRAWING ORDER (never sorted: several published
    curves double back -- the Fermi rise, a clipped re-entry -- and sorting
    by x would silently splice them)."""
    pts = []
    for it in drawing["items"]:
        if it[0] == "l":
            pts.append((it[1].x, it[1].y))
            pts.append((it[2].x, it[2].y))
        elif it[0] == "c":
            pts.append((it[1].x, it[1].y))
            pts.append((it[4].x, it[4].y))
        elif it[0] == "re":
            r = it[1]
            pts.append((r.x0, r.y0))
    return pts


def _overlaps(rect, frame, pad=1.0):
    x0, y0, x1, y1 = frame
    return (rect.x1 >= x0 - pad and rect.x0 <= x1 + pad
            and rect.y1 >= y0 - pad and rect.y0 <= y1 + pad)


def frames(page, min_side=40.0):
    """Axis-frame candidates: the closed boxes drawn as four long strokes."""
    hori, vert = {}, {}
    for d in page.get_drawings():
        for it in d["items"]:
            if it[0] != "l":
                continue
            (ax, ay), (bx, by) = (it[1].x, it[1].y), (it[2].x, it[2].y)
            if abs(ay - by) < 0.01 and abs(bx - ax) > min_side:
                hori.setdefault(round(ay, 2), set()).add(
                    (round(min(ax, bx), 2), round(max(ax, bx), 2)))
            elif abs(ax - bx) < 0.01 and abs(by - ay) > min_side:
                vert.setdefault(round(ax, 2), set()).add(
                    (round(min(ay, by), 2), round(max(ay, by), 2)))
    out = []
    for y0 in sorted(hori):
        for y1 in sorted(hori):
            if y1 - y0 < min_side:
                continue
            for (x0, x1) in sorted(hori[y0] & hori[y1]):
                if (y0, y1) in vert.get(x0, ()) and (y0, y1) in vert.get(x1, ()):
                    out.append((x0, y0, x1, y1))
    return out


def axis_words(page, frame, margin=45.0):
    """Numeric tick labels in the two margins, as (value, centre_pt)."""
    x0, y0, x1, y1 = frame
    below, left = [], []
    for wx0, wy0, wx1, wy1, txt, *_ in page.get_text("words"):
        if not _NUM.match(txt):
            continue
        val = float(txt.replace("−", "-"))
        cx, cy = 0.5 * (wx0 + wx1), 0.5 * (wy0 + wy1)
        if y1 < wy0 < y1 + margin and x0 - margin < cx < x1 + margin:
            below.append((val, cx))
        elif x0 - margin < wx1 < x0 and y0 - margin < cy < y1 + margin:
            left.append((val, cy))
    return below, left


def report_calibration(name, words, lo_pt, hi_pt, lo_val, hi_val, log=False):
    """Map each printed tick label through the calibration and show it."""
    if not words:
        print("  %s: no numeric tick labels found" % name)
        return
    print("  %s: printed -> read back through the calibration" % name)
    for val, pt in sorted(words, key=lambda vp: vp[1]):
        got = _to_value(np.array([pt]), lo_pt, hi_pt, lo_val, hi_val, log)[0]
        print("      %-8s %.4f   (%+.4f)" % (val, got, got - val))


def _to_value(pts, lo_pt, hi_pt, lo_val, hi_val, log=False):
    f = (np.asarray(pts, dtype=float) - lo_pt) / (hi_pt - lo_pt)
    if log:
        return 10.0 ** (np.log10(lo_val) + f * (np.log10(hi_val)
                                                - np.log10(lo_val)))
    return lo_val + f * (hi_val - lo_val)


def collect(page, frame, colour, dash, min_items, edge_tol):
    """All in-frame vertices of the paths matching (colour, dash), in order,
    with the frame-edge (clipped) vertices dropped."""
    x0, y0, x1, y1 = frame
    pts = []
    for d in page.get_drawings():
        k = _key(d)
        if colour is not None and k[0] != colour:
            continue
        if dash is not None and k[1] != dash:
            continue
        if len(d["items"]) < min_items or not _overlaps(d["rect"], frame):
            continue
        keep = []
        for px, py in _points(d):
            if not (x0 <= px <= x1 and y0 <= py <= y1):
                continue                      # outside the panel entirely
            if (abs(px - x0) < edge_tol or abs(px - x1) < edge_tol
                    or abs(py - y0) < edge_tol or abs(py - y1) < edge_tol):
                continue                      # clipped at the frame
            keep.append((px, py))
        if len(keep) >= min_items:            # a curve of THIS panel
            pts.extend(keep)
    return pts


def legend_label(page, rect, frame):
    """Text to the right of a legend handle, on the handle's own line."""
    cx, cy = 0.5 * (rect.x0 + rect.x1), 0.5 * (rect.y0 + rect.y1)
    near = []
    for wx0, wy0, wx1, wy1, txt, *_ in page.get_text("words"):
        wcy = 0.5 * (wy0 + wy1)
        if abs(wcy - cy) < 5.0 and wx0 > cx and wx0 - rect.x1 < 60.0:
            near.append((wx0, txt))
    return " ".join(t for _, t in sorted(near))


def inspect(page, frame):
    """List the (colour, dash, width) path groups inside the frame, with the
    legend text of every short path -- this is where the curve -> label map
    comes from: it is read off the legend handles, not assigned by eye."""
    groups = {}
    for d in page.get_drawings():
        if frame and not _overlaps(d["rect"], frame):
            continue
        groups.setdefault(_key(d), []).append(d)
    print("frame %s" % (frame,))
    for k, ds in sorted(groups.items(), key=lambda kv: -sum(len(d["items"])
                                                            for d in kv[1])):
        nitem = sum(len(d["items"]) for d in ds)
        print("  colour=%s dash=%r width=%.2f  paths=%d items=%d"
              % (k[0], k[1], k[2], len(ds), nitem))
        for d in ds:
            if len(d["items"]) <= 4:
                lab = legend_label(page, d["rect"], frame)
                if lab:
                    print("      legend handle at (%.1f, %.1f) -> %s"
                          % (d["rect"].x0, d["rect"].y0, lab))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdf", required=True)
    ap.add_argument("--page", type=int, required=True, help="1-based")
    ap.add_argument("--frame", nargs=4, type=float,
                    metavar=("X0", "Y0", "X1", "Y1"))
    ap.add_argument("--xrange", nargs=2, type=float, metavar=("LEFT", "RIGHT"))
    ap.add_argument("--yrange", nargs=2, type=float, metavar=("BOTTOM", "TOP"))
    ap.add_argument("--xlog", action="store_true")
    ap.add_argument("--ylog", action="store_true")
    ap.add_argument("--curve", action="append", default=[],
                    metavar="NAME:R,G,B[:DASH]",
                    help="one output column; DASH is the PDF dash string")
    ap.add_argument("--min-items", type=int, default=10,
                    help="skip paths shorter than this (legend handles)")
    ap.add_argument("--edge-tol", type=float, default=0.25,
                    help="drop vertices this close to the frame (clipped)")
    ap.add_argument("--grid", type=int, default=200,
                    help="points of the common output grid")
    ap.add_argument("--out")
    ap.add_argument("--inspect", action="store_true",
                    help="list frames, path groups and legend labels; no CSV")
    args = ap.parse_args(argv)

    page = _open(args.pdf, args.page)
    if args.inspect:
        if not args.frame:
            print("axis-frame candidates on page %d:" % args.page)
            for f in frames(page):
                print("  --frame %.2f %.2f %.2f %.2f" % f)
            return 0
        inspect(page, tuple(args.frame))
        return 0

    frame = tuple(args.frame)
    x0, y0, x1, y1 = frame
    below, left = axis_words(page, frame)
    print("calibration check (%s p.%d)" % (args.pdf, args.page))
    report_calibration("x", below, x0, x1, *args.xrange, log=args.xlog)
    report_calibration("y", left, y1, y0, *args.yrange, log=args.ylog)

    curves, meta = {}, {}
    for spec in args.curve:
        name, colour, dash = (spec.split(":", 2) + [None])[:3] \
            if spec.count(":") >= 2 else (spec.split(":") + [None])[:3]
        col = tuple(round(float(v), 3) for v in colour.split(",")) if colour else None
        pts = collect(page, frame, col, dash, args.min_items, args.edge_tol)
        if not pts:
            sys.exit("no path matched %r" % spec)
        px = np.array([p[0] for p in pts])
        py = np.array([p[1] for p in pts])
        xv = _to_value(px, x0, x1, *args.xrange, log=args.xlog)
        yv = _to_value(py, y1, y0, *args.yrange, log=args.ylog)
        order = np.argsort(xv, kind="stable")
        xs, ys = xv[order], yv[order]
        keep = np.concatenate(([True], np.diff(xs) > 0))
        curves[name] = (xs[keep], ys[keep])
        meta[name] = "colour=%s dash=%r points=%d x=[%.4f, %.4f]" % (
            col, dash, keep.sum(), xs[0], xs[-1])
        print("  %-14s %s" % (name, meta[name]))

    lo = max(c[0][0] for c in curves.values())
    hi = min(c[0][-1] for c in curves.values())
    grid = (np.logspace(np.log10(lo), np.log10(hi), args.grid) if args.xlog
            else np.linspace(lo, hi, args.grid))
    names = list(curves)
    cols = [np.interp(grid, *curves[n]) for n in names]
    head = ["source: %s page %d, frame %s" % (args.pdf, args.page, frame),
            "axes: x %s%s, y %s%s" % (args.xrange, " (log)" if args.xlog else "",
                                      args.yrange, " (log)" if args.ylog else ""),
            "extracted by tools/digitize_figure.py; provenance in SOURCES.md"]
    head += ["curve %s: %s" % (n, meta[n]) for n in names]
    if args.out:
        out = pathlib.Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(out, np.column_stack([grid] + cols), delimiter=",",
                   fmt="%.6g", header="\n".join(head) + "\nx," + ",".join(names))
        print("wrote %s (%d rows, x = %.4f .. %.4f)" % (out, len(grid), lo, hi))
    return 0


if __name__ == "__main__":
    sys.exit(main())
