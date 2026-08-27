#!/usr/bin/env python3
"""Assemble the self-contained reports/ pages from their templates.

Two build steps per page:

1. Figures: __TAG__ placeholders are replaced with base64-embedded PNGs
   produced by the evgen scripts (run those first; evgen/README.md).
2. Display math: every <img class="tex" data-tex="..."> element is
   typeset with matplotlib mathtext (STIX fonts, glyphs converted to
   paths) and embedded as an SVG data URI, sized in pt so screen and the
   A4 print layout agree.  No JavaScript, no external fonts -- the
   output stays fully self-contained.  Inline math stays as HTML.

Pages built (template -> html), in report-number order:
  0  polarized_li_primer.template.html          educational physics primer
  1  cos2phi_money_plots_report.template.html   the projected measurements
  2  reconstruction_chain_report.template.html  the reconstruction chain
  3  eic_epic_reference.template.html           the machine and the detector
  4  nanowire_far_forward.template.html         the near-beam detector study

With --pdf each page is also print-rendered through headless
Edge/Chrome (the templates carry the @page A4 setup).

Usage:  python reports/build_report.py [--pdf]
"""

import argparse
import base64
import glob
import io
import os
import pathlib
import re
import shutil
import subprocess
import sys

REPORTS = pathlib.Path(__file__).resolve().parent
REPO = REPORTS.parent

MONEY_FIGS = {
    "__PS__": "evgen/phase_space_bins_6Li.png",
    "__M5__": "evgen/money_cos2phi_6Li.png",
    "__M7__": "evgen/money_delta_extracted_6Li.png",
    "__M6__": "evgen/money_cos2phi_coherent_6Li.png",
}

EIC_FIGS = {
    "__EIC1__": "evgen/eic_energy_menu.png",
    "__EIC2__": "evgen/eic_beam_divergence.png",
}

NANOWIRE_FIGS = {
    "__NB1__": "evgen/nearbeam_aperture_6Li.png",
    "__NB2__": "evgen/nearbeam_reach_gain_6Li.png",
    "__NB3__": "evgen/nearbeam_zid_threshold.png",
    "__NB4__": "evgen/nearbeam_zid_power.png",
}

RECO_FIGS = {
    "__RC1__": "evgen/reco_chain_inclusive_6Li.png",
    "__RC2__": "evgen/reco_chain_coherent_6Li.png",
    "__RC3__": "evgen/money_cos2phi_reco_6Li.png",
    "__RC4__": "evgen/money_delta_extracted_reco_6Li.png",
    "__RC5__": "evgen/money_cos2phi_coherent_reco_6Li.png",
    "__RC6__": "evgen/hfs_resolution_6Li.png",
}

# In report-number order, which is also the order reports/index.html lists
# them and the order they are meant to be read: primer -> projections ->
# reconstruction -> detector study.
PAGES = (
    {"stem": "polarized_li_primer", "number": 0,
     "figures": {k: MONEY_FIGS[k] for k in ("__PS__", "__M5__")}},
    {"stem": "cos2phi_money_plots_report", "number": 1, "figures": MONEY_FIGS},
    {"stem": "reconstruction_chain_report", "number": 2, "figures": RECO_FIGS},
    {"stem": "eic_epic_reference", "number": 3, "figures": EIC_FIGS},
    {"stem": "nanowire_far_forward", "number": 4, "figures": NANOWIRE_FIGS},
)

BROWSERS = (
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    "msedge", "chrome", "chromium", "google-chrome",
    # Linux boxes that carry no system browser: playwright downloads one
    # into the user cache without root (pip install playwright &&
    # python3 -m playwright install chromium), which is how the PDFs are
    # built on the analysis machine.
    *sorted(glob.glob(os.path.expanduser(
        "~/.cache/ms-playwright/chromium-*/chrome-linux64/chrome"))),
)


INK = "#1A2530"  # report body ink; math is typeset in the same color
TEX_RE = re.compile(r'<img class="tex(?P<cls>[^"]*)" data-tex="(?P<tex>[^"]+)"\s*/?>')


def render_tex(tex, fontsize=10.5):
    """Typeset one display formula with matplotlib mathtext -> (data URI,
    width pt, height pt).  svg.fonttype=path outlines every glyph, so the
    result renders identically everywhere with no font dependencies."""
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams["mathtext.fontset"] = "stix"
    matplotlib.rcParams["svg.fonttype"] = "path"
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0.0, 0.0, "$%s$" % tex, fontsize=fontsize, color=INK)
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", bbox_inches="tight", pad_inches=0.015,
                transparent=True)
    plt.close(fig)
    svg = buf.getvalue().decode("utf-8")
    m = re.search(r'width="([\d.]+)pt" height="([\d.]+)pt"', svg)
    if not m:
        sys.exit("could not read SVG size for: %s" % tex)
    uri = ("data:image/svg+xml;base64,"
           + base64.b64encode(svg.encode("utf-8")).decode("ascii"))
    return uri, float(m.group(1)), float(m.group(2))


def typeset_math(html):
    """Replace every <img class="tex..." data-tex="..."> with the rendered
    SVG, carrying the width in pt so print and screen sizes agree."""
    def sub(m):
        tex = m.group("tex").replace("&amp;", "&").replace("&quot;", '"')
        uri, w, h = render_tex(tex)
        return ('<img class="tex%s" src="%s" style="width:%.1fpt" '
                'alt="%s">' % (m.group("cls"), uri, w, m.group("tex")))
    html, n = TEX_RE.subn(sub, html)
    print("typeset %d display formulas" % n)
    return html


def build_html(page):
    template = REPORTS / ("%s.template.html" % page["stem"])
    out_html = REPORTS / ("%s.html" % page["stem"])
    html = template.read_text(encoding="utf-8")
    html = typeset_math(html)
    for tag, rel in page["figures"].items():
        png = REPO / rel
        if not png.is_file():
            sys.exit("missing figure %s (for %s) -- run the evgen "
                     "scripts first" % (png, tag))
        if tag not in html:
            sys.exit("placeholder %s not found in %s" % (tag, template))
        html = html.replace(tag, base64.b64encode(png.read_bytes())
                            .decode("ascii"))
    leftover = re.findall(r"__[A-Z0-9]+__", html)
    if leftover:
        sys.exit("unreplaced placeholders in %s: %s" % (template, leftover))
    out_html.write_text(html, encoding="utf-8")
    print("wrote %s (%.1f MB)" % (out_html, out_html.stat().st_size / 1e6))
    return out_html


def build_pdf(out_html):
    exe = next((b for b in BROWSERS
                if pathlib.Path(b).is_file() or shutil.which(b)), None)
    if exe is None:
        sys.exit("no headless-capable browser found for --pdf")
    out_pdf = out_html.with_suffix(".pdf")
    subprocess.run(
        [exe, "--headless", "--disable-gpu", "--no-sandbox",
         "--no-pdf-header-footer",
         "--print-to-pdf=%s" % out_pdf, out_html.as_uri()],
        check=True)
    print("wrote %s (%.1f MB)" % (out_pdf, out_pdf.stat().st_size / 1e6))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", action="store_true",
                    help="also print-render the PDFs (headless Edge/Chrome)")
    args = ap.parse_args()
    for page in PAGES:
        out = build_html(page)
        if args.pdf:
            build_pdf(out)
