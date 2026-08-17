#!/usr/bin/env python3
"""Assemble the self-contained cos2phi report from its template.

Replaces the __TAG__ placeholders in
reports/cos2phi_money_plots_report.template.html with base64-embedded
PNGs produced by the evgen scripts (run those first; evgen/README.md):

    __PS__  evgen/phase_space_bins_6Li.png        phase space & binning
    __M5__  evgen/money_cos2phi_6Li.png           money plot 5
    __M7__  evgen/money_delta_extracted_6Li.png   money plot 7
    __M6__  evgen/money_cos2phi_coherent_6Li.png  money plot 6

and writes reports/cos2phi_money_plots_report.html.  With --pdf it also
print-renders the PDF through headless Edge/Chrome (the template carries
the @page A4 setup).

Usage:  python reports/build_report.py [--pdf]
"""

import argparse
import base64
import pathlib
import shutil
import subprocess
import sys

REPORTS = pathlib.Path(__file__).resolve().parent
REPO = REPORTS.parent
TEMPLATE = REPORTS / "cos2phi_money_plots_report.template.html"
OUT_HTML = REPORTS / "cos2phi_money_plots_report.html"
OUT_PDF = REPORTS / "cos2phi_money_plots_report.pdf"

FIGURES = {
    "__PS__": "evgen/phase_space_bins_6Li.png",
    "__M5__": "evgen/money_cos2phi_6Li.png",
    "__M7__": "evgen/money_delta_extracted_6Li.png",
    "__M6__": "evgen/money_cos2phi_coherent_6Li.png",
}

BROWSERS = (
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    "msedge", "chrome", "chromium", "google-chrome",
)


def build_html():
    html = TEMPLATE.read_text(encoding="utf-8")
    for tag, rel in FIGURES.items():
        png = REPO / rel
        if not png.is_file():
            sys.exit("missing figure %s (for %s) -- run the evgen "
                     "scripts first" % (png, tag))
        if tag not in html:
            sys.exit("placeholder %s not found in template" % tag)
        html = html.replace(tag, base64.b64encode(png.read_bytes())
                            .decode("ascii"))
    leftover = [t for t in ("__PS__", "__M5__", "__M6__", "__M7__")
                if t in html]
    if leftover:
        sys.exit("unreplaced placeholders: %s" % leftover)
    OUT_HTML.write_text(html, encoding="utf-8")
    print("wrote %s (%.1f MB)" % (OUT_HTML,
                                  OUT_HTML.stat().st_size / 1e6))


def build_pdf():
    exe = next((b for b in BROWSERS
                if pathlib.Path(b).is_file() or shutil.which(b)), None)
    if exe is None:
        sys.exit("no headless-capable browser found for --pdf")
    subprocess.run(
        [exe, "--headless", "--disable-gpu", "--no-pdf-header-footer",
         "--print-to-pdf=%s" % OUT_PDF, OUT_HTML.as_uri()],
        check=True)
    print("wrote %s (%.1f MB)" % (OUT_PDF, OUT_PDF.stat().st_size / 1e6))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", action="store_true",
                    help="also print-render the PDF (headless Edge/Chrome)")
    args = ap.parse_args()
    build_html()
    if args.pdf:
        build_pdf()
