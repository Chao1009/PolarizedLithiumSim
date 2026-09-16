#!/usr/bin/env bash
# paper/build.sh -- build the letter: four figures, Table 1, and the PDF.
#
# plans/07 SS 7.8 requires that ONE script regenerate the four letter figures
# and Table 1 (the run-19 ledger's T25).  This is that script, and it is the
# only supported way to rebuild any of them:
#
#   bash paper/build.sh              figures, captions, Table 1, then main.pdf
#   bash paper/build.sh --figures    figures, captions and Table 1 only
#
# What it makes, in order:
#
#   figs/fig1.pdf .. fig4.pdf   the letter figures, at the journal column
#                               widths, plus PNG previews
#   figs/figN.json              every number the driver measured
#   figs/figN_caption.tex       the caption of that figure, written from the
#                               same arrays, as a \newcommand main.tex inputs
#   table1.tex                  the projections table, generated whole from
#                               figs/fig2.json and figs/fig4.json
#   main.pdf                    the document
#
# The drivers never write into evgen/ or fastsim/: they import the producing
# scripts' functions and draw new stems under paper/figs/, so the published
# figures those scripts own stay bit for bit what they are.  Each takes a few
# seconds and a few hundred MB.
#
# TeX engine: this repository's machines have no pdflatex/xelatex/lualatex.
# The engine is tectonic, which fetches whatever packages main.tex asks for
# (elsarticle among them) from its own bundle on first use and caches them.
# Point $TECTONIC at the binary, or put `tectonic` on $PATH:
#
#   TECTONIC=/path/to/tectonic bash paper/build.sh
#
# The binary is deliberately NOT committed to this repository.  Get one from
# https://github.com/tectonic-typesetting/tectonic/releases (the
# x86_64-unknown-linux-musl tarball on a Linux box).
#
# Exit codes: 0 on a PDF, 1 on a missing engine or a failed driver, whatever
# tectonic returns otherwise.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$here"

python="${PYTHON:-python3}"

figures_only=0
if [[ "${1:-}" == "--figures" ]]; then
    figures_only=1
    shift
fi

engine=""
if [[ $figures_only -eq 0 ]]; then
    engine="${TECTONIC:-}"
    if [[ -z "$engine" ]]; then
        engine="$(command -v tectonic || true)"
    fi
    if [[ -z "$engine" ]]; then
        echo "build.sh: no TeX engine." >&2
        echo "  set \$TECTONIC to a tectonic binary, or put tectonic on \$PATH." >&2
        echo "  releases: https://github.com/tectonic-typesetting/tectonic/releases" >&2
        echo "  (bash build.sh --figures builds the figures and Table 1 without one)" >&2
        exit 1
    fi
    echo "build.sh: engine = $engine ($("$engine" --version))"
fi

# --- the four letter figures, their captions and their numbers -------------
for driver in fig1_phase_space fig2_inclusive fig3_extraction fig4_coherent; do
    echo "build.sh: $driver.py"
    "$python" "$driver.py"
done

# --- Table 1, from the numbers those drivers just wrote --------------------
echo "build.sh: table1.py"
"$python" table1.py

for f in figs/fig1.pdf figs/fig2.pdf figs/fig3.pdf figs/fig4.pdf table1.tex; do
    if [[ ! -f "$f" ]]; then
        echo "build.sh: $f was not produced" >&2
        exit 1
    fi
done

if [[ $figures_only -eq 1 ]]; then
    echo "build.sh: four figures, four captions and table1.tex are up to date"
    exit 0
fi

# --- the document ----------------------------------------------------------
"$engine" -X compile main.tex --outdir . --keep-logs --keep-intermediates

if [[ ! -f main.pdf ]]; then
    echo "build.sh: tectonic exited 0 but produced no main.pdf" >&2
    exit 1
fi
echo "build.sh: wrote $here/main.pdf"
