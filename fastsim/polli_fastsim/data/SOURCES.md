# Digitized theory curves — what they are and how to re-derive them

Every table here is a published figure read back from the PDF's own path
operators by `tools/digitize_figure.py`, not a pixel pick and not a
transcription.  All five figures are pure vector with text tick labels, so
the extraction is exact up to the drawing resolution of the curve itself
(200–600 vertices per curve).  The PDFs are the committed copies in
`refs/`, fetched with `python3 refs/find_ref.py --fetch` from the arXiv ids
in `refs/refs_dict.json`; the loader is `polli_fastsim.polarized._load_curve`.

`tools/digitize_figure.py` needs PyMuPDF, which is a **dev-time dependency
only** — install it into a scratch prefix with `pip install --target <dir>
pymupdf` and put that on `PYTHONPATH`.  Nothing at runtime reads a PDF; the
committed artefact is the CSV.

Two conventions run through the files.  The first column is always `x`, on
a uniform grid spanning the INTERSECTION of the curves' own x ranges, and
the header comment of each file records each curve's native range and its
(colour, dash) key.  Values are the plotted quantity: an EMC ratio, `x*b1`
where the paper plots `x*b1`, `b1` where it plots `100 b1` (divided by 100
at extraction, so the CSV is `b1`).

## Calibration, and the one systematic to know about

The tool takes the plot frame in PDF points and the axis values at its
edges, then maps every numeric word in the axis margins back through that
calibration and prints the residual, so a mis-set range shows up at once.
On the x axis the residuals are zero to four decimals at every figure: the
tick-label words are centred on their ticks.  On the y axis they are a
CONSTANT offset of one tick-label height fraction — +0.0047 for CBT,
+0.0016 for TMT, 0.0000 for CDKS — because a text bounding box spans the
font's ascent-to-descent box and sits about 0.9 pt above the tick it
labels, while the frame edge does not.  The frame-based calibration is
therefore the authority and the word check is the cross-check; for CBT it
is confirmed independently by the dotted reference line the paper draws at
1.0, which reads back as 1.0001.

The strongest validation is external: CDKS Fig. 5's Q² = 2.5 GeV² theory-1
sum, extracted from page 10 with its own frame and its own axis ranges,
agrees with the same curve extracted from Fig. 4 on page 9 to better than
1e-7 absolute (0.1% relative) at x = 0.05, 0.1, 0.3, 0.5 and 0.8 — two
independent figures, two independent calibrations.  Internally, CDKS's
solid curves reproduce the sum of their own dashed (SD) and dotted (DD)
curves to 1e-7 everywhere.

## `cbt_polemc_7Li_Q5.csv`

I. C. Cloët, W. Bentz, A. W. Thomas, *Spin-dependent structure functions in
nuclear matter and the polarized EMC effect*, PLB 642 (2006) 210,
arXiv:nucl-th/0605061 — **page 7, FIG. 6, upper-left panel: ⁷Li at
Q² = 5 GeV²**.  Frame (82.43, 55.92)–(296.11, 187.99) pt, x = 0…1,
ratio = 0.55…1.20 (the frame bottom is half a major tick below the 0.6
label; the major ticks are 20.32 pt apart and the paper's dotted 1.0 line
sits at y = 96.55, which this calibration returns as 1.0001).

Curve → label is read off the legend handles, not assigned by eye
(`--inspect` prints it): the handle at y = 151.0 in colour (0.182, 0.19,
0.573) dashed carries the text *Unpolarized EMC effect*; the one at
y = 163.2 in (0.93, 0.111, 0.141) solid *Polarized EMC effect:
R^(3/2 1)_As*, their Eq. (26); the one at y = 175.4 in the same red dotted
*Polarized EMC effect: R^{3/2 3/2}_As*, their Eq. (23) — the one
`plans/01_findings_physics_case.md` defines the programme's ΔR_A to be, and
the one the money plot draws.  Columns `R_unpol`, `R_pol_eq26`,
`R_pol_eq23`; covered range x = 0.0275–0.8710.

```bash
python3 tools/digitize_figure.py --pdf refs/nucl-th_0605061.pdf --page 7 \
  --frame 82.43 55.92 296.11 187.99 --xrange 0 1 --yrange 0.55 1.2 \
  --curve "R_unpol:0.182,0.19,0.573:[ 1.50769 ] 0" \
  --curve "R_pol_eq26:0.93,0.111,0.141:[] 0" \
  --curve "R_pol_eq23:0.93,0.111,0.141:[ 0 1.50769 ] 0" \
  --grid 300 --out fastsim/polli_fastsim/data/cbt_polemc_7Li_Q5.csv
```

## `tmt_polemc_nm_Q10.csv`

S. Tronchin, H. H. Matevosyan, A. W. Thomas, *Polarized EMC effect in the
QMC model*, PLB 783 (2018) 247, arXiv:1806.00481 — **page 9, Figure 4:
isospin-symmetric nuclear matter at Q² = 10 GeV²**.  Frame
(207.03, 116.67)–(403.28, 258.42) pt, x = 0…1, ratio = 0.6…1.2 (the
paper's own horizontal reference line at y = 163.94 reads back as 0.9999).
Legend handles: (0, 0, 1) solid → *Unpolarized*, (0.75, 0, 0.75) dashed →
*Polarized*.  Columns `R_unpol`, `R_pol`; covered range x = 0.0015–0.7392.

```bash
python3 tools/digitize_figure.py --pdf refs/1806.00481.pdf --page 9 \
  --frame 207.03 116.67 403.28 258.42 --xrange 0 1 --yrange 0.6 1.2 \
  --curve "R_unpol:0.0,0.0,1.0:[] 0" \
  --curve "R_pol:0.75,0.0,0.75:[ 2.4 2.4 ] 0" \
  --grid 300 --out fastsim/polli_fastsim/data/tmt_polemc_nm_Q10.csv
```

## `b1_cdks_q2p5.csv` and `b1_cdks_q2set.csv`

W. Cosyn, Yu-Bing Dong, S. Kumano, M. Sargsian, *Deuteron tensor structure
function b1*, PRD 95 (2017) 074036, arXiv:1702.05337.

`b1_cdks_q2p5.csv` — **page 9, FIG. 4: x·b1 at Q² = 2.5 GeV²**.  Frame
(371.31, 568.87)–(526.90, 668.51) pt, x = 0…1.6 (nucleon-scaled, so the
deuteron's x runs to 2), x·b1 = −0.001…0.0015 (six major ticks, the outer
two on the frame edges).  Blue (0, 0, 1) is theory 1 (their Eq. 16), green
(0.037, 0.833, 0) theory 2 (Eq. 44); the caption and the three black legend
handles at y = 644.0 / 652.2 / 659.3 give solid → SD+DD, dashed → SD,
dotted → DD.  The solid curves are drawn as 320 (160) one-segment paths, so
`--min-items 1` is needed; the black legend handles do not collide because
they are a different colour.  Covered range x = 0.0100–1.5900.

`b1_cdks_q2set.csv` — **page 10, FIG. 5: the Q² handle**, same quantity at
Q² = 1.0, 2.5 and 5.0 GeV² for both theories.  Frame
(106.07, 82.89)–(263.60, 183.77) pt, x = 0…1.6, x·b1 = −0.002…0.004.  Here
the Q² values are keyed by DASH, not colour: the black legend handles at
y = 160.4 / 168.1 / 175.6 carry *Q²=1.0 GeV²* (dashed), *Q²=2.5 GeV²*
(solid) and *Q²=5.0 GeV²* (dotted).  Covered range x = 0.0100–1.5850.

```bash
python3 tools/digitize_figure.py --pdf refs/1702.05337.pdf --page 9 \
  --frame 371.31 568.87 526.90 668.51 --xrange 0 1.6 --yrange -0.001 0.0015 \
  --min-items 1 \
  --curve "xb1_theory1_sum:0.0,0.0,1.0:[] 0" \
  --curve "xb1_theory1_SD:0.0,0.0,1.0:[ 2.84259 2.84259 ] 0" \
  --curve "xb1_theory1_DD:0.0,0.0,1.0:[ 1.13704 2.27407 ] 0" \
  --curve "xb1_theory2_sum:0.037,0.833,0.0:[] 0" \
  --curve "xb1_theory2_SD:0.037,0.833,0.0:[ 2.84259 2.84259 ] 0" \
  --curve "xb1_theory2_DD:0.037,0.833,0.0:[ 1.13704 2.27407 ] 0" \
  --grid 300 --out fastsim/polli_fastsim/data/b1_cdks_q2p5.csv

python3 tools/digitize_figure.py --pdf refs/1702.05337.pdf --page 10 \
  --frame 106.07 82.89 263.60 183.77 --xrange 0 1.6 --yrange -0.002 0.004 \
  --min-items 1 \
  --curve "xb1_theory1_q2_1p0:0.0,0.0,1.0:[ 2.8781203 2.8781203 ] 0" \
  --curve "xb1_theory1_q2_2p5:0.0,0.0,1.0:[] 0" \
  --curve "xb1_theory1_q2_5p0:0.0,0.0,1.0:[ 1.15125 2.3025 ] 0" \
  --curve "xb1_theory2_q2_1p0:0.037,0.833,0.0:[ 2.8781203 2.8781203 ] 0" \
  --curve "xb1_theory2_q2_2p5:0.037,0.833,0.0:[] 0" \
  --curve "xb1_theory2_q2_5p0:0.037,0.833,0.0:[ 1.15125 2.3025 ] 0" \
  --grid 300 --out fastsim/polli_fastsim/data/b1_cdks_q2set.csv
```

## `b1_miller.csv` and `b1_miller_q2set.csv`

G. A. Miller, *Pionic and hidden-color, six-quark contributions to the
deuteron b1 structure function*, PRC 89 (2014) 045203, arXiv:1311.4561.

These two are Mathematica plots with only a left and a bottom axis — no
closed frame — so `--inspect` finds no frame candidate and the box was read
off the axis lines and their major tick marks instead.  The y axis is
`100 b1(x)`, and the `--yrange` given below is already divided by 100, so
the CSV holds b1.  The tool's word check therefore prints the printed tick
labels against values a hundred times smaller; that is the intended
conversion, not a calibration error.

`b1_miller.csv` — **page 10, FIG. 5: b1 = b1(pion) + b1(6q)**, the total
that reproduces HERMES with a hidden-colour probability of 0.15%.  Axis
lines at x = 164.88 and y = 291.61 pt; major x ticks at 164.88, 226.47,
288.06, 349.65, 411.24 for x = 0, 0.2, 0.4, 0.6, 0.8 (61.59 pt per 0.2, so
the right edge of the extraction box at 447.81 is x = 0.918753); major y
ticks at 237.18, 182.76, 128.33 for 100 b1 = 5, 10, 15 (54.42 pt per 5,
which puts 100 b1 = 0 at y = 291.60, i.e. on the x axis, to 0.01 pt).  One
curve, colour (0.247, 0.24, 0.6), 578 vertices; covered range
x = 0.0100–0.9000.

`b1_miller_q2set.csv` — **page 11, FIG. 6: 100 b1 at Q² = 1.17, 1.76, 2.12
and 3.25 GeV²**.  Axis lines at x = 167.02 and y = 268.00 pt; major x ticks
45.845 pt per 0.1 with 212.86 ↔ 0.2, major y ticks 30.895 pt per 0.5 with
268.00 ↔ 0.  Six colour groups are present but only FIVE distinct traces:
(0.6, 0.547, 0.24) and (0.24, 0.6, 0.337) are the same 495-vertex path
drawn twice, one colour over the other, so the duplicate is not extracted.
This figure carries no legend, so the curve → label map comes from the
caption: *"for values of Q² = 1.17, 1.76, 2.12 and 3.25 GeV² [29]
distributions and for [35] (lowest curve at x = 0.15).  For the other
curves, b1(pion) increases as Q² increases for small values of x."*  At
x = 0.15 the five traces read 0.00993, 0.01074, 0.01155, 0.01197, 0.01301,
so the lowest is the [35] curve and the remaining four are Q² = 1.17, 1.76,
2.12, 3.25 in that order.  Columns are named accordingly; the assignment is
an inference from the caption, unlike every other file here, where it comes
off a legend handle.  Covered range x = 0.1007–0.7000.

```bash
python3 tools/digitize_figure.py --pdf refs/1311.4561.pdf --page 10 \
  --frame 164.88 95.67 447.81 305.99 --xrange 0 0.918753 \
  --yrange -0.0132213 0.1800165 --min-items 100 \
  --curve "b1:0.247,0.24,0.6:[] 0" \
  --grid 300 --out fastsim/polli_fastsim/data/b1_miller.csv

python3 tools/digitize_figure.py --pdf refs/1311.4561.pdf --page 11 \
  --frame 167.02 95.76 447.81 305.77 --xrange 0.1 0.712473 \
  --yrange -0.0061126 0.027875 --min-items 100 \
  --curve "b1_ref35:0.6,0.24,0.563:[] 0" \
  --curve "b1_q2_1p17:0.24,0.353,0.6:[] 0" \
  --curve "b1_q2_1p76:0.6,0.547,0.24:[] 0" \
  --curve "b1_q2_2p12:0.6,0.24,0.443:[] 0" \
  --curve "b1_q2_3p25:0.247,0.24,0.6:[] 0" \
  --grid 300 --out fastsim/polli_fastsim/data/b1_miller_q2set.csv
```

## Not digitized

Wang–Bentz–Cloët–Thomas, *Polarized gluon EMC effect*, J. Phys. G 49 (2022)
03LT01, arXiv:2109.03591, Figure 3 (page 8) carries g₁A/g₁p and Δg_A/Δg_p —
the dg₁/dlnQ² observable of `plans/02` step 1.2.2, which has no money plot
yet.  The PDF is in `refs/` and the entry in `refs_dict.json` records where
the curves are.

Cosyn–Weiss arXiv:2603.23700 FIG. 13 (page 36) is the tagged tensor
asymmetry A_T∥ in light-front variables (α_p, p_pT) the generator does not
carry, so it is a comparison rather than a drop-in.  It did not need
digitizing: page 35 gives the closed form (Eq. 6.12), its extrema
(Eqs. 6.13–6.14) and Table II, which is what
`evgen/tests/test_tagged.py::test_cosyn_weiss_tensor_gate` is pinned
against.
