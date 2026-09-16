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

## `wbct_emc_nm_Q5.csv` and `wbct_polemc_nm_Q5.csv`

X. G. Wang, W. Bentz, I. C. Cloët, A. W. Thomas, *Polarized gluon EMC
effect*, J. Phys. G 49 (2022) 03LT01, arXiv:2109.03591 — **page 8,
Figure 3: isospin-symmetric nuclear matter at Q² = 5 GeV²**, the gluon
sector of `plans/02` step 1.2.2.  The two panels sit on separate frames
with identical axes, x = 0…1 and ratio = 0.6…1.2.  Each drawn curve stops
where it leaves the top of the box, at a different x for each; the
tabulated grid is the range the curves of a panel have in COMMON (the
per-curve spans are in the CSV header comments, the solid ones running on
to x = 0.911), so the last rows carry every column of the panel rather
than one of them alone.  At the other end the paper does not draw the
small-x region at all — "because of the neglect of shadowing effects …
this region is not shown in Fig. 3", their p. 8 — so both tables start
near x = 0.05, and a caller must treat anything below that as the paper's
silence rather than extrapolate the tables into it.

`wbct_emc_nm_Q5.csv` — the LEFT panel, frame (117.25, 92.06)–(287.23,
202.45) pt.  Legend, read off the four handles: red solid *EMC effect,
NLO*, red dashed *Gluon EMC effect, NLO*, blue solid *EMC effect, NNLO*,
blue dashed *Gluon EMC effect, NNLO* — i.e. F_2A/F_2N and g_A/g_p of the
caption at the two orders.  Columns `R_F2_nlo`, `R_gluon_nlo`,
`R_F2_nnlo`, `R_gluon_nnlo`; covered range x = 0.0467–0.8062.  The black
paths in this panel are the empirical nuclear-matter points of their
Ref. [53] with their error bars, and are not extracted.

`wbct_polemc_nm_Q5.csv` — the RIGHT panel, frame (319.30, 92.45)–(488.83,
202.55) pt.  Two red handles: solid *Polarized EMC effect* = g_1A/g_1p,
dashed *Polarized gluon EMC effect* = Δg_A/Δg_p.  Columns `R_g1`,
`R_deltag`; covered range x = 0.0500–0.8259.  Read back, the trough of
g_1A/g_1p is 0.763 at x = 0.73 against 0.803 for F_2A/F_2N, and that of
Δg_A/Δg_p is 0.795 at x = 0.50 against 0.855 for g_A/g_p — the Letter's
two claims about this figure, which
`fastsim/tests/test_digitized_curves.py::test_wbct_figure_3_range_and_monotonicity`
pins.

```bash
python3 tools/digitize_figure.py --pdf refs/2109.03591.pdf --page 8 \
  --frame 117.25 92.06 287.23 202.45 --xrange 0 1 --yrange 0.6 1.2 \
  --curve "R_F2_nlo:1.0,0.0,0.0:[] 0" \
  --curve "R_gluon_nlo:1.0,0.0,0.0:[ 2.8176797 1.218456 ] 0" \
  --curve "R_F2_nnlo:0.0,0.0,1.0:[] 0" \
  --curve "R_gluon_nnlo:0.0,0.0,1.0:[ 2.8176797 1.218456 ] 0" \
  --grid 300 --out fastsim/polli_fastsim/data/wbct_emc_nm_Q5.csv

python3 tools/digitize_figure.py --pdf refs/2109.03591.pdf --page 8 \
  --frame 319.30 92.45 488.83 202.55 --xrange 0 1 --yrange 0.6 1.2 \
  --curve "R_g1:1.0,0.0,0.0:[] 0" \
  --curve "R_deltag:1.0,0.0,0.0:[ 2.8102425 1.21524 ] 0" \
  --grid 300 --out fastsim/polli_fastsim/data/wbct_polemc_nm_Q5.csv
```

## `av18/fdeut.av18` — the one file here that is not a digitized curve

R. B. Wiringa, Argonne National Laboratory, **the Argonne v18 deuteron wave
function**: the solver's own printout, not a figure and not a
transcription.  Live URL
`https://www.phy.anl.gov/theory/research/deuteron/fdeut.av18`; that host
sits behind Cloudflare bot mitigation which answers `HTTP/2 403
cf-mitigated: challenge` to every non-interactive client, so the copy here
came through the Internet Archive,
`http://web.archive.org/web/20250529230933/https://www.phy.anl.gov/theory/research/deuteron/fdeut.av18`,
which serves the original file bytes (both routes are recorded in
`LiPolGen/data/vmc/README.md`, fetched there 2026-08-29).  This file is
those **raw served bytes**, copied byte for byte from
`LiPolGen/data/vmc/deuteron/fdeut.av18` on 2026-09-15 — 816 110 bytes,
md5 `7f4a361cd833de83b6a1912f8d86bf58`.

The file carries a header, an r-space block (`r u du/dr w dw/dr`, 10 000
rows to r = 100 fm), a momentum-space block, and EM form factors.  The
loader — `polligen.tagged._av18_deuteron_tables`, which finds it with
`importlib.resources` exactly as `_load_curve` finds the CSVs — reads
**only the momentum-space block**: the three columns `k  u(k)  w(k)`,
**201 rows**, k = 0 … 20 **fm⁻¹** in steps of 0.1 fm⁻¹, converted to GeV/c
with ħc = 0.197327 GeV·fm.  Sign convention, which is the whole reason the
file is here: `u(k), w(k) = √(2/π) ×` the **plain** Bessel transforms
∫ j_L(kr) u_L(r) r dr — no i^L, no phase — and **both are positive at low
k**, the S wave having its node at k = 0.4128 GeV/c (2.092 fm⁻¹) and the D
wave staying positive through 1.2 GeV/c.  That is the convention
`TaggedModel._amp2_table` assumes when it applies φ_L = i^L ψ_L, so the
table is loaded with its signs untouched.

Header cross-checks against the block actually read: `dstate = 0.057599`
is the D-state probability `tagged.P_D_AV18_DEUTERON`, and the block
reproduces it, ∫w²k²dk / ∫(u²+w²)k²dk = 0.057600 with
∫(u²+w²)k²dk = 0.999976 (the table's own normalization, to its printed
precision); `qm = 0.269673 fm²` is the positive quadrupole moment that
fixes the sign of w.  The two Cosyn–Weiss landmarks fall where the paper
puts them: w/u = +√2 at k = 0.2988 GeV/c (their "k = 0.30 GeV") and
−1/√2 at 1.0257 GeV/c (their "≈ 1 GeV").

Nothing at run time reads this file; it is a gate.
`evgen/tests/test_tagged.py::test_cosyn_weiss_table_ii_on_av18` builds
`tagged.av18_deuteron_channel()` from it and reproduces Cosyn–Weiss II
TABLE II — the axial node of n_{±1}, A_T∥ = −2 along the spin axis at
k ≈ 0.30 GeV/c and +1 at θ_k = 90° — which the generator's own
Hulthén-type toy radials cannot, their w/u never reaching √2.  Measured
2026-09-15 at the w/u = √2 crossing k = 0.298 GeV/c: −1.937 at the grid
cell nearest θ_k = 0 and −1.998 on a fine (n_c = 4001) near-axis grid,
where n_{+1}/n_0 = 2.7×10⁻⁴ against 0.16 and 0.19 at k = 0.2 and
0.4 GeV/c; +0.999 at the 90° cell; and +0.967 at the k = 1.0 GeV/c row.

## `vmc/li6_ad1.momentum` and `vmc/li6.ad` — the ANL α+d overlap

Two more files that are not digitized curves: R. B. Wiringa *et al.*'s
variational Monte Carlo ⁶Li → α + d cluster overlap, the ab initio input
the generator's two-parameter Hulthén pair stands in for. They are read by
`polligen.tagged.li6_vmc_tables` and selected with
`tagged.li6_alpha_channel(wave='vmc')` / `scripts/money_tagged_azz.py
--cluster-wave vmc`; **nothing published reads them** — every shipped
number is the Hulthén default, bit for bit.

Provenance, quoted from `LiPolGen/data/vmc/README.md`, which fetched them
on 2026-08-29 and is where the URLs live. The live pages
`https://www.phy.anl.gov/theory/research/momenta/li6_ad1.momentum` and
`…/theory/research/overlap_old/li6.ad` sit behind Cloudflare bot
mitigation which answers `HTTP/2 403 cf-mitigated: challenge` to every
non-interactive client, so both copies came through the Internet Archive,

    web.archive.org/web/20250606203105id_/…/momenta/li6_ad1.momentum
    web.archive.org/web/20250617050329/…/overlap_old/li6.ad

the momentum file through the Wayback Machine's raw `id_` form, which
serves the original bytes rather than a rendered page. README.md states the
point explicitly: "Every file here is the **raw text table**, not a
re-typed or figure-digitized copy — the Wayback Machine serves the
original file bytes." The files here are those **raw served bytes**,
copied byte for byte from `LiPolGen/data/vmc/momenta/li6_ad1.momentum`
and `LiPolGen/data/vmc/li6_alpha_d/li6.ad` on 2026-09-16 — 5 167 bytes,
md5 `d4cc62c13e142208e6ca6a36d500fdfb`, and 23 293 bytes, md5
`452ddf07eef2e8751ff69850fc274bbe`.

Physics provenance: AV18+UX (momentum file, 22-Mar-14, VMC 1M samples)
and AV18+UIX (overlap file, Apr-2004, per its own banner) variational
Monte Carlo wave functions. B. S. Pudliner, V. R. Pandharipande,
J. Carlson, S. C. Pieper, R. B. Wiringa, *Quantum Monte Carlo calculations
of nuclei with A ≤ 7*, Phys. Rev. C **56**, 1720 (1997) for the A = 6 wave
function; J. L. Forest, V. R. Pandharipande, S. C. Pieper, R. B. Wiringa,
R. Schiavilla, A. Arriaga, *Femtometer toroidal structures in nuclei*,
Phys. Rev. C **54**, 646 (1996) for the cluster-overlap method.

**Which block the loader reads, and why there are two files.** A momentum
distribution is |ψ_L|² and carries no phase at all, so it cannot supply
the relative S–D sign that the tensor observables read; a 2004 amplitude
table can, but has 200 k samples against 1 M and no printed per-wave
normalization. The loader therefore takes the **magnitude** from
`li6_ad1.momentum`'s **second** block — the S/D split
`K RHOKA0 DRHOKA0 RHOKA2 DRHOKA2`, 51 rows, K = 0.001 … 5 fm⁻¹ — as
ψ̂_L = √ρ_L, and the **sign** from `li6.ad`'s k-space block
`k(fm-1) Aad00(k) Aad22(k)`, 51 rows, k = 0 … 5 fm⁻¹, converted to GeV/c
with ħc = 0.197327 GeV·fm. `li6.ad`'s r-space block, its `ndx / s-wave /
d-wave` line and the momentum file's first (total) block are not read.
Two of the three are cross-checked in `test_tagged.py` all the same: the
total block against the file's own printed norm, and the k-block's own
normalization against what the `ndx` line prints (0.856 / 0.838 / 0.017).
The r-space block is not used at any point.

The sign is taken from the reference's zero **crossings** below 3 fm⁻¹
rather than point by point: past ~3 fm⁻¹ both overlap columns are at the
Monte Carlo noise floor and wander while carrying ~10⁻⁴ of the norm. The
phase is anchored where the reference is largest and stepped across the
crossings from there.

**Sign convention**, which is the whole reason `li6.ad` is here: the
stored ψ̂_L are **plain** Bessel transforms — no i^L, no phase — in the
file's own global phase, ψ̂₂ positive at low k, exactly the convention the
AV18 deuteron's `u(k), w(k)` above are in. The unobservable global phase
is fixed at load by ψ̂₀(k → 0) > 0 (`li6.ad` prints `Aad00 < 0` at low k,
so both columns are negated), and the i^L phase of the momentum-space
amplitude is applied in `TaggedModel._amp2_table`, for every wave alike,
tabulated or analytic — the same treatment `av18_deuteron_channel` gets.

Header cross-checks against the blocks actually read, all measured on
these bytes by `evgen/tests/test_tagged.py` (2026-09-16). The file's own
printed normalizations `4π∫ρK²dK/(2π)³` are 0.81971 (total), 0.80362 (S)
and 0.015861 (D); the trapezoid of the tabulated columns reproduces them
to 7.8×10⁻⁶, 1.6×10⁻⁵ and 1.9×10⁻⁴ relative. P_D is taken from the printed
S/D pair, 0.015861 / (0.80362 + 0.015861) = **0.0193549**, against which
the tabulated columns integrate to 0.0193516. The 2004 overlap file, on
its own normalization S_ad = (2π)⁻³∫k²A²dk, gives 0.85463 and
P_D = 0.0201 — a different Hamiltonian and a fifth of the samples, which
is the size of the spread. The two sign-region boundaries come out of the
signed columns: the α-d **S node at 0.6779 fm⁻¹ = 0.1338 GeV/c** and the
α-d **D node at 2.2498 fm⁻¹ = 0.4439 GeV/c**, each confirmed in the same
0.1 fm⁻¹ bin by a minimum of the momentum file's own ρ_L, so
sign(ψ̂₂/ψ̂₀) = −1, +1, −1 across the three regions — a structure no
node-free analytic form can carry.

Every one of those numbers is also a check against an independent
implementation: LiPolGen's C++ `ClusterWaveSource::VmcAV18` reads the same
two files, and its published acceptance-weighted A_zz^tag at the Yellow
Report high-acceptance optics of 10 × 99.5 GeV/u — −0.5191, −0.2899,
−0.1993, −0.0822, +0.1170 at k = 0.1979, 0.2495, 0.3012, 0.4001,
0.4990 GeV/c — is reproduced here to 3.7×10⁻⁵, its printed rounding, and
its spin-blind accepted-fraction model integral 0.024675932148828
(Hulthén) to every one of fifteen digits.

## Not digitized

Cosyn–Weiss arXiv:2603.23700 FIG. 13 (page 36) is the tagged tensor
asymmetry A_T∥ in light-front variables (α_p, p_pT) the generator does not
carry, so it is a comparison rather than a drop-in.  It did not need
digitizing: page 35 gives the closed form (Eq. 6.12), its extrema
(Eqs. 6.13–6.14) and Table II.  Eq. (6.12) is pinned as an **identity** —
A_T∥ = +1 × A_zz^wf, to machine precision on the whole (k, cosθ_k) grid:
max |A_zz^wf − Eq. (6.12)| = 8.9×10⁻¹⁶ on the deuteron control and
1.1×10⁻¹⁵ on ⁶Li, with the (1 − 3cos²θ_k) factorization holding to
6.0×10⁻¹⁴ — by `evgen/tests/test_tagged.py::test_cosyn_weiss_tensor_gate`,
and TABLE II's own numbers by `…::test_cosyn_weiss_table_ii_on_av18` on the AV18
table above.  Before 2026-09-15 the first of those pinned a mapping
A_T∥ = −2 A_zz^wf, which double-counted the (1 − 3cos²θ_k) factor that
A_zz^wf already carries, and read Eq. (6.14)'s minimum −1/√2 as
Eq. (6.13)'s maximum +√2 on a toy wave function that never reaches √2.
