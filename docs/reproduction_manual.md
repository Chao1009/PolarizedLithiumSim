# Reproducing the simulation results

Every number and figure this program publishes comes out of a script in
this repository.  This manual is the map from result to command: what to
install, what to run, how long it takes, what the answer should be, and
what you cannot reproduce here and why.

It is written from a clean run on the Linux analysis box on 2026-08-26.
Every command below was executed; the runtimes are measured, not
estimated.  Where a number is quoted as an expected output it is the
number that run produced.

**Read this first if you only want one thing:**

| I want to … | Go to |
|---|---|
| check the code is sane in five minutes | §2 |
| redraw the money plots | §4 |
| regenerate the PYTHIA hadronic final state | §5.1 |
| rerun the July fast-sim production | §3.3 |
| rebuild the circulate-able reports | §6 |
| know whether my numbers agree with ours | §7 |
| understand why a number moved | §8 |
| fix an error I hit | §9 |
| know how long something will take | §10 |

---

## 0 · What is reproducible here, and what is not

**Fully reproducible on one Linux box, no accounts, no licences.**
Everything in `fastsim/` and `evgen/`; the PYTHIA 8 hadronic final state;
the BeAGLE e+d control study (the samples stream from a public xrootd
door); the ePIC far-forward acceptance scans; all five reports.  That is
every published number except the four cases below.

**Not reproducible here, and honestly so.**

| what | why | what stands in for it |
|---|---|---|
| BeAGLE breakup of A = 6, 7 | BeAGLE links FLUKA, whose licence is personal and per-user (`tools/beagle/README.md`) | the cluster-IA model of `polli_fastsim/spectator.py`, whose tail is *calibrated* against the official BeAGLE e+d sample (§5.2) — and found wanting, which is the result |
| ePIC calorimeter noise/threshold floor at Σ_h ≈ 0.2–0.5 GeV | an ePIC design number nobody outside the collaboration has (plans/04 #21) | a noise scan, 0 → 25 → 50 → 100 MeV, which brackets it (§4.3) |
| the backward-disk angular resolution | same (plans/04, F3) | the repository's placeholder table, flagged as such everywhere it is used |
| the TENSOR-sector radiative correction | never calculated by anybody (plans/05 §5.5, plans/04 #10); an unpolarized QED study does not bound a tensor one | nothing, and the affected claims say so.  The *unpolarized* collinear-ISR migration IS reproducible here — `--isr`, §4.3 and §8 |

**Reproducible but not bit-for-bit.**  Monte-Carlo results carry a seed
and repeat exactly on the same machine; across machines and NumPy
versions the last digits move.  §8 says how much.

---

## 1 · Environment

### 1.1 What needs what

Install only the rows you need.

| component | needed for | cost |
|---|---|---|
| Python 3.11 + numpy/scipy/matplotlib/pytest | everything | minutes |
| `parton` + four PDF grids | the fast-sim grid backends, the dated money-Δ line, 3 tests | 436 MB, minutes |
| PYTHIA 8 with Python bindings | the hadronic final state (§5.1) | 2.5 min to build |
| `jug_xl-nightly.sif` (eic-shell) + Singularity | BeAGLE samples and ePIC full simulation (§5.2, §5.3) | 3.2 GB container |
| a headless Chromium | the report PDFs only (§6) | 115 MB |

### 1.2 Python

```bash
python3 -V                      # 3.11.4 here
python3 -m pip install numpy scipy matplotlib pytest
```

Reference versions of the run this manual describes: numpy 1.25.1,
scipy 1.11.1, matplotlib 3.7.2, pytest 9.1.1.  Nothing in the repository
needs numpy ≥ 2; `delta_models.py` and `money_delta_20260729.py` carry a
two-sided `np.trapezoid`/`np.trapz` shim so either works.

No package needs installing: both suites and every script run from a
checkout with the repository root on the path, and the scripts put
`fastsim/` and `evgen/` on `sys.path` themselves.  Run them from the
package directory (`cd evgen && python3 scripts/…`).

### 1.3 PDF grids

The pure-Python `parton` package reads LHAPDF grids without LHAPDF:

```bash
python3 -m pip install parton
python3 -m parton update
for s in CT18NLO EPPS21nlo_CT18Anlo_Li6 NNPDFpol11_100 nNNPDF30_nlo_as_0118_A6_Z3; do
    yes | python3 -m parton install $s
done
```

68 + 124 + 83 + 161 MB into `~/.local/share/parton`.  The first three are
what the current code needs; the fourth is used only as the alternative
nuclear PDF in the 2026-07-20 and -07-21 dated scripts, and without it
they exit 2 with an install hint.  Without any of them the three grid
tests skip and the dated money-Δ line cannot run; everything else uses
the toy structure functions and is unaffected.

Check:

```bash
cd fastsim && python3 -c "
from polli_fastsim.structure import PartonF2
print('%.4f' % PartonF2('CT18NLO').f2p(0.1, 10.0))"     # 0.4274
```

### 1.4 PYTHIA 8

**PYTHIA's Python bindings do not ship in the eic-shell container** — it
carries the C++ library and headers only, and `import pythia8` fails
there.  Build them against your own interpreter instead; it takes two and
a half minutes:

```bash
curl -sSL -o pythia8311.tar.gz \
  https://gitlab.com/Pythia8/releases/-/archive/pythia8311/releases-pythia8311.tar.gz
tar xzf pythia8311.tar.gz && cd releases-pythia8311
./configure --prefix=$HOME/Apps/pythia8311 \
            --with-python-config=$(which python3-config) --with-gzip
make -j8 && make install
```

pybind11 is inside the release, so the only prerequisites are a C++
compiler and the Python development headers.  Then, in every shell that
generates events:

```bash
export PYTHONPATH=$HOME/Apps/pythia8311/lib
export PYTHIA8DATA=$HOME/Apps/pythia8311/share/Pythia8/xmldoc
python3 -c "import pythia8; print(pythia8.Pythia('', False).settings.parm('Pythia:versionNumber'))"
# 8.311
```

### 1.5 The eic-shell container

Needed only for §5.2 (BeAGLE samples: `xrdcp`, `pyHepMC3`) and §5.3
(ePIC: `npsim`, `uproot`).  This manual was run against
`~/Projects/eic/local/lib/jug_xl-nightly.sif` (3.2 GB, frozen ≈ September
2024), which carries detector configurations `epic-24.05…24.08` and an
`epic-main` at git 5a7dd057, uproot 5.0.5 and a working
`pyHepMC3.rootIO.ReaderRootTree`.  Fetch a current one with
`eic-shell --upgrade`, and see §5.3 for why the geometry version matters.

```bash
export SIF=~/Projects/eic/local/lib/jug_xl-nightly.sif
singularity exec $SIF bash -lc 'source /opt/detector/epic-main/bin/thisepic.sh; echo $DETECTOR_PATH'
```

### 1.6 A headless browser (reports only)

`reports/build_report.py --pdf` prints through headless Chrome/Edge.  If
the machine has none, playwright downloads one into the user cache with
no root:

```bash
python3 -m pip install playwright && python3 -m playwright install chromium
```

`build_report.py` looks in `~/.cache/ms-playwright/chromium-*/chrome-linux64/`
automatically.

---

## 2 · The five-minute check: the test suites

```bash
cd evgen   && python3 -m pytest tests/ -q     # 305 passed, ~60 s
cd fastsim && python3 -m pytest tests/ -q     # 111 passed, ~13 s
python3 tools/consistency_check.py --verbose  # 25 checks, whole repository
```

416 tests, all of which run without the PDF grids except three of the four
in `fastsim/tests/test_grids.py`, which skip.  These are not smoke tests:
they pin physics identities against independent constructions — the
spin-1 cross section against an explicit density-matrix trace, the
covariant azimuth against a boost-and-rotate construction of the
collinear frame, R1998 against the published fit's own worked values, the
nuclear masses against CODATA, the tensor sign against Cosyn Eq. (27).
If these pass, the machinery is sound and the rest of this manual is
about numbers, not correctness.

The consistency sweep is the other half of those five minutes: 25 checks
in five groups — PHYSICS invariants the simulation must satisfy, SOURCES
against the Yellow Report tables, DRIFT (superseded values a correction
should have removed, statements a rewrite must not drop), ARTEFACTS
(figures, report numbering, and the test counts quoted just above) and
REFERENCES.  Three of the ARTEFACTS checks guard the figures, and since
2026-08-29 the last of them is dependency-aware: the first asks that every
figure a report embeds exists and is registered in
`reports/build_report.py`, the second that it is newer than the script
that draws it, and the third that it is newer than every library module
that script imports — the `polligen` and `polli_fastsim` trees resolved
transitively with `ast` (a package `__init__` counts, because importing a
submodule executes it), plus the digitized CSVs under
`polli_fastsim/data/` that `polarized.py` reads.  A change to a library
module therefore marks every figure downstream of it stale even though no
script was touched, which is exactly the case a hand-run reproduction
misses; the message names the figure, the offending module and both
timestamps, and the fix is to rerun that figure with the command this
manual gives for it (§3, §4).

---

## 3 · Fast simulation (`fastsim/`)

Analytic rates, figures of merit, tagging acceptance, and the dated
money-Δ production line.  Run everything from `fastsim/`; every script
that draws a figure takes `--outdir`; the three that only print
(`validate_inputs`, `diag_sig2_grid`, `_check_reco_mask_invariants`) do
not.

### 3.1 Phase space, rates, and the three early money plots

```bash
cd fastsim
python3 scripts/validate_inputs.py                       # beam/scenario sanity
python3 scripts/phase_space_map.py     --outdir out
python3 scripts/money_delta.py         --outdir out      # money plot 1 (toy Δ)
python3 scripts/money_b1.py            --outdir out      # money plot 2 (b₁/A_zz)
python3 scripts/money_polemc.py        --outdir out      # money plot 3 (polarized EMC)
```

All four, and `coverage_and_stat_maps.py`, take `--run-share` (default
1.0).  The published figures and every number below are the share-1 ones:
the full programme year to this observable in this configuration.  §3.5
prices a share and states what it does and does not move.

Since 2026-08-28 both of the last two draw **digitized theory curves**
(§3.1b), and both keep the pre-digitization shapes behind a flag that also
changes the output stem, so an exploratory run cannot overwrite a
published PNG:

```bash
python3 scripts/money_polemc.py --ion 7Li --pdf grid  --outdir out   # 23 s
python3 scripts/money_polemc.py --ion 7Li --pdf toy   --outdir .     # published
python3 scripts/money_polemc.py --emc-mode constant   --outdir out   # the old 2x / 1x
python3 scripts/money_b1.py     --outdir .                           # published
python3 scripts/money_b1.py     --transfer legacy     --outdir out   # the old 0.87 x 1
```

`money_polemc.py` prints, at 10 fb⁻¹/u on grid inputs, δΔR = 0.0423 /
0.0405 / 0.0600 / 0.1893 at x = 0.09 / 0.28 / 0.45 / 0.71 against a CBT–TMT
separation of 0.044 / 0.040 / 0.034 / 0.059, i.e. 1.04 / 0.98 / 0.57 /
0.31 σ per bin and 3.3 / 3.1 / 1.8 / 1.0 σ at 100 fb⁻¹/u, best in the bin at
x = 0.141 (1.16 σ, 3.65 σ) with 5 of 23 bins above 1 σ.  On toy inputs the
errors are 0.0477 / 0.0511 / 0.0619 / 0.1237 and the best bin is x = 0.089 at
0.92 σ.  Neither unrestricted best bin is the number to quote — both sit
below the window in which the transfer is defined, and §3.1b explains why —
so the same run prints the window-restricted line as well: x = 0.355 at
0.84 σ and 2.66 σ on grid inputs, 0.72 σ and 2.27 σ on toy.  Both figure
panels shade that window, and the upper one draws TMT's published
nuclear-matter curve untransferred as a dotted line, where it lies on top
of CBT's below x ≈ 0.3.  The separations were 0.007 / 0.016 / 0.049 /
0.118 with the old constants: the digitization does not move the
discriminating region, it halves the separation there and takes its x
dependence from the calculations instead of from the hand-written EMC
table.

`money_b1.py` prints the transfer it used, the two signal curves and the
per-bin significance.  At the rank-2 default, 0.921947 × 2/6 = 0.307316:
|A_zz| = 2.42×10⁻⁴ / 4.83×10⁻⁴ / 7.79×10⁻⁴ / 1.42×10⁻³ / 1.37×10⁻³ /
3.34×10⁻³ at x = 0.005 / 0.01 / 0.03 / 0.07 / 0.2 / 0.5 for the Miller
curve and 1.15×10⁻⁵ / 1.15×10⁻⁵ / 6.89×10⁻⁶ / 1.32×10⁻⁶ / 6.94×10⁻⁵ /
4.13×10⁻⁴ for the CDKS convolution, against δA_zz = 9.31×10⁻⁵ / 8.84×10⁻⁵ /
9.73×10⁻⁵ / 3.55×10⁻⁴ / 1.17×10⁻³ at x = 0.0035 / 0.0089 / 0.0282 / 0.2818 /
0.5623 for P_zz = 0.60 (6.98×10⁻⁵ / 6.63×10⁻⁵ / 7.30×10⁻⁵ / 2.66×10⁻⁴ /
8.77×10⁻⁴ at 0.80).  That is 1.7 / 4.8 / 7.6 / 10.9 / 5.6 / 5.8 σ per bin
for Miller at P_zz = 0.6 and ≤ 0.2 σ everywhere for CDKS.  `--transfer
legacy` gives 5.0 / 13.7 / 21.6 / 31.0 / 15.9 / 16.5 σ.  The error rows are
untouched by any of this; only the signal moved.

### 3.1b The digitized theory curves

The polarized-EMC and b₁ curves the two money plots draw are the published
figures, read back from the PDFs' own path operators — not pixel picks.
`fastsim/polli_fastsim/data/SOURCES.md` carries paper, page, figure, frame
box, axis ranges, the legend handle each curve was identified by, and the
exact command; the tool is `tools/digitize_figure.py`, which needs PyMuPDF
as a **dev-time-only** dependency (`pip install --target <dir> pymupdf`,
then put `<dir>` on `PYTHONPATH`; nothing at runtime opens a PDF).  The
PDFs come from `python3 refs/find_ref.py --fetch`.

```bash
export PYTHONPATH=<dir-with-pymupdf>
python3 tools/digitize_figure.py --pdf refs/nucl-th_0605061.pdf --page 7 --inspect
python3 tools/digitize_figure.py --pdf refs/nucl-th_0605061.pdf --page 7 \
  --frame 82.43 55.92 296.11 187.99 --inspect      # groups + legend handles
python3 tools/digitize_figure.py --pdf refs/nucl-th_0605061.pdf --page 7 \
  --frame 82.43 55.92 296.11 187.99 --xrange 0 1 --yrange 0.55 1.2 \
  --curve "R_unpol:0.182,0.19,0.573:[ 1.50769 ] 0" \
  --curve "R_pol_eq26:0.93,0.111,0.141:[] 0" \
  --curve "R_pol_eq23:0.93,0.111,0.141:[ 0 1.50769 ] 0" \
  --grid 300 --out fastsim/polli_fastsim/data/cbt_polemc_7Li_Q5.csv
```

The other five commands are in `SOURCES.md` verbatim.  Two checks come free
and both are in `fastsim/tests/test_digitized_curves.py`: CDKS's Fig. 5 at
Q² = 2.5 GeV² reproduces the same curve extracted from Fig. 4 on a different
page with a different calibration to better than 2×10⁻⁶, and their solid
curves equal the sum of their own dashed and dotted ones to 2×10⁻⁷.  The
Close–Kumano integral ∫b₁ dx over the digitized range is 4.59×10⁻⁴ for the
CDKS convolution and 5.92×10⁻³ for Miller's total — the model that violates
the sum rule violates it, and the one that respects it does, which is the
sharpest statement the tables support given that they cover only what the
papers plot.

The two polarized-EMC camps are compared through the ratio of effects
(1 − R_pol)/(1 − R_unpol), which the report and `plans/01` quote:

```bash
cd fastsim && python3 - <<'PY'
import numpy as np
from polli_fastsim import polarized as P
print("  CBT", " / ".join("%.2f" % P.cbt_ratio_of_effects(x)
                          for x in (0.40, 0.45, 0.50, 0.60)))
print("  TMT", " / ".join("%.2f" % P.tmt_ratio_of_effects(x)
                          for x in (0.40, 0.45, 0.50, 0.60)))
xs = np.linspace(0.30, 0.83, 10001)
r = np.array([P.cbt_ratio_of_effects(x) for x in xs])
print("  CBT minimum %.2f at x = %.3f" % (r.min(), xs[r.argmin()]))
xs = np.linspace(*P.curve_x_range("cbt_polemc_7Li_Q5"), num=20001)
d = np.array([P.cbt_unpolarized_emc_ratio(x) for x in xs]) - 1.0
i = np.where(np.diff(np.sign(d)) != 0)[0]
print("  CBT unpolarized = 1 at x =", " / ".join("%.3f" % xs[j] for j in i))
for x in (0.36, 0.40, 0.45, 0.65):
    c, t = P.cbt_polarized_emc_ratio(x), P.tmt_polarized_emc_ratio(x)
    u = P.cbt_unpolarized_emc_ratio(x)
    print("  x=%.2f  depletion CBT %.3f  TMT->7Li %.3f  7Li unpol %.3f"
          "  sep %.4f" % (x, 1 - c, 1 - t, 1 - u, abs(c - t)))
PY
```

It prints 2.25 / 1.69 / 1.41 / 1.14 for CBT against 1.01 / 0.98 / 1.00 /
1.08 for TMT at x = 0.40 / 0.45 / 0.50 / 0.60, a CBT minimum of 1.06 at
x = 0.696, and the two zeros of the denominator at x = 0.280 and 0.840
(TMT's is at 0.721, from `tmt_polemc_nm_Q10.csv` the same way).  Those
zeros are why the two camps are put on a common baseline by a valence-window
strength factor rather than by the pointwise ratio, and CBT's Eq.-23 curve
never meets its unpolarized one in between; the Eq.-26 curve does, at
x = 0.651.  The last block prints the depletions the report reads the
comparison off: at x = 0.40 / 0.45 / 0.65 the transferred nuclear-matter
depletion is 0.039 / 0.048 / 0.083 against ⁷Li's own unpolarized
0.034 / 0.048 / 0.087 — within 0.005, which is what "polarized ≈
unpolarized" means — and against CBT's 0.077 / 0.082 / 0.094, so ΔR
separates by 0.040 at x = 0.36, 0.038 at 0.40, 0.034 at 0.45 and 0.011 at
0.65.

That strength factor is one constant, fitted over 0.35 < x < 0.65 and
applied at every x, so it is worth knowing what it does outside its window
before reading a separation off the money plot.  `money_polemc.py` prints
the answer beside every number it quotes — the transferred separation and,
in brackets, the separation between the two *published* curves at the same
x, from `polarized.tmt_published_emc_ratio`:

```
  x=0.09: dDR(10/fb)=0.0423  CBT-TMT sep=0.044  (1.04 sigma here, 3.28 at 100 fb^-1/u)  [published curves differ by 0.002 here]
  x=0.45: dDR(10/fb)=0.0600  CBT-TMT sep=0.034  (0.57 sigma here, 1.80 at 100 fb^-1/u)  [published curves differ by 0.038 here]
  valence window 0.35<x<0.65: best bin x=0.355, 0.84 sigma at 10 fb^-1/u, 2.66 at 100 fb^-1/u; 0 of 3 bins above 1 sigma
  published-curve separation: max 0.0078 over 0.028 < x < 0.3, where the transferred pair shows 0.0398-0.0456; max 0.1048 above x = 0.35
```

The two calculations agree to better than 0.008 over the whole decade in
which both tables have data below x = 0.3, so the ≈ 0.04 the projection
plots there is the ⁷Li ← nuclear-matter rescaling and not a disagreement
between the papers; the honest reach is the window-restricted line, 0.84 σ
per bin at 10 fb⁻¹/u and 2.66 σ at 100 (0.72 and 2.27 σ with `--pdf toy`,
which is what the published PNG draws).  Report 0 §5.3, `fastsim/README.md`
and `plans/01` all say so, and `tests/test_digitized_curves.py::
test_the_two_published_curves_agree_at_low_x_and_part_in_the_valence_window`
pins both ends of it.

### 3.2 Spectator tagging acceptance

```bash
python3 scripts/tagging_acceptance.py
```

The number the proposal asks for.  Since 2026-08-28 the script evaluates
every channel at each beam configuration and at four optics — the Yellow
Report high-acceptance and high-divergence rows of that configuration
(plans/10; the envelope is a rectangle 10(σ_h, σ_v) in angle applied to
each fragment's azimuth), the lithium tagging optics of Report 1 §6.1,
and the legacy proton-derived 73 μrad for reproduction.  Expected at the
β = 0.30 central short-range scale, 5 × 41 / 10 × 100 / 18 × 275:

```
6Li α-tag   YR high-acceptance  0.0186 / 0.0168 / 0.0261  tagging optics  0.3541 / 0.2754 / 0.2913  legacy 73 μrad  0.1679 / 0.0262 / 0.0276
6Li d-tag   YR high-acceptance  0.0930 / 0.0814 / 0.1389  tagging optics  0.5693 / 0.5108 / 0.5665
7Li α-tag   YR high-acceptance  0.9690 / 0.9683 / 0.9748  tagging optics  0.9877 / 0.9919 / 0.9941
7Li t-tag   YR high-acceptance  0.7777 / 0.9185 / 0.9386  tagging optics  0.7797 / 0.9218 / 0.9406
```

Every row above rose on 2026-08-28, and the ⁷Li triton row by two orders
of magnitude, because `farforward.route_charged` gained an **over-rigid
branch**: a fragment at R > 1.05 bends less than the beam, and the pot
dispersion carries it onto the Roman-Pot silicon on the inner side of the
bend wherever the displacement clears that configuration's blind block,
48 / 32 / 16 mm.  The destination is measured, not assumed — an R = 1.286
triton lands at dx = +66 mm in 60 of 60 events at every configuration
with a ZDC deposit behind it (plans/09 B1, `tools/fullsim/README.md`) —
so the ⁷Li α + t double tag exists and the t-tag row is 78–94% rather
than the 0.033 / 0.004 / 0.005 the routing-as-lost picture gave.  Almost
all of that row is route 6 (0.745 / 0.915 / 0.933), which is why the
tagging optics moves it by under 0.3 points: this tag is dispersive, not
angular.  The block is per configuration, so a caller that leaves
`pot_config` at its 18 × 275 default over-accepts at the two lower
configurations — 0.0277 / 0.0249 instead of 0.0186 / 0.0168 on the ⁶Li α
row.

The ⁶Li α is beam-blind (rigidity ratio 0.99813, inside the ±5% near-beam
band): at the Yellow Report optics its near-beam tail is inside the
envelope at every configuration and what survives is the 1.5% slice that
falls below R = 0.95 into the Roman-Pot window, plus 0.1 to 1.0 points of
over-rigid inner branch; the tagging optics (1/7–1/13 of the luminosity)
recovers the tail.  §5.2 still matters — the tail is the
least trustworthy part of the cluster model.

**Two samplers, one channel — read the ⁶Li α numbers with this.**  The
table above is `spectator.spectator_lab_kinematics`, which carries one
partial wave per channel (`ClusterChannel.l_wave` = 0 for α–d).  The
tagged generator behind money plot 4 (`polligen.tagged`) carries the full
S + D expansion, and the D wave is the whole difference: its S-wave radial
*is* `spectator.momentum_density`, ⟨k⟩ = 0.1071 GeV/c in both, while the
D wave is hard at ⟨k⟩ = 0.2778 and at P_D = 0.0867 pulls the channel mean
to 0.1219 and the off-rigidity R < 0.95 slice from 1.5% to 2.5%.  Set
`p_d = 0` and the two agree quantile by quantile
(`test_boost_matches_fastsim_spectator`).  So the ⁶Li α tag at the Yellow
Report optics is a 1.5% Roman-Pot slice on the pure spectator model and
2.5–2.9% on the tagged one, and neither is wrong — the D wave is not an optional tail, it
*is* the tensor observable (with P_D = 0 the α–d density is m-independent
and A_zz^tag vanishes identically), so the tagged observables are quoted
on the S + D spectrum and this table on the S-wave one.  For ⁷Li the two
densities are identical (both pure P wave), and the 0.7 points by which
the table's 0.9690 at 5 × 40.8 exceeds `tagged_polarimetry_7li.py`'s
0.9620 is **two different acceptance definitions, not two different
densities**.  This table is `1 - lost`, i.e. any far-forward system; the
tagged script's mask is the Roman Pots alone (route 1 | 4), and at
5 × 40.8 the B0 carries 1.1% of the ⁷Li α and the over-rigid branch a
further 0.2%.  Like for like, the pure
model's Roman-Pot tag is 0.9568 / 0.9668 / 0.9721, and restricted to
k ≤ 1.2 GeV/c where the tagged model's momentum grid ends it is
0.9626 / 0.9683 / 0.9736 against the tagged sampler's
0.9620 / 0.9678 / 0.9728 — 0.1 point at every configuration.  (The grid
truncation alone is worth +0.6 point at 5 × 40.8 and +0.2 at the other
two on the Roman-Pot mask, +0.1 uniformly on `1 - lost`; the B0 fraction
is 0.0000 at 10 × 99.5 and 18 × 117.9, where the two definitions
coincide.)  `tagged_polarimetry_7li.py` now prints both definitions.

### 3.3 The July money-Δ production, and the R that changes it

This is the line that produces the discovery-luminosity numbers.  It has
two switches that matter, and a third — `--run-share` on
`money_delta_realistic.py`, default 1.0 — that leaves every `L_5σ` below
exactly where it is (§3.5) and is here only so a run plan can quote the
years.

```bash
# reproduce the dated notes exactly (the default)
python3 scripts/money_delta_20260729.py --outdir out         # ~16 s, 16 PNGs
python3 scripts/money_delta_realistic.py --outdir out        # ~36 s

# the same with a corrected R = σ_L/σ_T
python3 scripts/money_delta_20260729.py  --r-model r1998 --outdir out
python3 scripts/money_delta_realistic.py --r-model r1998 --configs low,mid,top --outdir out
                                                             # ~54 s (3 configs)
```

`--r-model simplified` (the default) is the frozen July form, whose Θ
factor is multiplied into all three terms of the SLAC fit and therefore
returns R = 1 over 38% of the region that carries the sensitivity — the
defect of code review item S1.  `--r-model r1998` is the published
three-form E143 fit from `polli_fastsim.structure`.  **Both are correct
to run**: the default reproduces the dated notes, the switch produces
today's answer.  A non-default model writes to `out/money_delta_r_<model>`
so it cannot overwrite the July figures, and every figure caption carries
which R produced it.

Expected (§7 has the table), at the full programme year (`--run-share 1`,
the default; the reach is invariant under the share in any case): `L_5σ`
at Δ/F₁ = 10⁻³ goes from 135.3 / 131.3 / 274.6 fb⁻¹/u (frozen) to
**67.5 / 65.8 / 155.1** (published R) for LOW / MID / TOP.  **Caveat (2026-08-27):** this
script carries its own `ALL_CONFIGS` at 27.5 / 50 / 137.5 GeV/u, the
pre-correction rigidity-scaled energies, so only its TOP row is at a
machine configuration; the low and mid rows reproduce the dated notes
and are not current results (plans/10).  The current-energy toy reach
is `scripts/money_delta.py`: L_5σ(Δ/F₁ = 10⁻³, P_zz = 0.8) = 16.7 /
16.3 / 21.8 fb⁻¹/u at 5 × 40.8 / 10 × 99.5 / 18 × 137.5.

The six dated notes in `fastsim/notes/` carry a banner recording this;
their numbers are still exactly reproducible with the default.

Structural invariants of the reco path are checked separately:

```bash
python3 scripts/_check_reco_mask_invariants.py     # S1–S6, exit 0
```

### 3.4 The whole dated line

`fastsim/`scripts/money_delta_<date>.py`` are the dated production
scripts, each frozen against the note of the same date in
`fastsim/notes/`.  **All seven run** (verified 2026-08-26); together they
are about two minutes:

```bash
for d in 20260715 20260720 20260721 20260724 20260725 20260728 20260729; do
    python3 scripts/money_delta_$d.py --outdir out/$d
done
```

Two of them need the fourth PDF grid of §1.3 and exit 2 with an install
hint without it (-07-20, -07-21).  `money_delta_pdfgrid.py` is the
slowest thing in the repository at ~5 minutes, because it evaluates the
grids point by point.

A dated script is a *reproduction*, not a current result: where its
answer has since been superseded — as the R = σ_L/σ_T defect superseded
every `L_5σ` in the July line — the note carries a banner saying by how
much, and the script still returns its original number.

### 3.5 The run-plan share (plans/07 WP2)

Every published number in this manual is quoted at the **full 10 fb⁻¹/u
programme year in its own configuration** — its own spin states, its own
far-forward optics, its own isotope.  They are therefore alternatives and
not a programme: nothing here says how a real year would be divided
between them.  That division is parameterized, not chosen, by one flag:
`--run-share` in `fastsim/` and `--lumi-fraction` in `evgen/`, default
1.0 in both, carried into the library as `fom.Scenario.run_share` and
applied at its single luminosity line.

It is deliberately a third factor, distinct from the two luminosity
fractions that were already there: `farforward.Optics.lumi_fraction` is
what a de-squeezed β*_x costs at fixed wall time (1/7–1/13 for the
lithium tagging optics), and `bookkeeping.SpinCategory.lumi_fraction` is
how one measurement's own luminosity divides between its fill states
(0.5/0.5 for the tensor flip plan).  The three multiply; only the first
is a programme decision.  Every script prints all of them it uses, and a
non-default share appends its key (`share0p25`) to the output stem so it
cannot overwrite a published figure.

The law it obeys, verified end to end on 2026-08-28:

```bash
cd fastsim
python3 scripts/money_delta.py --ion 6Li --pdf toy --run-share 1    --outdir out
python3 scripts/money_delta.py --ion 6Li --pdf toy --run-share 0.25 --outdir out
```

`L_5σ(Δ/F₁ = 10⁻³, P_zz = 0.8)` = **16.7 / 16.3 / 21.8 fb⁻¹/u at both
shares** — it is the luminosity the measurement must accumulate, and the
share buys wall-clock time rather than physics — while the error on the
Δ/F₁ scale after one programme year doubles, 2.586 / 2.556 / 2.954 ×10⁻⁴
→ 5.172 / 5.112 / 5.907 ×10⁻⁴, and the programme years to 5σ go from
1.67 / 1.63 / 2.18 to 6.69 / 6.53 / 8.72.  The invariance is exact, not
approximate: `money_delta_realistic.py --configs mid,top` returns
131.26 / 274.64 fb⁻¹/u at `--run-share 1` and at `--run-share 0.25`
alike.

The same law on the reconstructed coherent chain:

```bash
cd evgen
python3 scripts/money_cos2phi_coherent_reco.py --config 1 --optics tagging \
        --exact --n-mc 600000 --lumi-fraction 1    --outdir .
python3 scripts/money_cos2phi_coherent_reco.py --config 1 --optics tagging \
        --exact --n-mc 600000 --lumi-fraction 0.25 --outdir .
```

The combined one-year δa_e over the seven |t| bins goes **0.00104 →
0.00208**, and per bin 0.00202 / 0.00230 / 0.00268 / 0.00221 / 0.00322 /
0.00537 / 0.00958 → 0.00405 / 0.00460 / 0.00535 / 0.00441 / 0.00644 /
0.01074 / 0.01916 — a factor two, bin by bin, with the recovered central
values unmoved.  The second run writes
`money_cos2phi_coherent_reco_6Li_c1_tagging_share0p25.png`, not the
published stem.

The coherent line of Plan A in plans/07 WP2 is the same arithmetic on the
full-year output of `tagging_optics.py` (§4.5) rather than a rerun, the
share entering the yield linearly, the 5σ floor as 1/√f and the years to a
target amplitude as 1/f.

On the inclusive reconstructed side the same check found a real defect
and fixed it.  `InclusiveSampler` recovers a per-cell cross section by
dividing `proj.n_events` by the scenario luminosity, and it was dividing
by the programme figure while `n_events` already carried the share — so
every cross section came out scaled by the share and the caller's
`lumi_pb` applied it a second time.  Dividing by the effective luminosity
(`polligen/sample.py`, 2026-08-28) is the fix; at share 1 nothing moves,
so no published number changes.  What it looked like:

```bash
cd evgen
python3 scripts/money_cos2phi_reco.py --n-mc-per-cell 60 --outdir . \
        --lumi-fraction 0.25
```

reported N(1 yr) = 1.17×10⁷ at sweet spot 1 where 4.67×10⁷ is right (the
share-1 value is 1.87×10⁸), and δÂ = 4.66×10⁻⁴ where 2.33×10⁻⁴ is right
against 1.16×10⁻⁴ at share 1.  It now writes
`money_cos2phi_reco_6Li_share0p25.png` and
`money_delta_extracted_reco_6Li_share0p25.png`, leaving 5R and 7R alone.

`fastsim/tests/test_run_share.py` and `evgen/tests/test_run_share.py` pin
both halves of this (errors × 1/√share, `L_5σ` invariant, the stem
guards, and the three fractions kept apart).  What must **not** be
rescaled this way is anything that is not pure statistics: the full bars
of Report 2 money plot 7R, which carry the shape fit, the response Monte
Carlo and the unfolding prior spread, do not shrink with luminosity and
their statistical part has to be separated first; the Table 4 bars of
that report are statistical only and do obey the law.

---

## 4 · Event generator (`evgen/`)

`polligen` is the doubly polarized e + ⁶,⁷Li generator.  Run from
`evgen/`; every script takes `--outdir`, and every script that draws a
Monte-Carlo sample takes `--seed`.

### 4.1 Phase space, closure, tagged observables

```bash
cd evgen
python3 scripts/phase_space_bins.py        --outdir .
python3 scripts/closure_fom.py             --outdir .    # estimator closure vs the analytic FOM
python3 scripts/money_tagged_azz.py        --outdir . --events 400000   # 3 s, money plot 4
python3 scripts/tagged_polarimetry_7li.py  --outdir .                  # 3 s
python3 scripts/target_mass_bound.py                                   # 2 s, prints only
```

`closure_fom.py` is the one to run if you doubt the generator: it
compares the spread of pseudo-experiments against the analytic error
formulas over ~65 x bins.

`target_mass_bound.py` answers, in four printed blocks, how much the
massless A_∥ = D·g₁/F₁ of the master formula costs: the W² ≥ 10 GeV² cap
γ² ≤ M²/(W²_min − M²) = 0.0965 and the measured maximum per
configuration; the shift at the four money-plot-5 sweet spots of each
⁶Li setting, beside the other O(γ²) effect at the same spots — the error
the lab-angle shortcut would make on the cos 2φ azimuth, which is
exactly zero for a massless target and reaches 1.38 and 1.21 mrad at the
mid and top configurations and 5.94 mrad at the low one (Report 2 §4.1);
the inverse-variance-weighted bias on the polarized-EMC ΔR — high,
because the exact A_∥ is the larger one and the massless inversion
passes that on — beside the δΔR those same weights give; and the
σ-weighted shift of the tagged-triton A_∥ overlay.  `--pdf grid` repeats
all four on CT18/NNPDFpol in about seven minutes and reproduces the
published spot's shift to two digits; there the two lowest-Q² sweet
spots of each configuration (Q² = 1.135 GeV² at the two lower ones,
1.462 at the top) sit below the grids' Q² floor and are reported, with
their Q², as not evaluable.  It is the source of every A_∥ γ² number in
Reports 0–2 — the b₁–b₄ cos 2φ leakage bound of Report 1 §2 and
Report 2 §4.1 is a separate calculation (plans/08 D2) — and of the γ²
statements in the docstrings of `asymmetries.a_parallel` and
`polligen/xsec.py` and in the comment above the g₁/F₁ division in
`fom.project_observables`.

The two tagged scripts take `--config {0,1,2}` and `--optics
{menu,legacy,high-acceptance,high-divergence,tagging}` since 2026-08-28
(plans/09 B2, B3); `menu`, the default, is the Yellow Report
high-acceptance optics of that configuration plus the lithium tagging
optics with its luminosity fraction, and `legacy` reproduces the retired
proton-derived 73/164 μrad pair.  The published figures are the DEFAULT
combination, `--config 1 --optics menu`; any other combination writes
`money_tagged_azz_6Li_<key>_<optics>.png` /
`tagged_polarimetry_7Li_<key>_<optics>.png`, so no exploratory run —
`--optics legacy` at the default configuration included — can overwrite
them.

**Money plot 4 at 10 × 99.5** (`--events 400000`, seed 20260713, and the
whole programme year to this observable — `--lumi-fraction 1`, the
default).  The α
tag is 0.0247 at the Yellow Report optics against 0.3061 at the tagging
optics — but the tagging optics runs at L/L_HA = 1/13.3, so in tagged
events per year the tagging optics costs 7% here (acc × L = 0.0247 vs
0.0231), which is the worst of the three configurations.  What the 7%
buys is the REACH, and that is the result: the median accepted
spectator momentum is 0.322 GeV/c at the Yellow Report optics with
**nothing at all below k = 0.15 GeV/c**, against 0.162 GeV/c with 44%
below 0.15 at the tagging optics (0.348/0.323/0.334 and
0.146/0.162/0.161 GeV/c at 5 × 41 / 10 × 100 / 18 × 275; the unfolded
sample has median 0.093 and 74% below 0.15).  The rate trade is even or
better at the other two: acc × L/L_HA is 0.0284 → 0.0541 at 5 × 41,
0.0247 → 0.0231 at 10 × 100 and 0.0265 → 0.0331 at 18 × 275, so at two
of the three configurations the tagging optics gains the reach and the
rate together, and pays for it only at 10 × 100.  At every published optics
the ⁶Li α tag admits only the high-k tail; the tagging optics is what
turns money plot 4 from a one-point measurement into a curve.

The right panel's coloured curves are the truth weighted by each optics'
own θ_k acceptance (`tagged.acceptance_weights`,
`tagged.azz_tensor_curve`), and that is what the markers measure: at
k ≈ 0.325 GeV/c the folded A_zz is +0.491 against a weighted truth of
+0.455 at the Yellow Report optics and −0.070 against −0.098 at the
tagging optics, while the θ_k = 90° curve — drawn as a grey reference —
says −0.482 at both.  The published version before this round overlaid
that 90° curve alone, on samples whose ⟨|cos θ_k|⟩ is 0.79 (the
off-rigidity window slice, longitudinal) or 0.40 (the near-beam tail,
transverse); its ±0.5 swing between two optics at k ≈ 0.3 GeV/c was the
envelope, not the wave function.  Closure: the residual is |ΔA_zz| ≤
0.018 in every populated bin at both optics, which at the figure's own
4 × 10⁵ events is within 1.4σ (Yellow Report) and 2.2σ (tagging).  At
`--events 8e6` the errors fall to 0.003–0.03 and one bin, k = 0.175 at
the tagging optics, reads 3.3σ on the same residual: that is the bin in
which the acceptance turns on — nothing at all is accepted below
k = 0.189 GeV/c at the Yellow Report optics — so the truth at the bin
*centre* is not what the bin-averaged marker measures.  Averaged over the
bin as the marker is, every populated bin is within 1.4σ at 8 × 10⁶ too.

**⁷Li polarimetry at 10 × 99.5** (`--events 300000`).  Nothing moves, and
that is the deliverable: the Roman-Pot tag is 0.9678 at the Yellow Report
optics against 0.9917 at the tagging optics, the folded ⟨P₂⟩ slope
−0.1947 and −0.1964 against the analytic −0.2000 (the retired 73 μrad
gave −0.1929), and the median δA_∥ 0.01152 against 0.01140 *at equal
generated statistics*.  At equal luminosity the
tagging optics multiplies every ⁷Li error bar by 2.83 / 3.87 / 3.15 at
5 × 41 / 10 × 100 / 18 × 275: ×1.02 in acceptance for ×1/8–1/15 in
luminosity is a factor 8–15 net loss, the exact inverse of ⁶Li.  ⁶Li and
⁷Li want different machine optics and are different runs.

Several of those numbers moved in the fourth decimal on 2026-08-28, from
two independent causes, and both are worth recording because neither is an
optics effect.

The first is plans/08 D7, which made every `Ion` slot per-nucleon, `TRITON`
included, so the struck triton's two neutrons both contribute to g₁ where
only one used to; the accepted sample is cross-section-weighted and the
cross section carries that g₁, so the tag follows.  Roman-Pot tag
0.9614 / 0.9676 / 0.9726 → **0.9620 / 0.9678 / 0.9728** (Yellow Report) and
0.9807 / 0.9916 / 0.9920 → **0.9807 / 0.9917 / 0.9919** (tagging); median
δA_∥ 0.01150 / 0.01138 → **0.01152 / 0.01140**; the error-bar penalty is
**2.83 / 3.87 / 3.15** either way, the fourth-decimal shift leaving
√(acc × L ratio) where it was.  The ⟨P₂⟩ slopes are unchanged.  What moves visibly is the analytic D(y) g₁ₜ/F₁ₜ overlay the
right panel draws, by +8.4% at x = 0.005, +5.3% at 0.01, +2.4% at 0.03,
+0.8% at 0.1, −0.5% at 0.5 and −0.7% at 0.7.

The second is the over-rigid routing branch of the same day
(`farforward.over_rigid_route`, plans/09 B1), which stopped treating an
R > 1.05 fragment as lost.  It does not touch the Roman-Pot tag, which
masks on the pots directly, but it moves `acc(any far-fwd)` — the
`1 − lost` definition — wherever the α's Fermi smearing carries it over
rigidity: 0.9683 / 0.9676 / 0.9726 → **0.9717 / 0.9703 / 0.9752** (Yellow
Report) and 0.9876 / 0.9916 / 0.9920 → **0.9904 / 0.9941 / 0.9942**
(tagging).  A review later the same day found that the script was passing
no `pot_config` at either of its two `route_charged` call sites, so it was
testing the over-rigid displacement against the *narrowest* blind block,
16 mm, at all three configurations; with the per-configuration 48 / 32 /
16 mm it publishes now, `acc(any far-fwd)` settles at **0.9699 / 0.9690 /
0.9751** and **0.9886 / 0.9929 / 0.9942**.  The Roman-Pot mask is
`pot_config`-blind by construction — routes 1 and 4 end at R = 1.05 where
route 6 begins — so no `acc(RP)` digit moves with it.  The isolation probe
that separates the causes forces
`farforward.over_rigid_route` to `False` and returns `acc(any far-fwd)` to
`acc(RP)` exactly, 0.9678 and 0.9917 at 10 × 99.5:

```python
# saved as probe.py and run from evgen/ as: python3 probe.py
import pathlib, runpy, sys
sys.path.insert(0, str(pathlib.Path("../fastsim").resolve()))
sys.argv = ["tagged_polarimetry_7li.py", "--config", "1", "--outdir", "/tmp"]
import polli_fastsim.farforward as ff
ff.over_rigid_route = lambda *a, **k: False
runpy.run_path("scripts/tagged_polarimetry_7li.py", run_name="__main__")
```  Report 3's
Table 6, its caption and its change log carry the post-2026-08-28 tags,
and `plans/09` §B3 the same triples.

The **Cosyn–Weiss deuteron limit** (`plans/05` §5.4) is a gate rather than a
figure, and it runs as a test:

```bash
cd evgen && python3 -m pytest tests/test_tagged.py::test_cosyn_weiss_tensor_gate -q
```

It pins the tagged sampler's wave-function tensor asymmetry against
Cosyn–Weiss II (`refs/2603.23700.pdf` p. 35, Eqs. (6.12)–(6.14) and
TABLE II), and needed no digitization because the paper gives the result in
closed form.  At fixed spectator momentum k = 0.3012 GeV/c the ratio
A_zz^wf / P₂(cos θ_k) is 0.99940, with a spread below 10⁻⁵, in every
angular cell away from the zero of P₂ at cos θ_k = 1/√3, where the ratio is
unbounded and the test excludes the cells with |P₂| ≤ 10⁻³; so their
angular factorization holds; the radial envelope peaks at
k = 0.3098 GeV/c against the 0.30 GeV/c they quote for AV18; and under
A_T∥ = −2 A_zz^wf the extremes are +0.99967 at the cell nearest θ_k = 90°
and −1.93782 at the outermost cell (|cos θ_k| = 0.98958), which extrapolate
through that factorization to +1.000 and −2.000 against their tabulated +1
and −2, the whole curve staying inside their [−2, 1].

### 4.2 The cos 2φ money plots, truth level

```bash
python3 scripts/money_cos2phi.py            --outdir .   # money plot 5
python3 scripts/money_cos2phi_coherent.py   --outdir .   # money plot 6
python3 scripts/money_delta_extraction.py   --outdir .   # money plot 7
```

Panel (b) of money plot 6 shades the reconstructed analysis window
0.017–0.25 GeV² of §4.3 and marks |t| = 0.05 GeV², where the digitized
deuteron anchor begins and below which the band is the linear model
extrapolated.  Both are annotations: no printed number of this script
moves with the window (A_hat = 0.0589 ± 0.0018 against a truth of 0.0600,
⟨a₂⟩_tag = 0.0360, 5σ floor 0.0090).

### 4.3 The reconstructed level — 5R, 7R, 6R

```bash
python3 scripts/reco_chain_figures.py       --outdir .   # chain schematics + acceptance curves
python3 scripts/money_cos2phi_reco.py       --outdir .   # 5R and 7R
python3 scripts/money_cos2phi_coherent_reco.py --config 0 --optics tagging --n-mc 6000000 --outdir .  # 6R
```

The 6R command is the published one (2026-08-28): the lithium tagging
optics of Report 1 §6.1 at the low configuration, with the pots following
the 10σ envelope in both planes (`--optics tagging` sets the divergence,
the cutout and the luminosity from `reco.tagging_optics_point`).  Without
it the script runs the pre-2026-08-27 legacy geometry, which is kept for
reproduction only.  The published stem `money_cos2phi_coherent_reco_6Li.png`
belongs to that one combination — `--config 0 --optics tagging` with the
default ratio fit, the assumed (u₁, u₂) and the published |t| edges; any
other run appends its keys, so none of the exploratory commands below can
overwrite money plot 6R.  `nearbeam_reach_gain.py` carries the same guard
on `--fit` and on `--t-edges` (all since 2026-08-28), and both — with
`money_cos2phi_reco.py` and `money_tagged_azz.py` — guard
`--lumi-fraction`, the programme share of §3.5, which is 1.0 in every
published run.  `--ensemble 20` repeats the one-year pseudo-experiment
and prints the bias test of Table 5; `--exact` switches the Poisson draw
off for the systematic scans; `--config 1/2` gives the other two
configurations.

The published |t| edges of both scripts are the seven bins of
`recopseudo.T_EDGES_PUBLISHED` — 0.017, 0.028, 0.039, 0.05, 0.08, 0.12,
0.17, 0.25 GeV² — adopted 2026-08-28 (Report 2 §7, plans/08 item 4).  The
window they replace is `recopseudo.T_EDGES_LEGACY`, reproducible in either
script by passing it explicitly:

```bash
python3 scripts/money_cos2phi_coherent_reco.py --config 0 --optics tagging \
        --n-mc 6000000 --t-edges 0.05,0.08,0.12,0.17,0.25 --outdir .
python3 scripts/nearbeam_reach_gain.py --n-mc 2000000 \
        --t-edges 0.05,0.08,0.12,0.17,0.25 --outdir .
```

Neither can overwrite a published figure: the stems gain `_tedges`
(`money_cos2phi_coherent_reco_6Li_c0_tagging_tedges.png`,
`nearbeam_reach_gain_6Li_tedges.png`), the sentinel for "published" being
an unset `--t-edges` rather than a value equal to the default.

The fraction of the tagged sample the window keeps — the number Report 2
§7 and plans/08 quote for the adoption — is the sum of the per-bin `N` the
run prints, over the `N_tag` of its header: 78.0 / 88.0 / 85.1% for the
seven published bins at 5 × 40.8 / 10 × 99.5 / 18 × 137.5, against
26.5 / 33.5 / 26.8% for the four of the legacy window.

Each bin's t_ref — the rate-weighted mean TRUE |t| at which a_t is
quoted, and the abscissa of panel (d) — is read off the same lines
without a second run: the published scenario has
c_2(t) = −(P_zz/2) ε_B0 B |t| = 2|t| per unit P_zz (ε_B0 = −0.08,
B = 50 GeV⁻²), so t_ref is half the printed `a_t truth`.  At 5 × 40.8
that gives 0.023, 0.029, 0.036, 0.045, 0.059, 0.074 and 0.090 GeV²: FOUR
of the seven bins, not the three added ones, sit below the 0.05 GeV² at
which `coherent.MANTYSAARI_A2_DEUTERON` begins, the fourth being the
0.05–0.08 bin that resolution pulls down to 0.045 — which is why Report 2
§4.5 says the wider window deepens the anchor extrapolation rather than
introducing it.  (`python3 -c "from polligen import coherent as c;
print(c.CoherentScenario().cos2phi_coefficient_deformation(0.05, 1.0))"`
prints 0.1, i.e. the factor two, and
`print(c.MANTYSAARI_A2_DEUTERON)` prints the anchor's four |t| keys.)

Every 6R run now prints, under each bin's summary line, a per-bin design
diagnostic:

```
design: 46 counts per (alpha, beta) cell, 24 of 24 beta bins populated, cond 1.91, rank 7/7
```

`cond` is √cond of the returned covariance, which is the σ_max/σ_min of
the weighted design that plans/08 item 4 defines, and `rank` is the rank
`reco._harmonic_rank_guard` MEASURED on that weighted design, returned by
both fitters as `design_rank` beside its column count `n_par` (until the
2026-08-28 review the field printed the parameter count twice and could
never have reported a deficiency; a deficient design does not reach this
line at all, because the guard raises a named `LinAlgError` first).  This
line is where the conditioning and β-occupancy statements of Report 2
§4.5 and §7 are read off.  It is unconditional stdout; it changes no
figure and no return value.  On exact counts
(`--exact`, the three configurations) `cond` runs 1.80–6.66 over the
seven published bins and 1.80–2.85 over the four of the legacy window —
2.33 / 1.94 / 1.80 / 2.01 at 5 × 40.8, 2.85 / 2.26 / 1.99 / 2.08 at
10 × 99.5 and 2.52 / 2.04 / 1.85 / 1.93 at 18 × 137.5 — which is the
range plans/08 item 4 quotes; the Poisson draw of the published runs
moves the last digit (1.79–6.65).

How far outside the beam envelope the measured silicon aperture sits —
the number Report 4 §2.1 quotes below Table 2 — is the ratio
of the two cutouts `nearbeam_reach_gain.py` prints per configuration, or,
without rounding:

```bash
python3 - <<'EOF'
import importlib.util
spec = importlib.util.spec_from_file_location(
    "nrg", "scripts/nearbeam_reach_gain.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
from polligen import reco
for i, cfg in enumerate(m._CFG):
    ap, top = reco.rp_aperture_for(cfg), m.ff.tagging_optics_point(cfg)
    print(m.NAMES[i], "%.2f x %.2f" % (ap[0] / top["env_x"],
                                       ap[1] / top["env_y"]))
EOF
```

which gives 7.63 × 2.33, 9.13 × 1.18 and 4.51 × 1.00 envelope widths at
5 × 41, 10 × 100 and 18 × 275: the vertical ratio of 1.00 at the top
configuration is why its silicon row is not identically dead.

**`--fit likelihood`** replaces the bin-wise spin-state ratio by the
acceptance-profiled Poisson likelihood (`reco.harmonic_likelihood_fit_2d`,
plans/08 A12): the same model and the same template basis, with the
per-bin acceptance profiled out, which makes the estimator exactly the
conditional multinomial given the bin totals and therefore unbiased at any
occupancy.  `ratio` is the default and reproduces every published number
bit for bit.  Run it with `--ensemble` to see the difference:

```bash
python3 scripts/money_cos2phi_coherent_reco.py --config 0 --optics tagging \
        --n-mc 6000000 --ensemble 20 --fit likelihood --outdir .
```

**`--u-in-situ`** measures (u₁, u₂) from the spin-averaged counts of the
same data against the response's acceptance shape instead of assuming the
ZEUS values, and propagates their covariance into the harmonic errors
(plans/08 A3).  With the per-bin acceptance free, u is not identifiable at
all, so this necessarily uses the acceptance MC — the point is that ZEUS
becomes the prior it is.  The fitted pair, its error and the term it
propagates into a_e are printed for **every** |t| bin, under that bin's
summary line.

**`--t-edges`** replaces the published reconstructed |t| binning
`0.017,0.028,0.039,0.05,0.08,0.12,0.17,0.25` by any increasing list.  The
window below 0.05 GeV² it once priced is now inside that binning, so what
the flag is for is reproducing the run-13 window (above) and pricing the
one bin that stays out of the published one:

```bash
python3 scripts/money_cos2phi_coherent_reco.py --config 0 --optics tagging \
        --exact --n-mc 6000000 \
        --t-edges 0.006,0.017,0.028,0.039,0.05,0.08,0.12,0.17,0.25 --outdir .
```

That bin, 0.006–0.017 GeV², is excluded on what this run prints (§7 and
plans/08 item 4): 16 of 24 β bins populated at 5 × 40.8 and 12 at the
other two, cond 9.79 / 20.39 / 14.79, and a 10⁻³ envelope split moving its
a_t by −26.5 / −76.2 / −56.4%.  The summary ends with the quadrature
combination of δa_e over the bins that survived — a_e is one constant
across |t|, so its errors combine, and a_t is not.

Three further switches change what the inclusive commands mean.

**`--y-source hfs --hfs-sample <p.npz> <n.npz>`** replaces the 25%
Gaussian stand-in for the hadronic y by a real hadronic final state
through the hadron-side detector response.  Needs §5.1.  Without
`--hfs-calibrate` the pseudo-events carry the response's capture bias on
y_Σ uncalibrated (purities 0.42 / 0.53 / 0.49 / 0.68); with it the
transferred sums are divided by the library's ⟨Σ_reco⟩/⟨Σ_true⟩ in bins
of the *reconstructed* (x_mixed, Q²_e) — the hadronic-scale calibration an
analysis derives from its own simulation and applies at the measured
point (until 2026-08-28 the factor was looked up at the event's true
cell, which no experiment can do) — and the purities are
0.56 / 0.59 / 0.64 / 0.75 at unchanged errors.  `--hfs-scale` then models
a residual scale error: +2% moves Â by 0.2–1.2%.  The published 5R/7R
figures are this run with `--unfold folded --tag _hfscal`.

**`--unfold folded`** replaces the model bin-centering factor K by a Δ(x)
shape fitted *through the response* per Q² slice, and puts the shape-fit,
response-MC and prior-spread terms into the 7R error bars.  `--unfold
model` is the default and reproduces every published number bit for bit.

**`--unfold-scan`** is the closure test behind the claim: it corrects the
same pseudo-data with each of the three Δ priors and prints the residual
bias per bin for both methods.  Run it if you want to see why the folded
fit is worth having.

```bash
python3 scripts/money_cos2phi_reco.py --unfold-scan --outdir .
```

**`--syst-scan`** runs the detector nuisances (electron energy scale,
hadronic-resolution mismatch with common random numbers, ε_eID η tilt,
the Yellow Report EMCal table) and prints their effect on Δ̂.

**`--isr`** is the radiative-correction bound (plans/07 WP4, Report 2 §7).
It rebuilds the response with collinear leading-log initial-state
radiation switched on — `polligen/radiative.py` draws the radiated
fraction z from its own random stream, so the ISR-on and ISR-off
responses sit on identical pseudo-events — and prints purity, efficiency
and the shift of Δ̂ an ISR-free bin-centering would leave, against the 5%
gate, plus the z spectrum and the five-method comparison table.  The hook
is off everywhere else, so no published number moves.

```bash
S=20260824,20260925,20261026,20261127,20261228,20270129,20270302,20270403
python3 scripts/money_cos2phi_reco.py --isr --isr-seeds $S \
        --n-mc-per-cell 1600 --outdir .
python3 scripts/money_cos2phi_reco.py --isr --isr-seeds $S --isr-gen-q2min 0.05 \
        --n-mc-per-cell 1600 --outdir .
python3 scripts/money_cos2phi_reco.py --isr --isr-seeds $S --isr-empz 0.85 \
        --n-mc-per-cell 1600 --outdir .
python3 scripts/money_cos2phi_reco.py --isr --isr-seeds $S --isr-empz 0.85 \
        --y-source hfs --hfs-calibrate --n-mc-per-cell 800 --outdir . \
        --hfs-sample samples/pythia8_e10_p99.5_dis.npz \
                     samples/pythia8_e10_n99.5_dis.npz
```

The two windowed runs (third and fourth) print the retention three ways:
globally, in bands of nominal y, and per analysis bin, since the global
number is dominated by the high-y bulk and is not what a sweet spot pays.

**`--isr-seeds` is not optional for a quoted number.**  One draw of the
response scatters by 4–14% of the bound: over the eight seeds above the
seed-to-seed standard deviation is 0.087 / 0.048 / 0.073 / 0.051
percentage points and the individual draws span 0.51–0.75, 0.43–0.57,
0.80–1.01 and 1.15–1.31% at the four sweet spots, so the default seed's
third spot sits 1.8 σ below the mean.  The published bound is therefore
the mean over that list and the printed `+-` is the standard error of the
mean; the run also prints the spread.  Plain `--isr` at the default seed
prints a single draw, which is what §8 records as the single-command
value.

The first run is the bound with the published generator window: Δ̂ biased
by +0.62 ± 0.03, +0.50 ± 0.02, +0.94 ± 0.03 and +1.22 ± 0.02% at the four
sweet spots.  `--isr-gen-q2min` opens the generator window so that events
below it can radiate into the analysis bins — the feed-in the published
window truncates — and the worst spot then reads 1.87, 2.81, 2.26 and
2.34% at 0.35, 0.15, 0.05 and 0.02 GeV², i.e. the bound saturates in a
1.8–2.8% band rather than at one value, and ≤ 2.9% is what the gate
should be read against.  `--isr-empz 0.85` adds the HERA-style E − p_z
window on both members of the pair and brings it to ≤ 0.25%: the bias
reads +0.23 / +0.16 / +0.22 / +0.18% on the 25% Gaussian y stand-in and
+0.17 / +0.19 / +0.20 / +0.16% through the PYTHIA final state with the
calibrated Σ scale (fourth command), and the window keeps 0.869 / 0.824
(stand-in) and 0.979 / 0.929 (PYTHIA, calibrated Σ) of the non-radiative
/ radiative selected rate; the loss is almost entirely above y = 0.2
(99.5% and 97.2% of the discarded non-radiative rate) and costs
0.01–0.06% at the four sweet spots.  The window is a documented
contingency, not a default: apply it if a published analysis opens the
generator window below Q² = 0.15 GeV², or if the 5% gate tightens.  It is
not the default because the gate already passes at ≤ 2.9% without it,
because it is free where the letter's numbers live, and because what it
does remove is the y > 0.2 rate that carries the low-x end of every Q²
slice in money plots 5R and 7R.

### 4.4 The hadronic-final-state resolution

```bash
python3 scripts/hfs_resolution.py --config 1 \
    --sample samples/pythia8_e10_p99.5_dis.npz samples/pythia8_e10_n99.5_dis.npz --outdir .
```

`--config 0/1/2` selects the low / mid / top beam configuration
(5 × 40.8, 10 × 99.5, 18 × 137.5 GeV/u for ⁶Li — γ-matched at the two
lower configurations, plans/10).  Without `--sample` the
toy string-fragmentation generator is used and every output is labelled
"toy"; its numbers are illustrative and, measured against PYTHIA,
optimistic (§7).

### 4.4b Where the hadronic E − p_z sum goes (Report 2 §3, Figure 2)

```bash
python3 scripts/hfs_acceptance.py --config 1 \
    --sample samples/pythia8_e10_p99.5_dis.npz samples/pythia8_e10_n99.5_dis.npz --outdir .
```

The true Σ_h at the four sweet spots by fate (captured by tracks / photons /
HCal objects; lost forward beyond the calorimeters, below thresholds, or
backward) for the response's calorimeter reach |η| ≤ 3.7 and the ePIC
nominal 4.0, with the Σ-method δy/y through the full response for both.
Expected (16 s): within acceptance and above threshold 0.80 / 0.87 /
0.83 / 0.92 at 3.7 (forward loss 0.17 / 0.08 / 0.15 / 0.07, thresholds
0.03 / 0.06 / 0.02 / 0.02) and 0.85 / 0.89 / 0.88 / 0.94 at 4.0; captured
through the full response 0.70 / 0.74 / 0.74 / 0.85 (median y_Σ/y
0.73 / 0.78 / 0.77 / 0.87); δy/y 0.32 / 0.22 / 0.29 / 0.14 at both reaches
— the forward escape is a scale bias, not the resolution driver.  The
acceptance is applied in the detector frame (25 mrad crossing) since
2026-08-28, worth +0.01 on the captured fractions.

### 4.5 The coherent optics scan (WP5)

```bash
python3 scripts/coherent_optics_scan.py --outdir .
```

Tagged fraction and the one-year error on the deformation coefficient
against the near-beam envelope, for the circular, square and slot
cutouts.  Compare its answer with the *measured* ePIC aperture of §5.3 —
they disagree, and that disagreement is an open item (plans/04 #20).

The last point of `--err-cuts` (0.60 GeV) is set by the response Monte
Carlo, not by data: it tags 18–27 recoils a year, and at `--n-mc 150000`
the harmonic design in its |t| window loses rank.  Since 2026-08-28 that
point is reported and dropped — with the recoil count and the guard's own
message, and a note on panel (c) — instead of aborting the figure; raise
`--n-mc` to recover it.

### 4.5b The tagging optics, priced (Report 1 §6.1)

```bash
python3 scripts/tagging_optics.py --outdir .
```

Per configuration, the tagged fraction and the luminosity against
r = β*_x/β*_x,HA with the horizontal plane de-squeezed alone and the
vertical held at high acceptance (σ_θ ∝ 1/√β*, L ∝ 1/√r; dispersion at
the pots added; parallel-to-point transport assumed), with the pots
following the 10σ envelope or fixed at the measured aperture, the optimum
of their product, and the tagged events, best-bin 5σ floor and shape-term
significance per year at the 10 fb⁻¹/u placeholder.  Expected: optimum
r_h = 49.7 / 175.6 / 89.3 at 5 × 40.8 / 10 × 99.5 / 18 × 137.5 with
ε = 0.423 / 0.323 / 0.332 and L/L_HA = 1/7.1 / 1/13.3 / 1/9.5,
N_tag/yr = 2.6×10⁶ / 3.0×10⁶ / 6.1×10⁶ (2.3 / 5.6 / 3.9× below the 0.20 GeV
reference at 10 fb⁻¹/u), best-super-bin 5σ floors of 1.7 / 2.1 / 1.6% per
unit P_zz, the shape term in the optics' own window 0.031–0.035 per unit
P_zz = 9.3 / 8.3 / 10.7σ per year and 2.8 / 4.4 / 2.6 years to 5σ on a 1%
exotic-glue term — both at the whole 10 fb⁻¹/u year given to this
observable at this optics, so a share f of the year divides the
significance by √f and multiplies the years by 1/f (§3.5, plans/07 WP2)
— and IR-8's ≈ 20% worth 3.3 / 8.2 / 5.7× the optimum at
equal luminosity; both planes de-squeezed gives a fifth of the yield at
1/24–1/70; with the pots fixed nothing is recovered at any β*.

### 4.6 The near-beam study (plans/09)

```bash
python3 scripts/eic_beam_figures.py       --outdir .   # 2 s
python3 scripts/nearbeam_aperture_scan.py --outdir .   # 7 s
python3 scripts/nearbeam_aperture_scan.py --isotope 7Li --outdir .   # 7 s
python3 scripts/nearbeam_reach_gain.py    --outdir .   # 4 s (9 s at the published --n-mc 2000000)
python3 scripts/nearbeam_sensor_budget.py --outdir .   # 1 s
python3 scripts/nearbeam_zid_power.py     --outdir .   # 48 s
python3 scripts/nearbeam_two_hit.py       --outdir .   # 27 s
```

Four questions, four scripts, all on the per-configuration Yellow
Report divergences of plans/10 since 2026-08-28. The first prices *every*
near-beam aperture — coherent tagged fraction and α-tag acceptance against
the horizontal half-width in angle, per configuration, with three markers:
the measured ePIC silicon aperture, the 10σ envelope of the Yellow Report
high-acceptance optics, and the envelope of the lithium tagging optics
(Report 1 §6.1). The second runs the full coherent chain at the tagging
optics of each configuration with the pots fixed at the silicon aperture
and with the pots following the envelope, and reports what the
*measurement* does; since 2026-08-28 it drops a |t| bin before fitting it
whenever the row's own N_tag puts fewer than `MIN_TAGGED_PER_BIN` = 1000
recoils in it, printing `DROPPED` with the count, which is what retires
the one 231-recoil bin the silicon row at 18 × 275 used to fit. The third asks whether a superconducting nanowire can
be the thing that delivers a closer approach: energy deposits in a 12 nm
NbN film, the hot-spot firing-threshold model of charge identification
(Figure 3 of the report), the sizing strip at the tagging envelope, and
the channel count at each available granularity. The fourth (plans/09 B4,
2026-08-28) samples BOTH fragments of a ⁶Li → α + d breakup from one
relative momentum (`spectator.breakup_lab_kinematics`) and asks what hit
multiplicity is worth against the coherent tag: the 2 × 2 topology per
optics, the separation at the pot plane in millimetres and in 500 μm
pixels, and — the result — how often the partner deuteron is there to veto
an α that faked an intact ⁶Li, scanned against the pot's outer edge; the
millimetres are on the per-configuration levers measured in plans/09 B1
and the scan should be read against the measured outer edge, 2.85 / 3.85 /
4.00 mrad, which the banner prints. Its
quantiles come from `--events` breakups (4 × 10⁵) held in memory and every
COUNT from a chunked pass over `--veto-events` (1.2 × 10⁷), because at the
Yellow Report optics the α fakes a coherent tag only once in 10³–10⁴
breakups and a small sample conditions on a handful of events. It
writes `nearbeam_two_hit_6Li.png`, which no report embeds; Report 4
Table 5 and plans/09 §9.2 quote its numbers, and
`tools/consistency_check.py` recomputes the three medians and compares
them with both documents.

The coherent script also grew `--near-beam-mrad`, which replaces the
measured *horizontal* aperture and keeps the measured vertical:

```bash
python3 scripts/money_cos2phi_coherent_reco.py --config 1 --optics high-acceptance \
        --rp-aperture measured --cut-scale-x 1.0 --near-beam-mrad 1.8
```

1.8 mrad is the Yellow Report high-acceptance envelope at 10 × 100 and
0.17 mrad the tagging one.  Asking for the first, as above, returns
nothing: all seven |t| bins are dropped for zero accepted recoils and the
script exits non-zero.  Asking for the second is what `--optics tagging`
already does — acc 0.316, N_tag 2.94 × 10⁶/yr, all seven bins — and it
has to be asked for that way, because `--optics tagging` takes its cutout
from `reco.tagging_optics_point` and therefore forces `--rp-aperture none`
and ignores `--near-beam-mrad`; the switch bites on the legacy and
high-acceptance optics only.  That is not a defect of the switch but the
result of §2 of
Report 4 restated at the reconstructed level — at the Yellow Report
optics the coherent channel has no acceptance at all — and it is why the
published 6R is `--optics tagging` (§4.3).  `--cut-scale-x 1.0` matters
in either case: the default 2.5 comes from the pre-measurement belief in
a wide horizontal slot, and on top of a measured geometric aperture it
imposes a 25σ horizontal retraction that binds *before* the geometry
does — hiding the whole effect.

The ⁷Li scoping memo (plans/09 B3a, 2026-08-28) is arithmetic on the
existing modules rather than a script of its own.  This is that
arithmetic, and it prints every number the memo quotes:

```bash
python3 - <<'PY'
from polligen import coherent as coh, reco, spin, xsec
from polli_fastsim import beams
from polli_fastsim.spectator import nucleus_mass
Q = {"d": 0.2858, "6Li": 0.0806, "7Li": 4.00}      # |Q| [fm^2], plans/01
R = {"d": 2.1305, "6Li": 2.589, "7Li": 2.444}      # R_ch [fm], plans/01
d = {k: Q[k] / R[k] ** 2 for k in Q}; B = coh.gaussian_slope(R["7Li"])
print(d, B, [d["7Li"] / d[k] for k in ("d", "6Li")])
for e in (0.21 * d["7Li"] / d["d"], 0.08 * d["7Li"] / d["6Li"]):
    s = coh.CoherentScenario(slope_b=B, eps_b0=-e)
    print(e, e * B, [float(s.cos2phi_coefficient_deformation(t, 1.0))
                     for t in (0.05, 0.25)])
print([coh.fragment_rigidity(a, z, beam_a=7, beam_z=3)
       for a, z in ((4, 2), (3, 1), (7, 3))],
      [coh.fragment_route_label(7, 3, config=c)
       for c in ("5x41", "10x100", "18x275")],
      coh.fragment_route_label(4, 2),
      [coh.fragment_route_label(a, z, beam_a=7, beam_z=3)
       for a, z in ((4, 2), (3, 1), (7, 3))])
for iso, s in (("6Li", coh.CoherentScenario()),
               ("7Li", coh.CoherentScenario(slope_b=B))):
    print(iso, [float(s.tag_acceptance_angular(reco.sigma_theta_for(c)[0],
                                               c.ion_momentum_per_nucleon,
                                               a_beam=int(iso[0])))
                for c in beams.default_configs(iso)])
m7 = nucleus_mass(3, 7)
print(m7, [2 * (7 * p / m7) * 477.6e-6 * 1e3 for p in (40.8, 99.5, 117.9)])
k = xsec.InclusiveKernel(beams.LI7)
print([k._tensor_moments(m) for m in spin.m_values(1.5)],
      spin.moments_along_axis(1.5, [0.4, 0.1, 0.1, 0.4]))
PY
```

In order: the shape dial Q/R_ch² = 0.0630 / 0.01202 / 0.6697 for
d / ⁶Li / ⁷Li, a factor 10.6 and 55.7 between ⁷Li and the two anchors;
B(⁷Li) = 51.13 GeV⁻²; the two linear rescalings of `eps_b0`, −2.23 from
the deuteron anchor and −4.46 from the ⁶Li default, whose |ΔB₀| = 114
and 228 GeV⁻² exceed B itself and whose cos 2φ coefficients are 2.86 and
14.3, or 5.70 and 28.5, at |t| = 0.05 and 0.25 GeV² — above unity, which
is why the memo declines to quote a ⁷Li deformation at all; the ⁷Li
fragment rigidities 0.856 (α), 1.290 (t) and 1.000 (intact), then what
`fragment_route_label` makes of them while it is left on its ⁶Li
defaults — an intact ⁷Li read as R = 1.166 and routed "lost
(over-rigid)" at 5 × 41, where the pots' blind block is widest, or
"RP-inner (over-rigid)" at the other two, and the ⁶Li α's "RP p_T-tail
only" handed to the ⁷Li α — against the three answers `beam_a=7,
beam_z=3` gives: "RomanPots" for the α, "RP-inner (over-rigid)" for the
triton and beam-blind for the intact recoil; the
intact tagged fractions at the Yellow Report high-acceptance divergences
and 10σ, 5.0×10⁻⁷ / 8.4×10⁻²⁶ / 3.7×10⁻¹³ for ⁶Li — the plans/10 row —
and 1.7×10⁻⁹ / 1.3×10⁻³⁵ / 1.5×10⁻¹⁵ for ⁷Li; m(⁷Li) = 6.5338 GeV and
the maximum lab energy 2γE* = 41.8 / 101.8 / 120.7 MeV of the 477.6 keV
M1 photon; and the J = 3/2 rank-2 moments (Q_NN, c_eff) = (±1, ±3),
whose fill average is 3T exactly (T = 0.600, ⟨c_eff⟩ = 1.800 for
populations 0.4 / 0.1 / 0.1 / 0.4).

---

## 5 · Third-party generators

*Run every command in this section from the repository root.*

### 5.1 PYTHIA 8 — the hadronic final state

`polligen` needs a library of hadronic final states: for a pseudo-event
at (x, Q²) it draws a generator event from the same cell, applies the
hadron-side response, and transfers the captured fraction of
Σ = Σ(E − p_z) and the p_T ratio onto the pseudo-event's exact
kinematics.  The tensor polarization enters only the inclusive azimuthal
weight, so an unpolarized generator is the right source.

With §1.4 done:

```bash
export PYTHONPATH=$HOME/Apps/pythia8311/lib
export PYTHIA8DATA=$HOME/Apps/pythia8311/share/Pythia8/xmldoc
python3 tools/pythia8/gen_dis_hfs.py --target p --n-events 2000000 \
    --electron-energy 10 --p-per-nucleon 99.5 --seed 101 --quiet \
    --out evgen/samples/pythia8_e10_p99.5_dis.npz
python3 tools/pythia8/gen_dis_hfs.py --target n --n-events 2000000 \
    --electron-energy 10 --p-per-nucleon 99.5 --seed 102 --quiet \
    --out evgen/samples/pythia8_e10_n99.5_dis.npz
```

7 000–8 500 events/s, 450 MB and ~140 s per million events at 10 × 99.5.
The standing production is six files (p and n at each of the three
configurations, 8 M events, 3.7 GB; regenerated 2026-08-28 with
`mHatMin = 0.5`, below); `evgen/samples/README.md` is the manifest with
each file's cross section and seed, and the whole set takes five minutes
with six jobs in parallel on eight cores.

The files are git-ignored.  `HFSSample.concatenate` merges the p and n
files with per-event weights σ_gen/n_events, so the library samples the
two targets in the ratio of their cross sections (0.947 against 0.855 μb
at 10 × 99.5) — the per-nucleon luminosity weighting Z = N = 3 asks for
(until 2026-08-28 the merge was by event count).

**Three traps, all handled in the script, all worth knowing.**

*`PhaseSpace:mHatMin` defaults to 4 GeV and applies to DIS.*  The hard
2 → 2 system's invariant mass is m̂² = x s, so the default silently removed
everything below x = 16/s — 0.004 at 10 × 99.5 — which is 39% of the
selected rate of the pseudo-experiments (the sweet spots at x ≥ 0.011 sit
above it; the low-x half of the Q² = 1.14 GeV² slice does not).  The
script now sets `mHatMin = 0.5` (`--mhat-min`); the cross section rose
from 0.666 to 0.947 μb at 10 × 99.5 (p) and the sample's x spectrum
follows the generator's rate map to ±25% down to x = 3×10⁻⁴.

*`PhaseSpace:Q2Min` does nothing on its own.*  PYTHIA applies it only
when `Q2Min ≥ pTHatMinDiverge²`, and `pTHatMinDiverge` defaults to 1 GeV,
so a requested 0.7 GeV² is silently ignored and the divergence cut sets
the floor.  Measured at 10 × 50 over 4000 events (a diagnostic run):
with the default,
min Q² = 1.002 GeV² and σ = 0.381 μb *whatever* is asked for below 1;
with `pTHatMinDiverge = 0.5` (PYTHIA's own lower limit, which the script
sets), min Q² = 0.697 and σ = 0.551 μb.  The 0.7–1.0 GeV² band is 31% of
the sample and is exactly the band the reconstructed-level
pseudo-experiments loosen their generator window to populate.

*The massless Σ identity does not hold for a massive target.*
`hfs.truth_kinematics_check` pins Σ over the final state = 2 E_e y, which
the toy satisfies exactly.  PYTHIA does not: the target nucleon carries
E − p_z = m²/(E + p) = 4.4 MeV at 99.5 GeV/u, so Σ − 2 E_e y = 4.5 MeV —
2% of Σ at y = 0.01.  `HFSLibrary` transfers the *ratio* of the library
event's sums, and since 2026-08-28 `HFSResponse` applies it to the
pseudo-event's Σ *including* the same mass term, so the transfer is
consistent; pass `target_mass=0.938` to the check for a PYTHIA sample.
The p_T identity is exact (≤ 5×10⁻¹¹).

### 5.2 BeAGLE — nuclear breakup, and the e+d control

BeAGLE itself links FLUKA and cannot be built without a personal licence
(`tools/beagle/README.md` has the full dependency recipe for when you
have one).  **The control study does not need it**: the official EIC
BeAGLE samples stream from a public xrootd door, and that is what
calibrates the cluster model's p_T tail — the single most important
model input to the ⁶Li α tag.

```bash
export SIF=~/Projects/eic/local/lib/jug_xl-nightly.sif
singularity exec $SIF xrdfs root://dtn-eic.jlab.org ls /volatile/eic/EPIC/EVGEN/DIS
#   BeAGLE1.03.02-{1.0,1.2,1.3,2.0,2.1,3.0,3.1}, plus Djangoh/pythia6/pythia8 trees
```

`eH2` (deuteron), `eHe3` and `eAu` are available; there is no A = 6 or 7
sample anywhere, which is the gap FLUKA gates.

```bash
B=root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/DIS/BeAGLE1.03.02-3.1/eH2/en/9x130/q2_1to1000
singularity exec $SIF python3 tools/analysis/dump_spectators.py \
    $B/BeAGLE1.03.02-3.1_DIS_eH2_en_9x130_q2_1to1000_ab_run001.hepmc3.tree.root \
    ed.csv --nevents 20000                    # ~2 min, streamed, no download
python3 tools/analysis/ed_control_analysis.py ed.csv --beta-scan --outdir out
```

`en` means DIS on the neutron, so the proton is the spectator.  Expected:
spectator protons in 100.0% of events, routed 93.1% to the off-momentum
detectors (rigidity ratio ≈ 0.5) — the routing logic agrees with BeAGLE
to better than two points in every window.  The p_T tail does not:

```
P(p_T >)                       0.10     0.20     0.30     0.45 GeV
BeAGLE e+d (x_L ∈ [0.9,1.1))  0.1882   0.0494   0.0261   0.0144
Hulthén β = 0.26 (default)    0.1522   0.0201   0.0037   0.0005
Hulthén β = 0.40 (best fit)   0.1934   0.0360   0.0089   0.0015
```

and the scan's conclusion is that **no β reproduces the shape** — BeAGLE
has a narrower core and a harder tail than a two-parameter Hulthén can
have at once.  Since the ⁶Li α tag is entirely a p_T-tail measurement,
its model uncertainty is one-sided *upward*.

Reading the container's HepMC3 tree files needs the *legacy*
`jug_xl-nightly` container; the newer `eic_xl-nightly` pyHepMC3
`rootIO.ReaderRootTree` segfaults on the same files.

### 5.3 ePIC full simulation — far-forward acceptance

**npsim cannot shoot a nucleus.**  `--gun.particle Li6`, `ion(3,6)` and
every other spelling return *"Geant4ParticleGenerator: Bad particle
type"*: DD4hep looks the name up in the G4 particle table, where generic
ions are not pre-instantiated.  Feed it a one-particle HepMC event
instead, which also buys a scan in p_T **and** azimuth:

The container is `~/Projects/eic-2026/local/lib/eic_xl-nightly.sif`
(present on the analysis box since 2026-08-28), whose `epic-main` is git
`9aaa2969`.  **npsim writes `calibrations/`, `fieldmaps/` and `gdml/`
into the current directory**, so `cd` to a scratch directory inside the
`bash -lc` string and never run it from the repository.

```bash
SIF=~/Projects/eic-2026/local/lib/eic_xl-nightly.sif
S=<scratch>                      # NOT the repository

# (a) the azimuthal acceptance map: 18 angles x 12 azimuths, one event each
python3 tools/fullsim/ion_gun_hepmc.py --out $S/gun.hepmc --pdg 1000030060 \
    --a 6 --p-per-nucleon 137.5 \
    --pt 0.2475 0.33 0.4125 0.495 0.5775 0.66 0.7425 0.825 0.9075 0.99 \
         1.11375 1.2375 1.65 2.0625 2.475 2.8875 3.30 4.125 \
    --nphi 12 > $S/gun.index      # theta = 0.3 ... 5.0 mrad; 216 events

# (b) the edges: theta 0.20-6.00 mrad in 0.05 steps at phi = 0/90/180/270
#     (117 angles x 4 azimuths = 468 events).  p_T = theta * A * p_u, so the
#     --pt list is per configuration; A p_u = 825 GeV at 18 x 275 (and
#     244.8 / 597.0 GeV at 5 x 41 / 10 x 100).
python3 tools/fullsim/ion_gun_hepmc.py --out $S/ladder.hepmc --pdg 1000030060 \
    --a 6 --p-per-nucleon 137.5 \
    --pt $(python3 -c "print(' '.join('%.5f' % (825e-3*(0.20+0.05*i)) \
                                      for i in range(117)))") \
    --nphi 4 > $S/ladder.index

# both files go through the same two steps; N is 216 for the map, 468 for
# the ladders, and the npsim caches land in $S
for F in gun:216 ladder:468; do
  NAME=${F%%:*}; N=${F##*:}
  singularity exec $SIF bash -lc "source /opt/detector/epic-main/bin/thisepic.sh
    cd $S                        # <-- the caches land here
    npsim --compactFile \$DETECTOR_PATH/epic_craterlake_18x275.xml \
          --inputFiles $S/$NAME.hepmc --numberOfEvents $N \
          --physics.list FTFP_BERT --part.minimalKineticEnergy '100*MeV' \
          --outputFile $S/$NAME.edm4hep.root"
                                 # ~47 s geometry + ~1.3 s/event

  singularity exec $SIF python3 tools/fullsim/ff_gun_hits.py \
      $S/$NAME.edm4hep.root --per-event --index $S/$NAME.index --positions
done
```

**Read the hits in the rotated frame.**  Both Roman-Pot stations write
into the *same* readout, `ForwardRomanPotHits`, so `--positions` used to
average two planes 1.7 m apart and call it a position; it groups by
station plane since 2026-08-28.  The stations are also rotated by
−0.04545 rad about y, so a plane is not a surface of constant global z:
`ff_gun_hits.py` assigns each hit to the plane whose normal distance is
|w| ≤ 5 mm and reports the in-plane offset, where a ±15 mm window on
global z put every hit beyond |dx| ≈ 110 mm into **both** layers of its
station.

Use the **matching optics file** for the beam energy.  The number here
is a *lattice-matching momentum*, not a fill energy: the gun must put the
⁶Li at the reference rigidity the chosen beamline file was built for, so
that the transport is the one the file describes.  For the two lower
files that rigidity corresponds to 20.5 and 50 GeV/u — the
rigidity-scaled pair that plans/10 A0 retired as ⁶Li *beam* energies in
favour of the γ-matched 40.8 and 99.5 GeV/u — and for the top file to
137.5, where the two conventions coincide because the ring rigidity caps
⁶Li there.  So: `epic_craterlake_5x41.xml` with `--p-per-nucleon 20.5`,
`epic_craterlake_10x100.xml` with 50, `epic_craterlake_18x275.xml` with
137.5.  Plain `epic_craterlake.xml` loads the 5×41 beamline fields, which
sends 275-optics momenta straight to the ZDC.  The angular edges this
scan measures are geometric and independent of that choice; plans/04 #20
converts them into p_T at the γ-matched 40.8 / 99.5 / 137.5 GeV/u, which
is where the lithium fill actually sits (0.61 / 0.90 / 0.44 GeV for the
2.50 / 1.51 / 0.53 mrad edges measured on 2026-08-28).

Expected: no Roman-Pot hits below |θ_x| = 0.55 mrad at 18 × 275 and none
from a *vertical* kick below 0.80 mrad.  An α at the same rigidity
(550 GeV) gives the same edge at the same *angle*, 0.55 mrad, which is
the check that the scan measures optics and not species.
`tools/fullsim/README.md` has the full azimuthal map, the transport
levers regressed from the hit positions, the outer edge, the over-rigid
triton and what all of it means for the coherent channel.

The single-particle routing scan for the other fragments is
`tools/fullsim/ff_gun_scan.sh` (α, triton, deuteron, proton, neutron at
their rigidity ratios).  It is what puts the over-rigid ⁷Li triton on the
Roman-Pot silicon at dx = +66 mm in 60 of 60 events, with a ZDC deposit
behind it.

**Which 5 × 41 lattice.**  `epic-main` ships two beamline files for that
ring setting and they are not a field scale of one another:
`compact/fields/beamline_5x41.xml` (41 GeV proton) and
`beamline_5x41_He4.xml` (Z/A = 0.5, 82 GV, the ⁶Li fill point).  No
`epic_craterlake_5x41_He4.xml` is built, but every include in the compact
file is an absolute `${DETECTOR_PATH}/…`, so the alternative needs one
`sed` and no bind mount:

```bash
singularity exec $SIF bash -lc 'source /opt/detector/epic-main/bin/thisepic.sh
  sed "s|compact/fields/beamline_5x41.xml|compact/fields/beamline_5x41_He4.xml|" \
      $DETECTOR_PATH/epic_craterlake_5x41.xml > '$S'/epic_craterlake_5x41_He4.xml'
```

It is worth R₁₂ = 29.81 m against 19.24 and a 1.61 mrad edge against
2.50 — 48 mm / 29.81 m, the threshold the code applies, the measured first
hit being 1.60 mrad on +x and 1.70 on −x; ×1.55 and ×0.64 — which is the systematic on every 5 × 41
millimetre in the reports (plans/09 D3).

---

## 6 · The reports

```bash
python3 reports/build_report.py          # five self-contained HTML pages,
                                         # numbered 0-4 in reading order
python3 reports/build_report.py --pdf    # and their PDFs (needs §1.6)
```

The builder embeds the figures from `evgen/*.png` as base64 and typesets
the display mathematics with matplotlib mathtext, so the pages have no
external dependency of any kind.  **Regenerate the figures first** (§4)
if you have changed anything that feeds them; the builder fails loudly on
a missing figure rather than silently using a stale one.

---

## 7 · Numbers to check against

Each row is what this machine produced on 2026-08-26.  §8 says how much
each class of number is allowed to move; a disagreement larger than that
is an environment difference, not noise, and worth chasing before you
trust anything downstream of it.

### Fast simulation

| what | command (from `fastsim/`) | expected |
|---|---|---|
| CT18NLO F₂ᵖ(0.1, 10) | see §1.3 | 0.4274 |
| R1998 at (0.1, 5) | `python3 -c "from polli_fastsim.structure import r1998; print(r1998(0.1,5))"` | 0.1844 |
| ⁶Li α-tag, β = 0.30 | `scripts/tagging_acceptance.py` | YR high-acceptance 0.0186 / 0.0168 / 0.0261, tagging optics 0.3541 / 0.2754 / 0.2913, legacy 73 μrad 0.1679 / 0.0262 / 0.0276 (5 × 41 / 10 × 100 / 18 × 275).  The YR high-acceptance column decomposes (the `acc` dict of the same run, in points) as Roman-Pot window 1.51 / 1.51 / 1.54, over-rigid inner branch 0.13 / 0.16 / 1.00, near-beam tail 0.20 / 0.01 / 0.08 and B0 0.02 / — / — — Report 3 Table 6's caption.  Every column rose on 2026-08-28 with the over-rigid branch of `farforward.route_charged`, which routes the R > 1.05 tail to the pots' inner half against the per-configuration blind block 48 / 32 / 16 mm; the same run before the branch existed gave 0.017 / 0.015 / 0.016, 0.353 / 0.273 / 0.281 and 0.167 / 0.025 / 0.018 |
| ⁶Li d-tag, β = 0.30 | same | 0.0930 / 0.0814 / 0.1389 (YR HA), 0.5693 / 0.5108 / 0.5665 (tagging) |
| ⁷Li α-tag, β = 0.30 | same | 0.9690 / 0.9683 / 0.9748 (YR HA), 0.9877 / 0.9919 / 0.9941 (tagging) |
| ⁷Li t-tag, β = 0.30 | same | 0.7777 / 0.9185 / 0.9386 (YR HA), 0.7797 / 0.9218 / 0.9406 (tagging).  This is the row the B1 measurement transformed: the over-rigid triton was routed as lost and read 0.033 / 0.004 / 0.005.  Almost all of it is route 6, RP-inner (0.745 / 0.915 / 0.933), so the tagging optics moves it by under 0.3 points — the tag is dispersive, not angular |
| … on the tagged generator instead | `evgen/scripts/tagged_polarimetry_7li.py` | Roman-Pot tag 0.9620 / 0.9678 / 0.9728 (YR HA), 0.9807 / 0.9917 / 0.9919 (tagging); the script also prints `acc(any far-fwd)` = 0.9699 / 0.9690 / 0.9751 and 0.9886 / 0.9929 / 0.9942, which is the definition the row above tabulates.  The two ⁷Li densities are identical (both pure P wave): the 0.7-point gap at 5 × 41 is the B0, which takes 1.1% of the α there and enters `1 - lost` but not the Roman-Pot mask, and the rest of the `1 - lost` margin is the over-rigid tail `farforward.over_rigid_route` began routing on 2026-08-28, evaluated against the per-configuration blind block since the review of the same day gave the script its `pot_config` (§4.1).  Like for like and inside the tagged model's k ≤ 1.2 GeV/c grid the pure model gives 0.9626 / 0.9683 / 0.9736 — 0.1 point.  For ⁶Li the two differ by the D wave — §4.1 |
| polarized-EMC reach and the CBT–TMT separation | `scripts/money_polemc.py --ion 7Li --pdf grid` | δΔR 0.0423 / 0.0405 / 0.0600 / 0.1893 at x = 0.09 / 0.28 / 0.45 / 0.71 against a separation of 0.044 / 0.040 / 0.034 / 0.059 — 1.04 / 0.98 / 0.57 / 0.31 σ per bin at 10 fb⁻¹/u and 3.3 / 3.1 / 1.8 / 1.0 σ at 100, best in the bin at x = 0.141 (1.16 σ, 3.65 σ) with 5 of 23 bins above 1 σ — but read the window-restricted line instead: best bin x = 0.355 at 0.84 σ and 2.66 σ, 0 of 3 bins above 1 σ, the unrestricted best bin being carried by the transfer (§3.1b).  The same run prints the published-curve separation, max 0.0078 over 0.028 < x < 0.3 against the 0.0398–0.0456 the transferred pair shows there, and max 0.1048 above x = 0.35.  `--pdf toy`, which writes the published PNG: errors 0.0477 / 0.0511 / 0.0619 / 0.1237, best bin x = 0.089 at 0.92 σ, valence window x = 0.355 at 0.72 σ and 2.27 σ.  `--emc-mode constant` returns the pre-digitization separations 0.007 / 0.016 / 0.049 / 0.118 and the old best bin x = 0.562 |
| b₁ signal and reach | `scripts/money_b1.py` | transfer 0.921947 × 2/6 = 0.307316; \|A_zz\| = 2.42×10⁻⁴ / 4.83×10⁻⁴ / 7.79×10⁻⁴ / 1.42×10⁻³ / 1.37×10⁻³ / 3.34×10⁻³ (digitized Miller) and 1.15×10⁻⁵ / 1.15×10⁻⁵ / 6.89×10⁻⁶ / 1.32×10⁻⁶ / 6.94×10⁻⁵ / 4.13×10⁻⁴ (CDKS convolution) at x = 0.005 / 0.01 / 0.03 / 0.07 / 0.2 / 0.5, against δA_zz = 9.31×10⁻⁵ / 8.84×10⁻⁵ / 9.73×10⁻⁵ / 3.55×10⁻⁴ / 1.17×10⁻³ at P_zz = 0.6: 1.7 / 4.8 / 7.6 / 10.9 / 5.6 / 5.8 σ per bin for Miller and ≤ 0.2 σ everywhere for CDKS.  `--transfer legacy` gives 5.0 / 13.7 / 21.6 / 31.0 / 15.9 / 16.5 σ |
| the digitized curves' ratio of effects | the §3.1b snippet | CBT 2.25 / 1.69 / 1.41 / 1.14 against TMT 1.01 / 0.98 / 1.00 / 1.08 at x = 0.40 / 0.45 / 0.50 / 0.60; CBT's minimum 1.06 at x = 0.696; the unpolarized ratios pass through 1 at x = 0.280 and 0.840 (CBT) and 0.721 (TMT), which is why the transfer is a valence-window strength factor (1 for CBT, 0.397009 for TMT); the two published polarized curves agree to better than 0.008 over 0.028 < x < 0.30, and inside the valence window the transferred TMT depletion tracks ⁷Li's own unpolarized 0.034 / 0.048 / 0.087 at x = 0.40 / 0.45 / 0.65 to within 0.005 against CBT's 0.077 / 0.082 / 0.094, so ΔR separates by 0.040 at x = 0.36 and 0.011 at 0.65 |
| A_bag triple (frozen R) | `scripts/money_delta_20260729.py --emit-a-bag-reference` | −0.317767 / −0.310041 / −0.296750 |
| A_bag triple (published R) | `… --r-model r1998 --emit-a-bag-reference` | −0.237040 / −0.235825 / −0.234926 |
| L₅σ (frozen R), any `--run-share` | `scripts/money_delta_realistic.py --configs low,mid,top` | 135.31 / 131.26 / 274.64 fb⁻¹/u (script-internal pre-2026-08-27 configs at 27.5 / 50 / 137.5 GeV/u, superseded by plans/10; only TOP is a machine configuration) |
| L₅σ (published R), any `--run-share` | `… --r-model r1998 --configs low,mid,top` | 67.51 / 65.80 / 155.12 fb⁻¹/u (same caveat) |
| L₅σ toy, current energies | `scripts/money_delta.py` | 16.7 / 16.3 / 21.8 fb⁻¹/u at Δ/F₁ = 10⁻³, P_zz = 0.8 (16.719 / 16.332 / 21.811 to the precision the 2026-08-28 min-events correction is pinned at, `tests/test_money_delta_mask.py`).  Unchanged at `--run-share 0.25`, where the error on the Δ/F₁ scale after one programme year doubles instead: 2.586 / 2.556 / 2.954 → 5.172 / 5.112 / 5.907 ×10⁻⁴ (§3.5) |

### Event generator

| what | command (from `evgen/`) | expected |
|---|---|---|
| 5R sweet spots (x, Q²) | `scripts/money_cos2phi.py` | (0.028, 1.14), (0.011, 1.14), (0.071, 3.13), (0.141, 14.3); A = 7.4 / 4.4 / 9.5 / 9.5 ×10⁻³, δA = 1.7 / 1.4 / 2.8 / 4.5 ×10⁻⁴ (1 yr) |
| 5R sweet-spot purity, 25% stand-in | `scripts/money_cos2phi_reco.py` | 0.66 / 0.63 / 0.69 / 0.69 (D = 0.92 / 0.99 / 0.90 / 0.96); δÂ = 1.2 / 0.9 / 1.6 / 2.9 ×10⁻⁴ |
| 5R sweet-spot purity, PYTHIA HFS, uncalibrated | `… --y-source hfs --hfs-sample …` | 0.42 / 0.53 / 0.49 / 0.68 |
| 5R amplitude dilution D, PYTHIA HFS, uncalibrated | same | 0.79 / 0.84 / 0.83 / 0.95 |
| 5R with the hadronic scale calibrated in reconstructed bins | `… --y-source hfs --hfs-sample … --hfs-calibrate` | purity 0.56 / 0.59 / 0.64 / 0.75, D = 1.00 / 1.06 / 0.95 / 0.98, δÂ = 1.2 / 0.9 / 1.7 / 3.0 ×10⁻⁴ |
| 7R folded-fit bars at the best bins, PYTHIA calibrated | `… --hfs-calibrate --unfold folded` | 4.2 / 2.6 / 7.0% (1 yr), 3.8 / 2.2 / 2.6% (10 yr) at Q² = 1.14 / 3.13 / 14.3; prior spread 3.7 / 2.1 / 1.3% |
| residual hadronic scale +2% (calibrated) | `… --hfs-calibrate --hfs-scale 1.02` | Â moves +0.2 to +1.2% |
| Σ at the mid sweet spots | `scripts/hfs_acceptance.py --config 1 --sample …` | within acceptance and above threshold 0.80 / 0.87 / 0.83 / 0.92 at |η| ≤ 3.7 (forward loss 0.17 / 0.08 / 0.15 / 0.07), 0.85 / 0.89 / 0.88 / 0.94 at 4.0; captured through the response 0.70 / 0.74 / 0.74 / 0.85; δy/y unchanged between reaches |
| Σ-method δy/y, PYTHIA, mid config (its own sweet spots) | `scripts/hfs_resolution.py --config 1 --sample …` | 0.32 / 0.22 / 0.29 / 0.15 |
| … low config | `--config 0` | 0.39 / 0.23 / 0.24 / 0.11 |
| … top config | `--config 2` | 0.23 / 0.19 / 0.21 / 0.18 |
| unfolding model dependence (moment_B prior) | `scripts/money_cos2phi_reco.py --unfold-scan` | bin-by-bin (−4.1, +8.0, −5.5, +4.9)% → folded (−1.5, −1.8, −0.3, +0.5)% |
| the γ² (target-mass) bound on A_∥ | `scripts/target_mass_bound.py` | cap 0.0965; grid maxima 0.0854 / 0.0577 / 0.0258 (⁶Li) and 0.0854 / 0.0577 / 0.0332 (⁷Li); sweet spots max γ² 0.00564 and max A_∥ shift 0.56% at 10 × 99.5 (0.49% at 18 × 137.5, 2.46% at 5 × 40.8), where the lab-angle azimuth shortcut errs by at most 1.38 / 1.21 / 5.94 mrad; polarized-EMC ΔR biased high by 0.118 / 0.444 / 0.709 / 1.055% at x = 0.089 / 0.282 / 0.447 / 0.708 against the same weights' δΔR = 0.0477 / 0.0511 / 0.0619 / 0.1237 (the toy inputs the block prints; the grid-input δΔR of Report 0's Table 3 is 0.042 / 0.041 / 0.060 / 0.189), the rate-weighted ⟨γ²⟩ over the same window running 0.001279 / 0.004406 / 0.006850 / 0.009903, i.e. ≤ 0.010; tagged-triton overlay ≤ 2.09% at the published configuration and 5.04% at 5 × 40.8 |
| collinear-ISR migration bound (plans/07 WP4) | `scripts/money_cos2phi_reco.py --isr --n-mc-per-cell 1600` | one draw at the default seed 20260824: Δ̂ biased by +0.616 / +0.501 / +0.804 / +1.240% at the four sweet spots; ⟨z⟩ = 0.0245, 26.6% of the rate radiates above z = 10⁻⁴, ⟨z\|z>10⁻⁴⟩ = 0.0923, t = 0.0700 at ⟨Q²⟩ = 4.37 GeV²; the covariant azimuth's residual under k → (1−z)k is 2.6×10⁻² rad and fakes cos 2φ′ at 9×10⁻⁸ |
| … the PUBLISHED bound, averaged over eight response seeds (`$S` below is that list) | `… --isr --isr-seeds $S --n-mc-per-cell 1600`, `S=20260824,20260925,20261026,20261127,20261228,20270129,20270302,20270403` | +0.62 ± 0.03 / +0.50 ± 0.02 / +0.94 ± 0.03 / +1.22 ± 0.02%; purity 0.653 → 0.638, 0.633 → 0.613, 0.679 → 0.659, 0.684 → 0.640; efficiency 0.414 → 0.404, 0.590 → 0.572, 0.374 → 0.369, 0.653 → 0.634 |
| … with the low-Q² feed-in opened up | `… --isr --isr-seeds $S --isr-gen-q2min 0.05` | +2.26 ± 0.03 / +2.24 ± 0.10 / +1.88 ± 0.07 / +1.88 ± 0.05%.  Worst spot over the window sequence 0.7 → 0.35 → 0.15 → 0.05 → 0.02 GeV²: 1.22 → 1.87 → 2.81 → 2.26 → 2.34%, so the bound saturates in a 1.8–2.8% band and **≤ 2.9%** is what the ≤5% gate is read against |
| … behind a HERA-style E − p_z window | `… --isr --isr-seeds $S --isr-empz 0.85` | +0.23 ± 0.02 / +0.16 ± 0.02 / +0.22 ± 0.02 / +0.18 ± 0.01%, independent of the generator window; the window keeps 0.8691 / 0.8235 of the non-radiative / radiative selected rate on the 25% Gaussian y stand-in and 0.9786 / 0.9287 through the PYTHIA final state with the calibrated Σ scale; per sweet spot the non-radiative retention is 0.9999 / 0.9996 / 0.9999 / 0.9995 and 0.9997 / 0.9995 / 0.9999 / 0.9994, i.e. a cost of 0.01–0.06%; retention by nominal y (0–0.05 / 0.05–0.2 / 0.2–0.5 / 0.5–1) 1.0000 / 0.9979 / 0.7800 / 0.2656 and 1.0000 / 0.9982 / 0.9727 / 0.8795 |
| … through the PYTHIA hadronic final state | `… --isr --isr-seeds $S --y-source hfs --hfs-sample … --hfs-calibrate --n-mc-per-cell 800` | +0.38 ± 0.03 / +0.44 ± 0.02 / +0.46 ± 0.03 / +0.88 ± 0.04% — the bound does not come from the Gaussian y stand-in |
| … that final state behind the E − p_z window | `… --isr --isr-seeds $S --isr-empz 0.85 --y-source hfs --hfs-sample … --hfs-calibrate --n-mc-per-cell 800` | +0.17 ± 0.03 / +0.19 ± 0.03 / +0.20 ± 0.03 / +0.16 ± 0.02% (retentions in the E − p_z row above) |
| coherent tagged fraction | `scripts/coherent_optics_scan.py` | 32% / 3.0% / 4×10⁻⁵ / 2×10⁻⁷ at 0.10 / 0.22 / 0.45 / 0.60 GeV.  The 0.60 GeV point is response-MC-limited (18–27 tagged/yr): at `--n-mc 150000` its harmonic design loses rank, and the point is now reported and dropped rather than aborting the figure |
| coherent tag at the three half-widths | `scripts/nearbeam_aperture_scan.py` | silicon / YR HA envelope / tagging envelope: 9.4×10⁻¹⁰ / 7.2×10⁻⁸ / 0.42 (5 × 41), 2.0×10⁻¹⁹ / 6.2×10⁻²⁷ / 0.32 (10 × 100), 1.2×10⁻⁵ / 7.1×10⁻¹⁴ / 0.33 (18 × 275).  The silicon column is the aperture re-measured on 2026-08-28 (2.50 / 1.51 / 0.53 mrad against 2.00 / 1.35 / 1.03), and the 18 × 275 YR column moved with it because that column takes the vertical at max(measured, 10σ_v) and the measured c_y fell from 2.3 to 0.92 mrad |
| ⁶Li α tag (routed) at the three half-widths | same | 0.0170 / 0.0177 / 0.3556, 0.0163 / 0.0162 / 0.2771, 0.0289 / 0.0247 / 0.2920 |
| ⁷Li α tag at the same three | `… --isotope 7Li` | 0.9668 / 0.9684 / 0.9883 (5 × 41), 0.9682 / 0.9672 / 0.9922 (10 × 100), 0.9812 / 0.9747 / 0.9943 (18 × 275) — flat to three points across the whole 0.05–3 mrad axis, at 1/8.2 / 1/15.3 / 1/10.1 of the luminosity for the tagging point.  Panel (a) stays the coherent intact ⁶Li at either setting: `polligen.coherent` is ⁶Li-specific and a ⁷Li coherent channel is a different amplitude (plans/09 B3, open) |
| money plot 4: ⁶Li α-tag reach, 10 × 99.5 | `scripts/money_tagged_azz.py --outdir . --events 400000` | acc 0.0247 (YR HA) vs 0.3061 (tagging, L/L_HA = 1/13.3) — acc × L 0.0247 vs 0.0231, a 7% cost — but median accepted k 0.323 vs 0.162 GeV/c and frac(k < 0.15) 0.000 vs 0.438.  At k ≈ 0.325 GeV/c A_zz = +0.491 (acceptance-weighted truth +0.455) and −0.070 (−0.098); the θ_k = 90° curve says −0.482 at both.  `acc` is the unbinned accepted fraction: 9–11% of the accepted α lie above the 0.6 GeV/c right edge of the panel.  The kernel forms b₁/F₁, and the digitized Miller b₁ that `polli_fastsim.polarized.toy_b1` returns by default is tapered above its last point (x = 0.9) with the (1 − x)³ falloff of F₁ (`polarized._interp_tapered`, 2026-08-29): frozen at its x = 0.9 value it made the ratio diverge at the generator grid's x = 0.9550 cell and raised `ValueError: negative phi-averaged density` (min(1 + w_avg) = −2.32 in 22 of 2750 cells).  b₁ enters this figure through w_avg alone, so the row's numbers do not depend on it at the digits quoted; a ⁶Li b₁ is read off `money_b1.py`, not off this generator |
| … the same at the other two configurations | `… --config 0` / `--config 2` | acc 0.0284 / 0.0265 (YR HA) vs 0.3816 / 0.3127 (tagging); median accepted k 0.348 / 0.334 GeV/c (YR HA, frac below 0.15 GeV/c = 0.000) against 0.146 / 0.161 with 0.522 / 0.446 (tagging); acc × L/L_HA 0.0284 → 0.0541 and 0.0265 → 0.0331, i.e. the tagging optics gains reach *and* rate at both, and costs rate only at 10 × 100 |
| ⁷Li polarimetry and tagged EMC, 10 × 99.5 | `scripts/tagged_polarimetry_7li.py` | acc(RP) 0.9678 (YR HA) vs 0.9917 (tagging); ⟨P₂⟩ slope −0.1947 vs −0.1964 against the analytic −0.2000 (legacy 73 μrad: −0.1929); median δA_∥ 0.01152 vs 0.01140 *at equal generated statistics* — the plotted bars are drawn there, and the figure says so.  At equal luminosity the tagging optics multiplies every ⁷Li error bar by 2.83 / 3.87 / 3.15 — a factor 8–15 net loss |
| Cosyn–Weiss deuteron limit (`plans/05` §5.4) | `python3 -m pytest tests/test_tagged.py::test_cosyn_weiss_tensor_gate -q` | A_zz^wf / P₂(cos θ_k) = 0.99940 at k = 0.3012 GeV/c in every angular cell away from the P₂ zero (cells with \|P₂\| ≤ 10⁻³ excluded), spread < 10⁻⁵; radial envelope peak at k = 0.3098 GeV/c against Cosyn–Weiss's 0.30 GeV/c for AV18; A_T∥ = −2 A_zz^wf extremes +0.99967 and −1.93782 at the outermost cells, +1.000 / −2.000 extrapolated, against their TABLE II's +1 and −2 |
| the chain at the tagging optics | `scripts/nearbeam_reach_gain.py --n-mc 2000000` | pots at the silicon: acc 0 and 0 of 7 bins at 5 × 41 and 10 × 100, acc 1.45×10⁻⁵ with 268 tagged/yr and 0 of 7 at 18 × 275 — its one populated bin, 231 expected tagged recoils in 0.17–0.25 GeV², is printed as `DROPPED` by the minimum-count guard `nearbeam_reach_gain.MIN_TAGGED_PER_BIN` = 1000 (2026-08-28), which is below two counts per (α, β) cell of the 12 × 24 design in each of the two fills; the fit that bin used to get returned a_t = −1.56 ± 2.22 at a t_ref outside the linear model, which is a count and not a measurement and alone set the vertical range of panel (a).  A dropped bin appears in neither panel and is counted out of the header's surviving-bin tally, and the empty-bin `RuntimeWarning` from `truth_reference` goes with it.  Nothing else in the run moves: pots following: acc 0.411 / 0.316 / 0.325, N_tag 2.52 / 2.94 / 6.00 ×10⁶/yr, 7 of 7 \|t\| bins, δa_t 0.0053 / 0.0047 / 0.0047 / 0.0035 / 0.0048 / 0.0086 / 0.0180 (5 × 41), 0.0074 / 0.0054 / 0.0051 / 0.0034 / 0.0042 / 0.0070 / 0.0135 (10 × 100), 0.0039 / 0.0031 / 0.0032 / 0.0023 / 0.0033 / 0.0064 / 0.0153 (18 × 275).  `--fit likelihood` leaves the errors alone (0.0035–0.0188 / 0.0034–0.0138 / 0.0023–0.0155) and moves only the recovered a_t of the sparsest bin, 0.1230 → 0.1932 against 0.1804 injected at 5 × 41 |
| … the run-13 \|t\| window on the same transport | `scripts/nearbeam_reach_gain.py --n-mc 2000000 --t-edges 0.05,0.08,0.12,0.17,0.25` | writes the guarded stem `nearbeam_reach_gain_6Li_tedges.png`; δa_t 0.0035 / 0.0048 / 0.0086 / 0.0181 (5 × 41), 0.0034 / 0.0042 / 0.0069 / 0.0134 (10 × 100), 0.0023 / 0.0033 / 0.0064 / 0.0149 (18 × 275) |
| sizing strip at the tagging envelope | `scripts/nearbeam_sensor_budget.py` | d50 / d90 / d99 = 183 / 504 / 843 μrad (5 × 41), 69 / 194 / 328 (10 × 100), 50 / 141 / 239 (18 × 275); α at 137.5 GeV/u 77 / 270 / 624.  The strip does not move with the aperture and that is the point of it: it measures the gap between the silicon and the *tagging* envelope, and the tagging envelope is the near edge.  What moves is how much strip there is to build — 12.4 mm at 18 × 275 (0.53 → 0.12 mrad at R₁₂ = 30.0 m) → 496 mm², 3968 microwire channels, 551 111 nanowire devices, against the 9.3 mm the retired 72.7 μrad envelope gave |
| tagging optics, priced | `scripts/tagging_optics.py` | horizontal-only optimum β*_x/β*_x,HA = 49.7 / 175.6 / 89.3, ε = 0.423 / 0.323 / 0.332, L/L_HA = 1/7.1 / 1/13.3 / 1/9.5, N_tag/yr = 2.6×10⁶ / 3.0×10⁶ / 6.1×10⁶, 5σ floor/yr = 1.7 / 2.1 / 1.6% per unit P_zz, shape term 9.3 / 8.3 / 10.7σ/yr.  None of it moved with the 2026-08-28 re-measurement, because this quantity reaches the pot dispersion only through the ratio D/R₁₂ and the re-measured 18 × 275 pair, 0.292 / 29.97, lands within 0.6% of the 0.30 / 30.6 it was priced with.  The banner's aperture line does move, to 2.50/8.84 (the vertical plane shut), 1.51/2.12 and 0.53/0.92 mrad, and with it the "pots FIXED at the aperture" ε, 9.36×10⁻¹⁰ / 1.97×10⁻¹⁹ / 1.23×10⁻⁵ |
| … on the per-configuration levers | `… --per-config-levers` | D/R₁₂ is 1.62×10⁻² / 1.35×10⁻² / 9.74×10⁻³ m/m, the two lower rows 66% and 39% larger than the single 18 × 275 ratio, and the optimum moves: β*_x/β*_x,HA = 46.5 / 164.1 / 89.3, ε = 0.374 / 0.251 / 0.332, L/L_HA = 1/6.8 / 1/12.8 / 1/9.5, N_tag/yr = 2.37×10⁶ / 2.42×10⁶ / 6.15×10⁶, 3.0 / 5.5 / 2.6 years to 5σ on a 1% exotic-glue term.  Priced and **opt-in**: the 22% drop in ε at 10 × 100 would propagate through `money_cos2phi_coherent_reco.py` into Reports 1, 3 and 4, and the coherent chain has not been re-run on it |
| 6R at the tagging optics, 5 × 40.8 | `scripts/money_cos2phi_coherent_reco.py --config 0 --optics tagging --n-mc 6000000 [--ensemble 20]` | σ_θ = 33/380 μrad, cutout 0.33 × 3.8 mrad, acc 0.4111, N_tag 2.52×10⁶/yr at L/L_HA = 1/7.1, ⟨cos 2β⟩ = −0.274; a_t 0.0466 ± 0.0053 / 0.0563 ± 0.0047 / 0.0672 ± 0.0048 / 0.0915 ± 0.0035 / 0.1099 ± 0.0048 / 0.1384 ± 0.0085 / 0.1215 ± 0.0174 (1 yr; inj. 0.045 / 0.059 / 0.071 / 0.090 / 0.118 / 0.147 / 0.180), 0.0441 / 0.0594 / 0.0713 / 0.0885 / 0.1190 / 0.1457 / 0.1763 at 10 yr; ensemble means 0.0440 / 0.0587 / 0.0699 / 0.0893 / 0.1174 / 0.1413 / 0.1157; combined one-year δa_e **0.00119**.  `--fit likelihood` on the same single draw returns the sparsest bin as **0.1842 ± 0.0183** against 0.1803 injected (+0.2σ) where the ratio estimator gives 0.1215 ± 0.0174 (−3.4σ), the combined δa_e being unchanged at 0.00119 (Report 1 §6.3).  The `design:` lines read 2011 / 1472 / 1031 / 1462 / 619 / 189 / 46 counts per (α, β) cell, 20 / 20 / 20 / 20 / 24 / 24 / 24 of 24 β bins populated, cond 4.24 / 3.24 / 2.76 / 2.33 / 1.94 / 1.79 / 1.91, rank 7/7 throughout |
| … the same with the likelihood estimator | `… --ensemble 20 --fit likelihood` | ensemble means 0.0432 / 0.0585 / 0.0700 / 0.0896 / 0.1179 / 0.1476 / 0.1755 (σ 0.0046 / 0.0049 / 0.0052 / 0.0031 / 0.0043 / 0.0072 / 0.0180) on the same twenty draws.  With `--ensemble 200` (~14 min on the seven-bin window): likelihood 0.0448 / 0.0588 / 0.0712 / 0.0901 / 0.1182 / 0.1476 / 0.1824, pulls of the mean −0.8 / +0.5 / +0.1 / +1.4 / +0.1 / +0.5 / +1.6, spreads 0.0055 / 0.0044 / 0.0048 / 0.0031 / 0.0049 / 0.0089 / 0.0193 against quoted errors 0.0053 / 0.0047 / 0.0048 / 0.0035 / 0.0048 / 0.0086 / 0.0187; ratio 0.0455 / 0.0590 / 0.0710 / 0.0898 / 0.1176 / 0.1413 / 0.1191, pulls +1.2 / +1.1 / −0.4 / −0.0 / −1.4 / −9.3 / −41.5 |
| … coarser (α, β) bins, same 200 draws | `… --ensemble 200 --n-alpha 8 --n-beta 16` (and `--n-alpha 6 --n-beta 12`) | the sparsest bin's ratio mean rises 0.1191 → 0.1549 → 0.1666 against 0.1803 injected, i.e. −33.9% → −14.1% → −7.6%, at δa_e 0.0139 → 0.0148 → 0.0162 in that bin and 0.00119 → 0.00126 → 0.00138 combined (+6%, +16%).  Coarser binning attenuates the bias, it does not remove it |
| 6R at the tagging optics, 18 × 137.5 | `… --config 2 --optics tagging --n-mc 6000000` | acc 0.3247, N_tag 6.00×10⁶/yr at L/L_HA = 1/9.5; a_t 0.0434 ± 0.0039 / 0.0637 ± 0.0031 / 0.0742 ± 0.0032 / 0.0985 ± 0.0023 / 0.1403 ± 0.0033 / 0.1666 ± 0.0063 / 0.1902 ± 0.0148 (inj. 0.047 / 0.062 / 0.077 / 0.100 / 0.137 / 0.179 / 0.228); ensemble means 0.0463 / 0.0610 / 0.0771 / 0.1001 / 0.1368 / 0.1776 / 0.2060 (ratio) and 0.0461 / 0.0606 / 0.0770 / 0.1002 / 0.1372 / 0.1787 / 0.2289 (`--fit likelihood`, same draws); combined δa_e **0.00074** |
| 6R at the tagging optics, 10 × 99.5 | `… --config 1 --optics tagging --n-mc 6000000 --ensemble 20 --fit likelihood` | acc 0.3158, N_tag 2.94×10⁶/yr at L/L_HA = 1/13.3; ensemble means 0.0468 / 0.0607 / 0.0732 / 0.0905 / 0.1138 / 0.1348 / 0.1582 against inj. 0.0481 / 0.0611 / 0.0728 / 0.0897 / 0.1135 / 0.1360 / 0.1594; combined δa_e **0.00104** |
| in-situ (u₁, u₂), tagging optics | `… --config 0 --optics tagging --exact --u-in-situ --n-mc 6000000` | δu₂ = 0.0029 / 0.0031 / 0.0035 / 0.0028 / 0.0040 / 0.0070 / 0.0136 per \|t\| bin at 5 × 40.8 in one year (0.0016 / 0.0017 / 0.0020 / 0.0016 / 0.0026 / 0.0049 / 0.0111 at 18 × 137.5, `--config 2`), i.e. 1.8–15× inside the ZEUS 1σ of 0.024; the propagated a_e term is 0.00002–0.00063 across the fourteen bins, negligible against the 0.00134–0.01406 statistical error |
| the 0.006–0.017 bin that stays out | `… --config {0,1,2} --optics tagging --exact --n-mc 6000000 --t-edges 0.006,0.017,0.028,0.039,0.05,0.08,0.12,0.17,0.25` (and the same with `--no-sin --envelope-split 1e-3 --split-axis x`) | that bin has 16 / 12 / 12 of 24 β bins populated, cond 9.79 / 20.39 / 14.79 and δa_t 0.0132 / 0.0328 / 0.0138, and a 10⁻³ envelope split moves its a_t by −26.5 / −76.2 / −56.4%; the printed combination over the eight bins is 0.00105 / 0.00097 / 0.00068 against 0.00119 / 0.00104 / 0.00074 over the seven published ones |
| 6R systematics at 5 × 40.8, exact counts | `… --exact --no-sin --envelope-split 1e-3 / --u2-assumed 0.044 / --rel-lumi-offset 1e-3` | a_t −5.5 / −2.2 / −1.4 / −0.8 / −0.5 / −0.3 / −0.4% at 1e-3 (−54.1 / −25.3 / −13.8 / −7.7 / −4.6 / −3.6 / −3.1% at 1e-2), against a worst bin of −12.5% (−131.0%) at 10 × 99.5 and −10.1% (−99.4%) at 18 × 137.5; a_e +0.0003 to +0.0009 over the seven bins at 5 × 40.8 (3–9%) and, with `--config 2 --optics tagging --exact --no-sin --u2-assumed 0.044 --n-mc 6000000`, +0.0004 to +0.0013 (4–13%) at 18 × 137.5 — fits 0.0104 / 0.0105 / 0.0106 / 0.0108 / 0.0109 / 0.0111 / 0.0113 against 0.0100 injected; a_t ≤ 0.2%, its printed resolution |
| α + d separation at the pots | `scripts/nearbeam_two_hit.py` | median 25.8 / 10.7 / 10.9 mm (16–84% 11.4–53.4 / 3.4–26.4 / 3.2–26.9) = 52 / 21 / 22 pixels of 500 μm at 5 × 41 / 10 × 100 / 18 × 275; the angular lever alone gives 23.1 / 6.2 / 6.2, the rest is the per-configuration pot dispersion D = 0.311 / 0.287 / 0.292 m acting on rigidities that move apart with k_z, a median 9.5 / 8.8 / 8.9 mm; θ_d/θ_α → 1.987 as k → 0.  The three medians fell by up to 42% on 2026-08-28 because R₃₄ was measured — 2.9–3.4 m against the R₁₂ it had been taken equal to — so the vertical half of the separation collapses.  The `--beta 0.20 / 0.40` scan has not been re-run on the per-configuration levers; the −12% to +9% spread is expected to survive but the 33.6 / 16.2 / 13.2 and 41.7 / 20.1 / 16.4 mm it gave on the single lever are superseded, and the veto still moves by less than 0.03 |
| α + d topology per breakup | same | both fragments recorded 0.0002 / 0.0000 / 0.0002 (YR high acceptance) and 0.2874 / 0.2170 / 0.2234 (tagging); d alone 0.075 / 0.061 / 0.068 and α alone 0.017 / 0.015 / 0.016 (YR); a recorded pair lands inside one 500 μm pixel in 5.2 × 10⁻⁵ / 2.7 × 10⁻³ / 4.2 × 10⁻³ of cases at the tagging optics and never at the YR or legacy envelopes — a factor 10–25 more than the single-lever arithmetic gave, through the same collapse of the vertical lever |
| partner-fragment veto | same | P(α fakes a coherent tag) 0.0020 / 0.0001 / 0.0007 (YR) and 0.338 / 0.259 / 0.266 (tagging); of those the partner d is recorded 0.120 ± 0.002 / 0.017 ± 0.004 / 0.245 ± 0.005 (YR, on 2.4 × 10⁴ / 1.4 × 10³ / 8.7 × 10³ fakes) and 0.850 / 0.839 / 0.840 (tagging), falling to 0.00 / 0.31 / 0.56 against a 0.5 mrad outer edge — which is 10 / 11 / 15 mm at the per-configuration R₁₂ and not one number.  The MEASURED outer edge is 2.85 / 3.85 / 4.00 mrad (`farforward.THETA_RP_OUTER_MEASURED`), not the 5 mrad the acceptance tables assume, and at it the veto is 0.80 / 0.84 / 0.84 — the script's `OUTER_SCAN` grid does not contain it, so that triple comes from the same counting pass run with `theta_outer = ff.theta_rp_outer_for(key)`: `route_charged` on both fragments of 1.2 × 10⁷ `spectator.breakup_lab_kinematics` breakups per configuration at seed 7 and β = 0.30, tagging optics, veto = P(route(d) ∈ {1, 4} \| route(α) = 4) |
| hot-spot Z-ID thresholds | `scripts/nearbeam_sensor_budget.py` | r_s = 134 / 268 / 402 nm for p,d / α / ⁶Li; at w = 1 µm, I_th/I_c = 0.73 / 0.46 / 0.20 |
| Z-ID fake rate, 4 planes at 95% eff | `scripts/nearbeam_zid_power.py` | 2.3×10⁻⁵ (8-bit LLR) / 3.1×10⁻⁵ (one bit) / 2.7×10⁻³ (truncated mean) / 5.3×10⁻² (plain sum); 50% fill cannot reach 95% |

### Third-party

| what | expected |
|---|---|
| PYTHIA σ_gen, e+p at 10 × 99.5 | 0.9473 μb (n: 0.8551), with mHatMin = 0.5 |
| PYTHIA sample, 2 M events | 25.3 M particles, 910 MB, ~270 s |
| BeAGLE e+d, P(p_T > 0.3) in the x_L peak | 0.0261, against the Hulthén model's 0.0037 |
| ⁶Li Roman-Pot edge, 18 × 275 optics | 0.55 mrad on both horizontal sides and 0.80–0.95 mrad vertical (0.05 mrad ladder, `epic-main` git 9aaa2969, 2026-08-28) |

---

## 8 · Determinism, seeds, and how much a number may move

Every Monte-Carlo script takes `--seed` and defaults to a fixed one, so a
rerun on the same machine reproduces its own output exactly.  Across
machines and library versions:

- **Analytic results** (`fastsim` rates, FOMs, `rp_hole_acceptance`,
  R1998, the structure functions) are deterministic to floating-point
  round-off.  A disagreement is a bug, not noise.
- **Pseudo-experiments** (`money_cos2phi_reco`, `money_tagged_azz`,
  `closure_fom`) carry the Poisson draw.  Fitted amplitudes move by their
  quoted statistical error; the *truth* references printed beside them do
  not move at all.
- **The response MC** (`RecoResponse`, 400 pseudo-events per sampler
  cell) is the floor on everything reconstructed: 0.13–0.21% on Δ̂ at the
  sweet spots, and 0.10–1.00% per plotted 7R bin.  Two runs at different
  response seeds differ by that much; anything smaller is not a result.
  This is why the systematic scans use **common random numbers** across
  the two responses they compare.
- **PYTHIA** is seeded and reproducible per run; different seeds move a
  resolution number in the third digit.

Where a comparison must be sharper than the MC floor, the scripts hold
the response fixed and vary one thing — that is what `--syst-scan` and
`--unfold-scan` do.

---

## 9 · Troubleshooting

Every entry below is an error we actually hit, with its cause.

| symptom | cause and fix |
|---|---|
| `ModuleNotFoundError: pythia8` inside eic-shell | the container has the C++ library only; build the bindings (§1.4) |
| `AttributeError: module 'numpy' has no attribute 'trapezoid'` | numpy < 2; the scripts carry a shim, so this means an old copy — pull |
| `UnicodeEncodeError: 'gbk' codec` at the third banner line | a non-UTF-8 console; `money_delta_20260729.py` switches the stream error handler itself, so this means an old copy. Reproduce deliberately with `PYTHONIOENCODING=gbk` |
| `FileNotFoundError: fastsim/scripts/money_delta_20260729.py` from another directory | `_check_reco_mask_invariants.py` now defaults to its own sibling; an old copy assumed the repository root |
| PYTHIA `sigmaGen()` returns 0 | `infoPython()` returns the `Info` object **by value**; a copy taken before the event loop is a stale snapshot. Re-fetch it at the end |
| PYTHIA segfaults after many `infoPython()` calls | same cause; do not call it per event |
| min Q² is 1.0 whatever `Q2Min` says | `pTHatMinDiverge` defaults to 1 GeV and gates the cut (§5.1) |
| `Geant4ParticleGenerator: Bad particle type: Li6` | npsim has no generic-ion gun; use `ion_gun_hepmc.py` (§5.3) |
| `Error when moving to event - EOF` from npsim on a `.hepmc` file | DD4hep routes `.hepmc` to HepMC3's `ReaderAscii`; write Asciiv3, not HepMC2 `IO_GenEvent` |
| `ERROR::ReaderAscii: event parsing failed … V -1 0 [] @ 0 0 0 0` | a vertex at the origin with nothing incoming must be omitted entirely, with the E-line vertex count 0 |
| everything flies to the ZDC in a gun scan | `epic_craterlake.xml` loads the 5×41 beamline fields; use the optics file matching the momenta (§5.3) |
| npsim leaves `calibrations/`, `fieldmaps/`, `gdml/` in the working directory | `cd` to a scratch directory inside the `bash -lc` string; never run it from the repository (§5.3) |
| the `--positions` output mixes the two Roman-Pot stations | both write into `ForwardRomanPotHits`; group by station plane, which `ff_gun_hits.py` has done since 2026-08-28 |
| a far-out hit appears in BOTH layers of a station | the stations are tilted 45 mrad, so a plane is not a surface of constant z; assign in the rotated frame (§5.3) |
| `pyHepMC3 rootIO.ReaderRootTree` segfaults | use the legacy `jug_xl-nightly` container for HepMC3 tree reading (§5.2) |
| `build_report.py` exits "no headless-capable browser" | install one into the user cache (§1.6) |
| a money plot is missing when building a report | run the `evgen/scripts` that produces it first (§4); the builder refuses to use a stale figure |

---

## 10 · Runtimes

Measured 2026-08-26 on eight cores (Python 3.11.4, numpy 1.25.1), single
job, warm page cache; the PYTHIA-consuming rows, the `fastsim` suite and
the rows added since were re-measured 2026-08-28, after the sample
regeneration of §5.1 raised the particle multiplicity 27–36%.  The two
coherent rows were re-measured the same day on the seven-bin |t| window,
and are quoted at their *published* `--n-mc` rather than the 6×10⁵ default
(4.0 s each at the default): 17 s for `money_cos2phi_coherent_reco` at
6×10⁶ recoils and 9 s for `nearbeam_reach_gain` at 2×10⁶, the latter
running all six configuration × pot combinations.  Everything
in `fastsim/` and `evgen/` together is about sixteen minutes, half of it
`money_delta_pdfgrid`, `nearbeam_zid_power` and the two PYTHIA rows.

| script | s | script | s |
|---|---:|---|---:|
| `fastsim/validate_inputs` | 1 | `evgen/phase_space_bins` | 3 |
| `fastsim/phase_space_map` | 6 | `evgen/closure_fom` | 16 |
| `fastsim/money_delta` | 2 | `evgen/money_cos2phi` | 1 |
| `fastsim/money_b1` | 1 | `evgen/money_cos2phi_coherent` | 3 |
| `fastsim/money_polemc` | 1 | `evgen/money_delta_extraction` | 2 |
| `fastsim/money_polemc --pdf grid` | 23 | | |
| `fastsim/tagging_acceptance` | 3 | `evgen/money_tagged_azz` | 3 |
| `tools/digitize_figure` (per curve set) | 1 | `evgen/target_mass_bound` | 2 |
| `fastsim/diag_sig2_grid` | 1 | `evgen/tagged_polarimetry_7li` | 3 |
| `fastsim/coverage_and_stat_maps` | 7 | `evgen/coherent_optics_scan` | 3 |
| `evgen/nearbeam_aperture_scan` | 7 | `evgen/nearbeam_reach_gain` | 9 |
| `evgen/tagging_optics` | 3 | `evgen/hfs_acceptance` (PYTHIA) | 19 |
| `evgen/nearbeam_sensor_budget` | 1 | `evgen/nearbeam_zid_power` | 48 |
| `evgen/nearbeam_two_hit` | 27 | `evgen/eic_beam_figures` | 2 |
| `fastsim/_check_reco_mask_invariants` | <1 | `evgen/reco_chain_figures` | 13 |
| `fastsim/money_delta_realistic` | 36 | `evgen/money_cos2phi_reco` | 4 |
| `fastsim/money_delta_pdfgrid` | **293** | `evgen/money_cos2phi_coherent_reco` | 17 |
| `fastsim/money_delta_20260715` | 5 | `evgen/hfs_resolution` (toy) | 10 |
| `fastsim/money_delta_20260720` | 17 | `evgen/hfs_resolution` (PYTHIA) | **98** |
| `fastsim/money_delta_20260721` | 40 | `evgen/money_cos2phi_reco --y-source hfs` | 20 |
| `fastsim/money_delta_20260724` | 24 | `evgen/money_cos2phi_reco --unfold-scan` | 4 |
| `fastsim/money_delta_20260725` | 21 | `evgen/money_cos2phi_reco --syst-scan` | 10 |
| | | `evgen/money_cos2phi_reco --isr --n-mc-per-cell 1600` | 14 |
| | | `… --isr --isr-seeds` (8 seeds) | 78 |
| `fastsim/money_delta_20260728` | 15 | `evgen` test suite | 44 |
| `fastsim/money_delta_20260729` | 15 | `fastsim` test suite | 12 |

The long poles are elsewhere:

| step | time |
|---|---|
| build PYTHIA 8 (§1.4) | 2.5 min |
| install the four PDF grids (§1.3) | minutes, network-bound |
| generate the 8 M-event PYTHIA production (§5.1) | 12 min, 3 jobs in parallel |
| stream 20 k BeAGLE e+d events (§5.2) | 2 min |
| one npsim far-forward scan (§5.3) | ~47 s geometry + ~1.3 s/event in `eic_xl-nightly`; the 468-event ladder is ~11 min |
| build the five reports with PDFs (§6) | 90 s when this was measured at three pages; the two added since are not timed |

---

## 11 · Provenance

- What each result means, and its caveats: `reports/` (five pages,
  numbered 0–4 in reading order: the educational primer, the cos 2φ
  projection report, the reconstruction-chain analysis note, the
  EIC/ePIC parameter reference, and the near-beam far-forward study).
- What was done when, and why a number changed: the development-run log
  in `plans/00_README.md`.
- What is still missing from the chain, ordered:
  `plans/08_simulation_chain_completion.md`.
- External dependencies with owners and default assumptions:
  `plans/04_open_questions.md`.
- The audit the current state was built against:
  `docs/code_review_2026-08-28.md` (its predecessor,
  `docs/code_review_2026-08-25.md`, is kept as the dated record).
- Per-tool detail this manual compresses: `tools/pythia8/README.md`,
  `tools/beagle/README.md`, `tools/fullsim/README.md`,
  `evgen/README.md`, `evgen/samples/README.md`, `fastsim/README.md`.
