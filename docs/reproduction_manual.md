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
cd evgen   && python3 -m pytest tests/ -q     # 276 passed, ~46 s
cd fastsim && python3 -m pytest tests/ -q     # 81 passed, ~11 s
python3 tools/consistency_check.py --verbose  # 23 checks, whole repository
```

357 tests, all of which run without the PDF grids except three of the four
in `fastsim/tests/test_grids.py`, which skip.  These are not smoke tests:
they pin physics identities against independent constructions — the
spin-1 cross section against an explicit density-matrix trace, the
covariant azimuth against a boost-and-rotate construction of the
collinear frame, R1998 against the published fit's own worked values, the
nuclear masses against CODATA, the tensor sign against Cosyn Eq. (27).
If these pass, the machinery is sound and the rest of this manual is
about numbers, not correctness.

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
6Li α-tag   YR high-acceptance  0.017 / 0.015 / 0.016    tagging optics  0.353 / 0.273 / 0.281    legacy 73 μrad  0.167 / 0.025 / 0.018
6Li d-tag   YR high-acceptance  0.078 / 0.061 / 0.068    tagging optics  0.55 / 0.49 / 0.50
7Li α-tag   YR high-acceptance  0.967 / 0.966 / 0.971    tagging optics  0.986 / 0.990 / 0.991
7Li t-tag   YR high-acceptance  0.033 / 0.004 / 0.005
```

The ⁶Li α is beam-blind (rigidity ratio 0.99813, inside the ±5% near-beam
band): at the Yellow Report optics its near-beam tail is inside the
envelope at every configuration and what survives is the slice that falls
below R = 0.95 into the Roman-Pot window; the tagging optics (1/7–1/13 of
the luminosity) recovers the tail.  §5.2 still matters — the tail is the
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
Report optics is 1.5–1.7% on the pure spectator model and 2.5–2.9% on the
tagged one, and neither is wrong — the D wave is not an optional tail, it
*is* the tensor observable (with P_D = 0 the α–d density is m-independent
and A_zz^tag vanishes identically), so the tagged observables are quoted
on the S + D spectrum and this table on the S-wave one.  For ⁷Li the two
densities are identical (both pure P wave), and the 0.6 points by which
the table's 0.9674 at 5 × 40.8 exceeds `tagged_polarimetry_7li.py`'s
0.9614 is **two different acceptance definitions, not two different
densities**.  This table is `1 - lost`, i.e. any far-forward system; the
tagged script's mask is the Roman Pots alone (route 1 | 4), and at
5 × 40.8 the B0 carries 1.1% of the ⁷Li α.  Like for like, the pure
model's Roman-Pot tag is 0.9568 / 0.9668 / 0.9721, and restricted to
k ≤ 1.2 GeV/c where the tagged model's momentum grid ends it is
0.9626 / 0.9683 / 0.9736 against the tagged sampler's
0.9614 / 0.9676 / 0.9726 — 0.1 point at every configuration.  (The grid
truncation alone is worth +0.6 point at 5 × 40.8 and +0.2 at the other
two on the Roman-Pot mask, +0.1 uniformly on `1 - lost`; the B0 fraction
is 0.0000 at 10 × 99.5 and 18 × 117.9, where the two definitions
coincide.)  `tagged_polarimetry_7li.py` now prints both definitions.

### 3.3 The July money-Δ production, and the R that changes it

This is the line that produces the discovery-luminosity numbers.  It has
two switches that matter.

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

Expected (§7 has the table): `L_5σ` at Δ/F₁ = 10⁻³ goes from
135.3 / 131.3 / 274.6 fb⁻¹/u (frozen) to **67.5 / 65.8 / 155.1**
(published R) for LOW / MID / TOP.  **Caveat (2026-08-27):** this
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
```

`closure_fom.py` is the one to run if you doubt the generator: it
compares the spread of pseudo-experiments against the analytic error
formulas over ~65 x bins.

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

**Money plot 4 at 10 × 99.5** (`--events 400000`, seed 20260713).  The α
tag is 0.0250 at the Yellow Report optics against 0.3046 at the tagging
optics — but the tagging optics runs at L/L_HA = 1/13.3, so in tagged
events per year the tagging optics costs 8% here (acc × L = 0.0250 vs
0.0230), which is the worst of the three configurations.  What the 8%
buys is the REACH, and that is the result: the median accepted
spectator momentum is 0.322 GeV/c at the Yellow Report optics with
**nothing at all below k = 0.15 GeV/c**, against 0.162 GeV/c with 44%
below 0.15 at the tagging optics (0.348/0.322/0.333 and
0.145/0.162/0.161 GeV/c at 5 × 41 / 10 × 100 / 18 × 275; the unfolded
sample has median 0.093 and 74% below 0.15).  The rate trade is even or
better at the other two: acc × L/L_HA is 0.0285 → 0.0540 at 5 × 41,
0.0250 → 0.0230 at 10 × 100 and 0.0267 → 0.0330 at 18 × 275, so at two
of the three configurations the tagging optics gains the reach and the
rate together, and pays for it only at 10 × 100.  At every published optics
the ⁶Li α tag admits only the high-k tail; the tagging optics is what
turns money plot 4 from a one-point measurement into a curve.

The right panel's coloured curves are the truth weighted by each optics'
own θ_k acceptance (`tagged.acceptance_weights`,
`tagged.azz_tensor_curve`), and that is what the markers measure: at
k ≈ 0.325 GeV/c the folded A_zz is +0.428 against a weighted truth of
+0.455 at the Yellow Report optics and −0.083 against −0.098 at the
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
that is the deliverable: the Roman-Pot tag is 0.9676 at the Yellow Report
optics against 0.9916 at the tagging optics, the folded ⟨P₂⟩ slope
−0.1947 and −0.1964 against the analytic −0.2000 (the retired 73 μrad
gave −0.1929), and the median δA_∥ 0.01150 against 0.01138 *at equal
generated statistics*.  At equal luminosity the
tagging optics multiplies every ⁷Li error bar by 2.83 / 3.87 / 3.15 at
5 × 41 / 10 × 100 / 18 × 275: ×1.02 in acceptance for ×1/8–1/15 in
luminosity is a factor 8–15 net loss, the exact inverse of ⁶Li.  ⁶Li and
⁷Li want different machine optics and are different runs.

### 4.2 The cos 2φ money plots, truth level

```bash
python3 scripts/money_cos2phi.py            --outdir .   # money plot 5
python3 scripts/money_cos2phi_coherent.py   --outdir .   # money plot 6
python3 scripts/money_delta_extraction.py   --outdir .   # money plot 7
```

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
on `--fit` (both since 2026-08-28).  `--ensemble 20` repeats the one-year pseudo-experiment
and prints the bias test of Table 5; `--exact` switches the Poisson draw
off for the systematic scans; `--config 1/2` gives the other two
configurations.

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
`0.05,0.08,0.12,0.17,0.25` by any increasing list, which is how the
window below 0.05 GeV² that the published binning discards is priced
(plans/08 §8.4).  The summary ends with the quadrature combination of
δa_e over the bins that survived — a_e is one constant across |t|, so its
errors combine, and a_t is not.

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
```

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
window on both members of the pair and brings it to ≤ 0.25% while keeping
87% of the non-radiative rate; the chain does not apply that cut, so it
is a mitigation and not the bound.

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

### 4.4b Where the hadronic E − p_z sum goes (Report 2 §3, Figure 4)

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
ε = 0.422 / 0.322 / 0.332 and L/L_HA = 1/7.1 / 1/13.3 / 1/9.5,
N_tag/yr = 2.6×10⁶ / 3.0×10⁶ / 6.1×10⁶ (2.3 / 5.6 / 3.9× below the 0.20 GeV
reference at 10 fb⁻¹/u), best-super-bin 5σ floors of 1.7 / 2.1 / 1.6% per
unit P_zz, the shape term in the optics' own window 0.031–0.035 per unit
P_zz = 9.3 / 8.3 / 10.7σ per year, 2.8 / 4.4 / 2.6 years to 5σ on a 1%
exotic-glue term, and IR-8's ≈ 20% worth 3.3 / 8.2 / 5.7× the optimum at
equal luminosity; both planes de-squeezed gives a fifth of the yield at
1/24–1/70; with the pots fixed nothing is recovered at any β*.

### 4.6 The near-beam study (plans/09)

```bash
python3 scripts/eic_beam_figures.py       --outdir .   # 2 s
python3 scripts/nearbeam_aperture_scan.py --outdir .   # 7 s
python3 scripts/nearbeam_aperture_scan.py --isotope 7Li --outdir .   # 7 s
python3 scripts/nearbeam_reach_gain.py    --outdir .   # 2 s
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
*measurement* does. The third asks whether a superconducting nanowire can
be the thing that delivers a closer approach: energy deposits in a 12 nm
NbN film, the hot-spot firing-threshold model of charge identification
(Figure 3 of the report), the sizing strip at the tagging envelope, and
the channel count at each available granularity. The fourth (plans/09 B4,
2026-08-28) samples BOTH fragments of a ⁶Li → α + d breakup from one
relative momentum (`spectator.breakup_lab_kinematics`) and asks what hit
multiplicity is worth against the coherent tag: the 2 × 2 topology per
optics, the separation at the pot plane in millimetres and in 500 μm
pixels, and — the result — how often the partner deuteron is there to veto
an α that faked an intact ⁶Li, scanned against the pot's outer edge. Its
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
python3 scripts/money_cos2phi_coherent_reco.py --config 1 --optics tagging \
        --rp-aperture measured --cut-scale-x 1.0 --near-beam-mrad 0.17
```

0.17 mrad is the tagging envelope at 10 × 100, and this run returns all
four |t| bins (acc 0.315, N_tag 2.93 × 10⁶/yr).  Asking for the Yellow
Report envelope instead — `--config 1 --rp-aperture measured
--cut-scale-x 1.0 --near-beam-mrad 1.8`, at 1.8 mrad — returns nothing:
all four bins are dropped for zero accepted recoils and the script exits
non-zero.  That is not a defect of the switch but the result of §2 of
Report 4 restated at the reconstructed level — at the Yellow Report
optics the coherent channel has no acceptance at all — and it is why the
published 6R is `--optics tagging` (§4.3).  `--cut-scale-x 1.0` matters
in either case: the default 2.5 comes from the pre-measurement belief in
a wide horizontal slot, and on top of a measured geometric aperture it
imposes a 25σ horizontal retraction that binds *before* the geometry
does — hiding the whole effect.

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

```bash
python3 tools/fullsim/ion_gun_hepmc.py --out li6.hepmc --pdg 1000030060 \
    --a 6 --p-per-nucleon 137.5 \
    --pt 0.78 0.87 1.24 1.65 2.06 2.48 3.30 --nphi 12 > li6.index
    # p_T = theta * p, so at p = 6 x 137.5 = 825 GeV this is
    # theta = 0.95 / 1.05 / 1.5 / 2.0 / 2.5 / 3.0 / 4.0 mrad; 7 x 12 = 84 events

singularity exec $SIF bash -lc '
  source /opt/detector/epic-main/bin/thisepic.sh
  npsim --compactFile $DETECTOR_PATH/epic_craterlake_18x275.xml \
        --inputFiles li6.hepmc --numberOfEvents 84 --physics.list FTFP_BERT \
        --part.minimalKineticEnergy "100*MeV" --outputFile gun_li6.edm4hep.root'
                                              # ~70 s geometry + ~1 s/event

singularity exec $SIF python3 tools/fullsim/ff_gun_hits.py gun_li6.edm4hep.root \
    --per-event --index li6.index --positions
```

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
is where the lithium fill actually sits (0.49 / 0.81 / 0.85 GeV rather
than 0.25 / 0.41 / 0.85).

Expected: no Roman-Pot hits below |θ_x| ≈ 1.03 mrad at 18 × 275, hits
from 1.09 mrad at φ = 0 and 180°, and nothing from a *vertical* kick
until ≈ 2 mrad.  An α at the same rigidity (550 GeV) gives the same edge
at the same *angle*, which is the check that the scan measures optics and
not species.  `tools/fullsim/README.md` has the full azimuthal map, the
optics levers read off the hit positions, and what it means for the
coherent channel.

The single-particle routing scan for the other fragments is
`tools/fullsim/ff_gun_scan.sh` (α, triton, deuteron, proton, neutron at
their rigidity ratios).

**The geometry version matters.**  The results above are the September
2024 `epic-main` inside `jug_xl-nightly`; that file already moved the
Roman-Pot stations from the published 26/28 m to 32.5/34.3 m.  Re-run in
a current container before quoting an acceptance in a paper.

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
| ⁶Li α-tag, β = 0.30 | `scripts/tagging_acceptance.py` | YR high-acceptance 0.017 / 0.015 / 0.016, tagging optics 0.353 / 0.273 / 0.281, legacy 73 μrad 0.167 / 0.025 / 0.018 (5 × 41 / 10 × 100 / 18 × 275) |
| ⁷Li α-tag, β = 0.30 | same | 0.967 / 0.966 / 0.971 (YR HA), 0.986 / 0.990 / 0.991 (tagging) |
| … on the tagged generator instead | `evgen/scripts/tagged_polarimetry_7li.py` | Roman-Pot tag 0.9614 / 0.9676 / 0.9726 (YR HA), 0.9807 / 0.9916 / 0.9920 (tagging); the script also prints `acc(any far-fwd)` = 0.9683 / 0.9676 / 0.9726 and 0.9876 / 0.9916 / 0.9920, which is the definition the row above tabulates.  The two ⁷Li densities are identical (both pure P wave): the 0.6-point gap at 5 × 41 is the B0, which takes 1.1% of the α there and enters `1 - lost` but not the Roman-Pot mask.  Like for like and inside the tagged model's k ≤ 1.2 GeV/c grid the pure model gives 0.9626 / 0.9683 / 0.9736 — 0.1 point.  For ⁶Li the two differ by the D wave — §4.1 |
| A_bag triple (frozen R) | `scripts/money_delta_20260729.py --emit-a-bag-reference` | −0.317767 / −0.310041 / −0.296750 |
| A_bag triple (published R) | `… --r-model r1998 --emit-a-bag-reference` | −0.237040 / −0.235825 / −0.234926 |
| L₅σ (frozen R) | `scripts/money_delta_realistic.py --configs low,mid,top` | 135.31 / 131.26 / 274.64 fb⁻¹/u (script-internal pre-2026-08-27 configs at 27.5 / 50 / 137.5 GeV/u, superseded by plans/10; only TOP is a machine configuration) |
| L₅σ (published R) | `… --r-model r1998 --configs low,mid,top` | 67.51 / 65.80 / 155.12 fb⁻¹/u (same caveat) |
| L₅σ toy, current energies | `scripts/money_delta.py` | 16.7 / 16.3 / 21.8 fb⁻¹/u at Δ/F₁ = 10⁻³, P_zz = 0.8 (16.719 / 16.332 / 21.811 to the precision the 2026-08-28 min-events correction is pinned at, `tests/test_money_delta_mask.py`) |

### Event generator

| what | command (from `evgen/`) | expected |
|---|---|---|
| 5R sweet spots (x, Q²) | `scripts/money_cos2phi.py` | (0.028, 1.14), (0.011, 1.14), (0.071, 3.13), (0.141, 14.3); A = 7.4 / 4.4 / 9.5 / 9.5 ×10⁻³, δA = 1.7 / 1.4 / 2.8 / 4.5 ×10⁻⁴ (1 yr) |
| 5R sweet-spot purity, 25% stand-in | `scripts/money_cos2phi_reco.py` | 0.66 / 0.63 / 0.69 / 0.69 (D = 0.92 / 0.99 / 0.90 / 0.96); δÂ = 1.2 / 0.9 / 1.6 / 3.0 ×10⁻⁴ |
| 5R sweet-spot purity, PYTHIA HFS, uncalibrated | `… --y-source hfs --hfs-sample …` | 0.42 / 0.53 / 0.49 / 0.68 |
| 5R amplitude dilution D, PYTHIA HFS, uncalibrated | same | 0.79 / 0.84 / 0.83 / 0.95 |
| 5R with the hadronic scale calibrated in reconstructed bins | `… --y-source hfs --hfs-sample … --hfs-calibrate` | purity 0.56 / 0.59 / 0.64 / 0.75, D = 1.00 / 1.06 / 0.95 / 0.98, δÂ = 1.3 / 0.9 / 1.7 / 3.0 ×10⁻⁴ |
| 7R folded-fit bars at the best bins, PYTHIA calibrated | `… --hfs-calibrate --unfold folded` | 4.2 / 2.6 / 7.0% (1 yr), 3.8 / 2.2 / 2.6% (10 yr) at Q² = 1.14 / 3.13 / 14.3; prior spread 3.7 / 2.1 / 1.3% |
| residual hadronic scale +2% (calibrated) | `… --hfs-calibrate --hfs-scale 1.02` | Â moves +0.2 to +1.2% |
| Σ at the mid sweet spots | `scripts/hfs_acceptance.py --config 1 --sample …` | within acceptance and above threshold 0.80 / 0.87 / 0.83 / 0.92 at |η| ≤ 3.7 (forward loss 0.17 / 0.08 / 0.15 / 0.07), 0.85 / 0.89 / 0.88 / 0.94 at 4.0; captured through the response 0.70 / 0.74 / 0.74 / 0.85; δy/y unchanged between reaches |
| Σ-method δy/y, PYTHIA, mid config (its own sweet spots) | `scripts/hfs_resolution.py --config 1 --sample …` | 0.32 / 0.22 / 0.29 / 0.15 |
| … low config | `--config 0` | 0.39 / 0.23 / 0.24 / 0.11 |
| … top config | `--config 2` | 0.23 / 0.19 / 0.21 / 0.18 |
| unfolding model dependence (moment_B prior) | `scripts/money_cos2phi_reco.py --unfold-scan` | bin-by-bin (−4.1, +8.0, −5.5, +4.9)% → folded (−1.5, −1.8, −0.3, +0.5)% |
| collinear-ISR migration bound (plans/07 WP4) | `scripts/money_cos2phi_reco.py --isr --n-mc-per-cell 1600` | one draw at the default seed 20260824: Δ̂ biased by +0.616 / +0.501 / +0.804 / +1.240% at the four sweet spots; ⟨z⟩ = 0.0245, 26.6% of the rate radiates above z = 10⁻⁴, ⟨z\|z>10⁻⁴⟩ = 0.0923, t = 0.0700 at ⟨Q²⟩ = 4.37 GeV²; the covariant azimuth's residual under k → (1−z)k is 2.6×10⁻² rad and fakes cos 2φ′ at 9×10⁻⁸ |
| … the PUBLISHED bound, averaged over eight response seeds (`$S` below is that list) | `… --isr --isr-seeds $S --n-mc-per-cell 1600`, `S=20260824,20260925,20261026,20261127,20261228,20270129,20270302,20270403` | +0.62 ± 0.03 / +0.50 ± 0.02 / +0.94 ± 0.03 / +1.22 ± 0.02%; purity 0.653 → 0.638, 0.633 → 0.613, 0.679 → 0.659, 0.684 → 0.640; efficiency 0.414 → 0.404, 0.590 → 0.572, 0.374 → 0.369, 0.653 → 0.634 |
| … with the low-Q² feed-in opened up | `… --isr --isr-seeds $S --isr-gen-q2min 0.05` | +2.26 ± 0.03 / +2.24 ± 0.10 / +1.88 ± 0.07 / +1.88 ± 0.05%.  Worst spot over the window sequence 0.7 → 0.35 → 0.15 → 0.05 → 0.02 GeV²: 1.22 → 1.87 → 2.81 → 2.26 → 2.34%, so the bound saturates in a 1.8–2.8% band and **≤ 2.9%** is what the ≤5% gate is read against |
| … behind a HERA-style E − p_z window | `… --isr --isr-seeds $S --isr-empz 0.85` | +0.23 ± 0.02 / +0.16 ± 0.02 / +0.22 ± 0.02 / +0.18 ± 0.01%, independent of the generator window; the window keeps 86.9% of the non-radiative selected rate and 82.4% of the radiative one |
| … through the PYTHIA hadronic final state | `… --isr --isr-seeds $S --y-source hfs --hfs-sample … --hfs-calibrate --n-mc-per-cell 800` | +0.38 ± 0.03 / +0.44 ± 0.02 / +0.46 ± 0.03 / +0.88 ± 0.04% — the bound does not come from the Gaussian y stand-in |
| coherent tagged fraction | `scripts/coherent_optics_scan.py` | 32% / 3.0% / 4×10⁻⁵ / 2×10⁻⁷ at 0.10 / 0.22 / 0.45 / 0.60 GeV.  The 0.60 GeV point is response-MC-limited (18–27 tagged/yr): at `--n-mc 150000` its harmonic design loses rank, and the point is now reported and dropped rather than aborting the figure |
| coherent tag at the three half-widths | `scripts/nearbeam_aperture_scan.py` | silicon / YR HA envelope / tagging envelope: 9.8×10⁻⁷ / 7.2×10⁻⁸ / 0.42 (5 × 41), 7.7×10⁻¹⁶ / 6.2×10⁻²⁷ / 0.32 (10 × 100), 1.9×10⁻¹⁷ / 3.9×10⁻¹⁴ / 0.33 (18 × 275) |
| ⁶Li α tag (routed) at the three half-widths | same | 0.018 / 0.017 / 0.35, 0.015 / 0.015 / 0.27, 0.015 / 0.015 / 0.28 |
| ⁷Li α tag at the same three | `… --isotope 7Li` | 0.969 / 0.968 / 0.987 (5 × 41), 0.967 / 0.965 / 0.990 (10 × 100), 0.968 / 0.968 / 0.991 (18 × 275) — flat to three points across the whole 0.05–3 mrad axis, at 1/8.2 / 1/15.3 / 1/10.1 of the luminosity for the tagging point.  Panel (a) stays the coherent intact ⁶Li at either setting: `polligen.coherent` is ⁶Li-specific and a ⁷Li coherent channel is a different amplitude (plans/09 B3, open) |
| money plot 4: ⁶Li α-tag reach, 10 × 99.5 | `scripts/money_tagged_azz.py --events 400000` | acc 0.0250 (YR HA) vs 0.3046 (tagging, L/L_HA = 1/13.3) — acc × L 0.0250 vs 0.0230, an 8% cost — but median accepted k 0.322 vs 0.162 GeV/c and frac(k < 0.15) 0.000 vs 0.438.  At k ≈ 0.325 GeV/c A_zz = +0.428 (acceptance-weighted truth +0.455) and −0.083 (−0.098); the θ_k = 90° curve says −0.482 at both.  `acc` is the unbinned accepted fraction: 9–11% of the accepted α lie above the 0.6 GeV/c right edge of the panel |
| … the same at the other two configurations | `… --config 0` / `--config 2` | median accepted k 0.348 / 0.333 GeV/c (YR HA, frac below 0.15 GeV/c = 0.000) against 0.145 / 0.161 with 0.522 / 0.446 (tagging); acc × L/L_HA 0.0285 → 0.0540 and 0.0267 → 0.0330, i.e. the tagging optics gains reach *and* rate at both, and costs rate only at 10 × 100 |
| ⁷Li polarimetry and tagged EMC, 10 × 99.5 | `scripts/tagged_polarimetry_7li.py` | acc(RP) 0.9676 (YR HA) vs 0.9916 (tagging); ⟨P₂⟩ slope −0.1947 vs −0.1964 against the analytic −0.2000 (legacy 73 μrad: −0.1929); median δA_∥ 0.01150 vs 0.01138 *at equal generated statistics* — the plotted bars are drawn there, and the figure says so.  At equal luminosity the tagging optics multiplies every ⁷Li error bar by 2.83 / 3.87 / 3.15 — a factor 8–15 net loss |
| the chain at the tagging optics | `scripts/nearbeam_reach_gain.py --n-mc 2000000` | pots at the silicon: acc 0, 0 bins at every configuration; pots following: acc 0.41 / 0.31 / 0.32, N_tag 2.5 / 2.9 / 6.0 ×10⁶/yr, 4 of 4 \|t\| bins, δa_t 0.0035–0.018 / 0.0034–0.0135 / 0.0023–0.0151.  `--fit likelihood` leaves the errors alone (0.0035–0.0187 / 0.0034–0.0137 / 0.0023–0.0152) and moves only the recovered a_t of the sparsest bin, 0.118 → 0.197 against 0.181 injected at 5 × 41 |
| sizing strip at the tagging envelope | `scripts/nearbeam_sensor_budget.py` | d50 / d90 / d99 = 183 / 504 / 842 μrad (5 × 41), 69 / 194 / 328 (10 × 100), 50 / 141 / 239 (18 × 275); α at 137.5 GeV/u 77 / 269 / 624 |
| tagging optics, priced | `scripts/tagging_optics.py` | horizontal-only optimum β*_x/β*_x,HA = 49.7 / 175.6 / 89.3, ε = 0.422 / 0.322 / 0.332, L/L_HA = 1/7.1 / 1/13.3 / 1/9.5, N_tag/yr = 2.6×10⁶ / 3.0×10⁶ / 6.1×10⁶, 5σ floor/yr = 1.7 / 2.1 / 1.6% per unit P_zz, shape term 9.3 / 8.3 / 10.7σ/yr |
| 6R at the tagging optics, 5 × 40.8 | `scripts/money_cos2phi_coherent_reco.py --config 0 --optics tagging --n-mc 6000000 [--ensemble 20]` | σ_θ = 33/380 μrad, cutout 0.33 × 3.8 mrad, acc 0.411, N_tag 2.52×10⁶/yr at L/L_HA = 1/7.1, ⟨cos 2β⟩ = −0.27; a_t 0.0899 ± 0.0035 / 0.118 ± 0.005 / 0.131 ± 0.009 / 0.104 ± 0.018 (1 yr; inj. 0.090 / 0.118 / 0.147 / 0.180), 0.0892 / 0.1202 / 0.1465 / 0.183 at 10 yr; ensemble means 0.0899 / 0.1170 / 0.139 / 0.114 |
| … the same with the likelihood estimator | `… --ensemble 20 --fit likelihood` | ensemble means 0.0902 / 0.1176 / 0.146 / 0.173 (σ 0.0033 / 0.0038 / 0.0085 / 0.0207) on the same twenty draws.  With `--ensemble 200` (8 min): likelihood 0.0902 / 0.1185 / 0.1477 / 0.1816, pulls of the mean +2.0 / +1.0 / +0.5 / +0.9, spreads 0.0031 / 0.0048 / 0.0091 / 0.0193 against quoted errors 0.0035 / 0.0048 / 0.0086 / 0.0187; ratio 0.0899 / 0.1179 / 0.1414 / 0.1184, pulls +0.4 / −0.7 / −9.0 / −42.1 |
| … coarser (α, β) bins, same 200 draws | `… --ensemble 200 --n-alpha 8 --n-beta 16` (and `--n-alpha 6 --n-beta 12`) | the sparsest bin's ratio mean rises 0.1184 → 0.1552 → 0.1662 against 0.1803 injected, i.e. −34.3% → −13.9% → −7.8%, at δa_e 0.0139 → 0.0148 → 0.0162 (+17%).  Coarser binning attenuates the bias, it does not remove it |
| 6R at the tagging optics, 18 × 137.5 | `… --config 2 --optics tagging --n-mc 6000000` | acc 0.324, N_tag 5.99×10⁶/yr; a_t 0.1005 ± 0.0023 / 0.133 ± 0.003 / 0.177 ± 0.006 / 0.221 ± 0.015 (inj. 0.100 / 0.137 / 0.179 / 0.228); ensemble means 0.0998 / 0.1368 / 0.1776 / 0.207 (ratio) and 0.0998 / 0.1372 / 0.1788 / 0.230 (`--fit likelihood`, same draws) |
| 6R at the tagging optics, 10 × 99.5 | `… --config 1 --optics tagging --n-mc 6000000 --ensemble 20 --fit likelihood` | acc 0.315, N_tag 2.93×10⁶/yr at L/L_HA = 1/13.3; ensemble means 0.0896 / 0.1131 / 0.1369 / 0.1576 against inj. 0.0897 / 0.1136 / 0.1360 / 0.1595 |
| in-situ (u₁, u₂), tagging optics | `… --config 0 --optics tagging --exact --u-in-situ --n-mc 6000000` | δu₂ = 0.0028 / 0.0040 / 0.0070 / 0.0136 per \|t\| bin at 5 × 40.8 in one year (0.0016 / 0.0026 / 0.0049 / 0.0111 at 18 × 137.5, `--config 2`), i.e. 1.8–15× the ZEUS 1σ of 0.024; the propagated a_e term is 0.00005–0.00063 across the eight bins, negligible against the 0.0016–0.0141 statistical error |
| \|t\| re-binned below the window | `… --config 0 --optics tagging --exact --n-mc 6000000 --t-edges 0.006,0.017,0.028,0.039,0.05,0.08,0.12,0.17,0.25` | δa_e = 0.0022 / 0.0022 / 0.0026 / 0.0031 in the four added bins, and the printed combination falls from 0.00205 (the four published bins) to 0.00105 over the eight |
| 6R systematics at 5 × 40.8, exact counts | `… --exact --no-sin --envelope-split 1e-3 / --u2-assumed 0.044 / --rel-lumi-offset 1e-3` | a_t −0.8 / −0.5 / −0.1 / −0.3%; a_e +0.0007 to +0.0009; a_t 0.1% |
| α + d separation at the pots | `scripts/nearbeam_two_hit.py` | median 38.4 / 18.5 / 15.1 mm (16–84% 17.1–79.4 / 8.2–38.7 / 6.6–32.1) = 77 / 37 / 30 pixels of 500 μm at 5 × 41 / 10 × 100 / 18 × 275; the angular lever alone gives 36.7 / 15.1 / 10.9, the rest is the pot dispersion D = 0.30 m acting on rigidities that move apart with k_z; θ_d/θ_α → 1.987 as k → 0.  `--beta 0.20 / 0.40` moves the medians to 33.6 / 16.2 / 13.2 and 41.7 / 20.1 / 16.4 mm (−12% to +9%) and the veto by less than 0.03 |
| α + d topology per breakup | same | both fragments recorded 0.0002 / 0.0000 / 0.0003 (YR high acceptance) and 0.287 / 0.216 / 0.223 (tagging); d alone 0.075 / 0.062 / 0.072 and α alone 0.017 / 0.015 / 0.016 (YR); a recorded pair lands inside one 500 μm pixel in 2 × 10⁻⁶ / 1.7 × 10⁻⁴ / 4.1 × 10⁻⁴ of cases at the tagging optics and never at the YR or legacy envelopes |
| partner-fragment veto | same | P(α fakes a coherent tag) 0.0020 / 0.0001 / 0.0007 (YR) and 0.338 / 0.258 / 0.265 (tagging); of those the partner d is recorded 0.120 ± 0.002 / 0.017 ± 0.004 / 0.245 ± 0.005 (YR, on 2.4 × 10⁴ / 1.4 × 10³ / 8.7 × 10³ fakes) and 0.850 / 0.838 / 0.840 (tagging), falling to 0.00 / 0.31 / 0.56 against a 0.5 mrad outer edge |
| hot-spot Z-ID thresholds | `scripts/nearbeam_sensor_budget.py` | r_s = 134 / 268 / 402 nm for p,d / α / ⁶Li; at w = 1 µm, I_th/I_c = 0.73 / 0.46 / 0.20 |
| Z-ID fake rate, 4 planes at 95% eff | `scripts/nearbeam_zid_power.py` | 2.3×10⁻⁵ (8-bit LLR) / 3.1×10⁻⁵ (one bit) / 2.7×10⁻³ (truncated mean) / 5.3×10⁻² (plain sum); 50% fill cannot reach 95% |

### Third-party

| what | expected |
|---|---|
| PYTHIA σ_gen, e+p at 10 × 99.5 | 0.9473 μb (n: 0.8551), with mHatMin = 0.5 |
| PYTHIA sample, 2 M events | 25.3 M particles, 910 MB, ~270 s |
| BeAGLE e+d, P(p_T > 0.3) in the x_L peak | 0.0261, against the Hulthén model's 0.0037 |
| ⁶Li Roman-Pot edge, 18 × 275 optics | between 0.970 and 1.091 mrad, horizontal |

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
| `pyHepMC3 rootIO.ReaderRootTree` segfaults | use the legacy `jug_xl-nightly` container for HepMC3 tree reading (§5.2) |
| `build_report.py` exits "no headless-capable browser" | install one into the user cache (§1.6) |
| a money plot is missing when building a report | run the `evgen/scripts` that produces it first (§4); the builder refuses to use a stale figure |

---

## 10 · Runtimes

Measured 2026-08-26 on eight cores (Python 3.11.4, numpy 1.25.1), single
job, warm page cache; the PYTHIA-consuming rows, the `fastsim` suite and
the rows added since were re-measured 2026-08-28, after the sample
regeneration of §5.1 raised the particle multiplicity 27–36%.  Everything
in `fastsim/` and `evgen/` together is about sixteen minutes, half of it
`money_delta_pdfgrid`, `nearbeam_zid_power` and the two PYTHIA rows.

| script | s | script | s |
|---|---:|---|---:|
| `fastsim/validate_inputs` | 1 | `evgen/phase_space_bins` | 3 |
| `fastsim/phase_space_map` | 6 | `evgen/closure_fom` | 16 |
| `fastsim/money_delta` | 2 | `evgen/money_cos2phi` | 1 |
| `fastsim/money_b1` | 1 | `evgen/money_cos2phi_coherent` | 3 |
| `fastsim/money_polemc` | 1 | `evgen/money_delta_extraction` | 2 |
| `fastsim/tagging_acceptance` | 3 | `evgen/money_tagged_azz` | 3 |
| `fastsim/diag_sig2_grid` | 1 | `evgen/tagged_polarimetry_7li` | 3 |
| `fastsim/coverage_and_stat_maps` | 7 | `evgen/coherent_optics_scan` | 3 |
| `evgen/nearbeam_aperture_scan` | 7 | `evgen/nearbeam_reach_gain` | 2 |
| `evgen/tagging_optics` | 3 | `evgen/hfs_acceptance` (PYTHIA) | 19 |
| `evgen/nearbeam_sensor_budget` | 1 | `evgen/nearbeam_zid_power` | 48 |
| `evgen/nearbeam_two_hit` | 27 | `evgen/eic_beam_figures` | 2 |
| `fastsim/_check_reco_mask_invariants` | <1 | `evgen/reco_chain_figures` | 13 |
| `fastsim/money_delta_realistic` | 36 | `evgen/money_cos2phi_reco` | 4 |
| `fastsim/money_delta_pdfgrid` | **293** | `evgen/money_cos2phi_coherent_reco` | 3 |
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
| one npsim far-forward scan, 84 events (§5.3) | 70 s geometry + ~1 s/event |
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
