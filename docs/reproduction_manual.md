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
door); the ePIC far-forward acceptance scans; all three reports.  That is
every published number except the four cases below.

**Not reproducible here, and honestly so.**

| what | why | what stands in for it |
|---|---|---|
| BeAGLE breakup of A = 6, 7 | BeAGLE links FLUKA, whose licence is personal and per-user (`tools/beagle/README.md`) | the cluster-IA model of `polli_fastsim/spectator.py`, whose tail is *calibrated* against the official BeAGLE e+d sample (§5.2) — and found wanting, which is the result |
| ePIC calorimeter noise/threshold floor at Σ_h ≈ 0.2–0.5 GeV | an ePIC design number nobody outside the collaboration has (plans/04 #21) | a noise scan, 0 → 25 → 50 → 100 MeV, which brackets it (§4.3) |
| the backward-disk angular resolution | same (plans/04, F3) | the repository's placeholder table, flagged as such everywhere it is used |
| radiative corrections | work package WP4, not written | nothing; the affected claims say so |

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
cd evgen   && python3 -m pytest tests/ -q     # 206 passed, ~35 s
cd fastsim && python3 -m pytest tests/ -q     # 56 passed, ~3 s
```

262 tests, all of which run without the PDF grids except the two in
`fastsim/tests/test_grids.py`, which skip.  These are not smoke tests:
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

The number the proposal asks for.  Expected, at the β = 0.30 central
short-range scale, high-acceptance / high-divergence optics:

```
6Li α-tag (137.5 GeV/u)   tagged 0.019 / 0.015
6Li d-tag                 tagged 0.091 / 0.060
7Li α-tag (117.9 GeV/u)   tagged 0.979 / 0.966
7Li t-tag                 tagged 0.006 / 0.004
```

The ⁶Li α is beam-blind (rigidity ratio 0.99813, inside the ±5% near-beam
band) and is recovered only through the p_T tail, which is why its number
is small and why §5.2 matters: the tail is the least trustworthy part of
the cluster model.

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
(published R) for LOW / MID / TOP.

The six dated notes in `fastsim/notes/` carry a banner recording this;
their numbers are still exactly reproducible with the default.

Structural invariants of the reco path are checked separately:

```bash
python3 scripts/_check_reco_mask_invariants.py     # S1–S6, exit 0
```

### 3.4 The whole dated line

`fastsim/scripts/money_delta_2026MMDD.py` are the dated production
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
python3 scripts/money_tagged_azz.py        --outdir . --events 400000
python3 scripts/tagged_polarimetry_7li.py  --outdir .
```

`closure_fom.py` is the one to run if you doubt the generator: it
compares the spread of pseudo-experiments against the analytic error
formulas over ~65 x bins.

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
python3 scripts/money_cos2phi_coherent_reco.py --outdir . # 6R
```

Three switches change what these mean.

**`--y-source hfs --hfs-sample <p.npz> <n.npz>`** replaces the 25%
Gaussian stand-in for the hadronic y by a real hadronic final state
through the hadron-side detector response.  Needs §5.1.

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

### 4.4 The hadronic-final-state resolution

```bash
python3 scripts/hfs_resolution.py --config 1 \
    --sample samples/pythia8_e10_p50_dis.npz samples/pythia8_e10_n50_dis.npz --outdir .
```

`--config 0/1/2` selects the low / mid / top beam configuration
(5 × 40.8, 10 × 99.5, 18 × 137.5 GeV/u for ⁶Li — γ-matched at the two
lower configurations, plans/10).  Without `--sample` the
toy string-fragmentation generator is used and every output is labelled
"toy"; its numbers are illustrative and, measured against PYTHIA,
optimistic (§7).

### 4.5 The coherent optics scan (WP5)

```bash
python3 scripts/coherent_optics_scan.py --outdir .
```

Tagged fraction and the one-year error on the deformation coefficient
against the near-beam envelope, for the circular, square and slot
cutouts.  Compare its answer with the *measured* ePIC aperture of §5.3 —
they disagree, and that disagreement is an open item (plans/04 #20).

### 4.6 The near-beam study (plans/09)

```bash
python3 scripts/eic_beam_figures.py       --outdir .   # 2 s
python3 scripts/nearbeam_aperture_scan.py --outdir .   # 3 s
python3 scripts/nearbeam_reach_gain.py    --outdir .   # 2 s
python3 scripts/nearbeam_sensor_budget.py --outdir .   # 1 s
python3 scripts/nearbeam_zid_power.py     --outdir .   # 48 s
```

Three questions, three scripts. The first prices *every* near-beam
aperture — coherent tagged fraction and α-tag acceptance against the
half-width in angle, per optics, with the measured ePIC aperture and the
10σ envelope marked. The second runs the full coherent chain at both
apertures and reports what the *measurement* does, not just the
acceptance. The third asks whether a superconducting nanowire can be the
thing that delivers a closer approach: energy deposits in a 12 nm NbN
film, the hot-spot firing-threshold model of charge identification
(Figure 3 of the report), the sizing, and the channel count at each
available granularity.

The coherent script also grew `--near-beam-mrad`, which replaces the
measured *horizontal* aperture and keeps the measured vertical:

```bash
python3 scripts/money_cos2phi_coherent_reco.py --config 1 \
        --rp-aperture measured --cut-scale-x 1.0 --near-beam-mrad 0.727
```

`--cut-scale-x 1.0` matters. The default 2.5 comes from the
pre-measurement belief in a wide horizontal slot, and on top of a
measured geometric aperture it imposes a 25σ horizontal retraction that
binds *before* the geometry does — hiding the whole effect.

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
    --electron-energy 10 --p-per-nucleon 50 --seed 101 --quiet \
    --out evgen/samples/pythia8_e10_p50_dis.npz
python3 tools/pythia8/gen_dis_hfs.py --target n --n-events 2000000 \
    --electron-energy 10 --p-per-nucleon 50 --seed 102 --quiet \
    --out evgen/samples/pythia8_e10_n50_dis.npz
```

9 000–13 000 events/s, 330 MB and ~115 s per million events at 10 × 99.5.
The standing production is six files (p and n at each of the three
configurations, 8 M events, 2.7 GB); `evgen/samples/README.md` is the
manifest with each file's cross section and seed, and the whole set takes
twelve minutes with three jobs in parallel on eight cores.

The files are git-ignored.  The p and n files are merged by event count,
which is what Z = N = 3 asks for in ⁶Li — but the two cross sections are
not equal (0.6661 against 0.5761 μb at 10 × 99.5), so equal event counts
are a *choice of weighting*, right for a library of hadronic shapes and
wrong for anything that reads the sample as a rate.

**Two traps, both already handled in the script, both worth knowing.**

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
E − p_z = 8.8 MeV at 50 GeV/u, so Σ − 2 E_e y = 9.4 MeV — 0.14% of Σ at
y = 1 but 2% at y = 0.01 and 7.7% at y = 0.004.  Nothing downstream is
affected, because `HFSLibrary` transfers the *ratio* of the library
event's own sums, in which the offset cancels; but do not apply the check
to a PYTHIA sample below y ≈ 0.2.  The p_T identity is exact (≤ 5×10⁻¹¹).

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

Use the **matching optics file** for the beam energy — the ⁶Li momentum
must put it at that ring's reference rigidity:
`epic_craterlake_5x41.xml` with `--p-per-nucleon 20.5`,
`epic_craterlake_10x100.xml` with 50, `epic_craterlake_18x275.xml` with
137.5.  Plain `epic_craterlake.xml` loads the 5×41 beamline fields, which
sends 275-optics momenta straight to the ZDC.

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
| ⁶Li α-tag, β = 0.30 | `scripts/tagging_acceptance.py` | 0.019 (HA) / 0.015 (HD) |
| ⁷Li α-tag, β = 0.30 | same | 0.979 / 0.966 |
| A_bag triple (frozen R) | `scripts/money_delta_20260729.py --emit-a-bag-reference` | −0.317767 / −0.310041 / −0.296750 |
| A_bag triple (published R) | `… --r-model r1998 --emit-a-bag-reference` | −0.237040 / −0.235825 / −0.234926 |
| L₅σ (frozen R) | `scripts/money_delta_realistic.py --configs low,mid,top` | 135.31 / 131.26 / 274.64 fb⁻¹/u |
| L₅σ (published R) | `… --r-model r1998 --configs low,mid,top` | 67.51 / 65.80 / 155.12 fb⁻¹/u |

### Event generator

| what | command (from `evgen/`) | expected |
|---|---|---|
| 5R sweet-spot purity, 25% stand-in | `scripts/money_cos2phi_reco.py` | 0.65 / 0.64 / 0.66 / 0.68 (D = 0.91 / 0.99 / 0.91 / 0.97) |
| 5R sweet-spot purity, PYTHIA HFS | `… --y-source hfs --hfs-sample …` | 0.43 / 0.54 / 0.47 / 0.69 |
| 5R amplitude dilution D, PYTHIA HFS | same | 0.79 / 0.85 / 0.82 / 0.95 |
| Σ-method δy/y, PYTHIA, mid config | `scripts/hfs_resolution.py --config 1 --sample …` | 0.55 / 0.28 / 0.50 / 0.15 |
| … low config | `--config 0` | 0.28 / 0.21 / 0.24 / 0.11 |
| … top config | `--config 2` | 0.74 / 0.34 / 0.69 / 0.18 |
| unfolding model dependence | `scripts/money_cos2phi_reco.py --unfold-scan` | bin-by-bin (−5.0, +8.5, −9.3, +6.3)% → folded (−0.9, −1.2, −0.2, +0.3)% |
| coherent tagged fraction | `scripts/coherent_optics_scan.py` | 32% / 3.0% / 4×10⁻⁵ / 2×10⁻⁷ at 0.10 / 0.22 / 0.45 / 0.60 GeV |
| near-beam gain, coherent | `scripts/nearbeam_aperture_scan.py` | silicon → 10σ: 1.41×10⁻² → 3.71×10⁻¹ (×26), 5.12×10⁻⁵ → 2.91×10⁻² (×569), 1.94×10⁻¹⁷ → 1.97×10⁻⁹ (dead either way) |
| near-beam gain, α tag | same | 0.103 → 0.550, 0.024 → 0.137, 0.0012 → 0.0054 |
| near-beam gain, through the chain | `scripts/nearbeam_reach_gain.py` | 5 × 41: acc 1.43×10⁻² → 3.62×10⁻¹, 2 → 3 \|t\| bins; 10 × 100: 1.02×10⁻⁴ → 3.25×10⁻², 2 → 4 bins and δa_t 1.45 → 0.0069 |
| hot-spot Z-ID thresholds | `scripts/nearbeam_sensor_budget.py` | r_s = 134 / 268 / 402 nm for p,d / α / ⁶Li; at w = 1 µm, I_th/I_c = 0.73 / 0.46 / 0.20 |
| Z-ID fake rate, 4 planes at 95% eff | `scripts/nearbeam_zid_power.py` | 2.3×10⁻⁵ (8-bit LLR) / 3.1×10⁻⁵ (one bit) / 2.7×10⁻³ (truncated mean) / 5.3×10⁻² (plain sum); 50% fill cannot reach 95% |

### Third-party

| what | expected |
|---|---|
| PYTHIA σ_gen, e+p at 10 × 99.5 | 0.6661 μb (n: 0.5761) |
| PYTHIA sample, 2 M events | 17.5 M particles, 662 MB, ~190 s |
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
job, warm page cache.  Everything in `fastsim/` and `evgen/` together is
about eleven minutes.

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
| `evgen/nearbeam_aperture_scan` | 3 | `evgen/nearbeam_reach_gain` | 2 |
| `evgen/nearbeam_sensor_budget` | 1 | `evgen/nearbeam_zid_power` | 48 |
| `fastsim/_check_reco_mask_invariants` | <1 | `evgen/reco_chain_figures` | 13 |
| `fastsim/money_delta_realistic` | 36 | `evgen/money_cos2phi_reco` | 4 |
| `fastsim/money_delta_pdfgrid` | **293** | `evgen/money_cos2phi_coherent_reco` | 3 |
| `fastsim/money_delta_20260715` | 5 | `evgen/hfs_resolution` (toy) | 10 |
| `fastsim/money_delta_20260720` | 17 | `evgen/hfs_resolution` (PYTHIA) | **64** |
| `fastsim/money_delta_20260721` | 40 | `evgen/money_cos2phi_reco --y-source hfs` | 20 |
| `fastsim/money_delta_20260724` | 24 | `evgen/money_cos2phi_reco --unfold-scan` | 4 |
| `fastsim/money_delta_20260725` | 21 | `evgen/money_cos2phi_reco --syst-scan` | 10 |
| `fastsim/money_delta_20260728` | 15 | `evgen` test suite | 30 |
| `fastsim/money_delta_20260729` | 15 | `fastsim` test suite | 3 |

The long poles are elsewhere:

| step | time |
|---|---|
| build PYTHIA 8 (§1.4) | 2.5 min |
| install the four PDF grids (§1.3) | minutes, network-bound |
| generate the 8 M-event PYTHIA production (§5.1) | 12 min, 3 jobs in parallel |
| stream 20 k BeAGLE e+d events (§5.2) | 2 min |
| one npsim far-forward scan, 84 events (§5.3) | 70 s geometry + ~1 s/event |
| build the three reports with PDFs (§6) | 90 s |

---

## 11 · Provenance

- What each result means, and its caveats: `reports/` (three pages:
  the cos 2φ projection report, the educational primer, and the
  reconstruction-chain analysis note).
- What was done when, and why a number changed: the development-run log
  in `plans/00_README.md`.
- What is still missing from the chain, ordered:
  `plans/08_simulation_chain_completion.md`.
- External dependencies with owners and default assumptions:
  `plans/04_open_questions.md`.
- The audit the current state was built against:
  `docs/code_review_2026-08-25.md`.
- Per-tool detail this manual compresses: `tools/pythia8/README.md`,
  `tools/beagle/README.md`, `tools/fullsim/README.md`,
  `evgen/README.md`, `evgen/samples/README.md`, `fastsim/README.md`.
