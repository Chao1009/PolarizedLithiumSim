# Running BeAGLE locally — status & instructions

*Reproducing the e+d control rather than building BeAGLE?  It needs no
licence and no local build:
[docs/reproduction_manual.md](../../docs/reproduction_manual.md) §5.2.*

Goal: build `../../../BeAGLE` (eic/BeAGLE, master) on this machine for the
e+⁶Li/⁷Li breakup & tagging study (plans/02 step 1.5).

## Where this stands on the analysis box (checked 2026-08-26, container line updated 2026-08-28)

**Nothing of the local build is present here**: no `~/Projects/eic/beagle_deps`,
no BeAGLE checkout, no FLUKA, and under `~/Projects/eic/local/lib` only the
legacy `jug_xl-nightly.sif` container.  That statement is about
`~/Projects/eic/local/lib` alone: an `eic_xl-nightly.sif` image (4.5 GB)
arrived on 2026-08-28 under `~/Projects/eic-2026/local/lib`, and it is the
container the far-forward aperture measurement was made in
(`tools/fullsim/README.md`, plans/09 B1).  The dependency table below is the
recipe, not a description of this machine.

**What is verified working here, and is the reason the local build is not
on the critical path**, is the no-FLUKA route: the container ships `xrdcp`,
`xrdfs` and a working `pyHepMC3.rootIO.ReaderRootTree`, and
`root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/DIS/` is reachable and
lists `BeAGLE1.03.02-{1.0,1.2,1.3,2.0,2.1,3.0,3.1}` with `eH2` (deuteron),
`eHe3` and `eAu` — the e+d control of plans/02 step 1.5.3 in full, no
licence required.  It has now been run; the result is below.

## The e+d control, run 2026-08-26 — the cluster model's tail is too soft

20 000 events of
`BeAGLE1.03.02-3.1/eH2/en/9x130/q2_1to1000/..._ab_run001.hepmc3.tree.root`
streamed through `../analysis/dump_spectators.py` and analysed with
`../analysis/ed_control_analysis.py --beta-scan`.  Spectator protons are
found in 100.0% of events and route to the OMD (93.1%, R ≈ 0.5) as the
rigidity argument says they must — the routing logic agrees with BeAGLE to
better than 2 points in every window.

The **p_T tail does not**.  Restricted to the spectator peak
x_L ∈ [0.9, 1.1) so that target fragmentation cannot contribute:

| P(p_T >) | 0.10 | 0.20 | 0.30 | 0.45 GeV |
|---|---|---|---|---|
| BeAGLE e+d | 0.1882 | 0.0494 | 0.0261 | 0.0144 |
| Hulthén, β = 0.26 (the script's default) | 0.1522 | 0.0201 | 0.0037 | 0.0005 |
| Hulthén, β = 0.40 (best fit) | 0.1934 | 0.0360 | 0.0089 | 0.0015 |
| Hulthén, β = 1.50 | 0.2505 | 0.0737 | 0.0303 | 0.0106 |

No β reproduces the shape: BeAGLE has a *narrower core* and a *harder
tail* than the two-parameter Hulthén can make simultaneously, and the best
shape fit, β = 0.40 GeV, sits at the top of the 0.20–0.40 scan range
`spectator.py` documents while still falling a factor 3 short at 0.30 GeV
and 10 at 0.45 GeV.  This matters because the ⁶Li α tag is *entirely* a
p_T-tail measurement — the angular near-beam envelope is 0.40 GeV of p_T
for a 550 GeV α under the high-acceptance optics — so the published
α-tag acceptance is model-limited **from below**, and its model
uncertainty is one-sided upward.  Caveats before this is quoted as a
correction rather than a warning: BeAGLE's deuteron is not its A = 6
treatment (A > 4 uses the C-12 Fermi parameterization with no α+d
clustering), and the external anchor plans/02 step 1.5.3 actually asks for
is Tu et al., not BeAGLE itself.

## Dependency recipe (built 2026-06-12 into `~/Projects/eic/beagle_deps/install`)

| dependency | status | notes |
|---|---|---|
| gfortran 11.4 | ✅ system | legacy code needs `-std=legacy` (handled by build script) |
| LHAPDF 5.9.1 | ✅ built | includes `eps09.f`; grids installed: cteq6ll.LHpdf (=CTEQ6L1, set 10042), EPS09LOR_{4,6,12}.LHgrid — **A=6 is a valid EPS09 nucleus**, and BeAGLE's `anear.f` maps both ⁶Li and ⁷Li → A=6 |
| RAPGAP 3.302 libs | ✅ built | `librapgap33.a` (sfecfe.o removed), `libar4.a`, `libbases.a` copied into `BeAGLE/RAPGAP-3.302/lib/` |
| PYTHIA 6.4.28 | ✅ built | standalone `libpythia6.a` (for RAPGAP configure); BeAGLE builds its own internally |
| CERNLIB 2024 (free) | ✅ built | mathlib/kernlib/packlib static, graphics/PAW tree patched out (no Motif on this box) |
| **FLUKA** | ❌ **needs your action** | license registration is personal — see below |
| nuclear.bin | ✅ in BeAGLE repo | FLUKA evaporation data file (link into run dir) |

## The one missing piece: FLUKA (~10 minutes of your time)

BeAGLE links `libflukahp.a` and compiles against `$FLUPRO/flukapro/`
includes — both come only with the **INFN FLUKA** distribution:

1. Register (free for research) at **https://www.fluka.org** → Download.
2. Download the 64-bit binary release matching **gfortran 11**
   (releases are tagged by gfortran major version, e.g.
   `fluka2024.x-linux-gfor64bit-11_amd64.tar.gz` — pick the gfortran-11
   variant; BeAGLE was developed against the older `fluka2011.2x` line,
   so if the final link complains we install that instead).
3. Unpack and set:
   ```bash
   mkdir -p ~/Projects/eic/beagle_deps/fluka && cd ~/Projects/eic/beagle_deps/fluka
   tar xzf ~/Downloads/fluka20XX...tar.gz
   export FLUPRO=$PWD            # dir containing libflukahp.a and flukapro/
   ```
4. Then: `cd ~/Projects/eic/pol_li/tools/beagle && ./build_beagle.sh`

## Build & run

```bash
./build_beagle.sh            # needs $FLUPRO; builds ./BeAGLE in the BeAGLE repo
source env.sh                # BEAGLESYS, LHAPATH, LD_LIBRARY_PATH
mkdir -p ~/scratch/beagle_test && cd ~/scratch/beagle_test
ln -sf $BEAGLESYS/nuclear.bin .
cp $BEAGLESYS/Examples/eAt1dfJn .                       # PYTHIA card
beagle_run() { $BEAGLESYS/BeAGLE < "$1" > "$1.log"; }
beagle_run ~/Projects/eic/pol_li/tools/beagle/cards/eD_18x135_test.inp   # control first!
beagle_run ~/Projects/eic/pol_li/tools/beagle/cards/eLi7_10x117_draft.inp
```

Validation order (plans/02 step 1.5.3): e+d control vs Tu et al. → fragment
yields for A=6,7 → cluster-IA comparison. Output: PYTHIA-style text on
fort.92? (per `OUTPUT` card) → eic-smear `BuildTree`/`TreeToHepMC` (run
inside eic-shell) → HepMC3.

## No-FLUKA alternatives (available right now)

- **Official samples via xrootd** (re-verified 2026-08-26; the e+d control
  above was produced this way): e+d (`eH2`), e+³He and e+Au EVGEN under
  `root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/DIS/BeAGLE1.03.02-*/`
  — enough to develop the full analysis/conversion chain and do the e+d
  control study before any local BeAGLE run.  One command:
  ```bash
  SIF=~/Projects/eic/local/lib/jug_xl-nightly.sif
  B=root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/DIS/BeAGLE1.03.02-3.1/eH2/en/9x130/q2_1to1000
  singularity exec $SIF python3 tools/analysis/dump_spectators.py \
      $B/BeAGLE1.03.02-3.1_DIS_eH2_en_9x130_q2_1to1000_ab_run001.hepmc3.tree.root \
      ed.csv --nevents 20000                       # ~2 min, streamed
  python3 tools/analysis/ed_control_analysis.py ed.csv --beta-scan
  ```
- **BNL SDCC / JLab ifarm**: BeAGLE prebuilt; submit the e+Li cards there.

## Caveats for e+Li physics (from plans/02 step 1.5)

A>4 uses the C-12 Fermi-momentum parameterization; Woods–Saxon geometry
(no α+d/α+t clustering) — fragment spectra come from FLUKA statistical
de-excitation. Validate, don't trust; the cluster-IA model in
`fastsim/polli_fastsim/spectator.py` brackets the model uncertainty.
