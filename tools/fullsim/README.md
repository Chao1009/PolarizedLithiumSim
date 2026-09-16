# Full-simulation runbook (Phase 2)

*Reproducing the far-forward scans rather than extending them?
[docs/reproduction_manual.md](../../docs/reproduction_manual.md) §5.3 has the
commands and what they should return.*

Local ePIC chain status and recipes. Versions matter — record them in
every output directory.

## Containers

| container | state |
|---|---|
| `~/Projects/eic/local/lib/jug_xl-nightly.sif` | legacy (frozen ~Sep 2024); detector configs `/opt/detector/epic-24.0{5,6,7,8}.0` **and** `epic-main` at git 5a7dd057; ships npsim, eicrecon, uproot, pyHepMC3, xrdcp — the container everything below was run in |
| `~/Projects/eic-2026/local/lib/eic_xl-nightly.sif` | **present since 2026-08-28** (the row said "not present on the analysis box" until then, and the whole aperture measurement below was blocked on it). `epic-main` at git `9aaa296976d3ad9de404f775ae89fc17a068c07c` (2026-08-22), EICrecon main `77d0cef8`; ships npsim, uproot, pyHepMC3. Source `/opt/detector/epic-main/bin/thisepic.sh`, which sets `$DETECTOR_PATH` to the `epic-git.9aaa2969…_main-ur6von3d…/share/epic` prefix; the compact files are `$DETECTOR_PATH/epic_craterlake_{5x41,10x100,18x275}.xml`. **npsim writes `calibrations/`, `fieldmaps/` and `gdml/` caches into the CURRENT DIRECTORY**, so `cd` to a scratch directory inside the `bash -lc` string before running it, never the repository. Cost here: ~47 s of geometry per job plus ~1.3 s/event. |

`./eic-shell` starts `jug_xl`. The *fragment routing* scan below is still
the **September-2024 `epic-main`** geometry; the intact-⁶Li aperture
measurement was repeated in the current container on 2026-08-28 and that
section carries both.

**Known issue (found 2026-06-12, diagnosed 2026-09-15):** the
`eic_xl-nightly` pyHepMC3 `rootIO.ReaderRootTree` **segfaults** where the
legacy container reads the same file — but the bindings are not at fault.
The crash follows *fatal error: module map file
`../../include/root/ROOT.modulemap` not found*: ROOT resolves its module
map through that relative path, so cling finds it only when the process's
working directory sits two levels below `$ROOTSYS = /opt/local`.
`singularity exec --pwd /opt/local/lib/root $SIF python3 ...`, with
absolute paths, reads the file; measured on the BeAGLE `eH2` xrootd tree
of `tools/beagle/README.md` — 20 000 events through
`tools/analysis/dump_spectators.py` — against a segfault from
`--pwd /opt/local` and from the default host directory.  Nothing to
report upstream beyond the ergonomics; npsim/eicrecon were never
affected.

Inside the legacy container, source the geometry with:
```bash
source /opt/detector/epic-24.07.0/setup.sh   # sets $DETECTOR_PATH
```
(current containers use `/opt/detector/epic-main/bin/thisepic.sh`).

## The reconstruction leg — the first EICrecon output (2026-09-15)

The runbook's third clause, `eicrecon -Ppodio:output_file=… sim.edm4hep.root`,
had never been run: until this date no EICrecon output existed anywhere in the
tree. It runs, and what it returns is worth more than the file.

```bash
SIF=~/Projects/eic-2026/local/lib/eic_xl-nightly.sif
S=<scratch>                      # npsim writes its caches into $PWD
python3 tools/fullsim/ion_gun_hepmc.py --out $S/li6_100.hepmc \
    --pdg 1000030060 --a 6 --p-per-nucleon 137.5 \
    --pt $(python3 -c "print(' '.join('%.5f' % (825e-3*(0.2+0.2*i)) for i in range(25)))") \
    --nphi 4 > $S/li6_100.index                       # 100 events
singularity exec --env S=$S $SIF bash -lc '
  source /opt/detector/epic-main/bin/thisepic.sh
  cd $S
  npsim --compactFile $DETECTOR_PATH/epic_craterlake_18x275.xml \
        -N 100 --inputFiles $S/li6_100.hepmc \
        --physics.list FTFP_BERT --part.minimalKineticEnergy "100*MeV" \
        --outputFile $S/sim.edm4hep.root
  eicrecon -Pdd4hep:xml_files=$DETECTOR_PATH/epic_craterlake_18x275.xml \
           -Ppodio:output_file=$S/reco.edm4eic.root $S/sim.edm4hep.root'
```

**Versions, recorded as the runbook asks.** `npsim` **1.8.0**
(`/opt/software/linux-x86_64_v2/npsim-1.8.0-s2yzj4e7osdlvcun2xu6saily2ilpfbr`);
`EICrecon version: git.77d0cef8f035e0f7c78d9eb240f4116a18f5430f=main`; geometry
`epic-git.9aaa296976d3ad9de404f775ae89fc17a068c07c_main-ur6von3dpmaimo6hibxvji7rseapbgc7`
(the `$DETECTOR_PATH` prefix — the epic hash is the one the container row above
carries). Geant4 runs single-threaded (`Multi-threaded mode requested, but not
supported by this compilation`). Cost on the 8-core box: **npsim 6 min 10 s**
for 100 ⁶Li events, **eicrecon 47 s**. Output `reco.edm4eic.root`, 37.5 MB,
**100 events, 1335 collections**.

**`eicrecon` does not inherit npsim's compact file.** `thisepic.sh` exports
`DETECTOR_CONFIG=epic`, so the literal command of plans/03 §2.0 loads
`$DETECTOR_PATH/epic.xml`, whose line 40 includes
`compact/fields/beamline_5x41.xml` — the same trap `ff_gun_scan.sh` encodes on
the simulation side, now on the reconstruction side: 18 × 275 hits reconstructed
against 5 × 41 beamline fields. The `-Pdd4hep:xml_files=…` above is the fix, and
it belongs in the runbook. For *this* input it changed nothing (both runs agree
collection by collection), which is the second finding, below.

**What the reconstruction does with a lithium ion.**

| collection | ⁶Li, 825 GeV | proton control, 275 GeV |
|---|---|---|
| `ForwardRomanPotHits` | 1179 hits, 76/100 ev | 474, 87/100 |
| `ForwardRomanPotRecHits` | 1031, 76 | 443, 87 |
| **`ForwardRomanPotRecParticles`** | **0, 0** | **72, 72** |
| `ForwardRomanPotStaticRecParticles` | 0, 0 | 72, 72 |
| `ForwardOffMRecParticles` | 0, 0 | 0, 0 |
| `B0TrackerHits` / `B0TrackerCKFTracks` | 1402, 17 / 0, 0 | — |
| `HcalFarForwardZDCClusters` | 19, 19 | 6, 6 |
| `ReconstructedHcalFarForwardZDCNeutrals` | 19, 19 | — |

**Digitisation works and reconstruction does not.** 1031 of the 1179 Roman-Pot
sim hits become `RecHits`, and **not one becomes a `RecParticle`**. The control
is the same chain with a 275 GeV proton at the same ring setting, 100 events,
θ = 0.3–3.9 mrad: 72 reconstructed particles at 275.0 ± 1.7 GeV (rms). So the zero is
not a broken chain, a wrong geometry or an empty window — it is the species.
This is `MatrixTransferStaticConfig.h`'s `partMass = 0.938272, partCharge = 1`
(plans/03 §2.2 (2)) with a measured consequence rather than a code reading: the
ePIC far-forward Roman-Pot reconstruction returns **nothing at all** for a
lithium ion, so `fastsim`'s own transport is not an interim stand-in for an
EICrecon path that exists — it is the only path there is.  *Bounded from the
other side on 2026-09-16* (§"The Mode-W chain" below): where the far-forward
species IS a proton — a deuteron beam's spectator, R ≈ 0.47 — the
off-momentum detector reconstructs it in 80 of 100 events at a 5% momentum
and a 0.25 mrad angle, so the matrices are proton-tuned rather than broken. The ZDC neutral chain,
by contrast, reconstructs in every event that reaches it, and the B0 tracker
digitises but seeds no track at these angles.

## Far-forward gun scan (plan 2.2.1)

```bash
SIF=~/Projects/eic/local/lib/jug_xl-nightly.sif
singularity exec $SIF bash ff_gun_scan.sh ~/Projects/eic/data/ff_gun 60
singularity exec $SIF python3 ff_gun_hits.py ~/Projects/eic/data/ff_gun
```

(`ff_gun_hits.py --b0-zdc` adds the B0 tracker/ECal and the two ZDC
collections' hit positions to `--per-event`, which the counting table below
has reported as a yes/no since 2026-06-12 and never as a place. For B0 the
place *is* the IP angle: IP6 → B0 layer 1 is a pure 5.900 m drift at every
ring setting — see the transfer-matrix section below.)

The scan shoots fragments along the ion axis (−25 mrad in x) with total
momenta chosen so the rigidity ratio R relative to the default
craterlake fields (18×275 ep optics) reproduces the Li-fragment cases:

| config | particle | p [GeV] | R | expectation (plans/03 §2.2) |
|---|---|---|---|---|
| alpha_R0857 | α | 471.4 | 0.857 | Roman Pots (⁷Li α-tag) |
| alpha_R100 | α | 550.0 | 1.000 | lost in beam pipe (⁶Li α, pT≈0) |
| alpha_R100pt3 | α | 550.0 | 1.000 (+0.55 mrad ≙ pT=0.3) | RP pT-tail |
| triton_R1286 | t | 353.6 | 1.286 | no coverage |
| proton_R050 | p | 137.5 | 0.500 | OMD |
| deuteron_R0857 | d | 235.7 | 0.857 | Roman Pots |
| neutron_ZDC | n | 117.9 | — | ZDC |

**Results (2026-06-12, epic-main geometry, 18×275 fields, 12 ev/config —
fraction of events with ≥1 hit):**

| config | RP | OMD | B0 | ZDC | verdict vs prediction |
|---|---|---|---|---|---|
| alpha_R0857 | **1.00** | 0 | 0 | 0.17 | ✓ RP (⁷Li α-tag works) |
| alpha_R100 | 0.08 | 0 | 0 | 0 | ✓ invisible (⁶Li α beam-blind) |
| alpha_R100pt3 | **1.00** | 0 | 0 | 0 | ✓ pT-tail recovers it |
| deuteron_R0857 | **1.00** | 0 | 0 | 0.25 | ✓ RP |
| proton_R050 | 0.25 | **1.00** | 0 | 0.42 | ✓ OMD (RP/ZDC = pass-through/splash) |
| triton_R1286 | **1.00** | 0 | 0 | **1.00** | ✗ **surprise: taggable!** |
| neutron_ZDC | 0 | 0 | 0 | **1.00** | ✓ ZDC |

Findings beyond the table:
- Current epic-main RP stations sit at **z ≈ 32.5 / 34.3 m** (papers
  quote 26/28 m — geometry has moved; update plans/03 numbers when the
  preTDR layout is confirmed).
- Offsets from the *measured* 275-GeV orbit (reference-proton gun,
  station-2 plane z ≈ 34.2 m, orbit x = −1212 mm): α(R=0.857) at
  **−48 mm** (dispersion side — outside a 10σ ≈ 18–36 mm exclusion for
  both optics ⇒ ⁷Li α-tag is envelope-safe); **triton (R=1.286) at
  +77 mm on the over-rigid side** of the same RP planes, then deposits
  in the ZDC (bends less than the beam) — far outside any beam
  envelope ⇒ the "no IP6 triton coverage" assumption looks wrong and
  the ⁷Li α+t double-tag may work via RP-inner-side tracking + ZDC
  energy. α(R=1.0) grazes at −15 mm ⇒ inside/marginal vs 10σ ⇒
  beam-blind, as predicted. Pending: divergence/vertex spread, whether
  the inner half is actually instrumented in the real pots, and
  reconstruction-level confirmation.
- Gotcha encoded in the script: plain `epic_craterlake.xml` loads the
  **5×41 beamline fields** — with the 275-optics momenta everything
  flies straight to the ZDC. Use `epic_craterlake_18x275.xml`.
- `epic_craterlake_18x110_Au.xml` (Z/A ≈ 0.40) is the closest existing
  optics proxy for a ⁷Li beam (Z/A = 0.429).

Caveats: single-particle, fixed vertex/direction (no beam envelope —
real pots retract to 10σ), hit-level only. This validates *routing
logic*, not final acceptances (plans/03 step 2.2 proper).

## The intact-⁶Li far-forward scan (2026-08-28) — the pot aperture, measured

The coherent channel's whole acceptance rests on one number: the angle at
which the Roman Pots stop seeing a nucleus at beam rigidity (plans/04 #20).
It was assumed until 2026-08-26, measured then in a **superseded**
geometry, and re-measured on **2026-08-28** in the current one — `epic-main`
at git `9aaa296976d3ad9de404f775ae89fc17a068c07c` inside
`eic_xl-nightly.sif`, which arrived on the analysis box the same day and
on which plans/09 B1 had been blocked. The 2026-08-26 numbers are kept at
the end of this section, dated and marked superseded.

**npsim cannot shoot a nucleus.** `--gun.particle Li6`, `ion(3,6)` and
every other spelling come back *"Geant4ParticleGenerator: Bad particle
type"* — DD4hep looks the name up in the G4 particle table, where generic
ions are not pre-instantiated.  `tools/fullsim/ion_gun_hepmc.py` writes a
one-particle HepMC3 event per scan point instead, which also buys what a
gun cannot do: a scan in **both** p_T and azimuth about the ion axis.
(Two format traps, both recorded in the script: DD4hep sends a `.hepmc`
file to HepMC3's `ReaderAscii`, so a HepMC2 `IO_GenEvent` file gives an
immediate EOF; and a vertex at the origin with nothing incoming must be
left out of the record entirely, with the E-line vertex count 0.)

**The reader had to be fixed first.** Both Roman-Pot stations write into
the *same* readout, `ForwardRomanPotHits`, so `ff_gun_hits.py --positions`
— which averaged over every hit of the event — was reporting the mean of
two planes 1.7 m apart and calling it a position. It now groups by station
plane (S1L1 32547.3, S1L2 32567.3, S2L1 34245.5, S2L2 34265.5 mm, ±15 mm)
and prints n, mean, min and max of x and y per plane, keeping the old
station-mixing row as `ALL`. The min/max are the useful pair in y, where
the 16 mm module tiles the coordinate and a mean is meaningless.

### The commands, exactly as run

All of it ran from a scratch directory, never the repository: **npsim
writes `calibrations/`, `fieldmaps/` and `gdml/` into the current
directory.**

```bash
SIF=~/Projects/eic-2026/local/lib/eic_xl-nightly.sif
S=<scratch>            # e.g. /tmp/…/b1 ; NOT the repo

# (a) the azimuthal acceptance map: 18 angles x 12 azimuths, one event each
python3 tools/fullsim/ion_gun_hepmc.py --out $S/gun.hepmc --pdg 1000030060 \
    --a 6 --p-per-nucleon 137.5 \
    --pt 0.2475 0.33 0.4125 0.495 0.5775 0.66 0.7425 0.825 0.9075 0.99 \
         1.11375 1.2375 1.65 2.0625 2.475 2.8875 3.30 4.125 \
    --nphi 12 > $S/gun.index          # theta = 0.3 … 5.0 mrad; 216 events

# (b) the ladders: theta 0.20-6.00 mrad in 0.05 steps at phi = 0/90/180/270
#     (117 angles x 4 azimuths = 468 events).  p_T = theta * A * p_u, so the
#     --pt list is per configuration.  The ladder below is 18 x 275, where
#     A p_u = 825 GeV IS the fill energy.  At 5 x 41 and 10 x 100 the ladder
#     is shot at the ring's rigidity-scaled reference momentum and not at the
#     fill energy: A p_u = 123.0 and 300.0 GeV, i.e. the rigidity-scaled
#     20.5 and 50 GeV/u every scan here used.  The gamma-matched fill
#     energies at those two settings, 244.8 and 597.0 GeV (plans/10), are
#     not used by any scan in this file.
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
    cd $S            # <-- the caches land here
    npsim --compactFile \$DETECTOR_PATH/epic_craterlake_18x275.xml \
          --inputFiles $S/$NAME.hepmc --numberOfEvents $N \
          --physics.list FTFP_BERT --part.minimalKineticEnergy '100*MeV' \
          --outputFile $S/$NAME.edm4hep.root"

  singularity exec $SIF python3 tools/fullsim/ff_gun_hits.py \
      $S/$NAME.edm4hep.root --per-event --index $S/$NAME.index --positions
done
```

These are **lattice-matching momenta, not fill energies**: the γ-matched
⁶Li menu is 40.8 / 99.5 / 137.5 GeV/u (plans/10 A0), and the
rigidity-scaled pair below is what puts the ion on each ring's *reference
rigidity* so that it follows the reference orbit — which is what makes the
result an aperture measurement rather than a statement about one fill (see
*Which lattice* below).

The rigidity-scaled matching momenta (plans/10, gamma-matched menu apart):
run `epic_craterlake_5x41.xml` at `--p-per-nucleon 20.5`,
`epic_craterlake_10x100.xml` at 50 and `epic_craterlake_18x275.xml` at
137.5, where the two conventions coincide because the ring rigidity caps
6Li there. The single-fragment runs (α at R = 0.857 and 1.0, triton
at R = 1.286, deuteron, proton, neutron) use `--enableGun` at
`--gun.direction "(-0.025,0,1.0)"`, the ion axis; `ff_gun_scan.sh` is that
loop. Cost: ~47 s of geometry per job plus ~1.3 s/event.

### What was measured

Positions are the **in-plane** offset from each plane's own centre. The
station centres are −1131.19 mm (station 1) and −1208.43 mm (station 2)
from `compact/far_forward/roman_pots_eRD24_design.xml`, and both stations
are **rotated by −0.04545 rad about y**, so a plane is not a surface of
constant global z: |dz| = 0.0454|dx| along it, 6.5 mm at the 144 mm module
edge, and layer 2 sits 20 mm along the *normal*, i.e. (−0.91, +19.98) mm
in global (x, z). `ff_gun_hits.py` assigns a hit to the plane whose normal
distance is |w| ≤ 5 mm and reports the in-plane u; grouping on global z
with a ±15 mm window instead, as it did until 2026-08-28, put every hit
beyond |dx| ≈ 110 mm into **both** layers of its station, which inflated
the multiplicity, tripped the analysis' ≤ 3-hit cleanliness cut and
dropped the outermost rows. φ = 0 is +x, the **inner** side of the bend,
where an over-rigid fragment goes.

| configuration | R₁₂ st.1 [m] | R₁₂ st.2 | R₃₄ st.1 [m] | R₃₄ st.2 | D [m] | D₂ [m] |
|---|---|---|---|---|---|---|
| 5 × 41 | **19.24** | 18.54 | **4.56** ‡ | 4.57 ‡ | **0.311** † | −0.190 † |
| 10 × 100 | **21.25** | 20.69 | **3.35** | 3.23 | **0.287** | −0.206 |
| 18 × 275 | **29.97** | 30.31 | **2.93** | 2.62 | **0.292** | −0.215 |

R₁₂, R₃₄ and D are all **station 1**, the plane the aperture thresholds
are divided by. † at 5 × 41 the R = 0.857 α reaches no station-1 silicon
at all (the 32–48 mm band is held off to |y| ≥ 18 mm there), so only
station 2 has the three rigidities a quadratic needs and D / D₂ are
station-2 numbers. The two stations differ by 8% in D and up to 14% in D₂
— 0.292 / −0.215 at station 1 against 0.315 / −0.249 at station 2 at
18 × 275 — which is real, station 2 having the longer lever, and is why
the pair must be read from one station and not mixed.

**R₁₂ = dx/dθ_x** is regressed on the φ = 0/180° ladders over the whole
accepted range, 46–138 primary rows per plane, by a Theil–Sen slope
followed by least squares on the 2 mm inliers (the outliers are showers,
not tracks). 42–141 primary rows per plane, residual rms 0.67–0.82 mm,
worst 2.4 mm. The fitted intercept is +0.17 / +0.50 / +0.48 mm at station
1 and within 0.6 mm of zero in **every** plane: the array is centred on
the reference orbit to under a millimetre, which is the check that the
regression measures transport and not a mis-set reference. (Before the
rotation-aware plane assignment the layer-2 intercepts sat at −0.4 to
−0.7 mm, which was the unmodelled −0.91 mm layer-2 x offset and not the
orbit. The levers themselves moved by 0.04–0.1%, an order of magnitude
inside the fit residual.)

**R₃₄ = dy/dθ_y** comes from the φ = 90/270° ladders restricted to rows
with |dx| < 16 mm, so that the vertical threshold is the single central
one. Intercepts −0.81 mm (10 × 100) and −0.11 mm (18 × 275), rms 0.9 mm.

‡ **The 5 × 41 vertical lever needed a second geometry** (2026-08-29).
The per-energy insertion there holds the central silicon off to |y| ≥ 29.6
mm and the intermediate bands to 27.5 and 18.0 mm, so the whole 0.20–6.00
mrad vertical ladder returned **one** accepted row (θ = 3.35 mrad, y =
+30.94 mm, a single hit) at φ = 90° and none at φ = 270°, and one point is
not a regression. The transport, however, belongs to the magnets and not
to where the pots sit, so it can be read off a geometry whose pots are not
retracted. Copy `compact/fields/beamline_5x41.xml` and
`epic_craterlake_5x41.xml` into a scratch directory, set the four
insertion constants of the copied beamline file to zero —

```xml
<constant name="offset_central_RP_section"        value="0.0*cm"/>
<constant name="offset_intermediate_1_RP_section" value="0.0*cm"/>
<constant name="offset_intermediate_2_RP_section" value="0.0*cm"/>
<constant name="offset_outer_RP_section"          value="0.0*cm"/>
```

— point the copied `epic_craterlake_5x41.xml` at it,

```bash
singularity exec $SIF bash -lc 'source /opt/detector/epic-main/bin/thisepic.sh
  mkdir -p '$S'/r34
  sed "s#\(offset_[a-z0-9_]*_RP_section\" value *=\)\"[^\"]*\"#\1\"0.0*cm\"#" \
      $DETECTOR_PATH/compact/fields/beamline_5x41.xml > '$S'/r34/beamline_5x41_zero.xml
  sed "s|\${DETECTOR_PATH}/compact/fields/beamline_5x41.xml|'$S'/r34/beamline_5x41_zero.xml|" \
      $DETECTOR_PATH/epic_craterlake_5x41.xml > '$S'/r34/epic_craterlake_5x41_zero.xml'
```

(`offset_OMD_x` does not match `_RP_section` and is left alone; the four
Roman-pot section offsets are the only constants touched.  Both
substitutions are written to be exact: the value pattern is `"[^"]*"`
rather than `.*`, so the two `<!-- rough extrapolation -->` comments on
lines 63–64 survive, and the include pattern carries its `${DETECTOR_PATH}/`
prefix, without which the rewritten path is appended to the detector
directory and no such file exists.  Both files then come out byte-identical
to the stored ones.)

and run the same θ = 0.20–6.00 mrad ladder in 0.05 mrad steps at
φ = 0/90/180/270° through it with an intact ⁶Li at the ring's reference
rigidity — the rigidity-scaled matching momentum of the paragraph above,
not a fill energy — i.e. `ion_gun_hepmc.py --a 6 --p-per-nucleon 20.5 --pt
$(python3 -c "print(' '.join('%.5f' % (123e-3*(0.20+0.05*i)) for i in
range(117)))") --nphi 4`, then `npsim --compactFile
$S/r34/epic_craterlake_5x41_zero.xml`, then `ff_gun_hits.py --per-event
--index … --positions`; 123 GeV total, 468 events. Every
field is untouched, so the transport is the baseline one; only the silicon
has slid onto the beam axis. Regressing y on θ_y over the φ = 90/270°
rows then gives **R₃₄ = 4.56 m** at station 1 (4.564 at layer 1, 4.558 at
layer 2, 149 of 154 rows on the line, residual rms 0.79 mm, stable to 1.5%
when the fit is restricted to |θ_y| < 3 mrad) and 4.57 m at station 2.
The control is R₁₂ in the same file: **19.18 m** against the 19.24
measured through the real insertion, i.e. zeroing the offsets moved the
transport by 0.3%, which is the fit scatter and not a change.

The plane is nonetheless *shut*, and now for a reason with a number
attached: reaching the 29.6 mm insertion needs θ_y = 29.6 mm / 4.56 m =
**6.49 mrad**, which is past the 5 mrad at which the far-forward routing
ends. It is *not* shut by the 2.85 mrad outer edge of the row above —
that is a **horizontal** edge, 54.8 mm of x on the 19.24 m horizontal
lever, while 6.49 mrad of y is only 29.6 mm on the 4.56 m vertical one;
the levers are 4.2× apart and the two angles do not compare. This very
ladder carries clean on-the-fit-line vertical rows out to θ_y = 4.65
mrad. What it does not carry is a row *at* the insertion: it stops at
6.00 mrad, i.e. y = 4.56 × 6.00 − 2.4 = 25 mm against 29.6, so the
constant is an extrapolation of a measured **lever**, not a measured
**edge**. The −2.42 mm fitted vertical orbit offset splits the two sides
to 7.02 mrad upward and 5.96 downward; 6.49 is the offset-free value, and
nothing this programme produces reaches either. That 6.49 mrad is what
`reco.RP_APERTURE_MEASURED["5x41"]` carries as c_y, replacing the 8.84
mrad bound (29.6 mm over the largest vertical lever seen anywhere,
R₃₄ = 3.35 m at 10 × 100) it held while the lever was unmeasured.

**D = dx/dR at θ = 0** from a quadratic through the three rigidities each
configuration provides — α at R = 0.857 and triton at R = 1.286 from the
θ = 0 frag guns, intact ⁶Li at R = 1.000 from the ladder's fitted
intercept. `mkdisp.py` in the scan directory is that fit, and D is the
linear coefficient at R = 1, which is what `farforward.POT_DISPERSION`
means. The second-order term is large and may not be dropped at the
triton's rigidity: a straight line through the same three points reads
0.25–0.28 m and misses the R = 1 slope by 12–13%. The result **confirms
the repository's 0.30 m to 4% at every configuration**, and it reproduces
the measured triton displacement to a tenth of a millimetre at 18 × 275:
0.292 × 286 − 0.215 × 286²/1000 = 65.9 mm against 65.95 measured (the
+66.43 mm hit minus the +0.48 mm intercept). An independent check comes
from the He4 lattice run below, where a ⁶Li at 81.6 GV in an 82 GV lattice
is off rigidity by −0.49% and its ladder intercept is −1.13 mm, i.e.
D ≈ 0.23 m — the same quantity to 20% from a different measurement.

The transport is **4.2× stiffer in x than in y at 5 × 41, 6.3× at
10 × 100 and 10.2× at 18 × 275** — the ratio grows with energy. That, and
not the pot orientation, is why the open aperture is a horizontal slot
although the pots insert vertically.

### The edges

The silicon is |x| ≤ 144 mm in four bands — central ≤ 16, intermediate
16–32, intermediate_541 32–48, outer 48–144 — each starting at |y| = its
per-energy offset (2.96 / 2.75 / 1.80 / 0 cm at 5 × 41, 0.71 / 0.55 / 0 /
0 at 10 × 100, 0.27 / 0 / 0 / 0 at 18 × 275). A track at y ≈ 0 therefore
needs |dx| > 48 / 32 / 16 mm.

| configuration | horizontal edge +x / −x | threshold/R₁₂ | vertical edge +y / −y | offset/R₃₄ | outer edge +x / −x | 144 mm/R₁₂ |
|---|---|---|---|---|---|---|
| 5 × 41 | **2.50 / 2.50 mrad** | 48/19.24 = 2.50 | none in 0.2–6.0 | 29.6/4.56 = 6.49 | **4.30 / 2.85 mrad** | 7.48 |
| 10 × 100 | **1.50 / 1.55** | 32/21.25 = 1.51 | 2.50 / 1.95 | 7.1/3.35 = 2.12 | **4.30 / 3.85** | 6.78 |
| 18 × 275 | **0.55 / 0.55** | 16/29.97 = 0.53 | 0.80 / 0.95 | 2.7/2.93 = 0.92 | **4.25 / 4.00** | 4.80 |

**The two readings close on each other at all three configurations at
once, which is what makes this a measurement.** The first angle with a
Roman-Pot hit equals the band threshold divided by the independently
regressed lever, to 4% at 18 × 275 and 1% at the other two, horizontally;
and to the same at 18 × 275 vertically. At 10 × 100 the two vertical sides
straddle the geometric 2.12 mrad (2.50 up, 1.95 down) and the −0.81 mm
fitted vertical orbit offset accounts for the split exactly:
(7.1 ± 0.81)/3.35 = 2.36 and 1.88 mrad. The 5 × 41 vertical row is the
one place the two readings cannot be compared: there is no measured
vertical edge to compare the 6.49 mrad against, because the ladder stops
at 6.00 mrad (y = 25 mm) before the 29.6 mm insertion is reached. The
2.85 mrad of the outer-edge column is **not** that comparison — it is a
horizontal edge on a 4.2× longer lever. What makes the plane *shut* is
that 6.49 mrad is past the 5 mrad at which the far-forward routing ends,
and the zero-insertion scan of 2026-08-29 is what turns that from an
inference into an arithmetic one.

**The outer edge is not the module edge.** At every configuration the
contiguous primary track stops between 2.85 and 4.30 mrad, at |dx| =
54–127 mm, well inside the last module, and beyond it the rows carry B0
and off-momentum hits and off-line positions — the ion has struck the pipe
or the magnet aperture and what the pots see is debris. The −x
(under-rigid, outer-bend) side is cut earliest: at 5 × 41 the primary
survives 0.35 mrad past the inner edge and no further. Both the "θ < 5
mrad" of the acceptance tables and the "144 mm / R₁₂" of the module
arithmetic are too generous, by 1.1–2.6 mrad. `farforward.
THETA_RP_OUTER_MEASURED` carries the smaller side per configuration.

**What that constant is, exactly.** It is the last angle of a *debris-free
contiguous run*: the run breaks at a gap of more than 0.30 mrad, and every
row in it must lie within 2 mm of the fitted transport line **and** carry
at most three hits in the plane. It is **not** the last row with a
primary, and the difference is the cleanliness cut rather than the
transport. At 5 × 41, φ = 180°, θ = 2.90 mrad the station-1 hit is at
dx = −55.0 mm against a fitted −55.6 — a primary — but comes with 13 hits
in the plane, so the run stops at 2.85. At 18 × 275, φ = 180°, θ = 4.35
mrad all four planes give a consistent on-line row (−130.8 / −131.2 /
−132.6 / −132.8 mm) beyond the 4.00 mrad quoted. Isolated on-line rows
survive 0.05–0.35 mrad past the constant, so it is the last angle at which
a *single clean track* is reconstructed through debris — which is the
quantity a partner-fragment veto needs — and it is conservative by up to
one 0.05 mrad step plus one gap.

**ZDC-only rows.** At 5 × 41 the twelve consecutive rows from θ = 1.80 to
2.45 mrad on both horizontal sides — the whole band between the beam and
the silicon edge — show ZDC hits and no Roman-Pot hit: the 123 GeV ⁶Li is
breaking up upstream and its neutrons reach the calorimeter. At the other
two configurations the ZDC-only rows begin only past the outer edge.

### The acceptance map at 18 × 275

Roman-Pot hit (●) for an intact ⁶Li at the ring reference rigidity, by
transverse angle at the IP and azimuth about the ion axis (φ = 0 is the
horizontal bend plane, 90° is up); `z` = ZDC only.

| θ [mrad] | 0 | 30 | 60 | 90 | 120 | 150 | 180 | 210 | 240 | 270 | 300 | 330 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.30 | | | | | | | | | | | | |
| 0.40 | | | | | | ● | | ● | | | | |
| 0.50 | | | | | | | ● | | | | | |
| 0.60 | ● | ● | | | | ● | ● | ● | ● | | | ● |
| 0.80 | ● | ● | | | | ● | ● | ● | | ● | | ● |
| 1.00 | ● | ● | ● | ● | | ● | ● | ● | ● | ● | | ● |
| 1.10 | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● |
| 1.50 | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● |
| 3.00 | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● |
| 4.00 | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● |
| 5.00 | z | z | z | ● | z | z | ● | z | z | z | z | ● |

The 12-azimuth map is coarse; the 0.05 mrad ladders are what the edges
above are read from. At 10 × 100 the same map opens at 1.5–2.0 mrad and at
5 × 41 at 2.5–3.0 mrad, and neither has vertical acceptance below 2 and
6 mrad respectively.

### The aperture, for the code

`reco.RP_APERTURE_MEASURED` is the cutout half-widths (c_x, c_y) — a
recoil is seen when |θ_x| > c_x **or** |θ_y| > c_y — taken as the band
threshold over the station-1 lever:

| configuration | c_x [mrad] | c_y [mrad] | superseded (Sep 2024) |
|---|---|---|---|
| 5 × 41 | **2.50** | **6.49** — measured, and the plane is *shut* | 2.0 × 3.0 |
| 10 × 100 | **1.51** | **2.12** | 1.35 × 3.0 |
| 18 × 275 | **0.53** | **0.92** | 1.03 × 2.3 |

As p_T at the γ-matched ⁶Li fill (40.8 / 99.5 / 137.5 GeV/u, i.e. 244.8 /
597 / 825 GeV — the fill is γ-matched, not rigidity-scaled, plans/10 A0),
and against the horizontal half-width of the 10σ Yellow Report
high-acceptance envelope, which is 2.2 / 1.8 / 0.9169 mrad:

| configuration | c_x | as p_T | vs 10σ_h envelope |
|---|---|---|---|
| 5 × 41 | 2.50 mrad | 0.61 GeV | **1.14×** — the silicon binds |
| 10 × 100 | 1.51 | 0.90 | 0.84× — the beam binds |
| 18 × 275 | 0.53 | 0.44 | 0.58× — the beam binds, by a lot |

The ordering the September-2024 numbers gave (0.91× / 0.75× / 1.12×, the
beam binding at the two lower configurations and the silicon marginally at
the top) is **reversed at both ends**.

At 18 × 275 both half-widths shrink — x by 1.9×, y by 2.5× — which is the
16 × 16 mm module and the 2.7 mm insertion replacing the 32 mm module and
7 mm insertion of September 2024, exactly as plans/09 §9.4 predicted from
reading the files (it named 0.52 mrad from 16 mm / 30.6 m; the scan gives
0.53). At 5 × 41 the aperture gets **worse**, also as predicted. The
cutout is **taller than it is wide at every configuration**, by 1.4–2.6×,
so `reco.rp_measure`'s `cut_scale_xy = (2.5, 1.0)` still has the aspect
the wrong way round — by 3.5–6.5× rather than the 5.8× the September-2024
numbers gave. The magnitude moved; the **sign** of the acceptance-induced
⟨cos 2φ_t⟩, which is the half of this that does not depend on the geometry
version, did not.

### The triton: the "no coverage above R = 1.05" assumption is wrong

Report 3 Table 6 records "no coverage" for an over-rigid fragment as a
*routing assumption*: `farforward.route_charged` carried no R > 1 branch
at all. It has now been measured at every configuration.

| configuration | fragment | R | RP | ZDC | dx st.1 | dx st.2 |
|---|---|---|---|---|---|---|
| 5 × 41 | triton | 1.286 | **1.00** | 0.80 | +66.2 mm | +73.2 mm |
| 10 × 100 | triton | 1.286 | **1.00** | 0.83 | +65.8 | +73.1 |
| 18 × 275 | triton | 1.286 | **1.00** | 0.98 | +66.4 | +70.1 |
| 5 × 41 | α | 0.857 | **0.57** | 0.33 | no hit | −48.7 |
| 10 × 100 | α | 0.857 | 1.00 | 0.38 | −44.8 | −49.1 |
| 18 × 275 | α | 0.857 | 1.00 | 0.43 | −45.8 | −49.8 |
| 18 × 275 | α | 1.000 | 0.03 | 0.00 | — | −3.4 (2/60) |
| 18 × 275 | α, +0.3 GeV p_T | 1.000 | 1.00 | 0.05 | +18.7 | +18.5 |
| any | proton | 0.500 | 0.15–0.38 | 0.40–0.58 | OMD 1.00 | — |
| 18 × 275 | neutron | — | 0 | 1.00 | — | — |

An **R = 1.286 triton is on the silicon in 60 of 60 events at every
configuration**, on the inner side of the bend, at dx = +66 mm at station
1 and +70 to +72 at station 2 — inside the 48–144 mm outer band, which
carries no vertical insertion anywhere, so the hit is at y ≈ 0 and needs
no p_T at all. It then deposits in the ZDC in 80 / 83 / 98% of the same
events. The ePIC signature of a ⁷Li α + t double tag is therefore an
RP-inner track plus calorimeter energy, and `route_charged` has an
over-rigid branch (code 6, `RP-inner (over-rigid)`) since 2026-08-28.

Its angular reach, from the θ = 0–3 mrad triton ladder at 18 × 275: on the
+x side it runs to 2.25 mrad (dx = +142.6 mm) and is lost at 2.50, which
*is* the 144 mm module edge; on the −x side the angular term eats the
dispersive one at 32.5 mm/mrad — the triton's own R₁₂, 8% larger than the
beam's because a stiffer particle is focused less — so it crosses into the
16 mm blind block at θ_x = −1.55 mrad and reappears at −2.53 mrad. There
is a real **acceptance hole** between them, which no purely angular
routing model can see.

**The ⁷Li α tag at 5 × 41 is a knife edge.** At 10 × 100 and 18 × 275 the
α at R = 0.857 is a clean 100% hit at dx = −45 to −50 mm, in the 32–48 mm
band whose vertical offset is zero at both — "mid-window" as Report 3 says.
At 5 × 41 that band is held off to |y| ≥ 18.0 mm, so an α at y ≈ 0 sees no
silicon at station 1 at all, and the tag survives only because station 2's
slightly larger lever carries it to −48.3 mm, a third of a millimetre past
the boundary into the offset-free outer band: **57%, station 2 only, on a
module edge**. The deuteron at the same rigidity behaves identically (58%).
The fast simulation routes in angle and knows nothing of the per-band
vertical offsets, so it cannot see this.

**Species check.** At 18 × 275 the intact ⁶Li (825 GeV) and an α at the
same rigidity (550 GeV) give their first Roman-Pot primary at **θ = 0.55
mrad on both horizontal sides**, dx = +17.7 / −16.9 mm for the ⁶Li and
+18.4 / −17.3 for the α. Same angle, different species, same rigidity: the
scan measures rigidity-driven optics and not the nucleus. (September 2024
read 0.909–1.091 mrad for the α against 0.970–1.091 for the ⁶Li on a
coarser grid; the two now agree to one 0.05 mrad step.)

### Which lattice do these numbers stand on?

`epic-main` ships **two different lattices** for the 5 × 41 ring setting,
and they are not a field scale of one another.
`compact/fields/beamline_5x41.xml` is the proton one — 41 GeV,
`FieldScaleFactor` 1.0, its own comment *"the values from Scott's tables
for the 50cm shifted hadron FF 41 GeV proton lattice. Notice they mostly
do not scale with beam energy in the straightforward way from the previous
lattice."* `beamline_5x41_He4.xml` is the Z/A = 0.5 one at exactly the ⁶Li
fill point — `IonBeamEnergy` 164 GeV (82 GV), `FieldScaleFactor` 82/275,
its own comment *"All magnet values are input by hand and represent the
implementation in EICRoot used for the Yellow Report."* Their Roman-Pot
offsets are byte-identical, so the difference is pure transport;
`Q1APF_GradientMax` is −15.38 T/m in one and −72.61 T/m in the other,
whose 82/275 scale would be −21.65.

No `epic_craterlake_5x41_He4.xml` is built. Every include of
`epic_craterlake_5x41.xml` is an absolute `${DETECTOR_PATH}/…`, so the
compact file can be copied anywhere and one line changed — **no bind mount
and no overlay are needed**:

```bash
singularity exec $SIF bash -lc 'source /opt/detector/epic-main/bin/thisepic.sh
  mkdir -p '$S'/he4
  sed "s|compact/fields/beamline_5x41.xml|compact/fields/beamline_5x41_He4.xml|" \
      $DETECTOR_PATH/epic_craterlake_5x41.xml > '$S'/he4/epic_craterlake_5x41_He4.xml'
# then npsim --compactFile $S/he4/epic_craterlake_5x41_He4.xml …
# (docs/reproduction_manual.md drops the he4/ component and writes into $S)
```

Shot with ⁶Li at 40.8 GeV/u = 244.8 GeV = 81.6 GV, the γ-matched fill
point, θ = 0.2–6.0 mrad in 0.1 steps at φ = 0/90/180/270:

| quantity | proton 5 × 41 (41 GV) | He4 5 × 41 (82 GV) | 18 × 275 |
|---|---|---|---|
| R₁₂ station 1 | 19.24 m | **29.81 m** | 29.97 m |
| R₁₂ station 2 | 18.54 m | **30.11 m** | 30.31 m |
| horizontal edge +x / −x | 2.50 / 2.50 mrad | **1.60 / 1.70 mrad** | 0.55 / 0.55 |
| 48 mm / R₁₂ | 2.50 mrad | 1.61 mrad | (16 mm band) |
| outer edge +x / −x | 4.30 / 2.85 mrad | 4.30 / 4.10 mrad | 4.25 / 4.00 |
| vertical acceptance | none in 0.2–6.0 | none in 0.2–6.0 | 0.80 / 0.95 mrad |

**The transport is not a pure field scale between the two files and the
edges do not agree: 2.50 mrad against the 1.61 the He4 lever implies
(48 mm / 29.81 m; first hit at 1.60 on +x, 1.70 on −x), a factor 1.55.** What the He4
lattice reproduces instead is the **18 × 275** lattice, to 1% in R₁₂ —
which is what it is, the Yellow Report 275 GeV magnet set scaled by 82/275,
and a pure field scale *does* preserve the transfer matrix. The ePIC
5 × 41 and 10 × 100 files are the newer, separately re-matched lattices
(R₁₂ = 19–21 m) and are a different optics, not a scaled one.

So the lattice has to be stated with the number, because it is worth a
factor 1.55 in the horizontal aperture at 5 × 41. Everything above is on
the ePIC baseline compact files, which is what every earlier scan here
used and what `docs/reproduction_manual.md` documents.

**Which file a ⁶Li fill at 81.6 GV would run in is not settled, and since
2026-08-29 this repository states a working assumption instead of leaving
the question open.** The assumption — an educated guess of ours, not a
machine-design statement — is the **ePIC baseline lattice of that ring
setting run at twice the field** for the Z/A = ½ fill: the same magnets,
the same orbit, and therefore the transport measured above. A lithium
beam at 40.8 GeV/u carries 81.6 GV against the 41 GeV proton lattice's
41 GV, so every gradient doubles, the optics is unchanged, and R₁₂, R₃₄
and D are the numbers in the table. The Yellow-Report-scaled Z/A = 0.5
file is carried as the **alternative**, and the ratio of the two is the
systematic: R₁₂ ×1.55, the horizontal edge ×0.64 at 5 × 41. The
far-forward working group may settle it the other way, and plans/09 B1
records what would then move.

**The alternative is reachable in code, not only here.**
`farforward.POT_LEVERS_LIGHT_ION_LATTICE` carries (R₁₂, R₃₄, D) =
(29.81, None, 0.23) at 5 × 41, `None` at 10 × 100 (it cannot be checked)
and the baseline triple at 18 × 275 (the 18 × 275 compact file *is* the
Yellow Report magnet set the He4 file scales);
`reco.RP_APERTURE_MEASURED_LIGHT_ION_LATTICE` carries the aperture that
follows, 1.61 mrad at 5 × 41, and is reachable through
`reco.rp_aperture_for(cfg, table=…)` exactly as `RP_APERTURE_SEP2024` is.
Both are pinned in `evgen/tests/test_nearbeam.py`, so a downstream caller
can price the factor 1.55 rather than read about it.

**10 × 100 cannot be checked the same way.** The only Z/A = 0.5 file near
that setting is `beamline_10x110_H2.xml` at 220 GV, and the ⁶Li fill at
10 × 100 is at 199 GV — a different *energy*, not a different lattice at
the same energy, so a scan there would confound the two. No 199 GV
light-ion file exists in `epic-main`. The check is not meaningful at
10 × 100 and was not run; carry the 5 × 41 factor 1.55 as the size of the
effect.

### Caveats

One event per scan point, no beam divergence and no vertex spread, hit
level only (no reconstruction), and the high-divergence insertion is
present in the XML but commented out, so this is the **high-acceptance**
configuration alone. The over-rigid branch's own criterion is the y ≈ 0
one: a fragment displaced in *both* x and y can reach silicon inside the
blind block where the model says it cannot, and `over_rigid_route` uses the
BEAM's R₁₂ for the angular term where a stiffer fragment's is up to 8%
larger. What a closer approach is *worth* does not depend on any of it —
`evgen/scripts/nearbeam_aperture_scan.py` prices every aperture.

### Superseded: the 2026-08-26 measurement

> Run in `jug_xl-nightly.sif`, whose `epic-main` is git 5a7dd057
> (≈ September 2024): **32 × 32 mm** modules, a single energy-independent
> insertion (3.2 cm outer, +0.7 cm central), and 1 mm aluminium RF shields
> active per module face. It gave a horizontal edge of **2.00 / 1.35 /
> 1.03 mrad** at 5 × 41 / 10 × 100 / 18 × 275 against a vertical
> **3.0 / 3.0 / 2.3 mrad**, one lever **R₁₂ = 30.6 m** read off four hit
> positions at 18 × 275 and applied at every configuration, no R₃₄ at all,
> and D = 0.30 m. Its own internal check was the 32 mm block over that
> lever, 1.046 mrad against 1.03 measured — 1.5%. Every one of those
> numbers is superseded by the table above; the 30.6 m survives to 2% at
> 18 × 275 and is 60% high at 5 × 41, the 2.3 mrad vertical at 18 × 275 is
> 2.5× the truth, and the 1.03 mrad horizontal is 1.9×. What survived
> intact: D = 0.30 m, the horizontal aspect of the slot, and the
> observation that a horizontal IP angle is imaged with a far larger lever
> than a vertical one — the September reading put that at "about ten
> times" and the measurement makes it 9.
>
> Two of its conclusions are withdrawn rather than refined. "At the
> published optics the machine binds at every configuration" is no longer
> true at 5 × 41, where the aperture widened to 2.50 mrad and the silicon
> is now the binding constraint (9.4×10⁻¹⁰ against the envelope's
> 7.2×10⁻⁸); and it is no longer true at 18 × 275 in the other direction,
> where the aperture halved to 0.53 mrad and the silicon is now nine
> orders of magnitude *more* permissive than the beam (1.2×10⁻⁵ against
> 7.1×10⁻¹⁴). The middle configuration is unchanged in kind. The
> superseded coherent-tag table it carried — 9.8×10⁻⁷ / 7.7×10⁻¹⁶ /
> 1.9×10⁻¹⁷ at the measured aperture — is replaced by
> **9.4×10⁻¹⁰ / 2.0×10⁻¹⁹ / 1.2×10⁻⁵** (`nearbeam_aperture_scan.py`).

## The transfer matrix (2026-09-15) — the second row, measured

Everything in the section above is the **first row** of the IP6 → Roman-Pot
transfer: R₁₂ = dx/dθ_x, R₃₄ = dy/dθ_y and D = dx/dR are all
*d(position at the pot)/d(something at the IP)*. Writing the ion line to one
pot plane in the usual form,

    x_pot   =  R₁₁ x_IP + R₁₂ x′_IP + D  δ
    x′_pot  =  R₂₁ x_IP + R₂₂ x′_IP + D′ δ

what was never measured was anything that needs either a **displacement at the
IP** (R₁₁, R₂₁) or an **angle at the pot** (R₂₁, R₂₂, D′). plans/03 §2.2 (2)
recorded them as outstanding. They are measured here, with the vertical
R₃₃/R₄₃/R₄₄ alongside, in the same geometry and with the same
Theil–Sen-then-least-squares fit the published levers use.

**Two things made it possible that August did not have.** `ion_gun_hepmc.py`
now writes a HepMC3 **vertex**, so (x₀, y₀) at the IP is a scan coordinate like
(p_T, φ) — the record DD4hep's reader accepts is a status-4 beam copy of the ion
with production vertex 0, then `V -1 0 [1] @ x y z t`, then the status-1 ion
(a *bare* `V` line with nothing incoming fails the parse; see `write_hepmc3`).
And `ForwardRomanPotHits` is an edm4hep `SimTrackerHit` collection carrying
**`momentum` as well as `position`**, so x′_pot is *read* at station 1 — the
plane R₁₂, R₃₄ and D are quoted at — instead of being differenced between
stations. The same branch identifies the primary: `|p| > 0.5 p_beam` separates
an 825 GeV ⁶Li from everything it knocks out, which is a cleaner cut than the
hit-count cleanliness one the position-only fits needed.

**Through a zero-insertion geometry.** A transfer matrix belongs to the magnets,
and the per-energy insertion holds the silicon 16 / 32 / 48 mm off the axis
horizontally at y ≈ 0, so an ion at θ = 0 with a millimetre of IP offset is
never seen. `tools/fullsim/ff_zero_insertion.sh` is the recipe of the R₃₄
paragraph above, now executable: the four `offset_*_RP_section` constants set to
`0.0*cm`, every field untouched.

### The commands, exactly as run

```bash
SIF=~/Projects/eic-2026/local/lib/eic_xl-nightly.sif
S=<scratch>                      # NOT the repo: npsim caches land in $PWD
R=<this repository>              # absolute; the container sees the host paths
for C in 18x275 10x100 5x41; do
  python3 $R/tools/fullsim/ff_transfer_scan.py build --config $C --out $S/x$C
  python3 $R/tools/fullsim/ff_transfer_scan.py build --config $C --out $S/x$C \
      --stem scan_topup --legs x0,y0 --topup          # 16 more rows for R21
done
singularity exec --env S=$S --env R=$R $SIF bash -lc '
  source /opt/detector/epic-main/bin/thisepic.sh
  cd $S
  for C in 18x275 10x100 5x41; do
    bash $R/tools/fullsim/ff_zero_insertion.sh $C $S/zero
    for STEM in scan scan_topup; do
      npsim --compactFile $S/zero/epic_craterlake_${C}_zero.xml \
            --inputFiles $S/x$C/$STEM.hepmc --numberOfEvents $(grep -c "^E " $S/x$C/$STEM.hepmc) \
            --physics.list FTFP_BERT --part.minimalKineticEnergy "100*MeV" \
            --outputFile $S/x$C/$STEM.edm4hep.root
    done
  done'
for C in 18x275 10x100 5x41; do          # the fit runs on the HOST (uproot)
  python3 $R/tools/fullsim/ff_transfer_scan.py fit $S/x$C --config $C
done
```

158 events per configuration (142 + 16), ~7 min of npsim each. The scan legs are
`thx`/`thy` (θ = 0.2–4.0 mrad in 0.2 steps at φ = 0/180 and 90/270 — the
published ladders, and the controls), `x0`/`y0` (θ = 0, vertex walked to
±32 mm), `dlt` (θ = 0, δ = 0 … ±0.15) and `b0x`/`b0y` (θ = 6–20 mrad, the B0
window). ⁶Li at each ring's **reference rigidity** — 123 / 300 / 825 GeV total,
the rigidity-scaled matching momenta, not the γ-matched fill menu.

### What came out

Station 1 layer 1. Every entry is a slope ± its standard error; the percentage
is that error over the slope, which is the quantity to read as the fit residual.

| configuration | R₁₁ | R₂₁ [rad/m] | R₂₂ | D′ [rad] |
|---|---|---|---|---|
| 5 × 41 | **1.148** ± 1.8% | **−0.0837** ± 2.8% | **−0.4955** ± 1.2% | **0.0175** ± 2.6% |
| 10 × 100 | **1.227** ± 0.8% | **−0.0651** ± 0.8% | **−0.3060** ± 0.9% | **0.0179** ± 2.4% |
| 18 × 275 | **1.852** ± 1.1% | **−0.0209** ± 4.5% | **0.1944** ± 1.5% | **0.0182** ± 4.0% |

| configuration | R₃₃ | R₄₃ [rad/m] | R₄₄ |
|---|---|---|---|
| 5 × 41 | **−4.199** ± 0.7% | **−0.2186** ± 1.7% | **0.0055** ± 0.0040 |
| 10 × 100 | **−3.711** ± 0.3% | **−0.2028** ± 0.5% | **−0.0856** ± 3.8% |
| 18 × 275 | **−3.031** ± 0.4% | **−0.1778** ± 0.4% | **−0.1512** ± 2.2% |

**The controls close, which is what makes it a measurement.** The three
published levers come out of the *same files* and reproduce `POT_LEVERS`:

| configuration | R₁₂ measured | carried | R₃₄ measured | carried | D measured | carried |
|---|---|---|---|---|---|---|
| 5 × 41 | 19.186 ± 0.061 | 19.24 (−0.3%) | 4.617 ± 0.055 | 4.56 (+1.3%) | 0.2845 ± 0.0078 | 0.311 ‡ |
| 10 × 100 | 21.341 ± 0.052 | 21.25 (+0.4%) | 3.255 ± 0.064 | 3.35 (−2.8%) | 0.2818 ± 0.0074 | 0.287 (−1.8%) |
| 18 × 275 | 29.961 ± 0.054 | 29.97 (−0.03%) | 3.081 ± 0.064 | 2.93 (+5.2%) | 0.2875 ± 0.0115 | 0.292 (−1.5%) |

‡ the one 8% outlier, and the expected one: the carried 5 × 41 D is a
**station-2** fit, because through the real insertion the R = 0.857 α reaches no
station-1 silicon there. The zero-insertion file has station 1, and the two
stations differ by 8% in D — exactly the spread `POT_DISPERSION_2`'s note
records. The 19.186 m against 19.24 is the same 0.3% zero-insertion control the
2026-08-29 R₃₄ run reported, reproduced independently.

**And the blocks come out symplectic without having been fitted to be.** A
linear transfer at fixed rigidity has determinant exactly 1:

| configuration | R₁₁R₂₂ − R₁₂R₂₁ | R₃₃R₄₄ − R₃₄R₄₃ |
|---|---|---|
| 5 × 41 | 1.038 | 0.986 |
| 10 × 100 | 1.013 | 0.978 |
| 18 × 275 | 0.985 | 1.006 |

(in-file values; with the *carried* R₁₂/R₃₄ substituted — which is what
`farforward.pot_transfer_for` returns — 1.042 / 1.008 / 0.986 and
0.974 / 0.997 / 0.979). Four numbers measured on four different legs of the
scan, closing on 1 to 4% at three ring settings, is a stronger statement than
any one of the residuals.

**Three things the numbers say.**

*R₂₂ changes sign* between 10 × 100 and 18 × 275: −0.496, −0.306, **+0.194**.
The horizontal phase advance from IP6 to the pots crosses a node as the energy
rises, so the top configuration is the only one where a positive angle at the IP
is still a positive angle at the pot.

*R₄₄ at 5 × 41 is zero to the measurement*, 0.0055 ± 0.0040. The vertical plane
there is point-to-parallel at the pots and the outgoing vertical angle carries
**no memory of the IP angle at all**. The symplectic check survives it because
R₃₄R₄₃ = −1.009 carries the whole determinant on its own.

*D′ is flat*, 0.0175 / 0.0179 / 0.0182 rad across a factor 6.7 in beam energy,
where R₁₁ grows by 1.6× and R₂₁ falls by 4×.

**The dispersion is quadratic over this range too.** Fitting the `dlt` leg's
nine rigidities with a quadratic rather than a line gives D₂ = −0.195 / −0.196 /
−0.311 m against the carried −0.190 / −0.206 / −0.215. The two lower
configurations close to 3–5%; 18 × 275 does not, and the reason is the window:
the carried fit is three points spanning δ ∈ [−0.143, +0.286] from *three
species* (α, ⁶Li, triton) and is anchored on the triton, while this one is nine
points of one ⁶Li over |δ| ≤ 0.15. They are curvatures over different ranges and
the disagreement is not a correction to either — `POT_DISPERSION_2` is unchanged.

**B0 and ZDC, on the same scan.** The `b0x`/`b0y` legs are the first ladder in
this repository to reach the B0 window (5.5–20 mrad) at all; the routing scan at
the top of this file found `B0 = 0` everywhere because nothing it shot went that
wide. The B0 tracker layers sit at z = **5896 / 6166 / 6436 / 6706 mm**, and in
layer 1

| configuration | dx/dθ_x [m] | dy/dθ_y [m] |
|---|---|---|
| 5 × 41 | 5.9003 ± 0.0007 | 5.9006 ± 0.0010 |
| 10 × 100 | 5.9010 ± 0.0003 | 5.9016 ± 0.0004 |
| 18 × 275 | 5.9001 ± 0.0001 | 5.9013 ± 0.0002 |

— **the lever is the distance**, 5.900 m against the layer's 5.896 m, at every
ring setting to 0.02%. IP6 → B0 layer 1 is a pure drift: a B0 hit measures the
IP angle directly and needs no per-configuration optics, which is why
`farforward.B0_DRIFT_M` is one number where `POT_LEVERS` needs three. The ZDC is
a calorimeter with no momentum branch, so what it contributes is where the
breakup lands: at 18 × 275, 41 events with `HcalFarForwardZDC` hits and 50 with
`EcalFarForwardZDC`, concentrated in the `thx` leg (20 and 25) with ⟨E⟩ = 1.8
and 11.9 GeV at ⟨x⟩ ≈ −1.07 m — the ⁶Li breaking up upstream, not the primary.

**What moved in the code.** Nothing did. `farforward.POT_SECOND_ROW`,
`POT_SECOND_ROW_VERTICAL`, `B0_DRIFT_M`, `B0_LAYER_Z_MM`, `pot_transfer_for`
and `propagate_to_pot` are **new names**; `POT_LEVERS`, `POT_DISPERSION_2` and
the scalar aliases are untouched, and `pot_transfer_for` returns the *carried*
R₁₂/R₃₄/D/D₂ rather than the values re-measured here. Every published figure and
acceptance number is bit-for-bit what it was.

### Caveats (this scan)

One event per scan point, no beam divergence and no vertex spread; hit level,
not reconstruction (the section above shows why reconstruction is not available
for a lithium ion, and §"The Mode-W chain" shows the same matrices returning a
deuteron beam's spectator proton at a 5% momentum). The IP offsets run to ±32 mm, far outside any real beam —
they are a lever for a linear coefficient, not a beam condition, and the
outermost rows drop out of the fit because the displaced ion no longer reaches
station-1 silicon — two of the nineteen horizontal rows at 5 × 41 (|x₀| = 32 mm)
and four at 10 × 100 and 18 × 275, asymmetrically on the −x side (|x₀| ≥ 16 and
≥ 12 mm respectively). R₃₃ and R₄₃
are fitted on 7–10 rows because the vertical top-up rows land outside the
silicon; their standard errors are nonetheless the smallest in the table. Every
element is quoted at **station 1 layer 1** and the station-to-station spread of
the first row (up to 8% in D, 2% in R₁₂) should be assumed for the second.

## The Mode-W chain (2026-09-16) — a weighted BeAGLE sample through EICrecon

The counterpart of the gun scans above on a physics sample, and the first
time a polarized weight of `polligen.reweight` reaches a reconstructed
file.  `../analysis/modew_beagle_hepmc.py` puts a NAMED Mode-W weight on
official BeAGLE `eH2/en/9x130` events and writes HepMC3 ASCII;
`modew_chain.sh` runs npsim → eicrecon on it at a FIXED SEED;
`../analysis/modew_reco_readback.py` reads the weight, the reco-level x and
Q², and the far forward back out.  docs/reproduction_manual.md §5.2 and
§5.3 carry the commands and the expected output.

```bash
singularity exec --env G=$S --env R=$R $SIF bash -lc \
  'source /opt/detector/epic-main/bin/thisepic.sh
   bash $R/tools/fullsim/modew_chain.sh $G/modew_100.hepmc $G/sim 18x275 100'
python3 ../analysis/modew_reco_readback.py $S/sim/reco.edm4eic.root \
    --sim $S/sim/sim.edm4hep.root --csv $S/modew_100.csv
```

Four things it measured, on 100 events, `epic_craterlake_18x275.xml` on
both legs (the nearest configuration to a 130 GeV/u × 9 GeV fill: 6% off in
ion rigidity, a factor two in electron energy — so nothing here is an
acceptance or a resolution).

* **npsim keeps only `weights()[0]`.**  It copies the nominal HepMC3 weight
  into the scalar `EventHeader.weight`, leaves the `EventHeader.weights`
  vector member empty and drops the weight NAMES entirely.  A weight
  appended after the sample's own `default` therefore arrives as a constant
  1.0; written first it arrives on 100/100 events with 100 distinct values,
  to the generator table's own 10-digit printing precision.
* **The event order survives** to a median 5.1 × 10⁻⁷ in x against BeAGLE's
  `trueX` read in file order, so a weight some tool drops can always be
  re-attached by index.
* **Every `InclusiveKinematics*` collection is empty** — Electron, Sigma,
  DA, JB, ESigma, ML all 0/100 — because `MCBeamProtons` is empty: the
  `eH2/en` files record the struck NEUTRON as the status-4 ion beam and
  EICrecon's beam finder wants a proton.  `InclusiveKinematicsTruth` is
  100/100 and `ScatteredElectronsTruth` 100/100, so it is the beam row and
  not the electron.  This is the inclusive-kinematics counterpart of the
  species blocker the section above measures in the far forward, and a ⁶Li
  beam row will fail it the same way.
* **Give npsim a seed.**  Without `--random.seed` it seeds from the clock,
  and the same 100 primaries shot twice gave 335 against 547
  `ForwardRomanPotHits`.  `modew_chain.sh` takes the seed as its fifth
  argument and defaults it to 20260916; `0` restores the clock.  The npsim
  recipes elsewhere in this file predate that and do not fix a seed.

**And the far forward, for a species that IS a proton.**  A deuteron at
130 GeV/u puts its spectator proton at ~128 GeV, R ≈ 0.47 against the
275 GV lattice — the OFF-MOMENTUM detector by the 2026-06-12 routing scan,
and the reconstruction agrees:

| collection | e+d spectator proton | ⁶Li gun ladder (2026-09-15) |
|---|---|---|
| `ForwardOffMTrackerRecHits` | 630 hits, 98 events | — |
| **`ForwardOffMRecParticles`** | **80 in 80 events** | 0 |
| `ForwardRomanPotRecHits` | 420, 37 | 1031, 76 |
| `ForwardRomanPotRecParticles` | 2, 2 | **0, 0** |

Matched event by event to the leading truth forward proton, the 80 OMD
particles come back at a weighted median |p| residual of −5.2% (68th
percentile 5.6%) and +0.05 mrad in polar angle (0.25 mrad).  So the
far-forward reconstruction is **not broken for a light ion — it is broken
for a light ion that is not a proton**: `MatrixTransferStaticConfig.h`'s
`partMass = 0.938272, partCharge = 1` is the right species here, at the
wrong rigidity.  *Labelled inference, untested:* the −5.2% is within 0.1%
of 260.4/275 = 0.9469, the fill's ion rigidity over the lattice's.  One
convention before comparing anything to truth: the tracker-based
far-forward reconstruction (RP, OMD) writes momenta in the ROTATED BEAM
frame, with the crossing angle taken out, while the calorimeter-based
neutrals (ZDC, B0 ECal) write theirs in the LAB frame.

## e+d control inputs (plan 1.5.3)

Official BeAGLE samples streamed via xrootd (no FLUKA needed):
```bash
singularity exec $SIF python3 ../analysis/dump_spectators.py \
  'root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/DIS/BeAGLE1.03.02-1.3/eH2/en/10x130/q2_1to1000/..._run001.hepmc3.tree.root' \
  out.csv --nevents 50000
python3 ../analysis/ed_control_analysis.py out.csv
```
