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
| `~/Projects/eic/local/lib/eic_xl-nightly.sif` | **not present on the analysis box** (checked 2026-08-26). The row is the recipe: EICrecon v1.38.0, detector configs epic-26.03.0…26.06.0 + `epic-main`; source `/opt/detector/epic-main/bin/thisepic.sh` |

`./eic-shell` starts `jug_xl`. Everything below therefore uses the
**September-2024 `epic-main`** far-forward geometry; re-run in the current
container (epic-26.06) before quoting acceptance numbers in the letter.

**Known issue (found 2026-06-12):** the `eic_xl-nightly` pyHepMC3
`rootIO.ReaderRootTree` **segfaults** on files the legacy container reads
fine — use `jug_xl-nightly.sif` for all HepMC3 tree.root *reading*
(dump_spectators.py); npsim/eicrecon in the new container are unaffected.
Consider reporting upstream (eic-shell / pyHepMC3 bindings).

Inside the legacy container, source the geometry with:
```bash
source /opt/detector/epic-24.07.0/setup.sh   # sets $DETECTOR_PATH
```
(current containers use `/opt/detector/epic-main/bin/thisepic.sh`).

## Far-forward gun scan (plan 2.2.1)

```bash
SIF=~/Projects/eic/local/lib/jug_xl-nightly.sif
singularity exec $SIF bash ff_gun_scan.sh ~/Projects/eic/data/ff_gun 60
singularity exec $SIF python3 ff_gun_hits.py ~/Projects/eic/data/ff_gun
```

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

## The intact-⁶Li far-forward scan (2026-08-26) — the pot aperture, measured

> **⚠ The geometry these numbers were measured in has been superseded.**
> Reading the current `main` of `eic/epic` directly (2026-08-26): modules
> went 32 × 32 mm → **16 × 16 mm**; the single energy-independent
> insertion became **per-energy 10σ offsets** in `compact/fields/beamline_*.xml`
> ("These are the ten-sigma cuts for the Roman pots, translated to the
> physical layout we currently have. They are not perfectly ten-sigma for
> reasons of physical geometry."); and the 1 mm aluminium RF shields in
> `roman_pots_eRD24_design.xml` are **commented out** ("we don't know if
> we will even need it — it depends on the beam impedence + RF heating,
> etc. Oct. 2025"). The old 32 mm block gives 32 mm / 30.6 m = 1.046 mrad
> against the 1.03 mrad measured below — agreement to 1.5%, which is the
> check that this scan and the file reading confirm each other. The
> current 16 mm block implies roughly half that, 0.52 mrad, *below* the
> 0.92 mrad 10σ envelope of the Yellow Report high-acceptance optics at
> 18 × 275 (plans/10) and 4× the 0.12 mrad envelope of the lithium tagging
> optics; at 5 × 41 the per-energy insertion moves the other way
> (29.6 mm inner edge). **Re-run this scan in the current container
> before quoting any of the numbers below** — plans/09 B1, priority. What
> a closer approach is *worth* does not depend on it:
> `evgen/scripts/nearbeam_aperture_scan.py` prices every aperture.

The coherent channel's whole acceptance rests on one assumption: where the
Roman Pots stop seeing a nucleus at beam rigidity (plans/04 #20).  It is
now measured rather than assumed.

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

```bash
python3 tools/fullsim/ion_gun_hepmc.py --out li6.hepmc --pdg 1000030060 \
    --a 6 --p-per-nucleon 137.5 --pt 0.78 0.87 1.24 1.65 2.06 2.48 3.30 --nphi 12 > li6.index
singularity exec $SIF bash -lc 'source /opt/detector/epic-main/bin/thisepic.sh
  npsim --compactFile $DETECTOR_PATH/epic_craterlake_18x275.xml --inputFiles li6.hepmc \
        --numberOfEvents 84 --physics.list FTFP_BERT --part.minimalKineticEnergy "100*MeV" \
        --outputFile gun_li6_map.edm4hep.root'
singularity exec $SIF python3 tools/fullsim/ff_gun_hits.py gun_li6_map.edm4hep.root \
        --per-event --index li6.index --positions
```

**Result — the acceptance boundary is HORIZONTAL, and the ⁶Li and the α
agree at the same angle.** Roman-Pot hit / no hit for an intact ⁶Li at
137.5 GeV/u (825 GeV, rigidity ratio 1.000 — beam-blind by construction),
by transverse angle at the IP and azimuth about the ion axis (φ = 0 is
the horizontal bend plane, 90° is up):

| θ [mrad] | 0 | 30 | 60 | 90 | 120 | 150 | 180 | 210 | 240 | 270 | 300 | 330 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.95 | ● | | | | | | | | | | | |
| 1.05 | ● | | | | | | ● | | | | | |
| 1.50 | ● | ● | | | | ● | ● | ● | | | | ● |
| 2.00 | ● | ● | | ● | ● | ● | ● | ● | | | ● | ● |
| 2.50 | ● | ● | ● | ● | ● | ● | ● | ● | ● | | ● | ● |
| 3.00 | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● | ● |
| 4.00 | | | ● | ● | ● | ● | ● | ● | ● | ● | | ● |

Reading it: the boundary is **|θ_x| ≳ 1.0 mrad**, or |θ_y| ≳ 1.8 mrad up
and ≳ 2.8 mrad down; the last row is the far edge, where a horizontal
track overshoots the outermost module (|Δx| ≈ 112 mm) and is lost again,
so the horizontal window is roughly 1.0 < θ_x < 3.5 mrad.  A separate
1-D scan puts the horizontal edge between 0.970 and 1.091 mrad for the
⁶Li and between 0.909 and 1.091 mrad for an α at the same rigidity
(550 GeV) — the same *angle*, which is the check that the scan is
measuring rigidity-driven optics and not the species.

**Why horizontal, when the geometry inserts the pots vertically.**
`compact/far_forward/roman_pots_eRD24_design.xml` builds top and bottom
assemblies at y = ±insertion, each 2×2-module array spanning 64 mm from
its centre, so the silicon covers |Δx| ≤ 112 mm continuously and reaches
down to y = 0 except under the central |Δx| < 32 mm block, which the
high-acceptance insertion holds off to |y| ≥ 7 mm.  The open aperture is
therefore a **horizontal slot ≈ 64 mm wide and ≈ 14 mm tall** at the pot
plane — and the reason a 2 mrad *vertical* kick still fits through it is
the optics.  Reading the first-crossing hit positions at station 1
(z = 32.55 m) straight off the scan: a horizontal IP angle arrives with a
lever of **R₁₂ ≈ 30.6 m** (Δx = 77.1 mm at θ = 2.5 mrad, 91.5 mm at 3.0,
−112 mm at 2.5 on the other side), while a vertical one arrives at only
**y = 8.0 / 8.6 / 9.0 / 10.6 mm for θ_y = 2.0 / 2.5 / 3.0 / 4.0 mrad**.
The far-forward transport is about ten times stiffer in x than in y, and
that, not the pot orientation, is what sets the near-beam boundary.  The
same reading gives the outer edge: horizontal coverage runs to Δx ≈ −120
and +112 mm, which is why θ = 4 mrad at φ = 0 is lost again.

The two readings close on each other, which is the check that this is a
measurement and not an artefact: the slot's half-widths are 32 mm in x
and 7 mm in y, and the observed edges are θ_x ≈ 1.03 mrad → Δx = 31 mm
with R₁₂ = 30.6 m, and θ_y ≈ 2.0 mrad → y = 8.0 mm, the first value above
7 mm.  Aperture and optics reproduce each other to a millimetre.

**Each energy configuration has its own optics, and its own edge.**  The
scan was repeated with `epic_craterlake_5x41.xml` and
`epic_craterlake_10x100.xml`, each with the ⁶Li momentum that puts it at
that ring's reference rigidity (123 and 300 GeV):

| configuration | ⁶Li p_u [GeV/u] | p [GeV] | horizontal edge | as p_T |
|---|---|---|---|---|
| 5 × 41 | 20.5 | 123 | ≈ 2.0 mrad | 0.25 GeV |
| 10 × 100 | 50 | 300 | ≈ 1.35 mrad | 0.41 GeV |
| 18 × 275 | 137.5 | 825 | ≈ 1.03 mrad | 0.85 GeV |

The angle tightens with energy but not as fast as the momentum grows, so
in p_T — which is what the physics delivers — the low-energy
configuration is 3.4× more permissive.  Against the fast simulation's
angular envelope (10 σ_θ, p_T = 0.089 / 0.218 / 0.60 GeV for high
acceptance and 0.18 / 0.45 / 1.23 for high divergence) the measured
aperture is 2.8× / 1.9× / 1.4× the high-acceptance number, i.e. it sits
between the two optics assumptions and closer to high divergence.

**What it changes.**  `polligen.reco.rp_measure` models the cutout as a
slot *wide in x and tight in y* — its docstring says so and
`cut_scale_xy = (2.5, 1.0)` encodes it, i.e. acceptance opened by the
VERTICAL displacement.  The measurement inverts that: (c_x, c_y) ≈ (1.0,
2.3) mrad at 18 × 275, an aspect the wrong way round by a factor ≈ 5.8.
Two consequences, both computed with the repository's own
`reco.rp_hole_acceptance` at B = 50 GeV⁻² with the per-optics edges above
(c_y taken as 3.0 mrad at 5 × 41 and 10 × 100, where the 30° scan only
bounds it between 2 and 3):

| p_u [GeV/u] | tagged fraction, measured aperture | assumed HA slot | ⟨cos 2φ_t⟩ measured | assumed |
|---|---|---|---|---|
| 20.5 | 1.4×10⁻² | 0.824 | **+0.78** | −0.09 |
| 50 | 5.1×10⁻⁵ | 0.401 | **+0.90** | −0.42 |
| 137.5 | 1.9×10⁻¹⁷ | 1.6×10⁻² | **+0.97** | −0.80 |

(Both columns are `rp_hole_acceptance`'s own convention, φ measured from
**x**; the analysis's β is measured from the *vertical* spin axis, so
⟨cos 2β⟩ is the negative of each — the assumed slot gives the "large and
positive" ⟨cos 2β⟩ that plans/04 #20 records, and the measured aperture
gives a large NEGATIVE one.)

The coherent tag is far more restrictive than assumed and survives only
at the low-energy configuration, and the acceptance-induced fake cos 2φ_t
**changes sign** — which is exactly the harmonic the coherent fit
templates (report §7, `reco.basis_2d`).  The magnitude is what the
template absorbs; the sign is what a mis-specified template gets wrong.

**Carried through the coherent chain.**  `reco.RP_APERTURE_MEASURED` is
the table above, `rp_measure(cut_theta_xy=…)` takes the larger of the
envelope and the aperture per axis (they are separate constraints), and
`money_cos2phi_coherent_reco.py --rp-aperture measured` runs 6R on it.
At the **low** configuration the measurement survives: acceptance
37.7% → 1.42%, N_tag 8.3×10⁶ → 3.1×10⁵, the acceptance-induced
⟨cos 2β⟩ **+0.426 → −0.772** — the sign flip, in the chain rather than in
a formula — two of the four |t| bins instead of four, δa_t worse by
6–34×, and a_e still recovered (0.0073 ± 0.0045 and 0.0091 ± 0.0045
against an injected 0.0100).  At **mid** and **top** the aperture leaves
no accepted recoil in the binned |t| window at all.  The coherent
programme is therefore a low-energy programme for a second and stronger
reason than the angular envelope already gave.

**Caveats before this is quoted as a correction.** One event per scan
point, 30° azimuthal steps, no beam divergence or vertex spread, and the
`jug_xl-nightly` geometry is `epic-main` at git 5a7dd057 (≈ September
2024) — the same file this README already flags as having moved the pot
stations from the published 26/28 m to 32.5/34.3 m.  The high-divergence
insertion is present in the XML but commented out, so the scan measures
the **high-acceptance** configuration only.  Repeat on the current
release, with a beam envelope, before the letter quotes it; the item is
raised in plans/04 #20.

## e+d control inputs (plan 1.5.3)

Official BeAGLE samples streamed via xrootd (no FLUKA needed):
```bash
singularity exec $SIF python3 ../analysis/dump_spectators.py \
  'root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/DIS/BeAGLE1.03.02-1.3/eH2/en/10x130/q2_1to1000/..._run001.hepmc3.tree.root' \
  out.csv --nevents 50000
python3 ../analysis/ed_control_analysis.py out.csv
```
