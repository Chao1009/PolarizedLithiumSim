# Phase 2 — Full EIC (ePIC) Simulation

**Goal.** End-to-end demonstration with the real detector: generator events →
beam-effects afterburner → Geant4 (npsim + ePIC DD4hep geometry) →
reconstruction (EICrecon) → physics analysis. Two questions only full sim
can answer: (1) far-forward acceptance/efficiency/PID for Li fragments with
real magnets, beam pipe, and optics; (2) reconstructed-level asymmetry
extraction with realistic resolutions, acceptances, and backgrounds
(closure test of the Phase-1 FOM).

*Tool-chain facts below were fetch-verified against eic.github.io, the
eic/* GitHub repos, and the cited arXiv papers on 2026-06-12.*

**Status legend:** ☐ todo ◐ started ☑ done

---

## Step 2.0 ☑ Local environment refresh (days)

Current local state (checked 2026-06-12): eic-shell at `~/Projects/eic`
with `jug_xl-nightly.sif` (4 GB, ~Sep 2024); `npsim`/`eicrecon` confirmed
runnable inside. The container has since been **renamed `eic_xl`**; current
releases: epic geometry **26.06.0**, EICrecon **v1.38.0**, npsim **v1.6.1**.
Local `epic` @ 24.08.0 and `EICrecon` @ v1.6.0 checkouts are ~2 years stale.

1. Reinstall fresh (old installer pulls the renamed image):
   ```bash
   cd ~/Projects/eic && curl --location https://get.epic-eic.org | bash
   ./eic-shell        # pin a release: curl -L https://get.epic-eic.org | bash -s -- -v 26.06-stable
   ```
2. Use the container-shipped geometry; source it the current way:
   ```bash
   source /opt/detector/epic-main/bin/thisepic.sh   # sets $DETECTOR_PATH, $DETECTOR_CONFIG
   ```
   `epic_craterlake.xml` is the flagship config; energy/species variants
   exist as `epic_craterlake_{5x41,10x100,18x275}.xml` and ion configs
   (`_18x110_Au`, `_10x110_He3`, `_10x130_H2`, …). Far-forward detectors
   are in craterlake by default. Build `epic`/`EICrecon` from source only
   when modifying geometry (needed in step 2.1 for Li beamline files).
   ☑ *2026-08-28: done in `tools/fullsim/ff_gun_scan.sh` — it sources
   `thisepic.sh` with fallbacks and selects `epic_craterlake_18x275.xml`
   explicitly (plain `epic_craterlake.xml` loads the 5×41 beamline fields);
   every `--compactFile` in the tree is `$DETECTOR_PATH`-relative, and the
   recipe is in `tools/fullsim/README.md` and docs/reproduction_manual.md §5.3.
   The build-from-source clause stays with step 2.1.*
3. Smoke test per tutorials (eic.github.io/documentation/tutorials.html;
   esp. tutorial-simulations-using-npsim-and-geant4, tutorial-analysis,
   tutorial-jana2):
   ```bash
   npsim --compactFile $DETECTOR_PATH/epic_craterlake.xml -N 100 \
         --inputFiles input.hepmc --outputFile sim.edm4hep.root
   eicrecon -Ppodio:output_file=reco.edm4eic.root sim.edm4hep.root
   ```
   Record exact commands + versions in `fullsim/README.md` as the runbook.
   *(2026-09-15: the block above is the tutorial's, and the runbook's form
   differs in two places the first real run found — `epic_craterlake.xml`
   loads the 5 × 41 beamline fields, and `eicrecon` does NOT inherit npsim's
   compact file because `thisepic.sh` exports `DETECTOR_CONFIG=epic`, so both
   legs need the configuration named explicitly, `--compactFile
   $DETECTOR_PATH/epic_craterlake_18x275.xml` and
   `-Pdd4hep:xml_files=$DETECTOR_PATH/epic_craterlake_18x275.xml`.)*

☑ *2026-09-06: the step header is set to ☑. Item 1 is done — the current
`eic_xl-nightly` container (`epic-main` 9aaa2969, 2026-08-22) is installed at
`~/Projects/eic-2026/` (run 14); item 2 carries its own ☑ above; item 3's runbook
is written, in `tools/fullsim/README.md` and `docs/reproduction_manual.md` §5.3,
and npsim is exercised against the container geometry by
`tools/fullsim/ff_gun_scan.sh`. The one clause never exercised is the literal
`eicrecon` reconstruction command — no EICrecon output exists anywhere in the
tree — and it closes with the first reconstructed run under step 2.3, not here.*

## Step 2.1 ☐ Beam configuration for Li species (2–3 weeks, w/ FF group)

Three Li-specific customizations are required (none exist today — verified
by absence in eic/afterburner, eic/epic, eic/BeAGLEsamples):

1. **Afterburner preset.** `abconv` applies the 25 mrad (+100 µrad vertical)
   crossing angle, divergence, crab kick, and vertex smearing *before*
   npsim (campaign convention; never combine with ddsim's
   `--crossingAngleBoost`). It requires exactly two status-4 beam
   particles (one PDG 11) and snaps the ion energy to brackets
   {275,250,166,130,115,110,100,41}×{18,10,9,5} GeV — ⁶Li @137.5 or
   ⁷Li @117.9 GeV/u **fall in no bracket** → abconv throws. Options:
   (a) run 110 GeV/u points with `-p ip6_eAu_110x10`-style presets as an
   optics proxy (β*/divergence will be Au's, not Li's — fine for first
   acceptance look); (b) add `ip6_eLi6_*`/`ip6_eLi7_*` preset functions in
   `cpp/afterburner/EicConfigurator.cc` (copy the eAu/eHe3 pattern) with
   C-AD-blessed parameters, and PR upstream.
2. **Beamline field maps in eic/epic.** Species configs differ only in
   `compact/fields/beamline_{E}x{P}_{species}.xml` (rigidity-scaled
   magnets) + a `configurations/craterlake_*.yml`. For ⁶Li (Z/A = 1/2) the
   existing `beamline_5x41_He4.xml` is a direct template (same Z/A); ⁷Li
   (Z/A = 3/7) needs new scale factors. Coordinate with the ePIC
   far-forward WG (A. Jentsch) — the Au files carry explicit "ASK FF
   EXPERTS" warnings.
3. **Generator runcard + conversion.** BeAGLE accepts A,Z via runcard
   (`eic/BeAGLEsamples` has d and ³He cards; write the Li ones). Chain:
   BeAGLE output → eic-smear `TreeToHepMC` → `abconv` → npsim. HepMC3
   conventions: status 1 = transported, 4 = beam (kept in
   `MCParticles.generatorStatus` but not transported); units read from the
   file header; fragments as final-state particles with 10-digit nuclear
   PDG codes (α = 1000020040, t = 1000010030). Watch-item: excited-ion
   codes (10LZZZAAAI, I≠0) from BeAGLE had DD4hep fixes — sanitize to
   ground-state codes if npsim drops them (`npsim/scripts/
   sanitize_hepmc3.py` exists for header fixes).
4. Validate: afterburned vertex/angle distributions reproduce the preset
   tables; α/t/d primaries survive HepMC3 → npsim → Geant4 ion transport
   (100-event gun jobs).

## Step 2.2 ◐ Far-forward acceptance for Li fragments (4–8 weeks, core)

**This is the novel deliverable.** Verified far-forward suite (YR matrix;
arXiv:2108.08314 Table I; arXiv:2409.02811; arXiv:2406.12877):

| detector | z | θ window | rigidity/x_L window | key resolution |
|---|---|---|---|---|
| B0 tracker+EMCal | 5.4–6.4 m | 5.5–20 mrad | any charged in window | δp/p ≈ 2–4%; γ: 6–7%/√E, E>50 MeV |
| Roman Pots ×2 | 26 & 28 m (10σ from beam) | 0–5 mrad | R ≈ 0.6–0.95 | pT cutoff ≈ 0.2 (high-acc optics) – 0.45 (high-div) GeV/c |
| Off-Momentum ×2 | 22.5 & 24.5 m | 0–5 (uniform <2) mrad | R ≈ 0.4–0.6 | σ_pT ≈ 20% @100 MeV/c; no low-pT cutoff |
| ZDC | 35–37.5 m | 0–4/4.5 mrad | neutrals | n: 50%/√E ⊕ 5%; γ crystal: 5%/√E ⊕ 3% |

(θ ≈ 5–5.5 mrad RP/OMD↔B0 gap is a known hole; beam-pipe material costs
5–20%.)

— *superseded (2026-08-28): the z positions and the p_T cutoffs of this table
have been overtaken by the geometry read and by plans/10. `farforward.py` (from
current `eic/epic` main) puts the Roman Pots at 32.55/34.25 m and the OMD at
25.50/27.00 m with R ∈ [0.45, 0.65], and the p_T cutoffs — a 275 GeV proton
number — are replaced by an angular 10(σ_h, σ_v) envelope that varies per
configuration (`farforward.sigma_theta_for`, Report 3 Table 7). Rewrite against
the code when the preTDR layout is confirmed (`tools/fullsim/README.md`).*

**Rigidity routing, R = (A_f·Z_beam)/(A_beam·Z_f)** — corrected mapping
that drives the whole tagging program:

| fragment | ⁶Li beam (Z/A=1/2) | ⁷Li beam (Z/A=3/7) |
|---|---|---|
| n | ZDC | ZDC |
| p | R=0.50 → OMD center | R=0.43 → OMD low edge |
| d | R=1.00 → **beam-blind** below RP pT cutoff | R=0.86 → Roman Pots |
| t | R=1.50 → **lost** (152 mm at the pot plane, past the 144 mm module edge) | R=1.29 → **RP-inner + ZDC** (measured 2026-08-28, dx = +66 mm) |
| ³He | R=0.75 → Roman Pots | R=0.64 → RP/OMD boundary |
| α | R=1.00 → **beam-blind** (pT tail only) | R=0.86 → **Roman Pots** |

— *superseded (2026-08-28): the table's own formula was replaced by
`spectator.py`'s mass-based R = (m_spec/Z_spec)/(m_beam/Z_beam) (plans/08 C1),
which separates the ⁶Li α (0.99813, under the orbit) from the d (1.00452, over
it) where A·Z arithmetic puts both at exactly 1; and the ⁷Li proton cell is
falsified by the corrected OMD window R ∈ [0.45, 0.65], which loses R = 0.43
(`route_charged`, `test_spectator.py`). The triton rows are the ◐ item below.*

Consequences to quantify (no literature exists for e+Li breakup tagging —
the only Li far-forward datapoint is coherent J/ψ at IR-8, arXiv:2511.05638,
⁷Li tagging eff. 17.75%):
- **⁷Li α-tag works at IP6** (R=0.86 sits mid-RP-window): e+⁷Li → e′+X+α
  selects DIS on the *triton cluster* — the polarized-EMC companion channel.
  Quantify acceptance × purity vs x_L, pT.
- **⁶Li α-tag is the hard case** (R=1.0): only the pT ≳ 0.2–0.45 GeV/c tail
  is visible at IP6, but the α–d cluster Fermi momentum is soft → fold the
  BeAGLE/cluster-model pT(α) distribution with the 10σ optics cutoff for
  both optics settings; this number decides whether the d-cluster tagging
  argument holds at IP6, needs the IR-8 secondary focus (RPs at 44–45.5 m
  recover R≈1 down to pT ~ 0), or relies on the ⁶Li → p/³He channels.
☑ *2026-08-28: the fold is done — `fastsim/scripts/tagging_acceptance.py`
folds the cluster-model p_T(α) against a 10σ angular envelope for four optics
per configuration, giving a ⁶Li α tag of 1.9 / 1.7 / 2.6% at the Yellow Report
optics and 31 / 22 / 29% at a lithium tagging optics costing 1/6.8–1/12.8 of the
luminosity (35 / 28 / 29% until the per-configuration transport became the
baseline on 2026-08-29) (`fastsim/out/tagging_acceptance.txt`, Report 3 Table 6). The
argument therefore does not hold at any published IP6 optics: it needs the
tagging optics or IR-8. The third branch, the ⁶Li → p/³He channels, was never
quantified as a tag — those fragments are treated as breakup backgrounds.*
- **Tritons at IP6 — revisit!** ☑ *2026-08-28 (tools/fullsim/README.md,
  plans/09 B1): in the current epic-main geometry (git 9aaa2969) the
  over-rigid triton (R = 1.286) crosses the Roman-Pot planes at
  dx = +66 mm at station 1 and +70 to +72 mm at station 2 on the inner
  (over-rigid) side of the beam orbit — 60 of 60 events at all three
  configurations, in the 48–144 mm outer band, which carries no vertical
  insertion, so the hit is at y ≈ 0 and needs no pT — and then deposits
  in the ZDC in 80 / 83 / 98% of the same events. The ~36 mm of the
  2026-06-12 scan was the September-2024 geometry. The "no coverage"
  assumption is withdrawn: `farforward.route_charged` has an over-rigid
  branch (code 6, RP-inner) and the ⁷Li t tag is 78 / 92 / 94%. What
  remains is reconstruction: the fast simulation routes in angle and
  cannot see the acceptance hole between θx = −1.55 and −2.53 mrad that
  the triton's own 8% larger R12 opens against the 16 mm blind block.*
- ZDC neutron tagging (evaporation n from both isotopes) + de-excitation γ.

Work plan: (1) particle-gun (x_L, pT, vertex) scans per fragment per beam
setting → acceptance maps; validate against published e+d maps
(arXiv:2108.08314 App. B) as control. (2) Check EICrecon FF reconstruction
(`ForwardRomanPotRecParticles`, `ForwardOffMRecParticles`, B0, ZDC
collections) for Li-rigidity transfer matrices; request/derive if
proton-tuned. (3) Fold BeAGLE e+Li events → tagging efficiency × purity ×
mis-tag matrix at reconstructed level. (4) Publish "Li far-forward tagging
performance" note; feed parameterizations back into the Phase-1 fast-sim.

☐ *2026-08-28 on (2): only the RP/OMD path was checked, and only for the beam-PDG
matrix selection (Report 3 item 13). B0 and ZDC were never examined, and
`MatrixTransferStaticConfig.h` hardcodes `partMass = 0.938272`, `partCharge = 1`
with no override in `RPOTS.cc`/`FOFFMTRK.cc` — every far-forward particle is
reconstructed as a proton. The "derive" half is ours, not the FF group's, and
the repo holds R12, R34 and D per configuration (`farforward.POT_LEVERS`:
19.24 / 4.56 / 0.311, 21.25 / 3.35 / 0.287 and 29.97 / 2.93 / 0.292 m, the 5 × 41 vertical lever read off a zero-insertion scratch geometry on 2026-08-29) plus the
second-order dispersion D2 = −0.190 / −0.206 / −0.215 m, which is what makes
the over-rigid arithmetic reproduce the measurement to a tenth of a
millimetre; R11, R21, R22 and D′ are still unmeasured.*

◐ *2026-09-15 on (2): the "derive" half is done and the B0/ZDC half is examined.
R11, R21, R22 and D′ are measured, per configuration, at pot station 1, by
`tools/fullsim/ff_transfer_scan.py` through the zero-insertion geometry of
`ff_zero_insertion.sh`: 1.148 / 1.227 / 1.852, −0.0837 / −0.0651 / −0.0209 rad m⁻¹,
−0.4955 / −0.3060 / +0.1944 and 0.0175 / 0.0179 / 0.0182 rad at 5 × 41 / 10 × 100 /
18 × 275, with 0.8–4.5 % standard errors, plus the vertical R33/R43/R44. They sit
in `farforward.POT_SECOND_ROW` and `POT_SECOND_ROW_VERTICAL` as NEW names —
POT_LEVERS and D2 are untouched and came out of the same files as controls (R12
19.186 / 21.341 / 29.961 m against the carried 19.24 / 21.25 / 29.97). The blocks
are symplectic to 4 % without having been fitted to be. B0: the 6–20 mrad leg is
the first ladder here to reach the B0 window, and IP6 → B0 layer 1 is a PURE DRIFT
— dx/dθ_x = 5.900 m at all three settings against the layer's z = 5896 mm, so a B0
hit measures the IP angle with no per-configuration optics (`farforward.B0_DRIFT_M`).
ZDC is a calorimeter with no momentum branch and is reported as hit counts and
energy-weighted centroids per leg. The "request" half is now a MEASURED external
blocker rather than a code reading: with the first EICrecon run (Step 2.3,
2026-09-15) the Roman-Pot reconstruction returns zero `RecParticles` for a lithium
ion where a proton control gives 72, so `MatrixTransferStaticConfig.h`'s
`partMass = 0.938272, partCharge = 1` does not merely mislabel the species, it
drops it. Item (3), the reco-level efficiency × purity × mis-tag matrix, is
unchanged and still externally gated.*

☐ *2026-08-28 on (4): feeding parameterizations back is done (the measured pot
aperture and the per-configuration optics are in the fast sim), and the
acceptance half of the note is published (Report 3 §5/Table 6, Report 4). What
remains is the reco-level efficiency × purity × mis-tag matrix of work-plan item
(3), which is externally gated — it needs BeAGLE A = 6,7 (FLUKA licence) and
step 2.1's Li presets.*

◐ *2026-09-16 (author decision A1, on (2) and on the step box): the step box
moves ☐ → ◐. The transfer half of work-plan item (2) is measured rather than
requested — R₁₁ = 1.148 / 1.227 / 1.852, R₂₁ = −0.0837 / −0.0651 / −0.0209
rad m⁻¹, R₂₂ = −0.4955 / −0.3060 / +0.1944 and D′ = 0.0175 / 0.0179 / 0.0182 rad
at 5 × 41 / 10 × 100 / 18 × 275, each to 0.8–4.5%, with the vertical block and
the B0 pure drift beside them (the 2026-09-15 note above) — so this step now
holds the whole first-order IP6 → pot transfer per configuration. What keeps it
short of ☑ is named and unchanged: EICrecon reconstructs no lithium ion at all
(zero `ForwardRomanPotRecParticles` where a 275 GeV proton control gives 72,
Step 2.3's 2026-09-15 note), and the ⁷Li beamline field maps of Step 2.1 item 2
do not exist, so item (3)'s reco-level efficiency × purity × mis-tag matrix
cannot be built here whatever this repository measures. Items (1) and (4) are
untouched by this decision.*

## Step 2.3 ☐ Central-detector physics performance (4 weeks)

☐ *2026-09-15: the residue Step 2.0 deferred here is closed — the literal
`eicrecon` command has been run and the tree has its first reconstructed sample.
100 ⁶Li events at 137.5 GeV/u through `epic_craterlake_18x275.xml` gave
`sim.edm4hep.root` (62.6 MB) and `reco.edm4eic.root` (37.5 MB, 100 events, 1335
collections); npsim 1.8.0, `EICrecon git.77d0cef8…=main`, epic `9aaa2969…`, 6 min
10 s + 47 s on the analysis box. Commands and versions are in
`tools/fullsim/README.md` §"The reconstruction leg" and docs/reproduction_manual.md
§5.3. Two things it taught. `thisepic.sh` exports `DETECTOR_CONFIG=epic`, so the
command as written loads `epic.xml`, which includes `beamline_5x41.xml` — the
`ff_gun_scan.sh` gotcha on the reconstruction side; pass
`-Pdd4hep:xml_files=$DETECTOR_PATH/epic_craterlake_18x275.xml`. And the
far-forward reconstruction returns NOTHING for a lithium ion: 1031 of 1179
Roman-Pot sim hits become `RecHits` and zero become `ForwardRomanPotRecParticles`,
where a 275 GeV proton control through the same chain gives 72 at 275.0 ± 1.7 GeV.
The central-detector items 1–4 below are untouched by this and stay ☐.*

1. Scattered electron: e-ID efficiency/purity vs (x,Q²); kinematic
   reconstruction (electron vs JB vs Σ/DA — tutorial-kinematic-
   reconstruction) on e+Li with Fermi smearing.
2. x–Q² migration matrices at reco level → re-derive per-bin FOM; confirm
   Phase-1 binning (purity ≳ 0.8) or rebin.
   — *superseded (2026-09-16, author decision A8): the purity ≳ 0.8 gate is
   retired with plans/02 Step 1.6's; the requirement in its place is that the
   folded shape fit be unbiased at the measured purity. See the ☑ note below.*
3. φ resolution & acceptance × crossing-angle correlations vs the cos 2φ
   gluonometry amplitude: inject known Δ-modulation by weighting, fit at
   reco level, quantify dilution + fake modulation. Make-or-break for the
   gluonometry case — do early with small samples.

☐ *2026-08-28: the parametric chain of plans/02 §1.6 and plans/07 WP3 answers the
questions of items 2 and 3 — migration matrices, per-bin purity/efficiency,
the φ′ dilution and the fake cos 2φ′ from a split acceptance — but it is the
Phase-1 FOM this step exists to close against, not this step. No central-detector
npsim/EICrecon run exists anywhere in the tree. Item 1 is additionally short of
its own words: `reco.eps_eid(η)` carries no (x, Q²) dependence and no purity or
background model, and there is no Fermi smearing in the reconstruction chain.
Item 2's "purity ≳ 0.8 or rebin" gate is unmet even parametrically (0.56–0.75
calibrated, 0.42–0.68 uncalibrated) — see the note under plans/02 Step 1.6.*

☐ *2026-09-15: item 1's Fermi half is now in the chain and its e-ID half is
closed as far as the reference set allows. `RecoModel.fermi_smear` (default
off) gives every pseudo-event's struck cluster an internal momentum from the
cluster momentum densities of `polli_fastsim.spectator`, so the vertex sits at
ss = α s (`reco.fermi_alpha`) and the analysis, reconstructing with the nominal
per-nucleon beam, reports x_meas = α x_vertex at the same Q²; ⟨|k|⟩ = 0.1071 GeV
and ⟨k²⟩^½ = 0.1345 GeV reproduce the density's own moments to 0.04% and 0.02%,
α = 0.9985 ± 0.0297, and the effect on the four 5R sweet spots is
+0.21 / +0.07 / −0.06 / −0.01% against statistical errors of 1.7–3.2%, with the
bin-centering ratio moving by 0.01% and 7R unmoved at the printed precision
(`money_cos2phi_reco.py --fermi-smear`, own stem). The (x, Q²) axis on ε_eID is
NOT closed and cannot be from the Yellow Report: its only electron efficiency is
the working point its pion-suppression numbers are quoted at — 95% for an E/p
cut, 92% with the shower shape (Table 11.29, Fig. 11.48) — flat in η, momentum,
x and Q²; Fig. 11.49 is purity, Fig. 8.4 is yields. The curve stays η-only, the
constant is recorded as `reco.EPS_EID_YR_WORKING_POINT`, and the gap is written
into the docstring, Report 2 Table 2 and Report 2 §7 rather than filled with
invented numbers. The purity ≳ 0.8 half of item 2 is untouched; so is the npsim/
EICrecon gate.*

☑ *2026-09-16 (author decision A8, on item 2): the "purity ≳ 0.8 or rebin" gate
is RETIRED here as it is in plans/02 Step 1.6 — superseded, not met. It was
written for a bin-by-bin correction, where a bin whose reconstructed events are
mostly migrants cannot be corrected by its own factor; the chain no longer works
that way. `recopseudo.fold_shape_fit` (plans/08 A6) fits the Δ(x) shape through
the response itself, so migration is folded rather than divided out and the
per-bin purity enters as a statistical cost, not as a bias. The MEASURED
purities stand as the state of the chain and are not a failure of it: 0.56–0.75
calibrated and 0.42–0.68 uncalibrated at the published spots, 0.68–0.70 at the
low configuration, and 0.73 / 0.80 / 0.76 / 0.77 for the electron method alone
above y = 0.05 against 0.65 / 0.64 / 0.70 / 0.69 for the mixed method in the same
bins. The requirement in their place is that THE FOLDED FIT IS UNBIASED AT THE
MEASURED PURITY, which the closure tests already pin:
`evgen/tests/test_recopseudo.py::test_folded_fit_recovers_the_truth_from_a_wrong_prior`
(a prior with the wrong F₁ factor biases the bin-by-bin K by −22%/+16% and the
folded fit by −6.3%/+3.5%, ~3.5× better on the worst bin),
`::test_folded_fit_is_the_bin_by_bin_K_when_the_prior_is_the_truth` (χ² < 10⁻¹²
and zero tilt when the prior is the truth, so opting in moves nothing),
`::test_prior_spread_covers_the_residual_the_tilt_cannot_absorb`,
`::test_folded_fit_does_not_oscillate` and
`::test_fold_reproduces_the_reco_bin_amplitude_of_any_model`; and the
pseudo-experiment pulls are unbiased in the runs themselves — the twenty pulls
of the electron-method window, five response seeds at purities 0.73–0.80, have
mean +0.10 and sd 0.99. Nothing is rebinned and no published number moves. The npsim/EICrecon gate on this step
is untouched.*

## Step 2.4 ☐ Pseudo-experiment closure tests (4 weeks)

1. Weight reconstructed unpolarized samples with Phase-1 asymmetry models
   (A∥, A_zz, A_cos2φ) incl. spin-state bookkeeping (bunch patterns;
   λ = +,0,− luminosity shares for tensor running).
2. Extract observables with the analysis estimators; verify unbiased pulls;
   quote stat ⊕ reco systematics.
3. Systematics: relative luminosity between spin states, δP/P polarimetry
   (HJET-for-Li is itself R&D), φ-acceptance stability.

☐ *2026-08-28: the estimators, the spin-state bookkeeping (`bookkeeping.py` run
plans, relative-luminosity offsets, δP/P smearing) and a parametric reco-level
closure all exist, but they are the plans/05 §5.A and plans/07 WP3 deliverables:
no reconstructed sample has ever come out of EICrecon (`eicrecon` is invoked in
zero scripts), and the generator has no HepMC3 writer (plans/05 step 5.D). Two
coverage gaps survive even at the parametric level: only A_cos2φ has been run
through a detector response — there is no reco-level A∥ or A_zz — and
`RunPlan.delta_p_over_p` is read by nothing outside the tests, so no money plot
yet carries a polarimetry-scale band.*

☐ *2026-09-15: both coverage gaps of that note are closed at the parametric
level; the EICrecon gate is not. `recopseudo.measure_azz` / `measure_apar` run
the tensor-thirds and helicity-flip estimators of `polligen.estimators` on the
exact expected per-fill counts of a RECONSTRUCTED (x, Q²) bin of the same
response 5R and 7R use — `RecoResponse.expected_rates` — with the analytic
errors of `asymmetries.err_azz` / `err_a_parallel`. The thirds estimator is an
identity there and is pinned as one: across the three fills the φ-averaged rate
factors sum to exactly 3, so read at true kinematics with neither selection nor
ε_eID it returns the σ-weighted analytic A_zz to 2×10⁻¹⁴, pinned at 1e-9 (A∥
to 0.29%, the residue
being the finite-γ longitudinal form the kernel carries by default). Over 240
pseudo-experiments of a bin at x = 0.03–0.10, Q² = 3–10 GeV², the pulls against
the reconstructed-bin truth have means −0.056 ± 0.066 (A_zz) and +0.103 ± 0.059
(A∥) with spreads 1.03 and 0.91 of the analytic errors, and the
relative-luminosity biases of item 3 are reproduced through the response and
removed by the luminosity-corrected estimators. `RunPlan.delta_p_over_p` is now
read by a money plot: `money_azz_reco.py --delta-p-over-p` draws it as the pure
scale it is — both estimators divide by a MEASURED polarization, and the band
is measured rather than assumed: each edge is the same drawn counts re-read by
the same estimator at P(1 ± δP/P), whose half-width comes out at
δP/P/(1 − (δP/P)²) of the central value at every bin, δP/P to 0.09%, and one
polarization draw scales all seven bins by the same 1/(1 + δP/P) to 2×10⁻¹⁶ — on
its own stem `money_azz_reco_pol0p03_6Li.png`. At the 3% ring value that band is
0.18–0.53 of the one-year statistical error and 0.56–1.67 of the ten-year one,
so polarimetry and not statistics limits a ten-year A_zz over the middle of the
Q² = 3–10 GeV² slice. Item 3's φ-acceptance stability stays where plans/07 WP3
left it.*

☐ *2026-09-15 on the EICrecon clause above: `eicrecon` now runs in the tree
and a reconstructed sample exists (Step 2.3, 100 ⁶Li events at 137.5 GeV/u;
`tools/fullsim/README.md` §"The reconstruction leg"), so "invoked in zero
scripts" no longer holds. What this step needs is unchanged: that sample is a
far-forward gun ladder, not a weighted physics sample, and nothing has been fed
back to the estimators. The two parametric coverage gaps of the 2026-08-28 note
are closed by the clause above, at the parametric level; the reconstructed
sample this step is written around is what stands.*

## Step 2.5 ☐ Campaign-scale production & write-up (ongoing)

1. Existing eA campaign datasets to reuse for technique validation
   (verified on `root://dtn-eic.jlab.org//volatile/eic/EPIC/`, browsable
   with `xrdfs ls`): e+d BeAGLE 10×130 (tagging samples), e+³He BeAGLE
   5×41/10×110/18×110/10×166 (A₁ⁿ double-tagging), e+Au BeAGLE DIS.
   Campaign infra: `eic/simulation_campaign_hepmc3` + condor; the beam pair
   and the species enter via the geometry filename
   (`${DETECTOR_CONFIG}_${EBEAM}x${PBEAM}.xml`, the species riding inside
   `PBEAM` — `100_Au197`, `110_He3`, `130_H2` — which is an inference from the
   compact and RECO listings, not from a submit file; see the note below),
   while the dataset path itself mirrors the generator sample's own path
   (`RECO/<version>/<DETECTOR_CONFIG>/<class>/…/<generator>/<species>/<EBEAMxPBEAM>/<Q2 bin>/`,
   with the middle depth varying by dataset);
   afterburning is done upstream of campaigns (`*_ABCONV` datasets).
☐ *2026-08-28: only the EVGEN e+d and e+³He (10×166) samples were streamed, for
the plans/02 §1.5.3 control. The reconstructed campaign subtree was never listed,
e+Au never opened, no podio/edm4eic reader exists, and the campaign-infrastructure
claims of this bullet (`simulation_campaign_hepmc3` + condor, geometry-filename
species routing, the `*_ABCONV` convention) appear nowhere but here and have never
been checked against the endpoint.*
☑ *2026-09-15: checked against the live endpoint and the container. The RECO
subtree is now on the record — `tools/analysis/list_campaign_tree.py` walks
`xrdfs root://dtn-eic.jlab.org ls` to a fixed depth and writes a dated listing
(`--path /volatile/eic/EPIC/RECO --depth 3`); `/volatile/eic/EPIC/RECO` holds 64
entries, 56 versioned campaigns `22.10.0`…`26.07.2` plus `main`,
`alternative-geometries-img-ecal`, `mktest`, four `25.10.{0,1,2,3}_manifest.txt`
and `transfer_manifest.py`. **Claim 1 (`simulation_campaign_hepmc3` + condor) —
cited**: the repository exists and is live (`scripts/run.sh`, pushed 2026-08-31),
and that script is the campaign job, keying off `_CONDOR_CREDS`,
`_CONDOR_SCRATCH_DIR`, `OSG_WN_TMP` and `GLIDEIN_Site` and writing its output with
`xrdcp` to `RECO/${DETECTOR_VERSION}/${DETECTOR_CONFIG}/…`, which is exactly the
tree listed above; the OSG half is visible from the endpoint alone in
`RECO/transfer_manifest.py`, which pulls each manifest line to
`/ospool/uc-shared/project/ePIC/RECO/` over Pelican/OSDF. **Claim 2
(`${DETECTOR_CONFIG}_${EBEAM}x${PBEAM}.xml`) — cited, with the species half
sharpened**: `run.sh:235` (npsim `--compactFile`) and `:275` (eicrecon
`-Pdd4hep:xml_files`) both build
`${DETECTOR_PATH}/${DETECTOR_CONFIG}${EBEAM:+${PBEAM:+_${EBEAM}x${PBEAM}}}.xml`,
and the container ships 44 `epic_craterlake*` compacts, 36 of them carrying an
`_${EBEAM}x${PBEAM}` tag and 25 of those a further species tag (`_H2`, `_He3`,
`_Au197` and one bare `_Au`, `_Cu63`, `_Ru96`, `_Pb207/8`).
The template has no separate species slot, and the 25.10.3 e+Au campaign ran at
`10x100`/`5x41` under the RECO directory `epic_craterlake_without_zdc`, for which
the only compacts that exist are `epic_craterlake_without_zdc_{10x100,5x41}_Au197.xml`
(there is no `epic_craterlake_without_zdc.xml` and no `…_without_zdc_10x100.xml`)
— so the species must travel inside `PBEAM`. That last step is an inference from
the two lists, not a reading of a submit file; the submit files are not on the
endpoint. **No `epic_craterlake_*_Li6/Li7*.xml` exists**, which is the
beamline-XML artifact the risk table below already lists.
**Claim 3 (`*_ABCONV` upstream) — cited**: the suffix is a dataset-level directory
in EVGEN (`EVGEN/EXCLUSIVE/{DVCS,DIFFRACTIVE_RHO,MESON_SF,…}_ABCONV`) that the
campaign mirrors into RECO (`RECO/25.12.0/epic_craterlake/EXCLUSIVE/DVCS_ABCONV`,
and under `SIDIS/` 103 609 `D0_ABCONV` files in the 25.10.3 manifest, each named
`…_q2_1to10000_ab_run368.0024.eicrecon.edm4eic.root` — that dataset's path,
`…/epic_craterlake_without_zdc/SIDIS/D0_ABCONV/HFsim-BeAGLE/BeAGLE1.03.01-2.0/eAu/10x100/q2_1to10000/`,
is the deeper shape the bullet's `…` stands for), and `run.sh` contains no
afterburner step at all — the campaign never afterburns, it consumes afterburnt
input. Two gaps from 2026-08-28 survive: no podio/edm4eic reader was written (a
separate ~8 h, out of scope here), and the BeAGLE e+d (`eH2/10x130`) and e+³He
(`eHe3/{5x41,10x110,18x110,10x166}`) RECO directories exist in 25.12.0 but their
Q² bins list zero files — the volatile copies have been swept, so a reco-level
e+d control would have to be re-produced or pulled from OSG.*
2. Estimate compute (Geant4 e+A ~ min/event → 10⁷ events needs farm/OSG);
   prepare configs to campaign standards so an official e+Li request can go
   through ePIC once endorsed.
3. Final note: "Feasibility of polarized ⁶,⁷Li physics with ePIC at the
   EIC" — reco-level FOMs, tagging performance, requirements flow-down,
   proposed run plan. Reference baseline: ePIC preTDR v3.1
   (Zenodo 10.5281/zenodo.19496158).

---

## Sequencing and gates

2.0 → 2.1 → 2.2 (start with guns — can begin before BeAGLE samples exist)
→ 2.3 → 2.4 → 2.5. Gate each step on validation against e+d / e+³He
references. Calendar: 4–6 months part-time once Phase-1 step 1.5 supplies
samples.

## Risks

| risk | mitigation |
|---|---|
| Li optics/afterburner configs don't exist | verified: 3 concrete artifacts to add (preset, beamline XML, runcard); ⁶Li can proxy d/He4 (same Z/A); engage FF WG early for ⁷Li |
| ⁶Li α-tag fails at IP6 (R=1 beam-blind) | quantify pT-tail acceptance; document IR-8 secondary-focus case; pivot ⁶Li tagging to p/³He channels — *superseded (2026-08-28): the risk materialised (1.7–2.6% at every published optics) and was answered by a mitigation this row predates — the one-plane β\* de-squeeze of Report 1 §6.1, 22–31% at 1/6.8–1/12.8 of the luminosity, with IR-8 (≈20%) as the fallback; the p/³He fragments were evaluated and classified as vetoable breakup backgrounds, not tags (plans/06 §6.2)* |
| Geant4/DD4hep mishandles light-ion or excited-ion primaries | gun tests in 2.1.4; sanitize PDG codes to ground states |
| EICrecon FF matrices proton-tuned | derive Li-rigidity matrices with FF WG — *sharpened (2026-09-15): the "derive" half is done in this repository (`farforward.POT_LEVERS` + `POT_SECOND_ROW`, the full first-order IP6 → pot transfer per configuration), and the request half is now measured rather than read off the source: the reconstruction returns ZERO `ForwardRomanPotRecParticles` for a lithium ion where a 275 GeV proton control gives 72, so `MatrixTransferStaticConfig.h` drops the species rather than mislabelling it.  **Sharpened again (2026-09-16)** on a BeAGLE e+d sample: the same assumption reaches the INCLUSIVE kinematics, where `MCBeamProtons` is empty for any status-4 beam row that is not a proton and every `InclusiveKinematics*` collection is therefore empty in 100/100 events; and where the species IS a proton — a deuteron beam's spectator, R ≈ 0.47 — the off-momentum detector reconstructs it in 80 of 100 events at a 5% momentum and a 0.25 mrad angle, so the matrices are proton-tuned rather than broken* |
| Compute exceeds local resources | guns + small samples locally; campaign production via collaboration |
| Container/geometry churn | pin container per study; record versions in every output dir |
