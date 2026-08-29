# polli_fastsim — Phase-1 fast simulation

*Reproducing these results?  [../docs/reproduction_manual.md](../docs/reproduction_manual.md)
§3 has every command with its runtime and expected output, including the
dated money-Δ line and the R switch that changes its reach.*

Fast (analytic + sampling) simulation for the polarized ⁶Li/⁷Li @ EIC
feasibility study. Companion to `../plans/02_phase1_event_generation.md`.

## Quick start

```bash
cd fastsim
python3 -m pytest tests/ -q                 # grid tests auto-skip
python3 scripts/phase_space_map.py --ion 7Li --lumi 10 --outdir out
python3 scripts/tagging_acceptance.py --outdir out      # spectator tagging
python3 scripts/money_polemc.py --outdir out            # polarized EMC FOM
python3 scripts/money_b1.py --outdir out                # tensor b1 FOM
python3 scripts/money_delta.py --outdir out             # gluonometry reach
python3 scripts/validate_inputs.py                      # toy vs PDF grids
python3 scripts/money_delta_realistic.py                # L_5σ reach (frozen July R)
python3 scripts/money_delta_20260729.py                 # dated production, 16 PNGs
python3 scripts/_check_reco_mask_invariants.py          # its S1–S6 static guard
```

`phase_space_map.py`, `money_polemc.py`, `money_b1.py`,
`money_delta.py`, `money_delta_realistic.py` and
`coverage_and_stat_maps.py` all take `--run-share` (default 1.0), this
observable's share of the programme year; every figure and number here is
the share-1 one, and a non-default share writes to its own file name.

The two dated `money_delta_*` scripts default to the R = σ_L/σ_T form
their notes were written with, defect and all; `--r-model theta-log` and
`--r-model r1998` re-run them with the corrected and the published fit
(plans/08 C3). `money_delta_20260729.py` additionally needs the EPPS21
⁶Li grid; `money_delta_realistic.py` runs on CT18NLO alone.

One-time PDF-grid setup (optional; toys are the default backend):
```bash
pip3 install --user parton && python3 -m parton update
yes | python3 -m parton install CT18NLO
yes | python3 -m parton install NNPDFpol11_100
yes | python3 -m parton install EPPS21nlo_CT18Anlo_Li6   # money_delta_20260729
```

## Modules

| module | content |
|---|---|
| `beams.py` | species (d, ³He, ⁶Li, ⁷Li), γ-matched energies with rigidity-capped tops (plans/10 A0), verified ⁷Li P_p/P_n |
| `kinematics.py` | DIS variables, scattered-electron lab kinematics, acceptance masks |
| `structure.py` | **TOY** F2 (±40% vs CT18, see validate_inputs) + `PartonF2` grid backend (five flavours; caveat 5), nuclear builder, NC cross section |
| `polarized.py` | **TOY** g1 + `PartonG1` (NNPDFpol11, three flavours); **digitized** theory curves in `data/` (`data/SOURCES.md`): CBT PLB 642:210 Fig. 6 and TMT PLB 783:247 Fig. 4 polarized EMC, Miller PRC 89:045203 Fig. 5 and CDKS PRD 95:074036 Fig. 4 b₁; Δ scenarios |
| `delta_models.py` | **unified Δ(x,Q²) model registry** — single home for all double-helicity-flip models: `toy`, sum-rule-constrained `moment_A`/`moment_B` (ported from the `money_delta` suite; ∫xΔdx = −0.012·α_s), shape variants, per-nucleon dilution convention (plans/04 #6); every consumer switches by name |
| `asymmetries.py` | spin-1 master-formula asymmetries (A∥, A_zz, A_cos2φ) + error estimators (toy-MC validated, `tests/test_closure.py`) |
| `fom.py` | luminosity scenarios → events/bin → δ(observable); Q²-combination helper |
| `spectator.py` | α+d / α+t cluster momentum densities (S/P-wave), lab boost → (pT, θ, R); both fragments of one breakup jointly (`breakup_lab_kinematics`, plans/09 B4) |
| `farforward.py` | verified far-forward windows (RP/OMD/B0/ZDC) + rigidity routing |

## First results (TOY inputs, statistical only — headline numbers)

Every number below assumes the **full 10 fb⁻¹/u programme year in its own
configuration** — its own spin states, its own far-forward optics, its own
isotope — which is `--run-share 1`, the default of every script here.  They
are alternatives, not a programme: a run plan that gives an observable a
share *f* of the year multiplies its statistical errors by 1/√*f* and
leaves any luminosity quoted as a reach (`L_5σ`) unchanged, so the years to
that reach go as 1/*f*.  The options are priced in plans/07 WP2; the
arithmetic is pinned by `tests/test_run_share.py`.

- **Tagging acceptance at IP6** (cluster model, β=0.3 GeV central):
  ⁷Li α-tag ≈ **97% into the Roman Pots** at every configuration and
  optics; ⁶Li α-tag **1.9 / 1.7 / 2.6%** at the Yellow Report optics of
  each configuration (the near-beam tail is inside the 10σ rectangle
  10(σ_h, σ_v) = 2.2×3.8 / 1.8×1.8 / 0.92×0.92 mrad, so what survives is
  the 1.5% off-rigidity slice below R = 0.95 plus 0.1–1.0 points of
  over-rigid inner branch) and **28–35%** at the lithium tagging
  optics of Report 1 §6.1 at 1/7–1/13 of the luminosity; ⁷Li t-tag
  **78 / 92 / 94%**
  (`scripts/tagging_acceptance.py`, 2026-08-28, plans/10). The t-tag was
  ~ 0 while an over-rigid fragment was lost by construction; the full
  simulation of plans/09 B1 puts the R = 1.286 triton on Roman-Pot silicon
  at dx = +66 mm in 60 of 60 events with a ZDC deposit behind it, and
  `farforward.route_charged` routes it there against the per-configuration
  blind block, 48 / 32 / 16 mm. The 3–9% and
  1.85% quoted earlier applied a proton p_T threshold, then one
  proton-derived 73 μrad divergence, at every configuration and
  azimuth-blind — a circle rather than the 10σ rectangle.
- **Gluonometry**: 5σ on Δ/F₁ = 10⁻³ (Sather–Schmidt scale) at
  **~66–155 fb⁻¹/u** (LOW/MID/TOP = 67.5 / 65.8 / 155.1 fb⁻¹/u; caveat
  2026-08-27 — `money_delta_realistic.py` carries its own pre-correction
  27.5 / 50 GeV/u low and mid configurations, so only the TOP number is
  at a machine configuration; `money_delta.py` at the γ-matched energies
  gives 16.7 / 16.3 / 21.8 fb⁻¹/u for the toy at P_zz = 0.8) with
  CT18 grid inputs, the 2-of-6-diluted P_zz = 0.267 and the published
  R₁₉₉₈ — inside the 1–100 fb⁻¹/u plausible-program band at 5 × 27.5 and
  10 × 50, above it at 18 × 137.5.
  Re-derived 2026-08-26 (plans/08 C3):
  `python3 scripts/money_delta_realistic.py --r-model r1998 --configs low,mid,top`.
  The two numbers this replaces were both wrong in the same direction:
  the old **~25–37 fb⁻¹/u** headline used P_zz = 0.8 with no 2-of-6
  dilution (code review 2026-08-25, S4), and the 131/275 fb⁻¹/u that
  corrected it (still reproducible with `--r-model simplified`, the
  default) used an `r1998` whose Θ multiplied all three terms of Abe
  et al. Eq. (2) and so returned R = 1.000 for x ≲ 0.1 (S1). Restricting
  Θ to the log term alone moves MID 131.26 → 67.11 and TOP
  274.64 → 156.91; the published three-form fit
  (`structure.r1998`) gives 65.80 / 155.12 — inside the 63–69 / 152–164
  the code review predicted.
  A defect found on 2026-08-28 and corrected in all three scripts left
  every number above unchanged: the ≥ 10-events-per-bin floor was applied
  at the 1 fb⁻¹/u the σ² is normalised to rather than at the luminosity
  the reach is quoted at, so bins that would clear ten events at a 17–275
  fb⁻¹/u reach were discarded and the reach scaled up from the truncated
  sum.  Applying it at the reach (`money_delta.reach_from_terms`, now
  imported by `diag_sig2_grid.py` and `money_delta_realistic.py` instead
  of copied) changes the mask by −15 to +3 bins of the 311–509 accepted,
  and those bins carry so little of σ² that L_5σ moves by at most
  3×10⁻⁹ relative anywhere on the plotted Δ/F₁ range: 16.719 / 16.332 /
  21.811 and 131.26 / 274.64 fb⁻¹/u are unchanged to every digit printed.
- **Polarized EMC**: δΔR ≈ 4.2% per x-bin at x = 0.09 and 4.1% at 0.28 at
  10 fb⁻¹/u (grid inputs, 3 energies combined; 6.0% at x = 0.45, 18.9% at
  0.71).  On the **digitized** curves (2026-08-28) the CBT–TMT separation
  is 0.044 / 0.040 / 0.034 / 0.059 at x = 0.09 / 0.28 / 0.45 / 0.71, i.e.
  1.04 / 0.98 / 0.57 / 0.31 σ per bin at 10 fb⁻¹/u on the grid inputs and
  3.3 / 3.1 / 1.8 / 1.0 σ at 100 (0.92 / 0.78 / 0.55 / 0.48 σ on toy).
  **Read only the valence window.** The ⁷Li ← nuclear-matter transfer is
  one constant fitted over 0.35 < x < 0.65 and applied at every x, so
  below that window the plotted separation is the transfer and not the
  papers: the two PUBLISHED polarized curves agree to better than 0.008
  over 0.028 < x < 0.30 (0.002 at x = 0.09, 0.0006 at 0.14).  Inside the
  window the comparison is real and fades across it: the transferred TMT
  depletion tracks ⁷Li's own unpolarized 0.034 / 0.048 / 0.087 at
  x = 0.40 / 0.45 / 0.65 to within 0.005 against CBT's 0.077 / 0.082 /
  0.094, so ΔR separates by 0.040 at x = 0.36 and 0.011 at 0.65.
  `money_polemc.py`
  shades that window in both panels, draws TMT's published nuclear-matter
  curve untransferred beside the transferred one, and prints both
  separations and the window-restricted reach: best bin
  x = 0.355 at 0.84 σ (10 fb⁻¹/u) and 2.66 σ (100), 0.72 / 2.27 σ on the
  toy inputs the published PNG draws.  The unrestricted best bin
  (x = 0.141, 1.16 / 3.65 σ, 5 of 23 bins above 1 σ) is the transfer's,
  not a prediction.  This retires the "≈ 5σ at x ≈ 0.5–0.7" headline all
  the same: with the constants 2 and 1 the two camps separated as
  |1 − R_EMC(x)|, which grows monotonically with x, while the published
  curves converge above x ≈ 0.6 — the ratio of effects bottoms out at 1.06
  near x = 0.70.  `python3 scripts/money_polemc.py --ion 7Li
  --pdf grid` (`--emc-mode constant` reproduces the old figure).
- **b₁**: δA_zz ≈ 0.9×10⁻⁴ (x < 0.05) to 1.2×10⁻³ (x ≈ 0.56) per x-bin,
  combined over Q² and the three energies (`money_b1.py`), at 10 fb⁻¹/u,
  P_zz = 0.6–0.8 — unchanged.  What changed is the **signal**: with the
  digitized Miller b₁ and the rank-2 transfer 0.921947 × 2/6 = 0.3073,
  |A_zz| is 2.4×10⁻⁴ at x = 0.005 rising to 1.4×10⁻³ at x = 0.07 and
  3.3×10⁻³ at 0.5, i.e. 1.7 / 4.8 / 7.6 / 10.9 / 5.6 / 5.8 σ per bin at
  x = 0.0035 / 0.009 / 0.028 / 0.071 / 0.18 / 0.45 (P_zz = 0.6).  The CDKS
  convolution scenario is |A_zz| = 10⁻⁵ at low x rising to 4×10⁻⁴ at
  x = 0.5, below 0.2 σ everywhere: the two
  camps are now "measurable" and "not measurable", not a factor ten.  The
  signal curve is drawn at Q² = 4 GeV² while the errors combine every
  accessible Q², so the low-x end understates what the low-Q² bins there
  would give.  `--transfer legacy` restores the pre-2026-08-28 0.87 × 1
  (31 σ at x = 0.07).

## Big caveats (by design — see plans/02 steps 1.2–1.5)

1. TOY inputs remain for F2, g1 and Δ.  The polarized-EMC and b₁ curves
   are no longer among them: since 2026-08-28 they are the published
   figures, digitized into `polli_fastsim/data/` (provenance and the
   extraction commands in `data/SOURCES.md`; `tools/digitize_figure.py`
   re-derives them).  Three caveats travel with the tables — the CBT curve
   covers x = 0.028–0.871 and the Miller b₁ x = 0.010–0.900, so both money
   plots draw a constant extrapolation below their leftmost bins; TMT is
   nuclear matter at Q² = 10 GeV² against CBT's ⁷Li at 5, which is why the
   two camps are compared through the ratio of EMC effects on one common
   ⁷Li baseline rather than by subtracting published R values — a single
   valence-window factor, so the comparison is only meaningful inside that
   window (see the polarized-EMC bullet above); and no ⁶Li
   b₁ prediction exists from anyone, so the embedded-deuteron transfer is
   still our own inference (plans/04 #9).  That transfer is also applied
   inconsistently across the repository: `money_b1.py` carries it, while
   the twelve `evgen` scripts that pass `toy_b1` straight to
   `InclusiveKernel(beams.LI6, …)` give their ⁶Li kernel the deuteron's b₁
   undiluted — a factor 3.3 too large — which is tolerated only because b₁
   enters the kernel through w_avg alone, at the 10⁻³ level, and never
   through the cos 2φ amplitude.  Read a ⁶Li b₁ off `money_b1.py`, not off
   the generator.  The unpolarized EMC ratio is
   still the hand-written 12-point table awaiting EPPS21 (plans/02 step
   1.2.1).
2. Ion in-ring polarizations are placeholders (source targets: P_z ≥ 0.90,
   P_zz ≥ 0.80); survival through EBIS+ring is open (plans/04 #1).
3. Cluster-spectator model: two-parameter wave functions; the ⁶Li tail
   (hence its IP6 acceptance) is genuinely model-dominated → VMC densities
   / BeAGLE needed (plans/02 step 1.5.3).
4. Far-forward windows are Phase-1 parameterizations; the near-beam band
   and dispersion assumptions need Phase-2 optics (plans/04 #11).
5. Grid backends mix flavour schemes: `PartonF2` is five-flavour (d u s c
   b) and `PartonG1` three (d u s), so every g1/F1 ratio is the physical
   one — NNPDFpol1.1 sets Δc = Δb = 0 — and is smaller than a
   light-flavour-only ratio by the heavy-quark share of F2.  Measured on
   the bins `money_polemc.py` combines: **7.8% of F2A** event-weighted
   over all of them, 10% at x = 10⁻³–10⁻² and up to 25% at the highest Q²
   there, but **0.65% at x = 0.3–0.5** and 0.23% at 0.5–0.7 — the valence
   window in which the two CBT/TMT curves actually differ, where it is one
   to two orders of magnitude below the 0.03–0.04 separation and costs the
   headline nothing (0.31% at x = 0.35 to 0.03% at 0.65 at Q² = 4 GeV²).
   At the low-x bins `money_polemc.py` also plots it is **2.6% at x = 0.09
   and 1.4% at 0.14** (4.3% and 2.4% at Q² = 10), i.e. comparable to the
   relative separation drawn there — a second reason, besides the transfer
   of the polarized-EMC bullet above, not to read those bins as the
   discriminating ones (2026-08-28, `tests/test_grids.py`).
