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
| `beams.py` | species (d, ³He, ⁶Li, ⁷Li), rigidity-scaled top energies, verified ⁷Li P_p/P_n |
| `kinematics.py` | DIS variables, scattered-electron lab kinematics, acceptance masks |
| `structure.py` | **TOY** F2 (±40% vs CT18, see validate_inputs) + `PartonF2` grid backend, nuclear builder, NC cross section |
| `polarized.py` | **TOY** g1 + `PartonG1` (NNPDFpol11); scenario curves: CBT 2× / TMT 1× polarized EMC, HERMES-like vs convolution b1, Δ scenarios |
| `delta_models.py` | **unified Δ(x,Q²) model registry** — single home for all double-helicity-flip models: `toy`, sum-rule-constrained `moment_A`/`moment_B` (ported from the `money_delta` suite; ∫xΔdx = −0.012·α_s), shape variants, per-nucleon dilution convention (plans/04 #6); every consumer switches by name |
| `asymmetries.py` | spin-1 master-formula asymmetries (A∥, A_zz, A_cos2φ) + error estimators (toy-MC validated, `tests/test_closure.py`) |
| `fom.py` | luminosity scenarios → events/bin → δ(observable); Q²-combination helper |
| `spectator.py` | α+d / α+t cluster momentum densities (S/P-wave), lab boost → (pT, θ, R) |
| `farforward.py` | verified far-forward windows (RP/OMD/B0/ZDC) + rigidity routing |

## First results (TOY inputs, statistical only — headline numbers)

- **Tagging acceptance at IP6** (cluster model, β=0.3 GeV central):
  ⁷Li α-tag ≈ **97% into the Roman Pots** at every configuration and
  optics; ⁶Li α-tag ≈ **1.5–1.7%** at the Yellow Report optics of each
  configuration (the near-beam tail is inside the 10σ rectangle
  10(σ_h, σ_v) = 2.2×3.8 / 1.8×1.8 / 0.92×0.92 mrad; what survives is the
  off-rigidity slice below R = 0.95) and **27–35%** at the lithium tagging
  optics of Report 1 §6.1 at 1/7–1/13 of the luminosity; ⁷Li t-tag ~ 0
  (`scripts/tagging_acceptance.py`, 2026-08-28, plans/10). The 3–9% and
  1.85% quoted earlier applied a proton p_T threshold, then one
  proton-derived 73 μrad divergence, at every configuration.
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
- **Polarized EMC**: δΔR ≈ 2.6–4% per x-bin at x = 0.3–0.5 at 10 fb⁻¹/u
  (grid inputs, 3 energies combined; 12% at x = 0.7); CBT-vs-TMT
  discrimination ≈ 5σ at x ≈ 0.5–0.7 with 100 fb⁻¹/u.
- **b₁**: δA_zz ≈ 10⁻⁴ (x < 0.05) to 10⁻³ (x > 0.3) per x-bin, combined over
  Q² and the three energies (`money_b1.py`), at 10 fb⁻¹/u, P_zz = 0.6–0.8.

## Big caveats (by design — see plans/02 steps 1.2–1.5)

1. TOY/scenario inputs everywhere until step 1.2 completes (CBT/TMT and
   b1 curves are qualitative shapes, not digitized theory).
2. Ion in-ring polarizations are placeholders (source targets: P_z ≥ 0.90,
   P_zz ≥ 0.80); survival through EBIS+ring is open (plans/04 #1).
3. Cluster-spectator model: two-parameter wave functions; the ⁶Li tail
   (hence its IP6 acceptance) is genuinely model-dominated → VMC densities
   / BeAGLE needed (plans/02 step 1.5.3).
4. Far-forward windows are Phase-1 parameterizations; the near-beam band
   and dispersion assumptions need Phase-2 optics (plans/04 #11).
