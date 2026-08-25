# polli_fastsim — Phase-1 fast simulation

Fast (analytic + sampling) simulation for the polarized ⁶Li/⁷Li @ EIC
feasibility study. Companion to `../plans/02_phase1_event_generation.md`.

## Quick start

```bash
cd fastsim
python3 -m pytest tests/ -q                 # 24 tests (grid tests auto-skip)
python3 scripts/phase_space_map.py --ion 7Li --lumi 10 --outdir out
python3 scripts/tagging_acceptance.py --outdir out      # spectator tagging
python3 scripts/money_polemc.py --outdir out            # polarized EMC FOM
python3 scripts/money_b1.py --outdir out                # tensor b1 FOM
python3 scripts/money_delta.py --outdir out             # gluonometry reach
python3 scripts/validate_inputs.py                      # toy vs PDF grids
```

One-time PDF-grid setup (optional; toys are the default backend):
```bash
pip3 install --user parton && python3 -m parton update
yes | python3 -m parton install CT18NLO
yes | python3 -m parton install NNPDFpol11_100
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
  ⁷Li α-tag ≈ **96–99% into the Roman Pots** (optics-robust);
  ⁶Li α-tag ≈ **3–9%** (high-acceptance optics; tail-dominated, ×3 model
  spread) collapsing to 1–2% with high-divergence optics; ⁷Li t-tag ~ 0.
  Caveat (code review 2026-08-25, S8): these apply the 0.20 GeV *proton*
  p_T threshold to the α; the 10σ envelope is angular, i.e. 0.40 GeV for
  a 4 × 137.5 GeV α, which gives ≈ 0.7–2.7%.
- **Gluonometry**: 5σ on Δ/F₁ = 10⁻³ (Sather–Schmidt scale) at
  **~25–37 fb⁻¹/u** with CT18 grid inputs (15–22 with the toy F2) —
  inside the plausible program either way. Caveat (code review
  2026-08-25, S1/S4): this headline uses P_zz = 0.8 without the 2-of-6
  dilution; the scripts' central P_zz = 0.267 gives 131/275 fb⁻¹/u, and
  the money_delta scripts' `r1998` returns R = 1 at low x (≈ 63–69 /
  152–164 fb⁻¹/u once corrected) — to be re-derived.
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
