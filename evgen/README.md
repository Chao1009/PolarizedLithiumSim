# polligen — doubly polarized e+⁶,⁷Li event generation

Event-level Monte Carlo for the *doubly polarized* eA process (plans/05):
polarized electrons on vector/tensor-polarized ⁶,⁷Li, spin-labeled events,
run-plan bookkeeping. Step 5.A (physics kernel + inclusive sampler) and
Step 5.B (tagged mode — spin ⊗ cluster-spectator correlation) are
implemented; the BeAGLE reweighter (5.C) and HepMC3 output (5.D) come
next. Imports `../fastsim/polli_fastsim` — nothing there is duplicated.

## Quick start

```bash
cd evgen
python3 -m pytest tests/ -q            # 66 tests (grid tests auto-skip)
python3 scripts/closure_fom.py --ion 6Li --events 200000 --trials 200
python3 scripts/closure_fom.py --ion 7Li --events 200000 --trials 200
python3 scripts/money_tagged_azz.py --events 400000       # money plot 4
python3 scripts/tagged_polarimetry_7li.py --events 300000
python3 scripts/money_cos2phi.py                          # money plot 5
python3 scripts/money_cos2phi_coherent.py                 # money plot 6
```

## Modules

| module | content |
|---|---|
| `polligen/spin.py` | ρ(m,m′) for J = 1, 3/2: Wigner-d/CG, populations ↔ normalized (vector, tensor, octupole) moments, arbitrary quantization axis, spin-temperature (max-entropy) fills |
| `polligen/xsec.py` | doubly polarized inclusive master formula: HJM spin-1 tensor sector (b₁/b₂, Δ cos 2φ), vector sector A∥ = D·g1/F1 (+ γ-suppressed g_T term, g₂ = g₂^WW), spin-3/2 rank-0/1 exact + rank-2 scenario slots |
| `polligen/bookkeeping.py` | run plans (helicity flips, tensor thirds, transverse fills), relative-luminosity offsets + first-order bias formulas, polarimetry smearing, per-(run,bunch) rng streams |
| `polligen/sample.py` | grid inverse-CDF sampler: per-spin-state Poisson rates (φ-averaged modulation shifts counting rates), φ accept-reject, Mode-W weight matrices |
| `polligen/estimators.py` | analysis-side estimators: helicity-flip, tensor thirds, cos 2φ moment + binned LSQ fit (holey-φ robust), luminosity-corrected yields |
| `polligen/tagged.py` | two-cluster spin model: CG-coupled L-waves → m-dependent spectator densities n_M(k,k̂), embedded-cluster spin populations, pair decomposition; TaggedSampler = spin ⊗ spectator ⊗ DIS ⊗ far-forward routing |
| `polligen/coherent.py` | coherent (intact-g.s.) e+⁶Li channel: scenario coherent fraction/slope, intact-recoil kinematics (R = 1.000 → RP pT-tail only), analytic tag acceptance exp(−B pT_cut²), per-bin tagged-rate projections, breakup veto table (plans/06) |

## Money plots 5–6 (2026-08-10): cos 2φ as projected data points

`money_cos2phi.py`: inclusive gluonometry pseudo-data (φ′ modulation
with statistical error bars) in the four sweet-spot (x, Q²) super-bins,
picked per Q² band from the significance map, plus amplitude-vs-x with
1-year AND 10-year error bars. The Δ(x,Q²) input comes from the unified
registry `polli_fastsim.delta_models` (default `moment_A`: the
sum-rule-constrained ansatz of the merged money_delta suite,
∫xΔdx = −0.012·α_s, with the ⁶Li per-nucleon dilution 1/3;
`moment_B`/`toy` one flag away). Full-luminosity statistics via
per-φ-bin Poisson draws (`sample.phi_histogram_pseudo` — exact for
binned estimators). At P_zz = 0.6 (mid energy) the moment-constrained
amplitudes are (0.6–1.2)×10⁻² at the sweet spots (table-α_s, the
money_delta production convention) with per-bin δA ≈ (1.5–4.5)×10⁻⁴
already in year 1 (10 fb⁻¹/u) — the
measurement discriminates interpretation A vs B, not just zero.

`money_cos2phi_coherent.py`: the coherent channel e⁶Li → e′X ⁶Li(g.s.)
with the intact recoil tagged in the Roman-Pot near-beam pT tail
(R = 1.000 exactly — beam-blind below the pT cut): acceptance
exp(−B pT_cut²) = 13.5% [9–20%] with high-acceptance optics vs 4×10⁻⁵
with high-divergence; ~10⁸ tagged events at 100 fb⁻¹/u (scenario
f₀ = 0.04), best-bin δA ≈ 6×10⁻⁴ (5σ floor at a 0.3% modulation).
The modulation is anchored (plans/06 §6.4b): a t-linear deformation
term scaled from the polarized-deuteron CGC calculation (PLB
858:139053, digitized in `coherent.MANTYSAARI_A2_DEUTERON`) giving
⟨a₂⟩_tag ≈ 0.036 [0.018–0.059] at P_zz = 0.6 with a predicted sign
flip vs the deuteron, plus a flat gluon-transversity scenario
(3×10⁻³–10⁻²) — separable by t-shape and x_P dependence.
Backgrounds and the RP charge-ID question: plans/06.

## Step-5.A validation gates (plans/05 §5.4) — all passing

1. **ρ moments, all axes** — exact against analytic Wigner rotations
   (`test_spin.py`).
2. **Master formula ↔ `asymmetries.py`** — bin-wise identity (rtol 1e-12)
   for A∥, A_zz, A_cos2φ, toy *and* CT18/NNPDFpol grid backends
   (`test_xsec_identity.py`).
3. **Pseudo-experiment estimator closure** — pulls unbiased, spreads match
   `err_a_parallel`/`err_azz`/`err_cos2phi_amplitude` within 15%
   (`test_pseudoexp.py`); full FOM-map closure over ~65 x-bins in
   `closure_fom_{6,7}Li.png` (spread/analytic compatible with 1 across
   2.5 decades in x for all three observables).
4. **φ-modulation recovery** — cos 2φ amplitude unbiased with uniform and
   holey acceptance (binned fit; the naive moment demonstrably breaks).

**First systematics numbers** (relative-luminosity offset, validated at
δ = 2% by pseudo-experiments, quoted at the plans/05 reference δ = 10⁻⁴,
in-ring P_e = 0.7, P_z = 0.7, P_zz = 0.6):
bias(A_zz) = −(2/3)δ/P_zz ≈ −1.1×10⁻⁴ and bias(A∥) = δ/(2P_eP_z) ≈
+1.0×10⁻⁴ for the naive estimators — below the ~10⁻³ per-x-bin statistical
floor at 10 fb⁻¹/u but not negligible after combining bins; the
luminosity-corrected estimators remove both exactly.

## Conventions and caveats

- Populations are ordered m = +J … −J everywhere.
- All SFs per nucleon (F2A/A, as in `fom.py`); luminosity per nucleon.
- Vector-sector y-factors use the fastsim approximation set:
  A∥ = D·A1 (γ²-suppressed η term dropped); A⊥ = d·(A2 − ξA1)
  = d·γ·((y/2)g₁ + g₂)/F₁ at leading γ (both O(γ) pieces kept — the
  ξA1 term is the same order as A2). The exact Cosyn–Weiss factors
  replace two functions when step 5.B lands.
- Spin-1 vector fills at |P_z| > 2/3 require tensor polarization
  (positivity); `helicity_flip_plan` defaults to spin-temperature
  populations, which record their implied P_zz in `pzz_true`.
- Modulation amplitudes are evaluated at cell centers of the sampling
  grid (default 100×72, ~4× finer than the 40×30 FOM analysis binning);
  adequate for binned estimators, not for unbinned in-cell shapes.
- Spin-3/2 rank-2 structure functions are scenario inputs (zero by
  default) pending the theory note (plans/04 #14); rank-3 dropped.
