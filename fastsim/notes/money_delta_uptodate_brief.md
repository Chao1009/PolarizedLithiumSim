# Money-Δ: Internal Working Summary (2026-07-31)

**Source:** distilled from `money_delta_uptodate.md` (master reference),
incorporating the 2026-07-31 reco-selection consistency fix audit
(`money_delta_note_2026-07-29_fix.md`). All numbers are sourced from those
documents; nothing is invented here.

---

## Purpose and Current Status

Goal: project EIC sensitivity to the tensor structure function Δ(x, Q²) of
⁶Li via the tensor asymmetry A_cos2φ. Four study phases complete (2026-07-16
through 2026-07-29); a reco-selection consistency fix to the Phase IV script
was applied and production-rerun validated on 2026-07-31
(`money_delta_note_2026-07-29_fix.md`). Discovery is feasible within a
1-year EIC program under all tested model combinations. The key unresolved
systematic is the choice of Δ ansatz (Interpretation A vs B, factor 2–8 in
reach).

---

## Core Observable

**Primary:** 5σ discovery luminosity L_5σ (fb⁻¹/nucleon) and fractional
asymmetry uncertainty δA/A at fixed luminosity.

**Secondary:** φ-differential yield modulation

    y(φ) = (N_φ − N_flat) / N_flat = P_zz · <A_cos2φ> · cos(2φ)

The asymmetry amplitude per (x, Q², y) bin (Hoodbhoy–Jaffe–Manohar):

    A_cos2φ = −[(1−y)/y²] · Δ / [F₁ + ((1−y)/(xy²))·F₂]

Combined significance across bins (Fisher sum):

    σ² = Σ_bins  A_cos2φ,bin² · P_zz² · N_bin / 2

    δA/A = 1/σ ;   L_5σ requires σ = 5  (i.e., δA/A = 0.20)

**Warning:** per-bin median |δA/A| is ~35–65× larger than the combined δA/A.
Do not quote per-bin medians as combined sensitivity.

Luminosity scale: 10 fb⁻¹/nucleon ≈ 1 year at nominal EIC luminosity.

---

## Beam Configurations

| Label | E_e [GeV] | p_ion [GeV/u] | √s [GeV/u] | <Q²> [GeV²] |
|-------|-----------|---------------|------------|-------------|
| LOW   | 5         | 27.5          | ~23        | 5.92        |
| MID   | 10        | 50            | ~45        | 7.38        |
| TOP   | 18        | 137.5         | ~100       | 10.26       |

Ion is ⁶Li throughout. <Q²> is the rate-weighted mean over accepted bins.

---

## Models Considered

**Δ ansatz — two interpretations:**

- **Interpretation A:** Δ(x,Q²) = s · α_s(Q²) · F₁(x,Q²) · x^α · (1−x)^β
  Template applied to Δ/F₁; sum-rule solver is a numerical integral over F₁;
  A varies ~5% across configs. More optimistic.

- **Interpretation B:** Δ(x,Q²) = s · α_s(Q²) · x^α · (1−x)^β
  Template applied to Δ directly; analytic Beta-function solver; A is strictly
  config-independent. Reach is 2–8× worse in L_5σ.

**x-shape variants** (applied identically under both interpretations):

| Label  | α   | β | Shape character             | x_peak |
|--------|-----|---|-----------------------------|--------|
| low_x  | 0.3 | 4 | Soft, Regge-like at small x | ~0.07  |
| mid_x  | 0.7 | 3 | Moderate peak               | ~0.19  |
| high_x | 1.5 | 2 | Suppressed at small x       | ~0.43  |

Shape systematic: factor ~2–3 in L_5σ within each interpretation.

**Integration cases (Phase III/IV):**

| Case | Description                              |
|------|------------------------------------------|
| 1    | Single peak-rate (x, Q²) bin            |
| 2    | Q² = 2.43 GeV² slice, all accepted x   |
| 3    | All accepted (x, Q²) bins — **preferred** |

---

## Inputs / Ingredients

| Ingredient       | Source / Value                                      |
|------------------|-----------------------------------------------------|
| Nuclear PDF      | EPPS21nlo_CT18Anlo_Li6                              |
| R = σ_L/σ_T     | R1998 parametrization → F₁ = F₂/(2x(1+R))         |
| α_s(Q²)         | `parton` package table                              |
| P_zz             | 0.267 (Cloet tensor polarization for ⁶Li)           |
| Sum-rule c (bag) | −0.012 (Sather–Schmidt 1993)                        |
| Sum-rule c (lat) | −0.009 (Detmold–Shanahan)                           |
| Detector model   | ePIC tracking-only σ_p/p, σ_θ (7 η regions) +      |
|                  | η-dependent ε_eID (ATHENA + ECCE anchors, 9 points) |

|A_bag| / |A_lat| = 4/3 exactly → L_5σ(lattice) ≈ 1.78× larger than bag.

---

## Stable Conclusions

1. **Discovery within 1 year is feasible** (pre-detector, ideal simulation)
   under both interpretations, all three x-shapes, and both sum-rule
   constraints. Worst case (Interpretation B, lattice, TOP, mid_x):
   L_5σ ≈ 0.82 fb⁻¹/u, still < 1 year. Detector validation has been
   performed under one model only (Interpretation A · mid_x · bag); see
   "Detector-Tested Conclusion Scope" below.

2. **At 100 fb⁻¹, statistics are no longer limiting** (δA/A ~ 0.5–1% under
   Interpretation A, bag); systematics become the bottleneck.

3. **Use the fully-integrated analysis (Case 3).** Pre-detector, Case 3 gives
   the best S/N for all three configs: about 2–3.5× over Case 2 and about
   5–25× over Case 1 in the tested (Interp. A, mid_x, bag) phi study.
   Under the tested detector model, Case 3 retains ~92% of pre-detector S/N.

4. **|A_bag|/|A_lat| = 4/3 exactly** under both interpretations.

---

## Model-Dependent Findings

| Finding                               | Magnitude            | Condition                              |
|---------------------------------------|----------------------|----------------------------------------|
| Interp. A vs B reach penalty          | 2–8× in L_5σ        | Leading modeling systematic            |
| Within-interpretation shape spread    | ~2–3× in L_5σ       | Across low_x / mid_x / high_x         |
| TOP Case 1 kinematic trough           | 13× suppression      | mid_x shape (α=0.7) at x≈0.002        |
| TOP Case 1 trough severity            | ~2× suppression      | low_x shape (α=0.3) instead            |
| P_zz uncertainty (~1.5×)             | ~2.3× in σ², ~1.5× in L_5σ | If P_zz_eff deviates from 0.267 |
| ε_eID model uncertainty               | ±5–10% efficiency    | ATHENA+ECCE anchors, forward region    |

---

## Detector-Tested Conclusion Scope

All detector conclusions apply **only** to:
> Interpretation A · mid_x shape · bag constraint · ePIC tracking-only smearing
> + ATHENA/ECCE η-dependent ε_eID · 10 fb⁻¹

Under that model (2026-07-29; conclusions confirmed by the 2026-07-31
production rerun after the reco-selection consistency fix):

- **Case 3 (integrated):** retains ~92% of pre-detector S/N across all configs.
  The 7–8% loss is from ε_eID < 1, not bin migration. Post-fix Case 3
  ⟨A_cos2φ⟩ values shifted by 2–4% relative to pre-fix values, consistent
  with the prior selection inconsistency being corrected; the integrated
  Case 3 recommendation is unchanged and reinforced.

- **Case 1 (single peak-rate bin):** changes by factors of 2–3 in either
  direction; direction is config-specific (MID drops to 42%, TOP rises to 214%
  — the latter is an artifact of migration moving the reco peak away from the
  worst available bin, not a real gain). Post-fix Case 1 values are
  unchanged (<0.1% shift).

- **TOP Case 2 ~5.1% shift** (post-fix vs pre-fix): the largest Case 2/3
  shift in the production rerun, marginally above the ≤5% plan tolerance.
  Recorded as a pass-with-note detail; it does not represent a headline
  reversal of conclusions.

Not yet tested under Interpretation B, other x-shapes, or the lattice
constraint.

---

## Reach Tables — Selected Highlights (bag constraint; lattice scales L_5σ × 1.78)

**Interpretation A** (Phase II, 2026-07-24):

| Config | Shape  | δA/A @ 10 fb⁻¹ | δA/A @ 100 fb⁻¹ | L_5σ [fb⁻¹/u] |
|--------|--------|-----------------|-----------------|----------------|
| LOW    | low_x  | 0.0059          | 0.0019          | 0.01           |
| MID    | mid_x  | 0.0086          | 0.0027          | 0.02           |
| TOP    | high_x | 0.0343          | 0.0108          | 0.29           |

Worst case Interp. A: TOP + high_x, L_5σ ≈ 0.29 fb⁻¹/u (~80 hr at EIC).

**Interpretation B** (Phase II, 2026-07-27):

| Config | Shape  | δA/A @ 10 fb⁻¹ | δA/A @ 100 fb⁻¹ | L_5σ [fb⁻¹/u] |
|--------|--------|-----------------|-----------------|----------------|
| LOW    | low_x  | 0.0109          | 0.0034          | 0.03           |
| MID    | mid_x  | 0.0175          | 0.0055          | 0.08           |
| TOP    | mid_x  | 0.0428          | 0.0135          | 0.46           |

Worst case Interp. B: TOP + mid_x, L_5σ ≈ 0.46 fb⁻¹/u (~128 hr at EIC).

Full 9-row tables for each interpretation are in §6 of the master note.

---

## Open Uncertainties

**Tier A — physics critical path:**

| Priority | Item                                        | Impact                             |
|----------|---------------------------------------------|------------------------------------|
| A1 (top) | P_zz = 0.267 convention for ⁶Li             | Factor ~1.5 in P_zz → factor ~1.5 in L_5σ |
| A2       | Δ x-shape from first principles             | Resolves Interp. A vs B (factor 2–8)       |
| A3       | Theorist consultation on ansatz choice      | Qualitative question about operator structure |
| A4       | Reach quotes must span both interpretations | Required for publication                   |
| A5       | Systematics study at 100 fb⁻¹              | P_zz calibration, F₂A normalization, radiative corrections |

**Tier B — detector realism:**

| Item | Description                                              |
|------|----------------------------------------------------------|
| B1   | ECal smearing as alternative/combined electron reco      |
| B2   | Replace ATHENA+ECCE ε_eID with official ePIC curve       |
| B3   | QED radiative corrections (~5–15% on reco x, Q², y)     |
| B4   | Non-uniform azimuthal acceptance (ePIC φ gaps)           |
| B5   | η-weighted ε_eID per bin for precision δA_bin            |
| B6   | Increase MC per bin from 1000 to ≥5000 events            |

**Tier C (long-term):** NNLO QCD, target-mass corrections, polarized-electron
observables (A_LT, A_LL), two-gluon operator contribution to Δ at small x.

---

## Next Steps

1. Resolve P_zz convention for ⁶Li (consult I. Cloet or nuclear-structure
   literature). **Highest priority.**
2. Obtain a motivated first-principles Δ x-shape (lattice moment, NJL, or CDKS
   convolution) to distinguish Interpretation A from B.
3. Replace ATHENA+ECCE ε_eID anchors with official ePIC parameterization (B2).
4. Extend Phase IV detector study to Interpretation B and the lattice constraint.
5. Begin systematics study for the 100 fb⁻¹ regime (A5).

---

## Package / File Structure

```
PolarizedLithiumSim/
└── fastsim/
    ├── notes/                              Documentation
    │   ├── money_delta_note_2026-07-16.md  Phase I: infrastructure, first reach
    │   ├── money_delta_note_2026-07-20.md  Phase I: PDF comparison, lum. ref
    │   ├── money_delta_note_2026-07-21.md  Phase I: three-tier scenario
    │   ├── money_delta_note_2026-07-24.md  Phase II: Interpretation A, sum-rule
    │   ├── money_delta_note_2026-07-27.md  Phase II: Interpretation B, A vs B
    │   ├── money_delta_note_2026-07-28.md  Phase III: phi-modulation, Cases 1-3
    │   ├── money_delta_note_2026-07-29.md       Phase IV: detector smearing
    │   ├── money_delta_note_2026-07-29_fix.md  Fix audit: reco-mask consistency, production rerun
    │   ├── money_delta_uptodate.md              Master synthesis (source of this brief)
    │   └── money_delta_uptodate_brief.md        This file
    ├── scripts/                            Run scripts (in order for each phase)
    │   ├── money_delta_20260715.py         Phase I (earlier)
    │   ├── money_delta_20260720.py         Phase I (earlier)
    │   ├── money_delta_20260721.py         Phase I (earlier)
    │   ├── money_delta_20260724.py         Phase II — Interpretation A reach
    │   ├── money_delta_20260725.py         Phase II — Interpretation B reach
    │   ├── money_delta_20260728.py         Phase III — phi-modulation projection
    │   ├── money_delta_20260729.py          Phase IV — detector-realistic (post-fix 2026-07-31)
    │   ├── _check_reco_mask_invariants.py  S1–S6 AST static checker (created 2026-07-31)
    │   ├── money_delta_pdfgrid.py          Utility: PDF grid diagnostics
    │   └── money_delta_realistic.py        Utility (see individual note)
    ├── polli_fastsim/                      Core library
    │   ├── structure.py                    Generic backend: ToyF2, NuclearF2, _safe_xfx,
    │   │                                   r_sigma_lt, dsigma_dx_dq2. EPPS21+R1998
    │   │                                   production setup is script-local (NuclearF2FromGrid).
    │   ├── polarized.py                    ToyG1, EMC ratio scenarios, toy_delta_gluon
    │   │                                   placeholder. A_cos2phi and P_zz are NOT here.
    │   ├── kinematics.py                   (x, Q2, y) grid, acceptance cuts
    │   ├── beams.py                        Ion dataclass, BeamConfig, LOW/MID/TOP configs.
    │   │                                   P_zz is a script/scenario input, not a field here.
    │   ├── asymmetries.py                  a_cos2phi (HJM Eq. 4), err_cos2phi_amplitude,
    │   │                                   azz, a_parallel
    │   └── fom.py                          Scenario dataclass, project_rates (per-bin
    │                                       event counts), project_observables (per-bin
    │                                       asymmetries/errors). Scripts compute sum-rule
    │                                       Fisher sums, L_5sigma, and delta-A/A locally.
    └── out/money_delta/                    Output plots (~50 PNGs across 4 phases)
```

---

## ASCII Flow Chart

```
Physics assumptions
  |  Hoodbhoy-Jaffe-Manohar operator structure of Delta
  |  Sum-rule: int x*Delta dx = c * alpha_s(Q^2)
  |    c = -0.012 (bag, Sather-Schmidt)
  |    c = -0.009 (lattice, Detmold-Shanahan)
  v
Structure inputs
  |  EPPS21 nuclear PDF  ->  F1(x,Q^2), F2(x,Q^2)
  |  R1998               ->  F1 = F2 / (2x(1+R))
  |  parton pkg          ->  alpha_s(Q^2)
  |  P_zz = 0.267        ->  tensor polarization of 6Li
  v
Delta ansatz / sum-rule normalization
  |  Choose interpretation:
  |    A: Delta = s * alpha_s * F1 * x^alpha * (1-x)^beta
  |       solver: s = c / INT[ x*F1*x^a*(1-x)^b dx ]  (numerical)
  |    B: Delta = s * alpha_s * x^alpha * (1-x)^beta
  |       solver: s = c * P_shape / B(alpha+2, beta+1)  (analytic)
  |  Choose x-shape: low_x / mid_x / high_x
  v
Rates and observables
  |  A_cos2phi(x,Q^2,y) = -[(1-y)/y^2] * Delta / [F1 + ((1-y)/(xy^2))*F2]
  |  N_bin = L * dsigma/dx dQ^2 * (bin area)
  |  Optionally: MC smear electron kinematics (ePIC tracking + eID efficiency)
  v
Statistical figures of merit
  |  sigma^2 = SUM_bins  A_bin^2 * P_zz^2 * N_bin / 2   (Fisher sum)
  |  delta-A/A = 1/sigma
  |  L_5sigma = 25 / (sigma^2 / L_ref)
  |
  |  OR  y(phi) = P_zz * <A_cos2phi> * cos(2phi)  [phi-modulation observable]
  |      S/N = P_zz * |<A>| * sqrt(N_total / 2)
  v
Integration strategy
  |  Case 1: single peak-rate bin   -- FRAGILE under detector smearing
  |  Case 2: Q^2 slice              -- intermediate
  |  Case 3: all accepted bins      -- PREFERRED (~92% S/N retained post-detector)
  v
Plots / conclusions
     Reach heatmaps, delta-A/A vs scale, L_5sigma vs scale
     phi-modulation curves (3 signal scales overlaid, Poisson error bars)
     Pre- vs post-detector S/N comparison tables
```

---

## Reproducibility Pointers

**Prerequisites (run once from repo root):**
```bash
python3 -m parton install EPPS21nlo_CT18Anlo_Li6
python3 -m parton install nNNPDF30_nlo_as_0118_A6_Z3
python3 -m parton install CT18NLO
```
Also apply: NumPy 2.0 patch in `parton/pdf.py:231` and `_safe_xfx` helper in
`fastsim/polli_fastsim/structure.py` and `polarized.py`
(documented in `money_delta_note_2026-07-16.md` §5).

**Run the four active phases in order:**
```bash
python3 fastsim/scripts/money_delta_20260724.py --outdir fastsim/out/money_delta  # ~15-20 min
python3 fastsim/scripts/money_delta_20260725.py --outdir fastsim/out/money_delta  # ~15-20 min
python3 fastsim/scripts/money_delta_20260728.py --outdir fastsim/out/money_delta  # ~5-10 min
python3 fastsim/scripts/money_delta_20260729.py --outdir fastsim/out/money_delta  # ~5-10 min  seed=42
```

Output plots land in `fastsim/out/money_delta/` (~50 PNGs total across 4 phases).
Full plot inventory is in master note §13.

---

*Brief last updated: 2026-07-31 (reco-selection consistency fix incorporated).*
*Sources: `money_delta_uptodate.md`; `money_delta_note_2026-07-29_fix.md`.*
*Detector conclusions (§ "Detector-Tested Conclusion Scope") are limited to*
*Interpretation A · mid_x · bag · ePIC tracking-only + ATHENA/ECCE ε_eID · 10 fb⁻¹.*
