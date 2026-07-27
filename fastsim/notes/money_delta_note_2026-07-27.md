# Money-Δ Plot: Interpretation A vs B — Dropping F₁ from Δ — 2026-07-27

**Observable.** 5σ discovery luminosity `L_5σ` (fb⁻¹/nucleon) and fractional
asymmetry uncertainty `δA/A` for the `cos(2φ)` double-helicity-flip tensor
asymmetry on transversely tensor-polarized ⁶Li. Today's work introduces
**Interpretation B**: the Δ ansatz drops the F₁ factor from the numerator,
making Δ = s · α_s(Q²) · x^α(1-x)^β. This change makes the sum-rule solver
analytic (Beta function, config-independent) and gives more conservative reach
estimates by factors of 2–8 compared to Interpretation A (2026-07-24).

**Revision history**

- 2026-07-27 (initial): Cloned `money_delta_20260724.py` → `money_delta_20260725.py`
  (note: script name carries the 20260725 date stamp; this note is dated 2026-07-27).
  Single physics change: F₁ removed from the Δ ansatz (Interpretation B). A becomes
  analytic and config-independent. 8 PNGs produced with `_20260725_` prefix.
  See `money_delta_note_2026-07-24.md` for the full prior derivation and formalism.

---

## 1. What was produced today

Script: `fastsim/scripts/money_delta_20260725.py` (~1400 lines, cloned from
`money_delta_20260724.py`). All figures land in `fastsim/out/money_delta/`
with `_20260725_` filename prefix. Ion is ⁶Li throughout. Three beam configs ×
three Δ x-shape variants = 9 (config, shape) combinations; A is evaluated once
per variant (not per config, since it is now config-independent).

### Beam configurations

Same as `money_delta_note_2026-07-24.md` §1:

| Label | E_e [GeV] | p_ion [GeV/u] | √s [GeV/u] |
|---|---|---|---|
| LOW | 5 | 27.5 | ~23 |
| MID | 10 | 50 | ~45 |
| TOP | 18 | 137.5 | ~100 |

### Δ x-shape variants

Same (α, β) grid as before:

| Label | α | β | Character |
|---|---|---|---|
| low_x | 0.3 | 4 | Soft, Regge-like rise at small x |
| mid_x | 0.7 | 3 | Moderate peak |
| high_x | 1.5 | 2 | Suppressed at small x, peaks at larger x |

### Output figures (8 PNGs)

| Filename | Content |
|---|---|
| `money_delta_20260725_p1_bag_10fb.png` | δA/A vs x-shape variant at L = 10 fb⁻¹, bag prediction |
| `money_delta_20260725_p1_bag_100fb.png` | Same at L = 100 fb⁻¹ |
| `money_delta_20260725_p2_bag.png` | L_5σ vs x-shape variant, bag prediction, 3 configs |
| `money_delta_20260725_p3_baglat.png` | δA/A, both bag and lattice predictions overlaid |
| `money_delta_20260725_p4_baglat.png` | L_5σ, both bag and lattice overlaid |
| `money_delta_20260725_perbin_mid_10fb.png` | Per-bin heatmap, MID config, L = 10 fb⁻¹ |
| `money_delta_20260725_perbin_mid_100fb.png` | Per-bin heatmap, MID config, L = 100 fb⁻¹ |
| `money_delta_20260725_perbin_top_10fb.png` | Per-bin heatmap, TOP config, L = 10 fb⁻¹ |
| `money_delta_20260725_perbin_top_100fb.png` | Per-bin heatmap, TOP config, L = 100 fb⁻¹ |

*Note: there are 4 per-bin files (2 configs × 2 luminosities), making 8 PNGs total.*

---

## 2. Formulas updated relative to 2026-07-24

Cross-reference `money_delta_note_2026-07-24.md` §2 for the full equation list
(Eqs. 2, 4–11). Only two equations change under Interpretation B.

### Eq. 1 — Δ ansatz

**Interpretation A (money_delta_20260724.py):**

    Δ(x, Q²) = s · α_s(Q²) · F₁(x, Q²) · x^α · (1−x)^β          (1A)

Ansatz on Δ/F₁. F₁ appears in the numerator of the asymmetry amplitude.

**Interpretation B (money_delta_20260725.py — today):**

    Δ(x, Q²) = s · α_s(Q²) · x^α · (1−x)^β                        (1B)

Ansatz on Δ directly. F₁ no longer boosts the numerator.

### Eq. 3 — Sum-rule solver for A ≡ s

The first-moment sum rule (Eq. 2, unchanged):

    ∫₀¹ x · Δ(x, Q²) dx = c · α_s(Q²)          c = −0.012 (bag), −0.009 (lattice)

**Interpretation A:** Substituting (1A) gives a numerical integral over x·F₁:

    A = c / ∫₀¹ x · F₁(x, ⟨Q²⟩) · x^α · (1−x)^β dx               (3A)

α_s cancels; the integral over F₁ must be evaluated numerically per config
(via the `parton` package α_s and PDF table). A varies ~5% across configs.

**Interpretation B:** Substituting (1B), F₁ drops out entirely:

    A = c · _PEAK_VALS[variant] / B(α+2, β+1)                       (3B)

where B(α+2, β+1) is the Euler Beta function ∫₀¹ x^(α+1)(1-x)^β dx and
`_PEAK_VALS[variant]` is the normalization factor needed to relate the shape
integral to the first moment. A is now **analytic and config-independent** —
no numerical PDF integration is required. α_s cancels between the two sides
of the sum rule as in Interpretation A.

All other formulas in `money_delta_note_2026-07-24.md` §2 (Eqs. 2, 4–11)
are unchanged and apply equally under Interpretation B.

---

## 3. A-value table (Interpretation B)

Under Interpretation B, A is determined by the Beta function and the sum-rule
constant c alone — F₁ evolution plays no role. A is therefore identical across
LOW, MID, and TOP configs.

| Variant | α | β | B(α+2, β+1) | A_bag (c = −0.012) | A_lat (c = −0.009) |
|---|---|---|---|---|---|
| low_x | 0.3 | 4 | 0.02202 | −0.18355 | −0.13766 |
| mid_x | 0.7 | 3 | 0.02242 | −0.08895 | −0.06671 |
| high_x | 1.5 | 2 | 0.02309 | −0.04762 | −0.03571 |

`|A_bag| / |A_lat| = 0.012 / 0.009 = 4/3` exactly — same as Interpretation A
(§3 of `money_delta_note_2026-07-24.md`), since the c-coefficient ratio is
unchanged.

*Contrast with Interpretation A:* the A values there ranged from −0.297 to
−0.391 (low_x/high_x shape, MID config). Interpretation B gives |A| values
3–8× smaller, because removing F₁ from the Δ numerator means the shape integral
that sets A no longer includes the low-x enhancement from F₁.

---

## 4. Reach summary at A_bag (Interpretation B)

Luminosity-scaling conventions are unchanged from `money_delta_note_2026-07-24.md`
§4 (Eqs. 6–11):

    10  fb⁻¹/nucleon  ≈ 1 year at nominal EIC instantaneous luminosity
    100 fb⁻¹/nucleon  ≈ 10 years

σ² scales as A²; lattice entries give σ² × 9/16 of bag values and
L_5σ × 16/9 ≈ 1.78× larger.

| Config | Variant | σ²/(fb⁻¹) | δA/A @ 10 fb⁻¹ | δA/A @ 100 fb⁻¹ | L_5σ [fb⁻¹/u] |
|---|---|---|---|---|---|
| LOW | low_x | 841.03 | 0.0109 | 0.0034 | 0.03 |
| LOW | mid_x | 464.42 | 0.0147 | 0.0046 | 0.05 |
| LOW | high_x | 317.16 | 0.0178 | 0.0056 | 0.08 |
| MID | low_x | 603.89 | 0.0129 | 0.0041 | 0.04 |
| MID | mid_x | 326.04 | 0.0175 | 0.0055 | 0.08 |
| MID | high_x | 309.04 | 0.0180 | 0.0057 | 0.08 |
| TOP | low_x | 98.07 | 0.0319 | 0.0101 | 0.25 |
| TOP | mid_x | 54.58 | 0.0428 | 0.0135 | 0.46 |
| TOP | high_x | 75.28 | 0.0364 | 0.0115 | 0.33 |

**Worst case under Interpretation B:** TOP + mid_x, L_5σ ≈ 0.46 fb⁻¹/u. At
nominal EIC instantaneous luminosity (~3.6 pb⁻¹/nucleon/hour) that is roughly
128 hours — still well within a 1-year program. Under the lattice sum rule
(c = −0.009), the worst case scales to ~0.82 fb⁻¹/u (still < 1 year).

---

## 5. Comparison: Interpretation A vs Interpretation B

The table below uses MID config + mid_x shape as the reference case; the first
four rows compare L_5σ across config/shape combinations. See §3 of
`money_delta_note_2026-07-24.md` for the full Interpretation A A-value table,
and §4 for the Interpretation A reach table.

| Quantity | Interp. A (20260724) | Interp. B (20260725) | B/A ratio |
|---|---|---|---|
| \|A_bag\| MID + mid_x | 0.310 | 0.0890 | 0.29 |
| L_5σ MID + mid_x + bag | 0.02 fb⁻¹/u | 0.08 fb⁻¹/u | 4× worse |
| L_5σ MID + low_x + bag | 0.01 fb⁻¹/u | 0.04 fb⁻¹/u | 4× worse |
| L_5σ TOP + mid_x + bag | 0.06 fb⁻¹/u | 0.46 fb⁻¹/u | 8× worse |
| δA/A MID + mid_x + 10 fb⁻¹ | 0.0086 | 0.0175 | 2× worse |
| δA/A MID + mid_x + 100 fb⁻¹ | 0.0027 | 0.0055 | 2× worse |

**The 4× penalty in L_5σ comes from two factors acting simultaneously:**

1. |A| drops by ~3.5× (from 0.310 to 0.089) because F₁ is no longer
   amplifying the signal in high-F₁, high-rate bins.
2. Those same high-rate bins (low x, moderate Q²) no longer receive the
   proportional signal boost, so the Fisher sum (Eq. 6) is also less
   concentrated in the best bins.

The combined effect is a ~12× drop in σ² per fb⁻¹, which translates to ~4×
worse L_5σ (since L_5σ ∝ 1/σ²) and ~3.5× worse δA/A (since δA/A ∝ 1/√σ²).

---

## 6. Per-bin heatmap diagnostics (Interpretation B)

The three-subpanel layout (|A_bin|, δA_bin, |δA_bin/A_bin|) is unchanged from
`money_delta_note_2026-07-24.md` §5. Numbers below are for the mid_x variant.

| Config | L [fb⁻¹] | \|A_bin\| min/med/max | δA_bin min/med/max | \|δA/A\|_bin min/med/max | Max \|δA/A\| at |
|---|---|---|---|---|---|
| MID | 10 | 2.42e-6 / 3.93e-3 / 0.157 | 1.07e-3 / 3.54e-3 / 0.839 | 8.68e-2 / 1.46 / 822 | x=3.57e-3, Q²=6.74 GeV² |
| MID | 100 | 2.42e-6 / 3.93e-3 / 0.157 | 3.38e-4 / 1.12e-3 / 0.265 | 2.74e-2 / 0.461 / 260 | x=3.57e-3, Q²=6.74 GeV² |
| TOP | 10 | 1.13e-7 / 7.70e-4 / 0.170 | 1.07e-3 / 4.27e-3 / 0.778 | 0.210 / 5.77 / 1.50e4 | x=5.66e-4, Q²=5.23 GeV² |
| TOP | 100 | 1.13e-7 / 7.70e-4 / 0.170 | 3.39e-4 / 1.35e-3 / 0.246 | 6.62e-2 / 1.82 / 4750 | x=5.66e-4, Q²=5.23 GeV² |

**The failure mode moves to small x.** Under Interpretation A (2026-07-24 §5),
the maximum |δA/A| bin was at x = 0.897, Q² = 1379–1776 GeV² (large x, large
Q²). Under Interpretation B, the worst bin is at small x, moderate Q² — the
opposite corner of phase space. The reason: without F₁ in the numerator of Δ,
low-x bins (which have the largest F₁ values) no longer get a signal boost
proportional to their rate. The signal is now more uniform across x, and the
worst-constrained bins are those at the smallest x where phase space is first
opening up but statistics are still sparse per bin.

**Physics implication.** Under Interpretation B, the tensor observable becomes
a moderate-x measurement rather than a small-x one. The combined significance
is still dominated by moderate-x, moderate-Q² bins where A_bin and N_bin are
both non-negligible, but the per-bin failure-mode diagnostic has shifted.

---

## 7. Physics discussion: two interpretations

### 7.1 What dropping F₁ from Δ means

Under Interpretation A, F₁ appears in the Δ formula (Eq. 1A) and also in the
denominator of the asymmetry amplitude A_cos2φ (Eq. 4, unchanged). The numerator
factor F₁ in Δ effectively amplifies signal in bins where F₁ is large — which
are the low-x, low-Q² bins that also have the highest event rate. This double
alignment (high signal × high rate) produces the small L_5σ values of
Interpretation A.

Under Interpretation B, this alignment is absent. Δ has a pure shape x^α(1-x)^β;
A_cos2φ,bin is still divided by a combination of F₁ and F₂ (Eq. 4), so signal
is suppressed in high-F₁ bins (low x) relative to Interpretation A. Result:
Interpretation B is more conservative by factors of 2–8 in reach.

### 7.2 Sum-rule compensation is incomplete

The analytic Beta-function solver (Eq. 3B) sets A to enforce the first-moment
sum rule by design. However, the sum-rule constraint fixes only the integrated
first moment; it cannot restore the per-bin signal boost that F₁ provided in
Interpretation A. The denominator integral changes when F₁ is dropped, so A
adjusts, but A_cos2φ,bin (Eq. 4) has F₁ in the denominator regardless of
interpretation. The net effect is that |A| under Interpretation B is 3–3.5×
smaller than under A, and the preferentially-boosted bins (low x) are no longer
preferentially boosted.

### 7.3 Which interpretation is physically appropriate?

Both are shape-modeling ansätze motivated by unpolarized-PDF templates. The
`x^α(1-x)^β` template describes the x-shape of a parton distribution (or ratio)
without predicting the overall normalization. The question is whether the template
applies to Δ itself (Interpretation B) or to the ratio Δ/F₁ (Interpretation A).

- *Argument for A:* Δ and F₁ share the same quark content at leading order.
  A ratio ansatz `Δ/F₁ ∝ x^α(1-x)^β` is analogous to the ratio `g₁/F₁` used
  in spin-structure modeling, where the asymmetry has a simpler shape than the
  individual structure functions.

- *Argument for B:* The x^α(1-x)^β template is most naturally applied to the
  structure function itself (not a ratio), matching the standard PDF fit ansatz.
  The form is also simpler, and the resulting A being config-independent is
  physically cleaner.

No current calculation determines which is correct. This is a genuine modeling
ambiguity that only a first-principles computation of the x-shape of Δ can
resolve.

---

## 8. Key physics conclusions

1. **Interpretation choice matters quantitatively (factor 2–8) but not
   qualitatively.** Under both A and B, the physics case is discovery in hours
   to days at 10 fb⁻¹ and systematics-limited within a year. The overall
   conclusion — that ⁶Li tensor gluonometry is feasible at the EIC — is robust.

2. **Interpretation B is the more conservative reach estimate.** Dropping F₁
   from Δ removes the accidental double-alignment of signal and rate at low x.
   Honest reach quotes should span both interpretations.

3. **A is truly config-independent under Interpretation B** — a parametric
   simplification. The analytic Beta-function solver produces a single A value
   per (variant, theory input) pair, with no residual F₁ evolution systematic.

4. **The x-region carrying sensitivity shifts.** Interpretation A: the worst
   per-bin failure mode and the dominant sensitivity are at large x (high-rate
   bins boosted by F₁). Interpretation B: worst bins shift to small x
   (moderate-Q², sparse-per-bin); combined sensitivity is carried by moderate-x
   bins where both shape and rate are non-trivial.

5. **|A_bag|/|A_lat| = 4/3 exactly** under both interpretations — fixed by the
   c-coefficient ratio (0.012/0.009) regardless of which sum-rule solver is used.

6. **The modeling ambiguity is the leading shape systematic.** Today's
   Interpretation A vs B comparison (factor 2–8) exceeds the within-interpretation
   shape systematic (factor ~2 across low_x/mid_x/high_x variants within each
   interpretation). Resolving the ansatz question is now a priority.

---

## 9. Assumptions and caveats

All assumptions from `money_delta_note_2026-07-24.md` §6 carry over unchanged
(bag model origin of Δ, gluonic probe language, status of existing Δ
calculations, small-modulation approximation, uniform φ acceptance). Two new
caveats specific to Interpretation B:

- **F₁ absence from Δ numerator.** Interpretation B does not simply remove a
  correction factor — it changes the physics model for the x-shape of Δ.
  Whether this is more or less physical than Interpretation A is unresolved.

- **Config independence as a diagnostic.** A being config-independent under B
  is analytically exact (α_s cancels, F₁ drops out). Under A, the ~5% residual
  config dependence was from F₁ evolution. The absence of that variation in B
  does not mean B is more accurate — it means one source of variation has been
  modeled away.

---

## 10. Updated next-steps ranking

Cross-reference `money_delta_note_2026-07-24.md` §8 for full Tier A/B/C
structure. The tier structure and existing items are unchanged. New and updated
items:

### Tier A additions

- **Tier A #1 (P_zz convention)** — unchanged; still the highest-priority
  physics input. See 2026-07-24 §8 for full discussion.

- **Tier A #2 (Δ shape anchoring) — sharpened.** Today's comparison shows that
  the ansatz choice (A vs B) alone changes reach by factors of 2–8 — comparable
  to or exceeding the within-ansatz shape systematic. A first-principles
  calculation of the x-shape of Δ would remove this leading ambiguity. Priority
  elevated.

- **Tier A #3 (new): Consult theorists on which ansatz is physically appropriate
  for Δ.** Specifically: is the `x^α(1-x)^β` template better applied to Δ
  directly (Interpretation B) or to Δ/F₁ (Interpretation A)? This is a
  qualitative question about operator structure that experienced tensor-structure
  theorists (Hoodbhoy, Jaffe, Cloet, or the Detmold–Shanahan lattice group)
  could likely answer quickly.

- **Tier A #4 (new): Honest reach quotes should span both interpretations.**
  Publication-ready projections should present Interpretation A and B reach tables
  side-by-side and characterize the factor 2–8 as a physics-modeling systematic,
  not a numerical uncertainty.

### Tiers B and C

Unchanged from `money_delta_note_2026-07-24.md` §8.

---

## 11. Reproducibility

All commands run from the repo root
(`/Users/L00338853/work/Polarized_Li/PolarizedLithiumSim`).

```bash
python3 fastsim/scripts/money_delta_20260725.py --outdir out/money_delta
```

Prerequisites (PDF sets, NumPy patch, `_safe_xfx` helper) are identical to
those documented in `money_delta_note_2026-07-24.md` §9. The script is a clone
of `money_delta_20260724.py` with the single physics change in Eqs. 1 and 3;
all infrastructure (grid backend, `parton` α_s table, P_zz = 0.267) is
unchanged.

---

## 12. File inventory

### Script

- `fastsim/scripts/money_delta_20260725.py` (~1400 lines)

### Plots in `fastsim/out/money_delta/` (new as of 2026-07-27, 8 files)

```
money_delta_20260725_p1_bag_10fb.png       (δA/A vs variant, L=10 fb⁻¹, bag)
money_delta_20260725_p1_bag_100fb.png      (δA/A vs variant, L=100 fb⁻¹, bag)
money_delta_20260725_p2_bag.png            (L_5σ vs variant, bag, 3 configs)
money_delta_20260725_p3_baglat.png         (δA/A, bag + lattice overlaid)
money_delta_20260725_p4_baglat.png         (L_5σ, bag + lattice overlaid)
money_delta_20260725_perbin_mid_10fb.png   (per-bin heatmap, MID, L=10 fb⁻¹, 3 subpanels)
money_delta_20260725_perbin_mid_100fb.png  (per-bin heatmap, MID, L=100 fb⁻¹, 3 subpanels)
money_delta_20260725_perbin_top_10fb.png   (per-bin heatmap, TOP, L=10 fb⁻¹, 3 subpanels)
money_delta_20260725_perbin_top_100fb.png  (per-bin heatmap, TOP, L=100 fb⁻¹, 3 subpanels)
```

### Notes

- `fastsim/notes/money_delta_note_2026-07-24.md` — Interpretation A; full
  derivation of formalism, statistical framework, and physics discussions.
- `fastsim/notes/money_delta_note_2026-07-27.md` (this file)

---

## 13. Bottom line

Interpretation B (Δ = s · α_s · x^α(1-x)^β, no F₁ factor) is the physically
simpler ansatz: Δ has a valence-like shape directly, α_s cancels analytically in
the sum rule, and A is fully config-independent. It gives more conservative reach
by factors of 2–8 compared to Interpretation A. The choice between A and B is a
genuine physics-modeling ambiguity that no current data or calculation resolves;
under both interpretations, ⁶Li tensor gluonometry at the EIC is discoverable in
hours to days at 10 fb⁻¹/nucleon. Honest reach projections must span both
interpretations, and the next priority is to consult theorists on which ansatz is
physically appropriate for Δ.
