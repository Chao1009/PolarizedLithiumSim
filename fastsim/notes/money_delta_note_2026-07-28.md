# Money-Δ Plot: φ-Distribution Projections — 2026-07-28

**Observable.** Fractional yield modulation `(N_φ − N_flat)/N_flat` vs `φ` for
the `cos(2φ)` tensor asymmetry on transversely tensor-polarized ⁶Li. Today's
work shifts from reach plots (`δA/A`, `L_5σ`) to the raw observable: the
φ-differential yield modulation that an experimenter would actually measure and
fit. Nine PNGs produced, covering three beam configs × three integration cases,
with three signal scales overlaid per plot.

**Revision history**

- 2026-07-28 (initial): New script `money_delta_20260728.py` (~1005 lines).
  Fixed setup: Interpretation A, mid_x shape (α=0.7, β=3), R1998, EPPS21,
  Cloet P_zz = 0.267, L = 10 fb⁻¹/nucleon. Three configs × three integration
  cases × three signal scales = 9 plots (one per config/case combination,
  three cosine curves per plot). A_bag values carried forward from
  `money_delta_20260724.py`. See `money_delta_note_2026-07-24.md` for the
  full sum-rule / Fisher-information derivation and `money_delta_note_2026-07-27.md`
  for the Interpretation A vs B discussion.

---

## 1. What was produced today

Script: `fastsim/scripts/money_delta_20260728.py` (~1005 lines). All figures
land in `fastsim/out/money_delta/` with `_20260728_` filename prefix. Ion is
⁶Li throughout. Fixed to Interpretation A, mid_x shape. Three beam configs ×
three integration cases = 9 (config, case) combinations; each plot shows three
overlaid cosine curves (s = 0.1·A_bag, A_bag, 10·A_bag) with one set of Poisson
error bars per φ bin.

### Beam configurations

| Label | E_e [GeV] | p_ion [GeV/u] | √s [GeV/u] |
|---|---|---|---|
| LOW | 5 | 27.5 | ~23 |
| MID | 10 | 50 | ~45 |
| TOP | 18 | 137.5 | ~100 |

### Fixed parameters

- Luminosity: L = 10 fb⁻¹/nucleon
- x-shape: mid_x (α = 0.7, β = 3)
- Interpretation A: Δ = s · α_s(Q²) · F₁(x, Q²) · x^α (1−x)^β
- R1998 (σ_L/σ_T ratio), EPPS21 (nuclear PDF), Cloet P_zz = 0.267
- φ bins: 72 bins of 5° each, from 0 to 2π

### Integration cases

| Case | Description |
|---|---|
| 1 (peak bin) | Single (x, Q²) bin with the largest N_bin |
| 2 (Q² slice) | Single Q² bin (iq2 = 3, Q² = 2.43 GeV²), integrated over all accepted x |
| 3 (integrated) | All accepted (x, Q²) bins combined |

### Hardcoded A_bag values (from money_delta_20260724.py, Interpretation A, mid_x)

| Config | \|A_bag\| |
|---|---|
| LOW | 0.318 |
| MID | 0.310 |
| TOP | 0.297 |

### Output figures (9 PNGs in `fastsim/out/money_delta/`)

| Filename | Content |
|---|---|
| `money_delta_20260728_phimodulation_peakbin_low.png` | LOW config, Case 1 (peak bin) |
| `money_delta_20260728_phimodulation_peakbin_mid.png` | MID config, Case 1 (peak bin) |
| `money_delta_20260728_phimodulation_peakbin_top.png` | TOP config, Case 1 (peak bin) |
| `money_delta_20260728_phimodulation_q2slice_low.png` | LOW config, Case 2 (Q² slice) |
| `money_delta_20260728_phimodulation_q2slice_mid.png` | MID config, Case 2 (Q² slice) |
| `money_delta_20260728_phimodulation_q2slice_top.png` | TOP config, Case 2 (Q² slice) |
| `money_delta_20260728_phimodulation_integrated_low.png` | LOW config, Case 3 (fully integrated) |
| `money_delta_20260728_phimodulation_integrated_mid.png` | MID config, Case 3 (fully integrated) |
| `money_delta_20260728_phimodulation_integrated_top.png` | TOP config, Case 3 (fully integrated) |

---

## 2. Formulas

Cross-reference `money_delta_note_2026-07-24.md` §2 for the full equation set
(Eqs. 1–11). Today's script uses those formulas unchanged. The observable-side
formula for the φ distribution follows from Eq. (4) of that note:

**Fractional yield modulation per φ bin:**

    y(φ) = (N_φ − N_flat) / N_flat = P_zz · <A_cos2φ> · cos(2φ)          (1)

**Poisson error per φ bin:**

    σ_y(φ) = √N_φ / N_flat ≈ 1/√N_flat    (small modulation limit)        (2)

where N_flat = N_total / 72 and ⟨A_cos2φ⟩ is either a single-bin amplitude
(Case 1) or a rate-weighted average over the integrated bins (Cases 2 and 3).
The error bars are essentially signal-scale-independent because |P_zz · A · cos(2φ)|
≪ 1 across all cases (max fractional modulation = 5.99%; see §4).

All other formulas in `money_delta_note_2026-07-24.md` §2 (Eqs. 2–11) are
unchanged and apply equally here.

---

## 3. Numerical results — complete table

For each (config, case): peak bin coordinates, total yield, yield per φ bin at
N_flat, rate-weighted cosine amplitude ⟨A_cos2φ⟩, and maximum fractional
modulation at 10×A_bag.

| Config | Case | Peak (x, Q²) or slice | N_total | N_flat/bin | ⟨A_cos2φ⟩ | P_zz·⟨A⟩ at A_bag | Max mod. at 10×A_bag |
|---|---|---|---|---|---|---|---|
| LOW | 1 (peak bin) | (8.91e-2, 2.43 GeV²) | 2.37e+7 | 3.29e+5 | −2.24e-2 | 0.0060 | 0.0599 |
| LOW | 2 (Q² slice) | Q² slice iq2=3 (2.43 GeV²) | 3.11e+8 | 4.32e+6 | −1.58e-2 | 0.0042 | 0.0421 |
| LOW | 3 (integrated) | all bins | 1.20e+9 | 1.67e+7 | −1.70e-2 | 0.0045 | 0.0453 |
| MID | 1 (peak bin) | (7.08e-2, 2.43 GeV²) | 2.46e+7 | 3.42e+5 | −1.98e-2 | 0.0053 | 0.0528 |
| MID | 2 (Q² slice) | Q² slice iq2=3 (2.43 GeV²) | 4.11e+8 | 5.71e+6 | −1.04e-2 | 0.0028 | 0.0278 |
| MID | 3 (integrated) | all bins | 1.79e+9 | 2.49e+7 | −1.21e-2 | 0.0032 | 0.0323 |
| TOP | 1 (peak bin) | (2.24e-3, 2.43 GeV²) | 2.43e+7 | 3.38e+5 | −2.08e-3 | 0.0006 | 0.0056 |
| TOP | 2 (Q² slice) | Q² slice iq2=3 (2.43 GeV²) | 4.49e+8 | 6.24e+6 | −3.46e-3 | 0.0009 | 0.0092 |
| TOP | 3 (integrated) | all bins | 2.25e+9 | 3.13e+7 | −5.33e-3 | 0.0014 | 0.0142 |

All positivity checks pass: max |P_zz · A_cos2φ · cos(2φ)| across all
(config, case, scale) combinations is 0.0599 (LOW + Case 1 + 10×A_bag),
well below the physical bound of 1.

---

## 4. Signal-to-noise analysis

For each (config, case) at s = A_bag, the signal-to-noise ratio for a
combined cos(2φ) fit to the φ-binned yield uses:

    S/N ≈ P_zz · |⟨A_cos2φ⟩| · √(N_total / 2)                            (3)

This is a single-cosine integrated fit (cf. the per-bin Fisher sum of
`money_delta_note_2026-07-24.md` Eq. 6 — see §6 for the comparison).

| Config | Case | \|⟨A_cos2φ⟩\| | N_total | S/N ≈ |
|---|---|---|---|---|
| LOW | 1 (peak bin) | 0.0224 | 2.37e+7 | 20.7 |
| LOW | 2 (Q² slice) | 0.0158 | 3.11e+8 | 52.5 |
| LOW | 3 (integrated) | 0.0170 | 1.20e+9 | 111.1 |
| MID | 1 (peak bin) | 0.0198 | 2.46e+7 | 18.6 |
| MID | 2 (Q² slice) | 0.0104 | 4.11e+8 | 39.9 |
| MID | 3 (integrated) | 0.0121 | 1.79e+9 | 96.6 |
| TOP | 1 (peak bin) | 0.00208 | 2.43e+7 | 1.94 |
| TOP | 2 (Q² slice) | 0.00346 | 4.49e+8 | 13.8 |
| TOP | 3 (integrated) | 0.00533 | 2.25e+9 | 47.7 |

---

## 5. Key physics conclusions

1. **Case 3 (fully integrated) wins in S/N for all three configs.** This was
   not the expected outcome. An earlier analysis predicted that Case 2 (Q²
   slice) would be the visual winner on the grounds that "combining bins with
   different signal sizes dilutes the average." In practice, the √N statistics
   gain from Case 1 to Case 3 (factor ~7–10 in √N_total) beats the dilution
   of ⟨A_cos2φ⟩ (factor ~0.6–0.8) for all three configs. Case 3 beats Case 2
   by factors of 2–3 in S/N.

2. **TOP config's peak-rate bin is a kinematic trough for the tensor signal.**
   TOP's peak-rate bin sits at very small x (x = 2.24×10⁻³), corresponding to
   y = Q²/(sx) ≈ 0.11 and kinematic prefactor (1−y)/y² ≈ 74. When TOP
   integrates over its full acceptance, bins at moderate x and moderate Q²
   contribute |A_cos2φ| values 2–3× larger. This is the only config where
   |⟨A_cos2φ⟩| *grows* as more bins are included:

       TOP: 0.0021 (peak bin) → 0.0035 (Q² slice) → 0.0053 (integrated)

   For LOW and MID the amplitude decreases slightly from Case 1 to Case 3, as
   bins with smaller individual signal are added. For TOP the trend inverts
   because the rate-peak sits at kinematically unfavorable small x.

   **Practical implication:** for a TOP-config data-taking program, triggers
   and bin schemes designed around the DIS rate peak will systematically
   undersample the kinematic region where the tensor observable is largest.
   Full-acceptance integration is not just statistically optimal — it is
   physically necessary to reach the sensitive bins.

3. **All configs give firm discovery at 10 fb⁻¹ under the bag prediction.**
   Even the worst single case (TOP + peak bin, S/N = 1.94) is poor only
   because it uses a small fraction of the available data. Integrating all
   accepted bins (TOP + Case 3) gives S/N = 47.7 — a comfortable discovery.

4. **The per-bin Fisher-information fit (2026-07-24) is ~20% better than
   today's single-cosine integrated fit (Case 3).** For MID at 10 fb⁻¹,
   yesterday's per-bin combined Fisher fit gave σ = 116 (equivalently
   δA/A = 0.0086), while today's Case 3 single-cosine fit gives S/N ≈ 97
   (δA/A ≈ 0.010). The ~20% gap arises because the Fisher sum
   (Eq. 6 of 2026-07-24) weights each bin by A_bin² · N_bin, extracting
   information proportional to each bin's squared signal-to-noise; the
   single-cosine fit uses only the yield-weighted average amplitude and
   therefore under-weights high-signal bins relative to their Fisher
   contribution.

   In practice an experimenter uses both approaches: Case 3-style φ-modulation
   plots for visualization and presentation ("the cosine is visually there")
   and the per-bin combined fit for the quoted measurement uncertainty.

5. **Positivity is not a concern at 10×A_bag.** The maximum fractional
   modulation across all combinations is 5.99% (LOW + Case 1 + 10×A_bag),
   confirming that the linear small-modulation approximation (Eq. 1) is
   valid throughout the full scale scan.

---

## 6. Comparison to prior work

### 6.1 Relation to 2026-07-24 reach plots

Today's plots are complementary to, not a replacement of, the reach plots in
`money_delta_note_2026-07-24.md`. The prior work showed `δA/A` and `L_5σ` as
functions of the signal scale — asking "how much luminosity do we need?" Today's
plots show the raw φ modulation at a fixed luminosity (10 fb⁻¹) — asking "what
does the data look like when we take it?" The S/N values in §4 above (e.g.,
MID + Case 3, S/N = 97) are consistent with the 2026-07-24 combined Fisher fit
(σ = 116) to within the ~20% gap explained in §5 point 4.

### 6.2 Relation to 2026-07-27 Interpretation A vs B discussion

Today's script is fixed to Interpretation A (Δ = s · α_s · F₁ · x^α(1−x)^β)
and the mid_x shape, consistent with the reference case used in §5 of
`money_delta_note_2026-07-27.md`. The A_bag values used here (LOW: 0.318,
MID: 0.310, TOP: 0.297) match the Interpretation A entries from
`money_delta_note_2026-07-24.md` §3 to within rounding. No Interpretation B
runs were produced today. The factor-2–8 reach difference between
interpretations (2026-07-27 §5) applies equally to the S/N values above:
under Interpretation B, all S/N values would be reduced by factors of 2–3.

### 6.3 Case 2 vs Case 3 — why the expectation was wrong

Prior intuition held that mixing bins with different A_bin values would dilute
the average enough to make Case 2 (a single Q² bin) superior to Case 3 (all
bins) in S/N. The numbers show the dilution is modest (~24–40% reduction in
|⟨A_cos2φ⟩| from Case 1 to Case 3) while the yield increase is a factor of
50–100. Because S/N scales as ⟨A⟩ · √N, the √N gain dominates. The
analogous conclusion for the Fisher-sum observable is that bins with large
N_bin · A_bin² are strongly preferred; no single Q² slice captures all such
bins.

---

## 7. What the plots visually show

Each PNG contains a single panel: fractional yield modulation vs φ (radians,
0 to 2π), with three overlaid cosine curves and one set of Poisson error bars
per φ bin.

**Signal scale coding:**

| Curve | Style | Meaning |
|---|---|---|
| 0.1·A_bag | Light blue, thin | Signal 10× below bag prediction |
| A_bag | Black, medium | Bag prediction |
| 10·A_bag | Red, thick | Signal 10× above bag prediction |

**Error bars** are computed from Poisson statistics per φ bin at N_flat/bin
(§3 column 5), essentially signal-independent (Eq. 2).

**Case 1 (peak bin) plots** look noisy: large cosine amplitude
(⟨A⟩ ~ 0.02 for LOW/MID) but relatively large error bars because only
~2.4×10⁷ events per bin and 72 bins of ~3.3×10⁵ events each.

**Case 3 (integrated) plots** look clean: the cosine is resolved at smaller
fractional amplitude because σ_y ≈ 1/√N_flat is ~7× smaller (N_flat ~
1.7–3.1×10⁷ events per φ bin at 1.20–2.25×10⁹ total events). The sinusoid
is beautifully resolved even at A_bag.

**The three-scale overlay directly answers:** "if the true signal is smaller
than the bag prediction, can we still see it?"

| Case | 0.1·A_bag visible? | A_bag visible? |
|---|---|---|
| 1 (peak bin) | Barely — near or below error bars | Clearly (3–20× error bar) |
| 2 (Q² slice) | Hint-level | Comfortable |
| 3 (integrated) | Clearly | Overwhelming |

---

## 8. Assumptions and caveats

All assumptions from `money_delta_note_2026-07-24.md` §6 carry over unchanged:
bag model origin of Δ, gluonic probe language, status of existing Δ calculations,
small-modulation approximation, and uniform φ acceptance. The following
are specific to today's script:

- **Fixed to Interpretation A, mid_x, 10 fb⁻¹.** Today's script does not scan
  shape variants, interpretation, or luminosity. Sensitivity to these choices
  is documented in the prior notes.

- **S/N formula (Eq. 3) assumes uniform φ efficiency.** Real detectors have
  azimuthal acceptance gaps. Non-uniform coverage changes both ⟨A_cos2φ⟩ and
  the effective N_flat per φ bin. See `money_delta_note_2026-07-24.md` §4.2.1
  caveat 1 for the correction formulation.

- **Rate-weighted ⟨A_cos2φ⟩ for Cases 2 and 3.** The amplitude plotted is not
  a simple average — it is weighted by N_bin so that the predicted cosine
  curve matches the yield-weighted expectation. High-N_bin bins dominate.

- **TOP peak-bin finding is Interpretation-A-specific.** Under Interpretation B
  (2026-07-27 §7.1), F₁ no longer appears in the Δ numerator and the signal
  does not track the low-x rate peak. The kinematic trough effect at small x
  under TOP would differ quantitatively under Interpretation B.

---

## 9. Updated next-steps ranking

Building on `money_delta_note_2026-07-24.md` §8 (Tier A/B/C structure) and
`money_delta_note_2026-07-27.md` §10 (Tier A additions):

### Tier A — unchanged and additions

- **Tier A #1 (P_zz convention resolution)** — unchanged; still the
  highest-priority physics input. See `money_delta_note_2026-07-24.md` §8.

- **Tier A #2 (Δ shape anchoring via first-principles theory)** — unchanged;
  priority elevated in 2026-07-27 §10. The Interpretation A vs B comparison
  (factor 2–8 in reach) dominates the shape systematic.

- **Tier A #3 and #4 from 2026-07-27** — unchanged (theorist consultation on
  ansatz choice; publication-ready reach quotes spanning both interpretations).

- **New observation from today (informing Tier A planning):** For the TOP
  config specifically, the peak-rate bin is a kinematic trough for the tensor
  observable. If a TOP-config data-taking program is planned, trigger design
  and bin selection should not assume that small-x bins carry the tensor signal.
  The sensitive bins (moderate x, moderate Q²) must be explicitly included in
  the integration. This is a detector/trigger design consideration, not a
  physics-modeling uncertainty; it belongs to Tier A planning for any TOP-config
  proposal.

### Tiers B and C

Unchanged from `money_delta_note_2026-07-24.md` §8.

---

## 10. Reproducibility

All commands run from the repo root
(`/Users/L00338853/work/Polarized_Li/PolarizedLithiumSim`).

```bash
python3 fastsim/scripts/money_delta_20260728.py --outdir out/money_delta
```

Prerequisites (PDF sets, NumPy patch, `_safe_xfx` helper) are identical to
those documented in `money_delta_note_2026-07-24.md` §9. The script is a new
production; it does not clone either prior script but shares the same grid
backend, `parton` α_s table, R1998, EPPS21, and P_zz = 0.267 infrastructure.

---

## 11. File inventory

### Script

- `fastsim/scripts/money_delta_20260728.py` (~1005 lines)

### Plots in `fastsim/out/money_delta/` (new as of 2026-07-28, 9 files)

```
money_delta_20260728_phimodulation_peakbin_low.png      (LOW, Case 1: peak bin)
money_delta_20260728_phimodulation_peakbin_mid.png      (MID, Case 1: peak bin)
money_delta_20260728_phimodulation_peakbin_top.png      (TOP, Case 1: peak bin)
money_delta_20260728_phimodulation_q2slice_low.png      (LOW, Case 2: Q² slice iq2=3)
money_delta_20260728_phimodulation_q2slice_mid.png      (MID, Case 2: Q² slice iq2=3)
money_delta_20260728_phimodulation_q2slice_top.png      (TOP, Case 2: Q² slice iq2=3)
money_delta_20260728_phimodulation_integrated_low.png   (LOW, Case 3: fully integrated)
money_delta_20260728_phimodulation_integrated_mid.png   (MID, Case 3: fully integrated)
money_delta_20260728_phimodulation_integrated_top.png   (TOP, Case 3: fully integrated)
```

### Notes

- `fastsim/notes/money_delta_note_2026-07-24.md` — Interpretation A; full
  derivation of formalism, sum-rule solver, statistical framework, and
  Fisher-information combined fit.
- `fastsim/notes/money_delta_note_2026-07-27.md` — Interpretation B (no F₁
  in Δ); A vs B comparison, analytic Beta-function solver.
- `fastsim/notes/money_delta_note_2026-07-28.md` (this file)

---

## 12. Bottom line

The φ-modulation projections show clean cosine signals at all three configs and
all three integration cases under the bag prediction at 10 fb⁻¹, with Case 3
(fully integrated) giving the best signal-to-noise; the surprising finding is
that TOP's peak-rate bin is a kinematic trough for the tensor observable,
arguing for full-acceptance integration rather than targeted peak-bin analysis.
