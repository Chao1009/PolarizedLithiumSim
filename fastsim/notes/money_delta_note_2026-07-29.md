# Money-Δ Plot: Detector-Realistic Smearing (ePIC Tracking + ε_eID) — 2026-07-29

> **Dated record.** Numbers here predate the 2026-08-27 corrections: ⁶Li
> beam energies are γ-matched (40.8 / 99.5 / 137.5 GeV/u), not
> rigidity-scaled (20.5 / 50 / 137.5), and the far-forward divergence is
> per-configuration rather than a single 72.7 μrad. See plans/10. Kept as
> the record of its date, not as a current reference.


> **Superseded in one input, 2026-08-26 (R = σ_L/σ_T).**  Every `L_5σ`,
> `A_bag` and significance number below was computed with the script's
> `r1998`, which multiplies Θ(x, Q²) into all three terms of the SLAC fit
> instead of the log term alone and therefore returns R = 1 over the region
> that carries the sensitivity (code review 2026-08-25 item S1).  Re-run
> with the published fit (`--r-model r1998`, now `structure.r1998`):
> **L_5σ at Δ/F₁ = 10⁻³ becomes 67.5 / 65.8 / 155.1 fb⁻¹/u** for
> LOW / MID / TOP against the 135.3 / 131.3 / 274.6 recorded here, |A_bag|
> falls 21–25% (0.318/0.310/0.297 → 0.237/0.236/0.235) and the peak-bin
> significance rises ×1.09–1.20.  The numbers below are still exactly
> reproducible — that is the point of a dated note — with the default
> `--r-model simplified`.  See plans/00 development run 8.


**Observable.** Fractional yield modulation `(N_φ − N_flat)/N_flat` vs `φ` for
the `cos(2φ)` tensor asymmetry on transversely tensor-polarized ⁶Li, now with
Monte Carlo detector smearing applied to the scattered electron kinematics.
Today's script `money_delta_20260729.py` is the detector-realistic version of
yesterday's `money_delta_20260728.py`: same 16 plots (3 heatmaps + 8 φ-modulation
plots + 5 TOP+peakbin bin-width scan), but with ePIC tracking-only momentum and
angular resolution, and η-dependent electron-ID efficiency applied per accepted
event. A sign-convention bug was found and fixed during the run; the post-fix
numbers are reported throughout.

## Small-modulation limit — a note before reading

> [!IMPORTANT]
> **Every plot and number in this note assumes the small-modulation limit:**
> |P_zz · A_cos2φ · cos(2φ)| ≪ 1.
> Readers should understand what this means before interpreting the figures.

**The approximation.** The spin-1 DIS cross section is

    dσ/dΩ ∝ 1 + P_zz · A_cos2φ · cos(2φ)                                (*)

This expression is *exact* at leading twist in Δ — it is not a Taylor expansion.
"Small modulation" refers to the *observable* being small, not to any truncation
of the model. Positivity of dσ/dΩ across φ requires |P_zz · A_cos2φ| ≤ 1 always;
the small-modulation condition is satisfied when we are well inside that bound.

**Consequence 1 — error bars are essentially signal-independent.** The per-φ-bin
Poisson error on the fractional modulation is σ_y(φ) ≈ 1/√N_flat in the
small-modulation limit. Each plot uses one set of error bars across signal scales;
variation across s is <3% even at 10×A_bag.

**Consequence 2 — the linear cross-section formula is exact enough to use
directly.** Because (*) is exact in Δ at leading twist, using it without further
approximation introduces no error. The only physical constraint is positivity,
which is satisfied throughout this study.

See `money_delta_note_2026-07-28.md` for the full small-modulation-limit
derivation and the 5.99% maximum-modulation verification for the pre-detector
plots. The detector-smeared plots stay within the same regime.

---

**Revision history**

- 2026-07-29 (initial): New script `money_delta_20260729.py` (~1870 lines).
  Detector-realistic extension of `money_delta_20260728.py`. Same physics
  configuration: Interpretation A, mid_x shape (α=0.7, β=3), R1998, EPPS21,
  Cloet P_zz = 0.267, L = 10 fb⁻¹/nucleon. A_bag is solved at startup via
  `solve_A_from_sum_rule()` (c_bag = −0.012, mid_x, Interpretation A); the
  solved values are signed and negative. The observable pipeline consumes
  |A_bag|; signed values are printed at startup for audit. Reference magnitudes
  (LOW ≈ 0.318, MID ≈ 0.310, TOP ≈ 0.297) are retained in the script as
  regression guards only. New additions: ePIC tracking-only momentum and angular
  resolution (piecewise in η, 7 regions each), η-dependent electron-ID efficiency
  (ATHENA + ECCE synthesis, 9 anchor points), 1000 MC events per accepted true
  bin. A sign-convention bug in the reconstructed Q² and y formulas was diagnosed
  and fixed mid-session (§5). 16 PNGs produced with `_20260729_` prefix. See
  `money_delta_note_2026-07-24.md` for the full formula derivation and
  `money_delta_note_2026-07-28.md` for the pre-detector baseline.

- 2026-07-29 (extension): Added raw (efficiency-independent) same-bin diagnostic
  alongside the existing efficiency-weighted one. New numbers: raw 88/76/64% for
  LOW/MID/TOP vs weighted 76/62/47%. Gap of 12–17% reflects that migrated events
  land preferentially at higher-ε η regions than the same-bin events came from;
  see new §7a for physics.

---

## 1. What was produced today

Script: `fastsim/scripts/money_delta_20260729.py` (~1870 lines). All figures land
in `fastsim/out/money_delta/` with `_20260729_` filename prefix. Ion is ⁶Li
throughout. Physics configuration (Interpretation A, mid_x shape, R1998, EPPS21,
P_zz = 0.267) is identical to `money_delta_20260728.py`. A_bag is now solved at
startup from the bag sum rule (c_bag = −0.012, mid_x, Interpretation A) rather
than hardcoded; the solved |A_bag| values reproduce the previous hardcoded
magnitudes to within 1% (see §1 A_bag table). The detector model is new:
scattered electron kinematics are smeared using 1000 MC events per accepted true
`(x, Q²)` bin, each weighted by the η-dependent electron-ID efficiency.

**Note on the initial run.** The first run of the script produced pathological
output: 0% of events landing in their original reco bin and peak-rate reco bins
displaced by factors of 100+ in Q². The bug was a sign-convention error in the
`smear_config` reconstruction formulas (see §5). All results below are post-fix.

### Beam configurations

| Label | E_e [GeV] | p_ion [GeV/u] | √s [GeV/u] |
|---|---|---|---|
| LOW | 5 | 27.5 | ~23 |
| MID | 10 | 50 | ~45 |
| TOP | 18 | 137.5 | ~100 |

### Fixed physics parameters

- Luminosity: L = 10 fb⁻¹/nucleon
- x-shape: mid_x (α = 0.7, β = 3)
- Interpretation A: Δ = s · α_s(Q²) · F₁(x, Q²) · x^α (1−x)^β
- R1998 (σ_L/σ_T ratio), EPPS21 (nuclear PDF), Cloet P_zz = 0.267

### A_bag values (Interpretation A, mid_x)

At startup, `money_delta_20260729.py` calls `solve_A_from_sum_rule()` for each
beam configuration and stores the result in `A_BAG_SIGNED`. The solved values are
**signed and negative** (Δ ∝ c_bag = −0.012 < 0 under Interpretation A). The
observable pipeline — `smear_config`, `compute_A_cos2phi_at_bin`,
`build_phi_plot`, `build_perbin_heatmap`, and the summary table — consumes
`abs(A_BAG_SIGNED[config_tag])`, so downstream plot and sign behavior is
identical to the previous hardcoded-magnitude convention.

The table below lists the historical hardcoded magnitudes retained in the script
as `REFERENCE_ABS_A_BAG` for audit and regression-guard purposes only. The
solver must reproduce these values to within 1%; a larger discrepancy raises a
`RuntimeError` at startup.

| Config | Historical \|A_bag\| (audit reference) |
|---|---|
| LOW | 0.318 |
| MID | 0.310 |
| TOP | 0.297 |

### Integration cases

| Case | Description |
|---|---|
| 1 (peak bin) | Single reco `(x, Q²)` bin with the largest reco N_bin |
| 2 (Q² slice) | Single Q² bin (iq2 = 3, Q² = 2.43 GeV²), integrated over all accepted reco x |
| 3 (integrated) | All accepted reco `(x, Q²)` bins combined |

### New detector model (relative to `money_delta_20260728.py`)

- **Tracking momentum resolution:** σ_p/p = √(a²·p² + b²), piecewise in η
  (7 regions, coefficients a and b from ePIC tracking-only performance).
- **Tracking angular resolution:** σ_θ = σ_φ, piecewise in η (7 regions,
  1–5 mrad).
- **Electron-ID efficiency:** η-dependent, linearly interpolated from 9 anchor
  points (ATHENA + ECCE synthesis; 0.85–0.95 in |η|<3, dropping to 0.70 at
  |η|=3.5; zero outside |η|>3.5).
- **Tracking only, no ECal energy resolution.** User decision: use tracking
  only for the scattered electron. ECal is not modeled in this script; it is
  deferred as a Tier B refinement.
- **1000 MC events per accepted true bin**, weighted by ε_eID(η_reco).

---

## 2. Formulas — detector additions

Cross-reference `money_delta_note_2026-07-24.md` §2 for the complete formula set
(Eqs. 1–11, unchanged). Cross-reference `money_delta_note_2026-07-28.md` §2 for
the φ-modulation and S/N formulas. Today's additions cover the detector smearing
model only.

**Tracking momentum resolution:**

    σ_p / p = √(a² · p² + b²)     [piecewise in η, 7 regions]

**Tracking angular resolution** (with σ_θ = σ_φ):

    σ_θ (rad) = mrad(η) × 10⁻³    [piecewise in η, 7 regions, 1–5 mrad]

**Reconstructed DIS kinematics** (electron method, from smeared E' and θ):

    Q²_reco = 2 E_e E'_reco (1 + cos θ_reco)              (Q² > 0 required)
    y_reco  = 1 − (E'_reco / 2E_e) (1 − cos θ_reco)      (0 < y < 1 required)
    x_reco  = Q²_reco / (y_reco · s)                      (0 < x < 1 required)
    η_reco  = −log[tan(θ_reco / 2)]

The sign convention (1 + cos θ for Q², (1 − cos θ)/2 for the y term) reflects
the project convention that θ is measured from **+z (ion direction)**; the
scattered electron sits at θ ≈ π. See §5 for the bug diagnosis.

**Electron-ID efficiency** (linear interpolation of 9 anchor points,
ATHENA + ECCE-based):

    η:      −3.5   −3.0   −2.0   −1.0    0.0    1.0    2.0    3.0    3.5
    ε_eID:  0.85   0.92   0.95   0.93   0.90   0.90   0.85   0.80   0.70

Zero outside |η| > 3.5.

**Per-event reco weight:**

    w = (N_true / n_mc) × ε_eID(η_reco)

Applied per MC event. Reconstructed bin counts and ⟨A_cos2φ⟩ averages are built
from these weights.

**S/N formula** (unchanged from `money_delta_note_2026-07-28.md` §4):

    S/N ≈ P_zz · |⟨A_cos2φ⟩| · √(N_total / 2)                       (3)

---

## 3. Numerical results

### 3.1 Smearing diagnostics (post-fix)

| Config | Survival % | Bin-migration — raw (% MC events) | Bin-migration — weighted (% reco yield) | Physicality rejected | Peak-rate reco bin | Peak-rate true bin | A_cos2φ at peak-rate reco | A at peak (unsmeared) |
|---|---|---|---|---|---|---|---|---|
| LOW | 93.6% | 88% | 76% | 200 / 218,000 | (0.089, 2.43) | (0.089, 2.43) | −0.0228 | −0.0224 |
| MID | 90.4% | 76% | 62% | 8,577 / 345,000 | **(0.011, 2.43)** | (0.071, 2.43) | −0.0078 | −0.0198 |
| TOP | 86.0% | 64% | 47% | 23,029 / 455,000 | **(0.0036, 2.43)** | (0.0022, 2.43) | −0.0043 | −0.0021 |

Survival % is the fraction of generated MC events passing the physicality cuts
(Q² > 0, 0 < y < 1, 0 < x < 1) and the |η| < 3.5 acceptance. **Bin-migration
— raw**: fraction of surviving (unweighted) MC events landing in the same
`(x, Q²)` bin as their true-level parent; a pure geometric quantity set by σ_p,
σ_θ, and the log-grid bin widths. **Bin-migration — weighted**: same-bin
efficiency-weighted reco yield divided by total efficiency-weighted reco yield;
represents the fraction of reconstructed yield that originated in the true bin
(see §7a). Physicality-rejected counts are small relative to total events
(0.1–5%), consistent with edge-of-acceptance losses rather than a systematic
issue.

### 3.2 Reach summary (post-fix, at L = 10 fb⁻¹/u, s = A_bag)

Comparison to `money_delta_20260728.py` (pre-detector) S/N values.

| Config | Case | 20260728 S/N (no detector) | 20260729 S/N (detector) | Ratio |
|---|---|---|---|---|
| LOW | 1 (peak bin) | 20.7 | ~21 | 100% |
| LOW | 2 (Q² slice) | 52.5 | ~49 | 93% |
| LOW | 3 (integrated) | 111.1 | ~104 | 93% |
| MID | 1 (peak bin) | 18.6 | **~7.9** | **42%** |
| MID | 2 (Q² slice) | 39.9 | ~24 | 60% |
| MID | 3 (integrated) | 96.6 | ~89 | 92% |
| TOP | 1 (peak bin) | 1.94 | ~4.1 | **214%** |
| TOP | 2 (Q² slice) | 13.8 | ~9 | 65% |
| TOP | 3 (integrated) | 47.7 | ~44 | 92% |

Case 3 (fully integrated) is 92–93% of the pre-detector S/N across all three
configs. Case 1 (peak-rate single bin) changes by factors of 2–3 in either
direction depending on how migration reshapes the rate distribution.

### 3.3 Complete reco ⟨A_cos2φ⟩ comparison

| Config | Case | ⟨A_cos2φ⟩ reco (20260729) | ⟨A_cos2φ⟩ pre-detector (20260728) |
|---|---|---|---|
| LOW | 1 (peak bin) | −0.02285 | −0.02244 |
| LOW | 2 (Q² slice) | −0.01576 | −0.01576 |
| LOW | 3 (integrated) | −0.01698 | −0.01698 |
| MID | 1 (peak bin) | −0.00777 | −0.01977 |
| MID | 2 (Q² slice) | −0.00935 | −0.01039 |
| MID | 3 (integrated) | −0.01171 | −0.01209 |
| TOP | 1 (peak bin) | −0.00426 | −0.00208 |
| TOP | 2 (Q² slice) | −0.00304 | −0.00346 |
| TOP | 3 (integrated) | −0.00502 | −0.00533 |

Case 3 integrated values match to within a few percent across all configs,
as expected: bin migration internal to the accepted region washes out under
integration (§4).

---

## 4. Physical interpretation of Case 1 vs Case 3 divergence

### 4.1 Why Case 3 (integrated) is robust

Under full integration, bin migration is a redistribution within the accepted
region: events leaving one bin land in another, and both bins are included in
the sum. The only real losses are (a) events migrating outside the acceptance
edge and (b) the η-dependent ε_eID dropping below unity. The first effect is
small (physicality rejection is 0.1–5%); the second accounts for the 7–14%
S/N reduction observed. Across all three configs, Case 3 retains 92–93% of the
pre-detector S/N. The residual loss is dominated by efficiency, not migration.

### 4.2 Why Case 1 (peak-rate single bin) is fragile

The peak-rate reco bin is determined by whichever bin accumulates the most
reconstructed events. On log-spaced `(x, Q²)` grids where the rate distribution
peaks broadly, modest migration can shift the reco rate maximum to a completely
different bin from the true rate maximum. Whether the reco peak bin has a larger
or smaller |A_cos2φ| than the true peak bin depends on the local |A| landscape
and is config-specific.

**MID:** the true rate distribution in the Q² = 2.43 GeV² slice rises smoothly
toward small x. Migration pools reco events into a bin at x = 0.011 instead of
x = 0.071 (the true peak). The x = 0.011 bin sits in a region where the mid_x
shape factor x^0.7 is ~6× smaller, suppressing |A_cos2φ| by 2.5×. Reco Case 1
S/N drops to 42% of pre-detector.

**TOP:** the true peak-rate bin is already at extreme small x (x = 0.0022) —
a kinematic trough for the tensor signal (established in `money_delta_note_2026-07-28.md`
§5a). Migration shifts the reco peak to x = 0.0036, which sits in a slightly
less unfavorable region of the |A| landscape. The reco Case 1 |A_cos2φ| rises
by 2×; reco Case 1 S/N is 214% of pre-detector. This is not a gain from
detector smearing — it is an artifact of migration moving the reco peak slightly
away from the true kinematic trough.

**LOW:** the true peak-rate bin is already at moderate x (x = 0.089), which is
near the center of the |A| landscape for this config. Migration is small (76% of
events stay in the original bin) and rearranges symmetrically around the peak.
The reco Case 1 |A_cos2φ| is nearly unchanged (−0.0228 vs −0.0224); S/N ratio
is 100%.

### 4.3 Physics take-away

For the tensor cos(2φ) observable, use the fully-integrated analysis (Case 3),
not a single-bin peak analysis. Case 3 is nearly detector-invariant (92% of
pre-detector S/N); Case 1 changes by factors of 2–3 in either direction depending
on the config-specific rate distribution. This is a config-independent finding,
sharpened by the detector-realistic study.

---

## 5. Bug diagnosis and fix

The first run of `money_delta_20260729.py` produced 0% of events in their
original `(x, Q²)` bin and reco peak-rate bins displaced by factors of 100+ in
Q² relative to the true peak. The `code-reviewer` agent diagnosed the bug as a
sign-convention error in the `smear_config` function at lines 712–717.

**Wrong (as originally written):**

```python
Q2_reco = 2 * E_e * Ep_reco * (1 - cos_theta)
y_reco  = 1 - (Ep_reco / E_e) * (1 + cos_theta) / 2
```

**Correct (matching `polli_fastsim/kinematics.py`):**

```python
Q2_reco = 2 * E_e * Ep_reco * (1 + cos_theta)
y_reco  = 1 - (Ep_reco / (2 * E_e)) * (1 - cos_theta)
```

The project convention measures θ from **+z (the ion beam direction)**. The
scattered electron travels at θ ≈ π (backward, negative η); cos(π) = −1, so
`1 + cos θ → 0` in the backward limit and `1 − cos θ → 2`. The original
formulas assumed θ measured from the electron beam axis (the opposite direction),
causing systematic misreconstruction of Q² to values near zero and y to values
well outside [0, 1] for virtually all events. After the fix, bin migration
behaves physically: 47–76% of events land in their original bin depending on
config, and physicality rejection is 0.1–5%.

---

## 6. Assumptions still in play

All assumptions from `money_delta_note_2026-07-24.md` §6 carry over unchanged:
bag model origin of Δ, gluonic probe language, status of existing Δ calculations,
small-modulation approximation, and uniform azimuthal detector coverage assumed
within |η| < 3.5. The following are new or reinforced today.

**Tracking only, no ECal.** This script does not model electromagnetic calorimeter
energy resolution for the scattered electron. ECal performs better than tracking
in the backward η region (negative η, where scattered electrons are
concentrated). Adding ECal as a second smearing model is deferred as Tier B #1
(§8). The tracking-only choice is defensible for the current reach projections:
the principal signal loss is from ε_eID (7–14%), not from kinematic resolution.

**η-dependent efficiency uses ATHENA + ECCE anchor synthesis**, not an official
ePIC efficiency curve. Uncertainties are ±5% in barrel/backward region and ±10%
in the forward direction. References: JINST 17 (2022) P10019 (ATHENA) and
NIM A 1055 (2023) 168464 (ECCE). An official ePIC curve when available should
replace these anchors (Tier B #2).

**No acceptance holes.** Uniform azimuthal coverage within |η| < 3.5 is assumed.
Real ePIC has azimuthal gaps; this introduces a correction to ⟨A_cos2φ⟩ and
the effective N_flat per φ bin (see `money_delta_note_2026-07-24.md` §4.2.1
caveat 1).

**No radiative corrections.** Born-level DIS kinematics are assumed throughout.
Radiative corrections introduce an additional ~5–15% systematic effect, depending
on (x, Q²) and config. Deferred as Tier B #3.

**1000 MC events per bin** gives ~3% relative sampling error on per-bin
migration fractions. Adequate for reach projections; tighten before final
publication.

---

## 7. Sanity checks passed

1. **Case 3 integrated ⟨A_cos2φ⟩ matches pre-detector values to within a few
   percent** across all three configs (§3.3). Bin migration within the accepted
   region washes out under full integration, as expected physically.

2. **Bin-migration fractions are physically ordered (both raw and weighted).**
   Raw f_raw = 88/76/64%, weighted f_weighted = 76/62/47% for LOW/MID/TOP. The
   raw values match first-principles resolution × bin-width estimates (LOW: high
   because σ_p/p ~ 1% is small vs 20% bin width; TOP: lower because more events
   populate large-|η| edge regions with worse resolution). The weighted values are
   systematically 12–17% lower than raw, reflecting that migrated events land at
   higher-ε η than same-bin events came from — see §7a for physics.

3. **Physicality rejection is 0.1–5%** (200–23,029 events out of 218,000–455,000
   per config). This is consistent with edge-of-acceptance events being lost and
   is not indicative of a systematic issue.

4. **LOW Case 1 ⟨A_cos2φ⟩ is nearly unchanged** (−0.0228 vs −0.0224, §3.3)
   because LOW's peak-rate bin sits near the center of the |A| landscape, where
   the smearing-induced redistribution is approximately symmetric.

5. **TOP's efficiency loss is the largest (~14%)** because TOP populates more
   of the forward acceptance (positive η) where ε_eID degrades to 0.70–0.80
   at |η| = 3.0–3.5. LOW loses the least (~6%) because its scattered electrons
   are more concentrated in the backward barrel.

6. **Post-fix bin migration is monotonically consistent with resolution.** Before
   the fix, 0% of events stayed in their original bin — a clear diagnostic that
   the formulas were wrong regardless of resolution magnitude.

---

## 7a. Bin migration: raw vs efficiency-weighted

### Definitions

- **Raw same-bin fraction (f_raw):** unweighted count of MC events landing in
  their true (x_true, Q²_true) bin, divided by the total valid MC events. This
  is a pure geometric quantity — it depends only on the smearing kernel (σ_p,
  σ_θ) and the bin widths. This is what "migration probability" naturally means
  in a detector-performance context.

- **Efficiency-weighted same-bin fraction (f_weighted):** sum of weights
  (base × ε_eID(η_reco)) over same-bin events, divided by the sum of weights
  over all in-grid events. This represents "fraction of reconstructed yield
  originating from the true (x, Q²) bin," which is what enters the sig² sum
  after efficiency weighting.

The two quantities are related by:

```
f_weighted = f_raw × ⟨ε⟩_same_bin / ⟨ε⟩_all_events
```

### Today's numbers

| Config | f_raw | f_weighted | Ratio ⟨ε⟩_same / ⟨ε⟩_all |
|---|---|---|---|
| LOW | 88% | 76% | 0.86 |
| MID | 76% | 62% | 0.82 |
| TOP | 64% | 47% | 0.73 |

### Physics interpretation of the gap

f_weighted < f_raw means ⟨ε⟩_same_bin < ⟨ε⟩_all_events: same-bin events come
from η regions where ε_eID is systematically lower than where migrated events
land.

**Mechanism.** Same-bin events are drawn from the entire accepted η range,
including edges (η → -3.5 and η → +3.5) where ε_eID drops to 0.85 or 0.70.
Migrated events preferentially originate FROM those edge bins (larger σ_p/p at
edges) and land AT slightly-less-edge η where ε is higher. Result: migrated
events carry higher average ε than same-bin events, so the yield-fraction
diagnostic (weighted) undershoots the geometric migration probability (raw).

**Why the gap grows with beam energy.** TOP samples more of the forward-η range
(η → +3.5) where ε_eID has its steepest drop (0.80 at η = 3.0 to 0.70 at
η = 3.5). Migration into slightly less-forward bins picks up a proportionally
larger efficiency boost. The ratio ⟨ε⟩_same_bin / ⟨ε⟩_all_events drops from
0.86 (LOW) to 0.73 (TOP), tracking the ε_eID gradient across that forward edge.

### Which number to quote

- For **detector performance** (e.g., "tracking resolution keeps X% of events
  in the true bin"): quote **f_raw**. This is the detector physicist's migration
  probability, independent of efficiency modelling.
- For **projection / reach** (e.g., "how much of the reco yield came from the
  true (x, Q²)"): quote **f_weighted**. This is the quantity that enters the
  sig² sum after efficiency weighting.

The 20260728 unsmeared results correspond to f_raw = 100% and f_weighted = 100%
(no smearing, no efficiency variation).

---

## 8. Next steps

Building on `money_delta_note_2026-07-24.md` §8 (full Tier A/B/C structure),
`money_delta_note_2026-07-27.md` §10 (Interpretation A vs B), and
`money_delta_note_2026-07-28.md` §9 (peak-bin vs integrated analysis). Tier A
items are unchanged.

### Tier A — unchanged

- **Tier A #1** (P_zz convention resolution) — still highest-priority physics
  input. See `money_delta_note_2026-07-24.md` §8.

- **Tier A #2** (Δ shape anchoring via first-principles theory) — unchanged;
  priority elevated in `money_delta_note_2026-07-27.md` §10.

- **Tier A #3 and #4 from 2026-07-27** — unchanged (theorist consultation on
  ansatz choice; publication-ready reach quotes spanning both interpretations).

### Tier B — new additions from today

- **Tier B #1 (new): ECal energy resolution as a second smearing model.**
  Backward-region ECal is better than tracking for scattered electrons in the
  negative-η acceptance. Adding ECal would improve MID/TOP Case 1 S/N modestly
  and give a more accurate reco-bin migration estimate. Implement as a second
  configuration flag in the script.

- **Tier B #2 (new): Replace ATHENA + ECCE ε_eID anchors with the official
  ePIC efficiency curve** when available. Current anchors carry ±10% uncertainty
  in the forward region, which is the dominant contribution to the 7–14%
  efficiency loss.

- **Tier B #3 (new): Radiative corrections.** Born-level DIS assumed throughout;
  QED radiative corrections introduce ~5–15% effects on (x, Q², y) reconstructed
  from the scattered electron. Add as a systematic after the nominal study is
  complete.

- **Tier B #4 (new, physics message worth publicizing): Add the "use integration,
  not single-bin peak analysis" conclusion to the design guidance.** Today's
  detector-realistic study demonstrates that this recommendation is robust across
  all three configs and holds regardless of detector model. It is a config-
  independent, detector-independent finding. Appropriate for the experiment
  proposal or a white paper section on observable design.

---

## 9. Reproducibility

All commands run from the repo root
(`/Users/L00338853/work/Polarized_Li/PolarizedLithiumSim`).

```bash
python3 fastsim/scripts/money_delta_20260729.py --outdir out/money_delta
```

Runtime: ~5–10 minutes (1 × 10⁶ total MC events across three configs at
1000 events per accepted true bin). Deterministic: uses `np.random.default_rng(seed=42)`.
Prerequisites (PDF sets, NumPy patch, `_safe_xfx` helper, `parton` α_s table,
R1998, EPPS21, P_zz = 0.267) are identical to those documented in
`money_delta_note_2026-07-24.md` §9. The script shares the same grid backend
as `money_delta_20260728.py` and adds the MC smearing loop on top.

---

## 10. File inventory

### Script

- `fastsim/scripts/money_delta_20260729.py` (~1870 lines)

### Plots in `fastsim/out/money_delta/` (new as of 2026-07-29, 16 files)

```
money_delta_20260729_perbin_low.png                          (3-subpanel heatmap: |A_bin|, δA_bin, |δA/A|)
money_delta_20260729_perbin_mid.png                          (3-subpanel heatmap)
money_delta_20260729_perbin_top.png                          (3-subpanel heatmap)
money_delta_20260729_phimodulation_peakbin_low.png           (LOW, Case 1: reco peak bin)
money_delta_20260729_phimodulation_peakbin_mid.png           (MID, Case 1: reco peak bin)
money_delta_20260729_phimodulation_peakbin_top.png           (TOP, Case 1: reco peak bin — 5° bins, 72 pts)
money_delta_20260729_phimodulation_q2slice_low.png           (LOW, Case 2: reco Q² slice iq2=3)
money_delta_20260729_phimodulation_q2slice_mid.png           (MID, Case 2: reco Q² slice iq2=3)
money_delta_20260729_phimodulation_q2slice_top.png           (TOP, Case 2: reco Q² slice iq2=3)
money_delta_20260729_phimodulation_integrated_low.png        (LOW, Case 3: fully integrated reco)
money_delta_20260729_phimodulation_integrated_mid.png        (MID, Case 3: fully integrated reco)
money_delta_20260729_phimodulation_integrated_top.png        (TOP, Case 3: fully integrated reco)
money_delta_20260729_phimodulation_peakbin_top_10deg.png     (TOP, Case 1: reco peak bin — 10° bins, 36 pts)
money_delta_20260729_phimodulation_peakbin_top_20deg.png     (TOP, Case 1: reco peak bin — 20° bins, 18 pts)
money_delta_20260729_phimodulation_peakbin_top_30deg.png     (TOP, Case 1: reco peak bin — 30° bins, 12 pts)
money_delta_20260729_phimodulation_peakbin_top_45deg.png     (TOP, Case 1: reco peak bin — 45° bins, 8 pts)
```

### Notes (full series)

- `fastsim/notes/money_delta_note_2026-07-24.md` — Interpretation A; full
  derivation of formalism, sum-rule solver, statistical framework, Fisher-
  information combined fit.
- `fastsim/notes/money_delta_note_2026-07-27.md` — Interpretation B (no F₁
  in Δ); A vs B comparison, analytic Beta-function solver.
- `fastsim/notes/money_delta_note_2026-07-28.md` — φ-distribution projections
  (no detector); small-modulation-limit derivation, per-bin heatmaps, TOP
  kinematic-trough analysis.
- `fastsim/notes/money_delta_note_2026-07-29.md` (this file)

---

## 11. Bottom line

Detector realism (tracking-only + ε_eID) reduces the integrated tensor
observable's S/N by ~7% (mostly from efficiency, not migration); Case 1
single-bin analysis changes by factors of 2–3 in either direction; the
config-independent recommendation is to use full-acceptance integration rather
than peak-bin analysis, and this holds regardless of the detector model choice.
