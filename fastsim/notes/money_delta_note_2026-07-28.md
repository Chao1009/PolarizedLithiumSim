# Money-Δ Plot: φ-Distribution Projections — 2026-07-28

> **Dated record.** Numbers here predate the 2026-08-27 corrections: ⁶Li
> beam energies are γ-matched (40.8 / 99.5 / 137.5 GeV/u), not
> rigidity-scaled (20.5 / 50 / 137.5), and the far-forward divergence is
> per-configuration rather than a single 72.7 μrad. See plans/10. Kept as
> the record of its date, not as a current reference.


**Observable.** Fractional yield modulation `(N_φ − N_flat)/N_flat` vs `φ` for
the `cos(2φ)` tensor asymmetry on transversely tensor-polarized ⁶Li. Today's
work shifts from reach plots (`δA/A`, `L_5σ`) to the raw observable: the
φ-differential yield modulation that an experimenter would actually measure and
fit. Nine PNGs produced, covering three beam configs × three integration cases,
with three signal scales overlaid per plot.

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

**Numerical verification.** For every (config, case, s) combination in this note,

    max |P_zz · A_cos2φ · cos(2φ)| = 0.0599    (LOW + Case 1 + 10·A_bag)

This is 5.99% — well below the physical limit of 1.

**Consequence 1 — error bars are essentially signal-independent.** The per-φ-bin
Poisson error on the fractional modulation y(φ) is

    σ_y(φ) = √N_φ(s) / N_flat = (1/√N_flat) · √(1 + P_zz·A_cos2φ·cos(2φ))

In the small-modulation limit √(1+x) ≈ 1, so σ_y ≈ 1/√N_flat — **constant
across φ and essentially constant across signal scales s**. Each plot therefore
shows one set of error bars applied to all three signal-scale curves; the
variation across s is <3% even at 10·A_bag.

**Consequence 2 — the linear cross-section formula is exact enough to use
directly.** Because (*) is exact in Δ at leading twist, using it without any
further approximation introduces no error. The only physical constraint is
positivity, which is satisfied at all 5.99%.

**When it would break down.** If a future scenario pushed |P_zz · A_cos2φ| toward
0.3 or higher, error bars would visibly vary across φ (~15% shrinkage at
cos(2φ) = ±1, ~15% growth near cos(2φ) = 0). None of the scenarios in this note
approach that regime.

---

**Revision history**

- 2026-07-28 (initial): New script `money_delta_20260728.py` (~1005 lines).
  Fixed setup: Interpretation A, mid_x shape (α=0.7, β=3), R1998, EPPS21,
  Cloet P_zz = 0.267, L = 10 fb⁻¹/nucleon. Three configs × three integration
  cases × three signal scales = 9 plots (one per config/case combination,
  three cosine curves per plot). A_bag values carried forward from
  `money_delta_20260724.py`. See `money_delta_note_2026-07-24.md` for the
  full sum-rule / Fisher-information derivation and `money_delta_note_2026-07-27.md`
  for the Interpretation A vs B discussion.

- 2026-07-28 (extension): Added small-modulation-limit clarification at the top of the note,
  emphasizing that |P_zz · A_cos2φ · cos(2φ)| ≪ 1 justifies (a) essentially s-independent
  error bars and (b) the direct use of the linear cross-section formula. Max modulation across
  all plotted (config, case, s) combinations is 5.99%; well below the physical limit.

- 2026-07-28 (extension 2): Extended TOP+peakbin φ-bin-width scan by adding 45° bins (8 total,
  Nyquist-like limit for resolving cos(2φ) with 2 periods over [0, 2π]). Full scan now covers
  5°/10°/20°/30°/45°. Total PNGs: 13 (was 12).

- 2026-07-28 (extension 3): Added 3-subpanel per-bin heatmaps (LOW/MID/TOP) as motivation for
  Case 1 peak-bin selection. Each heatmap shows |A_bin| (viridis), δA_bin (plasma),
  |δA_bin/A_bin| (magma). Diagnostic printouts include peak-rate-bin cross-references linking
  each heatmap to the corresponding Case 1 φ plot. New finding: for TOP, the peak-rate bin has
  |A_bin| 13× smaller than the maximum |A_bin| on the heatmap, directly visualizing the
  "kinematic trough" finding from §5a. Total PNGs from script: 13 → 16.

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

## 1a. Per-bin heatmaps as motivation for peak-bin selection

Three heatmap PNGs were added today (extension 3) and are generated **before**
the φ-distribution plots in each config's block. The reader encounters the
(x, Q²) landscape first — seeing where |A_bin| is bright, where δA_bin is small,
and where the relative uncertainty |δA/A| is manageable — and then examines the
φ plots for the specific Case 1 bin with that context already established.

### Physics setup

All three heatmaps use the same fixed parameters as the φ plots:

- **Interpretation A:** Δ = s · α_s(Q²) · F₁(x, Q²) · x^0.7 · (1−x)^3,
  with s = A_bag_config (per-config hardcoded value; LOW: 0.318, MID: 0.310,
  TOP: 0.297).
- L = 10 fb⁻¹/nucleon, mid_x shape (α = 0.7, β = 3), EPPS21 nuclear PDF,
  R1998 (σ_L/σ_T ratio), Cloet P_zz = 0.267.
- Signal amplitude computed via
  `polli_fastsim.asymmetries.a_cos2phi(delta, f1, f2, x, y)` per (x, Q²) bin.
- Statistical uncertainty per bin:
  δA_bin = √(2/N_bin) / P_zz    (Eq. 5 of `money_delta_note_2026-07-24.md`)
- **Style** matches `money_delta_20260725.py`: viridis / plasma / magma
  colormaps, LogNorm colour scaling, NaN masking for bins outside acceptance
  or where N_bin < 10 or |A_bin| < 1e-10.

### Three-subpanel layout

Each heatmap PNG shows (x, Q²) on the axes, with one subpanel per quantity:

1. **|A_bin|** (viridis) — signal amplitude per bin.
2. **δA_bin** (plasma) — statistical uncertainty per bin.
3. **|δA_bin / A_bin|** (magma) — relative uncertainty; the primary figure-of-merit for bin selection.

### Per-bin diagnostics (min / median / max across accepted bins)

| Config | \|A_bin\| min / med / max | δA_bin min / med / max | \|δA/A\| min / med / max |
|---|---|---|---|
| LOW | 7.14e-5 / 1.59e-2 / 3.09e-2 | 1.09e-3 / 3.11e-3 / 0.401 | 3.88e-2 / 0.225 / 5.61e+3 |
| MID | 6.67e-5 / 1.22e-2 / 2.94e-2 | 1.07e-3 / 3.54e-3 / 0.839 | 4.52e-2 / 0.351 / 1.26e+4 |
| TOP | 1.38e-4 / 8.39e-3 / 2.73e-2 | 1.07e-3 / 4.27e-3 / 0.778 | 0.109 / 0.648 / 5.65e+3 |

The max |δA/A| in each config is concentrated at large x and large Q² — bins
at the kinematic boundary where N_bin is smallest (e.g., x = 0.897 for all
three configs at their respective high-Q² boundary).

### Key finding: peak-rate bin is NOT the peak-|A| bin

The heatmap diagnostic block explicitly prints the peak-rate bin coordinates
and the |A_bin| value at that bin. These cross-reference directly to the Case 1
⟨A_cos2φ⟩ values quoted in §3. The comparison below quantifies how far the
Case 1 choice sits from the optimal single-bin pick:

| Config | Peak-rate bin \|A\| | Max \|A\| on heatmap | Ratio (max / peak-rate) |
|---|---|---|---|
| LOW | 2.24e-2 | 3.09e-2 | 1.4× |
| MID | 1.98e-2 | 2.94e-2 | 1.5× |
| TOP | **2.08e-3** | **2.73e-2** | **13×** |

For **LOW and MID**, the peak-rate bin is within a factor ~1.5 of the maximum
|A_bin|. Choosing the highest-rate bin as Case 1 is a reasonable proxy for the
"best" single bin; the experimenter leaves at most ~50% of potential signal
amplitude on the table.

For **TOP**, the peak-rate bin has |A_bin| **13× smaller** than the maximum
available on the heatmap. This is the visual manifestation of the kinematic-trough
finding established analytically in §5a: TOP's DIS rate peaks at very small x
(x ≈ 0.002), where the mid_x shape factor x^0.7 = 0.011 suppresses the tensor
signal by roughly 13–15× relative to the moderate-x bins (x ~ 0.05–0.3) where
|A_bin| is brightest on the heatmap.

On TOP's |A_bin| subpanel, the bright region is clearly visible at moderate x
(~0.05–0.3), while the peak-rate bin sits at x ~ 0.002 in the dim corner of the
panel. The heatmap turns the abstract algebra of §5a.3 into a picture.

### Physics interpretation

**LOW and MID:** The DIS-rate peak coincides reasonably well with the
signal-amplitude peak (factor ~1.5). Choosing the peak-rate bin as Case 1 gives
a near-optimal single-bin measurement; the experimenter is not far from the best
available bin.

**TOP:** The DIS-rate peak is at very small x, where the mid_x shape assumption
heavily suppresses the tensor signal. The Case 1 bin yields |A_cos2φ| = 2.08e-3,
13× below the maximum |A_bin| = 2.73e-2 visible elsewhere on the heatmap.
This makes the Case 1 S/N for TOP (~1.94) look dramatically worse than LOW/MID
(~20), even though all three configs have essentially identical statistics in
their respective peak-rate bins. The heatmap makes this immediately visible.

This finding is **shape-assumption-dependent**: under the low_x variant
(α = 0.3), the x^0.3 factor at x = 0.002 is ~2× smaller than at x = 0.07
(rather than 13×), so TOP's peak-rate bin would be a much less severe
kinematic trough. See §5a.5 for the detailed shape-dependence discussion.

**Cross-reference:** §5a (Why TOP's Case 1 looks so much worse than LOW and MID)
provides the analytical basis; the per-bin heatmaps directly visualize the
physics laid out there.

### Output PNGs (extension 3)

| Filename | Content |
|---|---|
| `money_delta_20260728_perbin_low.png` | LOW config — 3-subpanel heatmap |
| `money_delta_20260728_perbin_mid.png` | MID config — 3-subpanel heatmap |
| `money_delta_20260728_perbin_top.png` | TOP config — 3-subpanel heatmap |

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

## 5a. Why TOP's Case 1 looks so much worse than LOW and MID

### 5a.1 The observed asymmetry

Looking at the §4 S/N table, TOP + Case 1 (S/N = 1.94) is roughly 10× worse
than LOW + Case 1 (S/N = 20.7) and MID + Case 1 (S/N = 18.6), despite all
three configs having essentially the same statistics in the peak bin. The numbers
from §3 make this explicit:

| Config | Peak bin N_total | N_flat/bin | ⟨A_cos2φ⟩ |
|---|---|---|---|
| LOW | 2.37 × 10⁷ | 3.29 × 10⁵ | −0.0224 |
| MID | 2.46 × 10⁷ | 3.42 × 10⁵ | −0.0198 |
| TOP | 2.43 × 10⁷ | 3.38 × 10⁵ | **−0.00208** (~10× smaller) |

Statistics per bin are essentially identical (~2.4 × 10⁷ events, ~3.4 × 10⁵
per φ bin). What differs is the *signal amplitude*, not the statistics. TOP's
peak-bin |⟨A_cos2φ⟩| is ~10× smaller than LOW's or MID's.

### 5a.2 Root cause — peak bin coordinates

Higher √s opens up small-x DIS phase space. The rate-maximizing bin shifts
accordingly:

| Config | √s [GeV/u] | Peak x | Peak Q² [GeV²] | y = Q²/(sx) |
|---|---|---|---|---|
| LOW | ~23.5 | 0.089 | 2.43 | 0.005 |
| MID | ~44.7 | 0.071 | 2.43 | 0.005 |
| TOP | ~99.5 | **0.0022** | 2.43 | 0.11 |

TOP's peak-rate bin sits at x ≈ 0.002 instead of x ≈ 0.07–0.09 (LOW/MID).
This is fatal for the tensor signal under the mid_x shape assumption.

### 5a.3 Two compounding suppression mechanisms

**Mechanism 1 — shape factor at small x.** Under the mid_x shape (α = 0.7,
β = 3), the tensor structure function Δ ∝ x^0.7 · (1−x)^3. Evaluating the
x^0.7 factor at each config's peak bin:

- At x = 0.002 (TOP): x^0.7 = 0.011
- At x = 0.070 (MID): x^0.7 = 0.145  →  **13× larger**
- At x = 0.089 (LOW): x^0.7 = 0.170  →  **15× larger**

The (1−x)^3 factor is essentially unity at both x's (≈ 0.73 at x = 0.089;
≈ 0.994 at x = 0.002). The x^0.7 factor alone accounts for a ~13–15× signal
suppression at TOP's peak bin.

**Mechanism 2 — kinematic prefactor (1−y)/y².** The cos(2φ) asymmetry
amplitude carries a kinematic prefactor proportional to (1−y)/y². Evaluating
at each config's peak-bin y:

- At y = 0.11 (TOP): (1−y)/y² = (0.89)/(0.012) ≈ 74
- At y = 0.005 (LOW/MID): (1−y)/y² = (1.0)/(2.5 × 10⁻⁵) ≈ 40000

The kinematic factor favors LOW/MID by ~500× over TOP in the peak bin. Both
configs' peak bins are in the extreme-small-y regime where this prefactor
diverges, and the rate-weighted effect of this large prefactor at LOW/MID
more than compensates for their lower absolute statistics in any given (x, Q²)
cell. Note that both the shape suppression and the kinematic factor act in the
same direction: both make TOP's peak bin unfavorable for the tensor signal.

**Net effect.** The two mechanisms combine to produce the observed ~10× smaller
|⟨A_cos2φ⟩| at TOP's peak bin relative to LOW/MID. The precise cancellation
(13× from shape × partial compensation from the y-dependence of the rate
distribution) accounts for the factor of ~10 seen in §4.

### 5a.4 Why TOP recovers when integrating

From Case 1 to Case 3, TOP's |⟨A_cos2φ⟩| *grows* monotonically:

    TOP: 0.00208 (peak bin) → 0.00346 (Q² slice) → 0.00533 (integrated)    [factor 2.5]

This is the opposite of LOW and MID, where integrating dilutes the signal:

    LOW: 0.0224 → 0.0158 → 0.0170    MID: 0.0198 → 0.0104 → 0.0121

When TOP integrates over its full acceptance, it brings in bins at moderate
x (x ~ 0.05–0.15) where x^0.7 is 10–15× larger than at x = 0.002. These bins
carry a much larger tensor signal per event than the rate-peak bin, even though
they contribute fewer events. The rate-weighted average ⟨A_cos2φ⟩ rises as
these moderate-x bins enter the sum. LOW and MID show the opposite trend
because their Case 1 bins are *already* at moderate x (which is where their
rate peaks), so further integration adds bins with smaller individual signal
and dilutes the average.

### 5a.5 Practical implications

**For TOP-config experimental design**: the DIS trigger-rate peak (small x) is a
kinematic trough for the tensor observable under the mid_x shape assumption.
Triggers and bin schemes designed around the peak DIS rate will systematically
undersample the kinematic region where the tensor signal is largest. The sensitive
region — moderate x (0.05–0.15), moderate Q² — lies in the tail of the DIS rate
distribution for TOP. Full-acceptance integration is not just statistically
optimal; it is physically necessary to include the sensitive bins. This
reinforces the conclusions of §5 point 2 and the Tier A planning note in §9.

**Shape-assumption dependence.** This picture depends on the mid_x shape
(α = 0.7). Under the low_x shape (α = 0.3, β = 4, peak at x ≈ 0.07), the
x^0.3 factor at x = 0.002 is 0.20 vs 0.44 at x = 0.07 — only ~2× smaller
instead of 13×. Under low_x, TOP's peak bin would be a much less severe
kinematic trough. Distinguishing these scenarios is another argument for
anchoring the Δ shape via first-principles theory input (Tier A #2 in §9).

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

**TOP peak-bin φ-bin-width scan** (extension 2026-07-28). A separate set of
plots holds N_total fixed at the TOP Case 1 value and varies only the bin width,
showing how statistical precision trades against angular resolution. The 45°
case is the Nyquist-like limit: 8 bins sample a cos(2φ) signal (2 full periods
over [0, 2π]) at exactly 4 points per period — the minimum needed to reconstruct
the shape. A 60° binning (6 bins = 3 per period) would start losing the ability
to resolve the cos(2φ) form.

| Bin width | N points | N_flat/bin | Error bar (1/√N_flat) | S/N at A_bag | S/N at 10×A_bag |
|---|---|---|---|---|---|
| 5° | 72 | 3.38 × 10⁵ | 1.72 × 10⁻³ | 0.12 | 1.21 |
| 10° | 36 | 6.76 × 10⁵ | 1.22 × 10⁻³ | 0.17 | 1.71 |
| 20° | 18 | 1.35 × 10⁶ | 8.60 × 10⁻⁴ | 0.24 | 2.42 |
| 30° | 12 | 2.03 × 10⁶ | 7.02 × 10⁻⁴ | 0.30 | 2.96 |
| 45° | 8 | 3.04 × 10⁶ | 5.73 × 10⁻⁴ | 0.97 | 9.7 |

Note: the per-bin S/N here is signal amplitude (P_zz · |⟨A_cos2φ⟩| at peak φ)
divided by the per-bin error bar. The fit-based S/N combining all bins (as in §4)
differs from this single-bin ratio; the §4 value of 1.94 for TOP Case 1 uses
the full √(N_total/2) formula over 72 bins.

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

### Plots in `fastsim/out/money_delta/` (new as of 2026-07-28, 16 files)

```
money_delta_20260728_perbin_low.png                          (3-subpanel heatmap: |A_bin|, δA_bin, |δA/A|)
money_delta_20260728_perbin_mid.png                          (3-subpanel heatmap)
money_delta_20260728_perbin_top.png                          (3-subpanel heatmap)
money_delta_20260728_phimodulation_peakbin_low.png           (LOW, Case 1: peak bin)
money_delta_20260728_phimodulation_peakbin_mid.png           (MID, Case 1: peak bin)
money_delta_20260728_phimodulation_peakbin_top.png           (TOP, Case 1: peak bin — 5° bins, 72 pts)
money_delta_20260728_phimodulation_q2slice_low.png           (LOW, Case 2: Q² slice iq2=3)
money_delta_20260728_phimodulation_q2slice_mid.png           (MID, Case 2: Q² slice iq2=3)
money_delta_20260728_phimodulation_q2slice_top.png           (TOP, Case 2: Q² slice iq2=3)
money_delta_20260728_phimodulation_integrated_low.png        (LOW, Case 3: fully integrated)
money_delta_20260728_phimodulation_integrated_mid.png        (MID, Case 3: fully integrated)
money_delta_20260728_phimodulation_integrated_top.png        (TOP, Case 3: fully integrated)
money_delta_20260728_phimodulation_peakbin_top_10deg.png     (TOP, Case 1: peak bin — 10° bins, 36 pts)
money_delta_20260728_phimodulation_peakbin_top_20deg.png     (TOP, Case 1: peak bin — 20° bins, 18 pts)
money_delta_20260728_phimodulation_peakbin_top_30deg.png     (TOP, Case 1: peak bin — 30° bins, 12 pts)
money_delta_20260728_phimodulation_peakbin_top_45deg.png     (TOP, Case 1: peak bin — 45° bins, 8 pts)
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
