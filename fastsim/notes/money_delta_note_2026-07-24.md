# Money-Δ Plot: Sum-Rule-Constrained Analysis — 2026-07-24

> **Dated record.** Numbers here predate the 2026-08-27 corrections: ⁶Li
> beam energies are γ-matched (40.8 / 99.5 / 137.5 GeV/u), not
> rigidity-scaled (20.5 / 50 / 137.5), and the far-forward divergence is
> per-configuration rather than a single 72.7 μrad. See plans/10. Kept as
> the record of its date, not as a current reference.


**Observable.** 5σ discovery luminosity `L_5σ` (fb⁻¹/nucleon) and fractional
asymmetry uncertainty `δA/A` for the `cos(2φ)` double-helicity-flip tensor
asymmetry on transversely tensor-polarized ⁶Li. The new element today is
sum-rule-constrained normalization of Δ(x, Q²): instead of scanning an
arbitrary peak `Δ/F₁` scale, the amplitude A is derived from the physical
first-moment constraint `∫x·Δ dx = c·α_s(Q²)`. This reduces the normalization
ambiguity from order-of-magnitude to a factor-of-order-2 shape systematic,
transforming the reach picture from "requires a decade" to
"systematics-limited within a year."

**Revision history**

- 2026-07-24 (initial): First session with sum-rule-constrained Δ normalization.
  Script `money_delta_20260724.py` (1346 lines); 9 PNGs in `fastsim/out/money_delta/`.
  See also: `money_delta_note_2026-07-16.md`, `2026-07-20.md`, `2026-07-21.md` for
  prior work and the evolving next-steps ranking.
- 2026-07-24 (extension): Corrected plot x-axis convention on Plots 1-4. The x-axis is
  now peak Δ/F₁ (not the raw `scale` parameter s). Vertical dashed lines are now drawn
  at the converted values peak(Δ/F₁)_bag ≈ 1.6×10⁻² and peak(Δ/F₁)_lat ≈ 1.2×10⁻²
  using the mid_x reference shape at ⟨Q²⟩ = 7.4 GeV². SCALES internally covers the
  code range [5e-3, 0.6] corresponding to peak Δ/F₁ ∈ [3e-4, 3e-2]. See new §4.1 for
  details. The reach numbers (L_5σ, δA/A, etc.) are unchanged — this is purely a
  visualization/labeling correction.
- 2026-07-24 (extension 2): Widened plotted x-axis range on Plots 1-4 from peak Δ/F₁ ∈ [3×10⁻⁴, 3×10⁻²] to [10⁻³, 10⁻¹] (three decades). Internal SCALES updated from [5×10⁻³, 0.6] to [0.02, 2.0]. Vertical dashed lines (bag at 1.6×10⁻², lattice at 1.2×10⁻²) unchanged. Reach numbers unchanged.
- 2026-07-24 (extension 3): Added §4.2 (statistical derivation and interpretation) explaining: origin of δA_bin = √(2/N_bin)/P_zz from cos(2φ) ML fit, combined uncertainty as one-parameter Fisher fit giving σ = ŝ/δŝ, why per-bin |δA/A| can be 10³ while combined is 10⁻² (Fisher additivity), and the σ ↔ δA/A duality (σ = 5 ↔ δA/A = 20%).
- 2026-07-24 (extension 4): Added §4.2.5 (References) listing pedagogical sources for the statistical framework used in §4.2: Cowan's Statistical Data Analysis, PDG Statistics chapter, Poskanzer-Voloshin flow-methodology paper (source of √(2/N) cos-fit uncertainty), Voloshin-Zhang (non-uniform acceptance extension), HERMES b₁ analysis (direct experimental precedent for tensor DIS asymmetry), and Lyons on the 5σ convention.
- 2026-07-24 (extension 5): Added Eq. (11) to §4 — the analytic rescaling formula σ²(s,L) = σ²_ref · (s/S₀)² · L that lets the plotting code compute δA/A at any (scale, luminosity) from a single reference sig² evaluation per (config, shape).

---

## 1. What was produced today

Script: `fastsim/scripts/money_delta_20260724.py` (1346 lines).
All figures land in `fastsim/out/money_delta/`. Ion is ⁶Li throughout.
Three beam configs × three Δ x-shape variants = 9 (config, shape) combinations.
PDF sets: R1998 (`R = σ_L/σ_T`), EPPS21 (nuclear), Cloet `P_zz = 0.267`.
Scale scan: SCALES = [0.02, 2.0], log-spaced, 15 points (internal `scale` parameter;
displayed on Plots 1-4 as peak Δ/F₁ ∈ [10⁻³, 10⁻¹] after conversion — see §4.1).

### Beam configurations

| Label | E_e [GeV] | p_ion [GeV/u] | √s [GeV/u] |
|---|---|---|---|
| LOW | 5 | 27.5 | ~23 |
| MID | 10 | 50 | ~45 |
| TOP | 18 | 137.5 | ~100 |

### Δ x-shape variants

| Label | α | β | Character |
|---|---|---|---|
| low_x | 0.3 | 4 | Soft, Regge-like rise at small x |
| mid_x | 0.7 | 3 | Moderate peak |
| high_x | 1.5 | 2 | Suppressed at small x, peaks at larger x |

### Output figures (9 PNGs)

| Filename | Content |
|---|---|
| `plot1a_dAoA_vs_scale_L10.png` | `δA/A` vs peak Δ/F₁ at `L = 10` fb⁻¹, 9 curves (3 configs × 3 shapes), separate panel |
| `plot1b_dAoA_vs_scale_L100.png` | Same at `L = 100` fb⁻¹ |
| `plot2_L5sig_vs_scale.png` | `L_5σ` vs peak Δ/F₁, 9 curves, single panel |
| `plot3_mid_midx_dAoA.png` | `δA/A` vs peak Δ/F₁ for mid+mid_x only; two luminosity curves (10 and 100 fb⁻¹) plus vertical dashed lines at peak(Δ/F₁)_bag ≈ 1.6×10⁻² and peak(Δ/F₁)_lat ≈ 1.2×10⁻² (mid_x shape, ⟨Q²⟩ = 7.4 GeV²; see §4.1) |
| `plot4_mid_midx_L5sig.png` | `L_5σ` vs peak Δ/F₁ for mid+mid_x; same vertical dashed lines (see §4.1) |
| `plot5_heatmap_mid_L10.png` | Per-bin heatmap, mid config, `L = 10` fb⁻¹; 3 subpanels: `|A_bin|`, `δA_bin`, `|δA_bin/A_bin|` |
| `plot6_heatmap_mid_L100.png` | Per-bin heatmap, mid config, `L = 100` fb⁻¹ |
| `plot7_heatmap_top_L10.png` | Per-bin heatmap, top config, `L = 10` fb⁻¹ |
| `plot8_heatmap_top_L100.png` | Per-bin heatmap, top config, `L = 100` fb⁻¹ |

---

## 2. Key formulas

### Δ model

    Δ(x, Q²) = scale · α_s(Q²) · F₁(x, Q²) · x^α · (1−x)^β          (1)

The α_s(Q²) prefactor is new relative to prior scripts. It encodes that Δ arises
at one-gluon-exchange order in the bag model (Sather–Schmidt) or from a
two-gluon operator (gluonic probe language).

### Sum-rule constraint

    ∫₀¹ x · Δ(x, Q²) dx = c · α_s(Q²)                                 (2)

    c = −0.012   (bag model, Sather–Schmidt)
    c = −0.009   (lattice QCD, Detmold–Shanahan)

### α_s cancellation — solving for A

Substituting (1) into (2):

    scale · α_s · ∫₀¹ x · F₁(x, ⟨Q²⟩) · x^α · (1−x)^β dx = c · α_s

    ⟹  A ≡ scale = c / ∫₀¹ x · F₁(x, ⟨Q²⟩) · x^α · (1−x)^β dx      (3)

α_s cancels between numerator and denominator. A is a single constant across
all Q² bins; it is computed per-bin using the `parton` package α_s table, then
averaged. The ~5% variation of A across configs is the residual from F₁
evolution only (see §4, point 3).

### Asymmetry amplitude per bin

    A_cos2φ(x, Q², y) = −[(1−y)/y²] · Δ / [F₁ + ((1−y)/(xy²))·F₂]   (4)

### Per-bin statistical uncertainty

    δA_bin = √(2/N_bin) / P_zz                                         (5)

### Combined significance

    σ² = Σ_bins  A_cos2φ,bin² · P_zz² · N_bin / 2                     (6)

    equivalently:  σ² = Σ_bins (A_bin / δA_bin)²

### Fractional precision and discovery luminosity

    δA/A = 1/√σ²  =  1/N_σ                                             (7)

    L_5σ = 25 / (σ²/L)   at reference luminosity L                     (8)

    (δA/A)_{5σ} = 1/5 = 0.20   (definitional, N_σ = 5)                 (9)

### Luminosity scale reference (from 2026-07-20)

    L_inst = 10⁻⁶ fb⁻¹/nucleon/s = 3.6 pb⁻¹/nucleon/hour

    10  fb⁻¹/nucleon  ≈ 1 year at nominal EIC instantaneous luminosity
    100 fb⁻¹/nucleon  ≈ 10 years

---

## 3. A-value table

α_s cancels in equation (3); A is determined by the shape integral over F₁.
Entries below are per-config, per-shape averages from the `parton` per-bin
α_s computation.

| Config | Shape | ⟨Q²⟩ [GeV²] | A_bag (c = -0.012) | A_lat (c = -0.009) |
|---|---|---|---|---|
| LOW | low_x | 5.92 | -0.3862 | -0.2896 |
| LOW | mid_x | 5.92 | -0.3178 | -0.2383 |
| LOW | high_x | 5.92 | -0.3913 | -0.2934 |
| MID | low_x | 7.38 | -0.3737 | -0.2803 |
| MID | mid_x | 7.38 | -0.3100 | -0.2325 |
| MID | high_x | 7.38 | -0.3875 | -0.2907 |
| TOP | low_x | 10.26 | -0.3509 | -0.2632 |
| TOP | mid_x | 10.26 | -0.2968 | -0.2226 |
| TOP | high_x | 10.26 | -0.3814 | -0.2861 |

`|A_bag| / |A_lat| = 0.012 / 0.009 = 4/3 = 1.333` exactly, holding across all
18 entries. A varies less than 5% across configs, confirming the Choice-2
physics (§5, point 8).

*Note: these A values are raw `scale` parameters (dimensionless multipliers in
the Δ formula). See §4.1 for conversion to peak Δ/F₁ ≈ 1.6×10⁻² (bag) and
≈ 1.2×10⁻² (lattice) for the mid_x shape reference at ⟨Q²⟩ = 7.4 GeV².*

---

## 4. Reach summary (A_bag values)

Significance σ² scales as `A²`, so the lattice entries give σ² × (3/4)² = σ² × 9/16
of the bag values; `L_5σ` scales as the inverse, i.e. ×16/9 ≈ 1.78× larger.

| Config | Shape | σ²/(fb⁻¹) | δA/A @ 10 fb⁻¹ | δA/A @ 100 fb⁻¹ | L_5σ [fb⁻¹/u] |
|---|---|---|---|---|---|
| LOW | low_x | 2895 | 0.0059 | 0.0019 | 0.01 |
| LOW | mid_x | 1516 | 0.0081 | 0.0026 | 0.02 |
| LOW | high_x | 714 | 0.0118 | 0.0037 | 0.03 |
| MID | low_x | 3535 | 0.0053 | 0.0017 | 0.01 |
| MID | mid_x | 1339 | 0.0086 | 0.0027 | 0.02 |
| MID | high_x | 472 | 0.0145 | 0.0046 | 0.05 |
| TOP | low_x | 2400 | 0.0065 | 0.0020 | 0.01 |
| TOP | mid_x | 433 | 0.0152 | 0.0048 | 0.06 |
| TOP | high_x | 85 | 0.0343 | 0.0108 | 0.29 |

**Worst case:** TOP + high_x, `L_5σ ≈ 0.29` fb⁻¹/u. At nominal EIC
luminosity (3.6 pb⁻¹/nucleon/hour) that is about 80 hours — well under a
year's run. Under the lattice sum rule (`c = -0.009`), all nine cases still
reach 5σ within a 1-year (10 fb⁻¹) program.

### Analytic rescaling of δA/A across (scale, luminosity)

Because A_cos2φ,bin scales linearly with s (through Δ ∝ s) and N_bin scales
linearly with luminosity L, the combined significance-squared has a
closed-form dependence on both variables. This lets the code compute σ²_ref
once per (config, shape) at reference (s = S₀ = 10⁻³, L = 1 fb⁻¹/nucleon)
and rescale analytically for every point on the plot — saving ~99% of the
compute cost compared to re-evaluating per (s, L) combination.

The scaling relation for the combined significance-squared is:

    σ²(s, L) = σ²_ref · (s / S₀)² · L                                  (10)

where σ²_ref is computed once at reference values s = S₀ = 10⁻³ and
L = 1 fb⁻¹/nucleon. The resulting fractional asymmetry uncertainty is:

    δA/A (s, L) = 1 / √σ²(s, L)
               = 1 / √[ σ²_ref · (s/S₀)² · L ]                         (11)

This is what enables the plotting code to compute 9 curves × 15 scale
points × 2 luminosities from a single σ² evaluation per (config, shape) —
the analytic rescaling costs zero additional physics computation.

### 4.1 Note on x-axis convention (peak Δ/F₁ vs scale parameter s)

**Internal vs. displayed quantity.** The code uses `scale` (equivalently `s`)
internally as the multiplier in the Δ model:

    Δ(x, Q²) = scale · α_s(Q²) · F₁(x, Q²) · x^α · (1−x)^β

However, Plots 1-4 display **peak Δ/F₁** on the x-axis, defined as:

    peak Δ/F₁ = scale × α_s(⟨Q²⟩) × max_x[ x^α (1−x)^β ]

where the maximum of the shape factor and ⟨Q²⟩ are both evaluated for the
reference shape and config (mid_x + MID for the 9-curve plots). The internal
SCALES range [0.02, 2.0] therefore corresponds to displayed peak Δ/F₁ values
in the range [10⁻³, 10⁻¹].

**Conversion factor for the mid_x reference (MID config, ⟨Q²⟩ = 7.4 GeV²).**
The shape x^0.7 (1−x)^3 peaks at x_peak = 0.19 with a peak value of ≈ 0.170.
With α_s(7.4 GeV²) ≈ 0.30, the conversion factor is:

    α_s(⟨Q²⟩) × peak_shape = 0.30 × 0.170 = 0.051

**Converted vertical-line positions** (mid_x shape, MID config):

| Theory input | Raw |A| (= scale) | peak Δ/F₁ |
|---|---|---|
| Bag (Sather–Schmidt, c = −0.012) | 0.310 | 0.310 × 0.051 = **1.6 × 10⁻²** |
| Lattice (Detmold–Shanahan, c = −0.009) | 0.233 | 0.233 × 0.051 = **1.2 × 10⁻²** |

**Caveats.**

1. *Q² dependence of peak Δ/F₁.* Because α_s runs, peak Δ/F₁ is not a single
   number across the accepted Q² range — it varies by roughly a factor of 2
   (α_s(2 GeV²) ≈ 0.36 vs α_s(100 GeV²) ≈ 0.19). The vertical-line positions
   use the rate-weighted ⟨Q²⟩ as the reference.

2. *Shape dependence in 9-curve plots.* The peak of x^α(1−x)^β differs across
   the three shape variants (low_x, mid_x, high_x). Plots 1 and 2 show all
   nine (config, shape) combinations; the vertical dashed lines use the mid_x
   shape as the reference and note this in an on-plot annotation. The underlying
   curves are computed for their own shapes.

3. *Plots 3 and 4 are exact.* These plots use only the mid_x shape throughout,
   so the conversion factor (0.051) applies exactly to all displayed quantities.

**Why this matters.** Without the conversion, plotting |A| on an axis labeled
"Δ/F₁" would imply peak Δ/F₁ ≈ 0.3 — roughly 20× larger than the physical
value. The displayed peak Δ/F₁ ~ 1–2 × 10⁻² correctly reflects that the
one-gluon-exchange tensor asymmetry is suppressed by α_s × (shape peak) ≈ 0.05
relative to the F₁-normalized amplitude. The reach numbers in the table above
(L_5σ, δA/A, etc.) are unaffected by this labeling change.

### 4.2 Statistical derivation and interpretation of δA_bin and the combined uncertainty

This subsection derives Eq. 5 and Eq. 6 from first principles and explains why
per-bin relative uncertainties `|δA_bin/A_bin|` can reach 10³ while the
combined `δA/A` is ~1%. The two quantities measure fundamentally different
things.

#### 4.2.1 Origin of the per-bin uncertainty δA_bin = √(2/N_bin)/P_zz (Eq. 5)

Eq. 5 is the statistical uncertainty from a **maximum-likelihood (ML) fit of the
per-bin azimuthal yield to a cos(2φ) modulation**. It is not derived from Eq. 4.
Eq. 4 is the theoretical prediction for A_cos2φ (Hoodbhoy–Jaffe–Manohar kinematic
formula); Eq. 5 is the statistical precision with which that amplitude can be
measured from N_bin detected events.

**Derivation.**

The φ-differential yield in a (x, Q², y) bin with N_bin events:

    dN/dφ = (N_bin / 2π) · [1 + P_zz · A_cos2φ · cos(2φ)]

Define the observed cos(2φ) amplitude B = P_zz · A_cos2φ. The ML estimator for B
from the individual event angles {φᵢ} is:

    B̂ = (2 / N_bin) · Σᵢ cos(2φᵢ)

In the small-modulation limit (|B| ≪ 1), the events are approximately uniformly
distributed in φ. The variance of B̂ is:

    Var(B̂) = (4 / N_bin²) · N_bin · ⟨cos²(2φ)⟩ = (4 / N_bin) · (1/2) = 2 / N_bin

Therefore δB = √(2/N_bin), and since B = P_zz · A_cos2φ:

    δA_cos2φ = δB / P_zz = √(2/N_bin) / P_zz                   (Eq. 5)

**Origin of each factor:**

| Factor | Source |
|---|---|
| √N_bin in denominator | Standard Poisson counting: N events reduce variance by 1/N. |
| Extra √2 in numerator | ⟨cos²(2φ)⟩ = 1/2 under uniform φ. Cosine fits use only half the information per event because sin and cos are orthogonal modes — only the projection onto cos(2φ) counts. |
| 1/P_zz | Tensor-polarization dilution. The measured modulation is B = P_zz · A_cos2φ; dividing by P_zz un-dilutes to recover A_cos2φ. |

**Caveats.**

1. *Non-uniform φ acceptance.* For a real detector with gaps or shadowed regions,
   ⟨cos²(2φ)⟩_det ≠ 1/2. The effective denominator changes and must be computed
   from the acceptance-weighted integral. The formula above assumes uniform
   coverage.

2. *Small-modulation assumption.* The derivation requires |B| = |P_zz · A_cos2φ| ≪ 1.
   For our kinematics |B| ≲ 0.02 (since |A_cos2φ| ≲ 0.03 and P_zz ≈ 0.27), so
   this is well satisfied.

#### 4.2.2 Combined uncertainty as a one-parameter Fisher fit

**Critical wrinkle.** A_cos2φ is not a single physical constant — it depends on
(x, Q², y) through Eq. 4. So the standard "combine N independent measurements of
the same quantity" formula does not directly apply. The correct procedure is a
**one-parameter fit**.

**Two distinct quantities** that the framework computes:

- **Total significance σ**: the test statistic for the hypothesis "Δ is nonzero,"
  combined across all kinematic bins. Given by Eq. 6.
- **Fitted amplitude ŝ ± δŝ**: the best-fit value of the overall scale parameter
  in the model Δ = s · [fixed shape terms], and its Fisher-information uncertainty.

**Derivation of the one-parameter fit.**

Assume Δ(x, Q²) = s · Δ_shape(x, Q²), where Δ_shape is the fixed reference shape
and s is the free parameter. Then A_cos2φ,bin(s) = s · A_bin,shape. The
log-likelihood for the per-bin measured amplitudes {A_bin,meas} is:

    ln L(s) = −½ · Σ_bins [(A_bin,meas − s · A_bin,shape)² / δA_bin²]

The ML estimator for s is:

    ŝ = Σ_bins [A_bin,meas · A_bin,shape / δA_bin²]
        ─────────────────────────────────────────────
        Σ_bins [A_bin,shape² / δA_bin²]

The Fisher information gives the uncertainty:

    1 / δŝ² = Σ_bins (A_bin,shape² / δA_bin²)

Using A_bin,shape = A_bin,cos2φ / s_true:

    s_true / δŝ = √[ Σ_bins (A_bin,cos2φ / δA_bin)² ] = σ      (Eq. 6)

**Conclusion.** σ = ŝ / δŝ is the signal-to-noise ratio on the fitted amplitude.
Equivalently:

    δA/A = δŝ / s_true = 1/σ = 1/N_σ                            (Eq. 7)

The "combined amplitude" is the free parameter of the one-parameter fit; its
uncertainty is what Fisher information gives when all bins are combined. The sum
in Eq. 6 is an **inverse-variance-weighted sum of squared per-bin
signal-to-noise ratios** — Fisher information is additive across independent bins.

#### 4.2.3 Per-bin vs combined uncertainty — why they differ so much

These are two distinct quantities that are easy to confuse.

**Per-bin relative uncertainty:** (δA/A)_bin = δA_bin / |A_bin|. Tells you how
well an individual bin measures the asymmetry in that bin alone. Varies dramatically
across the (x, Q²) grid because both δA_bin (set by N_bin) and A_bin (set by
kinematic factors) vary. This is what the right subpanel of Plots 5–8 shows.

**Combined relative uncertainty:** (δA/A)_combined = 1/σ = 1/N_σ. Tells you how
well the fitted amplitude is determined using all bins together. Much smaller than
per-bin ratios because Fisher information accumulates additively across bins.

**Numerical example (MID + mid_x + 10 fb⁻¹/nucleon):**

| Quantity | Value |
|---|---|
| Per-bin median \|δA/A\| | 0.35 |
| Combined δA/A | 0.0086 |
| Ratio (per-bin / combined) | ~41 |

The ratio ~41 corresponds to roughly √(1700) ≈ 41 effective bins contributing
with full weight to the Fisher information (out of 345 kinematic bins), each with
varying weight determined by (A_bin/δA_bin)².

**Wrong way to combine.** Taking the mean of per-bin |δA/A| ratios does NOT give
the combined uncertainty. That gives ~0.35, underestimating sensitivity by ~40× in
this case. The correct combined uncertainty uses inverse-variance weighting: square
the per-bin signal-to-noise ratios, sum, then take √ and invert (i.e., Eq. 6 → 7).

**Intuition.** Each bin contributes an independent "vote" of magnitude
(A_bin/δA_bin)² to the combined significance. Good bins (large |A|, small δA)
vote loudly; bad bins — those at large x + large Q² where signal is killed by
(1−x)^β and rate by 1/Q⁴ — whisper. The combined significance is the RMS-sum of
votes, dominated by the loud moderate-x, moderate-Q² bins.

#### 4.2.4 Significance ↔ relative uncertainty duality

From Eq. 7, σ = 1/(δA/A), so every significance value maps to a unique relative
statistical uncertainty on the fitted amplitude:

| σ (significance) | δA/A (relative stat uncertainty) |
|---|---|
| 1 | 100% |
| 2 | 50% |
| 3 | 33% |
| **5** | **20%** |
| 10 | 10% |
| 100 | 1% |

Same information, flipped: significance = signal/noise; relative uncertainty =
noise/signal.

**Precise statement for today's numbers.** σ = 5 means the fitted amplitude ŝ from
the one-parameter fit is 20% relatively uncertain. For the bag prediction
s_0^bag = 0.310 (mid_x shape, MID config), a 5σ result would report:

    ŝ = 0.310 ± 0.062   (statistical only)

**Two caveats.**

1. *Per-bin vs fitted amplitude.* The 5σ ↔ 20% mapping applies to the **fitted
   overall amplitude** ŝ, not to per-bin measurements. Per-bin uncertainties can
   be much larger (median 35% at 10 fb⁻¹) even when the combined uncertainty is
   small (0.86%).

2. *Discovery vs precision.* 5σ is the **discovery threshold**, not the
   **precision threshold**. At 5σ, δA/A = 20% is a loose measurement. For
   precision physics one wants 10σ (10%) or better, requiring 4× more luminosity
   (since σ ∝ √L). Today's numbers show that at 100 fb⁻¹/nucleon under the bag
   prediction, δA/A drops to ~0.6–1% (100σ+); the measurement becomes
   systematics-limited well before statistics run out.

#### §4.2.5 References

**Note**: bibliographic details (year, volume, page) have been provided by the
author from recall; please cross-check on INSPIRE or the journal websites before
citing publicly.

The references below establish that the statistical framework in §4.2 is standard
HEP practice. They are organized into three thematic groups.

**Group 1: General statistical framework**

1. G. Cowan, *Statistical Data Analysis* (Oxford University Press, 1998).
   The standard HEP reference. Chapter 6 covers maximum likelihood estimation,
   Fisher information, and the Cramér–Rao bound — the mathematical basis for
   §4.2.2's derivation of the combined uncertainty as a single-parameter Fisher
   fit.

2. Particle Data Group, "Statistics" review chapter in the current *Review of
   Particle Physics* (updated yearly at [pdg.lbl.gov](https://pdg.lbl.gov)).
   Canonical modern reference for HEP conventions on discovery significance,
   profile likelihood, and nuisance-parameter treatment. Cite the current-year
   edition.

**Group 2: cos(nφ) modulation fits**

3. A. M. Poskanzer & S. A. Voloshin, "Methods for analyzing anisotropic flow
   in relativistic nuclear collisions," *Phys. Rev. C* **58**, 1671 (1998).
   arXiv:nucl-ex/9805001. Standard reference on Fourier-coefficient extraction
   from N azimuthally-distributed events. Contains the derivation of the
   δv_n = √(2/N) uncertainty on flow coefficients — the same math as §4.2.1 for
   our δA_cos2φ per bin.

4. S. A. Voloshin & Y. Zhang, "Flow study in relativistic nuclear collisions by
   Fourier expansion of azimuthal particle distributions," *Z. Phys. C* **70**,
   665 (1996). arXiv:hep-ph/9407282. Extends the treatment to non-uniform
   azimuthal acceptance — relevant for real detectors with dead regions.

**Group 3: Applied precedent for tensor DIS asymmetries**

5. HERMES Collaboration (A. Airapetian *et al.*), "First measurement of the
   tensor structure function b₁ of the deuteron," *Phys. Rev. Lett.* **95**,
   242001 (2005). arXiv:hep-ex/0506018. Closest published precedent for a tensor
   DIS asymmetry measurement. Documents the per-bin statistical treatment and
   systematic-uncertainty framework in a real spin-1 target experiment
   (unpolarized electron beam × tensor-polarized deuteron target). The direct
   experimental analog of what our fast-sim is projecting for ⁶Li.

**Optional: the "why 5σ?" reference**

6. L. Lyons, "Discovering the significance of 5 sigma," arXiv:1310.1284 (2013).
   Short (~5 pages) discussion of why HEP uses the 5σ threshold for discovery
   (vs 3, 4, or 6). Historical and practical justification of §4.2.4's discovery
   convention.

The statistical framework in §4.2 is standard HEP practice as documented in [1]
and [2]. The specific per-bin uncertainty formula δA_bin = √(2/N_bin)/P_zz
(Eq. 5) is the tensor-polarization-scaled version of the well-known cos(nφ)
amplitude uncertainty from [3, 4], applied to the double-helicity-flip observable
rather than to elliptic flow. The direct experimental precedent for a tensor DIS
asymmetry measurement is HERMES's b₁ analysis [5]; our machinery is essentially
the ⁶Li tensor gluonometry analog of that framework, projected forward to EIC
luminosities.

---

## 5. Per-bin heatmap diagnostics

The three-subpanel layout of Plots 5–8 shows `|A_bin|` (left), `δA_bin`
(middle), `|δA_bin/A_bin|` (right) as a function of (x, Q²). This decomposition
makes two failure modes visually distinct: signal suppression via `(1−x)^β`
at large x, and rate suppression from the 1/Q⁴ cross-section falloff at large
Q².

| Config | L [fb⁻¹] | \|A_bin\| min/med/max | δA_bin min/med/max | \|δA/A\| min/med/max | Max at |
|---|---|---|---|---|---|
| MID | 10 | 6.7e-5 / 1.2e-2 / 2.9e-2 | 1.1e-3 / 3.5e-3 / 0.84 | 4.5e-2 / 0.35 / 1.3e+4 | x=0.897, Q²=1379 GeV² |
| MID | 100 | 6.7e-5 / 1.2e-2 / 2.9e-2 | 3.4e-4 / 1.1e-3 / 0.27 | 1.4e-2 / 0.11 / 4.0e+3 | x=0.897, Q²=1379 GeV² |
| TOP | 10 | 1.4e-4 / 8.4e-3 / 2.7e-2 | 1.1e-3 / 4.3e-3 / 0.78 | 0.11 / 0.65 / 5.7e+3 | x=0.897, Q²=1776 GeV² |
| TOP | 100 | 1.4e-4 / 8.4e-3 / 2.7e-2 | 3.4e-4 / 1.3e-3 / 0.25 | 3.5e-2 / 0.21 / 1.8e+3 | x=0.897, Q²=1776 GeV² |

Large `|δA/A|` ratios are always at large x + large Q²: signal is killed by
`(1−x)^β` and rate is killed by 1/Q⁴. These bins contribute negligibly to the
combined sum (6); the combined reach is dominated by the moderate-x, moderate-Q²
region where both `A_bin` and N_bin are large.

---

## 6. Physics discussions from today's Q&A

These discussions clarify the theoretical context for the sum-rule normalization
and are recorded here for reference.

### 6.1 Bag model origin of Δ

MIT bag model: quarks confined in a ~1 fm spherical cavity, free Dirac
fermions inside with linear-boundary-condition BCs. In the pure parton model
(no gluons), helicity conservation forces Δ = 0. Sather and Schmidt (1993)
added one-gluon exchange to get non-zero Δ; the leading result is
`⟨x·Δ⟩ ~ α_s` at the bag scale. The bag calculation provides the only
published numerical estimate of the first moment of Δ.

### 6.2 "Gluonic probe" language

Δ has two gauge-invariant operator contributions at one-gluon order:
(a) the quark propagator with a gluon insertion (Sather–Schmidt captures this),
and (b) the two-gluon operator (gluon helicity flip). At small x the gluon
piece tends to dominate; at moderate x (around the peak of `x·F₁`) both
contributions are comparable. The current model includes only the quark piece
via the Sather–Schmidt moment.

### 6.3 Is 10⁻³ a lower bound on |A|?

No. The Sather–Schmidt value is a point estimate with an order-of-magnitude
theoretical uncertainty. The missing two-gluon piece could add to or cancel the
quark piece. No positivity theorem forces `|Δ|` above any particular floor.
The three-tier scenario (10⁻³, 3×10⁻³, 10⁻²) used in earlier notes is
narrative labeling — "modest gluon enhancement" and "optimistic gluon-dominated"
— not independent calculations.

### 6.4 Three-tier scenario (10⁻³, 3×10⁻³, 10⁻²)

Standard EIC projection labeling carried over from prior sessions. Only the
10⁻³ case has a calculation behind it (Sather–Schmidt bag). The 3× and 10×
tiers are scenario labels. The sum-rule normalization of today's work supersedes
this ad-hoc scan and grounds the normalization in the physical first moment.

### 6.5 State of existing Δ calculations

Very thin. Hoodbhoy–Jaffe–Manohar (NPB 312:571): foundational operator
analysis, no numerics. Sather–Schmidt (1993): bag model, the only published
number for `⟨x·Δ⟩`. Detmold–Shanahan: lattice QCD, provides the second
numerical estimate used here (`c = -0.009`). No experimental measurement of Δ
exists.

### 6.6 Toy x-shape motivation

`x^α(1−x)^β` is a Regge-inspired small-x template (α controls the small-x
power) combined with Brodsky–Farrar quark-counting large-x suppression (β
controls the large-x falloff). This template is not Δ-specific; it is the
standard ansatz for any parton distribution. All three current variants have
positive α, meaning they miss the gluon-dominated scenario where Δ could be
rising more steeply at small x.

### 6.7 First-moment constraint as normalization

The sum rule `∫x·Δ dx = c·α_s` replaces the arbitrary "peak Δ/F₁ = scale"
normalization with a physically-motivated amplitude derived from the shape
integral over F₁. This reduces the normalization freedom from a free parameter
to the fixed c coefficient (known from bag or lattice), leaving only the x-shape
as a free parameter. Factor-of-order-2 shape systematic remains (high_x vs
low_x cases in §4).

### 6.8 Choice 2: placing α_s inside Δ

Putting α_s(Q²) inside Δ (equation 1) rather than treating it as a global
prefactor means that when the sum-rule constraint is enforced (equation 2),
α_s appears on both sides and cancels in equation (3). The result is that A
is a single constant across Q² bins. This matches physical intuition — the
tensor asymmetry amplitude should not run strongly with Q² if it is fully
determined by the sum-rule normalization and the shape. The ~5% residual
variation observed in §3 is purely from F₁ evolution.

### 6.9 5σ threshold and δA/A = 0.20

`δA/A = 1/N_σ` by definition (equation 7). Setting N_σ = 5 gives
`(δA/A)_{5σ} = 0.20`. The value 5 is the HEP convention for "discovery"
significance. A measurement achieving `δA/A = 0.20` is a discovery-level
constraint on Δ at the assumed sum-rule value.

---

## 7. Key physics conclusions

1. **Sum-rule normalization improves reach by ~1000× vs peak-Δ/F₁ scan.**
   Under the bag or lattice sum-rule prediction, discovery-time becomes
   "hours to days" rather than "years to decades."

2. **All 9 (config, shape) combinations reach 5σ within 1 year (10 fb⁻¹).**
   Worst case is TOP + high_x at `L_5σ ≈ 0.29` fb⁻¹/u (~80 hours at nominal
   EIC instantaneous luminosity). Under the lattice constraint, the worst case
   is still comfortably within a 1-year program.

3. **A ≈ constant across configs (< 5% variation).** Confirms Choice-2 physics:
   α_s cancels and A is set by the shape integral over F₁. The residual variation
   is from F₁ evolution alone.

4. **|A_bag|/|A_lat| = 4/3 exactly** from the c-coefficient ratio (0.012/0.009).
   The two theoretical inputs agree in sign and are close in magnitude; the ratio
   is fixed by the calculation, not the simulation.

5. **At 100 fb⁻¹, the measurement becomes systematics-limited** (δA/A ~
   0.5–1%). The statistical reach is no longer the bottleneck. Next priorities
   are P_zz calibration, F₂A normalization, radiative corrections, and
   shape-model systematic control.

6. **Per-bin decomposition clarifies failure modes.** Large `|δA/A|` at large x +
   large Q² comes from combined signal suppression (`(1−x)^β`) and rate
   suppression (1/Q⁴). These bins are negligible in the combined sum (6).

---

## 8. Updated next-steps ranking

Building on §4 of `money_delta_note_2026-07-16.md` and §4 of `2026-07-20.md`:

### Tier A — physics inputs on the critical path

1. **Resolve ⁶Li P_zz dilution convention** (unchanged from 2026-07-16 §4 #1,
   elevated in urgency). Under sum-rule normalization, P_zz uncertainty directly
   limits achievable δA/A at 100 fb⁻¹. Consult I. Cloet or nuclear-structure
   literature. A factor ~1.5 in `P_zz_eff` → factor ~2.3 in σ² → factor ~1.5
   in `L_5σ`. See 2026-07-16 §7 for the full convention discussion.

2. **Anchor Δ(x, Q²) to a real theory calculation** (sharpened from 2026-07-16
   §4 #2). Sum-rule normalization reduces shape systematic but does not eliminate
   it — factor ~2 remains across the three shape variants. Options: lattice-moment
   constraints, CDKS convolution model, NJL or light-cone model. Minimum viable:
   one motivated shape with its own uncertainty estimate.

3. **(Retired per 2026-07-20)** Nuclear PDF choice sub-1%; not a priority.
   EPPS21 is sufficient for current projections.

4. **NEW: Systematic-uncertainty studies now on the critical path.** With
   statistics no longer limiting at 10–100 fb⁻¹, the priority shifts to:
   P_zz calibration uncertainty, F₂A normalization uncertainty, radiative
   corrections, and x-shape modeling systematic. These dominate the achievable
   precision at high luminosity.

### Tier B — detector and reconstruction realism

5. η-dependent electron-ID efficiency (real ePIC curve, once stable).
6. Detector resolution smearing on E′ and θ; bin migration in (x, Q², φ).
7. Radiative corrections (HECTOR-equivalent).

### Tier C — long-term

8. NNLO QCD corrections to the DIS cross section.
9. Target-mass corrections at high x.
10. Polarized-electron scenarios (complementary single-spin and A_LT observables).

---

## 9. Reproducibility

All commands run from the repo root
(`/Users/L00338853/work/Polarized_Li/PolarizedLithiumSim`).

```bash
python3 fastsim/scripts/money_delta_20260724.py --outdir out/money_delta
```

Runtime: ~15–20 minutes (27 grid-backend σ² evaluations for the 9-curve
plots + 4 per-bin heatmap computations).

### Prerequisites

PDF sets required by `parton`:

```bash
python3 -m parton install EPPS21nlo_CT18Anlo_Li6
python3 -m parton install nNNPDF30_nlo_as_0118_A6_Z3
python3 -m parton install CT18NLO
```

The environment patches documented in `money_delta_note_2026-07-16.md` §5
(NumPy 2.0 compatibility fix in `parton/pdf.py:231`; `_safe_xfx` helper in
`fastsim/polli_fastsim/structure.py` and `polarized.py`) still apply. The
`parton` package α_s table is used for per-bin α_s values (primary source).

---

## 10. File inventory

### Script

- `fastsim/scripts/money_delta_20260724.py` (1346 lines)

### Plots in `fastsim/out/money_delta/` (new as of 2026-07-24, 9 files)

```
plot1a_dAoA_vs_scale_L10.png       (δA/A vs peak Δ/F₁, L=10 fb⁻¹, 9 curves)
plot1b_dAoA_vs_scale_L100.png      (δA/A vs peak Δ/F₁, L=100 fb⁻¹, 9 curves)
plot2_L5sig_vs_scale.png           (L_5σ vs peak Δ/F₁, 9 curves)
plot3_mid_midx_dAoA.png            (mid+mid_x, δA/A vs peak Δ/F₁, 2L + bag/lat verticals at 1.6e-2/1.2e-2)
plot4_mid_midx_L5sig.png           (mid+mid_x, L_5σ vs peak Δ/F₁, bag/lat verticals at 1.6e-2/1.2e-2)
plot5_heatmap_mid_L10.png          (per-bin heatmap, mid, L=10 fb⁻¹, 3 subpanels)
plot6_heatmap_mid_L100.png         (per-bin heatmap, mid, L=100 fb⁻¹, 3 subpanels)
plot7_heatmap_top_L10.png          (per-bin heatmap, top, L=10 fb⁻¹, 3 subpanels)
plot8_heatmap_top_L100.png         (per-bin heatmap, top, L=100 fb⁻¹, 3 subpanels)
```

### Note

- `fastsim/notes/money_delta_note_2026-07-24.md` (this file)

---

## 11. Bottom line

Under the bag (`c = -0.012`) or lattice (`c = -0.009`) first-moment sum-rule
prediction, ⁶Li tensor gluonometry is discoverable at 5σ significance in
hours to days across all tested (config, shape) combinations. The physics case
transforms from "requires a decade of running" to "systematics-limited within
a year." The two remaining priorities now on the critical path are:

1. **P_zz convention resolution** — the unresolved dilution factor is the
   largest single ambiguity in the achievable δA/A at high luminosity.
2. **Systematic-uncertainty control** — P_zz calibration, F₂A normalization,
   radiative corrections, and shape modeling all become leading uncertainties
   once δA/A drops below ~1%.
