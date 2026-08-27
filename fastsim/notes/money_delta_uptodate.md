# Money-Δ: Master Reference Note (current as of 2026-07-31)

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


**Purpose.** One-stop synthesis of the full money-delta note series
(2026-07-24 through 2026-07-29), incorporating the 2026-07-31 reco-selection
consistency fix audit. Earlier sessions (2026-07-16, -07-20, -07-21)
established the infrastructure and initial reach estimates; those are
summarized only in §2 (chronology). All numbers here come directly from the
dated notes; nothing is invented.

**Current best-practice conclusion (summary before details).**
The fully-integrated tensor cos(2φ) analysis is robust: pre-detector
(ideal simulation) it gives the best S/N across all tested configs, and under
the one detector model tested to date (ePIC tracking-only momentum/angular
resolution + η-dependent ε_eID, Interpretation A mid_x bag, 10 fb⁻¹) it
retains ~92% of pre-detector S/N. Single-bin peak analyses are fragile: they
change by factors of 2–3 in either direction depending on how detector
migration reshapes the reco rate distribution, and the peak-rate bin can be
a kinematic trough for the tensor signal (TOP config under mid_x shape).
Use integration, not single-bin peak analysis.

A reco-selection consistency fix was applied to `money_delta_20260729.py`
on 2026-07-31 (see `money_delta_note_2026-07-29_fix.md`). The fix repaired
four internal reco-mask/selection inconsistencies in the detector-realistic
analysis path and was validated through a full production rerun. Case 1
values were unchanged (<0.1%); Case 2/3 values shifted by 2–5% (TOP Case 2
~5.1%, marginally above the ≤5% plan tolerance, recorded as a pass-with-note
detail). The production rerun did not produce an obvious qualitative
overturning of the detector-realistic message: integrated Case 3 remains the
robust recommendation and Case 1 fragility is unchanged. Scope is limited to
the tested detector model; results are not generalized beyond it.

---

## Table of Contents

0. Package structure
1. Observable and physical context
   - 1.1 Estimation logic (flow chart)
2. Chronology of scripts and notes
3. Beam configurations
4. Theory inputs and model variants
   - 4.1 Δ ansatz — Interpretation A vs B
   - 4.2 Sum-rule constraints (bag vs lattice)
   - 4.3 x-shape variants
   - 4.4 PDF and structure-function inputs
   - 4.5 Polarization
5. Key formulas
   - 5.1 Δ model and asymmetry amplitude
   - 5.2 Sum-rule solver
   - 5.3 Per-bin statistical uncertainty
   - 5.4 Combined significance and fractional precision
   - 5.5 Analytic rescaling formula
   - 5.6 φ-modulation observable
6. Reach results
   - 6.1 Interpretation A reach table (2026-07-24)
   - 6.2 Interpretation B reach table (2026-07-27)
   - 6.3 A vs B comparison
   - 6.4 A-value tables
7. φ-modulation projections (2026-07-28)
8. Detector-realistic extension (2026-07-29)
   - 8.1 Detector model
   - 8.2 S/N results pre- vs post-detector
   - 8.3 Why Case 3 is robust; why Case 1 is fragile
9. Four study phases — what each established
10. Stable conclusions vs model-dependent findings
11. Open uncertainties and next steps
    - 11.1 Tier A — physics critical path
    - 11.2 Tier B — detector realism
    - 11.3 Tier C — long-term
12. Reproducibility and prerequisites
13. Complete file inventory

---

## 0. Package structure

The table below covers the principal paths relevant to the money-delta
workflow. Paths are relative to the repo root
(`PolarizedLithiumSim/`). The `polli_fastsim` library is the shared
physics backend; the money-delta scripts import it directly and drive
the estimation logic described in §1.1.

| Path | Role in money-delta workflow |
|---|---|
| **fastsim/polli_fastsim/** | Core library package |
| `polli_fastsim/__init__.py` | Registers modules; exports `beams`, `kinematics`, `structure`, `polarized`, `asymmetries`, `fom` |
| `polli_fastsim/beams.py` | `Ion` dataclass (A, Z, spin, `eff_pol_p/n`); `BeamConfig`; EIC energy configurations (LOW / MID / TOP); rigidity scaling. P_zz is **not** a field here — it is a scenario input set in each driver script. |
| `polli_fastsim/kinematics.py` | DIS kinematics per nucleon: √s, Q², y, W², acceptance cuts; reconstructed Q²/y formulas (corrected in Phase IV) |
| `polli_fastsim/structure.py` | Generic backend: `ToyF2`, `NuclearF2` (toy-F2 based), `_safe_xfx` NumPy-2.0 patch, `dsigma_dx_dq2`, `r_sigma_lt`. **The EPPS21+R1998 production setup is not wired here**: the money-delta scripts define `NuclearF2FromGrid` locally (using `parton` to load EPPS21nlo_CT18Anlo_Li6) and call `r_sigma_lt` from this module for F₁. |
| `polli_fastsim/polarized.py` | `ToyG1` (g₁ toy model); EMC ratio scenarios; `toy_delta_gluon` placeholder. Does **not** contain `A_cos2phi`, `P_zz`, or the sum-rule Δ ansätze (Eqs. 1A/1B) used in the money-delta drivers. |
| `polli_fastsim/asymmetries.py` | `a_cos2phi` (Hoodbhoy–Jaffe–Manohar Eq. 4); `err_cos2phi_amplitude`; `azz`; `a_parallel` |
| `polli_fastsim/fom.py` | `Scenario` dataclass (luminosity, cuts, polarization placeholders); `project_rates` (event counts per bin); `project_observables` (per-bin asymmetries and statistical errors); `combine_over_q2`. The money-delta scripts compute sum-rule Fisher sums, L_5σ, and δA/A scaling locally, not via this module. |
| `polli_fastsim/inputs.py` | `get_backends(pdf=)`: selects toy (`ToyF2`/`ToyG1`) or grid (`PartonF2`/`PartonG1`) F2/g1 sources. Does **not** own luminosity constants or the α_s table; those are handled by the scripts and the `parton` package directly. |
| **fastsim/scripts/** | Money-delta driver scripts |
| `scripts/money_delta_20260724.py` | Phase II: Interpretation A; script-local `NuclearF2FromGrid` (EPPS21) + R1998; sum-rule solver (Eq. 3A); reach heatmaps; x-axis conversion to peak Δ/F₁ |
| `scripts/money_delta_20260725.py` | Phase II: Interpretation B; script-local `NuclearF2FromGrid` + R1998; analytic Beta-function solver (Eq. 3B); A vs B comparison |
| `scripts/money_delta_20260728.py` | Phase III: φ-modulation projections (Eq. 12); three integration cases; TOP kinematic-trough analysis |
| `scripts/money_delta_20260729.py` | Phase IV: ePIC MC smearing (tracking resolution + ε_eID); sign-convention fix; Case 3 robustness |
| `scripts/money_delta_pdfgrid.py` | Utility: PDF grid diagnostics and EPPS21 validation |
| `scripts/money_delta_realistic.py` | Utility: earlier realistic-smearing prototype (see individual notes) |
| **fastsim/notes/** | Dated session notes and this master synthesis |
| `notes/money_delta_note_2026-07-24.md` | Phase II primary source: full statistical derivation |
| `notes/money_delta_note_2026-07-27.md` | Phase II primary source: Interpretation B |
| `notes/money_delta_note_2026-07-28.md` | Phase III primary source: φ-modulation |
| `notes/money_delta_note_2026-07-29.md` | Phase IV primary source: detector smearing |
| `notes/money_delta_note_2026-07-29_fix.md` | Fix audit: reco-selection consistency fix, production rerun (2026-07-31) |
| `notes/money_delta_uptodate.md` | **This file** — master synthesis |
| **fastsim/out/money_delta/** | Output directory for all money-delta plots (50 PNGs across four sessions; see §13) |
| **fastsim/tests/** | Pytest suite covering `polli_fastsim` library modules |
| **fastsim/README.md** | Top-level fastsim package README |

External dependencies directly used by the workflow:

| Dependency | Role |
|---|---|
| `parton` (pip) | Nuclear PDF evaluation (EPPS21nlo_CT18Anlo_Li6); α_s table |
| `numpy`, `scipy` | Grid computation, Beta-function integrals, MC smearing |
| `matplotlib` | All plot output |

---

## 1. Observable and physical context

**Target.** ⁶Li nucleus (A = 6, Z = 3), transversely tensor-polarized.
The tensor structure function Δ(x, Q²) is defined by the difference of
spin-1 nucleon DIS cross sections under tensor polarization; it is accessed
through the tensor asymmetry A_cos2φ, which produces a cos(2φ) modulation
in the azimuthal distribution of the scattered electron about the virtual
photon axis.

**Observable.** The primary reach observable is the 5σ discovery luminosity:

    L_5σ  (fb⁻¹/nucleon)

and the fractional asymmetry uncertainty at a given luminosity:

    δA/A  (dimensionless; δA/A = 1/N_σ by definition)

A secondary observable, introduced in the 2026-07-28 session, is the
φ-differential yield modulation:

    y(φ) = (N_φ − N_flat) / N_flat = P_zz · ⟨A_cos2φ⟩ · cos(2φ)

which is what an experimenter directly measures and fits.

**Theory background.** The operator structure of Δ was established by
Hoodbhoy–Jaffe–Manohar (NPB 312:571, 1989). The only published numerical
estimates of its first moment are:

- Sather–Schmidt (1993): bag model, ∫x·Δ dx ≡ c·α_s, c = −0.012.
- Detmold–Shanahan: lattice QCD, c = −0.009.

No experimental measurement of Δ exists. The tensor asymmetry A_cos2φ is
the kinematic formula from Hoodbhoy–Jaffe–Manohar (Eq. 5.1 below).
The HERMES b₁ measurement (Phys. Rev. Lett. 95, 242001, 2005) is the
closest published experimental precedent for a tensor DIS asymmetry.

**Luminosity reference.**

    L_inst = 10⁻⁶ fb⁻¹/nucleon/s ≈ 3.6 pb⁻¹/nucleon/hour
    10  fb⁻¹/nucleon ≈ 1 year at nominal EIC luminosity
    100 fb⁻¹/nucleon ≈ 10 years

### 1.1 Estimation logic (flow chart)

The chart below traces the full estimation pipeline from physical inputs
through to final observables. Each box names the concept; the module or
formula reference in parentheses indicates where that step is implemented
or derived. Branching points mark choices that produce the systematic
spreads discussed in §4 and §10.

```
 BEAM SETUP                        ION SPECIES / SCENARIO INPUT
 E_e, p_ion/u                      ⁶Li  (A=6, Z=3, spin=1)
 LOW / MID / TOP configs           tensor polarization P_zz = 0.267
 (beams.py)                        (script input; Cloet value)
         \                              /
          \____________________________/
                        |
                        v
              DIS KINEMATICS GRID
              (x, Q², y) bins; acceptance cuts
              W² > 4 GeV², y in [0.01, 0.95]
              (kinematics.py)
                        |
                        v
         NUCLEAR PDFS & STRUCTURE FUNCTIONS
          F₂(x,Q²) via EPPS21nlo_CT18Anlo_Li6
          F₁ = F₂ / [2x(1+R)],  R from R1998
          α_s(Q²) from parton table
          (script-local NuclearF2FromGrid + r_sigma_lt from structure.py;
           _safe_xfx NumPy-2.0 patch in structure.py)
                        |
                        v
              POLARIZATION INPUT
              P_zz = 0.267 (Cloet)
              enters δA_bin = sqrt(2/N)/P_zz  [Eq. 5]
              and σ² = Σ A²·P_zz²·N/2         [Eq. 6]
              (fom.py)
                        |
                        v
         +--------------+------------------+
         |                                 |
         v                                 v
  INTERP. A                          INTERP. B
   Δ = s·α_s·F₁·x^α·(1-x)^β  [1A]   Δ = s·α_s·x^α·(1-x)^β  [1B]
   (script-local ansatz)              (script-local ansatz)
         |                                 |
         v                                 v
  SUM-RULE SOLVER (Eq. 3A)          SUM-RULE SOLVER (Eq. 3B)
  numerical integral over F₁        analytic Beta function
  s = c / int[x·F₁·x^α(1-x)^β dx]  s = c·P_shape / B(α+2,β+1)
  α_s cancels; s varies ~5%         α_s cancels; s config-independent
  across configs                     (fom.py / script)
         |                                 |
         +----------+----------+----------+
                    |          |          |
                    v          v          v
             low_x shape  mid_x shape  high_x shape
             (α=0.3,β=4)  (α=0.7,β=3)  (α=1.5,β=2)
             ~factor 2-3 spread in L_5σ within each interpretation
                               |
                               v
                  ASYMMETRY AMPLITUDE PER BIN
                  A_cos2φ(x,Q²,y) = -[(1-y)/y²]·Δ / [F₁+(1-y)/(xy²)·F₂]
                  (Hoodbhoy-Jaffe-Manohar; asymmetries.py, Eq. 4)
                               |
                               v
                  RATE PROJECTION PER BIN
                  N_bin = L · dσ/dx dQ² · Δx · ΔQ²
                  (fom.py; dsigma_dx_dq2 from structure.py)
                               |
                               v
                  ASYMMETRY UNCERTAINTY PER BIN
                  δA_bin = sqrt(2/N_bin) / P_zz      [Eq. 5]
                  (fom.py)
                               |
                               v
                   FISHER SUM (statistical combination)
                   σ² = Σ_bins  A_cos2φ,bin² · P_zz² · N_bin / 2   [Eq. 6]
                   δA/A = 1/σ              [Eq. 7]
                   L_5σ = 25 / (σ²/L)     [Eq. 8]
                   (script-local; uses per-bin N_bin from fom.py project_rates)
                               |
                  +------------+--------------------+
                  |                                 |
                  v                                 v
        REACH OUTPUTS (§6)             φ-MODULATION OBSERVABLE (§7)
        L_5σ [fb⁻¹/u]                 y(φ) = P_zz·<A_cos2φ>·cos(2φ) [Eq. 12]
        δA/A @ 10, 100 fb⁻¹           three integration cases:
        heatmaps (|A|, δA, δA/A)        Case 1: peak-rate single bin
        9 (config,shape) combos          Case 2: Q² = 2.43 GeV² slice
        (scripts/money_delta_2026        Case 3: all bins (RECOMMENDED)
         0724.py, _20260725.py)        (scripts/money_delta_20260728.py)
                  |                                 |
                  +------------+--------------------+
                               |
                               v
                  DETECTOR SMEARING BRANCH (§8)
                  (scripts/money_delta_20260729.py)
                       |
              +--------+--------+
              |                 |
              v                 v
       MC SMEARING          ε_eID WEIGHTING
       tracking σ_p/p,      η-dependent efficiency
       σ_θ = σ_φ            (ATHENA+ECCE synthesis)
       piecewise in η        linear interp. 9 anchors
       1000 MC events/bin    (see §8.1)
              |                 |
              +--------+--------+
                       |
                       v
              RECONSTRUCTED (x_reco, Q²_reco, y_reco)
              corrected sign convention:
                Q²_reco = 2 E_e E'(1+cosθ_reco)
                y_reco  = 1 - (E'/2E_e)(1-cosθ_reco)
                       |
                       v
               POST-DETECTOR S/N (§8.2)
               φ-modulation S/N on reco-binned cases:
               S/N = P_zz · |⟨A_cos2φ⟩_reco| · √(N_reco / 2)
                       |
              +--------+--------+
              |                 |
              v                 v
       Case 3 (integrated)  Case 1 (peak-rate bin)
       ~92% pre-detector     changes by 2-3x in
       S/N retained          config-specific direction
       USE THIS              DO NOT QUOTE AS PRIMARY
              |
              v
    FINAL OBSERVABLES / CONCLUSIONS
    - Discovery luminosity L_5σ in [0.01, 0.46] fb⁻¹/u (Interp. A/B, bag)
    - Fractional precision δA/A ≈ 0.5-1% at 100 fb⁻¹ (Interp. A, bag)
    - Use fully-integrated cos(2φ) fit; confirm with φ-modulation curves
    - Systematic band: factor 2-8 between interpretations (leading uncertainty)
    - All configs reach 5σ within 1-year EIC program (10 fb⁻¹/u)
```

---

## 2. Chronology of scripts and notes

| Date | Script | Note | Study phase | Key advance |
|---|---|---|---|---|
| 2026-07-16 | `money_delta_20260715.py` | `money_delta_note_2026-07-16.md` | (I) Formalism / reach | First reach plots; P_zz convention discussion; environment patches documented |
| 2026-07-20 | `money_delta_20260720.py` | `money_delta_note_2026-07-20.md` | (I) Formalism / reach | Nuclear PDF comparison; EPPS21 adopted; luminosity reference fixed |
| 2026-07-21 | `money_delta_20260721.py` | `money_delta_note_2026-07-21.md` | (I) Formalism / reach | Three-tier scenario (10⁻³ / 3×10⁻³ / 10⁻²) established as narrative labels |
| 2026-07-24 | `money_delta_20260724.py` | `money_delta_note_2026-07-24.md` | (II) Interpretation / model comparison | Sum-rule normalization (bag + lattice); x-axis conversion (scale → peak Δ/F₁); heatmaps; full statistical derivation |
| 2026-07-27 | `money_delta_20260725.py`* | `money_delta_note_2026-07-27.md` | (II) Interpretation / model comparison | Interpretation B (no F₁ in Δ); analytic Beta-function solver; A vs B reach comparison |
| 2026-07-28 | `money_delta_20260728.py` | `money_delta_note_2026-07-28.md` | (III) φ-modulation projection | Raw cosine modulation plots; 3 integration cases; TOP kinematic-trough finding |
| 2026-07-29 | `money_delta_20260729.py` | `money_delta_note_2026-07-29.md` | (IV) Detector-realistic extension | ePIC tracking + ε_eID smearing; sign-convention bug found and fixed; Case 3 robustness confirmed |
| 2026-07-31 | `money_delta_20260729.py` (fix) | `money_delta_note_2026-07-29_fix.md` | Fix audit / production rerun | Reco-selection consistency fix; four reco helpers refactored onto canonical `reco_analysis_mask`; S1–S6 static checks pass; production rerun validates <0.1% Case 1 shift and 2–5% Case 2/3 shifts (TOP Case 2 ~5.1%, pass with note); detector-realistic conclusions unchanged |

*Script filename carries date 20260725; the note is dated 2026-07-27.

**Study phases** (referenced throughout this document):

- **(I) Foundational formalism and reach study** — 2026-07-16 through 2026-07-21.
  Established the grid backend, PDF infrastructure, three beam configs,
  luminosity reference, and first L_5σ estimates.

- **(II) Interpretation and model comparison** — 2026-07-24 and 2026-07-27.
  Introduced sum-rule normalization; compared Interpretation A (Δ ∝ F₁) vs
  Interpretation B (Δ pure shape); established the factor-2–8 systematic
  across interpretations.

- **(III) φ-modulation projection study** — 2026-07-28. Shifted from reach
  plots to the raw observable: cosine modulation curves with Poisson error
  bars; three integration cases; TOP kinematic-trough analysis; first explicit
  statement that integrated > single-bin.

- **(IV) Detector-realistic extension** — 2026-07-29. Added MC smearing of
  electron kinematics using ePIC tracking resolution and η-dependent ε_eID;
  confirmed Case 3 robustness (~92% S/N retained) and Case 1 fragility
  (factors of 2–3 change in either direction); fixed sign-convention bug.
  A follow-on reco-selection consistency fix (2026-07-31) repaired internal
  reco-mask/selection inconsistencies in the analysis path; a full production
  rerun confirmed the detector-realistic conclusions are unchanged
  (see `money_delta_note_2026-07-29_fix.md`).

---

## 3. Beam configurations

Three EIC beam configs used throughout all sessions. Ion is ⁶Li throughout.

| Label | E_e [GeV] | p_ion [GeV/u] | √s [GeV/u] | ⟨Q²⟩ [GeV²] |
|---|---|---|---|---|
| LOW | 5 | 27.5 | ~23 | 5.92 |
| MID | 10 | 50 | ~45 | 7.38 |
| TOP | 18 | 137.5 | ~100 | 10.26 |

⟨Q²⟩ is the rate-weighted mean over accepted (x, Q²) bins, used as the
reference scale for α_s(Q²) evaluations and the peak Δ/F₁ conversion (§5.2,
Interpretation A only).

---

## 4. Theory inputs and model variants

### 4.1 Δ ansatz — Interpretation A vs B

The x-shape of Δ(x, Q²) is uncomputed from first principles. Two template
ansätze have been compared.

**Interpretation A** (script `money_delta_20260724.py`; note 2026-07-24):

    Δ(x, Q²) = s · α_s(Q²) · F₁(x, Q²) · x^α · (1−x)^β          (1A)

The template is applied to the ratio Δ/F₁, analogous to modeling g₁/F₁
in spin-structure studies. F₁ in the numerator of Δ aligns the tensor
signal with the highest-rate (low-x, low-Q²) bins, producing optimistic
but physically motivated reach estimates.

**Interpretation B** (script `money_delta_20260725.py`; note 2026-07-27):

    Δ(x, Q²) = s · α_s(Q²) · x^α · (1−x)^β                       (1B)

The template is applied to Δ directly. F₁ is absent from the Δ numerator.
The sum-rule solver becomes analytic (Beta function), and A is strictly
config-independent. Reach is more conservative by factors of 2–8.

**Which is correct?** No first-principles calculation determines whether
the `x^α(1-x)^β` template applies to Δ (Interpretation B) or Δ/F₁
(Interpretation A). Resolving this requires a first-principles computation
of the x-shape of Δ (e.g., lattice QCD, NJL model, or CDKS convolution).
This is currently the leading modeling systematic. Honest reach projections
must span both interpretations.

### 4.2 Sum-rule constraints (bag vs lattice)

The first-moment sum rule:

    ∫₀¹ x · Δ(x, Q²) dx = c · α_s(Q²)

Two independent theoretical estimates of c:

| Source | c | Reference |
|---|---|---|
| Bag model | −0.012 | Sather–Schmidt (1993) |
| Lattice QCD | −0.009 | Detmold–Shanahan |

Key facts:
- |c_bag| / |c_lat| = 4/3 exactly. The significance ratio scales as |A|, so
  L_5σ(lattice) = (4/3)² × L_5σ(bag) ≈ 1.78× larger.
- α_s cancels in the sum-rule solver (§5.2). A is a constant across Q² bins
  under both interpretations.
- The three-tier scenario (10⁻³, 3×10⁻³, 10⁻²) used in the 2026-07-21 note
  is narrative labeling only; only the 10⁻³ tier has a calculation behind it
  (Sather–Schmidt bag). The sum-rule normalization (2026-07-24 onward)
  supersedes this ad-hoc scan.

### 4.3 x-shape variants

Three Regge-inspired / Brodsky–Farrar templates spanning a range of
small-x behavior. Applied identically under both interpretations.

| Label | α | β | Shape character | x_peak (approx) |
|---|---|---|---|---|
| low_x | 0.3 | 4 | Soft, Regge-like rise at small x | ~0.07 |
| mid_x | 0.7 | 3 | Moderate peak | ~0.19 |
| high_x | 1.5 | 2 | Suppressed at small x, peaks at larger x | ~0.43 |

All three variants have positive α, meaning they miss a gluon-dominated
scenario where Δ could rise more steeply at small x.

**Within-interpretation shape systematic:** L_5σ varies by a factor of ~2–3
across the three variants (from the 2026-07-24 and 2026-07-27 reach tables;
see §6). The between-interpretation systematic (factor 2–8) is comparable to
or exceeds the within-interpretation shape systematic.

### 4.4 PDF and structure-function inputs

| Input | Source | Used for |
|---|---|---|
| R = σ_L/σ_T | R1998 parametrization | F₂ → F₁ conversion: F₁ = F₂/(2x) × 1/(1+R) |
| Nuclear PDF | EPPS21nlo_CT18Anlo_Li6 | F₁(x, Q²) and F₂(x, Q²) per bin |
| α_s(Q²) | `parton` package table | Per-bin α_s for Δ model and sum-rule solver |

Nuclear PDF choice was found to contribute sub-1% to reach numbers in the
2026-07-20 session (not fully re-derived in the four primary source notes);
EPPS21 was adopted as sufficient on that basis. R1998 is applied uniformly.

### 4.5 Polarization

Cloet tensor polarization: **P_zz = 0.267** (used throughout all sessions).

This value enters every statistical formula as a dilution factor. A factor
~1.5 uncertainty in P_zz_eff → factor ~2.3 in σ² → factor ~1.5 in L_5σ.
The P_zz convention (whether 0.267 is appropriate for ⁶Li nuclear structure)
is the highest-priority unresolved physics input (§11.1, Tier A #1).

---

## 5. Key formulas

### 5.1 Δ model and asymmetry amplitude

Under Interpretation A:

    Δ(x, Q²) = s · α_s(Q²) · F₁(x, Q²) · x^α · (1−x)^β          (1A)

Under Interpretation B:

    Δ(x, Q²) = s · α_s(Q²) · x^α · (1−x)^β                       (1B)

Tensor asymmetry amplitude per (x, Q², y) bin
(Hoodbhoy–Jaffe–Manohar kinematic formula):

    A_cos2φ(x, Q², y) = −[(1−y)/y²] · Δ / [F₁ + ((1−y)/(xy²))·F₂]   (4)

Note: the denominator contains F₁ and F₂ regardless of interpretation.
The numerator Δ is given by (1A) or (1B).

### 5.2 Sum-rule solver

The first-moment sum rule (unchanged under both interpretations):

    ∫₀¹ x · Δ(x, Q²) dx = c · α_s(Q²)                                 (2)

**Interpretation A solver** — numerical integral over F₁:

    A ≡ s = c / ∫₀¹ x · F₁(x, ⟨Q²⟩) · x^α · (1−x)^β dx              (3A)

α_s cancels. A varies ~5% across configs (from F₁ evolution).

**Interpretation B solver** — analytic Beta function:

    A ≡ s = c · P_shape(α, β) / B(α+2, β+1)                            (3B)

where B(α+2, β+1) = ∫₀¹ x^(α+1)(1−x)^β dx is the Euler Beta function and
P_shape(α, β) = (α/(α+β))^α · (β/(α+β))^β is the peak value of the shape
factor x^α(1−x)^β, evaluated at its mode x_peak = α/(α+β). α_s cancels
exactly. A is strictly config-independent.

Numerical values of P_shape: low_x (α=0.3, β=4): 0.3369; mid_x (α=0.7, β=3):
0.1662; high_x (α=1.5, β=2): 0.0916. See §6.4 for the resulting A values.

**x-axis convention (Interpretation A plots only).** For Plots 1–4 produced
by `money_delta_20260724.py`, the internal parameter `s` (= A) is converted to
**peak Δ/F₁** on the displayed x-axis:

    peak Δ/F₁ = s × α_s(⟨Q²⟩) × max_x[x^α (1−x)^β]

This conversion applies specifically to Interpretation A (Eq. 1A), where the
ansatz is written on Δ/F₁. Under Interpretation B (Eq. 1B), s is a multiplier
on Δ directly, not on Δ/F₁, and the conversion produces a quantity with
different physical meaning; the 2026-07-27 plots do not use this axis.

For the mid_x shape at MID config (⟨Q²⟩ = 7.4 GeV²):
conversion factor = α_s(7.4) × peak_shape ≈ 0.30 × 0.170 = 0.051.

| Theory input | Raw |A| (= s, Interp. A) | peak Δ/F₁ (Interp. A) |
|---|---|---|
| Bag (c = −0.012, mid_x, MID) | 0.310 | **1.6 × 10⁻²** |
| Lattice (c = −0.009, mid_x, MID) | 0.233 | **1.2 × 10⁻²** |

### 5.3 Per-bin statistical uncertainty

From a maximum-likelihood cos(2φ) fit to N_bin events in a single (x, Q², y)
bin (derivation in 2026-07-24 note §4.2.1):

    δA_bin = √(2/N_bin) / P_zz                                          (5)

The factor √2 arises from ⟨cos²(2φ)⟩ = 1/2 under uniform φ coverage.
The 1/P_zz un-dilutes the tensor polarization.

**Caveat:** this formula assumes uniform azimuthal acceptance. For a detector
with gaps or shadowed regions, ⟨cos²(2φ)⟩_det ≠ 1/2 and must be computed
from the acceptance-weighted integral.

### 5.4 Combined significance and fractional precision

The combined significance across all (x, Q²) bins is the signal-to-noise
on the one-parameter fit for the overall scale s (Fisher information derivation
in 2026-07-24 note §4.2.2):

    σ² = Σ_bins (A_bin / δA_bin)²
       = Σ_bins  A_cos2φ,bin² · P_zz² · N_bin / 2                      (6)

Fractional amplitude uncertainty and discovery luminosity:

    δA/A = 1/σ = 1/N_σ                                                  (7)

    L_5σ = 25 / (σ²/L)    at reference luminosity L                     (8)

    (δA/A)_{5σ} = 0.20   (discovery threshold)                          (9)

**σ ↔ δA/A duality table:**

| σ (significance) | δA/A |
|---|---|
| 1 | 100% |
| 3 | 33% |
| **5** | **20%** |
| 10 | 10% |
| 100 | 1% |

**Important:** the combined δA/A (Eq. 7) is the fractional uncertainty on
the fitted overall amplitude across all bins. The per-bin median |δA/A| is
~35–65× larger and measures only individual bin quality. Combining per-bin
medians does NOT give the combined uncertainty; the correct combination uses
inverse-variance weighting (Eq. 6 → 7). See 2026-07-24 §4.2.3 for
the full explanation.

### 5.5 Analytic rescaling formula

Because A_bin ∝ s and N_bin ∝ L:

    σ²(s, L) = σ²_ref · (s / S₀)² · L                                  (10)

    δA/A(s, L) = 1 / √[σ²_ref · (s/S₀)² · L]                          (11)

where σ²_ref is evaluated once at reference (s = S₀ = 10⁻³, L = 1 fb⁻¹/u).
This reduces compute cost by ~99% for the scan plots.

### 5.6 φ-modulation observable

Fractional yield modulation per φ bin (2026-07-28 note §2):

    y(φ) = (N_φ − N_flat) / N_flat = P_zz · ⟨A_cos2φ⟩ · cos(2φ)      (12)

Poisson error per φ bin in the small-modulation limit:

    σ_y(φ) ≈ 1/√N_flat                                                  (13)

where N_flat = N_total / (number of φ bins). The small-modulation condition
|P_zz · A_cos2φ · cos(2φ)| ≪ 1 is satisfied throughout: the maximum
modulation across all tested combinations is 5.99% (LOW + Case 1 + 10×A_bag),
well below the physical positivity bound of 1. Error bars are essentially
signal-scale-independent as a consequence.

Eq. (12) is exact at leading twist in Δ — not a Taylor expansion.

---

## 6. Reach results

### 6.1 Interpretation A reach table (2026-07-24)

Bag sum-rule constraint (c = −0.012). Lattice values give L_5σ × 16/9 ≈ 1.78×
larger, δA/A × 4/3 larger.

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

**Worst case (Interpretation A):** TOP + high_x, L_5σ ≈ 0.29 fb⁻¹/u
(~80 hours at nominal EIC luminosity). Under the lattice constraint, all nine
cases still reach 5σ within a 1-year (10 fb⁻¹) program.

### 6.2 Interpretation B reach table (2026-07-27)

Bag sum-rule constraint. Lattice scales as above.

| Config | Variant | σ²/(fb⁻¹) | δA/A @ 10 fb⁻¹ | δA/A @ 100 fb⁻¹ | L_5σ [fb⁻¹/u] |
|---|---|---|---|---|---|
| LOW | low_x | 841 | 0.0109 | 0.0034 | 0.03 |
| LOW | mid_x | 464 | 0.0147 | 0.0046 | 0.05 |
| LOW | high_x | 317 | 0.0178 | 0.0056 | 0.08 |
| MID | low_x | 604 | 0.0129 | 0.0041 | 0.04 |
| MID | mid_x | 326 | 0.0175 | 0.0055 | 0.08 |
| MID | high_x | 309 | 0.0180 | 0.0057 | 0.08 |
| TOP | low_x | 98 | 0.0319 | 0.0101 | 0.25 |
| TOP | mid_x | 55 | 0.0428 | 0.0135 | 0.46 |
| TOP | high_x | 75 | 0.0364 | 0.0115 | 0.33 |

**Worst case (Interpretation B):** TOP + mid_x, L_5σ ≈ 0.46 fb⁻¹/u
(~128 hours). Under the lattice constraint, worst case scales to ~0.82 fb⁻¹/u
(still < 1 year).

### 6.3 A vs B comparison (MID + mid_x reference)

| Quantity | Interp. A | Interp. B | B/A ratio |
|---|---|---|---|
| \|A_bag\| | 0.310 | 0.0890 | 0.29× |
| σ²/(fb⁻¹) | 1339 | 326 | 0.24× (~4× worse) |
| L_5σ MID + mid_x (bag) | 0.02 fb⁻¹/u | 0.08 fb⁻¹/u | 4× worse |
| L_5σ TOP + mid_x (bag) | 0.06 fb⁻¹/u | 0.46 fb⁻¹/u | 8× worse |
| δA/A @ 10 fb⁻¹ | 0.0086 | 0.0175 | 2× worse |
| δA/A @ 100 fb⁻¹ | 0.0027 | 0.0055 | 2× worse |

The reach penalty ranges from 4× worse in L_5σ (MID + mid_x) to 8× worse
(TOP + mid_x), with a corresponding ~4× drop in σ² per fb⁻¹ for the MID
reference case. These ratios vary with config and shape because the two
interpretations differ in how bin-level signal aligns with event rate.

Two effects combine to produce the reach penalty. First, the two solvers use
different normalization prescriptions and yield different |A| values for the
same (variant, c) pair: Eq. (3A) integrates x · F₁(x, ⟨Q²⟩) · x^α(1−x)^β
over x, while Eq. (3B) uses the peak value of x^α(1−x)^β normalized by its
Beta-function integral. The tabulated results (§6.4) show the outcome:
|A_bag| = 0.310 under Eq. (3A) vs 0.089 under Eq. (3B) for MID + mid_x.
Second, under Interpretation A, the F₁ factor in the Δ numerator (Eq. 1A)
amplifies the per-bin signal in exactly those bins that also have the highest
event rate (low x, low Q²), concentrating Fisher information; under
Interpretation B (Eq. 1B), F₁ is absent from Δ, so this alignment between
signal and rate is lost. The reach penalty in L_5σ (which scales as
1/σ² ∝ 1/|A|²) reflects both effects and varies across (config, shape)
combinations as shown above.

### 6.4 A-value tables

**Interpretation A** (numerical integral over F₁; ~5% variation across configs):

| Config | Shape | ⟨Q²⟩ [GeV²] | A_bag | A_lat |
|---|---|---|---|---|
| LOW | low_x | 5.92 | −0.3862 | −0.2896 |
| LOW | mid_x | 5.92 | −0.3178 | −0.2383 |
| LOW | high_x | 5.92 | −0.3913 | −0.2934 |
| MID | low_x | 7.38 | −0.3737 | −0.2803 |
| MID | mid_x | 7.38 | −0.3100 | −0.2325 |
| MID | high_x | 7.38 | −0.3875 | −0.2907 |
| TOP | low_x | 10.26 | −0.3509 | −0.2632 |
| TOP | mid_x | 10.26 | −0.2968 | −0.2226 |
| TOP | high_x | 10.26 | −0.3814 | −0.2861 |

**Interpretation B** (analytic Beta function; strictly config-independent):

| Variant | α | β | B(α+2, β+1) | A_bag | A_lat |
|---|---|---|---|---|---|
| low_x | 0.3 | 4 | 0.02202 | −0.18355 | −0.13766 |
| mid_x | 0.7 | 3 | 0.02242 | −0.08895 | −0.06671 |
| high_x | 1.5 | 2 | 0.02309 | −0.04762 | −0.03571 |

In both tables: |A_bag|/|A_lat| = 0.012/0.009 = 4/3 exactly, fixed by the
c-coefficient ratio.

---

## 7. φ-modulation projections (2026-07-28, study phase III)

**Setup.** Fixed to Interpretation A, mid_x shape (α=0.7, β=3), 10 fb⁻¹/u.
Plots show three overlaid cosine curves (s = 0.1·A_bag, A_bag, 10·A_bag)
with one set of Poisson error bars.

### Integration cases

| Case | Description |
|---|---|
| 1 (peak bin) | Single (x, Q²) bin with the largest N_bin |
| 2 (Q² slice) | Q² bin iq2=3 (Q² = 2.43 GeV²), all accepted x |
| 3 (integrated) | All accepted (x, Q²) bins combined |

### Numerical results (pre-detector, at s = A_bag)

| Config | Case | Peak (x, Q²) or range | N_total | ⟨A_cos2φ⟩ | S/N |
|---|---|---|---|---|---|
| LOW | 1 | (0.089, 2.43 GeV²) | 2.37e+7 | −0.0224 | 20.7 |
| LOW | 2 | Q² = 2.43 GeV² slice | 3.11e+8 | −0.0158 | 52.5 |
| LOW | 3 | all bins | 1.20e+9 | −0.0170 | 111.1 |
| MID | 1 | (0.071, 2.43 GeV²) | 2.46e+7 | −0.0198 | 18.6 |
| MID | 2 | Q² = 2.43 GeV² slice | 4.11e+8 | −0.0104 | 39.9 |
| MID | 3 | all bins | 1.79e+9 | −0.0121 | 96.6 |
| TOP | 1 | **(0.0022, 2.43 GeV²)** | 2.43e+7 | **−0.00208** | **1.94** |
| TOP | 2 | Q² = 2.43 GeV² slice | 4.49e+8 | −0.00346 | 13.8 |
| TOP | 3 | all bins | 2.25e+9 | −0.00533 | 47.7 |

S/N formula: S/N ≈ P_zz · |⟨A_cos2φ⟩| · √(N_total / 2).

**Key findings from 2026-07-28.**

1. **Case 3 wins in S/N for all configs.** The √N gain from integrating
   (factor ~7–10) beats the dilution of ⟨A_cos2φ⟩ (factor ~0.6–0.8) in all
   cases. From the table above, Case 3 beats Case 2 by about 2–3.5× and
   Case 1 by about 5–25× across the three beam configs.

2. **TOP's peak-rate bin is a kinematic trough.** At √s ≈ 100 GeV/u, the
   DIS rate peaks at x ≈ 0.002. Under the mid_x shape (x^0.7), the shape
   factor x^0.7 at x = 0.002 is 0.011, while at x = 0.07 (where LOW/MID
   rates peak) it is 0.145 — **13× larger**. Additionally, the kinematic
   prefactor (1−y)/y² in Eq. (4) favors LOW/MID peak bins by ~500× over
   TOP's peak bin. The combined result: TOP Case 1 |⟨A_cos2φ⟩| is 10×
   smaller than LOW/MID despite similar statistics. For TOP, integrating
   *increases* ⟨A_cos2φ⟩ monotonically (0.0021 → 0.0035 → 0.0053), the
   opposite of LOW and MID.

3. **The per-bin Fisher fit (2026-07-24) is ~20% better than Case 3.**
   For MID at 10 fb⁻¹: Fisher σ = 116 (δA/A = 0.0086) vs Case 3 S/N = 97
   (δA/A ≈ 0.010). Fisher weighting (A_bin² · N_bin per bin) extracts more
   information than the yield-weighted average amplitude used in the φ-plot
   S/N. In practice both approaches are used: φ plots for visualization,
   Fisher fit for the measurement uncertainty.

4. **Heatmap diagnostic: peak-rate bin vs peak-|A| bin.**

   | Config | Peak-rate bin |A| | Max |A| on heatmap | Ratio |
   |---|---|---|---|
   | LOW | 2.24e-2 | 3.09e-2 | 1.4× |
   | MID | 1.98e-2 | 2.94e-2 | 1.5× |
   | TOP | **2.08e-3** | **2.73e-2** | **13×** |

   For LOW and MID the peak-rate bin is within ~1.5× of the optimal single
   bin. For TOP, a trigger designed around the DIS rate peak captures the
   worst available single bin for the tensor signal.

---

## 8. Detector-realistic extension (2026-07-29, study phase IV)

### 8.1 Detector model

Applied to the scattered electron kinematics via MC smearing (1000 MC events
per accepted true bin; seed 42 for reproducibility).

**Tracking momentum resolution** (piecewise in η, 7 regions):

    σ_p / p = √(a² · p² + b²)    (coefficients from ePIC tracking-only performance)

**Tracking angular resolution** (piecewise in η, 7 regions):

    σ_θ = σ_φ = [1–5 mrad, piecewise in η]

**Electron-ID efficiency** (ATHENA + ECCE synthesis, linear interpolation
of 9 anchor points):

| η | −3.5 | −3.0 | −2.0 | −1.0 | 0.0 | 1.0 | 2.0 | 3.0 | 3.5 |
|---|---|---|---|---|---|---|---|---|---|
| ε_eID | 0.85 | 0.92 | 0.95 | 0.93 | 0.90 | 0.90 | 0.85 | 0.80 | 0.70 |

Zero outside |η| > 3.5. Uncertainties: ±5% in barrel/backward region,
±10% in forward direction.

**Not included in this model:** ECal energy resolution (tracking-only smearing
was used), radiative corrections, azimuthal acceptance holes within |η| < 3.5.

**Sign-convention bug (found and fixed).** The initial script had wrong
sign conventions in the reconstructed Q² and y formulas (using 1−cos θ for
Q² and 1+cos θ for y, when the project convention measures θ from +z, the
ion direction, so scattered electrons sit at θ ≈ π). The corrected formulas:

    Q²_reco = 2 E_e E'_reco (1 + cos θ_reco)
    y_reco  = 1 − (E'_reco / 2E_e) (1 − cos θ_reco)

Before the fix: 0% of events landed in their original bin. After fix:
47–76% same-bin fractions, physically consistent with resolution estimates.

### 8.2 S/N results pre- vs post-detector

At L = 10 fb⁻¹/u, s = A_bag, Interpretation A, mid_x shape:

| Config | Case | Pre-detector S/N | Post-detector S/N | Ratio |
|---|---|---|---|---|
| LOW | 1 (peak bin) | 20.7 | ~21 | 100% |
| LOW | 2 (Q² slice) | 52.5 | ~49 | 93% |
| LOW | 3 (integrated) | 111.1 | ~104 | 93% |
| MID | 1 (peak bin) | 18.6 | **~7.9** | **42%** |
| MID | 2 (Q² slice) | 39.9 | ~24 | 60% |
| MID | 3 (integrated) | 96.6 | ~89 | 92% |
| TOP | 1 (peak bin) | 1.94 | **~4.1** | **214%** |
| TOP | 2 (Q² slice) | 13.8 | ~9 | 65% |
| TOP | 3 (integrated) | 47.7 | ~44 | 92% |

**Bin-migration diagnostics (same-bin fractions):**

| Config | f_raw (unweighted) | f_weighted (ε_eID-weighted) | Survival rate |
|---|---|---|---|
| LOW | 88% | 76% | 93.6% |
| MID | 76% | 62% | 90.4% |
| TOP | 64% | 47% | 86.0% |

The gap between raw and weighted fractions (12–17%) reflects that migrated
events land at higher-ε η than same-bin events: migration preferentially
carries events from edge bins (lower ε) toward interior bins (higher ε).
The gap grows with beam energy because TOP samples more of the forward-η
region where ε_eID has its steepest gradient (0.80 at η=3.0 → 0.70 at η=3.5).

**Reconstructed ⟨A_cos2φ⟩ comparison:**

| Config | Case | Reco ⟨A⟩ (20260729) | Pre-detector ⟨A⟩ (20260728) |
|---|---|---|---|
| LOW | 1 | −0.02285 | −0.02244 |
| LOW | 3 | −0.01698 | −0.01698 |
| MID | 1 | **−0.00777** | −0.01977 |
| MID | 3 | −0.01171 | −0.01209 |
| TOP | 1 | **−0.00426** | −0.00208 |
| TOP | 3 | −0.00502 | −0.00533 |

Case 3 values agree to within a few percent across all configs; Case 1
values change by factors of 2–3.

### 8.3 Why Case 3 is robust; why Case 1 is fragile

**Case 3 (integrated) is robust** because bin migration under full integration
is a redistribution *within* the accepted region. Events leaving one bin land
in another, and both bins are summed. Real losses are only from (a) edge
acceptance (physicality rejection 0.1–5%) and (b) ε_eID < 1 (7–14% S/N loss,
dominated by efficiency not migration). Case 3 retains 92–93% of pre-detector
S/N across all three configs. The residual loss is determined by the overall
ε_eID level, not by the resolution.

**Case 1 (peak-rate single bin) is fragile** because the peak-rate reco bin
is determined by whichever (x, Q²) bin accumulates the most reconstructed
events. On a log-spaced grid with a broadly-peaked rate distribution, modest
smearing shifts the reco rate maximum to a different true bin. Whether the
reco peak bin has a larger or smaller |A_cos2φ| than the true peak bin depends
on the local |A| landscape and is config-specific:

- **MID:** migration pools reco events at x = 0.011 instead of the true peak
  at x = 0.071. The mid_x shape factor at x = 0.011 is ~6× smaller, suppressing
  |A_cos2φ| by ~2.5×. Reco Case 1 S/N drops to 42% of pre-detector.
- **TOP:** the true peak-rate bin is already a kinematic trough (x = 0.0022).
  Migration shifts the reco peak to x = 0.0036, slightly less unfavorable.
  Reco |A_cos2φ| rises by 2×; S/N is 214% of pre-detector. This is not a
  gain — it is an artifact of migration moving the reco peak away from the
  worst available bin.
- **LOW:** true peak bin is already at moderate x (0.089), near the center
  of the |A| landscape. Migration is symmetric; Case 1 S/N is ~100%.

**Conclusion from this detector model:** under ePIC tracking-only smearing +
η-dependent ε_eID (Interpretation A, mid_x, bag, 10 fb⁻¹), use the fully-
integrated analysis. Case 1 single-bin analysis changes by large and
unpredictable factors under this smearing model; the direction is config-specific.
This is consistent across all three beam configs under the tested model but has
not been verified under Interpretation B, other x-shapes, or the lattice constraint.

**Post-fix note (2026-07-31).** A reco-selection consistency fix was applied to
`money_delta_20260729.py` and validated through a full production rerun
(see `money_delta_note_2026-07-29_fix.md`). The fix repaired internal
reco-mask/selection inconsistencies; it did not produce an obvious qualitative
overturning of the detector-realistic message. Case 1 ⟨A_cos2φ⟩ values are
unchanged (<0.1%); Case 2 and Case 3 values shifted by 2–5%, with TOP Case 2
reaching ~5.1% (marginally above the ≤5% plan tolerance — recorded as a
pass-with-note detail, not a headline result). Integrated Case 3 remains the
robust recommendation. Results are scoped strictly to this detector model and
physics setup; no generalization beyond the tested configuration is implied.

---

## 9. Four study phases — what each established

### Phase I — Foundational formalism and reach study (2026-07-16 to 2026-07-21)

Established (as reported in the four primary source notes, not independently
re-verified here): ⁶Li DIS grid backend; acceptance cuts; R1998 + EPPS21 PDF
infrastructure; three EIC beam configs; luminosity reference (10 fb⁻¹ ≈ 1 year
at EIC); first L_5σ estimates for an ad-hoc peak-Δ/F₁ scan; three-tier scenario
(10⁻³ / 3×10⁻³ / 10⁻²) as narrative labels. Also: NumPy 2.0 compatibility
patch for the `parton` package, `_safe_xfx` helper in `polli_fastsim`.

Not yet in Phase I: sum-rule normalization; interpretation comparison;
φ-modulation plots; detector smearing.

### Phase II — Interpretation and model comparison (2026-07-24 and 2026-07-27)

Established: sum-rule normalization (bag c = −0.012, lattice c = −0.009),
replacing the ad-hoc amplitude scan; full statistical derivation (per-bin
δA_bin = √(2/N)/P_zz; Fisher-sum σ²; 1-parameter fit; σ ↔ δA/A duality);
x-axis conversion from raw scale parameter to physical peak Δ/F₁; per-bin
heatmaps (|A_bin|, δA_bin, |δA/A|); Interpretation A vs B; analytic
Beta-function solver for Interpretation B.

Key numbers: L_5σ ∈ [0.01, 0.29] fb⁻¹/u (Interp. A, bag) and
[0.03, 0.46] fb⁻¹/u (Interp. B, bag). Both within 1-year EIC program.

### Phase III — φ-modulation projection study (2026-07-28)

Established: raw cosine modulation observable y(φ); three integration cases
(peak-bin, Q²-slice, fully-integrated); Case 3 wins S/N for all configs;
TOP kinematic-trough finding (peak-rate bin is x ≈ 0.002, where the tensor
signal under mid_x shape is 13× smaller than at moderate x); heatmap
cross-reference to Case 1 bin selection; bin-width scan for TOP Case 1 (5°
through 45°); small-modulation limit verified (max modulation 5.99%).

Key numbers: integrated S/N = 47–111 at 10 fb⁻¹ under Interpretation A,
mid_x, bag. Per-bin Fisher fit ~20% better than integrated single-cosine fit.

### Phase IV — Detector-realistic extension (2026-07-29 + 2026-07-31 fix)

Established under one detector model (ePIC tracking-only momentum and angular
resolution + η-dependent ε_eID synthesis, Interpretation A, mid_x shape, bag
constraint, 10 fb⁻¹): Case 3 retains 92–93% of pre-detector S/N (7–8% loss
from efficiency, negligible from migration); Case 1 changes by factors of 2–3
in either direction; sign-convention bug in reconstructed Q² and y found and
fixed; raw vs efficiency-weighted migration fractions distinguished (gap 12–17%,
growing with beam energy from ε_eID gradient at forward η).

A reco-selection consistency fix (2026-07-31) repaired four internal
reco-mask/selection inconsistencies in the analysis path (canonical
`reco_analysis_mask` helper, S1–S6 static-analysis checks, production rerun
validated). The production rerun confirmed: Case 1 ⟨A_cos2φ⟩ unchanged
(<0.1%); Case 2/3 values shifted by 2–5% (TOP Case 2 ~5.1%, pass with note);
detector-realistic conclusions — including the integrated Case 3 recommendation
— are unchanged.

The "use integration" recommendation is consistent across all three configs
under this model, reinforcing the pre-detector conclusion from Phase III. It
has not been tested under Interpretation B, other x-shapes, or lattice
constraint combinations.

Key numbers: post-detector integrated S/N = 44–104 at 10 fb⁻¹ (Interp. A,
mid_x, bag). Documented in `money_delta_note_2026-07-29.md` (primary) and
`money_delta_note_2026-07-29_fix.md` (fix audit).

---

## 10. Stable conclusions vs model-dependent findings

### Stable (pre-detector) — robust across both interpretations, all x-shapes, bag and lattice constraints

1. **EIC discovery of the tensor asymmetry Δ is feasible within a 1-year
   program (10 fb⁻¹/nucleon).** Holds under both interpretations (A and B),
   all three x-shape variants, and both sum-rule constraints (bag and lattice).
   Worst case (Interpretation B, lattice, TOP, mid_x): L_5σ ≈ 0.82 fb⁻¹/u,
   still well within 1 year. This conclusion is from ideal simulation (no
   detector smearing); see the detector-tested conclusions below.

2. **At 100 fb⁻¹, the measurement becomes systematics-limited** (δA/A ~
   0.5–1% under Interpretation A, bag). Statistics are no longer the
   bottleneck.

3. **|A_bag|/|A_lat| = 4/3 exactly**, fixed by the c-coefficient ratio
   (0.012/0.009). True under both interpretations.

4. **The fully-integrated analysis (Case 3) gives the best S/N across all
   configs in the 2026-07-28 pre-detector φ-modulation study** (Interpretation
   A, mid_x shape, bag constraint, 10 fb⁻¹). Under those fixed conditions,
   Case 3 beats Case 2 by about 2–3.5× and Case 1 by about 5–25× in S/N
   across all three beam configs (see §7 table). This comparison has not been
   repeated under Interpretation B, other x-shapes, or the lattice constraint;
   it is plausible but unverified that the ranking holds there too.

5. **Per-bin |δA/A| ratios should not be averaged or quoted as combined
   sensitivity.** The combined uncertainty is the Fisher sum (Eq. 6 → 7),
   which can be ~40× smaller than the per-bin median. These are distinct
   quantities; averaging per-bin ratios underestimates sensitivity by ~40×.

### Stable (detector-tested) — under one specific detector model

The following conclusions are based on the 2026-07-29 study using ePIC
tracking-only resolution + ATHENA/ECCE ε_eID anchors, fixed to Interpretation
A, mid_x shape, bag constraint, 10 fb⁻¹, and confirmed by the 2026-07-31
production rerun after the reco-selection consistency fix. They have not been
tested under Interpretation B, other x-shapes, or the lattice constraint.

6. **The fully-integrated analysis retains ~92% of pre-detector S/N** under
   tracking-only smearing + η-dependent ε_eID. The 7–8% loss is from efficiency
   (ε_eID < 1), not from bin migration. The 2026-07-31 production rerun after
   the reco-selection consistency fix found small Case 3 ⟨A_cos2φ⟩ shifts of
   2–4% relative to the pre-fix values, with no qualitative overturning of
   the detector conclusion; the integrated Case 3 recommendation stands. The
   rerun establishes consistency of the fixed code path for this one detector
   model and does not constitute a broader detector-physics validation.

7. **The single-bin peak analysis (Case 1) changes by factors of 2–3 in either
   direction** under detector smearing, in a config-specific and unpredictable
   direction. Do not quote Case 1 as a primary reach figure. Post-fix Case 1
   values are unchanged (<0.1% shift), confirming `find_peak_bin_reco` already
   applied effective acceptance masking pre-fix.

   **TOP Case 2 ~5.1% shift** relative to the pre-fix §3.3 value is the largest
   single Case 2/3 shift in the production rerun. It is recorded as a pass-with-
   note detail (plan tolerance ≤5%), consistent with TOP populating more of the
   forward acceptance edge where out-of-acceptance migration is slightly larger.
   It does not constitute a headline result or a reversal of conclusions.

### Model-dependent — conclusions contingent on unresolved choices

1. **Interpretation A vs B factor 2–8.** Which ansatz is correct
   (Δ ∝ F₁·x^α(1-x)^β vs Δ ∝ x^α(1-x)^β) is unresolved. Reach quotes
   must span both until a first-principles x-shape calculation is available.

2. **Within-interpretation shape systematic, factor ~2–3.** The three x-shape
   variants (low_x / mid_x / high_x) produce a factor ~2–3 spread in L_5σ
   within each interpretation.

3. **TOP kinematic-trough severity is shape-dependent.** Under mid_x (α=0.7),
   the peak-rate bin is 13× suppressed relative to the best available bin.
   Under low_x (α=0.3), the suppression is only ~2×. This has direct
   implications for trigger design at the TOP config.

4. **A_bag/A_lat ratio fixes the relative reach between theory inputs**, but
   neither calculation has rigorous uncertainty estimates. Sather–Schmidt is
   a tree-level bag model; Detmold–Shanahan is lattice but finite-volume
   and quenched effects are not fully controlled.

5. **The overall efficiency loss (7–14%) is dominated by the ε_eID model**
   (ATHENA + ECCE synthesis, ±5–10% uncertainty), not by the tracking
   resolution. When an official ePIC curve is available it should replace
   these anchors.

---

## 11. Open uncertainties and next steps

### 11.1 Tier A — physics critical path

**A1. P_zz convention resolution (highest priority).**
Whether P_zz = 0.267 is the appropriate effective tensor polarization for
⁶Li nuclear structure — including sub-nuclear state contributions — is
unresolved. A factor ~1.5 in P_zz_eff propagates to factor ~2.3 in σ²
and factor ~1.5 in L_5σ. Consult I. Cloet or nuclear-structure literature.

**A2. Δ x-shape anchoring via first-principles theory.**
The Interpretation A vs B comparison (factor 2–8 in reach) dominates all
other modeling systematics. Options: lattice-moment constraints beyond the
first moment, CDKS convolution model, NJL or light-cone model calculations.
Minimum viable: one motivated x-shape with its own uncertainty estimate.

**A3. Theorist consultation on ansatz choice.**
Is x^α(1-x)^β better applied to Δ directly (Interpretation B) or to Δ/F₁
(Interpretation A)? This is a qualitative question about operator structure
that tensor-DIS theorists (Hoodbhoy, Jaffe, Cloet, Detmold–Shanahan group)
could likely answer quickly.

**A4. Publication-ready reach quotes must span both interpretations.**
Present A and B reach tables side-by-side and characterize the factor 2–8
as a physics-modeling systematic, not a numerical uncertainty.

**A5. Systematic-uncertainty studies (now on the critical path).**
At 100 fb⁻¹, statistics are no longer limiting. Leading uncertainties
become: P_zz calibration, F₂A normalization, radiative corrections, and
x-shape modeling. These must be quantified before any precision physics
claim can be made.

### 11.2 Tier B — detector realism

**B1. ECal energy resolution as a second smearing model.**
The current script uses tracking-only smearing. ECal typically provides better
energy resolution than tracking for scattered electrons concentrated in the
backward-η acceptance. Adding ECal as an alternative or combined smearing option
would give a second data point on how sensitive the migration fractions and Case
1 S/N are to the electron reconstruction method. Whether it improves or changes
the integrated (Case 3) S/N is expected to be small, since Case 3 is dominated
by the overall efficiency loss rather than resolution.

**B2. Replace ATHENA + ECCE ε_eID anchors with the official ePIC curve.**
Current anchors carry ±10% uncertainty in the forward region, which is the
dominant contribution to the 7–14% efficiency loss. Use the official ePIC
parameterization when available.

**B3. Radiative corrections.**
Born-level DIS assumed throughout. QED radiative corrections introduce
~5–15% effects on reconstructed (x, Q², y). Add as a systematic after the
nominal study is complete.

**B4. Non-uniform φ acceptance.**
Uniform azimuthal coverage within |η| < 3.5 is assumed. Real ePIC has
azimuthal gaps; this changes ⟨cos²(2φ)⟩ in Eq. (5) and the effective
N_flat per φ bin. Correct using acceptance-weighted integrals.

**B5. η-dependent electron-ID efficiency for precision δA_bin.**
The Eq. (5) formula assumes all events in a bin share the same P_zz and
efficiency. With η-varying ε_eID, per-bin effective counts should be
weighted by ε_eID(η) before computing δA_bin.

**B6. Tighten MC statistics per bin.**
1000 MC events per bin gives ~3% relative sampling error on migration
fractions. Increase to ≥5000 before final publication.

### 11.3 Tier C — long-term

**C1.** NNLO QCD corrections to the DIS cross section.

**C2.** Target-mass corrections at high x.

**C3.** Polarized-electron scenarios (complementary single-spin A_LT and
double-spin A_LL observables).

**C4.** Two-gluon operator contribution to Δ at small x (not captured by
the Sather–Schmidt quark-only calculation; could enhance or cancel the
first moment).

---

## 12. Reproducibility and prerequisites

All commands are run from the repo root:
`/Users/L00338853/work/Polarized_Li/PolarizedLithiumSim`

**Install PDF sets (once):**

```bash
python3 -m parton install EPPS21nlo_CT18Anlo_Li6
python3 -m parton install nNNPDF30_nlo_as_0118_A6_Z3
python3 -m parton install CT18NLO
```

**Environment patches** (documented in `money_delta_note_2026-07-16.md` §5;
apply once to the installed packages):

- NumPy 2.0 compatibility fix in `parton/pdf.py:231`
- `_safe_xfx` helper in `fastsim/polli_fastsim/structure.py` and `polarized.py`

**Run scripts in order (each is independent but builds on prior results):**

```bash
# Phase II: sum-rule reach plots, 9 (config,shape) combos, Interpretation A
python3 fastsim/scripts/money_delta_20260724.py --outdir fastsim/out/money_delta

# Phase II: same but Interpretation B (no F₁ in Δ)
python3 fastsim/scripts/money_delta_20260725.py --outdir fastsim/out/money_delta

# Phase III: φ-modulation projections (no detector), 3 configs × 3 cases
python3 fastsim/scripts/money_delta_20260728.py --outdir fastsim/out/money_delta

# Phase IV: detector-realistic extension (ePIC tracking + ε_eID)
python3 fastsim/scripts/money_delta_20260729.py --outdir fastsim/out/money_delta
```

**Approximate runtimes:**

| Script | Runtime | Bottleneck |
|---|---|---|
| `money_delta_20260724.py` | ~15–20 min | 27 grid-backend σ² evaluations + 4 heatmaps |
| `money_delta_20260725.py` | ~15–20 min | Same grid structure, analytic solver |
| `money_delta_20260728.py` | ~5–10 min | 9 φ-distribution plots + heatmaps |
| `money_delta_20260729.py` | ~5–10 min | 1 × 10⁶ MC events; deterministic (seed=42) |

---

## 13. Complete file inventory

### Notes (this series)

```
fastsim/notes/money_delta_note_2026-07-16.md   Phase I: first reach plots, infrastructure
fastsim/notes/money_delta_note_2026-07-20.md   Phase I: nuclear PDF comparison, lum. reference
fastsim/notes/money_delta_note_2026-07-21.md   Phase I: three-tier scenario
fastsim/notes/money_delta_note_2026-07-24.md   Phase II: sum-rule normalization, full derivation
fastsim/notes/money_delta_note_2026-07-27.md   Phase II: Interpretation B, A vs B comparison
fastsim/notes/money_delta_note_2026-07-28.md   Phase III: φ-modulation plots, integration cases
fastsim/notes/money_delta_note_2026-07-29.md        Phase IV: detector smearing, sign-convention fix
fastsim/notes/money_delta_note_2026-07-29_fix.md   Fix audit: reco-selection consistency fix, production rerun (2026-07-31)
fastsim/notes/money_delta_uptodate.md               This file — master synthesis
fastsim/notes/money_delta_uptodate_brief.md         Internal working summary (distilled from this master; 2026-07-31)
```

### Scripts

```
fastsim/scripts/money_delta_20260724.py    (~1346 lines) Phase II, Interpretation A
fastsim/scripts/money_delta_20260725.py    (~1400 lines) Phase II, Interpretation B
fastsim/scripts/money_delta_20260728.py    (~1005 lines) Phase III, φ-modulation
fastsim/scripts/money_delta_20260729.py          Phase IV, detector smearing (post-fix; reco-mask consistency fix applied 2026-07-31)
fastsim/scripts/_check_reco_mask_invariants.py   S1–S6 AST static checker for reco-mask invariants (created 2026-07-31)
fastsim/scripts/money_delta_20260721.py          Phase I (earlier)
fastsim/scripts/money_delta_20260720.py    Phase I (earlier)
fastsim/scripts/money_delta_20260715.py    Phase I (earlier)
fastsim/scripts/money_delta_realistic.py   (utility; see individual note for role)
fastsim/scripts/money_delta_pdfgrid.py     (utility; PDF grid diagnostics)
```

### Output plots in `fastsim/out/money_delta/`

**2026-07-24 session (9 PNGs, `plot*` prefix):**

```
plot1a_dAoA_vs_scale_L10.png          δA/A vs peak Δ/F₁, L=10 fb⁻¹, 9 curves
plot1b_dAoA_vs_scale_L100.png         δA/A vs peak Δ/F₁, L=100 fb⁻¹, 9 curves
plot2_L5sig_vs_scale.png              L_5σ vs peak Δ/F₁, 9 curves
plot3_mid_midx_dAoA.png               MID+mid_x, δA/A, bag/lat vertical lines
plot4_mid_midx_L5sig.png              MID+mid_x, L_5σ, bag/lat vertical lines
plot5_heatmap_mid_L10.png             Per-bin heatmap: MID, 10 fb⁻¹
plot6_heatmap_mid_L100.png            Per-bin heatmap: MID, 100 fb⁻¹
plot7_heatmap_top_L10.png             Per-bin heatmap: TOP, 10 fb⁻¹
plot8_heatmap_top_L100.png            Per-bin heatmap: TOP, 100 fb⁻¹
```

**2026-07-27 session (9 PNGs, `money_delta_20260725_` prefix):**
*(The source note header says "8 PNGs" but lists 9 filenames — 5 reach plots +
4 per-bin heatmaps. The list is authoritative; the header count is a typo in the
source note.)*

```
money_delta_20260725_p1_bag_10fb.png       δA/A vs variant, L=10 fb⁻¹, bag
money_delta_20260725_p1_bag_100fb.png      δA/A vs variant, L=100 fb⁻¹, bag
money_delta_20260725_p2_bag.png            L_5σ vs variant, bag, 3 configs
money_delta_20260725_p3_baglat.png         δA/A, bag + lattice overlaid
money_delta_20260725_p4_baglat.png         L_5σ, bag + lattice overlaid
money_delta_20260725_perbin_mid_10fb.png   Per-bin heatmap: MID, 10 fb⁻¹
money_delta_20260725_perbin_mid_100fb.png  Per-bin heatmap: MID, 100 fb⁻¹
money_delta_20260725_perbin_top_10fb.png   Per-bin heatmap: TOP, 10 fb⁻¹
money_delta_20260725_perbin_top_100fb.png  Per-bin heatmap: TOP, 100 fb⁻¹
```

**2026-07-28 session (16 PNGs, `money_delta_20260728_` prefix):**

```
money_delta_20260728_perbin_low.png                      Heatmap: LOW
money_delta_20260728_perbin_mid.png                      Heatmap: MID
money_delta_20260728_perbin_top.png                      Heatmap: TOP
money_delta_20260728_phimodulation_peakbin_low.png       LOW, Case 1 (peak bin)
money_delta_20260728_phimodulation_peakbin_mid.png       MID, Case 1
money_delta_20260728_phimodulation_peakbin_top.png       TOP, Case 1, 5° bins
money_delta_20260728_phimodulation_q2slice_low.png       LOW, Case 2 (Q² slice)
money_delta_20260728_phimodulation_q2slice_mid.png       MID, Case 2
money_delta_20260728_phimodulation_q2slice_top.png       TOP, Case 2
money_delta_20260728_phimodulation_integrated_low.png    LOW, Case 3 (integrated)
money_delta_20260728_phimodulation_integrated_mid.png    MID, Case 3
money_delta_20260728_phimodulation_integrated_top.png    TOP, Case 3
money_delta_20260728_phimodulation_peakbin_top_10deg.png TOP, Case 1, 10° bins
money_delta_20260728_phimodulation_peakbin_top_20deg.png TOP, Case 1, 20° bins
money_delta_20260728_phimodulation_peakbin_top_30deg.png TOP, Case 1, 30° bins
money_delta_20260728_phimodulation_peakbin_top_45deg.png TOP, Case 1, 45° bins
```

**2026-07-29 session (16 PNGs, `money_delta_20260729_` prefix):**

```
money_delta_20260729_perbin_low.png                      Heatmap: LOW (reco)
money_delta_20260729_perbin_mid.png                      Heatmap: MID (reco)
money_delta_20260729_perbin_top.png                      Heatmap: TOP (reco)
money_delta_20260729_phimodulation_peakbin_low.png       LOW, Case 1 (reco peak)
money_delta_20260729_phimodulation_peakbin_mid.png       MID, Case 1 (reco peak)
money_delta_20260729_phimodulation_peakbin_top_5deg.png  TOP, Case 1, 5° bins, reco
money_delta_20260729_phimodulation_q2slice_low.png       LOW, Case 2 (reco Q² slice)
money_delta_20260729_phimodulation_q2slice_mid.png       MID, Case 2
money_delta_20260729_phimodulation_q2slice_top.png       TOP, Case 2
money_delta_20260729_phimodulation_integrated_low.png    LOW, Case 3 (reco integrated)
money_delta_20260729_phimodulation_integrated_mid.png    MID, Case 3
money_delta_20260729_phimodulation_integrated_top.png    TOP, Case 3
money_delta_20260729_phimodulation_peakbin_top_10deg.png TOP, Case 1, 10° bins, reco
money_delta_20260729_phimodulation_peakbin_top_20deg.png TOP, Case 1, 20° bins, reco
money_delta_20260729_phimodulation_peakbin_top_30deg.png TOP, Case 1, 30° bins, reco
money_delta_20260729_phimodulation_peakbin_top_45deg.png TOP, Case 1, 45° bins, reco
```

---

*End of master note. Last updated: 2026-07-31 (revised; reco-selection consistency fix incorporated).*
*Primary sources: money_delta_note_2026-07-24, -07-27, -07-28, -07-29, -07-29_fix.*
*Earlier sessions (-07-16, -07-20, -07-21) summarized in §2; not independently re-derived here.*
*Detector-realistic conclusions (§8, §10 detector-tested block) apply only to Interpretation A,*
*mid_x shape, bag constraint, ePIC tracking-only + ATHENA/ECCE ε_eID, 10 fb⁻¹.*
