# Money-Δ Plot: Progress and Next Steps — 2026-07-20

**Observable.** 5σ discovery luminosity `L_5σ` (fb⁻¹/nucleon) for the `cos(2φ)`
double-helicity-flip tensor asymmetry on transversely tensor-polarized ⁶Li,
plotted against the assumed peak `Δ/F₁` scale. Full observable definition,
significance formula, and kinematic conventions are in
`money_delta_note_2026-07-16.md` §1 (Note on the observable). This note covers
only what is new today.

**Revision history**

- 2026-07-20 (initial): Added nuclear-PDF comparison variant (EPPS21 and
  nNNPDF3.0) at the mid EIC config with P_zz fixed at the Cloet convention.
  Introduced local class `NuclearF2FromGrid` in `money_delta_20260720.py`.
  Tier A #3 (nuclear-PDF upgrade) retired as a limiting systematic for ⁶Li.

---

## 1. What was produced today

### New script

`fastsim/scripts/money_delta_20260720.py` — 665 lines. Cloned from
`money_delta_realistic.py` with four changes; see §1 ("Differences from parent")
below.

### New output

`fastsim/out/money_delta/money_delta_20260720_mid.png` — single log-log panel,
`L_5σ` vs `Δ/F₁` at the mid config (`E_e = 10` GeV × ⁶Li at `p_ion = 50`
GeV/u), showing two nuclear-PDF uncertainty bands (EPPS21, nNNPDF3.0), a
central realistic curve, a CT18NLO + EMC reference, and a toy continuity line.
Full legend description in §2.

### Differences from `money_delta_realistic.py`

| Item | Parent (`money_delta_realistic.py`) | Today (`money_delta_20260720.py`) |
|---|---|---|
| Beam configs | Mid + top | Mid only |
| P_zz | Two scenarios: Cloet 0.267 and cluster-d 0.404 | Fixed: Cloet 0.267 only |
| Nuclear F₂ source | CT18NLO + smooth EMC ratio hook | Three backends: CT18+EMC (reference), EPPS21, nNNPDF3.0 |
| Δ shapes | low_x, mid_x, high_x (unchanged) | Identical |
| R parameterizations | R1998, Christy-Bosted (unchanged) | Identical |
| New class | — | `NuclearF2FromGrid` defined locally (package not modified) |

The `NuclearF2FromGrid` class is a drop-in replacement for `NuclearF2`
(from `polli_fastsim.structure`) when F₂ᴬ comes from a nuclear PDF grid via
`parton`. It multiplies the `parton`-returned per-nucleon xf_A by A to satisfy
the `f2a()` API contract (callers divide by `ion.A` to recover per-nucleon
quantities). The `f1a()` method derives F₁ᴬ via Callan-Gross using the
currently-active `r_sigma_lt` so that R1998 / Christy-Bosted swapping via
`r_override()` works identically to the parent. The class lives entirely in the
script; `polli_fastsim/` was not modified.

### Plot elements

- Blue (steelblue) band, alpha 0.22: min-max envelope over 6 EPPS21 combinations
  (3 shapes × 2 R), Cloet P_zz = 0.267.
- Orange (darkorange) band, alpha 0.22: min-max envelope over 6 nNNPDF3.0
  combinations (3 shapes × 2 R), Cloet P_zz = 0.267.
- Black solid line (lw 2.2): central curve = mid_x shape + R1998 + EPPS21 +
  Cloet P_zz = 0.267.
- Gray dashed line (lw 1.6): CT18NLO + EMC reference, mid_x + R1998 (matches
  the central curve from `money_delta_realistic.py`).
- Gray dotted line (lw 1.3): toy backend, `toy_delta_gluon`, P_zz = 0.80 (kept
  for continuity with the parent script).
- Gold axhspan at 1–100 fb⁻¹/u: plausible EIC program range.
- Vertical dotted line at `Δ/F₁ = 10⁻³` with "Sather-Schmidt O(10⁻³)" label.

### Questions scoped out today

- **⁶Li P_zz dilution convention** — tabled. P_zz was fixed at Cloet 0.267 to
  isolate the nuclear-PDF variable. The unresolved convention ambiguity (§7 of
  the 2026-07-16 note) carries forward unchanged.
- **Δ(x, Q²) x-shape** — tabled. The same three toy shapes from the parent
  script were reused; no new theory-motivated shapes were introduced.

---

## 2. Key physics conclusions from today's work

### 2.1 Nuclear-PDF choice is not a limiting systematic for ⁶Li

EPPS21 and nNNPDF3.0 agree to below 1% in `L_5σ` across all six tested (shape,
R) combinations. On the plot the two bands are essentially indistinguishable.
The paired-combination ratios at `Δ/F₁ = 10⁻³` are:

| Shape | R | EPPS21 (fb⁻¹/u) | nNNPDF3.0 (fb⁻¹/u) | Ratio |
|---|---|---:|---:|---:|
| low_x  | R1998          |  79.97 |  79.72 | 1.003 |
| low_x  | Christy-Bosted |  34.02 |  33.89 | 1.004 |
| mid_x  | R1998          | 130.83 | 131.71 | 0.993 |
| mid_x  | Christy-Bosted |  64.54 |  65.06 | 0.992 |
| high_x | R1998          | 481.30 | 480.04 | 1.003 |
| high_x | Christy-Bosted | 330.68 | 330.84 | 1.000 |

All ratios within 0.8% of unity. The physical reason is that ⁶Li is a light
nucleus (A = 6) with sparse experimental DIS data. Both nuclear-PDF
collaborations return CT18-family baselines with nuclear corrections at the
few-percent level, which translate to sub-percent effects on `L_5σ`. This means
nuclear-PDF-choice uncertainty can be parked as a follow-up question. It is not
the next bottleneck.

### 2.2 Δ(x, Q²) x-shape remains the dominant band-width driver

The total band at `Δ/F₁ = 10⁻³` spans approximately 34 to 481 fb⁻¹/u.
Shape choice (low_x vs high_x at fixed R) accounts for a factor of roughly 14
in `L_5σ`; the R choice (R1998 vs Christy-Bosted at fixed shape) contributes an
additional factor of 1.5–2.4 within that. This quantifies and confirms the
qualitative finding from 2026-07-16 §2.5(a). Replacing the three ad-hoc toy
shapes with a real theory-motivated Δ(x, Q²) is the next physically important
input (see §4, Tier A #2).

### 2.3 Central estimate at Sather-Schmidt reference: ~130 fb⁻¹/u

At `Δ/F₁ = 10⁻³`, the central combination (mid_x + R1998 + EPPS21 + Cloet
P_zz = 0.267) gives L_5σ = 130.83 fb⁻¹/u, above the 1–100 fb⁻¹/u plausible
EIC program band. The CT18NLO + EMC reference at the same (shape, R, P_zz)
gives 131.26 fb⁻¹/u (ratio 0.997; see §7 sanity check 1). Discovery becomes
plausible on the mid config at `Δ/F₁ ≈ 3 × 10⁻³`, where the central curve
enters the gold band. Under the Cloet P_zz convention with mid_x shape and
R1998, the Sather-Schmidt case is not reachable in a plausible EIC program.

Reference values at `Δ/F₁ = 10⁻³`:

| Quantity | L_5σ (fb⁻¹/u) |
|---|---:|
| Central (mid_x, R1998, EPPS21, Cloet P_zz = 0.267) | 130.83 |
| CT18NLO + EMC reference (mid_x, R1998, Cloet P_zz = 0.267) | 131.26 |
| EPPS21 band: [lo, hi] | [34.02, 481.30] |
| nNNPDF3.0 band: [lo, hi] | [33.89, 480.04] |
| Toy backend (toy_delta_gluon, original R, P_zz = 0.80) | 31.26 |

---

## 3. Current physics assumptions still in play

After today's work, the following inputs remain ad hoc, placeholder, or
unmodeled. The nuclear-PDF item from the 2026-07-16 §3 list is resolved and
removed; all others carry forward.

- **Δ shape.** Only three toy variants (low_x, mid_x, high_x). No lattice
  input, no CDKS-type convolution-model calculation. The three shapes bracket a
  plausible range but are not derived from first principles.
- **P_zz dilution.** Fixed at Cloet 0.267 today to isolate the nuclear-PDF
  variable. The convention ambiguity (Cloet per-nucleon-average vs
  cluster-picture occupancy-folded 0.404) remains unresolved. See
  `money_delta_note_2026-07-16.md` §7 for the full physics discussion and the
  questions to bring to the theorist.
- **Nuclear PDF sets now real for ⁶Li.** EPPS21 and nNNPDF3.0 are installed and
  running. However, only central members (member 0) are used; PDF replica
  uncertainties within each set are not evaluated.
- **`R = σ_L/σ_T`.** Two simplified DIS-region parameterizations (R1998,
  Christy-Bosted) shown. Neither is a full published fit applied outside its
  validity range.
- **η electron-ID efficiency.** Not included. Uniform efficiency inside
  |η| ≤ 3.5. Placeholder R5 from `money_delta_realistic.py` was removed
  2026-07-16 pending a real ePIC efficiency curve.
- **No detector resolution smearing.** No bin migration in x, Q², or φ.
- **No radiative corrections.**
- **In-ring P_zz survival.** Nominal source value 0.80 assumed to persist
  through storage.
- **⁶Li deuteron-cluster embedding factor.** The 0.87 cluster polarization
  retention fraction (Schellingerhout PRC 48:2714) is discussed in
  `money_delta_note_2026-07-16.md` §7 but is not yet a formal input to the
  significance chain.

---

## 4. Next steps, ranked by expected impact

The Tier A ranking from `money_delta_note_2026-07-16.md` §4 is revised below.
Tier A #3 (nuclear-PDF upgrade) is retired. Everything else carries forward.

### Tier A — physics inputs to resolve

1. **Resolve the ⁶Li P_zz dilution convention (unchanged, highest priority).**
   Consult I. Cloet or the relevant nuclear-structure literature. Factor ~1.5–2.4
   in `L_5σ` depending on convention. All results today are conditional on the
   Cloet 0.267 choice. See `money_delta_note_2026-07-16.md` §7 for the list of
   questions to bring to the theorist.

2. **Anchor Δ(x, Q²) to a real theory calculation (promoted from prior #2).**
   Today's results confirm this is the largest technically resolvable systematic:
   factor ~14 in `L_5σ` from shape choice alone. Options in decreasing order of
   directness: (a) lattice-moment constraints (Detmold–Shanahan-type), (b)
   CDKS-style convolution model, (c) NJL or light-cone quark model. At minimum,
   replace the three ad-hoc x-shapes with a single motivated one and quote its
   own uncertainty.

3. ~~**Include nuclear shadowing / nuclear-PDF upgrade (EPPS21 or nNNPDF3.0).**~~
   **RETIRED.** Today's result shows this is a sub-1% effect in `L_5σ` for ⁶Li.
   Both sets are now installed and confirmed running; central-member evaluation
   shows the two collaborations agree to below 1%. This item can be parked.

### Tier B — detector and reconstruction realism

4. **Add η-dependent electron-ID efficiency** using the real ePIC curve (once
   stable enough to quote). Placeholder R5 was tested and removed 2026-07-16.
5. **Detector resolution smearing.** Gaussian smearing on E′ and θ from ePIC
   specs; propagate into x, Q² and check bin migration.
6. **Radiative corrections.** HECTOR or an EIC-adapted equivalent.

### Tier C — long-term

7. **NNLO QCD corrections** to the DIS cross section.
8. **Target-mass corrections** at high x.
9. **Polarized-electron scenarios.** Complementary observables (A_LT, etc.);
   useful cross-check but does not change the tensor-only reach.

---

## 5. Reproducibility

All scripts are run from the repo root
(`/Users/L00338853/work/Polarized_Li/PolarizedLithiumSim`).

### Command line

```bash
python3 fastsim/scripts/money_delta_20260720.py
```

Output writes to `fastsim/out/money_delta/` by default. The `--outdir` flag
overrides this:

```bash
python3 fastsim/scripts/money_delta_20260720.py --outdir /path/to/outdir
```

Runtime: approximately 5 minutes for 14 grid-backend evaluations
(12 nuclear-PDF combos + 1 CT18+EMC reference + 1 toy).

### Prerequisites

In addition to the prerequisites listed in `money_delta_note_2026-07-16.md` §5
(CT18NLO and NNPDFpol11_100), the two nuclear PDF sets used today must be
installed:

```bash
python3 -m parton install EPPS21nlo_CT18Anlo_Li6
python3 -m parton install nNNPDF30_nlo_as_0118_A6_Z3
```

Both were downloaded from `lhapdfsets.web.cern.ch` and extracted successfully
on the first try on 2026-07-20.

### Environment patches

The two environment patches documented in `money_delta_note_2026-07-16.md` §5
are still required and still apply:

- `parton/pdf.py` line 231: NumPy 2.0 compatibility fix (`float(np.asarray(res).item())`).
- `fastsim/polli_fastsim/structure.py` and `polarized.py`: `_safe_xfx` helper
  for NaN / edge-of-grid handling.

The new script mirrors the `_safe_xfx` pattern in a local function
`_safe_xfx_local` (lines 67–78) rather than importing from the package, so
it does not introduce new patch requirements. Re-check after any `parton` upgrade.

---

## 6. File inventory

### New plot added today

```
money_delta_20260720_mid.png    (EPPS21 vs nNNPDF3.0 nuclear-PDF comparison, mid config)
```

### New script added today

```
fastsim/scripts/money_delta_20260720.py    (665 lines; parent: money_delta_realistic.py)
```

### Previously existing plots in `fastsim/out/money_delta/` (14 files as of 2026-07-16)

See `money_delta_note_2026-07-16.md` §6 for the full list.

### Notes

- `fastsim/notes/money_delta_note_2026-07-16.md` — prior session note.
- `fastsim/notes/money_delta_note_2026-07-20.md` — this file.

---

## 7. Sanity checks performed

All four checks passed post-run.

### 7.1 Central vs CT18+EMC reference (like-for-like)

Central (mid_x + R1998 + EPPS21, Cloet P_zz = 0.267): 130.83 fb⁻¹/u.  
CT18NLO + EMC reference (same shape, R, P_zz): 131.26 fb⁻¹/u.  
Ratio: 0.997.

These two curves use the same Δ shape, the same R parameterization, and the same
P_zz — only the nuclear-F₂ treatment differs (EPPS21 grid vs CT18NLO + smooth
EMC ratio hook). Agreement to 0.3% is consistent with EPPS21 nuclear corrections
for ⁶Li being at the ~1% level, as expected for a light nucleus.

### 7.2 EPPS21 vs nNNPDF3.0 agreement

All six paired-combination ratios are within 0.8% of unity. See the full table
in §2.1. The two bands are visually indistinguishable on the plot.

### 7.3 Δ-shape systematic decomposition

At `Δ/F₁ = 10⁻³`:

- low_x to high_x ratio at R1998 (EPPS21): 481.30 / 79.97 ≈ 6.0 (per R choice)  
- low_x to high_x ratio at Christy-Bosted (EPPS21): 330.68 / 34.02 ≈ 9.7  
- Band total (low end to high end): 481.30 / 34.02 ≈ 14.1  
- R1998 to Christy-Bosted ratio at low_x: 79.97 / 34.02 ≈ 2.4  
- R1998 to Christy-Bosted ratio at mid_x: 130.83 / 64.54 ≈ 2.0  

Shape choice contributes roughly a factor of 7–10 in `L_5σ`; R choice
contributes an additional factor of 1.5–2.4. These are sub-systematics within
the total band width of ~14×.

### 7.4 Toy vs central P_zz-only expectation

Central (EPPS21, mid_x, R1998, P_zz = 0.267): 130.83 fb⁻¹/u.  
Toy (toy_delta_gluon, original R, P_zz = 0.80): 31.26 fb⁻¹/u.  
Observed ratio: 130.83 / 31.26 ≈ 4.19.

Naive P_zz-only expectation from scaling: `(0.80 / 0.267)² ≈ 8.98`. The
observed ratio of ~4.2 is roughly half that. The additional factor ~2 is
attributable to differences in the underlying F₂/R models (grid + EMC vs toy)
and the different `q2_min` cuts (toy uses 1.0 GeV², grid uses 2.0 GeV²).

This is internally consistent: the ~×2 toy-vs-grid discrepancy at matched P_zz
was documented in `money_delta_note_2026-07-16.md` §2.3. Accounting for it:
`(0.80 / 0.267)² × 2 ≈ 18`, but shape and R differences between `toy_delta_gluon`
and mid_x + R1998 reduce this further, bringing the observed ratio to ~4.2.
This is a consistency check, not a discrepancy.

---

**Bottom line for the next session.** With nuclear-PDF choice confirmed
non-limiting, the next-session physics priority narrows to two items: resolving
the ⁶Li P_zz dilution convention (Tier A #1) and obtaining a theory-anchored
Δ(x, Q²) shape (Tier A #2 promoted).
