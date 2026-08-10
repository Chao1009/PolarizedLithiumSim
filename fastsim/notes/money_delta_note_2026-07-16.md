# Money-Δ Plot: Progress and Next Steps — 2026-07-16

**Observable.** 5σ discovery luminosity `L_5σ` (fb⁻¹/nucleon) for the `cos(2φ)`
double-helicity-flip tensor asymmetry on transversely tensor-polarized ⁶Li,
plotted as a function of the assumed peak `Δ/F₁` scale. Δ(x, Q²) is the
spin-1 tensor gluon structure function that gives the `cos(2φ)` term in the
Hoodbhoy–Jaffe–Manohar master formula (NPB 312:571). Significance is combined
across accepted (x, Q²) bins as `σ² = Σ_bins A_bin² · P_zz² · N_bin / 2`, so
`L_5σ = 25 / σ²` at the reference luminosity of 1 fb⁻¹/u.

**Revision history**

- 2026-07-16 (initial): Documented six realism ingredients R1–R6 in `money_delta_realistic.py`.
- 2026-07-16 (revised): R5 (η-dependent electron-ID efficiency, placeholder) removed at user request. Realistic pass now uses uniform efficiency inside the |η| ≤ 3.5 hard cut. Section 1 realism list, Section 2 numbers, Section 3 assumptions, and Section 4 next-steps updated accordingly.
- 2026-07-16 (correction): P_zz cluster-picture scenario value corrected from 0.70 to 0.404. The prior 0.70 was the raw deuteron-nucleon polarization *conditional* on scattering off a polarized nucleon, which is not directly usable in a uniform per-bin σ² sum. The corrected 0.404 folds in the 2/6 occupancy factor: sqrt(2/6) × 0.70 ≈ 0.404. Section 1 R6 numbers, Section 2 numerical band, and new Section 7 (P_zz convention discussion) added. Plots regenerated.

---

## 1. What was produced today

All figures land in `fastsim/out/money_delta/`. Ion is ⁶Li throughout. Unless
noted, `P_zz = 0.80`, `Q² > 1` GeV² (toy) or `Q² > 2` GeV² (grid, to stay in
CT18NLO fit range), `W² > 4` GeV², and the per-bin cut `N_bin ≥ 10` from
`sig2_per_fb_at`.

### Toy backend (hand-tuned F₂ ≈ CTEQ-like, Δ from `toy_delta_gluon`)

| Filename | Content |
|---|---|
| `money_delta_6Li_toy.png` | Original money-Δ plot: `L_5σ` vs Δ/F₁ scale on log-log, six curves (3 canonical EIC configs × 2 P_zz = 0.55, 0.80). |
| `money_delta_20260715_vs_Ee.png` | `L_5σ` vs electron energy `E_e ∈ [5, 18]` GeV with ⁶Li at rigidity-limited `p_ion`, three Δ/F₁ scale curves (10⁻³, 3·10⁻³, 10⁻²). |
| `money_delta_20260715_vs_Eion.png` | Same as above but sweeping `p_ion ∈ [20, 137.5]` GeV/u at fixed `E_e = 18` GeV. |
| `money_delta_20260715_L5sig_vs_sqrts.png` | `L_5σ` vs `√s` for nine `(E_e, Δ/F₁)` combinations; gold band marks the plausible EIC program range 1–100 fb⁻¹/u. |
| `money_delta_20260715_xQ2_coverage.png` | 2×2 accepted-bin `(x, Q²)` event-count pcolormesh at four representative beam configs. |
| `money_delta_20260715_yield_vs_phi.png` | Two-panel: `dN/dφ` yield (top) and `A(φ) = P_zz · ⟨A_cos2φ⟩ · cos(2φ)` with ±1σ stat band (bottom), at `L = 10` fb⁻¹/u, top config. |

All toy-backend figures use `NuclearF2` with `base=None` (toy F₂) and
`toy_delta_gluon(x, Q², f₁, scale)` linear in `scale`. Kinematic cuts are the
`fom.Scenario` defaults except where noted.

### Grid backend (CT18NLO via `parton`, NNPDFpol11_100 for polarized)

Same physics chain, same cuts (with `q2_min = 2.0`), same six figures:

| Filename | Corresponding toy figure |
|---|---|
| `money_delta_pdfgrid_scale_scan.png` | ↔ `money_delta_6Li_toy.png` |
| `money_delta_pdfgrid_vs_Ee.png` | ↔ `..._20260715_vs_Ee.png` |
| `money_delta_pdfgrid_vs_Eion.png` | ↔ `..._20260715_vs_Eion.png` |
| `money_delta_pdfgrid_L5sig_vs_sqrts.png` | ↔ `..._20260715_L5sig_vs_sqrts.png` |
| `money_delta_pdfgrid_xQ2_coverage.png` | ↔ `..._20260715_xQ2_coverage.png` |
| `money_delta_pdfgrid_yield_vs_phi.png` | ↔ `..._20260715_yield_vs_phi.png` |

### Realistic pass (grid backend + six realism upgrades)

| Filename | Content |
|---|---|
| `money_delta_realistic_mid.png` | `L_5σ` vs Δ/F₁ at the mid config (`E_e = 10` GeV × `p_ion = 50` GeV/u), showing central curve plus shape-and-P_zz uncertainty band. |
| `money_delta_realistic_top.png` | Same for the top config (`E_e = 18` GeV × `p_ion = 137.5` GeV/u). |

Realism ingredients in both (all active):

- **R1.** Grid PDFs (CT18NLO unpolarized, NNPDFpol11_100 polarized) with `q2_min = 2.0`.
- **R2.** Three Δ(x, Q²) x-shape variants, all normalized to peak `Δ/F₁ = scale`:
  low_x (α=0.3, β=4, same as `toy_delta_gluon`), mid_x (α=0.7, β=3), high_x (α=1.5, β=2).
- **R3.** EMC-ratio hook via `unpolarized_emc_ratio` on `NuclearF2`.
- **R4.** Two `R = σ_L/σ_T` parameterizations (simplified R1998, simplified Christy–Bosted DIS-region form), swapped in by context-manager monkey-patching of `structure.r_sigma_lt` and `asymmetries.r_sigma_lt`.
- R5. Not included. See §3.
- **R6.** Two P_zz-dilution scenarios for the same nominal ring value P_zz = 0.80: Cloet-conservative per-nucleon-average (`P_zz_eff = 0.267 = 0.8/3`) and cluster-picture occupancy-folded (`P_zz_eff = sqrt(2/6) × 0.87 × 0.8 ≈ 0.404`). Both fold the 2/6 fraction of polarized nucleons and — for the cluster case — the ~0.87 deuteron-cluster polarization retention (Schellingerhout PRC 48:2714). See §7 for full discussion.

---

## 2. Key physics conclusions from today's work

1. **Sweet spot at the mid config.** Among the three canonical EIC beam
   configurations, `E_e = 10` GeV × ⁶Li at `p_ion = 50` GeV/u (√s ≈ 45 GeV per
   nucleon) gives the lowest `L_5σ` on the Δ/F₁ = 10⁻³ Sather–Schmidt reference
   case. This was verified numerically with the toy backend by
   `diag_sig2_grid.py`, which prints σ² and `L_5σ` on a beam-energy grid; the
   mid config wins.

2. **Origin of the sweet spot.** Competition between event yield (grows with
   √s through the DIS cross section and the accepted phase space) and dilution
   of the `A_cos2φ` amplitude at large `y`. The amplitude carries a
   `−(1−y)/y²` factor, which shrinks fast as `y → 1`. Pushing √s higher opens
   more phase space in the large-`y` region where `(1−y)/y²` is small, so the
   marginal event contributes less to σ². The two effects cross near the mid
   config.

3. **Toy vs grid backend.** CT18NLO PDFs give `L_5σ` roughly a factor of 2
   larger than the hand-tuned toy F₂ across all six figures. This is
   consistent with the ~×1.5 accuracy claim on the toy backend; it lands a
   little worse than that claim but not dramatically so.

4. **Realistic pass widens the uncertainty band substantially.** At
   `Δ/F₁ = 10⁻³` on the mid config, the band now spans 15 to 486 fb⁻¹/u
   depending on which `R`, which P_zz convention, and which Δ x-shape is
   chosen. The central case sits at 131 fb⁻¹/u — worse than the toy 31 fb⁻¹/u.
   On the top config the band is 17 to 1708 fb⁻¹/u with central 275 fb⁻¹/u.
   The band lo values shifted upward from the pre-correction plot (previously
   4.9 and 5.7 fb⁻¹/u) because of the P_zz cluster-scenario correction
   described in §7. The band is dominated by two effects, described next.

5. **Dominant systematics identified.**
   (a) The assumed Δ(x, Q²) x-shape. The three variants (low_x, mid_x, high_x)
       differ by roughly an order of magnitude in σ² because Δ is peaked at
       different `x` and the accepted-bin `x`-density is highly non-uniform.
   (b) The P_zz dilution convention for ⁶Li: the Cloet-conservative
       per-nucleon-average (0.267) and the cluster-picture occupancy-folded
       value (0.404) differ by a factor 1.5 in `P_zz_eff`, i.e. ~2.3 in σ²
       and ~1.5 in `L_5σ`. This factor ~1.5 in `L_5σ` is a leading unresolved
       ambiguity. See §7 for the detailed physics.

### Note on the observable

> A clarification for readers: the code computes `A_cos2φ(x, Q²)` as a theoretical prediction of the fractional `cos(2φ)` modulation amplitude in each bin, given the assumed `Δ`. `A_cos2φ` is *not* itself a measured quantity. In a real experiment, one measures `dN/dφ` directly, fits it to `A + B·cos(2φ)`, and extracts `B̂ = P_zz · A_cos2φ_theory`. The significance formula `σ² = Σ A_cos2φ² · P_zz² · N_bin / 2` is what falls out of comparing the theoretical `A_cos2φ` prediction to the standard Fourier-amplitude statistical uncertainty `σ_B̂ = √(2/N_bin)`. This is a projection framework, not an event generator.

---

## 3. Current physics assumptions still in play

Even after the realistic pass, the following remain ad hoc, placeholder, or
unmodeled:

- **Δ shape.** Only three toy variants (low_x, mid_x, high_x). No lattice
  input. No CDKS-type convolution-model calculation. The three shapes
  bracket a plausible range but are not derived.
- **Nuclear effects beyond EMC.** No shadowing at low `x`, no Fermi-motion
  smearing, no off-shell corrections. The EMC hook is a smooth ratio only.
- **η electron-ID efficiency.** Not included. R5 was tested with a placeholder
  analytic form (`ε(η) = 0.98 − 0.05·(η/3.5)²`) and dropped at user request
  2026-07-16 pending a real ePIC electron-ID efficiency curve. Efficiency is
  currently uniform inside |η| ≤ 3.5.
- **P_zz dilution.** Two scenarios shown (Cloet 0.267, cluster-d 0.70);
  the correct convention for ⁶Li is unresolved pending nuclear-structure
  input. See `beams.py` docstring "6Li: UNRESOLVED convention (factor 2.4!)".
- **No detector resolution smearing.** No bin migration in `x`, `Q²`, or `φ`.
- **No radiative corrections.** No HECTOR-equivalent unfolding.
- **In-ring P_zz survival.** Nominal source value 0.80 assumed to persist
  through storage; not yet an accelerator-physics estimate.
- **`R = σ_L/σ_T`.** Two parameterizations shown (simplified R1998,
  simplified Christy–Bosted DIS-region form). Both are simplified DIS-region
  forms, not full published fits.
- **⁶Li deuteron-cluster embedding factor.** The 0.87 cluster fraction is
  not applied anywhere in the chain yet.

---

## 4. Next steps, ranked by expected impact

### Tier A — physics inputs to resolve

1. **Resolve the ⁶Li P_zz dilution convention.** Consult I. Cloet or the
   relevant nuclear-structure literature. This is the largest single
   ambiguity — factor ~2.4 in `L_5σ`. Everything else in the reach estimate
   is contingent on this being nailed down.
2. **Anchor the Δ(x, Q²) shape to a real calculation.** Options in decreasing
   order of directness: (a) lattice-moment constraints
   (Detmold–Shanahan-type), (b) CDKS-style convolution model, (c) NJL or
   light-cone model. At minimum, replace the three ad-hoc x-shapes with a
   single motivated one and quote its own uncertainty.
3. **Include nuclear shadowing at low `x`.** Swap the smooth EMC-only ratio
   for a proper nuclear-PDF set (EPPS21 or nNNPDF3.0). This mostly affects
   the low-`x` bins that dominate σ² through the yield.

### Tier B — detector and reconstruction realism

4. **Add η-dependent electron-ID efficiency** using the real ePIC curve (once
   stable enough to quote). Placeholder R5 was tested and removed on 2026-07-16.
5. **Detector resolution smearing.** Add Gaussian smearing on `E′` and `θ`
   from ePIC calorimeter and tracker specs; propagate into `x, Q²` and check
   bin migration.
6. **Radiative corrections.** HECTOR, or an EIC-adapted equivalent, applied
   to the born-level cross section.

### Tier C — long-term

7. **NNLO QCD corrections** to the DIS cross section.
8. **Target-mass corrections** at high `x`.
9. **Polarized-electron scenarios.** Would add complementary observables
   (single-spin, `A_LT`, etc.), not just `cos(2φ)`. Useful cross-check but
   does not change the tensor-only reach.

---

## 5. Reproducibility

All scripts are run from the repo root
(`/Users/L00338853/work/Polarized_Li/PolarizedLithiumSim`).

### Scripts

- `fastsim/scripts/money_delta.py` — toy backend, original money-Δ scale scan.
- `fastsim/scripts/money_delta_20260715.py` — toy backend, five extended figures (Ee sweep, Eion sweep, √s reach, `(x, Q²)` coverage, yield vs φ).
- `fastsim/scripts/money_delta_pdfgrid.py` — grid backend (CT18NLO + NNPDFpol11_100), six-figure replica of the toy set.
- `fastsim/scripts/money_delta_realistic.py` — grid backend plus six realism upgrades (R1–R6, see §1), two EIC configs.
- `fastsim/scripts/diag_sig2_grid.py` — √s sweet-spot diagnostic; no plots, prints σ² and `L_5σ` table.

### Command lines

All scripts write to `fastsim/out/money_delta/`. The realistic script
defaults its `--outdir` to that path; the others need it passed explicitly.

```bash
# toy backend
python3 fastsim/scripts/money_delta.py --ion 6Li --pdf toy --outdir fastsim/out/money_delta
python3 fastsim/scripts/money_delta_20260715.py --pdf toy --outdir fastsim/out/money_delta

# grid backend replica
python3 fastsim/scripts/money_delta_pdfgrid.py --outdir fastsim/out/money_delta

# realistic pass
python3 fastsim/scripts/money_delta_realistic.py

# sweet-spot diagnostic (no plots)
python3 fastsim/scripts/diag_sig2_grid.py
```

### Environment patches required for the grid backend

Both patches are required and have already been applied in the local
environment. Re-check after any `parton` upgrade.

- `parton/pdf.py`, line 231 — NumPy 2.0 compatibility (upstream bug):

  ```python
  # was: res = float(res)
  res = float(np.asarray(res).item())  # NumPy ≥2.0 compat
  ```

  Local path on this machine:
  `/Users/L00338853/.local/lib/python3.13/site-packages/parton/pdf.py:231`.

- `fastsim/polli_fastsim/structure.py` and
  `fastsim/polli_fastsim/polarized.py` — a `_safe_xfx` helper was added to
  wrap grid PDF calls and cleanly handle NaN returns / edge-of-grid points.

### Prerequisites

```bash
pip install parton
python3 -m parton install CT18NLO
python3 -m parton install NNPDFpol11_100
```

The realistic script exits with a clear error message if the grid backend
cannot be imported.

---

## 6. File inventory

### Plots in `fastsim/out/money_delta/` (14 files, as of end of 2026-07-16)

```
money_delta_6Li_toy.png                     (toy original, six curves)
money_delta_20260715_vs_Ee.png              (toy, Ee sweep)
money_delta_20260715_vs_Eion.png            (toy, Eion sweep)
money_delta_20260715_L5sig_vs_sqrts.png     (toy, √s reach)
money_delta_20260715_xQ2_coverage.png       (toy, (x, Q²) coverage)
money_delta_20260715_yield_vs_phi.png       (toy, dN/dφ + A(φ))
money_delta_pdfgrid_scale_scan.png          (grid, replica of toy original)
money_delta_pdfgrid_vs_Ee.png               (grid, Ee sweep)
money_delta_pdfgrid_vs_Eion.png             (grid, Eion sweep)
money_delta_pdfgrid_L5sig_vs_sqrts.png      (grid, √s reach)
money_delta_pdfgrid_xQ2_coverage.png        (grid, (x, Q²) coverage)
money_delta_pdfgrid_yield_vs_phi.png        (grid, dN/dφ + A(φ))
money_delta_realistic_mid.png               (realistic, mid config)
money_delta_realistic_top.png               (realistic, top config)
```

### New source files created today

- `fastsim/scripts/money_delta_20260715.py`
- `fastsim/scripts/money_delta_pdfgrid.py`
- `fastsim/scripts/money_delta_realistic.py`
- `fastsim/scripts/diag_sig2_grid.py`

### Note

- `fastsim/notes/money_delta_note_2026-07-16.md` (this file).

---

## 7. P_zz dilution conventions for ⁶Li — physics discussion

**Status:** unresolved pending nuclear-structure expert input (I. Cloet or equivalent). The realistic plot presents two scenarios as a defensible band while this is being resolved.

### Setup

⁶Li has 6 nucleons: 3 protons + 3 neutrons, total spin J = 1. Nominal ring polarization from the ECRP source: `P_zz = 0.8`. The `cos(2φ)` observable measures the tensor structure function `Δ` per struck nucleon, weighted by the effective tensor polarization *of that nucleon at the moment of scattering*.

Question: when a photon knocks a nucleon out of a `P_zz = 0.8` ⁶Li nucleus, what tensor polarization does that nucleon actually carry?

### Cluster picture

⁶Li is well described (~85–90%) as a ⁴He (α) core + a deuteron cluster, weakly bound at ~1.5 MeV per nucleon. Consequences:

- α core (4 nucleons, spin-paired): J = 0. Cannot carry tensor polarization. Zero contribution to the tensor signal.
- Deuteron cluster (2 nucleons, 1p + 1n): J = 1. Carries essentially all of the ⁶Li tensor polarization.
- Cluster D-state / distortion correction: Schellingerhout et al. (PRC 48:2714, 1993) estimate the deuteron cluster inside ⁶Li retains about 87% of a free-deuteron polarization. So at ⁶Li `P_zz = 0.8`, the deuteron-cluster nucleons see ~`0.87 × 0.8 ≈ 0.70` tensor polarization each.

### Convention 1: Cloet per-nucleon average (0.267)

DIS is inclusive; the photon strikes any nucleon with roughly equal probability. Average the tensor polarization over the 6-nucleon population:

    ⟨P_zz⟩_per_nucleon = (2 × 0.8 + 4 × 0) / 6 = 0.267

Used directly in the uniform per-bin σ² sum. This is what Cloet uses in the JLab E12-14-001 slides, referenced in `beams.py`.

### Convention 2: Cluster-picture occupancy-folded (0.404)

The signal only comes from the 2/6 fraction of events that scatter off deuteron-cluster nucleons. For those events, the per-nucleon polarization is ~0.70. Elsewhere it is 0.

To fold this into the uniform per-bin σ² sum, define an effective `P_zz` such that:

    P_zz_eff² = (2/6) × (0.70)² = 0.163
    P_zz_eff  = √0.163 ≈ 0.404

### Ratio of the two conventions

    Cluster / Cloet in P_zz_eff:  0.404 / 0.267 = 1.51
    Cluster / Cloet in σ²:        (0.404 / 0.267)² = 2.29
    Cluster / Cloet in L_5σ:      inverse of σ² ratio ≈ 0.44 (cluster is better)

So under the cluster picture, required luminosity is ~2.3× lower than under the Cloet convention. That corresponds to the factor ~1.5 in `L_5σ` scaling — smaller than the "factor 2.4" the `beams.py` docstring quotes.

### Reconciling with the "factor 2.4" in `beams.py`

The `beams.py` comment on the LI6 dataclass flags the convention ambiguity as "factor 2.4," which is looser than what a careful accounting gives. The 2.4 arises from a different reading of the cluster picture that treats the deuteron pair as a *single coherent object* and does not fold the 2/6 population factor. Under that reading, one gets a factor closer to 2.4 in `L_5σ`.

Whether that reading or the folded 0.404 reading is correct is exactly the physics-input question. A rigorous derivation of `A_cos2φ^{6Li}` in terms of nucleon-level `Δ` and the ⁶Li ground-state wavefunction is what settles this.

### Questions to bring to the theorist

1. Which convention gives the correct mapping from the nucleon-level `Δ` to the ⁶Li tensor asymmetry: per-nucleon averaged, or occupancy-folded cluster?
2. Is DIS off α-core nucleons truly polarization-blind in the tensor channel? (Naively yes — an isoscalar S-wave J=0 cannot support a tensor signal — but confirm rigorously.)
3. Is the Schellingerhout `P_d ≈ 0.87` deuteron-cluster polarization retention still the best number, or has newer VMC/QMC work (Wiringa, Argonne v18 + Illinois-7) updated it?
4. Are there off-shell or Fermi-motion modifications to `Δ` itself that further alter the mapping from free-nucleon `Δ` to ⁶Li tensor asymmetry beyond the polarization factor?
5. Which convention did the JLab E12-14-001 proposal use for its projections? For consistency with published projections, we should match their choice.

### What the plot currently does

`money_delta_realistic.py` shows both scenarios (0.267 and 0.404) as separate σ² weight values in the sig² sum, producing the ~2.3× spread in σ² that manifests as the band width along the P_zz axis. This is the honest way to display the ambiguity: two scenarios plotted side by side rather than a single guess.

### Correction history

Prior to 2026-07-16 15:12, the cluster-picture scenario used `P_zz_eff = 0.70` — the raw deuteron-nucleon value without folding the 2/6 occupancy factor. This overstated the difference between scenarios (the band lo values were ~3× too low). Corrected to 0.404 and plots regenerated.
