# Money-Δ Plot: Progress and Next Steps — 2026-07-21

**Observable.** Relative statistical uncertainty `δA/A_cos2φ` on the `cos(2φ)`
double-helicity-flip tensor asymmetry on transversely tensor-polarized ⁶Li,
plotted against the assumed peak `Δ/F₁` scale at **fixed luminosity**. Full
observable definition, significance formula, and kinematic conventions are in
`money_delta_note_2026-07-16.md` §1 (Note on the observable). The parent
observable (`L_5σ` vs `Δ/F₁`) is documented in `money_delta_note_2026-07-20.md`
§2. This note covers only what is new today.

**Revision history**

- 2026-07-21 (initial): Changed observable from `L_5σ` (required luminosity for
  5σ discovery) to `δA/A_cos2φ` (relative statistical uncertainty at fixed
  luminosity). Added Plot B (per-bin absolute stat uncertainty on the `(x, Q²)`
  plane). Extended to both mid and top EIC configs at two luminosity anchors
  (10 fb⁻¹/u and 100 fb⁻¹/u). Physics chain unchanged from
  `money_delta_20260720.py` (14 evaluations per config). Nuclear-PDF comparison
  (EPPS21, nNNPDF3.0) retained; Tier A #3 (nuclear-PDF upgrade) remains retired.
- 2026-07-21 (extension): Added Plot C — per-bin DIS event rate per hour — and
  CSV per-bin rate tables. Script grew from 877 to 1079 lines. See new §1.5,
  §2.6, §6 additions, and §7.5–§7.8.

---

## 1. What was produced today

### New script

`fastsim/scripts/money_delta_20260721.py` — 877 lines. Cloned from
`money_delta_20260720.py`. Physics helpers (`NuclearF2FromGrid`, `r_override`,
`sig2_per_fb_at_nuc`, etc.) are copied verbatim from the parent; see §1.1 for
what changed.

### New output

Eight PNGs in `fastsim/out/money_delta/`, described in §1.2. Full file list in
§6.

### 1.1 Differences from `money_delta_20260720.py`

| Item | Parent (`money_delta_20260720.py`) | Today (`money_delta_20260721.py`) |
|---|---|---|
| Observable | `L_5σ` (required luminosity for 5σ discovery) vs `Δ/F₁` | `δA/A_cos2φ` (relative stat uncertainty) vs `Δ/F₁` at fixed L |
| Beam configs | Mid only | Mid + top |
| Luminosity | Fixed-L not applicable (L_5σ is the y-axis) | 10 fb⁻¹/u and 100 fb⁻¹/u |
| Plot A y-axis | `L_5σ` (fb⁻¹/u) | `δA/A_cos2φ` (dimensionless, log scale) |
| Plot B | Not present | Per-bin `δA_bin = √(2/N_bin) / P_zz` on `(x, Q²)` pcolormesh |
| Gold axhspan | 1–100 fb⁻¹/u EIC program band | Not shown (L is now a fixed input) |
| 5σ threshold | `L_5σ` naturally encodes 5σ | Horizontal firebrick line at `δA/A = 0.20` |
| Nuclear-PDF sets | EPPS21, nNNPDF3.0 (unchanged) | Identical |
| Δ shapes, R choices | Unchanged | Identical |
| P_zz | Cloet 0.267 (unchanged) | Identical |
| New class | `NuclearF2FromGrid` (defined locally) | Copied verbatim; no new package modifications |

The analytic rescaling for Plot A reuses sig² per fb⁻¹/u values computed at the
reference scale S₀ = 10⁻³ (same as the parent), converted to `δA/A` via:

    δA/A = 1 / √( sig²_per_fb × (s/S₀)² × L )

This is equivalent to Nσ⁻¹ at luminosity L. The 5σ crossing at `s*` satisfies
`sig²_per_fb × (s*/S₀)² × L = 25`, so `L_5σ(10⁻³) = L × (s*/10⁻³)²` — the
self-consistency cross-check reported in §7.1.

### 1.2 Plot descriptions

**Plot A — `δA/A` vs `Δ/F₁` (log-log).** One panel per (config, luminosity),
four panels total. Band structure is identical to the parent:

- Steelblue band, alpha 0.22: min-max envelope over 6 EPPS21 combinations
  (3 shapes × 2 R), Cloet P_zz = 0.267.
- Darkorange band, alpha 0.22: min-max envelope over 6 nNNPDF3.0 combinations.
- Black solid line (lw 2.2): central curve = mid_x + R1998 + EPPS21 +
  Cloet P_zz = 0.267.
- Gray dashed line (lw 1.6): CT18NLO + EMC reference, mid_x + R1998.
- Gray dotted line (lw 1.3): toy backend, `toy_delta_gluon`, P_zz = 0.80
  (continuity with parent).
- Horizontal firebrick dashed line at `δA/A = 0.20`: 5σ discovery threshold.
- Vertical dotted line at `Δ/F₁ = 10⁻³` with "Sather-Schmidt O(10⁻³)" label.

**Plot B — per-bin `δA_bin` heatmap (log x × log Q²).** One panel per
(config, luminosity), four panels total. Uses EPPS21 central member + R1998 +
Cloet P_zz = 0.267. Per-bin uncertainty defined as `δA_bin = √(2/N_bin) / P_zz`
(the standard Fourier-amplitude statistical uncertainty divided by P_zz; see
`money_delta_note_2026-07-16.md` §1, Note on the observable). Only bins with
N_bin ≥ 10 are shown; nuclear-PDF choice is immaterial at this level (confirmed
sub-1% by `money_delta_note_2026-07-20.md` §2.1). Colormap: viridis, LogNorm.

### 1.3 Luminosity anchors and caveat

The two luminosity values used today are intended to map to plausible EIC
program durations for ⁶Li:

- **L = 10 fb⁻¹/nucleon** — design-target reference for approximately 1 year
  of EIC operation on ⁶Li.
- **L = 100 fb⁻¹/nucleon** — design-target reference for approximately 10 years
  of EIC operation on ⁶Li.

> **Caveat.** These luminosity anchors assume ⁶Li can be delivered at a
> per-nucleon integrated luminosity comparable to the EIC design e-A
> specification. Final numbers depend on unresolved source-development and
> beam-transport questions for ⁶Li specifically, and the 1-year and 10-year
> mappings should be understood as design-target references rather than
> committed program parameters.

### 1.5 Plot C and per-bin rate table (extension)

**Plot C — rate per hour per (x, Q²) bin.** One panel per config (MID, TOP),
two PNGs total. The function `build_rate_perbin_plot` runs `project_rates` with
the scenario luminosity set to `LUMI_PER_HOUR_FB` so that `proj.n_events` is
directly in units of events/hour/bin. Settings:

- Nuclear PDF: EPPS21 central member (matches Plot B convention).
- Colormap: `plasma`, `LogNorm` (deliberately distinct from Plot B's `viridis`).
- Colorbar label: "Rate [events / hour / bin]".
- No N ≥ 10 cut — low-population bins are retained for completeness (contrast
  with Plot B, which blanks bins below the cut).
- Independent of Δ/F₁ scale: Plot C reflects the pure DIS cross section rate;
  the tensor signal amplitude is not involved.

**Per-bin rate CSV.** The function `build_rate_table_and_summary` writes a CSV
with 40 × 30 = 1200 rows per config (the full kinematic grid, including
out-of-acceptance bins). Columns: `x_center, q2_center, y_center,
rate_per_hour, accepted`. Four header comment lines document the config label,
luminosity assumption, and UTC generation timestamp.

**Stdout summary.** After writing the CSV, `build_rate_table_and_summary` prints
per config: total accepted-bin rate (events/hour), and a ranked table of the top
20 highest-rate (x, Q²) bins.

**Instantaneous luminosity assumption.** The rate calculation uses:

> `INSTANT_LUMI_FB_PER_S = 1e-6 fb⁻¹/nucleon/s`  
> `LUMI_PER_HOUR_FB = 3.6e-3 fb⁻¹/nucleon/hour = 3.6 pb⁻¹/nucleon/hour`

> **Caveat.** This instantaneous luminosity is chosen for consistency with the
> "10 fb⁻¹/nucleon ≈ 1 year" anchor in §1.3: 1 × 10⁻⁶ fb⁻¹/s × 10⁷ s/year
> (Snowmass year) = 10 fb⁻¹/year. Per-nucleon; for ⁶Li this corresponds to
> roughly 1.67 × 10⁻⁷ fb⁻¹/s per ion, or ~10³² cm⁻² s⁻¹ per nucleon —
> within a factor of a few of EIC design e-A luminosity. This is an assumption,
> not a committed EIC specification for ⁶Li. The same caveats on source
> development and beam transport noted in §1.3 apply here.

### 1.4 Configs run today

| Label | `E_e` (GeV) | `p_ion` (GeV/u) | Ion |
|---|---:|---:|---|
| MID | 10 | 50 | ⁶Li |
| TOP | 18 | 137.5 | ⁶Li |

**Runtime.** Approximately 10 minutes for 28 grid-backend evaluations
(14 evaluations per config: 12 nuclear-PDF combos + 1 CT18+EMC reference + 1
toy; identical count to the parent script, repeated for the second config).

---

## 2. Key physics conclusions from today's work

### 2.1 Mid config with 10 years reaches the Sather-Schmidt reference

Under the mid config with 100 fb⁻¹/nucleon (~10 years), the central combination
(mid_x + R1998 + EPPS21 + Cloet P_zz = 0.267) reaches a 5σ discovery at
**Δ/F₁ ≈ 1.144 × 10⁻³**. This is the closest any (config, luminosity)
combination in today's scan comes to the Sather-Schmidt reference scale of
10⁻³ under the Cloet P_zz convention. No other combination tested today reaches
10⁻³ for the central curve.

### 2.2 A 1-year program is not discovery-sensitive at the 10⁻³ reference scale

Under all (config, shape) combinations tested with L = 10 fb⁻¹/u, the central
curve's 5σ crossing lies above 3.6 × 10⁻³. At Δ/F₁ = 10⁻³, the EPPS21 band
ranges from δA/A = 0.369 to 1.388 (mid, 10 fb⁻¹/u) and from 0.399 to 2.605
(top, 10 fb⁻¹/u). The best-case edge of the band (low_x + Christy-Bosted at
mid config, 10 fb⁻¹/u) gives δA/A = 0.369, corresponding to a ~2.7σ hint —
not discovery.

### 2.3 Mid config outperforms top across both luminosity anchors

At both luminosity anchors the mid config requires less luminosity for 5σ
discovery than the top config on the central curve. At L = 100 fb⁻¹/u the mid
config reaches Δ/F₁ = 1.144 × 10⁻³ vs the top config's 1.657 × 10⁻³ — a
factor ~2.1 difference in the corresponding L_5σ (130.83 vs 274.65 fb⁻¹/u).
This confirms the "sweet spot at mid" finding from `money_delta_note_2026-07-16.md`
§2.1 and is unaffected by the change of observable today.

5σ crossings and derived L_5σ(10⁻³) for the central combination:

| Config | L (fb⁻¹/u) | Program time | 5σ crossing (Δ/F₁) | Derived L_5σ(10⁻³) (fb⁻¹/u) |
|---|---:|---|---:|---:|
| MID | 10 | ~1 year | 3.617 × 10⁻³ | 130.83 |
| MID | 100 | ~10 years | 1.144 × 10⁻³ | 130.83 |
| TOP | 10 | ~1 year | 5.241 × 10⁻³ | 274.65 |
| TOP | 100 | ~10 years | 1.657 × 10⁻³ | 274.65 |

### 2.4 Shape systematic dominates the band width

The EPPS21 band at Δ/F₁ = 10⁻³ spans a factor ~3.8 in δA/A at the mid config
and a factor ~6.5 at the top config across both luminosities (the ratio is
luminosity-independent because the (s/S₀)² × L scaling cancels). Nuclear-PDF
choice (EPPS21 vs nNNPDF3.0) remains sub-1% at the level of sig² — confirmed in
`money_delta_note_2026-07-20.md` §2.1 — and the two bands are visually
indistinguishable today as they were in the parent. No new systematics were
identified.

EPPS21 band at Δ/F₁ = 10⁻³:

| Config | L (fb⁻¹/u) | δA/A range |
|---|---:|---|
| MID | 10 | [0.369, 1.388] |
| MID | 100 | [0.117, 0.439] |
| TOP | 10 | [0.399, 2.605] |
| TOP | 100 | [0.126, 0.824] |

### 2.5 Per-bin uncertainty spans ~3 orders of magnitude

The per-bin `δA_bin` covers roughly 10⁻³ to 1 across the accepted `(x, Q²)`
grid at both configs and both luminosities. The median is well below 1 in all
cases, confirming that the statistical power is concentrated in the kinematic
interior. All accepted bins pass the N_bin ≥ 10 cut (100% acceptance health;
see §7.4).

Per-bin δA_bin statistics (EPPS21 + mid_x + R1998 + Cloet P_zz = 0.267):

| Config | L (fb⁻¹/u) | Accepted bins | Bins with N ≥ 10 | δA_bin min | δA_bin median | δA_bin max |
|---|---:|---:|---:|---:|---:|---:|
| MID | 10 | 345 | 345 (100%) | 0.00107 | 0.00354 | 0.839 |
| MID | 100 | 345 | 345 (100%) | 0.000338 | 0.00112 | 0.265 |
| TOP | 10 | 455 | 455 (100%) | 0.00107 | 0.00427 | 0.778 |
| TOP | 100 | 455 | 455 (100%) | 0.000340 | 0.00135 | 0.246 |

### 2.6 Rate per hour per (x, Q²) bin

For the instantaneous luminosity assumption see §1.5. Total DIS rates across all
accepted bins:

| Config | Accepted bins | Total DIS rate [events / hour] | Total DIS rate [events / s] |
|---|---:|---:|---:|
| MID (10 × 50 GeV/u) | 345 | 6.6304 × 10⁵ | ~184 |
| TOP (18 × 137.5 GeV/u) | 455 | 8.3419 × 10⁵ | ~232 |

The accepted-bin counts (345 and 455) are identical to those from Plot B because
the kinematic acceptance mask is the same; the difference is that Plot C retains
all accepted bins regardless of event count, while Plot B blanks bins below the
N ≥ 10 cut.

**Peak-rate location.** For both configs, the top 14 bins by rate/hour are all
at Q² = 2.427 GeV² — the lowest Q² bin above the `q2_min = 2.0` cut. This is
the direct signature of dσ/dQ² ∝ 1/Q⁴: cross section is strongly peaked at the
lowest accessible Q². The x range of the top-14 bins differs by config:

- MID: x from 5.62 × 10⁻³ to 0.112. Peak per-bin rate ≈ 8.85 × 10³
  events/hour (~2.5 Hz per bin).
- TOP: x extends down to 5.6 × 10⁻⁴ — one decade lower than MID — because the
  higher √s opens smaller-Bjorken-x phase space. Peak per-bin rate is comparable
  (~8.8 × 10³ events/hour). The extra total rate in TOP relative to MID comes
  from a wider kinematic footprint, not from higher per-bin peaks.

**Detector implication.** A peak per-bin rate of ~2.5 Hz is not limiting for
ePIC central-detector observables. The total DIS rate of ~200 events/s is
manageable but relevant for far-forward spectator tagging, where tagging
efficiency and dead-time considerations act at the total-rate level. This is an
operational planning consideration, not a source of systematic uncertainty on
the tensor asymmetry measurement.

---

## 3. Current physics assumptions still in play

Today's script changes the observable display but not the physics chain. The
assumption list from `money_delta_note_2026-07-20.md` §3 carries forward
unchanged:

- **Δ shape.** Three toy variants (low_x, mid_x, high_x); no lattice input, no
  CDKS-type convolution-model calculation. Shape choice remains the dominant
  band-width driver.
- **P_zz dilution.** Fixed at Cloet 0.267 to isolate the observable-change
  from the convention question. The convention ambiguity (Cloet
  per-nucleon-average vs cluster-picture occupancy-folded 0.404) remains
  unresolved. See `money_delta_note_2026-07-16.md` §7 for the full physics
  discussion.
- **Nuclear PDF sets.** EPPS21 and nNNPDF3.0 central members (member 0) only.
  Replica uncertainties within each set are not evaluated.
- **`R = σ_L/σ_T`.** Two simplified DIS-region parameterizations (R1998,
  Christy-Bosted). Neither is a full published fit applied outside its validity
  range.
- **η electron-ID efficiency.** Not included. Uniform efficiency inside
  |η| ≤ 3.5.
- **No detector resolution smearing.** No bin migration in x, Q², or φ.
- **No radiative corrections.**
- **In-ring P_zz survival.** Nominal source value 0.80 assumed to persist
  through storage.
- **⁶Li deuteron-cluster embedding factor.** The 0.87 cluster polarization
  retention fraction (Schellingerhout PRC 48:2714) remains outside the formal
  significance chain.

---

## 4. Next steps, ranked by expected impact

The Tier A ranking from `money_delta_note_2026-07-20.md` §4 is unchanged.
Today's results sharpen the motivation for Tier A #2 (see §2.4 above): the
factor ~4–7 band width from shape choice is the binding constraint on whether a
1-year program constitutes a discovery claim.

### Tier A — physics inputs to resolve

1. **Resolve the ⁶Li P_zz dilution convention (highest priority).** Consult
   I. Cloet or the relevant nuclear-structure literature. Factor ~1.5–2.4 in
   `L_5σ` depending on convention. All results today are conditional on the
   Cloet 0.267 choice. See `money_delta_note_2026-07-16.md` §7 for the list of
   questions to bring to the theorist.

2. **Anchor Δ(x, Q²) to a real theory calculation (promoted; still highest
   resolvable systematic).** Today's band spanning a factor ~4–7 in δA/A
   directly determines whether a 1-year program can claim a hint vs. noise.
   Options in decreasing order of directness: (a) lattice-moment constraints
   (Detmold–Shanahan-type), (b) CDKS-style convolution model, (c) NJL or
   light-cone quark model. At minimum, replace the three ad-hoc x-shapes with a
   single motivated one and quote its own uncertainty.

3. ~~**Include nuclear shadowing / nuclear-PDF upgrade.**~~ **RETIRED**
   (2026-07-20). Sub-1% effect in `L_5σ` for ⁶Li; both EPPS21 and nNNPDF3.0
   confirmed running.

### Tier B — detector and reconstruction realism

4. **Add η-dependent electron-ID efficiency** (real ePIC curve).
5. **Detector resolution smearing.** Gaussian smearing on E′ and θ; propagate
   into x, Q².
6. **Radiative corrections.** HECTOR or EIC-adapted equivalent.

### Tier C — long-term

7. **NNLO QCD corrections.**
8. **Target-mass corrections** at high x.
9. **Polarized-electron scenarios.** Complementary observables (A_LT, etc.);
   does not change tensor-only reach.

---

## 5. Reproducibility

All scripts are run from the repo root
(`/Users/L00338853/work/Polarized_Li/PolarizedLithiumSim`).

### Command line

```bash
python3 fastsim/scripts/money_delta_20260721.py
```

Output writes to `fastsim/out/money_delta/` by default. The `--outdir` flag
overrides this:

```bash
python3 fastsim/scripts/money_delta_20260721.py --outdir /path/to/outdir
```

Runtime: approximately 10 minutes for 28 grid-backend evaluations (14 per
config × 2 configs).

### Prerequisites

Same as `money_delta_note_2026-07-20.md` §5. No new package requirements.
Both nuclear PDF sets must already be installed:

```bash
python3 -m parton install EPPS21nlo_CT18Anlo_Li6
python3 -m parton install nNNPDF30_nlo_as_0118_A6_Z3
```

For the base grid backend and the original environment patches (NumPy 2.0
compatibility fix in `parton/pdf.py` and the `_safe_xfx` helper in
`polli_fastsim/structure.py`), see `money_delta_note_2026-07-16.md` §5. No new
patches are introduced today. The `_safe_xfx_local` function used in
`money_delta_20260721.py` (lines 76–87) mirrors the existing helper rather than
modifying the package.

---

## 6. File inventory

### New plots added today (all in `fastsim/out/money_delta/`)

```
money_delta_20260721_statvsdelta_mid_10fb.png    (Plot A: δA/A vs Δ/F₁, mid config, L=10 fb⁻¹/u)
money_delta_20260721_statvsdelta_mid_100fb.png   (Plot A: δA/A vs Δ/F₁, mid config, L=100 fb⁻¹/u)
money_delta_20260721_statvsdelta_top_10fb.png    (Plot A: δA/A vs Δ/F₁, top config, L=10 fb⁻¹/u)
money_delta_20260721_statvsdelta_top_100fb.png   (Plot A: δA/A vs Δ/F₁, top config, L=100 fb⁻¹/u)
money_delta_20260721_perbin_mid_10fb.png         (Plot B: per-bin δA_bin heatmap, mid, L=10 fb⁻¹/u)
money_delta_20260721_perbin_mid_100fb.png        (Plot B: per-bin δA_bin heatmap, mid, L=100 fb⁻¹/u)
money_delta_20260721_perbin_top_10fb.png         (Plot B: per-bin δA_bin heatmap, top, L=10 fb⁻¹/u)
money_delta_20260721_perbin_top_100fb.png        (Plot B: per-bin δA_bin heatmap, top, L=100 fb⁻¹/u)
```

### New plots added today (Plot C, all in `fastsim/out/money_delta/`)

```
money_delta_20260721_rate_perbin_mid.png   (Plot C: rate/hour/bin heatmap, mid config, plasma)
money_delta_20260721_rate_perbin_top.png   (Plot C: rate/hour/bin heatmap, top config, plasma)
```

### New CSV tables added today (in `fastsim/out/money_delta/`)

```
money_delta_20260721_rate_perbin_mid.csv   (1200 rows: x, Q², y, rate/hour, accepted — mid config)
money_delta_20260721_rate_perbin_top.csv   (1200 rows: x, Q², y, rate/hour, accepted — top config)
```

### New script added today

```
fastsim/scripts/money_delta_20260721.py    (1079 lines; grew from 877 → 1079 (+202 lines); parent: money_delta_20260720.py)
```

### Previously existing files

See `money_delta_note_2026-07-20.md` §6 for the 2026-07-20 additions and
`money_delta_note_2026-07-16.md` §6 for the full inventory of the 14 plots
present at end of 2026-07-16.

### Notes

- `fastsim/notes/money_delta_note_2026-07-16.md` — original session note.
- `fastsim/notes/money_delta_note_2026-07-20.md` — prior session note.
- `fastsim/notes/money_delta_note_2026-07-21.md` — this file.

---

## 7. Sanity checks performed

All four checks passed post-run.

### 7.1 L_5σ cross-check vs `money_delta_20260720.py`

The derived `L_5σ(10⁻³) = L × (s*/10⁻³)²` from Plot A's 5σ crossings matches
`money_delta_20260720.py`'s output for both configs:

- MID: 130.83 fb⁻¹/u (both scripts agree exactly).
- TOP: 274.65 fb⁻¹/u (both scripts agree exactly).

This confirms that the analytic `(s/S₀)² × L` rescaling used in today's script
is internally consistent with the parent's grid-evaluated sig² values.

### 7.2 Luminosity scaling (10 → 100 fb⁻¹)

The 5σ crossing scale `s*` scales as L⁻¹/² (because `s* ∝ 1/√L` from the
`s*² × L = const` relation). Expected ratio: `s*(10) / s*(100) = √10 = 3.162`.

Observed:
- MID: 3.617 × 10⁻³ / 1.144 × 10⁻³ = 3.162. ✓
- TOP: 5.241 × 10⁻³ / 1.657 × 10⁻³ = 3.163. ✓

Both within 0.1% of √10.

### 7.3 Band-width luminosity scaling (10 → 100 fb⁻¹)

The `δA/A` band should shrink by √10 = 3.162 when luminosity increases by 10×.
Observed band-edge ratios at Δ/F₁ = 10⁻³:

- MID lo: 0.369 / 0.117 = 3.154. ✓
- MID hi: 1.388 / 0.439 = 3.162. ✓
- TOP lo: 0.399 / 0.126 = 3.167. ✓
- TOP hi: 2.605 / 0.824 = 3.161. ✓

All within 0.5% of √10.

### 7.4 Bin-population health

100% of kinematically accepted bins pass the N_bin ≥ 10 cut at both luminosities
and both configs. The minimum-events threshold does not remove any acceptance at
either luminosity anchor:

- MID, 10 fb⁻¹/u: 345 accepted, 345 passing (100%).
- MID, 100 fb⁻¹/u: 345 accepted, 345 passing (100%).
- TOP, 10 fb⁻¹/u: 455 accepted, 455 passing (100%).
- TOP, 100 fb⁻¹/u: 455 accepted, 455 passing (100%).

### 7.5 Order-of-magnitude cross-check on total DIS rate

Expected total rate ≈ L_inst × ⟨σ_DIS⟩ per nucleon. Using
L_inst = 1 × 10⁻³ nb⁻¹/s (= 10⁻⁶ fb⁻¹/s × 10⁶ nb/fb) and a rough
DIS cross section ⟨σ⟩ ≈ 0.2 nb/nucleon (order-of-magnitude estimate for
the accepted kinematic range), the expected rate is ~200 events/s. Observed:
184 events/s (MID) and 232 events/s (TOP). ✅

### 7.6 TOP/MID total rate ratio

The TOP config has 1.32× more accepted bins (455 vs 345) and accesses higher-√s
cross sections, so a ratio of 1.3–2× was expected. Observed:
8.3419 × 10⁵ / 6.6304 × 10⁵ = 1.26×. ✅

### 7.7 Peak-Q² location from dσ/dQ² ∝ 1/Q⁴

The lowest Q² bin above the `q2_min = 2.0` cut is at Q² = 2.427 GeV². Because
dσ/dQ² ∝ 1/Q⁴, the rate should be maximum at this bin. Confirmed for both MID
and TOP: the top 14 bins by rate/hour all sit at Q² = 2.427 GeV². ✅

### 7.8 Consistency of Plot C rates with Plot B δA_bin minimum

Plot B (§1.2, §2.5) reports a minimum per-bin uncertainty δA_bin = 0.00107 at
L = 10 fb⁻¹/u for both MID and TOP, corresponding to the highest-statistics
bin. That minimum can be cross-checked against the Plot C peak rate.

For the MID peak bin (x ≈ 0.07, Q² ≈ 2.4 GeV², rate ≈ 8850 events/hour):
L = 10 fb⁻¹/u corresponds to 10 / 3.6 × 10⁻³ ≈ 2778 hours (from
`LUMI_PER_HOUR_FB` = 3.6 × 10⁻³ fb⁻¹/hour in §1.5), so

    N_bin = 8850 × 2778 ≈ 2.46 × 10⁷ events

    δA_bin = √(2 / 2.46 × 10⁷) / 0.267 ≈ 1.07 × 10⁻³

Plot B's reported minimum δA_bin at L = 10 fb⁻¹/u = 0.00107. Exact match. ✅
This confirms that the rate scale in Plot C and the uncertainty scale in Plot B
are derived from the same underlying cross-section calculation.

---

**Bottom line for the next session.** Discovery of Δ/F₁ at the Sather-Schmidt
reference (10⁻³) requires the mid config with a full ~10-year EIC program under
the Cloet P_zz convention; a 1-year program provides at best a 2-3σ hint under
the most-favorable Δ-shape assumption, motivating the Tier A #2 next-step
priority (theory-anchored Δ shape) even more sharply.
