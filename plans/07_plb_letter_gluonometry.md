# Plan 07 — Simulation Letter (PLB-class): Nuclear Gluonometry with Tensor-Polarized ⁶Li at the EIC

**Goal.** Turn the development-run-5 results (money plots 5–7, the coherent
intact-recoil channel, the fact-checked report) into a submitted letter in a
Physics Letters B–class journal, on the arXiv before the INT program
"Towards Realizing the Program with Polarized Ion Beams at EIC"
(March 22 – April 2, 2027).

**Status legend:** ☐ todo ◐ started ☑ done · **Decision points are marked D1–D5.**

---

## 7.0 Scope decision (D1)

One paper, one message. Candidates, ranked:

| candidate | novelty | risk | verdict |
|---|---|---|---|
| **Gluonometry Δ(x,Q²): inclusive + coherent** | first numerical EIC projection for any target (verified gap); first intact-A=6 tag projection; the ⁶Li null test (Q = −0.0806 fm²) is a quotable idea | coherent inputs are scenario bands; optics undocumented for Li | **RECOMMENDED** |
| Gluonometry inclusive only | same first, smaller surface | loses the null test and the imaging connection | fallback if referees or co-authors balk at the coherent bands |
| b₁(⁶Li) via A_zz | first A > 2 b₁ projection | no theory prediction for ⁶Li b₁ → interpretation section is thin; HERMES/JLab context makes it "another projection" | second paper, after theory input |
| Polarized EMC g₁(⁷Li) | strong physics case | E12-14-001 exists; message is incremental precision → better as PRC/EPJ A with reco-level detail | third paper / long-form companion |
| polligen (the generator itself) | first polarized e+A generator | software papers fit CPC/EPJ C, not PLB | companion software paper, later |

The recommended letter carries three claims:

1. The EIC with a transversely tensor-polarized ⁶Li beam measures the
   double-helicity-flip structure function Δ(x,Q²) — never measured for any
   target — with per-bin δA = O(10⁻⁴) in one year: not merely detection, but
   discrimination between moment-constrained interpretations (A vs B) by the
   x dependence, and an extraction of xΔ(x,Q²) itself (money plot 7).
2. The same observable exists in coherent diffraction with the ⁶Li detected
   intact; ⁶Li's anomalously small quadrupole moment makes it a null test that
   separates nuclear-shape modulations from exotic glue — with a predicted
   sign flip of the deformation term relative to the polarized-deuteron CGC
   benchmark.
3. The far-forward analysis (beam-rigidity blindness, the pT-tail tag, the
   IR-8 secondary focus) defines what the machine must provide — the letter
   doubles as the physics case for Li beam optics and the second IR.

Excluded from this letter (and where it goes): tagged A_zz and ⁷Li
polarimetry (paper 2 with b₁ theory); polarized EMC (paper 3); generator
methodology beyond one validation paragraph (software companion).

## 7.1 Headline numbers the letter will quote

All to be re-derived under WP1–WP3 (grid structure functions, reco level);
current values from the 2026-08-17 report shown as the baseline:

| quantity | baseline (toy SFs, generator level) |
|---|---|
| sweet-spot amplitudes (moment_A, P_zz = 0.6, dilution 1/3) | (0.7–1.2)×10⁻² |
| per-bin δA, 1 yr / 10 yr (10/100 fb⁻¹/u) | (1.5–4.5)×10⁻⁴ / ÷√10 |
| significance per bin, year 1 | ≥ 21σ |
| extracted δΔ per bin, year 1 | (0.7–1.4)×10⁻³ on Δ ≈ −0.03…−0.09 |
| coherent tag acceptance (IP6, 0.20 GeV cut, B = 50) | 13.5% [9–20% for B ∈ 40–60] |
| N_tag, 1 yr / 10 yr (f₀ = 0.04 [0.02–0.08]) | 1.1×10⁷ / 1.1×10⁸ |
| coherent best-bin δA, 1 yr / 10 yr | 1.9×10⁻³ / 6×10⁻⁴ (5σ floors 1.0%/0.3%) |
| deformation-anchored ⟨a₂⟩ at P_zz = 0.6 | 0.036 [0.018–0.059], sign flip vs d |

**Superseded 2026-08-27.** The baseline above was computed at the
rigidity-scaled energies. Report 1 as rewritten on 2026-08-27 carries the
current values at e10 × ⁶Li 99.5 GeV/u: sweet spots x = 0.011–0.14,
amplitudes (0.44–0.95)×10⁻² vs δA = (1.4–4.5)×10⁻⁴ (21–44σ), δΔ 2–9%
relative in the best bins, N_tag = 1.7×10⁷ / 1.7×10⁸ for a 0.20 GeV
envelope with best-bin δA = 1.8×10⁻³ / 0.6×10⁻³ — and the finding that no
published IP6 optics delivers that envelope (plans/10), which turns claim 3
into the letter's far-forward requirement.

#### §7.1 re-derived on the two-fill estimator (2026-09-15)

The three truth-level drivers (`money_cos2phi.py`,
`money_cos2phi_coherent.py`, `money_delta_extraction.py`) now carry the
two-fill path WP3 asked for: `--pzz-plus/--pzz-zero`, both defaulting to
`None`, which split the same total luminosity between an m = ±1-rich fill
and an m = 0-rich one and read them with `reco.harmonic_ratio_fit` — the
spin-state-sorted ratio the reconstructed-level chain has used since
2026-08-28 — instead of the single-fill binned fit. With the flags unset
every printed line and every published PNG is bit for bit what it was; a
run with them on writes its own `_twofill` stem. Equal luminosity in the
two fills is not a flag, because the 0.67 rests on it and a free share
under one stem key would put two numbers behind one file name.

**The factor, measured three ways.** The closed form is
δÂ(two-fill)/δÂ(single) = P₊/σ_P = 0.6/0.9 = 0.6667 at (+0.6, −1.2). At the
four sweet spots the drivers print the ratio of the *fitted* errors as
**0.6665 / 0.6666 / 0.6666 / 0.6668** — the estimator delivers its own
algebra with no reconstruction-level loss, which is the statement the
generator level can make and the reconstructed level cannot (there
0.60–0.72, Report 1 §5.2, because efficiency and purity ride on top).
Â stays unbiased: the four pulls against the two-fill truth are
−0.4, +0.9, −1.2 and −0.8 σ.

**What moves in the table above, and by how much** (mid configuration,
toy backend, the published seed, so that only the estimator changes):

| §7.1 row | single fill (published) | two fills (+0.6, −1.2) | factor |
|---|---|---|---|
| sweet-spot amplitudes | (0.44–0.95)×10⁻² | unchanged | 1 |
| per-bin δA, 1 yr | 1.7 / 1.4 / 2.7 / 4.5 ×10⁻⁴ | 1.2 / 0.95 / 1.8 / 3.0 ×10⁻⁴ | 0.667 |
| per-bin δA, 10 yr | 5.5 / 4.5 / 8.7 / 14 ×10⁻⁵ | 3.6 / 3.0 / 5.8 / 9.4 ×10⁻⁵ | 0.667 |
| significance per bin, year 1 | 21–43σ | **32–64σ** | 1.5 |
| relative δΔ, best bin per Q² slice, year 1 | 2.5 / 2.4 / 9.5% | **1.7 / 1.6 / 6.3%** | 0.667 |
| coherent best-bin δÂ, 1 yr / 10 yr | 1.8 / 0.57 ×10⁻³ | **1.2 / 0.38 ×10⁻³** | 0.667 |
| coherent 5σ floors, 1 yr / 10 yr | 0.90% / 0.28% | **0.60% / 0.19%** | 0.667 |
| coherent tag acceptance, N_tag, ⟨a₂⟩ | — | unchanged | 1 |

Everything that moves moves by the one factor, because the estimator
changes the error and nothing else: the amplitudes, the tag acceptance,
N_tag and ⟨a₂⟩ are properties of the beam and the model, not of how the
fills are sorted. The 1.5× gain is real **only** if the source delivers
m = 0-rich bunches at |P_zz| = 1.2 at the same purity as the ±1 fill at
0.6 — the condition WP3 attached to it and the one the run plan below
prices. `two_fill_err_ratio` makes the dependence explicit: at (+0.6,
−0.6), equal purity with no m = 0 enrichment, the two-fill error equals
the single-fill one, and at (+0.6, 0) it is twice as large.

Which footing §7.1 and Report 1 quote is not settled here, for the same
reason the WP1 grid column is not: it is an author call, and the two calls
are coupled (the drift table above is on single-fill errors). Report 1
§3.3 records that the two paths now come from one code path; its published
single-state numbers are untouched.

#### The Yellow-Report binning, priced (2026-09-15)

plans/02 Step 1.1 item 3 asked for the YR inclusive convention, ~5
logarithmic bins per decade, in place of the generic 40 × 30 log grid.
`beams.YR_GRID` / `YR_X_EDGES` / `YR_Q2_EDGES` hard-code it as a
decade-anchored lattice — edges at 10^(k/5), so the bins are the published
ones and not a five-way split of this project's span — and
`phase_space_bins.py` and `coverage_and_stat_maps.py` reach it through
`--binning {log40x30,yr}`, default `log40x30`, which is bit for bit the
published grid (the four money-plot stems and all twelve coverage PNGs
verified by md5 before and after). A YR run writes its own `_yr` stem.

The x span 1e−4…1 is four exact decades, 20 bins. The Q² span 1…2×10³ is
3.301 decades, so the lattice edge above it is 10^3.4 = 2512 GeV² and the
YR range is rounded **out** to it (17 bins) rather than in, so that no part
of the published span is dropped. The extra row costs nothing: it holds
zero accepted events at every configuration of both isotopes, the electron
acceptance cutting off well below Q² = 2×10³ (pinned in
`evgen/tests/test_money_scripts.py`).

**The drift, mid configuration, one EIC year.** Accepted cells 461 → 124;
N_DIS 5.193 → 5.236 ×10⁹ (+0.8%, entirely the coarser grid's cell-centre
acceptance test, none of it the extra Q² row); median events per bin
3.00×10⁶ → 1.24×10⁷ (×4.1); max per-bin cos 2φ significance 24.8 → 46.8
(×1.89) and min per-bin δA 2.98 → 1.62 ×10⁻⁴ (×0.55, i.e. 1/1.83) — both of
them pooling and not a gain, of the order of the √4.1 the median cell gains
and short of it because the extremal cells are not the median one. The rate drift is configuration-dependent and grows with beam energy:
N_DIS ×1.005 / ×1.008 / ×1.079 at 5 × 40.8, 10 × 99.5 and 18 × 137.5 ⁶Li
(×1.005 / ×1.008 / ×1.094 for ⁷Li), because the coarse lattice straddles
the acceptance boundary more crudely where the boundary is longest. The
coherent map drifts the other way — N_coh 1.234 → 1.180 ×10⁸ and N_tag
1.669 → 1.598 ×10⁷, −4.4 and −4.3% — because f_coh(x) is steep across a
5-per-decade x bin and is evaluated at its centre.

**What the analysis binning does to the money plots.** The sweet spots
move to the nearest YR cells, (0.028, 1.14) → (0.032, 1.26), (0.011, 1.14)
→ (0.013, 1.26), (0.071, 3.13) → (0.079, 3.16) and (0.141, 14.3) →
(0.126, 12.6); the Δ-extraction combs drop from 10 / 9 / 10 merged x-bin
pairs per Q² slice to 5 / 5 / 5; and the tagged coherent super-bin, 3 × 2
cells on both grids (the pad-1 window is clipped against the Q² = 1 floor),
grows by ×3.02 in yield (1.75 → 5.28 ×10⁶). The last two are the
reason the YR grid is a comparison convention rather than a replacement:
it is the grid on which a YR projection can be read off bin for bin, and
the grid on which this measurement has the fewest Δ(x) points. Adopting it
for the published figures is therefore a separate call from providing it,
which is why `--binning` defaults to the grid the figures were made on.

## 7.2 Venue

- **Primary: Physics Letters B.** Precedents for exactly this genre:
  Friščić et al., PLB 823 (2021) 136726 (e+³He tagging projections);
  Mäntysaari et al., PLB 858 (2024) 139053 (polarized-d imaging proposal).
  Format: elsarticle, no hard length limit but letters run 6–8 published
  pages; we budget ~4,200 words + 4 figures + 1 table + ~40 references.
- Fallbacks, in order: EPJ A (EIC projection papers routine; no length
  pressure), PRD (if the coherent section grows), PRL only if a co-author
  push and a compressed 4-page cut emerge (not planned).
- arXiv: hep-ex primary, nucl-ex + hep-ph cross-list. Repository citation +
  Zenodo DOI minted at submission; the GitHub Pages report becomes the
  "extended companion note" link.

## 7.3 Gap analysis — what referees will check vs what we have

| # | item | current state | required for submission | WP | blocking? |
|---|---|---|---|---|---|
| 1 | structure functions | TOY in money plots 5–7 (grids exist behind interface; money_delta line ran EPPS21) | production plots on EPPS21nlo_CT18Anlo_Li6 + R1998; toy kept as cross-check band | WP1 | **yes** |
| 2 | polarization placeholder | single P_zz = 0.6, full luminosity in one fill | δA vs P_zz ∈ {0.4, 0.6, 0.8} and fill-share ∈ {0.5, 1}; one table | WP2 | **yes** |
| 3 | dilution convention (plans/04 #6, closed 2026-08-29) | ⅓ per-nucleon baseline; the ×2.4 cluster upside is withdrawn — the cluster value settled at 0.81123 is the VECTOR polarization of the pair, and Δ is rank 2, where the same wave function transfers 0.92 (plans/08 D9) | state the ⅓ and its rank in one paragraph, with no upside claimed; whether Δ should carry the rank-2 0.92 as b₁ does (⅓ → 0.3073) is a question for the Δ registry, not for the convention | WP2 | no longer blocking (one paragraph) |
| 4 | detector level | generator-level φ′; smearing machinery exists (money_delta_20260729: ePIC tracking + ε_eID, Case-3 retains ~92% S/N) | reco-level amplitude dilution + δA inflation per super-bin; 2-D φ-hole closure (plans/03 2.3.3) | WP3 | **yes** |
| 5 | radiative corrections | uncalculated for tensor observables (flagged); **the collinear-ISR migration is bounded — Δ̂ by +0.5 to +1.2% in the published generator window, ≤ 2.9% with the low-Q² feed-in opened up, ≤ 0.25% behind an E − p_z window (a documented contingency, not applied), against the 5% gate (WP4 closed 2026-08-28)** | quantified migration bound + spin-state-ratio cancellation argument (common acceptance, bunch-by-bunch) + explicit open-theory statement | WP4 | yes (as a bound, not a calculation) |
| 6 | coherent optics cut | two-point (0.20/0.45 GeV); Li optics undocumented (report assumption #1) | acceptance, N_tag, δA as curves vs pT_cut ∈ [0.1, 0.7] GeV; IR-8 alternative with published efficiencies; t_min and ×0.73 rate-weighting folded into the central curve | WP5 | **yes** |
| 7 | far-forward geometry | RP z = 26/28 m (YR-era) in text; run-2 Geant4 scan found 32.5/34.3 m in epic-main | quote current ePIC geometry, windows unchanged in θ/R | WP5 | minor |
| 8 | coherent fraction f₀ | scenario 0.04 [×2/÷2] | unchanged (no light-A prediction exists — stated); theory ask on record (IP-Glasma α–d, #18) | — | no |
| 9 | BeAGLE incoherent shapes | blocked on FLUKA license | not blocking: |t|-fit purity argument stays band-level with e+Pb benchmark [80–99%] | WP7 | no |
| 10 | bibliography | 40 entries, individually verified 2026-08-17 | BibTeX from INSPIRE keys; re-verify volume/pages at submission | WP6 | no |

## 7.4 Work packages

### WP1 — Grid structure functions in production (◐ infrastructure exists; R done)
- ☑ **R is done** (2026-08-26, plans/08 C2): the published SLAC/E143 R1998
  world fit lives in `polli_fastsim/structure.py` as `r1998` (all three of
  the paper's forms, so their spread is the fit's functional-form
  systematic) and reaches every consumer through one `r_func` hook.  The
  defect it replaces clipped R to exactly 1.000 over 38% of the sensitivity
  box.  Measured at the four sweet spots: **F₂ cancels exactly** in the
  cos 2φ amplitude (2×10⁻¹⁶) while Δ/F₁ moves +16.6 / +18.0 / +4.7 / −4.4%,
  so of the two structure-function inputs only R moves the physics result —
  which reorders this work package.
- ☑ The grids themselves are installed (CT18NLO, EPPS21nlo_CT18Anlo_Li6,
  NNPDFpol11_100 via `parton`); the fast-sim grid tests no longer skip.
- ☑ Promote the money_delta script-local `NuclearF2FromGrid`
  (EPPS21nlo_CT18Anlo_Li6) into `polli_fastsim/structure.py`
  behind the existing `NuclearF2` interface; wire through
  `inputs.get_backends` (pass `r_func=structure.r1998` there, not a
  monkey-patch).
  ☑ *2026-09-15: `structure.NuclearF2FromGrid` carries the `NuclearF2`
  signature `(ion, setname=None, member=0, r_func=None)` — the set
  resolved off (Z, A) through `NUCLEAR_F2_SET_FOR_ION`, the module's own
  `_safe_xfx`, and the call-time `r_sigma_lt` lookup that keeps the dated
  scripts' `r_override` working. `get_backends(pdf="toy", nuclear=None,
  r_func=None)` gains a `"nuclear"` key and threads R into the g1 model
  as well; both new arguments default `None`, at which the dict is the
  three keys it always was. `tests/test_grids.py` pins the promoted class
  against the still-frozen copy in `money_delta_20260729.py` at seven
  (x, Q²) points on both f2a and f1a: max |new/frozen − 1| = **0**
  (tolerance 1e-12).*
- ☑ Rerun money plots 5/6/7 and `phase_space_bins.py` with `--pdf grid`;
  re-solve moment_A at the grid ⟨Q²⟩; record sweet-spot drift and the
  toy-vs-grid ratio per headline number (expect ≲ ×1.5 by run-2/3 checks:
  toy-vs-CT18 F₂ within ±37%).
  ☑ *2026-09-15: done and tabulated in the addendum below. The flag is
  `--pdf {toy,grid}`, not the `--backend grid` this line asked for until
  today: no `--backend` flag has ever existed in this repository, while
  six scripts — five in `fastsim/scripts` and
  `evgen/scripts/target_mass_bound.py` — already spell it `--pdf`. The
  drift is inside the ≲ ×1.5 anticipated here: every level-quantity ratio
  of the table below lies in 0.56–1.25, and the widest excursion anywhere
  in it is the ×1.56 of a RELATIVE error (δΔ/Δ̂ at Q² = 1.14 GeV², where
  the denominator itself falls by 0.564). But the drift's largest single
  cause is not the one this work package expected.*

#### WP1 addendum — the toy → grid drift table (2026-09-15)

`--pdf grid` on the four truth-level drivers means: F₂ᴬ from the nuclear
set EPPS21nlo_CT18Anlo_Li6 (not Z·F₂ᵖ + N·F₂ⁿ on a free proton), g₁ from
NNPDFpol11_100, F₂ᵖ from CT18NLO for the g₁ denominator, and the SLAC/E143
R1998 fit in place of the toy R — one flag moving all four R consumers at
once, which is what `inputs.get_backends(pdf, nuclear, r_func)` is for.
Commands (from `evgen/`, adding `--pdf grid` to the published ones of
`docs/reproduction_manual.md` §4.1–4.2):

```bash
python3 scripts/phase_space_bins.py        --pdf grid --outdir .
python3 scripts/money_cos2phi.py           --pdf grid --outdir .
python3 scripts/money_cos2phi_coherent.py  --pdf grid --outdir .
python3 scripts/money_delta_extraction.py  --pdf grid --outdir .
```

Each writes a `_grid` stem beside the published one, which is bit-for-bit
unchanged (md5 before and after the run-19 edit: `ad907a83…`, `aae59967…`,
`2d1eae23…`, `b9a05445…`). The bin selection does not move: all four
super-bins, all three Q² slices and the tagged super-bin come out at the
same edges on both backends, so the table below compares the same bins.

**THE Q₀ FLOOR, and why it is the headline of this addendum.**
EPPS21nlo_CT18Anlo_Li6 begins at Q = 1.3 GeV, i.e. **Q² = 1.69 GeV²**, and
`parton` returns NaN below it rather than freezing or extrapolating. The
money maps start at Q² = 1 GeV², and **two of the four published sweet
spots sit at Q² = 1.14 GeV²** — below the set's support. `NuclearF2FromGrid`
therefore freezes F₂ᴬ at Q₀² below the floor (the treatment `r1998` already
gives its own fit support) and exposes `q2_min` and `q2_frozen_fraction`;
a bare backend without that guard poisons the low-Q² half of the map with
NaN and the sampler dies in `rng.poisson`. Measured at 10 × 99.5 GeV/u:
**8.7% of the accepted cells but 36.3% of the accepted one-year rate** lie
below Q₀². So the grid column at Q² = 1.14 — spots 1 and 2, the Q² = 1.14
Δ slice, and the whole coherent best super-bin, which is Q² ∈ [1, 1.66] —
is a **frozen-Q₀ continuation, not a grid evaluation**, and must be quoted
as such. Closing that gap needs a nuclear set with a lower Q₀, or the
sensitivity box raised to Q² > 1.69 GeV².

| quantity (10 × ⁶Li 99.5 GeV/u, P_zz = 0.6, 10 fb⁻¹/u) | toy | grid | ratio |
|---|---|---|---|
| moment_A ⟨Q²⟩ [GeV²] | 4.474 | 5.069 | 1.133 |
| moment_A bag amplitude A | −0.2924 | −0.2404 | 0.822 |
| N_DIS accepted, 1 yr | 5.193×10⁹ | 3.980×10⁹ | 0.766 |
| spot 1 (x 0.0282, Q² 1.14†) N | 1.903×10⁸ | 1.605×10⁸ | 0.843 |
| … A_truth | 7.416×10⁻³ | 5.217×10⁻³ | 0.703 |
| … δA 1 yr / 10 yr | 1.728 / 0.547 ×10⁻⁴ | 1.882 / 0.595 ×10⁻⁴ | 1.089 |
| … significance, 1 yr | 42.9σ | 27.7σ | 0.646 |
| spot 2 (x 0.0112, Q² 1.14†) N | 2.812×10⁸ | 2.080×10⁸ | 0.740 |
| … A_truth | 4.346×10⁻³ | 3.040×10⁻³ | 0.699 |
| … δA 1 yr / 10 yr | 1.422 / 0.450 ×10⁻⁴ | 1.653 / 0.523 ×10⁻⁴ | 1.163 |
| … significance, 1 yr | 30.6σ | 18.4σ | 0.602 |
| spot 3 (x 0.0708, Q² 3.13) N | 7.519×10⁷ | 7.798×10⁷ | 1.037 |
| … A_truth | 9.491×10⁻³ | 7.320×10⁻³ | 0.771 |
| … δA 1 yr / 10 yr | 2.750 / 0.869 ×10⁻⁴ | 2.700 / 0.854 ×10⁻⁴ | 0.982 |
| … significance, 1 yr | 34.5σ | 27.1σ | 0.785 |
| spot 4 (x 0.141, Q² 14.3) N | 2.852×10⁷ | 3.190×10⁷ | 1.119 |
| … A_truth | 9.517×10⁻³ | 8.197×10⁻³ | 0.861 |
| … δA 1 yr / 10 yr | 4.464 / 1.412 ×10⁻⁴ | 4.221 / 1.335 ×10⁻⁴ | 0.946 |
| … significance, 1 yr | 21.3σ | 19.4σ | 0.911 |
| Δ̂ at Q² 1.14†, x 0.02 | −0.1352 ± 0.0034 | −0.0763 ± 0.0030 | 0.564 / 0.88 |
| Δ̂ at Q² 3.13, x 0.0501 | −0.0696 ± 0.0016 | −0.0519 ± 0.0016 | 0.746 / 1.00 |
| Δ̂ at Q² 14.3, x 0.316 | −0.0047 ± 0.0004 | −0.0050 ± 0.0005 | 1.064 / 1.25 |
| relative δΔ in the three best bins | 2.5 / 2.3 / 8.5 % | 3.9 / 3.1 / 10.0 % | 1.56 / 1.34 / 1.18 |
| N_coh produced, 1 yr | 1.234×10⁸ | 8.604×10⁷ | 0.697 |
| N_tag, 1 yr (0.20 GeV envelope) | 1.669×10⁷ | 1.164×10⁷ | 0.697 |
| coherent best super-bin N† | 1.75×10⁶ | 1.14×10⁶ | 0.651 |
| coherent best-bin δÂ, 1 yr / 10 yr† | 1.8 / 0.6 ×10⁻³ | 2.2 / 0.7 ×10⁻³ | 1.22 / 1.17 |
| coherent 5σ floor† | 0.0090 | 0.0112 | 1.24 |
| ⟨a₂⟩_tag; coherent tag acceptance | 0.036 | unchanged by construction (both are model numbers with no structure function in them) | 1.000 |

† frozen-Q₀ rows: Q² < 1.69 GeV², where F₂ᴬ is held at Q₀².

**Which half of the drift is F₂ and which is R.** A third run with the grid
F₂ᴬ but the toy R (a scratch diagnostic, not a flag) separates them on the
four sweet-spot amplitudes. F₂ alone: 0.829 / 0.831 / 0.829 / 0.829 —
flat, because the moment constraint re-solves A against ∫x Δ dx and the
whole amplitude follows the solved bag amplitude (−0.2924 → −0.242, 0.828),
not the local F₂. R alone: 0.849 / 0.842 / 0.930 / 1.039, and after
dividing out the residual 0.993 of the bag re-solve, 0.855 / 0.848 / 0.936 /
1.046, against the reciprocals 0.858 / 0.847 / 0.955 / 1.046 of the
+16.6 / +18.0 / +4.7 / −4.4% Δ/F₁ shifts the ☑ R row above already
measured: the same numbers to 0.02–2%, reached on a different backend and
through the whole extraction rather than at the four points. **The ☑ R
row's conclusion survives and sharpens:** the R swap is the x-dependent
half of the drift and the one that moves the physics shape, while the
nuclear grid contributes an almost x-independent overall 0.83 that a
re-solved normalisation absorbs. What the grids do move on their own is
the RATE — N_DIS ×0.766, N_coh and N_tag ×0.697 — and through it every δA
and every significance.

**Independent check.** The sibling generator `LiPolGen` measures the same
switch in C++ on its own observables (`README.md` selector table,
`--unpol-sf {toy,mstw,ct18nlo}` at ⁶Li config 1): accepted σ ×0.7985
(ct18nlo) / ×0.7934 (mstw). This addendum's N_DIS ratio is **×0.766** on a
nuclear set with R1998 folded in — the same direction and within 4% of an
independently written generator's number.

- Acceptance: all §7.1 numbers re-derived on grids; drift table in the plan
  addendum. Effort: 2–4 days.
  ◐ *2026-09-15: the drift table is above and every §7.1 row it feeds is
  measured on both backends. What is NOT done here is the re-quote itself:
  §7.1 and Report 1 still carry the toy column, deliberately — which
  column becomes the published one is an author call, and it is coupled to
  the Q₀ floor above, since two of the four sweet spots cannot be
  evaluated on the nuclear grid at all. The §7.1 rows the grid column
  would change: amplitudes (0.44–0.95) → (0.30–0.82)×10⁻²; δA
  (1.4–4.5) → (1.7–4.2)×10⁻⁴; significance 21–43σ → 18–28σ; relative δΔ in
  the best bins 2–9% → 3–10%; N_tag 1.7×10⁷/1.7×10⁸ → 1.16×10⁷/1.16×10⁸;
  coherent best-bin δÂ 1.8/0.6 → 2.2/0.7 ×10⁻³ with the 5σ floors 0.9%/0.3%
  → 1.1%/0.4%. ⟨a₂⟩ and the tag acceptance are model numbers and do not
  move.*

### WP2 — Polarization, run-plan, and dilution bands
- ☑ δA scaling table vs P_zz and fill share (analytic 1/(P_zz√(fN)) checked
  against one sampler rerun).
  ◐ *2026-08-28: the fill-share half is done and priced below — the share is a
  flag (`--run-share` / `--lumi-fraction`, default 1.0), the 1/√f law is verified
  against direct reruns rather than assumed, and Plans A, B and A×B are tabulated.
  The P_zz half of the same table is still to write.*
  ☑ *2026-09-15: the P_zz half is written — `fastsim/scripts/wp2_pzz_table.py`,
  tabulated below. It scans P_zz on the grid the fill-share half used (the three
  ⁶Li configurations combined over Q², the five x bins Plan B quotes δA_zz at,
  the shares f = 1, ½, ⅓) and divides each projected δA_zz by the analytic law;
  the ratio is 1.0000 in all 120 cells, |ratio − 1| ≤ 2.2×10⁻¹⁶. The table is
  a table, not a figure: the script writes no PNG and nothing published moves.*
- ☑ Dilution paragraph: ⅓ baseline, no upside claimed.
  ☑ *2026-08-28: the paragraph is written (Report 1 §3.2 and assumption row 3)
  and ⅓ is the code default. 2026-08-29: the ×2.4 upside is withdrawn there and
  in Report 1 §4. plans/04 #6 closed on the cluster picture, but at 0.81123 and
  as a VECTOR polarization: the gap the upside was quoted from had already
  shrunk to 1.23 with the per-nucleon convention, and Δ is rank 2, where the
  α–d wave function transfers 0.92 rather than 0.87 (plans/08 D9). No Δ
  amplitude moves.*
- Acceptance: Table 1 of the letter exists. Effort: 1 day.

#### The P_zz half of the table (2026-09-15)

δA_zz per x bin, three ⁶Li configurations and every accepted Q² cell above
the 100-event floor combined, 10 fb⁻¹/u at f = 1, toy inputs, statistical only
(`python3 fastsim/scripts/wp2_pzz_table.py`):

| P_zz | x = 0.0035 | 0.0089 | 0.0282 | 0.2818 | 0.5623 |
|---|---|---|---|---|---|
| 0.30 | 1.86×10⁻⁴ | 1.77×10⁻⁴ | 1.95×10⁻⁴ | 7.09×10⁻⁴ | 2.34×10⁻³ |
| 0.40 | 1.40×10⁻⁴ | 1.33×10⁻⁴ | 1.46×10⁻⁴ | 5.32×10⁻⁴ | 1.75×10⁻³ |
| 0.50 | 1.12×10⁻⁴ | 1.06×10⁻⁴ | 1.17×10⁻⁴ | 4.26×10⁻⁴ | 1.40×10⁻³ |
| 0.60 | 9.31×10⁻⁵ | 8.84×10⁻⁵ | 9.73×10⁻⁵ | 3.55×10⁻⁴ | 1.17×10⁻³ |
| 0.70 | 7.98×10⁻⁵ | 7.57×10⁻⁵ | 8.34×10⁻⁵ | 3.04×10⁻⁴ | 1.00×10⁻³ |
| 0.80 | 6.98×10⁻⁵ | 6.63×10⁻⁵ | 7.30×10⁻⁵ | 2.66×10⁻⁴ | 8.77×10⁻⁴ |
| 0.90 | 6.21×10⁻⁵ | 5.89×10⁻⁵ | 6.49×10⁻⁵ | 2.36×10⁻⁴ | 7.79×10⁻⁴ |
| 1.00 | 5.59×10⁻⁵ | 5.30×10⁻⁵ | 5.84×10⁻⁵ | 2.13×10⁻⁴ | 7.02×10⁻⁴ |

Any other share multiplies the whole table by 1/√f; the accepted counts the
errors are built on are N = 6.41 / 7.12 / 5.86 ×10⁸, 4.42×10⁷ and 4.06×10⁶ in
the five bins at f = 1.  The law the work package names is met exactly:
writing D and A for the two amplitude attenuations of `fom.Scenario`
(plans/05:370, both 1.0 here),

δA_zz = √2 / (P_zz · D · A · √(f N)),

and measured/analytic = 1.0000 at every one of the 120 cells, to 2.2×10⁻¹⁶
— the ratio column the acceptance asks for, and it is exact rather than
approximate because nothing in the projection except the three error
functions sees P_zz at all.  Three things the table settles rather than
assumes:

* **The P_zz = 0.6 and 0.8 rows are the published ones.**  They reproduce
  `money_b1.py`'s two δA_zz rows digit for digit at all five x, so the new
  scan is the same projection the published figure is drawn from and not a
  second implementation of it.
* **One table serves both tensor observables.**  `err_azz` and
  `err_cos2phi_amplitude` are the same function of (N, P_zz), √(2/N)/P_zz,
  so the cos 2φ amplitude of this letter scales with P_zz exactly as A_zz
  does; the script asserts the two agree cell by cell rather than trusting
  it, and stops if they ever part.
* **The min-events floor is not part of the law.**  A cell enters the
  combination at ≥ 100 events, which at f < 1 is a different set of cells.
  Over the grid here the floor is inert — `--floor cell` and the default
  `--floor reference` agree to 2.2×10⁻¹⁶ down to f = ⅓ — but it does bite
  eventually: at f = 10⁻³ the x = 0.5623 bin departs from the law by 15% and
  x = 0.2818 by 0.75%, entirely through the selection.  A run-plan table that
  went below a few percent of the year would have to quote the floor with the
  share.

#### The run plan, priced (2026-08-28)

Today every projection in every report gives its observable the **whole**
10 fb⁻¹/u year in its own spin configuration, far-forward optics and
isotope.  That is stated once (Report 0 §6, Report 1 §3.1, Report 2
Table 2, Report 3 Table 9 row 16) and it means the reaches are
alternatives, not a programme.  The share is now a flag — `--run-share`
in `fastsim/`, `--lumi-fraction` in `evgen/`, both default 1.0, both
carried into `fom.Scenario.run_share` — so the options below are
arithmetic, not a rerun of the physics.  **This section prices them; it
does not choose one, and none of these numbers belongs in a report until
the programme decides.**

Three laws, verified exactly (§3.5 of the manual,
`fastsim/tests/test_run_share.py`, `evgen/tests/test_run_share.py`): a
share *f* multiplies every statistical error by 1/√*f*, leaves any
luminosity quoted as a reach (`L_5σ` = 16.7 / 16.3 / 21.8 fb⁻¹/u for the
toy inclusive Δ) exactly where it is, and multiplies the years to a
target significance by 1/*f*.

**Plan A — one ⁶Li year per configuration, split between the coherent and
the inclusive channel.**  The coherent channel is given exactly the share
that buys 5σ on the deformation (shape) term in the year:
*f*<sub>coh</sub> = (5/9.4)², (5/8.3)², (5/10.7)² = 0.283 / 0.363 / 0.218
from the 9.4 / 8.3 / 10.7 σ per year of `tagging_optics.py` at the
tagging optics.  The rest runs at the Yellow Report high-acceptance
optics and is shared *s* = 1, ½, ⅓ among the three inclusive ion fills.

| configuration | *f*<sub>coh</sub> | δa_e coherent, 1 yr at *f*<sub>coh</sub> | inclusive share, *s* = 1 | δ(Δ/F₁) | yr to 5σ | *s* = ½ | δ(Δ/F₁) | yr to 5σ | *s* = ⅓ | δ(Δ/F₁) | yr to 5σ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 × 40.8 | 0.283 | 0.00227 | 0.717 | 3.054×10⁻⁴ | 2.33 | 0.358 | 4.319×10⁻⁴ | 4.66 | 0.239 | 5.290×10⁻⁴ | 6.99 |
| 10 × 99.5 | 0.363 | 0.00184 | 0.637 | 3.202×10⁻⁴ | 2.56 | 0.319 | 4.529×10⁻⁴ | 5.13 | 0.212 | 5.546×10⁻⁴ | 7.69 |
| 18 × 137.5 | 0.218 | 0.00158 | 0.782 | 3.341×10⁻⁴ | 2.79 | 0.391 | 4.725×10⁻⁴ | 5.58 | 0.261 | 5.786×10⁻⁴ | 8.37 |

δa_e is the combined one-year error on the exotic-glue coefficient over
the seven |t| bins, 0.00121 / 0.00111 / 0.00074 at the full year
(`money_cos2phi_coherent_reco.py --optics tagging --exact --n-mc 600000`,
per configuration), divided by √*f*<sub>coh</sub>.  δ(Δ/F₁) is the error
on the Δ/F₁ scale after one programme year, 2.586 / 2.556 / 2.954 ×10⁻⁴
at the full year, and "yr to 5σ" the programme years to 5σ on
Δ/F₁ = 10⁻³, 1.67 / 1.63 / 2.18 at the full year
(`money_delta.py --ion 6Li --pdf toy`).  Every entry was **also produced
directly** by rerunning `money_delta.py --run-share` at the share in its
own row, taken unrounded as (1 − *f*<sub>coh</sub>)·*s* — the share
column is rounded to three decimals for display and rerunning at the
rounded value reproduces the entries only to two significant figures.
The direct and the rescaled values agree to the last digit printed
(5.786 against 5.787 in the last cell is rounding, nothing else).

What *f*<sub>coh</sub> does to the rest of the coherent line is the same
arithmetic on the numbers `tagging_optics.py` publishes at the full year
(rerun 2026-08-29 on the per-configuration levers: N_tag/yr = 2.37 / 2.42 /
6.15 ×10⁶, best-bin 5σ floor 1.74 / 2.34 / 1.62 % per unit P_zz, 3.0 / 5.5 /
2.6 years to 5σ on a 1% exotic-glue term, from ε = 0.374 / 0.251 / 0.332 at
r_h = 46.5 / 164.1 / 89.3 and L/L_HA = 1/6.8 / 1/12.8 / 1/9.5):

| at *f*<sub>coh</sub> = 0.283 / 0.363 / 0.218 | 5 × 40.8 | 10 × 99.5 | 18 × 137.5 |
|---|---|---|---|
| tagged events in the year (× *f*<sub>coh</sub>) | 6.7×10⁵ | 8.8×10⁵ | 1.34×10⁶ |
| best-bin 5σ floor, /√*f*<sub>coh</sub> | 3.27 % | 3.88 % | 3.47 % |
| years to 5σ on a 1% exotic-glue term, /*f*<sub>coh</sub> | 10.6 | 15.2 | 11.9 |

The last row is the reason "what the coherent channel needs" has two
answers and they differ by an order of magnitude.  Read as *5σ on the
deformation (shape) term within the year* — the null test the deuteron
cannot offer, and the thing the sign flip is a prediction about — the
channel needs *f*<sub>coh</sub> = 0.283 / 0.363 / 0.218 and Plan A is
affordable.  Read as *5σ on a 1% exotic-glue term in the best (x, Q²) super-bin alone* — the x-resolved requirement, priced on `tagging_optics.py`'s counting floor over the 19 / 10 / 9% of the tagged sample in that bin — it needs 3.0 / 5.5 /
2.6 **full** years at the full year and 10.6 / 15.2 / 11.9 programme years
at *f*<sub>coh</sub>, so no share of one year buys it and the question
becomes how many years the programme runs, not how one is divided (read instead as the x-integrated flat term of the seven-|t|-bin fit, the δa_e column of the table above, 0.00227 / 0.00184 / 0.00158 at *f*<sub>coh</sub>, gives 4.4 / 5.4 / 6.3σ on 1% per programme year, 5σ within 1.3 / 0.85 / 0.62 programme years).  The
tables here are built on the first reading; the second is the ten-year
column of the reports and is unaffected by any of this arithmetic except
through the count of years.

**Plan B — ⁶Li and ⁷Li each half a year.**  Every share above is halved,
so every error grows by √2 and every year count doubles.  Measured at
`--run-share 0.5`: the inclusive δ(Δ/F₁) is 3.657 / 3.615 / 4.177 ×10⁻⁴
with 3.34 / 3.27 / 4.36 years to 5σ; δA_zz per x-bin at P_zz = 0.8 goes
6.98 / 6.63 / 7.30 ×10⁻⁵ and 2.66 / 8.77 ×10⁻⁴ (x = 0.0035 / 0.0089 /
0.0282 / 0.2818 / 0.5623) to 9.88 / 9.37 / 10.3 ×10⁻⁵ and 3.76 / 12.4
×10⁻⁴ (`money_b1.py`); and the ⁷Li polarized-EMC δΔR goes 0.0477 /
0.0509 / 0.0615 / 0.1224 to 0.0674 / 0.0719 / 0.0870 / 0.1731 at
x = 0.09 / 0.28 / 0.45 / 0.71 (`money_polemc.py --ion 7Li`, toy inputs).

**Plan A × B — ⁶Li takes half a year and splits it, ⁷Li takes the
other half.**  The coherent share becomes ½ *f*<sub>coh</sub> = 0.141 /
0.181 / 0.109, δa_e becomes 0.00322 / 0.00261 / 0.00224 and the shape
term falls to 5/√2 = **3.5σ at every configuration** — by construction,
*f*<sub>coh</sub> having been defined to put it at exactly 5σ — so it
needs two years rather than one; the ⁶Li inclusive share at *s* = 1
is 0.358 / 0.319 / 0.391, which is numerically the *s* = ½ column of
Plan A and gives the same 4.319 / 4.529 / 4.725 ×10⁻⁴ and 4.66 / 5.13 /
5.58 years; ⁷Li is Plan B's column.

Two second-order effects the tables above ignore, in the conservative
direction, and one thing that must not be rescaled at all.

The inclusive measurement does not stop during the coherent share: the
tagging optics de-squeezes the horizontal β* only and inclusive DIS is
blind to the far-forward acceptance, so those stores still deliver
1/6.8, 1/12.8 and 1/9.5 of the high-acceptance luminosity to it.  The
effective inclusive share is (1 − *f*<sub>coh</sub>) +
*f*<sub>coh</sub> L/L<sub>HA</sub> = 0.759 / 0.665 / 0.805 rather than
0.717 / 0.637 / 0.782, which is 2.8 / 2.2 / 1.4 % off every inclusive
error in Plan A.

The longitudinal fills are shared, not additive: one set of
longitudinally polarized ⁶Li stores with electron-helicity flips delivers
A_∥ (and so g₁, the polarized EMC effect) and, in its tensor states, A_zz
(and so b₁).  b₁ and g₁ therefore cost one share between them, not two,
and Plan B's halving is the real cost of the second isotope rather than
of the second observable.

What must **not** be rescaled by 1/√*f*: anything that is not pure
statistics.  The full bars of Report 2 money plot 7R (§5.1) carry the
shape fit, the response Monte Carlo and the unfolding-prior spread, and
that spread does not shrink with luminosity; only their statistical part
obeys the law above, and it has to be separated before the share is
applied.  (The Table 4 bars of the same report are statistical only and
do obey it.)  The same holds for the polarimetry scale (a 1:1
normalisation on Δ) and for the model bands of Report 1 — the Δ/F₁
scenario band of §4 and the a₂ deformation band of §6.3.

### WP3 — Reconstructed-level closure for cos 2φ
- ☑ Port the 20260729 smearing model (tracking σ_p/p, σ_θ, η-dependent
  ε_eID) onto the polligen φ′ pseudo-experiments: migration in (x, Q², φ′),
  reco-binned refit per super-bin.
  ☑ *2026-08-28: done — `polligen/reco.py` carries the 20260729 tables verbatim
  (`tracking_resolution`, `tracking_angular_resolution`, `eps_eid`, plus the
  η-dependent `emcal_resolution`); `recopseudo.RecoResponse` generates on a
  loosened scenario so events migrate in, smears in the lab frame and refits per
  super-bin (`money_cos2phi_reco.py` = money plot 5R). φ′ migration enters as the
  exact second-harmonic dilution cos 2(φ′_reco − φ′_true).*
- ☐ 2-D acceptance-hole closure: (φ, η) hole map → binned-fit bias < stat
  error at 10 fb⁻¹/u (the generator-level version already passes).
  — *superseded (2026-08-28): the WP3 addendum below replaces it with a
  smooth-ε(φ′) closure under the ratio estimator, demonstrated in Report 2
  Figure 1(d) (400 pseudo-experiments under ε = 1 + 0.03 cos 2φ′ + 0.02 cos φ′:
  the single-fill fit biased by 0.03/P_zz, the two-state ratio unbiased and 1.5×
  narrower). No (φ, η) hole map exists or is wanted.*
- ☐ D2: if reco dilution < ~10% and unbiased, keep main figures at generator
  level with reco factors quoted per bin; otherwise switch figures to reco
  level.
- Acceptance: per-bin dilution factors; statement "Case-3-style retention X%"
  reproduced for the cos 2φ observable. Effort: ~1 week.
  ☑ *2026-08-28: met — `recopseudo.RecoResponse.bin_summary` gives purity,
  efficiency and D per super-bin, published in Report 1 Table 1 and Report 2
  Table 4 for both the Gaussian stand-in and the PYTHIA final state; the
  retention statement survives as the reco-over-truth error ratio, 0.60–0.72,
  quoted in both reports and in the WP3-results block below. The literal
  "Case-3-style retention X%" phrasing was never restated, and no fully
  integrated (rate-weighted, all-bins) observable exists — the chain closes per
  super-bin.*

**WP3 addendum (2026-08-24, `reports/reconstruction_chain_report`).** The
reconstructed-level closure is now specified end to end and seeded in
`polligen/reco.py`; the note changes three inputs of this WP:
- ☑ **Binning must be reconstructed (x, Q²) with the mixed (eΣ) method**:
  spots 1–3 of money plot 5 sit at y = 0.010–0.025 where the electron
  alone gives δy/y = 50–120% (Table 1 of the note). Use
  `reco.hadronic_y` with a 15–30% band (purity 0.75–0.83, reco-bin
  amplitude 0.96–0.99 of truth at 15–20%); quote the e′-only variant at
  y ≥ 0.05 and the low-energy configuration for x ≈ 0.1 (open question #21).
  ☐ *2026-08-28: the mixed-method reco binning is done and is the default
  (`reco.mixed_method`/`hadronic_y`, `recopseudo.RecoModel`,
  `money_cos2phi_reco.py --y-method mixed --y-source hfs`). The two trailing
  clauses are not: nothing quotes the e′-only variant at y ≥ 0.05 (the published
  e′-only panels sit at y = 0.010–0.025, where the conclusion is that the bins are
  not reconstructible), and the x ≈ 0.1 case is argued from a δy/y comparison
  rather than run as a low-configuration projection — every 5R/7R number is
  mid-config.*
  ☑ *2026-09-06: the bullet is set to ☑ — the mixed (eΣ) reconstructed binning is
  the default of the reco chain (`reco.mixed_method`/`hadronic_y`,
  `recopseudo.RecoModel`, `money_cos2phi_reco.py --y-method mixed --y-source hfs`).
  The residue stays where it is: the ☐ note above keeps the two trailing clauses
  open — the e′-only variant at y ≥ 0.05 and the x ≈ 0.1 low-configuration
  projection.*
  ☑ *2026-09-15: both trailing clauses are closed, run rather than argued, and
  are Report 2 §5.1 Figures 5 and 6 with their command lines. The x ≈ 0.1 case
  is now a low-configuration projection: at e5 × ⁶Li40.8 (s = 816 GeV²/u) the
  same figure-of-merit selection puts two of the four sweet spots at
  x = 0.089 and 0.141 with Q² = 1.14 and 3.13 GeV², where the middle
  configuration cannot go at all (y = Q²/(sx) = 0.0032 at x = 0.089,
  Q² = 1.14, below the y ≥ 0.01 floor), at purities 0.68–0.70 and
  Â = 13.38 ± 0.11 ×10⁻³ against a reconstructed-bin truth of 13.35 — 0.8% in
  one year (`money_cos2phi_reco.py --config 0`, own stem). The e′-only variant
  is quoted where the method is usable: at y ≥ 0.05 it is the BETTER of the
  two, purity 0.73 / 0.80 / 0.76 / 0.77 against 0.65 / 0.64 / 0.70 / 0.69 for
  the mixed method in the same bins and δÂ 0.98–0.99 of it, unbiased over
  twenty pulls from five response seeds (mean +0.10, sd 0.99); what the window
  costs is reach, x ≤ 0.072 at this s (`--y-method electron --y-min 0.05`, own
  stem). The `--y-min` flag and the filename keys for the beam configuration,
  the y method and the y floor were added with them, and the published 5R/7R
  stems are reproduced md5-identically after the change.*
- ☑ **Estimator: spin-state-sorted ratio** (`reco.harmonic_ratio_fit`) of
  m = ±1-rich (P_zz = +0.6) and m = 0-rich (−1.2) fills, with a sin 2φ′
  term; the single-fill fit is biased by the detector's cos 2φ′ acceptance
  harmonic ÷ P_zz. δA becomes 2√(2/N)/(P₊ − P₀) = 0.67× the current
  values — re-derive the §7.1 numbers with it (the 1.5× gain is real only
  if the source delivers m = 0-rich bunches at the same purity).
  ☐ *2026-08-28: the estimator is done and drives the reconstructed-level
  analysis (`reco.harmonic_ratio_fit` with the sin term, `err_harmonic_ratio`,
  the 0.67 factor derived and published). The second clause is not: §7.1 still
  quotes single-state values, and the three scripts that produce them
  (`money_cos2phi.py`, `money_cos2phi_coherent.py`, `money_delta_extraction.py`)
  have no two-fill path at all.*
  ☑ *2026-09-06: the bullet is set to ☑ — the spin-state-sorted ratio estimator
  is implemented and drives the reconstructed-level analysis
  (`reco.harmonic_ratio_fit` with the sin 2φ′ term, `err_harmonic_ratio`), and the
  0.67 factor is derived and published (commits 8a98f91, bb636b5). The residue
  stays where it is: the ☐ note above keeps the second clause open — §7.1 is still
  quoted from single-state fits and the three generator-level drivers have no
  two-fill path.*
  ☑ *2026-09-15: the residue is now half closed. The three generator-level
  drivers do carry a two-fill path — `--pzz-plus/--pzz-zero`, default off and
  bit-for-bit, on their own `_twofill` stem — and §7.1 is re-derived with it in
  the dated addendum under §7.1 above: δA and δΔ ×0.667 (measured 0.6665–0.6668
  at the four sweet spots), significance 21–43σ → 32–64σ, the amplitudes and the
  coherent model numbers unmoved. What stays open is only which footing the
  published §7.1 and Report 1 quote, which is an author call coupled to the WP1
  grid column (T05).*
- ☑ **Angles from four-vectors**: `reco.azimuth_wrt_lepton_plane` (covariant
  φ_S, = φ_e − φ_S to O(γ²)); head-on transformation applied (e′ odd
  harmonics only). The 2-D φ-hole closure of the original bullet becomes
  a smooth-ε(φ′) closure with the ratio estimator (already demonstrated
  at the super-bin level, Fig. 1d of the note).
  ☑ *2026-08-28: done and in use — the covariant azimuth is built from the
  four-vectors with the Bacchetta transverse projector and the ε-tensor sign
  convention, `recopseudo` forms φ′_true and φ′_reco from them, the head-on ↔ lab
  transformation is applied on both the electron and the hadron side, and the
  smooth-ε(φ′) closure exists with `harmonic_ratio_fit`/`_2d` and the
  `--eff-cos2` split. `lab_azimuth_shortcut_error` quantifies the shortcut it
  replaces.*
- ☑ **Coherent (feeds WP5)**: present the anchored a₂ as a modulation of the
  recoil azimuth φ_t − φ_S and the exotic-glue term of the electron azimuth;
  fit R(α, β) in 2-D; replace the constant 0.20 GeV cut in `coherent.py` by
  the angular cut `reco.tag_pt_cut` (5.0×10⁻⁷ / 8.4×10⁻²⁶ / 3.7×10⁻¹³ at 40.8 / 99.5 /
  137.5 GeV/u) and state the cutout geometry (open question #20); settle
  the a_n normalization convention of arXiv:2408.13213 (1 + 2Σ vs 1 + Σ).
  ☑ *2026-08-28: all five clauses discharged — the two-azimuth presentation and
  the 2-D R(α, β) fit (`reco.basis_2d`, `harmonic_ratio_fit_2d`,
  `recopseudo.CoherentResponse`) are money plot 6R and Report 2 §4.5/Table 5; the
  angular cut (`reco.tag_pt_cut`, `coherent.tag_acceptance_angular`) is what the
  reconstructed chain uses; the cutout geometry is no longer assumed but measured
  through the ePIC geometry (`reco.RP_APERTURE_MEASURED`, plans/04 #20); and
  Eq. (9) settles the convention as 1 + 2Σ. The generator-level money plot 6
  deliberately keeps the constant proton-referenced cut as the reference tag,
  labelled as such in Report 1's Figure 5(a) caption.*
- D2 input: the reco-level dilution is small (1–9% at the 25% default, 1–4% at 15–20%) and unbiased with the
  mixed method, so the main figures can stay at generator level with
  reco factors quoted — provided the binning and the estimator above are
  adopted in the generator-level plots themselves.

**WP3 results (2026-08-24, second pass — `polligen/recopseudo.py`, money
plots 5R/7R/6R, report §7).** ☑ Reconstructed-level pseudo-experiments
exist and close:
- ☑ Inclusive (mixed method, 25% hadronic y — the ePIC kinematic-fit
  study's smearing and ATHENA Fig. 22 at y ≈ 0.01, refs/ 2026-08-25 —
  EMCal E′, track angles, ε_eID, reco cuts, covariant φ′, two-fill ratio
  with a 3%/2% φ′-efficiency harmonic and a 10⁻³ rel-lumi offset on):
  sweet spots 1–4 in reco bins — purity 0.65/0.64/0.66/0.68, efficiency
  0.42/0.60/0.37/0.64, D = 0.907/0.986/0.910/0.969; Â unbiased vs the
  reco-bin truth; **δÂ = 1.2 / 1.0 / 1.8 / 3.0 ×10⁻⁴ (1 yr)** = 0.65–0.70 ×
  the §7.1 baseline (the m = 0-rich fill gain beats the efficiency loss).
  The y = 0.01 edge bins (1, 3) lose ~60% of their events to the reco y
  cut; D = 0.91 there.
- ☑ 7R: δΔ best bins 1.0×10⁻³ (Q² = 1.14) / 0.5×10⁻³ (3.13 GeV²) in year 1
  (vs 0.7–1.4 ×10⁻³ baseline), purities ≈ 0.55.
- ☑ Coherent 6R (feeds WP5): angular envelope 10σ_θ·6p_u (0.22 GeV at
  50 GeV/u) + slot-like cutout |p_x| < 0.55, |p_y| < 0.22 GeV (ePIC pots
  surround a horizontal slot; 2026-08-25 refs pass) → N_tag = 2.7×10⁶
  (1 yr) vs 1.1×10⁷; fake ⟨cos 2β⟩ = +0.77 cancelled by the spin-state
  ratio; template fit (acceptance-weighted MC basis, a_t ∝ |t|) recovers
  a_t(t_ref) in four |t| bins (0.117 ± 0.007 vs 0.119 injected at 1 yr;
  ± 0.002 at 10 yr) and the flat a_e = 0.010 ± 0.0013 (1 yr; ± 0.0004 at
  10 yr). The slot geometry costs ~3× in the deformation-term error vs a
  square cutout (a_e untouched) — a pot-design lever for WP5. u₁/u₂ at
  the ZEUS LPS values (u₁ its central A_LT at 0.01 < x_P < 0.1, u₂ within one error of zero; A_LT, A_TT consistent with zero, NPB 816:1). **Convention verified in
  arXiv:2408.13213 Eq. (9): 1 + 2Σ a_n e^{inΦ} → the deformation
  modulation is 2a₂ = 0.072 at P_zz = 0.6, ⟨|t|⟩ = 0.06 (money plot 6
  injects a₂: conservative ×2).** Update §7.1 / plans/06 / the projection
  report when the WP5 curves are redone.
- D2 decision input: reco-level figures are ready; the recommended
  presentation is reco-level 5R/7R with the ratio estimator (errors improve,
  not degrade), and 6R with the angular cut curves vs σ_θ and p_u.

### WP3-HFS — Hadronic final state and hadron-side detection (☑ chain built, sample generated)
Decision 2026-08-25 (option 1 of the reconstruction-note discussion): replace
the 25% Gaussian stand-in for the hadronic y by a real hadronic final state
through a hadron-side detector response, and treat the hadron-side detection
efficiency explicitly.
- ☑ `polligen/hfs.py`: `HFSSample` (generator-independent .npz format:
  scattered electron + final-state hadrons, head-on frame, per-nucleon
  kinematics), exact hadronic sums, Σ / Jacquet–Blondel / double-angle /
  mixed methods (`hadronic_kinematics`, identities tested against the truth
  with a perfect response), `HadronResponse` (tracker |η| ≤ 3.5, p_T > 0.2 GeV
  with a 95% plateau efficiency and the repository's momentum/angle tables;
  calorimeters |η| ≤ 3.7: photons in the EMCal above 0.1 GeV, neutral and
  untracked charged hadrons in the HCal above 0.5 GeV at 90% efficiency, Yellow
  Report resolution requirements per region; Gaussian noise on Σ and on each
  p_T component, 50 MeV default), `ToyHFS` (vectorized string-fragmentation
  stand-in with exact four-momentum closure and π⁰ → γγ), `HFSLibrary` +
  `HFSResponse` ((x, Q²)-cell library transferring the captured Σ fraction
  and p_T ratio onto the pseudo-events; noise per event).  7 tests.
- ☑ `recopseudo.RecoModel(y_source="hfs", hadronic_method=...)` and
  `RecoResponse(..., hfs=...)`; `money_cos2phi_reco.py --y-source hfs
  [--hfs-sample ...] [--hfs-noise ...]`; `scripts/hfs_resolution.py`
  (Figure 3 of the reconstruction report).
- ☑ `tools/pythia8/gen_dis_hfs.py` + README: PYTHIA 8 e+p / e+n DIS at the
  per-nucleon beam energies (head-on frame, dipole recoil, Q² > 0.7, lepton
  radiation off) → HFS .npz.  **Runs natively — PYTHIA 8.311 builds its own
  Python bindings against the analysis interpreter; the container has the C++
  library but no bindings, which is what made this look container-gated.**
  It also found that `PhaseSpace:Q2Min` had never been applied: PYTHIA
  honours it only for `Q2Min ≥ pTHatMinDiverge²` and that defaults to 1 GeV,
  so the 0.7–1.0 GeV² band — 31% of the sample, and the band the loosened
  generator window exists for — was missing until `pTHatMinDiverge = 0.5`.
- ☑ Samples generated (regenerated 2026-08-27 at the γ-matched energies,
  plans/10): 8 M events, p and n at 5 × 40.8, 10 × 99.5 and 18 × 137.5,
  80.4 M particles, 2.8 GB, manifest in `evgen/samples/README.md`.
- Results **with PYTHIA** (the toy's, kept for contrast, in brackets):
  captured Σ fraction 0.90 (tracks 0.51 [0.60], EMCal 0.28 [0.30], neutral
  hadrons 0.11 [0.03]) *(superseded 2026-08-27 — pre-correction energies and
  spots; at the corrected mid spots 69–85% is captured through the response
  and HCal objects — neutral hadrons plus the charged particles the tracker
  does not see — carry 0.09–0.10 of Σ within acceptance, Report 2 §3 Figure 2
  and Table 1b)*; Σ-method δy/y at the four sweet spots with 50 MeV
  noise **0.55 / 0.28 / 0.50 / 0.15** [0.28 / 0.17 / 0.24 / 0.07], i.e. the
  toy was optimistic by 0.04–0.05 absolute everywhere; 0.21 / 0.12 / 0.17 /
  0.10 at LOW and 0.74 / 0.34 / 0.69 / 0.18 at TOP, so **the x ≈ 0.1 bins
  belong to the low-energy configuration**.  5R with the HFS-based y:
  purities 0.43 / 0.54 / 0.47 / 0.69 [0.60 / 0.68 / 0.68 / 0.86], amplitude
  dilution D = 0.79 / 0.85 / 0.82 / 0.95, δÂ = 1.2 / 1.0 / 1.8 / 2.9 × 10⁻⁴
  (Table 2 errors unchanged — the loss is in purity, not in statistics).
- ☑ Quote the resolution table by y and Q² in the reconstruction report and
  decide whether 5R/7R are published on the HFS y or the stand-in.
  ☑ *2026-08-28: Report 2 Table 3 is that table — the four sweet spots with y, W,
  e′ kinematics, δy/y for the electron alone and for Σ/JB/DA, with the acceptance
  and threshold columns, backed by Figure 3(a), (c) and the 0/25/50/100 MeV noise
  scan of 3(d) (`hfs_resolution.py`, `hfs_acceptance.py`). The decision resolved as
  publish both: Table 4 carries the stand-in and the calibrated PYTHIA columns side
  by side and §5.1 quotes δΔ both ways, with the figures made on the PYTHIA final
  state.*
- ☐ Replace the Yellow-Report response magnitudes by the ePIC design values
  (calorimeter noise/threshold floor at Σ_h ≈ 0.2–0.5 GeV is the decisive
  input — plans/04 #21 narrowed to it); add the HFS energy-scale
  calibration and the Σ-method ISR test once WP4 exists.
- Effort: sample generation ≈ 1 h machine time; the rerun and the report
  update 1 day.

### WP4 — Radiative-correction bound (not a calculation) — ☑ **closed 2026-08-28**
- ☑ Leading-log unpolarized RC weights (plans/02 step 1.4 route) applied as
  kinematic migration. `evgen/polligen/radiative.py`: the exponentiated
  leading-log electron structure function D(z, Q²) = (t/2)z^(t/2−1)S(t) −
  (t/4)(2−z), t = (2α/π)[ln(Q²/mₑ²) − 1] (Kuraev–Fadin / Nicrosini–Trentadue;
  ∫D = 1 + O(t²), residual 7×10⁻⁴ at t = 0.070), a per-event sampler with its
  own random stream, the closed-form observed kinematics of all five
  reconstruction methods, the (1 − z) rescalings of the hard rate and of a₂,
  and `migration_bound` / `migration_bound_seeds`. Hook:
  `recopseudo.RecoResponse(isr=…)`, **default off and bit-for-bit inert**
  (pinned against a stored digest of the response arrays, against the closed
  form of the generator weight, and against the state of the response's own
  random stream); driver `money_cos2phi_reco.py --isr [--isr-seeds]
  [--isr-gen-q2min] [--isr-empz]`. 30 tests.
  - **The fake-modulation term is identically zero.** The covariant azimuth is
    invariant under k → (1 − z)k for a massless target: cos φ′ and sin φ′ carry
    the same factor [2ac((1−z)a−b)]^(−1/2) and the arctan divides it out.
    Measured **3.6×10⁻¹⁵ rad** over the 2×10⁴-event flat sample of
    `test_covariant_azimuth_is_invariant_under_a_collinear_photon` (z ≤ 0.9);
    over the 1.84×10⁶ events of the response, where the physical ⁶Li mass
    leaves the O(γ²) residual, max |Δφ′| = 2.6×10⁻² rad and the fake cos 2φ′ is
    **9×10⁻⁸** rate-weighted (`RecoResponse.isr_dphi`, printed by `--isr`).
    The two samples must not be conflated.
  - **x is exact, the Q²ₑ label migrates by 1/(1 − z)** (plans/08 D3, code review
    R16), pinned against a four-vector construction through `hfs.hadronic_kinematics`.
- ☑ Spin-state-ratio cancellation argument written up **with the bound**. The
  cancellation, the non-cancelling residual (ε₂⁺ − ε₂⁰)/(P₊ − P₀) = 5.6×10⁻⁴ and
  the bunch-by-bunch requirement are in Report 1 §3.3 / assumption row 8, Report 2
  §4.3/§6 and plans/06; the bound they were missing is now Report 2 §7 and its
  Table 2 row. An unpolarized-lepton QED correction is common to the fills by
  construction, so the whole rate effect (−0.3% to +3.7% per bin) cancels and only
  the amplitude migration survives.
- ☑ Gate: **PASS** → one paragraph (Report 2 §7) + assumptions row (Report 2
  Table 2); no appendix needed. Mid configuration, four sweet spots, 1600
  pseudo-events per cell, ISR seed 20260828, common random numbers,
  Δ̂ = Â × K(ISR-free). One response draw scatters by 4–14% of the bound
  (seed-to-seed sd 0.087 / 0.048 / 0.073 / 0.051 points; the eight draws span
  0.51–0.75, 0.43–0.57, 0.80–1.01 and 1.15–1.31%), so **every number below is the
  mean ± sem over the eight response seeds** of `--isr-seeds
  20260824,20260925,20261026,20261127,20261228,20270129,20270302,20270403`
  (the plain `--isr` at the default seed 20260824 prints one draw of the same
  quantity: +0.62 / +0.50 / +0.80 / +1.24%):
  - published generator window (Q² ≥ 0.7 GeV²): **+0.62 ± 0.03 / +0.50 ± 0.02 /
    +0.94 ± 0.03 / +1.22 ± 0.02%**; purity 0.653 → 0.638, 0.633 → 0.613,
    0.679 → 0.659, 0.684 → 0.640; efficiency 0.414 → 0.404, 0.590 → 0.572,
    0.374 → 0.369, 0.653 → 0.634; selected rate ×0.997, ×0.998, ×1.017, ×1.037.
  - the window truncates the feed-in (an event below it cannot radiate into an
    analysis bin). Opening it, the worst spot rises **1.22 → 1.87 → 2.81 →
    2.26 → 2.34%** at Q²_gen = 0.7 → 0.35 → 0.15 → 0.05 → 0.02 GeV², i.e. it
    saturates in a 1.8–2.8% band below Q²_gen ≈ 0.15 rather than at one value.
    At Q²_gen = 0.05 the four spots are **+2.26 ± 0.03 / +2.24 ± 0.10 /
    +1.88 ± 0.07 / +1.88 ± 0.05%**; the largest value found at any window is
    +2.81 ± 0.08% (Q²_gen = 0.15, second spot). **≤ 2.9% is the number the ≤5%
    gate is read against**, and it passes by a factor 1.7.
  - not an artefact of the 25% Gaussian y stand-in: through the PYTHIA hadronic
    final state with the calibrated scale (`--y-source hfs --hfs-sample …
    --hfs-calibrate`, 800/cell) the same bound is **+0.38 / +0.44 / +0.46 /
    +0.88%** (± 0.02–0.04).
  - **mitigation the chain does not use:** the HERA E − p_z window. The
    visible sum is 2(1 − z)E_e and is already reconstructed
    as Σ_h + E′(1 − cos θ) = E′(1 − cos θ)/(1 − y_Σ); requiring it within 15% of
    2E_e brings the bias to **+0.23 / +0.16 / +0.22 / +0.18%** on the 25%
    Gaussian y stand-in and **+0.17 / +0.19 / +0.20 / +0.16%** through the PYTHIA
    hadronic final state with the calibrated Σ scale (± 0.01–0.03 on both),
    independent of the generator window, and keeps 0.869 / 0.824 (stand-in) and
    0.979 / 0.929 (PYTHIA, calibrated Σ) of the non-radiative / radiative
    selected rate; the loss is almost entirely above y = 0.2 — 99.5% and 97.2%
    of the discarded non-radiative rate — and costs 0.01–0.06% at the four
    sweet spots. The window is a documented contingency, not a default: apply it
    if a published analysis opens the generator window below Q² = 0.15 GeV², or
    if the 5% gate tightens. It is not the default because the gate already
    passes at ≤ 2.9% without it, because it is free where the letter's numbers
    live, and because what it does remove is the y > 0.2 rate that carries the
    low-x end of every Q² slice in money plots 5R and 7R
    (`radiative.empz_fraction`, `apply_empz_cut`,
    `empz_bin_retention`, `empz_y_retention`).
  - **method comparison** at z = 0.092 (the mean radiated fraction of the emitting
    events), observed/hard for (Q², y, x). The electron-method rows are a strong
    function of y — (y + z)(1 − z)/y — so the y has to be stated. At the
    rate-weighted ⟨y⟩ = 0.189 of the whole selected sample: electron
    (1.102, 1.351, 0.740), Σ (1.000, 1.000, 0.908), Jacquet–Blondel
    (0.976, 0.908, 0.976), double angle (1.214, 1.000, 1.102),
    **mixed (1.102, 1.000, 1.000)**. At the four sweet spots themselves
    (y = 0.0101, 0.0254, 0.0111, 0.0254) the electron method is far worse —
    y is off by 9.2, 4.2, 8.5, 4.2 and x by 0.109, 0.239, 0.118, 0.238 — while
    every other row is unchanged except Jacquet–Blondel's Q² (0.998). That
    factor 3–7 between the two y is why the chain uses the mixed method.
    Q²_Σ = p_T,e²/(1 − y_Σ) is ISR-exact too, so an e-Σ *label* would have no
    migration at all — a chain change, not made.
- **Still open and outside this bound:** the TENSOR-sector radiative correction
  (plans/05 §5.5) — never calculated, and no unpolarized study stands in for it;
  the polarized-lepton correction (irrelevant, unpolarized beam); wide-angle real
  emission; FSR; the elastic and quasi-elastic radiative tails (removed by
  W² ≥ 10 GeV²).
- Effort spent: ~1 day (vs 3–5 estimated; the DJANGOH route of plans/02 step 1.4
  was not needed).

### WP5 — Coherent channel presentation
**◐ 2026-08-25: the scan exists** — `evgen/scripts/coherent_optics_scan.py`
(plans/08 A4) gives all four panels: analytic acceptance vs the envelope
over B = 40–60 for the slot / square / circular cutouts with the three
beam configurations marked; tagged yield with exp(−B t_min) (×0.85) and
the 0.73 rate weighting folded in, against the IR-8 published band; the
fitted δa_t and δa_e from the full response (importance-sampled above the
cut — the plain sampler leaves *zero* accepted recoils above 0.3 GeV);
and acceptance vs beam momentum. Numbers: tagged fraction
32% / 3.0% / 4×10⁻⁵ / 2×10⁻⁷ and δa_t/a_t = 1.2% / 4.6% / 79% / 392% at
an envelope of 0.10 / 0.22 / 0.45 / 0.60 GeV. **The coherent measurement
lives at the low- and mid-energy configurations and is dead at the top
energy.** Remaining: fold the curve into the letter figure and quote the
cutout geometry as the assumption it is (#20).
- ☑ Replace two-point optics with curves vs pT_cut (0.1–0.7 GeV):
  acceptance, N_tag, best-bin δA; mark 0.20 (documented top-rigidity scale,
  ³He precedent) and 0.45 (our derivation) on the curves; state that Li
  optics are undocumented and the physics case constrains them.
  ☑ *2026-08-28: done in `evgen/scripts/coherent_optics_scan.py` — acceptance,
  N_tag and the fitted δa_t/a_t, δa_e/a_e over a 0.05–0.70 GeV cut scan for the
  slot / square / circular cutouts across B = 40–60, with the undocumented-optics
  statement in the module docstring and Report 1 §6.1. The "mark 0.20 and 0.45"
  clause is half superseded: 0.45 came from the legacy 164 µrad high-divergence
  optics retired in b9d2e82, and the script now marks the per-configuration
  Yellow Report envelopes, the tagging optics and the measured pot aperture
  instead.*
- ☐ Fold exp(−B t_min) and the ×0.73 rate-weighting option into the central
  tagged-yield curve (small changes in `polligen/coherent.py` + tests).
  ☐ *2026-08-28: the figure half is done (panel (b) draws N_tag, ×exp(−B t_min)
  and ×0.73), but the library half is not: `t_min` is hard-coded in the script —
  with its own slope, so the suppression curve cannot follow the B band the panel
  draws — `coherent.tag_acceptance`/`mean_t_tagged`/`tag_acceptance_angular` take
  no t_min, `a2_tagged` still asks the caller to apply `RATE_WEIGHT_SYST` by hand,
  and no test mentions it. Money plot 6 still quotes the unweighted value.*
  ◐ *2026-09-15: the library half is done and the fold is OFF by default —
  whether it becomes the central curve is an open author call, raised to the
  supervisor in run 19 with the numbers below. `coherent.py`
  gains `t_min_coherent(x_P, M_A) = (x_P M_A)²/(1 − x_P)` and
  `CoherentScenario.t_min_suppression(t_min) = exp(−B t_min)`, and a `t_min=None`
  keyword on `tag_acceptance`, `mean_t_tagged`, `tag_acceptance_angular`,
  `a2_tagged`, `recoil_lab` and `project_coherent`; `a2_tagged` gains
  `rate_weighted=False`, which applies `RATE_WEIGHT_SYST` in place of asking the
  caller to multiply by hand. Because the suppression is now a method of the
  scenario it follows the B band: at |t_min| = 3.1×10⁻³ GeV² it is −11.7% /
  −14.4% / −17.0% at B = 40 / 50 / 60, where the script's hard-coded constant
  gave one number. Every default is None/False, so money plot 6
  (`money_cos2phi_coherent_6Li.png`, md5 aae59967…), the WP5 scan
  (`coherent_optics_scan_6Li.png`, md5 e30729d5…) and `phase_space_bins_6Li.png`
  regenerate bit-for-bit against the committed PNGs. The 2026-08-10 audit comment
  in `recoil_lab` is now eight tests in `test_coherent.py` rather than prose
  (evgen 334 → 342): its −14% is measured at −14.36% (and −14.66% with the
  1/(1 − x_P) the note's own formula carries but its number dropped — that 1% of
  t_min is the whole gap between this note's −14% and the 2026-08-25 code
  review's −15%). Its second number does not survive as quoted: "rate-weighted
  over f_coh ≈ 10%" is a statement about the x_P weight, not a number — f_coh with
  x_P flat over 0 < x_P < 0.02 gives −11.0%, but the DIS-like dN/dx_P ~ 1/x_P
  gives only −2.5%. What is robust is only the one-sidedness. Review added the
  missing caveat: x_P = 0.01 is `x_coh`, the coherence half-point, not the top
  of the window — |t_min| ∝ x_P², so at the x_P ≈ 0.02 edge this file's own
  module docstring quotes, the B = 50 suppression is **−47%**, three times
  deeper. Both measured weights land above the 0.01 value only because f_coh
  kills the large-x_P end, so 0.01 is a reference point, not a bound, and the
  docstrings of `t_min_coherent` and `project_coherent` now say so. Still open:
  the script-side
  `t_min = 3.2e-3` constant of `coherent_optics_scan.py:214` should call
  `t_min_coherent(0.01)` (= 3.169×10⁻³) and `t_min_suppression` so panel (b)'s
  suppression curve follows the band — not done here, the figure is bit-for-bit
  and the re-quote belongs with money plot 6.*
- ☑ IR-8 panel/inset: published efficiencies d 47% / ³He 32% / ⁴He 29% /
  ⁷Li 17.8% (no ⁶Li — interpolation labeled ours), pT ≈ 0 reach.
  ☑ *2026-08-28: delivered as an overlay on panel (b) rather than a separate
  panel — `IR8_PUBLISHED` with `IR8_LI6_INTERPOLATED = 0.20` drawn as a
  cut-independent line (the p_T ≈ 0 reach) over an axhspan of the published
  min/max, labelled "interpolated, ours", and repeated in Report 1 §6.4 and
  plans/06 §6.5. The annotation rounds ⁷Li to 18%; the exact 17.8% lives in
  Report 1 §6.1 and refs/README.md.*
- ☑ Geometry note: quote current RP z (32.5/34.3 m scan) alongside the
  YR-era 26/28 m, windows in θ/R unchanged.
  ☑ *2026-08-28: done in `farforward.py`'s module header (z = 32.55/34.25 m read
  from current `eic/epic` main, θ < 5 mrad, R ∈ [0.60, 0.95]) and in Report 3
  Table 7, whose caption makes the alongside-quote explicit; `tools/fullsim/README.md`
  carries the reason for the move. Residual elsewhere: plans/03 §2.2 still lists
  26/28 m and 22.5/24.5 m.*
- Effort: 2–3 days.

### WP6 — Paper production
- ☐ `paper/` directory: `main.tex` (elsarticle, two-column), `refs.bib`
  (INSPIRE keys for the verified list), `figs/`, build script.
- ☐ Letter-quality figure pass: a shared matplotlib style (column widths,
  8–9 pt fonts, consistent Okabe–Ito accents); condense money plot 5 →
  Fig. 2 (two φ′ panels + amplitude vs x), money plot 7 → Fig. 3 (two Q²
  slices), money plot 6 → Fig. 4 (a₂ anchor/band + tagged φ′, acceptance
  inset); phase space → Fig. 1 (single inclusive panel with bins, coherent
  support as contour or inset).
- ☐ Text: draft from the fact-checked report + `docs/note_cos2phi_coherent_6Li.md`;
  budget below. Cover letter: first-projection claims + the two verified
  literature gaps.
- Effort: ~1 week drafting + iteration.

### WP7 — Non-blocking parallel items
- ☐ FLUKA license → local BeAGLE incoherent ⁶Li shapes (upgrade the |t|-fit
  purity band if it lands in time).
- ☐ Theory contacts on record: IP-Glasma α–d ask (#18); Cloët convention
  (#6); tensor RC (#10). None block submission; each strengthens revision.

## 7.5 Letter skeleton (word budget ≈ 4,200)

| § | content | words | assets |
|---|---|---|---|
| 1 | Introduction: exotic glue, never measured, the two verified gaps, why tensor-polarized ⁶Li at the EIC | 600 | — |
| 2 | Observable: master formula, amplitude, Δ–δG relation, moment estimates (Δ⁺⁺ bag; lattice φ/deuteron with caveats) | 450 | typeset eqs from report |
| 3 | Simulation framework: beams (EPIOS), polligen validation (one paragraph), grid SFs, acceptance, statistics method | 500 | Fig. 1 |
| 4 | Inclusive projections: φ′ pseudo-data, amplitude vs x, xΔ extraction, A-vs-B discrimination | 700 | Figs. 2–3, Table 1 |
| 5 | Coherent channel: detection at IP6/IR-8, rate model bands, anchored a₂, sign flip, null test, two-component fit | 700 | Fig. 4 |
| 6 | Systematics and assumptions: reco dilution, RC bound (collinear-ISR migration ≤ 2.9% of Δ̂, ≤ 0.25% behind an E − p_z window (a documented contingency, not applied); tensor-sector RC still uncalculated), polarimetry, purity via |t| fit, acceptance stability between the spin-state samples (bunch-by-bunch requirement: a 10⁻³ difference of the cos 2φ′ efficiency harmonic fakes half the signal), K model dependence (3–11% between Δ shapes with the bin-by-bin factor, ≤ 1.2% with the folded shape fit), O(γ²) b₁ leakage, measured since 2026-08-29 rather than bounded (`evgen/scripts/tensor_gamma_leakage.py` on the exact finite-γ tensor kernel, plans/08 D2): Δ_fake = (0.14–0.16) γ²b₁ — the full Cosyn combination, whose leading-twist and twist-3 channels stand as T_LL : T_LT : T_TT = 3 : −3 : 1 and cancel, leaving the twist-4 Eq. (17e) term almost alone — worth **0.109%** of the cos 2φ amplitude at the worst of the twelve sweet spots (5 × 41, x = 0.089, Q² = 1.14), ≤ 0.033% at 10 × 99.5 and ≤ 0.027% at 18 × 137.5; negative at every spot in the adopted tensor convention, so it cancels part of the amplitude of the moment-constrained models rather than faking one, and subtractable in situ from the same transverse fills because 99.96% of it is proportional to b₂ = 2x b₁ at leading twist — the combination the fit's own constant κ̂ measures, once the bin-independent luminosity pedestal on κ̂ is calibrated out — with no longitudinal fill required (b₃ = b₄ = 0.1 b₂, which breaks the cancellation, would move it to 0.175% and is the 0.6 of the leakage the subtraction cannot reach) — R model in Δ/F₁ = −2(1+R)Â, now the published R1998 and worth +16.6 / +18.0 / +4.7 / −4.4% at the sweet spots — code review 2026-08-25 — condensed assumptions | 450 | — |
| 7 | Summary and outlook (theory and machine asks) | 250 | — |

Title candidates (D5): (a) "Nuclear gluonometry with a tensor-polarized ⁶Li
beam at the Electron-Ion Collider"; (b) "Projections for the
double-helicity-flip structure function of ⁶Li at the EIC"; (c) variant of
(a) with "exotic glue" in the subtitle. Recommend (a).

## 7.6 Referee risk register

| objection | response | where |
|---|---|---|
| "The Δ model is arbitrary" | moment-constrained (Δ⁺⁺ bag provenance stated); A/B interpretations bracket; the measurement discriminates them — that is the point | §2, §4 |
| "Bag moment ported to a nucleus" | labeled scenario; literature brackets both directions (NPLQCD suppression vs binding enhancement); x-shape ours by necessity (S–S give none) | §2, §6 |
| "The b₁ sector contaminates the cos 2φ amplitude at finite γ" | measured, not bounded: the exact Cosyn kernel (b₃/b₄ slots, Eqs. 9/10/14/17/24), anchored on both finite-γ rows of that paper's own Table 1, puts it at 0.109% of the amplitude at the worst of the twelve sweet spots and ≤ 0.033% at the published configuration, with the sign that cancels rather than fakes, and 99.96% of it is proportional to b₂ = 2x b₁ at leading twist — the combination the fit's own constant κ̂ measures, so it is subtracted in situ from the same transverse fills, with no longitudinal fill required, once the bin-independent luminosity pedestal on κ̂ is calibrated out; the b₃/b₄ band leaves 0.6 of it | §6, plans/08 D2 |
| "No tensor radiative corrections" | quantified migration bound (≤ 2.9% of Δ̂ with the low-Q² feed-in opened up, ≤ 0.25% with an E − p_z window (a documented contingency, not applied); Report 2 §7); collinear ISR fakes no cos 2φ′ at all and cancels in the spin-state ratio; open-theory statement with citation trail | §6 (WP4) |
| "Li beams do not exist; no luminosity" | EPIOS PRC 113:060501 feasibility; stated 10 fb⁻¹/u with linear scaling; P_zz band quoted | §3, Table 1 |
| "Coherent fraction is invented" | explicit f₀ band bracketing HERA ep and heavy-A saturation; first-of-kind labeled; IP-Glasma ask on record | §5 |
| "pT cut undocumented for Li" | curves vs cut, not a point estimate; documented anchors marked; IR-8 alternative with published numbers | §5 (WP5) |
| "Generator-level only" | reco-level dilution factors from ePIC-parameter smearing; φ-hole closure demonstrated | §6 (WP3) |
| "α+d background fakes the tag" | m-state-blind → dilutes, cannot fake; |t|-shape purity with e+Pb benchmark; Z-ID question stated as open | §5–6 |

## 7.7 Authorship, circulation, timeline (D3, D4)

- **D3 authorship** (user's call): lead C. Peng, second author J. Zhou (ANL
  Physics Division; *2026-08-28: Reports 0–4 now carry "C. Peng and
  J. Zhou" with a "Writing assisted by Claude (Anthropic)" line under the
  affiliation — the letter keeps both*); candidates to invite —
  I. Cloët (⁶Li structure/conventions), J. Maxwell (gluonometry lineage,
  LOI), EPIOS accelerator contact (one, for the beam paragraph), optionally
  W. Chang/A. Jentsch (far-forward blessing). A short-author-list projection
  letter is also viable (Friščić precedent had six).
- **D4 circulation**: 2–3 week comment window on a frozen v1 before
  submission; the INT organizers' orbit is the natural reviewer pool.

| date (2026–27) | milestone |
|---|---|
| Sep 5 | WP1 + WP2 done; §7.1 numbers re-derived on grids (addendum here) |
| Sep 26 | WP3 (reco closure) + WP5 (coherent curves) done; D2 decided |
| Oct 10 | WP4 bound done; gate passed or appendix planned — ☑ **done 2026-08-28, gate passed, no appendix** |
| Oct 31 | WP6: full draft v1 + letter figures |
| Nov | circulation (D3 list), revisions |
| Dec 19 | v2 frozen; co-author sign-off |
| Jan 2027 | arXiv + PLB submission; Zenodo DOI |
| Mar 22 | INT program talk with the paper on arXiv (referee reports likely in hand) |

Total new effort ≈ 4–5 working weeks spread over ~4 months; the writing
draws on already fact-checked text.

## 7.8 Definition of done

- ☐ All §7.1 numbers grid-based, reco-informed, with P_zz/dilution bands.
- ☐ Four letter figures + Table 1 regenerated by one `paper/` build script.
- ☐ Every citation BibTeX-verified against INSPIRE; the three "do not say"
  lists from the 2026-08-17 fact checks enforced in the text.
- ☐ Assumptions table of the report reduced to a §6 paragraph with no claim
  the fact checks flagged as unsourceable.
- ☐ Cover letter states the two verified literature firsts.
- ☐ Submitted; arXiv number recorded here.
