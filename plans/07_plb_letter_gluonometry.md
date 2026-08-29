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
| 3 | dilution convention (plans/04 #6) | 1/3 baseline (Cloët convention); cluster picture gives ≈0.81 (×2.4) | keep 1/3 as conservative baseline, quote the cluster value as labeled upside; ask Cloët (D3) | WP2 | yes (one paragraph) |
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
- ☐ Promote the money_delta script-local `NuclearF2FromGrid`
  (EPPS21nlo_CT18Anlo_Li6) into `polli_fastsim/structure.py`
  behind the existing `NuclearF2` interface; wire through
  `inputs.get_backends` (pass `r_func=structure.r1998` there, not a
  monkey-patch).
- ☐ Rerun money plots 5/6/7 and `phase_space_bins.py` with `--backend grid`;
  re-solve moment_A at the grid ⟨Q²⟩; record sweet-spot drift and the
  toy-vs-grid ratio per headline number (expect ≲ ×1.5 by run-2/3 checks:
  toy-vs-CT18 F₂ within ±37%).
- Acceptance: all §7.1 numbers re-derived on grids; drift table in the plan
  addendum. Effort: 2–4 days.

### WP2 — Polarization, run-plan, and dilution bands
- ◐ δA scaling table vs P_zz and fill share (analytic 1/(P_zz√(fN)) checked
  against one sampler rerun).
  ◐ *2026-08-28: the fill-share half is done and priced below — the share is a
  flag (`--run-share` / `--lumi-fraction`, default 1.0), the 1/√f law is verified
  against direct reruns rather than assumed, and Plans A, B and A×B are tabulated.
  The P_zz half of the same table is still to write.*
- ☐ Dilution paragraph: 1/3 baseline; 0.81-cluster upside as the labeled
  alternative (amplitudes ×2.4); one email to I. Cloët (D3) — do not block.
  ☐ *2026-08-28: the paragraph is written (Report 1 §3.2 and assumption row 3)
  and 1/3 is the code default, with the ×2.4 exact because the bag-moment
  normalization is solved before the dilution multiplies Δ. The email is unsent —
  external, plans/04 #6.*
- Acceptance: Table 1 of the letter exists. Effort: 1 day.

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
*f*<sub>coh</sub> = (5/9.3)², (5/8.3)², (5/10.7)² = 0.289 / 0.363 / 0.218
from the 9.3 / 8.3 / 10.7 σ per year of `tagging_optics.py` at the
tagging optics.  The rest runs at the Yellow Report high-acceptance
optics and is shared *s* = 1, ½, ⅓ among the three inclusive ion fills.

| configuration | *f*<sub>coh</sub> | δa_e coherent, 1 yr at *f*<sub>coh</sub> | inclusive share, *s* = 1 | δ(Δ/F₁) | yr to 5σ | *s* = ½ | δ(Δ/F₁) | yr to 5σ | *s* = ⅓ | δ(Δ/F₁) | yr to 5σ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 × 40.8 | 0.289 | 0.00221 | 0.711 | 3.067×10⁻⁴ | 2.35 | 0.355 | 4.337×10⁻⁴ | 4.70 | 0.237 | 5.312×10⁻⁴ | 7.05 |
| 10 × 99.5 | 0.363 | 0.00173 | 0.637 | 3.202×10⁻⁴ | 2.56 | 0.319 | 4.529×10⁻⁴ | 5.13 | 0.212 | 5.546×10⁻⁴ | 7.69 |
| 18 × 137.5 | 0.218 | 0.00158 | 0.782 | 3.341×10⁻⁴ | 2.79 | 0.391 | 4.725×10⁻⁴ | 5.58 | 0.261 | 5.786×10⁻⁴ | 8.37 |

δa_e is the combined one-year error on the exotic-glue coefficient over
the seven |t| bins, 0.00119 / 0.00104 / 0.00074 at the full year
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
(rerun 2026-08-28: N_tag/yr = 2.59 / 3.01 / 6.15 ×10⁶, best-bin 5σ floor
1.67 / 2.10 / 1.62 % per unit P_zz, 2.8 / 4.4 / 2.6 years to 5σ on a 1%
exotic-glue term, from ε = 0.423 / 0.323 / 0.332 at r_h = 49.7 / 175.6 /
89.3 and L/L_HA = 1/7.1 / 1/13.3 / 1/9.5):

| at *f*<sub>coh</sub> = 0.289 / 0.363 / 0.218 | 5 × 40.8 | 10 × 99.5 | 18 × 137.5 |
|---|---|---|---|
| tagged events in the year (× *f*<sub>coh</sub>) | 7.5×10⁵ | 1.09×10⁶ | 1.34×10⁶ |
| best-bin 5σ floor, /√*f*<sub>coh</sub> | 3.11 % | 3.49 % | 3.47 % |
| years to 5σ on a 1% exotic-glue term, /*f*<sub>coh</sub> | 9.7 | 12.1 | 11.9 |

The last row is the reason "what the coherent channel needs" has two
answers and they differ by an order of magnitude.  Read as *5σ on the
deformation (shape) term within the year* — the null test the deuteron
cannot offer, and the thing the sign flip is a prediction about — the
channel needs *f*<sub>coh</sub> = 0.289 / 0.363 / 0.218 and Plan A is
affordable.  Read as *5σ on a 1% exotic-glue term* it needs 2.8 / 4.4 /
2.6 **full** years at the full year and 9.7 / 12.1 / 11.9 programme years
at *f*<sub>coh</sub>, so no share of one year buys it and the question
becomes how many years the programme runs, not how one is divided.  The
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
0.0511 / 0.0619 / 0.1237 to 0.0675 / 0.0722 / 0.0876 / 0.1749 at
x = 0.09 / 0.28 / 0.45 / 0.71 (`money_polemc.py --ion 7Li`, toy inputs).

**Plan A × B — ⁶Li takes half a year and splits it, ⁷Li takes the
other half.**  The coherent share becomes ½ *f*<sub>coh</sub> = 0.145 /
0.181 / 0.109, δa_e becomes 0.00313 / 0.00244 / 0.00224 and the shape
term falls to 5/√2 = **3.5σ at every configuration** — by construction,
*f*<sub>coh</sub> having been defined to put it at exactly 5σ — so it
needs two years rather than one; the ⁶Li inclusive share at *s* = 1
is 0.355 / 0.319 / 0.391, which is numerically the *s* = ½ column of
Plan A and gives the same 4.337 / 4.529 / 4.725 ×10⁻⁴ and 4.70 / 5.13 /
5.58 years; ⁷Li is Plan B's column.

Two second-order effects the tables above ignore, in the conservative
direction, and one thing that must not be rescaled at all.

The inclusive measurement does not stop during the coherent share: the
tagging optics de-squeezes the horizontal β* only and inclusive DIS is
blind to the far-forward acceptance, so those stores still deliver
1/7.1, 1/13.3 and 1/9.5 of the high-acceptance luminosity to it.  The
effective inclusive share is (1 − *f*<sub>coh</sub>) +
*f*<sub>coh</sub> L/L<sub>HA</sub> = 0.752 / 0.664 / 0.805 rather than
0.711 / 0.637 / 0.782, which is 2.7 / 2.1 / 1.4 % off every inclusive
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
  retention statement survives as the reco-over-truth error ratio, 0.59–0.69,
  quoted in both reports and in the WP3-results block below. The literal
  "Case-3-style retention X%" phrasing was never restated, and no fully
  integrated (rate-weighted, all-bins) observable exists — the chain closes per
  super-bin.*

**WP3 addendum (2026-08-24, `reports/reconstruction_chain_report`).** The
reconstructed-level closure is now specified end to end and seeded in
`polligen/reco.py`; the note changes three inputs of this WP:
- ☐ **Binning must be reconstructed (x, Q²) with the mixed (eΣ) method**:
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
- ☐ **Estimator: spin-state-sorted ratio** (`reco.harmonic_ratio_fit`) of
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
  the ZEUS LPS 1σ bounds (A_LT, A_TT consistent with zero, NPB 816:1). **Convention verified in
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
    (y = 0.0102, 0.0256, 0.0111, 0.0255) the electron method is far worse —
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
| 6 | Systematics and assumptions: reco dilution, RC bound (collinear-ISR migration ≤ 2.9% of Δ̂, ≤ 0.25% behind an E − p_z window (a documented contingency, not applied); tensor-sector RC still uncalculated), polarimetry, purity via |t| fit, acceptance stability between the spin-state samples (bunch-by-bunch requirement: a 10⁻³ difference of the cos 2φ′ efficiency harmonic fakes half the signal), K model dependence (3–11% between Δ shapes with the bin-by-bin factor, ≤ 1.2% with the folded shape fit), O(γ²) b₁ leakage — **not** γ²b₁/6: the full Cosyn Eq. (17e) combination is ≈ 6.9× that, so the bound is γ²b₁ × 1.15, still ≤ 0.15% of everything published and exposed only at Δ/F₁ ~ 10⁻³ and x ≳ 0.2, Q² ≈ 1 (plans/08 D2) — R model in Δ/F₁ = −2(1+R)Â, now the published R1998 and worth +16.6 / +18.0 / +4.7 / −4.4% at the sweet spots — code review 2026-08-25 — condensed assumptions | 450 | — |
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
