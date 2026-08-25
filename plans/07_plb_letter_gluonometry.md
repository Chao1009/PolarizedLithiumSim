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
| 5 | radiative corrections | uncalculated for tensor observables (flagged) | quantified migration bound + self-normalization argument + explicit open-theory statement | WP4 | yes (as a bound, not a calculation) |
| 6 | coherent optics cut | two-point (0.20/0.45 GeV); Li optics undocumented (report assumption #1) | acceptance, N_tag, δA as curves vs pT_cut ∈ [0.1, 0.7] GeV; IR-8 alternative with published efficiencies; t_min and ×0.73 rate-weighting folded into the central curve | WP5 | **yes** |
| 7 | far-forward geometry | RP z = 26/28 m (YR-era) in text; run-2 Geant4 scan found 32.5/34.3 m in epic-main | quote current ePIC geometry, windows unchanged in θ/R | WP5 | minor |
| 8 | coherent fraction f₀ | scenario 0.04 [×2/÷2] | unchanged (no light-A prediction exists — stated); theory ask on record (IP-Glasma α–d, #18) | — | no |
| 9 | BeAGLE incoherent shapes | blocked on FLUKA license | not blocking: |t|-fit purity argument stays band-level with e+Pb benchmark [80–99%] | WP7 | no |
| 10 | bibliography | 40 entries, individually verified 2026-08-17 | BibTeX from INSPIRE keys; re-verify volume/pages at submission | WP6 | no |

## 7.4 Work packages

### WP1 — Grid structure functions in production (◐ infrastructure exists)
- ☐ Promote the money_delta script-local `NuclearF2FromGrid`
  (EPPS21nlo_CT18Anlo_Li6, R1998 for F₁) into `polli_fastsim/structure.py`
  behind the existing `NuclearF2` interface; wire through
  `inputs.get_backends`.
- ☐ Rerun money plots 5/6/7 and `phase_space_bins.py` with `--backend grid`;
  re-solve moment_A at the grid ⟨Q²⟩; record sweet-spot drift and the
  toy-vs-grid ratio per headline number (expect ≲ ×1.5 by run-2/3 checks:
  toy-vs-CT18 F₂ within ±37%).
- Acceptance: all §7.1 numbers re-derived on grids; drift table in the plan
  addendum. Effort: 2–4 days.

### WP2 — Polarization, run-plan, and dilution bands
- ☐ δA scaling table vs P_zz and fill share (analytic 1/(P_zz√(fN)) checked
  against one sampler rerun).
- ☐ Dilution paragraph: 1/3 baseline; 0.81-cluster upside as the labeled
  alternative (amplitudes ×2.4); one email to I. Cloët (D3) — do not block.
- Acceptance: Table 1 of the letter exists. Effort: 1 day.

### WP3 — Reconstructed-level closure for cos 2φ
- ☐ Port the 20260729 smearing model (tracking σ_p/p, σ_θ, η-dependent
  ε_eID) onto the polligen φ′ pseudo-experiments: migration in (x, Q², φ′),
  reco-binned refit per super-bin.
- ☐ 2-D acceptance-hole closure: (φ, η) hole map → binned-fit bias < stat
  error at 10 fb⁻¹/u (the generator-level version already passes).
- ☐ D2: if reco dilution < ~10% and unbiased, keep main figures at generator
  level with reco factors quoted per bin; otherwise switch figures to reco
  level.
- Acceptance: per-bin dilution factors; statement "Case-3-style retention X%"
  reproduced for the cos 2φ observable. Effort: ~1 week.

**WP3 addendum (2026-08-24, `reports/reconstruction_chain_report`).** The
reconstructed-level closure is now specified end to end and seeded in
`polligen/reco.py`; the note changes three inputs of this WP:
- ☐ **Binning must be reconstructed (x, Q²) with the mixed (eΣ) method**:
  spots 1–3 of money plot 5 sit at y = 0.010–0.025 where the electron
  alone gives δy/y = 50–120% (Table 1 of the note). Use
  `reco.hadronic_y` with a 15–30% band (purity 0.75–0.83, reco-bin
  amplitude 0.96–0.99 of truth at 15–20%); quote the e′-only variant at
  y ≥ 0.05 and the low-energy configuration for x ≈ 0.1 (open question #21).
- ☐ **Estimator: spin-state-sorted ratio** (`reco.harmonic_ratio_fit`) of
  m = ±1-rich (P_zz = +0.6) and m = 0-rich (−1.2) fills, with a sin 2φ′
  term; the single-fill fit is biased by the detector's cos 2φ′ acceptance
  harmonic ÷ P_zz. δA becomes 2√(2/N)/(P₊ − P₀) = 0.67× the current
  values — re-derive the §7.1 numbers with it (the 1.5× gain is real only
  if the source delivers m = 0-rich bunches at the same purity).
- ☐ **Angles from four-vectors**: `reco.azimuth_wrt_lepton_plane` (covariant
  φ_S, = φ_e − φ_S to O(γ²)); head-on transformation applied (e′ odd
  harmonics only). The 2-D φ-hole closure of the original bullet becomes
  a smooth-ε(φ′) closure with the ratio estimator (already demonstrated
  at the super-bin level, Fig. 1d of the note).
- ☐ **Coherent (feeds WP5)**: present the anchored a₂ as a modulation of the
  recoil azimuth φ_t − φ_S and the exotic-glue term of the electron azimuth;
  fit R(α, β) in 2-D; replace the constant 0.20 GeV cut in `coherent.py` by
  the angular cut `reco.tag_pt_cut` (67% / 9% / 10⁻⁸ at 20.5 / 50 /
  137.5 GeV/u) and state the cutout geometry (open question #20); settle
  the a_n normalization convention of arXiv:2408.13213 (1 + 2Σ vs 1 + Σ).
- D2 input: the reco-level dilution is small (1–4%) and unbiased with the
  mixed method, so the main figures can stay at generator level with
  reco factors quoted — provided the binning and the estimator above are
  adopted in the generator-level plots themselves.

**WP3 results (2026-08-24, second pass — `polligen/recopseudo.py`, money
plots 5R/7R/6R, report §7).** ☑ Reconstructed-level pseudo-experiments
exist and close:
- ☑ Inclusive (mixed method, 20% hadronic y, EMCal E′, track angles,
  ε_eID, reco cuts, covariant φ′, two-fill ratio with a 3%/2% φ′-efficiency
  harmonic and a 10⁻³ rel-lumi offset on): sweet spots 1–4 in reco bins —
  purity 0.70/0.69/0.72/0.74, efficiency 0.45/0.64/0.40/0.70, D = 0.912/
  0.991/0.919/0.980; Â unbiased vs the reco-bin truth; **δÂ = 1.2 / 1.0 /
  1.8 / 3.0 ×10⁻⁴ (1 yr)** = 0.65–0.70 × the §7.1 baseline (the m = 0-rich
  fill gain beats the efficiency loss). The y = 0.01 edge bins (1, 3) lose
  55–60% of their events to the reco y cut; D = 0.91 there.
- ☑ 7R: δΔ best bins 1.0×10⁻³ (Q² = 1.14) / 0.5×10⁻³ (3.13 GeV²) in year 1
  (vs 0.7–1.4 ×10⁻³ baseline), purities ≈ 0.6.
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

### WP4 — Radiative-correction bound (not a calculation)
- ☐ Leading-log unpolarized RC weights (plans/02 step 1.4 route) applied as
  (x, Q², φ) kinematic migration on the modulated cross section → bound on
  amplitude dilution and on a fake-a₁ term.
- ☐ Self-normalization argument (single-fill modulation; no relative-lumi
  lever) written up with the bound.
- ☐ Gate: bound ≤ ~5% of amplitude → one paragraph + assumptions row;
  larger → appendix in the arXiv version and a flagged systematic band.
- Effort: 3–5 days.

### WP5 — Coherent channel presentation
- ☐ Replace two-point optics with curves vs pT_cut (0.1–0.7 GeV):
  acceptance, N_tag, best-bin δA; mark 0.20 (documented top-rigidity scale,
  ³He precedent) and 0.45 (our derivation) on the curves; state that Li
  optics are undocumented and the physics case constrains them.
- ☐ Fold exp(−B t_min) and the ×0.73 rate-weighting option into the central
  tagged-yield curve (small changes in `polligen/coherent.py` + tests).
- ☐ IR-8 panel/inset: published efficiencies d 47% / ³He 32% / ⁴He 29% /
  ⁷Li 17.8% (no ⁶Li — interpolation labeled ours), pT ≈ 0 reach.
- ☐ Geometry note: quote current RP z (32.5/34.3 m scan) alongside the
  YR-era 26/28 m, windows in θ/R unchanged.
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
| 6 | Systematics and assumptions: reco dilution, RC bound, polarimetry, purity via |t| fit, condensed assumptions | 450 | — |
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
| "No tensor radiative corrections" | quantified migration bound; self-normalization; open-theory statement with citation trail | §6 (WP4) |
| "Li beams do not exist; no luminosity" | EPIOS PRC 113:060501 feasibility; stated 10 fb⁻¹/u with linear scaling; P_zz band quoted | §3, Table 1 |
| "Coherent fraction is invented" | explicit f₀ band bracketing HERA ep and heavy-A saturation; first-of-kind labeled; IP-Glasma ask on record | §5 |
| "pT cut undocumented for Li" | curves vs cut, not a point estimate; documented anchors marked; IR-8 alternative with published numbers | §5 (WP5) |
| "Generator-level only" | reco-level dilution factors from ePIC-parameter smearing; φ-hole closure demonstrated | §6 (WP3) |
| "α+d background fakes the tag" | m-state-blind → dilutes, cannot fake; |t|-shape purity with e+Pb benchmark; Z-ID question stated as open | §5–6 |

## 7.7 Authorship, circulation, timeline (D3, D4)

- **D3 authorship** (user's call): lead C. Peng; candidates to invite —
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
| Oct 10 | WP4 bound done; gate passed or appendix planned |
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
