# Polarized ⁶,⁷Li @ EIC — Simulation Program Plans

Demonstrate feasibility and showcase the physics a polarized lithium source
brings to the EIC, in support of the ANL polarized ⁶,⁷Li ion-source program
(`docs/ecrp_2026_proposal.pdf`, NOFO DE-FOA-0003602).

## Document map

| doc | content |
|---|---|
| [01_findings_physics_case.md](01_findings_physics_case.md) | What the source docs + literature establish: the three flagship observables (gluonometry Δ, tensor b₁, polarized EMC), formulas, beam parameters, polarization inputs |
| [02_phase1_event_generation.md](02_phase1_event_generation.md) | Phase 1: fast simulation + BeAGLE — phase space, rates, FOMs, spectator-tagging purity. Step-by-step with effort estimates |
| [03_phase2_full_simulation.md](03_phase2_full_simulation.md) | Phase 2: full ePIC chain (eic-shell → npsim → EICrecon) — far-forward acceptance for Li fragments, reconstructed-level closure tests |
| [04_open_questions.md](04_open_questions.md) | External dependencies: ring spin dynamics, optics for Li, BeAGLE validity, theory curves — each with owner and default assumption |
| [05_doubly_polarized_generator.md](05_doubly_polarized_generator.md) | "polligen": doubly polarized e+⁶,⁷Li event generator — spin-density-matrix ⊗ cluster-IA kernel, reweighting + native modes, tagged spin observables, HepMC3 for ePIC |
| [06_cos2phi_coherent_backgrounds.md](06_cos2phi_coherent_backgrounds.md) | cos 2φ money plots as projected data points (sweet-spot bins); the coherent intact-⁶Li channel (RP tag, 13.5% HA acceptance) and its full background budget (α+d beam-blindness, T=1 veto pattern, Z-ID question #19) |
| [07_plb_letter_gluonometry.md](07_plb_letter_gluonometry.md) | The PLB-class simulation letter: scope decision (gluonometry, inclusive + coherent), gap analysis vs referees, work packages WP1–WP7 (grid SFs, reco closure, RC bound, coherent curves, paper production), skeleton, risk register, timeline to the INT program (submit Jan 2027) |

## The physics in three lines

1. **⁶Li (spin-1)** = practical replacement for polarized deuterons (which
   the EIC ring likely cannot keep polarized): tensor structure function b₁
   and the purely gluonic double-helicity-flip Δ(x,Q²) ("nuclear
   gluonometry"), plus medium modification via ⁶Li-vs-d comparison.
2. **⁷Li (spin-3/2)** = effective polarized proton in-medium (P_p ≈ 0.86):
   the polarized EMC effect over an order of magnitude more x–Q² than JLab,
   discriminating mean-field vs SRC origins of the EMC effect.
3. **Collider-mode spectator tagging** (α, t, n in Roman Pots/OMD/B0/ZDC)
   selects the struck cluster — impossible in fixed-target; its purity is
   the simulation deliverable the ECRP proposal explicitly calls for.

## Strategic findings from the verified research sweep (2026-06-12)

1. **Four confirmed literature gaps we can publish first**: b₁ for any
   A > 2 nucleus; any numerical EIC gluonometry (Δ) projection; α-tagged
   DIS on ⁶Li; a Li-beam luminosity/polarization parameter set.
2. **Tagging inverts between isotopes at IP6** (rigidity-verified): ⁷Li
   α-tag lands mid-Roman-Pot window (works); ⁶Li α-tag is beam-blind below
   the RP pT cutoff (needs the pT tail, p/³He channels, or IR-8 secondary
   focus); tritons have no IP6 coverage at all. ⁷Li — the isotope the
   source commissions first — is also the tagging-friendly one.
3. **No generator anywhere does polarized nuclei** — asymmetry reweighting
   on unpolarized samples (BeAGLE for breakup) is the established route;
   BeAGLE runs A = 6,7 but untuned (C-12 Fermi momentum, no cluster
   geometry) → validation + cluster-IA cross-check is mandatory.
4. **Machine feasibility is already argued in EPIOS** (arXiv:2510.10794,
   PRC 113:060501): G(⁶Li) = −0.178 handled deuteron-like, G(⁷Li) = +1.532
   with partial snakes; ~138/~117 GeV/u top energies.
5. **Calendar anchor**: INT program on polarized ion beams at EIC,
   March 22 – April 2, 2027 — target for Phase-1 money plots.

## Development run 6 (2026-08-24): the reconstruction chain, measured quantities, and the audit

- ☑ **Analysis note `reports/reconstruction_chain_report`** (HTML/PDF,
  three schematics, two computed figures): what the inclusive and the
  coherent cos 2φ measurements record; the covariant azimuth of the
  observable (= φ_e − φ_S exactly for a massless target, O(γ²) < 5 mrad
  for ⁶Li; verified against an explicit collinear-frame construction);
  the head-on frame (the 25 mrad crossing changes only odd harmonics of
  the e′ azimuth, ≈ 10⁻³, because the ePIC axis is the electron beam);
  kinematic reconstruction; the spin-state-sorted estimator; the
  Roman-Pot measurement; Δ̂ = −Â y²D_φ/(1−y) → Δ/F₁ = −2(1+R)Â at low y;
  a step-by-step audit table of the code.
- ☑ **`polligen/reco.py`** (15 tests, 83 collected in evgen): frames, covariant
  azimuths, electron/hadronic/mixed reconstruction with resolution
  models, the harmonic ratio estimator, recoil four-vector (exact
  light-cone t), Roman-Pot emulation (divergence smearing, rectangular/
  elliptical cutout, angular cut) — the seed of plans/07 WP3.
- **Findings that change the analysis** (open questions #20, #21):
  (1) three of four inclusive sweet spots sit at y = 0.010–0.025 where
  e′ alone gives δy/y = 50–120% → x from the hadronic (Σ) y, which
  polligen does not generate; at 15–20% hadronic resolution the
  super-bins keep 75–83% purity, reco-bin amplitude within 1–4% of
  truth; (2) the single-fill fit of money plots 5–7 is biased by the
  detector's cos 2φ′ acceptance harmonic ÷ P_zz (3% → 4× the signal);
  the ratio of m = ±1-rich / m = 0-rich fills cancels it exactly and is
  1.5× more precise (δA = 2√(2/N)/(P₊ − P₀)); (3) coherent: the
  anchored a₂ lives in the recoil azimuth φ_t − φ_S, the RP cutout fakes
  ⟨cos 2φ_t⟩ ≈ 0.5 at a 25% aspect ratio, and the near-beam cut is
  angular (pT_cut = 10σ_θ·6p_u): 67% / 9% / 10⁻⁸ acceptance at
  20.5 / 50 / 137.5 GeV/u for the same optics — the low-energy
  configuration is the coherent program's home; the a_n normalization
  convention of the deuteron anchor (factor 2) is flagged for a source
  check.
- ☑ **Second pass (same day): WP3 reconstructed-level pseudo-experiments**
  — `polligen/recopseudo.py` (+10 tests, 93 collected), money plots
  5R/7R (`money_cos2phi_reco.py`) and 6R (`money_cos2phi_coherent_reco.py`),
  report §7. Inclusive, mid energy, mixed method (20% hadronic y), EMCal
  + track angles + ε_eID, spin-state ratio with the acceptance harmonic
  and a rel-lumi offset on: sweet-spot purities 0.69–0.74, efficiencies
  0.40–0.70, D = 0.91–0.99, Â unbiased vs the reco-bin truth, δÂ = 1.2 /
  1.0 / 1.8 / 3.0 ×10⁻⁴ in year 1 (0.65–0.70 × money plot 5); 7R best
  bins δΔ = 1.0 / 0.5 ×10⁻³. Coherent: angular cut + r = 1.25 cutout →
  N_tag = 3.4×10⁶ (1 yr), fake ⟨cos 2β⟩ = −0.44 cancelled by the ratio,
  template fit (acceptance-weighted MC basis, a_t ∝ |t| — the plans/06
  two-component fit in response form) recovers a_t and a_e in all four
  |t| bins. **Convention verified from arXiv:2408.13213 Eq. (9):
  1 + 2Σ a_n e^{inΦ}, Φ = vector-meson (recoil) azimuth to the
  polarization axis → the deformation modulation is 2a₂; money plot 6
  injects a₂ (conservative ×2)** — `coherent.cos2phi_coefficient_deformation`.
  Estimator upgrades: luminosity-weighted ratio (exact R = σ_P²T/(1+P̄T)),
  denominator correction, empty-bin handling, 2-D fit with template basis.

## Development run 5 addendum (2026-08-10, later): merge + unified Δ models + 1/10-year projections

- ☑ **Merged the parallel `money_delta` line** (TooLate0800/
  PolarizedLithiumSim, 12 commits, July 2026): moment-constrained Δ
  ansätze, dated study scripts + notes in `fastsim/`, NumPy-2.0 grid
  fix. Disjoint from the evgen work; merge commit ffd6e52; all suites
  pass.
- ☑ **Unified Δ-model registry** `polli_fastsim/delta_models.py` —
  single home for every Δ(x,Q²) model behind one (x,q2,f1) interface:
  `toy`, `moment_A` (Δ = A·α_s·F₁·x^a(1−x)^b, A solved from
  ∫xΔdx = −0.012·α_s; ported from `money_delta_20260729.py`),
  `moment_B` (no F₁, analytic Beta solver, the conservative reading);
  shape variants; explicit per-nucleon `dilution` (plans/04 #6 mapping:
  their P_zz = 0.267 ≡ our P_zz = 0.8 × dilution 1/3). 7 new tests
  (23 fastsim total).
- ☑ **Money plots 5–6 re-cut with the moment constraint and separate
  1-year (10 fb⁻¹/u) / 10-year (100 fb⁻¹/u) projections**: inclusive
  amplitudes now (0.7–1.2)×10⁻² at the sweet spots (moment_A,
  table-α_s, dilution 1/3, P_zz 0.6) vs per-bin δA ≈ (1.5–4.5)×10⁻⁴ in year 1 — the
  A-vs-B interpretation spread (the dominant ansatz systematic) is
  itself resolvable; coherent tag: N_tag = 1.1×10⁷ (1 yr) / 1.1×10⁸
  (10 yr), best-bin δA = 1.9×10⁻³ / 6×10⁻⁴.
- ☑ Top-level README; physics note `docs/note_cos2phi_coherent_6Li.md`
  (50 verified references); PDF report in `reports/`.
- ☑ **Money plot 7** (`money_delta_extraction.py`): the extracted
  Δ(x,Q²) itself as xΔ data points at the three sweet-spot Q² slices
  (Δ̂ = Â·y²D_φ/(1−y), model bin-centering, independent 1-/10-yr
  draws) — best bins δΔ ≈ (0.7–1.4)×10⁻³ (1 yr) on Δ ≈ −0.03…−0.09.

## Development run 5 (2026-08-10): cos 2φ money plots + coherent channel

- ☑ **Money plots 5–6 — cos 2φ as projected data points** (user request):
  `evgen/scripts/money_cos2phi.py` (inclusive gluonometry: φ′ pseudo-data
  with stat error bars in 4 sweet-spot (x,Q²) super-bins, Q² = 1.1→14
  GeV² at x ≈ 0.02–0.06, + amplitude-vs-x scenario curves; per-bin
  δA ≈ (0.4–1.2)×10⁻⁴ at 100 fb⁻¹/u, P_zz = 0.6) and
  `evgen/scripts/money_cos2phi_coherent.py` (coherent e⁶Li→e′X⁶Li(g.s.):
  intact recoil is exactly beam-blind (R = 1.000) ⇒ RP pT-tail only —
  **acc = exp(−B pT_cut²) = 13.5% with high-acceptance optics, 4×10⁻⁵
  with high-divergence: the coherent program fixes the optics**; ~10⁸
  tagged events at 100 fb⁻¹/u, best-bin δA ≈ 6×10⁻⁴, 5σ floor at a
  0.3% modulation — apparently the first intact-tag projection for any
  A = 6). New `polligen/coherent.py` (+10 tests, 64 total),
  `phi_histogram_pseudo` full-luminosity binned pseudo-experiments,
  `cos2phi_fit_binned`, `effective_modulation`.
- ☑ **plans/06 background budget** for the coherent tag, anchored on a
  32-claim verified research sweep: the α+d channel is the killer (both
  fragments R = 1.000, same velocity — only dE/dx Z² separates); every
  T = 1 ⁶Li* state is γ- or nucleon-vetoable while T = 0 states feed
  the blind channel; B0 EMCal (not ZDC) catches the 3.56 MeV
  de-excitation γ; ⁶Li is a deformation null test (Q = −0.0806 fm²).
  New open questions #18 (coherent model) and #19 (RP Z-ID for
  A/Z = 2 — no EIC document addresses it).
- ☑ **Second sweep — modulation anchored** (same day): the polarized-d
  CGC calculation (Mäntysaari et al. PLB 858:139053) digitized and
  scaled to ⁶Li (only such calculation for any polarized nucleus; all
  8 forward citations checked): deformation term a₂(t) = −(P_zz/4)ε_B0
  B|t| with ε_B0 ∈ −(0.04–0.13) → ⟨a₂⟩_tag ≈ 0.036 [0.018–0.059] at
  P_zz = 0.6, sign flip vs d predicted; flat gluonic term 3×10⁻³–10⁻²
  (lattice/NPLQCD + Kumano–Song). Money plot 6 now shows the d-anchor
  vs ⁶Li band and the two-component a₂(t) decomposition strategy
  (plans/06 §6.4b; 66 tests).

## Development run 4 (2026-07-13): polligen Step 5.A

- ☑ **`evgen/polligen` created — plans/05 step 5.A complete** (kernel +
  inclusive sampler): spin-density machinery (J = 1, 3/2, any axis,
  Wigner-d/CG, population solvers incl. spin-temperature fills), doubly
  polarized inclusive master formula (HJM tensor sector + vector sector
  on the fastsim SF backends, spin-3/2 rank-2 scenario slots), run-plan
  bookkeeping (helicity flips, tensor thirds, rel-lumi offsets,
  polarimetry smear), Poisson spin-labeled sampler + Mode-W weight hook,
  analysis estimators. 35 new tests; all §5.4 inclusive gates pass
  (master formula ↔ asymmetries.py at rtol 1e-12 on toy AND PDF-grid
  backends; pseudo-experiment spreads close on every analytic FOM map,
  `evgen/closure_fom_{6,7}Li.png`; cos 2φ recovery unbiased with holey
  φ acceptance where the naive moment demonstrably fails).
- ☑ First systematics numbers: naive-estimator biases from a 10⁻⁴
  relative-luminosity offset are ≈1×10⁻⁴ on both A_zz and A∥
  (analytic formulas validated by MC at δ = 2%); lumi-corrected
  estimators remove them exactly.
- Physics guard discovered by the machinery: spin-3/2 vector fills with
  zero tensor/octupole moments violate positivity for P_z > ~0.55 —
  spin-temperature populations are now the default fill model.

## Development run 3 (2026-07-07): fast-sim evaluation + plans/05

- ☑ Fast-sim re-evaluated end-to-end on this machine: 15/15 tests pass
  (grid tests included after installing `parton` + CT18NLO/NNPDFpol11
  grids); all scripts rerun and reproduce the README headline numbers
  (⁷Li α-tag 96–99%, ⁶Li 3–9%/HA optics, gluonometry L_5σ = 15–22 fb⁻¹/u
  toy, δΔR ≈ 3.5% at x = 0.3); toy-vs-CT18 F2 within ±37% as documented.
  Findings (kept as report, not yet fixed): two inert test assertions
  (`test_smoke.py:17`, `test_spectator.py:68`); FOMs implicitly give each
  observable the full luminosity in its own spin configuration (run-plan
  split not modeled); no γ²/A₂ target-mass terms at high x, low Q²;
  F2 uses 5 flavors vs g1's 3; `money_delta` applies its min-events cut at
  the 1 fb⁻¹ reference (mildly conservative); RP z in `farforward.py`
  docstring (26/28 m) predates the 32.5/34.3 m geometry (windows in θ/R
  unchanged); DIS kinematics and spectator kinematics are sampled
  independently (fixed by plans/05 tagged mode).
- ☑ **plans/05 added — doubly polarized e+Li event generation
  ("polligen")**: reuse-vs-reinvent decision (reinvent only the
  polarized-nucleus vertex + spin-correlated cluster spectator; reuse
  PYTHIA/DJANGOH/BeAGLE for hadronization/RC/backgrounds), three-mode
  architecture (reweight | native cluster-IA | RC bands), Cosyn–Weiss
  (arXiv:2006.03033, 2603.23699/23700) as the tagged formalism, step plan
  5.A–5.E (~7–9 weeks) with validation gates, money plot 4 (first tagged
  tensor asymmetry for A > 2) and a tagged-α alignment-polarimetry bonus;
  open questions #14–17 added to plans/04.

## Development run 2 (2026-06-12, autonomous; commits abd2bce…)

- ☑ Money plots on real PDF grids (CT18/NNPDFpol): gluonometry 5σ at
  Δ/F₁=10⁻³ within 25–37 fb⁻¹/u; conclusions stable vs toys.
- ☑ Control studies on **official BeAGLE samples** (xrootd-streamed):
  routing validated at all three regions — e+d p→OMD 96.6%, e+d n→ZDC
  99.2%, e+³He p→RP 99.8%; BeAGLE pT tails 2–13× harder than the cluster
  model ⇒ R≈1 tag acceptances are model-dominated (the ⁶Li α-tag number).
- ☑ Far-forward gun scan in **current epic-main geometry** (18×275
  fields): routing table confirmed by Geant4; RP stations now at
  z=32.5/34.3 m; **discovery — the "no-coverage" triton (R=1.286)
  crosses the RP planes on the over-rigid side + showers in the ZDC ⇒
  the ⁷Li α+t double-tag may work at IP6** (10σ/reco check pending).
- ☑ Container refreshed (eic_xl, EICrecon v1.38.0); known issue: its
  pyHepMC3 rootIO segfaults — HepMC3 reading stays on the legacy image.
- ❌ Showstopper for local BeAGLE generation: **FLUKA license** (user
  registration at fluka.org) — everything else is built and validated.

## Current status (2026-06-12, end of first development sprint)

- ☑ Source docs digested; verified findings + benchmarks table in 01.
- ☑ Phase-1 fast-sim developed and tested (15 tests): rates/FOM maps,
  **first tagging-acceptance numbers** (⁷Li α→RP 96–99%; ⁶Li α 3–9% at
  IP6 — quantitative beam-blindness), **first money plots** for all three
  observables (gluonometry 5σ at Δ/F₁=10⁻³ within ~15–40 fb⁻¹/u;
  CBT-vs-TMT ≈5σ at high x with 100 fb⁻¹/u), estimator closure tests,
  and PDF-grid backends (CT18NLO, NNPDFpol11) behind the toy interfaces.
  See `../fastsim/README.md` for the headline numbers and caveats.
- ☑ Local stack surveyed: BeAGLE present but needs FLUKA (use BNL/JLab
  prebuilds); eic-shell container runnable but stale (image renamed
  `eic_xl`; current releases epic 26.06.0, EICrecon v1.38.0) → reinstall
  at Phase-2 start.
- ☐ Next actions (cheapest-first): ask Cloët about the ⁶Li
  effective-polarization convention (factor 2.4 in the g₁ FOM, plans/04
  item 6); request SDCC/ifarm access + email BeAGLE authors (long pole,
  plans/02 step 1.5); adopt EPIOS scenario numbers and digitized
  CBT/TMT/b₁ theory curves (steps 1.1–1.2); then rerun money plots on
  grid inputs.
