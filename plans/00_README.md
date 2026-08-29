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
| [../docs/reproduction_manual.md](../docs/reproduction_manual.md) | **How to reproduce any of it**: environment (Python, PDF grids, PYTHIA 8, eic-shell, a headless browser), every command with its measured runtime and expected output, the third-party generators (PYTHIA 8, BeAGLE over xrootd, ePIC/npsim), the numbers to check against, and what cannot be reproduced here and why |
| [08_simulation_chain_completion.md](08_simulation_chain_completion.md) | Completing the simulation chain (2026-08-25): the 23 gaps that survived an adversarial audit of the reconstruction chain, the kernel, the hadronic final state and the fast simulation — ordered, with the convention items reserved for the author (D1, D7, D9) and the externally blocked tail (D2–D8) separated out |
| [09_nearbeam_nanowire_far_forward.md](09_nearbeam_nanowire_far_forward.md) | A near-beam layer for the far-forward lithium tags (2026-08-26; §9.0 re-derived 2026-08-28 on the Yellow Report divergences — at the published optics the machine binds everywhere and a closer approach buys nothing; at the tagging optics a layer that follows the 0.12–0.33 mrad envelope is the difference between no tag and a 32–42% tag), whether a superconducting nanowire can deliver it, the hot-spot firing-threshold answer to open question #19, the obstacle table, and the correction that the ePIC pot geometry has moved since the snapshot `tools/fullsim` measured — report in `reports/nanowire_far_forward.html` |
| [10_beam_divergence_light_ions.md](10_beam_divergence_light_ions.md) | **The beam energies and the divergence the whole far-forward programme rests on** (2026-08-27): every far-forward acceptance is exp(−B(10σ_θ·A·p_u)²). Two corrections. **The energies**: EIC ions are γ-matched, not rigidity-scaled — the rings must share a revolution period, so ⁶Li sits at **40.8 / 99.5 / 137.5 GeV/u**, not 20.5 / 50 / 137.5. **The divergence**: σ_θ was one energy-independent, isotropic, proton-derived 72.7 µrad; YR Tables 10.1/10.2 give it per configuration and optics, and the species step applies *only* where rigidity binds — so ⁶Li carries the proton's **220/380 and 180/180 µrad** at the two lower configurations and pays √2 only at the top (**92/92**). Together these cost the coherent tag six to twenty-four orders of magnitude, and the recovery needs a two-ring β* de-squeeze the machine may not have |

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
   focus); the over-rigid triton IS taggable — it lands on the inner-side
   outer band of the pots at +66 mm at every configuration and then in
   the ZDC (measured 2026-08-28, run 14; `farforward.over_rigid_route`),
   which retires the "no IP6 coverage" of the earlier routing. ⁷Li — the isotope the
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

## Development run 15 (2026-08-29): the target-mass term on, the light-ion lattice and the tagging optics decided, the coherent chain on its own transport, the theory curves on a data-driven baseline, the ⁷Li theory note

The author's decisions of 2026-08-29 — switch the finite-γ term on if it
is cheap, choose the ⁶Li lattice and the tagging optics as educated
guesses that maximise the physics and minimise the machine work, write
the ⁷Li theory questions down, regenerate the inclusive figures once at
the end — were carried out in six streams, reviewed, verified through the
five lenses of run 14 and repaired.  Before that, every registered figure
had been regenerated on the run-14 libraries (only two report numbers had
gone stale; ~30 pre-existing mis-roundings were corrected) and a
twenty-fifth check now fails whenever a figure is older than a module its
script imports.

- ☑ **The exact target-mass kinematics is in the chain by default.**
  A∥ = D_γ(A₁ + ηA₂) with A₁ = (g₁ − γ²g₂)/F₁, A₂ = γ(g₁ + g₂)/F₁ and
  g₂ = g₂^WW from one shared quadrature (`asymmetries.a_parallel_exact`,
  cached per grid, so the term costs nothing per call);
  `InclusiveKernel(target_mass=True)` is the default and `fom` extracts
  g₁/F₁ with the same kinematics, so the 0.12–1.06 % bias on ΔR is gone
  (δΔR 4.2 / 4.0 / 5.9 / 18.7 % on the grids).  What is left is the
  twist-3 uncertainty on g₂, ≤ 0.16 % on ΔR at g₂ = 0 and ≤ 0.08 % at
  1.5 g₂^WW; the tagged-triton overlay carries its 2.1 / 5.0 % as a term
  with a residual of 0.015 / 0.073 %.  `target_mass=False` stays
  bit-for-bit reachable.
- ☑ **The lattice and the tagging optics are decided, as educated
  guesses stated as such.**  A ⁶Li fill at 5 × 41 or 10 × 100 runs, in
  our assumption, in the ePIC baseline lattice of that ring setting at
  twice the field — the same magnets and orbit, the transport the scan
  measured (R₁₂ = 19.24 / 21.25 m) — with the Yellow-Report-scaled
  Z/A = ½ file (R₁₂ = 29.8 m, edge 1.60 mrad) carried as the alternative;
  consequently the per-configuration levers are the baseline everywhere,
  including `tagging_optics_point` (the 18 × 275 lever reachable).  The
  5 × 41 vertical lever was measured on a zero-insertion copy of the
  geometry (`tools/fullsim/README.md`): R₃₄ = 4.56 m, so the real
  29.6 mm insertion is 6.49 mrad, beyond the 5 mrad outer bound of the
  pot plane — the plane is shut by the pipe (two of 2 × 10⁵ α spectators
  clear it, both past 5 mrad), and the α + d separations become 17.3 /
  10.7 / 10.9 mm.
  The tagging optics is now a requirement list (Report 1 §6.1, plans/10):
  β*_x raised ×46.5 / ×164 / ×89 (≈ 42 / 100 / 71 m), a horizontal
  envelope of 0.36 / 0.19 / 0.12 mrad, the pots inserted to 7.0 / 4.1 /
  3.5 mm at station 1 (the vertical needs nothing: 17.3 / 6.0 / 2.7 mm
  against the current 29.6 / 7.1 / 2.7), the IP-to-pot transport kept at
  the baseline levers, luminosity 1/6.8 / 1/12.8 / 1/9.5 — achievable with
  moderate effort, since the Yellow Report already carries two optics per
  configuration.
- ☑ **The coherent chain re-run on its own transport.**  The dispersive
  term with each configuration's (R₁₂, D) widens the tagging envelope at
  the two lower configurations (0.33 → 0.36, 0.17 → 0.19 mrad): ε 0.423 /
  0.323 / 0.332 → 0.374 / 0.251 / 0.332, N_tag 2.59 / 3.01 / 6.15 × 10⁶ →
  2.37 / 2.42 / 6.15 × 10⁶ per year, the 5σ floors 1.74 / 2.34 / 1.62 %,
  the years to 5σ on a 1 % term 3.0 / 5.5 / 2.6, the combined one-year
  δa_e 0.00121 / 0.00111 / 0.00074; Report 2's Table 5 moves in its
  5 × 40.8 block and not at all at 18 × 137.5; the ⁶Li α tag at the
  tagging optics is 0.315 / 0.223 / 0.291 and the ⁷Li α tag 0.987 /
  0.991 / 0.994.  The δa_t/a_t "move" of the WP5 scan that the refresh
  reported was a stale run-13 quote: the exact Jacobian of the ratio
  inversion (8a98f91) had correctly reduced the error by up to 12 %, and
  Table 5 was on it already.
- ☑ **The unpolarized EMC baseline is data-driven, and the b₁ signal is
  binned in Q².**  `unpolarized_emc_ratio` defaults to the per-nucleon F₂
  of EPPS21 (Li-6) over CT18ANLO's free isoscalar nucleon — the proton
  fit the modification is defined against — at Q² = 5 GeV² (nNNPDF3.0,
  the digitized CBT ⁷Li curve and the legacy table reachable), with an
  EPPS21 Hessian band; on it the two camps' valence strengths are 0.53
  (CBT) and 0.21 (TMT), the CBT–TMT separation halves to 0.023 / 0.021 /
  0.018 / 0.032 and the discrimination to 0.55 / 0.53 / 0.31 / 0.17σ per
  bin at 10 fb⁻¹/u (1.8 / 1.7 / 1.0 / 0.5σ at 100; the
  isotope caveat, Li-6 grid against a ⁷Li observable, is quantified against
  the CBT curve).  `money_b1` draws the signal at each bin's Q² from
  Miller's Q² set: |A_zz| 2.2 × 10⁻⁴ … 3.3 × 10⁻³, 1.6 / 4.7 / 7.3 / 10.4 /
  5.5 / 5.8σ per bin (the move against the fixed slice is the bins' own
  rate-weighted y and ⟨Q²⟩ in F₁, not Miller's Q² dependence, which is
  flat in every accepted bin).
- ☑ **The ⁷Li theory questions are written down**
  (`docs/note_7li_theory_questions.md`): the spin-3/2 rank-2 basis and its
  normalisation, the coherent cos 2φ amplitude beyond linear order for a
  deformation fifty times the deuteron's, the α + t P-wave overlap with its
  m-dependence, the coherent excitation of the 477.6 keV state and whether
  its photon can be vetoed, FSI and two-body currents in the tagged α + t
  channel, the tensor-sector radiative corrections, and the ⁶Li asks that
  go with them, each with what the simulation needs, what exists, and to
  whom the ask goes.
- ☐ **Left**: D1 (the b₁-sector sign); the ⁷Li theory asks themselves;
  whether the far-forward group's light-ion lattices confirm the ×2-field
  assumption; a measured R₁₂ for a de-squeezed lattice.

Tests: 312 evgen + 119 fastsim, 25 consistency checks.

## Development run 14 (2026-08-28): the pot aperture re-measured in the current ePIC geometry, the seven-bin |t| window, the run-plan share, the target-mass bound, the theory curves digitized

Run 13 left nine items.  The current `eic_xl-nightly` container
(`epic-main` 9aaa2969, 2026-08-22) was installed at
`~/Projects/eic-2026/`, which unblocked the first and the two that
depended on it; the rest were closable in code.  Six streams were
implemented, each reviewed adversarially and repaired, the far-forward
measurement reviewed twice, and the whole round passed through a
five-lens verification (numbers, stale values, physics, register,
reproducibility) before the rebuild.

- ☑ **The Roman-pot aperture is measured in the current geometry, per
  configuration, and the over-rigid triton is taggable** (plans/09 B1,
  D3, the R > 1 branch, R₁₂/R₃₄).  The pots are now 1.6 cm modules in
  four horizontal bands (|x| ≤ 16, 16–32, 32–48, 48–144 mm) held off the
  beam by per-configuration insertions (2.7 mm at 18 × 275 for the
  central band only; 7.1 / 5.5 mm at 10 × 100; 29.6 / 27.5 / 18.0 mm at
  5 × 41), and the hit reader (`tools/fullsim/ff_gun_hits.py`) now
  resolves the four station planes instead of averaging them.  Ladders at
  0.05 mrad steps on the ring reference rigidity give R₁₂ = 19.24 /
  21.25 / 29.97 m (station 1; the array is centred on the reference orbit
  to under a millimetre), R₃₄ = — / 3.35 / 2.93 m, and the dispersion
  D = 0.31 / 0.29 / 0.29 m from the α (R = 0.857), the intact ion and the
  triton (R = 1.286) — the 0.30 m the code carried is confirmed at every
  configuration, with a quadratic term that the triton needs.  The edges
  close on the band thresholds divided by the lever to 1–4 %: horizontal
  2.50 / 1.51 / 0.53 mrad, vertical (shut) / 2.12 / 0.92 mrad — at 5 × 41
  the whole 0.2–6 mrad vertical ladder reaches silicon once, so the plane
  is encoded as shut; the outer edge is 2.85–4.30 mrad, not 5 mrad, the
  ion breaking up on the pipe beyond it.  Against the 10σ envelope the
  silicon is 1.14× / 0.84× / 0.58×: it binds at 5 × 41 and the machine
  binds at the other two (the 0.91× / 0.75× / 1.12× of run 13 withdrawn).
  The R = 1.286 triton lands on the inner-side outer band at +66 mm in
  60 of 60 events at all three configurations and then in the ZDC in
  80–98 % of them, so the repository's "lost (over-rigid)" — Report 3's
  Table 6 assumption — was wrong: `farforward.over_rigid_route`,
  `route_charged(pot_config=…)` and a beam-generic
  `coherent.fragment_route_label` with `LI7_BREAKUP` carry it, and the
  ⁷Li α + t double tag is possible at IP6, with a triton acceptance hole
  at θ_x ≈ −1.6 to −2.5 mrad where the dispersion is cancelled.  The ⁷Li
  α tag at R = 0.857 is clean at 10 × 100 and 18 × 275 and a knife-edge at
  5 × 41 (57 %, station 2 only, on a module boundary).  With the
  per-configuration blind block the ⁶Li α tag at the published optics is
  0.0170 / 0.0163 / 0.0289 and the α + d separations at the pots are
  25.7 / 10.8 / 11.1 mm (Report 4 Table 5, pinned).  `epic-main` also
  ships a Z/A = 0.5 lattice at the γ-matched 82 GV for 5 × 41
  (`beamline_5x41_He4.xml`, the Yellow Report 275 GeV magnet set scaled),
  which gives R₁₂ = 29.8 m and a 1.60 mrad edge against the re-matched
  proton lattice's 19.2 m and 2.50 mrad: the published numbers stand on
  the ePIC baseline compact files and carry the alternative as
  `RP_APERTURE_MEASURED_LIGHT_ION_LATTICE` / `POT_LEVERS_LIGHT_ION_LATTICE`
  (×0.64 on the edge); which lattice a ⁶Li fill runs in was left as a
  question for the far-forward group *(decided as an educated guess in
  run 15: the baseline lattice at twice the field)*.  The tagging-optics
  envelope kept the 18 × 275 lever with the per-configuration levers
  opt-in *(the default since run 15)*.  The September-2024 table survives as
  `RP_APERTURE_SEP2024`; `tools/fullsim/README.md` carries the scans as
  the reproduction recipe.
- ☑ **The seven-bin |t| window 0.017–0.25 GeV² is the published coherent
  default** (`recopseudo.T_EDGES_PUBLISHED`; the four-bin 0.05–0.25
  window reachable as `T_EDGES_LEGACY` via `--t-edges`).  The combined
  one-year δa_e falls 0.00205 / 0.00169 / 0.00132 → 0.00119 / 0.00104 /
  0.00074 at the three configurations, the tagged sample inside the
  window rising from 27–34 % to 78–88 %; the 0.006–0.017 bin stays out
  (8–12 of 24 β cells empty, condition number 10–21, −29 % on a_t for a
  10⁻³ envelope split).  The 200-draw ensembles at 5 × 40.8 give ratio
  pulls +1.2 … −41.5 and likelihood pulls within ±1.6 over the seven
  bins; the cutout-split sensitivity is −5.5 / −2.2 / −1.4 % in the three
  added bins against −0.8 … −0.3 % in the published four; four of the
  seven bins have t_ref below the deformation anchor's lowest digitized
  point (0.05 GeV²), the fourth being the retired window's own lowest bin,
  which is now said wherever a_t is quoted there.  A per-bin design
  diagnostic (counts per cell, populated β bins, condition number,
  measured rank) is printed by every run.
- ☑ **The run-plan share is parameterized and stated, not chosen.**
  `fom.Scenario.run_share`, `--run-share` on the six fast-sim scripts and
  `--lumi-fraction` on the four evgen projection scripts (a separate,
  printed factor beside the optics penalty and the spin-state share), with
  output-stem guards and tests of the 1/√share law and the invariance of
  L_5σ in fb⁻¹/u; every report now states that each projection assumes
  the full luminosity in its own spin configuration, optics and isotope,
  and a consistency check keeps the sentence from being dropped.
  plans/07 WP2 tabulates the options — a ⁶Li year split f_coh = 0.29 /
  0.36 / 0.22 to the tagging optics (5σ on the coherent shape term in the
  year) against the high-acceptance remainder shared one, two or three
  ways; ⁶Li and ⁷Li half a year each; both — and the two second-order
  effects.
- ☑ **The target-mass term is bounded and available.**  The exact
  finite-γ A∥ = D_γ(A₁ + ηA₂) with g₂ = g₂^WW is in the kernel as
  `InclusiveKernel(target_mass=…)`, default off (the fast simulation is
  unchanged); at small y the correction collapses to (1 + γ²) independent
  of g₂.  γ² ≤ M²/(W²_min − M²) = 0.097 by the W² cut; at the sweet spots
  0.0056 (mid, Report 1's ≤ 0.006 exact) and 0.0049 (top), 0.025 at the
  low configuration; on the polarized-EMC ΔR the weighted bias is
  0.12–1.06 % against 5–12 % statistics, on the tagged-triton A∥ overlay
  2.1 % at the published configuration and 5.0 % at 5 × 40.8.  plans/06's
  "x ≤ 0.06, γ² < 10⁻²" was false at the γ-matched low configuration and
  is corrected.
- ☑ **The theory curves are digitized, and D7 / D9 are closed.**
  `tools/digitize_figure.py` reads the vector paths of the published
  figures (CBT ⁷Li Fig. 6 at Q² = 5, both Eqs. 23 and 26; TMT Fig. 4;
  CDKS Fig. 4/5; Miller Fig. 5/6) into `fastsim/polli_fastsim/data/`
  with the recipe in `SOURCES.md`.  The two-camp comparison is made on one
  ⁷Li baseline with TMT's nuclear matter carried by a valence strength
  factor 0.397 on the legacy table *(0.21 on the EPPS21 baseline of run 15,
  where the CBT strength is 0.53)* (the pointwise ratio of effects is singular where the
  unpolarized ratio crosses one): the CBT ratio of effects is 2.25 /
  1.69 / 1.41 / 1.14 at x = 0.40–0.60 against TMT's 1.01 / 0.98 / 1.00 /
  1.08, so the discrimination sits at x ≈ 0.35–0.45 and not in the
  highest-x bins — 0.84σ per bin at 10 fb⁻¹/u and 2.7σ at 100 inside the
  valence window (the README's "≈ 5σ" retired).  The b₁ signal now uses
  Miller's absolute b₁ with the rank-2 transfer 0.9219 × ⅓ (D9: the
  vector 0.87 was the wrong rank and the 2/6 was missing, a factor 2.8
  between signal and error bars): |A_zz| = 2.4 × 10⁻⁴ at x = 0.005 to
  1.4 × 10⁻³ at 0.07, i.e. 1.7 / 4.8 / 7.6 / 10.9σ per bin at P_zz = 0.6,
  the convolution camp below 0.2σ — the b₁ case is weaker than Report 0
  said and now says so.  D7: `g1_nucleus` weights by Z and N as `f2a`
  does, the ⁷Li constants rescaled bit-for-bit, ⁶Li's double dilution gone
  (g₁(⁶Li)/g₁(d) 0.119 → 0.358 per nucleon; the value 1/3 vs 0.81 stays
  plans/04 #6), the triton mirrored on ³He.  The Cosyn–Weiss deuteron gate
  is closed analytically: the P₂(cos θ_k) factor to 10⁻⁵, the envelope
  peak at k = 0.310 GeV/c against their 0.30, A_T∥ = −2A_zz^wf reaching
  +1.000 / −2.000.  The digitized b₁ also exposed two defects: frozen at
  its last point it made b₁/F₁ diverge at the generator's x = 0.955 cell
  (now tapered with F₁'s (1 − x)³), and in the tagged deuteron kernel an
  inclusive b₁ double-counts the S/D interference the sampler already
  carries — money plot 4 now runs without it, and its 8 × 10⁶-event
  closure on the wave-function truth improves from −0.008 to +0.005 at
  k = 0.325 GeV/c (Report 4 §2.1: +0.49 / −0.07 at 4 × 10⁵ events).
- ☑ **The E − p_z window stays a documented contingency; the coherent
  ⁷Li channel has its scoping memo.**  Seed-averaged over the eight
  published seeds, the window keeps 86.9 % of the non-radiative rate on
  the Gaussian y stand-in but 97.9 % through the PYTHIA final state with
  the calibrated Σ scale, 99.5 / 97.2 % of the loss lying above y = 0.2
  and 0.01–0.06 % at the four sweet spots, with the mitigated bias
  +0.17 … +0.20 %; it is free where the numbers are and costs the low-x
  lever of money plot 7R, and the gate passes without it.  plans/09 §B3a
  records why the ⁶Li deformation model cannot be rescaled to ⁷Li
  (Q/⟨r²⟩ 56× larger gives |ΔB₀| > B and c₂ > 1 by |t| = 0.018 GeV²),
  the J = 3/2 polarization variable and the missing flip plan, the optics
  conflict with the ⁷Li run plan, the 477.6 keV level whose photon boosts
  to 42 / 102 / 121 MeV against the B0 EMCal threshold, and the theory
  asks.
- ☐ **Left**: D1 (the b₁-sector sign, the author's); whether to switch
  `target_mass` on and whether `fom` should divide by D_γ; the
  unpolarized EMC baseline (EPPS21, plans/02 step 1.2.1) and a bin-by-bin
  b₁ signal in `money_b1`; which 5 × 41 lattice a ⁶Li fill runs in, and
  R₃₄ there; the tagging optics' own lever; the ⁷Li channel's theory asks
  (plans/09 §B3a); the evgen cos 2φ figures were not regenerated after
  the b₁ default changed (the amplitude is exactly unchanged and the
  normalisation moves by ≤ 0.05 %).

Tests: 305 evgen + 111 fastsim, 25 consistency checks *(at the end of run 14)*.

## Development run 13 (2026-08-28): the to-do audit, the low-count estimator, the α+d veto, the tagged chain on the real optics, the radiative bound

Every ☐/◐ item in plans/02–10 and this log — 319 of them — was classified
against the code, the reports and the git history, and each classification
handed to an adversarial verifier: 133 were open and closable here, 100
external, 43 done but unticked, 28 superseded, 15 the author's.  The
closable ones were then worked in three rounds, each deliverable reviewed
adversarially and repaired.  Reports 0–4 now carry the author line
"C. Peng and J. Zhou" and a "Writing assisted by Claude (Anthropic)" line.

- ☑ **The low-count bias of the coherent estimator is removed** (plans/08
  A12).  The bin-wise ratio inverts R = σ_P²T/(1 + u + P̄T) bin by bin, and
  the inversion is strongly curved (g″ ∝ P̄ < 0 for the flip plan), so each
  (α, β) bin carries a Jensen offset ∝ 1/ν_b that is flat in α and lands on
  a_t: over 200 one-year pseudo-experiments at the tagging optics (5 × 40.8)
  the ratio's a_t is off by +0.1 / −0.2 / −4.0 / −34.3 % in the four |t|
  bins (pulls −9.0 and −42.1 in the last two), and at one count per cell it
  changes sign.  `reco.harmonic_likelihood_fit_2d` — the Poisson likelihood
  with the per-bin acceptance profiled out, which is exactly the
  conditional multinomial given the bin totals, whose score has zero mean
  at any count — returns pulls +0.5…+2.0 on the same draws, +0.181 where the
  ratio gives −0.085, errors that reproduce the ensemble spread, Asimov
  errors within 0.24 % of the ratio's and the blind systematics unchanged;
  on Table 5's twenty draws its means are 0.0902 / 0.1176 / 0.146 / 0.173
  against 0.0898 / 0.1181 / 0.1473 / 0.1803 injected.  Coarser bins only
  attenuate (−34.3 → −13.9 → −7.8 % at +17 % on δa_e; 4 × 8 is
  rank-deficient because ⟨cos 2α⟩ vanishes).  `--fit likelihood`; the ratio
  stays the default and the published record.  With it, the in-situ
  (u₁, u₂) fit of A3 (`--u-in-situ`, against the response's acceptance —
  the only way u is identifiable): δu₂ 1.8–15× tighter than the ZEUS 1σ,
  the u₂ leakage into a_e down from 7–9 % to 0.5–6.3 %.  §8.4's
  "anisotropic basis" item is superseded (the design is near-orthogonal at
  the tagging optics, cond 1.8–2.9); what survives is that the fixed |t|
  window discards 67–74 % of the tagged sample below 0.05 GeV² — four
  added bins would halve the one-year δa_e, 0.00205 → 0.00105, recorded as
  a follow-up *(adopted in run 14 as the seven-bin 0.017–0.25 GeV²
  window)*.  The WP5 optics scan no longer loses its figure to a
  rank-deficient last point, and a stalled Newton step is named.
- ☑ **The α + d partner-fragment veto is quantified, and Report 4's
  Table 5 was 2× wrong** (plans/09 B4).  Nothing sampled both fragments;
  `spectator.breakup_lab_kinematics` boosts one relative momentum twice
  (θ_d → 2θ_α at the opposite azimuth) and `farforward.separation_at_pots`
  carries the dispersive term the earlier estimate dropped.  The separations
  at the pots are 38.4 / 18.5 / 15.1 mm (5 × 41 / 10 × 100 / 18 × 275); the
  published 73.4 / 30.1 mm carried the retired 20.5 / 50 GeV/u through a
  derived millimetre no drift check could see — one now does (23 checks).
  Given an α that fakes a coherent tag, the deuteron is in acceptance
  0.85 / 0.84 / 0.84 of the time at the tagging optics (0.12 / 0.02 / 0.25
  at the Yellow Report optics), a recorded pair merges only rarely (≤ 4 × 10⁻⁴ in one
  pixel), and the veto collapses at 5 × 41 once the pot's outer edge is
  below 2 mrad — a station-layout question for B1.  plans/06 handle #3 is
  upgraded from "rare, not worth relying on" to the strongest handle after
  the |t| shape.  `nearbeam_two_hit.py`.
- ☑ **The tagged chain runs on the real optics** (plans/09 B2, B3; the
  last two published figures on the retired 73 µrad).  `polligen.tagged`
  now hands the spectator's lab azimuth to the rectangular envelope (without
  it the tagging-optics cut degenerated to a circle, ×1.7 too generous) and
  defaults to `yr_optics`; both scripts take `--optics/--config`.  ⁶Li α tag
  on the tagged sampler 0.0285 / 0.0250 / 0.0267 (YR high acceptance) →
  0.381 / 0.305 / 0.312 (tagging) at 1/7.1, 1/13.3, 1/9.5 of the
  luminosity — ×1.9, ×0.92, ×1.2 in tagged events (0.0284 / 0.0247 /
  0.0265 → 0.382 / 0.306 / 0.313 after run 14's per-configuration blind
  block and b₁-free tagged kernel; 0.315 / 0.223 / 0.291 on run 15's
  per-configuration transport) — and the reach is the
  result: at every published optics not one accepted spectator lies below
  k = 0.15 GeV/c (minimum 0.184), while at the tagging optics the median is
  0.145–0.162 with 44–52 % below 0.15, so money plot 4 becomes a curve.
  The published ±0.5 swing of the folded A_zz at k ≈ 0.3 was θ_k
  sculpting by the envelope (⟨|cos θ_k|⟩ 0.71–0.79 against 0.40), not the
  wave function; the overlay is now the acceptance-weighted prediction.
  The two ⁶Li samplers differ by the D wave alone (P_D = 0.0867).  ⁷Li:
  0.961–0.974 → 0.981–0.992 for 1/8–1/15 of the luminosity — a strict net
  loss (×2.8–3.9 on the error bars at equal time), so ⁶Li and ⁷Li want
  different optics and are different runs; `rp_aperture_for` is keyed by
  configuration for both isotopes.
- ☑ **The radiative-correction bound exists** (plans/07 WP4, plans/08 D3,
  plans/02 step 1.4).  `polligen/radiative.py` radiates a collinear photon
  off the incoming electron (exponentiated leading log) and lets the
  analysis reconstruct with the nominal beam: the covariant azimuth is
  invariant (3.6 × 10⁻¹⁵ rad over the massless flat sample, pinned; 2.6 × 10⁻² rad over the response itself, faking cos 2φ′ at 9 × 10⁻⁸), the mixed x is exact and the Q²_e label
  migrates by 1/(1 − z), Q²_DA by 1/(1 − z)².  With common random numbers
  the ISR-free corrections leave Δ̂ high by +0.62 / +0.50 / +0.94 / +1.22 % (mean ± sem over eight response seeds)
  at the four mid-configuration sweet spots in the published generator
  window, ≤ 2.9 % once the low-Q² feed-in is generated; a HERA-style
  E − p_z window would bring it to ≤ 0.25 % while keeping 87 % of the non-radiative rate
  on the Gaussian y stand-in (98 % through the PYTHIA final state; run 14).  The
  ≤ 5 % gate passes; the tensor-sector RC stays unmodelled and is said so
  wherever the bound is quoted.  `--isr`.
- ☑ **Far-forward and fast-sim hygiene.**  `yr_divergence_for` and its
  blanket √2 deleted (the test that pinned it now pins `sigma_theta_for`);
  the IBS normalisations stated once (Z⁴/A² per particle, Z³/A² at fixed
  current); the 149 µrad comments retired; `coherent.fragment_rigidity` on
  the physical masses (R = 0.998 / 1.005); two inert `or True` assertions
  made real; the 10σ envelope made a rectangle at every aspect ratio
  (`Optics.clears` fell back on its inscribed circle wherever σ_h = σ_v,
  discarding the azimuth at 10 × 100 and 18 × 275: the ⁶Li α tag at the
  Yellow Report optics is 1.7 / 1.5 / 1.6 %, the α + d fake rate
  0.0020 / 0.0001 / 0.0007 and its veto 0.12 / 0.02 / 0.25); the measured
  aperture priced against the per-configuration
  envelopes (0.91× / 0.75× / 1.12× — the envelope binds at the two lower
  configurations, the "never binding" of run 8 withdrawn); `tools/pythia8`
  defaults and commands on the γ-matched menu with the stale-energy guard
  extended to `tools/`; F₂'s five flavours against g₁'s three documented
  (charm + bottom are 7.8 % of F₂ event-weighted, 0.65 / 0.23 % in the EMC
  window); `money_delta`'s min-events mask applied at the evaluated
  luminosity (nothing moves: 16.7 / 16.3 / 21.8 and 66–155 fb⁻¹/u stand).
- ☑ **Documents.**  plans/10 §10.1 marked superseded, β* = 13 m → 45 m,
  §10.5 restated; plans/04 #3 / #11 / #19 / #20, plans/06, plans/09
  restated; 28 markers ticked and 11 superseded in plans/02, 03, 05, 07;
  Report 0's α tag on the current optics; Report 3's triton "no coverage"
  qualified (a routing assumption the gun scan contradicts); β = 0.30
  defined at its first use.  Reports rebuilt.
- ☑ **Left** at the time, in order of value: B1 (the aperture on the
  current `epic-main`); the fast-sim run-plan split; the γ²A₂ target-mass
  term; digitized CBT/TMT/b₁ theory curves; the R > 1 triton branch;
  per-configuration R₁₂/R₃₄; a coherent ⁷Li channel; the |t| window below
  0.05 GeV² and the E − p_z cut; D1 / D7 / D9 — all taken up in run 14
  (D1 and the ⁷Li channel's theory input remain the author's).

Tests: 276 evgen + 81 fastsim, 23 consistency checks.

## Development run 12 (2026-08-28): the divergence fixed everywhere, Reports 3 and 4 rewritten as papers

The last place the single proton-derived 72.7 µrad still drove a published
number was the fast simulation's spectator routing — the 1.85% azimuth-blind ⁶Li α-tag
of the README — and the near-beam scripts of plans/09.  Both now run on
the per-configuration optics (plans/10 A3, A4 ☑):

- `farforward.Optics` carries (σ_h, σ_v) and a luminosity fraction;
  `route_charged` / `acceptance_summary` take each fragment's azimuth and
  apply the rectangular 10(σ_h, σ_v) envelope; `yr_optics(config)` and
  `tagging_optics(config)` build the Yellow Report rows and the Report 1
  §6.1 optimum; `sigma_theta_for`, `hole_acceptance`, `tagging_optics_point`
  moved from `polligen.reco` into `farforward` and reco delegates
  (`fastsim/tests/test_optics_20260828.py`, 5 tests).
- **⁶Li α-tag** at the Yellow Report optics: 1.7 / 1.5 / 1.6% at
  5 × 41 / 10 × 100 / 18 × 275 — the near-beam tail is inside the envelope
  at every configuration and only the slice below R = 0.95 survives; the
  tagging optics recovers 35 / 27 / 28% at 1/7–1/13 of the luminosity
  *(31 / 22 / 29 % at 1/6.8–1/12.8 on run 15's transport)*.
  The legacy 73 µrad gave 17 / 2.5 / 1.8% (20 / 2.9 / 1.9% as it was
  applied, azimuth-blind); the old 1.85% was close at the
  top energy by coincidence (92 ≈ 73 µrad) and 12× wrong at the bottom.
  ⁷Li α: 97% at every optics.
- **Near-beam study re-derived** (`nearbeam_aperture_scan.py`,
  `nearbeam_reach_gain.py`, `nearbeam_sensor_budget.py`): the ×26 / ×569
  gains are withdrawn — at the published optics the machine binds at
  every configuration (coherent 7×10⁻⁸ / 6×10⁻²⁷ / 4×10⁻¹⁴ at the
  envelope); at the tagging optics the silicon at 1.0–2.0 mrad tags zero
  and pots following the 0.33 / 0.17 / 0.12 mrad envelope tag 0.41 / 0.31
  / 0.32 with four clean |t| bins each (δa_t 0.0035 / 0.0034 / 0.0023 in
  the lowest).  The gain lives in a strip d50 = 183 / 69 / 50 µrad wide.
  The energy ordering reverses: the top configuration leads.
- A fact-check workflow (51 agents) on the two reports then caught: the
  Yellow Report high-divergence rows were read from the wrong columns
  (e5 × p100 and the 1160-bunch 275 GeV; now 220/220 and 150/150, so ⁶Li
  212/212 at the top), the two near-beam scripts computed the α tag
  differently (both now route through `farforward.acceptance_summary`:
  0.35 / 0.27 / 0.28 at the tagging optics), the β* ask at 5 × 41 is ≈ 45 m
  (0.9 m × 49.7) and not the superseded isotropic 13 m, the ⁶Li deuteron
  spectator sits at R ≈ 1 (not 0.5), and the α + d merge column of Report 4
  was reversed.  Report 1 §6.3 now quotes the tagging-optics closure.
- `coherent_optics_scan.py` panel (d) is per-configuration YR curves with
  the tagging points; `eic_beam_figures.py` marks the legacy line.
- **Report 3** (`reports/eic_epic_reference.html`) rewritten as a paper
  with a new §5 / Table 6 — the lithium tags at each optics — and §4.2
  stating the tagging optics as the programme's one machine request;
  **Report 4** (`reports/nanowire_far_forward.html`) rewritten as a paper
  on the corrected divergence: the layer and the optics are multiplicative
  levers, the nanowire verdict unchanged, the charge-ID and geometry
  sections condensed.
- Docs: README α-tag line, `fastsim/README`, `docs/reproduction_manual.md`
  §3.2 / §4.6 / §7, plans/04 #20, plans/09 §9.0, `tools/fullsim/README`.

## Development run 11 (2026-08-28): the toolchain reviewed, and Report 2 rewritten as a paper

Nine independent readings of the reconstruction chain (truth leaks in the
inclusive and coherent analyses, the hadronic-final-state layer, frames
and kinematics, estimator statistics, the kernel and sampler, code-versus-
document consistency, experimental realism, tests and scripts) produced 65
findings; each code-level finding was re-verified numerically before it
was acted on.  Record: `docs/code_review_2026-08-28.md`.  What changed:

- ☑ **The hadronic-scale calibration read the truth.**
  `HFSResponse(calibrate=True)` divided each pseudo-event's measured sums
  by the library's captured fraction of the event's *true* (x, Q²) cell.
  It now uses a map in bins of the library's own *reconstructed*
  (x_mixed, Q²_e), looked up at the pseudo-event's uncalibrated
  reconstructed point — what an analysis does.  Purities
  0.56 / 0.59 / 0.64 / 0.75 against 0.52 / 0.59 / 0.59 / 0.76: not
  optimistic, but not measurable.  A test forbids the truth lookup.
- ☑ **The PYTHIA library was truncated at x = 16/s.**  PYTHIA's default
  `PhaseSpace:mHatMin = 4 GeV` applies to DIS (m̂² = x s), so nothing was
  generated below x = 0.004 at 10 × 99.5 — 39% of the selected generator
  rate, none of it at the sweet spots, all of the low-x half of the
  Q² = 1.14 slice.  Regenerated with `mHatMin = 0.5` (8 M events, five
  minutes): σ_gen +37–48%, the x spectrum follows the rate map to ±25%.
- ☑ **Two statistical defects of the estimator.**  A bin populated by one
  spin state had variance zero (R = w_f exactly) and an infinite weight —
  the "rank-deficient" |t| bins of the earlier Table 3; the variance now
  uses expected per-state counts.  The ratio inversion's error Jacobian
  dropped its second-order denominator (0.3%).
- ☑ **Smaller consistency fixes.**  Hadron acceptance in the detector
  frame (+0.7–1.5% of Σ_h); the library's target-mass term carried onto
  the pseudo-events (2% of Σ_h at y = 0.01); p + n merged by cross section
  rather than count; the beam-energy spread with the analysis sign on
  Q²_e against 1 − y_e; the physical nuclear mass in `beam_fourvectors`;
  the two legacy high-divergence constants unified (164 μrad); the 5R
  panel bins selected on measurable criteria only; the LO α_s fallback
  announced.
- ☑ **The coherent closure re-derived at the tagging optics** of Report 1
  §6.1 (`--optics tagging`: 0.33 × 3.8 mrad at 5 × 40.8, σ_θ = 33/380 μrad,
  6×10⁶ response recoils), at all three configurations, with a
  20-experiment ensemble: unbiased in every |t| bin with > 10⁵ recoils,
  biased low by the bin-wise Poisson ratio below that (−5%, −37% at one
  year; closes at ten).  The "residual template bias" on a_e was a 2σ
  fluctuation.  Systematics with exact counts: a 10⁻³ cutout change moves
  a_t by ≤ 0.8% (it was +19% under the slot), δu₂ = 0.024 moves a_e by
  7–9% (was 20%).  The per-fill perturbation now acts on the binding cut.
- ☑ **Report 2 rewritten** as a paper: the chain stage by stage with the
  truth/measured boundary explicit, the assumptions in one table with
  provenance and leverage, the methodology, the closure, the measured
  systematics and the specifications, what is not modelled.  Report 1's
  Table 1 and §7 restated; Report 4 carries a scope banner (its
  acceptances are on the legacy divergence; plans/10 A4 stays open).
- ☑ 16 new tests (225 evgen); the reproduction manual, the READMEs and
  the sample manifest restated; `plans/08` D3 corrected (the mixed x is
  ISR-invariant, the Q² label migrates).

What is still open from the review: a likelihood (or coarser-binned) fit
for the low-count coherent bins *(closed in run 13)*; Report 4's absolute acceptances on the
Yellow Report divergences (plans/10 A4); the ePIC inputs of plans/04 #20
and #21.

## Development run 10 (2026-08-27): the P_zz propagation was 2x wrong, and T1 gets a value

Report 2's §8 listed the two specifications that are ours to set but left
one of them without a number and got both of their propagations wrong.
Both are now measured rather than argued, and pinned.

- ☑ **A P_zz scale error propagates 1:1, not quadratically.**  Commit
  `f05d026` published (1+d)² − 1 = **4.0 / 10.3 / 21.0%** at
  δP_zz/P_zz = 2 / 5 / 10%.  The correct costs are **2.0 / 5.0 / 10.0%**.
  The estimator's weights wᶠ = Pᶠ − P̄ are built from the *assumed*
  polarizations, so the ratio R carries one power of the assumed scale
  while σₚ² carries two, and one cancels: **Â/A = P_zz(true)/P_zz(assumed)
  exactly** — measured 1.020000 / 1.050000 / 1.100000 / 1.300000 at
  d = 2 / 5 / 10 / 30%, on both aₑ and aₜ.  The quadratic is real but
  belongs to *reach*: δA ∝ 1/σₚ, so luminosity at fixed δA goes as 1/P_zz².
  The test that was meant to pin this checked only that σₚ² ∝ (1+d)² —
  true, irrelevant, and it never ran the estimator on rescaled data.
- ☑ **The relative-luminosity coefficient is 1/3 per unit ratio error**,
  not the ≈1.4 × δ that was published.  It is the analytic
  −(P₁+P₂)/(P₁−P₂), exactly 1/3 for the flip plan's (+0.6, −1.2), and
  measured 0.333 out to δ = 10⁻².  The old number conflated two
  conventions: the scripts' `--rel-lumi-offset d` sets the assumed shares
  to [0.5(1±d)] against equal truth, which is a *ratio* error of ≈2d, so
  that convention reads ≈2/3.  Both tests now run the estimator, and both
  conventions are named where the number is quoted.
- ☑ **T1 has a value: δP_zz/P_zz ≤ 5%** (new §8.4).  No EIC document
  states one and **nothing exists for ⁶Li at all** — EPIOS's proposed Li–Li
  CNI polarimeter is explicitly vector, and it puts tensor polarimetry on
  Lamb-shift or BRP-type source-level devices.  So the target is set from
  precedent, with stored beams separated from targets: **2.1%** (JINR
  Nuclotron, absolutely calibrated on ¹²C(d,α)¹⁰B* at 0°) and **≈4%**
  (COSY/ANKE, excluding the EDDA absolute scale) on stored deuterons;
  **3.0%** on the HERMES storage cell, the best tensor number in the
  literature, whose systematic is the *cell* — a mechanism a stored beam
  does not have; 4.9% (NIKHEF), 8.0% (VEPP-3), 4.7–9.7% (UVa solid) and
  8.0% (the JLab b₁ budget) elsewhere.  5% is the working value, 3% the
  optimistic case two experiments have reached, 8–10% the conservative one.
- ☑ **Relative luminosity ≤ 10⁻³ on the ratio**, costing 0.03%.  It is not
  a leading systematic at any plausible value, and bunch-by-bunch tensor
  alternation — what makes RHIC's 10⁻⁴ achievable for helicity — removes
  it.  That is a design choice to state early, not a number to measure late.
- ☑ **Two new DRIFT checks** (22 total) refuse the 2×-pessimistic table and
  an unqualified 1.4 anywhere in `reports/`, `plans/` or `docs/`.  266
  tests and 22 checks pass.

The consequence for what the programme claims is unchanged and worth
repeating: a P_zz scale error multiplies every bin identically, so **every
shape claim is immune to it** — the x and Q² dependence, the sign, the
ratio between sweet spots, and the discrimination between the two Δ ansätze
all survive a wrong P_zz scale untouched.  It moves only the absolute
magnitude of Δ.

## Development run 9 (2026-08-27): Reports 0 and 1 rewritten as papers

The two circulate-able physics reports were rewritten to a publishable
state — abstract, author line, one logical thread, references numbered
in citation order, no repeated facts — and every number in them restated
from the current scripts at the γ-matched energies (plans/10).  Doing so
found that the report text had only been search-replaced when the
energies changed: the figures were regenerated at 99.5 GeV/u on the
morning of 2026-08-27 but the prose still carried the 50 GeV/u numbers.

- ☑ **Report 0** (`polarized_li_primer`): the physics case as a paper —
  atom → nucleus → ion; polarization → collider; formalism → observables →
  a single requirements-and-reach table.  Reach numbers from
  `tagging_acceptance.py`, `money_b1.py`, `money_polemc.py`,
  `money_delta.py` at the current energies.
- ☑ **Report 1** (`cos2phi_money_plots_report`): the projections as a
  paper.  Sweet spots now x = 0.011–0.14 (were 0.02–0.14), amplitudes
  (0.44–0.95)×10⁻² against δA = (1.4–4.5)×10⁻⁴ (21–44σ), xΔ to 2–9%;
  reconstructed δÂ = (0.9–3.0)×10⁻⁴ with purities 0.43–0.69 on the PYTHIA
  final state.  The coherent channel is presented as the far-forward
  requirement it is: with the Yellow Report divergences and the measured
  aperture no published IP6 optics tags the recoil (Table 2 of the
  paper), the figures stand for a 0.20 GeV envelope (1.7×10⁷ tagged in
  year 1), and the analysis closure is shown at a 0.73 mrad approach at
  5 × 41 (`nearbeam_reach_gain.py`).  The R statement corrected: every
  projection uses the toy R; R1998 sits behind the hook and would move
  Δ/F₁ by +17.8 / +18.5 / +7.6 / −4.4% at the current spots.
- ☑ **`hfs_resolution.py` selected the sweet spots from a constant** at
  the pre-correction energies; it now takes them per configuration from
  the money-plot-5 selection.  Σ-method δy/y at the mid spots is
  0.32 / 0.22 / 0.29 / 0.15 (low 0.39 / 0.23 / 0.24 / 0.11, top
  0.23 / 0.19 / 0.21 / 0.18, each at its own spots; the low-configuration
  first spot re-read 0.39 on 2026-08-29, not the 0.38 first logged).
- ☑ **Report 2 patched, not rewritten**: Table 2 (5R) and the text on the
  current spots; money plot 6R re-derived at the low configuration with a
  0.727 mrad near-beam approach (acceptance 7.8%, N_tag 3.4×10⁶, three
  |t| bins, a_t recovered to ±0.003–0.008; a_e carries a residual
  template bias of the order of the one-year error at this aperture —
  0.0151 ± 0.0014 and 0.0102 ± 0.0014 in two pseudo-experiments against
  0.010), since with the measured aperture nothing survives at any
  configuration.  README, evgen/README, the reproduction manual §7 and
  fastsim/README carry the same numbers; `money_delta_realistic.py`
  still runs its own superseded (pre-2026-08-27, rigidity-scaled)
  27.5 / 50 / 137.5 GeV/u configurations and is flagged as such rather
  than changed.
- ☑ **The tagging optics is priced** (`evgen/scripts/tagging_optics.py`,  Report 1 Figure 4, same day, at the user's request).  σ_θ ∝ 1/√β* but
  L ∝ 1/β* per plane, and a recoil escapes through the horizontal pot gap
  OR the vertical one, so only the horizontal β* needs de-squeezing (the
  adversarial review's point): the yield ε × L peaks at
  β*_x/β*_x,HA = 50 / 180 / 90 with the vertical plane at high acceptance
  (planar pots, rectangle envelope 10σ_h × 10σ_v with dispersion at the
  pots), ε = 0.42 / 0.32 / 0.33, L/L_HA = 1/7 / 1/13 / 1/9.5, horizontal
  envelope 0.33 / 0.17 / 0.12 mrad *(0.37 / 0.25 / 0.33, 1/6.8 / 1/12.8 /
  1/9.5 and 0.36 / 0.19 / 0.12 mrad on run 15's per-configuration
  transport)*: 2.6×10⁶ / 3.0×10⁶ / 6.1×10⁶ tagged
  events per year at the 10 fb⁻¹/u placeholder, 2–6× below the 0.20 GeV
  reference; the shape term in the optics' own window (3.1–3.5% per unit
  P_zz) is a 9 / 8 / 11σ/yr measurement and a 1% exotic-glue term reaches
  5σ in 2.8 / 4.4 / 2.6 years.  De-squeezing both planes (the naive
  reading) gives a fifth of that at 1/24–1/70; with the pots fixed at the
  measured aperture nothing is recovered at any β*.  IR-8's secondary
  focus (≈ 20%, our interpolation) is still worth 3 / 8 / 6× the IP6
  optimum if the second IR delivers the same luminosity.  Assumptions
  stated in the paper: the electron β* raised in step; a parallel-to-point
  far-forward transport for the de-squeezed lattice.  The proposed IP6
  configuration (18 × 137.5, β*_x × 90, pots at 0.12 mrad; or 5 × 40.8,
  × 50, 0.33 mrad; a separate running configuration) is stated in Report 1
  §6.1 with that cost.
- ☑ **Is the hadronic final state in acceptance at the sweet spots?**
  (`evgen/scripts/hfs_acceptance.py`, Report 2 §3 Figure 4 / Table 1b,
  same day, at the user's question.)  Mostly: the hadronic system sits at
  η_h = 0.7–2.1, but at y ≈ 0.01 it is a W ≈ 6 GeV, five-charged-particle
  system whose target-fragmentation side reaches the calorimeter edge —
  78 / 86 / 82 / 91% of Σ_h is captured at the four mid spots, 19 / 8 / 16 /
  7% escapes forward beyond |η| = 3.7 (geometry and thresholds; 69 / 74 /
  73 / 85% captured through the full response) — lab-frame figures,
  superseded 2026-08-28 by the detector frame's 80 / 87 / 83 / 92, 17 / 8 /
  15 / 7 and 70 / 74 / 74 / 85% (evgen/README) — 1–6% is below threshold,
  nothing is lost backward, and the lithium fragments carry
  E − p_z = m²/(2p) — 4.4 MeV per nucleon, 26 MeV for an intact recoil — at
  η ≈ 8, never in the measured Σ_h whether tagged or not.  The escape is a
  13–28% scale bias on y_Σ: the library *reproduces* it in the
  pseudo-events (the adversarial review's catch — the reports had said
  "corrects"), the bin-centering factor absorbs it, and an analysis would
  calibrate it — `HFSResponse(calibrate=True)` now does, per cell, taking
  the 5R PYTHIA purities from 0.43 / 0.54 / 0.47 / 0.69 to
  0.52 / 0.59 / 0.59 / 0.76 and D from 0.79–0.95 to 0.95–1.03 at unchanged
  errors; a residual 1% scale error moves Δ̂ by 0.2–0.7%.  It is not the
  resolution driver: the ePIC nominal reach of 4.0 recovers a quarter of it
  (19 → 14% in the lab frame of this entry) with δy/y unchanged (0.32 / 0.22 / 0.29 / 0.14) — the width is
  the noise at the y ≈ 0.01 spots and the within-acceptance capture
  fluctuation at y ≈ 0.025.  Report 1 §5.2 and Table 1 now quote the
  calibrated numbers with the uncalibrated ones alongside; the purity loss
  is the partial, fluctuating capture plus the noise floor, not "neutral
  hadrons" (the 0.11 neutral-hadron share was a pre-correction number; it
  is 0.09–0.10 of Σ as HCal objects including untracked charged particles).
- ☐ Left: Report 2 as prose is still the 2026-08-24 analysis note with
  patched numbers, not a paper; the coherent a_e template bias at a
  strongly anisotropic aperture is unexplained.

## Development run 8 (2026-08-26): the generators arrive, and the aperture is measured

The work moved to the Linux box, and three of the things the plans had
filed as *external* turned out not to be.  PYTHIA 8 builds its own Python
bindings against the analysis interpreter in two and a half minutes; the ePIC
far-forward geometry will accept an intact ⁶Li if you hand it a HepMC
record instead of asking the particle gun for a nucleus; and the EPPS21,
CT18 and NNPDFpol grids install with one `pip` command, which is what had
been keeping the July production un-runnable.  The remaining tail of
[08_simulation_chain_completion.md](08_simulation_chain_completion.md)
(A6, B4, C1, C2, C3) was closed in parallel, each item implemented, then
handed to an adversarial reviewer, then repaired.  Both suites are green
at **48 fastsim + 183 evgen**, from 25 + 143.

- ☑ **The hadronic final state is PYTHIA 8, not a toy** (D4).  8 M
  events over the three beam configurations, 2.7 GB, twelve minutes of
  wall time (`tools/pythia8`, manifest in `evgen/samples/README.md`).
  Building it found a defect in the settings this repository had
  documented for two months: PYTHIA applies `PhaseSpace:Q2Min` **only**
  when `Q2Min ≥ pTHatMinDiverge²`, and that parameter defaults to 1 GeV,
  so the requested 0.7 GeV² was silently ignored — min Q² = 1.002 and
  σ = 0.381 μb whatever was asked for.  With `pTHatMinDiverge = 0.5`:
  min Q² = 0.697 and σ = 0.551 μb.  The missing 0.7–1.0 GeV² band is 31%
  of the sample and is exactly the band the loosened generator window
  exists to populate.
- ☑ **The hadronic resolution is measured, and the toy was optimistic.**
  Σ-method δy/y at the four sweet spots, 50 MeV of calorimeter noise:
  **0.55 / 0.28 / 0.50 / 0.15** at the mid configuration against the
  toy's 0.28 / 0.17 / 0.24 / 0.07 — optimistic by 0.04–0.05 absolute at
  every spot, because the toy put 0.60 of Σ in tracks and 0.03 in neutral
  hadrons where PYTHIA puts 0.51 and 0.11, and neutral hadrons are the
  HCal's problem.  *(Superseded 2026-08-27, run 9: at the corrected spots
  HCal objects — neutral hadrons plus the charged particles the tracker does
  not see — carry 0.09–0.10 of Σ within acceptance, and the purity loss is
  the partial, uncalibrated capture of Σ_h, not neutral hadrons — Report 2
  §3 Figure 4 and Table 1b.)*  It degrades monotonically with beam energy
  (0.28 / 0.21 / 0.24 / 0.11 low, 0.74 / 0.34 / 0.69 / 0.18 top), which
  settles the x ≈ 0.1 bins on the **low-energy** configuration — the open
  item the reconstruction note lists.  Sweet-spot purity falls from
  64–68% to 40–73% and the reco-bin amplitude sits 4–21% below the
  true-bin value instead of 1–9%: the response is a bigger correction
  than the stand-in implied.
- ☑ **The Roman-Pot aperture for an intact ⁶Li is measured, and it opens
  horizontally.**  npsim cannot shoot a nucleus (`Bad particle type` for
  every spelling), so `tools/fullsim/ion_gun_hepmc.py` feeds one through
  as HepMC and scans p_T *and* azimuth.  The boundary is
  |θ_x| ≳ 2.0 / 1.35 / 1.03 mrad in the 5×41 / 10×100 / 18×275 optics
  against |θ_y| ≳ 1.8–3 mrad; the ⁶Li and an α at the same rigidity agree
  at the same *angle*, which is the check that this is optics and not
  species.  The cause is the transport, not the pot orientation: the
  horizontal lever is R₁₂ ≈ 30.6 m against a vertical few metres.
  *(Superseded 2026-08-28, run 14: in the current geometry the edges are
  2.50 / 1.51 / 0.53 mrad against (shut) / 2.12 / 0.92 mrad, R₁₂ =
  19.24 / 21.25 / 29.97 m against R₃₄ = — / 3.35 / 2.93 m, and the ratio
  to the 10σ envelope is 1.14× / 0.84× / 0.58× — the silicon binds at
  5 × 41; the chain numbers below were computed against the
  September-2024 aperture.)*  This
  **inverts** `rp_measure`'s assumed slot (`cut_scale_xy = (2.5, 1)`
  against a measured ≈ (1, 2.3)) and therefore flips the sign of the
  acceptance-induced ⟨cos 2β⟩.  Carried through the chain the same day:
  at the low configuration the coherent measurement survives — acceptance
  37.7% → 1.42%, the fake ⟨cos 2β⟩ **+0.426 → −0.772**, two of the four
  |t| bins, δa_t worse by 6–34× and a_e still recovered — and at mid and
  top the aperture leaves nothing in the binned window, so the coherent
  programme is a low-energy programme.  The measured edge is
  2.8× / 1.9× / 1.4× the 10σ envelope, so the envelope is never the
  binding constraint; it is marked on all three WP5 curves.  *[Both
  statements superseded 2026-08-28 (run 13): against the per-configuration
  Yellow Report envelopes the edge is 0.91× / 0.75× / 1.12× the horizontal
  10σ half-width, binding at the two lower configurations, and at the
  tagging optics the top configuration is the best covered — plans/04
  #20, `tools/fullsim/README.md`.]*  Two latent
  defects surfaced with it — a tight cutout made the two-azimuth fit
  raise a bare "Singular matrix", and one dead |t| bin aborted the whole
  figure.  Caveats and what is left in plans/08 §8.4 and plans/04 #20.
- ☑ **The e+d control, second pass** (plans/02 step 1.5.3), on the
  official BeAGLE 1.03.02-3.1 sample over xrootd.  Routing agrees to
  better than two points; the p_T tail does not, and restricting to the
  spectator peak makes it worse: BeAGLE is 2.5× / 6.9× / 28× the Hulthén
  model at p_T > 0.2 / 0.3 / 0.45 GeV.  **No β reproduces the shape** —
  BeAGLE has a narrower core and a harder tail than a two-parameter
  Hulthén can have at once, and the best fit sits at β = 0.40 GeV, the
  top of the documented scan range.  Since the ⁶Li α tag is entirely a
  p_T-tail measurement, its model uncertainty is one-sided **upward**.
- ☑ **The bin-centering factor has a fitted alternative** (A6).  The
  amplitude is *exactly* linear in Δ, so `RecoResponse.fold` is the
  response as an exact linear operator and `fold_shape_fit` fits a Δ(x)
  shape through it per Q² slice on a bounded grid — bounded because that,
  not a matrix inverse, is the conditioning guard.  Correcting with the
  wrong prior, the residual bias over the plotted bins falls from
  22% / 24% / 99% to at most 3.5% / 7.1% / 6.9%, and the sweet-spot model
  dependence from (−5.0, +8.5, −9.3, +6.3)% to (−0.9, −1.2, −0.2, +0.3)%.
  The 7R bar now carries four terms: statistics, the shape fit, the
  response MC (0.10–1.00% per plotted bin) and the spread over the priors
  the tilt family cannot absorb — the last of which exists only because
  the adversarial pass caught the bar being blind to the biggest of them.
  `--unfold model` stays the default and reproduces every published 7R
  number bit for bit.
- ☑ **R = σ_L/σ_T is the published fit, behind one hook** (C2).  The
  E143 R1998 world fit is transcribed from the paper (all three forms,
  refs/hep-ex_9808028.pdf), and `r_func` threads it through every
  consumer — including `polarized.ToyG1._f1` and `fom.py`, which the
  dated scripts' monkey-patch misses.  The defect it replaces clipped R
  to exactly 1.000 over 38% of the sensitivity box.  Measured
  consequence at the four sweet spots: **F₂ cancels exactly** (2×10⁻¹⁶,
  checked with two independent levers) while Δ/F₁ moves
  **+16.6 / +18.0 / +4.7 / −4.4%**.  One correction to the code review:
  its five R values are a property of the script's mongrel function, not
  of R1998, and should not be quoted.
- ☑ **The July production runs, and its reach was half of what it should
  be** (C3).  It reproduces the dated note to all four recorded digits
  (A_bag −0.3178 / −0.3100 / −0.2967; the 0.07% residual is the |A_bag|
  provenance the fix note itself flags).  With the published R, at
  Δ/F₁ = 10⁻³: **L_5σ = 67.5 / 65.8 / 155.1 fb⁻¹/u** for LOW / MID / TOP
  against the recorded 135.3 / 131.3 / 274.6 — the code review predicted
  63–69 and 152–164 and both land inside.  |A_bag| moves 21–25%, not the
  ~10% the review estimated.  The frozen numbers stay reproducible with
  `--r-model simplified`; the six July notes carry a superseded banner.
- ☑ **Nuclear masses and the last untested path** (C1, B4).  A·M_U is not
  the mass of a nucleus; the beam boost now uses AME2020-derived nuclear
  masses in both `spectator.py` and `tagged.py` (the same beam had two
  masses in two modules).  The kinematics move 2×10⁻³, the ⁶Li α tag
  moves 1.65% → **1.85%** because the R = 0.95 window edge is hard.
  `helicity_flip_plan(...).pzz_true` — the quantity the money scripts
  divide by — is pinned for every branch.
- ☑ **A reproduction manual, and the whole dated line runs.**
  [../docs/reproduction_manual.md](../docs/reproduction_manual.md) is the
  map from result to command: environment, every script with its measured
  runtime, the numbers to check against, the third-party generators, and
  what cannot be reproduced here and why.  Writing it meant running
  everything, which found the last two things that did not: the
  2026-07-20 and -07-21 scripts want a fourth PDF grid
  (`nNNPDF30_nlo_as_0118_A6_Z3`) and exited 2 without it, and
  `money_delta_20260724.py` still carried the bare `np.trapezoid` that
  code review S11 flagged — the same one-line shim as its successor fixes
  it, and renaming a function moves no number.  **All seven dated scripts
  and all 30 scripts in the repository now run**, about eleven minutes for
  the lot.
- ☑ **Housekeeping the box unlocked**: the CT18NLO / EPPS21nlo_CT18Anlo_Li6
  / NNPDFpol11_100 grids are installed (the fast-sim grid tests no longer
  skip), and `reports/build_report.py` renders the PDFs on Linux through
  a user-cache chromium, so the whole report pipeline runs here.

## Development run 7 (2026-08-25): the systematics the estimator cannot cancel

The reconstruction-chain note was rewritten as current content: the
archaeology of the code went (§5 and §6 now state today's status, the old
§5.1 implementation recipe is gone, and the annotated recommendations are
reduced to what is open), which cut the prose by 11%; the results below
then put half of that back, so the note is 3% shorter and says
considerably more.  The gaps it left were audited and closed.  The audit
itself is [08_simulation_chain_completion.md](08_simulation_chain_completion.md):
six parallel passes proposed 67 gaps, an adversarial reviewer refuted 44
against the code, and the 23 that survived are the plan.  Nine of them are
now done; the convention items that would change a documented choice are
left as author decisions (D1, D7, D9).

- ☑ **The one systematic the spin-state ratio cannot cancel is modelled**
  (F1).  `phi_eff` accepts one efficiency per fill and
  `reco.fill_acceptance_bias` gives the analytic bias
  Σᶠ lᶠ (Pᶠ − P̄) eᶠ / σₚ²; for two fills it reduces to
  (e₊ − e₀)/(P₊ − P₀) at **any** luminosity split, so a lopsided run plan
  buys no protection.  A 10⁻³ difference of the cos 2φ′ harmonic between
  the samples fakes 5.6×10⁻⁴ — 5% of the sweet-spot amplitudes, but 4.6
  one-year statistical errors.
- ☑ **Two documented numbers were wrong and are corrected against the
  code.**  (i) A 1% fill-to-fill change of the Roman-Pot vertical envelope
  was estimated at 1.3% on aₜ; the template fit gives **+169%**, and 10⁻³
  already gives +19%.  The amplification is in the fit — the slot leaves
  12 of 24 β bins live, where the t-template is 99% anti-correlated with
  the constant — so the requirement is bunch-by-bunch alternation or
  ≈10⁻⁴ envelope stability, and the higher-|t| bins are the fallback
  (+3.9%, +0.04%).  (ii) The unpolarized u₂ was described as a
  second-order nuisance; an error in it reaches aₑ at **first** order as
  aₜ·δu₂·⟨cos 2β⟩, and a ZEUS-1σ δu₂ moves aₑ by 20%.
- ☑ **The coherent channel has a null test.**  sin 2α, sin 2β and
  sin(α+β) are exactly forbidden by reflection symmetry; fitted alongside,
  a spin-axis error gives tan 2δ in both ratios while a Roman-Pot
  azimuthal roll gives it in β alone — the recoil null resolves a pot roll
  to 5–8 mrad in year one.  The columns are orthogonal to the cos ones, so
  Table 3 is reproduced digit for digit with them on.
- ☑ **WP5 exists as a curve.**  `coherent_optics_scan.py`: tagged fraction
  32% / 3.0% / 4×10⁻⁵ / 2×10⁻⁷ and δaₜ/aₜ = 1.6% / 5.7% / 104% / 540% at
  a near-beam envelope of 0.10 / 0.22 / 0.45 / 0.60 GeV (re-run 2026-08-29
  on the run-14 chain: 1.2% / 4.6% / 79% / 392%, the tagged fractions
  unchanged) — **the coherent
  measurement lives at the low- and mid-energy configurations and is dead
  at the top energy.**  Needed importance sampling above the cut: the
  plain sampler leaves zero accepted recoils above 0.3 GeV.
- ☑ **The far-forward near-beam cut is angular** (S8).  0.20 GeV on a
  4 × 137.5 GeV α is 5σ, not 10σ: the ⁶Li α tag falls to 1.7% / 1.3%
  from the 3–9% the README quoted (which already flagged it).  ⁷Li is off
  rigidity and untouched at 96–98%.  OMD window → ζ = 0.45–0.65 (S7).
- ☑ **Detector nuisances, measured** rather than argued
  (`money_cos2phi_reco.py --syst-scan`, MC noise floor 0.13–0.21%):
  electron energy scale ±1% → 0.2–1.4% on Δ̂; hadronic-resolution
  mismatch (generate 0.30, correct 0.25) → 0.5–1.3%; ε_eID η tilt → 0.02%;
  the Yellow Report EMCal η table → 0.000% at spots 1–3 and 0.5% at spot 4.
  All below the 3–11% model dependence of K.  A correction to the audit
  that proposed this: the electron lever is 2 − y only for the Σ method —
  the Gaussian y stand-in never sees E′ and gives exactly 1.
- ☑ **The tensor sign has an external anchor** (G1).  The suite compared
  code to code, so flipping the sign in both files left all 125 tests
  green.  `asymmetries.TENSOR_LL_SIGN` is now the single constant and
  A_zz(θ_S = 0)(1 + ε(y)R) = SIGN·(2/3)b₁/F₁ — exact, F₂-free — is pinned
  against Cosyn Eq. (27).  Verified non-blind: flipping it fails only that
  file.  The decision itself stays the author's (D1).
- ☑ **Kernel hygiene**: an exact positivity guard on the φ density (G3 —
  the accept–reject silently sampled max(W, 0)); one rank-2 geometry for
  J = 1 and J = 3/2 (the spin-3/2 branch's rate and cos 2φ channels
  disagreed with each other by 3; latent); coverage for the
  spin-temperature populations, a non-default b₂, and θ_S between 0 and
  π/2.
- ☐ **Left for the author** (plans/08 §8.2): the b₁ sign flip (D1, now one
  line), the ⁶Li effective-polarization convention and the double
  dilution it produces in g₁ (D7 — the structural fix is unambiguous, the
  *value* is plans/04 #6), and the b₁ money-plot dilution (D9).  Also
  deliberately not unified: the high-divergence envelope, quoted at 0.41
  GeV in `polligen.reco` and the rounded-up 0.45 in the fast sim, where
  every published fast-sim number comes from.
- ☐ **Still blocked externally** (plans/08 §8.3): the ePIC calorimeter
  noise floor at Σ_h ≈ 0.2–0.5 GeV (#21, the one number the letter cannot
  do without), the backward-disk angular resolution, PYTHIA 8 samples,
  BeAGLE breakup shapes, a coherent diffractive model for ⁶Li.

Tests: 143 evgen + 25 fastsim (from 94 + 24 at the start of run 6).

## Development run 6 addendum (2026-08-25): placeholders filled from the references in refs/

- ☑ Seven papers added by the user to `refs/` (index: `refs/README.md`;
  PDFs tracked in git since 2026-08-25) read and folded into the reconstruction-chain report,
  the code defaults and plans/04: **ZEUS LPS** (NPB 816:1) bounds the
  diffractive azimuthal harmonics u₁, u₂ (A_LT, A_TT consistent with zero
  within ±0.03–0.05) and documents a beam pT spread of 45 (x) vs 100 (y)
  MeV at the HERA IP, Φ resolution 0.2 rad, σ(t)/t = 0.14/√|t|; the
  **Nikolaev–Pronyaev–Zakharov / Pronyaev** papers give the y-factors and
  the LT/T model (A_LT ≈ 0.03 at our β ≈ 0.5); **Cosyn et al.** EPJ A 61:83
  give the exact spin-1 harmonics — b₂–b₄ feed cos φ_TL at O(γ) and
  cos 2φ_TT at O(γ²) (Eqs. 17d–e): ≤ 3% of the Δ amplitude at the worst
  sweet spot, computable from A_zz; **Jentsch DIS 2023** fixes the RP
  geometry (planes around a horizontal slot → slot-like cutout, new
  default of money plot 6R: acceptance 3.2%, N_tag 2.7×10⁶, fake ⟨cos 2β⟩
  = +0.77, deformation-term error ×3 vs a square cutout, a_e unaffected)
  and the beam-effect-dominated RP resolution (ΔpT ≈ 40 MeV at 275 GeV,
  detector ≤ 1.5%); **Mäntysaari et al.** full text confirms Eq. (9) and
  the γ*d-frame definition of Φ. 94 evgen tests.
- ☑ **Reference dictionary + web search for the missing item** (same
  day): `refs/refs_dict.json` (machine-readable: identifiers, key content
  with equation/figure/slide pointers, where used) + `refs/find_ref.py`
  lookup; the ePIC inclusive-DIS reconstruction performance was found in
  S. Maple's Dec-2024 ePIC seminar (σ(δ_h)/δ_h = 25% smearing; Δy/y per
  method in y bins; 0.01 < y < 0.95 coverage rule; DA best at y ≈ 0.01),
  the ATHENA proposal Sec. 3.1/Fig. 22 (e−Σ/DA for y ≲ 0.1; ≈ 25% at
  y ≈ 0.01, ≈ 10% at y ≈ 0.1; JB 20–30%) and Arratia et al. NIM A
  1025:166164 Fig. 5 (IΣ/DA/JB 13–17% at y ≈ 0.05–0.2) → the reco-level
  default hadronic-y resolution is now 25% (money plots 5R/7R rerun);
  #21 reduced to "no published ePIC full-sim number at Q² ≈ 1–3 GeV²,
  y ≈ 0.01, for light ions".
- ☑ **Reference copies committed** (2026-08-25, second pass): every
  dictionary entry with a free copy is now in `refs/` (13 arXiv PDFs incl.
  the Yellow Report, split into four < 50 MB parts, plus the TUNL A = 6
  evaluation from nucldata.tunl.duke.edu) and the PDFs are tracked in git
  at the user's request (≈ 240 MB; `refs/find_ref.py --check` / `--fetch`).
  First pages and the cited numbers were verified against the dictionary:
  Jentsch–Tu–Weiss Table I, the IR-8 efficiencies (17.75% ⁷Li, no ⁶Li),
  the EPIOS G-factors, the TUNL widths (24 keV, 541 keV) all confirmed;
  two citations corrected — arXiv:2509.18558 is Hamwi–Hoffstaetter,
  *Polarization transmission in the EIC's Hadron Storage Ring* (Devlin is
  acknowledged, not an author; Table I with 81/573 resonances confirmed),
  and arXiv:2603.23699/23700 are *SIDIS on a polarized spin-1 target*,
  Parts I/II (the “order unity at ≳ 300 MeV” statement is verbatim in the
  Part II abstract). Only HJM/JM 1989, Sather–Schmidt, Jacquet–Blondel
  1979 and Li–Sick 1971 have no free copy.
- ☑ **Simulation code review and reconstruction audit** (2026-08-25,
  `docs/code_review_2026-08-25.md`): whole code base read under three
  lenses (consistency with the references, truth-information use,
  measurability with the ePIC design), three independent reviewers plus
  numerical checks. Verdict: no truth-only quantity enters the 5R/7R/6R
  selections, azimuths, binning or fits (truth only in the response and
  the closure references). Findings: the spin-state ratio cancels only an
  acceptance common to both spin states — a 10⁻³ fill-to-fill
  difference of the cos 2φ′ efficiency harmonic fakes half the signal,
  so bunch-by-bunch alternation is a requirement (plans/04 #3); the
  hadronic-y stand-in is an absolute 50–255 MeV (50–125 MeV at spots 1–3) resolution on
  Σ_h = 0.2–0.5 GeV at the sweet spots (no published ePIC number, #21);
  report §3 numbers corrected (δQ²/Q² = 5% with the unsourced 3 mrad
  angular table, φ′ dilution 0.99); K model dependence 3–11%;
  EMCal term applied at all η. Outside the chain: the b₁-sector sign is
  opposite to HJM/HERMES/Cosyn Eq. (27) (A_zz = −(2/3)b₁/F₁; Δ sector
  untouched) and the fast-sim `r1998` returns R = 1 at low x (July
  L_5σ values 1.5–2× too pessimistic; polligen unaffected) — both
  left for an author decision. Documentation-only changes in this pass.
- ◐ **WP3-HFS: hadronic final state and hadron-side detection**
  (2026-08-25, option 1 of the reconstruction-note discussion):
  `polligen/hfs.py` (sample format, exact hadronic sums, Σ/JB/DA/mixed
  methods, hadron-side response with tracker/calorimeter coverage,
  thresholds, efficiencies, resolutions and noise, toy string-fragmentation
  stand-in, (x, Q²) event library → `recopseudo` hook, 7 tests, 101 evgen
  tests); `tools/pythia8/gen_dis_hfs.py` for the PYTHIA 8 sample (eic-shell);
  `scripts/hfs_resolution.py` → Figure 3 of the reconstruction report. Toy
  result: Σ-method δy/y 0.28 / 0.17 / 0.24 / 0.07 at the sweet spots with a
  50 MeV noise floor (9–12% without noise) — the 25% stand-in is the noise
  floor on Σ_h ≈ 0.2 GeV; 5R with the HFS-based y reproduces Table 2's
  errors. Pending: the PYTHIA sample and the ePIC noise/threshold floor.

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
  40.8 / 99.5 / 137.5 GeV/u for the same optics — the low-energy
  configuration is the coherent program's home; the a_n normalization
  convention of the deuteron anchor (factor 2) is flagged for a source
  check.
- ☑ **Second pass (same day): WP3 reconstructed-level pseudo-experiments**
  — `polligen/recopseudo.py` (+11 tests, 94 collected), money plots
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
  (Δ̂ = −Â·y²D_φ/(1−y), model bin-centering, independent 1-/10-yr
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
  fragments R ≈ 1 — 0.998 and 1.005 with the physical nuclear masses,
  run 13 — same velocity; dE/dx Z² separates, and since run 13 so does
  the partner-fragment two-hit veto, plans/09 B4); every
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
  (⁷Li α-tag 96–99%, ⁶Li 3–9%/HA optics — 1.7% since the angular cut of
  2026-08-25, gluonometry L_5σ = 15–22 fb⁻¹/u
  toy, δΔR ≈ 3.5% at x = 0.3); toy-vs-CT18 F2 within ±37% as documented.
  Findings (kept as report, not yet fixed): two inert test assertions
  (`test_smoke.py:17`, `test_spectator.py:68`); FOMs implicitly give each
  observable the full luminosity in its own spin configuration (run-plan
  split not modeled — parameterized and stated in run 14); no γ²/A₂
  target-mass terms at high x, low Q² (bounded and implemented default-off
  in run 14);
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
