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
| [09_nearbeam_nanowire_far_forward.md](09_nearbeam_nanowire_far_forward.md) | A near-beam layer for the far-forward lithium tags (2026-08-26): what a closer approach is worth (×26 / ×569 / still-dead per optics, and mid-energy going from unusable to the best-covered configuration), whether a superconducting nanowire can deliver it, the hot-spot firing-threshold answer to open question #19, the obstacle table, and the correction that the ePIC pot geometry has moved since the snapshot `tools/fullsim` measured — report in `reports/nanowire_far_forward.html` |
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
  0.32 / 0.22 / 0.29 / 0.15 (low 0.38 / 0.23 / 0.24 / 0.11, top
  0.23 / 0.19 / 0.21 / 0.18, each at its own spots).
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
  HCal's problem.  It degrades monotonically with beam energy
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
  horizontal lever is R₁₂ ≈ 30.6 m against a vertical few metres.  This
  **inverts** `rp_measure`'s assumed slot (`cut_scale_xy = (2.5, 1)`
  against a measured ≈ (1, 2.3)) and therefore flips the sign of the
  acceptance-induced ⟨cos 2β⟩.  Carried through the chain the same day:
  at the low configuration the coherent measurement survives — acceptance
  37.7% → 1.42%, the fake ⟨cos 2β⟩ **+0.426 → −0.772**, two of the four
  |t| bins, δa_t worse by 6–34× and a_e still recovered — and at mid and
  top the aperture leaves nothing in the binned window, so the coherent
  programme is a low-energy programme.  The measured edge is
  2.8× / 1.9× / 1.4× the 10σ envelope, so the envelope is never the
  binding constraint; it is marked on all three WP5 curves.  Two latent
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
  a near-beam envelope of 0.10 / 0.22 / 0.45 / 0.60 GeV — **the coherent
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
  (⁷Li α-tag 96–99%, ⁶Li 3–9%/HA optics — 1.7% since the angular cut of
  2026-08-25, gluonometry L_5σ = 15–22 fb⁻¹/u
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
