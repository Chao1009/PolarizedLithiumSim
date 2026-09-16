# Open Questions & External Dependencies

Items that gate or shape the simulation program but are *not* solvable
inside it. Each has an owner-to-engage and a default assumption we proceed
with until answered. Updated 2026-09-15 (first version 2026-06-12): run 19 added
items 25–29, the five author decisions its measurements opened, and closed
the definition half of #17 and the ask half of #15 against the sibling
generator's delivered work.  The 2026-09-06 revision recorded the
consistency review of Reports 0–4, which added items 22–24, the three
passages it could not settle from the repository. The 2026-08-25 revision followed the fetch-verified
literature sweep — several items moved from "unknown" to "answered, needs
adoption".

## Answered by the literature sweep (adopt, then verify with owners)

A1. **Ring spin dynamics for Li** — EPIOS white paper (arXiv:2510.10794,
    PRC 113:060501, Table 3): G(⁶Li) = −0.178 with 81 linear resonances to
    top energy, handled like the deuteron (15 Tm solenoid partial snake +
    jump quads; AGS crossing only imperfection resonances); G(⁷Li) = +1.532,
    573 resonances, partial snakes keeping spin tune 0.2–0.8 + jump quads.
    Top energies ~138 (⁶Li) / ~117 (⁷Li) GeV/u. Deuteron antecedent: Huang
    et al. PRAB 23:021001. *Remaining question:* polarization-survival
    fraction through the full chain (not in any paper) → keep FOMs
    parameterized in P.
A2. **Why ³He is quoted at 166 GeV/u in old documents** — eRHIC-era
    250-GeV-proton rigidity legacy (Milner arXiv:1809.05626). Current design
    number is 183 GeV/u = ⅔ × 275 with 6-snake spin preservation (CDR
    Sec. 5.5; eic.jlab.org/Requirements). Li rigidity scaling (137.5/117.9
    GeV/u) therefore stands, pending C-AD blessing of Li specifically.
A3. **Far-forward routing of Li fragments** — verified windows (plans/03
    §2.2): ⁷Li α → Roman Pots (R = 0.86); the ⁷Li triton at R = 1.29 →
    **RP-inner + ZDC** (measured 2026-08-28), the ⁶Li ³He+t triton at
    R = 1.50 still **lost**, 152 mm at the pot plane against a 144 mm
    module edge; ⁶Li α/d → beam-blind below RP pT cutoff
    (R = 0.998 and 1.005 from the physical nuclear masses since
    2026-08-28, not the 1.000 the A·Z ratio gives — the two fragments are
    separated by 0.7% of rigidity, which is inside the ±5% near-beam band
    and so still undispersed, but they are not the same trajectory);
    p → OMD; n/γ → ZDC. IR-8 secondary focus (RPs 44–45.5 m) recovers R ≈ 1
    at pT → 0. *Remaining question:* none at concept level — quantitative
    acceptance is exactly Phase-2 step 2.2.  *The tritons, settled
    2026-08-28:* "no coverage" was the routing's own answer — `route_charged`
    carried no R > 1 branch, so an over-rigid fragment was lost by
    construction — and the particle-gun scan of `tools/fullsim`
    contradicted it. `farforward.over_rigid_route` is that scan's answer:
    the pot dispersion carries a fragment that bends less than the beam to
    +x and the ⁷Li triton lands on Roman-Pot silicon at every
    configuration, then in the ZDC, so `route_charged` returns route 6,
    "RP-inner (over-rigid)". The re-measurement in the current geometry
    (git 9aaa2969, 2026-08-28) settles it: dx = +66 mm at station 1 and
    +70 to +72 at station 2, 60 of 60 events at all three configurations,
    with a ZDC deposit in 80–98% of them, and the ⁷Li t tag is
    78 / 92 / 94% against the 0.033 / 0.004 / 0.005 the routing-as-lost
    picture gave. plans/03 §2.2 "Tritons at IP6 — revisit" is closed ☑.
    What the scan still does not carry is a beam envelope and a
    reconstruction, and the triton's own R12 is 8% larger than the beam's,
    which opens an acceptance hole between θx = −1.55 and −2.53 mrad that
    no purely angular routing can see; plans/09 B1 records it.
A4. **BeAGLE status for light nuclei** — runs any (A,Z) but A>4 uses the
    C-12 Fermi-momentum parameterization, Woods–Saxon geometry without
    α+d/α+t clustering, FLUKA evaporation untuned for A<12, code frozen
    since 2023, FLUKA license required (prebuilt at BNL/JLab/CVMFS).
    *Remaining question (Q7 below):* is the collaboration maintaining it /
    can we get light-ion guidance?

## Still open — machine / accelerator

1. **Polarization survival through EBIS charge-breeding + ring for Li.**
   Explicitly "a goal of the study rather than a promised outcome" in the
   ECRP proposal; no number exists anywhere.
   *Engage:* EPIOS (Raparia, Rathmann), MIT ³He group (Milner).
   *Default:* P_z = 0.7 at IP, band {0.5, 0.9}.
2. **Transverse ion polarization at IP6 for Li.** Gluonometry needs
   transverse spin with unpolarized electrons. HSR stable direction is
   vertical (transverse) in the arcs — possibly the *easy* orientation —
   but the IP6 rotator/snake configuration for Li species is undefined.
   Also: clean b₁ extraction prefers polarization along the momentum
   transfer (Cosyn et al. arXiv:2410.12764) — a spin-direction systematic
   to design for.
   *Engage:* C-AD spin group via EPIOS; INT program (Mar 22–Apr 2, 2027).
   *Default:* transverse running available with P = P_z value.
3. **Tensor (λ=0) bunch operations.** Source RF transitions can prepare
   m = 0 (proposal Sec. 3.5); unknowns: survival through acceleration,
   bunch-by-bunch spin patterns, relative-luminosity control at 10⁻⁴.
   *Code review 2026-08-25 (docs/code_review_2026-08-25.md, F1):* the
   spin-state ratio cancels the φ′ acceptance only when both spin states
   see the same acceptance — a 10⁻³ difference of its cos 2φ′ harmonic
   between the m = ±1-rich and m = 0-rich samples fakes
   (ε₂⁺ − ε₂⁰)/(P₊ − P₀) = 5.6×10⁻⁴ — half a Δ/F₁ ~ 10⁻³ signal, 5% of the
   sweet-spot amplitudes, and 4.6 one-year statistical errors (modelled
   since 2026-08-25 by `reco.fill_acceptance_bias` and
   `money_cos2phi_reco.py --eff-cos2-split`; plans/08 A1).  The coherent
   counterpart is a property of the cutout, so it moved when the cutout
   did.  *Dated record, 2026-08-27 (plans/08 A1b, the assumed 2.5 : 1
   slot):* half the β bins were blind there and the t-template was 99%
   anti-correlated with the constant, so a shape perturbation was
   amplified ~100× over δ⟨cos 2β⟩/(P₊ − P₀) = 1.3% — 10⁻³ of the
   Roman-Pot vertical envelope between the samples biased a_t by **19%**
   and 1% by **169%**, the higher-|t| bins being the fallback at +3.9%
   and +0.04%.  That is what made 10⁻⁴ envelope stability the stated
   requirement.  *2026-08-28 (run 11), re-measured at the tagging optics
   of Report 1 §6.1, where the coherent channel is now measured:* a 10⁻³
   change of the binding (horizontal) half-width between the fills moves
   a_t by **−9.1 / −3.2 / −1.8 / −1.1 / −0.5 / −0.5 / −0.5 %** in the
   seven |t| bins of the window adopted 2026-08-28, and 10⁻² by −83.3 /
   −30.1 / −18.0 / −9.8 / −4.9 / −4.5 / −3.5 %, with a_e untouched
   (Report 2 Table 6); the worst bin at the other two configurations is
   −31.2 % (10 × 99.5) and −10.1 % (18 × 137.5), always the lowest.  The
   tagging cutout's horizontal edge sits in a shallow part of the recoil
   spectrum and every bin stays live, so the ×100 amplification is gone,
   but the requirement is now edge-of-window dependent: per-mille
   stability costs 0.2–1.4 % in the four bins above 0.05 GeV² and
   1.4–12.5 % in the three below it, so 10⁻⁴ is what those three would
   need.
   Bunch-by-bunch alternation remains a requirement of the measurement
   — it is what makes the acceptance cancel at all — but fill-by-fill
   running is no longer excluded by this systematic at the tagging
   optics; the φ′ efficiency of the inclusive channel still needs 10⁻⁴.
   *Default:* equal thirds (+,0,−), δ(rel-lumi) as a Phase-2 systematic.
4. **Li luminosity.** Confirmed gap — no Li number exists in any document
   (EPIOS included). Space charge, IBS, cooling for Li bunches unstudied.
   *Engage:* EPIOS/C-AD. *Default:* 10 fb⁻¹/nucleon per setting, quoted
   ∈ {1, 10, 100}.
5. **Li ring polarimetry.** EPIOS concept: Li–Li elastic CNI vs polarized
   Li jet (HJET analog) + Breit–Rabi absolute; analyzing-power theory
   flagged as needing work. R&D scale in EPIOS: ~26 FTE-yr/12 yr.
   *Default:* δP_z/P_z = 3% (vector) systematic in FOMs; the tensor scale is T1, δP_zz/P_zz ≤ 5% with 3% the optimistic case (run 10, plans/00).

## Still open — generator / theory

6. **⁶Li effective-polarization convention** — ☑ **closed 2026-08-29**
   (author decision), on the **cluster picture**.  The ⁶Li spin is carried
   by the α–d relative motion and by the deuteron inside it; the α is
   J = 0 and contributes nothing, so the polarized proton and neutron are
   each polarized along the ⁶Li spin by the product of the two vector
   dilutions, (1 − 1.5 P_D^{α−d})(1 − 1.5 P_D^{d}) = 0.86995 × 0.9325 =
   **0.81123** whole-nucleus, which in the per-nucleon slot convention of
   D7 is P_p = P_n = 0.81123/3 = 0.27041 (`beams.LI6`).  The alternative
   it retires is Cloët's slides' P_p = P_n = 1/3, i.e. a whole-nucleus
   1.0 — 1.233 times as large, and the optimistic end of the pair since
   D7 — which survives as `beams.LI6_NAIVE_ONE_THIRD` and is pinned as
   the pre-2026-08-29 value.  ⁷Li was already settled: P_p = +0.866,
   P_n = −0.037.
   The other reading of the same quantity is ab initio and is not
   retired but not adopted: six-body VMC (Wiringa et al. PRC 89:024305
   Table I, 1.924 spin-up against 1.076 spin-down protons and neutrons
   in M = 1 — the table ⁷Li's own slots are read from, and plans/01's
   verified fact-check entry) gives **0.848** whole-nucleus, 4.5% above
   the cluster product, and read backwards implies an α–d vector factor
   0.848/0.9325 = 0.909 rather than E155's 0.870.  The cluster
   construction is preferred anyway because it shares one wave function
   with the tagged sector instead of transcribing a constant that no
   other observable here would then constrain; **0.81–0.85 is the band**
   recorded with the adopted value (plans/02 step 1.1 item 2), and #15
   below — a VMC α–d overlap — is what would collapse it rather than
   bracket it.
   The decisive property is not the number but where it comes from: the
   two D-state probabilities now live in `polli_fastsim.beams`
   (`P_D_LI6` = 0.0867, `P_D_DEUTERON` = 0.045) and `polligen.tagged`
   re-exports them, so the inclusive effective polarization and the
   tagged α–d S/D interference are the same wave function seen in two
   experiments and cannot drift apart; the deuteron's own slot is the
   expression 1 − 1.5 P_D^{d} verbatim, which makes per-nucleon
   g₁(⁶Li)/g₁(d) = (1 − 1.5 P_D^{α−d})/3 = **0.290** exactly — the
   deuteron's D state cancels between the two isoscalar ions — where the
   naive constant gave 0.358.  What moved with it: g₁(⁶Li) is multiplied
   by 0.81123 and nothing else in the repository changes.  `target_mass_
   bound.py` is byte-identical (its shifts are ratios linear in g₁),
   `money_tagged_azz.py` is byte-identical (unpolarized electrons), and
   `closure_fom.py`'s ⁶Li A_∥ panel moves within its Monte-Carlo band —
   worst |spread/analytic − 1| 0.20 → 0.18 in the same x = 0.00105 bin,
   ⟨N⟩ 2970 → 2974, against a 95% band of 0.15.  The one unpublished
   number that scales is `money_polemc.py --ion 6Li`, whose δΔR rises by
   the full 1.233 (0.0496 → 0.0612 at x = 0.09); the published
   polarized-EMC reach is ⁷Li and does not move.  This is the **vector**
   (g₁) polarization only: the tensor sector's ⅓ dilution and the rank-2
   transfer 0.9219 are a different object and are untouched (plans/08 D9).
   *Engage (no longer blocking):* I. Cloët (ANL, local) — a VMC α–d
   overlap would replace the α–d D-state scenario, which is #15.
   *2026-08-28: the STRUCTURAL half of this was closed first (plans/08 D7).
   `ToyG1.g1_nucleus` now weights by Z and N exactly as `NuclearF2.f2a`
   does, every `Ion` slot is per-nucleon, and `beams.LI7` holds the verified
   VMC sums divided by Z and N so that the published ⁷Li path is
   bit-for-bit unchanged (`fastsim/tests/test_polarized_normalisation.py`).
   The ⁶Li 1/3 was being diluted a second time by the callers' 1/A; with
   that gone, per-nucleon g₁(⁶Li)/g₁(d) is 0.358 instead of 0.119, against
   the cluster picture's 0.29.  The same change gives ³He's proton term the
   ×2 the Bissey numbers intend and the struck triton of `polligen.tagged`
   its second neutron.  The VALUE — 1/3 against 0.81 — was untouched by it,
   but the gap it spans became 1.23 rather than 2.4 and the default sat at
   the optimistic end of it, which is what the decision above turned
   round.  Nothing published moved with the structural half:
   `fom.project_observables`' err_azz and err_g1_over_f1 are counting
   errors that carry no eff_pol, `phase_space_map.py` defaults to ⁷Li, and
   every published cos 2φ / Δ figure runs the transverse categories at
   θ_S = π/2 with unpolarized electrons, where the g₁ term enters only
   through cos θ_S = 0 — w_avg and a₂ are bit-for-bit unchanged there.  The
   tripling showed up only in the longitudinal vector-L term, i.e. in
   `closure_fom.py`'s A_∥ panel, whose estimator variance is
   1/(P_e P_z)²N to O(A_∥²) ≤ 2×10⁻⁵ — below its own Monte-Carlo band.
   The same term is where the 2026-08-29 value change shows up, and the
   current measurement of it is recorded in
   `fastsim/polli_fastsim/beams.py` beside the constant.*
   *2026-09-15, on the BAND rather than the value:* the sibling generator
   adopts the same central 0.811228 and the same per-nucleon convention,
   and disagrees only on the top edge — it forbids quoting an inclusive
   ⁶Li polarization without the band **0.81 … 0.91**
   (`LiPolGen/README.md` Conventions; `PHYSICS_CHANNELS.md:83` gives the
   cluster product's own span as 0.811 … 0.905, with Wiringa's ab-initio
   0.848 *inside* it rather than at its top).  The two trees have made
   different choices of what the band is *for*: 0.81–0.85 brackets the
   adopted cluster product against the ab-initio reading of the same
   quantity, which is the comparison this item made; 0.81–0.905 is the
   spread of the cluster product itself over the inputs of its two factors.
   Nothing here computes with the top edge — it is quoted, in
   `beams.py` and in Report 0's 2026-08-29 appendix row — so this is a
   labelling choice and not a live inconsistency, so it is recorded here as
   a clause on a closed item rather than reopened as one of its own.  *The
   call, for the authors:* keep 0.81–0.85 and say in the same breath that
   it is the cluster-product-against-ab-initio bracket, or widen to
   0.81–0.905 and say it is the cluster product's own spread.  Either is
   defensible; quoting the bare pair without saying which is not, because
   the two brackets answer different questions and the sibling tree already
   quotes the other one.  *Default until then:* 0.81–0.85 stands, as
   published.
7. **BeAGLE light-ion guidance / maintenance.** See A4.
   *Engage:* M. Baker, A. Jentsch, Z. Tu, W. Chang.
   *Fallback:* cluster-IA toy fragmenter (Phase-1 step 1.5.3).
   *Update 2026-06-12:* local build prepared up to the FLUKA wall
   (tools/beagle/); remaining user action = fluka.org registration.
   e+d/e+³He official EVGEN confirmed downloadable via xrootd for the
   control study.
8. **Nuclear (n)PDFs at A = 6,7.** EPPS21/nNNPDF coverage to confirm;
   polarized nuclear PDFs don't exist → effective-polarization convolution
   + CBT/TMT medium curves as scenarios.
   *2026-08-28: the CBT and TMT curves are no longer scenarios — both are
   digitized from the published figures (plans/02 step 1.2.2,
   `fastsim/polli_fastsim/data/SOURCES.md`).  Two limits remain and are
   this item: CBT computes ⁷Li at Q² = 5 GeV² while TMT computes nuclear
   matter at Q² = 10, so the comparison needs a target/scale transfer (a
   single valence strength factor).  *2026-08-29: the second half of this
   item is closed.  The unpolarized baseline both camps are transferred
   onto is EPPS21's ⁶Li F₂ per nucleon over CT18ANLO's free isoscalar
   nucleon — EPPS21's own proton baseline, so the fit cancels — not the
   hand-written 12-point table; its valence depletion is
   0.03105 over 0.35 < x < 0.65, giving s_CBT = 0.5322 and s_TMT = 0.2113
   where CBT's own model curve gave 1 and 0.397.  The transfer, and with
   it the whole reach, halves.  What the closure leaves is the baseline
   SPREAD — 0.01372 (nNNPDF3.0), 0.03105 (EPPS21), 0.05835 (CBT), with a
   90% CL Hessian band on EPPS21 alone of +0.039 / −0.041 — which is now
   the leading uncertainty on this figure of merit, wider than the
   statistics.*  The scale transfer itself remains open.  A
   second ⁷Li-specific polarized-EMC calculation, from either camp, would
   remove both — that is what to ask Cloët for.*
9. **b₁ and Δ theory for ⁶Li specifically.** Confirmed literature gaps:
   no b₁ prediction for any A > 2; no EIC Δ projection for any target.
   First-mover opportunity — co-author with theory.
   *Engage:* Cloët (ANL), Cosyn, Miller; lattice: Detmold/Shanahan.
   *2026-08-29:* written up as an ask, and packaged with the ⁷Li ones
   so that one message reaches the same people, in
   `docs/note_7li_theory_questions.md` §7(a).
10. **Radiative corrections on tensor observables** (A_zz, cos 2φ).
    Vector-case tools exist (DJANGOH/HERACLES); tensor RC uncharted.
    *Default:* **no band.**  An unpolarized QED study does not bound a
    tensor one, so the affected claims carry no correction and say so
    (plans/07 WP4, Report 2 §7).  What *is* measured is the unpolarized
    collinear-ISR migration, in `polligen/radiative.py`: +0.62 / +0.50 /
    +0.94 / +1.22% of Δ̂ at the four sweet spots in the published
    generator window, ≤ 2.9% once the low-Q² feed-in is opened
    (Phase-1 step 1.4, plans/08 D3).  *2026-08-29:* carried as a ⁷Li ask
    in `docs/note_7li_theory_questions.md` §6, where the spin-3/2 case
    has not even the spin-1 formal starting point of #14 to build on.
14. **Complete inclusive structure-function basis for spin-3/2** (⁷Li):
    rank-2 (b₁-analog) and rank-3 functions are not classified anywhere we
    can adopt; needed by the doubly polarized generator (plans/05 §5.2).
    *Engage:* Cloët, Cosyn — co-author opportunity.
    *Default:* rank ≤ 2 truncation, rank-2 shapes as scenarios.
    *2026-08-28:* the rank-2 *kernel* is already isotope-generic
    (Q_NN = (1,−1,−1,1) and c_eff = 3T for J = 3/2, verified), so what
    this item blocks is the structure function and not the machinery —
    scoped, with the rest of a ⁷Li channel, in plans/09 B3a.
    *2026-08-29:* the ask itself — the basis, its normalization against a
    named alignment tensor, the frame the alignment axis is defined in,
    and the map from the source literature's P_zz to `spin`'s T — is
    `docs/note_7li_theory_questions.md` §1.
15. **VMC two-cluster overlaps with m-dependence** (α+d S/D for ⁶Li,
    α+t P for ⁷Li) to replace the two-parameter radial forms whose tail
    dominates tagged acceptances (e+d control: BeAGLE tails 2–13× Hulthén).
    *Engage:* R.B. Wiringa (ANL, local). *Default:* β = 0.20–0.40 band.
    *2026-08-29, moments corrected 2026-09-15:* the cost of that band is
    measured on the ⁷Li P wave — ⟨k⟩ = 0.2362 / 0.3059 / 0.3717 GeV/c and
    P(k > 0.3 GeV/c) = 0.2209 / 0.3630 / 0.4803 at β = 0.20 / 0.30 / 0.40,
    a factor 2.17 in exactly the region where the tagged tensor
    asymmetries are O(1), while the angular moment ⟨P₂⟩ = −T/5 is fixed by
    Clebsch–Gordan and does not move — and the ask is
    `docs/note_7li_theory_questions.md` §3.  (The triples first entered here
    were `momentum_density(…, l_wave=0)`, the S-wave form, evaluated for a
    channel whose wave is `l_wave = 1`; they are retired in
    `tools/retired_numbers.json`.  On the P wave the band is wider in ⟨k⟩
    and far heavier in the tagged region — a fifth of the ⁷Li α spectrum
    is already above 0.3 GeV/c at β = 0.20 — while the direction of the
    argument, that this tail is what a VMC overlap must settle, is
    unchanged.)
    *2026-09-15:* the ask acquired a second, sharper part — the **sign**
    of the α–d D radial is now a live physics input of the tagged ⁶Li
    channel, not a convention.  `polligen/tagged.py` applies the i^L phase
    of the momentum-space partial-wave expansion (plans/00 run 19), so the
    sign of the S–D interference — and with it the sign of every tagged
    A_zz — follows the sign of ψ₂/ψ₀ directly.  The model takes it
    deuteron-like (ψ₂/ψ₀ > 0, as AV18's w(k)/u(k) is at low k).  The
    sibling generator's VMC α+d overlap supports that over 0.134–0.444 GeV/c
    and nowhere else: it measures sign(ψ₂/ψ₀) = −1 below the α–d S node at
    0.134 GeV/c, +1 between the nodes, and −1 again above the α–d D node at
    0.444 GeV/c, and the node-free Hulthén forms used here can represent
    neither reversal.  Nor is the accepted sample confined to the supported
    window: 27 / 21 / 26 % of the Yellow-Report-accepted α lie above
    0.444 GeV/c at the three configurations — where LiPolGen's own
    acceptance-weighted VMC A_zz^tag changes sign between k = 0.40 and 0.50
    (−0.08 → +0.12) against a Hulthén −0.73 / −0.59 — and 41 / 28 / 37 %
    of the tagging-optics sample lies below 0.134 GeV/c.  The k = 0.325 GeV/c
    headline bin of money plot 4 is inside the supported window; the tails
    on either side of it are not.  A VMC overlap would settle
    both the tail and the sign in one object; the statement of the adopted
    sign lives next to `P_D_LI6` in `tagged.py`.
    *2026-09-15, the overlaps are in hand — this item is delivered, and
    the default is superseded by measurement.*  The ANL VMC α + d and
    α + t two-cluster overlaps this item asks for exist and are shipped in
    the sibling generator, `LiPolGen/data/vmc` (ANL's own tables — the
    2004 `overlap_old/` and the 2024 `momenta/` sets, retrieved through the
    Wayback mirror; the `README` there carries the URLs), so they can be
    obtained here independently of that tree.  Measured against them, the
    β = 0.20–0.40 band this item proceeds on **does not bracket the VMC
    α–d density at either end**: in the Roman-Pot window VMC's
    P(k > 0.20 GeV/c) = 0.2458 sits *above* the whole band
    (0.0836 / 0.1508 / 0.2005 at β = 0.20 / 0.30 / 0.40), while in the far
    tail its P(k > 0.45 GeV/c) = 0.0020 sits *below* it
    (0.0046 / 0.0156 / 0.0321) — the Hulthén form is too soft where the
    pots look and too hard beyond, and P_D itself is 0.0193–0.0207 against
    the 0.0867 scenario `tagged.P_D_LI6` adopts to reproduce the 0.87
    vector dilution.  The consequence for the tag is not a single
    direction: on the Yellow-Report high-acceptance optics the ⁶Li α-tag
    fraction moves **0.0264 → 0.0348** (×1.32) at 10 × 99.5, while on the
    tagging optics it moves **0.2551 → 0.2486** (×0.975) — opposite signs,
    which is exactly what the one-sided band plans/05 asks for would get
    wrong.  Those two ratios, ×1.32 and ×0.975, are what
    transfers here: the pair Report 0 §5.4 and Table 3 publish at the same
    configuration is 2.4% and 25% — 2.5% before this run's S–D phase
    correction, which is the value LiPolGen's 0.0264 corresponds to — so on
    the VMC density it reads 3.2% and 24.4%; whether it is re-quoted as a
    band is #29.  *Superseded 2026-09-16: the ratios above are the
    sibling's, and the α + d tables are now read by this generator
    (`tagged.li6_alpha_channel(wave='vmc')`), so the pair is measured
    rather than transferred — ×1.394 and ×0.974, i.e. 3.4% and 24.8%.*
    The band
    clause of plans/05 §169–182 reads "until VMC lands"; it has landed,
    and what remains of this item is the adoption, not the ask — with the
    **sign** of ψ₂/ψ₀ above the live half of it.
16. **FSI for cluster spectators beyond IA** (α rescattering off DIS
    debris; deuteron case solved in PRC 97:035209, pole extrapolation).
    *Engage:* Cosyn, Sargsian. *Default:* IA, quoted at small |t′|.
    *2026-08-29:* the ⁷Li form of the same ask — a triton spectator, and
    the two-body currents an A = 7 cluster pair carries — is
    `docs/note_7li_theory_questions.md` §5.
17. **HepMC3 convention for ion spin states** — ☑ **the schema is written
    and delivered 2026-09-15**; what is left is the proposal, not the
    definition.  None exists upstream in HepMC3 or in the EIC stack;
    plans/05 step 5.D defines named attributes and proposes them upstream.
    *Engage:* ePIC MC/software group. *Default:* our attribute schema.
    *2026-09-15:* the sibling generator has written that schema down and
    implemented it — `LiPolGen/docs/HEPMC3_CONVENTION.md`, against
    `lipolgen::HepMC3Writer`, and the document names this item as the thing
    it answers.  The schema is a block of **`GenEvent`-level attributes
    (id 0)** carrying the spin state beside the kinematics —
    `spin_J`, `spin_M`, `struck_cluster_m` (NaN when inclusive), `lam_e`,
    `P_e`, `P_z`, `P_zz`, `spin_axis_theta`, `spin_axis_phi`,
    `spin_category`, `run`, `bunch`, plus `channel` and the DIS /
    spectator / coherent kinematics (`dis_x`, `dis_Q2`, `dis_y`, `dis_phi`,
    `spectator_k`, `spectator_cos_theta`, `spectator_phi`, `alpha_s`,
    `pt_s`, `t`, `x_pom`) — with a per-particle `pol` attribute written
    only when it is not "unknown", weight slot 0 named `nominal` and one
    `spin_weight_<k>` per spin category, and the ion written as a 10-digit
    nuclear PDG code on a status-4 beam particle in Asciiv3 (GeV/mm).  Two
    of its constraints are this repository's own lessons, sourced to
    `tools/fullsim/ion_gun_hepmc.py`: HepMC2 `IO_GenEvent` handed to
    DD4hep's `HEPMC3FileReader` fails as a bare EOF rather than a parse
    error, so Asciiv3 is a requirement and not a style choice, and a vertex
    with zero incoming particles does not survive `ReaderAscii`.  The item
    stays on the tracking board until the ePIC MC group has been asked.
18. **Coherent diffractive model for (tensor-polarized) ⁶Li** — no
    light-nucleus coherent-fraction prediction exists (lightest
    published is Ca); the tensor cos 2φ of the coherent yield has no
    calculation for any A > 2 (all forward citations of the deuteron
    template checked, 2026-08-10). Template: Mäntysaari et al.
    arXiv:2408.13213 / PLB 858:139053 — now digitized into
    `coherent.MANTYSAARI_A2_DEUTERON` and scaled to ⁶Li (plans/06
    §6.4b: deformation term ε_B0 ∈ −(0.04–0.13), sign flip predicted).
    *Engage:* the Mäntysaari–Schenke group — the concrete ask is
    rerunning their IP-Glasma polarized-deuteron setup with an α–d
    cluster density; the ⁶Li case is a clean null test
    (Q(⁶Li) = −0.0806 fm²). *Default:*
    `polligen/coherent.CoherentScenario` bands (f₀ = 0.04 ×2÷2,
    B = 50 ± 10 GeV⁻², deformation + flat-gluonic modulation) —
    plans/06.  *2026-08-28:* the ⁷Li version of the ask (an α + t
    density, and the rank-2 slope amplitude beyond linear order) is the
    same conversation with the same group: ε_B0 cannot be rescaled to
    ⁷Li — the linear form gives −2.2 to −4.5, |ΔB₀| > B and c₂ > 1
    inside the tagged window — so ⁷Li has no scenario until the
    amplitude exists.  Scoped in plans/09 B3a, and written out as one
    ask with two densities — α + t and α + d — in
    `docs/note_7li_theory_questions.md` §2 and §7(b).

## Still open — detector / software

11. **Far-forward transfer matrices & optics at Li rigidities** (RP/OMD
    reconstruction in EICrecon is tuned per beam setting).
    *Engage:* ePIC FF WG (A. Jentsch). Phase-2 step 2.2.
    *2026-08-28:* the item has split into its two halves, each with its
    own owner and its own written-down ask.  The **optics** half was
    promoted to a plan of its own — plans/10 **D1–D3**: σ_θ(h, v) and
    Δp/p at the IP for a ⁶Li/⁷Li fill (D1), which cooling scenario is the
    baseline for ion running (D2), and β* for light-ion running (D3), for
    C-AD with the ePIC FF WG.  plans/10 §10.3 answers D1 provisionally by
    scaling the Yellow Report's own proton tables with an equal-emittance
    assumption calibrated on gold, and the sharp question inside it is
    whether a light-ion *tagging* optics can exist at all, since nothing
    else recovers the coherent channel.  The **transfer-matrix** half
    stays with the ePIC FF WG and is written down as plans/09 **D3**:
    R₁₂, R₃₄ and the pot dispersion D were measured at all three
    configurations on 2026-08-28 (`farforward.POT_LEVERS`: 19.24 / 21.25 /
    29.97 m, 4.56 / 3.35 / 2.93 m and 0.311 / 0.287 / 0.292 m, the 5 × 41 vertical lever read off a zero-insertion scratch geometry on 2026-08-29), so plans/09
    quotes millimetres everywhere; R₁₁, R₂₁, R₂₂ and D′ remain unmeasured
    — *superseded 2026-09-16: the angle row is measured in the same
    zero-insertion geometry, R₁₁ = 1.148 / 1.227 / 1.852, R₂₁ = −0.0837 /
    −0.0651 / −0.0209 rad m⁻¹, R₂₂ = −0.4955 / −0.3060 / +0.1944 and
    D′ = 0.0175 / 0.0179 / 0.0182 rad (`farforward.POT_SECOND_ROW`,
    plans/09 B1 and D3), so what stays with the FF WG is the lattice
    question alone* —
    the 10σ offsets in `beamline_*.xml` are still marked a *"rough
    extrapolation"* at 5 × 41, and what the FF WG is now asked is the
    lattice question: `beamline_5x41.xml` (41 GeV proton) and
    `beamline_5x41_He4.xml` (Z/A = 0.5, 82 GV) give R₁₂ = 19.24 and 29.81 m
    at the same ring setting, a factor 1.55 in every millimetre and 0.64 in
    the 5 × 41 aperture.
12. **Geant4/DD4hep light-ion & excited-ion primaries** (10LZZZAAAI codes
    from BeAGLE; DD4hep had fixes ~PR #920; `sanitize_hepmc3.py` exists).
    Verify in Phase-2 step 2.1.4.
13. **Afterburner + beamline configs for Li** — verified absent; three
    concrete artifacts to add (EicConfigurator.cc preset; beamline_*.xml;
    BeAGLE runcard). ⁶Li can start from He-4/deuteron files (same Z/A).
19. **Roman-Pot charge discrimination for A/Z = 2 fragments** — an
    intact ⁶Li, an α, and a d from breakup have identical rigidity AND
    velocity; only dE/dx (∝ Z²: 9/4/1) separates them. The EICROC
    AC-LGAD chain records pulse amplitude (for charge-sharing) but no
    EIC document addresses Z-ID in the RPs; the documented concept is a
    Z² Cherenkov behind the IR-8 secondary-focus pots (arXiv:2211.15746
    §VIII, 2602.04636). Make-or-break for the coherent intact-⁶Li tag
    (plans/06 §6.2). *Engage:* ePIC FF WG (A. Jentsch). *Default:*
    assume no event-by-event Z-ID at IP6 → two-component |t| fit.
    *2026-08-26 (plans/09, `reports/nanowire_far_forward`):* **the
    question has been asked of the wrong technology.** A nanowire
    supplies a candidate mechanism (below), but the INCUMBENT already
    carries more information: EICROC provides per channel an 8-bit
    40 MHz SAR ADC for charge (the ToT of ALTIROC was replaced by it for
    dynamic range) behind an AC-LGAD with a 30 µm active thickness, over
    the **four** planes ePIC already has (2 stations × 2 layers).
    *Retracted 2026-08-27 (plans/09 §9.2):* this entry priced that
    advantage as "**4.8σ per plane** against a nanowire's one bit", and a
    σ is the wrong figure of merit — a gap over the quadrature sum of two
    Landau core widths is neither a separation power nor a fake rate, and
    what puts an α inside a ⁶Li's window is the Landau *upper tail*.
    Restated as an α fake rate at a matched 95% ⁶Li efficiency over the
    four planes (`nearbeam_zid_power.py`, a sampled Landau, 1.5×10⁶
    events): **2.3×10⁻⁵** for the 8-bit per-plane likelihood ratio, the
    optimum, against **3.1×10⁻⁵** for one bit per plane with a
    majority-of-k — **a factor 1.4, not orders of magnitude**, because
    the power comes from requiring coincidence across planes rather than
    precision within one, and the two species are far apart (MPV 31.7
    against 75.2 keV). More bits is not automatically better Z-ID: a
    truncated mean, the standard analogue dE/dx estimator, gives
    2.7×10⁻³ and a plain sum of the four planes 5.3×10⁻². Where the
    nanowire actually loses is **geometric fill factor** — the
    coincidence needs every plane to record the track, which silicon does
    ~99% of the time and a wire comb only over its fill (25–50% in the
    published devices), capping the reachable ⁶Li efficiency at 0.68–0.94
    over four planes, so 95% is out of reach at any working point. That
    is a fabrication number rather than an information-theoretic one, and
    it is the actionable thing to put to the MEP group.
    **Ask the incumbent first** (plans/09 D1): a Geant4 study through the
    four layers, plus EICROC's input charge dynamic range in fC and the
    sensor's gain-suppression curve at ~9 MIP. One person-month, no
    hardware, closes #19 either way.
    *Better still, and free:* the background #19 exists to reject is
    ⁶Li → α + d, and that is **two hits**. The relative momentum
    (κ = 60.7 MeV/c) is transverse and unboosted, so the α (4p_u) and the
    d (2p_u) take opposite kicks and land a median 10.9 / 10.7 / 17.3 mm
    apart at 18×275 / 10×100 / 5×41 — **21 to 35 pixels of the existing
    500 µm pitch**, and 6 to 82 over the 16–84% band (3.2–41.2 mm across
    the three configurations). An intact ⁶Li is
    one hit. And the second fragment is a *veto*: conditioned on an α
    that fakes a coherent tag, the partner deuteron is recorded in **84%**
    of events at the tagging optics (0.02–0.25 at the published ones,
    where the fake rate is 10⁻³–10⁻⁴ anyway), so topology beats dE/dx
    here in sensors that already exist — subject to how far the pot
    stations extend, which is B1. Measured 2026-08-28 with both fragments
    sampled from one relative momentum (plans/09 B4 §9.2); the
    6.7 / 18.4 / 44.8 mm this entry carried until then were a single
    k = 40 MeV/c at the retired rigidity-scaled energies, with no
    dispersion.
    *The nanowire mechanism, for the record:* A
    superconducting nanowire latches — its pulse amplitude is the
    diverted bias current and is identical for hadrons, muons, pions and
    showering electrons (arXiv:2510.11725, 2410.00251) — so any
    pulse-height scheme is dead. What the deposit sets is the FIRING
    THRESHOLD, I_th/I_c = 1 − 2 r_s/w with r_s = √(Q/(eπcρ(T_c−T_0)))
    (Argonne's own Eqs. 1–2, arXiv:2312.13405). Since dE/dx ∝ z² at
    fixed β and ⁶Li at 137.5 GeV/u has βγ = 148 against 128 for their
    calibration proton, **r_s ∝ z**: 134 / 268 / 402 nm for d,p / α /
    ⁶Li, anchored on their *measured* 134 nm. At w = 1 µm — the
    microwire width that already exists — the three turn-ons are 0.73 /
    0.46 / 0.20 I_c, and **two bias points at 0.33 and 0.60 I_c tag Z by
    the firing pattern alone**, both below the 0.80 I_c dark-count wall.
    The scheme is a granted patent (US 8,872,109) demonstrated on singly-
    vs doubly-charged lysozyme. One plane is a one-bit measurement (α →
    ⁶Li confusion 20–25%), so 3–5 planes are needed for sub-percent.
    *Blocking measurement, and only if the LGAD route fails* (plans/09
    D5): the α turn-on curve on the same wires as the 120 GeV proton —
    their α analysis is "underway". Note it tests the √Q law and **not**
    z² at fixed β: their ²⁴¹Am α differs from the 120 GeV proton almost
    entirely through 1/β² (β = 0.054 vs ≈1). Nobody has ever varied Z at
    fixed velocity on one of these devices, so the "interpolation, not
    extrapolation" argument is only half true.
20. **Roman-Pot cutout geometry and Li beam divergence** (reconstruction-
    chain note, 2026-08-24). The near-beam cut is angular,
    pT_cut = 10σ_θ·A·p_u, so the same optics gives tag acceptances of
    67% / 9% / 10⁻⁸ at the then-assumed 20.5 / 50 / 137.5 GeV/u for ⁶Li
    (superseded: the energies are γ-matched, plans/10) (constant 0.20 GeV
    in the code → 13.5%); and the pots' rectangular cutout has its sides
    parallel to the vertical spin axis, faking ⟨cos 2φ_t⟩ ≈ 0.5 for a 25%
    aspect ratio (physics a₂ ≈ 0.036) unless the spin-state ratio is
    used. *Needed:* σ_θx, σ_θy (β*, emittance) at the light-ion energies
    and the pot geometry (`reco.rp_measure` takes both). *Engage:* ePIC
    FF WG / C-AD optics. *Default (2026-08-28, plans/10):* the
    **per-configuration** Yellow Report divergences, anisotropic, from
    `farforward.sigma_theta_for` — 10σ_θ = 2.20 × 3.80, 1.80 × 1.80 and
    0.92 × 0.92 mrad for ⁶Li at 5 × 41 / 10 × 100 / 18 × 275 high
    acceptance — with the measured pot aperture as a second constraint
    per axis, and the tagging optics of Report 1 §6.1 (0.36 × 3.80,
    0.19 × 1.80, 0.12 × 0.92 mrad at 1/6.8, 1/12.8, 1/9.5 of the
    luminosity) as the setting at which the coherent channel is
    published.  The single proton-derived pair the code carried until
    then — 73 μrad high acceptance and 164 μrad high divergence at every
    configuration, from a 0.20 / 0.45 GeV p_T at 275 GeV — is **retired**
    (`--optics legacy` keeps it as the dated record; the "149 μrad" this
    entry used to give was `reco.SIGMA_THETA_HD`, an unreconciled second
    copy of the same constant, unified on 164 by plans/10 A1b).  Angular envelope
    10σ_θ·A·p_u throughout; ratio estimator.
    *2026-08-25 (refs/):* the ePIC pots are sensor planes around a
    horizontal SLOT (Jentsch DIS 2023, slide 15) — wide in x for the
    beam's momentum spread and dispersion, tight in y — so the cutout
    aspect ratio is < 1 and the fake ⟨cos 2β⟩ w.r.t. the vertical axis
    is large and positive; HERA's proton beam had a 45 (x) vs 100 (y)
    MeV transverse-momentum spread at the IP (ZEUS NPB 816:1), i.e.
    factor-2 anisotropies are the norm; ePIC's full-beam-effects
    simulation gives ΔpT ≈ 40 MeV at 275 GeV with the detector alone at
    ≤ 1.5% ("beam effects the dominant source", slide 20) —
    `reco.rp_measure(cut_scale_xy=(2.5, 1))` is the new default of
    money plot 6R. Still needed: the slot dimensions and the light-ion
    optics.
    *2026-08-26 (measured, `tools/fullsim`):* an intact ⁶Li shot through
    the ePIC geometry (`ion_gun_hepmc.py` → npsim, epic-main of
    jug_xl-nightly, 84 points in p_T × azimuth) **inverts the aspect
    ratio**. The pot silicon does surround a horizontal slot, as the
    2026-08-25 entry says; but the far-forward optics image an IP angle
    onto the pot plane with R₁₂ = 19.24 / 21.25 / 29.97 m horizontally
    against R₃₄ = 4.56 / 3.35 / 2.93 m vertically (re-measured 2026-08-28; the 5 × 41 entry on a zero-insertion scratch geometry, 2026-08-29), a
    factor 4.2 at 5 × 41, 6.3 at 10 × 100 and 10.2 at 18 × 275, so what clears the slot is
    the HORIZONTAL angle: the boundary
    is |θ_x| ≳ 2.50 / 1.51 / 0.53 mrad in the 5×41 / 10×100 /
    18×275 optics (p_T = A p_u |θ_x| = 0.61 / 0.90 / 0.44 GeV for the ⁶Li
    at the γ-matched 40.8 / 99.5 / 137.5 GeV/u; the 0.25 / 0.41 / 0.85
    this entry carried until 2026-08-28 priced the September-2024 edges at
    the retired rigidity-scaled momenta), against |θ_y| ≳ 0.92–2.12 mrad
    where the vertical plane is open at all, and nothing at 5 × 41.
    In `rp_measure` terms
    that is `cut_scale_xy ≈ (1, 1.4–1.7)` where both axes are open, not
    (2.5, 1) — a factor 3.5–4.4 the
    wrong way — and the fake ⟨cos 2β⟩ about the vertical spin axis is
    therefore large and **negative**, not positive. The tagged fraction
    falls with it (B = 50 GeV⁻², `evgen/scripts/nearbeam_aperture_scan.py`,
    re-run 2026-08-28 at the γ-matched momenta and on the re-measured
    aperture): **9.4×10⁻¹⁰ / 2.0×10⁻¹⁹ / 1.2×10⁻⁵** through the measured
    aperture at 40.8 / 99.5 / 137.5
    GeV/u, against **7.2×10⁻⁸ / 6.2×10⁻²⁷ / 7.1×10⁻¹⁴** through the
    Yellow Report high-acceptance envelope on the scan's own convention
    (the envelope horizontally, the larger of silicon and envelope
    vertically; the pure 10σ envelope in both planes gives 7.2×10⁻⁸ /
    1.2×10⁻²⁶ / 7.8×10⁻¹⁴ — the two are different quantities sharing a
    name, `tools/fullsim/README.md`).  **Which of the two binds changed twice
    on 2026-08-28.** Priced against the retired single 73 μrad, plans/08
    §8.4 read 2.8× / 1.9× / 1.4× and concluded the envelope is never
    binding; the per-configuration envelopes turned that into
    0.91× / 0.75× / 1.12× on the September-2024 aperture; and the
    re-measurement in the current geometry makes it **1.14× / 0.84× /
    0.58×**, so the SILICON binds at 5 × 41 and the machine at the other
    two, by eight orders of magnitude at the top (1.2×10⁻⁵ at the silicon
    against 7.1×10⁻¹⁴ at the envelope).  The 1.4×10⁻² and 5×10⁻⁵
    this entry carried until then were the September-2024 aperture at the
    retired rigidity-scaled 20.5 and 50 GeV/u.
    *2026-08-26 (the chain, run on it — `money_cos2phi_coherent_reco.py
    --rp-aperture measured`):* at the LOW configuration the measurement
    survives. Acceptance 37.7% → 1.42%, N_tag 8.3×10⁶ → 3.1×10⁵, the
    acceptance-induced ⟨cos 2β⟩ **+0.426 → −0.772**, two of the four |t|
    bins instead of four (the cutout leaves |t| = 0.061–0.273 GeV², and
    the higher bins no longer separate the seven harmonic columns), and
    δa_t 0.0014 → 0.0482 and 0.0026 → 0.0150 in the two survivors — a
    factor 6–34. a_e is still recovered: 0.0073 ± 0.0045 and
    0.0091 ± 0.0045 against an injected 0.0100. At MID and TOP the
    aperture leaves no accepted recoil in the binned window at all.
    (This chain entry is itself at the retired rigidity-scaled menu —
    "LOW" is 20.5 GeV/u, superseded by plans/10 — and at the
    September-2024 aperture: 37.7% → 1.42% and +0.426 → −0.772 were
    computed against a 2.0 mrad edge that measures 2.50 mrad in the
    current geometry, so the whole chain entry must be re-run before it is
    quoted again, plans/09 B1.)
    **The conclusion drawn here — "so the coherent programme is a
    low-energy programme … for a second and stronger reason" — is
    withdrawn (2026-08-28).** It rested on the measured aperture binding
    everywhere, which it does not, and on the rigidity-scaled energies.
    Against the per-configuration envelopes the aperture is the binding
    constraint only at the top; and at the tagging optics, where the
    channel is now measured, the ordering reverses outright — 0.37 / 0.25
    / 0.33 tagged at 1/6.8 / 1/12.8 / 1/9.5 of the luminosity, with the
    TOP configuration, which has four times the coherent rate at equal
    luminosity, the best covered (Report 4 §3).  What survives is the
    weaker statement the angular envelope alone makes at the *published*
    optics, where no configuration is usable.
    *Caveats:* one event per scan point, 30° azimuthal steps,
    no beam envelope, and a September-2024 `epic-main`. *Action:* repeat
    on the current release with beam effects, and put the aspect ratio to
    the FF WG as a question with a number attached.
    *2026-08-26 (plans/09 §9.4 — **the geometry has moved**):* reading
    the current `main` of `eic/epic` directly, the pot layout changed
    after the September-2024 snapshot `tools/fullsim` ran in. Modules
    went 32 × 32 mm → **16 × 16 mm**; the single energy-independent
    insertion became **per-energy 10σ offsets** in `beamline_*.xml`
    ("These are the ten-sigma cuts for the Roman pots, translated to the
    physical layout we currently have. They are not perfectly ten-sigma
    for reasons of physical geometry."); and the 1 mm aluminium RF
    shields are **commented out** ("we don't know if we will even need it
    … Oct. 2025"). The old 32 mm block gives 32/30.6 m = 1.046 mrad
    against the 1.03 then measured — agreement to 1.5%, so measurement and
    file reading confirm each other. The current 16 mm block **measures
    0.53 mrad** against the 0.52 predicted here, *below* the 0.9169 mrad
    Yellow Report high-acceptance envelope at 18 × 275 (plans/10) but 4×
    the 0.12 mrad tagging-optics envelope; at 5 × 41 the per-energy
    insertion moves the other way (29.6 mm inner edge) and the aperture
    **measures 2.50 mrad**, outside the 2.20 mrad envelope.
    **The re-measurement is done** (2026-08-28, plans/09 B1,
    `tools/fullsim/README.md`): #20 is now measured rather than assumed,
    and every aperture-conditional number in this file has been moved onto
    it except the reconstructed-chain entry below, which is flagged. What a closer approach
    is WORTH is unaffected and now curved rather than tabulated:
    `nearbeam_aperture_scan.py` prices every aperture per configuration.
    On the Yellow Report divergences (2026-08-28, plans/10 A4) the
    machine envelope binds at 10 × 100 and 18 × 275 and the re-measured
    silicon at 5 × 41 (2.50 against 2.20 mrad), so a closer approach alone
    buys a factor 77 there and nothing at the other two; under the tagging
    optics of Report 1 §6.1 a layer
    that follows the 0.36 / 0.19 / 0.12 mrad envelope is the difference
    between no tag and 0.36 / 0.25 / 0.32 with seven populated |t| bins each
    (`nearbeam_reach_gain.py`; the earlier ×26 / ×569 were artefacts of the
    73 μrad divergence).
21. **Hadronic-method y resolution at y = 0.01–0.05 for e + light ions**
    (same note). Three of the four inclusive sweet spots sit at
    y = 0.010–0.025 where the electron alone gives δy/y = 50–120%; the
    mixed (eΣ) method needs the hadronic final state, which polligen does
    not generate. *Needed:* δy_Σ/y vs y from the ePIC inclusive WG (or a
    BeAGLE/PYTHIA e+Li sample through eic-smear). *Default:* 15–30%
    band (`reco.hadronic_y`); purity 0.75–0.83 per super-bin at 15–20%
    (0.64–0.68 at the 25% default).
    *2026-08-25 (WP3-HFS, plans/07):* the chain now takes the hadronic y from
    a hadronic final state through a hadron-side detector response
    (`polligen/hfs.py`; PYTHIA 8 sample via `tools/pythia8`, toy stand-in
    locally). Toy result: Σ-method δy/y = 0.28 / 0.17 / 0.24 / 0.07 at the
    sweet spots with a 50 MeV calorimeter noise floor (9–12% without noise),
    i.e. the 25% default is the noise floor acting on Σ_h ≈ 0.2 GeV. *Needed
    now:* the PYTHIA sample (one eic-shell command) and the ePIC calorimeter
    noise/threshold floor at Σ_h ≈ 0.2–0.5 GeV.
    *2026-08-26 (the PYTHIA sample exists):* 8 M events over the three beam
    configurations, generated natively (`tools/pythia8`). Σ-method
    δy/y = **0.55 / 0.28 / 0.50 / 0.15** at the sweet spots with the 50 MeV
    floor — the toy was optimistic by 0.04–0.05 absolute at every one,
    because it put 0.03 of Σ into neutral hadrons where PYTHIA puts 0.11 —
    *(superseded 2026-08-27: at the corrected spots HCal objects — neutral
    hadrons plus untracked charged particles — carry 0.09–0.10 of Σ within
    acceptance; see the acceptance entry below)* —
    and 0.28 / 0.21 / 0.24 / 0.11 at LOW against 0.74 / 0.34 / 0.69 / 0.18
    at TOP. **Half of this question is therefore answered**: what remains
    is only the ePIC noise/threshold floor, which is what the scan
    0 → 25 → 50 → 100 MeV turns into 0.20 → 0.32 → 0.54 → 1.01 at y = 0.005.
    Reco purity at the sweet spots falls from 0.64–0.68 to 0.40–0.73 with
    the real final state.
    *2026-08-27 (acceptance, `hfs_acceptance.py`, Report 2 §3 Figure 2):*
    80 / 87 / 83 / 92% of Σ_h is captured at the four mid sweet spots;
    17 / 8 / 15 / 7% escapes forward beyond |η| = 3.7 (the target-
    fragmentation side of a W ≈ 6–10 GeV system), 1–6% is below threshold.
    Through the full response 70 / 74 / 74 / 85% is captured: a 13–28% scale
    bias on y_Σ that the library *reproduces* in the pseudo-events and the
    bin-centering factor absorbs; `--hfs-calibrate` (per-cell mean captured
    fraction, the analysis's own calibration) takes the 5R purities to
    0.52–0.76 at unchanged errors, and a residual 1% scale error moves Δ̂ by
    0.2–0.7%.  Not the resolution driver: the ePIC nominal reach of 4.0
    recovers a third of the escape with δy/y unchanged.  What the ePIC
    inclusive WG could settle: the forward calorimeter reach and thresholds
    at Σ_h ≈ 0.2–0.5 GeV, alongside the noise floor.
    *2026-08-28 (review):* the calibration is now keyed on the reconstructed
    (x_mixed, Q²_e) rather than the true cell (purities 0.56 / 0.59 / 0.64 /
    0.75), the acceptance is applied in the detector frame (+0.01 on the
    captured fractions: 0.80 / 0.87 / 0.83 / 0.92), and the PYTHIA library
    was regenerated without the m̂ ≥ 4 GeV floor that had removed x < 16/s.
    The 50 MeV noise is labelled as this programme's stand-in; the ePIC
    floor at Σ_h ≈ 0.2–0.5 GeV remains the one number to obtain.
    *2026-08-25 (web search, refs/README.md):* now bracketed by documents —
    ATHENA proposal JINST 17:P10019 Sec. 3.1/Fig. 22 (e−Σ or DA for
    y ≲ 0.1; ≈ 25% y resolution at y ≈ 0.01 → ≈ 10% at y ≈ 0.1; JB 20–30%);
    the ePIC kinematic-fit study (S. Maple, Dec 2024 seminar) smears
    σ(δ_h)/δ_h = 25% and shows Σ/JB/DA widths of 0.2–0.3 in
    0.01 < y < 0.05 with the electron method flat; Arratia et al. NIM A
    1025:166164 Fig. 5 (ATHENA fast sim, Q² > 200): RMS(y)/y ≈ 0.13–0.17
    for IΣ/DA/JB at y ≈ 0.05–0.2. The 15–30% band stands (25% default = ePIC's
    own value). Remaining gap: no published ePIC full-simulation
    number at Q² ≈ 1–3 GeV², y ≈ 0.01, for e + light ions specifically.
    *2026-09-15 (plans/07 WP3 addendum, Report 2 §5.1):* the question is
    bounded from the other side as well. Run through the full chain at a
    reconstructed y ≥ 0.05, where the bin-centre δy/y of the electron method
    alone falls to 0.08–0.22 from the 0.46–1.18 of the published spots, the
    electron method is the better of the two and needs no hadronic final
    state: purity 0.73 / 0.80 / 0.76 / 0.77 against 0.65 / 0.64 / 0.70 / 0.69
    for the mixed method in the same bins, D = 0.95–1.01 against 0.90–0.99 and
    δÂ 0.98–0.99 of it. The missing number stays the one named above — it is
    needed only below y ≈ 0.05, which is where every published sweet spot
    sits, because above it the hadronic sum can be dropped.

## Still open — sourcing (raised by the consistency review of Reports 0–4)

Three passages the review of 2026-09-02
(`docs/consistency_review_2026-09-02.md` §2.5) could not settle from the
repository: each states a number whose producer or citation exists
nowhere in it, so only the authors can say what it is. The review left a
placeholder edit for each, unapplied. Two of them were settled on
2026-09-06 by locating the sources (items 22 and 23); item 24 remains an
author decision.

22. **The source of ξ_p ≤ 0.015** (review finding F226) — ☑ **closed
    2026-09-06.** The bound is the EIC Conceptual Design Report's hadron
    beam-beam design limit and not this programme's own: §1.2 p. 5 lists it
    among the design parameters ("hadrons: ξ_p ≤ 0.015; electrons:
    ξ_e ≤ 0.1"), §1.4 p. 11 calls the value "suggested by previous
    performance and confirmed by extensive beam-beam simulations", §3.1.1
    p. 95 adopts RHIC's achieved ξ_p = 0.015 in proton-proton collisions as
    the EIC design value, and §4.6.1 pp. 392 and 395 carry it into the
    interaction-region design; Tables 3.3 and 3.4 run the proton at 3/3,
    12/12, 12/12, 14/14 and 15/9 × 10⁻³ (h/v) at 18 × 275, 10 × 275,
    10 × 100, 5 × 100 and 5 × 41 GeV, at or under it. For ions the CDR does
    not restate the number and does not leave it open either: §1.4 p. 14
    scales the hadron tune shift with Z and inversely with A, Table 3.5
    runs the Au beam at 1–5 × 10⁻³, and Appendix A §A.2.2 caps the two
    interaction regions together at 0.03, i.e. 0.015 apiece — so applying
    the proton limit to a lithium beam (Z/A = ½) is conservative. Cited as
    Report 3 [11] and Report 4 [15] on 2026-09-06; the Yellow Report
    carries no beam-beam row (confirmed).
    *The question was:* "Where does ξ_p ≤ 0.015 come from: a published EIC
    beam-beam limit for the hadron beam (which document and table?), or is
    it this programme's own assumed bound, to be labelled 'ours to
    specify' in Table 9 as rows 11, 12 and 16 are?"
    *Engage:* none needed.
23. **The upper end of the ePIC momentum resolution at 1 GeV/c**
    (review finding F237) — ☑ **closed 2026-09-06.** The η-resolved ePIC
    full simulation is in the document Report 3 [7] already named: the
    Preliminary Design Report of September 2024, Figure 8.9 ("ePIC 24.08",
    single pions, five η ranges) — public snapshot
    doi:10.5281/zenodo.13866213, the same figure standing as Fig. 3.56 of
    the Preliminary Technical Design Report v3.1,
    doi:10.5281/zenodo.18271601. Digitised 2026-09-06 by two independent
    agents, it reads 3.00 / 1.08 / 0.38 / 1.13 / 2.55% at p = 1 GeV/c in
    −3.5 < η < −2.5, −2.5 < η < −1, −1 < η < 1, 1 < η < 2.5 and
    2.5 < η < 3.5. No slice gives 0.6%: the 0.6% end came from the 1.4 T
    BaBar-magnet reference design of arXiv:2305.15593 (X. Li, 2023), not
    from ePIC. Report 3 Table 8 and the `reco.py` docstring are restated on
    the figure — the barrel model's 0.50% is ≈ 30% pessimistic against
    0.38%, while the two backward slices, which hold eleven of the twelve
    sweet spots, are right to ≈ 8% — and Maple's 0.45% in 0 ≤ η ≤ 0.5 is
    Report 3 [10].
    *The question was:* which ePIC document gives the upper end of the
    resolution band Table 8 quoted, Maple's seminar slide 47 being the only
    ePIC full-simulation momentum resolution then held in `refs/` and
    giving 0.45% at 1 GeV/c in the single slice 0 ≤ η ≤ 0.5. The answer is
    the Preliminary Design Report's own figure, and Table 8 now carries it
    η by η.
    *Engage:* none needed.
24. **Report 4 §7's ×1.75 relocation lever and its 30–60× dispersive
    shortfall** (review finding F254). The lever list of §7
    (`reports/nanowire_far_forward.template.html`) states that "detector
    relocation along the far-forward line is worth ≤ ×1.75" and that
    "dispersive tagging at IP6 is 30–60× short of the 0.1 rigidity
    threshold of Gamage et al. [12]". Neither has a producer anywhere in
    the repository. The neighbouring ×230 does re-derive
    (`farforward.hole_acceptance` at 5 × 41, p_ion = 244.8 GeV, B = 50:
    acc(1.76 mrad)/acc(2.20 mrad) = 227.7); the ×1.75 reproduces in no
    unit — read as a tagged-fraction factor it corresponds to a 1.9%
    relocation and contradicts §2's ×77 for the same move at 5 × 41, and
    read as a transport lever `farforward.POT_LEVERS` (19.24 / 21.25 /
    29.97 m) gives 1.56 over the whole instrumentable stretch. For the
    band, the report's own §1 rigidity for the α, R = 0.99813, fixes
    0.1/(1 − 0.99813) = 53×, the middle of it; the intact ⁶Li sits at
    R = 1 exactly and the ⁷Li α at R = 0.856 is 0.144 off rigidity, so the
    band is not a spread over the fragments, and arXiv:2105.13564 is not
    in `refs/`, so its Table 1 cannot be checked here. Appendix A dates
    the 30–60× to the 2026-08-27 lever table; the ×1.75 has no
    revision-history trace at all.
    *The question:* "Two: (1) What quantity is '≤ ×1.75' a factor in, and
    what produces it? As a tagged-fraction factor it corresponds to a
    1.5–1.9% relocation and contradicts §2's ×77 for the same move at
    5 × 41; as a transport-lever factor the measured R₁₂ (19.24 / 21.25 /
    29.97 m) give 1.56 for the whole instrumentable stretch, not 1.75 —
    and R₁₂ is a matrix element, not a distance, so a relocation lever
    needs the twiss rather than a ratio of z. Give it a producer (a manual
    row under §4.6 with the command, or a plans/09 derivation) or drop it.
    (2) Which entries of Gamage et al. (arXiv:2105.13564, not in refs/)
    Table 1 set the 30 and the 60? The report's own §1 rigidity,
    R = 0.99813 for the α, fixes 0.1/(1 − 0.99813) = 53×, and neither
    endpoint of the band appears anywhere in the repository."
    *Engage:* the authors; the ePIC FF WG if the ×1.75 is a transport
    bound somebody measured.
    *Default:* the text stands as written until the source is located.
    *2026-09-06:* arXiv:2105.13564 read in the arXiv v2 and the JACoW
    version of record (IPAC'21 TUPAB040, pp. 1435–1438): its title is
    "Design Concept for the Second Interaction Region for Electron-Ion
    Collider" (Report 4 [12] corrected), and its Table 1 row "Minimum
    Δ(Bρ)/(Bρ) allowing for detection of p_T = 0 fragments" reads 0.1 at
    IR 1 (IP6) and 0.003–0.01 at IR 2, which supports a 10–33× band, not
    30–60×; the paper carries no per-species efficiency. Of the
    repository's fragments, 0.1/|1 − R| is 53.6 (⁶Li α), 22.1 (⁶Li d), ∞
    (intact ⁶Li), 0.69 (⁷Li α), 0.35 (⁷Li t) — no pair spans 30–60. The
    ×1.75 remains without a producer (candidates computed: R₁₂
    29.97/19.24 = 1.56, R₃₄ 4.56/2.93 = 1.56, station z 1.05, RP/OMD z
    1.28, 10σ envelope over silicon at 18 × 275 0.92/0.53 = 1.74; the ×230
    and ×77 neighbours reproduce). Both numbers stand as written until the
    authors decide.
    *Recommendation (2026-09-16):* drop the ×1.75 relocation lever from Report 4
    §7 rather than give it a producer — no quantity in the repository produces
    it as a relocation lever (the nearest arithmetic match, the 18 × 275
    10σ-envelope-over-silicon ratio 0.92/0.53 = 1.74, is an aperture ratio and
    not a lever), and read as a tagged-fraction factor it contradicts §2's ×77
    for the same move at 5 × 41 — and replace the 30–60× by the sourced 10–33× band,
    which is what Table 1 of arXiv:2105.13564 supports (0.1 at IR 1 against
    0.003–0.01 at IR 2). The neighbouring ×230 and ×77 reproduce and stay.

## Still open — author decisions raised by run 19 (2026-09-15)

Five questions the run-19 work opened and could not answer for itself.  Each
is a published number that is *not wrong* — every one of them is labelled
with the footing it stands on — but whose footing the run has now priced, so
the authors can choose the lead.  Each carries the question, whom to engage,
and the default the repository proceeds on meanwhile.

*2026-09-16: dispositioned. One of the five is DECIDED (#26); the other four
carry a Recommendation line of their own below, as does #24 of the section
above. A recommendation is not a decision — the defaults still stand and no
published number moves on one.*

25. **Which column of §7.1 leads: toy or grid.**  Run 19 put the nuclear
    grid into production (`inputs.get_backends(pdf, nuclear, r_func)`,
    `structure.NuclearF2FromGrid` on EPPS21nlo_CT18Anlo_Li6) and measured
    the whole toy → grid drift of the four truth-level drivers (plans/07
    **WP1 addendum**, this run).  The obstacle is not the drift but the
    support: **EPPS21nlo_CT18Anlo_Li6 begins at Q = 1.3 GeV, Q₀² = 1.69 GeV²,
    and returns NaN below it**, while the money maps start at Q² = 1 GeV²
    and **two of the four published sweet spots sit at Q² = 1.14 GeV²**.
    `NuclearF2FromGrid` therefore freezes F₂ᴬ at Q₀² below the floor and
    exposes `q2_min` and `q2_frozen_fraction`; at 10 × 99.5 GeV/u, 8.7% of
    the accepted cells and **36.3% of the accepted one-year rate** lie
    below Q₀².  So spots 1 and 2, the Q² = 1.14 Δ slice and the **whole
    coherent best super-bin** (Q² ∈ [1, 1.66]) are a frozen-Q₀
    continuation, not a grid evaluation.
    *The question:* adopt the grid column as the lead of Report 1 §7.1 and
    mark its frozen-Q₀ rows; keep the toy as the lead with the grid beside
    it; or restrict the grid claim to Q² > Q₀² and re-site the sensitivity
    box above 1.69 GeV², which costs the two low-Q² sweet spots and the
    coherent channel's best bin.
    *Engage:* the authors; the EPPS21 authors or nCTEQ if a nuclear set
    with a lower Q₀ for A = 6 can be had — that, and nothing else, is what
    would close the gap rather than label it.
    *Default:* the toy column leads, the grid column is published beside it
    with its frozen-Q₀ rows daggered, as the addendum table has them.
    *Recommendation (2026-09-16):* lead with the toy column and restrict any
    statement made on the grid to Q² > 1.69 GeV², until a nuclear set with a
    lower Q₀ for A = 6 exists — the frozen-Q₀ continuation is a labelled
    extrapolation and should not carry a claim, and 36.3% of the accepted
    one-year rate, both low-Q² sweet spots and the whole coherent best
    super-bin sit under the floor.
26. **Whether the coherent tagged-yield curve is suppressed by default.**
    `polligen/coherent.py` gained `t_min_coherent(x_P, M_A)` and
    `CoherentScenario.t_min_suppression` this run, so exp(−B t_min) now
    follows the B band instead of a hard-coded constant, and `a2_tagged`
    gained `rate_weighted=False`, which applies `RATE_WEIGHT_SYST` (×0.73)
    in place of asking the caller to multiply by hand (plans/07 **WP5**
    clause, this run).  **Both default to off**, so money plot 6 and the
    WP5 scan regenerate bit-for-bit.  The reason the fold is not simply
    switched on is that **a scalar cannot represent it**: at x_P = 0.01,
    the coherence half-point, the suppression is −11.7 / −14.4 / −17.0% at
    B = 40 / 50 / 60, but |t_min| ∝ x_P², so at the x_P ≈ 0.02 edge of the
    window the same B = 50 gives **−47%**.  Rate-weighting it over the
    window does not fix that either: f_coh with x_P flat gives −11.0% and
    the DIS-like dN/dx_P ~ 1/x_P gives −2.5%, so the weight is a statement
    about the x_P spectrum and not a number.  What is robust is only the
    one-sidedness — the true yield is below the unsuppressed one
    everywhere.
    *The question:* does the central tagged-yield curve become (i) the
    unsuppressed one with a one-sided −2.5 … −47% systematic stated, (ii)
    the exp(−B t_min) one at a named reference x_P, with the band drawn,
    or (iii) a curve in x_P rather than a number, which is a figure change
    in money plot 6 and in the WP5 panel?  And separately: is the ×0.73
    rate weighting part of the central value or a systematic beside it?
    *Engage:* the authors.  (The script-side constant
    `coherent_optics_scan.py:214` `t_min = 3.2e-3` should call
    `t_min_coherent(0.01)` = 3.169×10⁻³ whichever way this goes; that
    re-quote belongs with money plot 6 and is not done.)
    *Default:* both folds stay off and the published curves are the
    unsuppressed ones, as they are today.
    ☑ **Decided (2026-09-16).** The coherent tagged-yield curve stays
    UNSUPPRESSED by default, and exp(−B t_min) and the ×0.73 rate weighting are
    quoted as a BAND beside it — −2.5 … −47% across the window — rather than
    folded into the central value. The reason is the one the item establishes:
    the suppression is not a scalar, |t_min| ∝ x_P² running it from
    −11.7 … −17.0% at the x_P = 0.01 coherence half-point to −47% at the
    x_P ≈ 0.02 edge of the window, while rate-weighting is a statement about the
    x_P spectrum (−11.0% with x_P flat, −2.5% for a DIS-like 1/x_P) and not a
    number. Until a diffractive model supplies that x_P spectrum — #18, the same
    ask — a central curve folded with either factor would publish an assumed
    spectrum as a measurement, where the band publishes what is robust: the true
    yield is below the unsuppressed one everywhere. Both flags therefore stay
    default off, money plot 6 and the WP5 scan stay bit-for-bit, and the band is
    what the text quotes beside the curve. The one re-quote this does not settle
    stands as the item records it: `coherent_optics_scan.py:214`'s hard-coded
    `t_min = 3.2e-3` should call `t_min_coherent(0.01)` = 3.169×10⁻³, which
    belongs with money plot 6.
27. **Whether the coherent 8–11σ is re-quoted over the charge-sector
    band.**  `coherent.CoherentScenario` ships `eps_b0` = −0.08 with the
    band −(0.04 … 0.13), anchored on Mäntysaari et al.'s deuteron
    ΔB₀/B ≈ 0.21 and argued in the `coherent.py` docstring as a deliberate
    *no-cancellation* scenario: the gluonic deformation need not share the
    charge sector's cancellations, and whether it does **is the
    measurement**.  The sibling generator has now priced that scenario:
    −0.08 at `slope_b` = 50 implies a ⁶Li charge quadrupole of
    **−0.9345 fm², 11.4× the measured −0.0818** — **11.6× the −0.0806
    this tree carries**, the quadrupole being the only difference —
    and run backwards the measured Q gives `eps_b0` = −0.0070 there and
    **−0.0069** here, so the charge-sector band is
    **−(0.0069 … 0.0527)** in this tree (−(0.0070 … 0.0527) in the
    sibling) — of which this
    repository's floor of 0.04 contains only the α + d model row 0.0527,
    with the GFMC 0.0171 and the measured 0.0069 both beneath it.  The
    a₂ amplitude is linear in `eps_b0`, so at the measured value the
    coherent shape term is 11.6× smaller and the significance falls by the
    same factor at fixed luminosity.
    *The question:* do the two sites that carry the figure — Report 1's
    abstract and Report 0 §5.2 (it is Report 0, not Report 1, whose §5.2
    states it; Report 1 §5.2 is *Reconstructed level* and Report 1's own
    band lives in §6.3) — keep 8–11σ as the headline with the
    no-cancellation premise stated, or quote it over the charge-sector
    band as well — i.e. "8–11σ on the gluonic scenario, under 1σ/yr if the
    gluonic deformation tracks the charge quadrupole" — which is the honest
    span of what is known?  (Report 0's Table 3 carries the coherent row's
    5σ-at-0.9% floor rather than this figure, and would not move with it.)
    Note that this 8–11σ and the sibling's O5 J/ψ a₂ reach, 2.84–3.29σ,
    are different quantities and must not be conflated, though both ride
    this constant.
    *Engage:* the Mäntysaari–Schenke group — #18's ask, which is the only
    thing that would replace the scenario with a calculation.
    *Default:* `eps_b0` = −0.08 stays the default, as it does in the
    sibling tree and for the same reason, and both σ statements stand as
    written.  The half of this that is independent of the call was done on
    2026-09-15: the arithmetic — −0.08 implies Q_charge(⁶Li) = −0.93 fm²,
    11.6× the measured value, and inverting the moment gives the band — is
    now stated in the `eps_b0` docstring of `coherent.py` and in Report 1
    §6.3, where the band it prices is stated, so the distance from the
    charge sector is a number rather than an adjective wherever the band
    appears.  What is left is the σ re-quote itself.
    *Recommendation (2026-09-16):* keep the 8–11σ as the SCENARIO's number, on
    the no-cancellation premise stated where it appears, and state the charge
    sector's counterpart beside it — that if the gluonic deformation tracked the
    measured charge quadrupole the same year would be an UNDER-1σ measurement:
    the a₂ amplitude is linear in `eps_b0` and the measured value is 11.6× below
    the scenario, so 8–11σ becomes 0.7–0.9σ at the same luminosity. The pair is
    the honest span of what is known, and neither number replaces the other.
28. **Whether Report 0's polarized-EMC reach is re-led with the grid leg.**
    `polli_fastsim.polarized.ToyG1.a1n(x) = −0.07(1 − x)² + 0.8x^2.2`
    crosses zero near x ≈ 0.25 and is positive above it, while
    NNPDFpol1.1's g₁ⁿ stays negative to x ≈ 0.6: **the toy neutron carries
    the wrong sign over roughly 0.25 < x < 0.6** (measured at Q² = 10:
    toy / NNPDFpol g₁ⁿ = +0.0052 / −0.0274 at x = 0.30 and
    +0.0078 / −0.0004 at x = 0.50).  Report 0 Table 3 publishes both legs
    of the polarized-EMC row — δΔR = 4.2 / 4.0 / 6.0 / 18.7% on the grids
    and 4.8 / 5.1 / 6.2 / 12.2% on the toy inputs, best-window bin
    0.45σ / 1.43σ against 0.38σ / 1.22σ — and the **published**
    `money_polemc` PNG is the `--pdf toy` one.  Two of the four x points,
    0.28 and 0.45, are inside the flagged window, and the target-mass
    systematic of the same row is inverse-variance weighted on the
    toy-input δΔR.  Nothing is unlabelled, so this is a choice of lead and
    not a correction; the labelling half was applied on 2026-09-15 (every
    toy-input site in Report 0 now carries the one-line fact).
    *The question:* does Table 3's row lead with the grid δΔR and carry the
    toy as the cross-check, and is the published PNG redrawn on
    `--pdf grid` under its own stem — or does the toy stay the lead on the
    ground that the grid leg mixes NNPDFpol1.1's three flavours against
    `PartonF2`'s five in every g₁/F₁ it forms?  (That mixing is the
    standard treatment and is measured, not a defect — see the `PartonF2`
    docstring — but it is the one argument for the toy lead that is about
    physics rather than habit.)
    *Engage:* none needed; both legs exist and both are cheap to run.
    *Default:* the toy leads and the grid is quoted beside it, as today,
    with the sign caveat now stated at every toy-input site.
    *Recommendation (2026-09-16):* re-lead Report 0 Table 3's polarized-EMC row
    with the GRID leg and keep the toy leg published beside it, labelled. The
    toy neutron carries the wrong sign of g₁ⁿ over roughly 0.25 < x < 0.6 and two
    of the four x points sit inside that window, so the grid leg is the one whose
    inputs are right where the row is read; the flavour-mixing argument for the
    toy lead is about a standard and measured treatment, not a defect. The
    published `money_polemc` PNG follows the lead under its own `_grid` stem, and
    the toy stem is not overwritten.
    *Note:* the exposure is isotope-dependent and Report 0's row is the
    small case — ⁷Li weights the neutron 23× less than the proton, ⁶Li
    weights the two equally.  On a ⁶Li map the swap moves the rate-weighted
    |A_∥| by ×0.61 / ×0.23 / ×0.13 and flips its sign in 13.3 / 9.1 / 2.5%
    of the accepted rate (plans/02 step 1.2, this run).
29. **Whether Report 0's ⁶Li α-tag pair is quoted as a VMC band.**  §5.4
    and Table 3 publish the S + D generator's tagged fractions as
    **2.4%** (Yellow-Report high-acceptance optics) and **25%** (tagging
    optics) at 10 × 99.5, point values on the Hulthén α–d density at
    β = 0.30.  Measured against the ANL VMC α + d overlap (#15, delivered
    2026-09-15), the same pair moves to **3.2%** (×1.32) and **24.4%**
    (×0.975) — in **opposite directions**, because the VMC density is
    harder than the Hulthén form where the Roman Pots look and softer in
    the far tail.
    *Measured in this tree (2026-09-16):* the ×1.32 above is the sibling
    generator's ratio applied to the published point.  The ANL α + d tables
    are now read by this generator (`fastsim/polli_fastsim/data/vmc`,
    `tagged.li6_alpha_channel(wave='vmc')`, `money_tagged_azz.py
    --cluster-wave vmc`), so the other end of the pair is measured rather
    than transferred: at `--config 1 --events 400000 --seed 20260713` the
    α-tag fraction reads **0.0336** against 0.0241 at the Yellow Report
    optics (**×1.394**) and **0.2476** against 0.2542 at the tagging optics
    (**×0.974**), and the sample-free model integrals move 0.024675 →
    0.033810 (×1.370) and 0.254028 → 0.247701 (×0.975).  Quoted as the
    published pair is quoted, that is **3.4%** and **24.8%**, not the 3.2%
    and 24.4% of the transfer — the tagging end moves because 24.4% was the
    ratio applied to the rounded 25%, where 24.8% is the measured point
    itself (0.2476 against the published 0.2542 = 25.4%, quoted as 25%).
    *The question:* are the two published tags re-quoted as bands with the
    VMC values as their other end, and if so, is the band written as the
    two-sided pair measured above, or does §5.4 adopt the VMC
    point?  The plans/05 rule as written — "always run the
    β = 0.20/0.30/0.40 band until VMC lands" — cannot be the answer here,
    since the β band brackets the VMC density at neither end.
    *Engage:* R.B. Wiringa (ANL, local), as #15; the tables themselves are
    already obtainable.
    *Default:* the point values stand, as published, with #15's dated
    clause carrying the in-house measurement.
    *Recommendation (2026-09-16, measured):* quote the pair as a two-sided
    band — **2.4–3.4% and 24.8–25%** — on the in-house measurement above,
    not on the transfer.  The two optics move in opposite directions, so
    neither the VMC point nor the β = 0.20–0.40 band can stand for both:
    the β band spans ×3.08 in tag acceptance and holds the k = 0.325 GeV/c
    asymmetry inside −0.889 … −0.824, where the VMC wave function leaves
    the sign alone and divides the magnitude by five.  A band whose other
    end is a different wave function is the only honest form, and the
    condition the earlier wording attached to it — that the waves be in
    this tree rather than reached through the sibling — is now met.

## Tracking

| # | item | status | next action |
|---|---|---|---|
| A1–A4 | answered | adopt | encode in fastsim + plans (done); verify with owners |
| 1–2 | P survival, transverse@IP | open | EPIOS contacts; aim INT 2027 |
| 3–5 | tensor ops, lumi, polarimetry | open | EPIOS; state assumptions in every plot |
| 6 | ⁶Li P convention | ☑ **closed 2026-08-29** | cluster picture, 0.81123 whole-nucleus from the tagged sector's own wave function (`beams.LI6_CLUSTER_POLARIZATION`), band 0.81–0.85 against the ab-initio VMC reading 0.848; the naive 1/3 stays reachable and pinned.  A VMC α–d overlap would sharpen P_D^{α−d} — that is #15; 2026-09-15: the sibling tree quotes the same central value with the wider band 0.81–0.905 (the cluster product's own spread, against this one's cluster-vs-ab-initio bracket) — a labelling call recorded in the item, nothing computes with the top edge |
| 7 | BeAGLE access+guidance | open | SDCC/ifarm accounts + email authors — **long pole, start now** |
| 8–10 | theory inputs | open | Cloët/Cosyn/Miller engagement |
| 11 | FF transfer matrices + optics at Li rigidities | **split** | optics half → plans/10 D1–D3 (C-AD, provisionally answered in §10.3); transfer matrices → plans/09 D3 (ePIC FF WG) |
| 12–13 | software checks | scheduled | inside Phase-1/2 steps |
| 14–16 | generator theory inputs (plans/05) | open; **15 delivered 2026-09-15, α + d adopted 2026-09-16** | Cloët/Cosyn (14), Wiringa (15), Cosyn/Sargsian (16).  #15's ANL VMC α + d / α + t overlaps exist and ship in `LiPolGen/data/vmc` (ANL's own tables); **2026-09-16: the α + d pair is adopted into this tree** — `fastsim/polli_fastsim/data/vmc/li6_ad1.momentum` and `li6.ad`, the same raw served bytes, read by `tagged.li6_vmc_tables` and selected with `tagged.li6_alpha_channel(wave='vmc')` or `money_tagged_azz.py --cluster-wave vmc`, opt-in and with every published number still on the Hulthén default, bit for bit.  The sign of ψ₂/ψ₀ is now MEASURED here and no longer quoted: −1 below the α–d S node at 0.1338 GeV/c, +1 between, −1 above the α–d D node at 0.4439 GeV/c, both nodes from `li6.ad`'s own signed columns.  The β = 0.20–0.40 band still brackets the VMC α–d density at neither end and the tagged fractions still move in opposite directions on the two optics, now measured in-house (×1.394 and ×0.974) — the quoting call is #29.  No α + t table is committed: that channel is a lone L = 1 wave with no observable relative phase, and `li7_alpha_channel(wave='vmc')` refuses |
| 17 | HepMC3 ion-spin convention | **schema delivered 2026-09-15**; proposal outstanding | the attribute schema is written and implemented in `LiPolGen/docs/HEPMC3_CONVENTION.md` (GenEvent-level `spin_*` / `P_z` / `P_zz` / `lam_e` / `spin_axis_*` attributes, per-particle `pol`, Asciiv3 GeV/mm, 10-digit nuclear PDG on a status-4 beam particle); plans/05 step 5.D → ePIC MC group is now the ask |
| 18 | coherent-⁶Li diffractive model (plans/06) | open | small-x theory engagement; scenario bands until then |
| 19 | RP Z-ID for A/Z = 2 (plans/06) | **redirected** | ask the incumbent: EICROC's 8-bit charge over 4 AC-LGAD planes, where one bit costs only ×1.4 in α fake rate (3.1 vs 2.3 × 10⁻⁵ at 95% ⁶Li efficiency) and the nanowire loses on fill factor instead (plans/09 D1, §9.2); and α + d is two hits 30–77 pixels apart whose second fragment vetoes 84% of the α fakes at the tagging optics (plans/09 B4) |
| 20 | RP cutout geometry + Li divergence (reco note) | **measured** (2026-08-28, plans/09 B1, `epic-main` 9aaa2969); optics half → plans/10 | re-run `tools/fullsim` if the pot geometry moves again; ePIC FF WG / C-AD for the light-ion optics |
| 21 | hadronic-method δy/y at y ≈ 0.01–0.05 (reco note) | open | ePIC inclusive WG; e+Li sample through eic-smear |
| 22 | ξ_p ≤ 0.015: published beam-beam limit or our own bound? (2026-09-02 review, F226) | ☑ **closed 2026-09-06** | cited (Report 3 [11], Report 4 [15]) — the EIC CDR's hadron beam-beam design limit, conservative for Z/A = ½ |
| 23 | the 0.6% end of the ePIC σ_p/p at 1 GeV/c (2026-09-02 review, F237) | ☑ **closed 2026-09-06** | restated (Report 3 Table 8, reco.py) on the ePIC Preliminary Design Report Fig. 8.9; the 0.6% was a 2023 non-ePIC design |
| 24 | the ×1.75 relocation lever and the 30–60× band (2026-09-02 review, F254) | **author decision** | name the quantity each factor multiplies and what produced it, or withdraw them; the text stands as written until then; 2026-09-06: Table 1 of the source supports 10–33×, the ×1.75 has no producer — see the item; **recommendation 2026-09-16**: drop the ×1.75 from Report 4 §7 and replace the 30–60× by the sourced 10–33× band |
| 25 | toy vs grid as the lead of §7.1, and the EPPS21 Q₀² = 1.69 GeV² floor (run 19, plans/07 WP1 addendum) | **author decision** | pick the lead; 36.3% of the accepted one-year rate and the whole coherent best super-bin are a frozen-Q₀ continuation and must be quoted as such either way; **recommendation 2026-09-16**: lead with the toy column and restrict any grid statement to Q² > 1.69 GeV² until a nuclear set with a lower Q₀ exists |
| 26 | exp(−B t_min) and the ×0.73 rate weighting as the CENTRAL coherent tagged yield (run 19, plans/07 WP5) | ☑ **decided 2026-09-16** | the curve stays UNSUPPRESSED by default and the two folds are quoted as a band beside it, −2.5 … −47% across the window, until a diffractive model supplies an x_P spectrum (#18); both flags stay default off and every published curve is bit-for-bit. Residue: `coherent_optics_scan.py:214`'s `t_min = 3.2e-3` should call `t_min_coherent(0.01)` = 3.169×10⁻³, with money plot 6 |
| 27 | `eps_b0` = −0.08 is 11.6× the ⁶Li charge quadrupole; re-quote the coherent 8–11σ (Report 1's abstract, Report 0 §5.2) over the charge-sector band −(0.0069…0.0527)? | **author decision** | the default stays, as in the sibling tree; the arithmetic (−0.08 ⇒ Q_charge = −0.93 fm²) is stated in `coherent.py` and Report 1 §6.3 since 2026-09-15, the σ re-quote is the open half; **recommendation 2026-09-16**: keep the 8–11σ as the scenario's number and state the charge-sector band's 1σ-level counterpart beside it |
| 28 | the toy g₁ⁿ's wrong sign over ≈ 0.25 < x < 0.6; re-lead Report 0 Table 3 and the `money_polemc` PNG with the grid leg? | **author decision** | labelling applied 2026-09-15 at every toy-input site in Report 0; both legs published, toy still the lead; **recommendation 2026-09-16**: re-lead Report 0 Table 3 with the grid leg and keep the toy leg labelled beside it |
| 29 | Report 0 §5.4 / Table 3's ⁶Li α-tag pair (2.4% / 25%) as a VMC band? | **author decision** | **2026-09-16: measured in this tree**, the α + d tables having landed here — 3.4% and 24.8% (×1.394 and ×0.974), opposite directions on the two optics, against the 3.2% / 24.4% the sibling transfer estimated; the point values stand until the authors choose; **recommendation 2026-09-16**: quote the pair as the two-sided band 2.4–3.4% and 24.8–25% |
