# Open Questions & External Dependencies

Items that gate or shape the simulation program but are *not* solvable
inside it. Each has an owner-to-engage and a default assumption we proceed
with until answered. Updated 2026-09-02 (first version 2026-06-12): the consistency review of
Reports 0–4 added items 22–24, the three passages it could not settle from
the repository. The 2026-08-25 revision followed the fetch-verified
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
    *2026-08-29:* the cost of that band is now measured on the ⁷Li P
    wave — ⟨k⟩ = 0.1113 / 0.1333 / 0.1505 GeV/c and P(k > 0.3 GeV/c) =
    0.0231 / 0.0522 / 0.0836 at β = 0.20 / 0.30 / 0.40, a factor 3.6 in
    exactly the region where the tagged tensor asymmetries are O(1),
    while the angular moment ⟨P₂⟩ = −T/5 is fixed by Clebsch–Gordan and
    does not move — and the ask is `docs/note_7li_theory_questions.md`
    §3.
16. **FSI for cluster spectators beyond IA** (α rescattering off DIS
    debris; deuteron case solved in PRC 97:035209, pole extrapolation).
    *Engage:* Cosyn, Sargsian. *Default:* IA, quoted at small |t′|.
    *2026-08-29:* the ⁷Li form of the same ask — a triton spectator, and
    the two-body currents an A = 7 cluster pair carries — is
    `docs/note_7li_theory_questions.md` §5.
17. **HepMC3 convention for ion spin states** — none exists; plans/05
    step 5.D defines named attributes and proposes them upstream.
    *Engage:* ePIC MC/software group. *Default:* our attribute schema.
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
    quotes millimetres everywhere; R₁₁, R₂₁, R₂₂ and D′ remain unmeasured,
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

## Still open — sourcing (raised by the consistency review of Reports 0–4)

Three passages the review of 2026-09-02
(`docs/consistency_review_2026-09-02.md` §2.5) could not settle from the
repository: each states a number whose producer or citation exists
nowhere in it, so only the authors can say what it is. The review left a
placeholder edit for each, unapplied.

22. **The source of ξ_p ≤ 0.015** (review finding F226). Report 3 Table 9
    row 2 (`reports/eic_epic_reference.template.html`) asks C-AD whether
    IR6 can be matched with the horizontal hadron β* raised 46–164× "and
    the electron β* co-de-squeezed so ξ_p stays ≤ 0.015"; Report 4 §7
    (`reports/nanowire_far_forward.template.html`) repeats the clause from
    there. The bound lives in those two sentences and nowhere else — no
    derivation in Report 3 §4.2, no entry in either reference list or in
    `refs/refs_dict.json`, nothing in the four Yellow-Report parts (Tables
    10.1/10.2 carry no beam-beam row), and, before this entry, no other
    occurrence of "beam-beam" anywhere in the repository — while the
    Report 3 dateline
    promises that every Table 9 entry carries its source and its status.
    *The question:* "Where does ξ_p ≤ 0.015 come from: a published EIC
    beam-beam limit for the hadron beam (which document and table?), or is
    it this programme's own assumed bound, to be labelled 'ours to
    specify' in Table 9 as rows 11, 12 and 16 are?"
    *Engage:* the authors; C-AD if it is a design limit to be cited.
    *Default:* the text stands as written until the source is located.
23. **The upper end of the ePIC momentum resolution, 0.45–0.6% at
    1 GeV/c** (review finding F237). Report 3 Table 8's tracking σ_p/p row
    reads "ePIC full simulation reaches 0.45–0.6% at 1 GeV/c, so realistic
    to ≈ ±20%", uncited in a report whose masthead promises a source for
    every entry, and the claim originated in the docstring of
    `evgen/polligen/reco.py` (line 413), the only other live copy. The one
    ePIC full-simulation momentum resolution held on disk is S. Maple,
    Birmingham seminar 11 Dec 2024, slide 47
    (`refs/EIC_Seminar_SMaple_2024.pdf` p. 47, indexed in
    `refs/refs_dict.json`): PYTHIA 8 NC DIS, Craterlake 23.12.0,
    Q² > 100 GeV², electron from tracking, 0 ≤ η ≤ 0.5, fit
    Δp/p [%] = 0.055 p ⊕ 0.45, i.e. 0.453% at p = 1 GeV/c against the
    barrel model's 0.503% — a ≈ 10% comparison, in one η slice. Nothing in
    the repository gives 0.6%.
    *The question:* "Which ePIC document gives the upper end of 'ePIC full
    simulation reaches 0.45–0.6% at 1 GeV/c depending on η slice'? Maple's
    Dec-2024 seminar slide 47 — the only ePIC full-simulation momentum
    resolution in refs/ — fits Δp/p = 0.055 p ⊕ 0.45% for the single slice
    0 ≤ η ≤ 0.5, i.e. 0.45% at 1 GeV/c and nothing at 0.6%. If a second,
    η-resolved source exists it should be added to refs/ and cited in
    Table 8; if not, Table 8 and reco.py must narrow the claim to 0.45% in
    0 ≤ η ≤ 0.5, which turns the row's ±20% into ≈ 10% there."
    *Engage:* the authors; the ePIC tracking/inclusive WG if a second,
    η-resolved source exists.
    *Default:* the text stands as written until the source is located.
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

## Tracking

| # | item | status | next action |
|---|---|---|---|
| A1–A4 | answered | adopt | encode in fastsim + plans (done); verify with owners |
| 1–2 | P survival, transverse@IP | open | EPIOS contacts; aim INT 2027 |
| 3–5 | tensor ops, lumi, polarimetry | open | EPIOS; state assumptions in every plot |
| 6 | ⁶Li P convention | ☑ **closed 2026-08-29** | cluster picture, 0.81123 whole-nucleus from the tagged sector's own wave function (`beams.LI6_CLUSTER_POLARIZATION`), band 0.81–0.85 against the ab-initio VMC reading 0.848; the naive 1/3 stays reachable and pinned.  A VMC α–d overlap would sharpen P_D^{α−d} — that is #15 |
| 7 | BeAGLE access+guidance | open | SDCC/ifarm accounts + email authors — **long pole, start now** |
| 8–10 | theory inputs | open | Cloët/Cosyn/Miller engagement |
| 11 | FF transfer matrices + optics at Li rigidities | **split** | optics half → plans/10 D1–D3 (C-AD, provisionally answered in §10.3); transfer matrices → plans/09 D3 (ePIC FF WG) |
| 12–13 | software checks | scheduled | inside Phase-1/2 steps |
| 14–16 | generator theory inputs (plans/05) | open | Cloët/Cosyn (14), Wiringa (15), Cosyn/Sargsian (16) |
| 17 | HepMC3 ion-spin convention | scheduled | plans/05 step 5.D → ePIC MC group |
| 18 | coherent-⁶Li diffractive model (plans/06) | open | small-x theory engagement; scenario bands until then |
| 19 | RP Z-ID for A/Z = 2 (plans/06) | **redirected** | ask the incumbent: EICROC's 8-bit charge over 4 AC-LGAD planes, where one bit costs only ×1.4 in α fake rate (3.1 vs 2.3 × 10⁻⁵ at 95% ⁶Li efficiency) and the nanowire loses on fill factor instead (plans/09 D1, §9.2); and α + d is two hits 30–77 pixels apart whose second fragment vetoes 84% of the α fakes at the tagging optics (plans/09 B4) |
| 20 | RP cutout geometry + Li divergence (reco note) | **measured** (2026-08-28, plans/09 B1, `epic-main` 9aaa2969); optics half → plans/10 | re-run `tools/fullsim` if the pot geometry moves again; ePIC FF WG / C-AD for the light-ion optics |
| 21 | hadronic-method δy/y at y ≈ 0.01–0.05 (reco note) | open | ePIC inclusive WG; e+Li sample through eic-smear |
| 22 | ξ_p ≤ 0.015: published beam-beam limit or our own bound? (2026-09-02 review, F226) | **author decision** | locate the document and table, or relabel the ask as ours to assume; the text stands as written until then |
| 23 | the 0.6% end of the ePIC σ_p/p at 1 GeV/c (2026-09-02 review, F237) | **author decision** | name the ePIC document the upper end comes from, or mark the range as our reading of one slide; the text stands as written until then |
| 24 | the ×1.75 relocation lever and the 30–60× band (2026-09-02 review, F254) | **author decision** | name the quantity each factor multiplies and what produced it, or withdraw them; the text stands as written until then |
