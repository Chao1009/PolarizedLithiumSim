# Open Questions & External Dependencies

Items that gate or shape the simulation program but are *not* solvable
inside it. Each has an owner-to-engage and a default assumption we proceed
with until answered. Updated 2026-08-25 (first version 2026-06-12) after the fetch-verified literature
sweep — several items moved from "unknown" to "answered, needs adoption".

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
    §2.2): ⁷Li α → Roman Pots (R = 0.86); tritons → **no IP6 coverage**
    (R > 1, both beams); ⁶Li α/d → beam-blind below RP pT cutoff (R = 1.0);
    p → OMD; n/γ → ZDC. IR-8 secondary focus (RPs 44–45.5 m) recovers R ≈ 1
    at pT → 0. *Remaining question:* none at concept level — quantitative
    acceptance is exactly Phase-2 step 2.2.
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
   `money_cos2phi_reco.py --eff-cos2-split`; plans/08 A1).  The coherent counterpart is
   tighter than the naive estimate: 10⁻³ of the Roman-Pot vertical
   envelope between the samples biases a_t by **19%** and 1% by 169%
   (measured, plans/08 A1b — under the slot half the β bins are blind and
   the t-template is 99% anti-correlated with the constant, so a shape
   perturbation is amplified ~100× over δ⟨cos 2β⟩/(P₊ − P₀) = 1.3%).
   Bunch-by-bunch alternation is therefore a requirement of the
   measurement; fill-by-fill needs 10⁻⁴ stability of both the φ′
   efficiency and the pot envelope.  The higher-|t| bins are the
   fallback: +3.9% and +0.04% at the same 10⁻³.
   *Default:* equal thirds (+,0,−), δ(rel-lumi) as a Phase-2 systematic.
4. **Li luminosity.** Confirmed gap — no Li number exists in any document
   (EPIOS included). Space charge, IBS, cooling for Li bunches unstudied.
   *Engage:* EPIOS/C-AD. *Default:* 10 fb⁻¹/nucleon per setting, quoted
   ∈ {1, 10, 100}.
5. **Li ring polarimetry.** EPIOS concept: Li–Li elastic CNI vs polarized
   Li jet (HJET analog) + Breit–Rabi absolute; analyzing-power theory
   flagged as needing work. R&D scale in EPIOS: ~26 FTE-yr/12 yr.
   *Default:* δP/P = 3% systematic in FOMs.

## Still open — generator / theory

6. **⁶Li effective-polarization convention** (factor 2.4 in the g₁ FOM):
   slides' P_p = P_n = 1/3 (per-nucleon-normalized) vs cluster-model
   whole-nucleus ≈ 0.81 (= 0.87 P_d × 0.93 D-state). ⁷Li is settled:
   P_p = +0.866, P_n = −0.037.
   *Engage:* I. Cloët (ANL, local). *Default in code:* 1/3 (conservative),
   switch after resolution.
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
9. **b₁ and Δ theory for ⁶Li specifically.** Confirmed literature gaps:
   no b₁ prediction for any A > 2; no EIC Δ projection for any target.
   First-mover opportunity — co-author with theory.
   *Engage:* Cloët (ANL), Cosyn, Miller; lattice: Detmold/Shanahan.
10. **Radiative corrections on tensor observables** (A_zz, cos 2φ).
    Vector-case tools exist (DJANGOH/HERACLES); tensor RC uncharted.
    *Default:* RC band from vector-case studies (Phase-1 step 1.4).
14. **Complete inclusive structure-function basis for spin-3/2** (⁷Li):
    rank-2 (b₁-analog) and rank-3 functions are not classified anywhere we
    can adopt; needed by the doubly polarized generator (plans/05 §5.2).
    *Engage:* Cloët, Cosyn — co-author opportunity.
    *Default:* rank ≤ 2 truncation, rank-2 shapes as scenarios.
15. **VMC two-cluster overlaps with m-dependence** (α+d S/D for ⁶Li,
    α+t P for ⁷Li) to replace the two-parameter radial forms whose tail
    dominates tagged acceptances (e+d control: BeAGLE tails 2–13× Hulthén).
    *Engage:* R.B. Wiringa (ANL, local). *Default:* β = 0.20–0.40 band.
16. **FSI for cluster spectators beyond IA** (α rescattering off DIS
    debris; deuteron case solved in PRC 97:035209, pole extrapolation).
    *Engage:* Cosyn, Sargsian. *Default:* IA, quoted at small |t′|.
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
    plans/06.

## Still open — detector / software

11. **Far-forward transfer matrices & optics at Li rigidities** (RP/OMD
    reconstruction in EICrecon is tuned per beam setting).
    *Engage:* ePIC FF WG (A. Jentsch). Phase-2 step 2.2.
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
    *2026-08-26 (plans/09, `reports/nanowire_far_forward`):* **a
    candidate answer exists, and it is not the dE/dx one.** A
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
    *Blocking measurement:* the α turn-on curve on the same wires as the
    120 GeV proton (plans/09 D1) — their α analysis is "underway", and a
    relativistic ⁶Li is an INTERPOLATION between their two measured
    points (134 nm proton, ≈1 µm for a 5.5 MeV ²⁴¹Am α).
20. **Roman-Pot cutout geometry and Li beam divergence** (reconstruction-
    chain note, 2026-08-24). The near-beam cut is angular,
    pT_cut = 10σ_θ·A·p_u, so the same optics gives tag acceptances of
    67% / 9% / 10⁻⁸ at 20.5 / 50 / 137.5 GeV/u for ⁶Li (constant 0.20 GeV
    in the code → 13.5%); and the pots' rectangular cutout has its sides
    parallel to the vertical spin axis, faking ⟨cos 2φ_t⟩ ≈ 0.5 for a 25%
    aspect ratio (physics a₂ ≈ 0.036) unless the spin-state ratio is
    used. *Needed:* σ_θx, σ_θy (β*, emittance) at the light-ion energies
    and the pot geometry (`reco.rp_measure` takes both). *Engage:* ePIC
    FF WG / C-AD optics. *Default:* proton-derived 73/149 μrad,
    slot-like cutout (2.5 : 1, `rp_measure(cut_scale_xy=(2.5, 1))`), angular
    envelope 10σ_θ·A·p_u; ratio estimator.
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
    onto the pot plane with R₁₂ ≈ 30.6 m horizontally against a few
    metres vertically, so what clears the slot is the HORIZONTAL angle:
    the boundary is |θ_x| ≳ 2.0 / 1.35 / 1.03 mrad in the 5×41 / 10×100 /
    18×275 optics (p_T = 0.25 / 0.41 / 0.85 GeV for the ⁶Li at that
    ring's rigidity), against |θ_y| ≳ 1.8–3 mrad. In `rp_measure` terms
    that is `cut_scale_xy ≈ (1, 2.3)`, not (2.5, 1) — a factor ≈ 5.8 the
    wrong way — and the fake ⟨cos 2β⟩ about the vertical spin axis is
    therefore large and **negative**, not positive. The tagged fraction
    falls with it (B = 50 GeV⁻²): 1.4×10⁻² against 0.82 at 20.5 GeV/u,
    5×10⁻⁵ against 0.40 at 50.
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
    **So the coherent programme is a low-energy programme**, which is what
    the angular envelope already implied, but for a second and stronger
    reason. *Caveats:* one event per scan point, 30° azimuthal steps,
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
    against the 1.03 measured — agreement to 1.5%, so measurement and
    file reading confirm each other. The current 16 mm block implies
    roughly half that, *below* the 0.727 mrad 10σ envelope; at 5 × 41 the
    per-energy insertion moves the other way (29.6 mm inner edge).
    **Every aperture-conditional number — the tagged fractions, the
    sign-flipped ⟨cos 2β⟩, "the coherent programme is a low-energy
    programme" — is therefore conditioned on a superseded geometry and
    must be re-measured** (plans/09 B1, priority). What a closer approach
    is WORTH is unaffected and now curved rather than tabulated:
    `nearbeam_aperture_scan.py` prices every aperture, and the chain says
    the gain is ×26 at 5 × 41, ×569 at 10 × 100, and at mid energy the
    difference between a_t = 0.48 ± 1.45 and four clean bins.
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
    δy/y = **0.32 / 0.22 / 0.28 / 0.11** at the sweet spots with the 50 MeV
    floor — the toy was optimistic by 0.04–0.05 absolute at every one,
    because it put 0.03 of Σ into neutral hadrons where PYTHIA puts 0.11 —
    and 0.21 / 0.12 / 0.17 / 0.10 at LOW against 0.74 / 0.34 / 0.69 / 0.18
    at TOP. **Half of this question is therefore answered**: what remains
    is only the ePIC noise/threshold floor, which is what the scan
    0 → 25 → 50 → 100 MeV turns into 0.20 → 0.32 → 0.54 → 1.01 at y = 0.005.
    Reco purity at the sweet spots falls from 0.64–0.68 to 0.40–0.73 with
    the real final state.
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

## Tracking

| # | item | status | next action |
|---|---|---|---|
| A1–A4 | answered | adopt | encode in fastsim + plans (done); verify with owners |
| 1–2 | P survival, transverse@IP | open | EPIOS contacts; aim INT 2027 |
| 3–5 | tensor ops, lumi, polarimetry | open | EPIOS; state assumptions in every plot |
| 6 | ⁶Li P convention | open | ask Cloët (local) — **first, cheapest, factor-2.4 impact** |
| 7 | BeAGLE access+guidance | open | SDCC/ifarm accounts + email authors — **long pole, start now** |
| 8–10 | theory inputs | open | Cloët/Cosyn/Miller engagement |
| 11–13 | software checks | scheduled | inside Phase-1/2 steps |
| 14–16 | generator theory inputs (plans/05) | open | Cloët/Cosyn (14), Wiringa (15), Cosyn/Sargsian (16) |
| 17 | HepMC3 ion-spin convention | scheduled | plans/05 step 5.D → ePIC MC group |
| 18 | coherent-⁶Li diffractive model (plans/06) | open | small-x theory engagement; scenario bands until then |
| 19 | RP Z-ID for A/Z = 2 (plans/06) | candidate answer | nanowire firing threshold, r_s ∝ z (plans/09 §9.2); blocked on the α turn-on curve (plans/09 D1) — ask ANL MEP |
| 20 | RP cutout geometry + Li divergence (reco note) | **re-measure** | the ePIC pot geometry moved after the Sep-2024 snapshot (plans/09 §9.4, B1); then ePIC FF WG / C-AD for the light-ion optics |
| 21 | hadronic-method δy/y at y ≈ 0.01–0.05 (reco note) | open | ePIC inclusive WG; e+Li sample through eic-smear |
