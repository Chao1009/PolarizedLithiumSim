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
    (R > 1, both beams); ⁶Li α/d → beam-blind below RP pT cutoff
    (R = 0.998 and 1.005 from the physical nuclear masses since
    2026-08-28, not the 1.000 the A·Z ratio gives — the two fragments are
    separated by 0.7% of rigidity, which is inside the ±5% near-beam band
    and so still undispersed, but they are not the same trajectory);
    p → OMD; n/γ → ZDC. IR-8 secondary focus (RPs 44–45.5 m) recovers R ≈ 1
    at pT → 0. *Remaining question:* none at concept level — quantitative
    acceptance is exactly Phase-2 step 2.2.  *Caveat on the tritons
    (2026-08-28):* "no coverage" is the routing's own answer —
    `farforward.route_charged` has no R > 1 branch, so an over-rigid
    fragment is lost by construction — and the 2026-06-12 particle-gun
    scan of `tools/fullsim` contradicts it, putting the ⁷Li triton on the
    inner side of the Roman-Pot planes and then in the ZDC. Unverified
    against a beam envelope, a divergence or a reconstruction; open as
    plans/03 §2.2 "Tritons at IP6 — revisit" ◐.
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
   a_t by **−0.8 / −0.5 / −0.1 / −0.3 %** in the four |t| bins, and 10⁻²
   by −7.6 / −4.5 / −3.6 / −2.9 %, with a_e untouched (Report 2
   Table 6).  The tagging cutout's horizontal edge sits in a shallow part
   of the recoil spectrum and all four bins stay live, so the
   amplification is gone: per-mille stability now costs per-mille.
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
    *Default:* **no band.**  An unpolarized QED study does not bound a
    tensor one, so the affected claims carry no correction and say so
    (plans/07 WP4, Report 2 §7).  What *is* measured is the unpolarized
    collinear-ISR migration, in `polligen/radiative.py`: +0.62 / +0.50 /
    +0.94 / +1.22% of Δ̂ at the four sweet spots in the published
    generator window, ≤ 2.9% once the low-Q² feed-in is opened
    (Phase-1 step 1.4, plans/08 D3).
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
    R₁₂ = 30.6 m and the pot dispersion D = 0.30 m were measured at
    18 × 275 only and are carried to the other configurations for want of
    per-optics values, R₃₄ is unmeasured, and the 10σ offsets in
    `beamline_*.xml` are marked a *"rough extrapolation"* at 5 × 41 —
    which is why plans/09 works in angle and quotes millimetres for one
    optics alone.
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
    d (2p_u) take opposite kicks and land a median 15.1 / 18.5 / 38.4 mm
    apart at 18×275 / 10×100 / 5×41 — **30 to 77 pixels of the existing
    500 µm pitch**, and 13 to 159 over the 16–84% band. An intact ⁶Li is
    one hit. And the second fragment is a *veto*: conditioned on an α
    that fakes a coherent tag, the partner deuteron is recorded in **84%**
    of events at the tagging optics (0.04–0.29 at the published ones,
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
    per axis, and the tagging optics of Report 1 §6.1 (0.33 × 3.80,
    0.17 × 1.80, 0.12 × 0.92 mrad at 1/7.1, 1/13.3, 1/9.5 of the
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
    onto the pot plane with R₁₂ ≈ 30.6 m horizontally against a few
    metres vertically, so what clears the slot is the HORIZONTAL angle:
    the boundary is |θ_x| ≳ 2.0 / 1.35 / 1.03 mrad in the 5×41 / 10×100 /
    18×275 optics (p_T = A p_u |θ_x| = 0.49 / 0.81 / 0.85 GeV for the ⁶Li
    at the γ-matched 40.8 / 99.5 / 137.5 GeV/u; the 0.25 / 0.41 / 0.85
    this entry carried until 2026-08-28 priced the same edges at the
    retired rigidity-scaled momenta and made the low configuration look
    3.4× rather than 1.7× more permissive), against |θ_y| ≳ 1.8–3 mrad.
    In `rp_measure` terms
    that is `cut_scale_xy ≈ (1, 2.3)`, not (2.5, 1) — a factor ≈ 5.8 the
    wrong way — and the fake ⟨cos 2β⟩ about the vertical spin axis is
    therefore large and **negative**, not positive. The tagged fraction
    falls with it (B = 50 GeV⁻², `evgen/scripts/nearbeam_aperture_scan.py`,
    re-run 2026-08-28 at the γ-matched momenta): **9.8×10⁻⁷ / 7.7×10⁻¹⁶ /
    1.9×10⁻¹⁷** through the measured aperture at 40.8 / 99.5 / 137.5
    GeV/u, against **7.2×10⁻⁸ / 6.2×10⁻²⁷ / 3.9×10⁻¹⁴** through the
    Yellow Report high-acceptance envelope on the scan's own convention
    (the envelope horizontally, the larger of silicon and envelope
    vertically; the pure 10σ envelope in both planes gives 7.2×10⁻⁸ /
    1.2×10⁻²⁶ / 7.8×10⁻¹⁴ — the two are different quantities sharing a
    name, `tools/fullsim/README.md`).  **Which of the two binds is the
    correction of 2026-08-28:** the measured edge is 0.91× / 0.75× /
    1.12× the envelope's horizontal half-width, so the ENVELOPE binds at
    the two lower configurations and the silicon only marginally at the
    top — the opposite of the "2.8× / 1.9× / 1.4×, the envelope is never
    binding" that plans/08 §8.4 obtained by pricing the same edges
    against the retired single 73 μrad.  The 1.4×10⁻² and 5×10⁻⁵ this
    entry carried until then were the measured aperture at the retired
    rigidity-scaled 20.5 and 50 GeV/u.
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
    (This chain entry is itself at the retired rigidity-scaled menu:
    "LOW" is 20.5 GeV/u, superseded by plans/10.)
    **The conclusion drawn here — "so the coherent programme is a
    low-energy programme … for a second and stronger reason" — is
    withdrawn (2026-08-28).** It rested on the measured aperture binding
    everywhere, which it does not, and on the rigidity-scaled energies.
    Against the per-configuration envelopes the aperture is the binding
    constraint only at the top; and at the tagging optics, where the
    channel is now measured, the ordering reverses outright — 0.42 / 0.32
    / 0.33 tagged at 1/7.1 / 1/13.3 / 1/9.5 of the luminosity, with the
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
    against the 1.03 measured — agreement to 1.5%, so measurement and
    file reading confirm each other. The current 16 mm block implies
    roughly half that, 0.52 mrad, *below* the 0.92 mrad Yellow Report
    high-acceptance envelope at 18 × 275 (plans/10) but 4× the 0.12 mrad
    tagging-optics envelope; at 5 × 41 the per-energy insertion moves the
    other way (29.6 mm inner edge).
    **Every aperture-conditional number — the tagged fractions, the
    sign-flipped ⟨cos 2β⟩, "the coherent programme is a low-energy
    programme" — is therefore conditioned on a superseded geometry and
    must be re-measured** (plans/09 B1, priority). What a closer approach
    is WORTH is unaffected and now curved rather than tabulated:
    `nearbeam_aperture_scan.py` prices every aperture per configuration.
    On the Yellow Report divergences (2026-08-28, plans/10 A4) the
    machine envelope binds at every configuration and a closer approach
    alone buys nothing; under the tagging optics of Report 1 §6.1 a layer
    that follows the 0.33 / 0.17 / 0.12 mrad envelope is the difference
    between no tag and 0.41 / 0.31 / 0.32 with four clean |t| bins each
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
    *2026-08-27 (acceptance, `hfs_acceptance.py`, Report 2 §3 Figure 4):*
    78 / 86 / 82 / 91% of Σ_h is captured at the four mid sweet spots;
    19 / 8 / 16 / 7% escapes forward beyond |η| = 3.7 (the target-
    fragmentation side of a W ≈ 6–10 GeV system), 1–6% is below threshold.
    Through the full response 69 / 74 / 73 / 85% is captured: a 13–28% scale
    bias on y_Σ that the library *reproduces* in the pseudo-events and the
    bin-centering factor absorbs; `--hfs-calibrate` (per-cell mean captured
    fraction, the analysis's own calibration) takes the 5R purities to
    0.52–0.76 at unchanged errors, and a residual 1% scale error moves Δ̂ by
    0.2–0.7%.  Not the resolution driver: the ePIC nominal reach of 4.0
    recovers a quarter of the escape with δy/y unchanged.  What the ePIC
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

## Tracking

| # | item | status | next action |
|---|---|---|---|
| A1–A4 | answered | adopt | encode in fastsim + plans (done); verify with owners |
| 1–2 | P survival, transverse@IP | open | EPIOS contacts; aim INT 2027 |
| 3–5 | tensor ops, lumi, polarimetry | open | EPIOS; state assumptions in every plot |
| 6 | ⁶Li P convention | open | ask Cloët (local) — **first, cheapest, factor-2.4 impact** |
| 7 | BeAGLE access+guidance | open | SDCC/ifarm accounts + email authors — **long pole, start now** |
| 8–10 | theory inputs | open | Cloët/Cosyn/Miller engagement |
| 11 | FF transfer matrices + optics at Li rigidities | **split** | optics half → plans/10 D1–D3 (C-AD, provisionally answered in §10.3); transfer matrices → plans/09 D3 (ePIC FF WG) |
| 12–13 | software checks | scheduled | inside Phase-1/2 steps |
| 14–16 | generator theory inputs (plans/05) | open | Cloët/Cosyn (14), Wiringa (15), Cosyn/Sargsian (16) |
| 17 | HepMC3 ion-spin convention | scheduled | plans/05 step 5.D → ePIC MC group |
| 18 | coherent-⁶Li diffractive model (plans/06) | open | small-x theory engagement; scenario bands until then |
| 19 | RP Z-ID for A/Z = 2 (plans/06) | **redirected** | ask the incumbent: EICROC's 8-bit charge over 4 AC-LGAD planes, where one bit costs only ×1.4 in α fake rate (3.1 vs 2.3 × 10⁻⁵ at 95% ⁶Li efficiency) and the nanowire loses on fill factor instead (plans/09 D1, §9.2); and α + d is two hits 30–77 pixels apart whose second fragment vetoes 84% of the α fakes at the tagging optics (plans/09 B4) |
| 20 | RP cutout geometry + Li divergence (reco note) | **re-measure**; optics half → plans/10 | the ePIC pot geometry moved after the Sep-2024 snapshot (plans/09 §9.4, B1); then ePIC FF WG / C-AD for the light-ion optics |
| 21 | hadronic-method δy/y at y ≈ 0.01–0.05 (reco note) | open | ePIC inclusive WG; e+Li sample through eic-smear |
