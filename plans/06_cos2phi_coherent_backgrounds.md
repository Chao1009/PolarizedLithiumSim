# Plan 06 — cos 2φ Money Plots: the Coherent (Intact-⁶Li) Channel and the Background Budget

**Context (2026-08-10).** Money plots 5 and 6 present the cos 2φ tensor
modulation as *projected data points with statistical error bars at the
sweet-spot (x, Q²) bins* — the presentation the INT-2027 audience will
want. Detection assumption: the scattered electron (central detector,
`fom.Scenario` cuts) plus whatever part of the lithium remnant the
far-forward systems can see; the channel of primary interest keeps the
⁶Li **intact in its 1⁺ ground state — the coherent process**.

| artifact | content |
|---|---|
| `evgen/scripts/money_cos2phi.py` → `money_cos2phi_6Li.png` | inclusive gluonometry: φ′ pseudo-data in 4 sweet-spot super-bins (Q² = 1.1→14 GeV² at x ≈ 0.02–0.06) + amplitude vs x with Δ/F₁ scenario curves |
| `evgen/scripts/money_cos2phi_coherent.py` → `money_cos2phi_coherent_6Li.png` | coherent channel: recoil t-spectrum vs RP cuts, tagged yield vs x, φ′ pseudo-data in the best low-x bin |
| `evgen/polligen/coherent.py` | scenario coherent model + intact-recoil routing + breakup veto table (tested) |

Headline numbers (mid energy e10×⁶Li50, L = 100 fb⁻¹/u, P_zz = 0.6,
scenario inputs):

- **Inclusive:** per-super-bin δA ≈ (0.4–1.2)×10⁻⁴ → the Δ/F₁ = 10⁻³
  scenario amplitude (≈1.4×10⁻⁴ at the sweet spots) is a 2.5–2.7σ
  effect per bin in the three low/mid-Q² panels (1.2σ at Q² ≈ 14), 4.6σ
  combining the four displayed bins, and 9.7σ combining the full
  accepted (x,Q²) map — consistent with the `money_delta` reach curves
  (numbers verified by the 2026-08-10 adversarial review).
- **Coherent:** the intact ⁶Li has exactly beam rigidity (A/Z = 2,
  R = 1.000), so at IP6 it is visible **only** in the Roman-Pot
  near-beam pT tail: acceptance exp(−B·pT_cut²) = **13.5%**
  [9–20% for B ∈ 40–60 GeV⁻²] with high-acceptance optics —
  and **4×10⁻⁵ with high-divergence optics: the coherent program
  fixes the optics choice**. Scenario rates give ~10⁸ tagged events at
  100 fb⁻¹/u; the best super-bin (x ≈ 0.002, Q² ≈ 1.3) reaches
  δA ≈ 6×10⁻⁴, i.e. a 5σ floor at a 0.3% modulation amplitude.

Everything below is the requested background think-through, ordered by
how dangerous the background is to the coherent tag.

---

## 6.1 What each channel actually measures

- **Inclusive (breakup) DIS**, transverse tensor fill, unpolarized e:
  the HJM double-helicity-flip Δ(x,Q²) — "exotic glue". Only e′ is
  measured; φ is the lepton-plane azimuth relative to the transverse
  alignment axis. The nucleus fragments; tagging is optional (and mostly
  beam-blind for ⁶Li anyway, plans/00).
- **Coherent diffractive DIS** e ⁶Li → e′ X ⁶Li(g.s.), same fill: the
  cos 2φ modulation of the *coherent* yield — the double-helicity-flip /
  gluonic-quadrupole analog in the nucleus-intact channel (the small-x
  "elliptic"-gluon side of the same physics; scenario amplitude until a
  real model is adopted, open question #18). The intact recoil is the
  event-defining tag, so background = anything that puts a
  beam-rigidity track into the RP pT tail with a quiet central detector.

## 6.2 Coherent-channel backgrounds, ranked

### (1) Incoherent diffraction / quasi-elastic breakup with a beam-blind fragment — THE background

Diffractive events where the ⁶Li breaks up but the surviving fragment is
indistinguishable from an intact ⁶Li in the Roman Pots. The rigidity
table from `coherent.veto_table()` (beam A/Z = 2):

| breakup channel | threshold | fragment routing |
|---|---|---|
| **α + d** | 1.474 MeV | α: R = 1.000 → RP pT tail only (beam-blind); d: R = 1.000 → same |
| α + p + n | 3.70 MeV | α beam-blind; p: R = 0.50 → OMD; n → ZDC |
| ³He + t | 15.79 MeV | ³He: R = 0.75 → RP main window; t: R = 1.50 → lost (over-rigid) |

The **α+d channel is the killer**: *both* fragments sit at exactly beam
rigidity, exactly like the intact ⁶Li — same trajectory for the same
rigidity, same velocity (so time-of-flight cannot help). The
lowest-lying strength (direct α+d continuum from 1.47 MeV and the 2.186
MeV 3⁺ resonance, Γ ≈ 24 keV, which decays ~100% to α+d with only
Q ≈ 0.7 MeV) produces an α with ~2/3 of the recoil pT plus a small
breakup smear (k_rel ≈ √(2μQ) ≈ 40 MeV/c) — precisely the events whose
α can wander above the 0.2 GeV cut while the partner d (∼1/3 of the pT)
stays invisible.

Handles, in order of power:

1. **Charge discrimination in the RPs.** dE/dx ∝ Z²: ⁶Li (Z=3)
   deposits 2.25× an α and 9× a deuteron. Research-sweep finding: **no
   EIC document addresses this A/Z = 2 degeneracy** — the RP AC-LGAD
   chain (EICROC) does record pulse amplitude, but only for
   charge-sharing position interpolation; the only *documented* Z-ID
   concept is a Cherenkov counter (N_γ ∝ Z²) behind the Roman Pots at
   the IR-8 secondary focus (YR §7.5; arXiv:2211.15746 §VIII;
   arXiv:2602.04636). Whether the IP6 sensor dynamic range separates
   Z = 3/2/1 is a genuinely open hardware question — **new open
   question #19** for the ePIC FF WG. If yes, this background is
   largely eliminated; if no, the coherent tag rests on shape fits.
2. **t-shape separation.** Coherent: e^(−B|t|), B ≈ 50 GeV⁻²;
   incoherent: slope set by the cluster wave function / nucleon size,
   an order of magnitude flatter. In the tagged window |t| ≳ 0.04 GeV²
   (pT > 0.2 GeV; identical to the Yellow-Report coherent-helium DVCS
   window, YR §8.4.3) and below the first ⁶Li form-factor minimum at
   |t| ≈ 0.31 GeV² (q² = 8 fm⁻², Li–Sick–Whitney–Yearian 1971), the two
   components have very different |t| dependences — a two-component fit
   in |t| extracts the coherent yield statistically even without
   event-by-event ID. Requires modeling the incoherent slope: BeAGLE
   Mode-W study (plans/02 step 1.5.4) plus the cluster-IA spectra of
   `tagged.py`. Benchmark for what vetoes+fits achieve on a heavy ion:
   the e+Pb coherent-J/ψ study (arXiv:2108.01694, PRD 104:114030)
   rejects 80–99% of incoherent events vs |t| (~2% residual after all
   far-forward vetoes) — light nuclei are *harder* (fewer emitted
   neutrons/photons), which is exactly why the α+d channel needs Z-ID.
3. **Partner-fragment veto.** Vetoes kill the channels with a
   non-beam-rigidity fragment: p → OMD, n → ZDC, ³He → RP main window.
   These channels are therefore *not* the problem (their veto
   inefficiency is a second-order correction, quantifiable with the
   same BeAGLE study). The pure α+d channel is vetoable only when the
   *second* fragment also fluctuates above the pT cut (rare —
   correlated small pT) — worth quantifying, not worth relying on.
4. **Central/B0 activity.** Genuine coherent diffraction has a rapidity
   gap and an empty B0/forward region; incoherent DIS breakup (the
   non-diffractive feed-down) mostly does not. A gap/quietness cut
   suppresses the non-diffractive part of the α+d feed by a large
   factor before the far-forward question even arises.
5. **IR-8 secondary focus (structural fix).** RPs at the z ≈ 44/45.5 m
   secondary focus resolve ~1% rigidity changes and reach pT ≈ 0; the
   eSTARlight coherent-J/ψ study (arXiv:2511.05638) quotes tagging
   efficiency × acceptance of 47% (d), 32% (³He), 29% (⁴He), **17.8%
   (⁷Li)** at top energies — and the proposed Z² Cherenkov lives there.
   The coherent-⁶Li program is one more physics case for the second IR.

### (2) ⁶Li excited states that come back to an intact-looking ⁶Li

- **3.5629 MeV 0⁺, T = 1** (Γ = 8.2 eV; isobaric analog of the ⁶He
  ground state): particle-*stable* in practice — the α+d decay is
  strictly forbidden by angular momentum + parity (0⁺ → 0⁺⊗1⁺ needs
  L = 1, wrong parity; measured Γ_α ≤ 6×10⁻⁷ eV) on top of being
  isospin-forbidden — so it M1-γ-decays to the ground state (~100%,
  TUNL A=6). A coherent excitation delivers a real A = 6, Z = 3 track
  at R ≈ 1.000 — indistinguishable in the RP — plus a de-excitation
  photon: E* = 3.56 MeV boosts to E_lab ~ γE*(1+cosθ*) ≈ 0.1–0.4 GeV
  at mid beam energy, emitted within ~1/γ ≈ 10–20 mrad — i.e. mostly
  in the **B0 EMCal photon acceptance (≥ 50–100 MeV, 5.5–20 mrad,
  arXiv:2409.02811)** with the small-angle tail in the ZDC (0.1–1 GeV
  photon requirement). So this background is *vetoable by design*, and
  doubly suppressed at production: **isoscalar diffractive exchange
  cannot drive T = 0 → T = 1**. Treat as small and modelable.
- **2.186 MeV 3⁺, T = 0** (Γ = 24 keV, ~100% α+d) and the broad
  4.31 MeV 2⁺ and 5.65 MeV 1⁺ (T = 0): particle-unstable → α+d; these
  are exactly background (1), listed here to note that a "quasi-elastic
  to low-lying resonances" component sits right on top of the coherent
  kinematics with only ΔpT ~ 40–90 MeV of breakup smearing (Q = 0.71
  MeV for the 3⁺).
- **5.366 MeV 2⁺, T = 1** (Γ = 541 keV): its α+d branch is
  isospin-forbidden (< 1%); it decays ~35% p+⁵He / ~65% p+n+α (TUNL) —
  i.e. it *always* emits a proton or neutron and is caught by the
  OMD/ZDC vetoes like the other nucleon-emission channels.

The pattern worth stating on any veto-strategy slide: **every T = 1
state is either γ-vetoable (3.56) or nucleon-vetoable (5.37), while the
T = 0 states feed the beam-blind α+d channel** — the veto budget is set
entirely by the T = 0 α+d continuum and resonances.

### (3) QED with an intact nucleus: radiative elastic / Bethe–Heitler / QED Compton

e ⁶Li → e γ ⁶Li keeps the nucleus intact *by construction* and passes
the coherent selection whenever the photon escapes (small-angle ISR/FSR
peaks). It is fully calculable, its t-dependence is the *charge* form
factor squared — so it doubles as an **in-situ calibration candle** for
the RP acceptance and (with tensor polarization) a monitor of the
alignment axis; the charge-quadrupole cos 2φ it induces is bounded by
the anomalously small Q(⁶Li) = −0.0806(6) fm² (the famous
near-cancellation),
i.e. negligible against percent-level gluonic scenarios — but it must be
*subtracted*, not ignored, at per-mille precision. Standard kill: E−p_z
balance of (e′, γ, recoil) and the exclusivity cut on M_X.

### (4) Machine and accidental backgrounds

- **Beam-gas**: ion-beam + residual-gas collisions produce beam-rigidity
  fragments (including α at R = 1) into the RP tail, in accidental
  coincidence with a DIS e′. Measured directly with non-colliding
  bunches and killed by vertex/timing association (AC-LGAD ~30 ps
  timing) — the standard EIC procedure; flagged because the tail
  acceptance is exactly where beam-halo lives.
- **Off-momentum halo / satellite bunches**: populate |R−1| ≲ 10⁻³
  trajectories; sets how close to 10σ the pots can actually sit — feeds
  back into the pT_cut value, i.e. directly into the 13.5% acceptance.
- **Pileup of a breakup fragment from a *different* crossing**:
  suppressed by timing, same machinery as beam-gas.

### (5) Electron-side backgrounds (both channels)

Photoproduction π⁻ faking e′ and pair-symmetric electrons at high y —
standard DIS practice applies (E/p + calorimeter ID, charge-symmetric
subtraction with wrong-sign tracks / positron running). For the
*coherent* selection these are further suppressed by requiring the
(e′, recoil) kinematic consistency (x_P, t balance).

## 6.3 Systematics specific to the cos 2φ extraction

- **φ-acceptance holes** couple the cos 2φ moment to the acceptance
  Fourier modes: demonstrated in step 5.A — the binned-LSQ estimator
  (`cos2phi_fit`, now also `cos2phi_fit_binned`) is unbiased where the
  naive moment fails. For the coherent channel the *recoil-side*
  φ_t acceptance of the RPs is also non-uniform; the measurement
  azimuth φ comes from the electron, but t⃗ and the lepton plane are
  correlated, so a 2-D (φ, φ_t) acceptance map is the safe upgrade
  (Phase-2 closure test, plans/03 step 2.3).
- **Relative luminosity: immune.** The modulation is measured *within a
  single fill category* (transverse tensor fill), self-normalized in φ —
  the 10⁻⁴-level rel-lumi systematics of A_zz/A∥ (step 5.A) do not
  enter. This is a real advantage of the cos 2φ observable and worth
  stating on the plot.
- **Spin-axis control**: dilution ∝ sin²θ_S (a 5° polar error is a
  0.8% dilution — negligible); a *drifting* φ_S rotates the modulation
  phase — bunch-by-bunch φ_S bookkeeping is already in the run-plan
  machinery. Polarimetry scale δP_zz/P_zz ≈ 3% is an overall amplitude
  scale, common with A_zz (plans/04 #5).
- **Vector-polarization leakage**: with unpolarized electrons the
  vector terms enter only through ⟨λ_e⟩ ≠ 0 residuals (cos φ′, not
  cos 2φ′ — orthogonal in a full-φ fit; relevant only combined with
  acceptance holes; covered by the binned fit + bookkeeper).
- **Tensor-observable radiative corrections**: still uncharted
  (plans/04 #10); for the coherent channel the QED radiative tail is
  item 6.2(3) and is calculable — flag both on every plot.
- **Low-Q², high-x nuclear effects**: target-mass/γ² terms are absent
  from the master formula (known fastsim caveat); the sweet spots sit
  at x ≤ 0.06 where γ² = 4M²x²/Q² < 10⁻² even at Q² ≈ 1 GeV², so this
  is a labeling caveat, not a bias, for money plot 5.

## 6.3b Interpretational backgrounds to "exotic glue" (theory, verified)

- **Leading-twist protection**: Δ(x,Q²) is defined in Jaffe–Manohar
  PLB 223:218 ("Nuclear gluonometry"; NPB 312:571 is the b₁ companion —
  attribution nuance for our docstrings): spin-½ nucleons and pions
  cannot transfer two units of helicity, so conventional nuclear
  physics cannot generate Δ *at leading twist*. Power-suppressed
  (target-mass / elastic-intermediate-state / quadrupole) contributions
  must die as powers of 1/Q² — but the sweep found **no published
  quantification** of this low-Q² coherent background for Δ; our Q²
  lever arm across the sweet-spot panels (1.1 → 14 GeV²) is the
  built-in consistency check, and quantifying the power correction is
  publishable new work, not a citation.
- **Non-nucleonic but non-gluonic mechanisms**: Nzar–Hoodbhoy
  (PRD 45:2264, 1992) showed a ΔΔ isobar admixture (spin-3/2
  constituents *can* double-flip) generates a deuteron Δ — an
  interpretational background even at leading twist. For ⁶Li this
  motivates quoting Δ alongside the b₁ program rather than alone.
- **The b₁ cautionary analog**: at x < 0.1, conventional D-state +
  coherent double-scattering/shadowing generate *large* b₁
  (Bora–Jaffe PRD 57:6906; Nikolaev–Schäfer PLB 398; Edelmann–Piller–
  Weise PRC 57:3392) — the reason the HERMES low-x b₁ is not exotic.
  No equivalent calculation exists for the cos 2φ Δ channel; the same
  mechanisms should be checked there before any discovery claim.
- **The ⁶Li null-test advantage**: coherent azimuthal modulations from
  *nuclear deformation* — exactly what Mantysaari–Salazar–Schenke–
  Shen–Zhao (arXiv:2408.13213) predict for coherent VM production off
  polarized deuterons — are bounded for ⁶Li by its anomalously small
  quadrupole moment, Q(⁶Li) = −0.0806(6) fm² (≈ 50× smaller than ⁷Li's
  −4.00(3) fm²). A sizable coherent cos 2φ on ⁶Li therefore cannot be
  a shape effect — the deformation-driven and gluon-driven mechanisms
  that are entangled for the deuteron are cleanly separated by the
  ⁶Li/d comparison. This is a genuine selling point of the ⁶Li beam.

## 6.4 Scenario parameters (to be replaced, with owners)

| parameter | default | band | grounding (2026-08-10 sweep) |
|---|---|---|---|
| coherent fraction f₀ (x→0) | 0.04 | 0.02–0.08 | **no light-A prediction exists** (lightest published is Ca); band brackets the HERA ep diffractive fraction (10–15% of DIS, of which coherent-elastic is part) and heavy-A saturation estimates (20–25%); replace via #18 |
| coherence falloff x_coh | 0.01 | l_c = 0.105 fm/x_P vs R–2R ≈ 2.4–5 fm → bare bound x_P ≲ 0.04; EIC practice x_P < 0.01 | diffractive model, #18 |
| slope B | 50 GeV⁻² | 40–60 | B = ⟨r²⟩/3: matter radius 2.32–2.45 fm → 45–51; charge radius 2.589(39) fm → ≈57; STAR d+Au found the *gluon* distribution wider than charge, so the gluonic B may sit high. Exponential valid below the first FF minimum, \|t\| ≈ 0.31 GeV² (q² = 8 fm⁻²) |
| modulation amp (P_zz = 1) | 0.02 | 10⁻³ (Sather–Schmidt-scale) … few×10⁻² (max Δ_T g = Δg, Kumano–Song-style; flagged overestimate) | scenario until #18; best-bin 5σ floor 0.003 at 100 fb⁻¹/u — the pessimistic 10⁻³ needs bin combination (~×3) and/or IR-8 acceptance |

**New open questions for plans/04**: **#18** a coherent diffractive
model for (tensor-polarized) ⁶Li — engage the small-x/diffraction
theory community (elliptic-gluon / spin-1 GPD angle) — and **#19** RP
AC-LGAD charge discrimination (Z = 3 vs Z = 2 vs Z = 1 by amplitude)
— engage the ePIC far-forward WG. The IR-8 secondary focus remains the
structural fix for everything at R ≈ 1 (same as the ⁶Li α-tag).

## 6.5 Verification ledger (2026-08-10 research sweep, 4 agents, 32 claims)

All external numbers above were checked against primary sources
(TUNL A=6/A=7 evaluations, EIC Yellow Report + ePIC FF documents, HERA/
EIC diffraction papers, the gluonometry theory chain). Corrections
applied in place:

- 3.5629 MeV 0⁺ α-decay is *parity/angular-momentum* forbidden exactly
  (Γ_α ≤ 6×10⁻⁷ eV), not merely isospin-forbidden; Γ = 8.2(2) eV.
- De-excitation-photon veto belongs mainly to the **B0 EMCal**
  (γ ≥ 50–100 MeV at 5.5–20 mrad), with ZDC (≤ 4–5.5 mrad, 0.1–1 GeV
  requirement) taking the small-angle tail.
- Modern radius pair: R_ch(⁶Li) = 2.589(39) fm > R_ch(⁷Li) = 2.444(42)
  fm; matter radius 2.32–2.45 fm; Q(⁶Li) = −0.0806(6) fm²,
  Q(⁷Li) = −4.00(3) fm².
- The 0.45 GeV high-divergence pT cut is *derived* from YR divergence
  tables (≈ 0.41 GeV at 275 GeV), not a documented spec; 0.2 GeV
  (high-acceptance) is documented and is the value our routing uses.
- Coherent-slope literature check: B_d ≈ 35–40 GeV⁻², B_⁴He ≈ 24
  (LEPS: 23.8(1.0)) — larger nucleus ⇒ steeper slope, consistent with
  our B(⁶Li) ≈ 50; light-nucleus coherent/incoherent crossover
  ~0.05–0.2 GeV², i.e. *inside* our tagged window — reinforcing the
  Z-ID/t-fit requirements of §6.2(1).
- arXiv:2108.01694 is e+**Pb** (not e+Au); veto numbers as cited.
- Kumano–Song citations: PRD 101:054011 + PRD 101:**094013**.
- No published study tags an intact nucleus in EIC RPs quantitatively;
  closest are the YR coherent-helium DVCS (pT ≳ 0.2 GeV ⇒ \|t\| ≳ 0.04
  GeV²), ECCE ⁴He DVCS (arXiv:2208.14575), polarized ³He/⁴He DVCS
  (arXiv:2606.11491), and IR-8 eSTARlight (arXiv:2511.05638) — our
  money plot 6 appears to be the first A = 6 projection of this type,
  extending the plans/00 first-publication list.
