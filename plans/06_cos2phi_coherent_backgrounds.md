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
| `evgen/scripts/money_cos2phi.py` → `money_cos2phi_6Li.png` | inclusive gluonometry: φ′ pseudo-data in 4 sweet-spot super-bins (Q² = 1.1→14 GeV² at x ≈ 0.02–0.14) + amplitude vs x with Δ/F₁ scenario curves |
| `evgen/scripts/money_cos2phi_coherent.py` → `money_cos2phi_coherent_6Li.png` | coherent channel: recoil t-spectrum vs RP cuts, tagged yield vs x, φ′ pseudo-data in the best low-x bin |
| `evgen/polligen/coherent.py` | scenario coherent model + intact-recoil routing + breakup veto table (tested) |

Headline numbers (mid energy e10×⁶Li50, L = 100 fb⁻¹/u, P_zz = 0.6,
scenario inputs):

- **Inclusive:** with the unified Δ registry
  (`polli_fastsim/delta_models`, default `moment_A`: sum-rule-
  constrained, ∫xΔdx = −0.012·α_s, dilution 1/3) the sweet-spot
  amplitudes are (0.7–1.2)×10⁻² (table-α_s, run below the table edge) against per-bin δA ≈ (1.5–4.5)×10⁻⁴
  already at the 1-year program (10 fb⁻¹/u; ÷√10 at 10 years) — the
  measurement resolves the moment_A-vs-moment_B ansatz spread, not just
  a null. [Historical toy-model numbers, kept for reference and
  verified by the 2026-08-10 adversarial review: at Δ/F₁ = 10⁻³ flat
  scale, per-bin δA (0.4–1.2)×10⁻⁴ at 100 fb⁻¹/u gave 2.5–2.7σ per
  low/mid-Q² bin, 4.6σ over the four displayed bins, 9.7σ over the
  full map.]
- **Coherent:** the intact ⁶Li has exactly beam rigidity (A/Z = 2,
  R = 1.000), so at IP6 it is visible **only** in the Roman-Pot
  near-beam pT tail: acceptance exp(−B·pT_cut²) = **13.5%**
  — but see the 2026-08-26 measurement below: the ePIC pot aperture is
  tighter than this envelope and opens horizontally, not vertically
  [9–20% for B ∈ 40–60 GeV⁻²] with high-acceptance optics —
  and **4×10⁻⁵ with high-divergence optics: the coherent program
  fixes the optics choice**. Scenario rates give ~10⁸ tagged events at
  100 fb⁻¹/u; the best super-bin (x ≈ 0.002, Q² ≈ 1.3) reaches
  δA ≈ 6×10⁻⁴ at 10 years (1.9×10⁻³ in year 1), i.e. 5σ floors at
  0.3% / 1.0% modulation amplitudes — with the whole luminosity in the
  transverse-tensor fill, as for all FOMs. The modulation itself is now
  *anchored* (second sweep, §6.4b): the deformation term scaled from
  the polarized-deuteron CGC calculation gives ⟨a₂⟩_tag ≈ 0.036
  [0.018–0.059] at P_zz = 0.6 (one-sided −27% rate-weighting
  systematic, §6.4b audit note) with a predicted sign flip vs the
  deuteron; the gluon-transversity term sits at 3×10⁻³–10⁻² — both
  well above the statistical floors, so the measurement becomes a
  two-component a₂(t) decomposition.

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
table from `coherent.veto_table()`, in mass-to-charge ratios against the
beam's (physical nuclear masses since 2026-08-28, not A/Z: the naive
ratio said both fragments of the α + d channel sat at exactly 1.000,
which was also the value the fast simulation's own spectator model
disagreed with):

| breakup channel | threshold | fragment routing |
|---|---|---|
| **α + d** | 1.474 MeV | α: R = 0.998 → RP pT tail only (beam-blind); d: R = 1.005 → same |
| α + p + n | 3.70 MeV | α beam-blind; p: R = 0.503 → OMD; n → ZDC |
| ³He + t | 15.79 MeV | ³He: R = 0.752 → RP main window; t: R = 1.504 → lost (over-rigid) |

The **α+d channel is the killer**: both fragments sit within 0.5% of beam
rigidity, deep inside the ±5% near-beam band and so undispersed, exactly
like the intact ⁶Li — same trajectory, same velocity (so time-of-flight
cannot help). The lowest-lying strength (direct α+d continuum from 1.47
MeV and the 2.186 MeV 3⁺ resonance, Γ ≈ 24 keV, which decays ~100% to α+d
with only Q ≈ 0.7 MeV) shares the recoil pT in proportion to mass, ~2/3 to
the α, and adds a breakup smear k_rel ≈ √(2μQ) ≈ 40 MeV/c — precisely the
events whose α can wander outside the near-beam envelope and be counted as
a coherent recoil. What the first version of this plan then assumed, that
the partner deuteron stays invisible, is the opposite of what the same
k_rel does: it is *transverse* and unboosted, and the deuteron carries
half the α's longitudinal momentum, so it opens twice the angle. Handle 3
below rests on that.

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
3. **Partner-fragment veto — quantified, and stronger than this plan
   assumed.** Vetoes kill the channels with a non-beam-rigidity
   fragment: p → OMD, n → ZDC, ³He → RP main window. These channels are
   therefore *not* the problem (their veto inefficiency is a
   second-order correction, quantifiable with the same BeAGLE study).
   For the pure α+d channel the two fragments were measured together
   for the first time on 2026-08-28 (plans/09 B4,
   `spectator.breakup_lab_kinematics`, `evgen/scripts/nearbeam_two_hit.py`):
   they are back to back in azimuth with θ_d = 1.987 θ_α, so the
   deuteron is *systematically* the wider fragment rather than a rare
   fluctuation, and the veto is not rare at all where it is needed.
   Conditioned on the α landing in the near-beam tail — the only
   configuration in which this channel is a background — the partner
   deuteron is on a pot in **84–85%** of events at the tagging optics of
   Report 1 §6.1. At the published Yellow Report optics it is 4–29%, but
   there the α fakes a coherent tag in only 10⁻³–10⁻⁴ of breakups: the
   optics that creates the background supplies the veto with it, the
   same conditional structure as every other number in this study.
   Two limits: the 84% assumes pot acceptance out to θ = 5 mrad and
   falls to 0.00 / 0.31 / 0.56 (5 × 41 / 10 × 100 / 18 × 275) against a
   single 16 mm module, so it is a station-layout question (plans/09 B1,
   D3); and it is a *veto*, so its inefficiency enters the coherent
   normalisation and must be folded rather than assumed. Rank: below the
   |t|-shape fit only because of those two, and above everything after
   it.
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
  stating on the plot. *Superseded 2026-08-24/25:* the single-fill fit is
  biased by the detector's cos 2φ acceptance harmonic ÷ P_zz; the
  estimator is the spin-state ratio of m = ±1-rich / m = 0-rich bunches
  (acceptance and relative luminosity cancel for a common acceptance,
  bunch-by-bunch alternation required) — reconstruction-chain report §3,
  code review F1.
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
  at x ≤ 0.15; the low-Q² spots have x ≤ 0.06 (γ² < 10⁻² at Q² ≈ 1
  GeV²) and the x ≈ 0.14 spots sit at Q² ≥ 3 GeV² (γ² ≲ 3×10⁻²), so
  this is a labeling caveat, not a bias, for money plot 5.

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
| modulation, deformation term ε_B0 | −0.08 | −0.04 … −0.13 | scaled from the polarized-d CGC calculation (§6.4b); a₂(t) = −(P_zz/4)·ε_B0·B·\|t\|, i.e. ⟨a₂⟩_tag ≈ 0.036 [0.018–0.059] at P_zz = 0.6 |
| modulation, gluon-transversity term (flat, P_zz = 1) | 0.01 | 3×10⁻³ … 10⁻² | lattice φ-meson transversity ÷ NPLQCD nuclear suppression; Kumano–Song few-% only at max Δ_T g (§6.4b); best-bin 5σ floor 0.003 at 100 fb⁻¹/u |

## 6.4b Literature anchors for the modulation (second research sweep, 2026-08-10)

The question "are there estimates like Mäntysaari et al. 2408.13213 we
can use?" has a sharp answer: **that paper is the only small-x
calculation of polarization-dependent azimuthal coefficients in coherent
diffraction on any polarized nucleus** — all 8 forward citations were
checked; none extends it to A > 2, and the group's own follow-ups are
unpolarized. So we anchored by scaling, and digitized their figures from
the arXiv source (now in `coherent.MANTYSAARI_A2_DEUTERON`):

- **The deuteron anchor** (PLB 858:139053; coherent J/ψ, x_P = 1.7×10⁻³,
  transverse polarization): a₂(m=0) = +0.08/+0.15/+0.30/+0.43 at
  \|t\| = 0.05/0.1/0.2/0.3 GeV², with a₂(m=±1) ≈ −a₂(m=0)/2 — a factor-2
  cross-section difference between m-states at \|t\| = 0.1 for full
  alignment. Mechanism: the D-wave-deformed (6% P_D) gluon density
  modulates the diffractive slope, a₂ ≈ (ΔB_m/2)·\|t\| with
  ΔB₀/B ≈ 0.21. Ensemble relation: a₂ ∝ P_zz (exact cancellation for an
  unpolarized fill). Two more transferable facts: only even harmonics
  appear; and the incoherent (breakup) cross section is
  **m-state-blind**, so incoherent contamination dilutes but does not
  fake the modulation — one more reason the tagged-coherent selection is
  the right one.
- **⁶Li scaling, two dials that disagree usefully**: relative quadrupole
  deformation Q/⟨r²⟩ is 5× smaller than the deuteron's (and opposite in
  sign); the α–d asymptotic D/S ratio is itself
  DISPUTED (2026-08-11 reference audit): the (⁶Li,d) transfer
  measurement of Veal et al. PRC 60:064003 gives η = +0.0003(9) —
  consistent with zero, mirroring the Q cancellation — while older
  sub-Coulomb analyses favored η ≈ −0.01 (vs η_d = +0.0256). Hence
  ε_B0(⁶Li) ∈ −(0.04–0.13) spans the Q/⟨r²⟩ dial up to a
  no-cancellation deuteron-like value, with near-zero gluonic
  deformation an open possibility covered by the quoted 5σ floors;
  the **predicted sign flip relative to the deuteron** (from Q < 0)
  remains a testable statement. The Wiringa–Schiavilla VMC result (even exact
  few-body theory misses the ⁶Li quadrupole cancellation by ~3×) is why
  the band stays wide, and why the *gluonic* deformation need not cancel
  like the charge one — measuring it is the physics.
- **Existence proofs at our exact \|t\|**: elastic e-d T₂₀ is measured
  at −(0.1–0.3) for Q² = 0.04–0.3 GeV² (monopole–quadrupole
  interference); HERMES measured A_zz(d) ≈ −1×10⁻² at x ≈ 0.01, the
  low-x rise being coherent-double-scattering in origin
  (Nikolaev–Schäfer predicted the ~1% shadowing-driven tensor scale) —
  the same diffractive physics as our channel. The only
  coherent-exclusive tensor bound is HERMES A_Lzz < 0.2 — everything
  below is unexplored.
- **Gluonic term**: lattice gives near-maximal gluon transversity in the
  compact φ meson but a ~10× smaller transversity/unpolarized ratio in
  the deuteron (NPLQCD) → flat scenario 3×10⁻³–10⁻² at P_zz = 1,
  corroborated by Kumano–Song pd Drell-Yan (few % only for Δ_T g = Δg).
- **The analysis strategy this buys**: the two mechanisms separate by
  *shape* — deformation is t-linear, sign-flipped vs the deuteron, and
  should weaken toward small x_P (JIMWLK washes out deformation,
  arXiv:2411.14934); gluon transversity is flat-ish in t and survives at
  small x. A two-component fit a₂(t) = c_def·\|t\| + a_g over the tagged
  window, repeated in x_P bins, is the measurement — and both components
  sit well above the best-bin statistical floors (δa₂ = P_zz·δA ≈
  1.1×10⁻³ / 3.6×10⁻⁴ at the 1-/10-year programs). Audit note
  (2026-08-10): the simple population average behind a₂ ∝ P_zz assumes
  equal m-state rates — exact only as |t| → 0; the anchor's own m-state
  rate differences imply a one-sided ≈ −27% rate-weighting systematic
  on ⟨a₂⟩_tag (0.026 vs 0.036 at the defaults,
  `coherent.RATE_WEIGHT_SYST`), covered by the ε_B0 band edge.
  Validity limit: quote for \|t\| < 0.2 GeV²; the
  C₀ zero at ≈ 0.31 GeV² locally enhances and sign-flips a₂ (the
  deuteron analog of 2408.13213's dip region) — flagged, not modeled.
- **Realism checks**: solid-target programs reach P_zz ≈ 0.5–0.7 (JLab
  Hall D tensor-deuteron, arXiv:2504.21177); ABS-type Li sources have
  demonstrated ~95% of maximal tensor polarization, so the in-ring 0.6
  default stands. The eSTARlight IR-8 study's intact-recoil efficiencies
  (d 47%, ³He 32%, ⁴He 29%, ⁷Li 17.8%) interpolate to ≈ 20% for ⁶Li at
  IR-8 — versus our 13.5% RP pT-tail estimate at IP6 (constant 0.20 GeV cut;
  and versus the **measured** aperture of 2026-08-26, |θ_x| ≳ 2.0 / 1.35 /
  1.03 mrad in the 5×41 / 10×100 / 18×275 optics, which at the γ-matched
  40.8 / 99.5 / 137.5 GeV/u gives a tagged fraction of
  **9.8×10⁻⁷ / 7.7×10⁻¹⁶ / 1.9×10⁻¹⁷** at B = 50 GeV⁻² and inverts
  the cutout aspect this section assumes — `evgen/scripts/nearbeam_aperture_scan.py`,
  `tools/fullsim/README.md`, plans/04 #20.  The measured edge is
  0.91× / 0.75× / 1.12× the horizontal half-width of the 10σ Yellow
  Report high-acceptance envelope, so at the two lower configurations the
  ENVELOPE binds and the tagged fraction is 7.2×10⁻⁸ / 6.2×10⁻²⁷ /
  3.9×10⁻¹⁴ instead; the 1.4×10⁻² / 2.7×10⁻³ this bullet carried until
  2026-08-28 was the measured aperture at the retired rigidity-scaled
  momenta, and the "envelope is never binding" pricing of plans/08 §8.4
  went with it.  Nothing published survives at either aperture: the
  channel exists only under the lithium tagging optics of Report 1 §6.1,
  which recovers 0.42 / 0.32 / 0.33 at 1/7.1 / 1/13.3 / 1/9.5 of the
  luminosity), an independent cross-check of the tagging scale.

**Convention check (2026-08-24, full text of arXiv:2408.13213).** Its
Eq. (9) expands d²σ/(dΦ d|t|) = dσ/d|t| · (1/2π)[1 + 2 Σ_n a_n(|t|) e^{inΦ}]
with Φ "the angle between the produced vector meson and the polarization
direction of the deuteron" (polarizations defined in the γ*d c.m.
frame), i.e. the RECOIL azimuth. The cos 2Φ modulation coefficient of a
pure state is therefore 2a₂(m), and the ensemble coefficient is
c₂ = 2·a₂(t; P_zz) = −(P_zz/2)·ε_B0·B·|t| — twice the ⟨a₂⟩_tag quoted
above and injected by money plot 6 (which is thus conservative by a
factor two: 0.072 instead of 0.036 at P_zz = 0.6, ⟨|t|⟩ = 0.06). The
reconstructed-level plot 6R (`money_cos2phi_coherent_reco.py`) injects
c₂ (`coherent.cos2phi_coefficient_deformation`); the numbers in this
section and in the projection report keep the a₂ convention until the
WP5 curves are redone. The modulation lives in the recoil azimuth
β = φ_t − φ_S; the gluon-transversity term lives in the electron azimuth
α = φ_e − φ_S (reconstruction-chain report §4).

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
