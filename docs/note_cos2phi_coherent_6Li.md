# The cos 2φ Modulation on Tensor-Polarized ⁶Li at the EIC: Projections, the Coherent Intact-Recoil Channel, and Anchored Amplitudes

*Working note, 2026-08-10 — polarized ⁶,⁷Li @ EIC simulation program
(companions: `plans/06` for the running background budget,
`evgen/README.md` for the code). Prepared for circulation ahead of the
INT program "Towards Realizing the Program with Polarized Ion Beams at
EIC", March 22 – April 2, 2027. All rates use TOY/scenario structure
functions and a scenario coherent model — bands, not predictions; every
external number below was verified against primary sources in two
research sweeps (2026-08-10), and the two literature-gap claims
(§4, §7) were checked against complete forward-citation lists.*

---

## 1. The observable and the two channels

For an unpolarized electron beam on a **transversely tensor-polarized**
spin-1 target, the inclusive cross section carries a cos 2φ modulation
proportional to the double-helicity-flip structure function Δ(x,Q²)
[1,2]:

    dσ/(dx dy dφ) = (2yα²/Q²) [ F₁ + (2/3)a_m b₁
                     + (1−y)/(xy²) (F₂ + (2/3)a_m b₂)
                     − (1−y)/y² · c_m sin²θ_m · Δ(x,Q²) cos 2φ ]

Photon helicity flips by two units, so no quark — and no bound spin-½
nucleon — can contribute at leading twist: Δ is purely gluonic, "exotic
glue" [1]. It has never been measured for any target; the only prior
experimental document is a JLab letter of intent [3]. Magnitude
guidance: the Sather–Schmidt bag-model moment is −0.012 α_s [4], i.e.
per-mille scale, while maximal-Δ_T g scenarios reach a few percent [5];
we use Δ/F₁ ∈ {10⁻³ … 10⁻²} as the scenario band.

Two measurement channels, with different detection demands:

- **Inclusive (breakup) DIS** — only the scattered electron is needed;
  φ is the lepton-plane azimuth relative to the transverse alignment
  axis. This is the classic Jaffe–Manohar measurement [1]; the full
  spin-1 DIS formalism including b₁ is [2].
- **Coherent diffractive DIS**, e ⁶Li → e′ X ⁶Li(g.s.) — the nucleus
  stays in its 1⁺ ground state and is detected in the far-forward
  spectrometers. The cos 2φ modulation of the *coherent* yield probes
  the gluonic quadrupole / double-helicity-flip structure of the intact
  nucleus — the nucleus-intact analog of the same physics, connected at
  small x to the "elliptic" gluon distribution program [6]. No
  published calculation of this observable exists for any A > 2 nucleus
  (§4), and no projection of an intact-A = 6 tag exists at the EIC
  (§5) — both are firsts available to this program.

## 2. Money plot 5 — inclusive gluonometry as projected data points

`evgen/scripts/money_cos2phi.py` (→ `money_cos2phi_6Li.png`).
Machinery: polligen pseudo-experiments; per-φ-bin Poisson statistics
drawn from exact expected yields, so any luminosity is exact for the
binned estimator; the whole luminosity is assigned to the transverse
tensor fill; scattered-electron acceptance cuts applied.

Setup: mid energy e(10 GeV) × ⁶Li(50 GeV/u), P_zz = 0.60 (in-ring
placeholder); projections quoted separately for the 1-year (10 fb⁻¹/u)
and 10-year (100 fb⁻¹/u) EIC programs. The Δ(x,Q²) input comes from the
**unified model registry** `polli_fastsim/delta_models.py` — the single
home for all Δ models and constraints, every consumer switching by
name. Default: `moment_A`, the sum-rule-constrained ansatz of the
merged money_delta suite, Δ = A·α_s(Q²)·F₁·x^a(1−x)^b with A solved
from the bag-model moment ∫x Δ dx = −0.012·α_s (Sather–Schmidt), times
the ⁶Li per-nucleon dilution 1/3 (Cloët convention; the whole-nucleus
P_zz is a separate factor — plans/04 #6). `moment_B` (no F₁ factor,
conservative; factor ~2–8 below A) and the legacy `toy` shape are one
flag away. Sweet-spot super-bins per Q² band: Q² = 1.1 → 14 GeV² at
x ≈ 0.02–0.14.

Results (moment_A, mid_x variant, solved A = −0.294 at ⟨Q²⟩ = 3.9 GeV²,
table-α_s — the money_delta production convention; the bag-target moment
applied to the per-nucleon Δ of ⁶Li is itself a scenario choice
inherited from that line, made explicit in the 2026-08-10 audit):

| super-bin (x, Q²) | amplitude | δA, 1 yr | δA, 10 yr |
|---|---|---|---|
| 0.056, 1.14 GeV² | +1.07×10⁻² | 1.8×10⁻⁴ | 5.7×10⁻⁵ |
| 0.022, 1.14 GeV² | +0.65×10⁻² | 1.5×10⁻⁴ | 4.6×10⁻⁵ |
| 0.14, 3.1 GeV²   | +1.24×10⁻² | 2.6×10⁻⁴ | 8.1×10⁻⁵ |
| 0.14, 14.3 GeV²  | +0.96×10⁻² | 4.5×10⁻⁴ | 1.4×10⁻⁴ |

(α_s from the CT18NLO table, continued running below the table edge —
2026-08-11 audit fix; a frozen edge value had suppressed the two
low-Q² amplitudes by ~10%.)

Under the moment constraint the measurement is no longer a null test:
every sweet-spot bin is a 21–44σ effect in year 1, and the
moment_A-vs-moment_B spread — the dominant ansatz systematic (their
notes: factor 2–8 in reach) — is itself resolved by the amplitude-vs-x
shape. At 10-year precision the statistical errors (≈ 5×10⁻⁵) sit
well below the polarimetry-scale systematic (3% of A ≈ 10⁻² is
3×10⁻⁴): the 10-year measurement is systematics-limited. The
extraction here divides by the true P_zz; polarimetry smearing exists
in the run-plan bookkeeper but is not exercised in these plots.
(Historical flat-toy numbers at Δ/F₁ = 10⁻³, verified by the
adversarial review: 2.5–2.7σ per bin, 9.7σ full-map at 100 fb⁻¹/u.)
The Q² lever arm across the panels remains the built-in consistency
check against power-suppressed (non-partonic) contributions (§6).

**Money plot 7** (`money_delta_extraction.py` →
`money_delta_extracted_6Li.png`) unfolds the same pseudo-measurements
to the structure function itself, Δ̂ = Â·y²D_φ/(1−y) per bin (model
bin-centering; errors δΔ = δA·|∂Δ/∂A|), presented as xΔ(x,Q²) at the
three sweet-spot Q² slices with independent 1-/10-year draws — the
world-data-style figure a reader expects for a structure-function
program, with the area under each curve equal to the Sather–Schmidt
moment (÷3 dilution). Best bins reach δΔ ≈ (0.7–1.4)×10⁻³ in year 1
on Δ ≈ −0.03…−0.09 per nucleon.

## 3. Money plot 6 — the coherent intact-⁶Li channel

`evgen/scripts/money_cos2phi_coherent.py`
(→ `money_cos2phi_coherent_6Li.png`), model in `polligen/coherent.py`.

**The tagging problem.** The intact ⁶Li recoil has A/Z = 2 — exactly
the beam rigidity (R = 1.000; the diffractive momentum loss x_P < 0.01
changes rigidity by < 1%). At IP6 it is therefore visible **only** in
the Roman-Pot near-beam pT tail above the 10σ beam-envelope cut [7,8]:

- acceptance = exp(−B·pT_cut²) = **13.5%** [9–20% for B ∈ 40–60 GeV⁻²]
  with high-acceptance optics (pT > 0.20 GeV, documented [7,9]);
- **4×10⁻⁵** with high-divergence optics (pT ≳ 0.45 GeV, derived from
  Yellow-Report divergence tables, not a documented spec) —
  **the coherent program fixes the optics choice**;
- the IR-8 secondary focus reaches pT ≈ 0; the eSTARlight study [10]
  quotes intact-recoil tagging efficiency × acceptance of 47% (d),
  32% (³He), 29% (⁴He), 17.8% (⁷Li) at top energies — interpolating to
  ≈ 20% for ⁶Li, an independent cross-check of the tagging scale.

**Rates** (scenario coherent fraction f₀ = 0.04 of DIS at x → 0, ×2÷2
band; coherence falloff x_coh = 0.01; slope B = 50 GeV⁻² from the
matter radius, exponential valid below the ⁶Li form-factor minimum at
|t| ≈ 0.31 GeV² [11,12]): ~1.1×10⁷ RP-tagged events in year 1
(10 fb⁻¹/u), ~1.1×10⁸ over the 10-year program. Best super-bin
(x ≈ 0.002, Q² ≈ 1.3 GeV²): δA = 1.9×10⁻³ (1 yr) / 6×10⁻⁴ (10 yr),
i.e. 5σ floors at 1.0% / 0.3% modulation amplitudes — both below the
deformation-anchored signal of §4.

## 4. Anchoring the modulation amplitude

**The literature situation** (verified 2026-08-10 against the complete
forward-citation list): the polarized-deuteron CGC calculation of
Mäntysaari, Salazar, Schenke, Shen, and Zhao [13] is the **only**
published calculation of polarization-dependent azimuthal coefficients
in coherent diffraction on any polarized nucleus. None of its citing
papers extends it to A > 2; the group's own follow-ups are unpolarized
[14] or concern deformation evolution [15]. A ⁶Li projection is
therefore an extrapolation by construction — ours is anchored as
follows.

**The deuteron anchor** [13] (coherent J/ψ photoproduction,
x_P = 1.7×10⁻³, transverse polarization; figures digitized from the
arXiv source, tabulated in `coherent.MANTYSAARI_A2_DEUTERON`): with
d²σ/dΦd|t| ∝ 1 + 2Σ a_n cos nΦ, Φ the angle between the recoil
t-vector and the transverse polarization axis,

| \|t\| [GeV²] | a₂(m = 0) | a₂(m = ±1) |
|---|---|---|
| 0.05 | +0.08 | −0.04 |
| 0.10 | +0.15 | −0.08 |
| 0.20 | +0.30 | −0.17 |
| 0.30 | +0.43 | −0.28 |

Only even harmonics appear; a₂(0) ≈ −2a₂(±1); at |t| = 0.1 GeV² the
m = ±1 coherent cross section is ~2× the m = 0 one. The mechanism is
geometric: the D-wave-deformed (P_D ≈ 6%) gluon density modulates the
diffractive slope, a₂(m) ≈ (ΔB_m/2)·|t| with ΔB₀/B ≈ 0.21. For a fill
with populations p_m the ensemble coefficient is exactly ∝ P_zz:

    a₂(t; P_zz) = −(P_zz/4) · ε_B0 · B · |t|,   ε_B0 ≡ ΔB₀/B.

Two more facts transfer: the incoherent (breakup) cross section is
**m-state-blind**, so incoherent contamination dilutes but cannot fake
the modulation; and near the diffractive dip a₂ is locally enhanced and
sign-flipping (the deuteron analog of our |t| ≈ 0.31 GeV² caveat).

**Scaling d → ⁶Li**, two dials that usefully disagree:

1. relative quadrupole deformation Q/⟨r²⟩: deuteron 0.286 fm²/(2.13 fm)²
   ≈ 0.063 vs ⁶Li 0.0806 fm²/(2.59 fm)² ≈ 0.012 — a factor ~5 smaller,
   **opposite sign** (Q(⁶Li) < 0);
2. the α–d asymptotic D/S ratio — itself DISPUTED: the (⁶Li,d)
   transfer measurement of Veal et al. [16] gives η = +0.0003(9),
   consistent with zero and mirroring the Q cancellation, while older
   sub-Coulomb analyses favored η ≈ −0.01 (vs η_d = +0.0256). The
   geometric dial therefore spans zero to ~0.5× the deuteron's.

Hence **ε_B0(⁶Li) ∈ −(0.04 … 0.13), default −0.08** (near-zero
deformation remains open if the gluonic sector shares the charge-sector
cancellations — detecting or excluding it is covered by the quoted 5σ
floors), giving

    ⟨a₂⟩_tag = a₂(⟨|t|⟩_tag; P_zz) ≈ 0.036  [0.018 – 0.059]

at P_zz = 0.6 with ⟨|t|⟩_tag = pT_cut² + 1/B = 0.06 GeV² (equal-rate
population average; the anchor's m-state rate differences imply a
one-sided ≈ −27% rate-weighting systematic, `coherent.RATE_WEIGHT_SYST`,
covered by the band edge — 2026-08-10 audit), and a
**predicted sign flip relative to the deuteron** — itself a testable
wave-function statement. The width of the band is physics, not
sloppiness: even exact few-body theory (VMC, AV18+UIX) misses the ⁶Li
quadrupole cancellation by a factor ~3 [17], and the *gluonic*
deformation need not cancel like the charge one — measuring it is the
point. A caveat inherited from [15]: JIMWLK evolution drives deformed
densities toward sphericity at small x, so the deformation term may
weaken toward the smallest x_P.

**Existence proofs at exactly our |t|**: the measured elastic e-d
tensor analyzing power T₂₀ = −(0.1–0.3) over Q² = 0.04–0.3 GeV²
(monopole–quadrupole interference) [18,19]; and HERMES measured
A_zz(d) ≈ −1×10⁻² at x ≈ 0.01 [20], whose low-x rise is
coherent-double-scattering in origin — the ~1% shadowing-driven tensor
scale predicted by Nikolaev–Schäfer [21] (cf. [22,23]). The only bound
on a coherent-*exclusive* tensor asymmetry is HERMES A_Lzz < 0.2 [24]:
everything below is unexplored territory.

**The gluon-transversity (exotic glue) term**, kept as a flat scenario:
lattice QCD finds a near-maximal gluon transversity in the compact φ
meson [25] (with off-forward matrix elements in [26]) but a ~10×
smaller transversity/unpolarized ratio in the deuteron [27]; the
Kumano–Song pd Drell-Yan asymmetry reaches a few percent only for the
maximal assumption Δ_T g = Δg [5,28]. Band: **3×10⁻³ – 10⁻² at
P_zz = 1**.

**The resulting measurement strategy**: the two mechanisms separate by
shape —

| | deformation term | gluon-transversity term |
|---|---|---|
| t-dependence | linear, a₂ ∝ \|t\| | ~flat |
| sign vs deuteron | flipped (Q < 0) | independent |
| x_P dependence | weakens toward small x_P [15] | survives (leading-twist gluonic) |

A two-component fit a₂(t) = c_def·|t| + a_g over the tagged window,
repeated in x_P bins, is the measurement; both components sit ≥ 5σ
above the 3×10⁻⁴-level per-bin statistical floor at 100 fb⁻¹/u. Quote
for |t| < 0.2 GeV² only; the C₀ zero at ≈ 0.31 GeV² [11,12] is flagged,
not modeled. The formal spin-1 language (S_LL/S_TT structure functions;
9 gluon transversity GPDs at leading twist, with sum rules tying their
moments to the tensor form factors) is in [29–31]; deuteron DVCS
formalism in [32,33]; a non-exotic ΔΔ-isobar mechanism that can mimic
Δ at leading twist is estimated in [34].

**The ⁶Li null-test argument**: deformation-driven azimuthal
modulations — exactly what [13] predicts for the deuteron — are bounded
on ⁶Li by its anomalously small quadrupole moment,
Q(⁶Li) = −0.0806(6) fm² ≈ Q(⁷Li)/50 [35]. A sizable flat-in-t coherent
cos 2φ on ⁶Li cannot be a shape effect; the deformation and gluon
mechanisms that are entangled for the deuteron separate cleanly in the
⁶Li/d comparison. Tensor-polarized ⁶Li beams are explicitly motivated in the EPIOS white
paper [50] ("the best chance for discovery may come from larger
nuclei… spin-1 ⁶Li"); polarized internal-target options at the EIC are
discussed in [36]; funded tensor-target programs at JLab
(P_zz = 0.5–0.7) show the observable class is considered
experimentally realistic [37,38].

## 5. Background budget for the coherent tag (summary; running version in plans/06)

Background = anything putting a beam-rigidity track into the RP pT tail
with a quiet central detector. Ranked:

1. **α+d breakup — the killer** (threshold 1.474 MeV; the 2.186 MeV 3⁺,
   Γ = 24 keV, decays ~100% to α+d [39]): both fragments have
   A/Z = 2 → R = 1.000, same trajectory *and* same velocity as an
   intact ⁶Li — only dE/dx (∝ Z²: 9/4/1 for Li/α/d) can separate.
   **No EIC document addresses this A/Z = 2 degeneracy**; the RP
   AC-LGAD chain records pulse amplitude but only for position
   interpolation, and the only documented Z-ID concept is a Z²
   Cherenkov behind the IR-8 secondary-focus pots [8,40]. Open question
   #19 (ePIC FF WG). Fallback: two-component |t|-shape fit (coherent
   e^{−50|t|} vs the ~10× flatter incoherent slope; the light-nucleus
   coherent/incoherent crossover sits at ~0.05–0.2 GeV², *inside* the
   tagged window [41,42]). Heavy-ion benchmark for veto power: the
   e+Pb coherent-J/ψ study rejects 80–99% of incoherent events vs |t|,
   ~2% residual [43] — light nuclei are harder (fewer emitted
   neutrons/photons).
2. **⁶Li* states** [39]: the 3.5629 MeV 0⁺ T=1 (Γ = 8.2 eV) is
   particle-stable — α+d is strictly parity/angular-momentum forbidden
   (Γ_α ≤ 6×10⁻⁷ eV) — and γ-decays to an intact-looking ⁶Li; but
   isoscalar diffractive exchange cannot drive T=0 → T=1, and the
   boosted 0.1–0.4 GeV photon lands mostly in the B0 EMCal acceptance
   (γ ≥ 50–100 MeV, 5.5–20 mrad [9]) with the small-angle tail in the
   ZDC. The 5.366 MeV 2⁺ T=1 always emits a nucleon (35% p+⁵He,
   65% p+n+α; α+d < 1%, isospin-forbidden [39]) → OMD/ZDC-vetoable.
   **Pattern: every T=1 state is γ- or nucleon-vetoable; the T=0 states
   (2.19, 4.31, 5.65 MeV) feed the beam-blind α+d channel** — the veto
   budget is set entirely by them.
3. **Nucleon-emission breakup** (α+p+n at 3.70 MeV, ³He+t at
   15.79 MeV [39]): vetoable — p → OMD (R = 0.50), n → ZDC,
   ³He → RP main window (R = 0.75); t is lost (R = 1.50, over-rigid).
4. **Bethe–Heitler / QED Compton** e ⁶Li → e γ ⁶Li: intact by
   construction; fully calculable (t-dependence = charge FF² [11,12]);
   subtract, and reuse as the RP-acceptance calibration candle. The
   charge-quadrupole cos 2φ it induces is < 10⁻³ for |t| < 0.15 GeV²
   (the ⁶Li C2 form factor is unobservably small below q ≈ 3 fm⁻¹
   [11,12,17]).
5. **Machine backgrounds**: beam-gas/halo fragments at R ≈ 1 in
   accidental coincidence — non-colliding bunches, ~30 ps AC-LGAD
   timing, vertex association; halo determines how close the pots sit,
   i.e. feeds the pT cut itself.
6. **e′-side**: photoproduction π and pair-symmetric electrons at high
   y — standard DIS practice plus (e′, recoil) kinematic consistency.

Systematics of the extraction: the modulation is measured within a
single fill and self-normalized in φ — relative-luminosity systematics
drop out; φ-acceptance holes are handled by the binned-LSQ estimator
(demonstrated unbiased with holey acceptance); polarimetry scale
δP_zz/P_zz ≈ 3% is common with A_zz; tensor-observable radiative
corrections remain uncharted and are flagged on every plot.

## 6. Scenario parameters

| parameter | default | band | grounding |
|---|---|---|---|
| coherent fraction f₀ (x→0) | 0.04 | 0.02–0.08 | no light-A prediction exists (lightest published: Ca); brackets HERA ep diffraction (10–15% of DIS [44]) and heavy-A saturation estimates (20–25% [45,46]) |
| coherence falloff x_coh | 0.01 | l_c = 0.105 fm/x_P vs 2R ≈ 5 fm | EIC practice x_P < 0.01 |
| slope B | 50 GeV⁻² | 40–60 | B = ⟨r²⟩/3: matter radius 2.32–2.45 fm → 45–51; charge radius 2.589(39) fm [47] → ≈57; gluonic B may sit high [42] |
| deformation ε_B0 | −0.08 | −0.04 … −0.13 | §4 scaling of [13]; sign flip vs d predicted |
| gluonic amp (flat, P_zz=1) | 0.01 | 3×10⁻³–10⁻² | [25,27,5] |

## 7. Firsts available to this program (verified gaps)

1. First Δ(x,Q²) sensitivity projection for any target (extends the
   plans/00 gap list; the LOI [3] has no projection).
2. First tagged intact-recoil coherent projection for any A = 6 — no
   published study tags an intact nucleus in EIC Roman Pots
   quantitatively; closest are coherent-helium DVCS in the Yellow
   Report [7], ECCE ⁴He DVCS [48], polarized ³He/⁴He DVCS [49], and the
   IR-8 eSTARlight study [10].
3. First extension of the polarized-nucleus coherent-diffraction
   program [13] beyond the deuteron (verified against all forward
   citations, 2026-08-10).
4. The ⁶Li deformation null test (§4) — unavailable to the deuteron.

## References

[1] R. L. Jaffe, A. Manohar, *Nuclear gluonometry*, Phys. Lett. B 223
(1989) 218.
[2] P. Hoodbhoy, R. L. Jaffe, A. Manohar, *Novel effects in deep
inelastic scattering from spin-one hadrons*, Nucl. Phys. B 312 (1989)
571.
[3] J. Maxwell et al., *Search for Exotic Gluonic States in the
Nucleus*, JLab LOI12-16-006, arXiv:1803.11206.
[4] E. Sather, C. Schmidt, *Size and scaling of the double-helicity-flip
hadronic structure function*, Phys. Rev. D 42 (1990) 1424.
[5] S. Kumano, Q.-T. Song, *Gluon transversity in polarized
proton-deuteron Drell-Yan process*, Phys. Rev. D 101 (2020) 054011,
arXiv:1910.12523.
[6] Y. Hatta, B.-W. Xiao, F. Yuan, *Probing the small-x gluon tomography
in correlated hard diffractive dijet production in DIS*, Phys. Rev.
Lett. 116 (2016) 202301.
[7] R. Abdul Khalek et al., *Science requirements and detector concepts
for the Electron-Ion Collider: EIC Yellow Report*, Nucl. Phys. A 1026
(2022) 122447, arXiv:2103.05419.
[8] Second EIC detector / IR-8 far-forward studies: arXiv:2211.15746
(§VIII, secondary focus and Z² Cherenkov concept); BNL LDRD 23-050
closeout, arXiv:2602.04636.
[9] M. Pitt, *Physics perspectives with the ePIC far-forward and
far-backward detectors*, PoS DIS2024 (2025) 259, arXiv:2409.02811;
E.-C. Aschenauer et al., *Study of deeply virtual Compton scattering at
the future Electron-Ion Collider*, Phys. Rev. D 112 (2025) 036010,
arXiv:2503.05908; M. Boër et al., 3D-imaging white paper,
arXiv:2512.15064 (current RP station positions).
[10] W. Chang, E.-C. Aschenauer, A. Jentsch, A. Kumar, Z. Tu, Z. Yin,
*Opportunities for imaging light nuclei with a second interaction
region at the Electron-Ion Collider*, Phys. Rev. D 113 (2026) 032018,
arXiv:2511.05638.
[11] L. R. Suelzle, M. R. Yearian, H. Crannell, *Elastic electron
scattering from ⁶Li and ⁷Li*, Phys. Rev. 162 (1967) 992.
[12] G. C. Li, I. Sick, R. R. Whitney, M. R. Yearian, *High-energy
electron scattering from ⁶Li*, Nucl. Phys. A 162 (1971) 583.
[13] H. Mäntysaari, F. Salazar, B. Schenke, C. Shen, W. Zhao, *Spatial
imaging of polarized deuterons at the Electron-Ion Collider*, Phys.
Lett. B 858 (2024) 139053, arXiv:2408.13213.
[14] H. Mäntysaari, H. Roch, B. Schenke, C. Shen, W. Zhao, *Nuclear
structure and saturation effects from diffractive vector meson
production*, Phys. Rev. D 114 (2026) 014068, arXiv:2605.00454
(same framework, unpolarized light nuclei).
[15] H. Mäntysaari, P. Singh, *Energy dependence of the deformed
nuclear structure at small-x*, Eur. Phys. J. C 85 (2025) 1449,
arXiv:2411.14934.
[16] K. D. Veal, C. R. Brune, W. H. Geist et al., *Determination of the
asymptotic D- to S-state ratio for ⁶Li via (⁶Li,d) transfer reactions*,
Phys. Rev. C 60 (1999) 064003 (η = +0.0003(9), consistent with zero).
[17] R. B. Wiringa, R. Schiavilla, *Microscopic calculation of ⁶Li
elastic and transition form factors*, Phys. Rev. Lett. 81 (1998) 4317,
nucl-th/9807037.
[18] M. Garçon, J. W. Van Orden, *The deuteron: structure and form
factors*, Adv. Nucl. Phys. 26 (2001) 293, nucl-th/0102049.
[19] Elastic e-d T₂₀ measurements: M. Bouwhuis et al. (NIKHEF),
Phys. Rev. Lett. 82 (1999) 3755, nucl-ex/9810004; D. Abbott et al.
(JLab t20 Collaboration), Phys. Rev. Lett. 84 (2000) 5053,
nucl-ex/0001006; D. M. Nikolenko et al. (VEPP-3), Phys. Rev. Lett. 90
(2003) 072501.
[20] A. Airapetian et al. (HERMES), *First measurement of the tensor
structure function b₁ of the deuteron*, Phys. Rev. Lett. 95 (2005)
242001, hep-ex/0506018.
[21] N. N. Nikolaev, W. Schäfer, *Nonvanishing tensor polarization of
sea quarks in polarized deuterons*, Phys. Lett. B 398 (1997) 245,
hep-ph/9611460 (and erratum, Phys. Lett. B 407 (1997)).
[22] K. Bora, R. L. Jaffe, *The double scattering contribution to
b₁(x,Q²) in the deuteron*, Phys. Rev. D 57 (1998) 6906, hep-ph/9711323.
[23] J. Edelmann, G. Piller, W. Weise, *Deuteron spin structure
functions at small Bjorken x*, Phys. Rev. C 57 (1998) 3392,
hep-ph/9709455.
[24] A. Airapetian et al. (HERMES), *Measurement of azimuthal
asymmetries associated with deeply virtual Compton scattering on a
longitudinally polarized deuterium target*, Nucl. Phys. B 842 (2011)
265, arXiv:1008.3996 (A_Lzz compatible with zero; tensor contribution
to coherent scattering small at low −t).
[25] W. Detmold, P. E. Shanahan, *Gluonic transversity from lattice
QCD*, Phys. Rev. D 94 (2016) 014507, arXiv:1606.04505.
[26] W. Detmold, D. Pefkou, P. E. Shanahan, *Off-forward gluonic
structure of vector mesons*, Phys. Rev. D 95 (2017) 114515,
arXiv:1703.08220.
[27] F. Winter, W. Detmold, A. S. Gambhir, K. Orginos, M. J. Savage,
P. E. Shanahan, M. L. Wagman (NPLQCD), *First lattice QCD study of the
gluonic structure of light nuclei*, Phys. Rev. D 96 (2017) 094512,
arXiv:1709.00395.
[28] S. Kumano, Q.-T. Song, *Deuteron polarizations in the
proton-deuteron Drell-Yan process for finding the gluon transversity*,
Phys. Rev. D 101 (2020) 094013, arXiv:2003.06623.
[29] J. Zhao, A. Bacchetta, S. Kumano, T. Liu, Y.-J. Zhou,
*Semi-inclusive deep inelastic scattering off a tensor-polarized spin-1
target*, JHEP 12 (2025) 067, arXiv:2508.06134 (23 tensor-polarized TMD
structure functions; S_LL/S_TT nomenclature).
[30] W. Cosyn, B. Pire, *Transversity generalized parton distributions
for the deuteron*, Phys. Rev. D 98 (2018) 074020, arXiv:1806.01177;
W. Cosyn et al., polynomiality sum rules for spin-1 GPDs, Phys. Rev. D
99 (2019) 094035; W. Cosyn, B. Pire, L. Szymanowski, PoS DIS2019
(2019) 254, arXiv:1907.08662.
[31] D. Boer, S. Cotogno, T. van Daal, P. J. Mulders, A. Signori,
Y.-J. Zhou, *Gluon and Wilson loop TMDs for hadrons of spin ≤ 1*,
JHEP 10 (2016) 013, arXiv:1607.01654; S. Cotogno, T. van Daal,
P. J. Mulders, *Positivity bounds on gluon TMDs for hadrons of
spin ≤ 1*, JHEP 11 (2017) 185, arXiv:1709.07827.
[32] E. R. Berger, F. Cano, M. Diehl, B. Pire, *Generalized parton
distributions in the deuteron*, Phys. Rev. Lett. 87 (2001) 142302.
[33] F. Cano, B. Pire, *Deep electroproduction of photons and mesons on
the deuteron*, Eur. Phys. J. A 19 (2004) 423, hep-ph/0307231;
A. Kirchner, D. Müller, *Deeply virtual Compton scattering off
nuclei*, Eur. Phys. J. C 32 (2003) 347, hep-ph/0302007.
[34] M. Nzar, P. Hoodbhoy, *Estimation of the double-helicity-flip
deuteron structure function*, Phys. Rev. D 45 (1992) 2264.
[35] Moments and radii: I. Angeli, K. P. Marinova, At. Data Nucl. Data
Tables 99 (2013) 69; recommended Q(⁶Li) = −0.0806(6) fm²,
Q(⁷Li) = −4.00(3) fm² (Pyykkö compilations, via arXiv:2403.06384;
N. J. Stone, At. Data Nucl. Data Tables 111-112 (2016) 1).
[36] B. Wojtsekhowski, *Polarized internal target experiments based on
EIC beams*, arXiv:2406.11480.
[37] M. M. Dalton, A. Deur, C. D. Keith, *Potential for tensor
polarized deuterons in Hall D at Jefferson Lab*, Eur. Phys. J. A 61
(2025) 111, arXiv:2504.21177.
[38] M. M. Dalton, A. Deur, C. D. Keith, N. Fomin, M. Sargsian, *High
precision measurement of φ-nucleon cross section using a tensor
polarized deuteron target*, JLab PAC 53 proposal, arXiv:2508.06481.
[39] D. R. Tilley, C. M. Cheves, J. L. Godwin, G. M. Hale, H. M.
Hofmann, J. H. Kelley, C. G. Sheu, H. R. Weller, *Energy levels of
light nuclei A = 5, 6, 7*, Nucl. Phys. A 708 (2002) 3 (TUNL
evaluation).
[40] EIC second-detector / IR-8 rare-isotope instrumentation:
arXiv:2211.15746 §VIII; arXiv:2602.04636.
[41] LEPS collaboration, coherent φ photoproduction on ⁴He
(b = 23.81 ± 0.95 (stat) GeV⁻², with a sizable one-sided systematic),
arXiv:1711.01095.
[42] STAR collaboration, *Probing the gluonic structure of the deuteron
with J/ψ photoproduction in d+Au ultra-peripheral collisions*,
Phys. Rev. Lett. 128 (2022) 122303, arXiv:2109.07625.
[43] W. Chang, E.-C. Aschenauer, M. D. Baker, A. Jentsch, J.-H. Lee,
Z. Tu, Z. Yin, L. Zheng, *Investigation of the background in coherent
J/ψ production at the EIC*, Phys. Rev. D 104 (2021) 114030,
arXiv:2108.01694.
[44] H1 collaboration, *Inclusive measurement of diffractive
deep-inelastic scattering at HERA*, Eur. Phys. J. C 72 (2012) 2074,
arXiv:1203.4495.
[45] T. Lappi, *Understanding saturation and AA collisions with an eA
collider*, Nucl. Phys. A 830 (2009) 403c, arXiv:0907.4588 ("20–25% of
the events at an EIC could be diffractive"); H. Kowalski, T. Lappi,
C. Marquet, R. Venugopalan, *Nuclear enhancement and suppression of
diffractive structure functions at high energies*, Phys. Rev. C 78
(2008) 045201, arXiv:0805.4071.
[46] A. Accardi et al., *Electron Ion Collider: the next QCD frontier*,
Eur. Phys. J. A 52 (2016) 268, arXiv:1212.1701.
[47] I. Sick, *Form factors and radii of light nuclei*,
arXiv:1505.06924 (R_ch(⁶Li) = 2.589(39) fm).
[48] ECCE consortium, exclusive/diffractive/tagging physics studies,
arXiv:2208.14575.
[49] Polarized ³He/⁴He coherent DVCS at the EIC, arXiv:2606.11491.
[50] EPIOS white paper: G. Atoian et al., *Realizing the scientific
program with polarized ion beams at EIC*, Phys. Rev. C 113 (2026)
060501, arXiv:2510.10794.
