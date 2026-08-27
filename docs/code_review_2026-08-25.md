# Simulation code review and reconstruction audit (2026-08-25)

> **Two of this review's own numbers were corrected when its items were
> implemented (2026-08-26, plans/08 C2 and C3).**  (i) §7 S1 estimates
> that fixing the fast-sim `r1998` moves the moment-solved |A_bag| "only
> ~10%"; measured, it moves **21–25%** (0.318/0.310/0.297 →
> 0.237/0.236/0.235).  (ii) The five R values S1 quotes for "Θ only in
> the log term" (0.34, 0.31, 0.41, 0.20, 0.12) are a property of the
> script's own function, which is not R1998: its third term has no
> counterpart in Abe et al. Eq. (2) and contributes 23–26% of R at the
> two high-x points.  The published three-form R1998 gives 0.350, 0.318,
> 0.374, 0.184, 0.097 there, and at the high-x points the difference
> exceeds the fit's own R_a/R_b/R_c spread — so it is the R1998 numbers,
> not S1's, that belong in the letter.  Everything else in S1 stands, and
> its central prediction (L_5σ 131 → 63–69 and 275 → 152–164 fb⁻¹/u) is
> confirmed at 65.8 and 155.1.  (iii) S1's "Θ … — 6.5–9.5 at x ≤ 0.05" is
> also measured wrong: over 1 ≤ Q² ≤ 130 at x ≤ 0.05 the factor runs
> **6.2–12.9**, falling to 4.45 at the fit's Q² = 0.5 support edge.  The
> qualitative point — Θ is large where the sensitivity is, so multiplying
> it into every term saturates the sum — is exactly right, which is why
> the item was worth implementing.


Scope: the whole simulation code (`evgen/polligen`, `evgen/scripts`,
`fastsim/polli_fastsim`, `fastsim/scripts`, both test suites), audited
under three questions: (1) is each step of the reconstruction chain
consistent with the reference it cites (the reconstruction-chain report,
HJM/JM, Bacchetta et al., the Yellow Report, Jentsch–Tu–Weiss, Maple/ATHENA);
(2) does any reconstruction or analysis path use information a real
experiment cannot measure (Monte-Carlo truth); (3) is every measured
quantity one the expected ePIC design delivers. Method: line-by-line
reading of the reconstruction modules (`reco.py`, `recopseudo.py`,
`coherent.py`, the 5R/7R/6R scripts, the report text), three independent
reviewers over the generator core, the fast-sim package and the
truth-leak question, numerical checks run against the code (§5), and the
full test suites (118 passed: 94 evgen + 24 fastsim).

## 0. Verdict

* **No truth-only quantity enters the reconstructed-level analyses.**
  In money plots 5R, 7R and 6R the event selection, the binning, the
  azimuths, the spin-state ratio and the fits consume only quantities an
  experiment records: the smeared scattered-electron four-vector, the
  nominal beams and crossing angle, the per-fill tensor polarization and
  luminosity share, the Roman-Pot angle pair after the near-beam cut.
  Truth enters only where a real analysis also uses its simulation: the
  response (purity/efficiency, the bin-centering factor K, the
  acceptance-weighted β basis) and the closure references (§1).
* **Where the chain is not yet the experiment's:** the hadronic-method y
  is a Gaussian stand-in smeared around the true y (no hadronic final
  state is generated, plans/04 #21) — at the sweet spots that means an
  absolute resolution of 50–255 MeV on an E − p_z sum of 0.2–1.0 GeV
  (50–125 MeV at spots 1–3),
  which no published ePIC study covers (§2, F2); radiative events,
  backgrounds and the coherent exclusivity selection (Z-ID of the intact
  ⁶Li, vetoes, M_X) are absent, as the report states.
* **One systematic the chain does not model and the report does not
  state:** the spin-state ratio cancels the φ′ acceptance only if both
  spin states see the *same* acceptance. A difference of the cos 2φ′
  harmonic of the efficiency between the m = ±1-rich and m = 0-rich
  samples leaks into Â as (ε₂⁺ − ε₂⁰)/(P₊ − P₀): 10⁻² → 5.6× the
  signal, 10⁻³ → 0.56×, 10⁻⁴ → 6% (§5.1). Fill-by-fill alternation
  therefore needs 10⁻⁴ stability of the acceptance harmonic; bunch-by-bunch
  alternation within a fill (plans/04 #3) makes the cancellation exact.
  The coherent cutout has the same property, and far more sharply than
  this pass estimated: 10⁻³ of the vertical envelope between the samples
  biases a_t by 19% (§5.2, superseded note; plans/08 A1b).
* **Consistency with the references:** the generator kernel, the
  covariant azimuth, the estimator algebra, the Mäntysaari 2a₂
  convention, the Jentsch far-forward windows, the ZEUS u₁/u₂ bounds and
  the Maple/ATHENA 25% resolution are used as the sources state. Two
  numbers in the report's §3 text disagree with the code that produced
  its figures (δQ²/Q² "1–2%" against 5% measured with the 3 mrad
  angular table; the φ′ dilution "1 − 2×10⁻⁵" against 0.99), and the
  tracking angular-resolution table itself has no published source (§2,
  F3). The migration/bin-centering factor K of money plot 7R is
  model-dependent at the 3–11% level between the repository's own Δ
  shapes (F5), larger than the 10-year statistical error.
* **Two items outside the reconstruction chain proper, found by the
  reviewers and confirmed:** (i) the sign of the b₁ (tensor-rate)
  sector is opposite to HJM/HERMES/Cosyn — the code's A_zz = +(2/3)b₁/F₁
  where the reference (Cosyn Eq. 27, the HERMES extraction relation, and a
  direct contraction of HJM's tensor) gives −(2/3)b₁/F₁; magnitudes and
  the Δ sector are unaffected, the sign of Δ itself cannot be checked
  without JM 1989 (§6, G1); (ii) the fast-sim money-delta scripts'
  "R1998" returns R = σ_L/σ_T = 1 at x ≲ 0.1, Q² ≲ 5 GeV² (Θ multiplies
  all three terms), making the July fast-sim L_5σ values ≈ 1.5–2× too
  pessimistic; polligen and the R-plots use a different (toy) R and are
  unaffected (§7, S1). Neither was changed in this pass — both change
  published numbers/signs and are decisions for the author.
* Reviewer findings on the generator core and the fast-sim package are
  in §6–§7; what this pass changed is listed in §9.

## 1. Measurability audit of the reconstructed-level analyses

Every quantity that flows into the 5R/7R (inclusive) and 6R (coherent)
analyses, with where it is produced and whether an experiment has it.
"Response" = built from the simulation and applied to data, as a real
analysis applies its Monte-Carlo corrections.

### 1.1 Inclusive (money plots 5R, 7R)

| Step | Quantity used | Code | Measurable? | Note |
|---|---|---|---|---|
| Beams | E_e, E_ion/u, crossing angle, spin axis φ_S | `recopseudo.RecoResponse.__init__` (`reco.beam_fourvectors`, `spin_fourvector(m.phi_s)`) | yes | nominal beams; the axis is the stable (vertical) spin direction, known from the machine/polarimeter |
| e′ measurement | E′ (EMCal 2%/√E ⊕ 1% or tracking), θ, φ (track direction table) in the lab | `RecoResponse.__init__` lines 126–147 | yes | smearing applied in the lab frame; azimuth resolution σ_θ/sin θ′ |
| Head-on frame | boost + rotation with the nominal crossing angle | `reco.lab_to_head_on` | yes | the standard ePIC transformation |
| Q² | electron method on the smeared k′ | `reco.electron_method` | yes | |
| y | hadronic y = y_true·(1 + 0.25·N(0,1)) | `reco.hadronic_y` | **stand-in** | measurable in principle (Σ/JB/DA from the HFS) but the HFS is not simulated; see F2 |
| x | mixed method Q²_e/(s·y_had) | `reco.mixed_method` | yes | the HERA eΣ method |
| Selection | Q²_e ≥ 1, 0.01 ≤ y_had ≤ 0.95, W²(x_reco, Q²_e) ≥ 10, |η_lab| ≤ 3.5, E′ ≥ 0.5 GeV, ε_eID(η_lab) | `RecoResponse.__init__` lines 169–176 | yes | all on reconstructed quantities |
| Azimuth | φ′ = covariant φ_S from (k, k′_smeared, P, S) | `reco.azimuth_wrt_lepton_plane` on `kp_m` | yes | the true φ′ is used only for the per-event dilution of the *expected* counts, i.e. for generating the pseudo-data |
| Binning | reconstructed (x, Q²) super-bins | `RecoResponse.mask_reco` | yes | `mask_true` is used only in `bin_summary` (purity, efficiency, D) |
| Counts per fill | expected φ′ counts per spin state, Poisson-fluctuated | `RecoResponse.expected_counts`, `measure_inclusive` | data | exact bin integrals of den_f + P_f A cos 2φ′ at the events' true kinematics — the data-generation model |
| P_zz per fill | plan values (+0.6 / −1.2) | `measure_inclusive` (`c.moments()[1]`) | yes (polarimeter) | the 3% polarimetry scale of the bookkeeping is not propagated into 5R/7R (statistical errors only, as labelled) |
| Luminosity shares | analysis assumes 0.5/0.5; truth carries a 10⁻³ offset | `money_cos2phi_reco.py` line 138, `bookkeeping.tensor_flip_plan(rel_lumi_offset)` | yes | the relative-luminosity error is emulated, not assumed away |
| Fit | weighted LSQ of the inverted ratio, cos 2φ′ + sin 2φ′ | `reco.harmonic_ratio_fit` | yes | needs only counts, P_f, L_f |
| Amplitude → Δ | K = Δ_model(x_c, Q²_c)/A_reco-bin(model) | `recopseudo.delta_from_amplitude`, `bin_summary` | response | model-derived migration + bin-centering correction; see F5 |
| F₁, F₂, R | model at the bin centre (inside K) | `money_cos2phi_reco.py` lines 313–316 | response | a real analysis takes them from the same data or from nuclear PDFs |

### 1.2 Coherent (money plot 6R)

| Step | Quantity used | Code | Measurable? | Note |
|---|---|---|---|---|
| Recoil generation | t, φ_t, x_P sampled | `CoherentResponse.__init__` | truth (generation) | x_P is never used downstream: the pots do not resolve it and the analysis bins in reconstructed |t| only |
| Roman-Pot measurement | θ_x, θ_y smeared by the beam divergence; rectangular slot-like cutout; p_T = A p_u·|θ|; φ_t; t = −p_T² (x_L = 1) | `reco.rp_measure` | yes | detector resolution neglected against the beam divergence (Jentsch DIS 2023) — stated |
| Acceptance | `accepted` mask of the smeared angles | `rp_measure`, `CoherentResponse` | yes | events inside the cutout are never used |
| Azimuths | β = φ_t,reco − φ_S; α = φ_e − φ_S | `CoherentResponse` (`beta_reco`), `expected_counts_2d` | yes | α is integrated analytically per event (the e′ azimuth smearing is negligible for it) |
| t bins | reconstructed |t| | `expected_counts_2d(sel = t_reco in bin)` | yes | |
| Counts per fill | exact (α, β) integrals with the true β, t inside, Poisson-fluctuated | `expected_counts_2d`, `measure_coherent` | data | data-generation model |
| Fit basis | acceptance-weighted in-bin means of the true β and the t-template per reco β bin | `CoherentResponse.basis_means` → `reco.basis_2d` | response | the MC-corrected basis a real analysis builds from its simulation |
| u₁, u₂ | taken as known (0.05, 0.02) | `measure_coherent(u_coeffs)` | assumption | in reality measured from the spin-averaged sample of the same data; enters the ratio only through the denominator |
| Truth reference | a_t(t_ref), a_e, a_m | `CoherentResponse.truth_reference` | closure only | plotted, never fitted |
| Exclusivity | intact-⁶Li identification, vetoes, M_X | — | **absent** | 6R assumes a pure coherent sample (plans/06, open #18/#19) |

An independent reviewer traced the same paths (every assignment of
`x, y, q2, phi_true, t_true, beta_true, mask_true` and every input of
`harmonic_ratio_fit` / `harmonic_ratio_fit_2d`) and reached the same
result: `phi_true` is assigned and never read; `mask_true` is read only
in `bin_summary` and the closure test; `dil`, `t_true`, `beta_true` are
read only in the expected-count generation, the response and the
display; per-event spin states are never drawn (the fill enters only
through the population sums of `_fill_arrays`).

**Conclusion of §1.** Run on real data, the 5R/7R/6R analysis code would
need only: per-bunch-crossing spin state, per-fill P_zz and luminosity
share, the reconstructed e′ four-vector, a hadronic-method y, the
Roman-Pot angle pair, and the simulation-derived corrections (purity,
K, β basis). Nothing in the fit paths reads `x`, `y`, `q2`, `phi_true`,
`t_true` or `beta_true` (they appear only in `bin_summary`,
`basis_means`, `truth_reference`, `t_mean_true` and the expected-count
generation).

## 2. Findings, ranked

**F1 — HIGH (systematic not modelled, not stated): acceptance stability
between the spin-state samples.** `reco.spin_state_ratio` cancels a
φ′-dependent efficiency ε(φ′) exactly *because* the same ε multiplies
every fill (docstring, report §3). `RecoResponse.expected_counts` and
`reco.expected_counts_by_fill` accept a single `phi_eff` for all fills,
so the pseudo-experiments can only demonstrate the common-acceptance
case. With a fill-dependent harmonic the bias is analytic and large
(§5.1): Â_bias = (ε₂⁺ − ε₂⁰)/(P₊ − P₀). For the m = ±1-rich / m = 0-rich
pattern (P₊ − P₀ = 1.8) a 10⁻³ difference of the cos 2φ′ harmonic of the
efficiency between the two samples fakes 56% of the signal (Â = 10⁻³),
so fill-by-fill alternation requires the harmonic to be stable to ~10⁻⁴
— tighter than any calibration the report discusses. The bookkeeping
docstring (`tensor_flip_plan`) describes *bunches*, the plan and the
category names say *fills*; plans/04 #3 already lists "bunch-by-bunch
spin patterns" as an unknown of tensor operations. The same holds for
the coherent cutout (§5.2): the fake ⟨cos 2β⟩ of the slot changes by
0.0027 per 1% change of the vertical envelope (10σ_y), i.e. a naive bias
of 0.0015 on a_t ≈ 0.12 (1.3%) — an estimate the template fit refutes by
two orders of magnitude (§5.2 superseded note; plans/08 A1b). *Action:*
state the requirement (bunch-by-bunch alternation, or 10⁻⁴ acceptance
stability and 1% envelope stability) in the report, plans/04 #3 and the
plans/07 systematics list; add a per-fill efficiency option to the
pseudo-experiments so the number is produced by the code, not by hand.

**F2 — HIGH (realism, known #21, now quantified): the hadronic-y
stand-in.** `reco.hadronic_y` smears the *true* y by 25% (Gaussian,
floor 10⁻⁴), independently of the electron measurement. For the sweet
spots of the mid-energy configuration the hadronic E − p_z sum is
Σ_h = 2E_e·y = 0.20, 0.50, 0.22 and 1.02 GeV (spots 1–4), so "25%" is an
*absolute* resolution of 50–255 MeV on Σ_h; the documented 25% comes
from pseudodata smearing at Q² > 1 GeV² without a stated y range
(Maple) and from full simulation at Q² > 100 GeV² (Maple slide 47,
Arratia), where Σ_h is tens of GeV. Arratia et al. Sec. 5 show that a
50 MeV-per-component calorimeter noise already changes the low-y
resolution qualitatively. Everything in 5R/7R (purities 0.63–0.68,
efficiencies 0.38–0.65, the δÂ) rests on this number. *Action:* keep #21
as the top external request (an ePIC full-simulation Σ_h/DA resolution
at Q² ≈ 1–3 GeV², Σ_h ≈ 0.2–1 GeV, e + light ions); quote 5R/7R for the
15–30% band *and* for an absolute-noise model (e.g. σ(Σ_h) = 25% ⊕
50 MeV) so the reader sees the sensitivity; the low-energy configuration
(y = 0.05–0.25 at the same (x, Q²)) is the fallback the report already
names.

**F3 — MEDIUM (report ↔ code inconsistency; unsourced table).** Report
§3 step 2 and the Schematic 2 box say δQ²/Q² ≈ 1–2% for the electron
method. With the angular table the code uses
(`reco.tracking_angular_resolution`: 3 mrad for −3.5 < η < −2.5, 2 mrad
for −2.5 < η < −1), the response gives rms δQ²/Q² = 5.2% at spots 1–2
(θ′ = 107 mrad → cot(θ′/2)·δθ′ = 18.7 × 3 mrad), 2.9% at spot 3, 1.6% at
spot 4 (§5.3); with the table divided by three it is 2.1% and the
spot-1 purity rises from 0.64 to 0.70. The table is copied from
`fastsim/scripts/money_delta_20260729.py`, whose docstring says the
tracking numbers were "provided in parent communication" — no
publication. Likewise §3 step 4 says the track angular resolution
dilutes cos 2φ′ by "exp(−2σ_φ²) ≈ 1 − 2×10⁻⁵": the code correctly uses
σ_φ = σ_θ/sin θ′ ≈ 20–30 mrad at the sweet spots and the measured
dilution is 0.985–0.999 (`bin_summary["dilution_phi"]`, §5.3) — still
negligible, but 100× the quoted number. *Action:* correct the two
statements (done in this pass), and either cite an ePIC angular
resolution for the backward disks or label the table as a placeholder in
`reco.py` (done).

**F4 — MEDIUM (realism): the EMCal resolution is applied at all η.**
`RecoModel(energy="emcal")` uses 2%/√E ⊕ 1% (the backward PbWO₄
specification) for every electron. The Yellow Report requirements are
2%/√E ⊕ 1–3% for −3.5 < η < −2, 7%/√E ⊕ 1–3% for −2 < η < −1,
(10–12)%/√E ⊕ 1–3% in the barrel and 7%/√E ⊕ 1–3% forward; the ePIC
barrel imaging calorimeter and the forward endcap are coarser than the
backward crystals. The sweet-spot electrons sit at η = −2.9 to −1.6
(Table 1 of the report), so the headline numbers are unaffected, but the
amplitude-vs-x panels (5R right, 7R) include bins whose electrons reach
the barrel, where 2%/√E ⊕ 1% is optimistic by a factor 3–5 in the
stochastic term; there tracking (0.05%·p ⊕ 0.5%) is the better
measurement, which `energy="best"` selects. *Action:* make the
calorimeter term η-dependent per the Yellow Report table and use "best"
by default (recommended; not changed in this pass because it changes
published figures). *Done 2026-08-25 (plans/08 A7):*
`emcal_resolution(..., eta=)` carries the Yellow Report table and
`RecoModel.emcal_eta_table` switches it on, default off so nothing
published moves. Measured with `--syst-scan`: the table changes Δ̂ by
0.000% at spots 1–3 (backward electrons) and 0.5% at spot 4
(η = −1.64, the transition region), against an MC noise floor of
0.13–0.21%.

**F5 — MEDIUM (method): the 7R migration correction is model-dependent
at the 3–11% level and the closure is by construction.** The factor
K = Δ_model(x_c, Q²_c)/A_reco-bin(model) is evaluated with the same Δ
model that generated the pseudo-data, so 7R points sit on the injected
curve by construction (the report says so). Recomputing K with the
repository's alternative shapes (§5.4): K_A/K_B − 1 = −6%, +7%, −11%,
+4% (moment_B) and +10%, +3%, +11%, +4% (toy) at the four sweet spots —
larger than the 10-year statistical errors (0.4–1%) and comparable to
the 1-year ones (1–3%). This is the usual model dependence of a
bin-by-bin unfolding with purity 0.64 (36% of the events come from
neighbouring true bins). *Action:* in the letter, either iterate K with
the fitted shape (two iterations suffice for a smooth Δ(x)) or quote the
spread over shapes as a systematic; a migration-matrix unfolding across
the x bins of a Q² slice removes most of it. Cheap with the existing
code (`bin_summary` under a second kernel).

**F6 — MEDIUM (realism, coherent): the 6R sample is background-free and
Z-identified by assumption.** `CoherentResponse` tags every intact ⁶Li
above the envelope; the α + d breakup fragments have the same rigidity
(A/Z = 2) and land in the same pots, and no Z-ID, ZDC/OMD/B0 veto or M_X
selection is simulated (report §4 step 3, plans/06 #18/#19). Because the
incoherent yield is spin-blind it dilutes rather than fakes the tensor
ratio, but the dilution (and its |t| dependence, the breakup slope being
~10× flatter) is not in the 6R error bars. u₁, u₂ are taken as known
and the relative-luminosity offset is not exercised in 6R (it is in 5R).
*Action:* label the 6R figure "background-free, Z-identified sample";
carry the plans/06 background budget into WP5.

**F7 — LOW (report ↔ code wording): fills versus bunches.** See F1. The
report's Schematic 2 says "interleaved bunches", `tensor_flip_plan`
names "fills"; the physics of the cancellation depends on which. *Action:*
say "bunch-by-bunch pattern (required for the acceptance cancellation)"
consistently (done in the report and the plan).

**F8 — LOW (documentation): systematic errors are not propagated in the
R-plots, as labelled.** The polarimetry scale (3%, `bookkeeping`), the
axis tilt (sin²θ_S), radiative corrections and the F₂/R model are
outside the 5R/7R/6R error bars; the report lists them in §3 step 7.
No action beyond keeping the "statistical errors only" label.

**F10 — LOW (report ↔ code): the coherent pseudo-experiment has no
electron side.** `CoherentResponse` generates no electron: α is
integrated analytically per event and β = φ_t,reco − φ_S is the lab
difference, not `azimuth_wrt_lepton_plane`; there is no (x, Q², x_P)
reconstruction, only reconstructed |t| bins integrated over
x_P ∈ [10⁻³, 10⁻²]. Numerically harmless (the recoil shortcut error is
< 0.2 mrad, report §4; test_reco bounds it below 1 mrad), but the
report's §5.1 said "compute α and β covariantly" (independent reviewer,
M5). *Action:* wording corrected in §5.1 (done).

**F11 — LOW (report formula): the ratio was written with Y_f = N_f/L_f
and a luminosity-weighted P̄.** That mixture is exact only for equal
luminosities; the code's form — raw counts N_f with
w_f = P_f − P̄, P̄ = ΣL_fP_f/ΣL_f, so that ΣL_fw_f = 0 — is exact for any
shares (derivation in the `spin_state_ratio` docstring; independent
reviewer). *Action:* report §3 formula and Schematic 2 corrected to the
code's form, T(φ′) = c + d cos 2φ′ fitted after the exact inversion
(done).

**F12 — LOW (bookkeeping): `tensor_flip_plan` defaulted to φ_S = 0
while `RecoModel.phi_s = π/2`.** Harmless — the kernel's a₂ depends on
θ_S only and φ_S enters `dsigma`, which the reco layer never calls — but
inconsistent. *Action:* default changed to π/2 (vertical stable-spin
direction), docstring states the bunch-by-bunch requirement (done; 94
evgen tests pass).

**F13 — LOW (report wording, several):** "weight σ_c/400" is σ_c/n_kept
(the cell rate is conserved after the generator-acceptance rejection);
the Table 2 caption's D is ε_eID-weighted and φ′-diluted in the
numerator only (10⁻³ effect); §3 step 7's "D_φ from F₂, R of the same
data" describes the experiment while the code uses the model inside K
(stated in §7.1). *Action:* first two corrected (done).

**F9 — LOW (numerical hygiene, verified OK).** Empty φ bins inside the
cutout get zero weight (`spin_state_ratio`: var = ∞); the ratio
inversion `_ratio_to_modulation` converges in 4 iterations for |P̄T| ≪ 1;
the beam-energy-spread emulation moves Q²_e and (1 − y_e) coherently;
the generator scenario is loose enough for migration into every analysis
cut (y ≥ 0.004 against the 0.01 cut: a +43% fluctuation of 0.007 is
1.7σ at 25%); `eps_eid` is applied on the reconstructed η and vanishes
beyond |η| = 3.5, consistent with the selection.

## 3. Consistency with the references (chain ↔ report ↔ sources)

| Item | Source | Code | Status |
|---|---|---|---|
| Master formula: cos 2φ term −(1−y)/y²·c_m sin²θ_m Δ against F₁ + (1−y)/(xy²)F₂ | HJM NPB 312:571, JM PLB 223:218 (no local copy) | `polligen.xsec`, `polli_fastsim.asymmetries.a_cos2phi` | y-structure consistent: with F₂ = 2xF₁ the ratio is −(1−y)/[1+(1−y)²]·Δ/F₁ = −(ε/2)Δ/F₁, the standard TT-interference depolarization (numerically exact incl. R); m-dependence c_m sin²θ_S cos 2(φ−φ_S) verified from the density matrix; **overall sign of Δ unverifiable locally** (§6 G1) |
| Tensor-rate (b₁) sector: A_zz = (2/3)[b₁ + (1−y)/(xy²)b₂]/D_φ | Cosyn et al. Eq. (2), (27) (`refs/2410.12764v1.pdf`); HERMES; HJM tensor contraction | `asymmetries.azz`, `xsec.amplitudes` (t_geo) | **sign opposite to the reference** (§6 G1); magnitude, b₂ = 2xb₁ (Cosyn Eq. 18), geometry (3cos²θ − 1)/2 consistent |
| R = σ_L/σ_T | SLAC R1990/R1998 (Θ only in the log term) | `structure.r_sigma_lt` (toy, polligen); `money_delta_realistic.r1998` (fast-sim scripts) | toy: placeholder of the right magnitude; **fast-sim `r1998` wrong (R = 1 at low x)** (§7 S1) |
| Far-forward OMD window | Jentsch–Tu–Weiss Table I: ξ 0.45–0.65 | `farforward.OMD_R_WINDOW = (0.40, 0.60)` | **inconsistent** (§7 S7), low impact |
| α-tag near-beam cut | YR: 10σ beam envelope (angular) | `farforward` 0.20 GeV proton p_T threshold applied to the α | **not scaled to the α momentum** (§7 S8): README 3–9% → 0.7–2.7% |
| Covariant azimuth φ_S (Eqs. 2.4–2.5) | Bacchetta et al. JHEP 02 (2007) 093 (refs/hep-ph_0611265.pdf) | `reco.azimuth_wrt_lepton_plane`, `transverse_part` | consistent; = φ_e − φ_S exactly for a massless target (test), O(γ²) for ⁶Li |
| Head-on frame, 25 mrad crossing, ePIC axis = electron beam | ePIC convention | `reco.head_on_to_lab / lab_to_head_on` | consistent; only odd harmonics of the e′ azimuth change (test) |
| Electron-method resolution (1−y)/y·[δE′/E′ ⊕ tan(θ′/2)δθ′ ⊕ δE_e/E_e] | standard | `reco.electron_method_resolution` | consistent |
| Σ / JB / mixed method | Bassler–Bernardi (refs/hep-ex_9412004.pdf) | `reco.hadronic_y`, `mixed_method` | consistent in form; resolution is a stand-in (F2) |
| Hadronic-y resolution 25% | Maple ePIC seminar slide 44/47; ATHENA Fig. 22; Arratia Fig. 5 | `RecoModel.y_had_res = 0.25` | number as documented; regime mismatch (F2) |
| EMCal 2%/√E ⊕ 1% | YR backward requirement | `reco.emcal_resolution` | consistent backward only (F4) |
| Tracking momentum/angle tables | none ("parent communication") | `reco.tracking_resolution`, `tracking_angular_resolution` | **unsourced** (F3) |
| ε_eID(η) anchors | ATHENA Table 5, ECCE Sec. 3.5.2 | `reco.eps_eid` | as documented; no ePIC curve exists |
| Spin-state ratio R = σ_P²T/(1 + P̄T), δÂ = √(2/N)/σ_P/dilution | derived in the report | `reco.spin_state_ratio`, `harmonic_ratio_fit`, `err_harmonic_ratio` | consistent; verified numerically (§5.1: unbiased under a common ε(φ′), analytic error reproduced) |
| Roman-Pot 10σ cut, pT_cut = 10σ_θ·A·p_u | YR (0.20/0.45 GeV for 275 GeV protons); Jentsch DIS 2023 slide 20 (beam effects dominate) | `reco.tag_pt_cut`, `rp_measure` | consistent; σ_θ proton-derived placeholder (plans/04 #11/#20) |
| Far-forward windows (B0 5.5–20 mrad, OMD 0–5 mrad ξ 0.45–0.65, RP 0–5 mrad ξ 0.6–0.95, ZDC 0–4 mrad) | Jentsch–Tu–Weiss Table I (refs/2108.08314.pdf, verified 2026-08-25) | `polli_fastsim.farforward` | consistent |
| Slot-like cutout (wide in x, tight in y) | Jentsch DIS 2023 slide 15 | `rp_measure(cut_scale_xy=(2.5, 1))` | shape as documented; dimensions illustrative (#20) |
| Deformation coefficient c₂ = 2a₂ | Mäntysaari et al. Eq. (9) (refs/2408.13213v1.pdf) | `coherent.cos2phi_coefficient_deformation` | consistent (money plot 6 injects a₂: conservative ×2, stated) |
| u₁, u₂ y-factors and bounds | Nikolaev–Pronyaev–Zakharov Eq. (1); ZEUS NPB 816:1 Sec. 10.2 | `money_cos2phi_coherent_reco.py` defaults | consistent; defaults at the 1σ edge |
| t = −(p_T² + M²x_P²)/(1 − x_P), t_min | light-cone kinematics | `reco.recoil_fourvector`, `t_from_fourvectors` | consistent (test) |
| Tensor-SF leakage: cos φ at O(γ), cos 2φ at O(γ²) | Cosyn et al. Eqs. 17d–e (refs/2410.12764v1.pdf) | report §2 (not in the code; ≤ 3% at the sweet spots) | consistent |

## 4. What each script is (truth-level projection vs reconstructed level)

| Script | Level | Estimator | Notes |
|---|---|---|---|
| `evgen/scripts/money_cos2phi.py` (5), `money_cos2phi_coherent.py` (6), `money_delta_extraction.py` (7) | truth (x, Q²) bins, exact statistics | single-fill cos 2φ′ fit | statistically faithful; acceptance/relative-luminosity not modelled (report §5 verdict); 6 injects a₂ (×2 conservative) |
| `evgen/scripts/money_cos2phi_reco.py` (5R, 7R), `money_cos2phi_coherent_reco.py` (6R) | reconstructed bins; smeared e′; RP emulation | two-fill spin-state ratio | §1 of this review |
| `evgen/scripts/reco_chain_figures.py` | diagnostics | both | Figure 1–2 of the report |
| `evgen/scripts/money_tagged_azz.py`, `tagged_polarimetry_7li.py`, `closure_fom.py`, `phase_space_bins.py` | truth | — | generator closure and tagged mode (§6) |
| `fastsim/scripts/money_delta_20260729.py` and predecessors | electron-method smearing + (x, Q²) migration, truth φ′ | single-fill | the "reco mask" is on reconstructed bins (§7); no hadronic method, so the y ≈ 0.01 bins are electron-method only |
| `fastsim/scripts/money_b1.py`, `money_polemc.py`, `tagging_acceptance.py`, `coverage_and_stat_maps.py`, `phase_space_map.py` | truth | — | §7 |

## 5. Numerical checks performed for this review

All scripts are in the session scratchpad; the numbers are reproducible
from the repository code with the defaults of the R-scripts
(e10 × ⁶Li50, P_zz = +0.6/−1.2, moment_A, 25% hadronic y).

### 5.1 Fill-dependent φ′ efficiency (F1)

`reco.harmonic_ratio_fit` on exact expected counts (N = 10⁹, no Poisson
noise), A_true = 10⁻³ per unit P_zz, 24 bins:

| Efficiency | Â | bias / signal |
|---|---|---|
| common 1 + 0.03 cos 2φ′ | 1.000×10⁻³ | 0 (−7×10⁻⁹) |
| 0.03 (m = ±1-rich) vs 0.02 (m = 0-rich) | 6.56×10⁻³ | 5.6 |
| 0.03 vs 0.029 | 1.56×10⁻³ | 0.56 |
| 0.03 vs 0.0299 | 1.06×10⁻³ | 0.06 |

Analytic: bias = (ε₂⁺ − ε₂⁰)/(P₊ − P₀) = 0.01/1.8 = 5.6×10⁻³ (exact
match). The single-fill fit under a common efficiency is biased by
ε₂/P_zz = 0.05 (report Figure 1d) — the ratio removes that entirely, and
replaces it by the *difference* between the samples.

### 5.2 Cutout stability (F1, coherent)

`reco.rp_hole_acceptance(B = 50, c_x = 2.5 c, c_y = c)` for the mid-energy
envelope c = 0.219 GeV (acceptance 2.85% analytic and unsmeared — the smeared MC
of 6R gives 3.2% —, fake ⟨cos 2β⟩ = −0.772 about
the horizontal, +0.772 about the vertical spin axis): a +1% change of
c_y moves the fake coefficient by −0.0027 (→ a_t bias −0.0015 for
P₊ − P₀ = 1.8, i.e. 1.3% of a_t = 0.117); +1% of c_x: no change (the
slot's x-sides are outside the spectrum); +5% of both: −0.0132 (→ −0.007,
6%). The envelope is set by the beam size at the pots (β-function,
emittance, orbit), which can differ between fills.

*Superseded 2026-08-25 (plans/08 A1b).* This analytic estimate is the
naive one and understates the effect by about two orders of magnitude.
Feeding a per-fill cutout through the template fit
(`CoherentResponse.with_cut`, `measure_coherent(responses=...)`) gives
+19% on a_t for a **10⁻³** relative change of the vertical half-height
and +169% for 1%, with exact closure at zero. The cause is in the fit,
not in ⟨cos 2β⟩: the slot leaves only 12 of 24 β bins live, where the
t-shape template is 99% anti-correlated with the constant. The
requirement is bunch-by-bunch alternation or ≈10⁻⁴ envelope stability;
the higher-|t| bins are far less exposed (+3.9%, +0.04%).

### 5.3 Resolutions and purities of the response (F3)

`RecoResponse` with the script defaults except 200 events per cell (the
script uses 400), events of
the true super-bin that pass the selection:

| Spot (x, Q², y) | rms δQ²/Q² | rms δy/y | rms δx/x | purity | efficiency | D | ⟨cos 2Δφ′⟩ |
|---|---|---|---|---|---|---|---|
| 1 (0.056, 1.14, 0.010) | 0.052 | 0.23 | 0.20 | 0.64 | 0.40 | 0.905 | 0.991 |
| 2 (0.022, 1.14, 0.025) | 0.052 | 0.24 | 0.34 | 0.63 | 0.60 | 0.979 | 0.997 |
| 3 (0.141, 3.13, 0.011) | 0.029 | 0.24 | 0.20 | 0.66 | 0.38 | 0.910 | 0.985 |
| 4 (0.141, 14.3, 0.051) | 0.016 | 0.24 | 0.36 | 0.68 | 0.65 | 0.965 | 0.999 |

With the angular table divided by three (1 mrad at η ≈ −2.9): δQ²/Q² =
0.021, 0.022, 0.015, 0.012 and purities 0.70, 0.69, 0.67, 0.69. The
Table 2 numbers of the report (purity 0.65/0.64/0.66/0.68, D
0.907/0.986/0.910/0.969) are reproduced within the MC statistics of this
check.

### 5.4 Model dependence of K (F5)

K = Δ_model(x_c, Q²_c)/A_reco-bin(model) at the four spots, same events:

| Model | spot 1 | spot 2 | spot 3 | spot 4 |
|---|---|---|---|---|
| moment_A (injected) | −7.68 | −19.67 | −2.41 | −2.39 |
| moment_B | −8.16 | −18.43 | −2.70 | −2.29 |
| toy (Δ/F₁ = 10⁻³) | −6.96 | −19.09 | −2.18 | −2.31 |
| K_A/K_B − 1 | −0.059 | +0.068 | −0.106 | +0.043 |
| K_A/K_toy − 1 | +0.104 | +0.030 | +0.109 | +0.035 |

### 5.5 Test suites

`python -m pytest evgen/tests fastsim/tests -q`: 118 passed (2026-08-25).

## 6. Generator core (`evgen/polligen`: xsec, spin, sample, bookkeeping, estimators, coherent, tagged)

Independent reviewer (83 evgen tests run, three numpy check scripts, three
local PDFs text-extracted); every item below was re-checked by the lead
reviewer against the code lines and, where stated, re-derived.

**G1 — HIGH (convention against the reference): the sign of the b₁
sector is opposite to HJM/HERMES/Cosyn.** `xsec.py:213,217` implement
F₁,eff(m) = F₁ + (2/3)a_m b₁ with a_m = ¼c_m(3cos²θ_S − 1), c_m = (1, −2, 1),
and `asymmetries.azz` (lines 41–47) returns +(2/3)K/D_φ, so for the axis
along the photon the code has F₁(±1) = F₁ + b₁/3, F₁(0) = F₁ − 2b₁/3 and
A_zz = +(2/3)b₁/F₁. Contracting HJM's tensor term −b₁r_μν,
r_μν ∝ [|q̂·E_m|² − ⅓]g_μν, with the explicit polarization vectors gives
|q̂·E|² − ⅓ = −(2/3)a_m for every m and θ (lead reviewer's derivation,
numerically confirmed by the reviewer for θ ∈ {0, 0.7, π/2}), i.e.
F₁(±1) = F₁ − b₁/3, F₁(0) = F₁ + 2b₁/3 and A_zz = −(2/3)b₁/F₁. That is
exactly Cosyn et al. EPJ A 61 (2025) 83 (`refs/2410.12764v1.pdf`),
Eq. (2) A_T = [dσ(+1) + dσ(−1) − 2dσ(0)]/Σ — the code's estimator — and
Eq. (27) "A_T = −(2/3) b₁/F₁. This is the relation that was used in the
b₁ extraction of the HERMES result" (verified verbatim in the local
PDF), with HERMES' b₁ > 0 at low x measured through A_zz < 0. The code
follows the transcription of docs/Discussions.pptx p. 5 recorded in the
`asymmetries.py` header; the tests cannot see it because
`test_xsec_identity.py` compares `xsec` to `asymmetries`, which share the
transcription. *Impact:* the sign of every A_zz curve against a b₁ model
(money_b1, the tagged A_zz plots, the transverse-fill rate shift κ);
magnitudes, error bars and the whole Δ/cos 2φ sector are untouched. The
m-dependence of the cos 2φ term (c_m sin²θ_S cos 2(φ − φ_S)) is verified
independently from the lab density matrix, ρ_{+1,−1} = (c_m/4)sin²θ_S
e^{−2iφ_S} (reviewer, 3×10⁻¹⁶); the overall sign of Δ relative to the
modulation could not be checked because HJM NPB 312:571 / JM PLB 223:218
are not available locally. *Action (decision for the author):* flip the
sign of `t_geo` in `xsec.amplitudes` and of `asymmetries.azz` (or
document b₁ ≡ −b₁^HJM everywhere), add a test against Cosyn Eq. (27),
and obtain JM 1989 to fix the Δ sign convention before any sign of an
extracted Δ is quoted. Not changed in this pass.

*Half done 2026-08-25 (plans/08 B2).* The sign is now the single
constant `asymmetries.TENSOR_LL_SIGN`, imported by `xsec`, and
`evgen/tests/test_tensor_convention.py` is the repository's first
EXTERNAL tensor test:

    A_zz(θ_S = 0) · (1 + ε(y) R) = TENSOR_LL_SIGN · (2/3) b₁/F₁,
    ε(y) = (1−y)/(1−y+y²/2)

exact at every y, F₂-free, and pinned at six (x, Q²) points.  (The
widely-quoted form that also divides by [1+2(1−y)/y²] and multiplies by
[1+2(1−y)(1+R)/y²] double-counts R — those brackets are identically
1+εR — and misses by 1.17.)  Demonstrated to be non-blind: flipping the
constant leaves the other 125 tests green and fails this one with the
decision spelled out in the assertion message.  The decision itself
(and the Δ sign, which needs JM 1989) is still the author's.

**G2 — MEDIUM (known, now quantified): the kernel is massless — the
O(γ²) cos 2φ and O(γ) cos φ contributions of b₁–b₄ to a transverse fill
are absent.** From Cosyn Eqs. (10), (17d–e), (18) with T_TT = Q/2 for a
transverse axis, the b₁-induced amplitude has the form of the Δ term
with Δ_fake = γ²b₁/6 (b₃ = b₄ = 0). At Q² = 1.14 GeV²: γ² = 0.0077
(x = 0.05), 0.031 (x = 0.1), 0.12 (x = 0.2); with the toy b₁/F₁ ≈ 0.02
and Δ/F₁ = 10⁻³ this is 3% at x = 0.05 and 40% at x = 0.2 (0.1–2% for
Δ/F₁ ≈ 0.02). The report's §2 statement (≤ 3% at the sweet spots) is
correct for the sweet spots only. *Action:* add the Eq. (17e) term to
the kernel and subtract it with the A_zz measured in the same data (κ);
list it in the plans/07 systematics.

**G3 — LOW:** the φ accept–reject (`sample.py:229–238`) checks only
1 + w_avg > 0, not |a₁| + |a₂| ≤ 1 + w_avg; for unphysical amplitudes it
silently samples max(W, 0) (reviewer: Δ = 3F₁ gives 2⟨cos 2φ′⟩ = −1.117
vs a₂ = −1.232). Harmless at physical 10⁻³–10⁻² amplitudes; add the
guard. *Done 2026-08-25 (plans/08 B1):* `InclusiveSampler._density_min`
computes the EXACT minimum over φ of 1 + A cos φ + B cos 2φ (the vertex
of the quadratic in cos φ, not the envelope 1 + |A| + |B|), and
`_amplitudes` raises with the offending cell.  `phi_histogram_pseudo`
carries its own guard: it bypasses `_amplitudes` and is the path behind
the headline figures.  Every repository scenario keeps a margin > 0.9,
so nothing existing fires.

**G4 — LOW:** the single-fill estimators Â = a₂/P_zz are biased by the
transverse-fill tensor rate shift when b₁ ≠ 0 (Â/A = 1.0015 at
P_zz = +0.6, 0.997 at −1.2 for the toy b₁); the two-fill ratio absorbs κ
exactly (reviewer's MC: const −0.0150 recovered for −0.015 injected).

**G5 — LOW:** cell rates are σ(centre) × full cell area while in-cell
positions are redrawn inside the acceptance: boundary cells are
overstated (rate-weighted 1.6% on a 16×12 grid, less on the default
grids).

**G6 — LOW:** fixed default seeds in `sample.py`, `tagged.py`,
`coherent.py` — repeated calls without `rng` return identical events;
all tests and scripts pass an rng.

**G7 — documented:** `tagged.py` is an instant-form impulse
approximation (non-relativistic |A|²k² sampling, on-shell spectator
boost), not the Cosyn–Weiss light-front formalism (no α_p, no
1/(2 − α_p)² flux factor, no x → x_n rescaling, no pole extrapolation,
no FSI) — as its docstring says; nothing unmeasurable is assumed beyond
the model (the struck-cluster spin populations are model quantities,
as in CW).

**G8 — documented scenario:** coherent `f₀ = 0.04`, `x_coh = 0.01`
written against Bjorken x rather than x_P, the constant proton p_T cut
coexisting with the angular cut (both reported by `project_coherent`;
0.668/0.0909/1.5×10⁻⁸ at 20.5/50/137.5 GeV/u reproduced), t_min
neglected in `recoil_lab` (−15% at x_P = 0.01, as the docstring states),
`a2_deformation` a population average (RATE_WEIGHT_SYST = 0.73). The
2a₂ convention is verified verbatim from Mäntysaari Eq. (9).

**G9 — documented:** `tensor_flip_plan` carries one `pzz_true` for both
samples (the m = 0-rich sample is assumed to have exactly −2× the P_zz of
the ±1-rich one: same source purity).

**G10 — LOW (doc):** `evgen/README.md` said 93 tests; 94 are collected
(corrected).

*Verified correct (reviewer, with method):* master-formula transcription
modulo G1 (reading + identity tests); unpolarized normalisation
4πα²/(xQ⁴)[(1 − y + y²/2)F₂ − y²F_L/2] ≡ (4πα²y/Q²)D_φ/(sx) (by hand);
cos 2φ y-dependence = −(ε/2)(Δ/F₁)/(1 + εR) to machine precision
(numerical); per-nucleon bookkeeping; spin conventions
(P_zz = ⟨3J_z² − 2⟩, Wigner-d, rotation, ρ coherences); two-fill
populations (0.433/0.133/0.433 and 0.133/0.733/0.133); sampler closure at
10⁶ events (2⟨cos 2φ′⟩ = −0.0360 ± 0.0014 vs −0.0342 injected, P_zz = 0.6;
+0.0677 ± 0.0014 vs +0.0681 at −1.2; sin 2φ′, cos 4φ′ consistent with 0);
`phi_histogram_pseudo` exact; estimator variances (std/√(2/N) = 1.02 over
2000 trials; binned-fit dilution sin(2w)/(2w)); the two-fill ratio
(1500 trials, N = 4×10⁶: amp 0.009987 ± 0.000785 for 0.01 injected,
analytic 0.000795; single-fill fake +0.046 under a 3% harmonic);
coherent formulas (B = R²/3, truncated-exponential t sampling,
exp(−Bp_T²), ⟨t⟩_tag = p_T² + 1/B); tagged internals (θ_lm, CG
normalisation, boost = fastsim's). *Untested formulas:* the b₁/Δ signs
against any external source (G1); the O(γ) tensor terms (absent);
`effective_modulation` with b₁ ≠ 0; the positivity bound;
`with_perp=True` sampling; `tag_acceptance_angular`;
`recoil_lab(x_pom ≠ 0)`; `populations_maxent` against the
spin-temperature relation.

## 7. Fast-simulation package (`fastsim/polli_fastsim`, `fastsim/scripts`)

Independent reviewer (24 fastsim tests run; numerical checks with the
installed CT18NLO/NNPDFpol11 grids; the EPPS21 ⁶Li grid used by the
production script is not installed here); items S1, S2, S7 re-verified
by the lead reviewer.

**S1 — HIGH (bug): the "R1998" parameterisation of the fast-sim
money-delta scripts returns R = σ_L/σ_T = 1 over the region that drives
the cos 2φ sensitivity.** `money_delta_realistic.py:121–136` (copied
verbatim into `money_delta_20260729.py:238–253`) multiplies the factor
Θ(x, Q²) = 1 + 12[Q²/(Q² + 1)][0.125²/(0.125² + x²)] — 6.5–9.5 at
x ≤ 0.05 — into all three terms, whereas in the SLAC R1990/R1998 fits it
multiplies only the 1/ln(Q²/0.04) term. Re-evaluated by the lead
reviewer: x = 0.005, Q² = 2 → 2.17 (clipped to 1.000); 0.01, 2.5 → 2.04
(1.000); 0.05, 1.14 → 2.17 (1.000); 0.1, 5 → 1.00; 0.3, 7.4 → 0.27; with Θ
only in the log term the same coefficients give 0.34, 0.31, 0.41, 0.20,
0.12 (world data 0.2–0.35). Since F₁ = F₂/(2x(1 + R)) and Δ ∝ F₁ in the
scenarios, the per-bin amplitude is ×0.63 too small where R = 1, and the
reviewer's re-run of `sig2_per_fb_at_realistic` moves the July central
L_5σ from 131 to 63–69 fb⁻¹/u (mid) and 275 to 152–164 (top) with any
sensible R. The fast-sim notes' own "R1998 vs Christy–Bosted ratio
≈ 2.0" (notes 2026-07-20) is this bug's signature. *Scope:* the
fast-sim notes of July 16–29 and any L_5σ / peak-bin significance
quoted from them; the moment-solved |A_bag| only ~10% (R ≈ 0.27 at
x ≈ 0.3). *Not affected:* polligen (its kernel uses
`structure.r_sigma_lt` = 0.18/(1 + Q²/50), a placeholder of the right
magnitude), hence money plots 5/6/7 and 5R/7R/6R. *Action:* restrict Θ
to the log term (or implement the full published R1998 with its three
functions), re-run the 20260729 production (needs the EPPS21 grid and the
`np.trapezoid` shim, S11), and correct the notes. Not changed in this
pass.

**S2 — HIGH (convention):** the A_zz sign, as G1
(`asymmetries.py:41–47`, docstring lines 13–17); `money_b1.py` plots
|A_zz|, so the projections are sign-blind.

**S3 — MEDIUM (realism, the reason for the R-chain):** 62%/72% (mid/top)
of the analytic significance Σa²P²N/2 comes from y < 0.1 (40/54% from
y < 0.05), where the electron method used for the acceptance has
δy/y ≈ 1 (reviewer: 1.2 at y = 0.028 with the script's own tracking
model) and no hadronic method exists in the fast-sim; the fast-sim
L_5σ values presuppose hadronic-method reconstruction at y ≈ 0.01,
which is what §1–§2 of this review examine.

**S4 — MEDIUM (stale documentation):** `fastsim/README.md:46–48`
"5σ on Δ/F₁ = 10⁻³ … at ~25–37 fb⁻¹/u" is reproduced with P_zz = 0.80
and no 2-of-6 dilution (`sig2_per_fb_toy(pzz=0.8)` → 31.3/35.9); the
scripts' central P_zz = 0.267 gives 131/275 fb⁻¹/u (63–69/152–164 after
S1), above the 1–100 fb⁻¹/u "plausible program" band the plot draws.

**S5 — MEDIUM (assumption):** the "Cluster d, folded" band edge
P_zz = √(2/6)·0.70 = 0.40 in `money_delta_realistic.py:175–185` is the
polarization of a perfectly tagged deuteron-cluster sample; for an
inclusive measurement a fraction f of polarized nucleons gives
P_eff = f·P = 0.23 (below the "Cloet" 0.267), since sig² ∝ (fP)²N, while
√f·P holds only for f·N tagged events — inconsistent with the package's
own ≤ 3% α-tag acceptance.

**S6 — LOW:** `money_b1.py` draws b₁ curves with the ×0.87 deuteron
transfer factor and no 2/6 factor against errors computed with
P_zz = 0.6/0.8 and no 2/6 factor, while the Δ line folds 1/3 into P_zz
and omits the 0.87 (`delta_models.py:39–50` acknowledges the mismatch):
the b₁ signal-to-error ratio is ×3 inconsistent with the Δ convention.

**S7 — LOW (reference):** `farforward.py:27` OMD window R ∈ [0.40, 0.60]
against Jentsch–Tu–Weiss Table I ζ = 0.45–0.65 (`refs/2108.08314.pdf`,
verified; the paper's symbol is ζ, Eq. 58); RP 0.6–0.95, B0 5.5–20 mrad,
ZDC 4 mrad agree. Affects only R ≈ 0.5 fragments. *Done 2026-08-25
(plans/08 C1):* window corrected to [0.45, 0.65). The ⁷Li spectator
proton (R = 3/7 = 0.4286) now falls below the low edge — optics-limited,
not "lost by construction".

**S8 — MEDIUM (realism):** the fast-sim α-tag applies the 0.20 GeV
*proton* p_T threshold (`farforward.py:37,59`) to a 4 × 137.5 GeV α; the
YR cut is an angular envelope (0.73 mrad), i.e. 0.40 GeV (0.90 GeV
high-divergence) for the α — the same scaling the reconstruction report
applied to the coherent ⁶Li. Reviewer's re-run: α-tag acceptance
3.2/6.3/8.9% (β = 0.2/0.3/0.4 GeV) → 0.7/1.7/2.7% with 0.40 GeV, so the
README's "3–9%" is ≈ 3× optimistic; no σ_θ exists in the fast-sim.
*Done 2026-08-25 (plans/08 C1):* `Optics` now carries `sigma_theta` and
`route_charged` cuts on the angle, which is exact per fragment and needs
no extra argument at any of the eight call sites; `pt_cut_near_beam`
remains as the 275 GeV proton reference and `pt_cut_for(p)` gives the
threshold at any momentum. Measured: ⁶Li α-tag 1.7% (high-acceptance) /
1.3% (high-divergence) at 137.5 GeV/u; ⁷Li α-tag unchanged at 96–98%
(off rigidity, no near-beam cut). `coherent.tag_acceptance_sampled` now
closes on `tag_acceptance_angular` instead of the constant cut.

**S9 — documented:** spectator routing uses the true fragment identity
and momentum; all Z/A = ½ fragments share the beam rigidity; no Z-ID or
misidentification model (plans/06 #19) — the quoted ⁷Li "96–99% into
the Roman Pots" is an acceptance, not a purity.

**S10 — LOW:** the 20260729 smearing loop generates only truth-accepted
bins at their bin centres (outflow only, 4.4% of reconstructed weight
dropped by `reco_mask`, migration underestimated); the reco selection
itself is on reconstructed bin centres (measurable). Superseded by
`RecoResponse` (in-cell positions, loosened generator cuts).

**S11 — LOW (portability):** `money_delta_20260729.py:478` uses
`np.trapezoid` (absent in NumPy < 2.0; the installed 1.23 fails —
`delta_models.py` carries the shim, the script does not);
`_check_reco_mask_invariants.py:86` opens the source without an
encoding (fails on a GBK console; passes S1–S6 with PYTHONUTF8=1); the
EPPS21 ⁶Li grid is not installed; the script aborts by design if the
solved |A_bag| moves > 1% from hard-coded values (so fixing S1 requires
updating those).

**S12 — LOW:** the 20260729 LOW configuration is 5 × 27.5 GeV/u;
`beams.default_configs` gives 20.5 GeV/u (41 × Z/A).

**S13 — LOW (documented convention):** `beams.LI6` effective
polarizations 1/3 per nucleon combined with /A in `fom.py` give
g₁(⁶Li)/g₁(d) = 0.119 per nucleon — the 2-of-6 dilution applied twice
in the A_∥ projections.

**S14 — LOW:** no crossing angle in the fast-sim (the |η| = 3.5 edge
moves by 0.15–0.5 units); the reco layer has it.

**S15 — known:** vacuous assertions in `test_smoke.py:17` and
`test_spectator.py:64` (reported earlier, deliberately left).

*Verified correct (reviewer):* kinematics (4-vector re-derivation to
10⁻¹², s massless to ≤ 2×10⁻³); electron-method resolution by finite
differences; cos 2φ amplitude and its low-y limit −Δ/(2(1 + R)F₁);
error estimators; rates N = σ·Δx·ΔQ²·L with L per nucleon (10 fb⁻¹/u =
1.67 fb⁻¹ per nucleus; N_DIS(18 × 137.5, Q² > 1) = 5.7×10⁹);
sig² = Σa²P²N/2; the moment constraint of `delta_models` (10⁻³);
spectator κ values and boost; far-forward RP/B0/ZDC windows; the
20260729 reco pipeline invariants (reco_mask ⊆ accepted, S1–S6).

## 8. Recommendations (ordered)

0. **Two author decisions before the next production numbers:** (a) the
   b₁-sector sign (G1/S2) — flip `t_geo` and `azz` to the HJM/HERMES
   convention or document b₁ ≡ −b₁^HJM, add a test against Cosyn
   Eq. (27), and obtain JM PLB 223:218 to settle the Δ sign; (b) the
   fast-sim `r1998` (S1) — restrict Θ to the log term, re-run the July
   production (S11 first), correct the notes and the fast-sim README
   headline (S4) and the "folded" P_zz band (S5).
1. **State and model the spin-state-sample acceptance stability (F1):**
   bunch-by-bunch alternation as the baseline of the run plan; a
   per-fill efficiency option in `RecoResponse.expected_counts` /
   `expected_counts_by_fill` and a printed stability requirement in the
   R-scripts; a line in plans/07's systematics table.
2. **Get the low-Σ_h hadronic resolution (F2):** the one external number
   the letter cannot do without; meanwhile quote 5R/7R for the 15–30%
   band and for an absolute-noise variant.
3. **Fix the report's δQ²/Q² and φ′-dilution statements and source or
   label the angular table (F3)** — done in this pass.
4. **η-dependent calorimeter resolution + "best" energy (F4).**
5. **Quote or iterate the K model dependence (F5)**; consider a
   migration-matrix unfolding along x for the letter.
6. **Label 6R as background-free/Z-identified and carry the plans/06
   budget into WP5 (F6).**
7. **Kernel completeness (G2):** add the O(γ²) b₁–b₃ cos 2φ term (Cosyn
   Eq. 17e) so the pseudo-experiments can subtract it with κ.
8. **Fast-sim hygiene (S7, S8, S11–S13):** OMD window 0.45–0.65; angular
   α-tag cut 10σ_θ·A·p_u; `np.trapezoid` shim and file encoding in the
   20260729 script and the invariant checker; one LOW configuration;
   resolve the ⁶Li effective-polarization convention.
9. **Positivity guard in the φ sampler (G3)** and an external-convention
   test for the tensor sector (the identity tests cannot see G1).
   *Both done 2026-08-25* (plans/08 B1, B2), together with one rank-2
   geometry for J = 1 and J = 3/2 (B3: the spin-3/2 branch's rate and
   cos 2φ channels disagreed with each other by 3; latent, since the ⁷Li
   rank-2 slots default to None) and coverage for the spin-temperature
   populations, a non-default b₂, and θ_S between 0 and π/2 (B4).

## 9. What this pass changed (documentation and consistency only)

* `reports/reconstruction_chain_report.template.html` (+ rebuilt
  HTML/PDF): δQ²/Q² statement (F3), φ′-dilution statement (F3), the ratio
  formula and its fit description (F11), Schematic 2 box text, §5.1
  coherent wording (F10), σ_c/n_kept (F13), Table 2 caption (F13), the
  acceptance-stability requirement in §3 step 6 (F1), and the §7.1
  known-limits paragraph (F1, F4, F5, F6).
* `evgen/polligen/reco.py`: docstrings of `emcal_resolution`,
  `tracking_resolution`, `tracking_angular_resolution` label the
  placeholders (F3, F4).
* `evgen/polligen/bookkeeping.py`: `tensor_flip_plan` default
  φ_S = π/2 (F12) and the bunch-by-bunch requirement in its docstring.
* `plans/04_open_questions.md` #3: the quantitative acceptance-stability
  requirement (F1). `plans/07`: systematics row extended. `plans/00`:
  addendum. `evgen/README.md`: test count (G10). `README.md`: pointer.
* No physics code path was changed; all 118 tests pass.
