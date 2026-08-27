# Plan 08 — Completing the simulation chain (2026-08-25)

**Goal.** Close the pieces of the simulation chain that the
reconstruction-chain note (`reports/reconstruction_chain_report`) and the
code review (`docs/code_review_2026-08-25.md`) leave open, so that every
number the PLB letter quotes comes out of code with a test behind it.

**Method of this plan.** Six independent audits (inclusive reco chain,
coherent chain, kernel, hadronic final state, fast simulation, and the
report's own prose) proposed 67 gaps; each was then handed to an
adversarial reviewer instructed to refute it against the code. 44 were
refuted — already implemented, explicitly justified in the documents, or
wrong physics — and 23 survived. This plan is the 23, ordered, with the
corrections the verification pass produced folded in.

**Status legend:** ☐ todo ◐ started ☑ done · **D-items are decisions or
external inputs, not work.**

---

## 8.0 What is actually missing

The measurement layer exists and closes (report §7). What is missing
falls into four groups:

| group | what is missing | why it matters |
|---|---|---|
| **A** | the systematics the estimator *cannot* cancel: a fill-dependent acceptance, an assumed relative luminosity, assumed nuisance harmonics, an energy scale, a wrong-MC correction | the spin-state ratio cancels a *common* acceptance; every non-common effect is unmodelled today, and the largest of them is far bigger than the report states |
| **B** | the analysis's own null tests and its optics scan: sin-harmonic columns, the WP5 curves versus the near-beam cut | the letter's answer to "the p_T cut is undocumented for Li" does not exist as a curve |
| **C** | kernel hygiene and completeness: a positivity guard, an external convention anchor, one rank-2 geometry for both spins | the tensor sector is tested only against itself |
| **D** | the fast simulation's sourced defects: the near-beam cut applied as a momentum cut to a nucleus, the OMD window, R₁₉₉₈ | published acceptance numbers |

Nothing in A–D needs PYTHIA 8, FLUKA or an ePIC-official number. Those
are §8.3.

---

## 8.1 Ordered work

### A1 — per-fill acceptance and assumed luminosity ☑
*Files:* `polligen/reco.py` (`spin_state_ratio` neighbourhood,
`expected_counts_by_fill`), `polligen/recopseudo.py`
(`RecoResponse.expected_counts`, `measure_inclusive`, `CoherentResponse`,
`measure_coherent`), both R-scripts.

The ratio estimator cancels ε(φ′) *because the same ε multiplies every
fill*. Today `phi_eff` is a single callable and `u1, u2` are handed to
the fit as exactly known, so the pseudo-experiments can only exhibit the
case that works. Add:

- a per-fill `phi_eff` (sequence of callables) on both
  `reco.expected_counts_by_fill` and `RecoResponse.expected_counts`;
- `reco.fill_acceptance_bias(eps2_per_fill, pzz, lumi_fractions)` =
  Σ_f l_f (P_f − P̄) e_f / σ_P² — for two equal-luminosity fills this is
  (ε₂⁺ − ε₂⁰)/(P₊ − P₀), the analytic form the code review derived by
  hand;
- `measure_coherent(lumi_assumed=…)`, mirroring `measure_inclusive`;
- `CoherentResponse.with_cut(cut_scale_xy)` — a re-masking *view* of the
  same recoils, so a per-fill cutout perturbation is not confounded with
  MC noise (rebuilding the response would only cancel the common mode to
  MC statistics).

*Correction this produces.* The audit measured the coherent case: with
the slot only half the β bins are live, the t-shape template spans
0.44–0.91 there and is 99% anti-correlated with the constant, so a
**10⁻³** relative change of the vertical envelope between fills biases
a_t by ~19% — not the ~1% that `docs/code_review_2026-08-25.md` §5.2,
report §7.1, `plans/04` #3 and the `plans/07` systematics row all quote
for a **1%** change. The requirement is bunch-by-bunch alternation or
~10⁻⁴ envelope stability. Verify with the code before rewriting those
four places.

*Tests:* exact closure when the fills share an acceptance; the analytic
inclusive bias reproduced to a few %; monotone, near-linear coherent
bias versus the envelope perturbation.

### A2 — sin-harmonic columns: the azimuthal-alignment null test ☑
*Files:* `reco.basis_2d`, `reco.harmonic_ratio_fit_2d`,
`recopseudo.CoherentResponse.expected_counts_2d` / `basis_means`.

The inclusive fit already carries sin 2φ′ (`with_sin=True`) as the
spin-axis calibration. The coherent fit is cos-only, so the coherent
channel has no null test at all. With unpolarized leptons and a headless
alignment axis, σ(α, β) = σ(−α, −β): sin 2α, sin 2β and sin(α+β) are
exactly forbidden, which makes all three a clean diagnostic.

Add, behind `with_sin=False` so nothing existing moves: `"e_s"` =
⟨sin 2α⟩, `"t_s"` = the *t-weighted* ⟨(t/t_ref) sin 2β⟩ (the plain
⟨sin 2β⟩ fails its own closure by 7%, exactly as `c2t` does for the
cosine), `"m_s"` = ⟨sin(α+β)⟩, and the two missing analytic α-integrals
in `expected_counts_2d`.

*Signature separation:* a spin-axis error δ gives (tan 2δ, tan 2δ,
tan δ) on the three ratios; a Roman-Pot azimuthal **roll** gives
(0, tan 2δ_t, tan δ_t). Frame it in the report as the pot-roll null —
the inclusive fit calibrates φ_S about an order of magnitude better.

### A3 — u₁, u₂ as a fit argument ☑
`measure_coherent` passes the same `u1, u2` to the generator and to the
fit. Add `u_coeffs_assumed`; the fit side already takes an independent
value. Quote the band from **u₂**, not u₁: the leakage into a_e is
first order, δa_e ≈ a_t·Δu₂·⟨cos 2β⟩, and with the slot's ⟨cos 2β⟩ =
+0.77 and a_t/a_e ≈ 12–34 a ZEUS-1σ Δu₂ = ±0.024 moves a_e by ~20% —
several times the 10-year error. Report §4's "u feeds the pure harmonics
only at second order" is wrong for u₂ and must be corrected.

The honest systematic is the in-situ one: fit (u₁, u₂) from the
spin-averaged counts of the same data and propagate; ZEUS is then the
prior it is.

### A4 — WP5: the optics scan ☑
*New:* `evgen/scripts/coherent_optics_scan.py`.

Take acceptance and N_tag **analytically** from `reco.rp_hole_acceptance`
(exact for any cutout, no MC) over p_T^cut ∈ [0.1, 0.7] GeV, with
exp(−B t_min) and the 0.73 rate weighting folded into the central curve
and the IR-8 published efficiencies overlaid. For the δa_t axis,
importance-sample the recoil above the cutout inside `CoherentResponse`
— the present sampler leaves *zero* accepted events above 0.6 GeV, so
configuration 2 cannot be marked on any MC curve today.

*Fixes an inconsistency:* `reco_chain_figures.py` panel (c) plots the
circular acceptance exp(−B c²) while the chain uses the slot — a factor
≈ 2.9 on the tagged yield, which also feeds `plans/04` #20.

*Done 2026-08-25* — panel (a) of the new script plots the circular,
square and slot curves together so the difference is visible rather than
implicit.  Measured (slot, B = 50): tagged fraction
32% / 3.0% / 4×10⁻⁵ / 2×10⁻⁷ and δa_t/a_t = 1.6% / 5.7% / 104% / 540%
at 0.10 / 0.22 / 0.45 / 0.60 GeV.  `t_floor` must stay well below the
cut: it is a hard lower bound on the TRUE |t|, so with divergence
smearing it biases the acceptance downward if pushed too close (the
script uses 0.25 cut²).

### A5 — energy scales, wrong-MC closure, HFS beam guard ☑
- `RecoModel.e_scale` (applied to the smeared E′) and an η-slope on
  ε_eID. The **electron** scale is the bigger lever:
  d ln x_mixed/d ln E′ = 2 − y ≈ 2, twice the hadronic one. A flat
  ε_eID scale is identically null (it cancels in the ratio) — do not add
  that row.
- One post-calibration multiplier in `HFSResponse.hadronic`, not in
  `HadronResponse` (whose parameters are reused to build the quiet
  library, so a scale there would be baked into f_Σ and double-counted).
- Wrong-MC systematic: generate at one `y_had_res`, correct with
  another. No new argument is needed — `delta_from_amplitude` already
  takes `summary` freely. Use **common random numbers** across the two
  responses, or the effect is indistinguishable from seed noise.
  Measured: generate at 0.30, correct with 0.25 → Δ̂ shifts 0.4–1.5%,
  against 10-year statistical errors of 0.4–1.0%.
- `HFSSample.concatenate` and `HFSResponse` must refuse to merge or
  transfer a response across beam energies.

*Done 2026-08-25.*  Measured at the four sweet spots with `--syst-scan`
(common random numbers; MC noise floor 0.13–0.21%): electron scale ±1% →
Δ̂ shifts 0.2–1.4%; hadronic resolution mismatch (generate 0.30, correct
0.25) → 0.5–1.3%; ε_eID η tilt 0.05/unit → 0.02% (a smooth efficiency
shape is almost exactly null, a flat one exactly null — so the audit's
advice to drop the flat row was right and the η row is nearly null too);
Yellow Report EMCal table → 0.000% at spots 1–3, 0.5% at spot 4.  All
four sit below the 3–11% model dependence of K.

One correction to the audit: the electron lever is 2 − y only for the Σ
method.  The Gaussian y stand-in smears the *true* y and never sees E′,
so it gives exactly 1 and understates the electron energy scale by about
two.  Both levers are pinned by tests.

### A6 — forward-folded response replaces the bin-by-bin K ☑
`RecoResponse.fold(delta_func, mask)` = the response-weighted
⟨(∂A/∂Δ)·Δ(x_true, Q²_true)⟩, then either two K iterations with the
shape refitted to the data, or a 2–3-parameter Δ(x) shape fitted
**folded through the response** per Q² slice.

*Two design constraints from the audit.* (i) This is not a yield
migration matrix: Â is a weighted *mean* of per-event amplitudes and
∂A/∂Δ changes by ~2× between adjacent x bins, so a w·eff matrix
mis-weights every off-diagonal element. (ii) Do not invert it: with
C_ii/row-sum ≈ 0.5 the matrix is ill-conditioned and plain least squares
on noise-free pseudo-data oscillates to ±30% — worse than the
+1.6…+27% bias it is meant to cure. Unfold along x only; include
feed-in bins beyond the plotted range; propagate the response MC error
(0.5–1.2% per bin, currently unpropagated).

*Done 2026-08-26.*  The amplitude is **exactly linear** in Δ at fixed
kinematics — Δ enters the master formula in one place — so
`RecoResponse.fold(delta_func, mask, category)` is the response as an
exact linear operator, not a linearisation, and it reproduces
`bin_summary["a_reco_bin"]` to 3×10⁻¹⁶ when handed the model the response
was built from (which is what makes it a drop-in for K).
`fold_shape_fit` then fits a 2-parameter tilt of a prior shape *through*
the response per Q² slice, on a bounded grid with the normalisation
profiled out — bounded because that, not a matrix inverse, is the
conditioning guard the audit asks for.  Measured on money plot 7R with
the wrong prior: the residual bias over the plotted bins falls from
22% / 24% / 99% (bin-by-bin K) to at most 3.5% / 7.1% / 6.9% (moment_B
prior) and 4.7% / 5.5% / 5.0% (toy prior), and at the four sweet spots
the model dependence K(moment_A)/K(prior) − 1 falls from (−5.0, +8.5,
−9.3, +6.3)% and (+9.7, +2.5, +10.6, +2.8)% to (−0.9, −1.2, −0.2, +0.3)%
and (+4.7, +2.8, +2.2, −0.9)%.  The 7R error bar now carries four terms, not one:
statistics, the shape fit, the response MC (`fold_mc_error`, a stratified
within-cell estimator validated against eight independent response
seeds — 0.10–1.00% per plotted bin, median 0.29%, and
0.32 / 0.27 / 0.21 / 0.11% at the four sweet spots at
`n_mc_per_cell = 400`) and, dominating them, the spread over the priors
the tilt family cannot absorb.  That last one exists because the
adversarial pass caught the bar advertising three sources while being
blind to the biggest: it was as small as 0.070× the residual it did not
cover.  Refitting from every other registry shape bounds it, and the
bound is measured to hold — ≥ 1.037× the residual at every plotted bin
of every (slice, prior) pair.  Both paths stay: `--unfold model` is the
default and reproduces every published 7R number bit for bit.

### A7 — η-dependent calorimeter resolution (F4) ☑
`emcal_resolution(..., eta=None)` with the Yellow Report table, `eta=None`
reproducing today. Keep the backward constant at the documented 1% ePIC
PbWO₄ spec. Do **not** deprecate `energy="best"`: at η = −1.6 tracking
gives 1.1% against a YR EMCal 3.0%, so `min(cal, trk)` *is* the switch
recommendation 4 endorses — make it the default. No published bin moves;
the higher-value neighbour is the *angular* placeholder table (F3),
which dominates δQ²/Q² at every sweet spot and is D6.

### B1 — φ-sampler positivity guard (G3) ☑
`InclusiveSampler._amplitudes` guards only 1 + w_avg ≤ 0, so wherever the
φ density dips negative the accept–reject silently samples max(W, 0) and
`weights_for` does not clip at all. Guard at that chokepoint with the
**exact** minimum over φ (for W/(1+w) = 2Bc² + Ac + (1 − B) the minimum
is at c = −A/(4B) when B > 0 and |A| ≤ 4|B|, else at c = ±1), and extend
it to `phi_histogram_pseudo`, which bypasses `_amplitudes` and is the
path behind the headline figures. Margins in every existing
configuration are ≥ 0.97, so nothing existing fires.

### B2 — an external anchor for the tensor sign (G1 gate) ☑
Every tensor assertion in the suite compares code to code: patching the
sign in both `xsec.py` and `asymmetries.py` leaves all tests green. Add
`TENSOR_LL_SIGN` as one module constant and one *external* test against
Cosyn Eq. (27), A_zz = −(2/3) b₁/F₁, in the F₂/R-free form

    A_T(θ_S = 0) · (1 + ε(y) R(x, Q²)) == SIGN · (−2/3) b₁/F₁

(the widely quoted form that also divides by [1 + 2(1−y)/y²] and
multiplies by [1 + 2(1−y)(1+R)/y²] double-counts R — those two brackets
are identically 1 + εR given F₁ = F₂/(2x(1+R)) — and misses by 1.17).
Land it green at the current value with a comment that Cosyn requires
the opposite, so **D1 becomes a one-line change**.

*Done 2026-08-25.*  Verified non-blind: with `TENSOR_LL_SIGN = -1` the
other 125 tests still pass and only `test_tensor_convention.py` fails,
with the decision in the assertion message.

### B3 — one rank-2 geometry for J = 1 and J = 3/2 ☑
`_tensor_moments` returns (Q_NN, Q_NN); Cosyn Eq. (9) gives
T_LL = Q_NN P₂(cos Θ) and T_TT = (3/2) Q_NN sin²Θ, so the cos 2φ channel
wants 3 Q_NN. Spin-1 already satisfies this; the spin-3/2 branch is off
by 3 between its own rate and cos 2φ channels. Latent (the rank-2 ⁷Li
slots default to None) — which is why it needs a test, not a comment.
Anchor the test on the density matrix built from `spin.py`'s own
operators, not on one branch against the other. For J = 3/2 write a
*characterization* test citing `plans/04` #14: asserting cross-J equality
would freeze a convention the documents declare open.

### B4 — pin the untested paths ☑
`populations_maxent` (spin-temperature populations), `b2_func` ≠ 2xb₁,
the θ_S ≠ π/2 thirds identity, and `helicity_flip_plan(...).pzz_true` —
the quantity the money scripts divide by — have no test.

*Done 2026-08-25* for the first three: the spin-temperature populations
are pinned for both spins against P_zz = 2 − √(4 − 3p_z²) and the
geometric-ratio property; a `b2_func` override is pinned against `azz`
with an explicit "not the default" guard; and the tensor rate is checked
across θ_S including the magic angle, where it vanishes exactly.  The
spin-3/2 rank-2 sector is a characterization test (plans/04 #14).

*Done 2026-08-26* for the fourth: `evgen/tests/test_bookkeeping.py` pins
`helicity_flip_plan(...).pzz_true` for every branch — exactly 0 for
j = 1/2, the spin-temperature moment for `pzz=None` (against closed forms
of the geometric population ladder with exact rational anchors at t = 3,
and against the density matrix built from `spin.py`'s own operators),
closure for an explicit `pzz`, and evenness under pz → −pz.

### C1 — far forward: angular near-beam cut, OMD window, nucleus mass ☑
`route_charged` already receives θ, so the 10σ envelope becomes one line
plus an `Optics` field — better than threading (A_frag, p_u) through
eight call sites, and exact per event. 0.20 GeV on a 4 × 137.5 GeV α is
5σ, not 10σ; the ⁶Li α tag falls 6.26% → 1.66%, which is the 1–3% the
README already quotes from the code review. OMD window → (0.45, 0.65)
per Jentsch–Tu–Weiss Table I (the symbol there is ζ, not ξ). Add a
physical `NUCLEUS_MASS` table beside `MASSES`; A·M_U shifts the rigidity
ratio and the off-rigidity tail.

*Done 2026-08-25* except the nucleus mass. Measured: ⁶Li α-tag 1.7% /
1.3% (high-acceptance / high-divergence) at 137.5 GeV/u, 13.3% at
50 GeV/u; ⁷Li α-tag unchanged at 96–98%. One discrepancy is left
deliberately: the high-divergence envelope has no documented spec —
plans/06 §6.5 derives ≈ 0.41 GeV at 275 GeV and the fast sim has always
quoted the rounded-up 0.45, which is where every published fast-sim
number and the money-plot report come from, while `polligen.reco` quotes
0.41. Both are now cross-referenced in the source; unifying them changes
the published coherent high-divergence acceptance by 5× (4×10⁻⁵ →
2.2×10⁻⁴) and is a documentation decision for the author, not a code
one.

*Done 2026-08-26* for the nucleus mass.  `spectator.NUCLEUS_MASS` carries
the ground-state **nuclear** masses (AME2020 atomic masses less the
electrons; the recipe reproduces the six CODATA light nuclear masses, and
the separation energies of the channel table close against it), and the
beam boost and the two-body reduced mass in `kappa` both use it.  A·M_U
is not the mass of a nucleus: it is low by 12.6 MeV for ⁶Li, and the
boost γ, γβ, p_lab, θ and R of every spectator inherited that.  The
kinematics move by 2×10⁻³; the *acceptance* moves more, because the R
window edges are hard — the ⁶Li α sits at R = 0.99813, not 1.00000, and
sliding it 0.0022 across the RP edge at R = 0.95 takes the α tag from
1.65% to **1.85%** (high acceptance) and 1.31% to **1.51%** (high
divergence) at 137.5 GeV/u, and 13.2% to 13.5% at 50 GeV/u.  ⁷Li, off
rigidity at R = 0.85571 with no nearby edge, moves +0.13% relative.
`tagged.boost_spectator` was carrying the old A·M_U and is now the same
mass — the same beam had two masses in two modules until this run.

### C2 — R₁₉₉₈ (S1) and one R hook ☑
Θ must multiply only the 0.0485/ln(Q²/0.04) term; the present form
clips R to 1.000 at x ≲ 0.1. Thread `r_func` through `NuclearF2`,
`dsigma_dx_dq2`, `depolarization_d`, `InclusiveKernel` **and**
`polarized.ToyG1._f1` — the fourth consumer the dated scripts'
monkey-patch misses. Keep the monkey-patch in the dated scripts: they
are frozen reproductions of dated notes.

R is the only backend systematic that moves the physics result: F₂
cancels **exactly** in the cos 2φ amplitude (scenario Δ is defined as
(Δ/F₁)·F₁ and D_φ/F₁ depends only on x, y, R), while Δ/F₁ = −2(1+R)Â
carries R directly — toy R = 0.157–0.171 across the acceptance against
a corrected R₁₉₉₈ = 0.12–0.34. Report the F₂ and R drifts separately.

### C3 — make the July production runnable (S11) ☑
`np.trapezoid` shim (it already exists in `delta_models`), `encoding=
"utf-8"`, `Path(__file__).with_name(...)`; and the two blockers that
fire first: a `UnicodeEncodeError` on "fb⁻¹" under cp936 before the
trapezoid line is ever reached, and the EPPS21 grid preflight (the grid
installs cleanly). Keep the `REFERENCE_ABS_A_BAG` hard guard — demoting
a designed hard stop is the opposite of what the review asks — and add
`--emit-a-bag-reference` so the R fix can paste the new triple.

---

## 8.2 Convention items reserved for the author

These are real defects, but each changes a documented convention rather
than fixing an unambiguous error. They are stated here with the one-line
fix and left for the author.

| # | item | state |
|---|---|---|
| **D1** | **b₁-sector sign.** The code is opposite to Cosyn Eq. (27) / HERMES. B2 lands the anchor and the single constant, so the flip is one line. What flips is the sign of κ — and therefore the O(γ²) subtraction, which must not be coded before the decision. `money_b1.py` plots \|A_zz\|; `money_tagged_azz` is dominated by the O(1) S/D interference. | open |
| **D7** | **⁶Li effective polarization.** `beams.LI6` uses the Cloët 1/3 *per-nucleon* convention in a slot that `ToyG1.g1_nucleus` treats as whole-nucleus and callers then divide by A — the 2-of-6 dilution applied twice (code review S13; g₁(⁶Li)/g₁(d) = 0.119 per nucleon). The structural fix (carry Z and N in `g1_nucleus`, as `NuclearF2.f2a` already does) is unambiguous and also restores ³He's missing ×2 on the proton term; the *value* — 1/3 versus the 0.81 cluster picture — is `plans/04` #6. No headline number uses g₁(⁶Li). | structural fix safe; value open |
| **D9** | **b₁ money plot dilution (S6).** `money_b1.py` draws the signal with the ×0.87 deuteron transfer and no 2/6 against errors with neither — a factor 3 between signal and error. The correct rank-2 transfer is `TaggedModel.tensor_dilution()`, not the vector 0.87. Blocked with D7 because they are the same convention question. | open |

## 8.3 Blocked on external input

| # | item | blocker |
|---|---|---|
| D2 | exact finite-γ kernel (Cosyn Eqs. 9/10/14/24), b₃/b₄ slots, Eqs. 17d/17e | gated on D1 for the subtraction sign. **Correct now, without code:** report §2, code review G2 and the `plans/07` systematics row all quote the leakage as γ²b₁/6 and are ≈ 7× low; the exact ratio a₂(full)/a₂(17e alone) ≈ 6.9. Impact on everything published is ≤ 0.15%; the exposure is the Δ/F₁ ~ 10⁻³ scenario and x ≳ 0.2 at Q² ≈ 1 |
| D3 | radiative corrections (WP4) | **two of the four planned deliverables are void**: collinear ISR generates exactly zero fake cos φ′/cos 2φ′ (the covariant azimuth is invariant under k → (1−z)k to 2×10⁻¹⁶ rad in the massless-target limit), and the ratio cancellation is already demonstrated. What survives is (a) the migration bound on purity/efficiency/K — the number the ≤5% gate should apply to — and (b) the method comparison: y_Σ and y_DA use no beam energy and are ISR-robust, but the chain's own x = Q²_e/(s y_Σ) is biased by exactly 1/(1−z), so report §3's robustness claim is true for y and false for the x the analysis uses |
| D4 | PYTHIA 8 HFS samples | ☑ **closed 2026-08-26** — it was never external. PYTHIA 8.311 builds its own Python bindings against the analysis machine's interpreter (`tools/pythia8/README.md`); the eic-shell container has the C++ library but no bindings, which is what had made this look like a container problem. 12 M events over the three beam configurations now stand in `evgen/samples/`, and A5's beam guard makes the merge safe |
| D5 | incoherent breakup shapes, veto efficiencies, event-level Z-ID | FLUKA licence → BeAGLE. `plans/07` already rules this non-blocking; the m-state-blind dilution argument is algebra, not a simulation gap. *2026-08-26:* the no-FLUKA half is done — the official BeAGLE e+d sample streams over xrootd and the control study is run (`tools/beagle/README.md`), which is what calibrates the cluster model's p_T tail; what FLUKA still gates is A = 6, 7 breakup itself |
| D6 | ePIC numbers: the calorimeter noise/threshold floor at Σ_h ≈ 0.2–0.5 GeV (#21, the one number the letter cannot do without); the **backward-disk angular resolution** (F3 — worth more than A7); Li ring σ_θ (#20); the RP slot geometry; EICROC Z-ID (#19) | external — except the RP slot geometry, which is now measured rather than assumed: an intact ⁶Li through the ePIC far-forward geometry puts the aperture at \|θ_x\| ≳ 2.0 / 1.35 / 1.03 mrad in the 5×41 / 10×100 / 18×275 optics against \|θ_y\| ≳ 1.8–3 mrad, i.e. **open horizontally**, the opposite aspect to `rp_measure`'s slot (`tools/fullsim/README.md`, plans/04 #20) |
| D8 | a coherent diffractive model for ⁶Li (#18) | blocks a coherent electron side and x_P/M_X binning. Pairing the inclusive map with the log-uniform x_P now in `CoherentResponse` gives β = x/x_P > 1 for 58% of events, so this cannot be faked locally |

---

## 8.4 What development run 8 found beyond the plan (2026-08-26)

Moving the work to the Linux box turned three of the "external" blockers
into measurements.  Each is written up where it belongs; this is the
index, and the numbers are in `plans/00` run 8.

1. **PYTHIA 8 was never external** (D4 above).  Building it also found a
   defect in the settings this repository had documented for two months:
   `PhaseSpace:Q2Min` is applied **only** when `Q2Min ≥
   pTHatMinDiverge²`, and `pTHatMinDiverge` defaults to 1 GeV, so the
   requested 0.7 GeV² was silently ignored and the 0.7–1.0 GeV² band —
   31% of the sample, and precisely the band the loosened generator
   window exists to populate — was missing.  `tools/pythia8`.
2. **The hadronic-resolution band is measured, and the toy was
   optimistic**: δy/y = 0.32 / 0.22 / 0.28 / 0.11 at the four sweet spots
   against the toy's 0.28 / 0.17 / 0.24 / 0.07, and the ordering across
   beam energies (0.21–0.10 low, 0.74–0.18 top) settles the energy
   strategy for the x ≈ 0.1 bins, which the reconstruction note lists as
   open.  Purity at the sweet spots falls from 64–68% (stand-in) to
   40–73%, so the response is a bigger correction than the note implies —
   which is the argument for A6.
3. **Every script in the repository runs** (2026-08-26).  Writing the
   reproduction manual meant executing all 30 of them, which found the
   two survivors of code review S11 outside the 20260729 script that C3
   covered: the 2026-07-20 and -07-21 scripts need a fourth PDF grid
   (`nNNPDF30_nlo_as_0118_A6_Z3`, 161 MB) and exit 2 with an install hint
   without it, and `money_delta_20260724.py` still called `np.trapezoid`
   bare.  Both fixed; the dated line is now reproducible end to end.
4. **The Roman-Pot aperture is measured** (D6 above), by shooting an
   intact ⁶Li through the ePIC geometry.  npsim cannot do that with a
   gun — DD4hep has no generic-ion particle type — so the scan feeds one
   through as HepMC (`tools/fullsim/ion_gun_hepmc.py`).  It inverts the
   assumed cutout aspect and therefore the SIGN of the acceptance-induced
   ⟨cos 2β⟩.  **Carried through the chain the same day**
   (`reco.RP_APERTURE_MEASURED`, `rp_measure(cut_theta_xy=…)`,
   `money_cos2phi_coherent_reco.py --rp-aperture measured`): the envelope
   and the aperture are separate constraints and the cutout is the larger
   of the two per axis.  At the LOW configuration the measurement
   survives — acceptance 37.7% → 1.42%, fake ⟨cos 2β⟩ +0.426 → −0.772,
   two |t| bins instead of four, δa_t worse by 6–34× and a_e still
   recovered — and at MID and TOP the aperture leaves nothing in the
   binned window.  Two defects surfaced on the way: a cutout tight enough
   to empty part of the circle made `harmonic_ratio_fit_2d` raise a bare
   "Singular matrix" (it now names the cause and how many bins carry
   weight), and one dead |t| bin aborted the whole figure (it is now
   reported and skipped).  **What is left:** a harmonic basis that works
   under a strongly anisotropic acceptance — fewer columns, wider bins,
   or |t| re-binned inside the window the cutout leaves.  The WP5 optics
   scan now carries the measured edge as a marked line on all three of
   its curves: it is 2.8× / 1.9× / 1.4× the 10σ envelope at
   5×41 / 10×100 / 18×275, so the envelope is never the binding
   constraint.
5. **The e+d control calibrates the cluster tail** — and says no β in
   the two-parameter Hulthén reproduces BeAGLE's shape (plans/02 step
   1.5.3).  Since the ⁶Li α tag is entirely a p_T-tail measurement, its
   model uncertainty is one-sided **upward**.

---

## 8.5 Commit sequence

Each commit is self-contained and leaves both suites green (baseline
125 = 101 evgen + 24 fastsim).

```
A1  reco: per-fill phi' acceptance and the analytic spin-state bias
A1b coherent: per-fill cutout and assumed luminosity in measure_coherent
A3  coherent: assumed u1, u2 as a fit argument
A2  coherent: sin-harmonic columns as an azimuthal-alignment null test
B1  kernel: exact positivity guard on the phi density
B2  kernel: external anchor for the tensor sign behind one constant
B3  kernel: one rank-2 geometry for J = 1 and J = 3/2
B4  tests: pin the untested spin-population and b2 paths
C1  farforward: angular near-beam cut, OMD window, physical nucleus mass
A7  reco: eta-dependent calorimeter resolution
A5  reco/hfs: energy-scale nuisances, wrong-MC closure, beam guard
A4  coherent: WP5 optics scan
C2  structure: correct r1998 behind one r_func hook
C3  scripts: make the July production runnable
A6  reco: forward-folded amplitude response replaces the K factor
```

Commits that change published numbers land with their report/plan edits:
C1 (α-tag acceptance), A1b/A3/A2 (report §7, Table 3, `plans/04` #3 and
the `plans/07` systematics row), A4 (`plans/04` #20 and report §7),
C2 (Δ/F₁), A6 (money plot 7R error bars).
