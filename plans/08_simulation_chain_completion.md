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
prior it is.  **Done 2026-08-28** (`reco.unpolarized_insitu_fit_2d`,
`measure_coherent(u_coeffs_assumed="in-situ")`, `--u-in-situ`).  The
identifiability had to be faced first: with the per-bin acceptance free —
which is what makes both estimators acceptance-free — u is *not*
identifiable at all, because u_b multiplies ε_b in exactly the same way.
So the in-situ measurement necessarily uses the acceptance MC, and the
implementation says so: it fits the spin-averaged counts against the
response's own ε^MC shape (u = 0, P_zz = 0) with the normalisation free,
by Poisson likelihood, alternating with the harmonic fit (five rounds
reach machine precision on exact counts) and propagating cov(u₁, u₂) into
the harmonic covariance through the numerical Jacobian dA/du.  Measured at
the tagging optics, one year, per |t| bin (`money_cos2phi_coherent_reco.py
--config 0 --optics tagging --exact --u-in-situ --n-mc 6000000`, which
now prints the pair and its error for every |t| bin): δu₂ = 0.0028 /
0.0040 / 0.0070 / 0.0136 at 5 × 40.8 and 0.0016 / 0.0026 / 0.0049 /
0.0111 at 18 × 137.5 (`--config 2`), i.e. **1.8–15× better than the ZEUS
1σ band of 0.024** that Table 6's systematic assumes — 15× where the bin
is full, and only 1.8–2.2× in the sparsest bin of either configuration,
which is where the leakage matters most.  The leakage into a_e falls with
it all the same, from the +0.00066–0.00093 (6.6–9.3%) of an assumed wrong
u₂ to a propagated 0.00005–0.00063 (0.5–6.3% of a_e) across the eight
bins, which adds nothing measurable to the statistical error in
quadrature.  What is left is the acceptance-shape uncertainty this
transfers the problem to.

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
| D3 | radiative corrections (WP4) | ☑ **closed 2026-08-28 — it was never external.** Two of the four planned deliverables were void by algebra and the other two are now measured, in `evgen/polligen/radiative.py` (exponentiated leading-log D(z, Q²), Kuraev–Fadin / Nicrosini–Trentadue; ∫D = 1 + O(t²) with a 7×10⁻⁴ residual at t = 0.070), the default-off `recopseudo.RecoResponse(isr=…)` hook and `money_cos2phi_reco.py --isr`. **Void:** collinear ISR generates exactly zero fake cos φ′/cos 2φ′ — the covariant azimuth is invariant under k → (1−z)k in the massless-target limit, now *pinned* at 3.6×10⁻¹⁵ rad over the 2×10⁴-event flat sample of `test_covariant_azimuth_is_invariant_under_a_collinear_photon` (z ≤ 0.9; the 2×10⁻¹⁶ quoted here before was a gentler sample). Over the 1.84×10⁶ events of the response, where the physical ⁶Li mass leaves the O(γ²) residual, max |Δφ′| = 2.6×10⁻² rad and the fake cos 2φ′ is 9×10⁻⁸ rate-weighted (`RecoResponse.isr_dphi`) — the two samples are not the same and their residuals differ by two orders — and the ratio cancellation is demonstrated. **Measured:** (a) the migration bound on purity/efficiency/K with common random numbers at the four mid-configuration sweet spots — purity 0.653→0.638, 0.633→0.613, 0.679→0.659, 0.684→0.640, efficiency 0.414→0.404, 0.590→0.572, 0.374→0.369, 0.653→0.634, and Δ̂ corrected with an ISR-free K biased by +0.62±0.03/+0.50±0.02/+0.94±0.03/+1.22±0.02% (mean ± sem over eight response seeds; one draw scatters by 4–14% of the bound, seed-to-seed sd 0.05–0.09 points, so a single-seed number is not publishable), rising to a 1.8–2.8% band as the generator window is opened to Q² ≥ 0.15–0.02 GeV², where the truncated feed-in saturates — ≤2.9% is what the 5% gate is read against, and it passes; a HERA-style E − p_z window brings it to ≤0.25% while keeping 87% of the non-radiative rate; (b) the method comparison at z = 0.092, observed/hard for (Q², y, x), at the rate-weighted ⟨y⟩ = 0.189 of the whole selected sample: electron (1.102, 1.351, 0.740), Σ (1.000, 1.000, 0.908), JB (0.976, 0.908, 0.976), DA (1.214, 1.000, 1.102), mixed (1.102, 1.000, **1.000**). The electron rows go as (y + z)(1 − z)/y and are therefore far worse at the sweet spots themselves (y = 0.010–0.026), where y is off by 4.2–9.2 and x by 0.109–0.239; the mixed row is 1.000 at every y, which is the whole argument for the chain's choice. The 2026-08-28 correction stands and is now pinned against a four-vector construction: x = Q²_e/(s y_Σ) is exact under a collinear photon and it is the Q²_e *label* that migrates by 1/(1−z). Q²_Σ = p_T,e²/(1 − y_Σ) uses no beam energy either, so an e-Σ label would carry no migration at all — a chain change, not made. Still external / uncalculated: the TENSOR-sector RC (plans/05 §5.5) |
| D4 | PYTHIA 8 HFS samples | ☑ **closed 2026-08-26** — it was never external. PYTHIA 8.311 builds its own Python bindings against the analysis machine's interpreter (`tools/pythia8/README.md`); the eic-shell container has the C++ library but no bindings, which is what had made this look like a container problem. 8 M events over the three beam configurations now stand in `evgen/samples/`, and A5's beam guard makes the merge safe |
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
   optimistic**: δy/y = 0.55 / 0.28 / 0.50 / 0.15 at the four sweet spots
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
   reported and skipped).  ~~**What is left:** a harmonic basis that works
   under a strongly anisotropic acceptance — fewer columns, wider bins,
   or |t| re-binned inside the window the cutout leaves.~~
   **SUPERSEDED 2026-08-28** (verified numerically; the disposition is
   below).  The WP5 optics
   scan now carries the measured edge as a marked line on all three of
   its curves.  Priced against the per-configuration Yellow Report
   envelopes it is 0.91× / 0.75× / 1.12× the horizontal 10σ half-width at
   5×41 / 10×100 / 18×275, so the envelope binds at the two lower
   configurations and the silicon only marginally at the top (the
   "2.8× / 1.9× / 1.4×, never binding" quoted here until 2026-08-28 was
   the same edge against the retired single 73 µrad;
   `tools/fullsim/README.md`).

   *Why item 4's "what is left" is superseded.*  The problem it names does
   not exist at the configurations the programme now quotes, and where the
   anisotropy is extreme it has become an **empty-sample** problem rather
   than an ill-conditioned-fit one.  At the tagging optics the acceptance
   *is* strongly anisotropic — ⟨cos 2β⟩ of the tagged sample runs −0.50 to
   +0.14 across the four |t| bins, and four of twenty-four β bins are hard
   empty in the lowest ones — and the design is nonetheless near
   orthogonal: condition number 1.80–2.86, rank 7/7, corr(a_e, a_t) <
   0.005 (worst 0.0046) and no non-constant parameter pair correlated
   above 0.031, at all three configurations, and the fit closes on the
   injected coefficients to every printed digit on exact counts even where
   those four β bins are empty (an empty bin already carries weight zero).
   *Measured on* `money_cos2phi_coherent_reco.py --config {0,1,2}
   --optics tagging --exact --n-mc 6000000`: the condition number is
   σ_max/σ_min of the weighted design the fit actually solves (`basis_2d`
   with the response's `beta_means`, scaled by 1/√var_t — the array
   `_harmonic_rank_guard` receives), the correlations come from the
   returned `cov`, and ⟨cos 2β⟩ is the tagged sample's mean over the true
   β of the recoils reconstructed into the bin.  The item was written on
   the pre-2026-08-27 rigidity-scaled energies, where 2.0 mrad of silicon
   sat at |t| ≈ 0.06
   GeV², inside the fitted window; at the γ-matched energies
   p_ion = 244.8 GeV, so 2.0 mrad is |t| > 0.24 GeV² and the
   measured-aperture case has
   **zero** accepted recoils, not an awkward fit.  Of the three levers:
   *fewer columns* buys nothing, because the three sin columns are already
   orthogonal to the cos ones — dropping them (`--no-sin`) leaves the
   condition number unchanged to every printed digit, 2.33 / 1.94 / 1.80 /
   2.01 in the four bins at 5 × 40.8 either way; *wider bins* buys nothing
   either (2.33 → 2.40 → 2.46 in the sparsest bin going 12 × 24 → 8 × 16 →
   6 × 12, i.e. slightly worse) and hits the rank guard at 4 × 8, where
   ⟨cos 2α⟩ vanishes identically — it belongs to A12, which closed the
   low-count problem with a likelihood instead; and *|t| re-binned inside
   the window the cutout leaves* is the only surviving content — and it is
   an unclaimed gain, not a repair.  **Follow-up, measurable but not
   adopted:** the fixed 0.05–0.25 GeV² window discards 67–74% of the
   tagged sample below 0.05 GeV² (reco |t| quantiles 0.010 / 0.018 /
   0.031 / 0.052 / 0.101 at 5 × 40.8); four added bins over 0.006–0.05
   fit cleanly (rank 7/7, exact-count closure exact) with δa_e = 0.0022 /
   0.0022 / 0.0026 / 0.0031 at edges 0.006 / 0.017 / 0.028 / 0.039 / 0.05,
   halving the combined one-year δa_e from 0.00205 to 0.00105 — the script
   takes `--t-edges` and prints the combination, so
   `money_cos2phi_coherent_reco.py --config 0 --optics tagging --exact
   --n-mc 6000000 --t-edges 0.006,0.017,0.028,0.039,0.05,0.08,0.12,0.17,0.25`
   is the whole measurement.  Those bins are
   resolution-dominated (δp_T,y = 93 MeV; the 0.17–0.25 reco bin already
   has t_ref = 0.0895), so they buy statistics for the flat a_e far more
   than shape information for a_t(t), and t_min = 3 × 10⁻³ GeV² at
   x_P = 0.01 bounds the window below.
5. **The e+d control calibrates the cluster tail** — and says no β in
   the two-parameter Hulthén reproduces BeAGLE's shape (plans/02 step
   1.5.3).  Since the ⁶Li α tag is entirely a p_T-tail measurement, its
   model uncertainty is one-sided **upward**.

---

## 8.6 The 2026-08-28 review (development run 11)

A second adversarial pass over the chain, this time with the questions
"is each step right", "does any analysis step read what an experiment
cannot", "is every assumption stated".  Findings and their disposition are
in `docs/code_review_2026-08-28.md`; the repository-level summary is
plans/00 run 11.  Items that change the plan above:

| # | item | state |
|---|---|---|
| A5b | hadronic-scale calibration keyed on the true cell → reconstructed-bin map (`HFSLibrary.cell_means_reco`) | ☑ |
| A8 | PYTHIA library truncated at x = 16/s (`mHatMin`) → regenerated | ☑ |
| A9 | `spin_state_ratio` variance for single-state bins; `_ratio_to_modulation` Jacobian | ☑ |
| A10 | hadron acceptance in the lab frame; target-mass term; σ-weighted p + n merge | ☑ |
| A11 | 6R at the tagging optics with ensembles; per-fill perturbation on the binding cut; `--exact` systematics | ☑ |
| A12 | low-count bias of the bin-wise 2-D ratio below ~30 counts per (α, β) bin: a likelihood fit, or adaptive binning | ☑ **closed 2026-08-28**, see below |
| D3 | the ISR statement corrected (above), then **closed in code** 2026-08-28: `polligen/radiative.py`, the default-off `RecoResponse(isr=…)` hook, `money_cos2phi_reco.py --isr`, 25 tests; bound +0.5 to +1.2% of Δ̂ in the published generator window (≤2.9% with the low-Q² feed-in opened up, ≤0.25% behind an E − p_z window) against the 5% gate | ☑ |

### A12 — what the low-count bias was, and what closed it

The defect is the *inversion*, not the ratio.  R_b itself is unbiased at
any count (conditional on the bin total the per-fill split is multinomial
and R_b is linear in it), but T̂ = R(1 + u)/(σ_P² − P̄R) is strongly curved
over the range |R| ≤ max_f |P_f − P̄| that R actually explores when the bin
holds a handful of counts — a bin with one count has R = ±0.9 exactly.  Its
second derivative 2(1 + u)σ_P²P̄/(σ_P² − P̄R)³ has the sign of P̄, so the
flip plan's P₀ = −2P₊ (P̄ = −0.3) makes it strictly *concave*, each bin
carries a *negative* Jensen offset (1 + u_b)P̄/(σ_P² ν_b), and the
data-driven 1/var weights add an opposite-sign term of the same order.
ε_b is flat in α and modulated ×25 in β by the cutout, so the offset is
orthogonal to ⟨cos 2α⟩ and ⟨cos(α+β)⟩ and lands entirely on the constant
and on a_t.  That is exactly the pattern Table 5 shows.

The fix is `reco.harmonic_likelihood_fit_2d`: the same model, the same
`basis_2d` columns, with ε_b carried as a free nuisance per bin and
profiled out.  ε_b is a linear scale, so the profile is closed-form,
ε̂_b = n_b/D_b, and substituting it back gives **exactly** the conditional
multinomial given the bin totals, p_{f,b} = l_f q_{f,b}/D_b.  Profiling B
nuisances that grow with the data would normally raise the Neyman–Scott
worry; it does not here, because the conditional score has zero mean bin
by bin at any count — Σ_f p_f (P_f/q_f − P̄/D) = 0 identically — so there
is no 1/ν_b term anywhere.  Newton from A = 0, six iterations, a
step-halving guard that keeps every populated bin's density positive,
covariance from the conditional Fisher information at the observed
totals, the same rank guard as the ratio (the Fisher information is the
same Gram matrix as the weighted design, so the two fail on the same
acceptances and say so in the same words), 1.5 ms per fit.  The loop
stops on the step, which cannot on its own tell convergence from a
step-halving *stall* on the positivity boundary, so
`_check_profile_convergence` tests the gradient at the solution in units
of the error bar it would move and raises a named `LinAlgError` instead
of returning a boundary point with a meaningless covariance (2026-08-28
review, finding 5; the margin is twelve orders of magnitude at every
occupancy this chain reaches, and the guard first fires around a tenth of
a count per (α, β) cell).  Empty bins
contribute nothing and single-fill bins contribute a finite n log p term,
which removes the R3 variance pathology structurally rather than by a
patch.

Measured at the tagging optics, `--n-mc 6000000`, twenty one-year
pseudo-experiments per |t| bin at 5 × 40.8 on **the same Poisson draws**
(`money_cos2phi_coherent_reco.py --optics tagging --ensemble 20
[--fit likelihood]`):

| \|t\| | N_tag/yr | counts/(α,β) bin | injected a_t | ratio mean | likelihood mean |
|---|---|---|---|---|---|
| 0.05–0.08 | 4.2×10⁵ | 1451 | 0.0898 | 0.0899 (+0.1%) | 0.0902 |
| 0.08–0.12 | 1.8×10⁵ | 620 | 0.1181 | 0.1170 (−0.9%) | 0.1176 |
| 0.12–0.17 | 5.4×10⁴ | 192 | 0.1473 | 0.1394 (−5.4%) | 0.1463 |
| 0.17–0.25 | 1.3×10⁴ | 48 | 0.1803 | 0.1143 (−37%) | 0.173 |

Two hundred experiments turn that comparison from an indication into a
measurement (the same command with `--config 0 --n-mc 6000000
--ensemble 200`, 8 min).  Ratio: 0.0899 / 0.1179 / 0.1414 / 0.1184,
i.e. +0.1 / −0.2 / −4.0 / −34.3 % with pulls of the mean +0.4 / −0.7 /
**−9.0** / **−42.1**.
Likelihood: 0.0902 / 0.1185 / 0.1477 / 0.1816, i.e. +0.4 / +0.3 / +0.3 /
+0.7 % with pulls +2.0 / +1.0 / +0.5 / +0.9 — and the first bin's residual
is the response Monte Carlo's own floor at 6 × 10⁶ recoils, not the
estimator.  At 18 × 137.5 the likelihood gives 0.0998 / 0.1373 / 0.1793 /
0.2298 against 0.0997 / 0.1372 / 0.1791 / 0.2283, pulls +0.8 / +0.6 /
+0.4 / +1.5.  The error bars become honest with it: the likelihood's
spreads are 0.0031 / 0.0048 / 0.0091 / 0.0193 against quoted errors
0.0035 / 0.0048 / 0.0086 / 0.0187, while the ratio's quoted error at 48
counts per bin is 13% *below* its own spread (0.0181 against 0.0208) and
at one count per bin 62% *above* it — a compressed estimator's error bar
is not its spread in either direction.  Push the same bin to ONE count
per (α, β) cell and the ratio changes sign, mean −0.085 against +0.182
injected, while the likelihood still returns +0.181.  On exact counts the
likelihood returns the injected coefficients to machine precision (rel
4 × 10⁻¹⁶ against the ratio's 2 × 10⁻⁶) and the two quote the same Asimov
errors to 0.24%, so **no published error bar moves**.  The default stays
`fit="ratio"`, and a test asserts the default reproduces
`harmonic_ratio_fit_2d` bit for bit.

*Adaptive binning is not the fix.*  In the sparsest bin it only
attenuates — over the same two hundred experiments, −34.3% at 12 × 24,
−13.9% at 8 × 16, −7.8% at 6 × 12 (`--n-alpha 8 --n-beta 16` and
`--n-alpha 6 --n-beta 12` on that same `--ensemble 200` command) —
because the offset falls as 1/ν_b while ν_b grows only as the bin count
falls; it inflates err(a_e) by 17% through the wider-bin dilutions
(0.0139 → 0.0148 → 0.0162); and it hits a hard floor at K_α = 4, where
⟨cos 2α⟩ vanishes identically and the design is rank-deficient at 430
counts per bin.  It survives as a cross-check, not as the remedy.

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
