# Review of the simulation and reconstruction toolchain (2026-08-28)

Scope: `evgen/polligen/{reco,hfs,recopseudo,xsec,sample,spin,bookkeeping,
coherent}.py`, the reconstructed-level scripts (`money_cos2phi_reco.py`,
`money_cos2phi_coherent_reco.py`, `reco_chain_figures.py`,
`hfs_resolution.py`, `hfs_acceptance.py`, `nearbeam_reach_gain.py`,
`coherent_optics_scan.py`, `tagging_optics.py`), the PYTHIA production
(`tools/pythia8/gen_dis_hfs.py`, `evgen/samples`), the fast-simulation
inputs the chain consumes (`beams`, `farforward`, `structure`,
`delta_models`, `fom`), the three test files of the chain, and the
documents that quote its numbers (Report 2 above all, Reports 1 and 4, the
READMEs, the reproduction manual, plans/04, 08, 10).

Three questions were put to the code: **(1)** is each step right —
kinematics, frames, angles, estimator algebra, statistics, the response
transfer, importance-sampling weights, conventions between modules;
**(2)** does any *analysis-side* step consume information an experiment
cannot have — generator truth, the injected model, exact nuisances —
outside the three legitimate places (the data-generation model, corrections
built from the simulation and applied as functions of reconstructed
variables, closure references that are plotted but never fitted);
**(3)** is every measured quantity one ePIC can deliver, and is every
idealisation labelled, in the code and in the report.

Method: nine independent readings, one lens each (inclusive truth leaks,
coherent truth leaks, the hadronic-final-state layer, kinematics and
frames, estimator statistics, kernel and sampler, code-versus-document
consistency, experimental realism, tests and scripts), returning 65
findings with the code lines and the numerical check behind each.  The
first was verified by three adversarial reviewers and reproduced
independently; the session limit interrupted the rest of the automated
verification, and every remaining code-level finding was then verified
numerically by hand before being acted on (the scripts are in the session
scratchpad; the reproduced numbers are quoted below).  Documentation
findings were resolved by rewriting Report 2 and restating the other
documents.  Baseline: 209 evgen + 57 fastsim tests green, 22/22
consistency checks; after the pass 225 + 57 and 22/22.

## 1. Verdict

* **One analysis step read generator truth.**  `HFSResponse(calibrate=True)`
  divided each pseudo-event's measured hadronic sums by the library's
  captured fraction of the event's *true* (x, Q²) cell.  Fixed: the
  calibration is a map of ⟨Σ_reco⟩/⟨Σ_true⟩ in bins of the library's own
  *reconstructed* (x_mixed, Q²_e), built with the library electron smeared
  as the chain smears it and the noise in, and looked up at the
  pseudo-event's uncalibrated reconstructed point.  Nothing else on the
  paths from counts to Â, a_t, a_e reads truth (§3).
* **One input was wrong by construction.**  The PYTHIA library carried
  PYTHIA's default `PhaseSpace:mHatMin = 4 GeV`, which for DIS (m̂² = x s)
  removes x < 16/s: 0.004 at 10 × 99.5, 39% of the selected rate of the
  pseudo-experiments (0.2–0.7% of the four sweet spots' rate, but the
  whole low-x half of the Q² = 1.14 GeV² slice).  Regenerated with
  `mHatMin = 0.5`.
* **Two statistical defects of the estimator**, both in low-count bins:
  the bin variance vanished for a bin populated by one spin state
  (R = w_f exactly), giving it an infinite weight — the cause of the
  "rank-deficient" |t| bins of the earlier Table 3 — and the ratio
  inversion's error Jacobian dropped its second-order denominator.  A
  third, deeper than either, was found and closed the same day: the ratio
  inversion is a nonlinear function of Poisson counts, so it carries an
  O(1/ν_b) Jensen offset that biases a_t **low** below ~10² counts per
  (α, β) bin, no matter how the variance is written; the
  acceptance-profiled conditional-multinomial likelihood
  `reco.harmonic_likelihood_fit_2d` removes it, returning 0.0902 / 0.1176
  / 0.146 / 0.173 on the same draws where the ratio gives −4% and −34% in
  the two sparsest |t| bins (plans/08 A12, Report 2 §4.3 and Table 5;
  `--fit likelihood`, with the ratio kept as the default and as the
  published record).
* **The coherent closure was on the wrong geometry.**  6R ran on the
  legacy 73 μrad isotropic divergence with a 0.73 mrad "near-beam
  approach" while §8.1 of the same report listed the Yellow Report values;
  the "residual template bias" it reported on a_e was a 2σ fluctuation.
  Re-derived at the tagging optics of Report 1 §6.1 at all three
  configurations, with ensembles and exact-count systematics.
* **Consistency**: eight smaller code inconsistencies fixed (§2);
  Report 2 rewritten; Report 1, Report 4, the READMEs, the manual, the
  sample manifest and the plans restated.  No headline number moved by
  more than its statistical error.

## 2. Findings and their disposition

Severity as verified.  "Impact" is on the numbers published before this
pass.

| # | finding | severity | verified by | disposition |
|---|---|---|---|---|
| R1 | `hfs.HFSResponse.hadronic` calibrated at `library.cell_means(x_true, q2_true)` | HIGH, truth leak | code + three adversarial reviewers + prototype: reco-keyed map gives purities 0.56/0.59/0.64/0.75, D 1.00/1.06/0.95/0.98 vs truth-keyed 0.52/0.59/0.59/0.76, 1.01/1.03/0.95/0.98 | `HFSLibrary.cell_means_reco`, `HFSResponse.hadronic` two-pass; test forbids the truth lookup |
| R2 | PYTHIA `mHatMin` default removed x < 16/s | HIGH, input | sample x min vs 16/s (1.8% of events below, spectrum ratio 0.00–0.01 below x = 0.003); 39.2% of selected rate; 40k-event test with `mHatMin = 0.5` follows the rate map to ±25% | `--mhat-min 0.5`; library regenerated (σ +37–48%) |
| R3 | `reco.spin_state_ratio` var = 0 for single-state bins | HIGH, statistics | toy counts `[[100,50,0,3],[90,60,5,0]]` → var `[.., 2e-33, 0]`, fit raised `LinAlgError` | expected-count variance with the first-order T clipped; tests (finite variance, Poisson spread reproduced to 6%) |
| R4 | `_ratio_to_modulation` Jacobian (1+u+P̄T)/σ² instead of (1+u+P̄T)/(σ² − P̄R) | LOW, statistics | analytic; 0.3% on errors at P̄T ~ 10⁻³ | fixed; test |
| R4b | the ratio inversion's O(1/ν_b) Jensen offset — a bias of a_t low that R3 and R4 do not touch, since it is in the estimator's form and not in its variance | MEDIUM, statistics | 200 pseudo-experiments per \|t\| bin at the tagging optics, 5 × 40.8: ratio −4.0% and −34.3% (pulls −9.0, −42.1) at 192 and 48 counts per (α, β) cell, and a sign change at 1 count; adaptive binning only attenuates it (−34.3% → −13.9% → −7.8% from 12 × 24 through 8 × 16 to 6 × 12, at +17% on δa_e) | `reco.harmonic_likelihood_fit_2d`, the acceptance-profiled conditional-multinomial likelihood, `--fit likelihood` (plans/08 A12); pulls +0.5…+2.0 at every bin, Asimov closure 4×10⁻¹⁶, errors within 0.24% of the ratio's; 16 tests in `evgen/tests/test_likelihood_fit.py` (14 for the likelihood, 2 for the in-situ u fit) |
| R5 | hadron acceptance in the head-on frame, electron in the lab | MEDIUM, inconsistency | captured Σ within \|η\| ≤ 3.7 at the four spots: head-on 0.816/0.918/0.837/0.927, lab 0.831/0.925/0.850/0.934 | `HadronResponse(xing=)`: acceptance and smearing in the lab, sums in the head-on frame; test |
| R6 | library f = Σ_reco/Σ_true(massive) applied to a massless 2E_e y | MEDIUM, physics | Σ_true − 2E_e y = 4.5 MeV in the library (2.3% at y ≈ 0.01) | the pseudo-event Σ carries the same mass term; test |
| R7 | e + p / e + n merged by event count | LOW, inconsistency | σ_p/σ_n = 1.11–1.16 | `concatenate` weights by σ_gen/n; `HFSLibrary` draws by weight; test |
| R8 | beam-energy spread with the same sign on Q²_e and 1 − y_e | LOW, physics | the mixed x should be exactly beam-energy independent | fixed; test pins Q²_e(1 − y_e) = truth with the spread on |
| R9 | `reco.beam_fourvectors` used A·M_U, `beams`/`spectator` the physical mass | LOW, inconsistency | 12.6 MeV | physical mass; test |
| R10 | `reco.SIGMA_THETA_HD` 149 μrad vs `farforward.HIGH_DIVERGENCE` 164 μrad | LOW, inconsistency | plans/10 A1b | unified; test |
| R11 | 5R/7R bins selected partly on the truth amplitude (`abs(a_reco_bin) < 1e-5`) | LOW, truth leak (inert) | never fires for a registry model | removed; test |
| R12 | `--envelope-split` a no-op wherever the aperture binds | MEDIUM, inconsistency | `_apply_cut` took max(envelope, aperture) before scaling | perturbation on the binding half-widths, on the axis the recoils escape through |
| R13 | 6R on the legacy divergence; "residual template bias" on a_e | HIGH, inconsistency / doc | ensemble of 20: a_e means 0.0092 ± 0.0007, 0.0097 ± 0.0008 vs 0.010; a_t unbiased above 10⁵ recoils per bin | `--optics tagging`, `--ensemble`, `--exact`, `--n-alpha/--n-beta`; Table 5 of Report 2 |
| R14 | ten-year a_t errors omitted the response's own MC statistics | LOW, statistics | 6×10⁵ recoils vs 10-yr errors | published runs at 6×10⁶ |
| R15 | `delta_models` fell back silently to LO α_s without `parton` (+14% on Δ) | MEDIUM, inconsistency | — | warning printed once; manual states it |
| R16 | ISR statement "x → x/(1 − z)" in Report 2 §3 and plans/08 D3 | MEDIUM, physics (doc) | (1 − z) cancels between Q²_e and s for y_Σ | corrected: the Q² label migrates, x does not; *closed 2026-08-28 (run 13)*: pinned by `test_mixed_x_is_exact_and_only_the_q2_label_migrates`, and the migration bound itself measured by `polligen/radiative.py` (plans/07 WP4: +0.50…+1.22% of Δ̂ at the sweet spots in the published window — the eight-seed mean ± sem +0.62 ± 0.03 / +0.50 ± 0.02 / +0.94 ± 0.03 / +1.22 ± 0.02% — and ≤ 2.9% with the low-Q² feed-in; gate ≤ 5% passed) |
| R17 | Report 2 Table 1 at the pre-correction sweet spots; callout (3) with 50 GeV/u numbers; §4 resolution numbers at 50 GeV/u; the 73/149 μrad footer; §8.1 "PYTHIA" as the value used against Table 2's stand-in; §8.2 without the acceptance-stability requirement; 3% vs 5% polarimetry | HIGH–LOW, doc | read against the scripts | Report 2 rewritten |
| R18 | Report 4 near-beam table on the legacy divergence | MEDIUM, doc | plans/10 A4 | scope banner; A4 stays open |
| R19 | stale sample names (e10_p50), test counts (183/209 vs 225), 6R command in the manual | LOW, doc | — | restated |
| R20 | docstrings: beams 40.7/99.3; delta_models ⟨Q²⟩ = 3.9; hfs "YR requirement values" for stand-ins; the 3 mrad defence citing an unverifiable slide | LOW, doc | — | rewritten |
| R21 | no test of a spin-axis tilt in the inclusive chain; no external test of the hadronic methods; no test of the tagging helper against the pricing script | LOW, test gap | — | tests added (tan 2δ recovered to 10⁻³; DA/JB/Σ against hand formulas to 10⁻⁹; optimum r_h = 49.7/175.6/89.3 reproduced) |

Findings examined and **not** acted on, with the reason: the cell-centre
amplitude and the σ(centre) × area cell rates of the importance sampling
(≤ 0.2% on the amplitude, +0.2–0.9% on N, documented in Report 2 §4.4);
the per-event sin 2δφ′ term dropped from the expected counts (averages to
zero); the Gaussian noise stand-in (labelled as such); the two EMCal tables
in `hfs.py` and `reco.py` (different purposes, both documented); P_zz taken
as exactly known (its scale is a 1:1 normalisation, pinned by a test, quoted
as a systematic); u₁, u₂ taken as known (fit arguments, the in-situ
measurement stated as the requirement); the coherent pseudo-experiment
without an electron side (stated); the dispersion degeneracy of the
off-momentum recoil (below every divergence, documented in
`rp_measure`); the template's dependence on the divergence model (a
systematic to be varied with the optics, stated).

## 3. Measurability, re-audited after the fixes

Paths traced: `RecoResponse.__init__` → `measure_inclusive` →
`harmonic_ratio_fit` → `delta_from_amplitude` / `fold_shape_fit`;
`CoherentResponse` → `measure_coherent` → `harmonic_ratio_fit_2d`.

| quantity | used by | source | measurable |
|---|---|---|---|
| spin state, P_zz, L_f per state | the ratio weights | run-plan bookkeeping | yes (bunch pattern, polarimeter, luminosity monitor); the offset the plan carries is unknown to the analysis |
| k′ smeared, in the lab, transformed to head-on | Q²_e, θ_e, φ′ | detector response | yes |
| Σ_h, p_T,h | y_Σ, x_mixed | library transfer + noise (data generation) | yes |
| calibration factor | y_Σ | `cell_means_reco(x_mixed^uncal, Q²_e)` | yes — reconstructed bins of the analysis's MC |
| analysis cuts, bin assignment | selection | reconstructed quantities | yes |
| φ′ | the histogram | covariant formula on k, k′_measured, P, S | yes |
| purity, efficiency, D, K, fold | corrections and closure references | the response | yes (MC-derived), never fitted |
| Roman-Pot angle pair | β, \|t\| | `rp_measure` (data generation) | yes |
| β basis, t_ref | the 2-D fit | the response per reconstructed (β, \|t\|) bin | yes (MC-derived) |
| u₁, u₂, luminosity shares | the 2-D fit | assumed, as fit arguments | assumption, exercised |
| truth references a_t(t_ref), a_e | plotted | `truth_reference` | closure only |

The bins shown in the figures are selected on their expected count and
statistical error.  Generator truth is read to generate the pseudo-data,
to build the corrections above, and for the closure references.
