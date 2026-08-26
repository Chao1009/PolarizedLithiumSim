# polligen — doubly polarized e+⁶,⁷Li event generation

Event-level Monte Carlo for the *doubly polarized* eA process (plans/05):
polarized electrons on vector/tensor-polarized ⁶,⁷Li, spin-labeled events,
run-plan bookkeeping. Step 5.A (physics kernel + inclusive sampler) and
Step 5.B (tagged mode — spin ⊗ cluster-spectator correlation) are
implemented; the BeAGLE reweighter (5.C) and HepMC3 output (5.D) come
next. Imports `../fastsim/polli_fastsim` — nothing there is duplicated.

## Quick start

```bash
cd evgen
python3 -m pytest tests/ -q            # 136 tests (grid tests auto-skip)
python3 scripts/closure_fom.py --ion 6Li --events 200000 --trials 200
python3 scripts/closure_fom.py --ion 7Li --events 200000 --trials 200
python3 scripts/money_tagged_azz.py --events 400000       # money plot 4
python3 scripts/tagged_polarimetry_7li.py --events 300000
python3 scripts/money_cos2phi.py                          # money plot 5
python3 scripts/money_cos2phi_coherent.py                 # money plot 6
python3 scripts/money_delta_extraction.py                 # money plot 7
python3 scripts/phase_space_bins.py       # (x,Q2) rate maps + binning
python3 scripts/reco_chain_figures.py     # reconstruction-chain figures
python3 scripts/money_cos2phi_reco.py          # money plots 5R + 7R (reco level)
python3 scripts/money_cos2phi_coherent_reco.py # money plot 6R (reco level)
python3 scripts/coherent_optics_scan.py   # WP5: the coherent tag vs the near-beam envelope
python3 ../reports/build_report.py --pdf  # assemble reports/ pages
#   -> cos2phi_money_plots_report.html/pdf (projection report)
#   -> polarized_li_primer.html/pdf (educational physics primer)
#   -> reconstruction_chain_report.html/pdf (measurement + reco audit)
#   figures embedded from the PNGs above; display math typeset at
#   build time (matplotlib mathtext, no JS/fonts in the output)
```

## Modules

| module | content |
|---|---|
| `polligen/spin.py` | ρ(m,m′) for J = 1, 3/2: Wigner-d/CG, populations ↔ normalized (vector, tensor, octupole) moments, arbitrary quantization axis, spin-temperature (max-entropy) fills |
| `polligen/xsec.py` | doubly polarized inclusive master formula: HJM spin-1 tensor sector (b₁/b₂, Δ cos 2φ), vector sector A∥ = D·g1/F1 (+ γ-suppressed g_T term, g₂ = g₂^WW), spin-3/2 rank-0/1 exact + rank-2 scenario slots |
| `polligen/bookkeeping.py` | run plans (helicity flips, tensor thirds, transverse fills), relative-luminosity offsets + first-order bias formulas, polarimetry smearing, per-(run,bunch) rng streams |
| `polligen/sample.py` | grid inverse-CDF sampler: per-spin-state Poisson rates (φ-averaged modulation shifts counting rates), φ accept-reject, Mode-W weight matrices |
| `polligen/estimators.py` | analysis-side estimators: helicity-flip, tensor thirds, cos 2φ moment + binned LSQ fit (holey-φ robust), luminosity-corrected yields |
| `polligen/tagged.py` | two-cluster spin model: CG-coupled L-waves → m-dependent spectator densities n_M(k,k̂), embedded-cluster spin populations, pair decomposition; TaggedSampler = spin ⊗ spectator ⊗ DIS ⊗ far-forward routing |
| `polligen/coherent.py` | coherent (intact-g.s.) e+⁶Li channel: scenario coherent fraction/slope, intact-recoil kinematics (R = 1.000 → RP pT-tail only), analytic tag acceptance exp(−B pT_cut²), per-bin tagged-rate projections, breakup veto table (plans/06) |
| `polligen/recopseudo.py` | **reconstructed-level pseudo-experiments** (plans/07 WP3, 2026-08-24): importance-sampled inclusive response (`RecoResponse`), exact spin-sorted expected counts per reco bin, ratio-fit measurement, MC bin-centering to Δ; coherent Roman-Pot response (`CoherentResponse`), exact two-azimuth counts with the deformation/gluonic/unpolarized harmonics, acceptance-weighted template basis for the 2-D fit |
| `polligen/reco.py` | **measured quantities and reconstruction** (plans/07 WP3 seed, 2026-08-24): head-on ↔ lab frames with the 25 mrad crossing angle; covariant azimuths φ_S / φ_t from four-vectors (Bacchetta et al. conventions, verified against an explicit collinear-frame construction); electron / hadronic (Σ, parametrized) / mixed kinematic reconstruction with EMCal and tracking resolution models; the spin-state-sorted harmonic ratio estimator (acceptance- and relative-luminosity-immune, 1.5× more precise with m = 0-rich fills); coherent recoil four-vector (exact light-cone t) and Roman-Pot emulation (divergence smearing, rectangular/elliptical 10σ cutout, angular near-beam cut pT_cut = 10σ_θ A p_u) |

## Money plots 5–6 (2026-08-10): cos 2φ as projected data points

`money_cos2phi.py`: inclusive gluonometry pseudo-data (φ′ modulation
with statistical error bars) in the four sweet-spot (x, Q²) super-bins,
picked per Q² band from the significance map, plus amplitude-vs-x with
1-year AND 10-year error bars. The Δ(x,Q²) input comes from the unified
registry `polli_fastsim.delta_models` (default `moment_A`: the
sum-rule-constrained ansatz of the merged money_delta suite,
∫xΔdx = −0.012·α_s, with the ⁶Li per-nucleon dilution 1/3;
`moment_B`/`toy` one flag away). Full-luminosity statistics via
per-φ-bin Poisson draws (`sample.phi_histogram_pseudo` — exact for
binned estimators). At P_zz = 0.6 (mid energy) the moment-constrained
amplitudes are (0.7–1.2)×10⁻² at the sweet spots (table-α_s, the
money_delta production convention) with per-bin δA ≈ (1.5–4.5)×10⁻⁴
already in year 1 (10 fb⁻¹/u) — the
measurement discriminates interpretation A vs B, not just zero.

`money_delta_extraction.py` (money plot 7): the same pseudo-
measurements unfolded to the structure function itself —
Δ̂ = Â·y²D_φ/(1−y) per bin with model bin-centering — presented as
xΔ(x,Q²) data points at the three sweet-spot Q² slices with
independent 1-yr/10-yr draws against the moment_A/moment_B curves
(the area under each curve is the Sather–Schmidt moment, ÷3 dilution).

`phase_space_bins.py`: the companion phase-space figure — (x, Q²)
event-rate maps on the 40×30 analysis grid (1-yr program) for the
inclusive DIS and the RP-tagged coherent channel, with every bin of
the detailed plots overlaid: the four sweet-spot φ′ super-bins
(plot 5), the three Q² slices with their merged x-bin pairs (the
Δ-extraction points of plots 5/7, exact money-plot selection), and
the tagged-count-maximal super-bin of plot 6(d). 98% of the tagged
sample sits at x < 0.03 (median 3.5×10⁻³).

`money_cos2phi_coherent.py`: the coherent channel e⁶Li → e′X ⁶Li(g.s.)
with the intact recoil tagged in the Roman-Pot near-beam pT tail
(R = 1.000 exactly — beam-blind below the pT cut): acceptance
exp(−B pT_cut²) = 13.5% [9–20%] with high-acceptance optics vs 4×10⁻⁵
with high-divergence; ~10⁸ tagged events at 100 fb⁻¹/u (scenario
f₀ = 0.04), best-bin δA ≈ 6×10⁻⁴ (5σ floor at a 0.3% modulation).
The modulation is anchored (plans/06 §6.4b): a t-linear deformation
term scaled from the polarized-deuteron CGC calculation (PLB
858:139053, digitized in `coherent.MANTYSAARI_A2_DEUTERON`) giving
⟨a₂⟩_tag ≈ 0.036 [0.018–0.059] at P_zz = 0.6 with a predicted sign
flip vs the deuteron, plus a flat gluon-transversity scenario
(3×10⁻³–10⁻²) — separable by t-shape and x_P dependence.
Backgrounds and the RP charge-ID question: plans/06.

## Reconstruction chain (2026-08-24): what is measured, and does the simulation do it?

`reports/reconstruction_chain_report` describes, for the inclusive and
the coherent cos 2φ measurements, what the detector records, how the
azimuth of the observable is built from the reconstructed four-vectors,
the kinematic reconstruction, the spin-state-sorted estimator, the
Roman-Pot measurement of the recoil, and the conversion to Δ — and
audits each step against the code. `scripts/reco_chain_figures.py`
(→ `reco_chain_inclusive_6Li.png`, `reco_chain_coherent_6Li.png`)
quantifies the findings with `polligen/reco.py`:

- three of the four sweet spots of money plot 5 sit at y = 0.010–0.025,
  where e′ alone gives δy/y = 50–120%: x needs the hadronic (Σ/JB) y;
  with a 15–20% hadronic resolution the super-bins keep 75–83% purity
  and the reco-bin amplitude is within 1–4% of the true-bin value
  (0.64–0.68 and 1–9% at the 25% default of the R-plots);
- the single-fill cos 2φ′ fit measures the detector's own cos 2φ′
  acceptance harmonic ÷ P_zz (a 3% harmonic aligned with the vertical
  spin axis fakes 4× the signal); the spin-state ratio of m = ±1-rich
  (P_zz = +0.6) and m = 0-rich (−1.2) bunches cancels any acceptance
  common to both states bin by bin (bunch-by-bunch alternation required:
  a 10⁻³ difference of the harmonic between the states fakes half the
  signal — code review 2026-08-25), turns relative luminosity into a
  constant, and is 1.5× more precise (`reco.harmonic_ratio_fit`);
- the crossing angle changes only odd harmonics of the e′ azimuth
  (≈ 10⁻³ cos φ, cos 2φ untouched) because the ePIC axis is the electron
  beam; the covariant φ_S equals φ_e − φ_S exactly (massless target),
  O(γ²) < 5 mrad for the nucleus;
- coherent: the anchored deformation term modulates the recoil azimuth
  φ_t − φ_S (not the electron azimuth of money plot 6); the RP cutout
  fakes ⟨cos 2φ_t⟩ = 0.49 (0.71) for aspect ratio 1.25 (1.5) against a
  physics a₂ ≈ 0.036; the near-beam cut is angular, pT_cut = 10σ_θ·6p_u,
  giving 67% / 9% / 10⁻⁸ tag acceptance at 20.5 / 50 / 137.5 GeV/u.

### Reconstructed-level closure (money plots 5R / 7R / 6R, same day)

`polligen/recopseudo.py` composes the pieces into the WP3 pseudo-
experiments: `RecoResponse` (importance-sampled response: 400
pseudo-events per sampler cell through lab frame → smearing → head-on →
electron/mixed reconstruction → reco cuts → ε_eID → covariant φ′),
`measure_inclusive` (exact spin-sorted φ′ counts per reco bin, Poisson,
`reco.harmonic_ratio_fit`), `delta_from_amplitude` (MC bin-centering with
migration), `CoherentResponse` (RP emulation + exact two-azimuth counts)
and `measure_coherent` (template fit `reco.harmonic_ratio_fit_2d` with
the acceptance-weighted basis). Results (mid energy, P_zz = 0.6, moment_A,
mixed method 25%, ε(φ′) harmonic + 10⁻³ rel-lumi offset on):

- 5R (25% hadronic y — the ePIC kinematic-fit study's smearing and the
  ATHENA Fig. 22 value at y ≈ 0.01, `refs/README.md`): sweet spots 1–4 in
  reco bins: purity 0.64–0.68, efficiency 0.37–0.64, D = 0.91–0.99;
  Â unbiased vs the reco-bin truth; δÂ = 1.2 / 1.0 / 1.8 / 3.0 ×10⁻⁴
  (1 yr) — 0.65–0.70 of money plot 5's single-fill errors (the m = 0-rich
  fill gain beats the efficiency loss);
- 7R: best bins δΔ = 1.0×10⁻³ (Q² = 1.14) and 0.5×10⁻³ (3.13 GeV²) in
  year 1, purities ≈ 0.55;
- 6R: angular envelope (0.22 GeV at 50 GeV/u) + slot-like cutout
  (|p_x| < 0.55, |p_y| < 0.22 GeV — the ePIC sensor planes surround a
  horizontal slot, Jentsch DIS 2023) → N_tag = 2.7×10⁶ (1 yr); the cutout
  fakes ⟨cos 2β⟩ = +0.77; the template fit recovers a_t(t_ref) =
  0.117 ± 0.007 (inj. 0.119) … 0.331 ± 0.018 (0.343) and a_e = 0.0101 ±
  0.0013 (inj. 0.010) in year 1 — the slot concentrates the recoils where
  cos 2β hardly varies, costing ~3× in the deformation-term error relative
  to a square cutout (a_e untouched). u₁ = 0.05, u₂ = 0.02 sit at the ZEUS
  LPS 1σ bounds (NPB 816:1); see `refs/README.md` for the sources.
  **Convention (verified in the paper's Eq. 9):** arXiv:2408.13213 expands
  1 + 2Σ a_n cos nΦ with Φ the vector-meson (recoil) azimuth relative to
  the polarization axis, so the deformation modulation coefficient is
  2a₂ — money plot 6 injects a₂ and is conservative by ×2
  (`coherent.cos2phi_coefficient_deformation`).

### Hadronic final state and hadron-side detection (WP3-HFS, 2026-08-25)

`polligen/hfs.py` replaces the 25% Gaussian stand-in for the hadronic y by a
hadronic final state through a hadron-side detector response: `HFSSample`
(generator-independent .npz), exact hadronic sums and the Σ / Jacquet–Blondel /
double-angle / mixed methods, `HadronResponse` (tracker |η| ≤ 3.5, p_T > 0.2 GeV,
95% efficiency; calorimeters |η| ≤ 3.7 with photon/neutral-hadron thresholds,
Yellow-Report resolutions, 50 MeV noise), `ToyHFS` (string-fragmentation
stand-in, exact four-momentum closure, π⁰ → γγ) and `HFSLibrary`/`HFSResponse`
((x, Q²) library transferring a generator's response onto the pseudo-events).
The production sample is PYTHIA 8 (`tools/pythia8/gen_dis_hfs.py`, runs in
eic-shell); locally the toy is used and labelled.

```bash
python3 scripts/hfs_resolution.py --outdir .                 # Figure 3: resolution vs y, captured Σ, sweet spots, noise scan
python3 scripts/money_cos2phi_reco.py --y-source hfs --outdir .   # 5R/7R with the HFS-based y (…_hfs.png)
python3 scripts/hfs_resolution.py --sample samples/pythia8_e10_p50_dis.npz samples/pythia8_e10_n50_dis.npz
```

Toy result (illustrative): captured Σ fraction 0.90; Σ-method δy/y at the four
sweet spots 0.28 / 0.17 / 0.24 / 0.07 with 50 MeV noise (9–12% without noise) —
the 25% stand-in is the noise floor acting on Σ_h ≈ 0.2 GeV; 5R rerun with the
HFS-based y reproduces the Table 2 errors (purities 0.60 / 0.68 / 0.68 / 0.86).

## Step-5.A validation gates (plans/05 §5.4) — all passing

1. **ρ moments, all axes** — exact against analytic Wigner rotations
   (`test_spin.py`).
2. **Master formula ↔ `asymmetries.py`** — bin-wise identity (rtol 1e-12)
   for A∥, A_zz, A_cos2φ, toy *and* CT18/NNPDFpol grid backends
   (`test_xsec_identity.py`).
3. **Pseudo-experiment estimator closure** — pulls unbiased, spreads match
   `err_a_parallel`/`err_azz`/`err_cos2phi_amplitude` within 15%
   (`test_pseudoexp.py`); full FOM-map closure over ~65 x-bins in
   `closure_fom_{6,7}Li.png` (spread/analytic compatible with 1 across
   2.5 decades in x for all three observables).
4. **φ-modulation recovery** — cos 2φ amplitude unbiased with uniform and
   holey acceptance (binned fit; the naive moment demonstrably breaks).

**First systematics numbers** (relative-luminosity offset, validated at
δ = 2% by pseudo-experiments, quoted at the plans/05 reference δ = 10⁻⁴,
in-ring P_e = 0.7, P_z = 0.7, P_zz = 0.6):
bias(A_zz) = −(2/3)δ/P_zz ≈ −1.1×10⁻⁴ and bias(A∥) = δ/(2P_eP_z) ≈
+1.0×10⁻⁴ for the naive estimators — below the ~10⁻³ per-x-bin statistical
floor at 10 fb⁻¹/u but not negligible after combining bins; the
luminosity-corrected estimators remove both exactly.

## Conventions and caveats

- Populations are ordered m = +J … −J everywhere.
- All SFs per nucleon (F2A/A, as in `fom.py`); luminosity per nucleon.
- Vector-sector y-factors use the fastsim approximation set:
  A∥ = D·A1 (γ²-suppressed η term dropped); A⊥ = d·(A2 − ξA1)
  = d·γ·((y/2)g₁ + g₂)/F₁ at leading γ (both O(γ) pieces kept — the
  ξA1 term is the same order as A2). The exact Cosyn–Weiss factors
  can replace two functions behind the same interface when adopted.
- Spin-1 vector fills at |P_z| > 2/3 require tensor polarization
  (positivity); `helicity_flip_plan` defaults to spin-temperature
  populations, which record their implied P_zz in `pzz_true`.
- Modulation amplitudes are evaluated at cell centers of the sampling
  grid (default 100×72, ~4× finer than the 40×30 FOM analysis binning);
  adequate for binned estimators, not for unbinned in-cell shapes.
- Spin-3/2 rank-2 structure functions are scenario inputs (zero by
  default) pending the theory note (plans/04 #14); rank-3 dropped.
