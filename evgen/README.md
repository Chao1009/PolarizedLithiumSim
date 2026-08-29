# polligen — doubly polarized e+⁶,⁷Li event generation

*Reproducing these results?  [../docs/reproduction_manual.md](../docs/reproduction_manual.md)
§4 has every command with its runtime and expected output, and §5.1 the
PYTHIA 8 sample they consume.*

Event-level Monte Carlo for the *doubly polarized* eA process (plans/05):
polarized electrons on vector/tensor-polarized ⁶,⁷Li, spin-labeled events,
run-plan bookkeeping. Step 5.A (physics kernel + inclusive sampler) and
Step 5.B (tagged mode — spin ⊗ cluster-spectator correlation) are
implemented; the BeAGLE reweighter (5.C) and HepMC3 output (5.D) come
next. Imports `../fastsim/polli_fastsim` — nothing there is duplicated.

## Quick start

```bash
cd evgen
python3 -m pytest tests/ -q            # 305 tests
python3 scripts/closure_fom.py --ion 6Li --events 200000 --trials 200
python3 scripts/closure_fom.py --ion 7Li --events 200000 --trials 200
python3 scripts/money_tagged_azz.py --events 400000       # money plot 4 (--config, --optics, --lumi-fraction)
python3 scripts/tagged_polarimetry_7li.py --events 300000 # 7Li polarimetry + tagged EMC (--config, --optics)
python3 scripts/money_cos2phi.py                          # money plot 5
python3 scripts/money_cos2phi_coherent.py                 # money plot 6
python3 scripts/money_delta_extraction.py                 # money plot 7
python3 scripts/phase_space_bins.py       # (x,Q2) rate maps + binning
python3 scripts/reco_chain_figures.py     # reconstruction-chain figures
python3 scripts/money_cos2phi_reco.py          # money plots 5R + 7R (reco level; --lumi-fraction)
python3 scripts/money_cos2phi_coherent_reco.py --config 0 --optics tagging --n-mc 6000000 # money plot 6R (reco level, tagging optics; --lumi-fraction)
python3 scripts/coherent_optics_scan.py   # WP5: the coherent tag vs the near-beam envelope
python3 scripts/tagging_optics.py         # report 1 §6.1: a lithium tagging optics at IP6, priced in luminosity
python3 scripts/hfs_acceptance.py --config 1 --sample samples/pythia8_e10_p99.5_dis.npz samples/pythia8_e10_n99.5_dis.npz  # report 2 §3: where the hadronic E - p_z sum goes
python3 scripts/target_mass_bound.py       # how big is the γ² term A∥ drops? (prints, no figure)
python3 scripts/eic_beam_figures.py       # report 3: the ion energy menu and the divergence
python3 scripts/nearbeam_aperture_scan.py # plans/09: what every near-beam aperture is worth (--isotope 7Li for the 7Li alpha panel)
python3 scripts/nearbeam_reach_gain.py    # plans/09: the coherent chain at both apertures (--lumi-fraction)
python3 scripts/nearbeam_sensor_budget.py # plans/09: hot-spot Z-ID, sizing, channel count
python3 scripts/nearbeam_zid_power.py     # plans/09: how much charge information Z-ID needs
python3 scripts/nearbeam_two_hit.py       # plans/09 B4: two-hit topology of 6Li -> alpha + d, and the partner-fragment veto
python3 ../reports/build_report.py --pdf  # assemble reports/ pages
python3 ../tools/consistency_check.py     # 24 checks: does everything still agree?
#   reports are numbered in reading order:
#   -> 0 polarized_li_primer.html/pdf (educational physics primer)
#   -> 1 cos2phi_money_plots_report.html/pdf (the projected measurements)
#   -> 2 reconstruction_chain_report.html/pdf (measurement + reco audit)
#   -> 3 eic_epic_reference.html/pdf (the machine and the detector)
#   -> 4 nanowire_far_forward.html/pdf (the near-beam detector study)
#   figures embedded from the PNGs above; display math typeset at
#   build time (matplotlib mathtext, no JS/fonts in the output)
```

Every per-year number below assumes the full 10 fb⁻¹/u programme year in
its own configuration — spin states, far-forward optics and isotope — which
is `--lumi-fraction 1`, the default of the scripts that take it.  A run plan
that gives an observable a share *f* multiplies its statistical errors by
1/√*f* and the years to a target significance by 1/*f*, leaving any
luminosity quoted as a reach where it is; the options are priced in
plans/07 WP2 and the arithmetic is pinned by `tests/test_run_share.py`.

## Modules

| module | content |
|---|---|
| `polligen/nearbeam.py` | plans/09: thin-film energy deposit; the hot-spot firing threshold I_th/I_c = 1 − 2r_s/w with r_s ∝ z (anchored on Argonne's *extrapolated* 134 nm, arXiv:2312.13405); and a sampled Landau with `zid_fake_rate`, which says how much charge information the ⁶Li/α separation actually needs — one bit per plane costs a factor 1.4 against the Neyman-Pearson optimum, and the nanowire loses on fill factor instead |
| `polligen/spin.py` | ρ(m,m′) for J = 1, 3/2: Wigner-d/CG, populations ↔ normalized (vector, tensor, octupole) moments, arbitrary quantization axis, spin-temperature (max-entropy) fills |
| `polligen/xsec.py` | doubly polarized inclusive master formula: HJM spin-1 tensor sector (b₁/b₂, Δ cos 2φ), vector sector A∥ = D·g1/F1 (+ γ-suppressed g_T term, g₂ = g₂^WW), spin-3/2 rank-0/1 exact + rank-2 scenario slots; the exact finite-γ A∥ = D_γ(A₁ + ηA₂) behind `target_mass=True` (default off; what it removes is ≤ 0.6% at the published sweet spots — `scripts/target_mass_bound.py`) |
| `polligen/bookkeeping.py` | run plans (helicity flips, tensor thirds, transverse fills), relative-luminosity offsets + first-order bias formulas, polarimetry smearing, per-(run,bunch) rng streams |
| `polligen/sample.py` | grid inverse-CDF sampler: per-spin-state Poisson rates (φ-averaged modulation shifts counting rates), φ accept-reject, Mode-W weight matrices |
| `polligen/estimators.py` | analysis-side estimators: helicity-flip, tensor thirds, cos 2φ moment + binned LSQ fit (holey-φ robust), luminosity-corrected yields |
| `polligen/tagged.py` | two-cluster spin model: CG-coupled L-waves → m-dependent spectator densities n_M(k,k̂), embedded-cluster spin populations, pair decomposition; TaggedSampler = spin ⊗ spectator ⊗ DIS ⊗ far-forward routing |
| `polligen/coherent.py` | coherent (intact-g.s.) e+⁶Li channel: scenario coherent fraction/slope, intact-recoil kinematics (R = 1.000 → RP pT-tail only), analytic tag acceptance exp(−B pT_cut²), per-bin tagged-rate projections, breakup veto table (plans/06) |
| `polligen/recopseudo.py` | **reconstructed-level pseudo-experiments** (plans/07 WP3, 2026-08-24): importance-sampled inclusive response (`RecoResponse`), exact spin-sorted expected counts per reco bin, ratio-fit measurement, MC bin-centering to Δ; coherent Roman-Pot response (`CoherentResponse`), exact two-azimuth counts with the deformation/gluonic/unpolarized harmonics, acceptance-weighted template basis for the 2-D fit |
| `polligen/reco.py` | **measured quantities and reconstruction** (plans/07 WP3 seed, 2026-08-24): head-on ↔ lab frames with the 25 mrad crossing angle; covariant azimuths φ_S / φ_t from four-vectors (Bacchetta et al. conventions, verified against an explicit collinear-frame construction); electron / hadronic (Σ, parametrized) / mixed kinematic reconstruction with EMCal and tracking resolution models; the spin-state-sorted harmonic ratio estimator (acceptance- and relative-luminosity-immune, 1.5× more precise with m = 0-rich fills); coherent recoil four-vector (exact light-cone t) and Roman-Pot emulation (divergence smearing, rectangular/elliptical 10σ cutout, angular near-beam cut pT_cut = 10σ_θ A p_u) |
| `polligen/radiative.py` | **collinear initial-state radiation as a migration bound** (plans/07 WP4, 2026-08-28): the exponentiated leading-log electron structure function, a per-event sampler with its own random stream, the closed-form observed kinematics of the electron / Σ / JB / double-angle / mixed methods (only the mixed x is ISR-exact; the Q²_e label migrates by 1/(1 − z)), the E − p_z handle, and `migration_bound` / `migration_bound_seeds` — the default-off `RecoResponse(isr=…)` hook and `money_cos2phi_reco.py --isr --isr-seeds …` measure it at +0.5 to +1.2% of Δ̂ in the published generator window and ≤ 2.9% with the low-Q² feed-in opened up, against the letter's 5% gate; a single response seed is not enough for a quoted number, so the bound is a seed average |

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
binned estimators). At P_zz = 0.6 (mid energy, e10 × ⁶Li 99.5 GeV/u) the moment-constrained amplitudes are (0.44–0.95)×10⁻² at the sweet spots x = 0.011–0.14 (table-α_s, the money_delta production convention) with per-bin δA ≈ (1.4–4.5)×10⁻⁴ already in year 1 (10 fb⁻¹/u) — the
measurement discriminates interpretation A vs B, not just zero.

`money_delta_extraction.py` (money plot 7): the same pseudo-
measurements unfolded to the structure function itself —
Δ̂ = Â·y²D_φ/(1−y) per bin with model bin-centering — presented as
xΔ(x,Q²) data points at the three sweet-spot Q² slices with
independent 1-yr/10-yr draws against the moment_A/moment_B curves
(the area under each curve is the Sather–Schmidt moment, ÷3 dilution).

`money_tagged_azz.py` (money plot 4): the tagged tensor asymmetry
A_zz^tag(k) of the α-tagged embedded deuteron, the first tagged spin
observable projected for any A > 2.  Since 2026-08-28 it routes at the
per-configuration optics (`--optics menu` = Yellow Report high-acceptance
plus the lithium tagging optics, `--config {0,1,2}`) with the spectator's
own lab azimuth against the rectangular 10(σ_h, σ_v) envelope, and it
overlays the truth **weighted by each optics' θ_k acceptance** beside the
θ_k = 90° curve.  Both matter.  At 10 × 99.5 the α tag is 2.5% at the
Yellow Report optics against 30% at the tagging optics — a 9% cost in
tagged events per year, since the tagging optics runs at L/L_HA = 1/13.3,
and a gain of 1.9× and 1.2× at the other two configurations — but
the median accepted spectator momentum is 0.32 GeV/c with **nothing below
k = 0.15 GeV/c** against 0.16 GeV/c with 44% below it: at every published
optics the ⁶Li α tag admits only the high-k tail, and the tagging optics
is what turns money plot 4 from a one-point measurement into a curve.  And
the accepted sample is not at θ_k = 90°: the surviving off-rigidity window
slice is longitudinal (⟨|cos θ_k|⟩ = 0.79) while the near-beam tail the
tagging optics opens is transverse (0.40), so the folded A_zz reads +0.49
and −0.07 at k ≈ 0.33 GeV/c where the 90° curve says −0.48.  The ±0.5
swing the pre-2026-08-28 figure showed between its two optics was that
θ_k sculpting, not the wave function (plans/09 B2).

`tagged_polarimetry_7li.py`: the ⁷Li α-tag pair — the in-situ alignment
polarimeter ⟨P₂(cos θ_k)⟩ = −T/5 and the tagged polarized-EMC companion
A_∥^tag(x) on the quasi-free triton.  Both are optics-blind: the ⁷Li α is
off rigidity at R = 0.856 and sits in the middle of the Roman-Pot momentum
window, so it never has to clear the near-beam envelope, and the folded
⟨P₂⟩ slope moves from −0.193 (legacy 73 μrad) to −0.194 (Yellow Report) to
−0.196 (tagging) against the analytic −0.200.  The consequence is a
programme statement: for ⁷Li the tagging optics buys ×1.02 in acceptance
for ×1/8–1/15 in luminosity, a factor 8–15 net loss and the exact inverse
of ⁶Li, so the two isotopes want different machine optics and are
different runs (plans/09 B3).

`phase_space_bins.py`: the companion phase-space figure — (x, Q²)
event-rate maps on the 40×30 analysis grid (1-yr program) for the
inclusive DIS and the RP-tagged coherent channel, with every bin of
the detailed plots overlaid: the four sweet-spot φ′ super-bins
(plot 5), the three Q² slices with their merged x-bin pairs (the
Δ-extraction points of plots 5/7, exact money-plot selection), and
the tagged-count-maximal super-bin of plot 6(d). The tagged sample is concentrated at x < 0.03.

`money_cos2phi_coherent.py`: the coherent channel e⁶Li → e′X ⁶Li(g.s.)
with the intact recoil tagged in the Roman-Pot near-beam pT tail
(R = 1.000 exactly — beam-blind below the pT cut): acceptance exp(−B pT_cut²) = 13.5% [9–20%] for a 0.20 GeV near-beam envelope vs 4×10⁻⁵ at 0.45 GeV; with the Yellow Report divergences and the measured pot aperture no published IP6 optics delivers either (plans/10; Report 1 §6.1). For the 0.20 GeV envelope: 1.7×10⁷ / 1.7×10⁸ tagged events at 10 / 100 fb⁻¹/u (scenario f₀ = 0.04), best-bin δA ≈ 1.8×10⁻³ / 0.6×10⁻³ (5σ floors at 0.9% / 0.3% modulations).
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
  φ_t − φ_S (not the electron azimuth of money plot 6); an anisotropic RP
  cutout fakes ⟨cos 2φ_t⟩ = 0.49 (0.71) at aspect ratio 1.25 (1.5)
  against a physics a₂ ≈ 0.036, and its SIGN follows the aspect — the
  slot-like (2.5 : 1, wide in x) cutout assumed before 2026-08-26 gives a
  large positive fake ⟨cos 2β⟩ about the vertical spin axis, the measured
  aperture and the real, taller-than-wide envelope a large negative one
  (`tools/fullsim/README.md`); the near-beam cut is angular,
  pT_cut = 10σ_θ·6p_u, giving 5.0×10⁻⁷ / 8.4×10⁻²⁶ / 3.7×10⁻¹³ tag
  acceptance at the γ-matched 40.8 / 99.5 / 137.5 GeV/u with the
  per-configuration divergence (plans/10) — but that is the analytic
  CIRCULAR cut at the horizontal 10σ.  The envelope is a rectangle
  10(σ_h, σ_v) and the measured silicon aperture is a second constraint
  per axis: taking the larger of the two per axis gives 9.4×10⁻¹⁰ /
  6.2×10⁻²⁷ / 7.1×10⁻¹⁴, and the aperture scan's envelope-only column
  7.2×10⁻⁸ / 6.2×10⁻²⁷ / 7.1×10⁻¹⁴ against 9.4×10⁻¹⁰ / 2.0×10⁻¹⁹ /
  1.2×10⁻⁵ for the silicon alone, on the aperture re-measured 2026-08-28
  (2.50 / 1.51 / 0.53 mrad against the September-2024 2.00 / 1.35 / 1.03;
  `scripts/nearbeam_aperture_scan.py`, Report 3 Table 6, plans/09 B1).

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
  reco bins (x = 0.028 / 0.011 / 0.071 / 0.141): purity 0.63–0.69, efficiency 0.38–0.66, D = 0.90–0.99; Â within ≈ 2σ of the reco-bin truth; δÂ = 1.2 / 0.9 / 1.6 / 3.0 ×10⁻⁴ (1 yr) — 0.59–0.69 of money plot 5's single-fill errors (the m = 0-rich
  fill gain beats the efficiency loss);
- 7R: best bins δΔ = 2.4×10⁻³ (Q² = 1.14) and 1.2×10⁻³ (3.13 GeV²) in year 1, purities 0.54–0.57;
- 6R (re-derived 2026-08-28 at the lithium tagging optics of Report 1 §6.1,
  `--optics tagging`: horizontal β* × 50 at 5 × 40.8, σ_θ = 33/380 μrad,
  pots following the 10σ envelope 0.33 × 3.8 mrad, L/L_HA = 1/7.1 — with
  the Yellow Report divergences and the measured pot aperture no recoil
  survives at any configuration): acceptance 41%, N_tag = 2.5×10⁶/yr, the
  cutout fakes ⟨cos 2β⟩ = −0.27.  The published binning is the seven bins
  of `recopseudo.T_EDGES_PUBLISHED`, 0.017–0.25 GeV² (adopted 2026-08-28;
  Report 2 §7, plans/08 item 4), of which five carry more than 10⁵ recoils
  (5.8 / 4.2 / 3.0 / 4.2 / 1.8 × 10⁵), and they combine to a one-year
  δa_e = 0.00119.  The template fit recovers a_t(t_ref) = 0.0466 ± 0.0053 /
  0.0563 ± 0.0047 / 0.0672 ± 0.0048 / 0.0915 ± 0.0035 / 0.1099 ± 0.0048 /
  0.1384 ± 0.0085 / 0.1215 ± 0.0174 against 0.045 / 0.059 / 0.071 / 0.090 /
  0.118 / 0.147 / 0.180 injected, and a_e = 0.010 on average: in the lowest
  bin one seed gives 0.0130 ± 0.0022 and the 20-experiment ensemble a mean
  of 0.0092 with a spread of 0.0020, so the single-seed high value is a
  fluctuation and not the "residual template bias" earlier text called it.
  Below 10⁵ recoils the bin-wise ratio of Poisson counts is biased low
  (−4.1% and −35.8% on twenty one-year draws, at 189 and 46 counts per
  (α, β) cell) and closes at ten years.  **`--fit likelihood`** — the
  acceptance-profiled
  conditional-multinomial Poisson likelihood of `reco.harmonic_likelihood_fit_2d`
  (plans/08 A12) — removes that bias outright: on the same twenty draws its
  means are 0.0432 / 0.0585 / 0.0700 / 0.0896 / 0.1179 / 0.1476 / 0.1755
  against 0.0451 / 0.0586 / 0.0711 / 0.0898 / 0.1181 / 0.1473 / 0.1803
  injected.  **`--t-edges`** takes any increasing list and appends `tedges`
  to the output stem, so the pre-2026-08-28 window
  (`recopseudo.T_EDGES_LEGACY`, `0.05,0.08,0.12,0.17,0.25`) stays
  reproducible without overwriting a published figure.
  **`--u-in-situ`** measures (u₁, u₂) from the data
  (`reco.unpolarized_insitu_fit_2d`, plans/08 A3) instead of assuming the
  ZEUS values: δu₂ = 0.0016–0.0136 per |t| bin in one year, 1.8–15× inside
  the ZEUS 1σ of 0.024, which takes the u₂ systematic on a_e from 3–9% to
  negligible in quadrature.  18 × 137.5 gives 6.0×10⁶ tagged per year and
  seven bins, six of them unbiased.  u₁ = 0.05, u₂ = 0.02 sit at the ZEUS LPS 1σ bounds
  (NPB 816:1); see `refs/README.md` for the sources.
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
The production sample is **PYTHIA 8, and it now exists**: 8 M events over the
three beam configurations (`tools/pythia8/gen_dis_hfs.py` builds natively —
no container — and `evgen/samples/README.md` is the manifest).  Without a
sample both scripts fall back to the toy and label it.

```bash
python3 scripts/hfs_resolution.py --config 1     --sample samples/pythia8_e10_p99.5_dis.npz samples/pythia8_e10_n99.5_dis.npz --outdir .
python3 scripts/money_cos2phi_reco.py --y-source hfs --hfs-sample <the same pair> --hfs-calibrate --unfold folded --tag _hfscal --outdir .   # the published 5R/7R
```

**Where Σ_h goes (2026-08-27, `hfs_acceptance.py`, Report 2 §3 Figure 2):** at the mid
sweet spots 78 / 86 / 82 / 91% of the true Σ_h is within the acceptance and
above threshold, 69 / 74 / 73 / 85% is captured through the full response, 19 /
8 / 16 / 7% escapes forward beyond |η| = 3.7 (the target-fragmentation side of
a W ≈ 6–10 GeV system; the lithium fragments at η ≈ 8 never enter).  The
library *reproduces* that capture bias in the pseudo-events;
`HFSResponse(calibrate=True)` / `money_cos2phi_reco.py --hfs-calibrate`
divides by the library's ⟨Σ_reco⟩/⟨Σ_true⟩ in bins of the RECONSTRUCTED
(x_mixed, Q²_e) — the analysis's own hadronic-scale calibration, applied at
the measured point (the 2026-08-28 review found the first version keyed on
the event's true cell) — and the 5R purities go from 0.42 / 0.53 / 0.49 / 0.68
to 0.56 / 0.59 / 0.64 / 0.75 at unchanged errors.

**Measured with PYTHIA (2026-08-28, γ-matched energies), Σ method, 50 MeV calorimeter noise, at the four money-plot-5 sweet spots — which each configuration selects for itself (`hfs_resolution.py` takes them from the same selection as `money_cos2phi.py`):**

| configuration | spot 1 (x, Q²) → δy/y | spot 2 | spot 3 | spot 4 |
|---|---|---|---|---|
| low, 5 × 40.8 | (0.089, 1.14) → 0.39 | (0.035, 1.14) → 0.23 | (0.141, 3.13) → 0.24 | (0.141, 14.3) → 0.11 |
| **mid, 10 × 99.5** | (0.028, 1.14) → **0.32** | (0.011, 1.14) → **0.22** | (0.071, 3.13) → **0.29** | (0.141, 14.3) → **0.15** |
| top, 18 × 137.5 | (0.014, 1.46) → 0.23 | (0.0056, 1.46) → 0.19 | (0.028, 3.13) → 0.21 | (0.141, 14.3) → 0.18 |

Three things the toy could not say.  (i) The 25% stand-in is right in magnitude at the mid configuration but the toy was **optimistic everywhere** (at the pre-correction spots by 0.04–0.05 absolute).  (ii) The electron method alone gives 0.46–1.18 at the mid spots and the sweet spots move to lower y as s grows, so **the low-energy configuration is where the x ≈ 0.14 bins are best measured** (0.24 and 0.11).  (iii) The noise floor is the whole story below y ≈ 0.02: in the lowest y bin (y = 0.005) the Σ resolution runs 0.22 → 0.33 → 0.54 → 1.01 for 0 → 25 → 50 → 100 MeV at the mid configuration, and at y = 0.0125 it runs 0.19 → 0.21 → 0.28 → 0.44, so plans/04 #21 (the ePIC noise and threshold floor at
Σ_h ≈ 0.2–0.5 GeV) sets it, exactly as the note says.

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
