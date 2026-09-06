# PolarizedLithiumSim — Polarized ⁶,⁷Li at the EIC

Simulation program building the physics case for polarized ⁶Li and ⁷Li
ion beams at the Electron-Ion Collider, in support of the ANL polarized
Li ion-source development (`docs/ecrp_2026_proposal.pdf`,
NOFO DE-FOA-0003602). Target venue for Phase-1 results: the INT program
on polarized ion beams at the EIC, March 22 – April 2, 2027.

## The physics in three lines

1. **⁶Li (spin-1)** — practical replacement for polarized deuterons:
   tensor structure function b₁, the purely gluonic double-helicity-flip
   Δ(x,Q²) ("nuclear gluonometry", cos 2φ modulation), and — new — the
   coherent intact-recoil channel with its deformation/exotic-glue
   decomposition.
2. **⁷Li (spin-3/2)** — an effective polarized proton in-medium
   (P_p ≈ 0.87): the polarized EMC effect far beyond JLab reach.
3. **Collider-mode spectator tagging** — α/t/n fragments in the
   far-forward detectors select the struck cluster; tagging purity is
   the simulation deliverable the proposal explicitly calls for.

## Repository map

| directory | content |
|---|---|
| `plans/` | the program plan and findings — start at [plans/00_README.md](plans/00_README.md) (document map + development-run log) |
| `fastsim/` | `polli_fastsim`: analytic fast simulation — rates, FOM maps, tagging acceptance, money plots 1–3, and the dated `money_delta_*` study suite (moment-constrained Δ ansatz, detector efficiency, reco selection) with working notes in `fastsim/notes/` |
| `evgen/` | `polligen`: the doubly polarized e+⁶,⁷Li event generator (the first of its kind) — spin-density kernel, inclusive + tagged samplers, run-plan bookkeeping, estimators, the coherent intact-⁶Li channel, money plots 4–7 and their reconstructed-level versions 5R/7R/6R, the hadronic-final-state layer with the hadron-side detector response (`polligen/hfs.py`) on the standing **PYTHIA 8 production** (`tools/pythia8`, manifest in `evgen/samples/README.md`), the systematics the spin-state ratio cannot cancel (per-fill acceptance, assumed luminosity and nuisance harmonics, energy scales), the forward-folded amplitude response that replaces the bin-by-bin K, and the WP5 coherent optics scan |
| `docs/` | the **[reproduction manual](docs/reproduction_manual.md)** — environment, every command, expected numbers, and the third-party generators — plus source documents, the self-contained physics notes [docs/note_cos2phi_coherent_6Li.md](docs/note_cos2phi_coherent_6Li.md) (verified 50-entry bibliography; the state of 2026-08-10 — its projections are superseded by Report 1 §6–7) and [docs/note_7li_theory_questions.md](docs/note_7li_theory_questions.md) (the six theory inputs a ⁷Li program waits on and the two ⁶Li asks that travel with them, each with its ask and its interim assumption), and the audit the current state was built against, [docs/code_review_2026-08-28.md](docs/code_review_2026-08-28.md) (nine readings of the reconstruction chain, 65 findings, the code-level ones re-verified by hand and the documentation ones by rewriting); its predecessor [docs/code_review_2026-08-25.md](docs/code_review_2026-08-25.md) (measurability audit of 5R/7R/6R, findings F1-F13) is kept as the dated record, and the consistency review of Reports 0–4 [docs/consistency_review_2026-09-02.md](docs/consistency_review_2026-09-02.md) (2026-09-02; 260 findings — the fixes and the 53 wording suggestions applied in run 17, the ten mechanical checks it specifies now modules under `tools/checks/`, and its three author decisions registered as items 22–24 of [plans/04_open_questions.md](plans/04_open_questions.md)) |
| `reports/` | the five circulate-able reports (self-contained HTML + rendered PDF), numbered 0–4 in reading order: **0** the educational primer, **1** the cos 2φ money-plot report, **2** the reconstruction-chain analysis note (what is measured, how the azimuth and Δ are reconstructed, audit of the simulation), **3** the EIC/ePIC machine-and-detector reference, **4** the near-beam far-forward study — all restated 2026-08-29 on the per-configuration transport and corrected 2026-09-02 on the cross-report consistency review, with Reports 0, 2, 3 and 4 re-sourced 2026-09-06; served as a website by the GitHub Pages workflow (`.github/workflows/pages.yml` + `reports/index.html`) — activate once via Settings → Pages → Source: "GitHub Actions" (public repo required on free plans) |
| `tools/` | PYTHIA 8 hadronic-final-state generation (`tools/pythia8`), the BeAGLE route and the e+d control that calibrates the cluster model's tail (`tools/beagle`, `tools/analysis`), and the full ePIC chain — including the ion gun that shoots an intact ⁶Li at the far-forward detectors (`tools/fullsim`) |

## Quick start

Everything below, and everything else this program publishes, is
reproducible from a clean Linux box: **[docs/reproduction_manual.md](docs/reproduction_manual.md)**
is the map from result to command — what to install, what to run, how
long it takes, what the answer should be, and what cannot be reproduced
here and why.

```bash
# fast simulation (121 tests; PDF-grid tests need the `parton` package
# with CT18NLO, CT18ANLO, EPPS21nlo_CT18Anlo_Li6, nNNPDF30_nlo_as_0118_A6_Z3
# and NNPDFpol11_100 installed)
cd fastsim && python -m pytest tests/ -q

# event generator (323 tests)
cd evgen && python -m pytest tests/ -q

# money plots (outputs land next to the scripts' working directory)
cd evgen
python scripts/money_cos2phi.py                       # inclusive gluonometry
python scripts/money_cos2phi_coherent.py              # coherent intact-6Li tag
python scripts/money_tagged_azz.py --events 400000    # tagged tensor asymmetry
```

## Headline results so far

Every projection below is quoted at the full 10 fb⁻¹/u programme year in
its own configuration — its own spin states, its own far-forward optics,
its own isotope — so the reaches are alternatives and not a programme.
How a real year divides between them is unspecified by any EIC document
and ours to propose (priced in `plans/07` WP2); a share *f* multiplies
every statistical error by 1/√*f*, leaves a luminosity quoted as a reach
where it is, and multiplies the years to it by 1/*f*.  It is a flag
(`--run-share` in `fastsim/`, `--lumi-fraction` in `evgen/`), default 1.0.

- **Gluonometry**: with the sum-rule-constrained Δ model
  (`delta_models`, moment_A, table-α_s) the sweet-spot cos 2φ amplitudes are (0.44–0.95)×10⁻² against per-bin δA = (1.4–4.5)×10⁻⁴ already at the 1-year program (10 fb⁻¹/u), across Q² = 1.1–14 GeV² at x = 0.011–0.14 (e10 × ⁶Li 99.5 GeV/u, 2026-08-27) — the measurement resolves the interpretation-A-vs-B ansatz spread, not just a null.
- **Coherent intact-⁶Li tag**: the recoil is exactly beam-blind
  (A/Z = 2), so the only handle is its angle against the 10σ envelope and the pot aperture. With the Yellow Report divergences (a rectangular 10(σ_h, σ_v) envelope) and the measured ePIC aperture the tag does not survive at IP6 at any configuration (7×10⁻⁸ / 6×10⁻²⁷ / 7×10⁻¹⁴ at the envelope, plans/10, Report 4); it needs the lithium tagging optics of Report 1 §6.1 — the horizontal β* raised 46 / 164 / 89× with pots that follow the 0.36 / 0.19 / 0.12 mrad envelope, tagging 0.37 / 0.25 / 0.33 analytically at 1/6.8–1/12.8 of the luminosity, which the reconstructed chain returns as 2.3–6.0×10⁶ recoils per year (Report 4 Table 2) — or the IR-8 secondary focus. For a 0.20 GeV near-beam envelope (13.5%): 1.7×10⁷ / 1.7×10⁸ tagged events in the 1-/10-year programs; modulation amplitude anchored on the polarized-deuteron CGC calculation (sign flip vs d predicted).
- **Tagging inverts between isotopes at IP6**: ⁷Li α-tag works
  (97% to the Roman Pots at every configuration and optics in angle; at 5 × 41 the measured per-band insertion of the pot silicon leaves 57%, station 2 only, on a module edge — plans/00 run 14, `tools/fullsim/README.md`); ⁶Li α-tag is beam-blind — at the Yellow Report optics of each configuration (220/380, 180/180, 92/92 μrad h/v, plans/10) almost all of its near-beam tail is inside the 10σ envelope and what survives is the ≈ 1.5% that falls below R = 0.95 into the Roman-Pot window, an over-rigid branch on the inner side of the bend and 0.01–0.20 points of tail (1.9 / 1.7 / 2.6% in all at 5 × 41 / 10 × 100 / 18 × 275); the lithium tagging optics of Report 1 §6.1 recovers 31 / 22 / 29% at 1/6.8–1/12.8 of the luminosity (2026-08-29 on the per-configuration transport, `fastsim/scripts/tagging_acceptance.py`; the 1.85% quoted before applied one proton-derived 73 μrad at every configuration);
  the over-rigid ⁷Li triton (R = 1.29) is taggable at IP6 as well — it lands on the Roman-Pot inner side at every configuration and then in the ZDC, 78 / 92 / 94% at either optics (measured 2026-08-28, run 14; the triton bullet below) — so only the ⁶Li side needs the tagging optics or the IR-8 secondary focus.
- **The hadronic final state is PYTHIA 8, not a toy** (2026-08-26): 8 M
  events over the three beam configurations, generated natively
  (`tools/pythia8`).  Σ-method δy/y at the four sweet spots is 0.32 / 0.22 / 0.29 / 0.15 at the mid configuration (2026-08-27, corrected spots; library regenerated 2026-08-28 without PYTHIA's default m̂ ≥ 4 GeV floor, which had removed x < 16/s) — the toy was optimistic everywhere — and each configuration's own spots resolve to 0.39–0.11 (low) and 0.23–0.18 (top), so the x ≈ 0.14 bins belong to the **low-energy** configuration.
- **The Roman-Pot aperture for an intact ⁶Li is measured, and it is
  horizontal** (re-measured 2026-08-28 in the current ePIC geometry,
  `epic-main` git 9aaa2969): npsim cannot shoot a nucleus, so
  `tools/fullsim/ion_gun_hepmc.py` feeds one through the geometry as
  HepMC.  The boundary is |θ_x| ≳ 2.50 / 1.51 / 0.53 mrad in the
  5×41 / 10×100 / 18×275 optics against |θ_y| ≳ (shut) / 2.12 / 0.92 mrad
  — the opposite aspect to the slot the coherent chain assumes, which
  flips the sign of the acceptance-induced ⟨cos 2β⟩ (plans/04 #20).  The
  transport is measured with it: (R₁₂, R₃₄, D) = 19.24 / 4.56 / 0.311,
  21.25 / 3.35 / 0.287 and 29.97 / 2.93 / 0.292 m, a horizontal lever
  4.2× the vertical at 5 × 41 (R₃₄ there read off a zero-insertion scratch geometry on 2026-08-29, the real 29.6 mm insertion sitting at 6.49 mrad on it), 6.3× at 10 × 100 and 10.2× at 18 × 275, and the pot's outer
  edge at 2.85 / 3.85 / 4.00 mrad rather
  than the 5 mrad the acceptance tables assume.  Which of the two
  constraints binds is now per configuration: the silicon at 5 × 41
  (2.50 against a 2.20 mrad envelope), the machine at the other two, by
  eight orders of magnitude at 18 × 275 (0.53 against 0.92).  The 5 × 41
  row stands on the ePIC baseline proton lattice; the Z/A = 0.5 file at
  the ⁶Li fill point gives 1.61 mrad instead (48 mm / 29.81 m, the
  threshold the code applies; first hit at 1.60 on +x, 1.70 on −x), a
  ×0.64 systematic (plans/09 B1, D3).
- **The over-rigid ⁷Li triton is taggable** (2026-08-28): it lands on
  Roman-Pot silicon at dx = +66 mm in 60 of 60 events at every
  configuration, with a zero-degree-calorimeter deposit behind it in
  80–98%, so "no IP6 coverage" is withdrawn, `farforward.route_charged`
  has an over-rigid branch, the ⁷Li α + t double tag exists, and the t tag
  runs 78 / 92 / 94% against the 3.3 / 0.4 / 0.5% it was published at.
- **The coherent estimator is unbiased at any count** (2026-08-28): the
  bin-wise spin-state ratio carries an O(1/ν) offset that biased a_t by
  −34% in the sparsest |t| bin at one year; the acceptance-profiled
  likelihood (`reco.harmonic_likelihood_fit_2d`, `--fit likelihood`)
  removes it with the same errors and the same blind spots (plans/08 A12).
- **The α + d background has a second handle** (2026-08-28): given an α
  that fakes a coherent tag, its partner deuteron is in the pots 82–84% of
  the time at the tagging optics (12 / 2 / 25% at the Yellow Report
  optics), 11–17 mm away on the measured per-configuration levers and
  merged into one 500 μm pixel in at most
  4 × 10⁻³ of recorded pairs, none at all at the Yellow Report envelope
  (`nearbeam_two_hit.py`, plans/09 B4).  The tagged observables now run
  on the same optics: the ⁶Li α tag admits nothing below k = 0.15 GeV/c
  at any published optics and half its sample there at the tagging
  optics; for ⁷Li the tagging optics is a net loss.
- **Radiative corrections are bounded** (2026-08-28): collinear ISR
  leaves the mixed x exact and the azimuth invariant; the ISR-free
  corrections mis-state Δ̂ by 0.5–1.2% at the sweet spots (≤ 2.9% with
  the low-Q² feed-in), inside the 5% gate (`polligen/radiative.py`,
  plans/07 WP4).
- **First systematics numbers**: relative-luminosity bias formulas
  (removed exactly by lumi-corrected estimators); the α+d breakup
  background is controlled by the two-component |t| fit, with no
  event-by-event Z identification assumed, and by the partner-deuteron
  veto above where the background is made; Roman-Pot charge
  discrimination is a redirected question — the incumbent readout's
  charge dynamic range and gain suppression at ≈ 9 MIP, for the ePIC
  far-forward WG and the AC-LGAD group (plans/04 #19, plans/09 D1).

Status, caveats, and the full development-run log live in
[plans/00_README.md](plans/00_README.md); what is still missing from the
simulation chain, ordered, in
[plans/08_simulation_chain_completion.md](plans/08_simulation_chain_completion.md);
open external dependencies in
[plans/04_open_questions.md](plans/04_open_questions.md). Everything is
scenario-input-driven — bands, not predictions — with every external
number verified against primary sources.
