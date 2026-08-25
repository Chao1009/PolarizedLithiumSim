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
| `evgen/` | `polligen`: the doubly polarized e+⁶,⁷Li event generator (the first of its kind) — spin-density kernel, inclusive + tagged samplers, run-plan bookkeeping, estimators, the coherent intact-⁶Li channel, money plots 4–7 and their reconstructed-level versions 5R/7R/6R |
| `docs/` | source documents + the self-contained physics note [docs/note_cos2phi_coherent_6Li.md](docs/note_cos2phi_coherent_6Li.md) (verified 50-entry bibliography) + the simulation code review and reconstruction audit [docs/code_review_2026-08-25.md](docs/code_review_2026-08-25.md) (measurability audit of 5R/7R/6R, findings F1-F13, numerical checks) |
| `reports/` | circulate-able reports (self-contained HTML + rendered PDF): the cos 2φ money-plot report, the educational primer, and the reconstruction-chain analysis note (what is measured, how the azimuth and Δ are reconstructed, audit of the simulation — 2026-08-24); served as a website by the GitHub Pages workflow (`.github/workflows/pages.yml` + `reports/index.html`) — activate once via Settings → Pages → Source: "GitHub Actions" (public repo required on free plans) |
| `tools/` | BeAGLE and full-simulation (eic-shell/ePIC) setup notes |

## Quick start

```bash
# fast simulation (24 tests; PDF-grid tests need the `parton` package)
cd fastsim && python -m pytest tests/ -q

# event generator (94 tests)
cd evgen && python -m pytest tests/ -q

# money plots (outputs land next to the scripts' working directory)
cd evgen
python scripts/money_cos2phi.py                       # inclusive gluonometry
python scripts/money_cos2phi_coherent.py              # coherent intact-6Li tag
python scripts/money_tagged_azz.py --events 400000    # tagged tensor asymmetry
```

## Headline results so far

- **Gluonometry**: with the sum-rule-constrained Δ model
  (`delta_models`, moment_A, table-α_s) the sweet-spot cos 2φ
  amplitudes are (0.7–1.2)×10⁻² against per-bin δA = (1.5–4.5)×10⁻⁴
  already at the 1-year program (10 fb⁻¹/u), across Q² = 1.1–14 GeV²
  at x ≈ 0.02–0.14 — the measurement resolves the interpretation-A-vs-B
  ansatz spread, not just a null.
- **Coherent intact-⁶Li tag**: the recoil is exactly beam-blind
  (A/Z = 2) — 13.5% acceptance via the Roman-Pot pT tail with
  high-acceptance optics, 4×10⁻⁵ with high-divergence: the coherent
  program fixes the optics choice. 1.1×10⁷ / 1.1×10⁸ tagged events in
  the 1-/10-year programs; modulation amplitude anchored on the
  polarized-deuteron CGC calculation (sign flip vs d predicted).
- **Tagging inverts between isotopes at IP6**: ⁷Li α-tag works
  (96–99% to the Roman Pots); ⁶Li α-tag is beam-blind (3–9% with the constant 0.20 GeV proton cut; ≈ 1–3% once the 10σ envelope is scaled to the α momentum — code review 2026-08-25);
  tritons need IR-8.
- **First systematics numbers**: relative-luminosity bias formulas
  (removed exactly by lumi-corrected estimators); the α+d breakup
  background requires Roman-Pot charge discrimination — an open
  hardware question posed to the ePIC far-forward WG
  (plans/04 #19).

Status, caveats, and the full development-run log live in
[plans/00_README.md](plans/00_README.md); open external dependencies in
[plans/04_open_questions.md](plans/04_open_questions.md). Everything is
scenario-input-driven — bands, not predictions — with every external
number verified against primary sources.
