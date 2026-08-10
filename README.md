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
| `evgen/` | `polligen`: the doubly polarized e+⁶,⁷Li event generator (the first of its kind) — spin-density kernel, inclusive + tagged samplers, run-plan bookkeeping, estimators, the coherent intact-⁶Li channel, money plots 4–6 |
| `docs/` | source documents + the self-contained physics note [docs/note_cos2phi_coherent_6Li.md](docs/note_cos2phi_coherent_6Li.md) (verified 50-entry bibliography) |
| `reports/` | circulate-able reports (HTML source + rendered PDF), currently the cos 2φ money-plot report |
| `tools/` | BeAGLE and full-simulation (eic-shell/ePIC) setup notes |

## Quick start

```bash
# fast simulation (15 tests; PDF-grid tests need the `parton` package)
cd fastsim && python -m pytest tests/ -q

# event generator (66 tests)
cd evgen && python -m pytest tests/ -q

# money plots (outputs land next to the scripts' working directory)
cd evgen
python scripts/money_cos2phi.py --lumi 100            # inclusive gluonometry
python scripts/money_cos2phi_coherent.py --lumi 100   # coherent intact-6Li tag
python scripts/money_tagged_azz.py --events 400000    # tagged tensor asymmetry
```

## Headline results so far

- **Gluonometry reach**: 5σ on Δ/F₁ = 10⁻³ within ~15–40 fb⁻¹/u
  (combined bins); projected-data money plots with per-bin
  δA ≈ (0.4–1.2)×10⁻⁴ at 100 fb⁻¹/u across Q² = 1.1–14 GeV².
- **Coherent intact-⁶Li tag**: the recoil is exactly beam-blind
  (A/Z = 2) — 13.5% acceptance via the Roman-Pot pT tail with
  high-acceptance optics, ~10⁻⁴ with high-divergence: the coherent
  program fixes the optics choice. ~10⁸ tagged events at 100 fb⁻¹/u;
  modulation amplitude anchored on the polarized-deuteron CGC
  calculation (sign flip vs d predicted).
- **Tagging inverts between isotopes at IP6**: ⁷Li α-tag works
  (96–99% to the Roman Pots); ⁶Li α-tag is beam-blind (3–9%);
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
