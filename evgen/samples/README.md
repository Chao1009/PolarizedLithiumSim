# PYTHIA 8 HFS samples

Git-ignored `.npz` files in the polligen HFS format (`polligen.hfs`,
`HFSSample`); regenerate with `tools/pythia8/gen_dis_hfs.py`, whose README
carries the build recipe and the two generator subtleties that matter
(`PhaseSpace:Q2Min` needs `pTHatMinDiverge = 0.5` to have any effect; the
massless Σ = 2 E_e y identity does not hold for a massive target).

## Standing production — 2026-08-27, PYTHIA 8.311, `Q2Min = 0.7 GeV²`

One p and one n file per beam configuration; ⁶Li (Z = N = 3) is the pair
merged by event count.  Configurations are `beams.default_configs("6Li")`,
i.e. `--config 0/1/2` of the consuming scripts.

| file | config | events | particles | σ_gen [μb] | seed | size |
|---|---|---|---|---|---|---|
| `pythia8_e5_p40.8_dis.npz` | 0 low, 5 × 40.8 | 1 000 000 | 7 169 605 | 0.4202 | 201 | 285 MB |
| `pythia8_e5_n40.8_dis.npz` | 0 low, 5 × 40.8 | 1 000 000 | 7 715 892 | 0.3390 | 202 | 299 MB |
| `pythia8_e10_p99.5_dis.npz` | 1 mid, 10 × 99.5 | 2 000 000 | 19 980 470 | 0.6661 | 101 | 558 MB |
| `pythia8_e10_n99.5_dis.npz` | 1 mid, 10 × 99.5 | 2 000 000 | 21 346 737 | 0.5761 | 102 | 312 MB |
| `pythia8_e18_p137.5_dis.npz` | 2 top, 18 × 137.5 | 1 000 000 | 11 762 519 | 0.8367 | 301 | 426 MB |
| `pythia8_e18_n137.5_dis.npz` | 2 top, 18 × 137.5 | 1 000 000 | 12 477 453 | 0.7448 | 302 | 448 MB |

**8 000 000 events, 80.4 M final-state particles, 2.8 GB**, in 4 minutes of
wall time with four jobs in parallel on eight cores (95 / 101 / 229 / 240 s).

**Regenerated 2026-08-27** at the corrected beam energies. The two lower
configurations were produced at 20.5 and 50 GeV/u until then, which are
rigidity-scaled and are not machine configurations: EIC ions are
gamma-matched, so 6Li sits at 40.8 and 99.5 GeV/u (plans/10). The top
configuration is rigidity-capped and did not move, which is why its files
and every number derived from them are unchanged -- a useful check that
nothing else drifted.  The mid configuration is doubled because it is the one money plots
5R/7R are published at; the low configuration is what the reconstruction
note lists as not yet run for the x ≈ 0.1 bins, and the top configuration
is where the coherent channel dies and the inclusive one does not.

Reconstructed Q² reaches below the 0.7 GeV² generator cut in the two
largest files (0.166 and 0.560 GeV²): the cut is on the hard process's
−t̂, while the recorded Q² is rebuilt from the final scattered electron,
which the shower and the primordial k_T move.  It is a physical tail of a
few events per million, not a violated cut.

## What consumes them

```bash
python3 evgen/scripts/hfs_resolution.py --config 1 \
    --sample evgen/samples/pythia8_e10_p50_dis.npz \
             evgen/samples/pythia8_e10_n50_dis.npz
python3 evgen/scripts/money_cos2phi_reco.py --y-source hfs --hfs-sample <the same pair>
```

Both fall back to `hfs.ToyHFS` and label their output "toy" when no sample
is given.
