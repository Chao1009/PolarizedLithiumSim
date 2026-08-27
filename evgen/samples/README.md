# PYTHIA 8 HFS samples

Git-ignored `.npz` files in the polligen HFS format (`polligen.hfs`,
`HFSSample`); regenerate with `tools/pythia8/gen_dis_hfs.py`, whose README
carries the build recipe and the two generator subtleties that matter
(`PhaseSpace:Q2Min` needs `pTHatMinDiverge = 0.5` to have any effect; the
massless Σ = 2 E_e y identity does not hold for a massive target).

## Standing production — 2026-08-26, PYTHIA 8.311, `Q2Min = 0.7 GeV²`

One p and one n file per beam configuration; ⁶Li (Z = N = 3) is the pair
merged by event count.  Configurations are `beams.default_configs("6Li")`,
i.e. `--config 0/1/2` of the consuming scripts.

| file | config | events | particles | σ_gen [μb] | seed | size |
|---|---|---|---|---|---|---|
| `pythia8_e5_p20.5_dis.npz` | 0 low, 5 × 20.5 | 1 000 000 | 6 026 571 | 0.3255 | 201 | 246 MB |
| `pythia8_e5_n20.5_dis.npz` | 0 low, 5 × 20.5 | 1 000 000 | 6 482 737 | 0.2514 | 202 | 260 MB |
| `pythia8_e10_p50_dis.npz` | 1 mid, 10 × 50 | 2 000 000 | 17 464 095 | 0.5527 | 101 | 662 MB |
| `pythia8_e10_n50_dis.npz` | 1 mid, 10 × 50 | 2 000 000 | 18 719 553 | 0.4659 | 102 | 702 MB |
| `pythia8_e18_p137.5_dis.npz` | 2 top, 18 × 137.5 | 1 000 000 | 11 762 519 | 0.8367 | 301 | 426 MB |
| `pythia8_e18_n137.5_dis.npz` | 2 top, 18 × 137.5 | 1 000 000 | 12 477 453 | 0.7448 | 302 | 448 MB |

Total 2.7 GB, 12 minutes of wall time with three jobs in parallel on eight
cores.  The mid configuration is doubled because it is the one money plots
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
