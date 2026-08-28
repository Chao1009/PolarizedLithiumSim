# PYTHIA 8 HFS samples

Git-ignored `.npz` files in the polligen HFS format (`polligen.hfs`,
`HFSSample`); regenerate with `tools/pythia8/gen_dis_hfs.py`, whose README
carries the build recipe and the three generator subtleties that matter
(`PhaseSpace:Q2Min` needs `pTHatMinDiverge = 0.5` to have any effect;
`PhaseSpace:mHatMin` must be lowered from its 4 GeV default or nothing is
generated below x = 16/s; the massless Σ = 2 E_e y identity does not hold
for a massive target).

## Standing production — 2026-08-28, PYTHIA 8.311, `Q2Min = 0.7 GeV²`, `mHatMin = 0.5 GeV`

One p and one n file per beam configuration; ⁶Li (Z = N = 3) is the pair
merged by `HFSSample.concatenate`, which weights the events of each file by
σ_gen/n_events so that the targets enter in the ratio of their cross
sections.  Configurations are `beams.default_configs("6Li")`, i.e.
`--config 0/1/2` of the consuming scripts.

| file | config | events | particles | σ_gen [μb] | seed | size |
|---|---|---|---|---|---|---|
| `pythia8_e5_p40.8_dis.npz` | 0 low, 5 × 40.8 | 1 000 000 | 9 698 916 | 0.6423 | 201 | 362 MB |
| `pythia8_e5_n40.8_dis.npz` | 0 low, 5 × 40.8 | 1 000 000 | 10 383 046 | 0.5551 | 202 | 383 MB |
| `pythia8_e10_p99.5_dis.npz` | 1 mid, 10 × 99.5 | 2 000 000 | 25 314 441 | 0.9473 | 101 | 910 MB |
| `pythia8_e10_n99.5_dis.npz` | 1 mid, 10 × 99.5 | 2 000 000 | 26 749 460 | 0.8551 | 102 | 955 MB |
| `pythia8_e18_p137.5_dis.npz` | 2 top, 18 × 137.5 | 1 000 000 | 14 482 707 | 1.164 | 301 | 512 MB |
| `pythia8_e18_n137.5_dis.npz` | 2 top, 18 × 137.5 | 1 000 000 | 15 225 574 | 1.070 | 302 | 536 MB |

**8 000 000 events, 101.8 M final-state particles, 3.7 GB**, in five minutes
of wall time with six jobs in parallel on eight cores (117–283 s).

**Regenerated 2026-08-28** with `PhaseSpace:mHatMin = 0.5`.  The previous
production (2026-08-27, same seeds and event counts) had PYTHIA's default
m̂ ≥ 4 GeV, which for DIS is x ≥ 16/s: nothing was generated below
x = 0.004 at 10 × 99.5 (0.020 at 5 × 40.8, 0.0016 at 18 × 137.5) except a
1.8% tail from shower migration, while 39% of the selected rate of the
reconstructed-level pseudo-experiments lies there.  The cross sections
rose by 37–48% and the particle multiplicities by 27–36% with the low-x
region in; the sweet spots (x ≥ 0.011) were not affected, the low-x half of
the Q² = 1.14 GeV² slice was.  The mid configuration is doubled because it
is the one money plots 5R/7R are published at.

Reconstructed Q² reaches below the 0.7 GeV² generator cut in every file
(down to 0.58 GeV²): the cut is on the hard process's −t̂, while the
recorded Q² is rebuilt from the final scattered electron, which the shower
and the primordial k_T move.  It is a physical tail of a few events per
million, not a violated cut.

## What consumes them

```bash
python3 evgen/scripts/hfs_resolution.py --config 1 \
    --sample evgen/samples/pythia8_e10_p99.5_dis.npz \
             evgen/samples/pythia8_e10_n99.5_dis.npz
python3 evgen/scripts/hfs_acceptance.py --config 1 --sample <the same pair>
python3 evgen/scripts/money_cos2phi_reco.py --y-source hfs --hfs-sample <the same pair> \
    --hfs-calibrate --unfold folded --tag _hfscal
```

`hfs_resolution.py` and `money_cos2phi_reco.py` fall back to `hfs.ToyHFS`
and label their output "toy" when no sample is given.
