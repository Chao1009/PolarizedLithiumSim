# tools/pythia8 — hadronic final state for the reconstruction chain

`gen_dis_hfs.py` generates unpolarized e + p / e + n DIS events with
PYTHIA 8 and writes the scattered electron plus the final-state hadrons of
every event in the polligen HFS format (`evgen/polligen/hfs.py`,
`HFSSample`, one compressed `.npz`).  polligen then uses these events as a
*library* of hadronic final states: for each importance-sampled
pseudo-event at (x, Q²) it takes a library event from the same (x, Q²)
cell, applies the hadron-side detector response, and transfers the
captured fraction of Σ = Σ(E − p_z) and the p_T ratio onto the
pseudo-event's exact kinematics.  The tensor polarization enters only the
inclusive azimuthal weight, so an unpolarized generator is the right source
of the HFS (plans/05, tier T2).

## Where it runs

PYTHIA 8 with Python bindings is not available on the Windows analysis
machine (no compiler, no WSL, no container).  The eic-shell container on
the Linux box ships `import pythia8`; run there:

```bash
./eic-shell
cd PolarizedLithiumSim
python3 tools/pythia8/gen_dis_hfs.py --target p --n-events 1000000 \
    --electron-energy 10 --p-per-nucleon 50 --seed 1 --out pythia8_e10_p50_dis.npz
python3 tools/pythia8/gen_dis_hfs.py --target n --n-events 1000000 \
    --electron-energy 10 --p-per-nucleon 50 --seed 2 --out pythia8_e10_n50_dis.npz
```

Then copy the files into `evgen/samples/` (git-ignored if large; ≈ 200 MB
per million events) and run on any machine:

```bash
python evgen/scripts/hfs_resolution.py --sample evgen/samples/pythia8_e10_p50_dis.npz \
    evgen/samples/pythia8_e10_n50_dis.npz            # resolution figure + table
python evgen/scripts/money_cos2phi_reco.py --y-source hfs \
    --hfs-sample evgen/samples/pythia8_e10_p50_dis.npz evgen/samples/pythia8_e10_n50_dis.npz
```

Until the PYTHIA sample exists, both scripts fall back to the toy
string-fragmentation generator `hfs.ToyHFS` (flagged "toy" in every
output); its numbers are illustrative.

## Settings

Head-on frame (ion along +z at p_u GeV/u, electron along −z at E_e),
`WeakBosonExchange:ff2ff(t:gmZ)`, `SpaceShower:dipoleRecoil = on`,
`PhaseSpace:Q2Min = 0.7`, lepton QED radiation off (radiative effects are
plans/07 WP4).  6Li is represented per nucleon by a p/n mixture; nuclear
effects on the hadronic final state (Fermi motion, nuclear breakup —
BeAGLE, plans/02 step 1.5) are Phase-2 items.  The scattered electron is
taken as the most energetic final-state e⁻; x, Q², y are computed from the
beam and electron four-vectors (exact without radiation).
