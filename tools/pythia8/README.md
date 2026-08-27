# tools/pythia8 — hadronic final state for the reconstruction chain

*Reproducing a result rather than changing one?  Start at
[docs/reproduction_manual.md](../../docs/reproduction_manual.md) §5.1, which
compresses this file into the commands and their expected outputs.*

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

## Where it runs — natively, since 2026-08-26

PYTHIA 8 no longer needs a container.  The Python bindings build from the
PYTHIA source against the analysis machine's own interpreter, which is what
makes the samples reproducible next to the code that consumes them:

```bash
curl -sSL -o pythia8311.tar.gz \
  https://gitlab.com/Pythia8/releases/-/archive/pythia8311/releases-pythia8311.tar.gz
tar xzf pythia8311.tar.gz && cd releases-pythia8311
./configure --prefix=$HOME/Apps/pythia8311 \
            --with-python-config=$(which python3-config) --with-gzip
make -j8 && make install
export PYTHONPATH=$HOME/Apps/pythia8311/lib
export PYTHIA8DATA=$HOME/Apps/pythia8311/share/Pythia8/xmldoc
python3 -c "import pythia8; print(pythia8.Pythia('', False).settings.parm('Pythia:versionNumber'))"
```

pybind11 ships inside the release (`plugins/python/include/`), so the only
prerequisites are a C++ compiler and the Python development headers; the
whole build — configure, library, bindings, install — took **2.5 minutes**
on eight cores (14:44:00 to 14:46:27, 2026-08-26).  The eic-shell container
(`jug_xl-nightly.sif`) carries the C++ library and headers but **no Python
bindings** — `import pythia8` fails there, which is why the previous
version of this file routed everything through it.

`infoPython()` in the 8.3 bindings returns the `Info` object **by value**:
a copy taken before the event loop still reports `sigmaGen() = 0` after a
million events, and re-fetching it per event is expensive and, repeated
enough times, unstable.  The writer therefore fetches it once at the end
for the cross section, and checks once (after the first event) that the run
is unweighted rather than assuming it.

## Settings, and the cut that was not being applied

Head-on frame (ion along +z at p_u GeV/u, electron along −z at E_e),
`WeakBosonExchange:ff2ff(t:gmZ)`, `SpaceShower:dipoleRecoil = on`,
`PhaseSpace:Q2Min = 0.7`, lepton QED radiation off (radiative effects are
plans/07 WP4).  ⁶Li is represented per nucleon by a p/n mixture; nuclear
effects on the hadronic final state (Fermi motion, nuclear breakup —
BeAGLE, plans/02 step 1.5) are Phase-2 items.  The scattered electron is
taken as the most energetic final-state e⁻; x, Q², y are computed from the
beam and electron four-vectors (exact without radiation).

**`PhaseSpace:Q2Min` alone does nothing.**  PYTHIA applies it only when
`Q2Min ≥ pTHatMinDiverge²` (`PhaseSpace.cc`: `hasQ2Min = (Q2GlobalMin >=
pow2(pTHatMinDiverge))`), and `pTHatMinDiverge` defaults to 1 GeV — so a
requested 0.7 GeV² was silently ignored and the divergence cut set the
floor instead.  Measured at 10 × 50 GeV, 4000 events (a diagnostic run, not a
production configuration): with the default,
min Q² = 1.002 GeV² and σ = 0.381 μb *whatever* `Q2Min` asks for below 1;
with `PhaseSpace:pTHatMinDiverge = 0.5` (PYTHIA's own lower limit for that
parameter, now set by the script), min Q² = 0.697 GeV² and σ = 0.551 μb.
The 0.7–1.0 GeV² band is 31% of the sample and is exactly the band the
reconstructed-level pseudo-experiments loosen their generator window to,
so that events can migrate *into* the analysis bins — it was missing
entirely from every earlier sample plan.

**The massless Σ identity does not hold for a massive target.**
`hfs.truth_kinematics_check` pins Σ over the final state = 2 E_e y, which
the toy generator satisfies exactly.  A PYTHIA sample does not: the target
nucleon carries E − p_z = 8.8 MeV at 50 GeV/u, so Σ − 2 E_e y = 9.4 MeV
(measured median), which is 0.14% of Σ at y = 1 but 2% at y = 0.01 and
7.7% at y = 0.004.  Nothing downstream is affected — `HFSLibrary`
transfers the *ratio* Σ_reco/Σ_true of the library event's own sums, in
which the offset cancels — but the check itself must not be applied to a
PYTHIA sample below y ≈ 0.2.  The p_T identity is exact (≤ 5×10⁻¹¹).

## Production

```bash
export PYTHONPATH=$HOME/Apps/pythia8311/lib
export PYTHIA8DATA=$HOME/Apps/pythia8311/share/Pythia8/xmldoc
python3 tools/pythia8/gen_dis_hfs.py --target p --n-events 2000000 \
    --electron-energy 10 --p-per-nucleon 50 --seed 101 --quiet \
    --out evgen/samples/pythia8_e10_p50_dis.npz
python3 tools/pythia8/gen_dis_hfs.py --target n ... --seed 102 \
    --out evgen/samples/pythia8_e10_n50_dis.npz
```

Measured 2026-08-26 on eight cores, three jobs in parallel: **9 000–13 000
events/s**, 330 MB and ~115 s per million events at 10 × 99.5.  The six files
of the standing production and their cross sections are the manifest in
[evgen/samples/README.md](../../evgen/samples/README.md).

Then, on any machine that has the files:

```bash
python3 evgen/scripts/hfs_resolution.py --config 1 \
    --sample evgen/samples/pythia8_e10_p50_dis.npz \
             evgen/samples/pythia8_e10_n50_dis.npz     # resolution figure + table
python3 evgen/scripts/money_cos2phi_reco.py --y-source hfs \
    --hfs-sample evgen/samples/pythia8_e10_p50_dis.npz \
                 evgen/samples/pythia8_e10_n50_dis.npz
```

Without a sample both scripts fall back to the toy string-fragmentation
generator `hfs.ToyHFS` (flagged "toy" in every output); its numbers are
illustrative.

The p and n files are merged by event count (`HFSSample.concatenate`),
which is what Z = N = 3 asks for in ⁶Li — but the two cross sections are
not equal (0.666 vs 0.576 μb at 10 × 99.5, 0.420 vs 0.339 at 5 × 40.8), so
generating the same number of events for both is a *choice of weighting*,
not a neutral merge.  It is the right one for a library that supplies
hadronic shapes per (x, Q²) cell and the wrong one for anything that reads
the sample as a rate.
