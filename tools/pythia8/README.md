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

**And a second silent floor: `PhaseSpace:mHatMin`** (found by the
2026-08-28 code review).  PYTHIA's default lower limit on the invariant
mass of the hard 2 → 2 system is 4 GeV, and it applies to the DIS process
too, where m̂² = x s: nothing was generated below x = 16/s — 0.004 at
10 × 99.5, 0.020 at 5 × 40.8 — which is 39% of the selected rate of the
pseudo-experiments (the sweet spots at x ≥ 0.011 sit above it, the low-x
half of the Q² = 1.14 GeV² slice does not; the 1.8% of events that did
appear below the floor came from shower and primordial-k_T migration).
The loosened window needs m̂ ≥ √(Q²_min/y_max) = 0.84 GeV; the script now
sets `PhaseSpace:mHatMin = 0.5` (`--mhat-min`) and the standing production
was regenerated with it — the generated cross section rose from 0.666 to
0.947 μb at 10 × 99.5 (p), and the sample's x spectrum now follows the
generator's rate map to ±25% down to x = 3×10⁻⁴.

**The massless Σ identity does not hold for a massive target.**
`hfs.truth_kinematics_check` pins Σ over the final state = 2 E_e y, which
the toy generator satisfies exactly.  A PYTHIA sample does not: the target
nucleon enters with E − p_z ≈ m_N²/(2 p_u), which is 8.8 MeV at the
50 GeV/u of the diagnostic run above and 10.8 / 4.4 / 3.2 MeV at the
standing production's 40.8 / 99.5 / 137.5 GeV/u.  On that diagnostic run
Σ − 2 E_e y = 9.4 MeV (measured median), which is 0.14% of Σ at y = 1 but
2% at y = 0.01 and 7.7% at y = 0.004; the offset scales as 1/p_u, so it is
smaller in every file of the standing production except the low one.
Nothing downstream is affected — `HFSLibrary` transfers the *ratio*
Σ_reco/Σ_true of the library event's own sums, in which the offset cancels
— but the check itself must not be applied to a PYTHIA sample below
y ≈ 0.2.  The p_T identity is exact (≤ 5×10⁻¹¹).

## Production

The beam energies are the three γ-matched machine configurations,
`beams.default_configs("6Li")` = 5 × 40.8, 10 × 99.5 and 18 × 137.5 GeV/u
(plans/10 A0); the rigidity-scaled 5 × 20.5 and 10 × 50 this file carried
before 2026-08-27 are not machine configurations.  The mid one, doubled
because money plots 5R/7R are published at it:

```bash
export PYTHONPATH=$HOME/Apps/pythia8311/lib
export PYTHIA8DATA=$HOME/Apps/pythia8311/share/Pythia8/xmldoc
python3 tools/pythia8/gen_dis_hfs.py --target p --n-events 2000000 \
    --electron-energy 10 --p-per-nucleon 99.5 --seed 101 --quiet \
    --out evgen/samples/pythia8_e10_p99.5_dis.npz
python3 tools/pythia8/gen_dis_hfs.py --target n ... --seed 102 \
    --out evgen/samples/pythia8_e10_n99.5_dis.npz
```

Measured 2026-08-28 on eight cores, three jobs in parallel, with
`PhaseSpace:mHatMin = 0.5`: **7 000–8 500 events/s**, 455 MB and ~140 s per
million events at 10 × 99.5 (the 9 000–13 000 events/s, 330 MB and 115 s
quoted here before are the pre-2026-08-28 production, whose default m̂ floor
removed the low-x half of the sample and 27–36% of the particles with it).  The six files
of the standing production and their cross sections are the manifest in
[evgen/samples/README.md](../../evgen/samples/README.md).

Then, on any machine that has the files:

```bash
python3 evgen/scripts/hfs_resolution.py --config 1 \
    --sample evgen/samples/pythia8_e10_p99.5_dis.npz \
             evgen/samples/pythia8_e10_n99.5_dis.npz   # resolution figure + table
python3 evgen/scripts/money_cos2phi_reco.py --y-source hfs \
    --hfs-sample evgen/samples/pythia8_e10_p99.5_dis.npz \
                 evgen/samples/pythia8_e10_n99.5_dis.npz
```

Without a sample both scripts fall back to the toy string-fragmentation
generator `hfs.ToyHFS` (flagged "toy" in every output); its numbers are
illustrative.

The p and n files are merged by `HFSSample.concatenate`, which since
2026-08-28 weights each file by σ_gen/n_events (plans/08 A10).  That
matters because the two cross sections are not equal — 0.947 vs 0.855 μb
at 10 × 99.5, 0.642 vs 0.555 at 5 × 40.8, 1.164 vs 1.070 at 18 × 137.5 —
so the equal-count merge the manifest still generates would have mixed the
targets 1 : 1 where ⁶Li (Z = N = 3) asks for σ_p : σ_n, 1.11 : 1 at the mid
configuration.  With the weighting in place the event counts of the two
files no longer have to match.
