#!/usr/bin/env python3
"""Generate unpolarized e + N deep-inelastic events with PYTHIA 8 and write
the hadronic final state in the polligen HFS sample format
(evgen/polligen/hfs.py: HFSSample, compressed .npz).

Needs the PYTHIA 8 Python bindings (`import pythia8`).  They are built from
the PYTHIA source with `./configure --with-python-config=<python3-config>`
(see tools/pythia8/README.md); the eic-shell container ships the C++ library
only.  polligen then uses these events as a *library* of hadronic final
states: for each importance-sampled pseudo-event at (x, Q2) it takes a
library event from the same (x, Q2) cell, applies the hadron-side detector
response, and transfers the captured fraction of Sigma = Sigma(E - p_z) and
the p_T ratio onto the pseudo-event's exact kinematics.  The tensor
polarization enters only the inclusive azimuthal weight, so an unpolarized
generator is the right source of the HFS (plans/05, tier T2).

Physics setup (PYTHIA 8.3 DIS, cf. examples/main36 and the EIC tutorials):
  Beams:frameType = 2 with idA = 2212 (proton) or 2112 (neutron) along +z at
  eA = p_u GeV/u and idB = 11 (electron) along -z at eB = E_e -- the head-on
  frame of polligen; WeakBosonExchange:ff2ff(t:gmZ) = on (photon/Z
  exchange), SpaceShower:dipoleRecoil = on (the recommended DIS setting),
  PhaseSpace:Q2Min = 0.7 GeV^2 (the loosened generator window of the
  reconstructed-level pseudo-experiments), lepton QED radiation switched
  off (PDF:lepton = off, TimeShower:QEDshowerByL = off) so that the
  electron-method kinematics of the record are exact and radiative effects
  stay a separate work package (plans/07 WP4).  6Li is represented per
  nucleon by an equal mixture of e+p and e+n runs (nuclear effects on the
  HFS are a Phase-2 item).

  PhaseSpace:mHatMin = 0.5 (--mhat-min) removes PYTHIA's default 4 GeV
  floor on the hard-process invariant mass, which for DIS is mHat^2 = x s
  and had silently removed everything below x = 16/s (2026-08-28; the
  standing production of evgen/samples was regenerated with it).

  PhaseSpace:pTHatMinDiverge = 0.5 is what MAKES the Q2 cut real, and it
  has to be set explicitly.  PYTHIA applies PhaseSpace:Q2Min only when
  Q2Min >= pTHatMinDiverge^2 (PhaseSpace.cc: hasQ2Min = (Q2GlobalMin >=
  pow2(pTHatMinDiverge))), and pTHatMinDiverge defaults to 1 GeV -- so
  with the default the requested 0.7 GeV^2 is silently ignored and the
  divergence cut itself sets the floor.  Measured at 10 x 50 GeV -- a
  diagnostic run, not a machine configuration -- with 4000
  events: default -> min Q2 = 1.002 GeV^2 and sigma = 0.381 ub, whatever
  Q2Min asks for below 1; with pTHatMinDiverge = 0.5 -> min Q2 = 0.697
  GeV^2 and sigma = 0.551 ub, i.e. the 0.7-1.0 GeV^2 band that the
  loosened window exists to populate is 31% of the sample and was missing
  entirely (2026-08-26; 0.5 is PYTHIA's own lower limit for the parameter).

Neutrinos are KEPT in the record: the exact relation the sample is tested
against, Sigma over the whole final state = 2 E_e, holds only with them in
(hfs.truth_kinematics_check).  The detector response drops them
(hfs.HadronResponse, NEUTRINOS).

The beam energies are the gamma-matched machine configurations of
`beams.default_configs("6Li")` -- 5 x 40.8, 10 x 99.5 and 18 x 137.5 GeV/u
(plans/10 A0).  The rigidity-scaled 5 x 20.5 and 10 x 50 that this script
defaulted to before 2026-08-27 are not machine configurations and the
samples named after them no longer exist.

Usage (the mid configuration of the standing production; the manifest with
all six files and their seeds is evgen/samples/README.md):
  export PYTHONPATH=$HOME/Apps/pythia8311/lib
  export PYTHIA8DATA=$HOME/Apps/pythia8311/share/Pythia8/xmldoc
  python3 tools/pythia8/gen_dis_hfs.py --target p --n-events 2000000 \
      --electron-energy 10 --p-per-nucleon 99.5 --seed 101 \
      --out evgen/samples/pythia8_e10_p99.5_dis.npz
  python3 tools/pythia8/gen_dis_hfs.py --target n ... --seed 102 \
      --out evgen/samples/pythia8_e10_n99.5_dis.npz
The p and n files are merged by the consuming scripts
(HFSSample.concatenate), which since 2026-08-28 weights each file by
sigma_gen/n_events so the two targets enter in the ratio of their cross
sections rather than of their event counts -- for 6Li (Z = N = 3) that is
the p : n mix the nucleus asks for, and it no longer matters whether the
two files hold the same number of events.
"""

import argparse
import json
import sys
import time

import numpy as np

try:
    import pythia8
except ImportError:  # pragma: no cover
    sys.exit("pythia8 python module not found -- set PYTHONPATH to the "
             "directory holding pythia8.so (tools/pythia8/README.md)")


def _info(pythia):
    """The Info object.  The 8.3 python plugin exposes it as a method
    (infoPython()); older bindings expose the C++ member directly."""
    return pythia.infoPython() if hasattr(pythia, "infoPython") else pythia.info


class _Buffer:
    """Flat per-particle and per-event accumulators that convert to numpy in
    chunks -- a million events hold ~3x10^7 particles, which is several GB
    as Python lists but under a GB as float64 arrays."""

    def __init__(self, chunk=20000):
        self.chunk = chunk
        self._pid, self._charge, self._p4 = [], [], []
        self._counts = []
        self._ev = {k: [] for k in ("x", "q2", "y", "weight")}
        self._kp = []
        self.parts = {"pid": [], "charge": [], "p4": []}
        self.evs = {k: [] for k in ("x", "q2", "y", "weight", "kp", "counts")}

    def add(self, pid, charge, p4, x, q2, y, kp, weight):
        self._pid.extend(pid)
        self._charge.extend(charge)
        self._p4.extend(p4)
        self._counts.append(len(pid))
        self._ev["x"].append(x)
        self._ev["q2"].append(q2)
        self._ev["y"].append(y)
        self._ev["weight"].append(weight)
        self._kp.append(kp)
        if len(self._counts) >= self.chunk:
            self.flush()

    def flush(self):
        if not self._counts:
            return
        self.parts["pid"].append(np.array(self._pid, dtype=np.int64))
        self.parts["charge"].append(np.array(self._charge, dtype=float))
        self.parts["p4"].append(np.array(self._p4, dtype=float).reshape(-1, 4))
        self.evs["counts"].append(np.array(self._counts, dtype=np.int64))
        for k in ("x", "q2", "y", "weight"):
            self.evs[k].append(np.array(self._ev[k], dtype=float))
        self.evs["kp"].append(np.array(self._kp, dtype=float).reshape(-1, 4))
        self._pid, self._charge, self._p4 = [], [], []
        self._counts = []
        self._ev = {k: [] for k in ("x", "q2", "y", "weight")}
        self._kp = []

    def arrays(self):
        self.flush()
        cat = lambda d, k: (np.concatenate(d[k]) if d[k]
                            else np.zeros(0, dtype=float))
        counts = cat(self.evs, "counts")
        offsets = np.zeros(counts.size + 1, dtype=np.int64)
        np.cumsum(counts, out=offsets[1:])
        return {
            "offsets": offsets,
            "pid": cat(self.parts, "pid").astype(np.int64),
            "charge": cat(self.parts, "charge"),
            "p4": (np.concatenate(self.parts["p4"]) if self.parts["p4"]
                   else np.zeros((0, 4))),
            "x": cat(self.evs, "x"), "q2": cat(self.evs, "q2"),
            "y": cat(self.evs, "y"), "weight": cat(self.evs, "weight"),
            "kp": (np.concatenate(self.evs["kp"]) if self.evs["kp"]
                   else np.zeros((0, 4))),
        }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="p", choices=("p", "n"))
    ap.add_argument("--n-events", type=int, default=100000)
    ap.add_argument("--electron-energy", type=float, default=10.0)
    ap.add_argument("--p-per-nucleon", type=float, default=99.5,
                    help="GeV/u; the standing production runs 5 x 40.8, "
                         "10 x 99.5 and 18 x 137.5, i.e. "
                         "beams.default_configs('6Li') -- the default is the "
                         "mid one, which pairs with --electron-energy 10")
    ap.add_argument("--q2min", type=float, default=0.7)
    ap.add_argument("--mhat-min", type=float, default=0.5,
                    help="PhaseSpace:mHatMin [GeV], the invariant mass of the "
                         "hard 2 -> 2 system; for DIS mHat^2 = x s, so PYTHIA's "
                         "default of 4 GeV silently removes x < 16/s (2026-08-28)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="pythia8_dis_hfs.npz")
    ap.add_argument("--report-every", type=int, default=50000)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    p = pythia8.Pythia("", not args.quiet)
    settings = [
        "Beams:frameType = 2",
        "Beams:idA = %d" % (2212 if args.target == "p" else 2112),
        "Beams:idB = 11",
        "Beams:eA = %g" % args.p_per_nucleon,
        "Beams:eB = %g" % args.electron_energy,
        "WeakBosonExchange:ff2ff(t:gmZ) = on",
        "PhaseSpace:Q2Min = %g" % args.q2min,
        # without this the line above does nothing: the cut is applied only
        # for Q2Min >= pTHatMinDiverge^2, and the default is 1 GeV
        "PhaseSpace:pTHatMinDiverge = 0.5",
        # and this is the SECOND silent floor: mHatMin defaults to 4 GeV for
        # every 2 -> 2 process, DIS included, where mHat^2 = x s -- so with
        # the default nothing is generated below x = 16/s (0.004 at
        # 10 x 99.5), 39% of the selected rate of the pseudo-experiments
        # (code review 2026-08-28).  The loosened window needs
        # mHat >= sqrt(Q2min / y_max) = 0.84 GeV; 0.5 leaves a margin.
        "PhaseSpace:mHatMin = %g" % args.mhat_min,
        "SpaceShower:dipoleRecoil = on",
        "PDF:lepton = off",
        "TimeShower:QEDshowerByL = off",
        "Random:setSeed = on",
        "Random:seed = %d" % args.seed,
        # per-event printing off: a million events must not print a million
        # listings (the three counters default to 1, not 0)
        "Next:numberCount = 0",
        "Next:numberShowInfo = 0",
        "Next:numberShowProcess = 0",
        "Next:numberShowEvent = 0",
    ]
    if args.quiet:
        settings.append("Print:quiet = on")
    for s in settings:
        p.readString(s)
    if args.q2min < 0.25:
        sys.exit("--q2min below 0.25 GeV^2 cannot be enforced: PYTHIA needs "
                 "Q2Min >= pTHatMinDiverge^2 and 0.5 GeV is its lower limit "
                 "for that parameter")
    if not p.init():
        sys.exit("pythia init failed")

    buf = _Buffer()
    t0 = time.time()
    n_done = 0
    n_tried = 0
    while n_done < args.n_events:
        n_tried += 1
        if not p.next():
            continue
        ev = p.event
        # beams: record 1 = A (nucleon, +z), record 2 = B (electron, -z)
        k = np.array([ev[2].e(), ev[2].px(), ev[2].py(), ev[2].pz()])
        P = np.array([ev[1].e(), ev[1].px(), ev[1].py(), ev[1].pz()])
        pid, charge, p4 = [], [], []
        best, ebest = -1, -1.0
        for i in range(ev.size()):
            pt = ev[i]
            if not pt.isFinal():
                continue
            e = pt.e()
            if pt.id() == 11 and e > ebest:
                if best >= 0:                      # demote the previous best
                    pid.append(11)
                    charge.append(-1.0)
                    p4.append(kp_best)
                best, ebest = i, e
                kp_best = (e, pt.px(), pt.py(), pt.pz())
                continue
            pid.append(pt.id())
            charge.append(pt.charge())
            p4.append((e, pt.px(), pt.py(), pt.pz()))
        if best < 0:
            continue
        kp = np.array(kp_best)
        q = k - kp
        q2 = -(q[0] ** 2 - q[1] ** 2 - q[2] ** 2 - q[3] ** 2)
        pq = P[0] * q[0] - P[1] * q[1] - P[2] * q[2] - P[3] * q[3]
        pk = P[0] * k[0] - P[1] * k[1] - P[2] * k[2] - P[3] * k[3]
        y = pq / pk
        x = q2 / (2.0 * pq)
        buf.add(pid, charge, p4, x, q2, y, kp_best, 1.0)
        n_done += 1
        if n_done == 1:
            # Info is returned BY VALUE by the 8.3 python plugin
            # (infoPython() is a snapshot, not a reference: a copy taken
            # before the loop still reports sigmaGen = 0 afterwards), and it
            # is far too expensive -- and, repeated a million times,
            # unstable -- to re-fetch per event.  This configuration is
            # unweighted, so the weight column is exactly 1; check that
            # once rather than assume it, so a future biased run fails loudly
            # instead of silently writing 1.
            i = _info(p)
            if i.nWeights() != 1 or abs(i.weight() - 1.0) > 1e-12:
                sys.exit("weighted generation (nWeights = %d, weight = %g): "
                         "the per-event weight must be recorded per event, "
                         "which this writer does not do"
                         % (i.nWeights(), i.weight()))
        if not args.quiet and args.report_every > 0 \
                and n_done % args.report_every == 0:
            dt = time.time() - t0
            print("%d events, %.0f s (%.0f ev/s)"
                  % (n_done, dt, n_done / max(dt, 1e-9)), flush=True)

    a = buf.arrays()
    info = _info(p)                      # fresh: see the note in the loop
    sigma_gen_mb = info.sigmaGen()
    version = p.settings.parm("Pythia:versionNumber")
    meta = {"generator": "pythia8", "version": "%.3f" % version,
            "settings": settings, "target": args.target,
            "sigma_gen_mb": sigma_gen_mb,
            "sigma_err_mb": info.sigmaErr(), "n_events": n_done,
            "n_tried": n_tried, "seconds": time.time() - t0,
            "frame": "head-on: ion +z, electron -z; per-nucleon x,y,Q2"}
    np.savez_compressed(args.out, offsets=a["offsets"], pid=a["pid"],
                        charge=a["charge"], p4=a["p4"], x=a["x"], q2=a["q2"],
                        y=a["y"], kp=a["kp"], weight=a["weight"],
                        e_energy=np.array(args.electron_energy),
                        p_per_nucleon=np.array(args.p_per_nucleon),
                        meta=np.array(json.dumps(meta)))
    print("wrote %s: %d events, %d particles, sigma_gen = %.4g mb, %.0f s"
          % (args.out, n_done, a["pid"].size, sigma_gen_mb,
             time.time() - t0))
    if not args.quiet:
        p.stat()


if __name__ == "__main__":
    main()
