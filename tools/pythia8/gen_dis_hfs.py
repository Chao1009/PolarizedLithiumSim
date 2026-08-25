#!/usr/bin/env python3
"""Generate unpolarized e + N deep-inelastic events with PYTHIA 8 and write
the hadronic final state in the polligen HFS sample format
(evgen/polligen/hfs.py: HFSSample, compressed .npz).

Runs where the PYTHIA 8 Python bindings exist -- the eic-shell container
ships them (`import pythia8`).  It does NOT run on the Windows analysis
machine (no PYTHIA there); the resulting .npz is copied into the
repository (evgen/samples/) and consumed by evgen/scripts/hfs_resolution.py
and money_cos2phi_reco.py --y-source hfs --hfs-sample <file>.

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

Usage (inside eic-shell):
  python3 tools/pythia8/gen_dis_hfs.py --target p --n-events 1000000 \
      --electron-energy 10 --p-per-nucleon 50 --seed 1 \
      --out pythia8_e10_p50_dis.npz
  python3 tools/pythia8/gen_dis_hfs.py --target n ... --out pythia8_e10_n50_dis.npz
Merge p and n files with evgen/scripts/hfs_resolution.py --sample a.npz b.npz
(HFSSample.concatenate).  ~1 M events take of order an hour and ~200 MB.
"""

import argparse
import json
import sys
import time

import numpy as np

try:
    import pythia8
except ImportError:  # pragma: no cover
    sys.exit("pythia8 python module not found -- run inside eic-shell "
             "(or any environment with the PYTHIA 8 bindings)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="p", choices=("p", "n"))
    ap.add_argument("--n-events", type=int, default=100000)
    ap.add_argument("--electron-energy", type=float, default=10.0)
    ap.add_argument("--p-per-nucleon", type=float, default=50.0)
    ap.add_argument("--q2min", type=float, default=0.7)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="pythia8_dis_hfs.npz")
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
        "SpaceShower:dipoleRecoil = on",
        "PDF:lepton = off",
        "TimeShower:QEDshowerByL = off",
        "Random:setSeed = on",
        "Random:seed = %d" % args.seed,
        "Next:numberCount = 0",
        "Next:numberShowEvent = 0",
    ]
    for s in settings:
        p.readString(s)
    if not p.init():
        sys.exit("pythia init failed")

    offsets = [0]
    pid, charge, p4 = [], [], []
    xs, q2s, ys, kps, wts = [], [], [], [], []
    t0 = time.time()
    n_done = 0
    while n_done < args.n_events:
        if not p.next():
            continue
        ev = p.event
        # beams (records 1, 2), scattered electron = most energetic final e-
        k = np.array([ev[2].e(), ev[2].px(), ev[2].py(), ev[2].pz()])
        P = np.array([ev[1].e(), ev[1].px(), ev[1].py(), ev[1].pz()])
        best, ebest = -1, -1.0
        for i in range(ev.size()):
            pt = ev[i]
            if pt.isFinal() and pt.id() == 11 and pt.e() > ebest:
                best, ebest = i, pt.e()
        if best < 0:
            continue
        kp = np.array([ev[best].e(), ev[best].px(), ev[best].py(), ev[best].pz()])
        q = k - kp
        q2 = -(q[0] ** 2 - q[1] ** 2 - q[2] ** 2 - q[3] ** 2)
        pq = P[0] * q[0] - P[1] * q[1] - P[2] * q[2] - P[3] * q[3]
        pk = P[0] * k[0] - P[1] * k[1] - P[2] * k[2] - P[3] * k[3]
        y = pq / pk
        x = q2 / (2.0 * pq)
        for i in range(ev.size()):
            pt = ev[i]
            if not pt.isFinal() or i == best:
                continue
            pid.append(pt.id())
            charge.append(pt.charge())
            p4.append((pt.e(), pt.px(), pt.py(), pt.pz()))
        offsets.append(len(pid))
        xs.append(x); q2s.append(q2); ys.append(y); kps.append(kp)
        wts.append(p.info.weight())
        n_done += 1
        if not args.quiet and n_done % 10000 == 0:
            print("%d events, %.0f s" % (n_done, time.time() - t0), flush=True)

    sigma_gen_mb = p.info.sigmaGen()
    meta = {"generator": "pythia8", "version": pythia8.__version__
            if hasattr(pythia8, "__version__") else "unknown",
            "settings": settings, "target": args.target,
            "sigma_gen_mb": sigma_gen_mb, "n_events": n_done,
            "frame": "head-on: ion +z, electron -z; per-nucleon x,y,Q2"}
    np.savez_compressed(args.out, offsets=np.array(offsets, dtype=np.int64),
                        pid=np.array(pid, dtype=np.int64),
                        charge=np.array(charge, dtype=float),
                        p4=np.array(p4, dtype=float), x=np.array(xs),
                        q2=np.array(q2s), y=np.array(ys), kp=np.array(kps),
                        weight=np.array(wts),
                        e_energy=np.array(args.electron_energy),
                        p_per_nucleon=np.array(args.p_per_nucleon),
                        meta=np.array(json.dumps(meta)))
    print("wrote %s: %d events, sigma_gen = %.4g mb" % (args.out, n_done, sigma_gen_mb))
    if not args.quiet:
        p.stat()


if __name__ == "__main__":
    main()
