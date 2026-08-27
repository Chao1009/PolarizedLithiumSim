#!/usr/bin/env python3
"""Write a HepMC2 ASCII file of single ions for the far-forward scan.

npsim's particle gun cannot shoot a nucleus: DD4hep's
Geant4ParticleGenerator looks the name up in the G4 particle table and
"Li6", "ion(3,6)" and every other spelling come back "Bad particle type"
(checked 2026-08-26 in jug_xl-nightly).  A one-particle HepMC event does
work, and it also buys the thing a gun cannot do -- a scan of the
transverse kick in BOTH p_T and azimuth, which is what the Roman-Pot
cutout is sensitive to (plans/04 #20: the cutout geometry is still an
explicit assumption of the coherent channel).

Each event carries one ion of total momentum A * p_u along the ion axis
(the crossing angle, -25 mrad in x at IP6) plus a transverse kick
(p_T, phi) measured about that axis.  Event i of the output corresponds
to grid point i of the (p_T, phi) scan, in the order this script prints,
so the hit counting downstream (ff_gun_hits.py --per-event) can turn
event index straight back into (p_T, phi).

Usage (no container needed to write the file):
  python3 tools/fullsim/ion_gun_hepmc.py --out li6.hepmc \
      --pt 0 0.2 0.4 0.6 0.8 --nphi 8 --repeat 2
  singularity exec $SIF npsim --compactFile $DETECTOR_PATH/epic_craterlake_18x275.xml \
      --inputFiles li6.hepmc --numberOfEvents <n> ...
"""

import argparse

import numpy as np

# Ground-state nuclear masses [GeV] (fastsim/polli_fastsim/spectator.py
# NUCLEUS_MASS, AME2020 atomic masses less the electrons).
MASS = {1000030060: 5.601518702, 1000030070: 6.533833028,
        1000020040: 3.727379407, 1000010020: 1.875612942,
        1000010030: 2.808921133, 2212: 0.938272088}

XING_IP6 = -0.025          # ion beam direction in x at IP6 [rad]


def events(pdg, p_total, pts, phis, repeat, xing=XING_IP6):
    """(px, py, pz, E) of every scan point, ion axis = rotate z by xing."""
    m = MASS[pdg]
    axis = np.array([np.sin(xing), 0.0, np.cos(xing)])
    # transverse basis about the ion axis: x' in the bend plane, y' vertical
    ex = np.array([np.cos(xing), 0.0, -np.sin(xing)])
    ey = np.array([0.0, 1.0, 0.0])
    out = []
    for pt in pts:
        for phi in phis:
            pl = np.sqrt(max(p_total ** 2 - pt ** 2, 0.0))
            p = pl * axis + pt * (np.cos(phi) * ex + np.sin(phi) * ey)
            e = np.sqrt(p @ p + m * m)
            for _ in range(repeat):
                out.append((pt, phi, p[0], p[1], p[2], e, m))
    return out


def write_hepmc3(path, rows, pdg):
    """HepMC3 Asciiv3 -- one vertex at the origin, one status-1 ion.

    DD4hep routes a .hepmc/.hepmc3 input to its HEPMC3FileReader, which
    opens the file with HepMC3's ReaderAscii: a HepMC2 IO_GenEvent file
    hands it an immediate EOF ("Error when moving to event - EOF"), so
    the v3 record layout is the one to write.  Fields: E <event>
    <n_vertices> <n_particles>; P <id> <production vertex id> <pdg>
    px py pz e m <status>.  A vertex at the origin with nothing incoming
    is NOT written as a V line -- HepMC3's own WriterAscii drops it and
    marks the particle's production vertex 0, and a V line for it fails
    the reader's parse ("event parsing failed ... V -1 0 [] @ 0 0 0 0").
    The vertex COUNT in the E line is then 0, not the 1 that HepMC3's own
    writer emits: with 1 the reader still delivers the event but warns
    "not enough implicit vertices" on every one of them.  Both checked
    against pyHepMC3 inside the container.
    """
    with open(path, "w") as f:
        f.write("HepMC::Version 3.02.05\n")
        f.write("HepMC::Asciiv3-START_EVENT_LISTING\n")
        for i, (_pt, _phi, px, py, pz, e, m) in enumerate(rows):
            f.write("E %d 0 1\n" % i)
            f.write("U GEV MM\n")
            f.write("P 1 0 %d %.9g %.9g %.9g %.9g %.9g 1\n"
                    % (pdg, px, py, pz, e, m))
        f.write("HepMC::Asciiv3-END_EVENT_LISTING\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="ion_gun.hepmc")
    ap.add_argument("--pdg", type=int, default=1000030060,
                    help="1000030060 = 6Li, 1000020040 = alpha")
    ap.add_argument("--a", type=int, default=6)
    ap.add_argument("--p-per-nucleon", type=float, default=137.5)
    ap.add_argument("--pt", type=float, nargs="*",
                    default=(0.0, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0))
    ap.add_argument("--nphi", type=int, default=8)
    ap.add_argument("--repeat", type=int, default=1)
    args = ap.parse_args()

    phis = 2.0 * np.pi * np.arange(args.nphi) / max(args.nphi, 1)
    rows = events(args.pdg, args.a * args.p_per_nucleon, args.pt, phis,
                  args.repeat)
    write_hepmc3(args.out, rows, args.pdg)
    print("# %s: %d events, pdg %d, p = %g GeV (A = %d, %g GeV/u)"
          % (args.out, len(rows), args.pdg, args.a * args.p_per_nucleon,
             args.a, args.p_per_nucleon))
    print("# ievt pT[GeV] phi[deg]")
    for i, (pt, phi, *_r) in enumerate(rows):
        print("%d %g %.1f" % (i, pt, np.degrees(phi)))


if __name__ == "__main__":
    main()
