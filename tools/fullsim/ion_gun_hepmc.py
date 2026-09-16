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

Each event carries one ion of total momentum (1 + delta) * A * p_u along
the ion axis (the crossing angle, -25 mrad in x at IP6) plus a transverse
kick (p_T, phi) measured about that axis, from a vertex offset (x0, y0)
perpendicular to that axis.  Event i of the output corresponds to grid
point i of the scan, in the order this script prints, so the hit counting
downstream (ff_gun_hits.py --per-event) can turn event index straight
back into the scan point.

(x0, y0) and delta were added on 2026-09-15 for the transfer-matrix scan
of `ff_transfer_scan.py`: R12 = dx/dtheta_x and D = dx/dR need only the
(p_T, phi) grid, but R11 = dx/dx_IP and R21 = dtheta_x/dx_IP need the ion
moved off the reference orbit AT THE IP, and D-prime = dtheta_x/dR needs
the rigidity walked at fixed angle.  Both default to zero and the written
file is then byte-identical to what this script wrote before, V lines and
all (there are none).

Usage (no container needed to write the file):
  python3 tools/fullsim/ion_gun_hepmc.py --out li6.hepmc \
      --pt 0 0.2 0.4 0.6 0.8 --nphi 8 --repeat 2
  singularity exec $SIF npsim --compactFile $DETECTOR_PATH/epic_craterlake_18x275.xml \
      --inputFiles li6.hepmc --numberOfEvents <n> ...
"""

import argparse
import collections
import sys

import numpy as np

# Ground-state nuclear masses [GeV] (fastsim/polli_fastsim/spectator.py
# NUCLEUS_MASS, AME2020 atomic masses less the electrons).
MASS = {1000030060: 5.601518702, 1000030070: 6.533833028,
        1000020040: 3.727379407, 1000010020: 1.875612942,
        1000010030: 2.808921133, 2212: 0.938272088}

XING_IP6 = -0.025          # ion beam direction in x at IP6 [rad]

#: One scan point.  `pt`, `phi`, `x0`, `y0`, `delta` are the SCAN
#: coordinates (x0, y0 in mm, perpendicular to the ion axis; delta =
#: dp/p, i.e. R - 1 for a fixed species); `vx, vy, vz` is the vertex in
#: global mm and `px, py, pz, e, m` the four-momentum in GeV.
Point = collections.namedtuple(
    "Point", "pt phi x0 y0 delta vx vy vz px py pz e m")


def axis_basis(xing=XING_IP6):
    """(axis, ex, ey): the ion axis and the transverse basis about it.

    `ex` lies in the bend plane and `ey` is vertical, so a (p_T, phi)
    kick and an (x0, y0) vertex offset are both measured about the SAME
    axis -- which is what makes the fitted levers transfer-matrix
    elements of the ion line rather than of the global frame.
    """
    axis = np.array([np.sin(xing), 0.0, np.cos(xing)])
    ex = np.array([np.cos(xing), 0.0, -np.sin(xing)])
    ey = np.array([0.0, 1.0, 0.0])
    return axis, ex, ey


def point(pdg, p_total, pt, phi, x0=0.0, y0=0.0, delta=0.0, xing=XING_IP6):
    """One `Point` at momentum (1 + delta) * p_total, kick (pt, phi) and
    vertex (x0, y0) [mm] about the ion axis."""
    m = MASS[pdg]
    axis, ex, ey = axis_basis(xing)
    ptot = p_total * (1.0 + delta)
    pl = np.sqrt(max(ptot ** 2 - pt ** 2, 0.0))
    p = pl * axis + pt * (np.cos(phi) * ex + np.sin(phi) * ey)
    e = np.sqrt(p @ p + m * m)
    v = x0 * ex + y0 * ey
    return Point(pt, phi, x0, y0, delta, v[0], v[1], v[2],
                 p[0], p[1], p[2], e, m)


def events(pdg, p_total, pts, phis, repeat, xing=XING_IP6,
           offsets=((0.0, 0.0),), deltas=(0.0,)):
    """Every scan point of the (p_T, phi) x (x0, y0) x delta product.

    `offsets` and `deltas` default to the single zero entry, so the
    (p_T, phi) grid and its event ORDER are what they were before
    2026-09-15.
    """
    out = []
    for pt in pts:
        for phi in phis:
            for x0, y0 in offsets:
                for delta in deltas:
                    for _ in range(repeat):
                        out.append(point(pdg, p_total, pt, phi,
                                         x0, y0, delta, xing))
    return out


def write_hepmc3(path, rows, pdg):
    """HepMC3 Asciiv3 -- one status-1 ion, at the origin or at (x0, y0).

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

    AN OFF-ORIGIN VERTEX therefore cannot be written as a bare `V` line
    either, and the transfer-matrix scan needs one.  What works, measured
    2026-09-15 against npsim itself (`PrimaryHandler INFO +++++
    G4PrimaryVertex at (...)` echoes the vertex of every event it reads),
    is to give the vertex something INCOMING: a status-4 beam copy of the
    same ion with production vertex 0, then `V -1 0 [1] @ x y z t`, then
    the status-1 ion with production vertex -1.  Geant4 tracks status 1
    only, so the beam copy costs nothing.  Events whose vertex is exactly
    the origin keep the old two-line record, which is why a file written
    without --x0/--y0 is byte-identical to what this script wrote before.
    """
    with open(path, "w") as f:
        f.write("HepMC::Version 3.02.05\n")
        f.write("HepMC::Asciiv3-START_EVENT_LISTING\n")
        for i, r in enumerate(rows):
            px, py, pz, e, m = r.px, r.py, r.pz, r.e, r.m
            if r.vx == 0.0 and r.vy == 0.0 and r.vz == 0.0:
                f.write("E %d 0 1\n" % i)
                f.write("U GEV MM\n")
                f.write("P 1 0 %d %.9g %.9g %.9g %.9g %.9g 1\n"
                        % (pdg, px, py, pz, e, m))
                continue
            f.write("E %d 1 2\n" % i)
            f.write("U GEV MM\n")
            f.write("P 1 0 %d %.9g %.9g %.9g %.9g %.9g 4\n"
                    % (pdg, px, py, pz, e, m))
            f.write("V -1 0 [1] @ %.9g %.9g %.9g 0\n" % (r.vx, r.vy, r.vz))
            f.write("P 2 -1 %d %.9g %.9g %.9g %.9g %.9g 1\n"
                    % (pdg, px, py, pz, e, m))
        f.write("HepMC::Asciiv3-END_EVENT_LISTING\n")


#: Header of the index file, and what `ff_transfer_scan.py` parses.  The
#: three-column form is what every scan before 2026-09-15 wrote and what
#: `ff_gun_hits.py --index` reads (it uses columns 1 and 2 only), so the
#: extra columns appear only when the scan actually moves the vertex or
#: the rigidity.
INDEX_HEADER_3 = "# ievt pT[GeV] phi[deg]"
INDEX_HEADER_6 = "# ievt pT[GeV] phi[deg] x0[mm] y0[mm] delta"


def write_index(rows, out):
    """The event index, three columns or six (see INDEX_HEADER_*)."""
    wide = any(r.x0 or r.y0 or r.delta for r in rows)
    out.write((INDEX_HEADER_6 if wide else INDEX_HEADER_3) + "\n")
    for i, r in enumerate(rows):
        line = "%d %g %.1f" % (i, r.pt, np.degrees(r.phi))
        if wide:
            line += " %g %g %g" % (r.x0, r.y0, r.delta)
        out.write(line + "\n")


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
    ap.add_argument("--x0", type=float, nargs="*", default=(0.0,),
                    help="vertex offsets [mm] in the bend plane, "
                         "perpendicular to the ion axis")
    ap.add_argument("--y0", type=float, nargs="*", default=(0.0,),
                    help="vertical vertex offsets [mm]")
    ap.add_argument("--delta", type=float, nargs="*", default=(0.0,),
                    help="dp/p about --p-per-nucleon (= R - 1 at fixed A, Z)")
    args = ap.parse_args()

    phis = 2.0 * np.pi * np.arange(args.nphi) / max(args.nphi, 1)
    offsets = tuple((x0, y0) for x0 in args.x0 for y0 in args.y0)
    rows = events(args.pdg, args.a * args.p_per_nucleon, args.pt, phis,
                  args.repeat, offsets=offsets, deltas=tuple(args.delta))
    write_hepmc3(args.out, rows, args.pdg)
    print("# %s: %d events, pdg %d, p = %g GeV (A = %d, %g GeV/u)"
          % (args.out, len(rows), args.pdg, args.a * args.p_per_nucleon,
             args.a, args.p_per_nucleon))
    write_index(rows, sys.stdout)


if __name__ == "__main__":
    main()
