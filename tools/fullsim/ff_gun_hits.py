#!/usr/bin/env python3
"""Count far-forward hits per subsystem for the gun-scan outputs.
Run INSIDE eic-shell (needs uproot):
  singularity exec $SIF python3 ff_gun_hits.py <outdir-from-ff_gun_scan>

With --per-event the same counting is reported event by event instead of
averaged over the file, which is what the ion scan of ion_gun_hepmc.py
needs: there every event is a different (p_T, phi) of the same nucleus,
so the per-event table IS the acceptance map.  Pass the index file that
ion_gun_hepmc.py prints as --index to label the rows, and --positions to
print where in the Roman-Pot plane each accepted ion landed (the cutout
geometry itself, plans/04 #20).
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import uproot

PATTERNS = {
    "RomanPots": re.compile(r"ForwardRomanPot.*Hits$"),
    "OMD": re.compile(r"ForwardOffM.*Hits$"),
    "B0": re.compile(r"B0Tracker.*Hits$|B0ECal.*Hits$"),
    "ZDC": re.compile(r"(Hcal|Ecal)FarForwardZDC.*Hits$|ZDC.*Hits$"),
}

#: Roman-Pot sensor planes of the current ePIC geometry, read from
#: compact/far_forward/roman_pots_eRD24_design.xml (epic-main
#: 9aaa296976d3ad9de404f775ae89fc17a068c07c): two stations of two layers,
#: 20 mm apart along the station normal.  Both stations write into the
#: SAME readout (ForwardRomanPotHits), so an event-level mean over the
#: collection mixes them and is not a transport lever -- group by plane
#: instead (plans/09 B1).
#:
#: THE STATIONS ARE ROTATED about y by -0.04545 rad (the crossing angle),
#: so a plane is NOT a surface of constant global z: |dz| = 0.0454 |dx|
#: along it, which is 5 mm at |dx| = 110 mm and 6.5 mm at the 144 mm
#: module edge.  Grouping on global z with the +-15 mm window this file
#: used until 2026-08-28 therefore put every hit beyond |dx| ~ 110 mm into
#: BOTH layers of its station (69 of the 358 Roman-pot events of
#: lad_18x275 had per-plane counts summing above the event hit count), and
#: the inflated multiplicity then tripped the analysis' cleanliness cut
#: and dropped the outermost rows.  Hits are assigned in the ROTATED frame
#: instead: with (x0, z0) the plane centre and th the station rotation,
#:
#:     u = (x - x0) cos th - (z - z0) sin th      (in-plane, "dx")
#:     w = (x - x0) sin th + (z - z0) cos th      (along the normal)
#:
#: and a hit belongs to the plane when |w| <= PLANE_TOL_MM.  Layer 2 sits
#: at local z = +20 mm, i.e. global (+20 sin th, +20 cos th) = (-0.91,
#: +19.98) mm from layer 1, so each plane carries its own x0.
STATION_ROTATION_RAD = -0.04545
_S1 = (-1131.19, 32547.3)
_S2 = (-1208.43, 34245.5)
_DX_L2 = 20.0 * np.sin(STATION_ROTATION_RAD)
_DZ_L2 = 20.0 * np.cos(STATION_ROTATION_RAD)
#: (name, plane-centre x [mm], plane-centre z [mm])
STATION_PLANES = (("S1L1", _S1[0], _S1[1]),
                  ("S1L2", _S1[0] + _DX_L2, _S1[1] + _DZ_L2),
                  ("S2L1", _S2[0], _S2[1]),
                  ("S2L2", _S2[0] + _DX_L2, _S2[1] + _DZ_L2))
PLANE_TOL_MM = 5.0


def plane_coords(x, z, x0, z0, th=STATION_ROTATION_RAD):
    """(u, w): the in-plane offset from a plane's centre and the distance
    along its normal, for hits at global (x, z).  See STATION_PLANES."""
    c, s = np.cos(th), np.sin(th)
    dx, dz = np.asarray(x, float) - x0, np.asarray(z, float) - z0
    return c * dx - s * dz, s * dx + c * dz


def hit_branches(tree, pattern):
    """Every <collection>.cellID branch whose collection matches."""
    out = []
    for b in tree.keys():
        base = b.split("/")[0].split(".")[0]
        if pattern.match(base) and b.endswith(".cellID"):
            out.append(b)
    return out


def per_event(path, index=None, positions=False):
    """One row per event: which far-forward systems saw >= 1 hit."""
    t = uproot.open(path)["events"]
    nev = t.num_entries
    hit = {}
    for label, pat in PATTERNS.items():
        seen = np.zeros(nev, dtype=bool)
        for b in hit_branches(t, pat):
            counts = t[b].array(library="np")
            seen |= np.array([len(c) > 0 for c in counts])
        hit[label] = seen
    rows = None
    if index:
        rows = [ln.split() for ln in open(index)
                if ln.strip() and not ln.startswith("#")]
    head = "ievt " + ("  pT[GeV]  phi[deg] " if rows else "")
    print(head + " ".join(f"{k:>10s}" for k in PATTERNS))
    for i in range(nev):
        lbl = ""
        if rows and i < len(rows):
            lbl = f"  {float(rows[i][1]):7.2f}  {float(rows[i][2]):7.1f} "
        print(f"{i:4d} " + lbl
              + " ".join(f"{int(hit[k][i]):10d}" for k in PATTERNS))
    if positions:
        print()
        print("Roman-Pot hit positions [mm], GROUPED BY STATION PLANE.")
        print("Planes are assigned in the ROTATED station frame (rotation "
              "%.5f rad about y): a hit belongs to a plane when its "
              "distance ALONG THE NORMAL is |w| <= %.0f mm.  Grouping on "
              "global z instead double-counts every hit beyond |dx| ~ 110 "
              "mm into both layers of its station."
              % (STATION_ROTATION_RAD, PLANE_TOL_MM))
        print("Plane centres (x, z) [mm]: " + ", ".join(
            "%s (%.2f, %.2f)" % p for p in STATION_PLANES))
        print("'dx' columns are the IN-PLANE offset u from the plane's own "
              "centre -- the quantity R12 = du/dtheta_x is regressed on. "
              "'ALL' is the old, station-mixing row kept for continuity, "
              "and its dx is taken against station 1.")
        print("The silicon y is tiled by the 16 mm module, so mean y is not "
              "a lever -- read miny/maxy.")
        print("ievt " + ("  pT[GeV]  phi[deg] " if rows else "")
              + " plane  nhit    meandx    meany     mindx    maxdx"
                "      miny     maxy       z")
        for b in hit_branches(t, PATTERNS["RomanPots"]):
            coll = b.split(".")[0]
            xs = t[f"{coll}.position.x"].array(library="np")
            ys = t[f"{coll}.position.y"].array(library="np")
            zs = t[f"{coll}.position.z"].array(library="np")
            for i in range(nev):
                if len(xs[i]) == 0:
                    continue
                lbl = ""
                if rows and i < len(rows):
                    lbl = (f"  {float(rows[i][1]):7.2f}  "
                           f"{float(rows[i][2]):7.1f} ")
                x, y, z = (np.asarray(xs[i], dtype=float),
                           np.asarray(ys[i], dtype=float),
                           np.asarray(zs[i], dtype=float))
                groups = []
                us = {}
                for name, x0, z0 in STATION_PLANES:
                    u, w = plane_coords(x, z, x0, z0)
                    us[name] = u
                    groups.append((name, np.abs(w) <= PLANE_TOL_MM))
                assigned = np.zeros(len(z), dtype=bool)
                for _n, m in groups:
                    assigned |= m
                if not assigned.all():
                    groups.append(("other", ~assigned))
                u_all = plane_coords(x, z, *STATION_PLANES[0][1:])[0]
                for nm in ("other", "ALL"):
                    us[nm] = u_all
                groups.append(("ALL", np.ones(len(z), dtype=bool)))
                for name, m in groups:
                    if not m.any():
                        continue
                    u = us[name]
                    print(f"{i:4d} " + lbl
                          + f" {name:5s} {int(m.sum()):5d} "
                            f"{u[m].mean():9.2f} {y[m].mean():8.2f} "
                            f"{u[m].min():9.2f} {u[m].max():8.2f} "
                            f"{y[m].min():9.2f} {y[m].max():8.2f} "
                            f"{z[m].mean():9.1f}   {coll}")
    return hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", default="/tmp/ff_gun",
                    help="directory of gun_*.edm4hep.root, or one file")
    ap.add_argument("--per-event", action="store_true")
    ap.add_argument("--index", default=None,
                    help="index file printed by ion_gun_hepmc.py")
    ap.add_argument("--positions", action="store_true")
    args = ap.parse_args()
    outdir = args.target
    if args.per_event:
        files = ([outdir] if os.path.isfile(outdir)
                 else sorted(glob.glob(f"{outdir}/gun_*.edm4hep.root")))
        if not files:
            sys.exit(f"no gun_*.edm4hep.root in {outdir}")
        for f in files:
            print(f"# {f}")
            per_event(f, args.index, args.positions)
        return
    files = ([outdir] if os.path.isfile(outdir)
             else sorted(glob.glob(f"{outdir}/gun_*.edm4hep.root")))
    if not files:
        sys.exit(f"no gun_*.edm4hep.root in {outdir}")
    print(f"{'config':16s} {'nev':>4s} " +
          " ".join(f"{k:>10s}" for k in PATTERNS) + "   (fraction of events with >=1 hit)")
    for f in files:
        tag = f.split("gun_")[-1].replace(".edm4hep.root", "")
        try:
            t = uproot.open(f)["events"]
        except Exception as e:
            print(f"{tag:16s} OPEN-FAIL {e}")
            continue
        nev = t.num_entries
        cols = {}
        for label, pat in PATTERNS.items():
            frac = 0.0
            for b in t.keys():
                base = b.split("/")[0].split(".")[0]
                if pat.match(base) and b.endswith(".cellID"):
                    counts = t[b].array(library="np")
                    frac = max(frac,
                               float(sum(1 for c in counts if len(c) > 0))
                               / max(nev, 1))
            cols[label] = frac
        print(f"{tag:16s} {nev:4d} " +
              " ".join(f"{cols[k]:10.2f}" for k in PATTERNS))
    # list available FF-ish collections once, for reference
    t = uproot.open(files[0])["events"]
    ff = sorted({b.split(".")[0] for b in t.keys()
                 if re.search(r"RomanPot|OffM|B0|ZDC", b)})
    print("\ncollections present:", ", ".join(ff[:12]))


if __name__ == "__main__":
    main()
