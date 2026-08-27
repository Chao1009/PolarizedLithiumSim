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
        print("Roman-Pot hit positions [mm] (mean over the hits of the event)")
        print("ievt " + ("  pT[GeV]  phi[deg] " if rows else "")
              + "     x        y        z      nhit")
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
                print(f"{i:4d} " + lbl
                      + f" {np.mean(xs[i]):8.1f} {np.mean(ys[i]):8.1f} "
                        f"{np.mean(zs[i]):8.1f} {len(xs[i]):6d}   {coll}")
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
