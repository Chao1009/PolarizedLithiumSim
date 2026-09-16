#!/usr/bin/env python3
"""The far-forward TRANSFER MATRIX, measured (plans/03 step 2.2 (2)).

The 2026-08-28/29 ion scans measured three elements of the ion line
between IP6 and the Roman Pots -- R12 = dx/dtheta_x, R34 = dy/dtheta_y
and D = dx/dR (with its second-order term D2).  All three are `dPOT/dIP`
columns of the FIRST row of the transfer matrix: the ladders walk the
angle at the IP and the dispersion walks the rigidity, and every one of
them is read off a POSITION at the pot plane.

The other half of the matrix was never measured.  Writing the ion line
between the IP and one pot plane in the usual (x, x') form,

    x_pot    =  R11 x_IP + R12 x'_IP + D  delta
    x'_pot   =  R21 x_IP + R22 x'_IP + D' delta

what was missing was everything that needs either a DISPLACEMENT at the
IP (R11, R21) or an ANGLE at the pot (R21, R22, D').  This script
measures those four, and the vertical R33/R43/R44 alongside them, in the
same geometry and with the same regression the three published levers
use, so that the two halves of the matrix are one measurement.

Two things make it possible that were not available in August.

*The IP displacement.*  `ion_gun_hepmc.py` now writes a vertex, so
(x0, y0) at the IP is a scan coordinate like (p_T, phi).  See its
`write_hepmc3` for the HepMC3 record that DD4hep's reader accepts.

*The angle at the pot.*  `ForwardRomanPotHits` is an edm4hep
SimTrackerHit collection and carries `momentum` as well as `position`:
the outgoing angle is READ, not differenced between stations, so R21,
R22 and D' come from station 1 alone -- the same plane R12, R34 and D
are quoted at -- and the 1.7 m station lever is left as a cross-check
rather than being the measurement.  The same `momentum` branch also
identifies the primary: a 6Li at the ring reference rigidity carries
123 / 300 / 825 GeV and every secondary in the plane carries far less,
so `|p| > PRIMARY_P_FRACTION * p_beam` isolates it without the hit-count
cleanliness cut the position-only fits needed.

*The pots must not be retracted.*  A transfer matrix is a property of
the magnets, and the per-energy insertion holds the silicon off the beam
axis -- 16 / 32 / 48 mm horizontally at y ~ 0 -- so an ion at theta = 0
with a millimetre of IP offset lands in the blind band and is never
seen.  The scan therefore runs through the ZERO-INSERTION geometry, the
recipe `tools/fullsim/README.md` already documents and validated for the
5 x 41 vertical lever: the four `offset_*_RP_section` constants set to
0.0 cm, every field untouched.  R12, R34 and D come out of the same file
as controls, and the published numbers are what they must reproduce.

Usage (build needs no container; the npsim leg does):

    python3 tools/fullsim/ff_transfer_scan.py build --config 18x275 --out $S
    # prints the npsim command to run inside eic-shell, then:
    python3 tools/fullsim/ff_transfer_scan.py fit $S --config 18x275

`fit` runs on the host: uproot reads the edm4hep file directly.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ion_gun_hepmc as gun                                      # noqa: E402
from ff_gun_hits import (STATION_PLANES, PLANE_TOL_MM,           # noqa: E402
                         STATION_ROTATION_RAD, plane_coords)

#: 6Li total momentum [GeV] that puts the ion on each ring's REFERENCE
#: rigidity -- the rigidity-scaled matching momenta of the aperture scan
#: (20.5 / 50 / 137.5 GeV/u x A = 6), not the gamma-matched fill menu.
#: The transfer matrix is a property of the lattice, so the reference
#: rigidity is the momentum to measure it at.
REFERENCE_MOMENTUM = {"5x41": 123.0, "10x100": 300.0, "18x275": 825.0}

#: A hit is the primary when it carries this fraction of the beam
#: momentum.  0.5 separates an 825 GeV 6Li from anything it knocks out.
PRIMARY_P_FRACTION = 0.5

#: Scan legs.  `thx`/`thy` reproduce the published ladders (and are the
#: controls for R12 and R34); `x0`/`y0` and `dlt` are the new ones;
#: `b0x`/`b0y` walk out to the B0 window, 5.5-20 mrad, which no ladder
#: in this repository had ever reached -- that is the B0 collection
#: joining the scan rather than reporting zeros.
LADDER_MRAD = tuple(0.2 * (i + 1) for i in range(20))       # 0.2 .. 4.0
B0_MRAD = (6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0)
OFFSETS_MM = (0.0, 0.5, -0.5, 1.0, -1.0, 2.0, -2.0, 4.0, -4.0, 8.0, -8.0)
#: The top-up offsets of the second pass.  R21 = dtheta_x/dx_IP is the
#: weakest of the four -- the angle a millimetre at the IP buys is tens
#: of microradians -- so its lever arm is doubled in a second, 34-event
#: file rather than by re-running the whole scan; `fit` merges every
#: scan*.index / scan*.edm4hep.root pair in the directory.
OFFSETS_MM_TOPUP = (12.0, -12.0, 16.0, -16.0, 24.0, -24.0, 32.0, -32.0)
DELTAS = (0.0, 0.02, -0.02, 0.05, -0.05, 0.10, -0.10, 0.15, -0.15)

INDEX_HEADER = "# ievt pT[GeV] phi[deg] x0[mm] y0[mm] delta leg"


def scan_points(config, pdg=1000030060, legs=None, topup=False):
    """Every (Point, leg) of the transfer scan at `config`."""
    p0 = REFERENCE_MOMENTUM[config]
    out = []
    offs = OFFSETS_MM_TOPUP if topup else OFFSETS_MM

    def add(leg, pt, phi, x0=0.0, y0=0.0, delta=0.0):
        if legs and leg not in legs:
            return
        out.append((gun.point(pdg, p0, pt, phi, x0, y0, delta), leg))

    for th in LADDER_MRAD:
        for phi in (0.0, np.pi):
            add("thx", 1e-3 * th * p0, phi)
    for th in LADDER_MRAD:
        for phi in (0.5 * np.pi, 1.5 * np.pi):
            add("thy", 1e-3 * th * p0, phi)
    for x0 in offs:
        add("x0", 0.0, 0.0, x0=x0)
    for y0 in (offs if topup else offs[1:]):
        add("y0", 0.0, 0.0, y0=y0)
    for d in DELTAS:
        add("dlt", 0.0, 0.0, delta=d)
    for th in B0_MRAD:
        for phi in (0.0, np.pi):
            add("b0x", 1e-3 * th * p0, phi)
    for th in B0_MRAD:
        for phi in (0.5 * np.pi, 1.5 * np.pi):
            add("b0y", 1e-3 * th * p0, phi)
    return out


def write_scan(outdir, config, pdg=1000030060, stem="scan", legs=None,
               topup=False):
    """`<stem>.hepmc` + `<stem>.index`; returns the event count."""
    os.makedirs(outdir, exist_ok=True)
    pts = scan_points(config, pdg, legs, topup)
    rows = [p for p, _leg in pts]
    gun.write_hepmc3(os.path.join(outdir, stem + ".hepmc"), rows, pdg)
    with open(os.path.join(outdir, stem + ".index"), "w") as f:
        f.write(INDEX_HEADER + "\n")
        for i, (r, leg) in enumerate(pts):
            f.write("%d %g %.1f %g %g %g %s\n"
                    % (i, r.pt, np.degrees(r.phi), r.x0, r.y0, r.delta, leg))
    return len(pts)


def read_index(path):
    """The rows of `scan.index` as dicts, in event order."""
    rows = []
    for ln in open(path):
        if not ln.strip() or ln.startswith("#"):
            continue
        f = ln.split()
        rows.append(dict(ievt=int(f[0]), pt=float(f[1]),
                         phi=np.radians(float(f[2])), x0=float(f[3]),
                         y0=float(f[4]), delta=float(f[5]),
                         leg=f[6] if len(f) > 6 else "?"))
    return rows


def scan_files(outdir, root=None):
    """Every (index, edm4hep) pair of a scan directory, `scan*` first.

    A scan may be spread over several npsim jobs -- the top-up of a weak
    leg, or a repeat -- and the fits are pooled over all of them.
    """
    import glob
    if root:
        return [(os.path.join(outdir, "scan.index"), root)]
    pairs = []
    for idx in sorted(glob.glob(os.path.join(outdir, "scan*.index"))):
        rt = idx[:-len(".index")] + ".edm4hep.root"
        if os.path.exists(rt):
            pairs.append((idx, rt))
    if not pairs:
        raise SystemExit("no scan*.index / scan*.edm4hep.root pair in "
                         + outdir)
    return pairs


# ---------------------------------------------------------------- fits

def robust_fit(x, y):
    """Theil-Sen slope, then least squares on the inliers.

    The same two-step the published R12/R34 fits use: the median of
    pairwise slopes is immune to the shower rows, and the inlier band is
    4 x 1.4826 x MAD of the Theil-Sen residual rather than the fixed
    2 mm of the position-only fit, so it carries over unchanged to an
    ANGLE regression whose natural scale is microradians.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    if n < 3:
        return None
    sl = [(y[j] - y[i]) / (x[j] - x[i])
          for i in range(n) for j in range(i + 1, n) if x[j] != x[i]]
    if not sl:
        return None
    s0 = float(np.median(sl))
    b0 = float(np.median(y - s0 * x))
    r = y - (s0 * x + b0)
    mad = float(np.median(np.abs(r - np.median(r))))
    tol = max(4.0 * 1.4826 * mad, 1e-15)
    keep = np.abs(r) <= tol
    if keep.sum() >= 3:
        s, b = np.polyfit(x[keep], y[keep], 1)
    else:
        keep, s, b = np.ones(n, bool), s0, b0
    res = y[keep] - (s * x[keep] + b)
    m = int(keep.sum())
    rms = float(np.sqrt(np.mean(res ** 2)))
    spread = float(np.std(y[keep]))
    # standard error of the slope: rms x sqrt(m / (m - 2)) over
    # sqrt(m) x the rms lever arm.  THIS is the "fit residual" the
    # acceptance is quoted on -- how well the LEVER is determined, not
    # how much of the dependent variable's spread the line explains.
    lever = float(np.std(x[keep]))
    sig = (rms * np.sqrt(m / max(m - 2, 1)) / (np.sqrt(m) * lever)
           if lever > 0 and m > 2 else float("nan"))
    return dict(slope=float(s), intercept=float(b), n=m, ntot=n, rms=rms,
                sigma=float(sig),
                rel=(abs(sig / s) if s != 0 else float("nan")),
                frac=(rms / spread if spread > 0 else float("nan")),
                x=x[keep], y=y[keep], res=res)


def _fmt(tag, unit, fit, scale=1.0):
    if fit is None:
        return "%-6s %-9s   --  (no rows)" % (tag, unit)
    return ("%-6s %-9s %10.4f +- %-8.4f n %3d/%-3d  resid %5.2f%%  "
            "rms %.4g  (unexplained %4.1f%%)"
            % (tag, unit, fit["slope"] * scale, fit["sigma"] * scale,
               fit["n"], fit["ntot"], 100.0 * fit["rel"], fit["rms"],
               100.0 * fit["frac"]))


# ------------------------------------------------------------- reading

def primary_rp(path, p_beam, planes=STATION_PLANES):
    """Per event and per Roman-Pot plane, the primary's (u, y, thx, thy).

    `u` is the in-plane offset from the plane centre in the ROTATED
    station frame, and the angles are the hit MOMENTUM rotated into the
    same frame, so position and angle share one convention.
    """
    import uproot
    t = uproot.open(path)["events"]
    c = "ForwardRomanPotHits"
    px = t[f"{c}.momentum.x"].array(library="np")
    py = t[f"{c}.momentum.y"].array(library="np")
    pz = t[f"{c}.momentum.z"].array(library="np")
    hx = t[f"{c}.position.x"].array(library="np")
    hy = t[f"{c}.position.y"].array(library="np")
    hz = t[f"{c}.position.z"].array(library="np")
    th = STATION_ROTATION_RAD
    cth, sth = np.cos(th), np.sin(th)
    out = {}
    for name, x0, z0 in planes:
        rows = {}
        for i in range(t.num_entries):
            if len(hx[i]) == 0:
                continue
            x = np.asarray(hx[i], float)
            y = np.asarray(hy[i], float)
            z = np.asarray(hz[i], float)
            qx = np.asarray(px[i], float)
            qy = np.asarray(py[i], float)
            qz = np.asarray(pz[i], float)
            u, w = plane_coords(x, z, x0, z0)
            p = np.sqrt(qx ** 2 + qy ** 2 + qz ** 2)
            m = (np.abs(w) <= PLANE_TOL_MM) & (p > PRIMARY_P_FRACTION * p_beam)
            if not m.any():
                continue
            pu = cth * qx[m] - sth * qz[m]
            pw = sth * qx[m] + cth * qz[m]
            rows[i] = dict(u=float(u[m].mean()), y=float(y[m].mean()),
                           thx=float(np.mean(pu / pw)),
                           thy=float(np.mean(qy[m] / pw)),
                           n=int(m.sum()))
        out[name] = rows
    return out


def b0_zdc_summary(path, p_beam):
    """The B0 and ZDC rows of the scan, by leg.

    B0: the primary's mean (x, y) and the fitted dx/dtheta_x, dy/dtheta_y
    in the first B0 tracker layer.  ZDC: the fraction of events with a
    hit and the energy-weighted centroid -- the ZDC is a calorimeter with
    no momentum branch, so what it contributes to a transfer measurement
    is where the breakup lands, not a lever.
    """
    import uproot
    t = uproot.open(path)["events"]
    keys = set(k.split("/")[0] for k in t.keys())
    out = {}
    if "B0TrackerHits" in keys:
        c = "B0TrackerHits"
        hx = t[f"{c}.position.x"].array(library="np")
        hy = t[f"{c}.position.y"].array(library="np")
        hz = t[f"{c}.position.z"].array(library="np")
        qx = t[f"{c}.momentum.x"].array(library="np")
        qy = t[f"{c}.momentum.y"].array(library="np")
        qz = t[f"{c}.momentum.z"].array(library="np")
        zs = np.concatenate([np.asarray(z, float) for z in hz if len(z)]) \
            if any(len(z) for z in hz) else np.array([])
        layers = []
        for z in np.sort(np.unique(np.round(zs, 0))):
            if not layers or z - layers[-1] > 20.0:
                layers.append(float(z))
        rows = {}
        for i in range(t.num_entries):
            if len(hx[i]) == 0:
                continue
            p = np.sqrt(np.asarray(qx[i], float) ** 2
                        + np.asarray(qy[i], float) ** 2
                        + np.asarray(qz[i], float) ** 2)
            m = p > PRIMARY_P_FRACTION * p_beam
            if not m.any() or not layers:
                continue
            z = np.asarray(hz[i], float)[m]
            first = np.abs(z - layers[0]) < 20.0
            if not first.any():
                continue
            rows[i] = dict(x=float(np.asarray(hx[i], float)[m][first].mean()),
                           y=float(np.asarray(hy[i], float)[m][first].mean()),
                           n=int(first.sum()))
        out["B0"] = dict(layers=layers, rows=rows)
    zdc = {}
    for c in ("HcalFarForwardZDCHits", "EcalFarForwardZDCHits"):
        if c not in keys:
            continue
        e = t[f"{c}.energy"].array(library="np")
        zx = t[f"{c}.position.x"].array(library="np")
        zy = t[f"{c}.position.y"].array(library="np")
        rows = {}
        for i in range(t.num_entries):
            w = np.asarray(e[i], float)
            if w.sum() <= 0:
                continue
            rows[i] = dict(e=float(w.sum()),
                           x=float(np.average(np.asarray(zx[i], float),
                                              weights=w)),
                           y=float(np.average(np.asarray(zy[i], float),
                                              weights=w)))
        zdc[c] = rows
    out["ZDC"] = zdc
    return out


# -------------------------------------------------------------- report

def fit_config(outdir, config, root=None, plane="S1L1"):
    pairs = scan_files(outdir, root)
    p_beam = REFERENCE_MOMENTUM[config]
    idx, files, nev = [], [], 0
    for ipath, rpath in pairs:
        rows = primary_rp(rpath, p_beam)[plane]
        rix = read_index(ipath)
        idx.append((rix, rows))
        files.append(os.path.basename(rpath))
        nev += len(rix)

    def leg(name):
        return [(r, rows[r["ievt"]]) for rix, rows in idx for r in rix
                if r["leg"] == name and r["ievt"] in rows]

    def th_ip(r):
        """(theta_x, theta_y) at the IP, about the ion axis."""
        p = p_beam * (1.0 + r["delta"])
        th = np.arcsin(min(r["pt"] / p, 1.0))
        return th * np.cos(r["phi"]), th * np.sin(r["phi"])

    print("=" * 74)
    print("far-forward transfer matrix, %s, plane %s, %d/%d events with a "
          "primary   [%s]"
          % (config, plane, sum(len(r) for _x, r in idx), nev,
             ", ".join(files)))
    print("  (slope +- its standard error, inliers/rows, RESID = the "
          "error as a fraction of the slope --")
    print("   the acceptance metric -- then the residual rms in the "
          "dependent variable's own unit)")
    print("=" * 74)

    out = {}
    lx = leg("thx")
    out["R12"] = robust_fit([th_ip(r)[0] for r, _h in lx],
                            [h["u"] for _r, h in lx])
    out["R22"] = robust_fit([th_ip(r)[0] for r, _h in lx],
                            [h["thx"] for _r, h in lx])
    ly = leg("thy")
    out["R34"] = robust_fit([th_ip(r)[1] for r, _h in ly],
                            [h["y"] for _r, h in ly])
    out["R44"] = robust_fit([th_ip(r)[1] for r, _h in ly],
                            [h["thy"] for _r, h in ly])
    lx0 = leg("x0")
    out["R11"] = robust_fit([r["x0"] for r, _h in lx0],
                            [h["u"] for _r, h in lx0])
    out["R21"] = robust_fit([r["x0"] for r, _h in lx0],
                            [h["thx"] for _r, h in lx0])
    ly0 = leg("y0")
    out["R33"] = robust_fit([r["y0"] for r, _h in ly0],
                            [h["y"] for _r, h in ly0])
    out["R43"] = robust_fit([r["y0"] for r, _h in ly0],
                            [h["thy"] for _r, h in ly0])
    ld = leg("dlt")
    out["D"] = robust_fit([r["delta"] for r, _h in ld],
                          [h["u"] for _r, h in ld])
    out["Dp"] = robust_fit([r["delta"] for r, _h in ld],
                           [h["thx"] for _r, h in ld])

    print("-- CONTROLS (must reproduce POT_LEVERS / POT_DISPERSION)")
    print(_fmt("R12", "[m]", out["R12"], 1e-3))
    print(_fmt("R34", "[m]", out["R34"], 1e-3))
    print(_fmt("D", "[m]", out["D"], 1e-3))
    # the dispersion leg spans |delta| <= 0.15, where the SECOND-order
    # term is 5 mm at 18 x 275 -- most of the linear fit's residual.
    # Fitting the quadratic here is the control on POT_DISPERSION_2.
    if len(ld) >= 4:
        c2, c1, _c0 = np.polyfit([r["delta"] for r, _h in ld],
                                 [h["u"] for _r, h in ld], 2)
        out["D_quad"], out["D2"] = c1 * 1e-3, c2 * 1e-3
        print("D quad [m]     %10.4f   with D2 %8.4f [m]   "
              "(control on POT_DISPERSION_2)" % (c1 * 1e-3, c2 * 1e-3))
    print("-- NEW: the second row of the matrix, and R11")
    print(_fmt("R11", "[1]", out["R11"]))
    print(_fmt("R21", "[rad/m]", out["R21"], 1e3))
    print(_fmt("R22", "[1]", out["R22"]))
    print(_fmt("D'", "[rad]", out["Dp"]))
    print("-- NEW: the vertical block")
    print(_fmt("R33", "[1]", out["R33"]))
    print(_fmt("R43", "[rad/m]", out["R43"], 1e3))
    print(_fmt("R44", "[1]", out["R44"]))
    if all(out[k] for k in ("R11", "R22", "R12", "R21")):
        det = (out["R11"]["slope"] * out["R22"]["slope"]
               - (out["R12"]["slope"] * 1e-3) * (out["R21"]["slope"] * 1e3))
        print("-- symplectic check   R11 R22 - R12 R21 = %.4f   (1 exactly "
              "for a linear line at fixed rigidity)" % det)
    if all(out[k] for k in ("R33", "R44", "R34", "R43")):
        det = (out["R33"]["slope"] * out["R44"]["slope"]
               - (out["R34"]["slope"] * 1e-3) * (out["R43"]["slope"] * 1e3))
        print("-- symplectic check   R33 R44 - R34 R43 = %.4f" % det)

    extras = [(b0_zdc_summary(rpath, p_beam), rix)
              for (ipath, rpath), (rix, _rows) in zip(pairs, idx)]
    print("-- B0 and ZDC on the same scan")
    layers = next((e["B0"]["layers"] for e, _r in extras
                   if e.get("B0") and e["B0"]["rows"]), None)
    if layers:
        nb0 = sum(len(e["B0"]["rows"]) for e, _r in extras if e.get("B0"))
        print("   B0 tracker layers at z = %s mm; primary in layer 1 in "
              "%d events" % (", ".join("%.0f" % z for z in layers), nb0))
        for nm, lg, coord, key in (("B0 R12", "b0x", "x", 0),
                                   ("B0 R34", "b0y", "y", 1)):
            pr = [(r, e["B0"]["rows"][r["ievt"]]) for e, rix in extras
                  if e.get("B0") for r in rix
                  if r["leg"] == lg and r["ievt"] in e["B0"]["rows"]]
            f = robust_fit([th_ip(r)[key] for r, _h in pr],
                           [h[coord] for _r, h in pr])
            print("   " + _fmt(nm, "[m]", f, 1e-3))
            out[nm.replace(" ", "_")] = f
    else:
        print("   B0: no primary reached the tracker")
    for c in ("HcalFarForwardZDCHits", "EcalFarForwardZDCHits"):
        legs, tot = {}, 0
        for e, rix in extras:
            rws = e["ZDC"].get(c, {})
            tot += len(rws)
            for r in rix:
                if r["ievt"] in rws:
                    legs.setdefault(r["leg"], []).append(rws[r["ievt"]])
        if not tot:
            continue
        print("   %-24s %3d events; by leg: %s" % (
            c, tot, ", ".join(
                "%s %d (<E> %.1f GeV, <x> %.0f mm)"
                % (k, len(v), np.mean([q["e"] for q in v]),
                   np.mean([q["x"] for q in v]))
                for k, v in sorted(legs.items()))))
    return out


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--config", default="18x275",
                   choices=sorted(REFERENCE_MOMENTUM))
    b.add_argument("--out", required=True)
    b.add_argument("--pdg", type=int, default=1000030060)
    b.add_argument("--stem", default="scan",
                   help="output stem; `fit` pools every scan*.index / "
                        "scan*.edm4hep.root pair in the directory")
    b.add_argument("--legs", default=None,
                   help="comma list, e.g. x0,y0 (default: every leg)")
    b.add_argument("--topup", action="store_true",
                   help="use OFFSETS_MM_TOPUP for the x0/y0 legs")
    f = sub.add_parser("fit")
    f.add_argument("outdir")
    f.add_argument("--config", default="18x275",
                   choices=sorted(REFERENCE_MOMENTUM))
    f.add_argument("--root", default=None)
    f.add_argument("--plane", default="S1L1",
                   choices=[p[0] for p in STATION_PLANES])
    args = ap.parse_args()
    if args.cmd == "build":
        legs = set(args.legs.split(",")) if args.legs else None
        n = write_scan(args.out, args.config, args.pdg, args.stem, legs,
                       args.topup)
        print("# %s/%s.hepmc: %d events, %s, 6Li at %g GeV "
              "(ring reference rigidity)"
              % (args.out, args.stem, n, args.config,
                 REFERENCE_MOMENTUM[args.config]))
        print("# next, INSIDE the container and through a ZERO-INSERTION "
              "geometry (tools/fullsim/ff_zero_insertion.sh):")
        print("#   npsim --compactFile $S/zero/epic_craterlake_%s_zero.xml \\"
              % args.config)
        print("#         --inputFiles %s/%s.hepmc --numberOfEvents %d \\"
              % (args.out, args.stem, n))
        print("#         --physics.list FTFP_BERT "
              "--part.minimalKineticEnergy '100*MeV' \\")
        print("#         --outputFile %s/%s.edm4hep.root"
              % (args.out, args.stem))
        return
    fit_config(args.outdir, args.config, args.root, args.plane)


if __name__ == "__main__":
    main()
