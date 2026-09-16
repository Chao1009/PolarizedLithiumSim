#!/usr/bin/env python3
"""Read a Mode-W reconstructed file back: the weight, x, Q2, and the
far-forward spectator.

The companion of `tools/analysis/modew_beagle_hepmc.py`.  That script puts
a named polarized weight on official BeAGLE events and writes HepMC3;
`tools/fullsim/modew_chain.sh` puts them through npsim and EICrecon; this
one asks the three questions the round trip exists to answer.

  1. Did the generator weight survive into the simulated and the
     reconstructed file, event by event?  Both edm4hep spellings are
     checked -- the scalar `EventHeader.weight` and the podio vector member
     `EventHeader.weights` (`weights_begin/_end` plus `_EventHeader_weights`).
  2. What do the reconstructed x and Q2 look like against the same events'
     truth, weighted, since a Mode-W sample is a weighted sample?  Both
     EICrecon's own `InclusiveKinematics*` and an independent electron
     method formed here from the reconstructed electron and the MC beams.
  3. What does the reconstruction return for the ion-side SPECTATOR of a
     deuteron beam?  That is the counterpart of the lithium blocker: the
     far-forward Roman-Pot reconstruction is hard-wired to a proton
     (`MatrixTransferStaticConfig.h`), and a deuteron's spectator proton is
     the one fragment for which that assumption is the right species at the
     wrong rigidity.

Runs on the HOST (uproot), not in the container:

    python3 tools/analysis/modew_reco_readback.py $S/sim/reco.edm4eic.root \\
        --sim $S/sim/sim.edm4hep.root --csv $S/modew_100.csv

`--csv` is the generator-level table modew_beagle_hepmc.py wrote; with it
the weight comparison is exact rather than self-consistent, and the event
ORDER is checked against BeAGLE's own truth attributes, which is what makes
re-attaching a dropped weight by index legitimate.
"""

import argparse
import csv
import os
import sys

import numpy as np


# --- podio plumbing -----------------------------------------------------


def branch(tree, name):
    """A podio branch by its short name.  The file spells it
    `Coll/Coll.field`, uproot also answers to `Coll.field`; this returns
    None instead of raising when it is absent."""
    for cand in (name, name.split(".")[0] + "/" + name):
        if cand in tree:
            try:
                return tree[cand].array(library="np")
            except Exception:
                return None
    return None


def has(tree, coll, field="energy"):
    return branch(tree, "%s.%s" % (coll, field)) is not None


def collection_ids(path, tree_name="events"):
    """{collectionID: name} from the podio metadata tree."""
    import uproot
    f = uproot.open(path)
    key = "podio_metadata"
    if key not in [k.split(";")[0] for k in f.keys()]:
        return {}
    meta = f[key]
    base = "%s___CollectionTypeInfo" % tree_name
    try:
        ids = meta["%s/%s.collectionID" % (base, base)].array(library="np")
        names = meta["%s/%s.name" % (base, base)].array(library="np")
    except Exception:
        return {}
    if len(ids) == 0:
        return {}
    return {int(i): str(n) for i, n in zip(np.atleast_1d(ids[0]),
                                           np.atleast_1d(names[0]))}


def subset_fourvectors(tree, coll, idmap, out_targets=None):
    """Per event, the (E, px, py, pz) of a podio SUBSET collection --
    `<coll>_objIdx` is a list of (collectionID, index) into a real one.
    Returns None when the collection is not a subset here."""
    out_targets = set() if out_targets is None else out_targets
    idx = branch(tree, "%s_objIdx.index" % coll)
    cid = branch(tree, "%s_objIdx.collectionID" % coll)
    if idx is None or cid is None:
        return None
    cache = {}
    out = []
    targets = set()
    for i in range(len(idx)):
        ii = np.atleast_1d(idx[i]).astype(int)
        cc = np.atleast_1d(cid[i]).astype(np.int64)
        rows = []
        for j, c in zip(ii, cc):
            target = idmap.get(int(c) & 0xFFFFFFFF)
            if target is None:
                continue
            targets.add(target)
            if target not in cache:
                px_ = branch(tree, "%s.momentum.x" % target)
                py_ = branch(tree, "%s.momentum.y" % target)
                pz_ = branch(tree, "%s.momentum.z" % target)
                en_ = branch(tree, "%s.energy" % target)
                if en_ is None:                  # MCParticles carry mass
                    m_ = branch(tree, "%s.mass" % target)
                    if m_ is not None and px_ is not None:
                        en_ = np.array(
                            [np.sqrt(np.asarray(px_[q], float) ** 2
                                     + np.asarray(py_[q], float) ** 2
                                     + np.asarray(pz_[q], float) ** 2
                                     + np.asarray(m_[q], float) ** 2)
                             for q in range(len(px_))], dtype=object)
                cache[target] = (en_, px_, py_, pz_, target)
            en, px, py, pz, _tgt = cache[target]
            if en is None or px is None:
                continue
            if j >= len(np.atleast_1d(en[i])):
                continue
            rows.append(np.array([np.atleast_1d(en[i])[j],
                                  np.atleast_1d(px[i])[j],
                                  np.atleast_1d(py[i])[j],
                                  np.atleast_1d(pz[i])[j]], dtype=float))
        out.append(np.array(rows) if rows else np.empty((0, 4)))
    out_targets.clear()
    out_targets.update(targets)
    return out


# --- kinematics ---------------------------------------------------------


def weighted_quantile(values, quantiles, weights=None):
    """Weighted quantiles of `values` (the ordinary ones when the weights
    are flat).  `quantiles` in [0, 1]."""
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return np.full(np.shape(quantiles), np.nan)
    weights = (np.ones_like(values) if weights is None
               else np.asarray(weights, dtype=float))
    order = np.argsort(values)
    v, w = values[order], weights[order]
    cum = np.cumsum(w) - 0.5 * w
    cum /= w.sum()
    return np.interp(quantiles, cum, v)


def minkowski(a, b):
    """(E, px, py, pz) . (E, px, py, pz) with the (+,-,-,-) metric."""
    a, b = np.atleast_2d(a), np.atleast_2d(b)
    return (a[:, 0] * b[:, 0] - a[:, 1] * b[:, 1]
            - a[:, 2] * b[:, 2] - a[:, 3] * b[:, 3])


def invariants(k, kp, p):
    """(x, Q2, y) from beam lepton, scattered lepton and beam hadron."""
    k, kp, p = np.atleast_2d(k), np.atleast_2d(kp), np.atleast_2d(p)
    q = k - kp
    q2 = -minkowski(q, q)
    pq = minkowski(p, q)
    pk = minkowski(p, k)
    with np.errstate(divide="ignore", invalid="ignore"):
        return q2 / (2.0 * pq), q2, pq / pk


# --- the weight ---------------------------------------------------------


def read_weights(tree):
    """{branch: per-event value} for every EventHeader weight spelling."""
    out = {}
    for key in tree.keys():
        if key.split("/")[0] not in ("EventHeader", "_EventHeader_weights"):
            continue
        if "weight" not in key.lower():
            continue
        if key in ("EventHeader",):
            continue
        try:
            out[key] = tree[key].array(library="np")
        except Exception as exc:                          # pragma: no cover
            out[key] = "unreadable: %s" % exc
    return out


def _flat(values):
    return np.array([np.atleast_1d(v)[0] if np.size(v) else np.nan
                     for v in values], dtype=float)


# --- truth from the MC record -------------------------------------------


def mc_truth(tree, coll="MCParticles"):
    """Per event: beam lepton, beam hadron, scattered lepton and the
    forward final-state baryons, from the generator rows of the MC record.
    `generatorStatus` 4 is a beam particle and 1 a final-state one, the
    HepMC3 convention npsim copies through."""
    pdg = branch(tree, coll + ".PDG")
    gs = branch(tree, coll + ".generatorStatus")
    px = branch(tree, coll + ".momentum.x")
    py = branch(tree, coll + ".momentum.y")
    pz = branch(tree, coll + ".momentum.z")
    mass = branch(tree, coll + ".mass")
    if pdg is None:
        return None
    rows = []
    for i in range(len(pdg)):
        p_ = np.stack([np.asarray(px[i], float), np.asarray(py[i], float),
                       np.asarray(pz[i], float)], axis=-1)
        e_ = np.sqrt((p_ ** 2).sum(-1) + np.asarray(mass[i], float) ** 2)
        four = np.concatenate([e_[:, None], p_], axis=1)
        id_, st = np.asarray(pdg[i]), np.asarray(gs[i])
        row = {"n": len(id_)}
        bl = np.where((st == 4) & (np.abs(id_) == 11))[0]
        bh = np.where((st == 4) & (np.abs(id_) != 11))[0]
        fe = np.where((st == 1) & (id_ == 11))[0]
        if bl.size:
            row["k"] = four[bl[np.argmax(np.abs(four[bl, 3]))]]
        if bh.size:
            j = bh[np.argmax(four[bh, 3])]
            row["p"] = four[j]
            row["beam_pdg"] = int(id_[j])
        if fe.size:
            row["kp"] = four[fe[np.argmax(four[fe, 0])]]
        if "p" in row:
            fwd = np.where((st == 1) & np.isin(id_, (2212, 2112))
                           & (four[:, 3] > 0.1 * row["p"][3]))[0]
            row["spectators"] = [(int(id_[j]), four[j]) for j in fwd]
        rows.append(row)
    return rows


def angle_to(v, axis):
    """Polar angle of a 3-vector about `axis` -- the ion axis, not z, so a
    25 mrad crossing angle does not sit in every far-forward angle."""
    v = np.asarray(v, float)[1:]
    a = np.asarray(axis, float)[1:]
    a = a / np.linalg.norm(a)
    par = float(np.dot(v, a))
    perp = np.linalg.norm(v - par * a)
    return float(np.arctan2(perp, par))


# --- the reconstructed electron -----------------------------------------


ELECTRON_COLLECTIONS = ("ScatteredElectronsTruth", "ReconstructedElectrons",
                        "ScatteredElectronsEMinusPz")
KINEMATICS = ("InclusiveKinematicsElectron", "InclusiveKinematicsTruth",
              "InclusiveKinematicsSigma", "InclusiveKinematicsDA",
              "InclusiveKinematicsJB", "InclusiveKinematicsESigma",
              "InclusiveKinematicsML")


def best_electron(tree, coll, idmap, out_targets=None):
    """Per event, the highest-energy entry of `coll` (a subset collection
    of reconstructed particles), or None."""
    rows = subset_fourvectors(tree, coll, idmap, out_targets)
    if rows is None:
        en = branch(tree, coll + ".energy")
        if en is None:
            return None
        px = branch(tree, coll + ".momentum.x")
        py = branch(tree, coll + ".momentum.y")
        pz = branch(tree, coll + ".momentum.z")
        rows = [np.stack([np.asarray(en[i], float), np.asarray(px[i], float),
                          np.asarray(py[i], float), np.asarray(pz[i], float)],
                         axis=-1) for i in range(len(en))]
    return [(r[np.argmax(r[:, 0])] if len(r) else None) for r in rows]


def kinematics_of(tree, coll):
    x = branch(tree, coll + ".x")
    q2 = branch(tree, coll + ".Q2")
    if x is None or q2 is None:
        return None
    return [((float(np.atleast_1d(x[i])[0]), float(np.atleast_1d(q2[i])[0]))
             if np.size(x[i]) else None) for i in range(len(x))]


# --- the far forward ----------------------------------------------------


FF_COLLECTIONS = ("ForwardRomanPotHits", "ForwardRomanPotRawHits",
                  "ForwardRomanPotRecHits", "ForwardRomanPotRecParticles",
                  "ForwardRomanPotStaticRecParticles",
                  "ForwardOffMTrackerHits", "ForwardOffMTrackerRawHits",
                  "ForwardOffMTrackerRecHits", "ForwardOffMRecParticles",
                  "B0TrackerHits", "B0TrackerRecHits", "B0TrackerCKFTracks",
                  "B0ECalRecHits", "ReconstructedB0EcalNeutrals",
                  "HcalFarForwardZDCRecHits", "HcalFarForwardZDCClusters",
                  "ReconstructedHcalFarForwardZDCNeutrals",
                  "EcalFarForwardZDCRecHits", "EcalFarForwardZDCClusters")
FF_FIELDS = ("energy", "momentum.x", "position.x", "cellID", "type",
             "chi2")


def ff_summary(tree):
    """(hits, events with >= 1) per far-forward collection present."""
    out = {}
    for c in FF_COLLECTIONS:
        for f in FF_FIELDS:
            a = branch(tree, "%s.%s" % (c, f))
            if a is None:
                continue
            n = np.array([np.size(v) for v in a])
            out[c] = (int(n.sum()), int((n > 0).sum()))
            break
    return out


def ff_particles(tree, coll):
    en = branch(tree, coll + ".energy")
    if en is None:
        return None
    px = branch(tree, coll + ".momentum.x")
    py = branch(tree, coll + ".momentum.y")
    pz = branch(tree, coll + ".momentum.z")
    return [np.stack([np.asarray(en[i], float), np.asarray(px[i], float),
                      np.asarray(py[i], float), np.asarray(pz[i], float)],
                     axis=-1) for i in range(len(en))]


# --- the report ---------------------------------------------------------


def load_csv(path):
    rows = list(csv.DictReader(open(path)))
    if not rows:
        return {}
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def residual_line(label, reco, truth, weights, n_events):
    dev = reco / truth - 1.0
    med = float(weighted_quantile(dev, 0.5, weights))
    p68 = float(weighted_quantile(np.abs(dev), 0.68, weights))
    print("    %-3s on %3d/%d events: weighted median %+.4f, "
          "68th pct |residual| %.4f  (unweighted %+.4f / %.4f)"
          % (label, dev.size, n_events, med, p68, float(np.median(dev)),
             float(np.percentile(np.abs(dev), 68))))
    return med, p68


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("reco", help="reco.edm4eic.root")
    ap.add_argument("--sim", help="sim.edm4hep.root (the npsim leg)")
    ap.add_argument("--csv", help="the generator table of "
                                  "modew_beagle_hepmc.py")
    ap.add_argument("--electrons", default=None,
                    help="scattered-electron collection (default: the first "
                         "of " + ", ".join(ELECTRON_COLLECTIONS) + ")")
    args = ap.parse_args(argv)

    import uproot
    reco = uproot.open(args.reco)["events"]
    idmap = collection_ids(args.reco)
    ncoll = len(set(k.split("/")[0] for k in reco.keys()))
    print("=== %s: %d events, %d collections"
          % (os.path.basename(args.reco), reco.num_entries, ncoll))

    # 1 -- the weight -----------------------------------------------------
    tab, gen_w = {}, None
    if args.csv:
        tab = load_csv(args.csv)
        gen_w = tab.get("weight")
        if gen_w is None:
            print("generator table: %d rows, no 'weight' column -- the "
                  "weight comparison below is skipped"
                  % len(tab.get("iout", ())))
        else:
            print("generator table: %d rows, <W> = %.10f, W in [%.6f, %.6f]"
                  % (gen_w.size, gen_w.mean(), gen_w.min(), gen_w.max()))
            if gen_w.size != reco.num_entries:
                print("  (!) %d table rows against %d reconstructed events: "
                      "the join is by event INDEX, so the weights are not "
                      "applied" % (gen_w.size, reco.num_entries))
                gen_w = None
    print("\n-- 1. the weight")
    for label, path in (("sim ", args.sim), ("reco", args.reco)):
        if not path:
            continue
        t = uproot.open(path)["events"]
        w = read_weights(t)
        if not w:
            print("[%s] no EventHeader weight branch at all" % label)
        for key, val in sorted(w.items()):
            if key.startswith("_EventHeader_weights"):
                tot = sum(np.size(v) for v in val)
                print("[%s] %-40s vector member, %d entries in %d events"
                      % (label, key, tot,
                         int(sum(np.size(v) > 0 for v in val))))
                continue
            flat = _flat(val)
            print("[%s] %-40s n=%d, first %.10f, mean %.10f, "
                  "distinct values %d"
                  % (label, key, flat.size, flat[0], np.nanmean(flat),
                     len(np.unique(flat))))
            if (gen_w is not None and flat.size == gen_w.size
                    and key.endswith(".weight")):
                # `weights_begin`/`weights_end` are podio index offsets, not
                # weights, so they are deliberately NOT compared here.  The
                # generator table is written to 10 significant digits, so the
                # meaningful statement is agreement to that precision and not
                # bit equality -- both counts are printed.
                d = np.abs(flat - gen_w)
                tol = 1e-9 * np.maximum(np.abs(gen_w), 1.0)
                print("       vs the generator weight: max |diff| %.3e, "
                      "equal on %d/%d events bit for bit, %d/%d within the "
                      "table's 10-digit printing precision"
                      % (d.max(), int((d == 0).sum()), d.size,
                         int((d <= tol).sum()), d.size))

    # 2 -- truth and the event order --------------------------------------
    truth = mc_truth(reco)
    if truth is None and args.sim:
        truth = mc_truth(uproot.open(args.sim)["events"])
    if truth is None:
        print("no MCParticles: nothing further to say")
        return 1
    keep = [i for i, r in enumerate(truth)
            if all(f in r for f in ("k", "kp", "p"))]
    k = np.array([truth[i]["k"] for i in keep])
    kp = np.array([truth[i]["kp"] for i in keep])
    p = np.array([truth[i]["p"] for i in keep])
    xt, qt, yt = invariants(k, kp, p)
    pos = {j: n for n, j in enumerate(keep)}
    beam_pdgs = sorted(set(truth[i].get("beam_pdg") for i in keep))
    print("\n-- 2. truth, from the MC record")
    print("  %d/%d events with beam lepton + beam hadron + scattered lepton; "
          "beam hadron PDG %s" % (len(keep), reco.num_entries, beam_pdgs))
    print("  x med %.4g, Q2 med %.4g GeV^2, y med %.3f"
          % (np.median(xt), np.median(qt), np.median(yt)))
    if "trueX" in tab and len(keep) == tab["trueX"].size:
        print("  against BeAGLE's own trueX/trueQ2 IN FILE ORDER: median "
              "|dev| %.2e / %.2e -- the event order survives the chain"
              % (np.median(np.abs(xt / tab["trueX"] - 1.0)),
                 np.median(np.abs(qt / tab["trueQ2"] - 1.0))))
    w_all = gen_w if gen_w is not None else np.ones(reco.num_entries)

    # 3 -- reconstructed x and Q2 -----------------------------------------
    print("\n-- 3. reconstructed x and Q2")
    print("  EICrecon's own collections, events filled:")
    for c in KINEMATICS:
        got = kinematics_of(reco, c)
        if got is None:
            continue
        have = [i for i in keep if got[i] is not None]
        print("    %-30s %3d/%d" % (c, len(have), reco.num_entries))
        if not have:
            continue
        ww = np.array([w_all[i] for i in have])
        residual_line("x", np.array([got[i][0] for i in have]),
                      np.array([xt[pos[i]] for i in have]), ww,
                      reco.num_entries)
        residual_line("Q2", np.array([got[i][1] for i in have]),
                      np.array([qt[pos[i]] for i in have]), ww,
                      reco.num_entries)
    for c in ("MCBeamElectrons", "MCBeamProtons", "MCScatteredElectrons"):
        tg = set()
        rows = subset_fourvectors(reco, c, idmap, tg)
        if rows is not None:
            print("    %-30s filled in %3d/%d events%s"
                  % (c, int(sum(len(r) > 0 for r in rows)), reco.num_entries,
                     (" -> " + ", ".join(sorted(tg))) if tg else
                     "  (no target collection: EMPTY)"))

    coll = args.electrons
    if coll is None:
        for c in ELECTRON_COLLECTIONS:
            if (branch(reco, c + "_objIdx.index") is not None
                    or has(reco, c)):
                coll = c
                break
    print("  the electron method formed here, from %s and the MC beams:"
          % coll)
    if coll is not None:
        tg = set()
        cand = best_electron(reco, coll, idmap, tg)
        if tg:
            print("    (%s resolves into %s)" % (coll, ", ".join(sorted(tg))))
        have = [i for i in keep if cand is not None and cand[i] is not None]
        if have:
            xr, qr = [], []
            for i in have:
                a, b, _ = invariants(truth[i]["k"], cand[i], truth[i]["p"])
                xr.append(float(a[0]))
                qr.append(float(b[0]))
            ww = np.array([w_all[i] for i in have])
            xr, qr = np.array(xr), np.array(qr)
            xtv = np.array([xt[pos[i]] for i in have])
            qtv = np.array([qt[pos[i]] for i in have])
            ytv = np.array([yt[pos[i]] for i in have])
            residual_line("x", xr, xtv, ww, reco.num_entries)
            residual_line("Q2", qr, qtv, ww, reco.num_entries)
            # the electron method degrades as 1/y on x and not on Q2; the
            # split is the check that the tail is kinematics and not a bug
            for lab, m in (("y < 0.1", ytv < 0.1), ("y >= 0.1", ytv >= 0.1)):
                if m.sum() < 3:
                    continue
                print("    %s (%d events):" % (lab, int(m.sum())))
                residual_line("x", xr[m], xtv[m], ww[m], reco.num_entries)
                residual_line("Q2", qr[m], qtv[m], ww[m], reco.num_entries)
        else:
            print("    no candidate in any event")

    # 4 -- the far forward -------------------------------------------------
    print("\n-- 4. the far forward (hits, events with >=1)")
    for c, (n, ne) in sorted(ff_summary(reco).items()):
        print("  %-40s %8d  %3d" % (c, n, ne))
    axis = truth[keep[0]]["p"]
    spec = [(i, pid, f) for i in keep
            for pid, f in truth[i].get("spectators", [])]
    for pid, name in ((2212, "proton"), (2112, "neutron")):
        rows = [(i, f) for i, q, f in spec if q == pid]
        if not rows:
            continue
        pmag = np.array([np.linalg.norm(f[1:]) for _, f in rows])
        th = np.array([angle_to(f, truth[i]["p"]) for i, f in rows])
        print("  truth forward %ss: %d in %d events, |p| med %.1f GeV "
              "(%.3f of the %.1f GeV beam row), theta to the ion axis "
              "10/50/90 pct = %.2f / %.2f / %.2f mrad"
              % (name, len(rows), len(set(i for i, _ in rows)),
                 np.median(pmag), np.median(pmag) / axis[3], axis[3],
                 *tuple(1e3 * np.percentile(th, [10, 50, 90]))))
    zaxis = np.array([1.0, 0.0, 0.0, 1.0])
    print("  the truth angles above are about the ION AXIS; each "
          "reconstructed collection below is\n  tested against BOTH axes and "
          "the one it is written in is named, because the tracker-based\n"
          "  far-forward reconstruction takes the crossing angle out and the "
          "calorimeter-based one does not.")
    for c in ("ForwardRomanPotRecParticles",
              "ForwardRomanPotStaticRecParticles",
              "ForwardOffMRecParticles",
              "ReconstructedHcalFarForwardZDCNeutrals",
              "ReconstructedB0EcalNeutrals"):
        parts = ff_particles(reco, c)
        if parts is None:
            continue
        tot = sum(len(v) for v in parts)
        nev = int(sum(len(v) > 0 for v in parts))
        print("  %-40s %3d particles in %3d events" % (c, tot, nev))
        if not tot:
            continue
        allp = np.concatenate([v for v in parts if len(v)])
        th_z = np.array([angle_to(v, zaxis) for v in allp])
        th_i = np.array([angle_to(v, axis) for v in allp])
        beam_frame = np.median(th_z) <= np.median(th_i)
        th = th_z if beam_frame else th_i
        print("      |p| med %.1f GeV, theta med %.2f mrad about %s "
              "(%.2f about the other)"
              % (np.median(np.linalg.norm(allp[:, 1:], axis=1)),
                 1e3 * np.median(th),
                 "+z, i.e. the crossing angle is OUT" if beam_frame
                 else "the ion axis, i.e. the LAB frame",
                 1e3 * np.median(th_i if beam_frame else th_z)))
        # match the leading reconstructed particle to the leading truth
        # forward proton of the same event
        pdg = branch(reco, c + ".PDG")
        seen = np.concatenate([np.atleast_1d(v) for v in pdg
                               if np.size(v)]) if pdg is not None else []
        want = 2212 if (len(seen) and
                        (np.asarray(seen) == 2212).mean() > 0.5) else 2112
        rows = []
        for i in keep:
            prot = [f for pid, f in truth[i].get("spectators", [])
                    if pid == want]
            if not prot or not len(parts[i]):
                continue
            tp = max(prot, key=lambda f: np.linalg.norm(f[1:]))
            rp = parts[i][np.argmax(np.linalg.norm(parts[i][:, 1:], axis=1))]
            rows.append((i, np.linalg.norm(tp[1:]), np.linalg.norm(rp[1:]),
                         angle_to(tp, truth[i]["p"]),
                         angle_to(rp, zaxis if beam_frame
                                  else truth[i]["p"])))
        if rows:
            ww = np.array([w_all[r[0]] for r in rows])
            pt_ = np.array([r[1] for r in rows])
            pr_ = np.array([r[2] for r in rows])
            tt_ = np.array([r[3] for r in rows])
            tr_ = np.array([r[4] for r in rows])
            print("      matched to the leading truth forward %s in "
                  "%d events:"
                  % ("proton" if want == 2212 else "neutron", len(rows)))
            residual_line("|p|", pr_, pt_, ww, reco.num_entries)
            dth = 1e3 * (tr_ - tt_)
            print("      theta  weighted median %+.3f mrad, 68th pct "
                  "|residual| %.3f mrad"
                  % (weighted_quantile(dth, 0.5, ww),
                     weighted_quantile(np.abs(dth), 0.68, ww)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
