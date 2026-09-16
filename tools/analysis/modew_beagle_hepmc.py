#!/usr/bin/env python3
"""Mode-W weights on an official BeAGLE sample, written back as HepMC3.

Reads an (optionally remote, xrootd) BeAGLE `*.hepmc3.tree.root` file,
forms the DIS invariants of every event from the record itself, evaluates
`polligen.reweight.ModeWReweighter` for ONE spin category of a run plan,
and writes the same events out as HepMC3 ASCII with the weight NAMED in the
event's weights vector -- first in it by default, which is the placement
npsim preserves (see below).  That file is what
abconv/npsim/EICrecon consume, so the Mode-W injection of plans/05 step 5.C
travels with the event into the reconstruction (step 5.D).

Runs INSIDE the eic-shell image, whose working directory must sit two
levels below $ROOTSYS or ROOT never finds its module map (see
tools/fullsim/README.md, "Known issue"):

    SIF=~/Projects/eic-2026/local/lib/eic_xl-nightly.sif
    B=root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/DIS/BeAGLE1.03.02-3.1
    singularity exec --pwd /opt/local/lib/root $SIF python3 \\
        <repo>/tools/analysis/modew_beagle_hepmc.py \\
        $B/eH2/en/9x130/q2_1to1000/<file>.hepmc3.tree.root \\
        $S/modew_100.hepmc --nevents 100 --csv $S/modew_100.csv

Every path must be absolute, because of the `--pwd`.

Two things about this pyHepMC3 build, both measured 2026-09-16:

* `GenEvent::weights()` is bound BY COPY -- `ev.weights().append(w)` is
  silently dropped.  The weights vector is reachable only through the
  `GenEventData` round trip, `ev.write_data(d)` -> `d.weights = ...` ->
  `ev.read_data(d)`, which is what `_with_weight_data` below does.
* `GenRunInfo::set_weight_names` takes a `std::vector<std::string>`, not a
  Python list; build it with `pyHepMC3.pyHepMC3.std.vector_std_string`.

And one about the chain, measured the same day: npsim (DD4hep 1.37) copies
**only `weights()[0]`** into the edm4hep `EventHeader.weight` and leaves the
`EventHeader.weights` vector member empty, so a weight APPENDED to a
BeAGLE file's `default` weight reaches the reconstructed file as 1.0 and a
weight written FIRST reaches it bit for bit.  That is why
`--weight-placement` defaults to `first`: the Mode-W weight becomes the
sample's nominal weight (correct for a reweighted sample -- the BeAGLE
files are unweighted, `default` = 1) and the original moves to index 1.

The kinematics are the ones `reweight.sample_from_beagle_rows` forms from
the CSV of `dump_spectators.py`, with one improvement the CSV cannot carry:
the beam LEPTON is taken from the event's own status-4 row, so the
afterburned electron direction and energy enter q = k - k' per event rather
than as an assumed (0, 0, -E, E).
"""

import argparse
import csv
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))   # tools/analysis -> repo

from pyHepMC3 import HepMC3                      # noqa: E402
from pyHepMC3 import pyHepMC3 as _p              # noqa: E402
from pyHepMC3.rootIO import HepMC3 as rootIO     # noqa: E402

# BeAGLE's own truth attributes, carried through to the CSV so the
# reco-level residuals of modew_reco_readback.py have a reference that does
# not depend on this script's own reconstruction of the invariants.
TRUTH_ATTRS = ("trueX", "trueQ2", "trueY", "trueW2", "x", "y", "QSquared")

# `x`, `y` and `QSquared` are also names this script writes for its OWN
# reconstruction of the invariants, and a CSV with a repeated header silently
# loses one of the two to `csv.DictReader` (it keeps the last).  The BeAGLE
# attributes are therefore written under a prefixed column name; `trueX` and
# the rest, which modew_reco_readback.py reads, keep theirs.
TRUTH_COLUMNS = tuple(("beagle_" + a if a in ("x", "y", "QSquared") else a)
                      for a in TRUTH_ATTRS)


def _sys_path(repo):
    for sub in ("evgen", "fastsim"):
        p = os.path.join(repo, sub)
        if p not in sys.path:
            sys.path.insert(0, p)


def _minkowski(a, b):
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]


def _four(mom):
    return np.array([mom.e(), mom.px(), mom.py(), mom.pz()], dtype=float)


def event_kinematics(evt, e_beam_fallback=None):
    """(x, Q2, y, phi, beam dict) from one BeAGLE GenEvent, or None.

    The ion-side beam is the forward-going status-4 hadron (the BeAGLE eH2
    files record the STRUCK NUCLEON there, so x is the nucleon's Bjorken x,
    which is the variable `polli_fastsim`'s per-nucleon structure functions
    take).  The scattered electron is the highest-energy status-1 e-.
    """
    kbeam = pbeam = None
    scat = None
    for p in evt.particles():
        st, pid = p.status(), p.pid()
        if st == 4:
            if abs(pid) == 11:
                if kbeam is None or abs(p.momentum().pz()) > abs(kbeam.pz()):
                    kbeam = p.momentum()
            elif pbeam is None or p.momentum().pz() > pbeam.pz():
                pbeam = p.momentum()
        elif st == 1 and pid == 11:
            if scat is None or p.momentum().e() > scat.e():
                scat = p.momentum()
    if pbeam is None or scat is None:
        return None
    if kbeam is None:
        if e_beam_fallback is None:
            return None
        k = np.array([e_beam_fallback, 0.0, 0.0, -e_beam_fallback])
    else:
        k = _four(kbeam)
    p4 = _four(pbeam)
    kp = _four(scat)
    q = k - kp
    q2 = -_minkowski(q, q)
    pq = _minkowski(p4, q)
    pk = _minkowski(p4, k)
    if not (q2 > 0.0 and pq > 0.0 and pk > 0.0):
        return None
    return dict(x=q2 / (2.0 * pq), q2=q2, y=pq / pk,
                phi=math.atan2(kp[2], kp[1]) % (2.0 * math.pi),
                e_beam=k[0], pz_ion=p4[3], e_prime=kp[0],
                theta_e=math.atan2(math.hypot(kp[1], kp[2]), kp[3]))


def _attr(evt, name):
    try:
        s = evt.attribute_as_string(name)
    except Exception:
        return float("nan")
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("input", help="BeAGLE hepmc3.tree.root (local or xrootd)")
    ap.add_argument("output", help="HepMC3 ASCII to write")
    ap.add_argument("--csv", help="per-event kinematics + weight + truth")
    ap.add_argument("--nevents", type=int, default=100,
                    help="events to WRITE (default 100)")
    ap.add_argument("--scan", type=int, default=4000,
                    help="max events to READ while filling --nevents")
    ap.add_argument("--e-beam", type=float, default=None,
                    help="beam electron energy [GeV] if the record has no "
                         "status-4 lepton (it normally does)")
    ap.add_argument("--nucleus", default="DEUTERON",
                    help="polli_fastsim.beams name (default DEUTERON)")
    ap.add_argument("--plan", default="tensor_thirds",
                    choices=("tensor_thirds", "helicity_flip"))
    ap.add_argument("--category", default="azz0",
                    help="which category of the plan to weight (default "
                         "azz0, the m0-enriched tensor third)")
    ap.add_argument("--pz", type=float, default=0.6)
    ap.add_argument("--pzz", type=float, default=0.6)
    ap.add_argument("--pe", type=float, default=0.7)
    ap.add_argument("--weight-name", default=None,
                    help="name of the Mode-W weight (default modew_<cat>)")
    ap.add_argument("--weight-placement", default="first",
                    choices=("first", "append"),
                    help="first (default): the Mode-W weight becomes weights"
                         "[0], the nominal one, and the sample's own weight "
                         "moves to index 1.  append: it goes last.  npsim "
                         "keeps ONLY weights[0] (see below), so `first` is "
                         "the placement that survives the chain.")
    ap.add_argument("--no-scenario", action="store_true",
                    help="keep every event instead of the fom.Scenario window")
    ap.add_argument("--repo", default=_REPO)
    args = ap.parse_args(argv)

    _sys_path(args.repo)
    from polligen import bookkeeping as bk
    from polligen import reweight as rw
    from polligen.xsec import InclusiveKernel
    from polli_fastsim import beams, fom
    from polli_fastsim.polarized import toy_b1

    reader = rootIO.ReaderRootTree(args.input)
    if reader.failed():
        sys.exit("cannot open %s" % args.input)

    # --- pass: read, keep the events inside the window --------------------
    kin, data, truth = [], [], []
    evt = HepMC3.GenEvent()
    run_info = None
    nread = 0
    scenario = fom.Scenario()
    while len(kin) < args.nevents and nread < args.scan:
        if not reader.read_event(evt) or reader.failed():
            break
        nread += 1
        if run_info is None:
            run_info = evt.run_info()
        k = event_kinematics(evt, args.e_beam)
        if k is not None:
            s_evt = 4.0 * k["e_beam"] * abs(k["pz_ion"])
            if args.no_scenario or _in_window(k, s_evt, scenario, rw):
                kin.append(k)
                truth.append({a: _attr(evt, a) for a in TRUTH_ATTRS})
                d = HepMC3.GenEventData()
                evt.write_data(d)
                data.append(d)
        evt.clear()
    reader.close()
    if not kin:
        sys.exit("no usable events in the first %d read" % nread)
    print("read %d events, kept %d" % (nread, len(kin)), flush=True)

    x = np.array([k["x"] for k in kin])
    q2 = np.array([k["q2"] for k in kin])
    phi = np.array([k["phi"] for k in kin])
    e_beam = float(np.mean([k["e_beam"] for k in kin]))
    pz_ion = float(np.mean([abs(k["pz_ion"]) for k in kin]))
    s = 4.0 * e_beam * pz_ion
    print("s = %.2f GeV^2 (E_e = %.3f, p_z,ion = %.3f)" % (s, e_beam, pz_ion))

    pool = rw.ExternalSample(x, q2, phi, s)
    kernel = InclusiveKernel(getattr(beams, args.nucleus), b1_func=toy_b1)
    if args.plan == "tensor_thirds":
        plan = bk.tensor_thirds_plan(args.pz, args.pzz)
    else:
        plan = bk.helicity_flip_plan(1.0, args.pz, args.pe)
    cats = {c.name: c for c in plan.categories}
    if args.category not in cats:
        sys.exit("category %s not in plan (%s)"
                 % (args.category, ", ".join(sorted(cats))))
    cat = cats[args.category]
    rwt = rw.ModeWReweighter(kernel, pool)
    w = np.asarray(rwt.event_weights(cat), dtype=float)
    mean_w = rwt.mean_weight(cat)
    print("category %s: <W> = %.8f, W in [%.6f, %.6f]"
          % (cat.name, mean_w, w.min(), w.max()))
    print("injected pool truth: A_zz = %.6g, A_par = %.6g"
          % (rwt.truth_azz(), rwt.truth_a_parallel()))

    # --- write ------------------------------------------------------------
    wname = args.weight_name or ("modew_" + cat.name)
    if run_info is None:
        run_info = HepMC3.GenRunInfo()
    names = [str(n) for n in run_info.weight_names()] or ["default"]
    insert = False
    if wname in names:
        idx = names.index(wname)
    elif args.weight_placement == "first":
        idx, insert = 0, True
        names = [wname] + names
    else:
        idx = len(names)
        names = names + [wname]
    run_info.set_weight_names(_p.std.vector_std_string(names))
    for key, val in (("modew_category", cat.name),
                     ("modew_weight_name", wname),
                     ("modew_nucleus", args.nucleus),
                     ("modew_plan", args.plan),
                     ("modew_pz", "%.6f" % args.pz),
                     ("modew_pzz", "%.6f" % args.pzz),
                     ("modew_mean_weight", "%.10f" % mean_w),
                     ("modew_s_gev2", "%.6f" % s),
                     ("modew_source", os.path.basename(args.input))):
        run_info.add_attribute(key, HepMC3.StringAttribute(val))

    writer = HepMC3.WriterAscii(args.output, run_info)
    out = HepMC3.GenEvent()
    out.set_run_info(run_info)
    for i, d in enumerate(data):
        _with_weight_data(d, float(w[i]), idx, insert=insert)
        out.read_data(d)
        out.set_run_info(run_info)
        writer.write_event(out)
        out.clear()
    writer.close()
    print("wrote %d events to %s (weight '%s' at index %d of %s)"
          % (len(data), args.output, wname, idx, names))

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            cw = csv.writer(fh)
            cw.writerow(["iout", "x", "q2", "y", "phi", "e_prime", "theta_e",
                         "weight"] + list(TRUTH_COLUMNS))
            for i, k in enumerate(kin):
                cw.writerow(["%d" % i] +
                            ["%.10g" % k[f] for f in ("x", "q2", "y", "phi",
                                                      "e_prime", "theta_e")] +
                            ["%.10g" % w[i]] +
                            ["%.10g" % truth[i][a] for a in TRUTH_ATTRS])
        print("wrote %s" % args.csv)
    return 0


def _with_weight_data(d, weight, index, insert=False):
    """Put `weight` at `index` of a GenEventData's weights vector.  The
    GenEventData round trip is the only writable path: `GenEvent::weights()`
    is bound by COPY in this pyHepMC3 build, so appending to it is dropped."""
    w = list(d.weights)
    if insert:
        w.insert(index, float(weight))
    else:
        if len(w) <= index:
            w = w + [1.0] * (index + 1 - len(w))
        w[index] = float(weight)
    d.weights = _p.std.vector_double(w)
    return d


def _in_window(k, s, scenario, rw):
    """The `fom.Scenario` window, event by event -- the same cuts
    `ExternalSample.in_scenario` applies, on a single event."""
    pool = rw.ExternalSample(np.array([k["x"]]), np.array([k["q2"]]),
                             np.array([k["phi"]]), s)
    return bool(pool.in_scenario(scenario, k["e_beam"])[0])


if __name__ == "__main__":
    sys.exit(main())
