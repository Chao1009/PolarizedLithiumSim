#!/usr/bin/env python3
"""WP2, the P_zz half of the delta-A scaling table (plans/07 WP2).

The fill-share half of that table is already written: `--run-share` /
`--lumi-fraction` scan the programme share f, the 1/sqrt(f) law is
verified against direct reruns rather than assumed, and Plans A, B and
A x B are priced in plans/07.  This is the other axis.  It scans the
tensor polarization over the same grid -- the three 6Li configurations
combined, the five x bins the fill-share half quotes delta-A_zz at, and
the same set of shares f = 1, 1/2, 1/3 -- and checks the projection
against the analytic law the work package names,

    delta A  =  sqrt(2) / (P_zz * D * A * sqrt(f * N)),

where N is the accepted event count at f = 1 summed over the cells that
enter the combination, D and A are the two amplitude attenuations of
`fom.Scenario` (both 1 in every published figure), and the sqrt(2) is the
equal-thirds three-state estimator of `asymmetries.err_azz`.  The ratio
column is the check: it is the measured delta-A divided by that
expression, and it must be 1.000 at every cell.

Nothing published moves here -- the script writes no figure and no file,
only this table.

TWO TENSOR OBSERVABLES, ONE TABLE.  `err_azz` and
`err_cos2phi_amplitude` are the same function of (N, P_zz) -- sqrt(2/N)
/ P_zz, the three-state estimator and the two-state cos 2phi fit landing
on the same variance -- so the table below serves the gluonometry
amplitude as well as A_zz.  That is asserted cell by cell rather than
assumed, and the assertion is printed; if the two ever part company the
line says so and the table stops being the whole story.

THE MIN-EVENTS FLOOR IS NOT THE LAW.  A cell enters the combination when
it holds at least `--min-events` events, and at f < 1 it holds f times
fewer, so a floor applied at the share being scanned SELECTS a different
set of cells and the measured ratio then departs from 1 through the
selection, not through the polarization.  `--floor reference` (the
default) fixes the cell set once, on the f = 1 projection, so the table
measures the law; `--floor cell` applies the floor at each share and
prints what the selection alone does.  The P_zz axis is untouched by
either: no event count depends on P_zz.

Usage:
    python3 scripts/wp2_pzz_table.py
    python3 scripts/wp2_pzz_table.py --pzz 0.6,0.8 --shares 1 --floor cell
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.inputs import get_backends
from polli_fastsim.polarized import toy_b1, toy_delta_gluon
from polli_fastsim.structure import NuclearF2

# the x bins the fill-share half of the table quotes delta-A_zz at
# (plans/07 WP2, Plan B): bin centres of the 40-bin log x grid
X_QUOTED = (0.0035, 0.0089, 0.0282, 0.2818, 0.5623)
# the shares that half tabulates: the whole year, half of it, a third
SHARES = (1.0, 0.5, 1.0 / 3.0)
# P_zz: the in-ring placeholder (0.60), the ECRP source target (0.80) and
# the span either side that the requirement is read off
PZZ = (0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)


def _floats(text):
    return tuple(float(v) for v in str(text).replace(",", " ").split())


def combined_errors(lumi, share, pzz, dilution, acceptance, base, g1m,
                    used=None, min_events=100.0):
    """delta A_zz per x bin, combined over Q2 and the three 6Li
    configurations, and the event count that goes with it.

    Returns (x centres, delta A per bin, N per bin, used masks).  `used`
    replays a cell set measured at another share (`--floor reference`);
    None applies the floor to this share's own counts."""
    inv2_tot = n_tot = x_c = None
    masks = []
    for i, cfg in enumerate(beams.default_configs("6Li")):
        sc = fom.Scenario(lumi_fb_per_nucleon=lumi, run_share=share,
                          pol_ion_tensor=pzz, dilution=dilution,
                          acceptance=acceptance)
        proj = fom.project_rates(cfg, sc,
                                 nuclear_f2=NuclearF2(cfg.ion, base=base))
        obs = fom.project_observables(cfg, sc, proj, g1m, toy_b1,
                                      toy_delta_gluon)
        use = (proj.accepted & (proj.n_events >= min_events) if used is None
               else used[i])
        masks.append(use)
        # err_azz and err_cos2phi_amplitude are the same function of
        # (N, P_zz); say so out loud rather than trusting it
        both = obs["err_azz"][use], obs["err_a_cos2phi"][use]
        if both[0].size and not np.array_equal(*both):
            raise AssertionError("err_azz and err_cos2phi_amplitude have "
                                 "parted company: this table no longer "
                                 "serves both tensor observables")
        inv2 = np.where(use, 1.0 / obs["err_azz"]**2, 0.0).sum(axis=1)
        n_b = np.where(use, proj.n_events, 0.0).sum(axis=1)
        x_c = proj.x[:, 0]
        inv2_tot = inv2 if inv2_tot is None else inv2_tot + inv2
        n_tot = n_b if n_tot is None else n_tot + n_b
    err = np.full(x_c.shape, np.inf)
    np.divide(1.0, np.sqrt(inv2_tot), out=err, where=inv2_tot > 0)
    return x_c, err, n_tot, masks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lumi", type=float, default=10.0,
                    help="programme luminosity [fb^-1/nucleon] (default: "
                         "%(default)s, one EIC year)")
    ap.add_argument("--pzz", default=",".join("%g" % p for p in PZZ),
                    help="tensor polarizations to scan (default: %(default)s)")
    ap.add_argument("--shares", default=",".join("%g" % s for s in SHARES),
                    help="run-plan shares to scan, the same grid the "
                         "fill-share half used (default: %(default)s)")
    ap.add_argument("--dilution", type=float, default=1.0,
                    help="measured amplitude dilution D (default 1.0, which "
                         "every published figure assumes)")
    ap.add_argument("--acceptance", type=float, default=1.0,
                    help="measured amplitude acceptance A (default 1.0)")
    ap.add_argument("--min-events", type=float, default=100.0,
                    dest="min_events",
                    help="cells below this count do not enter the "
                         "combination (default: %(default)s)")
    ap.add_argument("--floor", default="reference",
                    choices=["reference", "cell"],
                    help="where the min-events floor is applied: at the "
                         "share-1 projection ('reference', the default, so "
                         "the table measures the law) or at each share")
    ap.add_argument("--pdf", default="toy", choices=["toy", "grid"])
    ap.add_argument("--tol", type=float, default=2e-3,
                    help="the acceptance criterion on |ratio - 1| "
                         "(default: %(default)s)")
    args = ap.parse_args()
    pzz_grid, shares = _floats(args.pzz), _floats(args.shares)
    for name, value in (("--dilution", args.dilution),
                        ("--acceptance", args.acceptance)):
        if not value > 0:
            ap.error("%s must be positive" % name)
    if not all(s > 0 for s in shares) or not all(p > 0 for p in pzz_grid):
        ap.error("--shares and --pzz must be positive")

    backends = get_backends(args.pdf)
    base, g1m = backends["base"], backends["g1"]
    da = args.dilution * args.acceptance

    print("# wp2_pzz_table.py  pdf=%s  %s"
          % (args.pdf, fom.run_share_header(args.lumi, 1.0)))
    print("# shares scanned: %s (the compact table at the end is the "
          "share-1 one)" % ", ".join("%g" % s for s in shares))
    print("# delta A_zz per x bin, 3 configurations x Q2 combined, "
          "stat. only; D x A = %g; floor %g events applied at the %s"
          % (da, args.min_events, args.floor))
    print("# the same delta applies to the cos 2phi amplitude: "
          "err_azz and err_cos2phi_amplitude agree cell by cell (asserted)")
    print("# analytic = sqrt(2) / (P_zz D A sqrt(f N)), N the f = 1 "
          "accepted count summed over the cells that enter the bin")

    # the reference: f = 1, P_zz = 1, which fixes both the cell set and N
    x_c, _, n_ref, masks = combined_errors(
        args.lumi, 1.0, 1.0, 1.0, 1.0, base, g1m,
        min_events=args.min_events)
    replay = masks if args.floor == "reference" else None
    idx = [int(np.argmin(np.abs(x_c - xq))) for xq in X_QUOTED]

    print("\n# f      P_zz    x        N(f=1)      dA_zz measured   "
          "dA_zz analytic   ratio")
    worst, rows = 0.0, []
    for share in shares:
        for pzz in pzz_grid:
            _, err, _, _ = combined_errors(
                args.lumi, share, pzz, args.dilution, args.acceptance,
                base, g1m, used=replay, min_events=args.min_events)
            for i in idx:
                analytic = np.sqrt(2.0 / (share * n_ref[i])) / (pzz * da)
                ratio = err[i] / analytic
                worst = max(worst, abs(ratio - 1.0))
                rows.append((share, pzz, x_c[i], err[i], ratio))
                print("  %-6.4f %-6.2f  %-8.4f %-11.4g %-16.6g %-16.6g "
                      "%.4f" % (share, pzz, x_c[i], n_ref[i], err[i],
                                analytic, ratio))

    print("\n# the compact table for plans/07 WP2 -- delta A_zz at f = 1, "
          "x 1/sqrt(f) at any other share")
    print("  P_zz  " + "  ".join("x=%-9.4f" % x_c[i] for i in idx))
    for pzz in pzz_grid:
        _, err, _, _ = combined_errors(
            args.lumi, 1.0, pzz, args.dilution, args.acceptance, base, g1m,
            used=replay, min_events=args.min_events)
        print("  %-5.2f " % pzz
              + "  ".join("%-11.3g" % err[i] for i in idx))

    print("\n# %d cells scanned (%d shares x %d P_zz x %d x bins); "
          "max |ratio - 1| = %.3g (criterion %g)"
          % (len(rows), len(shares), len(pzz_grid), len(idx), worst,
             args.tol))
    ok = worst <= args.tol
    print("# %s: every cell reproduces sqrt(2)/(P_zz D A sqrt(f N))"
          % ("PASS" if ok else "FAIL"))
    if not ok and args.floor == "cell":
        print("#   with --floor cell a departure can be the min-events "
              "selection rather than the law: the cells the floor drops at "
              "this share are not in N.  Rerun with --floor reference to "
              "separate the two.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
