#!/usr/bin/env python3
"""Diagnostic: sig^2 and L_5sigma grid for 6Li gluonometry cos(2phi).

Verifies whether the mid EIC config (E_e=10 GeV, p_ion=50 GeV/u) really
gives the lowest L_5sigma on the Delta/F1=1e-3 Sather-Schmidt case, or
whether that is a plot-reading artifact.  Prints a table; no plots.
"""

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim.beams import LI6, BeamConfig
from polli_fastsim.fom import Scenario, project_rates
from polli_fastsim.asymmetries import a_cos2phi
from polli_fastsim.polarized import toy_delta_gluon
from polli_fastsim.structure import NuclearF2
from polli_fastsim.inputs import get_backends


# ---------------------------------------------------------------------------
# Exact copy of sig2_per_fb_at from money_delta.py — do not modify.
# ---------------------------------------------------------------------------
def sig2_per_fb_at(cfg, scale, pzz, base=None, min_events=10):
    """Significance^2 per fb^-1/nucleon at a reference Delta/F1 scale.
    The toy Delta shape is linear in `scale`, so reach curves follow
    L_5sig(s) = 25 / (sig2 * (s/scale)^2) analytically."""
    sc = Scenario(lumi_fb_per_nucleon=1.0, pol_ion_tensor=pzz)
    nf2_in = NuclearF2(cfg.ion, base=base) if base is not None else None
    proj = project_rates(cfg, sc, nuclear_f2=nf2_in)
    nf2 = proj.extras["nf2"]
    f1 = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2 = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y = proj.extras["y"]
    delta = toy_delta_gluon(proj.x, proj.q2, f1, scale=scale)
    amp = a_cos2phi(delta, f1, f2, proj.x, y)
    use = proj.accepted & (proj.n_events >= min_events)
    return np.where(use, amp**2 * pzz**2 * proj.n_events / 2.0, 0.0).sum()


# ---------------------------------------------------------------------------
# Helper: total accepted events (at 1 fb^-1/u)
# ---------------------------------------------------------------------------
def n_accepted_at(cfg, base=None, min_events=10):
    """Total accepted event count at 1 fb^-1/nucleon."""
    sc = Scenario(lumi_fb_per_nucleon=1.0)
    nf2_in = NuclearF2(cfg.ion, base=base) if base is not None else None
    proj = project_rates(cfg, sc, nuclear_f2=nf2_in)
    use = proj.accepted & (proj.n_events >= min_events)
    return proj.n_events[use].sum()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    P_ZZ = 0.8
    SCALE = 1e-3
    MIN_EVENTS = 10

    backends = get_backends("toy")
    base = backends["base"]

    E_e_vals   = [5.0, 7.0, 10.0, 14.0, 18.0]          # GeV
    p_ion_vals = [20.0, 35.0, 50.0, 75.0, 100.0, 137.5]  # GeV/u

    # Canonical EIC configs for 6Li (from beams.default_configs)
    canonical = [(5.0, 20.5), (10.0, 50.0), (18.0, 137.5)]

    # -----------------------------------------------------------------------
    # Full grid
    # -----------------------------------------------------------------------
    rows = []
    for Ee in E_e_vals:
        for p in p_ion_vals:
            cfg = BeamConfig(electron_energy=Ee, ion=LI6,
                             ion_momentum_per_nucleon=p)
            sqrts = np.sqrt(4.0 * Ee * p)
            sig2 = sig2_per_fb_at(cfg, scale=SCALE, pzz=P_ZZ, base=base,
                                  min_events=MIN_EVENTS)
            l5sig = 25.0 / max(sig2, 1e-30)
            n_acc = n_accepted_at(cfg, base=base, min_events=MIN_EVENTS)
            rows.append((Ee, p, sqrts, n_acc, sig2, l5sig))

    # -----------------------------------------------------------------------
    # Print full grid table
    # -----------------------------------------------------------------------
    hdr = (f"{'E_e':>5}  {'p_ion':>7}  {'sqrt(s)':>8}  "
           f"{'N_acc (fb^-1/u)':>18}  {'sig^2':>12}  {'L_5sig (fb^-1/u)':>18}")
    sep = "-" * len(hdr)
    print()
    print("=" * len(hdr))
    print("  6Li gluonometry cos(2phi) grid diagnostic")
    print(f"  Delta/F1 scale = {SCALE:.0e}  |  P_zz = {P_ZZ}  |  backend = toy")
    print("=" * len(hdr))
    print(hdr)
    print(sep)
    for (Ee, p, sqrts, n_acc, sig2, l5sig) in rows:
        print(f"{Ee:>5.1f}  {p:>7.1f}  {sqrts:>8.2f}  "
              f"{n_acc:>18.3e}  {sig2:>12.4e}  {l5sig:>18.3f}")
    print(sep)

    # -----------------------------------------------------------------------
    # Best config
    # -----------------------------------------------------------------------
    best = min(rows, key=lambda r: r[5])
    print()
    print(">>> BEST CONFIG (minimum L_5sig on this grid):")
    print(f"    E_e = {best[0]:.1f} GeV,  p_ion = {best[1]:.1f} GeV/u,  "
          f"sqrt(s) = {best[2]:.2f} GeV")
    print(f"    sig^2 = {best[4]:.4e},  L_5sig = {best[5]:.3f} fb^-1/u")

    # -----------------------------------------------------------------------
    # Canonical configs subsection
    # -----------------------------------------------------------------------
    print()
    print("=" * len(hdr))
    print("  Canonical EIC configs for 6Li")
    print("=" * len(hdr))
    print(hdr)
    print(sep)
    for (Ee_c, p_c) in canonical:
        cfg = BeamConfig(electron_energy=Ee_c, ion=LI6,
                         ion_momentum_per_nucleon=p_c)
        sqrts = np.sqrt(4.0 * Ee_c * p_c)
        sig2 = sig2_per_fb_at(cfg, scale=SCALE, pzz=P_ZZ, base=base,
                              min_events=MIN_EVENTS)
        l5sig = 25.0 / max(sig2, 1e-30)
        n_acc = n_accepted_at(cfg, base=base, min_events=MIN_EVENTS)
        label = f"e({Ee_c:g})x6Li({p_c:g}/u)"
        print(f"{Ee_c:>5.1f}  {p_c:>7.1f}  {sqrts:>8.2f}  "
              f"{n_acc:>18.3e}  {sig2:>12.4e}  {l5sig:>18.3f}"
              f"   <- {label}")
    print(sep)
    print()


if __name__ == "__main__":
    main()
