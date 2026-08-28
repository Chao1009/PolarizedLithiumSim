#!/usr/bin/env python3
"""Diagnostic: sig^2 and L_5sigma grid for 6Li gluonometry cos(2phi).

Verifies whether the mid EIC config (E_e = 10 GeV, p_ion = 99.5 GeV/u)
really gives the lowest L_5sigma on the Delta/F1 = 1e-3 Sather-Schmidt
case, or whether that is a plot-reading artifact.  Prints a table; no
plots.  The scanned grid deliberately runs off the machine menu -- the
question is which sqrt(s) the figure of merit likes, and the three
canonical configurations are printed underneath it for comparison.
"""

import importlib.util
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams
from polli_fastsim.beams import LI6, BeamConfig
from polli_fastsim.fom import Scenario, project_rates
from polli_fastsim.structure import NuclearF2
from polli_fastsim.inputs import get_backends


# ---------------------------------------------------------------------------
# The significance is money_delta.py's own, IMPORTED rather than copied: it
# used to be a verbatim copy under a "do not modify" comment, which is a
# convention and not a mechanism -- and when the min-events floor moved to
# the luminosity the reach is quoted at (2026-08-28) the copy would have
# been left behind.  `reach_fb` replaces the old `25 / sig2_per_fb_at`.
# ---------------------------------------------------------------------------
def _load_money_delta():
    path = pathlib.Path(__file__).resolve().parent / "money_delta.py"
    spec = importlib.util.spec_from_file_location("_money_delta", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_md = _load_money_delta()
bin_terms = _md.bin_terms
reach_from_terms = _md.reach_from_terms
sig2_per_fb_at = _md.sig2_per_fb_at
reach_fb = _md.reach_fb


# ---------------------------------------------------------------------------
# Helper: total accepted events (rate per fb^-1/u, over the bins the
# significance beside it actually sums)
# ---------------------------------------------------------------------------
def n_accepted_at(cfg, base=None, min_events=10, lumi_fb=1.0):
    """Total accepted event RATE, per fb^-1/nucleon, over the bins that
    hold `min_events` at `lumi_fb`.

    The rate is per fb^-1/u whatever `lumi_fb` is; only the min-events
    mask moves with it.  Passing the solved L_5sigma keeps this column and
    the sig^2 column beside it describing the same set of bins -- until
    2026-08-28 the mask sat at 1 fb^-1/u here while the significance had
    already moved to the reach, which at the canonical configurations
    (16.7-21.8 fb^-1/u) differs by up to three bins."""
    sc = Scenario(lumi_fb_per_nucleon=1.0)
    nf2_in = NuclearF2(cfg.ion, base=base) if base is not None else None
    proj = project_rates(cfg, sc, nuclear_f2=nf2_in)
    use = proj.accepted & (proj.n_events * lumi_fb >= min_events)
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
    p_ion_vals = [20.0, 40.8, 60.0, 80.0, 99.5, 137.5]  # GeV/u

    # Canonical EIC configs for 6Li -- DERIVED, never hard-coded (plans/10:
    # the two lower ones are gamma-matched, not rigidity-scaled)
    canonical = [(c.electron_energy, c.ion_momentum_per_nucleon)
                 for c in beams.default_configs("6Li")]

    # -----------------------------------------------------------------------
    # Full grid
    # -----------------------------------------------------------------------
    rows = []
    for Ee in E_e_vals:
        for p in p_ion_vals:
            cfg = BeamConfig(electron_energy=Ee, ion=LI6,
                             ion_momentum_per_nucleon=p)
            sqrts = np.sqrt(4.0 * Ee * p)
            terms, n_ev = bin_terms(cfg, scale=SCALE, pzz=P_ZZ, base=base)
            l5sig = reach_from_terms(terms, n_ev, min_events=MIN_EVENTS)
            sig2 = float(terms[n_ev * l5sig >= MIN_EVENTS].sum())
            n_acc = n_accepted_at(cfg, base=base, min_events=MIN_EVENTS,
                                  lumi_fb=l5sig)
            rows.append((Ee, p, sqrts, n_acc, sig2, l5sig))

    # -----------------------------------------------------------------------
    # Print full grid table
    # -----------------------------------------------------------------------
    hdr = (f"{'E_e':>5}  {'p_ion':>7}  {'sqrt(s)':>8}  "
           f"{'N_acc/fb^-1/u @L5s':>18}  {'sig^2':>12}  "
           f"{'L_5sig (fb^-1/u)':>18}")
    sep = "-" * len(hdr)
    print()
    print("=" * len(hdr))
    print("  6Li gluonometry cos(2phi) grid diagnostic")
    print(f"  Delta/F1 scale = {SCALE:.0e}  |  P_zz = {P_ZZ}  |  backend = toy")
    print(f"  N_acc and sig^2 are both summed over the bins holding "
          f"{MIN_EVENTS} events AT L_5sig, not at 1 fb^-1/u")
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
        terms, n_ev = bin_terms(cfg, scale=SCALE, pzz=P_ZZ, base=base)
        l5sig = reach_from_terms(terms, n_ev, min_events=MIN_EVENTS)
        sig2 = float(terms[n_ev * l5sig >= MIN_EVENTS].sum())
        n_acc = n_accepted_at(cfg, base=base, min_events=MIN_EVENTS,
                              lumi_fb=l5sig)
        label = f"e({Ee_c:g})x6Li({p_c:g}/u)"
        print(f"{Ee_c:>5.1f}  {p_c:>7.1f}  {sqrts:>8.2f}  "
              f"{n_acc:>18.3e}  {sig2:>12.4e}  {l5sig:>18.3f}"
              f"   <- {label}")
    print(sep)
    print()


if __name__ == "__main__":
    main()
