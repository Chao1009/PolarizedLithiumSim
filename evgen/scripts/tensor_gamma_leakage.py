#!/usr/bin/env python3
"""How much of the tensor (b1-b4) sector leaks into the cos 2phi amplitude.

The gluonometry observable of this programme is the cos 2phi modulation of
a TRANSVERSELY tensor-polarized fill, and the structure function it
extracts is Delta.  In the massless master formula the b1/b2 rate sector
is purely phi-independent and cannot contaminate it.  At finite gamma it
can: the virtual photon is not along the beam (Cosyn Eq. 24 gives the
rest-frame angle theta_q between them), so a spin axis transverse to the
BEAM is not transverse to q, and the rank-2 alignment tensor acquires
photon-frame components that depend on the lepton-plane azimuth.

Until 2026-08-29 the size of that leakage was quoted as a bound,
Delta_fake = 1.15 gamma^2 b1 -- the Eq. (17e) term alone, gamma^2 b1/6,
times the 6.9 the full combination was estimated to be.  This script
replaces the bound with the measurement, at the same twelve sweet spots
`target_mass_bound.py` uses (the four money-plot-5 spots of each of the
three 6Li configurations, picked by `money_cos2phi.pick_sweet_spots_banded`
exactly as the money plot picks them), with the kernel's own exact
finite-gamma tensor path (`InclusiveKernel(tensor_gamma=True)`,
plans/08 D2) so that what is printed is the code and not a transcription.

Three contributions make up the leakage, all O(gamma^2), and the two
large ones very nearly cancel:

  T_LL  the b1/b2 RATE term seen through the tilted axis (sin^2 theta_q
        times the leading-twist combination) -- +3x the Eq. (17e) term;
  T_LT  Eq. (17d), higher twist but O(gamma) in the structure function and
        O(gamma) in the geometry -- -3x, opposite to the T_LL one;
  T_TT  Eq. (17e) itself, the only one the earlier bound counted.

The 3 : -3 : 1 ratio is exact as gamma^2 and y go to zero, so what
survives at the sweet spots is the twist-4 Eq. (17e) term almost alone
and the coefficient Delta_fake/(gamma^2 b1) sits at 0.14-0.16, i.e. at
the 1/6 of Eq. (17e).  The guess the 1.15 came from had the T_LT term
ADDING; the relative sign is fixed by both finite-gamma rows of the
paper's own Table 1, which the kernel reproduces exactly (see
`tests/test_tensor_gamma.py`).  Because the residual is the small
difference of two large terms, it is the higher-twist b3/b4 and not
gamma^2 that now dominates its uncertainty -- see below.

Columns: gamma^2; the three contributions' sum as an equivalent Delta,
Delta_fake/F1 and the coefficient Delta_fake/(gamma^2 b1) that REPLACES
the 1.15; the leakage as a fraction of the published Delta amplitude at
that spot (moment_A at dilution 1/3, the model the money plot draws); and
the ratio of the full combination to the Eq. (17e) term alone.

Sign: with the literature tensor convention adopted 2026-08-29
(`TENSOR_LL_SIGN = -1`, plans/08 D1) and a positive b1, the leakage is
NEGATIVE at every spot, i.e. OPPOSITE in sign to the cos 2phi amplitude
of the moment-constrained (negative-Delta) models, so it cancels part of
the measured amplitude instead of faking one.  Flipping the constant back
flips it, which is why D2 was gated on D1.  It is subtractable in the
same data: it is proportional to b1, which the A_zz of the longitudinal
fills measures at the same (x, Q^2).

b3 and b4 are the unmeasured higher-twist tensor structure functions;
they default to zero, and `--b3-frac`/`--b4-frac` set them to a fraction
of b2 to exercise the slots.  Since they break the T_LL-T_LT
cancellation, they matter at the level of the whole effect rather than
at ten per cent of it: b3 = b4 = 0.1 b2 moves the coefficient below from
0.14-0.16 to 0.23-0.26 and the worst fraction from 0.11% to 0.17%.  That
is the honest width of the systematic, and it is still four times smaller
than the bound it replaces.

Usage:  python3 scripts/tensor_gamma_leakage.py
        python3 scripts/tensor_gamma_leakage.py --b3-frac 0.1 --b4-frac 0.1
"""

import argparse
import pathlib
import sys

import numpy as np

_SCRIPTS = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS.parent))
sys.path.insert(0, str(_SCRIPTS))

from polligen.xsec import (EventSpinState, InclusiveKernel,  # noqa: E402
                           cosyn_tensor_sfs, cosyn_unpolarized_sfs,
                           gamma_squared, theta_q_cos_sin)

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.asymmetries import epsilon_gamma  # noqa: E402
from polli_fastsim.polarized import toy_b1  # noqa: E402

from target_mass_bound import published_spots  # noqa: E402


def delta_model_of(cfg):
    """The published Delta model of the money plot (moment_A, dilution 1/3)."""
    import types
    from money_cos2phi import build_delta_model
    args = types.SimpleNamespace(delta_model="moment_A", variant="mid_x",
                                 dilution=1.0 / 3.0, scale=1e-2)
    return build_delta_model(args, cfg, fom.Scenario())[0]


def contributions(t, x, q2, y, b3f, b4f):
    """(T_LL, T_LT, T_TT) pieces of the cos 2phi harmonic at theta_S = 90 deg,
    for the m = +1 state -- the same decomposition
    `InclusiveKernel._tensor_harmonics_gamma` sums, split by channel."""
    g2 = gamma_squared(x, q2)
    eps = epsilon_gamma(y, g2)
    c, sn = theta_q_cos_sin(y, g2)
    b1, b2 = t["b1"], t["b2"]
    b3, b4 = b3f * b2, b4f * b2
    f_t, f_l, f_lt, f_tt = cosyn_tensor_sfs(b1, b2, b3, b4, x, g2)
    fu_t, fu_l = cosyn_unpolarized_sfs(t["f1"], t["f2"], x, g2)
    den = fu_t + eps * fu_l
    # cos 2phi' coefficients of t_zz, t_xz, t_xx - t_yy at theta_S = 90 deg,
    # with the pure m = +1 prefactor 3 Q_NN/2 = 1/2
    half = 0.5
    return (half * (0.5 * sn * sn) * (f_t + eps * f_l) / den,
            -half * (0.5 * c * sn) * np.sqrt(2.0 * eps * (1.0 + eps)) * f_lt
            / den,
            half * (0.5 * (c * c + 1.0)) * eps * f_tt / den)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--b3-frac", type=float, default=0.0,
                    help="b3 as a fraction of b2 (default 0)")
    ap.add_argument("--b4-frac", type=float, default=0.0,
                    help="b4 as a fraction of b2 (default 0)")
    args = ap.parse_args()

    print("cos 2phi leakage of the tensor sector at the twelve money-plot-5"
          " sweet spots")
    print("  kernel: InclusiveKernel(tensor_gamma=True), toy b1, b3 = %.3g b2,"
          " b4 = %.3g b2" % (args.b3_frac, args.b4_frac))
    print("  %-8s %-8s %-9s %-10s %-9s %-11s %-11s %-9s %-9s %s"
          % ("x", "Q2", "y", "gamma^2", "b1/F1", "a2(tensor)", "Dfake/F1",
             "Dfake/", "a2tens/", "full/"))
    print("  %-8s %-8s %-9s %-10s %-9s %-11s %-11s %-9s %-9s %s"
          % ("", "", "", "", "", "(m=+1)", "", "(g^2 b1)", "a2(Delta)",
             "17e"))
    worst_frac = worst_coef = 0.0
    for ci in range(3):
        cfg = beams.default_configs("6Li")[ci]
        s = cfg.sqrt_s_per_nucleon ** 2
        delta_func = delta_model_of(cfg)
        spots = published_spots(cfg)
        x = np.array([sp[0] for sp in spots])
        q2 = np.array([sp[1] for sp in spots])
        y = q2 / (s * x)
        b3 = ((lambda xx, qq, f1: args.b3_frac * 2.0 * xx
               * toy_b1(xx, qq, f1)) if args.b3_frac else None)
        b4 = ((lambda xx, qq, f1: args.b4_frac * 2.0 * xx
               * toy_b1(xx, qq, f1)) if args.b4_frac else None)
        kern = InclusiveKernel(beams.LI6, b1_func=toy_b1,
                               delta_func=delta_func, b3_func=b3, b4_func=b4,
                               tensor_gamma=True)
        kern0 = InclusiveKernel(beams.LI6, b1_func=toy_b1,
                                delta_func=delta_func)
        t = kern.tables(x, q2)
        state = EventSpinState(0, 0.0, 1.0, 1.0, theta_s=0.5 * np.pi)
        # tensor-only kernel: the Delta amplitude is added on top of it
        kern_nod = InclusiveKernel(beams.LI6, b1_func=toy_b1, b3_func=b3,
                                   b4_func=b4, tensor_gamma=True)
        a2_tensor = kern_nod.amplitudes(kern_nod.tables(x, q2), x, q2, s,
                                        state)[2]
        a2_delta = kern0.amplitudes(kern0.tables(x, q2), x, q2, s, state)[2]
        g2v = gamma_squared(x, q2)
        dphi = t["f1"] + (1.0 - y) / (x * y * y) * t["f2"]
        d_fake = -a2_tensor * y * y * dphi / (1.0 - y)
        ll, lt, tt = contributions(t, x, q2, y, args.b3_frac, args.b4_frac)
        print("  %s" % cfg.label())
        for i in range(x.size):
            print("  %-8.4g %-8.4g %-9.4g %-10.5g %-9.4g %-11.4g %-11.4g"
                  " %-9.4g %-9.4g %.4g"
                  % (x[i], q2[i], y[i], g2v[i], t["b1"][i] / t["f1"][i],
                     a2_tensor[i], d_fake[i] / t["f1"][i],
                     d_fake[i] / (g2v[i] * t["b1"][i]),
                     a2_tensor[i] / a2_delta[i],
                     (ll[i] + lt[i] + tt[i]) / tt[i]))
            worst_frac = max(worst_frac, abs(a2_tensor[i] / a2_delta[i]))
            worst_coef = max(worst_coef, abs(d_fake[i]
                                             / (g2v[i] * t["b1"][i])))
        print("    channels (T_LL : T_LT : T_TT) of the cos 2phi harmonic: "
              + ", ".join("%.3f:%.3f:%.3f"
                          % (ll[i] / tt[i], lt[i] / tt[i], 1.0)
                          for i in range(x.size)))
    print("\n  worst |a2(tensor)/a2(Delta)| = %.4g (%.3g%%);"
          " worst Delta_fake/(gamma^2 b1) = %.4g"
          % (worst_frac, 100.0 * worst_frac, worst_coef))


if __name__ == "__main__":
    main()
