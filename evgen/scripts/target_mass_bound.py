#!/usr/bin/env python3
"""What the target-mass (gamma^2) term moved, and what residual is left.

Since 2026-08-29 every longitudinal double-spin number in this repository
carries the exact E143 form (PRD 58:112003; Anselmino-Efremov-Leader
Phys.Rept. 261:1)

  A_par = D_gamma (A1 + eta A2),  A1 = (g1 - gamma^2 g2)/F1,
  A2 = gamma (g1 + g2)/F1,        gamma^2 = 4 M^2 x^2/Q^2,
  eps = (1-y-gamma^2 y^2/4)/(1-y+y^2/2+gamma^2 y^2/4),
  D_gamma = [1-(1-y) eps]/(1+eps R),  eta = eps gamma y/[1-(1-y) eps],

by default -- `polligen.xsec.InclusiveKernel(target_mass=True)`,
`polli_fastsim.asymmetries.a_parallel_exact` -- and inverts it with the
matching divisor D_eff = D_gamma [1 - gamma^2 rho + eta gamma (1 + rho)],
rho = g2/g1 (`asymmetries.depolarization_effective`).  Before that flip
the massless A_par = D(y) g1/F1 was used on both sides, which left every
extracted g1/F1 -- and the polarized-EMC Delta-R built from it -- high by
(1 + gamma^2) + O(gamma^2 y).  This script measures BOTH halves of that
change, with the kernel's own implementation so that what is printed is
the code and not a transcription of it:

  * how far the default has moved from the retired massless path
    ("full", the A_par ratio minus one), which is the bias that is now
    gone rather than a bound on an approximation; and
  * the RESIDUAL that survives it.  The term is exact GIVEN g2, and g2 is
    Wandzura-Wilczek here, so what is left is twist-3.  The residual is
    quoted as the mis-extraction C(rho)/C(s rho) - 1 an assumed
    g2 = s g2^WW makes when the truth is g2^WW, at the two variations
    s = 0 (no g2 at all) and s = 1.5 (half again as much).

Four blocks:

  (1) the kinematic cap.  The W^2 >= 10 GeV^2 cut of `fom.Scenario`
      forces Q^2 >= (W2_min - M^2) x/(1-x), hence
      gamma^2 <= 4 M^2 x (1-x)/(W2_min - M^2) <= M^2/(W2_min - M^2),
      maximal at x = 1/2 -- printed against the measured maximum over
      the accepted 40x30 analysis grid of all six ion/energy settings.
  (2) the twelve cos 2phi sweet spots of `money_cos2phi.py` (the same
      `pick_sweet_spots_banded` call the money plot makes), per
      configuration -- and, in the same table, the other O(gamma^2)
      effect at those spots: the error the lab-angle shortcut would
      make on the cos 2phi azimuth, |phi_S(covariant) - (phi_e -
      phi_s)|, which is exactly zero for a massless target.
  (3) the polarized-EMC window of `fastsim/scripts/money_polemc.py`:
      the shift on the extracted g1A/F1A -- and so on DR -- combined
      over Q2 and the three energies with the SAME inverse-variance
      weights the published error bars use, against which the g2
      residual is the number to compare.
  (4) the tagged polarized-EMC companion of
      `tagged_polarimetry_7li.py`: the sigma-weighted shift of the
      analytic A_par overlay in its ten x bins.

Blocks (2)-(4) report the shift both as the full A_par ratio and as the
eta A2 piece alone; the ratio to A1 is printed for continuity with the
literature but is a bad figure of merit wherever g1 crosses zero (with
the grid backend at x ~ 0.03 it reads per cent on an absolute shift of
1e-4), so the number quoted in the reports is the full multiplicative
one.  At the small y of every sweet spot the whole correction collapses
to a factor (1 + gamma^2) on A_par, independent of g2: the column
full/gamma^2 shows it.

Usage:  python3 scripts/target_mass_bound.py            # toy backends
        python3 scripts/target_mass_bound.py --pdf grid # CT18/NNPDFpol
"""

import argparse
import pathlib
import sys

import numpy as np

_SCRIPTS = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS.parent))
sys.path.insert(0, str(_SCRIPTS))

from polligen import reco, tagged  # noqa: E402
from polligen.xsec import (InclusiveKernel, M_NUCLEON,  # noqa: E402
                           eta_gamma, gamma_squared)

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.asymmetries import (  # noqa: E402
    depolarization_effective)
from polli_fastsim.inputs import get_backends  # noqa: E402
from polli_fastsim.polarized import toy_b1, toy_delta_gluon  # noqa: E402
from polli_fastsim.structure import NuclearF2  # noqa: E402

W2_MIN = fom.Scenario().w2_min


# the two g2 variations the residual systematic is quoted from
G2_SCALES = (0.0, 1.5)


def kernel_pair(ion, backends, **kw):
    """(massless, finite-gamma) kernels on ONE backend pair.

    The second is the DEFAULT kernel since 2026-08-29; the first has to
    ask for `target_mass=False` explicitly."""
    common = dict(f2_source=backends["base"], g1_model=backends["g1"], **kw)
    return (InclusiveKernel(ion, target_mass=False, **common),
            InclusiveKernel(ion, **common))


def shifts(kern_tm, x, q2, y):
    """(gamma^2, full A_par ratio - 1, eta A2 / (g1/F1)) per cell."""
    t = kern_tm.tables(x, q2)
    g2v = gamma_squared(x, q2)
    f1 = np.maximum(t["f1"], 1e-30)
    a2 = np.sqrt(g2v) * (t["g1"] + t["g2"]) / f1
    eta_a2 = eta_gamma(y, g2v) * a2 / (t["g1"] / f1)
    kern_0 = InclusiveKernel(kern_tm.ion, f2_source=kern_tm.nf2.base,
                             g1_model=kern_tm.g1_model, target_mass=False)
    full = (kern_tm.a_parallel(t, x, q2, y)
            / kern_0.a_parallel(t, x, q2, y) - 1.0)
    return g2v, full, eta_a2


def g2_residual(kern_tm, x, q2, y, scales=G2_SCALES):
    """The twist-3 residual, per cell and per variation.

    An extraction assumes g2 = s g2^WW while the truth is g2^WW; the
    recovered g1/F1 is then wrong by C(rho)/C(s rho) - 1 with
    C = 1 - gamma^2 rho + eta gamma (1 + rho) and rho = g2^WW/g1.  This
    is what is left of the target-mass correction now that its
    kinematic half is computed exactly, and it is the systematic the
    reports quote in place of the removed bias.  It vanishes identically
    at s = 1 and is linear in (s - 1) up to O(gamma^4)."""
    t = kern_tm.tables(x, q2)
    g1 = t["g1"]
    rho = t["g2"] / np.where(np.abs(g1) > 1e-300, g1, 1e-300)
    truth = depolarization_effective(y, x, q2, g2_over_g1=rho,
                                     r_func=kern_tm.r_func)
    return {s: truth / depolarization_effective(y, x, q2,
                                                g2_over_g1=s * rho,
                                                r_func=kern_tm.r_func) - 1.0
            for s in scales}


def azimuth_shortcut_error(cfg, x, q2, y, n_phi=720):
    """max over the lepton azimuth of |phi_S(covariant) - (phi_e - phi_s)|
    in mrad, for a transverse spin axis: the O(gamma^2) error of the
    lab-angle shortcut that the cos 2phi observable would make if it used
    lab angles instead of `reco.azimuth_wrt_lepton_plane`.  Exactly zero
    for a massless target (`tests/test_reco.py`), so it is the same
    target-mass effect as the rest of this script, seen in the azimuth
    rather than in A_par."""
    phi_e = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    phi_s = 0.5 * np.pi
    k, p = reco.beam_fourvectors(cfg)
    s_vec = reco.spin_fourvector(phi_s)
    out = np.empty(np.shape(x))
    s_nn = cfg.sqrt_s_per_nucleon**2
    for i in range(out.size):
        kp = reco.electron_fourvector(x[i], y[i], s_nn, cfg.electron_energy,
                                      phi_e)
        err = reco.lab_azimuth_shortcut_error(k, kp, p, s_vec, phi_s)
        out[i] = 1e3 * np.abs(err).max()
    return out


def wmean(values, weights):
    w = np.sum(weights)
    return np.sum(values * weights) / w if w > 0 else np.nan


# --- (1) the kinematic cap ------------------------------------------------

def block_cap():
    cap = M_NUCLEON**2 / (W2_MIN - M_NUCLEON**2)
    print("(1) kinematic cap from W^2 >= %.4g GeV^2 (fom.Scenario):" % W2_MIN)
    print("    gamma^2 <= 4 M^2 x(1-x)/(W2min - M^2) <= M^2/(W2min - M^2)"
          " = %.5f  (at x = 1/2, M = %.4f GeV)" % (cap, M_NUCLEON))
    print("    %-26s %-28s %s"
          % ("configuration", "max gamma^2 on the grid", "with >= 100 events"))
    for ion in ("6Li", "7Li"):
        for cfg in beams.default_configs(ion):
            proj = fom.project_rates(cfg, fom.Scenario())
            g2v = gamma_squared(proj.x, proj.q2)
            acc = np.where(proj.accepted, g2v, 0.0)
            n100 = np.where(proj.accepted & (proj.n_events >= 100), g2v, 0.0)
            i, j = np.unravel_index(np.argmax(acc), acc.shape)
            k, m = np.unravel_index(np.argmax(n100), n100.shape)
            print("    %-26s %.4f at (x=%.4g, Q2=%.4g)   %.4f at (x=%.4g,"
                  " Q2=%.4g)" % (cfg.label(), acc[i, j], proj.x[i, j],
                                 proj.q2[i, j], n100[k, m], proj.x[k, m],
                                 proj.q2[k, m]))


# --- (2) the cos 2phi sweet spots ----------------------------------------

def published_spots(cfg):
    """The four money-plot-5 sweet spots of `cfg`, picked exactly as
    `money_cos2phi.main` picks them: default Scenario, toy structure
    functions and the moment_A Delta at dilution 1/3.  The selection is
    held at the published inputs whatever --pdf asks for, because these
    are the bins the report names."""
    import types
    from money_cos2phi import build_delta_model, pick_sweet_spots_banded
    args = types.SimpleNamespace(delta_model="moment_A", variant="mid_x",
                                 dilution=1.0 / 3.0, scale=1e-2)
    scenario = fom.Scenario()
    delta_func, _ = build_delta_model(args, cfg, scenario)
    proj = fom.project_rates(cfg, scenario)
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1, delta_func=delta_func)
    obs = fom.project_observables(cfg, scenario, proj, kern.g1_model,
                                  toy_b1, delta_func)
    return pick_sweet_spots_banded(proj, obs["sig_a_cos2phi"])[:4]


def block_spots(backends):
    print("\n(2) cos 2phi sweet spots (money_cos2phi.pick_sweet_spots_banded,"
          " 6Li, moment_A Delta; shifts on %s):" % backends["tag"])
    print("    %-9s %-9s %-9s %-11s %-11s %-11s %-13s %-11s %-11s %s"
          % ("x", "Q2", "y", "gamma^2", "etaA2/A1", "full", "full/gamma^2",
             "res(g2=0)", "res(1.5WW)", "dphi_S(mrad)"))
    worst = worst_res = 0.0
    for ci in range(3):
        cfg = beams.default_configs("6Li")[ci]
        s = cfg.sqrt_s_per_nucleon**2
        spots = published_spots(cfg)
        _, kern_tm = kernel_pair(cfg.ion, backends)
        x = np.array([sp[0] for sp in spots])
        q2 = np.array([sp[1] for sp in spots])
        y = q2 / (s * x)
        g2v, full, eta_a2 = shifts(kern_tm, x, q2, y)
        res = g2_residual(kern_tm, x, q2, y)
        dphi = azimuth_shortcut_error(cfg, x, q2, y)
        print("  %s" % cfg.label())
        for i in range(x.size):
            print("    %-9.4g %-9.4g %-9.4g %-11.5g %-11.4g %-11.4g %-13.3g"
                  " %-11.4g %-11.4g %.4g"
                  % (x[i], q2[i], y[i], g2v[i], eta_a2[i], full[i],
                     full[i] / g2v[i], res[0.0][i], res[1.5][i], dphi[i]))
        bad = ~np.isfinite(full)
        if np.any(bad):
            print("    (%d spot(s) not evaluable on these backends -- below"
                  " the PDF grids' Q^2 floor; here Q^2 = %s GeV^2)"
                  % (np.sum(bad),
                     ", ".join("%.3f" % v for v in np.unique(q2[bad]))))
        worst = max(worst, float(np.nanmax(np.abs(full))))
        worst_res = max(worst_res, float(np.nanmax(np.abs(res[0.0]))))
        print("    -> %s: max |full| = %.4g, max gamma^2 = %.4g,"
              " max |residual| = %.4g (g2 = 0) / %.4g (1.5 g2^WW),"
              " max dphi_S = %.3g mrad"
              % ("published configuration" if ci == 1 else "this"
                 " configuration", np.nanmax(np.abs(full)), np.max(g2v),
                 np.nanmax(np.abs(res[0.0])), np.nanmax(np.abs(res[1.5])),
                 dphi.max()))
    print("    -> over all twelve spots: max |full| = %.4g,"
          " max |residual| = %.4g" % (worst, worst_res))


# --- (3) the polarized-EMC window ----------------------------------------

def block_polemc(backends, ion="7Li", min_events=100, lumi=10.0):
    print("\n(3) polarized-EMC window (money_polemc.py combination, %s, %s,"
          " >= %d events):" % (ion, backends["tag"], min_events))
    print("    %-9s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s"
          % ("x", "<etaA2>_w", "<full>_w", "<full/g^2>_w", "worst |full|",
             "<gamma^2>_w", "<res g2=0>_w", "<res 1.5WW>_w",
             "dDR(%g/fb)" % lumi))
    rows = {}
    for cfg in beams.default_configs(ion):
        scenario = fom.Scenario(lumi_fb_per_nucleon=lumi)
        proj = fom.project_rates(
            cfg, scenario, nuclear_f2=NuclearF2(cfg.ion,
                                                base=backends["base"]))
        obs = fom.project_observables(cfg, scenario, proj, backends["g1"],
                                      toy_b1, toy_delta_gluon)
        X, Q2, Y = proj.x, proj.q2, proj.extras["y"]
        f1 = proj.extras["nf2"].f1a(X, Q2) / cfg.ion.A
        g1n = backends["g1"].g1_nucleus(cfg.ion, X, Q2) / cfg.ion.A
        a1_naive = np.abs(g1n / np.maximum(f1, 1e-30))
        err_dr = obs["err_g1_over_f1"] / np.maximum(a1_naive, 1e-12)
        use = (proj.accepted & (proj.n_events >= min_events) & (err_dr > 0)
               & np.isfinite(err_dr))
        _, kern_tm = kernel_pair(cfg.ion, backends)
        g2v, full, eta_a2 = shifts(kern_tm, X, Q2, Y)
        res = g2_residual(kern_tm, X, Q2, Y)
        good = (np.isfinite(full) & np.isfinite(eta_a2)
                & np.isfinite(res[0.0]) & np.isfinite(res[1.5]))
        if not np.all(good[use]):
            print("    (%s: %d of %d accepted cells dropped, backend"
                  " not evaluable there)"
                  % (cfg.label(), np.sum(use & ~good), np.sum(use)))
        use = use & good
        full = np.where(good, full, 0.0)
        eta_a2 = np.where(good, eta_a2, 0.0)
        r0 = np.where(good, res[0.0], 0.0)
        r15 = np.where(good, res[1.5], 0.0)
        w = np.zeros_like(err_dr)
        np.divide(1.0, err_dr**2, out=w, where=use)
        for i in range(X.shape[0]):
            if not use[i].any():
                continue
            acc = rows.setdefault(X[i, 0], [0.0] * 8)
            acc[0] += np.sum(w[i] * eta_a2[i])
            acc[1] += np.sum(w[i] * full[i])
            acc[2] += np.sum(w[i] * np.where(use[i], full[i]
                                             / np.maximum(g2v[i], 1e-300),
                                             0.0))
            acc[3] += np.sum(w[i] * g2v[i])
            acc[4] += np.sum(w[i])
            here = np.abs(np.where(use[i], full[i], 0.0))
            acc[5] = max(acc[5], np.max(here))
            acc[6] += np.sum(w[i] * r0[i])
            acc[7] += np.sum(w[i] * r15[i])
    xs = np.array(sorted(rows))
    for xq in (0.1, 0.3, 0.5, 0.7):
        xi = xs[np.argmin(np.abs(xs - xq))]
        a = rows[xi]
        # the published statistical error on Delta R is 1/sqrt(sum w)
        # over exactly these cells (money_polemc.delta_dr_per_x)
        print("    %-9.4f %-13.4g %-13.4g %-13.4g %-13.4g %-13.4g %-13.4g"
              " %-13.4g %.4g"
              % (xi, a[0] / a[4], a[1] / a[4], a[2] / a[4], a[5],
                 a[3] / a[4], a[6] / a[4], a[7] / a[4],
                 1.0 / np.sqrt(a[4])))


# --- (4) the tagged polarized-EMC companion ------------------------------

def block_tagged(backends):
    print("\n(4) tagged A_par overlay (tagged_polarimetry_7li.py, quasi-free"
          " triton, %s, sigma-weighted):" % backends["tag"])
    x_edges = np.logspace(-3, np.log10(0.7), 11)
    model = tagged.TaggedModel(tagged.li7_alpha_channel())
    for ci in range(3):
        cfg = beams.default_configs("7Li")[ci]
        _, kern_tm = kernel_pair(tagged.TRITON, backends)
        smp = tagged.TaggedSampler(model, kern_tm, cfg, fom.Scenario())
        inner = smp.inner
        x, q2, w = inner.x_cells, inner.q2_cells, inner.xsec_flat
        y = q2 / (inner.s * x)
        g2v, full, eta_a2 = shifts(kern_tm, x, q2, y)
        res = g2_residual(kern_tm, x, q2, y)
        good = (np.isfinite(full) & np.isfinite(eta_a2)
                & np.isfinite(res[0.0]) & np.isfinite(res[1.5]))
        if not np.all(good):
            print("    (%d of %d cells dropped, backend not evaluable"
                  " there)" % (np.sum(~good), good.size))
        x, q2, w = x[good], q2[good], w[good]
        g2v, full, eta_a2 = g2v[good], full[good], eta_a2[good]
        r0, r15 = res[0.0][good], res[1.5][good]
        xb = np.digitize(x, x_edges) - 1
        best_full = best_eta = best_res = 0.0
        rows = []
        for b in range(x_edges.size - 1):
            m = xb == b
            if not m.any():
                continue
            rows.append((np.sqrt(x_edges[b] * x_edges[b + 1]),
                         wmean(g2v[m], w[m]), wmean(eta_a2[m], w[m]),
                         wmean(full[m], w[m]),
                         full[m][np.argmax(np.abs(full[m]))],
                         wmean(r0[m], w[m]), wmean(r15[m], w[m])))
            best_full = max(best_full, abs(rows[-1][3]))
            best_eta = max(best_eta, abs(rows[-1][2]))
            best_res = max(best_res, abs(rows[-1][5]))
        print("  %s%s" % (cfg.label(),
                          "   <- published (--config 1)" if ci == 1 else ""))
        print("    %-9s %-12s %-12s %-12s %-12s %-12s %s"
              % ("x_bin", "<gamma^2>", "<etaA2/A1>", "<full>", "worst cell",
                 "<res g2=0>", "<res 1.5WW>"))
        for r in rows:
            print("    %-9.4g %-12.4g %-12.4g %-12.4g %-12.4g %-12.4g %.4g"
                  % r)
        print("    -> max over x bins: |<full>| = %.4g, |<etaA2/A1>| = %.4g,"
              " |<residual>| = %.4g (g2 = 0)" % (best_full, best_eta,
                                                 best_res))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pdf", default="toy", choices=("toy", "grid"),
                    help="structure-function backends; toy is what every "
                         "published number uses")
    ap.add_argument("--blocks", default="1234",
                    help="which blocks to run, e.g. '13'")
    args = ap.parse_args()
    backends = get_backends(args.pdf)
    print("target-mass (gamma^2) term of the inclusive A_par chain, ON by "
          "default since 2026-08-29 -- backends: %s, M = %.4f GeV, "
          "g2 = g2_WW, residual quoted at g2 = %s x g2_WW"
          % (backends["tag"], M_NUCLEON,
             " and ".join("%g" % v for v in G2_SCALES)))
    if "1" in args.blocks:
        block_cap()
    if "2" in args.blocks:
        block_spots(backends)
    if "3" in args.blocks:
        block_polemc(backends)
    if "4" in args.blocks:
        block_tagged(backends)


if __name__ == "__main__":
    main()
