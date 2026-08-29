#!/usr/bin/env python3
"""How big is the target-mass (gamma^2) term the A_par chain drops?

Every published longitudinal double-spin number in this repository is
built from the massless approximation A_par = D(y) g1/F1
(`polli_fastsim.asymmetries.a_parallel`, `polligen.xsec.InclusiveKernel`
with target_mass=False).  The exact E143 form (PRD 58:112003;
Anselmino-Efremov-Leader Phys.Rept. 261:1) is

  A_par = D_gamma (A1 + eta A2),  A1 = (g1 - gamma^2 g2)/F1,
  A2 = gamma (g1 + g2)/F1,        gamma^2 = 4 M^2 x^2/Q^2,
  eps = (1-y-gamma^2 y^2/4)/(1-y+y^2/2+gamma^2 y^2/4),
  D_gamma = [1-(1-y) eps]/(1+eps R),  eta = eps gamma y/[1-(1-y) eps],

and this script measures the difference where it could matter, with the
kernel's own implementation (`InclusiveKernel(..., target_mass=True)`)
and g2 = g2_WW, so that what is printed is the code, not a transcription
of it.  Four blocks:

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
      weights the published error bars use.
  (4) the tagged polarized-EMC companion of
      `tagged_polarimetry_7li.py`: the sigma-weighted shift of the
      analytic D g1t/F1t overlay in its ten x bins.

Blocks (2)-(4) report the shift both as the full A_par ratio and as the
eta A2 piece alone; the ratio to A1 is printed for continuity with the
literature but is a bad figure of merit wherever g1 crosses zero (with
the grid backend at x ~ 0.03 it reads per cent on an absolute shift of
1e-4), so the bound quoted in the reports is the full multiplicative
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
from polli_fastsim.inputs import get_backends  # noqa: E402
from polli_fastsim.polarized import toy_b1, toy_delta_gluon  # noqa: E402
from polli_fastsim.structure import NuclearF2  # noqa: E402

W2_MIN = fom.Scenario().w2_min


def kernel_pair(ion, backends, **kw):
    """(massless, finite-gamma) kernels on ONE backend pair."""
    common = dict(f2_source=backends["base"], g1_model=backends["g1"], **kw)
    return (InclusiveKernel(ion, **common),
            InclusiveKernel(ion, target_mass=True, **common))


def shifts(kern_tm, x, q2, y):
    """(gamma^2, full A_par ratio - 1, eta A2 / (g1/F1)) per cell."""
    t = kern_tm.tables(x, q2)
    g2v = gamma_squared(x, q2)
    f1 = np.maximum(t["f1"], 1e-30)
    a2 = np.sqrt(g2v) * (t["g1"] + t["g2"]) / f1
    eta_a2 = eta_gamma(y, g2v) * a2 / (t["g1"] / f1)
    kern_0 = InclusiveKernel(kern_tm.ion, f2_source=kern_tm.nf2.base,
                             g1_model=kern_tm.g1_model)
    full = (kern_tm.a_parallel(t, x, q2, y)
            / kern_0.a_parallel(t, x, q2, y) - 1.0)
    return g2v, full, eta_a2


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
    print("    %-9s %-9s %-9s %-11s %-11s %-11s %-13s %s"
          % ("x", "Q2", "y", "gamma^2", "etaA2/A1", "full", "full/gamma^2",
             "dphi_S(mrad)"))
    worst = 0.0
    for ci in range(3):
        cfg = beams.default_configs("6Li")[ci]
        s = cfg.sqrt_s_per_nucleon**2
        spots = published_spots(cfg)
        _, kern_tm = kernel_pair(cfg.ion, backends)
        x = np.array([sp[0] for sp in spots])
        q2 = np.array([sp[1] for sp in spots])
        y = q2 / (s * x)
        g2v, full, eta_a2 = shifts(kern_tm, x, q2, y)
        dphi = azimuth_shortcut_error(cfg, x, q2, y)
        print("  %s" % cfg.label())
        for i in range(x.size):
            print("    %-9.4g %-9.4g %-9.4g %-11.5g %-11.4g %-11.4g %-13.3g"
                  " %.4g"
                  % (x[i], q2[i], y[i], g2v[i], eta_a2[i], full[i],
                     full[i] / g2v[i], dphi[i]))
        bad = ~np.isfinite(full)
        if np.any(bad):
            print("    (%d spot(s) not evaluable on these backends -- below"
                  " the PDF grids' Q^2 floor; here Q^2 = %s GeV^2)"
                  % (np.sum(bad),
                     ", ".join("%.3f" % v for v in np.unique(q2[bad]))))
        worst = max(worst, float(np.nanmax(np.abs(full))))
        print("    -> %s: max |full| = %.4g, max gamma^2 = %.4g,"
              " max dphi_S = %.3g mrad"
              % ("published configuration" if ci == 1 else "this"
                 " configuration", np.nanmax(np.abs(full)), np.max(g2v),
                 dphi.max()))
    print("    -> over all twelve spots: max |full| = %.4g" % worst)


# --- (3) the polarized-EMC window ----------------------------------------

def block_polemc(backends, ion="7Li", min_events=100, lumi=10.0):
    print("\n(3) polarized-EMC window (money_polemc.py combination, %s, %s,"
          " >= %d events):" % (ion, backends["tag"], min_events))
    print("    %-9s %-13s %-13s %-13s %-13s %-13s %s"
          % ("x", "<etaA2>_w", "<full>_w", "<full/g^2>_w", "worst |full|",
             "<gamma^2>_w", "dDR(%g/fb)" % lumi))
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
        good = np.isfinite(full) & np.isfinite(eta_a2)
        if not np.all(good[use]):
            print("    (%s: %d of %d accepted cells dropped, backend"
                  " not evaluable there)"
                  % (cfg.label(), np.sum(use & ~good), np.sum(use)))
        use = use & good
        full = np.where(good, full, 0.0)
        eta_a2 = np.where(good, eta_a2, 0.0)
        w = np.zeros_like(err_dr)
        np.divide(1.0, err_dr**2, out=w, where=use)
        for i in range(X.shape[0]):
            if not use[i].any():
                continue
            acc = rows.setdefault(X[i, 0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            acc[0] += np.sum(w[i] * eta_a2[i])
            acc[1] += np.sum(w[i] * full[i])
            acc[2] += np.sum(w[i] * np.where(use[i], full[i]
                                             / np.maximum(g2v[i], 1e-300),
                                             0.0))
            acc[3] += np.sum(w[i] * g2v[i])
            acc[4] += np.sum(w[i])
            here = np.abs(np.where(use[i], full[i], 0.0))
            acc[5] = max(acc[5], np.max(here))
    xs = np.array(sorted(rows))
    for xq in (0.1, 0.3, 0.5, 0.7):
        xi = xs[np.argmin(np.abs(xs - xq))]
        a = rows[xi]
        # the published statistical error on Delta R is 1/sqrt(sum w)
        # over exactly these cells (money_polemc.delta_dr_per_x)
        print("    %-9.4f %-13.4g %-13.4g %-13.4g %-13.4g %-13.4g %.4g"
              % (xi, a[0] / a[4], a[1] / a[4], a[2] / a[4], a[5],
                 a[3] / a[4], 1.0 / np.sqrt(a[4])))


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
        good = np.isfinite(full) & np.isfinite(eta_a2)
        if not np.all(good):
            print("    (%d of %d cells dropped, backend not evaluable"
                  " there)" % (np.sum(~good), good.size))
        x, q2, w = x[good], q2[good], w[good]
        g2v, full, eta_a2 = g2v[good], full[good], eta_a2[good]
        xb = np.digitize(x, x_edges) - 1
        best_full = best_eta = 0.0
        rows = []
        for b in range(x_edges.size - 1):
            m = xb == b
            if not m.any():
                continue
            rows.append((np.sqrt(x_edges[b] * x_edges[b + 1]),
                         wmean(g2v[m], w[m]), wmean(eta_a2[m], w[m]),
                         wmean(full[m], w[m]),
                         full[m][np.argmax(np.abs(full[m]))]))
            best_full = max(best_full, abs(rows[-1][3]))
            best_eta = max(best_eta, abs(rows[-1][2]))
        print("  %s%s" % (cfg.label(),
                          "   <- published (--config 1)" if ci == 1 else ""))
        print("    %-9s %-12s %-12s %-12s %s"
              % ("x_bin", "<gamma^2>", "<etaA2/A1>", "<full>", "worst cell"))
        for r in rows:
            print("    %-9.4g %-12.4g %-12.4g %-12.4g %.4g" % r)
        print("    -> max over x bins: |<full>| = %.4g, |<etaA2/A1>| = %.4g"
              % (best_full, best_eta))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pdf", default="toy", choices=("toy", "grid"),
                    help="structure-function backends; toy is what every "
                         "published number uses")
    ap.add_argument("--blocks", default="1234",
                    help="which blocks to run, e.g. '13'")
    args = ap.parse_args()
    backends = get_backends(args.pdf)
    print("target-mass (gamma^2) bound on the inclusive A_par chain "
          "-- backends: %s, M = %.4f GeV, g2 = g2_WW"
          % (backends["tag"], M_NUCLEON))
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
