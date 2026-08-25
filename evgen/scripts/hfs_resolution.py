#!/usr/bin/env python3
"""Hadronic-method resolution with a hadronic final state and the
hadron-side detector response (plans/07 WP3-HFS).

For an HFS sample (PYTHIA 8 via tools/pythia8/gen_dis_hfs.py, or the toy
string fragmentation of polligen.hfs when no sample is given) at the
mid-energy configuration, the script applies hfs.HadronResponse and shows

  (a) delta y/y of the Sigma, Jacquet-Blondel and double-angle methods
      versus y at Q2 < 3 GeV2 (the sweet-spot slice), for the default
      response and for a noise scan, against the 25% Gaussian stand-in;
  (b) the captured fraction of Sigma = sum(E - p_z) versus y and its
      decomposition (tracks / photons / neutral hadrons / lost);
  (c) the same resolutions at the four sweet spots of money plot 5,
      with the electron-method delta y/y for comparison.

Outputs evgen/hfs_resolution_6Li.png and a printed table.  With the toy
generator the numbers are ILLUSTRATIVE (labelled "toy"); the PYTHIA sample
is the deliverable of tools/pythia8/README.md.

Usage:  python3 scripts/hfs_resolution.py [--sample a.npz b.npz] [--noise 0 0.025 0.05 0.1]
"""

import argparse
import pathlib
import sys

import numpy as np

_SCRIPTS = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS.parent))
sys.path.insert(0, str(_SCRIPTS))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from polligen import hfs, reco, recopseudo as rp  # noqa: E402
from polligen.sample import InclusiveSampler  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.polarized import toy_b1  # noqa: E402

SPOTS = ((0.0562, 1.14), (0.0224, 1.14), (0.141, 3.13), (0.141, 14.3))
C = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9")


def resolution_by_y(y_true, y_reco, edges):
    """Half 68% width of y_reco/y_true - 1 (robust) and rms per y bin."""
    hw, rms, mu, cen = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (y_true >= lo) & (y_true < hi) & np.isfinite(y_reco)
        d = y_reco[m] / y_true[m] - 1.0
        d = d[np.isfinite(d)]
        if d.size < 30:
            hw.append(np.nan); rms.append(np.nan); mu.append(np.nan)
        else:
            q16, q50, q84 = np.percentile(d, [16, 50, 84])
            hw.append(0.5 * (q84 - q16)); rms.append(np.sqrt(np.mean(d ** 2))); mu.append(q50)
        cen.append(np.sqrt(lo * hi))
    return np.array(cen), np.array(hw), np.array(rms), np.array(mu)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2))
    ap.add_argument("--sample", nargs="*", default=None)
    ap.add_argument("--n-toy", type=int, default=400000)
    ap.add_argument("--noise", type=float, nargs="*", default=(0.0, 0.025, 0.05, 0.1))
    ap.add_argument("--q2max", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=20260826)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    config = beams.default_configs("6Li")[args.config]
    analysis = fom.Scenario(lumi_fb_per_nucleon=10.0, pol_ion_tensor=0.6)
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    sampler = InclusiveSampler(kern, config, rp.generator_scenario(analysis),
                               nx=60, nq2=45, q2_range=(0.7, 2e3))
    if args.sample:
        smp = hfs.HFSSample.concatenate([hfs.HFSSample.load(f) for f in args.sample])
        src = "%s, %d events" % (smp.meta.get("generator", "sample"), smp.n_events)
        toy = False
    else:
        smp = hfs.toy_library_sample(sampler, args.n_toy, rng, per_cell_uniform=False)
        src = "TOY string fragmentation, %d events" % smp.n_events
        toy = True
    if abs(smp.e_energy - config.electron_energy) > 1e-6 or \
            abs(smp.p_per_nucleon - config.ion_momentum_per_nucleon) > 1e-6:
        sys.exit("sample beams (%g, %g) do not match config %s"
                 % (smp.e_energy, smp.p_per_nucleon, config.label()))
    s = smp.s
    e_e = config.electron_energy
    theta_e = np.arccos(smp.kp[:, 3] / np.sqrt((smp.kp[:, 1:] ** 2).sum(axis=1)))
    e_prime = smp.kp[:, 0]
    # reconstructed electron (EMCal + track angles) for the hadronic methods
    eta_lab = -np.log(np.tan(np.minimum(theta_e, np.pi - 1e-9) / 2.0))
    de = reco.emcal_resolution(e_prime)
    sig_ang = reco.tracking_angular_resolution(eta_lab)
    e_r = e_prime * (1.0 + de * rng.standard_normal(e_prime.size))
    th_r = np.clip(theta_e + sig_ang * rng.standard_normal(e_prime.size), 1e-9, np.pi - 1e-9)
    q2_e = 2.0 * e_e * e_r * (1.0 + np.cos(th_r))
    y_e = 1.0 - e_r * (1.0 - np.cos(th_r)) / (2.0 * e_e)

    slice_mask = smp.q2 < args.q2max
    y_edges = np.logspace(np.log10(0.004), np.log10(0.95), 13)
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.4))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    # (a) methods at the default response, and (d) noise scan for Sigma
    table = []
    results = {}
    for k, noise in enumerate(args.noise):
        resp = hfs.HadronResponse(noise_sigma=noise)
        r = resp.reconstruct_sums(smp, rng=np.random.default_rng(args.seed + k))
        kin = hfs.hadronic_kinematics(r["sigma"], r["ptx"], r["pty"], e_r, th_r, e_e, s)
        results[noise] = (r, kin)
        cen, hw, rms, mu = resolution_by_y(smp.y[slice_mask], kin["y_sigma"][slice_mask], y_edges)
        ax_d.plot(cen, hw, "o-", color=C[k % len(C)], ms=3.5, lw=1.2,
                  label=r"$\Sigma$ method, noise %.0f MeV" % (1e3 * noise))
        table.append((noise, cen, hw, rms, mu))
    noise0 = args.noise[min(2, len(args.noise) - 1)]
    r, kin = results[noise0]
    for k, (key, lab) in enumerate((("y_sigma", r"$\Sigma$"), ("y_jb", "Jacquet–Blondel"),
                                    ("y_da", "double angle"))):
        cen, hw, rms, mu = resolution_by_y(smp.y[slice_mask], kin[key][slice_mask], y_edges)
        ax_a.plot(cen, hw, "o-", color=C[k], ms=3.5, lw=1.2, label="%s: 68%% half-width" % lab)
        ax_a.plot(cen, rms, ":", color=C[k], lw=1.0)
    cen, hw, rms, mu = resolution_by_y(smp.y[slice_mask], y_e[slice_mask], y_edges)
    ax_a.plot(cen, hw, "s--", color="0.4", ms=3, lw=1.0, label="electron method (EMCal 2%/√E ⊕ 1%)")
    ax_a.axhline(0.25, color="0.7", lw=0.9, ls="-.", label="25% Gaussian stand-in (5R/7R)")
    ax_a.set_xscale("log"); ax_a.set_yscale("log")
    ax_a.set_ylim(0.02, 3.0)
    ax_a.set_xlabel(r"$y$ (true)"); ax_a.set_ylabel(r"$\delta y/y$")
    ax_a.set_title(r"(a) hadronic methods at $Q^2 < %g$ GeV$^2$, noise %.0f MeV "
                   r"(dotted: rms)" % (args.q2max, 1e3 * noise0), fontsize=9)
    ax_a.legend(fontsize=7, loc="upper right")
    ax_d.axhline(0.25, color="0.7", lw=0.9, ls="-.")
    ax_d.set_xscale("log"); ax_d.set_yscale("log"); ax_d.set_ylim(0.02, 3.0)
    ax_d.set_xlabel(r"$y$ (true)"); ax_d.set_ylabel(r"$\delta y_\Sigma/y$ (68% half-width)")
    ax_d.set_title(r"(d) noise scan: $\Sigma_h = 2E_e y$ = %.2f GeV at $y$ = 0.01"
                   % (2 * e_e * 0.01), fontsize=9)
    ax_d.legend(fontsize=7)

    # (b) captured fraction of Sigma and its decomposition by detector object
    resp = hfs.HadronResponse(noise_sigma=0.0)
    rng_b = np.random.default_rng(args.seed + 99)
    p4m, w = resp.reconstruct_particles(smp.p4, smp.pid, smp.charge, rng_b)
    ev = smp.event_index()
    empz_true = smp.p4[:, 0] - smp.p4[:, 3]
    is_nu = np.isin(np.abs(smp.pid), hfs.NEUTRINOS)
    sig_true = np.bincount(ev, np.where(is_nu, 0.0, empz_true), smp.n_events)
    empz_reco = (p4m[:, 0] - p4m[:, 3]) * w
    charged = (smp.charge != 0)
    photon = smp.pid == 22
    parts = {"tracks / charged": charged, "photons": photon,
             "neutral hadrons": ~charged & ~photon & ~is_nu}
    bottom = np.zeros(y_edges.size - 1)
    ycen = np.sqrt(y_edges[:-1] * y_edges[1:])
    for k, (lab, sel) in enumerate(parts.items()):
        frac = []
        for lo, hi in zip(y_edges[:-1], y_edges[1:]):
            m = (smp.y >= lo) & (smp.y < hi) & slice_mask
            if m.sum() < 30:
                frac.append(np.nan); continue
            num = np.bincount(ev, np.where(sel, empz_reco, 0.0), smp.n_events)[m].sum()
            frac.append(num / sig_true[m].sum())
        frac = np.array(frac)
        ax_b.bar(np.arange(ycen.size), np.nan_to_num(frac), bottom=bottom, color=C[k],
                 label=lab, width=0.8)
        bottom = bottom + np.nan_to_num(frac)
    ax_b.set_xticks(np.arange(ycen.size))
    ax_b.set_xticklabels(["%.3g" % v for v in ycen], fontsize=7, rotation=45)
    ax_b.axhline(1.0, color="0.5", lw=0.8)
    ax_b.set_xlabel(r"$y$ (true)"); ax_b.set_ylabel(r"reconstructed / true $\Sigma$ (rate-weighted)")
    ax_b.set_title("(b) captured fraction of Σ by detector object (no noise)", fontsize=9)
    ax_b.legend(fontsize=7, loc="lower right")
    ax_b.set_ylim(0, 1.15)

    # (c) resolutions at the four sweet spots (cells around the spot)
    r, kin = results[noise0]
    labels, vals = [], {"y_sigma": [], "y_jb": [], "y_da": [], "electron": []}
    for (xs, qs) in SPOTS:
        m = ((smp.x > xs / 1.41) & (smp.x < xs * 1.41) & (smp.q2 > qs / 1.29) & (smp.q2 < qs * 1.29))
        yspot = qs / (s * xs)
        labels.append("x=%.3g\nQ²=%.3g\ny=%.3f" % (xs, qs, yspot))
        for key in ("y_sigma", "y_jb", "y_da"):
            d = kin[key][m] / smp.y[m] - 1.0
            d = d[np.isfinite(d)]
            vals[key].append(0.5 * np.diff(np.percentile(d, [16, 84]))[0] if d.size > 30 else np.nan)
        d = y_e[m] / smp.y[m] - 1.0
        vals["electron"].append(0.5 * np.diff(np.percentile(d, [16, 84]))[0] if d.size > 30 else np.nan)
    xpos = np.arange(len(SPOTS))
    for k, (key, lab) in enumerate((("y_sigma", "Σ"), ("y_jb", "JB"), ("y_da", "DA"),
                                    ("electron", "electron"))):
        ax_c.bar(xpos + (k - 1.5) * 0.2, vals[key], width=0.2, color=C[k] if k < 3 else "0.5",
                 label=lab)
    ax_c.axhline(0.25, color="0.7", lw=0.9, ls="-.")
    ax_c.set_xticks(xpos); ax_c.set_xticklabels(labels, fontsize=7)
    ax_c.set_yscale("log"); ax_c.set_ylim(0.02, 3.0)
    ax_c.set_ylabel(r"$\delta y/y$ (68% half-width)")
    ax_c.set_title("(c) at the four sweet spots of money plot 5 (noise %.0f MeV)" % (1e3 * noise0),
                   fontsize=9)
    ax_c.legend(fontsize=7)

    fig.suptitle("Hadronic-method resolution with a hadronic final state — %s; %s\n%s"
                 % (config.label(), src, hfs.HadronResponse(noise_sigma=noise0).describe()),
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "hfs_resolution_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print("source:", src, "| toy" if toy else "")
    print("response:", hfs.HadronResponse(noise_sigma=noise0).describe())
    print("delta y/y (68%% half-width) of the Sigma method at Q2 < %g by y bin, per noise level:" % args.q2max)
    print("  y-bin centre: " + " ".join("%7.4f" % v for v in table[0][1]))
    for noise, cen, hw, rms, mu in table:
        print("  noise %3.0f MeV: " % (1e3 * noise) + " ".join("%7.3f" % v for v in hw))
    print("sweet spots (noise %.0f MeV): " % (1e3 * noise0)
          + "; ".join("x=%.3g Q2=%.3g: Sigma %.2f JB %.2f DA %.2f e %.2f"
                      % (xs, qs, a, b, c, d) for (xs, qs), a, b, c, d in
                      zip(SPOTS, vals["y_sigma"], vals["y_jb"], vals["y_da"], vals["electron"])))


if __name__ == "__main__":
    main()
