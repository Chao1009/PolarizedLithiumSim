#!/usr/bin/env python3
"""Where the hadronic E - p_z sum goes at the four sweet spots (Report 2 §3).

The Sigma method reconstructs y from Sigma_h = sum(E - p_z) over the
hadronic final state, so what the central detector does NOT see biases
y_Sigma.  This script takes the PYTHIA 8 sample at one beam configuration,
selects the events around the four sweet spots of money plot 5 (the same
selection the money-plot-5 script makes, via reco_chain_figures.sweet_spots) and
splits the TRUE Sigma_h of those events into

  captured   -- tracks (|eta| <= eta_track, p_T > p_T,min), photons and
                neutral hadrons in the calorimeters (|eta| <= eta_cal,
                above threshold)
  forward    -- beyond the calorimeters on the hadron-going side
  thresholds -- charged below p_T,min inside the tracker, photons and
                neutral hadrons below their cluster thresholds
  backward   -- beyond the calorimeters on the electron-going side

for the response's default calorimeter reach (|eta| <= 3.7) and for the
ePIC nominal forward reach (4.0), and runs the full hadron-side response
(HadronResponse, 50 MeV noise) at both to give the Sigma-method delta y/y
at each spot, the captured fraction through the full response (Sigma-
weighted mean and per-event median, noise off) and the median y_Sigma/y.
It also reports the double-angle hadronic angle gamma_h = 2 atan(Sigma/p_T)
as a pseudorapidity -- the struck-quark direction -- and the Sigma-weighted
eta quantiles of the particles, which is where Sigma is actually carried.

Outputs evgen/hfs_acceptance_6Li.png and a printed table.

Usage:  python3 scripts/hfs_acceptance.py --config 1 \\
            --sample samples/pythia8_e10_p99.5_dis.npz samples/pythia8_e10_n99.5_dis.npz
"""

import argparse
import pathlib
import sys

import numpy as np

_SCRIPTS = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS.parent))
sys.path.insert(0, str(_SCRIPTS))
sys.path.insert(0, str(_SCRIPTS.parent.parent / "fastsim"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from polligen import hfs, reco  # noqa: E402
from polli_fastsim import beams  # noqa: E402
from reco_chain_figures import sweet_spots  # noqa: E402

C_TRK, C_GAM, C_NH = "#0072B2", "#E69F00", "#009E73"
C_FWD, C_THR, C_BWD = "#D55E00", "#CC79A7", "0.6"


def event_slice(smp, sel_events):
    """Particle arrays and local offsets for a subset of events."""
    counts = np.diff(smp.offsets)
    ev_of = np.repeat(np.arange(smp.n_events), counts)
    pm = sel_events[ev_of]
    ev_ids = np.flatnonzero(sel_events)
    local = np.searchsorted(ev_ids, ev_of[pm])
    offsets = np.concatenate([[0], np.cumsum(counts[ev_ids])])
    return pm, local, offsets, ev_ids


def decompose(p4, pid, charge, local, n_ev, resp):
    """Fractions of the true Sigma by fate, for the acceptance/threshold
    parameters of `resp` (no efficiencies, no smearing: pure geometry and
    thresholds)."""
    E, px, py, pz = p4[:, 0], p4[:, 1], p4[:, 2], p4[:, 3]
    # acceptance edges live in the detector frame (25 mrad crossing), the
    # E - p_z shares in the head-on frame -- as in HadronResponse
    eta = resp.lab_eta(p4)
    lab = reco.head_on_to_lab(p4, resp.xing) if resp.xing else p4
    pt = np.hypot(lab[:, 1], lab[:, 2])
    emz = E - pz
    apid = np.abs(pid)
    nu = np.isin(apid, hfs.NEUTRINOS)
    is_ch = (charge != 0) & ~nu
    is_gam = apid == 22
    is_nh = ~is_ch & ~is_gam & ~nu
    tot = emz[~nu].sum()
    f = lambda m: float(emz[m & ~nu].sum() / tot)   # noqa: E731
    cap_trk = f(is_ch & (np.abs(eta) <= resp.eta_track) & (pt > resp.pt_min_track))
    # charged particles inside the calorimeters but not tracked are HCal objects
    cap_gam = f(is_gam & (np.abs(eta) <= resp.eta_cal) & (E >= resp.e_min_photon))
    cap_nh = f((is_nh | (is_ch & ((np.abs(eta) > resp.eta_track) | (pt <= resp.pt_min_track))))
               & (np.abs(eta) <= resp.eta_cal) & (E >= resp.e_min_nhad))
    fwd = f(eta > resp.eta_cal)
    bwd = f(eta < -resp.eta_cal)
    thr = 1.0 - cap_trk - cap_gam - cap_nh - fwd - bwd
    # hadronic-system angle and Sigma-weighted eta quantiles
    S = np.bincount(local, weights=emz * (~nu), minlength=n_ev)
    PX = np.bincount(local, weights=px * (~nu), minlength=n_ev)
    PY = np.bincount(local, weights=py * (~nu), minlength=n_ev)
    eta_h = float(-np.log(np.median(S / np.maximum(np.hypot(PX, PY), 1e-9))))
    order = np.argsort(eta[~nu]); c = np.cumsum(emz[~nu][order]) / tot
    q = [float(eta[~nu][order][min(np.searchsorted(c, v), c.size - 1)]) for v in (0.16, 0.5, 0.84)]
    n_ch = float(is_ch.sum() / n_ev)
    return dict(trk=cap_trk, gam=cap_gam, nh=cap_nh, fwd=fwd, thr=thr, bwd=bwd,
                captured=cap_trk + cap_gam + cap_nh, eta_h=eta_h, q=q, n_ch=n_ch)


def sigma_resolution(smp, pm, local, offsets, n_ev, ev_ids, resp, rng, e_e, s, rng_e):
    """68% half-width of y_Sigma / y - 1 through the full response, the
    captured Sigma fraction through it (Sigma-weighted mean, per-event
    median; noise off) and the median y_Sigma / y."""
    p4, pid, ch = smp.p4[pm], smp.pid[pm], smp.charge[pm]
    p4m, w = resp.reconstruct_particles(p4, pid, ch, rng)
    sig, ptx, pty = hfs.hadronic_sums(p4m, offsets, w)
    sig_true, _, _ = hfs.hadronic_sums(p4, offsets, (~np.isin(np.abs(pid), hfs.NEUTRINOS)).astype(float))
    cap_mean = float(sig.sum() / sig_true.sum())
    cap_median = float(np.median(sig / np.maximum(sig_true, 1e-9)))
    if resp.noise_sigma > 0:
        sig = sig + resp.noise_sigma * rng.standard_normal(n_ev)
        ptx = ptx + resp.noise_sigma * rng.standard_normal(n_ev)
        pty = pty + resp.noise_sigma * rng.standard_normal(n_ev)
    kp = smp.kp[ev_ids]
    theta_e = np.arccos(kp[:, 3] / np.sqrt((kp[:, 1:] ** 2).sum(axis=1)))
    e_prime = kp[:, 0]
    eta_lab = -np.log(np.tan(np.minimum(theta_e, np.pi - 1e-9) / 2.0))
    e_r = e_prime * (1.0 + reco.emcal_resolution(e_prime) * rng_e.standard_normal(n_ev))
    th_r = np.clip(theta_e + reco.tracking_angular_resolution(eta_lab)
                   * rng_e.standard_normal(n_ev), 1e-9, np.pi - 1e-9)
    kin = hfs.hadronic_kinematics(sig, ptx, pty, e_r, th_r, e_e, s)
    ratio = kin["y_sigma"] / smp.y[ev_ids]
    d = ratio[np.isfinite(ratio)] - 1.0
    y_e = 1.0 - e_r * (1.0 - np.cos(th_r)) / (2.0 * e_e)
    de = y_e / smp.y[ev_ids] - 1.0
    hw = lambda a: float(0.5 * np.diff(np.percentile(a, [16, 84]))[0])   # noqa: E731
    return dict(dy_sigma=hw(d), dy_e=hw(de[np.isfinite(de)]), cap_mean=cap_mean,
                cap_median=cap_median, y_ratio_median=float(np.median(ratio[np.isfinite(ratio)])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2))
    ap.add_argument("--sample", nargs="+", required=True)
    ap.add_argument("--eta-cal", type=float, nargs="+", default=(3.7, 4.0),
                    help="calorimeter reaches to compare (default: the response's 3.7 "
                         "and the ePIC nominal 4.0)")
    ap.add_argument("--noise", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=20260827)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs("6Li")[args.config]
    smp = hfs.HFSSample.concatenate([hfs.HFSSample.load(f) for f in args.sample])
    if abs(smp.e_energy - config.electron_energy) > 1e-6 or \
            abs(smp.p_per_nucleon - config.ion_momentum_per_nucleon) > 1e-6:
        sys.exit("sample beams (%g, %g) do not match config %s"
                 % (smp.e_energy, smp.p_per_nucleon, config.label()))
    spots = [(float(x), float(q)) for x, q, _i, _j in sweet_spots(args.config)[4]]
    s, e_e = smp.s, config.electron_energy
    print("configuration %s, %d events; spots %s" % (config.label(), smp.n_events,
          ", ".join("(%.3g, %.3g)" % sp for sp in spots)))

    rows = []
    for xs, qs in spots:
        sel = (smp.x > xs / 1.41) & (smp.x < xs * 1.41) & (smp.q2 > qs / 1.29) & (smp.q2 < qs * 1.29)
        pm, local, offsets, ev_ids = event_slice(smp, sel)
        n_ev = ev_ids.size
        y_spot = qs / (xs * s)
        w_spot = (qs * (1 - xs) / xs + 0.938272 ** 2) ** 0.5
        row = dict(x=xs, q2=qs, y=y_spot, W=w_spot, n_ev=n_ev, scen=[])
        for k, eta_cal in enumerate(args.eta_cal):
            resp = hfs.HadronResponse(eta_cal=eta_cal, noise_sigma=args.noise)
            dec = decompose(smp.p4[pm], smp.pid[pm], smp.charge[pm], local, n_ev, resp)
            res = sigma_resolution(smp, pm, local, offsets, n_ev, ev_ids, resp,
                                   np.random.default_rng(args.seed + k), e_e, s,
                                   np.random.default_rng(args.seed + 1000))
            dec.update(eta_cal=eta_cal, **res)
            row["scen"].append(dec)
        rows.append(row)
        d0 = row["scen"][0]
        print("spot x=%.4f Q2=%.3g y=%.3f W=%.1f GeV <n_ch>=%.1f (%d events): DA angle gamma_h as eta = %.2f, "
              "Sigma-weighted particle eta 16/50/84%% = %.1f/%.1f/%.1f"
              % (xs, qs, y_spot, w_spot, d0["n_ch"], n_ev, d0["eta_h"], *d0["q"]))
        for d in row["scen"]:
            print("   cal |eta| <= %.1f: geometry+thresholds captured %.2f (tracks %.2f, photons %.2f, HCal objects %.2f) | "
                  "lost forward %.3f, thresholds %.3f, backward %.3f | full response: captured %.2f (mean) %.2f (median), "
                  "median y_Sigma/y %.2f, dy/y %.2f (electron %.2f)"
                  % (d["eta_cal"], d["captured"], d["trk"], d["gam"], d["nh"], d["fwd"], d["thr"],
                     d["bwd"], d["cap_mean"], d["cap_median"], d["y_ratio_median"], d["dy_sigma"], d["dy_e"]))

    # --- figure ---------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 4.6), gridspec_kw=dict(width_ratios=(1.6, 1)))
    nsc = len(args.eta_cal)
    width = 0.8 / nsc
    labels = []
    for i, row in enumerate(rows):
        labels.append("x=%.3g\nQ²=%.3g\ny=%.3f" % (row["x"], row["q2"], row["y"]))
        for k, d in enumerate(row["scen"]):
            xpos = i + (k - (nsc - 1) / 2) * width
            bottom = 0.0
            for key, col, lab in (("trk", C_TRK, "captured: tracks"), ("gam", C_GAM, "captured: photons"),
                                  ("nh", C_NH, "captured: neutral hadrons"),
                                  ("fwd", C_FWD, "lost forward, beyond the calorimeters"),
                                  ("thr", C_THR, "lost below thresholds"), ("bwd", C_BWD, "lost backward")):
                ax1.bar(xpos, d[key], width * 0.92, bottom=bottom, color=col,
                        alpha=1.0 if k == 0 else 0.55, label=lab if (i == 0 and k == 0) else None)
                bottom += d[key]
            ax1.text(xpos, 1.01, "|η|≤%.1f" % d["eta_cal"], ha="center", va="bottom", fontsize=6.5)
    ax1.set_xticks(range(len(rows))); ax1.set_xticklabels(labels, fontsize=7.5)
    ax1.set_ylim(0, 1.08); ax1.set_ylabel("fraction of the true Σ$_h$ = Σ(E − p$_z$)")
    ax1.set_title("(a) where the hadronic E − p$_z$ sum goes: geometry and thresholds "
                  "(calorimeter reach as labelled; lighter bars: 4.0)", fontsize=9)
    handles, labels_ = ax1.get_legend_handles_labels()
    for k, eta_cal in enumerate(args.eta_cal):
        xpos = np.arange(len(rows)) + (k - (nsc - 1) / 2) * width
        ax2.bar(xpos, [r["scen"][k]["dy_sigma"] for r in rows], width * 0.92,
                color=C_TRK, alpha=1.0 if k == 0 else 0.55, label="Σ method, cal |η| ≤ %.1f" % eta_cal)
    ax2.plot(range(len(rows)), [r["scen"][0]["dy_e"] for r in rows], "s", color="0.4", ms=5,
             label="electron method")
    ax2.axhline(0.25, color="0.7", lw=0.9, ls="-.")
    ax2.set_xticks(range(len(rows))); ax2.set_xticklabels(labels, fontsize=7.5)
    ax2.set_yscale("log"); ax2.set_ylim(0.05, 3.0)
    ax2.set_ylabel(r"$\delta y/y$ (68% half-width)")
    ax2.set_title("(b) Σ-method resolution, %.0f MeV noise" % (1e3 * args.noise), fontsize=9)
    ax2.legend(fontsize=7)
    fig.suptitle("Hadronic-final-state acceptance at the four sweet spots — %s, PYTHIA 8; tracker |η| ≤ 3.5, "
                 "p$_T$ > 0.2 GeV; photons > 0.1, neutral hadrons > 0.5 GeV" % config.label(), fontsize=9.5)
    fig.legend(handles, labels_, loc="lower center", ncol=3, fontsize=7.2, frameon=False)
    fig.tight_layout(rect=(0, 0.1, 1, 0.94))
    out = pathlib.Path(args.outdir) / "hfs_acceptance_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)


if __name__ == "__main__":
    main()
