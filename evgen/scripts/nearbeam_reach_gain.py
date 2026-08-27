#!/usr/bin/env python3
"""What a near-beam layer buys the coherent measurement (plans/09).

nearbeam_aperture_scan.py gives the acceptance versus aperture; this
script runs the FULL coherent chain -- Roman-Pot emulation, the two
azimuths, the spin-state-sorted 2-D harmonic fit -- at two apertures per
configuration and reports what the measurement itself does:

  silicon    the aperture measured by shooting an intact 6Li through the
             ePIC geometry (reco.RP_APERTURE_MEASURED, tools/fullsim)
  near beam  the same, with the HORIZONTAL half-width replaced by the
             10 sigma beam envelope -- the closest a detector is allowed
             to approach, and therefore the ceiling on what any
             technology can buy

The envelope cutout is held at 10 sigma in BOTH axes here (cut_scale_xy
= (1, 1)), unlike money plot 6R, whose 2.5 in x came from the pre-
measurement belief that the pots surround a wide horizontal slot.  With
the geometric aperture measured, carrying 2.5 as well would impose a 25
sigma horizontal retraction that no machine requirement asks for, and it
binds before the geometry does -- hiding the whole effect.

Panels: (a) delta a_t per |t| bin, 1 yr, both apertures, both live
configurations; (b) the recovered a_t against the injected curve.

Usage:  python3 scripts/nearbeam_reach_gain.py [--outdir .]
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from polligen import bookkeeping as bk  # noqa: E402
from polligen import coherent as coh  # noqa: E402
from polligen import reco, recopseudo as rp  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.farforward import HIGH_ACCEPTANCE  # noqa: E402

C_SI, C_NB = "#C0392B", "#1F4E79"
T_EDGES = (0.05, 0.08, 0.12, 0.17, 0.25)
LIVE = ((0, "5 x 41", 20.5), (1, "10 x 100", 50.0))


def run(cfg_index, aperture, args, seed):
    """One configuration at one aperture -> per-bin fit results."""
    config = beams.default_configs("6Li")[cfg_index]
    rng = np.random.default_rng(seed)
    sc = coh.CoherentScenario(amp=args.amp, eps_b0=args.eps_b0)
    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    _, n_coh, _ = coh.project_coherent(config, scenario, sc,
                                       optics_list=(HIGH_ACCEPTANCE,),
                                       sigma_theta_list=(args.sigma_theta,))
    cresp = rp.CoherentResponse(sc, config, args.sigma_theta, n_mc=args.n_mc,
                                rng=rng, cut_scale_xy=(1.0, 1.0),
                                cut_theta_xy=aperture)
    plan = bk.tensor_flip_plan(args.pzz)

    def a_t_func(t):
        return sc.cos2phi_coefficient_deformation(t, 1.0)

    out = {"acc": cresp.acceptance, "bins": []}
    for tlo, thi in zip(T_EDGES[:-1], T_EDGES[1:]):
        try:
            f = rp.measure_coherent(cresp, float(n_coh.sum()), plan, tlo, thi,
                                    args.amp, a_t_func, u1=args.u1,
                                    u2=args.u2, rng=rng, with_sin=True)
        except np.linalg.LinAlgError:
            continue
        out["bins"].append({"t": f["truth"]["t_ref"], "tlo": tlo, "thi": thi,
                            "a_t": f["a_t"], "d_a_t": f["err_t"],
                            "a_t_truth": f["truth"]["a_t"], "a_e": f["a_e"],
                            "d_a_e": f["err_e"], "n": f["n"]})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lumi-1yr", type=float, default=10.0)
    ap.add_argument("--pzz", type=float, default=0.60)
    ap.add_argument("--eps-b0", type=float, default=-0.08)
    ap.add_argument("--amp", type=float, default=0.01)
    ap.add_argument("--u1", type=float, default=0.05)
    ap.add_argument("--u2", type=float, default=0.02)
    ap.add_argument("--sigma-theta", type=float, default=reco.SIGMA_THETA_HA)
    ap.add_argument("--near-beam-mrad", type=float,
                    default=1e3 * 10.0 * reco.SIGMA_THETA_HA)
    ap.add_argument("--n-mc", type=int, default=600000)
    ap.add_argument("--seed", type=int, default=20260826)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.5))
    res = {}
    for cfg_index, name, pu in LIVE:
        key = name.replace(" ", "")
        meas = reco.RP_APERTURE_MEASURED[key]
        near = (1e-3 * args.near_beam_mrad, meas[1])
        for tag, apert, col, ls, mk in (("silicon", meas, C_SI, "--", "s"),
                                        ("near beam", near, C_NB, "-", "o")):
            r = run(cfg_index, apert, args, args.seed)
            res[(name, tag)] = r
            print("%-9s %-10s aperture %.3f mrad -> acc %.3e, %d of %d "
                  "|t| bins survive"
                  % (name, tag, 1e3 * apert[0], r["acc"], len(r["bins"]),
                     len(T_EDGES) - 1))
            if not r["bins"]:
                continue
            t = [b["t"] for b in r["bins"]]
            ax1.plot(t, [b["d_a_t"] for b in r["bins"]], ls, marker=mk,
                     color=col, ms=5, lw=1.5, alpha=0.55 if pu > 30 else 1.0,
                     label="%s, %s" % (name, tag))
            ax2.errorbar(t, [b["a_t"] for b in r["bins"]],
                         yerr=[b["d_a_t"] for b in r["bins"]], fmt=mk,
                         color=col, ms=5, lw=1.2, capsize=2,
                         alpha=0.55 if pu > 30 else 1.0,
                         label="%s, %s" % (name, tag))
            for b in r["bins"]:
                print("   |t| %.2f-%.2f  N %8.3g  a_t %+.4f +- %.4f "
                      "(truth %+.4f)  a_e %+.4f +- %.4f"
                      % (b["tlo"], b["thi"], b["n"], b["a_t"], b["d_a_t"],
                         b["a_t_truth"], b["a_e"], b["d_a_e"]))

    tt = np.linspace(T_EDGES[0], T_EDGES[-1], 100)
    sc = coh.CoherentScenario(amp=args.amp, eps_b0=args.eps_b0)
    ax2.plot(tt, [sc.cos2phi_coefficient_deformation(t, 1.0) for t in tt],
             "-", color="0.35", lw=1.2, label="injected a$_t$(|t|)")
    ax1.set_yscale("log")
    ax1.set_ylabel(r"$\delta a_t$ per $|t|$ bin, 1 yr")
    ax1.set_title("(a) what the aperture costs the measurement", fontsize=9.5)
    ax2.set_ylim(0.0, 0.55)
    for (name, tag), r in res.items():
        for b in r["bins"]:
            if b["a_t"] + b["d_a_t"] > 0.55 or b["a_t"] - b["d_a_t"] < 0.0:
                ax2.annotate("%.2f ± %.2f" % (b["a_t"], b["d_a_t"]),
                             xy=(b["t"], 0.52), fontsize=6.5, color=C_SI,
                             ha="center")
    ax2.set_ylabel(r"recovered $a_t$")
    ax2.set_title("(b) recovery against the injected deformation term",
                  fontsize=9.5)
    for ax in (ax1, ax2):
        ax.set_xlabel(r"$|t|$ [GeV$^2$]")
        ax.grid(alpha=0.25, lw=0.5)
        ax.legend(fontsize=7.5)
    fig.suptitle("Coherent intact-⁶Li chain at the measured Roman-Pot silicon "
                 "aperture and at the 10σ beam envelope", fontsize=9.5)
    fig.tight_layout()
    out = outdir / "nearbeam_reach_gain_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
