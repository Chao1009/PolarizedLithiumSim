#!/usr/bin/env python3
"""What a near-beam layer buys the coherent measurement (plans/09;
re-derived 2026-08-28 at the tagging optics of Report 1 Section 6.1).

nearbeam_aperture_scan.py gives the acceptance versus aperture; this
script runs the FULL coherent chain -- Roman-Pot emulation with the
anisotropic divergence, the two azimuths, the spin-state-sorted 2-D
harmonic fit -- at two apertures per configuration and reports what the
measurement itself does:

  silicon    the pots fixed at the aperture measured by shooting an intact
             6Li through the ePIC geometry (reco.RP_APERTURE_MEASURED,
             tools/fullsim), with the beam de-squeezed to the tagging optics
  near beam  the pots following the 10 sigma envelope of the tagging optics
             in both planes -- the closest a detector is allowed to
             approach, and therefore the ceiling on what any technology
             can buy

Both run at the tagging optics' divergence and luminosity (1/7 - 1/13 of
the high-acceptance value).  With the Yellow Report optics and the
measured aperture no recoil is tagged at any configuration (Report 2
Figure 2c), and the de-squeeze alone does not change that: the envelope
shrinks to 0.12-0.33 mrad but the silicon still starts at 1.0-2.0 mrad.
The near-beam layer and the optics are strictly multiplicative levers,
and this script prices the layer GIVEN the optics.

Panels: (a) delta a_t per |t| bin, 1 yr, both apertures, all three
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
from polli_fastsim import farforward as ff  # noqa: E402
from polli_fastsim.farforward import HIGH_ACCEPTANCE  # noqa: E402

C_SI, C_NB = "#C0392B", "#1F4E79"
T_EDGES = (0.05, 0.08, 0.12, 0.17, 0.25)
_CFG = beams.default_configs("6Li")
NAMES = ("5 x 41", "10 x 100", "18 x 275")


def run(cfg_index, aperture, args, seed):
    """One configuration at one aperture, at the tagging optics -> per-bin
    fit results."""
    config = beams.default_configs("6Li")[cfg_index]
    rng = np.random.default_rng(seed)
    sc = coh.CoherentScenario(amp=args.amp, eps_b0=args.eps_b0)
    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    top = ff.tagging_optics_point(config, slope_b=sc.slope_b)
    _, n_coh, _ = coh.project_coherent(config, scenario, sc,
                                       optics_list=(HIGH_ACCEPTANCE,))
    n_produced = float(n_coh.sum()) * top["lumi_fraction"]
    cresp = rp.CoherentResponse(sc, config, (top["sigma_x_eff"], top["sigma_y"]),
                                n_mc=args.n_mc, rng=rng, cut_scale_xy=(1.0, 1.0),
                                cut_theta_xy=aperture)
    plan = bk.tensor_flip_plan(args.pzz)

    def a_t_func(t):
        return sc.cos2phi_coefficient_deformation(t, 1.0)

    out = {"acc": cresp.acceptance, "n_tag": n_produced * cresp.acceptance,
           "cut": cresp.cut_theta_eff, "bins": [], "top": top}
    for tlo, thi in zip(T_EDGES[:-1], T_EDGES[1:]):
        try:
            f = rp.measure_coherent(cresp, n_produced, plan, tlo, thi,
                                    args.amp, a_t_func, u1=args.u1,
                                    u2=args.u2, rng=rng, with_sin=True)
        except (np.linalg.LinAlgError, ValueError, ZeroDivisionError):
            continue
        if not np.isfinite(f["err_t"]):
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
    ap.add_argument("--n-mc", type=int, default=600000)
    ap.add_argument("--seed", type=int, default=20260826)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.5))
    res = {}
    alphas = (1.0, 0.75, 0.5)
    for cfg_index, (name, alpha) in enumerate(zip(NAMES, alphas)):
        cfg = _CFG[cfg_index]
        key = name.replace(" ", "")
        meas = reco.RP_APERTURE_MEASURED[key]
        top = ff.tagging_optics_point(cfg)
        near = (top["env_x"], top["env_y"])
        for tag, apert, col, ls, mk in (("silicon", meas, C_SI, "--", "s"),
                                        ("pots follow", near, C_NB, "-", "o")):
            r = run(cfg_index, apert, args, args.seed)
            res[(name, tag)] = r
            print("%-9s %-12s cutout %.2f x %.2f mrad at sigma %.0f/%.0f urad, "
                  "L/L_HA = 1/%.1f -> acc %.3e, N_tag %.3g/yr, %d of %d |t| bins survive"
                  % (name, tag, 1e3 * r["cut"][0], 1e3 * r["cut"][1],
                     1e6 * r["top"]["sigma_x_eff"], 1e6 * r["top"]["sigma_y"],
                     1.0 / r["top"]["lumi_fraction"], r["acc"], r["n_tag"],
                     len(r["bins"]), len(T_EDGES) - 1))
            if not r["bins"]:
                continue
            t = [b["t"] for b in r["bins"]]
            ax1.plot(t, [b["d_a_t"] for b in r["bins"]], ls, marker=mk,
                     color=col, ms=5, lw=1.5, alpha=alpha,
                     label="%s, %s" % (name, tag))
            ax2.errorbar(t, [b["a_t"] for b in r["bins"]],
                         yerr=[b["d_a_t"] for b in r["bins"]], fmt=mk,
                         color=col, ms=5, lw=1.2, capsize=2, alpha=alpha,
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
    ax1.set_title("(a) what the aperture costs the measurement, at the tagging optics",
                  fontsize=9.5)
    ax2.set_ylim(0.0, 0.4)
    ax2.set_ylabel(r"recovered $a_t$")
    ax2.set_title("(b) recovery against the injected deformation term", fontsize=9.5)
    for ax in (ax1, ax2):
        ax.set_xlabel(r"$|t|$ [GeV$^2$]")
        ax.grid(alpha=0.25, lw=0.5)
        ax.legend(fontsize=7.0)
    fig.suptitle("Coherent intact-⁶Li chain at the tagging optics: pots fixed at the "
                 "measured silicon aperture against pots following the 10σ envelope",
                 fontsize=9.5)
    fig.tight_layout()
    out = outdir / "nearbeam_reach_gain_6Li.png"
    fig.savefig(out, dpi=140)
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
