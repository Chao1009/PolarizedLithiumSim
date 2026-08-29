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

The |t| binning is the published seven-bin window of
`recopseudo.T_EDGES_PUBLISHED`, 0.017-0.25 GeV^2 (adopted 2026-08-28);
`--t-edges` takes any increasing list, and the run-13 four-bin window is
`--t-edges 0.05,0.08,0.12,0.17,0.25`.  Below 0.05 GeV^2 the injected a_t
is the linear-in-|t| deformation model extrapolated below the lowest
digitized anchor point (coherent.MANTYSAARI_A2_DEUTERON starts at 0.05).

Panels: (a) delta a_t per |t| bin, 1 yr, both apertures, all three
configurations; (b) the recovered a_t against the injected curve.
`--fit likelihood` swaps the bin-wise ratio for the acceptance-profiled
Poisson likelihood (plans/08 A12); the errors are the same, the recovered
a_t of panel (b) is unbiased in the sparse bins.

A |t| bin is only plotted if it expects at least MIN_TAGGED_PER_BIN
tagged recoils in one year at this run's luminosity share; the rest are
printed as DROPPED, with their expected count, and appear in neither
panel.  The guard exists because the silicon row is a POINT, not a
curve: at 18 x 275 the measured aperture tags 268 recoils a year, 231 of
them in one |t| bin, and the fit of that bin returned
a_t = -1.56 +- 2.22 -- a number with the shape of a measurement and none
of the content, which alone set the vertical range of panel (a) over two
extra decades.  With the guard the silicon aperture is what it is at the
other two configurations: no bins, and therefore no curve.

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
# the published reconstructed |t| binning, shared with
# money_cos2phi_coherent_reco.py (adopted 2026-08-28; plans/08 8.4)
T_EDGES = rp.T_EDGES_PUBLISHED

#: Expected tagged recoils a |t| bin must hold, at ONE YEAR and at this
#: run's luminosity share, before its fit is reported as a measurement
#: (2026-08-28).  The two-azimuth harmonic fit solves seven harmonics on
#: an n_alpha x n_beta = 12 x 24 design in each of the two fills, so a
#: thousand recoils is already under two counts per (alpha, beta) cell;
#: below it the ratio estimator's linearization has nothing to stand on
#: and the number it returns is an artefact of the empty design, not an
#: error bar.  The bin this retires is the only one the measured silicon
#: aperture ever produced -- 0.170 < |t| < 0.250 at 18 x 275, 231 of the
#: 268 tagged/yr there, which returned a_t = -1.56 +- 2.22 against an
#: injected +0.42 and stretched panel (a) over two extra decades.  The
#: guard's count is the acceptance share alone, which is why it reads 231
#: where the fit that used to run on the bin reported N = 258: the fit
#: also carries the model's |t| weight across the bin.  A threshold is
#: not a measurement, so the cheaper count is the right one.  Dropped
#: bins are printed, with their expected count, and are absent from both
#: panels; the surviving-bin counter in the header line counts them out.
MIN_TAGGED_PER_BIN = 1000.0
_CFG = beams.default_configs("6Li")
NAMES = ("5 x 41", "10 x 100", "18 x 275")


def t_edges_for(args):
    """The |t| edges of one run: the published seven by default, any
    increasing list through `--t-edges` (the run-13 window is
    0.05,0.08,0.12,0.17,0.25)."""
    if getattr(args, "t_edges", None) is None:
        return list(T_EDGES)
    edges = [float(v) for v in args.t_edges.split(",")]
    if len(edges) < 2 or any(hi <= lo for lo, hi in zip(edges[:-1],
                                                        edges[1:])):
        raise SystemExit("--t-edges wants at least two increasing edges, "
                         "not %r" % (args.t_edges,))
    return edges


def expected_tagged_in_bin(n_tag, t_reco, tlo, thi):
    """Tagged recoils a reconstructed |t| bin expects in one year.

    The accepted recoils of a `CoherentResponse` carry ONE weight each
    (`CoherentResponse.w` is constant for the shifted-exponential
    importance sampling), so the bin's share of them times the row's
    N_tag is the expected count, and no fit has to be attempted to know
    it.  Returns 0.0 when nothing was accepted at all."""
    t_reco = np.asarray(t_reco, dtype=float)
    if t_reco.size == 0:
        return 0.0
    share = float(((t_reco >= tlo) & (t_reco < thi)).sum()) / t_reco.size
    return float(n_tag) * share


def run(cfg_index, aperture, args, seed):
    """One configuration at one aperture, at the tagging optics -> per-bin
    fit results.  `args.fit` selects the estimator (see
    money_cos2phi_coherent_reco.py --fit)."""
    config = beams.default_configs("6Li")[cfg_index]
    rng = np.random.default_rng(seed)
    sc = coh.CoherentScenario(amp=args.amp, eps_b0=args.eps_b0)
    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            run_share=args.lumi_fraction,
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

    edges = t_edges_for(args)
    out = {"acc": cresp.acceptance, "n_tag": n_produced * cresp.acceptance,
           "cut": cresp.cut_theta_eff, "bins": [], "dropped": [], "top": top,
           "n_bins": len(edges) - 1}
    for tlo, thi in zip(edges[:-1], edges[1:]):
        # MIN_TAGGED_PER_BIN, applied BEFORE the fit is attempted
        n_exp = expected_tagged_in_bin(out["n_tag"], cresp.t_reco, tlo, thi)
        if n_exp < MIN_TAGGED_PER_BIN:
            if cresp.t_reco.size:
                out["dropped"].append((tlo, thi, n_exp))
            continue
        try:
            f = rp.measure_coherent(cresp, n_produced, plan, tlo, thi,
                                    args.amp, a_t_func, u1=args.u1,
                                    u2=args.u2, rng=rng, with_sin=True,
                                    fit=args.fit)
        except (np.linalg.LinAlgError, ValueError, ZeroDivisionError):
            continue
        if not np.isfinite(f["err_t"]):
            continue
        out["bins"].append({"t": f["truth"]["t_ref"], "tlo": tlo, "thi": thi,
                            "a_t": f["a_t"], "d_a_t": f["err_t"],
                            "a_t_truth": f["truth"]["a_t"], "a_e": f["a_e"],
                            "d_a_e": f["err_e"], "n": f["n"]})
    return out


def output_stem(args):
    """File stem for one run.  Report 4's __NB2__ is the DEFAULT run --
    the ratio fit at the published (u1, u2) and the published |t| edges.
    A non-default `--fit` or `--t-edges` gets its key appended rather
    than overwriting the published PNG (the same guard as
    `money_tagged_azz.output_stem`, 2026-08-28; the `--t-edges` half was
    added when the seven-bin window replaced the four-bin one)."""
    base = "nearbeam_reach_gain_6Li"
    keys = []
    if args.fit != "ratio":
        keys.append(args.fit)
    if getattr(args, "t_edges", None) is not None:
        keys.append("tedges")
    share_key = fom.run_share_tag(getattr(args, "lumi_fraction", 1.0))
    if share_key:                      # a run-plan share is a non-default run
        keys.append(share_key)
    return base if not keys else "%s_%s" % (base, "_".join(keys))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lumi-1yr", type=float, default=10.0)
    ap.add_argument("--lumi-fraction", type=float, default=1.0,
                    dest="lumi_fraction",
                    help="this observable's share of the PROGRAMME "
                         "luminosity (plans/07 WP2; default 1.0, which the "
                         "published figure assumes).  A THIRD factor "
                         "beside the tagging optics' own L/L_HA = 1/7-1/13, "
                         "printed beside it, and the 0.5/0.5 spin-state "
                         "share of the flip plan")
    ap.add_argument("--pzz", type=float, default=0.60)
    ap.add_argument("--eps-b0", type=float, default=-0.08)
    ap.add_argument("--amp", type=float, default=0.01)
    ap.add_argument("--u1", type=float, default=0.05)
    ap.add_argument("--u2", type=float, default=0.02)
    ap.add_argument("--n-mc", type=int, default=600000)
    ap.add_argument("--fit", default="ratio", choices=("ratio", "likelihood"),
                    help="estimator of the two-azimuth harmonics: 'ratio' "
                         "(default, the published curve) or 'likelihood', "
                         "the acceptance-profiled Poisson likelihood that "
                         "is unbiased at any count (plans/08 A12).  The two "
                         "give the same ERRORS -- this figure plots errors "
                         "-- and differ in the recovered a_t of panel (b) "
                         "wherever a bin is sparse")
    ap.add_argument("--t-edges", default=None,
                    help="comma-separated reconstructed |t| bin edges in "
                         "GeV^2, replacing the published seven-bin "
                         "0.017,0.028,0.039,0.05,0.08,0.12,0.17,0.25 "
                         "(recopseudo.T_EDGES_PUBLISHED, adopted "
                         "2026-08-28, and the same default as "
                         "money_cos2phi_coherent_reco.py).  The run-13 "
                         "window is --t-edges 0.05,0.08,0.12,0.17,0.25; "
                         "a non-default list appends its key to the "
                         "output stem")
    ap.add_argument("--seed", type=int, default=20260826)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    if not args.lumi_fraction > 0:
        ap.error("--lumi-fraction must be positive")
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    print("run plan: programme %g fb^-1/u/yr x share %g -> %g fb^-1/u "
          "delivered, before the tagging optics' own L/L_HA printed per "
          "configuration below; spin-state share 0.5 / 0.5 within it"
          % (args.lumi_1yr, args.lumi_fraction,
             args.lumi_1yr * args.lumi_fraction))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.5))
    res = {}
    alphas = (1.0, 0.75, 0.5)
    for cfg_index, (name, alpha) in enumerate(zip(NAMES, alphas)):
        cfg = _CFG[cfg_index]
        # keyed on the CONFIGURATION, not on a momentum or on the panel
        # label: the aperture is a property of the ring, and the momentum
        # path of rp_aperture_for cannot resolve 7Li (plans/09 B3)
        meas = reco.rp_aperture_for(cfg)
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
                     len(r["bins"]), r["n_bins"]))
            for tlo, thi, n_exp in r["dropped"]:
                print("   |t| %.3f-%.3f: DROPPED -- %.3g expected tagged "
                      "recoils at 1 yr, below MIN_TAGGED_PER_BIN = %g"
                      % (tlo, thi, n_exp, MIN_TAGGED_PER_BIN))
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
                print("   |t| %.3f-%.3f  N %8.3g  a_t %+.4f +- %.4f "
                      "(truth %+.4f)  a_e %+.4f +- %.4f"
                      % (b["tlo"], b["thi"], b["n"], b["a_t"], b["d_a_t"],
                         b["a_t_truth"], b["a_e"], b["d_a_e"]))

    edges = t_edges_for(args)
    tt = np.linspace(edges[0], edges[-1], 100)
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
    out = outdir / ("%s.png" % output_stem(args))
    fig.savefig(out, dpi=140)
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
