#!/usr/bin/env python3
"""What the far-forward near-beam aperture is worth (plans/09; re-derived
2026-08-28 on the Yellow Report divergences of plans/10).

The coherent intact-6Li tag and the 6Li alpha tag are both beam-blind:
A/Z = 2 exactly for the intact nucleus and a rigidity ratio of 0.99813
for the alpha, so neither is separated from the beam by dispersion.  What
separates them is the transverse ANGLE at the interaction point, and the
acceptance is therefore exp(-B t) evaluated at t = (theta p)^2 -- an
exponential in the SQUARE of the angle times the square of the momentum.
Approach distance is the whole measurement.

Three horizontal half-widths matter at each configuration, and they are
not the same number:

  * the 10 sigma beam envelope of the Yellow Report optics -- an
    operational retraction, a machine-protection question rather than a
    detector one: 2.2 / 1.8 / 0.92 mrad at 5 x 41 / 10 x 100 / 18 x 275
    for a 6Li at the high-acceptance optics (farforward.sigma_theta_for);
  * the ePIC Roman-Pot silicon APERTURE, measured by shooting an intact
    6Li through the geometry: 2.50 / 1.51 / 0.53 mrad, re-measured
    2026-08-28 in the CURRENT ePIC geometry (tools/fullsim, epic-main
    9aaa2969, plans/09 B1) -- the sensor package, not the beam.  The
    September-2024 2.0 / 1.35 / 1.03 mrad every number published before
    that date was priced at are still reachable as
    `reco.RP_APERTURE_SEP2024`;
  * the envelope of the lithium TAGGING optics of Report 1 Section 6.1 --
    the horizontal beta* raised to the optimum of acceptance x luminosity,
    0.33 / 0.17 / 0.12 mrad, at 1/7 - 1/13 of the high-acceptance
    luminosity (farforward.tagging_optics_point).

Which of the three binds is now a PER-CONFIGURATION question, and the
re-measurement changed the answer.  At 5 x 41 the silicon binds: 2.50 mrad
against the 2.20 mrad envelope, 1.14x outside it, and the coherent tag
reads 9.4e-10 at the silicon against 7.2e-8 at the envelope.  At 10 x 100
the machine still binds, 1.51 against 1.80 mrad.  At 18 x 275 the machine
binds by nine orders of magnitude -- the aperture halved to 0.53 mrad
against a 0.92 mrad envelope, so the sensor that used to be the limit
there no longer is.  "At the published optics the machine binds
everywhere" was true of the September-2024 table and is not true of this
one.  Under the tagging optics the envelope shrinks by 7-11x and the
sensor package IS the limit at every configuration: a near-beam layer
that reaches the envelope is what makes the optics worth having.  This
script prices every aperture per configuration and marks the three.

`--isotope` (plans/09 B3) selects the species of PANEL (b) alone.  For 7Li
the same scan is flat: the 7Li alpha is off rigidity at R = 0.856, so it
is accepted by the Roman-Pot momentum window and never has to clear the
near-beam cutout -- the curve sits at 0.96-0.99 across the whole
0.05-3 mrad axis and the three marked half-widths coincide.  That is the
picture of a tag that does not care about the aperture, and it is why a
near-beam layer is worth nothing to 7Li.

PANEL (a) IS ALWAYS 6Li, whatever `--isotope` says.  The coherent channel
is 6Li-specific by construction (`polligen.coherent`: a J = 1 nucleus
whose quadrupole deformation is scaled from the deuteron's).  7Li has
J = 3/2 and a quadrupole moment ~50x larger, so its coherent cos 2phi
amplitude is a different physics case with its own model, not a re-run of
this one -- recorded as an open item in plans/09 B3 rather than
half-built here.

Usage:  python3 scripts/nearbeam_aperture_scan.py [--isotope 7Li] [--outdir .]
"""

import argparse
import pathlib
import sys

import numpy as np

_SCRIPTS = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS.parent))
sys.path.insert(0, str(_SCRIPTS.parent.parent / "fastsim"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from polligen import reco  # noqa: E402
from polli_fastsim import beams  # noqa: E402
from polli_fastsim import farforward as ff  # noqa: E402
from polli_fastsim import spectator as sp  # noqa: E402

C_TRUTH, C_FIT, C_ALT, C_GREY = "#1F4E79", "#C0392B", "#B8860B", "#8A8A8A"
_SUP = {"6Li": "⁶Li", "7Li": "⁷Li"}
R12 = ff.POT_R12       # m, IP angle -> pot-plane x at 18 x 275 (measured)


def coherent_acceptance(theta_x, p_ion, slope_b, theta_y):
    """Tagged fraction of an exp(-B t) recoil outside a rectangular
    cutout of half-widths (theta_x, theta_y) in ANGLE."""
    return np.array([reco.rp_hole_acceptance(slope_b, t * p_ion,
                                             theta_y * p_ion)["acc"]
                     for t in np.atleast_1d(theta_x)])


ALPHA_CHANNEL = {"6Li": sp.LI6_ALPHA_TAG, "7Li": sp.LI7_ALPHA_TAG}


def alpha_acceptance(theta_x, theta_y, p_per_nucleon, channel=None,
                     n=200000, beta=0.30, seed=7, pot_config="18x275"):
    """Tagged fraction of alpha spectators in ANY far-forward system,
    with the near-beam band (|R - 1| < 0.05) cut by a rectangular envelope
    of half-widths (theta_x, theta_y) in angle -- the same routing
    (farforward.acceptance_summary) as fastsim/scripts/tagging_acceptance.py
    (Report 3 Table 6), but NOT at the same vertical half-width: this
    script's callers pass max(measured aperture, 10 sigma_v), which since
    the 2026-08-28 re-measurement is the ENVELOPE at 18 x 275 (0.917
    against the aperture's 0.92 mrad, so the two now coincide there) and
    the shut vertical plane at 5 x 41 (8.84 against a 3.80 mrad
    envelope).  The two agree wherever the envelope is the binding
    vertical constraint.
    For 6Li the off-rigidity slice below R = 0.95 (~1.5%) is tagged at any
    envelope and the rest is the R ~ 1 tail outside the rectangle; for 7Li
    the alpha is at R = 0.856, entirely inside the momentum window, and the
    curve is flat.

    `pot_config` is the machine configuration whose blind block the
    OVER-RIGID branch of the routing tests against (48 / 32 / 16 mm at
    5 x 41 / 10 x 100 / 18 x 275).  It has to be passed: leaving it at
    the 18 x 275 default while sweeping configurations lets the R > 1.05
    tail onto silicon that is retracted at the lower two, and inflates
    the 6Li alpha tag from 0.0171 to 0.0257 at 5 x 41 and from 0.0163 to
    0.0241 at 10 x 100 (a factor 8 on the RP-inner share alone)."""
    channel = sp.LI6_ALPHA_TAG if channel is None else channel
    k = sp.spectator_lab_kinematics(channel, p_per_nucleon, n,
                                    beta=beta, rng=np.random.default_rng(seed))
    out = []
    for t in np.atleast_1d(theta_x):
        cut = ff.Optics("cut", t / 10.0, 10.0, theta_y / 10.0)
        acc = ff.acceptance_summary(k["R"], k["theta"], k["pT"], cut,
                                    phi=k["phi"], pot_config=pot_config)
        out.append(1.0 - acc["lost"])
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slope-b", type=float, default=50.0,
                    help="coherent t-slope [GeV^-2]")
    ap.add_argument("--n-spectator", type=int, default=200000)
    ap.add_argument("--isotope", default="6Li", choices=("6Li", "7Li"),
                    help="species of the alpha-tag panel (b).  Panel (a), "
                         "the coherent intact-nucleus recoil, is 6Li at "
                         "either setting -- see the module docstring")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    th = np.logspace(np.log10(0.05e-3), np.log10(3.0e-3), 300)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.5))

    lines = ["# near-beam aperture scan (plans/09), Yellow Report divergences "
             "(plans/10) and the tagging optics (Report 1 Section 6.1)",
             "# R12 = %s m at 5 x 41 / 10 x 100 / 18 x 275, so 1 mrad is "
             "that many mm at the pots -- measured per configuration on "
             "2026-08-28 (farforward.POT_LEVERS, plans/09 B1), where one "
             "18 x 275 number, 30.6 m, used to serve all three"
             % " / ".join("%.1f" % v[0] for v in ff.POT_LEVERS.values()),
             "# the vertical half-width is the larger of the measured "
             "aperture and the 10 sigma_v envelope of the high-acceptance "
             "optics, per configuration"]
    iso = args.isotope
    channel = ALPHA_CHANNEL[iso]
    lines.append("# panel (a) coherent: 6Li always; panel (b) alpha tag: %s"
                 % iso)
    # the machine configuration is shared; the BEAM in it is not (7Li's top
    # per-nucleon momentum is 117.9 GeV/u against 6Li's 137.5), so the
    # coherent panel keeps the 6Li configuration and the alpha panel takes
    # the isotope's own.  The measured aperture is a property of the ring
    # and is the same for both (reco.rp_aperture_for, keyed on the
    # configuration since 2026-08-28).
    for cfg, cfg_a, name, key, col in zip(beams.default_configs("6Li"),
                                          beams.default_configs(iso),
                                          ("5 x 41", "10 x 100", "18 x 275"),
                                          ("5x41", "10x100", "18x275"),
                                          (C_TRUTH, C_FIT, C_ALT)):
        pu = cfg.ion_momentum_per_nucleon
        pu_a = cfg_a.ion_momentum_per_nucleon
        p_ion = cfg.ion.A * pu
        meas_x, meas_y = reco.rp_aperture_for(cfg)
        yr = ff.yr_optics(cfg, "high-acceptance")
        env_x, env_y = yr.envelope
        top = ff.tagging_optics_point(cfg, slope_b=args.slope_b)
        tag_x = top["env_x"]
        ty = max(meas_y, env_y)
        # panel (b) rides the isotope's own optics
        env_a_x, env_a_y = ff.yr_optics(cfg_a, "high-acceptance").envelope
        top_a = ff.tagging_optics_point(cfg_a, slope_b=args.slope_b)
        ty_a = max(meas_y, env_a_y)

        acc = coherent_acceptance(th, p_ion, args.slope_b, ty)
        ax1.plot(1e3 * th, acc, "-", color=col, lw=1.6,
                 label="%s (⁶Li %g GeV/u), vertical %.1f mrad%s"
                   % (name, pu, 1e3 * ty, " (shut)" if ty > 5e-3 else ""))
        pts = {}
        for ckey, tx, mk, fill in (("silicon", meas_x, "o", "white"),
                                   ("YR envelope", env_x, "s", col),
                                   ("tagging", tag_x, "^", col)):
            pts[ckey] = coherent_acceptance(tx, p_ion, args.slope_b, ty)[0]
            ax1.plot([1e3 * tx], [pts[ckey]], mk, color=col, ms=7, mfc=fill, mew=1.6)
        lines.append("coherent %-9s silicon %.2f mrad -> %.3e ; YR HA envelope %.2f mrad "
                     "-> %.3e (gain x%.3g) ; tagging optics %.2f mrad -> %.3e "
                     "(x%.3g over silicon at L/L_HA = 1/%.1f)"
                     % (name, 1e3 * meas_x, pts["silicon"], 1e3 * env_x,
                        pts["YR envelope"], pts["YR envelope"] / max(pts["silicon"], 1e-300),
                        1e3 * tag_x, pts["tagging"],
                        pts["tagging"] / max(pts["silicon"], 1e-300),
                        1.0 / top["lumi_fraction"]))

        acc_a = alpha_acceptance(th, ty_a, pu_a, channel, n=args.n_spectator,
                                 pot_config=key)
        ax2.plot(1e3 * th, acc_a, "-", color=col, lw=1.6, label="%s" % name)
        apts = {}
        for akey, tx, mk, fill in (("silicon", meas_x, "o", "white"),
                                   ("YR envelope", env_a_x, "s", col),
                                   ("tagging", top_a["env_x"], "^", col)):
            apts[akey] = alpha_acceptance(tx, ty_a, pu_a, channel,
                                          n=args.n_spectator,
                                          pot_config=key)[0]
            ax2.plot([1e3 * tx], [apts[akey]], mk, color=col, ms=7, mfc=fill, mew=1.6)
        lines.append("alpha %s %-9s silicon -> %.4f ; YR HA envelope -> %.4f (x%.3g) ; "
                     "tagging optics -> %.4f (x%.3g at L/L_HA = 1/%.1f)"
                     % (iso, name, apts["silicon"], apts["YR envelope"],
                        apts["YR envelope"] / max(apts["silicon"], 1e-12),
                        apts["tagging"], apts["tagging"] / max(apts["silicon"], 1e-12),
                        1.0 / top_a["lumi_fraction"]))

    for ax, ylab, title in (
            (ax1, "tagged fraction of the coherent recoil",
             "(a) coherent intact ⁶Li, exp(−B|t|), B = %g GeV$^{-2}$" % args.slope_b),
            (ax2, "%s α-tag, any far-forward system (routed)" % _SUP[iso],
             "(b) %s α spectator, cluster model (β = 0.30)" % _SUP[iso])):
        ax.set_xscale("log")
        # 7Li's alpha panel spans 0.96-0.99, not fourteen decades: a log
        # axis would render the whole B3 result as one flat line
        if not (ax is ax2 and iso == "7Li"):
            ax.set_yscale("log")
        ax.set_xlabel("horizontal half-width of the near-beam cutout at the IP [mrad]")
        ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=9.5)
        ax.grid(alpha=0.25, lw=0.5)
        ax.plot([], [], "o", color="0.3", mfc="white", mew=1.6, label="silicon aperture (measured 2026-08-28)")
        ax.plot([], [], "s", color="0.3", label="10σ envelope, YR high-acceptance optics")
        ax.plot([], [], "^", color="0.3", label="10σ envelope, tagging optics (pots follow)")
        ax.legend(fontsize=7.0, loc="lower left")
    ax1.set_ylim(1e-18, 2.0)
    ax2.set_ylim((1e-3, 1.5) if iso == "6Li" else (0.94, 1.005))
    fig.suptitle("What the near-beam aperture is worth, per configuration: the tag is an "
                 "angle, and the three half-widths that compete for it"
                 + ("" if iso == "6Li" else
                    "\n(b) is ⁷Li: off rigidity at R = 0.856, so the α is "
                    "inside the momentum window and the aperture never "
                    "enters"), fontsize=9.5)
    fig.tight_layout()
    out = outdir / ("nearbeam_aperture_%s.png" % iso)
    fig.savefig(out, dpi=140)
    lines.append("wrote %s" % out)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
