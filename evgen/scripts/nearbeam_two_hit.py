#!/usr/bin/env python3
"""Two-hit topology of 6Li -> alpha + d in the Roman Pots (plans/09 B4).

The background the coherent cos 2phi measurement exists to reject is
6Li -> alpha + d: both fragments sit within 0.5% of beam rigidity
(alpha R = 0.99813, d R = 1.00452), so neither is separated from the beam
by dispersion and either can fake an intact 6Li in the near-beam tail.
Charge identification is one answer (Report 4 SS4); TOPOLOGY is another,
free one -- an intact 6Li is ONE hit and a breakup is TWO.

What this script measures, and what no module could until the joint
sampler of `spectator.breakup_lab_kinematics` existed: the two fragments
of ONE breakup, from one relative momentum k, boosted with +k and -k.
Both take transverse kicks of equal size and opposite sign, but the alpha
is carried by 4 p_u and the deuteron by 2 p_u, so

    theta_d / theta_alpha -> m_alpha / m_d = 1.987   at   phi_alpha + pi

in the k -> 0 limit (the ratio is not fixed away from it -- it runs from
0.70 to 5.3 over the sampled k -- but the deuteron is the wider fragment
in all but 3e-5 of breakups, so it is the easier one to catch and the
alpha is the binding one).  Three questions follow:

  1. the 2 x 2 topology per breakup (both fragments in acceptance, one,
     neither) at each optics;
  2. how far apart the two hits land, in millimetres and in 500 um
     pixels, and how often they merge;
  3. the number that matters for the coherent measurement: GIVEN that the
     alpha faked a coherent tag (near-beam route), how often is the
     partner deuteron there to veto it -- and how that collapses as the
     pot's OUTER edge is brought in from the 5 mrad the acceptance tables
     assume to the ~16 mm of a single sensor module.

Two sample sizes, because the two halves of the study need different
statistics.  The separation quantiles and the figure come from --events
breakups held in memory; every COUNT -- the topology fractions, the fake
rate, the conditional veto, the merge fraction -- comes from a chunked
pass over --veto-events, because at the Yellow Report optics the alpha
fakes a coherent tag in only 1e-4 to 2e-3 of breakups and a 4e5 sample
conditions on a few dozen events.  (It did: the 0.05 published for the
10 x 100 veto on 2026-08-28 rested on four successes out of 77.)

Millimetres carry one caveat since 2026-08-29 where they carried two.
(R12, R34, D) are measured PER CONFIGURATION -- 19.2 / 4.56 / 0.31,
21.2 / 3.35 / 0.29 and 30.0 / 2.93 / 0.30 m (plans/09 B1,
`farforward.POT_LEVERS`) -- in place of the single 18 x 275 pair
R12 = 30.6, D = 0.30 m with R34 taken equal to R12, which made the two
lower configurations 45-60% too wide.  The 5 x 41 vertical lever was the
last to arrive: the 29.6 mm insertion shuts that plane, so it had to be
read off a geometry with the pots slid onto the axis, and until it was
this script fell back on R12 = 19.24 m there and printed a 5 x 41 median
of 25.8 mm against the 17.3 it prints now.
And the millimetres inherit the
cluster density's short-range scale beta, which is 0.30 GeV here and
uncertain to 0.20-0.40: that band moves the medians by -12% to +9%
-- but not the veto fractions, which move by less than 0.03.  Nothing in the
ACCEPTANCE depends on any of it: the routing is done in angle.

Usage:  python3 scripts/nearbeam_two_hit.py [--events N] [--veto-events N]
                                            [--seed S] [--beta B]
                                            [--outdir .]
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

from polli_fastsim import beams  # noqa: E402
from polli_fastsim import farforward as ff  # noqa: E402
from polli_fastsim import spectator as sp  # noqa: E402

C_TRUTH, C_FIT, C_ALT, C_GREY = "#1F4E79", "#C0392B", "#B8860B", "#8A8A8A"

PIXEL_PITCH_MM = 0.5      # AC-LGAD readout pitch assumed in Report 4
MODULE_MM = 16.0          # ePIC Roman-Pot sensor module, 16 x 16 mm
OUTER_SCAN = (5.0e-3, 2.0e-3, 1.0e-3, 0.5e-3)   # pot outer edge [rad]
CHUNK = 1_000_000         # events per block of the counting pass


def optics_menu(config):
    """The four envelopes every near-beam number in this programme is
    quoted at: the two Yellow Report settings of this configuration, the
    lithium tagging optics of Report 1 SS6.1, and the single proton-derived
    73 urad the fast simulation used at every configuration before
    2026-08-28 (kept as the reference the published tables were made at)."""
    return (("YR high-acceptance", ff.yr_optics(config, "high-acceptance")),
            ("YR high-divergence", ff.yr_optics(config, "high-divergence")),
            ("tagging optics", ff.tagging_optics(config)),
            ("legacy 73 urad", ff.HIGH_ACCEPTANCE))


def route(frag, optics, theta_outer=None, pot_config="18x275"):
    """`pot_config` is the machine configuration whose blind block the
    over-rigid branch tests against (48 / 32 / 16 mm); it must be passed
    per configuration, since 18 x 275 is the most permissive of the
    three."""
    return ff.route_charged(frag["R"], frag["theta"], frag["pT"], optics,
                            phi=frag["phi"], theta_outer=theta_outer,
                            pot_config=pot_config)


def seen(frag, optics, theta_outer=None, pot_config="18x275"):
    """Is the fragment in a Roman-Pot acceptance -- the main window below
    R = 0.95 or the near-beam tail outside the envelope?"""
    r = route(frag, optics, theta_outer, pot_config)
    return (r == 1) | (r == 4)


def fakes_coherent_tag(frag, optics, theta_outer=None, pot_config="18x275"):
    """Is the fragment a NEAR-BEAM tag, i.e. does it look like the intact
    6Li recoil the coherent channel counts?  Route 4 alone: the main
    window below R = 0.95 is off rigidity and is not a coherent fake."""
    return route(frag, optics, theta_outer, pot_config) == 4


def wilson(k, n):
    """Binomial uncertainty on k/n, as the plain sqrt(p(1-p)/n) where that
    is defined and 1/n where k = 0 (so a zero column carries its own
    resolution rather than a spurious 0.000)."""
    if n == 0:
        return float("nan")
    p = k / n
    return max((p * (1.0 - p) / n) ** 0.5, 1.0 / n)


def _merge_scale(config, optics, a):
    """The scale a recorded pair's angular separation sits at: three times
    the alpha's own pot-plane displacement at the envelope, on the lever of
    the axis it ESCAPED through.  Both axes were the same number while R34
    was taken equal to R12; they are 4 to 10x apart now, and an axis no
    recorded alpha clears must not enter the minimum."""
    r12, r34, _d = ff.pot_levers_for(config)
    lev = ([r12 * optics.envelope[0]] if a["esc_x"] else []
           ) + ([r34 * optics.envelope[1]] if a["esc_y"] else [])
    return "%.0f mm" % (3e3 * min(lev)) if lev else "no pair"


def counting_pass(config, n_total, seed, beta):
    """Chunked pass returning, per optics label, the counts every fraction
    in the table is built from: the 2 x 2 topology, the near-beam fake rate
    of the alpha, the partner-seen count given that fake, the number of
    recorded pairs and how many of those land inside one pixel; plus the
    veto against each outer bound at the tagging optics."""
    menu = optics_menu(config)
    key = ff.yr_config_key(config)
    acc = {label: dict(n=0, both=0, d_only=0, a_only=0, fake=0, veto=0,
                       rec=0, merge=0, min_sep=np.inf, min_ang=np.inf,
                       esc_x=0, esc_y=0)
           for label, _ in menu}
    scan = {t: [0, 0] for t in OUTER_SCAN}
    tag = ff.tagging_optics(config)
    rng = np.random.default_rng(seed)
    done = 0
    while done < n_total:
        n = min(CHUNK, n_total - done)
        done += n
        ev = sp.breakup_lab_kinematics(sp.LI6_ALPHA_TAG,
                                       config.ion_momentum_per_nucleon, n,
                                       beta=beta, rng=rng)
        alpha, deut = ev["spectator"], ev["partner"]
        r12, r34, _d = ff.pot_levers_for(config)
        sep_mm = 1e3 * ff.separation_at_pots(alpha, deut, config=config)
        ang_mm = 1e3 * ff.separation_at_pots(alpha, deut, r12=r12, r34=r34,
                                             dispersion=0.0)
        for label, opt in menu:
            ha = seen(alpha, opt, pot_config=key)
            hd = seen(deut, opt, pot_config=key)
            fake = fakes_coherent_tag(alpha, opt, pot_config=key)
            a = acc[label]
            a["n"] += n
            a["both"] += int(np.sum(ha & hd))
            a["d_only"] += int(np.sum(~ha & hd))
            a["a_only"] += int(np.sum(ha & ~hd))
            a["fake"] += int(np.sum(fake))
            a["veto"] += int(np.sum(fake & hd))
            rec = ha & hd
            # which axis a RECORDED alpha escaped through: the merge scale
            # below needs the lever of that axis, and since 2026-08-29 the
            # two levers differ by 4 to 10x at every configuration
            cx, cy = opt.envelope
            a["esc_x"] += int(np.sum(rec & (np.abs(
                alpha["theta"] * np.cos(alpha["phi"])) > cx)))
            a["esc_y"] += int(np.sum(rec & (np.abs(
                alpha["theta"] * np.sin(alpha["phi"])) > cy)))
            pair = sep_mm[rec]
            a["rec"] += int(pair.size)
            a["merge"] += int(np.sum(pair < PIXEL_PITCH_MM))
            if pair.size:
                a["min_sep"] = min(a["min_sep"], float(pair.min()))
                a["min_ang"] = min(a["min_ang"], float(ang_mm[rec].min()))
        for th in OUTER_SCAN:
            f = fakes_coherent_tag(alpha, tag, theta_outer=th, pot_config=key)
            scan[th][0] += int(np.sum(f))
            scan[th][1] += int(np.sum(f & seen(deut, tag, theta_outer=th,
                                               pot_config=key)))
    return acc, scan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", type=int, default=400000,
                    help="breakups held in memory for the quantiles/figure")
    ap.add_argument("--veto-events", type=int, default=12000000,
                    help="breakups in the chunked counting pass")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--beta", type=float, default=0.30,
                    help="short-range scale of the cluster density [GeV]")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    lines = ["# two-hit topology of 6Li -> alpha + d (plans/09 B4), seed %d, "
             "beta = %.2f GeV" % (args.seed, args.beta),
             "# quantiles from %d breakups per configuration; every COUNT "
             "(topology, fake rate, veto, merge) from %d, chunked -- the "
             "Yellow Report fake rate is 1e-4 to 2e-3 and a small sample "
             "conditions on a handful of events"
             % (args.events, args.veto_events),
             "# both fragments from ONE relative momentum k "
             "(spectator.breakup_lab_kinematics); routed with their azimuth",
             "# millimetres use the PER-CONFIGURATION levers measured on "
             "2026-08-28 (farforward.POT_LEVERS, tools/fullsim, plans/09 "
             "B1): (R12, R34, D) = " + " ; ".join(
                 "%s %.1f / %s / %.2f m"
                 % (k, v[0], "--" if v[1] is None else "%.2f" % v[1], v[2])
                 for k, v in ff.POT_LEVERS.items())
             + ".  Until then one 18 x 275 pair, R12 = 30.6 and D = 0.30 m, "
             "was applied at every configuration and R34 was taken equal to "
             "R12; the lower configurations' millimetres were 45-60% high. "
             " R34 at 5 x 41 was measured a day later, off a geometry whose "
             "Roman-pot insertions are zeroed (tools/fullsim): 4.56 m, "
             "against the R12 = 19.24 m this script fell back on until "
             "then, which cost the 5 x 41 median 8.5 mm.  "
             "They also inherit beta (0.20-0.40 moves the medians by -12% "
             "to +9%).  The acceptance itself is angular and does not "
             "depend on any of them.",
             "# the MEASURED outer edge of the pot acceptance is "
             + " / ".join("%.2f" % (1e3 * ff.theta_rp_outer_for(k))
                          for k in ("5x41", "10x100", "18x275"))
             + " mrad (5 x 41 / 10 x 100 / 18 x 275), not the 5 mrad the "
             "acceptance tables assume and not the 144 mm / R12 = 7.5 / "
             "6.8 / 4.8 mrad the module arithmetic gives: past it the ion "
             "has struck the pipe.  Read the veto scan below against it."]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.5))
    cfgs = beams.default_configs("6Li")
    names = ("5 x 41", "10 x 100", "18 x 275")
    for cfg, name, col in zip(cfgs, names, (C_TRUTH, C_FIT, C_ALT)):
        pu = cfg.ion_momentum_per_nucleon
        ev = sp.breakup_lab_kinematics(sp.LI6_ALPHA_TAG, pu, args.events,
                                       beta=args.beta,
                                       rng=np.random.default_rng(args.seed))
        alpha, deut = ev["spectator"], ev["partner"]
        tag_opt = ff.tagging_optics(cfg)
        r12, r34, disp = ff.pot_levers_for(cfg)
        sep_mm = 1e3 * ff.separation_at_pots(alpha, deut, config=cfg)
        ang_mm = 1e3 * ff.separation_at_pots(alpha, deut, r12=r12,
                                             r34=r34, dispersion=0.0)
        q16, q50, q84 = np.percentile(sep_mm, [16, 50, 84])
        a16, a50, a84 = np.percentile(ang_mm, [16, 50, 84])
        acc, scan = counting_pass(cfg, args.veto_events, args.seed, args.beta)

        lines.append("")
        lines.append("== %s (6Li %.1f GeV/u) ==  R_alpha = %.5f, R_d = %.5f "
                     "at k = 0 ; theta_d/theta_alpha -> %.3f as k -> 0"
                     % (name, pu, sp.MASSES["alpha"] / 2.0
                        / (sp.nucleus_mass(3, 6) / 3.0),
                        sp.nucleus_mass(1, 2) / (sp.nucleus_mass(3, 6) / 3.0),
                        sp.MASSES["alpha"] / sp.nucleus_mass(1, 2)))
        lines.append("   separation at the pots: 16/50/84%% = %.1f / %.1f / "
                     "%.1f mm = %.0f / %.0f / %.0f pixels of %g um"
                     % (q16, q50, q84, q16 / PIXEL_PITCH_MM,
                        q50 / PIXEL_PITCH_MM, q84 / PIXEL_PITCH_MM,
                        1e3 * PIXEL_PITCH_MM))
        d_mm = np.abs(1e3 * disp * (np.asarray(alpha["R"])
                                    - np.asarray(deut["R"])))
        lines.append("      angular lever alone (no dispersion): %.1f / %.1f "
                     "/ %.1f mm -- the D (R_a - R_d) term is %.1f / %.1f / "
                     "%.1f mm at this configuration's D = %.2f m and does "
                     "NOT cancel between the fragments"
                     % ((a16, a50, a84)
                        + tuple(np.percentile(d_mm, [16, 50, 84])) + (disp,)))
        lines.append("   %-20s %-16s %8s %8s %8s %8s | %10s %14s %9s | %8s %9s"
                     % ("optics", "envelope [mrad]", "both", "d only",
                        "alpha only", "neither", "P(a fakes)",
                        "partner seen", "n_fake", "min sep", "merged"))
        for label, opt in optics_menu(cfg):
            a = acc[label]
            n = a["n"]
            veto = ("%.3f +- %.3f" % (a["veto"] / a["fake"],
                                      wilson(a["veto"], a["fake"]))
                    if a["fake"] else "-")
            lines.append("   %-20s %6.2f x %-7.2f %8.4f %8.4f %8.4f %8.4f | "
                         "%10.5f %14s %9d | %8s %9s"
                         % (label, 1e3 * opt.envelope[0], 1e3 * opt.envelope[1],
                            a["both"] / n, a["d_only"] / n, a["a_only"] / n,
                            (n - a["both"] - a["d_only"] - a["a_only"]) / n,
                            a["fake"] / n, veto, a["fake"],
                            "%.2f mm" % a["min_sep"] if a["rec"] else "-",
                            "%.1e" % (a["merge"] / a["rec"])
                            if a["rec"] else "-"))
        lines.append("   merge = recorded pair inside one %g um pixel; RARE, "
                     "not impossible.  In the ANGULAR lever alone a recorded "
                     "pair stays about 3x the alpha's own displacement apart "
                     "-- back to back with theta_d ~ 2 theta_alpha -- and no "
                     "recorded pair could merge at all: %s (a scale, not a "
                     "bound; the measured minima sit within 10%% of it either "
                     "way, because theta_d/theta_alpha is 1.987 only as "
                     "k -> 0 and the sparser rows are sampling their own "
                     "tail).  The "
                     "scale is 3 min(R12 env_x, R34 env_y) over the axes the "
                     "recorded alphas are SEEN to escape through, which is "
                     "not the same as over both: at 5 x 41 the vertical "
                     "envelope is 3.80 mrad and no recorded alpha reaches it "
                     "at the Yellow Report or the tagging optics, so the "
                     "4.56 m vertical lever is not a scale of anything there. "
                     " The dispersive term is free to cancel the angular one, "
                     "and does.  Over ALL breakups the minimum is %.2f mm and "
                     "%.4f fall inside a pixel; most of those are collinear "
                     "pairs inside the envelope that no pot records, but not "
                     "all of them, which is the merge column."
                     % (1e3 * PIXEL_PITCH_MM,
                        " ; ".join("%s scale %s, min %s"
                                   % (label, _merge_scale(cfg, opt, acc[label]),
                                      "%.1f" % acc[label]["min_ang"]
                                      if acc[label]["rec"] else "no pair")
                                   for label, opt in optics_menu(cfg)),
                        sep_mm.min(),
                        float(np.mean(sep_mm < PIXEL_PITCH_MM))))

        # the veto, and what the pot's outer edge does to it
        key = ff.yr_config_key(cfg)
        fake = fakes_coherent_tag(alpha, tag_opt, pot_config=key)
        s = sep_mm[fake & seen(deut, tag_opt, pot_config=key)]
        lines.append("   given an alpha near-beam tag at the tagging optics: "
                     "separation 16/50/84%% = %.0f / %.0f / %.0f mm"
                     % tuple(np.percentile(s, [16, 50, 84])))
        curve = [scan[t][1] / scan[t][0] if scan[t][0] else np.nan
                 for t in OUTER_SCAN]
        lines.append("   partner-veto efficiency vs the pot outer edge: "
                     + " ".join("%.1f mrad (%.0f mm): %.2f"
                                % (1e3 * t, 1e3 * ff.pot_levers_for(cfg)[0] * t,
                                   v)
                                for t, v in zip(OUTER_SCAN, curve)))

        bins = np.logspace(-2.0, 2.7, 110)
        ax1.hist(sep_mm, bins=bins, histtype="step", lw=1.0, color=col,
                 alpha=0.45, density=True,
                 label="%s, all breakups (median %.1f mm)" % (name, q50))
        rec = sep_mm[seen(alpha, tag_opt, pot_config=key)
                     & seen(deut, tag_opt, pot_config=key)]
        ax1.hist(rec, bins=bins, histtype="stepfilled", lw=1.6, color=col,
                 alpha=0.20, density=True)
        ax1.hist(rec, bins=bins, histtype="step", lw=1.7, color=col,
                 density=True, label="        both fragments recorded")
        ax2.plot([1e3 * t for t in OUTER_SCAN], curve, "o-", color=col, lw=1.6,
                 label="%s" % name)

    ax1.axvline(PIXEL_PITCH_MM, color=C_GREY, ls=":", lw=1.4)
    ax1.text(PIXEL_PITCH_MM * 1.15, 0.60, "500 μm pixel", rotation=90,
             fontsize=7.5, color=C_GREY, va="top",
             bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.0),
             transform=ax1.get_xaxis_transform())
    ax1.axvline(MODULE_MM, color=C_GREY, ls="--", lw=1.4)
    ax1.text(MODULE_MM * 1.15, 0.60, "16 mm module", rotation=90, fontsize=7.5,
             color=C_GREY, va="top",
             bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.0),
             transform=ax1.get_xaxis_transform())
    ax1.set_xscale("log")
    ax1.set_xlabel("α–d separation at the pot plane [mm], per-configuration "
                   "R₁₂ = %s m, D = %s m"
                   % (" / ".join("%.1f" % v[0] for v in ff.POT_LEVERS.values()),
                      " / ".join("%.2f" % v[2] for v in ff.POT_LEVERS.values())))
    ax1.set_ylabel("breakups per unit log separation")
    ax1.set_title("(a) a recorded pair is two hits tens of pixels apart; the "
                  "pairs that could merge\nare mostly the collinear ones "
                  "inside the envelope, which no pot records (tagging optics)",
                  fontsize=9.0)
    ax2.set_xlabel("outer edge of the pot acceptance [mrad]")
    ax2.set_ylabel("partner d seen | α faked a coherent tag")
    ax2.set_title("(b) the veto is a station-layout question", fontsize=9.5)
    ax2.set_ylim(0.0, 1.0)
    for ax in (ax1, ax2):
        ax.grid(alpha=0.25, lw=0.5)
        ax.legend(fontsize=6.8 if ax is ax1 else 7.5,
                  loc="upper left" if ax is ax1 else "lower right",
                  framealpha=0.85)
    fig.suptitle("⁶Li → α + d in the Roman Pots: what hit multiplicity is "
                 "worth against the coherent tag (tagging optics)",
                 fontsize=9.5)
    fig.tight_layout()
    out = outdir / "nearbeam_two_hit_6Li.png"
    fig.savefig(out, dpi=140)
    lines.append("")
    lines.append("wrote %s" % out)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
