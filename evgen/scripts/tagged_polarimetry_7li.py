"""7Li alpha-tag bonus observables (plans/05 step 5.B deliverables):

Left: the in-situ alignment polarimeter -- the tagged-alpha polar
anisotropy <P2(cos theta_k)> = -T/5 vs the fill tensor moment T, with
Roman-Pot-folded pseudo-experiment markers.

Right: the tagged polarized-EMC companion -- A_par^tag(x) on the
quasi-free triton (helicity-flip fills, RP-accepted alpha tag), analytic
D(y) g1_t/F1_t overlay.  P_p(t) = 0.86 makes the tagged triton an
effective polarized proton at 86% polarization with an alpha in the RP
as the event tag.

WHY THE 7Li FIGURE BARELY MOVES, AND WHAT THAT MEANS (plans/09 B3,
2026-08-28).  Both panels used to be routed through the retired
proton-derived 73 / 164 microrad envelope, applied as a circle because
the spectator's lab azimuth was never passed to `route_charged`.  They
are now routed at the per-configuration Yellow Report optics and at the
lithium tagging optics, with the rectangular 10 (sigma_h, sigma_v)
envelope and the true azimuth (`--optics`, `--config`).  Nothing in the
physics moves, and that is the deliverable: the 7Li alpha sits at
rigidity ratio R = 0.856, far off the beam and in the MIDDLE of the
Roman-Pot momentum window, so it is accepted by the window and never
has to clear the near-beam envelope at all.  Its acceptance is therefore
optics-blind: this script's own Roman-Pot tag reads 0.962 / 0.968 / 0.973
at the Yellow Report optics against 0.980 / 0.991 / 0.992 at the tagging
optics, and the folded <P2> markers and A_par bars move by ~1%.

(Report 3 Table 6 and `fastsim/out/tagging_acceptance.txt` quote 0.9690 /
0.9683 / 0.9748 and 0.9872 / 0.9909 / 0.9941 for the same tag.  Those are
the pure spectator model, and they are `1 - lost`, i.e. ANY far-forward
system, while `accepted()` below is the Roman Pots alone.  At 5 x 40.8
the B0 carries 1.1% of the 7Li alpha, so the like-for-like pure-model
Roman-Pot numbers are 0.9568 / 0.9668 / 0.9721; restricted further to
k <= 1.2 GeV/c, where the tagged model's momentum grid ends, they are
0.9626 / 0.9683 / 0.9736 against this script's 0.9617 / 0.9678 / 0.9728
-- agreement to 0.1 point at every configuration.  The headline block
prints both definitions so the comparison can be made directly.  Those
three numbers were 0.9614 / 0.9676 / 0.9726 before 2026-08-28 and
0.9620 / 0.9678 / 0.9728 until the 5 x 41 R34 measurement of 2026-08-29: the
accepted sample is cross-section-weighted and the cross section carries
g1 of the struck triton, so making TRITON per-nucleon like every other
Ion slot -- plans/08 D7, which gives its N = 2 neutrons their second
term -- moves the tag in the fourth decimal.  The g1_t/F1_t overlay this
figure draws moves by +8.4% at x = 0.005, +5.3% at 0.01, +2.4% at 0.03,
+0.8% at 0.1, -0.5% at 0.5 and -0.7% at 0.7.)

The consequence is a programme statement, not a plot.  For 7Li the
tagging optics buys x1.02 in acceptance for x1/8 to x1/15 in luminosity:
it is a factor 8-15 NET LOSS, the exact inverse of 6Li, where the same
optics is the difference between a 1.5% tag and a 30% one.  6Li and 7Li
therefore want different machine optics and are different runs, not one
fill plan (plans/09 B3, Report 4 conclusion).

Output: `tagged_polarimetry_7Li.png` at the published configuration
(--config 1) and `tagged_polarimetry_7Li_<key>.png` at the other two.

Usage:  python3 scripts/tagged_polarimetry_7li.py --events 300000
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
from polligen import estimators as est  # noqa: E402
from polligen import tagged  # noqa: E402
from polligen.spin import spin32_populations  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.asymmetries import a_parallel, err_a_parallel  # noqa: E402
from polli_fastsim import farforward as ff  # noqa: E402
from polli_fastsim.farforward import route_charged  # noqa: E402

PE, PZ = 0.7, 0.7
C_HA, C_TAG, C_LEG, C_HD = "#1F4E79", "#C0392B", "#8A8A8A", "#B8860B"


def optics_menu(config, which):
    """(label, Optics, colour, marker) per optics -- the same menu as
    money_tagged_azz.py: 'menu' is the Yellow Report high-acceptance optics
    of this configuration plus the tagging optics with its luminosity
    fraction; 'legacy' the retired proton-derived pair."""
    ha = (ff.yr_optics(config, "high-acceptance"), C_HA, "o")
    hd = (ff.yr_optics(config, "high-divergence"), C_HD, "v")
    tag = (ff.tagging_optics(config), C_TAG, "^")
    table = {
        "menu": [ha, tag],
        "high-acceptance": [ha],
        "high-divergence": [hd],
        "tagging": [tag],
        "legacy": [(ff.HIGH_ACCEPTANCE, C_LEG, "s"),
                   (ff.HIGH_DIVERGENCE, C_HD, "D")],
    }
    return [(o.name, o, c, m) for o, c, m in table[which]]


def output_stem(base, key, config, optics):
    """As in money_tagged_azz.py: the published artefact is the default
    combination alone, so no `--optics` at the default `--config` can
    overwrite it."""
    if config == 1 and optics == "menu":
        return base
    return "%s_%s_%s" % (base, key, optics)


def accepted(ev, optics, pot_config):
    """RP mask (main window + near-beam tail) with the SPECTATOR's lab
    azimuth, so the rectangular envelope is applied as a rectangle.

    `pot_config` is the machine configuration whose blind block the
    over-rigid branch tests against.  Routes 1 and 4 cannot themselves
    move with it -- RP_R_WINDOW ends at R = 0.95 and the near-beam tail
    at |R - 1| < 0.05, while the over-rigid branch starts at R > 1.05,
    so the three are disjoint -- but it is passed because the call is
    meant to be the configuration's routing and not the default one, and
    because `1 - lost` in the headline block below DOES move with it."""
    route = route_charged(ev["R"], ev["theta"], ev["pT"], optics,
                          phi=ev["phi_spec"], pot_config=pot_config)
    return (route == 1) | (route == 4)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", type=float, default=3e5)
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2),
                    help="beam configuration index (5 x 40.8, 10 x 99.5, "
                         "18 x 117.9 GeV/u for 7Li)")
    ap.add_argument("--optics", default="menu",
                    choices=("menu", "legacy", "high-acceptance",
                             "high-divergence", "tagging"),
                    help="which envelope(s) to route the alpha through; "
                         "'menu' (default) is the Yellow Report "
                         "high-acceptance optics of this configuration and "
                         "the tagging optics, 'legacy' the retired "
                         "proton-derived 73/164 microrad pair")
    ap.add_argument("--seed", type=int, default=20260713)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs("7Li")[args.config]
    key = ff.yr_config_key(config)
    menu = optics_menu(config, args.optics)
    model = tagged.TaggedModel(tagged.li7_alpha_channel())
    kern = InclusiveKernel(tagged.TRITON)
    sampler = tagged.TaggedSampler(model, kern, config, fom.Scenario())
    rng = np.random.default_rng(args.seed)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # --- left: alignment polarimetry ----------------------------------------
    # one sample per T, re-routed at each optics: the optics comparison is
    # then paired and any difference between the markers is acceptance,
    # not Monte-Carlo noise.
    t_scan = np.linspace(-1.0, 1.0, 9)
    ax1.plot(t_scan, -t_scan / 5.0, "k-", lw=1.5,
             label=r"analytic $\langle P_2\rangle = -T/5$")
    n_per_point = int(args.events / t_scan.size)
    p2 = {name: [] for name, _o, _c, _m in menu}
    for t in t_scan:
        pops = spin32_populations(0.0, t, 0.0)
        cat = bk.SpinCategory("pol T=%.2f" % t, 1.5, pops)
        ev = sampler.sample_category(cat, n=n_per_point, rng=rng)
        for name, optics, _c, _m in menu:
            c = ev["cos_theta_k"][accepted(ev, optics, key)]
            p2[name].append(float(np.mean(0.5 * (3 * c * c - 1.0))))
    slopes = {}
    for name, optics, colour, marker in menu:
        vals = np.array(p2[name])
        slopes[name] = float(np.polyfit(t_scan, vals, 1)[0])
        ax1.plot(t_scan, vals, marker, ms=4, color=colour,
                 label="RP-folded, %s" % name)
    ax1.set_xlabel(r"fill tensor moment $T$")
    ax1.set_ylabel(r"$\langle P_2(\cos\theta_k)\rangle$ of the tagged "
                   r"$\alpha$")
    ax1.set_title("in-situ alignment polarimeter")
    ax1.legend(fontsize=8)
    ax1.axhline(0, color="0.6", lw=0.5)

    # --- right: tagged A_par on the quasi-free triton -----------------------
    plan = bk.helicity_flip_plan(1.5, PZ, PE)
    x_edges = np.logspace(-3, np.log10(0.7), 11)
    xc = np.sqrt(x_edges[:-1] * x_edges[1:])
    evs = {cat.name: sampler.sample_category(
        cat, n=int(args.events * cat.lumi_fraction), rng=rng)
        for cat in plan.categories}
    n_gen = sum(len(ev["x"]) for ev in evs.values())
    apar = {}
    for i, (name, optics, colour, marker) in enumerate(menu):
        counts = {cname: np.histogram(ev["x"][accepted(ev, optics, key)],
                                      bins=x_edges)[0]
                  for cname, ev in evs.items()}
        n_tot = counts["apar+"] + counts["apar-"]
        with np.errstate(divide="ignore", invalid="ignore"):
            a = np.where(n_tot > 0, est.apar_flip(
                counts["apar+"], counts["apar-"], PE, PZ), np.nan)
        errs = np.where(n_tot > 1,
                        err_a_parallel(np.maximum(n_tot, 1), PE, PZ), np.nan)
        use = n_tot > 100
        apar[name] = (a, errs, n_tot)
        # a small multiplicative x offset per optics: the two series would
        # otherwise sit exactly on top of each other, which IS the result
        ax2.errorbar(xc[use] * (1.0 + 0.04 * i), a[use], yerr=errs[use],
                     fmt=marker, ms=4, capsize=2, color=colour,
                     mfc="none" if i else colour,
                     label="tagged pseudo-exp, %s" % name)
    # Analytic overlay: the sigma-weighted A_par of the inner grid, per
    # x bin.  It is the FINITE-GAMMA A_par (default since 2026-08-29), the
    # same one the sampler above drew the pseudo-experiment from, so the
    # curve and the points test the generator rather than the difference
    # between two conventions.  Taking the massless D(y) g1/F1 here
    # instead -- what this line did before -- would put the curve below
    # the points by the target-mass term, up to 2.1% in the top x bin at
    # this configuration and 5.0% at 5 x 40.8
    # (`evgen/scripts/target_mass_bound.py` block 4).
    inner = sampler.inner
    t = inner.tables
    y = inner.q2_cells / (inner.s * inner.x_cells)
    apar_cells = a_parallel(t["g1"], t["f1"], y, inner.x_cells,
                            inner.q2_cells, g2=t.get("g2"))
    xb = np.digitize(inner.x_cells, x_edges) - 1
    analytic = np.array([
        np.average(apar_cells[xb == b], weights=inner.xsec_flat[xb == b])
        if np.any(xb == b) else np.nan for b in range(xc.size)])
    ax2.plot(xc, analytic, "k-", lw=1.5,
             label=r"$D_\gamma(A_1+\eta A_2)$ ($\sigma$-weighted)")
    ax2.set_xscale("log")
    ax2.set_xlabel("x")
    ax2.set_ylabel(r"$A_\parallel^{\rm tag}$ (quasi-free $t$)")
    ax2.set_title(r"tagged polarized-EMC companion, $P_p(t)=0.86$")
    ax2.legend(fontsize=8)
    ax2.axhline(0, color="0.6", lw=0.5)

    # --- headline numbers, computed before the figure is titled -------------
    # `acc` is the Roman-Pot tag this figure folds with (route 1 | 4);
    # `acc_ff` is 1 - lost, ANY far-forward system, which is the definition
    # Report 3 Table 6 and fastsim/out/tagging_acceptance.txt tabulate.
    # The two differ by the B0, which takes 1.1% of the alpha at 5 x 40.8
    # and nothing at the other two, and -- since the per-configuration
    # blind block reached this script on 2026-08-28 -- by the over-rigid
    # RP-inner branch, which the high-k tail of the same distribution
    # enters at a few times 1e-3.  Hence 0.9699 / 0.9690 / 0.9751 against
    # 0.9620 / 0.9678 / 0.9728 at the Yellow Report optics.
    cat = bk.SpinCategory("acc", 1.5, (0.25, 0.25, 0.25, 0.25))
    ev_acc = sampler.sample_category(cat, n=100_000, rng=rng)
    head, ref = [], None
    for name, optics, _c, _m in menu:
        route = route_charged(ev_acc["R"], ev_acc["theta"], ev_acc["pT"],
                              optics, phi=ev_acc["phi_spec"],
                              pot_config=key)
        acc = float(np.mean((route == 1) | (route == 4)))
        acc_ff = float(np.mean(route != 0))
        fom_ = acc * optics.lumi_fraction
        ref = fom_ if ref is None else ref
        head.append((name, optics, acc, acc_ff, fom_,
                     float(np.sqrt(ref / fom_))))

    penalty = "x%s" % " / x".join("%.1f" % h[-1] for h in head[1:]) \
        if len(head) > 1 else "none"
    fig.suptitle(r"$^7$Li($e,e^\prime\alpha$)X, %s (TOY inputs): the "
                 r"$\alpha$ tag is optics-blind -- both observables are "
                 r"free at the published optics" % config.label()
                 + "\n" + "error bars at equal GENERATED statistics; at "
                 "equal luminosity the tagging optics costs %s in every "
                 "bar, which is why it is a net loss for $^7$Li" % penalty,
                 fontsize=9.5)
    fig.tight_layout()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = output_stem("tagged_polarimetry_7Li", key, args.config,
                       args.optics)
    out = outdir / (stem + ".png")
    fig.savefig(out, dpi=140)
    print("wrote", out)

    # --- headline numbers ---------------------------------------------------
    print("7Li alpha-tag, %s, %d events in the A_par panel:"
          % (config.label(), n_gen))
    for name, optics, acc, acc_ff, fom_, mult in head:
        a, errs, n_tot = apar[name]
        use = n_tot > 100
        print("  %-30s acc(RP) %.4f  acc(any far-fwd) %.4f  L/L_HA %.3f  "
              "acc x L %.4f  -> error bars x%.2f at equal luminosity; "
              "<P2> slope %+.4f (analytic -0.2000); median dA_par %.5f "
              "(at equal generated statistics)"
              % (name, acc, acc_ff, optics.lumi_fraction, fom_, mult,
                 slopes[name], float(np.median(errs[use]))))


if __name__ == "__main__":
    main()
