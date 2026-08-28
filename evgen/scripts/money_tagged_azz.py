"""Money plot 4: tagged tensor asymmetry A_zz^tag for the alpha-tagged
embedded deuteron in 6Li -- the first tagged spin observable projected for
any A > 2 (plans/05 step 5.B; extends the four-gap list of plans/00).

Left: the analytic wave-function asymmetry A_zz^wf(k, theta_k = 90 deg)
from the alpha-d S/D interference, with the model band (beta and P_D
scans -- the same tail systematics that dominate the 6Li RP acceptance).
Right: acceptance-folded pseudo-experiment -- tensor-thirds fills, events
routed through the far-forward windows, Roman-Pot-accepted (main window +
pT tail) thirds estimator vs the spectator momentum k, at each optics of
the menu.  The plan's question -- does the Cosyn-Weiss-style O(1)
asymmetry at p_s ~ 300 MeV/c survive the RP pT-tail acceptance that
dominates 6Li? -- is answered by the right panel: yes, and at the
published optics ONLY there.

Two corrections of 2026-08-28 (plans/09 B2) change what this figure says.

1.  THE OPTICS ARE PER CONFIGURATION.  Until this round the routing used
    the retired proton-derived pair (73 / 164 microrad, isotropic, the
    same at every energy) and applied it as a CIRCLE, because the
    spectator's lab azimuth was never passed to `route_charged`.  The
    envelope is a rectangle of half-widths 10 (sigma_h, sigma_v) that
    depends on the beam configuration (`farforward.yr_optics`,
    `tagging_optics`), and the azimuth decides which half-width a
    fragment has to clear: ignoring it overstates the tag by 1.7x at the
    tagging optics.  `--optics` selects the menu.

2.  THE TRUTH OVERLAY IS ACCEPTANCE-WEIGHTED.  The analytic curve is
    A_zz^wf at theta_k = 90 deg, and the ACCEPTED sample is not at 90
    deg: the near-beam tail is transverse (<|cos theta_k|> = 0.40 at the
    tagging optics) while the off-rigidity R < 0.95 window slice, which
    is all that survives at the Yellow Report optics, is longitudinal
    (0.71-0.79).  The two read the S/D interference at opposite ends of
    theta_k, and the +-0.5 swing the published version showed between its
    two optics at k ~ 0.3 GeV/c was that, not the wave function.  The
    right panel therefore carries, per optics, the acceptance-weighted
    prediction

        A_zz^acc(k) = sum_c eps(k,c) [n_+1 + n_-1 - 2 n_0]
                    / sum_c eps(k,c) [n_+1 + n_-1 +   n_0]

    with eps(k, c) the routed acceptance of the model grid, which is what
    the markers should be compared against; the 90 deg curve is kept as a
    labelled reference.  eps concentrated at c = 0 returns it exactly
    (pinned in test_tagged.py).  Closure, measured with `--events 8e6
    --config 1` at the k = 0.325 GeV/c bin: +0.4471 +- 0.0129 against the
    weighted prediction +0.4548 at the Yellow Report optics (0.6 sigma)
    and -0.1019 +- 0.0069 against -0.0984 at the tagging optics (0.5
    sigma).  Over all populated bins the residual is |dA_zz| <= 0.018 at
    both optics -- within 1.4 sigma (YR) and 2.2 sigma (tagging) at the
    default 4e5 events, where the errors are 0.03-0.09.  At 8e6 the errors
    fall to 0.003-0.03 and one bin, k = 0.175 at the tagging optics, then
    reads 3.3 sigma: that bin is where eps(k) turns on (nothing at all is
    accepted below k = 0.189 GeV/c at the Yellow Report optics), so the
    truth at the bin CENTRE is not the truth the bin-averaged marker
    measures.  Averaged over each bin the way the marker is, every
    populated bin is within 1.4 sigma at 8e6 as well.

Output: `money_tagged_azz_6Li.png` for the published combination
(--config 1 --optics menu) and `money_tagged_azz_6Li_<key>_<optics>.png`
otherwise, so no exploratory run can overwrite the published artefact.

Usage:  python3 scripts/money_tagged_azz.py --events 400000
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
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.asymmetries import err_azz  # noqa: E402
from polli_fastsim import farforward as ff  # noqa: E402
from polli_fastsim.farforward import route_charged  # noqa: E402
from polli_fastsim.polarized import toy_b1  # noqa: E402

PZ, PZZ = 0.7, 0.6
C_HA, C_TAG, C_LEG, C_HD = "#1F4E79", "#C0392B", "#8A8A8A", "#B8860B"


def optics_menu(config, which):
    """The optics to route at, as (label, Optics, colour, marker) tuples.

    'menu' is the pair every per-configuration figure of this round shows:
    the Yellow Report high-acceptance optics of THIS configuration and the
    lithium tagging optics with its luminosity fraction.  'legacy' is the
    retired proton-derived 73 microrad, kept so a pre-2026-08-28 number
    can be reproduced."""
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
    """File stem for one run.  The published artefact is the DEFAULT
    combination alone -- `--config 1 --optics menu`.  Keying the guard on
    the configuration alone would let `--optics legacy` at the default
    configuration overwrite the published PNG with the retired 73/164
    microrad figure, which is exactly the run the manual documents."""
    if config == 1 and optics == "menu":
        return base
    return "%s_%s_%s" % (base, key, optics)


# the two live in polligen.tagged so a test can pin them (and so Report 4
# can quote the weighted curve without importing a script)
azz_wf_curve = tagged.azz_tensor_curve
acceptance_table = tagged.acceptance_weights


def folded_asymmetry(sampler, plan, n_events, k_edges, optics_list, rng):
    """RP-accepted thirds estimator vs k at each optics.

    The sample is drawn ONCE and re-routed per optics, so the optics
    curves are paired and their difference is not statistical."""
    evs = {cat.name: sampler.sample_category(
        cat, n=int(n_events * cat.lumi_fraction), rng=rng)
        for cat in plan.categories}
    out = {}
    for name, optics, _c, _m in optics_list:
        counts, kacc = {}, []
        for cname, ev in evs.items():
            route = route_charged(ev["R"], ev["theta"], ev["pT"], optics,
                                  phi=ev["phi_spec"])
            acc = (route == 1) | (route == 4)
            counts[cname] = np.histogram(ev["k"][acc], bins=k_edges)[0]
            kacc.append(ev["k"][acc])
        n_tot = sum(counts.values())
        with np.errstate(divide="ignore", invalid="ignore"):
            a = np.where(n_tot > 0, est.azz_thirds(
                counts["azz+"], counts["azz-"], counts["azz0"], PZZ), np.nan)
        out[name] = (a, n_tot, np.concatenate(kacc))
    n_gen = sum(len(ev["k"]) for ev in evs.values())
    return out, n_gen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", type=float, default=4e5,
                    help="generated events per pseudo-experiment")
    ap.add_argument("--lumi-ref", type=float, default=10.0,
                    help="reference lumi [fb^-1/u] for projected errors, "
                         "AT THE HIGH-ACCEPTANCE OPTICS: each optics' own "
                         "lumi_fraction multiplies it")
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2),
                    help="beam configuration index (5 x 40.8, 10 x 99.5, "
                         "18 x 137.5 GeV/u)")
    ap.add_argument("--optics", default="menu",
                    choices=("menu", "legacy", "high-acceptance",
                             "high-divergence", "tagging"),
                    help="which envelope(s) to route the alpha through. "
                         "'menu' (default) is the published pair -- the "
                         "Yellow Report high-acceptance optics of this "
                         "configuration and the tagging optics of Report 1 "
                         "Section 6.1 with its luminosity fraction; "
                         "'legacy' the retired proton-derived 73/164 "
                         "microrad pair, for reproduction only")
    ap.add_argument("--seed", type=int, default=20260713)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    config = beams.default_configs("6Li")[args.config]
    key = ff.yr_config_key(config)
    menu = optics_menu(config, args.optics)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # --- left: analytic band ------------------------------------------------
    central = tagged.TaggedModel(tagged.li6_alpha_channel())
    ic = np.argmin(np.abs(central.c))
    ax1.plot(central.k, azz_wf_curve(central, ic), "k-", lw=2,
             label=r"$\beta=0.3$, $P_D=%.3f$" % tagged.P_D_LI6)
    for beta in (0.20, 0.40):
        m = tagged.TaggedModel(tagged.li6_alpha_channel(beta=beta))
        ax1.plot(m.k, azz_wf_curve(m, np.argmin(np.abs(m.c))), "--", lw=1,
                 label=r"$\beta=%.2f$" % beta)
    for p_d in (0.04, 0.13):
        m = tagged.TaggedModel(tagged.li6_alpha_channel(p_d=p_d))
        ax1.plot(m.k, azz_wf_curve(m, np.argmin(np.abs(m.c))), ":", lw=1,
                 label=r"$P_D=%.2f$" % p_d)
    ax1.set_xlabel(r"$k$ [GeV/$c$]")
    ax1.set_ylabel(r"$A_{zz}^{\rm tag}(k,\ \theta_k=90^\circ)$")
    ax1.set_title(r"analytic: $\alpha$-$d$ S/D interference")
    ax1.set_xlim(0, 0.6)
    ax1.legend(fontsize=7)
    ax1.axhline(0, color="0.6", lw=0.5)

    # --- right: acceptance-folded pseudo-experiment -------------------------
    kern = InclusiveKernel(beams.DEUTERON, b1_func=toy_b1)
    sampler = tagged.TaggedSampler(central, kern, config, fom.Scenario())
    plan = bk.tensor_thirds_plan(PZ, PZZ)
    k_edges = np.linspace(0.0, 0.6, 13)
    kc = 0.5 * (k_edges[:-1] + k_edges[1:])
    rng = np.random.default_rng(args.seed)
    folded, n_gen = folded_asymmetry(sampler, plan, args.events, k_edges,
                                     menu, rng)
    sigma_pb = sum(sampler.sigma_tot_pb(c) * c.lumi_fraction
                   for c in plan.categories)
    summary = ["configuration %s (%s), %d events generated"
               % (args.config, config.label(), n_gen)]
    ax2.plot(central.k, azz_wf_curve(central, ic), "k-", lw=1, alpha=0.45,
             label=r"unfolded, $\theta_k=90^\circ$")
    for name, optics, colour, marker in menu:
        a, n, k_acc = folded[name]
        errs = np.where(n > 1, err_azz(np.maximum(n, 1), PZZ), np.nan)
        use = n > 50
        ax2.errorbar(kc[use], a[use], yerr=errs[use], fmt=marker, ms=4,
                     capsize=2, color=colour,
                     label="%s (%.2f x %.2f mrad)"
                           % (name, 1e3 * optics.envelope[0],
                              1e3 * optics.envelope[1]))
        eps = acceptance_table(central, config, optics)
        wf_acc = azz_wf_curve(central, weights=eps)
        ax2.plot(central.k, wf_acc, "-", lw=1.4, color=colour, alpha=0.85,
                 label="acceptance-weighted truth, " + name)
        # projected errors at the reference luminosity, at THIS optics'
        # luminosity: the tagging optics buys acceptance with luminosity
        scale = (args.lumi_ref * optics.lumi_fraction * 1e3 * sigma_pb) / n_gen
        n_ref = n * scale
        i300 = np.argmin(np.abs(kc - 0.3))
        # the tag is the UNBINNED accepted fraction: 9-11 % of the accepted
        # alpha spectators are above the 0.6 GeV/c right edge of the panel's
        # histogram, so n.sum() / n_gen would understate it by that much
        acc = k_acc.size / n_gen
        summary.append(
            "%-30s acc %.4f, L/L_HA %.3f, acc x L %.4f; median accepted k "
            "%.3f GeV/c, frac(k < 0.15) %.3f; at k = %.3f GeV/c A_zz = "
            "%+.3f (acceptance-weighted truth %+.3f, 90 deg curve %+.3f), "
            "delta(%g fb^-1/u) = %s"
            % (name, acc, optics.lumi_fraction, acc * optics.lumi_fraction,
               np.median(k_acc) if k_acc.size else np.nan,
               float(np.mean(k_acc < 0.15)) if k_acc.size else np.nan,
               kc[i300], a[i300] if n[i300] > 0 else np.nan,
               float(np.interp(kc[i300], central.k, wf_acc)),
               float(np.interp(kc[i300], central.k,
                               azz_wf_curve(central, ic))),
               args.lumi_ref,
               "%.4f" % err_azz(max(n_ref[i300], 1), PZZ)
               if n[i300] > 0 else "n/a"))
    ax2.set_xlabel(r"$k_{\rm rec}$ [GeV/$c$]")
    ax2.set_ylabel(r"$A_{zz}^{\rm tag}$ (RP-accepted)")
    ax2.set_title("acceptance-folded, tensor-thirds fills")
    ax2.set_xlim(0, 0.6)
    ax2.legend(fontsize=6.5)
    ax2.axhline(0, color="0.6", lw=0.5)

    fig.suptitle(r"$^6$Li($e,e^\prime\alpha$)X: tagged tensor asymmetry of the "
                 r"embedded deuteron, %s (TOY/scenario inputs)" % config.label()
                 + "\n" + r"coloured curves: the truth weighted by each "
                 r"optics' own $\theta_k$ acceptance -- what the markers "
                 r"measure", fontsize=9.5)
    fig.tight_layout()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = output_stem("money_tagged_azz_6Li", key, args.config, args.optics)
    out = outdir / (stem + ".png")
    fig.savefig(out, dpi=140)
    print("wrote", out)
    for line in summary:
        print(line)


if __name__ == "__main__":
    main()
