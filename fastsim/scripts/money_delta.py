#!/usr/bin/env python3
"""Money plot 3 (plan step 1.3.3): gluonometry discovery reach.

5-sigma luminosity for the cos(2phi) double-helicity-flip amplitude as a
function of the Delta/F1 scale (scenario range 1e-3 .. 1e-2 spans the
Sather-Schmidt bag estimate). Assumes transversely tensor-polarized 6Li
running with unpolarized electrons; significance combines all accepted
(x, Q2) bins:  sig^2 = sum_bins A_bin^2 * P_zz^2 * N_bin / 2.
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.asymmetries import a_cos2phi
from polli_fastsim.inputs import get_backends
from polli_fastsim.polarized import toy_delta_gluon
from polli_fastsim.structure import NuclearF2

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def bin_terms(cfg, scale, pzz, base=None):
    """Per-bin (sig^2 contribution, event count), both at 1 fb^-1/nucleon,
    over the accepted bins.

    Bins whose structure functions come back non-finite are dropped here,
    explicitly: a PDF grid outside its fit range returns NaN, and until
    2026-08-28 those bins were removed only as a side effect of the
    min-events comparison (NaN >= 10 is False), which stopped being true
    the moment the threshold was allowed to fall below one event."""
    sc = fom.Scenario(lumi_fb_per_nucleon=1.0, pol_ion_tensor=pzz)
    nf2_in = NuclearF2(cfg.ion, base=base) if base is not None else None
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_in)
    nf2 = proj.extras["nf2"]
    f1 = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2 = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y = proj.extras["y"]
    delta = toy_delta_gluon(proj.x, proj.q2, f1, scale=scale)
    amp = a_cos2phi(delta, f1, f2, proj.x, y)
    terms = amp**2 * pzz**2 * proj.n_events / 2.0
    ok = (proj.accepted & np.isfinite(terms) & np.isfinite(proj.n_events))
    return terms[ok].ravel(), proj.n_events[ok].ravel()


def reach_from_terms(terms, n_events, min_events=10, target=25.0):
    """L_5sigma [fb^-1/nucleon] from the per-bin sig^2 at 1 fb^-1/u, with
    the MIN-EVENTS floor applied at the luminosity the reach is quoted at.

    The floor is there to keep a Gaussian counting significance honest, so
    it belongs at the luminosity of the measurement, not at the 1 fb^-1/u
    the sig^2 is normalised to.  Until 2026-08-28 it was applied at
    1 fb^-1/u and the reach was then scaled up from the truncated sum: a
    bin holding two events per fb^-1 -- forty at a 20 fb^-1 reach, and a
    perfectly good bin there -- was discarded.

    L sig2(L) is non-decreasing in L (raising L both scales the
    significance and admits bins), so `L sig2(L) = target` has one
    solution and it is found exactly rather than iterated: bins enter in
    order of decreasing rate, bin k at L = min_events / n_k, so a walk
    over the cumulative sum lands on the interval containing its own
    solution.  Where the target is already met at the left edge of that
    interval the answer is the edge itself -- the first luminosity at
    which the bins carrying it have min_events each."""
    n_events = np.asarray(n_events, dtype=float)
    terms = np.asarray(terms, dtype=float)
    if n_events.size == 0:
        return np.inf
    order = np.argsort(-n_events)
    n_sorted, cum = n_events[order], np.cumsum(terms[order])
    if not cum[-1] > 0:            # no bin carries any significance
        return np.inf
    enter = min_events / np.maximum(n_sorted, 1e-300)
    for k in range(n_sorted.size):
        lumi = target / cum[k] if cum[k] > 0 else np.inf
        hi = enter[k + 1] if k + 1 < n_sorted.size else np.inf
        if lumi < hi:
            return float(max(lumi, enter[k]))
    return np.inf


def sig2_per_fb_at(cfg, scale, pzz, base=None, min_events=10, lumi_fb=None):
    """Significance^2 per fb^-1/nucleon at a reference Delta/F1 scale.
    The toy Delta shape is linear in `scale`, so reach curves follow
    L_5sig(s) = 25 / (sig2 * (s/scale)^2) analytically -- but the
    min-events floor does not scale with it, so it is applied at
    `lumi_fb`, the luminosity the reach is evaluated at (None = the
    self-consistent one, `reach_fb`; see `reach_from_terms`)."""
    terms, n_events = bin_terms(cfg, scale, pzz, base=base)
    if lumi_fb is None:
        lumi_fb = reach_from_terms(terms, n_events, min_events=min_events)
    return float(terms[n_events * lumi_fb >= min_events].sum())


def reach_fb(cfg, scale, pzz, base=None, min_events=10, target=25.0):
    """L_5sigma [fb^-1/nucleon] at this Delta/F1 scale (target = 25 is
    5 sigma), with the min-events floor at the luminosity it is quoted
    at."""
    terms, n_events = bin_terms(cfg, scale, pzz, base=base)
    return reach_from_terms(terms, n_events, min_events=min_events,
                            target=target)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="6Li", choices=["6Li", "7Li"])
    ap.add_argument("--pdf", default="toy", choices=["toy", "grid"])
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    backends = get_backends(args.pdf)
    scales = np.logspace(-3.3, -1.7, 15)
    s0 = 1e-3
    fig, ax = plt.subplots(figsize=(7, 5))
    reach_ref = {}
    for cfg, color in zip(beams.default_configs(args.ion),
                          ("crimson", "seagreen", "navy")):
        for pzz, ls in ((0.60, "-"), (0.80, "--")):
            # the amplitude is linear in the scale, so the per-bin sig^2
            # terms scale as (s/s0)^2 and only the min-events floor has to
            # be re-solved per point
            terms, n_events = bin_terms(cfg, s0, pzz, base=backends["base"])
            reach = np.array([reach_from_terms(terms * (s / s0) ** 2,
                                               n_events) for s in scales])
            reach_ref[(cfg.label(), pzz)] = reach_from_terms(terms, n_events)
            ax.plot(scales, reach, ls, color=color, lw=1.5,
                    label=f"{cfg.label()}, $P_{{zz}}$={pzz:g}")
    ax.axhspan(1, 100, color="gold", alpha=0.12,
               label="1-100 fb$^{-1}$/u (plausible program)")
    ax.axvline(1e-3, color="gray", ls=":", lw=1)
    ax.text(1.05e-3, 0.93, "Sather-Schmidt\n$O(10^{-3})$", fontsize=7,
            transform=ax.get_xaxis_transform(), va="top")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta/F_1$ scale (peak of scenario shape)")
    ax.set_ylabel(r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]")
    ax.set_title(f"Nuclear gluonometry reach, transversely polarized "
                 f"{args.ion}\n(cos 2$\\phi$ amplitude, all bins combined; "
                 f"{backends['tag'].upper()} inputs)", fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    path = outdir / f"money_delta_{args.ion}_{backends['tag']}.png"
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")
    for cfg in beams.default_configs(args.ion):
        l5 = reach_ref[(cfg.label(), 0.8)]
        l5_hi = reach_fb(cfg, 1e-2, 0.8, base=backends["base"])
        print(f"  {cfg.label():26s} L_5sig(Delta/F1=1e-3, Pzz=0.8) = "
              f"{l5:9.1f} fb^-1/u ; (1e-2) = {l5_hi:7.3f}")


if __name__ == "__main__":
    main()
