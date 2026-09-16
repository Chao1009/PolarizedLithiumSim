#!/usr/bin/env python3
"""Money plot 5: the cos 2phi gluonometry measurement as projected data.

Left 2x2: phi'-modulation pseudo-data (points with statistical error bars
at the 1-YEAR EIC luminosity) in the four sweet-spot (x, Q2) super-bins --
the bins with the highest projected cos 2phi significance, spread in
(x, Q2).  Right: the extracted amplitude vs x along the sweet-spot Q2
slice with 1-year AND 10-year error bars and the Delta-model curves.

The Delta(x, Q2) input comes from the unified model registry
`polli_fastsim.delta_models` (single home for all Delta models and
constraints): default `moment_A` -- the sum-rule-constrained ansatz of
the money_delta suite, Delta = A alpha_s F1 x^a(1-x)^b with A solved
from int x Delta dx = -0.012 alpha_s (Sather-Schmidt bag moment) --
with the 6Li per-nucleon dilution 1/3 (two polarized nucleons of six;
the whole-nucleus P_zz stays a separate factor).  It is a counting
fraction, not the vector polarization of `beams.LI6` -- see the
convention note in `delta_models`.  `moment_B` (no F1 factor,
the conservative reading) and the legacy `toy` shape are one flag away.

Detection assumption: scattered electron only (Scenario e' cuts).
Statistics via per-phi-bin Poisson draws from exact expected yields
(sample.phi_histogram_pseudo) -- identical to event-level sampling for
the binned estimator at any luminosity.  The whole luminosity is
assigned to the transverse-tensor fill (FOM-map convention).

Estimator: ONE fill at --pzz and the binned cos 2phi' fit, which is what
every published number here is.  `--pzz-plus/--pzz-zero` instead splits
the same luminosity between an m = +-1-rich and an m = 0-rich fill and
reads them with the spin-state-sorted ratio reco.harmonic_ratio_fit --
the estimator of the measurement (plans/07 WP3), worth 0.67x on dA at
(+0.6, -1.2) -- on its own '_twofill' stem.

Usage:  python3 scripts/money_cos2phi.py --delta-model moment_A
        python3 scripts/money_cos2phi.py --pzz-plus 0.6 --pzz-zero -1.2
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

from polligen import bookkeeping as bk  # noqa: E402
from polligen import reco  # noqa: E402
from polligen import spin as spinmod  # noqa: E402
from polligen.estimators import (cos2phi_fit_binned,  # noqa: E402
                                 cos2phi_fit_err)
from polligen.sample import InclusiveSampler, phi_histogram_pseudo  # noqa: E402
from polligen.xsec import (InclusiveKernel,  # noqa: E402
                           tensor_leakage_amplitude)

from polli_fastsim import beams, delta_models as dm, fom  # noqa: E402
PRETTY = {"moment_A": "moment-constrained $\\Delta$", "moment_B": "no-$F_1$ variant", "toy": "flat toy"}  # figure labels
from polli_fastsim import structure  # noqa: E402
from polli_fastsim.asymmetries import a_cos2phi  # noqa: E402
from polli_fastsim.inputs import get_backends  # noqa: E402
from polli_fastsim.kinematics import kinematic_mask, y_from_xq2  # noqa: E402
from polli_fastsim.polarized import toy_b1  # noqa: E402

# Okabe-Ito (colorblind-safe): blue = injected, vermillion = fit/1-yr,
# green = alternative interpretation
C_TRUTH, C_FIT, C_ALT = "#0072B2", "#D55E00", "#009E73"

# --- the O(gamma^2) tensor-leakage switches, shared by the three money
# scripts (plans/08 D2).  EVERY DEFAULT IS OFF: with none of these flags
# given, the kernel is the massless Hoodbhoy-Jaffe-Manohar one and every
# published number is bit for bit what it has always been.

# What the subtraction cannot reach.  The correction removes the
# b3 = b4 = 0 leakage; at b3 = b4 = 0.1 b2 the true leakage is 1.60-1.65x
# that (xsec.tensor_leakage_ratio, measured at all twelve sweet spots), so
# ~0.6 of the subtracted value survives as an irreducible model band.  It
# is the reason the subtraction shrinks the systematic by 1.6 and not by
# an order of magnitude.
#
# It is a PRIOR width on the unmeasured b3, b4 and not the residual
# against the b3, b4 a given run assumes, so it is a fixed constant: with
# --b3-frac 0.1 --b4-frac 0.1 the generator and the subtraction share the
# same higher twist and `tensor_gamma_leakage.py` prints a zero residual
# against the reference, while this band stays open, because an
# experiment does not know that it guessed right.  The in-situ route
# carries a second, MEASURED width on top of it (the residual mis-scaling
# of its rate sector); money_cos2phi_reco.py adds the two in quadrature.
LEAKAGE_B34_BAND = 0.6


def add_tensor_leakage_args(ap, with_scan=False):
    """The four (or five) tensor-leakage flags, identical in all three
    money scripts."""
    ap.add_argument("--tensor-gamma", action="store_true",
                    dest="tensor_gamma",
                    help="generate with the EXACT finite-gamma tensor "
                         "b-sector (Cosyn Eqs. 9/10/14/16/17/24, "
                         "xsec.InclusiveKernel(tensor_gamma=True)) instead "
                         "of the massless Hoodbhoy-Jaffe-Manohar one.  OFF "
                         "by default: every published number is on the "
                         "massless path, and switching this on re-draws "
                         "every pseudo-measurement (the leakage moves the "
                         "phi-independent harmonic too, hence the expected "
                         "counts).  A run with it on writes its own PNG")
    ap.add_argument("--subtract-tensor-leakage", default="none",
                    choices=("none", "kappa", "model"),
                    dest="subtract_tensor_leakage",
                    help="subtract the O(gamma^2) b1-b4 leakage from the "
                         "extracted cos 2phi amplitude (and so from "
                         "Delta-hat).  'none' (default) carries it as a "
                         "systematic, which is what every published number "
                         "does.  'kappa' takes the rate sector from the "
                         "data -- the constant term of the same ratio fit, "
                         "with the bin-independent pedestal a "
                         "relative-luminosity error leaves on it measured "
                         "across bins and removed -- and only the kinematic "
                         "ratio L from the model; 'model' takes the whole "
                         "leakage from the "
                         "b1 model, which is a closure test and not a "
                         "publishable route.  Requires --tensor-gamma: "
                         "subtracting a leakage the pseudo-data do not "
                         "carry would bias them")
    ap.add_argument("--b3-frac", type=float, default=0.0, dest="b3_frac",
                    help="b3 as a fraction of b2 (default 0).  The "
                         "unmeasured higher-twist slot that breaks the "
                         "3 : -3 : 1 cancellation the small leakage rests "
                         "on; --tensor-gamma only")
    ap.add_argument("--b4-frac", type=float, default=0.0, dest="b4_frac",
                    help="b4 as a fraction of b2 (default 0); "
                         "--tensor-gamma only")
    if with_scan:
        ap.add_argument("--leakage-scan", action="store_true",
                        dest="leakage_scan",
                        help="print the leakage budget per sweet spot -- "
                             "A_leak, L = A_leak/kappa, the fitted kappa, "
                             "the shift of Delta-hat the subtraction makes "
                             "and the b3/b4 band it leaves.  A pure "
                             "addition to the output: no figure changes")


def check_tensor_leakage_args(ap, args):
    """Refuse the combinations that would silently bias a number."""
    if args.subtract_tensor_leakage != "none" and not args.tensor_gamma:
        ap.error("--subtract-tensor-leakage needs --tensor-gamma: on the "
                 "massless path the pseudo-data carry no leakage, so "
                 "subtracting one would bias the extraction by its full "
                 "size")
    if (args.b3_frac or args.b4_frac) and not args.tensor_gamma:
        ap.error("--b3-frac/--b4-frac need --tensor-gamma: only the exact "
                 "finite-gamma b-sector reads b3 and b4")


def b34_funcs(args):
    """(b3_func, b4_func) as fractions of b2 = 2 x b1 on the toy b1 --
    the convention of `scripts/tensor_gamma_leakage.py`."""
    def frac(f):
        return ((lambda xx, qq, f1: f * 2.0 * xx * toy_b1(xx, qq, f1))
                if f else None)
    return frac(args.b3_frac), frac(args.b4_frac)


def tensor_leakage_tag(args):
    """Filename key for a non-default tensor-leakage setting ('' at the
    published defaults).

    The published PNGs belong to the massless path with no subtraction; a
    run with any of these switches appends its key so that it cannot
    overwrite them -- the same guard as `fom.run_share_tag` and
    `money_tagged_azz.output_stem`."""
    keys = []
    if getattr(args, "tensor_gamma", False):
        keys.append("tgamma")
    sub = getattr(args, "subtract_tensor_leakage", "none")
    if sub != "none":
        keys.append("sub" + sub)
    for name, val in (("b3f", getattr(args, "b3_frac", 0.0)),
                      ("b4f", getattr(args, "b4_frac", 0.0))):
        if val:
            keys.append(("%s%g" % (name, val)).replace(".", "p"))
    return "_".join(keys)


def add_pdf_arg(ap):
    """`--pdf {toy,grid}` -- the structure-function backend.

    The spelling is the one the six scripts that already take a
    backend use -- five in `fastsim/scripts` (`money_delta_20260715.py`,
    `money_b1.py`, `money_delta.py`, `coverage_and_stat_maps.py`,
    `money_polemc.py`) and `evgen/scripts/target_mass_bound.py`;
    plans/07 WP1 wrote `--backend grid`, a flag that has never existed
    anywhere in this repository.
    """
    ap.add_argument("--pdf", default="toy", choices=("toy", "grid"),
                    help="structure-function backend.  'toy' (the "
                         "default, and every published PNG) is the toy F2 "
                         "and the toy R; 'grid' is the published-input "
                         "path -- F2A from the EPPS21nlo_CT18Anlo_Li6 "
                         "NUCLEAR set, g1 from NNPDFpol11_100, and the "
                         "SLAC/E143 R1998 fit in place of the toy R.  "
                         "Needs `parton` with those grids installed, and "
                         "writes its own '_grid' PNG stem")


def pdf_backends(args, ion):
    """The backend set `--pdf` names, for `ion`.

    One call gives the three objects that have to move together -- the
    free-nucleon F2 the g1 denominator is built on, the g1 model, and the
    whole-nucleus F2A -- plus, on 'grid', the one R they all use.  On
    'toy' with r_func None it is bit-for-bit the bare `NuclearF2(ion)` /
    `ToyG1()` these scripts built inline before run 19.
    """
    pdf = getattr(args, "pdf", "toy")
    return get_backends(pdf, nuclear=ion,
                        r_func=structure.r1998 if pdf == "grid" else None)


def pdf_tag(args):
    """Filename key for a non-default `--pdf` ('' at the published toy)."""
    return "grid" if getattr(args, "pdf", "toy") == "grid" else ""


def add_binning_arg(ap):
    """`--binning {log40x30,yr}` -- the (x, Q2) ANALYSIS grid.

    The published maps run this project's own 40 x 30 log grid, which is
    10.0 x 9.1 bins per decade and matches no published projection bin for
    bin.  `yr` is the Yellow-Report inclusive convention, five logarithmic
    bins per decade in both variables on the decade-anchored lattice
    (`beams.YR_GRID` / `beams.YR_X_EDGES` / `beams.YR_Q2_EDGES`), which is
    what makes a comparison with published YR projections one-to-one --
    plans/02 Step 1.1 item 3.  `generic` is an accepted alias of the
    default spelling.
    """
    ap.add_argument("--binning", default="log40x30",
                    choices=("log40x30", "generic", "yr"),
                    help="(x, Q2) analysis binning.  'log40x30' (the "
                         "default, and every published PNG; 'generic' is "
                         "the same grid) is this project's 40 x 30 log "
                         "grid over x in [1e-4, 1], Q2 in [1, 2e3]; 'yr' "
                         "is the Yellow-Report five-bins-per-decade "
                         "lattice (20 x 17 over x in [1e-4, 1], Q2 in "
                         "[1, 2512]) and writes its own '_yr' PNG stem")


def binning_tag(args):
    """Filename key for a non-default `--binning` ('' on the published
    40 x 30 grid)."""
    return beams.binning_tag(getattr(args, "binning", "log40x30"))


# --- the two-fill (spin-state-sorted) switches ----------------------------
#
# The published truth-level numbers come from a SINGLE tensor-polarized
# fill at P_zz = +0.6 and the binned cos 2phi' fit
# (`estimators.cos2phi_fit_binned`), whose statistical error is
# sqrt(2/N)/P_zz / dilution.  The reconstructed-level chain instead runs
# the acceptance-cancelling pattern of plans/07 WP3 -- an m = +-1-rich fill
# at P_zz = +0.6 alternating bunch by bunch with an m = 0-rich fill at
# -1.2 -- and the spin-state-sorted ratio estimator
# `reco.harmonic_ratio_fit`, whose error is 2 sqrt(2/N)/(P+ - P0) /
# dilution at equal luminosity.  At (+0.6, -1.2) that is 0.6/0.9 = 2/3 of
# the single-fill error on the SAME total luminosity: the m = 0-rich fill
# carries twice the tensor lever arm of the fill it replaces.
#
# BOTH FLAGS DEFAULT TO None: with neither given the drivers take exactly
# today's single-fill path and every published line and PNG is bit for
# bit what it has always been.  Given together they switch the estimator
# and append '_twofill' to the stem.  The luminosity split is fixed at
# 50/50, because the 0.67 rests on it and a free share under the same
# stem key would put two different numbers behind one file name.


def add_two_fill_args(ap):
    """`--pzz-plus/--pzz-zero` -- the two-fill spin-state-sorted path."""
    ap.add_argument("--pzz-plus", type=float, default=None, dest="pzz_plus",
                    help="tensor polarization of the m = +-1-rich fill.  "
                         "Given together with --pzz-zero it switches the "
                         "driver from the single-fill binned cos 2phi' fit "
                         "to the spin-state-sorted ratio estimator "
                         "reco.harmonic_ratio_fit at equal luminosity in "
                         "the two fills, and writes its own '_twofill' "
                         "PNG stem.  The published run is (0.6, -1.2)")
    ap.add_argument("--pzz-zero", type=float, default=None, dest="pzz_zero",
                    help="tensor polarization of the m = 0-rich fill "
                         "(-2 x --pzz-plus at the same source purity); "
                         "--pzz-plus only")


def check_two_fill_args(ap, args):
    """Refuse a half-specified or degenerate two-fill run."""
    if (args.pzz_plus is None) != (args.pzz_zero is None):
        ap.error("--pzz-plus and --pzz-zero must be given together: the "
                 "ratio estimator needs both fills, and one alone would "
                 "silently fall back to the single-fill path under a "
                 "stem that claims otherwise")
    if args.pzz_plus is None:
        return
    if args.pzz_plus == args.pzz_zero:
        ap.error("--pzz-plus and --pzz-zero must differ: the estimator "
                 "divides by (P+ - P0)")


def two_fill_pzz(args):
    """(P+, P0) when the two-fill path is on, else None."""
    if getattr(args, "pzz_plus", None) is None:
        return None
    return (float(args.pzz_plus), float(args.pzz_zero))


def two_fill_tag(args):
    """Filename key for the two-fill path ('' on the published
    single-fill one)."""
    return "twofill" if two_fill_pzz(args) is not None else ""


def two_fill_plan(pzz_pair, phi_s=0.0, name="cos2phi"):
    """Equal-luminosity two-fill run plan at the two tensor polarizations.

    `bookkeeping.tensor_flip_plan` is the same object at the fixed
    purity-matched pair (P, -2P); this builder takes the two values
    independently so that a run can price a source that does not deliver
    the m = 0-rich bunches at the same purity."""
    cats = [bk.SpinCategory("%s%+g" % (name, p), 1.0,
                            spinmod.spin1_populations(0.0, p),
                            theta_s=np.pi / 2.0, phi_s=phi_s,
                            lumi_fraction=0.5)
            for p in pzz_pair]
    return bk.RunPlan(cats, pzz_true=pzz_pair[0])


def describe_backends(args, backends):
    """One line naming what `--pdf` actually loaded, for the run record."""
    nf2 = backends["nuclear"]
    r = "R1998" if getattr(nf2, "r_func", None) is not None else "toy R"
    if getattr(args, "pdf", "toy") != "grid":
        return "toy F2/g1, %s (published path)" % r
    return ("F2A = %s, g1 = NNPDFpol11_100, F2p = CT18NLO, %s; F2A FROZEN "
            "below Q2 = %.4g GeV^2 (the set's Q0^2 -- parton returns NaN "
            "there, and the money maps start at Q2 = 1)"
            % (nf2.setname, r, nf2.q2_min))


def output_stem_tag(args):
    """The whole non-default-settings key of a truth-level money PNG.

    `--pdf grid` first, then the binning, then the two-fill estimator,
    then the tensor-leakage keys, so that the published toy massless-path
    single-fill stem stays the bare one and any other run writes beside
    it rather than over it -- the guard convention of
    `tensor_leakage_tag` below, extended to the backend, the analysis
    grid and the estimator.
    """
    return "_".join(k for k in (pdf_tag(args), binning_tag(args),
                                two_fill_tag(args),
                                tensor_leakage_tag(args)) if k)


def truth_leakage_route(args):
    """How the TRUTH-LEVEL scripts should name the route they ran.

    There is no fitted constant at truth level, so 'kappa' and 'model'
    perform the same subtraction here -- the model's own folded rate
    sector.  The distinction is a reconstructed-level one
    (money_cos2phi_reco.py, where kappa_hat is a measured quantity), and
    a run record that said 'kappa' unqualified would claim an in-situ
    route this script cannot take."""
    route = args.subtract_tensor_leakage
    if route == "kappa":
        return ("kappa (identical to model at truth level: there is no "
                "fitted constant here)")
    return route


def leakage_amplitude(sampler, cat, mask, sigma_pb):
    """The cos 2phi' tensor leakage of a truth-level super-bin, per unit
    P_zz, on the SAME rate weights `measure` averages the amplitude over.

    `sampler.effective_modulation` returns a2_eff = sum_m p_m sum_cells
    xsec a_2 divided by the total sigma; the b-sector part of that
    numerator is P_zz sum_cells xsec h2 because every harmonic is linear
    in the rank-2 moment, so this is a2_eff's leakage half divided by
    P_zz.  It is what `--subtract-tensor-leakage` removes from both the
    fitted amplitude and the truth reference."""
    sel = np.asarray(mask, dtype=bool)
    y = sampler.q2_cells / (sampler.s * sampler.x_cells)
    _h0, h2 = tensor_leakage_amplitude(sampler.tables, sampler.x_cells,
                                       sampler.q2_cells, y,
                                       theta_s=cat.theta_s, j=cat.j)
    if sigma_pb <= 0.0:
        return 0.0
    return float((sampler.xsec_flat[sel] * np.asarray(h2)[sel]).sum()
                 / sigma_pb)


def pick_sweet_spots(proj, sig, n=4, min_dx=0.35):
    """Greedy top-significance bins, separated by min_dx in log10(x)."""
    order = np.argsort(sig, axis=None)[::-1]
    picked = []
    for flat in order:
        i, j = np.unravel_index(flat, sig.shape)
        if sig[i, j] <= 0:
            break
        x, q2 = proj.x[i, j], proj.q2[i, j]
        if not any(abs(np.log10(x / xs)) < min_dx
                   for xs, *_ in picked):
            picked.append((x, q2, i, j))
        if len(picked) == n:
            break
    return picked


def pick_sweet_spots_banded(proj, sig,
                            bands=((1.0, 3.0, 2), (3.0, 12.0, 1),
                                   (12.0, 200.0, 1))):
    """Best bins per Q2 band: the raw significance map clusters at the
    lowest Q2 (rate-dominated), but the money plot must also show the
    Q2 lever arm, so one spot is forced into each higher band.

    Audit note (2026-08-10): the ranking runs on the 40x30 analysis map
    with cell-center acceptance; counts of bins touching the y_min edge
    carry ~20% grid bias there.  All PLOTTED numbers instead come from
    the finer 60x45 sampler grid (within ~3% of a 240x180 reference),
    so the bias can only influence which near-edge bin gets picked."""
    picked = []
    for qlo, qhi, npick in bands:
        in_band = np.where((proj.q2 >= qlo) & (proj.q2 < qhi), sig, 0.0)
        picked += pick_sweet_spots(proj, in_band, n=npick)
    return picked


def superbin_edges(proj, i, j, pad=1):
    """(x_lo, x_hi, q2_lo, q2_hi) of the FOM bin padded by `pad` bins."""
    xe, qe = proj.x_edges, proj.q2_edges
    return (xe[max(i - pad, 0)], xe[min(i + pad + 1, xe.size - 1)],
            qe[max(j - pad, 0)], qe[min(j + pad + 1, qe.size - 1)])


def superbin_mask(sampler, xlo, xhi, q2lo, q2hi):
    return ((sampler.x_cells >= xlo) & (sampler.x_cells < xhi)
            & (sampler.q2_cells >= q2lo) & (sampler.q2_cells < q2hi))


def measure(sampler, cat, mask, lumi_pb, pzz, rng, nbins=24):
    """One full-luminosity binned pseudo-measurement on a super-bin."""
    sigma_pb, _a1, a2_eff = sampler.effective_modulation(cat, mask=mask)
    n_exp = lumi_pb * sigma_pb
    counts, edges = phi_histogram_pseudo(n_exp, a2_eff, nbins=nbins, rng=rng)
    amp_hat = cos2phi_fit_binned(counts, edges, pzz)
    return {"n": n_exp, "truth": a2_eff / pzz, "amp": amp_hat,
            "err": cos2phi_fit_err(n_exp, pzz, nbins),
            "counts": counts, "edges": edges, "a2_eff": a2_eff,
            "sigma_pb": sigma_pb}


def measure_two_fill(sampler, cats, mask, lumi_pb, pzz_pair, rng, nbins=24):
    """One two-fill spin-state-sorted pseudo-measurement on a super-bin.

    The same TOTAL luminosity as `measure`, split 50/50 between the two
    fills.  Each fill is drawn exactly as the single-fill path draws its
    one sample -- the category's own accepted cross section sigma_f and
    its own effective modulation a2_f from
    `sampler.effective_modulation`, Poisson per phi' bin -- and the two
    count rows go into `reco.harmonic_ratio_fit`, which forms the
    acceptance-free spin-state ratio, inverts it for the modulation
    T_i = kappa + A cos 2phi' and fits A by weighted LSQ.

    The estimator's target is EXACT here, not linearized.  Both sigma_f
    and the numerator of a2_f are linear in P_zz (the populations are),
    so sigma_f = sigma_0 (1 + P_f kappa) and a2_f = P_f A/(1 + P_f kappa),
    and the drawn yields are L_f sigma_0 (1 + P_f kappa +
    P_f A cos 2phi') bin by bin -- the estimator's model with no
    remainder.  Two fills therefore determine (sigma_0, kappa) and hence
    the amplitude A the fit should return; `truth` below is that A, and
    `truth_single` is the single-fill convention a2(P+)/P+ that the
    published numbers quote, so the two can be compared without confusing
    a change of estimator with a change of definition.
    """
    lumis = [0.5 * lumi_pb, 0.5 * lumi_pb]
    p_plus, p_zero = float(pzz_pair[0]), float(pzz_pair[1])
    sig, a2 = [], []
    counts, edges = [], None
    for cat, lf in zip(cats, lumis):
        sigma_pb, _a1, a2_eff = sampler.effective_modulation(cat, mask=mask)
        sig.append(sigma_pb)
        a2.append(a2_eff)
        c, edges = phi_histogram_pseudo(lf * sigma_pb, a2_eff,
                                        nbins=nbins, rng=rng)
        counts.append(c)
    counts = np.asarray(counts, dtype=float)
    n_exp = float(sum(l * s for l, s in zip(lumis, sig)))
    # (sigma_0, kappa) from the two fills, then the exact amplitude
    denom = p_plus - p_zero
    sigma_0 = (sig[1] * p_plus - sig[0] * p_zero) / denom
    if sigma_0 > 0.0:
        amp_true = a2[0] * sig[0] / (sigma_0 * p_plus) if p_plus else 0.0
        kappa_true = (sig[0] - sig[1]) / (sigma_0 * denom)
    else:
        amp_true, kappa_true = 0.0, 0.0
    fit = reco.harmonic_ratio_fit(counts, lumis, [p_plus, p_zero], edges)
    # the per-bin modulation the fit runs on, for the figure: the exact
    # inversion T = R/(sigma_P^2 - Pbar R) of the same spin-state ratio
    r, var_r, sig2, pbar = reco.spin_state_ratio(counts, lumis,
                                                 [p_plus, p_zero])
    t_bin = r / (sig2 - pbar * r)
    jac = (1.0 + pbar * t_bin) / (sig2 - pbar * r)
    return {"n": n_exp, "truth": amp_true,
            "truth_single": (a2[0] / p_plus) if p_plus else 0.0,
            "const_truth": kappa_true,
            "amp": fit["amp"], "err": fit["err"], "const": fit["const"],
            "t": t_bin, "t_err": np.sqrt(np.maximum(var_r, 0.0)) * np.abs(jac),
            "counts": counts, "edges": edges,
            "a2_eff": a2[0], "sigma_pb": sigma_0,
            "sigma_fills": tuple(sig),
            "err_analytic": float(reco.err_harmonic_ratio(
                max(n_exp, 1e-12), [p_plus, p_zero], nbins=nbins))}


def two_fill_err_ratio(pzz_single, pzz_pair):
    """dA(two-fill) / dA(single-fill) at the same total N.

    sqrt(2/N)/sigma_P vs sqrt(2/N)/P_single, i.e. P_single/sigma_P; the
    finite-bin dilution is the same factor in both and cancels.  At
    (0.6; +0.6, -1.2) it is 0.6/0.9 = 2/3."""
    p = np.asarray(pzz_pair, dtype=float)
    pbar = p.mean()
    sigma_p = float(np.sqrt(((p - pbar) ** 2).mean()))
    return abs(pzz_single) / sigma_p


def build_delta_model(args, config, scenario, backends=None):
    """Delta model from the unified registry; moment_A needs the
    rate-weighted <Q2> of the accepted phase space and a per-nucleon
    F1 handle.

    Both come off the `--pdf` backend since run 19: on `--pdf grid` the
    sum rule is re-solved against the grid F1 at the GRID <Q2>, which is
    the whole of what "re-solve moment_A at the grid <Q2>" means
    (plans/07 WP1).  `backends` is the already-built set when the caller
    has one, so the grids are opened once per run.
    """
    if args.delta_model == "toy":
        return dm.make("toy", scale=args.scale), None
    if args.delta_model == "moment_B":
        return dm.make("moment_B", variant=args.variant,
                       dilution=args.dilution), None
    nf2 = (backends or pdf_backends(args, config.ion))["nuclear"]
    proj0 = fom.project_rates(config, scenario, nuclear_f2=nf2)
    acc = proj0.accepted
    q2_ref = float((proj0.n_events[acc] * proj0.q2[acc]).sum()
                   / proj0.n_events[acc].sum())
    model = dm.make(
        "moment_A", f1_func=lambda x, q2: nf2.f1a(x, q2) / config.ion.A,
        q2_ref=q2_ref, variant=args.variant, dilution=args.dilution)
    return model, q2_ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=int, default=1, choices=(0, 1, 2),
                    help="beam-energy point (default mid, e10 x 6Li 99.5 GeV/u)")
    ap.add_argument("--delta-model", default="moment_A",
                    choices=dm.available())
    ap.add_argument("--variant", default="mid_x",
                    choices=sorted(dm.VARIANTS))
    ap.add_argument("--dilution", type=float, default=1.0 / 3.0,
                    help="6Li per-nucleon dilution folded into Delta "
                         "(2 of 6 nucleons); P_zz stays whole-nucleus")
    ap.add_argument("--scale", type=float, default=1e-2,
                    help="peak Delta/F1 (toy model only)")
    ap.add_argument("--lumi-1yr", type=float, default=10.0,
                    help="1-year EIC program [fb^-1/nucleon]")
    ap.add_argument("--lumi-10yr", type=float, default=100.0,
                    help="10-year EIC program [fb^-1/nucleon]")
    ap.add_argument("--pzz", type=float, default=0.60)
    ap.add_argument("--nspots", type=int, default=4)
    add_pdf_arg(ap)
    add_two_fill_args(ap)
    add_tensor_leakage_args(ap)
    ap.add_argument("--seed", type=int, default=20260810)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    check_tensor_leakage_args(ap, args)
    check_two_fill_args(ap, args)

    config = beams.default_configs("6Li")[args.config]
    lumi1_pb = args.lumi_1yr * 1e3
    lumi10_pb = args.lumi_10yr * 1e3
    rng = np.random.default_rng(args.seed)

    scenario = fom.Scenario(lumi_fb_per_nucleon=args.lumi_1yr,
                            pol_ion_tensor=args.pzz)
    backends = pdf_backends(args, config.ion)
    model, q2_ref = build_delta_model(args, config, scenario,
                                      backends=backends)
    delta_func = model  # DeltaModel is (x, q2, f1)-callable

    # --- sweet spots from the analytic significance map -------------------
    proj = fom.project_rates(config, scenario,
                             nuclear_f2=backends["nuclear"])
    b3_func, b4_func = b34_funcs(args)
    kern = InclusiveKernel(config.ion, b1_func=toy_b1, delta_func=delta_func,
                           b3_func=b3_func, b4_func=b4_func,
                           tensor_gamma=args.tensor_gamma,
                           f2_source=backends["base"],
                           g1_model=backends["g1"],
                           nuclear_f2=backends["nuclear"],
                           r_func=backends["nuclear"].r_func)
    obs = fom.project_observables(config, scenario, proj,
                                  kern.g1_model, toy_b1, delta_func)
    spots = pick_sweet_spots_banded(proj, obs["sig_a_cos2phi"])[:args.nspots]

    # --- one sampler over the full range; super-bins are cell masks ------
    sampler = InclusiveSampler(kern, config, scenario, nx=60, nq2=45)
    plan = bk.transverse_tensor_plan(args.pzz)
    cat = plan.categories[0]
    subtract = args.subtract_tensor_leakage != "none"

    # the two-fill path: the same super-bins and the same total
    # luminosity, drawn in two fills and read with the spin-state-sorted
    # ratio estimator instead of the single-fill binned fit
    pzz_pair = two_fill_pzz(args)
    two_fill = pzz_pair is not None
    plan2 = two_fill_plan(pzz_pair) if two_fill else None

    def do_measure(mask, lumi_pb):
        if two_fill:
            return measure_two_fill(sampler, plan2.categories, mask,
                                    lumi_pb, pzz_pair, rng)
        return measure(sampler, cat, mask, lumi_pb, plan.pzz_true, rng)

    def single_fill_err(mask, lumi_pb):
        """dA the published single-fill estimator would give on the same
        super-bin and the same total luminosity -- the denominator of the
        0.67 the two-fill path is priced against."""
        sigma_pb, _a1, _a2 = sampler.effective_modulation(cat, mask=mask)
        return cos2phi_fit_err(lumi_pb * sigma_pb, plan.pzz_true)

    def leak_of(mask, m):
        """The leakage to remove from one super-bin's amplitude and truth,
        per unit P_zz.  At the truth level the 'kappa' and 'model' routes
        coincide -- there is no fitted kappa here, only the model's own
        rate sector -- so both take the folded leakage itself; the
        distinction is a reconstructed-level one (money_cos2phi_reco.py).
        """
        if not subtract:
            return 0.0
        return leakage_amplitude(sampler, cat, mask, m["sigma_pb"])

    fig = plt.figure(figsize=(12.5, 6.8))
    gs = GridSpec(2, 3, figure=fig, width_ratios=(1, 1, 1.35),
                  hspace=0.42, wspace=0.30)
    summary = []

    for k, (xs, qs, i, j) in enumerate(spots):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        xlo, xhi, q2lo, q2hi = superbin_edges(proj, i, j)
        mask = superbin_mask(sampler, xlo, xhi, q2lo, q2hi)
        m = do_measure(mask, lumi1_pb)
        leak = leak_of(mask, m)
        if leak:
            # additive on the amplitude; the truth reference and the drawn
            # curve lose the same leakage, so what is plotted and what is
            # fitted stay the same quantity (plans/08 D2 risk R3)
            m["amp"] -= leak
            m["truth"] -= leak
            m["a2_eff"] = m["a2_eff"] - leak * plan.pzz_true
        err10 = (float(reco.err_harmonic_ratio(lumi10_pb * m["n"]
                                               / max(lumi1_pb, 1e-30),
                                               pzz_pair))
                 if two_fill else
                 cos2phi_fit_err(lumi10_pb * m["sigma_pb"], plan.pzz_true))
        centers = 0.5 * (m["edges"][:-1] + m["edges"][1:])
        phi = np.linspace(0, 2 * np.pi, 200)
        if two_fill:
            # the observable of the two-fill analysis is the spin-state
            # ratio inverted for the modulation, T = kappa + A cos 2phi';
            # the fitted pedestal is removed so that the panel shows the
            # same thing the single-fill one does
            ax.errorbar(centers, 1e3 * (m["t"] - m["const"]),
                        yerr=1e3 * m["t_err"], fmt="o", color="black",
                        ms=3.5, capsize=2, lw=1, zorder=3)
            ax.plot(phi, 1e3 * m["truth"] * np.cos(2 * phi), "-",
                    color=C_TRUTH, lw=1.6)
            ax.plot(phi, 1e3 * m["amp"] * np.cos(2 * phi), "--",
                    color=C_FIT, lw=1.4)
        else:
            nbar = m["counts"].mean()
            mod = 1e3 * (m["counts"] / nbar - 1.0)
            mod_err = 1e3 * np.sqrt(np.maximum(m["counts"], 1.0)) / nbar
            ax.errorbar(centers, mod, yerr=mod_err, fmt="o", color="black",
                        ms=3.5, capsize=2, lw=1, zorder=3)
            ax.plot(phi, 1e3 * m["a2_eff"] * np.cos(2 * phi), "-",
                    color=C_TRUTH, lw=1.6)
            fit_amp = m["amp"] * plan.pzz_true
            ax.plot(phi, 1e3 * fit_amp * np.cos(2 * phi), "--", color=C_FIT,
                    lw=1.4)
        ax.set_xlim(0, 2 * np.pi)
        ax.set_xticks([0, np.pi, 2 * np.pi])
        ax.set_xticklabels(["0", r"$\pi$", r"$2\pi$"])
        ax.set_xlabel(r"$\phi' = \phi - \phi_S$", fontsize=9, labelpad=1)
        if k % 2 == 0:
            ax.set_ylabel((r"$T(\phi') - \hat\kappa$"
                           r"  $[\times 10^{-3}]$") if two_fill else
                          (r"$N(\phi')/\langle N\rangle - 1$"
                           r"  $[\times 10^{-3}]$"), fontsize=9)
        ax.tick_params(labelsize=8)
        ax.axhline(0.0, color="0.85", lw=0.6, zorder=0)
        ax.annotate(
            (r"$x\approx%.3g$, $Q^2\approx%.3g$ GeV$^2$" "\n"
             r"$\hat A=(%.2f\pm%.2f)\times10^{-3}$ [1 yr]; "
             r"$\pm%.2f$ [10 yr]")
            % (xs, qs, 1e3 * m["amp"], 1e3 * m["err"], 1e3 * err10),
            xy=(0.5, 1.02), xycoords="axes fraction", ha="center",
            fontsize=7.5)
        if k == 0:
            ax.annotate("injected (%s)" % PRETTY.get(args.delta_model, args.delta_model),
                        xy=(0.03, 0.93), xycoords="axes fraction",
                        color=C_TRUTH, fontsize=7)
            ax.annotate("spin-state ratio fit" if two_fill else "binned fit",
                        xy=(0.03, 0.84),
                        xycoords="axes fraction", color=C_FIT, fontsize=7)
        summary.append(
            "spot %d: x=%.3g Q2=%.3g  N_1yr=%.2e  A_truth=%+.2e  "
            "A_hat=%+.2e +- %.1e (1yr) +- %.1e (10yr)%s"
            % (k + 1, xs, qs, m["n"], m["truth"], m["amp"], m["err"],
               err10,
               ("  dA/dA_single=%.4f  sig=%.1f"
                % (m["err"] / single_fill_err(mask, lumi1_pb),
                   abs(m["truth"]) / m["err"])) if two_fill else ""))

    # --- right: amplitude vs x along the sweet-spot Q2 slice --------------
    ax = fig.add_subplot(gs[:, 2])
    q2_spot = spots[0][1]
    q2lo, q2hi = q2_spot / 1.6, q2_spot * 1.6
    xe = proj.x_edges
    pts = []
    for i0 in range(0, xe.size - 2, 2):  # merge pairs of x bins
        xc = np.sqrt(xe[i0] * xe[i0 + 2])
        if not kinematic_mask(xc, q2_spot, sampler.s):
            continue
        mask = superbin_mask(sampler, xe[i0], xe[i0 + 2], q2lo, q2hi)
        if not mask.any():
            continue
        m = do_measure(mask, lumi1_pb)
        if m["n"] < 1e3 or m["err"] > 8e-3:
            continue
        # independent 10-yr pseudo-measurement (2026-08-11 audit: reusing
        # the 1-yr draw with 10-yr bars scattered points ~sqrt(10) sigma)
        m10 = do_measure(mask, lumi10_pb)
        leak = leak_of(mask, m)
        if leak:
            m["amp"] -= leak
            m["truth"] -= leak
            m10["amp"] -= leak
        m["amp10"] = m10["amp"]
        m["err10"] = m10["err"]
        pts.append((xc, m))
    xoff = 1.045  # slight offset so the two series stay legible
    ax.errorbar([p[0] for p in pts], [1e3 * p[1]["amp"] for p in pts],
                yerr=[1e3 * p[1]["err"] for p in pts], fmt="s",
                mfc="none", color=C_FIT, ms=4, capsize=2, lw=1, zorder=3,
                label=r"1-year EIC (%g fb$^{-1}$/u)" % args.lumi_1yr)
    ax.errorbar([p[0] * xoff for p in pts],
                [1e3 * p[1]["amp10"] for p in pts],
                yerr=[1e3 * p[1]["err10"] for p in pts], fmt="o",
                color="black", ms=4, capsize=2, lw=1, zorder=4,
                label=r"10-year EIC (%g fb$^{-1}$/u)" % args.lumi_10yr)

    xg = np.logspace(np.log10(2e-4), np.log10(0.8), 250)
    q2g = np.full_like(xg, q2_spot)
    ok = kinematic_mask(xg, q2g, sampler.s)
    f2 = kern.nf2.f2a(xg, q2g) / kern.ion.A
    f1 = kern.nf2.f1a(xg, q2g) / kern.ion.A
    y = y_from_xq2(xg, q2g, sampler.s)
    curves = [(delta_func, C_TRUTH, "-",
               "%s (injected)" % PRETTY.get(args.delta_model, args.delta_model))]
    if args.delta_model != "moment_B":
        alt = dm.make("moment_B", variant=args.variant,
                      dilution=args.dilution)
        curves.append((alt, C_ALT, "-", "no-$F_1$ variant (conservative)"))
    if args.delta_model != "toy":
        curves.append((dm.make("toy", scale=1e-3), "0.45", "--",
                       r"flat toy, $\Delta/F_1=10^{-3}$"))
    for dfunc, color, ls, lab in curves:
        amp = a_cos2phi(dfunc(xg, q2g, f1), f1, f2, xg, y)
        ax.plot(xg[ok], 1e3 * amp[ok], ls, color=color, lw=1.5, label=lab)
    ax.set_xscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$A^{\cos 2\phi}$  $[\times 10^{-3}]$")
    ax.axhline(0, color="0.85", lw=0.6, zorder=0)
    ax.tick_params(labelsize=8)
    ax.set_title(r"amplitude vs $x$,  $Q^2\approx%.3g$ GeV$^2$" % q2_spot,
                 fontsize=9)
    ax.legend(fontsize=7, loc="upper left")

    fig.suptitle(
        r"Nuclear gluonometry, transversely tensor-polarized $^6$Li, "
        r"%s, %s""\n"
        r"$\Delta$: %s;  $\phi'$ pseudo-data at 1 yr (%g fb$^{-1}$/u); "
        "statistical errors only, no backgrounds or tensor radiative corrections"
        % (config.label(),
           (r"two fills $P_{zz}=%+.2f / %+.2f$ (spin-state ratio)"
            % pzz_pair) if two_fill else r"$P_{zz}=%.2f$" % plan.pzz_true,
           PRETTY.get(args.delta_model, args.delta_model)
           + (" (bag moment, dilution 1/3)" if args.delta_model == "moment_A" else ""),
           args.lumi_1yr),
        fontsize=10)
    fig.subplots_adjust(top=0.86, bottom=0.09, left=0.06, right=0.985)
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    # a non-default backend or tensor-leakage setting writes its own
    # stem, so it cannot overwrite the published massless-path figure
    tag = output_stem_tag(args)
    out = outdir / ("money_cos2phi_6Li%s.png" % ("_" + tag if tag else ""))
    fig.savefig(out, dpi=140)
    print("wrote", out)
    print("backend:", describe_backends(args, backends))
    if two_fill:
        print("estimator: two-fill spin-state ratio "
              "(reco.harmonic_ratio_fit), P_zz = %+.3g / %+.3g at equal "
              "luminosity; dA(two-fill)/dA(single-fill at P_zz = %.3g) "
              "= %.4f" % (pzz_pair[0], pzz_pair[1], args.pzz,
                          two_fill_err_ratio(args.pzz, pzz_pair)))
    print("delta model:", model.info(),
          "" if q2_ref is None else "(<Q2> = %.3g GeV^2)" % q2_ref)
    if args.tensor_gamma:
        print("tensor sector: EXACT finite-gamma kernel, b3 = %.3g b2, "
              "b4 = %.3g b2; leakage subtraction: %s"
              % (args.b3_frac, args.b4_frac, truth_leakage_route(args)))
    for line in summary:
        print(line)


if __name__ == "__main__":
    main()
