"""Polarized / tensor / gluonometric inputs.

The polarized-EMC and b1 curves are no longer toy shapes: they are the
published theory curves, digitized from the papers' own vector figures by
`tools/digitize_figure.py` and committed as CSV in `data/` (provenance in
`data/SOURCES.md`).  What is still a scenario is named as one --
`toy_delta_gluon`, and the `mode='constant'` / `mode='toy'` legacy
branches that reproduce every figure published before 2026-08-28.

  - g1p, g1n: JAM / DSSV grids (LHAPDF) -- `PartonG1` wires NNPDFpol11
  - g2: Wandzura-Wilczek from whichever g1 backend is in use (`g2_ww`,
    `ToyG1.g2_nucleus`), the one quadrature both packages share; the
    twist-3 departure from it is varied at the caller and is the residual
    systematic of the finite-gamma A_par (`asymmetries.a_parallel_exact`)
  - unpolarized EMC ratio: the A = 6 NUCLEAR PDF GRIDS (EPPS21, nNNPDF3.0)
    over CT18NLO's free isoscalar nucleon -- data-driven since 2026-08-29,
    and the common baseline the two polarized camps are transferred onto
  - polarized EMC ratio: DIGITIZED Cloet-Bentz-Thomas (PLB 642:210, Fig. 6,
    7Li at Q2 = 5 GeV2) and Tronchin-Matevosyan-Thomas (PLB 783:247, Fig. 4,
    nuclear matter at Q2 = 10 GeV2)
  - b1: DIGITIZED Cosyn-Dong-Kumano-Sargsian (PRD 95:074036, Fig. 4) and
    Miller (PRC 89:045203, Fig. 5); the 6Li transfer is still our own
    inference (plans/04 #9)
  - Delta(x,Q2): lattice-moment-normalized model (`delta_models.py`) or
    scenario values Delta/F1 in {1e-3 ... 1e-2}
"""

import numpy as np

from .structure import ToyF2, r_sigma_lt, _safe_xfx

# NumPy compat: np.trapz removed in NumPy >= 2.4, np.trapezoid absent < 2.0
_trapezoid = getattr(np, "trapezoid", None) or np.trapz


def g2_ww(g1_func, x, q2, npts=96):
    """Wandzura-Wilczek g2(x) = -g1(x) + int_x^1 du g1(u)/u at fixed Q2.

    Substitution u = x^(1-t) maps the integral to -ln(x) int_0^1 g1(x^(1-t)) dt,
    evaluated by trapezoid; accepts any mutually broadcastable (x, q2).

    This is the ONE quadrature both packages use: `polligen.xsec` re-exports
    it under its old name, and `ToyG1.g2_nucleus` below wraps it with a
    per-grid cache.  It is linear in g1, which is what makes the extraction
    of section (c) legitimate -- scaling g1 by a constant scales g2^WW by the
    same constant, so the ratio g2/g1 that the finite-gamma A_par needs is a
    property of the SHAPE and survives the multiplicative medium
    modification the polarized-EMC ratio measures.
    """
    xb, qb = np.broadcast_arrays(np.asarray(x, dtype=float),
                                 np.asarray(q2, dtype=float))
    shape = xb.shape
    xf = np.atleast_1d(xb).ravel()
    qf = np.atleast_1d(qb).ravel()
    t = np.linspace(0.0, 1.0, npts)
    xu = np.power(xf[:, None], 1.0 - t)
    qu = np.broadcast_to(qf[:, None], xu.shape)
    g1u = np.asarray(g1_func(xu, qu), dtype=float)
    integral = -np.log(np.maximum(xf, 1e-12)) * _trapezoid(g1u, t, axis=-1)
    out = -np.asarray(g1_func(xf, qf), dtype=float) + integral
    return out.reshape(shape)


class ToyG1:
    """g1 = A1 * F1 with toy A1(x) shapes (TOY; replace with JAM/DSSV).

    `r_func(x, q2) -> R = sigma_L/sigma_T` is the R used to turn F2 into
    F1; None (the default) uses this module's `r_sigma_lt` global.  This
    is the fourth R consumer, and the one the dated money scripts'
    `r_override` context manager does NOT reach: it rebinds
    `structure.r_sigma_lt` and `asymmetries.r_sigma_lt` only, so every
    g1/F1 those scripts compute still carries the toy R (code review
    2026-08-25 S1, plans/08 C2).  The hook is the fix for new code; the
    scripts stay as they are.
    """

    def __init__(self, base=None, r_func=None):
        self.base = base or ToyF2()
        self.r_func = r_func

    def a1p(self, x):
        # x^0.7 tracks moderate/high-x world data; saturates below 1
        return np.clip(np.power(np.asarray(x, dtype=float), 0.7), 0.0, 1.0)

    def a1n(self, x):
        x = np.asarray(x, dtype=float)
        # small negative at low/mid x, positive rise at high x
        return -0.07 * np.power(1.0 - x, 2.0) + 0.8 * np.power(x, 2.2)

    def _f1(self, f2, x, q2):
        r = (r_sigma_lt if self.r_func is None else self.r_func)(x, q2)
        return f2 / (2.0 * x * (1.0 + r))

    def g1p(self, x, q2):
        return self.a1p(x) * self._f1(self.base.f2p(x, q2), x, q2)

    def g1n(self, x, q2):
        return self.a1n(x) * self._f1(self.base.f2n(x, q2), x, q2)

    def g1_nucleus(self, ion, x, q2, medium_ratio=None):
        """g1A = Z P_p g1p + N P_n g1n, times an optional medium ratio.

        P_p and P_n are PER-NUCLEON effective polarizations, so the
        nucleon counts multiply them exactly as `NuclearF2.f2a` multiplies
        f2p and f2n -- the two builders now agree, and a caller dividing
        by A gets a per-nucleon g1 (plans/08 D7).  Before 2026-08-28 the
        Z and N were missing, so 6Li's per-nucleon 1/3 was diluted a
        second time by the caller's 1/A and 3He's proton term lost its
        factor 2; `beams.LI7` was rescaled to 0.866/3 and -0.037/4 in the
        same change so that every published 7Li number is bit-for-bit
        what it was (`tests/test_polarized_normalisation.py`).

        `medium_ratio(x)` is the polarized-EMC modification DR(x); the
        observable of interest is precisely deviations of DR from 1.
        """
        g1 = (ion.Z * ion.eff_pol_p * self.g1p(x, q2)
              + ion.N * ion.eff_pol_n * self.g1n(x, q2))
        if medium_ratio is not None:
            g1 = g1 * medium_ratio(x)
        return g1

    # --- g2 -------------------------------------------------------------
    #
    # A backend is asked for g2 by everything that uses the finite-gamma
    # A_par -- `polligen.xsec.InclusiveKernel` (default since 2026-08-29),
    # `fom.project_observables`, and the g2 systematic of
    # `evgen/scripts/target_mass_bound.py`.  It is Wandzura-Wilczek by
    # construction here: the twist-3 piece is the model's own uncertainty
    # and is varied at the CALLER (g2 = 0 and 1.5 g2^WW), not silently
    # inside the backend.
    #
    # The cache is what makes the term free.  `g2_ww` costs one g1
    # evaluation per quadrature node -- 96 of them -- so on a PDF grid it
    # is the whole cost of a projection, while every consumer asks for the
    # same (x, Q2) grid over and over: `money_polemc` runs three energy
    # configurations that SHARE one 40 x 30 grid, and its `main()` calls
    # `delta_dr_per_x` once per luminosity, so the six asks of a run
    # collapse to a single 96-node build.  Keyed on the ion and the raw
    # bytes of the two arrays, so it is exact rather than approximate;
    # bounded, so a script sweeping grids cannot grow it without limit.
    _G2_CACHE_MAX = 12

    def g2_nucleus(self, ion, x, q2, medium_ratio=None, npts=96):
        """g2A^WW built from THIS backend's g1A, cached per (x, Q2) grid.

        Same normalization as `g1_nucleus`: whole-nucleus, so a caller
        dividing by A gets the per-nucleon g2 that pairs with a
        per-nucleon F1.  g2^WW is linear in g1, so the division commutes
        with the quadrature."""
        x = np.asarray(x, dtype=float)
        q2 = np.asarray(q2, dtype=float)
        cache = getattr(self, "_g2_cache", None)
        if cache is None:
            cache = self._g2_cache = {}
        key = None
        if medium_ratio is None:
            # `ion` itself, NOT id(ion): beams.Ion is a frozen dataclass
            # and hashable, while an id is only unique for as long as the
            # object lives -- a transient Ion (dataclasses.replace, a test
            # fixture) can be allocated at a recycled address and would
            # collect another nucleus's table.
            key = (ion, x.shape, q2.shape, npts,
                   x.tobytes(), q2.tobytes())
            hit = cache.get(key)
            if hit is not None:
                return hit

        def g1f(xx, qq):
            return self.g1_nucleus(ion, xx, qq, medium_ratio=medium_ratio)

        out = g2_ww(g1f, x, q2, npts=npts)
        if key is not None:
            if len(cache) >= self._G2_CACHE_MAX:
                cache.pop(next(iter(cache)))
            # a hit hands back the stored array itself, so it is frozen:
            # a caller writing into a result would otherwise poison every
            # later hit.  Every caller in the repository divides by A and
            # gets a fresh array from that.
            out.setflags(write=False)
            cache[key] = out
        return out


class PartonG1(ToyG1):
    """LO g1 from a polarized LHAPDF grid via `parton` (default
    NNPDFpol11_100): g1 = (1/2) sum_q e_q^2 [Dq + Dqbar]; neutron by
    isospin. Unpolarized denominator (F1) still comes from `base`.
    Install:  yes | python3 -m parton install NNPDFpol11_100

    THREE flavours (d u s) against `structure.PartonF2`'s five (d u s c b),
    so every g1/F1 ratio built here mixes the two schemes.  That is the
    standard treatment and not a defect: NNPDFpol1.1 sets Delta c =
    Delta b = 0, so a three-flavour g1 is the whole of what it predicts,
    while the F2 a measurement sees does carry charm.  The size of the
    difference from a light-flavour-only denominator is measured and
    recorded in the `PartonF2` docstring (7.8% of F2A event-weighted over
    the money_polemc.py bins, 0.31%-0.03% across the valence window
    0.35 < x < 0.65 where the two polarized-EMC camps actually differ, but
    2.6% at x = 0.09 and 1.4% at 0.14 at Q2 = 4 GeV2, where the projection
    also plots); `tests/test_grids.py` pins the subset relation.
    """

    _E2 = {1: 1 / 9, 2: 4 / 9, 3: 1 / 9}

    def __init__(self, base=None, setname="NNPDFpol11_100", member=0,
                 r_func=None):
        # r_func is inherited for interface uniformity only: g1p/g1n come
        # straight from the polarized grid here and never divide by F1.
        super().__init__(base=base, r_func=r_func)
        from parton import mkPDF  # lazy: optional dependency
        self._pol = mkPDF(setname, member)

    def _g1_scalar(self, x, q2, swap_ud):
        if not (0.0 < x < 1.0):
            return 0.0
        tot = 0.0
        for pid, e2 in self._E2.items():
            if swap_ud and pid in (1, 2):
                e2 = self._E2[3 - pid]
            tot += e2 * (_safe_xfx(self._pol, pid, x, q2)
                         + _safe_xfx(self._pol, -pid, x, q2))
        return 0.5 * tot / x  # xfxQ2 returns x*Dq

    def g1p(self, x, q2):
        return np.vectorize(lambda a, b: self._g1_scalar(a, b, False))(
            np.asarray(x, dtype=float), np.asarray(q2, dtype=float))

    def g1n(self, x, q2):
        return np.vectorize(lambda a, b: self._g1_scalar(a, b, True))(
            np.asarray(x, dtype=float), np.asarray(q2, dtype=float))


# The pre-2026-08-29 hand-written shape, kept as mode='table' below:
# shadowing dip, anti-shadowing bump ~1.01 at x ~ 0.1, valence dip 0.88 at
# x ~ 0.7, Fermi rise beyond.  Every figure published before 2026-08-29
# that carried an unpolarized EMC ratio carried this one.
LEGACY_EMC_TABLE_X = np.array([1e-4, 0.01, 0.06, 0.10, 0.20, 0.30, 0.45,
                               0.60, 0.70, 0.80, 0.88, 0.95])
LEGACY_EMC_TABLE_R = np.array([0.96, 0.98, 1.00, 1.01, 1.00, 0.98, 0.95,
                               0.91, 0.88, 0.90, 1.00, 1.15])

UNPOL_EMC_MODES = ("epps21", "nnnpdf", "cbt", "table")
UNPOL_EMC_MODE = "epps21"        # the default of unpolarized_emc_ratio
UNPOL_EMC_Q2 = 5.0               # reference Q2 [GeV^2] of the grid modes


def unpolarized_emc_ratio(x, q2=None, mode=None):
    """Unpolarized EMC ratio R(x, Q2) for a light nucleus, per nucleon.

    mode='epps21' (the default, `UNPOL_EMC_MODE`) is DATA-DRIVEN: the
    per-nucleon F2 of EPPS21nlo_CT18Anlo_Li6 over the free isoscalar
    nucleon (F2p + F2n)/2 of CT18ANLO -- EPPS21's OWN proton baseline, so
    the fit cancels and the ratio is the nuclear modification alone (it
    was CT18NLO until 2026-08-29, 4.2% shallower) -- both through
    `structure.NuclearF2Ratio` and so through the same five-flavour charge
    sum, at `q2` or at the reference `UNPOL_EMC_Q2` = 5 GeV2 (CBT's scale,
    so that the two-camp comparison below is built at one scale).  This
    closes plans/02 step 1.2.1.

    mode='nnnpdf' is the same construction on nNNPDF30_nlo_as_0118_A6_Z3.
    Its A = 6 EMC effect is much the shallower of the two -- 0.0137 mean
    valence depletion against EPPS21's 0.0311 -- because nNNPDF3.0 has
    almost no A = 6 data to constrain it there; the spread between the two
    fits is the leading uncertainty on everything built on this baseline
    and is reported as such, not averaged away.

    mode='cbt' is the digitized CBT 7Li unpolarized curve (the right
    ISOTOPE, but a model rather than a fit to data: mean valence depletion
    0.0584, 1.9x EPPS21's).  mode='table' is the legacy 12-point shape
    (0.0646), which every figure published before 2026-08-29 used and
    which the `mode='constant'` branches below and the dated
    `money_delta_*.py` notes stay pinned to.

    ISOTOPE CAVEAT.  LHAPDF has A = 6 grids and nothing for A = 7, while
    the polarized-EMC observable of `money_polemc.py` is 7Li.  Over
    0.1 < x < 0.7 the EPPS21 Li-6 ratio and the CBT 7Li curve differ by up
    to 0.056 (mean 0.022), but that difference is a fit against a model,
    not 6Li against 7Li: the EMC slope grows roughly as ln A, so 7Li's
    depletion should exceed 6Li's by about ln7/ln6 = 1.086, i.e. the
    isotope alone moves the valence depletion from 0.0298 to about 0.032
    -- a twelfth of the gap to CBT and a fifth of the gap to nNNPDF3.0.
    The isotope is therefore the smallest of the three uncertainties here
    (plans/04 #8).

    `q2` is ignored by the two digitized/tabulated modes, which are curves
    at one scale (CBT at Q2 = 5, the legacy table at none); passing it is
    not an error, because callers switch modes without switching
    signatures.
    """
    mode = UNPOL_EMC_MODE if mode is None else mode
    if mode == "table":
        return np.interp(np.asarray(x, dtype=float),
                         LEGACY_EMC_TABLE_X, LEGACY_EMC_TABLE_R)
    if mode == "cbt":
        return cbt_unpolarized_emc_ratio(x)
    if mode not in ("epps21", "nnnpdf"):
        raise ValueError("mode must be one of %s" % (UNPOL_EMC_MODES,))
    q2 = UNPOL_EMC_Q2 if q2 is None else q2
    return _nuclear_ratio(mode)(x, q2)


def legacy_emc_ratio(x):
    """`unpolarized_emc_ratio` frozen at mode='table', as a bare r(x).

    The `NuclearF2(emc_ratio=)` hook takes one argument, and the three
    scripts that use it -- `money_delta_realistic.py` and the dated
    `money_delta_20260720/21.py` -- are reproductions of dated notes
    (`notes/money_delta_note_2026-07-16.md` R3), so they must keep the
    12-point shape they were written against when the default moved to
    the nuclear grids on 2026-08-29."""
    return unpolarized_emc_ratio(x, mode="table")


_NUCLEAR_RATIOS = {}


def _nuclear_ratio(mode):
    """Cached `structure.NuclearF2Ratio`; opening a grid is not cheap."""
    if mode not in _NUCLEAR_RATIOS:
        from .structure import NuclearF2Ratio
        _NUCLEAR_RATIOS[mode] = NuclearF2Ratio(mode)
    return _NUCLEAR_RATIOS[mode]


# --------------------------------------------------------------------------
# Digitized theory curves (plans/02 step 1.2.2 / 1.2.3, closed 2026-08-28)
# --------------------------------------------------------------------------

_CURVES = {}


def _load_curve(name):
    """Columns of a digitized CSV in `polli_fastsim/data`, cached.

    The tables are committed data, not a generated artefact: they are read
    with `importlib.resources` so they travel with the package however it
    is put on the path (the fast sim is used from a source checkout via
    `sys.path`, so there is no installer and no `package_data` to keep in
    step).  Column names come from the last comment line of the file, the
    provenance from the ones above it and from `data/SOURCES.md`."""
    if name not in _CURVES:
        try:
            from importlib.resources import files
            text = (files(__package__) / "data" / (name + ".csv")).read_text()
        except Exception:                      # pragma: no cover - fallback
            import os
            here = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "data", name + ".csv")
            with open(here, encoding="utf-8") as fh:
                text = fh.read()
        header, rows = None, []
        for line in text.splitlines():
            if line.startswith("#"):
                header = line
            elif line.strip():
                rows.append([float(v) for v in line.split(",")])
        cols = [c.strip() for c in header.lstrip("# ").split(",")]
        arr = np.asarray(rows, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != len(cols):
            raise ValueError("malformed digitized table %s" % name)
        if np.any(np.diff(arr[:, 0]) <= 0):
            raise ValueError("non-monotone x in digitized table %s" % name)
        _CURVES[name] = dict(zip(cols, arr.T))
    return _CURVES[name]


def curve_x_range(name):
    """(x_min, x_max) actually covered by a digitized table."""
    xs = _load_curve(name)["x"]
    return float(xs[0]), float(xs[-1])


def _interp(name, column, x):
    """Table lookup, CONSTANT-EXTRAPOLATED outside the digitized range.

    np.interp clamps to the end values; every caller's docstring says
    where its table stops and what that means, because the money plots
    sample x from 0.003-0.005 while the published figures start at
    0.0015-0.028."""
    t = _load_curve(name)
    return np.interp(np.asarray(x, dtype=float), t["x"], t[column])


def _interp_tapered(name, column, x, x_end=1.0, power=3):
    """Table lookup for a STRUCTURE FUNCTION: constant below the table,
    tapered as ((x_end - x)/(x_end - x_max))**power above it.

    A constant extrapolation to x -> 1 is wrong for anything that must
    vanish at the elastic edge: the kernel forms b1/F1, and with b1 frozen
    at its x = 0.9 value while F1 falls as (1 - x)^3 the ratio diverges by
    the generator grid's x = 0.955 cell (found 2026-08-28 in
    money_tagged_azz.py, whose positivity guard raised).  The taper uses
    the same counting-rule power as F1, so b1/F1 stays frozen at its value
    at the table's end instead of growing.  Below x_min the constant
    extrapolation is kept (see _interp)."""
    t = _load_curve(name)
    xa = np.asarray(x, dtype=float)
    x_max = float(t["x"][-1])
    val = np.interp(xa, t["x"], t[column])
    if x_max < x_end:
        frac = np.clip((x_end - xa) / (x_end - x_max), 0.0, 1.0)
        val = np.where(xa > x_max, val * frac ** power, val)
    return val


CBT_TABLE = "cbt_polemc_7Li_Q5"      # Cloet-Bentz-Thomas Fig. 6, 7Li, Q2 = 5
TMT_TABLE = "tmt_polemc_nm_Q10"      # Tronchin-Matevosyan-Thomas Fig. 4, NM
MILLER_TABLE = "b1_miller"           # Miller Fig. 5, b1(pion) + b1(6q)
CDKS_TABLE = "b1_cdks_q2p5"          # Cosyn-Dong-Kumano-Sargsian Fig. 4

# The two camps are computed for DIFFERENT targets at DIFFERENT scales --
# CBT for 7Li at Q2 = 5 GeV2, TMT for isospin-symmetric nuclear matter at
# Q2 = 10 GeV2 -- so subtracting their published R values would compare
# two nuclei, not two models.  What each model actually predicts is the
# STRENGTH of the polarized EMC effect relative to the unpolarized one in
# its own calculation ("about twice" vs "about equal"), and that ratio is
# what transfers.  Both camps are therefore put on one common unpolarized
# baseline through
#
#     DR_model(x) = 1 - s_model * [1 - R_pol,model(x)],
#     s_model     = <1 - R_unpol,baseline> / <1 - R_unpol,model>  over the
#                   valence window 0.35 < x < 0.65.
#
# WHICH BASELINE, changed 2026-08-29 (plans/02 step 1.2.1).  It used to be
# CBT's own digitized 7Li unpolarized curve, so s was 1 for CBT by
# construction and 0.397 for TMT, whose nuclear matter is depleted 0.1470
# against 7Li's 0.0583 over the window.  It is now the DATA-DRIVEN
# `unpolarized_emc_ratio` default, EPPS21's Li-6 over the free isoscalar
# nucleon, whose valence depletion is 0.03105 -- just over half of CBT's
# model value.  Both camps are therefore scaled: s = 0.5322 for CBT and 0.2113
# for TMT, their RATIO unchanged at 2.52 because that ratio is a property
# of the two calculations and not of the baseline.
#
# WHAT THAT COSTS THE HEADLINE, stated plainly because it is the single
# largest number this module moved: the separation the money plot draws is
# proportional to the baseline depletion, so it halves, from 0.0395 to
# 0.0202 at x = 0.36 and from 0.0110 to 0.0056 at 0.65.  The reach halves
# with it.  That is not a loss of information but the removal of a
# borrowed one: the old figure asserted CBT's own model value for the
# unpolarized EMC effect of a nucleus nobody has measured, and EPPS21's
# fit says half of it.  `baseline='cbt'` reproduces the old figure
# exactly, `baseline='nnnpdf'` gives 0.213 and 0.0846 (nNNPDF3.0's A = 6
# EMC effect is shallower again), and the spread between those three is
# the honest uncertainty on the reach -- larger than any of the detector
# or luminosity terms in the same projection.  EPPS21's own 90%-CL
# Hessian band on the valence depletion is 0.0297 +0.0393 -0.0406, i.e.
# compatible with no EMC effect at A = 6 at all, so the spread is not a
# choice between three well-determined numbers.
#
# The strength ratio is taken over the valence window rather than
# pointwise for a physical reason, not a numerical one: pointwise,
# r(x) = (1 - R_pol)/(1 - R_unpol) is singular wherever a model's
# unpolarized ratio crosses 1 -- at x = 0.280 and 0.840 for CBT's 7Li
# curve and at x = 0.721 for TMT's nuclear matter -- and in the
# shadowing / anti-shadowing region below x ~ 0.3 the "scale the medium
# effect by the unpolarized one" logic has no content anyway, because the
# unpolarized effect there is an ENHANCEMENT while both models keep a 7%
# polarized depletion.  `cbt_ratio_of_effects` / `tmt_ratio_of_effects`
# expose the pointwise ratio for the valence region, where it is the
# number the "2x" and "1x" slogans refer to.
#
# WHAT THIS COSTS, stated because a projection reads it straight off the
# figure: s is CONSTANT in x, so it is applied everywhere the money plot
# draws, including two decades below the window in which it is defined.
# The two PUBLISHED polarized curves agree to better than 0.007 over
# 0.028 < x < 0.30 (0.0016 at x = 0.09, 0.0006 at x = 0.14), so the 0.020
# separation the transferred pair shows there is the baseline rescaling of
# two different targets and not a disagreement between the papers.
# Inside the window, where the comparison does mean something, the
# separation is 0.0202 at x = 0.36 falling to 0.0056 at 0.65 as CBT's
# ratio of effects falls towards TMT's.  `tmt_published_emc_ratio`
# returns the untransferred curve so that both separations can be quoted
# side by side, `money_polemc.py` prints them side by side, and Report 0
# section 5.3 says so in words.
POLEMC_VALENCE_WINDOW = (0.35, 0.65)


POLEMC_BASELINE = "epps21"        # default baseline of the transfer
POLEMC_BASELINE_MODES = UNPOL_EMC_MODES


def valence_depletion(mode=None, q2=None, n=301):
    """<1 - R_unpol> of one baseline over POLEMC_VALENCE_WINDOW.

    0.03105 for 'epps21' and 0.01372 for 'nnnpdf' at the reference
    Q2 = 5 GeV2, 0.05835 for 'cbt' and 0.06459 for 'table'.  The two grid
    values are almost Q2-independent over the range the projection covers
    (EPPS21: 0.03108 at Q2 = 4, 0.03090 at 10, 0.02952 at 100).
    """
    lo, hi = POLEMC_VALENCE_WINDOW
    g = np.linspace(lo, hi, n)
    return float((1.0 - unpolarized_emc_ratio(g, q2=q2, mode=mode)).mean())


_VALENCE_SCALES = {}


def valence_scale(table, baseline=None, q2=None):
    """<1 - R_unpol,baseline> / <1 - R_unpol,table> over the valence window.

    The strength with which one camp's POLARIZED depletion is carried onto
    the common unpolarized baseline.  On the default 'epps21' baseline it
    is 0.5322 for CBT and 0.2113 for TMT; on the legacy 'cbt' baseline it
    is exactly 1 for CBT (which computed 7Li itself) and 0.397009 for TMT,
    which is what every figure published before 2026-08-29 used.
    """
    baseline = POLEMC_BASELINE if baseline is None else baseline
    key = (table, baseline, q2)
    if key not in _VALENCE_SCALES:
        lo, hi = POLEMC_VALENCE_WINDOW
        g = np.linspace(lo, hi, 301)
        d_ref = valence_depletion(mode=baseline, q2=q2)
        d_mod = float((1.0 - _interp(table, "R_unpol", g)).mean())
        _VALENCE_SCALES[key] = d_ref / d_mod
    return _VALENCE_SCALES[key]


# The legacy pair, evaluated at import because neither needs a grid.  They
# are the transfer on the 'cbt' baseline, i.e. the published-before-
# 2026-08-29 numbers, and are what `baseline='cbt'` reproduces.
CBT_VALENCE_SCALE_ON_CBT = 1.0   # CBT computed 7Li itself -- exact identity
_VW = np.linspace(POLEMC_VALENCE_WINDOW[0], POLEMC_VALENCE_WINDOW[1], 301)
TMT_VALENCE_SCALE_ON_CBT = float(
    (1.0 - _interp(CBT_TABLE, "R_unpol", _VW)).mean()
    / (1.0 - _interp(TMT_TABLE, "R_unpol", _VW)).mean())     # = 0.397009


def cbt_unpolarized_emc_ratio(x):
    """Digitized 7Li UNPOLARIZED EMC ratio (CBT Fig. 6, blue dashed).

    The common baseline of the two-camp comparison.  Covers
    x = 0.028-0.871; constant-extrapolated outside, so the value below
    x = 0.028 is frozen at 1.006 and the Fermi rise above 0.871 is not
    followed."""
    return _interp(CBT_TABLE, "R_unpol", x)


def cbt_polarized_emc_ratio(x, mode="digitized", eq=23, baseline=None):
    """Cloet-Bentz-Thomas polarized EMC ratio for 7Li, on the baseline.

    mode='digitized' (default) takes the published curve of PLB 642:210
    Fig. 6 -- `eq=23` the red dotted R^{3/2 3/2}_{As} of their Eq. (23),
    which is what plans/01 defines the programme's DR_A to be, and
    `eq=26` the red solid R^{(3/2 1)}_{As} of their Eq. (26) -- and
    carries its depletion onto the common unpolarized baseline with the
    valence factor `valence_scale(CBT_TABLE, baseline)`.

    Until 2026-08-29 the baseline WAS CBT's own 7Li unpolarized curve, so
    that factor was exactly 1 and this returned the published curve
    untouched; `baseline='cbt'` still does.  On the default data-driven
    'epps21' baseline it is 0.5106, because EPPS21's Li-6 valence
    depletion (0.0298) is half of CBT's model 7Li one (0.0584): what
    transfers is the model's RATIO of polarized to unpolarized effect,
    and the size of the unpolarized effect is then taken from the fit.

    mode='constant' is the pre-2026-08-28 stand-in, 1 - 2(1 - R_EMC) on
    the legacy 12-point table, kept so the figures published before the
    digitization can be reproduced; it ignores `baseline`.

    DIGITIZED RANGE x = 0.028-0.871.  Below 0.028 the value is frozen at
    the endpoint 0.919 and above 0.871 at 1.007; `money_polemc.py` draws
    from x = 0.005, so the leftmost decade of that figure is a flat
    continuation, not a prediction."""
    if mode == "constant":
        return 1.0 - 2.0 * (1.0 - unpolarized_emc_ratio(x, mode="table"))
    if mode != "digitized":
        raise ValueError("mode must be 'digitized' or 'constant'")
    if eq not in (23, 26):
        raise ValueError("eq must be 23 (R^{3/2 3/2}) or 26 (R^{(3/2 1)})")
    s = valence_scale(CBT_TABLE, baseline)
    return 1.0 - s * (1.0 - _interp(CBT_TABLE, "R_pol_eq%d" % eq, x))


def tmt_polarized_emc_ratio(x, mode="digitized", baseline=None):
    """Tronchin-Matevosyan-Thomas polarized EMC ratio, on the baseline.

    mode='digitized' (default) rescales the published nuclear-matter
    polarized depletion of PLB 783:247 Fig. 4 with the valence factor
    `valence_scale(TMT_TABLE, baseline)` -- 0.2113 on the default
    data-driven 'epps21' baseline, 0.397009 on the legacy 'cbt' one --
    so that it is comparable with `cbt_polarized_emc_ratio` bin by bin.
    mode='constant' is the pre-2026-08-28 stand-in, the legacy 12-point
    table itself, and ignores `baseline`.

    DIGITIZED RANGE x = 0.0015-0.739 (the published curve stops there);
    above 0.739 the value is frozen at the endpoint, which also freezes
    the nuclear-matter Fermi rise at a different x than 7Li's, so the
    x > 0.74 end of the money plot is not a like-for-like comparison."""
    if mode == "constant":
        return unpolarized_emc_ratio(x, mode="table")
    if mode != "digitized":
        raise ValueError("mode must be 'digitized' or 'constant'")
    s = valence_scale(TMT_TABLE, baseline)
    return 1.0 - s * (1.0 - _interp(TMT_TABLE, "R_pol", x))


def cbt_published_emc_ratio(x, eq=23):
    """CBT's PUBLISHED 7Li polarized ratio, BEFORE the baseline transfer.

    Until 2026-08-29 this was `cbt_polarized_emc_ratio` itself, because
    the baseline was CBT's own unpolarized curve and the transfer factor
    was 1.  It no longer is, so the two accessors have parted, and the
    question "how far apart are the two CALCULATIONS?" is answered by this
    one against `tmt_published_emc_ratio` -- not by the transferred pair,
    which answers "how far apart does the projection draw them on one
    common baseline?".

    DIGITIZED RANGE x = 0.028-0.871, constant-extrapolated."""
    if eq not in (23, 26):
        raise ValueError("eq must be 23 (R^{3/2 3/2}) or 26 (R^{(3/2 1)})")
    return _interp(CBT_TABLE, "R_pol_eq%d" % eq, x)


def tmt_published_emc_ratio(x):
    """TMT's PUBLISHED nuclear-matter polarized ratio, BEFORE the transfer.

    `tmt_polarized_emc_ratio` returns the same curve rescaled to 7Li
    valence strength, which is what a projection against CBT's 7Li needs;
    this one is what PLB 783:247 Fig. 4 actually draws.  The difference
    matters when the question is how far apart the two CALCULATIONS are
    rather than how far apart the projection plots them, and the answer is
    strongly x-dependent: over 0.028 < x < 0.30 the two published
    polarized curves agree to better than 0.007 (0.0016 at x = 0.09,
    0.0006 at x = 0.14), while at x = 0.45, 0.55 and 0.65 they differ by
    0.038, 0.097 and 0.115 -- but most of THAT is nuclear matter against
    7Li rather than one model against the other, which is the whole reason
    for the transfer.  The comparison is meaningful only inside
    POLEMC_VALENCE_WINDOW, the only window in which the transfer factor is
    defined: the separation the transferred pair shows below x ~ 0.3 is the
    strength rescaling extrapolated out of that window, not a disagreement
    between the papers, and inside it the transferred pair separates by
    0.0202 at x = 0.36 and 0.0056 at 0.65 on the default EPPS21 baseline
    (0.0395 and 0.0110 on the legacy CBT one) as CBT's ratio of effects
    falls towards TMT's.  `money_polemc.py` prints both separations for
    this reason, and Report 0 section 5.3 says it in words.

    DIGITIZED RANGE x = 0.0015-0.739, constant-extrapolated."""
    return _interp(TMT_TABLE, "R_pol", x)


def cbt_ratio_of_effects(x, eq=23):
    """Pointwise (1 - R_pol)/(1 - R_unpol) of CBT's own 7Li figure.

    2.25 / 1.69 / 1.41 / 1.14 at x = 0.40 / 0.45 / 0.50 / 0.60 -- the "about
    twice the unpolarized EMC effect" of PLB 642:210 holds at x ~ 0.4 and
    decays to a minimum of 1.06 at x ~ 0.70 without reaching 1 (the eq=26
    curve does cross the unpolarized one, at x = 0.651).  DIVERGES at the
    x = 0.280 and 0.840 crossings of the unpolarized curve; the caller is
    responsible for staying in the valence region."""
    d_unpol = 1.0 - _interp(CBT_TABLE, "R_unpol", x)
    return (1.0 - _interp(CBT_TABLE, "R_pol_eq%d" % eq, x)) / d_unpol


def tmt_ratio_of_effects(x):
    """Pointwise (1 - R_pol)/(1 - R_unpol) of TMT's nuclear-matter figure:
    1.01 / 0.98 / 1.00 / 1.08 at x = 0.40 / 0.45 / 0.50 / 0.60 -- the
    "polarized ~ unpolarized" of PLB 783:247.  Diverges at the x = 0.721
    crossing of their unpolarized curve."""
    return ((1.0 - _interp(TMT_TABLE, "R_pol", x))
            / (1.0 - _interp(TMT_TABLE, "R_unpol", x)))


def toy_polarized_emc_ratio(x, **kw):
    """Backward-compatible alias for the CBT curve."""
    return cbt_polarized_emc_ratio(x, **kw)


# --------------------------------------------------------------------------
# b1
# --------------------------------------------------------------------------
#
# CHANGE OF MEANING, 2026-08-28.  Before the digitization both b1 hooks
# returned a SHAPE times the F1 they were handed, so their magnitude was
# set by the caller's structure function.  They now return the published
# ABSOLUTE b1 of the deuteron and ignore `f1`, which is kept in the
# signature because `fom.project_observables` and `xsec.InclusiveKernel`
# pass it positionally.  The published figures are b1^d in the
# Hoodbhoy-Jaffe-Manohar per-deuteron normalization (CDKS plot x*b1 out to
# x = 1.6, i.e. nucleon-scaled x); every consumer here pairs b1 with a
# PER-NUCLEON F1, so the tables are halved on the way out -- see
# B1_PER_DEUTERON_TO_PER_NUCLEON.  `mode='toy'` restores the old shapes.
#
# WHOSE b1 IS IT?  These functions return the DEUTERON's, and the two
# halves of the repository do different things with it.  `money_b1.py`
# wraps them in `b1_li6_from_deuteron`, so its 6Li signal carries the
# rank-2 transfer and the 2/6, a factor 0.307.  The twelve evgen scripts
# that hand `toy_b1` straight to `xsec.InclusiveKernel(beams.LI6, ...)`
# do not, so their 6Li kernel carries the deuteron's b1 undiluted -- a
# factor 3.3 too large for 6Li.  Before the digitization both were an
# explicit scenario SHAPE and the mismatch was only labelling; now it is
# numerical, and it is tolerated only because b1 enters the kernel through
# w_avg alone (never through the cos 2phi amplitude a_2, which depends on
# `delta`), where it is O(1e-3) of the cross section, so 3.3x of it is
# below those figures' Monte-Carlo noise.  The proper fix is a transfer
# argument on InclusiveKernel for spin-1 nuclei other than the deuteron;
# until then, read a 6Li b1 off `money_b1.py` and not off the generator.
B1_PER_DEUTERON_TO_PER_NUCLEON = 0.5


# --------------------------------------------------------------------------
# The Q2 handle on b1 (plans/02 step 1.2.1, wired in 2026-08-29)
# --------------------------------------------------------------------------
#
# Both b1 papers publish a Q2 SET alongside the curve we digitize as the
# nominal: Miller's Fig. 6 at Q2 = 1.17, 1.76, 2.12 and 3.25 GeV2, CDKS's
# Fig. 5 at 1.0, 2.5 and 5.0.  Neither set covers the nominal curve's full
# x range, so the Q2 dependence is carried as a RATIO to the set's own
# member at the nominal scale, which cancels the digitization of the
# common shape and leaves only the Q2 lever.
#
# WHICH SCALE THE NOMINAL CURVE IS AT is not printed in Miller's Fig. 5
# caption; it is READ OFF the two digitizations, which is possible because
# Fig. 6's Q2 = 3.25 curve and Fig. 5's total coincide: 1.001 / 1.000 /
# 1.000 / 1.001 / 1.003 at x = 0.1007 / 0.15 / 0.20 / 0.40 / 0.50, and
# within 1.3% at every point of the overlap except the x ~ 0.3 zero
# crossing.  So `b1_miller.csv` is Miller at Q2 = 3.25 GeV2, and the
# hidden-colour six-quark piece Fig. 5 adds to Fig. 6's pion is invisible
# on that scale.  CDKS name theirs: `b1_cdks_q2p5.csv` is Q2 = 2.5.
#
# WHAT IT COSTS.  The scale is frozen in x outside the set's range (Miller
# 0.1007-0.7000, CDKS 0.0050-1.5850) and in ln Q2 outside the set's nodes,
# so every bin above Q2 = 3.25 gets Miller's 3.25 curve -- i.e. the
# published nominal, ratio exactly 1 -- and the lever bites only on the
# low-Q2 bins, where it reaches 0.826 at x = 0.15 and Q2 = 1.17.  Below
# x = 0.1007 the ratio is frozen at its value there, 0.889 at Q2 = 1.17,
# and that is where most of the b1 reach sits, so the low-x bins carry a
# Q2 correction inherited from x = 0.1 rather than one Miller plots.  Near
# a zero crossing the ratio is not a useful object; it is clipped to
# `B1_Q2_SCALE_CLIP`, which for Miller bites only inside the two windows
# 0.2991 < x < 0.3027 and 0.5753 < x < 0.5789 around his two zero
# crossings, where |b1| stays below 5.9e-4 and 3.7e-4 per deuteron and
# |A_zz| is at the 1e-6 level.
#
# ONLY MILLER'S SET IS WIRED IN.  CDKS's Fig. 5 exists as
# `b1_cdks_q2set.csv` and is not used, because their theory-1 curve
# crosses zero twice (x = 0.06 and 0.42) and the crossings MOVE with Q2,
# so the ratio swings over the whole clip band across 0.1 < x < 0.7 and
# would be a digitization artefact rather than a Q2 lever.  The CDKS camp
# is drawn at its published Q2 = 2.5 GeV2 everywhere, which costs nothing
# at the digits that matter: it is below 0.2 sigma in every bin.
MILLER_Q2SET_TABLE = "b1_miller_q2set"
MILLER_Q2SET_COLS = ("b1_q2_1p17", "b1_q2_1p76", "b1_q2_2p12", "b1_q2_3p25")
MILLER_Q2SET_NODES = (1.17, 1.76, 2.12, 3.25)
MILLER_B1_Q2_REF = 3.25          # the scale of b1_miller.csv, read off above
B1_Q2_SCALE_CLIP = (0.5, 2.0)


def _q2_scale(table, cols, nodes, q2_ref, x, q2, clip=B1_Q2_SCALE_CLIP):
    """b1(x, Q2) / b1(x, q2_ref) from a digitized Q2 set.

    Linear in ln Q2 between the set's nodes and CONSTANT outside them;
    constant in x outside the set's own range (`_interp`'s convention);
    clipped to `clip` so that the model's zero crossings cannot turn the
    ratio into a pole."""
    t = _load_curve(table)
    xa, qa = np.broadcast_arrays(np.asarray(x, dtype=float),
                                 np.asarray(q2, dtype=float))
    shape = xa.shape
    xf = xa.reshape(-1)
    lq = np.log(np.maximum(qa.reshape(-1), 1e-6))
    vals = np.stack([np.interp(xf, t["x"], t[c]) for c in cols])
    ln_nodes = np.log(np.asarray(nodes, dtype=float))
    cell = np.arange(vals.shape[1])

    def at(ln):
        j = np.clip(np.searchsorted(ln_nodes, ln) - 1, 0, len(nodes) - 2)
        w = np.clip((ln - ln_nodes[j]) / (ln_nodes[j + 1] - ln_nodes[j]),
                    0.0, 1.0)
        v0, v1 = vals[j, cell], vals[j + 1, cell]
        return v0 + w * (v1 - v0)

    num = at(lq)
    den = at(np.full(lq.shape, np.log(float(q2_ref))))
    safe = np.abs(den) > 0
    r = np.where(safe, num / np.where(safe, den, 1.0), 1.0)
    r = np.clip(np.nan_to_num(r, nan=1.0), clip[0], clip[1])
    return r.reshape(shape)


def miller_b1_q2_scale(x, q2):
    """Miller's Fig.-6 Q2 lever on the Fig.-5 total, relative to Q2 = 3.25."""
    return _q2_scale(MILLER_Q2SET_TABLE, MILLER_Q2SET_COLS,
                     MILLER_Q2SET_NODES, MILLER_B1_Q2_REF, x, q2)


def _toy_b1_shape(x, q2, f1):
    """The pre-2026-08-28 'HERMES-like' scenario shape times F1."""
    x = np.asarray(x, dtype=float)
    shape = 0.01 * np.power(np.maximum(x, 1e-6), -0.2) * (1.0 - x / 0.20)
    return shape * np.exp(-3.0 * x) * f1


def toy_b1(x, q2, f1, mode="digitized", q2_evolve=False):
    """b1 of the deuteron, per nucleon: Miller's pion + hidden-colour total.

    mode='digitized' (default) is PRC 89:045203 Fig. 5, the curve that
    reproduces the HERMES measurement with a hidden-colour probability of
    0.15% -- the 'HERMES-like' camp.  mode='toy' is the old analytic shape
    times F1; `f1` is ignored otherwise.

    `q2_evolve=True` multiplies by `miller_b1_q2_scale(x, q2)`, Miller's
    own Fig.-6 Q2 set relative to the Q2 = 3.25 GeV2 of Fig. 5 (see the
    block above), so the curve is evaluated at the caller's scale rather
    than at one slice.  It stays OFF by default because the twelve evgen
    scripts that hand this function to `xsec.InclusiveKernel` want one
    fixed b1 per sample and already carry a much larger normalisation
    caveat (the 6Li transfer below); `money_b1.py --signal-q2 binned`, the
    projection that reads a number off b1, turns it on.  Above
    Q2 = 3.25 GeV2 the scale is exactly 1, so a bin at the old Q2 = 4
    slice is bit-for-bit what it was.

    DIGITIZED RANGE x = 0.010-0.900: b1 is frozen at 0.0647 per nucleon
    below x = 0.01, where the generator does sample, and above x = 0.900 it
    is tapered to zero at x = 1 with the (1 - x)^3 falloff of F1
    (_interp_tapered), so that b1/F1 -- what the kernel forms -- stays at
    its x = 0.9 value instead of diverging at the generator's x = 0.955
    cell."""
    if mode == "toy":
        return _toy_b1_shape(x, q2, f1)
    if mode != "digitized":
        raise ValueError("mode must be 'digitized' or 'toy'")
    b1 = B1_PER_DEUTERON_TO_PER_NUCLEON * _interp_tapered(MILLER_TABLE, "b1", x)
    if q2_evolve:
        b1 = b1 * miller_b1_q2_scale(x, q2)
    return b1


def b1_convolution(x, q2, f1, mode="digitized"):
    """b1 of the deuteron, per nucleon: the standard convolution.

    mode='digitized' (default) is the SD + DD sum of Cosyn-Dong-Kumano-
    Sargsian PRD 95:074036 Fig. 4, theory 1 (their Eq. 16), at
    Q2 = 2.5 GeV2, converted from the plotted x*b1.  It is the camp that
    finds |b1| < 1e-3 at x >~ 0.2 -- an order of magnitude below the data --
    and it is NOT simply a tenth of the HERMES-like curve: the old
    `0.1 * toy_b1` had the same shape and one zero crossing, while the real
    convolution changes sign at x = 0.06 and again at x = 0.42.
    mode='toy' restores that tenth.

    DIGITIZED RANGE x = 0.010-1.590, constant-extrapolated below 0.01."""
    if mode == "toy":
        return 0.1 * _toy_b1_shape(x, q2, f1)
    if mode != "digitized":
        raise ValueError("mode must be 'digitized' or 'toy'")
    x = np.maximum(np.asarray(x, dtype=float), 1e-6)
    xb1 = _interp(CDKS_TABLE, "xb1_theory1_sum", x)
    return B1_PER_DEUTERON_TO_PER_NUCLEON * xb1 / x


def close_kumano_integral(mode="cdks"):
    """Integral of b1 over the digitized range -- the Close-Kumano sum rule
    int b1 dx = 0, which holds when the sea is not tensor polarized.

    Reported, not enforced: the tables cover only what the papers plot
    (CDKS x = 0.01-1.59, Miller 0.01-0.90), so a nonzero value is as much a
    statement about the missing ends as about the model."""
    if mode == "cdks":
        t = _load_curve(CDKS_TABLE)
        return float(np.trapz(t["xb1_theory1_sum"] / t["x"], t["x"]))
    t = _load_curve(MILLER_TABLE)
    return float(np.trapz(t["b1"], t["x"]))


# The rank-2 transfer of the embedded deuteron's tensor polarization to
# 6Li, from the tagged two-cluster spin model:
# `polligen.tagged.TaggedModel(li6_alpha_channel()).tensor_dilution()` =
# 1 - (9/10) P_D at P_D = 0.0867 (pinned in evgen/tests/test_tagged.py).
# The 0.87 it replaces is the VECTOR dilution 1 - (3/2) P_D -- the wrong
# rank for b1, and the whole of plans/08 D9.
LI6_B1_RANK2_TRANSFER = 0.921947
LI6_B1_LEGACY_TRANSFER = 0.87
# 6Li carries two polarized nucleons out of six.  The factor puts the
# signal on the same per-nucleon footing as the F1 it is divided by and as
# the Delta sector's `delta_models.py` dilution=1/3, which is what D9 asks
# for: before 2026-08-28 the b1 money plot drew a signal with 0.87 and no
# 2/6 against errors with neither.
LI6_B1_PER_NUCLEON = 2.0 / 6.0


def b1_li6_from_deuteron(b1_d, transfer=LI6_B1_RANK2_TRANSFER,
                         per_nucleon=LI6_B1_PER_NUCLEON):
    """Whole-nucleus 6Li b1 from the embedded deuteron's, per nucleon.

    b1(6Li)/nucleon = transfer * (2/6) * b1(d)/nucleon.  No published 6Li
    b1 exists (plans/04 #9), so the embedded-deuteron scaling is our own
    inference whatever b1(d) is used.  `transfer=LI6_B1_LEGACY_TRANSFER,
    per_nucleon=1.0` reproduces the pre-2026-08-28 curve exactly."""
    return transfer * per_nucleon * b1_d


def toy_delta_gluon(x, q2, f1, scale=1e-3):
    """TOY double-helicity-flip Delta(x,Q2) = scale * F1 * (1-x)^4 * x^0.3.

    Delta is alpha_s-suppressed and unknown; `scale` sets Delta/F1 at its
    maximum. Use scale in {1e-3, 3e-3, 1e-2} as discovery scenarios, to be
    replaced by lattice-moment-anchored models.
    """
    x = np.asarray(x, dtype=float)
    return scale * f1 * np.power(x, 0.3) * np.power(1.0 - x, 4.0)
