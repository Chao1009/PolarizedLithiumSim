"""Polarized / tensor / gluonometric inputs.

The polarized-EMC and b1 curves are no longer toy shapes: they are the
published theory curves, digitized from the papers' own vector figures by
`tools/digitize_figure.py` and committed as CSV in `data/` (provenance in
`data/SOURCES.md`).  What is still a scenario is named as one --
`toy_delta_gluon`, and the `mode='constant'` / `mode='toy'` legacy
branches that reproduce every figure published before 2026-08-28.

  - g1p, g1n: JAM / DSSV grids (LHAPDF) -- `PartonG1` wires NNPDFpol11
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


def unpolarized_emc_ratio(x):
    """Qualitative unpolarized EMC ratio R_EMC(x) for a light nucleus:
    shadowing dip, anti-shadowing bump ~1.01 at x~0.1, valence dip ~0.88
    at x~0.7, Fermi rise beyond. Smooth interpolation of the canonical
    shape (SCENARIO; EPPS21 / data fits are its own replacement, plans/02
    step 1.2.1 -- deliberately NOT touched by the 2026-08-28
    digitization, which replaced the POLARIZED curves).

    `cbt_unpolarized_emc_ratio` is the digitized 7Li unpolarized curve and
    is the better object; this one stays because `NuclearF2(emc_ratio=)`
    consumers (`money_delta_realistic.py`) and the `mode='constant'`
    legacy branches below are pinned to it."""
    xs = np.array([1e-4, 0.01, 0.06, 0.10, 0.20, 0.30, 0.45, 0.60,
                   0.70, 0.80, 0.88, 0.95])
    rs = np.array([0.96, 0.98, 1.00, 1.01, 1.00, 0.98, 0.95, 0.91,
                   0.88, 0.90, 1.00, 1.15])
    return np.interp(np.asarray(x, dtype=float), xs, rs)


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
# what transfers.  Both camps are therefore put on one common 7Li
# unpolarized baseline -- CBT's own digitized 7Li unpolarized curve --
# through
#
#     DR_model(x) = 1 - s_model * [1 - R_pol,model(x)],
#     s_model     = <1 - R_unpol,7Li> / <1 - R_unpol,model>   over the
#                   valence window 0.35 < x < 0.65.
#
# For CBT s = 1 identically (it IS the 7Li calculation), so the money plot
# draws CBT's published Eq.-23 curve untouched.  For TMT s = 0.397: the
# 7Li unpolarized depletion averages 0.0583 against nuclear matter's
# 0.1470 over that window.
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
# 0.028 < x < 0.30 (0.0016 at x = 0.09, 0.0006 at x = 0.14), so the 0.04
# separation the transferred pair shows there is the 7Li-vs-nuclear-matter
# strength rescaling and not a disagreement between the papers.  Inside the
# window, where the comparison does mean something, the transferred
# depletion tracks 7Li's own unpolarized 0.034 / 0.048 / 0.087 at
# x = 0.40 / 0.45 / 0.65 to within 0.005 against CBT's 0.077 / 0.082 /
# 0.094 -- "polarized ~ unpolarized" against "about twice" -- so DR
# separates by 0.040 at x = 0.36 and 0.011 at 0.65 as CBT's ratio of
# effects falls to 1.  `tmt_published_emc_ratio` returns the untransferred
# curve so that both separations can be quoted side by side, `money_polemc.py`
# prints them side by side, and Report 0 section 5.3 says so in words.
POLEMC_VALENCE_WINDOW = (0.35, 0.65)


def _valence_scale(table):
    """<1 - R_unpol,7Li> / <1 - R_unpol,table> over the valence window."""
    lo, hi = POLEMC_VALENCE_WINDOW
    g = np.linspace(lo, hi, 301)
    d_ref = 1.0 - _interp(CBT_TABLE, "R_unpol", g)
    d_mod = 1.0 - _interp(table, "R_unpol", g)
    return float(d_ref.mean() / d_mod.mean())


# Evaluated once at import (the tables are small and the lookup cached).
CBT_VALENCE_SCALE = 1.0          # CBT computed 7Li itself -- exact identity
TMT_VALENCE_SCALE = _valence_scale(TMT_TABLE)                    # = 0.397009


def cbt_unpolarized_emc_ratio(x):
    """Digitized 7Li UNPOLARIZED EMC ratio (CBT Fig. 6, blue dashed).

    The common baseline of the two-camp comparison.  Covers
    x = 0.028-0.871; constant-extrapolated outside, so the value below
    x = 0.028 is frozen at 1.006 and the Fermi rise above 0.871 is not
    followed."""
    return _interp(CBT_TABLE, "R_unpol", x)


def cbt_polarized_emc_ratio(x, mode="digitized", eq=23):
    """Cloet-Bentz-Thomas polarized EMC ratio for 7Li.

    mode='digitized' (default) returns the published curve of PLB 642:210
    Fig. 6 -- `eq=23` the red dotted R^{3/2 3/2}_{As} of their Eq. (23),
    which is what plans/01 defines the programme's DR_A to be, and
    `eq=26` the red solid R^{(3/2 1)}_{As} of their Eq. (26).  No transfer
    is applied because CBT computed 7Li itself (s = 1 above).

    mode='constant' is the pre-2026-08-28 stand-in, 1 - 2(1 - R_EMC) on
    the 12-point `unpolarized_emc_ratio` table, kept so the figures
    published before the digitization can be reproduced.

    DIGITIZED RANGE x = 0.028-0.871.  Below 0.028 the value is frozen at
    the endpoint 0.919 and above 0.871 at 1.007; `money_polemc.py` draws
    from x = 0.005, so the leftmost decade of that figure is a flat
    continuation, not a prediction."""
    if mode == "constant":
        return 1.0 - 2.0 * (1.0 - unpolarized_emc_ratio(x))
    if mode != "digitized":
        raise ValueError("mode must be 'digitized' or 'constant'")
    if eq not in (23, 26):
        raise ValueError("eq must be 23 (R^{3/2 3/2}) or 26 (R^{(3/2 1)})")
    return _interp(CBT_TABLE, "R_pol_eq%d" % eq, x)


def tmt_polarized_emc_ratio(x, mode="digitized"):
    """Tronchin-Matevosyan-Thomas polarized EMC ratio, transferred to 7Li.

    mode='digitized' (default) rescales the published nuclear-matter
    polarized depletion of PLB 783:247 Fig. 4 to 7Li strength with the
    valence factor s = 0.397 defined above, so that it is comparable with
    `cbt_polarized_emc_ratio` bin by bin.  mode='constant' is the
    pre-2026-08-28 stand-in, the unpolarized 12-point table itself.

    DIGITIZED RANGE x = 0.0015-0.739 (the published curve stops there);
    above 0.739 the value is frozen at the endpoint, which also freezes
    the nuclear-matter Fermi rise at a different x than 7Li's, so the
    x > 0.74 end of the money plot is not a like-for-like comparison."""
    if mode == "constant":
        return unpolarized_emc_ratio(x)
    if mode != "digitized":
        raise ValueError("mode must be 'digitized' or 'constant'")
    return 1.0 - TMT_VALENCE_SCALE * (1.0 - _interp(TMT_TABLE, "R_pol", x))


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
    between the papers, and inside it the transferred depletion tracks
    7Li's own unpolarized 0.034 / 0.048 / 0.087 at x = 0.40 / 0.45 / 0.65
    to within 0.005 against CBT's 0.077 / 0.082 / 0.094, so DR separates by
    0.040 at x = 0.36 and 0.011 at 0.65.  `money_polemc.py` prints both separations for this
    reason, and Report 0 section 5.3 says it in words.

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


def _toy_b1_shape(x, q2, f1):
    """The pre-2026-08-28 'HERMES-like' scenario shape times F1."""
    x = np.asarray(x, dtype=float)
    shape = 0.01 * np.power(np.maximum(x, 1e-6), -0.2) * (1.0 - x / 0.20)
    return shape * np.exp(-3.0 * x) * f1


def toy_b1(x, q2, f1, mode="digitized"):
    """b1 of the deuteron, per nucleon: Miller's pion + hidden-colour total.

    mode='digitized' (default) is PRC 89:045203 Fig. 5, the curve that
    reproduces the HERMES measurement with a hidden-colour probability of
    0.15% -- the 'HERMES-like' camp.  `q2` and `f1` are ignored (the Q2
    handle is `b1_miller_q2set.csv`, not wired in).  mode='toy' is the old
    analytic shape times F1.

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
    return B1_PER_DEUTERON_TO_PER_NUCLEON * _interp_tapered(MILLER_TABLE, "b1", x)


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
