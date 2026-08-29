"""Unpolarized structure functions: toy parameterizations + nuclear builder.

The TOY F2 below is a two-component (valence + sea) form, anchored by eye
to world ep data: F2p(1e-3, 10) ~ 1.2, F2p(0.1, 10) ~ 0.36, F2p(0.4, 10) ~ 0.12.
It is adequate for phase-space maps and ~factor-1.5 rate estimates ONLY.

Step 1.3 of the Phase-1 plan replaces this with LHAPDF (CT18 + EPPS21 or
nNNPDF) evaluated inside eic-shell, or the `parton` pip package locally.
Keep the F2Source interface so the swap is one line in the scripts.

R = sigma_L/sigma_T lives here too, in two flavours that are NOT
interchangeable: `r_sigma_lt` (the toy placeholder that is the default of
every published number) and `r1998` (the published SLAC/E143 world fit,
K. Abe et al., PLB 452 (1999) 194, arXiv:hep-ex/9808028).  Consumers take
an optional `r_func` hook defaulting to None = `r_sigma_lt`, so switching
R is a per-call-site decision and never a silent one.
"""

import numpy as np

ALPHA_EM = 1.0 / 137.036


def _safe_xfx(pdf, pid, x, q2):
    """Coerce xfxQ2 return value to a plain Python float.

    Guards against the NumPy ≥ 2.0 breakage where float() on a 0-d ndarray
    raises TypeError instead of extracting the scalar.
    """
    v = pdf.xfxQ2(pid, x, q2)
    try:
        return float(v)
    except TypeError:
        # NumPy ≥ 2 0-d array path
        return float(np.asarray(v).item())


GEV2_TO_PB = 0.3894e9  # (hbar c)^2 in GeV^2 * pb


class ToyF2:
    """Crude but smooth F2p / F2n with mild log Q2 evolution (TOY)."""

    def f2p(self, x, q2):
        x = np.asarray(x, dtype=float)
        q2 = np.asarray(q2, dtype=float)
        lq = np.log(np.maximum(q2, 1.1) / 0.04)  # ln(Q2/Lambda^2), Lambda=0.2
        # sea: low-x rise steepens slowly with Q2 (HERA-like lambda ~ 0.05 lnQ2)
        lam = 0.045 * lq
        sea = 0.20 * np.power(x, -lam) * np.power(1.0 - x, 7.0)
        val = 1.05 * np.power(x, 0.55) * np.power(1.0 - x, 3.0)
        return sea + val

    def f2n_over_f2p(self, x):
        # simple fit to NMC-like d/p ratio behavior
        return np.clip(1.0 - 0.75 * np.asarray(x, dtype=float), 0.25, 1.0)

    def f2n(self, x, q2):
        return self.f2p(x, q2) * self.f2n_over_f2p(x)


def r_sigma_lt(x, q2):
    """R = sigma_L/sigma_T, simplified R1990-like magnitude (TOY).

    This is the DEFAULT R of the whole program and every published
    polligen number (money plots 5/6/7 and 5R/7R/6R) was produced with
    it; it is a placeholder of the right magnitude (0.14-0.18 over the
    acceptance), not a fit.  `r1998` below is the published fit; switch
    to it per call site with the `r_func` hooks, never by editing here.
    """
    q2 = np.asarray(q2, dtype=float)
    return 0.18 / (1.0 + q2 / 50.0)


# --- R1998: the SLAC/E143 world fit to R = sigma_L/sigma_T ------------------
#
# K. Abe et al. (E143 Collaboration), "Measurements of R = sigma_L/sigma_T
# for 0.03 < x < 0.1 and Fit to World Data", Phys. Lett. B 452 (1999) 194
# (SLAC-PUB-7927, arXiv:hep-ex/9808028; refs/hep-ex_9808028.pdf).
#
# Table II of that paper, the coefficients of the three six-parameter forms
# of its Eq. (2).  Transcribed 2026-08-26 from the arXiv PDF, together with
# the chi2/df the paper quotes for 231 degrees of freedom.
R1998_PARAMS = {
    "a": (0.0485, 0.5470, 2.0621, -0.3804, 0.5090, -0.0285),   # chi2/df 0.9
    "b": (0.0481, 0.6114, -0.3509, -0.4611, 0.7172, -0.0317),  # chi2/df 0.9
    "c": (0.0577, 0.4644, 1.8288, 12.3708, -43.1043, 41.7415),  # chi2/df 1.0
}

# Support of the fit: "There were 237 points with a kinematic range of
# 0.005 <= x <= 0.86 and 0.5 <= Q2 <= 130" (paper, paragraph above Eq. 2).
R1998_X_RANGE = (0.005, 0.86)
R1998_Q2_RANGE = (0.5, 130.0)

# With clip=False the Q2 floor is still enforced here: the leading term of
# every form is 1/log(Q2/0.04), whose pole sits at Q2 = Lambda^2 = 0.04
# GeV2, so the parameterization changes sign just below the fit's own
# lower edge and means nothing there.
_R1998_Q2_POLE_FLOOR = 0.05


def _r1998_theta(x, q2):
    """Theta(x, Q2) of Abe et al. Eq. (3).

    Theta = 1 + 12 [Q2/(Q2+1)] [0.125^2/(0.125^2 + x^2)], which is
    6.2-12.9 over x <= 0.05, 1 <= Q2 <= 130 (4.5 at the Q2 = 0.5 support
    edge) -- measured, the code review's "6.5-9.5" is a little low and a
    little narrow.  It multiplies the log term of Eq. (2) and
    NOTHING ELSE -- multiplying it into the other terms as well saturates
    the sum and drives R to 1 over the whole low-x region (code review
    2026-08-25 item S1, the defect this function replaces).
    """
    return 1.0 + 12.0 * (q2 / (q2 + 1.0)) * (0.125 ** 2 / (0.125 ** 2 + x * x))


def _r1998_prepare(x, q2, clip):
    x = np.asarray(x, dtype=float)
    q2 = np.asarray(q2, dtype=float)
    if clip:
        x = np.clip(x, *R1998_X_RANGE)
        q2 = np.clip(q2, *R1998_Q2_RANGE)
    else:
        q2 = np.maximum(q2, _R1998_Q2_POLE_FLOOR)
        x = np.maximum(x, 1e-12)  # x^a6 with a6 < 0 diverges at x = 0
    return x, q2


def r1998_forms(x, q2, clip=True):
    """(R_a, R_b, R_c) of Abe et al. Eq. (2); R1998 is their average.

        R_a = a1 Theta/log(Q2/0.04) + a2 (Q2^4 + a3^4)^(-1/4)
                                         [1 + a4 x + a5 x^2] x^a6
        R_b = b1 Theta/log(Q2/0.04) + [b2/Q2 + b3/(Q2^2 + 0.3^2)]
                                         [1 + b4 x + b5 x^2] x^b6
        R_c = c1 Theta/log(Q2/0.04) + c2 [(Q2 - Q2thr)^2 + c3^2]^(-1/2),
              Q2thr = c4 x + c5 x^2 + c6 x^3

    (the paper writes Q^4 for (Q2)^2 and Q^8 for (Q2)^4; log is natural,
    0.04 GeV2 = Lambda^2).  The spread of the three IS the paper's
    functional-form systematic: "A systematic error associated with the
    functional form can be assigned from the spread of the three fits."
    """
    x, q2 = _r1998_prepare(x, q2, clip)
    rlog = _r1998_theta(x, q2) / np.log(q2 / 0.04)

    a1, a2, a3, a4, a5, a6 = R1998_PARAMS["a"]
    r_a = a1 * rlog + (a2 * np.power(q2 ** 4 + a3 ** 4, -0.25)
                       * (1.0 + a4 * x + a5 * x * x) * np.power(x, a6))

    b1, b2, b3, b4, b5, b6 = R1998_PARAMS["b"]
    r_b = b1 * rlog + ((b2 / q2 + b3 / (q2 * q2 + 0.3 ** 2))
                       * (1.0 + b4 * x + b5 * x * x) * np.power(x, b6))

    c1, c2, c3, c4, c5, c6 = R1998_PARAMS["c"]
    q2thr = c4 * x + c5 * x * x + c6 * x ** 3
    r_c = c1 * rlog + c2 * np.power((q2 - q2thr) ** 2 + c3 * c3, -0.5)

    return r_a, r_b, r_c


def r1998(x, q2, form="average", clip=True):
    """R = sigma_L/sigma_T from the SLAC/E143 world fit "R1998".

    K. Abe et al. (E143), Phys. Lett. B 452 (1999) 194
    (arXiv:hep-ex/9808028), Eqs. (2)-(3) with Table II; R1998 is defined
    there as the AVERAGE of the three forms R_a, R_b, R_c ("As in the
    case of R1990, we define R1998 to be the average of the three fits").

    Parameters
    ----------
    form : {"average", "a", "b", "c"}
        "average" is R1998 itself; the individual forms are exposed so
        that the fit's own functional-form systematic (their spread,
        `r1998_spread`) is visible to the caller rather than hidden.
    clip : bool
        The fit's support is 0.005 <= x <= 0.86, 0.5 <= Q2 <= 130 GeV2.
        clip=True (default) FREEZES the fit at the nearest boundary
        outside that box, so no extrapolation beyond the world data is
        ever implied.  clip=False evaluates the analytic form as written:
        it is smooth below x = 0.005 (Theta has saturated there -- the
        raw form is only 6% above the frozen one at x = 1e-4, Q2 = 2
        GeV2, measured) and falls slowly above Q2 = 130, but Q2 is still
        floored just above the log pole at Lambda^2 = 0.04 GeV2.

    Not clipped in R: over the whole support the three forms stay in
    0.012 <= R <= 0.45 (R_b is the largest, 0.4462 at x = 0.005, Q2 =
    0.97; R1998 itself, their average, in 0.012-0.403) -- measured on a
    60 x 60 log grid and confirmed on an 800 x 800 one.  So a clip in R
    would only ever mask a coding error, which is exactly how the defect
    this replaces stayed invisible (it was clipped to R = 1.000).
    """
    r_a, r_b, r_c = r1998_forms(x, q2, clip=clip)
    if form == "average":
        return (r_a + r_b + r_c) / 3.0
    try:
        return {"a": r_a, "b": r_b, "c": r_c}[form]
    except KeyError:
        raise ValueError("form must be 'average', 'a', 'b' or 'c'") from None


def r1998_spread(x, q2, clip=True):
    """max - min of (R_a, R_b, R_c): the fit's functional-form systematic.

    Abe et al. (arXiv:hep-ex/9808028), sentence below the delta-R formula.
    """
    forms = np.stack(np.broadcast_arrays(*r1998_forms(x, q2, clip=clip)))
    return forms.max(axis=0) - forms.min(axis=0)


def r1998_fit_error(x, q2):
    """delta-R of the R1998 fit, Abe et al. (arXiv:hep-ex/9808028).

    The unnumbered equation below their Eq. (3):
    delta-R = 0.0078 - 0.013 x + (0.070 - 0.39 x + 0.70 x^2)/(1.7 + Q2).
    This is the FITTING error only; the functional-form systematic is
    `r1998_spread` and long-range correlated errors enlarge both.
    """
    x = np.asarray(x, dtype=float)
    q2 = np.asarray(q2, dtype=float)
    return (0.0078 - 0.013 * x
            + (0.070 - 0.39 * x + 0.70 * x * x) / (1.7 + q2))


class PartonF2:
    """LO F2 from an LHAPDF grid via the pure-python `parton` package.

    F2 = sum_q e_q^2 (x q + x qbar); neutron by isospin u<->d swap.
    Install grids once:  yes | python3 -m parton install CT18NLO
    Drop-in replacement for ToyF2 (same interface).

    FLAVOUR SCHEME (documented 2026-08-28).  `_E2` runs over FIVE flavours,
    d u s c b, while `polarized.PartonG1._E2` runs over three, d u s.  That
    asymmetry is deliberate and physical, not an oversight: the unpolarized
    F2 a measurement sees contains charm and bottom, and NNPDFpol1.1 -- the
    polarized set behind PartonG1 -- assumes Delta c = Delta b = 0, so a
    three-flavour g1 IS the polarized structure function that set predicts.
    Every g1/F1 ratio built from the two is therefore the physical ratio,
    and it is smaller than a light-flavour-only ratio by the heavy-quark
    share of F2.  Measured on the (x, Q2) bins money_polemc.py combines
    (CT18NLO, 7Li, >= 100 events at 10 fb^-1/u): charm + bottom carry 7.8%
    of F2A event-weighted over all of them, 10% at x = 1e-3 to 1e-2 and up
    to 25% at the highest Q2 there, but only 0.65% at x = 0.3-0.5 and 0.23%
    at 0.5-0.7.

    WHETHER THAT COSTS THE HEADLINE ANYTHING depends on which x the
    polarized-EMC discrimination is read at, and after the 2026-08-28
    digitization those are two different windows.  Where the two published
    calculations genuinely differ is the valence window 0.35 < x < 0.65
    (`polarized.POLEMC_VALENCE_WINDOW`), and there the share is 0.31% at
    x = 0.35 falling to 0.03% at 0.65 at Q2 = 4 GeV2 (0.51% to 0.05% at
    Q2 = 10) -- one to two orders of magnitude below the 0.03-0.04
    separation in DR, so the mixed scheme costs that headline nothing.
    The same projection also plots low-x bins, where the share is 2.60% at
    x = 0.09 and 1.43% at 0.14 at Q2 = 4 (4.26% and 2.39% at Q2 = 10) --
    comparable to the 4.7% relative separation drawn there, which is a
    second reason (besides the valence-window transfer of
    `polarized.tmt_published_emc_ratio`) not to read the low-x bins as the
    discriminating ones.  `tests/test_grids.py` pins the two flavour sets
    against each other and both ends of this statement.
    """

    _E2 = {1: 1 / 9, 2: 4 / 9, 3: 1 / 9, 4: 4 / 9, 5: 1 / 9}

    def __init__(self, setname="CT18NLO", member=0):
        from parton import mkPDF  # lazy: optional dependency
        self._pdf = mkPDF(setname, member)
        self._f2p_vec = np.vectorize(self._f2p_scalar)

    def _f2p_scalar(self, x, q2):
        if not (0.0 < x < 1.0):
            return 0.0
        tot = 0.0
        for pid, e2 in self._E2.items():
            tot += e2 * (_safe_xfx(self._pdf, pid, x, q2)
                         + _safe_xfx(self._pdf, -pid, x, q2))
        return max(tot, 0.0)

    def f2p(self, x, q2):
        return self._f2p_vec(np.asarray(x, dtype=float),
                             np.asarray(q2, dtype=float))

    def f2n(self, x, q2):
        # isospin: swap the u and d charges
        x = np.asarray(x, dtype=float)
        q2 = np.asarray(q2, dtype=float)

        def scalar(xx, qq):
            if not (0.0 < xx < 1.0):
                return 0.0
            e2n = {1: 4 / 9, 2: 1 / 9, 3: 1 / 9, 4: 4 / 9, 5: 1 / 9}
            return max(sum(e2 * (_safe_xfx(self._pdf, p, xx, qq)
                                 + _safe_xfx(self._pdf, -p, xx, qq))
                           for p, e2 in e2n.items()), 0.0)

        return np.vectorize(scalar)(x, q2)

    def f2n_over_f2p(self, x, q2=10.0):
        f2p = self.f2p(x, q2)
        return np.where(f2p > 0, self.f2n(x, q2) / np.maximum(f2p, 1e-30),
                        1.0)


class NuclearF2:
    """Whole-nucleus F2A = Z*F2p + N*F2n, with an optional EMC-ratio hook.

    Per-nucleon F2 is F2A/A. The emc_ratio callable r(x) multiplies the
    isoscalar combination; default None means r = 1 (no medium effect),
    fine for rates (the EMC effect is a <20% shape effect here).

    `r_func(x, q2) -> R = sigma_L/sigma_T` selects the R that converts F2
    to F1 in `f1a`; None (the default) keeps the module-level
    `r_sigma_lt`, so every existing call site is bit-for-bit unchanged.
    Pass `r1998` for the published fit.
    """

    def __init__(self, ion, base=None, emc_ratio=None, r_func=None):
        self.ion = ion
        self.base = base or ToyF2()
        self.emc_ratio = emc_ratio
        self.r_func = r_func

    def f2a(self, x, q2):
        f2 = self.ion.Z * self.base.f2p(x, q2) + self.ion.N * self.base.f2n(x, q2)
        if self.emc_ratio is not None:
            f2 = f2 * self.emc_ratio(x)
        return f2

    def f1a(self, x, q2):
        """F1 via the Callan-Gross relation modified by R (massless).

        The default R is looked up as the MODULE global `r_sigma_lt` at
        call time, not captured in __init__, because the dated money
        scripts rebind `structure.r_sigma_lt` inside their `r_override`
        context manager (fastsim/scripts/money_delta_*.py) and those
        scripts are frozen reproductions of dated notes.
        """
        r = (r_sigma_lt if self.r_func is None else self.r_func)(x, q2)
        return self.f2a(x, q2) / (2.0 * np.asarray(x, dtype=float) * (1.0 + r))


def dsigma_dx_dq2(x, q2, s, f2, fl=None, r_func=None):
    """NC DIS double-differential cross section [pb/GeV^2].

    d2sigma/dxdQ2 = 4 pi alpha^2 / (x Q^4) [ (1 - y + y^2/2) F2 - y^2/2 FL ]
    FL defaults to F2 * R/(1+R), with R = `r_func(x, q2)`; r_func=None
    keeps the module-level `r_sigma_lt` (looked up at call time so the
    dated scripts' `r_override` monkey-patch still bites) and is
    therefore bit-for-bit today's result.  An explicit `fl` overrides R
    entirely, as before.
    """
    x = np.asarray(x, dtype=float)
    q2 = np.asarray(q2, dtype=float)
    y = q2 / (s * x)
    if fl is None:
        r = (r_sigma_lt if r_func is None else r_func)(x, q2)
        fl = f2 * r / (1.0 + r)
    bracket = (1.0 - y + 0.5 * y * y) * f2 - 0.5 * y * y * fl
    xsec = 4.0 * np.pi * ALPHA_EM**2 / (x * q2**2) * bracket
    return np.maximum(xsec, 0.0) * GEV2_TO_PB
