"""Polarized / tensor / gluonometric inputs (TOY models + hooks).

Every quantity here is a *scenario input*, not a prediction. The toy shapes
have roughly the right magnitude and x-dependence so that FOM machinery can
be exercised end-to-end; Phase-1 steps 1.3-1.5 replace them with:
  - g1p, g1n: JAM / DSSV grids (LHAPDF)
  - polarized EMC ratio: digitized Cloet-Bentz-Thomas curves (PLB 642:210)
  - b1: Cosyn-Dong-Kumano-Song deuteron convolution + 6Li cluster model
  - Delta(x,Q2): lattice-moment-normalized model (Detmold-Shanahan) or
    scenario values Delta/F1 in {1e-3 ... 1e-2}
"""

import numpy as np

from .structure import ToyF2, r_sigma_lt, _safe_xfx


class ToyG1:
    """g1 = A1 * F1 with toy A1(x) shapes (TOY; replace with JAM/DSSV)."""

    def __init__(self, base=None):
        self.base = base or ToyF2()

    def a1p(self, x):
        # x^0.7 tracks moderate/high-x world data; saturates below 1
        return np.clip(np.power(np.asarray(x, dtype=float), 0.7), 0.0, 1.0)

    def a1n(self, x):
        x = np.asarray(x, dtype=float)
        # small negative at low/mid x, positive rise at high x
        return -0.07 * np.power(1.0 - x, 2.0) + 0.8 * np.power(x, 2.2)

    def _f1(self, f2, x, q2):
        return f2 / (2.0 * x * (1.0 + r_sigma_lt(x, q2)))

    def g1p(self, x, q2):
        return self.a1p(x) * self._f1(self.base.f2p(x, q2), x, q2)

    def g1n(self, x, q2):
        return self.a1n(x) * self._f1(self.base.f2n(x, q2), x, q2)

    def g1_nucleus(self, ion, x, q2, medium_ratio=None):
        """g1A = P_p g1p + P_n g1n (naive), times optional medium ratio.

        `medium_ratio(x)` is the polarized-EMC modification DR(x); the
        observable of interest is precisely deviations of DR from 1.
        """
        g1 = ion.eff_pol_p * self.g1p(x, q2) + ion.eff_pol_n * self.g1n(x, q2)
        if medium_ratio is not None:
            g1 = g1 * medium_ratio(x)
        return g1


class PartonG1(ToyG1):
    """LO g1 from a polarized LHAPDF grid via `parton` (default
    NNPDFpol11_100): g1 = (1/2) sum_q e_q^2 [Dq + Dqbar]; neutron by
    isospin. Unpolarized denominator (F1) still comes from `base`.
    Install:  yes | python3 -m parton install NNPDFpol11_100
    """

    _E2 = {1: 1 / 9, 2: 4 / 9, 3: 1 / 9}

    def __init__(self, base=None, setname="NNPDFpol11_100", member=0):
        super().__init__(base=base)
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
    shape (SCENARIO; replace with EPPS21 / data fits in step 1.2)."""
    xs = np.array([1e-4, 0.01, 0.06, 0.10, 0.20, 0.30, 0.45, 0.60,
                   0.70, 0.80, 0.88, 0.95])
    rs = np.array([0.96, 0.98, 1.00, 1.01, 1.00, 0.98, 0.95, 0.91,
                   0.88, 0.90, 1.00, 1.15])
    return np.interp(np.asarray(x, dtype=float), xs, rs)


def cbt_polarized_emc_ratio(x):
    """Cloet-Bentz-Thomas-like scenario: polarized EMC effect ~2x the
    unpolarized one (verbatim claim of PLB 642:210; curve shape here is
    qualitative -- digitize the published figure before publication)."""
    return 1.0 - 2.0 * (1.0 - unpolarized_emc_ratio(x))


def tmt_polarized_emc_ratio(x):
    """Tronchin-Matevosyan-Thomas-like scenario (PLB 783:247): polarized
    EMC effect comparable to the unpolarized one."""
    return unpolarized_emc_ratio(x)


def toy_polarized_emc_ratio(x):
    """Backward-compatible alias for the CBT-like scenario."""
    return cbt_polarized_emc_ratio(x)


def toy_b1(x, q2, f1):
    """SCENARIO b1 'HERMES-like': large at low x (b1 ~ 0.1 at x ~ 0.01,
    cf. HERMES PRL 95:242001), zero crossing near x ~ 0.2, small negative
    at high x (Miller pion + hidden-color shape). Returns b1 (absolute).
    """
    x = np.asarray(x, dtype=float)
    shape = 0.01 * np.power(np.maximum(x, 1e-6), -0.2) * (1.0 - x / 0.20)
    return shape * np.exp(-3.0 * x) * f1


def b1_convolution(x, q2, f1):
    """SCENARIO b1 'standard convolution': an order of magnitude below the
    HERMES-like curve at x > 0.2 (Cosyn-Dong-Kumano-Sargsian PRD 95:074036
    find |b1| < 1e-3 there). Same zero-crossing structure."""
    return 0.1 * toy_b1(x, q2, f1)


def b1_li6_from_deuteron(b1_d, p_d_cluster=0.87):
    """Embedded-deuteron scaling for whole-nucleus 6Li b1:
    b1(6Li) ~ P_d * b1(d) for the d-cluster contribution (our inference --
    no published 6Li b1 exists; see plans/04 item 9). Per-nucleon
    normalization adds the 2/6 dilution where needed."""
    return p_d_cluster * b1_d


def toy_delta_gluon(x, q2, f1, scale=1e-3):
    """TOY double-helicity-flip Delta(x,Q2) = scale * F1 * (1-x)^4 * x^0.3.

    Delta is alpha_s-suppressed and unknown; `scale` sets Delta/F1 at its
    maximum. Use scale in {1e-3, 3e-3, 1e-2} as discovery scenarios, to be
    replaced by lattice-moment-anchored models.
    """
    x = np.asarray(x, dtype=float)
    return scale * f1 * np.power(x, 0.3) * np.power(1.0 - x, 4.0)
