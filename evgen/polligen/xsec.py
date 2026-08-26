"""Inclusive doubly polarized master cross section (the polligen "wheel").

Per-nucleon differential cross section for a polarized electron
(helicity lam_e = +-1, magnitude P_e) on a spin-J ion in the definite
projection state m along axis n(theta_S, phi_S):

  dsigma/(dx dQ2 dphi) = sigma_unpol(x,Q2)/(2 pi) * W(phi')
  W = 1 + w_avg + a_1 cos(phi') + a_2 cos(2 phi'),      phi' = phi - phi_S

with sigma_unpol = `polli_fastsim.structure.dsigma_dx_dq2` (per-nucleon F2A/A,
FL from R) and, writing D_phi = F1 + (1-y)/(x y^2) F2 (the phi-averaged
density of `asymmetries.py`), K = b1 + (1-y)/(x y^2) b2:

  w_avg = t_geo(m, theta_S) * K / D_phi                       [tensor, J>=1]
        + lam_e P_e (m/J) cos(theta_S) * A_par(x,y)           [vector-L]
  a_1   = lam_e P_e (m/J) sin(theta_S) * A_perp(x,y)          [vector-T, gT]
  a_2   = -(1-y)/y^2 * c_eff(m) sin^2(theta_S) * Delta / D_phi  [gluonometry]

Spin-1 (Hoodbhoy-Jaffe-Manohar NPB 312:571, exactly the master formula of
`polli_fastsim.asymmetries`): c_m = 3 m^2 - 2 -> (1,-2,1);
t_geo = (2/3) a_m with a_m = (1/4) c_m (3 cos^2 theta_S - 1); c_eff = c_m.
The overall sign of the b1/b2 (tensor RATE) sector is the single constant
`asymmetries.TENSOR_LL_SIGN`; it is opposite to Cosyn et al. Eq. (27) and
is an open author decision (plans/08 D1).  The Delta sector does not
depend on it.

Spin-3/2 (7Li): rank-0/1 exact via (F1, F2, g1, g2); rank-2 is a SCENARIO
slot (plans/04 #14) with normalized T_m = (3 m^2 - 15/4)/3 -> (1,-1,-1,1):
t_geo = T_m P2(cos theta_S), c_eff = T_m; both default to zero SFs.
Rank-3 dropped (plans/05 truncation).

Vector-sector y-factors follow the fastsim approximation set:
  A_par  = D(y) g1/F1                       (`asymmetries.a_parallel`)
  A_perp = D(y) sqrt(2 eps/(1+eps)) * gamma ((y/2) g1 + g2)/F1,
           eps = (1-y)/(1-y+y^2/2),  gamma = 2 M x / Q
A_par = D*A1 drops the gamma^2-suppressed eta*A2 term (consistent with the
fastsim); A_perp = d*(A2 - xi*A1) keeps BOTH O(gamma) pieces -- in the
massless-kinematics limit gamma - xi = gamma*y/2 exactly, giving the
textbook transverse combination (y/2) g1 + g2 (E143 conventions,
PRD 58:112003).  g2 defaults to Wandzura-Wilczek from the g1 backend; the
exact Cosyn-Weiss y-factors can replace these two functions behind the
same interface (plans/05 step 5.B+).

Consistency gates (tested bin-wise in evgen/tests/test_xsec_identity.py):
with vector-only polarization reproduce `asymmetries.a_parallel`; with
tensor-only reproduce `asymmetries.azz` and `asymmetries.a_cos2phi`.
"""

from dataclasses import dataclass

import numpy as np

from polli_fastsim.asymmetries import TENSOR_LL_SIGN, depolarization_d
from polli_fastsim.structure import NuclearF2, ToyF2, dsigma_dx_dq2
from polli_fastsim.polarized import ToyG1

M_NUCLEON = 0.9383  # GeV (per-nucleon target mass for gamma = 2Mx/Q)

# NumPy compat: np.trapz removed in NumPy >= 2.4, np.trapezoid absent < 2.0
_trapezoid = getattr(np, "trapezoid", None) or np.trapz


def g2_ww(g1_func, x, q2, npts=96):
    """Wandzura-Wilczek g2(x) = -g1(x) + int_x^1 du g1(u)/u at fixed Q2.

    Substitution u = x^(1-t) maps the integral to -ln(x) int_0^1 g1(x^(1-t)) dt,
    evaluated by trapezoid; accepts any mutually broadcastable (x, q2).
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


@dataclass(frozen=True)
class EventSpinState:
    """Spin configuration of one bunch crossing category / event."""
    lam_e: int          # electron helicity sign (+1/-1); 0 = unpolarized
    pe: float           # electron polarization magnitude in [0, 1]
    j: float            # ion spin (1 or 3/2; 1/2 supported vector-only)
    m: float            # projection along the quantization axis
    theta_s: float = 0.0  # axis polar angle in the lab [rad]
    phi_s: float = 0.0    # axis azimuth in the lab [rad]


class InclusiveKernel:
    """Doubly polarized inclusive e+A cross section on polli_fastsim SFs.

    All structure functions are per-nucleon (F2A/A convention of fom.py).
    b1/b2/Delta enter as callables (x, q2, f1) -> SF like polarized.py's
    scenario curves; b2 defaults to 2*x*b1 (asymmetries.azz convention).
    Spin-3/2 rank-2 slots (b1_32/delta_32) default to None -> zero.
    """

    def __init__(self, ion, f2_source=None, g1_model=None,
                 b1_func=None, b2_func=None, delta_func=None,
                 b1_32_func=None, b2_32_func=None, delta_32_func=None,
                 g2_mode="ww", emc_ratio=None):
        self.ion = ion
        self.nf2 = NuclearF2(ion, base=f2_source or ToyF2(),
                             emc_ratio=emc_ratio)
        self.g1_model = g1_model or ToyG1(base=self.nf2.base)
        self.b1_func = b1_func
        self.b2_func = b2_func
        self.delta_func = delta_func
        self.b1_32_func = b1_32_func
        self.b2_32_func = b2_32_func
        self.delta_32_func = delta_32_func
        if g2_mode not in ("ww", "zero"):
            raise ValueError("g2_mode must be 'ww' or 'zero'")
        self.g2_mode = g2_mode

    # --- structure-function tables -------------------------------------

    def _g1a(self, x, q2):
        return self.g1_model.g1_nucleus(self.ion, x, q2) / self.ion.A

    def tables(self, x, q2, with_g2=False):
        """Per-nucleon SF tables on arrays (x, q2)."""
        x = np.asarray(x, dtype=float)
        q2 = np.asarray(q2, dtype=float)
        f2 = self.nf2.f2a(x, q2) / self.ion.A
        f1 = self.nf2.f1a(x, q2) / self.ion.A
        out = {"f1": f1, "f2": f2, "g1": self._g1a(x, q2)}
        zeros = np.zeros(np.broadcast(x, q2).shape)
        j = self.ion.spin
        if abs(j - 1.0) < 1e-9:
            b1 = self.b1_func(x, q2, f1) if self.b1_func else zeros
            b2 = self.b2_func(x, q2, f1) if self.b2_func else 2.0 * x * b1
            delta = self.delta_func(x, q2, f1) if self.delta_func else zeros
        elif abs(j - 1.5) < 1e-9:
            b1 = self.b1_32_func(x, q2, f1) if self.b1_32_func else zeros
            b2 = (self.b2_32_func(x, q2, f1) if self.b2_32_func
                  else 2.0 * x * b1)
            delta = (self.delta_32_func(x, q2, f1) if self.delta_32_func
                     else zeros)
        else:
            b1, b2, delta = zeros, zeros, zeros
        out.update(b1=b1, b2=b2, delta=delta)
        if with_g2:
            out["g2"] = (g2_ww(self._g1a, x, q2)
                         if self.g2_mode == "ww" else zeros)
        return out

    # --- kinematic factors ----------------------------------------------

    @staticmethod
    def _dphi(t, x, y):
        return t["f1"] + (1.0 - y) / (x * y * y) * t["f2"]

    @staticmethod
    def _tensor_kernel(t, x, y):
        return t["b1"] + (1.0 - y) / (x * y * y) * t["b2"]

    def a_parallel(self, t, x, q2, y):
        """A_par(x,y) = D(y) g1/F1 (identical to asymmetries.a_parallel)."""
        return depolarization_d(y, x, q2) * t["g1"] / np.maximum(t["f1"],
                                                                 1e-30)

    def a_perp(self, t, x, q2, y):
        """gamma-suppressed transverse-vector amplitude d(y)*(A2 - xi*A1).

        Both O(gamma) pieces are kept: in massless kinematics
        gamma - xi = gamma*y/2, so the amplitude reduces to
        d(y) * gamma * ((y/2) g1 + g2) / F1.
        """
        if "g2" not in t:
            raise KeyError("tables(..., with_g2=True) required for a_perp")
        eps = (1.0 - y) / (1.0 - y + 0.5 * y * y)
        d = depolarization_d(y, x, q2) * np.sqrt(
            np.maximum(2.0 * eps / (1.0 + eps), 0.0))
        gamma = 2.0 * M_NUCLEON * x / np.sqrt(q2)
        amp = gamma * (0.5 * y * t["g1"] + t["g2"]) / np.maximum(t["f1"],
                                                                 1e-30)
        return d * amp

    def _tensor_moments(self, m):
        """(t_geo prefactor c-like, c_eff) for a pure state m."""
        j = self.ion.spin
        if abs(j - 1.0) < 1e-9:
            c_m = 3.0 * m * m - 2.0
            return c_m, c_m
        if abs(j - 1.5) < 1e-9:
            t_m = (3.0 * m * m - j * (j + 1)) / 3.0
            return t_m, t_m
        return 0.0, 0.0

    # --- modulation amplitudes -------------------------------------------

    def amplitudes(self, t, x, q2, s, state, with_perp=False):
        """(w_avg, a_1, a_2) of W = 1 + w_avg + a_1 cos phi' + a_2 cos 2phi'.

        `t` are tables() on the same (x, q2) arrays; `state` an
        EventSpinState.  a_1 is only computed when with_perp (needs g2).
        """
        x = np.asarray(x, dtype=float)
        q2 = np.asarray(q2, dtype=float)
        y = q2 / (s * x)
        dphi = self._dphi(t, x, y)
        j = state.j
        ct, st = np.cos(state.theta_s), np.sin(state.theta_s)

        w_avg = np.zeros(np.broadcast(x, q2).shape)
        a_2 = np.zeros_like(w_avg)
        if j >= 1.0 - 1e-9:
            c_like, c_eff = self._tensor_moments(state.m)
            if abs(j - 1.0) < 1e-9:
                # HJM: (2/3) a_m K / D_phi, a_m = (1/4) c_m (3 cos^2 - 1)
                t_geo = (c_like / 6.0) * (3.0 * ct * ct - 1.0)
            else:
                t_geo = c_like * 0.5 * (3.0 * ct * ct - 1.0)
            t_geo = TENSOR_LL_SIGN * t_geo   # asymmetries.TENSOR_LL_SIGN
            kern = self._tensor_kernel(t, x, y)
            w_avg = w_avg + t_geo * kern / np.maximum(dphi, 1e-30)
            a_2 = (-(1.0 - y) / (y * y) * c_eff * st * st
                   * t["delta"] / np.maximum(dphi, 1e-30))

        helicity = state.lam_e * state.pe
        v = state.m / j if j > 0 else 0.0
        if helicity != 0.0 and v != 0.0:
            w_avg = w_avg + helicity * v * ct * self.a_parallel(t, x, q2, y)

        a_1 = np.zeros_like(w_avg)
        if with_perp and helicity != 0.0 and v != 0.0 and abs(st) > 1e-12:
            a_1 = helicity * v * st * self.a_perp(t, x, q2, y)
        return w_avg, a_1, a_2

    # --- differential cross sections --------------------------------------

    def dsigma_unpol(self, x, q2, s):
        """Unpolarized per-nucleon d2sigma/dxdQ2 [pb/GeV^2] (fastsim's)."""
        f2 = self.nf2.f2a(x, q2) / self.ion.A
        return dsigma_dx_dq2(x, q2, s, f2)

    def dsigma(self, x, q2, phi, s, state, with_perp=False):
        """Doubly polarized d3sigma/dxdQ2dphi [pb/GeV^2/rad]."""
        t = self.tables(x, q2, with_g2=with_perp)
        w_avg, a_1, a_2 = self.amplitudes(t, x, q2, s, state,
                                          with_perp=with_perp)
        phip = np.asarray(phi, dtype=float) - state.phi_s
        w = 1.0 + w_avg + a_1 * np.cos(phip) + a_2 * np.cos(2.0 * phip)
        return self.dsigma_unpol(x, q2, s) / (2.0 * np.pi) * np.maximum(w, 0.0)
