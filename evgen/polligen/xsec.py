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
slot (plans/04 #14) sharing the SAME geometry as spin 1 --
Q_NN = [3m^2 - J(J+1)]/3 -> (1,-1,-1,1) x 1 for J = 3/2, t_geo = Q_NN
P2(cos theta_S), c_eff = 3 Q_NN (Cosyn Eq. 9); both default to zero SFs.
Rank-3 dropped (plans/05 truncation).

Vector-sector y-factors:
  A_par  = D_gamma (A1 + eta A2)            (`asymmetries.a_parallel_exact`)
           A1 = (g1 - gamma^2 g2)/F1,  A2 = gamma (g1 + g2)/F1
  A_perp = D(y) sqrt(2 eps/(1+eps)) * gamma ((y/2) g1 + g2)/F1,
           eps = (1-y)/(1-y+y^2/2),  gamma = 2 M x / Q
A_perp = d*(A2 - xi*A1) keeps BOTH O(gamma) pieces -- in the
massless-kinematics limit gamma - xi = gamma*y/2 exactly, giving the
textbook transverse combination (y/2) g1 + g2 (E143 conventions,
PRD 58:112003).  g2 is Wandzura-Wilczek from the g1 backend; the exact
Cosyn-Weiss y-factors can replace these two functions behind the same
interface (plans/05 step 5.B+).

A_par carries its target-mass term by DEFAULT (target_mass=True, author
decision of 2026-08-29); target_mass=False is the massless limit
D(y) g1/F1 that every figure published before that date was made on, kept
reachable and pinned bit-for-bit.  The default was flipped because the
term is exact given g2 and costs one g2^WW table per grid: leaving it out
would have meant carrying its size as a systematic rather than computing
it.  That size is bounded by the W^2 >= 10 GeV^2 cut of `fom.Scenario`,
which caps gamma^2 <= M^2/(W2_min - M^2) = 0.0965 over the entire
accepted phase space (grid maxima 0.085 / 0.058 / 0.026 at the three 6Li
settings, 0.085 / 0.058 / 0.033 for 7Li), and at the small y of this
programme the whole correction collapses to a multiplicative
(1 + gamma^2) on A_par, independent of g2.  What is left over is the
twist-3 uncertainty on g2 itself, which `evgen/scripts/target_mass_bound.py`
measures by repeating the extraction with g2 = 0 and g2 = 1.5 g2^WW.
The tensor sector's own O(gamma^2) terms (Cosyn Eqs. 17d/17e) are a
separate open item (plans/08 D2) and target_mass does not touch them.

Consistency gates (tested bin-wise in evgen/tests/test_xsec_identity.py):
with vector-only polarization reproduce `asymmetries.a_parallel`; with
tensor-only reproduce `asymmetries.azz` and `asymmetries.a_cos2phi`.
"""

from dataclasses import dataclass

import numpy as np

from polli_fastsim import asymmetries as _asym
from polli_fastsim.asymmetries import (M_NUCLEON, TENSOR_LL_SIGN,
                                       a_parallel_exact,
                                       depolarization_d,
                                       depolarization_d_gamma,
                                       epsilon_gamma, eta_gamma,
                                       gamma_squared)
from polli_fastsim.structure import NuclearF2, ToyF2, dsigma_dx_dq2
from polli_fastsim.polarized import ToyG1, g2_ww

# The finite-gamma kinematics (eps, D_gamma, eta, gamma^2) and the
# Wandzura-Wilczek quadrature now live in `polli_fastsim` -- ONE
# implementation for both packages, so the generator kernel and the fast
# simulation cannot drift apart in the O(gamma^2) term -- and are
# re-exported here under the names this module has always used, so that
# `from polligen.xsec import gamma_squared, g2_ww, ...` keeps working.
depolarization_gamma = depolarization_d_gamma

__all__ = ["EventSpinState", "InclusiveKernel", "M_NUCLEON",
           "depolarization_gamma", "epsilon_gamma", "eta_gamma",
           "g2_ww", "gamma_squared"]


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

    `r_func(x, q2) -> R = sigma_L/sigma_T` is threaded to ALL FOUR places
    the kernel needs R -- F1 = F2/(2x(1+R)) in `NuclearF2.f1a`, F_L in
    `dsigma_dx_dq2`, D(y) in `depolarization_d`, and the g1/F1 of the
    default `ToyG1` -- so that one argument moves R consistently instead
    of the three-of-four monkey-patch of the dated fast-sim scripts
    (plans/08 C2).  None (the default) is `structure.r_sigma_lt`
    everywhere, i.e. bit-for-bit every published polligen number.  An
    EXPLICIT `g1_model` keeps whatever R it was built with: r_func only
    fills the default one.

    `g2_mode` chooses the g2 model ("ww", the default Wandzura-Wilczek
    table of the g1 backend, or "zero") and `g2_scale` multiplies it.
    Together they are the twist-3 handle on the finite-gamma A_par.  The
    residual systematic that replaced the target-mass bias is NOT taken by
    re-running an extraction through this knob: `target_mass_bound.
    g2_residual` builds one table and forms D_eff(rho)/D_eff(s rho) - 1
    from it, which is the same variation at a hundredth of the cost
    because g2^WW is linear in g1 and rho therefore scales with s
    exactly.  `g2_scale` is the kernel-level equivalent of that s -- what
    a sampler or a written-out event set has to vary through, having no
    rho of its own to rescale -- and it is pinned against the analytic
    form in `tests/test_target_mass.py`.  A kernel with target_mass=True
    and no g2 at all (g2_mode="zero") is legitimate and is exactly the
    g2_scale = 0 variation, so it is no longer refused.

    `target_mass` (DEFAULT True since 2026-08-29) selects the longitudinal
    vector kernel: True is the exact finite-gamma D_gamma (A1 + eta A2) of
    `asymmetries.a_parallel_exact`, which needs a g2 column and refuses a
    table built without one (the KeyError in `a_parallel`, not a
    restriction on g2_mode); False is the massless A_par = D g1/F1 of
    the fast simulation, bit-for-bit the path every figure published
    before that date was made on, kept reachable and pinned by
    `tests/test_target_mass.py`.  The default was flipped because the term
    is exact and free -- one g2^WW table per grid, nothing per call -- so
    leaving it off would have meant carrying its O(gamma^2) effect as a
    systematic on Delta-R instead of simply computing it.  What survives
    the flip is the twist-3 uncertainty on g2 itself, measured by
    `evgen/scripts/target_mass_bound.py`.
    """

    def __init__(self, ion, f2_source=None, g1_model=None,
                 b1_func=None, b2_func=None, delta_func=None,
                 b1_32_func=None, b2_32_func=None, delta_32_func=None,
                 g2_mode="ww", g2_scale=1.0, emc_ratio=None, r_func=None,
                 target_mass=True):
        self.ion = ion
        self.r_func = r_func
        self.nf2 = NuclearF2(ion, base=f2_source or ToyF2(),
                             emc_ratio=emc_ratio, r_func=r_func)
        self.g1_model = g1_model or ToyG1(base=self.nf2.base, r_func=r_func)
        self.b1_func = b1_func
        self.b2_func = b2_func
        self.delta_func = delta_func
        self.b1_32_func = b1_32_func
        self.b2_32_func = b2_32_func
        self.delta_32_func = delta_32_func
        if g2_mode not in ("ww", "zero"):
            raise ValueError("g2_mode must be 'ww' or 'zero'")
        self.g2_mode = g2_mode
        self.g2_scale = float(g2_scale)
        self.target_mass = bool(target_mass)

    # --- structure-function tables -------------------------------------

    def _g1a(self, x, q2):
        return self.g1_model.g1_nucleus(self.ion, x, q2) / self.ion.A

    def tables(self, x, q2, with_g2=False):
        """Per-nucleon SF tables on arrays (x, q2).

        g2 is included when `with_g2` (a_perp needs it) or whenever
        `target_mass` is on -- i.e. by default -- so that callers which
        never ask for the transverse sector (`sample.InclusiveSampler`,
        `tagged`) get a table the finite-gamma `a_parallel` can use.

        It is fetched from the g1 backend's own `g2_nucleus`, which caches
        the Wandzura-Wilczek quadrature per (x, Q2) grid; that cache is
        what makes the default-on target mass free on a PDF grid, where
        the 96-node quadrature would otherwise dominate a projection.  A
        backend without the method (any duck-typed g1 model) falls back to
        `g2_ww` on the per-nucleon g1 directly.
        """
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
        if with_g2 or self.target_mass:
            out["g2"] = self._g2a(x, q2) if self.g2_mode == "ww" else zeros
        return out

    def _g2a(self, x, q2):
        """Per-nucleon g2, `g2_scale` times Wandzura-Wilczek.

        The table comes from the backend's cached `g2_nucleus` when it has
        one; g2^WW is linear in g1, so dividing by A before or after the
        quadrature is the same number up to the last bit.  `g2_scale` is
        the twist-3 handle: 0 and 1.5 are the endpoints of the variation
        `evgen/scripts/target_mass_bound.py` quotes the residual
        systematic between -- that script forms the same variation
        analytically from rho, and this is its kernel-level equivalent --
        and 1 is the model."""
        cached = getattr(self.g1_model, "g2_nucleus", None)
        if cached is not None:
            out = cached(self.ion, x, q2) / self.ion.A
        else:
            out = g2_ww(self._g1a, x, q2)
        return out if self.g2_scale == 1.0 else self.g2_scale * out

    # --- kinematic factors ----------------------------------------------

    @staticmethod
    def _dphi(t, x, y):
        return t["f1"] + (1.0 - y) / (x * y * y) * t["f2"]

    @staticmethod
    def _tensor_kernel(t, x, y):
        return t["b1"] + (1.0 - y) / (x * y * y) * t["b2"]

    def a_parallel(self, t, x, q2, y):
        """Longitudinal double-spin asymmetry A_par(x, y).

        Default (`target_mass=True`): the exact E143 form
        (PRD 58:112003) computed by `asymmetries.a_parallel_exact`,

            A_par = D_gamma (A1 + eta A2),
            A1 = (g1 - gamma^2 g2)/F1,   A2 = gamma (g1 + g2)/F1,

        with gamma^2 = 4 M^2 x^2/Q^2 and the finite-gamma eps, D_gamma,
        eta of `polli_fastsim.asymmetries`.  With `target_mass=False` it
        is the massless limit A_par = D(y) g1/F1, bit-for-bit
        `asymmetries.a_parallel`, which is how every figure published
        before 2026-08-29 was made.

        R is resolved exactly as `depolarization_d` resolves it --
        self.r_func, else the `asymmetries` module global at call time,
        which the dated money scripts rebind.
        """
        if not self.target_mass:
            return (depolarization_d(y, x, q2, r_func=self.r_func)
                    * t["g1"] / np.maximum(t["f1"], 1e-30))
        if "g2" not in t:
            raise KeyError("target_mass=True needs g2 in tables(); build "
                           "them with this kernel, not a massless one")
        r_func = (_asym.r_sigma_lt if self.r_func is None else self.r_func)
        return a_parallel_exact(t["g1"], t["g2"], t["f1"], y, x, q2,
                                r_func=r_func)

    def a_perp(self, t, x, q2, y):
        """gamma-suppressed transverse-vector amplitude d(y)*(A2 - xi*A1).

        Both O(gamma) pieces are kept: in massless kinematics
        gamma - xi = gamma*y/2, so the amplitude reduces to
        d(y) * gamma * ((y/2) g1 + g2) / F1.
        """
        if "g2" not in t:
            raise KeyError("tables(..., with_g2=True) required for a_perp")
        eps = (1.0 - y) / (1.0 - y + 0.5 * y * y)
        d = depolarization_d(y, x, q2, r_func=self.r_func) * np.sqrt(
            np.maximum(2.0 * eps / (1.0 + eps), 0.0))
        gamma = 2.0 * M_NUCLEON * x / np.sqrt(q2)
        amp = gamma * (0.5 * y * t["g1"] + t["g2"]) / np.maximum(t["f1"],
                                                                 1e-30)
        return d * amp

    def _tensor_moments(self, m):
        """Rank-2 alignment of a pure state m, in ONE form for every spin.

        Cosyn et al. (arXiv:2410.12764) Eq. (9) writes the alignment
        tensor of any spin-J state as t_ij = (Q_NN/2)(3 n_i n_j - d_ij)
        with the single scalar

            Q_NN(m) = [3 m^2 - J(J+1)] / 3,

        so the longitudinal-longitudinal and transverse-transverse
        projections are T_LL = Q_NN P_2(cos Theta) and
        T_TT = (3/2) Q_NN sin^2 Theta.  Returning (Q_NN, 3 Q_NN)
        therefore gives ONE geometry for both spins:

            t_geo = Q_NN P_2(cos theta_S)     (the b1/b2 rate shift)
            c_eff = 3 Q_NN                    (the cos 2phi coefficient)

        For J = 1, Q_NN = c_m/3 with c_m = 3m^2 - 2, so this reproduces
        the Hoodbhoy-Jaffe-Manohar transcription exactly, digit for digit.
        For J = 3/2 it CHANGES the cos 2phi channel by a factor 3: the
        previous code returned (Q_NN, Q_NN), i.e. its rate and cos 2phi
        channels were inconsistent with each other by that factor.  The
        7Li rank-2 structure functions default to None, so nothing
        published moves; the overall rank-2 normalization convention for
        spin 3/2 is plans/04 #14 and this is the adopted one.
        """
        j = self.ion.spin
        if j < 1.0 - 1e-9:
            return 0.0, 0.0
        q_nn = (3.0 * m * m - j * (j + 1.0)) / 3.0
        return q_nn, 3.0 * q_nn

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
            q_nn, c_eff = self._tensor_moments(state.m)
            # T_LL = Q_NN P_2(cos theta_S), one line for every spin; for
            # J = 1 this is the HJM (2/3) a_m with a_m = (1/4) c_m
            # (3 cos^2 - 1), digit for digit
            t_geo = TENSOR_LL_SIGN * q_nn * 0.5 * (3.0 * ct * ct - 1.0)
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
        return dsigma_dx_dq2(x, q2, s, f2, r_func=self.r_func)

    def dsigma(self, x, q2, phi, s, state, with_perp=False):
        """Doubly polarized d3sigma/dxdQ2dphi [pb/GeV^2/rad]."""
        t = self.tables(x, q2, with_g2=with_perp)
        w_avg, a_1, a_2 = self.amplitudes(t, x, q2, s, state,
                                          with_perp=with_perp)
        phip = np.asarray(phi, dtype=float) - state.phi_s
        w = 1.0 + w_avg + a_1 * np.cos(phip) + a_2 * np.cos(2.0 * phip)
        return self.dsigma_unpol(x, q2, s) / (2.0 * np.pi) * np.maximum(w, 0.0)
