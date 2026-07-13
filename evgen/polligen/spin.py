"""Spin-density matrices and multipole moments for spin-1 and spin-3/2 ions.

The beam is an incoherent mixture of |J m> states along a quantization axis
n(theta_S, phi_S): populations p_m with sum(p_m) = 1, ordered m = +J ... -J
throughout this package.  The lab-frame density matrix is

    rho_lab = U rho_diag U^dagger,   U_{m'm} = D^J_{m'm}(phi_S, theta_S, 0)
                                            = exp(-i m' phi_S) d^J_{m'm}(theta_S)

Normalized moments (dimensionless; conventions used across polligen):

  vector          P    = <J_z>/J                      any J, in [-1, 1]
  tensor  J=1     Pzz  = <3 J_z^2 - 2>  = p+ + p- - 2 p0,   in [-2, +1]
  tensor  J=3/2   T    = <3 J_z^2 - J(J+1)>/3
                       -> +1 for |m|=3/2 pure, -1 for |m|=1/2 pure
  octupole J=3/2  O    = <J_z^3 - (41/20) J_z>/(3/10)
                       -> +1 for pure m=+3/2   (rank-3; carried but unused
                          by the cross section, plans/05 truncation)

The J=1 conventions match Hoodbhoy-Jaffe-Manohar NPB 312:571 as used in
`polli_fastsim.asymmetries` (c_m = 3 m^2 - 2).  Wigner-d / Clebsch-Gordan
follow the standard (Varshalovich / Wikipedia) conventions.
"""

from dataclasses import dataclass
from math import factorial

import numpy as np


def _fact(n):
    """Factorial of a value that must be a non-negative integer."""
    k = int(round(n))
    if abs(n - k) > 1e-9 or k < 0:
        raise ValueError("factorial argument %r is not a non-negative integer"
                         % (n,))
    return factorial(k)


def m_values(j):
    """Spin projections ordered +J ... -J (package-wide convention)."""
    n = int(round(2 * j + 1))
    return np.array([j - i for i in range(n)])


def clebsch_gordan(j1, m1, j2, m2, j, m):
    """<j1 m1 j2 m2 | j m> via the Racah formula (exact for small spins)."""
    if abs(m1 + m2 - m) > 1e-9:
        return 0.0
    if not (abs(j1 - j2) - 1e-9 <= j <= j1 + j2 + 1e-9):
        return 0.0
    pre = ((2 * j + 1)
           * _fact(j1 + j2 - j) * _fact(j1 - j2 + j) * _fact(-j1 + j2 + j)
           / _fact(j1 + j2 + j + 1))
    pre *= (_fact(j1 + m1) * _fact(j1 - m1) * _fact(j2 + m2) * _fact(j2 - m2)
            * _fact(j + m) * _fact(j - m))
    k_min = int(round(max(0.0, j2 - j - m1, j1 - j + m2)))
    k_max = int(round(min(j1 + j2 - j, j1 - m1, j2 + m2)))
    tot = 0.0
    for k in range(k_min, k_max + 1):
        tot += ((-1.0) ** k
                / (_fact(k) * _fact(j1 + j2 - j - k) * _fact(j1 - m1 - k)
                   * _fact(j2 + m2 - k) * _fact(j - j2 + m1 + k)
                   * _fact(j - j1 - m2 + k)))
    return np.sqrt(pre) * tot


def wigner_d(j, beta):
    """Small Wigner matrix d^j_{m'm}(beta), indexed [m', m] with the
    package m-ordering (+J ... -J)."""
    ms = m_values(j)
    n = ms.size
    d = np.zeros((n, n))
    c, s = np.cos(beta / 2.0), np.sin(beta / 2.0)
    for a, mp in enumerate(ms):
        for b, m in enumerate(ms):
            pre = np.sqrt(_fact(j + mp) * _fact(j - mp)
                          * _fact(j + m) * _fact(j - m))
            k_min = int(round(max(0.0, m - mp)))
            k_max = int(round(min(j + m, j - mp)))
            tot = 0.0
            for k in range(k_min, k_max + 1):
                num = (-1.0) ** (mp - m + k)
                den = (_fact(j + m - k) * _fact(k)
                       * _fact(j - mp - k) * _fact(mp - m + k))
                tot += (num / den
                        * c ** (2 * j + m - mp - 2 * k)
                        * s ** (mp - m + 2 * k))
            d[a, b] = pre * tot
    return d


def rotation_matrix(j, theta, phi):
    """U_{m'm} = exp(-i m' phi) d^j_{m'm}(theta): |jm>_axis = sum U|jm'>_z."""
    ms = m_values(j)
    return np.exp(-1j * ms[:, None] * phi) * wigner_d(j, theta)


def angular_momentum_ops(j):
    """(J_x, J_y, J_z) matrices in the package m-ordering."""
    ms = m_values(j)
    n = ms.size
    jz = np.diag(ms).astype(complex)
    jp = np.zeros((n, n), dtype=complex)  # raising: <m+1|J+|m>
    for b, m in enumerate(ms):
        if b > 0:  # m+1 is at index b-1 (descending order)
            jp[b - 1, b] = np.sqrt(j * (j + 1) - m * (m + 1))
    jm = jp.conj().T
    jx = 0.5 * (jp + jm)
    jy = -0.5j * (jp - jm)
    return jx, jy, jz


def multipole_operator(j, k, q):
    """Irreducible tensor operator T_kq with
    <j m'|T_kq|j m> = <j m; k q|j m'> sqrt((2k+1)/(2j+1))."""
    ms = m_values(j)
    n = ms.size
    t = np.zeros((n, n), dtype=complex)
    norm = np.sqrt((2 * k + 1) / (2 * j + 1))
    for a, mp in enumerate(ms):
        for b, m in enumerate(ms):
            t[a, b] = clebsch_gordan(j, m, k, q, j, mp) * norm
    return t


def rho_from_populations(j, populations, theta=0.0, phi=0.0):
    """Lab-frame density matrix for populations p_m along axis n(theta, phi)."""
    p = np.asarray(populations, dtype=float)
    if p.size != int(round(2 * j + 1)):
        raise ValueError("populations must have 2j+1 entries (m=+J..-J)")
    if np.any(p < -1e-12) or abs(p.sum() - 1.0) > 1e-9:
        raise ValueError("populations must be >= 0 and sum to 1")
    u = rotation_matrix(j, theta, phi)
    return u @ np.diag(p).astype(complex) @ u.conj().T


def vector_polarization(rho, j):
    """<J>/J as a real 3-vector (lab frame)."""
    jx, jy, jz = angular_momentum_ops(j)
    return np.array([np.trace(rho @ op).real for op in (jx, jy, jz)]) / j


def tensor_polarization(rho, j):
    """Normalized J_z tensor moment: J=1 -> Pzz = <3Jz^2-2>;
    J=3/2 -> T = <3Jz^2 - J(J+1)>/3."""
    _, _, jz = angular_momentum_ops(j)
    q = np.trace(rho @ (3.0 * jz @ jz - j * (j + 1) * np.eye(jz.shape[0]))).real
    if abs(j - 1.0) < 1e-9:
        return q          # <3Jz^2> - 2
    if abs(j - 1.5) < 1e-9:
        return q / 3.0
    raise ValueError("tensor_polarization defined for j = 1, 3/2 only")


def octupole_moment(rho, j=1.5):
    """Normalized rank-3 moment for J=3/2: <Jz^3 - (41/20) Jz>/(3/10)."""
    if abs(j - 1.5) > 1e-9:
        raise ValueError("octupole moment defined for j = 3/2 only")
    _, _, jz = angular_momentum_ops(j)
    op = jz @ jz @ jz - (41.0 / 20.0) * jz
    return np.trace(rho @ op).real / 0.3


def moments_along_axis(j, populations):
    """(vector, tensor[, octupole]) moments of p_m along its own axis."""
    p = np.asarray(populations, dtype=float)
    ms = m_values(j)
    vec = float((ms * p).sum() / j)
    if abs(j - 1.0) < 1e-9:
        ten = float(((3.0 * ms**2 - 2.0) * p).sum())
        return vec, ten
    if abs(j - 1.5) < 1e-9:
        ten = float(((3.0 * ms**2 - j * (j + 1)) * p).sum() / 3.0)
        octu = float(((ms**3 - (41.0 / 20.0) * ms) * p).sum() / 0.3)
        return vec, ten, octu
    raise ValueError("moments defined for j = 1, 3/2 only")


def spin1_populations(pz, pzz):
    """(p+, p0, p-) with <Jz> = pz and <3Jz^2-2> = pzz.
    Physical domain: pzz in [-2, 1], and all p_m >= 0."""
    p_plus = (2.0 + 3.0 * pz + pzz) / 6.0
    p_zero = (1.0 - pzz) / 3.0
    p_minus = (2.0 - 3.0 * pz + pzz) / 6.0
    pops = np.array([p_plus, p_zero, p_minus])
    if np.any(pops < -1e-12):
        raise ValueError("unphysical (pz=%g, pzz=%g): populations %s"
                         % (pz, pzz, pops))
    return tuple(np.clip(pops, 0.0, 1.0))


def spin32_populations(pz, t=0.0, o=0.0):
    """(p_{3/2}, p_{1/2}, p_{-1/2}, p_{-3/2}) with normalized moments
    (vector pz, tensor t, octupole o) as defined in this module."""
    ms = m_values(1.5)
    a = np.vstack([
        np.ones(4),
        ms / 1.5,                       # vector, normalized
        (3.0 * ms**2 - 3.75) / 3.0,     # tensor, normalized
        (ms**3 - (41.0 / 20.0) * ms) / 0.3,  # octupole, normalized
    ])
    pops = np.linalg.solve(a, np.array([1.0, pz, t, o]))
    if np.any(pops < -1e-12):
        raise ValueError("unphysical (pz=%g, t=%g, o=%g): populations %s"
                         % (pz, t, o, pops))
    return tuple(np.clip(pops, 0.0, 1.0))


def populations_maxent(j, pz, beta_max=60.0, tol=1e-13):
    """Spin-temperature (maximum-entropy) populations p_m ~ exp(beta m)
    with <J_z>/J = pz -- the standard population model of polarized ion
    sources.  Works for any j; solves for beta by bisection.  Note that
    for j >= 1 this implies non-zero higher moments (e.g. the spin-1
    spin-temperature relation between pz and pzz)."""
    if not -1.0 < pz < 1.0:
        raise ValueError("maxent populations need |pz| < 1")
    ms = m_values(j)

    def mean(beta):
        w = np.exp(beta * (ms - ms.max()))  # stable for large beta
        return float((ms * w).sum() / w.sum() / j)

    lo, hi = -beta_max, beta_max
    if not mean(lo) < pz < mean(hi):
        raise ValueError("pz=%g out of reach for beta in [%g, %g]"
                         % (pz, lo, hi))
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if mean(mid) < pz:
            lo = mid
        else:
            hi = mid
    beta = 0.5 * (lo + hi)
    w = np.exp(beta * (ms - ms.max()))
    return tuple(w / w.sum())


@dataclass(frozen=True)
class SpinDensity:
    """Populations along an axis + the axis orientation, with lab moments."""
    j: float
    populations: tuple        # m = +J ... -J along the axis
    theta: float = 0.0        # axis polar angle in the lab [rad]
    phi: float = 0.0          # axis azimuth in the lab [rad]

    def rho_lab(self):
        return rho_from_populations(self.j, self.populations,
                                    self.theta, self.phi)

    def lab_moments(self):
        rho = self.rho_lab()
        out = {"vector": vector_polarization(rho, self.j),
               "tensor_zz": tensor_polarization(rho, self.j)}
        if abs(self.j - 1.5) < 1e-9:
            out["octupole_z"] = octupole_moment(rho, self.j)
        return out
