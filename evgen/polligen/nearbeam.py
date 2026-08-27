"""Near-beam sensor physics for the far-forward lithium tags (plans/09).

Two independent pieces, both used by the near-beam study and its report
(`reports/nanowire_far_forward`) and both anchored on published numbers:

1.  How much energy a relativistic ion leaves in a thin superconducting
    film.  The Bethe mean is quoted and the Landau most-probable value is
    NOT: a 12 nm NbN film has xi/I ~ 2e-3, far outside the xi >> I regime
    the most-probable-value formalism requires, so only the mean and the
    z^2 scaling are defensible here.

2.  The hot-spot firing threshold, which is how a superconducting
    nanowire carries charge information at all.  The device latches --
    its pulse amplitude is the diverted BIAS current and is the same for
    120 GeV hadrons, 120 GeV muons, 8 GeV pions and 8 GeV showering
    electrons (measured: arXiv:2510.11725 Sec. 4, arXiv:2410.00251
    Sec. 4.2) -- so no pulse-height scheme can work.  What the deposit
    sets is the bias current at which the wire fires:

        r_s = sqrt(Q / (e pi c rho (T_c - T_0)))
        I_th / I_c = 1 - 2 r_s / w

    Argonne's own Eqs. 1-2 (Lee et al., arXiv:2312.13405, NIM A 1069
    (2024) 169956).  Two cautions the paper states about itself and one
    the literature states about it:

      * r_s = 134 nm is the EXTRAPOLATED zero-crossing of a four-point
        fit, not a measurement.  The authors write "While the physical
        validity of this simple model is a question of future work".
        Inverting their four Fig. 7 points individually gives 102-118 nm
        (mean 113) -- see tests/test_nearbeam.py.
      * The VOLUME in which Q is counted is unresolved by a factor 5000.
        The 2024 paper defines Q as "the energy that the proton has
        deposited into the thin film" (~20 eV by the Bethe mean below);
        the 2026 paper, with overlapping authorship, applies the same
        sqrt(Q) scaling to the SUBSTRATE deposit (0.1 MeV for the same
        proton, arXiv:2601.03158).  Both scale as z^2, so r_s ~ z
        survives; the absolute radii do not.
      * Their 5.5 MeV 241Am alpha is NOT a z-scaling calibration point.
        It differs from the 120 GeV proton almost entirely through
        1/beta^2 (beta = 0.054 against ~1), so it tests the sqrt(Q) law
        across energy and says nothing about z^2 at fixed velocity.
        Nobody has varied Z at fixed beta on one of these devices.  The same relation appears on the photon side as
    E = (w/C)^2 (1 - I_b/I_c)^2 (Renema et al., arXiv:1301.3337) and in
    the ion literature as I_th/I_c = 1 - (zeV)^(1/2) C/w (Cristiano et
    al., Supercond. Sci. Technol. 28 (2015) 124004).

    Since dE/dx goes as z^2 at fixed velocity, r_s goes as z LINEARLY.
    A 6Li at 137.5 GeV/u has beta*gamma = 148 against 128 for Argonne's
    120 GeV calibration proton -- the same velocity to 15% -- so the step
    from their r_s = 134 nm to z = 3 needs no velocity correction.

The consequence that decides the design: below w = 2 r_s the threshold is
zero, the wire fires at any bias, and it carries no charge information at
all.  Argonne's stated efficiency optimum for a MIP, w ~ 250 nm, is
exactly 2 r_s for a proton -- maximal efficiency and zero charge
information -- so Z identification needs deliberately WIDE wires, which
is the microwire (SMSPD) geometry that already exists at w = 1 um.

WHERE THIS ENDS UP (reports/nanowire_far_forward, 2026-08-26): a
threshold nanowire delivers ONE BIT per plane, while the incumbent ePIC
AC-LGAD digitises an 8-bit charge over a 30 um active layer -- which
separates 6Li from alpha at ~4.8 sigma per plane, over the four planes
ePIC already has.  A nanowire cannot beat that observable.  This module
is kept because the threshold physics is correct, testable and worth
recording; it is not kept because the design won.
"""

import math

__all__ = ["NBN", "Film", "bethe_mean_ev", "hot_spot_nm", "threshold_ratio",
           "saturation_width_nm", "R_S_PROTON_NM", "W_MIP_OPTIMAL_NM",
           "DARK_COUNT_WALL", "SNSPD_THRESHOLD_EV", "landau_xi_ev"]

K_MEV = 0.307075                 # MeV mol^-1 cm^2
M_E_EV = 0.510998e6

R_S_PROTON_NM = 134.0            # EXTRAPOLATED, 120 GeV p in 12 nm NbN
                                 # -- a fit parameter, not a datum: see the
                                 # module docstring
W_MIP_OPTIMAL_NM = 250.0         # ANL's stated optimum for a MIP (= 2 r_s)
DARK_COUNT_WALL = 0.80           # I_b/I_c above which background rises
                                 # exponentially (arXiv:2312.13405 Sec. 3.1)
SNSPD_THRESHOLD_EV = 0.8         # one 1.55 um photon


class Film:
    """A thin absorber: thickness, density, <Z/A> and mean excitation
    energy.  `delta` is the density-effect correction, held flat on the
    Fermi plateau; it enters every species as the same xi * delta, so the
    species RATIOS are insensitive to it."""

    def __init__(self, thickness_nm, density, z_over_a, i_ev, delta=5.7,
                 name=""):
        self.thickness_nm = float(thickness_nm)
        self.density = float(density)
        self.z_over_a = float(z_over_a)
        self.i_ev = float(i_ev)
        self.delta = float(delta)
        self.name = name

    @property
    def mass_thickness(self):
        """g/cm^2."""
        return 1e-7 * self.thickness_nm * self.density

    def xi_over_i(self, z=1, beta=1.0):
        """The Landau validity ratio.  The most-probable-value formalism
        needs xi >> I; a nanometre film is nowhere near it."""
        return landau_xi_ev(z, beta, self) / self.i_ev


def _bragg_i_ev(pairs):
    """Mean excitation energy by Bragg additivity over electron fraction."""
    ztot = sum(z for z, _ in pairs)
    return math.exp(sum(z * math.log(i) for z, i in pairs) / ztot)


# The beam-tested Argonne device: 12 nm NbN (arXiv:2312.13405).
NBN = Film(thickness_nm=12.0, density=8.47,
           z_over_a=(41.0 + 7.0) / (92.906 + 14.007),
           i_ev=_bragg_i_ev(((41.0, 417.0), (7.0, 82.0))),
           name="NbN 12 nm")


def landau_xi_ev(z, beta, film=NBN):
    """Landau width parameter xi = (K/2)(Z/A)(z/beta)^2 x, in eV.  Kept
    for the validity check only -- see the module docstring."""
    return (1e6 * 0.5 * K_MEV * film.z_over_a * (z / beta) ** 2
            * film.mass_thickness)


def bethe_mean_ev(z, beta_gamma, film=NBN):
    """Mean energy loss of a charge-z ion in the film [eV].

    The unrestricted Bethe mean with the density-effect correction of
    `film`.  Delta rays escape a nanometre film, so this OVERSTATES what
    is deposited in it; what the design argument uses is the z^2 ratio,
    which is unaffected."""
    beta_gamma = float(beta_gamma)
    beta2 = beta_gamma ** 2 / (1.0 + beta_gamma ** 2)
    t_max = 2.0 * M_E_EV * beta_gamma ** 2          # eV, no mass correction
    coeff = 1e6 * K_MEV * film.z_over_a * z ** 2 / beta2 * film.mass_thickness
    return coeff * (0.5 * math.log(2.0 * M_E_EV * beta_gamma ** 2 * t_max
                                   / film.i_ev ** 2)
                    - beta2 - 0.5 * film.delta)


def hot_spot_nm(z, r_s_proton=R_S_PROTON_NM):
    """Normal-core radius r_s = sqrt(Q / (e pi c rho (T_c - T_0))).

    Q goes as z^2 at fixed velocity, so r_s goes as z.  Anchored on
    Argonne's MEASURED 134 nm for a 120 GeV proton, which is why this is
    a one-parameter scaling and not a thermodynamic calculation."""
    return float(r_s_proton) * abs(float(z))


def threshold_ratio(r_s_nm, w_nm):
    """I_th / I_c = 1 - 2 r_s / w, clipped at zero.

    Zero means the hot spot already spans the wire at zero bias: the wire
    fires on anything and carries NO charge information."""
    if w_nm <= 0.0:
        raise ValueError("wire width must be positive")
    return max(0.0, 1.0 - 2.0 * float(r_s_nm) / float(w_nm))


def saturation_width_nm(z, r_s_proton=R_S_PROTON_NM):
    """The wire width below which charge z carries no information,
    w = 2 r_s.  For a proton this is ~268 nm -- which is why Argonne's
    MIP-efficiency optimum of ~250 nm has zero Z discrimination."""
    return 2.0 * hot_spot_nm(z, r_s_proton)
