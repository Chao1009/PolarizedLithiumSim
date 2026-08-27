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

import numpy as np

__all__ = ["NBN", "SI_ACLGAD", "Film", "bethe_mean_ev", "landau_mpv_ev",
           "landau_core_sigma_ev", "hot_spot_nm", "threshold_ratio",
           "saturation_width_nm", "R_S_PROTON_NM", "W_MIP_OPTIMAL_NM",
           "DARK_COUNT_WALL", "SNSPD_THRESHOLD_EV", "landau_xi_ev",
           "landau_density", "landau_sample", "zid_fake_rate",
           "LANDAU_MODE"]

K_MEV = 0.307075                 # MeV mol^-1 cm^2
M_E_EV = 0.510998e6
LANDAU_FWHM_OVER_XI = 4.02       # exact for the Landau distribution

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


# The ePIC Roman-Pot incumbent: an AC-LGAD with a 30 um active layer at
# 500 um pitch, read out by EICROC's 8-bit charge ADC.  It is the thing a
# near-beam nanowire would have to beat, and it does not (plans/09 9.2).
SI_ACLGAD = Film(thickness_nm=30e3, density=2.33, z_over_a=0.4993,
                 i_ev=173.0, delta=5.604,   # Sternheimer, Si at beta*gamma 148
                 name="Si 30 um (ePIC AC-LGAD active)")


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


def landau_mpv_ev(z, beta_gamma, film=NBN):
    """Most probable energy loss in the film [eV], PDG Eq. 34.12.

    ONLY valid where xi >> I.  A nanometre superconducting film is not
    (NBN.xi_over_i() ~ 2e-3), which is why the near-beam study quotes the
    Bethe mean there; the 30 um active layer of an ePIC AC-LGAD is
    (SI_ACLGAD.xi_over_i() ~ 3), which is why the incumbent comparison
    quotes the MPV.  Raises if asked outside the regime."""
    if film.xi_over_i(1, 1.0) < 1.0:
        raise ValueError(
            "Landau MPV is not defined for %s: xi/I = %.2e, far below the "
            "xi >> I regime the formalism needs.  Use bethe_mean_ev."
            % (film.name or "this film", film.xi_over_i(1, 1.0)))
    beta_gamma = float(beta_gamma)
    beta2 = beta_gamma ** 2 / (1.0 + beta_gamma ** 2)
    xi = landau_xi_ev(z, math.sqrt(beta2), film)
    return xi * (math.log(2.0 * M_E_EV * beta_gamma ** 2 / film.i_ev)
                 + math.log(xi / film.i_ev) + 0.200 - beta2 - film.delta)


def landau_core_sigma_ev(z, film=NBN):
    """Gaussian-equivalent core width of the Landau: FWHM/2.355 with the
    exact FWHM = 4.02 xi.  The Landau's upper tail is heavier than
    Gaussian, so a separation quoted in these sigmas OVERSTATES how well
    a light species is kept out of a heavy one's peak."""
    return LANDAU_FWHM_OVER_XI * landau_xi_ev(z, 1.0, film) / 2.355


# --- the Landau itself, because a sigma is the wrong figure of merit -----
#
# A separation quoted as "gap / quadrature sum of core widths" is not a PID
# separation power and is not a fake rate: the Landau's upper tail is what
# puts an alpha inside a 6Li's acceptance window, and the tail falls as
# 1/lambda, not as a Gaussian.  What follows samples the real distribution
# so that the comparison between a THRESHOLD device (one bit per plane) and
# an ADC (many bits) can be stated as a fake rate at a matched efficiency,
# which is unambiguous.

LANDAU_MODE = -0.2228            # argmax of the Landau density, reference


def landau_density(lam):
    """Landau density phi(lambda) = (1/pi) Int_0^inf exp(-u ln u - lambda u)
    sin(pi u) du, integrated half-period by half-period so the oscillating
    integrand converges.  Reproduces the reference mode (-0.2228), FWHM
    (4.02) and phi(0) = 0.17805."""
    from scipy import integrate
    lam = np.atleast_1d(np.asarray(lam, dtype=float))
    out = np.empty_like(lam)
    for i, l in enumerate(lam):
        tot = 0.0
        for a in range(0, 60):
            v, _ = integrate.quad(
                lambda u: math.exp(-u * math.log(u) - l * u) * math.sin(math.pi * u)
                if u > 0 else 0.0, a, a + 1, limit=200)
            tot += v
            if a > 5 and abs(v) < 1e-15:
                break
        out[i] = tot / math.pi
    return out


def _landau_cdf_grid(lam_lo=-4.0, lam_hi=400.0, n=3000):
    """(lambda, CDF) on a log-stretched grid.  Cached: building it costs a
    few seconds of quadrature and the result is universal."""
    global _LANDAU_GRID
    try:
        return _LANDAU_GRID
    except NameError:
        pass
    lam = np.concatenate([np.linspace(lam_lo, 20.0, n),
                          np.geomspace(20.0, lam_hi, n // 4)[1:]])
    dens = landau_density(lam)
    dens = np.clip(dens, 0.0, None)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1])
                                           * np.diff(lam))])
    cdf /= cdf[-1]
    _LANDAU_GRID = (lam, cdf)
    return _LANDAU_GRID


def landau_pdf(lam):
    """Landau density, interpolated off the cached grid.  Fast enough for a
    per-event likelihood; landau_density() is the slow exact version."""
    g_lam, _ = _landau_cdf_grid()
    global _LANDAU_PDF
    try:
        dens = _LANDAU_PDF
    except NameError:
        dens = np.clip(landau_density(g_lam), 0.0, None)
        _LANDAU_PDF = dens
    return np.interp(np.asarray(lam, dtype=float), g_lam, dens,
                     left=0.0, right=0.0)


def landau_sample(n, rng=None):
    """n draws of the Landau variable lambda, by inverse CDF."""
    rng = rng or np.random.default_rng(20260827)
    lam, cdf = _landau_cdf_grid()
    return np.interp(rng.uniform(size=n), cdf, lam)


def zid_fake_rate(z_signal, z_background, film=SI_ACLGAD, beta_gamma=147.6,
                  n_planes=4, efficiency=0.95, readout="llr", n_mc=300000,
                  plane_efficiency=1.0, rng=None):
    """Fake rate of `z_background` faking `z_signal`, at a MATCHED signal
    efficiency, over `n_planes` independent samples of `film`.

    `readout` picks the estimator, which matters far more here than the
    number of bits:
      "llr"       per-plane Landau log-likelihood ratio, summed.  The
                  Neyman-Pearson optimum, and therefore the ceiling any
                  readout can reach.
      "threshold" ONE BIT per plane -- each plane fires or does not -- with
                  a majority-of-k decision, k and the threshold chosen
                  together to minimise the fake at the stated efficiency.
                  This is what a superconducting nanowire can deliver.
      "sum"       the plain sum of the plane deposits.  Included because it
                  is the naive choice and it is BAD: a Landau has no mean,
                  so one delta ray in one plane drags the sum, and the sum
                  loses to a one-bit coincidence.
      "trunc"     truncated mean, dropping the largest plane.  The standard
                  dE/dx estimator, and the reason the naive comparison
                  "8 bits beats 1 bit" does not hold.

    `plane_efficiency` is the probability that a plane records the track at
    all -- ~0.99 for a silicon pixel plane, but only the GEOMETRIC FILL
    FACTOR for a superconducting wire comb: 0.25 (arXiv:2510.11725), 0.40
    (arXiv:2410.00251), 0.50 (the ANL EIC-targeted device, arXiv:2312.13405).
    This is where the two technologies actually part company, and it is a
    fabrication number rather than an information-theoretic one.

    Returns (fake_rate, efficiency_achieved).  If the device cannot reach
    the requested efficiency at any working point, returns (nan, best_eff).

    CAVEAT: the Landau is used untruncated, so the background's upper tail
    is OVERSTATED -- a delta ray above ~50 keV escapes a 30 um layer and
    never deposits.  Every fake rate here is therefore a conservative upper
    bound, and the more tail-sensitive the estimator, the more conservative.
    """
    rng = rng or np.random.default_rng(20260827)
    par = {z: (landau_xi_ev(z, 1.0, film),
               landau_mpv_ev(z, beta_gamma, film))
           for z in (z_signal, z_background)}

    def draw(z):
        xi, mpv = par[z]
        lam = landau_sample(n_mc * n_planes, rng).reshape(n_mc, n_planes)
        return mpv + xi * (lam - LANDAU_MODE)

    sig, bkg = draw(z_signal), draw(z_background)
    if plane_efficiency < 1.0:
        # a plane that does not record the track can never be above threshold
        sig = np.where(rng.uniform(size=sig.shape) < plane_efficiency,
                       sig, -np.inf)
        bkg = np.where(rng.uniform(size=bkg.shape) < plane_efficiency,
                       bkg, -np.inf)

    if readout == "threshold":
        best, reach = (1.0, 0.0), 0.0
        for k in range(1, n_planes + 1):
            lo, hi = 0.0, float(sig.max())
            for _ in range(60):
                t = 0.5 * (lo + hi)
                if float(((sig > t).sum(axis=1) >= k).mean()) >= efficiency:
                    lo = t
                else:
                    hi = t
            eff = float(((sig > lo).sum(axis=1) >= k).mean())
            reach = max(reach, eff)
            if eff + 1e-9 < efficiency:
                continue
            fake = float(((bkg > lo).sum(axis=1) >= k).mean())
            if fake < best[0]:
                best = (fake, eff)
        if best[1] == 0.0:
            return float("nan"), reach
        return best

    def statistic(e):
        if readout == "sum":
            return e.sum(axis=1)
        if readout == "trunc":
            return np.sort(e, axis=1)[:, :-1].sum(axis=1)
        if readout == "llr":
            tot = np.zeros(e.shape[0])
            for z, sign in ((z_signal, 1.0), (z_background, -1.0)):
                xi, mpv = par[z]
                pdf = landau_pdf((e - mpv) / xi + LANDAU_MODE) / xi
                tot += sign * np.log(np.clip(pdf, 1e-300, None)).sum(axis=1)
            return tot
        raise ValueError("unknown readout %r" % (readout,))

    s_stat, b_stat = statistic(sig), statistic(bkg)
    cut = float(np.quantile(s_stat, 1.0 - efficiency))
    return float((b_stat > cut).mean()), float((s_stat > cut).mean())


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
