"""Leading-log initial-state radiation as a kinematic migration
(plans/07 WP4, plans/08 D3, plans/02 step 1.4).

WHAT THIS MODULE IS.  Not a radiative-correction calculation: a *bound*.
The letter's inclusive cos 2phi' measurement needs to know how much a
QED photon radiated off the incoming electron can move the numbers the
chain quotes -- purity, efficiency, the bin-centering factor K and,
through K, Delta.  Two of the four things WP4 originally planned are
void by algebra rather than by simulation (plans/08 D3):

  * a collinear photon generates EXACTLY ZERO fake cos phi' or cos 2phi'.
    The covariant azimuth of `reco.azimuth_wrt_lepton_plane` is invariant
    under k -> (1 - z) k for a massless target: writing a = P.k,
    b = P.k', c = k.k' and l = q = (1-z)k - k', both the cosine and the
    sine of the azimuth pick up the SAME factor [2ac((1-z)a - b)]^(-1/2),
    which the arctan2 divides out (proof in the module test
    `test_covariant_azimuth_is_invariant_under_a_collinear_photon`);
  * the spin-state ratio cancels any radiative effect COMMON to the
    fills, which an unpolarized-lepton QED correction is by construction
    (the fills differ only in the ion's tensor population).

What survives is the migration, and that is what this module measures.

THE SPECTRUM.  The leading-log electron structure function with the
soft-photon part exponentiated (Kuraev-Fadin, Sov. J. Nucl. Phys. 41
(1985) 466; Nicrosini-Trentadue, Phys. Lett. B196 (1987) 551; the form
quoted as Eq. (2) of Skrzypek, Acta Phys. Polon. B23 (1992) 135 and used
by every HERA radiative-correction study), written for the RADIATED
fraction z = 1 - x_e:

    D(z, Q2) = (t/2) z^(t/2 - 1) S(t)  -  (t/4) (2 - z),
    S(t)     = exp[(t/2)(3/4 - gamma_E)] / Gamma(1 + t/2),
    t        = (2 alpha / pi) [ ln(Q2 / m_e^2) - 1 ].

The first term resums the soft/virtual series into an integrable
z^(t/2-1) endpoint singularity; the second is the O(alpha) hard
collinear remainder.  int_0^1 D dz = S - 3t/8 = 1 + O(t^2) --
`total_weight` returns exactly that residual, and the module test pins
it below t^2.  t = 0.071 at Q2 = 3 GeV^2, so the normalisation is good
to 7e-4 and <z> = 2.3e-2.

THE MIGRATION.  A pseudo-event of the response carries the drawn
(x, Q2) and y = Q2/(s x).  Reading it as a radiative event means
reading (x, Q2) as the HARD variables at the reduced invariant
ss = (1 - z) s, so that

    y_hard = Q2 / (x ss) = y / (1 - z),   E_hard = (1 - z) E_e,

and the scattered electron is the one that reduced beam makes, which is
NOT the one the same (x, Q2) makes at the nominal beam.  Q2 is invariant
under this reading, which is why it is the reading the response uses:
the sampler cell, its cross-section weight and its cos 2phi' amplitude
table are all keyed on (x, Q2) and do not move.  The kinematic ceiling
y_hard <= 1 caps z at 1 - y.

The analysis does not know z.  It reconstructs with the NOMINAL beam and
the nominal s, and `observed_kinematics` below gives what each method
then returns; every entry is an identity, tested against a full
four-vector construction through `hfs.hadronic_kinematics`:

  method              Q2_obs / Q2_hard    y_obs / y_hard   x_obs / x_hard
  electron            1/(1-z)             (y+z)/y_hard     (see below)
  Sigma               1                   1                1 - z
  Jacquet-Blondel     (1-y_hard)/(1-y)    1 - z            ...
  double angle        1/(1-z)^2           1                1/(1-z)
  MIXED (the chain)   1/(1-z)  [Q2_e]     1  [y_Sigma]     1   EXACT

So the chain's x = Q2_e/(s y_Sigma) is exact under a collinear photon --
the (1 - z) cancels between Q2_e and the Sigma-method y, both of which
are off by it in the same direction -- and the ONLY thing that migrates
is the Q2 LABEL of the bin, upward by 1/(1 - z).  (This corrects the
"x -> x/(1-z)" statement that stood in Report 2 section 3 and in
plans/08 D3 until 2026-08-28; code review R16.)  The Sigma method's own
Q2 = pT_e^2/(1 - y_Sigma) uses no beam energy either and is exact, so a
chain that labelled its bins with Q2_Sigma instead of Q2_e would have no
ISR migration at all -- at the price of the Sigma resolution on the
label.  That is a chain change, not made here; the bound below is for
the chain as published.

WHAT IS AND IS NOT MODELLED.  Modelled: collinear ISR off the incoming
electron, its kinematic migration, the (1 - z) shift of the hard
unpolarized rate (`rate_scale`) and of the cos 2phi' amplitude
(`amplitude_scale`).  NOT modelled, and stated wherever the bound is
quoted: the TENSOR-SECTOR radiative correction, which has never been
calculated (plans/05 section 5.5, plans/07 systematics row 5) -- this
module bounds the migration of an unpolarized QED correction, not the
size of an unknown tensor one; the polarized-LEPTON correction (the
beam is unpolarized in this measurement); non-collinear (wide-angle)
real emission, which is O(alpha) without the log and does not factorise
into a beam-energy shift; final-state radiation off the scattered
electron, which a calorimeter cluster largely reabsorbs; and the elastic
and quasi-elastic radiative tails, which the W2 >= 10 GeV^2 cut removes
from the inclusive sample.
"""

import numpy as np
from scipy.special import gamma as _gamma_fn

from polli_fastsim.structure import ALPHA_EM

#: electron mass [GeV]
M_E = 0.51099895e-3
#: Euler-Mascheroni constant
EULER_GAMMA = 0.5772156649015329


# --- the leading-log spectrum -------------------------------------------------

def beta_ll(q2):
    """t = (2 alpha/pi) [ln(Q2/m_e^2) - 1], the LL radiator exponent."""
    q2 = np.asarray(q2, dtype=float)
    return (2.0 * ALPHA_EM / np.pi) * (np.log(np.maximum(q2, 4.0 * M_E ** 2)
                                              / M_E ** 2) - 1.0)


def soft_factor(t):
    """S(t) = exp[(t/2)(3/4 - gamma_E)] / Gamma(1 + t/2)."""
    t = np.asarray(t, dtype=float)
    return np.exp(0.5 * t * (0.75 - EULER_GAMMA)) / _gamma_fn(1.0 + 0.5 * t)


def d_electron(z, q2):
    """D(z, Q2), the density of the RADIATED fraction z (module docstring).

    Positive on (0, 1] for the t <~ 0.1 of interest; the z -> 0 endpoint
    singularity is integrable (z^(t/2-1))."""
    z = np.asarray(z, dtype=float)
    t = beta_ll(q2)
    return (0.5 * t * np.power(np.maximum(z, 1e-300), 0.5 * t - 1.0)
            * soft_factor(t) - 0.25 * t * (2.0 - z))


def cdf_z(z, q2):
    """int_0^z D(z', Q2) dz' = S z^(t/2) - (t/4)(2z - z^2/2)."""
    z = np.asarray(z, dtype=float)
    t = beta_ll(q2)
    return (soft_factor(t) * np.power(np.maximum(z, 0.0), 0.5 * t)
            - 0.25 * t * (2.0 * z - 0.5 * z * z))


def moment_z(q2, z_max=1.0, z_min=0.0, order=1):
    """int_{z_min}^{z_max} z^order D dz in closed form (order 0, 1 or 2)."""
    t = np.asarray(beta_ll(q2), dtype=float)
    s_t = soft_factor(t)

    def prim(z):
        z = np.asarray(z, dtype=float)
        if order == 0:
            return s_t * z ** (0.5 * t) - 0.25 * t * (2.0 * z - 0.5 * z ** 2)
        if order == 1:
            return (s_t * (0.5 * t) * z ** (1.0 + 0.5 * t) / (1.0 + 0.5 * t)
                    - 0.25 * t * (z ** 2 - z ** 3 / 3.0))
        if order == 2:
            return (s_t * (0.5 * t) * z ** (2.0 + 0.5 * t) / (2.0 + 0.5 * t)
                    - 0.25 * t * (2.0 * z ** 3 / 3.0 - 0.25 * z ** 4))
        raise ValueError("order must be 0, 1 or 2")

    return prim(z_max) - prim(z_min)


def total_weight(q2):
    """int_0^1 D dz - 1 = S(t) - 3t/8 - 1 = O(t^2): the normalisation
    residual of the truncated exponentiation, the module's own error."""
    return moment_z(q2, 1.0, 0.0, order=0) - 1.0


# --- the sampler --------------------------------------------------------------

class ISRModel:
    """Per-event collinear-ISR sampler and migration hook.

    Its own random stream.  `RecoResponse` draws z from THIS generator,
    never from the response's rng, so that an ISR-on and an ISR-off
    response built from the same seed share every other random number:
    smearing, hadronic resolution and in-cell positions are identical and
    the difference between the two is the radiation alone, not seed noise
    (the `--syst-scan` convention of `money_cos2phi_reco.py`).

    Parameters
    ----------
    z_min : float
        photons below this radiated fraction are treated as no emission.
        The migration they cause is a relative shift 1/(1 - z) - 1 of the
        Q2 label, i.e. 1e-4 at the default -- three orders below the
        bin widths and below the MC floor.
    y_max_hard : float
        ceiling on y_hard = y/(1 - z); caps z at 1 - y/y_max_hard.  This
        is the kinematic limit of the emission, not a cut.
    reweight : bool
        also apply the (1 - z) shift of the HARD cross section
        (`rate_scale`) to the event weight and of the cos 2phi'
        amplitude (`amplitude_scale`).  Both are exact closed forms in
        the reduced ss; switching them off leaves the pure migration.
    seed : int
        seed of the module's own stream.
    """

    def __init__(self, z_min=1e-4, y_max_hard=0.999, reweight=True,
                 seed=20260828, n_reject=12):
        self.z_min = float(z_min)
        self.y_max_hard = float(y_max_hard)
        self.reweight = bool(reweight)
        self.seed = int(seed)
        self.n_reject = int(n_reject)
        self._rng = np.random.default_rng(self.seed)

    def describe(self):
        return ("exponentiated LL ISR: z_min = %.0e, y_hard <= %.3f, "
                "reweight = %s, seed = %d"
                % (self.z_min, self.y_max_hard, self.reweight, self.seed))

    def reset(self):
        """Rewind the stream (so two responses radiate identically)."""
        self._rng = np.random.default_rng(self.seed)
        return self

    def z_max(self, y):
        """Kinematic ceiling on z from y_hard = y/(1 - z) <= y_max_hard."""
        y = np.asarray(y, dtype=float)
        return np.clip(1.0 - y / self.y_max_hard, 0.0, 1.0)

    def p_emit(self, q2, y):
        """P(z > z_min | z < z_max), the fraction of events that radiate
        a photon this module resolves.

        The spectrum is RENORMALISED over the kinematically allowed
        range, not truncated at it: an event whose hard configuration
        would need y_hard > y_max_hard cannot exist, and the missing mass
        int_{z_max}^1 D dz is put back into (z_min, z_max) instead of
        being removed from the rate.  That over-counts the resolved
        emission by f = [F(1) - F(z_max)] / F(z_max), which is
        conservative for a bound and negligible where the bound is
        quoted: f = 1.3e-3 at y = 0.03 and 9e-4 at y = 0.01 (the sweet
        spots), rising to 1.5e-2 at y = 0.5 and 6.3e-2 at y = 0.9, where
        the rate is already gone (`cdf_z`, t = 0.073 at <Q2> = 4.4)."""
        zmax = self.z_max(y)
        f_hi = cdf_z(zmax, q2)
        f_lo = cdf_z(np.minimum(self.z_min, zmax), q2)
        with np.errstate(divide="ignore", invalid="ignore"):
            p = np.where(f_hi > 0.0, 1.0 - f_lo / np.maximum(f_hi, 1e-300), 0.0)
        return np.clip(np.where(zmax > self.z_min, p, 0.0), 0.0, 1.0)

    def sample(self, q2, y, rng=None):
        """Radiated fraction z per event (0 = no resolved emission).

        Rejection off the exponentiated soft term g(z) = (t/2) z^(t/2-1),
        which majorises D by the constant S(t) because the hard remainder
        enters with a minus sign; g is inverted analytically on
        (z_min, z_max), so the sampler is exact up to the `n_reject`
        cap (the acceptance is > 0.99 wherever the density has weight)."""
        rng = self._rng if rng is None else rng
        q2 = np.asarray(q2, dtype=float)
        y = np.asarray(y, dtype=float)
        q2, y = np.broadcast_arrays(q2, y)
        n = y.size
        t = beta_ll(q2)
        s_t = soft_factor(t)
        zmax = self.z_max(y)
        z = np.zeros(n)
        emit = rng.random(n) < self.p_emit(q2, y)
        if not np.any(emit):
            return z
        half = 0.5 * t[emit]
        g_lo = np.power(np.minimum(self.z_min, zmax[emit]), half)
        g_hi = np.power(zmax[emit], half)
        todo = np.ones(int(emit.sum()), dtype=bool)
        out = np.zeros(todo.size)
        for _ in range(self.n_reject):
            m = todo
            u = rng.random(int(m.sum()))
            cand = np.power(g_lo[m] + u * (g_hi[m] - g_lo[m]),
                            1.0 / half[m])
            # accept with D/(S g) = 1 - (2 - z) z^(1 - t/2) / (2 S)
            acc = 1.0 - ((2.0 - cand) * np.power(cand, 1.0 - half[m])
                         / (2.0 * s_t[emit][m]))
            take = rng.random(cand.size) < acc
            idx = np.flatnonzero(m)
            out[idx[take]] = cand[take]
            todo[idx[take]] = False
            if not todo.any():
                break
        out[todo] = np.minimum(self.z_min, zmax[emit][todo])
        z[emit] = out
        return z

    # --- the two hard-side rescalings ------------------------------------

    def rate_scale(self, kernel, x, q2, y, z):
        return rate_scale(kernel, x, q2, y, z) if self.reweight \
            else np.ones_like(np.asarray(z, dtype=float))

    def amplitude_scale(self, kernel, x, q2, y, z):
        return amplitude_scale(kernel, x, q2, y, z) if self.reweight \
            else np.ones_like(np.asarray(z, dtype=float))


# --- the hard-side (1 - z) rescalings ----------------------------------------

def rate_scale(kernel, x, q2, y, z):
    """d2sigma_unpol(x, Q2, ss)/d2sigma_unpol(x, Q2, s) at ss = (1-z)s.

    At fixed (x, Q2) the Born cross section depends on the beam only
    through y, so this is the ratio of the [(1 - y + y^2/2) F2 -
    (y^2/2) F_L] brackets at y_hard = y/(1-z) and at y.  It is the RC
    weight of the convolution
    dsigma^rad/dxdQ2 = int dz D(z,Q2) dsigma^0(x, Q2, ss(z))/dxdQ2."""
    x = np.asarray(x, dtype=float)
    q2 = np.asarray(q2, dtype=float)
    z = np.asarray(z, dtype=float)
    s = q2 / (x * np.asarray(y, dtype=float))
    den = kernel.dsigma_unpol(x, q2, s)
    num = kernel.dsigma_unpol(x, q2, s * np.maximum(1.0 - z, 1e-12))
    return np.where(den > 0.0, num / np.where(den > 0.0, den, 1.0), 0.0)


def amplitude_scale(kernel, x, q2, y, z):
    """a_2(y_hard) / a_2(y) at fixed (x, Q2), ss = (1-z)s.

    a_2 = -[(1-y)/y^2] c_eff sin^2(theta_S) Delta / D_phi with
    D_phi = F1 + [(1-y)/(x y^2)] F2 (xsec.InclusiveKernel.amplitudes), so
    the whole m-, spin- and Delta-dependence cancels in the ratio and
    what is left is

        g(y) = x / (F2 + F1 x y^2/(1 - y)),

    a pure function of y at fixed (x, Q2).  g is flat at small y (where
    the F2 term dominates D_phi), which is why this rescaling is a
    sub-per-mille effect at the sweet spots and the migration of the Q2
    label is the whole story."""
    x = np.asarray(x, dtype=float)
    q2 = np.asarray(q2, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    y_hard = np.clip(y / np.maximum(1.0 - z, 1e-12), 0.0, 1.0 - 1e-12)
    f2 = kernel.nf2.f2a(x, q2) / kernel.ion.A
    f1 = kernel.nf2.f1a(x, q2) / kernel.ion.A

    def g(yy):
        return x / np.maximum(f2 + f1 * x * yy * yy / np.maximum(1.0 - yy,
                                                                 1e-12), 1e-300)

    return g(y_hard) / np.maximum(g(y), 1e-300)


# --- what each reconstruction method sees ------------------------------------

def observed_kinematics(x, q2, y, z, s):
    """Closed-form observed kinematics of a radiative event.

    `(x, q2, y)` are the drawn nominal-s variables of a pseudo-event and
    `z` the radiated fraction; the HARD event is (x, q2, y_hard) at
    ss = (1-z)s (module docstring).  Everything is returned as the value
    the analysis obtains with the NOMINAL beam and the nominal s, next to
    the hard truth, so that a caller can read the bias of each method
    straight off.  Verified against a four-vector construction through
    `hfs.hadronic_kinematics` in `tests/test_radiative.py`."""
    x = np.asarray(x, dtype=float)
    q2 = np.asarray(q2, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    one_z = 1.0 - z
    y_hard = y / one_z
    out = {
        # the hard truth
        "z": z, "s_hard": s * one_z, "y_hard": y_hard, "x_hard": x,
        "q2_hard": q2 * np.ones_like(one_z),
        # electron method with the nominal beam
        "q2_e": q2 / one_z, "y_e": y + z,
        "x_e": q2 / (one_z * s * (y + z)),
        # Sigma method (no beam energy anywhere)
        "q2_sigma": q2 * np.ones_like(one_z), "y_sigma": y_hard,
        "x_sigma": q2 / (s * y_hard),
        # Jacquet-Blondel (y uses the nominal 2 E_e)
        "y_jb": y * np.ones_like(one_z),
        "q2_jb": q2 * (1.0 - y_hard) / (1.0 - y),
        "x_jb": q2 * (1.0 - y_hard) / ((1.0 - y) * s * y),
        # double angle (angles only for y; 4 E_e^2 for Q2)
        "y_da": y_hard, "q2_da": q2 / (one_z * one_z),
        "x_da": q2 / (one_z * one_z * s * y_hard),
        # the chain's mixed method: Q2 from the electron, y from Sigma
        "x_mixed": q2 / (one_z * s * y_hard),
        "x_eda": q2 / (one_z * s * y_hard),
        "x_ejb": q2 / (one_z * s * y),
    }
    return out


def method_bias_table(y, z, s=None):
    """Rows (label, Q2_obs/Q2_hard, y_obs/y_hard, x_obs/x_hard) for the
    five reconstruction methods at one (y, z) -- the method comparison of
    plans/08 D3.  `s` cancels in every ratio and is only carried so the
    caller can pass its own.

    `y` is the NOMINAL y = Q2/(x s) of the bin, i.e. the y BEFORE the
    radiative reading; the hard y is y/(1 - z) (`observed_kinematics`).
    Passing a hard y instead divides by (1 - z) twice.  The rows are
    strong functions of y -- the electron method's y ratio is
    (y + z)(1 - z)/y, which is 1.34 at y = 0.19 but 4.2 at y = 0.026 and
    9.2 at y = 0.010 for the same z = 0.092 -- so a caller that prints
    the table must print the y with it."""
    s = 1.0 if s is None else float(s)
    x = 1.0
    q2 = x * y * s
    k = observed_kinematics(x, q2, y, z, s)
    rows = []
    for lab, q2k, yk, xk in (("electron", "q2_e", "y_e", "x_e"),
                             ("Sigma", "q2_sigma", "y_sigma", "x_sigma"),
                             ("Jacquet-Blondel", "q2_jb", "y_jb", "x_jb"),
                             ("double angle", "q2_da", "y_da", "x_da"),
                             ("mixed (Q2_e, y_Sigma)", "q2_e", "y_sigma",
                              "x_mixed")):
        rows.append((lab, float(k[q2k] / k["q2_hard"]),
                     float(k[yk] / k["y_hard"]),
                     float(k[xk] / k["x_hard"])))
    return rows


# --- the E - p_z handle -------------------------------------------------------

def empz_fraction(resp):
    """(E - p_z)_visible / 2E_e of every pseudo-event of a RecoResponse.

    The incoming electron carries E - p_z = 2E_e; a collinear photon takes
    2 z E_e down the beam pipe with it, so the visible sum of the event is
    2(1 - z)E_e and this returns 1 - z.  It is MEASURABLE: the Sigma method
    already gives the hadronic part through
    Sigma_h + E'(1 - cos theta) = E'(1 - cos theta)/(1 - y_Sigma), so no
    quantity beyond the ones the chain reconstructs is needed.

    This is the cut every HERA analysis used against radiative events
    (H1 and ZEUS: 35-38 < E - p_z < 65-70 GeV against 2E_e = 55 GeV), and
    the chain does NOT apply it -- which is why the bound below is quoted
    both with and without it.  Two caveats when it is applied for real:
    the resolution here is that of y (with the Gaussian stand-in,
    delta(1-z)/(1-z) = y dy/y /(1-y) -- 0.8% at y = 0.03 and 25%), and
    with a real hadronic final state the sum is low by whatever escapes
    the acceptance (Report 2 Table 3: 15-30% of Sigma_h), so the cut has
    to sit after the hadronic-scale calibration and be correspondingly
    loose."""
    if getattr(resp.model, "y_method", None) != "mixed":
        # with y_method = "electron" the chain's own 1 - y_e IS
        # E'(1 - cos theta)/2E_e, so the ratio below is identically 1 and
        # the cut would silently keep everything.  E - p_z is only a
        # handle when y carries an INDEPENDENT hadronic measurement.
        raise ValueError("empz_fraction needs model.y_method == 'mixed' "
                         "(an independent hadronic y); got %r"
                         % (getattr(resp.model, "y_method", None),))
    e_p = np.asarray(resp.e_prime_reco, dtype=float)
    ct = np.cos(np.asarray(resp.theta_reco, dtype=float))
    y_r = np.asarray(resp.y_reco, dtype=float)
    two_ee = 2.0 * resp.sampler.config.electron_energy
    with np.errstate(divide="ignore", invalid="ignore"):
        return e_p * (1.0 - ct) / np.maximum(1.0 - y_r, 1e-12) / two_ee


def empz_keep_mask(resp, lo=0.85, hi=1.15):
    """The boolean the window keeps, WITHOUT touching `resp.eff`.

    Separated from `apply_empz_cut` because the per-bin retention has to
    be taken against the pre-cut selection: `RecoResponse.mask_reco`
    carries an `eff > 0` factor, so a bin mask evaluated after the cut
    already has the rejected events removed and every retention would
    come back 1."""
    f = empz_fraction(resp)
    return np.isfinite(f) & (f >= lo) & (f <= hi)


def empz_bin_retention(resp, edges, lo=0.85, hi=1.15, keep=None):
    """Per-bin (w eff mask keep).sum()/(w eff mask).sum() at the pre-cut
    `resp.eff`, one entry per (xlo, xhi, q2lo, q2hi) of `edges`.

    The global retention `apply_empz_cut` returns is dominated by the
    high-y bulk of the selected sample, where the radiative tail lives;
    what an analysis bin actually pays is this, and at the sweet spots
    (y = 0.01-0.03) it is a different number by two orders."""
    keep = empz_keep_mask(resp, lo, hi) if keep is None else keep
    we = resp.w * resp.eff
    out = []
    for e in edges:
        mr = resp.mask_reco(*e)
        den = float(we[mr].sum())
        out.append(float(we[mr & keep].sum()) / den if den > 0 else np.nan)
    return out


EMPZ_Y_BANDS = (0.0, 0.05, 0.2, 0.5, 1.0)


def empz_y_retention(resp, bands=EMPZ_Y_BANDS, lo=0.85, hi=1.15, keep=None):
    """Retention and rate share in bands of the NOMINAL y, pre-cut.

    Where the window's loss sits is the whole question of whether it can
    be adopted as a default: the reconstructed E - p_z is
    E'(1 - cos theta)/(1 - y_Sigma), so its resolution carries 1/(1 - y)
    and the window rejects unradiated events at high y purely on the
    hadronic y resolution -- 25% for the Gaussian stand-in, a few per
    cent for the calibrated Sigma of the PYTHIA chain.  The analysis
    bins of the tensor programme live at y = 0.01-0.03, which is why the
    global retention badly overstates what they pay.

    `share` is normalised to the BANDED rate and therefore sums to one by
    construction; `covered` is the banded rate as a fraction of the whole
    selected rate, so that the identity the published numbers rest on,

        global retention = covered * sum(keep_b * share_b),

    stays exact whatever `bands` is asked for.  The default bands span
    y = 0 to 1 and `covered` is 1 for any physical selection; a narrower
    `bands` argument -- or, one day, a selection leaking outside [0, 1) --
    shows up there instead of silently rescaling the shares."""
    keep = empz_keep_mask(resp, lo, hi) if keep is None else keep
    y = np.asarray(resp.y_nominal, dtype=float)
    we = resp.w * resp.eff
    tot = float(we.sum())
    out_k, out_d = [], []
    for a, b in zip(bands[:-1], bands[1:]):
        m = (y >= a) & (y < b)
        den = float(we[m].sum())
        out_k.append(float(we[m & keep].sum()) / den if den > 0 else np.nan)
        out_d.append(den)
    banded = float(np.sum(out_d))
    out_s = [d / banded if banded > 0 else np.nan for d in out_d]
    return {"bands": tuple(bands), "keep": np.array(out_k),
            "share": np.array(out_s),
            "covered": banded / tot if tot > 0 else np.nan}


def apply_empz_cut(resp, lo=0.85, hi=1.15):
    """Zero `resp.eff` outside lo <= (E - p_z)/2E_e <= hi, IN PLACE.

    Every reco-level quantity of the response (`mask_reco`, `bin_summary`,
    `fold_kernel`, `expected_counts`) reads `eff` at call time, so this is
    the whole cut.  Returns the accepted fraction of the selected rate --
    the GLOBAL one, which the high-y bulk dominates and which is NOT what
    an analysis bin pays; for that see `empz_bin_retention` and
    `empz_y_retention`, both of which must be called before this."""
    keep = empz_keep_mask(resp, lo, hi)
    before = float((resp.w * resp.eff).sum())
    resp.eff = resp.eff * keep
    after = float((resp.w * resp.eff).sum())
    return after / before if before > 0 else np.nan


# --- the migration bound ------------------------------------------------------

def migration_bound(build_response, edges, category, isr=None,
                    labels=None, empz_cut=None):
    """Purity, efficiency, K and the reco-bin amplitude with ISR on and
    off, with COMMON RANDOM NUMBERS.

    `build_response(isr)` must return a fresh `recopseudo.RecoResponse`
    built from the SAME seed for both calls -- the ISR model draws z from
    its own stream, so the two responses differ by the radiation alone.
    `edges` is a list of (xlo, xhi, q2lo, q2hi).

    The quantity the gate applies to is `d_amp`: the analysis fits A on
    radiative data and converts it with the bin-centering factor
    K = Delta(x_c, Q2_c)/A_reco-bin of an ISR-FREE MC, so
    Delta_hat/Delta = A_reco(ISR)/A_reco(no ISR) -- the "wrong-MC" logic
    of the detector systematics scan, and the same convention as
    `money_cos2phi_reco.py --syst-scan`."""
    isr = (isr or ISRModel()).reset()   # reproducible for any caller,
    off = build_response(None)         # even a reused model instance
    on = build_response(isr)
    keep = bin_keep = None
    if empz_cut is not None:
        # the same cut on both members of the pair, so the comparison
        # stays a comparison of the radiation and not of the selection.
        # The per-bin retention is taken BEFORE the cut is applied, for
        # the reason in `empz_bin_retention`.
        keep, bin_keep, y_keep = [], [], []
        for resp in (off, on):
            km = empz_keep_mask(resp, *empz_cut)
            bin_keep.append(empz_bin_retention(resp, edges, keep=km))
            y_keep.append(empz_y_retention(resp, keep=km))
            keep.append(apply_empz_cut(resp, *empz_cut))
        keep = tuple(keep)
    rows = []
    for i, e in enumerate(edges):
        a = off.bin_summary(*e, category)
        b = on.bin_summary(*e, category)
        rows.append({
            "label": (labels[i] if labels else "bin %d" % (i + 1)),
            "edges": tuple(e),
            "purity_off": a["purity"], "purity_on": b["purity"],
            "eff_off": a["efficiency"], "eff_on": b["efficiency"],
            "amp_off": a["a_true_bin"], "amp_on": b["a_true_bin"],
            "a_reco_off": a["a_reco_bin"], "a_reco_on": b["a_reco_bin"],
            "d_purity": b["purity"] - a["purity"],
            "d_eff": b["efficiency"] - a["efficiency"],
            "d_amp": b["a_reco_bin"] / a["a_reco_bin"] - 1.0,
            "d_k": a["a_reco_bin"] / b["a_reco_bin"] - 1.0,
            "sigma_reco_ratio": b["sigma_reco_pb"] / a["sigma_reco_pb"],
        })
        if bin_keep is not None:
            rows[-1]["empz_keep_off"] = bin_keep[0][i]
            rows[-1]["empz_keep_on"] = bin_keep[1][i]
    return {"rows": rows, "off": off, "on": on, "isr": isr,
            "empz_cut": empz_cut, "empz_keep": keep,
            "empz_y": (tuple(y_keep) if empz_cut is not None else None)}


def migration_bound_seeds(build_pair, edges, category, seeds,
                          isr_seed=20260828, labels=None, empz_cut=None,
                          **isr_kw):
    """`migration_bound` repeated over several RESPONSE seeds.

    The single-seed bound carries a Monte-Carlo scatter of 5e-4 to 9e-4
    in absolute d_amp at the statistics the money script runs -- 4 to 14
    per cent of the bound, and up to two standard deviations for an
    individual seed; a published number therefore has to be a seed
    average with its error, not one draw.  `build_pair(seed, isr_model)` must return a fresh
    `RecoResponse` at that response seed.  Every seed gets its OWN fresh
    `ISRModel(seed=isr_seed)`, so the pair at each seed still sits on
    common random numbers and the seeds differ only in the pseudo-events.

    The returned dict has the same shape as `migration_bound`'s, with the
    per-bin entries replaced by the mean over seeds and `<key>_sem` added
    for `d_amp`, `d_purity`, `d_eff` and `sigma_reco_ratio`; `off`, `on`
    and `isr` are those of the FIRST seed, so a caller can still read the
    z spectrum off the pair.  With `empz_cut` the retentions are averaged
    too -- the global pair as `empz_keep` with `empz_keep_sem`, the
    per-bin ones as row keys, the y bands as `empz_y` -- because they are
    published numbers and a single draw of them is not."""
    seeds = [int(v) for v in seeds]
    runs = [migration_bound(lambda m, sd=sd: build_pair(sd, m), edges,
                            category, isr=ISRModel(seed=isr_seed, **isr_kw),
                            labels=labels, empz_cut=empz_cut)
            for sd in seeds]
    keys = ("purity_off", "purity_on", "eff_off", "eff_on", "amp_off",
            "amp_on", "a_reco_off", "a_reco_on", "d_purity", "d_eff",
            "d_amp", "d_k", "sigma_reco_ratio")
    if empz_cut is not None:
        keys = keys + ("empz_keep_off", "empz_keep_on")
    n = float(len(seeds))
    rows = []
    for i in range(len(edges)):
        r = dict(runs[0]["rows"][i])
        for k in keys:
            v = np.array([run["rows"][i][k] for run in runs], dtype=float)
            r[k] = float(v.mean())
            r[k + "_sem"] = (float(v.std(ddof=1) / np.sqrt(n)) if n > 1
                             else 0.0)
            r[k + "_vals"] = v
        rows.append(r)
    out = dict(runs[0])
    out["rows"] = rows
    out["seeds"] = seeds
    out["runs"] = runs
    if empz_cut is not None:
        # the GLOBAL retention pair is a per-seed number as well
        k = np.array([run["empz_keep"] for run in runs], dtype=float)
        out["empz_keep"] = tuple(k.mean(axis=0))
        out["empz_keep_sem"] = tuple(k.std(axis=0, ddof=1) / np.sqrt(n)
                                     if n > 1 else np.zeros(k.shape[1]))
        out["empz_y"] = tuple(
            {"bands": runs[0]["empz_y"][i]["bands"],
             "keep": np.mean([run["empz_y"][i]["keep"] for run in runs],
                             axis=0),
             "share": np.mean([run["empz_y"][i]["share"] for run in runs],
                              axis=0),
             "covered": float(np.mean([run["empz_y"][i]["covered"]
                                       for run in runs]))}
            for i in (0, 1))
    return out


def format_bound(bound, gate=0.05):
    """The printed table of `migration_bound` (or of
    `migration_bound_seeds`, whose extra +- sem columns it picks up)
    plus the WP4 gate verdict."""
    rows = bound["rows"]
    multi = "seeds" in bound
    out = [bound["isr"].describe()]
    if multi:
        out.append("  mean +- sem over %d response seeds: %s"
                   % (len(bound["seeds"]),
                      ", ".join(str(v) for v in bound["seeds"])))
    empz = bool(bound.get("empz_cut"))
    if empz:
        sem = bound.get("empz_keep_sem")
        fmt = ("%.4f +- %.4f (no ISR) / %.4f +- %.4f (ISR)"
               % (bound["empz_keep"][0], sem[0], bound["empz_keep"][1],
                  sem[1]) if sem else
               "%.4f (no ISR) / %.4f (ISR)" % (bound["empz_keep"][0],
                                               bound["empz_keep"][1]))
        out.append("  E - p_z cut %.2f-%.2f x 2E_e: keeps %s of the "
                   "selected rate GLOBALLY; the per-bin retentions below "
                   "are what the analysis bins pay"
                   % (bound["empz_cut"][0], bound["empz_cut"][1], fmt))
        yk = bound.get("empz_y")
        if yk:
            bands = " ".join("%g-%g" % (a, b) for a, b in
                             zip(yk[0]["bands"][:-1], yk[0]["bands"][1:]))
            out.append("    retention in bands of nominal y (%s): "
                       "no ISR %s / ISR %s; band share of the selected "
                       "rate %s"
                       % (bands,
                          " ".join("%.4f" % v for v in yk[0]["keep"]),
                          " ".join("%.4f" % v for v in yk[1]["keep"]),
                          " ".join("%.3f" % v for v in yk[0]["share"])))
    out.append("  %-22s %8s %8s %8s %8s %9s%s %16s"
               % ("bin", "P(off)", "P(on)", "eff off", "eff on",
                  "sig on/off",
                  " %8s %8s" % ("keep off", "keep on") if empz else "",
                  "dDelta [%]"))
    for r in rows:
        tail = ("%+9.3f +- %.3f" % (100.0 * r["d_amp"],
                                    100.0 * r["d_amp_sem"]) if multi
                else "%+16.3f" % (100.0 * r["d_amp"],))
        mid = (" %8.4f %8.4f" % (r["empz_keep_off"], r["empz_keep_on"])
               if empz else "")
        out.append("  %-22s %8.3f %8.3f %8.3f %8.3f %9.4f%s %s"
                   % (r["label"], r["purity_off"], r["purity_on"],
                      r["eff_off"], r["eff_on"], r["sigma_reco_ratio"],
                      mid, tail))
    if multi:
        out.append("  seed-to-seed spread of dDelta [%]: "
                   + ", ".join("%s %.3f (sd %.3f, min %.3f, max %.3f)"
                               % (r["label"], 100.0 * r["d_amp"],
                                  100.0 * float(np.std(r["d_amp_vals"],
                                                       ddof=1)),
                                  100.0 * float(r["d_amp_vals"].min()),
                                  100.0 * float(r["d_amp_vals"].max()))
                               for r in rows))
    worst_row = max(rows, key=lambda r: abs(r["d_amp"]))
    worst = abs(worst_row["d_amp"]) + (2.0 * worst_row["d_amp_sem"] if multi
                                       else 0.0)
    out.append("  gate (plans/07 WP4): max |dDelta/Delta|%s = %.2f%% vs %.0f%%"
               " -> %s" % (" + 2 sem" if multi else "", 100.0 * worst,
                           100.0 * gate, "PASS" if worst <= gate else "FAIL"))
    return "\n".join(out)
