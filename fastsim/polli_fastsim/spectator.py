"""Cluster-spectator kinematics for e+6Li and e+7Li tagging.

Model: the lithium ground state is a two-cluster system (6Li = alpha + d,
S-wave dominant; 7Li = alpha + t, P-wave, 3/2-). DIS strikes one cluster;
the other is a spectator carrying the internal relative momentum k. We
sample k from simple analytic momentum densities, boost to the lab, and
report the far-forward routing variables (pT, theta, rigidity ratio R).

Momentum-density models (two-parameter, crude BY DESIGN — the high-k tail
is the dominant model uncertainty and exactly what VMC densities/BeAGLE
must eventually pin down; see plans/02 step 1.5.3):

  S-wave (6Li):  psi(k) ∝ 1/(k^2+kappa^2) - 1/(k^2+beta^2)   (Hulthen form)
  P-wave (7Li):  psi(k) ∝ k / [(k^2+kappa^2)(k^2+beta^2)]    (vanishes at 0)

kappa = sqrt(2*mu*S_alpha) from the verified separation energies, with
mu the reduced mass of the two SEPARATED clusters (see m_partner):
  6Li -> alpha+d: S = 1.474 MeV  => kappa = 60.66 MeV
  7Li -> alpha+t: S = 2.467 MeV  => kappa = 88.90 MeV
beta is the short-range scale (default 0.30 GeV, scan 0.20-0.40 to span
the tail uncertainty; deuteron analogy: beta/kappa ~ 5.7).

Masses are the physical nuclear ones (NUCLEUS_MASS), not A * M_U.  The
visible consequence is the rigidity: a spectator at rest in the beam
frame has p_lab = (p_beam / m_beam) * m_spec, so
    R(k=0) = (m_spec / Z_spec) / (m_beam / Z_beam),
a ratio of mass-to-charge ratios that A * M_U collapses onto the naive
(A_spec Z_beam) / (A_beam Z_spec).  The alpha is more bound per nucleon
than the lithium it comes from, so it sits BELOW beam rigidity:
R = 0.99813 (not 1.00000) from 6Li and 0.85571 (not 0.85714) from 7Li
-- 0.19% and 0.17%, entirely inside the +-5% near-beam band, but enough
to slide the 6Li alpha across the RP window edge at R = 0.95 (plans/08
C1 measurement, recorded in spectator_lab_kinematics).

Upgrade path: replace n(k) with VMC two-cluster overlap densities
(R.B. Wiringa et al., tables available from ANL) — same interface.
"""

from dataclasses import dataclass

import numpy as np

M_U = 0.93149  # amu [GeV]
MASSES = {"p": 0.93827, "n": 0.93957, "d": 1.87561, "t": 2.80892,
          "3He": 2.80839, "alpha": 3.72738}

# Ground-state NUCLEAR (not atomic) masses [GeV], keyed by (Z, A) so a
# ClusterChannel can look up its own beam.  A * M_U is NOT the mass of a
# nucleus: it drops the mass excess and keeps the electrons, so it is
# low by (excess - Z m_e) = 12.6 MeV for d, 14.5 for t, 13.9 for 3He,
# 1.4 for the alpha, 12.6 for 6Li and 13.4 for 7Li.  The beam mass sets
# the boost gamma and gamma*beta, and with them the lab momentum, theta
# and the rigidity ratio R of every spectator, so all of those were
# wrong by 2.0-2.2e-3 for the lithium beams and 6.7e-3 for the deuteron
# control channel.
#
# Built from the AME2020 ATOMIC masses (Wang, Huang, Kondev, Audi, Naimi,
# Chin. Phys. C 45 (2021) 030003, Table I) [u]
#     1H  1.00782503190   2H  2.01410177784   3H   3.01604928132
#     n   1.00866491590   3He 3.01602932197   4He  4.00260325413
#     6Li 6.01512288742   7Li 7.01600343426
# by removing the electrons,
#     M_nuc = M_atomic * u - Z * m_e + B_e(Z),
# with u = 931.49410242 MeV and m_e = 0.51099895000 MeV (CODATA 2018,
# Tiesinga, Mohr, Newell, Taylor, Rev. Mod. Phys. 93 (2021) 025010,
# Table I) and the total atomic electron binding energies B_e =
# 13.598 / 79.005 / 203.486 eV for Z = 1 / 2 / 3 (sums of the NIST ASD
# ionization energies; 2e-7 GeV, at the edge of the digits kept here).
# Before rounding, the recipe reproduces the six CODATA nuclear masses
# m_p, m_n, m_d, m_t, m_h, m_alpha to <= 1.9e-10 GeV (worst the
# deuteron).  The entries below are that recipe kept to 1 eV, so as
# WRITTEN they reproduce those six to <= 5.7e-10 GeV -- the rounding,
# not the evaluation.  test_spectator.py pins the table as written, at
# 2e-9 GeV: an external check of these digits, not of this arithmetic.
NUCLEUS_MASS = {
    (0, 1): 0.939565420,   # n
    (1, 1): 0.938272088,   # p
    (1, 2): 1.875612942,   # d
    (1, 3): 2.808921133,   # t
    (2, 3): 2.808391607,   # 3He
    (2, 4): 3.727379407,   # alpha (4He)
    (3, 6): 5.601518702,   # 6Li
    (3, 7): 6.533833028,   # 7Li
}

# MASSES is the same evaluation rounded to 10 keV: every entry agrees
# with NUCLEUS_MASS to <= 4.6 keV (worst the neutron, 4.9e-6 relative;
# the alpha is 0.6 keV, 1.6e-7), i.e. to its own last digit.  It is kept
# verbatim rather than re-derived from NUCLEUS_MASS because it sets the
# SPECTATOR mass in every published distribution and 5e-6 relative is
# far below anything the cluster model resolves.


def nucleus_mass(z, a):
    """Ground-state nuclear mass [GeV] of the nuclide (Z, A).

    Falls back to A * M_U for a nuclide the table does not carry.  That
    is the mass with the mass excess dropped and the Z electrons left
    in, so for everything in the table it is LOW -- by 0.04% for the
    alpha, 0.2% for 6,7Li, 0.7% for the deuteron -- because each of
    those is less bound per nucleon than the 12C that defines u by more
    than its electrons weigh.  The sign is not universal: u is one
    twelfth of the mass of a neutral 12C ATOM, so at (Z, A) = (6, 12)
    the fallback returns the atomic mass and is HIGH of the nuclear one
    by 6 m_e - B_e(6) = 3.016 MeV (test_spectator.py).  Extend
    NUCLEUS_MASS rather than lean on it.
    """
    return NUCLEUS_MASS.get((int(z), int(a)), a * M_U)


@dataclass(frozen=True)
class ClusterChannel:
    """Beam nucleus -> struck cluster + spectator cluster."""
    name: str            # e.g. "6Li -> d* + alpha(spec)"
    beam_A: int
    beam_Z: int
    spectator: str       # key into MASSES
    spectator_A: int
    spectator_Z: int
    separation_energy: float  # [GeV]
    l_wave: int          # 0 (S) or 1 (P) relative motion

    @property
    def m_spec(self):
        return MASSES[self.spectator]

    @property
    def m_beam(self):
        """Ground-state nuclear mass of the beam [GeV] (NUCLEUS_MASS)."""
        return nucleus_mass(self.beam_Z, self.beam_A)

    @property
    def partner_A(self):
        """Mass number of the struck (partner) cluster."""
        return self.beam_A - self.spectator_A

    @property
    def partner_Z(self):
        """Charge of the struck (partner) cluster."""
        return self.beam_Z - self.spectator_Z

    @property
    def m_partner(self):
        """Mass of the FREE struck cluster [GeV].

        kappa is a two-body bound-state scale, so the reduced mass it
        needs is the one of the two SEPARATED fragments.  Inside the
        nucleus the partner is worth only m_beam - m_spec, which is
        m_partner - S: the separation energy is exactly the difference.
        Either the table (used here) or m_beam - m_spec + S is the free
        mass; m_beam - m_spec alone is low by S, and that -- not the
        beam mass -- is what the pre-2026-08-26 kappa used.  The table
        is preferred because it does not propagate the rounding of
        `separation_energy`; the two routes agree to <= 3.6 keV, which
        test_spectator.py pins.  Every partner of every channel defined
        below is in the table, so the fallback is reached only by a
        channel a caller builds; it is written to keep the DEFINITION
        of the separation energy, S = m_spec + m_partner - m_beam,
        exact off the table too, and test_spectator.py exercises it on
        a 12C -> alpha + 8Be channel because nothing here does.
        """
        z = self.beam_Z - self.spectator_Z
        a = self.beam_A - self.spectator_A
        if (z, a) in NUCLEUS_MASS:
            return NUCLEUS_MASS[(z, a)]
        return self.m_beam - self.m_spec + self.separation_energy

    @property
    def kappa(self):
        """Bound-state momentum scale sqrt(2 mu S) [GeV].

        mu is the reduced mass of the two FREE clusters: kappa is the
        asymptotic decay constant of the relative wave function
        psi ~ exp(-kappa r) / r of two fragments separated by S, so the
        Schroedinger problem it comes from is the one of the separated
        fragments (m_spec, m_partner).  Before 2026-08-26 the beam mass
        was beam_A * M_U and the partner was m_beam - m_spec; those two
        errors partly cancelled in mu, leaving kappa 0.06-0.25% low and
        -- the symptom -- DIFFERENT for the two channels of the same
        beam (60.50 vs 60.62 MeV for 6Li), which the two-body scale
        cannot be (plans/08 C1).
        """
        mu = self.m_spec * self.m_partner / (self.m_spec + self.m_partner)
        return np.sqrt(2.0 * mu * self.separation_energy)


DEUTERON_P_TAG = ClusterChannel(
    "d: DIS on neutron, proton spectator (control)", 2, 1, "p", 1, 1,
    2.2246e-3, l_wave=0)
DEUTERON_N_TAG = ClusterChannel(
    "d: DIS on proton, neutron spectator (control)", 2, 1, "n", 1, 0,
    2.2246e-3, l_wave=0)
HE3_P_TAG = ClusterChannel(
    "3He: p spectator (p+d two-body approx., control)", 3, 2, "p", 1, 1,
    5.49e-3, l_wave=0)
LI6_ALPHA_TAG = ClusterChannel(
    "6Li: DIS on d-cluster, alpha spectator", 6, 3, "alpha", 4, 2,
    1.4743e-3, l_wave=0)
LI6_D_TAG = ClusterChannel(
    "6Li: DIS on alpha-cluster, d spectator", 6, 3, "d", 2, 1,
    1.4743e-3, l_wave=0)
LI7_ALPHA_TAG = ClusterChannel(
    "7Li: DIS on t-cluster, alpha spectator", 7, 3, "alpha", 4, 2,
    2.4670e-3, l_wave=1)
LI7_T_TAG = ClusterChannel(
    "7Li: DIS on alpha-cluster, triton spectator", 7, 3, "t", 3, 1,
    2.4670e-3, l_wave=1)

CHANNELS = (LI6_ALPHA_TAG, LI6_D_TAG, LI7_ALPHA_TAG, LI7_T_TAG)


def momentum_density(k, kappa, beta, l_wave):
    """Unnormalized n(k) (NOT including the k^2 phase-space factor)."""
    k = np.asarray(k, dtype=float)
    if l_wave == 0:
        psi = 1.0 / (k * k + kappa * kappa) - 1.0 / (k * k + beta * beta)
    else:
        psi = k / ((k * k + kappa * kappa) * (k * k + beta * beta))
    return psi * psi


def sample_k(channel, n, beta=0.30, k_max=1.5, rng=None):
    """Sample |k| from k^2 * n(k) on [0, k_max] GeV by inverse-CDF on a grid,
    plus isotropic angles. Returns (kx, ky, kz) arrays."""
    rng = rng or np.random.default_rng(20260612)
    grid = np.linspace(1e-4, k_max, 30000)
    pdf = grid * grid * momentum_density(grid, channel.kappa, beta,
                                         channel.l_wave)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    k = np.interp(rng.uniform(size=n), cdf, grid)
    cos_t = rng.uniform(-1.0, 1.0, size=n)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n)
    sin_t = np.sqrt(1.0 - cos_t**2)
    return k * sin_t * np.cos(phi), k * sin_t * np.sin(phi), k * cos_t


def _boost_fragment(channel, p_per_nucleon, kx, ky, kz, m, frag_Z, frag_A):
    """Lab kinematics of ONE fragment of a breakup, from its rest-frame
    momentum (kx, ky, kz), mass `m`, charge `frag_Z` and mass number
    `frag_A`.  Returns the dict `spectator_lab_kinematics` documents.

    The single implementation of the boost, the rigidity ratio and the
    azimuth convention, shared by `spectator_lab_kinematics` (the spectator
    alone, +k) and `breakup_lab_kinematics` (both fragments, +-k) so that
    the two cannot drift apart -- the lesson of plans/08 C1, where one
    two-body scale computed in two places came out as two numbers.
    """
    e_rest = np.sqrt(m * m + kx * kx + ky * ky + kz * kz)
    # beam boost (per-nucleon momentum sets the velocity of the nucleus).
    # The beam mass is the physical nuclear mass: with beam_A * M_U the
    # boost gamma*beta = p_beam / m_beam, and with it p_lab, theta and R,
    # is high by the relative mass excess (2.2e-3 for 6Li, 2.1e-3 for
    # 7Li) -- see NUCLEUS_MASS.
    m_beam = channel.m_beam
    p_beam = channel.beam_A * p_per_nucleon
    e_beam = np.sqrt(p_beam**2 + m_beam**2)
    gamma = e_beam / m_beam
    gbeta = p_beam / m_beam
    pz_lab = gamma * kz + gbeta * e_rest
    pt = np.sqrt(kx * kx + ky * ky)
    p_lab = np.sqrt(pt * pt + pz_lab * pz_lab)
    theta = np.arctan2(pt, pz_lab)
    rigidity_beam = p_beam / channel.beam_Z
    if frag_Z > 0:
        rig_ratio = (p_lab / frag_Z) / rigidity_beam
    else:
        rig_ratio = np.full_like(p_lab, np.nan)  # neutral: no rigidity
    return {
        "pT": pt,
        "theta": theta,
        "phi": np.arctan2(ky, kx),        # lab azimuth, for a planar cut
        "p_lab": p_lab,
        "R": rig_ratio,
        "xL": p_lab / (frag_A * p_per_nucleon),
        "k": np.sqrt(kx**2 + ky**2 + kz**2),
    }


def spectator_lab_kinematics(channel, p_per_nucleon, n=200_000, beta=0.30,
                             rng=None):
    """Boost spectator (rest-frame momentum -k of the struck cluster) to the
    lab. Returns dict of arrays: pT, theta [rad], p_lab, R (rigidity ratio
    vs beam), xL (= p_lab / (A_spec * p_per_nucleon)).

    What the physical beam mass buys (measured 2026-08-26 against
    beam_A * M_U, 2e6 events, beta = 0.30, seed 7, plans/08 C1).  It is
    a 2e-3 effect on the kinematics -- R low by 2.0-2.3e-3, theta high
    by 2.2-4.2e-3, pT and |k| high by 1.3e-4 to 2.0e-3 (those two come
    from kappa, which the same masses fix) at every quantile from 5% to
    99% and at every documented energy (measured at the pre-2026-08-27
    6Li 20.5/50/137.5 and 7Li
    17.6/42.9/117.9 GeV/u) -- but NOT a 2e-3 effect on the acceptance,
    because the R window edges are hard.  Sliding the 6Li alpha, which
    sits at R ~ 1, down by 0.0022 moves events across the RP/near-beam
    boundary at R = 0.95, where the density is high: the Roman-Pot
    fraction rises 1.31% -> 1.51% (+15% relative) and the total alpha
    tag 1.65% -> 1.85% (legacy 73 urad "high-acceptance") and 1.31% ->
    1.51% (legacy 164 urad) at 137.5 GeV/u; 13.2% -> 13.5% at 50 GeV/u.
    (These legacy numbers are the AZIMUTH-BLIND circular cut the routing
    applied at the time; since 2026-08-28 the published tag is evaluated
    per configuration at the Yellow Report optics -- 1.7 / 1.5 / 1.6% --
    and at the tagging optics, farforward.yr_optics / tagging_optics,
    with the rectangle applied to each fragment's azimuth.)  The 7Li
    alpha, off rigidity at R = 0.856 with no nearby edge, moves by
    +0.13% relative (97.7% -> 97.9%).  Splitting the two changes shows
    the acceptance shift is the mass alone: at 137.5 GeV/u the new
    kappa moves the 6Li tag 1.645% -> 1.652% and the new beam mass
    1.645% -> 1.839%.
    """
    kx, ky, kz = sample_k(channel, n, beta=beta, rng=rng)
    return _boost_fragment(channel, p_per_nucleon, kx, ky, kz,
                           channel.m_spec, channel.spectator_Z,
                           channel.spectator_A)


def breakup_lab_kinematics(channel, p_per_nucleon, n=200_000, beta=0.30,
                           rng=None):
    """Both fragments of one two-body breakup, correlated.

    `spectator_lab_kinematics` boosts ONE fragment; the two channels of a
    beam (6Li -> alpha + d as LI6_ALPHA_TAG and LI6_D_TAG) sampled
    separately are two unrelated events.  Here a single relative momentum k
    is drawn and BOTH fragments are boosted from it: the spectator with +k
    and the partner with -k, which is all the two-body rest-frame condition
    says.  The two channels of a beam share kappa and carry complementary
    masses, so nothing new is assumed -- only the correlation that sampling
    them apart throws away.

    Returns {"spectator": {...}, "partner": {...}, "k": |k|, "kx", "ky",
    "kz"}, the two fragment dicts carrying exactly the keys
    `spectator_lab_kinematics` returns and the rest-frame components being
    the SPECTATOR's (the partner's are their negatives).

    The observable consequence, and the reason this is worth sampling: the
    two fragments take opposite transverse kicks of equal size but are
    carried by different longitudinal momenta, so their lab angles are in
    the inverse ratio of their masses and their azimuths differ by pi.  For
    6Li -> alpha + d that is theta_d / theta_alpha -> m_alpha / m_d = 1.987
    as k -> 0 (the naive mass-number ratio would say 2), and the deuteron
    is always the wider, easier fragment.
    """
    kx, ky, kz = sample_k(channel, n, beta=beta, rng=rng)
    spec = _boost_fragment(channel, p_per_nucleon, kx, ky, kz,
                           channel.m_spec, channel.spectator_Z,
                           channel.spectator_A)
    part = _boost_fragment(channel, p_per_nucleon, -kx, -ky, -kz,
                           channel.m_partner, channel.partner_Z,
                           channel.partner_A)
    return {"spectator": spec, "partner": part, "kx": kx, "ky": ky, "kz": kz,
            "k": spec["k"]}
