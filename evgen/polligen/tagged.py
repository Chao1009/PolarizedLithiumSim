"""Tagged mode: spin-correlated (e', spectator) events (plans/05 step 5.B).

Model: two-cluster light-front-flavored impulse approximation.  The ion
ground state is expanded in cluster relative partial waves,

  |J M> = sum_L a_L sum_{m_L, m_S} <L m_L S_c m_S | J M>
            psihat_L(k) Y_L^{m_L}(khat) |S_c m_S>,

with S_c the channel spin (the coupled cluster spins), psihat_L normalized
radial waves and a_L^2 = P_L the wave probabilities.  The radial tables
(`Wave.radial`, and the tabulated AV18 control) are in the PLAIN
Bessel-transform convention, psihat_L(k) ~ int j_L(kr) u_L(r) r dr, which
is positive at low k -- so the momentum-space amplitude carries the
plane-wave phase phi_L = i^L psihat_L: a coordinate-space
sum_L u_L(r)/r [Y_L x chi_S]_JM transforms into sum_L (-i)^L psihat_L(k)
[Y_L x chi_S]_JM.  Both phases appear because they belong to opposite
conventions -- (-i)^L is the expansion of e^{-ik.r}, i^L that of
e^{+ik.r} -- and they differ by (-1)^L, a GLOBAL sign within one parity;
only the RELATIVE phase between waves is observable, so nothing here
depends on which of the two is taken.  Everything the tagged
observables need follows from the joint amplitude

  A_{m_S}(M; k, khat) = sum_L i^L a_L psihat_L(k) C_L(M, m_S)
                              Y_L^{M-m_S}(khat):

For the same-parity mixtures this module allows (guarded in
TaggedChannel.__post_init__) the relative phase is real, (-1)^(L//2) --
+1 and -1 for L = 0 and 2, +1 for a lone L = 1 wave -- so A stays real.
Dropping it inverts the S-D interference: it moves the axial node of the
deuteron density from M = +-1 to M = 0 and flips the sign of A_zz^wf
(run 19, 2026-09-15; the pre-fix sign was published through 2026-09-06).

* spectator momentum density  n_M(k, khat) = sum_{m_S} |A_{m_S}|^2
  (L-interference gives the Cosyn-Weiss-style m-dependence: S/D for the
  deuteron and the alpha-d system, |Y_1^{m_L}|^2 for alpha-t);
* struck-cluster spin populations p(m_S | M, k, khat) = |A_{m_S}|^2 / n_M
  (diagonal truncation: coherences feed phi-dependent structures beyond
  the Step-5.A master formula and are dropped, documented);
* pair decomposition p(m_struck, m_spec | m_S) via a second CG factor
  (e.g. the struck-neutron polarization inside the deuteron's triplet).

The tagged cross section then factorizes (IA, no FSI -- quote at low
spectator virtuality; plans/04 #16):

  dsigma^tag(lam_e, M) = |A_{m_S}(M; k)|^2 (x) dsigma_struck(lam_e, m_S)

with dsigma_struck the Step-5.A inclusive master formula evaluated on the
struck cluster's own structure functions (embedded deuteron: F_{1,2}d,
b1_d, g1_d -- this is the embedded-b1 observable; quasi-free triton: g1_t
with (P_p, P_n) = (0.86, -0.028)).  DIS kinematics use per-nucleon x at
the ion's per-nucleon beam momentum; the light-cone alpha_s dependence of
the struck-cluster SFs is beyond this tier (Cosyn-Weiss upgrade path).

Angular conventions: khat is measured in the SPIN frame (quantization
axis = z).  All densities are even in cos(theta_k) (no parity-odd L
mixing), so the spectator-vs-struck momentum sign convention drops out.
For longitudinal fills the spin frame coincides with the lab frame up to
the (uniform) azimuth; tilted axes are rotated before the lab boost.

Verified consequences implemented here (tested in test_tagged.py):
* 7Li alpha-tag P-wave: n_{3/2} ~ sin^2(theta_k), <P2(cos theta_k)> =
  -T/5 for fill tensor moment T -- the in-situ alignment polarimeter;
  triton polarization P_t(M=3/2)=1, P_t(M=1/2)=1/3, and the khat-resolved
  P_t(theta; M=1/2) = (5c^2-1)/(3c^2+1);
* 6Li alpha-tag: vector dilution 1 - (3/2) P_D of the embedded deuteron
  (default P_D chosen to reproduce the 0.87 of
  `polarized.b1_li6_from_deuteron`); S-wave limit reduces exactly to the
  inclusive polarized-deuteron master formula;
* deuteron control: the S/D interference tensor structure of n_m(k, theta)
  (the Cosyn-Weiss tagged-Azz mechanism), O(1) at k ~ 300 MeV/c.

RELATION TO `spectator.spectator_lab_kinematics`, which is the OTHER
sampler of the same 6Li alpha tag (fastsim/scripts/tagging_acceptance.py,
Report 3 Table 6).  The two give different spectra and neither is wrong:
`spectator` carries ONE partial wave per channel (`ClusterChannel.l_wave`
= 0 for the alpha-d system), this module carries the full S+D expansion,
and the D wave is the whole difference.  Measured on the current tables
(beta = 0.30, kappa = 60.7 MeV): the S-wave radial here is EXACTLY
`spectator.momentum_density`, <k> = 0.1071 GeV/c in both; the D wave is
hard, <k> = 0.2778 GeV/c, and at P_D = 0.0867 it pulls the channel mean
to <k> = 0.1219 GeV/c and the off-rigidity R < 0.95 slice from 1.5% to
2.5%.  Set p_d = 0 and the two samplers agree quantile by quantile
(test_boost_matches_fastsim_spectator).  The D wave is not an optional
tail: it IS the tensor observable -- with P_D = 0 the alpha-d density is
m-independent and A_zz^tag vanishes identically -- so the tagged
observables must be quoted on the S+D spectrum and the pure-model
acceptance table on the S-wave one, with the 1.5 vs 2.5% gap stated
wherever both appear (plans/09 B2, docs/reproduction_manual SS7).
"""

from dataclasses import dataclass
from typing import Tuple

import numpy as np

_trapezoid = getattr(np, "trapezoid", None) or np.trapz

from polli_fastsim import beams, spectator
from polli_fastsim.farforward import route_charged, yr_config_key, yr_optics

from . import bookkeeping as bk
from .sample import InclusiveSampler
from .spin import clebsch_gordan, m_values

# Struck-cluster DIS targets that are not beam species.  TRITON mirrors
# beams.HE3 under p <-> n, and its constants are PER NUCLEON like every
# other Ion slot since plans/08 D7: with Z = 1 and N = 2, ToyG1.g1_nucleus
# gives g1(t) = 0.86 g1p + 2 x (-0.028) g1n.  The second neutron used to
# be dropped, which moved the tagged_polarimetry_7li g1t/F1t overlay by
# +8.4% at x = 0.005 falling to -0.7% at x = 0.7.  No published number
# reads that overlay, but the alpha-tag acceptances the same script
# tabulates DO move, in the fourth decimal, because the accepted sample is
# cross-section-weighted and the cross section carries g1 of the struck
# triton: Roman-Pot tag 0.9614 / 0.9676 / 0.9726 -> 0.9620 / 0.9678 /
# 0.9728 at the Yellow Report optics (docs/reproduction_manual.md 4.1,
# plans/08 D7).  The figure was regenerated with the change.
TRITON = beams.Ion("t", 3, 1, 0.5, eff_pol_p=0.86, eff_pol_n=-0.028)
NEUTRON = beams.Ion("n", 1, 0, 0.5, eff_pol_p=0.0, eff_pol_n=1.0)


def theta_lm(l, m, c):
    """|m|-azimuth-stripped spherical harmonic Theta_l^m(theta) with
    Condon-Shortley signs: Y_l^m = Theta_l^m(theta) exp(i m phi).

    These are the harmonic's own signs only.  The plane-wave phase i^L of
    the momentum-space amplitude is a separate factor, applied per wave in
    `TaggedModel._amp2_table`; the two never interact, because two waves
    of a channel share a cell only at m_l = 0, where Theta_L^0 is real,
    positive-normalized and free of any Condon-Shortley (-1)^m."""
    c = np.asarray(c, dtype=float)
    s = np.sqrt(np.maximum(1.0 - c * c, 0.0))
    if l == 0 and m == 0:
        return np.full_like(c, np.sqrt(1.0 / (4.0 * np.pi)))
    if l == 1:
        if m == 0:
            return np.sqrt(3.0 / (4.0 * np.pi)) * c
        if abs(m) == 1:
            return -np.sign(m) * np.sqrt(3.0 / (8.0 * np.pi)) * s
    if l == 2:
        if m == 0:
            return np.sqrt(5.0 / (16.0 * np.pi)) * (3.0 * c * c - 1.0)
        if abs(m) == 1:
            return -np.sign(m) * np.sqrt(15.0 / (8.0 * np.pi)) * s * c
        if abs(m) == 2:
            return np.sqrt(15.0 / (32.0 * np.pi)) * s * s
    raise ValueError("theta_lm implemented for l <= 2, got (l=%s, m=%s)"
                     % (l, m))


@dataclass(frozen=True)
class Wave:
    """One cluster relative partial wave: probability and radial shape.

    Radial forms (unnormalized; kappa from the channel's separation
    energy, beta the short-range scale scanned as the model band):
      L=0  Hulthen        1/(k^2+kappa^2) - 1/(k^2+beta^2)
      L=1  P-wave         k / ((k^2+kappa^2)(k^2+beta^2))
      L=2  D-wave         k^2 / ((k^2+kappa^2)(k^2+beta^2)^2)
    """
    l: int
    prob: float
    beta: float = 0.30

    def radial(self, k, kappa):
        k = np.asarray(k, dtype=float)
        k2 = k * k
        if self.l == 0:
            return 1.0 / (k2 + kappa**2) - 1.0 / (k2 + self.beta**2)
        if self.l == 1:
            return k / ((k2 + kappa**2) * (k2 + self.beta**2))
        if self.l == 2:
            return k2 / ((k2 + kappa**2) * (k2 + self.beta**2) ** 2)
        raise ValueError("Wave.l must be 0, 1, or 2")


@dataclass(frozen=True)
class TabulatedWave:
    """One partial wave whose radial shape is a TABLE, not a formula.

    Same duck type as `Wave` -- `l`, `prob`, `radial(k, kappa)` -- but the
    table already is the wave function, so `kappa` is accepted and
    ignored.  `k_table` is in GeV/c and must be increasing; values are
    linearly interpolated and ZERO beyond the table's last point (the AV18
    table runs to 20 fm^-1 = 3.95 GeV/c, far outside any grid used here).
    Signs are the table's own: no phase is applied at load, because the
    i^L phase of the momentum-space amplitude belongs to
    `TaggedModel._amp2_table` and is applied there for every wave alike.
    """
    l: int
    prob: float
    k_table: Tuple[float, ...]
    psi_table: Tuple[float, ...]
    label: str = ""

    def radial(self, k, kappa=None):
        kt = np.asarray(self.k_table, dtype=float)
        pt = np.asarray(self.psi_table, dtype=float)
        return np.interp(np.asarray(k, dtype=float), kt, pt,
                         left=pt[0], right=0.0)


@dataclass(frozen=True)
class TaggedChannel:
    """Spin structure on top of a kinematic spectator.ClusterChannel."""
    base: spectator.ClusterChannel   # masses, separation energy, boost
    j_ion: float
    s_struck: float                  # struck-cluster spin
    s_spec: float                    # spectator-cluster spin
    s_channel: float                 # coupled channel spin S_c
    waves: Tuple[Wave, ...]
    dis_target: beams.Ion            # SF target for the struck cluster
    label: str = ""

    def __post_init__(self):
        tot = sum(w.prob for w in self.waves)
        if abs(tot - 1.0) > 1e-9:
            raise ValueError("wave probabilities must sum to 1 (got %g)"
                             % tot)
        # All waves must share L mod 2.  The momentum-space amplitude
        # carries phi_L = i^L psihat_L (module docstring); for one parity
        # that phase is real up to an irrelevant overall factor, which is
        # what lets `TaggedModel._amp2_table` sum REAL amplitudes with the
        # relative factor (-1)^(L//2).  Mix parities and i^L is imaginary
        # for the odd waves, so the amplitude -- and the machinery here --
        # would have to become complex.  Parity also forbids the mixture
        # physically for a state of good parity, so this is a guard
        # against a mis-specified channel, not a missing feature.
        parities = sorted({w.l % 2 for w in self.waves})
        if len(parities) > 1:
            raise ValueError(
                "all waves of a channel must share L mod 2 (got L = %s): "
                "with mixed parity the i^L phase of the momentum-space "
                "amplitude is imaginary for the odd waves and the "
                "amplitude would have to be complex"
                % ([w.l for w in self.waves],))


# --- default channels (radial betas carry the model band, as in fastsim) --

# The two D-state probabilities are RE-EXPORTED from polli_fastsim.beams,
# which is where they live since 2026-08-29: the inclusive effective
# polarization of 6Li (beams.LI6_CLUSTER_POLARIZATION, the two dilutions
# multiplied) and the tagged S/D interference below are then the same
# wave function seen in two experiments, and cannot drift apart.  The
# names, values and meanings are unchanged -- alpha-d D state chosen so
# the embedded-deuteron vector dilution 1 - (3/2) P_D reproduces the 0.87
# of b1_li6_from_deuteron (SCENARIO -- VMC overlaps are the scheduled
# replacement, plans/04 #15), and the deuteron's own AV18-like D state.
P_D_LI6 = beams.P_D_LI6
P_D_DEUTERON = beams.P_D_DEUTERON

# SIGN of the alpha-d D wave -- a live physics input, not a convention.
# P_D_LI6 is a probability and fixes only the magnitude; since the i^L
# phase was restored (run 19) the SIGN of the alpha-d D radial relative to
# the S wave sets the sign of the 6Li tagged A_zz, and this model takes it
# deuteron-like (psihat_2/psihat_0 > 0, both Hulthen-type forms positive).
# LiPolGen's VMC alpha+d overlap measures THREE sign regions, not one:
# sign(psi_2/psi_0) = -1 below the alpha-d S node at 0.134 GeV/c, +1
# between the nodes, and -1 again above the alpha-d D node at 0.444 GeV/c.
# The Hulthen-type forms used here are node-free and positive-definite and
# can represent NEITHER reversal, so the adopted sign is the VMC's
# measurement over 0.134-0.444 GeV/c only and the model's assumption
# outside it -- and the accepted sample is not confined to that window.
# Acceptance-weighted over this model's density at the three 6Li
# configurations (5x40.8 / 10x99.5 / 18x137.5 GeV/u): 27 / 21 / 26 % of
# the Yellow-Report-accepted alpha lie ABOVE 0.444 GeV/c, where LiPolGen's
# own acceptance-weighted VMC A_zz^tag changes sign between k = 0.40 and
# 0.50 (-0.08 -> +0.12) while its Hulthen channel stays at -0.73 / -0.59;
# and 41 / 28 / 37 % of the tagging-optics sample lies BELOW 0.134 GeV/c.
# The k = 0.325 GeV/c headline bin of money plot 4 is inside the supported
# window; the tails on either side of it are not.  VMC overlaps are the
# scheduled replacement of this whole radial input (plans/04 #15).


def li6_alpha_channel(beta=0.30, p_d=P_D_LI6):
    """6Li: DIS on the embedded deuteron, alpha spectator (S+D waves)."""
    return TaggedChannel(spectator.LI6_ALPHA_TAG, 1.0, 1.0, 0.0, 1.0,
                         (Wave(0, 1.0 - p_d, beta), Wave(2, p_d, beta)),
                         beams.DEUTERON, "6Li alpha-tag (embedded d)")


def li7_alpha_channel(beta=0.30):
    """7Li: DIS on the quasi-free triton, alpha spectator (pure P-wave)."""
    return TaggedChannel(spectator.LI7_ALPHA_TAG, 1.5, 0.5, 0.0, 0.5,
                         (Wave(1, 1.0, beta),), TRITON,
                         "7Li alpha-tag (quasi-free t)")


def deuteron_channel(beta=0.30, p_d=P_D_DEUTERON):
    """Deuteron control: DIS on the neutron, proton spectator (S+D).
    The Cosyn-Weiss tagged limit of the machinery."""
    return TaggedChannel(spectator.DEUTERON_P_TAG, 1.0, 0.5, 0.5, 1.0,
                         (Wave(0, 1.0 - p_d, beta), Wave(2, p_d, beta)),
                         NEUTRON, "d control (n struck, p tagged)")


# --- AV18 deuteron control (the real wave function, not a toy) ---------

# `fastsim/polli_fastsim/data/av18/fdeut.av18`, R. B. Wiringa (ANL), the
# raw served bytes (provenance in `data/SOURCES.md`).  Its k-block is
# `k [fm^-1]  u(k)  w(k)`, 201 rows, k = 0 .. 20 fm^-1, with
# u(k), w(k) = sqrt(2/pi) x the PLAIN Bessel transforms and both positive
# at low k -- exactly the convention `_amp2_table`'s i^L phase assumes.
# The file header's `dstate` is the D-state probability.
HBARC_GEV_FM = 0.197327
P_D_AV18_DEUTERON = 0.057599
_AV18_CACHE = {}


def _av18_deuteron_tables():
    """(k [GeV/c], u(k), w(k)) of the committed AV18 deuteron file.

    Read with `importlib.resources` from `polli_fastsim/data/av18`, the
    same way `polarized._load_curve` reads the digitized CSVs, so the
    table travels with the package however it is put on the path.  Cached.
    """
    if "d" not in _AV18_CACHE:
        try:
            from importlib.resources import files
            text = (files("polli_fastsim") / "data" / "av18"
                    / "fdeut.av18").read_text()
        except Exception:                  # pragma: no cover - fallback
            import os
            here = os.path.join(os.path.dirname(
                os.path.abspath(beams.__file__)), "data", "av18",
                "fdeut.av18")
            with open(here, encoding="utf-8") as fh:
                text = fh.read()
        lines = text.splitlines()
        i = next(j for j, ln in enumerate(lines)
                 if ln.strip().startswith("k ") and "u(k)" in ln)
        rows = []
        for ln in lines[i + 1:]:
            parts = ln.split()
            if len(parts) != 3:
                if rows:
                    break
                continue
            rows.append([float(x) for x in parts])
        tab = np.array(rows)
        if tab.shape != (201, 3):
            raise ValueError("unexpected AV18 k-block shape %s" % (tab.shape,))
        _AV18_CACHE["d"] = (tab[:, 0] * HBARC_GEV_FM, tab[:, 1], tab[:, 2])
    return _AV18_CACHE["d"]


def av18_deuteron_channel():
    """Deuteron control on the REAL AV18 wave function.

    Same kinematics and spin structure as `deuteron_channel`, with the
    Hulthen-type toy radials replaced by the tabulated AV18 u(k), w(k) and
    the file's own D-state probability.  This is the channel that
    reproduces Cosyn-Weiss II TABLE II quantitatively -- the axial node of
    n_{+-1} at k = 0.2988 GeV/c where w/u = sqrt2, A_zz^wf = -2 there and
    +1 at theta_k = 90 degrees -- which the toy pair cannot, its w/u never
    reaching sqrt2 (test_cosyn_weiss_table_ii_on_av18).  It is a gate, not
    a production channel: nothing the generator ships samples from it.
    """
    k, u, w = _av18_deuteron_tables()
    p_d = P_D_AV18_DEUTERON
    return TaggedChannel(
        spectator.DEUTERON_P_TAG, 1.0, 0.5, 0.5, 1.0,
        (TabulatedWave(0, 1.0 - p_d, tuple(k), tuple(u), "AV18 u(k)"),
         TabulatedWave(2, p_d, tuple(k), tuple(w), "AV18 w(k)")),
        NEUTRON, "d control (AV18 u(k), w(k))")


class TaggedModel:
    """Grid tables of the joint amplitude |A_{m_S}(M; k, cos theta)|^2."""

    def __init__(self, channel, k_max=1.2, nk=280, nc=96):
        self.channel = channel
        self.k = np.linspace(1e-4, k_max, nk)
        self.c = np.linspace(-1.0 + 1.0 / nc, 1.0 - 1.0 / nc, nc)
        self.dk = self.k[1] - self.k[0]
        self.dc = 2.0 / nc
        kappa = channel.base.kappa
        # normalized radial tables: integral psihat^2 k^2 dk = 1
        self._rad = {}
        for w in channel.waves:
            psi = w.radial(self.k, kappa)
            norm = np.sqrt(_trapezoid(psi * psi * self.k**2, self.k))
            self._rad[w.l] = np.sqrt(w.prob) * psi / norm
        self._ms_struck = m_values(channel.s_channel)
        self._amp2 = {}   # (M) -> array (n_mS, nk, nc)

    # --- amplitudes and densities ------------------------------------

    def _amp2_table(self, M):
        """|A_{m_S}|^2 on the (k, c) grid for ion projection M."""
        if M in self._amp2:
            return self._amp2[M]
        ch = self.channel
        kk = self.k[:, None]
        cc = self.c[None, :]
        out = np.zeros((self._ms_struck.size, self.k.size, self.c.size))
        for i, m_s in enumerate(self._ms_struck):
            m_l = M - m_s
            amp = np.zeros((self.k.size, self.c.size))
            for w in ch.waves:
                if abs(m_l) > w.l:
                    continue
                cg = clebsch_gordan(w.l, m_l, ch.s_channel, m_s,
                                    ch.j_ion, M)
                if cg == 0.0:
                    continue
                # phi_L = i^L psihat_L: the momentum-space amplitude of
                # a coordinate-space wave function carries the plane-wave
                # phase, and the stored radials are the PLAIN Bessel
                # transforms (positive at low k, like the AV18 deuteron's
                # tabulated u(k), w(k)), not phased ones.  Only the
                # relative phase is observable; for the same-parity
                # mixtures TaggedChannel allows it is the real
                # (-1)^(L//2) = +1, -1 for L = 0, 2 (and +1 for a lone
                # L = 1 wave, so 7Li is bit-for-bit unaffected).
                phase = (-1.0) ** ((w.l // 2) % 2)
                amp = amp + (phase * self._rad[w.l][:, None] * cg
                             * theta_lm(w.l, int(round(m_l)), cc))
            out[i] = amp * amp
        self._amp2[M] = out
        return out

    def n_of_kc(self, M, k=None, c=None):
        """Spectator density n_M(k, cos theta_k) (spin frame), normalized
        so that integral n k^2 dk dOmega = 1.  Grid table if k, c None."""
        table = self._amp2_table(M).sum(axis=0)
        if k is None:
            return table
        k = np.asarray(k, dtype=float)
        c = np.asarray(c, dtype=float)
        ik = np.clip(np.searchsorted(self.k, k) - 1, 0, self.k.size - 2)
        ic = np.clip(np.searchsorted(self.c, c) - 1, 0, self.c.size - 2)
        return table[ik, ic]

    def struck_populations(self, M):
        """p(m_S | M, k, c) tables, shape (n_mS, nk, nc)."""
        a2 = self._amp2_table(M)
        n = a2.sum(axis=0)
        return a2 / np.maximum(n, 1e-300)

    def population_integrated(self, M):
        """P(m_S | M): khat- and k-integrated channel-spin populations."""
        a2 = self._amp2_table(M)
        raw = (a2 * self.k[None, :, None] ** 2).sum(axis=(1, 2))
        return raw / raw.sum()

    def norm(self, M):
        """Grid quadrature of integral n_M k^2 dk dOmega (should be ~1)."""
        a2 = self._amp2_table(M).sum(axis=0)
        return float((a2 * self.k[:, None] ** 2).sum()
                     * self.dk * 2.0 * np.pi * self.dc)

    def pair_populations(self, m_s):
        """p(m_struck, m_spec | m_S) from the cluster-spin CG factor."""
        ch = self.channel
        m1s = m_values(ch.s_struck)
        m2s = m_values(ch.s_spec)
        out = np.zeros((m1s.size, m2s.size))
        for i, m1 in enumerate(m1s):
            for j, m2 in enumerate(m2s):
                out[i, j] = clebsch_gordan(ch.s_struck, m1, ch.s_spec, m2,
                                           ch.s_channel, m_s) ** 2
        return out

    # --- analytic moments (gates + polarimetry) -----------------------

    def vector_dilution(self, M=None):
        """khat,k-integrated <m_S>/S_c for the stretched state M = J
        (e.g. 1 - (3/2) P_D for the embedded deuteron)."""
        M = self.channel.j_ion if M is None else M
        p = self.population_integrated(M)
        return float((self._ms_struck * p).sum() / self.channel.s_channel)

    def tensor_dilution(self, M=None):
        """khat,k-integrated <3 m_S^2 - 2> for S_c = 1, stretched state."""
        if abs(self.channel.s_channel - 1.0) > 1e-9:
            raise ValueError("tensor dilution defined for S_c = 1")
        M = self.channel.j_ion if M is None else M
        p = self.population_integrated(M)
        return float(((3.0 * self._ms_struck**2 - 2.0) * p).sum())

    def p2_moment(self, M):
        """<P2(cos theta_k)> of the spectator direction for ion state M."""
        n = self.n_of_kc(M)
        w = (n * self.k[:, None] ** 2).sum(axis=0)  # -> c distribution
        p2 = 0.5 * (3.0 * self.c**2 - 1.0)
        return float((w * p2).sum() / w.sum())

    def p2_moment_mixture(self, populations):
        """<P2> for a fill with populations over M (linearity in rho)."""
        ms = m_values(self.channel.j_ion)
        return float(sum(p * self.p2_moment(m)
                         for p, m in zip(populations, ms) if p > 0))

    # --- sampling -------------------------------------------------------

    def sample_kc(self, M, m_s, n, rng):
        """Sample (k, cos theta_k, phi_k) from |A_{m_S}(M)|^2 k^2."""
        i = int(np.flatnonzero(np.isclose(self._ms_struck, m_s))[0])
        dens = self._amp2_table(M)[i] * self.k[:, None] ** 2
        prob = (dens / dens.sum()).ravel()
        cell = rng.choice(prob.size, size=n, p=prob)
        ik, ic = np.unravel_index(cell, dens.shape)
        k = self.k[ik] + (rng.uniform(size=n) - 0.5) * self.dk
        c = self.c[ic] + (rng.uniform(size=n) - 0.5) * self.dc
        phi = rng.uniform(0.0, 2.0 * np.pi, size=n)
        return np.abs(k), np.clip(c, -1.0, 1.0), phi


def boost_spectator(channel, k, c, phi_k, p_per_nucleon,
                    theta_s=0.0, phi_s=0.0):
    """Lab kinematics of the spectator cluster.

    (k, c, phi_k) are spherical components in the SPIN frame; for a tilted
    quantization axis they are rotated to the lab before the longitudinal
    beam boost (same boost algebra as spectator.spectator_lab_kinematics).
    Returns dict with pT, theta, p_lab, R, xL, the lab k components and
    `phi_spec`, the spectator's LAB azimuth.

    `phi_spec` is arctan2(ky, kx) evaluated AFTER the spin-frame rotation
    and is what a planar Roman Pot sees: the boost is longitudinal, so it
    scales pz and leaves (kx, ky) -- hence the azimuth -- untouched.  It
    is the argument `farforward.Optics.clears` needs for the rectangular
    10(sigma_h, sigma_v) envelope of the per-configuration optics; without
    it the cut degenerates to a circle at n sigma_h and overstates the tag
    by 1.7x at the tagging optics (plans/09 B2).  The key is deliberately
    NOT "phi": `TaggedSampler.sample_category` merges this dict into an
    event record whose "phi" is the DIS azimuth.
    """
    base = channel.base
    s = np.sqrt(np.maximum(1.0 - c * c, 0.0))
    kx, ky, kz = (k * s * np.cos(phi_k), k * s * np.sin(phi_k), k * c)
    if theta_s != 0.0 or phi_s != 0.0:
        ct, st = np.cos(theta_s), np.sin(theta_s)
        cp, sp = np.cos(phi_s), np.sin(phi_s)
        # R_z(phi_s) R_y(theta_s): spin frame -> lab
        kx, kz = ct * kx + st * kz, -st * kx + ct * kz
        kx, ky = cp * kx - sp * ky, sp * kx + cp * ky
    m = base.m_spec
    e_rest = np.sqrt(m * m + kx * kx + ky * ky + kz * kz)
    # the physical nuclear mass, not beam_A * M_U: the docstring promises
    # the same boost algebra as spectator.spectator_lab_kinematics, and
    # that function moved to NUCLEUS_MASS on 2026-08-26 (plans/08 C1).
    # Leaving A * M_U here would have given the same beam two masses in
    # two modules, 2.2e-3 apart for 6Li.
    m_beam = base.m_beam
    p_beam = base.beam_A * p_per_nucleon
    e_beam = np.sqrt(p_beam**2 + m_beam**2)
    gamma = e_beam / m_beam
    gbeta = p_beam / m_beam
    pz_lab = gamma * kz + gbeta * e_rest
    pt = np.sqrt(kx * kx + ky * ky)
    p_lab = np.sqrt(pt * pt + pz_lab * pz_lab)
    theta = np.arctan2(pt, pz_lab)
    rigidity_beam = p_beam / base.beam_Z
    if base.spectator_Z > 0:
        rig = (p_lab / base.spectator_Z) / rigidity_beam
    else:
        rig = np.full_like(p_lab, np.nan)
    return {"pT": pt, "theta": theta, "p_lab": p_lab, "R": rig,
            "xL": p_lab / (base.spectator_A * p_per_nucleon),
            "kx": kx, "ky": ky, "kz": kz,
            "phi_spec": np.arctan2(ky, kx)}


def acceptance_weights(model, config, optics, n_phi=64, theta_s=0.0,
                       phi_s=0.0):
    """eps(k, cos theta_k) on a TaggedModel's grid: the fraction of lab
    azimuths at which a spectator of that (k, c) is Roman-Pot accepted
    (main window or near-beam tail).

    theta, R and pT depend on (k, c) alone -- the beam boost is
    longitudinal and does not touch (kx, ky) -- so the azimuth enters only
    through the rectangular envelope, and a uniform phi average is the
    exact marginal (the spectator's azimuth is uniform and independent of
    the DIS azimuth).  Returned shape is (model.k.size, model.c.size).
    """
    phi = (np.arange(n_phi) + 0.5) * 2.0 * np.pi / n_phi
    kk, cc, pp = np.meshgrid(model.k, model.c, phi, indexing="ij")
    lab = boost_spectator(model.channel, kk.ravel(), cc.ravel(), pp.ravel(),
                          config.ion_momentum_per_nucleon, theta_s, phi_s)
    route = route_charged(lab["R"], lab["theta"], lab["pT"], optics,
                          phi=lab["phi_spec"],
                          pot_config=yr_config_key(config))
    return ((route == 1) | (route == 4)).reshape(kk.shape).mean(axis=2)


def azz_tensor_curve(model, ic=None, weights=None):
    """The wave-function tensor asymmetry of a spin-1 channel vs k,

        A_zz^wf = (n_+1 + n_-1 - 2 n_0) / (n_+1 + n_-1 + n_0),

    at one cos theta_k cell (`ic`, e.g. the 90 degree slice the analytic
    panel of money plot 4 draws) or ACCEPTANCE-WEIGHTED: with `weights` an
    (nk, nc) table -- `acceptance_weights` -- both sums are integrated over
    cos theta_k against it, which is the prediction for a sample the
    far-forward acceptance sculpts in theta_k.

    The distinction is not cosmetic.  The accepted 6Li alpha sample is not
    at theta_k = 90 degrees: at the Yellow Report optics only the
    off-rigidity R < 0.95 window slice survives and it is longitudinal
    (<|cos theta_k|> = 0.71-0.79), while the near-beam tail that the
    tagging optics opens is transverse (0.40).  Comparing markers from
    either against the 90 degree curve reads the S/D interference at the
    wrong angle (plans/09 B2).  A `weights` concentrated in a single c cell
    returns the `ic` curve exactly.
    """
    n = {m: model.n_of_kc(m) for m in (1.0, 0.0, -1.0)}
    if weights is None:
        n = {m: v[:, ic] for m, v in n.items()}
    else:
        n = {m: (v * weights).sum(axis=1) for m, v in n.items()}
    num = n[1.0] + n[-1.0] - 2.0 * n[0.0]
    den = n[1.0] + n[-1.0] + n[0.0]
    return np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)


class TaggedSampler:
    """Spin-correlated (e', spectator) events for one tagged channel.

    Composes the two-cluster spin model with a Step-5.A InclusiveSampler
    built on the STRUCK cluster's structure functions at the ion beam's
    per-nucleon momentum.  Every event carries DIS kinematics, ion and
    struck-cluster spin labels, the spectator momentum (spin frame and
    lab), and its far-forward route.
    """

    def __init__(self, model, kernel_struck, beam_config, scenario=None,
                 optics=None, **sampler_kw):
        """`optics` defaults to the Yellow Report HIGH-ACCEPTANCE optics OF
        THIS CONFIGURATION (`farforward.yr_optics`), not to the legacy
        proton-derived 73 microrad: since 2026-08-28 the envelope is
        per-configuration and anisotropic (plans/10 A3, plans/09 B2), and a
        module-level default cannot know the beam.  Pass
        `farforward.HIGH_ACCEPTANCE` / `HIGH_DIVERGENCE` explicitly to
        reproduce a pre-2026-08-28 number.  The route is only a convenience
        column on the event record -- every script re-routes with its own
        optics -- but the default is what a caller that does not think
        about it gets, so it is the current one."""
        self.model = model
        ch = model.channel
        if kernel_struck.ion is not ch.dis_target:
            raise ValueError("kernel target %r != channel DIS target %r"
                             % (kernel_struck.ion.name, ch.dis_target.name))
        self.config = beam_config
        self.optics = (yr_optics(beam_config, "high-acceptance")
                       if optics is None else optics)
        dis_config = beams.BeamConfig(beam_config.electron_energy,
                                      ch.dis_target,
                                      beam_config.ion_momentum_per_nucleon)
        self.inner = InclusiveSampler(kernel_struck, dis_config,
                                      scenario=scenario, **sampler_kw)
        self._pure = {}

    def _pure_category(self, m_s, lam_e, pe):
        key = (m_s, lam_e, pe)
        if key not in self._pure:
            ms = m_values(self.model.channel.s_channel)
            pops = tuple(1.0 if np.isclose(m, m_s) else 0.0 for m in ms)
            self._pure[key] = bk.SpinCategory(
                "struck m=%g lam=%d" % (m_s, lam_e),
                self.model.channel.s_channel, pops, lam_e=lam_e, pe=pe)
        return self._pure[key]

    def sigma_tot_pb(self, category):
        """Accepted tagged cross section [pb] for an ION spin category
        (bookkeeping.SpinCategory with j = j_ion).  Tagging probability is
        1 in the two-cluster model; purity/efficiency come from routing."""
        tot = 0.0
        ms_ion = m_values(self.model.channel.j_ion)
        for pop, M in zip(category.populations, ms_ion):
            if pop <= 0.0:
                continue
            p_ms = self.model.population_integrated(M)
            for p, m_s in zip(p_ms, m_values(self.model.channel.s_channel)):
                if p <= 0.0:
                    continue
                cat = self._pure_category(m_s, category.lam_e, category.pe)
                tot += pop * p * self.inner.sigma_tot_pb(cat)
        return tot

    def sample_category(self, category, lumi_pb=None, n=None, rng=None):
        """Tagged events for one ion spin category.

        Joint sampling: M from the fill populations weighted by the
        (m_S-summed) DIS rates; then (m_S, k, khat) from
        |A_{m_S}(M)|^2 * S_DIS(m_S); then (x, Q2, phi) from the inclusive
        sampler conditioned on m_S; then lab boost + far-forward route.
        """
        rng = rng or np.random.default_rng(20260713)
        ch = self.model.channel
        ms_ion = m_values(ch.j_ion)
        ms_c = m_values(ch.s_channel)
        # per-(M, m_S) rates in pb
        rates = np.zeros((ms_ion.size, ms_c.size))
        for a, M in enumerate(ms_ion):
            pop = category.populations[a]
            if pop <= 0.0:
                continue
            p_ms = self.model.population_integrated(M)
            for b, m_s in enumerate(ms_c):
                if p_ms[b] <= 0.0:
                    continue
                cat = self._pure_category(m_s, category.lam_e, category.pe)
                rates[a, b] = pop * p_ms[b] * self.inner.sigma_tot_pb(cat)
        if lumi_pb is not None:
            counts = rng.poisson(lumi_pb * rates)
        elif n is not None:
            counts = rng.multinomial(n, (rates / rates.sum()).ravel()
                                     ).reshape(rates.shape)
        else:
            raise ValueError("pass lumi_pb or n")

        chunks = []
        for a, M in enumerate(ms_ion):
            for b, m_s in enumerate(ms_c):
                cnt = int(counts[a, b])
                if cnt == 0:
                    continue
                cat = self._pure_category(m_s, category.lam_e, category.pe)
                dis = self.inner.sample_category(cat, n=cnt, rng=rng)
                k, c, phi_k = self.model.sample_kc(M, m_s, cnt, rng)
                lab = boost_spectator(ch, k, c, phi_k,
                                      self.config.ion_momentum_per_nucleon,
                                      category.theta_s, category.phi_s)
                chunk = {"x": dis["x"], "q2": dis["q2"], "y": dis["y"],
                         "phi": dis["phi"], "cell": dis["cell"],
                         "m_ion": np.full(cnt, float(M)),
                         "m_struck": np.full(cnt, float(m_s)),
                         "k": k, "cos_theta_k": c, "phi_k": phi_k}
                chunk.update(lab)
                chunks.append(chunk)
        keys = chunks[0].keys() if chunks else ()
        out = {kk: np.concatenate([c[kk] for c in chunks]) for kk in keys}
        if chunks:
            # the LAB azimuth of the spectator, not the DIS azimuth: the
            # near-beam envelope is a rectangle (Optics.clears), and
            # without phi it degenerates to a circle at n sigma_h
            out["route"] = route_charged(out["R"], out["theta"], out["pT"],
                                         self.optics, phi=out["phi_spec"],
                                         pot_config=yr_config_key(self.config))
        out["category"] = category.name
        out["lam_e"] = category.lam_e
        return out


def rp_accepted(events):
    """Roman-Pot mask: main window (route 1) + near-beam pT tail (4)."""
    return (events["route"] == 1) | (events["route"] == 4)
