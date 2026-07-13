"""Tagged mode: spin-correlated (e', spectator) events (plans/05 step 5.B).

Model: two-cluster light-front-flavored impulse approximation.  The ion
ground state is expanded in cluster relative partial waves,

  |J M> = sum_L a_L sum_{m_L, m_S} <L m_L S_c m_S | J M>
            psihat_L(k) Y_L^{m_L}(khat) |S_c m_S>,

with S_c the channel spin (the coupled cluster spins), psihat_L normalized
radial waves and a_L^2 = P_L the wave probabilities.  Everything the
tagged observables need follows from the joint amplitude

  A_{m_S}(M; k, khat) = sum_L a_L psihat_L(k) C_L(M, m_S) Y_L^{M-m_S}(khat):

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
"""

from dataclasses import dataclass
from typing import Tuple

import numpy as np

from polli_fastsim import beams, spectator
from polli_fastsim.farforward import HIGH_ACCEPTANCE, route_charged

from . import bookkeeping as bk
from .sample import InclusiveSampler
from .spin import clebsch_gordan, m_values

# struck-cluster DIS targets that are not beam species
TRITON = beams.Ion("t", 3, 1, 0.5, eff_pol_p=0.86, eff_pol_n=-0.028)
NEUTRON = beams.Ion("n", 1, 0, 0.5, eff_pol_p=0.0, eff_pol_n=1.0)


def theta_lm(l, m, c):
    """|m|-azimuth-stripped spherical harmonic Theta_l^m(theta) with
    Condon-Shortley signs: Y_l^m = Theta_l^m(theta) exp(i m phi)."""
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


# --- default channels (radial betas carry the model band, as in fastsim) --

# alpha-d D-state probability chosen so the embedded-deuteron vector
# dilution 1 - (3/2) P_D reproduces the 0.87 of b1_li6_from_deuteron
# (SCENARIO -- VMC overlaps are the scheduled replacement, plans/04 #15).
P_D_LI6 = 0.0867
# deuteron D-state probability (AV18-like)
P_D_DEUTERON = 0.045


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
            norm = np.sqrt(np.trapz(psi * psi * self.k**2, self.k))
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
                amp = amp + (self._rad[w.l][:, None] * cg
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
    Returns dict with pT, theta, p_lab, R, xL, and the lab k components.
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
    m_beam = base.beam_A * spectator.M_U
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
            "kx": kx, "ky": ky, "kz": kz}


class TaggedSampler:
    """Spin-correlated (e', spectator) events for one tagged channel.

    Composes the two-cluster spin model with a Step-5.A InclusiveSampler
    built on the STRUCK cluster's structure functions at the ion beam's
    per-nucleon momentum.  Every event carries DIS kinematics, ion and
    struck-cluster spin labels, the spectator momentum (spin frame and
    lab), and its far-forward route.
    """

    def __init__(self, model, kernel_struck, beam_config, scenario=None,
                 optics=HIGH_ACCEPTANCE, **sampler_kw):
        self.model = model
        ch = model.channel
        if kernel_struck.ion is not ch.dis_target:
            raise ValueError("kernel target %r != channel DIS target %r"
                             % (kernel_struck.ion.name, ch.dis_target.name))
        self.config = beam_config
        self.optics = optics
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
            out["route"] = route_charged(out["R"], out["theta"], out["pT"],
                                         self.optics)
        out["category"] = category.name
        out["lam_e"] = category.lam_e
        return out


def rp_accepted(events):
    """Roman-Pot mask: main window (route 1) + near-beam pT tail (4)."""
    return (events["route"] == 1) | (events["route"] == 4)
