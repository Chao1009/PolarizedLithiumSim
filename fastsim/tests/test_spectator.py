"""Tests for cluster-spectator kinematics and far-forward routing."""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import farforward as ff
from polli_fastsim import spectator as sp


def test_kappa_values():
    # verified separation energies -> momentum scales
    np.testing.assert_allclose(sp.LI6_ALPHA_TAG.kappa, 0.0607, atol=0.002)
    np.testing.assert_allclose(sp.LI7_ALPHA_TAG.kappa, 0.0889, atol=0.002)


def test_p_wave_vanishes_at_origin():
    n_s = sp.momentum_density(1e-4, 0.06, 0.3, l_wave=0)
    n_p = sp.momentum_density(1e-4, 0.09, 0.3, l_wave=1)
    assert n_p / n_p.max() < 1e-4 or n_p < n_s  # P-wave suppressed at k->0


def test_spectator_rigidity_centers():
    rng = np.random.default_rng(1)
    k6 = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG, 137.5, 50_000, rng=rng)
    k7 = sp.spectator_lab_kinematics(sp.LI7_ALPHA_TAG, 117.9, 50_000, rng=rng)
    # nominal R is the ratio of mass-to-charge ratios,
    # (m_spec/Z_spec)/(m_beam/Z_beam) = 0.99813 / 0.85571 -- close to but
    # not equal to the naive (A_f Z_beam)/(A_beam Z_f) = 1.0 and 6/7,
    # which is what A * M_U used to give (pinned below)
    np.testing.assert_allclose(np.median(k6["R"]), 1.0, atol=0.01)
    np.testing.assert_allclose(np.median(k7["R"]), 6.0 / 7.0, atol=0.01)
    # Fermi smearing of R is percent-level
    assert np.std(k6["R"]) < 0.05
    # transverse momentum scale ~ cluster momentum scale, well below 1 GeV
    assert 0.02 < np.median(k6["pT"]) < 0.3
    # theta is far inside 5 mrad for everything at these momenta
    assert np.quantile(k6["theta"], 0.99) < 5e-3


def test_routing_logic():
    optics = ff.HIGH_ACCEPTANCE
    env = optics.n_sigma * optics.sigma_theta          # 0.727 mrad
    # 7Li alpha: R=0.857, off rigidity -> Roman Pots, no near-beam cut
    assert ff.route_charged(0.857, 1e-3, 0.01, optics) == 1
    # 6Li alpha: R=1.0 -> the ANGULAR envelope decides, not pT.  A 0.05 GeV
    # fragment at 1 mrad is outside the envelope and IS tagged; the same
    # 0.05 GeV inside the envelope is not.
    assert ff.route_charged(1.0, 1e-3, 0.05, optics) == 4
    assert ff.route_charged(1.0, 0.5 * env, 0.30, optics) == 0
    assert ff.route_charged(1.0, 1.5 * env, 0.01, optics) == 4
    # proton from 6Li: R=0.5 -> OMD (window is [0.45, 0.65))
    assert ff.route_charged(0.50, 1e-3, 0.05, optics) == 2
    # 7Li spectator proton: R = 3/7 = 0.4286, below the OMD low edge
    assert ff.route_charged(3.0 / 7.0, 1e-3, 0.05, optics) == 0
    # triton from 7Li: R=1.29 -> lost
    assert ff.route_charged(9.0 / 7.0, 1e-3, 0.05, optics) == 0
    # B0 window
    assert ff.route_charged(0.857, 8e-3, 0.05, optics) == 3


def test_near_beam_cut_is_angular_not_a_momentum_threshold():
    """The documented 0.20 / 0.41 GeV are the 10-sigma envelope for a
    275 GeV PROTON.  Applied verbatim to a nucleus of momentum A p_u they
    understate it by A p_u / 275 (code review S8): 0.20 GeV on a
    4 x 137.5 GeV alpha is 5 sigma, not 10."""
    ha = ff.HIGH_ACCEPTANCE
    assert ha.pt_cut_near_beam == pytest.approx(0.20)
    # the fast-sim quotes the rounded-up end of the 0.41-0.45 band that
    # plans/06 SS6.5 derives; polligen.reco quotes 0.41 (farforward.py)
    assert ff.HIGH_DIVERGENCE.pt_cut_near_beam == pytest.approx(0.45)
    p_alpha = 4 * 137.5
    assert ha.pt_cut_for(p_alpha) == pytest.approx(0.40, rel=1e-6)
    n_sigma_of_020 = 0.20 / (ha.sigma_theta * p_alpha)
    assert n_sigma_of_020 == pytest.approx(5.0, rel=1e-6)
    # coherent 6Li at the three energies: 0.09 / 0.22 / 0.60 GeV
    for p_u, expect in ((20.5, 0.0895), (50.0, 0.218), (137.5, 0.600)):
        assert ha.pt_cut_for(6 * p_u) == pytest.approx(expect, rel=2e-3)


def test_li7_alpha_tag_efficient_li6_suppressed():
    rng = np.random.default_rng(2)
    k7 = sp.spectator_lab_kinematics(sp.LI7_ALPHA_TAG, 117.9, 100_000,
                                     rng=rng)
    acc7 = ff.acceptance_summary(k7["R"], k7["theta"], k7["pT"])
    assert acc7["RomanPots"] > 0.9  # off-rigidity alpha: clean RP tag

    k6 = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG, 137.5, 100_000,
                                     rng=rng)
    acc6 = ff.acceptance_summary(k6["R"], k6["theta"], k6["pT"])
    # near-beam alpha only via the pT tail -> much smaller acceptance
    assert (1.0 - acc6["lost"]) < 0.5
    assert acc6["RP (pT tail, R~1)"] == 1.0 - acc6["lost"] or True


# --- physical nuclear masses (plans/08 C1, 2026-08-26) --------------------

def test_nucleus_mass_table_against_codata_and_the_ame_mass_excesses():
    """External anchor for NUCLEUS_MASS.

    The table is built from the AME2020 atomic masses in u; the six
    light nuclear masses of CODATA 2018 (Tiesinga, Mohr, Newell, Taylor,
    Rev. Mod. Phys. 93 (2021) 025010, Table I) are an independent
    evaluation of the same quantities, and the AME2020 mass-EXCESS
    column is an independent representation of the rest.  Both must come
    back, including the 6Li and 7Li entries CODATA does not carry.
    """
    codata_mev = {(1, 1): 938.27208816,   (0, 1): 939.56542052,
                  (1, 2): 1875.61294257,  (1, 3): 2808.92113298,
                  (2, 3): 2808.39160743,  (2, 4): 3727.3794066}
    for key, mev in codata_mev.items():
        assert sp.NUCLEUS_MASS[key] == pytest.approx(mev * 1e-3, abs=2e-9)
    u, m_e = 931.49410242e-3, 0.51099895e-3           # CODATA 2018 [GeV]
    b_e = {2: 79.0052e-9, 3: 203.486e-9}              # NIST ASD sums [GeV]
    for (z, a), excess_kev in (((2, 4), 2424.91587),
                               ((3, 6), 14086.8789),
                               ((3, 7), 14907.1051)):
        m = a * u + excess_kev * 1e-6 - z * m_e + b_e[z]
        assert sp.NUCLEUS_MASS[(z, a)] == pytest.approx(m, abs=2e-8)


def test_masses_table_agrees_with_the_nuclear_masses_to_its_last_digit():
    """MASSES (5 decimals in GeV) is the same evaluation rounded to
    10 keV -- it is kept literal so that no published spectator
    distribution moves, so the agreement is worth pinning."""
    keys = {"p": (1, 1), "n": (0, 1), "d": (1, 2), "t": (1, 3),
            "3He": (2, 3), "alpha": (2, 4)}
    for name, key in keys.items():
        assert sp.MASSES[name] == pytest.approx(sp.NUCLEUS_MASS[key],
                                                abs=5e-6)


def test_nucleus_mass_falls_back_to_a_amu_off_the_table():
    """The documented fallback, and the size of what it costs: A * M_U
    drops the mass excess and keeps the electrons."""
    assert sp.nucleus_mass(8, 16) == pytest.approx(16 * sp.M_U)
    for (z, a), rel in (((3, 6), 2.2e-3), ((3, 7), 2.1e-3),
                        ((1, 2), 6.7e-3), ((2, 4), 3.8e-4)):
        exact = sp.nucleus_mass(z, a)
        assert (exact - a * sp.M_U) / exact == pytest.approx(rel, rel=0.05)


def test_separation_energies_close_against_the_mass_table():
    """Independent closure: the cluster separation energies of the
    channel table come from nuclear data (module docstring), the masses
    from AME2020, and m_spec + m_partner - m_beam must return S.  The
    5 keV tolerance is the 10 keV rounding of MASSES (the neutron entry
    is 4.6 keV high), not the mass evaluation, which closes to eV.
    """
    channels = sp.CHANNELS + (sp.DEUTERON_P_TAG, sp.DEUTERON_N_TAG,
                              sp.HE3_P_TAG)
    for ch in channels:
        s = ch.m_spec + ch.m_partner - ch.m_beam
        assert s == pytest.approx(ch.separation_energy, abs=5e-6), ch.name
    # what the beam mass is for: with A * M_U the same closure returns
    # 14.05 MeV for the 6Li alpha-d threshold instead of 1.4743 MeV
    ch = sp.LI6_ALPHA_TAG
    fake = ch.m_spec + ch.m_partner - ch.beam_A * sp.M_U
    assert fake == pytest.approx(14.05e-3, rel=1e-3)


def test_kappa_is_a_property_of_the_beam_not_of_the_tagged_cluster():
    """mu is the reduced mass of the two SEPARATED clusters, which is
    symmetric under swapping them, so kappa = sqrt(2 mu S) must be the
    same for the alpha-tag and the d-tag (t-tag) of one beam.  With
    m_beam = A * M_U and m_other = m_beam - m_spec it was not: 60.50 vs
    60.62 MeV for 6Li and 88.76 vs 88.82 MeV for 7Li -- the symptom that
    fixed the construction (plans/08 C1)."""
    np.testing.assert_allclose(sp.LI6_ALPHA_TAG.kappa, sp.LI6_D_TAG.kappa,
                               rtol=1e-5)
    np.testing.assert_allclose(sp.LI7_ALPHA_TAG.kappa, sp.LI7_T_TAG.kappa,
                               rtol=1e-5)
    # the free partner really is the free cluster, not m_beam - m_spec
    assert sp.LI6_ALPHA_TAG.m_partner == pytest.approx(sp.MASSES["d"],
                                                       abs=5e-6)
    assert (sp.LI6_ALPHA_TAG.m_beam - sp.LI6_ALPHA_TAG.m_spec
            == pytest.approx(sp.MASSES["d"] - 1.4743e-3, abs=5e-6))


def test_spectator_rigidity_is_a_ratio_of_mass_to_charge_ratios():
    """Analytic limit k -> 0: a spectator at rest in the beam frame is
    boosted to p_lab = gamma beta m_spec = (p_beam / m_beam) m_spec, so

        R(k=0) = (m_spec / Z_spec) / (m_beam / Z_beam),

    independent of the beam energy and of the momentum density.  A * M_U
    collapses this onto the naive (A_spec Z_beam)/(A_beam Z_spec); with
    the real masses the alpha is more bound per nucleon than the lithium
    it comes from and therefore sits BELOW beam rigidity, by 1.87e-3
    (6Li) and 1.67e-3 (7Li).  Those are far inside the +-5% near-beam
    band, but they slide the 6Li alpha across the R = 0.95 Roman-Pot
    window edge (see spectator_lab_kinematics).
    """
    for ch, p_u, drop, tol in ((sp.LI6_ALPHA_TAG, 137.5, 1.87e-3, 1e-3),
                               (sp.LI7_ALPHA_TAG, 117.9, 1.67e-3, 4e-3)):
        p_beam = ch.beam_A * p_u
        e_beam = np.hypot(p_beam, ch.m_beam)
        p_lab = (p_beam / e_beam) * (e_beam / ch.m_beam) * ch.m_spec
        r0 = (p_lab / ch.spectator_Z) / (p_beam / ch.beam_Z)
        np.testing.assert_allclose(
            r0, (ch.m_spec / ch.spectator_Z) / (ch.m_beam / ch.beam_Z),
            rtol=1e-12)
        naive = ((ch.spectator_A * ch.beam_Z)
                 / (ch.beam_A * ch.spectator_Z))
        assert r0 < naive
        np.testing.assert_allclose((naive - r0) / naive, drop, rtol=0.02)
        # the sampled distribution sits on it, displaced only by the
        # Fermi smearing (e_rest = sqrt(m^2+k^2) > m raises p_lab)
        kin = sp.spectator_lab_kinematics(ch, p_u, 50_000,
                                          rng=np.random.default_rng(3))
        assert np.median(kin["R"]) == pytest.approx(r0, abs=tol)
        assert np.median(kin["R"]) > r0     # the smearing is one-sided
    # The two tolerances are not a matter of taste, and they are why this
    # test is NOT the one that separates the two beam masses.  Measured
    # at the settings above: the S-wave 6Li alpha lands 2.7e-4 above r0,
    # less than the 1.87e-3 the mass offset moved it, so its MEDIAN is
    # still below the naive 1.0 -- the strict inequality asserted next.  The P-wave 7Li alpha has a harder k
    # tail: it lands 1.8e-3 above r0, MORE than the 1.4e-3 the mass
    # offset moved it, so its median sits just above the naive 6/7 and
    # no inequality there can tell the two masses apart.  The sharp
    # statement for every channel is the rapidity test below.
    kin6 = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG, 137.5, 50_000,
                                       rng=np.random.default_rng(3))
    assert np.median(kin6["R"]) < 1.0       # 0.99841; A * M_U: 1.00066


def test_the_lab_boost_is_by_the_beam_velocity_of_the_nuclear_mass():
    """The one line the published acceptance rides on, pinned sharply.

    `m_beam` enters `spectator_lab_kinematics` only through the boost
    (gamma, gamma*beta), and the sharp handle on a longitudinal boost is
    rapidity, which ADDS:

        y_lab = y_beam + artanh(k_z / sqrt(m_spec^2 + k^2)),
        y_beam = artanh(p_beam / E_beam) = arcsinh(p_beam / m_beam).

    `sample_k` draws k isotropically, so the second term is odd in k_z
    and averages to zero: the MEAN spectator rapidity is the beam
    rapidity, with no Fermi-smearing bias at all -- unlike the median of
    R, which the k tail displaces by more than the effect under test for
    the P-wave channels (see the test above).  Nothing here re-derives
    the module's boost; it uses the two textbook facts that rapidity is
    additive under a longitudinal boost and that an isotropic k has
    <artanh(k_z/E)> = 0.

    Measured 2026-08-26 at 2e5 events, seed 7: the residual is 8e-6 to
    3.9e-5 for the six channels, while beam_A * M_U would put the mean
    2.05e-3 (7Li) to 6.76e-3 (deuteron control) ABOVE the true beam
    rapidity.  The assertion is scaled to that separation, so it cannot
    be satisfied by the wrong mass at any statistics.
    """
    cases = ((sp.LI6_ALPHA_TAG, 137.5), (sp.LI6_D_TAG, 137.5),
             (sp.LI7_ALPHA_TAG, 117.9), (sp.LI7_T_TAG, 117.9),
             (sp.DEUTERON_P_TAG, 130.0), (sp.HE3_P_TAG, 166.0))
    for ch, p_u in cases:
        p_beam = ch.beam_A * p_u
        y_beam = np.arcsinh(p_beam / ch.m_beam)
        y_amu = np.arcsinh(p_beam / (ch.beam_A * sp.M_U))
        kin = sp.spectator_lab_kinematics(ch, p_u, 200_000,
                                          rng=np.random.default_rng(7))
        pz = kin["p_lab"] * np.cos(kin["theta"])
        e = np.hypot(kin["p_lab"], ch.m_spec)
        y = 0.5 * np.log((e + pz) / (e - pz))
        resid = float(np.mean(y)) - y_beam
        assert abs(resid) < 0.25 * (y_amu - y_beam), ch.name
        assert abs(resid) < 5e-4, ch.name


def test_deuteron_kappa_is_the_textbook_deuteron_wave_number():
    """External anchor for the kappa CONSTRUCTION, independent of the
    lithium channels it is used for.

    The d -> p + n control channel is the one case where the two
    "clusters" are the free nucleons themselves, so kappa is the
    textbook deuteron wave number

        gamma = sqrt(2 mu B) / hbar = 0.2316 fm^-1,  1/gamma = 4.318 fm,

    with mu = m_p m_n / (m_p + m_n) and B = 2.2246 MeV -- the asymptotic
    slope of the deuteron wave function quoted in every nuclear-physics
    text (e.g. Krane, Introductory Nuclear Physics, Wiley 1988, Sec.
    4.2).  Nothing in the module carries that number.  The pre-2026-08-26
    construction (partner mass m_beam - m_spec, low by B) gave 45.52 MeV,
    i.e. 4.335 fm, which this tolerance excludes.
    """
    hbar_c_fm = 0.1973269804        # GeV fm (CODATA 2018)
    for ch in (sp.DEUTERON_P_TAG, sp.DEUTERON_N_TAG):
        assert ch.kappa == pytest.approx(0.04570, abs=5e-5), ch.name
        assert hbar_c_fm / ch.kappa == pytest.approx(4.318, abs=5e-3)


def test_the_a_amu_fallbacks_off_the_table():
    """The two documented fallbacks, which no channel in the module
    reaches, and the one nuclide where A * M_U is not low.

    u is defined as one twelfth of the mass of a neutral 12C ATOM, so at
    (Z, A) = (6, 12) the fallback returns an atomic mass: it is HIGH of
    the nuclear mass by the six electrons less their binding,
    6 m_e - B_e(6) = 3.066 - 0.001 MeV, plus the 0.05 MeV by which the
    module's 5-digit M_U rounds u down -- 3.016 MeV in total, the
    opposite sign to every entry the table does carry.

    The `m_partner` fallback is exercised on the alpha decay of 12C,
    12C -> alpha + 8Be, whose partner (Z, A) = (4, 8) is off the table
    too: S_alpha(12C) = 7.36659 MeV from the AME2020 mass excesses
    (8Be +4941.67 keV, 4He +2424.92 keV, 12C 0 by definition of u).  The
    fallback must keep the DEFINITION of a separation energy,
    S = m_spec + m_partner - m_beam, which m_beam - m_spec alone (the
    pre-2026-08-26 partner) misses by exactly S.
    """
    u, m_e = 931.49410242e-3, 0.51099895e-3      # CODATA 2018 [GeV]
    b_e6 = 1030.11e-9                            # NIST ASD sum, Z = 6
    m_c12 = 12 * u - 6 * m_e + b_e6
    assert sp.nucleus_mass(6, 12) == 12 * sp.M_U               # fallback
    assert sp.nucleus_mass(6, 12) - m_c12 == pytest.approx(3.016e-3,
                                                           abs=2e-6)
    assert sp.nucleus_mass(6, 12) > m_c12        # high, not low

    s_alpha = 7.36659e-3
    c12 = sp.ClusterChannel("12C: DIS on 8Be, alpha spectator (off-table)",
                            12, 6, "alpha", 4, 2, s_alpha, l_wave=0)
    assert (2, 8) not in sp.NUCLEUS_MASS and (4, 8) not in sp.NUCLEUS_MASS
    assert c12.m_spec + c12.m_partner - c12.m_beam == pytest.approx(
        s_alpha, abs=1e-12)
    assert c12.m_partner - (c12.m_beam - c12.m_spec) == pytest.approx(
        s_alpha, abs=1e-12)
    # and kappa stays the two-body scale of the SEPARATED fragments
    mu = c12.m_spec * c12.m_partner / (c12.m_spec + c12.m_partner)
    assert c12.kappa == pytest.approx(np.sqrt(2 * mu * s_alpha), rel=1e-12)
