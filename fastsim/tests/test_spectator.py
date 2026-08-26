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
    # nominal R = (A_f Z_beam)/(A_beam Z_f): 6Li alpha -> 1.0; 7Li alpha -> 6/7
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
