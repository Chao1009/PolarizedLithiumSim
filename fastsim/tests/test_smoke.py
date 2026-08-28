"""Smoke tests: kinematic self-consistency, sane magnitudes, FOM scaling."""

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim import kinematics as kin
from polli_fastsim.polarized import ToyG1, toy_b1, toy_delta_gluon
from polli_fastsim.structure import NuclearF2, ToyF2, dsigma_dx_dq2


def test_beam_energies():
    """The rigidity cap, and the top configuration that it sets.

    `momentum_per_nucleon_max` is the ring's 917.3 T*m limit expressed per
    nucleon, Z/A times the 275 GeV proton: 117.86 GeV/u for 7Li (Z/A =
    3/7) and 137.5 for 6Li and the deuteron (1/2).  At the TOP
    configuration the cap binds before gamma-matching does -- a 275 GeV
    proton is gamma = 293, a 6Li at 137.5 GeV/u only 147 -- so the cap is
    also the published top energy, 117.9 GeV/u for 7Li and 137.5 for 6Li
    (plans/10 A0).  The assertion on the first line was inert until
    2026-08-28: `assert x == assert_allclose(...) or True` is `assert
    False or True`."""
    np.testing.assert_allclose(beams.LI7.momentum_per_nucleon_max,
                               275.0 * 3 / 7)
    np.testing.assert_allclose(beams.LI7.momentum_per_nucleon_max, 117.857,
                               atol=5e-4)
    np.testing.assert_allclose(beams.LI6.momentum_per_nucleon_max, 137.5)
    np.testing.assert_allclose(beams.DEUTERON.momentum_per_nucleon_max, 137.5)
    # and the cap is what the top configuration returns, rounded to the
    # 0.1 GeV/u the published menu is quoted at
    for name, top in (("7Li", 117.9), ("6Li", 137.5)):
        cap = beams.IONS[name].momentum_per_nucleon_max
        got = beams.default_configs(name)[2].ion_momentum_per_nucleon
        np.testing.assert_allclose(got, cap, atol=0.05)
        np.testing.assert_allclose(got, top, atol=1e-9)
    cfg = beams.BeamConfig(18.0, beams.LI6, 137.5)
    np.testing.assert_allclose(cfg.sqrt_s_per_nucleon, np.sqrt(4 * 18 * 137.5))


def test_kinematic_roundtrip():
    s = kin.sqrt_s(18.0, 117.9) ** 2
    x, y = 0.01, 0.5
    q2 = kin.q2_from_xy(x, y, s)
    e_p, theta, eta = kin.scattered_electron(x, y, s, 18.0)
    # recompute (x, Q2) from (E', theta)
    q2_back = 2.0 * 18.0 * e_p * (1.0 + np.cos(theta))
    y_back = 1.0 - (e_p / (2 * 18.0)) * (1.0 - np.cos(theta))
    np.testing.assert_allclose(q2_back, q2, rtol=1e-9)
    np.testing.assert_allclose(y_back, y, rtol=1e-9)
    assert eta < 0  # scattered electron goes electron-side at low x


def test_f2_magnitudes():
    f2 = ToyF2()
    # anchors good to ~30% -- this is a TOY model by design
    assert 0.8 < f2.f2p(1e-3, 10.0) < 1.6
    assert 0.25 < f2.f2p(0.1, 10.0) < 0.5
    assert 0.05 < f2.f2p(0.4, 10.0) < 0.2
    assert f2.f2n(0.5, 10.0) < f2.f2p(0.5, 10.0)


def test_cross_section_positive_and_scale():
    s = 4 * 18 * 117.9
    nf2 = NuclearF2(beams.LI7)
    x, q2 = 0.1, 10.0
    xs = dsigma_dx_dq2(x, q2, s, nf2.f2a(x, q2) / 7.0)
    assert xs > 0
    # per-nucleon ep-like DIS xsec at this point: O(1e3-1e5) pb/GeV^2
    assert 1e2 < xs < 1e6


def test_rate_projection_and_fom_scaling():
    cfg = beams.BeamConfig(10.0, beams.LI7, 100 * 3 / 7.0)
    sc1 = fom.Scenario(lumi_fb_per_nucleon=10.0)
    sc2 = fom.Scenario(lumi_fb_per_nucleon=40.0)
    p1 = fom.project_rates(cfg, sc1)
    p2 = fom.project_rates(cfg, sc2)
    assert p1.n_events.sum() > 1e6  # DIS at EIC lumi = many events
    np.testing.assert_allclose(p2.n_events.sum(), 4 * p1.n_events.sum(),
                               rtol=1e-12)
    g1 = ToyG1()
    o1 = fom.project_observables(cfg, sc1, p1, g1, toy_b1, toy_delta_gluon)
    o2 = fom.project_observables(cfg, sc2, p2, g1, toy_b1, toy_delta_gluon)
    m = p1.accepted & (p1.n_events > 10)
    # delta ~ 1/sqrt(L): 4x lumi -> 2x smaller errors
    np.testing.assert_allclose(o1["err_azz"][m], 2 * o2["err_azz"][m],
                               rtol=1e-9)
    # asymmetries are bounded and finite where accepted
    for key in ("a_par", "azz", "a_cos2phi"):
        vals = o1[key][m]
        assert np.all(np.isfinite(vals))
        assert np.all(np.abs(vals) < 1.0)


def test_acceptance_masks_make_sense():
    cfg = beams.BeamConfig(18.0, beams.LI6, 137.5)
    sc = fom.Scenario()
    p = fom.project_rates(cfg, sc)
    assert p.accepted.any()
    y = p.extras["y"]
    assert np.all(y[p.accepted] <= sc.y_max + 1e-9)
    assert np.all(y[p.accepted] >= sc.y_min - 1e-9)
    w2 = kin.w2(p.x[p.accepted], p.q2[p.accepted])
    assert np.all(w2 >= sc.w2_min - 1e-9)
