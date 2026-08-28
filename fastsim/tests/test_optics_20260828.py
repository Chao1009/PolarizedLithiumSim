"""The per-configuration far-forward optics of 2026-08-28 (plans/10 A2-A4):
the Yellow Report divergence with the gamma-matched species step, the
rectangular envelope in the spectator routing, and the tagging optics."""

import numpy as np
import pytest

from polli_fastsim import beams
from polli_fastsim import farforward as ff
from polli_fastsim import spectator as sp


def test_yr_optics_reproduces_the_yellow_report_rows():
    lo, mid, top = beams.default_configs("6Li")
    assert [round(1e6 * x) for x in ff.sigma_theta_for(lo)] == [220, 380]
    assert [round(1e6 * x) for x in ff.sigma_theta_for(mid)] == [180, 180]
    assert [round(1e6 * x) for x in ff.sigma_theta_for(top)] == [92, 92]
    # protons: the tables themselves
    for cfg, hv in zip(beams.default_configs("p"), ((220, 380), (180, 180), (65, 65))):
        assert tuple(round(1e6 * x) for x in ff.sigma_theta_for(cfg)) == hv
    o = ff.yr_optics(lo)
    assert o.envelope == pytest.approx((2.2e-3, 3.8e-3), rel=1e-3)
    assert not o.isotropic and ff.yr_optics(mid).isotropic


def test_rectangular_envelope_in_the_routing():
    o = ff.yr_optics(beams.default_configs("6Li")[0])      # 2.2 x 3.8 mrad
    # a beam-rigidity fragment at 3 mrad clears the envelope horizontally
    # but not vertically
    assert ff.route_charged(1.0, 3.0e-3, 0.1, o, phi=0.0) == 4
    assert ff.route_charged(1.0, 3.0e-3, 0.1, o, phi=np.pi / 2) == 0
    # without an azimuth the cut is the horizontal circle
    assert ff.route_charged(1.0, 3.0e-3, 0.1, o) == 4
    assert ff.route_charged(1.0, 2.0e-3, 0.1, o) == 0
    # the legacy isotropic optics is unchanged
    leg = ff.HIGH_ACCEPTANCE
    assert ff.route_charged(1.0, 1.0e-3, 0.05, leg, phi=1.0) == 4
    assert ff.route_charged(1.0, 0.5e-3, 0.30, leg, phi=1.0) == 0


def test_tagging_optics_reproduces_the_priced_optimum():
    expect = ((49.7, 0.422, 7.1, 0.33, 3.80), (175.6, 0.322, 13.3, 0.17, 1.80),
              (89.3, 0.332, 9.5, 0.12, 0.92))
    for cfg, (r_h, eps, one_over_l, env_x, env_y) in zip(beams.default_configs("6Li"), expect):
        t = ff.tagging_optics_point(cfg)
        assert abs(t["r_h"] / r_h - 1.0) < 0.02
        assert abs(t["acceptance"] - eps) < 0.005
        assert abs(1.0 / t["lumi_fraction"] - one_over_l) < 0.15
        o = t["optics"]
        assert o.envelope == pytest.approx((1e-3 * env_x, 1e-3 * env_y), abs=1e-5)
        assert o.lumi_fraction == pytest.approx(t["lumi_fraction"])


def test_alpha_tag_at_the_yellow_report_optics_is_the_off_rigidity_window_only():
    """With the Yellow Report divergence the 6Li alpha's near-beam tail is
    below the 10 sigma envelope at every configuration: what survives is
    the ~1.5% that falls below R = 0.95, and a tagging optics recovers
    ~30% at 1/7-1/13 of the luminosity."""
    rng = np.random.default_rng(3)
    for cfg in beams.default_configs("6Li"):
        k = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG, cfg.ion_momentum_per_nucleon,
                                        100_000, rng=rng)
        yr = ff.acceptance_summary(k["R"], k["theta"], k["pT"], ff.yr_optics(cfg),
                                   phi=k["phi"])
        tag = ff.acceptance_summary(k["R"], k["theta"], k["pT"], ff.tagging_optics(cfg),
                                    phi=k["phi"])
        assert 1.0 - yr["lost"] < 0.02
        assert yr["RP (pT tail, R~1)"] < 0.003
        assert 0.25 < 1.0 - tag["lost"] < 0.40


def test_the_reconstruction_chain_and_the_fast_simulation_share_one_divergence():
    import pathlib
    import sys
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "evgen"))
    reco = pytest.importorskip("polligen.reco")
    for cfg in beams.default_configs("6Li"):
        assert reco.sigma_theta_for(cfg) == ff.sigma_theta_for(cfg)
        assert reco.rp_hole_acceptance(50.0, 0.08, 0.93)["acc"] == \
            ff.hole_acceptance(50.0, 0.08, 0.93)["acc"]
