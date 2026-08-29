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
    # THE PRICED OPTIMUM, on the per-configuration levers that became the
    # default on 2026-08-29.  The pot dispersion enters here only as the
    # RATIO D / R12, which is 1.62e-2 / 1.35e-2 / 9.74e-3 m/m; the single
    # 18 x 275 pair this function used until then understated the two
    # lower configurations by 66% and 39%, and the top row is unmoved
    # because POT_R12 and POT_DISPERSION ARE the 18 x 275 triple's first
    # and third entries.  The single-lever numbers -- 49.7 / 175.6 / 89.3,
    # eps 0.423 / 0.323 / 0.332, 1/7.1 / 1/13.3 / 1/9.5, the envelope
    # 0.33 / 0.17 / 0.12 mrad, which is what every tagging number
    # published before that date was priced with -- are reproduced by
    # levers="18x275" and pinned below.
    expect = ((46.5, 0.374, 6.8, 0.363, 3.80), (164.1, 0.251, 12.8, 0.192, 1.80),
              (89.3, 0.332, 9.5, 0.117, 0.92))
    for cfg, (r_h, eps, one_over_l, env_x, env_y) in zip(beams.default_configs("6Li"), expect):
        t = ff.tagging_optics_point(cfg)
        assert abs(t["r_h"] / r_h - 1.0) < 0.02
        assert abs(t["acceptance"] - eps) < 0.005
        assert abs(1.0 / t["lumi_fraction"] - one_over_l) < 0.15
        o = t["optics"]
        assert o.envelope == pytest.approx((1e-3 * env_x, 1e-3 * env_y), abs=1e-5)
        assert o.lumi_fraction == pytest.approx(t["lumi_fraction"])
    # the retired single-lever behaviour, still reachable and still exact
    old = ((49.7, 0.423, 7.1, 0.328), (175.6, 0.323, 13.3, 0.166),
           (89.3, 0.332, 9.5, 0.117))
    for cfg, (r_h, eps, one_over_l, env_x) in zip(beams.default_configs("6Li"),
                                                  old):
        t = ff.tagging_optics_point(cfg, levers="18x275")
        assert abs(t["r_h"] / r_h - 1.0) < 0.02
        assert abs(t["acceptance"] - eps) < 0.005
        assert abs(1.0 / t["lumi_fraction"] - one_over_l) < 0.15
        assert abs(1e3 * t["env_x"] - env_x) < 0.005
    with pytest.raises(ValueError):
        ff.tagging_optics_point(beams.default_configs("6Li")[0],
                                levers="per-configuration")


def test_alpha_tag_at_the_yellow_report_optics_is_the_off_rigidity_window_only():
    """With the Yellow Report divergence the 6Li alpha's near-beam tail is
    below the 10 sigma envelope at every configuration: what survives is
    the ~1.5% that falls below R = 0.95, and a tagging optics recovers
    23-32% at 1/6.8-1/12.8 of the luminosity."""
    rng = np.random.default_rng(3)
    for cfg in beams.default_configs("6Li"):
        k = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG, cfg.ion_momentum_per_nucleon,
                                        100_000, rng=rng)
        yr = ff.acceptance_summary(k["R"], k["theta"], k["pT"], ff.yr_optics(cfg),
                                   phi=k["phi"])
        tag = ff.acceptance_summary(k["R"], k["theta"], k["pT"], ff.tagging_optics(cfg),
                                    phi=k["phi"])
        # 2026-08-28: the tag is no longer the low-R window alone.  The
        # over-rigid branch of route_charged (plans/09 B1) hands the HIGH-R
        # side of the same k distribution to the pots' inner half, where
        # the measured dispersion puts an R > 1.05 fragment past the
        # central blind block, and that slice is 1.2 points at 18 x 275.
        # It was routed "lost" by construction until the triton scan.
        assert 1.0 - yr["lost"] < 0.045
        assert yr["RP-inner (over-rigid)"] > 0.005
        assert yr["RomanPots"] < 0.02
        assert yr["RP (pT tail, R~1)"] < 0.003
        # 0.324 / 0.229 / 0.292 on the per-configuration levers of
        # 2026-08-29; on the single 18 x 275 pair the band was 0.27-0.40,
        # and the 10 x 100 row is where the flip is felt.
        assert 0.22 < 1.0 - tag["lost"] < 0.35


def test_the_reconstruction_chain_and_the_fast_simulation_share_one_divergence():
    import pathlib
    import sys
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "evgen"))
    reco = pytest.importorskip("polligen.reco")
    for cfg in beams.default_configs("6Li"):
        assert reco.sigma_theta_for(cfg) == ff.sigma_theta_for(cfg)
        assert reco.rp_hole_acceptance(50.0, 0.08, 0.93)["acc"] == \
            ff.hole_acceptance(50.0, 0.08, 0.93)["acc"]


def test_the_per_configuration_pot_levers_are_the_default_and_what_they_cost():
    """plans/09 B1's one open consequence, closed on 2026-08-29 and pinned.

    `tagging_optics_point` reaches the horizontal envelope through the pot
    dispersion, D dp/p / R12.  Both were measured PER CONFIGURATION on
    2026-08-28 and are 1.62e-2 / 1.35e-2 / 9.74e-3 m/m against the single
    18 x 275 pair 0.30 / 30.6 = 9.80e-3 that every tagging number
    published before 2026-08-29 was priced with; the top configuration
    lands within 0.6% of the historical scalar, so only the two lower rows
    move.  They are now the default, because the per-configuration
    transport is the working assumption everywhere (POT_LEVERS), and this
    test says what that costs: eps at 10 x 100 falls 22%, and
    `money_cos2phi_coherent_reco.py` runs on the optimum."""
    expect = ((46.5, 0.374, 6.8), (164.1, 0.251, 12.8), (89.3, 0.332, 9.5))
    for cfg, (r_h, eps, one_over_l) in zip(beams.default_configs("6Li"),
                                           expect):
        t = ff.tagging_optics_point(cfg)
        assert abs(t["r_h"] / r_h - 1.0) < 0.02
        assert abs(t["acceptance"] - eps) < 0.005
        assert abs(1.0 / t["lumi_fraction"] - one_over_l) < 0.15
    # the two lower configurations are where the flip is felt; at 18 x 275
    # the two paths are the same numbers and agree exactly.
    for cfg in beams.default_configs("6Li")[:2]:
        assert (ff.tagging_optics_point(cfg)["r_h"]
                != pytest.approx(
                    ff.tagging_optics_point(cfg, levers="18x275")["r_h"],
                    rel=1e-6))
    top = beams.default_configs("6Li")[2]
    assert ff.tagging_optics_point(top)["r_h"] == pytest.approx(
        ff.tagging_optics_point(top, levers="18x275")["r_h"], rel=1e-12)
    # the opt-in flag is gone, not silently ignored
    import inspect
    params = inspect.signature(ff.tagging_optics_point).parameters
    assert "per_config_levers" not in params
    assert params["levers"].default == "per-config"

