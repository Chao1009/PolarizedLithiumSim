"""plans/09: the near-beam sensor physics behind reports/nanowire_far_forward.

Two claims carry the report and both are pinned here:

  * the deposit ratio is z^2 EXACTLY at fixed velocity, so the hot-spot
    radius goes as z linearly.  This is what lets Argonne's single
    measured point (r_s = 134 nm for a 120 GeV proton) fix the 6Li and
    alpha radii without a thermodynamic calculation;
  * below w = 2 r_s the firing threshold is zero, the wire fires at any
    bias, and it carries NO charge information -- which is why Argonne's
    ~250 nm MIP-efficiency optimum has zero Z discrimination and Z-ID
    needs deliberately wide wires.

The published anchors themselves (134 nm; the ~250 nm optimum; the
0.8 I_c dark-count wall; the ~1 um hot spot of a 5.5 MeV alpha) are
checked against the numbers the papers state.
"""

import math

import pytest

from polligen import nearbeam as nb


BG_6LI_TOP = 137.5 / 0.9315          # 6Li at 137.5 GeV/u
BG_ANL_PROTON = 120.0 / 0.9383       # ANL's 120 GeV calibration proton


def test_same_velocity_so_deposit_ratio_is_exactly_z_squared():
    """The whole extrapolation to z = 3 rests on the two beams having the
    same velocity, so no velocity correction is needed."""
    assert BG_6LI_TOP == pytest.approx(148, abs=1)
    assert BG_ANL_PROTON == pytest.approx(128, abs=1)

    q1 = nb.bethe_mean_ev(1, BG_6LI_TOP)
    for z in (2, 3):
        assert nb.bethe_mean_ev(z, BG_6LI_TOP) / q1 == pytest.approx(z ** 2)


def test_deposit_is_far_above_the_single_photon_threshold():
    """Detection is not the question -- discrimination is."""
    for z in (1, 2, 3):
        q = nb.bethe_mean_ev(z, BG_6LI_TOP)
        assert q > 10.0 * nb.SNSPD_THRESHOLD_EV


def test_bethe_mean_is_insensitive_to_the_density_effect_in_RATIO():
    """delta enters every species the same way, so the z^2 ratio the
    design argument uses does not depend on a number we took flat."""
    soft = nb.Film(nb.NBN.thickness_nm, nb.NBN.density, nb.NBN.z_over_a,
                   nb.NBN.i_ev, delta=0.0)
    r_ref = (nb.bethe_mean_ev(3, BG_6LI_TOP) / nb.bethe_mean_ev(1, BG_6LI_TOP))
    r_soft = (nb.bethe_mean_ev(3, BG_6LI_TOP, soft)
              / nb.bethe_mean_ev(1, BG_6LI_TOP, soft))
    assert r_ref == pytest.approx(r_soft, rel=1e-12)


def test_landau_formalism_is_out_of_range_for_this_film():
    """The report quotes the Bethe mean and NOT the Landau MPV; this is
    the reason, and it must stay true if the film ever changes."""
    assert nb.NBN.xi_over_i() < 1e-2
    # a Landau-valid absorber would be millimetres of this material
    thick = nb.Film(1e6, nb.NBN.density, nb.NBN.z_over_a, nb.NBN.i_ev)
    assert thick.xi_over_i() > 1.0


def test_hot_spot_radius_is_linear_in_z_and_anchored_on_the_measurement():
    assert nb.hot_spot_nm(1) == pytest.approx(nb.R_S_PROTON_NM)
    assert nb.hot_spot_nm(2) == pytest.approx(268.0)
    assert nb.hot_spot_nm(3) == pytest.approx(402.0)
    # r_s ~ sqrt(Q) with Q ~ z^2 is r_s ~ z: the two routes must agree
    q1, q3 = (nb.bethe_mean_ev(1, BG_6LI_TOP), nb.bethe_mean_ev(3, BG_6LI_TOP))
    assert math.sqrt(q3 / q1) == pytest.approx(3.0)


def test_argonnes_mip_optimum_is_exactly_the_proton_saturation_width():
    """~250 nm is 2 r_s for a proton: maximal MIP efficiency and zero
    charge information.  This coincidence is the design tension."""
    assert nb.saturation_width_nm(1) == pytest.approx(268.0)
    assert nb.saturation_width_nm(1) == pytest.approx(nb.W_MIP_OPTIMAL_NM,
                                                      rel=0.08)
    assert nb.threshold_ratio(nb.hot_spot_nm(1), nb.W_MIP_OPTIMAL_NM) == 0.0


def test_threshold_is_zero_below_saturation_and_rises_above_it():
    for z in (1, 2, 3):
        w_sat = nb.saturation_width_nm(z)
        assert nb.threshold_ratio(nb.hot_spot_nm(z), 0.9 * w_sat) == 0.0
        assert nb.threshold_ratio(nb.hot_spot_nm(z), w_sat) == 0.0
        assert nb.threshold_ratio(nb.hot_spot_nm(z), 2.0 * w_sat) > 0.0
    with pytest.raises(ValueError):
        nb.threshold_ratio(100.0, 0.0)


def test_threshold_orders_the_species_and_is_monotone_in_width():
    """Heavier ion -> bigger hot spot -> fires at LOWER bias."""
    for w in (1000.0, 1500.0, 2000.0):
        t = [nb.threshold_ratio(nb.hot_spot_nm(z), w) for z in (1, 2, 3)]
        assert t[0] > t[1] > t[2]
    for z in (1, 2, 3):
        t = [nb.threshold_ratio(nb.hot_spot_nm(z), w)
             for w in (1000.0, 1500.0, 2000.0)]
        assert t[0] < t[1] < t[2]


def test_the_existing_microwire_separates_all_three_below_the_noise_wall():
    """The report's design point: w = 1 um, two bias points, all three
    turn-ons resolved and every one of them below the 0.8 I_c wall where
    Argonne measure the dark-count rate rising exponentially."""
    w = 1000.0                                  # arXiv:2510.11725
    t_p, t_al, t_li = (nb.threshold_ratio(nb.hot_spot_nm(z), w)
                       for z in (1, 2, 3))
    assert (t_li, t_al, t_p) == pytest.approx((0.196, 0.464, 0.732), abs=5e-3)
    assert t_p < nb.DARK_COUNT_WALL
    # the two bias points, and what each must and must not fire on
    for bias, should_fire in ((0.5 * (t_li + t_al), (3,)),
                              (0.5 * (t_al + t_p), (2, 3))):
        fired = tuple(z for z in (1, 2, 3)
                      if bias > nb.threshold_ratio(nb.hot_spot_nm(z), w))
        assert fired == should_fire
        assert bias < nb.DARK_COUNT_WALL


def test_anl_measured_proton_thresholds_invert_to_a_consistent_radius():
    """Their Fig. 7 points, digitised.  Inverting each one for r_s must
    give a CONSISTENT radius, or the one-parameter model is wrong and the
    step to z = 3 means nothing.  The inverted values cluster at ~113 nm,
    0.85 of the 134 nm Argonne quote -- 134 nm is their extrapolated
    'hard' hot spot, and the 15% offset moves no conclusion, since every
    design number here is a RATIO between species."""
    inverted = [0.5 * w * (1.0 - ratio)
                for w, ratio in ((300.0, 0.215), (400.0, 0.49),
                                 (600.0, 0.62), (800.0, 0.70))]
    mean = sum(inverted) / len(inverted)
    assert mean == pytest.approx(113.0, abs=6.0)
    for r in inverted:                       # consistent point to point
        assert r == pytest.approx(mean, rel=0.15)
    for r in inverted:                       # and consistent with the quote
        assert 0.7 * nb.R_S_PROTON_NM < r < 1.05 * nb.R_S_PROTON_NM


def test_species_ordering_survives_the_softer_inverted_radius():
    """Every conclusion is a ratio, so re-anchoring on the inverted
    113 nm instead of the quoted 134 nm must not change the design
    point: 1 um still separates all three below the noise wall."""
    w = 1000.0
    t = [nb.threshold_ratio(nb.hot_spot_nm(z, r_s_proton=113.0), w)
         for z in (1, 2, 3)]
    assert t[0] > t[1] > t[2] > 0.0
    assert t[0] < nb.DARK_COUNT_WALL


def test_a_relativistic_6Li_brackets_between_argonnes_two_measurements():
    """The Z-ID case is an interpolation: 134 nm (their 120 GeV proton)
    < ~400 nm (a relativistic 6Li) < ~1 um (their 5.5 MeV 241Am alpha,
    by the same sqrt(Q) scaling)."""
    r_li = nb.hot_spot_nm(3)
    # their own arithmetic: sqrt(5.5 MeV / 0.1 MeV) x 134 nm ~ 1 um
    r_am_alpha = nb.R_S_PROTON_NM * math.sqrt(5.5 / 0.1)
    assert r_am_alpha == pytest.approx(1000.0, rel=0.05)
    assert nb.R_S_PROTON_NM < r_li < r_am_alpha
