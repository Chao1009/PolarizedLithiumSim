"""plans/09: the near-beam sensor physics behind reports/nanowire_far_forward.

Two claims carry the report and both are pinned here:

  * the deposit ratio is z^2 EXACTLY at fixed velocity, so the hot-spot
    radius goes as z linearly.  This is what lets Argonne's single
    anchor (r_s = 134 nm for a 120 GeV proton -- an EXTRAPOLATED fit
    parameter, not a datum) fix the 6Li and alpha radii without a
    thermodynamic calculation;
  * below w = 2 r_s the firing threshold is zero, the wire fires at any
    bias, and it carries NO charge information -- which is why Argonne's
    ~250 nm MIP-efficiency optimum has zero Z discrimination and Z-ID
    needs deliberately wide wires.

The published anchors themselves (134 nm; the ~250 nm optimum; the
0.8 I_c dark-count wall; the ~1 um hot spot of a 5.5 MeV alpha) are
checked against the numbers the papers state.

Several tests here exist to pin NEGATIVE results, including one that
corrects an earlier version of this study: the ePIC AC-LGAD stack does
NOT beat a threshold device by orders of magnitude on charge information
-- one bit per plane is worth a factor 1.4 against the Neyman-Pearson
optimum.  What actually decides open question #19 against a nanowire is
its GEOMETRIC FILL FACTOR, which is a fabrication number.
"""

import math

import numpy as np
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


def test_the_241Am_alpha_is_a_beta_point_not_a_z_point():
    """The report's first version argued that a relativistic 6Li is an
    'interpolation' between Argonne's two measured points.  It is -- in Q.
    But their 241Am alpha differs from the 120 GeV proton almost entirely
    through 1/beta^2, not through z^2, so it calibrates sqrt(Q) across
    ENERGY and says nothing about z^2 at fixed velocity.  This test pins
    the distinction so the claim cannot quietly come back."""
    r_li = nb.hot_spot_nm(3)
    # their own arithmetic: sqrt(5.5 MeV / 0.1 MeV) x 134 nm ~ 1 um
    r_am_alpha = nb.R_S_PROTON_NM * math.sqrt(5.5 / 0.1)
    assert r_am_alpha == pytest.approx(1000.0, rel=0.05)
    assert nb.R_S_PROTON_NM < r_li < r_am_alpha          # in Q, yes

    # but the alpha's velocity is nothing like the proton's, so the
    # 1/beta^2 factor alone spans more than the whole z^2 range
    beta_am = math.sqrt(2.0 * 5.5 / 3727.4)              # 5.5 MeV alpha
    assert beta_am == pytest.approx(0.054, abs=0.005)
    assert 1.0 / beta_am ** 2 > 3.0 ** 2


def test_the_incumbent_ac_lgad_carries_more_information_per_plane():
    """The finding that decides open question #19 against this
    technology.  A threshold nanowire yields ONE BIT per plane; the
    ePIC AC-LGAD digitises an 8-bit charge over a 30 um active layer,
    which separates 6Li from alpha well enough that a nanowire can at
    best match it.  Landau in silicon, Gaussian-equivalent core width
    FWHM/2.355 = 4.02 xi / 2.355."""
    si = nb.SI_ACLGAD
    # the AC-LGAD active layer is thick enough for Landau to apply here,
    # unlike the 12 nm superconducting film -- and the module refuses to
    # quote an MPV where it does not
    assert si.xi_over_i() > 1.0
    assert nb.NBN.xi_over_i() < 1e-2
    with pytest.raises(ValueError):
        nb.landau_mpv_ev(3, BG_6LI_TOP, nb.NBN)

    mpv = {z: nb.landau_mpv_ev(z, BG_6LI_TOP, si) for z in (1, 2, 3)}
    sg = {z: nb.landau_core_sigma_ev(z, si) for z in (1, 2, 3)}
    assert mpv[1] == pytest.approx(7.1e3, rel=0.05)     # ~1 MIP in 30 um Si
    assert mpv[3] / mpv[1] == pytest.approx(10.5, rel=0.05)
    separation = (mpv[3] - mpv[2]) / math.hypot(sg[3], sg[2])
    assert separation > 4.0                      # against a nanowire's 1 bit

    # even a MERGED alpha + d pair in one 500 um pixel -- the only way the
    # breakup fakes a single hit -- stays separable, because charges add:
    # 4 + 1 = 5 MIP against the 6Li's 9
    xi_pair = nb.landau_xi_ev(2, 1.0, si) + nb.landau_xi_ev(1, 1.0, si)
    sig_pair = nb.LANDAU_FWHM_OVER_XI * xi_pair / 2.355
    mpv_pair = xi_pair * (math.log(2.0 * nb.M_E_EV * BG_6LI_TOP ** 2 / si.i_ev)
                          + math.log(xi_pair / si.i_ev) + 0.200 - 1.0
                          - si.delta)
    assert (mpv[3] - mpv_pair) / math.hypot(sg[3], sig_pair) > 3.0


def test_alpha_deuteron_breakup_is_two_hits_tens_of_pixels_apart():
    """The handle that beats any dE/dx scheme for the background #19 was
    written about.  Both fragments sit at beam rigidity so neither is
    dispersed, but k_rel is TRANSVERSE and therefore unboosted: the alpha
    at 4 p_u and the deuteron at 2 p_u take opposite kicks of the same
    k_rel and land far apart against a 500 um pitch."""
    k_gev, r12_m, pitch_mm = 0.040, 30.6, 0.5
    for p_u, expect_mm in ((137.5, 6.7), (50.0, 18.4), (20.5, 44.8)):
        sep_mm = 1e3 * r12_m * k_gev * (1.0 / (4.0 * p_u) + 1.0 / (2.0 * p_u))
        assert sep_mm == pytest.approx(expect_mm, rel=0.02)
        assert sep_mm / pitch_mm > 13.0          # resolvable at every optics


# --- the Landau itself, and what the 6Li/alpha separation really needs ---

def test_landau_density_reproduces_its_reference_values():
    """The sampler underneath every fake rate below.  If these move, the
    fake rates are meaningless."""
    lam = np.linspace(-1.0, 1.0, 81)
    dens = nb.landau_density(lam)
    assert lam[int(np.argmax(dens))] == pytest.approx(nb.LANDAU_MODE, abs=0.03)
    assert float(nb.landau_density(0.0)[0]) == pytest.approx(0.17805, rel=0.01)
    # exact FWHM of the Landau is 4.02 xi
    i = int(np.argmax(dens))
    wide = np.linspace(-2.0, 6.0, 400)
    d2 = nb.landau_density(wide)
    j = int(np.argmax(d2))
    half = 0.5 * d2[j]
    lo = float(np.interp(half, d2[:j + 1], wide[:j + 1]))
    hi = float(np.interp(-half, -d2[j:], wide[j:]))
    assert hi - lo == pytest.approx(nb.LANDAU_FWHM_OVER_XI, rel=0.02)


def test_landau_sampler_matches_the_density_it_came_from():
    lam = nb.landau_sample(200000, np.random.default_rng(3))
    assert float(np.median(lam)) == pytest.approx(1.3557, abs=0.05)
    # the tail falls as ~1/lambda, nothing like a Gaussian -- which is the
    # whole reason a "sigma" is the wrong figure of merit here
    assert 0.05 < float((lam > 10).mean()) < 0.20
    assert float((lam > 100).mean()) > 1e-3


def test_one_bit_per_plane_is_nearly_as_good_as_the_optimum():
    """The correction to the first version of reports/nanowire_far_forward.
    Against the Neyman-Pearson optimum (per-plane Landau log-likelihood
    ratio), a one-bit threshold with a majority-of-k decision loses a small
    factor, NOT the orders of magnitude an '8 bits versus 1 bit' framing
    implies.  The power comes from coincidence across planes, not from
    precision within one."""
    kw = dict(efficiency=0.95, n_planes=4, n_mc=400000, plane_efficiency=0.99)
    opt, _ = nb.zid_fake_rate(3, 2, readout="llr",
                              rng=np.random.default_rng(7), **kw)
    bit, _ = nb.zid_fake_rate(3, 2, readout="threshold",
                              rng=np.random.default_rng(7), **kw)
    assert opt > 0.0 and bit > 0.0
    assert 1.0 <= bit / opt < 4.0


def test_the_naive_analogue_estimators_lose_to_one_bit():
    """A Landau has no mean, so a plain sum is dragged by a single delta
    ray and a truncated mean is far behind a coincidence.  This is why
    'more bits' does not automatically mean 'better Z-ID'."""
    kw = dict(efficiency=0.95, n_planes=4, n_mc=400000, plane_efficiency=0.99)
    bit, _ = nb.zid_fake_rate(3, 2, readout="threshold",
                              rng=np.random.default_rng(7), **kw)
    trunc, _ = nb.zid_fake_rate(3, 2, readout="trunc",
                                rng=np.random.default_rng(7), **kw)
    tot, _ = nb.zid_fake_rate(3, 2, readout="sum",
                              rng=np.random.default_rng(7), **kw)
    assert trunc > 10.0 * bit
    assert tot > trunc


def test_the_nanowire_loses_on_fill_factor_not_on_bits():
    """Where the two technologies actually part company.  The coincidence
    that makes one bit sufficient needs every plane to record the track; a
    wire comb only does so over its fill factor.  At the ANL device's 50%
    a four-plane array cannot reach 95% 6Li efficiency at all."""
    kw = dict(efficiency=0.95, n_planes=4, n_mc=300000, readout="threshold")
    ok, e_ok = nb.zid_fake_rate(3, 2, plane_efficiency=0.99,
                                rng=np.random.default_rng(7), **kw)
    assert ok == ok and e_ok >= 0.95            # silicon reaches it
    for fill in (0.50, 0.40, 0.25):
        f, reach = nb.zid_fake_rate(3, 2, plane_efficiency=fill,
                                    rng=np.random.default_rng(7), **kw)
        assert math.isnan(f)                    # cannot reach 95% at all
        assert reach < 0.95
    # and the reach falls monotonically with fill
    reaches = [nb.zid_fake_rate(3, 2, plane_efficiency=fl, efficiency=0.999,
                                n_planes=4, n_mc=120000, readout="threshold",
                                rng=np.random.default_rng(7))[1]
               for fl in (0.25, 0.50, 0.99)]
    assert reaches[0] < reaches[1] < reaches[2]
