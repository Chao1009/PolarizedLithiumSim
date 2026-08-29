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
    written about.  Both fragments sit within 0.5% of beam rigidity, but
    k_rel is TRANSVERSE and therefore unboosted: the alpha at 4 p_u and the
    deuteron at 2 p_u take opposite kicks of the same k_rel and land far
    apart against a 500 um pitch.

    This is the back-of-envelope version, kept because it is the one line
    of arithmetic behind the claim: a single k = 40 MeV/c, the angular
    lever alone, no dispersion.  It is a LOWER bound and it is NOT what the
    documents quote.  The measurement draws both fragments from one
    relative momentum at kappa = 60.7 MeV/c and adds the pot dispersion,
    giving medians of 10.9 / 10.7 / 25.8 mm at 18 x 275 / 10 x 100 / 5 x 41,
    an ordering that stopped being monotone in the beam energy when R12
    was measured per configuration (fastsim/tests/test_two_hit.py,
    plans/09 B4).  Until 2026-08-28 this
    test was evaluated at the retired rigidity-scaled 50 and 20.5 GeV/u,
    and plans/04 #19 quoted its output as the answer."""
    k_gev, r12_m, pitch_mm = 0.040, 30.6, 0.5
    for p_u, expect_mm in ((137.5, 6.7), (99.5, 9.2), (40.8, 22.5)):
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


# --- the aperture lookup must serve both isotopes (plans/09 B3) ------------


def test_rp_aperture_resolves_for_both_isotopes_at_every_configuration():
    """`reco.rp_aperture_for` used to key on a MOMENTUM matched against the
    6Li configurations, so it returned None for 7Li at the top energy
    (117.9 GeV/u against 6Li's 137.5) and the near-beam scans could not be
    run for 7Li at all.  It is keyed on the CONFIGURATION now -- the pot
    slot is a property of the ring, not of the beam in it -- and the
    momentum path is kept for the callers written before 2026-08-28."""
    from polligen import reco
    from polli_fastsim import beams

    # re-measured 2026-08-28 in the current epic-main (plans/09 B1,
    # tools/fullsim): band threshold / measured lever, 48/19.24,
    # 32/21.25, 16/29.97 horizontally and 7.1/3.35, 2.7/2.93 vertically,
    # each reproduced by the first accepted angle of a 0.05 mrad ladder to
    # 1-4%.  The September-2024 triple (2.0, 3.0), (1.35, 3.0), (1.03,
    # 2.3) is kept reachable as reco.RP_APERTURE_SEP2024.
    expect = ((2.50e-3, 8.84e-3), (1.51e-3, 2.12e-3), (0.53e-3, 0.92e-3))
    for iso in ("6Li", "7Li"):
        got = [reco.rp_aperture_for(c) for c in beams.default_configs(iso)]
        assert got == list(expect), iso
    # the two isotopes' top configurations are DIFFERENT beams at the same
    # machine setting, which is exactly what the momentum path could not see
    top6, top7 = (beams.default_configs(i)[2] for i in ("6Li", "7Li"))
    assert top6.ion_momentum_per_nucleon != top7.ion_momentum_per_nucleon
    assert reco.rp_aperture_for(top6) == reco.rp_aperture_for(top7)
    # backward compatibility: 6Li momenta still resolve, everything else
    # returns None rather than interpolating
    for cfg, want in zip(beams.default_configs("6Li"), expect):
        assert reco.rp_aperture_for(cfg.ion_momentum_per_nucleon) == want
    assert reco.rp_aperture_for(20.5) is None          # the retired energy
    assert reco.rp_aperture_for(117.9) is None         # 7Li, momentum path

    # the superseded table is still reachable for reproducing a number
    # published against it, and is NOT what the default returns
    for cfg, want in zip(beams.default_configs("6Li"),
                         ((2.0e-3, 3.0e-3), (1.35e-3, 3.0e-3),
                          (1.03e-3, 2.3e-3))):
        assert reco.rp_aperture_for(
            cfg, table=reco.RP_APERTURE_SEP2024) == want
    # every half-width moved, and at 18 x 275 both shrank -- the 16 mm
    # module and the 2.7 mm insertion replacing 32 mm and 7 mm
    assert reco.RP_APERTURE_MEASURED["18x275"][0] < 0.6 * 1.03e-3
    assert reco.RP_APERTURE_MEASURED["18x275"][1] < 0.5 * 2.3e-3
    # and the cutout is TALLER THAN IT IS WIDE at every configuration,
    # which is the sign statement rp_measure's docstring rests on
    for cx, cy in reco.RP_APERTURE_MEASURED.values():
        assert cy > cx


def test_the_5x41_vertical_plane_is_shut_not_a_3mrad_edge():
    """The 5 x 41 c_y encodes a SHUT plane, and it has to, because
    `rp_measure` accepts on |theta_x| > c_x OR |theta_y| > c_y: a c_y that
    is too small over-accepts.

    The scan gave three Roman-pot rows in the whole 0.20-6.00 mrad vertical
    ladder at phi = 90 and none at phi = 270, and the only one above the
    29.6 mm module edge -- theta = 3.35 mrad, y = +30.94 mm -- is a single
    hit in a single plane with both its 0.05 mrad neighbours empty.  That
    is an isolated grazing hit, not an edge.  c_y is instead the smallest
    angle at which any 5 x 41 silicon could be reached at all: 29.6 mm over
    the largest vertical lever measured anywhere (R34 = 3.35 m at
    10 x 100)."""
    import numpy as np

    from polligen import reco
    from polli_fastsim import beams
    from polli_fastsim import farforward as ff
    from polli_fastsim import spectator as sp

    cx, cy = reco.RP_APERTURE_MEASURED["5x41"]
    assert cy > 8.0e-3                       # shut, not 3.35 mrad
    assert cy == pytest.approx(29.6e-3 / ff.POT_LEVERS["10x100"][1], rel=0.01)
    assert np.isfinite(cy)                   # finite, so callers can print it
    # nothing this programme produces reaches it: no 6Li alpha spectator at
    # 5 x 41 is tagged on the vertical alone
    cfg = beams.default_configs("6Li")[0]
    k = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG,
                                    cfg.ion_momentum_per_nucleon, 200000,
                                    beta=0.30, rng=np.random.default_rng(7))
    thx = np.abs(k["theta"] * np.cos(k["phi"]))
    thy = np.abs(k["theta"] * np.sin(k["phi"]))
    assert not np.any((thy > cy) & (thx <= cx))
    # the retired 3.35 mrad row WOULD have admitted some of them
    assert np.any((thy > 3.35e-3) & (thx <= cx))


def test_the_light_ion_lattice_alternative_is_reachable_and_priced():
    """The 5 x 41 aperture depends on WHICH compact file the ring is run
    from, by a factor 1.55, and the alternative has to be reachable in code
    rather than only in prose (`beamline_5x41_He4.xml`, the Z/A = 0.5
    lattice, against the 41 GeV proton lattice the ePIC baseline
    includes)."""
    from polligen import reco
    from polli_fastsim import farforward as ff

    base_r12 = ff.POT_LEVERS["5x41"][0]
    alt_r12 = ff.POT_LEVERS_LIGHT_ION_LATTICE["5x41"][0]
    assert alt_r12 == pytest.approx(29.81, abs=0.05)
    assert alt_r12 / base_r12 == pytest.approx(1.55, abs=0.02)
    # the vertical is shut in both, and 10 x 100 could not be checked
    assert ff.POT_LEVERS_LIGHT_ION_LATTICE["5x41"][1] is None
    assert ff.POT_LEVERS_LIGHT_ION_LATTICE["10x100"] is None
    # the He4 file is the Yellow Report 275 GeV magnet set scaled, so it
    # reproduces the 18 x 275 lattice to 1% -- and the 18 x 275 row of the
    # alternative IS the baseline row
    assert ff.POT_LEVERS_LIGHT_ION_LATTICE["18x275"] == ff.POT_LEVERS["18x275"]
    assert abs(alt_r12 / ff.POT_LEVERS["18x275"][0] - 1.0) < 0.01

    alt = reco.RP_APERTURE_MEASURED_LIGHT_ION_LATTICE
    assert alt["5x41"][0] == pytest.approx(48e-3 / alt_r12, rel=0.01)
    assert (reco.RP_APERTURE_MEASURED["5x41"][0] / alt["5x41"][0]
            == pytest.approx(1.55, abs=0.03))
    assert alt["10x100"] is None
    assert alt["18x275"] == reco.RP_APERTURE_MEASURED["18x275"]
    # and it is reachable through the lookup, like the September-2024 table
    from polli_fastsim import beams
    cfg = beams.default_configs("6Li")[0]
    assert reco.rp_aperture_for(cfg, table=alt) == alt["5x41"]


def test_li7_alpha_tag_is_optics_blind_and_the_tagging_optics_is_a_net_loss():
    """The B3 result in one assertion.  The 7Li alpha sits at rigidity
    ratio 0.856, inside the Roman-Pot momentum window, so it never has to
    clear the near-beam envelope: its tag is 0.96-0.97 at the Yellow Report
    optics and 0.98-0.99 at the tagging optics -- x1.02 -- bought at 1/8 to
    1/15 of the luminosity.  For 7Li the tagging optics is therefore a
    factor 8-15 NET LOSS, the exact inverse of 6Li, where the same optics
    turns a 1.7-2.6% tag into a 28-35% one.  6Li and 7Li want different machine
    optics and are different runs."""
    import numpy as np
    from polli_fastsim import beams, farforward as ff, spectator as sp

    for i, cfg in enumerate(beams.default_configs("7Li")):
        k = sp.spectator_lab_kinematics(sp.LI7_ALPHA_TAG,
                                        cfg.ion_momentum_per_nucleon,
                                        200_000, beta=0.30,
                                        rng=np.random.default_rng(7))
        assert np.median(k["R"]) == pytest.approx(0.856, abs=0.01)
        tag = {}
        for label, o in (("yr", ff.yr_optics(cfg, "high-acceptance")),
                         ("top", ff.tagging_optics(cfg))):
            s = ff.acceptance_summary(k["R"], k["theta"], k["pT"], o,
                                      phi=k["phi"])
            tag[label] = 1.0 - s["lost"]
            # the momentum window carries it, not the near-beam tail
            assert s["RomanPots"] > 0.94
        assert 0.96 < tag["yr"] < 0.98
        assert 0.98 < tag["top"] < 1.00
        gain = tag["top"] / tag["yr"]
        loss = 1.0 / ff.tagging_optics(cfg).lumi_fraction
        assert 1.00 < gain < 1.04                  # acceptance bought
        assert 8.0 < loss < 16.0                   # luminosity paid
        assert gain / loss < 0.15, i               # a net loss, decisively


# --- the 7Li breakup table and the over-rigid route (plans/09 B1, B3) ------


def test_the_breakup_tables_are_beam_generic_and_7li_has_one():
    """plans/09 B3 held the 7Li table back deliberately: "a routing table
    with no amplitude behind it is the half-built channel this section
    refuses ... the far-forward stream is measuring the triton's
    destination in this run".  It has been measured (B1), so the table
    exists, and `fragment_route_label` / `veto_table` take the beam.

    Before this, `fragment_route_label` hard-coded 6Li through
    `fragment_rigidity`'s defaults and answered the wrong question on a
    7Li fragment: it called an INTACT 7Li "lost (over-rigid)" on
    R = m(7Li)/m(6Li) = 1.166 and the 7Li alpha "RP pT-tail only", which
    is the 6Li alpha's destination."""
    from polligen import coherent as coh

    li7 = coh.veto_table(7, 3)
    assert set(li7) == {"alpha+t", "6Li+n", "alpha+d+n", "6He+p"}
    frag = {name: {f: (r, d) for f, r, d in rows}
            for name, rows in li7.items()}

    # the tag channel: the alpha lands mid-window and the triton, measured
    # 2026-08-28, lands on the pots' INNER half rather than nowhere
    r_a, dest_a = frag["alpha+t"]["alpha"]
    assert r_a == pytest.approx(0.85571, abs=5e-6) and dest_a == "RomanPots"
    r_t, dest_t = frag["alpha+t"]["t"]
    assert r_t == pytest.approx(1.28971, abs=5e-6)
    assert dest_t == "RP-inner (over-rigid)"
    assert frag["6Li+n"]["n"][1] == "ZDC"
    assert frag["6Li+n"]["6Li"][1] == "RomanPots"
    # 6He is not in spectator.NUCLEUS_MASS and falls back on A * M_U,
    # 0.3% low; the destination does not notice
    r_he6 = frag["6He+p"]["6He"][0]
    assert r_he6 == pytest.approx(1.283, abs=1e-3)
    assert frag["6He+p"]["6He"][1] == "RP-inner (over-rigid)"

    # separation energies from AME2020 through spectator.NUCLEUS_MASS
    thr = {name: t for name, t, _f in coh.LI7_BREAKUP}
    assert thr["alpha+t"] == pytest.approx(2.468, abs=0.002)   # plans/01
    assert thr["6Li+n"] == pytest.approx(7.251, abs=0.002)
    assert thr["alpha+d+n"] == pytest.approx(8.725, abs=0.002)
    assert thr["6He+p"] == pytest.approx(9.975, abs=0.002)
    # 7Li is 3p + 4n: there is no alpha + p + n channel, that is 6Li
    assert "alpha+p+n" not in thr
    assert "alpha+p+n" in dict((n, t) for n, t, _f in coh.LI6_BREAKUP)

    # the 6Li table is unchanged, and its 3He + t triton at R = 1.504 is
    # still lost: 152 mm at the pot plane, past the last module at 144
    li6 = coh.veto_table()
    assert li6 == coh.veto_table(6, 3)
    assert li6["3He+t"][1][2] == "lost (over-rigid)"

    # and an intact 7Li is at beam rigidity when the beam IS 7Li
    assert coh.fragment_rigidity(7, 3, beam_a=7, beam_z=3) == 1.0
    assert "beam-blind" in coh.fragment_route_label(7, 3, beam_a=7, beam_z=3)
    # on the 6Li defaults the same call answers the wrong question: an
    # intact 7Li reads R = m(7Li)/m(6Li) = 1.166 and is routed over-rigid
    assert "over-rigid" in coh.fragment_route_label(7, 3)

    with pytest.raises(KeyError):
        coh.breakup_table(4, 2)



def test_the_reach_gain_min_count_guard_drops_the_empty_bin_and_keeps_a_full_one():
    """`nearbeam_reach_gain.py` used to fit every |t| bin it was given.

    At 18 x 275 the measured silicon aperture tags 268 recoils a year and
    231 of them land in one bin, whose fit returned a_t = -1.56 +- 2.22
    against an injected +0.42 -- a point with the shape of a measurement
    and none of the content, and the only thing setting the vertical
    range of panel (a).  Since 2026-08-28 a bin is fitted only if it
    expects at least `MIN_TAGGED_PER_BIN` tagged recoils at one year, the
    expectation being the bin's share of the equally weighted accepted
    response recoils times the row's N_tag.  This pins both directions of
    that threshold and the arithmetic behind it.
    """
    import importlib
    import pathlib as _pl
    import sys as _sys
    _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]
                            / "scripts"))
    reach = importlib.import_module("nearbeam_reach_gain")

    assert reach.MIN_TAGGED_PER_BIN > 231.0     # the bin it exists to retire

    # nine accepted recoils, one of them in 0.17-0.25 and eight in
    # 0.05-0.08, against the 268 tagged/yr of the silicon row at 18 x 275
    t_reco = np.array([0.06] * 8 + [0.20])
    sparse = reach.expected_tagged_in_bin(268.0, t_reco, 0.17, 0.25)
    full = reach.expected_tagged_in_bin(268.0, t_reco, 0.05, 0.08)
    assert sparse == pytest.approx(268.0 / 9.0)
    assert full == pytest.approx(268.0 * 8.0 / 9.0)
    assert sparse < reach.MIN_TAGGED_PER_BIN     # dropped
    # the same shape at the tagging optics' 6.0e6 tagged/yr keeps both
    rich_sparse = reach.expected_tagged_in_bin(6.0e6, t_reco, 0.17, 0.25)
    assert rich_sparse >= reach.MIN_TAGGED_PER_BIN
    assert full * (6.0e6 / 268.0) >= reach.MIN_TAGGED_PER_BIN

    # half-open bins: the upper edge belongs to the next bin, and an
    # empty response expects nothing rather than dividing by zero
    edge = reach.expected_tagged_in_bin(1.0e6, np.array([0.08]), 0.05, 0.08)
    assert edge == 0.0
    assert reach.expected_tagged_in_bin(1.0e6, np.array([]), 0.05, 0.08) == 0.0
