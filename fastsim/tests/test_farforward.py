"""plans/10: what beam divergence a polarized lithium fill would actually have.

Every far-forward acceptance in this repository is exp(-B (10 sigma_theta
A p_u)^2), and `farforward.HIGH_ACCEPTANCE` is a single energy-independent,
isotropic, proton-derived 72.7 urad.  The Yellow Report's own beam tables
(10.1 for e+p, 10.2 for e+Au) are none of those things, and these tests pin
the estimate that follows from them:

  * the tables are being read correctly -- Table 11.48's angular-divergence
    column follows from Table 10.1 times the proton momentum;
  * the light-ion step is KINEMATIC and applies only where the ring rigidity
    CAPS the ion: a gamma-matched ion has the proton's beta*gamma and pays
    nothing, while 6Li at the top configuration sits at half of it and picks
    up sqrt(2) at equal normalised emittance;
  * that equal-emittance assumption is calibrated against gold, which is far
    more IBS-prone than lithium under either normalisation of the same law
    (446x per particle, 17x at fixed beam current) and still costs at most a
    factor 2.6 in eps_N -- and nothing at all horizontally.

The estimate is 1.3x the repo's single 72.7 urad at top energy and 3x at the
low-energy configuration the coherent programme calls home.
"""

import pytest



def test_yr_proton_divergence_table_reproduces_table_11_48():
    """YR Table 11.48's angular-divergence dpT column is Table 10.1's
    divergence times the proton momentum.  If that closes, the tables are
    being read correctly and the light-ion scaling can stand on them."""
    import math
    from polli_fastsim import farforward as ff
    for cfg, p_gev, tab_11_48 in (("10x100", 100.0, 22.0), ("5x41", 41.0, 14.0)):
        (hd_h, hd_v, _, _), _ = ff.YR_PROTON_DIVERGENCE[cfg]
        rms = math.sqrt(0.5 * (hd_h ** 2 + hd_v ** 2))
        assert 1e-3 * rms * p_gev == pytest.approx(tab_11_48, rel=0.12)


def test_light_ion_divergence_is_kinematic_not_ibs():
    """The species step is a beta*gamma ratio, not a blanket sqrt(A/Z).  At
    a machine configuration the lattice is set by RIGIDITY, so a
    rigidity-capped 6Li at 137.5 GeV/u sits at half the 275 GeV proton's
    beta*gamma and picks up sqrt(2) at equal normalised emittance, while a
    gamma-matched one pays nothing.  Calibrating equal eps_N against the
    published gold rows shows gold costs at most 2.6 in eps_N and nothing
    horizontally, and gold is far more IBS-prone than lithium under either
    normalisation of the same law."""
    import math
    from polli_fastsim import beams, farforward as ff
    u, m_p = 0.9315, 0.9383
    cfgs = beams.default_configs("6Li")
    top_factor = ff.sigma_theta_for(cfgs[2])[0] / (
        1e-6 * ff.YR_PROTON_DIVERGENCE["18x275"][0][2])
    assert top_factor == pytest.approx(math.sqrt((275.0 / m_p) / (137.5 / u)),
                                       rel=2e-3)
    for cfg, key in zip(cfgs[:2], ("5x41", "10x100")):    # gamma-matched
        assert ff.sigma_theta_for(cfg)[0] == pytest.approx(
            1e-6 * ff.YR_PROTON_DIVERGENCE[key][0][2], rel=1e-3)

    # the gold calibration, from Tables 10.1 and 10.2
    bg_p, bg_au = 275.0 / m_p, 110.0 / u
    (au_h, au_v), _ = ff.YR_GOLD_DIVERGENCE["110GeV/u"]
    expect = 150.0 * math.sqrt(bg_p / bg_au)          # 290-bunch HD proton
    assert (au_h / expect) ** 2 == pytest.approx(0.85, abs=0.10)
    assert (au_v / expect) ** 2 == pytest.approx(2.6, abs=0.4)

    # IBS: the per-particle growth rate goes as Z^4/A^2 (gold 446x lithium,
    # lithium 2.25x a proton); at fixed BEAM CURRENT one factor Z comes off,
    # Z^3/A^2 (gold 17x lithium, lithium 0.75x a proton -- the normalisation
    # plans/10 SS10.3 quotes).  Lithium is proton-class either way, which is
    # the only thing the equal-emittance assumption needs.
    per_particle = lambda z, a: z ** 4 / a ** 2
    at_fixed_current = lambda z, a: z ** 3 / a ** 2
    assert per_particle(3, 6) / per_particle(1, 1) == pytest.approx(2.25)
    assert per_particle(79, 197) / per_particle(3, 6) == pytest.approx(446,
                                                                      abs=2)
    assert at_fixed_current(3, 6) / at_fixed_current(1, 1) == pytest.approx(0.75)
    assert at_fixed_current(79, 197) / at_fixed_current(3, 6) == pytest.approx(
        17, abs=1)


def test_the_estimate_is_larger_than_what_the_repo_currently_assumes():
    """The point of plans/10: every far-forward acceptance in the repo used
    a single 72.7 urad, and the per-configuration estimate is 1.26x that at
    top energy and 3.0x at the low-energy configuration the coherent
    programme calls home.  Evaluated through `sigma_theta_for`, which
    applies the species step only where rigidity binds -- the blanket
    sqrt(A/Z) of the retired `yr_divergence_for` gave 4.3x at 5 x 41, where
    the ion is gamma-matched and pays nothing (it was pinned here until
    2026-08-28)."""
    from polli_fastsim import beams, farforward as ff
    ratios = {}
    for cfg in beams.default_configs("6Li"):
        h, v = ff.sigma_theta_for(cfg)
        ratios[ff.yr_config_key(cfg)] = h / ff.HIGH_ACCEPTANCE.sigma_theta
        dpp = 1e-4 * ff.YR_PROTON_DIVERGENCE[ff.yr_config_key(cfg)][1]
        assert 6e-4 < dpp < 1.1e-3            # dp/p is species-insensitive
    assert ratios["18x275"] == pytest.approx(1.26, abs=0.03)
    assert ratios["10x100"] == pytest.approx(2.48, abs=0.05)
    assert ratios["5x41"] == pytest.approx(3.03, abs=0.05)
    assert ratios["5x41"] > ratios["10x100"] > ratios["18x275"]
    assert not hasattr(ff, "yr_divergence_for")
    assert not hasattr(ff, "LIGHT_ION_DIVERGENCE_FACTOR")


# --- the pot-plane transport, measured per configuration (plans/09 B1) -----


def test_pot_levers_resolve_at_every_gamma_matched_configuration():
    """`POT_LEVERS` is keyed by machine configuration and reachable three
    ways: by BeamConfig, by key, and by the bare 6Li per-nucleon momentum
    the callers written before 2026-08-28 pass.  The momenta are the
    GAMMA-MATCHED 40.8 / 99.5 / 137.5 GeV/u -- the rigidity-scaled
    20.5 / 50 that this repository carried until 2026-08-27 are not beam
    energies and must not resolve."""
    import pytest
    from polli_fastsim import beams
    from polli_fastsim import farforward as ff

    expect = {"5x41": (19.24, 4.56, 0.311),
              "10x100": (21.25, 3.35, 0.287),
              "18x275": (29.97, 2.93, 0.292)}
    for pu, key in ((40.8, "5x41"), (99.5, "10x100"), (137.5, "18x275")):
        assert ff.pot_levers_for(pu) == expect[key]
        assert ff.pot_levers_for(key) == expect[key]
    for cfg in beams.default_configs("6Li"):
        assert ff.pot_levers_for(cfg) == expect[ff.yr_config_key(cfg)]
    # the transport is a property of the RING, so 7Li at the same machine
    # setting sees the same levers even at a different per-nucleon momentum
    for c6, c7 in zip(beams.default_configs("6Li"),
                      beams.default_configs("7Li")):
        assert ff.pot_levers_for(c6) == ff.pot_levers_for(c7)
    for stale in (20.5, 50.0, 117.9):
        with pytest.raises(KeyError):
            ff.pot_levers_for(stale)
    # R34 at 5 x 41 was None until 2026-08-29 -- the 29.6 mm per-energy
    # insertion shuts the vertical plane and the ladder had nothing to
    # regress on -- and is now 4.56 m, read off a scratch geometry whose
    # four `offset_*_RP_section` constants are zero (tools/fullsim).  The
    # table carries no None any more, which is what lets
    # `separation_at_pots` stop falling back on R12.
    assert all(t[1] is not None for t in ff.POT_LEVERS.values())
    assert ff.pot_levers_for("5x41")[1] == 4.56
    # the plane is still SHUT: reaching the insertion needs 29.6 mm /
    # 4.56 m = 6.49 mrad, past THETA_RP_MAX where the routing ends.  NOT
    # compared against THETA_RP_OUTER_MEASURED: that is a horizontal edge
    # on a 4.2x longer lever, so the two angles do not compare.
    assert (29.6e-3 / ff.pot_levers_for("5x41")[1] > ff.THETA_RP_MAX)
    # the scalar aliases are the 18 x 275 row, so nothing written before
    # the measurement changes meaning
    assert (ff.POT_R12, ff.POT_R34, ff.POT_DISPERSION) == expect["18x275"]
    # x is 4.2 to 10.2x stiffer than y at the three configurations --
    # which is why the pot aperture is a horizontal slot though the pots
    # insert vertically
    for key in ("5x41", "10x100", "18x275"):
        r12, r34, _d = ff.POT_LEVERS[key]
        assert 4.0 < r12 / r34 < 12.0
    # and the stiffness ratio GROWS with energy, 4.2 / 6.3 / 10.2
    ratios = [ff.POT_LEVERS[k][0] / ff.POT_LEVERS[k][1]
              for k in ("5x41", "10x100", "18x275")]
    assert ratios == sorted(ratios)


def test_the_over_rigid_fragment_is_routed_where_the_scan_put_it():
    """Report 3 Table 6 called "no coverage above R = 1.05" a routing
    assumption; the scan of 2026-08-28 measured it and it was wrong.  An
    R = 1.286 triton is on the Roman-Pot silicon in 60 of 60 events at
    every configuration, at dx = +66 mm on the INNER side of the bend,
    because the pot dispersion carries it past the central blind block
    into the offset-free outer band."""
    import numpy as np
    from polli_fastsim import farforward as ff

    for key in ("5x41", "10x100", "18x275"):
        # the 7Li triton, measured
        assert bool(ff.over_rigid_route(1.2897, 0.0, key))
        # the displacement the arithmetic gives is the one the scan saw,
        # +66 to +73 mm, and only because the second-order term is kept
        r12, _r34, d = ff.POT_LEVERS[key]
        d2 = ff.POT_DISPERSION_2[key]
        x_mm = 1e3 * (d * 0.286 + d2 * 0.286 ** 2)
        assert 60.0 < x_mm < 80.0
        assert 1e3 * d * 0.286 > 82.0          # linear alone is 15-25% high
        # the 6Li 3He + t triton at R = 1.5044 is past the last module and
        # stays lost -- unmeasured, and routed as it always was
        assert not bool(ff.over_rigid_route(1.5044, 0.0, key))

    # routing: code 6, its own label, and nothing else moved
    assert ff.route_charged(1.2897, 1e-3, 0.05) == 6
    assert ff.ROUTE_LABELS[6] == "RP-inner (over-rigid)"
    assert ff.route_charged(1.5044, 1e-3, 0.05) == 0
    assert ff.route_charged(0.857, 1e-3, 0.01) == 1     # unchanged
    assert ff.route_charged(0.50, 1e-3, 0.05) == 2      # unchanged

    # the measured HOLE: a horizontal angle of the opposite sign cancels
    # the dispersion and pushes the triton back into the 16 mm block.  At
    # 18 x 275 the scan loses it between -1.6 and -2.5 mrad and finds it
    # again beyond.
    assert bool(ff.over_rigid_route(1.286, 0.0, "18x275"))
    assert not bool(ff.over_rigid_route(1.286, -2.0e-3, "18x275"))
    assert bool(ff.over_rigid_route(1.286, -3.0e-3, "18x275"))
    # and the outer module edge closes it on the other side
    assert not bool(ff.over_rigid_route(1.286, +3.0e-3, "18x275"))


def test_the_measured_outer_edge_is_not_the_module_arithmetic():
    """The pot acceptance stops between 2.9 and 4.0 mrad, at |dx| = 54 to
    127 mm, because the ion strikes the pipe or the magnet aperture -- not
    at 144 mm / R12, which would be 7.5 / 6.8 / 4.8 mrad, and not at the
    5 mrad the acceptance tables assume.  The DEFAULT must stay at
    THETA_RP_MAX so no published acceptance moves under the measurement
    (test_two_hit.py pins that); this is the value to pass explicitly."""
    from polli_fastsim import beams
    from polli_fastsim import farforward as ff

    assert ff.THETA_RP_OUTER == ff.THETA_RP_MAX
    expect = {"5x41": 2.85e-3, "10x100": 3.85e-3, "18x275": 4.00e-3}
    for key, want in expect.items():
        assert ff.theta_rp_outer_for(key) == want
        assert want < ff.POT_OUTER_HALF_WIDTH / ff.POT_LEVERS[key][0]
        assert want < ff.THETA_RP_MAX
    for cfg in beams.default_configs("6Li"):
        assert ff.theta_rp_outer_for(cfg) == expect[ff.yr_config_key(cfg)]
    assert ff.theta_rp_outer_for(137.5) == expect["18x275"]


def test_acceptance_summary_forwards_the_pot_configuration():
    """The over-rigid branch tests the pot-plane displacement against the
    configuration's own blind block -- 48 / 32 / 16 mm -- and 18 x 275 is
    the most permissive of the three.  `acceptance_summary` carried no
    `pot_config` at all until 2026-08-28, so every per-configuration
    caller silently priced the two lower configurations at the 18 x 275
    block: the single-number-everywhere error this measurement exists to
    remove.  The parameter must exist, must be forwarded, and must matter.
    """
    import numpy as np
    from polli_fastsim import beams
    from polli_fastsim import farforward as ff
    from polli_fastsim import spectator as sp

    for cfg, key in zip(beams.default_configs("6Li"),
                        ("5x41", "10x100", "18x275")):
        k = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG,
                                        cfg.ion_momentum_per_nucleon, 100000,
                                        beta=0.30, rng=np.random.default_rng(7))
        cut = ff.Optics("cut", 2.5e-4, 10.0, 2.5e-4)
        wrong = ff.acceptance_summary(k["R"], k["theta"], k["pT"], cut,
                                      phi=k["phi"], pot_config="18x275")
        right = ff.acceptance_summary(k["R"], k["theta"], k["pT"], cut,
                                      phi=k["phi"], pot_config=key)
        label = ff.ROUTE_LABELS[6]
        if key == "18x275":
            assert right[label] == wrong[label]
        else:
            # the 48 and 32 mm blocks swallow the tail the 16 mm one passes
            assert right[label] < 0.25 * wrong[label]
        # and it is the same routing route_charged does
        r = ff.route_charged(k["R"], k["theta"], k["pT"], cut, phi=k["phi"],
                             pot_config=key)
        assert right[label] == pytest.approx(float(np.mean(r == 6)))


def test_the_second_row_of_the_pot_transfer_is_measured_and_symplectic():
    """plans/03 2.2 (2), measured 2026-09-15.

    Until this scan only the FIRST row of the IP6 -> Roman-Pot transfer
    existed: R12, R34 and D are all `d(position at the pot)/d(something
    at the IP)`.  `POT_SECOND_ROW` adds R11, R21, R22, D' and
    `POT_SECOND_ROW_VERTICAL` adds R33, R43, R44, from the
    `tools/fullsim/ff_transfer_scan.py` ladders through a zero-insertion
    geometry.

    The check that this is a MEASUREMENT and not four numbers is that the
    2x2 blocks come out symplectic without having been fitted to be: a
    linear transfer at fixed rigidity has determinant exactly 1, the
    horizontal block is built out of R11/R21/R22 measured on the new legs
    and R12 read off `POT_LEVERS` untouched, and the two halves close on
    each other to 4 %.  Nothing above them moves -- the published triples
    are asserted here bit-for-bit, because the whole addition is only
    safe if they do not.
    """
    from polli_fastsim import farforward as ff

    # (i) APPEND-ONLY.  The published levers are what they were.
    assert ff.POT_LEVERS == {"5x41": (19.24, 4.56, 0.311),
                             "10x100": (21.25, 3.35, 0.287),
                             "18x275": (29.97, 2.93, 0.292)}
    assert ff.POT_DISPERSION_2 == {"5x41": -0.190, "10x100": -0.206,
                                   "18x275": -0.215}
    assert (ff.POT_R12, ff.POT_R34, ff.POT_DISPERSION) == (29.97, 2.93, 0.292)

    keys = ("5x41", "10x100", "18x275")
    assert set(ff.POT_SECOND_ROW) == set(keys)
    assert set(ff.POT_SECOND_ROW_VERTICAL) == set(keys)

    for key in keys:
        m = ff.pot_transfer_for(key)
        # the first row IS POT_LEVERS, not a re-measurement of it
        assert (m["R12"], m["R34"], m["D"]) == ff.POT_LEVERS[key]
        assert m["D2"] == ff.POT_DISPERSION_2[key]
        # (ii) SYMPLECTIC: |R| = 1 for a linear transfer at fixed rigidity
        det_x = m["R11"] * m["R22"] - m["R12"] * m["R21"]
        det_y = m["R33"] * m["R44"] - m["R34"] * m["R43"]
        assert abs(det_x - 1.0) < 0.05, (key, det_x)
        assert abs(det_y - 1.0) < 0.05, (key, det_y)

    # (iii) the two physics facts the scan turned up, which are the
    # reason the second row is not a formality.
    # R22 CHANGES SIGN between 10 x 100 and 18 x 275: the horizontal
    # phase advance to the pots crosses a node, so only at the top
    # configuration is a positive IP angle still positive at the pot.
    r22 = [ff.POT_SECOND_ROW[k][2] for k in keys]
    assert r22[0] < 0 and r22[1] < 0 and r22[2] > 0
    # R44 at 5 x 41 is ZERO to the measurement (0.0055 +- 0.0040): the
    # vertical plane is point-to-parallel there and the outgoing angle
    # carries no memory of the IP angle.  It is an order of magnitude
    # below the other two.
    r44 = [abs(ff.POT_SECOND_ROW_VERTICAL[k][2]) for k in keys]
    assert r44[0] < 0.02
    assert r44[0] < 0.2 * r44[1] < r44[2]
    # D' is FLAT across a factor 6.7 in beam energy, where every other
    # element moves by 2-4x
    dp = [ff.POT_SECOND_ROW[k][3] for k in keys]
    assert max(dp) / min(dp) < 1.10
    # R11 and R21 move the OTHER way and monotonically: the point-to-point
    # magnification grows with energy (1.15 -> 1.23 -> 1.85) while the
    # angle a millimetre at the IP buys falls by 4x (84 -> 65 -> 21
    # microrad/mm), which is the same stiffening R12/R34 shows
    r11 = [ff.POT_SECOND_ROW[k][0] for k in keys]
    r21 = [abs(ff.POT_SECOND_ROW[k][1]) for k in keys]
    assert r11 == sorted(r11)
    assert r21 == sorted(r21, reverse=True)
    assert r21[0] / r21[2] > 3.0


def test_the_pot_transfer_resolves_and_propagates_like_the_levers():
    """`pot_transfer_for` takes what `pot_levers_for` takes -- a
    BeamConfig of any species, a key, or a bare gamma-matched 6Li GeV/u
    -- and `propagate_to_pot` is the first far-forward transport in this
    repository that carries an ANGLE.  Its position row must reduce to
    the published levers exactly, which is what keeps it from becoming a
    second, disagreeing copy of `separation_at_pots`."""
    import pytest
    from polli_fastsim import beams
    from polli_fastsim import farforward as ff

    for pu, key in ((40.8, "5x41"), (99.5, "10x100"), (137.5, "18x275")):
        assert ff.pot_transfer_for(pu) == ff.pot_transfer_for(key)
    for cfg in beams.default_configs("6Li"):
        assert ff.pot_transfer_for(cfg) == \
            ff.pot_transfer_for(ff.yr_config_key(cfg))
    for c6, c7 in zip(beams.default_configs("6Li"),
                      beams.default_configs("7Li")):
        assert ff.pot_transfer_for(c6) == ff.pot_transfer_for(c7)
    # A BARE FLOAT is read as a 6Li GeV/u and nothing else, exactly as in
    # `pot_levers_for`: 20.5 and 50 are the retired rigidity-scaled fill
    # energies, and 117.9 is live but is 7Li's top setting -- pass the
    # BeamConfig, not the number, and it resolves (asserted above).
    for off_menu in (20.5, 50.0, 117.9):
        with pytest.raises(KeyError):
            ff.pot_transfer_for(off_menu)
        with pytest.raises(KeyError):
            ff.pot_levers_for(off_menu)

    for key in ("5x41", "10x100", "18x275"):
        r12, r34, d = ff.POT_LEVERS[key]
        r11, r21, r22, dp = ff.POT_SECOND_ROW[key]
        r33, r43, r44 = ff.POT_SECOND_ROW_VERTICAL[key]
        # one unit at a time: each column of the matrix, read back
        assert ff.propagate_to_pot(key, xp=1e-3) == \
            pytest.approx((r12 * 1e-3, r22 * 1e-3, 0.0, 0.0))
        assert ff.propagate_to_pot(key, x=1e-3) == \
            pytest.approx((r11 * 1e-3, r21 * 1e-3, 0.0, 0.0))
        assert ff.propagate_to_pot(key, yp=1e-3) == \
            pytest.approx((0.0, 0.0, r34 * 1e-3, r44 * 1e-3))
        assert ff.propagate_to_pot(key, y=1e-3) == \
            pytest.approx((0.0, 0.0, r33 * 1e-3, r43 * 1e-3))
        assert ff.propagate_to_pot(key, delta=0.01) == \
            pytest.approx((d * 0.01, dp * 0.01, 0.0, 0.0))
        # and the second-order dispersion is the SAME quadratic
        # `over_rigid_route` uses, not a second fit
        x2 = ff.propagate_to_pot(key, delta=0.2, second_order=True)[0]
        assert x2 == pytest.approx(d * 0.2 + ff.POT_DISPERSION_2[key] * 0.04)

    # the displacement row is not a rounding of the angle row: at
    # 18 x 275 one millimetre at the IP is 1.85 mm at the pot and
    # -21 microrad of outgoing angle, and the pot lever R12 turns that
    # angle into nothing the position row already said
    x, xp, _y, _yp = ff.propagate_to_pot("18x275", x=1e-3)
    assert x == pytest.approx(1.852e-3)
    assert xp == pytest.approx(-2.09e-5)


def test_the_ip_to_b0_transport_is_a_pure_drift():
    """The 6.0-20.0 mrad leg of the 2026-09-15 scan is the first ladder
    in this repository to reach the B0 window at all, and what it found
    is that IP6 -> B0 layer 1 is a DRIFT: dx/dtheta_x = 5.900 m and
    dy/dtheta_y = 5.901 m at all three ring settings, to 0.02 %, and the
    layer sits at z = 5896 mm.  A B0 hit therefore measures the IP angle
    directly, with no per-configuration optics between -- which is why
    `B0_DRIFT_M` is one number where `POT_LEVERS` needs three."""
    import pytest
    from polli_fastsim import farforward as ff

    assert ff.B0_DRIFT_M == pytest.approx(ff.B0_LAYER_Z_MM[0] * 1e-3, rel=1e-3)
    assert len(ff.B0_LAYER_Z_MM) == 4
    assert ff.B0_LAYER_Z_MM == tuple(sorted(ff.B0_LAYER_Z_MM))
    # the B0 window is where the routing already puts it, and the drift
    # maps that window onto the tracker's own transverse size
    assert 5.5e-3 == ff.THETA_B0_MIN and 20.0e-3 == ff.THETA_B0_MAX
    assert 30.0 < ff.B0_DRIFT_M * ff.THETA_B0_MIN * 1e3 < 35.0
    assert 115.0 < ff.B0_DRIFT_M * ff.THETA_B0_MAX * 1e3 < 120.0
