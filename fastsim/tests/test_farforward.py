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
