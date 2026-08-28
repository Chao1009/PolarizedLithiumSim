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
