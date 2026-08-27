"""plans/10: what beam divergence a polarized lithium fill would actually have.

Every far-forward acceptance in this repository is exp(-B (10 sigma_theta
A p_u)^2), and `farforward.HIGH_ACCEPTANCE` is a single energy-independent,
isotropic, proton-derived 72.7 urad.  The Yellow Report's own beam tables
(10.1 for e+p, 10.2 for e+Au) are none of those things, and these tests pin
the estimate that follows from them:

  * the tables are being read correctly -- Table 11.48's angular-divergence
    column follows from Table 10.1 times the proton momentum;
  * the light-ion step is KINEMATIC.  At a machine configuration the lattice
    is set by rigidity, so an A/Z = 2 ion sits at half the proton's
    beta*gamma and picks up sqrt(2) of divergence at equal normalised
    emittance;
  * that equal-emittance assumption is calibrated against gold, which has
    ~450x lithium's intrabeam scattering and still costs at most a factor
    2.6 in eps_N -- and nothing at all horizontally.

The estimate is 1.3x the repo's assumption at top energy and 4x at the
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
    """At a machine configuration the lattice is set by RIGIDITY, so an
    A/Z = 2 ion sits at half the proton's beta*gamma and picks up sqrt(2)
    of divergence at equal normalised emittance.  Calibrating that against
    the published gold rows shows gold -- with ~450x lithium's intrabeam
    scattering -- costs at most a factor 2.6 in eps_N and nothing
    horizontally, so equal eps_N is the right assumption for lithium."""
    import math
    from polli_fastsim import farforward as ff
    u, m_p = 0.9315, 0.9383
    assert ff.LIGHT_ION_DIVERGENCE_FACTOR == pytest.approx(
        math.sqrt((275.0 / m_p) / (137.5 / u)), rel=1e-3)

    # the gold calibration, from Tables 10.1 and 10.2
    bg_p, bg_au = 275.0 / m_p, 110.0 / u
    (au_h, au_v), _ = ff.YR_GOLD_DIVERGENCE["110GeV/u"]
    expect = 150.0 * math.sqrt(bg_p / bg_au)          # 290-bunch HD proton
    assert (au_h / expect) ** 2 == pytest.approx(0.85, abs=0.10)
    assert (au_v / expect) ** 2 == pytest.approx(2.6, abs=0.4)

    # IBS goes as Z^4/A^2: lithium is far closer to a proton than to gold
    ibs = lambda z, a: z ** 4 / a ** 2
    assert ibs(3, 6) / ibs(1, 1) == pytest.approx(2.25)
    assert ibs(79, 197) / ibs(3, 6) > 400


def test_the_estimate_is_larger_than_what_the_repo_currently_assumes():
    """The point of plans/10: every far-forward acceptance in the repo uses
    72.7 urad, and the estimate is 1.3x that at top energy and 4x at the
    low-energy configuration the coherent programme calls home."""
    from polli_fastsim import farforward as ff
    ratios = {}
    for cfg in ("18x275", "10x100", "5x41"):
        (h, v), dpp = ff.yr_divergence_for(cfg)
        ratios[cfg] = h / ff.HIGH_ACCEPTANCE.sigma_theta
        assert 6e-4 < dpp < 1.1e-3            # dp/p is species-insensitive
    assert ratios["18x275"] == pytest.approx(1.27, abs=0.05)
    assert ratios["10x100"] == pytest.approx(3.5, abs=0.2)
    assert ratios["5x41"] == pytest.approx(4.3, abs=0.2)
    assert ratios["5x41"] > ratios["10x100"] > ratios["18x275"]
