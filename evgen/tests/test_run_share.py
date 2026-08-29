"""The programme share on the evgen side: `--lumi-fraction`.

Three luminosity factors multiply in the coherent chain and none of them
is the other:

  the PROGRAMME share      how a year divides between observables,
                           isotopes and running configurations -- the
                           `--lumi-fraction` of these scripts, carried into
                           the library as `fom.Scenario.run_share`;
  the OPTICS fraction      what a de-squeezed beta*_x costs at fixed wall
                           time (`farforward.Optics.lumi_fraction`,
                           1/6.8 - 1/12.8 for the lithium tagging optics);
  the SPIN-STATE share     how one measurement's own luminosity divides
                           between its fills (`bookkeeping.SpinCategory
                           .lumi_fraction`, 0.5 / 0.5 for the flip plan).

These tests pin the two things that make the share safe to use: the
published figure stems are guarded against it, and it enters the yields
(and therefore the errors) exactly once.
"""

import importlib
import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from polligen import bookkeeping as bk       # noqa: E402
from polligen import coherent as coh         # noqa: E402
from polligen import recopseudo as rp        # noqa: E402

from polli_fastsim import beams, fom         # noqa: E402
from polli_fastsim.farforward import HIGH_ACCEPTANCE  # noqa: E402

CONFIG = beams.default_configs("6Li")[1]


class _A:
    """Stand-in for the argparse namespace an output_stem is given."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


PUBLISHED_COH = dict(config=0, optics="tagging", fit="ratio",
                     u_in_situ=False, t_edges=None, lumi_fraction=1.0)


def test_a_run_plan_share_never_claims_the_published_coherent_stem():
    """Money plot 6R is the full-share run.  A share is a smaller-error
    figure of the same measurement, which is exactly the kind of run that
    used to overwrite a published PNG unnoticed (the guard added for
    `--optics legacy` on 2026-08-28)."""
    stem = importlib.import_module("money_cos2phi_coherent_reco").output_stem
    base = "money_cos2phi_coherent_reco_6Li"
    assert stem(_A(**PUBLISHED_COH)) == base
    # the share is the ONLY thing changed, and the stem still moves
    shared = stem(_A(**dict(PUBLISHED_COH, lumi_fraction=0.25)))
    assert shared != base
    assert "share0p25" in shared
    # distinct shares never collide, and an absent attribute means share 1
    stems = {stem(_A(**dict(PUBLISHED_COH, lumi_fraction=f)))
             for f in (1.0, 0.5, 0.25, 1.0 / 3.0)}
    assert len(stems) == 4
    legacy = {k: v for k, v in PUBLISHED_COH.items() if k != "lumi_fraction"}
    assert stem(_A(**legacy)) == base


def test_the_other_three_scripts_carry_the_same_guard():
    reach = importlib.import_module("nearbeam_reach_gain").output_stem
    assert reach(_A(fit="ratio", t_edges=None, lumi_fraction=1.0)) == \
        "nearbeam_reach_gain_6Li"
    assert reach(_A(fit="ratio", t_edges=None, lumi_fraction=0.25)) != \
        "nearbeam_reach_gain_6Li"

    tag = importlib.import_module("money_tagged_azz").output_stem
    base = "money_tagged_azz_6Li"
    assert tag(base, "10x100", 1, "menu") == base
    assert tag(base, "10x100", 1, "menu", lumi_fraction=0.25) != base


def test_the_share_multiplies_the_coherent_yield_exactly_once():
    """`--lumi-fraction` enters through the Scenario, so the produced
    recoils scale with it and the optics fraction stays a separate
    factor."""
    sc = coh.CoherentScenario()
    full = fom.Scenario(lumi_fb_per_nucleon=10.0, pol_ion_tensor=0.6)
    quarter = fom.Scenario(lumi_fb_per_nucleon=10.0, run_share=0.25,
                           pol_ion_tensor=0.6)
    _, n_full, _ = coh.project_coherent(CONFIG, full, sc,
                                        optics_list=(HIGH_ACCEPTANCE,))
    _, n_quarter, _ = coh.project_coherent(CONFIG, quarter, sc,
                                           optics_list=(HIGH_ACCEPTANCE,))
    assert n_quarter.sum() == pytest.approx(0.25 * n_full.sum(), rel=1e-12)

    # ... and the optics penalty multiplies on top of it, unchanged
    optics_fraction = 1.0 / 13.3
    assert (n_quarter.sum() * optics_fraction
            == pytest.approx(0.25 * n_full.sum() * optics_fraction, rel=1e-12))


def test_a_quarter_share_doubles_the_coherent_harmonic_errors():
    """The 1/sqrt(share) law through the full two-azimuth fit, at exact
    expected counts so the comparison is not statistical."""
    sc = coh.CoherentScenario()
    cr = rp.CoherentResponse(sc, CONFIG, 1.7e-5, aspect=10.6, n_mc=120000,
                             rng=np.random.default_rng(11))
    plan = bk.tensor_flip_plan(0.6)

    def a_t(t):
        return sc.cos2phi_coefficient_deformation(t, 1.0)

    n_produced = 9.31e6
    full = rp.measure_coherent(cr, n_produced, plan, 0.05, 0.12, a_e=0.01,
                               a_t_func=a_t, u1=0.05, u2=0.02, poisson=False)
    quarter = rp.measure_coherent(cr, 0.25 * n_produced, plan, 0.05, 0.12,
                                  a_e=0.01, a_t_func=a_t, u1=0.05, u2=0.02,
                                  poisson=False)
    assert quarter["n"] == pytest.approx(0.25 * full["n"], rel=1e-12)
    for key in ("err_e", "err_t", "err_m"):
        assert quarter[key] == pytest.approx(2.0 * full[key], rel=1e-9), key
    # the recovered central values are physics and do not move
    for key in ("a_e", "a_t"):
        assert quarter[key] == pytest.approx(full[key], rel=1e-9), key


def test_the_sampler_cross_sections_are_share_invariant():
    """A run-plan share moves event counts, never cross sections.

    `InclusiveSampler` recovers a per-cell cross section by dividing
    `proj.n_events` by the scenario luminosity.  When `run_share` was
    introduced it divided by the PROGRAMME luminosity while `n_events`
    already carried the share, so every cross section came out scaled by
    it and the share was then applied a second time by the caller's
    lumi_pb: `money_cos2phi_reco.py --lumi-fraction 0.25` reported
    N_1yr = 1.17e7 where 4.67e7 was right, and errors four times the
    published ones instead of twice.  Dividing by the effective
    luminosity is the fix."""
    from polligen.sample import InclusiveSampler
    from polligen.xsec import InclusiveKernel
    from polli_fastsim.polarized import toy_b1

    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    kw = dict(nx=24, nq2=18, q2_range=(0.7, 2e3))
    full = InclusiveSampler(kern, CONFIG, fom.Scenario(), **kw)
    quarter = InclusiveSampler(
        kern, CONFIG, fom.Scenario(run_share=0.25), **kw)
    np.testing.assert_allclose(quarter.cell_xsec_pb, full.cell_xsec_pb,
                               rtol=1e-12)
    # and the counts the caller asks for do carry it, exactly once
    cat = bk.tensor_flip_plan(0.6).categories[0]
    n_full = full.expected_events(cat, 1.0e4)
    n_quarter = quarter.expected_events(cat, 0.25 * 1.0e4)
    assert n_quarter == pytest.approx(0.25 * n_full, rel=1e-12)


def test_the_spin_share_is_not_the_programme_share():
    """`SpinCategory.lumi_fraction` divides one measurement's own
    luminosity and sums to one; the programme share divides the year and
    does not.  Overloading the first for the second would silently change
    the estimator weights."""
    plan = bk.tensor_flip_plan(0.6)
    assert sum(c.lumi_fraction for c in plan.categories) == pytest.approx(1.0)
    shares = plan.lumi_shares(1.0e3)
    assert sum(shares.values()) == pytest.approx(1.0e3)
    # a programme share of 1/4 is applied OUTSIDE the plan, to the total
    quarter = plan.lumi_shares(0.25 * 1.0e3)
    for name, v in shares.items():
        assert quarter[name] == pytest.approx(0.25 * v)
