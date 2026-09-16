"""The two switches run 19 added to the truth-level money drivers.

Both are guarded the same way: they are OFF by default, the default path
is bit for bit what it has always been, and a run that turns one on
writes its own PNG stem so that it cannot overwrite a published figure.

  --pzz-plus/--pzz-zero   the TWO-FILL spin-state-sorted path (plans/07
                          WP3).  The published truth-level numbers come
                          from one fill at P_zz = +0.6 read with
                          `estimators.cos2phi_fit_binned`; the two-fill
                          path splits the same luminosity between an
                          m = +-1-rich and an m = 0-rich fill and reads
                          them with `reco.harmonic_ratio_fit`, which is
                          what the reconstructed-level chain already does.
                          The gain it buys is the factor 0.67 of WP3, and
                          these tests pin that factor from three
                          directions: the closed form, the analytic error
                          of `reco.err_harmonic_ratio`, and the fitted
                          error of a pseudo-experiment.

  --binning yr            the Yellow-Report five-bins-per-decade analysis
                          grid (plans/02 Step 1.1 item 3) in place of this
                          project's 40 x 30 log grid.  The lattice is
                          anchored on the decade boundaries, so the tests
                          pin the density at exactly 5 per decade and pin
                          the written-out edge tuples of `beams` against
                          the grid `kinematics.log_grid` actually builds.
"""

import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from polligen import reco                       # noqa: E402

from polli_fastsim import beams                 # noqa: E402
from polli_fastsim import kinematics as kin     # noqa: E402

import money_cos2phi as m5                      # noqa: E402

PUBLISHED_PZZ = 0.60
PUBLISHED_PAIR = (0.60, -1.20)


class _A:
    """Stand-in for the argparse namespace a stem helper is given."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


# --- the stem guards ------------------------------------------------------

def test_the_published_stem_is_the_bare_one():
    """No switch on -> no key, i.e. `money_cos2phi_6Li.png`."""
    args = _A(pdf="toy", binning="log40x30", pzz_plus=None, pzz_zero=None,
              tensor_gamma=False, subtract_tensor_leakage="none",
              b3_frac=0.0, b4_frac=0.0)
    assert m5.output_stem_tag(args) == ""


@pytest.mark.parametrize("kw, key", [
    (dict(pzz_plus=0.6, pzz_zero=-1.2), "twofill"),
    (dict(binning="yr"), "yr"),
    (dict(pdf="grid"), "grid"),
])
def test_every_switch_appends_its_own_key(kw, key):
    base = dict(pdf="toy", binning="log40x30", pzz_plus=None,
                pzz_zero=None, tensor_gamma=False,
                subtract_tensor_leakage="none", b3_frac=0.0, b4_frac=0.0)
    base.update(kw)
    assert key in m5.output_stem_tag(_A(**base)).split("_")


def test_the_keys_compose_in_a_fixed_order():
    """A run with several switches must land on one predictable name, not
    on a permutation of the same keys."""
    args = _A(pdf="grid", binning="yr", pzz_plus=0.6, pzz_zero=-1.2,
              tensor_gamma=True, subtract_tensor_leakage="none",
              b3_frac=0.0, b4_frac=0.0)
    assert m5.output_stem_tag(args) == "grid_yr_twofill_tgamma"


# --- the two-fill estimator ----------------------------------------------

def test_two_fill_pzz_is_none_unless_both_flags_are_given():
    assert m5.two_fill_pzz(_A(pzz_plus=None, pzz_zero=None)) is None
    assert m5.two_fill_pzz(_A(pzz_plus=0.6, pzz_zero=-1.2)) == (0.6, -1.2)


def test_half_a_two_fill_run_is_refused():
    """One flag alone would silently take the single-fill path under a
    stem that does not say so."""
    import argparse

    ap = argparse.ArgumentParser()
    m5.add_two_fill_args(ap)
    with pytest.raises(SystemExit):
        m5.check_two_fill_args(ap, ap.parse_args(["--pzz-plus", "0.6"]))
    with pytest.raises(SystemExit):
        m5.check_two_fill_args(
            ap, ap.parse_args(["--pzz-plus", "0.6", "--pzz-zero", "0.6"]))
    # both given and distinct: accepted
    m5.check_two_fill_args(
        ap, ap.parse_args(["--pzz-plus", "0.6", "--pzz-zero", "-1.2"]))


def test_the_published_pair_buys_the_factor_0p67():
    """dA(two-fill)/dA(single-fill) = P_single/sigma_P = 0.6/0.9."""
    assert m5.two_fill_err_ratio(PUBLISHED_PZZ, PUBLISHED_PAIR) == \
        pytest.approx(2.0 / 3.0, abs=1e-12)


def test_the_closed_form_agrees_with_the_analytic_error():
    """The same 0.67, this time as the ratio of the two error formulas
    the scripts actually call, at one fixed N."""
    n = 1.0e8
    single = float(np.sqrt(2.0 / n) / PUBLISHED_PZZ
                   / (np.sin(2 * np.pi / 24) / (2 * np.pi / 24)))
    two = float(reco.err_harmonic_ratio(n, list(PUBLISHED_PAIR), nbins=24))
    assert two / single == pytest.approx(2.0 / 3.0, abs=1e-10)
    assert two / single == pytest.approx(
        m5.two_fill_err_ratio(PUBLISHED_PZZ, PUBLISHED_PAIR), abs=1e-10)


def test_the_m_zero_fill_must_carry_twice_the_lever_arm_for_the_gain():
    """The 1.5x is real only if the source delivers |P_zz| = 1.2 at the
    same purity (plans/07 WP3).  At (0.6, -0.6) -- equal purity, opposite
    sign, i.e. no m = 0 enrichment -- the two-fill error is WORSE than the
    single-fill one, which is the statement that costs the gain."""
    assert m5.two_fill_err_ratio(0.6, (0.6, -0.6)) == pytest.approx(1.0)
    assert m5.two_fill_err_ratio(0.6, (0.6, 0.0)) == pytest.approx(2.0)


def test_the_fitted_error_of_a_pseudo_experiment_lands_on_0p67():
    """End to end through `reco.harmonic_ratio_fit`: two equal-luminosity
    fills drawn at the published amplitude, against the single-fill binned
    fit on the same total N."""
    from polligen.estimators import cos2phi_fit_err
    from polligen.sample import phi_histogram_pseudo

    rng = np.random.default_rng(20260915)
    n_total, amp = 2.0e8, 7.4e-3
    rows, edges = [], None
    for p in PUBLISHED_PAIR:
        c, edges = phi_histogram_pseudo(0.5 * n_total, p * amp,
                                        nbins=24, rng=rng)
        rows.append(c)
    fit = reco.harmonic_ratio_fit(np.asarray(rows, dtype=float),
                                  [0.5, 0.5], list(PUBLISHED_PAIR), edges)
    single = float(cos2phi_fit_err(n_total, PUBLISHED_PZZ, 24))
    assert fit["err"] / single == pytest.approx(2.0 / 3.0, abs=0.01)
    assert fit["amp"] == pytest.approx(amp, abs=5.0 * fit["err"])


def test_two_fill_plan_is_the_purity_matched_flip_plan_at_the_default():
    """`two_fill_plan((P, -2P))` and `bookkeeping.tensor_flip_plan(P)`
    must be the same two fills, so that the truth-level path and the
    reconstructed-level one are not two different run plans."""
    from polligen import bookkeeping as bk

    mine = m5.two_fill_plan(PUBLISHED_PAIR)
    theirs = bk.tensor_flip_plan(PUBLISHED_PZZ)
    assert len(mine.categories) == len(theirs.categories) == 2
    for a, b in zip(mine.categories, theirs.categories):
        assert np.allclose(a.populations, b.populations)
        assert a.lumi_fraction == pytest.approx(b.lumi_fraction)
        assert a.theta_s == pytest.approx(b.theta_s)


# --- the Yellow-Report binning -------------------------------------------

def test_the_default_binning_is_the_published_grid():
    """`analysis_grid()` threaded into `fom.project_rates` must be the
    call that function already makes by itself."""
    assert beams.analysis_grid() == {
        "nx": 40, "nq2": 30, "x_range": (1e-4, 1.0), "q2_range": (1.0, 2e3)}
    assert beams.analysis_grid("generic") == beams.analysis_grid("log40x30")
    assert beams.binning_tag() == ""
    assert beams.binning_tag("generic") == ""
    assert beams.binning_tag("yr") == "yr"


def test_an_unknown_binning_is_refused_rather_than_silently_generic():
    with pytest.raises(ValueError):
        beams.analysis_grid("40x30")
    with pytest.raises(ValueError):
        beams.binning_tag("YR")


def test_the_yr_lattice_is_five_bins_per_decade_in_both_variables():
    g = beams.analysis_grid("yr")
    xe, qe, _xc, _qc = kin.log_grid(g["x_range"], g["q2_range"],
                                    g["nx"], g["nq2"])
    assert beams.bins_per_decade(xe) == pytest.approx(5.0, abs=1e-12)
    assert beams.bins_per_decade(qe) == pytest.approx(5.0, abs=1e-12)
    # and the published grid is not: 10.0 and 9.09
    gg = beams.analysis_grid()
    xe0, qe0, _, _ = kin.log_grid(gg["x_range"], gg["q2_range"],
                                  gg["nx"], gg["nq2"])
    assert beams.bins_per_decade(xe0) == pytest.approx(10.0, abs=1e-12)
    assert beams.bins_per_decade(qe0) == pytest.approx(9.088, abs=1e-3)


def test_the_written_out_yr_edges_are_the_grid_that_is_actually_built():
    """`YR_X_EDGES`/`YR_Q2_EDGES` are the human-readable copy; production
    code goes through `log_grid`.  If the two ever part, the copy is the
    one that lies."""
    g = beams.analysis_grid("yr")
    xe, qe, _, _ = kin.log_grid(g["x_range"], g["q2_range"],
                                g["nx"], g["nq2"])
    assert np.allclose(beams.YR_X_EDGES, xe, rtol=1e-12, atol=0.0)
    assert np.allclose(beams.YR_Q2_EDGES, qe, rtol=1e-12, atol=0.0)


def test_the_yr_lattice_sits_on_the_decade_boundaries():
    """Anchored, not an arbitrary five-way split of this project's span:
    every decade boundary inside the range is an edge."""
    for decade in (1e-4, 1e-3, 1e-2, 1e-1, 1.0):
        assert np.min(np.abs(np.asarray(beams.YR_X_EDGES) / decade - 1.0)) \
            < 1e-12
    for decade in (1.0, 10.0, 100.0, 1000.0):
        assert np.min(np.abs(np.asarray(beams.YR_Q2_EDGES) / decade - 1.0)) \
            < 1e-12


def test_the_yr_q2_range_is_rounded_out_and_costs_nothing():
    """The Q2 span 1-2e3 is 3.3 decades, so the lattice edge above it is
    10^3.4 = 2512 and the range is rounded OUT to it rather than in --
    nothing inside the published span is dropped.  The extra row is empty
    at every configuration of both isotopes: the electron acceptance cuts
    off well below Q2 = 2e3, so rounding out buys coverage that costs no
    events and loses none."""
    from polli_fastsim import fom
    from polli_fastsim.structure import NuclearF2

    g = beams.analysis_grid("yr")
    assert g["q2_range"][1] == pytest.approx(10.0 ** 3.4)
    assert g["q2_range"][1] > 2e3
    sc = fom.Scenario(lumi_fb_per_nucleon=10.0, pol_ion_tensor=0.6)
    for ion in ("6Li", "7Li"):
        for cfg in beams.default_configs(ion):
            proj = fom.project_rates(cfg, sc, nuclear_f2=NuclearF2(cfg.ion),
                                     **g)
            above = proj.accepted & (proj.q2 > 2e3)
            assert proj.n_events[above].sum() == 0.0
