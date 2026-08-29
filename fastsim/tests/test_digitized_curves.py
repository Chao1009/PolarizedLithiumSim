"""The digitized theory tables load, and read back what the papers print.

Pins the plans/02 step 1.2.2 / 1.2.3 curves against values read off the
published figures, so a re-extraction that silently changes a calibration
or a curve -> label assignment cannot pass.  Provenance and the extraction
commands are in `polli_fastsim/data/SOURCES.md`.
"""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import polarized as P

TABLES = (P.CBT_TABLE, P.TMT_TABLE, P.MILLER_TABLE, P.CDKS_TABLE,
          "b1_cdks_q2set", "b1_miller_q2set")


def test_every_table_loads_with_a_strictly_increasing_x():
    for name in TABLES:
        t = P._load_curve(name)
        assert "x" in t and len(t) >= 2, name
        assert np.all(np.diff(t["x"]) > 0), name
        assert all(np.all(np.isfinite(v)) for v in t.values()), name
        lo, hi = P.curve_x_range(name)
        assert 0.0 < lo < hi, name


def test_cbt_figure_6_values():
    """CBT PLB 642:210 Fig. 6, 7Li at Q2 = 5, at x = 0.10/0.30/0.50/0.70."""
    x = np.array([0.10, 0.30, 0.50, 0.70])
    assert np.allclose(P.cbt_unpolarized_emc_ratio(x),
                       [1.024, 0.994, 0.939, 0.910], atol=0.005)
    assert np.allclose(P.cbt_polarized_emc_ratio(x),
                       [0.927, 0.930, 0.914, 0.905], atol=0.005)
    assert np.allclose(P.cbt_polarized_emc_ratio(x, eq=26),
                       [0.928, 0.927, 0.909, 0.921], atol=0.005)


def test_tmt_figure_4_values():
    """TMT PLB 783:247 Fig. 4, nuclear matter at Q2 = 10, as published."""
    x = np.array([0.10, 0.30, 0.50, 0.70])
    t = P._load_curve(P.TMT_TABLE)
    assert np.allclose(np.interp(x, t["x"], t["R_unpol"]),
                       [0.985, 0.942, 0.850, 0.906], atol=0.006)
    assert np.allclose(np.interp(x, t["x"], t["R_pol"]),
                       [0.926, 0.925, 0.849, 0.881], atol=0.005)


def test_ratio_of_effects_is_the_two_camps_statement():
    """'about twice' vs 'about equal' -- 2.25/1.69/1.41/1.14 against
    1.01/0.98/1.00/1.08 over the valence region, the numbers Report 0 and
    fastsim/README.md quote."""
    x = np.array([0.40, 0.45, 0.50, 0.60])
    assert np.allclose(P.cbt_ratio_of_effects(x),
                       [2.25, 1.69, 1.41, 1.14], atol=0.05)
    assert np.allclose(P.tmt_ratio_of_effects(x),
                       [1.01, 0.98, 1.00, 1.08], atol=0.05)
    assert P.cbt_ratio_of_effects(0.5) == pytest.approx(1.41, abs=0.05)
    assert P.tmt_ratio_of_effects(0.5) == pytest.approx(1.00, abs=0.05)


def test_the_common_baseline_leaves_cbt_untouched():
    """CBT computed 7Li itself, so its valence scale is exactly 1 and the
    money plot draws its published curve; TMT's nuclear matter is scaled
    down to 7Li strength by 0.397."""
    assert P.CBT_VALENCE_SCALE == 1.0
    assert P.TMT_VALENCE_SCALE == pytest.approx(0.397, abs=0.005)
    x = np.linspace(0.05, 0.85, 41)
    assert np.allclose(P.cbt_polarized_emc_ratio(x),
                       P._interp(P.CBT_TABLE, "R_pol_eq23", x))
    # and the transferred TMT curve stays a physical ratio everywhere the
    # money plot draws it -- the pointwise ratio form does not (it diverges
    # where either model's unpolarized curve crosses 1)
    xs = np.logspace(np.log10(0.005), np.log10(0.9), 300)
    tmt = P.tmt_polarized_emc_ratio(xs)
    assert np.all(np.isfinite(tmt)) and np.all(tmt > 0.7) and np.all(tmt < 1.2)


def test_the_two_published_curves_agree_at_low_x_and_part_in_the_valence_window():
    """Where the separation the money plot draws is real and where it is not.

    The transfer factor is a single constant defined over
    POLEMC_VALENCE_WINDOW and applied everywhere, so below that window it
    manufactures a separation the two papers do not have.  This pins both
    ends of that statement, because the report reads its headline off it
    (Report 0 section 5.3; `money_polemc.py` prints both separations).
    """
    lo = np.linspace(P.curve_x_range(P.CBT_TABLE)[0], 0.30, 200)
    raw_lo = np.abs(P.cbt_polarized_emc_ratio(lo)
                    - P.tmt_published_emc_ratio(lo))
    assert raw_lo.max() < 0.009          # the two CALCULATIONS agree here
    assert abs(P.cbt_polarized_emc_ratio(0.09)
               - P.tmt_published_emc_ratio(0.09)) < 0.004
    assert abs(P.cbt_polarized_emc_ratio(0.14)
               - P.tmt_published_emc_ratio(0.14)) < 0.002
    # ... while the transferred pair is drawn ~0.04 apart there
    trans_lo = np.abs(P.cbt_polarized_emc_ratio(lo)
                      - P.tmt_polarized_emc_ratio(lo))
    assert trans_lo.min() > 0.035
    assert trans_lo.min() > 4.0 * raw_lo.max()
    # in the valence window the two published curves genuinely part
    for x, floor in ((0.45, 0.03), (0.55, 0.08), (0.65, 0.10)):
        assert abs(P.cbt_polarized_emc_ratio(x)
                   - P.tmt_published_emc_ratio(x)) > floor
    # and the accessor is the untransferred table, not the scaled one
    g = np.linspace(0.05, 0.7, 33)
    assert np.allclose(P.tmt_published_emc_ratio(g),
                       P._interp(P.TMT_TABLE, "R_pol", g))
    assert np.allclose(P.tmt_polarized_emc_ratio(g),
                       1.0 - P.TMT_VALENCE_SCALE
                       * (1.0 - P.tmt_published_emc_ratio(g)))


def test_legacy_constant_mode_is_the_old_curve_exactly():
    x = np.array([0.01, 0.1, 0.3, 0.45, 0.7])
    assert np.allclose(P.cbt_polarized_emc_ratio(x, mode="constant"),
                       1.0 - 2.0 * (1.0 - P.unpolarized_emc_ratio(x)))
    assert np.allclose(P.tmt_polarized_emc_ratio(x, mode="constant"),
                       P.unpolarized_emc_ratio(x))
    old_sep = np.abs(P.cbt_polarized_emc_ratio(x, mode="constant")
                     - P.tmt_polarized_emc_ratio(x, mode="constant"))
    assert np.allclose(old_sep, np.abs(1.0 - P.unpolarized_emc_ratio(x)))


def test_b1_tables_are_absolute_and_carry_the_right_shape():
    """Miller reproduces HERMES at low x; CDKS is an order of magnitude
    below it at x >~ 0.2 and changes sign twice, which `0.1 * toy_b1`
    never did."""
    # The tables hold the published PER-DEUTERON b1 and the accessors halve
    # it, because every consumer pairs b1 with a per-nucleon F1.  At
    # x = 0.012 the digitized Miller total is 0.114 per deuteron, against
    # HERMES's measured 0.112 +- 0.055 +- 0.028 and the 0.105 of his own
    # TABLE I -- i.e. the curve does reproduce the datum it was built for.
    assert 2.0 * P.toy_b1(0.012, 2.5, None) == pytest.approx(0.114, abs=0.01)
    assert abs(P.b1_convolution(0.3, 2.5, None)) < 1e-3
    assert abs(P.b1_convolution(0.25, 2.5, None)) < 1e-3
    # the CDKS sum is its own SD + DD
    t = P._load_curve(P.CDKS_TABLE)
    assert np.allclose(t["xb1_theory1_sum"],
                       t["xb1_theory1_SD"] + t["xb1_theory1_DD"], atol=2e-7)
    assert np.allclose(t["xb1_theory2_sum"],
                       t["xb1_theory2_SD"] + t["xb1_theory2_DD"], atol=2e-7)
    # two sign changes in 0.01 < x < 0.5 (0.06 and 0.42)
    m = (t["x"] > 0.02) & (t["x"] < 0.5)
    s = np.sign(t["xb1_theory1_sum"][m])
    assert int(np.sum(np.diff(s) != 0)) == 2
    # the constant extrapolation below Miller's x = 0.01, which the
    # generator does sample and whose value the docstring quotes
    assert P.toy_b1(0.001, 4.0, 1.0) == pytest.approx(0.0647, abs=5e-4)
    assert P.toy_b1(1e-5, 4.0, 1.0) == P.toy_b1(0.001, 4.0, 1.0)


def test_cdks_figure_5_reproduces_figure_4_independently():
    """Two figures, two frames, two calibrations, same curve."""
    a = P._load_curve("b1_cdks_q2set")
    b = P._load_curve(P.CDKS_TABLE)
    x = np.array([0.05, 0.1, 0.3, 0.5, 0.8])
    assert np.allclose(np.interp(x, a["x"], a["xb1_theory1_q2_2p5"]),
                       np.interp(x, b["x"], b["xb1_theory1_sum"]), atol=2e-6)


def test_close_kumano_integral_is_reported():
    """int b1 dx = 0 when the sea is not tensor polarized (Close-Kumano).
    Reported, not enforced: the tables cover only what the papers plot, so
    a nonzero value is partly a statement about the missing ends.  The
    convolution respects it to 5e-4 over x = 0.01-1.59; Miller's total,
    which is the model that violates it (his Sec. V), does not."""
    cdks = P.close_kumano_integral("cdks")
    miller = P.close_kumano_integral("miller")
    assert abs(cdks) < 1e-3, cdks
    assert abs(miller) > 5 * abs(cdks), (miller, cdks)
    print("Close-Kumano integral: CDKS %.3e, Miller %.3e" % (cdks, miller))


def test_b1_li6_transfer_constants():
    """rank-2 x 2/6 by default; the legacy pair reproduces the old curve."""
    assert P.LI6_B1_RANK2_TRANSFER == pytest.approx(0.921947, abs=1e-6)
    assert P.LI6_B1_PER_NUCLEON == pytest.approx(1.0 / 3.0)
    assert P.b1_li6_from_deuteron(1.0) == pytest.approx(0.307316, abs=1e-6)
    assert P.b1_li6_from_deuteron(3.0, P.LI6_B1_LEGACY_TRANSFER,
                                  1.0) == pytest.approx(2.61)


def test_toy_modes_reproduce_the_pre_digitization_shapes():
    x = np.array([0.01, 0.05, 0.2, 0.5])
    f1 = np.array([12.0, 6.0, 2.0, 0.4])
    shape = (0.01 * x ** -0.2 * (1.0 - x / 0.20)) * np.exp(-3.0 * x) * f1
    assert np.allclose(P.toy_b1(x, 4.0, f1, mode="toy"), shape)
    assert np.allclose(P.b1_convolution(x, 4.0, f1, mode="toy"), 0.1 * shape)


def test_bad_modes_are_rejected():
    for fn in (P.cbt_polarized_emc_ratio, P.tmt_polarized_emc_ratio):
        with pytest.raises(ValueError):
            fn(0.3, mode="nonsense")
    with pytest.raises(ValueError):
        P.toy_b1(0.3, 4.0, 1.0, mode="nonsense")
    with pytest.raises(ValueError):
        P.cbt_polarized_emc_ratio(0.3, eq=99)


def test_digitized_b1_is_tapered_to_zero_at_the_elastic_edge():
    """b1/F1 is what the kernel forms.  Frozen at its x = 0.9 value while
    F1 falls as (1 - x)^3, the ratio diverged at the generator grid's
    x = 0.955 cell and drove the phi-averaged density negative
    (money_tagged_azz.py, 2026-08-28).  Above the table's last point the
    digitized b1 now falls with the same counting-rule power, so the ratio
    against a (1 - x)^3 F1 is frozen, not growing, and b1 vanishes at
    x = 1."""
    x_lo, x_hi = P.curve_x_range(P.MILLER_TABLE)
    assert x_hi == pytest.approx(0.9, abs=1e-3)
    f1 = lambda x: (1.0 - np.asarray(x)) ** 3  # noqa: E731
    xs = np.array([0.85, x_hi, 0.93, 0.955, 0.99, 1.0])
    b1 = P.toy_b1(xs, 4.0, f1(xs))
    assert np.all(np.isfinite(b1))
    assert b1[-1] == pytest.approx(0.0, abs=1e-12)
    ratio = b1[1:-1] / f1(xs[1:-1])
    assert np.allclose(ratio, ratio[0], rtol=1e-9)   # frozen at the table end
    assert np.all(np.abs(b1[2:]) < abs(b1[1]))       # and falling beyond it
    # inside the table nothing changed
    assert np.allclose(P.toy_b1(np.array([0.1, 0.5, 0.85]), 4.0, 1.0),
                       P.B1_PER_DEUTERON_TO_PER_NUCLEON
                       * P._interp(P.MILLER_TABLE, "b1",
                                   np.array([0.1, 0.5, 0.85])))
