"""Unified Delta-model registry: moment constraints, shapes, interface."""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, delta_models as dm
from polli_fastsim.polarized import toy_delta_gluon
from polli_fastsim.structure import NuclearF2


@pytest.fixture(scope="module")
def f1_func():
    nf2 = NuclearF2(beams.LI6)
    return lambda x, q2: nf2.f1a(x, q2) / beams.LI6.A


def _moment(model, f1_func, q2):
    """Numerical int_0^1 x Delta(x, q2) dx."""
    x = np.unique(np.concatenate([np.logspace(-5, np.log10(0.5), 400),
                                  np.linspace(0.5, 0.9999, 200)]))
    q2a = np.full_like(x, q2)
    d = model(x, q2a, f1_func(x, q2a))
    return float(np.trapz(x * d, x))


def test_moment_a_satisfies_sum_rule(f1_func):
    q2_ref = 7.4
    model = dm.make("moment_A", f1_func=f1_func, q2_ref=q2_ref,
                    alphas=dm.alpha_s_lo)
    # at Q2 = q2_ref the constraint holds by construction:
    # int x Delta dx = c * alpha_s(q2_ref)
    expect = dm.C_BAG * dm.alpha_s_lo(q2_ref)
    assert _moment(model, f1_func, q2_ref) == pytest.approx(expect,
                                                            rel=1e-3)


def test_moment_b_analytic_sum_rule(f1_func):
    model = dm.make("moment_B", alphas=dm.alpha_s_lo)
    for q2 in (3.0, 10.0):
        expect = dm.C_BAG * dm.alpha_s_lo(q2)
        assert _moment(model, f1_func, q2) == pytest.approx(expect,
                                                            rel=1e-3)


def test_default_alphas_resolution(f1_func):
    # alphas=None resolves to the parton table when grids exist, else LO;
    # the sum rule must hold with whatever was resolved
    resolved = dm.default_alpha_s()
    model = dm.make("moment_B")
    q2 = 5.0
    assert _moment(model, f1_func, q2) == pytest.approx(
        dm.C_BAG * float(resolved(q2)), rel=1e-3)
    # table alpha_s, when present, sits below the LO analytic at Q2 ~ 5
    table = dm.alpha_s_table()
    if table is not None:
        assert float(table(5.0)) < float(dm.alpha_s_lo(5.0))


def test_dilution_scales_linearly(f1_func):
    full = dm.make("moment_B")
    third = dm.make("moment_B", dilution=1.0 / 3.0)
    x = np.array([0.05, 0.2])
    q2 = np.array([5.0, 5.0])
    f1 = f1_func(x, q2)
    assert np.allclose(third(x, q2, f1), full(x, q2, f1) / 3.0)


def test_variants_and_peak_normalization():
    for name, (a, b) in dm.VARIANTS.items():
        xp = a / (a + b)
        assert dm.shape_normalized(xp, name) == pytest.approx(1.0)
        assert dm.shape_normalized(0.999, name) < 0.05


def test_interp_a_magnitude_matches_money_delta(f1_func):
    # money_delta_20260729 quotes |A| ~ 0.30-0.32 (mid_x, EPPS21 grids);
    # with the toy F1 backend the same solver should land near that.
    a = dm.solve_A_interp_a(f1_func, q2_ref=7.4, variant="mid_x")
    assert -0.5 < a < -0.15


def test_toy_passthrough():
    model = dm.make("toy", scale=1e-2)
    x = np.array([0.05, 0.3])
    q2 = np.array([4.0, 4.0])
    f1 = np.array([2.0, 0.5])
    assert np.allclose(model(x, q2, f1),
                       toy_delta_gluon(x, q2, f1, scale=1e-2))


def test_registry_and_info(f1_func):
    assert dm.available() == ["moment_A", "moment_B", "toy"]
    with pytest.raises(KeyError):
        dm.make("nope")
    m = dm.make("moment_A", f1_func=f1_func, q2_ref=7.4,
                dilution=1.0 / 3.0)
    assert "moment_A" in m.info() and "mid_x" in m.info()
