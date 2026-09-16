"""The 95% CL exclusion contour of the gluonometry reach figure.

plans/02 Step 1.3 item 3 asks for the exclusion contour next to the 5 sigma
discovery one and says what it must be: "the same curve shifted by
(1.645/5)^2".  `money_delta.py --cl-band` draws it, and re-solves it from
the per-bin sigma^2 at target = 1.645^2 rather than multiplying the 5 sigma
answer by 0.108241 -- so that if the min-events floor ever bound at the
lower luminosity, the contour would show it instead of inheriting a shift
that is no longer true.

Over the plotted Delta/F1 range and all six published curves the floor
barely grazes -- the ratio is the flat 0.108241 to five digits everywhere,
the worst departure being 1.1e-5 relative at the top of the scale range
(e(18), P_zz = 0.8) -- and these tests pin that, together with the filename
guard that keeps the published single-contour PNG out of reach of a
two-contour run.
"""

import importlib.util
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.inputs import get_backends

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"


def _load(name):
    spec = importlib.util.spec_from_file_location("_" + name,
                                                  SCRIPTS / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


MD = _load("money_delta")
EXPECTED = (1.645 / 5.0) ** 2           # 0.108241


def test_targets_are_the_gaussian_points():
    assert MD.TARGET_5SIG == 25.0
    assert MD.TARGET_95 == 1.645 ** 2
    assert round(MD.TARGET_95 / MD.TARGET_5SIG, 5) == round(EXPECTED, 5)
    assert round(EXPECTED, 5) == 0.10824


def test_exclusion_is_the_discovery_curve_shifted_by_the_ratio():
    """The number the script prints, for every configuration and P_zz."""
    base = get_backends("toy")["base"]
    for cfg in beams.default_configs("6Li"):
        for pzz in (0.60, 0.80):
            terms, n_events = MD.bin_terms(cfg, 1e-3, pzz, base=base)
            l5 = MD.reach_from_terms(terms, n_events,
                                     target=MD.TARGET_5SIG)
            l95 = MD.reach_from_terms(terms, n_events, target=MD.TARGET_95)
            assert np.isfinite(l5) and l5 > 0
            assert round(l95 / l5, 5) == round(EXPECTED, 5)


def test_the_shift_is_flat_across_the_plotted_scale_range():
    """Not just at the reference point, and not just at one curve: every
    one of the six curves the figure draws, at all fifteen plotted
    scales."""
    base = get_backends("toy")["base"]
    scales = np.logspace(-3.3, -1.7, 15)
    worst = 0.0
    for cfg in beams.default_configs("6Li"):
        for pzz in (0.60, 0.80):
            terms, n_events = MD.bin_terms(cfg, 1e-3, pzz, base=base)
            ratio = np.array([
                MD.reach_from_terms(terms * (s / 1e-3) ** 2, n_events,
                                    target=MD.TARGET_95)
                / MD.reach_from_terms(terms * (s / 1e-3) ** 2, n_events,
                                      target=MD.TARGET_5SIG)
                for s in scales])
            # the floor grazes the top of the scale range: the largest
            # scale of e(18) at P_zz = 0.8 departs by 1.1e-5 relative, the
            # worst of the six, and every point still rounds to 0.10824
            np.testing.assert_allclose(ratio, EXPECTED, rtol=5e-5)
            assert set(np.round(ratio, 5)) == {round(EXPECTED, 5)}
            worst = max(worst, float(np.abs(ratio / EXPECTED - 1).max()))
    assert 1e-6 < worst < 2e-5, worst      # pins the grazing itself


def test_the_two_contour_figure_cannot_overwrite_the_published_one():
    assert MD.cl_band_tag(False) == ""
    assert MD.cl_band_tag(True) == "cl95"
    # and it composes with the run-share guard rather than replacing it
    assert fom.run_share_tag(1.0) == ""
    assert fom.run_share_tag(0.5) == "share0p5"
