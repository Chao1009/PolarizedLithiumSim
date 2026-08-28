"""The min-events floor of the gluonometry reach, and where it belongs.

`money_delta.py` normalises sigma^2 to 1 fb^-1/nucleon and scales the reach
analytically.  The min-events floor -- the cut that keeps the Gaussian
counting significance honest -- does NOT scale with it, and until
2026-08-28 it was applied at that 1 fb^-1/u: a bin holding two events per
fb^-1, and forty at a 20 fb^-1 reach, was discarded and the reach was then
scaled up from the truncated sum.  `reach_from_terms` applies it at the
luminosity the reach is quoted at.

Measured when the change was made: the mask changes by -15 to +3 bins of
the 311-509 accepted, and L_5sigma moves by at most 3e-9 relative anywhere
on the plotted Delta/F1 range, while money_delta_realistic.py reproduces
its frozen 131.26 / 274.64 fb^-1/u exactly.  The defect was real and its
numerical effect is nil; these tests pin both halves of that statement.
"""

import importlib.util
import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams
from polli_fastsim.inputs import get_backends

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"


def _load(name):
    spec = importlib.util.spec_from_file_location("_" + name,
                                                  SCRIPTS / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


md = _load("money_delta")


def test_reach_solves_l_times_sig2_equals_the_target():
    """Two bins, one rich and one poor, with the poor one below the floor
    at 1 fb^-1/u and above it at the reach."""
    n_events = np.array([1000.0, 3.0])
    terms = np.array([1.0, 0.5])
    lumi = md.reach_from_terms(terms, n_events, min_events=10, target=25.0)
    # both bins are in at the solution: 10/3 = 3.3 fb^-1/u admits the poor
    # one, and 25/1.5 = 16.7 is past that
    assert lumi == pytest.approx(25.0 / 1.5)
    assert np.all(n_events * lumi >= 10.0)
    # the floor at 1 fb^-1/u would have kept only the rich bin, and the
    # reach would have been scaled up from that truncated sum
    old = 25.0 / terms[n_events >= 10.0].sum()
    assert old == pytest.approx(25.0)          # the rich bin alone
    assert lumi < old


def test_reach_stops_at_the_floor_when_the_target_is_met_below_it():
    """One bin that is worth 5 sigma the moment it has ten events: the
    answer is the luminosity at which it does, not a smaller one."""
    lumi = md.reach_from_terms(np.array([100.0]), np.array([2.0]),
                               min_events=10, target=25.0)
    assert lumi == pytest.approx(5.0)          # 10 events / 2 per fb^-1


def test_reach_is_monotone_in_the_target():
    n_events = np.array([500.0, 40.0, 6.0, 0.4])
    terms = np.array([0.9, 0.4, 0.2, 0.1])
    reaches = [md.reach_from_terms(terms, n_events, target=t)
               for t in (1.0, 4.0, 25.0, 100.0)]
    assert reaches == sorted(reaches)


def test_empty_and_dead_inputs_give_an_infinite_reach():
    assert md.reach_from_terms(np.array([]), np.array([])) == np.inf
    assert md.reach_from_terms(np.zeros(3), np.ones(3)) == np.inf


def test_bin_terms_drops_non_finite_bins():
    """The floor used to remove NaN bins by accident (NaN >= 10 is False).
    Now that the threshold can fall below one event, the drop has to be
    explicit or a grid outside its fit range poisons the sum."""
    for cfg in beams.default_configs("6Li"):
        terms, n_events = md.bin_terms(cfg, 1e-3, 0.8,
                                       base=get_backends("toy")["base"])
        assert terms.size and terms.size == n_events.size
        assert np.all(np.isfinite(terms)) and np.all(np.isfinite(n_events))


def test_the_published_toy_reach_did_not_move():
    """What the correction cost at the published point: nothing.  The bins
    it recovers there hold 1.7 to 9 events per fb^-1/u out of 311-509
    accepted, and they sit where the cos 2phi amplitude is smallest."""
    base = get_backends("toy")["base"]
    published = {"e(5) x 6Li(40.8/u)": 16.7,
                 "e(10) x 6Li(99.5/u)": 16.3,
                 "e(18) x 6Li(137.5/u)": 21.8}
    for cfg in beams.default_configs("6Li"):
        terms, n_events = md.bin_terms(cfg, 1e-3, 0.8, base=base)
        new = md.reach_from_terms(terms, n_events, min_events=10)
        old = 25.0 / terms[n_events >= 10].sum()      # the floor at 1 fb^-1/u
        assert new == pytest.approx(old, rel=1e-8)
        assert new == pytest.approx(published[cfg.label()], abs=0.05)


def test_the_diagnostic_and_the_realistic_script_use_the_same_solver():
    """`diag_sig2_grid.py` used to carry a verbatim copy of the
    significance under a 'do not modify' comment; both it and
    `money_delta_realistic.py` now import this one, so the three cannot
    disagree about where the floor goes."""
    for name in ("diag_sig2_grid", "money_delta_realistic"):
        assert _load(name).reach_from_terms.__doc__ == \
            md.reach_from_terms.__doc__
