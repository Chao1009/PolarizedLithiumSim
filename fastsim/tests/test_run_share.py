"""The run-plan share: what it scales, and what it must leave alone.

Today every published projection gives its observable the whole of the
10 fb^-1/u year in its own spin configuration, optics and isotope, so the
reaches of one report cannot be added up into a single year.  `run_share`
(plans/07 WP2) is the knob that prices any other division of the year, and
the arithmetic it obeys is the only thing a run-plan table needs:

  statistical errors  x 1/sqrt(share)
  event counts        x share
  L_5sigma, quoted as the luminosity the measurement must ACCUMULATE,
                      invariant -- the share buys wall-clock time, not
                      physics, so the YEARS to that reach go as 1/share.

The last of the three is the one a rescaling can silently get wrong, so it
is pinned exactly rather than approximately.
"""

import importlib.util
import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.polarized import ToyG1, toy_b1, toy_delta_gluon

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"


def _load(name):
    spec = importlib.util.spec_from_file_location("_" + name,
                                                  SCRIPTS / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CFG = beams.BeamConfig(10.0, beams.LI6, 99.5)


def test_a_quarter_share_quarters_the_events_and_doubles_the_errors():
    full = fom.Scenario(lumi_fb_per_nucleon=10.0)
    quarter = fom.Scenario(lumi_fb_per_nucleon=10.0, run_share=0.25)
    p_full = fom.project_rates(CFG, full)
    p_quarter = fom.project_rates(CFG, quarter)
    np.testing.assert_allclose(p_quarter.n_events, 0.25 * p_full.n_events,
                               rtol=1e-12)

    g1 = ToyG1()
    o_full = fom.project_observables(CFG, full, p_full, g1, toy_b1,
                                     toy_delta_gluon)
    o_quarter = fom.project_observables(CFG, quarter, p_quarter, g1, toy_b1,
                                        toy_delta_gluon)
    m = p_full.accepted & (p_full.n_events > 10)
    for key in ("err_azz", "err_a_cos2phi", "err_a_par", "err_g1_over_f1"):
        np.testing.assert_allclose(o_quarter[key][m], 2.0 * o_full[key][m],
                                   rtol=1e-9)
    # the asymmetries themselves are physics and do not move
    for key in ("a_par", "azz", "a_cos2phi"):
        np.testing.assert_allclose(o_quarter[key][m], o_full[key][m],
                                   rtol=1e-12)


def test_the_share_and_the_programme_luminosity_stay_separate_knobs():
    """A report must be able to say '10 fb^-1/u of programme, 1/3 of it
    here': the two multiply, and only their product is observable."""
    a = fom.Scenario(lumi_fb_per_nucleon=10.0, run_share=0.25)
    b = fom.Scenario(lumi_fb_per_nucleon=2.5)
    assert a.lumi_effective_fb_per_nucleon == pytest.approx(2.5)
    np.testing.assert_allclose(fom.project_rates(CFG, a).n_events,
                               fom.project_rates(CFG, b).n_events, rtol=1e-12)
    # ... but they are not the same field: the programme figure is the
    # machine's and the share is the run plan's
    assert a.lumi_fb_per_nucleon == 10.0 and b.run_share == 1.0


def test_a_non_positive_share_is_refused():
    for bad in (0.0, -0.5):
        with pytest.raises(ValueError):
            fom.Scenario(run_share=bad)


def test_l5sigma_is_exactly_invariant_under_the_share():
    """money_delta.reach_fb quotes the luminosity the measurement must
    accumulate.  The share scales the per-bin sig^2 terms and the event
    counts together, so the solved reach scales as 1/share and multiplying
    it back leaves the answer bit-for-bit where it was -- including the
    min-events floor, which enters at the same luminosity in both."""
    md = _load("money_delta")
    for share in (1.0, 0.5, 0.25, 3.0):
        got = md.reach_fb(CFG, 1e-3, 0.8, run_share=share)
        ref = md.reach_fb(CFG, 1e-3, 0.8, run_share=1.0)
        assert got == pytest.approx(ref, rel=1e-12), share
    # the significance accumulated per fb^-1 of PROGRAMME luminosity, on
    # the other hand, is linear in the share
    s_full = md.sig2_per_fb_at(CFG, 1e-3, 0.8, lumi_fb=20.0)
    s_quarter = md.sig2_per_fb_at(CFG, 1e-3, 0.8, lumi_fb=20.0, run_share=0.25)
    assert s_quarter == pytest.approx(0.25 * s_full, rel=1e-9)


def test_the_share_key_guards_the_published_stems():
    """Published figures are the share-1 ones; anything else must write to
    its own file (the guard of money_tagged_azz.output_stem)."""
    assert fom.run_share_tag(1.0) == ""
    assert fom.run_share_tag(0.25) == "share0p25"
    assert fom.run_share_tag(1.0 / 3.0).startswith("share0p333")
    keys = {fom.run_share_tag(s) for s in (1.0, 0.5, 0.25, 1.0 / 3.0, 2.0)}
    assert len(keys) == 5


def test_the_header_states_both_numbers():
    """'10 fb^-1/u of programme, 1/3 of it here' has to be readable off the
    header a script prints, or a table built from a log cannot say which of
    the two a number assumed."""
    line = fom.run_share_header(10.0, 0.25)
    assert "0.25" in line and "10 fb^-1/u" in line and "2.5 fb^-1/u" in line
