"""The two measured attenuations of the amplitude, D and A.

`recopseudo` measures a phi' dilution bin by bin (`RecoResponse.bin_summary`
-> `dilution_phi`) and a beta dilution for the tagged channel
(`CoherentResponse.truth_reference` -> `dilution_beta`); until now neither
could reach the analytic layer, which had no field for them (plans/05:370).
`Scenario.dilution` and `Scenario.acceptance` are that field, and the
arithmetic they obey is one line: a fit whose amplitude is D*A times the
physics one has a statistical error 1/(D*A) times the ideal one.

Two things this must NOT do, both pinned below:

  * it must not touch the event count, and it must not enter as 1/sqrt --
    that is what a luminosity share does (`run_share`), and the two laws
    differ by a square root;
  * at the defaults (1.0, 1.0) it must be the identity BIT FOR BIT, since
    every published figure is made there.

The flow-back is by NUMBER, not by import: fastsim never imports evgen
(the discipline is test-guarded), so a script passes the measured D in.
"""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.asymmetries import (err_a_parallel, err_azz,
                                       err_cos2phi_amplitude)
from polli_fastsim.polarized import ToyG1, toy_b1, toy_delta_gluon

CFG = beams.BeamConfig(10.0, beams.LI6, 99.5)
ERR_KEYS = ("err_azz", "err_a_cos2phi", "err_a_par", "err_g1_over_f1")


def _observables(**kw):
    sc = fom.Scenario(lumi_fb_per_nucleon=10.0, **kw)
    proj = fom.project_rates(CFG, sc)
    obs = fom.project_observables(CFG, sc, proj, ToyG1(), toy_b1,
                                  toy_delta_gluon)
    return sc, proj, obs


def test_half_the_dilution_exactly_doubles_every_error():
    """The acceptance criterion of the flow-back: 1/(D*A), not 1/sqrt(D*A).

    Asserted with `==`, not a tolerance: 0.5 and 2 are exact in binary, so
    a law that is right is right in the last bit, and anything that has to
    be rounded into agreement is a different law."""
    _, proj, full = _observables()
    _, _, half = _observables(dilution=0.5)
    m = proj.accepted & (proj.n_events > 10)
    for key in ERR_KEYS:
        assert np.array_equal(half[key][m], 2.0 * full[key][m]), key


def test_the_acceptance_is_the_same_factor_and_the_two_multiply():
    _, proj, full = _observables()
    _, _, split = _observables(dilution=0.5, acceptance=0.25)
    m = proj.accepted & (proj.n_events > 10)
    for key in ERR_KEYS:
        assert np.array_equal(split[key][m], 8.0 * full[key][m]), key
    assert fom.Scenario(dilution=0.5, acceptance=0.25).analyzing_power == 0.125


def test_the_defaults_are_the_identity_bit_for_bit():
    """Every published figure is made at D = A = 1, so the new fields must
    leave the three error functions returning exactly what they returned
    without them -- compared against `asymmetries` directly, so the check
    does not merely compare the new code with itself."""
    sc, proj, obs = _observables()
    assert sc.analyzing_power == 1.0
    n = proj.n_events
    assert np.array_equal(obs["err_azz"], err_azz(n, sc.pol_ion_tensor))
    assert np.array_equal(obs["err_a_cos2phi"],
                          err_cos2phi_amplitude(n, sc.pol_ion_tensor))
    assert np.array_equal(obs["err_a_par"],
                          err_a_parallel(n, sc.pol_electron,
                                         sc.pol_ion_vector))
    _, _, explicit = _observables(dilution=1.0, acceptance=1.0)
    for key in ERR_KEYS:
        assert np.array_equal(explicit[key], obs[key]), key


def test_a_dilution_is_not_a_luminosity_share():
    """D attenuates the amplitude; f buys wall-clock time.  The event
    count sees only f, and the error law differs by the square root -- so
    D = 0.25 and f = 0.25 must NOT agree."""
    _, p_full, full = _observables()
    _, p_dil, dil = _observables(dilution=0.25)
    _, _, share = _observables(run_share=0.25)
    assert np.array_equal(p_dil.n_events, p_full.n_events)
    m = p_full.accepted & (p_full.n_events > 10)
    assert np.array_equal(dil["err_azz"][m], 4.0 * full["err_azz"][m])
    np.testing.assert_allclose(share["err_azz"][m], 2.0 * full["err_azz"][m],
                               rtol=1e-12)
    assert not np.allclose(dil["err_azz"][m], share["err_azz"][m])
    # and the asymmetries themselves are physics: neither knob moves them
    for key in ("a_par", "azz", "a_cos2phi"):
        assert np.array_equal(dil[key][m], full[key][m]), key


def test_the_significance_maps_carry_the_same_factor():
    """sig = |A|/delta(A) is what the money plots read off, so a diluted
    projection must lose significance in step with its error."""
    _, proj, full = _observables()
    _, _, half = _observables(dilution=0.5)
    m = proj.accepted & (proj.n_events > 10)
    for key in ("sig_azz", "sig_a_cos2phi", "sig_a_par"):
        assert np.array_equal(half[key][m], 0.5 * full[key][m]), key


def test_a_non_positive_dilution_or_acceptance_is_refused():
    for bad in (0.0, -0.5):
        with pytest.raises(ValueError, match="dilution"):
            fom.Scenario(dilution=bad)
        with pytest.raises(ValueError, match="acceptance"):
            fom.Scenario(acceptance=bad)


def test_the_measured_dilution_travels_as_a_number_not_an_import():
    """`recopseudo.RecoResponse.bin_summary`'s `dilution_phi` is a float in
    a dict; the fast-sim takes it as `Scenario(dilution=...)`.  fastsim
    importing evgen would invert the dependency, so it stays out."""
    import ast

    import polli_fastsim
    root = pathlib.Path(polli_fastsim.__file__).resolve().parent
    offenders = []
    for mod in sorted(root.glob("*.py")):
        for node in ast.walk(ast.parse(mod.read_text())):
            names = ([a.name for a in node.names]
                     if isinstance(node, ast.Import) else
                     [node.module or ""] if isinstance(node, ast.ImportFrom)
                     else [])
            if any(n.split(".")[0] == "polligen" for n in names):
                offenders.append(mod.name)
    assert offenders == []
    measured = {"dilution_phi": 0.8312}          # shape of a bin_summary row
    assert fom.Scenario(dilution=measured["dilution_phi"]).analyzing_power \
        == 0.8312
