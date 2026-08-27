"""What `bookkeeping.helicity_flip_plan(...).pzz_true` actually is.

plans/08 B4: `pzz_true` is the number the money scripts divide by when
they turn a fitted harmonic into an asymmetry, and until now nothing
constrained it.  This file pins it for every branch of
`helicity_flip_plan`:

  * j = 1/2 -- exactly 0.0 (spin 1/2 has no rank-2 moment at all), and
    never a divisor anywhere: the only two callers of the flip plan
    (`scripts/closure_fom.py` "apar", `scripts/tagged_polarimetry_7li.py`)
    feed `estimators.apar_flip`, which divides by pe * pz alone, while
    every `plan.pzz_true` division in the tree (`money_cos2phi.py`,
    `money_cos2phi_reco.py`, `money_delta_extraction.py`,
    `phase_space_bins.py`) is fed by a spin-1 `transverse_tensor_plan`
    or `tensor_flip_plan` with an explicit pzz;
  * pzz = None -- the spin-temperature (max-entropy) fill, whose vector
    polarization DRAGS a rank-2 moment along with it.  Checked against
    the closed forms of the geometric population ladder p_m ~ t^m
    derived below, against the density matrix built from `spin.py`'s
    angular-momentum operators, and against the small-pz limit;
  * an explicit pzz -- closure, which is what the money scripts assume;
  * pz -> -pz -- the rank-2 moment is even in pz, the property that
    makes the helicity-flip pattern tensor-blind.

Closed forms used as the independent construction (spin temperature,
p_m proportional to t^m with t = exp(beta), populations ordered
m = +J ... -J).  Writing the ladder as t^k, k = 0 ... 2J:

  J = 1:    S = t^2 + t + 1
            P_z  = (t^2 - 1) / S,      P_zz = <3 J_z^2 - 2> = (t-1)^2 / S
            which satisfies the spin-temperature relation
            P_zz = 2 - sqrt(4 - 3 P_z^2) already used in test_spin.py.
  J = 3/2:  P_z  = (t-1)(3t^2 + 4t + 3) / (3 (1+t)(1+t^2))
            T    = <3 J_z^2 - J(J+1)>/3 = (t-1)^2 / (1 + t^2)

Both give rational anchors at t = 3: (P_z, P_zz) = (8/13, 4/13) for
J = 1 and (P_z, T) = (7/10, 2/5) for J = 3/2 -- exact numbers that owe
nothing to the module under test.  Small pz: expanding p_m ~ e^(beta m)
to second order gives beta = 3 P_z/(J+1) and a rank-2 moment
beta^2/3 (J = 1) and beta^2/2 (J = 3/2), i.e. (3/4) P_z^2 and
(18/25) P_z^2.
"""

import inspect
import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen import bookkeeping as bk  # noqa: E402
from polligen import estimators as est  # noqa: E402
from polligen import spin  # noqa: E402

PE = 0.7


def _spin_temperature_ladder(j, pz):
    """Populations p_m ~ t^m with <J_z>/J = pz, solved for t here rather
    than taken from `spin.populations_maxent` (bisection in t, not in
    beta, and no reuse of the module under test)."""
    ms = spin.m_values(j)

    def vector(t):
        w = t ** ms
        return float((ms * w).sum() / w.sum() / j)

    lo, hi = 1e-12, 1e12
    for _ in range(400):
        mid = np.sqrt(lo * hi)          # bisection in log t
        if vector(mid) < pz:
            lo = mid
        else:
            hi = mid
    t = np.sqrt(lo * hi)
    w = t ** ms
    return tuple(w / w.sum()), t


def test_spin_half_flip_plan_has_no_rank_two_moment():
    """A spin-1/2 fill has no rank-2 moment: pzz_true must be exactly
    zero, not merely small, and it must survive the polarimetry smear
    as zero (a smeared 1e-17 divisor would be worse than a zero one)."""
    plan = bk.helicity_flip_plan(0.5, 0.7, PE)
    assert plan.pzz_true == 0.0
    assert plan.measured["pzz"] == 0.0
    smeared = bk.RunPlan(plan.categories, pe_true=plan.pe_true,
                         pz_true=plan.pz_true, pzz_true=plan.pzz_true,
                         delta_p_over_p=0.03)
    assert smeared.measured["pzz"] == 0.0
    assert smeared.measured["pz"] != 0.7          # the smear does fire
    # populations are the plain (1+pz)/2, (1-pz)/2 doublet
    np.testing.assert_allclose(plan.categories[0].populations,
                               (0.85, 0.15), atol=1e-15)


def test_nothing_downstream_divides_a_flip_plan_by_pzz():
    """The estimator a helicity-flip plan feeds divides by pe * pz and
    takes no pzz at all; it recovers an injected A_par from counts built
    by hand, so the spin-1/2 zero can never reach a denominator."""
    assert "pzz" not in inspect.signature(est.apar_flip).parameters
    plan = bk.helicity_flip_plan(0.5, 0.6, PE)
    a_par = 0.037
    scale = plan.pe_true * plan.pz_true
    n_plus = 1.0e6 * (1.0 + scale * a_par)
    n_minus = 1.0e6 * (1.0 - scale * a_par)
    assert est.apar_flip(n_plus, n_minus, plan.pe_true,
                         plan.pz_true) == pytest.approx(a_par, rel=1e-12)


@pytest.mark.parametrize("j,pz,pzz", [(1.0, 8.0 / 13.0, 4.0 / 13.0),
                                      (1.5, 0.7, 0.4)])
def test_spin_temperature_pzz_at_the_rational_anchor(j, pz, pzz):
    """t = 3 makes both ladders rational: p_m ~ (9, 3, 1)/13 for J = 1
    and (27, 9, 3, 1)/40 for J = 3/2, giving (P_z, rank-2) = (8/13,
    4/13) and (7/10, 2/5).  The rank-2 value is re-derived here from
    the density matrix and `spin.py`'s angular-momentum operators, not
    from `moments_along_axis` (the route the branch under test uses)."""
    plan = bk.helicity_flip_plan(j, pz, PE)
    assert plan.pzz_true == pytest.approx(pzz, abs=1e-12)
    ladder = np.array([3.0 ** k for k in range(int(round(2 * j)), -1, -1)])
    pops = tuple(ladder / ladder.sum())
    np.testing.assert_allclose(plan.categories[0].populations, pops,
                               atol=1e-12)
    rho = spin.rho_from_populations(j, pops)
    assert spin.tensor_polarization(rho, j) == pytest.approx(pzz, abs=1e-12)


@pytest.mark.parametrize("j", [1.0, 1.5])
@pytest.mark.parametrize("pz", [0.05, 0.3, 0.7, 0.9])
def test_spin_temperature_pzz_across_pz(j, pz):
    """pzz = None is the physical vector fill, and its rank-2 moment is
    not a free parameter.  Independent construction: solve the geometric
    ladder here, then take the moment through the operator trace."""
    plan = bk.helicity_flip_plan(j, pz, PE)
    pops, _ = _spin_temperature_ladder(j, pz)
    rho = spin.rho_from_populations(j, pops)
    assert plan.pzz_true == pytest.approx(spin.tensor_polarization(rho, j),
                                          abs=1e-10)
    if abs(j - 1.0) < 1e-9:
        # the published spin-1 relation, used elsewhere in the suite
        assert plan.pzz_true == pytest.approx(
            2.0 - np.sqrt(4.0 - 3.0 * pz * pz), abs=1e-12)
    assert plan.pzz_true > 0.0      # a vector fill is prolate, never zero


@pytest.mark.parametrize("j,coeff", [(1.0, 0.75), (1.5, 0.72)])
def test_spin_temperature_pzz_small_pz_limit(j, coeff):
    """Second-order limit of p_m ~ exp(beta m): beta = 3 pz/(J+1) and the
    rank-2 moment is beta^2/3 (J = 1) / beta^2/2 (J = 3/2), i.e.
    (3/4) pz^2 and (18/25) pz^2.  This is where the tensor contamination
    of a vector run is smallest and the ratio is most diagnostic."""
    pz = 1e-4
    plan = bk.helicity_flip_plan(j, pz, PE)
    assert plan.pzz_true / pz ** 2 == pytest.approx(coeff, rel=1e-6)


@pytest.mark.parametrize("j,pz,pzz", [(1.0, 0.6, 0.35), (1.0, 0.5, -0.4),
                                      (1.0, 0.0, 1.0), (1.5, 0.5, 0.2),
                                      (1.5, 0.0, -1.0), (1.5, 0.3, -0.1)])
def test_explicit_pzz_closes(j, pz, pzz):
    """The money scripts divide by the pzz they asked for, so asking must
    be getting.  (The argument is P_zz in [-2, 1] for J = 1 and the
    normalized T in [-1, 1] for J = 3/2 -- different quantities behind
    one keyword, which is the reason to pin both.)"""
    plan = bk.helicity_flip_plan(j, pz, PE, pzz=pzz)
    assert plan.pzz_true == pytest.approx(pzz, abs=1e-12)
    assert plan.pz_true == pz
    rho = spin.rho_from_populations(j, plan.categories[0].populations)
    assert spin.tensor_polarization(rho, j) == pytest.approx(pzz, abs=1e-12)
    assert spin.vector_polarization(rho, j)[2] == pytest.approx(pz,
                                                                abs=1e-12)


@pytest.mark.parametrize("j,pz,kwargs", [(1.0, 0.7, {}), (1.5, 0.7, {}),
                                         (1.0, 0.6, {"pzz": 0.35}),
                                         (1.5, 0.5, {"pzz": 0.2})])
def test_pzz_true_is_even_in_pz(j, pz, kwargs):
    """The rank-2 moment cannot know the sign of the vector one, so the
    two halves of a helicity-flip pair -- and a plan built at -pz --
    carry the same pzz_true.  That is what makes the flip tensor-blind:
    the tensor term is common mode and cancels in (N+ - N-)/(N+ + N-)."""
    up = bk.helicity_flip_plan(j, pz, PE, **kwargs)
    down = bk.helicity_flip_plan(j, -pz, PE, **kwargs)
    assert up.pzz_true == pytest.approx(down.pzz_true, abs=1e-12)
    assert up.pz_true == -down.pz_true
    # both fills of one plan are the same ion state (only the electron
    # helicity flips), so each carries the full pzz_true
    for cat in up.categories:
        assert cat.moments()[1] == pytest.approx(up.pzz_true, abs=1e-12)
    assert [c.lam_e for c in up.categories] == [+1, -1]


def test_pzz_true_is_the_moment_along_the_fill_axis():
    """pzz_true is the moment along the quantization axis, not along the
    lab z: tilting theta_S leaves it alone while the lab moment picks up
    P2(cos theta_S).  The kernel works in the axis frame, so the axis
    moment is the right divisor -- and it means a theta_S sweep must not
    be read as a change of the beam's tensor polarization."""
    pz, theta_s = 0.7, 0.9
    flat = bk.helicity_flip_plan(1.0, pz, PE)
    tilted = bk.helicity_flip_plan(1.0, pz, PE, theta_s=theta_s)
    assert tilted.pzz_true == pytest.approx(flat.pzz_true, abs=1e-15)
    dens = spin.SpinDensity(1.0, tilted.categories[0].populations, theta_s,
                            0.0)
    p2 = 0.5 * (3.0 * np.cos(theta_s) ** 2 - 1.0)
    assert dens.lab_moments()["tensor_zz"] == pytest.approx(
        flat.pzz_true * p2, abs=1e-12)


def test_unsupported_spin_is_refused():
    with pytest.raises(ValueError):
        bk.helicity_flip_plan(2.0, 0.5, PE, pzz=0.3)
    with pytest.raises(ValueError):
        bk.helicity_flip_plan(1.0, 0.5, PE, pzz=-2.0)   # unphysical corner
