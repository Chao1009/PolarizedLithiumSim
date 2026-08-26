"""Gate 2 (plans/05 §5.4): master formula == asymmetries.py bin-wise.

The doubly polarized kernel must reproduce the fastsim's analytic
asymmetries exactly, sector by sector: vector -> a_parallel, tensor ->
azz, transverse tensor -> a_cos2phi.  Toy backends here; the same
identities on PDF-grid backends run in the grid-marked test (auto-skip).
"""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen.xsec import EventSpinState, InclusiveKernel, g2_ww  # noqa: E402

from polli_fastsim import beams  # noqa: E402
from polli_fastsim.asymmetries import (  # noqa: E402
    a_cos2phi, a_parallel, azz)
from polli_fastsim.polarized import toy_b1, toy_delta_gluon  # noqa: E402


def _grid(s):
    x = np.logspace(-3, -0.15, 14)
    q2 = np.logspace(0.1, 2.3, 11)
    xx, qq = np.meshgrid(x, q2, indexing="ij")
    y = qq / (s * xx)
    ok = (y > 0.01) & (y < 0.95)
    return xx[ok], qq[ok], y[ok]


@pytest.fixture(scope="module")
def setup6():
    config = beams.default_configs("6Li")[1]
    s = config.sqrt_s_per_nucleon**2
    return config, s


def test_vector_sector_matches_a_parallel(setup6):
    _, s = setup6
    kern = InclusiveKernel(beams.LI6)  # no tensor SFs -> pure vector sector
    x, q2, y = _grid(s)
    t = kern.tables(x, q2)
    pe = 0.7
    wp, _, _ = kern.amplitudes(t, x, q2, s,
                               EventSpinState(+1, pe, 1.0, 1.0), )
    wm, _, _ = kern.amplitudes(t, x, q2, s,
                               EventSpinState(-1, pe, 1.0, 1.0))
    measured = (wp - wm) / (2.0 + wp + wm)
    expected = pe * a_parallel(t["g1"], t["f1"], y, x, q2)
    np.testing.assert_allclose(measured, expected, rtol=1e-12)


def test_vector_sector_spin32(setup6):
    _, s = setup6
    kern = InclusiveKernel(beams.LI7)
    x, q2, y = _grid(s)
    t = kern.tables(x, q2)
    # m = +3/2: full vector polarization m/J = 1
    wp, _, _ = kern.amplitudes(t, x, q2, s, EventSpinState(+1, 1.0, 1.5, 1.5))
    wm, _, _ = kern.amplitudes(t, x, q2, s, EventSpinState(-1, 1.0, 1.5, 1.5))
    measured = (wp - wm) / (2.0 + wp + wm)
    expected = a_parallel(t["g1"], t["f1"], y, x, q2)
    np.testing.assert_allclose(measured, expected, rtol=1e-12)
    # rank-2 slots default to zero: unpolarized electron -> no modulation
    w0, a1, a2 = kern.amplitudes(t, x, q2, s, EventSpinState(0, 0.0, 1.5, 1.5))
    assert np.all(w0 == 0.0) and np.all(a1 == 0.0) and np.all(a2 == 0.0)


def test_tensor_sector_matches_azz(setup6):
    _, s = setup6
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    x, q2, y = _grid(s)
    t = kern.tables(x, q2)
    sigma = {}
    for m in (1.0, 0.0, -1.0):
        w, _, _ = kern.amplitudes(t, x, q2, s, EventSpinState(0, 0.0, 1.0, m))
        sigma[m] = 1.0 + w
    measured = ((sigma[1.0] + sigma[-1.0] - 2 * sigma[0.0])
                / (sigma[1.0] + sigma[-1.0] + sigma[0.0]))
    expected = azz(t["b1"], t["f1"], t["f2"], x, y)
    np.testing.assert_allclose(measured, expected, rtol=1e-12)


def test_transverse_tensor_matches_a_cos2phi(setup6):
    _, s = setup6
    kern = InclusiveKernel(beams.LI6,
                           delta_func=lambda x, q2, f1:
                           toy_delta_gluon(x, q2, f1, scale=1e-3))
    x, q2, y = _grid(s)
    t = kern.tables(x, q2)
    for m in (1.0, -1.0):
        w, _, a2 = kern.amplitudes(t, x, q2, s,
                                   EventSpinState(0, 0.0, 1.0, m,
                                                  theta_s=np.pi / 2))
        assert np.all(w == 0.0)  # no b1 -> no phi-averaged tensor shift
        expected = a_cos2phi(t["delta"], t["f1"], t["f2"], x, y)
        np.testing.assert_allclose(a2, expected, rtol=1e-12)


def test_dsigma_reduces_to_unpolarized(setup6):
    _, s = setup6
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1,
                           delta_func=toy_delta_gluon)
    x, q2, _ = _grid(s)
    phi = np.linspace(0, 2 * np.pi, 48, endpoint=False)
    # population-average over an unpolarized mixture = unpolarized xsec
    tot = 0.0
    for m in (1.0, 0.0, -1.0):
        state = EventSpinState(0, 0.0, 1.0, m, theta_s=0.7, phi_s=1.3)
        dsig = kern.dsigma(x[:, None], q2[:, None], phi[None, :], s, state)
        tot = tot + dsig.mean(axis=1) / 3.0
    np.testing.assert_allclose(tot * 2 * np.pi,
                               kern.dsigma_unpol(x, q2, s), rtol=1e-9)


def test_g2_ww_analytic_power_law():
    # g1 = x^a  ->  g2_WW = -x^a + (1 - x^a)/a
    a = 0.7
    x = np.array([0.01, 0.1, 0.4, 0.8])
    got = g2_ww(lambda xx, qq: np.power(xx, a), x, 10.0, npts=400)
    expected = -x**a + (1.0 - x**a) / a
    np.testing.assert_allclose(got, expected, rtol=2e-3)


def test_a_perp_transverse_amplitude(setup6):
    """a_perp against an independent transcription of the leading-gamma
    formula d(y) * gamma * ((y/2) g1 + g2) / F1 (E143 conventions)."""
    from polli_fastsim.asymmetries import depolarization_d
    _, s = setup6
    kern = InclusiveKernel(beams.LI6)  # g2 = g2_WW default
    x = np.full(3, 0.2)
    q2 = np.array([5.0, 20.0, 80.0])
    y = q2 / (s * x)
    t = kern.tables(x, q2, with_g2=True)
    state = EventSpinState(+1, 1.0, 1.0, 1.0, theta_s=np.pi / 2)
    _, a1, _ = kern.amplitudes(t, x, q2, s, state, with_perp=True)
    eps = (1.0 - y) / (1.0 - y + 0.5 * y * y)
    d = depolarization_d(y, x, q2) * np.sqrt(2.0 * eps / (1.0 + eps))
    gamma = 2.0 * 0.9383 * x / np.sqrt(q2)
    expected = d * gamma * (0.5 * y * t["g1"] + t["g2"]) / t["f1"]
    np.testing.assert_allclose(a1, expected, rtol=1e-12)
    assert np.all(expected != 0.0)
    # longitudinal axis -> no cos(phi) modulation
    _, a1_l, _ = kern.amplitudes(t, x, q2, s,
                                 EventSpinState(+1, 1.0, 1.0, 1.0),
                                 with_perp=True)
    assert np.all(a1_l == 0.0)


def test_g2_ww_broadcasting():
    def g1(xx, qq):
        return np.power(xx, 0.7)

    assert g2_ww(g1, 0.1, np.array([5.0, 10.0])).shape == (2,)
    grid = g2_ww(g1, np.array([[0.05], [0.2]]), np.array([[5.0, 20.0]]))
    assert grid.shape == (2, 2)
    np.testing.assert_allclose(
        grid[0], g2_ww(g1, np.array([0.05, 0.05]), np.array([5.0, 20.0])))


# --- grid backends (auto-skip when parton/grids absent) --------------------
# NOTE: must NOT use module-level importorskip -- it would silently skip the
# Gate-2 identity tests above on machines without the optional parton package.


def _have_grids():
    try:
        from parton import mkPDF
        mkPDF("CT18NLO", 0)
        mkPDF("NNPDFpol11_100", 0)
        return True
    except Exception:
        return False


@pytest.mark.skipif(not _have_grids(),
                    reason="parton or PDF grids not installed")
def test_identities_on_grid_backends(setup6):
    from polli_fastsim.polarized import PartonG1
    from polli_fastsim.structure import PartonF2
    _, s = setup6
    base = PartonF2()
    g1_model = PartonG1(base=base)
    x = np.array([3e-3, 0.03, 0.15, 0.45])
    q2 = np.array([4.0, 10.0, 30.0, 60.0])
    y = q2 / (s * x)
    # vector identity on a vector-only kernel (no tensor dilution)
    kern_v = InclusiveKernel(beams.LI6, f2_source=base, g1_model=g1_model)
    t = kern_v.tables(x, q2)
    wp, _, _ = kern_v.amplitudes(t, x, q2, s,
                                 EventSpinState(+1, 0.7, 1.0, 1.0))
    wm, _, _ = kern_v.amplitudes(t, x, q2, s,
                                 EventSpinState(-1, 0.7, 1.0, 1.0))
    np.testing.assert_allclose((wp - wm) / (2.0 + wp + wm),
                               0.7 * a_parallel(t["g1"], t["f1"], y, x, q2),
                               rtol=1e-12)
    # tensor identity with the scenario b1 on grid F1/F2
    kern = InclusiveKernel(beams.LI6, f2_source=base, g1_model=g1_model,
                           b1_func=toy_b1)
    t = kern.tables(x, q2)
    sigma = {}
    for m in (1.0, 0.0, -1.0):
        w, _, _ = kern.amplitudes(t, x, q2, s, EventSpinState(0, 0.0, 1.0, m))
        sigma[m] = 1.0 + w
    np.testing.assert_allclose(
        (sigma[1.0] + sigma[-1.0] - 2 * sigma[0.0])
        / (sigma[1.0] + sigma[-1.0] + sigma[0.0]),
        azz(t["b1"], t["f1"], t["f2"], x, y), rtol=1e-12)


# --- paths with no coverage until 2026-08-25 (plans/08 B4) -----------------

def test_b2_override_is_used_and_is_not_the_default():
    """`b2_func` has never been exercised: every test and script leaves it
    None, so the kernel silently uses b2 = 2x b1."""
    cfg = beams.default_configs("6Li")[1]
    s = cfg.sqrt_s_per_nucleon ** 2
    x = np.array([0.02, 0.056, 0.14])
    q2 = np.array([1.14, 1.14, 3.13])
    b1f = lambda x, q2, f1: 0.05 * f1                       # noqa: E731
    b2f = lambda x, q2, f1: 5.0 * x * 0.05 * f1             # noqa: E731
    default = InclusiveKernel(beams.LI6, b1_func=b1f)
    override = InclusiveKernel(beams.LI6, b1_func=b1f, b2_func=b2f)
    td, to = default.tables(x, q2), override.tables(x, q2)
    np.testing.assert_allclose(td["b2"], 2.0 * x * td["b1"], rtol=1e-12)
    np.testing.assert_allclose(to["b2"], 5.0 * x * to["b1"], rtol=1e-12)
    assert not np.allclose(td["b2"], to["b2"])
    state = EventSpinState(lam_e=0, pe=0.0, j=1.0, m=1.0, theta_s=0.0)
    y = q2 / (s * x)
    wo, _, _ = override.amplitudes(to, x, q2, s, state)
    np.testing.assert_allclose(
        wo, 0.5 * azz(to["b1"], to["f1"], to["f2"], x, y, b2=to["b2"]),
        rtol=1e-12)
    wd, _, _ = default.amplitudes(td, x, q2, s, state)
    assert not np.allclose(wd, wo)


def test_tensor_rate_at_an_intermediate_axis_angle():
    """theta_S is only ever 0 or pi/2 in the suite; the P_2(cos theta_S)
    geometry between them is untested."""
    cfg = beams.default_configs("6Li")[1]
    s = cfg.sqrt_s_per_nucleon ** 2
    x, q2 = np.array([0.056]), np.array([1.14])
    y = q2 / (s * x)
    kern = InclusiveKernel(beams.LI6, b1_func=lambda x, q2, f1: 0.05 * f1)
    t = kern.tables(x, q2)
    for th in (0.0, 0.4, 0.7, np.arccos(1.0 / np.sqrt(3.0)), np.pi / 2):
        state = EventSpinState(lam_e=0, pe=0.0, j=1.0, m=1.0, theta_s=th)
        w, _, _ = kern.amplitudes(t, x, q2, s, state)
        np.testing.assert_allclose(
            w, 0.5 * azz(t["b1"], t["f1"], t["f2"], x, y, b2=t["b2"],
                         theta_m=th), rtol=1e-12)
    # the magic angle kills the rate shift exactly
    magic = EventSpinState(lam_e=0, pe=0.0, j=1.0, m=1.0,
                           theta_s=np.arccos(1.0 / np.sqrt(3.0)))
    w, _, _ = kern.amplitudes(t, x, q2, s, magic)
    assert abs(float(w)) < 1e-15
