"""Gate 1 (plans/05 §5.4): rho-matrix moments, all axes — exact."""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen import spin  # noqa: E402


def test_wigner_d_spin1_known_matrix():
    beta = 0.7
    d = spin.wigner_d(1.0, beta)
    c, s = np.cos(beta), np.sin(beta)
    # rows/cols ordered m = +1, 0, -1
    expected = np.array([
        [(1 + c) / 2, -s / np.sqrt(2), (1 - c) / 2],
        [s / np.sqrt(2), c, -s / np.sqrt(2)],
        [(1 - c) / 2, s / np.sqrt(2), (1 + c) / 2],
    ])
    np.testing.assert_allclose(d, expected, atol=1e-12)


def test_wigner_d_spin32_diagonal_entries():
    beta = 1.1
    d = spin.wigner_d(1.5, beta)
    ch = np.cos(beta / 2)
    np.testing.assert_allclose(d[0, 0], ch**3, atol=1e-12)
    np.testing.assert_allclose(d[1, 1], (3 * np.cos(beta) - 1) / 2 * ch,
                               atol=1e-12)


@pytest.mark.parametrize("j", [0.5, 1.0, 1.5])
def test_rotation_matrix_unitary(j):
    u = spin.rotation_matrix(j, 0.9, 2.3)
    np.testing.assert_allclose(u @ u.conj().T, np.eye(u.shape[0]), atol=1e-12)


def test_clebsch_gordan_reference_values():
    assert np.isclose(spin.clebsch_gordan(1, 1, 1, -1, 0, 0), 1 / np.sqrt(3))
    assert np.isclose(spin.clebsch_gordan(0.5, 0.5, 0.5, -0.5, 1, 0),
                      1 / np.sqrt(2))
    assert np.isclose(spin.clebsch_gordan(1, 0, 1, 0, 2, 0), np.sqrt(2 / 3))


@pytest.mark.parametrize("j", [1.0, 1.5])
def test_multipole_operators_orthonormal(j):
    ks = range(0, int(round(2 * j)) + 1)
    ops = {(k, q): spin.multipole_operator(j, k, q)
           for k in ks for q in range(-k, k + 1)}
    for (k1, q1), t1 in ops.items():
        for (k2, q2), t2 in ops.items():
            val = np.trace(t1.conj().T @ t2)
            expected = 1.0 if (k1, q1) == (k2, q2) else 0.0
            np.testing.assert_allclose(val, expected, atol=1e-12)


def test_pure_and_unpolarized_moments():
    # pure m=+1 along z
    rho = spin.rho_from_populations(1.0, (1.0, 0.0, 0.0))
    np.testing.assert_allclose(spin.vector_polarization(rho, 1.0), [0, 0, 1],
                               atol=1e-12)
    assert np.isclose(spin.tensor_polarization(rho, 1.0), 1.0)
    # unpolarized: all moments vanish
    for j, n in ((1.0, 3), (1.5, 4)):
        rho = spin.rho_from_populations(j, tuple([1.0 / n] * n), 0.4, 1.2)
        np.testing.assert_allclose(spin.vector_polarization(rho, j),
                                   [0, 0, 0], atol=1e-12)
        assert np.isclose(spin.tensor_polarization(rho, j), 0.0, atol=1e-12)
    # HERMES-style m0-enriched fill: Pzz = -2
    rho = spin.rho_from_populations(1.0, (0.0, 1.0, 0.0))
    assert np.isclose(spin.tensor_polarization(rho, 1.0), -2.0)


def test_spin1_population_round_trip():
    pz, pzz = 0.5, 0.3
    pops = spin.spin1_populations(pz, pzz)
    v, t = spin.moments_along_axis(1.0, pops)
    assert np.isclose(v, pz) and np.isclose(t, pzz)
    with pytest.raises(ValueError):
        spin.spin1_populations(1.0, -2.0)  # unphysical corner


def test_spin32_population_round_trip():
    pz, t, o = 0.6, 0.2, -0.1
    pops = spin.spin32_populations(pz, t, o)
    v, tt, oo = spin.moments_along_axis(1.5, pops)
    np.testing.assert_allclose([v, tt, oo], [pz, t, o], atol=1e-12)
    # pure m=+3/2 has all normalized moments = +1
    np.testing.assert_allclose(spin.spin32_populations(1.0, 1.0, 1.0),
                               (1.0, 0.0, 0.0, 0.0), atol=1e-12)


@pytest.mark.parametrize("j,pops", [
    (1.0, (0.6, 0.3, 0.1)),
    (1.5, (0.5, 0.25, 0.15, 0.1)),
])
def test_rotated_moments_analytic(j, pops):
    theta, phi = 0.8, 2.1
    dens = spin.SpinDensity(j, pops, theta, phi)
    mom = dens.lab_moments()
    axis_moments = spin.moments_along_axis(j, pops)
    v, t = axis_moments[0], axis_moments[1]
    n_hat = np.array([np.sin(theta) * np.cos(phi),
                      np.sin(theta) * np.sin(phi), np.cos(theta)])
    np.testing.assert_allclose(mom["vector"], v * n_hat, atol=1e-12)
    # rank-2 moment transforms with P2(cos theta) along z
    p2 = 0.5 * (3 * np.cos(theta) ** 2 - 1)
    np.testing.assert_allclose(mom["tensor_zz"], t * p2, atol=1e-12)


def test_unphysical_populations_raise():
    with pytest.raises(ValueError):
        spin.rho_from_populations(1.0, (0.7, 0.6, -0.3))
    with pytest.raises(ValueError):
        spin.spin32_populations(1.2, 0.0, 0.0)
