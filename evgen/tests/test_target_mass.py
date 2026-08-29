"""The gamma^2 (target-mass) term of the longitudinal double-spin kernel.

Everything published in this repository uses the massless
A_par = D(y) g1/F1; `InclusiveKernel(..., target_mass=True)` restores the
exact E143 form (PRD 58:112003) A_par = D_gamma (A1 + eta A2).  These
tests pin, in this order:

* the kinematic cap.  W^2 >= 10 GeV^2 (`fom.Scenario`) forces
  gamma^2 <= M^2/(W2_min - M^2) over the whole accepted phase space of
  every ion/energy configuration, so the term can never be large;
* the finite-gamma eps, D_gamma and eta against E143's OWN lab-frame
  definitions built from (E, E', theta) -- an independent construction,
  to double precision -- and the same for the assembled A_par;
* that switching the flag moves A_par by O(gamma^2) and by nothing else:
  it scales exactly with 1/Q^2 at fixed (x, y), it collapses to the
  g2-independent factor (1 + gamma^2) at small y, and with the flag off
  the kernel is bit-for-bit `asymmetries.a_parallel`;
* and the same target mass seen in the azimuth rather than in A_par --
  the mrad bounds Report 2 section 4.1 quotes for the lab-angle
  shortcut at the twelve published sweet spots.
"""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen.xsec import (EventSpinState, InclusiveKernel,  # noqa: E402
                           M_NUCLEON, depolarization_gamma, epsilon_gamma,
                           eta_gamma, gamma_squared)

from polli_fastsim import asymmetries, beams, fom  # noqa: E402
from polli_fastsim.asymmetries import (a_parallel,  # noqa: E402
                                       depolarization_d)
from polli_fastsim.polarized import toy_b1  # noqa: E402


def _kernels(ion=None):
    ion = ion or beams.LI6
    return InclusiveKernel(ion), InclusiveKernel(ion, target_mass=True)


def _dis_grid():
    """Flat (x, Q2) arrays over the DIS region W^2 >= 10 GeV^2."""
    x = np.logspace(-3, np.log10(0.7), 40)
    q2 = np.logspace(0.0, 2.0, 12)
    xx, qq = np.meshgrid(x, q2, indexing="ij")
    w2 = qq * (1.0 - xx) / xx + M_NUCLEON**2
    return xx[w2 >= 10.0], qq[w2 >= 10.0]


# --- (a) the kinematic cap ------------------------------------------------


def test_gamma2_capped_by_the_w2_cut_at_every_configuration():
    """gamma^2 = 4M^2x^2/Q^2 with Q^2 >= (W2min - M^2) x/(1-x) gives
    gamma^2 <= 4M^2 x(1-x)/(W2min - M^2) <= M^2/(W2min - M^2), maximal at
    x = 1/2.  Checked against the accepted 40x30 analysis grid of all six
    ion/energy settings; the numbers are quoted in Reports 0-2 and
    printed by evgen/scripts/target_mass_bound.py."""
    scenario = fom.Scenario()
    cap = M_NUCLEON**2 / (scenario.w2_min - M_NUCLEON**2)
    assert cap == pytest.approx(0.09654, rel=1e-3)
    seen = {}
    for ion in ("6Li", "7Li"):
        for cfg in beams.default_configs(ion):
            proj = fom.project_rates(cfg, scenario)
            g2v = np.where(proj.accepted, gamma_squared(proj.x, proj.q2), 0.0)
            # the cell-wise cap, tighter than the global one
            x = proj.x[proj.accepted]
            assert np.all(g2v[proj.accepted] <= 4.0 * M_NUCLEON**2 * x
                          * (1.0 - x) / (scenario.w2_min - M_NUCLEON**2)
                          * (1.0 + 1e-12))
            seen[cfg.label()] = g2v.max()
            assert g2v.max() < cap
    # the published per-configuration maxima (target_mass_bound.py block 1)
    got = sorted(round(v, 4) for v in seen.values())
    assert got == [0.0258, 0.0332, 0.0577, 0.0577, 0.0854, 0.0854]


# --- (b) the finite-gamma factors against E143's lab-frame forms ---------


def test_finite_gamma_factors_match_the_lab_frame_definitions():
    """eps = 1/[1 + 2(1 + nu^2/Q^2) tan^2(theta/2)],
    D = (1 - E' eps/E)/(1 + eps R), eta = eps sqrt(Q^2)/(E - E' eps)
    (E143 PRD 58:112003), from random (E, E', theta) in the target rest
    frame.  Nothing but double-precision agreement is acceptable: the
    two forms are the same object written in different variables."""
    rng = np.random.default_rng(20260828)
    r = 0.18
    for _ in range(64):
        e_in = rng.uniform(5.0, 50.0)
        e_out = e_in * rng.uniform(0.05, 0.95)
        theta = rng.uniform(0.05, 1.2)
        nu = e_in - e_out
        q2 = 4.0 * e_in * e_out * np.sin(theta / 2.0) ** 2
        x = q2 / (2.0 * M_NUCLEON * nu)
        y = nu / e_in
        g2v = gamma_squared(x, q2)
        assert g2v == pytest.approx(q2 / nu**2, rel=1e-12)
        eps_lab = 1.0 / (1.0 + 2.0 * (1.0 + nu * nu / q2)
                         * np.tan(theta / 2.0) ** 2)
        d_lab = (1.0 - e_out * eps_lab / e_in) / (1.0 + eps_lab * r)
        eta_lab = eps_lab * np.sqrt(q2) / (e_in - e_out * eps_lab)
        assert epsilon_gamma(y, g2v) == pytest.approx(eps_lab, rel=1e-12)
        assert depolarization_gamma(y, g2v, r) == pytest.approx(d_lab,
                                                                rel=1e-12)
        assert eta_gamma(y, g2v) == pytest.approx(eta_lab, rel=1e-12)


def test_target_mass_a_parallel_matches_a_lab_frame_construction():
    """The assembled A_par = D(A1 + eta A2) rebuilt from (E, E', theta)
    with gamma^2 = Q^2/nu^2 instead of 4M^2x^2/Q^2."""
    _, kern = _kernels()
    x, q2 = _dis_grid()
    x, q2 = x[::7], q2[::7]
    y = np.full(x.shape, 0.3)
    t = kern.tables(x, q2)
    nu = q2 / (2.0 * M_NUCLEON * x)
    e_in = nu / y
    e_out = e_in - nu
    sin2 = q2 / (4.0 * e_in * e_out)
    ok = (e_out > 0) & (sin2 > 0) & (sin2 < 1.0)
    assert ok.sum() > 20
    theta = 2.0 * np.arcsin(np.sqrt(sin2[ok]))
    eps = 1.0 / (1.0 + 2.0 * (1.0 + nu[ok] ** 2 / q2[ok])
                 * np.tan(theta / 2.0) ** 2)
    r = asymmetries.r_sigma_lt(x[ok], q2[ok])
    d = (1.0 - e_out[ok] * eps / e_in[ok]) / (1.0 + eps * r)
    eta = eps * np.sqrt(q2[ok]) / (e_in[ok] - e_out[ok] * eps)
    gam2 = q2[ok] / nu[ok] ** 2
    f1 = t["f1"][ok]
    a1 = (t["g1"][ok] - gam2 * t["g2"][ok]) / f1
    a2 = np.sqrt(gam2) * (t["g1"][ok] + t["g2"][ok]) / f1
    np.testing.assert_allclose(kern.a_parallel(t, x, q2, y)[ok],
                               d * (a1 + eta * a2), rtol=1e-12)


def test_massless_limit_of_the_finite_gamma_factors():
    """At gamma = 0 the finite-gamma set must BE the fast simulation's:
    eps -> (1-y)/(1-y+y^2/2), D -> asymmetries.depolarization_d, eta -> 0."""
    y = np.array([0.01, 0.05, 0.2, 0.5, 0.9])
    x, q2 = 0.1, 4.0
    r = asymmetries.r_sigma_lt(x, q2)
    np.testing.assert_allclose(epsilon_gamma(y, 0.0),
                               (1.0 - y) / (1.0 - y + 0.5 * y * y),
                               rtol=1e-14)
    np.testing.assert_allclose(depolarization_gamma(y, 0.0, r),
                               depolarization_d(y, x, q2), rtol=1e-14)
    assert np.all(eta_gamma(y, 0.0) == 0.0)


# --- (c) what the flag does, and does not, move -------------------------


def test_target_mass_off_is_bit_for_bit_the_published_kernel():
    kern0, kern1 = _kernels()
    x, q2 = _dis_grid()
    y = np.full(x.shape, 0.2)
    t0 = kern0.tables(x, q2)
    assert "g2" not in t0                      # nothing extra is computed
    t1 = kern1.tables(x, q2)
    assert "g2" in t1                          # ... and it IS when needed
    expected = a_parallel(t0["g1"], t0["f1"], y, x, q2)
    assert np.all(kern0.a_parallel(t0, x, q2, y) == expected)
    # the two kernels share every unpolarized/tensor number
    for key in ("f1", "f2", "g1", "b1", "b2", "delta"):
        np.testing.assert_array_equal(t0[key], t1[key])


def test_target_mass_needs_g2():
    with pytest.raises(ValueError):
        InclusiveKernel(beams.LI6, g2_mode="zero", target_mass=True)
    kern0, kern1 = _kernels()
    x, q2 = np.array([0.2]), np.array([5.0])
    with pytest.raises(KeyError):
        kern1.a_parallel(kern0.tables(x, q2), x, q2, np.array([0.3]))


def test_flag_moves_a_parallel_by_o_gamma_squared_only():
    """|A_par(tm)/A_par(0) - 1| stays within a few gamma^2 everywhere the
    W^2 cut allows, and scales EXACTLY as 1/Q^2 at fixed (x, y): a factor
    100 in Q^2 divides the shift by 100 to better than 5%."""
    kern0, kern1 = _kernels()
    x, q2 = _dis_grid()
    t = kern1.tables(x, q2)
    g2v = gamma_squared(x, q2)
    # away from the g1 zero crossing, where g2_WW/g1 is unbounded and any
    # ratio to A1 stops being a figure of merit
    tame = np.abs(t["g2"] / t["g1"]) <= 3.0
    for yy in np.linspace(0.02, 0.9, 12):
        y = np.full(x.shape, yy)
        shift = kern1.a_parallel(t, x, q2, y) / kern0.a_parallel(t, x, q2,
                                                                 y) - 1.0
        assert np.all(np.abs(shift[tame]) <= 3.0 * g2v[tame])
    for xx, yy in ((0.3, 0.1), (0.1, 0.05), (0.5, 0.3)):
        s = []
        for qq in (2.0, 200.0):
            xa, qa, ya = (np.array([xx]), np.array([qq]), np.array([yy]))
            tt = kern1.tables(xa, qa)
            ratio = (kern1.a_parallel(tt, xa, qa, ya)
                     / kern0.a_parallel(tt, xa, qa, ya) - 1.0)
            s.append(float(ratio[0]))
        assert s[0] / s[1] == pytest.approx(100.0, rel=0.05)


def test_small_y_collapse_to_one_plus_gamma_squared():
    """A_par(tm)/A_par(0) = (1 + gamma^2) + O(gamma^2 y), independent of
    g2: with kappa = eps y/[1-(1-y) eps] the exact ratio is
    (D_gamma/D)[1 + gamma^2 kappa + gamma^2 (kappa-1) g2/g1] and
    kappa - 1 = -0.5 y + O(y^2).  This is the form the reports quote."""
    kern0, kern1 = _kernels()
    x, q2 = _dis_grid()
    t = kern1.tables(x, q2)
    g2v = gamma_squared(x, q2)
    tame = np.abs(t["g2"] / t["g1"]) <= 3.0
    devs = []
    for yy in (0.05, 0.01):
        y = np.full(x.shape, yy)
        shift = kern1.a_parallel(t, x, q2, y) / kern0.a_parallel(t, x, q2,
                                                                 y) - 1.0
        devs.append(np.abs(shift[tame] / g2v[tame] - 1.0).max())
    assert devs[0] <= 0.08          # y <= 0.05: the collapse holds to 8%
    assert devs[1] <= 0.02          # and improves linearly in y
    assert devs[0] / devs[1] > 4.0
    # kappa -> 1 is the model-independent half of the statement
    for yy in (0.05, 0.01):
        for gg in (0.0, M_NUCLEON**2 / (10.0 - M_NUCLEON**2)):
            eps = epsilon_gamma(yy, gg)
            kappa = eps * yy / (1.0 - (1.0 - yy) * eps)
            assert abs(kappa - 1.0) < 0.6 * yy


def test_flag_leaves_the_tensor_and_unpolarized_sectors_alone():
    """Only the vector-L term moves: with an unpolarized electron the
    two kernels give the same b1 rate shift and the same cos 2phi."""
    cfg = beams.default_configs("6Li")[1]
    s = cfg.sqrt_s_per_nucleon**2
    kern0 = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    kern1 = InclusiveKernel(beams.LI6, b1_func=toy_b1, target_mass=True)
    x, q2 = _dis_grid()
    keep = q2 / (s * x) < 0.95
    x, q2 = x[keep], q2[keep]
    state = EventSpinState(0, 0.0, 1.0, 1.0, theta_s=0.4)
    got = [k.amplitudes(k.tables(x, q2), x, q2, s, state)
           for k in (kern0, kern1)]
    for a, b in zip(got[0], got[1]):
        np.testing.assert_array_equal(a, b)
    # ... and the vector-L term does move, by the a_parallel shift and
    # nothing else (the b1 rate term is common to w0 and w1)
    lam = EventSpinState(+1, 0.7, 1.0, 1.0)
    t0, t1 = kern0.tables(x, q2), kern1.tables(x, q2)
    w0 = kern0.amplitudes(t0, x, q2, s, lam)[0]
    w1 = kern1.amplitudes(t1, x, q2, s, lam)[0]
    assert np.any(w0 != w1)
    y = q2 / (s * x)
    # (rtol only: w carries the b1 rate term, so the O(1e-9) vector
    # difference is taken between numbers ~1e-2 and loses eight digits)
    np.testing.assert_allclose(
        w1 - w0, 0.7 * (kern1.a_parallel(t1, x, q2, y)
                        - kern0.a_parallel(t0, x, q2, y)),
        rtol=1e-6, atol=1e-15)


def test_azimuth_shortcut_error_at_the_published_sweet_spots():
    """The same target-mass effect seen in the azimuth rather than in
    A_par.  Report 2 section 4.1 quotes the error the lab-angle shortcut
    phi_S = phi_e - phi_s would make on the cos 2phi observable at the
    twelve money-plot-5 sweet spots; this pins those three numbers, and
    that the error is gamma^2/4 rather than something with its own
    kinematics.  (Exactness for a massless target is
    tests/test_reco.py's.)"""
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]
                           / "scripts"))
    import target_mass_bound as tmb
    quoted = {0: 6.0, 1: 1.4, 2: 1.3}          # mrad, low / mid / top
    for ci, bound in quoted.items():
        cfg = beams.default_configs("6Li")[ci]
        s = cfg.sqrt_s_per_nucleon**2
        spots = tmb.published_spots(cfg)
        x = np.array([sp[0] for sp in spots])
        q2 = np.array([sp[1] for sp in spots])
        dphi = tmb.azimuth_shortcut_error(cfg, x, q2, q2 / (s * x))
        assert dphi.max() <= bound             # the report's bound holds
        assert dphi.max() > 0.9 * bound        # ... and is not loose
        # O(gamma^2) with the coefficient 1/4, not a shape of its own
        ratio = 1e-3 * dphi / gamma_squared(x, q2)
        assert np.all((ratio > 0.20) & (ratio < 0.26))
