"""The gamma^2 (target-mass) term of the longitudinal double-spin kernel.

Since 2026-08-29 the exact E143 form (PRD 58:112003)
A_par = D_gamma (A1 + eta A2) is the DEFAULT
(`InclusiveKernel(..., target_mass=True)`) and the massless
A_par = D(y) g1/F1 every figure published before that date was made on is
`target_mass=False`, kept reachable and pinned here bit-for-bit.  These
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
* the extraction that inverts it: A_par = D_eff (g1/F1) with
  D_eff = `asymmetries.depolarization_effective`, so that
  `fom.project_observables` returns delta(g1/F1) with no O(gamma^2) bias
  left in it, and a closure on the toy that recovers an injected
  multiplicative modification of g1 exactly;
* the residual that survives -- the twist-3 uncertainty on g2, whose
  span between g2 = 0 and 1.5 g2^WW is what `target_mass_bound.py` now
  quotes as the systematic in place of the removed bias;
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
                           eta_gamma, g2_ww, gamma_squared)

from polli_fastsim import asymmetries, beams, fom  # noqa: E402
from polli_fastsim.asymmetries import (a_parallel,  # noqa: E402
                                       a_parallel_exact,
                                       depolarization_d,
                                       depolarization_effective)
from polli_fastsim.polarized import ToyG1, toy_b1, toy_delta_gluon  # noqa: E402


def _kernels(ion=None):
    """(massless, finite-gamma) -- the second one is the DEFAULT kernel."""
    ion = ion or beams.LI6
    return InclusiveKernel(ion, target_mass=False), InclusiveKernel(ion)


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
    """A_par and the unpolarized/tensor tables are bit-for-bit.

    g2 and a_perp are NOT, at the 1e-14 level, and deliberately so: g2
    now comes from the backend's cached `g2_nucleus(...)/A` rather than
    from `g2_ww` on the already-divided g1, so the division by A crossed
    the quadrature.  g2^WW is linear in g1, so the two orderings differ
    only by the rounding of the 96-node sum; the tolerance below is what
    pins that they differ by NOTHING ELSE."""
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
    # g2 against the OLD ordering, g2_ww on the per-nucleon g1
    old_order = g2_ww(kern1._g1a, x, q2)
    scale = np.abs(old_order).max()
    assert np.any(t1["g2"] != old_order)       # the reorder is visible ...
    assert np.abs(t1["g2"] - old_order).max() < 1e-13 * scale   # ... and tiny


def test_target_mass_needs_g2_in_the_tables_it_is_handed():
    """A finite-gamma kernel cannot use a massless kernel's tables.

    g2_mode="zero" is NOT refused any more: with target_mass on it is
    exactly the g2 = 0 end of the twist-3 variation, which is the
    residual systematic the reports now quote."""
    kern0, kern1 = _kernels()
    x, q2 = np.array([0.2]), np.array([5.0])
    with pytest.raises(KeyError):
        kern1.a_parallel(kern0.tables(x, q2), x, q2, np.array([0.3]))
    zero = InclusiveKernel(beams.LI6, g2_mode="zero")
    tz = zero.tables(x, q2)
    assert np.all(tz["g2"] == 0.0)
    assert np.isfinite(zero.a_parallel(tz, x, q2, np.array([0.3]))).all()


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
    kern0 = InclusiveKernel(beams.LI6, b1_func=toy_b1, target_mass=False)
    kern1 = InclusiveKernel(beams.LI6, b1_func=toy_b1)
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


# --- (d) the extraction that inverts the term ---------------------------


def test_depolarization_effective_is_exactly_the_a_parallel_divisor():
    """A_par = D_eff (g1/F1) identically, cell by cell, with
    D_eff = D_gamma [1 - gamma^2 rho + eta gamma (1 + rho)] and
    rho = g2/g1.  This is the whole algebra of the extraction: what
    multiplies g1/F1 in the asymmetry is what delta(A_par) must be divided
    by to give delta(g1/F1)."""
    _, kern = _kernels()
    x, q2 = _dis_grid()
    for yy in (0.05, 0.2, 0.6):
        y = np.full(x.shape, yy)
        t = kern.tables(x, q2)
        rho = t["g2"] / t["g1"]
        np.testing.assert_allclose(
            kern.a_parallel(t, x, q2, y),
            depolarization_effective(y, x, q2, g2_over_g1=rho)
            * t["g1"] / t["f1"], rtol=1e-11)
    # ... and with no g2 it is bit-for-bit the massless divisor
    y = np.full(x.shape, 0.3)
    assert np.all(depolarization_effective(y, x, q2)
                  == depolarization_d(y, x, q2))


def test_the_massless_inversion_was_high_by_one_plus_gamma_squared():
    """The bias the flip removed, stated as the ratio of the two
    divisors: dividing the exact A_par by the massless D returns g1/F1
    high by (1 + gamma^2) + O(gamma^2 y), so D_eff/D - 1 IS gamma^2 at
    small y.  Same statement as the small-y collapse above, seen from the
    extraction rather than from the asymmetry."""
    _, kern = _kernels()
    x, q2 = _dis_grid()
    t = kern.tables(x, q2)
    g2v = gamma_squared(x, q2)
    rho = t["g2"] / t["g1"]
    tame = np.abs(rho) <= 3.0
    for yy, tol in ((0.05, 0.08), (0.01, 0.02)):
        y = np.full(x.shape, yy)
        ratio = (depolarization_effective(y, x, q2, g2_over_g1=rho)
                 / depolarization_d(y, x, q2) - 1.0)
        assert np.abs(ratio[tame] / g2v[tame] - 1.0).max() <= tol


def test_g2_ww_is_linear_in_g1_so_rho_is_a_property_of_the_shape():
    """The step the extraction rests on: g2^WW[c g1] = c g2^WW[g1], so
    the multiplicative medium modification Delta-R measures leaves
    rho = g2/g1 alone and D_eff is the same for model and measurement."""
    from polli_fastsim.polarized import g2_ww
    x = np.logspace(-2.5, np.log10(0.8), 25)
    q2 = np.full(x.shape, 8.0)

    def g1(xx, qq):
        return ToyG1().g1_nucleus(beams.LI7, xx, qq)

    for c in (0.5, 1.3):
        np.testing.assert_allclose(
            g2_ww(lambda xx, qq: c * g1(xx, qq), x, q2),
            c * g2_ww(g1, x, q2), rtol=1e-13)


def test_extraction_closure_recovers_an_injected_modification():
    """Inject Delta-R = c on g1, measure the exact A_par, invert it with
    D_eff built from the model's own rho: the ratio comes back as c to
    machine precision at every cell, for every c.  With the massless D it
    comes back as c(1 + gamma^2) -- the bias this stream removed."""
    x, q2 = _dis_grid()
    y = np.full(x.shape, 0.15)
    model = ToyG1()
    g1 = model.g1_nucleus(beams.LI7, x, q2) / beams.LI7.A
    g2 = model.g2_nucleus(beams.LI7, x, q2) / beams.LI7.A
    f1 = 1.0 / (2.0 * x) * np.ones_like(x)     # any positive F1 will do
    rho = g2 / g1
    g2v = gamma_squared(x, q2)
    d_eff = depolarization_effective(y, x, q2, g2_over_g1=rho)
    for c in (0.85, 1.0, 1.12):
        apar = a_parallel_exact(c * g1, c * g2, f1, y, x, q2)
        np.testing.assert_allclose(apar / d_eff / (g1 / f1), c, rtol=1e-11)
        naive = apar / depolarization_d(y, x, q2) / (g1 / f1)
        np.testing.assert_allclose(naive / c - 1.0, g2v, rtol=0.10)


def test_fom_extraction_is_self_consistent_and_moved_the_published_error():
    """`fom.project_observables` must invert the very A_par it returns:
    err_g1_over_f1 * D_eff == err_a_par exactly wherever D_eff > 0.  And
    the change is not cosmetic -- the published delta(g1/F1) moved by
    (1 + gamma^2), i.e. down, cell by cell.

    D_eff changes sign off the physical region (y > 1 here, 525 of the
    1200 cells of this map, none of them accepted) and at a zero crossing
    of g1 on a polarized grid, and a negative statistical error is not a
    statement.  Since 2026-08-29 the divisor is guarded and the error is
    +inf there, which is the honest reading -- A_par carries no
    information about g1/F1 where its own divisor vanishes -- and which
    drops out of every inverse-variance weight on its own."""
    cfg = beams.default_configs("7Li")[1]
    scenario = fom.Scenario()
    proj = fom.project_rates(cfg, scenario)
    model = ToyG1()
    obs = fom.project_observables(cfg, scenario, proj, model, toy_b1,
                                  toy_delta_gluon)
    X, Q2 = proj.x, proj.q2
    y = proj.extras["y"]
    f1 = proj.extras["nf2"].f1a(X, Q2) / cfg.ion.A
    g1 = model.g1_nucleus(cfg.ion, X, Q2) / cfg.ion.A
    g2 = model.g2_nucleus(cfg.ion, X, Q2) / cfg.ion.A
    rho = g2 / g1
    d_eff = depolarization_effective(y, X, Q2, g2_over_g1=rho)
    ok = d_eff > 0.0
    np.testing.assert_allclose((obs["err_g1_over_f1"] * d_eff)[ok],
                               obs["err_a_par"][ok], rtol=1e-12)
    assert np.all(np.isinf(obs["err_g1_over_f1"][~ok]))    # guarded, not negative
    assert np.all(obs["err_g1_over_f1"] > 0.0)             # never negative
    assert not np.any(proj.accepted & ~ok)                 # and never accepted
    np.testing.assert_allclose(obs["a_par"],
                               a_parallel_exact(g1, g2, f1, y, X, Q2),
                               rtol=1e-12)
    old = obs["err_a_par"] / depolarization_d(y, X, Q2)
    acc = proj.accepted
    shift = (old / obs["err_g1_over_f1"] - 1.0)[acc]
    # O(gamma^2) everywhere on the accepted grid -- exactly gamma^2 only
    # at the small y where the collapse holds, which the accepted grid
    # does not confine itself to
    assert np.all(np.abs(shift) <= 3.0 * gamma_squared(X, Q2)[acc])
    assert shift.max() > 1e-3          # the move is real, not round-off


def test_the_residual_is_the_twist_three_uncertainty_on_g2():
    """What survives the flip.  Extracting with an assumed g2 = s g2^WW
    while the truth is g2^WW mis-scales g1/F1 by C(rho)/C(s rho), which
    vanishes at s = 1, is linear in (s - 1) for small gamma^2, and is
    bounded by 2 gamma^2 |rho| |s - 1|.  s = 0 and s = 1.5 are the two
    variations `target_mass_bound.py` quotes the systematic from."""
    _, kern = _kernels(beams.LI7)
    x, q2 = _dis_grid()
    y = np.full(x.shape, 0.2)
    t = kern.tables(x, q2)
    rho = t["g2"] / t["g1"]
    g2v = gamma_squared(x, q2)
    tame = np.abs(rho) <= 3.0
    truth = depolarization_effective(y, x, q2, g2_over_g1=rho)
    got = {}
    for s in (0.0, 1.0, 1.5):
        assumed = depolarization_effective(y, x, q2, g2_over_g1=s * rho)
        got[s] = (truth / assumed - 1.0)
        assert np.all(np.abs(got[s][tame])
                      <= 2.0 * g2v[tame] * np.abs(rho[tame]) * abs(s - 1.0)
                      + 1e-12)
    assert np.all(got[1.0] == 0.0)
    # linear in (s - 1): the 1.5 variation is half the 0 one, up to O(g^4)
    ratio = np.abs(got[1.5][tame] / got[0.0][tame])
    assert np.abs(ratio - 0.5).max() < 0.05


def test_g2_scale_is_the_kernel_level_twist_three_knob():
    """`InclusiveKernel(g2_scale=s)` is the s of `g2_residual`.

    The script forms the variation analytically, from rho = g2/g1; a
    sampler cannot, so the kernel carries the same knob.  Pinned here
    because nothing else in the repository exercises it: the table scales
    exactly, s = 0 is g2_mode="zero", and the A_par it produces is the
    one the analytic residual predicts, to double precision."""
    x, q2 = _dis_grid()
    y = np.full(x.shape, 0.2)
    base = InclusiveKernel(beams.LI6)
    t1 = base.tables(x, q2)
    for s in (0.0, 0.5, 1.5):
        kern = InclusiveKernel(beams.LI6, g2_scale=s)
        ts = kern.tables(x, q2)
        np.testing.assert_allclose(ts["g2"], s * t1["g2"],
                                   rtol=1e-15, atol=0.0)
        # A_par at scale s, against D_eff(s rho) (g1/F1) -- the divisor
        # `target_mass_bound.g2_residual` varies instead of the kernel
        rho = t1["g2"] / t1["g1"]
        want = (depolarization_effective(y, x, q2, g2_over_g1=s * rho)
                * t1["g1"] / t1["f1"])
        np.testing.assert_allclose(kern.a_parallel(ts, x, q2, y), want,
                                   rtol=1e-12, atol=0.0)
    zero = InclusiveKernel(beams.LI6, g2_mode="zero")
    np.testing.assert_array_equal(
        zero.tables(x, q2)["g2"],
        InclusiveKernel(beams.LI6, g2_scale=0.0).tables(x, q2)["g2"])


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
