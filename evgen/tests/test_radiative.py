"""Collinear initial-state radiation as a migration (polligen.radiative,
plans/07 WP4 / plans/08 D3): the leading-log spectrum against its closed
forms, the observed kinematics of every reconstruction method against a
four-vector construction, the RecoResponse hook's default-off
immutability, and the migration bound at reduced statistics."""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen import bookkeeping as bk  # noqa: E402
from polligen import hfs, radiative as rad, reco, recopseudo as rp  # noqa: E402
from polligen.sample import InclusiveSampler  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.kinematics import scattered_electron  # noqa: E402
from polli_fastsim.polarized import toy_b1, toy_delta_gluon  # noqa: E402

CONFIG = beams.default_configs("6Li")[1]
S_NN = CONFIG.sqrt_s_per_nucleon ** 2


# --- the leading-log spectrum -------------------------------------------------

@pytest.mark.parametrize("q2", [1.0, 3.13, 14.3, 200.0])
def test_spectrum_normalisation_is_second_order_in_t(q2):
    """int_0^1 D dz = 1 + O(t^2): the truncated exponentiation's own
    error, and the reason the module can call itself a bound."""
    t = rad.beta_ll(q2)
    resid = rad.total_weight(q2)
    assert abs(resid) < t ** 2
    assert abs(resid) < 2.0e-3


def test_beta_and_soft_factor_are_the_textbook_expansions():
    """t = (2a/pi)(L - 1) and S(t) = 1 + 3t/8 + O(t^2)."""
    q2 = 3.13
    t = rad.beta_ll(q2)
    lhs = (2.0 / 137.035999 / np.pi) * (np.log(q2 / rad.M_E ** 2) - 1.0)
    assert t == pytest.approx(lhs, rel=2e-5)
    assert rad.soft_factor(t) == pytest.approx(1.0 + 3.0 * t / 8.0, abs=t ** 2)


def test_cdf_is_the_integral_of_the_density():
    """F is the primitive of D.  Quadrature away from the z -> 0 endpoint,
    where the z^(t/2-1) singularity defeats the trapezoid rule but the
    closed form is exact by construction."""
    q2 = 3.13
    z = np.linspace(1e-2, 1.0, 200001)
    num = float(np.trapz(rad.d_electron(z, q2), z))
    assert num == pytest.approx(float(rad.cdf_z(1.0, q2)
                                      - rad.cdf_z(1e-2, q2)), rel=1e-6)
    # and moment_z(order=0) is the same primitive
    assert float(rad.moment_z(q2, 1.0, 1e-2, order=0)) == pytest.approx(num,
                                                                        rel=1e-6)


def test_sampler_reproduces_the_closed_form_moments():
    """The rejection sampler is the density it advertises: <z> and <z^2>
    over the emitting events match moment_z / (F(zmax) - F(zmin))."""
    q2 = np.full(400000, 3.13)
    y = np.full(400000, 0.05)
    model = rad.ISRModel(seed=4242)
    z = model.sample(q2, y)
    zmax = float(model.z_max(0.05))
    norm = float(rad.moment_z(3.13, zmax, model.z_min, order=0))
    emit = z > 0.0
    n = z.size
    for order, tol in ((1, 4e-2), (2, 1.2e-1)):
        want = float(rad.moment_z(3.13, zmax, model.z_min, order=order)) / norm
        got = float((z[emit] ** order).sum() / emit.sum())
        assert got == pytest.approx(want, rel=tol)
    # and the emission probability is the CDF ratio
    assert emit.mean() == pytest.approx(float(model.p_emit(3.13, 0.05)),
                                        abs=4.0 / np.sqrt(n))


def test_kinematic_ceiling_caps_the_radiated_fraction():
    """y_hard = y/(1 - z) <= y_max_hard is a limit, not a cut."""
    model = rad.ISRModel(seed=7)
    y = np.repeat([0.02, 0.4, 0.9, 0.999], 20000)
    z = model.sample(np.full(y.size, 5.0), y)
    assert np.all(y / (1.0 - z) <= model.y_max_hard + 1e-12)
    assert np.all(z[y >= 0.999] == 0.0)
    # more phase space at low y -> more emission
    assert (z[y == 0.02] > 0).mean() > (z[y == 0.9] > 0).mean()


# --- what each reconstruction method sees ------------------------------------

def _four_vector_reference(x, y, z):
    """The observed kinematics of a radiative event built from actual
    four-vectors: the hard event at the REDUCED beam, then reconstructed
    with the NOMINAL one (Sigma_h = 2 E_hard y_hard and pT_h = pT_e for a
    massless target, hfs.truth_kinematics_check)."""
    q2 = x * y * S_NN
    one_z = 1.0 - z
    e_hard = CONFIG.electron_energy * one_z
    e_p, theta, _eta = scattered_electron(x, y / one_z, S_NN * one_z, e_hard)
    q2_e, y_e, x_e = reco.electron_method(e_p, theta, CONFIG.electron_energy,
                                          S_NN)
    sigma_h = 2.0 * e_hard * (y / one_z)
    pt = e_p * np.sin(theta)
    had = hfs.hadronic_kinematics(sigma_h, pt, np.zeros_like(pt), e_p, theta,
                                  CONFIG.electron_energy, S_NN)
    had.update(q2_e=q2_e, y_e=y_e, x_e=x_e, q2=q2)
    return had


def test_observed_kinematics_match_a_four_vector_construction():
    """Every closed form of `observed_kinematics` is the number the full
    construction plus `hfs.hadronic_kinematics` returns."""
    x = np.array([0.005, 0.02, 0.05, 0.1, 0.005, 0.02])
    y = np.array([0.02, 0.05, 0.1, 0.3, 0.5, 0.8])
    z = np.array([0.0, 0.01, 0.05, 0.2, 0.4, 0.05])
    ref = _four_vector_reference(x, y, z)
    got = rad.observed_kinematics(x, x * y * S_NN, y, z, S_NN)
    for key in ("q2_e", "y_e", "x_e", "q2_sigma", "y_sigma", "x_sigma",
                "q2_jb", "y_jb", "x_jb", "q2_da", "y_da", "x_da", "x_mixed"):
        assert got[key] == pytest.approx(ref[key], rel=1e-10), key


def test_mixed_x_is_exact_and_only_the_q2_label_migrates():
    """The plans/08 D3 statement, pinned: with a collinear photon the
    chain's x = Q2_e/(s y_Sigma) is EXACT because the (1 - z) cancels
    between Q2_e and y_Sigma, and it is the Q2_e label that migrates by
    1/(1 - z).  (The 'x -> x/(1 - z)' reading this replaced was wrong;
    code review R16.)  The Sigma-method Q2 uses no beam energy and is
    exact too, so an e-Sigma label would not migrate at all; the
    electron-only and double-angle x do migrate."""
    x = np.array([0.005, 0.02, 0.05, 0.1])
    y = np.array([0.02, 0.05, 0.1, 0.3])
    z = np.array([0.02, 0.1, 0.3, 0.5])
    ref = _four_vector_reference(x, y, z)
    x_hard, one_z = x, 1.0 - z
    q2_hard = x * y * S_NN
    assert ref["x_mixed"] == pytest.approx(x_hard, rel=1e-12)
    assert ref["q2_e"] == pytest.approx(q2_hard / one_z, rel=1e-12)
    assert ref["q2_sigma"] == pytest.approx(q2_hard, rel=1e-12)
    assert ref["y_sigma"] == pytest.approx(y / one_z, rel=1e-12)
    assert ref["y_da"] == pytest.approx(y / one_z, rel=1e-12)
    assert ref["y_jb"] == pytest.approx(y, rel=1e-12)
    assert ref["q2_da"] == pytest.approx(q2_hard / one_z ** 2, rel=1e-12)
    assert ref["x_da"] == pytest.approx(x_hard / one_z, rel=1e-12)
    assert ref["x_sigma"] == pytest.approx(x_hard * one_z, rel=1e-12)
    assert ref["x_e"] / x_hard == pytest.approx(y / (one_z * (y + z)),
                                                rel=1e-12)


def test_method_bias_table_rows():
    rows = dict((r[0], r[1:]) for r in rad.method_bias_table(0.05, 0.1))
    assert rows["mixed (Q2_e, y_Sigma)"][2] == pytest.approx(1.0, abs=1e-12)
    assert rows["mixed (Q2_e, y_Sigma)"][0] == pytest.approx(1.0 / 0.9,
                                                             rel=1e-12)
    assert rows["Sigma"][0] == pytest.approx(1.0, abs=1e-12)
    assert rows["double angle"][0] == pytest.approx(1.0 / 0.81, rel=1e-12)


# --- the RecoResponse hook ---------------------------------------------------

@pytest.fixture(scope="module")
def sampler():
    analysis = fom.Scenario(q2_min=1.0, y_min=0.01, y_max=0.95)
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1,
                           delta_func=lambda x, q2, f1:
                           toy_delta_gluon(x, q2, f1, scale=3e-2))
    gen = rp.generator_scenario(analysis)
    return InclusiveSampler(kern, CONFIG, gen, nx=28, nq2=22,
                            q2_range=(0.7, 3e2))


def _response(sampler, isr=None, seed=20260824, n_mc=120):
    return rp.RecoResponse(sampler, rp.RecoModel(), n_mc_per_cell=n_mc,
                           rng=np.random.default_rng(seed), isr=isr)


def _grid_positions(sampler, seed=20260824, n_mc=120):
    """(cell, x, q2, y) of an ISR-OFF response, rebuilt HERE from the
    sampler grid and a fresh generator alone -- the first two uniforms of
    `RecoResponse.__init__`, the log-uniform map, the generator
    acceptance and y = Q2/(s x).  Nothing from the cross section and
    nothing from reco.py enters, so this is an independent construction
    of the four arrays and not a second call to the code under test.

    It replaces a stored sha256 over the raw IEEE bytes (until
    2026-08-28), which asserted bitwise float equality across machines
    and NumPy versions -- the one thing the manual's Section 1 says the
    repository does NOT promise -- and would have failed on a NumPy
    upgrade with a message telling the reader to rebuild every published
    figure."""
    rng = np.random.default_rng(seed)
    ncell = sampler.xsec_flat.size
    cell = np.repeat(np.arange(ncell), n_mc)
    u1 = rng.uniform(size=cell.size)
    u2 = rng.uniform(size=cell.size)
    x = np.exp(sampler.logx_lo[cell]
               + u1 * (sampler.logx_hi[cell] - sampler.logx_lo[cell]))
    q2 = np.exp(sampler.logq2_lo[cell]
                + u2 * (sampler.logq2_hi[cell] - sampler.logq2_lo[cell]))
    ok = sampler._in_acceptance(x, q2)
    cell, x, q2 = cell[ok], x[ok], q2[ok]
    return cell, x, q2, q2 / (sampler.s * x)


def test_isr_off_path_is_the_pre_hook_constructor(sampler):
    """Every published number reproduces bit for bit with the hook in
    place.  Nothing here compares two identical calls: each quantity the
    hook could have touched is pinned against an INDEPENDENT value --
    a stored digest for the arrays the response rng and the sampler grid
    fix between them, and closed forms for the rest."""
    a = _response(sampler)

    # (i) the hook draws nothing from the RESPONSE stream: an ISR-on
    # response leaves the very same generator in the very same state.
    rng_off = np.random.default_rng(20260824)
    off = rp.RecoResponse(sampler, rp.RecoModel(), n_mc_per_cell=120,
                          rng=rng_off)
    state_off = rng_off.bit_generator.state
    rng_on = np.random.default_rng(20260824)
    rp.RecoResponse(sampler, rp.RecoModel(), n_mc_per_cell=120, rng=rng_on,
                    isr=rad.ISRModel(seed=99))
    assert rng_on.bit_generator.state == state_off

    # (ii) the in-cell positions are the grid's own, rebuilt here from a
    # fresh generator: the hook consumes no draw and rescales no position
    ref_cell, ref_x, ref_q2, ref_y = _grid_positions(sampler)
    for resp in (a, off):
        assert np.array_equal(resp.cell, ref_cell)
        assert np.array_equal(resp.x, ref_x)
        assert np.array_equal(resp.q2, ref_q2)
        assert np.array_equal(resp.y, ref_y)

    # (iii) the weight carries NO radiative rescaling: it is exactly the
    # cell cross section over the accepted count, as before the hook
    n_kept = np.bincount(a.cell, minlength=sampler.xsec_flat.size)
    assert np.array_equal(a.w, sampler.xsec_flat[a.cell]
                          / np.maximum(n_kept[a.cell], 1))

    # (iv) y is the nominal q2/(s x), the hard-y reading is inert, and
    # nothing multiplies the amplitude or the dilution
    assert np.array_equal(a.y, a.q2 / (sampler.s * a.x))
    assert np.array_equal(a.y_nominal, a.y)
    assert a.isr is None and a._amp_isr is None
    assert np.all(a.isr_z == 0.0)
    assert a._dil_amp() is a.dil
    cat = bk.tensor_flip_plan(0.6).categories[0]
    _den, num = a._fill_arrays(cat)
    assert np.array_equal(a.amplitude_per_event(cat),
                          num[a.cell] / float(cat.moments()[1]))


def test_isr_draws_from_its_own_stream_so_the_pair_shares_random_numbers(
        sampler):
    """Common random numbers: the ISR model never touches the response's
    rng, so the two members of a bound pair sit on the SAME pseudo-events
    and their difference is the radiation, not seed noise."""
    off = _response(sampler)
    on = _response(sampler, isr=rad.ISRModel(seed=99))
    assert np.array_equal(off.cell, on.cell)
    # cell, x and q2 are computed BEFORE the ISR block, from the same two
    # uniforms of the same generator, so they are identical -- not close.
    # (Until 2026-08-28 this was an rtol=1e-9 comparison under a comment
    # about allocation history; at that tolerance an ISR path that
    # perturbed x or q2 at the 1e-10 level would have passed.)
    for name in ("x", "q2"):
        assert np.array_equal(getattr(off, name), getattr(on, name)), name
    assert np.any(on.isr_z > 0.0)
    assert np.array_equal(on.y, off.y / (1.0 - on.isr_z))
    # and the model is reproducible on its own
    again = _response(sampler, isr=rad.ISRModel(seed=99))
    assert np.array_equal(again.isr_z, on.isr_z)


def test_isr_moves_the_q2_label_and_leaves_the_mixed_x_alone(sampler):
    """The response's own reconstruction shows the D3 statement: at zero
    hadronic resolution the reconstructed x is the true x whatever the
    radiation, while Q2_e is high by 1/(1 - z)."""
    model = rp.RecoModel(y_had_res=0.0, emcal_stoch=0.0, emcal_const=0.0,
                         beam_e_spread=0.0)
    kw = dict(n_mc_per_cell=120, rng=np.random.default_rng(4))
    on = rp.RecoResponse(sampler, model, isr=rad.ISRModel(seed=5), **kw)
    ok = (on.isr_z > 0.05) & np.isfinite(on.x_reco)
    assert ok.sum() > 50
    # the angular smearing is the only survivor, so compare medians
    assert np.median(on.x_reco[ok] / on.x[ok]) == pytest.approx(1.0, abs=2e-3)
    assert np.median(on.q2_reco[ok] * (1.0 - on.isr_z[ok]) / on.q2[ok]) \
        == pytest.approx(1.0, abs=2e-3)


def test_empz_fraction_recovers_one_minus_z(sampler):
    """(E - p_z)/2E_e = 1 - z, from quantities the chain already
    reconstructs -- the HERA-style handle the bound quotes."""
    model = rp.RecoModel(y_had_res=0.0, emcal_stoch=0.0, emcal_const=0.0,
                         beam_e_spread=0.0)
    on = rp.RecoResponse(sampler, model, n_mc_per_cell=120,
                         rng=np.random.default_rng(4),
                         isr=rad.ISRModel(seed=5))
    f = rad.empz_fraction(on)
    ok = (on.eff > 0) & np.isfinite(f)
    assert np.median(f[ok] / (1.0 - on.isr_z[ok])) == pytest.approx(1.0,
                                                                    abs=3e-3)
    kept = rad.apply_empz_cut(on, 0.85, 1.15)
    assert 0.5 < kept < 1.0
    assert np.all(rad.empz_fraction(on)[on.eff > 0] >= 0.85)


# --- the bound ----------------------------------------------------------------

@pytest.fixture(scope="module")
def bound_inputs(sampler):
    analysis = fom.Scenario(q2_min=1.0, y_min=0.01, y_max=0.95)
    proj = fom.project_rates(CONFIG, analysis, nx=24, nq2=18)
    xe, qe = proj.x_edges, proj.q2_edges
    edges = []
    for i, j in ((8, 6), (11, 8)):
        edges.append((xe[i], xe[i + 2], qe[j], qe[j + 2]))
    cat = bk.tensor_flip_plan(0.6).categories[0]
    return edges, cat


def test_migration_bound_is_stable_at_reduced_statistics(sampler,
                                                         bound_inputs):
    """The bound is a property of the physics, not of the response
    sample: halving the MC statistics moves it by less than the spread
    the two samples themselves carry."""
    edges, cat = bound_inputs
    out = []
    for n_mc in (100, 300):
        def build(isr, n_mc=n_mc):
            return _response(sampler, isr=isr, n_mc=n_mc)
        b = rad.migration_bound(build, edges, cat,
                                isr=rad.ISRModel(seed=20260828))
        out.append(np.array([r["d_amp"] for r in b["rows"]]))
        assert all(0.0 < r["purity_on"] < 1.0 for r in b["rows"])
        assert all(r["purity_on"] <= r["purity_off"] + 0.02
                   for r in b["rows"])
    assert np.all(np.abs(out[0] - out[1]) < 0.02)
    assert np.all(np.abs(out[1]) < 0.05)


def test_empz_cut_removes_the_hard_radiation(sampler, bound_inputs):
    """The mitigation the report quotes next to the uncut bound: the
    HERA-style E - p_z window empties the sample of the hard photons that
    move the Q2 label, at a few per cent of the non-radiative rate."""
    edges, cat = bound_inputs
    on = _response(sampler, isr=rad.ISRModel(seed=20260828), n_mc=300)
    w = on.w * on.eff
    # z > 0.2 is well outside the 0.85 window; 0.10-0.15 survives by
    # construction, which is what sets where the cut can be placed
    before = float(w[on.isr_z > 0.2].sum() / w.sum())
    keep = rad.apply_empz_cut(on, 0.85, 1.15)
    w = on.w * on.eff
    after = float(w[on.isr_z > 0.2].sum() / w.sum())
    assert before > 0.01
    assert after < 0.15 * before
    assert 0.5 < keep < 1.0
    # and the formatter says so
    def build(isr):
        return _response(sampler, isr=isr, n_mc=150)
    txt = rad.format_bound(rad.migration_bound(
        build, edges, cat, isr=rad.ISRModel(seed=20260828),
        empz_cut=(0.85, 1.15)))
    assert "E - p_z cut" in txt and ("PASS" in txt or "FAIL" in txt)


def test_hard_rescalings_are_unity_without_radiation(sampler):
    kern = sampler.kernel
    x = np.array([0.01, 0.05, 0.2])
    q2 = np.array([1.5, 4.0, 20.0])
    y = q2 / (sampler.s * x)
    zero = np.zeros_like(x)
    assert rad.rate_scale(kern, x, q2, y, zero) == pytest.approx(1.0, rel=1e-12)
    assert rad.amplitude_scale(kern, x, q2, y, zero) == pytest.approx(1.0,
                                                                      rel=1e-12)
    # and at small y the amplitude rescaling stays a sub-per-mille effect
    small = rad.amplitude_scale(kern, np.array([0.03]), np.array([3.13]),
                                np.array([3.13 / (sampler.s * 0.03)]),
                                np.array([0.05]))
    assert abs(float(small[0]) - 1.0) < 1e-3


def test_p_emit_renormalises_over_the_allowed_range_by_a_stated_amount():
    """`p_emit` puts the kinematically forbidden mass back into
    (z_min, z_max) instead of dropping it, which over-counts the resolved
    emission -- conservative for a bound, and the docstring's numbers."""
    m = rad.ISRModel()
    q2 = 4.37   # <Q2> of the selected sample of the money script
    for y, quoted in ((0.01, 9e-4), (0.03, 1.3e-3), (0.5, 1.5e-2),
                      (0.9, 6.3e-2)):
        zmax = float(m.z_max(y))
        renorm = float(m.p_emit(q2, y))
        trunc = float(rad.cdf_z(zmax, q2) - rad.cdf_z(m.z_min, q2))
        assert renorm / trunc - 1.0 == pytest.approx(quoted, rel=0.06)


def test_method_bias_table_is_a_strong_function_of_y():
    """The comparison the reports quote must be read at a stated y: the
    electron method's failure at the sweet spots (y = 0.010-0.026) is
    times larger than at the rate-weighted <y> = 0.19 of the whole
    selected sample.  Pinned so the table cannot be re-attributed."""
    z = 0.0923   # <z | z > z_min> of the published-window run
    # the rate-weighted <y> of the whole selected sample, then the four
    # sweet spots of money_cos2phi_reco.py at s = 3980 GeV^2
    y_bar = 0.1889
    y_spots = [1.14 / (3980.0 * 0.0282), 1.14 / (3980.0 * 0.0112),
               3.13 / (3980.0 * 0.0708), 14.3 / (3980.0 * 0.141)]
    rows = {}
    for y in [y_bar] + y_spots:
        rows[y] = dict((r[0], r[1:]) for r in rad.method_bias_table(y, z))
        # the chain's own row does not move at all, at any y
        assert rows[y]["mixed (Q2_e, y_Sigma)"] == pytest.approx(
            (1.0 / (1.0 - z), 1.0, 1.0), rel=1e-12)
        # y + z over y, times (1 - z): the closed form the reports quote
        assert rows[y]["electron"][1] == pytest.approx(
            (y + z) * (1.0 - z) / y, rel=1e-12)
    assert rows[y_bar]["electron"][1] == pytest.approx(1.351, abs=5e-3)
    assert rows[y_bar]["electron"][2] == pytest.approx(0.740, abs=5e-3)
    spot_y = [rows[y]["electron"][1] for y in y_spots]
    assert spot_y == pytest.approx([9.156, 4.184, 8.450, 4.196], abs=5e-3)
    spot_x = [rows[y]["electron"][2] for y in y_spots]
    assert spot_x == pytest.approx([0.109, 0.239, 0.118, 0.238], abs=5e-3)
    # 3 to 7 times the failure the whole-sample <y> shows
    assert 3.0 < min(spot_y) / rows[y_bar]["electron"][1] < 3.2
    assert 6.7 < max(spot_y) / rows[y_bar]["electron"][1] < 6.9


def test_empz_needs_an_independent_hadronic_y(sampler):
    """With y_method = 'electron' the chain's 1 - y IS E'(1-cos)/2E_e, so
    E - p_z would be identically 1 and the cut a silent no-op; the module
    refuses rather than reporting a mitigation that did nothing."""
    on = rp.RecoResponse(sampler, rp.RecoModel(y_method="electron"),
                         n_mc_per_cell=60, rng=np.random.default_rng(3),
                         isr=rad.ISRModel(seed=5))
    with pytest.raises(ValueError, match="y_method"):
        rad.empz_fraction(on)
    with pytest.raises(ValueError, match="y_method"):
        rad.apply_empz_cut(on)


def test_migration_bound_rewinds_the_isr_stream(sampler, bound_inputs):
    """A reused ISRModel gives the same bound: `migration_bound` resets
    the stream, so the answer does not depend on how many bounds the
    caller has already taken with that instance."""
    edges, cat = bound_inputs

    def build(isr):
        return _response(sampler, isr=isr, n_mc=100)

    model = rad.ISRModel(seed=20260828)
    first = rad.migration_bound(build, edges, cat, isr=model)
    second = rad.migration_bound(build, edges, cat, isr=model)
    assert np.array_equal(first["on"].isr_z, second["on"].isr_z)
    assert [r["d_amp"] for r in first["rows"]] == \
           [r["d_amp"] for r in second["rows"]]


def test_migration_bound_seeds_averages_the_response_seeds(sampler,
                                                           bound_inputs):
    """The published bound is a seed average: one response draw carries a
    Monte-Carlo scatter of the size of the bound itself.  Each seed keeps
    its own common-random-number pair, and the mean/sem are those of the
    individual runs."""
    edges, cat = bound_inputs
    seeds = (11, 22, 33)

    def build(seed, isr):
        return _response(sampler, isr=isr, seed=seed, n_mc=100)

    out = rad.migration_bound_seeds(build, edges, cat, seeds,
                                    isr_seed=20260828)
    assert out["seeds"] == [11, 22, 33]
    for i, r in enumerate(out["rows"]):
        vals = np.array([run["rows"][i]["d_amp"] for run in out["runs"]])
        assert r["d_amp"] == pytest.approx(float(vals.mean()), rel=1e-12)
        assert r["d_amp_sem"] == pytest.approx(
            float(vals.std(ddof=1) / np.sqrt(3.0)), rel=1e-12)
        assert np.array_equal(r["d_amp_vals"], vals)
    # every pair still sits on common random numbers
    for run in out["runs"]:
        assert np.array_equal(run["off"].x, run["on"].x)
        assert np.array_equal(run["off"].cell, run["on"].cell)
    txt = rad.format_bound(out)
    assert "mean +- sem over 3 response seeds" in txt and "+-" in txt


def test_response_level_azimuth_residual_is_negligible(sampler):
    """The theorem of plans/08 D3 on the actual response: with the
    physical 6Li mass the covariant azimuth's residual rotation under
    k -> (1 - z)k leaves a fake cos 2phi' far below the amplitudes the
    letter measures (which are 4e-3 to 1e-2).  `isr_dphi` is None with
    the hook off, so the default path carries no extra array."""
    off = _response(sampler)
    assert off.isr_dphi is None
    on = _response(sampler, isr=rad.ISRModel(seed=20260828), n_mc=300)
    w = on.w * on.eff
    fake = 1.0 - float(np.average(np.cos(2.0 * on.isr_dphi), weights=w))
    assert np.abs(on.isr_dphi).max() < 0.2
    assert abs(fake) < 1e-5
