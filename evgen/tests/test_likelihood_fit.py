"""The acceptance-profiled likelihood estimator of the two-azimuth
harmonics (plans/08 A12) and the in-situ (u1, u2) fit (plans/08 A3).

`reco.harmonic_ratio_fit_2d` inverts the bin-wise spin-state ratio, and
that inversion is strongly curved over the range the ratio explores at low
counts -- concave for the flip plan, whose Pbar = -0.3 sets the sign of
its second derivative -- so it carries a negative per-bin Jensen offset of
order 1/nu_b that the weighted LSQ projects onto the constant and onto
a_t.
`reco.harmonic_likelihood_fit_2d` profiles the per-bin acceptance out of
the Poisson likelihood instead; the profile is EXACTLY the conditional
multinomial given the bin totals, whose score has zero mean bin by bin at
any count.  These tests pin that difference against the chain's own
sparsest published |t| bin, and pin that nothing else about the estimator
changed -- the same basis, the same assumptions, the same blindness to a
fill-dependent acceptance, the same rank diagnostics, and the same
default.
"""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen import bookkeeping as bk  # noqa: E402
from polligen import coherent as coh  # noqa: E402
from polligen import reco, recopseudo as rp  # noqa: E402

from polli_fastsim import beams  # noqa: E402

PLAN = bk.tensor_flip_plan(0.6)
PZZ = [float(c.moments()[1]) for c in PLAN.categories]
FRAC = [c.lumi_fraction for c in PLAN.categories]
U = (0.05, 0.02)
N_ALPHA, N_BETA = 12, 24


# --- the published working point, in the sparsest |t| bin ------------------

@pytest.fixture(scope="module")
def sparse_bin():
    """The coherent chain at the tagging optics of the low configuration,
    in the |t| bin that carries the fewest counts of the four published
    ones (0.17-0.25 GeV^2), with the produced flux scaled so that the
    expected occupancy is the 48 counts per (alpha, beta) cell of Report 2
    Table 5.  Returns (mu, basis, truth, fit_kwargs) with mu the EXACT
    expected counts -- every study below draws Poisson from this one mu,
    so the response Monte Carlo is common to both estimators and cancels
    from their comparison."""
    cfg = beams.default_configs("6Li")[0]
    sc = coh.CoherentScenario(amp=0.01, eps_b0=-0.08)
    top = reco.tagging_optics_point(cfg, slope_b=sc.slope_b)
    cr = rp.CoherentResponse(sc, cfg, (top["sigma_x_eff"], top["sigma_y"]),
                             n_mc=400000, rng=np.random.default_rng(20260828),
                             cut_theta_xy=(top["env_x"], top["env_y"]))
    tlo, thi = 0.17, 0.25
    n_produced = 48.0 * N_ALPHA * N_BETA / cr.t_bin_fraction(tlo, thi)
    fit = rp.measure_coherent(
        cr, n_produced, PLAN, tlo, thi, a_e=0.010,
        a_t_func=lambda t: sc.cos2phi_coefficient_deformation(t, 1.0),
        u1=U[0], u2=U[1], kappa=0.002, with_sin=True, poisson=False)
    kw = dict(u_coeffs=U, beta_means=fit["beta_means"], with_sin=True)
    return fit["expected"], (fit["alpha_edges"], fit["beta_edges"]), \
        fit["truth"], kw


def _both(counts, edges, kw):
    ae, be = edges
    return (reco.harmonic_ratio_fit_2d(counts, FRAC, PZZ, ae, be, **kw),
            reco.harmonic_likelihood_fit_2d(counts, FRAC, PZZ, ae, be, **kw))


# --- 1. exactness and high-count agreement --------------------------------

def test_likelihood_is_exact_on_asimov_counts(sparse_bin):
    """On exact expected counts the conditional score vanishes identically
    at the truth (sum_f mu_{f,b} s_{f,b} = eps_b [Pbar - Pbar D/D] = 0), so
    the fit must return the injected coefficients to machine precision --
    where the ratio fit's own closure is limited by the nonlinearity of
    its inversion (1e-6 here)."""
    mu, edges, truth, kw = sparse_bin
    ratio, lik = _both(mu, edges, kw)
    assert lik["a_t"] == pytest.approx(truth["a_t"], rel=1e-9)
    assert lik["a_e"] == pytest.approx(0.010, rel=1e-9)
    assert lik["const"] == pytest.approx(0.002, rel=1e-9)
    for key in ("a_e_s", "a_t_s", "a_m_s"):
        assert abs(lik[key]) < 1e-12
    assert lik["grad_max"] < 1e-6 and lik["n_iter"] < 20
    # the ratio fit closes far less well on the same input ...
    assert abs(ratio["a_t"] / truth["a_t"] - 1.0) > 1e-8
    # ... and yet the two quote essentially the same Asimov errors, so no
    # published error bar moves if the estimator is swapped.  The residual
    # is the difference between the delta-method variance of the inverted
    # ratio and the exact Fisher information: 0.24% on err_t and 0.009% on
    # err_e here, and INDEPENDENT of the luminosity (the same figure at 48
    # and at 2 x 10^4 counts per bin), with the likelihood the smaller of
    # the two, as an efficient estimator must be.
    for key in ("err_t", "err_e", "err_m"):
        assert lik[key] == pytest.approx(ratio[key], rel=3e-3)
        assert lik[key] <= ratio[key]


def test_likelihood_and_ratio_agree_at_high_counts(sparse_bin):
    """The two estimators converge as the counts grow, and the way they
    converge identifies what separates them: the discrepancy on a_t is the
    ratio's per-bin offset, which falls as 1/nu_b, while a_e -- flat in
    alpha, where the acceptance is flat too -- never differs at all.  At
    2 x 10^4 counts per (alpha, beta) bin the two agree to a fifth of the
    statistical error draw by draw."""
    mu, edges, _truth, kw = sparse_bin
    means = {}
    for target in (1500.0, 20000.0):
        rng = np.random.default_rng(3)
        d_t, d_e = [], []
        for _ in range(10):
            counts = rng.poisson(mu * (target / 48.0))
            ratio, lik = _both(counts, edges, kw)
            d_t.append(lik["a_t"] - ratio["a_t"])
            d_e.append((lik["a_e"] - ratio["a_e"]) / lik["err_e"])
        means[target] = float(np.mean(d_t))
        assert np.max(np.abs(d_e)) < 0.1
        tol = 0.7 if target < 2e3 else 0.2
        assert np.max(np.abs(d_t)) < tol * lik["err_t"]
    # the discrepancy falls like 1/counts -- it is the ratio's per-bin
    # offset, not noise (measured 16.5 against the 13.3 the counts predict)
    assert means[1500.0] / means[20000.0] == pytest.approx(20000.0 / 1500.0,
                                                           rel=0.35)


# --- 2. unbiasedness where the ratio fails --------------------------------

def test_likelihood_is_unbiased_where_the_ratio_fails(sparse_bin):
    """The defect and its fix, on the same 300 Poisson draws of the same
    expected counts: at the 48 counts per (alpha, beta) cell of the
    sparsest published |t| bin the ratio's a_t is low by more than 15%
    (measured -37%, reproducing Report 2 Table 5's -37% at 12 x 24), while
    the likelihood's mean sits within 4 standard errors of the injected
    value."""
    mu, edges, truth, kw = sparse_bin
    assert mu.sum() / (N_ALPHA * N_BETA) == pytest.approx(48.0, rel=0.05)
    rng = np.random.default_rng(101)
    ratio, lik = [], []
    for _ in range(300):
        counts = rng.poisson(mu)
        r, l = _both(counts, edges, kw)
        ratio.append(r["a_t"])
        lik.append(l["a_t"])
    ratio, lik = np.asarray(ratio), np.asarray(lik)
    assert ratio.mean() / truth["a_t"] - 1.0 < -0.15
    assert abs(lik.mean() - truth["a_t"]) < 4 * lik.std() / np.sqrt(lik.size)
    # and the fix costs nothing in variance -- it is the smaller of the two
    assert lik.std() <= 1.05 * ratio.std()


# --- 3. the quoted errors are the ensemble spread -------------------------

def test_likelihood_errors_reproduce_the_ensemble_spread(sparse_bin):
    """The compressed ratio estimator fluctuates LESS than it says it does
    at low counts; the likelihood's quoted error is the spread at 48 counts
    per bin and at one count per bin alike."""
    mu, edges, _truth, kw = sparse_bin
    ae, be = edges
    for scale, ratio_is_wrong in ((1.0, False), (1.0 / 48.0, True)):
        rng = np.random.default_rng(8)
        lik, ratio = [], []
        for _ in range(250):
            counts = rng.poisson(mu * scale)
            f = reco.harmonic_likelihood_fit_2d(counts, FRAC, PZZ, ae, be,
                                                **kw)
            lik.append((f["a_t"], f["err_t"]))
            g = reco.harmonic_ratio_fit_2d(counts, FRAC, PZZ, ae, be, **kw)
            ratio.append((g["a_t"], g["err_t"]))
        lik, ratio = np.asarray(lik), np.asarray(ratio)
        assert lik[:, 0].std() == pytest.approx(lik[:, 1].mean(), rel=0.15)
        if ratio_is_wrong:
            # at ~1 count per bin the ratio quotes an error 60% larger than
            # the spread it actually has
            assert ratio[:, 1].mean() > 1.3 * ratio[:, 0].std()


# --- 4. the parity-forbidden null columns ---------------------------------

@pytest.fixture(scope="module")
def slot():
    """The tight-cutout fixture the sin-null tests of test_recopseudo.py
    use: half the beta bins live, so the null test is genuinely loaded."""
    cfg = beams.default_configs("6Li")[1]
    sc = coh.CoherentScenario(amp=0.01, eps_b0=-0.08)
    sig = 0.218 / (10.0 * cfg.ion.A * cfg.ion_momentum_per_nucleon)
    return sc, rp.CoherentResponse(sc, cfg, sig, aspect=1.0,
                                   cut_scale_xy=(2.5, 1.0), n_mc=300000,
                                   rng=np.random.default_rng(11))


def _rotated(cr, sc, delta, roll_only=False, drop_s2t=False, n=8e7):
    base = (lambda t: sc.cos2phi_coefficient_deformation(t, 1.0))
    c, s_ = np.cos(2 * delta), np.sin(2 * delta)
    kw = dict(u1=U[0], u2=U[1], poisson=False, with_sin=True,
              a_t_s_func=lambda t: base(t) * s_, fit="likelihood")
    if not roll_only:
        kw["a_e_s"] = 0.010 * s_
    if drop_s2t:
        ae = np.linspace(0.0, 2.0 * np.pi, N_ALPHA + 1)
        be = np.linspace(0.0, 2.0 * np.pi, N_BETA + 1)
        mu = cr.expected_counts_2d(n, PZZ, FRAC, 0.05, 0.08, ae, be,
                                   0.010 if roll_only else 0.010 * c,
                                   lambda t: base(t) * c, u1=U[0], u2=U[1],
                                   a_t_s_func=lambda t: base(t) * s_)
        bm = dict(cr.basis_means(0.05, 0.08, be))
        bm.pop("s2t")
        return reco.harmonic_likelihood_fit_2d(mu, FRAC, PZZ, ae, be,
                                               u_coeffs=U, beta_means=bm,
                                               with_sin=True)
    return rp.measure_coherent(cr, n, PLAN, 0.05, 0.08,
                               0.010 if roll_only else 0.010 * c,
                               lambda t: base(t) * c, **kw)


def test_sin_nulls_survive_the_likelihood(slot):
    """All three parity-forbidden partners vanish, turning the columns on
    leaves the cos sector alone, an axis tilt gives tan 2delta in BOTH
    azimuths, a pot roll only in beta, and the beta null still needs the
    t-weighted sine template."""
    sc, cr = slot
    base = (lambda t: sc.cos2phi_coefficient_deformation(t, 1.0))
    common = dict(u1=U[0], u2=U[1], poisson=False, fit="likelihood")
    off = rp.measure_coherent(cr, 8e7, PLAN, 0.05, 0.08, 0.010, base,
                              **common)
    on = rp.measure_coherent(cr, 8e7, PLAN, 0.05, 0.08, 0.010, base,
                             with_sin=True, **common)
    for key in ("a_e_s", "a_t_s", "a_m_s"):
        assert abs(on[key]) < 1e-10
    assert on["a_e"] == pytest.approx(off["a_e"], abs=1e-12)
    assert on["a_t"] == pytest.approx(off["a_t"], abs=1e-9)

    for delta in (0.02, 0.05):
        f = _rotated(cr, sc, delta)
        assert f["a_e_s"] / f["a_e"] == pytest.approx(np.tan(2 * delta),
                                                      rel=1e-3)
        assert f["a_t_s"] / f["a_t"] == pytest.approx(np.tan(2 * delta),
                                                      rel=1e-3)
        g = _rotated(cr, sc, delta, roll_only=True)
        assert abs(g["a_e_s"] / g["a_e"]) < 1e-5
        assert g["a_t_s"] / g["a_t"] == pytest.approx(np.tan(2 * delta),
                                                      rel=1e-3)
    bad = _rotated(cr, sc, 0.05, roll_only=True, drop_s2t=True)
    assert abs(bad["a_t_s"] / bad["a_t"] / np.tan(0.10) - 1.0) > 0.03


def test_low_count_sin_nulls_are_consistent_with_zero(sparse_bin):
    """The null test must still be a null where the ratio's cos sector is
    already 37% low: the sin means over 200 draws sit within 4 standard
    errors of zero."""
    mu, edges, _truth, kw = sparse_bin
    ae, be = edges
    rng = np.random.default_rng(21)
    vals = []
    for _ in range(200):
        f = reco.harmonic_likelihood_fit_2d(rng.poisson(mu), FRAC, PZZ, ae,
                                            be, **kw)
        vals.append([f["a_e_s"], f["a_t_s"], f["a_m_s"]])
    vals = np.asarray(vals)
    se = vals.std(axis=0) / np.sqrt(vals.shape[0])
    assert np.all(np.abs(vals.mean(axis=0)) < 4 * se)


# --- 5. single-fill and empty bins ----------------------------------------

def test_single_fill_and_empty_bins_need_no_variance_floor():
    """The R3 pathology, removed structurally: a bin populated by ONE fill
    contributes a finite n ln p term and an empty bin contributes nothing
    at all, so no bin can acquire an infinite weight.  The counts are the
    R3 counter-example spread over eight phi' bins (two empty, two
    single-fill)."""
    counts = np.array([[100.0, 50.0, 0.0, 3.0, 80.0, 40.0, 0.0, 7.0],
                       [90.0, 60.0, 5.0, 0.0, 70.0, 55.0, 6.0, 0.0]])
    edges = np.linspace(0.0, 2.0 * np.pi, 9)
    fit = reco.harmonic_likelihood_fit(counts, [0.5, 0.5], [0.6, -1.2],
                                       edges)
    assert np.isfinite(fit["amp"]) and np.isfinite(fit["err"])
    assert fit["err"] > 0.0


def test_sparse_two_dimensional_counts_give_a_positive_definite_covariance(
        sparse_bin):
    """One count per (alpha, beta) cell: most bins empty, many populated by
    one fill.  Every draw must return finite parameters and a
    positive-definite covariance."""
    mu, edges, _truth, kw = sparse_bin
    ae, be = edges
    rng = np.random.default_rng(13)
    for _ in range(20):
        counts = rng.poisson(mu / 48.0)
        assert (counts.sum(axis=0) == 0).sum() > 50      # genuinely sparse
        assert (counts.min(axis=0) == 0).sum() > 100     # single-fill bins
        f = reco.harmonic_likelihood_fit_2d(counts, FRAC, PZZ, ae, be, **kw)
        assert np.all(np.isfinite(f["cov"]))
        assert np.all(np.linalg.eigvalsh(f["cov"]) > 0.0)
        assert np.all(np.isfinite([f["a_e"], f["a_t"], f["a_m"]]))


# --- 6. the systematics it is equally blind to ----------------------------

def test_blind_systematics_are_identical_to_the_ratio(slot):
    """The likelihood makes the SAME common-acceptance assumption as the
    ratio, so it must be blind in the same way and by the same amount: a
    per-fill cutout, a wrong u2 and a wrong luminosity share move both
    estimators identically on exact counts.  It must not be advertised as
    fixing a systematic it shares.  The leading shift of each -- a_t for
    the cutout and the luminosity, a_e for u2, which is what Report 2
    Table 6 quotes -- agrees to better than 1% relative (measured 0.4% on
    a_t, 0.7% on a_e, the residual being second order in a deliberately
    LARGE 12% shift), and no parameter's shift differs by as much as 2% of
    its own statistical error."""
    sc, cr = slot
    base = (lambda t: sc.cos2phi_coefficient_deformation(t, 1.0))
    common = dict(u1=U[0], u2=U[1], poisson=False)

    def pair(**kw):
        a = rp.measure_coherent(cr, 8e7, PLAN, 0.05, 0.08, 0.010, base,
                                fit="ratio", **common, **kw)
        b = rp.measure_coherent(cr, 8e7, PLAN, 0.05, 0.08, 0.010, base,
                                fit="likelihood", **common, **kw)
        return a, b

    ref_r, ref_l = pair()
    # the slot's BINDING half-width is the vertical one; perturbing the
    # horizontal 25-sigma edge moves no recoil at all
    cases = ((dict(responses=[cr.with_cut(eff_scale_xy=(1.0, 1.001)), cr]),
              "a_t"),
             (dict(u_coeffs_assumed=(0.05, 0.044)), "a_e"),
             (dict(lumi_assumed=[0.5005, 0.4995]), "a_t"))
    for kw, leading in cases:
        got_r, got_l = pair(**kw)
        assert abs(got_r[leading] - ref_r[leading]) > 1e-6   # it does move
        for key, err in (("a_t", "err_t"), ("a_e", "err_e"),
                         ("a_m", "err_m"), ("const", "err_t")):
            d_r = got_r[key] - ref_r[key]
            d_l = got_l[key] - ref_l[key]
            assert abs(d_l - d_r) < 0.02 * ref_l[err]
            if key == leading:
                assert d_l == pytest.approx(d_r, rel=1e-2)


# --- 7. the default has not moved -----------------------------------------

def test_measure_coherent_default_is_the_ratio_fit(sparse_bin, slot):
    """Every published coherent number is the ratio fit, so the default
    must reproduce `harmonic_ratio_fit_2d` BIT FOR BIT."""
    sc, cr = slot
    base = (lambda t: sc.cos2phi_coefficient_deformation(t, 1.0))
    kw = dict(u1=U[0], u2=U[1], with_sin=True)
    a = rp.measure_coherent(cr, 8e7, PLAN, 0.05, 0.08, 0.010, base,
                            poisson=False, **kw)
    b = rp.measure_coherent(cr, 8e7, PLAN, 0.05, 0.08, 0.010, base,
                            poisson=False, fit="ratio", **kw)
    direct = reco.harmonic_ratio_fit_2d(a["expected"], FRAC, PZZ,
                                        a["alpha_edges"], a["beta_edges"],
                                        u_coeffs=U,
                                        beta_means=a["beta_means"],
                                        with_sin=True)
    for key in ("a_e", "a_t", "a_m", "err_e", "err_t", "err_m", "const"):
        np.testing.assert_allclose(a[key], b[key], rtol=0.0, atol=0.0)
        np.testing.assert_allclose(a[key], direct[key], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(a["cov"], direct["cov"], rtol=0.0, atol=0.0)
    with pytest.raises(ValueError, match="ratio.*likelihood"):
        rp.measure_coherent(cr, 8e7, PLAN, 0.05, 0.08, 0.010, base,
                            poisson=False, fit="lsq", **kw)


# --- 8. the rank guard ----------------------------------------------------

def test_rank_guard_names_its_cause_from_both_estimators():
    """A design whose live bins no longer span the template space is a
    statement about the ACCEPTANCE, and both estimators must say so: the
    Fisher information of the profiled likelihood is the same Gram matrix
    as the weighted least-squares design, up to weights, so they fail
    together.  Here only two alpha bins carry counts, half a period apart,
    where <cos 2a> and <sin 2a> are each constant and therefore collinear
    with the constant column."""
    counts = np.zeros((2, N_ALPHA, N_BETA))
    counts[:, [0, N_ALPHA // 2], :] = 100.0
    ae = np.linspace(0.0, 2.0 * np.pi, N_ALPHA + 1)
    be = np.linspace(0.0, 2.0 * np.pi, N_BETA + 1)
    for fit in (reco.harmonic_ratio_fit_2d, reco.harmonic_likelihood_fit_2d):
        with pytest.raises(np.linalg.LinAlgError) as exc:
            fit(counts, FRAC, PZZ, ae, be, with_sin=True)
        msg = str(exc.value)
        assert "rank-deficient" in msg
        assert "%d of %d" % (2 * N_BETA, N_ALPHA * N_BETA) in msg
        assert "e/t/m/e_s/t_s/m_s" in msg
        assert "7-column" in msg


def test_a_stalled_newton_is_named_rather_than_returned(sparse_bin):
    """A fill with no counts anywhere leaves the profiled likelihood with
    no interior maximum: the positivity guard halves the Newton step
    towards the boundary 1 + u_b + P_f T_b = 0 until the step underflows,
    which the step test alone reads as convergence.  The gradient AT the
    returned point is what separates the two, and it must be checked --
    otherwise the fit returns a boundary point with an error bar a
    thousandth of the truth and nothing to say so (2026-08-28 review,
    finding 5).
    The guard costs nothing where the estimator is used: at the published
    48 counts per (alpha, beta) cell, and still at one count per cell, the
    gradient at the solution is 1e-15 error bars against a tolerance of
    1e-3."""
    counts = np.zeros((2, N_ALPHA, N_BETA))
    counts[0] = 30.0                       # the -1.2 fill collected nothing
    ae = np.linspace(0.0, 2.0 * np.pi, N_ALPHA + 1)
    be = np.linspace(0.0, 2.0 * np.pi, N_BETA + 1)
    with pytest.raises(np.linalg.LinAlgError) as exc:
        reco.harmonic_likelihood_fit_2d(counts, FRAC, PZZ, ae, be)
    assert "positivity boundary" in str(exc.value)
    # what it would have returned: the boundary of the EMPTY fill,
    # 1 + P_- T = 0 with P_- = -1.2, i.e. T = 5/6, and an error bar three
    # orders below the a_t error this occupancy really carries
    theta, _it, grad, _q, _d = reco._profile_likelihood_newton(
        counts.reshape(2, N_ALPHA * N_BETA), np.array(FRAC) / sum(FRAC),
        np.array(PZZ), np.vstack([np.ones(N_ALPHA * N_BETA)]
                                 + [reco.basis_2d(ae, be)[c]
                                    for c in ("e", "t", "m")]).T,
        np.zeros(N_ALPHA * N_BETA))
    assert theta[0] == pytest.approx(5.0 / 6.0, rel=1e-5)
    assert np.max(np.abs(grad)) > 1.0
    # and the good case passes it: the published sparse bin, drawn Poisson
    mu, edges, _truth, kw = sparse_bin
    rng = np.random.default_rng(4)
    for scale in (1.0, 1.0 / 48.0):
        for _ in range(20):
            f = reco.harmonic_likelihood_fit_2d(
                rng.poisson(mu * scale), FRAC, PZZ, *edges, **kw)
            assert f["grad_max"] * f["err_t"] < 1e-9


def test_four_phi_bins_are_named_rather_than_silently_infinite():
    """Four bins of width pi/2 have <cos 2phi> = 0 identically, so the
    amplitude is not measurable at all.  The ratio fit returns a finite
    but meaningless 10^14 error there; the likelihood names the cause."""
    counts = np.array([[100.0, 50.0, 0.0, 3.0], [90.0, 60.0, 5.0, 0.0]])
    edges = np.linspace(0.0, 2.0 * np.pi, 5)
    loose = reco.harmonic_ratio_fit(counts, [0.5, 0.5], [0.6, -1.2], edges)
    assert loose["err"] > 1e10
    with pytest.raises(np.linalg.LinAlgError, match="rank-deficient"):
        reco.harmonic_likelihood_fit(counts, [0.5, 0.5], [0.6, -1.2], edges)


# --- the one-azimuth twin -------------------------------------------------

def test_one_dimensional_likelihood_matches_the_ratio(slot):
    """The inclusive phi' bins carry >= 10^4 counts, where the two
    estimators must agree; on exact counts both return the injected
    amplitude with the finite-bin dilution divided out."""
    edges = np.linspace(0.0, 2.0 * np.pi, 25)
    centers = 0.5 * (edges[:-1] + edges[1:])
    dil = np.sin(2.0 * (np.pi / 24)) / (2.0 * np.pi / 24)
    amp, kappa = 0.007, 0.001
    eps = 1.0 + 0.03 * np.cos(2 * centers) + 0.02 * np.cos(centers)
    t = kappa + amp * dil * np.cos(2 * centers)
    mu = np.array([lf * 2e5 * eps * (1.0 + pf * t)
                   for pf, lf in zip(PZZ, FRAC)])
    exact = reco.harmonic_likelihood_fit(mu, FRAC, PZZ, edges)
    assert exact["amp"] == pytest.approx(amp, rel=1e-9)
    assert exact["const"] == pytest.approx(kappa, rel=1e-9)
    rng = np.random.default_rng(5)
    counts = rng.poisson(mu)
    lik = reco.harmonic_likelihood_fit(counts, FRAC, PZZ, edges,
                                       with_sin=True)
    ratio = reco.harmonic_ratio_fit(counts, FRAC, PZZ, edges, with_sin=True)
    assert abs(lik["amp"] - ratio["amp"]) < 0.05 * lik["err"]
    assert lik["err"] == pytest.approx(ratio["err"], rel=5e-3)


# --- the in-situ (u1, u2) fit (plans/08 A3) -------------------------------

def test_insitu_u_recovers_the_injected_value(slot):
    """With a free per-bin acceptance u is not identifiable at all, so the
    in-situ fit uses the response's own acceptance shape (what an analysis
    takes from its acceptance MC).  Given that shape it must return the
    generated (u1, u2) exactly on exact counts."""
    sc, cr = slot
    base = (lambda t: sc.cos2phi_coefficient_deformation(t, 1.0))
    f = rp.measure_coherent(cr, 8e7, PLAN, 0.05, 0.08, 0.010, base,
                            u1=U[0], u2=U[1], poisson=False,
                            u_coeffs_assumed="in-situ", fit="likelihood")
    assert f["u_insitu"][0] == pytest.approx(U[0], rel=1e-8)
    assert f["u_insitu"][1] == pytest.approx(U[1], rel=1e-8)
    # and the harmonics still close, now with u measured rather than known
    assert f["a_t"] == pytest.approx(f["truth"]["a_t"], rel=1e-8)
    assert f["a_e"] == pytest.approx(0.010, rel=1e-8)
    # one alternation of the two fits is NOT enough: the pair converges by
    # a factor ~40 per round, and the default of five reaches the counts
    one = rp.measure_coherent(cr, 8e7, PLAN, 0.05, 0.08, 0.010, base,
                              u1=U[0], u2=U[1], poisson=False,
                              u_coeffs_assumed="in-situ", fit="likelihood",
                              u_iter=1)
    assert abs(one["u_insitu"][0] / U[0] - 1.0) > 1e-3
    assert f["u_err"][0] > 0.0 and f["u_err"][1] > 0.0
    with pytest.raises(ValueError, match="in-situ"):
        rp.measure_coherent(cr, 8e7, PLAN, 0.05, 0.08, 0.010, base,
                            poisson=False, u_coeffs_assumed="fitted")


def test_insitu_u_propagates_into_the_a_e_error(slot):
    """The propagated error is the Jacobian of the u2 systematic Report 2
    Table 6 quotes: da_e/du2 times a ZEUS 1-sigma du2 = 0.024 must
    reproduce the shift measured by simply assuming a wrong u2, and the
    propagated term must enter the returned a_e error in quadrature."""
    sc, cr = slot
    base = (lambda t: sc.cos2phi_coefficient_deformation(t, 1.0))
    common = dict(u1=U[0], u2=U[1], poisson=False, fit="likelihood")
    ref = rp.measure_coherent(cr, 8e7, PLAN, 0.05, 0.08, 0.010, base,
                              **common)
    wrong = rp.measure_coherent(cr, 8e7, PLAN, 0.05, 0.08, 0.010, base,
                                u_coeffs_assumed=(U[0], U[1] + 0.024),
                                **common)
    insitu = rp.measure_coherent(cr, 8e7, PLAN, 0.05, 0.08, 0.010, base,
                                 u_coeffs_assumed="in-situ", **common)
    shift = wrong["a_e"] - ref["a_e"]
    assert shift == pytest.approx(0.024 * insitu["u_jacobian"][1, 1],
                                  rel=2e-3)
    stat = np.sqrt(insitu["cov_stat"][1, 1])
    extra = np.sqrt(insitu["cov_u"][1, 1])
    # the propagated term, rebuilt from OUTSIDE measure_coherent: the
    # Jacobian by refitting at u +- h and the u covariance the in-situ fit
    # quotes.  (That err_e is the quadrature sum of the two is an identity
    # of the implementation and would hold for any Jacobian at all, so it
    # is the Jacobian and the covariance that have to be checked.)
    jac = np.zeros(2)
    for j, u_hat in enumerate(insitu["u_insitu"]):
        h = 0.25 * insitu["u_err"][j]
        step = [0.0, 0.0]
        step[j] = h
        hi, lo = (rp.measure_coherent(
            cr, 8e7, PLAN, 0.05, 0.08, 0.010, base, u_coeffs_assumed=(
                insitu["u_insitu"][0] + sgn * step[0],
                insitu["u_insitu"][1] + sgn * step[1]), **common)["a_e"]
            for sgn in (+1.0, -1.0))
        jac[j] = (hi - lo) / (2.0 * h)
    assert np.sqrt(jac @ insitu["u_cov"] @ jac) == pytest.approx(extra,
                                                                 rel=1e-3)
    assert insitu["err_e"] >= stat
    # the in-situ u is far better than the ZEUS band at this luminosity:
    # the leakage it leaves is a small fraction of the assumed-u shift
    assert extra < 0.2 * abs(shift)
