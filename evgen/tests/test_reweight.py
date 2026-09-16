"""Gate: the Mode-W reweighter of plans/05 step 5.C.

Everything here runs on a SYNTHETIC in-memory pool -- no network, no
container -- except `test_streamed_beagle_pool_closure`, which reads a
cached `tools/analysis/dump_spectators.py` CSV named by the environment
variable POLLIGEN_BEAGLE_CSV and skips when it is not there.

The acceptance criterion of the task (plans/05:276) is the last two
tests: over >= 200 pseudo-experiments the pull mean is |mu| < 0.15 and the
pull width 1.00 +- 0.10 on both A_par and A_zz.
"""

import os
import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen import bookkeeping as bk  # noqa: E402
from polligen import estimators as est  # noqa: E402
from polligen import reweight as rw  # noqa: E402
from polligen.xsec import InclusiveKernel  # noqa: E402

from polli_fastsim import beams, fom  # noqa: E402
from polli_fastsim.asymmetries import (  # noqa: E402
    err_a_parallel, err_azz)
from polli_fastsim.polarized import toy_b1, toy_delta_gluon  # noqa: E402

CONFIG = beams.default_configs("6Li")[1]
S = CONFIG.sqrt_s_per_nucleon ** 2
NPOOL = 60_000
NTRIALS = 1000        # the gate asks >= 200; 1000 puts the |mu| < 0.15
                      # threshold at 4.7 sigma of the pull mean
NTRIALS_DRAW = 500    # the resampling variant, 4 times the cost per trial
PE, PZ, PZZ = 0.7, 0.6, 0.6

# the acceptance gate of plans/05 step 5.C
PULL_MEAN_MAX = 0.15
PULL_WIDTH_TOL = 0.10


def _pool(kern, seed, n=NPOOL):
    return rw.unpolarized_pool(kern, S, n, x_range=(0.03, 0.3),
                               q2_range=(4.0, 60.0),
                               rng=np.random.default_rng(seed),
                               scenario=fom.Scenario(),
                               electron_energy=CONFIG.electron_energy)


@pytest.fixture(scope="module")
def vector_rw():
    kern = InclusiveKernel(beams.LI6)          # vector sector only
    return rw.ModeWReweighter(kern, _pool(kern, 7))


@pytest.fixture(scope="module")
def tensor_rw():
    kern = InclusiveKernel(beams.LI6, b1_func=toy_b1)
    return rw.ModeWReweighter(kern, _pool(kern, 8))


# --- the sample container ------------------------------------------------


def test_external_sample_guards():
    x = np.array([0.1, 0.2])
    with pytest.raises(ValueError):
        rw.ExternalSample(x, x[:1], x, 1000.0)
    with pytest.raises(ValueError):
        rw.ExternalSample(x, x, x, -1.0)
    with pytest.raises(ValueError):
        rw.ExternalSample(x, x, x, 1000.0, weight=[1.0, -1.0])
    with pytest.raises(ValueError):
        rw.ExternalSample(x.reshape(2, 1), x.reshape(2, 1),
                          x.reshape(2, 1), 1000.0)
    s = rw.ExternalSample(x, np.array([4.0, 9.0]), np.array([0.0, 7.0]),
                          1000.0, extra={"tag": np.array([1, 2])})
    assert len(s) == 2
    assert s.phi[1] == pytest.approx(7.0 - 2.0 * np.pi)   # wrapped
    assert s.y == pytest.approx([4.0 / 100.0, 9.0 / 200.0])
    half = s.select([True, False])
    assert len(half) == 1 and half.extra["tag"][0] == 1


def test_unpolarized_pool_is_unpolarized():
    """The one property the closure rests on: phi flat, so the pool
    carries no cos phi or cos 2phi of its own."""
    kern = InclusiveKernel(beams.LI6)
    pool = _pool(kern, 3)
    n = len(pool)
    assert n > NPOOL // 4
    for k in (1, 2):
        moment = float(np.mean(np.cos(k * pool.phi)))
        assert abs(moment) < 5.0 / np.sqrt(2.0 * n)
    # and it follows the unpolarized cross section: the accepted pool's
    # x spectrum is far steeper than the log-uniform proposal
    assert np.median(pool.x) < 0.09


# --- BeAGLE CSV -> kinematics (no network) -------------------------------


def _fake_beagle_rows(x, q2, phi, e_beam=9.0, p_ion=130.0, m_n=0.9383):
    """A dump_spectators.py-shaped record array built BACKWARDS from known
    (x, Q2, phi) -- the exact inverse of `sample_from_beagle_rows`.

    Head-on beams, electron along -z.  With A = Q2/(2E) and
    B = (1 - y) E (E_i + p_i), the two invariants Q2 = 2 E E'(1 + cos th)
    and y = 1 - E'(E_i - p_i cos th)/(E (E_i + p_i)) invert exactly:
        cos th = (A E_i - B)/(A p_i + B),   E' = A/(1 + cos th).
    """
    dt = np.dtype([("ievt", int), ("kind", "U1"), ("pdg", int),
                   ("px", float), ("py", float), ("pz", float),
                   ("e", float)])
    e_ion = np.sqrt(p_ion ** 2 + m_n ** 2)
    pk = e_beam * (e_ion + p_ion)          # P.k, exact
    rows = []
    for i, (xi, qi, ph) in enumerate(zip(x, q2, phi)):
        y = qi / (2.0 * xi * pk)
        a = qi / (2.0 * e_beam)
        b = (1.0 - y) * e_beam * (e_ion + p_ion)
        cos_th = (a * e_ion - b) / (a * p_ion + b)
        e_p = a / (1.0 + cos_th)
        sin_th = np.sqrt(max(1.0 - cos_th ** 2, 0.0))
        rows.append((i, "B", 2112, 0.0, 0.0, p_ion, e_ion))
        rows.append((i, "E", 11, e_p * sin_th * np.cos(ph),
                     e_p * sin_th * np.sin(ph), e_p * cos_th, e_p))
    return np.array(rows, dtype=dt)


def test_beagle_row_reconstruction_round_trip():
    """(x, Q2, phi) -> four-vectors -> back, through the same invariants
    the streamed loader forms."""
    x = np.array([0.01, 0.05, 0.2])
    q2 = np.array([5.0, 20.0, 60.0])
    phi = np.array([0.3, 2.9, 5.5])
    rows = _fake_beagle_rows(x, q2, phi)
    ievt, xr, q2r, yr, phir = rw.sample_from_beagle_rows(
        rows, rows["kind"], electron_beam=(0.0, 0.0, -9.0, 9.0))
    assert ievt.tolist() == [0, 1, 2]
    assert q2r == pytest.approx(q2, rel=1e-9)
    assert xr == pytest.approx(x, rel=1e-9)
    assert phir == pytest.approx(phi, abs=1e-12)
    assert np.all((yr > 0.0) & (yr < 1.0))


def test_beagle_csv_round_trip(tmp_path):
    """The loader end to end: the same rows written as a CSV in the
    dump_spectators.py column order and read back."""
    x = np.array([0.01, 0.05])
    q2 = np.array([5.0, 20.0])
    phi = np.array([0.3, 2.9])
    rows = _fake_beagle_rows(x, q2, phi)
    csv = tmp_path / "ed.csv"
    with open(csv, "w") as fh:
        fh.write("ievt,kind,pdg,px,py,pz,e\n")
        for r in rows:
            fh.write("%d,%s,%d,%.17g,%.17g,%.17g,%.17g\n"
                     % (r["ievt"], r["kind"], r["pdg"], r["px"], r["py"],
                        r["pz"], r["e"]))
    pool = rw.load_beagle_spectator_csv(
        str(csv), electron_beam=(0.0, 0.0, -9.0, 9.0))
    assert pool.s == pytest.approx(4.0 * 9.0 * 130.0)
    assert pool.x == pytest.approx(x, rel=1e-9)
    assert pool.q2 == pytest.approx(q2, rel=1e-9)
    assert pool.extra["ievt"].tolist() == [0, 1]
    assert np.all(pool.weight == 1.0)


def test_beagle_rows_drop_events_without_a_scattered_electron():
    rows = _fake_beagle_rows(np.array([0.05]), np.array([20.0]),
                             np.array([1.0]))
    only_beam = rows[rows["kind"] == "B"]
    ievt, xr, _, _, _ = rw.sample_from_beagle_rows(
        only_beam, only_beam["kind"], electron_beam=(0.0, 0.0, -9.0, 9.0))
    assert ievt.size == 0 and xr.size == 0


# --- the weights ---------------------------------------------------------


def test_weights_carry_the_injected_asymmetry(vector_rw):
    """<W>_pool for a helicity-flip pair is 1 -+ P_e P_z <A_par>, which is
    what makes the flip estimator exact -- checked against the asymmetry
    computed independently through polli_fastsim.asymmetries."""
    plan = bk.helicity_flip_plan(1.0, PZ, PE)
    truth = vector_rw.truth_a_parallel()
    wp = vector_rw.mean_weight(plan.categories[0])
    wm = vector_rw.mean_weight(plan.categories[1])
    assert plan.categories[0].lam_e == +1
    assert (wp - wm) / (wp + wm) == pytest.approx(PE * PZ * truth, rel=1e-9)
    assert wp + wm == pytest.approx(2.0, rel=1e-9)   # no tensor sector


def test_tensor_thirds_weights_are_lumi_neutral(tensor_rw):
    """(+, +, -2) tensor fills: the three mean weights sum to 3, so the
    thirds denominator is the unpolarized rate and the estimator returns
    <A_zz> with no dilution."""
    plan = bk.tensor_thirds_plan(PZ, PZZ)
    ws = [tensor_rw.mean_weight(c) for c in plan.categories]
    assert sum(ws) == pytest.approx(3.0, rel=1e-9)
    truth = tensor_rw.truth_azz()
    assert est.azz_thirds(ws[0], ws[1], ws[2], PZZ) == pytest.approx(
        truth, rel=1e-9)


def test_weight_cache_is_share_blind(tensor_rw):
    """`bookkeeping.with_offset` moves luminosity, not physics: the same
    category under a different share must return the same weights, and a
    different population must not hit the cache."""
    plan = bk.tensor_thirds_plan(PZ, PZZ)
    offset = bk.with_offset(plan, "azz0", 0.05)
    for a, b in zip(plan.categories, offset.categories):
        assert tensor_rw.mean_weight(a) == tensor_rw.mean_weight(b)
    other = bk.tensor_thirds_plan(PZ, 0.3)
    assert tensor_rw.mean_weight(other.categories[0]) != pytest.approx(
        tensor_rw.mean_weight(plan.categories[0]))


def test_negative_weight_is_refused():
    """A weight below zero is not a probability; the reweighter must say
    so rather than resample max(W, 0) and dilute the modulation."""
    base = InclusiveKernel(beams.LI6)
    kern = InclusiveKernel(beams.LI6,
                           delta_func=lambda x, q2, f1: 50.0 * f1)
    rwt = rw.ModeWReweighter(kern, _pool(base, 5, n=4000))
    cat = bk.transverse_tensor_plan(PZZ).categories[0]
    with pytest.raises(ValueError, match="negative Mode-W weight"):
        rwt.event_weights(cat)


# --- resampling ----------------------------------------------------------


def test_resampled_events_reproduce_the_polarized_phi(tensor_rw):
    """The cos 2phi' channel: the pool's phi is flat, so the modulation in
    the resampled events is entirely the reweighter's.  Recovered with the
    ordinary `estimators.cos2phi_fit`."""
    kern = InclusiveKernel(beams.LI6,
                           delta_func=lambda x, q2, f1:
                           toy_delta_gluon(x, q2, f1, scale=3.0))
    pool = _pool(InclusiveKernel(beams.LI6), 9)
    rwt = rw.ModeWReweighter(kern, pool)
    plan = bk.transverse_tensor_plan(PZZ, phi_s=0.7)
    cat = plan.categories[0]
    truth = rwt.truth_cos2phi()
    assert abs(truth) > 0.2                       # a signal worth fitting
    rng = np.random.default_rng(31)
    ev = rwt.draw_category(cat, 400_000, rng)
    amp = est.cos2phi_fit(ev["phi"] - cat.phi_s, PZZ)
    # bootstrap resampling cannot beat the pool's OWN resolution on the
    # amplitude, sqrt(2/N_pool)/P_zz, however many events are drawn
    tol = 4.0 * np.sqrt(2.0 / len(pool)) / PZZ
    assert abs(amp - truth) < tol
    assert amp == pytest.approx(truth, rel=0.12)
    # and the unpolarized pool itself carries none of it
    flat = est.cos2phi_fit(pool.phi - cat.phi_s, PZZ)
    assert abs(flat) < tol and abs(flat) < 0.1 * abs(truth)


def test_draw_category_returns_pool_events(vector_rw):
    cat = bk.helicity_flip_plan(1.0, PZ, PE).categories[0]
    ev = vector_rw.draw_category(cat, 500, np.random.default_rng(4))
    assert ev["x"].size == 500 and ev["category"] == cat.name
    assert ev["lam_e"] == +1
    assert np.all(ev["x"] == vector_rw.sample.x[ev["idx"]])
    assert ev["y"] == pytest.approx(ev["q2"] / (vector_rw.s * ev["x"]))
    empty = vector_rw.draw_category(cat, 0, np.random.default_rng(4))
    assert empty["x"].size == 0


def test_pseudo_experiment_counts_track_the_shares(vector_rw):
    plan = bk.helicity_flip_plan(1.0, PZ, PE, rel_lumi_offset=0.2)
    mus = vector_rw.expected_counts(plan, 100_000)
    # shares 0.5(1+0.2) and 0.5 -> the offset is visible in the means
    assert mus["apar+"] / mus["apar-"] == pytest.approx(
        1.2 * vector_rw.mean_weight(plan.categories[0])
        / vector_rw.mean_weight(plan.categories[1]), rel=1e-9)
    ev, _ = vector_rw.run_pseudo_experiment(
        plan, 100_000, rng=np.random.default_rng(2), draw_events=False)
    counts = rw.counts_of(ev)
    assert set(counts) == {"apar+", "apar-"}
    assert abs(counts["apar+"] - mus["apar+"]) < 6.0 * np.sqrt(mus["apar+"])


# --- the acceptance criterion -------------------------------------------


def _apar_scan(rwt, ntrials=NTRIALS, seed=101, draw_events=False):
    plan = bk.helicity_flip_plan(1.0, PZ, PE)
    return rw.pull_scan(
        rwt, plan, 60_000,
        lambda c: est.apar_flip(c["apar+"], c["apar-"], PE, PZ),
        rwt.truth_a_parallel(), lambda n: err_a_parallel(n, PE, PZ),
        ntrials=ntrials, seed=seed, draw_events=draw_events)


def _azz_scan(rwt, ntrials=NTRIALS, seed=102, draw_events=False):
    plan = bk.tensor_thirds_plan(PZ, PZZ)
    return rw.pull_scan(
        rwt, plan, 50_000,
        lambda c: est.azz_thirds(c["azz+"], c["azz-"], c["azz0"], PZZ),
        rwt.truth_azz(), lambda n: err_azz(n, PZZ),
        ntrials=ntrials, seed=seed, draw_events=draw_events)


def test_apar_pull_closure(vector_rw):
    pulls, vals, _ = _apar_scan(vector_rw)
    assert pulls.size == NTRIALS
    assert abs(pulls.mean()) < PULL_MEAN_MAX
    assert abs(pulls.std(ddof=1) - 1.0) < PULL_WIDTH_TOL
    assert vals.mean() == pytest.approx(
        vector_rw.truth_a_parallel(),
        abs=4.0 * vals.std(ddof=1) / np.sqrt(vals.size))


def test_azz_pull_closure(tensor_rw):
    pulls, vals, _ = _azz_scan(tensor_rw)
    assert pulls.size == NTRIALS
    assert abs(pulls.mean()) < PULL_MEAN_MAX
    assert abs(pulls.std(ddof=1) - 1.0) < PULL_WIDTH_TOL
    assert vals.mean() == pytest.approx(
        tensor_rw.truth_azz(),
        abs=4.0 * vals.std(ddof=1) / np.sqrt(vals.size))


def test_pull_closure_survives_event_resampling(vector_rw, tensor_rw):
    """The same gate with the events actually drawn from the pool, not
    only their Poisson counts -- the path a reco-level closure would use."""
    for pulls, _, _ in (_apar_scan(vector_rw, ntrials=NTRIALS_DRAW, seed=11,
                                   draw_events=True),
                        _azz_scan(tensor_rw, ntrials=NTRIALS_DRAW, seed=12,
                                  draw_events=True)):
        assert abs(pulls.mean()) < PULL_MEAN_MAX
        assert abs(pulls.std(ddof=1) - 1.0) < PULL_WIDTH_TOL


# --- the streamed sample (skipped without the cached CSV) ----------------

_CSV = os.environ.get("POLLIGEN_BEAGLE_CSV", "")


@pytest.mark.skipif(not (_CSV and pathlib.Path(_CSV).exists()),
                    reason="set POLLIGEN_BEAGLE_CSV to a cached "
                           "dump_spectators.py CSV (tools/beagle/README.md)")
def test_streamed_beagle_pool_closure():
    """The same acceptance gate on the official BeAGLE e+d sample.

    Stream it once with the recipe in `load_beagle_spectator_csv` and
    point POLLIGEN_BEAGLE_CSV at the CSV; the test never touches the
    network itself.
    """
    pool = rw.load_beagle_spectator_csv(
        _CSV, electron_beam=(0.0, 0.0, -9.0, 9.0))
    assert len(pool) > 2000
    pool = pool.select(pool.in_scenario(fom.Scenario(), electron_energy=9.0))
    # unpolarized in phi, as an external sample must be
    assert abs(np.mean(np.cos(2.0 * pool.phi))) < 5.0 / np.sqrt(2 * len(pool))

    kern = InclusiveKernel(beams.DEUTERON)
    pulls, _, _ = _apar_scan(rw.ModeWReweighter(kern, pool),
                             ntrials=NTRIALS, seed=201)
    assert abs(pulls.mean()) < PULL_MEAN_MAX
    assert abs(pulls.std(ddof=1) - 1.0) < PULL_WIDTH_TOL

    kern2 = InclusiveKernel(beams.DEUTERON, b1_func=toy_b1)
    pulls2, _, _ = _azz_scan(rw.ModeWReweighter(kern2, pool),
                             ntrials=NTRIALS, seed=202)
    assert abs(pulls2.mean()) < PULL_MEAN_MAX
    assert abs(pulls2.std(ddof=1) - 1.0) < PULL_WIDTH_TOL
