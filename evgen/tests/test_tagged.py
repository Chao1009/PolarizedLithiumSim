"""Step-5.B gates (plans/05 §5.4): tagged-mode spin (x) spectator physics.

Structural identities that must hold for any two-cluster spin model
(normalization, isotropy sums, pure-wave CG limits), the analytic 7Li
polarimetry/forward-limit predictions, the 6Li embedded-deuteron dilution
and S-wave reduction to the inclusive master formula, the deuteron-limit
S/D tensor mechanism (Cosyn-Weiss style), and sampler/boost consistency
against the fastsim spectator module.
"""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polligen import bookkeeping as bk  # noqa: E402
from polligen import tagged  # noqa: E402
from polligen.spin import clebsch_gordan, m_values  # noqa: E402
from polligen.xsec import EventSpinState, InclusiveKernel  # noqa: E402

from polli_fastsim import beams, spectator  # noqa: E402
from polli_fastsim.asymmetries import azz  # noqa: E402
from polli_fastsim.polarized import toy_b1  # noqa: E402


@pytest.fixture(scope="module")
def li6():
    return tagged.TaggedModel(tagged.li6_alpha_channel())


@pytest.fixture(scope="module")
def li7():
    return tagged.TaggedModel(tagged.li7_alpha_channel())


@pytest.fixture(scope="module")
def deut():
    return tagged.TaggedModel(tagged.deuteron_channel())


# --- structural identities -------------------------------------------------


def test_normalization_all_channels(li6, li7, deut):
    for model in (li6, li7, deut):
        for m in m_values(model.channel.j_ion):
            assert model.norm(m) == pytest.approx(1.0, abs=2e-3)


def test_state_sum_isotropic(li6, li7, deut):
    """sum_M n_M(k, c) must be c-independent (unpolarized beam isotropy
    about the axis; L-interference cancels in the M sum)."""
    for model in (li6, li7, deut):
        tot = sum(model.n_of_kc(m) for m in m_values(model.channel.j_ion))
        spread = np.ptp(tot, axis=1) / np.maximum(tot.mean(axis=1), 1e-300)
        assert spread.max() < 1e-10


def test_pure_s_wave_is_m_independent():
    ch = tagged.li6_alpha_channel(p_d=0.0)
    model = tagged.TaggedModel(ch)
    n = {m: model.n_of_kc(m) for m in (1.0, 0.0, -1.0)}
    for m in (0.0, -1.0):
        np.testing.assert_allclose(n[m], n[1.0], rtol=1e-12)
    # S-wave: struck deuteron spin = ion spin exactly
    p = model.population_integrated(1.0)
    np.testing.assert_allclose(p, [1.0, 0.0, 0.0], atol=1e-12)


def test_pure_d_wave_forward_ratio():
    """theta_k = 0 limit of a pure L=2 wave: only m_L = 0 survives, so
    n_{+-1}/n_0 = CG(2 0 1 1|1 1)^2 / CG(2 0 1 0|1 0)^2 = (1/10)/(4/10)
    = 1/4, and the thirds combination gives A_zz^wf(theta=0) = -1 exactly
    (maximal tensor analyzing power of the D-wave forward direction)."""
    ch = tagged.TaggedChannel(spectator.LI6_ALPHA_TAG, 1.0, 1.0, 0.0, 1.0,
                              (tagged.Wave(2, 1.0),), beams.DEUTERON, "pureD")
    model = tagged.TaggedModel(ch, nc=4001)
    ic = np.argmax(model.c)  # closest grid point to c = 1 (sin^2 ~ 5e-4)
    n = {m: model.n_of_kc(m)[:, ic] for m in (1.0, 0.0, -1.0)}
    expected = (clebsch_gordan(2, 0, 1, 1, 1, 1) ** 2
                / clebsch_gordan(2, 0, 1, 0, 1, 0) ** 2)
    np.testing.assert_allclose(n[1.0] / n[0.0], expected, rtol=5e-3)
    assert expected == pytest.approx(0.25, abs=1e-12)
    a = (n[1.0] + n[-1.0] - 2 * n[0.0]) / (n[1.0] + n[-1.0] + n[0.0])
    np.testing.assert_allclose(a, -1.0, atol=3e-3)


def test_density_matches_independent_transcription(li6):
    """n_M and p(m_S|M) against a test-local, loop-based CG sum."""
    model = li6
    ch = model.channel
    kappa = ch.base.kappa
    k = np.array([0.05, 0.15, 0.30, 0.45])
    c = np.array([-0.7, 0.1, 0.9])
    # test-local normalized radials
    rad = {}
    kg = model.k
    for w in ch.waves:
        psi = w.radial(kg, kappa)
        norm = np.sqrt(getattr(np, "trapezoid", np.trapz)(psi**2 * kg**2, kg))
        rad[w.l] = lambda kk, w=w, norm=norm: (np.sqrt(w.prob)
                                               * w.radial(kk, kappa) / norm)
    for M in (1.0, 0.0, -1.0):
        for kk in k:
            for cc in c:
                # evaluate the transcription at the nearest cell center
                # (the module's lookup point) -> exact identity
                ik = np.clip(np.searchsorted(model.k, kk) - 1, 0,
                             model.k.size - 2)
                ic = np.clip(np.searchsorted(model.c, cc) - 1, 0,
                             model.c.size - 2)
                kc, ccc = model.k[ik], model.c[ic]
                dens = 0.0
                for m_s in (1.0, 0.0, -1.0):
                    amp = 0.0
                    for w in ch.waves:
                        m_l = M - m_s
                        if abs(m_l) > w.l:
                            continue
                        amp += (rad[w.l](kc)
                                * clebsch_gordan(w.l, m_l, 1.0, m_s, 1.0, M)
                                * tagged.theta_lm(w.l, int(m_l),
                                                  np.array([ccc]))[0])
                    dens += amp * amp
                got = model.n_of_kc(M, np.array([kk]), np.array([cc]))[0]
                assert got == pytest.approx(dens, rel=1e-9)


# --- 7Li: polarimetry + forward limit --------------------------------------


def test_li7_stretched_state_angular_shape(li7):
    """n_{3/2}(c) ~ sin^2(theta): pure |Y_1^1|^2."""
    n = li7.n_of_kc(1.5)
    prof = n[50] / n[50].max()
    expected = (1.0 - li7.c**2) / (1.0 - li7.c**2).max()
    np.testing.assert_allclose(prof, expected, atol=1e-12)


def test_li7_p2_moments_and_polarimeter(li7):
    # midpoint-rule c-grid: discretization ~ (dc)^2 ~ 1e-4
    tol = 3e-4
    assert li7.p2_moment(1.5) == pytest.approx(-0.2, abs=tol)
    assert li7.p2_moment(0.5) == pytest.approx(+0.2, abs=tol)
    assert li7.p2_moment(-0.5) == pytest.approx(+0.2, abs=tol)
    assert li7.p2_moment(-1.5) == pytest.approx(-0.2, abs=tol)
    # polarimeter: <P2> = -T/5 for any fill (T = normalized tensor moment)
    from polligen.spin import moments_along_axis, spin32_populations
    pops = spin32_populations(0.5, 0.4, 0.1)
    t = moments_along_axis(1.5, pops)[1]
    assert li7.p2_moment_mixture(pops) == pytest.approx(-t / 5.0, abs=tol)


def test_li7_triton_polarization_forward_limit(li7):
    """P_t(M=3/2) = 1, P_t(M=1/2) = 1/3; with the triton's own effective
    proton polarization 0.86 this gives P_p(7Li) ~ 0.86 for a stretched
    fill -- the plans/05 forward-limit gate (VMC: 0.866)."""
    p32 = li7.population_integrated(1.5)
    np.testing.assert_allclose(p32, [1.0, 0.0], atol=1e-12)
    p12 = li7.population_integrated(0.5)
    np.testing.assert_allclose(p12, [2.0 / 3.0, 1.0 / 3.0], atol=3e-4)
    p_t = (p12 * np.array([1.0, -1.0])).sum()
    assert p_t == pytest.approx(1.0 / 3.0, abs=6e-4)
    p_p_7li = 1.0 * tagged.TRITON.eff_pol_p  # stretched fill
    assert abs(p_p_7li - 0.866) < 0.02  # within the D-state band


def test_li7_khat_resolved_triton_polarization(li7):
    """P_t(theta; M=1/2) = (5c^2 - 1)/(3c^2 + 1): +1 forward, -1 at 90 deg."""
    p = li7.struck_populations(0.5)
    c = li7.c
    pol = (p[0] - p[1])  # m_t = +1/2 minus -1/2, any k row (P-wave only)
    expected = (5 * c**2 - 1.0) / (3 * c**2 + 1.0)
    np.testing.assert_allclose(pol[50], expected, atol=1e-10)


# --- 6Li: embedded-deuteron dilution + inclusive reduction -----------------


def test_li6_embedded_deuteron_dilutions(li6):
    """Vector dilution 1 - (3/2) P_D = 0.87 with the default P_D --
    the b1_li6_from_deuteron scaling recovered from the tagged model."""
    got = li6.vector_dilution()
    assert got == pytest.approx(1.0 - 1.5 * tagged.P_D_LI6, abs=1e-3)
    assert got == pytest.approx(0.87, abs=2e-3)
    # tensor dilution of the stretched state (embedded-d Pzz)
    tzz = li6.tensor_dilution()
    assert 0.7 < tzz < 1.0  # diluted but not destroyed
    # deuteron control: standard 1 - (3/2) w_D
    dm = tagged.TaggedModel(tagged.deuteron_channel())
    assert dm.vector_dilution() == pytest.approx(
        1.0 - 1.5 * tagged.P_D_DEUTERON, abs=1e-3)


def test_li6_s_wave_reduces_to_inclusive_deuteron():
    """P_D = 0: k-integrated tagged tensor asymmetry == inclusive
    polarized-deuteron azz bin-wise (tagged -> inclusive integration)."""
    model = tagged.TaggedModel(tagged.li6_alpha_channel(p_d=0.0))
    kern = InclusiveKernel(beams.DEUTERON, b1_func=toy_b1)
    config = beams.default_configs("6Li")[1]
    s = config.sqrt_s_per_nucleon**2
    x = np.array([0.05, 0.15, 0.35])
    q2 = np.array([5.0, 12.0, 30.0])
    y = q2 / (s * x)
    t = kern.tables(x, q2)
    # S-wave: p(m_d | M) = delta_{m_d, M} -> tagged sigma_m = inclusive
    sig = {}
    for m in (1.0, 0.0, -1.0):
        p = model.population_integrated(m)
        w_mix = 0.0
        for pj, m_d in zip(p, (1.0, 0.0, -1.0)):
            if pj <= 0:
                continue
            w, _, _ = kern.amplitudes(t, x, q2, s,
                                      EventSpinState(0, 0.0, 1.0, m_d))
            w_mix = w_mix + pj * (1.0 + w)
        sig[m] = w_mix
    measured = ((sig[1.0] + sig[-1.0] - 2 * sig[0.0])
                / (sig[1.0] + sig[-1.0] + sig[0.0]))
    np.testing.assert_allclose(measured,
                               azz(t["b1"], t["f1"], t["f2"], x, y),
                               rtol=1e-10)


# --- deuteron control: the Cosyn-Weiss tagged tensor mechanism -------------


def test_deuteron_tagged_azz_wf_shape(deut):
    """S/D interference: A_zz^wf(k, theta=90deg) vanishes at k -> 0 (D-wave
    threshold), grows to O(1) in the 0.25-0.5 GeV/c region."""
    ic = np.argmin(np.abs(deut.c))  # cos theta ~ 0
    n = {m: deut.n_of_kc(m)[:, ic] for m in (1.0, 0.0, -1.0)}
    a = (n[1.0] + n[-1.0] - 2 * n[0.0]) / (n[1.0] + n[-1.0] + n[0.0])
    k = deut.k
    assert abs(a[k < 0.02].max()) < 0.02          # threshold suppression
    # O(1) asymmetry peaking near 300 MeV/c (CW mechanism); peak location
    peak = k[np.argmax(np.abs(a))]
    assert np.abs(a[(k > 0.25) & (k < 0.5)]).max() > 0.45
    assert 0.2 < peak < 0.45
    # k-integrated over all angles the wf tensor asymmetry vanishes
    # (L-orthogonality; midpoint-rule residual ~ dc^2) -- the observable
    # lives in the angular structure
    for m in (1.0, 0.0):
        nk = (deut.n_of_kc(m) * 2 * np.pi * deut.dc).sum(axis=1)
        nk1 = (deut.n_of_kc(1.0) * 2 * np.pi * deut.dc).sum(axis=1)
        np.testing.assert_allclose(nk, nk1, rtol=3e-4)


def test_deuteron_struck_neutron_polarization(deut):
    """Pair decomposition: triplet m_S = +-1 -> fully polarized neutron,
    m_S = 0 -> unpolarized."""
    for m_s, expected in ((1.0, +1.0), (0.0, 0.0), (-1.0, -1.0)):
        pair = deut.pair_populations(m_s)     # (m_n, m_p) with m ordering
        p_n = pair.sum(axis=1)                # marginal over proton spin
        pol = 2.0 * (p_n * m_values(0.5)).sum()
        assert pol == pytest.approx(expected, abs=1e-12)


# --- sampler + boost + routing ---------------------------------------------


@pytest.fixture(scope="module")
def li6_sampler(li6):
    from polli_fastsim import fom
    kern = InclusiveKernel(beams.DEUTERON, b1_func=toy_b1)
    config = beams.default_configs("6Li")[1]
    return tagged.TaggedSampler(li6, kern, config, fom.Scenario(),
                                nx=24, nq2=18, x_range=(1e-3, 0.5),
                                q2_range=(1.5, 100.0))


def test_sampled_kc_matches_density(li6):
    rng = np.random.default_rng(21)
    n = 200_000
    k, c, _ = li6.sample_kc(1.0, 1.0, n, rng)
    # cos-theta marginal: histogram bins aligned with the model c-cells
    # (16 bins x 6 cells of the 96-cell grid)
    prob_c = (li6._amp2_table(1.0)[0] * li6.k[:, None] ** 2).sum(axis=0)
    prob_c = prob_c / prob_c.sum()
    edges = np.linspace(-1.0, 1.0, 17)
    counts, _ = np.histogram(c, bins=edges)
    cell_bin = np.digitize(li6.c, edges) - 1
    expected = np.array([prob_c[cell_bin == b].sum()
                         for b in range(16)]) * n
    z = (counts - expected) / np.sqrt(np.maximum(expected, 1.0))
    assert np.abs(z).max() < 5.0


def test_tagged_events_complete_and_routed(li6_sampler):
    cat = bk.SpinCategory("t+", 1.0, (1.0, 0.0, 0.0))
    ev = li6_sampler.sample_category(cat, n=20_000,
                                     rng=np.random.default_rng(22))
    assert np.all(ev["m_ion"] == 1.0)
    assert set(np.unique(ev["m_struck"])) <= {-1.0, 0.0, 1.0}
    # spectator kinematics present and physical
    assert np.all(ev["pT"] >= 0) and np.all(ev["p_lab"] > 0)
    assert np.all((ev["route"] >= 0) & (ev["route"] <= 4))
    # 6Li alpha spectator is BEAM-BLIND: R ~ 1 (rigidity 2 p_u vs beam
    # 2 p_u) -- only the pT tail reaches the Roman Pots (fastsim headline
    # 3-9% with high-acceptance optics, beta = 0.3 central)
    r_med = np.median(ev["R"])
    assert 0.97 < r_med < 1.03
    frac_tail = float(np.mean(ev["route"] == 4))
    frac_rp_main = float(np.mean(ev["route"] == 1))
    assert 0.005 < frac_tail < 0.20
    assert frac_rp_main < 0.05
    assert float(np.mean(ev["route"] == 0)) > 0.5  # mostly lost


def test_boost_matches_fastsim_spectator(li6_sampler):
    """S-wave limit: tagged pT spectrum == spectator.py's (same density,
    same boost), cross-checked statistically."""
    model = tagged.TaggedModel(tagged.li6_alpha_channel(p_d=0.0))
    rng = np.random.default_rng(23)
    k, c, phi = model.sample_kc(1.0, 1.0, 100_000, rng)
    lab = tagged.boost_spectator(model.channel, k, c, phi,
                                 li6_sampler.config.ion_momentum_per_nucleon)
    ref = spectator.spectator_lab_kinematics(
        spectator.LI6_ALPHA_TAG,
        li6_sampler.config.ion_momentum_per_nucleon, n=100_000,
        rng=np.random.default_rng(24))
    for key in ("pT", "R", "xL"):
        q_new = np.percentile(lab[key], [25, 50, 75])
        q_ref = np.percentile(ref[key], [25, 50, 75])
        np.testing.assert_allclose(q_new, q_ref, rtol=0.02)


def test_tagged_rate_asymmetry_matches_analytic(li6_sampler):
    """sigma_tot per pure ion state reproduces the diluted inclusive azz
    (S+D model: DIS-side b1 modulation with the k-integrated embedded-d
    populations; wf part integrates to zero over 4pi)."""
    model = li6_sampler.model
    sig = {}
    for m in (1.0, 0.0, -1.0):
        pops = tuple(1.0 if np.isclose(mm, m) else 0.0
                     for mm in (1.0, 0.0, -1.0))
        sig[m] = li6_sampler.sigma_tot_pb(
            bk.SpinCategory("s%g" % m, 1.0, pops))
    measured = ((sig[1.0] + sig[-1.0] - 2 * sig[0.0])
                / (sig[1.0] + sig[-1.0] + sig[0.0]))
    # analytic: dilution(tensor) * sigma-weighted inclusive azz of the
    # struck deuteron over the accepted DIS phase space
    inner = li6_sampler.inner
    t = inner.tables
    y = inner.q2_cells / (inner.s * inner.x_cells)
    azz_w = np.average(azz(t["b1"], t["f1"], t["f2"], inner.x_cells, y),
                       weights=inner.xsec_flat)
    expected = model.tensor_dilution() * azz_w
    assert measured == pytest.approx(expected, rel=2e-2)
