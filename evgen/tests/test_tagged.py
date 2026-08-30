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

_trapezoid = getattr(np, "trapezoid", None) or np.trapz

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
        norm = np.sqrt(_trapezoid(psi**2 * kg**2, kg))
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
    # 6 is the over-rigid RP-inner branch added on 2026-08-28 (plans/09
    # B1); 5 is the neutral ZDC code and never appears for a charged tag
    from polli_fastsim import farforward as _ff
    assert np.all(np.isin(ev["route"], list(_ff.ROUTE_LABELS)))
    assert np.all((ev["route"] >= 0) & (ev["route"] <= 6)
                  & (ev["route"] != 5))
    # 6Li alpha spectator is BEAM-BLIND: R ~ 1 (rigidity 2 p_u vs beam
    # 2 p_u), so it reaches the Roman Pots only through the pT tail
    # outside the near-beam envelope or through the off-rigidity slice
    # below R = 0.95.  Since 2026-08-28 the sampler's own route column is
    # at the CONFIGURATION's Yellow Report high-acceptance envelope
    # (1.80 x 1.80 mrad at 10 x 99.5), not at the legacy 73 microrad: the
    # tail collapses from ~2.5% to ~0.05% and the off-rigidity slice --
    # which no envelope touches -- is what is left (plans/09 B2).
    r_med = np.median(ev["R"])
    assert 0.97 < r_med < 1.03
    frac_tail = float(np.mean(ev["route"] == 4))
    frac_rp_main = float(np.mean(ev["route"] == 1))
    assert 0.0 < frac_tail < 0.005
    assert 0.01 < frac_rp_main < 0.05
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


# --- per-configuration optics and the lab azimuth (plans/09 B2) ------------
#
# Until 2026-08-28 the tagged sampler routed through the retired
# proton-derived HIGH_ACCEPTANCE (73 microrad, isotropic, the same at every
# energy) and passed no azimuth, so the rectangular near-beam envelope
# degenerated to a circle at n sigma_h.  Both are fixed; these pin the fix.


def test_phi_spec_is_the_lab_azimuth_of_the_transverse_momentum():
    """`boost_spectator` must expose the azimuth the pots see -- and it is
    NOT the DIS azimuth, which is why the key is not "phi"."""
    ch = tagged.li6_alpha_channel()
    rng = np.random.default_rng(31)
    k = rng.uniform(0.02, 0.6, 5000)
    c = rng.uniform(-1.0, 1.0, 5000)
    phi_k = rng.uniform(0.0, 2.0 * np.pi, 5000)
    lab = tagged.boost_spectator(ch, k, c, phi_k, 99.5)
    np.testing.assert_allclose(lab["phi_spec"],
                               np.arctan2(lab["ky"], lab["kx"]), atol=1e-12)
    # the longitudinal boost does not touch the transverse plane, so with
    # an untilted spin axis the lab azimuth IS the spin-frame one
    np.testing.assert_allclose(np.cos(lab["phi_spec"]), np.cos(phi_k),
                               atol=1e-12)
    # and pT / phi_spec reconstruct (kx, ky)
    np.testing.assert_allclose(lab["pT"] * np.cos(lab["phi_spec"]),
                               lab["kx"], atol=1e-12)
    # a tilted quantization axis is rotated BEFORE the boost, so the lab
    # azimuth is no longer phi_k
    tilt = tagged.boost_spectator(ch, k, c, phi_k, 99.5, theta_s=0.5)
    assert np.mean(np.abs(np.cos(tilt["phi_spec"]) - np.cos(phi_k))) > 0.05
    np.testing.assert_allclose(tilt["phi_spec"],
                               np.arctan2(tilt["ky"], tilt["kx"]), atol=1e-12)


def test_routing_uses_the_rectangular_envelope():
    """Mirror of fastsim/tests/test_optics_20260828.py::
    test_rectangular_envelope_in_the_routing, through the tagged sampler's
    own routing call: a beam-rigidity fragment at 3 mrad clears the
    2.2 x 3.8 mrad envelope of 5 x 40.8 horizontally and not vertically."""
    from polli_fastsim import farforward as ff
    o = ff.yr_optics(beams.default_configs("6Li")[0])
    assert o.envelope == pytest.approx((2.20e-3, 3.80e-3), abs=1e-5)
    assert ff.route_charged(1.0, 3.0e-3, 0.1, o, phi=0.0) == 4
    assert ff.route_charged(1.0, 3.0e-3, 0.1, o, phi=np.pi / 2) == 0
    # what dropping the azimuth would do: the circle at n sigma_h accepts
    # the vertical fragment too
    assert ff.route_charged(1.0, 3.0e-3, 0.1, o) == 4


def test_sampler_defaults_to_the_configuration_optics():
    """No module-level default can know the beam: the default envelope is
    this configuration's Yellow Report high-acceptance one, and the legacy
    constants stay reachable by passing them explicitly."""
    from polli_fastsim import farforward as ff, fom
    kern = InclusiveKernel(beams.DEUTERON, b1_func=toy_b1)
    model = tagged.TaggedModel(tagged.li6_alpha_channel())
    for cfg in beams.default_configs("6Li"):
        s = tagged.TaggedSampler(model, kern, cfg, fom.Scenario(),
                                 nx=8, nq2=6)
        assert s.optics.envelope == ff.yr_optics(cfg).envelope
        assert s.optics.sigma_theta_v is not None      # anisotropic-capable
    s = tagged.TaggedSampler(model, kern, beams.default_configs("6Li")[1],
                             fom.Scenario(), optics=ff.HIGH_ACCEPTANCE,
                             nx=8, nq2=6)
    assert s.optics is ff.HIGH_ACCEPTANCE


def test_alpha_tag_acceptance_per_optics_at_10x99():
    """The headline B2 numbers, on the tagged sampler at 10 x 99.5: 2.5% at
    the Yellow Report high-acceptance optics against 25% at the tagging
    optics, and the whole difference is the near-beam tail.  The 25% was
    30% until 2026-08-29, when `tagging_optics_point` began pricing the
    dispersive envelope term on each configuration's own (R12, D): at
    10 x 100 that term is 39% larger than the 18 x 275 pair it replaced,
    the optimum de-squeeze falls from 175.6 to 164.1, and the horizontal
    envelope opens from 0.166 to 0.192 mrad.  The circular cut the
    azimuth-less call applies would still read ~1.6x more."""
    from polli_fastsim import farforward as ff, fom
    kern = InclusiveKernel(beams.DEUTERON, b1_func=toy_b1)
    model = tagged.TaggedModel(tagged.li6_alpha_channel())
    cfg = beams.default_configs("6Li")[1]
    s = tagged.TaggedSampler(model, kern, cfg, fom.Scenario(),
                             nx=16, nq2=12)
    cat = bk.SpinCategory("flat", 1.0, (1 / 3., 1 / 3., 1 / 3.))
    ev = s.sample_category(cat, n=200_000, rng=np.random.default_rng(32))
    acc = {}
    for label, o in (("yr", ff.yr_optics(cfg, "high-acceptance")),
                     ("tag", ff.tagging_optics(cfg))):
        r = ff.route_charged(ev["R"], ev["theta"], ev["pT"], o,
                             phi=ev["phi_spec"])
        acc[label] = float(np.mean((r == 1) | (r == 4)))
    assert acc["yr"] == pytest.approx(0.025, rel=0.10)
    assert acc["tag"] == pytest.approx(0.254, rel=0.10)
    # the off-rigidity R < 0.95 window slice is optics-independent and is
    # all that survives at the Yellow Report optics
    r_yr = ff.route_charged(ev["R"], ev["theta"], ev["pT"],
                            ff.yr_optics(cfg), phi=ev["phi_spec"])
    assert float(np.mean(r_yr == 1)) == pytest.approx(0.024, rel=0.10)
    # dropping the azimuth overstates the tagging optics by ~1.7x
    r_circ = ff.route_charged(ev["R"], ev["theta"], ev["pT"],
                              ff.tagging_optics(cfg))
    assert float(np.mean((r_circ == 1) | (r_circ == 4))) > 1.5 * acc["tag"]


def test_acceptance_weighted_curve_reduces_to_the_90_degree_curve(li6):
    """The overlay of money plot 4.  Weights concentrated at cos theta_k =
    0 must return the analytic 90 degree curve exactly; a real acceptance
    table must not, and at the Yellow Report optics it must come out with
    the OPPOSITE sign at k ~ 0.3 GeV/c, which is the defect the weighting
    corrects."""
    from polli_fastsim import farforward as ff
    ic = int(np.argmin(np.abs(li6.c)))
    ref = tagged.azz_tensor_curve(li6, ic)
    w = np.zeros((li6.k.size, li6.c.size))
    w[:, ic] = 1.0
    np.testing.assert_allclose(tagged.azz_tensor_curve(li6, weights=w), ref,
                               rtol=1e-12)
    # a uniform weight is the 4pi average: the L cross terms integrate
    # away and CG completeness makes the c-integral of n_M the same for
    # every M, so the tensor combination vanishes at every k -- to the
    # accuracy of the midpoint rule on the 96-cell c grid, which is 3e-5
    flat = tagged.azz_tensor_curve(li6, weights=np.ones_like(w))
    assert np.abs(flat).max() < 1e-3

    cfg = beams.default_configs("6Li")[1]
    j = int(np.argmin(np.abs(li6.k - 0.30)))
    eps_yr = tagged.acceptance_weights(li6, cfg, ff.yr_optics(cfg))
    eps_tag = tagged.acceptance_weights(li6, cfg, ff.tagging_optics(cfg))
    a_yr = tagged.azz_tensor_curve(li6, weights=eps_yr)[j]
    a_tag = tagged.azz_tensor_curve(li6, weights=eps_tag)[j]
    assert ref[j] < -0.4                      # the 90 degree curve
    assert a_yr > +0.4                        # longitudinal acceptance
    assert -0.2 < a_tag < 0.0                 # transverse near-beam tail
    # eps is a probability, and the near-beam tail the tagging optics opens
    # is the transverse half of the sphere
    assert eps_tag.min() >= 0.0 and eps_tag.max() <= 1.0
    wt = eps_tag * li6.n_of_kc(1.0) * li6.k[:, None] ** 2
    wy = eps_yr * li6.n_of_kc(1.0) * li6.k[:, None] ** 2
    mean_abs_c = lambda ww: float((ww * np.abs(li6.c)).sum() / ww.sum())
    assert mean_abs_c(wt) < 0.5 < mean_abs_c(wy)


def test_the_alpha_tag_is_not_the_sub_0p6_histogram():
    """money_tagged_azz.py printed `acc` as the sum of its k < 0.6 GeV/c
    histogram, i.e. the accepted fraction TRUNCATED at the right edge of
    the panel, and six documents published that as the 6Li alpha tag.  The
    tail above 0.6 GeV/c is 9-11% of the accepted sample at the Yellow
    Report optics, so the two differ by that much and the truncated one is
    not an acceptance.  This pins the tail, and pins that the script now
    divides the UNBINNED accepted count by the generated one."""
    import importlib
    import pathlib as _pl
    import sys as _sys
    from polli_fastsim import farforward as ff, fom
    _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]
                            / "scripts"))
    mod = importlib.import_module("money_tagged_azz")

    kern = InclusiveKernel(beams.DEUTERON, b1_func=toy_b1)
    model = tagged.TaggedModel(tagged.li6_alpha_channel())
    cfg = beams.default_configs("6Li")[1]
    s = tagged.TaggedSampler(model, kern, cfg, fom.Scenario(), nx=16, nq2=12)
    cat = bk.SpinCategory("flat", 1.0, (1 / 3., 1 / 3., 1 / 3.))
    ev = s.sample_category(cat, n=200_000, rng=np.random.default_rng(32))
    r = ff.route_charged(ev["R"], ev["theta"], ev["pT"],
                         ff.yr_optics(cfg, "high-acceptance"),
                         phi=ev["phi_spec"])
    k_acc = ev["k"][(r == 1) | (r == 4)]
    assert 0.05 < float(np.mean(k_acc > 0.6)) < 0.15
    # the two candidate definitions, on the script's own machinery
    plan = bk.tensor_thirds_plan(0.7, 0.6)
    k_edges = np.linspace(0.0, 0.6, 13)
    menu = mod.optics_menu(cfg, "high-acceptance")
    folded, n_gen = mod.folded_asymmetry(s, plan, 60_000, k_edges, menu,
                                         np.random.default_rng(5),
                                         ff.yr_config_key(cfg))
    (_a, n, k_acc2), = folded.values()
    assert n.sum() < k_acc2.size            # the histogram truncates
    tag = k_acc2.size / n_gen
    assert tag == pytest.approx(0.025, rel=0.20)
    assert n.sum() / n_gen < 0.96 * tag     # by 4% or more, and here ~9%


def test_the_published_figure_stems_are_guarded_by_config_AND_optics():
    """The guard used to key on --config alone, so `--optics legacy` at the
    default configuration -- the reproduction command the manual documents
    -- silently overwrote the published PNG with the retired 73/164 microrad
    figure.  Only the default combination may claim the published stem."""
    import importlib
    import pathlib as _pl
    import sys as _sys
    _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]
                            / "scripts"))
    for name, base in (("money_tagged_azz", "money_tagged_azz_6Li"),
                       ("tagged_polarimetry_7li", "tagged_polarimetry_7Li")):
        stem = importlib.import_module(name).output_stem
        assert stem(base, "10x100", 1, "menu") == base
        for cfg, key, opt in ((1, "10x100", "legacy"),
                              (1, "10x100", "tagging"),
                              (1, "10x100", "high-acceptance"),
                              (0, "5x41", "menu"),
                              (2, "18x275", "legacy")):
            assert stem(base, key, cfg, opt) != base
        # and distinct runs never collide
        stems = {stem(base, k, c, o)
                 for c, k in ((0, "5x41"), (1, "10x100"), (2, "18x275"))
                 for o in ("menu", "legacy", "tagging", "high-acceptance",
                           "high-divergence")}
        assert len(stems) == 15

    # the two coherent scripts carry the same guard (2026-08-28 review):
    # money plot 6R is `--config 0 --optics tagging` with the ratio fit and
    # the published |t| edges, and Report 4's reach figure is the ratio fit
    class _A:
        def __init__(self, **kw):
            self.__dict__.update(kw)

    coh = importlib.import_module("money_cos2phi_coherent_reco").output_stem
    published = dict(config=0, optics="tagging", fit="ratio",
                     u_in_situ=False, t_edges=None)
    assert coh(_A(**published)) == "money_cos2phi_coherent_reco_6Li"
    variants = [dict(published, config=1), dict(published, optics="legacy"),
                dict(published, fit="likelihood"),
                dict(published, u_in_situ=True),
                dict(published, t_edges="0.006,0.05,0.25")]
    for kw in variants:
        assert coh(_A(**kw)) != "money_cos2phi_coherent_reco_6Li", kw
    assert len({coh(_A(**kw)) for kw in variants}) == len(variants)

    reach = importlib.import_module("nearbeam_reach_gain").output_stem
    assert reach(_A(fit="ratio", t_edges=None)) == "nearbeam_reach_gain_6Li"
    for kw in (dict(fit="likelihood", t_edges=None),
               dict(fit="ratio", t_edges="0.05,0.08,0.12,0.17,0.25"),
               dict(fit="likelihood", t_edges="0.05,0.25")):
        assert reach(_A(**kw)) != "nearbeam_reach_gain_6Li", kw
    assert (reach(_A(fit="ratio", t_edges="0.05,0.08,0.12,0.17,0.25"))
            == "nearbeam_reach_gain_6Li_tedges")
    assert (reach(_A(fit="likelihood", t_edges="0.05,0.25"))
            == "nearbeam_reach_gain_6Li_likelihood_tedges")


def test_the_published_coherent_t_window_is_the_seven_bin_one():
    """The reconstructed |t| window of the coherent intact-6Li cos 2beta
    channel became the seven bins 0.017-0.25 GeV^2 on 2026-08-28 (plans/08
    8.4): the three bins added below 0.05 carry most of the tagged sample
    and nearly halve the combined one-year delta(a_e), 0.00207 -> 0.00121
    at 5 x 40.8.  Both coherent scripts must default to the SAME list, the
    default must stay expressed as `--t-edges` unset (the sentinel the
    published stems key on), and the run-13 four-bin window must stay
    reproducible behind the flag -- under a stem of its own, so that it
    cannot overwrite either published PNG."""
    import importlib
    import pathlib as _pl
    import sys as _sys
    from polligen import recopseudo as _rp
    _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]
                            / "scripts"))

    class _A:
        def __init__(self, **kw):
            self.__dict__.update(kw)

    assert _rp.T_EDGES_PUBLISHED == (0.017, 0.028, 0.039, 0.05, 0.08,
                                     0.12, 0.17, 0.25)
    assert _rp.T_EDGES_LEGACY == (0.05, 0.08, 0.12, 0.17, 0.25)
    # the lowest published edge clears the tagging-optics aperture floor
    # |t|_min = (A p_u env_x)^2 = 0.0064 / 0.0098 / 0.0094 GeV^2, and the
    # window that would reach down to it (0.006) is the one plans/08 rules
    # out on empty beta cells, conditioning and envelope-split sensitivity
    assert _rp.T_EDGES_PUBLISHED[0] > 0.0098
    coh = importlib.import_module("money_cos2phi_coherent_reco")
    reach = importlib.import_module("nearbeam_reach_gain")
    assert reach.T_EDGES == _rp.T_EDGES_PUBLISHED
    legacy = ",".join("%g" % v for v in _rp.T_EDGES_LEGACY)
    for mod in (coh, reach):
        assert mod.t_edges_for(_A(t_edges=None)) == list(
            _rp.T_EDGES_PUBLISHED)
        assert mod.t_edges_for(_A(t_edges=legacy)) == list(
            _rp.T_EDGES_LEGACY)
        for bad in ("0.05", "0.08,0.05", "0.05,0.05"):
            with pytest.raises(SystemExit):
                mod.t_edges_for(_A(t_edges=bad))
    # the truth-level money plot 6 shades the same window
    truth = importlib.import_module("money_cos2phi_coherent")
    assert truth.T_WINDOW == (_rp.T_EDGES_PUBLISHED[0],
                              _rp.T_EDGES_PUBLISHED[-1])
    # ... and the legacy window keeps a stem of its own in BOTH scripts
    published = dict(config=0, optics="tagging", fit="ratio",
                     u_in_situ=False, t_edges=None)
    assert (coh.output_stem(_A(**dict(published, t_edges=legacy)))
            == "money_cos2phi_coherent_reco_6Li_c0_tagging_tedges")
    assert (reach.output_stem(_A(fit="ratio", t_edges=legacy))
            == "nearbeam_reach_gain_6Li_tedges")


def test_li7_two_samplers_agree_once_the_acceptance_definition_matches():
    """The published explanation of the 7Li two-sampler residual blamed the
    tagged model's momentum grid.  It is not that.  `tagging_acceptance.py`
    (Report 3 Table 6) reports 1 - lost, i.e. ANY far-forward system, while
    `tagged_polarimetry_7li.py` masks on the Roman Pots alone, and at
    5 x 40.8 the B0 carries 1.1% of the 7Li alpha -- which is the entire
    0.6-point gap at that configuration and nothing at the other two.  Like
    for like, and inside the k <= 1.2 GeV/c on which the tagged model's grid
    ends, the two samplers agree to a few tenths of a point."""
    from polli_fastsim import farforward as ff, spectator as sp

    # the tagged generator's Roman-Pot tag, per configuration (the numbers
    # tagged_polarimetry_7li.py prints as acc(RP)).  They moved by +0.0006
    # / +0.0002 / +0.0002 on 2026-08-28 when TRITON became per-nucleon like
    # every other Ion slot (plans/08 D7): the accepted sample is weighted by
    # a cross section that carries g1 of the struck triton, whose neutron
    # term gained its second neutron.
    published = (0.9620, 0.9678, 0.9728)
    b0_expected = (0.011, 0.000, 0.000)
    for i, cfg in enumerate(beams.default_configs("7Li")):
        k = sp.spectator_lab_kinematics(sp.LI7_ALPHA_TAG,
                                        cfg.ion_momentum_per_nucleon,
                                        200_000, beta=0.30,
                                        rng=np.random.default_rng(7))
        o = ff.yr_optics(cfg, "high-acceptance")
        r = ff.route_charged(k["R"], k["theta"], k["pT"], o, phi=k["phi"],
                             pot_config=ff.yr_config_key(cfg))
        rp = (r == 1) | (r == 4)
        assert float(np.mean(r == 3)) == pytest.approx(b0_expected[i],
                                                       abs=0.002), i
        # the two definitions differ by the B0 and, since 2026-08-28, by
        # the over-rigid RP-inner branch (route 6): the 7Li alpha at
        # R = 0.856 cannot reach it, but the high-k tail of the same
        # distribution can, at 3e-3 (plans/09 B1).  Up to a 3e-5
        # off-momentum sliver (route 2) that neither mask nor argument
        # turns on.
        assert (float(np.mean(r != 0)) - float(np.mean(rp))
                == pytest.approx(float(np.mean(r == 3))
                                 + float(np.mean(r == 6)), abs=1e-4)), i
        # like for like, inside the tagged model's k grid
        inside = k["k"] <= 1.2
        assert float(np.mean(inside)) > 0.99
        assert float(np.mean(rp[inside])) == pytest.approx(published[i],
                                                           abs=0.003), i


# --- plans/08 D9 and the plans/05 §5.4 deuteron-limit gate ------------------


def test_li6_b1_rank2_transfer_constant_is_pinned_to_the_model(li6):
    """`polarized.LI6_B1_RANK2_TRANSFER` is not a free number: it is the
    tagged model's own rank-2 dilution for the alpha-tagged embedded
    deuteron, a quadrature which the closed form 1 - (9/10) P_D at
    P_D = P_D_LI6 approximates to better than 1e-4 -- both tolerances
    below are that, and neither is an assertion of equality.  The
    0.87 it replaces on the b1 money plot is the VECTOR dilution
    1 - (3/2) P_D -- the wrong rank for a tensor structure function, which
    is the whole of plans/08 D9.  The test lives in evgen because it needs
    both packages; the constant lives in fastsim, which imports nothing
    from the generator."""
    from polli_fastsim import polarized

    tzz = li6.tensor_dilution()
    assert polarized.LI6_B1_RANK2_TRANSFER == pytest.approx(tzz, abs=1e-4)
    assert tzz == pytest.approx(1.0 - 0.9 * tagged.P_D_LI6, abs=1e-4)
    assert polarized.LI6_B1_LEGACY_TRANSFER == pytest.approx(
        li6.vector_dilution(), abs=2e-3)
    # what the money plot multiplies the deuteron b1 by, signal and error
    # now on the same per-nucleon footing as delta_models' dilution = 1/3
    assert polarized.b1_li6_from_deuteron(1.0) == pytest.approx(
        0.921947 / 3.0, abs=1e-6)
    assert polarized.b1_li6_from_deuteron(
        1.0, polarized.LI6_B1_LEGACY_TRANSFER, 1.0) == pytest.approx(0.87)


def test_cosyn_weiss_tensor_gate():
    """plans/05 §5.4, deuteron limit of tagged mode, quantitatively.

    Cosyn-Weiss II (arXiv:2603.23700) page 35 gives the closed form its
    FIG. 13 only illustrates.  Their Eq. (6.12) is a ratio of quadratic
    forms in the S- and D-wave radials times an angular factor
    (1 - 3 cos^2 theta_k) = -2 P2(cos theta_k); the ratio "takes values in
    [-2, 1]", the quadratic form peaks at +1 where f2/f0 = sqrt(2)
    (Eq. 6.13), which for AV18 is k = 0.30 GeV, and the angular factor runs
    from +1 at theta_k = pi/2 to -2 at theta_k = 0, pi.  TABLE II lists the
    extremal settings (A_T|| = -2 at k = 0.3 GeV, theta_k = 0; +1 at
    k = 0.3 GeV, theta_k = pi/2).

    Our `azz_tensor_curve` computes A_zz^wf = (n_+1 + n_-1 - 2 n_0) /
    (n_+1 + n_-1 + n_0) in c.m. variables, which maps onto theirs as
    A_T|| = -2 A_zz^wf.  This gate is the quantitative version of the
    plans/05 row that used to read "met qualitatively"; no digitization is
    needed, because the paper states the analytic result."""
    m = tagged.TaggedModel(tagged.deuteron_channel())
    p2 = 0.5 * (3.0 * m.c ** 2 - 1.0)

    # (a) the P2(cos theta_k) angular factor, exactly: A_zz^wf / P2 is
    #     independent of the angle bin at fixed k.  Cells within 1e-3 of
    #     the P2 zero at cos theta_k = 1/sqrt(3) are excluded, where the
    #     ratio is unbounded and says nothing; the report and the manual
    #     state the exclusion rather than claiming the full range.
    ik = int(np.argmin(np.abs(m.k - 0.30)))
    assert m.k[ik] == pytest.approx(0.3012, abs=5e-4)
    ratios = np.array([tagged.azz_tensor_curve(m, ic)[ik] / p2[ic]
                       for ic in range(m.c.size) if abs(p2[ic]) > 1e-3])
    assert ratios.max() - ratios.min() < 1e-5
    assert ratios.mean() == pytest.approx(0.99940, abs=1e-4)

    # (b) the k-envelope: the quadratic form reaches its maximum 1 at
    #     f2/f0 = sqrt(2), which CW put at k = 0.30 GeV for AV18
    ic0 = int(np.argmax(np.abs(m.c)))               # nearest cell to theta=0
    envelope = tagged.azz_tensor_curve(m, ic0) / p2[ic0]
    j = int(np.argmax(envelope))
    assert envelope[j] == pytest.approx(1.0, abs=1e-3)
    assert m.k[j] == pytest.approx(0.3098, abs=5e-4)
    assert m.k[j] == pytest.approx(0.30, abs=0.02)  # against CW's 0.30 GeV

    # (c) A_T|| = -2 A_zz^wf against CW TABLE II's +1 and -2.  The grid's
    #     outermost cos theta_k cell is 0.9896, not 1, so the exact extremes
    #     are recovered through the P2 factorization pinned in (a).
    ic90 = int(np.argmin(np.abs(m.c)))
    a_par_90 = -2.0 * tagged.azz_tensor_curve(m, ic90)[j]
    a_par_0 = -2.0 * tagged.azz_tensor_curve(m, ic0)[j]
    assert a_par_90 == pytest.approx(0.9997, abs=1e-3)     # CW: +1
    assert a_par_0 == pytest.approx(-1.9378, abs=2e-3)     # cell centre
    assert -2.0 * envelope[j] * 1.0 == pytest.approx(-2.0, abs=3e-3)
    assert -2.0 * envelope[j] * (-0.5) == pytest.approx(1.0, abs=3e-3)
    # and the whole curve stays inside CW's stated range [-2, 1]
    for ic in range(m.c.size):
        a = -2.0 * tagged.azz_tensor_curve(m, ic)
        assert np.nanmin(a) > -2.001 and np.nanmax(a) < 1.001
