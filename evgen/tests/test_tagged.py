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
    """n_M against Cosyn-Weiss Eq. (3.22b) built in CARTESIAN form.

    Genuinely independent of the module: no partial-wave sum, no
    Clebsch-Gordan coefficients, no spherical harmonics and no i^L.  CW
    write the momentum-space wave function of a spin-1 S+D pair as a
    Cartesian tensor contracted with polarization vectors,

      M^ij = (f0 - f2/sqrt2) delta^ij + 3 (f2/sqrt2) khat^i khat^j,
      n_Lambda = sum_{m_s} |eps^i_Lambda M^ij eps^j*_{m_s}|^2,

    and with f0, f2 the channel's own normalized radial tables that is
    this module's n_M times the 4pi of Y_00^2 -- exactly, to double
    precision -- but ONLY with the i^L phase in place: CW's Cartesian +f2
    is the spherical -psi_2 (S_12 = -2 T on the triplet), so a phase-less
    sum reproduces this construction at f2 -> -f2 and gets the S-D cross
    term backwards.  Until 2026-09-15 this test retyped `_amp2_table`'s
    own CG sum, which tracked whatever that function did, phase and all;
    it passed throughout the period the phase was missing.
    """
    model = li6
    ch = model.channel
    kappa = ch.base.kappa
    k = np.array([0.05, 0.15, 0.30, 0.45])
    c = np.array([-0.7, 0.1, 0.9])
    # test-local normalized radials f0, f2 (the same quadrature the model
    # uses, recomputed here so the test owns its own numbers)
    kg = model.k
    rad = {}
    for w in ch.waves:
        psi = w.radial(kg, kappa)
        norm = np.sqrt(_trapezoid(psi**2 * kg**2, kg))
        rad[w.l] = np.sqrt(w.prob) * psi / norm

    s2 = np.sqrt(2.0)
    eps = {1: -np.array([1.0, 1.0j, 0.0]) / s2,
           0: np.array([0.0, 0.0, 1.0], dtype=complex),
           -1: np.array([1.0, -1.0j, 0.0]) / s2}

    def n_cw(f0, f2, cth, lam, phi=0.0):
        sth = np.sqrt(max(0.0, 1.0 - cth * cth))
        kh = np.array([sth * np.cos(phi), sth * np.sin(phi), cth])
        mat = (f0 - f2 / s2) * np.eye(3) + 3.0 * (f2 / s2) * np.outer(kh, kh)
        return float(sum(abs(eps[lam] @ mat @ np.conj(eps[ms])) ** 2
                         for ms in (1, 0, -1)).real)

    for M in (1.0, 0.0, -1.0):
        for kk in k:
            for cc in c:
                # the module's lookup point: the nearest cell centre
                ik = np.clip(np.searchsorted(model.k, kk) - 1, 0,
                             model.k.size - 2)
                ic = np.clip(np.searchsorted(model.c, cc) - 1, 0,
                             model.c.size - 2)
                dens = n_cw(rad[0][ik], rad[2][ik], model.c[ic],
                            int(M)) / (4.0 * np.pi)
                got = model.n_of_kc(M, np.array([kk]), np.array([cc]))[0]
                assert got == pytest.approx(dens, rel=1e-9)
                # the phase-less sum is the SAME construction at f2 -> -f2,
                # and it is a different number: this is the defect the
                # retyped-CG-sum version of this test could not see
                wrong = n_cw(rad[0][ik], -rad[2][ik], model.c[ic],
                             int(M)) / (4.0 * np.pi)
                if rad[2][ik] > 1e-6:
                    assert abs(wrong - dens) > 1e-3 * abs(dens)


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
    fill -- the PROTON half of the plans/05 SS5.4 forward-limit gate.

    The gate's 0.866 is the WHOLE-NUCLEUS 7Li VMC sum (JLab PR12-14-001
    Eq. 29, from Wiringa PRC 89 (2014) 024305 Table I's 0.868), and
    TRITON.eff_pol_p is PER NUCLEON.  The two footings coincide here and
    only here, because the triton has Z = 1; the neutron half, where
    N = 2 makes them differ by a factor 2, is the test below."""
    p32 = li7.population_integrated(1.5)
    np.testing.assert_allclose(p32, [1.0, 0.0], atol=1e-12)
    p12 = li7.population_integrated(0.5)
    np.testing.assert_allclose(p12, [2.0 / 3.0, 1.0 / 3.0], atol=3e-4)
    p_t = (p12 * np.array([1.0, -1.0])).sum()
    assert p_t == pytest.approx(1.0 / 3.0, abs=6e-4)
    p_p_7li = 1.0 * tagged.TRITON.eff_pol_p  # stretched fill, Z_t = 1
    assert abs(p_p_7li - 0.866) < 0.02  # within the D-state band


#: Whole-nucleus 7Li neutron polarization of the ab initio calculation the
#: plans/05 SS5.4 gate quotes: JLab PR12-14-001 Eq. (29) rounds Wiringa
#: PRC 89 (2014) 024305 Table I (1.981 up - 2.019 down = -0.038) to -0.037.
VMC_P_N_7LI_WHOLE = -0.037
VMC_P_N_7LI_WHOLE_TABLE_I = -0.038


def test_li7_neutron_forward_limit_both_footings(li7, capsys):
    """The NEUTRON half of the gate, on both footings, with the gap to the
    ab initio value printed.

    The model's 7Li neutron spin lives entirely on the struck triton: the
    alpha spectator is spin-0 and contributes exactly zero, so with the
    triton fully polarized in the stretched fill (P_t(M = 3/2) = 1, the
    test above) the whole-nucleus sum is P_t x N_t x TRITON.eff_pol_n
    = 1 x 2 x (-0.028) = -0.056, per-nucleon -0.014 over 7Li's N = 4 and
    -0.028 over the triton's own two neutrons (which is TRITON's slot
    back again, by construction).

    -0.056 against -0.037 is a KNOWN MODEL DIFFERENCE, not a failure and
    not a band: the alpha + t decomposition puts all the neutron spin on
    two neutrons where VMC spreads it over four correlated ones, and this
    channel is a lone L = 1 wave with no D-state admixture to widen a
    tolerance around.  The tolerance below is therefore a PIN on the
    model's own value, not an agreement criterion."""
    p32 = li7.population_integrated(1.5)
    p_t = float((p32 * np.array([1.0, -1.0])).sum())
    n_t = tagged.TRITON.A - tagged.TRITON.Z          # 2 neutrons
    whole = p_t * n_t * tagged.TRITON.eff_pol_n      # whole-nucleus 7Li
    per_nucleon_a = whole / 4.0                      # over 7Li's N = 4
    per_nucleon_t = whole / n_t                      # over the triton's 2

    # the model's own values, pinned
    assert whole == pytest.approx(-0.056, abs=1e-6)
    assert per_nucleon_a == pytest.approx(-0.014, abs=1e-6)
    assert per_nucleon_t == pytest.approx(tagged.TRITON.eff_pol_n, abs=1e-12)

    # beams.LI7 holds the ab initio value on the SAME per-nucleon footing
    # as TRITON, divided by N = 4, so N * P_n returns the whole-nucleus one
    assert 4.0 * beams.LI7.eff_pol_n == pytest.approx(VMC_P_N_7LI_WHOLE,
                                                      abs=1e-12)

    gap = whole - VMC_P_N_7LI_WHOLE
    gap_table_i = whole - VMC_P_N_7LI_WHOLE_TABLE_I
    with capsys.disabled():
        print("\n7Li forward limit, neutron half (plans/05 SS5.4):"
              "\n  model, whole-nucleus   P_n = %+.6f"
              "\n  model, per nucleon/N=4 P_n = %+.6f"
              "\n  model, per triton n    P_n = %+.6f"
              "\n  VMC, whole-nucleus     P_n = %+.6f (JLab PR12-14-001)"
              "\n                             = %+.6f (Wiringa Table I)"
              "\n  GAP (known model difference, no D-state band exists):"
              " %+.6f (x%.4f) / %+.6f (x%.4f)"
              % (whole, per_nucleon_a, per_nucleon_t,
                 VMC_P_N_7LI_WHOLE, VMC_P_N_7LI_WHOLE_TABLE_I,
                 gap, whole / VMC_P_N_7LI_WHOLE,
                 gap_table_i, whole / VMC_P_N_7LI_WHOLE_TABLE_I))

    # the gap is recorded, with a tolerance measured on it (2026-09-15:
    # -0.019000, a factor 1.5135); it moves only if a channel constant or
    # the cluster decomposition moves, and then this row must be restated
    assert gap == pytest.approx(-0.019, abs=1e-6)
    assert whole / VMC_P_N_7LI_WHOLE == pytest.approx(1.5135, abs=1e-3)


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
    """S/D interference on the toy deuteron: A_zz^wf(theta_k = 90 deg) is
    POSITIVE where f2/f0 > 0, vanishes at k -> 0 (D-wave threshold) and
    grows monotonically to O(1) through the 0.25-0.5 GeV/c region.

    The sign is the physics; the monotonicity is a property of THIS toy.
    At theta_k = 90 degrees A_zz^wf is CW's quadratic form Q(f2/f0) itself
    (Eq. 6.12 at 1 - 3cos^2 = 1), which is monotone on 0 <= f2/f0 <= sqrt2
    and reaches +1 only at sqrt2 -- a ratio this node-free Hulthen pair
    never attains (max 1.287, at the grid edge) -- so the curve is monotone
    in k here only because this pair's f2/f0 is.  It is not a general
    statement: the toy 6Li pair turns over past its sqrt2 crossing at
    0.744 GeV/c.  The familiar PEAK near 300 MeV/c is not a property of
    this toy either: it is an AV18 feature, the S-wave node that sends
    f2/f0 through sqrt2 at k = 0.2988 GeV/c and puts the 90-degree
    maximum at 0.297, and it is pinned where it belongs, on the real
    wave function, in test_cosyn_weiss_table_ii_on_av18.  Before
    2026-09-15 this test pinned a peak at 0.2-0.45 GeV/c which the
    phase-less code produced by hitting Eq. (6.14)'s MINIMUM at
    f2/f0 = -1/sqrt2 there -- the right k for the wrong reason, and with
    |A_zz| in place of A_zz it could not see the sign.
    """
    ic = np.argmin(np.abs(deut.c))  # cos theta ~ 0
    n = {m: deut.n_of_kc(m)[:, ic] for m in (1.0, 0.0, -1.0)}
    a = (n[1.0] + n[-1.0] - 2 * n[0.0]) / (n[1.0] + n[-1.0] + n[0.0])
    k = deut.k
    assert a.min() >= 0.0                          # sign: f2/f0 >= 0 here
    assert abs(a[k < 0.02].max()) < 0.02           # threshold suppression
    assert np.diff(a).min() > 0.0                  # monotone in k
    assert np.abs(a[(k > 0.25) & (k < 0.5)]).max() > 0.45
    assert a.max() < 1.0                           # CW Eq. (6.13) ceiling
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
    # signs measured on the corrected wave function (2026-09-15): the
    # 90 degree curve is +0.897 here, the longitudinal Yellow Report
    # acceptance folds it to -0.953 and the transverse near-beam tail of
    # the tagging optics to +0.158.  Before the i^L phase was restored
    # these three read -0.482, +0.522 and -0.087; the STRUCTURE the test
    # exists for -- opposite sign to the 90 degree curve at the Yellow
    # Report optics, same sign and much smaller at the tagging optics --
    # is what survived, and it is what is asserted.
    assert ref[j] > +0.4                      # the 90 degree curve
    assert a_yr < -0.4                        # longitudinal acceptance
    assert 0.0 < a_tag < +0.4                 # transverse near-beam tail
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

    # `--cluster-wave vmc` redraws BOTH panels off a different wave
    # function, so it takes the `_vmc` suffix and can never land on the
    # published stem (2026-09-16, plans/04 #29)
    mta = importlib.import_module("money_tagged_azz").output_stem
    mb = "money_tagged_azz_6Li"
    assert mta(mb, "10x100", 1, "menu", cluster_wave="hulthen") == mb
    assert (mta(mb, "10x100", 1, "menu", cluster_wave="vmc")
            == mb + "_10x100_menu_vmc")
    assert (mta(mb, "10x100", 1, "menu", beta_band=True, cluster_wave="vmc")
            == mb + "_10x100_menu_betaband_vmc")
    assert (mta(mb, "5x41", 0, "tagging", cluster_wave="vmc")
            == mb + "_5x41_tagging_vmc")

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
        0.921949 / 3.0, abs=1e-6)
    assert polarized.b1_li6_from_deuteron(
        1.0, polarized.LI6_B1_LEGACY_TRANSFER, 1.0) == pytest.approx(0.87)


def test_cosyn_weiss_tensor_gate():
    """plans/05 SS5.4, deuteron limit of tagged mode, as an IDENTITY.

    Cosyn-Weiss II (arXiv:2603.23700) page 35 gives the closed form its
    FIG. 13 only illustrates.  Their Eq. (6.12),

      A_T|| = (2 f0 + f2/sqrt2)(f2/sqrt2)/(f0^2 + f2^2) x (1 - 3cos^2 th_k),

    is a statement about the spin algebra alone: it holds for ANY pair of
    radial functions, so this module's A_zz^wf must equal it to double
    precision on the whole (k, cos theta_k) grid, with f0, f2 the
    channel's own normalized S- and D-wave tables and the mapping
    A_T|| = +1 x A_zz^wf.  Eqs. (6.13)/(6.14) bound the quadratic form by
    +1 (at f2/f0 = +sqrt2) and -1/2 (at -1/sqrt2), so with the angular
    factor in [-2, +1] the whole curve lies in CW's stated [-2, 1].

    What the version of this gate retired on 2026-09-15 got wrong, twice.
    (i) It mapped A_T|| = -2 A_zz^wf; the -2 is the value of CW's angular
    factor at theta_k = 0, which A_zz^wf already carries, so the map
    double-counted it (the correct coefficient is +1, exactly, at every
    angle).  (ii) It read the toy deuteron's |f2/f0| = 1/sqrt2 at
    k = 0.3098 GeV/c -- Eq. (6.14)'s MINIMUM -- as Eq. (6.13)'s maximum at
    +sqrt2, a ratio this Hulthen pair never reaches (its f2/f0 tops out at
    1.287 on the grid), and so certified the phase-less code's M = 0 axial
    node as CW's M = +-1 node.  CW's TABLE II is an AV18 result and is
    pinned on AV18, in test_cosyn_weiss_table_ii_on_av18.
    """
    s2 = np.sqrt(2.0)
    for ch in (tagged.deuteron_channel(), tagged.li6_alpha_channel()):
        m = tagged.TaggedModel(ch)
        f0 = m._rad[0][:, None]
        f2 = m._rad[2][:, None]
        cw = ((2.0 * f0 + f2 / s2) * (f2 / s2) / (f0 ** 2 + f2 ** 2)
              * (1.0 - 3.0 * m.c[None, :] ** 2))
        got = np.stack([tagged.azz_tensor_curve(m, ic)
                        for ic in range(m.c.size)], axis=1)

        # (a) the identity itself, everywhere on the grid.  Measured
        #     residual 8.9e-16 (deuteron) and 1.1e-15 (6Li)
        assert np.abs(got - cw).max() < 1e-12

        # (b) the mapping is A_T|| = +1 x A_zz^wf, so CW's own range is
        #     the range of A_zz^wf: [-1.9319, +0.9966] for the deuteron
        #     pair and [-1.9378, +0.9997] for the alpha-d one -- both
        #     inside [-2, 1], neither attaining it
        a_par = got                                    # the +1 mapping
        assert a_par.min() > -2.0 and a_par.max() < 1.0
        assert a_par.min() < -1.9 and a_par.max() > 0.99

        # (c) the (1 - 3 cos^2 theta_k) factorization: A_zz^wf / P2 is
        #     independent of the angle bin at fixed k.  Cells within 1e-3
        #     of the P2 zero at cos theta_k = 1/sqrt(3) are excluded,
        #     where the ratio is unbounded and says nothing
        p2 = 0.5 * (3.0 * m.c ** 2 - 1.0)
        keep = np.abs(p2) > 1e-3
        ratios = got[:, keep] / p2[None, keep]
        spread = ratios.max(axis=1) - ratios.min(axis=1)
        assert spread.max() < 1e-10                    # measured 6.0e-14
        # and the factorization is -2 Q(k), Q the quadratic form
        q = (2.0 * m._rad[0] + m._rad[2] / s2) * (m._rad[2] / s2) \
            / (m._rad[0] ** 2 + m._rad[2] ** 2)
        np.testing.assert_allclose(ratios.mean(axis=1), -2.0 * q,
                                   atol=1e-12)
        # Eq. (6.13)'s ceiling Q <= 1 is ATTAINED at f2/f0 = sqrt2, which
        # the toy 6Li pair crosses at k = 0.744 GeV/c: q.max = 1 - 4.5e-8
        # there, and only because no grid cell lands on the crossing.  The
        # bound is <=, not <
        assert q.min() >= 0.0 and q.max() <= 1.0 + 1e-12

    # (d) the honest statement about these toy radials: the deuteron pair
    #     never reaches CW's f2/f0 = sqrt2, so it cannot produce TABLE II;
    #     the alpha-d pair does, at k = 0.743 GeV/c, far above the tagged
    #     window.  Both are node-free Hulthen-type forms, which is why
    #     neither has AV18's f2/f0 -> sqrt2 crossing at 0.30 GeV/c
    d = tagged.TaggedModel(tagged.deuteron_channel())
    rd = d._rad[2] / d._rad[0]
    assert rd.max() == pytest.approx(1.2866, abs=1e-3)
    assert d.k[int(np.argmax(rd))] == pytest.approx(d.k[-1], abs=1e-9)
    assert rd.max() < s2
    li6m = tagged.TaggedModel(tagged.li6_alpha_channel())
    r6 = li6m._rad[2] / li6m._rad[0]
    i = int(np.flatnonzero(np.diff(np.sign(r6 - s2)))[0])
    k_cross = li6m.k[i] + ((s2 - r6[i]) * (li6m.k[i + 1] - li6m.k[i])
                           / (r6[i + 1] - r6[i]))
    assert k_cross == pytest.approx(0.743, abs=5e-3)


def test_cosyn_weiss_table_ii_on_av18():
    """CW TABLE II, on the Argonne v18 deuteron rather than on a toy.

    `tagged.av18_deuteron_channel()` is the same spin structure as the
    deuteron control with the tabulated AV18 u(k), w(k) in place of the
    Hulthen forms (`polli_fastsim/data/av18/fdeut.av18`, provenance in
    `data/SOURCES.md`).  It is the only channel here that can reproduce
    the paper's numbers, because they are AV18 numbers: their f2/f0 passes
    +sqrt2 at k = 0.30 GeV, where CW note "the polarized neutron
    distributions have a node", and -1/sqrt2 near 1 GeV.

    Cell-centre caveat: the grid's outermost cos theta_k cell is 0.9896,
    not 1, so the theta_k = 0 row reads -1.937 rather than -2: the angular
    factor is 1 - 3 cos^2 th_k = -2 P2(0.9896) = -1.9378, not the -2 it
    takes on the axis.  The exact -2 is recovered through the P2
    factorization pinned in test_cosyn_weiss_tensor_gate, and directly by
    the fine-angle grid below (-1.998 at |cos theta_k| = 0.99975).
    """
    m = tagged.TaggedModel(tagged.av18_deuteron_channel())
    s2 = np.sqrt(2.0)
    ratio = m._rad[2] / m._rad[0]
    ic0 = int(np.argmax(np.abs(m.c)))          # nearest cell to theta_k = 0
    ic90 = int(np.argmin(np.abs(m.c)))         # nearest cell to 90 degrees
    assert abs(m.c[ic0]) == pytest.approx(0.9896, abs=1e-3)

    # row 1 and 2 of TABLE II: k = 0.3 GeV, f2/f0 = sqrt2, A_T|| = -2 at
    # theta_k = 0 and +1 at theta_k = pi/2 (mapping A_T|| = +1 A_zz^wf)
    ik = int(np.argmin(np.abs(m.k - 0.30)))
    assert m.k[ik] == pytest.approx(0.3012, abs=5e-4)
    a0 = tagged.azz_tensor_curve(m, ic0)[ik]
    a90 = tagged.azz_tensor_curve(m, ic90)[ik]
    assert a0 == pytest.approx(-1.937, abs=3e-3)        # CW: -2
    assert abs(a0 - (-2.0)) < 0.07                      # the cell centre
    assert a90 == pytest.approx(+1.0, abs=3e-3)         # CW: +1

    # the f2/f0 = sqrt2 crossing that puts CW's landmark at 0.30 GeV
    i = int(np.flatnonzero(np.diff(np.sign(ratio - s2)))[0])
    k_cross = m.k[i] + ((s2 - ratio[i]) * (m.k[i + 1] - m.k[i])
                        / (ratio[i + 1] - ratio[i]))
    assert k_cross == pytest.approx(0.299, abs=5e-3)

    # the node is in n_{+-1} ALONG THE SPIN AXIS, not in n_0: that is the
    # convention-free content of TABLE II and the reason the i^L phase
    # matters.  Resolved on a fine cos theta_k grid, where the outermost
    # cell is 0.99975 and the off-axis leakage no longer masks the zero
    fine = tagged.TaggedModel(tagged.av18_deuteron_channel(), nc=4001)
    jf = int(np.argmin(np.abs(fine._rad[2] / fine._rad[0] - s2)))
    icf = int(np.argmax(np.abs(fine.c)))
    n1 = fine.n_of_kc(1.0)[:, icf]
    n0 = fine.n_of_kc(0.0)[:, icf]
    assert fine.k[jf] == pytest.approx(0.2968, abs=5e-3)
    assert n1[jf] / n0[jf] < 1e-3                       # measured 2.7e-4
    for k_off in (0.2, 0.4):                            # a node, not smallness
        j = int(np.argmin(np.abs(fine.k - k_off)))
        assert n1[j] / n0[j] > 0.15
    assert tagged.azz_tensor_curve(fine, icf)[jf] == pytest.approx(
        -2.0, abs=3e-3)                                 # -1.9984 measured

    # row 3: k = 1 GeV, theta_k = 0, f2/f0 = -1/sqrt2 -> A_T|| = +1.  CW
    # call this row "not presumed to be a realistic prediction"; it pins
    # the OTHER extremum of the quadratic form, Eq. (6.14)
    j1 = int(np.argmin(np.abs(m.k - 1.0)))
    assert m.k[j1] == pytest.approx(0.9979, abs=5e-4)
    assert tagged.azz_tensor_curve(m, ic0)[j1] == pytest.approx(
        +0.967, abs=0.01)
    i2 = int(np.flatnonzero(np.diff(np.sign(ratio + 1.0 / s2)))[-1])
    k2 = m.k[i2] + ((-1.0 / s2 - ratio[i2]) * (m.k[i2 + 1] - m.k[i2])
                    / (ratio[i2 + 1] - ratio[i2]))
    assert k2 == pytest.approx(1.03, abs=0.01)


def test_mixed_parity_channel_is_rejected():
    """The i^L phase `_amp2_table` applies is real only while every wave
    of a channel has the same L parity; mix L = 0 with L = 1 and i^1 = i
    makes the amplitude complex, which this module's real arrays cannot
    carry.  Parity forbids the mixture for a state of good parity anyway,
    so the guard is against a mis-specified channel."""
    with pytest.raises(ValueError) as err:
        tagged.TaggedChannel(spectator.LI6_ALPHA_TAG, 1.0, 1.0, 0.0, 1.0,
                             (tagged.Wave(0, 0.6), tagged.Wave(1, 0.4)),
                             beams.DEUTERON, "mixed parity")
    assert "L mod 2" in str(err.value)
    # the shipped channels are all single-parity and still build
    for ch in (tagged.li6_alpha_channel(), tagged.li7_alpha_channel(),
               tagged.deuteron_channel(), tagged.av18_deuteron_channel()):
        assert len({w.l % 2 for w in ch.waves}) == 1


# --- the ANL alpha-d VMC tables (plans/04 #15, #29) --------------------
#
# The two files under `polli_fastsim/data/vmc` are the raw served bytes of
# R. B. Wiringa et al.'s ANL VMC output (provenance in `data/SOURCES.md`).
# They are also what LiPolGen reads for its own `ClusterWaveSource::VmcAV18`
# 6Li channel, so every number below is a two-implementation check: the
# reference values are LiPolGen's, quoted from its
# `docs/open_items/vmc_reconciliation.md`, and reproduced here by an
# independent Python loader over the same bytes.

#: `docs/open_items/vmc_reconciliation.md`, "Tagged tensor asymmetry
#: A_zz^tag(k), 6Li, YR high-acceptance", the VMC AV18 column -- printed to
#: four decimals, which is the tolerance below.
LIPOLGEN_AZZ_TAG_VMC = ((0.1979, -0.5191), (0.2495, -0.2899),
                        (0.3012, -0.1993), (0.4001, -0.0822),
                        (0.4990, +0.1170))
#: the same table's Hulthen beta = 0.30 column, the pin that the default
#: path did not move when the VMC one was added.
LIPOLGEN_AZZ_TAG_HULTHEN = ((0.1979, -1.2069), (0.2495, -1.0639),
                            (0.3012, -0.9533), (0.4001, -0.7267),
                            (0.4990, -0.5851))


@pytest.fixture(scope="module")
def li6_vmc():
    return tagged.TaggedModel(tagged.li6_alpha_channel(wave="vmc"))


def test_vmc_li6_normalization_and_p_d():
    """The S/D split block reproduces the file's OWN printed norms, and
    P_D comes out where LiPolGen documents it.

    `li6_ad1.momentum` prints `4*PI*TOTINT(RHO*K**2:K)/(2*PI)**3` for the
    total block and for each of the S and D blocks; the loader takes P_D
    from the printed S/D pair rather than from a re-integration, and this
    checks that the tabulated columns carry the same content.  The 2004
    `li6.ad` overlap file, committed beside it as the sign reference, is
    also the independent cross-check on the norm: a different Hamiltonian
    (AV18+UIX vs AV18+UX) and 200k samples instead of 1M, which is why the
    two P_D's differ in the third decimal.
    """
    mo = tagged._vmc_data_text(tagged.VMC_LI6_MOMENTUM)
    norms = tagged._parse_anl_momentum_norms(mo)
    assert norms == [0.81971, 0.80362, 0.015861]        # the printed lines

    blocks = tagged._parse_anl_momentum(mo)
    assert len(blocks) == 2
    x_tot, c_tot, _e = blocks[0]
    x_sd, c_sd, e_sd = blocks[1]
    assert len(c_tot) == 1 and len(c_sd) == 2 and len(e_sd) == 2
    assert x_sd.size == 51 and x_sd[0] == 0.001 and x_sd[-1] == 5.0
    fac = 4.0 * np.pi / (2.0 * np.pi) ** 3
    # the file's quadrature is not stated; the trapezoid reproduces every
    # printed norm to 2e-5 relative (measured 7e-6, 1.6e-5, 1.9e-4)
    assert fac * _trapezoid(c_tot[0] * x_tot**2, x_tot) == pytest.approx(
        norms[0], rel=2e-5)
    assert fac * _trapezoid(c_sd[0] * x_sd**2, x_sd) == pytest.approx(
        norms[1], rel=2e-5)
    assert fac * _trapezoid(c_sd[1] * x_sd**2, x_sd) == pytest.approx(
        norms[2], rel=2e-4)

    k, psi_s, psi_d, p_d, nodes = tagged.li6_vmc_tables()
    assert p_d == pytest.approx(0.015861 / (0.80362 + 0.015861), rel=0, abs=0)
    assert p_d == pytest.approx(0.0193549, abs=1e-7)    # LiPolGen VMC_P_D_LI6
    # LiPolGen's documented spread over its three estimators of the SAME
    # quantity (momentum file, overlap k-block, overlap r-block FT)
    assert 0.0193 <= p_d <= 0.0207
    n0 = _trapezoid(k**2 * psi_s**2, k)
    n2 = _trapezoid(k**2 * psi_d**2, k)
    # the table's own trapezoid against the printed split: 0.0193516 vs
    # 0.0193549, the same 1.7e-4 relative LiPolGen measures
    assert n2 / (n0 + n2) == pytest.approx(0.0193516, abs=1e-7)
    assert n2 / (n0 + n2) == pytest.approx(p_d, rel=2e-4)

    # the 2004 overlap file, the sign reference, on its own normalization
    # S_ad = (2 pi)^-3 int k^2 A^2 dk -- LiPolGen reports 0.85463 / 0.0201
    k_ov, a0, a2 = tagged._parse_anl_overlap_k(
        tagged._vmc_data_text(tagged.VMC_LI6_OVERLAP))
    assert k_ov.size == 51 and k_ov[0] == 0.0 and k_ov[-1] == 5.0
    s0 = _trapezoid(k_ov**2 * a0**2, k_ov) / (2.0 * np.pi) ** 3
    s2 = _trapezoid(k_ov**2 * a2**2, k_ov) / (2.0 * np.pi) ** 3
    assert s0 + s2 == pytest.approx(0.85463, abs=5e-5)
    assert s2 / (s0 + s2) == pytest.approx(0.0201, abs=5e-5)

    # and the channel carries the file's P_D, not P_D_LI6
    ch = tagged.li6_alpha_channel(wave="vmc")
    assert [w.l for w in ch.waves] == [0, 2]
    assert ch.waves[1].prob == p_d
    assert ch.waves[0].prob == 1.0 - p_d
    assert nodes == pytest.approx((0.13376, 0.44394), abs=1e-5)


def test_vmc_li6_sign_regions():
    """THE measurement the tables exist for: the relative S-D phase.

    A momentum density is |psi_L|^2 and carries no phase, so the sign
    comes from `li6.ad`'s signed k-space columns.  With the global phase
    fixed by psihat_0(k -> 0) > 0 there are exactly two zero crossings
    below the noise floor at 3 fm^-1 and three sign regions, and this pins
    both crossings within 0.01 GeV/c of the alpha-d S node at 0.134 and
    the alpha-d D node at 0.444 GeV/c.  Each is confirmed independently,
    in the same 0.1 fm^-1 bin, by a minimum of the momentum file's own
    rho_L -- a density does not know the phase but it does know where the
    amplitude vanishes.
    """
    k, psi_s, psi_d, _p_d, nodes = tagged.li6_vmc_tables()
    assert psi_s[0] > 0.0                       # the global-phase convention

    # psihat_L = s(k) sqrt(rho_L) with s a STEP function, so the ratio
    # jumps across a node rather than passing through zero: the crossings
    # are the signed reference amplitude's own, which is what the loader
    # returns, and the table can only bracket them to its 0.1 fm^-1 =
    # 0.0197 GeV/c spacing.
    ratio = psi_d / psi_s
    brackets = [i for i in range(k.size - 1)
                if ratio[i] * ratio[i + 1] < 0.0 and k[i + 1] < 0.6]
    assert len(brackets) == 2
    assert len(nodes) == 2
    assert nodes[0] == pytest.approx(0.134, abs=0.01)       # 0.13376
    assert nodes[1] == pytest.approx(0.444, abs=0.01)       # 0.44394
    for node, i in zip(nodes, brackets):
        assert k[i] < node < k[i + 1]
    crossings = nodes

    # the three regions, sampled away from the crossings
    assert np.all(ratio[k < crossings[0] - 0.01] < 0.0)
    mid = (k > crossings[0] + 0.01) & (k < crossings[1] - 0.01)
    assert np.all(ratio[mid] > 0.0)
    upper = (k > crossings[1] + 0.01) & (k < 0.59)
    assert np.all(ratio[upper] < 0.0)

    # each crossing is a node of ONE wave, confirmed by a minimum of the
    # momentum file's own |psi_L| in the same 0.1 fm^-1 = 0.0197 GeV bin
    for cross, psi in ((crossings[0], psi_s), (crossings[1], psi_d)):
        win = np.flatnonzero(np.abs(k - cross) < 0.06)   # +-3 table bins
        j = win[int(np.argmin(np.abs(psi[win])))]
        assert abs(k[j] - cross) < 0.0198

    # the Hulthen pair can carry none of this: it is positive-definite
    hul = tagged.li6_alpha_channel()
    kk = np.linspace(0.01, 0.6, 200)
    kappa = hul.base.kappa
    assert np.all(hul.waves[0].radial(kk, kappa) > 0.0)
    assert np.all(hul.waves[1].radial(kk, kappa) > 0.0)


def test_vmc_li6_azz_tag_matches_lipolgen(li6_vmc, li6):
    """The acceptance-weighted tagged tensor asymmetry, against the
    independent C++ implementation reading the same bytes.

    LiPolGen's table (`docs/open_items/vmc_reconciliation.md`) is
    40000-event-independent -- it is a MODEL integral, not a sample -- and
    is printed to four decimals at the YR high-acceptance optics of
    e x ion = 10 x 99.5 GeV/u, on the same (k, cos theta_k) grid this
    module uses.  Measured agreement is 3.7e-5 at worst, i.e. their
    printed rounding; the pin is 5e-5.

    THEIR n_phi IS 32 AND THIS MODULE'S DEFAULT IS 64, and at these cells
    that makes no difference AT ALL -- not a small one.  Measured: the two
    azimuthal grids give bit-identical curves for every k < 0.718 GeV/c,
    because at the Yellow Report optics the only alpha that survive are
    the off-rigidity R < 0.95 window slice and the R cut does not depend
    on the azimuth, so eps(k, c) is exactly 0 or 1 there (asserted below).
    n_phi first matters at k = 0.718 GeV/c, where the accepted sample has
    0.08 % of its weight in total.
    """
    from polli_fastsim import farforward as ff

    cfg = beams.default_configs("6Li")[1]
    assert "99.5" in cfg.label() and cfg.label().startswith("e(10)")
    optics = ff.yr_optics(cfg, "high-acceptance")

    eps64 = tagged.acceptance_weights(li6_vmc, cfg, optics)
    eps32 = tagged.acceptance_weights(li6_vmc, cfg, optics, n_phi=32)
    low = li6_vmc.k < 0.6
    assert set(np.unique(eps64[low])) == {0.0, 1.0}
    assert np.array_equal(eps64[low], eps32[low])

    a64 = tagged.azz_tensor_curve(li6_vmc, weights=eps64)
    a32 = tagged.azz_tensor_curve(li6_vmc, weights=eps32)
    for k_ref, azz_ref in LIPOLGEN_AZZ_TAG_VMC:
        i = int(np.argmin(np.abs(li6_vmc.k - k_ref)))
        assert li6_vmc.k[i] == pytest.approx(k_ref, abs=5e-4)
        assert a64[i] == pytest.approx(azz_ref, abs=5e-5)
        assert a32[i] == a64[i]                  # exactly, see the docstring

    # the same table's Hulthen column, on the DEFAULT channel: the pin
    # that adding the `wave` switch moved nothing
    a_hul = tagged.azz_tensor_curve(
        li6, weights=tagged.acceptance_weights(li6, cfg, optics))
    for k_ref, azz_ref in LIPOLGEN_AZZ_TAG_HULTHEN:
        i = int(np.argmin(np.abs(li6.k - k_ref)))
        assert a_hul[i] == pytest.approx(azz_ref, abs=5e-5)

    # nothing at all is accepted below k = 0.189 GeV/c at this optics
    assert np.isnan(a64[li6_vmc.k < 0.189]).all()


def test_vmc_li6_tensor_dilution_and_accepted_fraction(li6_vmc, li6):
    """khat,k-integrated <3 m_S^2 - 2> of the embedded deuteron, and the
    spin-blind accepted fraction, both against LiPolGen.

    LiPolGen prints 0.9825758723 for its VMC alpha-d channel and this
    module gives 0.9825758755.  The 3.2e-9 is not quadrature noise and it
    is worth naming, because it is ONE GRID POINT: the model grid starts
    at k = 1e-4 GeV/c and the S/D momentum table's first abscissa is
    K = 0.001 fm^-1 = 1.97e-4 GeV/c, so exactly one cell falls below the
    table.  `TabulatedWave` extrapolates it flat and LiPolGen's
    `VmcRadial` returns zero there.  Flat is the better reading -- the
    file's own total block prints rho(K = 0) = 1041.5, the same value as
    rho_0(K = 0.001), so the density really is flat there -- and adopting
    LiPolGen's rule instead reproduces its ten printed digits exactly
    (0.982575872345, measured).  Nothing else moves appreciably with it:
    the five A_zz^tag cells above shift by at most 4.7e-8 absolute
    between the two rules (largest at k = 0.1979: -0.519129029686 flat
    against -0.519128982429 zeroed), three orders below the 5e-5 at which
    they are pinned against LiPolGen.

    The Hulthen default is 0.9219 -- its P_D is 4.5x the VMC one, so its
    tensor dilution is correspondingly further from 1.
    """
    from polli_fastsim import farforward as ff

    assert li6_vmc.tensor_dilution() == pytest.approx(0.9825758723, abs=1e-8)
    assert li6_vmc.vector_dilution() == pytest.approx(
        1.0 - 1.5 * li6_vmc.channel.waves[1].prob, abs=2e-5)
    assert li6.tensor_dilution() == pytest.approx(0.9219489770, abs=1e-8)
    for m in (1.0, 0.0, -1.0):
        assert li6_vmc.norm(m) == pytest.approx(1.0, abs=1e-5)

    # the spin-blind (uniform-M) accepted fraction as a MODEL integral,
    # which LiPolGen publishes to fifteen digits at n_phi = 32: the
    # Hulthen value reproduces every one of them (0.024675932148827632
    # against a printed 0.024675932148828, i.e. to 1.5e-14 relative --
    # two independent implementations of the same integral, agreeing at
    # the last double-precision digit) and the VMC one agrees to 7.3e-8
    # relative (the one grid point of the docstring, plus the two
    # libraries' own table normalizations)
    cfg = beams.default_configs("6Li")[1]
    optics = ff.yr_optics(cfg, "high-acceptance")
    for model, ref, tol in ((li6, 0.024675932148828, 1e-13),
                            (li6_vmc, 0.033810227625842, 1e-7)):
        n_tot = sum(model.n_of_kc(m) for m in (1.0, 0.0, -1.0))
        eps = tagged.acceptance_weights(model, cfg, optics, n_phi=32)
        w = n_tot * model.k[:, None] ** 2
        assert (eps * w).sum() / w.sum() == pytest.approx(ref, rel=tol,
                                                          abs=0.0)


def test_vmc_wave_switch_is_opt_in_and_guarded(li6):
    """`wave` selects a radial input and nothing else; the default path is
    bit-for-bit the analytic pair, and every way of asking for something
    the tables cannot give is an error rather than a silent no-op."""
    default = tagged.li6_alpha_channel()
    explicit = tagged.li6_alpha_channel(wave="hulthen")
    assert default == explicit                        # frozen dataclasses
    assert all(isinstance(w, tagged.Wave) for w in default.waves)
    m_def = tagged.TaggedModel(default)
    m_exp = tagged.TaggedModel(explicit)
    for m in (1.0, 0.0, -1.0):
        assert np.array_equal(m_def._amp2_table(m), m_exp._amp2_table(m))
        assert np.array_equal(m_def._amp2_table(m), li6._amp2_table(m))

    vmc = tagged.li6_alpha_channel(wave="vmc")
    assert all(isinstance(w, tagged.TabulatedWave) for w in vmc.waves)
    assert vmc.base is default.base and vmc.dis_target is default.dis_target

    # beta and p_d have no meaning for a table
    with pytest.raises(ValueError, match="do not apply"):
        tagged.li6_alpha_channel(beta=0.40, wave="vmc")
    with pytest.raises(ValueError, match="do not apply"):
        tagged.li6_alpha_channel(p_d=0.04, wave="vmc")
    with pytest.raises(ValueError, match="hulthen"):
        tagged.li6_alpha_channel(wave="av18")

    # 7Li: refused, and for a stated reason -- a lone L = 1 wave has no
    # interference and no observable relative phase, which is the one
    # thing the ANL overlaps carry that the analytic forms cannot
    with pytest.raises(ValueError) as err:
        tagged.li7_alpha_channel(wave="vmc")
    assert "no VMC table is committed" in str(err.value)
    assert len(tagged.li7_alpha_channel().waves) == 1


def test_vmc_li6_accepted_sign_region_fractions(li6_vmc, li6):
    """Which sign region the accepted alpha actually sit in -- the numbers
    the `P_D_LI6` sign paragraph states, re-measured here.

    The answer depends on the density asked, and that is the point: the
    Hulthen pair puts 21-27 % of the Yellow-Report-accepted sample ABOVE
    the alpha-d D node, where its positive-definite forms carry the wrong
    sign, while the VMC pair -- whose far tail is 8x softer -- puts 2-3 %
    there.  Both numbers are on the uniform-M mixture.
    """
    from polli_fastsim import farforward as ff

    node_lo, node_hi = tagged.li6_vmc_tables()[4]
    expect = {                      # (below, between, above) per config
        "vmc": {"high-acceptance": [(0.000, 0.971, 0.029),
                                    (0.000, 0.981, 0.019),
                                    (0.000, 0.980, 0.020)],
                "tagging": [(0.255, 0.740, 0.005),
                            (0.113, 0.881, 0.006),
                            (0.206, 0.789, 0.005)]},
        "hulthen": {"high-acceptance": [(0.000, 0.732, 0.268),
                                        (0.000, 0.787, 0.213),
                                        (0.000, 0.741, 0.259)],
                    "tagging": [(0.405, 0.566, 0.029),
                                (0.281, 0.680, 0.039),
                                (0.367, 0.600, 0.033)]},
    }
    for name, model in (("vmc", li6_vmc), ("hulthen", li6)):
        n_tot = sum(model.n_of_kc(m) for m in (1.0, 0.0, -1.0))
        for ic, cfg in enumerate(beams.default_configs("6Li")):
            for which in ("high-acceptance", "tagging"):
                optics = (ff.yr_optics(cfg, which) if which != "tagging"
                          else ff.tagging_optics(cfg))
                eps = tagged.acceptance_weights(model, cfg, optics)
                w = (eps * n_tot * model.k[:, None] ** 2).sum(axis=1)
                got = (w[model.k < node_lo].sum() / w.sum(),
                       w[(model.k >= node_lo)
                         & (model.k <= node_hi)].sum() / w.sum(),
                       w[model.k > node_hi].sum() / w.sum())
                assert got == pytest.approx(expect[name][which][ic],
                                            abs=1e-3), (name, which, ic)

    # and the k = 0.325 GeV/c headline bin of money plot 4 has the SAME
    # sign on both wave functions -- it sits inside the supported window
    i = int(np.argmin(np.abs(li6_vmc.k - 0.325)))
    ic90 = int(np.argmin(np.abs(li6_vmc.c)))
    a_vmc = tagged.azz_tensor_curve(li6_vmc, ic90)[i]
    a_hul = tagged.azz_tensor_curve(li6, ic90)[i]
    assert a_vmc == pytest.approx(0.1698, abs=1e-3)
    assert a_hul == pytest.approx(0.9241, abs=1e-3)
    assert a_vmc > 0.0 and a_hul > 0.0
