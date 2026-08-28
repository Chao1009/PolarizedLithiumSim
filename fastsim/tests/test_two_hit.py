"""Two-hit topology of 6Li -> alpha + d (plans/09 B4, 2026-08-28).

`spectator.breakup_lab_kinematics` samples ONE relative momentum and boosts
both fragments from it, which is what no function in this repository could
do before: LI6_ALPHA_TAG and LI6_D_TAG sampled separately are two unrelated
events.  These tests pin the three things that follow -- the two-body rest
frame, the mass-ratio angle relation, and the millimetres that Report 4
Table 5 and plans/09 SS9.2 quote -- and the two guarantees the refactor
has to keep: the single-fragment path is unchanged bit for bit, and the new
outer angular bound defaults to where every published acceptance was
computed.
"""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams
from polli_fastsim import farforward as ff
from polli_fastsim import spectator as sp

#: median alpha-d separation at the pot plane [mm], R12 = 30.6 m, pot
#: dispersion D = 0.30 m, beta = 0.30, seed 7, in configuration order
#: (5 x 41, 10 x 100, 18 x 275).  These are the numbers Report 4 Table 5,
#: plans/09 SS9.2 and the reproduction manual carry, and
#: tools/consistency_check.py compares the documents against this table.
#: They have been wrong twice, in two different ways.  Table 5 published
#: 73.4 / 30.1 / 10.9 mm until 2026-08-28: the two lower rows were computed
#: at the pre-2026-08-27 rigidity-scaled 20.5 and 50 GeV/u and were a
#: factor 2 high, which no energy-drift check could see because the stale
#: energy survives only inside a derived millimetre.  The correction to
#: 36.7 / 15.1 / 10.9 mm the same day then dropped the DISPERSIVE
#: displacement, on the argument that two fragments within 0.7% of beam
#: rigidity share it -- they do not, because they take opposite k_z and
#: their rigidities move in opposite directions (see
#: `farforward.separation_at_pots`).  Restoring it adds 5 / 23 / 39%.
MEDIAN_SEPARATION_MM = {"5x41": 38.4, "10x100": 18.5, "18x275": 15.1}

#: median alpha-d separation with the ANGULAR lever alone [mm]: what the
#: same table said between the two corrections.  Pinned so that the size of
#: the dispersive term stays visible and cannot be dropped again in
#: silence.
ANGULAR_ONLY_SEPARATION_MM = {"5x41": 36.7, "10x100": 15.1, "18x275": 10.9}


def _breakup(config, n=200_000, seed=7):
    return sp.breakup_lab_kinematics(sp.LI6_ALPHA_TAG,
                                     config.ion_momentum_per_nucleon, n,
                                     rng=np.random.default_rng(seed))


def test_the_two_fragments_are_back_to_back_in_the_rest_frame():
    """The whole physics input: a two-body breakup at rest in the beam frame
    has p_1 + p_2 = 0.  Observably, the fragments carry equal and opposite
    transverse momentum -- so equal pT and azimuths exactly pi apart -- and
    equal |k|."""
    ev = _breakup(beams.default_configs("6Li")[2], n=20_000)
    a, d = ev["spectator"], ev["partner"]
    assert np.allclose(a["pT"], d["pT"], rtol=0, atol=0)
    assert np.allclose(a["k"], d["k"], rtol=0, atol=0)
    dphi = np.abs(np.abs(a["phi"] - d["phi"]) - np.pi)
    assert float(dphi.max()) < 1e-12
    # and the rest-frame components the sampler returns are the spectator's
    assert np.allclose(np.hypot(ev["kx"], ev["ky"]), a["pT"], atol=0)


def test_the_deuteron_opens_twice_the_alpha_angle_at_the_mass_ratio():
    """Equal transverse kicks carried by different longitudinal momenta: the
    lab angles go inversely as the fragment masses, theta_d / theta_alpha ->
    m_alpha / m_d = 1.987 as k -> 0.  The naive mass-number ratio says 2;
    the 0.6% is the same physical-mass effect that puts the alpha at
    R = 0.99813 rather than 1 (plans/08 C1)."""
    ev = _breakup(beams.default_configs("6Li")[2], n=200_000)
    a, d = ev["spectator"], ev["partner"]
    at_rest = ev["k"] < 0.01
    assert at_rest.sum() > 200
    ratio = float(np.median(d["theta"][at_rest] / a["theta"][at_rest]))
    expect = sp.MASSES["alpha"] / sp.nucleus_mass(1, 2)
    assert expect == pytest.approx(1.9873, abs=1e-3)
    assert ratio == pytest.approx(expect, rel=2e-3)
    assert ratio == pytest.approx(2.0, rel=0.01)      # 2 to within the masses
    # 1.987 is the LIMIT and not an identity: away from k = 0 the ratio runs
    # over the better part of a decade, and the deuteron is the wider
    # fragment in all but ~3e-5 of breakups rather than in all of them.  The
    # documents must say "-> as k -> 0", which is why this is pinned.
    full = d["theta"] / a["theta"]
    assert float(full.min()) < 0.9 and float(full.max()) > 4.0
    wider = float(np.mean(d["theta"] > a["theta"]))
    assert 0.9999 < wider < 1.0


def test_the_single_fragment_path_is_unchanged_bit_for_bit():
    """The refactor that made the joint sampler possible put the boost in
    one private helper used by both entry points.  `spectator_lab_kinematics`
    must return exactly what it returned before -- every published spectator
    distribution depends on it -- and the joint sampler's spectator arm must
    BE that function, not a copy of it."""
    cfg = beams.default_configs("6Li")[2]
    pu = cfg.ion_momentum_per_nucleon
    one = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG, pu, 50_000,
                                      rng=np.random.default_rng(7))
    both = _breakup(cfg, n=50_000)["spectator"]
    assert set(one) == set(both)
    for key in one:
        assert np.array_equal(one[key], both[key]), key


def test_the_partner_marginal_is_the_complementary_channel():
    """The deuteron of the joint sampler and the deuteron of LI6_D_TAG are
    the same distribution: the two channels share kappa and carry
    complementary masses, so the correlation is all the joint sampler adds.
    Event by event they differ (the partner takes -k), so this is a
    distributional identity, checked on the quantiles that matter for the
    routing."""
    cfg = beams.default_configs("6Li")[2]
    pu = cfg.ion_momentum_per_nucleon
    solo = sp.spectator_lab_kinematics(sp.LI6_D_TAG, pu, 400_000,
                                       rng=np.random.default_rng(11))
    joint = _breakup(cfg, n=400_000, seed=13)["partner"]
    for key in ("theta", "pT", "R", "p_lab"):
        q_solo = np.percentile(solo[key], [16, 50, 84])
        q_joint = np.percentile(joint[key], [16, 50, 84])
        assert q_joint == pytest.approx(q_solo, rel=0.02), key
    # the rigidities are the fragment's own, not the spectator's
    assert float(np.min(joint["R"])) < 1.00452
    assert np.percentile(joint["R"], 5) == pytest.approx(
        np.percentile(solo["R"], 5), rel=2e-3)


def test_the_separation_medians_report_4_table_5_quotes():
    """Pin the three medians so Table 5 cannot go stale again.  They are
    quoted to 0.1 mm in Report 4, plans/09 and the manual; 3% is the
    sampling spread at 200k events."""
    for cfg in beams.default_configs("6Li"):
        ev = _breakup(cfg)
        sep_mm = 1e3 * ff.separation_at_pots(ev["spectator"], ev["partner"])
        med = float(np.median(sep_mm))
        assert med == pytest.approx(MEDIAN_SEPARATION_MM[ff.yr_config_key(cfg)],
                                    rel=0.03), cfg.label()
    # the ordering is the boost: the separation is an ANGLE times R12, and
    # the angle goes as 1/p_u, so the low-energy configuration is the widest
    meds = [float(np.median(ff.separation_at_pots(_breakup(c)["spectator"],
                                                  _breakup(c)["partner"])))
            for c in beams.default_configs("6Li")]
    assert meds[0] > meds[1] > meds[2]


def test_the_dispersive_term_is_the_same_size_as_the_angular_one():
    """The correction of the 2026-08-28 review.  `separation_at_pots` may
    not drop the dispersion: the two fragments carry OPPOSITE k_z, so their
    rigidities move apart rather than together (dR/dk_z = +0.268 for the
    alpha, -0.536 for the deuteron), and D (R_a - R_d) is a median 9.2 mm
    against separations of 15 to 38 mm.  Pinning both columns keeps the
    size of the term visible."""
    for cfg in beams.default_configs("6Li"):
        ev = _breakup(cfg)
        a, d = ev["spectator"], ev["partner"]
        key = ff.yr_config_key(cfg)
        full = float(np.median(1e3 * ff.separation_at_pots(a, d)))
        ang = float(np.median(1e3 * ff.separation_at_pots(a, d,
                                                          dispersion=0.0)))
        assert ang == pytest.approx(ANGULAR_ONLY_SEPARATION_MM[key], rel=0.03)
        assert full == pytest.approx(MEDIAN_SEPARATION_MM[key], rel=0.03)
        assert full > ang
        # the rigidity difference itself, and its k = 0 offset
        dr = np.asarray(a["R"]) - np.asarray(d["R"])
        assert float(np.median(np.abs(1e3 * ff.POT_DISPERSION * dr))) == \
            pytest.approx(9.2, rel=0.1)
        slope = np.polyfit(ev["kz"], a["R"], 1)[0]
        assert slope == pytest.approx(0.268, rel=0.02)
        assert np.polyfit(ev["kz"], d["R"], 1)[0] == pytest.approx(-0.536,
                                                                   rel=0.02)
    # and the fragment order does not matter
    assert np.allclose(ff.separation_at_pots(a, d),
                       ff.separation_at_pots(d, a))


def test_a_recorded_pair_merges_only_through_the_dispersion():
    """The merge probability Report 4 once bounded at 0.006-0.26%, and
    briefly retired as identically zero, is neither: it is a few 1e-4 of
    RECORDED pairs at the tagging optics and zero everywhere else.

    Both halves matter.  In the angular lever alone the fragments are back
    to back with theta_d ~ 2 theta_alpha, so a pair whose alpha cleared the
    envelope stays about 3 R12 min(env_x, env_y) apart -- a scale and not a
    bound, since 1.987 is only the k -> 0 ratio, and the measured minima run
    3-5% under it -- and no recorded pair merges.  The dispersive term is
    free to cancel that, and does."""
    for cfg in beams.default_configs("6Li"):
        ev = _breakup(cfg)
        a, d = ev["spectator"], ev["partner"]
        sep_mm = 1e3 * ff.separation_at_pots(a, d)
        ang_mm = 1e3 * ff.separation_at_pots(a, d, dispersion=0.0)
        assert float(sep_mm.min()) < 0.5          # they do exist ...
        for optics in (ff.yr_optics(cfg), ff.tagging_optics(cfg),
                       ff.HIGH_ACCEPTANCE):
            ra = ff.route_charged(a["R"], a["theta"], a["pT"], optics,
                                  phi=a["phi"])
            rd = ff.route_charged(d["R"], d["theta"], d["pT"], optics,
                                  phi=d["phi"])
            rec = ((ra == 1) | (ra == 4)) & ((rd == 1) | (rd == 4))
            if rec.sum() < 50:
                continue
            scale = 3e3 * ff.POT_R12 * min(optics.envelope)
            assert 0.90 * scale < float(ang_mm[rec].min()) < scale
            assert float(ang_mm[rec].min()) > 20.0 * 0.5   # > 20 pixels
            # with the dispersion the same pairs can close up, but rarely
            assert float(np.mean(sep_mm[rec] < 0.5)) < 1e-3


def test_the_partner_veto_given_an_alpha_that_fakes_a_coherent_tag():
    """The headline of B4.  At the tagging optics 26-34% of breakups put the
    alpha in the near-beam tail, where it is indistinguishable from an
    intact 6Li -- and 84% of those events also put the deuteron on a pot, so
    hit multiplicity alone rejects them.  At the Yellow Report optics the
    fake rate is 10^-3 and the partner is rarely there, which is the same
    conditional structure as everything else in this study: the optics that
    creates the background also supplies the veto."""
    for cfg in beams.default_configs("6Li"):
        ev = _breakup(cfg, n=200_000)
        a, d = ev["spectator"], ev["partner"]
        tag = ff.tagging_optics(cfg)
        fake = ff.route_charged(a["R"], a["theta"], a["pT"], tag,
                                phi=a["phi"]) == 4
        rd = ff.route_charged(d["R"], d["theta"], d["pT"], tag, phi=d["phi"])
        assert 0.20 < float(fake.mean()) < 0.40
        assert float(((rd == 1) | (rd == 4))[fake].mean()) == pytest.approx(
            0.84, abs=0.02)
        # the Yellow Report column of the same table is a small-sample
        # number and is quoted from 1.2e7 breakups, not from this one: at
        # 2e5 the alpha fakes a tag in ~40 events at 10 x 100.  Assert only
        # what that supports.
        yr_fake = ff.route_charged(a["R"], a["theta"], a["pT"],
                                   ff.yr_optics(cfg), phi=a["phi"]) == 4
        assert float(yr_fake.mean()) < 5e-3
        # the deuteron is the commoner single hit, at every optics
        yr = ff.yr_optics(cfg)
        seen_a = ff.route_charged(a["R"], a["theta"], a["pT"], yr,
                                  phi=a["phi"])
        seen_d = ff.route_charged(d["R"], d["theta"], d["pT"], yr,
                                  phi=d["phi"])
        ha = (seen_a == 1) | (seen_a == 4)
        hd = (seen_d == 1) | (seen_d == 4)
        assert float(np.mean(~ha & hd)) > 3.0 * float(np.mean(ha & ~hd))
        assert float(np.mean(ha & hd)) < 1e-3


def test_the_outer_bound_defaults_to_where_every_published_number_is():
    """`theta_outer` is new and the veto depends on it more than on anything
    else, so its default must not move a single published acceptance."""
    assert ff.THETA_RP_OUTER == ff.THETA_RP_MAX
    cfg = beams.default_configs("6Li")[2]
    k = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG,
                                    cfg.ion_momentum_per_nucleon, 100_000,
                                    rng=np.random.default_rng(2))
    for optics in (ff.HIGH_ACCEPTANCE, ff.yr_optics(cfg),
                   ff.tagging_optics(cfg)):
        base = ff.acceptance_summary(k["R"], k["theta"], k["pT"], optics,
                                     phi=k["phi"])
        for theta_outer in (None, ff.THETA_RP_MAX):
            assert ff.acceptance_summary(k["R"], k["theta"], k["pT"], optics,
                                         phi=k["phi"],
                                         theta_outer=theta_outer) == base
    # and it bites when it is asked to: a 16 mm module is 0.5 mrad
    tight = ff.acceptance_summary(k["R"], k["theta"], k["pT"],
                                  ff.tagging_optics(cfg), phi=k["phi"],
                                  theta_outer=0.5e-3)
    wide = ff.acceptance_summary(k["R"], k["theta"], k["pT"],
                                 ff.tagging_optics(cfg), phi=k["phi"])
    assert tight["lost"] > wide["lost"]
