"""g1A normalisation: Z and N carried, 7Li held bit-for-bit (plans/08 D7),
and the 6Li slot on the cluster picture since 2026-08-29 (plans/04 #6).

`ToyG1.g1_nucleus` used to read the Ion's effective polarizations as
WHOLE-NUCLEUS sums while every caller divided by A, so 6Li's per-nucleon
1/3 was diluted a second time and 3He's proton term lost its factor 2.
Since 2026-08-28 the builder mirrors `NuclearF2.f2a` exactly, the slots
are per-nucleon throughout, and `beams.LI7` holds the verified VMC sums
divided by Z and N so that nothing published moves.
"""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, polarized
from polli_fastsim.structure import NuclearF2, ToyF2

X = np.array([0.005, 0.02, 0.08, 0.2, 0.35, 0.5, 0.7])
Q2 = np.full_like(X, 4.0)


def test_li7_whole_nucleus_sums_reproduced_bit_for_bit():
    """The published 7Li path is 0.866 g1p - 0.037 g1n and must not move."""
    g1 = polarized.ToyG1()
    got = g1.g1_nucleus(beams.LI7, X, Q2)
    want = 0.866 * g1.g1p(X, Q2) - 0.037 * g1.g1n(X, Q2)
    assert np.allclose(got, want, rtol=0, atol=1e-15)
    # and the slots really are the sums divided by the nucleon counts
    assert beams.LI7.Z * beams.LI7.eff_pol_p == pytest.approx(0.866, abs=1e-15)
    assert beams.LI7.N * beams.LI7.eff_pol_n == pytest.approx(-0.037, abs=1e-15)


def test_g1_nucleus_weights_by_z_and_n_exactly_as_f2a_does():
    """Structural pin: the polarized and unpolarized builders agree."""
    g1 = polarized.ToyG1()
    base = ToyF2()
    for ion in beams.IONS.values():
        f2a = NuclearF2(ion, base=base).f2a(X, Q2)
        assert np.allclose(f2a, ion.Z * base.f2p(X, Q2)
                           + ion.N * base.f2n(X, Q2), rtol=1e-14)
        assert np.allclose(g1.g1_nucleus(ion, X, Q2),
                           ion.Z * ion.eff_pol_p * g1.g1p(X, Q2)
                           + ion.N * ion.eff_pol_n * g1.g1n(X, Q2), rtol=1e-14)


def test_he3_proton_term_carries_its_factor_two():
    """D7's 'restores 3He's missing x2': Z = 2 now multiplies P_p."""
    g1 = polarized.ToyG1()
    got = g1.g1_nucleus(beams.HE3, X, Q2)
    want = 2 * (-0.028) * g1.g1p(X, Q2) + 1 * 0.86 * g1.g1n(X, Q2)
    assert np.allclose(got, want, rtol=1e-14)


def test_li6_per_nucleon_g1_ratio_to_the_deuteron():
    """Per-nucleon g1(6Li)/g1(d) = (1 - 1.5 P_D_LI6)/3 = 0.290.

    Both ions are isoscalar, so the ratio is the ratio of the slots,
    P(6Li)/P(d).  In the cluster picture adopted 2026-08-29 (plans/04 #6)
    P(6Li) = (1 - 1.5 P_D_LI6)(1 - 1.5 P_D_DEUTERON)/3 and P(d) is the
    second factor alone, so the DEUTERON'S OWN D state cancels and what
    is left is the alpha-d dilution over the three protons (neutrons) the
    per-nucleon slot convention divides by.  It was 0.3584 = (1/3)/0.93
    while the naive Cloet constant and the rounded deuteron slot were the
    default, i.e. 1.233 times as large, and 0.119 before the double
    dilution went at D7."""
    g1 = polarized.ToyG1()
    ratio = ((g1.g1_nucleus(beams.LI6, X, Q2) / beams.LI6.A)
             / (g1.g1_nucleus(beams.DEUTERON, X, Q2) / beams.DEUTERON.A))
    assert np.allclose(ratio, 0.290, atol=5e-4)
    assert ratio.mean() == pytest.approx(
        (1.0 - 1.5 * beams.P_D_LI6) / 3.0, rel=1e-12)


def test_li6_slots_are_the_cluster_constants_and_nothing_hard_coded():
    """The 6Li slot is built from the wave function the tagged sector
    uses, not from a transcribed 0.81 (author decision 2026-08-29)."""
    assert beams.P_D_LI6 == 0.0867 and beams.P_D_DEUTERON == 0.045
    assert beams.LI6_CLUSTER_POLARIZATION == pytest.approx(
        (1.0 - 1.5 * beams.P_D_LI6) * (1.0 - 1.5 * beams.P_D_DEUTERON),
        rel=0, abs=0)
    assert beams.LI6_CLUSTER_POLARIZATION == pytest.approx(0.81123, abs=5e-6)
    # whole-nucleus Z*P_p = N*P_n = the cluster polarization itself
    for slot in (beams.LI6.Z * beams.LI6.eff_pol_p,
                 beams.LI6.N * beams.LI6.eff_pol_n):
        assert slot == pytest.approx(beams.LI6_CLUSTER_POLARIZATION,
                                     rel=1e-15)
    # and the deuteron carries the SAME expression, bit for bit, which is
    # what makes the ratio above exact
    assert (beams.DEUTERON.eff_pol_p
            == 1.0 - 1.5 * beams.P_D_DEUTERON == beams.DEUTERON.eff_pol_n)
    # the tagged sector re-exports these, it does not keep its own copies
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]
                           / "evgen"))
    from polligen import tagged
    assert tagged.P_D_LI6 is beams.P_D_LI6
    assert tagged.P_D_DEUTERON is beams.P_D_DEUTERON


def test_naive_one_third_stays_reachable_as_the_pre_20260829_value():
    """The retired Cloet convention, pinned: a whole-nucleus 1.0, 1.233
    times the cluster picture, and the g1(6Li) every number published
    before 2026-08-29 was computed with."""
    import dataclasses

    assert beams.LI6_NAIVE_ONE_THIRD == 1.0 / 3.0
    naive = dataclasses.replace(beams.LI6,
                                eff_pol_p=beams.LI6_NAIVE_ONE_THIRD,
                                eff_pol_n=beams.LI6_NAIVE_ONE_THIRD)
    assert naive.Z * naive.eff_pol_p == pytest.approx(1.0, rel=1e-15)
    g1 = polarized.ToyG1()
    got = g1.g1_nucleus(beams.LI6, X, Q2)
    old = g1.g1_nucleus(naive, X, Q2)
    assert np.allclose(got, old * beams.LI6_CLUSTER_POLARIZATION, rtol=1e-14)
    assert beams.LI6_NAIVE_ONE_THIRD / (beams.LI6_CLUSTER_POLARIZATION / 3.0) \
        == pytest.approx(1.233, abs=5e-4)


def test_ion_slots_are_all_per_nucleon():
    """No slot may hold a whole-nucleus sum again: |P| <= 1 per nucleon."""
    for ion in beams.IONS.values():
        assert abs(ion.eff_pol_p) <= 1.0 + 1e-12, ion.name
        assert abs(ion.eff_pol_n) <= 1.0 + 1e-12, ion.name


def test_g2_cache_key_separates_ion_grid_npts_and_medium_ratio():
    """The g2^WW cache is what makes the default-on target-mass term free
    (`polarized.ToyG1.g2_nucleus`), so what it is keyed on is a
    correctness statement, not a performance one.

    It is keyed on the ION ITSELF -- a frozen, hashable dataclass -- and
    not on id(ion), which is unique only while the object lives: a
    transient Ion could be allocated at a recycled address and collect
    another nucleus's table.  A medium-modified call is not cached at all,
    and a hit hands back a FROZEN array, so a caller writing into a result
    cannot poison every later hit.
    """
    import dataclasses

    model = polarized.ToyG1()
    a6 = model.g2_nucleus(beams.LI6, X, Q2)
    a7 = model.g2_nucleus(beams.LI7, X, Q2)
    assert not np.allclose(a6, a7)             # the ion is in the key
    assert model.g2_nucleus(beams.LI6, X, Q2) is a6      # ... and it hits
    assert not a6.flags.writeable                        # frozen
    with pytest.raises(ValueError):
        a6[0] = 1234.0
    # a value-equal Ion at a different address is the SAME key (hashable,
    # not identity), and a genuinely different one is not
    twin = dataclasses.replace(beams.LI6)
    assert model.g2_nucleus(twin, X, Q2) is a6
    other = dataclasses.replace(beams.LI6, eff_pol_p=0.5)
    assert model.g2_nucleus(other, X, Q2) is not a6
    # npts is in the key, and a medium-modified call is not cached
    assert model.g2_nucleus(beams.LI6, X, Q2, npts=32) is not a6
    n = len(model._g2_cache)
    med = model.g2_nucleus(beams.LI6, X, Q2, medium_ratio=lambda x: 0.9)
    assert len(model._g2_cache) == n
    assert not np.allclose(med, a6)
    assert model.g2_nucleus(beams.LI6, X, Q2) is a6      # and did not evict
