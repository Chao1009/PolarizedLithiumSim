"""Grid-backend tests (skipped automatically if parton/grids are absent)."""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

parton = pytest.importorskip("parton")


def _have(setname):
    try:
        from parton import mkPDF
        mkPDF(setname, 0)
        return True
    except Exception:
        return False


@pytest.mark.skipif(not _have("CT18NLO"), reason="CT18NLO grid not installed")
def test_parton_f2_sane():
    from polli_fastsim.structure import PartonF2
    f2 = PartonF2()
    # HERA-anchored magnitudes and n/p < 1 at valence x
    assert 0.8 < f2.f2p(1e-3, 10.0) < 1.3
    assert 0.3 < f2.f2p(0.1, 10.0) < 0.55
    assert f2.f2n(0.5, 10.0) < f2.f2p(0.5, 10.0)
    # vectorized call works on arrays
    vals = f2.f2p(np.array([1e-3, 0.1, 0.5]), np.array([10.0, 10.0, 10.0]))
    assert vals.shape == (3,) and np.all(vals > 0)


@pytest.mark.skipif(not _have("NNPDFpol11_100"),
                    reason="NNPDFpol11 grid not installed")
def test_parton_g1_sane():
    from polli_fastsim.polarized import PartonG1
    g1 = PartonG1()
    # proton g1 positive at mid x; A1p grows with x and stays below 1
    assert g1.g1p(0.1, 10.0) > 0
    f1p = g1.base.f2p(0.3, 10.0) / (2 * 0.3 * 1.18)
    a1p = float(g1.g1p(0.3, 10.0)) / f1p
    assert 0.2 < a1p < 0.9
    # neutron g1 small/negative at low-mid x (isospin swap sanity)
    assert g1.g1n(0.05, 10.0) < g1.g1p(0.05, 10.0)


def test_flavour_schemes_of_f2_and_g1_are_the_documented_pair():
    """F2 is five-flavour, g1 three-flavour, and the second is a SUBSET of
    the first with the same charges.

    Not a style check: every g1/F1 ratio in the package mixes the two, and
    the mixture is deliberate (NNPDFpol1.1 sets Delta c = Delta b = 0, so
    a three-flavour g1 is all it predicts, while the measured F2 carries
    charm).  This pins the pair so that neither can be widened or narrowed
    without the docstrings that explain the ratio moving with it.
    """
    from polli_fastsim.polarized import PartonG1
    from polli_fastsim.structure import PartonF2
    assert set(PartonF2._E2) == {1, 2, 3, 4, 5}          # d u s c b
    assert set(PartonG1._E2) == {1, 2, 3}                # d u s
    assert set(PartonG1._E2) < set(PartonF2._E2)
    for pid, e2 in PartonG1._E2.items():
        assert PartonF2._E2[pid] == e2


@pytest.mark.skipif(not _have("CT18NLO"), reason="CT18NLO grid not installed")
def test_charm_share_of_f2_is_small_where_the_polarized_emc_lives():
    """The size of the mixed scheme, at the bins money_polemc.py combines.

    Measured 2026-08-28 over the accepted (x, Q2) bins of the three ⁷Li
    configurations at 10 fb^-1/u with >= 100 events: charm + bottom carry
    7.8% of F2A event-weighted over all of them, but only 0.65% at
    x = 0.3-0.5 and 0.23% at x = 0.5-0.7.  That window is where the two
    digitized polarized-EMC camps genuinely differ -- the valence window
    `polarized.POLEMC_VALENCE_WINDOW` in which their transfer factor is
    defined -- and the share there is 0.31% at x = 0.35 down to 0.03% at
    0.65 at Q2 = 4 GeV2, far below the 0.03-0.04 separation in DR.
    `money_polemc.py` also plots low-x bins, and its unrestricted best bin
    sits at x = 0.141; the share there is 2.6% at x = 0.09 and 1.4% at
    0.14 at Q2 = 4 GeV2 (4.3% and 2.4% at Q2 = 10), i.e. comparable to the
    relative separation those bins are drawn with, which is why the report
    reads the discrimination off the valence window instead.  This checks
    both ends of that statement at a representative Q2 rather than
    re-running the whole projection.
    """
    from polli_fastsim.structure import PartonF2, _safe_xfx
    f2 = PartonF2()
    pdf = f2._pdf

    def share(x, q2):
        tot = light = 0.0
        for pid, e2 in PartonF2._E2.items():
            v = e2 * (_safe_xfx(pdf, pid, x, q2) + _safe_xfx(pdf, -pid, x, q2))
            tot += v
            if pid <= 3:
                light += v
        return 1.0 - light / tot

    assert 0.10 < share(3e-3, 100.0) < 0.30     # low x, high Q2: the big end
    assert 0.002 < share(0.4, 10.0) < 0.02      # the polarized-EMC window
    assert share(0.65, 10.0) < share(0.4, 10.0)
    # the valence window, where the two digitized camps genuinely differ:
    # negligible against the 0.03-0.04 separation in DR
    from polli_fastsim.polarized import POLEMC_VALENCE_WINDOW
    vlo, vhi = POLEMC_VALENCE_WINDOW
    assert share(vlo, 4.0) < 0.005
    assert share(vhi, 4.0) < share(vlo, 4.0) < share(0.14, 4.0)
    assert share(vlo, 10.0) < 0.01
    # the low-x bins money_polemc.py also plots: comparable to the
    # separation drawn there, so the report does not read them as the
    # discriminating ones
    assert 0.02 < share(0.09, 4.0) < 0.04
    assert 0.01 < share(0.14, 4.0) < 0.02
    assert share(0.14, 10.0) < share(0.09, 10.0)


@pytest.mark.skipif(not (_have("CT18ANLO") and _have("EPPS21nlo_CT18Anlo_Li6")
                         and _have("nNNPDF30_nlo_as_0118_A6_Z3")),
                    reason="A = 6 nuclear grids not installed")
def test_the_a6_grids_are_per_nucleon_and_give_the_documented_emc_ratio():
    """`polarized.unpolarized_emc_ratio`'s data-driven default.

    Two statements, both load-bearing for the polarized-EMC baseline.
    FIRST, the A = 6 grids hold the AVERAGE BOUND NUCLEON and not the
    bound proton, which is what lets `PartonF2.f2p` on them be F2^A/A:
    Li-6 has Z = N = 3, so the average nucleon is isoscalar and x*u must
    equal x*d, while the free CT18ANLO proton is nowhere near that.
    SECOND, the ratio it builds against the free isoscalar nucleon is the
    shallow EMC effect the projection is now based on -- 1.012 at x = 0.1,
    0.967 at 0.5 and 0.967 at 0.7 for EPPS21 at Q2 = 5 GeV2, against the
    0.939 and 0.910 of CBT's model 7Li curve at the last two.  The
    denominator is CT18ANLO since 2026-08-29 -- EPPS21's OWN proton
    baseline, so the fit cancels; on CT18NLO the same ratio read 1.010 /
    0.969 / 0.965 and its valence depletion was 4.2% shallower.
    """
    from polli_fastsim.polarized import unpolarized_emc_ratio
    from polli_fastsim.structure import (NUCLEAR_F2_SETS, NuclearF2Ratio,
                                         f2_backend)
    for key, setname in NUCLEAR_F2_SETS.items():
        pdf = f2_backend(setname)._pdf
        for x in (0.01, 0.1, 0.5):
            assert pdf.xfxQ2(2, x, 10.0) == pytest.approx(
                pdf.xfxQ2(1, x, 10.0), rel=1e-9), (key, x)
    free = f2_backend("CT18ANLO")._pdf
    assert free.xfxQ2(1, 0.5, 10.0) < 0.5 * free.xfxQ2(2, 0.5, 10.0)

    x = np.array([0.1, 0.5, 0.7])
    assert np.allclose(unpolarized_emc_ratio(x, mode="epps21"),
                       [1.012, 0.967, 0.967], atol=0.004)
    assert np.allclose(unpolarized_emc_ratio(x, mode="nnnpdf"),
                       [0.992, 0.987, 0.989], atol=0.004)
    # the callable form is the same object the money plot uses
    assert np.allclose(NuclearF2Ratio("epps21")(x, 5.0),
                       unpolarized_emc_ratio(x, mode="epps21"))
    # the free isoscalar denominator, not the bound proton: the ratio is
    # within 4% of 1 everywhere here, which a proton-over-proton ratio
    # would not be
    assert np.all(np.abs(unpolarized_emc_ratio(x, mode="epps21") - 1) < 0.04)
