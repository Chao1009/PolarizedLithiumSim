"""plans/10: the EIC energy menu, and that nothing has drifted off it.

Ion energies are GAMMA-matched to their proton configuration and then
capped by the ring rigidity -- not rigidity-scaled.  Before 2026-08-27 this
repository scaled the two lower configurations as 41*Z/A and 100*Z/A,
putting 6Li at 20.5 and 50 GeV/u, which are not machine configurations.

The last test here is a guard rather than a physics check: it sweeps the
source tree for the stale energies so the correction cannot quietly rot
back in.
"""

import glob
import pathlib
import re

import pytest

from polli_fastsim import beams, spectator


def test_proton_reproduces_the_machine_configurations_exactly():
    """The strongest self-check available: run the ion machinery on a
    PROTON and it must return the configuration energies themselves."""
    got = [c.ion_momentum_per_nucleon for c in beams.default_configs("p")]
    assert got == [41.0, 100.0, 275.0]


def test_ions_are_gamma_matched_below_the_rigidity_cap():
    """Every species sits at ~41 GeV/u at the low configuration -- the same
    SPEED as the 41 GeV proton, not the same rigidity.  Yellow Report
    Table 10.2 lists gold at 41 GeV/u, which is what settles it."""
    for name in ("d", "3He", "6Li", "7Li"):
        lo = beams.default_configs(name)[0].ion_momentum_per_nucleon
        assert 40.5 <= lo <= 41.1, name
        # and NOT the rigidity-scaled value, which for 6Li is 20.5
        ion = beams.IONS[name]
        assert abs(lo - 41.0 * ion.Z / ion.A) > 1.0 or ion.Z == ion.A


def test_top_energies_are_rigidity_capped_and_match_the_published_menu():
    for name, top in (("d", 137.5), ("3He", 183.3), ("6Li", 137.5),
                      ("7Li", 117.9)):
        got = beams.default_configs(name)[2].ion_momentum_per_nucleon
        assert got == pytest.approx(top, abs=0.15), name


def test_mass_per_nucleon_uses_the_physical_nuclear_mass():
    """Binding matters at the 0.5% level, which is exactly what separates a
    41 GeV proton from a gamma-matched 6Li at 40.8 GeV/u."""
    for (name, a, z), m in beams.NUCLEUS_MASS.items():
        assert m == pytest.approx(spectator.NUCLEUS_MASS[(z, a)], abs=1e-9)
        assert beams.IONS[name].mass_per_nucleon == pytest.approx(m / a)
    # 6Li is bound by ~5.3 MeV/nucleon relative to a free proton, which is
    # the 0.5%-level effect that separates a 41 GeV proton from a
    # gamma-matched 6Li at 40.8 GeV/u
    assert beams.PROTON_MASS - beams.LI6.mass_per_nucleon == pytest.approx(
        4.7e-3, abs=1e-3)


def test_no_source_file_still_hard_codes_the_stale_energies():
    """A guard, not a physics test.  Lines that DISCUSS the correction are
    allowed; lines that use 20.5 or 50.0 as a 6Li energy are not."""
    root = pathlib.Path(__file__).resolve().parents[2]
    allow = re.compile(r"gamma-matched|2026-08-27|plans/10|before|It gave|"
                       r"pre-2026|slope_b|B = 50|default=50")
    stale = re.compile(r"(?<![\d.])20\.5(?![\d])")
    bad = []
    for pat in ("evgen/scripts/*.py", "evgen/polligen/*.py",
                "fastsim/polli_fastsim/*.py", "fastsim/scripts/*.py"):
        for f in glob.glob(str(root / pat)):
            if "money_delta_2026" in f:      # frozen dated productions
                continue
            lines = pathlib.Path(f).read_text().splitlines()
            for n, line in enumerate(lines, 1):
                # a note EXPLAINING the correction may span a few lines, so
                # look at the neighbourhood rather than the single line
                near = "\n".join(lines[max(0, n - 3):n + 2])
                if stale.search(line) and not allow.search(near):
                    bad.append("%s:%d %s" % (pathlib.Path(f).name, n,
                                             line.strip()[:70]))
    assert not bad, "stale 6Li energies still present:\n" + "\n".join(bad)
