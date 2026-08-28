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


#: Text a line may contain that makes a stale energy legitimate: a note
#: EXPLAINING the correction, a deliberately dated reproduction, or the
#: diagnostic runs of tools/pythia8 that were taken at 10 x 50 GeV and are
#: labelled as such.
STALE_ALLOWED = re.compile(r"gamma-matched|2026-08-27|plans/10|before|It gave|"
                           r"pre-2026|slope_b|B = 50|default=50000|diagnostic|"
                           r"superseded|retired|frozen|legacy|"
                           r"rigidity-scaled|stale")

#: The retired 6Li energies, in the forms they actually occur in: the
#: rigidity-scaled 20.5 and 50 GeV/u themselves, the sample files named
#: after them (pythia8_e10_p50_dis.npz and friends, which do not exist --
#: the standing production is e5_p40.8 / e10_p99.5 / e18_p137.5), and the
#: generator command that would make them.  Deliberately narrow: a bare
#: "50" is a slope, a bin count and a percentage all over this tree.
STALE_ENERGY = re.compile(
    r"(?<![\d.])20\.5(?![\d])"
    r"|[pn]50_dis|e10_[pn]50"
    r"|--p-per-nucleon\s+50(?:\.0+)?(?![\d.])"
    r"|(?<![\d.])50(?:\.0+)?\s*GeV/(?:u|nucleon)"
    r"|(?:ion_momentum_per_nucleon|p_per_nucleon)\s*=\s*50(?:\.0+)?(?![\d.])"
    r"|BeamConfig\([^)]*,\s*50(?:\.0+)?\s*\)")

#: An argparse DEFAULT is stale only in the option it belongs to: a bare
#: `default=50` is a report interval, a bin count and a grid size
#: elsewhere in this tree.  This pattern therefore fires only when the
#: option name is on the same line -- and it is the form the sweep exists
#: for, since what was actually wrong on 2026-08-28 was
#: `gen_dis_hfs.py`'s `--p-per-nucleon` DEFAULTING to 50.0 rather than any
#: command line that spelled it out.
STALE_PPN_DEFAULT = re.compile(
    r"p[-_]per[-_]nucleon.{0,80}?default\s*=\s*50(?:\.0+)?(?![\d.])"
    r"|default\s*=\s*50(?:\.0+)?(?![\d.]).{0,80}?p[-_]per[-_]nucleon")

#: Files exempt from the sweep, and why.  `money_delta_2026*.py` are dated
#: reproductions of published figures.  `money_delta_realistic.py` is the
#: same thing without a date in its name: it is frozen on the superseded
#: 5 x 27.5 / 10 x 50 / 18 x 137.5 menu because the L_5sigma numbers of the
#: 2026-07 notes are quoted against it, and its own docstring says so --
#: only its TOP configuration is at a machine configuration, which is the
#: caveat fastsim/README.md carries next to the 66-155 fb^-1/u headline.
STALE_EXEMPT = ("money_delta_2026", "money_delta_realistic")

#: Known offenders in files this test's own round could not edit, kept as
#: exact text rather than line numbers so that the fix silently retires
#: them.  `polligen/hfs.py` quotes the target-mass term E_N - p_N =
#: M^2/(2 p_u) at the retired 50 GeV/u; it is 10.8 / 4.4 / 3.2 MeV at
#: 40.8 / 99.5 / 137.5 GeV/u (reported 2026-08-28).
STALE_REPORTED = ("(8.8 MeV at 50 GeV/u,",)


def test_no_source_file_still_hard_codes_the_stale_energies():
    """A guard, not a physics test.  Lines that DISCUSS the correction are
    allowed; lines that use 20.5 or 50.0 as a 6Li energy are not.

    The sweep covers `tools/` as well as the two packages: the PYTHIA
    driver and its README name the sample files by their beam energies, so
    a stale energy there is a command that cannot be run rather than a
    number that is merely wrong.
    """
    root = pathlib.Path(__file__).resolve().parents[2]
    bad = []
    for pat in ("evgen/scripts/*.py", "evgen/polligen/*.py",
                "fastsim/polli_fastsim/*.py", "fastsim/scripts/*.py",
                "tools/**/*.py", "tools/**/*.md"):
        for f in glob.glob(str(root / pat), recursive=True):
            if any(k in f for k in STALE_EXEMPT):
                continue
            lines = pathlib.Path(f).read_text().splitlines()
            for n, line in enumerate(lines, 1):
                if any(k in line for k in STALE_REPORTED):
                    continue
                # a note EXPLAINING the correction may span a few lines, so
                # look at the neighbourhood rather than the single line
                near = "\n".join(lines[max(0, n - 3):n + 2])
                hit = (STALE_ENERGY.search(line)
                       or STALE_PPN_DEFAULT.search(line))
                if hit and not STALE_ALLOWED.search(near):
                    bad.append("%s:%d %s" % (pathlib.Path(f).name, n,
                                             line.strip()[:70]))
    assert not bad, "stale 6Li energies still present:\n" + "\n".join(bad)


def test_the_stale_energy_guard_would_actually_fire():
    """The guard is a regex over comments; a regex that matches nothing is
    indistinguishable from a clean tree.  These are the forms it exists to
    catch, and the neighbours it must not."""
    for line in ("    --out evgen/samples/pythia8_e10_p50_dis.npz",
                 "  python3 gen_dis_hfs.py --p-per-nucleon 50 --seed 1",
                 "| 5 x 41 | 20.5 | 123 | 2.0 mrad |",
                 "the mid configuration (10 GeV x 50 GeV/u)",
                 "ion_momentum_per_nucleon=50)",
                 "  python3 gen_dis_hfs.py --p-per-nucleon 50.0",
                 "    cfg = BeamConfig(10.0, LI6, 50.0)"):
        assert STALE_ENERGY.search(line), line
    # the argparse default that item C of the 2026-08-28 round actually
    # repaired: caught by name, so that reverting it fails the suite
    for line in ('    ap.add_argument("--p-per-nucleon", type=float, '
                 'default=50.0)',
                 "    ap.add_argument('--p-per-nucleon', default=50)"):
        assert not STALE_ENERGY.search(line), line     # not by number alone
        assert STALE_PPN_DEFAULT.search(line), line
    for line in ("    slope_b = 50.0", "  50 bins in log x", "n_grid=400",
                 "    ap.add_argument('--report-every', default=50000)",
                 "    ap.add_argument('--n-report', default=50)",
                 "  python3 gen_dis_hfs.py --p-per-nucleon 50.5",
                 "sigma = 0.550 ub", "beams at 5 x 40.8 and 10 x 99.5 GeV/u"):
        assert not STALE_ENERGY.search(line), line
        assert not STALE_PPN_DEFAULT.search(line), line
