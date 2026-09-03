#!/usr/bin/env python3
"""Whole-repository consistency check: does everything still agree?

This programme produces circulate-able reports whose numbers other people
act on, so a correct-but-stale value in one document is worse than the work
it takes to keep them aligned.  This script is the sweep that catches that.

It checks five kinds of agreement:

  PHYSICS      invariants the simulation must satisfy, independent of any
               published number -- the proton self-check on the ion energy
               menu, the tagging optics giving 1/e, the two mass tables
               agreeing, the divergence species step applying only where
               rigidity binds.
  SOURCES      published values the code carries, against the tables they
               come from (Yellow Report 10.1 / 10.2).
  DRIFT        stale values that a correction should have removed --
               principally the pre-2026-08-27 rigidity-scaled 6Li
               energies and the 2x-pessimistic P_zz propagation -- and
               statements a rewrite must not drop, such as the run-plan
               share every per-year reach is quoted at.
  ARTEFACTS    that every figure a report embeds exists, is registered in
               build_report.py, and is newer both than the script that makes
               it and than every library module that script imports; that
               every sample matches a configuration energy; that the report
               numbering agrees across templates, index and builder.
  REFERENCES   that every refs_dict entry with a local file has one, and
               that documents citing refs/ point at files that exist.

Exit status is 0 when everything agrees and 1 otherwise, so it can gate a
commit.  `--verbose` lists the checks that passed as well.

Usage:  python3 tools/consistency_check.py [--verbose] [--full]

Checks added after 2026-09-02 live one module each under tools/checks/ and are
loaded by name; see the loader at the end of this file.
"""

import argparse
import ast
import glob
import json
import math
import pathlib
import re
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "fastsim"))
sys.path.insert(0, str(ROOT / "evgen"))

ISSUES = []
PASSED = []
SKIPPED = []
# a check returns SKIP when it cannot run in this mode -- the two that
# re-execute the producing scripts run only under --full -- so that the
# default sweep counts them as skipped rather than reporting them passed
SKIP = object()
# --full also runs the checks that re-derive expected numbers by executing
# the producing scripts (minutes rather than seconds); the default sweep,
# which the fastsim test suite runs, must stay under about a minute.
FULL = "--full" in sys.argv


def check(name):
    """Decorator: run a check, collect what it returns as issues."""
    def deco(fn):
        try:
            bad = fn()
        except Exception as exc:                      # a broken check is an issue
            bad = ["check raised %s: %s" % (type(exc).__name__, exc)]
        if bad is SKIP:
            SKIPPED.append(name)
            return fn
        bad = bad or []
        if bad:
            ISSUES.append((name, bad))
        else:
            PASSED.append(name)
        return fn
    return deco


# --- PHYSICS ---------------------------------------------------------------

@check("physics: the ion energy menu is gamma-matched and rigidity-capped")
def _():
    from polli_fastsim import beams
    bad = []
    p = [c.ion_momentum_per_nucleon for c in beams.default_configs("p")]
    if p != [41.0, 100.0, 275.0]:
        bad.append("proton self-check failed: %s, expected [41.0, 100.0, 275.0]" % p)
    for name, top in (("d", 137.5), ("3He", 183.3), ("6Li", 137.5), ("7Li", 117.9)):
        cfg = beams.default_configs(name)
        if abs(cfg[2].ion_momentum_per_nucleon - top) > 0.15:
            bad.append("%s top = %.1f, expected %.1f (rigidity cap)"
                       % (name, cfg[2].ion_momentum_per_nucleon, top))
        if not 40.5 <= cfg[0].ion_momentum_per_nucleon <= 41.1:
            bad.append("%s low = %.1f, expected ~41 (gamma-matched)"
                       % (name, cfg[0].ion_momentum_per_nucleon))
    return bad


@check("physics: the two nuclear mass tables agree")
def _():
    from polli_fastsim import beams, spectator
    return ["%s: %.9f vs %.9f" % (n, m, spectator.NUCLEUS_MASS[(z, a)])
            for (n, a, z), m in beams.NUCLEUS_MASS.items()
            if abs(m - spectator.NUCLEUS_MASS[(z, a)]) > 1e-9]


@check("physics: divergence species step applies only where rigidity binds")
def _():
    from polli_fastsim import beams
    from polligen import reco
    bad = []
    for i, expect in ((0, 1.0), (1, 1.0), (2, math.sqrt(2.0))):
        f = (reco.sigma_theta_for(beams.default_configs("6Li")[i])[0]
             / reco.sigma_theta_for(beams.default_configs("p")[i])[0])
        if abs(f - expect) > 0.02 * expect:
            bad.append("6Li/p sigma ratio at config %d = %.3f, expected %.3f"
                       % (i, f, expect))
    return bad


@check("physics: the tagging optics gives acceptance 1/e at every configuration")
def _():
    from polli_fastsim import beams
    from polligen import reco, coherent as coh
    sc = coh.CoherentScenario()
    bad = []
    for c in beams.default_configs("6Li"):
        st = reco.sigma_theta_tagging(c, slope_b=sc.slope_b)
        acc = float(sc.tag_acceptance_angular(st, c.ion_momentum_per_nucleon))
        if abs(acc - math.exp(-1.0)) > 0.02:
            bad.append("acceptance at the tagging optimum = %.4f, expected 1/e" % acc)
    return bad


@check("physics: the aperture lookup is keyed off the current configurations")
def _():
    from polli_fastsim import beams
    from polligen import reco
    bad = []
    for c in beams.default_configs("6Li"):
        if reco.rp_aperture_for(c.ion_momentum_per_nucleon) is None:
            bad.append("no aperture for the current %.1f GeV/u configuration"
                       % c.ion_momentum_per_nucleon)
    for stale in (20.5, 50.0):
        if reco.rp_aperture_for(stale) is not None:
            bad.append("aperture still resolves for the stale %.1f GeV/u" % stale)
    return bad


# --- SOURCES ---------------------------------------------------------------

@check("sources: the divergence table reproduces Yellow Report Table 10.1")
def _():
    from polli_fastsim import beams
    from polligen import reco
    want = {0: (220, 380), 1: (180, 180), 2: (65, 65)}          # high acceptance
    want_hd = {0: (220, 380), 1: (220, 220), 2: (150, 150)}   # e10 x p100, e18 x p275 (290 b)
    bad = []
    for i, c in enumerate(beams.default_configs("p")):
        got = tuple(round(1e6 * x) for x in reco.sigma_theta_for(c))
        if got != want[i]:
            bad.append("proton HA[%d] = %s, YR Table 10.1 says %s" % (i, got, want[i]))
        got = tuple(round(1e6 * x)
                    for x in reco.sigma_theta_for(c, "high-divergence"))
        if got != want_hd[i]:
            bad.append("proton HD[%d] = %s, YR says %s" % (i, got, want_hd[i]))
    return bad


@check("sources: Table 11.48's divergence column follows from Table 10.1")
def _():
    from polli_fastsim import farforward as ff
    bad = []
    for cfg, p_gev, tab in (("10x100", 100.0, 22.0), ("5x41", 41.0, 14.0)):
        (hd_h, hd_v, _, _), _ = ff.YR_PROTON_DIVERGENCE[cfg]
        rms = math.sqrt(0.5 * (hd_h ** 2 + hd_v ** 2))
        got = 1e-3 * rms * p_gev
        if abs(got - tab) > 0.15 * tab:
            bad.append("%s: Table 10.1 implies dpT = %.1f MeV, Table 11.48 says %.0f"
                       % (cfg, got, tab))
    return bad


@check("sources: the configurations sit in EPIOS's synchronisation windows")
def _():
    """EPIOS pp. 12-13 gives the mechanism -- a +-20 mm radial shift covering
    118 < gamma < 293, plus a 'Blue' arc bypass at gamma = 43.5.  Two of the
    three Yellow Report configurations fall inside; the 100 GeV proton does
    not, and that is a KNOWN CONFLICT between two EIC documents rather than
    a repository error.  This check exists so the conflict stays visible and
    so a NEW configuration cannot be added without noticing."""
    from polli_fastsim import beams
    bad = []
    for e in beams.PROTON_CONFIG_ENERGIES:
        w = beams.epios_window_of(e)
        if w is None and e != 100.0:
            bad.append("the %g GeV configuration (gamma %.1f) is in neither "
                       "EPIOS window and is not the known 100 GeV conflict"
                       % (e, beams.gamma_of(e)))
        if w is not None and e == 100.0:
            bad.append("the 100 GeV conflict has resolved (gamma %.1f now in "
                       "%s) -- update the note in beams.py and plans/10"
                       % (beams.gamma_of(e), w))
    return bad


# --- DRIFT -----------------------------------------------------------------

STALE_ENERGY = re.compile(r"(?<![\d.])20\.5(?![\d])")
ALLOW = re.compile(r"gamma-matched|γ-matched|2026-08-27|plans/10|before|It gave|"
                   r"pre-2026|superseded|then-assumed|rigidity-scaled|not 20|"
                   r"Dated record|diagnostic|slope_b|B = 50|default=50")


@check("drift: no live source file carries the pre-correction 6Li energies")
def _():
    bad = []
    for pat in ("evgen/scripts/*.py", "evgen/polligen/*.py",
                "fastsim/polli_fastsim/*.py", "fastsim/scripts/*.py",
                "tools/**/*.py"):
        for f in glob.glob(str(ROOT / pat), recursive=True):
            if ("money_delta_2026" in f
                    or pathlib.Path(f).name == "consistency_check.py"):
                continue   # frozen dated productions, and this file
            lines = pathlib.Path(f).read_text(errors="ignore").splitlines()
            for n, line in enumerate(lines, 1):
                near = "\n".join(lines[max(0, n - 3):n + 2])
                if STALE_ENERGY.search(line) and not ALLOW.search(near):
                    bad.append("%s:%d %s"
                               % (pathlib.Path(f).relative_to(ROOT), n,
                                  line.strip()[:70]))
    return bad


@check("drift: no report template carries the pre-correction 6Li energies")
def _():
    bad = []
    for f in glob.glob(str(ROOT / "reports/*.template.html")):
        txt = pathlib.Path(f).read_text()
        for m in STALE_ENERGY.finditer(txt):
            near = txt[max(0, m.start() - 260):m.start() + 260]
            if not ALLOW.search(near):
                bad.append("%s: ...%s..."
                           % (pathlib.Path(f).name,
                              txt[max(0, m.start() - 45):m.start() + 45]
                              .replace("\n", " ")))
    return bad


@check("drift: the alpha+d separation table matches the current energies")
def _():
    """Report 4 Table 5 and plans/09 SS9.2 quote the median alpha-d
    separation at the pot plane in MILLIMETRES, which is an angle (~ 1/p_u)
    times R12.  A retired beam energy therefore survives inside a derived
    number that no energy-drift check can scan for: the two lower rows were
    published as 30.1 and 73.4 mm, computed at the rigidity-scaled 20.5 and
    50 GeV/u, and stayed there through the correction to the gamma-matched
    menu (a factor 2).  This check recomputes them from the current
    configurations and compares, so the class of drift is closed rather than
    the instance -- and it caught the second error in the same table on the
    same day, when the correction to the energies dropped the DISPERSIVE
    displacement between the two fragments (farforward.separation_at_pots)
    and came out 5-39% low.  Pinned in fastsim/tests/test_two_hit.py at the
    same numbers, together with the angular-only column; 400k events per
    configuration, seed 7, ~1 s.

    Third error, 2026-08-28 (plans/09 B1): the recomputation itself was
    calling `separation_at_pots` with its DEFAULTS, which are the 18 x 275
    levers, so it validated the documents against the same single-lever
    arithmetic the documents used and could not see that R12 and R34 are
    per configuration.  It now passes `config=cfg`.

    Fourth error, 2026-08-29: 5 x 41 was still on the fallback
    R34 = R12 = 19.24 m, because the vertical plane is shut and the B1
    ladder had nothing to regress there; the zero-insertion scratch
    geometry gives R34 = 4.56 m and the median falls 25.8 -> 17.3 mm.
    The current values are 17.3 / 10.7 / 10.9 mm."""
    import numpy as np
    from polli_fastsim import beams
    from polli_fastsim import farforward as ff
    from polli_fastsim import spectator as sp
    bad, computed = [], {}
    for cfg in beams.default_configs("6Li"):
        ev = sp.breakup_lab_kinematics(sp.LI6_ALPHA_TAG,
                                       cfg.ion_momentum_per_nucleon, 400_000,
                                       rng=np.random.default_rng(7))
        computed[ff.yr_config_key(cfg)] = float(np.median(
            1e3 * ff.separation_at_pots(ev["spectator"], ev["partner"],
                                        config=cfg)))
    docs = (("reports/nanowire_far_forward.template.html",
             r'<tr><td>%s</td><td class="mono">([\d.]+) mm</td>'),
            ("plans/09_nearbeam_nanowire_far_forward.md",
             r'\| %s \| \*\*([\d.]+) mm\*\*'))
    labels = {"5x41": "5 × 41", "10x100": "10 × 100", "18x275": "18 × 275"}
    for doc, pat in docs:
        txt = (ROOT / doc).read_text()
        for key, med in computed.items():
            m = re.search(pat % re.escape(labels[key]), txt)
            if not m:
                bad.append("%s no longer carries an alpha+d separation row "
                           "for %s" % (pathlib.Path(doc).name, labels[key]))
            elif abs(float(m.group(1)) - med) > 0.04 * med:
                bad.append("%s says %s mm for %s, the current energies give "
                           "%.1f mm" % (pathlib.Path(doc).name, m.group(1),
                                        labels[key], med))
    return bad


@check("drift: no report claims the P_zz scale propagates quadratically")
def _():
    """Commit f05d026 published (1+d)^2 - 1 = 4.0/10.3/21.0% at
    d = 2/5/10%.  The propagation is 1:1 (2.0/5.0/10.0%) -- the estimator's
    weights are built from the ASSUMED polarizations, so one power of the
    scale cancels between R and sigma_P^2.  Pinned by
    test_pzz_scale_error_propagates_one_to_one_not_quadratically.  This
    check stops the 2x-pessimistic table coming back."""
    bad, allow = [], re.compile(
        r"corrected|Corrected|wrong|too pessimistic|reach|row below|"
        r"was corrected|earlier version")
    stale = re.compile(r"10\.3\s*(/|%)|21\.0\s*%|4\.0\s*/\s*10\.3")
    for f in glob.glob(str(ROOT / "reports/*.template.html")) + \
             glob.glob(str(ROOT / "plans/*.md")) + \
             glob.glob(str(ROOT / "docs/*.md")):
        txt = pathlib.Path(f).read_text(errors="ignore")
        for m in stale.finditer(txt):
            near = txt[max(0, m.start() - 400):m.start() + 400]
            if "P_zz" in near or "Pzz" in near or "P<sub>zz" in near:
                if not allow.search(near):
                    bad.append("%s: ...%s..."
                               % (pathlib.Path(f).name,
                                  txt[max(0, m.start() - 60):m.start() + 60]
                                  .replace("\n", " ")))
    return bad


@check("drift: the relative-luminosity coefficient is quoted with its convention")
def _():
    """The bias is -[(P1+P2)/(P1-P2)] x delta_ratio = 1/3 for the flip
    plan's (+0.6, -1.2).  The scripts' --rel-lumi-offset d is a RATIO error
    of ~2d, so that convention shows ~2/3.  An unqualified 1.4 was the old
    published number and was simply wrong."""
    bad = []
    stale = re.compile(r"1\.4\s*(x|\u00d7)\s*(\u03b4|delta)")
    for f in glob.glob(str(ROOT / "reports/*.template.html")) + \
             glob.glob(str(ROOT / "plans/*.md")):
        txt = pathlib.Path(f).read_text(errors="ignore")
        for m in stale.finditer(txt):
            near = txt[max(0, m.start() - 400):m.start() + 400]
            if not re.search(r"corrected|Corrected|wrong|row below|"
                             r"not the|was published|conflated", near):
                bad.append("%s: ...%s..."
                           % (pathlib.Path(f).name,
                              txt[max(0, m.start() - 60):m.start() + 60]
                              .replace("\n", " ")))
    return bad


# The sentence every report that quotes a per-year reach has to carry.  A
# pending entry exempts a report a rewrite has not reached yet, and is
# REMOVED the moment the sentence lands: the check fails on an entry whose
# template already has it, so the exemption cannot quietly become
# permanent.  All five reports carry it as of 2026-08-28, so the map is
# empty; keep it, because the next report added starts out without.
RUNPLAN_STATEMENT = "each projection assumes the full luminosity"
RUNPLAN_PENDING = {}


@check("drift: every report quoting a per-year reach states the run-plan share")
def _():
    """Every projection in these reports gives its observable the whole of
    the 10 fb^-1/u year in its own spin configuration, far-forward optics
    and isotope, so the reaches of one table cannot be summed into a year.
    That is a statement about the RUN PLAN, not about any one measurement,
    and a rewrite that drops it turns a set of alternatives into a
    programme.  Any report that quotes a luminosity per nucleon or a
    per-year number therefore has to carry it (plans/07 WP2)."""
    triggers = ("fb⁻¹/u", "one EIC year", "per year")
    bad = []
    for f in sorted(glob.glob(str(ROOT / "reports/*.template.html"))):
        name = pathlib.Path(f).name
        txt = pathlib.Path(f).read_text(errors="ignore")
        has = RUNPLAN_STATEMENT in txt.lower()
        if name in RUNPLAN_PENDING:
            if has:
                bad.append("%s now carries the run-plan statement -- remove "
                           "its RUNPLAN_PENDING entry in %s (%s)"
                           % (name, pathlib.Path(__file__).name,
                              RUNPLAN_PENDING[name]))
            continue
        if any(s in txt for s in triggers) and not has:
            hit = next(s for s in triggers if s in txt)
            bad.append("%s quotes %r but does not say %r"
                       % (name, hit, RUNPLAN_STATEMENT))
    return bad


@check("drift: the PYTHIA samples match the current configuration energies")
def _():
    from polli_fastsim import beams
    pus = [c.ion_momentum_per_nucleon for c in beams.default_configs("6Li")]
    bad = []
    for f in glob.glob(str(ROOT / "evgen/samples/*.npz")):
        m = re.search(r"_[pn]([\d.]+)_", pathlib.Path(f).name)
        if m and not any(abs(float(m.group(1)) - p) < 0.05 for p in pus):
            bad.append("%s is at an energy no configuration uses"
                       % pathlib.Path(f).name)
    return bad


# --- ARTEFACTS -------------------------------------------------------------

# the ".png" is load-bearing: build_report.py's own comment about the
# registry spells out a "__TAG__": "path" pair that is not a figure
FIG_RE = re.compile(r'"(__[A-Z0-9_]+__)":\s*"([^"]+\.png)"')

# the two library trees every figure's numbers come out of, and the
# digitized curves polarized.py reads: data a figure depends on exactly
# as it depends on a module
PKG_ROOTS = (("polligen", ROOT / "evgen" / "polligen"),
             ("polli_fastsim", ROOT / "fastsim" / "polli_fastsim"))
POLARIZED = ROOT / "fastsim" / "polli_fastsim" / "polarized.py"
POLARIZED_DATA = tuple(sorted(
    pathlib.Path(p) for p in
    glob.glob(str(ROOT / "fastsim/polli_fastsim/data/*.csv"))))

_BODIES = {}
_FIG_SCRIPT = {}
_IMPORTS = {}
_DEPS = {}


def registered_figures():
    """Every ("__TAG__", "path.png") pair build_report.py registers, in file
    order -- the embedded ones and UNEMBEDDED_FIGS alike."""
    return FIG_RE.findall((ROOT / "reports/build_report.py").read_text())


def _script_bodies():
    """evgen/scripts/*.py read once, in a deterministic order."""
    if not _BODIES:
        for sc in sorted(glob.glob(str(ROOT / "evgen/scripts/*.py"))):
            _BODIES[sc] = pathlib.Path(sc).read_text()
    return _BODIES


def script_of_figure(rel):
    """Which evgen/scripts/*.py draws the figure at repository path `rel`,
    or None when no script names it.  Both staleness checks below use this
    one map.

    A "--tag" suffix (money_cos2phi_reco_6Li_hfscal.png) is not a literal
    in any script, so three candidates are tried, most specific first: the
    full file stem (which is what an output_stem guard spells out), the
    isotope-less one, then that without the suffix.  Each is matched as a
    figure-file NAME rather than as the prefix of a longer one, so
    money_cos2phi cannot claim money_cos2phi_reco_6Li.png; two spellings,
    in this order: the literal "...png", and (since 2026-08-28, when the
    published stems went behind an output_stem guard) the bare stem as a
    quoted string."""
    if rel in _FIG_SCRIPT:
        return _FIG_SCRIPT[rel]
    stem = pathlib.Path(rel).stem
    bare = stem.replace("_6Li", "").replace("_7Li", "")
    cands = [stem, bare]
    if re.sub(r"_[a-z0-9]+$", "", bare) != bare:
        cands.append(re.sub(r"_[a-z0-9]+$", "", bare))
    bodies = _script_bodies()
    found = None
    for cand in cands:
        short = re.escape(cand) + r"(?![A-Za-z0-9]|_[A-Za-z])"
        for pat in (short + r"[^\"'\s]*\.png", short + r"[\"']"):
            for sc in bodies:
                if re.search(pat, bodies[sc]):
                    found = sc
                    break
            if found:
                break
        if found:
            break
    _FIG_SCRIPT[rel] = found
    return found


def _own_package(path):
    """The package a file lives in, for resolving `from . import x`."""
    for pkg, root in PKG_ROOTS:
        if path.parent == root:
            return pkg
    return None


def _resolve(dotted):
    """A dotted import name -> the files of ours it reaches.  Importing a
    submodule executes the package __init__, which is why the package
    itself is a dependency of every submodule import."""
    for pkg, root in PKG_ROOTS:
        if dotted != pkg and not dotted.startswith(pkg + "."):
            continue
        out = [root / "__init__.py"]
        if dotted != pkg:
            f = root / (dotted[len(pkg) + 1:].split(".")[0] + ".py")
            if f.is_file():
                out.append(f)
        return [p for p in out if p.is_file()]
    return []


def _imports_of(path):
    """The dotted names `path` imports, read with ast (so a commented-out
    or merely quoted import does not count) and including the ones written
    inside functions.  Relative imports are resolved against the package
    the file lives in, and `from pkg import mod` yields `pkg.mod` as well
    as `pkg`."""
    if path in _IMPORTS:
        return _IMPORTS[path]
    names = set()
    try:
        tree = ast.parse(path.read_text(), filename=str(path))
    except (OSError, SyntaxError):
        tree = None
    own = _own_package(path)
    for node in ast.walk(tree) if tree else ():
        if isinstance(node, ast.Import):
            names.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            base = node.module or ""
            if node.level:
                base = ".".join(x for x in (own, base) if x)
            if base:
                names.add(base)
                names.update("%s.%s" % (base, a.name) for a in node.names)
    _IMPORTS[path] = names
    return names


def module_deps(script):
    """Every file under evgen/polligen or fastsim/polli_fastsim that
    running `script` imports, transitively, plus the digitized CSVs
    polarized.py reads at import time.  Sibling scripts are followed (one
    figure script imports another's helpers) but not returned: the
    script's own mtime is the other check's business."""
    script = pathlib.Path(script)
    if script in _DEPS:
        return _DEPS[script]
    seen, deps, queue = set(), set(), [script]
    while queue:
        cur = queue.pop()
        if cur in seen:
            continue
        seen.add(cur)
        for dotted in sorted(_imports_of(cur)):
            for f in _resolve(dotted):
                deps.add(f)
                queue.append(f)
                if f == POLARIZED:
                    deps.update(POLARIZED_DATA)
            sib = script.parent / (dotted + ".py")
            if "." not in dotted and sib.is_file():
                queue.append(sib)
    _DEPS[script] = sorted(deps)
    return _DEPS[script]


def _stamp(mtime):
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(mtime))


@check("artefacts: every figure a report embeds exists and is registered")
def _():
    bad = []
    tags = dict(registered_figures())
    for f in sorted(glob.glob(str(ROOT / "reports/*.template.html"))):
        for tag in sorted(set(re.findall(r"__[A-Z0-9_]+__",
                                         pathlib.Path(f).read_text()))):
            if tag not in tags:
                bad.append("%s uses %s, which build_report.py does not define"
                           % (pathlib.Path(f).name, tag))
            elif not (ROOT / tags[tag]).exists():
                bad.append("%s -> %s does not exist" % (tag, tags[tag]))
    return bad


@check("artefacts: every embedded figure is newer than the script that makes it")
def _():
    bad = []
    for _tag, rel in registered_figures():
        png = ROOT / rel
        script = script_of_figure(rel)
        if not png.exists() or script is None:
            continue
        if png.stat().st_mtime < pathlib.Path(script).stat().st_mtime:
            bad.append("%s is older than %s -- rerun with the manual's command"
                       % (rel, pathlib.Path(script).name))
    return bad


@check("artefacts: every embedded figure is newer than every module its script "
       "imports")
def _():
    """A figure goes stale when its own script moves on -- and just as
    surely when a library module that script imports does, which the check
    above cannot see.  Development run 14 changed polarized.py, farforward.py
    and half of polligen under figures whose scripts were never touched.
    The import graph is read with ast and followed transitively through
    polligen and polli_fastsim (a package __init__ counts: importing a
    submodule runs it), with polarized.py carrying the digitized CSVs it
    reads."""
    bad = []
    for _tag, rel in registered_figures():
        png = ROOT / rel
        script = script_of_figure(rel)
        if not png.exists() or script is None:
            continue
        t_png = png.stat().st_mtime
        newer = sorted((d.stat().st_mtime, d) for d in module_deps(script)
                       if d.stat().st_mtime > t_png)
        if not newer:
            continue
        t_mod, mod = newer[-1]
        bad.append("%s (%s) is older than %s (%s)%s -- rerun with the "
                   "manual's command"
                   % (rel, _stamp(t_png), mod.relative_to(ROOT),
                      _stamp(t_mod),
                      "" if len(newer) == 1
                      else ", and %d more module%s" % (len(newer) - 1,
                                                       "" if len(newer) == 2
                                                       else "s")))
    return bad


@check("artefacts: report numbering agrees across builder, templates and index")
def _():
    bad = []
    src = (ROOT / "reports/build_report.py").read_text()
    order = re.findall(r'"stem":\s*"([a-z_0-9]+)",\s*"number":\s*(\d+)', src)
    if [int(n) for _, n in order] != list(range(len(order))):
        bad.append("build_report numbers are %s, not 0..N"
                   % [n for _, n in order])
    index = (ROOT / "reports/index.html").read_text()
    idx_nums = [int(x) for x in re.findall(r'class="num">Report (\d+)<', index)]
    if idx_nums != [int(n) for _, n in order]:
        bad.append("index.html lists %s, builder has %s"
                   % (idx_nums, [int(n) for _, n in order]))
    for stem, num in order:
        tpl = ROOT / ("reports/%s.template.html" % stem)
        if not tpl.exists():
            bad.append("no template for %s" % stem)
            continue
        m = re.search(r"<b>Report (\d+)</b>", tpl.read_text())
        if not m:
            bad.append("%s carries no report number in its eyebrow" % stem)
        elif int(m.group(1)) != int(num):
            bad.append("%s says Report %s, builder says %s"
                       % (stem, m.group(1), num))
        for ext in ("html", "pdf"):
            if not (ROOT / ("reports/%s.%s" % (stem, ext))).exists():
                bad.append("reports/%s.%s missing -- rebuild" % (stem, ext))
    return bad


@check("artefacts: every built report is newer than its template")
def _():
    bad = []
    for tpl in glob.glob(str(ROOT / "reports/*.template.html")):
        t = pathlib.Path(tpl)
        for ext in ("html", "pdf"):
            out = t.with_name(t.name.replace(".template.html", "." + ext))
            if out.exists() and out.stat().st_mtime < t.stat().st_mtime:
                bad.append("%s is older than its template -- rebuild" % out.name)
    return bad


@check("artefacts: the manual's test counts match what pytest collects")
def _():
    """The reproduction manual quotes the size of each suite as a check a
    reader can run.  If it drifts, the first thing a new reader does fails."""
    import subprocess
    txt = (ROOT / "docs/reproduction_manual.md").read_text()
    bad = []
    collected = {}
    for pkg, pat in (("evgen", r"cd evgen\s+&& python3 -m pytest tests/ -q\s+# (\d+) passed"),
                     ("fastsim", r"cd fastsim && python3 -m pytest tests/ -q\s+# (\d+) passed")):
        r = subprocess.run([sys.executable, "-m", "pytest", "tests/",
                            "--collect-only", "-q"],
                           capture_output=True, text=True, cwd=str(ROOT / pkg))
        got = re.search(r"(\d+) tests? collected", r.stdout)
        if not got:
            bad.append("could not collect %s tests" % pkg)
        else:
            collected[pkg] = int(got.group(1))
        m = re.search(pat, txt)
        if not m:
            bad.append("the manual no longer states a %s test count" % pkg)
        elif pkg in collected and collected[pkg] != int(m.group(1)):
            bad.append("the manual says %s has %s tests, pytest collects %d"
                       % (pkg, m.group(1), collected[pkg]))
    # the same count is quoted in the package quick-starts, which drifted
    # unguarded until 2026-08-28 (evgen/README.md said 270 against 276)
    for pkg, pat in (("evgen", r"python3 -m pytest tests/ -q\s+# (\d+) tests"),
                     ("fastsim", r"python3 -m pytest tests/ -q\s+# (\d+) tests")):
        rd = ROOT / pkg / "README.md"
        if not rd.exists():
            continue
        m = re.search(pat, rd.read_text())
        if m and pkg in collected and collected[pkg] != int(m.group(1)):
            bad.append("%s/README.md says %s tests, pytest collects %d"
                       % (pkg, m.group(1), collected[pkg]))
    # and the manual's headline total must be the sum of the two suites
    m = re.search(r"^(\d+) tests, all of which run", txt, re.M)
    if not m:
        bad.append("the manual no longer states a total test count")
    elif len(collected) == 2 and sum(collected.values()) != int(m.group(1)):
        bad.append("the manual's total is %s, the two suites collect %d"
                   % (m.group(1), sum(collected.values())))
    return bad


@check("artefacts: the 6Li configuration energies appear consistently in the docs")
def _():
    """Every document that names the three 6Li energies must name the same
    three.  This is the check that a hand-edit missed a file."""
    from polli_fastsim import beams
    pus = [c.ion_momentum_per_nucleon for c in beams.default_configs("6Li")]
    triple = re.compile(r"(\d+\.?\d*)\s*/\s*(\d+\.?\d*)\s*/\s*(\d+\.?\d*)\s*GeV/u")
    bad = []
    for pat in ("reports/*.template.html", "plans/*.md", "docs/reproduction_manual.md",
                "evgen/README.md"):
        for f in glob.glob(str(ROOT / pat)):
            txt = pathlib.Path(f).read_text()
            for m in triple.finditer(txt):
                vals = [float(x) for x in m.groups()]
                if abs(vals[2] - pus[2]) > 0.2:       # not a 6Li energy triple
                    continue
                near = txt[max(0, m.start() - 300):m.start() + 120]
                if ALLOW.search(near):
                    continue
                if any(abs(v - p) > 0.2 for v, p in zip(vals, pus)):
                    bad.append("%s quotes %s GeV/u, configurations are %s"
                               % (pathlib.Path(f).name,
                                  " / ".join("%g" % v for v in vals),
                                  " / ".join("%g" % p for p in pus)))
    return bad


# --- REFERENCES ------------------------------------------------------------

@check("references: every refs_dict entry with a local file has one")
def _():
    d = json.loads((ROOT / "refs/refs_dict.json").read_text())
    bad = []
    for key, e in d["entries"].items():
        f = e.get("file")
        if not f:
            continue
        for part in str(f).split(";"):
            if not (ROOT / part.strip()).exists():
                bad.append("%s -> %s missing" % (key, part.strip()))
    return bad


@check("references: every refs/ file cited in a document exists")
def _():
    bad = []
    seen = set()
    for pat in ("reports/*.template.html", "plans/*.md", "docs/*.md",
                "evgen/README.md", "refs/README.md"):
        for f in glob.glob(str(ROOT / pat)):
            for m in re.findall(r"refs/([A-Za-z0-9_.\-]+\.pdf)",
                                pathlib.Path(f).read_text()):
                if m in seen:
                    continue
                seen.add(m)
                if not (ROOT / "refs" / m).exists():
                    bad.append("%s cites refs/%s, which does not exist"
                               % (pathlib.Path(f).name, m))
    return bad


@check("references: every script the manual names exists")
def _():
    txt = (ROOT / "docs/reproduction_manual.md").read_text()
    bad = []
    for m in sorted(set(re.findall(r"(?:scripts|tools)/[A-Za-z0-9_/]+\.py", txt))):
        hits = [p for p in (ROOT / "evgen" / m, ROOT / "fastsim" / m, ROOT / m)
                if p.exists()]
        if not hits:
            bad.append("the manual names %s, which does not exist" % m)
    return bad


# --- extension checks: tools/checks/*.py ------------------------------------
# Every module in tools/checks/ registers its checks with the same @check
# decorator at import time (the decorator runs the check as it is defined), so
# they run here, after the built-in checks above.  The loader injects the
# decorator and this module's namespace (as `checker`) into each module before
# executing it, so a module never re-imports this file; a module that fails to
# import is itself reported as an issue.
CHECKS_DIR = pathlib.Path(__file__).resolve().parent / "checks"


def _load_extension_checks():
    import importlib.util
    for path in sorted(CHECKS_DIR.glob("*.py")):
        if path.name.startswith("_"):
            continue
        try:
            spec = importlib.util.spec_from_file_location("polli_check_" + path.stem, path)
            mod = importlib.util.module_from_spec(spec)
            mod.__dict__["check"] = check
            mod.__dict__["checker"] = sys.modules[__name__]
            sys.modules[spec.name] = mod
            spec.loader.exec_module(mod)
        except Exception as exc:                  # a broken module is an issue
            ISSUES.append(("extension %s" % path.name,
                           ["failed to import: %s: %s" % (type(exc).__name__, exc)]))


_load_extension_checks()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--full", action="store_true",
                    help="also run the slow checks that re-execute the producing scripts")
    args = ap.parse_args()

    print("Consistency check -- %d checks\n"
          % (len(PASSED) + len(ISSUES) + len(SKIPPED)))
    if args.verbose:
        for name in PASSED:
            print("  ok    %s" % name)
        print()
    for name, bad in ISSUES:
        print("  FAIL  %s" % name)
        for b in bad[:12]:
            print("          - %s" % b)
        if len(bad) > 12:
            print("          ... and %d more" % (len(bad) - 12))
    if SKIPPED:
        for name in SKIPPED:
            print("  skip  %s -- runs with --full" % name)
        print()
    print("\n%d passed, %s%d failed"
          % (len(PASSED),
             ("%d skipped, " % len(SKIPPED)) if SKIPPED else "",
             len(ISSUES)))
    return 1 if ISSUES else 0


if __name__ == "__main__":
    sys.exit(main())
