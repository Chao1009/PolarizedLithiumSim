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
               energies and the 2x-pessimistic P_zz propagation.
  ARTEFACTS    that every figure a report embeds exists, is newer than the
               script that makes it, and is registered in build_report.py;
               that every sample matches a configuration energy; that the
               report numbering agrees across templates, index and builder.
  REFERENCES   that every refs_dict entry with a local file has one, and
               that documents citing refs/ point at files that exist.

Exit status is 0 when everything agrees and 1 otherwise, so it can gate a
commit.  `--verbose` lists the checks that passed as well.

Usage:  python3 tools/consistency_check.py [--verbose]
"""

import argparse
import glob
import json
import math
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "fastsim"))
sys.path.insert(0, str(ROOT / "evgen"))

ISSUES = []
PASSED = []


def check(name):
    """Decorator: run a check, collect what it returns as issues."""
    def deco(fn):
        try:
            bad = fn() or []
        except Exception as exc:                      # a broken check is an issue
            bad = ["check raised %s: %s" % (type(exc).__name__, exc)]
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
    want_hd = {0: (220, 380), 1: (206, 206), 2: (119, 119)}
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

@check("artefacts: every figure a report embeds exists and is registered")
def _():
    bad = []
    src = (ROOT / "reports/build_report.py").read_text()
    registered = set(re.findall(r'"(__[A-Z0-9_]+__)":\s*"([^"]+)"', src))
    tags = {t: p for t, p in registered}
    for f in glob.glob(str(ROOT / "reports/*.template.html")):
        for tag in set(re.findall(r"__[A-Z0-9_]+__", pathlib.Path(f).read_text())):
            if tag not in tags:
                bad.append("%s uses %s, which build_report.py does not define"
                           % (pathlib.Path(f).name, tag))
            elif not (ROOT / tags[tag]).exists():
                bad.append("%s -> %s does not exist" % (tag, tags[tag]))
    return bad


@check("artefacts: every embedded figure is newer than the script that makes it")
def _():
    bad = []
    src = (ROOT / "reports/build_report.py").read_text()
    for _, rel in re.findall(r'"(__[A-Z0-9_]+__)":\s*"([^"]+)"', src):
        png = ROOT / rel
        if not png.exists():
            continue
        stem = png.stem
        for script in glob.glob(str(ROOT / "evgen/scripts/*.py")):
            body = pathlib.Path(script).read_text()
            if stem in body or stem.replace("_6Li", "") in body:
                if png.stat().st_mtime < pathlib.Path(script).stat().st_mtime:
                    bad.append("%s is older than %s -- rerun it"
                               % (rel, pathlib.Path(script).name))
                break
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
    for pkg, pat in (("evgen", r"cd evgen\s+&& python3 -m pytest tests/ -q\s+# (\d+) passed"),
                     ("fastsim", r"cd fastsim && python3 -m pytest tests/ -q\s+# (\d+) passed")):
        m = re.search(pat, txt)
        if not m:
            bad.append("the manual no longer states a %s test count" % pkg)
            continue
        r = subprocess.run([sys.executable, "-m", "pytest", "tests/",
                            "--collect-only", "-q"],
                           capture_output=True, text=True, cwd=str(ROOT / pkg))
        got = re.search(r"(\d+) tests? collected", r.stdout)
        if not got:
            bad.append("could not collect %s tests" % pkg)
        elif int(got.group(1)) != int(m.group(1)):
            bad.append("the manual says %s has %s tests, pytest collects %s"
                       % (pkg, m.group(1), got.group(1)))
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    print("Consistency check -- %d checks\n" % (len(PASSED) + len(ISSUES)))
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
    print("\n%d passed, %d failed" % (len(PASSED), len(ISSUES)))
    return 1 if ISSUES else 0


if __name__ == "__main__":
    sys.exit(main())
