"""derived: the manual's expected numbers, re-derived rather than re-read.

Motivated by `docs/consistency_review_2026-09-02.md` section 5.2, item 1
("Expected numbers re-derived, not re-read").  The review's diagnosis is
that a number moves in the code, the owning report section and the
manual's own expected-output block are restated, and every other copy is
left behind; the manual is the document the rest of the programme reads
its expected values off, so it is the one place where "still true today"
can be settled by running the producing script.

TWO CHECKS
----------
`derived: every manual command names a script whose argparse takes its
flags` runs in the DEFAULT sweep.  It reads every command line out of the
manual's fenced ``bash`` blocks (backslash continuations joined) and out
of the command column of the section 7 table -- 85 lines naming 31
scripts as this is written -- resolves each script under `evgen/` or
`fastsim/`, and asserts that every long flag on the line is one the
script accepts.  The flag names are read from the script source with
`ast` (no import, so it costs nothing); only a script carrying a flag the
source scan cannot see is actually executed as `<script> --help`, and a
flag is reported missing only after `--help` has confirmed it.  No script
needs the fallback today, so the check costs about 0.2 s.

`derived: the manual's expected numbers are what the scripts print today`
runs only under `--full`.  It executes the fifteen runs listed below into
a scratch directory outside the repository (`$POLLI_REDERIVE_DIR`, else
`<tmp>/polli_rederive`; every script that has one is given `--outdir`
there, so no published figure is touched), and compares each of the 55
CLAIMS -- a span of the manual anchored on its own wording -- with the
stdout of the command that span is about: 448 numbers in all.  A manual
number fails when nothing the run printed agrees with it at the manual's
own precision, allowing +-1 in its last digit (equivalently: when nothing
printed rounds to the manual's figure or to either neighbour of it).  A
number is also read against a hundredth of itself, but only where the
manual writes it as a per cent -- last in its slash list before a literal
`%` -- since several scripts print as a fraction what the manual quotes
as per cent.  Runtime is about two minutes, of which 49 s is
`money_polemc --pdf grid`.

Sensitivity, measured on this repository: of 896 single-token mutations
of +-3 units in the last manual digit, 63% are caught; the rest land on
some other number the same run prints, which is the ceiling for "is this
value anywhere in the output" and the reason the comparison is scoped per
claim rather than per row.  A re-derivation that moves a whole triple or
row -- which is what section 5.1 shows actually happens -- is caught with
near certainty.

WHAT IS RUN, AND WHAT IS NOT
----------------------------
Run: `tagging_acceptance.py`; `money_polemc.py --ion 7Li` on both `--pdf
grid` (48 s) and `--pdf toy`; the section 3.1b snippet, taken verbatim
out of the manual's own fenced block rather than copied here;
`nearbeam_aperture_scan.py`; `tagging_optics.py`; `money_cos2phi.py`;
`money_cos2phi_reco.py --syst-scan` and `--unfold-scan`;
`target_mass_bound.py`; `nearbeam_reach_gain.py --n-mc 2000000`;
`tagged_polarimetry_7li.py --config 0/1/2`; `money_tagged_azz.py --events
400000`.

Excluded, with the reason:

  * everything section 10 times at over ~60 s: `money_polemc --emc-band`
    (22 s is fine but the manual quotes it only through the band, which
    the graded runs do not print), `target_mass_bound --pdf grid` (209 s),
    `money_delta_pdfgrid` (293 s), `hfs_resolution` on PYTHIA (98 s),
    `money_cos2phi_reco --y-source hfs` and the `_hfscal` 5R/7R rows
    (20 s + a PYTHIA sample), `nearbeam_zid_power` (48 s but no manual
    row of numbers), `nearbeam_two_hit` (27 s), the eight-seed ISR
    average (78 s) and every PYTHIA/BeAGLE/npsim row of sections 5.  The
    review itself proposes the `_hfscal` and `hfs_resolution` runs as a
    second, opt-in tier; this module is the first tier.
  * `target_mass_bound.py` is run at its published arguments (no flags,
    2 s, all four blocks) rather than at the review's "--blocks 2 3":
    the flag is spelled `--blocks 23`, and section 7 also quotes the
    block-1 kinematic cap and the block-4 tagged overlay, so the whole
    run is both cheaper to justify and stricter.
  * `money_cos2phi_reco.py --syst-scan` is run, but the manual publishes
    no expected numbers for the nuisance table it prints (section 4.3
    describes it in words only), so the only claim checked against it is
    the 5R block it reprints.

SCOPING, so that the check has no false positives on this repository
-------------------------------------------------------------------
A section 7 row is not a list of one run's numbers: it also carries the
values the same quantity had before a run changed it, the values other
flags of the same script return, ratios the reader is expected to form,
and library constants no script prints.  Comparing a whole row would
therefore fail on correct text.  Each CLAIM below is instead a span
anchored on the manual's own wording, chosen so that every number inside
it is one the command prints, and `NOISE` strips the document furniture
that is a name rather than a value -- ISO dates, beam-configuration
labels ("5 x 41", "10 x 99.5"), section/plan/report/table pointers,
isotope names, mixed identifiers such as CT18ANLO or EPPS21, and
word-hyphen-number forms such as "twist-3".  The stripping runs AFTER
"x10^-n" has been resolved (including an exponent shared by a whole slash
list, "7.4 / 4.4 / 9.5 / 9.5 x10^-3"), or the configuration-label rule
would eat the "5 x10" of a scientific value.  An
anchor that stops matching is itself reported, so a rewrite of the manual
cannot silently drop a claim.  Deliberately outside the spans, all of it
correct text: historical readings ("the same run before the branch
existed gave ...", "the fit that bin used to get returned a_t = -1.56 +-
2.22"); other commands' outputs (`--emc-baseline cbt`, `--emc-mode
constant`, `--emc-band`, `--isotope 7Li`, `--levers 18x275`, `--fit
likelihood`, `--events 8000000`, `--optics legacy`); numbers the manual
says outright are not on the output line (money plot 4's <|cos theta_k|>,
"from the probe of section 4.1"); ratios the reader forms from two
printed numbers (section 4.5b's "2.5 / 6.9 / 3.9x below the 0.20 GeV
reference" and "IR-8's ~20% worth 3.6 / 10.2 / 5.7x", both of which
section 7 quotes in the printed form 2x/7x/4x and 4x/10x/6x, which IS
checked); the section 3.2 decomposition in points (1.51 / 0.13 / 0.20),
which is quoted two digits finer than the `acc` dict is printed; and
library constants (D/R12, P_D, the pot blind block) that belong to the
canonical-triples check, not to this one.  Appendix A revision rows and
plans/00 run sections are not read by this module at all.
"""

import ast
import bisect
import os
import pathlib
import re
import subprocess
import sys
import tempfile

ROOT = checker.ROOT
MANUAL = ROOT / "docs/reproduction_manual.md"
MANUAL_TEXT = MANUAL.read_text()
MANUAL_LINE_STARTS = [0]
for _ln in MANUAL_TEXT.splitlines(True):
    MANUAL_LINE_STARTS.append(MANUAL_LINE_STARTS[-1] + len(_ln))


def _line_of(offset):
    return bisect.bisect_right(MANUAL_LINE_STARTS, offset)


# --- the commands ----------------------------------------------------------
# key -> (working directory, argv after "python3", takes --outdir)

COMMANDS = {
    "tagging_acceptance":
        ("fastsim", ["scripts/tagging_acceptance.py"], True),
    "money_polemc_grid":
        ("fastsim", ["scripts/money_polemc.py", "--ion", "7Li", "--pdf", "grid"], True),
    "money_polemc_toy":
        ("fastsim", ["scripts/money_polemc.py", "--ion", "7Li", "--pdf", "toy"], True),
    "sec31b":
        ("fastsim", None, False),                 # the manual's own snippet
    "nearbeam_aperture_scan":
        ("evgen", ["scripts/nearbeam_aperture_scan.py"], True),
    "tagging_optics":
        ("evgen", ["scripts/tagging_optics.py"], True),
    "money_cos2phi":
        ("evgen", ["scripts/money_cos2phi.py"], True),
    "money_cos2phi_reco_syst":
        ("evgen", ["scripts/money_cos2phi_reco.py", "--syst-scan"], True),
    "money_cos2phi_reco_unfold":
        ("evgen", ["scripts/money_cos2phi_reco.py", "--unfold-scan"], True),
    "target_mass_bound":
        ("evgen", ["scripts/target_mass_bound.py"], False),
    "nearbeam_reach_gain":
        ("evgen", ["scripts/nearbeam_reach_gain.py", "--n-mc", "2000000"], True),
    "tagged_polarimetry_7li_0":
        ("evgen", ["scripts/tagged_polarimetry_7li.py", "--config", "0"], True),
    "tagged_polarimetry_7li_1":
        ("evgen", ["scripts/tagged_polarimetry_7li.py", "--config", "1"], True),
    "tagged_polarimetry_7li_2":
        ("evgen", ["scripts/tagged_polarimetry_7li.py", "--config", "2"], True),
    "money_tagged_azz":
        ("evgen", ["scripts/money_tagged_azz.py", "--events", "400000"], True),
}


def _scratch():
    d = os.environ.get("POLLI_REDERIVE_DIR")
    base = pathlib.Path(d) if d else pathlib.Path(tempfile.gettempdir()) / "polli_rederive"
    base = base.resolve()
    if base == ROOT or ROOT in base.parents:      # never into the repository
        raise RuntimeError("POLLI_REDERIVE_DIR must be outside %s" % ROOT)
    return base


def _sec31b_snippet():
    """The section 3.1b block, read out of the manual so that it is the
    manual's own snippet that is run and not a copy of it here."""
    m = re.search(r"```bash\ncd fastsim && python3 - <<'PY'\n(.*?)\nPY\n```",
                  MANUAL_TEXT, re.S)
    return m.group(1) if m else None


# --- numbers ---------------------------------------------------------------

SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
       "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-"}

NUM_RE = re.compile(r"(?<![0-9.])[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
RUN_EXP_RE = re.compile(
    r"(?P<run>(?:[+-]?\d+(?:\.\d+)?\s*[/,]\s*)*[+-]?\d+(?:\.\d+)?)\s*Ê(?P<exp>-?\d+)")

# document furniture: names that happen to contain digits, never outputs
NOISE = [
    re.compile(r"\d{4}-\d\d-\d\d"),                      # ISO dates
    re.compile(r"§\s*\d+(?:\.\d+)*[a-z]?"),              # section pointers
    re.compile(r"\bplans/\d+"), re.compile(r"\bWP\d+"),
    re.compile(r"\b(?:Report|Table|Figure|item|run)\s+\d+"),
    re.compile(r"\bIR-\d"),
    re.compile(r"[⁶⁷]Li|\b\d+Li\b|\b\d+He\b"),           # isotope names
    re.compile(r"\b\d{1,2}\s*[×x]\s*\d{2,3}(?:\.\d)?\b"),  # 5 x 41, 10 x 99.5
    re.compile(r"\b[A-Za-z]+\d+[A-Za-z0-9]*\b"),         # CT18ANLO, EPPS21, F5
    re.compile(r"(?<![0-9][eE])(?<=[A-Za-z])-\d+(?![.\d])"),  # twist-3, not 7.4e-3
]


def _normalise(span):
    """Fold the unicode minus, resolve every "x10^-n" -- including one
    shared by a whole slash list ("7.4 / 4.4 / 9.5 / 9.5 x10^-3") -- and
    only then strip the furniture, so that the "5 x 41" rule cannot eat
    the "5 x10" of a scientific value."""
    span = span.replace("−", "-")                   # MINUS SIGN
    span = re.sub(r"[×x]\s*10([⁻]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+)",
                  lambda m: "Ê" + "".join(SUP[c] for c in m.group(1)), span)

    def _spread(m):
        exp = m.group("exp")
        return re.sub(r"([+-]?\d+(?:\.\d+)?)", r"\1e" + exp, m.group("run"))
    span = RUN_EXP_RE.sub(_spread, span)
    for pat in NOISE:
        span = pat.sub(" ", span)
    return span


PCT_AHEAD = re.compile(r"[\s/,\d.]*%")


def _manual_numbers(span):
    """[(text, value, tolerance, per-cent?)] for every number the span
    states.  A token counts as a per-cent form -- so that it may also be
    read against a hundredth of itself -- only when it is the last number
    before a literal '%', alone or at the end of a slash list."""
    out, norm = [], _normalise(span)
    for m in NUM_RE.finditer(norm):
        tok = m.group(0)
        try:
            val = float(tok)
        except ValueError:
            continue
        mant, _, exp = tok.lower().partition("e")
        dec = len(mant.split(".")[1]) if "." in mant else 0
        unit = 10.0 ** (-dec + (int(exp) if exp else 0))
        out.append((tok, val, 1.5 * unit,
                    bool(PCT_AHEAD.match(norm, m.end()))))
    return out


def _output_numbers(text):
    vals = set()
    for tok in NUM_RE.findall(text):
        try:
            vals.add(float(tok))
        except ValueError:
            pass
    return sorted(vals)


def _present(val, tol, pool, pct):
    """Is `val` among `pool` at the manual's precision?  `pct` also allows
    the manual's per-cent form of a fraction the script prints."""
    for target in ((val,) if not pct else (val, val / 100.0)):
        i = bisect.bisect_left(pool, target - tol)
        if i < len(pool) and pool[i] <= target + tol:
            return True
    return False


def _nearest(val, pool, pct):
    best = None
    for target in ((val,) if not pct else (val, val / 100.0)):
        for v in pool:
            d = abs(v - target)
            if best is None or d < best[0]:
                best = (d, v)
    return "none" if best is None else ("%g" % best[1])


# --- the claims ------------------------------------------------------------
# (command keys the claim may be verified against, anchor regex)
# The anchor is matched against the whole manual with re.S; the span it
# matches is the text whose numbers must all be in the run's stdout.

S = re.S
CLAIMS = [
    # --- tagging_acceptance.py (section 3.2 block, section 7 row) ---
    (["tagging_acceptance"],
     r"6Li α-tag +YR high-acceptance.*?7Li t-tag +YR high-acceptance[^\n]*"),
    (["tagging_acceptance"],
     r"Almost\s+all of that row is route 6 \(0\.745 / 0\.915 / 0\.933\)"),
    (["tagging_acceptance"],
     r"YR high-acceptance 0\.0186 / 0\.0168 / 0\.0261, tagging optics [\d. /]+,"
     r" legacy 73 μrad 0\.1679 / 0\.0262 / 0\.0276"),

    # --- money_polemc.py --ion 7Li --pdf grid ---
    (["money_polemc_grid"],
     r"  x=0\.09: dDR\(10/fb\).*?published-curve separation:[^\n]*"),
    (["money_polemc_grid"],
     r"δΔR 0\.0423 .*?best bin x = 0\.355 at 0\.45 σ and 1\.43 σ, 0 of 3 bins above 1 σ"),
    (["money_polemc_grid"],
     r"0\.03105 with s_CBT = \n?0?\.?5322 and s_TMT = 0\.2113"),
    (["money_polemc_grid"],
     r"0\.01372 with 0\.\n?2351 / 0\.0933 .*?0\.06459 with 1\.1070 / 0\.4395"),
    (["money_polemc_grid"],
     r"max 0\.0078 over 0\.028 < x < 0\.3 against the 0\.0212–0\.0243 the transferred"
     r" pair shows there, and m\n?ax 0\.1048 above x = 0\.35"),
    (["money_polemc_grid"],
     r"the grid-input δΔR of Report 0's Table 3 is 0\.042 / 0\.040 / 0\.060 / 0\.187"),

    # --- money_polemc.py --ion 7Li --pdf toy ---
    (["money_polemc_toy"],
     r"`--pdf toy`, which writes the published PNG: errors 0\.0477 / 0\.0509 / 0\.0615"
     r" / 0\.1224, best bin x = 0\.089 at 0\.49 σ, valen\n?ce window x = 0\.355 at"
     r" 0\.38 σ and 1\.22 σ"),
    (["money_polemc_toy"],
     r"\(0\.38 and 1\.22 σ with `--pdf toy`"),

    # --- the section 3.1b snippet ---
    (["sec31b"],                       # the four x are the snippet's input, not its output
     r"It prints 2\.25 / 1\.69 / 1\.41 / 1\.14 for CBT against 1\.01 / 0\.98 / 1\.00 /\s+1\.08 for TMT"),
    (["sec31b"],
     r"a CBT minimum of 1\.06 at\s+x = 0\.696, and the two zeros of the denominator"
     r" at x = 0\.280 and 0\.840"),
    (["sec31b"],
     r"the transferred nuclear-matter depletion is\s+0\.021 / 0\.026 / 0\.044"),
    (["sec31b"],
     r"against 0\.041 / 0\.043 / 0\.050 for the transferred\s+mean-field calculation,"
     r" so ΔR separates by 0\.021 at x = 0\.36, 0\.020 at\s+0\.40, 0\.018 at 0\.45"
     r" and 0\.006 at 0\.65"),
    (["sec31b"],
     r"unpolarized curve, 0\.034 / 0\.048 / 0\.087"),

    # --- nearbeam_aperture_scan.py ---
    (["nearbeam_aperture_scan"],
     r"silicon / YR HA envelope / tagging envelope: 9\.4×10⁻¹⁰.*?"
     r"1\.2×10⁻⁵ / 7\.1×10⁻¹⁴ / 0\.33"),
    (["nearbeam_aperture_scan"],
     r"horizontal 2\.50 / 1\.51 / 0\.53 mrad"),
    (["nearbeam_aperture_scan"],
     r"0\.0170 / 0\.0177 / 0\.3159, 0\.0163 / 0\.0162 / 0\.2235, 0\.0289 / 0\.0247 / 0\.2920"),

    # --- tagging_optics.py (section 4.5b paragraph, section 7 row) ---
    (["tagging_optics"],
     r"optimum\s+r_h = 46\.5 / 164\.1 / 89\.3 at.*?L/L_HA = 1/6\.8 / 1/12\.8 / 1/9\.5"),
    (["tagging_optics"],
     r"N_tag/yr = 2\.4×10⁶ / 2\.4×10⁶ / 6\.1×10⁶ \("),
    (["tagging_optics"],
     r"best-super-bin 5σ floors of 1\.7 / 2\.3 / 1\.6% per\s+unit P_zz"),
    (["tagging_optics"],
     r"the shape term in the optics' own window 0\.033–0\.039 per unit\s+P_zz"
     r" = 9\.4 / 8\.3 / 10\.7σ per year and 3\.0 / 5\.5 / 2\.6 years to 5σ"),
    (["tagging_optics"],
     r"the bin holding 0\.193 / 0\.105 / 0\.086 o"),
    (["tagging_optics"],
     r"both planes de-squeezed gives a fifth to a quarter of\s+the yield at 1/25–1/71"),
    (["tagging_optics"],
     r"Horizontal-only optimum β\*_x/β\*_x,HA = 46\.5.*?years to 5σ on a 1% exotic-glue term"),
    (["tagging_optics"],
     r"the banner prints that bin's fraction of the tagged sample, 0\.193 / 0\.105 / 0\.086\),"
     r" IR-8 at L_HA worth 4× / 10× / 6×"),
    (["tagging_optics"],
     r"The banner's aperture line reads 2\.50/6\.49.*?"
     r"ε, 9\.36×10⁻¹⁰ / 1\.97×10⁻¹⁹ / 1\.23×10⁻⁵"),

    # --- money_cos2phi.py ---
    (["money_cos2phi"],
     r"\(0\.028, 1\.14\), \(0\.011, 1\.14\), \(0\.071, 3\.13\), \(0\.141, 14\.3\);"
     r" A = 7\.4 / 4\.4 / 9\.5 / 9\.5 ×10⁻\n?³, δA = 1\.7 / 1\.4 / 2\.7 / 4\.5 ×10⁻⁴"),

    # --- money_cos2phi_reco.py (5R block, reprinted by both scans) ---
    (["money_cos2phi_reco_syst", "money_cos2phi_reco_unfold"],
     r"0\.66 / 0\.63 / 0\.69 / 0\.69 \(D = 0\.92 / 0\.99 / 0\.90 / 0\.96\); δÂ = 1\.2 / 0\.9 /\n?"
     r" ?1\.6 / 2\.9 ×10⁻⁴"),
    (["money_cos2phi_reco_unfold"],
     r"bin-by-bin \(−4\.2, \+8\.0, −5\.6, \+4\.9\)% → folded \(.*?\+0\.5\)%"),

    # --- target_mass_bound.py ---
    (["target_mass_bound"],            # section 4.1's prose on the same four blocks
     r"γ² ≤ M²/\(W²_min − M²\) = 0\.0965 and the measured maximum per", ),
    (["target_mass_bound"],
     r"reaches 1\.38 and 1\.20 mrad at the\s+mid and top configurations and 5\.93 mrad"
     r" at the low one", ),
    (["target_mass_bound"],
     r"at most 1\.54×10⁻⁴ over the twelve sweet spots, 1\.59×10⁻³ over the\s+"
     r"polarized-EMC window and 7\.27×10⁻⁴ over the tagged x bins", ),
    (["target_mass_bound"],
     r"cap 0\.0965; grid maxima 0\.0854 / 0\.0577 / 0\.0258 \(⁶Li\) and 0\.0854 / 0\.057\n?"
     r"7 / 0\.0332 \(⁷Li\)"),
    (["target_mass_bound"],
     r"sweet spots max γ² 0\.00564 and max A_∥ shift 0\.56%.*?"
     r"le azimuth shortcut errs by at most 1\.38 / 1\.20 / 5\.93 mrad"),
    (["target_mass_bound"],
     r"polarized-EMC ΔR carries a target-mass term of 0\.120 / 0\.456 / 0\.731 / 1\.084%"
     r" at x = 0\.0\n?89 / 0\.282 / 0\.447 / 0\.708"),
    (["target_mass_bound"],
     r"against the same weights' δΔR = 0\.0477 / 0\.0509 / 0\.0615 / 0\.1224"),
    (["target_mass_bound"],
     r"the rate-weighted ⟨γ²⟩ over the same wind\n?ow running"
     r" 0\.001291 / 0\.004525 / 0\.007067 / 0\.010190, i\.e\. 0\.0102 in the top bin"
     r" and below 0\.010 in the other three"),
    (["target_mass_bound"],
     r"the twist-3 residual at g₂ = 0 \n?is 1\.55×10⁻⁵ / 3\.43×10⁻⁴ / 7\.76×10⁻⁴"
     r" / 1\.59×10⁻³"),
    (["target_mass_bound"],
     r"tagged-triton overlay ≤ 2\n?\.09% at the published configuration and 5\.04% at"),

    # --- nearbeam_reach_gain.py --n-mc 2000000 ---
    (["nearbeam_reach_gain"],
     r"pots at the silicon \(2\.50 × 6\.49, 1\.51 × 2\.12 and 0\.53 × 0\.92\n?"
     r" ?mrad\): acc 0 and 0 of 7 bins at.*?acc 1\.45×10⁻⁵ with 268 tagged/yr and 0 of 7 at"),
    (["nearbeam_reach_gain"],
     r"231 e\n?xpected tagged recoils in 0\.17–0\.25 GeV²"),
    (["nearbeam_reach_gain"],
     r"`nearbeam_reach_gain\.MIN_TAGGED_PER_BIN` = 1000"),
    (["nearbeam_reach_gain"],
     r"pots following: acc 0\.3643 / 0\.2471 / 0\.3246, N_tag 2\.31 / 2\.38 / 6\.00 ×10⁶/yr,"
     r" 7 of 7 \\\|t\\\| bins, δa_t 0\.0063.*?0\.0153 \("),

    # --- tagged_polarimetry_7li.py --config 0/1/2 ---
    (["tagged_polarimetry_7li_1"],
     r"acc\(RP\) 0\.9678 \(YR HA\) vs 0\.9909 \(tagging\); acc\(any far-fwd\) 0\.\n?"
     r"9690 vs 0\.9921; ⟨P₂⟩ slope −0\.1947 vs −0\.1962 against the analytic −0\.2000"),
    (["tagged_polarimetry_7li_1"],
     r"median δA_∥ 0\.01152 vs 0\.01141"),
    (["tagged_polarimetry_7li_0", "tagged_polarimetry_7li_2"],
     r"the tags are 0\.9617 /\n? ?0\.9728 \(YR HA\) against 0\.9800 / 0\.9919 \(tagging\)"),
    (["tagged_polarimetry_7li_0", "tagged_polarimetry_7li_1", "tagged_polarimetry_7li_2"],
     r"the tagging optics multiplies every ⁷Li error bar by 2\.78 / 3\.81 / 3\.15"),
    (["tagged_polarimetry_7li_0", "tagged_polarimetry_7li_1", "tagged_polarimetry_7li_2"],
     r"Roman-Pot tag 0\.9617 / 0\.9678 / 0\.9728 \(YR HA\), 0\.9800 / 0\.9909 / \n?"
     r"0\.9919 \(tagging\)"),
    (["tagged_polarimetry_7li_0", "tagged_polarimetry_7li_1", "tagged_polarimetry_7li_2"],
     r"`acc\(any far-fwd\)` = 0\.9699 / 0\.9690 / 0\.9751 and 0\.9882 / 0\.9921 / 0\.9942"),

    # --- money_tagged_azz.py --events 400000 ---
    (["money_tagged_azz"],
     r"acc 0\.0247 \(YR HA\) vs 0\.2545 \(tagging"),
    (["money_tagged_azz"],
     r"acc × L 0\.0247 vs 0\.0199"),
    (["money_tagged_azz"],
     r"median accepted k 0\.323 vs 0\.177 GeV/c and frac\(k < 0\.15\) 0\.000 vs 0\.363"),
    (["money_tagged_azz"],
     r"A_zz = \+0\.491 \(acceptance-weighted truth \+0\.455\) and −0\.066 \(−0\.095\);"
     r" the θ_k = 90° curve says −0\.482 at both"),
]


# --- default mode: the manual's command lines still parse ------------------

CMD_RE = re.compile(r"python3\s+(?:\S*/)?(scripts/[A-Za-z0-9_]+\.py)([^\n]*)")
FLAG_RE = re.compile(r"(?<![\w-])--[a-z][a-z0-9-]*")
CELL_RE = re.compile(r"(?:python3\s+)?(?:\S*/)?(scripts/[A-Za-z0-9_]+\.py)(.*)$")


def _manual_command_lines():
    """(line number, script path as written, [flags]) for every command the
    manual publishes -- the fenced bash blocks (backslash continuations
    joined) and the command column of the section 7 table."""
    out, lines, inblock = [], MANUAL_TEXT.splitlines(), False
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("```"):
            inblock = line.startswith("```bash")
            i += 1
            continue
        start, joined = i + 1, line
        while inblock and joined.rstrip().endswith("\\") and i + 1 < len(lines):
            i += 1
            joined = joined.rstrip()[:-1] + " " + lines[i]
        if inblock:
            m = CMD_RE.search(joined)
            if m:
                out.append((start, m.group(1), FLAG_RE.findall(m.group(2))))
        elif line.startswith("|"):
            for cell in re.findall(r"`([^`]*)`", line):
                m = CELL_RE.match(cell.strip())
                if m:
                    out.append((start, m.group(1), FLAG_RE.findall(m.group(2))))
        i += 1
    return out


def _resolve(rel):
    for base in ("evgen", "fastsim", "."):
        p = ROOT / base / rel
        if p.is_file():
            return base, p
    return None, None


def _flags_in_source(path):
    """Long option strings the script's own add_argument calls spell out."""
    found = set()
    try:
        tree = ast.parse(path.read_text(), filename=str(path))
    except (OSError, SyntaxError):
        return None
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                and node.func.attr == "add_argument":
            for a in node.args:
                if isinstance(a, ast.Constant) and isinstance(a.value, str) \
                        and a.value.startswith("--"):
                    found.add(a.value)
    return found


def _flags_in_help(base, path):
    try:
        out = subprocess.run([sys.executable, str(path), "--help"],
                             cwd=str(ROOT / base), capture_output=True,
                             text=True, timeout=120)
    except (OSError, subprocess.SubprocessError):
        return set()
    return set(FLAG_RE.findall(out.stdout + out.stderr))


@check("derived: every manual command names a script whose argparse takes its flags")
def _():
    bad, cache = [], {}
    for line, rel, flags in _manual_command_lines():
        base, path = _resolve(rel)
        if path is None:
            bad.append("%s:%d — the manual runs %s, which does not exist under "
                       "evgen/ or fastsim/" % (MANUAL.name, line, rel))
            continue
        if not flags:
            continue
        if rel not in cache:
            cache[rel] = _flags_in_source(path)
        known = cache[rel]
        if known is None:
            bad.append("%s:%d — %s does not parse, so its flags cannot be checked"
                       % (MANUAL.name, line, rel))
            continue
        missing = [f for f in flags if f not in known]
        if missing:                               # only now pay for an import
            known = known | _flags_in_help(base, path)
            cache[rel] = known
            for f in [f for f in flags if f not in known]:
                bad.append("%s:%d — the manual runs `%s %s`, but %s/%s's argparse "
                           "does not accept %s" % (MANUAL.name, line, rel, f,
                                                   base, rel, f))
    return bad


# --- full mode: the numbers themselves --------------------------------------

def _run(key, scratch):
    where, argv, wants_outdir = COMMANDS[key]
    cwd = ROOT / where
    out = scratch / key
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, MPLBACKEND="Agg")
    if argv is None:                              # the manual's own snippet
        src = _sec31b_snippet()
        if src is None:
            return None, "section 3.1b no longer carries a `python3 - <<'PY'` block"
        proc = subprocess.run([sys.executable, "-"], input=src, cwd=str(cwd),
                              capture_output=True, text=True, timeout=900, env=env)
    else:
        cmd = [sys.executable] + argv + (["--outdir", str(out)] if wants_outdir else [])
        proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True,
                              text=True, timeout=900, env=env)
    if proc.returncode != 0:
        tail = (proc.stderr or "").strip().splitlines()
        return None, "exited %d: %s" % (proc.returncode, tail[-1] if tail else "no stderr")
    (out / "stdout.txt").write_text(proc.stdout)
    return proc.stdout, None


def _command_line(key):
    where, argv, wants_outdir = COMMANDS[key]
    if argv is None:
        return "cd %s && the section 3.1b snippet" % where
    return "cd %s && python3 %s%s" % (where, " ".join(argv),
                                      " --outdir <scratch>" if wants_outdir else "")


@check("derived: the manual's expected numbers are what the scripts print today")
def _():
    if not checker.FULL:
        return checker.SKIP
    scratch = _scratch()
    scratch.mkdir(parents=True, exist_ok=True)

    needed = sorted({k for keys, _ in CLAIMS for k in keys})
    bad, pools = [], {}
    for key in needed:
        text, err = _run(key, scratch)
        if err:
            bad.append("%s — %s" % (_command_line(key), err))
            pools[key] = None
        else:
            pools[key] = _output_numbers(text)

    for keys, pattern in CLAIMS:
        m = re.search(pattern, MANUAL_TEXT, S)
        if m is None:
            bad.append("%s — no line now matches the claim anchored on %r; the "
                       "manual was rewritten and this check must follow it"
                       % (MANUAL.name, pattern[:60]))
            continue
        line = _line_of(m.start())
        pool = sorted({v for k in keys if pools.get(k) for v in pools[k]})
        if not pool:
            continue                              # the run already failed above
        for tok, val, tol, pct in _manual_numbers(m.group(0)):
            if not _present(val, tol, pool, pct):
                bad.append("%s:%d — the manual states %s for `%s`, which prints "
                           "nothing within %.3g of it today (nearest %s)"
                           % (MANUAL.name, line, tok, _command_line(keys[0]),
                              tol, _nearest(val, pool, pct)))
    return bad
