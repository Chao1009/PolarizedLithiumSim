#!/usr/bin/env python3
"""Build `paper/table1.tex` -- the letter's projections table.

plans/07 SS 7.5 gives Table 1 to Section 4 and plans/07 SS 7.8 requires it to
come out of the same build script as the figures (the run-19 ledger's T25).
That is what this module is: it reads `figs/fig2.json` and `figs/fig4.json`,
which the figure drivers wrote from their own pseudo-experiments, and formats
them.  It measures nothing, so a number in Table 1 and the same number in
Figure 2 or Figure 4 cannot disagree -- they are one draw of one
pseudo-experiment, at one seed.

The columns are those Report 1 publishes: SS 5.1 gives the four sweet spots
their injected amplitude, their one-year and ten-year statistical errors and
the resulting significance; SS 6.2 gives the coherent best super-bin its yield,
its errors and its 5 sigma floors, and SS 6.3 the deformation-anchored
<a_2>.  Everything is quoted per unit P_zz, which is the footing the two
sections share once SS 6.3's "0.036 at P_zz = 0.6, i.e. 0.060 per unit P_zz"
is read; the caption says so, and says which of the two the figure draws.

`figs/table1.json` is written beside the .tex with the same numbers in
machine form, so that `tools/checks/paper_numbers.py` can assert that the
typeset cells are the measured values and not a transcription.

Usage:  python3 paper/table1.py                 (writes table1.tex)
"""

import pathlib
import re

import _sources as src
from _sources import beam_math, dec, sci, sig3

OUT = src.PAPER / "table1.tex"

# LaTeX label of each figure the table points at, as main.tex defines them.
FIG_INCL = "fig:inclusive"
FIG_COH = "fig:coherent"


def _model_constant(info):
    """The bag moment the Delta model is normalized to, out of its own
    `info()` string (`moment_A(variant=..., c=-0.012, ...)`)."""
    m = re.search(r"\bc=(-?\d+(?:\.\d+)?)", info)
    return None if m is None else float(m.group(1))


def build():
    f2 = src.load("fig2")
    f4 = src.load("fig4")

    rows = []
    for s in f2["spots"]:
        rows.append({
            "kind": "inclusive", "spot": s["spot"],
            "bin": r"$(%.3g,\ %.3g)$" % (s["x"], s["q2"]),
            "n": s["n_1yr"], "amp": s["a_truth"],
            "err1": s["err_1yr"], "err10": s["err_10yr"],
            "significance": s["significance_1yr"],
        })
    b = f4["best_bin"]
    rows.append({
        "kind": "coherent", "spot": None,
        "bin": (r"$(%.2f\text{--}%.2f)\times10^{-3},\ %.2f\text{--}%.2f$"
                % (1e3 * b["x_lo"], 1e3 * b["x_hi"],
                   b["q2_lo"], b["q2_hi"])),
        "n": b["n_1yr"], "amp": f4["a2_per_pzz"],
        "err1": f4["err_1yr"], "err10": f4["err_10yr"],
        "significance": abs(f4["a2_per_pzz"] / f4["err_1yr"]),
    })

    cells = []
    for r in rows:
        cells.append({
            "kind": r["kind"], "spot": r["spot"], "bin": r["bin"],
            "n": "$%s$" % sci(r["n"]),
            "amp": "$%s$" % sig3(1e3 * r["amp"]),
            "err1": "$%s$" % sig3(1e3 * r["err1"]),
            "err10": "$%s$" % sig3(1e3 * r["err10"]),
            "significance": "$%.0f$" % r["significance"],
            "values": {"n_1yr": r["n"], "amp": r["amp"],
                       "err_1yr": r["err1"], "err_10yr": r["err10"],
                       "significance_1yr": r["significance"]},
        })
    return f2, f4, cells


def caption(f2, f4):
    c = _model_constant(f2["delta_model"])
    dil = float(f2["settings"]["dilution"])
    lumi1 = float(f2["settings"]["lumi_1yr"])
    lumi10 = float(f2["settings"]["lumi_10yr"])
    return (
        r"Projected statistical precision for the double-helicity-flip "
        r"observable at %s with $P_{zz} = %.2f$, in one-year and ten-year "
        r"EIC programmes (%g and %g\,fb$^{-1}$ per nucleon). "
        r"$A$ is the injected $\cos 2\phi'$ amplitude per unit $P_{zz}$ and "
        r"$\delta\hat A$ the statistical error on the fitted amplitude on "
        r"the same footing. "
        r"The four inclusive rows are the analysis super-bins of Fig.~\ref{%s}, "
        r"each of $3\times3$ grid cells ($3\times2$ for the two clipped at the "
        r"grid\'s $Q^2 = 1$\,GeV$^2$ edge), with $\Delta$ the moment-constrained "
        r"ansatz normalized to $\int x\Delta\,\mathrm{d}x = %.3f\,\alpha_s$ at "
        r"$\langle Q^2\rangle = %.2f$\,GeV$^2$ and carrying the $^{6}$Li "
        r"per-nucleon dilution $1/%.0f$. "
        r"The coherent row is the tagged-count-maximal super-bin of the "
        r"intact-recoil channel behind a near-beam $p_T > %.2f$\,GeV tag, "
        r"which keeps $%.1f\%%$ of the coherent recoils for $f_0 = %.2f$ and "
        r"$B = %g$\,GeV$^{-2}$, a reference envelope that no published IP6 optics "
        r"delivers (Sec.~\ref{sec:coherent}); its amplitude is the deformation-anchored "
        r"$\langle a_2\rangle_{\mathrm{tag}} = %.3f$ per unit $P_{zz}$, "
        r"i.e.\ $%.3f$ at $P_{zz} = %.2f$ over a band $%s$--$%s$, and "
        r"its $5\sigma$ detection floors are $%s$ (one year) and $%s$ "
        r"(ten years) per unit $P_{zz}$ (Fig.~\ref{%s}). "
        r"Errors are statistical only; backgrounds, purity, polarimetry and "
        r"tensor radiative corrections are not included. "
        r"Every entry is produced by the drivers that build "
        r"Figs.~\ref{fig:phasespace}--\ref{%s}, from the same seeds as those "
        r"figures."
        % (beam_math(f2["beam"]), f2["pzz"], lumi1, lumi10,
           FIG_INCL, c, f2["q2_ref"], 1.0 / dil,
           f4["pt_cut"], 100 * f4["tag_acceptance"], f4["f0"], f4["slope_b"],
           f4["a2_per_pzz"], f4["a2_tagged"], f4["pzz"],
           dec(min(f4["a2_tagged_band"]), 3),
           dec(max(f4["a2_tagged_band"]), 3),
           "%s\\times10^{-3}" % sig3(1e3 * f4["floor5_1yr"]),
           "%s\\times10^{-3}" % sig3(1e3 * f4["floor5_10yr"]),
           FIG_COH, FIG_COH))


def render(f2, f4, cells):
    head = [
        r"%% paper/table1.tex -- GENERATED by paper/table1.py; do not edit.",
        r"%% Every number is read out of paper/figs/fig2.json and fig4.json,",
        r"%% which the figure drivers wrote from their own pseudo-experiments.",
        r"\begin{table*}[t]",
        r"  \centering",
        r"  \caption{%s}" % caption(f2, f4),
        r"  \label{tab:projections}",
        r"  \begin{tabular}{@{}l r r r r r@{}}",
        r"    \toprule",
        (r"    bin $(x,\ Q^{2}\,[\mathrm{GeV}^{2}])$ & $N$ (1\,yr) & "
         r"$A$ $[10^{-3}]$ & $\delta\hat A$ (1\,yr) $[10^{-3}]$ & "
         r"$\delta\hat A$ (10\,yr) $[10^{-3}]$ & $A/\delta\hat A$ \\"),
        r"    \midrule",
        (r"    \multicolumn{6}{@{}l}{\emph{inclusive} $\cos2\phi'$, "
         r"scattered electron only} \\"),
    ]
    body = []
    for c in cells:
        if c["kind"] == "coherent":
            body.append(r"    \midrule")
            body.append(
                r"    \multicolumn{6}{@{}l}{\emph{coherent} "
                r"$e\,^{6}\mathrm{Li}\to e'X\,^{6}\mathrm{Li}$(g.s.), "
                r"intact recoil tagged} \\")
        body.append("    %s & %s & %s & %s & %s & %s \\\\"
                    % (c["bin"], c["n"], c["amp"], c["err1"], c["err10"],
                       c["significance"]))
    tail = [
        r"    \bottomrule",
        r"  \end{tabular}",
        r"\end{table*}",
    ]
    return "\n".join(head + body + tail) + "\n"


def main():
    f2, f4, cells = build()
    OUT.write_text(render(f2, f4, cells), encoding="utf-8")
    print("wrote", OUT)

    payload = {
        "table": "table1",
        "built_from": ["figs/fig2.json", "figs/fig4.json"],
        "beam": f2["beam"], "pzz": f2["pzz"],
        "lumi_1yr": f2["settings"]["lumi_1yr"],
        "lumi_10yr": f2["settings"]["lumi_10yr"],
        "coherent": {
            "pt_cut": f4["pt_cut"], "tag_acceptance": f4["tag_acceptance"],
            "f0": f4["f0"], "slope_b": f4["slope_b"],
            "a2_per_pzz": f4["a2_per_pzz"], "a2_tagged": f4["a2_tagged"],
            "a2_tagged_band": f4["a2_tagged_band"],
            "floor5_1yr": f4["floor5_1yr"], "floor5_10yr": f4["floor5_10yr"],
        },
        "rows": cells,
    }
    src.dump("table1", payload)

    print("%-30s %12s %8s %8s %9s %6s"
          % ("bin", "N (1 yr)", "A[1e-3]", "dA 1yr", "dA 10yr", "sig"))
    for c in cells:
        v = c["values"]
        label = (c["bin"].replace("$", "").replace(r"\text{--}", "-")
                 .replace(r"\times", "x").replace(r"\ ", " "))
        print("%-30s %12.3e %8.3f %8.3f %9.4f %6.0f"
              % (label, v["n_1yr"], 1e3 * v["amp"], 1e3 * v["err_1yr"],
                 1e3 * v["err_10yr"], v["significance_1yr"]))


if __name__ == "__main__":
    main()
