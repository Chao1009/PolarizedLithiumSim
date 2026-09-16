# paper/figs — the four letter figures

This directory holds the four letter figures (`figN.pdf`, `figN.png`), the
numbers behind each (`figN.json`), the generated caption of each
(`figN_caption.tex`) and `table1.json`, all written by `bash paper/build.sh`
(T23–T25 of plans/07 WP6). The register below says for each letter figure which
*published* figure it condenses and which script produces that published
figure; `tools/checks/paper_numbers.py` re-reads the table and asserts that
every producing script it names exists and that each `figN.json` names the
same script and published figure.

## The register

| letter figure | condenses | producing script | published output | what the letter figure must show (plans/07 WP6) |
|---|---|---|---|---|
| Fig. 1 | phase space | `evgen/scripts/phase_space_bins.py` | `evgen/phase_space_bins_6Li.png` | one inclusive panel with the analysis bins, the coherent support as a contour or an inset |
| Fig. 2 | money plot 5 | `evgen/scripts/money_cos2phi.py` | `evgen/money_cos2phi_6Li.png` | two $\phi'$ panels plus the amplitude versus $x$ |
| Fig. 3 | money plot 7 | `evgen/scripts/money_delta_extraction.py` | `evgen/money_delta_extracted_6Li.png` | two $Q^2$ slices of the extracted $x\Delta(x,Q^2)$ |
| Fig. 4 | money plot 6 | `evgen/scripts/money_cos2phi_coherent.py` | `evgen/money_cos2phi_coherent_6Li.png` | the $a_2$ anchor and band, the tagged $\phi'$ distribution, an acceptance inset |

The four published outputs are the stems registered in `reports/build_report.py`
as `__PS__`, `__M5__`, `__M7__` and `__M6__`; the reproduction commands are in
`docs/reproduction_manual.md` §4.1–4.2. `tools/checks/paper_numbers.py` re-reads
the table above and asserts that every producing script it names exists.

## The rule the drivers must respect

The published figure stems are bit-for-bit reproducible and are **not** to be
changed to suit the letter. A letter figure is a *new* artifact drawn by a new
driver under `paper/figs/`; it never overwrites `evgen/*.png`, and any
non-default setting it needs must reach the producing script through that
script's own stem-tag flag (`money_cos2phi.py --output-stem-tag`,
`fom.run_share_tag`, `money_tagged_azz.py --output-stem`), so that a tagged
output can never be mistaken for a published one.

## Table 1

`paper/table1.tex` is generated whole by `paper/table1.py` from `fig2.json` and
`fig4.json`, so a cell of the table and the number beside it in Figure 2 or 4
are one draw of one pseudo-experiment; `table1.json` is the cell-by-cell record
the sweep checks the typeset table against.
