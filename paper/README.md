# paper/ — the PLB-class letter

This directory is the letter of `plans/07` (WP6), built on 2026-09-16 by T22 and
T27 (the scaffold and the bibliography), T23–T25 (the figures and Table 1) and
T26/T28/T29 (the text, the assumptions section and the cover letter).
`main.tex` and `cover_letter.md` are marked **DRAFT v0, 2026-09-16, for the
authors' review** in their header comments and nowhere in their body text.

## What is here, and what state it is in

| file | state |
|---|---|
| `main.tex` | **written.** elsarticle two-column, seven sections, ~3,700 words of body text against the ~4,200 of plans/07 §7.5, a 119-word abstract, four figures, Table 1 and 35 references. Every body line that quotes a number ends with a `% src:` comment naming the report section, plan section or build artifact it comes from; a displayed equation carries the same comment on the line above it. |
| `cover_letter.md` | **written.** States the two verified literature firsts of plans/07 §7.8 and what makes each one verified. Editor name, address and date are placeholders for the submitting author. |
| `refs.bib` | **complete and machine-verified.** The 35 records the letter cites, and no others: 31 fetched from INSPIRE and keyed by INSPIRE texkey, plus 4 hand-written `@misc` entries that INSPIRE has no record for, each flagged `HAND-WRITTEN` in an `annote` field (which the bibliography style does not print; the `note` field beside it is what the reference list shows). The correspondence with `main.tex`'s `\cite` commands is exact in both directions and is checked (2 below). Five entries were edited after the fetch for what the style printed; `refs_provenance.md` lists them. |
| `refs_reserve.bib` | the other 40 records of the same fetch: the entries of Report 1's and the coherent note's reference lists the letter does not cite (38 INSPIRE, 2 hand-written). Together the two files cover both source lists in full (4 below). `main.tex` compiles this one only under `\draftbibtrue`; to cite one of these works, move its block into `refs.bib`. |
| `refs_provenance.md` | the audit trail of both bibliography files: the query behind every entry, whether it resolved, the four second queries and what they caught, and the discrepancies found between the source reference lists and the INSPIRE records. |
| `figstyle.py`, `_sources.py`, `fig1..4_*.py`, `table1.py` | the figure drivers and the table generator. They import the producing scripts of `evgen/scripts/` and run at those scripts' published defaults; they never write into `evgen/` or `fastsim/`. |
| `figs/` | the four letter figures (`figN.pdf`, `figN.png`), the numbers behind each (`figN.json`), the generated caption of each (`figN_caption.tex`), `table1.json`, and the register `README.md` that ties each letter figure to the published figure it condenses. |
| `table1.tex` | generated whole by `table1.py` from `figs/fig2.json` and `figs/fig4.json`. Do not edit. |
| `build.sh` | one command for all of it: the four figures, their captions, Table 1 and the PDF. |

## How to build

There is no `pdflatex`, `xelatex`, `lualatex` or `latexmk` on the machines this
programme runs on. The engine is **tectonic**, which carries its own package
bundle and fetches `elsarticle` (and everything else `main.tex` asks for) on
first use, then caches it under `~/.cache/Tectonic`.

```bash
# get a binary once -- NOT into this repository
curl -L -o tectonic.tar.gz \
  https://github.com/tectonic-typesetting/tectonic/releases/download/tectonic%400.17.0/tectonic-0.17.0-x86_64-unknown-linux-musl.tar.gz
tar xzf tectonic.tar.gz

bash paper/build.sh --figures                 # figures, captions and Table 1 only
TECTONIC=$PWD/tectonic bash paper/build.sh    # the above, then paper/main.pdf
```

`build.sh` takes the engine from `$TECTONIC` if that is set and otherwise looks
for `tectonic` on `$PATH`; with neither it exits 1 and says so. The first build
downloads a few hundred files from the tectonic bundle and takes a couple of
minutes; later builds are seconds. The figure half is deterministic — two runs
on an unchanged tree give identical bytes, PDFs included.

Output: `paper/main.pdf`, with `main.aux`, `main.bbl`, `main.blg`, `main.log`
and `main.out` beside it. All six are build artifacts and are git-ignored.

## Which layout, and how long the letter is

`main.tex` carries both elsarticle two-column class lines, one active and one
commented, and the text is identical between them:

| class line | pages | float placement |
|---|---|---|
| `[final,5p,times,twocolumn]` (**active**) | **8** (four lines of the reference list on the last page) | figures and Table 1 beside the text they belong to |
| `[preprint,twocolumn,12pt]` | 16 | every double-column float deferred to the end, four pages of nothing but floats |

The journal layout is the default because it is the one a letter's length can be
judged in. Eight pages (the last carrying four lines of references) is at the
top of the six-to-eight published pages plans/07 §7.2 expects of a PLB letter,
and above the four the WP6 brief for the text asked for; roughly two and a half
of them are the four figures, Table 1's caption and the reference list, so the
text itself is not what would have to shrink. Whether to cut a figure, shorten
the generated captions, or accept eight pages is an author call.

Measured on the active layout (2026-09-16, after the review pass and the
number-guard pass): 0 TeX errors,
0 undefined citations or references, 2 BibTeX warnings (both "empty journal", for the two arXiv-only
records the letter cites: `Maxwell:2018gci`, `Sick:2015spa`), and one 1.9 pt
overfull box in the output routine, which is Table 1's tabular against the text
width and is generated rather than written.

## What is checked mechanically

`tools/checks/paper_numbers.py` runs inside `python3 tools/consistency_check.py`
and, standalone, prints the `0 unresolved` line that plans/07 §7.8 asks for:

```bash
python3 tools/checks/paper_numbers.py
```

Nine checks, none of which touches the network:

1. the scaffold has the files WP6 names;
2. `refs.bib` and the letter's `\cite` commands are an exact correspondence:
   every `\cite` key resolves there, **every entry is reached by a `\cite`**,
   no key is duplicated or sits in both bibliography files, and the letter
   stays inside its reference budget. A record the letter stops citing is not
   left in the file it ships: it moves to `refs_reserve.bib`;
3. every entry is either an INSPIRE-shaped record or flagged `HAND-WRITTEN`
   (in its `annote` field) and accounted for in `refs_provenance.md`;
4. **every numbered entry of Report 1's and the coherent note's reference lists
   is represented in `refs.bib` or `refs_reserve.bib`** — the mapping back is
   carried in the two files themselves, as a `% cited by R1[n]/N[m]` comment
   above every entry, so adding a reference to either source document without
   adding it here fails the sweep;
5. the four producing scripts named in `figs/README.md` exist;
6. the four figures and Table 1 are built, `table1.tex` matches cell for cell
   the cells the driver wrote into `figs/table1.json`, and those cells' raw
   values match `figs/fig2.json` and `figs/fig4.json`;
7. **every number in Table 1, in its caption and in the abstract is a number
   Report 1's template or the manual's expected-numbers table (§7) already
   prints** — each is declared in the `SOURCED` table of the checker against
   the string that document prints it in, the string must still be there
   verbatim, and the letter's value must agree with a number inside it to the
   precision that number is printed at. Five numbers are in neither document
   and are declared against the driver that measured them, with the reason:
   the four inclusive per-bin event counts (Report 1 quotes the programme
   totals) and the coherent row's $A/\delta\hat A$. Adding a number to the
   abstract or to Table 1 without declaring it fails the sweep, and so does
   editing either source document so that the string is gone;
8. `main.tex` carries the programme's author line, affiliation and
   writing-assistance line, and an elsarticle two-column class line;
9. **every number in the letter's text names where it came from** — see below —
   the abstract is inside its 120-word budget, `main.tex`'s header carries the
   draft marker while its body does not use the word, `cover_letter.md` states
   the two verified firsts, and no form retired in
   `tools/retired_numbers.json` reaches either file.

The last of those is the T26 leg, and it is the weaker of the two number
checks: it asks that a number be attributed, not that it be one already
published. `tools/checks/retired_strings.py` owns the
retired-forms list but its corpus is the reports, the plans, the docs and the
module docstrings of `evgen/` and `fastsim/`; `paper/` is outside it, and the
letter is the one document in this repository that will leave it.

### The sourcing convention

A "quoted quantity" is a decimal, a number of two digits or more, or a
percentage; a bare single digit is notation (`spin-1`, `IP6`, `cos 2phi`,
`A > 2`, `Section 3`) and carries nothing to source. Every body line that quotes
one ends with

```latex
... errors of $1.7$, $1.4$, $2.7$ and $4.5 \times 10^{-4}$ --- % src: R1 sec 5.1
```

where the comment names a Report section (`R0`–`R4`), a plan section
(`plans/07 sec 7.4`), the coherent note, a manual row, or the JSON a figure
driver wrote. Whether the source says what the sentence says is a reading and
belongs to the author; that a number has a source at all is a grep, and it is
this check.

## The footing the letter quotes

Single-fill estimator, toy structure functions with a placeholder `R`,
generator level — the footing Report 1 publishes. Section 6 of the letter bounds
the two alternatives rather than adopting either: the R1998 swap (+18, +18, +8,
−4% on Δ/F₁ at the four bins) and the nuclear-grid backend (rates ×0.77, tagged
coherent yield ×0.70, significances 18–28σ), the latter with the statement that
EPPS21nlo_CT18Anlo_Li6 has no support below Q² = 1.69 GeV² and that nothing is
claimed on that footing below that scale. The two-fill estimator of plans/07
§7.1 is not quoted anywhere in the letter.

## What is decided, and what is not

Decided and written into `main.tex`:

* **D1** (scope) — plans/07 §7.0's recommended inclusive + coherent gluonometry
  letter; the seven sections are that letter's.
* **D3** (authorship) — lead C. Peng, second author J. Zhou, Physics Division,
  Argonne National Laboratory. The co-authors §7.7 lists as *candidates to
  invite* are not in the file.
* **D5** (title) — candidate (a), the one §7.5 recommends. Candidates (b) and
  (c) are recorded there, not here.

Not decided:

* **D2** (generator-level or reconstructed-level main figures). All four letter
  figures are generator-level, because the four published money plots they
  condense are; the reconstructed-level results appear in the text of Section 4
  rather than in a figure. If D2 goes the other way the drivers point at the
  reconstructed producers and the structure does not change.
* The single-fill / two-fill footing of §7.1, and the toy / grid column of the
  WP1 addendum. Both are bounded in Section 6 and neither is adopted.
* **D4** (circulation), the funding line and the people thanked in the
  acknowledgments, and the editor, address and date of the cover letter.
* **D6 and the submission itself.** The text is a draft v0: the journal
  decision, the submission, the arXiv posting and the arXiv number, the
  Zenodo DOIs the four companion-report entries need in place of "companion
  report, in preparation", and the re-run of the forward-citation sweep behind
  the cover letter's second first (last run 2026-08-10) are all the authors'.
* The **three "do not say" lists** of §7.8. They are named in plans/07 and
  enumerated nowhere in this repository. What is enforced in their place is
  `tools/retired_numbers.json`, over `main.tex` and `cover_letter.md`, by
  check 8 above.

## Edit record

| date | what |
|---|---|
| 2026-09-16 | written: the scaffold and the bibliography (T22/T27), the four figures and Table 1 (T23–T25), the text, the assumptions section and the cover letter (T26/T28/T29). |
| 2026-09-16 | review pass (R-PAPER): the tensor-sign convention of Eq. (1), the bag normalization of Eq. (4), the null-test footing, the reference-tag labelling, the 3×2 clipped bins, six reference-list defects, and thirteen smaller corrections. Its open items are in that review. |
| 2026-09-16 | number-guard pass: `refs.bib` split into the letter's 35 records and `refs_reserve.bib`'s 40, `main.tex`'s bibliography switch made to compile the reserve only in draft mode, one `\cite{Nikolaev:1996jy}` added to §5 (no change to the count), and the Table 1 / caption / abstract numbers declared against Report 1 and the manual in `tools/checks/paper_numbers.py` (check 7). |
