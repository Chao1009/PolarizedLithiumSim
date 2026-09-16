# paper/refs_provenance.md — where every bibliography entry came from

*Written 2026-09-16 by the T27 pass of plans/07 WP6 (development run 19, W4).*

`paper/refs.bib` is machine-fetched. This file is its audit trail: one row per
work, the exact INSPIRE query, whether it resolved, and the texkey that came
back. Nothing in either bibliography file was typed by hand except the six
flagged `@misc` entries at the bottom of each.

The fetch produced 75 records, and they are held in two files. `refs.bib` is
the letter's own bibliography: the 35 works `main.tex` cites, and nothing else,
so that the file that goes to the journal has no record a reader cannot find in
the reference list. `refs_reserve.bib` holds the other 40 — the entries of the
two source reference lists the letter does not cite — in the same format, with
the same `% cited by` mapping comments, so that the coverage below is still the
coverage of both source lists. The split was made on 2026-09-16, after the
fetch and after the review pass, by moving whole blocks: every one of the 75
records is byte for byte the record the fetcher wrote.

## Scope

The list of works was built by walking, entry by entry, the two reference lists
plans/07 WP6 names as the sources of the letter's text:

* `reports/cos2phi_money_plots_report.template.html` — Report 1, entries `[1]`–`[43]`;
* `docs/note_cos2phi_coherent_6Li.md` — entries `[1]`–`[50]`.

The two lists overlap heavily (the `cited by` column below records both labels
where they do), and several entries of each bundle two or three distinct works
under one number — Report 1 `[24]` is Lappi *and* Kowalski *et al.*, the note's
`[30]` is three Cosyn papers. A bundled entry contributes one BibTeX record per
work, which is why 43 + 50 report entries reduce to **69 unique INSPIRE
records + 6 hand-written `@misc`** rather than to 93.

## Method

```
https://inspirehep.net/api/literature?q=<query>&format=bibtex
```

fetched with `curl --max-time 60`, the query URL-quoted. Query forms used, in
order of preference: `arxiv:<id>` where the source list gives an arXiv number;
`doi:<doi>` where it gives a journal reference whose DOI is unambiguous;
`t "<title>"` for the pre-arXiv literature (1967–1992) and for reports. The
**first** record of a multi-hit response is taken, and the hit count is recorded
in the notes column so that a reviewer can see which rows were ambiguous. The
INSPIRE texkey is kept verbatim as the citation key.

## Results

**All 69 works in the fetch list resolved on INSPIRE.** Four of them needed a
second query (recorded in the notes column): three title queries returned
nothing and one returned the wrong paper. Two further works — the ePIC seminar
of Report 1 `[17]` and the 2021 IAEA quadrupole report of Report 1 `[39]` — were
queried and came back empty; together with the four internal reports of this
series they are the six hand-written `@misc` entries listed near the end of this
file.

| texkey | cited by | query that produced it | status | notes |
|---|---|---|---|---|
| `Jaffe:1989xy` | R1[1]/N[1] | `t "Nuclear Gluonometry"` | resolved |  |
| `Hoodbhoy:1988am` | R1[2]/N[2] | `t "Novel Effects in Deep Inelastic Scattering from Spin 1 Hadrons"` | resolved |  |
| `Maxwell:2018gci` | R1[3]/N[3] | `arxiv:1803.11206` | resolved |  |
| `Sather:1990bq` | R1[4]/N[4] | `t "Size and scaling of the double helicity flip hadronic structure function"` | resolved |  |
| `Kumano:2019igu` | R1[5]/N[5] | `arxiv:1910.12523` | resolved |  |
| `Kumano:2020gfk` | R1[5]/N[28] | `arxiv:2003.06623` | resolved |  |
| `Detmold:2016gpy` | R1[6]/N[25] | `arxiv:1606.04505` | resolved |  |
| `Detmold:2017oqb` | R1[6]/N[26] | `arxiv:1703.08220` | resolved |  |
| `Winter:2017bfs` | R1[7]/N[27] | `arxiv:1709.00395` | resolved |  |
| `EPIOSScientificConsortium:2025dfc` | R1[8]/N[50] | `arxiv:2510.10794` | resolved |  |
| `Mantysaari:2024xmy` | R1[10]/N[13] | `arxiv:2408.13213` | resolved |  |
| `Cosyn:2024drt` | R1[11] | `arxiv:2410.12764` | resolved |  |
| `Nzar:1992ax` | R1[12]/N[34] | `t "Estimation of the double helicity flip deuteron structure function"` | resolved |  |
| `AbdulKhalek:2021gbh` | R1[13]/N[7] | `arxiv:2103.05419` | resolved |  |
| `E143:1998nvx` | R1[14] | `arxiv:hep-ex/9808028` | resolved |  |
| `Bacchetta:2006tn` | R1[15] | `arxiv:hep-ph/0611265` | resolved |  |
| `ATHENA:2022hxb` | R1[17] | `arxiv:2210.09048` | resolved |  |
| `Arratia:2021tsq` | R1[17] | `arxiv:2110.05505` | resolved |  |
| `Jentsch:2021qdp` | R1[18] | `arxiv:2108.08314` | resolved |  |
| `Burkert:2022hjz` | R1[19]/N[8] | `arxiv:2211.15746` | resolved |  |
| `Kim:2026ytt` | R1[19]/N[8] | `arxiv:2602.04636` | resolved |  |
| `Pitt:2024utg` | R1[20]/N[9] | `arxiv:2409.02811` | resolved |  |
| `Aschenauer:2025cdq` | R1[20]/N[9] | `arxiv:2503.05908` | resolved |  |
| `Boer:2025ixc` | N[9] | `arxiv:2512.15064` | resolved |  |
| `Chang:2025pgi` | R1[22]/N[10] | `arxiv:2511.05638` | resolved |  |
| `H1:2012pbl` | R1[23]/N[44] | `arxiv:1203.4495` | resolved |  |
| `Lappi:2009wz` | R1[24]/N[45] | `arxiv:0907.4588` | resolved |  |
| `Kowalski:2008sa` | R1[24]/N[45] | `arxiv:0805.4071` | resolved |  |
| `Sick:2015spa` | R1[25]/N[47] | `arxiv:1505.06924` | resolved |  |
| `STAR:2021wwq` | R1[26]/N[42] | `arxiv:2109.07625` | resolved |  |
| `Suelzle:1967zz` | R1[27]/N[11] | `t "Elastic Electron Scattering from Li-6 and Li-7"` | resolved | 2 hits, first taken |
| `Li:1971tk` | R1[28]/N[12] | `t "High-energy electron scattering from Li-6"` | resolved |  |
| `Mantysaari:2026dps` | R1[29]/N[14] | `arxiv:2605.00454` | resolved |  |
| `Mantysaari:2024qmt` | R1[30]/N[15] | `arxiv:2411.14934` | resolved |  |
| `Veal:1999yk` | R1[31]/N[16] | `doi:10.1103/PhysRevC.60.064003` | resolved | first query 't "Determination of the asymptotic D-state to S-state ratio for Li-6"' returned 0 hit(s); re-queried by DOI |
| `Hupin:2014iqa` | R1[31] | `doi:10.1103/PhysRevLett.114.212502` | resolved | first query 't "Ab initio predictions for polarized deuterium-tritium thermonuclear fusion"' returned 1 hit(s); re-queried by DOI |
| `Wiringa:1998hr` | R1[32]/N[17] | `arxiv:nucl-th/9807037` | resolved |  |
| `Bouwhuis:1998jj` | R1[33]/N[19] | `arxiv:nucl-ex/9810004` | resolved |  |
| `JLABt20:2000uor` | R1[33]/N[19] | `arxiv:nucl-ex/0001006` | resolved |  |
| `Nikolenko:2003zq` | R1[33]/N[19] | `doi:10.1103/PhysRevLett.90.072501` | resolved | first query 't "Measurement of the tensor analyzing power T20 in elastic electron deuteron scattering"' returned 0 hit(s); re-queried by DOI |
| `HERMES:2005pon` | R1[34]/N[20] | `arxiv:hep-ex/0506018` | resolved |  |
| `Nikolaev:1996jy` | R1[35]/N[21] | `arxiv:hep-ph/9611460` | resolved |  |
| `Bora:1997pi` | R1[36]/N[22] | `arxiv:hep-ph/9711323` | resolved |  |
| `Edelmann:1997ik` | R1[37]/N[23] | `arxiv:hep-ph/9709455` | resolved |  |
| `Hatta:2016dxp` | R1[38]/N[6] | `t "Probing the Small-x Gluon Tomography in Correlated Hard Diffractive Dijet Production"` | resolved |  |
| `Stone:2016bmk` | R1[39]/N[35] | `t "Table of nuclear electric quadrupole moments"` | resolved |  |
| `Tilley:2002vg` | R1[40]/N[39] | `doi:10.1016/S0375-9474(02)00597-3` | resolved | first query 't "Energy levels of light nuclei A=5, 6, 7"' returned 0 hit(s); re-queried by DOI |
| `LEPS:2017nqz` | R1[41]/N[41] | `arxiv:1711.01095` | resolved |  |
| `Chang:2021jnu` | R1[42]/N[43] | `arxiv:2108.01694` | resolved |  |
| `ZEUS:2008xhs` | R1[43] | `arxiv:0812.2003` | resolved |  |
| `Garcon:2001sz` | N[18] | `arxiv:nucl-th/0102049` | resolved |  |
| `HERMES:2010hnl` | N[24] | `arxiv:1008.3996` | resolved |  |
| `Zhao:2025vol` | N[29] | `arxiv:2508.06134` | resolved |  |
| `Cosyn:2018rdm` | N[30] | `arxiv:1806.01177` | resolved |  |
| `Cosyn:2018thq` | N[30] | `t "Polynomiality sum rules for generalized parton distributions of spin-1 targets"` | resolved |  |
| `Cosyn:2019eeg` | N[30] | `arxiv:1907.08662` | resolved |  |
| `Boer:2016xqr` | N[31] | `arxiv:1607.01654` | resolved |  |
| `Cotogno:2017puy` | N[31] | `arxiv:1709.07827` | resolved |  |
| `Berger:2001zb` | N[32] | `t "Generalized parton distributions in the deuteron"` | resolved |  |
| `Cano:2003ju` | N[33] | `arxiv:hep-ph/0307231` | resolved |  |
| `Kirchner:2003wt` | N[33] | `arxiv:hep-ph/0302007` | resolved |  |
| `Angeli:2013epw` | N[35] | `t "Table of experimental nuclear ground state charge radii: An update"` | resolved |  |
| `Guan:2024oyp` | N[35] | `arxiv:2403.06384` | resolved |  |
| `Wojtsekhowski:2024lzg` | N[36] | `arxiv:2406.11480` | resolved |  |
| `Dalton:2025jzt` | N[37] | `arxiv:2504.21177` | resolved |  |
| `Dalton:2025fbh` | N[38] | `arxiv:2508.06481` | resolved |  |
| `Accardi:2012qut` | N[46] | `arxiv:1212.1701` | resolved |  |
| `Bylinkin:2022rxd` | N[48] | `arxiv:2208.14575` | resolved |  |
| `Pybus:2026udh` | N[49] | `arxiv:2606.11491` | resolved |  |

## The four second queries, and what they caught

| work | first query | result | second query | what it changed |
|---|---|---|---|---|
| Veal *et al.*, PRC 60 (1999) 064003 | `t "Determination of the asymptotic D-state to S-state ratio for Li-6"` | 0 hits | `doi:10.1103/PhysRevC.60.064003` | INSPIRE spells it "asymptotic D- to S-state ratio"; the title query failed on the hyphenation, not on the record. `Veal:1999yk`. |
| Nikolenko *et al.*, PRL 90 (2003) 072501 | `t "Measurement of the tensor analyzing power T20 in elastic electron deuteron scattering"` | 0 hits | `doi:10.1103/PhysRevLett.90.072501` | INSPIRE's title is "the tensor analyzing power**s** T(20) **and T(21)**". `Nikolenko:2003zq`. |
| TUNL A = 5, 6, 7 evaluation | `t "Energy levels of light nuclei A=5, 6, 7"` | 0 hits | `doi:10.1016/S0375-9474(02)00597-3` | INSPIRE's title is "Energy levels of light nuclei **A=5, A=6, A=7**". `Tilley:2002vg`. |
| Hupin, Quaglioni, Navratil (Report 1 `[31]`, "NCSMC") | `t "Ab initio predictions for polarized deuterium-tritium thermonuclear fusion"` | 1 hit, **wrong paper** (`Hupin:2018biv`, Nature Commun. 10 (2019) 351) | `doi:10.1103/PhysRevLett.114.212502` | The PRL the repository cites is `Hupin:2014iqa`, *Unified description of ⁶Li structure and deuterium-⁴He dynamics with chiral two- and three-nucleon forces*, PRL 114 (2015) 212502. **The repository's citation was right and the title guess was wrong**; recorded here because a title query that returns exactly one plausible-looking hit is the failure mode this audit trail exists to catch. |

## Discrepancies between the source lists and the INSPIRE record

These are differences of wording or of bibliographic detail between what Report 1
or the coherent note prints and what INSPIRE returns. None of them is a wrong
citation; all of them are things the submission pass should decide about, since
`\bibliography` will print the INSPIRE form.

| entry | source list says | INSPIRE record says |
|---|---|---|
| `EPIOSScientificConsortium:2025dfc` | *Realizing the scientific program with polarized ion beams at EIC* (Report 1 `[8]`, note `[50]`, `refs/README.md`) | *Realizing the scientific program with polarized ion beams at **the future BNL** Electron Ion Collider*, and the author field is a collaboration ("EPIOS Scientific Consortium") with "Atoian, Grigor and others" |
| `Guan:2024oyp` | note `[35]`: "Q(⁶Li) = −0.0806(6) fm² … (Pyykkö compilations, via arXiv:2403.06384)" | that arXiv id is Guan *et al.*, *Direct extraction of nuclear structure information using precision lithium-ion spectroscopy*, Phys. Rev. A 112 (2025) L010801 — a lithium-ion spectroscopy paper, not a Pyykkö compilation. The note's "via" is literally true, but the bibliography will print the spectroscopy paper, so the letter must not attribute the compilation to it. |
| `Kim:2026ytt` | Report 1 `[19]`, note `[8]`: "BNL LDRD 23-050 closeout, arXiv:2602.04636" | *A Second EIC Detector: Physics Case and Conceptual Design*, Kim *et al.*, report number "BNL LDRD 23-050", no journal yet (Feb 2026) |
| `Suelzle:1967zz` | Report 1 `[27]`, note `[11]` | the title query returned **2** hits; the second is Bumiller *et al.*, PRC 5 (1972) 391, a different measurement. The first hit is the right one (Phys. Rev. 162 (1967) 992). |
| `Stone:2016bmk` vs `Stone:2021indc` | Report 1 `[39]` bundles "IAEA INDC(NDS)-0833 (2021) **and** At. Data Nucl. Data Tables 111–112 (2016) 1" | INSPIRE has the ADNDT table (`Stone:2016bmk`) but no record at all for the 2021 IAEA report, which is therefore hand-written |

## Hand-written entries

Six, all flagged with a `note = "HAND-WRITTEN: ..."` field inside `refs.bib`:

| key | what it is | query tried | why it is hand-written |
|---|---|---|---|
| `Peng:Report0` | Report 0 of this series | — | internal document |
| `Peng:Report2` | Report 2 of this series | — | internal document |
| `Peng:Report3` | Report 3 of this series | — | internal document |
| `Peng:Report4` | Report 4 of this series | — | internal document |
| `Maple:2024sem` | ePIC seminar, Birmingham, 11 Dec 2024 | `t "Tracking and inclusive DIS reconstruction with the ePIC detector"` | 0 hits — seminar slides |
| `Stone:2021indc` | IAEA INDC(NDS)-0833 (2021) | `r INDC(NDS)-0833`, then `t "Table of Recommended Nuclear Electric Quadrupole Moments"` | 0 hits both times — IAEA nuclear-data report |

The **EIC Yellow Report is not among them**. The T22/T27 brief anticipated that
it might have to be hand-written; `arxiv:2103.05419` resolves cleanly to
`AbdulKhalek:2021gbh` (Nucl. Phys. A 1026 (2022) 122447), so it is a machine
entry like the rest.

## Re-running this

The fetcher is not in the repository — it is a one-off under the run-19
scratchpad (`W4/G-PAPER-1/{works.py,fetch.py,gen_bib.py,write_bib.py}`), with the
raw INSPIRE responses in `inspire_raw.json`. What *is* in the repository and is
meant to be re-run is `tools/checks/paper_numbers.py`, which re-checks the
structural invariants of both files (every `\cite` key defined in `refs.bib`,
every entry of `refs.bib` reached by a `\cite`, no key duplicated or present
in both files, every entry either INSPIRE-shaped or flagged hand-written)
without touching the network. Promoting a reserve record is therefore a move
between the two files and not a copy: leaving it in both fails the sweep.

## Not done here, and why

**Volume and page re-verification at submission** (plans/07 §7.3 row 10) is
deliberately *not* done: INSPIRE's own record is taken as the verification, and
several entries are 2025–2026 papers whose journal data are still settling
(`Kim:2026ytt` and `Boer:2025ixc` have no journal field at all yet). The
re-verification is a submission-time pass, not a scaffolding-time one.

## Edits after the fetch (2026-09-16, review pass)

Five entries were edited by hand after the INSPIRE fetch, each for what the
bibliography style printed rather than for what the record says; the INSPIRE
texkeys are unchanged and no query above is affected.

| texkey | edit | why |
|---|---|---|
| `Jaffe:1989xy` | `title` recased from `{NUCLEAR GLUONOMETRY}` to `{Nuclear gluonometry}` | the INSPIRE record carries the 1989 all-capitals title, and `elsarticle-num` prints the braced title verbatim |
| `Maxwell:2018gci`, `Sick:2015spa` | `month` field removed | with no `journal`, `elsarticle-num` printed the numeric month as "(3 2018)" and "(5 2015)" |
| `Nikolaev:1996jy` | second author `Schafer` → `Sch{\"a}fer` | the umlaut the record dropped |
| the six hand-written `@misc` | the `HAND-WRITTEN: …` flag moved from `note` to `annote`; `note` now carries the text the reference list prints ("companion report, in preparation" for Reports 0/2/3/4, "seminar slides", "IAEA nuclear-data report") | `elsarticle-num` prints `note` and ignores `annote`, and the flag was appearing in the typeset reference list as "hAND-WRITTEN: internal document, not on INSPIRE" |

`tools/checks/paper_numbers.py` still finds the flag (it reads the whole entry)
and still requires a `note` field beside it. What the four companion-report
entries should eventually carry — a Zenodo DOI or arXiv number at submission
(plans/07 §7.2) — is an author item.
