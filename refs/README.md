# refs/ — reference papers consulted for the reconstruction-chain work

The PDFs are committed in this directory (tracked in git since
2026-08-25 at the user's request, about 240 MB in total; the EIC Yellow
Report is stored as four parts because GitHub rejects single files above
100 MB). This index records what each one was used for and which numbers
were taken from it (2026-08-25 read-through). A missing arXiv copy can be
re-fetched with `python refs/find_ref.py --fetch`; `--check` lists the
local-copy status of every dictionary entry.

**Machine-readable dictionary:** `refs/refs_dict.json` maps every local
file (and the external references the reports lean on) to identifiers,
title, key content with equation/figure/slide pointers, and where it is
used in the repository. Look things up with

```bash
python refs/find_ref.py slot            # keyword search over all fields
python refs/find_ref.py "Eq. (9)" a2    # several terms, all must match
python refs/find_ref.py --key maple     # full entry
python refs/find_ref.py --list
python refs/find_ref.py --check     # local-copy status per entry
python refs/find_ref.py --fetch     # re-download missing arXiv PDFs
```

| file | reference | used for |
|---|---|---|
| `2410.12764v1.pdf` | W. Cosyn, B. Roldan Tomei, A. Sosa, A. Zec, *Polarization options in inclusive DIS off tensor polarized deuteron*, EPJ A 61 (2025) 83 | Exact spin-1 inclusive decomposition for unpolarized electrons (their Eq. 10): the tensor SFs b₁–b₄ feed the cos φ_TL term at O(γ) (Eq. 17d) and the cos 2φ_TT term at O(γ²) (Eq. 17e); ε and γ definitions (Eqs. 11c, 13); axis along q kills T_LT, T_TT (Eq. 20); b₃, b₄ non-negligible at Q² = 2 GeV² for the photon-direction axis (their Fig. 6). Their decomposition carries no Δ (b₁ focus). |
| `2408.13213v1.pdf` | H. Mäntysaari, F. Salazar, B. Schenke, C. Shen, W. Zhao, *Spatial imaging of polarized deuterons at the EIC*, PLB 858 (2024) 139053 | Eq. (9): d²σ/dΦd\|t\| ∝ 1 + 2Σₙ aₙ e^{inΦ} → the cos 2Φ modulation coefficient is **2a₂**; Φ = angle between the vector-meson (recoil) momentum and the polarization axis, polarizations defined in the γ*d c.m. frame ("a Lorentz transformation is required from the lab frame ... which would mix the polarization states" at Q² > 0); Fig. 4 a₂, a₄ per m-state; effective radius R = √(2B_D) from a fit within \|t\| < 0.3 GeV². |
| `0812.2003v3.pdf` | ZEUS Collaboration, *Deep inelastic scattering with leading protons or large rapidity gaps at HERA*, NPB 816 (2009) 1 | Sec. 10.2 azimuthal asymmetries dσ/dΦ ∝ 1 + A_LT cos Φ + A_TT cos 2Φ: A_LT = −0.036 ± 0.036 (+0.016/−0.014), A_TT = −0.030 ± 0.037 (+0.022/−0.006) for 0.0002 < x_P < 0.01; A_LT = +0.051 ± 0.024, A_TT = −0.010 ± 0.024 for 0.01 < x_P < 0.1 → our u₁ = 0.05, u₂ = 0.02 sit at the 1σ edge. LPS: t-resolution σ(t)/t = 0.14 GeV/√\|t\| dominated by the beam angular spread; Φ resolution ≈ 0.2 rad; **beam transverse-momentum spread at the IP ≈ 45 MeV horizontal, 100 MeV vertical** (an anisotropic divergence, ×2); t-slope b = 7.0 ± 0.3 GeV⁻²; LPS acceptance ≈ 2% at x_L > 0.98; 0.09 < \|t\| < 0.55 GeV². |
| `9812212v1.pdf` | N. N. Nikolaev, A. V. Pronyaev, B. G. Zakharov, *Azimuthal asymmetry as a new handle on σ_L/σ_T in diffractive DIS*, hep-ph/9812212 | Eq. (1): the y-dependence of the cos φ and cos 2φ terms of diffractive DIS with a detected proton: u₁ = (2−y)√(1−y) A_LT/[2(1−y)+y²], u₂ = 2(1−y) A_TT/[2(1−y)+y²]; A_LT = F_LT/(F_T+F_L); LT/T ratio model-independent at large β. |
| `9808432v1.pdf` | A. V. Pronyaev, *The forward cone and L/T separation in diffractive DIS*, hep-ph/9808432 | Eq. (7): R_LT/T = (p_⊥/Q)·12β³(2−3β)/[(1−β)(3+4β+8β²)] — at our coherent kinematics (β ≈ 0.5, p_T/Q ≈ 0.2) this gives A_LT ≈ 0.03, consistent with the ZEUS bound; sign change at β = 2/3. |
| `0206031v1.pdf` | M. Ruspa (H1/ZEUS), *Inclusive diffraction at HERA*, hep-ex/0206031 | Diffractive-to-total ratio dσ_diff/dM_X/σ_tot vs W (Fig. 5) — the qualitative basis of the f₀ bracket in plans/06 (no single number quoted in the text). |
| `ePIC_far_forward_talk_DIS_2023_v2.pdf` | A. Jentsch (ePIC), *Far-Forward Detectors and Physics with ePIC @ the EIC*, DIS 2023 | Far-forward acceptance table (RP 0.0* < θ < 5.0 mrad, lower bound optics-dependent; OMD 0–5 mrad; B0 5.5–20 mrad; ZDC < 5.5 mrad); RP technology (500 μm AC-LGAD, "potless" RF-foil design, stations 2 m apart, 25.6 × 12.8 cm sensors); OMD implemented as horizontal RP-style sensors (protons 45% < x_L < 55%); "beam effects the dominant source of momentum smearing" (angular divergence, crossing angle, crab rotation / vertex smearing); RP p_T resolution plots with and without beam effects. |

## Found by web search (2026-08-25): inclusive-DIS kinematic reconstruction at the EIC

The item still missing after the first read-through — Σ / Jacquet–Blondel
δy/y at y ≈ 0.01 — is now bracketed by these documents (downloaded into
refs/):

| file | reference | used for |
|---|---|---|
| `EIC_Seminar_SMaple_2024.pdf` | S. Maple (ePIC), *Tracking and inclusive DIS reconstruction with the ePIC detector at the EIC*, seminar, Birmingham, 11 Dec 2024 (Indico Global 4787) | **The ePIC source for the hadronic-y resolution:** smeared EIC pseudodata (Djangoh 18×275, Q² > 1 GeV²) with σ(δ_h)/δ_h = 25%, σ(p_T,h)/p_T,h = 25%, σ(E_e)/E = 11%/√E ⊕ 2%, σ(θ_e) = 0.1 mrad (slide 44); the same 25% parametrization fitted to full ePIC simulation (Craterlake 23.12.0, Pythia8, Q² > 100 GeV²; slide 47); Δy/y per method in y bins — Σ/JB/DA widths 0.2–0.3 in 0.01 < y < 0.05, electron flat (slides 45, 48); the 18×275 "best method for y" map with 25/10/5/1% markers, DA at the y = 0.01 edge (slide 42); coverage 0.01 < y < 0.95, Q² > 1 GeV², five log bins per decade, "lower y accessible → easier to rely on overlap between data at different √s"; Bayesian kinematic fit incl. ISR energy (slides 43–56). |
| `2210.09048.pdf` | ATHENA Collaboration, JINST 17 (2022) P10019 (43 MB) | Sec. 3.1 / Fig. 22: "at the lowest y values, the electron method resolution degrades like 1/y ... e−Σ or Double Angle methods ... for y ≲ 0.1"; "the resolution with the JB method is at the 20–30% level throughout the kinematic range"; Fig. 22 marker sizes: ≈ 25% at y ≈ 0.01 (DA), ≈ 10% at low Q², y ≈ 0.02–0.1 (e−Σ); y > 0.01 cut motivated by reconstructibility. |
| `2110.05505v2.pdf` | M. Arratia, D. Britzger, O. Long, B. Nachman, NIM A 1025 (2022) 166164 | Table 1 with the formulas of all classical methods; Fig. 5 (ATHENA fast simulation, Q² > 200 GeV², 18×275): RMS(y)/y at y → 0.05 ≈ 0.3 (electron), 0.17 (JB), 0.15 (DA), 0.13 (IΣ); Sec. 5 on why Σ is fragile at low y (noise, acceptance) and the Delphes-vs-full-simulation comparison (Fig. 11). |
| `2206.04897.pdf` | R. Aggarwal, A. Caldwell, JINST 17 (2022) P09035 | The Bayesian kinematic-fit method (HERA kinematics, Q² > 400 GeV²) that the ePIC seminar applies. |
| `2209.14489.pdf` | C. Pecar, A. Vossen, DIS2022 proceedings | SIDIS reconstruction with a particle-flow network vs electron/JB/DA: HFS methods surpass the electron method at very low y. |

| `2108.11638.pdf` | M. Diefenthaler, A. Farhat, A. Verbytskyi, Y. Xu, EPJ C 82 (2022) 1064 | ZEUS-based DNN reconstruction of Q², x (HERA, not EIC) — search record only. |

The ePIC Inclusive WG wiki (open task on hadron treatment in JB/DA/(e)Σ),
the BNL "DIS Kinematics" wiki page and B. Schmookler's
`JeffersonLab/dis-reconstruction` repository (Yellow-Report kinematic
maps) are recorded in `refs_dict.json` as the places to look for newer
ePIC numbers. What is still not published anywhere we found: an ePIC
full-simulation δy/y at Q² ≈ 1–3 GeV², y ≈ 0.01, for e + light ions.

## Downloaded 2026-08-25 (second pass): the external references of the dictionary

Every dictionary entry with a free copy now has one here (`python
refs/find_ref.py --check`). First pages and the cited numbers were checked
against the dictionary; two entries had wrong titles/authors in the repo
text and were corrected (marked below).

| file | reference | used for / verified |
|---|---|---|
| `hep-ph_0611265.pdf` | A. Bacchetta, M. Diehl, K. Goeke, A. Metz, P. J. Mulders, M. Schlegel, *SIDIS at small transverse momentum*, JHEP 02 (2007) 093 | Eqs. (2.4)–(2.5): covariant azimuths cos φ_S, sin φ_S with g_⊥, ε_⊥ — `reco.azimuth_wrt_lepton_plane`. |
| `hep-ph_0410050.pdf` | A. Bacchetta, U. D'Alesio, M. Diehl, C. A. Miller, *Single-spin asymmetries: the Trento conventions*, PRD 70 (2004) 117504 | Sign/orientation conventions of the azimuthal angles about the virtual photon. |
| `hep-ph_0503023.pdf` | M. Diehl, S. Sapeta, *On the analysis of lepton scattering on longitudinally or transversely polarized protons*, EPJ C 41 (2005) 515 | O(γ) mixing of longitudinal/transverse target polarization between the lab and photon frames. |
| `hep-ex_9412004.pdf` | U. Bassler, G. Bernardi, *On the kinematic reconstruction of DIS at HERA: the Σ method*, NIM A 361 (1995) 197 (DESY 94-231) | y_Σ = Σ/(Σ + E′(1 − cos θ_e)), ISR insensitivity; the e−Σ mixed method of `reco.mixed_method`. (Jacquet–Blondel 1979, DESY 79/48, has no free copy.) |
| `2103.05419_part1..4.pdf` | R. Abdul Khalek et al., *EIC Yellow Report*, NPA 1026 (2022) 122447 — arXiv v3, 902 pages, split with PyMuPDF into pages 1–300 / 301–450 / 451–600 / 601–902 (the unsplit 124 MB file is git-ignored) | Sec. 8.1 kinematic-reconstruction comparison; far-forward 10σ Roman-pot cuts; luminosity accounting. |
| `2108.08314.pdf` | A. Jentsch, Z. Tu, C. Weiss, PRC 104 (2021) 065205 | Table I verified: B0 5.5–20 mrad; OMD 0–5 mrad, ξ 0.45–0.65; RP 0–5 mrad, ξ 0.6–0.95; ZDC 0–4 mrad — the Phase-1 far-forward model of `fastsim/polli_fastsim/farforward.py`. |
| `2511.05638.pdf` | W. Chang, E.-C. Aschenauer, A. Jentsch, A. Kumar, Z. Tu, Z. Yin, PRD 113 (2026) 032018 | Verified: IR-8 intact-recoil efficiencies 47.12 / 32.23 / 29.42 / **17.75** / 12.37 / 6.36 / 1.59% for d / ³He / ⁴He / ⁷Li / ⁹Be / ¹²C / ¹⁶O; no ⁶Li sample. |
| `2603.23699.pdf`, `2603.23700.pdf` | W. Cosyn, C. Weiss, *SIDIS on a polarized spin-1 target. I. Cross section and spin observables; II. Deuteron and spectator nucleon tagging* (JLAB-THY-26-4663/-4661) — **titles corrected** (the repo had called them "tagged DIS on spin-1 targets") | Part II abstract verbatim: "Tensor-polarized spin asymmetries of order unity are achieved for spectator momenta ≳ 300 MeV, which select configurations with large D-wave" (plans/01 §2.4, primer §4.4). |
| `2510.10794.pdf` | G. Atoian et al. (EPIOS), *Realizing the scientific program with polarized ion beams at EIC*, PRC 113 (2026) 060501 | Verified: G = 1.793 / −0.143 / −4.184 / −0.178 / 1.532 for p / d / ³He / ⁶Li / ⁷Li; "⁶Li, ¹⁴N which, together with the deuteron, are the only stable spin-one nuclei"; Δ(x,Q²) needs transversely polarized spin ≥ 1 nuclei. |
| `2509.18558.pdf` | E. Hamwi, G. H. Hoffstaetter, *Polarization transmission in the EIC's Hadron Storage Ring* (Cornell; PRAB 29 (2026) 073501 as cited in plans/01) — **title and author list corrected** (the repo had "Hamwi–Devlin–Hoffstaetter, *Spin dynamics of light polarized ions in the EIC hadron ring*"; J. Devlin is acknowledged, not an author) | Table I verified: G / max\|Gγ\| / resonances = p 1.7928 / 525 / 1575; d −0.1430 / 21 / 63; ³He −4.1842 / 819 / 2457; ⁶Li −0.1818 / 27 / 81; ⁷Li 1.5196 / 191 / 573; "particles with small anomalous magnetic moment (²H, ⁶Li) are not amenable to Siberian snake devices". |
| `TUNL_A6_2002.pdf` | D. R. Tilley et al. (TUNL), *Energy levels of light nuclei A = 6*, NPA 708 (2002) 3 — the 2017 revised manuscript from nucldata.tunl.duke.edu (`ourpubs/06_2002.pdf`) | Verified: 2.186 MeV 3⁺;0 Γ = 24 ± 2 keV (α, d); 5.366 MeV 2⁺;1 Γ = 541 ± 20 keV (γ, n, p, α); Q_m(⁴He(d,γ)⁶Li) = 1.4743 MeV; 3.5629 MeV 0⁺ α + d decay via (1996CS03) — `coherent.LI6_BREAKUP`, plans/06 §6.2. |

No free copy exists for: Hoodbhoy–Jaffe–Manohar NPB 312 (1989) 571 and
Jaffe–Manohar PLB 223 (1989) 218 (the Δ definition), Sather–Schmidt PRD 42
(1990) 1424, Jacquet–Blondel DESY 79/48, and Li–Sick–Whitney–Yearian NPA
162 (1971) 583 (`find_ref.py --check` lists them as "no free copy"); the
dictionary entries carry what the repo takes from each.

## Downloaded 2026-08-26: superconducting nanowire and microwire detectors

The far-forward near-beam study (`plans/09`, `reports/nanowire_far_forward.html`)
asks whether a superconducting nanowire layer can reach closer to the ⁶Li beam
than the ePIC Roman-Pot silicon does. These are its primary sources — three from
the Argonne MEP/Physics-Division group that develops the technology, two from
the Caltech/JPL/Fermilab group that develops the large-area variant.

| file | reference | used for / verified |
|---|---|---|
| `2312.13405.pdf` | S. Lee, T. Polakovic, W. Armstrong, A. Dibos, T. Draher, N. Pastika, Z.-E. Meziani, V. Novosad, *Beam Tests of SNSPDs with 120 GeV Protons*, NIM A 1069 (2024) 169956 | The Argonne device as actually built and beam-tested: **12 nm NbN** on 300 μm Si, **30 × 30 μm² active area**, fill factor 1/2, 8 devices per 8 × 8 mm² chip; wire widths 300–800 nm with ~250 nm named ideal and >400 nm inefficient at low bias; **T = 2.82 K** (GM cryocooler), Tc ≈ 7 K, I_c = 25.2 μA; background/signal (2.08 ± 0.89)% at FTBF MT6.2, with background rising exponentially above I_b/I_c > 0.8. **The hot-spot anchor**: their Eqs. 1–2 give the firing threshold I_th/I_c = 1 − 2r_s/w, and their thresholds versus wire width fit r_s = **134 nm** for a 120 GeV proton, crossing zero at w ≈ 268 nm (hence the ~250 nm optimum). **Read that 134 nm as a fit parameter, not a datum**: it is an *extrapolated* zero-crossing of four points, the authors write "While the physical validity of this simple model is a question of future work", and inverting the four points individually gives 102–120 nm. Note also that Q here is *"the energy that the proton has deposited into the thin film"* (≈20 eV), while `2601.03158` applies the same √Q scaling to the **substrate** deposit (0.1 MeV) — a factor 5000 in what Q means. Radiation hardness explicitly **not yet measured** ("planned in the near future"). Their own EIC motivation is the aperture argument: acceptance "is limited by the beamline magnets, which can be mitigated by operating SNSPDs within the frigid bore of superconducting magnets". |
| `1907.13059.pdf` | T. Polakovic, W. R. Armstrong, V. Yefremenko, J. E. Pearson, K. Hafidi, G. Karapetrov, Z.-E. Meziani, V. Novosad, *Superconducting nanowires as high-rate photon detectors in strong magnetic fields*, NIM A 959 (2020) 163543 | The field result the siting argument needs: saturated internal efficiency to **5 T parallel** to the device plane (setup limit; ~8 T extrapolated) but only **~0.5 T perpendicular** — orientation, not magnitude, is the constraint. NbN meander 13.5 nm × 80 nm on 110 nm pitch, 10 × 10 μm² pixel, Tc = 8 K, H_c2(0) = 32 T; saturation at 9 μA against I_c = 23 μA; **10⁷ counts/s measured, 10⁸ /s expected** from τ_F = 11.78 ns. |
| `2510.11725.pdf` | C. Wang et al. (Caltech / JPL / Fermilab), *Towards High-Efficiency Particle Detection Using Superconducting Microwire Arrays*, submitted to JINST | The large-area branch: 8-channel **1 × 1 mm²** WSi SMSPD on a 4.7 nm film, 1 μm wires on a 3 μm gap → **25% fill factor**, **0.8 K** operation (Tc = 1.85 K). At CERN SPS H6 with 120 GeV hadrons and muons: fill-factor-normalized efficiency **75%**, time resolution **130 ± 17 ps** (up from 60% and 1.15 ns on the earlier 3 nm / 1.5 μm device). mm² of area bought with a colder fridge and a quarter fill factor. |
| `2601.03158.pdf` | S. Lee, W. Armstrong, J. DiPreta, C. Dulya, V. Novosad, T. Polakovic, *Optimization of Cryogenic Detector Test Station by Rejecting Electromagnetic Interference*, NIM A 1093 (2027) 171953 | The **second anchor of the hot-spot scaling**: "the previous work with a 120 GeV proton determined the hot spot size to be 134 nm, and a 5.5 MeV α particle creates an approximately 1 μm hot spot by the same scaling". A relativistic ⁶Li lands at r_s ≈ 400 nm, between their two points — but **in Q, not in z**: this α differs from the 120 GeV proton almost entirely through 1/β² (β = 0.054 against ≈1), so it calibrates √Q across energy and says nothing about z² at fixed velocity. Nobody has varied Z at fixed β on one of these devices. A ~1 μm hot spot does **not latch** a 100–200 nm wire (clean count-rate plateaus). Low-bias running needed an EMI-rejection scheme to see counts at all — a caution for any "bias down until only Z = 3 fires" scheme. The α analysis itself is not published ("a detailed study of α detection is underway"). |
| `2410.00251.pdf` | C. Peña, C. Wang, S. Xie et al. (Caltech / JPL / Fermilab), *Characterization of a Superconducting Microwire Single-Photon Detector Array for Charged-Particle Detection*, JINST 20 (2025) P03001 | The largest SMSPD run in a GeV hadron beam: 8 channels over **2 × 2 mm²**, 3 nm WSi, 1.5 μm wires, 40% fill, 30 ± 1 μm spatial and 1.15 ns time resolution. The measurement that **kills pulse-height Z-ID**: waveforms and amplitude distributions from 120 GeV protons, 8 GeV pions and showering electrons are indistinguishable, and one amplitude threshold serves all — the amplitude is the diverted bias current, not the deposit. |

## Downloaded 2026-08-28: the theory curves the fast simulation draws

Phase-1 step 1.2 asked for the polarized-EMC and b₁ predictions "digitized
as `medium_ratio(x)`" rather than approximated by constants. These are the
five papers those curves come from. Four of the five figures are now
committed as CSV in `fastsim/polli_fastsim/data/`, read back from the PDFs'
own path operators by `tools/digitize_figure.py`; that directory's
`SOURCES.md` records page, figure, frame box, axis ranges, the legend
handle each curve was identified by, and the exact command.

| file | reference | used for / verified |
|---|---|---|
| `nucl-th_0605061.pdf` | I. C. Cloët, W. Bentz, A. W. Thomas, *Spin-dependent structure functions in nuclear matter and the polarized EMC effect*, PLB 642 (2006) 210 | **The CBT camp.** Page 7 FIG. 6 upper-left panel is ⁷Li at Q² = 5 GeV² — the only ⁷Li-specific polarized-EMC calculation in the literature, and therefore the common baseline of the two-camp comparison. Three curves digitized into `cbt_polemc_7Li_Q5.csv`: the unpolarized EMC ratio (blue dashed), R^{(3/2 1)}_As of their Eq. (26) (red solid) and R^{3/2 3/2}_As of their Eq. (23) (red dotted), which is what `plans/01` defines ΔR_A to be and what the money plot draws. **What the digitization corrects:** the "2× the unpolarized effect" is a valence-region statement, not a constant — (1 − R_pol)/(1 − R_unpol) is 2.25 / 1.69 / 1.41 / 1.14 at x = 0.40 / 0.45 / 0.50 / 0.60, the ratio bottoms out at 1.06 near x = 0.70 without reaching 1 (their Eq.-26 curve does cross the unpolarized one, at x = 0.651), and below x ≈ 0.28 the ratio has no meaning because ⁷Li's unpolarized ratio is above 1 there while the polarized one keeps a 7% depletion. Curves span x = 0.028–0.871. |
| `1806.00481.pdf` | S. Tronchin, H. H. Matevosyan, A. W. Thomas, *Polarized EMC effect in the QMC model*, PLB 783 (2018) 247 | **The TMT camp.** Page 9 Figure 4, isospin-symmetric nuclear matter at Q² = 10 GeV², unpolarized (blue solid) and polarized (purple dashed), into `tmt_polemc_nm_Q10.csv`; x = 0.0015–0.739. Their ratio of effects is 1.01 / 0.98 / 1.00 / 1.08 at x = 0.40 / 0.45 / 0.50 / 0.60 — "polarized ≈ unpolarized" holds pointwise across the valence region. Different target and scale from CBT, which is why the repository transfers each model's effect onto one common unpolarized baseline — EPPS21's ⁶Li per-nucleon F₂ over CT18ANLO since 2026-08-29 — with a single valence strength factor (0.5322 for CBT, 0.2113 for TMT; 1 and 0.397 on the legacy `cbt` baseline, which is what every figure published before that date used) instead of subtracting published R values.  That factor is fitted over 0.35 < x < 0.65 and applied everywhere, so it is worth knowing what it does outside its window: the two PUBLISHED polarized curves agree to better than 0.008 over 0.028 < x < 0.30 (0.002 at x = 0.09), and the ≈ 0.02 separation (0.0212–0.0243) the transferred pair shows there is the rescaling, not a disagreement between the papers.  Inside the window it is the comparison the two camps make: the transferred nuclear-matter depletion tracks the baseline's own unpolarized 0.021 / 0.027 / 0.041 at x = 0.40 / 0.45 / 0.65 to within 0.003 against the transferred CBT's 0.041 / 0.043 / 0.050, so ΔR separates by 0.021 at x = 0.36 and by 0.006 at 0.65 (⁷Li's own unpolarized curve, 0.034 / 0.048 / 0.087, and the untransferred pair 0.077 / 0.082 / 0.094 against 0.039 / 0.048 / 0.083, separating by 0.040 and 0.011, are the pre-2026-08-29 reading `--emc-baseline cbt` returns). |
| `1702.05337.pdf` | W. Cosyn, Yu-Bing Dong, S. Kumano, M. Sargsian, *Deuteron tensor structure function b1*, PRD 95 (2017) 074036 | **The convolution camp for b₁.** Page 9 FIG. 4 (x·b₁ at Q² = 2.5 GeV², SD / DD / sum for two convolution formalisms) into `b1_cdks_q2p5.csv` and page 10 FIG. 5 (the same sums at Q² = 1.0 / 2.5 / 5.0) into `b1_cdks_q2set.csv`. The extraction validates itself twice: the solid curves equal the sum of their own dashed and dotted ones to 2×10⁻⁷, and Fig. 5's Q² = 2.5 curve reproduces Fig. 4's to 2×10⁻⁶ from a different page with a different calibration. \|b₁\| < 10⁻³ at x ≳ 0.2 confirmed; the sum changes sign at x ≈ 0.06 and again at 0.42, which the `0.1 × toy_b1` stand-in it replaces never did. ∫b₁ dx over the digitized range is 4.6×10⁻⁴ — consistent with the Close–Kumano sum rule over the range the figure covers, which is the sharpest statement the table supports: it stops at x = 1.59 while the deuteron's x runs to 2, so `close_kumano_integral` reports the number and does not enforce it. |
| `1311.4561.pdf` | G. A. Miller, *Pionic and hidden-color, six-quark contributions to the deuteron b1 structure function*, PRC 89 (2014) 045203 | **The HERMES-like camp for b₁.** Page 10 FIG. 5, the total b₁ = b₁^π + b₁^{6q} that a hidden-colour probability of 0.15% suffices to reproduce the data with, into `b1_miller.csv` (x = 0.010–0.900); page 11 FIG. 6, 100 b₁ at Q² = 1.17 / 1.76 / 2.12 / 3.25 GeV², into `b1_miller_q2set.csv`. At x = 0.012 the digitized total is 0.114 per deuteron against HERMES's measured 0.112 ± 0.055 ± 0.028 and his own TABLE I. Both are Mathematica plots with no closed frame, so the axis box came from the axis lines and their major ticks; FIG. 6 has no legend either and the curve → Q² map is read from the caption's ordering rule, the only assignment in `data/` not taken off a legend handle. ∫b₁ dx = 5.9×10⁻³ — his Sec. V shows the pionic contribution violates Close–Kumano, and it does. |
| `2109.03591.pdf` | X. G. Wang, W. Bentz, I. C. Cloët, A. W. Thomas, *Polarized gluon EMC effect*, J. Phys. G 49 (2022) 03LT01 | The gluon-spin arm of `plans/02` step 1.2.2: page 8 Figure 3 carries g₁A/g₁p and Δg_A/Δg_p, with ⁷Li named the most promising case and the polarized gluon EMC effect larger than the unpolarized one. **Not digitized** — the dg₁/dlnQ² observable has no money plot yet; the entry records where the curves are. |
