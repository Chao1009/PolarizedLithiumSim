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
| `2511.05638.pdf` | W. Chang, E.-C. Aschenauer, A. Jentsch, A. Kumar, Z. Tu, Z. Yin, PRD 113 (2026) 032018 | Verified: IR-8 intact-recoil **acceptances** 47.12 / 32.23 / 29.42 / **17.75** / 12.37 / 6.36 / 1.59% for d / ³He / ⁴He / ⁷Li / ⁹Be / ¹²C / ¹⁶O; no ⁶Li sample.  The paper calls these “detection efficiency”, but its Sec. IV (Event Generator) states that “the current simulation only accounts for the acceptance effect and does not incorporate the efficiencies of the detector” and that it does not account for the efficiency and acceptance of the reconstructed distribution (verified from the PDF, 2026-09-15).  Read every number in this row as acceptance only — an UPPER BOUND on efficiency × acceptance. |
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

## Downloaded 2026-09-06 (run 18): the sources of the sourcing pass

Run 18 re-read, against the primary sources, the numbers Reports 0–4 take
from the ⁶Li/⁷Li structure calculation, the polarized-target literature,
the JLab proposals, the EIC machine limits and the ePIC tracking
performance. Four of those sources are small enough to commit:

| file | reference | used for / verified |
|---|---|---|
| `1309.3794.pdf` | R. B. Wiringa, R. Schiavilla, S. C. Pieper, J. Carlson, *Nucleon and nucleon-pair momentum distributions in A ≤ 12 nuclei*, PRC 89 (2014) 024305 | **The source of every spin-projected nucleon number in the programme.** Table I (p. 3, caption "Total number of spin-up/down and proton/neutron nucleons in J > 0 nuclei with M_J = J") verified: ⁶Li(1⁺) 1.924 / 1.076 for both protons and neutrons, whence the 0.848 the reports round to 0.85; ⁷Li(3/2⁻) 1.934 / 1.066 and 1.981 / 2.019, whence +0.868 and −0.038. The Hamiltonian is AV18 + **Urbana X**, not IL7. §II carries the cluster normalisations verbatim: "The integrated N_αd = 0.86 is a sum of S- and D-wave parts of 0.846 and 0.017, respectively" (relative 1s, node at k ≈ 0.7 fm⁻¹) and "with N_αt = 1.00" (relative 1p). The paper's own online tables at `phy.anl.gov/theory/research/momenta/` are later runs and give 0.820 and 1.0084 — the reports quote the paper, and `.../overlaps/`, cited until this run, is single-nucleon pickup/stripping and holds no α–d or α–t entry at all. Note that the +0.866 / −0.037 in `beams.py` are JLab PR12-14-001's numbers, not Table I's. |
| `hep-ex_9904002.pdf` | P. L. Anthony et al. (E155), *Measurement of the deuteron spin structure function g₁ᵈ(x)*, PLB 463 (1999) 339 | The ⁶LiD target paragraph behind Report 0 §§2.2 and 3.2: "The ⁶Li polarization was measured to be 97% of the free deuteron polarization"; "we conclude that the effective deuteron in ⁶Li has a net polarization of 87% of the Li polarization" — *effective deuteron* is E155's own phrase, while the companion target paper calls the same number a nucleon polarization P_n = 0.866 ± 0.012. Dilution: f = 0.18–0.20 for the free deuterons, and "C₁ f then gives an effective dilution factor of ∼ 0.36, as compared with ∼ 0.22 for ¹⁵ND₃" — both are *effective* factors, and the 0.22 the reports quote is exactly this one. Target composition 18% free D, 53% ⁶Li, 14% ⁴He, 11% Al, 3% O, 1% N by weight. |
| `hep-ex_0703049.pdf` | P. Abbon et al. (COMPASS), *The COMPASS experiment at CERN*, NIM A 577 (2007) 455 | §4.1: "the fraction of polarisable material f is of the order of 0.35, taking into account also the He content in the target region", against "(f ≈ 0.15)" for irradiated NH₃ — the two dilution numbers Report 0 §3.2 cites to [31], verified verbatim. Also the 4.2% ⁷Li isotopic admixture of the COMPASS ⁶LiD. What is **not** in this paper: the +54%/−47% deuteron polarizations (Ball, NIM A 498), the 51 ± 3% / 92 ± 4% Li polarizations (Koivuniemi, SPIN2004), and any ⁶LiD f above 0.35. |
| `2105.13564.pdf` | B. R. Gamage et al., *Design Concept for the Second Interaction Region for Electron-Ion Collider*, Proc. IPAC'21 TUPAB040, pp. 1435–1438 | Table 1, row "Minimum Δ(Bρ)/(Bρ) allowing for detection of p_T = 0 fragments": **0.1 for EIC IR 1, 0.003–0.01 for EIC IR 2** — the 0.1 that Report 4 §7 uses as the IP6 dispersive-tagging threshold, verified here and against the version of record. It is a detection threshold, not a resolution. The paper publishes **no** per-species efficiencies, and its two entries support a 10–33× ratio, not the 30–60× the report states. The title Report 4 [12] gave this arXiv id belongs to no paper on INSPIRE; the title above is the one on the arXiv, JACoW and INSPIRE records. |

Added to `refs_dict.json` in the same pass, with `"file": null` and a URL
or DOI, because they are paywalled, too large to commit, or not documents
at all: the **EIC Conceptual Design Report** (BNL-221006-2021-FORE,
doi:10.2172/1765663, 972 pp. — the ξ_p ≤ 0.015 hadron beam-beam limit at
§1.2 p. 5, §1.4 p. 11, §3.1.1 p. 95 and §4.6.1 pp. 392/395, Tables
3.3/3.4/3.5, the Z/A scaling at §1.4 p. 14 and the two-IR ion cap at
Appendix A §A.2.2 p. A-12); the **ePIC Preliminary Design Report** of
September 2024 (Zenodo 13866213 — Fig. 8.9 p. 36, Table 8.1 p. 30,
Roman-Pot rates p. 203) and the **Preliminary TDR v3.1** (concept
doi:10.5281/zenodo.18271601, record 18271602 — Fig. 3.56 p. 64, Table 3.22
p. 63, p. 270), whose full-simulation δp/p is 0.38% in the barrel at
1 GeV/c and has no η slice at 0.6%; the COMPASS and Saclay target papers
**J. Ball et al.**, NIM A 498 (2003) 101 (+54.2%/−47.1%), **J. Koivuniemi
et al.**, SPIN2004 p. 796 (51 ± 3% for ⁶Li, 92 ± 4% for ⁷Li) and **P.
Chaumette et al.**, AIP Conf. Proc. 187 (1989) 1275 (47%/56% in ⁷LiH),
together with **S. Bültmann et al.**, NIM A 425 (1999) 23 (SLAC-PUB-7904),
the E155 target paper; the JLab proposals **PR12-14-001** (65–80% for ⁷Li
at 5 T, P_z = 0.8, 15 nA) and **PR12-13-011** with its 2023 jeopardy
update (0.16 < x < 0.49, P_zz = 20% in the proposal and 26% in the update,
the 30% condition lifted in August 2022); **N. J. Stone**, INDC(NDS)-0833
(2021) and INDC(NDS)-0794 (2019), the quadrupole and dipole moment tables;
and **X. Li**, arXiv:2305.15593 — recorded, and flagged, as the 1.4 T
BaBar reference design measured with the ECCE MC, **not** ePIC, so that
its ≈ 0.59% at 1 GeV/c is never re-imported as an ePIC number.

## Downloaded 2026-09-16: the last two open-access gaps from the LiPolGen reference-research pass

The 2026-09-16 `docs/references/` pass across LiPolGen's domains (tensor
b₁, ⁶Li/⁷Li nuclear data, radiative corrections, coherent diffraction, the
generator chain, nucleon/deuteron data) named 148 candidate arXiv ids;
136 of them were already on disk (most fetched earlier the same day).
These are the only two gaps it left, both marked "yes (not fetched)" by
the domain write-ups:

| file | reference | used for / verified |
|---|---|---|
| `hep-ph_0309123.pdf` | G. I. Gakh, O. Shekhovtsova, *Radiative events in DIS of unpolarized electron by tensor polarized deuteron: Radiative corrections*, eConf C030626 (2003) 351 (PIC03 proceedings) | The 3-page conference precursor of the JETP paper already on disk (`hep-ph_0403262.pdf`, `gakh-shekhovtsova2004-tensor-rc`); kept only to date the calculation (`docs/references/01_radiative.md` R15). Verified: %PDF-1.4, 3 pages, 117975 bytes, pdftotext reproduces the title verbatim. |
| `2605.00502.pdf` | I. Helenius, J. O. Laulainen, C. T. Preuss, *Multiplicative matching of neutral current deep-inelastic scattering processes at next-to-leading order in PYTHIA 8*, submitted (no journal ref / DOI on the arXiv record as of this pass) | The NLO successor of `2410.20950` (Helenius–Laulainen–Preuss, JHEP 05 (2025) 153, already on disk): multiplicative matching of NC DIS with the PYTHIA 8 shower (default and Vincia) at O(α_s), validated against HERA reduced cross sections (`docs/references/01_generators-chain.md` GC-6). Verified: %PDF-1.7, 1985710 bytes, pdftotext reproduces the title verbatim. |

Both are open access on arXiv; both entries were appended to `refs_dict.json`.

## Downloaded 2026-09-16 (reconciliation pass): 101 unindexed files, plus six free-not-fetched items

The 2026-09-16 reconciliation counted 194 PDFs on disk against 89 files
named in `refs_dict.json`'s `file` fields and found **101** committed PDFs
that no dictionary entry named. Each was identified from the arXiv API
(`export.arxiv.org/api/query?id_list=<id>`) or, for the ten non-arXiv
tags, from INSPIRE and the PDF's own first page; the rationale for each
(what it is cited for, and where) was taken from the matching bullet in
`docs/references/REFERENCES.md` in the LiPolGen repository, which had
independently re-verified every one of these citations against arXiv,
INSPIRE, Crossref or the PDF text in the same pass. Every title below was
checked against `pdftotext -f 1 -l 1` of the committed PDF before being
written; two disagreed with the arXiv API's own `<title>` field (arXiv
returns a stale/mis-cased string) and two disagreed with
`REFERENCES.md`'s prose citation (a title mismatch in that document,
not in the PDF) — in all four cases the PDF's own first page is what was
kept. `hep-ph_9302320.pdf`'s and `JINR-P2-10113.pdf`'s OCR/font extraction
is degraded (old scans), but report number, author list and topic on the
first page(s) still confirm the citation.

| file | reference | used for / verified |
|---|---|---|
| `0712.2670.pdf` | G. Watt, H. Kowalski, *Impact parameter dependent colour glass condensate dipole model*, Phys. Rev. D 78 (2008) 014016 | Sartre's `DipoleModel_bCGC` -- the model-systematic partner of the two above. |
| `0901.0002.pdf` | A. D. Martin, W. J. Stirling, R. S. Thorne, G. Watt, *Parton distributions for the LHC (MSTW 2008)*, Eur. Phys. J. C 63 (2009) 189 | A gate-defining reference, not a convenience: `b1_nuclear.hpp` states the CDKS-comparison gate passes for the MSTW2008 LO configuration specifically, because MSTW2008 LO is the unpolarized PDF CDKS computed their b1 with. |
| `0909.1254.pdf` | A. Caldwell, H. Kowalski, *The J/psi Way to Nuclear Structure*, arXiv preprint (2009); the refereed version (Phys. Rev. C 81 (2010) 025203) carries no free e-print | The \|t\| -> b Fourier-Bessel inversion that gives `slope_b` its meaning as an imaging observable rather than a fit parameter -- the argument `cluster_config.hpp`'s `eps_b0_equivalent()` makes when it turns a configuration ensemble into a slope. |
| `1005.4524.pdf` | S. Kumano, *Tensor-polarized quark and antiquark distribution functions in a spin-one hadron*, Phys. Rev. D 82 (2010) 017501 | `01_tensor-b1.md` section 4.1 calls it "the only b1 source in this domain that publishes numbers": a closed-form delta_T w(x) with fitted parameters at Q2 = 2.5 GeV2, fitted to HERMES. |
| `1106.2091.pdf` | L. Frankfurt, V. Guzey, M. Strikman, *Leading twist nuclear shadowing phenomena in hard processes with nuclei*, Phys. Rept. 512 (2012) 255 | The missing `f0` calculation class. |
| `1202.2225.pdf` | G. I. Gakh, M. I. Konchatnij, N. P. Merenkov, *Radiative corrections to polarization observables in elastic electron-deuteron scattering in leptonic variables*, JETP 115 (2012) 212 | The only independent tensor-target radiative calculation besides Gakh-Shekhovtsova -- the one thing available to check POLRAD Eq. (A.4) against. |
| `1211.3048.pdf` | T. Toll, T. Ullrich, *Exclusive diffractive processes in electron-ion collisions*, Phys. Rev. C 87 (2013) 024913 | The method behind Sartre, and the statement that F(b) is recovered from the coherent dsigma/dt -- the imaging claim `slope_b` stands in for. |
| `1212.2974.pdf` | A. H. Rezaeian, M. Siddikov, M. Van de Klundert, R. Venugopalan, *Analysis of combined HERA data in the Impact-Parameter dependent Saturation model*, Phys. Rev. D 87 (2013) 034002 | The current IPsat parameters from combined H1+ZEUS data -- what a new 6Li dipole table would be generated with. |
| `1307.8059.pdf` | T. Toll, T. Ullrich, *The dipole model Monte Carlo generator Sartre 1*, Comput. Phys. Commun. 185 (2014) 1835 | Kept for its negative: section 3.3.4 enumerates the supported nuclei and A = 6 is not among them, and its nucleus model is Woods-Saxon only -- which is why `open_items/physics_literature.md` Item 7 rules Sartre out as the route to a tensor axis. |
| `1311.4835.pdf` | K. Slifer, E. Long, *Novel Physics with Tensor Polarized Deuteron Targets*, PoS PSTP2013 (2013) 008 | Programme-level motivation plus the target-technology context for 6LiD / ND3. |
| `1406.5539.pdf` | E. R. Nocera, R. D. Ball, S. Forte, G. Ridolfi, J. Rojo, *A first unbiased global determination of polarized PDFs and their uncertainties (NNPDF)*, Nucl. Phys. B 887 (2014) 276 | `LhapdfG1`'s default `NNPDFpol11_100`, including the replica-based uncertainties and the Deltab = 0 convention the tree's header note attributes to this fit. |
| `1407.1653.pdf` | W. Cosyn, M. Sargsian, *Final-state interactions in deep-inelastic scattering from a tensor polarized deuteron target*, J. Phys. Conf. Ser. 543 (2014) 012006 | The only calculation anywhere of A_zz with any non-Born effect in it. |
| `1407.3852.pdf` | S. Kumano, *Tensor-polarized structure functions: Tensor structure of deuteron in 2020's*, J. Phys. Conf. Ser. 543 (2014) 012001 | Carries closed-form projection operators extracting b1-b4 from the hadron tensor W_munu -- a free, independent check on `InclusiveKernel`'s transcription of the Cosyn et al. |
| `1412.7420.pdf` | A. Buckley, J. Ferrando, S. Lloyd, K. Nordstrom, B. Page, M. Ruefenacht, M. Schoenherr, G. Watt, *LHAPDF6: parton density access in the LHC precision era*, Eur. Phys. J. C 75 (2015) 132 | `[LHAPDF6]` -- the grid library behind `src/lhapdf/lhapdf_sf.cpp`, i.e. the whole optional LHAPDF tier and every PDF-set name quoted in section 2 and section 4c. |
| `1509.00792.pdf` | K. Kovarik et al., *nCTEQ15 - Global analysis of nuclear parton distributions with uncertainties in the CTEQ framework*, Phys. Rev. D 93 (2016) 085037 | Second fit of the `T-39` four-fit EMC spread, and the source of a real 7Li grid (`nCTEQ15_7_3`) -- `01_generators-chain.md` `T-40` records that 7Li is therefore not PDF-less, contrary to the tree's current assumption. |
| `1606.03149.pdf` | S. Kumano, Q. Song, *Theoretical estimate on tensor-polarization asymmetry in proton-deuteron Drell-Yan process*, Phys. Rev. D 94 (2016) 054022 | The Fermilab route to delta_T qbar using the 2010 PDFs above -- the non-DIS way the same distributions could be measured. |
| `1606.04505.pdf` | W. Detmold, P. E. Shanahan, *Gluonic Transversity from Lattice QCD*, Phys. Rev. D 94 (2016) 014507; Erratum PRD 95 (2017) 079902 | The lattice half of `CoherentScenario::amp`'s 3 x 10-3 ... |
| `1607.01711.pdf` | H. Mantysaari, B. Schenke, *Revealing proton shape fluctuations with incoherent diffraction at high energy*, Phys. Rev. D 94 (2016) 034042 | The paper behind `subnucleondiffraction`, the code `write_snd_configs()` writes configuration files for. |
| `1607.03838.pdf` | S. R. Klein, J. Nystrand, J. Seger, Y. Gorbunov, J. Butterworth, *STARlight: A Monte Carlo simulation program for ultra-peripheral collisions of relativistic ions*, Comput. Phys. Commun. 212 (2017) 258 | The Q2 -> 0 parent of eSTARlight, with a decade of RHIC/LHC UPC data behind it, and the documented Z <= 6 Gaussian form-factor branch -- the published basis both for "6Li runs today" and for `slope_b`'s Gaussian model. |
| `1709.00395.pdf` | F. Winter, W. Detmold, A. S. Gambhir, K. Orginos, M. J. Savage, P. E. Shanahan, M. L. Wagman, *First lattice QCD study of the gluonic structure of light nuclei (NPLQCD)*, Phys. Rev. D 96 (2017) 094512 | The first exotic-glue number in a nucleus. |
| `1709.07827.pdf` | S. Cotogno, T. van Daal, P. J. Mulders, *Positivity bounds on gluon TMDs for hadrons of spin <= 1*, JHEP 11 (2017) 185 | A free positivity gate on `toy_delta_gluon`: leading-twist gluon-distribution bounds including tensor polarization, with a small-x limit. |
| `1710.00391.pdf` | B. Cabouat, T. Sjostrand, *Some Dipole Shower Studies*, Eur. Phys. J. C 78 (2018) 226 | The DIS-specific pedigree of the shower the bridge inherits (dipole recoil) -- what `D-1` is actually testing when it tests "PYTHIA's DIS shower". |
| `1803.06420.pdf` | M. Lomnitz, S. Klein, *Exclusive vector meson production at an electron-ion collider (eSTARlight)*, Phys. Rev. C 99 (2019) 015203 | `estarlight_li6_coherent()` and `estarlight_li6_q2_floors()` -- the code every O5 coherent rate in the tree comes from (`D-4`). |
| `1803.11206.pdf` | J. Maxwell, D. Crabb, D. Day, W. Detmold, R. Jaffe, M. Jones, C. Keith, D. Keller, D. Meekins, R. Milner et al., *Search for Exotic Gluonic States in the Nucleus, A Letter of Intent to Jefferson Lab PAC 44*, Letter of Intent to Jefferson Lab PAC 44 (dated 6 June 2016, posted 2018) | The only experimental programme anywhere aiming at Delta: unpolarized electron beam on a transversely polarized spin-1 nuclear target, inclusive DIS below x = 0.3, via single-spin tensor asymmetries. |
| `1805.05877.pdf` | S. Fucini, S. Scopetta, M. Viviani, *Coherent deeply virtual Compton scattering off 4He*, Phys. Rev. C 98 (2018) 015203 | The impulse-approximation convolution that TOPEG implements -- the citable primary source for a generator that has no code paper. |
| `1806.10820.pdf` | C. Bierlich, G. Gustafson, L. Lonnblad, H. Shah, *The Angantyr model for Heavy-Ion Collisions in PYTHIA8*, JHEP 10 (2018) 134 | Documents why a nucleus id must never reach `Beams:idA/idB` -- the published basis for the negative result recorded in `benchmarking/01` section 2.2. |
| `1908.03355.pdf` | M. Walt, I. Helenius, W. Vogelsang, *Open-source QCD analysis of nuclear parton distribution functions at NLO and NNLO (TUJU19)*, Phys. Rev. D 100 (2019) 096015 | The source of the `TUJU19_*_7_3` 7Li grids named in `T-40`. |
| `1910.12523.pdf` | S. Kumano, Q. Song, *Gluon transversity in polarized proton-deuteron Drell-Yan process*, Phys. Rev. D 101 (2020) 054011 | The Drell-Yan half of the same band -- the "Drell-Yan estimates" `coherent.hpp`'s bound comment refers to. |
| `1912.08005.pdf` | A. Buckley, P. Ilten, D. Konstantinov, L. Lonnblad, J. Monk, W. Pokorski, T. Przedzinski, A. Verbytskyi, *The HepMC3 Event Record Library for Monte Carlo Event Generators*, Comput. Phys. Commun. 260 (2021) 107310 | The library LiPolGen writes with (3.3.0 in `deps/install`) and the format the whole ePIC chain consumes. |
| `1912.10053.pdf` | T. Hou et al. (CTEQ-TEA), *New CTEQ global analysis of quantum chromodynamics with high-precision data from the LHC*, Phys. Rev. D 103 (2021) 014013 | `LhapdfSF`'s default `CT18NLO`, and `CT18Anlo` is EPPS21's proton baseline -- so one paper stands behind both the free-nucleon F2 and the nuclear-ratio denominator. |
| `2003.06623.pdf` | S. Kumano, Q. Song, *Deuteron polarizations in the proton-deuteron Drell-Yan process for finding the gluon transversity*, Phys. Rev. D 101 (2020) 094013 | The same, re-expressed in conventional (measurable) deuteron polarizations. |
| `2005.14706.pdf` | Z. Tu, A. Jentsch, M. Baker, L. Zheng, J. Lee, R. Venugopalan, O. Hen, D. Higinbotham, E. Aschenauer, T. Ullrich, *Probing short-range correlations in the deuteron via incoherent diffractive J/psi production with spectator tagging at the EIC*, Phys. Lett. B 811 (2020) 135877 | The incoherent, tagged partner of the coherent channel -- what a tagged diffractive vector-meson event looks like at the EIC, and the BeAGLE precedent for `CoherentSampler` plus spectator tagging. |
| `2006.03033.pdf` | W. Cosyn, C. Weiss, *Polarized electron-deuteron deep-inelastic scattering with spectator nucleon tagging*, Phys. Rev. C 102 (2020) 065204 | 52 pages and 19 figures: the light-front polarized spectral function in closed form. |
| `2006.06206.pdf` | H. Xing, C. Zhang, J. Zhou, Y. Zhou, *The cos 2phi azimuthal asymmetry in rho0 meson production in ultraperipheral heavy ion collisions*, JHEP 10 (2020) 064 | The dipole-model calculation of the same linear-polarization cos 2phi -- the background amplitude `phase_C_numbers.md` section C2.5 bounds. |
| `2006.12099.pdf` | W. Zha, J. Daniel Brandenburg, L. Ruan, Z. Tang, Z. Xu, *Exploring the double-slit interference with linearly polarized photons*, Phys. Rev. D 103 (2021) 033007 | The interference reading of the same modulation -- the alternative that makes the two-nucleus part of the STAR result inapplicable to e+A, which is what keeps the O5 background estimate honest. |
| `2008.02895.pdf` | T. Liu, W. Melnitchouk, J. Qiu, N. Sato, *Factorized approach to radiative corrections for inelastic lepton-hadron collisions*, Phys. Rev. D 104 (2021) 094033 | Whether a fixed-target O(alpha) tail ports to collider kinematics at all -- the question that decides whether the tree's POLRAD transcription is the right object for the EIC. |
| `2008.11437.pdf` | S. Fucini, S. Scopetta, M. Viviani, *Incoherent deeply virtual Compton scattering off 4He*, Phys. Rev. C 102 (2020) 065205 | Its incoherent companion, and the paper that names the Orsay-Perugia generator in text. |
| `2011.02464.pdf` | H. Mantysaari, K. Roy, F. Salazar, B. Schenke, *Gluon imaging using azimuthal correlations in diffractive scattering at the Electron-Ion Collider*, Phys. Rev. D 103 (2021) 094026 | The unpolarized e-V azimuthal modulation -- a third cos-type background under `cos2phi_coefficient`, and the size any tensor signal has to beat. |
| `2011.08583.pdf` | S. Kumano, Q. Song, *Transverse-momentum-dependent parton distribution functions up to twist 4 for spin-1 hadrons*, Phys. Rev. D 103 (2021) 014025 | The complete spin-1 correlator decomposition -- 40 TMDs over twists 2-4 -- i.e. the map of everything a spin-1 target can carry, against which the tree implements a handful. |
| `2012.09970.pdf` | A. Afanasev et al., *CFNS Ad-Hoc meeting on Radiative Corrections Whitepaper*, CFNS ad-hoc workshop whitepaper, arXiv:2012.09970 (2020), unpublished | The EIC-specific radiative problem statement, including nuclear targets -- the document that says what the community thinks is unsolved. |
| `2106.13466.pdf` | Y. Hagiwara, C. Zhang, J. Zhou, Y. Zhou, *Probing the gluon tomography in photoproduction of di-pions*, Phys. Rev. D 104 (2021) 094021 | cos 4phi and the elliptic-gluon vs QED-radiation splitting -- the higher harmonic a rho0 control measurement has to survive. |
| `2106.15849.pdf` | S. Kumano, Q. Song, *Twist-2 relation and sum rule for tensor-polarized parton distribution functions of spin-1 hadrons*, JHEP 09 (2021) 141 | A second free, closed-form sum rule beside Close-Kumano -- integral dx f2LT = 0 with f2LT = (2/3)f_LT - f1LL, plus a Wandzura-Wilczek-type twist-2 relation. |
| `2112.12462.pdf` | K. J. Eskola, P. Paakkinen, H. Paukkunen, C. A. Salgado, *EPPS21: A global QCD analysis of nuclear PDFs*, Eur. Phys. J. C 82 (2022) 413 | `Epps21Ratio` / `EmcBaseline::Epps21` and the grid `EPPS21nlo_CT18Anlo_Li6` -- the shipped unpolarized 6Li EMC baseline (`D-5`), whose proton baseline is the CT18A fit of section 2. |
| `2201.12363.pdf` | R. Abdul Khalek, R. Gauld, T. Giani, E. R. Nocera et al., *nNNPDF3.0: Evidence for a modified partonic structure in heavy nuclei*, Eur. Phys. J. C 82 (2022) 507 | The third fit of the same spread; `nNNPDF30_nlo_as_0118_A6_Z3` is its 6Li member. |
| `2202.12200.pdf` | V. Guzey, M. Rinaldi, S. Scopetta, M. Strikman, M. Viviani, *Coherent J/psi electroproduction on 4He and 3He at the Electron-Ion Collider: probing nuclear shadowing one nucleon at a time*, Phys. Rev. Lett. 129 (2022) 242503 | The closest published thing to a 6Li coherent amplitude: a Gribov-Glauber expansion on 4He and 3He with one-, two- and three-body form factors built from realistic AV18 few-body wave functions -- the same class the tree already carries. |
| `2203.11601.pdf` | C. Bierlich, S. Chakraborty, N. Desai, L. Gellersen, I. Helenius, P. Ilten, L. Lonnblad, S. Mrenna, S. Prestel, C. T. Preuss, T. Sjostrand, P. Skands, M. Utheim, R. Verheyen, *A comprehensive guide to the physics and usage of PYTHIA 8.3*, SciPost Phys. Codebases 8 (2022) | The hadronizer itself and the manual for everything in `src/pythia/` and `docs/PYTHIA_BRIDGE.md`; gate D-1. |
| `2204.01625.pdf` | M. S. Abdallah et al. (STAR Collaboration), *Tomography of Ultra-relativistic Nuclei with Polarized Photon-gluon Collisions*, Sci. Adv. 9 (2023) eabq3903 | The measured size of the photon-polarization cos 2phi background -- mechanism (ii) in `physics_literature.md` section 7's warning, and the reason `tensor_flip_plan` exists at all. |
| `2204.11998.pdf` | W. Chang, E. Aschenauer, M. D. Baker, A. Jentsch, J. Lee, Z. Tu, Z. Yin, L. Zheng et al., *BeAGLE: Benchmark eA Generator for LEptoproduction in high energy lepton-nucleus collisions*, (BeAGLE), Phys. Rev. D 106 (2022) 012007 | The standard eA generator for the EIC (DPMJet-III + PyQM + FLUKA evaporation on a nuclear spectral function). |
| `2207.03712.pdf` | H. Mantysaari, F. Salazar, B. Schenke, *Nuclear geometry at high energy from exclusive vector meson production*, Phys. Rev. D 106 (2022) 074019 | Strong-interaction radius against charge radius -- a named systematic on `gaussian_slope(r_rms_fm)`, which takes the charge radius and uses it as a gluonic one. |
| `2209.12161.pdf` | D. Fu, B. Sun, Y. Dong, *Generalized parton distributions of spin-3/2 particles*, Phys. Rev. D 106 (2022) 116012 | `[Fu22]` -- the J = 3/2 forward-limit structure functions, i.e. the b1 definition 7Li would be normalised against. |
| `2210.03785.pdf` | D. Byer, V. Khachatryan, H. Gao, I. Akushevich, A. Ilyichev et al., *SIDIS-RC EvGen: a Monte-Carlo event generator of semi-inclusive deep inelastic scattering with the lowest-order QED radiative corrections*, Comput. Phys. Commun. 287 (2023) 108702 | The maintained modern implementation of the same O(alpha) formalism -- the "Byer" reference, and the living code to compare a transcription against. |
| `2212.04730.pdf` | A. Afanasev, I. Akushevich, A. Ilyichev, N. Merenkov, *ESFRAD. FORTRAN code for calculation of QED corrections to polarized ep-scattering by the electron structure function method*, arXiv:2212.04730 (2022), unpublished preprint | Higher-order QED against POLRAD and RADGEN -- the published measure of the exponentiation `TPeakPlusLL` does not have, and the bound on how much that costs. |
| `2303.04866.pdf` | H. Mantysaari, B. Schenke, C. Shen, W. Zhao, *Multi-scale Imaging of Nuclear Deformation at the Electron Ion Collider*, Phys. Rev. Lett. 131 (2023) 062301 | Deformation -> coherent/incoherent ratio: the unpolarized half of the `delta_b_m` geometry, and the precedent for reading a shape off a diffractive pattern. |
| `2306.14578.pdf` | A. Afanasev, J. C. Bernauer, P. Blunden, J. Blumlein et al., *Radiative Corrections: From Medium to High Energy Experiments*, (topical review), Eur. Phys. J. A 60 (2024) 91 | The modern baseline for "what 1.5 % on A_zz means", and the survey of which radiative codes are still alive. |
| `2310.13211.pdf` | S. Gardiner, J. Isaacson, L. Pickering, *NuHepMC: A standardized event record format for neutrino event generators*, SciPost Phys. Codebases 57 (2025) | The template for a community attribute convention layered on HepMC3 -- the published precedent for the ion-spin attribute convention (`F-4`, `CH-18`). |
| `2402.11561.pdf` | D. Fu, Y. Dong, S. Kumano, *Transversity generalized parton distributions in spin-3/2 particles*, Phys. Rev. D 109 (2024) 096006 | 16 transversity GPDs per parton -- the gluon-transversity analogue of Delta for J = 3/2. |
| `2404.15984.pdf` | C. Bierlich et al., *Robust Independent Validation of Experiment and Theory: Rivet version 4 release note*, SciPost Phys. Codebases 36 (2024) | Rivet 4.1.2 ships in `eic_xl-nightly.sif`; it is the route to the 46 H1/ZEUS DIS analyses behind `CH-1`...`CH-5` -- the cheapest external validation surface the chain has. |
| `2406.01180.pdf` | S. Kumano, *Parton distribution functions and fragmentation functions of spin-1 hadrons*, Eur. Phys. J. A 60 (2024) 205 | The single best entry point to the whole spin-1 sector -- the review to cite once instead of six primaries, when a primary is not what is wanted. |
| `2410.20950.pdf` | I. Helenius, J. O. Laulainen, C. T. Preuss, *Multi-Jet Production in Deep Inelastic Scattering with Pythia*, JHEP 05 (2025) 153 | PYTHIA's DIS shower validated against H1 data -- the published pedigree gate `D-1` borrows rather than re-deriving. |
| `2502.20044.pdf` | J. Poudel, A. Bacchetta, J. Chen, D. Keller, I. Fernando, E. Long, D. Ruth, N. Santiesteban, K. Slifer, *Spin 1 Transverse Momentum Dependent Tensor Structure Functions in CLAS12*, arXiv:2502.20044 (2025), unpublished preprint | The next b1 datum: an RG-C ND3 analysis extracting b1 inclusively and F_U(LL),T and F^{cos2phi}_U(LL) in SIDIS. |
| `2504.21177.pdf` | M. M. Dalton, A. Deur, C. Keith, *Potential for Tensor Polarized Deuterons in Hall D at Jefferson Lab*, Eur. Phys. J. A 61 (2025) 111 | The only other proposal anywhere to measure a tensor observable in coherent vector-meson production, and the m = 0 / intermediate-\|t\| sensitivity argument behind `a2_m_state`. |
| `2505.23487.pdf` | J. Cammarota, J. Qiu, K. Watanabe, J. Zhang, *Factorized QED and QCD Contribution to Deeply Inelastic Scattering*, Phys. Rev. D 112 (2025) 056007 | The successor framework to the entry above -- NLO joint QED x QCD factorization. |
| `2508.06134.pdf` | J. Zhao, A. Bacchetta, S. Kumano, T. Liu, Y. Zhou, *Semi-inclusive deep inelastic scattering off a tensor-polarized spin-1 target*, JHEP 12 (2025) 067 | 23 structure functions, 21 non-vanishing at tree level through twist 3 -- the SIDIS counterpart of the inclusive decomposition the tree implements, and the map of what a SIDIS extension would owe. |
| `2602.11587.pdf` | D. Fu, Y. Dong, S. Kumano, J. Xie, *Generalizing the Soffer Bound: Positivity Constraints on Parton Distributions of Spin-3/2 Particles*, Phys. Rev. D 113 (2026) L111901 | The complete set of spin-3/2 positivity bounds, for the first time. |
| `2605.00454.pdf` | H. Mantysaari, H. Roch, B. Schenke, C. Shen, W. Zhao, *Nuclear structure and saturation effects from diffractive vector meson production*, Phys. Rev. D 114 (2026) 014068 | Retires the tree's "lightest published nucleus is Ca" note: it generates initial-state nucleon configurations from VMC, PGCM (with clustering and uniform sampling), NLEFT and GFMC (used for 3He and 4He) inside a Good-Walker amplitude. |
| `2606.14633.pdf` | T. Toll, D. Ghosh, A. Srivastav, *Efficient calculation of exclusive diffractive cross sections at the EIC and LHeC with the Sartre event generator*, arXiv:2606.14633 (2026), unpublished preprint | The table-speedup paper -- the cost half of the "Sartre route" decision, and the paper whose abstract prompted (and whose source code then disproved) the Sartre-feasibility claim. |
| `2607.09237.pdf` | S. Kumano, K. Kuroki, *Tensor-polarized parton distribution functions of the deuteron by a convolution model*, arXiv:2607.09237 (2026), unpublished preprint | An independent second convolution calculation at exactly Q2 = 2.5 GeV2 -- the same point as `tables::kB1CdksQ2p5`, which makes it a directly comparable external check rather than another camp. |
| `2608.23445.pdf` | A. Mondal, A. Kumar, D. Sarkar, *Imprints of nuclear shell structure in exclusive vector meson production*, arXiv:2608.23445 (2026), unpublished preprint | How far a realistic light-nucleus density moves the coherent \|t\| shape -- a named systematic on `slope_b` beyond the Gaussian. |
| `JINR-P2-10113.pdf` | D. Yu. Bardin, N. M. Shumeiko, *An Exact Calculation of the Lowest Order Electromagnetic Correction to the Elastic Scattering*, JINR preprint P2-10113 (Dubna, September 1976), in Russian; English text published as Nucl. Phys. B 127 (1977) 242 | The free JINR preprint of the Bardin-Shumeiko lowest-order electromagnetic radiative-correction calculation to elastic scattering; the covariant infrared separation POLRAD Eq. (18) and `polrad_tpeak_quadrature` are built on. |
| `JLab_PR12-13-011.pdf` | K. Slifer et al., *The Deuteron Tensor Structure Function b1*, JLab proposal PR12-13-011 to PAC 40 (2013), update to PR12-11-110 | The approved b1 experiment: 0.16 < x < 0.49, 0.8 < Q2 < 5.0 GeV2, 30 days at 11 GeV, Pzz = 20% (read from the PDF). |
| `JLab_PR12-15-005.pdf` | E. Long, K. Slifer, P. Solvignon et al., *Measurements of the Quasi-Elastic and Elastic Deuteron Tensor Asymmetries*, JLab proposal PR12-15-005 to PAC 43 (2015), update to LOI12-14-002 | The Azz proposal: quasi-elastic Azz at x > 1 plus elastic T20 over 0.2 < Q2 < 1.8 GeV2 (read from the PDF). |
| `JPCS_543_012003.pdf` | K. Slifer, *The Deuteron Polarized Tensor Structure Function b1*, J. Phys. Conf. Ser. 543 (2014) 012003 (IOP, open access CC-BY) | The four-page public, citable summary of proposal PR12-13-011 -- the version to cite in place of the proposal PDF itself. |
| `JPCS_543_012008.pdf` | N. Kalantarians, *Tensor Polarized Deuteron at an Electron-Ion Collider*, J. Phys. Conf. Ser. 543 (2014) 012008 (IOP) | The earliest written case for tensor polarization at an EIC -- the direct ancestor of this project's premise. |
| `JPCS_543_012010.pdf` | E. Long, *Potential for a Tensor Asymmetry Azz Measurement in the x > 1 Region at Jefferson Lab*, J. Phys. Conf. Ser. 543 (2014) 012010 (IOP) | The citable version of the Azz quasi-elastic case; it names Frankfurt-Strikman 1988 as the first quasi-elastic Azz calculation. |
| `SLAC-PUB-0380.pdf` | L. W. Mo, Y. S. Tsai, *Radiative Corrections to Elastic and Inelastic ep and mu p Scattering*, SLAC-PUB-380 (January 1968); published as Rev. Mod. Phys. 41 (1969) 205 | [MT69] -- Section III's exact elastic radiative tail, Appendix B's exact bremsstrahlung expression, and the paper's own exact-vs-peaking comparison. |
| `SLAC-PUB-0848.pdf` | Y. S. Tsai, *Radiative Corrections to Electron Scatterings*, SLAC-PUB-848 (1971) (INSPIRE recid 67278, free preprint) | The same Mo-Tsai formalism restated at greater length, and what most of the field actually cites for the exact radiative-correction formulas. |
| `SLAC-PUB-1528.pdf` | S. Stein et al., *Electron Scattering at 4 Degrees with Energies of 4.5 GeV - 20 GeV*, SLAC-PUB-1528 (January 1975); published as Phys. Rev. D 12 (1975) 1884 | Reference [12] of the HERMES b1 paper: the inclusive elastic-plus-quasi-elastic radiative-tail recipe that HERMES's own radiative correction is built on. |
| `hep-ex_0606004.pdf` | A. Aktas et al. (H1 Collaboration), *Measurement and QCD Analysis of the Diffractive Deep-Inelastic Scattering Cross Section at HERA (the H1 2006 DPDF Fit A/B)*, Eur. Phys. J. C 48 (2006) 715 | `[H1FitB]` -- the diffractive PDF family behind `PDF:PomSet = 6`, i.e. the definition of the T2 Pomeron tier (`F-8`). |
| `hep-ph_0007120.pdf` | A. Bacchetta, P. J. Mulders, *Deep inelastic leptoproduction of spin-one hadrons*, Phys. Rev. D 62 (2000) 114004 | The spin-1 TMD basis (f1LL, hperp1LL, ...) and the Trento-convention names the modern spin-1 literature converts into -- the dictionary between `TensorSF`'s labels and everyone else's. |
| `hep-ph_0102086.pdf` | A. Afanasev, I. Akushevich, N. Merenkov, *Model independent radiative corrections in processes of polarized electron-nucleon elastic scattering (MASCARAD)*, Phys. Rev. D 64 (2001) 113009 | Polarized elastic radiative corrections with realistic acceptance -- the acceptance question buried inside POLRAD Eq. (44). |
| `hep-ph_0105032.pdf` | A. V. Afanas'ev, I. Akushevich, G. I. Gakh, N. P. Merenkov, *Radiative Corrections to Polarized Inelastic Scattering in Coincidence*, JETP 93 (2001) 449 | Radiative corrections in coincidence, polarized and model-independent -- the calculation against which the tagged channels' `rc_tail == 1` (the one place the generator asserts a correction is exactly 1.0) should be priced. |
| `hep-ph_0106180.pdf` | I. Akushevich, A. Ilyichev, N. Shumeiko, *Radiative effects in scattering of polarized leptons by polarized nucleons and light nuclei*, review, arXiv:hep-ph/0106180 (2001), unpublished | The code map -- POLRAD, RADGEN, HAPRAD, DIFFRAD, MASCARAD -- and the covariant formulae in one place; the document that says which of these codes does what. |
| `hep-ph_0106192.pdf` | E. R. Berger, F. Cano, M. Diehl, B. Pire, *Generalized parton distributions in the deuteron*, Phys. Rev. Lett. 87 (2001) 142302 | The spin-1 GPD basis -- 5 quark and 9 gluon distributions, including the gluon-transversity GPD -- i.e. the exclusive counterpart of the coherent 6Li channel and the formal home of the object `CoherentScenario::amp` parameterises. |
| `hep-ph_0109068.pdf` | E. Boos et al., *Generic User Process Interface for Event Generators (Les Houches Accord 1)*, arXiv:hep-ph/0109068 (2001), unpublished Les Houches Accord 1 note | Defines `HEPRUP`/`HEPEUP` and the strategy-3 weighting convention the bridge uses -- the normative document for the one interface LiPolGen must get exactly right. |
| `hep-ph_0304189.pdf` | H. Kowalski, D. Teaney, *An Impact Parameter Dipole Saturation Model*, Phys. Rev. D 68 (2003) 114005 | IPsat/bSat itself -- Sartre's `bSat`, and the b-dependence any 6Li amplitude would inherit. |
| `hep-ph_0403262.pdf` | G. I. Gakh, O. Shekhovtsova, *Radiative corrections to deep-inelastic ed- scattering. Case of tensor polarized deuteron*, JETP 99 (2004) 898 | The source of `RC_DELTA_LOW_X` = 0.30 and of the shape the tree examined and rejected -- one of only two tensor-target radiative calculations that exist. |
| `hep-ph_0507286.pdf` | L. Frankfurt, M. Strikman, C. Weiss, *Small-x Physics: From HERA to LHC and beyond*, Ann. Rev. Nucl. Part. Sci. 55 (2005) 403 | The review behind the framing that coherent diffraction images transverse gluon structure. |
| `hep-ph_0606272.pdf` | H. Kowalski, L. Motyka, G. Watt, *Exclusive diffractive processes at HERA within the dipole picture*, Phys. Rev. D 74 (2006) 074016 | The bSat fit Sartre cites, with the vector-meson wave functions and t-slopes. |
| `hep-ph_0609017.pdf` | J. Alwall et al., *A standard format for Les Houches Event Files*, Comput. Phys. Commun. 176 (2007) 300 | The LHEF record `LhaupDis::setEvent` writes (`PYTHIA_BRIDGE.md` section 6). |
| `hep-ph_9302320.pdf` | S. Kumano, *Tensor Structure Function b1(x) For Spin-One Hadrons*, MKPH-T-93-03 (1993), unpublished preprint | The earliest statement of the b1-b4 measurement case and the historical head of the Kumano series that runs through section 6c -- useful for dating the field's claims, not for numbers. |
| `hep-ph_9605286.pdf` | G. Ingelman, A. Edin, J. Rathsman, *LEPTO 6.5 - A Monte Carlo Generator for Deep Inelastic Lepton-Nucleon Scattering*, Comput. Phys. Commun. 101 (1997) 108 | The DIS hard-process + Lund-fragmentation layer under both DJANGOH and PEPSI, and the closest published precedent for LiPolGen's own architecture (parton-level cross section handed to a string hadronizer). |
| `hep-ph_9605291.pdf` | A. Yu. Umnikov, *Relativistic Calculation of Structure Functions b1,2(x) of the Deuteron*, Phys. Lett. B 391 (1997) 177 | The validity floor on `b1_convolution`: a non-relativistic convolution gets small x wrong and violates the exact sum rules. |
| `hep-ph_9611460.pdf` | N. N. Nikolaev, W. Schaefer, *Nonvanishing tensor polarization of sea quarks in polarized deuterons*, Phys. Lett. B 398 (1997) 245; Erratum PLB 407 (1997) 453 | The shadowing camp, and the reason the small-x end of the tree's b1 is not a settled number: `01_tensor-b1.md` section 3.1 records that it predicts b2 rising at small x and A2 ~ 1 %, two orders of magnitude above the impulse approximation, with an explicit claim that the Close-Kumano sum rule is broken. |
| `hep-ph_9706516.pdf` | I. Akushevich, A. Ilyichev, N. Shumeiko, A. Soroko, A. Tolkachev, *POLRAD 2.0. FORTRAN code for the Radiative Corrections Calculation to Deep Inelastic Scattering of Polarized Particles*, Comput. Phys. Commun. 104 (1997) 201 | The single most load-bearing reference in the radiative sector: `01_radiative.md` Table 1 lists Eqs. |
| `hep-ph_9709455.pdf` | J. Edelmann, G. Piller, W. Weise, *Deuteron Spin Structure Functions at Small Bjorken-x*, Phys. Rev. C 57 (1998) 3392 | Its long version -- spin-extended Glauber-Gribov multiple scattering worked out in full. |
| `hep-ph_9711323.pdf` | K. Bora, R. L. Jaffe, *The Double Scattering Contribution to b1(x,Q2) in the Deuteron*, Phys. Rev. D 57 (1998) 6906 | The VMD/double-scattering camp, with the opposite small-x verdict to Nikolaev-Schafer: b1 -> 0 as x -> 0 by rotational symmetry. |
| `hep-ph_9906408.pdf` | I. Akushevich, H. Boettcher, D. Ryckbosch, *RADGEN 1.0. Monte Carlo Generator for Radiative Events in DIS on Polarized and Unpolarized Targets*, arXiv:hep-ph/9906408 (1999), unpublished preprint | The POLRAD-2.0-derived radiative generator used by HERMES and COMPASS, and the code that physically ships inside both BeAGLE and PEPSI. |
| `jung1995_rapgap_desy93-182.pdf` | H. Jung, *Hard Diffractive Scattering in High Energy ep Collisions and the Monte Carlo Generator RAPGAP*, Comput. Phys. Commun. 86 (1995) 147; preprint DESY 93-182 | The diffractive-DIS generator behind the official ePIC `DDIS/rapgap3.310` samples, and the closest external analogue of the T2 gamma*-Pomeron tier's hadronic final state (`CH-9`). |
| `nucl-ex_9809002.pdf` | Z.-L. Zhou, M. Bouwhuis, M. Ferro-Luzzi, E. Passchier, R. Alarcon, M. Anghinolfi, H. Arenhoevel et al., *Tensor Analyzing Powers for Quasi-Elastic Electron Scattering from Deuterium*, Phys. Rev. Lett. 82 (1999) 687 | The sole citation under `RcOptions::qe_tensor_scale = 0.0`. |
| `nucl-th_0002058.pdf` | L. C. Maximon, J. A. Tjon, *Radiative Corrections to Electron-Proton Scattering*, Phys. Rev. C 62 (2000) 054320 | How much Mo-Tsai's soft-photon/peaking treatment leaves out -- i.e. a published bound on the error `RcTailModel::TPeak` carries. |
| `nucl-th_9701026.pdf` | J. Edelmann, G. Piller, W. Weise, *Polarized deuteron structure functions at small x*, Z. Phys. A 357 (1997) 129 | A third shadowing calculation: b1 "surprisingly large at x < 0.1", dominated by coherent double scattering. |

`REFERENCES.md` section 1a #16's Bardin–Shumeiko 1977 JINR preprint
(`JINR-P2-10113.pdf`, Russian-language, INSPIRE recid 111694) and its
`used_in`/companion notes for `JLab_PR12-13-011.pdf` are cross-referenced
against the pre-existing corpus entries `jlab-pr12-13-011-b1-and-2023-jeopardy`
(`file: null`) rather than merging into them, per the append-only rule on
this dictionary.

Separately, `REFERENCES.md` section 10.1/10.4 listed **13** references as
`free, not fetched` — open access and reachable, but no PDF had been
committed (12 bullets in sections 7a/8, plus item 27 in section 1, found
free in the 2026-09-16 reconciliation pass). Each was tried once with a
plain fetch (arXiv, an INSPIRE `documents` URL, or the publisher's
open-access link). Six resolved to a real PDF (`%PDF` magic, > 20 kB, and
a first page carrying the title) and are now on disk:

| file | reference | used for / verified |
|---|---|---|
| `FERMILAB-PUB-03-339.pdf` | S. Agostinelli et al. (GEANT4 Collaboration), *GEANT4 - A Simulation Toolkit*, Nucl. Instrum. Meth. A 506 (2003) 250-303, doi:10.1016/S0168-9002(03)01368-8; free preprint FERMILAB-Pub-03/339 (INSPIRE recid 593382) | The transport engine `npsim` runs (`F-1`). |
| `IAEA-INDC-NDS-0794.pdf` | N. J. Stone, *Table of recommended nuclear magnetic dipole moments: Part I, long-lived states*, IAEA INDC(NDS)-0794 (November 2019) | The free IAEA twin of the paywalled ADNDT table (REFERENCES.md section 1 item 42), fetched here alongside the electric-quadrupole twin above. |
| `IAEA-INDC-NDS-0833.pdf` | N. J. Stone, *Table of nuclear electric quadrupole moments*, IAEA INDC(NDS)-0833 (October 2021), doi:10.61092/iaea.a6te-dg7q | The free IAEA twin of the paywalled ADNDT table (REFERENCES.md section 1 item 42). |
| `JINR-E1-12727.pdf` | D. Albrecht, M. Csatlos, J. Ero, Z. Fodor et al., *Large-angle quasi-free scattering in 6Li(p,pd)4He at 670 MeV*, JINR preprint E1-12727 (Dubna, 1979); published as Nucl. Phys. A338 (1980) 477, doi:10.1016/0375-9474(80)90045-7 | A hadronic probe of the same alpha-d distribution the tagged channel samples, and where the factor-two width disagreement (FWHM ~ 70 vs 120 MeV/c) lives -- the honest external band on the distribution's width. |
| `JPCS_513_022010.pdf` | M. Frank, F. Gaede, C. Grefe, P. Mato, *DD4hep: A Detector Description Toolkit for High Energy Physics Experiments*, J. Phys. Conf. Ser. 513 (2014) 022010, doi:10.1088/1742-6596/513/2/022010 (IOP, open access) | The geometry layer under `npsim` / `eic/epic` -- the first thing downstream of the HepMC3 file (`F-1`, `CH-13`). |
| `SLAC-PUB-1365.pdf` | Y. S. Tsai, *Pair Production and Bremsstrahlung of Charged Leptons*, Rev. Mod. Phys. 46 (1974) 815; Erratum Rev. Mod. Phys. 49 (1977) 421; preprint SLAC-PUB-1365 | The equivalent-radiator whose single-z form `ll_radiator` implements. |

`IAEA-INDC-NDS-0833.pdf` and `IAEA-INDC-NDS-0794.pdf` (the Stone
quadrupole/dipole-moment tables) and `JINR-E1-12727.pdf` (the Albrecht
*et al.* preprint) are, likewise, new file-bearing entries standing beside
the pre-existing `file: null` corpus entries `stone2021-quadrupole-moments`,
`stone2019-magnetic-dipole-moments` and `albrecht1980-6li-ppd-670mev`
rather than edits to them.

The remaining eight of the 13 were not added, for three different reasons:

- **Too large to commit**, confirmed reachable by a plain fetch (no bot
  challenge) but not kept, consistent with this corpus's existing
  size policy for the same documents (see the 2026-09-06 section above):
  the **EIC Conceptual Design Report** — <https://www.osti.gov/servlets/purl/1765663>
  (a complete, valid 172 MB PDF; the corpus's existing `file: null` entry
  already notes 178 MB and PLEASE-DOWNLOAD-scale content); the **ePIC
  Preliminary Design Report** — <https://doi.org/10.5281/zenodo.13866213>
  (114 MB, not re-tried this pass, same size policy); the **ePIC
  Preliminary TDR v3.1** — <https://doi.org/10.5281/zenodo.18271601>
  (260 MB, not re-tried, same policy); and the **PDG Review of Particle
  Physics** — <https://doi.org/10.1103/PhysRevD.110.030001> (a complete,
  valid 161 MB / 2382-page PDF, fetched to confirm reachability then not
  committed — the corpus cites it for one convention only, the Monte
  Carlo particle numbering scheme, which does not need the whole book on
  disk).
- **Not a single document to fetch**: the **Nuclear Charge Density
  Archive** (Day *et al.*, UVA) — <https://discovery.phys.virginia.edu/research/groups/ncd/index.html>;
  the **Quasi-elastic Electron–Nucleus Scattering Archive** ⁶Li dataset —
  <https://discovery.phys.virginia.edu/research/groups/qes-archive/data/6Li.dat>;
  and the two **HEPData** records —
  <https://www.hepdata.net/record/ins394050> and
  <https://www.hepdata.net/record/ins393377>. These are data archives and
  machine-readable records, not papers with a PDF and a title page; the
  URLs are recorded here and in `refs_dict.json` for the user to fetch the
  data directly if a specific number is needed.
- **A volume, not a paper**: the *Proceedings, Tensor Polarized Solid
  Target Workshop* (JLab, 10–12 March 2014, INSPIRE recid 1324998) is a
  conference volume whose four individual papers this project actually
  uses are already separate, file-bearing entries in `docs/references/REFERENCES.md`
  sections 3a and 6f (`JPCS_543_012003.pdf`, `JPCS_543_012008.pdf`,
  `JPCS_543_012010.pdf`, and `1407.1653.pdf`'s sibling proceedings entry).
  There is no additional single PDF to add for the volume itself.

None of the eight above hit a bot challenge; all are simply out of scope
for a single committed PDF under this pass's size and single-document
criteria. Their exact URLs are recorded above (and, for the six that have a
`file: null` corpus entry -- all but the PDG review and the workshop volume --
in `refs_dict.json`) so
a future pass (or the user) can revisit the size question directly.
