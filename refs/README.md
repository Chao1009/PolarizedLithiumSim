# refs/ — reference papers consulted for the reconstruction-chain work

PDFs are kept locally only (`refs/*.pdf` is git-ignored, like the other
local paper copies); this index records what each one was used for and
which numbers were taken from it (2026-08-25 read-through).

**Machine-readable dictionary:** `refs/refs_dict.json` maps every local
file (and the external references the reports lean on) to identifiers,
title, key content with equation/figure/slide pointers, and where it is
used in the repository. Look things up with

```bash
python refs/find_ref.py slot            # keyword search over all fields
python refs/find_ref.py "Eq. (9)" a2    # several terms, all must match
python refs/find_ref.py --key maple     # full entry
python refs/find_ref.py --list
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
refs/; PDFs git-ignored):

| file | reference | used for |
|---|---|---|
| `EIC_Seminar_SMaple_2024.pdf` | S. Maple (ePIC), *Tracking and inclusive DIS reconstruction with the ePIC detector at the EIC*, seminar, Birmingham, 11 Dec 2024 (Indico Global 4787) | **The ePIC source for the hadronic-y resolution:** smeared EIC pseudodata (Djangoh 18×275, Q² > 1 GeV²) with σ(δ_h)/δ_h = 25%, σ(p_T,h)/p_T,h = 25%, σ(E_e)/E = 11%/√E ⊕ 2%, σ(θ_e) = 0.1 mrad (slide 44); the same 25% parametrization fitted to full ePIC simulation (Craterlake 23.12.0, Pythia8, Q² > 100 GeV²; slide 47); Δy/y per method in y bins — Σ/JB/DA widths 0.2–0.3 in 0.01 < y < 0.05, electron flat (slides 45, 48); the 18×275 "best method for y" map with 25/10/5/1% markers, DA at the y = 0.01 edge (slide 42); coverage 0.01 < y < 0.95, Q² > 1 GeV², five log bins per decade, "lower y accessible → easier to rely on overlap between data at different √s"; Bayesian kinematic fit incl. ISR energy (slides 43–56). |
| `2210.09048.pdf` | ATHENA Collaboration, JINST 17 (2022) P10019 (43 MB) | Sec. 3.1 / Fig. 22: "at the lowest y values, the electron method resolution degrades like 1/y ... e−Σ or Double Angle methods ... for y ≲ 0.1"; "the resolution with the JB method is at the 20–30% level throughout the kinematic range"; Fig. 22 marker sizes: ≈ 25% at y ≈ 0.01 (DA), ≈ 10% at low Q², y ≈ 0.02–0.1 (e−Σ); y > 0.01 cut motivated by reconstructibility. |
| `2110.05505v2.pdf` | M. Arratia, D. Britzger, O. Long, B. Nachman, NIM A 1025 (2022) 166164 | Table 1 with the formulas of all classical methods; Fig. 5 (ATHENA fast simulation, Q² > 200 GeV², 18×275): RMS(y)/y at y → 0.05 ≈ 0.3 (electron), 0.17 (JB), 0.15 (DA), 0.13 (IΣ); Sec. 5 on why Σ is fragile at low y (noise, acceptance) and the Delphes-vs-full-simulation comparison (Fig. 11). |
| `2206.04897.pdf` | R. Aggarwal, A. Caldwell, JINST 17 (2022) P09035 | The Bayesian kinematic-fit method (HERA kinematics, Q² > 400 GeV²) that the ePIC seminar applies. |
| `2209.14489.pdf` | C. Pecar, A. Vossen, DIS2022 proceedings | SIDIS reconstruction with a particle-flow network vs electron/JB/DA: HFS methods surpass the electron method at very low y. |

Not downloaded: M. Diefenthaler et al., EPJ C 82 (2022) 1064
(arXiv:2108.11638) — ZEUS-based DNN reconstruction; the ePIC Inclusive WG
wiki (open task on hadron treatment in JB/DA/(e)Σ), the BNL "DIS
Kinematics" wiki page and B. Schmookler's `JeffersonLab/dis-reconstruction`
repository (Yellow-Report kinematic maps) are recorded in
`refs_dict.json` as the places to look for newer ePIC numbers. What is
still not published anywhere we found: an ePIC full-simulation δy/y at
Q² ≈ 1–3 GeV², y ≈ 0.01, for e + light ions.
