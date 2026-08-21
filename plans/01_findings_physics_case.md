# Findings: Physics Case and Inputs for the Polarized ⁶,⁷Li Simulation

*Sources: `docs/Discussions.pptx` (slides from J. Maxwell/JLab and I. Cloët/ANL),
`docs/ecrp_2026_proposal.pdf` (C. Peng, ANL, NOFO DE-FOA-0003602), and a
fetch-verified literature sweep (2026-06-12). Citations below were verified
against primary sources unless marked (unverified).*

## 1. What the simulation must demonstrate

The ECRP proposal develops the **source** (polarized ⁶,⁷Li atomic beam →
EBIS/Tandem → EIC). The simulation effort builds the **physics case**: show
that polarized Li beams enable measurements no other configuration delivers,
quantified with figures of merit. The proposal explicitly defers one item to
simulation: *"Spectator tagging via forward Roman pots and the zero-degree
calorimeter offers an additional handle, though its achievable purity for
charged-α and neutral fragments must be quantified through BeAGLE-class
simulations within the broader EIC program."*

**Four confirmed literature gaps = first-publication opportunities:**
1. No b₁ prediction or measurement exists for any A > 2 nucleus.
2. No numerical EIC gluon-transversity (Δ) projection exists for *any* target.
3. No study of α-tagged DIS on ⁶Li (or ⁷Li) exists.
4. No Li-beam luminosity/polarization parameter set exists for the EIC.

## 2. The three flagship observables

### 2.1 Nuclear gluonometry — double-helicity-flip Δ(x,Q²)  [⁶Li spin-1; ⁷Li spin-3/2 also eligible]
- Jaffe–Manohar PLB 223:218 (1989); Hoodbhoy–Jaffe–Manohar NPB 312:571
  (1989). Photon helicity flip by 2 ⇒ no quark contribution at leading
  twist ⇒ **purely gluonic**; bound nucleons (spin-½) cannot contribute —
  "exotic glue". Requires J ≥ 1, **transverse** polarization, **unpolarized
  electrons**. Never measured.
- Master formula (slides p.5; target spin angle θ_m, projection λ_m,
  c_m = 3|λ_m|−2, a_m = ¼c_m(3cos²θ_m−1)):

  dσ/(dx dy dφ) = (2yα²/Q²)·[F₁ + ⅔a_m b₁ + (1−y)/(xy²)(F₂ + ⅔a_m b₂)
                   − (1−y)/y²·c_m sin²θ_m Δ(x,Q²) cos 2φ]

- **Expected size**: bag-model/Sather–Schmidt estimate ∫xΔdx = −0.012·α_s ≈
  O(10⁻³); Kumano–Song "few percent" Drell-Yan asymmetries assume Δ_Tg = Δg
  and are explicitly an overestimate. **Plan FOMs around A_cos2φ ~
  10⁻³–10⁻²** ⇒ needs ~100 fb⁻¹-class samples and tensor polarization.
- Lattice anchor: Detmold–Shanahan PRD 94:014507 (2016) — first moment of
  gluon transversity in the **φ meson**, clearly nonzero (~10σ bare) at
  m_π = 450 MeV. NPLQCD PRD 96:094512 (2017): unpolarized gluon structure
  of A ≤ 3 nuclei.
- Prior proposal art: JLab **LOI12-16-006** (arXiv:1803.11206; Maxwell,
  Detmold, Jaffe, Shanahan et al.) using ¹⁴N in NH₃, x < 0.3 — LOI only, no
  sensitivity projection. EIC advantage (slides p.8–9): cos 2φ modulation
  with vector transverse polarization, vastly larger kinematic reach, no
  out-of-plane detector limitation.
- EPIOS (arXiv:2510.10794): "the best chance for discovery may come from
  larger nuclei… spin-1 ⁶Li or ¹⁴N, or spin-3/2 ⁷Li or ²³Na."

### 2.2 Tensor structure function b₁  [⁶Li vs free d, spin-1]
- b₁ = ½Σe_q²[q⁰ − (q⁺¹+q⁻¹)/2]; vanishes for weakly bound S-wave pairs ⇒
  sensitive to D-state, pions, hidden color, tensor sea. Close–Kumano sum
  rule ∫b₁dx = 0 unless tensor-polarized sea.
- **HERMES (PRL 95:242001)**: A_zz ≤ 2×10⁻² over 0.01 < x < 0.45;
  b₁(x=0.012) = +0.112±0.055±0.028 (large at low x!), zero crossing at
  x ≈ 0.2; ∫b₁dx = (1.05±0.34±0.35)×10⁻² (~2σ tensor sea).
- Theory tension = science driver: standard convolution gives |b₁| < 10⁻³
  at x ≳ 0.2 (Cosyn–Dong–Kumano–**Sargsian** PRD 95:074036), an order below
  HERMES; Miller (PRC 89:045203) reproduces HERMES with pions + only
  **0.15% hidden-color** probability.
- JLab: E12-13-011 (b₁, Hall C, 0.16 < x < 0.49, 0.8 < Q² < 5; active per
  PAC51 jeopardy 2023); E12-15-005 is quasi-elastic A_zz (x = 0.8–1.75),
  not DIS. EIC adds low-x reach (where HERMES saw b₁ ~ 0.1) and Q² arm.
- **⁶Li**: no b₁ theory exists (gap #1). First model for FOM purposes
  (our inference, to be developed with theory colleagues):
  b₁(⁶Li)/nucleon ≈ ⅓·P_d(⁶Li)·b₁(d) ≈ 0.29·b₁(d) + α–d D-wave term.
  EPIOS explicitly motivates "b₁ of the deuteron, free or embedded in ⁶Li".
- Extraction systematics (Cosyn et al. arXiv:2410.12764, EPJ A 61:83):
  at collider kinematics the clean b₁ extraction prefers polarization
  **along the momentum transfer**, not the beam axis — build this into the
  simulation early (spin-direction systematics).

### 2.3 Polarized EMC effect — g₁ of ⁷Li (and ⁶Li)
- ΔR_A(x) = g₁^A/[P_p·g₁^p + P_n·g₁^n] (Cloët–Bentz–Thomas PLB 642:210,
  Eq. 23). CBT (NJL): polarized EMC ≈ **2× the unpolarized EMC** (~10–20%
  valence depletion). Counter-model Tronchin–Matevosyan–Thomas
  (PLB 783:247): polarized ≈ unpolarized. **A 5%-level measurement
  discriminates the two camps** — and they map onto mean-field vs SRC
  origins of the EMC effect.
- Gluon sector: Wang–Bentz–Cloët–Thomas (J.Phys.G 49:03LT01,
  arXiv:2109.03591): polarized *gluon* EMC larger than unpolarized; ⁷Li
  named "the most promising case" — EIC Q² lever arm via dg₁/dlnQ² is the
  unique handle.
- **Effective polarizations (verified)**:
  | nucleus | P_p | P_n | source |
  |---|---|---|---|
  | ⁷Li | **+0.866** | **−0.037** | QMC/VMC (Wiringa et al. 1309.3794, via E12-14-001) |
  | ⁶Li | d-cluster pol ≈ 0.87; valence pair carries >90% of spin | (same) | Schellingerhout PRC 48:2714 |
  | ³He | −0.028 | +0.86 | Bissey PRC 65:064317 |
  Convention warning: Cloët's slides quote P_p = P_n = 1/3 for ⁶Li — this
  is the per-nucleon-normalized dilution (2 polarized of 6 nucleons);
  whole-nucleus values are ≈ 0.87×0.93 ≈ 0.81 per cluster nucleon. Fix the
  convention in step 1.1 before quoting FOMs (factor 2.4 at stake).
- JLab benchmark to beat: **E12-14-001** (CLAS12, Brooks–Kuhn): ⁷LiD +
  ⁶LiH targets, 0.06 < x < 0.8, 1 < Q² < 15 GeV², point-to-point syst
  3.2–4.2% + 4% scale, sensitive to depletion > 5%. EIC adds x < 0.06,
  Q² arm, and collider kinematics (no target-fragment absorption).
- COMPASS ⁶LiD prior art: P_d ≈ +0.56/−0.53, dilution f ≈ 0.4–0.5.

### 2.4 Spectator-tagged DIS
- Formalism: Cosyn–Weiss PLB 799:135035 (method) + PRC 102:065204
  (formalism); FSI: PRC 97:035209. Pole extrapolation in spectator pT at
  fixed light-cone fraction → free-nucleon structure. Tagged *tensor*
  asymmetries are O(1) at p_s ~ 300 MeV/c (D-wave dominance). New
  (Mar 2026): arXiv:2603.23699/23700 extend to SIDIS on spin-1 with tagging.
- Detector feasibility template: Jentsch–Tu–Weiss PRC 104:065205
  (arXiv:2108.08314) — e+d tagging with full far-forward simulation.
- FOM-study template: **Friščić et al. PLB 823:136726 (arXiv:2106.08805)**
  — A₁ⁿ from e+³He with **double spectator-proton tagging** (5×41 GeV,
  P = 0.7, 100 fb⁻¹; transverse case survives at 10 fb⁻¹): tagged method
  beats effective-polarization unfolding ×2–3 (×10 at low x). Replicate
  with α-tagged ⁶Li → embedded polarized deuteron.
- JLab methodology antecedents: BAND (E12-11-003A) and LAD (E12-11-107)
  backward-spectator double ratios vs virtuality.
- Far-forward routing for Li fragments (verified vs detector windows; see
  plans/03 §2.2): ⁷Li α-tag → Roman Pots (works at IP6); ⁶Li α-tag is
  beam-blind below the RP pT cutoff (R = 1.0); tritons have **no IP6
  coverage** (R > 1). This asymmetry between isotopes is itself a finding —
  the tagging program leads with ⁷Li at IP6; ⁶Li α-tagging favors the
  IR-8/second-detector secondary focus.

## 3. Machine and beam inputs (EPIOS-verified)

EPIOS white paper: Atoian et al., **arXiv:2510.10794**, PRC 113:060501
(2026). Species: d, ³He, ⁶Li, ⁷Li (+ B, Na source development).

| | d | ³He | ⁶Li | ⁷Li |
|---|---|---|---|---|
| G-factor | −0.143 | −4.184 | **−0.178** | **+1.532** |
| linear resonances to top | 63 | 2457 | **81** | **573** |
| snake solution | 15 Tm solenoid partial snake + jump quads | 6 full snakes | as d | partial snakes (ν_s 0.2–0.8) + jump quads |
| top energy | ~137 GeV/u | ~183 GeV/u | **~138 GeV/u** | **~117 GeV/u** |

- ⁶Li is spin-dynamically deuteron-like (small G): few resonances, partial
  solenoid snake; in the AGS only imperfection resonances. Deuteron
  antecedent: Huang et al. PRAB 23:021001 (tune jumps to >100 GeV/u).
- **No Li luminosity number exists** (gap #4) — FOMs must state a scaling
  assumption (we use 10 fb⁻¹/nucleon baseline, quoted ∈ {1,10,100}).
- Li polarimetry concept: Li–Li elastic CNI vs polarized Li jet (HJET
  analog) + Breit–Rabi absolute; analyzing-power theory needs work.
- Source (ANL + UKY, Sec. V.4): oven + de Laval nozzle, optical pumping
  ~671 nm, RF transitions; target ~100 µA (2–3 orders above the historic
  FSU source); siting option BNL Tandem. ECRP deliverables: P_z ≥ 0.90,
  P_zz ≥ 0.80, ≥10¹³ Li⁰/s into EBIS acceptance.
- **Timeline hooks**: EIC CD-4 ≈ 2034; INT program "Towards Realizing the
  Program with Polarized Ion Beams at EIC" **March 22 – April 2, 2027** —
  natural target venue for Phase-1 results.

## 4. Source-side constraints carried into the physics study

- ⁷Li is commissioned first → physics projections cover both isotopes with
  ⁷Li (polarized EMC + α-tag) deliverable first. Conveniently, ⁷Li is also
  the isotope whose tagging works at IP6 (§2.4).
- Polarization through EBIS charge-breeding is "a goal of the study rather
  than a promised outcome" → quote FOMs vs P ∈ {0.5, 0.7, 0.9}, never at a
  single point.
- Tensor polarization P_zz is the scarcer resource → b₁/Δ FOMs shown vs
  P_zz; spin-direction (along q vs along beam) systematics from 2410.12764.

## 5. Benchmarks the FOM study must beat (verified numbers)

| observable | existing best | EIC must show |
|---|---|---|
| b₁(d) | HERMES: δb₁ ~ (1–5)×10⁻² stat, 0.01<x<0.45; E12-13-011: high-x precision | low-x reach where b₁ ~ 0.1, zero-crossing mapping, ⁶Li-vs-d difference |
| Δ(x,Q²) | nothing (LOI only) | first projection; 5σ reach vs Δ/F₁ ∈ 10⁻³–10⁻² |
| ΔR_pol-EMC(⁷Li) | E12-14-001 projected 4–5% on ratio, x ≤ 0.8 | x < 0.06 + Q² arm + 5%-level CBT-vs-QMC discrimination |
| tagged structure | Friščić e+³He ×2–3 improvement | α-tagged ⁶Li embedded-deuteron program |

## 6. Fact-check addendum (2026-08-21, primer revision "why Li / why polarization / fixed target vs collider")

Verified against primary sources while strengthening `reports/polarized_li_primer`;
corrections to entries above:

- **⁶Li cluster polarization (§2.3 table).** Schellingerhout–Kok–Coon–Adam
  (PRC 48:2714, erratum PRC 52:439) give a *valence-neutron polarization
  > 90%* (92–98% across models) with an α+d probability of 60–70%; the
  paper does not state "valence pair carries > 90% of the spin", and its
  "87%" refers to ³He. The ≈0.85–0.87 per-nucleon value is six-body QMC:
  Pudliner et al. PRC 56:1720 (n↑/n↓ = 1.93/1.07) and Wiringa et al.
  PRC 89:024305 Table I (1.924/1.076 → 0.848); SLAC E155 (PLB 463:339)
  adopted "87% of the Li polarization" for the effective deuteron.
- **⁷Li effective polarizations.** P_p = 0.866 / P_n = −0.037 are the values
  *quoted in PR12-14-001* from the VMC spin-projected nucleon numbers of
  Wiringa et al. (printed Table I: 1.934/1.066, 1.981/2.019 → 0.868/−0.038;
  online tables 0.8675/−0.0374). Shell model 13/15 = 0.867 (proposal),
  Cohen–Kurath via Cloët–Bentz–Thomas; earliest QMC 0.88/−0.04 (Pudliner).
  Say "VMC", not "GFMC".
- **Tagged tensor asymmetries O(1) at ~300 MeV/c** are from Cosyn–Weiss
  PRC 114:025202 (2026), arXiv:2603.23700 (Part I: 025201, 2603.23699) — not
  PRC 102:065204 (vector) or PLB 799:135035.
- **Resonance counts 81 (⁶Li) / 573 (⁷Li) / 63 (d) / 2457 (³He)** are from
  Hamwi–Devlin–Hoffstaetter, PRAB 29:073501 (2026), arXiv:2509.18558, Table I
  (G = −0.1818 / +1.5196 there; EPIOS Table III gives −0.178 / +1.532). Top
  energies ~138 / ~117 GeV/u are the rigidity limit 275 GeV × Z/A (137.5 /
  117.9), not a number quoted by EPIOS.
- **Source targets** (P_z ≥ 0.90, P_zz ≥ 0.80, ≥ 10¹³ atoms/s) are from the ANL
  ECRP proposal, not EPIOS (which gives no Li source numbers).
- **EIC baseline species.** CDR (doi:10.2172/1765663) and EIC Global
  Requirements (EIC-ORG-PLN-010 Rev. 03, 2026): polarized p and ³He in
  project scope; polarized d "with further development", explicitly outside
  scope. Lithium is not mentioned in the CDR or the Yellow Report; EPIOS calls
  ⁶,⁷Li "a desirable addition" and lists them as critical technology.
- **Yellow Report luminosity accounting:** 10 fb⁻¹ = 30 weeks at
  10³³ cm⁻² s⁻¹ with 60% efficiency (1.5 fb⁻¹/month) — "one year" is our
  rounding.
- **Fixed-target polarized Li** (new, §1.4 of the primer): COMPASS ⁶LiD
  D-polarization +54.2/−47.1% (Ball et al. NIM A 498:101), EST → ⁶Li 51 ± 3%,
  ⁷Li admixture 92 ± 4% (Koivuniemi SPIN2004), dilution f ≈ 0.35–0.45
  (NIM A 577:455; PLB 612:154); E155 ⁶LiD: P(⁶Li) = 97% of P(d), f ≈ 0.36 vs
  0.22 for ND₃ (PLB 463:339); ⁷LiH at Saclay: ⁷Li 47%, H 56% (Chaumette
  AIP Conf. Proc. 187:1275); PR12-14-001 expects ⁷Li 65–80% in ⁷LiD at 5 T,
  assumes P_z = 0.8 (P_zz = 0.55, P_zzz = 0.11 by EST). Tensor polarization
  in solids: EST P_zz = 2 − √(4 − 3P_z²) (0.20 at P_z = 0.5); RF enhancement
  to 29–37% in d-butanol (Keller–Crabb–Day NIM A 981:164504); E12-13-011
  C1-approved on demonstrating 30%. HERMES ABS target: P_zz = +0.89/−1.65
  pure-tensor states (effective 0.83). LOI12-16-006: ¹⁴N in NH₃, x < 0.3,
  2–12% alignment; "cos 2φ … difficult with a fixed target at JLab". No
  polarized-Li gas/jet target in lepton scattering exists (only Wittchow et
  al. PLB 59:29 (1975), low-energy reaction).
- **Spectator/recoil energetics:** α at 100/300 MeV/c → 1.3/12 MeV → CSDA
  range 7–12 / 170–240 μm in LiD (NIST ASTAR scaled); coherent ⁶Li recoil at
  |t| = 0.05 GeV² → 4.5 MeV. Target cells: JLab 2–5 cm, COMPASS 30–60 cm.
  BoNuS (p ≥ 70 MeV/c) and ALERT (α 230–400 MeV/c) use unpolarized gas.
- **Coherent recoil rigidity:** an intact coherent recoil keeps the beam
  rigidity for any A/Z; A/Z = 2 matters for the α/d *fragments* of ⁶Li.
  Chang et al. (PRD 113:032018) studied ⁷Li (17.75% IR-8 efficiency) but not
  ⁶Li.
- **Source history (primary numbers):** Heidelberg 1977 ⁶Li³⁺ ≤ 150 nA,
  |P| = 0.57–0.65 (Steffens NIM 143:409); 1987 atomic beam P_z, P_zz ≥ 0.85,
  ⁷Li accelerated at P_zz = 0.69 (Jänsch NIM A 254:7); FSU ⁶Li³⁺ > 150 nA,
  P_yy = −1.33 ± 0.02, P_y > 0.74, no loss in the linac (Myers NIM B 79:701;
  Mendez NIM A 329:37). Zelenski 2023: 2×10¹¹ ³He⁺⁺/pulse at ≥ 70% is the
  design goal; polarization transmission through the EBIS not yet measured
  (EPIOS: extraction measurement planned 2026).
- **Nuclear data:** Q(⁶Li) = −0.0806(6) fm², Q(⁷Li) = −4.00(3) fm² (Stone
  INDC(NDS)-0833, 2021; −4.06 is superseded), Q(d) = +0.2858 fm²;
  μ(⁶Li) = 0.8220, μ(⁷Li) = 3.2564 μ_N; S(α+d) = 1.474, S(α+t) = 2.468 MeV
  (AME2020); stable J = 1 nuclei are exactly ²H, ⁶Li, ¹⁴N; stable J = 3/2
  below A = 30: ⁷Li, ⁹Be, ¹¹B, ²¹Ne, ²³Na (IAEA ground-state table).
