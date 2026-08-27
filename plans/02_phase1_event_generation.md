# Phase 1 — Event-Generation-Level Study (fast simulation)

**Goal.** Establish the kinematic phase space, rates, backgrounds, and
statistical figures of merit for polarized e+⁶Li / e+⁷Li at the EIC, and
quantify spectator-tagging purity — the deliverable the ECRP proposal
explicitly assigns to "BeAGLE-class simulations". Everything here runs on a
laptop/workstation; no full detector simulation.

**Strategy in one paragraph.** No public generator produces *polarized*
e+A events — verified across BeAGLE, DJANGOH (HPOLAR ignored for A>1),
PEPSI/CLASDIS (nucleon-level only), eHIJING, Sartre, TOPEG, GCF (fetch-
verified survey, 2026-06-12). None needs to: the program separates into
(i) *rates and kinematics* from unpolarized tools (analytic fast-sim →
BeAGLE for anything involving nuclear breakup/fragments), and (ii)
*asymmetry injection* by reweighting with structure-function models — the
established EIC practice (DSSV route arXiv:2007.08300; JAM route
arXiv:2105.04434; event-level mechanics per the ECCE recipes NIM A 1056
(2023) 168563 and arXiv:2207.10890; no public reweighting package exists,
so our `polli_fastsim` grows into that role). For tensor/spin-1
observables the Cosyn–Weiss covariant spin-density-matrix formalism
(arXiv:2006.03033, includes tensor polarization) is the theory input.
Detector effects enter Phase 1 only as parameterized acceptance/smearing;
full simulation is Phase 2.

**Status legend:** ☐ todo ◐ started ☑ done

---

## Step 1.0 ☑ Bootstrap analytic fast-sim (rates + FOM skeleton)

Done in this session: `fastsim/polli_fastsim` (tested, 6/6 passing).
- x–Q² coverage, per-bin DIS rates, δ(g₁/F₁), δA_zz, δA_cos2φ maps for an
  energy scan of both isotopes; toy structure functions clearly labeled.
- First numbers (10 fb⁻¹/u, P_e = 0.7, P_z = 0.7, P_zz = 0.6): N_DIS ≈
  (2.5–5.7)×10⁹; δ(g₁/F₁) ≈ 2×10⁻² per bin at x ≈ 0.3 at √s_eN ≈ 19–20 GeV;
  δA_zz ≈ 9×10⁻⁴ per bin — to be compared with |A_zz(b₁)| ~ 10⁻³–10⁻².

## Step 1.1 ☐ Pin down scenario inputs (1–2 weeks, mostly reading)

1. Work through the EPIOS white paper (Atoian et al., arXiv:2510.10794,
   PRC 113:060501): verified there — G(⁶Li) = −0.178 / G(⁷Li) = +1.532,
   top energies ~138/~117 GeV/u, partial-snake schemes; **no Li luminosity
   number exists** → adopt and *state* a per-nucleon scaling assumption.
2. **Resolve the ⁶Li effective-polarization convention** (factor 2.4 in the
   g₁ FOM!): Cloët slides quote P_p = P_n = 1/3 (per-nucleon-normalized,
   2-of-6 dilution); cluster picture gives ≈ 0.87 (P_d) × 0.93 (D-state)
   ≈ 0.81 per cluster nucleon whole-nucleus. ⁷Li verified: P_p = +0.866,
   P_n = −0.037 (QMC via E12-14-001). Record adopted values + band in
   `beams.py`; cross-check against EPIOS Fig. 4 / I. Cloët directly (ANL).
3. Binning conventions: match HERMES b₁ x-points (zero crossing at x ≈ 0.2,
   b₁ ~ 0.1 at x ~ 0.01), E12-14-001 (0.06 < x < 0.8) and YR inclusive
   binning (~5 bins/decade) so every comparison is one-to-one.
4. Calendar anchor: aim Phase-1 money plots at the INT program "Towards
   Realizing the Program with Polarized Ion Beams at EIC",
   **March 22 – April 2, 2027**.

## Step 1.2 ◐ Replace toy structure functions (1–2 weeks)

*2026-06-12: grid backends wired and validated — `PartonF2` (CT18NLO) and
`PartonG1` (NNPDFpol11_100) via the `parton` package; toy F2 certified to
±40% (`scripts/validate_inputs.py`). Scenario curves added: CBT 2× / TMT 1×
polarized-EMC, HERMES-like vs convolution b1. Remaining below.*

1. Unpolarized: LHAPDF inside eic-shell (container ships LHAPDF6) or the
   pure-python `parton` package locally; CT18NNLO + EPPS21/nNNPDF3.0 nuclear
   ratios for ⁶,⁷Li (A=6,7 grids exist in EPPS21? if not, interpolate A or
   use light-nuclei convolution). Validate toy-vs-grid F₂ maps (expect ≤30%
   shifts, no FOM conclusions change).
2. Polarized: JAM (e.g., JAMpol via LHAPDF) or DSSV14 grids for g₁p/g₁n →
   A₁ maps; **two-camp** medium-modification curves digitized as
   `medium_ratio(x)`: CBT (PLB 642:210; ~2× unpolarized EMC) vs
   Tronchin–Matevosyan–Thomas (PLB 783:247; ≈ unpolarized) — the FOM
   question becomes "at what lumi/P do we discriminate the camps at 5σ".
   Gluon-spin EMC curves from Wang et al. (arXiv:2109.03591) for the
   dg₁/dlnQ² observable.
3. b₁ model: deuteron convolution (Cosyn–Dong–Kumano–Sargsian PRD 95:074036,
   |b₁| < 10⁻³ at x ≳ 0.2) vs Miller pion+hidden-color (PRC 89:045203,
   reproduces HERMES b₁ ~ 0.1 at x ~ 0.01) as the two scenario curves,
   rescaled by ⅓·P_d ≈ 0.29 for the ⁶Li embedded deuteron (our model — no
   ⁶Li b₁ theory exists; engage Cloët/Cosyn/Miller to publish one).
   Δ model: normalize shapes to ∫xΔdx = −0.012·α_s (Sather–Schmidt bag
   estimate, the LOI12-16-006 reference point) + flat Δ/F₁ ∈
   {10⁻³, 3×10⁻³, 10⁻²} scenarios; lattice φ-meson moment
   (Detmold–Shanahan PRD 94:014507) as the nonzero-existence argument.

**Deliverable:** updated FOM maps with credible central curves and scenario
bands; a short note fixing the input set.

## Step 1.3 ◐ Inclusive figure-of-merit study (2–3 weeks)

*2026-06-12: first versions of all three money plots exist
(`scripts/money_{polemc,b1,delta}.py`) with the error estimators validated
by toy-MC closure (`tests/test_closure.py`). First numbers: gluonometry
5σ at Δ/F₁=10⁻³ needs ~15–40 fb⁻¹/u; CBT-vs-TMT ≈5σ at x≈0.5–0.7 with
100 fb⁻¹/u; δA_zz ~ 10⁻³/x-bin at 10 fb⁻¹/u. To finalize: rerun on grid
inputs + adopted binning, add the items below.*

For each observable, produce the "money plots":
1. **Polarized EMC (⁷Li first):** projected δ(ΔR(x)) vs x for the energy
   scan; overlay CBT prediction and JLab E12-14-001 projected errors —
   demonstrate the order-of-magnitude x–Q² extension and Q²-lever arm
   (gluon-spin EMC via dg₁/dlnQ²). Quote vs P_z ∈ {0.5, 0.7, 0.9} and
   lumi ∈ {1, 10, 100} fb⁻¹/u.
2. **b₁(⁶Li) vs b₁(d):** δA_zz(x) per bin vs predictions; significance of a
   d-vs-⁶Li *difference* (the medium-modification signal) including the
   P_zz dependence — this sets the tensor-polarization requirement, feeding
   back to the source spec (P_zz ≥ 0.80).
3. **Gluonometry Δ:** 95% CL exclusion / 5σ discovery reach in Δ/F₁ vs
   integrated luminosity from the cos 2φ fit; required transverse-spin
   running time. Note: needs *transverse* ion polarization at IP6 → flag to
   04_open_questions (spin-rotator configuration for ions).
4. Cross-check δA estimators against a toy MC (generate Poisson counts in φ
   bins and spin states, fit, compare pulls) — validates the analytic FOM.

## Step 1.4 ☐ Radiative corrections sanity pass (1 week)

g₁/A∥ at low x and the φ-modulations are RC-sensitive. Run **DJANGOH
4.6.22** (github.com/spiesber/DJANGOH, maintained, HERACLES full O(α) EW;
EIC reference arXiv:1309.5327) on an *effective polarized nucleon* —
hadron polarization is ignored for A>1 and there is no Fermi motion, so
follow the arXiv:2406.05591 ³He workaround (effective nucleon + spectators
added by hand from spectral functions; CLASDIS is the alternative polarized
nucleon-level generator, used by Friščić et al. for e+³He). Outcome: a
multiplicative RC-uncertainty band on the FOMs, and a decision whether
Phase-2 needs full RC treatment (note: RC on tensor A_zz/cos 2φ is
uncharted — flag to theory colleagues).

## Step 1.5 ◐ BeAGLE e+Li breakup & tagging study (4–6 weeks, core novelty)

*2026-06-12: the cluster-IA seed (item 3 cross-check + item 4 routing) is
implemented (`polli_fastsim/spectator.py`, `farforward.py`,
`scripts/tagging_acceptance.py`): ⁷Li α-tag 96–99% into RP; ⁶Li α-tag
3–9% at IP6 (tail-dominated, the quantitative beam-blindness statement;
1.7% since the near-beam cut became angular, 2026-08-25 — plans/08 C1);
⁷Li t-tag ~0. BeAGLE itself (items 1–2, evaporation backgrounds, purity)
still todo — access is the long pole.*

The proposal's explicit ask: charged-α and neutral-fragment tagging purity.
BeAGLE (arXiv:2204.11998) is the only breakup-capable eA generator; it runs
arbitrary (A,Z) but is *untuned* for Li — known caveats (verified): A>4
inherits the **C-12 Fermi-momentum parameterization** (A=2,3,4 use
Ciofi-degli-Atti–Simula); geometry is Woods–Saxon with **no α+d / α+t
cluster structure** (fragments come from FLUKA statistical de-excitation,
not cluster knockout); no Li shadowing map ships (use `genShd=1`); code
frozen since v1.03.02 (2023). The deuteron needed special treatment
(`deutfix.f`, IA mode of arXiv:2005.14706 / 2108.08314) — expect Li to
need analogous care.

1. **Run environment.** ◐ *2026-06-12: local build is one step from done —
   LHAPDF 5.9.1 (+ CTEQ6L1 and EPS09 A=4,6,12 grids; `anear.f` maps both
   isotopes to the valid A=6 EPS09 nucleus), RAPGAP-3.302 libs, PYTHIA6,
   and CERNLIB 2024 core are built in `~/Projects/eic/beagle_deps/`;
   the gfortran-11 recipe is pre-validated (PYTHIA/RADGEN/PyQM libs
   compile; DPMJET stops only on FLUKA's `(DIMPAR)` includes). The single
   missing piece is **FLUKA (personal license — user registers at
   fluka.org, gfortran-11 64-bit build, export FLUPRO)**; then
   `tools/beagle/build_beagle.sh` finishes the build. See
   `tools/beagle/README.md`. Alternatives remain: BNL/JLab prebuilds, and
   the official e+d/e+³He EVGEN samples are downloadable from here via
   xrootd (verified) for the control study.
2. **Configuration.** Cards: `TARPAR 6. 3.` / `TARPAR 7. 3.` (mixed n/p
   mode), `MOMENTUM <Ee> <p/u>` (⁶Li ≤137.5, ⁷Li ≤117.9 GeV/u), `genShd=1`,
   radgen on, `L-TAG` cuts as in the local e+D example
   (`../BeAGLE/Examples/eD_18x135_*.inp`); PYTHIA card from `eAt1dfJn`.
   ≥10⁶ DIS events per isotope per energy.
3. **Validation for A=6,7 — mandatory, not optional.**
   - e+d control: reproduce the published BeAGLE e+d tagging distributions
     (Tu et al. 2005.14706) with our setup before touching Li.
     ◐ *2026-06-12 first pass on the official BeAGLE eH2 'en' 10×130
     sample (50k events streamed via xrootd; `tools/analysis/`):
     spectator-proton routing confirmed — 96.6% OMD (cluster model:
     98.1%), x_L and pT bulk agree with the Hulthén model, but BeAGLE's
     pT tail is 2.2× (pT>0.2) to 13× (pT>0.45) harder ⇒ near-beam (R≈1)
     tag acceptances are tail-model-dependent at the O(10) level — the
     ⁶Li α-tag number needs VMC/BeAGLE input, as suspected.*
     ☑ *2026-08-26 second pass, on BeAGLE **1.03.02-3.1** eH2 'en' 9×130
     (20k events, same xrootd route; `--beta-scan` added to
     `ed_control_analysis.py`): the first pass reproduces on the newer
     sample — routing 93.1% OMD, tail 2.24× / 4.86× / 12.95× at
     p_T > 0.2 / 0.3 / 0.45 GeV — and restricting to the spectator peak
     x_L ∈ [0.9, 1.1), which removes any target-fragmentation
     contamination, makes it worse, not better: 2.50× / 6.92× / 28×.
     The scan says **no β reproduces the shape**: BeAGLE has a narrower
     core and a harder tail than the two-parameter Hulthén can have at
     once, the best shape fit sits at β = 0.40 GeV (the top of the
     documented 0.20–0.40 range) and still falls 3× short at 0.3 GeV and
     10× at 0.45.  Since the ⁶Li α tag is entirely a p_T-tail measurement
     — 0.40 GeV of p_T for a 550 GeV α under the high-acceptance envelope
     — the published tag acceptance is model-limited **from below** and
     its model uncertainty is one-sided upward.  Still owed: the Tu et al.
     comparison the step actually asks for, and BeAGLE's A = 6 treatment
     is not its deuteron one.*
   - Fragment yields (α, t, d, ³He, p, n, residues with E*) vs
     photo-/electro-disintegration data and cluster-model expectations
     (⁷Li → α+t dominance at low E*); event-by-event momentum conservation.
   - **Cluster-IA cross-check:** build a light-front impulse-approximation
     toy (α+d for ⁶Li, α+t for ⁷Li, cluster momentum densities from
     few-body theory) in the spirit of the deuteron special treatment;
     compare spectator x_L–pT spectra against BeAGLE's evaporation picture.
     The difference spans the model uncertainty on tagging acceptance.
4. **Physics outputs.**
   - Spectator spectra: x_L, pT, θ per fragment species — the far-forward
     acceptance inputs.
   - **Rigidity routing** (corrected, R = (A_f·Z_beam)/(A_beam·Z_f); see
     plans/03 §2.2 table): ⁷Li beam → α at R = 0.86 lands **in the Roman
     Pots** (x_L window 0.6–0.95) — the IP6-friendly tag; p → OMD; n → ZDC;
     t → **no IP6 coverage** (R = 1.29). ⁶Li beam → α/d at R = 1.0 are
     **beam-blind** below the RP pT cutoff (0.2–0.45 GeV/c by optics):
     fold the soft cluster pT(α) spectrum with the 10σ cutoff — this single
     number decides whether ⁶Li d-cluster tagging works at IP6 or needs
     the IR-8 secondary focus; ³He → RP (R = 0.75); p → OMD.
   - **Tagging purity/efficiency:** P(tag α | DIS on d/t-cluster) vs
     α from evaporation/INC background; same for n tags (+ de-excitation γ
     in ZDC). Defines the tagged-sample dilution for the tagged FOM.
   - Backgrounds: low-Q² photoproduction leakage; multi-fragment
     combinatorics; coherent/diffractive e+⁶Li (intact-nucleus RP
     signature — also a future signal channel for coherent studies).
5. **Output format.** BeAGLE text → eic-smear `BuildTree`/`TreeToHepMC` →
   HepMC3 once, so identical samples feed Phase 2 (then `abconv` adds beam
   effects there).

## Step 1.6 ☐ Parameterized detector smearing (2 weeks)

Tool status (verified 2026-06): **eic-smear is alive** (v1.1.17, 05/2026;
reads BeAGLE natively) with YR-era cards including **"Matrix 0.1 +
Far-Forward"** — use that; delphes_EIC is dormant (no ePIC card) and no
official ePIC fast-sim exists. Apply to the fast-sim/BeAGLE samples:
scattered-electron resolution → x–Q² migration matrices → effect on per-bin
FOM and on the cos 2φ amplitude (φ resolution). Verify YR-style binning
keeps purity ≳ 0.8 per bin; else rebin. Real far-forward acceptances come
only from Phase-2 full sim — keep the FF parameterization swappable.

## Step 1.7 ☐ Synthesis & write-up (2 weeks)

- Note (10–15 pages): phase space, rates, FOM money plots, tagging purity
  tables, RC band, energy-scan recommendation, source-polarization
  requirement flow-down (what P_z/P_zz/lumi the physics actually needs).
- Feed figures back into: ECRP renewal material, EPIOS white-paper
  follow-ups, and an EIC user-group / DNP talk.

---

## Suggested order & effort

1.0 ☑ → 1.1 → 1.2 → 1.3 (first money plots, ~6 weeks in) → 1.5 in parallel
with 1.4 → 1.6 → 1.7. Total ≈ 3–4 months of focused effort; BeAGLE access
(1.5.1) is the long-pole external dependency — start it immediately.

## Risks specific to Phase 1

| risk | mitigation |
|---|---|
| BeAGLE invalid for A=6,7 breakup (C-12 n(k), no cluster geometry, frozen code) | validation step 1.5.3 incl. cluster-IA cross-check; fallback: cluster-model toy fragmenter (α+d / α+t momentum densities + flat E*) good enough for acceptance maps |
| BeAGLE access (FLUKA license, "no mere mortal" build) | prebuilt BNL/JLab/CVMFS installs; start access requests immediately |
| ⁶Li α-tag beam-blind at IP6 (R = 1.0 vs RP pT cutoff) | quantify pT-tail acceptance early (step 1.5.4); lead the tagging story with ⁷Li (α → RP works); document IR-8 secondary-focus case |
| No nuclear PDF grids at A=6,7 | interpolate EPPS21 in A; or convolution from d/³He/⁴He |
| Transverse ion polarization at IP unavailable | gluonometry FOM quoted conditional on rotator configuration; raise early with EPIOS/C-AD (04_open_questions) |
| Tensor (λ=0) bunches operationally undefined | source RF transitions can prepare m=0; needs machine fill-pattern concept — document requirement, don't solve |
| ⁶Li effective-polarization convention (1/3 vs 0.81, factor 2.4 in FOM) | resolve in step 1.1 with Cloët before any public plot |
