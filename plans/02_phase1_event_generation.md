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

## Step 1.1 ◐ Pin down scenario inputs (1–2 weeks, mostly reading)

*2026-08-28: items 1 and 4 are done (notes below); item 2's record is in place
but the value itself is external (plans/04 #6); item 3 is untouched — the
analysis still runs the generic 40×30 log grid, not the YR ~5 bins/decade. The
header glyph was ☐ and is set to ◐ here on the audit's own reading — the step
itself stays open, and no other open item's marker was moved.*

1. Work through the EPIOS white paper (Atoian et al., arXiv:2510.10794,
   PRC 113:060501): verified there — G(⁶Li) = −0.178 / G(⁷Li) = +1.532,
   top energies ~138/~117 GeV/u, partial-snake schemes; **no Li luminosity
   number exists** → adopt and *state* a per-nucleon scaling assumption.
   ☑ *2026-08-28: done — EPIOS pp. 12–13 are encoded in `beams.py`
   (`EPIOS_GAMMA_BYPASS`, `EPIOS_GAMMA_SHIFT_RANGE`, `epios_window_of`, guarded
   by `tools/consistency_check.py`) with the G factors and snake schemes in
   plans/04 A1; the per-nucleon assumption is adopted as
   `fom.Scenario.lumi_fb_per_nucleon = 10.0` and stated with its "no Li number
   exists" premise in Report 1 §3.1 and Report 3 §4 (Reports 0 and 2 carry
   the placeholder without the premise).*
2. **Resolve the ⁶Li effective-polarization convention** (factor 2.4 in the
   g₁ FOM!): Cloët slides quote P_p = P_n = 1/3 (per-nucleon-normalized,
   2-of-6 dilution); cluster picture gives ≈ 0.87 (P_d) × 0.93 (D-state)
   ≈ 0.81 per cluster nucleon whole-nucleus. ⁷Li verified: P_p = +0.866,
   P_n = −0.037 (QMC via E12-14-001). Record adopted values + band in
   `beams.py`; cross-check against EPIOS Fig. 4 / I. Cloët directly (ANL).
   ☐ *2026-08-28: the record exists (`beams.py` adopts 1/3 with the band and
   provenance) but it was written with this bullet, so what remains is the
   resolution itself — external, plans/04 #6. And the recorded value is not the
   operative one: `ToyG1.g1_nucleus` reads the per-nucleon 1/3 as whole-nucleus
   while callers divide by A, applying the 2-of-6 dilution twice (plans/08 D7);
   the fix changes a documented convention, so it is the author's call.*
3. Binning conventions: match HERMES b₁ x-points (zero crossing at x ≈ 0.2,
   b₁ ~ 0.1 at x ~ 0.01), E12-14-001 (0.06 < x < 0.8) and YR inclusive
   binning (~5 bins/decade) so every comparison is one-to-one.
4. Calendar anchor: aim Phase-1 money plots at the INT program "Towards
   Realizing the Program with Polarized Ion Beams at EIC",
   **March 22 – April 2, 2027**.
   ☑ *2026-08-28: adopted and propagated — plans/00, plans/04 #2 and the
   plans/07 milestone table all work back from this date.*

## Step 1.2 ◐ Replace toy structure functions (1–2 weeks)

*2026-06-12: grid backends wired and validated — `PartonF2` (CT18NLO) and
`PartonG1` (NNPDFpol11_100) via the `parton` package; toy F2 certified to
±40% (`scripts/validate_inputs.py`). Scenario curves added: CBT 2× / TMT 1×
polarized-EMC, HERMES-like vs convolution b1. Remaining below.*

*2026-08-28: items 2 and 3 below are done for the polarized-EMC and b₁
curves — all four are digitized from the published figures into
`fastsim/polli_fastsim/data/` by `tools/digitize_figure.py`
(`data/SOURCES.md`), the two money plots draw them by default and keep the
old shapes behind `--emc-mode constant` / `--transfer legacy`.  Item 1
(EPPS21 for the unpolarized ratio) closed on 2026-08-29 — see step 1.2.1
below — and the Wang gluon-spin curves are the remainder.*

1. Unpolarized: LHAPDF inside eic-shell (container ships LHAPDF6) or the
   pure-python `parton` package locally; CT18NNLO + EPPS21/nNNPDF3.0 nuclear
   ratios for ⁶,⁷Li (A=6,7 grids exist in EPPS21? if not, interpolate A or
   use light-nuclei convolution). Validate toy-vs-grid F₂ maps (expect ≤30%
   shifts, no FOM conclusions change).
   ☑ *2026-08-29: closed for ⁶Li.  `polarized.UNPOL_EMC_MODE = "epps21"` is
   EPPS21nlo_CT18Anlo_Li6's F₂ per nucleon over CT18ANLO's free isoscalar
   nucleon — EPPS21's own proton baseline, so the fit cancels (CT18NLO
   until 2026-08-29, 4.2% shallower) — with `nnnpdf` (nNNPDF30 A = 6) as the spread and the retired
   hand-written 12-point shape reachable as `mode="table"`.  Mean valence
   depletion over 0.35 < x < 0.65: 0.03105 (EPPS21), 0.01372 (nNNPDF3.0),
   0.06459 (table), against 0.05835 for CBT's own model curve; EPPS21's 90%
   CL Hessian band on its own value is +0.039 / −0.041.  It is the common
   baseline both polarized-EMC camps are now transferred onto (item 2), and
   the spread between the choices is the leading uncertainty on that figure
   of merit.  ⁷Li has no nuclear PDF grid from anyone, so the ⁷Li
   projection uses the A = 6 baseline; that substitution is the residue of
   this item.*
2. Polarized: JAM (e.g., JAMpol via LHAPDF) or DSSV14 grids for g₁p/g₁n →
   A₁ maps; **two-camp** medium-modification curves digitized as
   `medium_ratio(x)`: CBT (PLB 642:210; ~2× unpolarized EMC) vs
   Tronchin–Matevosyan–Thomas (PLB 783:247; ≈ unpolarized) — the FOM
   question becomes "at what lumi/P do we discriminate the camps at 5σ".
   Gluon-spin EMC curves from Wang et al. (arXiv:2109.03591) for the
   dg₁/dlnQ² observable.
   ☑ *2026-08-28 for the two camps: CBT Fig. 6 (⁷Li, Q² = 5 GeV², both the
   Eq.-23 and Eq.-26 curves and the unpolarized one) and TMT Fig. 4
   (nuclear matter, Q² = 10 GeV²) are in `cbt_polemc_7Li_Q5.csv` and
   `tmt_polemc_nm_Q10.csv`.  Because the two are different targets at
   different scales, they are compared as the ratio of EMC effects applied
   to one common baseline with a single valence strength factor.  That
   factor is fitted in 0.35 < x < 0.65 and applied at every x, so the
   comparison is a statement about the two papers only inside that window:
   the published polarized curves agree to better than 0.008 over
   0.028 < x < 0.30, while inside the window the transferred TMT depletion
   tracks the baseline's own unpolarized one against CBT's roughly twice it.
   *2026-08-29: the baseline is no longer CBT's own model curve but EPPS21's
   ⁶Li over CT18ANLO (step 1.2.1), whose valence depletion 0.0311 is just
   over half of CBT's 0.0584, so the factors are 0.5322 and 0.2113 in place
   of 1 and 0.397 and the transferred TMT depletion tracks it to within
   0.003 while
   ΔR separates by 0.021 at x = 0.36 falling to 0.006 at 0.65.  The answer
   to the FOM question is 1.43σ in the best bin of that window at
   100 fb⁻¹/u, not 5σ and not the 2.7σ the model baseline gave; the 1.95σ
   the unrestricted scan returns at x = 0.141 is carried by the transfer,
   not by the calculations, and `money_polemc.py` prints both, together with
   the baseline spread that is now the leading uncertainty on the answer.*
   Wang's Figure 3 is NOT digitized:
   the dg₁/dlnQ² observable has no money plot (`refs/2109.03591.pdf` and
   `data/SOURCES.md` record where the curves are).*
3. b₁ model: deuteron convolution (Cosyn–Dong–Kumano–Sargsian PRD 95:074036,
   |b₁| < 10⁻³ at x ≳ 0.2) vs Miller pion+hidden-color (PRC 89:045203,
   reproduces HERMES b₁ ~ 0.1 at x ~ 0.01) as the two scenario curves,
   rescaled by ⅓·P_d ≈ 0.29 for the ⁶Li embedded deuteron (our model — no
   ⁶Li b₁ theory exists; engage Cloët/Cosyn/Miller to publish one).
   ☑ *2026-08-28: both digitized (`b1_miller.csv` from PRC 89:045203 Fig. 5,
   `b1_cdks_q2p5.csv` from PRD 95:074036 Fig. 4 theory-1 SD+DD, plus the two
   Q² sets), and they are now ABSOLUTE b₁ rather than a shape times F₁, so
   `toy_b1`/`b1_convolution` ignore the F₁ they are handed (`mode='toy'`
   restores the shapes).  The ⅓·P_d rescaling became ⅓ × 0.921947 = 0.3073:
   the rank-2 tensor transfer, not the vector 0.87 (plans/08 D9).  The two
   camps are no longer a factor ten apart — Miller reaches 11σ per bin at
   x ≈ 0.07 while CDKS stays under 0.2σ everywhere.*
   Δ model: normalize shapes to ∫xΔdx = −0.012·α_s (Sather–Schmidt bag
   estimate, the LOI12-16-006 reference point) + flat Δ/F₁ ∈
   {10⁻³, 3×10⁻³, 10⁻²} scenarios; lattice φ-meson moment
   (Detmold–Shanahan PRD 94:014507) as the nonzero-existence argument.
   ☑ *2026-08-28: done in `fastsim/polli_fastsim/delta_models.py` — `C_BAG =
   −0.012` with `solve_A_interp_a/b` enforcing ∫xΔdx = c·α_s (unit-tested in
   `test_delta_models.py`), `make_toy(scale)` for the flat scenarios, and the
   lattice argument carried in plans/06 §6.4b and Report 2.*

**Deliverable:** updated FOM maps with credible central curves and scenario
bands; a short note fixing the input set.

## Step 1.3 ◐ Inclusive figure-of-merit study (2–3 weeks)

*2026-06-12: first versions of all three money plots exist
(`scripts/money_{polemc,b1,delta}.py`) with the error estimators validated
by toy-MC closure (`tests/test_closure.py`). First numbers: gluonometry
5σ at Δ/F₁=10⁻³ needs ~15–40 fb⁻¹/u; CBT-vs-TMT ≈5σ at x≈0.5–0.7 with
100 fb⁻¹/u; δA_zz ~ 10⁻³/x-bin at 10 fb⁻¹/u. To finalize: rerun on grid
inputs + adopted binning, add the items below.*

*2026-08-28, revised 2026-08-29: the CBT-vs-TMT line above is retired.  It
was the two camps as the constants 2 and 1 on a hand-written EMC table,
which made them separate as |1 − R_EMC(x)| and therefore grow with x.  On
the digitized curves, transferred onto the EPPS21 unpolarized baseline that
became the default on 2026-08-29, the separation is 0.023 / 0.021 / 0.018 /
0.032 at x = 0.09 / 0.28 / 0.45 / 0.71 against δΔR = 0.0423 / 0.0403 /
0.0595 / 0.1868 on grid inputs at 10 fb⁻¹/u — 1.75 / 1.66 / 0.97 / 0.53 σ
at 100 fb⁻¹/u, and nothing above 1σ per bin at 10.  The replacement headline
is the valence window in which the transfer is defined and the two published
curves genuinely differ: best bin x = 0.355, 0.45σ at 10 fb⁻¹/u and 1.43σ at
100.  The x ≈ 0.5–0.7 of the old line is inside that window; what is retired
is the 5σ, not the location.  The model baseline the day before gave twice
those separations, which `--emc-baseline cbt` reproduces.*

For each observable, produce the "money plots":
1. **Polarized EMC (⁷Li first):** projected δ(ΔR(x)) vs x for the energy
   scan; overlay CBT prediction and JLab E12-14-001 projected errors —
   *2026-08-28: the CBT overlay is now the digitized Eq.-23 curve of their
   Fig. 6 rather than a 2× stand-in (`money_polemc.py`); both panels shade
   the valence window 0.35 < x < 0.65 in which the TMT → ⁷Li transfer is
   defined, and the upper panel also draws TMT's published nuclear-matter
   curve untransferred as a dotted line, because the significance panel is
   otherwise read as claiming a low-x discrimination the two published
   curves do not have — untransferred, TMT lies on top of CBT below
   x ≈ 0.3; the E12-14-001 projected errors are still not overlaid.*
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
   ☑ *2026-08-28: done — `money_delta.py` gives the 5σ reach over the
   (Δ/F₁, luminosity) plane for all three configurations, the running time is
   stated in years for both channels — inclusively as L_5σ ÷ 10 fb⁻¹/u per year
   (`fastsim/notes/money_delta_note_2026-07-24.md`: worst case L_5σ ≈ 0.29 fb⁻¹/u,
   ≈ 80 hours, all nine cases inside a one-year program), and for the coherent tag
   by `tagging_optics.py`'s 2.8/4.4/2.6 yr (Report 1) — and the flag is
   plans/04 #2. The 95% CL band was never drawn: it is the same curve shifted
   by (1.645/5)².*
4. Cross-check δA estimators against a toy MC (generate Poisson counts in φ
   bins and spin states, fit, compare pulls) — validates the analytic FOM.
   ☑ *2026-08-28: done in `fastsim/tests/test_closure.py` and, through the full
   generator chain, `evgen/tests/test_pseudoexp.py` with
   `evgen/scripts/closure_fom.py` — pseudo-experiment spreads close on every
   analytic FOM map (plans/00 run 4).*

## Step 1.4 ☑ Radiative corrections sanity pass (closed 2026-08-28)

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

— *superseded and then closed (2026-08-28). plans/08 §8.3 D3 voids two of the
four deliverables — collinear ISR generates exactly zero fake cos φ′/cos 2φ′ and
the ratio cancellation is already demonstrated — so an external DJANGOH build on
an effective nucleon no longer buys what it was chosen for. The two that survive
are done in this repository instead:* `evgen/polligen/radiative.py` samples the
exponentiated leading-log photon spectrum D(z, Q²) = (t/2)z^(t/2−1)S(t) −
(t/4)(2−z), t = (2α/π)[ln(Q²/mₑ²) − 1] (Kuraev–Fadin; Nicrosini–Trentadue;
∫D = 1 + O(t²), residual 7×10⁻⁴ at t = 0.070) and applies it as the kinematic
migration this step asked for, through a **default-off** hook
`recopseudo.RecoResponse(isr=…)` and the driver `money_cos2phi_reco.py --isr`.
Outcome (mid configuration, four sweet spots, common random numbers, 1600
pseudo-events per cell, mean ± sem over eight response seeds — one draw of the
response scatters by 4–14% of the bound, so it is averaged with `--isr-seeds`): the multiplicative RC band on Δ̂ is **+0.62 ± 0.03 / +0.50 ± 0.02
/ +0.94 ± 0.03 / +1.22 ± 0.02%** with the published generator window, rising to
a 1.8–2.8% band as the window is opened to Q² ≥ 0.15–0.02 GeV², where the
truncated low-Q² feed-in saturates, i.e. **≤ 2.9%** — inside the ≤5% gate of
plans/07 WP4, so **Phase 2 does not need a full RC treatment for the unpolarized
sector**. A HERA-style E − p_z window would bring it to ≤ 0.25% while keeping
87% (25% Gaussian y stand-in) / 98% (PYTHIA, calibrated Σ) of the
non-radiative rate; it is a documented contingency and is not applied.
What remains uncharted and is flagged to
theory colleagues is unchanged: RC on the tensor observables A_zz / cos 2φ′
(plans/05 §5.5) — no unpolarized study stands in for it. Written up in Report 2
§7 and its Table 2 row; 30 tests in `evgen/tests/test_radiative.py` and
`test_reco.py`.

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
     ☐ *2026-08-28: the toy is built (`spectator.py`, `tagged.py`); the
     comparison this bullet asks for is not, and is blocked with the BeAGLE
     build — the e+d control runs against BeAGLE's special deuteron treatment,
     not against its A > 4 evaporation machinery.*
4. **Physics outputs.**
   - Spectator spectra: x_L, pT, θ per fragment species — the far-forward
     acceptance inputs.
     ☐ *2026-08-28: only the two-body cluster channels (α, d, t) exist
     (`spectator.py`, `fastsim/out/spectator_*.png`); ³He, p, n and excited
     residues come from BeAGLE's FLUKA de-excitation and are blocked with the
     A = 6,7 run itself.*
   - **Rigidity routing** (corrected, R = (A_f·Z_beam)/(A_beam·Z_f); see
     plans/03 §2.2 table): ⁷Li beam → α at R = 0.86 lands **in the Roman
     Pots** (x_L window 0.6–0.95) — the IP6-friendly tag; p → OMD; n → ZDC;
     t → **Roman-Pot silicon on the inner side of the bend** (R = 1.29;
     "no IP6 coverage" was a routing assumption until 2026-08-28, when
     `farforward.over_rigid_route` gave the fast simulation the R > 1
     branch the `tools/fullsim` gun scan had asked for — still without a
     beam envelope or a reconstruction, plans/03 §2.2, plans/09 B1). ⁶Li beam → α/d at R ≈ 1 (0.998 / 1.005) are
     **beam-blind** below the RP pT cutoff (0.2–0.45 GeV/c by optics):
     fold the soft cluster pT(α) spectrum with the 10σ cutoff — this single
     number decides whether ⁶Li d-cluster tagging works at IP6 or needs
     the IR-8 secondary focus; ³He → RP (R = 0.75); p → OMD.
     ☑ *2026-08-28: the number is in — the ⁶Li α tag is 1.9 / 1.7 / 2.6% at the
     Yellow Report optics of the three configurations and 31 / 22 / 29% at a
     lithium tagging optics costing 1/6.8–1/12.8 of the luminosity (the tagging
     triple read 35 / 28 / 29% until the per-configuration transport became the
     baseline on 2026-08-29)
     (`fastsim/out/tagging_acceptance.txt`, Report 3 Table 6); the near-beam cut
     is now the angular 10(σ_h, σ_v) envelope, not the 0.2–0.45 GeV/c p_T cutoff
     this bullet still words it as. The decision it forces — a tagging optics or
     the IR-8 secondary focus (≈20%), not any published IP6 optics — is the
     programme's one machine request (Report 1 §6.1).*
     — *superseded (2026-08-28): the bullet's formula and its ⁷Li p → OMD cell.
     Rigidity is `spectator.py`'s mass-based R = (m_spec/Z_spec)/(m_beam/Z_beam),
     which separates the ⁶Li α (0.99813, under the orbit) from the d (1.00452,
     over it) where the A·Z arithmetic puts both at exactly 1; and the ⁷Li proton
     at R = 0.43081 falls below the corrected OMD window R ∈ [0.45, 0.65], i.e.
     it is lost, not routed (`route_charged`, `fastsim/tests/test_spectator.py`).
     The plans/03 §2.2 table this bullet cites carries the same note and is no
     longer the authority.*
   - **Tagging purity/efficiency:** P(tag α | DIS on d/t-cluster) vs
     α from evaporation/INC background; same for n tags (+ de-excitation γ
     in ZDC). Defines the tagged-sample dilution for the tagged FOM.
   - Backgrounds: low-Q² photoproduction leakage; multi-fragment
     combinatorics; coherent/diffractive e+⁶Li (intact-nucleus RP
     signature — also a future signal channel for coherent studies).
5. **Output format.** BeAGLE text → eic-smear `BuildTree`/`TreeToHepMC` →
   HepMC3 once, so identical samples feed Phase 2 (then `abconv` adds beam
   effects there).
   — *superseded (2026-08-28): the chain went another way — the official BeAGLE
   samples are read as HepMC3 through `pyHepMC3.rootIO`, the hadronic final
   state arrives as PYTHIA 8 `.npz` (`polligen.hfs`), and Phase-2 input is
   written by `tools/fullsim/ion_gun_hepmc.py`; eic-smear appears nowhere in
   the tree.*

## Step 1.6 ◐ Parameterized detector smearing (2 weeks)

Tool status (verified 2026-06): **eic-smear is alive** (v1.1.17, 05/2026;
reads BeAGLE natively) with YR-era cards including **"Matrix 0.1 +
Far-Forward"** — use that; delphes_EIC is dormant (no ePIC card) and no
official ePIC fast-sim exists. Apply to the fast-sim/BeAGLE samples:
scattered-electron resolution → x–Q² migration matrices → effect on per-bin
FOM and on the cos 2φ amplitude (φ resolution). Verify YR-style binning
keeps purity ≳ 0.8 per bin; else rebin. Real far-forward acceptances come
only from Phase-2 full sim — keep the FF parameterization swappable.

☑ *2026-08-28: the smearing itself is delivered, by a home-grown chain rather
than eic-smear — `polligen/reco.py` (η-dependent EMCal, tracking σ_p/p and σ_θ,
ε_eID, `smear_electron`) on the electron side and `polligen/hfs.py` on an 8 M-event
PYTHIA 8 sample on the hadron side, with the (x, Q²) migration measured by
`recopseudo.RecoResponse.bin_summary` and its cost quoted per bin (reco
δÂ = (0.9–3.0)×10⁻⁴ against truth-level (1.4–4.5)×10⁻⁴; Report 1 §5.2, Report 2
§7). The FF parameterization stayed swappable: `farforward.Optics` is an argument
to `route_charged`/`acceptance_summary` with four interchangeable instances, and
the measured ePIC pot aperture already substitutes for it (`reco.rp_measure`,
`--rp-aperture measured`).*

— *superseded (2026-08-28): the eic-smear "Matrix 0.1 + Far-Forward" tool choice
— the response was built instead on Yellow Report requirements plus sourced ePIC
design numbers (refs/README.md), and the far-forward went past the fast-sim tier
to a measurement in the real geometry (`tools/fullsim/ion_gun_hepmc.py`).*

☐ *2026-08-28: still open — "purity ≳ 0.8 per bin, else rebin" is not met
(measured 0.56–0.75 calibrated, 0.42–0.68 uncalibrated) and nothing was rebinned;
the project changed the unfolding instead (`recopseudo.fold_shape_fit`,
plans/08 A6). The gate needs retiring or re-deciding — plans/03 §2.3 still
carries it as a live Phase-2 requirement.*

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
| BeAGLE invalid for A=6,7 breakup (C-12 n(k), no cluster geometry, frozen code) | validation step 1.5.3 incl. cluster-IA cross-check; fallback: cluster-model toy fragmenter (α+d / α+t momentum densities + flat E*) good enough for acceptance maps — *2026-08-28: the fallback carries every published tagging number, but of step 1.5.3 only the e+d control ran; still owed and unblocked are the Tu et al. 2005.14706 comparison and the flat-E\* fragmenter (no E\* code exists)* |
| BeAGLE access (FLUKA license, "no mere mortal" build) | prebuilt BNL/JLab/CVMFS installs; start access requests immediately |
| ⁶Li α-tag beam-blind at IP6 (R = 1.0 vs RP pT cutoff) | quantify pT-tail acceptance early (step 1.5.4); lead the tagging story with ⁷Li (α → RP works); document IR-8 secondary-focus case — *2026-08-28: quantified per configuration and per optics (1.7–2.6% at the YR optics, 22–31% at the tagging optics since the per-configuration transport of 2026-08-29, 28–35% before it) and IR-8 priced; Report 0 §5.4/Table 3 restated on it* |
| No nuclear PDF grids at A=6,7 | interpolate EPPS21 in A; or convolution from d/³He/⁴He — *superseded (2026-08-28): A = 6 grids exist and are in use (EPPS21nlo_CT18Anlo_Li6, nNNPDF30 A6, compared against each other in the money-Δ line); the risk survives for ⁷Li alone, where LHAPDF has nothing* |
| Transverse ion polarization at IP unavailable | gluonometry FOM quoted conditional on rotator configuration; raise early with EPIOS/C-AD (04_open_questions) — *2026-08-28: the conditional quoting is done (`money_delta.py` docstring, Report 1 §3.1 and Table 4 #8); the raise itself has not happened — plans/04 #2 records owner and default, no contact yet* |
| Tensor (λ=0) bunches operationally undefined | source RF transitions can prepare m=0; needs machine fill-pattern concept — document requirement, don't solve — ☑ *2026-08-28: documented past the ask — plans/04 #3 carries the requirement plus the measured consequence (a 10⁻³ inter-fill difference of the cos 2φ′ harmonic fakes 5.6×10⁻⁴), which turns bunch-by-bunch alternation into a requirement of the measurement (Reports 1 and 2)* |
| ⁶Li effective-polarization convention (1/3 vs 0.81, factor 2.4 in FOM) | resolve in step 1.1 with Cloët before any public plot |
