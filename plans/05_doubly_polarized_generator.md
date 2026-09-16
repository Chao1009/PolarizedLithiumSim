# Plan 05 — Doubly Polarized e+⁶,⁷Li Event Generation ("polligen")

**Goal.** An event-level Monte Carlo of the *doubly polarized* eA process —
polarized electrons (helicity λ_e, P_e) on vector/tensor-polarized ⁶,⁷Li
(populations p_m over spin projections, arbitrary quantization axis) —
producing spin-labeled events with correlated scattered-electron,
spectator-fragment, and (optionally) hadronic final states, in HepMC3 for
the ePIC chain. This is the layer the fast simulation (plans/02, `fastsim/`)
deliberately postponed: the fastsim computes asymmetries and δA maps
*analytically per bin*; only an event generator can deliver (i) tagged
observables where the asymmetry is correlated with the spectator momentum,
(ii) reconstructed-level closure with realistic φ acceptance (make-or-break
for the cos 2φ gluonometry amplitude, plans/03 step 2.3.3), (iii)
pseudo-experiments with bunch-pattern / relative-luminosity / polarimetry
systematics, and (iv) samples an official ePIC e+Li production request can
be built on (plans/03 step 2.5).

**The build-vs-reuse decision in one paragraph.** The 2026-06-12
fetch-verified survey (plans/02, re-checked 2026-07-07) stands: **no public
generator produces polarized e+A events for A > 1** — BeAGLE is unpolarized;
DJANGOH's hadron polarization is nucleon-level and *longitudinal only*;
PEPSI/CLASDIS are polarized nucleon-level LO (longitudinal only); eHIJING,
Sartre, TOPEG, GCF are unpolarized. Nothing does tensor polarization or
transverse target spin at any A, and nothing correlates nuclear spin with
fragment kinematics. So we reinvent the wheel — but *only the wheel that
does not exist*: the polarized-nucleus vertex (spin-density matrix ×
structure-function master formula) and the spin-correlated cluster-spectator
sampler. We deliberately do **not** reinvent hadronization (PYTHIA),
radiative corrections (HERACLES/DJANGOH, step 1.4 workaround per
arXiv:2406.05591), or nuclear evaporation/backgrounds (BeAGLE, step 1.5).
The theory backbone exists and is current: Hoodbhoy–Jaffe–Manohar
(NPB 312:571) for inclusive spin-1, and Cosyn–Weiss for polarized tagged
DIS — PLB 799:135035, PRC 102:065204 (arXiv:2006.03033), and the 2026
spin-1 tagged-SIDIS pair arXiv:2603.23699/23700, whose stated purpose is
exactly "simulations of spectator tagging … at the EIC".

**Status legend:** ☐ todo ◐ started ☑ done

---

## 5.0 What "doubly polarized" must cover (requirements)

Spin configurations per observable — the generator must express all rows
from one master cross section and one run-plan bookkeeper:

| observable | e beam | ion spin | axis | terms needed |
|---|---|---|---|---|
| A∥ → g₁ᴬ (pol. EMC, ⁷Li & ⁶Li) | λ_e = ±, P_e | vector P_z | longitudinal | F₁,F₂ + λ_e P_z (g₁, g₂) |
| A⊥ (g_T access, control) | λ_e = ± | vector P_T | transverse | λ_e P_T cos(φ−φ_S), γ-suppressed |
| A_zz → b₁ (⁶Li) | unpol (any) | tensor P_zz | longitudinal (and ‖q, arXiv:2410.12764) | (2/3)a_m (b₁, b₂) |
| A_cos2φ → Δ (gluonometry) | unpol | tensor P_zz | **transverse** | c_m sin²θ_m Δ cos 2φ |
| tagged versions of all three | as above | as above | as above | tagged SFs vs (α_s, p_sT); spin ⊗ cluster-wave-function correlation |

Plus run-plan realism: bunch-by-bunch λ_e flips; ion fills with population
patterns (p₊, p₀, p₋) [⁶Li] / (p₃/₂ … p₋₃/₂) [⁷Li]; relative-luminosity
offsets between spin states at the 10⁻⁴ level; per-fill P_e, P_z, P_zz
values with polarimetry uncertainty (δP/P ≈ 3%, plans/04 #5). Every event
carries its spin labels so analysis-side estimators see exactly what the
experiment would.

Kinematic scope: inclusive DIS 10⁻⁴ < x < 1, Q² > 1 GeV², the three
reference energies of `beams.default_configs` per isotope; spectator
channels of `spectator.py` (⁶Li → α+d*, d+α*; ⁷Li → α+t*, t+α*), plus
evaporation-n ZDC tags left to BeAGLE.

## 5.1 Architecture: one kernel, three modes

```
                    ┌──────────────────────────────────────┐
                    │  physics kernel (new, the "wheel")   │
                    │  spin.py    ρ(m,m'), moments, axes   │
                    │  xsec.py    master formula, all SFs  │◄── polli_fastsim
                    │  tagged.py  cluster LFIA ⊗ spin      │    structure/polarized
                    └───────┬──────────────┬───────────────┘    (toy↔grid backends)
                            │              │
        Mode W (weights)    │              │    Mode G (native generator)
   per-event w(λ_e,m|x,Q²,φ,p_s)     sample spin state → (x,Q²,φ) →
   attached to ANY unpolarized       spectator k from n_m(k,k̂) →
   sample (BeAGLE, campaigns)        [optional PYTHIA on struck cluster]
                            │              │
                            └──────┬───────┘
                                   ▼
                     io_hepmc.py  spin-labeled HepMC3
                     (status-4 beams, 10-digit fragment PDG,
                      spin attributes) → abconv → npsim
```

- **Mode W (reweighting)** — the established EIC practice (DSSV
  arXiv:2007.08300, JAM arXiv:2105.04434, ECCE recipes NIM A 1056:168563).
  Works on any unpolarized sample whose *unpolarized* kinematics are right;
  injects any asymmetry whose SFs we can evaluate. Cheapest route to
  reco-level closure on the already-streamable official BeAGLE e+d samples.
  Limitation, understood and documented: for **Li tagged** states BeAGLE's
  spectator spectrum is evaporation-driven (C-12 n(k), no clusters) — a
  weight cannot fix a wrong unpolarized p_s distribution. Hence:
- **Mode G (native cluster-IA generator)** — samples the full doubly
  polarized tagged cross section; the only route to correct spin ⊗ spectator
  correlations for Li. Final states in three fidelity tiers:
  T0 = (e′, spectator, X as one pseudo-particle) — enough for every Phase-1
  FOM; T1 = + struck-cluster internal nucleon and its partner spectator(s)
  (α + p from d*; α + d/nn from t*) — double-tagging studies; T2 = + PYTHIA
  hadronization of the γ*–nucleon system — detector occupancy, e′ isolation,
  ZDC/B0 backgrounds (Phase 2).
- **Mode R (radiative corrections)** — not event-level at first: DJANGOH
  4.6.22 runs on an effective polarized nucleon (the arXiv:2406.05591 ³He
  workaround, already plans/02 step 1.4) provide multiplicative RC bands
  (*2026-08-28: the unpolarized collinear-ISR half is done in-repo —
  `polligen/radiative.py`, plans/07 WP4 — without DJANGOH; only the
  tensor-sector RC stays external, plans/04 #10*)
  per (x, y) that Mode W/G attach as optional weights. Tensor-observable RC
  remains uncharted (plans/04 #10) — the hook is there, the numbers await
  theory.

Language: numpy-vectorized Python, same style/test discipline as
`polli_fastsim` (it imports the existing SF backends and `farforward.py`);
4-vector-tier throughput ≥ 10⁴ ev/s is ample for Phase-1. PYTHIA tier runs
**natively** — the bindings do NOT ship in the eic-shell container (it has
the C++ library and headers only); PYTHIA 8.311 builds its own against the
analysis interpreter in two and a half minutes (2026-08-26,
`tools/pythia8/README.md`), which is how the standing 8 M-event HFS
production was made and what removes the fifo/LHE hand-off from step 5.D. C++ port only if an official
campaign demands it.

Layout: `evgen/polligen/{spin,xsec,tagged,sample,bookkeeping,reweight,io_hepmc,hadronize}.py`
+ `evgen/tests/` + `evgen/scripts/` (money plots 4–5, closure drivers).

## 5.2 The physics kernel (what is actually new)

1. **Spin-density matrix** ρ(m,m′) for J = 1 and J = 3/2 from
   (P_z, P_zz, quantization axis n̂(θ_S, φ_S)); diagonal in the n̂ frame
   (populations), rotated analytically into the lab/photon frame. Exposes
   vector/tensor (and for J = 3/2, neglected rank-3) moments. Unit-tested
   against analytic limits (pure states, unpolarized, HERMES-style
   P_zz = −2 λ₀-enriched fills).
2. **Inclusive master cross section** dσ(λ_e, ρ)/dx dy dφ. Spin-1: the
   Hoodbhoy–Jaffe–Manohar set {F₁, F₂, b₁…b₄, Δ, g₁, g₂} in the covariant
   Cosyn–Weiss classification (their inclusive limit fixes every sign and
   y-factor, including the transverse-vector λ_e P_T cos(φ−φ_S) g_T term
   the fastsim never needed). Spin-3/2 (⁷Li): rank-0/1 exact
   (F₁, F₂, g₁, g₂ with effective polarizations P_p = 0.866, P_n = −0.037),
   rank-2 as scenario inputs (b₁-analog shapes ← plans/04 #9; the complete
   spin-3/2 basis is itself a publishable theory note — engage Cloët/Cosyn;
   the ask, with its conventions and the P_zz ↔ T map, is written out in
   `docs/note_7li_theory_questions.md` §1),
   rank-3 dropped. Backends: existing `ToyF2/PartonF2`, `ToyG1/PartonG1`,
   scenario b₁/Δ curves — zero duplication; g₂ = g₂^WW default.
   Consistency gate: with ρ → (P_z only), reproduce `asymmetries.a_parallel`;
   with tensor-only ρ, reproduce `asymmetries.azz`/`a_cos2phi` bin by bin.
3. **Tagged cross section (cluster light-front IA ⊗ spin).** Cosyn–Weiss
   deuteron formalism transplanted to two-cluster Li:
   - ⁶Li(1⁺) = α ⊗ d, relative L = 0 (+ small L = 2): to leading order the
     ⁶Li spin *is* the embedded-deuteron spin (S-wave CG is trivial), so the
     α-tagged vertex = Cosyn–Weiss polarized deuteron with (a) the α–d
     momentum density replacing |ψ_d|² for the spectator, (b) the struck
     "nucleon" replaced by the polarized deuteron (its own g₁/b₁ enter —
     this is precisely the embedded-b₁ observable), (c) P_d(⁶Li) = 0.87 and
     the L = 2 admixture as the correction band.
   - ⁷Li(3/2⁻) = α ⊗ t, relative L = 1: Clebsch–Gordan structure
     |3/2 m⟩ = Σ ⟨1 m_L ½ m_t|3/2 m⟩|1 m_L⟩|½ m_t⟩ delivers, with *no new
     parameters*, both the triton polarization per m-state (→ P_p ≈ 0.87
     forward limit; VMC refines to 0.866) and the m_L-dependent **angular
     distribution of the tagged α**: |Y₁^{m_L}(k̂)|² correlates k̂ with the
     spin axis. Two free consequences to exploit: tagged tensor-type
     asymmetries for a spin-3/2 nucleus (never computed by anyone), and the
     tagged-α polar/azimuthal moments as an **in-situ alignment polarimeter**
     for the ⁷Li beam (our inference — FSI caveats apply; polarimetry is
     open question plans/04 #5, so even a cross-check is valuable).
   - Radial inputs: the existing two-parameter Hulthén/P-wave forms of
     `spectator.py` now, **VMC two-cluster overlaps (Wiringa, ANL) as the
     upgrade** — same interface, and the S/D (⁶Li) and P (⁷Li) amplitudes
     with their m-dependence are exactly what VMC tables provide. The tail
     dominance found in the e+d control (BeAGLE pT tails 2–13× the Hulthén
     model) is the driving systematic → always run the β = 0.20/0.30/0.40
     band until VMC lands.
     ☑ *2026-09-15 (run 19b; restates the 2026-08-28 box it replaces, which
     counted four call sites where there are three). `--beta-band` is now on
     `tagged_polarimetry_7li.py`, `nearbeam_aperture_scan.py` and
     `money_tagged_azz.py`'s acceptance-folded panel, **default off** and on its
     own `_betaband` stem, so every published PNG regenerates to the same md5
     (a449a360… / c091b7c0… / 9f13e2ba… against `df32780`) and the printed text
     byte for byte. The fourth site of the old box, `nearbeam_reach_gain.py`,
     carries **no cluster β at all** — an AST walk finds no `beta` attribute or
     keyword and no spectator import; it runs the coherent intact-⁶Li chain, whose
     model band is the 40–60 GeV⁻² one around B = 50 — so a flag there would be a
     no-op and the measurement is recorded in its docstring instead. **What the
     band is worth, measured:** the ⁶Li α tag spans **×3.35–3.44** over
     β = 0.20–0.40 at the Yellow Report high-acceptance envelope (0.0082 / 0.0177 /
     0.0283 at 5 × 41, 0.0076 / 0.0162 / 0.0254 at 10 × 100, 0.0114 / 0.0247 /
     0.0391 at 18 × 275) and **×1.35–1.52** at the tagging optics (0.2591 / 0.3159 /
     0.3490, 0.1695 / 0.2235 / 0.2584, 0.2298 / 0.2920 / 0.3298); the ⁷Li α tag
     spans ×1.01–1.04 and moves the *other* way (0.9840 / 0.9684 / 0.9509 at
     5 × 41 on the Yellow Report envelope, ×1.03), because
     it is caught by the momentum window and a harder spectrum only spills a little
     of it out. On money plot 4 the band moves the **rate** by ×3.08 (YR) and ×1.45
     (tagging) and the **asymmetry** by ≤ 0.07 absolute — acceptance-weighted truth
     at k = 0.325 GeV/c −0.889 / −0.871 / −0.824 and +0.185 / +0.181 / +0.171 — so
     the tagged A_zz is the robust half of that figure and the α-tag fraction is
     not. **The band is one-sided IN β and not in the tag:** the 2026-08-26 e+d
     control found that no β in a two-parameter Hulthén form reproduces BeAGLE's
     p_T tail (2–13×), so the true short-range scale sits at or above 0.40 — read
     that end, not the middle — while what the tag then does is upward for ⁶Li and
     downward for ⁷Li. Independent evidence, from the VMC upgrade the band stands
     in for in the sibling generator LiPolGen (`README.md`,
     `--cluster-wave {hulthen,vmc}`): the ⁶Li α-tag fraction goes 0.0264 → 0.0348
     at 10 × 99.5 on the YR high-acceptance envelope and 0.2551 → 0.2486 at the
     tagging optics — **opposite directions at the two optics**, which is why the
     band is quoted per optics and never as one factor. Reports 3 §5 / Table 6 and
     4 §2.1 no longer publish the ⁶Li α tag bare.*
     *2026-09-15: the VMC upgrade now also carries a sign. Since `tagged.py`
     applies the i^L phase of the partial-wave expansion (plans/00 run 19), the
     sign of every tagged ⁶Li A_zz follows sign(ψ₂/ψ₀) of the α–d D radial,
     which the two-parameter form cannot predict; the model takes it
     deuteron-like, and the sibling generator's VMC α+d overlap supports that
     over 0.134–0.444 GeV/c only — it reads −1 below the α–d S node at
     0.134 and −1 again above the α–d D node at 0.444, two reversals a
     node-free Hulthén form cannot carry, with 27 / 21 / 26 % of the
     Yellow-Report-accepted sample above the upper node and 41 / 28 / 37 %
     of the tagging-optics sample below the lower one. plans/04 #15.*
   - FSI: not modeled at first (IA). Quote tagged asymmetries at low
     spectator virtuality / small |t′| where pole dominance holds
     (Cosyn–Weiss FSI: PRC 97:035209); flag α-spectator FSI as a theory
     engagement item (new plans/04 #16).
4. **Sampling.** Spin config from the bunch bookkeeper → (x, Q²) by 2-D
   inverse-CDF of the unpolarized σ (grid-backed), φ uniform → accept-reject
   on the (small, bounded) polarized modulation → unweighted spin-labeled
   events; tagged mode adds k ~ n_m(k, k̂) per spin state. Weighted mode
   (all spin states per event, vector of weights) supported for
   FOM-efficiency studies. Reproducibility: fixed rng streams per
   (run, bunch) — same discipline as `spectator.sample_k`.

## 5.3 Steps

### Step 5.A ☑ Kernel + inclusive sampler (done 2026-07-13)
`spin.py`, `xsec.py`, `sample.py`, `bookkeeping.py` (+ `estimators.py`);
tests: ρ-matrix moments, master-formula ↔ `asymmetries.py` bin-wise
identity, estimator closure (full pseudo-experiments with bunch
patterns and a 10⁻⁴ relative-luminosity offset — first systematics number).
**Deliverable:** inclusive doubly polarized e+⁷Li/⁶Li pseudo-experiments
reproducing every Phase-1 analytic FOM map (δA∥, δA_zz, δA_cos2φ) —
the fastsim↔generator closure that certifies both.
**Done (see `evgen/README.md`):** 35 tests, all four §5.4 inclusive gates
pass (identities at rtol 1e-12 on toy + CT18/NNPDFpol backends; FOM-map
closure over ~65 x-bins per isotope, `evgen/closure_fom_{6,7}Li.png`);
rel-lumi systematics: bias(A_zz) = −(2/3)δ/P_zz, bias(A∥) = δ/(2P_eP_z),
≈1×10⁻⁴ at the reference δ = 10⁻⁴ — removed exactly by lumi-corrected
estimators. Spin-temperature (max-ent) fills added after positivity ruled
out naive (P_z, 0-tensor) spin-3/2 fills.

### Step 5.B ☑ Tagged mode (2–3 weeks, the core novelty)
`tagged.py`: spin-correlated (e′, spectator) events for the four Li
channels + deuteron/³He controls. Validations in §5.4. **Deliverables:**
☑ *2026-08-28: landed 2026-07-13 in commit e1b7547 — `tagged.py`
(`TaggedChannel`/`TaggedModel`/`TaggedSampler`/`boost_spectator`/`rp_accepted`)
with 17 tests in `tests/test_tagged.py`; both deliverables below re-run today.
What remains of the step: `tagged.py` wraps three channels (⁶Li α, ⁷Li α, the
deuteron control), not "the four Li channels + deuteron/³He controls" — the ⁶Li d,
⁷Li t and ³He p tags have no `TaggedChannel`.  The §5.4 deuteron-limit gate is
now met quantitatively and needed no digitization: Cosyn–Weiss II page 35 gives
the closed form (their Eq. 6.12) that FIG. 13 only illustrates, and
`tests/test_tagged.py::test_cosyn_weiss_tensor_gate` pins our model against it.
**Restated 2026-09-15** (plans/00 run 19): the gate is now the identity itself,
not a five-digit ratio.  With the channel's own radial tables as (f₀, f₂),
max |A_zz^wf − Eq. (6.12)| is 8.9e-16 on the deuteron control and 1.1e-15 on
⁶Li over the whole (k, cos θ_k) grid, and the (1 − 3cos²θ_k) factorization holds
to 6.0e-14; the normalization map is **A_T∥ = +1 × A_zz^wf**, because CW's
factor −2 is the angular factor at θ_k = 0 that A_zz^wf already carries.  The
old −2 map goes, and with it the sentence that Report 4's −0.48 on the 90°
curve "is CW's +0.96": on the corrected wave function that curve reads +0.92.
So does the claimed envelope peak at k = 0.3098 GeV — that is where the toy
|f₂/f₀| = 1/√2, i.e. Eq. (6.14)'s *minimum*, on a Hulthén deuteron whose f₂/f₀
never reaches √2 at all (max 1.2866 at the k = 1.2 GeV/c grid edge).  TABLE II
is instead pinned where it can be met, on the AV18 deuteron
(`tagged.av18_deuteron_channel()`, `test_cosyn_weiss_table_ii_on_av18`):
A_zz^wf = −1.937 at the cell nearest θ_k = 0 and −1.998 on a fine near-axis grid
where n_{±1} has its node, +0.999 at θ_k = 90°, both at the f₂/f₀ = √2 crossing
k = 0.298 GeV/c (CW: "k = 0.30 GeV"), and +0.967 at k = 1.0 GeV/c with the
f₂/f₀ = −1/√2 crossing at 1.03 GeV/c.*
- **Money plot 4:** tagged tensor asymmetry A_zz^tag(p_s) for the α-tagged
  embedded deuteron in ⁶Li, folded with `farforward.py` acceptance and both
  optics — *first tagged spin observable for any A > 2* (extends the
  four-gap list of plans/00). Check the Cosyn–Weiss O(1) asymmetry at
  p_s ≈ 300 MeV/c survives the RP pT-tail acceptance that dominates ⁶Li.
  ☑ *2026-08-28 (plans/09 B2). The figure was answering the question on the
  legacy proton-derived 73/164 µrad optics that commit b9d2e82 retired —
  its printed 5.0% / 2.5% were the two totals of that pair — and it is now
  per configuration (`--config {0,1,2}`) and per optics (`--optics`,
  default the Yellow Report high-acceptance optics plus the tagging optics
  with its luminosity fraction), routed with the spectator's own lab
  azimuth (`tagged.boost_spectator` returns `phi_spec`) against the
  rectangular envelope; without that azimuth the tagging optics read 0.51
  against the 0.30 with it, on the levers in force when the omission was
  found (2026-08-28). The answer to the plan's question is **yes, and at the
  published optics only there**: the α tag is 0.0241 at the Yellow Report
  optics against 0.2542 at the tagging optics — at L/L_HA = 1/12.8 an 18%
  cost in tagged events per year here, and a 1.8× / 1.3× *gain* at the
  other two configurations — but the median accepted spectator momentum is
  0.32 GeV/c with nothing below k = 0.15 GeV/c, against 0.18 GeV/c with
  36% below it. The tagging optics turns a one-point measurement into a
  curve. Two by-products. The overlay was wrong: an analytic curve at
  θ_k = 90° drawn over a sample the acceptance sculpts to
  ⟨|cos θ_k|⟩ = 0.797 (the off-rigidity slice, longitudinal) or 0.395 (the
  near-beam tail, transverse), so the swing between the two optics at
  k ≈ 0.3 GeV/c — 1.06 wide on the corrected wave function, ±0.5 on the
  figure as it stood then — was the envelope; the right panel now carries the
  acceptance-weighted truth per optics, which tracks the markers to 1–2σ,
  with the 90° curve kept as a labelled reference. And the two k spectra
  are reconciled: `tagged.TaggedSampler`'s ⟨k⟩ = 0.122 GeV and 2.5% below
  R = 0.95 against `spectator.spectator_lab_kinematics`'s 0.107 GeV and
  1.5% is the **D wave and nothing else** — the S-wave radials are
  identical, the D wave has ⟨k⟩ = 0.278 GeV and P_D = 0.0867 — and since
  A_zz^tag vanishes identically at P_D = 0, the tagged observables must be
  quoted on the S + D spectrum and the acceptance table on the S-wave one.*
  *2026-09-15 (plans/00 run 19): the numbers above are the corrected ones.
  `TaggedModel._amp2_table` summed the partial waves without the i^L phase, so
  the α–d S/D interference carried the wrong relative sign, and every A_zz of
  this figure flipped with it: at k = 0.325 GeV/c the folded markers now read
  −0.843 (acceptance-weighted truth −0.871) at the Yellow Report optics and
  +0.215 (+0.181) at the tagging optics against the θ_k = 90° curve's +0.922,
  where the published figure read +0.491 (+0.455), −0.066 (−0.095) and −0.482.
  The tag acceptances moved only in the fourth digit, which is seed noise in
  either build; the k-marginals, ⟨k⟩ = 0.122 GeV and the D-wave reconciliation
  are exactly invariant, because the interference cancels in the angular
  integral.  `evgen/money_tagged_azz_6Li.png` was regenerated.*
- ⁷Li α-tag: tagged A∥ (polarized-EMC companion on the quasi-free triton)
  + the tagged-α angular-moment polarimetry curve vs P_zz.
  ☑ *2026-08-28: both panels are in `scripts/tagged_polarimetry_7li.py`
  (`evgen/tagged_polarimetry_7Li.png`, regenerated after the γ-matched energy
  fix) — ⟨P₂(cos θ_k)⟩ on the analytic −T/5 line with RP-folded markers, and
  A∥^tag(x) on the quasi-free triton against the D(y)·g₁ᵗ/F₁ᵗ overlay; the
  analytic backing is gated by `tests/test_tagged.py`. It takes `--config`
  and `--optics` and defaults to the configuration's Yellow Report
  high-acceptance envelope, with the tagging optics beside it (plans/09 B3);
  the retired legacy pair is reachable only by asking for it.*

### Step 5.C ◐ Mode W reweighter + reco-level closure (1–2 weeks)
`reweight.py` driving the official BeAGLE e+d sample (xrootd, already
streamed for the control study): inject A∥(g₁d) and A_zz(b₁d), extract with
the analysis estimators through the existing conversion chain, verify pulls.
This is the ECCE-style pedigree demonstration on the nucleus where BeAGLE
*is* right, and it exercises the identical machinery later pointed at e+Li
BeAGLE samples for purity studies (evaporation background stays BeAGLE's
job; plans/02 step 1.5.4).
  ☑ *2026-09-15: the reweighter and the generator-level closure are in
  `polligen/reweight.py`, gated by `tests/test_reweight.py` (16 tests).
  `ModeWReweighter` evaluates the kernel's own (w_avg, a₁, a₂) triple at
  each external event's (x, Q², φ) — not at a grid cell — and, because an
  unpolarized sample has flat φ, resampling it with probability ∝ W
  reproduces the polarized cross section including its φ modulation;
  `bookkeeping`'s run plans and `estimators`' counting estimators then read
  the spin-labelled categories unchanged. Injected and recovered on the
  official `BeAGLE1.03.02-3.1/eH2/en/9x130` sample, 20 000 events streamed
  through `tools/analysis/dump_spectators.py` (19 935 inside the
  `fom.Scenario` window): over 1000 pseudo-experiments the pull mean is
  +0.005 on A∥ and −0.040 on A_zz with widths 0.992 and 1.007, against the
  gate |μ| < 0.15 and 1.00 ± 0.10; the same gate on the synthetic in-memory
  pool the tests use gives −0.074 / +0.018 and 1.010 / 1.007, and it holds
  with the events themselves resampled and not only their Poisson counts.
  The φ channel is closed too: a cos 2φ′ amplitude injected into a flat
  pool comes back through `estimators.cos2phi_fit` to 0.5% of itself. Two
  by-products. The streamed sample's kinematics are re-formed from the
  per-event beam and scattered-electron rows against a supplied
  beam-electron four-vector — the dumper writes no status-4 lepton — and
  agree with BeAGLE's own `trueX`/`trueQ2`/`trueY`/`leptonphi` to a median
  −8.3×10⁻⁴ / +7.5×10⁻⁵ / +5.9×10⁻⁴ / 0 and a 68th percentile of
  1.9×10⁻³ / 9.8×10⁻⁴ / 1.4×10⁻³ / 1.3 mrad, with a 0.27% tail where the
  dumper's highest-energy electron is not the scattered one. And the
  container note of docs/reproduction_manual.md §5.2 was wrong: the newer
  `eic_xl-nightly` image reads these tree files, it is ROOT's relative
  module-map path that needs `singularity exec --pwd /opt/local/lib/root`.*
  ◐ *Not closed: the reco-level half of this step. The pulls above are
  truth-level — the external sample supplies the kinematics and the
  hadronic final state, but nothing between the generator and the
  estimator. Putting a reweighted sample through abconv → npsim → EICrecon
  needs the HepMC3 writer of step 5.D (`io_hepmc.py`), which does not
  exist; `reco.py`/`recopseudo.py` fold resolutions analytically rather
  than reading a reconstructed file. `draw_category` already returns the
  Mode-G event-dict shape so that a writer can consume it unchanged.*
  **Disposition 2026-09-15 (run 19b):** the step splits cleanly and the two
  halves have different owners. The **generator-side half is delivered** —
  `polligen/reweight.py` and its 16 tests are the ☑ clause above, and nothing
  in this repository is owed for it. The **reco-level half stays ◐**, and its
  blocker has moved rather than cleared: the writer it waits on is step 5.D's,
  which the sibling generator LiPolGen now ships (see 5.D below), so what is
  left here is not a writer but a *reader* — `abconv → npsim → EICrecon` output
  fed back to `estimators.py` — and the EICrecon leg is open on both sides.

  ◐ *2026-09-16: the reader exists and the leg has been RUN, once, bounded.
  100 events of the official `BeAGLE1.03.02-3.1/eH2/en/9x130` sample carry a
  Mode-W weight of `polligen.reweight` for the `azz0` tensor third
  (`tools/analysis/modew_beagle_hepmc.py`, HepMC3 ASCII with the weight NAMED
  in the weights vector), go through npsim and EICrecon
  (`tools/fullsim/modew_chain.sh`, `epic_craterlake_18x275.xml` on both legs)
  and come back as 100 reconstructed events read with uproot
  (`tools/analysis/modew_reco_readback.py`); docs/reproduction_manual.md §5.2
  and §5.3 carry the commands. Four things the leg taught, none of them a
  number this repository publishes. **(a) npsim keeps only the NOMINAL HepMC3
  weight.** A weight appended after the sample's own `default` reaches the
  reconstructed file as a constant 1.0 and the weight NAMES do not survive at
  all; written at index 0 it arrives on 100/100 events with 100 distinct
  values in `EventHeader.weight`, to the generator table's own 10-digit
  printing precision. The writer therefore makes the Mode-W weight the
  nominal one by default — correct for a reweighted sample, the BeAGLE files
  being unweighted. **(b) The event order survives the chain**, to a median
  5.1×10⁻⁷ in x against BeAGLE's own `trueX` read in file order, so a weight
  a tool drops can always be re-attached by index. **(c) EICrecon's inclusive
  kinematics need a status-4 PROTON.** Every reconstructed method —
  `InclusiveKinematicsElectron`, `Sigma`, `DA`, `JB`, `ESigma`, `ML` — is
  empty in 100/100 events because `MCBeamProtons` is empty: the `eH2/en`
  files record the struck NEUTRON as the status-4 ion beam.
  `InclusiveKinematicsTruth`, which does not use it, is filled in 100/100,
  and the scattered electron itself is found (`ScatteredElectronsTruth`
  100/100), so the reco-level x and Q² have to be formed by the reader, which
  they are — weighted, and at a median residual of −0.052 in x and −0.0050 in
  Q² over 99 events, the x residual being the electron method's 1/y tail
  (−0.44 below y = 0.1 against −0.013 above it). This is the
  inclusive-kinematics counterpart of the far-forward species blocker and it
  will bite a ⁶Li beam row the same way. **(d) npsim seeds from the clock
  unless told otherwise**, and at 100 events the shower seed moves every
  reconstructed number by more than the smoke measures, so the chain script
  fixes a seed by default. What is still NOT closed is the pull closure at
  the reconstructed level: 100 events is a smoke test, the estimators want
  O(10⁴–10⁵), and `estimators.py` has not been run on a reconstructed file.*

### Step 5.D ☑-by-sibling Final states + HepMC3 + chain smoke test (2 weeks)
Tier T1 (cluster-internal nucleon + partner spectators; t* remnant → d or
nn per the triton wave function — crude, flagged), `io_hepmc.py` (ASCII
HepMC3; status-4 beams so `abconv` accepts it; 10-digit fragment PDG;
spin labels as named attributes — no HepMC3 convention exists for ion spin
states, so define one and propose it to the ePIC MC group, new plans/04
#17), tier T2 PYTHIA attachment — natively, and the HFS half of it is
already done and in production use (`polligen/hfs.py` on the 8 M-event
sample of `tools/pythia8`); what T2 still owes is the attachment to the
*tagged* final state, not the generator. 100-event
abconv → npsim → EICrecon smoke per plans/03 step 2.1.4 (⁶Li can proxy the
existing d/⁴He beamline configs; ⁷Li needs the new field maps — already
plans/03 step 2.1.2).
  ☑-by-sibling *2026-09-15 (run 19b). **This repository will not duplicate the
  writer.** The three deliverables this box names — the HepMC3 writer, the
  ion-spin attribute convention and the abconv → npsim chain gate — are shipped
  by the sibling generator `LiPolGen`, and `io_hepmc.py` is therefore withdrawn
  rather than deferred. Evidence, read there: `include/lipolgen/hepmc_writer.hpp`
  / `src/hepmc/hepmc_writer.cpp` write `lipolgen::Event` records as HepMC3
  Asciiv3; `docs/HEPMC3_CONVENTION.md` is the named-attribute schema and says in
  its own first paragraph that it **is** `plans/04` open item **#17** ("no such
  convention exists upstream in HepMC3 or in the EIC software stack, so LiPolGen
  defines one here and is prepared to propose it to the ePIC MC group") — so #17
  is answered by a document, and what remains of it is the proposal to the ePIC
  MC group, not the schema; `README.md:130` "ePIC chain gate passed: HepMC3 →
  `npsim` (direct, and via `abconv -p ip6_hiacc_100x10`)", with
  `docs/OPEN_ITEMS_SOLUTIONS.md:38` recording 10/10 events through `npsim` both
  ways. **The EICrecon leg is still open on both sides** — the sibling's
  `docs/PHYSICS_CHANNELS.md:954` says "the `abconv` → `npsim` → EICrecon smoke
  test is not run from this repository" — and it is the same leg step 5.C's
  reco-level half waits on. The T1 remnant tier and the tagged-final-state half
  of T2 are a separate question from the writer and stay where §5.5's risk row
  puts them.*

  ☑ *2026-09-16 on the "100-event abconv → npsim → EICrecon smoke" this box
  names: it has been run, from THIS repository, on an external sample rather
  than on a generated one — 100 official BeAGLE e+d events with a Mode-W
  weight attached, through npsim and EICrecon, read back (step 5.C's
  2026-09-16 clause, docs/reproduction_manual.md §5.3). Two notes on the
  abconv leg, which is the one this box is explicit about. The official EVGEN
  files are ALREADY afterburned — `GenRunInfo ab_afterburner_is_used = 1`,
  `ab_crossing_angle = 0.025`, and the beam rows carry −25 mrad — so abconv is
  a CHECK there (`abconv -p <preset> --exit-ca` exits 0 saying so) and not a
  transform; and abconv's auto preset ABORTS on a 9 × 130 file ("9x130 is not
  a valid energy combination!!", SIGABRT), its own nearest being the
  approximate `eD 10x130 GeV/n`. A `--ab-off` pass-through does preserve the
  named weight, 100 events in and 100 out. The remaining gap in this box is
  unchanged: the smoke ran on a BeAGLE sample, not on the sibling generator's
  own HepMC3 output.*

### Step 5.E ◐ Physics production + write-up (2 weeks + ongoing)
Regenerate all money plots from generator pseudo-experiments (statistical
FOMs now include acceptance × estimator effects); tagged-FOM table
(efficiency × purity × dilution per channel); short generator note —
"first polarized eA event generator" is itself one of the publishable
firsts, aimed with the money plots at the INT program (Mar 22 –
Apr 2, 2027). Upgrades stay behind interfaces: VMC overlaps, spin-3/2 rank-2
SFs, RC tables.  *2026-08-28: the digitized CBT/TMT and b₁ theory curves of
step 1.2 landed — `fastsim/polli_fastsim/data/`, drawn by default, with the
old shapes behind `--emc-mode constant` and `--transfer legacy`.*

— *superseded (2026-08-28) for the generator note: plans/07 §7.0 D1 ruled that a
software paper fits CPC/EPJ C, not the letter, and deferred it to a companion
"later"; the INT vehicle is now the PLB letter with the money plots, and the
INT-facing circulation note already exists as docs/note_cos2phi_coherent_6Li.md;
since 2026-08-29 docs/note_7li_theory_questions.md carries the ⁷Li theory asks
to the same audience.*

**Disposition 2026-09-15 (run 19b).** The step has three deliverables and they
are now in three different states. The **write-up half is superseded** by
plans/07 §7.0 D1, above — there is no generator note to write. The **money-plot
regeneration half** is carried by the individual figure boxes of §5.3 and by
plans/07 WP3, not by this step. The **tagged-FOM table is delivered**, below.

#### 5.E tagged-FOM table (efficiency × purity × dilution per channel) ☑ *2026-09-15*

`fastsim/scripts/tagged_fom_table.py` (new). Efficiency is the
`tagging_acceptance.py` machinery — `spectator.spectator_lab_kinematics` folded
with `farforward.acceptance_summary`, tagged = 1 − lost, the definition Report 3
Table 6 tabulates — at 4 × 10⁵ spectators per cell over the β = 0.20/0.30/0.40
band. Dilution is the tagged spin model's own,
`TaggedModel.tensor_dilution()` where the channel spin allows a rank-2 moment
and `.vector_dilution()` for ⁷Li, whose S_c = ½ has none; the script **asserts**
that the ⁶Li row reproduces `polarized.b1_li6_from_deuteron(1.0)` rather than
transcribing 0.921949 a second time, which is the drift plans/08 D9 exists to
prevent.

| channel | rank | D_model | footing | **D_published** | observable |
|---|---|---|---|---|---|
| ⁶Li α (embedded d) | 2 | 0.9219490 | 2/6 per-nucleon (`LI6_B1_PER_NUCLEON`, D9) | **0.3073163** | tagged A_zz, money plot 4 |
| ⁷Li α (quasi-free t) | 1 | 1.0000000 | P_p(t) = 0.86 (`tagged.TRITON.eff_pol_p`, per nucleon = whole-triton, §5.4) | **0.8600000** | tagged A_∥ |
| d–p control | 2 | 0.9594889 | 1 (quoted on the deuteron itself) | **0.9594889** | Cosyn–Weiss tagged A_zz |

| channel | optics | ε (β = 0.30) | ε over the β band | purity, \|t\| window | purity, incoherent bkg | ε × D | ε × D² |
|---|---|---|---|---|---|---|---|
| ⁶Li α | YR high-acceptance | 0.0168–0.0255 | 0.0079–0.0395 | n/a | **unavailable (FLUKA)** | 0.0052–0.0078 | 0.0016–0.0024 |
| ⁶Li α | tagging | 0.2236–0.3150 | 0.1691–0.3486 | n/a | **unavailable (FLUKA)** | 0.0687–0.0968 | 0.0211–0.0298 |
| ⁷Li α | YR high-acceptance | 0.9688–0.9755 | 0.9511–0.9873 | n/a | **unavailable (FLUKA)** | 0.8332–0.8389 | 0.7165–0.7215 |
| ⁷Li α | tagging | 0.9874–0.9941 | 0.9794–0.9970 | n/a | **unavailable (FLUKA)** | 0.8492–0.8550 | 0.7303–0.7353 |
| d–p | YR high-acceptance | 0.9549–0.9564 | 0.9461–0.9713 | n/a | **unavailable (FLUKA)** | 0.9163–0.9176 | 0.8791–0.8805 |
| d–p | tagging | 0.9549–0.9564 | 0.9461–0.9713 | n/a | **unavailable (FLUKA)** | 0.9163–0.9176 | 0.8791–0.8805 |

Ranges are over the three beam configurations. `ε × D` is this step's literal
product with the purity factor left open; `ε × D²` is the statistical figure of
merit, since δA divides by D (`fom.Scenario.analyzing_power`) and the count goes
as ε.

**The two purity columns are the honest half of the table.** The kinematic one is
`n/a` by construction: a spectator tag is selected by a rigidity/angle window and
not by a \|t\| fit, and the only \|t\|-window purity this project holds — 80–99%
incoherent rejection from the e+Pb coherent-J/ψ study arXiv:2108.01694
(PRD **104** 114030), plans/06 and plans/07 risk row 9 — belongs to the
**coherent** intact-⁶Li recoil, a different channel, and is recorded here as the
reference it is. The incoherent-background one is **unavailable and FLUKA-gated**:
BeAGLE links FLUKA, whose licence is personal and per-user, so no A = 6, 7
breakup sample exists anywhere in this project (`docs/reproduction_manual.md`,
plans/08 D5). It is left empty rather than guessed, which is why the product is
reported as ε × D and not as one number that would look complete — and why **no
Report 3 row is proposed for this table**: a published row needs every column
sourced, and one of them is not.

Total ≈ 7–9 focused weeks to 5.E; 5.A+5.C alone (≈ 3–4 weeks) already
upgrade every Phase-1 FOM to pseudo-experiment grade.

— *superseded (2026-08-28): the FOM upgrade arrived without 5.C —
`polligen/recopseudo.py` (plans/07 WP3) took the cos 2φ FOMs to reconstructed-level
pseudo-experiment grade and `hfs.py` replaced the 25% hadronic-y stand-in with a
PYTHIA-backed response; 5.B is done, so what is left of the 7–9 weeks is
5.C + 5.D + 5.E.*

— *restated 2026-09-15 (run 19b), on the three dispositions above: **what is left
of the 7–9 weeks is the reco-level half of 5.C, and nothing else in this
repository.** 5.D is delivered by the sibling generator and its writer is
withdrawn here rather than deferred; 5.E's write-up half is superseded by
plans/07 §7.0 D1 and its tagged-FOM table is delivered above, leaving only the
FLUKA-gated purity column, which is external. The residue is one reader —
`abconv → npsim → EICrecon` output fed back to `estimators.py` — and the EICrecon
leg is open on both sides, so it is gated on the ePIC chain and not on effort
here.*

## 5.4 Validation matrix (gates, in order)

| gate | reference | pass criterion |
|---|---|---|
| ρ moments, all axes | analytic | exact (unit test) — ☑ *2026-08-28: `tests/test_spin.py`, 16 tests; rotations and population round trips for J = 1 and 3/2 at atol 1e-12* |
| master formula, vector/tensor sectors | `asymmetries.py` | bin-wise identity (toy + grid backends) — ☑ *2026-08-28: `tests/test_xsec_identity.py`, 15 tests at rtol 1e-12, the grid half actually running on CT18NLO + NNPDFpol11_100* |
| pseudo-experiment estimators | `fom.py` maps | δA agree within trial statistics; pulls unbiased — ☑ *2026-08-28: `tests/test_pseudoexp.py` and `scripts/closure_fom.py` (~65 x-bins per isotope); means unbiased against the σ-weighted truth, spreads within 15% of the three analytic errors* |
| φ-modulation recovery | injected Δ scenarios | amplitude unbiased with uniform *and* holey φ acceptance — ☑ *2026-08-28: `test_cos2phi_fit_unbiased_with_holey_acceptance` removes two asymmetric φ sectors; the fit is unbiased at 5×SE while the naive moment is biased by >10×SE, so the gate is not vacuous* |
| deuteron limit of tagged mode | Cosyn–Weiss arXiv:2603.23700 Eq. (6.12), Eqs. (6.13)–(6.14), TABLE II (p. 35) | ☑ *2026-09-15 (restates the 2026-08-28 row, which was met on the wrong normalization map): the gate is the **identity**, `test_cosyn_weiss_tensor_gate` — with the channel's own (f₀, f₂), max \|A_zz^wf − Eq. (6.12)\| = 8.9e-16 (deuteron) and 1.1e-15 (⁶Li) over the whole grid, the (1 − 3cos²θ_k) factorization to 6.0e-14, the map A_T∥ = +1 × A_zz^wf (CW's −2 is the θ_k = 0 angular factor A_zz^wf already carries), and the whole curve inside their [−2, 1], approached (< −1.9, > 0.99) but not attained on either Hulthén toy. TABLE II is a separate gate on the **AV18** deuteron, `test_cosyn_weiss_table_ii_on_av18`: −1.937 at the cell nearest θ_k = 0 (−1.998 on a fine near-axis grid, at the n_{±1} node) and +0.999 at 90°, both at the f₂/f₀ = √2 crossing k = 0.298 GeV/c against CW's 0.30, and +0.967 at k = 1.0 GeV/c. The old row's 0.99940 / 0.3098 / +0.9997 / −2.000 are retired: 0.3098 was Eq. (6.14)'s minimum read as Eq. (6.13)'s maximum on a toy whose f₂/f₀ never reaches √2. The FIG. 13 panels themselves are in light-front variables (α_p, p_pT) the sampler does not carry, so they are a comparison, not a gate* |
| unpolarized spectator spectra | official BeAGLE e+d via `ed_control_analysis.py` | bulk agreement; tail differences documented as the model band — ☑ *2026-08-28: run on the BeAGLE 1.03.02-3.1 eH2 9×130 sample; routing agrees to better than 2 points, but no β reproduces the p_T tail, so the difference is carried as a one-sided upward band rather than the symmetric one this row assumed* |
| forward limit of tagged ⁷Li | **whole-nucleus ⁷Li VMC sums**: P_p = +0.866, P_n = −0.037 (JLab PR12-14-001 Eq. 29, rounding Wiringa *et al.* PRC **89** (2014) 024305 Table I's +0.868 / −0.038; `beams.LI7` stores them divided by Z = 3 and N = 4) | ☑ *2026-09-15 (restates the 2026-08-28 row, which compared two different objects). Both halves are now asserted and the neutron gap is recorded rather than left open. The **proton** half, `test_li7_triton_polarization_forward_limit`: the model's whole-nucleus P_p is P_t(M = 3/2) × Z_t × `TRITON.eff_pol_p` = 1 × 1 × 0.86 = **+0.86**, within 0.02 of +0.866. The **neutron** half, the new `test_li7_neutron_forward_limit_both_footings`, which prints the gap: `TRITON`'s −0.028 is **per nucleon** — the isospin mirror of Bissey's ³He (PRC **65** 064317), per-nucleon like every `Ion` slot since plans/08 D7 — and is **not** the gate's whole-nucleus −0.037. On the whole-nucleus footing the model gives P_t(M = 3/2) × N_t × (−0.028) = 1 × 2 × (−0.028) = **−0.056**, the α spectator being spin-0 and contributing exactly zero; per nucleon that is **−0.014** over ⁷Li's N = 4 and **−0.028** over the triton's own two neutrons. The gap to the ab initio value is **−0.019 (×1.5135)** against −0.037 and −0.018 (×1.4737) against Table I's −0.038, and it is a **known model difference, not a band**: the α + t decomposition puts all the ⁷Li neutron spin on two neutrons where VMC spreads it over four correlated ones, and this channel is a lone L = 1 wave with no D-state admixture to widen a tolerance around — which is why the proton half hid the distinction (Z_t = 1 makes 0.86 the same number on either footing) and the neutron half does not. The convention statement is at `tagged.py`'s `TRITON`* |
| ⁶Li embedded-d b₁ scaling | `b1_li6_from_deuteron` (rank-2 0.921949 × 2/6, *2026-09-15*: 0.921947 before the i^L phase fix, which moved the model's own value 0.9219467 → 0.9219490 and no printed money_b1 number at all) | ☑ *2026-08-28 (plans/08 D9): the transfer is `TaggedModel(li6_alpha_channel()).tensor_dilution()` itself, pinned against it and against the closed form 1 − (9/10) P_D in `test_li6_b1_rank2_transfer_constant_is_pinned_to_the_model`. The 0.87 it replaces is the VECTOR dilution 1 − (3/2) P_D, the wrong rank for b₁, and is still reachable as `--transfer legacy`* |
| conservation & chain | HepMC3 → abconv → npsim | event-by-event 4-momentum/charge; 100-event smoke passes |

## 5.5 Risks

| risk | mitigation |
|---|---|
| spin-3/2 inclusive SF basis incomplete in the literature | rank ≤ 2 truncation + scenario shapes; co-author the formal note (turns a risk into a paper) |
| cluster-overlap tail dominates tagged acceptances (known from e+d control) | β-band always quoted; VMC overlaps as the scheduled fix; BeAGLE-vs-IA spread as the model systematic — ☑ *2026-09-15 (run 19b): the band is now quoted everywhere a β exists — `--beta-band` on the three spectator call sites, default off and bit-for-bit — and it is restated as one-sided **in β** rather than as a bracket: no β in a two-parameter Hulthén form reproduces BeAGLE's p_T tail, so the true scale is at or above 0.40. The tag's own direction is then channel-dependent, up for ⁶Li (×3.35–3.44 at the Yellow Report envelope) and down for ⁷Li (×1.01–1.04), and the sibling generator's VMC upgrade moves it in opposite directions at the two optics (0.0264 → 0.0348 at the YR envelope, 0.2551 → 0.2486 at the tagging optics), so it is a per-optics band and never one factor. §5.3 step 5.B carries the numbers* |
| α-spectator FSI beyond IA | quote at small |t′| (pole dominance); engage Cosyn/Sargsian (plans/04 #16); Mode W on BeAGLE brackets rescattering qualitatively |
| t* remnant treatment (d vs nn) too crude for double-tag studies | affects T1 tier only; gate double-tag claims on a ³He control (Friščić et al. PLB 823:136726 as template) — *2026-08-28: the gate is holding and nothing is due — no double-tag claim exists and the T1 remnant tier is unbuilt; the unpolarized ³He control has run, its polarized Friščić-template version waits on T1. What remains is the author's judgement of when a double-tag claim may be made* |
| tensor-observable RC unknown | RC hook + vector-case band (step 1.4 — the unpolarized ISR migration bound is measured, `polligen/radiative.py`, 2026-08-28); flag in every tensor plot |
| PYTHIA-tier integration friction (container-only) | tiers T0/T1 carry all Phase-1 physics; T2 only gates Phase-2 detector studies — *superseded (2026-08-28): plans/08 D4 closed it — PYTHIA 8.311 builds its own bindings against the analysis interpreter, so it was never a container problem, and 8 M events now sit in `evgen/samples/`* |
| relative luminosity at 10⁻⁴ unproven for tensor fills | bookkeeper makes it a knob; quote A_zz FOMs vs δ(rel-lumi) — feeds the machine requirement back to plans/04 #3 — ☑ *2026-08-28: knob in `bookkeeping.py` and on the CLI (`--rel-lumi-offset`, on by default in the 5R/7R runs), closed-form biases tested, quoted against the per-bin floors in Reports 0 and 2, and fed back to plans/04 #3. The two tolerances belong to different observables: 10⁻⁴ for the A_zz thirds estimator, 10⁻³ for the two-state cos 2φ ratio* |

## 5.6 Interfaces to the rest of the program

- **plans/02:** steps 1.3 (FOMs → pseudo-experiment grade), 1.4 (RC bands
  consumed as weights), 1.5 (BeAGLE: backgrounds/purity via Mode W; the
  cluster-IA cross-check of step 1.5.3 *is* Mode G's unpolarized limit).
- **plans/03:** step 2.1 consumes `io_hepmc` output; steps 2.3–2.4
  (φ-dilution, closure tests) run on Mode G samples; step 2.5's campaign
  request needs exactly this generator.
- **plans/04:** uses #1–#5 defaults as bookkeeper inputs; adds #14 (spin-3/2
  SF basis), #15 (VMC α+d/α+t overlaps incl. m-dependence), #16 (cluster
  FSI), #17 (HepMC3 ion-spin attribute convention).
  ☑ *2026-08-28: all four are in plans/04 with their tracking-table rows, and the
  #1–#5 defaults are consumed by the bookkeeper; the questions themselves stay
  open/external. Residual: `RunPlan.delta_p_over_p` (#5) is set only in tests, so
  no money plot yet carries a polarimetry-scale band.*
- **fastsim:** `polli_fastsim` is imported, not duplicated; anything the
  generator learns (acceptance-weighted dilutions, tail bands) flows back
  into the analytic FOM notebooks as parameterizations.
  ☐ *2026-08-28: the import discipline holds and is test-guarded; the flow-back
  does not — `fom.Scenario` and the three error functions take no dilution or
  acceptance argument, so `recopseudo`'s measured φ dilution and the β tail band
  never reach the analytic layer.*
  ☑ *2026-09-15: the flow-back now has its channel. `fom.Scenario` carries
  `dilution` (D) and `acceptance` (A), both defaulting to 1.0, both refused
  at ≤ 0, and all three δA paths divide by D·A — so `bin_summary`'s
  `dilution_phi` and `truth_reference`'s `dilution_beta` reach the analytic
  layer as numbers, which is how the import discipline stays one-way
  (fastsim still imports no `polligen`, now asserted by an AST scan in
  `fastsim/tests/test_fom_dilution.py`). δA(D = 0.5) = 2 δA(D = 1) exactly,
  and at the defaults every error is bit-for-bit what it was: the registered
  figure `phase_space_bins_6Li.png` regenerates to the same md5
  (b9a05445f26f2c1b7f89284133f1b2e4) across the change. What does NOT belong
  in these two fields is anything that removes events rather than attenuating
  the amplitude — that enters as 1/√ through the luminosity knobs, and the
  two laws are pinned apart in the tests.*
