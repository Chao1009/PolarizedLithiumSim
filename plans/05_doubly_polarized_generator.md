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
  per (x, y) that Mode W/G attach as optional weights. Tensor-observable RC
  remains uncharted (plans/04 #10) — the hook is there, the numbers await
  theory.

Language: numpy-vectorized Python, same style/test discipline as
`polli_fastsim` (it imports the existing SF backends and `farforward.py`);
4-vector-tier throughput ≥ 10⁴ ev/s is ample for Phase-1. PYTHIA tier runs
inside eic-shell (PYTHIA8 python bindings ship in the container) or via a
thin fifo/LHE hand-off — decided in step 5.D. C++ port only if an official
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
   spin-3/2 basis is itself a publishable theory note — engage Cloët/Cosyn),
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

### Step 5.B ☐ Tagged mode (2–3 weeks, the core novelty)
`tagged.py`: spin-correlated (e′, spectator) events for the four Li
channels + deuteron/³He controls. Validations in §5.4. **Deliverables:**
- **Money plot 4:** tagged tensor asymmetry A_zz^tag(p_s) for the α-tagged
  embedded deuteron in ⁶Li, folded with `farforward.py` acceptance and both
  optics — *first tagged spin observable for any A > 2* (extends the
  four-gap list of plans/00). Check the Cosyn–Weiss O(1) asymmetry at
  p_s ≈ 300 MeV/c survives the RP pT-tail acceptance that dominates ⁶Li.
- ⁷Li α-tag: tagged A∥ (polarized-EMC companion on the quasi-free triton)
  + the tagged-α angular-moment polarimetry curve vs P_zz.

### Step 5.C ☐ Mode W reweighter + reco-level closure (1–2 weeks)
`reweight.py` driving the official BeAGLE e+d sample (xrootd, already
streamed for the control study): inject A∥(g₁d) and A_zz(b₁d), extract with
the analysis estimators through the existing conversion chain, verify pulls.
This is the ECCE-style pedigree demonstration on the nucleus where BeAGLE
*is* right, and it exercises the identical machinery later pointed at e+Li
BeAGLE samples for purity studies (evaporation background stays BeAGLE's
job; plans/02 step 1.5.4).

### Step 5.D ☐ Final states + HepMC3 + chain smoke test (2 weeks)
Tier T1 (cluster-internal nucleon + partner spectators; t* remnant → d or
nn per the triton wave function — crude, flagged), `io_hepmc.py` (ASCII
HepMC3; status-4 beams so `abconv` accepts it; 10-digit fragment PDG;
spin labels as named attributes — no HepMC3 convention exists for ion spin
states, so define one and propose it to the ePIC MC group, new plans/04
#17), tier T2 PYTHIA attachment inside eic-shell. 100-event
abconv → npsim → EICrecon smoke per plans/03 step 2.1.4 (⁶Li can proxy the
existing d/⁴He beamline configs; ⁷Li needs the new field maps — already
plans/03 step 2.1.2).

### Step 5.E ☐ Physics production + write-up (2 weeks + ongoing)
Regenerate all money plots from generator pseudo-experiments (statistical
FOMs now include acceptance × estimator effects); tagged-FOM table
(efficiency × purity × dilution per channel); short generator note —
"first polarized eA event generator" is itself one of the publishable
firsts, aimed with the money plots at the INT program (Mar 22 –
Apr 2, 2027). Upgrades stay behind interfaces: VMC overlaps, digitized
CBT/TMT and b₁ theory curves (step 1.2), spin-3/2 rank-2 SFs, RC tables.

Total ≈ 7–9 focused weeks to 5.E; 5.A+5.C alone (≈ 3–4 weeks) already
upgrade every Phase-1 FOM to pseudo-experiment grade.

## 5.4 Validation matrix (gates, in order)

| gate | reference | pass criterion |
|---|---|---|
| ρ moments, all axes | analytic | exact (unit test) |
| master formula, vector/tensor sectors | `asymmetries.py` | bin-wise identity (toy + grid backends) |
| pseudo-experiment estimators | `fom.py` maps | δA agree within trial statistics; pulls unbiased |
| φ-modulation recovery | injected Δ scenarios | amplitude unbiased with uniform *and* holey φ acceptance |
| deuteron limit of tagged mode | Cosyn–Weiss arXiv:2603.23700 figures | shapes/magnitudes reproduced (digitized) |
| unpolarized spectator spectra | official BeAGLE e+d via `ed_control_analysis.py` | bulk agreement; tail differences documented as the model band |
| forward limit of tagged ⁷Li | P_p = 0.866, P_n = −0.037 | recovered within the D-state band |
| ⁶Li embedded-d b₁ scaling | `b1_li6_from_deuteron` (0.87) | recovered from tagged→inclusive integration |
| conservation & chain | HepMC3 → abconv → npsim | event-by-event 4-momentum/charge; 100-event smoke passes |

## 5.5 Risks

| risk | mitigation |
|---|---|
| spin-3/2 inclusive SF basis incomplete in the literature | rank ≤ 2 truncation + scenario shapes; co-author the formal note (turns a risk into a paper) |
| cluster-overlap tail dominates tagged acceptances (known from e+d control) | β-band always quoted; VMC overlaps as the scheduled fix; BeAGLE-vs-IA spread as the model systematic |
| α-spectator FSI beyond IA | quote at small |t′| (pole dominance); engage Cosyn/Sargsian (plans/04 #16); Mode W on BeAGLE brackets rescattering qualitatively |
| t* remnant treatment (d vs nn) too crude for double-tag studies | affects T1 tier only; gate double-tag claims on a ³He control (Friščić et al. PLB 823:136726 as template) |
| tensor-observable RC unknown | RC hook + vector-case band (step 1.4); flag in every tensor plot |
| PYTHIA-tier integration friction (container-only) | tiers T0/T1 carry all Phase-1 physics; T2 only gates Phase-2 detector studies |
| relative luminosity at 10⁻⁴ unproven for tensor fills | bookkeeper makes it a knob; quote A_zz FOMs vs δ(rel-lumi) — feeds the machine requirement back to plans/04 #3 |

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
- **fastsim:** `polli_fastsim` is imported, not duplicated; anything the
  generator learns (acceptance-weighted dilutions, tail bands) flows back
  into the analytic FOM notebooks as parameterizations.
