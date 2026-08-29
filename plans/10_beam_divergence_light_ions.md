# Plan 10 — The beam divergence the whole far-forward programme rests on (2026-08-27)

**Why this is its own document.** Every far-forward acceptance number in
this repository — the coherent intact-⁶Li tag, the α tag, the near-beam
study of plans/09, the money plots' 10σ envelope — is an exponential in
the *square* of the beam angular divergence:

    tagged fraction ~ exp(−B (10 σ_θ · A · p_u)²)

and `σ_θ` is currently **one energy-independent, proton-derived number**.
Nothing else in the chain has that much leverage on the headline result,
and nothing else is fixed by an input we have never checked for a
lithium beam.

---

## 10.1 What the repository assumed before 2026-08-27

> **Superseded 2026-08-28 (A1, A1b, A2 and A3 of §10.4).** All three
> assumptions below have been replaced in the code.  The section is kept
> because every figure published before 2026-08-27 rests on them, and
> because the legacy constants survive, labelled, for those
> reproductions.

| where | value | provenance |
|---|---|---|
| `fastsim/polli_fastsim/farforward.py` | `HIGH_ACCEPTANCE` 0.20/(10·275) = **72.7 µrad**, `HIGH_DIVERGENCE` 0.45/(10·275) = **164 µrad** | a 275 GeV **proton**, back-derived from the documented p_T cuts |
| `evgen/polligen/reco.py` | `SIGMA_THETA_HA` = 72.7 µrad, `SIGMA_THETA_HD` = **164 µrad**, the `farforward` constant itself | the two packages disagreed by 10% until 2026-08-28: `reco.py` derived the high-divergence value from the 0.41 GeV end of the band, 0.41/(10·275) = 149 µrad, against `farforward.py`'s rounded 0.45. A1b made `reco.SIGMA_THETA_HD` an alias, pinned by `evgen/tests/test_review_20260828.py` |
| everywhere | **isotropic** (σ_x = σ_y), `aspect = 1.0` | the anisotropy parameter existed and was unused by default; A2 made `sigma_theta_for` return the (h, v) pair and threaded it to the call sites |
| everywhere | **energy-independent** | the same 72.7 µrad at 5 × 41, 10 × 100 and 18 × 275; A3 replaced it with `farforward.yr_optics(config)` and `tagging_optics(config)` |

Three assumptions, each of which the published tables contradict — and
§10.3 replaces all three with an estimate. `farforward.YR_PROTON_DIVERGENCE`
and `YR_GOLD_DIVERGENCE` carry the tables and `farforward.sigma_theta_for()`
the species step, with `fastsim/tests/test_farforward.py` pinning them.  The
blanket √(A/Z) that the first correction applied through
`yr_divergence_for()` and `LIGHT_ION_DIVERGENCE_FACTOR` was **retired on
2026-08-28** together with those two names: a γ-matched ion pays no species
penalty at all, and only the rigidity-capped top configuration pays √2
(§10.3).

## 10.2 What the Yellow Report actually says

**Table 10.1** (e+p) and **10.2** (e+Au) of `refs/2103.05419_part2.pdf` give
RMS divergence h/v and RMS Δp/p per configuration and per optics. The hadron
rows, verbatim:

| protons | high divergence h/v | high acceptance h/v | Δp/p |
|---|---|---|---|
| 275 GeV, 290 bunches | 150 / 150 µrad | 65 / 65 | 6.8×10⁻⁴ |
| 275 GeV, 1160 bunches | 119 / 119 | 65 / 65 | 6.8×10⁻⁴ |
| 100 GeV | 220 / 220 | 180 / 180 | 9.7×10⁻⁴ |
| 100 GeV (vs e 5) | 206 / 206 | 180 / 180 | 9.7×10⁻⁴ |
| 41 GeV | 220 / 380 | 220 / 380 | 10.3×10⁻⁴ |

| gold | strong hadron cooling h/v | Δp/p | stochastic h/v | Δp/p |
|---|---|---|---|---|
| 110 GeV/u | 218 / 379 | 6.2×10⁻⁴ | 77 / 380 | 10×10⁻⁴ |
| 41 GeV/u | 275 / 377 | 10×10⁻⁴ | 174 / 302 | 13×10⁻⁴ |

Three things the repo's single number cannot represent: the divergence is
**energy-dependent** (65 → 180 → 220 µrad in high acceptance), it is
**anisotropic at 41 GeV** (220/380), and at 41 GeV there is **no separate
high-acceptance option at all** — the two rows are identical.

*Consistency check:* YR Table 11.48's angular-divergence column (40 / 22 / 14
MeV of Δp_T) is Table 10.1's divergence times the proton momentum, to within
6–9% at 100 and 41 GeV. The two tables agree, so they are being read right.

## 10.3 The estimate for lithium

**First, the energies were wrong, and that matters more than the divergence.**
The HSR and the ESR must have equal revolution period and the electrons are
ultrarelativistic, so the hadron **γ** is fixed by the ring circumference and
the magnets supply whatever rigidity that γ demands, up to the 917.3 T·m cap.
Ions are therefore **γ-matched**, not rigidity-scaled.

**The mechanism, verbatim** (EPIOS pp. 12–13, added 2026-08-27): *"their
revolution frequencies have to be equal. This synchronization is accomplished
by applying a radial shift of up to ±20 mm in the arcs, which facilitates a
range of the Lorentz factor of **118 < γ < 293**. To allow for even lower ion
energies, a 'Blue' arc between IR12 and IR2 will be utilized as a bypass. The
average radius of this arc is about 90 cm smaller than that of the
corresponding 'Yellow' arc, which reduces the circumference of the HSR by
roughly 90 cm. The resulting circumference then corresponds to a Lorentz
factor of **γ = 43.5**."*

Those two numbers are exactly the species menu at *top* energies — gold at
110 GeV/u is γ = 118.1 and a 275 GeV proton is γ = 293.1 — so the window
describes what the ion programme needs rather than a hard reachability bound,
and the shift is quoted as "up to" ±20 mm.

> **A conflict, recorded rather than resolved.** YR Table 10.1 runs 100 GeV
> protons, i.e. γ = 106.6, which falls in *neither* stated window. The 41 GeV
> (γ = 43.7) and 275 GeV (γ = 293.1) anchors are both inside them, and YR
> Table 10.2's gold at 41 GeV/u (γ = 44.0) sits in the bypass window
> alongside the 41 GeV proton — which is the independent confirmation that
> ions are γ-matched. `beams.py` stays anchored on the Yellow Report's three
> configurations because every beam-parameter table is indexed by them, and
> `tools/consistency_check.py` flags the 100 GeV point so the conflict stays
> visible. Resolving it is a question for C-AD.

*The decisive check:* YR Table 10.2 lists **gold at 41 GeV/u** (γ = 44.02)
against the 41 GeV proton's γ = 43.70 — equal to 0.7%, exactly the u vs m_p
mass difference. Rigidity scaling would have put gold at 41 × 79/197 = 16.4
GeV/u. It does not. The same rule reproduces the published ³He menu
("41, and 100–183 GeV/nucleon") exactly.

So the ⁶Li menu is **41, and 99–138 GeV/u**, with nothing in between, and
`beams.default_configs("6Li")` now returns **40.8 / 99.5 / 137.5** where it
returned 20.5 / 50 / 137.5 before 2026-08-27. The implementation
self-checks: run the same machinery on a *proton* and it returns 41 / 100 /
275 exactly.

**Second, the species step then applies only at the top.** A γ-matched ion
has the *same* βγ as its proton, hence the same geometric emittance at equal
ε_N, hence the **same divergence — no √2 penalty**. Only at the top
configuration is ⁶Li rigidity-capped (γ = 147.6 against the proton's 293),
and there it does pick up √2.

| ⁶Li | p_ion | high acceptance h/v | high divergence h/v | Δp/p |
|---|---|---|---|---|
| **40.8 GeV/u** (5 × 41) | 245 GeV | **220 / 380 µrad** | 220 / 380 | 10.3×10⁻⁴ |
| **99.5 GeV/u** (10 × 100) | 597 GeV | **180 / 180** | 220 / 220 | 9.7×10⁻⁴ |
| **137.5 GeV/u** (18 × 275) | 825 GeV | **92 / 92** | 212 / 212 | 6.8×10⁻⁴ |

**Is equal ε_N defensible?** Gold is the published test. Scaling the 275 GeV
proton to Au at 110 GeV/u by 1/√(βγ) alone predicts 236 µrad against an
observed 218 (h) and 379 (v) — ε_N(Au)/ε_N(p) = **0.85 horizontally, 2.6
vertically**. Intrabeam scattering separates gold from lithium by far more
than that, under either normalisation of the same law. Per particle the
growth rate goes as N Z⁴/A², so at equal bunch intensity gold is **446×**
lithium (Z⁴/A² = 2.25 for ⁶Li, 1 for a proton, 1004 for gold); at fixed
beam *current*, where N ∝ 1/Z because a fixed current is fewer ions when
each carries more charge, it goes as Z³/A² and gold is **17×** lithium
(0.75, 1, 12.7). `farforward.py` and `fastsim/tests/test_farforward.py`
carry both. RHIC deuterons showed no measurable IBS growth while gold
grew 20–45%/h. Lithium is a proton-class ion on either count.

**Δp/p is the easy one.** Protons give 6.8–10.3×10⁻⁴ and gold 6.2–13×10⁻⁴
across every energy and both cooling schemes — set by the RF bucket rather
than the species, so ⁶Li is **7–10×10⁻⁴** with good confidence.

### What it costs

Coherent tagged fraction, B = 50 GeV⁻², 10σ:

| config | as published | energy fix alone | + divergence | total cost |
|---|---|---|---|---|
| 5 × 41 | 6.70×10⁻¹ | 2.05×10⁻¹ | **5.0×10⁻⁷** | ×1.3×10⁶ |
| 10 × 100 | 9.25×10⁻² | 8.1×10⁻⁵ | **8.4×10⁻²⁶** | ×1.1×10²⁴ |
| 18 × 275 | 1.52×10⁻⁸ | 1.52×10⁻⁸ | **3.7×10⁻¹³** | ×4.1×10⁴ |

**The coherent intact-⁶Li channel does not survive at any configuration.**
Note the energy correction alone costs a factor 3.3 at 5 × 41 and 1100 at
10 × 100, before the divergence is touched.

### What recovers it

Acceptance goes as exp(−C/β\*) while luminosity goes as 1/β\*, so the figure
of merit L × acceptance is maximised where the 10σ cut sits at t = 1/B, i.e.
p_T = 1/√B = 0.141 GeV. `reco.sigma_theta_tagging` returns that working
point, and by construction the acceptance there is **1/e at every
configuration** — an invariance that is itself the derivation:

| config | tagging σ_θ | β\* factor over high acceptance | acceptance |
|---|---|---|---|
| 5 × 41 | 57.8 µrad | **×14.5** | 0.368 |
| 10 × 100 | 23.7 µrad | ×57.7 | 0.368 |
| 18 × 275 | 17.1 µrad | ×28.6 | 0.368 |

That table is the circular-isotropic derivation, and run 12 replaced it
with the horizontal-only de-squeeze priced immediately below; the absolute
ask that follows from the corrected optimum is **β*_x ≈ 42 m at 5 × 41**
(0.90 m × 46.5), with ≈ 100 m and ≈ 71 m at the other two configurations
(Report 0 §4.2). Forty-two metres is the smallest of the three and sits
well inside the LHC's demonstrated forward-physics optics — TOTEM/ALFA ran
β\* = 90 m and 2500 m against a nominal 0.55 m. It does not, however, relax
the IR aperture, as the superseded ×14.5 would have: that optimum shrank the
beam in the final-focus quadrupoles (β ≈ β\* + L²/β\* ≈ 15 m at the first
quad against 29 m today, for L ≈ 5 m), whereas a 42 m horizontal de-squeeze
leaves β_x ≈ 43 m there, larger than today. Aperture, matching and
chromaticity are all questions for C-AD.

**Priced (2026-08-27, `evgen/scripts/tagging_optics.py`, Report 1 §6.1):**
with the envelope as planar pots see it (a rectangle 10σ_h × 10σ_v) only
the horizontal plane needs the de-squeeze, and the optimum of ε × L sits
at β*_x/β*_x,HA = 46.5 / 164.1 / 89.3 with the vertical plane at high
acceptance (L/L_HA = 1/6.8 / 1/12.8 / 1/9.5, ε = 0.37 / 0.25 / 0.33,
`farforward.tagging_optics_point`) rather than at the
circular-isotropic 14.5 / 57.7 / 28.6 above: 2.4×10⁶ / 2.4×10⁶ / 6.2×10⁶
tagged events per year at the 10 fb⁻¹/u placeholder — still 4–10× below
what IR-8's secondary focus (≈ 20%, our interpolation) would give at the
same luminosity.  Two assumptions stated there: the electron β* raised in
step, and a parallel-to-point far-forward transport for the de-squeezed
lattice (R₁₁σ* ≪ R₁₂σ_θ at the pots), which is the first thing to ask C-AD.

**But the detector must follow.** At the tagging optics the horizontal
10σ envelope is 0.36 / 0.19 / 0.12 mrad while the measured silicon edge is
2.50 / 1.51 / 0.53 mrad (`tools/fullsim`, re-measured 2026-08-28 in the
current ePIC geometry), 6.9 / 7.9 / 4.5 times outside it, so
the geometry pins the acceptance and the whole β\* gain is wasted unless the
pots follow the envelope in. The two levers are
strictly multiplicative, and the second one is the near-beam granularity
question of plans/09 — the module-quantised insertion and the x-moving layer
ePIC is already designing.

**There is no high-acceptance optics at 41 GeV today.** CDR Tables 3.3 and
3.4 are identical in that column — the only energy where they are. Whether
that is a physics limit or an unstudied case is the question of D1 below,
and it decides whether the recovery exists.

## 10.4 Ordered work

### A0 — the energies ☑
`beams.default_configs` is γ-matched and rigidity-capped, `beams.NUCLEUS_MASS`
carries the physical masses per nucleon, and `fastsim/tests/test_beams.py`
pins the menu, the proton self-check and a guard that sweeps the source tree
for the stale 20.5 / 50 GeV/u.

### A1 — per-configuration divergence ☑
`reco.sigma_theta_for(config, optics)` returns the YR Table 10.1 divergence
with the species step applied only where rigidity binds;
`reco.sigma_theta_tagging(config)` returns the L × acceptance optimum.
`farforward.YR_PROTON_DIVERGENCE` / `YR_GOLD_DIVERGENCE` hold the tables.
The legacy `SIGMA_THETA_HA` / `_HD` remain, labelled, because every number
published before 2026-08-27 used them.

### A1b — reconcile the two legacy values ☑
`farforward.py` said 164 µrad, `reco.py` 149 µrad, for the same
"high-divergence" optics. `reco.SIGMA_THETA_HD` is now the `farforward`
constant, pinned by a test (2026-08-28); both are legacy placeholders that
only the pre-2026-08-27 figures used.

### A2 — make the call sites use the anisotropy ☑
`sigma_theta_for` returns (h, v) and the low configuration is genuinely
anisotropic (220/380). `recopseudo.CoherentResponse` accepts the pair,
`reco.tagging_optics_point` returns it, and `money_cos2phi_coherent_reco.py
--optics tagging` / `high-acceptance` thread it through (2026-08-28); the
legacy scalar remains the script's default for reproduction only.

### A3 — publish the sensitivity band, not a single number ☑
Done 2026-08-28 the other way round: instead of a band, every far-forward
figure now states its optics per configuration.  `farforward.yr_optics`
(the Yellow Report high-acceptance / high-divergence rows, a rectangular
10(σ_h, σ_v) envelope applied to each fragment's azimuth) and
`farforward.tagging_optics` (the Report 1 §6.1 optimum, with its
luminosity fraction) are the two optics the per-configuration figures
evaluate: `coherent_optics_scan.py` panel (d) is per-configuration YR
curves with the tagging-optics points, `tagging_acceptance.py` tabulates
the spectator tags at both plus the legacy 73 µrad for reproduction
(Report 3 Table 6), and `money_cos2phi_coherent_reco.py --optics`,
`tagging_optics.py`, `eic_beam_figures.py`, `reco_chain_figures.py` and
the five `nearbeam_*.py` take the per-configuration values.
`sigma_theta_for` and `hole_acceptance` live in `farforward` and
`polligen.reco` delegates to them, pinned by
`fastsim/tests/test_optics_20260828.py`.

**The two tagged scripts joined them the same evening** (plans/09 B2,
B3), which closes this item for the tagged observables.
`money_tagged_azz.py` and `tagged_polarimetry_7li.py` take `--config
{0,1,2}` and `--optics {menu,legacy,high-acceptance,high-divergence,
tagging}`, default `menu` = the configuration's Yellow Report
high-acceptance optics plus the tagging optics with its luminosity
fraction; `legacy` reproduces the retired 73 / 164 µrad pair and is the
only place it survives in a published figure.  Threading the flag was
not the whole fix: `polligen.tagged` never exposed the spectator's lab
azimuth, so the rectangular envelope degenerated to a circle at n σ_h
wherever these scripts had used it — ×1.7 too generous at the tagging
optics — and `TaggedSampler`'s default optics is now
`farforward.yr_optics(beam_config)` rather than a module constant that
cannot know the beam.

**Two scripts are still legacy, for reproduction, and say so.**
`money_cos2phi_coherent.py` and `phase_space_bins.py` price
the coherent tag with `Optics.pt_cut_near_beam`, the 0.20 and 0.45 GeV
cuts as a **275 GeV proton** sees them, which for a ⁶Li at 40.8–137.5
GeV/u is between 5σ and 2σ rather than 10σ; their figure annotations, which
said only "near-beam envelope", now carry "(legacy)".

### A4 — re-derive the near-beam study's gains on the band ☑
Re-derived 2026-08-28 (plans/09 §9.0, Report 4 §2).  The ×26 / ×569 gains
are withdrawn: at the Yellow Report optics, and on the aperture
re-measured in the current ePIC geometry on 2026-08-28, the silicon binds
at 5 × 41 (2.50 mrad horizontally against a 2.20 mrad envelope) and the
machine at the other two (1.51 against 1.80, and 0.53 against 0.92), so a
nearer sensor buys ×77 at the low configuration and exactly nothing at
the other two, none of it worth having (coherent
7.2×10⁻⁸ / 6.2×10⁻²⁷ / 7.1×10⁻¹⁴ at the envelope —
`nearbeam_aperture_scan.py`'s column, i.e. the envelope horizontally and
the larger of silicon and envelope vertically, which is how the pots
move; the pure envelope in both planes would read
7.2×10⁻⁸ / 1.2×10⁻²⁶ / 7.8×10⁻¹⁴, and the scan's silicon column,
9.4×10⁻¹⁰ / 2.0×10⁻¹⁹ / 1.2×10⁻⁵, is the aperture alone and therefore
inside the beam at the two upper configurations, where it is not a
reachable number).  At the tagging
optics the envelope is 0.36 / 0.19 / 0.12 mrad and the silicon at
0.53–2.50 mrad tags zero; pots that follow the envelope tag 0.37 / 0.25 /
0.33 of coherent recoils and 0.31 / 0.22 / 0.29 of α spectators at
1/6.8–1/12.8 of the luminosity, with all seven bins of the published |t|
window per configuration through the chain (`nearbeam_reach_gain.py`).  The layer and the optics
are multiplicative levers.

### D1 — σ_θ (h and v) and Δp/p at the IP for a ⁶Li / ⁷Li fill · **C-AD / ePIC FF WG**
**Provisionally answered in §10.3** by scaling the Yellow Report's own
proton tables, with the equal-emittance assumption calibrated against gold.
What remains to be asked is narrower and sharper than "what is σ_θ":
1. **Is there, or could there be, a light-ion high-acceptance optics** — a
   β* choice made for tagging rather than luminosity? This is the one lever
   that recovers the coherent channel, and §10.3 shows nothing else does.
2. Does the polarized ⁶Li source deliver a normalised emittance comparable
   to the proton beam's, or is the ring unable to cool it to equilibrium at
   the low intensity a polarized source provides? (§10.3 assumes it can;
   the gold calibration supports it, the source does not confirm it.)
3. Confirm Δp/p ≈ 7–10×10⁻⁴ for a light ion — expected to be the least
   species-sensitive number in the problem.

### D2 — which cooling scenario is the baseline for ion running · **C-AD**
The two rows of YR Table 10.2 differ by a factor 2.8 in horizontal gold
divergence. **Partly answered 2026-08-27:** EPIOS p. 13 states that the
baseline is conventional electron cooling **at injection**, and that *"at a
later stage, cooling at collision energies could be added"* on a Coherent
electron Cooling scheme, worth a factor two in average luminosity. So the
baseline has **no cooling at store** — the emittance grows under IBS through
the fill with nothing fighting it, and any argument that a low-intensity
lithium beam reaches a colder *equilibrium* has no equilibrium to reach.
What remains open is the injection-cooler's reach for a light ion: the Low
Energy Cooler is matched to the proton injection γ, and an A/Z = 2 ion at the
same injection *rigidity* sits at roughly half of it.

### D3 — β* for light-ion running · **C-AD**
The other half of σ_θ = √(ε_N/(βγ β*)), and the one that the
high-acceptance / high-divergence choice actually names.

---

## 10.5 What this changes elsewhere

* **plans/04 #20** — the optics half is promoted here and given a number.
  **#11** is not: its transfer-matrix half — what the far-forward
  transport looks like at the pot planes, and what it becomes under a
  de-squeezed lattice — stays with the ePIC FF WG, and nothing in this
  file answers it.
* **plans/09** — neither the absolute tagged fractions nor the gains
  survive.  The ×26 / ×569 are withdrawn with them (A4): the recoil must
  clear the larger of the envelope and the aperture per plane, and on the
  re-measured aperture the envelope is the larger at 10 × 100 (1.80 mrad
  against a 1.51 mrad silicon edge) and at 18 × 275 (0.92 against 0.53),
  so approaching closer buys exactly nothing there; at 5 × 41 it buys the
  ×77 the aperture scan prints (7.2×10⁻⁸ against 9.4×10⁻¹⁰) on a tagged
  fraction of 9×10⁻¹⁰, which is nothing that matters.  The
  near-beam layer pays only *under a tagging optics*, where it is the
  difference between no coherent tag at all and 25–37%.
* **reports/1 (cos 2φ projections) and reports/3 (detector study)** —
  both quote acceptances conditioned on the single proton-derived σ_θ.
  Neither should be circulated further without the band of §10.3 or an
  explicit statement of which row it assumes.
* **The honest summary for a collaboration meeting:** the coherent
  intact-⁶Li tag was projected on a proton-derived divergence that the
  Yellow Report's own beam tables do not support for a light ion. Scaling
  those tables properly costs three to four orders of magnitude in the
  tagged fraction, and **no sensor, aperture or reconstruction change
  recovers it** — the only lever is a light-ion optics with β* chosen for
  tagging rather than luminosity. That is the question to put to C-AD, and
  it is a sharper one than "what is σ_θ".
