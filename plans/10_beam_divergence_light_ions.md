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

## 10.1 What the repository assumes

| where | value | provenance |
|---|---|---|
| `fastsim/polli_fastsim/farforward.py` | `HIGH_ACCEPTANCE` 0.20/(10·275) = **72.7 µrad**, `HIGH_DIVERGENCE` 0.45/(10·275) = **164 µrad** | a 275 GeV **proton**, back-derived from the documented p_T cuts |
| `evgen/polligen/reco.py` | `SIGMA_THETA_HA` = 72.7 µrad, `SIGMA_THETA_HD` = 0.41/(10·275) = **149 µrad** | same, and the high-divergence value **disagrees with `farforward.py` by 10%** — the module comment records the 0.41/0.45 rounding but the two are not reconciled |
| everywhere | **isotropic** (σ_x = σ_y), `aspect = 1.0` | the anisotropy parameter exists and is unused by default |
| everywhere | **energy-independent** | the same 72.7 µrad at 5 × 41, 10 × 100 and 18 × 275 |

Three assumptions, each of which the published tables contradict — and
§10.3 replaces all three with an estimate. `farforward.YR_PROTON_DIVERGENCE`,
`YR_GOLD_DIVERGENCE` and `yr_divergence_for()` now carry the tables and the
scaling, with `fastsim/tests/test_farforward.py` pinning them.

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
vertically**. Intrabeam scattering at fixed beam current goes as Z³/A² =
**0.75 for ⁶Li against 1 for a proton and 12.7 for gold, 17× lithium's**, and RHIC
deuterons showed no measurable IBS growth while gold grew 20–45%/h. Lithium
is a proton-class ion.

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

β\* = 13 m at 5 × 41 (from 0.90 m) is the smallest absolute ask and sits
well inside the LHC's demonstrated forward-physics optics — TOTEM/ALFA ran
β\* = 90 m and 2500 m against a nominal 0.55 m. Raising β\* also *shrinks*
the beam in the final-focus quadrupoles (β = 14.9 m at the first quad against
28.7 m today), so the IR aperture is not the constraint; matching and
chromaticity would be.

**Priced (2026-08-27, `evgen/scripts/tagging_optics.py`, Report 1 §6.1):**
with the envelope as planar pots see it (a rectangle 10σ_h × 10σ_v) only
the horizontal plane needs the de-squeeze, and the optimum of ε × L sits
at β*_x/β*_x,HA = 50 / 180 / 90 with the vertical plane at high acceptance
(L/L_HA = 1/7 / 1/13 / 1/9.5, ε = 0.32–0.42) rather than at the
circular-isotropic 14.5 / 57.7 / 28.6 above: 2.6×10⁶ / 3.0×10⁶ / 6.1×10⁶
tagged events per year at the 10 fb⁻¹/u placeholder — still 3–8× below
what IR-8's secondary focus (≈ 20%, our interpolation) would give at the
same luminosity.  Two assumptions stated there: the electron β* raised in
step, and a parallel-to-point far-forward transport for the de-squeezed
lattice (R₁₁σ* ≪ R₁₂σ_θ at the pots), which is the first thing to ask C-AD.

**But the detector must follow.** At β\* = 13 m the 10σ envelope is 0.58
mrad while the silicon aperture at 5 × 41 is ≈2–3 mrad, so the geometry
pins the acceptance and the whole β\* gain is wasted. The two levers are
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
luminosity fraction) are the two optics every script evaluates;
`coherent_optics_scan.py` panel (d) is per-configuration YR curves with
the tagging-optics points, `tagging_acceptance.py` tabulates the spectator
tags at both plus the legacy 73 µrad for reproduction (Report 3 Table 6).
`sigma_theta_for` and `hole_acceptance` live in `farforward` and
`polligen.reco` delegates to them, pinned by
`fastsim/tests/test_optics_20260828.py`.

### A4 — re-derive the near-beam study's gains on the band ☑
Re-derived 2026-08-28 (plans/09 §9.0, Report 4 §2).  The ×26 / ×569 gains
are withdrawn: at the Yellow Report optics the envelope is at or inside
the silicon at every configuration and a closer approach buys nothing
(coherent 7×10⁻⁸ / 6×10⁻²⁷ / 4×10⁻¹⁴ at the envelope).  At the tagging
optics the envelope is 0.33 / 0.17 / 0.12 mrad and the silicon at
1.0–2.0 mrad tags zero; pots that follow the envelope tag 0.42 / 0.32 /
0.33 of coherent recoils and 0.35 / 0.27 / 0.28 of α spectators at
1/7–1/13 of the luminosity, with four clean |t| bins per configuration
through the chain (`nearbeam_reach_gain.py`).  The layer and the optics
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

* **plans/04 #11 and #20** — the optics half of both is promoted here and
  given a number.
* **plans/09** — the near-beam gains are ratios and survive; the absolute
  tagged fractions do not.
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
