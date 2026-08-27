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

**The species step is kinematic, and it is the whole story.** At a given
machine configuration the lattice is set by *rigidity*, so β* is common to
every species and a fully stripped A/Z = 2 ion sits at **half** the proton's
βγ. Since σ_θ = √(ε_N/(βγ β*)),

    sigma_theta(6Li) = sqrt(A/Z) x sigma_theta(proton, same configuration)
                     = 1.409 x  ... at equal normalised emittance

**Is equal ε_N defensible?** Gold is the published test. Scaling the 275 GeV
proton to Au at 110 GeV/u by 1/√(βγ) alone predicts 236 µrad against an
observed 218 (h) and 379 (v) — i.e. ε_N(Au)/ε_N(p) = **0.85 horizontally and
2.6 vertically**. Gold's intrabeam scattering, which goes as N·Z⁴/A², is
**~450× lithium's** (Z⁴/A² = 1004 for Au, 2.25 for ⁶Li, 1 for p), and it
still costs at most a factor 2.6. **Lithium sits far closer to a proton than
to gold on every axis that matters, so the lithium divergence is set by
kinematics, not by IBS.**

| ⁶Li | high acceptance h/v | high divergence h/v | Δp/p | repo uses |
|---|---|---|---|---|
| 20.5 GeV/u (5 × 41) | **310 / 535 µrad** | 310 / 535 | 10.3×10⁻⁴ | 72.7 / 163.6 |
| 50 GeV/u (10 × 100) | **254 / 254** | 290 / 290 | 9.7×10⁻⁴ | 72.7 / 163.6 |
| 137.5 GeV/u (18 × 275) | **92 / 92** | 168 / 168 | 6.8×10⁻⁴ | 72.7 / 163.6 |

**The repo's 72.7 µrad is low by ×1.3 at top energy, ×3.5 at mid, and ×4.3
at the low-energy configuration** — and low energy is exactly where the
coherent programme was said to survive.

**Δp/p is the easy one.** Protons give 6.8–10.3×10⁻⁴ and gold 6.2–13×10⁻⁴
across every energy and both cooling schemes. It is set by the RF bucket and
the cooling equilibrium rather than by species, so **⁶Li ≈ 7–10×10⁻⁴** with
good confidence. Note this is *larger* than the α/⁶Li rigidity difference of
1.87×10⁻³ divided by 3 — i.e. the beam's own momentum spread is within a
factor ~2 of the rigidity separation any dispersive Z-ID would need to
resolve (plans/09 §9.2).

### What it costs

Coherent tagged fraction, B = 50 GeV⁻², 10σ envelope, high acceptance:

| config | repo (72.7 / 163.6) | estimate | cost |
|---|---|---|---|
| 5 × 41 | 3.99×10⁻¹ | **1.4×10⁻⁴** | ×2900 |
| 10 × 100 | 2.92×10⁻² | **5.5×10⁻¹⁴** | ×5×10¹¹ |
| 18 × 275 | 2.0×10⁻⁹ | 8.3×10⁻¹⁴ | ×2.4×10⁴ |

**On this estimate the coherent intact-⁶Li channel does not survive at any
configuration**, best case 1.4×10⁻⁴ at 5 × 41 against the 0.40 assumed. The
near-beam study's ×26 and ×569 gains (plans/09) are ratios and still hold —
but they are gains on a base three to four orders smaller than quoted.

### What would move it

The estimate assumes ε_N(⁶Li) = ε_N(proton) and the *published proton*
β* per configuration. Two things could recover acceptance, and both are
machine choices rather than simulation ones:

* **A dedicated light-ion high-acceptance optics.** The high-acceptance
  configuration is a deliberate β* choice that costs luminosity (the YR is
  explicit: σ ∝ 1/√β*, β* ∝ 1/L). No light-ion optics is published. Raising
  β* by 4 would halve σ_θ and buy back orders of magnitude in the tag.
* **Running lithium at the lowest rigidity available.** The kinematic √(A/Z)
  penalty is unavoidable, but the acceptance depends on σ_θ·A·p_u, so the
  low-energy configuration remains the right home even after this correction.

## 10.4 Ordered work

### A1 — reconcile the two internal values ☐
`farforward.py` says 164 µrad, `reco.py` says 149 µrad, for the same
"high-divergence" optics. One of them is wrong. Fix, and add a test that
pins them together.

### A2 — make σ_θ per-optics and anisotropic ☐
`Optics` already carries a single scalar. It needs (σ_x, σ_y) per beam
configuration, and every call site needs to stop assuming isotropy. The
`aspect` parameter in `reco.rp_measure` and `recopseudo.CoherentResponse`
already exists for this and defaults to 1.0.

### A3 — publish the sensitivity band, not a single number ☐
Every acceptance and reach figure in the programme should carry the
divergence band of §10.3, or state explicitly which row it assumes.
`coherent_optics_scan.py` already scans a σ_θ list; extend it to the
published ion rows and make the band the headline rather than a footnote.

### A4 — re-derive the near-beam study's gains on the band ☐
plans/09 §9.0 quotes ×26 / ×569 gains from closing the aperture to 10σ.
Those are ratios at fixed σ_θ and survive a change of σ_θ better than the
absolute acceptances do — but the *absolute* numbers in that table move
with §10.3 and must be re-quoted.

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
The two rows of YR Table 8.9 differ by ×47 in coherent acceptance at
110 GeV/u. Without knowing which is planned, the programme cannot quote a
projection at all.

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
