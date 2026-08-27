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

Three assumptions, each of which the published tables contradict.

## 10.2 What the Yellow Report actually says

YR Table 8.9 (`refs/2103.05419_part2.pdf`) gives RMS Δθ h/v for **e+Au**,
the only heavy-ion bracket published, and it is neither
energy-independent nor isotropic nor close to the proton number:

| species / energy | strong hadron cooling | stochastic cooling |
|---|---|---|
| Au 110 GeV/u | **218 / 379 µrad** | 77 / 380 µrad |
| Au 41 GeV/u | **275 / 377 µrad** | 174 / 302 µrad |

Against the repo's 72.7 / 164 µrad that is **2.5–3.8× horizontally and
2.3–2.3× vertically** under strong hadron cooling, and the *spread
between the two cooling options in the same table* is itself a factor 2.8
horizontally.

Independent corroboration that σ_θ cannot be energy-independent: ePIC's
own per-energy 10σ Roman-Pot retractions in
`compact/fields/beamline_*.xml` are **2.67 / 3.11 / 5.36 cm** at
18 × 275 / 10 × 100 / 5 × 41. A single σ_θ cannot produce a retraction
that grows by a factor 2 toward *lower* energy.

## 10.3 What it is worth

`reco.rp_hole_acceptance` at B = 50 GeV⁻², 10σ envelope only (no
geometric aperture), coherent intact ⁶Li:

| σ_θ assumption | 5 × 41 | 10 × 100 | 18 × 275 |
|---|---|---|---|
| repo (72.7 / 164 µrad, energy-independent) | 3.99×10⁻¹ | 2.92×10⁻² | 2.0×10⁻⁹ |
| YR e+Au 110 GeV/u, stochastic cooling | 3.44×10⁻¹ | 2.09×10⁻² | 2.1×10⁻¹⁰ |
| YR e+Au 110 GeV/u, strong hadron cooling | 7.34×10⁻³ | 6.2×10⁻¹¹ | ~10⁻⁷² |
| YR e+Au 41 GeV/u, strong hadron cooling | 7.22×10⁻⁴ | 1.6×10⁻¹⁶ | ~10⁻¹¹⁴ |

**A factor 553 at 5 × 41 and fourteen orders of magnitude at 10 × 100,
between two rows of the same published table.** The coherent programme is
either healthy at low energy or does not exist, and which one is true is
decided by a number nobody in this repository has ever sourced for
lithium.

**Caveat, stated plainly:** Au is not Li. Divergence is
√(ε_N/(βγ β*)), so it depends on the species, the source and the cooling
scheme, and the ⁶Li values could fall anywhere in or outside this
bracket. These rows are the available *bracket*, not a prediction. That
is precisely the problem.

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
**The blocking input for the entire far-forward programme.** Ask for it
per beam-energy configuration and per cooling scenario. This subsumes the
optics half of plans/04 #11 and #20 and is the same conversation as
plans/09 D2.

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
  intact-⁶Li measurement is feasible at low energy *if* the lithium beam
  divergence resembles the proton-derived assumption, and not otherwise —
  and settling that is one question to C-AD, not a simulation task.
