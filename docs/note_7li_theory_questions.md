# The ⁷Li programme and what blocks it: critical theoretical questions, asks, and interim assumptions

*Working note, 2026-08-29 — polarized ⁶,⁷Li @ EIC simulation program (companions:
`plans/09` §B3a, the coherent ⁷Li scoping memo; `plans/04` #9/#10/#14/#15/#16/#18; `plans/05`
§5.2, the kernel these inputs feed), for circulation ahead of the INT program "Towards
Realizing the Program with Polarized Ion Beams at EIC", 22 March – 2 April 2027. Every
number is reproduced by one of two commands in `docs/reproduction_manual.md`: the
scoping memo’s block in §4.6, and §4.6b, "The ⁷Li theory-questions note".*

---

⁷Li is not a variant of the ⁶Li program but its complement, and it buys two things.
The first is deformation: Q(⁷Li) = −4.00(3) fm² against Q(⁶Li) = −0.0806(6) fm²
(Stone, INDC(NDS)-0833 2021), a relative dial |Q|/R_ch² of 0.6697 against 0.01202 at
R_ch = 2.444 and 2.589 fm — a factor 56 between isotopes measured with one estimator
in one window, so that ⁶Li is the null test and ⁷Li the positive control of the same
coherent cos 2φ measurement. The second is the far-forward tag: the α spectator sits
at rigidity R = 0.8557 — clear of the R ≈ 1 at which ⁶Li's α is beam-blind — and lands
mid-window in the Roman Pots at 10 × 99.5 and 18 × 117.9; at 5 × 40.8 the measured
per-band insertion holds the 32–48 mm silicon off to |y| ≥ 18 mm, so the same α misses
station 1 entirely and reaches silicon at station 2 only, in 57% of events. The
over-rigid triton at R = 1.2897 reaches Roman-Pot silicon in 60 of 60 events, opening
an α + t double tag (Report 3 Table 6 and its footnote; `plans/09` B1).

What ⁷Li does not come with is theory. The rank-2 machinery is already isotope-generic
— `xsec.InclusiveKernel._tensor_moments` returns Q_NN = [3m² − J(J+1)]/3 = (+1, −1,
−1, +1) over m = +3/2 … −3/2 with c_eff = 3 Q_NN, so a fill averages to ⟨c_eff⟩ = 3T
exactly, populations (0.4, 0.1, 0.1, 0.4) giving T = 0.600 and ⟨c_eff⟩ = 1.800 — but
the slots it multiplies (`b1_32_func`, `b2_32_func`, `delta_32_func`) default to
`None`, and a kernel with no structure function in it produces no number. A question
appears below only when a *published* number cannot be produced without its answer,
ordered by how hard it binds; items 1–6 are ⁷Li's own, item 7 the two ⁶Li asks that go
to the same people in one message.

## 1. The inclusive rank-2 (and rank-3) basis for J = 3/2, and its conventions

**Needed:** the decomposition of the inclusive unpolarized-lepton cross section on a
spin-3/2 target — the b₁-analogue set multiplying the rank-2 alignment in the rate,
the double-helicity-flip analogue of Δ multiplying cos 2φ — normalized against a named
alignment tensor over 0.01 ≲ x ≲ 0.5, 1 ≲ Q² ≲ 20 GeV², and two things usually left
implicit: whether rank 3 contributes to cos 2φ at leading twist, and in which frame
the alignment axis is defined (the spin-1 literature puts it along **q** to kill T_LT
and T_TT; we inherit that untested). **Blocks:** every ⁷Li inclusive figure — with
`delta_32_func = None` the cos 2φ amplitude vanishes identically, so ⁶Li's gluonometry
reach has no ⁷Li counterpart and the "⁷Li eligible" entry of Report 0 Table 3 is an
eligibility of the kernel, not of the observable. **Exists:** for spin 1, Cosyn–Roldan
Tomei–Sosa–Zec, EPJ A 61 (2025) 83 (`refs/2410.12764v1.pdf`), whose Eq. (9) is the
single-scalar alignment form our kernel uses for both spins, and
Hoodbhoy–Jaffe–Manohar NPB 312 (1989) 571 for Δ; for J = 3/2, nothing. **Ask:** I.
Cloët (ANL) and W. Cosyn — classify the basis and fix its normalization — a
co-authorship offer, not a favour. **Meanwhile:** rank ≤ 2 truncation, rank-2 shapes
as scenarios, and `_tensor_moments` declared the adopted convention, not a derived
result.

**Also needed:** an identity between the rank-2 moment used by the ion-source and
target literature at spin 3/2 and the one `polligen.spin` normalizes as T = ⟨3J_z² −
J(J+1)⟩/3, which is +1 for a pure |m| = 3/2 fill and −1 for a pure |m| = 1/2 one
against P_zz ∈ [−2, +1] at spin 1; and the same for rank 3, `spin.octupole_moment`
being +1 for pure m = +3/2 but −3.000 for pure m = +1/2 and so not a bounded
polarization fraction. **Blocks:** the ⁷Li run plan and every projection resting on
it. The ⁶Li plan alternates m = ±1-rich fills at +P_zz against m = 0-rich fills at
−2P_zz and `bookkeeping.tensor_flip_plan` is spin-1 by construction; J = 3/2 has no m
= 0 state, so the ⁷Li pattern is |m| = 3/2-rich against |m| = 1/2-rich — two fills
differing in the *sign* of T, whose acceptance cancellation must be re-derived. The
literature's "⁷Li accelerated at P_zz = 0.69" (Jänsch, NIM A 254:7) is in that paper's
convention, not this T. **Exists:** module-internal definitions only. **Ask:** the
EPIOS source authors (Atoian et al., PRC 113 (2026) 060501, `refs/2510.10794.pdf`)
with Cosyn — the achievable spin-3/2 population pattern in that convention, and the
flip pattern that cancels acceptance. **Meanwhile:** T as defined, the two-fill sign
flip assumed but not costed, and no ⁷Li tensor-luminosity share published.

## 2. The coherent cos 2φ amplitude beyond the linear-in-|t| expansion

**Needed:** the Fourier coefficient a₂(|t|) of the coherent yield per m-state for an α
+ t cluster density over 0.017 ≲ |t| ≲ 0.25 GeV² at the EIC coherent x_P, with its
sign relative to the deuteron and without the linear expansion. **Blocks:** any ⁷Li
coherent number, and arithmetic rather than taste says so.
`coherent.CoherentScenario.eps_b0` is ΔB₀/B, anchored at ≈ 0.21 on the
polarized-deuteron IP-Glasma calculation and scaled through the quadrupole dial to
−0.08 for ⁶Li; the same scaling gives 2.23 from the deuteron anchor and 4.46 from the
⁶Li default, and at B(⁷Li) = 51.13 GeV⁻² the implied |ΔB₀| is 114 or 228 GeV⁻², larger
than B itself — J = 3/2 has no m = 0 state, and the deuteron m = 0 slope the
parametrization rescales is driven negative. The coefficient −(P/2) ε_B0 B |t| then
reaches unity at |t| = 0.0175 and 0.0088 GeV² and 14.3 or 28.5 at the upper edge of
the tagged window, against 0.100 at 0.05 and 0.500 at 0.25 GeV² for ⁶Li, which sits
inside the regime its anchor validates; a coefficient above one is a negative
probability density. The conclusion is not that ⁷Li modulates fifty times more
strongly but that the expansion has no domain of validity here. **Exists:**
Mäntysaari–Salazar–Schenke–Shen–Zhao, PLB 858 (2024) 139053 (`refs/2408.13213v1.pdf`),
digitized as `coherent.MANTYSAARI_A2_DEUTERON` — the only calculation of this
observable for any nucleus, and it is the deuteron at 6% D-wave; below Ca no
light-nucleus coherent-fraction prediction exists. **Ask:** the Mäntysaari–Schenke
group — rerun their polarized-deuteron IP-Glasma setup with an α + t density,
reporting a₂(|t|) per m-state rather than a slope modulation, and say whether the
JIMWLK washout of deformation toward small x (arXiv:2411.14934) removes the effect
over the EIC x_P range. **Meanwhile:** no ⁷Li coherent scenario and no ⁷Li coherent
number anywhere; the deformation-independent gluon-transversity term of the same
scenario would survive, but alone it would misrepresent what the channel is for.

## 3. The α + t P-wave overlap with its m-dependence, and the intact form factor

**Needed:** the two-cluster overlap ⟨α t|⁷Li, m⟩ in relative momentum k and in the
angle between k̂ and the spin axis, out to k ≈ 0.6 GeV/c; separately, the intact-⁷Li
elastic form factor in the coherent window, presently a Gaussian of the charge radius.
**Blocks:** every tagged ⁷Li asymmetry, through its radial half. The angular structure
is free — the pure P wave gives ⟨P₂(cos θ_k)⟩ = −T/5 (verified as −0.1999 at M = 3/2,
+0.1998 at M = 1/2), the in-situ alignment polarimeter of `plans/05` §5.2 — but the
radial input is not, and it is where the signal lives: over the band β = 0.20 / 0.30 /
0.40 GeV of the two-parameter P-wave form the mean spectator momentum moves 0.1113 →
0.1333 → 0.1505 GeV/c and the fraction above 0.3 GeV/c by a factor 3.6, 0.0231 →
0.0522 → 0.0836, in exactly the region where Cosyn–Weiss find tagged tensor
asymmetries of order unity ("spectator momenta ≳ 300 MeV, which select configurations
with large D-wave", `refs/2603.23700.pdf`). The model band is wider than the effect.
**Exists:** the QMC/VMC effective polarizations we use, P_p = +0.866 and P_n = −0.037
in `polli_fastsim.beams.LI7`, are Wiringa et al. arXiv:1309.3794 via JLab E12-14-001 —
a reference with no `refs_dict.json` entry and no local copy, the first thing to fix
when this item opens. **Ask:** R. B. Wiringa (ANL, local) — the α + t P-wave overlap
with its m-dependence tabulated in k, and the α + d S and D overlaps for ⁶Li in the
same format; `spectator.sample_k` inverts a tabulated cdf and would take a VMC n(k)
unchanged, while `polligen.tagged.Wave` carries hard-coded analytic radials and needs
a one-field extension to accept one, so the upgrade is close to a data drop.
**Meanwhile:** the two-parameter forms with the β = 0.20–0.40 band quoted on every ⁷Li
tagged number.

## 4. Coherent excitation of the 477.6 keV 1/2⁻ state

**Needed:** the cross section for coherent excitation of ⁷Li's first excited state
relative to elastic scattering off the ground state, and its |t| dependence.
**Blocks:** the interpretation of any intact-⁷Li coherent yield. The level lies far
below the α + t threshold (S = 2.468 MeV, AME2020), so it is particle-stable and
M1-decays to the ground state: a coherent excitation delivers an intact-looking A = 7,
Z = 3 recoil at R = 1.000 plus one photon, and it carries the ground state's isospin,
so the isoscalar-exchange protection that shields the analogous ⁶Li background
(`plans/06` §6.2) does not apply. In the α + t picture the 3/2⁻ ground state and the
1/2⁻ level are the spin-orbit partners of one L = 1 spatial orbit, so a rank-2 spatial
operator connects them by a 6j recoupling — our inference, and one for the answer to
confirm or refute: if it holds, the excitation is driven by the *same* gluonic
quadrupole amplitude item 2 asks for and is not an independent background at all. The
photon is a weak handle: at γ = 43.72, 106.60 and 126.32 (m(⁷Li) = 6.5338 GeV at 40.8,
99.5 and 117.9 GeV/u) its maximum lab energy 2γE* is 41.8, 101.8 and 120.7 MeV against
a B0 calorimeter threshold quoted at 50–100 MeV over 5.5–20 mrad — below it at the low
configuration, straddling it at the other two, with the isotropic decay putting half
the photons below γE* anyway. One point is favourable: a J = 1/2 recoil carries no
rank-2 moment, so the contamination dilutes the modulation rather than faking one.
**Exists:** the TUNL A = 5, 6, 7 evaluation (Tilley et al., NPA 708 (2002) 3), except
that the local copy `refs/TUNL_A6_2002.pdf` is the A = 6 chapter alone and carries no
A = 7 table — it must be fetched before the widths quoted here count as repo-verified;
no calculation of coherent nuclear excitation off a light nucleus is known to us.
**Ask:** Cosyn and the small-x group of item 2 — the inelastic coherent channel
alongside the elastic one, at least as a ratio; and, as an ePIC far-forward question,
whether the B0 photon threshold is a floor or a lowerable design value. **Meanwhile:**
no veto assumed, no subtraction, and no ⁷Li coherent number published.

## 5. FSI and two-body currents in the tagged α + t channel

**Needed:** the correction to the impulse approximation for an α spectator recoiling
against DIS debris in ⁷Li, and the size of two-body-current contributions to the
*tensor* tagged observables, across the k range of item 3. **Blocks:** the systematic
on every tagged ⁷Li asymmetry. Those are quoted above 300 MeV/c of spectator momentum,
where both corrections are largest and the deuteron pole extrapolation is least
tested; without a bound they are impulse-approximation numbers with an uncontrolled
error, and the alignment polarimeter of item 3 — attractive because it needs no
external polarimetry — inherits it. **Exists:** the deuteron case is solved
(Cosyn–Weiss, PRC 97 (2018) 035209, pole extrapolation in spectator p_T at fixed
light-cone fraction), with the spin-1 SIDIS formalism in `refs/2603.23699.pdf` and
`refs/2603.23700.pdf`; nothing for a cluster spectator off A = 7. That GFMC needs
two-body currents to reproduce ⁷Li's magnetic moment (Report 0 §2.2) is why that
sector should matter here. **Ask:** W. Cosyn and M. Sargsian — extend the pole
extrapolation to a composite spectator, or bound the error of not extrapolating.
**Meanwhile:** the impulse approximation, quoted at small |t′| only, with the omission
stated wherever a tagged ⁷Li asymmetry appears.

## 6. Radiative corrections in the tensor sector

**Needed:** a radiative-correction treatment for tensor-polarized targets — the A_zz
rate difference and the cos 2φ modulation — at spin 3/2 as well as spin 1. **Blocks:**
not a central value but an uncertainty; the affected claims carry no correction band
and say so (`plans/07` WP4, Report 2 §7). **Exists:** vector-case tools (DJANGOH,
HERACLES) and, in `polligen.radiative`, the unpolarized collinear-ISR migration —
+0.62 / +0.50 / +0.94 / +1.22% of Δ̂ at the four sweet spots of the published
generator window, ≤ 2.9% once the low-Q² feed-in is opened. An unpolarized QED study
does not bound a tensor one, and at J = 3/2 there is not even the spin-1 starting
point of item 1 to build on. **Ask:** the EIC radiative-corrections community (the
DJANGOH/HERACLES maintainers, the ePIC effort) — whether the tensor sector differs
from the vector one at the accuracy these observables need, before an implementation
is attempted. **Meanwhile:** no band, stated explicitly at every appearance.

## 7. The ⁶Li asks that travel with these

**(a) b₁ and Δ for ⁶Li specifically** (`plans/04` #9). No b₁ prediction has been
published for any A > 2 and no EIC Δ projection exists for any target, so the ⁶Li
flagship rests on scenario shapes — and the ⁷Li version would rest on the same shapes
through a basis that does not exist (item 1). The deuteron b₁ literature it would be
built against is in hand (Cosyn et al., `refs/1702.05337.pdf`; Miller,
`refs/1311.4561.pdf`). *Engage:* Cloët, Cosyn, Miller; lattice, Detmold and Shanahan —
one message with item 1.

**(b) The IP-Glasma α + d run** (`plans/04` #18). Item 2's ask with the ⁶Li density in
place of the ⁷Li one, and the cleaner half of the pair: ⁶Li's near-sphericity makes it
the null test, so a small a₂ there is as informative as a large one at ⁷Li. One group,
one conversation, two densities.

---

Items 1 and 2 bind hardest and are co-authorship offers rather than service requests;
item 3 is local, cheap and would improve numbers we already publish; half of item 4 is
a literature question, and the A = 7 evaluation should be fetched into `refs/` before
anything in it is quoted. Items 5 and 6 move error bars rather than central values.

Two ⁷Li questions are absent because they are not theoretical. Whether a de-squeezed
⁷Li fill is a run at all is a scheduling question for C-AD: the tagging optics is a
net loss for ⁷Li, buying a couple of percent in acceptance against a factor eight to
fifteen in luminosity and inflating every ⁷Li error bar by about three, which
`plans/09` B3 prices per configuration, so ⁶Li and ⁷Li are different runs. And whether
the intact-⁷Li recoil can be tagged at IP6 is
answered and negative: the tagged fraction at the Yellow Report divergences is
1.7×10⁻⁹, 1.3×10⁻³⁵ and 1.5×10⁻¹⁵ at the three configurations against 5.0×10⁻⁷,
8.4×10⁻²⁶ and 3.7×10⁻¹³ for ⁶Li — the intact recoil is the beam at either isotope — so
the coherent channel needs the tagging optics whatever the theory says, and IR-8,
where Chang et al. quote a 17.75% intact-⁷Li efficiency (`refs/2511.05638.pdf`),
belongs to `plans/09`.
