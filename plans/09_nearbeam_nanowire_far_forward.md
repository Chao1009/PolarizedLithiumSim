# Plan 09 — A near-beam layer for the far-forward lithium tags (2026-08-26)

**Question.** The far-forward tags of the polarized-lithium programme — the
intact ⁶Li of the coherent cos 2φ measurement and the α spectator of the
tagged measurements — are *beam-blind*: A/Z = 2 exactly for the intact
nucleus and the α's rigidity is 0.99813 of the beam's, so neither is
separated from the beam by dispersion. What separates them is the
transverse **angle** at the interaction point, and the acceptance is
exp(−B|t|) evaluated at |t| = (θp)² — an exponential in the *square* of
the approach angle times the square of the momentum. Approach distance is
the whole measurement.

Superconducting nanowire detectors, developed for nuclear physics by the
Argonne MEP and Physics Divisions, are the candidate technology.
**Report:** `reports/nanowire_far_forward.html` / `.pdf`.

**Answer: no.** Neither candidate role survives, and both fail on
arithmetic rather than on unknowns (§9.2, §9.3). What survives is the
pricing curve — which is technology-independent — a correction to the
programme's own baseline (§9.4), and a redirected open question #19.

**Status legend:** ☐ todo ◐ started ☑ done · **D-items are decisions or
external inputs, not work.**

---

## 9.0 What a closer approach is worth

> **Superseded 2026-08-28 (plans/10 A4; Report 4 rewritten).** The table
> below was computed with a single energy-independent, isotropic,
> proton-derived σ_θ = 72.7 µrad at every configuration, and neither its
> absolute numbers nor its gains survive the Yellow Report divergences
> (220/380, 180/180, 92/92 µrad h/v for ⁶Li).  The re-derived table, from
> the same two scripts on the rectangular 10(σ_h, σ_v) envelope:
>
> | | 5 × 41 | 10 × 100 | 18 × 275 |
> |---|---|---|---|
> | coherent, silicon (2.50 / 1.51 / 0.53 mrad) | 9.4×10⁻¹⁰ | 2.0×10⁻¹⁹ | 1.2×10⁻⁵ |
> | coherent, YR high-acceptance envelope (2.2 / 1.8 / 0.92 mrad) | 7.2×10⁻⁸ | 6.2×10⁻²⁷ | 7.1×10⁻¹⁴ |
> | coherent, tagging-optics envelope (0.36 / 0.19 / 0.12 mrad, L/L_HA = 1/6.8 / 1/12.8 / 1/9.5) | **0.37** | **0.25** | **0.33** |
> | ⁶Li α tag (routed), silicon / YR envelope / tagging | 0.017 / 0.018 / 0.32 | 0.016 / 0.016 / 0.22 | 0.029 / 0.025 / 0.29 |
> | chain at the tagging optics, pots at silicon → pots follow | 0 → 7 bins, N_tag 2.3×10⁶/yr | 0 → 7, 2.4×10⁶ | 0 → 7, 6.0×10⁶ |
> | | (at 18 × 275 the silicon run's one populated bin is dropped by `nearbeam_reach_gain.MIN_TAGGED_PER_BIN` = 1000 against 231 recoils/yr) | | |
> | δa_t, best bin, pots follow | 0.0037 | 0.0040 | 0.0023 |
>
> At the published optics the *machine* binds at 10 × 100 and 18 × 275,
> where a closer approach buys nothing; at 5 × 41 the silicon re-measured
> on 2026-08-28 binds instead (2.50 against a 2.20 mrad envelope) and a
> closer approach is worth a factor 77.  At the tagging optics the silicon
> binds everywhere and a layer that follows the envelope is the whole
> difference between no tag and a 25–37% tag.  The energy ordering reverses: the top
> configuration is the most productive at equal luminosity.  The three
> effects listed below (statistics, a configuration, bias) were features
> of the 72.7 µrad envelope against the silicon and are kept as a record.
>
> Two readings of the ⁶Li α row, added 2026-08-28 with B2.  It is the
> **pure spectator model** (`spectator.spectator_lab_kinematics`, one
> partial wave per channel); the tagged generator behind money plot 4
> carries the S + D expansion of the same channel and reads 2.5–2.9%
> instead of the 1.5% Roman-Pot slice of this row at the Yellow Report
> optics — the two masks compared like for like — the D wave being the
> entire difference and also the tensor observable itself (B2).  And an
> acceptance is not a measurement: what the tagging optics buys the tagged
> observables is **reach**, not rate — at the published optics the α tag
> admits nothing below a spectator momentum of 0.15 GeV/c, under the
> tagging optics a third to a half of the accepted sample is there, for
> a 19% cost in tagged events per year at 10 × 100 and a 1.2–1.8× gain at
> the other two configurations (B2).  For ⁷Li the whole trade inverts and
> the tagging optics is a factor 7.9–14.8 net loss (B3).

`evgen/scripts/nearbeam_aperture_scan.py` (analytic) and
`nearbeam_reach_gain.py` (the full chain: Roman-Pot emulation, two
azimuths, spin-state-sorted 2-D harmonic fit).

| | 5 × 41 | 10 × 100 | 18 × 275 |
|---|---|---|---|
| coherent tagged fraction, silicon aperture | 1.41×10⁻² | 5.12×10⁻⁵ | 1.94×10⁻¹⁷ |
| coherent tagged fraction, 10σ envelope | 3.71×10⁻¹ | 2.91×10⁻² | 1.97×10⁻⁹ |
| gain | **×26** | **×569** | ×1.0×10⁸, still dead |
| ⁶Li α tag, silicon → 10σ | 0.103 → 0.550 | 0.024 → 0.137 | 0.0012 → 0.0054 |
| \|t\| bins that survive, silicon → 10σ | 2 → 3 | 2 → 4 | 0 → 0 |
| δa_t, best bin, silicon → 10σ | 0.015 → 0.0019 | 0.175 → 0.0044 | — |

Three effects, not one:

1. **Statistics.** δa_t improves ×8 to ×40 where the measurement exists.
2. **A configuration.** At 10 × 100 the silicon aperture returns
   a_t = 0.48 ± 1.45 against a truth of 0.32 — not a measurement. At 10σ
   it returns four clean bins. Mid-energy goes from dead to the
   best-covered configuration in the programme. The reconstruction-chain
   report's "the coherent programme is a low-energy programme" is
   **aperture-conditional, not physics-conditional**.
3. **Bias, not only error.** At 5 × 41 the silicon aperture returns the
   gluon-transversity coefficient as −0.0060 ± 0.0045 and +0.0197 ± 0.0044
   against an injected 0.0100 — both wrong by > 2σ, because the emptied
   azimuthal circle leaves the seven harmonic columns nearly degenerate.
   At 10σ every bin is consistent.

**Top energy is not rescued by any technology.** At 825 GeV even a 10σ
approach demands p_T > 0.60 GeV, |t| > 0.36 GeV², exp(−50 × 0.36) = 10⁻⁹.

**The `cut_scale_x` trap.** Money plot 6R carries `cut_scale_x = 2.5` from
the pre-measurement belief in a wide horizontal slot. With the geometric
aperture measured, carrying 2.5 as well imposes a 25σ horizontal
retraction no machine requirement asks for, and it binds *before* the
geometry does — hiding the entire effect. Both new scripts hold the
envelope at 10σ per axis and let the measured geometry be the only
aperture.

## 9.1 The idea is already on record — and neither its stated mechanism nor its fallback survives

Near-beam recoil tagging is not a new application to propose to the
Argonne group. **EIC Yellow Report §14.5.3** lists the four SNSPD
applications the EIC R&D committee identified, the first being *"a Roman
pot detector in the forward region about 35 meters or more from the
interaction point to tag low momentum transfer recoiling ions"*.
**Armstrong's FY22 EIC generic-detector-R&D proposal** — reviewed and
funded — states the mechanism verbatim: *"The so-called 10 σ rule-of-thumb
for Roman pots prevents the beam from damaging the detectors by limiting
its proximity to the beam. Radiation hard nanowire detectors will be able
to move closer to the beam, thus extending the acceptance reach at lower
t."* What this study adds is the number attached to *lower t*, for a
species exponentially more sensitive to it than the DVCS protons that
motivate the rest of the far-forward programme.

**But the stated mechanism does not hold.** The 10σ rule is machine
protection and beam halo, not sensor damage; a more radiation-hard sensor
does not relax it. LHC ran the CT-PPS pots at 25σ at 6.5 TeV.

**And the dead-edge fallback does not survive arithmetic.** A nanowire's
active area reaches within a few hundred nanometres of the substrate edge
(and the SMSPD's edge is measurably sharper than the 10 µm telescope used
to look at it) — a real and unusual property. But slim-edge planar and 3D
silicon already deployed in Roman Pots reaches 100–200 µm, so the
difference is ≈150 µm, or **0.005 mrad at the measured R₁₂ = 29.97 m** —
1% of the 0.53 mrad measured edge, 0.3% of the 10 × 100 gap. Unmeasurable
in |t|.

**The granularity half does survive — and it belongs to ePIC.** At
10 × 100 the blind block is 32 mm wide against a 10σ_x of order 10–17 mm,
and the measured edge there is 1.51 mrad = 32 mm / 21.25 m:
a factor ~330 in coherent yield surrendered to retracting a whole 32 mm
block. The fix is a staircase retraction near x = 0, narrower modules, or
the x-moving layer ePIC is **already designing** — Jentsch, July 2025:
*"Simplest option would be one 'x' moving layer and one 'y' moving layer
in each station … Expect to have this finalized in the next 6-8 weeks."*
This is the highest-value item in the study, and it is a layout change to
a baselined detector, not new technology.

## 9.2 Charge identification: right mechanism, wrong technology

Open question #19 has no answer at IP6. A nanowire supplies a *candidate*
mechanism — but not the one it is usually sold with, and not one that
beats the incumbent.

**Not by pulse height.** The device latches: the amplitude is the
diverted *bias* current. Both Fermilab/JPL beam papers show amplitude
distributions from 120 GeV hadrons, 120 GeV muons, 8 GeV pions and 8 GeV
showering electrons to be indistinguishable, with one threshold serving
all (arXiv:2510.11725, 2410.00251). Any pulse-height scheme is dead.

**It would be by firing threshold.** In the normal-core hot-spot regime,

    r_s = sqrt(Q / (e π c ρ (T_c − T_0)))     I_th/I_c = 1 − 2 r_s / w

— Argonne's Eqs. 1–2 (arXiv:2312.13405), the same relation Renema writes
as E = (w/C)²(1 − I_b/I_c)² and the ion literature as
I_th/I_c = 1 − (zeV)^½ C/w (Cristiano et al., SUST 28 (2015) 124004).
Since dE/dx ∝ z² at fixed β, and ⁶Li at 137.5 GeV/u has βγ = 148 against
128 for Argonne's calibration proton, **r_s ∝ z linearly**: 134 / 268 /
402 nm for d,p / α / ⁶Li.

| wire width | I_th/I_c: p,d | α | ⁶Li | separated turn-ons |
|---|---|---|---|---|
| 250 nm | 0.00 | 0.00 | 0.00 | none — **ANL's MIP optimum** |
| 800 nm | 0.67 | 0.33 | 0.00 | α from d |
| **1000 nm** | **0.73** | **0.46** | **0.20** | **all three — the existing microwire** |
| 1500 nm | 0.82 | 0.64 | 0.46 | all three |

Two bias points at 0.33 and 0.60 I_c would tag Z by the firing pattern
alone, both below the 0.80 I_c dark-count wall. In its favour: that
operating point is also the one that makes the device blind to a 300 K
beam pipe (a 300 K surface radiates ≈3×10¹⁸ photons cm⁻² s⁻¹ above the
0.04 eV threshold), and the wide wire Z-ID needs is the same microwire
the area problem independently needs.

### Why it fails anyway

**(1) The anchor is softer than it looks.** r_s = 134 nm is the
*extrapolated* zero-crossing of a four-point fit, and the authors write
*"While the physical validity of this simple model is a question of
future work"*. Inverting their four points individually gives 102–118 nm
(mean 113). And the volume in which Q is counted differs by a **factor
5000** between their own two papers: the 2024 paper says *"the energy
that the proton has deposited into the thin film"* (≈20 eV), the 2026
paper applies the same √Q scaling to the *substrate* deposit (0.1 MeV).
Both scale as z², so r_s ∝ z survives; the absolute radii do not.

**(2) "An interpolation, not an extrapolation" is only half true.**
Argonne's ²⁴¹Am α differs from their 120 GeV proton almost entirely
through 1/β² (β = 0.054 against ≈1), not through z². It tests the √Q law
across energy; it says **nothing about z² at fixed velocity**. Nobody has
ever varied Z at fixed β on one of these devices, and every published ion
demonstration varies kinetic energy through an acceleration voltage, so
charge enters as E = zeV — linear in z, not z² from dE/dx.

**(3) The incumbent carries more information — but far less decisively
than the first version of this plan claimed.** EICROC provides per channel
an **8-bit 40 MHz SAR ADC** for charge (the ToT of ALTIROC was replaced by
it precisely for dynamic range), behind an AC-LGAD with a **30 µm active
thickness** at 500 µm pitch.

*Correction, 2026-08-27.* The first version quoted "4.8σ per plane" and
treated 8 bits against 1 bit as decisive. **A σ is the wrong figure of
merit**: gap over the quadrature sum of two Landau core widths is neither
a PID separation power nor a fake rate, and what puts an α inside a ⁶Li's
window is the Landau *upper tail*, which falls as 1/λ. Restated as a fake
rate at a matched efficiency (`nearbeam_zid_power.py`, a sampled Landau,
4 planes, 95% ⁶Li efficiency, 1.5×10⁶ events):

| readout | α fake | efficiency reached |
|---|---|---|
| 8-bit charge, per-plane likelihood ratio (**the optimum**) | 2.3×10⁻⁵ | 0.950 |
| **one bit per plane, majority-of-k, 99% efficient** | **3.1×10⁻⁵** | 0.950 |
| truncated mean — the standard analogue dE/dx estimator | 2.7×10⁻³ | 0.950 |
| plain sum of the four planes | 5.3×10⁻² | 0.950 |

**One bit costs a factor 1.4, not orders of magnitude.** The two species
are far apart (MPV 31.7 against 75.2 keV) and the power comes from
requiring *coincidence* across planes, not from precision within one. Note
also that a truncated mean is **86× worse than one bit** and a plain sum
worse still — a Landau has no mean, so one δ ray drags it. More bits does
not automatically mean better Z-ID.

**Where the nanowire actually loses is geometric fill factor.** The
coincidence that makes one bit sufficient needs every plane to record the
track. A silicon pixel plane does so ~99% of the time; a superconducting
wire comb only over its fill — 25% (arXiv:2510.11725), 40%
(arXiv:2410.00251), 50% (the ANL EIC-targeted device, arXiv:2312.13405):

| per-plane efficiency | best reachable ⁶Li efficiency, 4 planes |
|---|---|
| 0.99, silicon | ≥ 0.95 ✓ |
| 0.50, ANL device | 0.937 — **cannot reach 95% at any working point** |
| 0.40 | 0.870 |
| 0.25, SMSPD | 0.683 |

That is a **fabrication** number, not an information-theoretic one, and it
is a fairer and more actionable thing to put to the MEP group than "you
only have one bit". It is also the one thing on this list they could
straightforwardly change.

### The handle that beats both, and is free

The background #19 exists to reject is ⁶Li → α + d. Both fragments sit
within 0.5% of beam rigidity (α at 0.99813, d at 1.00452 — mass-to-charge
ratios, not mass-number ratios), deep inside the ±5% near-beam band, so
neither is separated by dispersion. What does separate them is that the
breakup's relative momentum is *transverse* and is therefore not boosted:
the α and the d take opposite kicks of the same k_T.

**Two momentum scales, and the first version of this plan mixed them.**
√(2μQ) with μ = 1247.7 MeV gives **60.7 MeV/c** for Q = the α+d separation
energy 1.4743 MeV — the *bound-state* scale κ, which is what
`spectator.sample_k` samples — and **42.2 MeV/c** for Q = 0.712 MeV, the
decay momentum of the 2.186 MeV 3⁺ resonance. Both are physical, for
different production mechanisms.

Since 2026-08-28 both fragments of one breakup are sampled together
(`spectator.breakup_lab_kinematics`: one k, the spectator boosted with +k
and the partner with −k), which the two independent channels could not do.
The α is carried by 4p_u and the d by 2p_u, so that common kick opens lab
angles in the inverse ratio of the fragment **masses**, θ_d → 1.987 θ_α as
k → 0, at azimuth φ_α + π. The ratio is not fixed away from that limit (it
runs from 0.7 to 5 over the sampled k), but the deuteron is the wider
fragment in all but 3 × 10⁻⁵ of breakups, so it is the easier one to catch
and the α is the binding one. `evgen/scripts/nearbeam_two_hit.py`, 4 × 10⁵
breakups per configuration for the quantiles and 1.2 × 10⁷ for every
fraction:

| optics | median separation | 16–84% | in 500 µm pixels | both recorded (YR / tagging) | α fakes a tag (YR / tagging) | partner d seen (YR / tagging) | recorded pair in one pixel (tagging) |
|---|---|---|---|---|---|---|---|
| 18 × 275 | **10.9 mm** | 3.2–26.9 | 22 | 0.0002 / 0.22 | 0.0007 / 0.27 | 0.25 / **0.84** | 4.2 × 10⁻³ |
| 10 × 100 | **10.7 mm** | 3.4–26.4 | 21 | 0.0000 / 0.17 | 0.0001 / 0.21 | 0.02 / **0.82** | 2.2 × 10⁻³ |
| 5 × 41 | **17.3 mm** | 6.3–41.2 | 35 | 0.0002 / 0.25 | 0.0020 / 0.30 | 0.12 / **0.84** | 1.5 × 10⁻⁴ |

("YR" is the high-acceptance setting; the high-divergence one is the same
beam at 5 × 41 and a wider envelope at the other two, where the α fakes a
tag in only 4 × 10⁻⁶–4 × 10⁻⁵ of breakups and no partner was recorded at
all in 1.2 × 10⁷. The
veto column is conditioned on 8.7 × 10³ / 1.4 × 10³ / 2.4 × 10⁴ fakes at
the YR optics — 0.245 ± 0.005, 0.017 ± 0.004, 0.120 ± 0.002 — and on
≈ 3 × 10⁶ each at the tagging optics.)

(The millimetre is a lever plus a dispersion, x = R₁₂·θ·cos φ + D·(R − 1),
with (R₁₂, R₃₄, D) measured PER CONFIGURATION — 19.24 / 4.56 / 0.311,
21.25 / 3.35 / 0.287 and 29.97 / 2.93 / 0.292 m, `farforward.POT_LEVERS`,
D3 closed. They also inherit the cluster density's short-range scale:
β = 0.20–0.40 GeV moves the medians by −12% to +9%, and the veto fractions
by less than 0.03. This table has been wrong three times. The 30.1 / 73.4 mm
published here and in Report 4 until 2026-08-28 were computed at the
rigidity-scaled 20.5 and 50 GeV/u that A0 of plans/10 retired, and were a
factor two high: a stale energy inside a derived millimetre, which the
energy-drift checks of `tools/consistency_check.py` structurally could not
see. The correction the same day then dropped the dispersive term, on the
argument that two fragments within 0.7% of beam rigidity share it — they do
not, because they take opposite k_z and their rigidities move in opposite
directions, and restoring it adds 5 / 23 / 39%. There is now a check that
recomputes these medians and compares them with both documents. The third
error was the levers: 38.4 / 18.5 / 15.1 mm carried ONE 18 × 275 pair,
R₁₂ = 30.6 m with R₃₄ taken equal to it, at every configuration, and the
B1 scan measured both per configuration — R₃₄ is 2.9–4.6 m rather than 30,
so the vertical half of every separation collapses and the medians fall
again, to 25.8 / 10.7 / 10.9 mm. A fourth correction followed on
2026-08-29: 5 × 41 was still carrying the fallback R₃₄ = R₁₂ = 19.24 m,
because the 29.6 mm insertion shuts the vertical plane and the B1 ladder
had nothing to regress there. Repeating that ladder through a copy of the
compact file whose four Roman-pot section offsets are zero — the silicon
slid onto the beam axis, every field untouched, R₁₂ = 19.18 m against
19.24 as the control — gives R₃₄ = 4.56 m, and the 5 × 41 median falls to
**17.3 mm**. The ordering stops being monotone with all of it: 18 × 275's
larger R₁₂ cancels its smaller angles and the two upper rows tie to 1%.
The earlier 6.7 / 18.4 / 44.8 mm used a single k = 40 MeV/c.)

**The merge is rare rather than impossible, and the dispersion is why.**
In the angular lever alone the back-to-back fragments with θ_d ≃ 2θ_α keep
a *recorded* pair about 3·min(R₁₂ env_x, R₃₄ env_y) apart, the minimum
taken over the axes a recorded α is *seen* to use rather than over both —
the two levers differ by 4.2× at 5 × 41, 6.3× at 10 × 100 and 10.2× at
18 × 275, and at 5 × 41 no recorded α reaches the 3.80 mrad vertical
envelope at either optics, so its 4.56 m vertical lever is not a scale of
anything there — measured minima 20.0 mm at 5 × 41, 11.7 at 10 × 100 and
8.1 at 18 × 275 against a scale of 21 / 12 / 8.1 mm, within 10% either way
because 1.987 is only the k → 0 ratio — and nothing
merges. The dispersive displacement is free to cancel that, and does: a
recorded pair falls inside one 500 µm pixel in 4.2 × 10⁻³ of cases at
18 × 275, 2.2 × 10⁻³ at 10 × 100 and 1.5 × 10⁻⁴ at 5 × 41, and in none at
the Yellow Report or legacy envelopes, where over 1.2 × 10⁷ breakups per
configuration nothing merges at all — the closest recorded pair is 8.0 mm
apart at the Yellow Report envelope and 4.1 mm at the legacy one, both
tens of pixels wide. At that level the Σz² = 5 against
9 argument is a backstop and not a load-bearing one. The single-hit topology is **asymmetric**: at the
Yellow Report optics 6–7.5% of breakups put the deuteron alone on a pot
against 1.5–1.7% for the α alone.

**An intact ⁶Li is one hit; the breakup is two, tens of pixels apart, in
sensors that already exist.** Conditioned on the α faking a coherent tag —
the only way this channel is a background at all — the partner deuteron is
recorded in **84%** of events at the tagging optics, so hit multiplicity
alone rejects five sixths of it with no dE/dx and no Z-ID. At the Yellow
Report optics the veto is worth 0.02–0.25, but the fake rate there is
10⁻³–10⁻⁴: the optics that creates the background supplies the veto
with it.

**The caveat is the pot's outer edge, and it is a station-layout question
(B1, D3).** The 84% assumes acceptance out to θ = 5 mrad, which is
96 / 106 / 150 mm at the measured per-configuration R₁₂ and an
acceptance-table convention rather than a sensor. The MEASURED outer edge
is 2.85 / 3.85 / 4.00 mrad — past it the primary track stops and the pots
see debris from an ion that has struck the pipe — and at that edge the
veto is 0.80 / 0.84 / 0.84, so the convention costs five points at 5 × 41
and nothing at the other two. Against an
outer edge of 2.0 / 1.0 / 0.5 mrad — 38 / 19 / 10 mm at 5 × 41,
42 / 21 / 11 mm at 10 × 100 and 60 / 30 / 15 mm at 18 × 275, the last of
these one 16 mm module — the veto is 0.68 / 0.22 / 0.00 at 5 × 41,
0.82 / 0.69 / 0.21 at 10 × 100 and 0.84 / 0.81 / 0.56 at 18 × 275: the
fractions are decided in angle and no lever enters them, only the
millimetres beside them. At 5 × 41 the median
separation given an α tag is 35 mm, 19 mm at 10 × 100 and 17 mm at
18 × 275.

## 9.3 What stands in the way

None of it is a physics objection; all of it is load-bearing.

| obstacle | the number | status |
|---|---|---|
| **cold-stage power** | NbN needs 2.8–3 K, WSi microwire 0.8 K; SNSPD cold-stage loads are 100s of µW to a few mW. A 4 K pulse tube gives 1.5 W at 4.2 K and **zero at its 2.8 K no-load floor**, on a 190 kg water-cooled compressor drawing 10.7 kW. | **hard** |
| **beam-induced RF heating** | A ferrite-shielded TOTEM Roman Pot dissipates ≈ 40 W at the LHC; scaling Q²M/σ_z^{3/2} to EIC gives ≈ 16 / 8 / 0.8 W. Even 41 GeV exceeds a pulse tube's entire 4.2 K lift. | **hard** |
| **in-vacuum precedent** | None found for a biased superconducting detector inside an accelerator's *beam* vacuum. LHC CryoBLM (1.9 K) and the CERN AD SQUID CCC both sit outside it, and the AD beam is ~10⁷ antiprotons against the EIC's 8×10¹³. | open |
| **vibration** | Cryocooler cold-head displacement is specified at 7–9 µm; decoupling brings it to 10 nm but is heavy and stiff — the opposite of a pot on motorised rails. | mitigable |
| **radiation hardness** | Unmeasured. The FY22 proposal's central deliverable was *"an upper limit for the accumulated dose … This currently unknown limit"*; four years on nothing is published. The 2023 Quantum Sensors for HEP report: SNSPDs *"are expected to be radiation hard … Although published studies in this area are lacking"*. | **unmeasured** |
| **area** | The largest SNSPD/SMSPD run in a GeV beam is 2 × 2 mm², 8 channels. A 3 mm-deep strip over 20 mm vertical, both sides, is 120 mm² — 30× the state of the art. As 30 × 30 µm² nanowires that is 133,000 channels; as 1 × 1 mm² microwire tiles, 960. | funded R&D |
| **latching** | A hard failure, not a dead time. Argonne's ²⁴¹Am α run is the reassurance: a ≈ 1 µm hot spot did *not* latch a 100–200 nm wire. | addressed |
| **rate** | **Not a blocker, and now sourced.** The ePIC Preliminary Design Report puts pot rates at a few Hz/channel normally and *"30–50 kHz at ~10σ"* from halo, on STAR experience — four orders below the ~10⁸ /s a nanowire's fall time allows. | **closed** |

**Cryogenics inverts the packaging argument.** Windowless in-vacuum
operation and a sub-µm dead edge are exactly what a near-beam sensor
wants. But a device that cannot see a 300 K surface must be enclosed in a
cold shield, and a cold shield is material in front of the sensor. The
resolution — bias below the photon threshold — is available, and has never
been demonstrated next to a beam.

## 9.4 A baseline correction: the ePIC pot geometry has moved — and has been re-measured

Independent of nanowires, and the most immediately actionable finding
here. The aperture `tools/fullsim` first measured came from the
September-2024 `epic-main` inside `jug_xl-nightly`. Reading the current
`main` of `eic/epic` directly (2026-08-26), and then measuring in it
(2026-08-28, B1 below):

| | September-2024 (what was measured) | current `main` |
|---|---|---|
| module | 32 × 32 mm | **16 × 16 mm** |
| insertion | energy-*independent*: 3.2 cm outer, +0.7 cm central | **per-energy 10σ offsets** in `beamline_*.xml` |
| RF shield | 1 mm aluminium per module face, active | **commented out** — *"we don't know if we will even need it … Oct. 2025"* |
| implied blind block at 18 × 275 | 32 mm (x) × 7 mm (y) | 16 mm (x) × 2.7 mm (y) |

The old block gives 32 mm / 30.6 m = 1.046 mrad against the 1.03 mrad the
first scan measured — agreement to 1.5%, which is the check that the
measurement and the file reading confirm each other. **The correction runs
in opposite directions at the two ends**, which is why it could not be
guessed, and the re-measurement confirms the prediction in both
directions. At 18 × 275 the 16 mm block **measures 0.53 mrad** against the
≈0.52 predicted from the files — *below* the 0.9169 mrad 10σ envelope, so
the sensor package no longer binds and the beam does, by eight orders of
magnitude (1.2×10⁻⁵ against 7.1×10⁻¹⁴). At 5 × 41 the per-energy insertion moves the other way, to a
29.6 mm inner edge, because the pots now retract properly for the larger
low-energy beam: the aperture there **measures 2.50 mrad**, worse, and now
the binding constraint at 1.14× the 2.20 mrad envelope. That is precisely
the configuration at which the coherent programme was said to survive.

**The re-measurement is done** (B1 below, `tools/fullsim/README.md`), and
what it replaces is recorded in §9.0: at the γ-matched momenta the
September-2024 aperture gave **9.8×10⁻⁷** at 5 × 41 and the current one
gives **9.4×10⁻¹⁰**, against a 10σ envelope at 7.2×10⁻⁸ that no longer
binds there.  (The 1.4×10⁻² this paragraph named until 2026-08-28 was that
aperture at the retired rigidity-scaled momenta.  The low-energy conclusion
is withdrawn in plans/04 #20; the ⟨cos 2β⟩ = −0.77 and the 37.7% → 1.42%
of the reconstructed chain were computed against the 2.0 mrad edge and
still need a re-run.)  The *curves* of §9.0 are unaffected — they price
every aperture, which is why they are plotted that way; only the marker
labelled "you are here" moves.

## 9.5 Ordered work

### A1 — the aperture scan as a repository script ☑
`evgen/scripts/nearbeam_aperture_scan.py`.

### A2 — the chain, run at both apertures ☑
`evgen/scripts/nearbeam_reach_gain.py`; `--near-beam-mrad` added to
`money_cos2phi_coherent_reco.py`.

### A3 — the sensor budget ☑
`evgen/scripts/nearbeam_sensor_budget.py`: Bethe deposits in the real
film, the hot-spot threshold model and its figure, the sizing, the
channel count, the α tag per optics.

### A4 — the report ☑
`reports/nanowire_far_forward.html` / `.pdf`.

### B1 — re-measure the pot aperture on the current `epic-main` ☑
Done 2026-08-28 in `eic_xl-nightly.sif`, `epic-main` at git 9aaa2969, and
written up in `tools/fullsim/README.md`; 216 map points in p_T × azimuth
plus four 0.05 mrad ladders of 468 points per configuration.

**The aperture** (`reco.RP_APERTURE_MEASURED`, half-widths c_x / c_y):
2.50 / 6.49 mrad at 5 × 41 — the vertical plane is *shut*, the insertion
holding the central silicon off to |y| ≥ 29.6 mm, and 6.49 mrad is that
29.6 mm over the R₃₄ = 4.56 m the zero-insertion scratch geometry measured
on 2026-08-29, so c_y says by how much the plane is shut and not where an
edge is; it replaces the 8.84 mrad bound the entry carried while R₃₄ was
unmeasured, and 6.49 < 8.84 means the retired number was the
anti-conservative one — 1.51 / 2.12 at
10 × 100 and 0.53 / 0.92 at 18 × 275, against the September-2024
2.00 / 3.0, 1.35 / 3.0 and 1.03 / 2.3. The cutout is **taller than it is
wide wherever both axes are open**, by 1.4–1.7×, so `rp_measure`'s
`cut_scale_xy = (2.5, 1)` still has the aspect the wrong way round, now by
3.5–4.4× rather than 5.8×; the *sign* of the acceptance-induced ⟨cos 2φ_t⟩
is unchanged.

**The levers** (`farforward.POT_LEVERS`, station 1): R₁₂ = 19.24 / 21.25 /
29.97 m, R₃₄ = — / 3.35 / 2.93 m, D = 0.311 / 0.287 / 0.292 m, plus the
second-order dispersion D₂ = −0.190 / −0.206 / −0.215 m. The transport is
6.3× stiffer in x than in y at 10 × 100 and 10.2× at 18 × 275, which is why the open aperture is a horizontal
slot although the pots insert vertically, and R₃₄ is not measurable at
5 × 41 — a result, not a gap. D confirms the repository's 0.30 m to 4% at
every configuration.

**The outer edge** (`farforward.THETA_RP_OUTER_MEASURED`): the debris-free
contiguous primary track stops at 2.85 / 3.85 / 4.00 mrad, not the 5 mrad
the acceptance tables assume nor the 144 mm / R₁₂ = 7.5 / 6.8 / 4.8 mrad
the module arithmetic gives; past it the ion has struck the pipe.

**The over-rigid triton is taggable.** An R = 1.286 triton is on the
Roman-Pot silicon in 60 of 60 events at every configuration, inner side of
the bend, dx = +66 mm at station 1 and +70 to +72 at station 2, in the
48–144 mm outer band, which carries no vertical insertion, so the hit is
at y ≈ 0 and needs no p_T; a ZDC deposit follows in 80 / 83 / 98% of the
same events. `farforward.route_charged` has an over-rigid branch (code 6,
RP-inner) and `coherent.LI7_BREAKUP` its ⁷Li table, so the ⁷Li α + t
double tag is available at IP6 and the t tag runs 78 / 92 / 94% rather
than 0.033 / 0.004 / 0.005. The branch is per configuration through the
blind half-width 48 / 32 / 16 mm. Its limits: the model uses the BEAM's
R₁₂ where a stiffer fragment's is up to 8% larger, and there is a real
acceptance hole between θ_x = −1.55 and −2.53 mrad that no angular routing
can see.

**The ⁷Li α tag at 5 × 41 is a knife edge.** At 10 × 100 and 18 × 275 the
α at R = 0.857 is a clean 100% hit at dx = −45 to −50 mm; at 5 × 41 the
32–48 mm band is held off to |y| ≥ 18.0 mm, the α misses station 1
entirely, and the tag survives only through station 2, at −48.3 mm, a
third of a millimetre past the boundary into the offset-free outer band:
**57%, station 2 only, on a module edge**. The fast simulation routes in
angle and cannot see this.

**Which lattice.** `epic-main` ships two lattices for the 5 × 41 ring
setting and they are not a field scale of one another: the re-matched
proton file (41 GeV) gives R₁₂ = 19.24 m and a 2.50 mrad edge, the
Z/A = 0.5 `beamline_5x41_He4.xml` at exactly the ⁶Li fill point (82 GV)
gives 29.81 m and 1.61 mrad (48 mm / 29.81 m, the threshold the code
applies; the first hit lands at 1.60 mrad on +x and 1.70 on −x) — ×1.55 in
the lever and ×0.64 in the aperture. The alternative is carried in code as
`farforward.POT_LEVERS_LIGHT_ION_LATTICE` and
`reco.RP_APERTURE_MEASURED_LIGHT_ION_LATTICE`, both pinned in
`evgen/tests/test_nearbeam.py`, so a downstream caller can price it. Every
number above is on the ePIC baseline compact files, and that is the stated
choice: the ⁶Li fill is taken to run in the ePIC baseline lattice of that
ring setting at twice the field for the Z/A = ½ fill — the same magnets and
the same orbit, which is the transport this scan measured, R₁₂ = 19.24 and
21.25 m — with the Yellow-Report-scaled Z/A = 0.5 file (R₁₂ = 29.8 m, a
1.60 mrad edge at 5 × 41) carried as the alternative.  It is an educated
guess and is labelled as one; D3 asks the FF WG to confirm it rather than
to supply it.
The check is not meaningful at 10 × 100, where the only Z/A = 0.5 file is
at a different *energy* (220 against 199 GV), and was not run.

**Caveats.** One event per scan point, no beam divergence, no vertex
spread, hit level only, and the high-acceptance configuration alone (the
high-divergence insertion is present in the XML but commented out).

### B2 — the α-tag and inclusive-tagging gain, through the chain ☑
Done 2026-08-28. The question as written — "re-run the two scripts at the
near-beam aperture" — stopped being well posed when run 12 replaced the
single 72.7 µrad envelope with per-configuration optics: at the published
optics there is no aperture lever (§9.0). What was done instead is the
same three-optics menu the rest of the study uses, pushed through the two
tagged *observables*. `money_tagged_azz.py` and `tagged_polarimetry_7li.py`
were the last two published figures still routed through the retired
proton-derived 73/164 µrad pair, and they were also applying it as a
**circle**, because `tagged.boost_spectator` never exposed the
spectator's lab azimuth and `TaggedSampler.sample_category` never passed
one to `route_charged`. Both are fixed: `boost_spectator` returns
`phi_spec`, the sampler routes with it, `TaggedSampler(optics=None)`
defaults to `farforward.yr_optics(beam_config)` instead of a module
constant that cannot know the beam, and both scripts take `--config
{0,1,2}` and `--optics {menu,legacy,high-acceptance,high-divergence,
tagging}` with the luminosity fraction in the projection scale.

**The result is a reach statement, not an acceptance one.** At 10 × 99.5,
4×10⁵ events: the ⁶Li α tag is 0.0247 at the Yellow Report optics against
0.2545 at the tagging optics, which at L/L_HA = 1/12.8 costs 19% in tagged
events per year (acc × L = 0.0247 vs 0.0199). What that 19% buys is
where in k those events sit — median accepted spectator momentum
0.348 / 0.323 / 0.334 GeV/c at the three configurations with **0.000
below k = 0.15 GeV/c**, against 0.153 / 0.177 / 0.161 GeV/c with
0.48 / 0.36 / 0.45 below it at the tagging optics (unfolded: median
0.093, 74% below 0.15).  10 × 99.5 is the *worst* case for the trade and
the only one that costs anything: acc × L/L_HA goes 0.0284 → 0.0503 at
5 × 41 and 0.0265 → 0.0331 at 18 × 275, so at two of the three
configurations the reach comes with a 1.2–1.8× gain in tagged events per
year. At every published optics the ⁶Li α tag admits
only the off-rigidity high-k tail; the tagging optics is what turns money
plot 4 from a one-point measurement into a curve, and that answers the
plans/05 question ("does the Cosyn–Weiss O(1) asymmetry at
p_s ≈ 300 MeV/c survive?") as "yes, and at the published optics *only*
there". Dropping the azimuth would have read 0.51 at the tagging optics
instead of the 0.30 it gave with the azimuth on the levers then in force
— ×1.7 — pinned in `test_tagged.py`.

**A correctness defect of money plot 4 was fixed with it.** The figure
overlaid an analytic A_zz^wf curve drawn at θ_k = 90° on a sample the
acceptance sculpts in θ_k, and the two optics sculpt it oppositely: the
surviving off-rigidity window slice is longitudinal (⟨|cos θ_k|⟩ = 0.80)
and the near-beam tail the tagging optics opens is transverse (0.39). The
published ±0.5 swing at k ≈ 0.3 GeV/c was that, not the wave function.
The right panel now carries, per optics, the acceptance-weighted truth
(`tagged.acceptance_weights` × `tagged.azz_tensor_curve`), which
reproduces the folded markers to |ΔA_zz| ≤ 0.063 at the Yellow Report
optics and ≤ 0.085 at the tagging optics in every populated bin.  At the
k = 0.325 GeV/c bin of one `--events 8e6 --config 1` run:
+0.4594 ± 0.0129 against +0.4548 at the Yellow Report optics (0.4σ) and
−0.0974 ± 0.0070 against −0.0945 at the tagging optics (0.4σ), where the
90° curve says −0.482.  (The default 4×10⁵ events give +0.491 and −0.066,
with errors 0.058 and 0.031, and every populated bin within 1.1σ / 1.6σ
of errors running 0.053–0.167 and 0.014–0.100.)
At 8×10⁶ the errors fall to 0.012–0.037 and 0.003–0.023, the bin-centre
residuals to ≤ 0.025 and ≤ 0.022, and one bin — k = 0.175 at the
tagging optics — reads 2.7σ on a 0.009 residual: that is the bin in which
ε(k) turns on (nothing is accepted below k = 0.189 GeV/c at the Yellow
Report optics), so the truth at the bin *centre* is not the truth the
bin-averaged marker measures; averaged over the bin as the marker is,
every populated bin is within 1.6σ (Yellow Report) and 2.2σ (tagging) at
8×10⁶.  These are the numbers
of the tagged kernel *without* an inclusive b₁ (2026-08-29): in the
impulse approximation the deuteron's b₁ is the k-integral of the
m-dependent spectator distribution the sampler already draws from, so an
inclusive b₁ on top counts the same physics twice — a −0.008 offset with
the old toy shape (+0.4471 at 8×10⁶), +0.04 once the digitized Miller b₁
replaced it, which is how it was found.  The 90° curve stays as a
labelled reference.  `--optics legacy` re-routes the same sample through
the retired isotropic pair: 4.2% / 2.5% accepted, A_zz = +0.08 and +0.49
at k ≈ 0.33 GeV/c against acceptance-weighted truths of +0.01 and +0.46
(the pre-2026-08-28 figure's −0.08 / +0.43 was routed without the
per-configuration blind block and is not reproduced digit for digit).
The swing between the two optics is the θ_k sculpting, computed from the
wave function without a fit.

**And the two k spectra were reconciled.** `spectator.spectator_lab_kinematics`
gave ⟨k⟩ = 0.107 GeV/c and a 1.5% R < 0.95 slice where `tagged.TaggedSampler`
gave 0.122 and 2.5%. The difference is the **D wave and nothing else**:
`spectator` carries one partial wave per channel (`l_wave` = 0 for α–d),
`tagged` the full S + D expansion; the S-wave radials are identical
(⟨k⟩ = 0.1071 in both), the D wave is hard (⟨k⟩ = 0.2778) and at
P_D = 0.0867 pulls the mean to 0.1219. Neither is a defect — with P_D = 0
the α–d density is m-independent and A_zz^tag vanishes identically, so
the D wave *is* the observable — and the difference is now documented in
`tagged.py`, the manual §4.1 and Report 3 Table 6's caption. For ⁷Li the
two densities are the same pure P wave, and the 0.7 points by which the
B3 table's 0.9690 at 5 × 40.8 exceeds `tagged_polarimetry_7li.py`'s
0.9617 is **two acceptance definitions, not two densities**: the table is
`1 - lost` (any far-forward system), the tagged script's mask is the
Roman Pots alone, and at 5 × 40.8 the B0 takes 1.1% of the ⁷Li α and the
over-rigid inner branch a further 0.2%. Like
for like the pure model's Roman-Pot tag is 0.9568 / 0.9668 / 0.9721, and
restricted to k ≤ 1.2 GeV/c where the tagged model's grid ends it is
0.9626 / 0.9683 / 0.9736 against the tagged sampler's
0.9617 / 0.9678 / 0.9728 — 0.1 point everywhere.  (The grid truncation
alone is worth +0.6 point at 5 × 40.8 and +0.2 at the other two on the
Roman-Pot mask, +0.1 uniformly on `1 - lost`; the B0 fraction is 0.0000
at 10 × 99.5 and 18 × 117.9.)  The script prints both definitions.

### B3 — ⁷Li ☑
Done 2026-08-28 with B2. `tagged_polarimetry_7li.py` takes the same
`--config` / `--optics` and routes with the azimuth; `reco.rp_aperture_for`
is keyed on the **configuration** (`farforward.yr_config_key`) instead of
matching a bare momentum against the ⁶Li list, which is why it returned
`None` for ⁷Li at the top energy (117.9 GeV/u against ⁶Li's 137.5) and no
near-beam scan could be run for ⁷Li at all; the momentum path is kept for
the pre-2026-08-28 callers. `nearbeam_aperture_scan.py --isotope 7Li`
produces `evgen/nearbeam_aperture_7Li.png`.

**⁷Li's physics does not move, and that is the deliverable.** The α is
off rigidity at R = 0.856 — in the *middle* of the Roman-Pot momentum
window — so it is tagged by the window and never has to clear the
near-beam envelope: 0.969 / 0.968 / 0.975 at the Yellow Report optics
against 0.987 / 0.991 / 0.994 at the tagging optics on the pure spectator
model (`tagging_acceptance.py`, Report 3 Table 6, as `1 - lost`) and
0.9617 / 0.9678 / 0.9728 against 0.9800 / 0.9909 / 0.9919 on the tagged
generator's Roman-Pot mask, the two agreeing to 0.1 point once the same
acceptance definition and the same k range are used (B2's last paragraph;
the tagged script also prints `1 - lost`, now 0.9699 / 0.9690 / 0.9751 and
0.9882 / 0.9921 / 0.9942 — that pair rose on 2026-08-28 with B1's
over-rigid routing branch, which stopped discarding the α's over-rigid
Fermi tail and leaves the Roman-Pot mask untouched, and then fell back at
the two lower configurations when the same day's review gave the script
the per-configuration blind block, 48 / 32 / 16 mm, instead of the
18 × 275 default it had been using at all three; forcing
`farforward.over_rigid_route` to `False` returns it to `acc(RP)` exactly,
0.9678 and 0.9909 at 10 × 99.5). The scan is flat to three
points across the whole 0.05–3 mrad aperture axis. Folded ⟨P₂⟩ slope
−0.1947 (YR) and −0.1962 (tagging) against −0.1929 at the retired
73 µrad and the analytic −0.2000; median δA_∥ 0.01152 and 0.01141 at
equal generated statistics (at equal *luminosity* the tagging bars are
2.78 / 3.81 / 3.15 times larger — the figure now says so on the panel).
The ⁷Li tags and the median δA_∥ moved in the fourth decimal on
2026-08-28 with plans/08 D7, which gave the struck triton its second
neutron; the ⟨P₂⟩ slopes did not (`docs/reproduction_manual.md` §4.1).

**For ⁷Li the tagging optics is a strict net loss, and a near-beam layer
is worth nothing.** ×1.02 in acceptance for ×1/7.9, ×1/14.8, ×1/10.1 in
luminosity is a factor 7.9–14.8 net loss, which at equal running time
multiplies every ⁷Li error bar by 2.78 / 3.81 / 3.15. So the inversion of
plans/00 strategic finding 2 is optics-independent **in both
directions**: ⁶Li needs the de-squeeze and ⁷Li must not have it. **⁶Li and
⁷Li therefore want different machine optics and are different runs, not
one fill plan** — a scheduling item for C-AD alongside plans/10 D1
(Report 3 §5, Report 4 §2.1, §7 and the conclusion).

**Open, and deliberately not half-built: a coherent ⁷Li channel.**
`polligen/coherent.py` is ⁶Li-specific by construction — a J = 1 nucleus
whose quadrupole deformation `CoherentScenario.eps_b0` is scaled from the
deuteron's, with `fragment_rigidity` defaulting to beam_A = 6, beam_Z = 3
and a `LI6_BREAKUP` veto table. ⁷Li has J = 3/2 and a quadrupole moment
≈ 50× larger, so its coherent cos 2φ amplitude is a different physics
case with its own model and its own breakup channels (α + t, and the
over-rigid triton the routing sent to "lost" until `over_rigid_route`
put it on the pot planes, B1), not a re-run of the ⁶Li one. Scoping it
needs a theory input first; recorded here rather than in the ☑ above.
Panel (a) of `nearbeam_aperture_scan.py` stays ⁶Li at either `--isotope`
for the same reason, and says so. What that theory input is, and what else the
channel would need, is scoped in B3a below — a memo, not code.

### B3a — a coherent ⁷Li channel: scoping memo (2026-08-28)

**The observable, and why ⁷Li is the positive control.** The coherent
channel measures the cos 2φ modulation of the recoil azimuth, carried in
plans/06 §6.4b as a t-linear *deformation* term scaled from the
polarized-deuteron IP-Glasma calculation plus a flat *gluon-transversity*
term. On ⁶Li the first is bounded to near zero by the anomalously small
quadrupole moment — that is the null test, and a sizable flat cos 2φ
there cannot be nuclear shape. ⁷Li is the same measurement with the
deformation switched on: Q(⁷Li) = −4.00(3) fm² against
Q(⁶Li) = −0.0806(6) fm² (Stone INDC(NDS)-0833 2021, plans/01; −4.06 and
−0.0818 are superseded and neither belongs here). The dial the scenario
uses is Q/R_ch² = 0.0630 (d), 0.01202 (⁶Li), 0.6697 (⁷Li) at
R_ch = 2.1305, 2.589, 2.444 fm — a factor 56 between the isotopes on the
same detector, in the same window, with the same estimator. Null test and
positive control differ by one ion species.

**The polarization variable is not P_zz.** For J = 3/2, `polligen/
spin.py` normalizes the rank-2 moment as T = ⟨3J_z² − J(J+1)⟩/3 ∈
[−1, +1], +1 for a pure |m| = 3/2 fill and −1 for a pure |m| = 1/2 one,
against P_zz ∈ [−2, +1] for spin 1. There is no m = 0 state, so the lever
the ⁶Li run plan rests on — m = ±1-rich bunches at +P_zz alternating with
m = 0-rich bunches at −2P_zz — has no analogue, and
`bookkeeping.tensor_flip_plan` is spin-1 by construction (j = 1.0,
`spin.spin1_populations`). The ⁷Li pattern would be |m| = 3/2-rich
against |m| = 1/2-rich, two fills differing in the *sign* of T rather
than in an enrichment ratio, and the acceptance cancellation of the
spin-state ratio must be re-derived for it before any projection is
quoted. The literature's "⁷Li accelerated at P_zz = 0.69" (Jänsch NIM A
254:7, via plans/01) is in that paper's convention and is not this T.

**The rank-2 kernel is already isotope-generic.** Verified today:
`xsec.InclusiveKernel._tensor_moments` returns, for J = 3/2,
Q_NN = [3m² − J(J+1)]/3 = (+1, −1, −1, +1) over m = +3/2 … −3/2 and
c_eff = 3 Q_NN, so c_eff averages over a fill to exactly 3T —
populations (0.4, 0.1, 0.1, 0.4) give T = 0.600 and ⟨c_eff⟩ = 1.800. An
*inclusive* ⁷Li cos 2φ needs a Δ(x, Q²) scenario and a spin-3/2 run
plan, not new kernel code. What is ⁶Li-specific is
`polligen/coherent.py`, and above all its deformation parameter.

**ε_B0 cannot be rescaled, and the arithmetic is why.**
`CoherentScenario.eps_b0` is ΔB₀/B, the fractional slope modulation of
the m = 0 state, anchored at ≈ 0.21 for the deuteron and scaled through
the dial above to −0.08 for ⁶Li (band −0.04 to −0.13). The same linear
scaling gives 0.21 × (0.6697/0.0630) = 2.23 from the deuteron anchor and
0.08 × (0.6697/0.01202) = 4.46 from the ⁶Li default, i.e. ε_B0 ≈ −2.2 to
−4.5. With B(⁷Li) = 51.13 GeV⁻² from R_ch = 2.444 fm through
`coherent.gaussian_slope`, |ΔB₀| = 114 to 228 GeV⁻² exceeds B itself:
the m = 0 slope changes sign and its form factor grows with |t|. The
observable fails first — `cos2phi_coefficient_deformation` =
−(P/2) ε_B0 B |t| reaches unity at |t| = 0.0175 GeV² for ε_B0 = −2.23
and at 0.0088 for −4.46.  Against the 0.017–0.25 GeV² window the tagged
pseudo-experiments have used since 2026-08-28 the first of those
thresholds falls *inside* the lowest published bin and only the second is
still below the window, so the coefficient is already at or above unity
over essentially the whole of it and reaches 14.3 or 28.5 at its upper
edge, and a cos 2φ coefficient above one is a negative probability
density (the ⁶Li scenario reads 0.100 at 0.05 GeV² and 0.500 at 0.25,
inside the regime its anchor validates).  Widening the window has
strengthened the argument rather than weakened it. The conclusion is not that the ⁷Li modulation is fifty times
larger, but that the linear-in-|t| expansion has no domain of validity at
this deformation: a ⁷Li number would be a formula evaluated where its own
anchor does not hold. That is the plans/04 #18 ask, not a parameter edit.

**The intact tag, and the optics conflict.** An intact ⁷Li recoil sits at
R = 1.000 exactly — it is the beam — so it is beam-blind as the intact
⁶Li is, and slightly worse, the near-beam envelope translating into a
p_T cut ∝ A p_u that is 7/6 larger. With
`coherent.tag_acceptance_angular` at the Yellow Report high-acceptance
divergences and 10σ (the circular form of Report 1 §5.2), at the ⁷Li
slope B = 51.13 GeV⁻² above, the intact-⁷Li tagged fraction is
1.7×10⁻⁹, 1.3×10⁻³⁵ and 1.5×10⁻¹⁵ at 5 × 40.8, 10 × 99.5 and
18 × 117.9 GeV/u; the same call at the ⁶Li default returns 5.0×10⁻⁷,
8.4×10⁻²⁶ and 3.7×10⁻¹³, reproducing the published row and
validating the usage (docs/reproduction_manual.md §4.6). The channel therefore needs the tagging optics of
Report 1 §6.1 exactly as the ⁶Li one does — which collides with B3 above,
where that optics is a strict net loss for ⁷Li (×1.02 in acceptance
against ×1/7.9–×1/14.8 in luminosity, every ⁷Li bar ×2.78–3.81) and where
this plan concluded that ⁶Li and ⁷Li want different optics and are
different runs. A coherent ⁷Li channel is a *third* configuration, ⁷Li
de-squeezed, competing for beam time with the ⁷Li spectator programme it
would otherwise share a fill with: a scheduling question for C-AD
alongside plans/10 D1, settled on paper and not in simulation.

**The breakup side is not implemented, and is not a rescale either.**
With physical masses, `coherent.fragment_rigidity(4, 2, beam_a=7,
beam_z=3)` = 0.8557 and `(3, 1, beam_a=7, beam_z=3)` = 1.2897: the α
lands mid-window in the Roman Pots (R ∈ 0.60–0.95) at every
configuration, and the triton is over-rigid but no longer lost by
construction. `farforward.over_rigid_route`, added with the B1 scan of
2026-08-28, tests the dispersive displacement at the pot plane against
the configuration's blind half-width (48 / 32 / 16 mm) and the last
module at 144 mm, and finds the R = 1.286 triton on silicon in 60 of 60
events at every configuration, so
`fragment_route_label(3, 1, beam_a=7, beam_z=3)` answers "RP-inner
(over-rigid)". What the function cannot do is answer a ⁷Li question from
the ⁶Li defaults it is called with everywhere else: on those an intact
⁷Li is R = m(⁷Li)/m(⁶Li) = 1.166 and comes back "RP-inner (over-rigid)"
at `config="18x275"`, or "lost (over-rigid)" at 5 × 41 where the blind
block is widest, against the correct answer — beam-blind at R = 1.000
exactly — and the ⁷Li α comes back "RP p_T-tail only", the ⁶Li α's
destination, where its own 0.8557 is "RomanPots". Passing `beam_a=7,
beam_z=3` is one argument per call and the destinations are then settled.
What is missing is the amplitude behind them. The table itself now exists:
`coherent.LI7_BREAKUP` beside `LI6_BREAKUP`, reached through
`coherent.BREAKUP_TABLES` and the beam-generic `fragment_route_label` /
`veto_table` (added with B1, pinned in `evgen/tests/test_nearbeam.py`),
its leading channel α + t at S = 2.468 MeV (AME2020, plans/01) with the
triton on the Roman-Pot inner side rather than lost. What it has no
density to fill it from is the third of the theory asks below.

**The background with no ⁶Li analogue: the first excited state.** ⁷Li's
477.6 keV 1/2⁻ level lies far below the α + t threshold quoted above, so
it is particle-stable and M1-decays to the ground state: a coherent
excitation delivers an intact-looking A = 7, Z = 3 recoil at R = 1.000
plus one photon. It carries isospin 1/2, the *same* as the ground state
— T above is the rank-2 spin moment, not this — so
the protection the ⁶Li channel enjoys — isoscalar diffractive exchange
cannot drive T = 0 → T = 1, plans/06 §6.2 (2) — does not apply, and
nothing suppresses the excitation at production. The photon is the only
handle and a weak one: at γ = 43.71, 106.60 and 126.31 (from
m(⁷Li) = 6.5338 GeV and the beam momenta 40.8, 99.5, 117.9 GeV/u) its
maximum lab energy 2γE* is 41.8, 101.8 and 120.7 MeV against a B0 EMCal
photon acceptance quoted at ≥ 50–100 MeV over 5.5–20 mrad (plans/06
§6.5) — below threshold at the low configuration, straddling it at the
other two, with the isotropic decay putting half the photons below γE* in
any case. To be checked before the channel is costed: the TUNL A = 7
evaluation for this level and any other particle-stable one, whether the
B0 threshold is a floor or a lowerable design value, and the coherent
excitation cross section relative to the elastic one. One structural
point is already favourable — a J = 1/2 recoil carries no rank-2 moment,
so the contamination dilutes the modulation rather than faking one, the
structure of the incoherent argument of plans/06 §6.2 (1) — which makes
this a normalization question, not a fake-signal one.

**Buildable now, versus genuinely blocked.** With today's code: the
inclusive spin-3/2 rank-2 projection; the ⁷Li α tag, already simulated at
both optics (B3); and the routing arithmetic, now that
`fragment_route_label` takes the beam arguments and `over_rigid_route`
has settled the triton — an afternoon's work, deliberately not done,
because a routing table with no amplitude behind it is the half-built
channel this section refuses. Blocked on theory, in the order they bind:
the rank-2 slope amplitude *beyond linear order* for a deformation 56×
the ⁶Li one, which is plans/04 #18 with an α + t density in place of the
α + d one; the spin-3/2 rank-2 (and rank-3) structure-function basis,
plans/04 #14, without which the inclusive ⁷Li cos 2φ has a kernel and no
structure function to put in it; and the α + t P-wave overlap with
m-dependence, plans/04 #15, which the tagged and breakup sides need.
*2026-08-29 (author decision (4)):* those three, together with the
excited-state, FSI and tensor-radiative questions and the two ⁶Li asks
that travel with them, are now written out as a circulable note —
`docs/note_7li_theory_questions.md`, seven questions, each with what the
simulation needs, what it blocks, what exists, the ask and its
addressee, and what we assume until it is answered.

**Order of work.** (i) The optics conflict, on paper: it decides whether
the channel is a run at all and costs nothing to answer. (ii) The
excited-state background — a literature question and a B0-threshold
question, needing no new simulation. (iii) The three theory asks,
packaged with the ⁶Li ones rather than sent separately; #18 is one
conversation with one group about two densities.

**Published statements a ⁷Li coherent channel would touch**, listed so
that opening it is a bounded edit; none should be edited before (i) and
(ii) are settled, with one exception below. Report 0 §2.2 ("⁷Li … is the
deformed counterpart"), §4.2's closing paragraph (the spin-3/2 rank-2
sector has no published structure-function basis), §5.4 Table 2 and its
caption (α at R = 0.85571 in the Roman-Pot window, triton at 1.28971
routed as lost) and the Table 3 row "Δ via cos 2φ (⁶Li; ⁷Li eligible)",
whose eligibility is inclusive-only on today's code; Report 1's
introduction and §6.4 (⁶Li's small Q as the null test,
Q(⁶Li) ≈ Q(⁷Li)/50, "the null test unavailable to the deuteron") and §7
(the coherent channel needs the tagging optics); Report 3 Table 6 and its
caption (the ⁷Li α and triton tags, and the R > 1 branch that replaced the
routing assumption).

The exception is the routing, and it was B1's rather than this memo's.
Report 0 §5.4 Table 2 and its caption and Report 3 Table 6 and its caption
routed the over-rigid triton as lost against an "unverified" gun scan;
`over_rigid_route` and the re-measurement settled it, and both tables now
carry the recomputed ⁷Li triton tag, 78 / 92 / 94%.

### B4 — two-hit topology for α + d ☑
§9.2, done 2026-08-28. `spectator.breakup_lab_kinematics` samples both
fragments from one k; `farforward.separation_at_pots` and the new outer
bound `THETA_RP_OUTER` price the pot plane; `evgen/scripts/
nearbeam_two_hit.py` produces the table and the figure. Both fragments are
in acceptance in 22–29% of breakups at the tagging optics and in ≤ 0.03% at
the Yellow Report optics; given an α that fakes a coherent tag the partner
deuteron is recorded in **84%** of events at the tagging optics against
0.02–0.25 at the published ones, where the fake rate is 10⁻³–10⁻⁴
anyway. A
recorded pair merges into one 500 µm pixel in ≤ 4 × 10⁻³ of cases, through
the dispersion rather than the angle. The number the veto really depends on
is the pot's outer edge (B1), now measured at 2.85 / 3.85 / 4.00 mrad,
where the veto is 0.80 / 0.84 / 0.84; at one 16 mm module it falls to
0.00 / 0.31 / 0.56. The pass also corrected the §9.2 separation table,
which had carried the retired 20.5 / 50 GeV/u through a derived millimetre;
the review of the same evening corrected it again, for the dispersive
displacement the first correction had dropped, and re-measured the Yellow
Report veto column on 1.2 × 10⁷ breakups per configuration after its
10 × 100 entry turned out to rest on four events.

### D1 — can the existing AC-LGAD stack do Z-ID? · **ePIC FF WG / OMEGA-IJCLab / BNL**
**Ask this before anything about nanowires.** A Geant4/DD4hep study
through all four RP layers with an EICROC front-end model, plus two
numbers nobody has published: the chip's **input charge dynamic range in
fC** (does 9 × MIP clip the front end?) and the sensor's **gain-suppression
curve** at ~9 MIP. One person-month, no hardware, and it closes #19
either way.

### D2 — what actually sets the 10σ retraction, and σ_x for light ions · **C-AD / ePIC FF WG**
A conversation, not an experiment. The first part decides whether any
"rad-hard ⇒ closer" argument can work at all. The second is worth three
orders of magnitude in coherent yield at top energy by itself: the
repository's envelope implies 10σ_x ≈ 22 mm (0.73 mrad) while Jentsch's
"1σ ~ 1 mm" and the published "inner detector edge … 1 cm or less" imply
~10 mm (0.33 mrad).

### D3 — per-optics R₁₂ and R₃₄ ☑ (the lattice question remains) · **ePIC FF WG**
R₁₂, R₃₄ and D are measured at all three configurations since 2026-08-28
— 19.24 / 21.25 / 29.97 m, — / 3.35 / 2.93 m and 0.311 / 0.287 / 0.292 m
(`farforward.POT_LEVERS`, `tools/fullsim`) — so this plan quotes
millimetres everywhere; R₁₁, R₂₁, R₂₂ and D′ are still unmeasured. The
10σ offsets in `beamline_*.xml` exist only for proton optics and are
marked *"rough extrapolation"* at 5 × 41, and what that turns into is the
**lattice question** this item now puts to the FF WG: `beamline_5x41.xml`
(41 GeV proton) and `beamline_5x41_He4.xml` (Z/A = 0.5, 82 GV) are not a
field scale of one another and give R₁₂ = 19.24 and 29.81 m at the same
ring setting — a factor 1.55 in every 5 × 41 millimetre and 0.64 in the
5 × 41 aperture (`farforward.POT_LEVERS_LIGHT_ION_LATTICE`,
`reco.RP_APERTURE_MEASURED_LIGHT_ION_LATTICE`). We publish on the
baseline file at twice the field, as an educated guess; the ask is
confirmation of it, with the scaled file as the labelled alternative.

### D4 — the x-moving layer and a staircase retraction near x = 0 · **ePIC FF WG (A. Jentsch)**
Where the ×60–×330 actually lives, and already being designed.

### D5 — *only if D1 says the LGAD cannot do it:* the α turn-on curve · **ANL MEP**
Turn-on bias versus wire width for the ²⁴¹Am α, on the same wires as the
120 GeV proton. Costs them one plot from data they already hold, and
kills the concept at zero cost if the α curve is not where √Q predicts.
**Be precise that it tests the √Q law, not z² at fixed β.** The
experiment that would test that — a relativistic light-ion beam on a
1–2 µm microwire at FTBF or the SPS — is a beam-time request, and this
study should not be what motivates it.

---

## 9.6 What this changes in the existing documents

* **plans/04 #19** (RP Z-ID for A/Z = 2) — a candidate answer exists and
  is sharper than expected: the firing threshold, not dE/dx amplitude.
  It rests on one unmeasured curve (D1).
* **plans/04 #20** (RP cutout geometry) — the measurement stands for the
  geometry it was made in, and that geometry has moved (§9.4).
* **The coherent programme's energy reach** — aperture-conditional, not
  physics-conditional (§9.0).
* **The recommendation, corrected.** The first version of this plan
  concluded that charge identification was the strong half of the case
  and should be put to the MEP group ahead of anything about apertures.
  An adversarial review established the opposite: the incumbent AC-LGAD
  already digitises an 8-bit charge over a 30 µm active layer across four
  planes, so a nanowire is a downgrade — though not for the reason first
  given here, since one bit per plane costs only a factor 1.4 in α fake
  rate at matched efficiency and the nanowire loses on geometric fill
  factor instead (§9.2, which retracts the "≈4.8σ per plane" of the
  first version). **Ask the incumbent first (D1), and ask C-AD
  what the 10σ rule actually protects (D2).** Ranked list of what would
  help the lithium tags: (1) a tagging optics, the only lever that turns
  the tag on at all, and with it the TOP configuration, whose four times
  the coherent rate at equal luminosity outweighs its slightly smaller
  tagged fraction — the low-energy configuration led only while a single
  fixed envelope was assumed (§9.0, Report 4 §7); (2) re-tiling
  the near-beam silicon, which ePIC is already designing; (3) using the
  charge EICROC already digitises; (4) two-hit topology for α + d,
  measured at an 84% veto under the tagging optics and 4–29% under the
  published ones, where the α fakes a tag in 10⁻³–10⁻⁴ of breakups
  anyway, and limited by how far the pot stations
  extend (B4); (5)
  the IR-8 secondary focus and its z² Cherenkov. The same tagging optics
  buys the ⁶Li *tagged* observables reach rather than rate — nothing
  below k = 0.15 GeV/c at any published optics, a third to a half of the
  accepted sample there under the tagging one, at a 19% cost in tagged
  events per year at 10 × 100 and a 1.2–1.8× gain at the other two (B2) —
  and buys ⁷Li a factor 7.9–14.8 net loss (B3), so item (1) is a per-isotope
  recommendation, not a programme-wide one. Superconducting
  nanowires appear nowhere on it.
* **Where the technology *would* belong.** The Argonne programme's own
  four-application list points at the cold bore of a superconducting
  magnet, in front of the ZDC, and a Compton polarimeter. The result that
  makes the technology unique — saturated efficiency to 5 T *parallel* to
  the device plane, against 0.5 T perpendicular — is worth nothing at a
  pot in a field-free drift and everything inside a magnet.
