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
difference is ≈150 µm, or **0.005 mrad at R₁₂ = 30.6 m** — 1% of the
0.52 mrad edge, 0.3% of the 10 × 100 gap. Unmeasurable in |t|.

**The granularity half does survive — and it belongs to ePIC.** At
10 × 100 the blind block is 32 mm wide against a 10σ_x of order 10–17 mm:
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

**(3) The incumbent already carries more information — and this is
decisive.** EICROC provides per channel an **8-bit 40 MHz SAR ADC** for
charge alongside its 10-bit 25 ps ToA TDC (the ToT of ALTIROC was
replaced by that ADC precisely for dynamic range), behind an AC-LGAD with
a **30 µm active thickness** at 500 µm pitch.

| discriminant | MPV d / α / ⁶Li | α ↔ ⁶Li per plane | planes ePIC has |
|---|---|---|---|
| AC-LGAD, 30 µm Si, 8-bit charge | 7.2 / 31.7 / 75.2 keV | **4.8σ** | **4** |
| nanowire, 12 nm NbN, threshold | — | 1 bit | 0 |

A nanowire cannot beat this observable; at best it matches it. The limit
on either is the **δ-ray upper tail**, which falls as 1/λ and is nearly
independent of sample thickness (30 µm of Si gives 4.8σ, a 1 µm substrate
volume 3.8σ) — so thickening does not rescue it, but four planes with a
majority vote take a 5–6% single-plane fake to ~10⁻³. **ePIC already has
four planes.**

### The handle that beats both, and is free

The background #19 exists to reject is ⁶Li → α + d. Both fragments sit at
beam rigidity, so neither is separated by dispersion — but the breakup's
relative momentum, k_rel ≈ √(2μQ) ≈ 40 MeV/c (plans/06 §6.2), is
*transverse* and is therefore not boosted. The α at 4p_u and the d at
2p_u take opposite kicks of the same k_rel:

| optics | α | d | separation | in 500 µm pixels |
|---|---|---|---|---|
| 18 × 275 | 2.2 mm | 4.5 mm | 6.7 mm | **13** |
| 10 × 100 | 6.1 mm | 12.2 mm | 18.4 mm | **37** |
| 5 × 41 | 14.9 mm | 29.9 mm | 44.8 mm | **90** |

(k_rel = 40 MeV/c at R₁₂ = 30.6 m, measured for 18 × 275; the other rows
carry that lever arm for want of the per-optics value, so read them as
scaling. Very small k_rel puts both fragments near the beam and may lose
both — quantifying that tail is the study, not this table.)

**An intact ⁶Li is one hit; the breakup is two, tens of pixels apart, in
sensors that already exist.** Topology is a stronger discriminant than
dE/dx for exactly the background #19 was written about, and nobody in
this programme had looked at it.

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

## 9.4 A baseline correction: the ePIC pot geometry has moved

Independent of nanowires, and the most immediately actionable finding
here. The aperture `tools/fullsim` measured came from the September-2024
`epic-main` inside `jug_xl-nightly`. Reading the current `main` of
`eic/epic` directly (2026-08-26):

| | September-2024 (what was measured) | current `main` |
|---|---|---|
| module | 32 × 32 mm | **16 × 16 mm** |
| insertion | energy-*independent*: 3.2 cm outer, +0.7 cm central | **per-energy 10σ offsets** in `beamline_*.xml` |
| RF shield | 1 mm aluminium per module face, active | **commented out** — *"we don't know if we will even need it … Oct. 2025"* |
| implied blind block at 18 × 275 | 32 mm (x) × 7 mm (y) | 16 mm (x) × 2.7 mm (y) |

The old block gives 32 mm / 30.6 m = 1.046 mrad against the 1.03 mrad the
scan measured — agreement to 1.5%, which is the check that the measurement
and the file reading confirm each other. **The correction runs in opposite directions at the two ends**, which is
why it cannot be guessed. At 18 × 275 the 16 mm block implies ≈0.52 mrad
— *below* the 0.727 mrad 10σ envelope, so the sensor package no longer
binds and the beam does. At 5 × 41 the per-energy insertion moves the
other way, to a 29.6 mm inner edge, because the pots now retract properly
for the larger low-energy beam: the aperture there gets **worse**. That
is precisely the configuration at which the coherent programme was said
to survive.

**So the measured aperture, the ⟨cos 2β⟩ = −0.77, the 1.4×10⁻² tagged
fraction and "the coherent programme is a low-energy programme" are all
conditioned on a superseded geometry** and must be re-measured before
being quoted further. The *curves* of §9.0 are unaffected — they price
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

### B1 — re-measure the pot aperture on the current `epic-main` ☐ **priority**
§9.4. `tools/fullsim/ion_gun_hepmc.py` already shoots an intact ⁶Li;
what is needed is the current container (`eic_xl-nightly`, not on the
analysis box as of 2026-08-26) and a repeat per optics, with the RF
shields as they now stand. Everything aperture-conditional depends on it.

### B2 — the α-tag and inclusive-tagging gain, through the chain ☐
§9.0 quotes the α-tag acceptance analytically. Re-run
`money_tagged_azz.py` and `tagged_polarimetry_7li.py` at the near-beam
aperture to say what it buys those, not only the coherent channel.

### B3 — ⁷Li ☐
Every number here is ⁶Li. ⁷Li's α tag already lands inside the Roman-Pot
window (plans/00, strategic finding 2), so the gain is smaller and needs
its own pass before any claim is made.

### B4 — two-hit topology for α + d ☐ **new, and cheap**
§9.2. What fraction of breakup events puts both fragments in acceptance,
and how well does hit multiplicity separate them from an intact ⁶Li?
Free, already instrumented, unexamined.

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

### D3 — per-optics R₁₂ and R₃₄ · **ePIC FF WG**
R₁₂ = 30.6 m was measured for 18 × 275 only, which is why this plan works
in angle and quotes millimetres for that optics alone. The 10σ offsets in
`beamline_*.xml` exist only for proton optics and are marked *"rough
extrapolation"* at 5 × 41.

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
  planes, giving ≈4.8σ per plane against a nanowire's one bit, so a
  nanowire is a downgrade. **Ask the incumbent first (D1), and ask C-AD
  what the 10σ rule actually protects (D2).** Ranked list of what would
  help the lithium tags: (1) the low-energy configuration; (2) re-tiling
  the near-beam silicon, which ePIC is already designing; (3) using the
  charge EICROC already digitises; (4) two-hit topology for α + d; (5)
  the IR-8 secondary focus and its z² Cherenkov. Superconducting
  nanowires appear nowhere on it.
* **Where the technology *would* belong.** The Argonne programme's own
  four-application list points at the cold bore of a superconducting
  magnet, in front of the ZDC, and a Compton polarimeter. The result that
  makes the technology unique — saturated efficiency to 5 T *parallel* to
  the device plane, against 0.5 T perpendicular — is worth nothing at a
  pot in a field-free drift and everything inside a magnet.
