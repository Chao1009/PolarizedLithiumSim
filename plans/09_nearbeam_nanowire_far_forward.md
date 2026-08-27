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

## 9.1 The idea is already on record — and its stated mechanism is wrong

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
does not relax it. The honest version of the argument is **granularity and
dead edge**: the insertion is quantised by the module, the innermost step
is held back to protect a block one module wide, and a millimetre-scale
tile whose active area reaches within a few hundred nanometres of the
substrate edge is the object that de-quantises that staircase. That is a
pot-mechanics argument worth making — but it is not the one the proposal
makes.

## 9.2 Charge identification is where the technology delivers

Open question #19 has no answer at IP6. A nanowire supplies one — but not
the way one would guess.

**Not by pulse height.** The device latches: the amplitude is the diverted
*bias* current. Both Fermilab/JPL beam papers show amplitude distributions
from 120 GeV hadrons, 120 GeV muons, 8 GeV pions and 8 GeV showering
electrons to be indistinguishable, with one threshold serving all
(arXiv:2510.11725, 2410.00251). Any pulse-height scheme is dead on arrival.

**By firing threshold.** In the normal-core hot-spot regime,

    r_s = sqrt(Q / (e π c ρ (T_c − T_0)))     I_th/I_c = 1 − 2 r_s / w

— Argonne's own Eqs. 1–2 (arXiv:2312.13405), the same relation Renema
writes as E = (w/C)²(1 − I_b/I_c)² and the ion literature as
I_th/I_c = 1 − (zeV)^½ C/w (Cristiano et al., SUST 28 (2015) 124004).
Since dE/dx ∝ z² at fixed β, and ⁶Li at 137.5 GeV/u has βγ = 148 against
128 for Argonne's calibration proton — the same velocity to 15% —
**r_s ∝ z linearly**: 134 / 268 / 402 nm for d,p / α / ⁶Li, anchored on
Argonne's *measured* 134 nm.

| wire width | I_th/I_c: p,d | α | ⁶Li | separated turn-ons |
|---|---|---|---|---|
| 250 nm | 0.00 | 0.00 | 0.00 | none — **ANL's MIP optimum** |
| 400 nm | 0.33 | 0.00 | 0.00 | none |
| 800 nm | 0.67 | 0.33 | 0.00 | α from d |
| **1000 nm** | **0.73** | **0.46** | **0.20** | **all three — the existing microwire** |
| 1500 nm | 0.82 | 0.64 | 0.46 | all three |

**Two bias points at 0.33 and 0.60 I_c tag Z by the firing pattern alone**,
both well below the 0.80 I_c at which Argonne measure the dark-count rate
rising exponentially. Bias-point charge-state discrimination is a granted
patent (US 8,872,109, Ohkubo & Suzuki), demonstrated on singly- versus
doubly-charged lysozyme by subtracting two bias points.

Three things make this the strong part of the case:

* **It is an interpolation.** Argonne have measured a 120 GeV proton at
  r_s = 134 nm and a 5.5 MeV ²⁴¹Am α that their own √Q scaling puts at
  ≈ 1 µm (arXiv:2601.03158). A relativistic ⁶Li lands at ≈ 400 nm —
  between their two measured points, on devices from the same programme.
* **The Z-ID operating point is also the blackbody-immune one.** A 300 K
  surface radiates ≈ 3×10¹⁸ photons cm⁻² s⁻¹ above the 0.04 eV SNSPD
  threshold, so a photon-mode device cannot look at a warm beam pipe.
  Biasing at 0.33 and 0.60 I_c puts the wire below its photon threshold
  and leaves only the hard-hot-spot particle mode — which is what the
  Z-ID scheme wants anyway.
* **The area answer and the Z-ID answer are the same device.** Argonne's
  efficiency optimum, w ≈ 250 nm, is precisely 2r_s for a proton: maximal
  MIP efficiency and *zero* charge information. Z-ID needs deliberately
  *wide* wires — the microwire (SMSPD) geometry Caltech/JPL/Fermilab
  already fabricate at 1 µm.

**The honest limit.** One threshold plane is a one-bit dE/dx measurement
on a straggling distribution: α → ⁶Li confusion is 20–25% per plane and
improves only weakly with film thickness, so sub-percent mis-ID needs 3–5
independent planes. That is a limitation of any threshold dE/dx detector —
and still the difference between a stack with charge information and a
Roman Pot with none.

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
and the file reading confirm each other. The current block implies roughly
half that, *below* the 0.727 mrad 10σ envelope; at 5 × 41 the per-energy
insertion moves the other way, to a 29.6 mm inner edge.

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

### D1 — the α turn-on curve · **ANL MEP**
Turn-on bias versus wire width for the ²⁴¹Am α, on the same wires as the
120 GeV proton. If the α curve sits where r_s ∝ z predicts, the Z-ID case
is validated with data they already own. Their α analysis is explicitly
*"underway"*.

### D2 — has the ~250 nm optimum been checked against anything but z = 1? · **ANL MEP**
Nobody has put a relativistic z > 1 nucleus in front of an SNSPD. If Z-ID
is real this is the first measurement of it, and FTBF or the SPS could do
it with a light-ion beam.

### D3 — per-optics R₁₂, R₃₄ and the light-ion beam sizes · **ePIC FF WG / C-AD**
R₁₂ = 30.6 m was measured for 18 × 275 only, which is why this plan works
in angle and quotes millimetres for that optics alone. The 10σ offsets in
`beamline_*.xml` exist only for proton optics, and are marked *"rough
extrapolation"* at 5 × 41. This is plans/04 #20's original ask.

### D4 — area, wire width, channel count, and by when · **ANL MEP**
Polakovic's 2024 DOE Early Career award (to 2029) is explicitly to
*"expand their effective sensing areas and interface them with
semiconducting readout electronics"*. Armstrong's parallel award targets
JLab ⁴He, not the EIC pot — so an EIC-pot demonstrator is not currently
anyone's funded near-term deliverable.

### D5 — a dose limit for NbN, and the pot radiation environment · **ANL MEP / ePIC FF WG**
Neither is published. Without both, "more radiation hard" cannot become
"therefore N mm closer".

### D6 — pot beam-coupling impedance and RF heating · **C-AD**
No quantitative EIC pot impedance study is public; the geometry comment
says the shield's necessity is unknown for exactly this reason. It sets
the cryogenic load of §9.3.

---

## 9.6 What this changes in the existing documents

* **plans/04 #19** (RP Z-ID for A/Z = 2) — a candidate answer exists and
  is sharper than expected: the firing threshold, not dE/dx amplitude.
  It rests on one unmeasured curve (D1).
* **plans/04 #20** (RP cutout geometry) — the measurement stands for the
  geometry it was made in, and that geometry has moved (§9.4).
* **The coherent programme's energy reach** — aperture-conditional, not
  physics-conditional (§9.0).
* **The recommendation.** Put the charge-identification question to the
  MEP group ahead of anything about apertures. It is specific, testable,
  answerable with data they may already own, and it addresses a gap
  nothing else at IP6 addresses.
