"""Beam species and EIC energy configurations.

TWO DIFFERENT CONSTRAINTS SET AN ION ENERGY, and which one binds depends on
the configuration (plans/10; corrected 2026-08-27).

TOP ENERGY -- magnetic RIGIDITY.  The ring caps at the 275 GeV proton
rigidity, so p/nucleon <= 275 GeV * (Z/A).  This reproduces d (137.5, CDR
rounds to "135"), Au (110 GeV/u) and 3He (183 GeV/u in the current design,
CDR Sec. 5.5; the 166 GeV/u of older papers is eRHIC-era legacy), and the
Li values of the EPIOS white paper (Atoian et al., arXiv:2510.10794,
PRC 113:060501: ~138 GeV/u 6Li, ~117 GeV/u 7Li).

LOWER ENERGIES -- revolution period, i.e. SPEED.  The HSR and the ESR must
have equal revolution period and the electrons are ultrarelativistic, so
the hadron GAMMA is fixed by the ring circumference and the magnets supply
whatever rigidity that gamma needs.  Ions are therefore GAMMA-MATCHED to
the proton configuration, NOT rigidity-scaled.

  The decisive check: Yellow Report Table 10.2 lists gold at 41 GeV/u
  (gamma = 44.02) against the 41 GeV proton's gamma = 43.70 -- equal to
  0.7%, exactly the u vs m_p mass difference.  Rigidity scaling would have
  put gold at 41 * 79/197 = 16.4 GeV/u.  It does not.  The same rule
  reproduces the published 3He menu "41, and 100-183 GeV/nucleon" exactly.

So the accessible menu for an ion is the isolated gamma ~ 43.7 point plus a
continuous band from gamma ~ 106.6 up to its own rigidity cap -- for 6Li,
"41, and 99-138 GeV/u", with NOTHING in between.  Before 2026-08-27 this
module scaled the low and mid points as 41*Z/A and 100*Z/A, giving 20.5 and
50 GeV/u for 6Li; those are not machine configurations, and correcting them
roughly DOUBLES p_ion at the low configuration.
"""

from dataclasses import dataclass

PROTON_TOP_MOMENTUM = 275.0  # GeV
PROTON_MASS = 0.938272088    # GeV

#: Physical nuclear masses [GeV], the same AME2020-derived values as
#: polli_fastsim.spectator.NUCLEUS_MASS (kept here to avoid an import cycle;
#: fastsim/tests/test_farforward.py pins the two together).  Mass per
#: NUCLEON is what converts a gamma into a per-nucleon momentum, and it is
#: not A * M_U -- a 6Li is bound by 32 MeV/nucleon relative to free protons.
NUCLEUS_MASS = {
    ("p", 1, 1): 0.938272088,
    ("d", 2, 1): 1.875612942,
    ("3He", 3, 2): 2.808391607,
    ("6Li", 6, 3): 5.601518702,
    ("7Li", 7, 3): 6.533833028,
}

#: Proton beam energies of the three reference EIC configurations [GeV].
PROTON_CONFIG_ENERGIES = (41.0, 100.0, 275.0)

#: How the ring actually achieves the revolution-frequency match, and the
#: Lorentz-factor windows it reaches.  EPIOS (arXiv:2510.10794) pp. 12-13,
#: verbatim: "their revolution frequencies have to be equal.  This
#: synchronization is accomplished by applying a radial shift of up to
#: +-20 mm in the arcs, which facilitates a range of the Lorentz factor of
#: 118 < gamma < 293.  To allow for even lower ion energies, a 'Blue' arc
#: between IR12 and IR2 will be utilized as a bypass.  The average radius of
#: this arc is about 90 cm smaller than that of the corresponding 'Yellow'
#: arc, which reduces the circumference of the HSR by roughly 90 cm.  The
#: resulting circumference then corresponds to a Lorentz factor of
#: gamma = 43.5."
#:
#: The two numbers are exactly the species menu at TOP energies -- gold at
#: 110 GeV/u is gamma = 118.1 and a 275 GeV proton is gamma = 293.1 -- so the
#: window is what the ion programme needs, not a hard reachability bound.
#: The shift is quoted as "up to" +-20 mm and spending all of it covers
#: ~200 mm of circumference against the 115 mm that 118 -> 293 requires, so
#: there is headroom below 118 that EPIOS does not quantify.
#:
#: KNOWN CONFLICT: Yellow Report Table 10.1 runs 100 GeV protons, i.e.
#: gamma = 106.6, which falls in NEITHER stated window.  The 41 GeV
#: (gamma 43.7) and 275 GeV (gamma 293.1) anchors are both inside them, and
#: YR Table 10.2's gold at 41 GeV/u (gamma 44.0) sits in the bypass window
#: with the 41 GeV proton -- which is the independent confirmation that ions
#: are gamma-matched.  This module stays anchored on the Yellow Report's
#: three configurations because every beam parameter table is indexed by
#: them; `tools/consistency_check.py` flags the 100 GeV point rather than
#: silently choosing between the two documents.
EPIOS_GAMMA_BYPASS = 43.5              # "Blue" arc bypass, -90 cm
EPIOS_GAMMA_SHIFT_RANGE = (118.0, 293.0)   # +-20 mm radial shift in the arcs


def gamma_of(proton_energy):
    """Lorentz factor of the ring at the configuration whose proton energy
    is `proton_energy` [GeV].  Species-independent: it is a property of the
    revolution frequency, which is why ions are gamma-matched."""
    return ((proton_energy ** 2 + PROTON_MASS ** 2) ** 0.5) / PROTON_MASS


def epios_window_of(proton_energy, shift_gamma=2.7, tol=0.005):
    """Which EPIOS synchronisation window a configuration falls in:
    "bypass", "radial-shift", or None if neither (see the note above).

    `tol` is a relative tolerance on the window edges, because EPIOS quotes
    them rounded: the 275 GeV proton is gamma = 293.1 against a stated 293,
    and gold's 110 GeV/u top is 118.1 against a stated 118 -- both are the
    endpoints the window was written to describe."""
    g = gamma_of(proton_energy)
    lo, hi = EPIOS_GAMMA_SHIFT_RANGE
    if abs(g - EPIOS_GAMMA_BYPASS) <= shift_gamma:
        return "bypass"
    if lo * (1.0 - tol) <= g <= hi * (1.0 + tol):
        return "radial-shift"
    return None


# --- the 6Li cluster wave function ---------------------------------------
#
# ONE SOURCE OF TRUTH for the two D-state probabilities the 6Li cluster
# picture is built from.  They live here, in the module every spin
# consumer already imports, and `polligen.tagged` re-exports them under
# the same names, so the INCLUSIVE effective polarization below and the
# TAGGED S/D interference of `tagged.li6_alpha_channel` cannot drift
# apart: they are the same wave function seen in two experiments.

#: alpha-d relative D-state probability.  Chosen so that the embedded
#: deuteron's vector dilution 1 - (3/2) P_D reproduces the 0.87 of
#: `polarized.b1_li6_from_deuteron` (SCENARIO -- VMC overlaps are the
#: scheduled replacement, plans/04 #15).
P_D_LI6 = 0.0867
#: the deuteron's own D-state probability (AV18-like).
P_D_DEUTERON = 0.045

#: Vector depolarization 1 - (3/2) P_D of a spin-1 system with D-state
#: probability P_D.  The deuteron slot below carries the second of these
#: verbatim, so the two ions are built from one expression and their
#: ratio is exact rather than rounded.
ALPHA_D_VECTOR_POLARIZATION = 1.0 - 1.5 * P_D_LI6        # 0.86995
DEUTERON_VECTOR_POLARIZATION = 1.0 - 1.5 * P_D_DEUTERON  # 0.9325

#: WHOLE-NUCLEUS vector polarization of the two polarized nucleons of
#: 6Li in the cluster picture (author decision 2026-08-29, plans/04 #6).
#: The 6Li spin is carried by the alpha-d relative motion and by the
#: deuteron inside it, so a nucleon of that deuteron is polarized along
#: the 6Li spin by the PRODUCT of the two dilutions -- 0.86995 x 0.9325 =
#: 0.81123.  The alpha contributes nothing (J = 0).
LI6_CLUSTER_POLARIZATION = (ALPHA_D_VECTOR_POLARIZATION
                            * DEUTERON_VECTOR_POLARIZATION)

#: The retired alternative, kept reachable and pinned: Cloet's slides use
#: P_p = P_n = 1/3, i.e. a whole-nucleus Z*P_p = N*P_n = 1 -- one fully
#: polarized proton and neutron out of three each.  It was the default
#: before 2026-08-29 and is 1/0.81123 = 1.233 times the cluster value.
LI6_NAIVE_ONE_THIRD = 1.0 / 3.0


@dataclass(frozen=True)
class Ion:
    name: str
    A: int
    Z: int
    spin: float
    # PER-NUCLEON effective polarizations P_p, P_n, from which
    #   g1A = Z*P_p*g1p + N*P_n*g1n   (naive denominator of the EMC ratio DR)
    # exactly as NuclearF2.f2a builds F2A from Z*f2p + N*f2n.  The nucleon
    # counts moved into ToyG1.g1_nucleus on 2026-08-28 (plans/08 D7); the
    # two slots that held WHOLE-NUCLEUS sums were rescaled in the same
    # change so that every published number is bit-for-bit unchanged.
    # 7Li: the verified whole-nucleus sums are P_p=+0.866, P_n=-0.037
    #      (QMC/VMC, Wiringa 1309.3794 via JLab E12-14-001); stored here
    #      DIVIDED BY Z=3 and N=4, so Z*P_p and N*P_n return them exactly
    #      (pinned in fastsim/tests/test_polarized_normalisation.py).
    # 3He: Bissey PRC 65:064317 -- verified, and per-nucleon already, so
    #      the proton term now carries the factor Z=2 the literature
    #      intends and that the old whole-nucleus reading dropped.
    # 6Li: the CLUSTER PICTURE, since the author decision of 2026-08-29
    #      that closed plans/04_open_questions.md item 6.  The 6Li spin
    #      is carried by the deuteron cluster, so a nucleon of that
    #      deuteron is polarized along the 6Li spin by the product of
    #      the alpha-d and deuteron vector dilutions,
    #      LI6_CLUSTER_POLARIZATION = 0.86995 x 0.9325 = 0.81123
    #      whole-nucleus (Schellingerhout PRC 48:2714), and the slots
    #      hold a third of it each so that Z*P_p = N*P_n = 0.81123.
    #      Built from the SAME two D-state probabilities the tagged
    #      sector uses -- P_D_LI6 and P_D_DEUTERON above, re-exported by
    #      polligen.tagged -- so the inclusive and tagged 6Li share one
    #      wave function, and from the same expression the deuteron slot
    #      carries, so the ratio below is exact.  NOT adopted, and the
    #      upper end of the band: the six-body VMC of Wiringa PRC
    #      89:024305 Table I -- the table LI7's slots come from -- counts
    #      1.924 spin-up against 1.076 spin-down protons (and neutrons)
    #      in the M = 1 state, i.e. the same whole-nucleus quantity read
    #      ab initio as 0.848, 4.5% above this one, and implying an
    #      alpha-d factor 0.848/0.9325 = 0.909 in place of E155's 0.870.
    #      The cluster construction is adopted over it because it and the
    #      tagged sector are then one wave function rather than two
    #      transcribed constants; 0.81-0.85 is the band (plans/02 step
    #      1.1 item 2, plans/04 #6, and #15 for the VMC alpha-d overlap
    #      that would replace the scenario P_D_LI6 altogether).  The
    #      retired alternative, Cloet's per-nucleon 1/3, survives as
    #      LI6_NAIVE_ONE_THIRD: it is a whole-nucleus 1.0, i.e. 1.233
    #      times this one and above the band.
    #      Per-nucleon g1(6Li)/g1(d) is (1-1.5*P_D_LI6)/3 = 0.290 -- the
    #      deuteron's own dilution cancels between the two isoscalar
    #      ions -- against the 0.358 the naive constant gave.
    #      WHAT MOVED WITH IT, measured 2026-08-29: g1(6Li) is
    #      multiplied by 0.81123 and no published number in the
    #      repository changes.  `fom.project_observables`' err_azz and
    #      err_g1_over_f1 are counting errors and carry no eff_pol;
    #      `phase_space_map.py` defaults to 7Li; and every published
    #      cos 2phi / Delta figure runs the transverse categories of
    #      `bookkeeping.transverse_tensor_plan` at theta_S = pi/2 with
    #      unpolarized electrons, where the g1 term enters only through
    #      cos(theta_S) = 0 and through a_1's lam_e = 0 -- w_avg and a_2
    #      are bit-for-bit unchanged there.  The one place it does move
    #      is the longitudinal vector-L term of `xsec.InclusiveKernel`,
    #      w_avg at theta_S = 0, m = 1, lam_e P_e = 0.7 with no b1, which
    #      is exactly proportional to the constant: at 10 x 99.5 GeV/u it
    #      is -0.000247 / 0.000010 / 0.000129 / 0.000363 / 0.000500 /
    #      0.000862 against the naive constant's -0.000304 / 0.000012 /
    #      0.000159 / 0.000448 / 0.000616 / 0.001063, over the six
    #      accepted grid cells nearest (x, Q2) = (0.005, 1.1) ...
    #      (0.28, 25) -- the ratio is 0.811228 in every one of them (the
    #      probe is in docs/reproduction_manual.md).  That is the
    #      closure_fom.py A_par panel, whose estimator variance is
    #      1/(P_e P_z)^2 N up to O(A_par^2) <= 2e-5 -- below its own
    #      Monte-Carlo band.
    eff_pol_p: float = 0.0
    eff_pol_n: float = 0.0

    @property
    def N(self) -> int:
        return self.A - self.Z

    @property
    def mass_per_nucleon(self) -> float:
        """Physical nuclear mass / A [GeV].  Not M_U: binding matters at the
        0.5% level, which is what separates a 41 GeV proton (0.9383) from a
        gamma-matched 6Li at 40.8 GeV/u (0.9336)."""
        return NUCLEUS_MASS[(self.name, self.A, self.Z)] / self.A

    @property
    def momentum_per_nucleon_max(self) -> float:
        """Rigidity-limited top momentum per nucleon [GeV]."""
        return PROTON_TOP_MOMENTUM * self.Z / self.A

    def momentum_per_nucleon_at(self, proton_energy: float) -> float:
        """Per-nucleon momentum at the configuration whose proton energy is
        `proton_energy` [GeV], i.e. GAMMA-MATCHED to that proton and then
        capped by the ring's rigidity limit (plans/10).

        Gamma-matching is what the equal-revolution-period constraint
        imposes; the rigidity cap is what the magnets impose.  Whichever
        binds, binds -- which is why the top configuration comes out
        rigidity-limited and the lower ones do not."""
        gamma = (proton_energy ** 2 + PROTON_MASS ** 2) ** 0.5 / PROTON_MASS
        p_u = self.mass_per_nucleon * (gamma ** 2 - 1.0) ** 0.5
        return min(p_u, self.momentum_per_nucleon_max)


PROTON = Ion("p", 1, 1, 0.5, eff_pol_p=1.0, eff_pol_n=0.0)
# the deuteron's slot IS the expression 6Li's is built from, so the two
# agree bit for bit and the ratio of their g1 is (1-1.5*P_D_LI6)/3 exactly
DEUTERON = Ion("d", 2, 1, 1.0, eff_pol_p=DEUTERON_VECTOR_POLARIZATION,
               eff_pol_n=DEUTERON_VECTOR_POLARIZATION)
HE3 = Ion("3He", 3, 2, 0.5, eff_pol_p=-0.028, eff_pol_n=0.86)
# LI6_CLUSTER_POLARIZATION/3 per slot: whole-nucleus 0.81123 (2026-08-29)
LI6 = Ion("6Li", 6, 3, 1.0, eff_pol_p=LI6_CLUSTER_POLARIZATION / 3.0,
          eff_pol_n=LI6_CLUSTER_POLARIZATION / 3.0)
# 0.866/3 and -0.037/4: the whole-nucleus VMC sums per nucleon (D7).
LI7 = Ion("7Li", 7, 3, 1.5, eff_pol_p=0.866 / 3.0, eff_pol_n=-0.037 / 4.0)

IONS = {i.name: i for i in (PROTON, DEUTERON, HE3, LI6, LI7)}

ELECTRON_ENERGIES = (5.0, 10.0, 18.0)  # GeV


@dataclass(frozen=True)
class BeamConfig:
    electron_energy: float  # GeV
    ion: Ion
    ion_momentum_per_nucleon: float  # GeV

    @property
    def sqrt_s_per_nucleon(self) -> float:
        """sqrt(s) of the electron-nucleon system [GeV] (massless approx)."""
        return (4.0 * self.electron_energy * self.ion_momentum_per_nucleon) ** 0.5

    def label(self) -> str:
        return (
            f"e({self.electron_energy:g}) x {self.ion.name}"
            f"({self.ion_momentum_per_nucleon:g}/u)"
        )


def default_configs(ion_name: str = "7Li") -> list:
    """Reference energy scan: low/mid/top, mirroring ep 5x41, 10x100, 18x275.

    Each point is GAMMA-MATCHED to its proton configuration and then capped
    by the ring rigidity (see the module docstring).  For 6Li that gives
    40.8 / 99.5 / 137.5 GeV/u -- the top one rigidity-limited, the other two
    not.  It gave 20.5 / 50 / 137.5 before 2026-08-27, which put the two
    lower points at energies the machine cannot deliver."""
    ion = IONS[ion_name]
    return [
        BeamConfig(e, ion, round(ion.momentum_per_nucleon_at(p), 1))
        for e, p in zip(ELECTRON_ENERGIES, PROTON_CONFIG_ENERGIES)
    ]
