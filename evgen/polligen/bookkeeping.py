"""Run-plan bookkeeping: spin categories, bunch patterns, luminosity shares.

A run plan is a list of SpinCategory: each combines an ion fill state
(populations along an axis) with an electron helicity sign and carries its
absolute share of the integrated luminosity, including relative-luminosity
offsets.  Every generated event is labeled by its category, and the
analysis-side estimators (estimators.py) see exactly the per-category
luminosities the "experiment" would measure.

This module also owns the two Step-5.A systematics knobs (plans/05 §5.0):

* relative-luminosity offsets between spin states.  Analytic first-order
  biases of the naive (equal-share) estimators, validated by the
  pseudo-experiment tests at large offset:
      tensor thirds:   bias(Azz)   = -(2/3) * delta / P_zz
                       (offset delta on the m0-enriched share)
      helicity flips:  bias(A_par) = +delta / (2 P_e P_z)
                       (offset delta on the lam_e=+1 share)
  At the plans/05 reference delta = 1e-4 with in-ring (P_e, P_z, P_zz) =
  (0.7, 0.7, 0.6): |bias| ~ 1.1e-4 on Azz and ~1.0e-4 on A_par -- below
  the ~1e-3 per-bin statistical floor at 10 fb^-1/u, but not negligible
  against combined bins; the lumi-corrected estimators remove it.

* polarimetry scale: measured polarizations = true * (1 + delta_p * gauss),
  drawn once per fill from a fixed-seed stream (delta P/P ~ 3%, plans/04 #5).
"""

from dataclasses import dataclass, field, replace
from typing import Tuple

import numpy as np

from . import spin as spinmod


@dataclass(frozen=True)
class SpinCategory:
    """One (ion fill state, electron helicity) luminosity category."""
    name: str
    j: float
    populations: Tuple[float, ...]   # m = +J ... -J along the axis
    lam_e: int = 0                   # +1/-1 electron helicity; 0 = unpol
    pe: float = 0.0                  # electron polarization magnitude
    theta_s: float = 0.0             # quantization-axis polar angle [rad]
    phi_s: float = 0.0               # quantization-axis azimuth [rad]
    lumi_fraction: float = 1.0       # absolute share of the total lumi

    def moments(self):
        return spinmod.moments_along_axis(self.j, self.populations)


@dataclass
class RunPlan:
    """Spin categories + true/measured polarization bookkeeping."""
    categories: list
    pe_true: float = 0.0
    pz_true: float = 0.0
    pzz_true: float = 0.0
    # polarimetry: measured = true * (1 + delta_p_over_p * gauss)
    delta_p_over_p: float = 0.0
    polarimetry_seed: int = 20260713
    measured: dict = field(default_factory=dict)

    def __post_init__(self):
        tot = sum(c.lumi_fraction for c in self.categories)
        if tot <= 0:
            raise ValueError("run plan has no luminosity")
        rng = np.random.default_rng(self.polarimetry_seed)
        for key, true in (("pe", self.pe_true), ("pz", self.pz_true),
                          ("pzz", self.pzz_true)):
            smear = 1.0 + self.delta_p_over_p * rng.standard_normal()
            self.measured[key] = true * smear if self.delta_p_over_p else true

    def lumi_shares(self, total_lumi_pb):
        """{category name: integrated lumi [pb^-1]} (absolute fractions)."""
        return {c.name: total_lumi_pb * c.lumi_fraction
                for c in self.categories}


def helicity_flip_plan(j, pz, pe, pzz=None, theta_s=0.0, phi_s=0.0,
                       rel_lumi_offset=0.0, name="apar"):
    """A_par run plan: one vector fill, electron helicity flipped +-.

    `pzz=None` (default) uses spin-temperature (max-entropy) populations
    for the requested pz -- the physical vector-fill model, which for
    j >= 1 implies a non-zero tensor moment (recorded in pzz_true).
    Passing pzz selects the explicit (pz, pzz) [spin-1] or (pz, t)
    [spin-3/2] populations instead.  `rel_lumi_offset` boosts the
    lam_e = +1 share only: L+ = L/2 (1+offset), L- = L/2 -- the
    convention of apar_rel_lumi_bias.
    """
    if pzz is None:
        pops = spinmod.populations_maxent(j, pz)
    elif abs(j - 1.0) < 1e-9:
        pops = spinmod.spin1_populations(pz, pzz)
    elif abs(j - 1.5) < 1e-9:
        pops = spinmod.spin32_populations(pz, pzz)
    elif abs(j - 0.5) < 1e-9:
        pops = ((1.0 + pz) / 2.0, (1.0 - pz) / 2.0)
    else:
        raise ValueError("unsupported spin %g" % j)
    lp, lm = 0.5 * (1.0 + rel_lumi_offset), 0.5
    cats = [
        SpinCategory(name + "+", j, pops, lam_e=+1, pe=pe, theta_s=theta_s,
                     phi_s=phi_s, lumi_fraction=lp),
        SpinCategory(name + "-", j, pops, lam_e=-1, pe=pe, theta_s=theta_s,
                     phi_s=phi_s, lumi_fraction=lm),
    ]
    pzz_true = (float(spinmod.moments_along_axis(j, pops)[1])
                if j >= 1.0 - 1e-9 else 0.0)
    return RunPlan(cats, pe_true=pe, pz_true=pz, pzz_true=pzz_true)


def tensor_thirds_plan(pz, pzz, rel_lumi_offset=0.0, theta_s=0.0, phi_s=0.0,
                       name="azz"):
    """Spin-1 A_zz run plan: equal-thirds fills (+, -, 0-enriched).

    Fills: (pz, +pzz), (-pz, +pzz) and the m0-enriched (0, -2 pzz) --
    the HERMES-style pattern that makes the canonical thirds estimator
    (y+ + y- - 2 y0)/(y+ + y- + y0)/pzz exact.  Electron unpolarized.
    `rel_lumi_offset` boosts the 0-enriched share only:
    L0 = L/3 (1+offset), L+- = L/3 -- the convention of azz_rel_lumi_bias.
    """
    pops_p = spinmod.spin1_populations(pz, pzz)
    pops_m = spinmod.spin1_populations(-pz, pzz)
    pops_0 = spinmod.spin1_populations(0.0, -2.0 * pzz)
    f0 = (1.0 + rel_lumi_offset) / 3.0
    fpm = 1.0 / 3.0
    cats = [
        SpinCategory(name + "+", 1.0, pops_p, theta_s=theta_s, phi_s=phi_s,
                     lumi_fraction=fpm),
        SpinCategory(name + "-", 1.0, pops_m, theta_s=theta_s, phi_s=phi_s,
                     lumi_fraction=fpm),
        SpinCategory(name + "0", 1.0, pops_0, theta_s=theta_s, phi_s=phi_s,
                     lumi_fraction=f0),
    ]
    return RunPlan(cats, pz_true=pz, pzz_true=pzz)


def transverse_tensor_plan(pzz, phi_s=0.0, name="cos2phi"):
    """Gluonometry run plan: tensor-polarized fill, axis in the transverse
    plane (theta_S = 90 deg), electron unpolarized.  Single category; the
    cos(2 phi') amplitude is extracted from the phi distribution."""
    pops = spinmod.spin1_populations(0.0, pzz)
    cats = [SpinCategory(name, 1.0, pops, theta_s=np.pi / 2.0, phi_s=phi_s,
                         lumi_fraction=1.0)]
    return RunPlan(cats, pzz_true=pzz)


def tensor_flip_plan(pzz, phi_s=0.0, share_plus=0.5, rel_lumi_offset=0.0,
                     name="flip"):
    """Gluonometry run plan with the acceptance-cancelling spin-state
    pattern (reconstruction-chain note, 2026-08-24): m = +-1-rich bunches
    at (P_z, P_zz) = (0, +pzz) and m = 0-rich bunches at (0, -2 pzz) --
    the same source purity, as in tensor_thirds_plan -- both with the
    alignment axis transverse (theta_S = 90 deg) at azimuth phi_s and
    unpolarized electrons.  `share_plus` is the luminosity share of the
    +pzz fill; `rel_lumi_offset` boosts that share by (1+offset) without
    the analysis knowing (the bookkeeper's convention).  The analysis
    estimator is reco.harmonic_ratio_fit."""
    pops_p = spinmod.spin1_populations(0.0, pzz)
    pops_0 = spinmod.spin1_populations(0.0, -2.0 * pzz)
    cats = [
        SpinCategory(name + "+", 1.0, pops_p, theta_s=np.pi / 2.0,
                     phi_s=phi_s,
                     lumi_fraction=share_plus * (1.0 + rel_lumi_offset)),
        SpinCategory(name + "0", 1.0, pops_0, theta_s=np.pi / 2.0,
                     phi_s=phi_s, lumi_fraction=1.0 - share_plus),
    ]
    return RunPlan(cats, pzz_true=pzz)


def with_offset(plan, category_name, offset):
    """Copy of `plan` with one category's lumi share scaled by (1+offset),
    other categories untouched (the convention of the bias formulas)."""
    if category_name not in {c.name for c in plan.categories}:
        raise KeyError(category_name)
    cats = [replace(c, lumi_fraction=c.lumi_fraction * (1.0 + offset))
            if c.name == category_name else c
            for c in plan.categories]
    new = RunPlan(cats, pe_true=plan.pe_true, pz_true=plan.pz_true,
                  pzz_true=plan.pzz_true,
                  delta_p_over_p=plan.delta_p_over_p,
                  polarimetry_seed=plan.polarimetry_seed)
    return new


def azz_rel_lumi_bias(offset, pzz):
    """First-order naive-thirds-estimator bias from a relative-luminosity
    offset on the m0-enriched share: bias(Azz) = -(2/3) offset / pzz."""
    return -(2.0 / 3.0) * offset / pzz


def apar_rel_lumi_bias(offset, pe, pz):
    """First-order helicity-flip bias: bias(A_par) = offset / (2 pe pz)."""
    return offset / (2.0 * pe * pz)


def bunch_rng(seed, run, bunch):
    """Fixed rng stream per (run, bunch) -- reproducibility discipline of
    plans/05 §5.2.4 (same idea as spectator.sample_k's fixed seed)."""
    return np.random.default_rng(np.random.SeedSequence((seed, run, bunch)))
