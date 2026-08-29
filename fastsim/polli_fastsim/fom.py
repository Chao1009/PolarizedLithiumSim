"""Figure-of-merit projections: rates per (x,Q2) bin -> delta(observable).

Luminosity convention: we quote integrated luminosity PER NUCLEON
L_eN [fb^-1] (EIC convention for eA: L_eA * A). The per-nucleon DIS cross
section is built from F2A/A, so N_events(bin) = L_eN * sigma(bin, F2A/A).

Default scenario: 10 fb^-1/nucleon at each energy point with electron
polarization 0.7. Ion polarizations default to the ECRP source targets
(Pz >= 0.90, Pzz >= 0.80) degraded to in-ring placeholders (0.7/0.6) --
ring transport survival is an open accelerator-physics question tracked
in plans/03_open_questions.md.

Run-plan share: `lumi_fb_per_nucleon` is the PROGRAMME luminosity ("one
EIC year" = 10 fb^-1/u) and `run_share` is the fraction of it this
observable is given.  The two are kept apart, rather than folded into one
number, so that a projection can say "10 fb^-1/u of programme, 1/3 of it
here": the programme figure is the machine's, the share is ours to
propose (plans/07 WP2).  Every published number is at run_share = 1.

Three shares stack multiplicatively and are deliberately NOT the same
object:

  run_share              the programme share -- how a year divides between
                         observables, isotopes and running configurations
                         (this field).
  Optics.lumi_fraction   the optics penalty of a running configuration --
                         what a de-squeezed beta*_x costs at fixed wall
                         time (farforward.Optics; 1/7 - 1/13 for the
                         lithium tagging optics).
  SpinCategory.lumi_fraction
                         the INTRA-observable spin share -- how one
                         measurement's own luminosity divides between its
                         fill states (evgen.polligen.bookkeeping; 0.5/0.5
                         for a tensor flip plan).  It sums to one within a
                         run plan and is not a programme decision.
"""

from dataclasses import dataclass, field

import numpy as np

from . import kinematics as kin
from .asymmetries import (
    a_cos2phi,
    a_parallel,
    azz,
    err_a_parallel,
    err_azz,
    err_cos2phi_amplitude,
    depolarization_d,
)
from .structure import NuclearF2, dsigma_dx_dq2


@dataclass
class Scenario:
    lumi_fb_per_nucleon: float = 10.0
    # fraction of the programme luminosity given to THIS observable in
    # THIS configuration (see the module docstring).  Statistical errors
    # scale as 1/sqrt(run_share); any luminosity at which a target
    # significance is reached is invariant under it.
    run_share: float = 1.0
    pol_electron: float = 0.70
    pol_ion_vector: float = 0.70   # P_z in the ring (placeholder)
    pol_ion_tensor: float = 0.60   # P_zz in the ring (placeholder)
    q2_min: float = 1.0
    y_min: float = 0.01
    y_max: float = 0.95
    w2_min: float = 10.0
    # crude central-detector acceptance for the scattered electron
    eta_min: float = -3.5
    eta_max: float = 3.5
    e_prime_min: float = 0.5  # GeV

    def __post_init__(self):
        if not self.run_share > 0:
            raise ValueError("run_share must be positive (it multiplies the "
                             "luminosity); got %r" % (self.run_share,))

    @property
    def lumi_effective_fb_per_nucleon(self):
        """Luminosity this observable actually receives [fb^-1/nucleon]."""
        return self.lumi_fb_per_nucleon * self.run_share


def run_share_tag(run_share):
    """Filename key for a non-default run-plan share ('' at share 1).

    Published figures are made at run_share = 1; a run at any other share
    appends this key so it cannot overwrite them (the same guard as
    `money_tagged_azz.output_stem`)."""
    if abs(float(run_share) - 1.0) < 1e-12:
        return ""
    return ("share%g" % float(run_share)).replace(".", "p")


def run_share_header(lumi_fb_per_nucleon, run_share):
    """One-line description of the programme luminosity and this
    observable's share of it, for a script header."""
    return ("run share %g of the %g fb^-1/u programme -> %g fb^-1/u "
            "delivered (published numbers are at share 1; errors scale as "
            "1/sqrt(share), any luminosity to a target significance does "
            "not)" % (run_share, lumi_fb_per_nucleon,
                      lumi_fb_per_nucleon * run_share))


@dataclass
class BinnedProjection:
    x_edges: np.ndarray
    q2_edges: np.ndarray
    x: np.ndarray          # 2D centers (nx, nq2)
    q2: np.ndarray
    accepted: np.ndarray   # bool mask
    n_events: np.ndarray
    extras: dict = field(default_factory=dict)


def project_rates(config, scenario, nx=40, nq2=30, x_range=(1e-4, 1.0),
                  q2_range=(1.0, 2e3), nuclear_f2=None):
    """Event counts per log-log (x,Q2) bin for a BeamConfig + Scenario."""
    s = config.sqrt_s_per_nucleon**2
    x_edges, q2_edges, x_c, q2_c = kin.log_grid(x_range, q2_range, nx, nq2)
    X, Q2 = np.meshgrid(x_c, q2_c, indexing="ij")

    mask = kin.kinematic_mask(
        X, Q2, s,
        q2_min=scenario.q2_min, y_min=scenario.y_min,
        y_max=scenario.y_max, w2_min=scenario.w2_min,
    )
    y = kin.y_from_xq2(X, Q2, s)
    e_p, _theta, eta = kin.scattered_electron(X, np.clip(y, 1e-6, 1.0), s,
                                              config.electron_energy)
    mask &= (eta >= scenario.eta_min) & (eta <= scenario.eta_max)
    mask &= e_p >= scenario.e_prime_min

    nf2 = nuclear_f2 or NuclearF2(config.ion)
    f2_per_nucleon = nf2.f2a(X, Q2) / config.ion.A
    # F_L needs the same R that turned F2 into F1 in `nf2`: carrying it on
    # the structure-function object is what stops the two going out of
    # step, which is the whole point of the r_func hook (plans/08 C2).
    # getattr, not nf2.r_func, because `nuclear_f2` is an interface.
    xsec = dsigma_dx_dq2(X, Q2, s, f2_per_nucleon,
                         r_func=getattr(nf2, "r_func", None))

    dx = np.diff(x_edges)[:, None]
    dq2 = np.diff(q2_edges)[None, :]
    # the ONE luminosity line: programme luminosity x this observable's
    # share of the run plan (plans/07 WP2)
    lumi_pb = scenario.lumi_effective_fb_per_nucleon * 1e3  # fb^-1 -> pb^-1
    n_events = np.where(mask, xsec * dx * dq2 * lumi_pb, 0.0)

    return BinnedProjection(x_edges, q2_edges, X, Q2, mask, n_events,
                            extras={"y": y, "eta": eta, "s": s, "nf2": nf2})


def project_observables(config, scenario, proj, g1_model, b1_func, delta_func):
    """Attach asymmetries + statistical errors for the three observables."""
    X, Q2, N = proj.x, proj.q2, proj.n_events
    y = proj.extras["y"]
    nf2 = proj.extras["nf2"]
    f2 = nf2.f2a(X, Q2) / config.ion.A
    f1 = nf2.f1a(X, Q2) / config.ion.A
    r_func = getattr(nf2, "r_func", None)   # see project_rates

    out = {}
    # (1) polarized EMC: A_par and delta(g1A/F1A)
    g1 = g1_model.g1_nucleus(config.ion, X, Q2) / config.ion.A
    apar = a_parallel(g1, f1, y, X, Q2, r_func=r_func)
    dapar = err_a_parallel(N, scenario.pol_electron, scenario.pol_ion_vector)
    # Dividing by the MASSLESS D inverts the massless A_par, so the
    # extraction comes out HIGH by a factor (1 + gamma^2): the exact
    # A_par = D_gamma (A1 + eta A2) is (1 + gamma^2) times the massless
    # one at these y, and feeding that larger asymmetry to the massless
    # inversion returns a correspondingly larger g1A/F1A (see
    # a_parallel's docstring).
    # That leaves Delta-R high by 0.12 / 0.44 / 0.71 / 1.06 % at
    # x = 0.089 / 0.282 / 0.447 / 0.708 -- inverse-variance weighted with
    # exactly the weights below, which give the published statistical
    # 4.77 / 5.11 / 6.19 / 12.37 % at 10 fb^-1/u, so the bias is 9-40x
    # smaller -- measured by evgen/scripts/target_mass_bound.py, whose
    # block 3 reproduces both columns.  It is a bias, not a spread, so it
    # does not enter these error bars; removing it means
    # dividing by D_gamma (1 + gamma^2 kappa + ...) instead, which is a
    # change to a published pipeline and the author's call.
    d_g1f1 = dapar / depolarization_d(y, X, Q2, r_func=r_func)
    out["a_par"] = apar
    out["err_a_par"] = dapar
    out["err_g1_over_f1"] = d_g1f1

    # (2) tensor b1 via Azz (spin-1 only)
    b1 = b1_func(X, Q2, f1)
    out["azz"] = azz(b1, f1, f2, X, y)
    out["err_azz"] = err_azz(N, scenario.pol_ion_tensor)

    # (3) gluonometry: cos(2phi) amplitude from Delta
    delta = delta_func(X, Q2, f1)
    out["a_cos2phi"] = a_cos2phi(delta, f1, f2, X, y)
    out["err_a_cos2phi"] = err_cos2phi_amplitude(N, scenario.pol_ion_tensor)

    # significance maps (|asym| / stat error), zeroed outside acceptance
    for key in ("a_par", "azz", "a_cos2phi"):
        sig = np.where(proj.accepted, np.abs(out[key]) / out["err_" + key], 0.0)
        out["sig_" + key] = sig
    proj.extras.update(out)
    return out


def combine_over_q2(err_map, accepted, min_events=None, n_events=None):
    """Combine per-(x,Q2)-bin errors into a per-x error:
    delta_x = 1/sqrt(sum_q2 1/delta^2) over accepted bins."""
    err_map = np.asarray(err_map, dtype=float)
    use = np.asarray(accepted, dtype=bool)
    if min_events is not None and n_events is not None:
        use = use & (np.asarray(n_events) >= min_events)
    inv2 = np.where(use & (err_map > 0), 1.0 / err_map**2, 0.0)
    tot = inv2.sum(axis=1)
    out = np.full(tot.shape, np.inf)
    np.divide(1.0, np.sqrt(tot), out=out, where=tot > 0)
    return out
