"""Inclusive spin-labeled event sampler (Mode G, inclusive tier).

Sampling scheme (plans/05 §5.2.4):

1. per spin category (bookkeeping.SpinCategory) and per spin projection m
   (drawn from the fill populations), the event count is Poisson with mean
   L_cat * p_m * sum_cells sigma_cell * (1 + w_avg(m)) -- the phi-averaged
   polarized modulation shifts the *rate* per spin state, which is exactly
   what counting-asymmetry estimators measure;
2. (x, Q2) by inverse-CDF over a log-log cell grid of the unpolarized
   cross section times (1 + w_avg), log-uniform inside the cell;
3. phi by accept-reject on 1 + a1 cos(phi') + a2 cos(2 phi'), phi'=phi-phi_S.

The cell grid, acceptance cuts, and per-cell rates come from
`polli_fastsim.fom.project_rates` with the kernel's own NuclearF2 -- the
pseudo-experiments are therefore normalized identically to the analytic
FOM maps they must close against (test_pseudoexp.py, scripts/closure_fom.py).

Fidelity note: modulation amplitudes are evaluated at cell centers (the
default 100x72 grid is ~4x finer than the 40x30 FOM analysis binning);
in-cell (x, Q2) placement is log-uniform.  Adequate for binned Phase-1
estimators by construction; not for unbinned in-cell shapes.
"""

import numpy as np

from polli_fastsim import fom
from polli_fastsim.kinematics import y_from_xq2

from .xsec import EventSpinState


class InclusiveSampler:
    """Grid-backed sampler for one beam configuration + kernel."""

    def __init__(self, kernel, beam_config, scenario=None, nx=100, nq2=72,
                 x_range=(1e-4, 1.0), q2_range=(1.0, 2e3), with_perp=False):
        self.kernel = kernel
        self.config = beam_config
        self.scenario = scenario or fom.Scenario()
        self.with_perp = with_perp
        proj = fom.project_rates(beam_config, self.scenario, nx=nx, nq2=nq2,
                                 x_range=x_range, q2_range=q2_range,
                                 nuclear_f2=kernel.nf2)
        self.proj = proj
        self.s = proj.extras["s"]
        # per-cell unpolarized cross section [pb]: divide out the scenario lumi
        lumi_pb = self.scenario.lumi_fb_per_nucleon * 1e3
        self.cell_xsec_pb = proj.n_events / lumi_pb
        self.accept = proj.accepted & (self.cell_xsec_pb > 0)
        self.flat_idx = np.flatnonzero(self.accept.ravel())
        self.xsec_flat = self.cell_xsec_pb.ravel()[self.flat_idx]
        self.x_cells = proj.x.ravel()[self.flat_idx]
        self.q2_cells = proj.q2.ravel()[self.flat_idx]
        # cell edges for in-cell sampling
        nxc = proj.x_edges.size - 1
        nq2c = proj.q2_edges.size - 1
        ii, jj = np.unravel_index(self.flat_idx, (nxc, nq2c))
        self.logx_lo = np.log(proj.x_edges[ii])
        self.logx_hi = np.log(proj.x_edges[ii + 1])
        self.logq2_lo = np.log(proj.q2_edges[jj])
        self.logq2_hi = np.log(proj.q2_edges[jj + 1])
        self.tables = kernel.tables(self.x_cells, self.q2_cells,
                                    with_g2=with_perp)
        self._amp_cache = {}

    # --- per-(category, m) amplitude tables --------------------------------

    def _state(self, category, m):
        return EventSpinState(lam_e=category.lam_e, pe=category.pe,
                              j=category.j, m=m, theta_s=category.theta_s,
                              phi_s=category.phi_s)

    def _amplitudes(self, category, m):
        key = (category.name, category.lam_e, category.pe, category.j, m,
               category.theta_s, category.phi_s)
        if key not in self._amp_cache:
            w_avg, a1, a2 = self.kernel.amplitudes(
                self.tables, self.x_cells, self.q2_cells, self.s,
                self._state(category, m), with_perp=self.with_perp)
            if np.any(1.0 + w_avg <= 0.0):
                raise ValueError("negative phi-averaged density for %s m=%g"
                                 % (category.name, m))
            self._amp_cache[key] = (w_avg, a1, a2)
        return self._amp_cache[key]

    # --- expected rates -----------------------------------------------------

    def sigma_tot_pb(self, category, m=None):
        """Accepted cross section [pb] for the category (or one m state)."""
        ms = np.array([m]) if m is not None else _m_values(category.j)
        pops = ([1.0] if m is not None else list(category.populations))
        tot = 0.0
        for p_m, mm in zip(pops, ms):
            if p_m <= 0.0:
                continue
            w_avg, _, _ = self._amplitudes(category, mm)
            tot += p_m * float((self.xsec_flat * (1.0 + w_avg)).sum())
        return tot

    def expected_events(self, category, lumi_pb):
        return lumi_pb * self.sigma_tot_pb(category)

    # --- sampling -------------------------------------------------------------

    def sample_category(self, category, lumi_pb=None, n=None, rng=None,
                        poisson=True):
        """Spin-labeled events for one category.

        Either `lumi_pb` (event count ~ Poisson per spin state) or a fixed
        total `n` (multinomial across spin states).  Returns dict of arrays:
        x, q2, y, phi, m, cell (flat accepted-cell index) + scalar labels.
        """
        rng = rng or np.random.default_rng(20260713)
        ms = _m_values(category.j)
        pops = np.asarray(category.populations, dtype=float)
        rates = np.array([p * (self.sigma_tot_pb(category, m=mm) if p > 0
                               else 0.0)
                          for p, mm in zip(pops, ms)])
        if lumi_pb is not None:
            counts = (rng.poisson(lumi_pb * rates) if poisson
                      else np.round(lumi_pb * rates).astype(int))
        elif n is not None:
            counts = rng.multinomial(n, rates / rates.sum())
        else:
            raise ValueError("pass lumi_pb or n")

        chunks = []
        for c_m, mm in zip(counts, ms):
            if c_m == 0:
                continue
            chunks.append(self._sample_state(category, mm, int(c_m), rng))
        if not chunks:
            empty = np.empty(0)
            return {"x": empty, "q2": empty, "y": empty, "phi": empty,
                    "m": empty, "cell": np.empty(0, dtype=int),
                    "category": category.name, "lam_e": category.lam_e}
        out = {k: np.concatenate([ch[k] for ch in chunks])
               for k in ("x", "q2", "y", "phi", "m", "cell")}
        out["cell"] = out["cell"].astype(int)
        out["category"] = category.name
        out["lam_e"] = category.lam_e
        return out

    def _sample_state(self, category, m, n, rng):
        w_avg, a1, a2 = self._amplitudes(category, m)
        weights = self.xsec_flat * (1.0 + w_avg)
        prob = weights / weights.sum()
        cell = rng.choice(prob.size, size=n, p=prob)
        u1 = rng.uniform(size=n)
        u2 = rng.uniform(size=n)
        x = np.exp(self.logx_lo[cell]
                   + u1 * (self.logx_hi[cell] - self.logx_lo[cell]))
        q2 = np.exp(self.logq2_lo[cell]
                    + u2 * (self.logq2_hi[cell] - self.logq2_lo[cell]))
        # phi accept-reject on the cell-center amplitudes
        a1n = a1[cell] / (1.0 + w_avg[cell])
        a2n = a2[cell] / (1.0 + w_avg[cell])
        bound = 1.0 + np.abs(a1n) + np.abs(a2n)
        phi_p = np.empty(n)
        todo = np.arange(n)
        while todo.size:
            cand = rng.uniform(0.0, 2.0 * np.pi, size=todo.size)
            u = rng.uniform(size=todo.size) * bound[todo]
            ok = u < (1.0 + a1n[todo] * np.cos(cand)
                      + a2n[todo] * np.cos(2.0 * cand))
            phi_p[todo[ok]] = cand[ok]
            todo = todo[~ok]
        phi = np.mod(phi_p + category.phi_s, 2.0 * np.pi)
        return {"x": x, "q2": q2, "y": y_from_xq2(x, q2, self.s),
                "phi": phi, "m": np.full(n, float(m)),
                "cell": cell}

    # --- Mode-W style weights ---------------------------------------------

    def weights_for(self, events, categories):
        """Per-event weight matrix w[i, k] = W(event_i | category_k), the
        polarized/unpolarized density ratio at the event's cell (mean over
        the fill populations).  This is the reweighting kernel Mode W
        (step 5.C) applies to external unpolarized samples."""
        cell = np.asarray(events["cell"], dtype=int)
        phi = np.asarray(events["phi"], dtype=float)
        out = np.empty((cell.size, len(categories)))
        for k, cat in enumerate(categories):
            pops = np.asarray(cat.populations, dtype=float)
            w = np.zeros(cell.size)
            for p_m, mm in zip(pops, _m_values(cat.j)):
                if p_m <= 0.0:
                    continue
                w_avg, a1, a2 = self._amplitudes(cat, mm)
                phip = phi - cat.phi_s
                w += p_m * (1.0 + w_avg[cell] + a1[cell] * np.cos(phip)
                            + a2[cell] * np.cos(2.0 * phip))
            out[:, k] = w
        return out


def _m_values(j):
    n = int(round(2 * j + 1))
    return np.array([j - i for i in range(n)])


def run_pseudo_experiment(sampler, plan, total_lumi_pb, rng=None,
                          poisson=True):
    """One pseudo-experiment: {category name: events}, {name: lumi [pb^-1]}.

    Event counts include the per-category luminosity shares (with any
    relative-luminosity offsets baked into the plan).
    """
    rng = rng or np.random.default_rng(20260713)
    shares = plan.lumi_shares(total_lumi_pb)
    events = {}
    for cat in plan.categories:
        events[cat.name] = sampler.sample_category(
            cat, lumi_pb=shares[cat.name], rng=rng, poisson=poisson)
    return events, shares
