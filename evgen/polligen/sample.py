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
        # per-cell unpolarized cross section [pb]: divide out the scenario
        # lumi.  It is the EFFECTIVE luminosity that went into n_events --
        # programme luminosity x run-plan share -- and dividing by the
        # programme figure alone would leave the share multiplying a cross
        # section, which is not a quantity a run plan can move (2026-08-28).
        lumi_pb = self.scenario.lumi_effective_fb_per_nucleon * 1e3
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

    @staticmethod
    def _density_min(a1n, a2n):
        """Minimum over phi of W/(1 + w_avg) = 1 + A cos phi + B cos 2phi.

        With c = cos phi in [-1, 1] this is the quadratic
        f(c) = 2B c^2 + A c + (1 - B), so the minimum is at the vertex
        c* = -A/(4B) when B > 0 and |c*| <= 1, and at an endpoint
        otherwise -- exact, where the accept-reject envelope
        1 + |A| + |B| is only a bound."""
        a, b = np.asarray(a1n, dtype=float), np.asarray(a2n, dtype=float)
        ends = 1.0 - np.abs(a) + b                    # min(f(+1), f(-1))
        with np.errstate(divide="ignore", invalid="ignore"):
            cstar = np.where(b > 0, -a / (4.0 * b), np.inf)
        vertex = (1.0 - b) - a * a / np.where(b > 0, 8.0 * b, 1.0)
        return np.where((b > 0) & (np.abs(cstar) <= 1.0), vertex, ends)

    def _check_positive(self, w_avg, a1, a2, label):
        """Raise if the phi density goes negative anywhere on the grid.

        The accept-reject of `_sample_state` draws u in [0, bound) and
        accepts on u < density, so a negative density is silently sampled
        as max(W, 0): the modulation comes out diluted AND the (x, Q2)
        mixture is skewed, because the cell weights assume the per-cell
        phi integral is 2 pi (1 + w_avg)."""
        den = 1.0 + w_avg
        margin = self._density_min(a1 / den, a2 / den)
        bad = margin < 0.0
        if np.any(bad):
            k = int(np.argmin(margin))
            raise ValueError(
                "negative phi density for %s: %d of %d cells, worst "
                "min(W)/(1+w_avg) = %.4f at x = %.4g, Q2 = %.4g "
                "(a1/den = %.3f, a2/den = %.3f).  The sampler would "
                "silently draw max(W, 0), diluting the modulation and "
                "skewing the (x, Q2) mixture; reduce the scenario "
                "amplitude." % (label, int(bad.sum()), margin.size,
                                float(margin[k]), float(self.x_cells[k]),
                                float(self.q2_cells[k]), float(a1[k] / den[k]),
                                float(a2[k] / den[k])))
        return float(np.min(margin))

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
            self._check_positive(w_avg, a1, a2,
                                 "%s m=%g" % (category.name, m))
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

    def effective_modulation(self, category, mask=None):
        """Fill-averaged accepted cross section and phi amplitudes.

        Returns (sigma_pb, a1_eff, a2_eff) for the category's population
        mixture over its accepted cells (optionally restricted to a
        boolean `mask` over the accepted-cell arrays, e.g. a sweet-spot
        super-bin in (x, Q2)):
            sigma_pb = sum_m p_m sum_c sigma_c (1 + w_avg)
            a_k_eff  = [sum_m p_m sum_c sigma_c a_k] / sigma_pb
        so that the category's phi' distribution is
        proportional to 1 + a1_eff cos phi' + a2_eff cos 2phi'.
        For a tensor-only transverse fill (b1 = 0), a2_eff reduces to
        P_zz * asymmetries.a_cos2phi rate-averaged over the cells
        (tested in test_coherent.py).
        """
        ms = _m_values(category.j)
        pops = np.asarray(category.populations, dtype=float)
        sel = (slice(None) if mask is None
               else np.asarray(mask, dtype=bool))
        xsec = self.xsec_flat[sel]
        den = 0.0
        num1 = 0.0
        num2 = 0.0
        for p_m, mm in zip(pops, ms):
            if p_m <= 0.0:
                continue
            w_avg, a1, a2 = self._amplitudes(category, mm)
            den += p_m * float((xsec * (1.0 + w_avg[sel])).sum())
            num1 += p_m * float((xsec * a1[sel]).sum())
            num2 += p_m * float((xsec * a2[sel]).sum())
        if den <= 0.0:
            return 0.0, 0.0, 0.0
        return den, num1 / den, num2 / den

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

    def _in_acceptance(self, x, q2):
        """Event-level acceptance: the same kinematic + scattered-electron
        cuts fom.project_rates applies at cell centers.  Needed because
        log-uniform in-cell placement in boundary cells would otherwise
        emit events outside the window -- including unphysical y > 1
        (2026-08-11 audit: 0.2% of events on the 60x45 grid before this
        guard)."""
        from polli_fastsim.kinematics import (kinematic_mask,
                                              scattered_electron)
        sc = self.scenario
        ok = kinematic_mask(x, q2, self.s, q2_min=sc.q2_min,
                            y_min=sc.y_min, y_max=sc.y_max,
                            w2_min=sc.w2_min)
        y = np.clip(q2 / (self.s * x), 1e-9, 1.0 - 1e-12)
        e_p, _th, eta = scattered_electron(x, y, self.s,
                                           self.config.electron_energy)
        return (ok & (eta >= sc.eta_min) & (eta <= sc.eta_max)
                & (e_p >= sc.e_prime_min))

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
        # redraw in-cell positions that fell outside the acceptance
        # (boundary cells only); fall back to the accepted cell center
        # after `tries` attempts so termination is guaranteed
        bad = ~self._in_acceptance(x, q2)
        tries = 20
        while bad.any() and tries > 0:
            nb = int(bad.sum())
            u1b = rng.uniform(size=nb)
            u2b = rng.uniform(size=nb)
            cb = cell[bad]
            x[bad] = np.exp(self.logx_lo[cb]
                            + u1b * (self.logx_hi[cb] - self.logx_lo[cb]))
            q2[bad] = np.exp(self.logq2_lo[cb]
                             + u2b * (self.logq2_hi[cb] - self.logq2_lo[cb]))
            bad = ~self._in_acceptance(x, q2)
            tries -= 1
        if bad.any():
            x[bad] = self.x_cells[cell[bad]]
            q2[bad] = self.q2_cells[cell[bad]]
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


def phi_histogram_pseudo(n_expected, a2, nbins=36, rng=None, a1=0.0,
                         poisson=True):
    """Binned phi' pseudo-experiment at full projected statistics.

    Expected counts per bin are the EXACT integral of
    n_expected/(2 pi) * (1 + a1 cos phi' + a2 cos 2phi') over each of
    `nbins` uniform bins, Poisson-fluctuated.  For binned estimators this
    carries statistics identical to event-level sampling, so 1e8-event
    projections cost 36 Poisson draws instead of 1e8 rows.

    Returns (counts, edges) ready for estimators.cos2phi_fit_binned.

    Raises if 1 + a1 cos phi' + a2 cos 2phi' goes negative anywhere: this
    path bypasses the sampler's own guard and would otherwise hand
    negative bin means to the estimator.
    """
    rng = rng or np.random.default_rng(20260713)
    margin = float(InclusiveSampler._density_min(a1, a2))
    if margin < 0.0:
        raise ValueError("negative phi' density: min(1 + %.4g cos phi' + "
                         "%.4g cos 2phi') = %.4f" % (a1, a2, margin))
    edges = np.linspace(0.0, 2.0 * np.pi, nbins + 1)
    lo, hi = edges[:-1], edges[1:]
    mu = (n_expected / (2.0 * np.pi)) * (
        (hi - lo)
        + a1 * (np.sin(hi) - np.sin(lo))
        + 0.5 * a2 * (np.sin(2.0 * hi) - np.sin(2.0 * lo)))
    counts = rng.poisson(mu) if poisson else mu
    return counts, edges


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
