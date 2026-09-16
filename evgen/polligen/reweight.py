"""Mode W: polarize an external, unpolarized event sample by reweighting.

Mode G (`sample.InclusiveSampler`) draws its own (x, Q2, phi) from a
cross-section grid.  Mode W does not generate kinematics at all: it takes
an *existing* unpolarized sample -- the official BeAGLE e+d EVGEN files are
the reference one -- and multiplies every event by the doubly polarized
density ratio

    W(x, Q2, phi | category) = 1 + w_avg + a_1 cos phi' + a_2 cos 2phi',
    phi' = phi - phi_S,

the same `xsec.InclusiveKernel.amplitudes` triple Mode G samples from, but
evaluated at the event's own kinematics instead of at a grid-cell centre.
Because the external sample is unpolarized, its phi is uniform and its
(x, Q2) is already distributed as the unpolarized cross section, so
resampling it with probability proportional to W reproduces the polarized
cross section exactly -- including the phi modulation, which is why the
weight must be evaluated at the event's phi and not phi-averaged.

That makes Mode W the ECCE-style pedigree demonstration of plans/05 step
5.C: the asymmetries are *injected* by the kernel into a sample whose
kinematics, hadronic final state and (later) detector response come from
outside this programme, and the analysis-side estimators of
`estimators.py` have to give them back.  The deuteron is the nucleus to do
it on, because it is the one where BeAGLE's nuclear model is right
(plans/02 step 1.5.3), and the machinery is the same one a later e+Li
BeAGLE sample would be run through for purity studies.

What is closed here and what is not
-----------------------------------
Closed: generator-level injection -> weights -> spin categories ->
counting estimators -> pull closure, on both the streamed BeAGLE sample
and a synthetic one (`tests/test_reweight.py`).  NOT closed: the
reco-level leg -- the same pulls after the events have gone through
abconv -> npsim -> EICrecon.  That needs the HepMC3 writer of step 5.D
(`io_hepmc.py`, not written) to put the reweighted sample back on disk in
a form the conversion chain accepts; `reco.py`/`recopseudo.py` currently
fold resolutions analytically rather than reading a reconstructed file.

Statistics
----------
One pseudo-experiment draws, per spin category k,

    n_k ~ Poisson(mu_k),   mu_k = N_target * f_k * <W_k>,

with f_k the category's luminosity share and <W_k> the sample-weighted
mean weight, and then n_k events from the pool with probability
proportional to W_k.  The counts are therefore exactly Poisson, which is
what the analytic error formulas of `polli_fastsim.asymmetries` assume;
the resampling is with replacement, so the pool only has to be large
enough that the *shape*-sensitive estimators (the cos 2phi fit) do not see
their own bootstrap correlation -- the counting estimators (A_par, A_zz)
are insensitive to it by construction, because the count is drawn before
the events are.

Caveat on the external sample's own weights: a sample carrying per-event
generator weights must pass them as `ExternalSample.weight`; they enter
both mu_k and the resampling probability.  The BeAGLE EVGEN files are
unweighted (GenCrossSection weight 1), so the loader leaves them at one.
"""

import numpy as np

from polli_fastsim import asymmetries as _asym
from polli_fastsim.kinematics import (kinematic_mask, scattered_electron,
                                      y_from_xq2)

from .xsec import EventSpinState


# --- the external sample ------------------------------------------------


class ExternalSample:
    """An unpolarized (x, Q2, phi) event pool with a fixed sqrt(s)^2.

    `phi` is the azimuth of the lepton scattering plane in the lab, the
    same angle Mode G's `sample.InclusiveSampler` stores: the kernel only
    ever uses phi' = phi - phi_S.  `weight` is the generator weight (one
    for an unweighted sample).
    """

    def __init__(self, x, q2, phi, s, weight=None, extra=None):
        self.x = np.asarray(x, dtype=float)
        self.q2 = np.asarray(q2, dtype=float)
        self.phi = np.mod(np.asarray(phi, dtype=float), 2.0 * np.pi)
        if not (self.x.shape == self.q2.shape == self.phi.shape):
            raise ValueError("x, q2, phi must have the same shape")
        if self.x.ndim != 1:
            raise ValueError("ExternalSample takes flat arrays")
        self.s = float(s)
        if self.s <= 0.0:
            raise ValueError("s must be positive")
        self.weight = (np.ones_like(self.x) if weight is None
                       else np.asarray(weight, dtype=float))
        if self.weight.shape != self.x.shape:
            raise ValueError("weight must match x")
        if np.any(self.weight < 0.0):
            raise ValueError("negative generator weights are not supported")
        self.extra = dict(extra or {})

    def __len__(self):
        return int(self.x.size)

    @property
    def y(self):
        return y_from_xq2(self.x, self.q2, self.s)

    def select(self, mask):
        """A new sample keeping the events `mask` selects."""
        mask = np.asarray(mask, dtype=bool)
        return ExternalSample(
            self.x[mask], self.q2[mask], self.phi[mask], self.s,
            weight=self.weight[mask],
            extra={k: np.asarray(v)[mask] for k, v in self.extra.items()})

    def in_scenario(self, scenario, electron_energy):
        """The `fom.Scenario` kinematic + scattered-electron cuts, applied
        event by event -- the same cuts `fom.project_rates` applies at cell
        centres and `sample.InclusiveSampler._in_acceptance` re-applies to
        in-cell draws, so a Mode-W run and a Mode-G run of the same
        scenario see the same window."""
        y = np.clip(self.y, 1e-12, 1.0 - 1e-12)
        ok = kinematic_mask(self.x, self.q2, self.s, q2_min=scenario.q2_min,
                            y_min=scenario.y_min, y_max=scenario.y_max,
                            w2_min=scenario.w2_min)
        e_p, _th, eta = scattered_electron(self.x, y, self.s, electron_energy)
        return (ok & (eta >= scenario.eta_min) & (eta <= scenario.eta_max)
                & (e_p >= scenario.e_prime_min))


def sample_from_beagle_rows(rows, kinds, electron_beam):
    """(x, Q2, y, phi) per event from `tools/analysis/dump_spectators.py`
    rows: the beam-ion row 'B' and the scattered-electron row 'E'.

    The dumper writes no beam-electron row (it skips status-4 leptons), so
    the incoming lepton four-vector is passed in: `electron_beam` is
    (px, py, pz, E) in GeV.  Everything else is per event, which is what
    matters, because the afterburned ion beam carries the crossing angle
    and the divergence and the DIS invariants are formed against that
    exact vector:

        q = k - k',  Q2 = -q.q,  x = Q2/(2 P.q),  y = P.q/P.k

    with P the beam-ion row.  The BeAGLE 'eH2' files record the struck
    NUCLEON as the status-4 ion beam, so P is already per nucleon and x is
    the Bjorken x of the nucleon, the variable the per-nucleon structure
    functions of `polli_fastsim` are functions of.

    phi is the lab azimuth of the scattered electron, i.e. of the lepton
    scattering plane about the beam axis.  The kernel needs the azimuth
    about the virtual-photon direction; the two agree up to the O(25 mrad)
    crossing angle, which is why `ModeWReweighter` is documented as
    generator-level and the residual is quoted, not ignored.
    """
    rows = np.asarray(rows)
    kinds = np.asarray(kinds).astype(str)
    kb = electron_beam
    k = np.asarray([float(kb[3]), float(kb[0]), float(kb[1]), float(kb[2])])
    beams = rows[kinds == "B"]
    elecs = rows[kinds == "E"]
    by_evt = {}
    for b in beams:
        ent = by_evt.setdefault(int(b["ievt"]), [None, None])
        # more than one status-4 hadron row can be written for a nucleus;
        # the ion side is the forward-going one
        if ent[0] is None or b["pz"] > ent[0]["pz"]:
            ent[0] = b
    for e in elecs:
        ent = by_evt.get(int(e["ievt"]))
        if ent is not None:
            ent[1] = e
    ievt, xs, q2s, ys, phis = [], [], [], [], []
    for i in sorted(by_evt):
        b, e = by_evt[i]
        if b is None or e is None:
            continue
        p = np.array([b["e"], b["px"], b["py"], b["pz"]], dtype=float)
        kp = np.array([e["e"], e["px"], e["py"], e["pz"]], dtype=float)
        q = k - kp
        q2 = -_minkowski(q, q)
        pq = _minkowski(p, q)
        pk = _minkowski(p, k)
        if q2 <= 0.0 or pq <= 0.0 or pk <= 0.0:
            continue
        ievt.append(i)
        q2s.append(q2)
        xs.append(q2 / (2.0 * pq))
        ys.append(pq / pk)
        phis.append(np.arctan2(e["py"], e["px"]))
    return (np.asarray(ievt, dtype=int), np.asarray(xs), np.asarray(q2s),
            np.asarray(ys), np.mod(np.asarray(phis), 2.0 * np.pi))


def _minkowski(a, b):
    """(E, px, py, pz) dot product with the (+, -, -, -) metric."""
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]


def unpolarized_pool(kernel, s, n, x_range=(1e-3, 0.5), q2_range=(1.0, 100.0),
                     rng=None, scenario=None, electron_energy=None,
                     oversample=4):
    """A synthetic stand-in for an external unpolarized sample.

    Events are proposed log-uniform in (x, Q2) over the window, accepted
    against `kernel.dsigma_unpol` (so the pool carries the unpolarized
    cross section in its density, not in a weight, exactly as a generator
    file does), and given a uniform phi.  With a `scenario` and the
    electron beam energy the `fom.Scenario` cuts are applied too.

    This is what `tests/test_reweight.py` runs on: it needs no network and
    no container, and it is unpolarized by construction, which is the one
    property the reweighter's closure depends on.  Fewer than `n` events
    come back unless `oversample` is raised -- the proposal is drawn
    `oversample * n` times and truncated.
    """
    rng = rng or np.random.default_rng(20260915)
    m = int(oversample * n)
    lx = rng.uniform(np.log(x_range[0]), np.log(x_range[1]), size=m)
    lq = rng.uniform(np.log(q2_range[0]), np.log(q2_range[1]), size=m)
    x, q2 = np.exp(lx), np.exp(lq)
    y = q2 / (s * x)
    ok = (y > 0.0) & (y < 1.0)
    x, q2 = x[ok], q2[ok]
    # d2sigma/dln x dln Q2 = x Q2 d2sigma/dx dQ2 is the log-uniform weight
    dens = kernel.dsigma_unpol(x, q2, s) * x * q2
    dens = np.where(np.isfinite(dens) & (dens > 0.0), dens, 0.0)
    keep = rng.uniform(size=x.size) * dens.max() < dens
    x, q2 = x[keep][:n], q2[keep][:n]
    phi = rng.uniform(0.0, 2.0 * np.pi, size=x.size)
    pool = ExternalSample(x, q2, phi, s)
    if scenario is not None:
        if electron_energy is None:
            raise ValueError("a scenario needs the electron beam energy")
        pool = pool.select(pool.in_scenario(scenario, electron_energy))
    return pool


def load_beagle_spectator_csv(path, electron_beam, s=None):
    """`ExternalSample` from a dump_spectators.py CSV.

    Streamed with the recipe of `tools/beagle/README.md`, which needs an
    eic-shell image with pyHepMC3 and its rootIO bindings, run with its
    working directory inside $ROOTSYS so ROOT finds its module map:

        singularity exec --pwd /opt/local/lib/root $SIF python3 \\
            $PWD/tools/analysis/dump_spectators.py \\
            $B/<file>.hepmc3.tree.root $PWD/ed.csv --nevents 20000

    `s` defaults to 4 E_e <|p_z|> over the recorded beam rows -- the
    massless s the kernel's y = Q2/(s x) is built on, not the exact
    (P + k)^2; pass it to pin the value.  The exact per-event y = P.q/P.k
    is kept in `extra["y_beam"]` for anything that wants the nucleon-mass
    term back, alongside the BeAGLE event index in `extra["ievt"]`.
    """
    raw = np.genfromtxt(path, delimiter=",", names=True, dtype=None,
                        encoding="utf-8")
    kinds = raw["kind"].astype(str)
    ievt, x, q2, y, phi = sample_from_beagle_rows(raw, kinds, electron_beam)
    if s is None:
        beams = raw[kinds == "B"]
        pz = np.mean(np.abs(beams["pz"]))
        s = 4.0 * float(electron_beam[3]) * float(pz)
    return ExternalSample(x, q2, phi, s,
                          extra={"ievt": ievt, "y_beam": y})


# --- the reweighter -----------------------------------------------------


class ModeWReweighter:
    """Per-event polarized weights for one kernel on one external sample.

    The structure-function tables are built once on the pool's (x, Q2) --
    that is the expensive step -- and the (w_avg, a_1, a_2) triple is
    cached per (category, m).
    """

    def __init__(self, kernel, sample, with_perp=False):
        self.kernel = kernel
        self.sample = sample
        self.with_perp = bool(with_perp)
        self.s = sample.s
        self.tables = kernel.tables(sample.x, sample.q2, with_g2=with_perp)
        self._amp_cache = {}
        self._w_cache = {}
        self._prob_cache = {}

    # -- amplitudes ------------------------------------------------------

    @staticmethod
    def _key(category, m=None):
        """Everything about a category the weights depend on.  Not the
        object identity and not its name alone: `bookkeeping.with_offset`
        makes copies that differ only in the luminosity share, which the
        weights must NOT see, while two plans can reuse a name for
        different populations."""
        return (float(category.j), int(category.lam_e), float(category.pe),
                float(category.theta_s), float(category.phi_s),
                tuple(float(p) for p in category.populations),
                None if m is None else float(m))

    def _state(self, category, m):
        return EventSpinState(lam_e=category.lam_e, pe=category.pe,
                              j=category.j, m=m, theta_s=category.theta_s,
                              phi_s=category.phi_s)

    def amplitudes(self, category, m):
        key = self._key(category, m)
        if key not in self._amp_cache:
            self._amp_cache[key] = self.kernel.amplitudes(
                self.tables, self.sample.x, self.sample.q2, self.s,
                self._state(category, m), with_perp=self.with_perp)
        return self._amp_cache[key]

    def event_weights(self, category):
        """W_i for every pool event, averaged over the fill populations.

        Cached per category: a pull scan asks for the same three or four
        weight vectors a few hundred times and they cost O(pool) each.

        Raises if the density goes negative anywhere: a negative weight is
        not a probability, and clipping it would dilute the modulation and
        skew the (x, Q2) mixture exactly as `sample.InclusiveSampler.
        _check_positive` describes for Mode G.
        """
        key = self._key(category)
        cached = self._w_cache.get(key)
        if cached is not None:
            return cached
        pops = np.asarray(category.populations, dtype=float)
        phip = self.sample.phi - category.phi_s
        w = np.zeros_like(self.sample.x)
        for p_m, mm in zip(pops, _m_values(category.j)):
            if p_m <= 0.0:
                continue
            w_avg, a_1, a_2 = self.amplitudes(category, mm)
            w += p_m * (1.0 + w_avg + a_1 * np.cos(phip)
                        + a_2 * np.cos(2.0 * phip))
        worst = float(np.min(w)) if w.size else 1.0
        if worst < 0.0:
            k = int(np.argmin(w))
            raise ValueError(
                "negative Mode-W weight for category %s: %d of %d events, "
                "worst W = %.4f at x = %.4g, Q2 = %.4g, phi = %.3f -- the "
                "kernel amplitude exceeds unity there; reduce the injected "
                "asymmetry." % (category.name, int((w < 0.0).sum()), w.size,
                                worst, float(self.sample.x[k]),
                                float(self.sample.q2[k]),
                                float(self.sample.phi[k])))
        w.flags.writeable = False
        self._w_cache[key] = w
        return w

    def _draw_prob(self, category):
        """(p_i, <W>) with p_i ~ w_i W_i the resampling law -- cached, so a
        pull scan pays for it once per category rather than once a trial."""
        key = self._key(category)
        cached = self._prob_cache.get(key)
        if cached is None:
            sw = self.sample.weight
            w = self.event_weights(category) * sw
            tot = float(w.sum())
            if tot <= 0.0:
                raise ValueError("category %s has zero total weight"
                                 % category.name)
            cached = (w / tot, tot / float(sw.sum()))
            self._prob_cache[key] = cached
        return cached

    def mean_weight(self, category):
        """<W>_pool, the category's rate relative to the unpolarized one."""
        return self._draw_prob(category)[1]

    # -- the injected truth ----------------------------------------------

    def _pool_average(self, values):
        return float(np.average(np.asarray(values, dtype=float),
                                weights=self.sample.weight))

    def truth_a_parallel(self):
        """<A_par>_pool -- what the helicity-flip estimator must return.

        Computed through `polli_fastsim.asymmetries`, not through the
        kernel's own `a_parallel`, so that the closure tests the weight
        assembly (population average, geometry, resampling, estimator) and
        not just the formula.  It follows the kernel's `target_mass`
        switch: the exact finite-gamma form when it is on (the default),
        the massless D(y) g1/F1 when it is off.
        """
        t, x, q2 = self.tables, self.sample.x, self.sample.q2
        y = q2 / (self.s * x)
        g2 = t["g2"] if self.kernel.target_mass else None
        return self._pool_average(
            _asym.a_parallel(t["g1"], t["f1"], y, x, q2,
                             r_func=self.kernel.r_func, g2=g2))

    def truth_azz(self):
        """<A_zz>_pool -- what the tensor thirds estimator must return."""
        t, x, q2 = self.tables, self.sample.x, self.sample.q2
        y = q2 / (self.s * x)
        return self._pool_average(
            _asym.azz(t["b1"], t["f1"], t["f2"], x, y, b2=t["b2"]))

    def truth_cos2phi(self):
        """<a_cos2phi>_pool -- the Delta-sector cos 2phi' amplitude per
        unit P_zz at the transverse working point."""
        t, x, q2 = self.tables, self.sample.x, self.sample.q2
        y = q2 / (self.s * x)
        return self._pool_average(
            _asym.a_cos2phi(t["delta"], t["f1"], t["f2"], x, y))

    # -- pseudo-experiments ----------------------------------------------

    def expected_counts(self, plan, n_target):
        """{category name: mu_k} at `n_target` total events on an
        unpolarized fill.  mu_k = n_target * f_k * <W_k>, with f_k the
        plan's absolute luminosity share (offsets included)."""
        shares = plan.lumi_shares(float(n_target))
        return {c.name: shares[c.name] * self.mean_weight(c)
                for c in plan.categories}

    def draw_category(self, category, n, rng):
        """`n` events resampled from the pool with probability ~ w_i W_i.

        Returns the Mode-G event-dict shape (x, q2, y, phi, idx + labels)
        so the same downstream code reads a Mode-W and a Mode-G sample.
        """
        prob, _ = self._draw_prob(category)
        idx = (rng.choice(prob.size, size=int(n), p=prob) if n > 0
               else np.empty(0, dtype=int))
        x, q2 = self.sample.x[idx], self.sample.q2[idx]
        return {"x": x, "q2": q2, "y": y_from_xq2(x, q2, self.s),
                "phi": self.sample.phi[idx], "idx": idx,
                "category": category.name, "lam_e": category.lam_e}

    def run_pseudo_experiment(self, plan, n_target, rng=None, poisson=True,
                              draw_events=True):
        """One Mode-W pseudo-experiment.

        Returns ({name: events}, {name: mu}).  With `draw_events=False`
        only the counts are drawn -- enough for the counting estimators
        (A_par, A_zz), and the cheap path for a long pull scan, since no
        resampling is done.
        """
        rng = rng or np.random.default_rng(20260915)
        mus = self.expected_counts(plan, n_target)
        events = {}
        for cat in plan.categories:
            mu = mus[cat.name]
            n = int(rng.poisson(mu)) if poisson else int(round(mu))
            if draw_events:
                events[cat.name] = self.draw_category(cat, n, rng)
            else:
                events[cat.name] = {"n": n, "category": cat.name,
                                    "lam_e": cat.lam_e}
        return events, mus


def _m_values(j):
    n = int(round(2 * j + 1))
    return np.array([j - i for i in range(n)])


def counts_of(events):
    """{name: count} from either event-dict shape of run_pseudo_experiment."""
    return {k: (int(v["n"]) if "n" in v else int(np.asarray(v["x"]).size))
            for k, v in events.items()}


# --- the closure ---------------------------------------------------------


def pull_scan(reweighter, plan, n_target, estimator, truth, err_func,
              ntrials=200, seed=20260915, draw_events=False):
    """`ntrials` pseudo-experiments -> (pulls, estimates, counts).

    `estimator(counts_dict)` returns the asymmetry estimate and
    `err_func(n_total)` the expected error the pull divides by, so the
    same scan serves A_par and A_zz.  This is the acceptance measurement
    of plans/05 step 5.C: |mean(pull)| and std(pull).
    """
    rng = np.random.default_rng(seed)
    est, ntot = [], []
    for _ in range(int(ntrials)):
        events, _ = reweighter.run_pseudo_experiment(
            plan, n_target, rng=rng, draw_events=draw_events)
        counts = counts_of(events)
        est.append(estimator(counts))
        ntot.append(sum(counts.values()))
    est = np.asarray(est, dtype=float)
    ntot = np.asarray(ntot, dtype=float)
    pulls = (est - truth) / err_func(ntot)
    return pulls, est, ntot
