"""Hadronic final state (HFS) for the reconstruction chain (plans/07 WP3-HFS).

The mixed (e-Sigma) reconstruction of x at the y ~ 0.01 sweet spots needs
the hadronic E - p_z sum, which polligen's inclusive kernel does not
generate.  This module supplies

  1. a compact, generator-independent sample format (`HFSSample`: the
     scattered electron and the list of final-state hadrons of unpolarized
     DIS events in the head-on frame, per-nucleon kinematics);
  2. the exact hadronic sums and the classical kinematic methods on them
     (`hadronic_sums`, `hadronic_kinematics`: Jacquet-Blondel, Sigma,
     double-angle, mixed) -- HERA definitions (Bassler-Bernardi NIM A 361
     (1995) 197; Arratia et al. NIM A 1025 (2022) 166164 Table 1);
  3. the hadron-side detector response (`HadronResponse`): tracker
     coverage and thresholds, tracking efficiency and resolution, EMCal
     response for photons and electrons, HCal response for neutral and
     untracked charged hadrons, calorimeter noise on the sums.  All
     magnitudes are the Yellow Report requirement values (Nucl. Phys. A
     1026 (2022) 122447, calorimeter and tracking requirement tables) and
     the repository's tracking placeholders -- to be replaced by the ePIC
     design values when available;
  4. a vectorized toy string-fragmentation generator (`ToyHFS`) so that the
     chain runs and is tested on a machine without PYTHIA -- a stand-in
     whose numbers are illustrative;
  5. an event library (`HFSLibrary`) that transfers the response of real
     generator events (PYTHIA 8 via tools/pythia8/gen_dis_hfs.py, or the
     toy) onto the importance-sampled pseudo-experiments of recopseudo:
     for a pseudo-event at (x, Q2) the captured fraction of Sigma and the
     p_T ratio of a library event from the same (x, Q2) cell are applied to
     the pseudo-event's exact truth, and the noise is added per event.

Conventions: head-on frame (ion along +z with momentum p_u per nucleon,
electron along -z with energy E_e), four-vectors (E, px, py, pz), per-
nucleon x, y, Q2 and s = 4 E_e p_u; theta_e is the scattered-electron
polar angle from the ION direction (so E'(1 - cos theta_e) is the
electron's E - p_z).  Exact relations used for tests: the total E - p_z
of the final state is 2 E_e (massless), hence Sigma_true = 2 E_e y.
"""

import json

import numpy as np

from polli_fastsim.kinematics import scattered_electron

from . import reco

M_PI, M_PI0, M_K, M_K0, M_P, M_N = 0.13957, 0.13498, 0.49368, 0.49761, 0.93827, 0.93957
NEUTRINOS = (12, 14, 16)


# --- sample container --------------------------------------------------------

class HFSSample:
    """Flat arrays of final-state hadrons (the scattered electron excluded)
    for n events: `offsets` (n+1), `pid`, `charge`, `p4` (N, 4); per event
    `x`, `q2`, `y`, `kp` (n, 4) scattered-electron four-vector, `weight`."""

    def __init__(self, offsets, pid, charge, p4, x, q2, y, kp, weight,
                 e_energy, p_per_nucleon, meta=None):
        self.offsets = np.asarray(offsets, dtype=np.int64)
        self.pid = np.asarray(pid, dtype=np.int64)
        self.charge = np.asarray(charge, dtype=float)
        self.p4 = np.asarray(p4, dtype=float)
        self.x = np.asarray(x, dtype=float)
        self.q2 = np.asarray(q2, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.kp = np.asarray(kp, dtype=float)
        self.weight = np.asarray(weight, dtype=float)
        self.e_energy = float(e_energy)
        self.p_per_nucleon = float(p_per_nucleon)
        self.meta = dict(meta or {})

    @property
    def n_events(self):
        return self.offsets.size - 1

    @property
    def s(self):
        return 4.0 * self.e_energy * self.p_per_nucleon

    def event_index(self):
        """Event index of every particle (N,)."""
        return np.repeat(np.arange(self.n_events),
                         np.diff(self.offsets))

    def save(self, path):
        np.savez_compressed(path, offsets=self.offsets, pid=self.pid,
                            charge=self.charge, p4=self.p4, x=self.x,
                            q2=self.q2, y=self.y, kp=self.kp,
                            weight=self.weight,
                            e_energy=np.array(self.e_energy),
                            p_per_nucleon=np.array(self.p_per_nucleon),
                            meta=np.array(json.dumps(self.meta)))

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as f:
            meta = json.loads(str(f["meta"])) if "meta" in f else {}
            return cls(f["offsets"], f["pid"], f["charge"], f["p4"], f["x"],
                       f["q2"], f["y"], f["kp"], f["weight"],
                       float(f["e_energy"]), float(f["p_per_nucleon"]), meta)

    @classmethod
    def concatenate(cls, samples):
        s0 = samples[0]
        offs = [s0.offsets]
        for s in samples[1:]:
            offs.append(s.offsets[1:] + offs[-1][-1])
        cat = lambda k: np.concatenate([getattr(s, k) for s in samples])
        return cls(np.concatenate(offs), cat("pid"), cat("charge"),
                   cat("p4"), cat("x"), cat("q2"), cat("y"), cat("kp"),
                   cat("weight"), s0.e_energy, s0.p_per_nucleon, s0.meta)


# --- exact sums and kinematic methods -----------------------------------------

def hadronic_sums(p4, offsets, weights=None):
    """Per-event Sigma = sum(E - p_z), p_T components and the summed
    transverse-energy-like scalar sum |p_T| (used for nothing but
    diagnostics), from flat particle arrays.  `weights` (N,) multiply each
    particle (0 = lost)."""
    p4 = np.asarray(p4, dtype=float)
    n_ev = offsets.size - 1
    w = np.ones(p4.shape[0]) if weights is None else np.asarray(weights, float)
    empz = w * (p4[:, 0] - p4[:, 3])
    px = w * p4[:, 1]
    py = w * p4[:, 2]

    def rsum(a):
        out = np.zeros(n_ev)
        nonempty = np.diff(offsets) > 0
        if nonempty.any():
            res = np.add.reduceat(a, offsets[:-1][nonempty])
            out[nonempty] = res
        return out
    return rsum(empz), rsum(px), rsum(py)


def hadronic_kinematics(sigma, ptx, pty, e_prime, theta_e, e_energy, s):
    """Jacquet-Blondel, Sigma, double-angle and mixed kinematics from the
    hadronic sums and the reconstructed electron (theta_e from the ion
    direction).  Returns a dict of arrays; unphysical values are NaN."""
    sigma = np.asarray(sigma, dtype=float)
    pt_h = np.hypot(np.asarray(ptx, dtype=float), np.asarray(pty, dtype=float))
    e_prime = np.asarray(e_prime, dtype=float)
    ct = np.cos(np.asarray(theta_e, dtype=float))
    st = np.sin(np.asarray(theta_e, dtype=float))
    with np.errstate(divide="ignore", invalid="ignore"):
        q2_e = 2.0 * e_energy * e_prime * (1.0 + ct)
        e_empz = e_prime * (1.0 - ct)                      # electron E - p_z
        y_jb = sigma / (2.0 * e_energy)
        q2_jb = pt_h ** 2 / (1.0 - y_jb)
        y_sig = sigma / (sigma + e_empz)
        q2_sig = (e_prime * st) ** 2 / (1.0 - y_sig)
        # double angle: tan(gamma/2) = Sigma / pT_h ; tan(theta_e/2) = st/(1+ct)
        tg = sigma / pt_h
        te = st / (1.0 + ct)
        y_da = tg / (tg + te)
        q2_da = 4.0 * e_energy ** 2 / te / (tg + te)   # = 4E^2 cot(th/2)/(tg+te)
        out = {"q2_e": q2_e, "y_jb": y_jb, "q2_jb": q2_jb,
               "x_jb": q2_jb / (s * y_jb), "y_sigma": y_sig, "q2_sigma": q2_sig,
               "x_sigma": q2_sig / (s * y_sig), "y_da": y_da, "q2_da": q2_da,
               "x_da": q2_da / (s * y_da), "x_mixed": q2_e / (s * y_sig),
               "x_ejb": q2_e / (s * y_jb), "x_eda": q2_e / (s * y_da)}
    for k, v in out.items():
        v = np.asarray(v, dtype=float)
        v[~np.isfinite(v)] = np.nan
        out[k] = v
    return out


# --- hadron-side detector response -------------------------------------------

def _region_res(eta, table):
    """Piecewise (stochastic, constant) resolution terms by eta from a
    table of (eta_lo, eta_hi, stoch, const)."""
    eta = np.asarray(eta, dtype=float)
    a = np.full(eta.shape, np.nan)
    b = np.full(eta.shape, np.nan)
    for lo, hi, sa, sb in table:
        m = (eta >= lo) & (eta < hi)
        a[m] = sa
        b[m] = sb
    return a, b


class HadronResponse:
    """Simple ePIC-like response of the hadronic final state.

    Parameters (defaults = Yellow Report requirement magnitudes; tracking
    from the repository's placeholder tables in `reco`):
      eta_track      tracker coverage |eta| <= 3.5
      pt_min_track   track p_T threshold [GeV] (0.2 GeV; the turn-on width
                     pt_turn = 0.05 GeV)
      eff_track      plateau tracking efficiency
      eta_cal        calorimeter coverage |eta| <= 3.7 (photons, neutral
                     hadrons, and charged particles without a track)
      e_min_photon   EMCal cluster threshold [GeV]
      e_min_nhad     HCal cluster threshold [GeV] for neutral hadrons and
                     untracked charged hadrons
      eff_nhad       neutral-hadron detection efficiency above threshold
      emcal, hcal    resolution tables (eta_lo, eta_hi, stoch, const):
                     sigma_E/E = stoch/sqrt(E) (+) const
      noise_sigma    Gaussian noise [GeV] added per event to Sigma and to
                     each p_T component (the Arratia et al. Sec. 5 effect;
                     0 switches it off)
    """

    EMCAL_YR = ((-4.0, -2.0, 0.02, 0.02), (-2.0, -1.0, 0.07, 0.02),
                (-1.0, 1.0, 0.10, 0.02), (1.0, 4.0, 0.10, 0.02))
    HCAL_YR = ((-4.0, -1.0, 0.50, 0.10), (-1.0, 1.0, 1.00, 0.10),
               (1.0, 4.0, 0.50, 0.10))

    def __init__(self, eta_track=3.5, pt_min_track=0.2, pt_turn=0.05,
                 eff_track=0.95, eta_cal=3.7, e_min_photon=0.1,
                 e_min_nhad=0.5, eff_nhad=0.9, emcal=None, hcal=None,
                 noise_sigma=0.05, perfect=False):
        self.eta_track = eta_track
        self.pt_min_track = pt_min_track
        self.pt_turn = pt_turn
        self.eff_track = eff_track
        self.eta_cal = eta_cal
        self.e_min_photon = e_min_photon
        self.e_min_nhad = e_min_nhad
        self.eff_nhad = eff_nhad
        self.emcal = emcal or self.EMCAL_YR
        self.hcal = hcal or self.HCAL_YR
        self.noise_sigma = noise_sigma
        self.perfect = perfect

    def describe(self):
        if self.perfect:
            return "perfect hadron response (all final-state hadrons, no smearing)"
        return ("tracker |eta|<=%.1f pT>%.2f GeV eff %.2f; cal |eta|<=%.1f, "
                "photon E>%.2f, n-had E>%.1f eff %.2f; noise %.0f MeV"
                % (self.eta_track, self.pt_min_track, self.eff_track,
                   self.eta_cal, self.e_min_photon, self.e_min_nhad,
                   self.eff_nhad, 1e3 * self.noise_sigma))

    # per-particle reconstruction ---------------------------------------
    def reconstruct_particles(self, p4, pid, charge, rng):
        """Measured four-vectors (N, 4) and a detection weight (N,) in
        {0, 1}; lost particles get weight 0.  Neutrinos are always lost."""
        p4 = np.asarray(p4, dtype=float)
        pid = np.asarray(pid)
        charge = np.asarray(charge, dtype=float)
        n = p4.shape[0]
        if self.perfect:
            w = np.ones(n)
            w[np.isin(np.abs(pid), NEUTRINOS)] = 0.0
            return p4.copy(), w
        p = np.sqrt((p4[:, 1:] ** 2).sum(axis=1))
        p_safe = np.maximum(p, 1e-12)
        pt = np.hypot(p4[:, 1], p4[:, 2])
        cth = np.clip(p4[:, 3] / p_safe, -1.0, 1.0)
        eta = -np.log(np.tan(np.clip(np.arccos(cth), 1e-9, np.pi - 1e-9) / 2.0))
        e = p4[:, 0]
        out = np.zeros_like(p4)
        w = np.zeros(n)
        is_nu = np.isin(np.abs(pid), NEUTRINOS)
        is_photon = pid == 22
        is_e = np.abs(pid) == 11
        is_charged = (charge != 0) & ~is_nu
        in_track = is_charged & (np.abs(eta) <= self.eta_track)
        in_cal = (np.abs(eta) <= self.eta_cal) & ~is_nu
        # --- tracks: efficiency turn-on, momentum + direction smearing
        eff = self.eff_track / (1.0 + np.exp(-(pt - self.pt_min_track) / self.pt_turn))
        eff = np.where(pt < 0.5 * self.pt_min_track, 0.0, eff)
        tracked = in_track & (rng.uniform(size=n) < eff)
        if tracked.any():
            idx = np.flatnonzero(tracked)
            sig_p = reco.tracking_resolution(p[idx], eta[idx])
            sig_a = reco.tracking_angular_resolution(eta[idx])
            p_m = p[idx] * (1.0 + sig_p * rng.standard_normal(idx.size))
            th = np.arccos(cth[idx]) + sig_a * rng.standard_normal(idx.size)
            th = np.clip(th, 1e-9, np.pi - 1e-9)
            ph = np.arctan2(p4[idx, 2], p4[idx, 1]) + sig_a / np.maximum(
                np.sin(np.arccos(cth[idx])), 1e-6) * rng.standard_normal(idx.size)
            mass = np.where(is_e[idx], 0.000511, M_PI)   # pion hypothesis
            out[idx, 0] = np.sqrt(p_m ** 2 + mass ** 2)
            out[idx, 1] = p_m * np.sin(th) * np.cos(ph)
            out[idx, 2] = p_m * np.sin(th) * np.sin(ph)
            out[idx, 3] = p_m * np.cos(th)
            w[idx] = 1.0
        # --- calorimeter objects: photons/electrons in the EMCal, hadrons in
        #     the HCal; charged particles only if not tracked
        cal_em = in_cal & ~tracked & (is_photon | is_e) & (e >= self.e_min_photon)
        cal_had = in_cal & ~tracked & ~is_photon & ~is_e & (e >= self.e_min_nhad)
        cal_had = cal_had & (rng.uniform(size=n) < self.eff_nhad)
        for sel, table in ((cal_em, self.emcal), (cal_had, self.hcal)):
            if not sel.any():
                continue
            idx = np.flatnonzero(sel)
            a, b = _region_res(eta[idx], table)
            a = np.where(np.isnan(a), table[-1][2], a)
            b = np.where(np.isnan(b), table[-1][3], b)
            res = np.sqrt((a / np.sqrt(np.maximum(e[idx], 1e-3))) ** 2 + b ** 2)
            e_m = np.maximum(e[idx] * (1.0 + res * rng.standard_normal(idx.size)), 0.0)
            dirn = p4[idx, 1:] / p_safe[idx, None]
            out[idx, 0] = e_m
            out[idx, 1:] = e_m[:, None] * dirn       # massless cluster
            w[idx] = 1.0
        return out, w

    def reconstruct_sums(self, sample_or_p4, offsets=None, pid=None,
                         charge=None, rng=None):
        """Reconstructed (Sigma, pT_x, pT_y) per event including the noise,
        plus the true sums, for an HFSSample or flat arrays."""
        rng = rng or np.random.default_rng(20260826)
        if isinstance(sample_or_p4, HFSSample):
            s = sample_or_p4
            p4, offsets, pid, charge = s.p4, s.offsets, s.pid, s.charge
        else:
            p4 = np.asarray(sample_or_p4, dtype=float)
        p4m, w = self.reconstruct_particles(p4, pid, charge, rng)
        sig_r, px_r, py_r = hadronic_sums(p4m, offsets, w)
        sig_t, px_t, py_t = hadronic_sums(p4, offsets,
                                          (~np.isin(np.abs(np.asarray(pid)),
                                                    NEUTRINOS)).astype(float))
        n_ev = offsets.size - 1
        if self.noise_sigma > 0 and not self.perfect:
            sig_r = sig_r + self.noise_sigma * rng.standard_normal(n_ev)
            px_r = px_r + self.noise_sigma * rng.standard_normal(n_ev)
            py_r = py_r + self.noise_sigma * rng.standard_normal(n_ev)
        return {"sigma": sig_r, "ptx": px_r, "pty": py_r,
                "sigma_true": sig_t, "ptx_true": px_t, "pty_true": py_t}


# --- toy string-fragmentation generator ----------------------------------------

def _boost(p4, beta):
    """Boost four-vectors (N, 4) by velocity vectors beta (N, 3)."""
    p4 = np.asarray(p4, dtype=float)
    beta = np.asarray(beta, dtype=float)
    b2 = (beta ** 2).sum(axis=1)
    g = 1.0 / np.sqrt(np.maximum(1.0 - b2, 1e-30))
    bp = (beta * p4[:, 1:]).sum(axis=1)
    g2 = np.where(b2 > 0, (g - 1.0) / np.maximum(b2, 1e-30), 0.0)
    out = np.empty_like(p4)
    out[:, 0] = g * (p4[:, 0] + bp)
    out[:, 1:] = p4[:, 1:] + (g2 * bp + g * p4[:, 0])[:, None] * beta
    return out


def _rotate_z_to(v3, axis):
    """Rotate vectors v3 (N, 3) so that the z axis maps onto the unit vectors
    `axis` (N, 3)."""
    ax = np.asarray(axis, dtype=float)
    ax = ax / np.maximum(np.sqrt((ax ** 2).sum(axis=1)), 1e-30)[:, None]
    z = np.array([0.0, 0.0, 1.0])
    c = np.clip(ax[:, 2], -1.0, 1.0)
    k = np.cross(np.broadcast_to(z, ax.shape), ax)   # rotation axis
    kn = np.sqrt((k ** 2).sum(axis=1))
    s = kn
    kh = k / np.maximum(kn, 1e-30)[:, None]
    # Rodrigues: v cos + (k x v) sin + k (k.v)(1 - cos)
    kv = (kh * v3).sum(axis=1)
    out = (v3 * c[:, None] + np.cross(kh, v3) * s[:, None]
           + kh * (kv * (1.0 - c))[:, None])
    flip = c < -0.999999                 # antiparallel: rotate by pi about x
    if flip.any():
        out[flip] = v3[flip] * np.array([1.0, -1.0, -1.0])
    return out


class ToyHFS:
    """Vectorized toy hadronic final state: the hadronic system X (mass W)
    is fragmented as a string between the struck quark and the target
    remnant -- particles uniform in rapidity along the quark axis in the
    X rest frame, Gaussian p_T (sigma_pt per component), species pi+-/pi0/
    K/p, mean multiplicity <n> = n0 + n1 ln W^2, exact four-momentum
    closure by longitudinal rescaling, pi0 -> gamma gamma decayed.  A
    STAND-IN for PYTHIA (tools/pythia8/gen_dis_hfs.py) so that the chain
    runs and is tested locally; its resolutions are illustrative only."""

    SPECIES = ((211, 1.0, M_PI, 0.31), (-211, -1.0, M_PI, 0.31),
               (111, 0.0, M_PI0, 0.30), (321, 1.0, M_K, 0.03),
               (-321, -1.0, M_K, 0.03), (130, 0.0, M_K0, 0.02))

    def __init__(self, e_energy, p_per_nucleon, sigma_pt=0.30, n0=-0.5,
                 n1=1.3, target_mass=0.0):
        # target_mass = 0 keeps the package's massless per-nucleon kinematics
        # exactly (Sigma_true = 2 E_e y); M_P adds the E_N - p_N = M^2/(2 p_u)
        # term (8.8 MeV at 50 GeV/u, 1-9% of Sigma at the sweet spots)
        self.e_energy = float(e_energy)
        self.p_per_nucleon = float(p_per_nucleon)
        self.sigma_pt = sigma_pt
        self.n0, self.n1 = n0, n1
        self.m_target = target_mass

    @property
    def s(self):
        return 4.0 * self.e_energy * self.p_per_nucleon

    def generate(self, x, q2, phi_e, rng, weight=None):
        x = np.asarray(x, dtype=float)
        q2 = np.asarray(q2, dtype=float)
        phi_e = np.asarray(phi_e, dtype=float)
        n_ev = x.size
        s = self.s
        y = q2 / (s * x)
        e_e, p_u = self.e_energy, self.p_per_nucleon
        # beams and the scattered electron (head-on frame)
        k = np.zeros((n_ev, 4)); k[:, 0] = e_e; k[:, 3] = -e_e
        P = np.zeros((n_ev, 4)); P[:, 0] = np.hypot(p_u, self.m_target); P[:, 3] = p_u
        kp = reco.electron_fourvector(x, y, s, e_e, phi_e)
        q = k - kp
        px_tot = P + q                                    # hadronic system X
        w2 = np.maximum(reco.mdot(px_tot, px_tot), 1.0)
        w_mass = np.sqrt(w2)
        # struck-quark direction in the X rest frame
        beta_x = px_tot[:, 1:] / px_tot[:, 0:1]
        pq_out = x[:, None] * P + q                       # massless quark
        pq_cm = _boost(pq_out, -beta_x)
        axis = pq_cm[:, 1:] / np.maximum(np.sqrt((pq_cm[:, 1:] ** 2).sum(axis=1)), 1e-30)[:, None]
        # multiplicity
        mean_n = np.maximum(self.n0 + self.n1 * np.log(w2), 2.0)
        mult = np.maximum(rng.poisson(mean_n), 2)
        ev = np.repeat(np.arange(n_ev), mult)
        n_p = ev.size
        # species
        probs = np.array([sp[3] for sp in self.SPECIES])
        kind = rng.choice(len(self.SPECIES), size=n_p, p=probs / probs.sum())
        pid = np.array([sp[0] for sp in self.SPECIES])[kind]
        charge = np.array([sp[1] for sp in self.SPECIES])[kind]
        mass = np.array([sp[2] for sp in self.SPECIES])[kind]
        # transverse momenta (balanced per event) and rapidities
        px = self.sigma_pt * rng.standard_normal(n_p)
        py = self.sigma_pt * rng.standard_normal(n_p)
        px -= np.repeat(np.bincount(ev, px, n_ev) / mult, mult)
        py -= np.repeat(np.bincount(ev, py, n_ev) / mult, mult)
        mt = np.sqrt(mass ** 2 + px ** 2 + py ** 2)
        # guarantee sum(m_T) < W: shrink p_T of over-heavy events
        smt = np.bincount(ev, mt, n_ev)
        shrink = np.where(smt > 0.9 * w_mass, 0.9 * w_mass / smt, 1.0)
        px *= np.repeat(shrink, mult); py *= np.repeat(shrink, mult)
        mt = np.sqrt(mass ** 2 + px ** 2 + py ** 2)
        ymax = np.maximum(np.log(w_mass / 0.6), 0.3)
        rap = rng.uniform(-1.0, 1.0, size=n_p) * np.repeat(ymax, mult)
        pz = mt * np.sinh(rap)
        pz -= np.repeat(np.bincount(ev, pz, n_ev) / mult, mult)   # sum pz = 0
        # longitudinal rescaling lambda per event: sum sqrt(mt^2 + lam^2 pz^2) = W
        lo = np.zeros(n_ev); hi = np.full(n_ev, 1.0)
        def esum(lam):
            return np.bincount(ev, np.sqrt(mt ** 2 + (np.repeat(lam, mult) * pz) ** 2), n_ev)
        for _ in range(60):                               # bracket
            grow = esum(hi) < w_mass
            if not grow.any():
                break
            hi = np.where(grow, hi * 2.0, hi)
        for _ in range(50):                               # bisection
            mid = 0.5 * (lo + hi)
            below = esum(mid) < w_mass
            lo = np.where(below, mid, lo)
            hi = np.where(below, hi, mid)
        lam = np.repeat(0.5 * (lo + hi), mult)
        pz = lam * pz
        e = np.sqrt(mt ** 2 + pz ** 2)
        p3 = np.stack([px, py, pz], axis=1)
        p3 = _rotate_z_to(p3, axis[ev])
        p4 = np.concatenate([e[:, None], p3], axis=1)
        p4 = _boost(p4, beta_x[ev])
        # pi0 -> gamma gamma
        is_pi0 = pid == 111
        if is_pi0.any():
            idx = np.flatnonzero(is_pi0)
            nd = idx.size
            cth = rng.uniform(-1.0, 1.0, size=nd)
            ph = rng.uniform(0.0, 2.0 * np.pi, size=nd)
            sth = np.sqrt(1.0 - cth ** 2)
            eg = 0.5 * M_PI0
            g1 = np.stack([np.full(nd, eg), eg * sth * np.cos(ph),
                           eg * sth * np.sin(ph), eg * cth], axis=1)
            g2 = g1.copy(); g2[:, 1:] *= -1.0
            bpi = p4[idx, 1:] / p4[idx, 0:1]
            g1 = _boost(g1, bpi); g2 = _boost(g2, bpi)
            keep = ~is_pi0
            p4 = np.concatenate([p4[keep], g1, g2])
            pid = np.concatenate([pid[keep], np.full(nd, 22), np.full(nd, 22)])
            charge = np.concatenate([charge[keep], np.zeros(2 * nd)])
            ev = np.concatenate([ev[keep], ev[idx], ev[idx]])
            order = np.argsort(ev, kind="stable")
            p4, pid, charge, ev = p4[order], pid[order], charge[order], ev[order]
        counts = np.bincount(ev, minlength=n_ev)
        offsets = np.concatenate([[0], np.cumsum(counts)])
        wgt = np.ones(n_ev) if weight is None else np.asarray(weight, float)
        return HFSSample(offsets, pid, charge, p4, x, q2, y, kp, wgt, e_e, p_u,
                         meta={"generator": "ToyHFS", "sigma_pt": self.sigma_pt,
                               "n0": self.n0, "n1": self.n1})


def toy_library_sample(sampler, n_events, rng, per_cell_uniform=True):
    """Generate a ToyHFS sample over the accepted cells of an
    InclusiveSampler: n_events (x, Q2) points drawn uniformly over the
    cells (so every cell of the response has library events) and
    log-uniform inside each cell; phi_e uniform."""
    ncell = sampler.xsec_flat.size
    cell = (rng.integers(0, ncell, size=n_events) if per_cell_uniform
            else rng.choice(ncell, size=n_events,
                            p=sampler.xsec_flat / sampler.xsec_flat.sum()))
    u1 = rng.uniform(size=n_events)
    u2 = rng.uniform(size=n_events)
    x = np.exp(sampler.logx_lo[cell] + u1 * (sampler.logx_hi[cell] - sampler.logx_lo[cell]))
    q2 = np.exp(sampler.logq2_lo[cell] + u2 * (sampler.logq2_hi[cell] - sampler.logq2_lo[cell]))
    ok = sampler._in_acceptance(x, q2)
    x, q2 = x[ok], q2[ok]
    cfg = sampler.config
    toy = ToyHFS(cfg.electron_energy, cfg.ion_momentum_per_nucleon)
    return toy.generate(x, q2, rng.uniform(0, 2 * np.pi, size=x.size), rng)


# --- event library and the response transfer ----------------------------------

class HFSLibrary:
    """Library of generator events with their hadron-side response,
    binned in (log x, log Q2).  `transfer(x, q2, rng)` returns, for every
    request, the captured Sigma fraction f = Sigma_reco/Sigma_true (noise
    excluded), the p_T ratio r = |pT_reco|/|pT_true| and the azimuthal
    shift of the reconstructed p_T of a random library event from the
    same cell (nearest populated cell if empty)."""

    def __init__(self, sample, response, nx=48, nq2=36, x_range=None,
                 q2_range=None, rng=None, min_per_cell=1):
        rng = rng or np.random.default_rng(20260826)
        self.sample = sample
        self.response = response
        # response without noise: the noise is added per pseudo-event
        quiet = HadronResponse(**{**response.__dict__, "noise_sigma": 0.0})
        r = quiet.reconstruct_sums(sample, rng=rng)
        st = np.maximum(r["sigma_true"], 1e-9)
        self.f_sigma = r["sigma"] / st
        pt_t = np.hypot(r["ptx_true"], r["pty_true"])
        pt_r = np.hypot(r["ptx"], r["pty"])
        self.r_pt = np.where(pt_t > 1e-9, pt_r / np.maximum(pt_t, 1e-9), 1.0)
        self.dphi_pt = np.arctan2(r["pty"], r["ptx"]) - np.arctan2(r["pty_true"], r["ptx_true"])
        # grid
        x_range = x_range or (sample.x.min(), sample.x.max())
        q2_range = q2_range or (sample.q2.min(), sample.q2.max())
        self.lx = np.linspace(np.log(x_range[0]), np.log(x_range[1]), nx + 1)
        self.lq = np.linspace(np.log(q2_range[0]), np.log(q2_range[1]), nq2 + 1)
        self.nx, self.nq2 = nx, nq2
        ci = self._cell(sample.x, sample.q2)
        order = np.argsort(ci, kind="stable")
        self.order = order
        counts = np.bincount(ci, minlength=nx * nq2)
        self.counts = counts
        self.start = np.concatenate([[0], np.cumsum(counts)[:-1]])
        # nearest populated cell for every cell
        self.redirect = self._nearest_populated(counts >= min_per_cell)

    def _cell(self, x, q2):
        i = np.clip(np.searchsorted(self.lx, np.log(x), side="right") - 1, 0, self.nx - 1)
        j = np.clip(np.searchsorted(self.lq, np.log(q2), side="right") - 1, 0, self.nq2 - 1)
        return i * self.nq2 + j

    def _nearest_populated(self, populated):
        nx, nq2 = self.nx, self.nq2
        pop = populated.reshape(nx, nq2)
        red = np.arange(nx * nq2)
        pi, pj = np.nonzero(pop)
        if pi.size == 0:
            raise ValueError("empty HFS library")
        for i in range(nx):
            for j in range(nq2):
                if not pop[i, j]:
                    d = (pi - i) ** 2 + (pj - j) ** 2
                    k = int(np.argmin(d))
                    red[i * nq2 + j] = pi[k] * nq2 + pj[k]
        return red

    def draw(self, x, q2, rng):
        c = self.redirect[self._cell(np.asarray(x, float), np.asarray(q2, float))]
        u = rng.uniform(size=c.size)
        k = self.start[c] + np.floor(u * self.counts[c]).astype(np.int64)
        return self.order[k]

    def transfer(self, x, q2, rng):
        k = self.draw(x, q2, rng)
        return self.f_sigma[k], self.r_pt[k], self.dphi_pt[k]

    def coverage(self):
        """Fraction of grid cells populated and the median count."""
        return float((self.counts > 0).mean()), float(np.median(self.counts[self.counts > 0]))


class HFSResponse:
    """Front end used by recopseudo: hadronic y (and x) of a pseudo-event
    from the library transfer plus per-event noise."""

    def __init__(self, library, method="sigma"):
        if method not in ("sigma", "jb", "da"):
            raise ValueError("method must be 'sigma', 'jb' or 'da'")
        self.library = library
        self.method = method
        self.noise_sigma = library.response.noise_sigma

    def hadronic(self, x, q2, y_true, e_prime_reco, theta_e_reco, e_energy, s,
                 rng):
        """Returns hadronic_kinematics() of the transferred, noisy sums plus
        the captured fraction 'f_sigma'."""
        x = np.asarray(x, float); q2 = np.asarray(q2, float)
        y_true = np.asarray(y_true, float)
        f, r, dphi = self.library.transfer(x, q2, rng)
        sig_true = 2.0 * e_energy * y_true
        pt_true = np.sqrt(np.maximum(q2 * (1.0 - y_true), 0.0))
        # true hadronic p_T is opposite to the electron's; direction irrelevant
        # for the magnitudes used by the methods, keep a fixed axis
        n = x.size
        sig = f * sig_true
        ptx = r * pt_true * np.cos(dphi)
        pty = r * pt_true * np.sin(dphi)
        if self.noise_sigma > 0:
            sig = sig + self.noise_sigma * rng.standard_normal(n)
            ptx = ptx + self.noise_sigma * rng.standard_normal(n)
            pty = pty + self.noise_sigma * rng.standard_normal(n)
        out = hadronic_kinematics(sig, ptx, pty, e_prime_reco, theta_e_reco,
                                  e_energy, s)
        out["f_sigma"] = f
        out["sigma"] = sig
        return out

    def y_hadronic(self, x, q2, y_true, e_prime_reco, theta_e_reco, e_energy, s,
                   rng, floor=1e-4):
        out = self.hadronic(x, q2, y_true, e_prime_reco, theta_e_reco,
                            e_energy, s, rng)
        key = {"sigma": "y_sigma", "jb": "y_jb", "da": "y_da"}[self.method]
        y = out[key]
        y = np.where(np.isfinite(y), y, floor)
        return np.clip(y, floor, 1.0 - 1e-9), out


def truth_kinematics_check(sample, tol=1e-6):
    """Sigma_true = 2 E_e y and pT_h = pT_e for every event (massless
    kinematics); returns the maximal relative deviations."""
    sig, px, py = hadronic_sums(sample.p4, sample.offsets)
    expect = 2.0 * sample.e_energy * sample.y
    d_sig = np.max(np.abs(sig - expect) / np.maximum(expect, 1e-12))
    pt_e = np.hypot(sample.kp[:, 1], sample.kp[:, 2])
    d_pt = np.max(np.abs(np.hypot(px, py) - pt_e) / np.maximum(pt_e, 1e-12))
    return d_sig, d_pt
