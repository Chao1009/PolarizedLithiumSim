"""Reconstructed-level pseudo-experiments (plans/07 WP3).

Composes the measured-quantity layer of `reco.py` with the generator
kernel of `sample.InclusiveSampler` and the coherent scenario of
`coherent.py` into the two analyses an experiment runs:

INCLUSIVE (RecoResponse + measure_inclusive)
  1. per accepted sampler cell, n_mc pseudo-events log-uniform in (x, Q2)
     with weight sigma_c / n_mc [pb] -- an importance-sampled response
     that gives every analysis bin the same MC statistics;
  2. the scattered electron: head-on four-vector -> lab (25 mrad
     crossing) -> Gaussian smearing of E' (EMCal or tracking) and of the
     track direction -> back to the head-on frame;
  3. reconstruction: Q2 by the electron method, y hadronic (parametrized
     Sigma/JB resolution) -> x by the mixed method, or the electron method
     alone; reco-level selection cuts; eps_eID(eta);
  4. the physics azimuth phi' from the covariant formula on the true and
     on the smeared four-vectors -> per-event cos 2phi' dilution;
  5. EXACT expected phi' counts per spin state in any reconstructed bin
     (bin integrals of den_f + P_f A cos 2phi' per event, cf.
     sample.phi_histogram_pseudo), Poisson-fluctuated at any luminosity;
  6. the spin-state-sorted ratio fit (reco.harmonic_ratio_fit) and the
     conversion to Delta with an MC-derived bin-centering factor.

COHERENT (CoherentResponse + measure_coherent)
  recoil (t, phi_t, x_P) -> exact four-vector -> Roman-Pot emulation
  (divergence smearing, rectangular/elliptical cutout) -> reconstructed
  |t| = pT^2 and beta = phi_t - phi_S -> exact expected 2-D (alpha, beta)
  counts per spin state with the injected harmonics (a_e cos 2alpha,
  a_t(t) cos 2beta, a_m cos(alpha+beta), unpolarized u_1, u_2) -> the
  two-azimuth ratio fit (reco.harmonic_ratio_fit_2d).

Everything is per nucleon and in the head-on frame conventions of the
package; the hadronic final state is NOT generated -- its only entry is
the parametrized y resolution (open question #21).
"""

from dataclasses import dataclass, replace

import numpy as np

from polli_fastsim.kinematics import w2

from . import reco
from .spin import m_values


# --- inclusive ---------------------------------------------------------------

@dataclass
class RecoModel:
    """Simple-detector model + reconstruction choices + analysis cuts."""
    energy: str = "emcal"          # "emcal" | "tracking" | "best"
    emcal_stoch: float = 0.02
    emcal_const: float = 0.01
    y_method: str = "mixed"        # "mixed" (Q2_e, y_had) | "electron"
    y_had_res: float = 0.20        # hadronic-method dy/y (assumption)
    beam_e_spread: float = 1.0e-3
    eid: bool = True
    xing: float = reco.XING_IP6
    phi_s: float = np.pi / 2.0     # vertical alignment axis
    # reconstructed-level analysis selection
    q2_min: float = 1.0
    y_min: float = 0.01
    y_max: float = 0.95
    w2_min: float = 10.0
    eta_min: float = -3.5
    eta_max: float = 3.5
    e_prime_min: float = 0.5

    def energy_resolution(self, e_prime, eta):
        cal = reco.emcal_resolution(e_prime, self.emcal_stoch, self.emcal_const)
        if self.energy == "emcal":
            return cal
        trk = reco.tracking_resolution(e_prime, eta)
        if self.energy == "tracking":
            return trk
        if self.energy == "best":
            return np.minimum(cal, trk)
        raise ValueError("energy must be 'emcal', 'tracking' or 'best'")


def generator_scenario(analysis, q2_min=0.7, y_min=0.004, y_max=0.985,
                       w2_min=8.0, eta_pad=0.3, e_prime_min=0.3):
    """Looser copy of an analysis fom.Scenario for the GENERATOR: events
    outside the analysis cuts must exist so that they can migrate in."""
    return replace(analysis, q2_min=q2_min, y_min=y_min, y_max=y_max,
                   w2_min=w2_min, eta_min=analysis.eta_min - eta_pad,
                   eta_max=analysis.eta_max + eta_pad,
                   e_prime_min=e_prime_min)


class RecoResponse:
    """Importance-sampled reconstructed-level response of an
    InclusiveSampler under a RecoModel (see module docstring)."""

    def __init__(self, sampler, model=None, n_mc_per_cell=400, rng=None):
        self.sampler = sampler
        self.model = model or RecoModel()
        self.n_mc_per_cell = n_mc_per_cell
        rng = rng or np.random.default_rng(20260824)
        m = self.model
        cfg = sampler.config
        s = sampler.s
        e_e = cfg.electron_energy
        ncell = sampler.xsec_flat.size

        # 1. in-cell positions, generator acceptance, weights ---------------
        cell = np.repeat(np.arange(ncell), n_mc_per_cell)
        u1 = rng.uniform(size=cell.size)
        u2 = rng.uniform(size=cell.size)
        x = np.exp(sampler.logx_lo[cell]
                   + u1 * (sampler.logx_hi[cell] - sampler.logx_lo[cell]))
        q2 = np.exp(sampler.logq2_lo[cell]
                    + u2 * (sampler.logq2_hi[cell] - sampler.logq2_lo[cell]))
        ok = sampler._in_acceptance(x, q2)
        n_kept = np.bincount(cell[ok], minlength=ncell)
        cell, x, q2 = cell[ok], x[ok], q2[ok]
        w = sampler.xsec_flat[cell] / np.maximum(n_kept[cell], 1)
        y = q2 / (s * x)
        n = x.size

        # 2. scattered electron: true -> lab -> smeared -> head-on ---------
        phi_e = rng.uniform(0.0, 2.0 * np.pi, size=n)
        k, p_ion = reco.beam_fourvectors(cfg)
        s_vec = reco.spin_fourvector(m.phi_s)
        kp = reco.electron_fourvector(x, y, s, e_e, phi_e)
        kp_lab = reco.head_on_to_lab(kp, m.xing)
        e_lab = kp_lab[..., 0]
        pmag = np.sqrt((kp_lab[..., 1:] ** 2).sum(axis=-1))
        th_lab = np.arccos(np.clip(kp_lab[..., 3] / pmag, -1.0, 1.0))
        ph_lab = reco.azimuth_about_z(kp_lab)
        eta_lab = -np.log(np.tan(np.minimum(th_lab, np.pi - 1e-9) / 2.0))
        de = m.energy_resolution(e_lab, eta_lab)
        sig_ang = reco.tracking_angular_resolution(eta_lab)
        e_m = e_lab * (1.0 + de * rng.standard_normal(n))
        th_m = np.clip(th_lab + sig_ang * rng.standard_normal(n), 1e-9,
                       np.pi - 1e-9)
        # transverse direction resolution -> azimuth resolution sig/sin th'
        ph_m = ph_lab + sig_ang / np.maximum(np.sin(th_lab), 1e-6) \
            * rng.standard_normal(n)
        kp_m_lab = reco.fourvector(e_m, e_m * np.sin(th_m) * np.cos(ph_m),
                                   e_m * np.sin(th_m) * np.sin(ph_m),
                                   e_m * np.cos(th_m))
        kp_m = reco.lab_to_head_on(kp_m_lab, m.xing)
        e_r = kp_m[..., 0]
        pr = np.sqrt((kp_m[..., 1:] ** 2).sum(axis=-1))
        th_r = np.arccos(np.clip(kp_m[..., 3] / pr, -1.0, 1.0))
        eta_r = -np.log(np.tan(np.minimum(th_m, np.pi - 1e-9) / 2.0))

        # 3. kinematic reconstruction and selection -----------------------
        q2_e, y_e, x_e = reco.electron_method(e_r, th_r, e_e, s)
        spread = m.beam_e_spread * rng.standard_normal(n)
        q2_e = q2_e * (1.0 + spread)
        y_e = 1.0 - (1.0 - y_e) * (1.0 + spread)
        if m.y_method == "mixed":
            y_r = reco.hadronic_y(y, m.y_had_res, rng)
            x_r = reco.mixed_method(q2_e, y_r, s)
        elif m.y_method == "electron":
            y_r = y_e
            with np.errstate(invalid="ignore", divide="ignore"):
                x_r = np.where(y_e > 0, q2_e / (s * np.where(y_e > 0, y_e, 1.0)),
                               np.nan)
        else:
            raise ValueError("y_method must be 'mixed' or 'electron'")
        valid = np.isfinite(x_r) & (x_r > 0) & (x_r < 1) & (q2_e > 0)
        x_safe = np.where(valid, x_r, 0.5)
        sel = (valid & (q2_e >= m.q2_min) & (y_r >= m.y_min) & (y_r <= m.y_max)
               & (w2(x_safe, q2_e) >= m.w2_min) & (eta_r >= m.eta_min)
               & (eta_r <= m.eta_max) & (e_m >= m.e_prime_min))
        eff = sel.astype(float)
        if m.eid:
            eff = eff * reco.eps_eid(eta_r)

        # 4. azimuth dilution --------------------------------------------
        phip_true = reco.azimuth_wrt_lepton_plane(k, kp, p_ion, s_vec)
        phip_reco = reco.azimuth_wrt_lepton_plane(k, kp_m, p_ion, s_vec)
        dil = np.cos(2.0 * (phip_reco - phip_true))

        self.cell, self.w = cell, w
        self.x, self.q2, self.y = x, q2, y
        self.x_reco = np.where(valid, x_r, np.nan)
        self.q2_reco, self.y_reco, self.eta_reco = q2_e, y_r, eta_r
        self.eff, self.dil = eff, dil
        self.phi_true = np.mod(phip_true, 2 * np.pi)
        self._cellamp = {}

    # --- per-cell fill amplitudes ---------------------------------------

    def _fill_arrays(self, category):
        """(den_c, num_c) over sampler cells: the fill's phi-averaged rate
        factor 1 + sum p_m w_avg and cos 2phi' coefficient sum p_m a2."""
        key = (category.name, tuple(category.populations), category.theta_s,
               category.phi_s, category.lam_e, category.pe)
        if key not in self._cellamp:
            pops = np.asarray(category.populations, dtype=float)
            den = np.zeros(self.sampler.xsec_flat.size)
            num = np.zeros_like(den)
            for p_m, mm in zip(pops, m_values(category.j)):
                if p_m <= 0.0:
                    continue
                w_avg, _a1, a2 = self.sampler._amplitudes(category, mm)
                den += p_m * (1.0 + w_avg)
                num += p_m * a2
            self._cellamp[key] = (den, num)
        return self._cellamp[key]

    def amplitude_per_event(self, category):
        """The P_zz-normalized cos 2phi' amplitude A at each event's TRUE
        kinematics (fill-independent: num/P_zz)."""
        _den, num = self._fill_arrays(category)
        pzz = float(category.moments()[1])
        return num[self.cell] / pzz

    # --- masks and bin bookkeeping ---------------------------------------

    def mask_reco(self, xlo, xhi, q2lo, q2hi):
        xr = self.x_reco
        with np.errstate(invalid="ignore"):
            return ((xr >= xlo) & (xr < xhi) & (self.q2_reco >= q2lo)
                    & (self.q2_reco < q2hi) & (self.eff > 0))

    def mask_true(self, xlo, xhi, q2lo, q2hi):
        return ((self.x >= xlo) & (self.x < xhi) & (self.q2 >= q2lo)
                & (self.q2 < q2hi))

    def bin_summary(self, xlo, xhi, q2lo, q2hi, category):
        """Purity, efficiency, rate ratio and the true-bin / reco-bin
        amplitudes (reco-bin one includes the phi' dilution) for a bin."""
        mr = self.mask_reco(xlo, xhi, q2lo, q2hi)
        mt = self.mask_true(xlo, xhi, q2lo, q2hi)
        a = self.amplitude_per_event(category)
        we = self.w * self.eff
        n_reco = float(we[mr].sum())
        n_true = float(self.w[mt].sum())
        both = float(we[mr & mt].sum())
        return {
            "sigma_reco_pb": n_reco, "sigma_true_pb": n_true,
            "purity": both / n_reco if n_reco > 0 else np.nan,
            "efficiency": both / n_true if n_true > 0 else np.nan,
            "a_true_bin": (float((self.w[mt] * a[mt]).sum() / n_true)
                           if n_true > 0 else np.nan),
            "a_reco_bin": (float((we[mr] * self.dil[mr] * a[mr]).sum()
                                 / n_reco) if n_reco > 0 else np.nan),
            "dilution_phi": (float((we[mr] * self.dil[mr]).sum() / n_reco)
                             if n_reco > 0 else np.nan),
        }

    # --- expected counts ---------------------------------------------------

    def expected_counts(self, categories, lumi_pb, mask, edges,
                        phi_eff=None, nsub=32):
        """Exact expected phi' counts per category (F, K) for the events
        in `mask` (a reco-bin mask), at integrated luminosity lumi_pb
        [pb^-1] shared by the categories' lumi_fraction.  Optional
        `phi_eff(phi')`: a smooth phi'-dependent efficiency (the thing
        the ratio estimator must cancel)."""
        edges = np.asarray(edges, dtype=float)
        frac = (np.arange(nsub) + 0.5) / nsub
        sub = edges[:-1, None] + frac[None, :] * (edges[1:] - edges[:-1])[:, None]
        epsi = np.ones_like(sub) if phi_eff is None else phi_eff(sub)
        dphi = (edges[1:] - edges[:-1]) / (2.0 * np.pi)
        base = epsi.mean(axis=1) * dphi
        mod = (epsi * np.cos(2.0 * sub)).mean(axis=1) * dphi
        we = (self.w * self.eff)[mask]
        dil = self.dil[mask]
        cells = self.cell[mask]
        out = np.empty((len(categories), edges.size - 1))
        for f, cat in enumerate(categories):
            den, num = self._fill_arrays(cat)
            s_den = float((we * den[cells]).sum())
            s_num = float((we * dil * num[cells]).sum())
            out[f] = lumi_pb * cat.lumi_fraction * (s_den * base + s_num * mod)
        return out


def measure_inclusive(resp, plan, lumi_pb, mask, rng=None, nbins=24,
                      phi_eff=None, poisson=True, with_sin=True,
                      lumi_assumed=None):
    """One reconstructed-level pseudo-measurement of a reco bin: expected
    spin-sorted phi' counts -> Poisson draw -> spin-state ratio fit.
    `lumi_assumed`: luminosity fractions the ANALYSIS believes (default:
    the plan's true shares; pass biased values to emulate a relative-
    luminosity error).  Returns the fit dict plus counts and the truth
    references of `bin_summary` (computed by the caller)."""
    rng = rng or np.random.default_rng(20260824)
    edges = np.linspace(0.0, 2.0 * np.pi, nbins + 1)
    mu = resp.expected_counts(plan.categories, lumi_pb, mask, edges,
                              phi_eff=phi_eff)
    counts = rng.poisson(mu) if poisson else mu
    pzz = [float(c.moments()[1]) for c in plan.categories]
    lum = ([c.lumi_fraction for c in plan.categories] if lumi_assumed is None
           else list(lumi_assumed))
    fit = reco.harmonic_ratio_fit(counts, lum, pzz, edges, with_sin=with_sin)
    fit.update({"counts": counts, "expected": mu, "edges": edges,
                "n": float(mu.sum()), "pzz": pzz, "lumi_fractions": lum})
    return fit


def delta_from_amplitude(fit, summary, delta_center):
    """MC bin-centering with migration: K = Delta_model(x_c, Q2_c) /
    A_reco-bin(model), Delta_hat = A_hat K, dDelta = dA |K| -- the pure
    kinematic inversion plus the in-bin-shape and migration correction
    an experiment takes from the model it fits."""
    k_conv = delta_center / summary["a_reco_bin"]
    return {"delta": fit["amp"] * k_conv, "err": fit["err"] * abs(k_conv),
            "k": k_conv}


# --- coherent ------------------------------------------------------------

class CoherentResponse:
    """Roman-Pot emulated response of the coherent recoil spectrum."""

    def __init__(self, scenario, config, sigma_theta, aspect=1.0,
                 shape="rectangle", n_mc=400000, x_pom_range=(1e-3, 1e-2),
                 t_max=0.5, phi_s=np.pi / 2.0, n_sigma=10.0, rng=None,
                 cut_scale_xy=(1.0, 1.0)):
        """`aspect` = sigma_theta_y / sigma_theta_x (beam divergence
        anisotropy; HERA's proton beam had 45 vs 100 MeV horizontal vs
        vertical pT spread at the IP, ZEUS NPB 816:1); `cut_scale_xy`
        scales the cutout half-widths (n_sigma sigma_x, n_sigma sigma_y):
        the ePIC pots surround a horizontal slot, cut_scale_xy = (2.5, 1)
        (see reco.rp_measure)."""
        rng = rng or np.random.default_rng(20260824)
        self.scenario, self.config = scenario, config
        self.sigma_theta, self.aspect, self.shape = sigma_theta, aspect, shape
        self.cut_scale_xy = tuple(cut_scale_xy)
        _k, p_ion = reco.beam_fourvectors(config)
        t = -scenario.sample_t(n_mc, rng, t_max=t_max)
        lo, hi = np.log10(x_pom_range[0]), np.log10(x_pom_range[1])
        x_pom = 10 ** rng.uniform(lo, hi, size=n_mc)
        t_min = reco.mdot(p_ion, p_ion) * x_pom ** 2 / (1.0 - x_pom)
        ok = -t > t_min
        t, x_pom = t[ok], x_pom[ok]
        phi_t = rng.uniform(0.0, 2.0 * np.pi, size=t.size)
        pp = reco.recoil_fourvector(t, phi_t, x_pom, p_ion)
        m = reco.rp_measure(pp, p_ion, (sigma_theta, sigma_theta * aspect),
                            n_sigma=n_sigma, shape=shape, rng=rng,
                            cut_scale_xy=cut_scale_xy)
        acc = m["accepted"]
        self.cut_pt_xy = m["cut_pt_xy"]
        self.n_produced_mc = n_mc
        self.acceptance = float(acc.sum()) / n_mc
        self.t_true = -t[acc]
        self.t_reco = -m["t_reco"][acc]
        self.beta_true = np.mod(phi_t[acc] - phi_s, 2.0 * np.pi)
        self.beta_reco = np.mod(m["phi_t"][acc] - phi_s, 2.0 * np.pi)
        self.x_pom = x_pom[acc]
        self.w = np.full(self.t_true.size, 1.0 / n_mc)   # per produced recoil
        self.pt_cut = reco.tag_pt_cut(sigma_theta, config.ion_momentum_per_nucleon,
                                      a_beam=config.ion.A, n_sigma=n_sigma)

    def t_bin_fraction(self, tlo, thi):
        sel = (self.t_reco >= tlo) & (self.t_reco < thi)
        return float(self.w[sel].sum())

    def expected_counts_2d(self, n_produced, pzz_list, lumi_fractions,
                           tlo, thi, alpha_edges, beta_edges, a_e, a_t_func,
                           a_m=0.0, u1=0.0, u2=0.0, kappa=0.0):
        """Exact expected (F, Ka, Kb) counts for recoils reconstructed in
        [tlo, thi): per event the alpha integrals of
        1 + u1 cos(a-b) + u2 cos 2(a-b) + P_f[kappa + a_e cos 2a
        + a_t(t) cos 2b + a_m cos(a+b)] (TRUE angles/t inside, RECO beta
        for the binning) are analytic.  `n_produced` = coherent recoils
        produced in the sample (before tagging)."""
        sel = (self.t_reco >= tlo) & (self.t_reco < thi)
        wgt = self.w[sel] * n_produced
        bt, br, tt = self.beta_true[sel], self.beta_reco[sel], self.t_true[sel]
        ae = np.asarray(alpha_edges, dtype=float)
        be = np.asarray(beta_edges, dtype=float)
        lo, hi = ae[:-1][None, :], ae[1:][None, :]
        i0 = (hi - lo) / (2.0 * np.pi)
        i2 = (np.sin(2 * hi) - np.sin(2 * lo)) / (4.0 * np.pi)
        b = bt[:, None]
        j1 = (np.sin(hi - b) - np.sin(lo - b)) / (2.0 * np.pi)
        j2 = (np.sin(2 * (hi - b)) - np.sin(2 * (lo - b))) / (4.0 * np.pi)
        j1p = (np.sin(hi + b) - np.sin(lo + b)) / (2.0 * np.pi)
        a_t = np.asarray(a_t_func(tt), dtype=float)[:, None]
        unpol = i0 + u1 * j1 + u2 * j2
        tens = (kappa * i0 + a_e * i2 + a_t * np.cos(2.0 * b) * i0
                + a_m * j1p)
        kb = np.clip(np.digitize(br, be) - 1, 0, be.size - 2)
        out = np.zeros((len(pzz_list), ae.size - 1, be.size - 1))
        for f, (pf, lf) in enumerate(zip(pzz_list, lumi_fractions)):
            arr = wgt[:, None] * lf * (unpol + pf * tens)
            for j in range(ae.size - 1):
                out[f, j] = np.bincount(kb, weights=arr[:, j],
                                        minlength=be.size - 1)
        return out

    def t_mean_true(self, tlo, thi):
        """Rate-weighted mean TRUE |t| of the events reconstructed in
        [tlo, thi) -- the reference point of the t-shape template."""
        sel = (self.t_reco >= tlo) & (self.t_reco < thi)
        return float((self.w[sel] * self.t_true[sel]).sum() / self.w[sel].sum())

    def basis_means(self, tlo, thi, beta_edges, t_shape=None, t_ref=None):
        """Acceptance-weighted in-bin means of cos/sin of the TRUE beta for
        the events reconstructed into each reco beta bin (the response
        basis of reco.basis_2d): {"c1","s1","c2","s2","c2t"}, each (Kb,).

        "c2t" = <g(t) cos 2beta_true>_k is the template basis of the
        t-dependent deformation term a_t(t) = a_t(t_ref) g(t/t_ref)
        (default g = t/t_ref, the linear model of plans/06 SS6.4b): the
        cutout correlates beta with |t| inside a reco t bin (the blind
        directions hold the larger-|t| recoils), so the coefficient must
        be fitted against the template, not against <cos 2beta> alone --
        the two-component fit a_2(t) = c_def |t| + a_g of plans/06 in
        MC-response form.  The fitted coefficient is a_t at t_ref (the
        bin's rate-weighted mean true |t| by default)."""
        sel = (self.t_reco >= tlo) & (self.t_reco < thi)
        be = np.asarray(beta_edges, dtype=float)
        kb = np.clip(np.digitize(self.beta_reco[sel], be) - 1, 0, be.size - 2)
        wgt = self.w[sel]
        bt = self.beta_true[sel]
        tt = self.t_true[sel]
        if t_ref is None:
            t_ref = float((wgt * tt).sum() / wgt.sum())
        g = (tt / t_ref) if t_shape is None else np.asarray(t_shape(tt / t_ref),
                                                            dtype=float)
        tot = np.bincount(kb, weights=wgt, minlength=be.size - 1)
        bc = 0.5 * (be[:-1] + be[1:])
        d1, d2 = reco._bin_dilutions(be, be)[1], reco._bin_dilutions(be, be)[3]
        out = {"t_ref": t_ref}
        for key, fn, fallback in (("c1", np.cos(bt), np.cos(bc) * d1),
                                  ("s1", np.sin(bt), np.sin(bc) * d1),
                                  ("c2", np.cos(2 * bt), np.cos(2 * bc) * d2),
                                  ("s2", np.sin(2 * bt), np.sin(2 * bc) * d2),
                                  ("c2t", g * np.cos(2 * bt),
                                   np.cos(2 * bc) * d2)):
            num = np.bincount(kb, weights=wgt * fn, minlength=be.size - 1)
            out[key] = np.where(tot > 0, num / np.where(tot > 0, tot, 1.0),
                                fallback)
        return out

    def truth_reference(self, tlo, thi, a_e, a_t_func, a_m=0.0, t_ref=None):
        """Injected coefficients the template fit should return for the
        reco t bin: a_t at t_ref (rate-weighted mean true |t| by default),
        a_e and a_m as injected (the MC basis absorbs the beta smearing);
        `dilution_beta` = <cos 2(beta_reco - beta_true)> for information."""
        sel = (self.t_reco >= tlo) & (self.t_reco < thi)
        wgt = self.w[sel]
        dbeta = self.beta_reco[sel] - self.beta_true[sel]
        tot = wgt.sum()
        t_mean = float((wgt * self.t_true[sel]).sum() / tot)
        if t_ref is None:
            t_ref = t_mean
        return {"a_e": a_e, "a_t": float(np.asarray(a_t_func(t_ref))),
                "a_m": a_m, "t_ref": t_ref, "t_mean_true": t_mean,
                "dilution_beta": float((wgt * np.cos(2.0 * dbeta)).sum() / tot),
                "fraction": float(tot)}


def measure_coherent(cresp, n_produced, plan, tlo, thi, a_e, a_t_func,
                     a_m=0.0, u1=0.0, u2=0.0, kappa=0.0, n_alpha=12,
                     n_beta=24, rng=None, poisson=True):
    """One two-azimuth pseudo-measurement of a reco t bin: expected 2-D
    counts per spin state -> Poisson -> reco.harmonic_ratio_fit_2d with
    the unpolarized (u1, u2) taken as known and the acceptance-weighted
    beta basis from the same response (cresp.basis_means)."""
    rng = rng or np.random.default_rng(20260824)
    ae = np.linspace(0.0, 2.0 * np.pi, n_alpha + 1)
    be = np.linspace(0.0, 2.0 * np.pi, n_beta + 1)
    pzz = [float(c.moments()[1]) for c in plan.categories]
    frac = [c.lumi_fraction for c in plan.categories]
    mu = cresp.expected_counts_2d(n_produced, pzz, frac, tlo, thi, ae, be,
                                  a_e, a_t_func, a_m=a_m, u1=u1, u2=u2,
                                  kappa=kappa)
    counts = rng.poisson(mu) if poisson else mu
    bm = cresp.basis_means(tlo, thi, be)
    fit = reco.harmonic_ratio_fit_2d(counts, frac, pzz, ae, be,
                                     u_coeffs=(u1, u2), beta_means=bm)
    fit["beta_means"] = bm
    fit.update({"counts": counts, "expected": mu, "alpha_edges": ae,
                "beta_edges": be, "n": float(mu.sum()),
                "truth": cresp.truth_reference(tlo, thi, a_e, a_t_func, a_m)})
    return fit
