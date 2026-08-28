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
     conversion to Delta with an MC-derived bin-centering factor -- either
     the bin-by-bin K of an assumed model (delta_from_amplitude) or, since
     the amplitude is exactly LINEAR in Delta, the shape fitted through
     the forward-folded response itself (RecoResponse.fold,
     fold_shape_fit; plans/08 A6).

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

import copy
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
    emcal_eta_table: bool = False  # True: Yellow Report requirements per
    #   eta region (reco.EMCAL_YR_TABLE) instead of the backward-endcap
    #   specification everywhere (code review F4).  With energy="best"
    #   the tracker then takes over in the barrel, which is what an
    #   experiment does.
    y_method: str = "mixed"        # "mixed" (Q2_e, y_had) | "electron"
    y_source: str = "param"        # "param": reco.hadronic_y Gaussian stand-in
    #   | "hfs": hadronic final state through hfs.HFSResponse (library of
    #   generator events + hadron-side detector response; plans/07 WP3-HFS)
    hadronic_method: str = "sigma" # HFS method for y: "sigma" | "jb" | "da"
    y_had_res: float = 0.25        # hadronic-method dy/y: 0.15-0.30 band;
    #   0.25 = the ePIC kinematic-fit study's own smearing of delta_h and
    #   the ATHENA Fig. 22 value at y ~ 0.01 (-> ~0.10 at y ~ 0.1); see
    #   reco.hadronic_y for the sources (y_source = "param" only)
    beam_e_spread: float = 1.0e-3
    e_scale: float = 1.0           # multiplicative CALIBRATION error on the
    #   measured E', unknown to the analysis.  The bigger of the two energy
    #   levers: with x = Q2_e/(s y_Sigma) and dln y_Sigma/dln E' = -(1-y),
    #   dln x/dln E' = 2 - y ~ 2, twice the hadronic one.
    eid_tilt: float = 0.0          # linear eta slope on eps_eID, normalized
    #   at eta = 0: eps -> eps (1 + eid_tilt * eta).  A FLAT eps_eID error is
    #   identically null (it cancels in every ratio the analysis forms);
    #   only an eta shape moves a number.
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
        cal = reco.emcal_resolution(
            e_prime, self.emcal_stoch, self.emcal_const,
            eta=(eta if self.emcal_eta_table else None))
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

    def __init__(self, sampler, model=None, n_mc_per_cell=400, rng=None,
                 hfs=None):
        """`hfs`: an hfs.HFSResponse (required when model.y_source == "hfs"):
        the hadronic y of every pseudo-event then comes from a generator
        event of the same (x, Q2) cell passed through the hadron-side
        detector response, instead of the Gaussian stand-in."""
        if hfs is not None:
            hfs.check_beams(sampler.config.electron_energy,
                            sampler.config.ion_momentum_per_nucleon)
        self.sampler = sampler
        self.model = model or RecoModel()
        self.n_mc_per_cell = n_mc_per_cell
        rng = rng or np.random.default_rng(20260824)
        m = self.model
        if m.y_source == "hfs" and hfs is None:
            raise ValueError("y_source='hfs' needs an hfs.HFSResponse")
        self.hfs = hfs
        self.sigma_capture = None
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
        e_m = m.e_scale * e_lab * (1.0 + de * rng.standard_normal(n))
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
        # electron-beam energy spread: the event was made at E_e (1 + d) but
        # the analysis reconstructs with the nominal E_e, so
        # Q2_e = 2 E_e E'(1 + cos) is LOW by (1 + d) and 1 - y_e =
        # E'(1 - cos)/(2 E_e) HIGH by it.  (The two carried the same sign
        # until 2026-08-28; with the Sigma-method y the mixed x is then
        # exactly beam-energy independent, as it should be, instead of
        # picking up a fake 2d.)
        spread = m.beam_e_spread * rng.standard_normal(n)
        q2_e = q2_e / (1.0 + spread)
        y_e = 1.0 - (1.0 - y_e) * (1.0 + spread)
        if m.y_method == "mixed":
            if m.y_source == "hfs":
                # hadronic final state of a library event from the same
                # (x, Q2) cell, through the hadron-side response, combined
                # with the RECONSTRUCTED electron (Sigma / JB / DA methods)
                y_r, hk = self.hfs.y_hadronic(x, q2, y, e_r, th_r, e_e, s, rng)
                self.sigma_capture = hk["f_sigma"]
            elif m.y_source == "param":
                y_r = reco.hadronic_y(y, m.y_had_res, rng)
            else:
                raise ValueError("y_source must be 'param' or 'hfs'")
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
            eps = reco.eps_eid(eta_r)
            if m.eid_tilt:
                eps = np.clip(eps * (1.0 + m.eid_tilt * eta_r), 0.0, 1.0)
            eff = eff * eps

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
        self._cellder = {}

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

    # --- forward folding: the response as an operator on Delta ------------

    def delta_response(self, category):
        """Per-cell dA/dDelta -- the EXACT derivative of the P_zz-normalized
        cos 2phi' amplitude with respect to the double-helicity-flip
        structure function at that cell's kinematics.

        A is linear in Delta at fixed (x, Q2): the master formula carries
        Delta in exactly one place, a_2 = -[(1-y)/y^2] c_eff(m) sin^2
        theta_S Delta / D_phi (Hoodbhoy-Jaffe-Manohar NPB 312:571 (1989)
        Eq. (30) as transcribed in xsec.InclusiveKernel.amplitudes), and
        nothing else in W depends on it.  The derivative is therefore the
        same kernel call with the Delta table replaced by 1, and `fold`
        below is exact rather than a linearization
        (test_amplitude_is_exactly_linear_in_delta).
        """
        key = (category.name, tuple(category.populations), category.theta_s,
               category.phi_s, category.lam_e, category.pe)
        if key not in self._cellder:
            smp = self.sampler
            tab = dict(smp.tables)
            tab["delta"] = np.ones_like(smp.x_cells)
            out = np.zeros_like(smp.x_cells)
            for p_m, mm in zip(np.asarray(category.populations, dtype=float),
                               m_values(category.j)):
                if p_m <= 0.0:
                    continue
                _w, _a1, a2 = smp.kernel.amplitudes(
                    tab, smp.x_cells, smp.q2_cells, smp.s,
                    smp._state(category, mm), with_perp=smp.with_perp)
                out = out + p_m * a2
            self._cellder[key] = out / float(category.moments()[1])
        return self._cellder[key]

    def delta_cells(self, delta_func):
        """`delta_func(x, q2, f1)` on the sampler CELL CENTRES -- the same
        (x, Q2, F1) grid on which the modulation amplitudes themselves are
        evaluated (sample.InclusiveSampler fidelity note), so that folding
        a model reproduces its own reco-bin amplitude identically."""
        smp = self.sampler
        return np.asarray(delta_func(smp.x_cells, smp.q2_cells,
                                     smp.tables["f1"]), dtype=float)

    def f1_center(self, x, q2):
        """Per-nucleon F1 at a bin centre, from the response's OWN kernel
        (the same NuclearF2 the amplitudes use), so that a Delta model
        evaluated at the centre and one folded through the response are
        evaluated with the same F1."""
        krn = self.sampler.kernel
        return np.asarray(krn.nf2.f1a(x, q2), dtype=float) / krn.ion.A

    def fold_kernel(self, mask, category):
        """(cells, v, norm) with fold(Delta) = sum_c v_c Delta_c / norm --
        the response of one reco bin as a linear functional of Delta on the
        cell grid.  Sparse (only the cells that reach the bin), so a shape
        fit costs one dot product per bin per trial shape."""
        cells = self.cell[mask]
        we = (self.w * self.eff)[mask]
        v = np.bincount(cells, weights=we * self.dil[mask],
                        minlength=self.sampler.xsec_flat.size)
        nz = np.flatnonzero(v)
        return nz, v[nz] * self.delta_response(category)[nz], float(we.sum())

    def fold(self, delta_func, mask, category, kernel=None):
        """Forward-folded amplitude of a reco bin: the response-weighted
        mean of (dA/dDelta) Delta(x_true, Q2_true) over the events that
        reconstruct into it -- i.e. the cos 2phi' amplitude the experiment
        would measure there if the true structure function were
        `delta_func` (plans/08 A6).

        Weights are the ones bin_summary already uses (sigma * eps_eID *
        selection, times the per-event cos 2 dphi' dilution), so with
        `delta_func` = the model the response was built from this returns
        bin_summary["a_reco_bin"] to floating-point roundoff (measured
        3e-16; the sum is only regrouped by cell).  That identity is what
        makes it a drop-in for the bin-by-bin K of delta_from_amplitude:
        K = Delta_shape(x_c, Q2_c) / fold(Delta_shape, bin) with any shape,
        the model-dependent bin-centering with the migration included.
        """
        nz, v, norm = (self.fold_kernel(mask, category) if kernel is None
                       else kernel)
        if norm <= 0.0:
            return np.nan
        return float((v * self.delta_cells(delta_func)[nz]).sum() / norm)

    def fold_mc_error(self, delta_func, mask, category):
        """1 sigma error of `fold` from the FINITE response sample.

        fold = sum_i z_i / sum_i u_i over the MC events of the bin, with
        u_i = sigma_c/n_c * eff_i and z_i = u_i dil_i (dA/dDelta)_c Delta_c;
        the response is stratified (n_c independent in-cell draws per
        sampler cell), so with r_i = z_i - fold u_i the delta-method
        variance is the sum over cells of the unbiased within-cell variance
        of the n_c draws, n_c/(n_c-1) [sum r_i^2 - (sum r_i)^2/n_c],
        divided by (sum u_i)^2.  Cross-cell terms vanish because the cells
        are drawn independently; n_c is taken as fixed (it is the number of
        draws that passed the GENERATOR acceptance, the conditioning under
        which the weights sigma_c/n_c are defined).  Validated against the
        spread over independent response seeds
        (test_fold_mc_error_matches_the_spread_over_response_seeds).
        """
        nz, v, norm = self.fold_kernel(mask, category)
        if norm <= 0.0:
            return np.nan
        val = float((v * self.delta_cells(delta_func)[nz]).sum() / norm)
        cells = self.cell[mask]
        u = (self.w * self.eff)[mask]
        a = (self.delta_response(category) * self.delta_cells(delta_func))
        r = u * (self.dil[mask] * a[cells] - val)
        ncell = self.sampler.xsec_flat.size
        n_c = np.bincount(self.cell, minlength=ncell).astype(float)
        r_sum = np.bincount(cells, weights=r, minlength=ncell)
        r_sq = np.bincount(cells, weights=r * r, minlength=ncell)
        ok = n_c > 1.0
        var = float((n_c[ok] / (n_c[ok] - 1.0)
                     * (r_sq[ok] - r_sum[ok] ** 2 / n_c[ok])).sum())
        return float(np.sqrt(max(var, 0.0)) / abs(norm))

    # --- expected counts ---------------------------------------------------

    def expected_counts(self, categories, lumi_pb, mask, edges,
                        phi_eff=None, nsub=32):
        """Exact expected phi' counts per category (F, K) for the events
        in `mask` (a reco-bin mask), at integrated luminosity lumi_pb
        [pb^-1] shared by the categories' lumi_fraction.  `phi_eff` is
        None, one smooth eps(phi') common to every category -- which the
        ratio estimator cancels exactly -- or a SEQUENCE of one per
        category, the fill-dependent case it cannot (F1;
        reco.fill_acceptance_bias gives the analytic bias)."""
        edges = np.asarray(edges, dtype=float)
        frac = (np.arange(nsub) + 0.5) / nsub
        sub = edges[:-1, None] + frac[None, :] * (edges[1:] - edges[:-1])[:, None]
        effs = reco.per_fill_acceptance(phi_eff, len(categories))
        dphi = (edges[1:] - edges[:-1]) / (2.0 * np.pi)
        we = (self.w * self.eff)[mask]
        dil = self.dil[mask]
        cells = self.cell[mask]
        out = np.empty((len(categories), edges.size - 1))
        for f, (cat, eff) in enumerate(zip(categories, effs)):
            epsi = np.ones_like(sub) if eff is None else eff(sub)
            base = epsi.mean(axis=1) * dphi
            mod = (epsi * np.cos(2.0 * sub)).mean(axis=1) * dphi
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
    `phi_eff`: None, one efficiency common to every fill, or one per
    fill (reco.per_fill_acceptance).  `lumi_assumed`: luminosity
    fractions the ANALYSIS believes.  The default is the plan's shares,
    which is the right thing only for a plan WITHOUT a relative-luminosity
    offset; a plan built with `rel_lumi_offset != 0` carries the offset in
    its shares, and the published scripts pass the nominal [0.5, 0.5]
    explicitly so that the analysis does not know the truth.  The fill
    polarizations are taken from the plan exactly: the polarimetry scale
    is a normalisation on the amplitude (1:1, pinned by a test) and is
    quoted as a systematic rather than emulated here.  Returns the fit
    dict plus counts and the truth references of `bin_summary` (computed
    by the caller)."""
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


def delta_from_amplitude(fit, summary, delta_center, k_conv=None,
                         k_rel_err=0.0):
    """MC bin-centering with migration: K = Delta_model(x_c, Q2_c) /
    A_reco-bin(model), Delta_hat = A_hat K, dDelta = dA |K| -- the pure
    kinematic inversion plus the in-bin-shape and migration correction
    an experiment takes from the model it fits.

    `k_conv` replaces that K with one the caller derived otherwise -- in
    particular `fold_shape_fit`'s, whose Delta(x) shape is FITTED THROUGH
    THE RESPONSE instead of assumed (plans/08 A6, code review F5) -- and
    `k_rel_err` is the relative error of that factor, added to the
    statistical error in quadrature.  Pass `fold_shape_fit`'s
    "k_rel_all", not "k_rel": the shape-fit and response-MC errors alone
    are 3-10x smaller than the residual prior dependence the tilt family
    cannot absorb, so a bar built from them is smaller than a bias the
    closure scan already measures.  Both default to the published
    behaviour, bit for bit.

    The K error and the statistical error of A_hat come from the same
    data, so the quadrature sum is the conservative reading and not an
    exact treatment; the correlation is not modelled.
    """
    if k_conv is None:
        k_conv = delta_center / summary["a_reco_bin"]
    stat = fit["err"] * abs(k_conv)
    err_k = abs(fit["amp"] * k_conv) * k_rel_err
    return {"delta": fit["amp"] * k_conv,
            "err": float(np.hypot(stat, err_k)) if k_rel_err else stat,
            "err_stat": stat, "err_k": err_k, "k": k_conv,
            "k_rel_err": k_rel_err}


def tilted_shape(base, x_ref, params):
    """`base` deformed by the two-parameter tilt (x/x_ref)^c1 (1-x)^c2.

    The shape family `fold_shape_fit` floats: the same x^a (1-x)^b form
    the repository's own Delta ansaetze use (delta_models.VARIANTS), so
    the tilt is a free low-x power and a free large-x fall-off around the
    prior -- the two directions in which moment_A, moment_B and toy
    actually differ.  params = (norm, c1, c2); (1, 0, 0) returns `base`
    itself, so freezing the tilt recovers today's bin-by-bin K exactly.
    """
    norm, c1, c2 = params

    def shape(x, q2, f1):
        x = np.asarray(x, dtype=float)
        return (norm * np.asarray(base(x, q2, f1), dtype=float)
                * np.power(x / x_ref, c1)
                * np.power(np.maximum(1.0 - x, 1e-12), c2))
    return shape


def _profiled_chi2(fold_rows, amp, wgt):
    """chi2 of A_hat = norm * fold(shape) with the normalization profiled
    out analytically (the fold is LINEAR in Delta, so norm is a linear
    parameter and only the SHAPE parameters are searched)."""
    denom = float((fold_rows * fold_rows * wgt).sum())
    if denom <= 0.0:
        return np.inf, 0.0
    norm = float((amp * fold_rows * wgt).sum() / denom)
    return float((wgt * (amp - norm * fold_rows) ** 2).sum()), norm


def fold_shape_fit(resp, category, bins, base, x_ref=None,
                   bounds=((-3.0, 3.0), (-8.0, 8.0)), n_grid=17, n_zoom=5,
                   mc_error=True, alt_bases=()):
    """Data-driven bin-centering for one Q2 slice: fit a 3-parameter
    Delta(x) shape FOLDED THROUGH THE RESPONSE to the measured amplitudes,
    and return the K it implies per bin (plans/08 A6, code review F5).

    `bins` is a sequence of dicts with the reco-bin mask ("mask"), the bin
    centre ("x", "q2") and the measurement ("amp", "err"); the shape is
    tilted_shape(base, x_ref, (norm, c1, c2)) and the fitted quantity is
        chi2(c1, c2) = sum_b [A_hat_b - norm fold(shape, b)]^2
                             / [dA_b^2 + dA_b(MC)^2]
    with norm profiled out (the fold is linear in Delta) and dA(MC) the
    response MC error of that bin's fold.  The reported
        K_b = shape(x_b, Q2_b) / fold(shape, bin b)
    is independent of norm, so only the two SHAPE parameters matter; the
    bins passed here should extend BEYOND the plotted range, because the
    edge bins of a plot are fed by true x outside it.

    Why a bounded 2-parameter shape and not an unfolding matrix (the
    audit's two constraints): A_hat is a weighted MEAN of per-event
    amplitudes and dA/dDelta changes by ~2x between adjacent x bins, so a
    yield migration matrix mis-weights every off-diagonal element; and
    with a purity near 0.6 that matrix is ill-conditioned, so least
    squares on it oscillates.  A smooth two-parameter tilt searched on a
    bounded grid cannot oscillate bin to bin at all -- that is the
    conditioning guard, and it is pinned by
    test_folded_fit_does_not_oscillate.

    Measured on money plot 7R (n_mc_per_cell = 400, moment_A injected,
    noise-free pseudo-data, `--unfold-scan`): correcting with the WRONG
    shape leaves a residual bias of at most 3.5% / 7.1% / 6.9% over the
    plotted bins of the three Q2 slices with the moment_B prior and
    4.7% / 5.5% / 5.0% with the toy prior, against 22% / 24% / 99% for the
    bin-by-bin K.  At the four sweet spots the model dependence
    K(moment_A)/K(prior) - 1 falls from (-5.0, +8.5, -9.3, +6.3)% and
    (+9.7, +2.5, +10.6, +2.8)% to (-0.9, -1.2, -0.2, +0.3)% and
    (+4.7, +2.8, +2.2, -0.9)%.

    That residual is REAL and is bigger than the fit's own errors: the
    tilt (x/x_ref)^c1 (1-x)^c2 cannot absorb the F1(x, Q2) factor that
    separates the repository's interpretations A and B over three decades
    of x, so the letter cannot quote the shape-fit and response-MC errors
    alone.  `alt_bases` is the family the prior is drawn from -- the fit
    is repeated from each of those shapes on the SAME data, and
        k_prior_rel[b] = max_p |K_p[b] / K[b] - 1|
    is returned as the third error.  It BRACKETS the residual bias
    whenever the truth is a member of that family: with the injected
    shape among `alt_bases` one of the refits IS the closure fit, whose
    own bias is the chain floor (0.03% on money plot 7R), so
    K_injected/K_prior - 1 is the residual up to O(residual^2).  Default
    () -- no alternatives, no third error, today's numbers unchanged.

    No goodness-of-fit penalty is applied to the alternatives, so this is
    an UPPER bound: the spread is a bias and does not shrink with
    luminosity, while the chi2 excess that would let an experiment throw
    an alternative out grows linearly with it.  Measured on money plot 7R
    (noise-free data, chi2 excess of the alternative's folded fit over
    the injected one, ndof = 9): Q2 = 1.14 toy 3.5 at 1 year and 34.7 at
    10 years, moment_B 0.16 and 1.55; Q2 = 14.3 toy 0.03 and 0.27,
    moment_B 0.17 and 1.65.  So at the low-Q2 slice a 10-year measurement
    rejects the shape that sets the 4.7% systematic, while at Q2 = 14.3
    no member of the family is distinguishable even at 10 years and the
    6.1% spread there is irreducible.  The per-alternative chi2 comes
    back in "alt_chi2" (same order as `alt_bases`) so the caller can see
    which alternatives the data disfavour.

    Returns a dict: "params" (norm, c1, c2), "shape" (the callable),
    "k"/"k_err" per bin (the latter the shape-fit and response-MC errors
    of K in quadrature), "k_fit_rel"/"k_mc_rel"/"k_prior_rel" separately,
    "k_rel_all"/"k_err_all" (all three in quadrature -- the honest bar),
    "fold", the "cov" of (c1, c2), "chi2", "ndof" and "at_bound".  A
    minimum with no positive curvature (a flat or saddle chi2, which
    "at_bound" usually accompanies) returns "cov" NaN and "k_fit_rel"
    zero: there is no defensible shape error there, so check them before
    quoting one.
    """
    masks = [b["mask"] for b in bins]
    x_c = np.array([b["x"] for b in bins], dtype=float)
    q2_c = np.array([b["q2"] for b in bins], dtype=float)
    amp = np.array([b["amp"] for b in bins], dtype=float)
    err = np.array([b["err"] for b in bins], dtype=float)
    if x_c.size < 4:
        raise ValueError("fold_shape_fit needs at least 4 bins for a "
                         "3-parameter shape, got %d" % x_c.size)
    if x_ref is None:
        x_ref = float(np.exp(np.mean(np.log(x_c))))
    f1_c = resp.f1_center(x_c, q2_c)
    base_c = np.asarray(base(x_c, q2_c, f1_c), dtype=float)
    # the response MC error of each bin's fold belongs in the chi2 next to
    # the measurement error: the feed-in bins beyond the plotted range are
    # exactly the ones with few MC events, and without it they would pull
    # the shape with a weight their response does not deserve
    # sparse folding kernels, on the UNION of the cells the bins reach
    kernels = [resp.fold_kernel(m, category) for m in masks]
    mc_c = np.zeros_like(x_c)
    if mc_error:
        mc_c = np.abs(amp) * np.array(
            [resp.fold_mc_error(base, m, category)
             / max(abs(resp.fold(base, m, category, kernel=k)), 1e-300)
             for m, k in zip(masks, kernels)])
    wgt = 1.0 / np.maximum(err ** 2 + mc_c ** 2, 1e-300)

    cells = np.unique(np.concatenate([k[0] for k in kernels]))
    rows = [(np.searchsorted(cells, k[0]), k[1] / k[2]) for k in kernels]
    smp = resp.sampler
    base_u = np.asarray(base(smp.x_cells[cells], smp.q2_cells[cells],
                             smp.tables["f1"][cells]), dtype=float)
    lx_u = np.log(smp.x_cells[cells] / x_ref)
    l1x_u = np.log(np.maximum(1.0 - smp.x_cells[cells], 1e-12))
    lx_c = np.log(x_c / x_ref)
    l1x_c = np.log(np.maximum(1.0 - x_c, 1e-12))

    def folds(theta):
        g = base_u * np.exp(theta[0] * lx_u + theta[1] * l1x_u)
        return np.array([float(v.dot(g[idx])) for idx, v in rows])

    def k_of(theta):
        f = folds(theta)
        d = base_c * np.exp(theta[0] * lx_c + theta[1] * l1x_c)
        return d / f, f

    # grid + zoom: deterministic, and the box is the conditioning guard
    lo = np.array([bounds[0][0], bounds[1][0]], dtype=float)
    hi = np.array([bounds[0][1], bounds[1][1]], dtype=float)
    box = np.array([lo, hi])
    best = (np.inf, 0.0, np.array([0.0, 0.0]))
    for _ in range(n_zoom):
        g1 = np.linspace(box[0, 0], box[1, 0], n_grid)
        g2 = np.linspace(box[0, 1], box[1, 1], n_grid)
        for c1 in g1:
            for c2 in g2:
                th = np.array([c1, c2])
                chi2, norm = _profiled_chi2(folds(th), amp, wgt)
                if chi2 < best[0]:
                    best = (chi2, norm, th)
        half = 2.0 * np.array([g1[1] - g1[0], g2[1] - g2[0]])
        box = np.array([np.maximum(best[2] - half, lo),
                        np.minimum(best[2] + half, hi)])
    chi2, norm, theta = best
    at_bound = bool(np.any(np.abs(theta - lo) < 1e-9)
                    or np.any(np.abs(theta - hi) < 1e-9))

    # covariance of the shape parameters from the profiled-chi2 curvature
    # (chi2 = chi2_min + dtheta^T H dtheta / 2 -> C = 2 H^-1), and the
    # gradient of K, by the same finite difference
    step = 0.02 * (hi - lo)
    f0 = _profiled_chi2(folds(theta), amp, wgt)[0]
    hess = np.zeros((2, 2))
    kk0, fold0 = k_of(theta)
    grad = np.zeros((2, kk0.size))
    for i in range(2):
        e_i = np.zeros(2)
        e_i[i] = step[i]
        fp = _profiled_chi2(folds(theta + e_i), amp, wgt)[0]
        fm = _profiled_chi2(folds(theta - e_i), amp, wgt)[0]
        hess[i, i] = (fp - 2.0 * f0 + fm) / step[i] ** 2
        grad[i] = (k_of(theta + e_i)[0] - k_of(theta - e_i)[0]) \
            / (2.0 * step[i])
    e_01 = np.array(step)
    fpp = _profiled_chi2(folds(theta + e_01), amp, wgt)[0]
    fmm = _profiled_chi2(folds(theta - e_01), amp, wgt)[0]
    hess[0, 1] = hess[1, 0] = (fpp + fmm - 2.0 * f0
                               - hess[0, 0] * step[0] ** 2
                               - hess[1, 1] * step[1] ** 2) / (2.0 * step[0]
                                                               * step[1])
    if np.all(np.linalg.eigvalsh(hess) > 0.0):
        cov = 2.0 * np.linalg.inv(hess)
        k_fit = np.sqrt(np.maximum(
            np.einsum("ib,ij,jb->b", grad, cov, grad), 0.0)) / np.abs(kk0)
    else:                       # flat or saddle: no defensible shape error
        cov = np.full((2, 2), np.nan)
        k_fit = np.zeros_like(kk0)
    shape = tilted_shape(base, x_ref, (norm, theta[0], theta[1]))
    k_mc = np.zeros_like(kk0)
    if mc_error:
        k_mc = np.array([resp.fold_mc_error(shape, m, category)
                         for m in masks]) / np.abs(fold0 * norm)
    # the residual the tilt cannot absorb: refit from every other shape
    # of the family on the same data (K is what the correction multiplies,
    # so the spread of K is the spread of the corrected Delta)
    k_prior = np.zeros_like(kk0)
    alt_chi2 = []
    for alt in alt_bases:
        alt_fit = fold_shape_fit(resp, category, bins, alt, x_ref=x_ref,
                                 bounds=bounds, n_grid=n_grid,
                                 n_zoom=n_zoom, mc_error=mc_error)
        k_prior = np.maximum(k_prior, np.abs(alt_fit["k"] / kk0 - 1.0))
        alt_chi2.append(alt_fit["chi2"])
    k_rel_all = np.sqrt(k_fit ** 2 + k_mc ** 2 + k_prior ** 2)
    return {"params": (norm, float(theta[0]), float(theta[1])),
            "shape": shape, "x_ref": x_ref, "k": kk0,
            "k_fit_rel": k_fit, "k_mc_rel": k_mc, "k_prior_rel": k_prior,
            "k_rel": np.hypot(k_fit, k_mc),
            "k_err": np.abs(kk0) * np.hypot(k_fit, k_mc),
            "k_rel_all": k_rel_all, "k_err_all": np.abs(kk0) * k_rel_all,
            "alt_chi2": alt_chi2,
            "fold": fold0 * norm, "cov": cov, "chi2": chi2,
            "ndof": int(x_c.size - 3), "at_bound": at_bound}


# --- coherent ------------------------------------------------------------

class CoherentResponse:
    """Roman-Pot emulated response of the coherent recoil spectrum."""

    def __init__(self, scenario, config, sigma_theta, aspect=1.0,
                 shape="rectangle", n_mc=400000, x_pom_range=(1e-3, 1e-2),
                 t_max=0.5, phi_s=np.pi / 2.0, n_sigma=10.0, rng=None,
                 cut_scale_xy=(1.0, 1.0), t_floor=0.0,
                 cut_theta_xy=None):
        """`sigma_theta` is the horizontal RMS divergence [rad], or a pair
        (sigma_x, sigma_y) -- the per-configuration values of
        reco.sigma_theta_for / reco.tagging_optics_point are anisotropic
        (2026-08-28); with a scalar, `aspect` = sigma_theta_y /
        sigma_theta_x (HERA's proton beam had 45 vs 100 MeV horizontal vs
        vertical pT spread at the IP, ZEUS NPB 816:1).  `cut_scale_xy`
        scales the cutout half-widths (n_sigma sigma_x, n_sigma sigma_y).
        `cut_theta_xy` adds the pots' measured GEOMETRIC aperture in
        absolute angle and takes the larger of the two per axis
        (reco.RP_APERTURE_MEASURED); at the Yellow Report optics it
        dominates the envelope at every configuration and inverts its
        aspect, under a tagging optics with pots that follow the envelope
        it is the envelope that binds.

        `t_floor` [GeV^2] importance-samples the recoil above |t| = t_floor
        instead of from |t| = 0.  The spectrum is exponential, so a shifted
        exponential has a CONSTANT likelihood ratio exp(-B t_floor): the
        only change is the per-event weight, and a tight envelope stops
        being unsimulable.  At the 0.60 GeV envelope of the top-energy
        configuration the plain sampler leaves ZERO accepted recoils in
        6e5 draws (exp(-B pT^2) = 1.5e-8), which is why the WP5 error
        curve could not be produced."""
        rng = rng or np.random.default_rng(20260824)
        self.scenario, self.config = scenario, config
        if np.ndim(sigma_theta) == 1:
            sigma_theta, aspect = (float(sigma_theta[0]),
                                   float(sigma_theta[1]) / float(sigma_theta[0]))
        self.sigma_theta, self.aspect, self.shape = sigma_theta, aspect, shape
        self.n_sigma, self.phi_s = n_sigma, phi_s
        _k, p_ion = reco.beam_fourvectors(config)
        t_floor = float(t_floor)
        if t_floor > 0.0:
            t = -(t_floor + rng.exponential(1.0 / scenario.slope_b,
                                            size=n_mc))
        else:
            t = -scenario.sample_t(n_mc, rng, t_max=t_max)
        self.t_floor = t_floor
        self._t_weight = float(np.exp(-scenario.slope_b * t_floor))
        lo, hi = np.log10(x_pom_range[0]), np.log10(x_pom_range[1])
        x_pom = 10 ** rng.uniform(lo, hi, size=n_mc)
        t_min = reco.mdot(p_ion, p_ion) * x_pom ** 2 / (1.0 - x_pom)
        ok = -t > t_min
        t, x_pom = t[ok], x_pom[ok]
        phi_t = rng.uniform(0.0, 2.0 * np.pi, size=t.size)
        pp = reco.recoil_fourvector(t, phi_t, x_pom, p_ion)
        # the PRE-CUT measurement is kept so that a perturbed cutout is a
        # re-masking of the SAME recoils (with_cut): rebuilding the
        # response would cancel the common mode only to MC statistics,
        # which is far larger than the perturbation of interest.
        self._m = reco.rp_measure(pp, p_ion,
                                  (sigma_theta, sigma_theta * aspect),
                                  n_sigma=n_sigma, shape=shape, rng=rng,
                                  cut_scale_xy=(1.0, 1.0))
        self._t, self._phi_t, self._x_pom = t, phi_t, x_pom
        self.n_produced_mc = n_mc
        self.pt_cut = reco.tag_pt_cut(sigma_theta, config.ion_momentum_per_nucleon,
                                      a_beam=config.ion.A, n_sigma=n_sigma)
        self.cut_theta_xy = (None if cut_theta_xy is None
                             else tuple(float(v) for v in cut_theta_xy))
        self._apply_cut(cut_scale_xy)

    def _apply_cut(self, cut_scale_xy, eff_scale_xy=(1.0, 1.0)):
        """(Re)select the tagged sample for a cutout of half-widths
        (n_sigma sigma_x, n_sigma sigma_y) * cut_scale_xy, or the measured
        aperture where that is larger, from the already-smeared angle
        pair.  `eff_scale_xy` scales the EFFECTIVE (binding) half-widths
        after that comparison -- the per-fill perturbation of `with_cut`,
        which therefore acts whichever constraint binds (2026-08-28; scaling
        the envelope alone was a no-op wherever the aperture bound)."""
        self.cut_scale_xy = tuple(float(v) for v in cut_scale_xy)
        m, n_mc = self._m, self.n_produced_mc
        sx, sy = self.sigma_theta, self.sigma_theta * self.aspect
        cx = self.n_sigma * sx * self.cut_scale_xy[0]
        cy = self.n_sigma * sy * self.cut_scale_xy[1]
        if self.cut_theta_xy is not None:
            cx = max(cx, self.cut_theta_xy[0])
            cy = max(cy, self.cut_theta_xy[1])
        cx, cy = cx * float(eff_scale_xy[0]), cy * float(eff_scale_xy[1])
        self.cut_theta_eff = (cx, cy)
        thx, thy = m["theta_x"], m["theta_y"]
        if self.shape == "rectangle":
            acc = (np.abs(thx) > cx) | (np.abs(thy) > cy)
        else:
            acc = (thx / cx) ** 2 + (thy / cy) ** 2 > 1.0
        p_beam = self.config.ion.A * self.config.ion_momentum_per_nucleon
        self.cut_pt_xy = (p_beam * cx, p_beam * cy)
        self.acceptance = float(acc.sum()) * self._t_weight / n_mc
        self.t_true = -self._t[acc]
        self.t_reco = -m["t_reco"][acc]
        self.beta_true = np.mod(self._phi_t[acc] - self.phi_s, 2.0 * np.pi)
        self.beta_reco = np.mod(m["phi_t"][acc] - self.phi_s, 2.0 * np.pi)
        self.x_pom = self._x_pom[acc]
        # per produced recoil; the importance-sampling likelihood ratio
        # is constant for a shifted exponential
        self.w = np.full(self.t_true.size, self._t_weight / n_mc)
        return self

    def with_cut(self, cut_scale_xy=None, eff_scale_xy=(1.0, 1.0)):
        """A view of the SAME recoils behind a perturbed cutout -- the
        fill-dependent acceptance the spin-state ratio cannot cancel
        (code review F1, coherent half).  `cut_scale_xy` rescales the
        envelope (the pre-2026-08-28 interface, inert wherever the measured
        aperture binds); `eff_scale_xy` rescales the binding half-widths
        themselves, envelope or aperture.  Cheap: no resampling."""
        if cut_scale_xy is None:
            cut_scale_xy = self.cut_scale_xy
        return copy.copy(self)._apply_cut(cut_scale_xy, eff_scale_xy)

    def t_bin_fraction(self, tlo, thi):
        sel = (self.t_reco >= tlo) & (self.t_reco < thi)
        return float(self.w[sel].sum())

    def expected_counts_2d(self, n_produced, pzz_list, lumi_fractions,
                           tlo, thi, alpha_edges, beta_edges, a_e, a_t_func,
                           a_m=0.0, u1=0.0, u2=0.0, kappa=0.0,
                           a_e_s=0.0, a_t_s_func=None, a_m_s=0.0):
        """Exact expected (F, Ka, Kb) counts for recoils reconstructed in
        [tlo, thi): per event the alpha integrals of
        1 + u1 cos(a-b) + u2 cos 2(a-b) + P_f[kappa + a_e cos 2a
        + a_t(t) cos 2b + a_m cos(a+b)] (TRUE angles/t inside, RECO beta
        for the binning) are analytic.  `n_produced` = coherent recoils
        produced in the sample (before tagging).

        `a_e_s`, `a_t_s_func`, `a_m_s` inject the parity-FORBIDDEN sin
        partners, which is how an azimuthal misalignment is emulated: a
        spin-axis error delta turns a_t cos 2beta into
        a_t cos 2d cos 2beta + a_t sin 2d sin 2beta (rotating phi_s in
        the constructor injects nothing -- it only relabels beta)."""
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
        if a_e_s or a_m_s or a_t_s_func is not None:
            i2s = (np.cos(2 * lo) - np.cos(2 * hi)) / (4.0 * np.pi)
            j1p_s = (np.cos(lo + b) - np.cos(hi + b)) / (2.0 * np.pi)
            a_t_s = (0.0 if a_t_s_func is None
                     else np.asarray(a_t_s_func(tt), dtype=float)[:, None])
            tens = (tens + a_e_s * i2s + a_t_s * np.sin(2.0 * b) * i0
                    + a_m_s * j1p_s)
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
        basis of reco.basis_2d): {"c1","s1","c2","s2","c2t","s2t"}, each
        (Kb,).

        "c2t" = <g(t) cos 2beta_true>_k is the template basis of the
        t-dependent deformation term a_t(t) = a_t(t_ref) g(t/t_ref)
        (default g = t/t_ref, the linear model of plans/06 SS6.4b): the
        cutout correlates beta with |t| inside a reco t bin (the blind
        directions hold the larger-|t| recoils), so the coefficient must
        be fitted against the template, not against <cos 2beta> alone.
        Within one reco t bin this is a ONE-coefficient fit of the
        assumed shape (a_t strictly proportional to |t|); the
        two-component decomposition a_2(t) = c_def |t| + a_g of plans/06
        is what the bin-to-bin t dependence of the fitted a_t(t_ref)
        gives, together with the flat a_e of the electron azimuth -- not
        something this basis fits inside a bin.  The fitted coefficient is
        a_t at t_ref (the bin's rate-weighted mean true |t| by default),
        and the basis itself carries the response's Monte-Carlo
        statistics: keep n_mc large enough that the in-bin means are
        known better than the data resolve them (6e6 for the ten-year
        Table 3 of Report 2)."""
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
                                   np.cos(2 * bc) * d2),
                                  ("s2t", g * np.sin(2 * bt),
                                   np.sin(2 * bc) * d2)):
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
                     n_beta=24, rng=None, poisson=True, responses=None,
                     lumi_assumed=None, u_coeffs_assumed=None,
                     with_sin=False, a_e_s=0.0, a_t_s_func=None, a_m_s=0.0):
    """One two-azimuth pseudo-measurement of a reco t bin: expected 2-D
    counts per spin state -> Poisson -> reco.harmonic_ratio_fit_2d with
    the acceptance-weighted beta basis from the same response
    (cresp.basis_means).

    Three arguments separate what the analysis ASSUMES from what is true,
    so that the systematics the ratio does not cancel can be exercised
    (code review F1, F6):
      `responses`      one CoherentResponse per fill for GENERATION (e.g.
                       cresp.with_cut(...) with a perturbed envelope);
                       the basis and the truth reference stay on `cresp`,
                       which is what a real analysis has.
      `lumi_assumed`   luminosity fractions the analysis believes.
      `u_coeffs_assumed`  (u1, u2) the analysis subtracts, if not the
                       generated ones.
    """
    rng = rng or np.random.default_rng(20260824)
    ae = np.linspace(0.0, 2.0 * np.pi, n_alpha + 1)
    be = np.linspace(0.0, 2.0 * np.pi, n_beta + 1)
    pzz = [float(c.moments()[1]) for c in plan.categories]
    frac = [c.lumi_fraction for c in plan.categories]
    gens = list(responses) if responses is not None else [cresp] * len(pzz)
    if len(gens) != len(pzz):
        raise ValueError("responses must have one entry per fill")
    mu = np.concatenate([
        g.expected_counts_2d(n_produced, [pf], [lf], tlo, thi, ae, be,
                             a_e, a_t_func, a_m=a_m, u1=u1, u2=u2,
                             kappa=kappa, a_e_s=a_e_s,
                             a_t_s_func=a_t_s_func, a_m_s=a_m_s)
        for g, pf, lf in zip(gens, pzz, frac)], axis=0)
    counts = rng.poisson(mu) if poisson else mu
    bm = cresp.basis_means(tlo, thi, be)
    lum = frac if lumi_assumed is None else list(lumi_assumed)
    u_ass = (u1, u2) if u_coeffs_assumed is None else tuple(u_coeffs_assumed)
    fit = reco.harmonic_ratio_fit_2d(counts, lum, pzz, ae, be,
                                     u_coeffs=u_ass, beta_means=bm,
                                     with_sin=with_sin)
    fit["beta_means"] = bm
    fit.update({"counts": counts, "expected": mu, "alpha_edges": ae,
                "beta_edges": be, "n": float(mu.sum()),
                "truth": cresp.truth_reference(tlo, thi, a_e, a_t_func, a_m)})
    return fit
