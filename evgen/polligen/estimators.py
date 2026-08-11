"""Analysis-side estimators for spin-labeled pseudo-experiments.

These are the counting/fit estimators an experiment would apply; the
Step-5.A closure tests check that their spreads reproduce the analytic
error formulas of `polli_fastsim.asymmetries` and that their pulls are
unbiased -- including with luminosity-corrected yields when the run plan
carries relative-luminosity offsets.
"""

import numpy as np


def yields(counts, lumis=None):
    """Luminosity-normalized yields.  `counts` and `lumis` are sequences
    aligned by category; lumis=None means the naive equal-share assumption."""
    counts = np.asarray(counts, dtype=float)
    if lumis is None:
        return counts
    lumis = np.asarray(lumis, dtype=float)
    return counts / lumis * lumis.mean()


def apar_flip(n_plus, n_minus, pe, pz, l_plus=None, l_minus=None):
    """Helicity-flip estimator: (y+ - y-)/(y+ + y-) / (pe pz)."""
    lumis = None if l_plus is None else (l_plus, l_minus)
    y = yields((n_plus, n_minus), lumis)
    return (y[0] - y[1]) / (y[0] + y[1]) / (pe * pz)


def azz_thirds(n_plus, n_minus, n_zero, pzz, lumis=None):
    """Tensor thirds estimator (y+ + y- - 2 y0)/(y+ + y- + y0) / pzz for the
    (+pzz, +pzz, -2 pzz) fill pattern of bookkeeping.tensor_thirds_plan."""
    y = yields((n_plus, n_minus, n_zero), lumis)
    return (y[0] + y[1] - 2.0 * y[2]) / (y[0] + y[1] + y[2]) / pzz


def cos2phi_moment(phi_prime, pzz):
    """Moment estimator 2 <cos 2phi'> / pzz -- uniform acceptance only."""
    return 2.0 * np.mean(np.cos(2.0 * np.asarray(phi_prime))) / pzz


def cos2phi_fit_binned(counts, edges, pzz, acceptance=None, nsub=9):
    """Binned least-squares amplitude of N(phi') = C (1 + A cos 2phi') from
    a precomputed phi' histogram (counts on `edges`), robust to phi
    acceptance holes.

    This is the histogram-level core of `cos2phi_fit`; use it directly for
    full-luminosity projections where the per-phi-bin counts are drawn
    Poisson from expected yields (binned estimators see statistics
    identical to event-level sampling, without 1e8-event arrays).

    `acceptance`: optional boolean callable acc(phi) marking live phi values.
    Each bin is probed at `nsub` points across its full width and dropped
    unless ALL are live -- a bin merely straddling a hole edge would enter
    the fit at full weight with depleted counts and bias the amplitude
    (by more than the gluonometry signal for unlucky hole placements).
    Returns A / pzz, with the finite-bin-width dilution sin(w)/w
    (w = bin half-width in 2phi) corrected.

    Preconditions: uniform bins; the counts follow C(1 + A cos 2phi') --
    a residual cos phi' modulation is orthogonal only at full coverage
    and would leak into A through acceptance holes.  The gluonometry
    fills satisfy this by design (m-symmetric tensor populations give
    a1_eff = 0 identically; see sample.effective_modulation).
    """
    counts = np.asarray(counts, dtype=float)
    edges = np.asarray(edges, dtype=float)
    nbins = counts.size
    centers = 0.5 * (edges[:-1] + edges[1:])
    use = np.ones(nbins, dtype=bool)
    if acceptance is not None:
        frac = np.linspace(0.0, 1.0, nsub)
        probes = edges[:-1, None] + frac[None, :] * (edges[1] - edges[0])
        live = np.asarray(acceptance(probes.ravel()),
                          dtype=bool).reshape(probes.shape)
        use &= live.all(axis=1)
    if use.sum() < 3:
        raise ValueError(
            "cos2phi_fit_binned: only %d live phi bins after the "
            "acceptance probe -- amplitude unidentifiable" % int(use.sum()))
    c2 = np.cos(2.0 * centers[use])
    n = counts[use]
    # linear LSQ for N = C + B*cos2phi  (B = C*A)
    design = np.vstack([np.ones(c2.size), c2]).T
    coef, *_ = np.linalg.lstsq(design, n, rcond=None)
    amp = coef[1] / coef[0]
    w = 0.5 * (edges[1] - edges[0])  # half-width; dilution in 2phi
    dilution = np.sin(2.0 * w) / (2.0 * w)
    return amp / dilution / pzz


def cos2phi_fit_err(n, pzz, nbins=24):
    """Statistical error of the cos2phi_fit_binned amplitude on n events
    over nbins uniform full-coverage bins: sqrt(2/n)/pzz inflated by the
    same 1/dilution the fit divides out of the amplitude (audit
    2026-08-10: +1.1% at 24 bins, +11% at 8 -- the bare sqrt(2/n)/pzz of
    asymmetries.err_cos2phi_amplitude is the unbinned/analytic-FOM
    convention and understates the binned-fit spread)."""
    n = np.maximum(np.asarray(n, dtype=float), 1e-12)
    w = np.pi / nbins  # bin half-width; dilution in 2*phi
    dilution = np.sin(2.0 * w) / (2.0 * w)
    return np.sqrt(2.0 / n) / pzz / dilution


def cos2phi_fit(phi_prime, pzz, nbins=36, acceptance=None, nsub=9):
    """Event-level wrapper of `cos2phi_fit_binned`: histogram phi' on
    `nbins` uniform bins over [0, 2pi) and fit the cos 2phi' amplitude."""
    phi_prime = np.mod(np.asarray(phi_prime, dtype=float), 2.0 * np.pi)
    edges = np.linspace(0.0, 2.0 * np.pi, nbins + 1)
    counts, _ = np.histogram(phi_prime, bins=edges)
    return cos2phi_fit_binned(counts, edges, pzz, acceptance=acceptance,
                              nsub=nsub)


def binned_counts(events, x_edges, q2_edges):
    """Event counts on an (x, Q2) analysis grid (fom.py log_grid edges)."""
    h, _, _ = np.histogram2d(events["x"], events["q2"],
                             bins=(x_edges, q2_edges))
    return h


def pulls(estimates, truth, expected_err):
    """(estimate - truth)/expected_err for arrays of trial estimates."""
    return (np.asarray(estimates, dtype=float) - truth) / expected_err
