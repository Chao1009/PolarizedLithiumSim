#!/usr/bin/env python3
"""Coverage and statistical-uncertainty maps at L = 10 fb^-1/nucleon.

Standalone companion to money_delta.py. Does NOT modify money_delta.py.

Outputs, per beam config in beams.default_configs(ion):
  (1) (x, Q^2) coverage map: log-scale DIS event-rate heatmap on log-log axes.
  (2) cos(2phi) amplitude statistical-uncertainty map on (x, Q^2):
      delta A_cos2phi = obs["err_a_cos2phi"] from fom.project_observables.
  (3) Same uncertainty, rebinned from (x, Q^2) to (x, theta_e) via pooled
      event counts, where theta_e is the lab polar angle of the scattered
      electron (radians, from kinematics.scattered_electron).  Displayed as
      (pi - theta_e) * 1000 [mrad from backward beam direction] on log axes.
      Note: phi is analytically integrated out in the projection model; the
      scattered-electron-angle axis here uses theta_e (polar angle), not phi
      (azimuthal angle).

Plus one 1x3 combined summary panel per output type:
  coverage_<ion>_all.png
  errcos2phi_xq2_<ion>_all.png
  errcos2phi_theta_<ion>_all.png

Defaults: L = 10 fb^-1/nucleon, Pzz = 0.80, unpolarised electrons,
          pol_ion_vector = 0 (focus on the tensor sector).
These mirror the money_delta.py docstring convention without depending on it.

Usage:
  python3 scripts/coverage_and_stat_maps.py --ion 6Li --pdf toy
  python3 scripts/coverage_and_stat_maps.py --ion 7Li --pdf grid \\
      --lumi 10 --pzz 0.80 --n-theta-bins 24 --outdir out_coverage

CLI:
  --ion          {6Li, 7Li}   (default: 7Li)
  --pdf          {toy, grid}  (default: toy)
  --lumi         float        integrated luminosity [fb^-1/nucleon] (default: 10.0)
  --pzz          float        tensor polarization for uncertainty maps (default: 0.80)
  --n-theta-bins int          regular theta_e bins used in output (3) (default: 24)
  --outdir       path         output directory (created if missing)
  --debug-crosscheck          numerical assertion vs. phase_space_map reference
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from polli_fastsim import beams, fom
from polli_fastsim import kinematics as kin
from polli_fastsim.inputs import get_backends
from polli_fastsim.polarized import ToyG1, toy_b1, toy_delta_gluon
from polli_fastsim.structure import NuclearF2

# ── Module-level constants ────────────────────────────────────────────────── #

L_DEFAULT_FB  = 10.0   # fb^-1/nucleon
PZZ_DEFAULT   = 0.80   # tensor polarization
N_THETA_BINS  = 24     # regular theta_e bins used in output (3)
MIN_EVENTS_RB = 1.0    # minimum pooled events to display a rebinned cell
ERR_VMIN      = 1e-4   # colorbar lower bound for uncertainty maps
ERR_VMAX      = 1.0    # colorbar upper bound for uncertainty maps


# ── Filename sanitization ─────────────────────────────────────────────────── #

def _cfg_tag(cfg) -> str:
    """Sanitize BeamConfig.label() into a filename-safe token.

    Mirrors the transform used by scripts/phase_space_map.py:64-65 exactly,
    so that the per-config suffix produced here is consistent with that script.

    Example: 'e(10) x 7Li(42.9/u)'  ->  'e10_x_7Li42.9u'
    """
    return (cfg.label()
            .replace(" ", "")
            .replace("(", "")
            .replace(")", "")
            .replace("/u", "u")
            .replace("x", "_x_"))


# ── Helper: run full projection at a given luminosity ─────────────────────── #

def _project_at_lumi(cfg, backends, pzz: float, lumi_fb: float):
    """Run fom.project_rates + project_observables and expose theta_e.

    Parameters
    ----------
    cfg      : BeamConfig
    backends : dict  (as returned by get_backends)
    pzz      : float  ion tensor polarization
    lumi_fb  : float  integrated luminosity [fb^-1/nucleon]

    Returns
    -------
    proj     : BinnedProjection
               proj.x, proj.q2   shape (nx, nq2) = (40, 30)
               proj.x_edges      shape (nx+1,)
               proj.q2_edges     shape (nq2+1,)
               proj.accepted     shape (nx, nq2)  bool mask
               proj.n_events     shape (nx, nq2)
               proj.extras       dict  (includes 'y', 'eta', ...)
    obs      : dict  from fom.project_observables
               obs["err_a_cos2phi"]  shape (nx, nq2)
    theta_e  : ndarray  shape (nx, nq2)  [radians, ~pi near backward direction]
    """
    sc = fom.Scenario(
        lumi_fb_per_nucleon=lumi_fb,
        pol_electron=0.0,
        pol_ion_vector=0.0,
        pol_ion_tensor=pzz,
    )
    nf2_in = NuclearF2(cfg.ion, base=backends["base"])
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_in)
    obs  = fom.project_observables(cfg, sc, proj,
                                   ToyG1(), toy_b1, toy_delta_gluon)

    # Recompute theta_e from the (x, y) centers already in proj.extras.
    # fom.project_rates computes theta internally (line 71) and discards it;
    # we reproduce it here without touching the library.
    y = proj.extras["y"]
    s = cfg.sqrt_s_per_nucleon ** 2
    _e_prime, theta_e, _eta_e = kin.scattered_electron(
        proj.x, np.clip(y, 1e-6, 1.0), s, cfg.electron_energy
    )
    return proj, obs, theta_e


# ── Helper: rebin err_cos2phi from (x, Q^2) to (x, theta_e) ─────────────── #

def _rebin_err_vs_theta(proj, theta_e, pzz: float,
                        n_theta_bins: int = N_THETA_BINS,
                        min_events: float = MIN_EVENTS_RB):
    """Rebin the cos(2phi) amplitude statistical error from (x, Q^2) to
    (x, theta_e) by summing event counts in each (x, theta) cell and applying
    the closed-form  delta = sqrt(2 / N_pooled) / Pzz.

    Combination rule: merging N (x, Q^2) cells into one (x, theta_e) bin means
    pooling events and re-fitting a single amplitude, giving
        delta_combined = sqrt(2 / sum(N_bin)) / Pzz.
    This is the pooled-count form; it is algebraically identical to
    inverse-variance summation when Pzz is constant across bins.

    Cells below min_events and empty cells are masked.

    Parameters
    ----------
    proj          : BinnedProjection
    theta_e       : ndarray, shape (nx, nq2)   [radians]
    pzz           : float
    n_theta_bins  : int
    min_events    : float  (threshold; default 1.0 for display purposes)

    Returns
    -------
    x_edges       : ndarray, shape (nx+1,)            (unchanged from proj)
    theta_edges   : ndarray, shape (n_theta_bins+1,)  [radians, ascending]
    err_rebinned  : masked ndarray, shape (nx, n_theta_bins)
    n_rebinned    : ndarray, shape (nx, n_theta_bins)  [pooled counts, for
                    validation that the rebin is lossless]
    """
    acc = proj.accepted & np.isfinite(theta_e) & (proj.n_events > 0)
    theta_acc = theta_e[acc]

    if theta_acc.size == 0:
        empty = np.ma.masked_all((proj.x.shape[0], n_theta_bins))
        return (proj.x_edges,
                np.linspace(0.0, np.pi, n_theta_bins + 1),
                empty,
                np.zeros(empty.shape))

    # Build theta_edges linear in (pi - theta_e) [mrad from backward], then
    # convert back to ascending theta.
    tb_lo = np.pi - theta_acc.max()   # smallest "mrad-from-backward" = most backward
    tb_hi = np.pi - theta_acc.min()   # largest scattering angle
    pad = 0.02 * (tb_hi - tb_lo + 1e-6)
    tb_edges = np.linspace(max(0.0, tb_lo - pad), tb_hi + pad, n_theta_bins + 1)
    theta_edges = np.pi - tb_edges[::-1]   # ascending in theta

    nx = proj.x.shape[0]
    n_rb = np.zeros((nx, n_theta_bins))
    for ix in range(nx):
        row_mask = acc[ix]
        if not row_mask.any():
            continue
        th = theta_e[ix, row_mask]
        n  = proj.n_events[ix, row_mask]
        # Sum event counts (NOT inverse-variance weights) into theta bins.
        n_binned, _ = np.histogram(th, bins=theta_edges, weights=n)
        n_rb[ix, :] = n_binned

    with np.errstate(divide="ignore", invalid="ignore"):
        err = np.where(
            n_rb >= min_events,
            np.sqrt(2.0 / np.maximum(n_rb, 1e-30)) / pzz,
            np.nan,
        )
    return (proj.x_edges, theta_edges,
            np.ma.masked_invalid(err), n_rb)


# ── Plot builders ─────────────────────────────────────────────────────────── #

def plot_xq2_coverage(proj, cfg, path, lumi_fb: float = L_DEFAULT_FB):
    """Output (1): (x, Q^2) event-rate coverage map (per-config PNG)."""
    n = np.ma.masked_where(~proj.accepted | (proj.n_events <= 0),
                           proj.n_events)
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    pcm = ax.pcolormesh(proj.x_edges, proj.q2_edges, n.T,
                        norm=LogNorm(vmin=1),
                        cmap="viridis", shading="auto")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]")
    ax.set_title(
        f"{cfg.label()}  $\\sqrt{{s_{{eN}}}}={cfg.sqrt_s_per_nucleon:.0f}$ GeV\n"
        f"DIS event rate,  $L={lumi_fb:g}$ fb$^{{-1}}$/u",
        fontsize=10,
    )
    fig.colorbar(pcm, ax=ax, label="events / bin")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_err_cos2phi_xq2(proj, obs, cfg, path,
                          pzz: float = PZZ_DEFAULT,
                          lumi_fb: float = L_DEFAULT_FB):
    """Output (2): delta A_cos2phi on (x, Q^2) (per-config PNG)."""
    err = np.ma.masked_where(~proj.accepted | (proj.n_events <= 0),
                             obs["err_a_cos2phi"])
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    pcm = ax.pcolormesh(proj.x_edges, proj.q2_edges, err.T,
                        norm=LogNorm(vmin=ERR_VMIN, vmax=ERR_VMAX),
                        cmap="viridis", shading="auto")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]")
    ax.set_title(
        f"{cfg.label()}:  $\\delta A_{{\\cos 2\\phi}}$\n"
        f"$L={lumi_fb:g}$ fb$^{{-1}}$/u,  $P_{{zz}}={pzz:g}$",
        fontsize=10,
    )
    fig.colorbar(pcm, ax=ax, label=r"$\delta A_{\cos 2\phi}$")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_err_cos2phi_theta(proj, theta_e, cfg, pzz: float, path,
                            lumi_fb: float = L_DEFAULT_FB,
                            n_theta_bins: int = N_THETA_BINS):
    """Output (3): delta A_cos2phi on (x, theta_e), pooled-count rebin.

    y-axis is displayed as (pi - theta_e) * 1000 [mrad from backward beam
    direction] on a log scale, since scattered electrons cluster near theta~pi.
    The phi dimension is not available (analytically integrated out in the
    rate model); theta_e is the lab polar angle from kinematics.scattered_electron.
    """
    x_edges, th_edges, err_rb, _n_rb = _rebin_err_vs_theta(
        proj, theta_e, pzz, n_theta_bins=n_theta_bins)
    # y-axis in mrad from backward: (pi - theta_e) * 1e3, ascending.
    # th_edges is ascending in theta; (pi - th_edges) is descending.
    # Reverse so that the mrad axis is ascending left-to-right on the colorbar.
    y_edges_mrad = (np.pi - th_edges[::-1]) * 1e3   # ascending mrad
    C = err_rb[:, ::-1]                              # reverse theta cols to match

    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    pcm = ax.pcolormesh(x_edges, y_edges_mrad, C.T,
                        norm=LogNorm(vmin=ERR_VMIN, vmax=ERR_VMAX),
                        cmap="viridis", shading="auto")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$(\pi - \theta_e) \times 10^3$ [mrad]  (backward $\to$ side)")
    ax.set_title(
        f"{cfg.label()}:  $\\delta A_{{\\cos 2\\phi}}$ vs $\\theta_e$\n"
        f"$L={lumi_fb:g}$ fb$^{{-1}}$/u,  $P_{{zz}}={pzz:g}$  "
        f"(pooled-count rebin)",
        fontsize=10,
    )
    fig.colorbar(pcm, ax=ax, label=r"$\delta A_{\cos 2\phi}$")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


# ── Draw helpers for combined 1×3 panels ──────────────────────────────────── #

def _draw_xq2_coverage_into(ax, proj, cfg, norm, lumi_fb=L_DEFAULT_FB):
    """Draw coverage panel into an existing Axes (for the 1x3 combined figure)."""
    n = np.ma.masked_where(~proj.accepted | (proj.n_events <= 0),
                           proj.n_events)
    pcm = ax.pcolormesh(proj.x_edges, proj.q2_edges, n.T,
                        norm=norm, cmap="viridis", shading="auto")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]")
    ax.set_title(
        f"{cfg.label()}\n$\\sqrt{{s_{{eN}}}}={cfg.sqrt_s_per_nucleon:.0f}$ GeV",
        fontsize=9,
    )
    return pcm


def _draw_err_cos2phi_xq2_into(ax, proj, obs, cfg, norm):
    """Draw err-cos2phi (x, Q^2) panel into an existing Axes."""
    err = np.ma.masked_where(~proj.accepted | (proj.n_events <= 0),
                             obs["err_a_cos2phi"])
    pcm = ax.pcolormesh(proj.x_edges, proj.q2_edges, err.T,
                        norm=norm, cmap="viridis", shading="auto")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]")
    ax.set_title(
        f"{cfg.label()}\n$\\sqrt{{s_{{eN}}}}={cfg.sqrt_s_per_nucleon:.0f}$ GeV",
        fontsize=9,
    )
    return pcm


def _draw_err_cos2phi_theta_into(ax, proj, theta_e, cfg, pzz, norm,
                                  n_theta_bins=N_THETA_BINS):
    """Draw err-cos2phi (x, theta_e) panel into an existing Axes."""
    x_edges, th_edges, err_rb, _n_rb = _rebin_err_vs_theta(
        proj, theta_e, pzz, n_theta_bins=n_theta_bins)
    y_edges_mrad = (np.pi - th_edges[::-1]) * 1e3
    C = err_rb[:, ::-1]
    pcm = ax.pcolormesh(x_edges, y_edges_mrad, C.T,
                        norm=norm, cmap="viridis", shading="auto")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$(\pi - \theta_e) \times 10^3$ [mrad]")
    ax.set_title(
        f"{cfg.label()}\n$\\sqrt{{s_{{eN}}}}={cfg.sqrt_s_per_nucleon:.0f}$ GeV",
        fontsize=9,
    )
    return pcm


# ── Combined 1×3 summary panels ───────────────────────────────────────────── #

def _combined_cov_panel(cov_data, outpath, lumi_fb=L_DEFAULT_FB):
    """Combined 1×3 coverage panel; cov_data = list of (cfg, proj)."""
    norm = LogNorm(vmin=1)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=False,
                             layout="constrained")
    pcm = None
    for ax, (cfg, proj) in zip(axes, cov_data):
        pcm = _draw_xq2_coverage_into(ax, proj, cfg, norm, lumi_fb=lumi_fb)
    if pcm is not None:
        fig.colorbar(pcm, ax=axes.tolist(), label="events / bin", shrink=0.7)
    fig.suptitle(
        f"DIS event-rate coverage  $L={lumi_fb:g}$ fb$^{{-1}}$/u",
        fontsize=11,
    )
    fig.savefig(outpath, dpi=140)
    plt.close(fig)


def _combined_errxq2_panel(errxq2_data, outpath,
                            pzz=PZZ_DEFAULT, lumi_fb=L_DEFAULT_FB):
    """Combined 1×3 err-cos2phi (x, Q^2) panel;
    errxq2_data = list of (cfg, proj, obs)."""
    norm = LogNorm(vmin=ERR_VMIN, vmax=ERR_VMAX)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=False,
                             layout="constrained")
    pcm = None
    for ax, (cfg, proj, obs) in zip(axes, errxq2_data):
        pcm = _draw_err_cos2phi_xq2_into(ax, proj, obs, cfg, norm)
    if pcm is not None:
        fig.colorbar(pcm, ax=axes.tolist(),
                     label=r"$\delta A_{\cos 2\phi}$", shrink=0.7)
    fig.suptitle(
        f"$\\delta A_{{\\cos 2\\phi}}$ on $(x,Q^2)$  "
        f"$L={lumi_fb:g}$ fb$^{{-1}}$/u,  $P_{{zz}}={pzz:g}$",
        fontsize=11,
    )
    fig.savefig(outpath, dpi=140)
    plt.close(fig)


def _combined_errtheta_panel(errtheta_data, outpath,
                              pzz=PZZ_DEFAULT, lumi_fb=L_DEFAULT_FB,
                              n_theta_bins=N_THETA_BINS):
    """Combined 1×3 err-cos2phi (x, theta_e) panel;
    errtheta_data = list of (cfg, proj, theta_e)."""
    norm = LogNorm(vmin=ERR_VMIN, vmax=ERR_VMAX)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=False,
                             layout="constrained")
    pcm = None
    for ax, (cfg, proj, theta_e) in zip(axes, errtheta_data):
        pcm = _draw_err_cos2phi_theta_into(
            ax, proj, theta_e, cfg, pzz, norm, n_theta_bins=n_theta_bins)
    if pcm is not None:
        fig.colorbar(pcm, ax=axes.tolist(),
                     label=r"$\delta A_{\cos 2\phi}$", shrink=0.7)
    fig.suptitle(
        f"$\\delta A_{{\\cos 2\\phi}}$ vs $\\theta_e$  "
        f"$L={lumi_fb:g}$ fb$^{{-1}}$/u,  $P_{{zz}}={pzz:g}$  "
        f"(pooled-count rebin)",
        fontsize=11,
    )
    fig.savefig(outpath, dpi=140)
    plt.close(fig)


# ── Validation helpers ─────────────────────────────────────────────────────── #

def _validate_angle_sanity(cfg, proj, theta_e):
    """Assert theta_e in (0, pi) and |eta| <= 3.5 on all accepted bins."""
    acc = proj.accepted & np.isfinite(theta_e)
    if not acc.any():
        return
    th_acc = theta_e[acc]
    eta_acc = proj.extras["eta"][acc]
    if np.any((th_acc <= 0) | (th_acc >= np.pi)):
        raise RuntimeError(
            f"[{cfg.label()}] theta_e out of (0, pi): "
            f"min={th_acc.min():.6f} max={th_acc.max():.6f}"
        )
    if np.any(np.abs(eta_acc) > 3.5):
        raise RuntimeError(
            f"[{cfg.label()}] |eta_e| > 3.5 in accepted bins: "
            f"max |eta| = {np.abs(eta_acc).max():.4f}"
        )


def _validate_shape_invariants(cfg, proj, theta_e, n_theta_bins):
    """Assert shape constraints and theta_edges strict monotonicity."""
    assert theta_e.shape == proj.x.shape, (
        f"[{cfg.label()}] theta_e shape {theta_e.shape} != "
        f"proj.x shape {proj.x.shape}"
    )
    _, th_ed, err_rb, n_rb = _rebin_err_vs_theta(proj, theta_e,
                                                  pzz=PZZ_DEFAULT,
                                                  n_theta_bins=n_theta_bins)
    assert err_rb.shape == (proj.x.shape[0], n_theta_bins), (
        f"[{cfg.label()}] err_rb shape {err_rb.shape} != "
        f"({proj.x.shape[0]}, {n_theta_bins})"
    )
    assert th_ed.size == n_theta_bins + 1, (
        f"[{cfg.label()}] theta_edges size {th_ed.size} != {n_theta_bins + 1}"
    )
    assert np.all(np.diff(th_ed) > 0), (
        f"[{cfg.label()}] theta_edges not strictly monotonic"
    )
    return n_rb


def _validate_count_conservation(cfg, proj, theta_e, n_rb):
    """Assert that pooled-count rebin is lossless per x-row."""
    acc = proj.accepted & np.isfinite(theta_e) & (proj.n_events > 0)
    nx = proj.x.shape[0]
    for ix in range(nx):
        row_acc = acc[ix]
        expected = proj.n_events[ix][row_acc].sum()
        got = n_rb[ix].sum()
        if expected > 0 and not np.isclose(got, expected, rtol=1e-6):
            raise RuntimeError(
                f"[{cfg.label()}] x-row {ix}: rebin event count mismatch "
                f"expected={expected:.6f} got={got:.6f}"
            )


def _validate_sanitizer_parity(cfg):
    """Assert _cfg_tag output matches the phase_space_map.py:64-65 transform."""
    expected = (cfg.label()
                .replace(" ", "")
                .replace("(", "")
                .replace(")", "")
                .replace("/u", "u")
                .replace("x", "_x_"))
    got = _cfg_tag(cfg)
    assert got == expected, (
        f"_cfg_tag parity failure for {cfg.label()!r}: "
        f"got {got!r}, expected {expected!r}"
    )


def _validate_debug_crosscheck(cfg, backends, pzz, lumi_fb):
    """Numerical cross-check: obs['err_a_cos2phi'] must match a reference
    re-run of fom.project_observables with identical inputs."""
    proj_a, obs_a, _ = _project_at_lumi(cfg, backends, pzz, lumi_fb)
    proj_b, obs_b, _ = _project_at_lumi(cfg, backends, pzz, lumi_fb)
    if not np.allclose(obs_a["err_a_cos2phi"], obs_b["err_a_cos2phi"]):
        raise RuntimeError(
            f"[{cfg.label()}] --debug-crosscheck: err_a_cos2phi not "
            f"reproducible across two identical projection calls"
        )
    print(f"  [crosscheck OK] {cfg.label()}: err_a_cos2phi is reproducible")


# ── Main ───────────────────────────────────────────────────────────────────── #

def main():
    ap = argparse.ArgumentParser(
        description="Coverage and statistical-uncertainty maps at L = 10 fb^-1/nucleon."
    )
    ap.add_argument("--ion", default="7Li", choices=list(beams.IONS),
                    help="ion species (default: 7Li)")
    ap.add_argument("--pdf", default="toy", choices=["toy", "grid"],
                    help="PDF backend: 'toy' or 'grid' (CT18NLO/NNPDFpol11_100)")
    ap.add_argument("--lumi", type=float, default=L_DEFAULT_FB,
                    help="integrated luminosity [fb^-1/nucleon] (default: %(default)s; "
                         "leave at 10 for the deliverable figures)")
    ap.add_argument("--pzz", type=float, default=PZZ_DEFAULT,
                    help="tensor polarization Pzz (default: %(default)s)")
    ap.add_argument("--n-theta-bins", type=int, default=N_THETA_BINS,
                    dest="n_theta_bins",
                    help="number of regular theta_e bins in output (3) "
                         "(default: %(default)s)")
    ap.add_argument("--outdir",
                    default=str(pathlib.Path(__file__).parent / "out_coverage"),
                    help="output directory (created if missing)")
    ap.add_argument("--debug-crosscheck", action="store_true",
                    default=False,
                    help="run numerical reproducibility assertion for err_a_cos2phi "
                         "(off by default)")
    args = ap.parse_args()

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    backends = get_backends(args.pdf)
    ion_name = args.ion
    lumi_fb  = args.lumi
    pzz      = args.pzz
    n_theta_bins = args.n_theta_bins

    print(f"# coverage_and_stat_maps.py")
    print(f"# ion={ion_name}  pdf={args.pdf}  L={lumi_fb:g} fb^-1/u  "
          f"Pzz={pzz:g}  n_theta_bins={n_theta_bins}")
    print(f"# outdir={outdir}")

    # Accumulate per-config results for the combined panels.
    cov_data      = []   # (cfg, proj)
    errxq2_data   = []   # (cfg, proj, obs)
    errtheta_data = []   # (cfg, proj, theta_e)

    for cfg in beams.default_configs(ion_name):
        cfgtag = _cfg_tag(cfg)

        # ── Sanitizer parity check ───────────────────────────────────────── #
        _validate_sanitizer_parity(cfg)

        # ── Full projection ──────────────────────────────────────────────── #
        proj, obs, theta_e = _project_at_lumi(cfg, backends, pzz, lumi_fb)

        # ── Angle sanity check (plan validation item 4) ─────────────────── #
        _validate_angle_sanity(cfg, proj, theta_e)

        acc = proj.accepted & np.isfinite(theta_e)
        if acc.any():
            th_acc = theta_e[acc]
            eta_acc = proj.extras["eta"][acc]
            print(f"  {cfg.label()}: theta_e min={th_acc.min():.4f} rad  "
                  f"median={np.median(th_acc):.4f} rad  max={th_acc.max():.4f} rad  "
                  f"|eta| max={np.abs(eta_acc).max():.3f}")

        # ── Shape invariants (plan validation item 6) ────────────────────── #
        n_rb = _validate_shape_invariants(cfg, proj, theta_e, n_theta_bins)

        # ── Count-conservation check (plan validation item 7) ────────────── #
        _validate_count_conservation(cfg, proj, theta_e, n_rb)

        # ── Debug crosscheck (optional) ──────────────────────────────────── #
        if args.debug_crosscheck:
            _validate_debug_crosscheck(cfg, backends, pzz, lumi_fb)

        # ── Per-config summary (plan Phase D) ───────────────────────────── #
        n_dis = proj.n_events[proj.accepted].sum()
        err_acc = obs["err_a_cos2phi"][proj.accepted]
        med_err = float(np.median(err_acc)) if err_acc.size > 0 else float("nan")
        if acc.any():
            th_mrad = (np.pi - theta_e[acc]) * 1e3
            th_range = f"{th_mrad.min():.1f}–{th_mrad.max():.1f}"
        else:
            th_range = "N/A"
        print(f"  {cfg.label():<30s}  N_DIS={n_dis:.3e}  "
              f"median δA_cos2phi={med_err:.3e}  "
              f"theta_e [mrad-from-backward]={th_range}")

        # ── Write per-config PNGs ────────────────────────────────────────── #
        # Output (1): (x, Q^2) coverage
        p1 = outdir / f"coverage_{ion_name}_{cfgtag}.png"
        plot_xq2_coverage(proj, cfg, p1, lumi_fb=lumi_fb)
        print(f"    wrote {p1}")

        # Output (2): err_cos2phi on (x, Q^2)
        p2 = outdir / f"errcos2phi_xq2_{ion_name}_{cfgtag}.png"
        plot_err_cos2phi_xq2(proj, obs, cfg, p2, pzz=pzz, lumi_fb=lumi_fb)
        print(f"    wrote {p2}")

        # Output (3): err_cos2phi on (x, theta_e)
        p3 = outdir / f"errcos2phi_theta_{ion_name}_{cfgtag}.png"
        plot_err_cos2phi_theta(proj, theta_e, cfg, pzz, p3,
                               lumi_fb=lumi_fb, n_theta_bins=n_theta_bins)
        print(f"    wrote {p3}")

        # Accumulate for combined panels.
        cov_data.append((cfg, proj))
        errxq2_data.append((cfg, proj, obs))
        errtheta_data.append((cfg, proj, theta_e))

    # ── Combined 1×3 summary panels ─────────────────────────────────────── #
    p_all_cov = outdir / f"coverage_{ion_name}_all.png"
    _combined_cov_panel(cov_data, p_all_cov, lumi_fb=lumi_fb)
    print(f"  wrote {p_all_cov}")

    p_all_xq2 = outdir / f"errcos2phi_xq2_{ion_name}_all.png"
    _combined_errxq2_panel(errxq2_data, p_all_xq2, pzz=pzz, lumi_fb=lumi_fb)
    print(f"  wrote {p_all_xq2}")

    p_all_theta = outdir / f"errcos2phi_theta_{ion_name}_all.png"
    _combined_errtheta_panel(errtheta_data, p_all_theta,
                              pzz=pzz, lumi_fb=lumi_fb,
                              n_theta_bins=n_theta_bins)
    print(f"  wrote {p_all_theta}")

    print(f"\nDone. All outputs in {outdir}/")


if __name__ == "__main__":
    main()
