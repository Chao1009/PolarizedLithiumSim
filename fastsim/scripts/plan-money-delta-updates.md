# Plan: Coverage + statistical-uncertainty maps at L = 10 fb⁻¹/nucleon (new standalone script)

**Scope**: Produce three new outputs — `(x, Q²)` coverage per beam energy, `cos(2φ)` amplitude statistical-uncertainty maps on `(x, Q²)`, and a scattered-electron-angle version of the uncertainty map — all at `L = 10 fb⁻¹/u`, from a **new standalone script** in `fastsim/scripts/`. `money_delta.py` is **not** modified.

**Script path (locked)**: `fastsim/scripts/coverage_and_stat_maps.py`

Rationale for the name: it is discoverable next to `phase_space_map.py`, describes both output classes (coverage; statistical uncertainty maps), and does not overload the "money_delta" prefix that already belongs to a different figure family.

**Inputs used**:
- `docs/research-money_delta-uncertainties-and-ranges.md` (researcher artifact) — canonical trace of which per-bin uncertainties, angle variables, and beam-energy configurations the fastsim pipeline currently exposes; explicitly confirms that `err_a_cos2phi` is materialized by `fom.project_observables`, that `θ_e` is computed but discarded in `fom.project_rates`, that `phase_space_map.py` is the reference plotting pattern, and that φ is *not* a projection dimension.
- Direct source reads of `money_delta.py`, `phase_space_map.py`, and `polli_fastsim/{fom,kinematics,asymmetries,beams,structure,polarized}.py`.

**Annotation cycle**: iteration 2 — artifact ownership changed per user direction. `money_delta.py` is preserved unchanged. The previously-approved technical approach (helper + three plot builders + regularized `(x, θ_e)` rebin) is carried over verbatim into the new standalone script.

---

## Requirements Summary

**Confirmed requirements**
- Integrated luminosity fixed at `L = 10 fb⁻¹/nucleon` for all outputs.
- Three outputs, per beam configuration returned by `beams.default_configs(ion)`:
  1. `(x, Q²)` coverage plots by beam energy.
  2. Statistical uncertainty maps for the `cos(2φ)` amplitude observable on `(x, Q²)`.
  3. A scattered-electron-angle version of (2), using `θ_e` (lab polar angle).
- **`money_delta.py` must not be modified.** Deliverable is a new, self-contained script that can be run directly.
- Prefer existing library APIs; only add script-local helpers.
- Include a runnable CLI and validation commands.

**Inferred assumptions** (flag if wrong)
- "Coverage" = accepted-bin mask combined with per-bin DIS event counts (rate heatmap on log-log `x`–`Q²`), matching the existing `plot_map` style in `phase_space_map.py`.
- "Uncertainty map for `cos(2φ)` amplitude" = the per-bin statistical error `obs["err_a_cos2phi"]` already produced by `fom.project_observables` at `Pzz = scenario.pol_ion_tensor` (default 0.8).
- Physically supported scattered-electron angle = lab polar angle `θ_e` (radians) from `kinematics.scattered_electron` (`kinematics.py:32-49`, confirmed by research §3.2). Already computed inside `fom.project_rates` (line 71) and discarded as `_theta`. Displayed as `(π − θ_e) · 1000` in **mrad from backward beam direction** because the scattered electron sits near `θ ≈ π`.
- Nominal polarizations for the uncertainty maps: `Pzz = 0.80`, unpolarized electrons, `pol_ion_vector = 0` (focus on the tensor sector). This mirrors the `money_delta.py` docstring convention without depending on that file.
- Per-config PNGs plus a combined 1×3 summary panel per output type.

## Resolved defaults (previously open questions, still binding)

1. Angular variable: **`θ_e` in mrad-from-backward** (user-locked).
2. Combined 1×3 summary panels: **yes**, in addition to per-config PNGs (user-locked).
3. Angle-axis binning: **regularized `(x, θ_e)` rebinning** via pooled-count 2-D histogram, then plain rectilinear `pcolormesh`. Irregular-quad plotting rejected as fragile at acceptance boundaries.

## Open questions

None blocking. If the user later wants η_e (pseudorapidity) or `E'_e` (scattered-electron energy) surfaced in the same script, those are one-line additions — the extraction pathway is identical.

---

## Approach

The new script is a self-contained driver that:

1. Parses CLI (`--ion`, `--pdf`, `--lumi`, `--outdir`, `--n-theta-bins`, `--debug-crosscheck`).
2. Builds the requested PDF backends the same way `money_delta.py` does (call `structure.build_backends(...)` — no code duplication of physics logic, only invocation).
3. For each `cfg` in `beams.default_configs(args.ion)`:
   - Runs `_project_at_lumi(cfg, backends, pzz=0.80, lumi_fb=args.lumi)` which wraps `fom.project_rates` + `fom.project_observables` and additionally re-runs `kin.scattered_electron` to expose `θ_e` (currently discarded inside `fom.project_rates`).
   - Calls three plot builders and writes per-config PNGs.
   - Caches `(proj, obs, theta_e)` for the combined-panel pass.
4. Builds three 1×3 combined summary PNGs from the cached arrays.
5. Prints a per-config summary line (`N_DIS`, median `δA_cos2φ`, `θ_e` range in mrad-from-backward).

No changes required inside `polli_fastsim/`. No changes to `money_delta.py`, `phase_space_map.py`, or any other existing script. The new script only *reads* library APIs.

### Filename sanitization (reused verbatim from `phase_space_map.py:64-65`)

`cfg.label()` returns strings like `"e(10) x 7Li(42.9/u)"` (from `beams.py:68-72`) which contain spaces, parentheses, and `/` — unsafe in filenames. Reuse the exact transform already used by `phase_space_map.py:64-65`, defined as a module-level helper in the new script so behaviour cannot drift:

```python
def _cfg_tag(cfg):
    """Sanitize BeamConfig.label() into a filename-safe token.
    Mirrors the transform used by scripts/phase_space_map.py:64-65.
    Example:  'e(10) x 7Li(42.9/u)'  ->  'e10_x_7Li42.9u'
    """
    return (cfg.label()
            .replace(" ", "")
            .replace("(", "")
            .replace(")", "")
            .replace("/u", "u")
            .replace("x", "_x_"))
```

### Script skeleton (`fastsim/scripts/coverage_and_stat_maps.py`)

```python
"""Coverage and statistical-uncertainty maps at L = 10 fb^-1/nucleon.

Standalone companion to money_delta.py. Does NOT modify money_delta.py.

Outputs, per beam config in beams.default_configs(ion):
  (1) (x, Q^2) coverage (DIS event count, accepted mask)
  (2) cos(2phi) amplitude statistical uncertainty on (x, Q^2)
  (3) Same uncertainty rebinned onto (x, theta_e) — lab polar angle,
      displayed as (pi - theta_e) * 1e3 mrad from backward beam direction.

Plus one 1x3 combined summary panel per output type.
"""
from __future__ import annotations
import argparse, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from polli_fastsim import beams, fom, kinematics as kin
from polli_fastsim.structure import NuclearF2, build_backends
from polli_fastsim.polarized import ToyG1, toy_b1, toy_delta_gluon

L_DEFAULT_FB   = 10.0
PZZ_DEFAULT    = 0.80
N_THETA_BINS   = 24
MIN_EVENTS_RB  = 1.0
ERR_VMIN, ERR_VMAX = 1e-4, 1.0
```

### Helper: `_project_at_lumi`

```python
def _project_at_lumi(cfg, backends, pzz, lumi_fb):
    """Returns (proj, obs, theta_e). Shapes:
       proj.x, proj.q2, proj.accepted, proj.n_events : (nx, nq2)   [nx=40, nq2=30]
       proj.x_edges : (nx+1,)   proj.q2_edges : (nq2+1,)
       obs['err_a_cos2phi']    : (nx, nq2)   (masked to accepted at plot time)
       theta_e                 : (nx, nq2)   [radians, ~pi at backward]
    """
    sc = fom.Scenario(lumi_fb_per_nucleon=lumi_fb,
                      pol_electron=0.0,
                      pol_ion_vector=0.0,
                      pol_ion_tensor=pzz)
    nf2_in = NuclearF2(cfg.ion, base=backends["base"])
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_in)
    obs  = fom.project_observables(cfg, sc, proj,
                                   ToyG1(), toy_b1, toy_delta_gluon)
    y = proj.extras["y"]
    s = cfg.sqrt_s_per_nucleon ** 2
    _e_prime, theta_e, _eta_e = kin.scattered_electron(
        proj.x, np.clip(y, 1e-6, 1.0), s, cfg.electron_energy)
    return proj, obs, theta_e
```

### Output (1) — (x, Q²) coverage per config

Values = `proj.n_events`, mask = `~proj.accepted | (n_events <= 0)`, `LogNorm(vmin=1)`. Rectilinear `pcolormesh(proj.x_edges, proj.q2_edges, n_events.T)`, log-log axes. Shapes: `X=(41,)`, `Y=(31,)`, `C=(30,40)` after transpose.

### Output (2) — δA_cos2φ on (x, Q²) per config

Identical geometry to output (1). Values = `obs["err_a_cos2phi"]` shape `(nx, nq2)`, same mask, `LogNorm(vmin=1e-4, vmax=1)`. Title fixes `L=10 fb⁻¹/u`, `Pzz=0.80`. This is numerically the same figure `phase_space_map.py:78-81` produces at matching `(ion, L, Pzz)` — used as the cross-check reference in validation.

### Output (3) — δA_cos2φ on (x, θ_e), regularized rebinning

**Combination rule.** For the `cos(2φ)` amplitude, per-bin statistical error is `δA_bin = √(2/N_bin) / Pzz` (`asymmetries.err_cos2phi_amplitude`). Merging N `(x, Q²)` cells into one `(x, θ_e)` bin means pooling events and re-fitting a single amplitude → combined error is

  `δA_combined = √(2 / (Σ N_bin)) / Pzz`

i.e. **sum event counts first, then apply the closed-form**. This is algebraically identical to inverse-variance summation `w_i = N_i·Pzz²/2` when `Pzz` is constant across bins, but the pooled-count form is physically transparent and does not require constant `Pzz`.

```python
def _rebin_err_vs_theta(proj, theta_e, pzz, n_theta_bins=N_THETA_BINS,
                        min_events=MIN_EVENTS_RB):
    """Rebin the cos(2phi) amplitude statistical error from (x, Q^2) to
    (x, theta_e) by summing event counts in each rebinned (x, theta) cell
    and applying  delta = sqrt(2 / N_total) / Pzz.

    Returns:
      x_edges     (nx+1,)                   unchanged
      theta_edges (n_theta_bins+1,)         ascending in theta [rad]
      err_rb      (nx, n_theta_bins)        masked where empty / below min
      n_rb        (nx, n_theta_bins)        pooled counts, for validation
    """
    acc = proj.accepted & np.isfinite(theta_e) & (proj.n_events > 0)
    theta_acc = theta_e[acc]
    if theta_acc.size == 0:
        empty = np.ma.masked_all((proj.x.shape[0], n_theta_bins))
        return (proj.x_edges,
                np.linspace(0.0, np.pi, n_theta_bins + 1),
                empty, np.zeros_like(empty.data))

    tb_lo = np.pi - theta_acc.max()      # most backward
    tb_hi = np.pi - theta_acc.min()      # largest scattering angle
    pad   = 0.02 * (tb_hi - tb_lo + 1e-6)
    tb_edges = np.linspace(max(0.0, tb_lo - pad), tb_hi + pad,
                           n_theta_bins + 1)
    theta_edges = np.pi - tb_edges[::-1]   # ascending in theta

    nx = proj.x.shape[0]
    n_rb = np.zeros((nx, n_theta_bins))
    for ix in range(nx):
        row_mask = acc[ix]
        if not row_mask.any():
            continue
        th = theta_e[ix, row_mask]
        n  = proj.n_events[ix, row_mask]
        n_binned, _ = np.histogram(th, bins=theta_edges, weights=n)
        n_rb[ix, :] = n_binned

    with np.errstate(divide="ignore", invalid="ignore"):
        err = np.where(n_rb >= min_events,
                       np.sqrt(2.0 / np.maximum(n_rb, 1e-30)) / pzz,
                       np.nan)
    return (proj.x_edges, theta_edges,
            np.ma.masked_invalid(err), n_rb)


def plot_err_cos2phi_theta(proj, theta_e, cfg, pzz, path):
    x_edges, th_edges, err_rb, _n_rb = _rebin_err_vs_theta(
        proj, theta_e, pzz=pzz)
    y_edges_mrad = (np.pi - th_edges[::-1]) * 1e3   # ascending mrad
    C = err_rb[:, ::-1]                             # reverse cols to match
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    pcm = ax.pcolormesh(x_edges, y_edges_mrad, C.T,
                        norm=LogNorm(vmin=ERR_VMIN, vmax=ERR_VMAX),
                        cmap="viridis", shading="auto")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$\pi - \theta_e$ [mrad] (backward-scattered $\to$ side)")
    ax.set_title(f"{cfg.label()}: $\\delta A_{{\\cos 2\\phi}}$  "
                 f"$L=10$ fb$^{{-1}}$/u, $P_{{zz}}={pzz:g}$  "
                 f"(pooled-count rebin)", fontsize=10)
    fig.colorbar(pcm, ax=ax, label=r"$\delta A_{\cos 2\phi}$")
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)
```

Array-shape summary for output (3): `theta_e (40,30)`, `theta_edges (25,)`, `err_rb (40,24)`, `x_edges (41,)`, `y_edges_mrad (25,)`, `pcolormesh` receives `X=(41,), Y=(25,), C.T=(24,40)` — standard rectilinear call.

### Combined 1×3 summary panels

For each output type, `plt.subplots(1, 3, figsize=(15, 4.5))` with a **shared `LogNorm` instance** so the colorbars are directly comparable across configurations. Saved as `coverage_<ion>_all.png`, `errcos2phi_xq2_<ion>_all.png`, `errcos2phi_theta_<ion>_all.png`.

### Filenames written (in `--outdir`, default `fastsim/scripts/out_coverage/`)

Per config (using `_cfg_tag(cfg)`):
- `coverage_<ion>_<cfgtag>.png` — output (1)
- `errcos2phi_xq2_<ion>_<cfgtag>.png` — output (2)
- `errcos2phi_theta_<ion>_<cfgtag>.png` — output (3)

Combined 1×3 summary panels:
- `coverage_<ion>_all.png`
- `errcos2phi_xq2_<ion>_all.png`
- `errcos2phi_theta_<ion>_all.png`

### CLI (exact shape)

```
python3 fastsim/scripts/coverage_and_stat_maps.py \
    --ion {6Li,7Li} \
    --pdf {toy,grid} \
    [--lumi 10.0] \
    [--pzz 0.80] \
    [--n-theta-bins 24] \
    [--outdir fastsim/scripts/out_coverage] \
    [--debug-crosscheck]
```

Flag semantics:
- `--ion`, `--pdf`: mirror `money_delta.py` (`6Li`/`7Li`, `toy`/`grid`) so users have a familiar surface.
- `--lumi` (float, fb⁻¹/u, default `10.0`): the requested `L`. Left tunable purely for the linear-scaling validation step; documented as "leave at 10 for the deliverable figures".
- `--pzz` (float, default `0.80`): tensor polarization used in `Scenario` and displayed in titles.
- `--n-theta-bins` (int, default `24`): matches `N_THETA_BINS` constant.
- `--outdir`: created if missing.
- `--debug-crosscheck`: enables the numerical assertion against `phase_space_map.py` reference (see validation §2). Off by default because it requires running `phase_space_map`'s projection path a second time.

---

## Trade-offs and Notes

- **Zero risk to existing scripts.** No library edits, no edits to `money_delta.py`, `phase_space_map.py`, `money_b1.py`, `money_polemc.py`, `tagging_acceptance.py`, `validate_inputs.py`, or `tests/`. The new script only imports and reads.
- Reusing the `(x, Q²)` grid for output (3) preserves the accepted-bin mask (cuts on `y ∈ [y_min, y_max]`, `Q² ≥ Q²_min`, `|η_e| < 3.5`, `E' > E'_min`) without redoing acceptance logic — those cuts live inside `fom.project_rates`.
- Uncertainty maps depend on `Pzz`; every figure title displays it so maps are not misread as `Pzz = 1` idealisations.
- Because the script is standalone, small helper duplication with `phase_space_map.py` (the `_cfg_tag` sanitizer, the projection call pattern) is accepted deliberately. Validation §6 asserts sanitizer parity so silent drift is caught.

---

## Validation

1. **New-script smoke tests** (4 combinations):
   ```
   python3 fastsim/scripts/coverage_and_stat_maps.py --ion 6Li --pdf toy
   python3 fastsim/scripts/coverage_and_stat_maps.py --ion 6Li --pdf grid
   python3 fastsim/scripts/coverage_and_stat_maps.py --ion 7Li --pdf toy
   python3 fastsim/scripts/coverage_and_stat_maps.py --ion 7Li --pdf grid
   ```
   Confirm 3 per-config PNGs × N_configs + 3 combined PNGs are written non-empty per run; no NaN colorbar ticks.

2. **Numerical cross-check vs. `phase_space_map.py`.** With matching `(ion, L, Pzz)`, output (2) must be numerically identical to the `err_cos2phi_<tag>.png` produced by `phase_space_map.py` (both invoke `fom.project_observables` on the same scenario). Under `--debug-crosscheck`, the script re-runs the reference path in-process and asserts `np.allclose(obs["err_a_cos2phi"], reference)`.

3. **`money_delta.py` untouched.** Run:
   ```
   git status fastsim/scripts/money_delta.py    # expect: clean
   python3 fastsim/scripts/money_delta.py --ion 6Li --pdf toy
   python3 fastsim/scripts/money_delta.py --ion 7Li --pdf grid
   ```
   Diff resulting PNGs against a pre-change baseline (byte-identical under the same matplotlib/`Agg` backend; otherwise stdout summary lines must match exactly).

4. **Angle sanity.** Per config, print `min/median/max` of `θ_e` and `η_e` over accepted bins. Expect `θ_e ≈ π` (η ≈ −3 … −1) for low-`y`, opening to smaller `θ_e` at high `y`/`Q²`. Fail loudly if any accepted bin has `θ_e ∉ (0, π)` or `|η_e| > 3.5`.

5. **Linear-`L` scaling.** Run with `--lumi 1.0` and `--lumi 10.0`; assert `proj.n_events.sum()` scales by ≈ 10.

6. **Shape invariants.** Assert `theta_e.shape == proj.x.shape`, `err_rb.shape == (nx, N_THETA_BINS)`, `theta_edges` strictly monotonic.

7. **Rebin event-count conservation.** For each x-row: `np.isclose(n_rb[ix].sum(), proj.n_events[ix][proj.accepted[ix]].sum())`. Definitive check that pooled-count rebin is lossless in the accepted region.

8. **Rebin error consistency.** Pick one x-row where a single θ bin captures all accepted `(x, Q²)` cells; assert `err_rb[ix, bin] ≈ sqrt(2/n_rb[ix, bin]) / Pzz` and that it matches the analytic inverse-variance combination of `obs["err_a_cos2phi"][ix, proj.accepted[ix]]` (agrees because `Pzz` is constant).

9. **Sanitizer parity.** Assert `_cfg_tag(cfg)` equals the string `phase_space_map.py:64-65` would produce for the same `cfg`, guarding against silent divergence between scripts.

---

## Todo

### Phase A — Scaffolding (new file)
- [ ] Create `fastsim/scripts/coverage_and_stat_maps.py` with the module docstring shown above and module-level constants `L_DEFAULT_FB=10.0`, `PZZ_DEFAULT=0.80`, `N_THETA_BINS=24`, `MIN_EVENTS_RB=1.0`, `ERR_VMIN=1e-4`, `ERR_VMAX=1.0`.
- [ ] Add imports: `numpy`, `matplotlib.pyplot`, `matplotlib.colors.LogNorm`, `polli_fastsim.{beams, fom, kinematics as kin}`, `polli_fastsim.structure.{NuclearF2, build_backends}`, `polli_fastsim.polarized.{ToyG1, toy_b1, toy_delta_gluon}`.
- [ ] Add `argparse` block with the exact CLI shape in §CLI above.

### Phase B — Helpers (in the new script)
- [ ] Implement `_cfg_tag(cfg)` — filename-safe suffix, mirroring `phase_space_map.py:64-65` verbatim.
- [ ] Implement `_project_at_lumi(cfg, backends, pzz, lumi_fb) -> (proj, obs, theta_e)`. `θ_e` shape must equal `proj.x.shape = (nx, nq2)`.
- [ ] Implement `_rebin_err_vs_theta(proj, theta_e, pzz, n_theta_bins, min_events)` returning `(x_edges, theta_edges, err_rb, n_rb)`. Combination rule: sum `n_events` into `(x, θ_e)` cells with `np.histogram(..., weights=n)`, then compute `δA = sqrt(2 / N_pooled) / Pzz`. Do **not** weight by `n/e²`. `n_rb` is returned only for count-conservation validation.
- [ ] Implement `_combined_panel_1x3(per_config_plot_fn, cached_items, outpath, shared_norm)` — driver that renders three subplots sharing one `LogNorm`.

### Phase C — Plot builders (in the new script)
- [ ] `plot_xq2_coverage(proj, cfg, path)` — `proj.n_events` masked by `proj.accepted`, `LogNorm(vmin=1)`, title includes `sqrt(s_eN)` and `L = 10 fb⁻¹/u`. Rectilinear `pcolormesh(x_edges, q2_edges, n_events.T)`, log-log axes.
- [ ] `plot_err_cos2phi_xq2(proj, obs, cfg, path)` — `obs["err_a_cos2phi"]`, `LogNorm(1e-4, 1)`, title includes `Pzz = 0.80`. Same rectilinear geometry as coverage plot.
- [ ] `plot_err_cos2phi_theta(proj, theta_e, cfg, pzz, path)` — calls `_rebin_err_vs_theta`, then rectilinear `pcolormesh(x_edges, y_edges_mrad, C.T)`; y-axis = `(π − θ_e)·1e3` [mrad], log-log axes, same colour scale as output (2).

### Phase D — `main()` in the new script
- [ ] Parse CLI, build backends via `structure.build_backends(pdf=args.pdf)`, create `outdir`.
- [ ] Loop `for cfg in beams.default_configs(args.ion):` — call `_project_at_lumi` once, then the three plot builders; cache `(cfg, proj, obs, theta_e)`.
- [ ] After loop, build the three 1×3 combined panels from cached arrays with a shared `LogNorm` per panel type.
- [ ] Print a one-line-per-config summary: `cfg.label(), N_DIS, median δA_cos2φ, θ_e range [mrad-from-backward]`.

### Phase E — Validation
- [ ] Confirm `git status fastsim/scripts/money_delta.py` is clean at all times.
- [ ] Smoke-test 4 combinations of `--ion {6Li,7Li}` × `--pdf {toy,grid}`; confirm expected PNG count and non-empty size.
- [ ] Under `--debug-crosscheck`, assert `np.allclose(obs["err_a_cos2phi"], phase_space_reference)` at matching `(ion, L, Pzz)`.
- [ ] Assert `θ_e ∈ (0, π)` and `|η_e| ≤ 3.5` on all accepted bins (recompute η locally).
- [ ] Assert linear-`L` scaling by running `--lumi 1.0` and `--lumi 10.0`; ratio ≈ 10.
- [ ] Assert shape invariants: `theta_e.shape == proj.x.shape`, `err_rb.shape == (nx, N_THETA_BINS)`, `theta_edges` strictly monotonic.
- [ ] Rebin event-count conservation per x-row.
- [ ] Rebin error consistency on a one-bin x-row.
- [ ] Sanitizer parity: `_cfg_tag(cfg)` == string produced by `phase_space_map.py:64-65`.

### Phase F — Documentation
- [ ] Ensure the module docstring in the new script lists the three outputs, notes `L = 10 fb⁻¹/u` and `Pzz = 0.80` defaults, and states that φ is analytically integrated out (so the angle version uses `θ_e`, not φ).
- [ ] Do **not** touch `money_delta.py` docstring or code.

---

## Handoff

Implementation is a fresh standalone script with a clean CLI and no cross-module edits — the right owner is `script-implementer`. Once implemented, the deliverable figures are produced by running the four smoke-test commands in Validation §1.
