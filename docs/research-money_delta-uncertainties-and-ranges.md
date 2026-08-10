# Research: statistical uncertainties and (x, Q²) ranges vs. beam energy / angle in the `money_delta` / fastsim workflow

Scope: determine what the fastsim projection pipeline that backs `fastsim/scripts/money_delta.py` currently exposes with respect to
1. per-bin statistical uncertainties,
2. any "angle" variable (electron `theta/eta`, azimuthal `phi`, hadron/spectator angles) that could be plotted alongside them,
3. the exact beam-energy configurations available, and
4. the most direct way today to visualize (x, Q²) acceptance/range by beam configuration.

All findings below are traced to concrete files and symbols in the working tree; nothing was edited.

---

## 1. The projection pipeline behind `money_delta.py`

`money_delta.py` is a thin driver on top of the same projection primitives used by every other `money_*.py` / `phase_space_map.py` script:

- `fastsim/scripts/money_delta.py:33-37` imports `beams`, `fom`, `asymmetries.a_cos2phi`, `structure.NuclearF2`, `polarized.toy_delta_gluon`, `inputs.get_backends`.
- The workhorse call is `fom.project_rates(cfg, sc, nuclear_f2=...)` (`money_delta.py:51`), which returns a `BinnedProjection`.

### 1.1 `BinnedProjection` — arrays actually available downstream

Defined in `fastsim/polli_fastsim/fom.py:47-55`:

```python
@dataclass
class BinnedProjection:
    x_edges: np.ndarray
    q2_edges: np.ndarray
    x: np.ndarray          # 2D centers (nx, nq2)
    q2: np.ndarray
    accepted: np.ndarray   # bool mask
    n_events: np.ndarray
    extras: dict = field(default_factory=dict)
```

Built in `project_rates(...)` (`fom.py:58-86`) on a log-log grid from `kin.log_grid` (`kinematics.py:73-79`, default `nx=40`, `nq2=30`, `x_range=(1e-4, 1.0)`, `q2_range=(1.0, 2e3)` GeV²).

`extras` is populated with:

- `"y"` — the inelasticity `y = Q²/(s·x)` (`kinematics.py:23`, computed at `fom.py:70`),
- `"eta"` — pseudorapidity of the *scattered electron* (`fom.py:71`, from `kin.scattered_electron`),
- `"s"` — `sqrt(s_eN)²`,
- `"nf2"` — the `NuclearF2` object used for the cross section.

Notably absent from `extras`: `theta` (the scattered-electron polar angle) and `e_prime` (its lab energy). Both are computed inside `project_rates` at `fom.py:71` (`e_p, _theta, eta = kin.scattered_electron(...)`) but only `eta` is stored — `theta` is discarded via `_theta` and `e_p` is only used for the `e_prime_min` cut.

### 1.2 Acceptance cuts baked into `Scenario`

`Scenario` (`fom.py:31-44`) provides the acceptance envelope used to build `accepted`:

- `q2_min=1.0`, `y_min=0.01`, `y_max=0.95`, `w2_min=10.0` — inclusive DIS cuts applied via `kin.kinematic_mask` (`kinematics.py:52-70`).
- `eta_min=-3.5`, `eta_max=3.5` and `e_prime_min=0.5 GeV` — crude central-detector cut on the scattered electron applied at `fom.py:73-74`.

The `accepted` boolean mask on the (x, Q²) grid is the AND of all of the above.

### 1.3 Events per bin (the raw statistical driver)

`n_events` is `xsec * dx * dq2 * lumi_pb` masked by `accepted` (`fom.py:80-83`), with `xsec = dsigma_dx_dq2(X, Q2, s, F2/A)` in pb/GeV²/nucleon and `lumi = scenario.lumi_fb_per_nucleon * 1e3` pb⁻¹/nucleon. This is the sole per-bin quantity from which every statistical uncertainty downstream is derived.

---

## 2. Are per-bin statistical uncertainties available?

**Yes, fully — and they are already used by `money_delta.py` (implicitly) and by `scripts/phase_space_map.py` (explicitly, plotted).**

### 2.1 Formulas — `fastsim/polli_fastsim/asymmetries.py:56-75`

```python
def err_a_parallel(n, pe, pz):     # delta(A_par)
    return 1.0 / (pe * pz * np.sqrt(n))

def err_azz(n, pzz):               # delta(A_zz)
    return np.sqrt(2.0 / n) / pzz

def err_cos2phi_amplitude(n, pzz): # delta(A_{cos 2phi})
    return np.sqrt(2.0 / n) / pzz
```

All three are element-wise NumPy functions on `proj.n_events`, so they naturally return arrays with the same `(nx, nq2)` shape as the (x, Q²) grid.

### 2.2 The convenience wrapper that attaches them

`fom.project_observables(config, scenario, proj, g1_model, b1_func, delta_func)` (`fom.py:89-122`) — this is the function that materializes per-bin errors for all three observables in one call. It stores into `proj.extras`:

- `"a_par"`, `"err_a_par"`, `"err_g1_over_f1"` (the vector-spin channel),
- `"azz"`, `"err_azz"` (tensor b1 channel),
- `"a_cos2phi"`, `"err_a_cos2phi"` (the Δ / gluonometry channel — the one `money_delta.py` targets),
- `"sig_a_par"`, `"sig_azz"`, `"sig_a_cos2phi"` (|A|/δA, zeroed outside acceptance).

There is also a Q²-combination helper `fom.combine_over_q2(err_map, accepted, min_events, n_events)` (`fom.py:125-136`) that folds per-(x,Q²) errors into a per-x error via inverse-variance summation. This is how `money_polemc.py` collapses to δ(g1/F1) vs x; it works for *any* per-bin error map.

### 2.3 What `money_delta.py` actually does with them

Interestingly, `money_delta.py` **does not call `project_observables`**. It reimplements a compressed significance-squared per fb⁻¹ inline (`money_delta.py:45-59`):

```python
use = proj.accepted & (proj.n_events >= min_events)
return np.where(use, amp**2 * pzz**2 * proj.n_events / 2.0, 0.0).sum()
```

This is arithmetically equivalent to `sum_bins (A_bin/δA_bin)²` using `err_cos2phi_amplitude` — i.e., the per-bin statistical uncertainty on A_cos2φ is baked into the reach formula via `sqrt(2/N)/Pzz`, then summed over accepted (x, Q²) bins. So the machinery is present; it is simply not exposed as an array in `money_delta.py`'s current outputs.

**Bottom line**: showing per-bin δ maps for A_cos2φ (or A∥, A_zz) at each beam energy is a straight consumer of existing outputs — call `fom.project_observables(cfg, sc, proj, ...)` after `project_rates`, then read `proj.extras["err_a_cos2phi"]`. `scripts/phase_space_map.py:70-81` already does exactly this (with `err_g1_over_f1`, `err_azz`, `err_a_cos2phi`) and produces one PNG per beam configuration.

---

## 3. Angle representations in the codebase

This is where support is genuinely limited. Three angle families exist, at three very different levels of readiness:

### 3.1 Azimuthal angle `phi` of the scattered electron — NOT a projection dimension

`money_delta.py` is explicit about this in its module docstring (`money_delta.py:13-22`) and again inline (`money_delta.py:264-272`):

> "The azimuthal angle phi of the scattered electron around the beam axis is NOT a binned dimension in the projection model. The rate arrays are integrated over phi using the standard inclusive-DIS cross section (phi-averaged). The cos(2phi) amplitude A_cos2phi returned by asymmetries.a_cos2phi is the *coefficient* of the phi modulation that would be extracted from an azimuthal fit — it is a function of (x, Q2, y) only."

The rate integrand `dsigma_dx_dq2` (used at `fom.py:78`) is the φ-integrated inclusive DIS cross section. A true φ-binned significance map would require event-level generation — noted as out-of-scope in the same docstring.

The spin-1 master formula that *defines* the φ dependence lives in `asymmetries.py:1-18` (docstring) and is used only to build the amplitude at (x, Q², y). Target-spin polar angle `theta_m` is a parameter of `azz(...)` (`asymmetries.py:41`, defaulted to 0 = longitudinal) and appears via `sin²(theta_m)` in the master formula, but it is a *setup* parameter (target spin orientation), not a per-bin observable.

### 3.2 Scattered-electron polar angle `theta` / pseudorapidity `eta` — COMPUTED but only `eta` is stored

`kin.scattered_electron(x, y, s, electron_energy)` (`kinematics.py:32-49`) returns `(e_prime, theta, eta)` in the lab frame from the standard massless DIS kinematics:

```python
Q² = 2 Ee E'(1+cosθ),   y = 1 − (E'/2Ee)(1−cosθ)
```

`fom.project_rates` at `fom.py:71` throws away `theta` (`_theta`) and stores only `eta` in `extras["eta"]`. Both are trivially recoverable — either by binding the return value or re-calling `kin.scattered_electron(proj.x, proj.extras["y"], proj.extras["s"], cfg.electron_energy)` after the fact. No new modeling is needed.

This is the natural handle for "different (scattered-electron) angles". A per-bin θ or η map is one line of Python on top of what already exists.

### 3.3 Hadron / spectator angles — present, but in a completely separate pipeline

`fastsim/polli_fastsim/spectator.py` samples cluster-spectator kinematics (α+d for ⁶Li, α+t for ⁷Li) via a two-parameter momentum-density MC, boosts to the lab and returns a dict containing:

- `pT`, `theta` (rad), `p_lab`, `R` (rigidity ratio) — see `spectator.py:114-136`.

These are Monte-Carlo *samples*, not binned rates on the (x, Q²) grid, and they describe the *forward spectator cluster*, not the struck-nucleon hadron system. `farforward.py:41-70` consumes them to route into RP / OMD / B0 / ZDC based on `theta` and `R`. `scripts/tagging_acceptance.py` is the only driver that exercises this angle information; `money_delta.py` does not touch it.

There is no representation of the *struck hadron* (jet, current-fragmentation hadron) angle anywhere in the codebase — no SIDIS or semi-inclusive angular variables. That would be new physics + new event generation.

---

## 4. Beam-energy configurations available

Defined in `fastsim/polli_fastsim/beams.py:54-83`:

- `ELECTRON_ENERGIES = (5.0, 10.0, 18.0)` GeV (line 54) — the three canonical EIC electron-ring settings.
- `default_configs(ion_name)` (lines 75-83) returns three `BeamConfig` objects mirroring ep 5×41, 10×100, 18×275, with the ion momentum scaled by rigidity `Z/A`:

    | E_e (GeV) | proton-analog p (GeV) | ⁶Li p/u (Z/A=1/2) | ⁷Li p/u (Z/A=3/7) |
    |---|---|---|---|
    | 5.0 | 41 | 20.5 | 17.6 |
    | 10.0 | 100 | 50.0 | 42.9 |
    | 18.0 | 275·(Z/A) | 137.5 | 117.9 |

The 18-GeV row uses `top = ion.momentum_per_nucleon_max = 275·Z/A` (via `Ion.momentum_per_nucleon_max`, `beams.py:40-43`), which matches the EPIOS white-paper values quoted in the module docstring (⁶Li ~138 GeV/u, ⁷Li ~117 GeV/u).

`BeamConfig` exposes:

- `.sqrt_s_per_nucleon` = `sqrt(4·Ee·Ep_per_nucleon)` (`beams.py:63-66`),
- `.label()` → e.g. `"e(10) x 7Li(42.9/u)"` (`beams.py:68-72`).

`money_delta.py:225-232` iterates over `beams.default_configs(args.ion)` for `--ion` ∈ `{6Li, 7Li}` (argparse choices at `money_delta.py:209`). No hooks exist for arbitrary user-supplied energies; adding one is a matter of constructing a `BeamConfig` manually.

Ion species with full parameters (name, A, Z, spin, effective nucleon polarizations): `p, d, ³He, ⁶Li, ⁷Li` in `beams.py:46-52`.

---

## 5. Most direct way to visualize (x, Q²) acceptance / range by beam configuration

**Already implemented in `fastsim/scripts/phase_space_map.py`** — it is the canonical answer to this question in the repo.

Concretely (`phase_space_map.py:29-43, 60-81`):

```python
def plot_map(proj, values, title, label, path, vmin, vmax, log=True):
    vals = np.ma.masked_where(~proj.accepted | (values <= 0), values)
    ax.pcolormesh(proj.x_edges, proj.q2_edges, vals.T, norm=LogNorm(...), ...)
    ax.set_xscale("log"); ax.set_yscale("log")
```

Per beam configuration it produces:
- `rate_<label>.png` — `proj.n_events` on a log-log (x, Q²) plane, masked outside `accepted` (this *is* the acceptance-range visualization, weighted by rate).
- `err_g1f1_<label>.png`, `err_azz_<label>.png`, `err_cos2phi_<label>.png` — per-bin statistical precision maps for the three observables, using `obs["err_g1_over_f1"]`, `obs["err_azz"]`, `obs["err_a_cos2phi"]` from `fom.project_observables`.

The (x, Q²) *range* per configuration is fully determined by `proj.accepted` (which encodes Q² ≥ 1, y ∈ [0.01, 0.95], W² ≥ 10, and −3.5 ≤ η_e ≤ 3.5, E'_e ≥ 0.5 GeV). Overlaying `accepted` boundaries for the three energies onto a single (x, Q²) axis is a natural extension — the arrays required are already there via `beams.default_configs(ion_name)` → `fom.project_rates(cfg, sc)` → `proj.accepted`.

`money_delta.py` itself does not draw an (x, Q²) map; it only produces its Δ/F1-vs-luminosity curve (Figure 1) and the beam-energy × Δ/F1 heatmap (Figure 2). Adding an (x, Q²) subplot per configuration would consume the same objects.

---

## 6. Summary — what is supported, what is not, what needs new work

### Currently supported (no new code, only consumption of existing arrays)

- **Per-bin statistical uncertainties** for A_cos2φ (`err_cos2phi_amplitude`), A_zz (`err_azz`), A∥ (`err_a_parallel`), and δ(g1/F1) via `depolarization_d` — all available element-wise on the (x, Q²) grid through `fom.project_observables` and stored in `proj.extras` (`fom.py:97-121`). Q²-combined per-x errors via `fom.combine_over_q2`.
- **(x, Q²) acceptance/range at each of the three beam energies**, per ion, from `proj.accepted` (built with Yellow-Report-like DIS cuts plus a crude central-detector electron η/E' cut). Directly plotted by `scripts/phase_space_map.py`; `money_delta.py` computes the same `proj` but does not plot the coverage.
- **Beam-energy scan** = the fixed 3-point set from `beams.default_configs(ion_name)`: e(5)×41-analog, e(10)×100-analog, e(18)×top-rigidity, for `6Li` and `7Li` (and `p`, `d`, `3He` if the caller passes those names).
- **Scattered-electron pseudorapidity** on the (x, Q²) grid (`proj.extras["eta"]`); polar angle `theta` and lab energy `E'` are computed inside `project_rates` but not persisted — one-line change (rebind `_theta`, or re-call `kin.scattered_electron`) to expose them without any new physics.

### Not currently supported (needs modest new code, no new physics or events)

- **A δA_cos2φ map plotted per beam configuration.** Nothing computes/plots it in `money_delta.py`; add a `project_observables` call and a `pcolormesh` block modeled on `phase_space_map.py:78-81`.
- **Per-bin θ_e / E'_e maps.** Values are computed once in `project_rates` but only η is stored. Trivial to expose.
- **Multi-energy (x, Q²) coverage overlay** (one panel comparing the three `accepted` regions). Arrays exist; only a new plotting script is needed.
- **Q²-combined δ vs x** for A_cos2φ. `fom.combine_over_q2` already accepts arbitrary error maps; `money_polemc.py` does the analogous thing for δ(g1/F1). Would need to be added as an output of `money_delta.py`.

### Would require new modeling and/or event generation

- **A per-φ-bin statistical figure of merit.** The rate array `n_events` is integrated over φ by construction (`dsigma_dx_dq2` is the φ-averaged inclusive cross section). Splitting statistics into φ bins requires either an analytical extension of the master formula into a differential dσ/dφ and re-binning, or true event-level MC. This is called out explicitly in `money_delta.py:264-272`.
- **Angle-dependence in the sense of hadron-level (SIDIS / current-jet) angles.** There is no representation of struck-quark fragmentation, hadron azimuth around q, or hadron polar angle anywhere in the projection pipeline. Would require importing/interfacing an event generator (BeAGLE is available in `BeAGLE_jz/`, but not wired into fastsim).
- **Target-spin-angle scans (θ_m ≠ 0 for A_zz, θ_m ≠ 90° for A_cos2φ).** `azz(..., theta_m=0)` supports it as an argument (`asymmetries.py:41-47`), but there is no driver that sweeps it, and the tensor-polarization setup in `Scenario` is a single scalar `pol_ion_tensor`. Would require re-deriving errors for arbitrary θ_m (currently `err_azz` and `err_cos2phi_amplitude` bake in the two canonical setups).
- **Beam energies outside the three default points.** No `argparse` hook; requires constructing custom `BeamConfig` objects. Trivial in principle but not exposed at the script level.

### Risks / caveats

- `money_delta.py` uses a toy `polarized.toy_delta_gluon` for the Δ shape; per-bin significance maps built from it are only qualitatively meaningful until a real Δ model is plugged in.
- The central-detector acceptance is a hard η ∈ [−3.5, 3.5] and E'_e ≥ 0.5 GeV cut — realistic detector efficiencies (tracking, EMcal, PID) are not modeled. Any "acceptance" figure derived from `proj.accepted` inherits this crudeness.
- Ion polarizations in `Scenario` (`P_z = 0.70`, `P_zz = 0.60`) are placeholders. Per-bin *statistical* δ scales as `1/(P·√N)` so the shape of the maps is P-independent; only absolute normalization moves.
- `money_delta.py`'s inline reach formula assumes Δ is strictly linear in `scale` (`money_delta.py:47-59`); if that assumption is broken by a future Δ model, the analytic `L_5σ = 25 / (sig² · (s/s₀)²)` extrapolation used to fill the FOM heatmap must be replaced with a per-scale recomputation.

---

## 7. Key file/function index

| File | Symbol | Purpose |
|---|---|---|
| `fastsim/scripts/money_delta.py` | `sig2_per_fb_at`, `build_fom_heatmap`, `plot_fom_heatmap`, `main` | Δ/F1 reach; does *not* plot per-bin δ or (x,Q²) coverage |
| `fastsim/scripts/phase_space_map.py` | `plot_map`, `main` | Reference for (x,Q²) rate + δ-per-bin visualization |
| `fastsim/polli_fastsim/fom.py` | `Scenario`, `BinnedProjection`, `project_rates`, `project_observables`, `combine_over_q2` | Pipeline core; `extras` dict holds y, eta, s, nf2 (+ observables when `project_observables` is called) |
| `fastsim/polli_fastsim/kinematics.py` | `y_from_xq2`, `w2`, `scattered_electron`, `kinematic_mask`, `log_grid` | (x,Q²) grid, scattered-electron (E', θ, η), inclusive DIS cuts |
| `fastsim/polli_fastsim/beams.py` | `Ion`, `BeamConfig`, `IONS`, `ELECTRON_ENERGIES`, `default_configs` | Three-point energy scan per ion; rigidity-scaled ion momenta |
| `fastsim/polli_fastsim/asymmetries.py` | `a_parallel`, `azz`, `a_cos2phi`, `err_a_parallel`, `err_azz`, `err_cos2phi_amplitude`, `depolarization_d` | Per-bin asymmetries and their √N-based statistical errors |
| `fastsim/polli_fastsim/spectator.py` | cluster sampler → `{pT, theta, p_lab, R}` | Forward-spectator angles (separate pipeline; not used by `money_delta`) |
| `fastsim/polli_fastsim/farforward.py` | `route_charged`, `route_neutral` | Consumes spectator θ, R for RP/OMD/B0/ZDC routing |
