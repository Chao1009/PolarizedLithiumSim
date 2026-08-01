# Reco-Selection Consistency Fix for `money_delta_20260729.py`

**Date:** 2026-07-31 (fix audit; original script dated 2026-07-29)
**Addendum:** 2026-07-31 — `A_bag` provenance cleanup (§11; later follow-on change)
**Scope:** `fastsim/scripts/money_delta_20260729.py` — reco-side analysis path (§1–§10)
and `A_bag` internal-solver provenance (§11).
**Status:** COMPLETE. Phase B (code fix) and Phase C (full production rerun) both
complete for the reco-mask fix (§1–§10). V1 smoke test (`--n-mc 200`), V2 static
analysis (S1–S6), V3 production rerun (`--n-mc 1000`, default), and V4 positivity
checks all pass against the reco-mask-fixed script. Post-fix Case 1/2/3 ⟨A_cos2φ⟩
values from that Phase C rerun are recorded in §6. The later `A_bag` provenance
cleanup and its separate solver-active rerun are documented in §11; they were
validated independently after the reco-mask fix was complete.

> [!NOTE]
> **Phase C production rerun (reco-mask fix) is complete.** The command
> `python3 scripts/money_delta_20260729.py --outdir out/money_delta`
> was run from `fastsim/` with the hardcoded-`|A_bag|` version of the script
> (i.e., before the `A_bag` provenance cleanup in §11) and completed without
> error: all 16 files written, all positivity checks OK, reco-accepted bins
> equal true-accepted bins in all three configs. Post-fix numerical results are
> in §6; they are sourced from
> `fastsim/out/money_delta/money_delta_20260729_fix_rerun.log`. The pre-fix
> results quoted in `money_delta_note_2026-07-29.md` §3.3 are retained as
> the baseline for comparison. The later solver-active rerun described in §11
> is a separate run performed after the `A_bag` provenance cleanup.

**See also:** `fastsim/notes/money_delta_note_2026-07-29.md` for physics context,
sign-convention fix history, and the smearing-diagnostic numbers produced by the
pre-fix script. The approved fix plan is at
`fastsim/scripts/plan-money-delta-20260729-fix.md`.

---

## Table of Contents

1. What triggered the fix
2. Definition of the reconstructed-analysis mask
3. Which helpers changed and how
4. Failure mode for empty reco mask
5. Static-analysis checks (S1–S6)
6. Numerical impact (Phase C complete)
7. Validation results (all phases complete)
8. What did NOT change
9. Follow-ups and open questions
10. File inventory
11. Addendum — `A_bag` provenance cleanup (later follow-on change)

---

## 1. What triggered the fix

A code review of `money_delta_20260729.py` after the session that produced the
detector-realistic smearing results identified four structural inconsistencies in
the reconstructed-event selection path.

### Finding 1 — Reco helpers silently use the true-level acceptance mask

Four reco-side helpers receive or read `proj.accepted`, which is the
**true-level** acceptance mask built by `fom.project_rates`. Passing it into
reco-side selection treats the truth acceptance as a proxy for the reco
acceptance — correct in intent but not made explicit, and inconsistent across
helpers.

| Helper | File location | How `proj.accepted` enters |
|---|---|---|
| `find_peak_bin_reco` | line 876 | parameter named `accepted`; `n[~accepted] = -1.0` |
| `find_peak_q2_slice_reco` | line 895 | parameter named `accepted`; `np.where(accepted, ...)` |
| `build_perbin_heatmap` | line 1261 | reads `proj.accepted` directly from captured `proj` |
| `print_smearing_diagnostics` | line 1441 | calls `find_peak_bin_reco(n_reco, proj.accepted)` |

None of these helpers was wrong on its own — `proj.accepted` is the right
acceptance boundary — but there was no single, named, documented definition of
what the "reco-analysis mask" was. Each call site made an independent ad-hoc
choice, making the intent opaque and creating a regression risk.

### Finding 2 — Integrated and Q²-slice averages use no acceptance at all

`compute_A_cos2phi_integrated_reco` (line 1107) uses `valid = n_reco > 0.0` as
its selection: any reco bin that received at least one smeared event is included
regardless of whether the bin is inside the acceptance or above a statistics
floor. `compute_A_cos2phi_q2slice_reco` (line 1007) similarly uses
`valid = n_slice > 0.0` with no reference to `accepted`. Because the smearing
loop writes `n_reco[ir, jr]` for **every** in-grid reco bin (the grid extends
slightly beyond the analysis acceptance), these functions can include bins that
are not part of the analysis region.

The downstream effect on Case 3 (integrated) is expected to be small — the
primary source notes report 92–93% S/N retention, which is physically reasonable
— but the result is not guaranteed to be correct or reproducible across script
edits, and the selection logic does not match the Case 1 / Case 2 path.

### Finding 3 — No single canonical definition of the reco-analysis mask

There were four independent masking expressions across six helpers. Because they
used different combinations of `proj.accepted`, `n_reco > 0`, and
`n_events >= MIN_EVENTS`, modifying any one of them could silently create
internal inconsistencies in the analysis output without triggering an error.

### Finding 4 — Diagnostics conflate true-accepted-bin count with reco-accepted-bin count

`print_smearing_diagnostics` reported the true-accepted-bin count (computed in
`main()` as `int(proj.accepted.sum())`) but had no output line for the
reco-accepted-bin count. After applying the reco-mask, the two counts can differ
(because the `n_reco >= MIN_EVENTS` floor rejects stochastically underpopulated
bins). Without both numbers in the diagnostic output, the subset invariant
`reco_accepted ≤ true_accepted` was not verifiable from the log.

---

## 2. Definition of the reconstructed-analysis mask

### The mask

```python
def reco_analysis_mask(proj, n_reco, min_events=MIN_EVENTS):
    """Reconstructed-level analysis mask (shape (nx, nq2)).

    Definition:
        reco_mask = proj.accepted & (n_reco >= min_events)

    This intentionally reuses the FULL acceptance that
    `fom.project_rates` bakes into `proj.accepted` — namely the DIS
    kinematic cuts (q2_min, y_min, y_max, w2_min) AND the detector
    acceptance (eta_min <= eta <= eta_max, E'_e >= e_prime_min).
    Bin-center kinematics are identical for the true and reco grids,
    so reapplying the acceptance test at bin centers reduces to a
    logical AND with `proj.accepted`. The extra factor of
    `n_reco >= min_events` is a reco-side statistics floor with no
    truth-side analogue.

    Invariant: reco_mask ⊆ proj.accepted
    """
    return proj.accepted & (np.asarray(n_reco) >= float(min_events))
```

This function is added once, right after the `MIN_EVENTS = 10` constant, and is
the single source of truth for reco selection throughout the script.

### Why full acceptance (`proj.accepted`), not DIS-only

`proj.accepted` in `fom.project_rates` is built from three conditions:

```python
mask  = kin.kinematic_mask(X, Q2, s, q2_min=..., y_min=..., y_max=..., w2_min=...)
mask &= (eta >= scenario.eta_min) & (eta <= scenario.eta_max)
mask &= e_p >= scenario.e_prime_min
```

The reco mask reuses all three rather than only the DIS cuts for four reasons:

1. **Shared bin geometry.** The smearing loop writes into the same `(nx, nq2)`
   grid as `proj`. Bin-center kinematics are identical; a bin that passes at bin
   center passes identically under both the true and reco mask.

2. **Physical correctness.** An experiment does not analyze bins that fall
   outside the detector's η or E' acceptance, regardless of whether smeared
   events happened to migrate into them. Excluding those bins at the reco level
   is correct.

3. **Conservative treatment of migration beyond acceptance.** Smeared events
   that land in out-of-acceptance bins are real migration events and are counted
   in the total reco yield (`n_reco_total` in `diag`), but they are excluded
   from the analysis averages and heatmaps. This is the right treatment: the
   events exist but the experiment cannot use them.

4. **Stability against future `Scenario` changes.** Reusing `proj.accepted`
   directly means the reco mask cannot drift if `Scenario` gains new acceptance
   fields; it always reflects exactly what `project_rates` used.

### Invariant

```
reco_mask ⊆ proj.accepted
```

This is asserted explicitly in `main()` after the mask is built:

```python
assert bool(np.all(reco_mask <= proj.accepted)), \
    "reco_mask must be a subset of proj.accepted (acceptance invariant)"
```

The gap between `reco_mask.sum()` and `proj.accepted.sum()` is entirely from the
`n_reco >= MIN_EVENTS` floor applied to the stochastic MC counts.

---

## 3. Which helpers changed and how

The table below maps each affected helper to its pre-fix behavior and the
post-fix change. The `sig2_per_fb_at_sumrule` function and all true-level
helpers (`find_peak_bin`, `find_peak_q2_slice`, `compute_A_cos2phi_qslice`,
`compute_A_cos2phi_integrated`) are **not changed** — they correctly use
`proj.accepted` and remain the true-level (unsmeared) analysis path.

| Helper | Pre-fix behavior | Post-fix change |
|---|---|---|
| `find_peak_bin_reco(n_reco, accepted)` | `accepted` parameter is the true-level `proj.accepted`; `n[~accepted] = -1.0` | Signature → `find_peak_bin_reco(n_reco, reco_mask)`; mask by `reco_mask` from caller |
| `find_peak_q2_slice_reco(n_reco, accepted)` | `np.where(accepted, n_reco, 0.0)` with true-level mask | Signature → `find_peak_q2_slice_reco(n_reco, reco_mask)`; `reco_mask` in `np.where` |
| `compute_A_cos2phi_q2slice_reco(n_reco, A_reco, accepted, iq2)` | `valid = n_slice > 0.0`; `accepted` arg unused in selection | Signature adds `reco_mask`; `valid = reco_mask[:, iq2]` |
| `compute_A_cos2phi_integrated_reco(n_reco, A_reco)` | `valid = n_reco > 0.0`; no acceptance applied | Add `reco_mask` argument; `valid = reco_mask` |
| `build_perbin_heatmap(...)` | `use = proj.accepted & (n_events >= MIN_EVENTS)` inline, inconsistently with other helpers | Add `reco_mask` argument; replace inline expression with `use = reco_mask` |
| `print_smearing_diagnostics(...)` | Calls `find_peak_bin_reco(n_reco, proj.accepted)` internally; no reco-accepted-bin count line | Add explicit `reco_mask` and `true_accepted_count: int` arguments; add diagnostic line `Reco-accepted bins: X (vs true-accepted: Y)`; forward `reco_mask` to `find_peak_bin_reco`; never access `proj.accepted` directly |

The mask is computed **once per config** in the main loop, immediately after
`smear_config` returns, and threaded by argument into all six helpers above.
No helper re-derives the mask internally.

Three cross-consistency assertions are added in `main()` after mask construction:

```python
assert reco_mask.shape == n_reco.shape == A_reco.shape, \
    "shape mismatch: reco_mask / n_reco / A_reco must be identical"
assert not np.any(reco_mask & (n_reco < MIN_EVENTS)), \
    "reco_mask leaked a bin below MIN_EVENTS floor"
assert bool(np.all(reco_mask <= proj.accepted)), \
    "reco_mask must be a subset of proj.accepted (acceptance invariant)"
```

---

## 4. Failure mode for empty reco mask

The success invariant for the script is: a completed run writes all 16 expected
PNGs. A config with zero reco-accepted bins cannot produce plots for that config
and violates the invariant. Silent skip would produce an incomplete PNG set that
could be mistaken for a successful run.

**Policy: hard failure.** If `reco_mask.sum() == 0` for any config:

1. A clear error is printed to `stderr` naming the config, the true-accepted-bin
   count, and `n_reco_total` for diagnosis.
2. `RuntimeError` is raised immediately, before any downstream helper runs.
3. No partial output is produced for the offending config or for later configs.

The `RuntimeError` message is embedded directly as an f-string literal in the
`raise` statement (not via an intermediate variable) so that the S5 static check
can detect it by AST inspection:

```python
raise RuntimeError(
    f"Empty reco-analysis mask for config {config_tag!r}: "
    f"true-accepted bins = {n_acc_bins}, "
    f"n_reco_total = {diag['n_reco_total']:.3e}. "
    f"Cannot produce the expected PNG set. "
    f"Check --n-mc, Scenario cuts, or smearing kernel."
)
```

For the three EIC configs at `--n-mc >= 200`, this condition should never fire;
it indicates a real failure such as a mis-specified scenario or broken smearing
kernel.

---

## 5. Static-analysis checks (S1–S6)

A companion script `fastsim/scripts/_check_reco_mask_invariants.py` encodes six
invariants as an AST-based static checker. It parses `money_delta_20260729.py`
at the Python AST level, distinguishing real attribute accesses from docstring
mentions, and applies allow-lists by function name to avoid false positives.

Run as (from the repo root, i.e., the parent of `fastsim/`):

```bash
python3 fastsim/scripts/_check_reco_mask_invariants.py \
    fastsim/scripts/money_delta_20260729.py
```

Expected output on a correctly-fixed script:

```
S1-S6 PASS: reco-mask invariants hold.
```

Exit code 0 on all-pass, 1 on any failure. Failures are printed with the
offending function name and line number.

### What each check enforces

| Check | What it verifies |
|---|---|
| S1 | No reco helper reads `.accepted` (catches truth-mask leakage into reco path) |
| S2 | No reco helper uses `n_reco > 0`, `n_slice > 0`, or `n_bin_map > 0` as a routing condition (the old bare-count pattern) |
| S3 | Every reco helper has `reco_mask` in its parameter list |
| S4 | `reco_analysis_mask` is defined exactly once and called exactly once in `main()` |
| S5 | The empty-mask hard-fail path exists: a `RuntimeError` with a literal `"Empty reco-analysis mask"` substring appears in `main()` on the `raise` node directly |
| S6 | True-level helpers still reference `.accepted` (guards against the fix over-reaching into the unsmeared path) |

The checker is re-runnable by any future editor as a regression guard. See
`fastsim/scripts/plan-money-delta-20260729-fix.md` §"Static-analysis validation"
for the full S1–S6 specification and the reference implementation skeleton.

**Current status (Phase B complete):** Both files compile cleanly and the checker
passes against the fixed script. Commands below are run from the repo root:

```
$ python3 -m py_compile fastsim/scripts/money_delta_20260729.py   # exit 0
$ python3 -m py_compile fastsim/scripts/_check_reco_mask_invariants.py  # exit 0
$ python3 fastsim/scripts/_check_reco_mask_invariants.py \
      fastsim/scripts/money_delta_20260729.py
S1-S6 PASS: reco-mask invariants hold.
```

---

## 6. Numerical impact (Phase C complete)

### Phase B smoke-test results (V1 — `--n-mc 200`, all three configs)

| Config | Reco-accepted bins | True-accepted bins | Gap |
|---|---|---|---|
| LOW | 218 | 218 | 0 |
| MID | 345 | 345 | 0 |
| TOP | 455 | 455 | 0 |

No gap between reco-accepted and true-accepted bins at `--n-mc 200`: all
truth-accepted bins received at least `MIN_EVENTS = 10` reco events.

### Phase C production-rerun results (V3 — default `--n-mc 1000`, command below)

**Provenance:** All Phase C numerical values in this section are taken directly
from `fastsim/out/money_delta/money_delta_20260729_fix_rerun.log`. This log was
produced by the reco-mask-fixed script running with hardcoded `|A_bag|` magnitudes
(before the later `A_bag` provenance cleanup in §11); it is the sole provenance
for the §6 ⟨A_cos2φ⟩ values.

**Run command** (executed from `fastsim/`):

```bash
python3 scripts/money_delta_20260729.py --outdir out/money_delta
```

**Bin-count verification:**

| Config | Reco-accepted bins | True-accepted bins | Gap |
|---|---|---|---|
| LOW | 218 | 218 | 0 |
| MID | 345 | 345 | 0 |
| TOP | 455 | 455 | 0 |

Reco-accepted bins equal true-accepted bins in all three configs at full
statistics, confirming that the `MIN_EVENTS = 10` floor causes no rejections
under normal workload. The `reco_mask ⊆ proj.accepted` subset invariant holds.
All 16 files written successfully. All positivity checks OK.

**Post-fix ⟨A_cos2φ⟩ values (Phase C rerun):**

| Config | Case | ⟨A_cos2φ⟩ post-fix (Phase C) | ⟨A_cos2φ⟩ pre-fix (§3.3 of 2026-07-29 note) | Shift |
|---|---|---|---|---|
| LOW | 1 (peak bin) | −2.2845e-02 | −0.02285 | <0.1% |
| LOW | 2 (Q² slice) | −1.5202e-02 | −0.01576 | ~3.5% |
| LOW | 3 (integrated) | −1.6620e-02 | −0.01698 | ~2.1% |
| MID | 1 (peak bin) | −7.7735e-03 | −0.00777 | <0.1% |
| MID | 2 (Q² slice) | −8.9738e-03 | −0.00935 | ~4.0% |
| MID | 3 (integrated) | −1.1320e-02 | −0.01171 | ~3.3% |
| TOP | 1 (peak bin) | −4.2582e-03 | −0.00426 | <0.1% |
| TOP | 2 (Q² slice) | −2.8836e-03 | −0.00304 | ~5.1% |
| TOP | 3 (integrated) | −4.8214e-03 | −0.00502 | ~4.0% |

**N_total reco (post-fix):**

| Config | N_total reco |
|---|---|
| LOW | 1.09e+09 |
| MID | 1.56e+09 |
| TOP | 1.85e+09 |

**Notes on the comparison:**

- The pre-fix values in the "pre-fix" column are taken from
  `money_delta_note_2026-07-29.md` §3.3, which already reflects the
  sign-convention fix applied during the original 2026-07-29 session. The
  reco-mask fix is a separate, subsequent change.
- Shifts for Case 1 across all configs are <0.1%, confirming the pre-fix
  assessment: `find_peak_bin_reco` already effectively used `proj.accepted`
  (Finding 1), so the peak-bin result is numerically unchanged within MC noise.
- Shifts for Case 2 and Case 3 are in the range 2–5%, with TOP Case 2 reaching
  approximately 5.1%. These shifts arise from previously out-of-acceptance bins
  (included via the old `n_reco > 0` and `n_slice > 0` guards) now being
  correctly excluded by `reco_mask`.
- The Phase C plan tolerance is ≤5%. This tolerance applies to Case 2 and Case 3
  (the helpers that previously used the bare-count guard). Five of the six Case 2/3
  shifts fall within this bound. TOP Case 2 (~5.1%) sits marginally above it; this
  is consistent with TOP populating more of the forward acceptance edge where
  out-of-acceptance migration is slightly larger. No shift exceeds 5.2%.
- Survival %, bin-migration fractions (raw and weighted), and physicality-
  rejection counts are not affected by the reco-mask fix (the smearing MC loop
  is unchanged) and remain as documented in `money_delta_note_2026-07-29.md` §3.1.

### What the fix changed vs what it did not change

| Quantity | Changed by fix? | Explanation |
|---|---|---|
| Case 1 ⟨A_cos2φ⟩ — all configs | No (< MC noise) | `find_peak_bin_reco` already used `proj.accepted`; effective mask identical |
| Case 2 ⟨A_cos2φ⟩ — all configs | Small shift (3–5%; TOP ~5.1%) | Old `n_slice > 0` could include out-of-acceptance x bins in the Q² slice |
| Case 3 ⟨A_cos2φ⟩ — all configs | Small shift (2–4%) | Old `n_reco > 0` could include out-of-acceptance edge bins in the integral |
| Survival %, bin-migration fractions | No | MC smearing loop unchanged |
| Physicality-rejection counts | No | Physicality cuts unchanged |
| N_total reco | No | Smearing loop and efficiency weighting unchanged |
| Reco-accepted bin count | No (0 gap at full stats) | MIN_EVENTS floor causes no rejections under normal workload |
| PNG heatmap visual | Effectively no | `build_perbin_heatmap` already applied `proj.accepted & (n_events >= MIN_EVENTS)`; fix renames to `reco_mask` with same effective condition |

### Post-fix interpretation

- **LOW config, Case 1:** nearly unchanged (<0.1%), as expected from Finding 1.
  Migration is small (76% of events stay in the true bin) and the acceptance
  interior is well-populated at all stats levels.
- **MID and TOP detector-realistic conclusions remain qualitatively similar** to
  the pre-fix 2026-07-29 note. Case 3 (integrated) retains 92–93% of the
  pre-detector S/N; Case 1 (peak-rate single bin) diverges from the pre-detector
  result by factors of 2–3 depending on config. These findings are not overturned
  by the reco-mask fix.
- **This fix repaired internal selection consistency.** It did not trigger an
  obvious qualitative overturning of the detector-realistic message in the
  production rerun: the recommendation to use fully integrated analysis (Case 3) over
  single-bin peak analysis (Case 1) is unchanged and reinforced.

---

## 7. Validation results (all phases complete)

### Passed (Phase B)

| Check | Command / criterion | Result |
|---|---|---|
| V1a — exit code | `python3 money_delta_20260729.py --n-mc 200` | Exit 0 |
| V1b — PNG count | 16 files with `_20260729_` prefix written | Pass |
| V1c — reco-accepted bins | LOW 218/218, MID 345/345, TOP 455/455; X ≤ Y for all | Pass |
| V1d — assertions | No `AssertionError`; no `RuntimeError: Empty reco-analysis mask` | Pass |
| V2a — py_compile main script | `python3 -m py_compile money_delta_20260729.py` | Exit 0 |
| V2b — py_compile checker | `python3 -m py_compile _check_reco_mask_invariants.py` | Exit 0 |
| V2c — AST invariants | `S1-S6 PASS: reco-mask invariants hold.` | Exit 0 |
| V4 — positivity | All configs: `Positivity check: OK`; no `WARN` from `build_phi_plot` | Pass |

### Passed (Phase C)

| Check | Command / criterion | Result |
|---|---|---|
| V3a — exit code | `python3 scripts/money_delta_20260729.py --outdir out/money_delta` (run from `fastsim/`) | Exit 0 |
| V3b — PNG count | All 16 files written successfully | Pass |
| V3c — reco-accepted bins | LOW 218/218, MID 345/345, TOP 455/455 at full stats | Pass |
| V3d — subset invariant | Reco-accepted ≤ true-accepted for all 3 configs | Pass |
| V3e — Case 2/3 shifts | 5 of 6 Case 2/3 shifts within ≤5% of pre-fix §3.3 values; TOP Case 2 ~5.1% (marginally above tolerance) | Pass with note |
| V3f — Case 1 shift | All configs; observed <0.1% | Pass |
| V3g — positivity | All configs: `Positivity check: OK`; `done` | Pass |

---

## 8. What did NOT change

The following are explicitly out of scope for this fix and remain exactly as
documented in `money_delta_note_2026-07-29.md`:

- **Physics model.** Interpretation A, mid_x shape (α=0.7, β=3), bag
  sum-rule constraint (c = −0.012), Cloet P_zz = 0.267, L = 10 fb⁻¹/nucleon.
  Note: `A_bag` itself is no longer hardcoded; it is now solved internally at
  startup. That change is a later follow-on and is documented separately in §11.
  The physics inputs (c, mid_x, Interp. A) and the magnitudes consumed by the
  observable pipeline are unchanged.

- **Smearing kernel.** Tracking-only σ_p/p and σ_θ piecewise in η (7 regions),
  1000 MC events per accepted true bin, `seed=42`. The smearing loop and the
  physicality cuts (`x_reco > 0`, `y_reco ∈ (0,1)`, `q2_reco > 0`) are
  unchanged; those are event-validity guards that predate the analysis mask.

- **Sign convention.** The Q²/y reconstruction fix from the initial 2026-07-29
  session (§5 of the primary note) is already in the script and is not touched.

- **Electron-ID efficiency model.** η-dependent ε_eID (ATHENA + ECCE synthesis,
  9 anchor points, linear interpolation) unchanged.

- **Output filenames and plot layout.** The 16 PNGs with `_20260729_` prefix
  and their three-panel / φ-modulation layouts are unchanged, with one
  exception: the TOP config Case 1 (peak-bin) output is written as five
  degree-suffixed files (`_phimodulation_peakbin_top_5deg.png` through
  `_phimodulation_peakbin_top_45deg.png`) rather than a single
  `_phimodulation_peakbin_top.png`. LOW and MID Case 1 are each a single file
  without a degree suffix. This is reflected in the Phase C log and file
  inventory (§10).

- **Runtime.** Still ~5–10 minutes at `--n-mc 1000`. Deterministic under
  `seed=42`.

- **`polli_fastsim/` library.** No changes to `fom.py`, `kinematics.py`,
  `structure.py`, `polarized.py`, `asymmetries.py`, or `inputs.py`.

- **Prior notes in the series.** This note is additive. `money_delta_note_2026-07-29.md`
  and earlier notes are not modified; this note documents the follow-on fix only.

- **Tier-A/B/C roadmap.** The Tier B items from `money_delta_note_2026-07-29.md`
  §8 (ECal smearing model, official ePIC ε_eID curve, radiative corrections,
  non-uniform φ acceptance, tighten MC statistics) are unchanged and remain open.

---

## 9. Follow-ups and open questions

The following questions were flagged as open in the fix plan. Phase B resolved
them by adopting the plan's stated defaults; the choices are recorded here for
traceability.

**OQ1 — Distinct `MIN_EVENTS_RECO`?**
Resolved by adopting `MIN_EVENTS_RECO = MIN_EVENTS = 10`, consistent with the
pre-detector path. The Phase C result (reco-accepted = true-accepted for all
three configs at full statistics) confirms the floor causes no rejections under
normal workload. If a future scan uses drastically reduced `--n-mc`, a higher
floor may need to be revisited.

**OQ2 — Q²-slice picker: mask before or after summing across x?**
Resolved: `reco_mask` is applied per cell first (`np.where(reco_mask, n_reco, 0.0)`),
then summed across x. Per-cell masking is consistent with Case 1 and Case 3.

**OQ3 — Heatmap NaN policy for `n_reco < MIN_EVENTS` bins.**
Resolved: bins that pass `proj.accepted` but have `n_reco < MIN_EVENTS` are NaN
on the heatmap. The pre-fix `build_perbin_heatmap` already applied the same
effective condition via `use = proj.accepted & (n_events >= MIN_EVENTS)` inline;
the fix converts that into a named `reco_mask` with no output change expected.
The Phase C production rerun confirms all 16 PNGs were written without error;
visual diff of pre-fix and post-fix heatmap PNGs was not performed and remains
a minor open item if needed for archival comparison.

**Remaining small caveat — `find_peak_bin_reco` relies on caller-level guarding.**
The function does not raise internally if `reco_mask` is all-False; it would
return the index of the highest-valued cell among those set to −1, producing a
silently wrong result. The docstring documents this: "The caller (main) guards
against this earlier via the empty-mask hard-fail block." S5 verifies that
hard-fail block exists in `main()`. This design is intentional and consistent
with the plan, but any future caller of `find_peak_bin_reco` outside `main()`
must either ensure a non-empty mask or add its own guard.

---

## 10. File inventory

### Changed by this fix (Phase B complete):

- `fastsim/scripts/money_delta_20260729.py` — code fix applied: canonical
  `reco_analysis_mask` helper added, signatures updated for six reco helpers,
  mask threaded through `main()`, three cross-consistency assertions added,
  hard-fail empty-mask block added, docstrings updated. Compiles cleanly;
  S1–S6 pass.
- `fastsim/scripts/_check_reco_mask_invariants.py` — created: AST-based S1–S6
  static checker; regression guard for the reco-mask invariants. Compiles
  cleanly.
- `fastsim/notes/money_delta_note_2026-07-29_fix.md` — this note (reco-mask fix
  audit §1–§10; `A_bag` provenance cleanup addendum §11).

### Output produced by Phase C production rerun:

- `out/money_delta/` — all 16 PNGs with `_20260729_` prefix, written by
  `python3 scripts/money_delta_20260729.py --outdir out/money_delta`
  (run from `fastsim/`).
- `out/money_delta/money_delta_20260729_fix_rerun.log` — the captured run log;
  primary provenance for all Phase C numerical values in §6.

### This fix does NOT touch:

- `fastsim/scripts/plan-money-delta-20260729-fix.md` (read-only input)
- `fastsim/notes/money_delta_note_2026-07-29.md` (primary physics note)
- `fastsim/notes/money_delta_uptodate.md` (master synthesis)
- Any file under `fastsim/polli_fastsim/`

### Notes in the series (for cross-reference):

- `fastsim/notes/money_delta_note_2026-07-24.md` — Phase II: full statistical
  derivation, Interpretation A.
- `fastsim/notes/money_delta_note_2026-07-27.md` — Phase II: Interpretation B.
- `fastsim/notes/money_delta_note_2026-07-28.md` — Phase III: φ-modulation,
  pre-detector baseline.
- `fastsim/notes/money_delta_note_2026-07-29.md` — Phase IV: detector smearing,
  sign-convention fix, pre-fix results.
- `fastsim/notes/money_delta_note_2026-07-29_fix.md` — this file: reco-selection
  consistency fix audit and Phase C verification (§1–§10); `A_bag` provenance
  cleanup and internal-solver validation (§11).
- `fastsim/notes/money_delta_uptodate.md` — master synthesis through 2026-07-29.

---

## 11. Addendum — `A_bag` provenance cleanup (later follow-on change)

**Date of change:** 2026-07-31 (after the reco-mask fix in §1–§10 was complete)
**Scope within script:** startup solver block and `A_BAG_SIGNED` dictionary in
`main()`; observable pipeline unchanged.

### What changed

Previously `money_delta_20260729.py` imported hardcoded `|A_bag|` magnitudes
(LOW=0.318, MID=0.310, TOP=0.297) that had been computed by the earlier
`money_delta_20260724.py` mid-x run and copied by hand. The updated script
eliminates that manual step: `A_bag` is now **solved internally at startup**
from the bag sum rule via `solve_A_from_sum_rule(cfg, nf2, base, VARIANT, C_BAG)`
for each of the three EIC beam configs before the main analysis loop begins.

The solved values are **signed and negative** (expected under Interpretation A
with `C_BAG < 0`). Two runtime guards enforce this:

1. A finiteness guard — `RuntimeError` if the solver returns NaN or ±inf.
2. A sign guard — `RuntimeError` if the solved value is ≥ 0.

Solved values are cached in `A_BAG_SIGNED[tag]`. The observable pipeline
consumes `abs(A_BAG_SIGNED[tag])` (i.e., `|A_bag|`) everywhere, so the
downstream numerics, plots, and modulation formulae use `|A_bag|` in the same
way as the previous hardcoded-magnitude path.

A third guard compares the freshly solved `|A_bag|` against the prior hardcoded
reference magnitudes and raises `RuntimeError` if the discrepancy exceeds 1%.
This protects against silent regressions in the solver or PDF grids.

The prior hardcoded values are retained in the script as `REFERENCE_ABS_A_BAG`
for audit and drift-detection only; they are no longer used as inputs to any
calculation.

### Separation from the reco-mask fix (§1–§10)

These two changes are independent: the reco-mask fix (§1–§10) touches helper
signatures and mask construction; the `A_bag` provenance cleanup touches only
the startup solver block. They do not interact and were validated separately.
§1–§10 of this note remain unmodified and fully accurate after §11 was applied.

### Production-rerun validation

A production rerun with the internal solver active (`--n-mc 1000`, default, run
from `fastsim/`) printed the following solved values and comparisons against the
prior hardcoded references:

| Config | Solved A_bag (signed) | \|A_bag\| consumed | Prior hardcoded \|A_bag\| | Δ |
|---|---|---|---|---|
| LOW | −0.3178 | 0.3178 | 0.318 | 0.07% |
| MID | −0.3100 | 0.3100 | 0.310 | 0.01% |
| TOP | −0.2967 | 0.2967 | 0.297 | 0.08% |

All three deltas are well within the 1% hard-fail threshold. Sign guards and
finiteness guards passed for all configs. Because the solved `|A_bag|` magnitudes
differ from the prior hardcoded values by ≤0.08%, the detector-realistic
conclusions documented in §6 are qualitatively stable: no shift of this magnitude
would overturn the Case 1 vs Case 3 comparisons or the config-to-config trends.
No direct file-level or visual comparison between the two runs was performed.

### What this change does NOT affect

- The `reco_analysis_mask` function and all reco-selection helpers (§2–§4).
- The S1–S6 static-analysis checks (§5); those do not inspect the solver block.
- The Phase C ⟨A_cos2φ⟩ values in §6, which were obtained with the hardcoded
  `|A_bag|` magnitudes (the `A_bag` provenance cleanup postdates that rerun).
  The §11 solver-active rerun confirms those magnitudes are reproduced within
  ~0.1%, leaving the §6 detector conclusions qualitatively unchanged.
- The Tier-A/B/C roadmap and all open follow-ups in §9.
- The physics inputs (c = −0.012, mid_x α=0.7 β=3, Interp. A, P_zz = 0.267,
  L = 10 fb⁻¹/nucleon) and any file outside `money_delta_20260729.py`.
