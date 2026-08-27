# Plan: Detector-realistic money-delta script — reco-selection consistency fix

> **Dated record.** Numbers here predate the 2026-08-27 corrections: ⁶Li
> beam energies are γ-matched (40.8 / 99.5 / 137.5 GeV/u), not
> rigidity-scaled (20.5 / 50 / 137.5), and the far-forward divergence is
> per-configuration rather than a single 72.7 μrad. See plans/10. Kept as
> the record of its date, not as a current reference.


**Scope**: Fix reconstructed-event selection in `fastsim/scripts/money_delta_20260729.py` so that all reco-side outputs (peak-bin find, Q²-slice find, integrated averages, per-bin heatmaps, φ-modulation plots, diagnostics) use a single, well-defined **reconstructed-analysis mask** that reapplies the same acceptance used at the true level, and that the true-level `proj.accepted` mask is never silently used as a *reco* selection where a reco mask is meant. Document the diagnosis and fixes in a new note `fastsim/notes/money_delta_note_2026-07-29_fix.md`.

**Files touched**:
- `fastsim/scripts/money_delta_20260729.py` — code fix (delegated to `script-implementer`).
- `fastsim/scripts/_check_reco_mask_invariants.py` — new AST-based static checker; regression guard for the reco-mask invariants (delegated to `script-implementer`; specification in this plan).
- `fastsim/notes/money_delta_note_2026-07-29_fix.md` — new note summarizing bugs found, definition of the reco-analysis mask, and post-fix behavior (delegated to `doc-writer`).
- (This plan) `fastsim/scripts/plan-money-delta-20260729-fix.md`.

**Not touched**:
- `polli_fastsim/` library. No changes required; `fom.Scenario` and `kinematics.kinematic_mask` are used as-is.
- Prior notes (`money_delta_note_2026-07-29.md` and earlier). The new note is additive.
- Prior script snapshots (`money_delta_20260728.py`, etc.).

**Annotation cycle**: iteration 4. Prior fixes retained: (1) removed nonexistent `scenario.x_max`; (2) reco mask reuses full `proj.accepted` (DIS + η + E'); (3) hard-fail on empty reco mask; (4) AST-based static checker `_check_reco_mask_invariants.py` (S1–S6). This iteration resolves three internal inconsistencies: (a) the empty-mask `raise RuntimeError(...)` must embed the message string directly on the Raise node (not via a `msg` variable) so S5's AST scan can see it; the S5 spec and skeleton are simplified accordingly. (b) `print_smearing_diagnostics` cannot compute `true-accepted` count internally without accessing `proj.accepted`, which S1 forbids — the helper now takes `true_accepted_count: int` as an explicit caller-supplied argument. (c) One stale "grep-based" phrase in the Requirements Summary updated to "AST-based".

**Input sources** (read while drafting):
- Code review findings supplied by the user (four bullets, first pass + critique).
- Direct reads of `fastsim/scripts/money_delta_20260729.py` (whole file).
- Direct read of `fastsim/polli_fastsim/kinematics.py` (`kinematic_mask`, `scattered_electron`, `log_grid`).
- Direct read of `fastsim/polli_fastsim/fom.py` (`Scenario` dataclass fields, `project_rates` construction of `proj.accepted`).
- `fastsim/notes/money_delta_note_2026-07-29.md` (context, prior sign-fix, terminology).

---

## Requirements Summary

**Confirmed requirements**
- The script must produce reconstructed-level results (peak-bin, Q²-slice, integrated, per-bin heatmaps, φ plots, smearing diagnostics) that are **internally consistent**: every reco output is derived from a single, explicitly-defined reconstructed-analysis mask.
- The reco-analysis mask must reapply **exactly the same acceptance** used to build `proj.accepted` (see §"Acceptance decision" below), plus a reco-side statistics floor (`n_reco >= MIN_EVENTS`).
- The true-level `proj.accepted` mask must not be used as a *reco selection* anywhere in the analysis flow. Where truth-level acceptance is genuinely needed (driving the smearing loop over source bins, computing the truth-side comparison A_true), it stays; where reco selection is intended, the reco mask is used.
- A successful full run (with `--n-mc 1000`) must write **all 16 expected PNGs** — this is the success invariant. No config may be silently skipped.
- Case 3 (fully integrated) reco results must remain close to the pre-detector integrated results (§3.3 of the 2026-07-29 note reports 92–93% S/N retention). The fix should not break this cross-check.
- Runtime characteristics unchanged: still ~5–10 min at `--n-mc 1000`, deterministic under `seed=42`.
- Filenames and plot layout unchanged. Same 16 PNGs with `_20260729_` prefix.
- Diagnostics printed to stdout must clearly distinguish true-accepted bin counts from reco-accepted bin counts.
- A static-analysis check (AST-based, run at end of Phase B) confirms that no reco-side helper still references `proj.accepted` or `n_reco > 0` as a mask.
- A companion note `fastsim/notes/money_delta_note_2026-07-29_fix.md` documents (a) each bug found, (b) the reco-analysis-mask definition and the acceptance decision, (c) what changed numerically, (d) which sanity checks pass post-fix.

**Inferred assumptions** (flag if wrong)
- The intended reco-analysis mask is `proj.accepted & (n_reco >= MIN_EVENTS)` — see "Acceptance decision" below for the derivation.
- The `MIN_EVENTS = 10` floor from the parent script should be reused for the reco mask (matches the true-level convention already present in `sig2_per_fb_at_sumrule` at line 406 of the current script). If the user prefers a distinct `MIN_EVENTS_RECO`, flag it (Open Question 1).
- All kinematic-cut and acceptance parameters are drawn from the single `Scenario` object built in `main()` (currently `Scenario(lumi_fb_per_nucleon=LUMI_FB, pol_ion_tensor=PZZ, q2_min=2.0)`; all other cut fields at their `Scenario` defaults). No new CLI flags introduced.
- "Reconstructed-analysis mask" is a **bin-level** mask (shape `(nx, nq2)`). Per-MC-event physicality cuts inside the smearing loop (`x_reco > 0`, `y_reco ∈ (0,1)`, `q2_reco > 0`) remain unchanged — those are event-validity guards, not acceptance.

**Resolved defaults** (previously ambiguous, now locked)
1. Reco cut source = the `Scenario` object already constructed in `main()`. No new configuration surface.
2. `MIN_EVENTS_RECO = MIN_EVENTS = 10`. Consistent with the pre-detector path.
3. The mask is computed **once per config** in `main()` after `smear_config` returns and passed by argument into every reco-side helper. No re-derivation inside helpers.
4. Empty reco-mask handling = **hard failure**, not silent skip (see "Failure mode for empty reco mask" below).

## Acceptance decision (blocking fix #2 — explicit)

`proj.accepted` in `fom.project_rates` (`fastsim/polli_fastsim/fom.py:65-74`) is built as:

```python
mask  = kin.kinematic_mask(X, Q2, s,
            q2_min=scenario.q2_min, y_min=scenario.y_min,
            y_max=scenario.y_max, w2_min=scenario.w2_min)          # DIS cuts
mask &= (eta >= scenario.eta_min) & (eta <= scenario.eta_max)       # eta acceptance
mask &= e_p >= scenario.e_prime_min                                 # E' threshold
```

The `Scenario` dataclass (`fom.py:31-44`) exposes exactly these fields:
`q2_min`, `y_min`, `y_max`, `w2_min`, `eta_min`, `eta_max`, `e_prime_min`
(**no `x_max` field** — the prior draft incorrectly referenced `scenario.x_max`; `kinematic_mask` has an `x_max` kwarg but `project_rates` never sets it, so the default 1.0 is always in force and it is not part of the `Scenario` surface).

**Decision — the reco-analysis mask uses the FULL acceptance, not just DIS cuts.** Definition:

```
reco_mask = proj.accepted & (n_reco >= MIN_EVENTS)
```

Rationale (four points):
1. The bin geometry is shared: the smearing loop writes `n_reco[ir, jr]` into the same `(nx, nq2)` grid whose centers are `proj.x`, `proj.q2`. Bin-center kinematics are identical for both masks; any bin that passes the reco acceptance test at bin center is exactly a bin in `proj.accepted`.
2. Reapplying η and E' cuts at *reco* bin centers is the correct physical criterion for a bin-level analysis: an experiment does not analyze bins outside the detector's kinematic acceptance regardless of whether events happened to migrate into them.
3. Bins that received migrated reco events but were not in `proj.accepted` must therefore be **excluded** from the reco analysis. The `n_reco` counts in those bins are physically real (events did arrive), but the bin is not part of the analysis region. Those events contribute to the `n_reco_total` yield sum (already tracked in `diag`) but not to reco averages or heatmaps.
4. Building the mask as `proj.accepted & (n_reco >= MIN_EVENTS)` (rather than recomputing `kinematic_mask + η + E'` locally from `proj.extras`) has three concrete benefits: (a) it cannot drift if `Scenario` gains new acceptance fields, (b) it makes the "reco mask ⊆ true mask" invariant lexically obvious to the reader, (c) it removes any need for the reco path to know the internal composition of `proj.accepted`.

Consequence for stdout diagnostics: `n_reco_acc_bins ≤ n_true_acc_bins` is an invariant. When the two differ, the gap is entirely from the `MIN_EVENTS` floor operating on the stochastic `n_reco`.

## Failure mode for empty reco mask (blocking fix #3 — explicit)

The success invariant is: a successful full run writes all 16 PNGs. A config with zero reco-accepted bins violates this. Given the parameters in the current script (three well-populated EIC configs, `--n-mc >= 200`), this should never happen for the intended workload; if it does, it indicates a real problem (mis-specified scenario, broken smearing, or drastically reduced luminosity), and silent skip would hide the bug and produce an incomplete output set. Therefore:

**On empty reco mask (`reco_mask.sum() == 0`) for any config**, the script must:
1. Print a clear, prominent error to `stderr` naming the config, the true-accepted bin count, and the reco `n_reco_total` for diagnosis.
2. Raise `RuntimeError` with the same message, halting the run *before* any downstream helper attempts to index into an empty selection.
3. Not attempt to produce partial output for the offending config or for later configs — the invariant is all-or-nothing.

Rationale: the empty-mask condition is inconsistent with the physics of any realistic EIC beam config × `n_mc >= 200` and indicates a real failure the user must see. The existing `--n-mc 200` smoke test in the validation section will exercise the mask construction end-to-end for all three configs and detect any surprise.

If a future workflow legitimately needs a "skip empty configs" mode (e.g., a scan over many hypothetical scenarios), that becomes an explicit CLI flag (`--allow-empty-reco-mask`) and a separate change — out of scope for this fix.

## Open Questions

1. **Distinct `MIN_EVENTS_RECO`?** The pre-detector path uses `MIN_EVENTS = 10`. The reco path applies smearing + efficiency, so `n_reco` is stochastic. Default: reuse `MIN_EVENTS = 10`. Flag if the user wants `MIN_EVENTS_RECO = 20` or `MIN_EVENTS_RECO = 5 · (1/⟨ε⟩)` to preserve pre-fix behavior for LOW/MID.
2. **Q²-slice picker: mask before or after summing across x?** Default: apply the reco mask per cell first (`np.where(reco_mask, n_reco, 0.0)`), then sum across x. Flag if the user prefers requiring the *slice sum* to pass a floor instead.
3. **Heatmap NaN policy.** Bins that pass acceptance but have `n_reco < MIN_EVENTS` will be NaN under the new mask (previously they may have shown as very-low-count cells). Confirm default: hide them.

---

## Approach

### 1. Introduce a canonical reco-analysis-mask helper

Add a single module-level helper. This replaces the four ad-hoc uses of either `proj.accepted` or `n_reco > 0.0` scattered across the reco path.

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

    Note: `proj.accepted` alone is the TRUE-level acceptance; do NOT
    apply it to reco arrays without also enforcing the reco statistics
    floor. This function is the single source of truth for reco selection.

    Invariant: reco_mask ⊆ proj.accepted (asserted by callers if desired).
    """
    return proj.accepted & (np.asarray(n_reco) >= float(min_events))
```

**Placement.** Right after the `MIN_EVENTS = 10` constant (~line 262) to keep related definitions together.

### 2. Rewrite each affected helper to take the reco mask by argument

Every existing helper below currently accepts either `proj` (and uses `proj.accepted`) or `accepted` (silently the true-level mask). Change signatures to accept an explicit `reco_mask` argument. Callers pass the result of `reco_analysis_mask(...)`.

| Helper | Current (buggy) behavior | Fix |
|---|---|---|
| `find_peak_bin_reco(n_reco, accepted)` | `n[~accepted] = -1.0` where `accepted = proj.accepted` (true-level) | Signature → `find_peak_bin_reco(n_reco, reco_mask)`; mask by `reco_mask`. If no cell passes, raise `RuntimeError` (caller-level empty-mask guard catches this earlier — see step 4). |
| `find_peak_q2_slice_reco(n_reco, accepted)` | `np.where(accepted, n_reco, 0.0).sum(axis=0)` — true mask on reco counts | Signature → `find_peak_q2_slice_reco(n_reco, reco_mask)`; use `reco_mask` in the `np.where`. |
| `compute_A_cos2phi_q2slice_reco(n_reco, A_reco, accepted, iq2)` | `valid = n_slice > 0.0` (no analytic cuts; `accepted` arg unused/misleading) | Signature → `compute_A_cos2phi_q2slice_reco(n_reco, A_reco, reco_mask, iq2)`; `valid = reco_mask[:, iq2]`. |
| `compute_A_cos2phi_integrated_reco(n_reco, A_reco)` | `valid = n_reco > 0.0` (no analytic cuts, no floor) | Add `reco_mask` argument: `compute_A_cos2phi_integrated_reco(n_reco, A_reco, reco_mask)`; `valid = reco_mask`. |
| `build_perbin_heatmap(..., proj, n_reco, A_reco)` | `use = proj.accepted & (n_events >= MIN_EVENTS)` where `n_events = n_reco` — same idea as the correct mask, but expressed inline and inconsistently with the other helpers | Add `reco_mask` argument; replace inline expression with `use = reco_mask`. |
| `print_smearing_diagnostics(...)` | Calls `find_peak_bin_reco(n_reco, proj.accepted)` internally | Add two explicit arguments: `reco_mask` (forwarded to `find_peak_bin_reco`) and `true_accepted_count: int` (passed by caller). Add a line: `Reco-accepted bins: X (vs true-accepted: Y)` using `int(reco_mask.sum())` for X and the caller-supplied `true_accepted_count` for Y. **The helper must not access `proj.accepted` directly** — that would violate S1. Both counts are precomputed in `main()` and threaded in as ints. |

**Not changed**: `find_peak_bin(proj)`, `compute_A_cos2phi_at_bin(...)`, `find_peak_q2_slice(proj)`, `compute_A_cos2phi_q2slice(...)`, `compute_A_cos2phi_integrated(...)`, and the `sig2_per_fb_at_sumrule` helper. These are the true-level (unsmeared) analogues and correctly use `proj.accepted`. They are still called for the `A_true_peakbin` comparison in smearing diagnostics.

### 3. Compute the mask once per config in `main()` and thread it through

In the `for cfg, config_tag, config_label in all_configs:` loop, right after `smear_config` returns (~line 1624), before any reco-side helper is called:

```python
# Reconstructed-analysis mask: the SINGLE source of truth for reco selection.
# = proj.accepted (full DIS + eta + E' acceptance) AND n_reco >= MIN_EVENTS.
# See reco_analysis_mask docstring for the acceptance decision.
reco_mask = reco_analysis_mask(proj, n_reco, min_events=MIN_EVENTS)
n_reco_acc_bins = int(reco_mask.sum())
print(f"  Reco-accepted bins: {n_reco_acc_bins}  "
      f"(vs true-accepted: {n_acc_bins})")

# Fail loudly if the mask is empty — see plan §"Failure mode for empty reco mask".
# NOTE: the f-string must be the DIRECT argument of RuntimeError(...) — not
# assigned to a `msg` variable first — so that S5's AST check can see the
# literal "Empty reco-analysis mask" substring on the Raise node itself.
if n_reco_acc_bins == 0:
    print(
        f"ERROR: Empty reco-analysis mask for config {config_tag!r}: "
        f"true-accepted bins = {n_acc_bins}, "
        f"n_reco_total = {diag['n_reco_total']:.3e}.",
        file=sys.stderr,
    )
    raise RuntimeError(
        f"Empty reco-analysis mask for config {config_tag!r}: "
        f"true-accepted bins = {n_acc_bins}, "
        f"n_reco_total = {diag['n_reco_total']:.3e}. "
        f"Cannot produce the expected PNG set. "
        f"Check --n-mc, Scenario cuts, or smearing kernel."
    )
```

Then update every downstream call site inside the loop to pass `reco_mask`:
- `print_smearing_diagnostics(..., reco_mask=reco_mask, true_accepted_count=n_acc_bins)` — pass both counts explicitly; the helper is forbidden from reading `proj.accepted` (S1).
- `build_perbin_heatmap(..., reco_mask=reco_mask)`
- `find_peak_bin_reco(n_reco, reco_mask)` (Case 1)
- `find_peak_q2_slice_reco(n_reco, reco_mask)` (Case 2)
- `compute_A_cos2phi_q2slice_reco(n_reco, A_reco, reco_mask, iq2_peak_slice_r)` (Case 2)
- `compute_A_cos2phi_integrated_reco(n_reco, A_reco, reco_mask)` (Case 3)

### 4. Cross-consistency assertions inside `main()`

After the mask is built, add three cheap invariants (assertion messages are informational; failures indicate a real regression):

```python
assert reco_mask.shape == n_reco.shape == A_reco.shape, (
    "shape mismatch: reco_mask / n_reco / A_reco must be identical")
assert not np.any(reco_mask & (n_reco < MIN_EVENTS)), (
    "reco_mask leaked a bin below MIN_EVENTS floor")
assert bool(np.all(reco_mask <= proj.accepted)), (
    "reco_mask must be a subset of proj.accepted (acceptance invariant)")
```

### 5. Update docstrings

- `smear_config` docstring: clarify that `n_reco` is populated for **every** in-grid reco bin regardless of the analysis mask, and that the analysis mask is applied downstream by `reco_analysis_mask`.
- Update all six helper docstrings to explain the `reco_mask` argument and cross-reference `reco_analysis_mask`.
- Update the module docstring: add one paragraph explaining the reco-analysis-mask definition and pointing to the note.

### 6. Diagnostic-print updates

- `print_smearing_diagnostics` gains one line, `Reco-accepted bins: X (vs true-accepted: Y)`, where X = `int(reco_mask.sum())` and Y = the `true_accepted_count` int passed in by the caller (computed once in `main()` as `int(proj.accepted.sum())` and stored in `n_acc_bins`). The helper itself never touches `proj.accepted` — required by S1.
- Per-config summary block in `main()` also prints `n_reco_acc_bins` alongside `n_acc_bins` right after the mask is built (before calling `print_smearing_diagnostics`).
- Case 1 / 2 / 3 blocks: no layout change, only input mask changes.

---

## Static-analysis validation (blocking fix #4)

At the end of Phase B, run a single Python AST-based checker that parses `money_delta_20260729.py`, walks each named function's body, and verifies six invariants (S1–S6). AST-based checking is chosen over grep/awk because:
- `awk '/^def foo/,/^def /'` cannot reliably delimit a function body (comments, decorators, nested defs, and the trailing `def` line pollute the range).
- Regex cannot distinguish real attribute access from mentions inside docstrings or comments.
- Legitimate occurrences (the mask helper itself; the true-level `sig2_per_fb_at_sumrule` control; the main-loop line that builds the mask) need to be **allow-listed by function name**, not by text pattern.
- The checker must exit nonzero on failure so it can act as a Phase B gate and, later, a regression guard.

**Deliverable**: a standalone script `fastsim/scripts/_check_reco_mask_invariants.py`. Not a plan artifact; created by `script-implementer` as part of Phase B. Runs as:

```bash
python3 fastsim/scripts/_check_reco_mask_invariants.py \
    fastsim/scripts/money_delta_20260729.py
```

Exit code 0 on all-pass, 1 on any failure. Stdout lists each check with `PASS`/`FAIL` and, for failures, the offending line numbers.

### Checker specification (implementer contract)

The checker parses the target file with `ast.parse`, indexes top-level `FunctionDef` nodes by name, and applies the six checks below. **Docstrings** (the first `ast.Expr(Constant(str))` of a function body) are excluded from all attribute/comparison scans. **Comments** are naturally excluded because AST discards them.

Define these name sets at the top of the checker:

```python
RECO_HELPERS = {
    "find_peak_bin_reco",
    "find_peak_q2_slice_reco",
    "compute_A_cos2phi_q2slice_reco",
    "compute_A_cos2phi_integrated_reco",
    "build_perbin_heatmap",
    "print_smearing_diagnostics",
}
TRUE_LEVEL_HELPERS_KEEPING_PROJ_ACCEPTED = {
    "sig2_per_fb_at_sumrule",
    "find_peak_bin",
    "find_peak_q2_slice",
    "compute_A_cos2phi_q2slice",
    "compute_A_cos2phi_integrated",
}
MASK_BUILDER = "reco_analysis_mask"
```

#### S1 — no reco helper reads `proj.accepted` or `<any>.accepted` from a passed-in projection
For each function `fn` whose name is in `RECO_HELPERS`, walk `fn.body` (excluding the docstring node) and collect every `ast.Attribute` node with `attr == "accepted"`. **Pass** iff the list is empty. On fail, emit `S1 FAIL: <fn> line <N>: attribute '.accepted' accessed`.

#### S2 — no reco helper uses `n_reco > 0`, `n_reco > 0.0`, or `n_slice > 0[.0]` as a routing condition
For each `fn` in `RECO_HELPERS`, walk the body and find every `ast.Compare` node where:
- `left` is `ast.Name` with `id in {"n_reco", "n_slice", "n_bin_map"}` **or** an `ast.Subscript` of such a name,
- the sole `op` is `ast.Gt`,
- the sole `comparator` is `ast.Constant` with `value in {0, 0.0}`.

Pass iff no such node is found. On fail, emit `S2 FAIL: <fn> line <N>: bare "<expr> > 0" comparison`.

Rationale: `reco_mask[:, iq2]` and `reco_mask` return bools already; any surviving `n_* > 0` in a routing position is the old bug pattern.

#### S3 — every reco helper has `reco_mask` in its parameter list
For each `fn` in `RECO_HELPERS`, check that `"reco_mask"` appears in `[a.arg for a in fn.args.args + fn.args.kwonlyargs]`. Pass iff true for all six. On fail, emit `S3 FAIL: <fn> missing reco_mask parameter`.

#### S4 — `reco_analysis_mask` is defined once and called exactly once in `main`
- Verify exactly one top-level `FunctionDef` named `reco_analysis_mask`.
- Locate top-level `FunctionDef` named `main` (or `_main` fallback). Walk `main.body` and count `ast.Call` nodes whose `func` is `ast.Name(id="reco_analysis_mask")`. Pass iff the count is exactly 1.

On fail, emit `S4 FAIL: reco_analysis_mask defined N times, called M times in main() (expected 1 each)`.

#### S5 — the empty-mask failure path exists
Walk `main.body` and find every `ast.Raise` node whose `exc` is a `Call` of `Name(id="RuntimeError")`. Pass iff at least one such node has an argument that is either an `ast.Constant(str)` or an `ast.JoinedStr` (f-string) whose literal string parts, concatenated, contain the substring `"Empty reco-analysis mask"`. On fail, emit `S5 FAIL: no "Empty reco-analysis mask" RuntimeError in main()`.

**Implementation contract this enforces on the raise site**: the message must be embedded directly as a string literal or f-string in the `raise RuntimeError(...)` call, not assigned to an intermediate variable (e.g. `msg = ...; raise RuntimeError(msg)` is NOT sufficient because `msg` is an `ast.Name`, not a string node). See the code fragment in §"Failure mode for empty reco mask" for the required form.

#### S6 — true-level helpers still reference `proj.accepted` (deliberate control)
For each `fn` in `TRUE_LEVEL_HELPERS_KEEPING_PROJ_ACCEPTED`, verify at least one `ast.Attribute` with `attr == "accepted"` in the body. Pass iff true for every listed helper. On fail, emit `S6 FAIL: <fn> no longer references .accepted — fix over-reached`.

### Checker skeleton (reference implementation for `script-implementer`)

Approximately 120 lines. Uses only stdlib (`ast`, `sys`). No third-party dependencies.

```python
#!/usr/bin/env python3
"""Static-analysis regression guard for the reco-selection consistency fix.

Verifies six invariants (S1-S6) on money_delta_20260729.py.
Exit code 0 on all-pass, 1 on any failure. See
plan-money-delta-20260729-fix.md § "Static-analysis validation" for the
authoritative specification.
"""
import ast, sys

RECO_HELPERS = {
    "find_peak_bin_reco", "find_peak_q2_slice_reco",
    "compute_A_cos2phi_q2slice_reco", "compute_A_cos2phi_integrated_reco",
    "build_perbin_heatmap", "print_smearing_diagnostics",
}
TRUE_HELPERS = {
    "sig2_per_fb_at_sumrule", "find_peak_bin", "find_peak_q2_slice",
    "compute_A_cos2phi_q2slice", "compute_A_cos2phi_integrated",
}
BARE_NAMES = {"n_reco", "n_slice", "n_bin_map"}


def body_without_docstring(fn):
    body = fn.body
    if (body and isinstance(body[0], ast.Expr)
            and isinstance(body[0].value, ast.Constant)
            and isinstance(body[0].value.value, str)):
        return body[1:]
    return body


def walk_nodes(nodes):
    for n in nodes:
        yield from ast.walk(n)


def check(path):
    src = open(path).read()
    tree = ast.parse(src, filename=path)
    funcs = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    fails = []

    # S1
    for name in RECO_HELPERS:
        fn = funcs.get(name)
        if fn is None:
            fails.append(f"S1 FAIL: reco helper {name!r} not found")
            continue
        for node in walk_nodes(body_without_docstring(fn)):
            if isinstance(node, ast.Attribute) and node.attr == "accepted":
                fails.append(f"S1 FAIL: {name} line {node.lineno}: '.accepted' accessed")

    # S2
    for name in RECO_HELPERS:
        fn = funcs.get(name)
        if fn is None:
            continue
        for node in walk_nodes(body_without_docstring(fn)):
            if (isinstance(node, ast.Compare)
                    and len(node.ops) == 1 and isinstance(node.ops[0], ast.Gt)
                    and len(node.comparators) == 1):
                cmp = node.comparators[0]
                if isinstance(cmp, ast.Constant) and cmp.value in (0, 0.0):
                    left = node.left
                    if isinstance(left, ast.Name) and left.id in BARE_NAMES:
                        fails.append(f"S2 FAIL: {name} line {node.lineno}: {left.id} > {cmp.value}")
                    elif (isinstance(left, ast.Subscript)
                          and isinstance(left.value, ast.Name)
                          and left.value.id in BARE_NAMES):
                        fails.append(f"S2 FAIL: {name} line {node.lineno}: {left.value.id}[...] > {cmp.value}")

    # S3
    for name in RECO_HELPERS:
        fn = funcs.get(name)
        if fn is None:
            continue
        arg_names = [a.arg for a in fn.args.args + fn.args.kwonlyargs]
        if "reco_mask" not in arg_names:
            fails.append(f"S3 FAIL: {name} missing 'reco_mask' parameter")

    # S4
    defs = [n for n in tree.body
            if isinstance(n, ast.FunctionDef) and n.name == "reco_analysis_mask"]
    main = funcs.get("main")
    n_calls = 0
    if main is not None:
        for node in ast.walk(main):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Name)
                    and node.func.id == "reco_analysis_mask"):
                n_calls += 1
    if len(defs) != 1 or n_calls != 1:
        fails.append(f"S4 FAIL: reco_analysis_mask defs={len(defs)}, calls in main()={n_calls} (expected 1/1)")

    # S5: message must be a Constant(str) or JoinedStr directly on the Raise
    # (not an intermediate `msg` Name). See plan §"Failure mode" contract.
    def _literal_parts(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            yield node.value
        elif isinstance(node, ast.JoinedStr):
            for v in node.values:
                if isinstance(v, ast.Constant) and isinstance(v.value, str):
                    yield v.value

    ok_s5 = False
    if main is not None:
        for node in ast.walk(main):
            if (isinstance(node, ast.Raise)
                    and isinstance(node.exc, ast.Call)
                    and isinstance(node.exc.func, ast.Name)
                    and node.exc.func.id == "RuntimeError"):
                for arg in node.exc.args:
                    joined = "".join(_literal_parts(arg))
                    if "Empty reco-analysis mask" in joined:
                        ok_s5 = True
    if not ok_s5:
        fails.append("S5 FAIL: no RuntimeError with literal 'Empty reco-analysis mask' in main() "
                     "(the string must be embedded directly on the Raise, not via an intermediate variable)")

    # S6
    for name in TRUE_HELPERS:
        fn = funcs.get(name)
        if fn is None:
            continue  # helper may be absent; not a regression
        found = any(isinstance(n, ast.Attribute) and n.attr == "accepted"
                    for n in walk_nodes(body_without_docstring(fn)))
        if not found:
            fails.append(f"S6 FAIL: {name} no longer references .accepted (fix over-reached)")

    for f in fails:
        print(f)
    if fails:
        print(f"\n{len(fails)} check(s) failed.")
        return 1
    print("S1-S6 PASS: reco-mask invariants hold.")
    return 0


if __name__ == "__main__":
    sys.exit(check(sys.argv[1] if len(sys.argv) > 1 else
                   "fastsim/scripts/money_delta_20260729.py"))
```

Any nonzero exit from this checker blocks Phase B sign-off. Capture stdout to `/tmp/mdfix_static.log`; the doc-writer references it in §5 of the note.

---

## Trade-offs and considerations

- **Full-acceptance vs DIS-only reco mask.** Choosing `proj.accepted & (n_reco >= MIN_EVENTS)` rather than DIS-cuts-only means bins where events migrated *into* a region outside the detector's η or E' acceptance are excluded from reco averages. This is the physically correct choice for an experimental analysis; the alternative (DIS-only mask) would let migrated events into the analysis from bins the detector cannot actually measure. Documented explicitly in the note.
- **Numerical impact expected to be small for Case 3.** The reco mask is `proj.accepted ∩ {n_reco >= MIN_EVENTS}`. For truth-accepted bins with substantial rate, the second condition is nearly always met. Expect Case 3 numbers to shift by <1%.
- **Case 1 numerical impact is bounded but not zero.** The peak-bin is the largest-N reco cell within `reco_mask`. If pre-fix Case 1 selected a bin that was in `proj.accepted` (which most peak candidates were), the answer is unchanged. If migration created a larger-N bin *outside* `proj.accepted`, Case 1 pre-fix would have missed it anyway (because the pre-fix code already applied `proj.accepted` at line 837) — but the pre-fix code's use of the true-mask happened to coincide with the acceptance intent for Case 1 alone. For Case 2 and Case 3, the pre-fix code was silently wrong (Case 2 masked by `proj.accepted` incorrectly, Case 3 used no mask at all beyond `n_reco > 0`).
- **Q²-slice picker impact.** Same reasoning as Case 1, applied per-x column.
- **Heatmap changes.** Bins that are truth-accepted but received `< MIN_EVENTS` reco events will now be NaN (previously, per the current inline `use = accepted & (n_reco >= MIN_EVENTS)`, they were already NaN — so the *heatmap* is a pure refactor with no output change if the current code is correct there; the fix converts an inline mask into a named one).
- **No library changes.** `polli_fastsim/` untouched; `Scenario` fields consumed by name only.
- **Failure-mode strictness.** Hard failure on empty reco mask means the script is safer as an automation target (CI, cron): incomplete output = nonzero exit code = visible failure. Silent skip would produce partial output sets that could be mistaken for successful runs.

---

## Validation

After Phase B is implemented, run in order:

### V1 — Smoke test
```bash
python3 fastsim/scripts/money_delta_20260729.py \
    --outdir fastsim/out/money_delta --n-mc 200
```
Reduced `--n-mc 200` for a fast check (~1–2 min).

**Pass criteria**:
1. Script exits 0.
2. All **16 PNGs** exist under `--outdir`.
3. Per-config stdout shows `Reco-accepted bins: X (vs true-accepted: Y)` with `X > 0`, `X ≤ Y`, and `X` within ~30% of `Y` for a normal config at `--n-mc 200`.
4. Assertions in `main()` pass (no `AssertionError`).
5. No `RuntimeError: Empty reco-analysis mask` on any config.

### V2 — Static-analysis checks
```bash
python3 fastsim/scripts/_check_reco_mask_invariants.py \
    fastsim/scripts/money_delta_20260729.py
```
Pass criteria: exit code 0 and final line `S1-S6 PASS: reco-mask invariants hold.`

### V3 — Full-statistics regression
```bash
python3 fastsim/scripts/money_delta_20260729.py \
    --outdir fastsim/out/money_delta --n-mc 1000 2>&1 | tee /tmp/mdfix_run.log
```
Compare stdout against §3.1 and §3.2 of `money_delta_note_2026-07-29.md`.

| Quantity | Expectation | Tolerance |
|---|---|---|
| Case 3 (integrated) ⟨A_cos2φ⟩ for all 3 configs | Within a few % of pre-fix values | ≤5% |
| Case 3 S/N ratio vs pre-detector | Still 92–93% | ≤2 percentage points |
| Survival % | Unchanged (mask does not affect the MC loop) | Exact |
| Bin-migration (raw and weighted) | Unchanged | Exact |
| Case 1 for LOW | Nearly unchanged (peak bin already in `proj.accepted`) | ≤5% shift in ⟨A⟩ |
| Case 1 for MID / TOP | May shift; document the new values in the note | N/A |
| Case 2 (Q² slice) for all 3 configs | Peak-slice iq2 may or may not change; document | N/A |
| `Reco-accepted bins` vs `True-accepted bins` | reco ≤ true, gap only from MIN_EVENTS floor | Invariant |

Any excursion beyond tolerance for the "Unchanged" or "Invariant" rows blocks Phase D.

### V4 — Positivity checks (already present in script)
- All three `Positivity check: OK` lines still print `OK` for all cases.
- No `WARN: |y_model| >= 1` warnings in `build_phi_plot`.

---

## Documentation deliverable

`fastsim/notes/money_delta_note_2026-07-29_fix.md` — created by `doc-writer`, structured as follows:

1. **Header** — date `2026-07-29 (fix)`, one-line summary: "Reco-selection consistency fix for `money_delta_20260729.py`; introduces a canonical reconstructed-analysis mask and removes silent uses of the true-level `proj.accepted` mask on reco arrays. Adds hard-fail behavior for empty reco masks and static-analysis regression checks."
2. **§1 What triggered the fix** — code review found four issues; list them with file:line refs from the pre-fix version:
   - Reco selection reapplied truth-level `proj.accepted` to reco arrays (four sites: `find_peak_bin_reco`, `find_peak_q2_slice_reco`, `build_perbin_heatmap`, `print_smearing_diagnostics`).
   - Integrated / Q²-slice averages used only `n_reco > 0` without reapplying any acceptance.
   - No single definition of a "reconstructed-analysis mask" existed; each helper made an ad-hoc choice.
   - Diagnostics conflated truth-accepted-bin counts with reco-accepted-bin counts.
3. **§2 Definition of the reconstructed-analysis mask** — copy the `reco_analysis_mask` helper docstring. State the acceptance decision explicitly: the mask reuses the *full* `proj.accepted` (DIS cuts + η + E' acceptance) and adds `n_reco >= MIN_EVENTS`. Explain why full acceptance rather than DIS-only (4-point rationale from the plan). State the `reco_mask ⊆ proj.accepted` invariant.
4. **§3 Which helpers changed** — six-row table matching the "Approach → 2." table above, with columns: helper name / what changed / effect. Include the deliberate non-change of `sig2_per_fb_at_sumrule` as a footnote.
5. **§4 Failure mode for empty reco mask** — describe the hard-fail policy, the diagnostic message format, and the rationale (success invariant = all 16 PNGs written).
6. **§5 Static-analysis checks** — reproduce the six S1–S6 commands as code fences with expected `PASS` output. Note that these are meant to be re-runnable by any future editor as a regression guard.
7. **§6 Numerical impact** — insert the post-fix numbers from V3 (Phase C): reco-accepted bins per config, updated Case 1/2/3 ⟨A_cos2φ⟩ and S/N, delta vs pre-fix numbers from `money_delta_note_2026-07-29.md` §3.
8. **§7 Sanity checks passed post-fix** — list the checks from V1–V4 that passed with concrete numbers.
9. **§8 What did NOT change** — physics model, smearing kernel, ε_eID curve, MC seed, output filenames, plot layout, runtime, tier-B roadmap.
10. **§9 Follow-ups** — flag any Open Questions (1–3) still open at fix time.
11. **File inventory** — one bullet each for the fixed script, this plan, the new note.

Style: consistent with the existing `money_delta_note_2026-07-29.md` (Markdown, `>` callouts for important framing, code fences for formulas, tables for numeric deltas). No emojis. Cross-reference prior notes in the series.

---

## Todo

### Phase A — Planning & sign-off (this file)
- [x] Enumerate the four code-review findings.
- [x] Define `reco_analysis_mask` semantics and signature.
- [x] Correct the acceptance decision: reco mask = `proj.accepted & (n_reco >= MIN_EVENTS)` — reuses full acceptance including η and E'.
- [x] Remove nonexistent `scenario.x_max` reference from the mask definition.
- [x] Replace silent-skip empty-mask behavior with hard-fail `RuntimeError`.
- [x] Add static-analysis checks S1–S6 that must pass at end of Phase B.
- [x] Table each helper that must change and its new signature.
- [x] Specify main-loop threading of `reco_mask`.
- [x] List validation and regression checks V1–V4.
- [x] List note structure and section headings (now 11 sections).
- [ ] User review and sign-off (annotation cycle 2 → N as needed).

### Phase B — Script fix (delegate to `script-implementer` / `code-editor`)
- [ ] Add `reco_analysis_mask(proj, n_reco, min_events=MIN_EVENTS)` right after the `MIN_EVENTS` constant.
- [ ] Change signature and body of `find_peak_bin_reco(n_reco, reco_mask)`; raise `RuntimeError` if no cell passes (guarded upstream by main-loop empty-mask check).
- [ ] Change signature and body of `find_peak_q2_slice_reco(n_reco, reco_mask)`.
- [ ] Change signature and body of `compute_A_cos2phi_q2slice_reco(n_reco, A_reco, reco_mask, iq2)`.
- [ ] Change signature and body of `compute_A_cos2phi_integrated_reco(n_reco, A_reco, reco_mask)`.
- [ ] Update `build_perbin_heatmap` to accept `reco_mask` and replace `use = accepted & (n_events >= MIN_EVENTS)` with `use = reco_mask`.
- [ ] Update `print_smearing_diagnostics` to accept `reco_mask` and `true_accepted_count: int` as explicit arguments; forward `reco_mask` to `find_peak_bin_reco`; add the "Reco-accepted bins: X (vs true-accepted: Y)" print line using X = `int(reco_mask.sum())`, Y = the passed `true_accepted_count`. Do not access `proj.accepted` inside the helper (S1).
- [ ] In `main()`, right after each `smear_config(...)` call, compute `reco_mask = reco_analysis_mask(proj, n_reco, min_events=MIN_EVENTS)`, count `n_reco_acc_bins`, print the diagnostic line.
- [ ] Add the empty-mask hard-fail block: stderr print + `raise RuntimeError(f"Empty reco-analysis mask ...")` with the message string embedded **directly** in the `raise` (f-string literal, not via a `msg` variable) so S5's AST check can detect it. Include config tag, true-accepted count, and `n_reco_total` in the message.
- [ ] Add the three cross-consistency assertions (shape / floor / subset).
- [ ] Thread `reco_mask` through all six downstream call sites (heatmap, diagnostics, Case 1, Case 2 slice picker + averager, Case 3).
- [ ] Update the module docstring: add one paragraph explaining the reco-analysis mask and the full-acceptance choice.
- [ ] Update the six affected helper docstrings.
- [ ] Create `fastsim/scripts/_check_reco_mask_invariants.py` implementing S1–S6 per the checker specification and skeleton in the plan.
- [ ] Run `python3 fastsim/scripts/_check_reco_mask_invariants.py fastsim/scripts/money_delta_20260729.py`; capture output to `/tmp/mdfix_static.log`. Must exit 0 with final line `S1-S6 PASS: ...`.
- [ ] Run V1 smoke test (`--n-mc 200`); confirm 16 PNGs written and no exceptions.

### Phase C — Numerical regression (delegate to `script-implementer`, reports back to planner and doc-writer)
- [ ] Run V3 full regression at `--n-mc 1000`; capture stdout to `/tmp/mdfix_run.log`.
- [ ] Compare Case 3 ⟨A_cos2φ⟩ against pre-fix values from `money_delta_note_2026-07-29.md` §3.3; confirm within ≤5%.
- [ ] Confirm reco-accepted ≤ true-accepted bin count for all three configs (invariant).
- [ ] Confirm survival % and bin-migration numbers unchanged (exact).
- [ ] Record post-fix Case 1 / Case 2 / Case 3 numbers per config in a table for the note.
- [ ] Confirm no `WARN` on positivity checks.

### Phase D — Note (delegate to `doc-writer`)
- [ ] Create `fastsim/notes/money_delta_note_2026-07-29_fix.md` with the 11-section structure above.
- [ ] Reference `fastsim/scripts/_check_reco_mask_invariants.py` in §5 of the note (one code fence showing the invocation and expected `S1-S6 PASS` output; point at the plan for the S1–S6 specification rather than duplicating it).
- [ ] Insert the numerical-impact table from Phase C (§6 of the note).
- [ ] State the acceptance decision explicitly in §2 (full `proj.accepted`, not DIS-only), with the 4-point rationale.
- [ ] State the hard-fail empty-mask policy in §4 with the diagnostic-message format.
- [ ] Cross-reference `money_delta_note_2026-07-29.md` for physics context; do not re-derive.
- [ ] Match style/formatting of the existing 2026-07-29 note (Markdown, `>` callouts, tables). No emojis.
- [ ] Add a `See also:` bullet in the new note pointing to `money_delta_note_2026-07-29.md`.

### Phase E — Sync back to this plan
- [ ] After Phase B is merged, mark all Phase B checklist items complete here and paste the S1–S6 pass/fail summary into a new "Static-analysis results" subsection under Validation.
- [ ] After Phase C is merged, insert the observed numeric deltas into V3's table as a follow-up subsection.
- [ ] After Phase D is merged, mark the doc checklist complete and record the note path.
- [ ] Close remaining Open Questions or promote them to Tier B in `money_delta_note_2026-07-29_fix.md`.
