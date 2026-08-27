# Plan — Internally solve signed `A_bag`, but consume `abs(A_bag)` in the observable pipeline (`money_delta_20260729.py`)

> **Dated record.** Numbers here predate the 2026-08-27 corrections: ⁶Li
> beam energies are γ-matched (40.8 / 99.5 / 137.5 GeV/u), not
> rigidity-scaled (20.5 / 50 / 137.5), and the far-forward divergence is
> per-configuration rather than a single 72.7 μrad. See plans/10. Kept as
> the record of its date, not as a current reference.


**Plan file**: `fastsim/scripts/plan-money-delta-20260729-abag-solve.md`
**Target script (operative, unchanged status)**: `fastsim/scripts/money_delta_20260729.py`
**Solver source of truth**: `fastsim/scripts/money_delta_20260724.py`, function `solve_A_from_sum_rule`
**Physics scope (locked by user)**: Interpretation A · `mid_x` · bag constraint (`c_bag = −0.012`) · per-config solve.
**Sign convention (locked by user, iteration 2)**: solve `A_bag` **signed** internally (negative under Interpretation A with `c_bag < 0`); print the signed value in the banner and pre-solve diagnostics; but consume **`abs(A_bag)`** in every downstream call that feeds the A_cos2phi / phi-modulation workflow. This preserves current visual and sign behavior on every plot and stdout row while making the internal derivation physically honest and consistent.
**Annotation cycle**: iteration 2 (user chose the "solve signed / consume abs" hybrid; the earlier "propagate signed everywhere" variant is retired — see §4.3 and §12 for what changed).
**Research inputs**: no dedicated `docs/research-*.md` for this change; sourced from the docstrings of `money_delta_20260724.py` (§"Sign convention"), `money_delta_20260728.py`, `money_delta_20260729.py`, `fastsim/notes/money_delta_note_2026-07-24.md` §3, `money_delta_note_2026-07-27.md` §3, `money_delta_note_2026-07-29.md` §"Hardcoded A_bag values", and the prior planner file `plan-money-delta-medium-priority-cleanup.md` (which deliberately deferred this change).

---

## 1. Goal

Replace the hardcoded absolute `A_BAG` magnitudes in `money_delta_20260729.py`

```python
A_BAG = {"low": 0.318, "mid": 0.310, "top": 0.297}   # currently
```

with a **per-config, signed** value computed at runtime by the sum-rule solver already implemented in `money_delta_20260724.py`. Under Interpretation A with `c_bag = −0.012`, the solver yields **negative** `A_bag` for all three configurations (LOW / MID / TOP); the magnitudes must reproduce the current hardcoded values to within numerical rounding.

**Sign-flow contract (locked)**:

1. **Solve**: `A_bag_signed = solve_A_from_sum_rule(...)` — negative, per-config, cached in `A_BAG_SIGNED`.
2. **Print** (banner + pre-solve loop stdout): both `A_bag_signed` (with sign) and `|A_bag_signed|` (magnitude, matches historical hardcoded value).
3. **Use** (every call site feeding the observable pipeline — `smear_config`, `compute_A_cos2phi_at_bin`, `build_perbin_heatmap`, `build_phi_plot`, `sig2_per_fb_at_sumrule`, and the summary-table `|A_bag|` column): **`abs_a_bag = abs(A_bag_signed)`**.

This is the "solve signed / consume abs" convention. The internal derivation is honest (the physics knows Δ ∝ c_bag < 0), the printed banner tells the user which sign fell out of the sum rule, but every downstream numerical and visual output is byte-identical to the current hardcoded-magnitude behavior. Existing PNGs, existing summary-table signs, existing stdout `A_cos2φ` rows, and existing positivity checks are all preserved.

## 2. Requirements Summary

**Confirmed (user-supplied, iteration 2)**
- `money_delta_20260729.py` remains the operative / latest script; no rename, no new script.
- Physics assumptions frozen: Interpretation A · `mid_x` (α = 0.7, β = 3) · bag constraint · per-config solve.
- Signed convention must originate from the solver rather than from a hand-typed sign flip.
- **The solver returns a signed value (negative). The observable pipeline (A_cos2phi, phi-modulation plots, heatmaps, summary table) consumes `abs(A_bag_signed)`, so existing visual and sign behavior is preserved byte-for-byte.**
- The signed value is surfaced only in the pre-solve banner / diagnostic stdout, so the user can see what sign fell out of the sum rule.
- Scope is bounded to this change and its directly associated documentation (module docstring, argparse help, banner text). Do not extend into detector-model, reco-mask, or other refactors, and do not edit `fastsim/notes/*` (that belongs to `doc-writer` — see §11 F3).

**Inferred (marked explicitly; user may override in annotation)**
- I1. `A_bag_signed` (negative) is the solver-native runtime value. `abs_a_bag = abs(A_bag_signed)` is computed once per config at the top of the main loop and is the value passed to every downstream call. Consequently every stdout row, every plot, and every summary-table cell that today prints a positive number continues to print the identical positive number.
- I2. `build_phi_plot` receives `A_bag_config = abs_a_bag`. Its three-scale ladder `s ∈ {0.1, 1.0, 10.0} × abs_a_bag` and the ratio `y_model = PZZ · (s / abs_a_bag) · A_ref · cos(2φ)` are numerically identical to today's plots. Because `A_ref` in this convention is also computed from `abs_a_bag` (via `smear_config`, see I3), no sign leaks into the plotted curves.
- I3. `smear_config(..., A_bag_config=abs_a_bag, ...)` and `compute_A_cos2phi_at_bin(..., A_bag_config=abs_a_bag, ...)` are called with the **absolute** value. Since `delta_shape_with_alphas(..., scale=abs_a_bag, ...)` uses the positive scale, `A_true_2d`, `A_reco`, `A_peakbin`, `A_q2slice`, `A_integrated`, and every per-bin `A_bin` retain today's sign (positive where they were positive before).
- I4. `A_BAG` (dict of hardcoded floats) is removed. A new module-level cache dict `A_BAG_SIGNED: dict[str, float]` is populated once at the top of `main()` after `base` and `nf2_obj` are ready; the cache is keyed by `config_tag ∈ {"low","mid","top"}`. It stores the **signed** value; the abs is recomputed at each call site with `abs(A_BAG_SIGNED[config_tag])` (cheap; makes the sign-flow contract visible at every consumer).
- I5. The startup banner and the pre-solve loop are the **only** places the signed value is visible to the user. Both print the signed value alongside the magnitude and a comparison to the previous hardcoded reference.
- I6. `sig2_per_fb_at_sumrule` accepts `scale` as an argument; if any call site passes it, it receives `abs_a_bag` (matching I3). No `20260729.py` main-path caller passes it today; the docstring is updated for consistency but no call-site change is needed.
- I7. No changes to `polli_fastsim` library files, notes, or sibling scripts (`20260728.py`, `20260724.py`, `20260725.py`).

**Out of scope (explicitly)**
- Any detector-model change (tracking resolutions, ε_eID anchor points, ECal decision).
- Any reco-mask / `MIN_EVENTS` / `_check_reco_mask_invariants.py` change.
- Any change to `20260724.py` / `20260728.py` / `20260725.py`.
- Any change to `polli_fastsim/*.py` library modules.
- Any change to the notes under `fastsim/notes/`.
- Any change to `A_lat` or lattice-model handling (bag only, per user).
- Any change to shape variant (`mid_x` only, per user).

## 3. Reuse-vs-duplication decision (explicit)

The sum-rule solver lives in `money_delta_20260724.py` as `solve_A_from_sum_rule(cfg, nf2_obj, base, variant, c_coef)`. There are three plausible integration paths; I recommend **path C** below and reject **A** and **B**.

| Path | Description | Verdict |
|---|---|---|
| **A. Import from `money_delta_20260724`** | `from money_delta_20260724 import solve_A_from_sum_rule` in `20260729.py`. | **Rejected.** Sibling scripts in `fastsim/scripts/` are executable entry points, not a package. Cross-imports between dated entry-point scripts create a hidden coupling that makes future `money_delta_20260730.py` also depend on `20260724.py`, and breaks if `20260724.py` is ever archived or renamed. No other date-suffixed script in this family imports from another. |
| **B. Promote solver to `polli_fastsim`** | Move `solve_A_from_sum_rule` into `polli_fastsim/asymmetries.py` (or a new `polli_fastsim/sumrules.py`). | **Rejected for this plan** (out of scope). Refactoring a library API is a separate, larger workstream — it would touch `20260724.py` and `20260725.py` as well, and require a library-side test. Flag it as future work in §11. |
| **C. Copy the solver verbatim into `20260729.py`, add `C_BAG` constant** | Add a `# ── Sum-rule solver (copied verbatim from money_delta_20260724.py) ──` block, with a docstring pointer to the origin script and a comment "keep in sync". | **Recommended.** Consistent with how `NuclearF2FromGrid`, `delta_shape_with_alphas`, `r1998`, and `alpha_s` are already copied verbatim from parent scripts (see the header comments on lines 121, 401, 236, 328 of `20260729.py`). Zero cross-script import, matches existing convention. |

Duplication cost: one function (~75 lines including docstring), one constant (`C_BAG = -0.012`). Benefit: `20260729.py` remains a self-contained entry-point script, matching the pattern of every other file in the family.

## 4. Approach

### 4.1 Structural changes (single file, `money_delta_20260729.py` only)

1. **Add module-level constant** (new, near `PZZ` at line 275):
   ```python
   # Sum-rule coefficient (bag model; matches money_delta_20260724.py line 313)
   C_BAG = -0.012
   ```

2. **Add solver function** (new, copied verbatim from `money_delta_20260724.py:484–557`; place immediately after `alpha_s(...)` at ~line 397, before `delta_shape_with_alphas(...)`):
   ```python
   # ══════════════════════════════════════════════════════════════════════════════
   # Sum-rule solver: A = c / ∫₀¹ x F₁(x, <Q²>) shape(x) dx
   # (copied verbatim from money_delta_20260724.py — keep in sync if the parent
   # solver is ever revised. Interpretation A: Δ = A · α_s · F₁ · x^α (1-x)^β.)
   # ══════════════════════════════════════════════════════════════════════════════

   def solve_A_from_sum_rule(cfg, nf2_obj, base, variant, c_coef):
       """... (verbatim body from money_delta_20260724.py) ..."""
       sc   = fom.Scenario(lumi_fb_per_nucleon=1.0, pol_ion_tensor=PZZ, q2_min=2.0)
       proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_obj)
       # ... rate-weighted <Q²>, numerical ∫ x F₁ x^α (1-x)^β dx / _PEAK_VALS[variant]
       # A = c_coef / integral  →  negative for c_coef = -0.012
       return A, q2_mean
   ```
   (Full body is 75 lines; the plan does not re-paste it. The implementer copies the current text of `money_delta_20260724.py:484–557` unchanged, then confirms via unit-diff that the copy is byte-identical to the source.)

3. **Delete** the hardcoded dict (lines 315–320 of the current script):
   ```python
   # REMOVED:
   # A_BAG = {"low": 0.318, "mid": 0.310, "top": 0.297}
   ```

4. **Add** a runtime-populated cache dict (new, replaces the deleted block; docstring only — populated in `main()`):
   ```python
   # Populated once per run in main() via solve_A_from_sum_rule.
   # Values are SIGNED: for Interpretation A with C_BAG = -0.012, all entries
   # are negative. Reference magnitudes from the previous hardcoded convention:
   #   |A_bag| ≈ {"low": 0.318, "mid": 0.310, "top": 0.297}  (mid_x, verified
   #   in money_delta_20260724.py output; reproduced within rounding).
   # The observable pipeline (smear_config, compute_A_cos2phi_at_bin,
   # build_phi_plot, build_perbin_heatmap, summary table) consumes
   # abs(A_BAG_SIGNED[config_tag]) — see §4.2 for the sign-flow contract.
   A_BAG_SIGNED: dict[str, float] = {}
   ```

5. **Update `main()`** — insert a new block immediately after the "Verify nuclear PDF set is installed" try/except (~line 1655), before the config loop begins:
   ```python
   # ── Solve A_bag for each config (once; cached in A_BAG_SIGNED) ────────────
   # Solver returns the SIGNED value (negative under Interpretation A with
   # C_BAG < 0). The signed value is printed here for transparency; the main
   # config loop below consumes abs(A_BAG_SIGNED[tag]) so downstream numerics
   # and plots remain identical to the previous hardcoded-magnitude behavior.
   print("Solving A_bag from bag sum rule (c = {:+.3f}, mid_x, Interpretation A)"
         .format(C_BAG))
   print("(reference magnitudes previously hardcoded: LOW=0.318, MID=0.310, TOP=0.297)")
   print("(downstream pipeline uses |A_bag|; signed value shown here for audit)")
   REFERENCE_ABS_A_BAG = {"low": 0.318, "mid": 0.310, "top": 0.297}   # audit-only
   with r_override(r1998):
       for cfg_pre, cfg_tag_pre, _label_pre in all_configs:
           nf2_pre = NuclearF2FromGrid(cfg_pre.ion, EPPS21_SET)
           A_signed, q2_mean = solve_A_from_sum_rule(
               cfg_pre, nf2_pre, base, VARIANT, C_BAG
           )
           A_BAG_SIGNED[cfg_tag_pre] = float(A_signed)
           ref_abs = REFERENCE_ABS_A_BAG[cfg_tag_pre]
           delta_pct = 100.0 * abs(abs(A_signed) - ref_abs) / ref_abs
           print(
               f"  {cfg_tag_pre.upper()}: A_bag (signed) = {A_signed:+.4f}   "
               f"|A_bag| (consumed) = {abs(A_signed):.4f}   "
               f"<Q²> = {q2_mean:.2f} GeV²   "
               f"(vs prior hardcoded |A_bag| = {ref_abs:.3f}: Δ = {delta_pct:.2f}%)"
           )
           # Hard-fail if the solver disagrees with the historical value by more
           # than 1% — protects against silent regressions in the solver / PDFs.
           assert delta_pct < 1.0, (
               f"Solved |A_bag| for {cfg_tag_pre} disagrees with hardcoded "
               f"reference by {delta_pct:.2f}% (>1% tolerance). Investigate "
               f"before proceeding."
           )
   print()
   ```
   Note: this block moves `all_configs` construction (lines 1657–1673) to just above it, so both the pre-solve loop and the main loop use the same list. **Preferred**: wrap both the pre-solve loop and the main config loop under a single `with r_override(r1998):` block to avoid duplicated context management.

6. **Update the main config loop** — replace the current line 1698:
   ```python
   # OLD:
   a_bag_config = A_BAG[config_tag]
   # NEW: solver returns signed value; observable pipeline consumes |A_bag|.
   a_bag_signed = A_BAG_SIGNED[config_tag]
   abs_a_bag    = abs(a_bag_signed)
   a_bag_config = abs_a_bag       # backwards-compatible alias for existing call sites
   ```
   Then every existing call site that currently receives `a_bag_config` continues to receive the same positive value (via the alias), because the alias now equals `abs_a_bag`. Explicit sign-flow:
   - `smear_config(..., A_bag_config=abs_a_bag, ...)` — **absolute** value (per I3). Ensures `delta_shape_with_alphas(scale=abs_a_bag, ...)` is positive, so `A_true_2d`, `A_reco`, `A_peakbin`, `A_q2slice`, `A_integrated` retain today's signs.
   - `compute_A_cos2phi_at_bin(..., A_bag_config=abs_a_bag, ...)` — **absolute** value; used inside `print_smearing_diagnostics` and any diagnostic that prints `A_reco_pk` / `A_true_peakbin`.
   - `build_perbin_heatmap(..., A_bag_config=abs_a_bag, ...)` — **absolute** value. Heatmap suptitle `s = |A_bag| = {A_bag_config:.3f}` and all three panels (`|A_bin|`, `δA_bin`, `|δA/A|`) render pixel-identically to today.
   - `build_phi_plot(..., A_bag_config=abs_a_bag, ...)` — **absolute** value. The three-scale ladder `s ∈ {0.1, 1.0, 10.0} × abs_a_bag`, the annotation `|A_bag| = {A_bag_config:.3f}`, the annotation `⟨A_cos2φ⟩ = {A_ref:.4f}` (positive because `A_ref` came from `smear_config` with `abs_a_bag`), and the plotted `y_model` curves are all byte-identical to today.
   - `sig2_per_fb_at_sumrule(..., scale=abs_a_bag, ...)` **if called** — **absolute** value. No main-path caller passes `scale` today; documented for future use.
   - The summary-table block (lines 1993–2008) already uses `A_BAG[current_config.lower()]` to print `|A_bag|`; replace with `abs(A_BAG_SIGNED[current_config.lower()])` (equivalently `abs_a_bag` if the summary-table build lives in scope). Value printed is unchanged.

   **Why keep the `a_bag_config` alias?** Two reasons: (a) it minimises the diff footprint (every existing call site keeps its keyword argument), and (b) it makes the "we solve signed but consume abs" convention visible in one line at the top of the loop rather than scattered across every call site. Reviewers see `a_bag_config = abs_a_bag` and immediately understand what the pipeline uses.

7. **Update the startup banner** (lines 1618–1621):
   ```python
   # OLD:
   print("Hardcoded |A_bag| values:")
   for k, v in A_BAG.items():
       print(f"  {k.upper()}: {v}")
   # NEW:
   print("A_bag values (solved internally from bag sum rule, mid_x, Interp. A):")
   print("(banner values printed after solver runs below)")
   ```
   The pre-solve loop (§4.1 step 5) prints the actual values immediately after. Keep the banner order: banner header → alpha_s init → PDF check → solve loop → main config loop.

8. **Update the module docstring** (lines 47–50):
   ```python
   # OLD:
   # A_bag values (hardcoded from money_delta_20260724.py mid_x output):
   #   LOW (5 × 27.5 GeV/u): |A_bag| = 0.318
   #   MID (10 × 50 GeV/u):  |A_bag| = 0.310
   #   TOP (18 × 137.5 GeV/u): |A_bag| = 0.297
   # NEW:
   # A_bag values (solved internally at startup from bag sum rule):
   #   c_bag = -0.012, mid_x (α=0.7, β=3), Interpretation A.
   #   Reference magnitudes (previous hardcoded values, reproduced within 1%):
   #     LOW ≈ 0.318, MID ≈ 0.310, TOP ≈ 0.297. Signed values are negative
   #     (Δ ∝ c_bag < 0). See solve_A_from_sum_rule() below (copied verbatim
   #     from money_delta_20260724.py).
   ```

9. **Update the argparse help block** (lines 1575–1576):
   ```python
   # OLD:
   "A_bag hardcoded from money_delta_20260724.py mid_x output:\n"
   "  LOW=0.318, MID=0.310, TOP=0.297\n\n"
   # NEW:
   "A_bag solved internally from bag sum rule (c_bag=-0.012, mid_x, Interp. A);\n"
   "signed values are negative. Reference magnitudes (within 1% of previous\n"
   "hardcoded convention): LOW≈0.318, MID≈0.310, TOP≈0.297.\n\n"
   ```

### 4.2 Sign-flow map — what is solved, what is printed, what is used

Under the "solve signed / consume abs" convention, the sign is contained entirely between (i) the solver return value and (ii) the pre-solve banner stdout. Every other site sees the absolute value. The following table makes the contract explicit and shows that **no downstream output changes sign** relative to today.

| Site | Value passed / used | Change vs today's hardcoded-magnitude output |
|---|---|---|
| `solve_A_from_sum_rule(...)` return value | **signed** (negative) | new site; did not exist before |
| `A_BAG_SIGNED[tag]` (cache dict) | **signed** (negative) | new site; did not exist before |
| Pre-solve banner stdout (`A_bag (signed) = -0.318 …`) | prints **signed** and abs | new stdout lines; visible to user |
| `a_bag_config = abs_a_bag` at top of main loop | **absolute** (positive) | value identical to old `A_BAG[tag]` |
| `smear_config(..., A_bag_config=abs_a_bag, ...)` | **absolute** | unchanged |
| `delta_shape_with_alphas(..., scale=abs_a_bag, ...)` inside `smear_config` | **absolute** | unchanged |
| `A_true_2d`, `A_reco`, `A_peakbin`, `A_q2slice`, `A_integrated` | sign identical to today | unchanged |
| Stdout `A_cos2φ (reco) = {A_peakbin:.4e}` (line 1817) | sign identical to today | unchanged |
| Stdout `⟨A_cos2φ⟩ (reco) = {A_q2slice:.4e}` (line 1908) | sign identical to today | unchanged |
| Stdout `⟨A_cos2φ⟩ (reco) = {A_integrated:.4e}` (line 1953) | sign identical to today | unchanged |
| `print_smearing_diagnostics` `A_reco_pk / A_true_peakbin` (~line 1542) | signs identical to today | unchanged |
| Summary-table `⟨A_cos2φ⟩ = {row['A_ref']:.4e}` (line 2002) | sign identical to today | unchanged |
| Summary-table `|A_bag| = {a_bag_val:.3f}` column (line ~1996) | uses `abs(A_BAG_SIGNED[cfg])` | value unchanged |
| `build_phi_plot` annotation `⟨A_cos2φ⟩ = {A_ref:.4f}` (line 1281) | sign identical to today (A_ref computed from abs_a_bag) | unchanged |
| `build_phi_plot` annotation `|A_bag| = {A_bag_config:.3f}` (line 1283) | receives `abs_a_bag` | unchanged |
| `build_phi_plot` `y_model = PZZ · (s / A_bag_config) · A_ref · cos(2φ)` with `s ∈ {0.1, 1, 10} × abs_a_bag` | positive ladder, positive A_ref | plotted curves **byte-identical** to today |
| `build_perbin_heatmap` suptitle `s = |A_bag| = {A_bag_config:.3f}` (line 1433) | receives `abs_a_bag` | unchanged |
| `build_perbin_heatmap` panels (`|A_bin|`, `δA_bin`, `|δA/A|`) | positive inputs feed abs / sqrt panels | pixel-identical PNGs |
| `mod_at_abag = abs(PZZ * A_peakbin)` (line 1810 + analogues) | unchanged (abs already applied, but input is positive to begin with) | unchanged |
| Positivity check `mod_at_10abag < 1.0` (line 1821 etc.) | unchanged | unchanged |
| `sig2_per_fb_at_sumrule` (if any caller passes `scale`) | **absolute** per I3/I6 | unchanged (no main-path caller) |
| `_check_reco_mask_invariants.py` (S1–S6) | untouched by this plan | unchanged |

### 4.3 Why we chose "consume abs" (superseding the earlier "propagate signed everywhere" variant)

An earlier iteration of this plan proposed passing the signed value everywhere and letting the sign flow into every downstream display quantity (inverting the φ plots vertically, negating the summary-table `⟨A_cos2φ⟩` entries, and prefixing every stdout `A_cos2φ` line with `-`). The user chose instead the "solve signed / consume abs" hybrid: **the sign is derived honestly from the sum rule and shown in the pre-solve banner, but the observable pipeline continues to operate on magnitudes**. Reasons this is preferable in the current scope:

- **No PNG regressions.** All 16 PNGs render byte-identically. Any downstream consumer (paper, note, slide deck) that has already embedded the current φ-modulation and heatmap PNGs remains valid.
- **No downstream comparison drift.** Summary-table diffs and stdout diffs are limited to the 3 new pre-solve loop lines; every existing line is unchanged.
- **The physics claim is not obscured.** The negative sign is still surfaced — it appears in the banner labelled `A_bag (signed) = -0.318` — so any reader who wants to see what fell out of the sum rule finds it in the first ~20 lines of stdout.
- **Sign semantics become internally consistent, but only inside the solver → cache boundary.** The signed value is contained in `A_BAG_SIGNED[tag]`. Every consumer that needs a scale takes `abs(...)` explicitly at the call site (or via the `a_bag_config = abs_a_bag` alias). This makes the contract auditable in a single grep.

**No visual or sign-behavior regressions are expected** anywhere in the observable pipeline. Validation §5.3 (previously "expected sign flips") is retired accordingly; §5 now only asserts invariance and the correctness of the new pre-solve banner.

If a future revision wants to expose the sign further down the pipeline (e.g., a plot mode that shows the signed cos(2φ) amplitude to make the negative-Δ physics visible), that is a separate change and belongs in a new plan. See §11 F5.

## 5. Validation

Under the "solve signed / consume abs" convention, validation is dominated by **invariance** checks: nothing downstream of the pre-solve banner should change. The only newly-visible outputs are the three pre-solve banner lines (each showing the signed value and its magnitude) and the module-docstring / argparse-help wording. The earlier "expected sign flip" section of this plan (V9–V11 in iteration 1) is retired.

### 5.1 Recovered-magnitude validation (numeric)

- **V1. Solver-vs-hardcoded audit** (built into the pre-solve loop, §4.1 step 5): for each of `{low, mid, top}`, assert `abs(abs(A_signed) - ref_abs) / ref_abs < 0.01`. The 1% tolerance covers PDF-set updates, α_s-table refreshes, and any small differences in the accepted-phase-space `<Q²>`. If any config exceeds 1%, the script hard-fails with a clear message before producing any plots.
- **V2. Manual cross-run** (implementer, one-time): after making the code changes, run once and confirm the three pre-solve lines print in the form `LOW: A_bag (signed) = -0.318X   |A_bag| (consumed) = 0.318X   ...`. The `A_bag (signed)` field must carry a `-` sign; the `|A_bag| (consumed)` field must match `{0.318, 0.310, 0.297}` to at least the 3rd decimal.
- **V3. Reference comparison against `money_delta_20260724.py`** (implementer, one-time): run `python3 fastsim/scripts/money_delta_20260724.py --outdir /tmp/audit_20260724` and grep the stdout for the mid_x rows of the "A-value summary table". Confirm `A_bag` values match those from `20260729.py`'s new pre-solve block to at least 4 decimal places (both scripts share `solve_A_from_sum_rule`, `NuclearF2FromGrid`, `alpha_s`, `_VARIANTS`, `_PEAK_VALS`, and `r1998` verbatim; any residual difference would indicate a copy-paste error).

### 5.2 Numerical-invariance validation (observable pipeline)

Because the observable pipeline consumes `abs_a_bag` — and `abs(A_BAG_SIGNED[tag])` reproduces the previous hardcoded `A_BAG[tag]` to within 1% — every downstream quantity should either match the pre-change value bit-for-bit (when the change is smaller than the 1% tolerance) or drift by <1%. Two invariance regimes therefore apply:

- **V4. `n_reco` bitwise invariance**: `n_reco` arrays are set by the smearing RNG (`default_rng(seed=42)`) and by `scale = abs_a_bag`. Where `abs_a_bag` matches the old hardcoded value to floating-point precision, `n_reco_new == n_reco_old` bitwise. Where the solver output differs from the hardcoded reference by even one ULP, `n_reco` may differ but only within the level of that ULP-scale drift (draws are deterministic in the seed but sensitive to the exact `scale` value). **Expected outcome**: `n_reco_new == n_reco_old` bitwise if the solver reproduces `{0.318, 0.310, 0.297}` at ≥4 decimal places; otherwise, `abs(n_reco_new - n_reco_old) / n_reco_old < 0.01` element-wise.
- **V5. `A_reco` sign invariance**: `A_reco_new` has the same sign as `A_reco_old` in every bin (both positive where before, both negative where before). Magnitude differs by <1%.
- **V6. `mod_at_abag` invariance**: for each of the 9 (config, case) rows in the summary table, the printed `P_zz·⟨A⟩ at A_bag` value must match the pre-change value to <1%.
- **V7. Positivity-check invariance**: the `pos_ok_{1,2,3}` booleans in the summary must remain OK/WARN as they were pre-change (a <1% drift in `A_peakbin` cannot flip a positivity boolean whose threshold is 1.0).
- **V8. Heatmap invariance**: the three per-bin heatmap PNGs (`money_delta_20260729_perbin_{low,mid,top}.png`) should render pixel-identically to the pre-change version when `abs_a_bag` matches the old value at ≥4 decimals. If solver output drifts at the 4th decimal, pixel-diff up to a small tolerance (`imagemagick compare -metric AE` count <100 pixels per PNG) is acceptable and does not indicate a bug. **Expected outcome for the current PDFs (July 2026)**: exact match.
- **V9. Summary-table sign invariance**: confirm all 9 `⟨A_cos2φ⟩` rows print with the **same sign** as the pre-change baseline (positive if positive before, negative if negative before). This directly supersedes the retired "expect leading `-` signs" check from iteration 1.
- **V10. φ-plot visual invariance**: `money_delta_20260729_phimodulation_peakbin_{low,mid,top}.png` and the φ-bin-width scan PNGs must render pixel-identically (or within the <100-pixel `compare -metric AE` tolerance) to the pre-change baseline. **No vertical inversion is expected** — this directly supersedes the retired "curves invert" expectation from iteration 1.
- **V11. Stdout `A_cos2φ` sign invariance**: the 9 stdout `A_cos2φ (reco)` / `⟨A_cos2φ⟩ (reco)` lines and the 3 `print_smearing_diagnostics` lines must print with the same sign and (within <1%) the same magnitude as the pre-change baseline.
- **V12. sig² invariance**: `sig2_per_fb_at_sumrule` return value is `∝ amp²` and unaffected by sign choice; if any test invokes it, output is unchanged (within the <1% solver drift tolerance).

### 5.3 Sanity checks unrelated to sign

- **V13. Static check passes**: `python3 fastsim/scripts/_check_reco_mask_invariants.py` must still exit 0 (this change touches no mask-related code; S1–S6 should be untouched).
- **V14. Empty-mask hard-fail path unchanged**: the empty-reco-mask hard-fail block (lines 1751–1764) is not touched; assertion behavior is preserved.
- **V15. `--n-mc` default unchanged** (still 1000); runtime estimate (~5–10 min) unchanged since the solver adds at most 3 `project_rates` calls at unit luminosity (each ~1–2 s).
- **V16. New pre-solve banner lines present**: confirm the 3 solver stdout lines and the 3 leading print statements ("Solving A_bag …", "(reference magnitudes …)", "(downstream pipeline uses |A_bag| …)") appear in the expected order before any config-loop output.

## 6. Open Questions

1. **Solver-vs-hardcoded assertion tolerance** — I have chosen 1% in V1. If the user wants stricter (e.g. 0.5%) or looser (e.g. 2%), change the `assert delta_pct < 1.0` value.
2. **REFERENCE_ABS_A_BAG persistence** — the audit-only reference dict `{"low": 0.318, "mid": 0.310, "top": 0.297}` is retained inline in the pre-solve loop as a sanity guard. Alternative: delete it entirely once the solver is trusted; the assert becomes a comment-only note. Default = keep it for at least one release cycle to detect PDF/library drift.
3. **`sig2_per_fb_at_sumrule` scale argument** — it currently accepts `scale` as an argument; no call site in `20260729.py`'s main path passes it. The docstring will be updated to note that any future caller should pass `abs(A_BAG_SIGNED[tag])`, matching the observable-pipeline convention. Default = docstring update only, no call-site change.

*(The prior iteration's "Option φ-A vs φ-B" question is resolved by the user's chosen convention: neither the plot nor any downstream display shows a sign change; `abs_a_bag` is used everywhere in the pipeline. This question is retired.)*

## 7. Files touched

- **Modified**: `fastsim/scripts/money_delta_20260729.py` (single file — all edits in §4.1 steps 1–9).
- **Created**: none (no new scripts, no new plan sub-files, no new notes).
- **Deleted**: none.
- **Not touched (explicit)**: `money_delta_20260724.py`, `money_delta_20260728.py`, `money_delta_20260725.py`, `_check_reco_mask_invariants.py`, all files under `polli_fastsim/`, all files under `fastsim/notes/`, all existing plan `.md` files (`plan-money-delta-20260729-fix.md`, `plan-money-delta-medium-priority-cleanup.md`, `plan-money-delta-updates.md`).

## 8. Trade-offs

- **Duplication of `solve_A_from_sum_rule` (~75 lines)** vs cross-script import or library promotion. Chosen: verbatim copy (§3), matches every other verbatim-copy convention in `20260729.py`. Cost: two copies of the solver exist. Mitigation: docstring comment "keep in sync with `money_delta_20260724.py:484–557`" and one-time V3 validation.
- **"Solve signed / consume abs" hides the negative sign from every downstream visual.** Cost: a reader who only inspects the PNGs or the summary table will not see that the bag-sum-rule Δ is negative. Benefit: no PNG or downstream-comparison regressions; the observable pipeline is unchanged. Mitigation: the sign is explicit in the pre-solve banner (`A_bag (signed) = -0.318`), in the module docstring, in the argparse help, and in the comment on `A_BAG_SIGNED`. A future plan can add a signed-plot mode if desired (§11 F5).
- **1% assertion tolerance in V1.** Cost: a genuine PDF / α_s drift of >1% will fail the script. Benefit: catches silent regressions. Mitigation: the failure mode is loud and informative; the user can widen the tolerance or update the reference values if drift is legitimate.
- **Pre-solve loop uses `r_override(r1998)` before the main loop.** Cost: `r_override` is entered twice per run unless the single-block form is used (which is preferred). Benefit: solver receives the correct `R` parameterization matching every downstream computation. This matches `20260724.py`'s convention exactly.
- **The `a_bag_config = abs_a_bag` alias in the main loop preserves the existing keyword-argument surface** at the cost of one extra local variable. Alternative: rename every call-site keyword to `abs_A_bag_config`. Chosen: keep the alias. Rationale: minimal diff, and the alias line reads as the sign-flow contract in the code itself.

## 9. Failure modes and recovery

- **Solver diverges from hardcoded reference (>1%)**: V1 assert fires with a clear message. Investigate `parton` package version, EPPS21 grid, or `_VARIANTS['mid_x']` values. Recovery: either fix the drift or update the reference dict (which is a physics decision, not a code decision).
- **`project_rates` at unit luminosity returns zero accepted events** (impossible for the three configs at Interpretation A, but guarded in the solver body): the solver's own `raise RuntimeError(...)` fires. No new failure paths introduced.
- **Downstream PNG/stdout differs from baseline** despite `abs_a_bag` matching the hardcoded value to ≥4 decimals: indicates a copy-paste error in the consumer call sites (probably one call site missed the alias and still receives the signed value). Recovery: grep for `A_BAG_SIGNED` in the main loop and confirm every consumer receives `abs_a_bag` (or the `a_bag_config` alias).

## 10. Recommendation on what to implement now

Implement everything in §4.1 (steps 1–9) and validate per §5 (V1–V16). This is the minimum-diff change that:

(a) removes the hardcoded magnitudes;
(b) sources them from the operative solver, signed;
(c) surfaces the signed value to the user in the pre-solve banner so the "sum rule → negative Δ" derivation is auditable in stdout;
(d) preserves numerical and visual behavior on every downstream output by consuming `abs_a_bag` in every call site of the observable pipeline.

Delegate the code edit to `code-editor` (this is an edit inside an existing multi-thousand-line script, not a new standalone script, so `script-implementer` is not the right handoff).

**Do not** additionally: promote the solver into `polli_fastsim` (item B in §3), touch sibling scripts, edit `fastsim/notes/*`, change any detector-model / mask code, or add a signed-plot mode. Any of those is a separate plan.

## 11. Follow-up work (not part of this plan)

- **F1. Library promotion of `solve_A_from_sum_rule`** (§3 path B). Would deduplicate the solver across `20260724.py`, `20260725.py`, `20260729.py`. Separate plan.
- **F2. Add `C_LAT = -0.009` and a `A_LAT_SIGNED` dict** if a future variant of `20260729.py` needs lattice reach. Not needed for the current bag-only scope.
- **F3. Update `fastsim/notes/money_delta_note_2026-07-29.md`** §"Hardcoded A_bag values" section to reflect the new internal solve and the "solve signed / consume abs" convention. Belongs to `doc-writer`, not this plan.
- **F4. Historical audit** — if the historical hardcoded values (0.318, 0.310, 0.297) turn out to disagree with the current solver output at the 1–2% level (e.g., PDF update since July 2026), decide whether the notes' numeric tables should be regenerated or annotated. Not part of this plan.
- **F5. Optional signed-plot mode** — add a CLI flag (e.g. `--show-signed-amplitude`) that swaps `abs_a_bag → a_bag_signed` in `build_phi_plot` only, so a reader who wants to see the negative bag amplitude visually can regenerate PNGs on demand without disturbing the default output. Separate plan; not needed today.

## 12. What changed vs iteration 1 of this plan

Iteration 1 proposed passing the signed value into every consumer of `A_bag` and letting the negative sign flow to every downstream display quantity — inverting the φ-modulation PNGs, negating the summary-table `⟨A_cos2φ⟩` column, and prefixing every stdout `A_cos2φ` row with `-`. The user chose instead the hybrid "solve signed / consume abs" convention. Concrete deltas in this iteration:

- **§1 Goal**: added the three-step sign-flow contract (Solve / Print / Use) making the boundary between "signed" and "abs" explicit.
- **§2 Requirements Summary**: I1, I2, I3, I5, I6 rewritten to reflect that every consumer receives `abs_a_bag`. I7 (was I6) unchanged.
- **§4.1 step 4**: added a docstring line to `A_BAG_SIGNED` describing the sign-flow contract.
- **§4.1 step 5**: pre-solve loop stdout labels changed to `A_bag (signed) = ...` and `|A_bag| (consumed) = ...` and a new banner line clarifies which value is consumed downstream.
- **§4.1 step 6**: introduces `a_bag_config = abs_a_bag` alias; every call site consumes the absolute value. The prior split (some receive signed, some receive abs) is replaced with a uniform "all receive abs" convention.
- **§4.2**: sign-propagation table replaced with a sign-flow map showing that no downstream site changes sign or magnitude beyond the <1% solver-vs-hardcoded drift.
- **§4.3**: former "Option φ-A vs φ-B" choice retired; new §4.3 explains why "consume abs" was chosen and what future revisions would need to expose the sign visually.
- **§5**: retired iteration 1's §5.3 "expected sign flip" checks (V9–V11). Replaced with pure invariance checks (V9–V11 now assert sign identity, plot pixel-identity, and stdout sign identity). Renumbered V1–V16.
- **§6 Open Questions**: retired "Option φ-A vs φ-B" question; kept assertion tolerance, reference-dict persistence, `sig2_per_fb_at_sumrule` docstring.
- **§8 Trade-offs**: replaced the "signed A_ref inverts figures" trade-off with the "consume abs hides sign from visuals" trade-off, plus mitigation via banner / docstring / argparse text.
- **§9 Failure modes**: removed the "user annotation reverses φ-plot sign choice" recovery path; added a "downstream diverges despite abs_a_bag matching" diagnostic.
- **§10 Recommendation**: rewritten to state the "solve signed / consume abs" outcome (a)–(d).
- **§11 Follow-up**: added F5 (optional signed-plot CLI mode).
- **Phase 1 / Phase 2 todo**: renumbered validation steps and removed the Option φ-A annotation-line task.

Nothing outside `plan-money-delta-20260729-abag-solve.md` and the target script itself is affected by these changes.

---

## Todo

### Phase 0 — Pre-flight (implementer, before editing)

- [ ] Read the current `money_delta_20260729.py` end-to-end to confirm line numbers cited in §4.1 still match the working copy.
- [ ] Read `money_delta_20260724.py:484–557` to capture the exact text of `solve_A_from_sum_rule` for verbatim copy.
- [ ] Read `money_delta_20260724.py:312–314` to capture `C_BAG = -0.012` in context.
- [ ] Run `python3 fastsim/scripts/money_delta_20260729.py --outdir /tmp/baseline_20260729 --n-mc 100` and save baseline stdout to `/tmp/baseline_stdout_20260729.txt` and baseline PNGs for later comparison. (Using `--n-mc 100` for speed; the invariance checks in V4 assume the seed is fixed so exact reproducibility holds even at low n_mc.)

### Phase 1 — Code edits (implementer, in `money_delta_20260729.py` only)

- [ ] Step 1: Add module-level constant `C_BAG = -0.012` (near line 275, after `PZZ`).
- [ ] Step 2: Add `solve_A_from_sum_rule(...)` function verbatim from `money_delta_20260724.py:484–557` (place after `alpha_s(...)` at ~line 397, before `delta_shape_with_alphas(...)`; include the "copied verbatim — keep in sync" header comment).
- [ ] Step 3: Delete the hardcoded `A_BAG = {...}` dict (lines 315–320).
- [ ] Step 4: Add empty `A_BAG_SIGNED: dict[str, float] = {}` module-level cache with the reference-values docstring block.
- [ ] Step 5: Move `all_configs = [...]` construction above the config-loop `with r_override(r1998):` block, and insert the pre-solve loop that populates `A_BAG_SIGNED` (per §4.1 step 5, including the 1% assertion `V1`). Prefer wrapping both the pre-solve loop and the main config loop under a single `with r_override(r1998):` block.
- [ ] Step 6a: In the main config loop, replace `a_bag_config = A_BAG[config_tag]` with:
  ```python
  a_bag_signed = A_BAG_SIGNED[config_tag]
  abs_a_bag    = abs(a_bag_signed)
  a_bag_config = abs_a_bag   # observable pipeline consumes |A_bag|
  ```
- [ ] Step 6b: Verify every existing call in the main loop that passed `A_bag_config=a_bag_config` continues to work (no keyword change needed; the alias preserves the value). Confirm the following sites now effectively receive `abs_a_bag`:
  - `smear_config(..., A_bag_config=a_bag_config, ...)` — abs.
  - `compute_A_cos2phi_at_bin(..., A_bag_config=a_bag_config, ...)` — abs (used inside `print_smearing_diagnostics`).
  - `build_perbin_heatmap(..., A_bag_config=a_bag_config, ...)` — abs.
  - `build_phi_plot(..., A_bag_config=a_bag_config, ...)` — abs.
- [ ] Step 6c: Update the summary-table block (lines 1993–2008) to use `abs(A_BAG_SIGNED[current_config.lower()])` in the `|A_bag| = {a_bag_val:.3f}` line.
- [ ] Step 7: Update the startup banner "Hardcoded |A_bag| values" block (lines 1618–1621) per §4.1 step 7.
- [ ] Step 8: Update the module docstring (lines 47–50) per §4.1 step 8. Include the "solve signed / consume abs" convention statement.
- [ ] Step 9: Update the argparse help block (lines 1575–1576) per §4.1 step 9. Include the sentence "downstream pipeline consumes |A_bag|; signed value shown at startup for audit".
- [ ] Optionally update `sig2_per_fb_at_sumrule` docstring to note that any future caller should pass `abs(A_BAG_SIGNED[tag])` (matches Open Question 3).

### Phase 2 — Validation (implementer, before commit)

- [ ] V1: Confirm the pre-solve loop's built-in assertion passes for all three configs.
- [ ] V2: Verify stdout shows three lines of the form `LOW: A_bag (signed) = -0.318X   |A_bag| (consumed) = 0.318X   ...` (signed field carries `-`, consumed field matches historical value to ≥3 decimals).
- [ ] V3: Run `money_delta_20260724.py --outdir /tmp/audit_20260724` and grep the "A-value summary table" `mid_x` rows for the three configs. Confirm `A_bag` values match §4.1 step 5's output to at least 4 decimal places.
- [ ] V4: Instrument a quick numpy diff: `n_reco_new == n_reco_old` (bitwise if `abs_a_bag` matches hardcoded at ≥4 decimals; else element-wise <1%). Delete instrumentation before commit.
- [ ] V5: Confirm `A_reco_new` has the same sign as `A_reco_old` in every bin, and magnitude matches within <1%.
- [ ] V6: Diff `/tmp/baseline_stdout_20260729.txt` against the new stdout. Confirm every `P_zz·⟨A⟩ at A_bag = ...` line matches within <1% (byte-identical if solver reproduces hardcoded values at ≥4 decimals).
- [ ] V7: Confirm all `pos_ok_{1,2,3}` booleans in the summary print the same OK/WARN as baseline.
- [ ] V8: `md5sum /tmp/baseline_20260729/money_delta_20260729_perbin_{low,mid,top}.png` and compare against the new run. Must match exactly (or diverge by <100 pixels via `compare -metric AE` if solver drifts at 4th decimal).
- [ ] V9: Confirm summary-table `⟨A_cos2φ⟩` rows have the **same sign** as baseline (no sign flip expected).
- [ ] V10: Compare each φ-modulation PNG (`money_delta_20260729_phimodulation_peakbin_{low,mid,top}.png` and the φ-bin-width scan PNGs) against baseline. **No vertical inversion expected**; images should match pixel-for-pixel (or within <100-pixel AE tolerance).
- [ ] V11: Confirm the 9 stdout `A_cos2φ (reco)` / `⟨A_cos2φ⟩ (reco)` lines and the 3 `print_smearing_diagnostics` lines carry the **same signs** as baseline and match magnitude within <1%.
- [ ] V12: If any test invokes `sig2_per_fb_at_sumrule` in the repo, run it and confirm unchanged output within <1%.
- [ ] V13: Run `python3 fastsim/scripts/_check_reco_mask_invariants.py` and confirm exit 0.
- [ ] V14: Confirm the empty-mask hard-fail block still triggers on a synthetic empty-mask input (existing test path unchanged).
- [ ] V15: Confirm runtime is within ~5–10 min on a laptop at `--n-mc 1000` (baseline vs new: within seconds of each other; solver adds ~3 s total).
- [ ] V16: Confirm the 3 pre-solve banner header prints ("Solving A_bag …", "(reference magnitudes …)", "(downstream pipeline uses |A_bag| …)") appear in order before any config-loop stdout.

### Phase 3 — Handoff to `supervisor` (planner-owned close-out)

- [ ] Report to supervisor: this plan file path, the "recommendation" in §10, and the list of open questions in §6 for user review before code implementation begins.
- [ ] After the user's annotation pass (if any), update this plan file per the annotation-cycle rules (preserve inline notes, mark addressed, do not silently delete).
- [ ] After code-editor completes Phase 1 and validation Phase 2, sync the checklist state in this file to reflect actual execution status per the "Execution-Status Sync" rules.
