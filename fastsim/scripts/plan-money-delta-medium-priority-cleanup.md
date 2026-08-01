# Plan: Medium-priority money-delta cleanup — φ-sign convention, 20260724 x-axis mismatch, shape-normalization documentation drift

**Scope.** Address three remaining medium-priority "money-delta" issues surfaced during earlier audits:

1. **φ-modulation sign convention.** `money_delta_20260728.py` and `money_delta_20260729.py` build the φ-modulation as if `A_bag` were signed, but the hardcoded `A_BAG` dict already stores absolute values. The physics is correct (sig² ∝ A²), yet the code comments, variable naming, and stdout / plot annotations equivocate between "signed" and "|A|" language.
2. **20260724 raw-scale vs peak-Δ/F₁ mismatch.** `money_delta_20260724.py` still labels Plots 1–4 with `xlabel = "scale parameter s"` and prints `x-axis is raw scale s (universal across all curves)`, while `money_delta_note_2026-07-24.md` §4.1 documents those same plots as **peak Δ/F₁**, with explicit vertical-line positions (1.6×10⁻², 1.2×10⁻²) that only make sense on a peak-Δ/F₁ axis. The reach numbers are unaffected; the label/interpretation channel is broken.
3. **Shape-normalization documentation drift.** All formulas in the notes (Eq. 1 / 1A across `money_delta_note_2026-07-24.md`, `2026-07-27.md`, `2026-07-28.md`, `2026-07-29.md`, and `money_delta_uptodate.md`) write `Δ = s · α_s · F₁ · x^α · (1-x)^β`, which omits the `1/_PEAK_VALS[variant]` factor the implementation actually applies. The docstring of `delta_shape_with_alphas` in `money_delta_20260729.py:404-432` describes the peak-normalization; the notes do not.

**Files touched by this plan**
- `fastsim/scripts/money_delta_20260728.py` — surgical code changes.
- `fastsim/scripts/money_delta_20260729.py` — surgical code changes.
- `fastsim/scripts/money_delta_20260724.py` — code changes limited to axis-label and stdout strings.
- `fastsim/notes/money_delta_note_2026-07-24.md` — text-only clarifications (annotated as addenda; no numeric table rewrites).
- `fastsim/notes/money_delta_note_2026-07-27.md` — text-only formula clarification.
- `fastsim/notes/money_delta_note_2026-07-28.md` — text-only formula clarification.
- `fastsim/notes/money_delta_note_2026-07-29.md` — text-only formula clarification.
- `fastsim/notes/money_delta_uptodate.md` — text-only formula clarification (Eq. 1A block).
- `fastsim/notes/money_delta_note_2026-07-31_cleanup.md` — **new** consolidation note documenting the three fixes and their scope.
- `fastsim/scripts/plan-money-delta-medium-priority-cleanup.md` — this file.

**Not touched**
- `polli_fastsim/` library. All fixes are script-local or documentation.
- `money_delta_note_2026-07-29_fix.md`. Its content is orthogonal (reco-selection mask); no rewrite needed.
- Prior scripts (`money_delta_20260715.py` … `money_delta_20260725.py` except `20260724.py`). They are historical snapshots; no re-labeling.
- The reach numbers in any note. All three fixes are labeling / naming / documentation and cannot change any computed value.

**Annotation cycle.** Iteration 1 (initial draft). Awaiting user review.

**Input sources**
- Direct reads of `money_delta_20260724.py`, `money_delta_20260728.py`, `money_delta_20260729.py`.
- `money_delta_note_2026-07-24.md` §§3, 4, 4.1 (documents peak-Δ/F₁ axis and A-value table).
- `money_delta_note_2026-07-27.md`, `2026-07-28.md`, `2026-07-29.md` (formulas Eq. 1 / 1A).
- `money_delta_uptodate.md` (canonical formula index, Eq. 1A blocks at lines 370, 463).
- Existing planner artifacts `plan-money-delta-updates.md`, `plan-money-delta-20260729-fix.md` for style / structure.
- Prior planner review that surfaced these three items as "medium-priority" (user-supplied).

---

## Requirements Summary

**Confirmed requirements** (from the user's request)

- Bounded scope: only the three enumerated items. No unrelated refactors, no library changes, no new physics.
- Code fixes must be separated from documentation fixes.
- Candor: for each item, explicitly state what changes now vs. what is documented as a caveat.
- Include validation and rollout/documentation steps.
- The plan lives in a planner-owned file in the repo.

**Inferred assumptions** (flag if wrong)

- The `A_BAG` dict entries (0.318, 0.310, 0.297) are meant to be `|A_bag|` values, per the module docstring on line 47–50 and the comment on line 315 of `money_delta_20260729.py`. The audit finding is a naming/comment issue, not a physics-sign question.
- The 20260724 x-axis mismatch is a label-only fix: internal `SCALES` values, `sig2_per_fb_at_sumrule` inputs, and the analytic `L_5σ = 25 / (sig² · (s/s₀)²)` extrapolation stay in the raw-`s` convention. Only the displayed axis (and the two lines of stdout that describe it) get the peak-Δ/F₁ transform, matching what the note already claims in §4.1.
- The shape-normalization drift is a note-only fix — Eq. 1 / 1A is rewritten with the explicit `1/_PEAK_VALS[variant]` factor, and one sentence explaining that the code enforces "at x = x_peak (f₁ = 1, α_s = 1), Δ = s" is added. The implementation is already correct; no code touched for this item.
- The consolidation note (`money_delta_note_2026-07-31_cleanup.md`) is preferred over inline addendum sections in each of the four dated notes. Dated notes remain historical snapshots; the new note pulls the three fixes together and is cross-referenced from `money_delta_uptodate.md`.

## Resolved defaults (locked unless the user flags)

1. Consolidation note filename `money_delta_note_2026-07-31_cleanup.md` (matches the existing "date-suffix + descriptor" pattern already established by `money_delta_note_2026-07-29_fix.md`).
2. Formula-drift fix in the historical dated notes = one-line "Corrected: divide by `_PEAK_VALS[variant]`" annotation immediately below each Eq. 1 / 1A block, plus a cross-reference to the consolidation note. Do **not** rewrite the surrounding text — historical notes stay historical.
3. 20260724 axis-label change: keep the `SCALES` internal grid unchanged; recompute displayed x-axis via a `_scale_to_peak_delta_over_f1(scale, alpha, beta, alpha_s_at_qbar)` conversion helper local to the script, using the same numbers already documented in the note (α_s(⟨Q²⟩) × peak_shape = 0.30 × 0.170 = 0.051 for mid_x @ MID). No new sig² evaluations.
4. φ-sign fix in `20260728.py` and `20260729.py`: **rename the variable, tighten the comments, remove the "signed convention" wording**. Do not introduce a signed `A_BAG_SIGNED` dict; the current absolute-value dict stays authoritative. This is the least invasive fix consistent with the physics being sign-agnostic.

## Open Questions

1. **Should the 20260724 axis-label fix also update the printed A_bag row in Table §3** (which shows raw `-0.310`, matching the raw-scale convention)? Default: **no** — that table is explicitly annotated as "raw `scale` parameters" in its footnote (§3, note after the table); it is internally consistent. Flag if the user wants Table §3 to also carry a "peak Δ/F₁" column.
2. **Should the consolidation note re-run any of the 20260728 / 20260729 scripts** to confirm that renaming does not perturb any PNG? Default: **yes** — one byte-diff check per script under `matplotlib.use('Agg')`. Flag if the user wants stricter (numeric stdout diff) or looser (skip the re-run).
3. **Should `money_delta_20260725.py` also get the axis-label fix**? It has the same `xlabel = "scale parameter s"` and the same "Note: on Plots 1-4, x-axis is raw scale s (universal across all curves)" stdout. Default: **no** — the 20260724 fix is enough for the medium-priority audit; 20260725 is out of scope unless the user asks. This is flagged explicitly here as the most likely scope-expansion request.

---

## Candor: what to actually change now vs. what to leave documented

**Change in code now (small, low-risk, no numeric impact)**
- **Item 1 (φ-sign).** Rename `A_bag_config` → `abs_A_bag_config` **only in the two functions and the dict-key access site where the misleading "signed → abs" comment currently lives**. Update comments at `money_delta_20260728.py:429-432` and `money_delta_20260729.py:946-949` to read plainly: "A_BAG stores |A|; the physics is sign-agnostic because sig² ∝ A². φ-modulation amplitudes use |A_bag| explicitly." Update the two stdout blocks that print `"|A_bag| values"` (already correct) and the two docstring "signed convention" mentions. Do **not** thread the new name through every internal helper — that would be a refactor, not a cleanup. The three call sites already labeled `|A_bag|` in the docstring headers (lines 644, 932, 1008, 1089, 1172, 1327 of 20260729) can keep their `A_bag_config` names because they are already documented as absolute values.
- **Item 2 (20260724 axis).** Add one 8-line helper `_scale_to_peak_delta_over_f1(scale, alpha, beta, alpha_s_qbar)` at module scope, insert one `disp_x = _scale_to_peak_delta_over_f1(SCALES, ...)` call before each of the four `ax.set_xlabel(r"scale parameter $s$", ...)` sites, plot against `disp_x`, and change the label to `r"peak $\Delta/F_1$ (mid_x reference, $\langle Q^2 \rangle = 7.4$ GeV²)"`. Update the two `# Vertical lines at raw |A|` comments and the vertical-line `axvline` x-positions to use the same converted values (1.6×10⁻² and 1.2×10⁻²) already in the note. Update the one stdout print at line 1366 to say `x-axis is peak Δ/F₁ (converted from raw scale s via _scale_to_peak_delta_over_f1)`.

**Document as a caveat, do not change code**
- **Item 3 (shape-normalization formula).** The implementation is correct (`delta_shape_with_alphas` uses `norm = scale / _PEAK_VALS[variant]`). The formulas in the notes are the ones that need the `1/_PEAK_VALS[variant]` factor added. Do **not** modify the implementation; do **not** rewrite the surrounding narrative. Add one corrected form of Eq. 1 / 1A per note, with a short justification line ("normalization chosen so that at x = x_peak with F₁ = 1, α_s = 1, one recovers Δ = s"). This is a pure documentation fix.
- **Item 1 (φ-sign) beyond the two annotated sites.** The `s = A_bag_config` line at 20260729:949 (and analogue in 20260728) is the only place the "signed convention" wording currently misleads a reader. Once that comment is rewritten and the local variable is renamed, propagating the rename globally is refactor territory and out of scope. **Document in the consolidation note** that `A_BAG` is understood to be `|A_bag|` throughout the family of scripts; future scripts should follow the same convention.
- **Item 2 (20260725 sibling).** Do not fix now; document in the consolidation note that `money_delta_20260725.py` has the same axis-label issue but is not in scope for this cleanup. If the user later needs 20260725 outputs regenerated, apply the same 8-line helper pattern.
- **Historical dated notes (20260727 / 20260728 / 20260729).** The formula-drift correction is a single-line annotation appended to each Eq. 1 / 1A block, plus a cross-reference to the consolidation note. Do not restructure or rewrite these notes.

**Explicitly NOT touched**
- No numeric tables, no reach numbers, no plot values, no library files.
- No changes to `A_BAG` dict entries. They remain `{low: 0.318, mid: 0.310, top: 0.297}`.
- No renames of `A_bag_config` outside the two annotated sites per script.
- No new script snapshots (no `money_delta_20260731_*.py`). All code lives in the existing 20260724 / 20260728 / 20260729 files.

---

## Detailed approach — code fixes

### Item 1 code — φ-sign in `money_delta_20260728.py` and `money_delta_20260729.py`

**Sites to change per script** (line numbers below are for 20260729; 20260728 has near-identical structure with different line numbers).

1. `money_delta_20260729.py:315` (comment).
   Before:
   ```python
   # Hardcoded A_bag values (abs) from money_delta_20260724.py mid_x output
   ```
   After:
   ```python
   # Hardcoded |A_bag| values from money_delta_20260724.py mid_x output.
   # Physics is sign-agnostic (sig² ∝ A²); all downstream code assumes |A_bag|.
   ```

2. `money_delta_20260729.py:946-949` (misleading "signed convention" comment).
   Before:
   ```python
   # Use A_bag as the scale (signed convention: negative A_bag → negative Δ;
   # we use the abs value per the plan since A_bag is pre-negated in the
   # hardcoded constants)
   scale = A_bag_config   # already absolute value from A_BAG dict
   ```
   After:
   ```python
   # A_BAG stores |A_bag|; we pass it directly as `scale`. The observable
   # sig² ∝ A² is sign-agnostic, and the φ modulation uses the absolute
   # amplitude (see build_phi_plot, which multiplies |s / A_bag_config|).
   scale = A_bag_config
   ```

3. `money_delta_20260729.py:1208-1209` (φ-modulation amplitude).
   Before:
   ```python
   # Modulation: y_model(φ, s) = P_zz · (s/A_bag) · A_ref · cos(2φ)
   amplitude = PZZ * (s / A_bag_config) * A_ref
   ```
   After:
   ```python
   # Modulation: y_model(φ, s) = P_zz · (s/|A_bag|) · A_ref · cos(2φ).
   # s and A_bag_config are both absolute values (A_BAG dict comment).
   amplitude = PZZ * (s / A_bag_config) * A_ref
   ```
   (No numerical change: the division was already using the absolute stored value; only the comment gets rewritten.)

4. `money_delta_20260729.py:1603` (stdout).
   Before:
   ```python
   print("Observable: dN/dφ = N_flat · [1 + P_zz · (s/A_bag) · <A_cos2φ> · cos(2φ)]")
   ```
   After:
   ```python
   print("Observable: dN/dφ = N_flat · [1 + P_zz · (s/|A_bag|) · <A_cos2φ> · cos(2φ)]")
   ```

5. `money_delta_20260728.py:429-432` and `money_delta_20260728.py:610-611` — apply the analogous edits at the matching lines.

Do **not** rename the `A_bag_config` parameter across the six helpers that already carry an `|A_bag|` marker in their docstring headers (lines 644, 932, 1008, 1089, 1172, 1327 of 20260729 and the analogues in 20260728). The docstrings are already correct; only the ambiguous inline comments and the stdout string need to change.

### Item 2 code — `money_delta_20260724.py` axis-label fix

**New helper** (add after the `_PEAK_VALS` block near line ≈ 250):

```python
def _scale_to_peak_delta_over_f1(scale, alpha, beta, alpha_s_qbar):
    """Convert internal raw `scale` (= s) to displayed peak Δ/F₁.

    peak(Δ/F₁) = scale · α_s(⟨Q²⟩) · max_x[x^α · (1-x)^β]

    Matches the conversion documented in money_delta_note_2026-07-24.md §4.1.
    For the mid_x reference shape (α = 0.7, β = 3) at MID (⟨Q²⟩ = 7.4 GeV²),
    the conversion factor is 0.30 × 0.170 ≈ 0.051.
    """
    xp = alpha / (alpha + beta)
    peak_shape = (xp ** alpha) * ((1.0 - xp) ** beta)
    return np.asarray(scale, dtype=float) * alpha_s_qbar * peak_shape
```

**Site edits** (four `ax.set_xlabel(r"scale parameter $s$", ...)` sites at lines 668, 785, 872, 951, plus the four vertical-line blocks at 652, 769, 858, 932, plus the stdout at 1366). Pattern:

Before:
```python
ax.plot(SCALES, curve, ...)
...
ax.axvline(abs(A_BAG), ls="--", color="k", label=r"$|A_\mathrm{bag}|$")
ax.axvline(abs(A_LAT), ls="--", color="gray", label=r"$|A_\mathrm{lat}|$")
...
ax.set_xlabel(r"scale parameter $s$", fontsize=10)
```

After:
```python
# Displayed x-axis = peak Δ/F₁, computed from raw scales via the
# mid_x reference shape at ⟨Q²⟩ = 7.4 GeV² (see note §4.1).
disp_x = _scale_to_peak_delta_over_f1(
    SCALES, alpha=0.7, beta=3.0, alpha_s_qbar=0.30
)
ax.plot(disp_x, curve, ...)
...
ax.axvline(_scale_to_peak_delta_over_f1(abs(A_BAG), 0.7, 3.0, 0.30),
           ls="--", color="k", label=r"peak $\Delta/F_1$ (bag) $\approx 1.6\times10^{-2}$")
ax.axvline(_scale_to_peak_delta_over_f1(abs(A_LAT), 0.7, 3.0, 0.30),
           ls="--", color="gray", label=r"peak $\Delta/F_1$ (lat) $\approx 1.2\times10^{-2}$")
...
ax.set_xlabel(r"peak $\Delta/F_1$  (mid$_x$ reference, $\langle Q^2\rangle = 7.4$ GeV$^2$)",
              fontsize=10)
```

Stdout line 1366:
Before: `print("Note: on Plots 1-4, x-axis is raw scale s (universal across all curves).")`
After:  `print("Note: on Plots 1-4, x-axis is peak Δ/F₁ (converted from raw s via _scale_to_peak_delta_over_f1; see note §4.1).")`

**Do not** change:
- Internal `SCALES` grid.
- `sig2_per_fb_at_sumrule` and the analytic `L_5σ` extrapolation.
- Table §3 in the note (per Open Question 1 default).
- The 9-curve legends (they identify config/shape, not axis units).

---

## Detailed approach — documentation fixes

### Item 3 doc — shape-normalization drift

For each of the four notes that carry Eq. 1 / 1A (`money_delta_note_2026-07-24.md:78`, `2026-07-27.md:76`, `2026-07-28.md:104,152`, `2026-07-29.md:92`, `money_delta_uptodate.md:206,370,463`), insert a single line **immediately below** each Eq. 1 / 1A block:

```
[Implementation note — corrected form]

    Δ(x, Q²) = (s / P_shape) · α_s(Q²) · F₁(x, Q²) · x^α · (1−x)^β

where P_shape = max_x[x^α (1−x)^β] is the shape peak (mid_x: 0.170;
low_x: 0.406; high_x: 0.088). The 1/P_shape factor is the peak-normalization
that fixes the code's convention: at x = x_peak with F₁ = 1 and α_s = 1,
Δ = s. See `delta_shape_with_alphas` in money_delta_20260729.py:404-432,
and the consolidation note money_delta_note_2026-07-31_cleanup.md §3 for
scope.
```

Do **not** rewrite the surrounding narrative. Do not touch the numeric tables. This is a minimum-diff clarification.

### Item 1 & 2 doc — φ-sign and 20260724 axis in the consolidation note

Create `fastsim/notes/money_delta_note_2026-07-31_cleanup.md` with the structure below. Delegated to `doc-writer`.

1. **Header** — date `2026-07-31 (cleanup)`, one-line summary: "Three medium-priority cleanup items on the money-delta family: φ-modulation sign convention (|A_bag| explicit), 20260724 x-axis label matches note §4.1 (peak Δ/F₁), and shape-normalization documentation drift (Eq. 1 / 1A) corrected across the notes."
2. **§1 What triggered the cleanup** — cross-reference the planner audit that surfaced the three items; list each with file:line refs.
3. **§2 φ-modulation sign convention** — quote the `A_BAG` docstring, the two amended comments, and the stdout line. State that no numeric result changes.
4. **§3 Shape-normalization documentation drift** — reproduce the corrected Eq. 1 form. Table listing `_PEAK_VALS` for the three variants (low_x 0.406, mid_x 0.170, high_x 0.088). Explicitly note: implementation was always correct; the formula in the notes was missing the 1/P_shape factor.
5. **§4 20260724 axis-label fix** — quote the new `_scale_to_peak_delta_over_f1` helper, list the four `ax.set_xlabel` sites, and confirm that reach numbers (L_5σ, δA/A, σ²) are unaffected. Cross-reference `money_delta_note_2026-07-24.md` §4.1.
6. **§5 What is out of scope** — enumerate: `money_delta_20260725.py` (same axis issue, not fixed here), all A_bag_config renames outside the two annotated sites per script, `A_BAG` dict itself, library files, historical notes' numeric tables, and the reco-mask fix (already handled in `money_delta_note_2026-07-29_fix.md`).
7. **§6 Validation summary** — reproduce V1–V4 from this plan's Validation section with pass/fail marks, once Phase C completes.
8. **§7 File inventory** — one bullet per script and note touched.
9. **See also** — bullet list of the four dated notes with the corrected-Eq. addendum, plus `money_delta_note_2026-07-29_fix.md`.

### `money_delta_uptodate.md` update

Add one bullet in the "Revision log" table (around line 314):

```
| 2026-07-31 | 20260728/20260729/20260724 scripts, various notes | money_delta_note_2026-07-31_cleanup.md | Medium-priority cleanup | φ-sign |A_bag| explicit; 20260724 x-axis matches note §4.1; Eq. 1/1A corrected across notes; no numeric changes |
```

Update the corrected Eq. 1A blocks at lines 370 and 463 with the addendum in Item 3 above.

---

## Trade-offs and considerations

- **Minimum-diff bias.** Every fix is chosen to be the smallest change that resolves the audit finding. This prevents scope creep and keeps the risk of accidental numeric perturbation at zero. The cost is that `A_bag_config` continues to appear as a variable name; readers must rely on the `|A_bag|` docstring markers plus the amended inline comment.
- **Peak-normalization convention lock-in.** By documenting the 1/P_shape factor in every note that carries Eq. 1 / 1A, we lock the convention. Future scripts that adopt a different normalization (e.g. area-normalized shapes) will need a new equation number and a new consolidation note.
- **20260724 vertical-line accuracy.** The vertical lines use `_scale_to_peak_delta_over_f1(|A_BAG|, 0.7, 3.0, 0.30)`. For the LOW and TOP configs this is a mid_x-reference approximation — same as the note §4.1 already documents (Caveat 2). The 9-curve plots inherit that caveat; the mid+mid_x-only Plots 3–4 are exact.
- **Historical notes as snapshots.** By restricting the historical dated notes to a one-line Eq. 1 / 1A addendum + cross-reference, we preserve their "what we thought on that date" character. Readers who want the consolidated current view read the consolidation note.
- **Non-issue risk: A_bag_config semantics change.** None. `A_BAG` stays absolute; the "signed convention" comment we remove was already misleading (the code never actually applied a sign flip). The one-word cleanup is a no-op physically.

---

## Validation

### V1 — code-fix smoke tests

```bash
# 20260728 (parent script; no MC loop)
python3 fastsim/scripts/money_delta_20260728.py --outdir /tmp/mdcleanup_28

# 20260729 (detector-realistic; use --n-mc 200 for speed)
python3 fastsim/scripts/money_delta_20260729.py --outdir /tmp/mdcleanup_29 --n-mc 200

# 20260724 (sum-rule reference)
python3 fastsim/scripts/money_delta_20260724.py --outdir /tmp/mdcleanup_24
```

**Pass criteria**
1. All three scripts exit 0.
2. PNG counts match pre-fix (11 for 20260728, 16 for 20260729, 9 for 20260724).
3. No new warnings on stdout, no `WARN: |y_model| >= 1`, no `AssertionError`.
4. The `20260729_fix` reco-mask AST check (`_check_reco_mask_invariants.py`) still passes: exit 0 with `S1-S6 PASS`.

### V2 — byte-identical outputs (20260728, 20260729)

Under `matplotlib.use('Agg')` and pinned matplotlib version, the φ-sign fix must not perturb any pixel of any PNG. Compare against a pre-fix baseline:

```bash
sha256sum fastsim/out/money_delta/money_delta_20260728_*.png > /tmp/mdcleanup_28_after.sha
sha256sum fastsim/out/money_delta/money_delta_20260729_*.png > /tmp/mdcleanup_29_after.sha
diff /tmp/mdcleanup_28_before.sha /tmp/mdcleanup_28_after.sha
diff /tmp/mdcleanup_29_before.sha /tmp/mdcleanup_29_after.sha
```

Pass criteria: both diffs empty. (The φ-sign fix touches only comments and one stdout line, so this must hold exactly.)

### V3 — 20260724 axis-label numeric consistency

```bash
python3 fastsim/scripts/money_delta_20260724.py --outdir /tmp/mdcleanup_24 2>&1 | tee /tmp/mdcleanup_24.log
```

Pass criteria:
1. `Plots 1-4, x-axis is peak Δ/F₁` appears in stdout.
2. Vertical lines fall at 1.6×10⁻² (bag) and 1.2×10⁻² (lattice) on all four plots, matching `money_delta_note_2026-07-24.md` §4.1 numbers.
3. Reach numbers in stdout table (§4 of the note) are byte-identical to pre-fix — the axis change is purely display.

### V4 — documentation-fix rendering check

For each of the four dated notes plus `money_delta_uptodate.md`, confirm:
1. The corrected Eq. 1 / 1A addendum block renders as a code fence in Markdown preview.
2. The addendum sits immediately after the original Eq. 1 / 1A block, before any surrounding paragraph.
3. Cross-reference to `money_delta_note_2026-07-31_cleanup.md` resolves to the new file.
4. No accidental removal of the original Eq. 1 / 1A block, table entries, or surrounding narrative.

### V5 — consolidation note completeness

`money_delta_note_2026-07-31_cleanup.md` must:
1. Cover all three items (§2, §3, §4 in its structure).
2. Include an explicit "out of scope" section (§5) naming `money_delta_20260725.py`, the `A_BAG` dict, and library files.
3. Cross-reference `money_delta_note_2026-07-29_fix.md`, `money_delta_note_2026-07-24.md`, and this plan file.
4. Be listed in the revision log of `money_delta_uptodate.md`.

---

## Rollout

1. **Phase A** (this file) — user review and sign-off of scope, split, and defaults.
2. **Phase B — code fixes** (delegated to `script-implementer`): items 1 and 2 above, `money_delta_20260728.py`, `money_delta_20260729.py`, `money_delta_20260724.py`. Runs V1 + V2 + V3 at end.
3. **Phase C — documentation fixes** (delegated to `doc-writer`): Eq. 1 / 1A addendum in the four dated notes + `money_delta_uptodate.md`, plus authoring `money_delta_note_2026-07-31_cleanup.md`. Runs V4 + V5 at end.
4. **Phase D — sync back** to this plan file: mark checklist items complete, insert observed byte-diff results (V2), and close open questions.

Code fixes (Phase B) and doc fixes (Phase C) are independent; run in parallel if convenient. The consolidation note (Phase C) references Phase B outcomes, so finalize Phase C last.

---

## Recommendation summary — what to execute now

- **Execute now**: Items 1, 2, and 3 as scoped above. All three are low-risk minimum-diff changes.
- **Item 1 (φ-sign, code)**: yes — narrow comment/stdout cleanup in `money_delta_20260728.py` and `money_delta_20260729.py`. Zero numeric impact; byte-identical PNGs expected.
- **Item 2 (20260724 axis, code)**: yes — 8-line helper + four axvline/xlabel site edits + one stdout line. Reach numbers unchanged; plots re-labeled to match the note.
- **Item 3 (shape-normalization, docs only)**: yes — one addendum block per note. No code change. Fixes the reader-facing formula so it matches the implementation.
- **Do not execute now**: the `money_delta_20260725.py` sibling of item 2; global rename of `A_bag_config`; any change to the `A_BAG` dict; any rewrite of historical notes' narrative or tables.

---

## Todo

### Phase A — planning and sign-off (this file)
- [x] Enumerate the three items and confirm scope.
- [x] Separate code fixes vs documentation fixes.
- [x] Draft candor section (now vs. caveat).
- [x] Draft file inventory and "not touched" list.
- [x] Draft validation V1–V5.
- [ ] User review and sign-off (annotation cycle 1 → N as needed).

### Phase B — code fixes (delegate to `script-implementer`)
- [ ] `money_delta_20260729.py`: rewrite comment at line 315 (`A_BAG` dict header).
- [ ] `money_delta_20260729.py`: rewrite comment block at lines 946–949 (remove "signed convention" wording).
- [ ] `money_delta_20260729.py`: rewrite comment at lines 1208–1209 (φ-modulation formula).
- [ ] `money_delta_20260729.py`: rewrite stdout at line 1603 (`"Observable: dN/dφ …"`).
- [ ] `money_delta_20260728.py`: apply the matching edits at lines 236 (dict header), 429–432 (comment), 610–611 (φ-formula), 933 (stdout).
- [ ] `money_delta_20260724.py`: add `_scale_to_peak_delta_over_f1` helper after the `_PEAK_VALS` block.
- [ ] `money_delta_20260724.py`: update four `ax.set_xlabel(r"scale parameter $s$", ...)` sites (lines 668, 785, 872, 951) to plot against `disp_x` and use the new peak-Δ/F₁ label.
- [ ] `money_delta_20260724.py`: update four vertical-line blocks (lines 652, 769, 858, 932) to place `axvline` at converted peak-Δ/F₁ positions using the same helper.
- [ ] `money_delta_20260724.py`: update stdout at line 1366 to state "x-axis is peak Δ/F₁".
- [ ] Snapshot PNG SHA256 sums before Phase B start (`/tmp/mdcleanup_28_before.sha`, `/tmp/mdcleanup_29_before.sha`, `/tmp/mdcleanup_24_before.sha`).
- [ ] Run V1 smoke tests; confirm exit 0 and expected PNG counts.
- [ ] Run V2 byte-diff for 20260728 and 20260729; must be identical.
- [ ] Run V3 for 20260724; confirm stdout says peak Δ/F₁, verticals at 1.6×10⁻² and 1.2×10⁻², reach numbers unchanged.
- [ ] Run `python3 fastsim/scripts/_check_reco_mask_invariants.py fastsim/scripts/money_delta_20260729.py`; must still `S1-S6 PASS`.

### Phase C — documentation fixes (delegate to `doc-writer`)
- [ ] `money_delta_note_2026-07-24.md`: insert Eq. 1 addendum immediately after line 78.
- [ ] `money_delta_note_2026-07-27.md`: insert Eq. 1A addendum immediately after line 76.
- [ ] `money_delta_note_2026-07-28.md`: insert Eq. 1A addendum after each of lines 104 and 152.
- [ ] `money_delta_note_2026-07-29.md`: insert Eq. 1A addendum immediately after line 92.
- [ ] `money_delta_uptodate.md`: update Eq. 1A blocks at lines 370 and 463 with the addendum.
- [ ] `money_delta_uptodate.md`: add revision-log row for 2026-07-31 (near line 314).
- [ ] Create `fastsim/notes/money_delta_note_2026-07-31_cleanup.md` with the 9-section structure specified above.
- [ ] Cross-reference `money_delta_note_2026-07-29_fix.md`, `money_delta_note_2026-07-24.md`, and this plan file.
- [ ] Run V4 (Markdown preview) and V5 (completeness).

### Phase D — sync back to this plan (planner)
- [ ] Mark all Phase B and Phase C checklist items complete.
- [ ] Paste V2 byte-diff results verbatim into a new "Validation results" subsection under §Validation.
- [ ] Paste V3 stdout excerpt confirming peak-Δ/F₁ axis into the same subsection.
- [ ] Resolve the three open questions or promote them to follow-ups in the consolidation note.
- [ ] Confirm no scope creep occurred; if `20260725.py` or a global `A_bag_config` rename is added later, it becomes a new plan file.
