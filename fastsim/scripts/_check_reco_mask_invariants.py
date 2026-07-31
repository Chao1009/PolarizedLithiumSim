#!/usr/bin/env python3
"""Static-analysis regression guard for the reco-selection consistency fix.

Verifies six invariants (S1-S6) on money_delta_20260729.py using Python's
AST (abstract syntax tree) rather than grep or awk.  AST-based checking is
reliable because:
  - It cannot be confused by mentions inside docstrings or comments.
  - Legitimate occurrences (the mask helper itself; the true-level helpers;
    the main-loop line that builds the mask) are allow-listed by function
    name, not by text pattern.
  - It distinguishes real attribute access from string literals.

Exit code 0 on all-pass, 1 on any failure.  Stdout lists each check with
PASS/FAIL and, for failures, the offending line numbers.

Usage:
    python3 fastsim/scripts/_check_reco_mask_invariants.py \\
        fastsim/scripts/money_delta_20260729.py

See plan-money-delta-20260729-fix.md § "Static-analysis validation" for the
authoritative S1-S6 specification.
"""
import ast
import sys

# ── Name sets ────────────────────────────────────────────────────────────────

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

# Names treated as the old bare-count pattern in S2
BARE_NAMES = {"n_reco", "n_slice", "n_bin_map"}


# ── AST utilities ─────────────────────────────────────────────────────────────

def body_without_docstring(fn):
    """Return fn.body with the leading docstring node stripped (if present)."""
    body = fn.body
    if (body
            and isinstance(body[0], ast.Expr)
            and isinstance(body[0].value, ast.Constant)
            and isinstance(body[0].value.value, str)):
        return body[1:]
    return body


def walk_nodes(nodes):
    """Yield every AST node reachable from the given list of statement nodes."""
    for n in nodes:
        yield from ast.walk(n)


# ── S5 helper — extract literal string parts from a Constant or JoinedStr ────

def _literal_parts(node):
    """Yield the literal string parts of a Constant(str) or JoinedStr node."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        yield node.value
    elif isinstance(node, ast.JoinedStr):
        for v in node.values:
            if isinstance(v, ast.Constant) and isinstance(v.value, str):
                yield v.value


# ── Main checker ──────────────────────────────────────────────────────────────

def check(path):
    """Parse `path` and run S1-S6.  Return 0 on all-pass, 1 on any failure."""
    with open(path) as fh:
        src = fh.read()
    tree = ast.parse(src, filename=path)

    # Index top-level FunctionDef nodes by name
    funcs = {
        n.name: n
        for n in tree.body
        if isinstance(n, ast.FunctionDef)
    }

    fails = []

    # ── S1 — no reco helper reads `.accepted` ────────────────────────────────
    # For each function in RECO_HELPERS, walk the body (excluding docstring)
    # and check for any ast.Attribute node with attr == "accepted".
    for name in sorted(RECO_HELPERS):
        fn = funcs.get(name)
        if fn is None:
            fails.append(f"S1 FAIL: reco helper {name!r} not found in file")
            continue
        for node in walk_nodes(body_without_docstring(fn)):
            if isinstance(node, ast.Attribute) and node.attr == "accepted":
                fails.append(
                    f"S1 FAIL: {name} line {node.lineno}: "
                    f"attribute '.accepted' accessed (must not use proj.accepted "
                    f"inside a reco helper — pass reco_mask instead)"
                )

    # ── S2 — no reco helper uses `n_reco > 0` / `n_slice > 0` as routing ────
    # These are the old bug patterns; reco_mask already encodes the floor.
    for name in sorted(RECO_HELPERS):
        fn = funcs.get(name)
        if fn is None:
            continue
        for node in walk_nodes(body_without_docstring(fn)):
            if (isinstance(node, ast.Compare)
                    and len(node.ops) == 1
                    and isinstance(node.ops[0], ast.Gt)
                    and len(node.comparators) == 1):
                cmp = node.comparators[0]
                if isinstance(cmp, ast.Constant) and cmp.value in (0, 0.0):
                    left = node.left
                    if isinstance(left, ast.Name) and left.id in BARE_NAMES:
                        fails.append(
                            f"S2 FAIL: {name} line {node.lineno}: "
                            f"bare '{left.id} > {cmp.value}' comparison "
                            f"(old bug pattern — use reco_mask instead)"
                        )
                    elif (isinstance(left, ast.Subscript)
                          and isinstance(left.value, ast.Name)
                          and left.value.id in BARE_NAMES):
                        fails.append(
                            f"S2 FAIL: {name} line {node.lineno}: "
                            f"bare '{left.value.id}[...] > {cmp.value}' comparison "
                            f"(old bug pattern — use reco_mask instead)"
                        )

    # ── S3 — every reco helper has `reco_mask` in its parameter list ─────────
    for name in sorted(RECO_HELPERS):
        fn = funcs.get(name)
        if fn is None:
            continue
        arg_names = [a.arg for a in fn.args.args + fn.args.kwonlyargs]
        if "reco_mask" not in arg_names:
            fails.append(
                f"S3 FAIL: {name} missing 'reco_mask' parameter "
                f"(found: {arg_names})"
            )

    # ── S4 — `reco_analysis_mask` defined once, called exactly once in main() ─
    defs = [
        n for n in tree.body
        if isinstance(n, ast.FunctionDef) and n.name == MASK_BUILDER
    ]
    main_fn = funcs.get("main")
    n_calls_in_main = 0
    if main_fn is not None:
        for node in ast.walk(main_fn):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Name)
                    and node.func.id == MASK_BUILDER):
                n_calls_in_main += 1
    if len(defs) != 1 or n_calls_in_main != 1:
        fails.append(
            f"S4 FAIL: {MASK_BUILDER!r} defs={len(defs)}, "
            f"calls in main()={n_calls_in_main} (expected 1/1)"
        )

    # ── S5 — the empty-mask hard-fail RuntimeError exists in main() ──────────
    # The message string must appear DIRECTLY on the Raise node as a
    # Constant(str) or JoinedStr — not via an intermediate variable — so this
    # AST scan can see it.  See plan §"Failure mode for empty reco mask".
    ok_s5 = False
    if main_fn is not None:
        for node in ast.walk(main_fn):
            if (isinstance(node, ast.Raise)
                    and isinstance(node.exc, ast.Call)
                    and isinstance(node.exc.func, ast.Name)
                    and node.exc.func.id == "RuntimeError"):
                for arg in node.exc.args:
                    joined = "".join(_literal_parts(arg))
                    if "Empty reco-analysis mask" in joined:
                        ok_s5 = True
    if not ok_s5:
        fails.append(
            "S5 FAIL: no RuntimeError with literal 'Empty reco-analysis mask' "
            "in main() — the string must be embedded directly on the Raise "
            "(f-string literal), not via an intermediate variable"
        )

    # ── S6 — true-level helpers still reference `.accepted` (control) ────────
    # These helpers legitimately use proj.accepted; if they no longer do, the
    # fix over-reached into the true-level path.
    for name in sorted(TRUE_LEVEL_HELPERS_KEEPING_PROJ_ACCEPTED):
        fn = funcs.get(name)
        if fn is None:
            # Helper may be absent in a future refactor; not a regression here
            continue
        found = any(
            isinstance(n, ast.Attribute) and n.attr == "accepted"
            for n in walk_nodes(body_without_docstring(fn))
        )
        if not found:
            fails.append(
                f"S6 FAIL: {name} no longer references '.accepted' — "
                f"the fix may have over-reached into the true-level path"
            )

    # ── Report ────────────────────────────────────────────────────────────────
    for f in fails:
        print(f)
    if fails:
        print(f"\n{len(fails)} check(s) failed.")
        return 1
    print("S1-S6 PASS: reco-mask invariants hold.")
    return 0


if __name__ == "__main__":
    target = (
        sys.argv[1]
        if len(sys.argv) > 1
        else "fastsim/scripts/money_delta_20260729.py"
    )
    sys.exit(check(target))
