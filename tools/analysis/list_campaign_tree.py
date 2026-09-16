#!/usr/bin/env python3
"""Record a dated, depth-limited listing of the ePIC campaign tree on the
JLab XRootD endpoint, so a plan or report can cite the listing instead of
an unchecked claim.

This is a thin wrapper around

    singularity exec $SIF xrdfs root://dtn-eic.jlab.org ls <path>

run depth-first to a fixed depth.  It is NOT a podio/edm4eic reader: it
never opens a file, it only walks directory names.  Nothing here imports
polligen or polli_fastsim, and it writes no repository output -- the
listing goes wherever --out points.

Usage (from the analysis box, outside the container; the wrapper enters it):

    python3 tools/analysis/list_campaign_tree.py \
        --path /volatile/eic/EPIC/RECO --depth 2 \
        --out /tmp/reco_tree.txt

The written file carries the UTC date, the endpoint, the container image
and the exact command template in its header, which is what makes it
citable.  Re-running it on another day produces another dated file; the
old one is the record of what the endpoint held on ITS date.
"""

import argparse
import datetime as _dt
import os
import shlex
import subprocess
import sys

DEFAULT_SERVER = "root://dtn-eic.jlab.org"
DEFAULT_SIF = "~/Projects/eic-2026/local/lib/eic_xl-nightly.sif"


def xrdfs_ls(path, server, sif, timeout):
    """Return (entries, error) for one `xrdfs ls` of `path`.

    `entries` is the list of absolute paths printed by xrdfs (empty for an
    empty directory or for a plain file).  `error` is None on exit 0.
    """
    if sif:
        cmd = ["singularity", "exec", os.path.expanduser(sif),
               "xrdfs", server, "ls", path]
    else:  # already inside eic-shell
        cmd = ["xrdfs", server, "ls", path]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=timeout)
    except subprocess.TimeoutExpired:
        return [], "TIMEOUT after %g s" % timeout
    if proc.returncode != 0:
        return [], (proc.stderr.strip() or
                    "exit %d" % proc.returncode)
    return [ln for ln in proc.stdout.splitlines() if ln.strip()], None


# Leaf extensions seen in the campaign tree.  A dot alone is NOT a leaf
# marker: the campaign directories are themselves named `25.12.0`,
# `BeAGLE1.03.02-1.2`, `pythia6.428-1.0`.
LEAF_SUFFIXES = (".root", ".txt", ".py", ".json", ".log", ".csv", ".xml",
                 ".md", ".dat", ".hepmc3", ".hepmc", ".gz", ".sh", ".yaml",
                 ".yml", ".png", ".pdf")


def is_dirlike(entry):
    """Heuristic: an entry whose basename ends in a known file extension is
    a leaf, everything else is descended into.

    xrdfs ls does not mark directories, and descending into a 3 GB .root
    file costs a round trip for nothing.  Getting this wrong is cheap in
    one direction only (a mis-classified directory is reported as a file
    and not descended), so the suffix list is explicit rather than a bare
    "contains a dot" test.
    """
    base = entry.rsplit("/", 1)[-1].lower()
    return not base.endswith(LEAF_SUFFIXES)


def walk(root, depth, server, sif, timeout, files_per_dir):
    """Depth-first walk to `depth` levels below `root`, so the recorded rows
    read as a tree (children directly under their parent).

    Returns a list of (level, path, note) rows.  `note` is "" for a listed
    directory, "(empty)" for one that listed nothing, "ERROR: ..." when the
    ls failed, and "(not descended: depth limit)" at the frontier.
    """
    rows = []
    # stack of (level, path); pushed reversed so siblings come out in order
    stack = [(0, root)]
    while stack:
        level, path = stack.pop()
        entries, err = xrdfs_ls(path, server, sif, timeout)
        if err is not None:
            rows.append((level, path, "ERROR: %s" % err))
            continue
        if not entries:
            rows.append((level, path, "(empty)"))
            continue
        rows.append((level, path, ""))
        dirs = [e for e in entries if is_dirlike(e)]
        leaves = [e for e in entries if not is_dirlike(e)]
        for leaf in leaves[:files_per_dir]:
            rows.append((level + 1, leaf, "(file)"))
        if len(leaves) > files_per_dir:
            rows.append((level + 1,
                         "... %d more file(s) in %s"
                         % (len(leaves) - files_per_dir, path), "(file)"))
        if level + 1 < depth:
            for d in reversed(dirs):
                stack.append((level + 1, d))
        else:
            for d in dirs:
                rows.append((level + 1, d, "(not descended: depth limit)"))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--path", default="/volatile/eic/EPIC/RECO",
                    help="absolute path on the endpoint to walk")
    ap.add_argument("--depth", type=int, default=2,
                    help="levels of names below --path to record; --depth N "
                         "lists directories down to N-1 levels below the root "
                         "and marks the level-N names as not descended "
                         "(default 2)")
    ap.add_argument("--server", default=DEFAULT_SERVER)
    ap.add_argument("--sif", default=DEFAULT_SIF,
                    help="container image; pass '' if already in eic-shell")
    ap.add_argument("--timeout", type=float, default=120.0,
                    help="seconds per xrdfs ls call")
    ap.add_argument("--files-per-dir", type=int, default=3,
                    help="leaf files to record per directory (default 3)")
    ap.add_argument("--out", default="-",
                    help="text file to write ('-' for stdout)")
    args = ap.parse_args(argv)

    if args.depth < 1:
        ap.error("--depth must be >= 1")
    if args.files_per_dir < 0:
        ap.error("--files-per-dir must be >= 0")

    stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    template = " ".join(shlex.quote(t) for t in (
        (["singularity", "exec", args.sif] if args.sif else [])
        + ["xrdfs", args.server, "ls", "<path>"]))

    rows = walk(args.path, args.depth, args.server, args.sif or "",
                args.timeout, args.files_per_dir)

    head = [
        "# ePIC campaign tree listing",
        "# date        : %s" % stamp,
        "# endpoint    : %s" % args.server,
        "# root path   : %s" % args.path,
        "# depth       : %d level(s) below the root path" % args.depth,
        "# container   : %s" % (args.sif or "(none; run inside eic-shell)"),
        "# command     : %s" % template,
        "# generated by: tools/analysis/list_campaign_tree.py",
        "# directories only; no file was opened (no podio/edm4eic reader).",
        "#",
    ]
    body = []
    for level, path, note in rows:
        body.append("%s%s%s" % ("  " * level, path,
                                ("   %s" % note) if note else ""))
    text = "\n".join(head + body) + "\n"

    if args.out == "-":
        sys.stdout.write(text)
    else:
        with open(args.out, "w") as fh:
            fh.write(text)
        sys.stderr.write("wrote %s (%d rows, %s)\n"
                         % (args.out, len(rows), stamp))
    return 0


if __name__ == "__main__":
    sys.exit(main())
