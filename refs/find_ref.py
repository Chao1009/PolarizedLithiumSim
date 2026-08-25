#!/usr/bin/env python3
"""Look up the program's reference dictionary (refs/refs_dict.json).

    python refs/find_ref.py slot             # keyword search (case-insensitive)
    python refs/find_ref.py "Eq. (9)" a2     # several terms: all must match
    python refs/find_ref.py --list           # one line per entry
    python refs/find_ref.py --key mantysaari # print the full entry
    python refs/find_ref.py --fetch          # download the missing arXiv PDFs
    python refs/find_ref.py --check          # which entries lack a local file

Matches are searched in every string field of every entry (identifiers,
title, keywords, key_content, used_in).  The PDFs themselves are committed
in refs/ (the EIC Yellow Report as four parts, see refs/README.md); the
dictionary records what was taken from each, with equation/figure/slide
pointers, so the source can be reopened at the right place.  `--fetch`
downloads https://arxiv.org/pdf/<id> for every entry whose `arxiv` field
is set (several ids may be separated by ';') into
refs/<id with '/' -> '_'>.pdf, skipping entries whose recorded `file`
(';'-separated list allowed) or id-named file exists; it uses urllib and
falls back to curl.  `--check` reports the local-copy status per entry.
"""

import json
import pathlib
import subprocess
import sys
import urllib.request

HERE = pathlib.Path(__file__).resolve().parent
DICT = HERE / "refs_dict.json"


def load():
    with open(DICT, encoding="utf-8") as f:
        return json.load(f)["entries"]


def arxiv_ids(entry):
    ids = entry.get("arxiv")
    if not ids:
        return []
    return [i.strip() for i in str(ids).split(";") if i.strip()]


def local_name(arxiv_id):
    return "%s.pdf" % arxiv_id.replace("/", "_")


def fetch(arxiv_id, dest):
    url = "https://arxiv.org/pdf/%s" % arxiv_id
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "polli-refs/1.0"})
        with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as f:
            f.write(r.read())
    except Exception as exc:  # sandboxed python: try curl
        print("  urllib failed (%s); trying curl" % exc)
        subprocess.run(["curl", "-sL", "--max-time", "300", "-o", str(dest),
                        url], check=True)
    ok = dest.is_file() and dest.stat().st_size > 10000 \
        and dest.read_bytes()[:5] == b"%PDF-"
    if not ok:
        dest.unlink(missing_ok=True)
        raise RuntimeError("download of %s did not produce a PDF" % url)


def recorded_files(entry):
    names = entry.get("file")
    if not names:
        return []
    return [HERE.parent / n.strip() for n in str(names).split(";") if n.strip()]


def fetch_missing(entries):
    for key, entry in entries.items():
        if any(p.is_file() for p in recorded_files(entry)):
            continue               # a local copy under the recorded name(s)
        for aid in arxiv_ids(entry):
            dest = HERE / local_name(aid)
            if dest.is_file():
                continue
            print("fetching %s -> %s (%s)" % (aid, dest.name, key))
            try:
                fetch(aid, dest)
                print("  ok, %.1f MB" % (dest.stat().st_size / 1e6))
            except Exception as exc:
                print("  FAILED:", exc)


def check(entries):
    for key, entry in entries.items():
        named = entry.get("file")
        rec = recorded_files(entry)
        have = bool(rec) and all(p.is_file() for p in rec)
        alt = [HERE / local_name(a) for a in arxiv_ids(entry)]
        have_alt = any(p.is_file() for p in alt)
        status = ("local: %s" % named if have else
                  "local (by arXiv id): %s" % ", ".join(p.name for p in alt if p.is_file())
                  if have_alt else
                  "MISSING (arXiv %s)" % entry.get("arxiv") if arxiv_ids(entry) else
                  "no free copy (%s)" % entry.get("journal"))
        print("%-52s %s" % (key, status))


def flatten(value):
    if isinstance(value, dict):
        for v in value.values():
            yield from flatten(v)
    elif isinstance(value, (list, tuple)):
        for v in value:
            yield from flatten(v)
    elif value is not None:
        yield str(value)


def show(key, entry, full=False):
    head = "%s\n  %s (%s)\n  file: %s | arXiv: %s | %s" % (
        key, entry["title"], entry.get("year", "?"), entry.get("file"),
        entry.get("arxiv"), entry.get("journal"))
    print(head)
    if full:
        print("  authors:", "; ".join(entry.get("authors", [])))
        print("  keywords:", ", ".join(entry.get("keywords", [])))
        for line in entry.get("key_content", []):
            print("  -", line)
        if entry.get("repo_usage_note"):
            print("  note:", entry["repo_usage_note"])
        for u in entry.get("used_in", []):
            print("  used in:", u)
    print()


def main(argv):
    entries = load()
    if not argv or argv == ["--list"]:
        for key, entry in entries.items():
            show(key, entry)
        return 0
    if argv[0] == "--fetch":
        fetch_missing(entries)
        return 0
    if argv[0] == "--check":
        check(entries)
        return 0
    if argv[0] == "--key":
        needle = argv[1].lower()
        for key, entry in entries.items():
            if needle in key.lower():
                show(key, entry, full=True)
        return 0
    terms = [a.lower() for a in argv]
    hits = 0
    for key, entry in entries.items():
        text = "\n".join(flatten({key: key, **entry})).lower()
        if all(t in text for t in terms):
            hits += 1
            show(key, entry)
            for line in entry.get("key_content", []):
                if any(t in line.lower() for t in terms):
                    print("  >", line)
            print()
    if not hits:
        print("no entry matches", terms)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
