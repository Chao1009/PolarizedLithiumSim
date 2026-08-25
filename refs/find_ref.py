#!/usr/bin/env python3
"""Look up the program's reference dictionary (refs/refs_dict.json).

    python refs/find_ref.py slot             # keyword search (case-insensitive)
    python refs/find_ref.py "Eq. (9)" a2     # several terms: all must match
    python refs/find_ref.py --list           # one line per entry
    python refs/find_ref.py --key mantysaari # print the full entry

Matches are searched in every string field of every entry (identifiers,
title, keywords, key_content, used_in).  The PDFs themselves live in
refs/ (git-ignored); the dictionary records what was taken from each,
with equation/figure/slide pointers, so the source can be reopened at
the right place.
"""

import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
DICT = HERE / "refs_dict.json"


def load():
    with open(DICT, encoding="utf-8") as f:
        return json.load(f)["entries"]


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
