"""The whole-repository consistency sweep, run as a test.

`tools/consistency_check.py` verifies that the simulation, the figures, the
reports and the reference list still agree with each other -- the physics
invariants, the published values the code carries, the absence of superseded
numbers, that every embedded figure exists and is newer both than the script
that makes it and than every library module that script imports, and that the
report numbering agrees across builder, templates and index.

It is wired in here so that a correction which lands in one place and not
another fails a test rather than reaching a report.  The checks that depend
on built artefacts (figures, HTML, PDFs) are skipped when those have not been
generated, so a fresh clone still passes.
"""

import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
CHECKER = ROOT / "tools" / "consistency_check.py"


@pytest.mark.skipif(not CHECKER.exists(), reason="checker not present")
def test_repository_is_self_consistent():
    if not (ROOT / "reports" / "index.html").exists():
        pytest.skip("reports not built")
    r = subprocess.run([sys.executable, str(CHECKER), "--verbose"],
                       capture_output=True, text=True, cwd=str(ROOT))
    assert r.returncode == 0, (
        "tools/consistency_check.py reports disagreements:\n" + r.stdout + r.stderr)
