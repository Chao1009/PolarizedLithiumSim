"""polligen — doubly polarized e+(6,7)Li event generation (plans/05).

Step 5.A scope (plus reco/recopseudo/hfs, WP3): the physics kernel (spin-density matrices, inclusive
doubly polarized master cross section) and the inclusive spin-labeled
sampler with run-plan bookkeeping.  Tagged mode (5.B), reweighting (5.C)
and HepMC3 output (5.D) build on these modules.

polligen imports the existing fast-simulation package `polli_fastsim`
(structure functions, kinematics, FOM machinery) rather than duplicating
it; the sibling `fastsim/` directory is added to sys.path on import.
"""

import pathlib
import sys

_FASTSIM = pathlib.Path(__file__).resolve().parents[2] / "fastsim"
if str(_FASTSIM) not in sys.path:
    sys.path.insert(0, str(_FASTSIM))

import polli_fastsim  # noqa: E402,F401  (path check; polligen depends on it)

from . import (spin, xsec, bookkeeping, sample,  # noqa: E402
               estimators, tagged, coherent, reco, recopseudo, hfs)

__all__ = ["spin", "xsec", "bookkeeping", "sample", "estimators", "tagged",
           "coherent", "reco", "recopseudo", "hfs"]
__version__ = "0.1.0"
