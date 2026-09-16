#!/usr/bin/env python3
"""Tagged-FOM table: efficiency x purity x dilution per tagged channel
(plans/05 step 5.E; run 19b, 2026-09-15).

One row per tagged channel -- 6Li alpha, 7Li alpha and the d-p control --
with each factor taken from the place that owns it and nothing taken from
anywhere else:

  EFFICIENCY  the spectator-tagging acceptance, from the machinery of
              `fastsim/scripts/tagging_acceptance.py` (Report 3 Table 6):
              `spectator.spectator_lab_kinematics` folded with
              `farforward.acceptance_summary`, tagged = 1 - lost, i.e. ANY
              far-forward system.  Evaluated at the Yellow Report
              high-acceptance optics of each configuration and at the
              lithium tagging optics of Report 1 Section 6.1, and over the
              beta = 0.20/0.30/0.40 wave-function band (plans/05 step 5.B;
              the band is ONE-SIDED UPWARD -- see `--beta-band` there).

  DILUTION    the tagged two-cluster spin model's own dilution,
              `polligen.tagged.TaggedModel.tensor_dilution()` for the two
              S_c = 1 channels and `.vector_dilution()` for 7Li, whose
              channel spin is 1/2 and which therefore has no rank-2 moment
              at all.  Times the FOOTING factor the published observable
              is quoted on, named per row.  For 6Li that reproduces the
              published b1 scaling, (2/6) x 0.921949 = 0.307316
              (`polarized.b1_li6_from_deuteron(1.0)`, plans/08 D9), and
              the script asserts the agreement rather than restating it.

  PURITY      TWO columns, and both are honest about what does not exist.
              (a) the kinematic-window purity: the only |t|-window purity
              number this repository holds is the COHERENT channel's
              80-99% incoherent rejection, read off the e+Pb coherent-J/psi
              study arXiv:2108.01694 (PRD 104:114030) -- plans/06 SS, the
              band-level argument of plans/07 risk row 9.  It belongs to
              the intact-nucleus recoil and NOT to a spectator tag, whose
              selection is a rigidity/angle window and not a |t| fit, so it
              is printed as a labelled reference line under the table and
              the tagged rows carry `n/a`.
              (b) the incoherent-background purity per channel:
              UNAVAILABLE, and FLUKA-gated.  BeAGLE links FLUKA, whose
              licence is personal and per-user, so no A = 6, 7 breakup
              sample exists anywhere in this project
              (docs/reproduction_manual.md SS, plans/08 D5, plans/07
              row 9).  It is printed as `unavailable (FLUKA)` and is never
              filled with a guess -- which is why the FOM product below is
              reported as eps x D with the purity factor left explicitly
              open, and not as a single number that would look complete.

The two derived columns.  `eps x D` is the plan's literal product with the
purity factor omitted.  `eps x D^2` is the statistical figure of merit:
the projected error of every observable in `polli_fastsim.fom` goes as
1/(D x A) at fixed count (`Scenario.analyzing_power`) and the count itself
goes as eps, so the inverse-variance per unit luminosity is eps D^2.  Both
are printed because the plan asks for the first and the error bars are the
second.

WHY THIS SCRIPT MAY IMPORT `polligen` AND `polli_fastsim` MAY NOT.  The
package dependency runs one way, polligen -> polli_fastsim, and it is
test-guarded by an AST scan over `polli_fastsim/*.py`
(`fastsim/tests/test_fom_dilution.py`).  A SCRIPT is a consumer of both
packages and a member of neither, exactly as every `evgen/scripts/*.py`
is; importing the generator here inverts nothing.  The alternative --
copying 0.921949 into the fast simulation a second time -- is the drift
plans/08 D9 exists to prevent, and the assertion below is the guard.

Usage:  python3 scripts/tagged_fom_table.py [--nevents 400000] [--outdir out]
"""

import argparse
import pathlib
import sys

import numpy as np

_HERE = pathlib.Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[1]))                  # fastsim/
sys.path.insert(0, str(_HERE.parents[2] / "evgen"))        # evgen/

from polli_fastsim import beams  # noqa: E402
from polli_fastsim import farforward as ff  # noqa: E402
from polli_fastsim import polarized  # noqa: E402
from polli_fastsim import spectator as sp  # noqa: E402

from polligen import tagged  # noqa: E402

#: The beta band of plans/05 step 5.B, always quoted for a tagged number.
BETA_BAND = (0.20, 0.30, 0.40)
BETA_CENTRAL = 0.30

#: The one |t|-window purity number in this repository, and the channel it
#: belongs to: the e+Pb coherent-J/psi study arXiv:2108.01694
#: (PRD 104:114030) rejects 80-99% of incoherent events versus |t|
#: (plans/06, plans/07 risk row 9, docs/note_cos2phi_coherent_6Li.md).
#: It is the COHERENT intact-6Li recoil's, not a spectator tag's.
COHERENT_T_WINDOW_PURITY = (0.80, 0.99)

PURITY_UNAVAILABLE = "unavailable (FLUKA)"


class Row:
    """One tagged channel: its kinematic spectator, its spin model, the
    footing its published observable is quoted on."""

    def __init__(self, label, isotope, kin_channel, model_factory, rank,
                 footing, footing_note, observable):
        self.label = label
        self.isotope = isotope
        self.kin = kin_channel
        self.model_factory = model_factory
        self.rank = rank                  # 2 or 1
        self.footing = footing
        self.footing_note = footing_note
        self.observable = observable

    def dilution(self, beta=BETA_CENTRAL):
        """(model dilution, footing, product) at the stretched state."""
        m = tagged.TaggedModel(self.model_factory(beta))
        d = m.tensor_dilution() if self.rank == 2 else m.vector_dilution()
        return float(d), self.footing, float(d) * self.footing


ROWS = (
    Row("6Li alpha-tag (embedded d)", "6Li", sp.LI6_ALPHA_TAG,
        lambda b: tagged.li6_alpha_channel(beta=b), 2,
        polarized.LI6_B1_PER_NUCLEON,
        "2/6 per-nucleon (polarized.LI6_B1_PER_NUCLEON, plans/08 D9)",
        "tagged A_zz of the embedded deuteron (money plot 4)"),
    Row("7Li alpha-tag (quasi-free t)", "7Li", sp.LI7_ALPHA_TAG,
        lambda b: tagged.li7_alpha_channel(beta=b), 1,
        None,   # filled below: TRITON.eff_pol_p
        "P_p(t) = 0.86 (tagged.TRITON.eff_pol_p, PER NUCLEON and, since "
        "Z_t = 1, whole-triton too -- plans/05 SS5.4)",
        "tagged A_par on the quasi-free triton"),
    Row("d-p control (n struck, p tagged)", "d", sp.DEUTERON_P_TAG,
        lambda b: tagged.deuteron_channel(beta=b), 2,
        1.0, "1 (quoted on the deuteron itself)",
        "Cosyn-Weiss tagged A_zz (the deuteron limit)"),
)
ROWS[1].footing = tagged.TRITON.eff_pol_p


def efficiency(row, config, optics, n, beta, seed):
    """Tagged fraction = 1 - lost, ANY far-forward system -- the
    `tagging_acceptance.py` definition, on the same call."""
    kin = sp.spectator_lab_kinematics(row.kin,
                                      config.ion_momentum_per_nucleon, n,
                                      beta=beta,
                                      rng=np.random.default_rng(seed))
    acc = ff.acceptance_summary(kin["R"], kin["theta"], kin["pT"], optics,
                                phi=kin["phi"],
                                pot_config=ff.yr_config_key(config))
    return 1.0 - acc["lost"]


def optics_menu(config):
    return (("YR high-acceptance", ff.yr_optics(config, "high-acceptance")),
            ("tagging optics", ff.tagging_optics(config)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nevents", type=int, default=400_000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    L = ["# Tagged-FOM table (plans/05 step 5.E), %d spectators per cell, "
         "2026-09-15" % args.nevents,
         "# efficiency: 1 - lost (ANY far-forward system), the "
         "tagging_acceptance.py definition of Report 3 Table 6",
         "# dilution:   polligen.tagged.TaggedModel, rank 2 where the "
         "channel spin allows it and rank 1 for 7Li (S_c = 1/2)",
         "# purity:     see the two purity lines below the table -- one "
         "column is n/a by kinematics and the other is FLUKA-gated",
         ""]

    # --- the dilution block, and the guard against a second copy of 0.921949
    L.append("== dilution ==")
    L.append("%-34s %-6s %12s %12s %12s  %s"
             % ("channel", "rank", "D_model", "footing", "D_published",
                "observable / footing"))
    dil = {}
    for row in ROWS:
        d_model, foot, d_pub = row.dilution()
        dil[row.label] = (d_model, foot, d_pub)
        L.append("%-34s %-6d %12.7f %12.7f %12.7f  %s; footing = %s"
                 % (row.label, row.rank, d_model, foot, d_pub,
                    row.observable, row.footing_note))
    # the published 6Li b1 scaling must BE this row, not a transcription
    d6 = dil[ROWS[0].label]
    assert abs(d6[0] - polarized.LI6_B1_RANK2_TRANSFER) < 1e-4, d6
    assert abs(d6[2] - polarized.b1_li6_from_deuteron(1.0)) < 1e-6, d6
    L.append("  guard: D_published(6Li) = %.7f reproduces "
             "polarized.b1_li6_from_deuteron(1.0) = %.7f to 1e-6, and "
             "D_model = %.7f is polarized.LI6_B1_RANK2_TRANSFER = %.6f "
             "to 1e-4 (a quadrature against a stored constant)"
             % (d6[2], polarized.b1_li6_from_deuteron(1.0), d6[0],
                polarized.LI6_B1_RANK2_TRANSFER))
    L.append("")

    # --- efficiency, per configuration and optics, over the beta band
    L.append("== efficiency (beta band %s, central %.2f) =="
             % ("/".join("%.2f" % b for b in BETA_BAND), BETA_CENTRAL))
    L.append("%-34s %-22s %-22s %9s %9s %9s %9s %9s"
             % ("channel", "configuration", "optics", "eps(0.30)",
                "eps(0.20)", "eps(0.40)", "span", "L/L_HA"))
    eff = {}
    for row in ROWS:
        for config in beams.default_configs(row.isotope):
            for name, optics in optics_menu(config):
                vals = {b: efficiency(row, config, optics, args.nevents, b,
                                      args.seed) for b in BETA_BAND}
                lo, hi = min(vals.values()), max(vals.values())
                span = hi / lo if lo > 0 else float("nan")
                eff[(row.label, config.label(), name)] = vals
                L.append("%-34s %-22s %-22s %9.4f %9.4f %9.4f %9.2f %9.3f"
                         % (row.label, config.label(), name,
                            vals[0.30], vals[0.20], vals[0.40], span,
                            optics.lumi_fraction))
    L.append("")

    # --- the table itself: one row per channel
    L.append("== TAGGED-FOM TABLE (one row per channel) ==")
    L.append("%-34s %-22s %13s %11s %11s %13s %10s %10s"
             % ("channel", "optics", "eps (beta 0.30)", "eps span",
                "D_published", "purity |t|", "purity bkg", "eps x D"))
    for row in ROWS:
        for name, _o in optics_menu(beams.default_configs(row.isotope)[0]):
            per_cfg = [eff[(row.label, c.label(), name)]
                       for c in beams.default_configs(row.isotope)]
            cen = [v[BETA_CENTRAL] for v in per_cfg]
            allv = [x for v in per_cfg for x in v.values()]
            d_pub = dil[row.label][2]
            L.append("%-34s %-22s %13s %11s %11.6f %13s %10s %10s"
                     % (row.label, name,
                        "%.4f-%.4f" % (min(cen), max(cen)),
                        "%.4f-%.4f" % (min(allv), max(allv)),
                        d_pub, "n/a", PURITY_UNAVAILABLE,
                        "%.4f-%.4f" % (min(cen) * d_pub, max(cen) * d_pub)))
    L.append("")
    L.append("  eps x D^2 (statistical FOM per unit luminosity, "
             "delta A ~ 1/(D sqrt(eps N)); fom.Scenario.analyzing_power):")
    for row in ROWS:
        for name, _o in optics_menu(beams.default_configs(row.isotope)[0]):
            cen = [eff[(row.label, c.label(), name)][BETA_CENTRAL]
                   for c in beams.default_configs(row.isotope)]
            d_pub = dil[row.label][2]
            L.append("    %-34s %-22s %.6f - %.6f"
                     % (row.label, name, min(cen) * d_pub * d_pub,
                        max(cen) * d_pub * d_pub))
    L.append("")
    L.append("  purity |t| : n/a for a spectator tag -- the selection is a "
             "rigidity/angle window, not a |t| fit.  The only |t|-window "
             "purity this repository holds is the COHERENT intact-6Li "
             "recoil's %.2f-%.2f incoherent rejection, from the e+Pb "
             "coherent-J/psi study arXiv:2108.01694 (PRD 104:114030), "
             "plans/06 and plans/07 risk row 9 -- a different channel, "
             "quoted here as the reference and not as a column entry."
             % COHERENT_T_WINDOW_PURITY)
    L.append("  purity bkg : %s.  BeAGLE links FLUKA, whose licence is "
             "personal and per-user, so no A = 6, 7 breakup sample exists "
             "in this project (docs/reproduction_manual.md, plans/08 D5). "
             "The column is left empty rather than guessed, and the FOM "
             "product above is therefore eps x D with the purity factor "
             "explicitly open." % PURITY_UNAVAILABLE)

    text = "\n".join(L)
    (outdir / "tagged_fom_table.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
