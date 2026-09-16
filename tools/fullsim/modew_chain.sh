#!/bin/bash
# Run one weighted HepMC3 file through npsim -> EICrecon.
# INSIDE eic-shell:  bash modew_chain.sh <input.hepmc> <outdir> [config] [N] [seed]
#
#   singularity exec --env S=$S $SIF bash \
#       $R/tools/fullsim/modew_chain.sh $S/modew_100.hepmc $S/sim 18x275 100
#
# The seed is the fifth argument and defaults to 20260916.  It is not
# cosmetic: npsim without --random.seed seeds from the clock, and the SAME
# 100 primaries shot twice (files differing only in their two HepMC3 W
# lines, the particle records md5-identical) gave 335 against 547
# ForwardRomanPotHits, 26 against 30 reconstructed ZDC neutrals and a
# reconstructed-x median residual of -0.45 against -0.18 in the y < 0.1 bin
# -- i.e. at 100 events the shower seed moves the reco-level numbers by more
# than anything this smoke is trying to measure.  Pass seed=0 to let npsim
# seed itself.
#
# The input is the Mode-W file of tools/analysis/modew_beagle_hepmc.py:
# official BeAGLE events with the polarized weight named in the HepMC3
# weights vector and written FIRST in it -- npsim copies only weights()[0]
# into EventHeader.weight, so a weight appended after the sample's own
# `default` arrives here as a constant 1.0.  The two traps this encodes
# are the ones the first reconstructed run found
# (tools/fullsim/README.md, 2026-09-15):
#
#   * plain epic_craterlake.xml loads compact/fields/beamline_5x41.xml, so
#     the configuration is named explicitly on the SIMULATION side, and
#   * eicrecon does NOT inherit npsim's compact file -- thisepic.sh exports
#     DETECTOR_CONFIG=epic, so without -Pdd4hep:xml_files it reconstructs
#     against epic.xml's 5x41 beamline -- named explicitly on the
#     RECONSTRUCTION side too.
#
# npsim writes calibrations/, fieldmaps/ and gdml/ caches into $PWD, so the
# script cds into <outdir> and never runs from the repository.
#
# The afterburner leg is NOT here.  The official BeAGLE EVGEN files are
# already afterburned (GenRunInfo ab_afterburner_is_used = 1,
# ab_crossing_angle = 0.025, and the beam rows carry -25 mrad); abconv
# re-applied would be a second crossing angle.  `abconv -p <preset>
# --exit-ca` is the check, and it exits 0 saying so.
set -euo pipefail
IN="${1:?input .hepmc}"
OUT="${2:?outdir}"
CFG="${3:-18x275}"
N="${4:-100}"
SEED="${5:-20260916}"

# The script cds into $OUT (below), so both paths have to be absolute before
# it does, or a relative argument silently resolves against the new $PWD.
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
case "$IN" in /*) ;; *) IN="$PWD/$IN" ;; esac
[ -f "$IN" ] || { echo "no such input file: $IN" >&2; exit 1; }

: "${DETECTOR_PATH:?source /opt/detector/epic-main/bin/thisepic.sh first}"
COMPACT="$DETECTOR_PATH/epic_craterlake_$CFG.xml"
[ -f "$COMPACT" ] || { echo "no such compact file: $COMPACT" >&2; exit 1; }

cd "$OUT"                       # the npsim caches land here, not in the repo

{
  echo "=== modew_chain.sh $(date -Is)"
  echo "input     $IN"
  echo "compact   $COMPACT"
  echo "events    $N"
  if [ "$SEED" = "0" ]; then
    echo "seed      <none: npsim seeds from the clock>"
  else
    echo "seed      $SEED"
  fi
  npsim --version 2>&1 | head -1 || true
  eicrecon --version 2>&1 | head -1 || true
} | tee "$OUT/chain_versions.txt"

NPSIM_SEED_ARGS=()
if [ "$SEED" != "0" ]; then
  NPSIM_SEED_ARGS=(--random.seed "$SEED")
fi

npsim --compactFile "$COMPACT" \
      -N "$N" --inputFiles "$IN" \
      --physics.list FTFP_BERT --part.minimalKineticEnergy "100*MeV" \
      "${NPSIM_SEED_ARGS[@]}" \
      --outputFile "$OUT/sim.edm4hep.root"

eicrecon -Pdd4hep:xml_files="$COMPACT" \
         -Ppodio:output_file="$OUT/reco.edm4eic.root" \
         "$OUT/sim.edm4hep.root"

ls -l "$OUT/sim.edm4hep.root" "$OUT/reco.edm4eic.root"
