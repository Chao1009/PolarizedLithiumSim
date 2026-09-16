#!/bin/bash
# Write a ZERO-INSERTION copy of one craterlake compact file.
# Run INSIDE eic-shell:  bash ff_zero_insertion.sh <5x41|10x100|18x275> <outdir>
#
# The Roman-Pot silicon is held off the beam axis by four per-energy
# constants in compact/fields/beamline_<cfg>.xml (2.96 / 2.75 / 1.80 /
# 0.0 cm at 5x41, 0.71 / 0.55 / 0 / 0 at 10x100, 0.27 / 0 / 0 / 0 at
# 18x275).  Setting all four to zero slides the silicon onto the axis
# and leaves EVERY FIELD untouched, so the transport is the baseline
# one -- which is what makes it legitimate for measuring a transfer
# matrix element the insertion would otherwise hide.  Validated
# 2026-08-29 on R34 at 5x41: the control R12 came out 19.18 m in the
# zeroed file against 19.24 through the real insertion, a 0.3% shift,
# which is the fit scatter (tools/fullsim/README.md).
#
# Both substitutions are exact: the value pattern is "[^"]*" and not .*,
# so the two <!-- rough extrapolation --> comments on lines 63-64 of the
# 5x41 file survive; and the include pattern carries its ${DETECTOR_PATH}/
# prefix, without which the rewritten path is appended to the detector
# directory and no such file exists.  offset_OMD_x does not match
# _RP_section and is left alone.
set -euo pipefail
CFG="${1:?config: 5x41 | 10x100 | 18x275}"
OUT="${2:?outdir}"
mkdir -p "$OUT"
sed "s#\(offset_[a-z0-9_]*_RP_section\" value *=\)\"[^\"]*\"#\1\"0.0*cm\"#" \
    "$DETECTOR_PATH/compact/fields/beamline_$CFG.xml" \
    > "$OUT/beamline_${CFG}_zero.xml"
sed "s|\${DETECTOR_PATH}/compact/fields/beamline_$CFG.xml|$OUT/beamline_${CFG}_zero.xml|" \
    "$DETECTOR_PATH/epic_craterlake_$CFG.xml" \
    > "$OUT/epic_craterlake_${CFG}_zero.xml"
# grep -c exits 1 on zero matches, which under `set -e` would abort the
# script before it could say what it found; || true keeps the assertion
# below the thing that reports the failure.
n=$(grep -c 'RP_section" value *="0.0\*cm"' "$OUT/beamline_${CFG}_zero.xml" || true)
inc=$(grep -c "$OUT/beamline_${CFG}_zero.xml" "$OUT/epic_craterlake_${CFG}_zero.xml" || true)
echo "$CFG: $n/4 RP section offsets zeroed, $inc include rewritten -> $OUT/epic_craterlake_${CFG}_zero.xml"
if [ "$n" != 4 ] || [ "$inc" != 1 ]; then
  echo "ff_zero_insertion.sh: expected 4 zeroed offsets and 1 rewritten include" >&2
  exit 1
fi
