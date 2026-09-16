#!/bin/bash
# Far-forward particle-gun rigidity scan (plan 2.2.1 seed).
# Run INSIDE eic-shell:  singularity exec $SIF bash ff_gun_scan.sh <outdir> [nev]
# Shoots fragments along the ion axis (-25 mrad in x) with momenta chosen
# so the rigidity ratio R (w.r.t. the 275-GeV-proton optics of the default
# craterlake fields) matches the Li-fragment cases of plans/03 #2.2.
#
# THIS SCRIPT VARIES RIGIDITY, NOT THE PHASE-SPACE POINT.  Every config
# below is shot at ONE direction and from the origin, so what it measures
# is ROUTING -- which subsystem sees which fragment.  The scan that walks
# the IP coordinates themselves is ff_transfer_scan.py (2026-09-15): the
# angle x'/y' through the (p_T, phi) grid, the DISPLACEMENT x/y through a
# HepMC3 vertex, and the rigidity through --delta, which between them fit
# the full first-order transfer matrix R11/R12/R21/R22, R33/R34/R43/R44,
# D and D'.  npsim's gun can do the angle (--gun.direction) but not the
# vertex, which is why that scan feeds HepMC and this one does not.
#
# Read the output with ff_gun_hits.py; --b0-zdc prints the B0 and ZDC hit
# positions that the counting table below reports only as a yes/no.
set -uo pipefail
OUT="${1:-/tmp/ff_gun}"; NEV="${2:-100}"
mkdir -p "$OUT"

for s in /opt/detector/epic-main/bin/thisepic.sh /opt/detector/setup.sh \
         /opt/detector/epic-24.07.0/setup.sh; do
    [ -f "$s" ] && { source "$s"; break; }
done
echo "DETECTOR_PATH=$DETECTOR_PATH"
# IMPORTANT: plain epic_craterlake.xml defaults to the 5x41 beamline
# fields; the rigidity-scan momenta below assume the 275-GeV-proton
# optics, so use the 18x275 configuration explicitly.
XML="$DETECTOR_PATH/epic_craterlake_18x275.xml"
[ -f "$XML" ] || XML="$DETECTOR_PATH/epic_craterlake.xml"
echo "using $XML"

# tag  particle  p_total[GeV]  comment (R = p/(Z*275))
CONFIGS="
alpha_R0857   alpha    471.4  7Li-alpha R=0.857 -> expect RP
alpha_R100    alpha    550.0  6Li-alpha R=1.0 pT~0 -> expect lost
alpha_R100pt3 alpha    550.0  6Li-alpha R=1.0 pT=0.3 -> expect RP tail
triton_R1286  triton   353.6  7Li-triton R=1.286 -> expect lost
proton_R050   proton   137.5  spectator p R=0.5 -> expect OMD
deuteron_R0857 deuteron 235.7 7Li-d R=0.857 -> expect RP
neutron_ZDC   neutron  117.9  neutron -> expect ZDC
"
echo "$CONFIGS" | while read -r tag part ptot comment; do
    [ -z "${tag:-}" ] && continue
    dir="(-0.025,0,1.0)"
    if [ "$tag" = "alpha_R100pt3" ]; then
        # add ~0.545 mrad to emulate pT = 0.3 GeV at 550 GeV
        dir="(-0.0244,0,1.0)"
    fi
    n=$NEV; [ "$part" = "neutron" ] && n=30
    echo "=== $tag : $part p=$ptot GeV dir=$dir ($comment)"
    # minimalKineticEnergy kills soft shower secondaries: this is a
    # ROUTING check (where does the primary go), not a response study
    npsim --compactFile "$XML" \
          --enableGun \
          --gun.particle "$part" \
          --gun.momentumMin "$ptot*GeV" --gun.momentumMax "$ptot*GeV" \
          --gun.direction "$dir" \
          --numberOfEvents "$n" \
          --physics.list FTFP_BERT \
          --part.minimalKineticEnergy "100*MeV" \
          --outputFile "$OUT/gun_$tag.edm4hep.root" \
          > "$OUT/gun_$tag.log" 2>&1
    echo "    exit=$? $(ls -la $OUT/gun_$tag.edm4hep.root 2>/dev/null | awk '{print $5}') bytes"
done
echo DONE
