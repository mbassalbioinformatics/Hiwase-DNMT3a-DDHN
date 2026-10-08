#!/usr/bin/env bash
###############################################################################
# 05_homer_annotation.sh                                               (server)
#
# HOMER annotation of the DMRs:
#   annotatePeaks.pl <dmrs.bed> hg38 -annStats <stats> > <annotated>
# for the metilene DMRs (q < 0.05, the primary set) and for the metilene DMRs
# also called by Defiant ("both"). The BEDs are written on the analysis
# computer by 04_integration/06_dmr_master.R (integration/dmr/bed/<set>/) and
# copied to $DMR_BEDS; names carry the comparison and direction (gained / lost),
# the score column the change (g1 - g2).
#
# USAGE (from the repository root, HOMER on PATH; the analysis used HOMER v4.11.1):
#   DMR_BEDS=/path/to/dmr/bed nohup bash 01_methylation_server/05_homer_annotation.sh > annot.log 2>&1 &
# Output: $DMR_OUT/homer_annotation/<set>/<comparison>.annotated.txt and .annStats.txt
###############################################################################
set -uo pipefail
source config/paths.sh
OUT=$DMR_OUT
DMR_BEDS=${DMR_BEDS:-$OUT/dmr_beds}
# HOMER writes large temporary files into the working directory
mkdir -p "$OUT/homer_tmp" && cd "$OUT/homer_tmp" || exit 1
log() { echo "$(date '+%F %T') $*" | tee -a "$OUT/progress.txt"; }
command -v annotatePeaks.pl >/dev/null || { echo "ERROR: annotatePeaks.pl not on PATH" >&2; exit 1; }
for set in metilene both; do
  mkdir -p "$OUT/homer_annotation/$set"
  for bed in "$DMR_BEDS/$set"/*.bed; do
    id=$(basename "$bed" .bed)
    [[ -s "$bed" ]] || { log "SKIP  annotate $set/$id: no DMRs"; continue; }
    a="$OUT/homer_annotation/$set/$id.annotated.txt"; s="$OUT/homer_annotation/$set/$id.annStats.txt"
    if annotatePeaks.pl "$bed" hg38 -annStats "$s" > "$a" 2> "$OUT/homer_annotation/$set/$id.log"; then
      log "OK    annotate $set/$id"
    else
      log "FAIL  annotate $set/$id, see homer_annotation/$set/$id.log"
    fi
  done
done
log "DONE homer annotation"
