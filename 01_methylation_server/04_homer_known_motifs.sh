#!/usr/bin/env bash
###############################################################################
# 04_homer_known_motifs.sh                                             (server)
#
# Known-motif enrichment (HOMER) in the metilene DMRs (q < 0.05), separately for
# DMRs that lost and gained methylation, per comparison:
#   findMotifsGenome.pl <bed> hg38 <out> -size given -nomotif -cpg
# (known motifs only, region sizes as given, CpG-normalised background, as DMRs
# are CpG-rich). Sets with fewer than MIN_DMR DMRs are skipped.
#
# USAGE (from the repository root, HOMER on PATH; the analysis used HOMER v4.11.1):
#   nohup bash 01_methylation_server/04_homer_known_motifs.sh > homer.log 2>&1 &
# Output: $DMR_OUT/homer/<comparison>_<lost|gained>/knownResults.txt
###############################################################################
set -uo pipefail
source config/paths.sh
CPU=${CPU:-30}
MIN_DMR=${MIN_DMR:-5}
COMP="$(pwd)/config/comparisons.tsv"
OUT=$DMR_OUT
H="$OUT/homer"; mkdir -p "$H/bed"
# HOMER writes large temporary files into the working directory
mkdir -p "$OUT/homer_tmp" && cd "$OUT/homer_tmp" || exit 1
log() { echo "$(date '+%F %T') $*" | tee -a "$OUT/progress.txt"; }
command -v findMotifsGenome.pl >/dev/null || { echo "ERROR: findMotifsGenome.pl not on PATH" >&2; exit 1; }

while IFS=$'\t' read -r id g1 g2 _; do
  [[ "$id" == id || "$id" == \#* || -z "$id" ]] && continue
  f="$OUT/metilene/${id}_metilene_DMRs_147_cpn3.tsv"
  [[ -s "$f" ]] || { log "SKIP  homer $id: no metilene output"; continue; }
  for dir in lost gained; do
    bed="$H/bed/${id}_${dir}.bed"
    if [[ $dir == lost ]]; then cond='$5 < 0'; else cond='$5 > 0'; fi
    awk -F'\t' -v OFS='\t' -v id="$id" -v d="$dir" "\$4 < 0.05 && $cond {n++; print \$1, \$2, \$3, id\"_\"d\"_\"n, \".\", \"+\"}" "$f" > "$bed"
    n=$(wc -l < "$bed" | tr -d ' ')
    res="$H/${id}_${dir}"
    if (( n < MIN_DMR )); then log "SKIP  homer ${id}_${dir}: $n DMRs (< $MIN_DMR)"; continue; fi
    [[ -s "$res/knownResults.txt" ]] && continue
    log "START homer ${id}_${dir} ($n DMRs)"
    if findMotifsGenome.pl "$bed" "$HOMER_GENOME" "$res" -size given -nomotif -cpg -p "$CPU" > "$H/${id}_${dir}.log" 2>&1; then
      log "OK    homer ${id}_${dir}"
    else
      log "FAIL  homer ${id}_${dir}, see homer/${id}_${dir}.log"
    fi
  done
done < "$COMP"
log "DONE homer"
