#!/usr/bin/env bash
###############################################################################
# 06_site_counts.sh
#
# Background check for the pooled clone test of 04_integration/17_genomewide_snvs.R:
# read counts, split by read orientation, at the clone mutations of every case in
# every library (taps_site_counts.py, on the masked BAMs from 02). A person's
# libraries read at another person's clone sites give the error background,
# since those mutations cannot be in them.
#   sites   site_counts_sites.tsv, written by 04_integration/17_genomewide_snvs.R
#           (integration/variants/); copy it to the server and give its path
#   output  $VAR_OUT/site_counts/ (sites.tsv, site_counts.tsv, site_counts.log);
#           copy the folder to emseq/04_variants/site_counts/ under DATA_ROOT
#
# Run on the server:  bash 05_variants_server/06_site_counts.sh <sites.tsv>
###############################################################################
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/sheet.sh"     # config/paths.sh and the sample sheet
OUT=${OUT:-$VAR_OUT}
SITES=$(realpath "${1:?give the sites file written by 04_integration/17_genomewide_snvs.R}")
PYF="$HERE/taps_site_counts.py"
SAMPLES=("${ALL[@]}")
D="$OUT/site_counts"; mkdir -p "$D"
exec 9> "$OUT/.lock.06"; flock -n 9 || { echo "ERROR: 06 is already running on $OUT"; exit 3; }
DOCKER="docker run --platform linux/amd64 --rm $DOCKER_MOUNTS $GATK_IMAGE"
if python3 -c "import pysam" 2>/dev/null; then PY="python3"; else PY="$DOCKER python"; fi
for f in "$SITES" "$PYF"; do [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }; done
args=()
for s in "${SAMPLES[@]}"; do
  b="$OUT/bams/$s.taps.bam"; [[ -s "$b" ]] || { echo "ERROR: missing $b (run 02)"; exit 1; }
  args+=("$s=$b")
done
cp -p "$SITES" "$D/sites.tsv"
echo "$(date '+%F %T') START site counts: $(($(wc -l < "$D/sites.tsv") - 1)) sites x ${#SAMPLES[@]} libraries"
if $PY "$PYF" "$D/sites.tsv" "$D/site_counts.tsv" "${args[@]}" > "$D/site_counts.log" 2>&1; then
  echo "$(date '+%F %T') DONE -> $D/site_counts.tsv ($(($(wc -l < "$D/site_counts.tsv") - 1)) rows)"
  echo "Copy $D to emseq/04_variants/site_counts/ on the analysis computer"
else
  echo "$(date '+%F %T') FAIL (see $D/site_counts.log)"; exit 1
fi
