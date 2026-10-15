#!/usr/bin/env bash
###############################################################################
# 02_prepare_bams.sh
#
# Analysis-ready BAMs for the genome-wide variant steps (03-05), one per library:
#   1. taps_mask_methylation.py, one process per contig: the bases that report
#      methylation (the C of a CpG on top-strand reads, the G on bottom-strand
#      reads) become N with base quality 0; contigs and order of the GATK
#      reference dictionary; read group = sample. A BAM or contig part made by
#      an older version of the masking (version in its @PG line) is rebuilt.
#   2. samtools cat of the contigs in dictionary order (already coordinate-sorted)
#   3. GATK MarkDuplicates (skip with MARKDUP=0 if 01 shows the duplicates are
#      already removed), index
# Output: $VAR_OUT/bams/<sample>.taps.bam (+ .bai), duplicate metrics in $VAR_OUT/logs.
# Resumable: finished libraries and finished contigs are skipped.
#
# Run on the server (samtools on PATH; Docker with the GATK image of
# config/paths.sh; pysam on the host or, failing that, inside the GATK image):
#   bash 05_variants_server/02_prepare_bams.sh --test chr21     # rehearsal on one contig
#   bash 05_variants_server/02_prepare_bams.sh
###############################################################################
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/sheet.sh"     # config/paths.sh and the sample sheet
MS=${MS:-$MS_OUT}
OUT=${OUT:-$VAR_OUT}
DICT="${REF%.fa}.dict"
JOBS=${JOBS:-24}          # contigs masked at the same time
MARKDUP=${MARKDUP:-1}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MASK="$HERE/taps_mask_methylation.py"
SAMPLES=("${ALL[@]}")
DOCKER="docker run --platform linux/amd64 --rm $DOCKER_MOUNTS $GATK_IMAGE"
GATK="$DOCKER gatk --java-options -Xmx24g"

TEST=""
if [[ "${1:-}" == "--test" ]]; then TEST="${2:?--test needs a contig, e.g. chr21}"; OUT="$OUT/test_$TEST"; fi
B="$OUT/bams"; T="$OUT/tmp"; L="$OUT/logs"; mkdir -p "$B" "$T" "$L"
exec 9> "$OUT/.lock.02"; flock -n 9 || { echo "ERROR: 02 is already running on $OUT"; exit 3; }
rm -f "$OUT/.done.02"     # written again only when every library is finished
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH"; exit 1; }
[[ -s "$MASK" && -s "$DICT" && -s "$REF.fai" ]] || { echo "ERROR: $MASK, $DICT or $REF.fai missing"; exit 1; }
if python3 -c "import pysam, numpy" 2>/dev/null; then PY="python3"; else PY="$DOCKER python"; fi
echo "$(date '+%F %T') START prepare; python: $PY; out: $OUT"

mapfile -t CONTIGS < <(awk '/^@SQ/{for(i=1;i<=NF;i++) if($i~/^SN:/) print substr($i,4)}' "$DICT")
[[ -n "$TEST" ]] && CONTIGS=("$TEST")
MASK_VERSION=$(sed -n 's/^MASK_VERSION = "\([^"]*\)".*/\1/p' "$MASK")
current() {  # current <bam>: made by this version of the masking (fields of its @PG line, in any order)
  [[ -s "$1" ]] && samtools view -H "$1" 2>/dev/null | awk -F'\t' -v vn="VN:$MASK_VERSION" '
    $1 == "@PG" { id = v = 0; for (i = 2; i <= NF; i++) { if ($i == "ID:taps_mask_methylation") id = 1; if ($i == vn) v = 1 }; if (id && v) ok = 1 }
    END { exit !ok }'
}
echo "$(date '+%F %T') masking version $MASK_VERSION"
FAILED=()
fail() { echo "$(date '+%F %T') FAIL $*"; FAILED+=("$s"); }

for s in "${SAMPLES[@]}"; do
  bam="$MS/${s}_msuite2/Msuite2.final.bam"; final="$B/$s.taps.bam"
  if [[ -s "$final" && -s "$final.bai" ]]; then
    current "$final" && { echo "$(date '+%F %T') skip $s (done)"; continue; }
    echo "$(date '+%F %T') rebuild $s (masked by an older version)"
  fi
  [[ -s "$bam" ]] || { fail "$s: no $bam"; continue; }
  [[ -s "$bam.bai" || -s "${bam%.bam}.bai" ]] || samtools index -@ 8 "$bam"

  # 1. mask, contigs in parallel
  running=0
  for c in "${CONTIGS[@]}"; do
    part="$T/$s.$c.bam"
    current "$part" && continue
    ( $PY "$MASK" "$bam" "$REF" "$DICT" "$c" "$s" "$part.tmp" 2> "$L/mask.$s.$c.log" && mv "$part.tmp" "$part" \
      || { echo "$(date '+%F %T') FAIL mask $s $c (see $L/mask.$s.$c.log)"; rm -f "$part.tmp"; } ) &
    running=$((running + 1))
    if (( running >= JOBS )); then wait -n || true; running=$((running - 1)); fi
  done
  wait
  missing=0; for c in "${CONTIGS[@]}"; do [[ -s "$T/$s.$c.bam" ]] || missing=1; done
  (( missing )) && { fail "$s: some contigs not masked; re-run to retry"; continue; }

  # 2. concatenate in dictionary order (each part is sorted; the header order is the dictionary order)
  cat_bam="$T/$s.masked.bam"
  if ! current "$cat_bam"; then
    samtools cat -o "$cat_bam.tmp" $(for c in "${CONTIGS[@]}"; do echo "$T/$s.$c.bam"; done) && mv "$cat_bam.tmp" "$cat_bam" \
      || { fail "cat $s"; rm -f "$cat_bam.tmp"; continue; }
  fi

  # 3. duplicates, index
  if (( MARKDUP )); then   # spill files on the data disk, not in the container's /tmp
    mkdir -p "$T/markdup.$s"
    $GATK MarkDuplicates -I "$cat_bam" -O "$final.tmp.bam" -M "$L/$s.duplicate_metrics.txt" \
         --VALIDATION_STRINGENCY SILENT --CREATE_INDEX false --TMP_DIR "$T/markdup.$s" > "$L/markdup.$s.log" 2>&1 \
      && mv "$final.tmp.bam" "$final" || { fail "MarkDuplicates $s (see $L/markdup.$s.log)"; rm -f "$final.tmp.bam"; continue; }
    rm -rf "$T/markdup.$s"
  else
    mv "$cat_bam" "$final"
  fi
  samtools index -@ 8 "$final" || { fail "index $s"; continue; }
  masked=$(cat "$L"/mask.$s.*.log | awk -F'\t' '{split($4, a, " "); m += a[3]} END{print m}')
  echo "$(date '+%F %T') OK   $s -> $final (masked bases: $masked)"
  rm -f "$cat_bam" "$T/$s".*.bam
done
(( ${#FAILED[@]} )) && { echo "$(date '+%F %T') FAILED: ${FAILED[*]} (re-run to retry; finished libraries are kept)"; exit 1; }
date '+%F %T' > "$OUT/.done.02"
echo "$(date '+%F %T') DONE prepare"
