#!/usr/bin/env bash
###############################################################################
# 05_annotate_collect.sh
#
# After 03 and 04: keeps the PASS calls of each case, flags each call's
# mappability (UMAP100: uniquely mappable 100-mers; SEGDUP: segmental
# duplication; taps_annotate_mappability.py, the same tracks as 04), annotates
# them with GATK Funcotator (gene, transcript, protein change), and gathers
# everything the analysis computer needs into one folder:
#   $VAR_OUT/transfer/  <case>.filtered.vcf.gz, <case>.pass.vcf.gz,
#                   <case>.pass.mappability.vcf.gz, <case>.funcotated.vcf.gz,
#                   contamination, segmentation and pileup tables, CNV segments
#                   (with the heterozygous-site allele counts), plots, logs and
#                   denoised copy ratios (gzipped), duplicate metrics, masking logs
# Copy $VAR_OUT/transfer/ to emseq/04_variants/ under DATA_ROOT (config/paths.R).
# Outputs older than the filtered VCF from 03 are redone.
#
# Run on the server:  bash 05_variants_server/05_annotate_collect.sh [--test chr21]
###############################################################################
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/sheet.sh"     # config/paths.sh and the sample sheet
OUT=${OUT:-$VAR_OUT}
UMAP="$UMAP_BED"; SEGDUP="$SEGDUP_BED"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAPF="$HERE/taps_annotate_mappability.py"
[[ "${1:-}" == "--test" ]] && OUT="$OUT/test_${2:?--test needs a contig}"
V="$OUT/mutect2/vcf"; X="$OUT/transfer"; mkdir -p "$X"/{vcf,cnv,logs/mutect2}
exec 9> "$OUT/.lock.05"; flock -n 9 || { echo "ERROR: 05 is already running on $OUT"; exit 3; }
rm -f "$OUT/.done.05"     # written again only when both cases are annotated
DOCKER="docker run --platform linux/amd64 --rm $DOCKER_MOUNTS $GATK_IMAGE"
GATK="$DOCKER gatk --java-options -Xmx16g"
FUNCO_DIR=${FUNCO_DIR:-$FUNCOTATOR_DIR}
[[ -d "$FUNCO_DIR" ]] || { echo "ERROR: Funcotator data sources not found: $FUNCO_DIR (set FUNCOTATOR_DIR in config/paths.sh)"; exit 1; }
for f in "$UMAP" "$SEGDUP" "$MAPF"; do [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }; done
if python3 -c "import pysam" 2>/dev/null; then PY="python3"; else PY="$DOCKER python"; fi
echo "$(date '+%F %T') Funcotator data sources: $FUNCO_DIR"
fresh() {  # fresh <output> <input...>: the output exists and is newer than every input
  local o="$1" i; shift; [[ -s "$o" ]] || return 1
  for i in "$@"; do [[ "$o" -nt "$i" ]] || return 1; done
}
FAILED=()

for c in DNMT3A TET2; do
  f="$V/$c.filtered.vcf.gz"; [[ -s "$f" ]] || { echo "FAIL $c: no $f (run 03)"; FAILED+=("$c"); continue; }
  pass="$X/vcf/$c.pass.vcf.gz"; mapd="$X/vcf/$c.pass.mappability.vcf.gz"; fun="$X/vcf/$c.funcotated.vcf.gz"
  fresh "$pass" "$f" || { rm -f "$pass"*; $DOCKER bash -c "bcftools view -f PASS -Oz -o $pass $f && bcftools index -t $pass" \
                         && echo "OK   $c PASS: $($DOCKER bash -c "bcftools view -H $pass | wc -l")" \
                         || { echo "FAIL PASS $c"; rm -f "$pass"*; FAILED+=("$c"); continue; }; }
  fresh "$mapd" "$pass" || { rm -f "$mapd"*; $PY "$MAPF" "$pass" "$UMAP" "$SEGDUP" "$mapd" > "$X/logs/mappability.$c.log" 2>&1 \
                         && echo "OK   $c mappability: $(tail -1 "$X/logs/mappability.$c.log")" \
                         || { echo "FAIL mappability $c (see $X/logs/mappability.$c.log)"; rm -f "$mapd"*; FAILED+=("$c"); continue; }; }
  fresh "$fun" "$mapd" || { rm -f "$fun"*; $GATK Funcotator -R "$REF" -V "$mapd" -O "$fun" --ref-version hg38 --data-sources-path "$FUNCO_DIR" \
                         --output-file-format VCF --remove-filtered-variants true > "$X/logs/funcotator.$c.log" 2>&1 \
                         && echo "OK   $c Funcotator" || { echo "FAIL Funcotator $c (see $X/logs/funcotator.$c.log)"; rm -f "$fun"*; FAILED+=("$c"); }; }
  cp -p "$f" "$f.tbi" "$X/vcf/" 2>/dev/null
done
cp -p "$V"/*.contamination.table "$V"/*.segments.table "$V"/*.pileups.table "$X/vcf/" 2>/dev/null
find "$OUT/mutect2/logs" -maxdepth 1 -name '*.log' ! -name 'mutect2.*' -exec cp -p {} "$X/logs/mutect2/" \; 2>/dev/null   # not the per-shard logs
for d in segments plots logs; do [[ -d "$OUT/cnv/$d" ]] && cp -rp "$OUT/cnv/$d" "$X/cnv/"; done
mkdir -p "$X/cnv/denoised"   # copy ratios per 1-kb bin, compressed, for the report figures
for f in "$OUT"/cnv/denoised/*.denoisedCR.tsv; do
  gz="$X/cnv/denoised/$(basename "$f").gz"; [[ -s "$f" ]] && ! fresh "$gz" "$f" && gzip -c "$f" > "$gz"
done
cp -p "$OUT"/logs/*duplicate_metrics.txt "$OUT"/logs/mask.*.log "$X/logs/" 2>/dev/null
(( ${#FAILED[@]} )) && { echo "$(date '+%F %T') FAILED: ${FAILED[*]} (re-run to retry)"; exit 1; }
date '+%F %T' > "$OUT/.done.05"
echo "$(date '+%F %T') DONE. Copy $X to emseq/04_variants/ on the analysis computer ($(du -sh "$X" | cut -f1))."
