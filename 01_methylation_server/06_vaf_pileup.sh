#!/usr/bin/env bash
###############################################################################
# 06_vaf_pileup.sh                                                     (server)
#
# Base counts per original DNA strand over the coding exons (+/- 10 bp) of a
# myeloid and clonal-haematopoiesis driver panel (vaf_panel_hg38.bed: 85 genes,
# DNMT3A and TET2 included), from each library's Msuite2 alignment. Used on the
# analysis computer to find the driver mutations and their VAFs. Strand-aware,
# because TAPS reads a methylated C as T: at a CpG, C>T is read from the
# opposite strand.
#
# samtools mpileup without a reference (bases printed as letters), MAPQ >= 10,
# base quality >= 20, no BAQ, read 1 and read 2 piled up separately so that
# vaf_pileup_counts.py can assign each base to the original top (OT) or bottom
# (OB) strand: OT = read 1 forward + read 2 reverse; OB = read 1 reverse + read 2
# forward.
#
# USAGE (from the repository root, samtools and python3 on PATH):
#   nohup bash 01_methylation_server/06_vaf_pileup.sh > vaf.log 2>&1 &
# Output: $DMR_OUT/vaf/<sample>.strand_counts.tsv.gz, vaf/ALL.strand_counts.tsv.gz,
#         vaf/read1_orientation_DNMT3A_region.tsv (a check that libraries are directional)
###############################################################################
set -uo pipefail
source config/paths.sh
BED=01_methylation_server/vaf_panel_hg38.bed
PARSE=01_methylation_server/vaf_pileup_counts.py
OUT=$DMR_OUT
V="$OUT/vaf"; mkdir -p "$V"
SAMPLES=()
while IFS= read -r s; do SAMPLES+=("$s"); done < <(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next}
  $h["methylation"] == "analysed" {print $h["sample"]}' config/samples.tsv)
log() { echo "$(date '+%F %T') $*" | tee -a "$OUT/progress.txt"; }
for t in samtools python3; do command -v $t >/dev/null || { echo "ERROR: $t not on PATH" >&2; exit 1; }; done
[[ -s "$BED" && -s "$PARSE" ]] || { echo "ERROR: $BED or $PARSE missing" >&2; exit 1; }
log "START vaf pileup ($(samtools --version | head -1)): ${SAMPLES[*]}"

for s in "${SAMPLES[@]}"; do
  bam="$MS_OUT/${s}_msuite2/Msuite2.final.bam"
  res="$V/$s.strand_counts.tsv.gz"
  [[ -s "$res" ]] && continue
  [[ -s "$bam" ]] || { log "FAIL  vaf $s: no $bam"; continue; }
  [[ -s "$bam.bai" || -s "${bam%.bam}.bai" ]] || samtools index -@ 8 "$bam"
  for m in 1 2; do
    flag=$([[ $m == 1 ]] && echo 64 || echo 128)
    samtools mpileup -l "$BED" -q 10 -Q 20 -B -d 100000 --rf "$flag" "$bam" > "$V/$s.r$m.mpileup" 2> "$V/$s.r$m.err"
  done
  if python3 "$PARSE" "$s" "$V/$s.r1.mpileup" "$V/$s.r2.mpileup" | gzip > "$res"; then
    rm -f "$V/$s.r1.mpileup" "$V/$s.r2.mpileup"
    log "OK    vaf $s ($(zcat "$res" | wc -l) positions)"
  else
    log "FAIL  vaf $s, see $V/$s.r*.err"
  fi
done
zcat "$V"/*.strand_counts.tsv.gz | gzip > "$V/ALL.strand_counts.tsv.gz"
for s in "${SAMPLES[@]}"; do
  bam="$MS_OUT/${s}_msuite2/Msuite2.final.bam"
  f=$(samtools view -c -f 64 -F 16 "$bam" chr2:25200000-25300000); r=$(samtools view -c -f 80 "$bam" chr2:25200000-25300000)
  echo -e "$s\tread1_forward\t$f\tread1_reverse\t$r"
done > "$V/read1_orientation_DNMT3A_region.tsv"
log "DONE vaf pileup"
