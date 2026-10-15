#!/usr/bin/env bash
###############################################################################
# 01_inspect_bams.sh
#
# Read-only check of the Msuite2 alignments before the genome-wide variant steps
# (02-05): sort order, read groups, the aligner, duplicates, mapping qualities,
# and whether every contig of the GATK reference dictionary is in each BAM with
# the same length. Nothing is written except the report.
#
# Run on the server (samtools on PATH), from the repository root:
#   bash 05_variants_server/01_inspect_bams.sh 2>&1 | tee inspect_bams.txt
###############################################################################
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/sheet.sh"     # config/paths.sh and the sample sheet
MS=${MS:-$MS_OUT}
DICT="${REF%.fa}.dict"
SAMPLES=("${ALL[@]}")
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH"; exit 1; }
[[ -s "$DICT" && -s "$REF.fai" ]] || { echo "ERROR: $DICT or $REF.fai missing"; exit 1; }
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
sq() { awk '/^@SQ/{for(i=1;i<=NF;i++){if($i~/^SN:/)n=substr($i,4);if($i~/^LN:/)l=substr($i,4)}; print n"\t"l}' "$1"; }
sq "$DICT" > "$tmp/ref"
echo "reference: $REF ($(wc -l < "$tmp/ref") contigs: $(cut -f1 "$tmp/ref" | tr '\n' ' '))"
echo "samtools: $(samtools --version | head -1)"
for s in "${SAMPLES[@]}"; do
  bam="$MS/${s}_msuite2/Msuite2.final.bam"
  echo; echo "=== $s"
  [[ -s "$bam" ]] || { echo "  MISSING: $bam"; continue; }
  echo "  file: $bam ($(du -h "$bam" | cut -f1))"
  idx=no; [[ -s "$bam.bai" || -s "${bam%.bam}.bai" ]] && idx=yes; echo "  index: $idx"
  samtools view -H "$bam" > "$tmp/hdr"
  echo "  @HD: $(grep '^@HD' "$tmp/hdr" | cut -f2- | tr '\t' ' ')"
  echo "  @RG lines: $(grep -c '^@RG' "$tmp/hdr")"
  echo "  @PG: $(grep '^@PG' "$tmp/hdr" | awk -F'\t' '{for(i=2;i<=NF;i++) if($i~/^(ID|PN|VN):/) printf "%s ", $i; printf "| "}' | cut -c1-300)"
  sq "$tmp/hdr" > "$tmp/bam"
  echo "  contigs in BAM: $(wc -l < "$tmp/bam")"
  awk 'NR==FNR{ref[$1]=$2; next} ($1 in ref){seen[$1]=1; if (ref[$1]!=$2) print "  LENGTH MISMATCH:", $1, "reference", ref[$1], "BAM", $2}
       END{for (c in ref) if (!(c in seen)) print "  MISSING FROM BAM:", c}' "$tmp/ref" "$tmp/bam"
  echo "  BAM contigs not in the reference (first 10): $(awk 'NR==FNR{ref[$1]=1; next} !($1 in ref){print $1}' "$tmp/ref" "$tmp/bam" | head -10 | tr '\n' ' ')"
  if [[ $idx == yes ]]; then
    samtools idxstats "$bam" | awk 'NR==FNR{ref[$1]=1; next} {tot+=$3; if ($1 in ref) prim+=$3; if ($1=="chrL") lam+=$3}
         END{printf "  mapped reads: %d; on reference contigs %.2f%%; on chrL (lambda) %d\n", tot, 100*prim/(tot>0?tot:1), lam}' "$tmp/ref" -
  fi
  samtools flagstat -@ 8 "$bam" | grep -E "in total|secondary|supplementary|duplicates|properly paired" | sed 's/^/  flagstat: /'
  echo "  reads at chr21:30,000,000 (flag MAPQ CIGAR):"; samtools view "$bam" chr21:30000000-30000400 2>/dev/null | head -4 | cut -f2,5,6 | sed 's/^/    /'
  echo "  MAPQ on chr21:30-31 Mb: $(samtools view "$bam" chr21:30000000-31000000 2>/dev/null | awk '{n++; if($5>=20) q++} END{printf "%d reads, %.1f%% MAPQ >= 20", n, 100*q/(n>0?n:1)}')"
done
echo; echo "Done. Check this report before running 02_prepare_bams.sh."
