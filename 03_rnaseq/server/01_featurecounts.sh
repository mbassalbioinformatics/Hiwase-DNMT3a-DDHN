#!/usr/bin/env bash
###############################################################################
# 01_featurecounts.sh                                                  (server)
#
# Gene and 5'UTR counts with featureCounts (Subread) against Ensembl 107:
#   featureCounts -p --countReadPairs -M -O -J -f -d 20 -s 1 -g gene_id
#                 --extraAttributes gene_name,transcript_id
#   5'UTR level: the same with -t five_prime_utr
# -s 1: the libraries are stranded with read 1 in the sense orientation. Before
# counting, the script checks this on one library by counting it with -s 0, 1
# and 2: -s 2 (antisense) should assign only a small fraction of the reads.
# --countReadPairs counts fragments (Subread >= 2.0.2; in older versions -p
# alone does, and the flag is omitted).
# -f reports one row per exon feature; 03_rnaseq/01_load_counts.R removes the
# duplicate rows of exons shared between transcripts before summing per gene.
#
# USAGE (from the repository root): bash 03_rnaseq/server/01_featurecounts.sh
# Output: $RNA_OUT/GENES.counts_ensembl_FC and $RNA_OUT/5UTR.counts_ensembl_FC
###############################################################################
set -euo pipefail
source config/paths.sh
STRAND=1
SAMPLES=()
while IFS= read -r s; do SAMPLES+=("$s"); done < <(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next}
  $h["rna_seq"] == "yes" {print $h["sample"]}' config/samples.tsv)
BAMS=(); for s in "${SAMPLES[@]}"; do BAMS+=("${BAM_DIR}/${s}${BAM_SUFFIX}"); done

command -v featureCounts >/dev/null || { echo "ERROR: featureCounts not on PATH" >&2; exit 1; }
FC_VER=$(featureCounts -v 2>&1 | grep -oE 'v[0-9.]+' | head -1)
echo "featureCounts ${FC_VER}"
if featureCounts --help 2>&1 | grep -q -- "--countReadPairs"; then CRP="--countReadPairs"; else CRP=""; fi
[[ -f "$GTF" ]] || { echo "ERROR: GTF not found: $GTF" >&2; exit 1; }
for b in "${BAMS[@]}"; do [[ -f "$b" ]] || { echo "ERROR: missing $b" >&2; exit 1; }; done
mkdir -p "$RNA_OUT"

echo "--- strand check on ${SAMPLES[0]} (share of reads assigned) ---"
TMPD=$(mktemp -d)
for s in 0 1 2; do
  featureCounts -p -T "$THREADS" -s "$s" -a "$GTF" -o "${TMPD}/s${s}.txt" "${BAMS[0]}" > /dev/null 2>&1 || true
  if [[ -f "${TMPD}/s${s}.txt.summary" ]]; then
    A=$(awk '$1 == "Assigned" {print $2}' "${TMPD}/s${s}.txt.summary")
    T=$(awk 'NR > 1 {sum += $2} END {print sum}' "${TMPD}/s${s}.txt.summary")
    printf "  -s %d : %s / %s (%.1f%%)\n" "$s" "$A" "$T" "$(awk -v a="$A" -v t="$T" 'BEGIN {print 100 * a / t}')"
  fi
done
rm -rf "$TMPD"
echo "  expected for read-1-sense libraries: -s 1 close to -s 0, -s 2 small"

echo "== gene level =="
featureCounts -g gene_id --extraAttributes gene_name,transcript_id \
  -p $CRP -M -O -J -f -d 20 -s "$STRAND" -T "$THREADS" -a "$GTF" \
  -o "${RNA_OUT}/GENES.counts_ensembl_FC" "${BAMS[@]}"

echo "== 5'UTR level =="
featureCounts -t five_prime_utr -g gene_id --extraAttributes gene_name,transcript_id \
  -p $CRP -M -O -J -f -d 20 -s "$STRAND" -T "$THREADS" -a "$GTF" \
  -o "${RNA_OUT}/5UTR.counts_ensembl_FC" "${BAMS[@]}"

for f in "${RNA_OUT}"/{GENES,5UTR}.counts_ensembl_FC.summary; do echo "--- $(basename "$f") ---"; cat "$f"; done
