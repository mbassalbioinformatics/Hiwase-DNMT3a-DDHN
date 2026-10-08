#!/usr/bin/env bash
###############################################################################
# 02_rmats.sh                                                          (server)
#
# Alternative splicing with rMATS-turbo: the six pairwise contrasts of
# config/rnaseq_contrasts.tsv, plus one grouped run (all recipient timepoints
# against the donor). The grouped run is the only one with replicates; its
# "replicates" are different timepoints, so it is conservative.
#
# Settings: --libType fr-secondstrand (read 1 sense, as featureCounts -s 1),
# read length detected from the BAMs, --variable-read-length (the reads are
# quality-trimmed), --novelSS, --allow-clipping where the installed rMATS
# supports them.
# Direction: the contrast numerator is always --b1, so IncLevelDifference =
# PSI(numerator) - PSI(denominator); 03_rnaseq/07_rmats_downstream.R re-derives
# it from the per-sample PSI and stops if they disagree. Settings of every run
# are written to <contrast>/RUN_INFO.txt.
#
# USAGE (from the repository root, rmats.py and samtools on PATH):
#   bash 03_rnaseq/server/02_rmats.sh 2>&1 | tee rmats_run.log
# Output: $RNA_OUT/06_rmats/<contrast>/
###############################################################################
set -euo pipefail
source config/paths.sh
RMATS_CMD=${RMATS_CMD:-rmats.py}
LIBTYPE="fr-secondstrand"
NOVEL_SS=1
OUTROOT="$RNA_OUT/06_rmats"
SHEET="$(pwd)/config/samples.tsv"; CONTRASTS="$(pwd)/config/rnaseq_contrasts.tsv"
field() { awk -F'\t' -v want="$1" -v cond="$2" -v val="$3" 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next}
          $h[cond] == val {print $h[want]}' "$SHEET"; }
SAMPLES=(); while IFS= read -r s; do SAMPLES+=("$s"); done < <(field sample rna_seq yes)
PATIENTS=(); while IFS= read -r s; do PATIENTS+=("$s"); done < <(field sample subject patient)
DONOR=$(field sample subject donor | while read -r s; do [[ " ${SAMPLES[*]} " == *" $s "* ]] && echo "$s"; done | head -1)
mkdir -p "$OUTROOT"; cd "$OUTROOT"

command -v "$RMATS_CMD" >/dev/null || { echo "ERROR: $RMATS_CMD not on PATH" >&2; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH" >&2; exit 1; }
[[ -f "$GTF" ]] || { echo "ERROR: GTF not found: $GTF" >&2; exit 1; }
for s in "${SAMPLES[@]}"; do
  b="${BAM_DIR}/${s}${BAM_SUFFIX}"
  [[ -f "$b" ]] || { echo "ERROR: BAM not found: $b" >&2; exit 1; }
  [[ -f "${b}.bai" || -f "${b%.bam}.bai" ]] || samtools index -@ "$THREADS" "$b"
done

set +o pipefail                      # samtools | head ends with SIGPIPE
READLEN=$(samtools view "${BAM_DIR}/${SAMPLES[0]}${BAM_SUFFIX}" 2>/dev/null | head -200000 \
          | awk '{print length($10)}' | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')
set -o pipefail
[[ -n "${READLEN:-}" ]] || { echo "ERROR: could not determine read length" >&2; exit 1; }
if [[ $(samtools view -c -f 1 "${BAM_DIR}/${SAMPLES[0]}${BAM_SUFFIX}") -gt 0 ]]; then TFLAG=paired; else TFLAG=single; fi
echo "read length ${READLEN}; ${TFLAG}; libType ${LIBTYPE}"

RM_HELP=$("$RMATS_CMD" --help 2>&1 || true)
have_flag() { grep -q -- "$1" <<< "$RM_HELP"; }
VARLEN=""; have_flag "--variable-read-length" && VARLEN="--variable-read-length"
NOVEL_FLAG=""; [[ $NOVEL_SS -eq 1 ]] && have_flag "--novelSS" && NOVEL_FLAG="--novelSS"
CLIP_FLAG=""; have_flag "--allow-clipping" && CLIP_FLAG="--allow-clipping"

mkdir -p bamlists
for s in "${SAMPLES[@]}"; do echo "${BAM_DIR}/${s}${BAM_SUFFIX}" > "bamlists/${s}.txt"; done
(IFS=,; for s in "${PATIENTS[@]}"; do printf '%s\n' "${BAM_DIR}/${s}${BAM_SUFFIX}"; done | paste -sd, -) > bamlists/PATIENT_ALL.txt

run_rmats() {
  local cid="$1" b1="$2" b2="$3" od="${OUTROOT}/$1" td="${OUTROOT}/tmp_$1"
  [[ -s "${od}/SE.MATS.JC.txt" ]] && { echo "== ${cid}: done already =="; return 0; }
  echo "== ${cid}: $(cat "$b1") vs $(cat "$b2") =="
  mkdir -p "$od" "$td"
  local cmd=("$RMATS_CMD" --b1 "$b1" --b2 "$b2" --gtf "$GTF" -t "$TFLAG" --readLength "$READLEN"
             --nthread "$THREADS" --od "$od" --tmp "$td" --libType "$LIBTYPE")
  [[ -n "$VARLEN" ]] && cmd+=("$VARLEN"); [[ -n "$NOVEL_FLAG" ]] && cmd+=("$NOVEL_FLAG"); [[ -n "$CLIP_FLAG" ]] && cmd+=("$CLIP_FLAG")
  set +e; "${cmd[@]}" 2>&1 | tee "${od}/rmats_stdout.log"; local rc=${PIPESTATUS[0]}; set -e
  if [[ $rc -ne 0 ]]; then echo "${cid} FAILED (exit $rc)"; tail -30 "${od}/rmats_stdout.log"; exit $rc; fi
  cat > "${od}/RUN_INFO.txt" <<EOF
contrast_id=${cid}
b1_numerator=$(tr '\n' ',' < "$b1")
b2_denominator=$(tr '\n' ',' < "$b2")
readLength=${READLEN}
variable_read_length=$([[ -n "$VARLEN" ]] && echo yes || echo no)
libType=${LIBTYPE}
readType=${TFLAG}
novelSS=${NOVEL_SS}
gtf=${GTF}
rmats=$("$RMATS_CMD" --version 2>&1 | head -1)
date=$(date -Iseconds)
sign_convention=IncLevelDifference == PSI(b1) - PSI(b2) == PSI(numerator) - PSI(denominator)
EOF
  rm -rf "$td"
}

while IFS=$'\t' read -r id num den _; do
  [[ "$id" == id || -z "$id" ]] && continue
  run_rmats "$id" "bamlists/${num}.txt" "bamlists/${den}.txt"
done < "$CONTRASTS"
run_rmats G1_patient_vs_donor bamlists/PATIENT_ALL.txt "bamlists/${DONOR}.txt"

printf "%-24s %8s %8s %8s %8s %8s\n" contrast SE A5SS A3SS MXE RI
for d in "${OUTROOT}"/C*_* "${OUTROOT}"/G1_*; do
  [[ -d "$d" ]] || continue
  line=$(printf "%-24s" "$(basename "$d")")
  for ev in SE A5SS A3SS MXE RI; do
    f="${d}/${ev}.MATS.JC.txt"; line+=$(printf "%8s" "$([[ -s "$f" ]] && echo $(( $(wc -l < "$f") - 1 )) || echo NA)")
  done
  echo "$line"
done
