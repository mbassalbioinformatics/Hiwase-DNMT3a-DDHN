#!/usr/bin/env bash
###############################################################################
# 03_shiba.sh                                                          (server)
#
# Alternative splicing with Shiba: ONE run over all RNA-seq libraries, each its
# own group, so that every event is defined once and every pairwise contrast can
# be taken from the same per-sample PSI matrix (03_rnaseq/08_shiba_downstream.R).
# Shiba tests differential PSI with Fisher's exact test on pooled junction read
# counts. Its reported dPSI is alternative - reference; direction is
# recomputed downstream from the per-sample PSI.
#
# Settings: strand 1 (read 1 sense, as featureCounts -s 1), unannotated events
# on (Shiba's counterpart of rMATS --novelSS), minimum 10 junction reads,
# individual PSI on, t-test off (one library per group).
#
# USAGE (from the repository root, Shiba on PATH):
#   bash 03_rnaseq/server/03_shiba.sh 2>&1 | tee shiba_run.log
# Output: $RNA_OUT/07_shiba/
###############################################################################
set -euo pipefail
source config/paths.sh
SHIBA_CMD=${SHIBA_CMD:-}
NTHREAD=16
MINREAD=10; FDR=0.05; DPSI=0.10
OUTROOT="$RNA_OUT/07_shiba"
SHEET="$(pwd)/config/samples.tsv"
field() { awk -F'\t' -v want="$1" -v cond="$2" -v val="$3" 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next}
          $h[cond] == val {print $h[want]}' "$SHEET"; }
SAMPLES=(); while IFS= read -r s; do SAMPLES+=("$s"); done < <(field sample rna_seq yes)
REF=$(for s in $(field sample subject donor); do [[ " ${SAMPLES[*]} " == *" $s "* ]] && echo "$s"; done | head -1)
ALT=$(field sample key DX | head -1)
mkdir -p "$OUTROOT"; cd "$OUTROOT"

if [[ -z "$SHIBA_CMD" ]]; then
  for c in shiba Shiba shiba.py; do command -v "$c" >/dev/null && { SHIBA_CMD="$c"; break; }; done
fi
[[ -n "$SHIBA_CMD" ]] || { echo "ERROR: Shiba not found (tried shiba, Shiba, shiba.py); set SHIBA_CMD" >&2; exit 1; }
SHIBA_VER=$(cat "$(dirname "$(readlink -f "$(command -v "$SHIBA_CMD")")")/VERSION" 2>/dev/null || echo unknown)
echo "Shiba ${SHIBA_VER}: $(command -v "$SHIBA_CMD")"
[[ -f "$GTF" ]] || { echo "ERROR: GTF not found: $GTF" >&2; exit 1; }
for s in "${SAMPLES[@]}"; do
  b="${BAM_DIR}/${s}${BAM_SUFFIX}"
  [[ -f "$b" ]] || { echo "ERROR: BAM not found: $b" >&2; exit 1; }
  [[ -f "${b}.bai" ]] || samtools index -@ "$NTHREAD" "$b"
done

EXPT="${OUTROOT}/experiment.tsv"
{ printf "sample\tbam\tgroup\n"; for s in "${SAMPLES[@]}"; do printf "%s\t%s/%s%s\t%s\n" "$s" "$BAM_DIR" "$s" "$BAM_SUFFIX" "$s"; done; } > "$EXPT"

# Shiba 0.8.2 configuration. strand uses the featureCounts convention (0, 1, 2, XS);
# processes come from the -p flag, not from the config.
CONFIG="${OUTROOT}/config.yaml"
cat > "$CONFIG" <<EOF
workdir:
  ${OUTROOT}/shiba_out
gtf:
  ${GTF}
experiment_table:
  ${EXPT}
unannotated:
  True
minimum_anchor_length:
  6
minimum_intron_length:
  70
maximum_intron_length:
  500000
strand:
  "1"
only_psi:
  False
only_psi_group:
  False
fdr:
  ${FDR}
delta_psi:
  ${DPSI}
reference_group:
  ${REF}
alternative_group:
  ${ALT}
minimum_reads:
  ${MINREAD}
individual_psi:
  True
ttest:
  False
excel:
  False
EOF

mkdir -p "${OUTROOT}/shiba_out"
set +e
"$SHIBA_CMD" -p "$NTHREAD" -v "$CONFIG" 2>&1 | tee "${OUTROOT}/shiba_stdout.log"
RC=${PIPESTATUS[0]}
set -e
[[ $RC -eq 0 ]] || { echo "Shiba FAILED (exit $RC)"; tail -40 "${OUTROOT}/shiba_stdout.log"; exit $RC; }

cat > "${OUTROOT}/RUN_INFO.txt" <<EOF
shiba_version=${SHIBA_VER}
command=${SHIBA_CMD} -p ${NTHREAD} -v ${CONFIG}
experiment_table=${EXPT}
event_set=one run, all libraries
samples=$(IFS=,; echo "${SAMPLES[*]}")
gtf=${GTF}
minimum_reads=${MINREAD}
date=$(date -Iseconds)
EOF
echo "done: ${OUTROOT}"
