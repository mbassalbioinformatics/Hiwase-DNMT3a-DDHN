#!/usr/bin/env bash
###############################################################################
# run_rnaseq.sh
# Runs the RNA-seq analysis (03_rnaseq/) in order, from the repository root,
# after the counts (03_rnaseq/server/01_featurecounts.sh) and, for the
# splicing steps, the rMATS and Shiba output have been copied from the server.
#   bash run_rnaseq.sh                 every step
#   bash run_rnaseq.sh 02 04           only the listed steps
# Log: rnaseq_run.log
###############################################################################
set -uo pipefail
cd "$(dirname "$0")"
STEPS=("02_dge_gene_level" "03_dge_5utr_level" "04_enrichment genes" "04_enrichment 5utr" "04b_enrichment_ma_plots"
       "05_pathway_ruler genes" "05b_bcv_and_seed_robustness" "06_top_signatures genes"
       "07_rmats_downstream" "08_shiba_downstream")
LOG=rnaseq_run.log
want=("$@")
for step in "${STEPS[@]}"; do
  s=${step%% *}; args=${step#"$s"}
  if (( ${#want[@]} )); then
    hit=0; for w in "${want[@]}"; do [[ "$s" == "$w"* ]] && hit=1; done
    (( hit )) || continue
  fi
  echo "$(date '+%F %T') START $step" | tee -a "$LOG"
  # shellcheck disable=SC2086
  if Rscript "03_rnaseq/$s.R" $args >> "$LOG" 2>&1; then
    echo "$(date '+%F %T') OK    $step" | tee -a "$LOG"
  else
    echo "$(date '+%F %T') FAIL  $step (see $LOG)" | tee -a "$LOG"; exit 1
  fi
done
echo "$(date '+%F %T') DONE" | tee -a "$LOG"
