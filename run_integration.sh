#!/usr/bin/env bash
###############################################################################
# run_integration.sh
# Runs the methylation x RNA-seq integration (04_integration/) in order, from
# the repository root. Needs config/paths.R and the config tables, the Msuite2
# per-CpG files and DMR results copied from the server, and the RNA-seq results.
#   bash run_integration.sh            every step
#   bash run_integration.sh 06 10      only the listed steps
# The locus-plot galleries (05, 10) are the slow part (~20 min).
# Log: integration_run.log
###############################################################################
set -uo pipefail
cd "$(dirname "$0")"
STEPS=(01_region_methylation 02_tf_activity 03_promoter_meth_vs_expression 04_global_direction
       07_gene_panels 06_dmr_master 08_dmr_plots 09_hic_loops 11_homer_summary
       13_vaf_from_taps_reads 14_allelic_imbalance
       05_locus_plot_gallery 10_dmr_locus_gallery)
LOG=integration_run.log
want=("$@")
for s in "${STEPS[@]}"; do
  if (( ${#want[@]} )); then
    hit=0; for w in "${want[@]}"; do [[ "$s" == "$w"* ]] && hit=1; done
    (( hit )) || continue
  fi
  echo "$(date '+%F %T') START $s" | tee -a "$LOG"
  if Rscript "04_integration/$s.R" >> "$LOG" 2>&1; then
    echo "$(date '+%F %T') OK    $s" | tee -a "$LOG"
  else
    echo "$(date '+%F %T') FAIL  $s (see $LOG)" | tee -a "$LOG"; exit 1
  fi
done
echo "$(date '+%F %T') DONE" | tee -a "$LOG"
