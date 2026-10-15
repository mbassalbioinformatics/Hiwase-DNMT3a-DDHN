#!/usr/bin/env bash
###############################################################################
# run_variants.sh
#
# Runs the genome-wide variant steps one after the other,
#   02 prepare BAMs -> 03 Mutect2 -> 04 copy number / LOH -> 05 annotate and collect
# and stops at the first step that fails.
#
# Picking up after a stop: every step skips what is already finished and redoes
# only what a change upstream made out of date. Fix the cause of the stop and
# start this script again; it carries on from where it stopped.
#   - A failed step is retried once (RETRIES=1) before the run stops. Finished
#     parts are skipped, so the retry costs little and gets past one-off
#     failures (a shard out of memory, a Docker hiccup).
#   - A step counts as finished only when it exits 0 and leaves its marker
#     <run folder>/.done.<step>.
#   - One run at a time per output folder. If a step is already running (for
#     example started by hand), this waits for it to end, then runs the step
#     itself, which skips whatever that run finished.
#   - Output of each step: <run folder>/run_logs/<step>.log, copied into
#     <run folder>/transfer/run_logs/ whenever this script ends.
#
# Run on the server, from the repository root:
#   bash 05_variants_server/run_variants.sh --test chr21     # rehearsal on chr21
#   nohup bash 05_variants_server/run_variants.sh > variants_run.log 2>&1 &
# Options:  --test <contig>     rehearsal on one contig (output in $VAR_OUT/test_<contig>)
#           --stop-after <NN>   stop after that step, e.g. --stop-after 02
# To stop a run: pkill -f run_variants.sh; pkill -f '05_variants_server/'; then
#   docker ps, and docker stop any GATK container still running. Start again to resume.
###############################################################################
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/sheet.sh"     # config/paths.sh and the sample sheet
OUT=${OUT:-$VAR_OUT}
RETRIES=${RETRIES:-1}
WAIT=${WAIT:-300}      # seconds between checks while another copy of a step runs
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEPS=(02_prepare_bams.sh 03_somatic_calls.sh 04_cnv_loh.sh 05_annotate_collect.sh)

TEST_ARGS=(); STOP_AFTER=""
while (( $# )); do
  case "$1" in
    --test)       TEST_ARGS=(--test "${2:?--test needs a contig, e.g. chr21}"); shift 2 ;;
    --stop-after) STOP_AFTER="${2:?--stop-after needs a step number, e.g. 02}"; shift 2 ;;
    *) echo "unknown option: $1 (options: --test <contig>, --stop-after <step>)"; exit 2 ;;
  esac
done
RUN="$OUT"; (( ${#TEST_ARGS[@]} )) && RUN="$OUT/test_${TEST_ARGS[1]}"
mkdir -p "$RUN/run_logs"
log() { echo "$(date '+%F %T') $*"; }
exec 8> "$RUN/.lock.run_all"
flock -n 8 || { log "ERROR: this pipeline is already running on $RUN"; exit 3; }
for f in "${STEPS[@]}"; do [[ -s "$HERE/$f" ]] || { log "ERROR: $HERE/$f missing"; exit 1; }; done
save_logs() { mkdir -p "$RUN/transfer"; cp -rp "$RUN/run_logs" "$RUN/transfer/" 2>/dev/null; }
trap save_logs EXIT

others() {  # pids of other running copies of step script $1 in this mode (same --test contig, or a full run)
  ps -eo pid=,args= | awk -v f="$1" -v t="${TEST_ARGS[*]:-}" -v me="$$" '
    $1 != me && $2 ~ /(^|\/)bash$/ && index($0, f) {
      test = index($0, "--test ") > 0
      if ((t == "" && !test) || (t != "" && index($0 " ", t " ") > 0)) printf "%s ", $1 }'
}

log "START ${TEST_ARGS[*]:-full run} -> $RUN${STOP_AFTER:+ (stop after $STOP_AFTER)}"
for f in "${STEPS[@]}"; do n=${f%%_*}; [[ -f "$RUN/.done.$n" ]] && log "  last finished run of $n: $(cat "$RUN/.done.$n")"; done
T0=$SECONDS
for f in "${STEPS[@]}"; do
  n=${f%%_*}; steplog="$RUN/run_logs/$n.log"
  while pids=$(others "$f"); [[ -n "$pids" ]]; do
    log "waiting: $n is already running (pid $pids); checking again in $WAIT s"; sleep "$WAIT"
  done
  attempt=0
  while :; do
    attempt=$((attempt + 1)); t0=$SECONDS
    log "=== $n, attempt $attempt (log: $steplog)"
    echo "===== $(date '+%F %T') attempt $attempt: bash $f ${TEST_ARGS[*]:-}" >> "$steplog"
    bash "$HERE/$f" ${TEST_ARGS[@]+"${TEST_ARGS[@]}"} 2>&1 | tee -a "$steplog"
    rc=${PIPESTATUS[0]}
    if (( rc == 0 )) && [[ -f "$RUN/.done.$n" ]]; then log "=== $n finished ($(( (SECONDS - t0) / 60 )) min)"; break; fi
    if (( rc == 0 )); then rc=1; log "$n ended without its finish marker $RUN/.done.$n"; fi
    (( rc == 3 )) && { log "STOPPED: another run holds the lock of $n on $RUN"; exit 3; }
    (( rc >= 128 )) && { log "STOPPED: $n was interrupted (signal $((rc - 128))). Start this script again to resume."; exit "$rc"; }
    if (( attempt > RETRIES )); then
      log "STOPPED at $n (exit $rc). Fix the cause (see $steplog), then start this script again: finished work is kept."
      exit "$rc"
    fi
    log "$n failed (exit $rc); trying again (finished parts are skipped)"
  done
  [[ "$n" == "$STOP_AFTER" ]] && { log "stopping after $n, as asked"; exit 0; }
done
log "ALL DONE in $(( (SECONDS - T0) / 3600 )) h $(( (SECONDS - T0) % 3600 / 60 )) min. Copy $RUN/transfer/ to emseq/04_variants/ on the analysis computer."
