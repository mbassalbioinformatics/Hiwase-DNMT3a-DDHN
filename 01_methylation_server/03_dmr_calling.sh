#!/usr/bin/env bash
###############################################################################
# 03_dmr_calling.sh                                                    (server)
#
# DMRs for every comparison in config/comparisons.tsv:
#   1. TAPS inputs per sample (02_taps_inputs.py)
#   2. per sample: decompress the fraction bedgraph; sortBed it for metilene
#   3. metilene (primary): the two samples' bedgraphs joined with
#        bedtools unionbedg -filler NA -header -names g1_<s1> g2_<s2>
#      then  metilene -m 3 -d 0.1 -t <cpu> -M 147 -f 1
#      (>= 3 CpGs, mean difference >= 0.1, <= 147 bp between CpGs, de novo mode)
#   4. Defiant (second caller):
#        defiant -b -c 0 -d 1 -G 147 -CpN 5 -p 0.05 -R
#      output <id>_defiant_DMR_147_c0_CpN5_d1_G147_p0.05_P10.tsv
#   5. summary: metilene DMRs at q < 0.05, gained and lost
#
# Comparisons are named g1_vs_g2; a positive difference means more methylation
# in g1 (metilene reports g1 - g2; Defiant reports g2 - g1, re-signed later).
#
# USAGE (from the repository root, in the environment with metilene, Defiant
# and bedtools):
#   nohup bash 01_methylation_server/03_dmr_calling.sh > dmr_run.log 2>&1 &
# Progress: $DMR_OUT/progress.txt. A re-run skips finished steps.
###############################################################################
set -uo pipefail
source config/paths.sh
CPU=${CPU:-30}
COMP=config/comparisons.tsv
REFORMAT=01_methylation_server/02_taps_inputs.py
OUT=$DMR_OUT
IN="$OUT/inputs"

mkdir -p "$IN" "$OUT/metilene/prep" "$OUT/defiant" "$OUT/logs"
log() { echo "$(date '+%F %T') $*" | tee -a "$OUT/progress.txt"; }

for t in python3 xz sortBed bedtools metilene defiant; do
  command -v "$t" >/dev/null || { echo "ERROR: '$t' is not on PATH" >&2; exit 1; }
done
[[ -s "$COMP" && -s "$REFORMAT" ]] || { echo "ERROR: $COMP or $REFORMAT missing" >&2; exit 1; }
free_gb=$(df -Pk "$OUT" | awk 'NR == 2 {print int($4 / 1048576)}')
(( ${free_gb:-0} < 60 )) && log "WARNING only ${free_gb} GB free on $OUT; about 40 GB is needed"

ROWS=(); SAMPLES=()
while IFS= read -r r; do ROWS+=("$r"); done < <(awk -F'\t' 'NR > 1 && $1 !~ /^#/ && NF >= 3' "$COMP")
while IFS= read -r s; do SAMPLES+=("$s"); done < <(printf '%s\n' "${ROWS[@]}" | cut -f2,3 | tr '\t' '\n' | sort -u)
log "START run: ${#ROWS[@]} comparisons, samples ${SAMPLES[*]}, CPU $CPU"

# 1. TAPS inputs ---------------------------------------------------------------
todo=()
for s in "${SAMPLES[@]}"; do
  [[ -s "$IN/$s.CpG.meth.TAPS.DefiantFormat3.bedgraph.xz" ]] || todo+=("$s")
done
if (( ${#todo[@]} )); then
  log "START inputs ${todo[*]}"
  if python3 "$REFORMAT" "$MS_OUT" "${todo[@]}" --out "$IN" > "$OUT/logs/01_inputs.log" 2>&1; then
    log "OK    inputs"; sed 's/^/        /' "$OUT/logs/01_inputs.log" | tee -a "$OUT/progress.txt"
  else
    log "FAIL  inputs, see logs/01_inputs.log"; cat "$OUT/logs/01_inputs.log"; exit 1
  fi
fi

# 2. decompress for Defiant, sort for metilene ---------------------------------
for s in "${SAMPLES[@]}"; do
  u="$IN/$s.CpG.meth.TAPS.DefiantFormat3.bedgraph"
  srt="$IN/$s.CpG.meth.TAPS.DefiantFormat3.sorted.bedgraph"
  [[ -s "$u" ]] || { xz -dc "$u.xz" > "$u.tmp" && mv "$u.tmp" "$u"; }
  [[ -s "$srt" ]] || { sortBed -i "$u" > "$srt.tmp" && mv "$srt.tmp" "$srt"; }
  [[ -s "$u" && -s "$srt" ]] || { log "FAIL  inputs for $s"; exit 1; }
done
log "OK    inputs decompressed and sorted"

# 3. metilene ------------------------------------------------------------------
for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r id g1 g2 _ <<< "$row"
  res="$OUT/metilene/${id}_metilene_DMRs_147_cpn3.tsv"
  [[ -s "$res" ]] && continue
  inp="$OUT/metilene/prep/${id}_metileneInput.bg"
  log "START metilene $id"
  if bedtools unionbedg -filler NA -header -names "g1_$g1" "g2_$g2" \
       -i "$IN/$g1.CpG.meth.TAPS.DefiantFormat3.sorted.bedgraph" "$IN/$g2.CpG.meth.TAPS.DefiantFormat3.sorted.bedgraph" > "$inp" \
     && metilene -m 3 -d 0.1 -t "$CPU" -M 147 -f 1 "$inp" > "$res.tmp" 2> "$OUT/metilene/${id}_metilene_DMRs.err" \
     && mv "$res.tmp" "$res"; then
    log "OK    metilene $id ($(wc -l < "$res" | tr -d ' ') regions)"
  else
    log "FAIL  metilene $id, see metilene/${id}_metilene_DMRs.err"
  fi
done

# 4. Defiant -------------------------------------------------------------------
for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r id g1 g2 _ <<< "$row"
  d="$OUT/defiant/$id"; mkdir -p "$d"
  [[ -f "$d/.done" ]] && continue
  log "START defiant $id"
  if ( cd "$d" && defiant -b -c 0 -cpu "$CPU" -d 1 -G 147 -CpN 5 -l "${id}_defiant_DMR_147" \
         -L "$g1,$g2" -p 0.05 -R \
         -i "$IN/$g1.CpG.meth.TAPS.DefiantFormat3.bedgraph" "$IN/$g2.CpG.meth.TAPS.DefiantFormat3.bedgraph" ) \
       > "$OUT/logs/defiant_$id.log" 2>&1; then
    touch "$d/.done"; log "OK    defiant $id"
  else
    log "FAIL  defiant $id, see logs/defiant_$id.log"
  fi
done

# 5. metilene summary ----------------------------------------------------------
{
  printf 'id\tregions\tq<0.05\thyper_in_g1\thypo_in_g1\tmedian_abs_diff\n'
  for row in "${ROWS[@]}"; do
    IFS=$'\t' read -r id _ <<< "$row"
    f="$OUT/metilene/${id}_metilene_DMRs_147_cpn3.tsv"
    [[ -s "$f" ]] || { printf '%s\tNA\n' "$id"; continue; }
    # columns: chr start stop q mean_diff(g1-g2) nCpG p_MWU p_2DKS mean_g1 mean_g2
    counts=$(awk -F'\t' '{n++} $4 < 0.05 {s++; if ($5 > 0) up++; else dn++}
                         END {printf "%d\t%d\t%d\t%d", n, s, up, dn}' "$f")
    med=$(awk -F'\t' '$4 < 0.05 {print ($5 < 0 ? -$5 : $5)}' "$f" | sort -g \
          | awk '{a[NR] = $1} END {if (NR) printf "%.3f", (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2); else printf "NA"}')
    printf '%s\t%s\t%s\n' "$id" "$counts" "$med"
  done
} > "$OUT/metilene_summary.tsv"
column -t "$OUT/metilene_summary.tsv" | tee -a "$OUT/progress.txt"
log "DONE"
