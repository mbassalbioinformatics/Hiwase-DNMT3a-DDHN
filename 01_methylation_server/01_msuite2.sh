#!/usr/bin/env bash
###############################################################################
# 01_msuite2.sh                                                        (server)
#
# Alignment and per-CpG methylation calls for the TAPS libraries with Msuite2,
# one job per library:
#   msuite2 -x hg38 -m TAPS -3 -p <threads> --minsize 15 --maxins 1500
#           --cut-r1-head 8 --cut-r2-head 8
#   -m TAPS   TAPS chemistry: a methylated C is read as T
#   -3        three-letter alignment
# Every sequencing run of a library (initial and top-up) is passed together,
# comma-separated. Libraries are found by their SAGC library number (column
# `uln` of config/samples.tsv) in the FASTQ file names.
#
# Writes $MS_OUT/run_msuite2_batch.sh, which runs every job; run that next.
# Per library, $MS_OUT/<sample>_msuite2/ then holds Msuite2.final.bam,
# Msuite2.CpG.meth.call.xz (base counts per CpG), Msuite2.CpG.meth.bedgraph.xz
# and Msuite2.report/.
#
# USAGE (from the repository root): bash 01_methylation_server/01_msuite2.sh
###############################################################################
set -euo pipefail
source config/paths.sh
SHEET=config/samples.tsv
RUN_SCRIPT="$MS_OUT/run_msuite2_batch.sh"

for d in "${FASTQ_DIRS[@]}"; do [[ -d "$d" ]] || { echo "ERROR: FASTQ directory not found: $d" >&2; exit 1; }; done
mkdir -p "$MS_OUT"
for f in Msuite2.pUC19.size Msuite2.c.pUC19.size Msuite2.w.pUC19.size; do
  [[ -f "$PUC19_DIR/$f" ]] || echo "WARNING: $f is missing from $PUC19_DIR; copy it there before running" >&2
done

cat > "$RUN_SCRIPT" <<EOF
#!/usr/bin/env bash
set -euo pipefail
echo "=== Msuite2: \$(date) ==="
EOF

# sample <tab> uln, from the sample sheet
while IFS=$'\t' read -r sample uln; do
  r1=()
  while IFS= read -r -d '' f; do r1+=("$f"); done \
    < <(find "${FASTQ_DIRS[@]}" -maxdepth 1 -type f -name "*${uln}*R1*.f*q.gz" -print0 | sort -z)
  if (( ${#r1[@]} == 0 )); then echo "WARNING: no FASTQ for $sample ($uln)" >&2; continue; fi
  r2=()
  for f in "${r1[@]}"; do
    [[ -f "${f/_R1/_R2}" ]] || { echo "ERROR: no read-2 file for $f" >&2; exit 1; }
    r2+=("${f/_R1/_R2}")
  done
  r1j=$(IFS=,; echo "${r1[*]}"); r2j=$(IFS=,; echo "${r2[*]}")
  od="$MS_OUT/${sample}_msuite2"
  cat >> "$RUN_SCRIPT" <<EOF

echo "--- $sample ($uln): ${#r1[@]} run(s) ---"
$MSUITE2 -1 "$r1j" -2 "$r2j" -x hg38 -m TAPS -o "$od" -3 -p $THREADS --minsize 15 --maxins 1500 --cut-r1-head 8 --cut-r2-head 8
cd "$od"
cp "$PUC19_DIR"/Msuite2.pUC19.size "$PUC19_DIR"/Msuite2.c.pUC19.size "$PUC19_DIR"/Msuite2.w.pUC19.size ./
make || echo "WARNING: first make failed for $sample; retrying"
cp "$PUC19_DIR"/Msuite2.pUC19.size "$PUC19_DIR"/Msuite2.c.pUC19.size "$PUC19_DIR"/Msuite2.w.pUC19.size ./
make || echo "WARNING: second make failed for $sample"
make clean || true
xz -vfT0 *.call *.bedgraph || true
cd "$MS_OUT"
EOF
done < <(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $h["uln"] != "" {print $h["sample"] "\t" $h["uln"]}' "$SHEET")

echo 'echo "=== Msuite2 done: $(date) ==="' >> "$RUN_SCRIPT"
chmod +x "$RUN_SCRIPT"
echo "Written: $RUN_SCRIPT"
echo "Run it:  nohup bash \"$RUN_SCRIPT\" > \"$MS_OUT/run_msuite2_batch.log\" 2>&1 &"
