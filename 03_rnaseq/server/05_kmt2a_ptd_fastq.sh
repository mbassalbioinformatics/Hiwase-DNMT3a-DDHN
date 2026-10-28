#!/usr/bin/env bash
###############################################################################
# 05_kmt2a_ptd_fastq.sh  --  KMT2A-PTD junction reads searched in the FASTQs
#                            (every read, including those STAR discarded)
#
# From the repository root, on the server, after 04_kmt2a_ptd_junctions.sh (it
# reuses that run's probes):
#   bash 03_rnaseq/server/05_kmt2a_ptd_fastq.sh 2>&1 | tee kmt2a_ptd_fastq_run.log
#
# Then copy $RNA_OUT/09_kmt2a_ptd/fastq_search/ to the analysis computer as
#   rnaseq/09_kmt2a_ptd/fastq_search/
#
# WHY
#   04_ searched the aligned BAMs and found no PTD junction read, but the BAMs
#   are blind to part of them: they keep no unmapped reads, and STAR drops a
#   read unless at least 66% of it aligns in one piece, so a read that crosses a
#   PTD junction near its middle never reaches the BAM. The FASTQs given to
#   STAR hold every read. Searching them with the same probes counts PTD
#   junction reads and normal junction reads with the same sensitivity.
#
# HOW
#   The FASTQs are those on STAR's own command line in each BAM header
#   (--readFilesIn), unless FASTQ_PATTERN is set. Each file is decompressed and
#   searched for the 78 probes of 04_ (66 PTD + 12 normal junctions, 15 + 15
#   bases, both orientations) with grep -F; the matches are re-checked in the
#   read sequence. A fragment counts once per junction even if both mates
#   match. The files are searched in parallel (MAX_JOBS at a time).
#
# OUTPUT (in PTD_DIR/fastq_search)
#   fastq_ptd_summary.tsv      per library: reads searched, PTD reads, the top PTD
#                              junction, PTD_e10_e2 against normal exon 10 -> 11
#   fastq_junction_matrix.tsv  every probe with any read, by library
#   fastq_junction_counts.tsv  every probe, every library (fragments, distinct sequences)
#   <sample>/supporting_reads.tsv   each matching read: probe, mate, sequence
#   RUN_INFO.txt
###############################################################################
set -euo pipefail
export LC_ALL=C

# ============================ CONFIGURE ME ==================================
source config/paths.sh
SHEET="$(pwd)/config/samples.tsv"
field() { awk -F'\t' -v want="$1" -v cond="$2" -v val="$3" 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next}
          $h[cond] == val {print $h[want]}' "$SHEET"; }
SAMPLES=(); while IFS= read -r s; do SAMPLES+=("$s"); done < <(field sample rna_seq yes)
CONDA_ENV="${RNA_CONDA_ENV:-}"      # conda environment providing samtools; empty = use PATH
MAX_JOBS=8                          # FASTQ files searched at the same time
PTD_DIR="$RNA_OUT/09_kmt2a_ptd"     # output of 04_kmt2a_ptd_junctions.sh
# Empty: the FASTQs on STAR's command line in each BAM header. Otherwise a pattern with
# {S} for the sample, e.g. "/path/to/rnaseq/02_processed_fq/{S}_*R[12]*.fq.gz"
FASTQ_PATTERN="${RNA_FASTQ_PATTERN:-}"
PY="python3"
# ============================================================================

OUT="$PTD_DIR/fastq_search"
mkdir -p "$OUT"
cd "$OUT"

echo "############ PREFLIGHT ############"
if [[ -n "$CONDA_ENV" ]] && command -v conda >/dev/null 2>&1; then
  set +u   # conda's scripts read variables that may be unset
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV" || echo "WARNING: could not activate '$CONDA_ENV'; using PATH"
  set -u
fi
command -v samtools >/dev/null 2>&1 || { echo "FATAL: samtools not on PATH"; exit 1; }
command -v "$PY" >/dev/null 2>&1 || { echo "FATAL: $PY not on PATH"; exit 1; }
for f in probes.tsv patterns.txt; do
  [[ -s "$PTD_DIR/$f" ]] || { echo "FATAL: $PTD_DIR/$f missing; run 03_rnaseq/server/04_kmt2a_ptd_junctions.sh first"; exit 1; }
done
if command -v pigz >/dev/null 2>&1; then GZ="pigz -dc -p 4"; else GZ="gzip -dc"; fi
echo "samtools     : $(samtools --version | head -1)"
echo "decompress   : $GZ"
echo "probes       : $PTD_DIR/probes.tsv ($(($(wc -l < "$PTD_DIR/probes.tsv") - 1)) junctions)"

# FASTQs per sample: from STAR's command line in the BAM header, or FASTQ_PATTERN
fastqs_of() {   # prints one FASTQ path per line
  local s="$1"
  if [[ -n "$FASTQ_PATTERN" ]]; then
    local p="${FASTQ_PATTERN//\{S\}/$s}"
    compgen -G "$p" | sort
    return
  fi
  samtools view -H "${BAM_DIR}/${s}${BAM_SUFFIX}" | awk -F'\t' '
    $1 == "@PG" && $0 ~ /ID:STAR/ {
      for (i = 2; i <= NF; i++) if ($i ~ /^CL:/) cl = substr($i, 4)
      n = split(cl, w, /[ ]+/); take = 0
      for (j = 1; j <= n; j++) {
        if (w[j] == "--readFilesIn") { take = 1; continue }
        if (w[j] ~ /^--/) take = 0
        if (take && w[j] != "") { m = split(w[j], f, ","); for (k = 1; k <= m; k++) print f[k] }
      }
      exit
    }'
}
declare -a JOB_SAMPLE JOB_MATE JOB_FILE
for s in "${SAMPLES[@]}"; do
  b="${BAM_DIR}/${s}${BAM_SUFFIX}"
  [[ -f "$b" ]] || { echo "FATAL: BAM not found: $b"; exit 1; }
  FQ=()
  while IFS= read -r line; do [[ -n "$line" ]] && FQ+=("$line"); done < <(fastqs_of "$s")
  [[ "${#FQ[@]}" -gt 0 ]] || { echo "FATAL: no FASTQ found for $s (set FASTQ_PATTERN)"; exit 1; }
  mate=0
  for f in "${FQ[@]}"; do
    mate=$((mate + 1))
    [[ -f "$f" ]] || { echo "FATAL: FASTQ not found for $s: $f"; echo "       set FASTQ_PATTERN to where the FASTQs are now"; exit 1; }
    JOB_SAMPLE+=("$s"); JOB_MATE+=("$mate"); JOB_FILE+=("$f")
    echo "  $s file $mate: $f ($(du -h "$f" | cut -f1))"
  done
done

echo
echo "############ SEARCH (${#JOB_FILE[@]} files, $MAX_JOBS at a time) ############"
search_one() {   # search_one <sample> <mate> <fastq>
  local s="$1" m="$2" f="$3" dir="$OUT/$1"
  mkdir -p "$dir"
  local dec="cat"; [[ "$f" == *.gz ]] && dec="$GZ"
  set +e
  # one line per read (name, sequence); awk also counts the reads, inside the pipeline
  $dec "$f" | paste - - - - \
    | awk -F'\t' -v out="$dir/reads_file${m}.txt" '{ print $1 "\t" $2 } END { print NR > out }' \
    | grep -F -f "$PTD_DIR/patterns.txt" > "$dir/hits_file${m}.txt"
  local st=("${PIPESTATUS[@]}")
  set -e
  if [[ "${st[0]}" -ne 0 || "${st[2]}" -ne 0 || "${st[3]}" -gt 1 ]]; then
    echo "FAILED $s file $m (decompress ${st[0]}, awk ${st[2]}, grep ${st[3]})" > "$dir/status_file${m}.txt"
  else
    echo "OK" > "$dir/status_file${m}.txt"
  fi
}
for i in "${!JOB_FILE[@]}"; do
  while [[ $(jobs -rp | wc -l) -ge "$MAX_JOBS" ]]; do sleep 2; done
  search_one "${JOB_SAMPLE[$i]}" "${JOB_MATE[$i]}" "${JOB_FILE[$i]}" &
done
wait
fail=0
for i in "${!JOB_FILE[@]}"; do
  st=$(cat "$OUT/${JOB_SAMPLE[$i]}/status_file${JOB_MATE[$i]}.txt" 2>/dev/null || echo "MISSING")
  [[ "$st" == "OK" ]] || { echo "  $st (${JOB_FILE[$i]})"; fail=1; }
done
[[ "$fail" -eq 0 ]] || { echo "FATAL: a search failed"; exit 1; }
echo "  all files searched"

# ---------------------------------------------------------------------------
# counting (Python, standard library only)
# ---------------------------------------------------------------------------
HELPER="$OUT/fastq_count_helper.py"
cat > "$HELPER" <<'PYEOF'
"""Count KMT2A junction reads from the FASTQ hits of 05_kmt2a_ptd_fastq.sh."""
import collections
import glob
import os
import re
import sys

probes_tsv, out = sys.argv[1], sys.argv[2]
samples = sys.argv[3:]
P = [l.rstrip("\n").split("\t") for l in open(probes_tsv)][1:]
rows = []
summary = []
for s in samples:
    d = os.path.join(out, s)
    support = collections.defaultdict(dict)      # probe -> fragment -> sequence
    recs = []
    reads = 0
    for rf in sorted(glob.glob(os.path.join(d, "reads_file*.txt"))):
        reads += int(open(rf).read().split()[0])
    for hf in sorted(glob.glob(os.path.join(d, "hits_file*.txt"))):
        mate = "R" + re.sub(r"\D", "", os.path.basename(hf))
        for line in open(hf):
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            frag = re.sub(r"/[12]$", "", f[0].lstrip("@").split()[0])
            seq = f[1].upper()
            for name, t, i, j, k, krc in P:
                if k in seq or krc in seq:
                    support[name].setdefault(frag, seq)
                    recs.append((s, name, t, frag, mate, seq))
    with open(os.path.join(d, "supporting_reads.tsv"), "w") as fo:
        fo.write("sample\tprobe\ttype\tfragment\tfile\tsequence\n")
        fo.writelines("\t".join(r) + "\n" for r in recs)
    c = {}
    for name, t, i, j, k, krc in P:
        q = support.get(name, {})
        c[name] = (t, len(q))
        rows.append((s, name, t, i, j, len(q), len(set(q.values()))))
    ptd = {k: v[1] for k, v in c.items() if v[0] == "PTD"}
    tot = sum(ptd.values())
    top = max(ptd, key=ptd.get) if tot else "none"
    a, b = c["PTD_e10_e2"][1], c["normal_e10_e11"][1]
    share = "%.3f" % (a / (a + b)) if a + b else "NA"
    summary.append((s, reads, tot, top, ptd.get(top, 0), a, b, share, c["normal_e1_e2"][1]))
with open(os.path.join(out, "fastq_junction_counts.tsv"), "w") as fo:
    fo.write("sample\tprobe\ttype\tdonor_exon\tacceptor_exon\tfragments\tdistinct_sequences\n")
    fo.writelines("\t".join(map(str, r)) + "\n" for r in rows)
M = {(r[0], r[1]): r[5] for r in rows}
with open(os.path.join(out, "fastq_junction_matrix.tsv"), "w") as fo:
    fo.write("probe\t" + "\t".join(samples) + "\n")
    for p in [x[0] for x in P]:
        v = [M[(s, p)] for s in samples]
        if p.startswith("normal") or any(v):
            fo.write(p + "\t" + "\t".join(map(str, v)) + "\n")
with open(os.path.join(out, "fastq_ptd_summary.tsv"), "w") as fo:
    fo.write("sample\treads_searched\tPTD_fragments_all_probes\ttop_PTD_junction\ttop_PTD_fragments\t"
             "PTD_e10_e2\tnormal_e10_e11\tPTD_e10_e2_share\tnormal_e1_e2\n")
    fo.writelines("\t".join(map(str, r)) + "\n" for r in summary)
for fn in ("fastq_ptd_summary.tsv", "fastq_junction_matrix.tsv"):
    print("--- %s ---" % fn)
    print(open(os.path.join(out, fn)).read())
PYEOF

echo
echo "############ SUMMARY ############"
"$PY" "$HELPER" "$PTD_DIR/probes.tsv" "$OUT" "${SAMPLES[@]}"

{
  echo "run          : $(date '+%F %T %Z') on $(hostname)"
  echo "script       : $0"
  echo "FASTQs       : $( [[ -n "$FASTQ_PATTERN" ]] && echo "$FASTQ_PATTERN" || echo "STAR --readFilesIn from each BAM header")"
  for i in "${!JOB_FILE[@]}"; do echo "               ${JOB_SAMPLE[$i]} file ${JOB_MATE[$i]}: ${JOB_FILE[$i]}"; done
  echo "probes       : $PTD_DIR/probes.tsv (66 PTD + 12 normal junctions, 15 + 15 bases, both orientations, exact match)"
} > "$OUT/RUN_INFO.txt"
cat "$OUT/RUN_INFO.txt"
echo
echo "Done. Copy $OUT to the analysis computer as rnaseq/09_kmt2a_ptd/fastq_search/"
