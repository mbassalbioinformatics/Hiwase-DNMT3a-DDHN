#!/usr/bin/env bash
###############################################################################
# 04_kmt2a_ptd_junctions.sh  --  RNA-seq check for a KMT2A partial tandem
#                                duplication (KMT2A-PTD), on the server
#
# From the repository root, on the server that holds the RNA-seq alignments:
#   bash 03_rnaseq/server/04_kmt2a_ptd_junctions.sh 2>&1 | tee kmt2a_ptd_run.log
#
# Then copy $RNA_OUT/09_kmt2a_ptd/ to the analysis computer as rnaseq/09_kmt2a_ptd/
#
# WHY
#   The methylation libraries, read as whole-genome sequencing, show a 23-kb
#   block inside KMT2A (exons 2-10 of the MANE Select transcript
#   ENST00000534358) at about twice the normal copy number at 3.5 years and at
#   DDHN diagnosis, and at the normal level at 1 month and in the donor 5 years
#   after donation. A tandem duplication of exons 2-10 makes an mRNA in which
#   the end of exon 10 is spliced back onto the start of exon 2, so reads that
#   cross that junction in the RNA-seq are an independent check of the duplication.
#
# HOW
#   1. Probes. For every donor exon i (2-12) joined back to every acceptor exon
#      j <= i (66 possible PTD junctions; the copy-number breakpoints are only
#      known to 1 kb, so all of them are tested), and for the 12 normal
#      junctions exon i -> i+1 (an expression control), a probe of the last
#      15 bases of exon i and the first 15 bases of exon j. The exon ends are
#      embedded below (MANE v1.5, GRCh38, checked against UCSC hg38 on
#      19 Oct 2026); with REF_FASTA set they are re-checked against the FASTA
#      the reads were aligned to. Every probe was checked to occur nowhere in
#      the KMT2A gene or its normal mRNA, on either strand.
#   2. Reads. Primary alignments in the KMT2A locus (both mates), plus the
#      unmapped reads if the BAM keeps them. A read crossing a PTD junction
#      cannot be aligned in one piece, so STAR soft-clips it at the exon end,
#      splits it, or leaves it unmapped; matching the read sequence finds it
#      whatever the aligner did. FULL_SCAN=1 searches every read in the BAM
#      (slow; catches reads placed elsewhere in the genome).
#   3. Matching. Each probe and its reverse complement is searched for exactly
#      (grep -F, then re-checked in the read-sequence field). A fragment counts
#      once per junction even if both mates match.
#
# EXPECTED IF THE DUPLICATION IS EXONS 2-10
#   Reads for PTD_e10_e2 at 3.5 years and DDHN diagnosis, none at 1 month or in the donor; the
#   normal junctions in all four. A different PTD junction in the clone
#   libraries only would mean a duplication with slightly different ends.
#
# OUTPUT (in OUT_DIR)
#   ptd_summary.tsv        per library: PTD reads, the top PTD junction, PTD_e10_e2
#                          against the normal exon 10 -> 11 junction, region depth
#   junction_matrix.tsv    every probe with any read, by library
#   junction_counts.tsv    every probe, every library (fragments, distinct sequences)
#   <sample>/supporting_reads.tsv   each matching read: probe, mate, flag, position,
#                                   CIGAR, where it was found, sequence
#   <sample>/supporting_reads.bam   all alignments of those reads (samtools >= 1.12)
#   probes.tsv, exons.tsv, RUN_INFO.txt
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
NTHREAD="${THREADS:-16}"
OUT_DIR="$RNA_OUT/09_kmt2a_ptd"

# The genome FASTA the reads were aligned to (Cell Ranger GRCh38-3.0.0), used only to re-check
# the embedded exon sequences and probes; without it the run relies on the embedded sequences.
REF_FASTA="${RNA_REF_FASTA:-}"

OVERHANG=15        # bases each side of the junction a read must cover (12-30)
FULL_SCAN=0        # 1 = also search every read in each BAM (slow)
PY="python3"       # any Python 3; standard library only
# ============================================================================

# KMT2A exons 1-13, MANE Select ENST00000534358.8 (GRCh38, + strand, 1-based):
# exon, start, end, first 30 bases, last 30 bases.
EXONS_TSV=$(cat <<'EOF'
1	118436492	118436944	GCTTCACTTCACGGGGCGAACATGGCGCAC	GAGAGCGGCGGGGGAGGCGGCAGCGGAGAG
2	118468775	118468844	GATGAGCAATTCTTAGGTTTTGGCTCAGAT	GAGTGCGAAGTCCCACAAGGTCTCCTTCAG
3	118471662	118474315	TTAAAACTAGTCCTCGAAAACCTCGTGGGA	CTTAAACAAACCGACCAGCCCAAAGCACAG
4	118476805	118476982	GGTCAAGAAAGTGACTCATCAGAGACCTCT	AAAAGATTTTGTCTTCCATGGGGAATGATG
5	118477967	118478201	ACAAGTCATCAATTGCTGGCTCAGAAGATG	TGGTCGCAATATAAAGAAGCAGTGCTGCAA
6	118480174	118480238	GATGAGAAAATGTCAGAATCTACAATGGAT	CCAAAGCCTACCTGCAGAAGCAAGCTAAAG
7	118481715	118482092	CTGTGAAAAAGAAAGAGAAAAAGTCTAAGA	AGAAAAAGCAGCCTCCACCACCAGAATCAG
8	118482422	118482495	GTCCAGAGCAGAGCAAACAGAAAAAAGTGG	ATCCCTGTAAAACAAAAACCAAAAGAAAAG
9	118484183	118484314	GAAAAACCACCTCCGGTCAATAAGCAGGAG	GGAGTCCACAGGATCAGAGTGGACTTTAAG
10	118484862	118484975	GAGGATTGTGAAGCAGAAAATGTGTGGGAG	TTTCTCTGTGCCAGTAGTGGGCATGTAGAG
11	118488614	118488760	TTTGTGTATTGCCAAGTCTGTTGTGAGCCC	GTTTGTGGAAGGCAACATCAGGCTACAAAG
12	118489792	118489887	CAGCTGCTGGAGTGTAATAAGTGCCGAAAC	ACCAAACCCACAAAGAAGAAGAAAGTCTGG
13	118490129	118490249	ATCTGTACCAAGTGTGTTCGCTGTAAGAGC	GTCATGATTGCGCCAAGCTCTTTGCTAAAG
EOF
)
LOCUS_START=118436000   # the KMT2A locus searched (GRCh38), exon 1 to beyond exon 13
LOCUS_END=118527000
GRCH38_CHR11_LEN=135086622

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "############ PREFLIGHT ############"
if [[ -n "$CONDA_ENV" ]] && command -v conda >/dev/null 2>&1; then
  # conda's own scripts (and the activate.d hooks of packages such as MKL) read
  # variables that may not be set, so nounset is off while conda runs
  set +u
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV" || echo "WARNING: could not activate '$CONDA_ENV'; using PATH"
  set -u
fi
command -v samtools >/dev/null 2>&1 || { echo "FATAL: samtools not on PATH"; exit 1; }
command -v "$PY" >/dev/null 2>&1 || { echo "FATAL: $PY not on PATH"; exit 1; }
SAMTOOLS_VERSION=$(samtools --version | head -1)
echo "samtools     : $(command -v samtools) ($SAMTOOLS_VERSION)"
echo "python       : $(command -v "$PY") ($("$PY" --version 2>&1))"
[[ "$OVERHANG" -ge 12 && "$OVERHANG" -le 30 ]] || { echo "FATAL: OVERHANG must be 12-30"; exit 1; }

for s in "${SAMPLES[@]}"; do
  b="${BAM_DIR}/${s}${BAM_SUFFIX}"
  [[ -f "$b" ]] || { echo "FATAL: BAM not found: $b"; exit 1; }
  [[ -f "${b}.bai" || -f "${b%.bam}.bai" ]] || { echo "  indexing $s"; samtools index -@ "$NTHREAD" "$b"; }
done
echo "BAM dir      : $BAM_DIR"

# chromosome naming (11 or chr11) and assembly, from the first BAM's header
B0="${BAM_DIR}/${SAMPLES[0]}${BAM_SUFFIX}"
HDR=$(samtools view -H "$B0")
sq_len() {   # sq_len <header text> <sequence name>: its length from the @SQ line, empty if absent
  awk -F'\t' -v n="SN:$2" '$1 == "@SQ" { has = 0; len = ""; for (i = 2; i <= NF; i++) { if ($i == n) has = 1; if ($i ~ /^LN:/) len = substr($i, 4) }
                                         if (has) { print len; exit } }' <<<"$1"
}
if   [[ -n "$(sq_len "$HDR" chr11)" ]]; then CHR="chr11"
elif [[ -n "$(sq_len "$HDR" 11)" ]];    then CHR="11"
else echo "FATAL: no chromosome 11 (11 or chr11) in the BAM header"; exit 1; fi
LEN=$(sq_len "$HDR" "$CHR")
[[ "$LEN" == "$GRCH38_CHR11_LEN" ]] || { echo "FATAL: ${CHR} length $LEN is not GRCh38 ($GRCH38_CHR11_LEN)"; exit 1; }
REGION="${CHR}:${LOCUS_START}-${LOCUS_END}"
echo "chromosome   : $CHR (GRCh38, length $LEN); locus searched $REGION"
for s in "${SAMPLES[@]}"; do
  h=$(sq_len "$(samtools view -H "${BAM_DIR}/${s}${BAM_SUFFIX}")" "$CHR")
  [[ "$h" == "$LEN" ]] || { echo "FATAL: $s has a different ${CHR} ($h)"; exit 1; }
done

echo "--- aligner (@PG) ---"
awk -F'\t' '$1 == "@PG"' <<<"$HDR" | head -2 | cut -c1-200

# ---------------------------------------------------------------------------
# helper (Python, standard library only): probes, checks, counting, summary
# ---------------------------------------------------------------------------
HELPER="$OUT_DIR/kmt2a_ptd_helper.py"
cat > "$HELPER" <<'PYEOF'
"""Helper for 04_kmt2a_ptd_junctions.sh (standard library only)."""
import collections
import os
import sys

RC = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(RC)[::-1]


def read_exons(path):
    ex = {}
    for line in open(path):
        if not line.strip() or line.startswith("#"):
            continue
        e, a, b, first, last = line.rstrip("\n").split("\t")
        ex[int(e)] = {"start": int(a), "end": int(b), "first": first, "last": last}
    return ex


def read_probes(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path)]
    return rows[1:]


def probes(exons_tsv, overhang, out_tsv, out_pat):
    a = int(overhang)
    ex = read_exons(exons_tsv)
    rows = []
    for i in range(2, 13):                      # PTD: exon i spliced back onto exon j <= i
        for j in range(2, i + 1):
            rows.append(("PTD_e%d_e%d" % (i, j), "PTD", i, j, ex[i]["last"][-a:] + ex[j]["first"][:a]))
    for i in range(1, 13):                      # normal splicing, exon i -> i+1
        rows.append(("normal_e%d_e%d" % (i, i + 1), "normal", i, i + 1, ex[i]["last"][-a:] + ex[i + 1]["first"][:a]))
    seen = {}
    for name, _, _, _, k in rows:               # no two probes may share a sequence, on either strand
        for x in (k, rc(k)):
            if x in seen and seen[x] != name:
                sys.exit("FATAL: probes %s and %s share a sequence" % (name, seen[x]))
            seen[x] = name
    with open(out_tsv, "w") as f:
        f.write("probe\ttype\tdonor_exon\tacceptor_exon\tsequence\treverse_complement\n")
        for name, t, i, j, k in rows:
            f.write("%s\t%s\t%d\t%d\t%s\t%s\n" % (name, t, i, j, k, rc(k)))
    with open(out_pat, "w") as f:
        for _, _, _, _, k in rows:
            f.write(k + "\n" + rc(k) + "\n")
    print("probes       : %d PTD + %d normal junctions, %d + %d bases"
          % (sum(r[1] == "PTD" for r in rows), sum(r[1] == "normal" for r in rows), a, a))


def verify(exons_tsv, probes_tsv, region_fa, region_start):
    """Re-check the embedded exon ends and the probes against the alignment reference."""
    seq = "".join(l.strip() for l in open(region_fa) if not l.startswith(">")).upper()
    off = int(region_start)
    g = lambda a, b: seq[a - off: b - off + 1]
    bad = 0
    for e, x in sorted(read_exons(exons_tsv).items()):
        ok = g(x["start"], x["start"] + 29) == x["first"] and g(x["end"] - 29, x["end"]) == x["last"]
        don = g(x["end"] + 1, x["end"] + 2)
        acc = g(x["start"] - 2, x["start"] - 1) if e > 1 else "AG"
        if not ok or don != "GT" or acc != "AG":
            bad += 1
            print("  MISMATCH exon %d (sequence %s, acceptor %s, donor %s)" % (e, "ok" if ok else "differs", acc, don))
    for p in read_probes(probes_tsv):
        if p[1] == "PTD" and (p[4] in seq or p[5] in seq):
            bad += 1
            print("  probe %s occurs in the KMT2A gene itself" % p[0])
    if bad:
        sys.exit("FATAL: %d checks failed against the reference FASTA" % bad)
    print("reference    : embedded exon ends, splice sites and probes agree with the FASTA")


def count(probes_tsv, sample, outdir, *hits):
    P = read_probes(probes_tsv)
    support = collections.defaultdict(dict)     # probe -> read name -> sequence
    recs = []
    for h in hits:
        if not os.path.exists(h):
            continue
        src = os.path.basename(h)[len("hits_"):-len(".sam")]
        for line in open(h):
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            qn, flag, chrom, pos, cigar, seq = f[0], int(f[1]), f[2], f[3], f[5], f[9].upper()
            for name, t, i, j, k, krc in P:
                if k in seq or krc in seq:      # the match must be in the read sequence itself
                    mate = "R1" if flag & 0x40 else ("R2" if flag & 0x80 else "SE")
                    support[name].setdefault(qn, seq)
                    recs.append((sample, name, t, qn, mate, flag, chrom, pos, cigar, src, seq))
    with open(os.path.join(outdir, "supporting_reads.tsv"), "w") as f:
        f.write("sample\tprobe\ttype\tread\tmate\tflag\tchrom\tpos\tcigar\tfound_in\tsequence\n")
        for r in recs:
            f.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(outdir, "junction_counts.tsv"), "w") as f:
        f.write("sample\tprobe\ttype\tdonor_exon\tacceptor_exon\tfragments\tdistinct_sequences\n")
        for name, t, i, j, k, krc in P:
            q = support.get(name, {})
            f.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\n" % (sample, name, t, i, j, len(q), len(set(q.values()))))
    names = sorted({r[3] for r in recs if r[2] == "PTD"})
    with open(os.path.join(outdir, "ptd_read_names.txt"), "w") as f:
        f.write("".join(n + "\n" for n in names))
    ptd = sum(len(v) for k, v in support.items() if k.startswith("PTD"))
    print("  %s: %d fragments on PTD probes, %d on normal probes"
          % (sample, ptd, sum(len(v) for k, v in support.items() if k.startswith("normal"))))


def summary(probes_tsv, outroot, *samples):
    rows = []
    for s in samples:
        rows += [l.rstrip("\n").split("\t") for l in list(open(os.path.join(outroot, s, "junction_counts.tsv")))[1:]]
    with open(os.path.join(outroot, "junction_counts.tsv"), "w") as f:
        f.write("sample\tprobe\ttype\tdonor_exon\tacceptor_exon\tfragments\tdistinct_sequences\n")
        f.writelines("\t".join(r) + "\n" for r in rows)
    order = [p[0] for p in read_probes(probes_tsv)]
    M = {(r[0], r[1]): r[5] for r in rows}
    with open(os.path.join(outroot, "junction_matrix.tsv"), "w") as f:
        f.write("probe\t" + "\t".join(samples) + "\n")
        for p in order:
            v = [M.get((s, p), "0") for s in samples]
            if p.startswith("normal") or any(x != "0" for x in v):
                f.write(p + "\t" + "\t".join(v) + "\n")
    with open(os.path.join(outroot, "ptd_summary.tsv"), "w") as f:
        f.write("sample\tlocus_reads\tPTD_fragments_all_probes\ttop_PTD_junction\ttop_PTD_fragments\t"
                "PTD_e10_e2\tnormal_e10_e11\tPTD_e10_e2_share\tnormal_e1_e2\n")
        for s in samples:
            r = {x[1]: (int(x[5]), x[2]) for x in rows if x[0] == s}
            ptd = {k: v[0] for k, v in r.items() if v[1] == "PTD"}
            tot = sum(ptd.values())
            top = max(ptd, key=ptd.get) if tot else "none"
            a, b = r["PTD_e10_e2"][0], r["normal_e10_e11"][0]
            share = "%.3f" % (a / (a + b)) if a + b else "NA"
            depth = open(os.path.join(outroot, s, "locus_reads.txt")).read().strip()
            f.write("%s\t%s\t%d\t%s\t%d\t%d\t%d\t%s\t%d\n"
                    % (s, depth, tot, top, ptd.get(top, 0), a, b, share, r["normal_e1_e2"][0]))
    for fn in ("ptd_summary.tsv", "junction_matrix.tsv"):
        print("--- %s ---" % fn)
        print(open(os.path.join(outroot, fn)).read())


if __name__ == "__main__":
    {"probes": probes, "verify": verify, "count": count, "summary": summary}[sys.argv[1]](*sys.argv[2:])
PYEOF

printf '%s\n' "$EXONS_TSV" > "$OUT_DIR/exons.tsv"
"$PY" "$HELPER" probes "$OUT_DIR/exons.tsv" "$OVERHANG" "$OUT_DIR/probes.tsv" "$OUT_DIR/patterns.txt"

# re-check the embedded sequences against a GRCh38 FASTA. Optional: the FASTA may
# name the chromosome 11 or chr11 whatever the BAM uses (the sequence is the same),
# and if it cannot be used the run carries on with the embedded sequences.
REF_CHECK="not re-checked"
if [[ -f "$REF_FASTA" ]]; then
  [[ -f "${REF_FASTA}.fai" ]] || samtools faidx "$REF_FASTA" 2>/dev/null || true
  FCHR=""
  [[ -f "${REF_FASTA}.fai" ]] && FCHR=$(awk -v L="$GRCH38_CHR11_LEN" '($1 == "chr11" || $1 == "11") && $2 == L { print $1; exit }' "${REF_FASTA}.fai")
  if [[ -n "$FCHR" ]] && samtools faidx "$REF_FASTA" "${FCHR}:${LOCUS_START}-${LOCUS_END}" > "$OUT_DIR/kmt2a_locus.fa" 2>/dev/null \
     && [[ $(grep -vc '^>' "$OUT_DIR/kmt2a_locus.fa") -gt 0 ]]; then
    "$PY" "$HELPER" verify "$OUT_DIR/exons.tsv" "$OUT_DIR/probes.tsv" "$OUT_DIR/kmt2a_locus.fa" "$LOCUS_START"
    REF_CHECK="checked against $REF_FASTA (chromosome $FCHR)"
  else
    echo "WARNING: no usable GRCh38 chromosome 11 (chr11 or 11, length $GRCH38_CHR11_LEN) in $REF_FASTA;"
    echo "         using the embedded sequences (checked against UCSC hg38)"
    REF_CHECK="not re-checked (no usable chromosome 11 in $REF_FASTA)"
  fi
else
  echo "WARNING: REF_FASTA not found ($REF_FASTA); using the embedded sequences (checked against UCSC hg38)"
  REF_CHECK="not re-checked (REF_FASTA not found)"
fi

# does samtools support selecting reads by name (-N, samtools >= 1.12)?
VIEW_HELP=$(samtools view --help 2>&1 < /dev/null || true)
if grep -q 'qname-file' <<<"$VIEW_HELP"; then HAVE_N=1; else HAVE_N=0; fi

echo
echo "############ SEARCH ############"
search() {   # search <label> <output> <samtools view arguments...>
  local label="$1" out="$2"; shift 2
  set +e
  samtools view "$@" | grep -F -f "$OUT_DIR/patterns.txt" > "$out"
  local st=("${PIPESTATUS[@]}")
  set -e
  [[ "${st[0]}" -eq 0 ]] || { echo "FATAL: samtools view failed ($label)"; exit 1; }
  [[ "${st[1]}" -le 1 ]] || { echo "FATAL: grep failed ($label)"; exit 1; }
}
for s in "${SAMPLES[@]}"; do
  b="${BAM_DIR}/${s}${BAM_SUFFIX}"
  mkdir -p "$OUT_DIR/$s"
  samtools view -c -F 0x900 "$b" "$REGION" > "$OUT_DIR/$s/locus_reads.txt"
  search "$s locus" "$OUT_DIR/$s/hits_locus.sam" -F 0x900 "$b" "$REGION"
  UNMAPPED=$(samtools idxstats "$b" | awk '$1 == "*" {print $4}')
  if [[ "${UNMAPPED:-0}" -gt 0 ]]; then
    search "$s unmapped" "$OUT_DIR/$s/hits_unmapped.sam" -f 4 -F 0x900 "$b" '*'
  fi
  if [[ "$FULL_SCAN" -eq 1 ]]; then
    search "$s full" "$OUT_DIR/$s/hits_full.sam" -@ "$NTHREAD" -F 0x900 "$b"
  fi
  echo "  $s: $(cat "$OUT_DIR/$s/locus_reads.txt") reads in the locus; ${UNMAPPED:-0} unplaced unmapped reads in the BAM"
  "$PY" "$HELPER" count "$OUT_DIR/probes.tsv" "$s" "$OUT_DIR/$s" "$OUT_DIR/$s"/hits_*.sam
  # every alignment of the PTD-supporting reads, for viewing in IGV
  if [[ "$HAVE_N" -eq 1 && -s "$OUT_DIR/$s/ptd_read_names.txt" ]]; then
    samtools view -b -N "$OUT_DIR/$s/ptd_read_names.txt" -o "$OUT_DIR/$s/supporting_reads.bam" "$b" "$REGION"
    samtools index "$OUT_DIR/$s/supporting_reads.bam"
  fi
done

echo
echo "############ SUMMARY ############"
"$PY" "$HELPER" summary "$OUT_DIR/probes.tsv" "$OUT_DIR" "${SAMPLES[@]}"

{
  echo "run          : $(date '+%F %T %Z') on $(hostname)"
  echo "script       : $0"
  echo "samtools     : $SAMTOOLS_VERSION"
  echo "BAMs         : ${BAM_DIR}/{$(IFS=,; echo "${SAMPLES[*]}")}${BAM_SUFFIX}"
  echo "locus        : $REGION (primary alignments) + unplaced unmapped reads; full scan: $FULL_SCAN"
  echo "probes       : 66 PTD (exon i -> exon j <= i, i = 2-12) + 12 normal junctions; ${OVERHANG} + ${OVERHANG} bases; exact match"
  echo "exons        : MANE Select ENST00000534358.8 (MANE v1.5), GRCh38; $REF_CHECK"
} > "$OUT_DIR/RUN_INFO.txt"
cat "$OUT_DIR/RUN_INFO.txt"
echo
echo "Done. Copy $OUT_DIR to the analysis computer as rnaseq/09_kmt2a_ptd/"
