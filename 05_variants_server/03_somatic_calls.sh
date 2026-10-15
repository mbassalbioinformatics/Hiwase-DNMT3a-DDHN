#!/usr/bin/env bash
###############################################################################
# 03_somatic_calls.sh
#
# Genome-wide somatic SNVs and small indels with Mutect2 on the methylation-
# masked BAMs from 02, one multi-sample run per case.
#
# Design (config/samples.tsv, column variant_role): in each case the libraries
# without the clone are the normals (the first one listed is the baseline), and
# the other libraries of the same person are the tumours, called jointly so that
# the clone's variants are followed across the timeline. The 1000 Genomes panel
# of normals and gnomAD filter artefacts and germline variants.
# A second normal of the same material as the tumours is worth having: with the
# donor's graft (progenitor cells collected from blood) as the only normal,
# hundreds of calls per chromosome that every marrow or blood mononuclear-cell
# library shared at ~20% were differences of material, not of clone; a day-30
# marrow without the clone cancels them.
# Contamination of every library other than the baseline is measured against
# the baseline (another person's DNA; after a transplant, the recipient's).
#
# Phases (each resumable; finished outputs are skipped, just re-run):
#   A  Mutect2 per (case x 25-Mb chunk), bounded parallel pool
#   B  per case: MergeVcfs, MergeMutectStats, LearnReadOrientationModel
#   C  per library: GetPileupSummaries; per library other than the baseline:
#      CalculateContamination vs the baseline (only the tumours' tables go to D)
#   D  per case: FilterMutectCalls
# Mapping quality >= 10 everywhere (MMQ). Msuite2's aligner gives at most 42,
# and the GATK defaults of 50 drop every read: empty contamination tables and no
# data for the read-orientation model.
# Each step is redone when one of its inputs is newer than its output. A change
# of a case's design (normals, tumours, MMQ, resources; kept in
# mutect2/<case>.design) removes that case's earlier outputs first.
# Output: $VAR_OUT/mutect2/vcf/<case>.filtered.vcf.gz, <library>.contamination.table
#
# Run on the server after 02:
#   bash 05_variants_server/03_somatic_calls.sh --test chr21    # rehearsal on the 02 test BAMs
#   bash 05_variants_server/03_somatic_calls.sh
###############################################################################
set -uo pipefail   # not -e: one failed shard must not stop a long run
source "$(dirname "${BASH_SOURCE[0]}")/sheet.sh"     # config/paths.sh and the sample sheet
OUT=${OUT:-$VAR_OUT}
DICT="${REF%.fa}.dict"
PON="$MUTECT2_PON"; GNOMAD="$GERMLINE_RESOURCE"; COMMON="$COMMON_SITES"
CHUNK_MB=${CHUNK_MB:-25}
MAX_PARALLEL=${MAX_PARALLEL:-12}   # Mutect2 shards at a time
HMM_THREADS=${HMM_THREADS:-3}      # pair-HMM threads per shard: 12 x 3 = 36 threads, for about 32 free cores
# Minimum mapping quality for every cut-off below. Msuite2's aligner (bowtie2) gives at most 42, and
# GATK's defaults sit near or above that: Mutect2 read filter 20, Mutect2's F1R2 (read-orientation)
# collection 50, GetPileupSummaries 50, FilterMutectCalls median of the alt reads 30.
MMQ=${MMQ:-10}
JAVA_MEM="-Xmx8g"
# cases, normals (the first one is the baseline) and tumours: sample sheet, column variant_role (sheet.sh)
declare -A NORMALS TUMOURS
for c in "${CASES[@]}"; do NORMALS[$c]="${M2_NORMALS[$c]:-}"; TUMOURS[$c]="${M2_TUMOURS[$c]:-}"; done

TEST=""
if [[ "${1:-}" == "--test" ]]; then TEST="${2:?--test needs a contig, e.g. chr21}"; OUT="$OUT/test_$TEST"; fi
BAMS="$OUT/bams"; M="$OUT/mutect2"; SHARD="$M/shards"; VCF="$M/vcf"; LOG="$M/logs"
mkdir -p "$SHARD" "$VCF" "$LOG"
exec 9> "$OUT/.lock.03"; flock -n 9 || { echo "ERROR: 03 is already running on $OUT"; exit 3; }
rm -f "$OUT/.done.03"     # written again only when both cases are filtered
DOCKER="docker run --platform linux/amd64 --rm $DOCKER_MOUNTS $GATK_IMAGE"
GATK="$DOCKER gatk --java-options $JAVA_MEM"
for f in "$REF" "$DICT" "$PON" "$GNOMAD" "$COMMON"; do
  [[ -s "$f" ]] || { echo "ERROR: missing $f (set it in config/paths.sh)"; exit 1; }
done
for c in "${CASES[@]}"; do for s in ${NORMALS[$c]} ${TUMOURS[$c]}; do
  [[ -s "$BAMS/$s.taps.bam" && -s "$BAMS/$s.taps.bam.bai" ]] || { echo "ERROR: $BAMS/$s.taps.bam missing (run 02 first)"; exit 1; }
done; done

# --- design of each case: when it changed, the case's earlier outputs are removed ---
for c in "${CASES[@]}"; do
  d="normals=${NORMALS[$c]}; tumours=${TUMOURS[$c]}; mmq=$MMQ; pon=${PON##*/}; germline=${GNOMAD##*/}; common=${COMMON##*/}"
  st="$M/$c.design"
  if [[ "$(cat "$st" 2>/dev/null)" != "$d" ]]; then
    if compgen -G "$SHARD/$c.*" > /dev/null || compgen -G "$VCF/$c.*" > /dev/null; then
      echo "$(date '+%F %T') design of $c changed (was: $(cat "$st" 2>/dev/null || echo 'not recorded')); removing its earlier outputs"
      rm -f "$SHARD/$c".* "$VCF/$c".*
      for s in ${NORMALS[$c]} ${TUMOURS[$c]}; do rm -f "$VCF/$s".{pileups,contamination,segments}.table; done
    fi
    echo "$d" > "$st"
  fi
  echo "$(date '+%F %T') $c: $d"
done

# --- chunks in reference dictionary order --------------------------------------
mapfile -t CHUNKS < <(awk -v step=$((CHUNK_MB * 1000000)) -v only="$TEST" '
  /^@SQ/ { n=""; l=0; for (i=1;i<=NF;i++) { if ($i~/^SN:/) n=substr($i,4); if ($i~/^LN:/) l=substr($i,4)+0 }
           if (n=="" || l==0 || (only!="" && n!=only)) next
           for (s=1; s<=l; s+=step) { e=s+step-1; if (e>l) e=l; printf "%s:%d-%d\n", n, s, e } }' "$DICT")
echo "$(date '+%F %T') ${#CASES[@]} cases x ${#CHUNKS[@]} chunks"

fresh() {  # fresh <output> <input...>: the output exists and is newer than every input
  local o="$1" i; shift; [[ -s "$o" ]] || return 1
  for i in "$@"; do [[ "$o" -nt "$i" ]] || return 1; done
}

# ============================ A: scattered Mutect2 ==============================
run_shard() {   # redone when a BAM of the case is newer than the shard (02 rebuilt it)
  local c="$1" region="$2" tag out args=() bams=()
  tag="${region//[:-]/_}"; out="$SHARD/$c.$tag"
  for s in ${NORMALS[$c]} ${TUMOURS[$c]}; do args+=(-I "$BAMS/$s.taps.bam"); bams+=("$BAMS/$s.taps.bam"); done
  for s in ${NORMALS[$c]}; do args+=(-normal "$s"); done
  fresh "$out.vcf.gz" "${bams[@]}" && [[ -s "$out.vcf.gz.stats" && -s "$out.f1r2.tar.gz" ]] && return 0
  if $GATK Mutect2 -R "$REF" "${args[@]}" -L "$region" \
       -pon "$PON" --germline-resource "$GNOMAD" --minimum-mapping-quality "$MMQ" --native-pair-hmm-threads "$HMM_THREADS" \
       --f1r2-median-mq "$MMQ" --f1r2-tar-gz "$out.f1r2.tar.gz" \
       -O "$out.vcf.gz" > "$LOG/mutect2.$c.$tag.log" 2>&1 \
     && $DOCKER bgzip -t "$out.vcf.gz" 2>/dev/null; then
    echo "OK   $c $region"
  else
    echo "FAIL $c $region (see $LOG/mutect2.$c.$tag.log)"; echo "$c $region" >> "$M/.failures.txt"
    rm -f "$out.vcf.gz" "$out.vcf.gz.tbi" "$out.vcf.gz.stats" "$out.f1r2.tar.gz"
  fi
  return 0
}
echo "$(date '+%F %T') === A: Mutect2 ($MAX_PARALLEL at a time)"
rm -f "$M/.failures.txt"; running=0
for c in "${CASES[@]}"; do for region in "${CHUNKS[@]}"; do
  run_shard "$c" "$region" &
  running=$((running + 1))
  if (( running >= MAX_PARALLEL )); then wait -n || true; running=$((running - 1)); fi
done; done
wait
if [[ -f "$M/.failures.txt" ]]; then
  echo "$(wc -l < "$M/.failures.txt") shard(s) failed; re-run to retry them:"; cat "$M/.failures.txt"; exit 1
fi

# ======================= B-D: per case gather and filter =========================
step() {   # step <label> <output removed on failure> <command...>; tool output in $LOG/<label>.log
  local label="$1" out="$2" lg; shift 2
  lg="$LOG/${label// /.}.log"; echo "  $label"
  "$@" > "$lg" 2>&1 && return 0
  echo "  FAIL $label (see $lg)"; rm -f "$out" "$out.tbi" "$out.idx"; return 1
}
sites() { local n; n=$(grep -vc -e '^#' -e '^contig' "$1" 2>/dev/null); echo "${n:-0}"; }   # data rows of a pileup table
pileup() { # pileup <sample> <table>: allele counts at common SNPs; an empty table counts as a failure
  local s="$1" t="$2" n
  fresh "$t" "$BAMS/$s.taps.bam" && [[ $(sites "$t") -gt 0 ]] && return 0
  step "GetPileupSummaries $s" "$t" $GATK GetPileupSummaries -I "$BAMS/$s.taps.bam" -V "$COMMON" -L "$COMMON" "${pil[@]}" \
       --minimum-mapping-quality "$MMQ" -O "$t" || return 1
  n=$(sites "$t"); echo "  $s: $n common SNP sites covered"
  (( n > 0 )) || { echo "  FAIL GetPileupSummaries $s: no site covered"; rm -f "$t"; return 1; }
}
FAILED=()
for c in "${CASES[@]}"; do
  echo "$(date '+%F %T') === $c"
  RAW="$VCF/$c.raw.vcf.gz"; STATS="$VCF/$c.merged.stats"; ORIENT="$VCF/$c.read_orientation_model.tar.gz"; FINAL="$VCF/$c.filtered.vcf.gz"
  av=(); as=(); af=(); sv=(); ss=(); sf=()
  for region in "${CHUNKS[@]}"; do
    tag="${region//[:-]/_}"
    av+=(-I "$SHARD/$c.$tag.vcf.gz"); as+=(-stats "$SHARD/$c.$tag.vcf.gz.stats"); af+=(-I "$SHARD/$c.$tag.f1r2.tar.gz")
    sv+=("$SHARD/$c.$tag.vcf.gz"); ss+=("$SHARD/$c.$tag.vcf.gz.stats"); sf+=("$SHARD/$c.$tag.f1r2.tar.gz")
  done
  ok=1
  if [[ $ok == 1 ]] && ! fresh "$RAW" "${sv[@]}"; then step "MergeVcfs $c" "$RAW" $GATK MergeVcfs "${av[@]}" -O "$RAW" -D "$DICT" || ok=0; fi
  # MergeVcfs keeps a single ##normal_sample and ##tumor_sample line, while FilterMutectCalls takes every
  # sample not listed as a normal for a tumour; the full lists are written back from the design
  if [[ $ok == 1 ]]; then
    roles="$VCF/$c.sample_roles.txt"
    { for s in ${NORMALS[$c]}; do echo "##normal_sample=$s"; done; for s in ${TUMOURS[$c]}; do echo "##tumor_sample=$s"; done; } > "$roles"
    if [[ "$($DOCKER bcftools view -h "$RAW" 2>/dev/null | grep -E '^##(normal|tumor)_sample=' | sort)" != "$(sort "$roles")" ]]; then
      echo "  sample roles written back into the header of $c.raw.vcf.gz ($(grep -c normal_sample "$roles") normal, $(grep -c tumor_sample "$roles") tumour)"
      $DOCKER bcftools view -h "$RAW" | grep -vE '^##(normal|tumor)_sample=' \
        | awk -v f="$roles" '/^#CHROM/ { while ((getline l < f) > 0) print l } { print }' > "$VCF/$c.raw.header.txt" \
        && $DOCKER bash -c "bcftools reheader -h $VCF/$c.raw.header.txt -o $RAW.tmp $RAW && mv $RAW.tmp $RAW && bcftools index -f -t $RAW" \
        || { echo "  FAIL sample roles $c"; ok=0; }
    fi
  fi
  if [[ $ok == 1 ]] && ! fresh "$STATS" "${ss[@]}"; then step "MergeMutectStats $c" "$STATS" $GATK MergeMutectStats "${as[@]}" -O "$STATS" || ok=0; fi
  if [[ $ok == 1 ]] && ! fresh "$ORIENT" "${sf[@]}"; then
    step "LearnReadOrientationModel $c" "$ORIENT" $GATK LearnReadOrientationModel "${af[@]}" -O "$ORIENT" || ok=0
    # every context skipped = no F1R2 data, so the orientation filter would do nothing
    skipped=$(grep -c "Skipping the reference context" "$LOG/LearnReadOrientationModel.$c.log" 2>/dev/null)
    [[ $ok == 1 ]] && echo "  $c read-orientation model: ${skipped:-0} contexts skipped for lack of data"
  fi

  # contamination of every library against the baseline normal (redone, with
  # everything after it, when an input is newer than the output); FilterMutectCalls
  # takes the tumours' tables
  read -r -a nrm <<< "${NORMALS[$c]}"; n=${nrm[0]}; NP="$VCF/$n.pileups.table"
  pil=(); [[ -n "$TEST" ]] && pil=(-L "$TEST" --interval-set-rule INTERSECTION)
  [[ $ok == 1 ]] && { pileup "$n" "$NP" || ok=0; }
  fc=(); cts=()
  for t in "${nrm[@]:1}" ${TUMOURS[$c]}; do
    TP="$VCF/$t.pileups.table"; CT="$VCF/$t.contamination.table"; SG="$VCF/$t.segments.table"
    [[ $ok == 1 ]] && { pileup "$t" "$TP" || ok=0; }
    if [[ $ok == 1 ]] && ! fresh "$CT" "$TP" "$NP"; then
      step "CalculateContamination $t" "$CT" $GATK CalculateContamination -I "$TP" -matched "$NP" --tumor-segmentation "$SG" -O "$CT" || ok=0
    fi
    [[ $ok == 1 ]] && echo "  $t contamination vs $n (fraction, error): $(awk '!/^#/ && $1 != "sample" {print $2", "$3; exit}' "$CT")"
    [[ " ${TUMOURS[$c]} " == *" $t "* ]] && { fc+=(--contamination-table "$CT" --tumor-segmentation "$SG"); cts+=("$CT"); }
  done

  if [[ $ok == 1 ]] && ! fresh "$FINAL" "$RAW" "$STATS" "$ORIENT" "${cts[@]}"; then
    step "FilterMutectCalls $c" "$FINAL" $GATK FilterMutectCalls -R "$REF" -V "$RAW" --stats "$STATS" \
         --ob-priors "$ORIENT" "${fc[@]}" --min-median-mapping-quality "$MMQ" -O "$FINAL" || ok=0
  fi
  if [[ $ok == 1 && -s "$FINAL" ]]; then
    echo "OK   $c -> $FINAL ($($DOCKER bcftools view -f PASS -H "$FINAL" 2>/dev/null | wc -l) PASS)"
  else
    FAILED+=("$c")
  fi
done
(( ${#FAILED[@]} )) && { echo "FAILED: ${FAILED[*]} (re-run to retry)"; exit 1; }
date '+%F %T' > "$OUT/.done.03"
echo "$(date '+%F %T') DONE Mutect2. Next: 04_cnv_loh.sh, 05_annotate_collect.sh"
