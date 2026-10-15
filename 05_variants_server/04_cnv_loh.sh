#!/usr/bin/env bash
###############################################################################
# 04_cnv_loh.sh
#
# Genome-wide copy number and loss of heterozygosity (LOH) with the GATK somatic
# CNV workflow, on the methylation-masked BAMs from 02.
#   - read counts in 1-kb bins (CNV_BINS), kept when mappable (Umap k100 >= 0.9),
#     outside segmental duplications (<= 0.5) and of ordinary GC (0.1-0.9)
#   - panel of normals from the least clonal libraries (config/samples.tsv,
#     column cnv_panel), leave-one-out for the libraries in it
#   - allelic counts at common gnomAD SNPs (AF 0.1-0.9); masked bases are ignored
#     (base quality 0), so CpG SNPs are counted on the informative strand. SNPs
#     whose alternative allele creates a CpG are dropped (taps_snp_sites.py):
#     TAPS reads that allele's methylated C as T on one strand
#   - mapping quality >= 10 for read and allele counts (MMQ): the aligner's
#     maximum is 42 and every TAPS conversion counts as a mismatch, so the tools'
#     default of 30 would drop more reads in heavily methylated regions
#   - ModelSegments with the heterozygous sites of the person's first normal
#     library; CallCopyRatioSegments; PlotModeledSegments
# Output: $VAR_OUT/cnv/segments/<sample>.modelFinal.seg, .called.seg; plots in $VAR_OUT/cnv/plots.
# Resumable: finished outputs are skipped; an output older than one of its
# inputs is redone.
#
# Run on the server after 02 (can run alongside 03):
#   bash 05_variants_server/04_cnv_loh.sh --test chr21
#   bash 05_variants_server/04_cnv_loh.sh
###############################################################################
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/sheet.sh"     # config/paths.sh and the sample sheet
OUT=${OUT:-$VAR_OUT}
DICT="${REF%.fa}.dict"
BINS="$CNV_BINS"; UMAP="$UMAP_BED"; SEGDUP="$SEGDUP_BED"
GNOMAD="$GERMLINE_RESOURCE"
JOBS=${JOBS:-7}           # counting jobs at a time
SEG_JOBS=${SEG_JOBS:-4}   # libraries segmented at a time (ModelSegments is memory-hungry genome-wide)
MMQ=${MMQ:-10}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNPF="$HERE/taps_snp_sites.py"
# ALL, the panel of normals (column cnv_panel) and MATCHED (the person's first normal): sample sheet (sheet.sh)
NORMALS=("${CNV_PANEL[@]}")

TEST=""
if [[ "${1:-}" == "--test" ]]; then TEST="${2:?--test needs a contig}"; OUT="$OUT/test_$TEST"; fi
C="$OUT/cnv"; mkdir -p "$C"/{counts,pon,denoised,allelic,segments,plots,logs}
exec 9> "$OUT/.lock.04"; flock -n 9 || { echo "ERROR: 04 is already running on $OUT"; exit 3; }
rm -f "$OUT/.done.04"     # written again only when every library is segmented and plotted
DOCKER="docker run --platform linux/amd64 --rm $DOCKER_MOUNTS $GATK_IMAGE"
GATK="$DOCKER gatk --java-options -Xmx16g"
for f in "$REF" "$DICT" "$BINS" "$UMAP" "$SEGDUP" "$GNOMAD" "$SNPF"; do [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }; done
for s in "${ALL[@]}"; do [[ -s "$OUT/bams/$s.taps.bam" ]] || { echo "ERROR: $OUT/bams/$s.taps.bam missing (run 02)"; exit 1; }; done
if python3 -c "import pysam" 2>/dev/null; then PY="python3"; else PY="$DOCKER python"; fi
log() { echo "$(date '+%F %T') $*"; }
fresh() {  # fresh <output> <input...>: the output exists and is newer than every input
  local o="$1" i; shift; [[ -s "$o" ]] || return 1
  for i in "$@"; do [[ "$o" -nt "$i" ]] || return 1; done
}
run() {    # run <label> <output> "<inputs>" <command...>: skipped while the output is fresh
  local label="$1" out="$2" deps="$3"; shift 3
  fresh "$out" $deps && return 0
  if "$@" > "$C/logs/$label.log" 2>&1; then log "OK   $label"; else log "FAIL $label (see $C/logs/$label.log)"; rm -f "$out"; return 1; fi
}
need() {   # need <what> <file...>: stop when a parallel step left an output missing
  local what="$1" f; shift
  for f in "$@"; do [[ -s "$f" ]] || { log "FAIL: $what missing ($f); re-run to retry"; exit 1; }; done
}
pool() { while (( $(jobs -rp | wc -l) >= ${1:-$JOBS} )); do sleep 5; done; }

# 0. bins; a test uses only the test contig (the copy-number tools accept only the union of -L intervals)
if [[ -n "$TEST" ]]; then
  sub="$C/bins.$TEST.interval_list"
  [[ -s "$sub" ]] || awk -v c="$TEST" '/^@/ || $1 == c' "$BINS" > "$sub"
  BINS="$sub"
fi

# 1. bins annotated with GC, mappability, segmental duplications; only well-behaved bins are kept
#    (FilterIntervals: mappability >= 0.9, segmental-duplication content <= 0.5, GC 0.1-0.9), which
#    removes the unmappable bins (acrocentric arms, pericentromeres, false duplications such as
#    U2AF1 / U2AF1L5) that otherwise give zero-count bins and tiny spurious segments
ANN_ALL="$C/bins.all.annotated.tsv"
run annotate_all "$ANN_ALL" "$BINS" $GATK AnnotateIntervals -R "$REF" -L "$BINS" --interval-merging-rule OVERLAPPING_ONLY \
    --mappability-track "$UMAP" --segmental-duplication-track "$SEGDUP" -O "$ANN_ALL" || exit 1
KEEP="$C/bins.kept.interval_list"
run filter_bins "$KEEP" "$ANN_ALL" $GATK FilterIntervals -L "$BINS" --annotated-intervals "$ANN_ALL" \
    --minimum-mappability 0.9 --maximum-segmental-duplication-content 0.5 --minimum-gc-content 0.1 --maximum-gc-content 0.9 \
    --interval-merging-rule OVERLAPPING_ONLY -O "$KEEP" || exit 1
log "bins kept: $(grep -vc '^@' "$KEEP") of $(grep -vc '^@' "$BINS")"
ANN="$C/bins.annotated.tsv"
run annotate "$ANN" "$KEEP" $GATK AnnotateIntervals -R "$REF" -L "$KEEP" --interval-merging-rule OVERLAPPING_ONLY \
    --mappability-track "$UMAP" --segmental-duplication-track "$SEGDUP" -O "$ANN" || exit 1

# 2. read counts in the kept bins
for s in "${ALL[@]}"; do pool
  run "counts.$s" "$C/counts/$s.counts.hdf5" "$OUT/bams/$s.taps.bam $KEEP" $GATK CollectReadCounts -I "$OUT/bams/$s.taps.bam" -L "$KEEP" \
      --interval-merging-rule OVERLAPPING_ONLY --minimum-mapping-quality "$MMQ" -O "$C/counts/$s.counts.hdf5" &
done; wait
need "read counts" $(for s in "${ALL[@]}"; do echo "$C/counts/$s.counts.hdf5"; done)

# 3. panels of normals: every library of the panel, and leave-one-out for each of them
pon_args() { local skip="$1" a=(); for n in "${NORMALS[@]}"; do [[ "$n" == "$skip" ]] || a+=(-I "$C/counts/$n.counts.hdf5"); done; echo "${a[@]}"; }
pon_deps() { local skip="$1" a=("$ANN"); for n in "${NORMALS[@]}"; do [[ "$n" == "$skip" ]] || a+=("$C/counts/$n.counts.hdf5"); done; echo "${a[@]}"; }
run pon.all "$C/pon/pon.all.hdf5" "$(pon_deps none)" $GATK CreateReadCountPanelOfNormals $(pon_args none) \
    --annotated-intervals "$ANN" -O "$C/pon/pon.all.hdf5" || exit 1
for n in "${NORMALS[@]}"; do
  run "pon.without_$n" "$C/pon/pon.without_$n.hdf5" "$(pon_deps "$n")" $GATK CreateReadCountPanelOfNormals $(pon_args "$n") \
      --annotated-intervals "$ANN" -O "$C/pon/pon.without_$n.hdf5" || exit 1
done

# 4. denoised copy ratios (GC correction comes from the panel, which holds the annotations)
for s in "${ALL[@]}"; do pool
  pon="$C/pon/pon.all.hdf5"; [[ " ${NORMALS[*]} " == *" $s "* ]] && pon="$C/pon/pon.without_$s.hdf5"
  run "denoise.$s" "$C/denoised/$s.denoisedCR.tsv" "$C/counts/$s.counts.hdf5 $pon" $GATK DenoiseReadCounts \
      -I "$C/counts/$s.counts.hdf5" --count-panel-of-normals "$pon" \
      --standardized-copy-ratios "$C/denoised/$s.standardizedCR.tsv" --denoised-copy-ratios "$C/denoised/$s.denoisedCR.tsv" &
done; wait
need "denoised copy ratios" $(for s in "${ALL[@]}"; do echo "$C/denoised/$s.denoisedCR.tsv"; done)

# 5. common SNP sites: gnomAD AF 0.1-0.9, biallelic SNVs, chr1-22 and X in dictionary order,
#    without SNPs whose alternative allele creates a CpG
SITES="$C/common_snps.vcf.gz"
if ! fresh "$SITES" "$GNOMAD"; then
  regions=$(awk '/^@SQ/{for(i=1;i<=NF;i++) if($i~/^SN:/){n=substr($i,4); if (n ~ /^chr([0-9]+|X)$/) printf "%s%s", (c++?",":""), n}}' "$DICT")
  [[ -n "$TEST" ]] && regions="$TEST"
  rm -f "$SITES" "$SITES.tbi"
  $DOCKER bash -c "bcftools view -r $regions -v snps -m2 -M2 -i 'INFO/AF>=0.1 && INFO/AF<=0.9' -Oz -o $SITES.all.vcf.gz $GNOMAD" \
    && $PY "$SNPF" "$SITES.all.vcf.gz" "$REF" "$SITES.tmp.vcf.gz" > "$C/logs/sites_filter.log" 2>&1 \
    && mv "$SITES.tmp.vcf.gz" "$SITES" && rm -f "$SITES.all.vcf.gz" \
    && $GATK IndexFeatureFile -I "$SITES" > "$C/logs/sites_index.log" 2>&1 \
    && log "OK   common SNP sites: $(tail -1 "$C/logs/sites_filter.log")" \
    || { log "FAIL common SNP sites (see $C/logs/sites_filter.log)"; rm -f "$SITES"*; exit 1; }
fi

# 6. allelic counts (base quality >= 20, so the masked methylation bases are ignored)
for s in "${ALL[@]}"; do pool
  run "allelic.$s" "$C/allelic/$s.allelicCounts.tsv" "$OUT/bams/$s.taps.bam $SITES" $GATK CollectAllelicCounts \
      -I "$OUT/bams/$s.taps.bam" -R "$REF" -L "$SITES" --minimum-base-quality 20 --minimum-mapping-quality "$MMQ" \
      -O "$C/allelic/$s.allelicCounts.tsv" &
done; wait
need "allelic counts" $(for s in "${ALL[@]}"; do echo "$C/allelic/$s.allelicCounts.tsv"; done)

# 7-9. segments, calls and plots, SEG_JOBS libraries at a time; tumours with their matched normal's
#      heterozygous sites, normals alone
segment() {   # segment <sample>: ModelSegments -> CallCopyRatioSegments -> PlotModeledSegments
  local s="$1" nac seg="$C/segments/$1.modelFinal.seg" cr="$C/denoised/$1.denoisedCR.tsv" ac="$C/allelic/$1.allelicCounts.tsv"
  if [[ -n "${MATCHED[$s]:-}" ]]; then
    nac="$C/allelic/${MATCHED[$s]}.allelicCounts.tsv"
    run "segments.$s" "$seg" "$cr $ac $nac" $GATK ModelSegments --denoised-copy-ratios "$cr" --allelic-counts "$ac" \
        --normal-allelic-counts "$nac" --minimum-total-allele-count-normal 10 --output "$C/segments" --output-prefix "$s" || return 1
  else
    run "segments.$s" "$seg" "$cr $ac" $GATK ModelSegments --denoised-copy-ratios "$cr" --allelic-counts "$ac" \
        --minimum-total-allele-count-case 10 --output "$C/segments" --output-prefix "$s" || return 1
  fi
  run "called.$s" "$C/segments/$s.called.seg" "$seg" $GATK CallCopyRatioSegments -I "$C/segments/$s.cr.seg" \
      -O "$C/segments/$s.called.seg" || return 1
  run "plot.$s" "$C/plots/$s.modeled.png" "$seg" $GATK PlotModeledSegments --denoised-copy-ratios "$cr" \
      --allelic-counts "$C/segments/$s.hets.tsv" --segments "$seg" --sequence-dictionary "$DICT" \
      --minimum-contig-length 46709983 --output "$C/plots" --output-prefix "$s"
}
for s in "${ALL[@]}"; do pool "$SEG_JOBS"; segment "$s" & done; wait
FAILED=()
for s in "${ALL[@]}"; do   # a library is finished when its calls and plot are newer than its segments
  seg="$C/segments/$s.modelFinal.seg"
  fresh "$C/segments/$s.called.seg" "$seg" && fresh "$C/plots/$s.modeled.png" "$seg" || FAILED+=("$s")
done
(( ${#FAILED[@]} )) && { log "FAILED: ${FAILED[*]} (re-run to retry; finished work is kept)"; exit 1; }
date '+%F %T' > "$OUT/.done.04"
log "DONE CNV/LOH: $C/segments, $C/plots"
