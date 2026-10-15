###############################################################################
# paths.example.sh
#
# Copy to config/paths.sh and edit. Sourced by the server scripts in
# 01_methylation_server/, 03_rnaseq/server/ and 05_variants_server/.
# config/paths.sh is not committed.
###############################################################################
THREADS=31

# ---- methylation (TAPS) ------------------------------------------------------
FASTQ_DIRS=(/path/to/emseq/01_raw_fq_1 /path/to/emseq/01_raw_fq_2)   # every sequencing run
MSUITE2=/path/to/Msuite2/msuite2
MS_OUT=/path/to/emseq/02_msuite2_analysis       # one <sample>_msuite2/ folder per sample
PUC19_DIR=$MS_OUT                               # Msuite2.pUC19.size, Msuite2.c.pUC19.size, Msuite2.w.pUC19.size
DMR_OUT=/path/to/emseq/03_dmr_analysis
HOMER_GENOME=hg38                               # a genome installed in HOMER, or an hg38 FASTA

# ---- RNA-seq -------------------------------------------------------------------
BAM_DIR=/path/to/rnaseq/03_aligned_bams
BAM_SUFFIX=_qual_trimmed.sorted.bam
GTF=/path/to/Homo_sapiens.GRCh38.107.chr.gtf
RNA_OUT=/path/to/rnaseq                         # counts, rMATS and Shiba output

# ---- genome-wide variants (05_variants_server/) ---------------------------------
VAR_OUT=/path/to/emseq/04_variants              # masked BAMs, Mutect2 and copy-number runs; transfer/ goes to the analysis computer
REF=/path/to/hg38.fa                            # the FASTA Msuite2 aligned to, with .fai and a GATK .dict beside it
GATK_IMAGE=broadinstitute/gatk:4.6.2.0          # Docker image; its python (pysam) is used when the host has none
DOCKER_MOUNTS="-v /path:/path"                  # every folder in this file, mounted at the same path in the container
# GATK somatic best-practice resources (hg38), contigs in the order of the reference dictionary
MUTECT2_PON=/path/to/1000g_pon.hg38.vcf.gz
GERMLINE_RESOURCE=/path/to/af-only-gnomad.hg38.vcf.gz
COMMON_SITES=/path/to/small_exac_common_3.hg38.vcf.gz
CNV_BINS=/path/to/hg38_1000bp_bins.interval_list   # 1-kb bins, e.g. GATK PreprocessIntervals --bin-length 1000 --padding 0
UMAP_BED=/path/to/k100.umap.hg38.bed.gz         # Umap k100 single-read mappability (bismap.hoffmanlab.org)
SEGDUP_BED=/path/to/hg38.segdups.bed.gz         # segmental duplications (UCSC genomicSuperDups), overlaps merged
FUNCOTATOR_DIR=/path/to/funcotator_dataSources.v1.8.hg38.20230908s   # GATK Funcotator somatic data sources
