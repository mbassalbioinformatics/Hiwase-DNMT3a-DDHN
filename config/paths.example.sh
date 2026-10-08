###############################################################################
# paths.example.sh
#
# Copy to config/paths.sh and edit. Sourced by the server scripts in
# 01_methylation_server/ and 03_rnaseq/server/. config/paths.sh is not committed.
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
