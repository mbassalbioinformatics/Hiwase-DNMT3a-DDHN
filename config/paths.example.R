###############################################################################
# paths.example.R
#
# Copy to config/paths.R and edit. Every R script is run from the repository
# root and reads config/paths.R first. config/paths.R is not committed.
###############################################################################
DATA_ROOT <- "/path/to/project"                  # the project data folder

EMS  <- file.path(DATA_ROOT, "emseq")            # <sample>_msuite2/ per sample, 03_dmr_analysis/ (from the server)
RNA  <- file.path(DATA_ROOT, "rnaseq")           # count tables; RNA-seq results are written here
INT  <- file.path(DATA_ROOT, "integration")      # integration results and reports are written here
RES  <- file.path(INT, "resources")              # annotation downloads (see resources/README.md)

COUNTS_GENE <- file.path(RNA, "GENES.counts_ensembl_FC")   # featureCounts, gene level
COUNTS_5UTR <- file.path(RNA, "5UTR.counts_ensembl_FC")    # featureCounts, 5'UTR level

MANE    <- "/path/to/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz"
HIC_DIR <- "/path/to/hic"                        # *_Pimanda_*.bedpe.xz (HSC, CMP, GMP, MEP; two replicates each)

SAMPLE_SHEET    <- "config/samples.tsv"
COMPARISONS     <- "config/comparisons.tsv"
RNA_CONTRASTS   <- "config/rnaseq_contrasts.tsv"

# Optional: point INT or RES elsewhere without editing this file, e.g. for a test run
INT <- Sys.getenv("DDHN_INT", INT)
RES <- Sys.getenv("DDHN_RES", RES)
