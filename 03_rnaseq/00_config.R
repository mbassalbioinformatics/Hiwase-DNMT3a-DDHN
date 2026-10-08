################################################################################
# 00_config.R  --  shared configuration for the RNA-seq analysis
#
# Sourced by every RNA-seq script (run from the repository root). Nothing here
# runs an analysis: it defines the samples, the contrasts, thresholds and the
# gene sets used to check for confounding, so that each exists in one place.
# Samples and contrasts are read from config/ (samples.tsv, rnaseq_contrasts.tsv),
# which hold the sample details and are not committed.
#
# DESIGN
#   One library per timepoint from the recipient (day 30, 3.5 years, DDHN
#   diagnosis) and one from the donor. The timepoints give two separate axes:
#     clonal establishment   day 30 -> 3.5 yr   the clone expands, no disease
#     transformation         3.5 yr -> DDHN     clone size flat, disease appears
#   Anything moving on the second axis is not a function of clone size.
#
# n = 1 PER TIMEPOINT
#   There are no replicates. Consequently, throughout:
#     * edgeR runs with a FIXED dispersion (BCV grid), never an estimated one;
#     * p and FDR rank genes; they are not calibrated error rates;
#     * confidence comes from agreement between independent contrasts and from
#       threshold-free set-level methods (GSEA, ssGSEA).
#
# CONFOUNDING
#   Globin and mitochondrial fractions differ widely between libraries; both are
#   removed before normalisation and reported. The donor library is peripheral
#   blood and the recipient libraries bone marrow, so every donor contrast is
#   compartment-confounded and is read against the reference contrast (day 30 vs
#   donor: no disease, a small clone at both ends).
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(RColorBrewer)
})
source("config/paths.R")

# ---- paths -----------------------------------------------------------------
PROJ <- RNA
PATHS <- list(
  proj        = PROJ,
  counts_gene = COUNTS_GENE,
  counts_5utr = COUNTS_5UTR,
  out_gene    = file.path(PROJ, "05_processed_counts_genes"),
  out_5utr    = file.path(PROJ, "05_processed_counts_5utr"),
  out_rmats   = file.path(PROJ, "06_rmats"),
  out_shiba   = file.path(PROJ, "07_shiba"),
  out_reports = file.path(PROJ, "08_reports"),
  scripts     = "03_rnaseq"
)
for (p in PATHS[grep("^out_", names(PATHS))])
  dir.create(p, showWarnings = FALSE, recursive = TRUE)

# ---- library preparation ---------------------------------------------------
# Stranded poly(A) libraries, read 1 in the sense orientation:
#   featureCounts -s 1 | rMATS fr-secondstrand | RSeQC "1++,1--,2+-,2-+"
SEQ <- list(
  kit           = "Tecan/NuGEN Universal Plus mRNA-Seq with NuQuant (part 0508), protocol M01442 v2",
  platform      = "MGI DNBSEQ-G400, paired-end 100 bp",
  strandedness  = "stranded, forward (read 1 sense)",
  fc_strand     = 1,
  rmats_libtype = "fr-secondstrand"
)

# ---- samples ---------------------------------------------------------------
# From config/samples.tsv, rows with rna_seq == "yes", in the order of the BAM
# columns of the count tables. Downstream code checks the order by name.
.sheet <- fread(SAMPLE_SHEET, colClasses = "character")[rna_seq == "yes"]
SAMPLES <- .sheet[, .(sample, subject, tp_label, short,
                      sample_date = as.Date(sample_date), material,
                      vaf = as.numeric(fifelse(nzchar(vaf_pct), vaf_pct, NA_character_)),
                      vaf_source, disease, clone, clone_ann)]
SAMPLE_ORDER <- SAMPLES$sample
PATIENT_SAMPLES <- SAMPLES[subject == "patient", sample]
setkey(SAMPLES, sample)

check_sample_order <- function(mat) {
  # Never trust column position. Fails loudly rather than mislabel silently.
  if (!identical(colnames(mat), SAMPLE_ORDER))
    stop("Count matrix columns are ", paste(colnames(mat), collapse = ", "),
         " but SAMPLE_ORDER is ", paste(SAMPLE_ORDER, collapse = ", "),
         " -- refusing to continue with ambiguous sample identity.")
  invisible(TRUE)
}

# ---- contrasts -------------------------------------------------------------
# From config/rnaseq_contrasts.tsv: id, num, den, label, class, meaning.
# Direction is ALWAYS numerator - denominator: positive logFC / dPSI = higher
# in `num`, in the gene, 5'UTR, rMATS and Shiba arms alike.
CONTRASTS <- fread(RNA_CONTRASTS, colClasses = "character")
setkey(CONTRASTS, id)

PRIMARY_CONTRASTS <- CONTRASTS[class == "primary", id]
RULER_CONTRAST    <- CONTRASTS[class == "ruler",   id]
# Both libraries steady-state marrow from the same patient, clone size flat:
# the best-controlled contrast for regeneration and tissue (not for cell mix).
CLEANEST_CONTRAST <- "C2_transform"

# ---- confound signatures ---------------------------------------------------
# Scored per sample and reported next to every result, because each is a live
# alternative explanation in this dataset:
#   Proliferation   day-30 marrow is regenerating after transplant
#   OxPhos_nuclear  nuclear-encoded mitochondrial genes (MT- genes are removed)
#   Erythroid_prog  erythroid precursor content
#   Interferon      interferon response
CONFOUND_SIGNATURES <- list(
  Proliferation  = c("MKI67", "TOP2A", "PCNA", "CCNB1", "CDK1", "TYMS", "RRM2",
                     "AURKB", "BIRC5", "MCM2", "TK1", "CENPF", "UBE2C"),
  OxPhos_nuclear = c("MRPL12", "MRPS12", "NDUFA4", "NDUFB1", "COX7B", "COX5A",
                     "ATP5F1E", "ATP5MC2", "UQCRB", "SDHB", "COX6C", "NDUFS6"),
  Erythroid_prog = c("AHSP", "ALAS2", "SLC4A1", "GYPA", "KLF1", "EPB42", "CA1",
                     "GATA1", "RHAG", "SPTA1"),
  Interferon     = c("IFI27", "OASL", "IFIT3", "OAS1", "RSAD2", "USP18", "MX1",
                     "ISG20", "IFITM1", "GBP2", "GBP4", "DDX60", "EPSTI1",
                     "IFIT1", "ISG15", "STAT1")
)

# ---- thresholds ------------------------------------------------------------
THRESH <- list(
  # expression filter: a feature must clear this CPM in at least this many
  # libraries to be tested
  cpm            = 1.0,
  min_samples    = 2,
  # edgeR fixed-dispersion grid. 0.4 is the edgeR User's Guide value for human
  # data and is the primary setting; 0.2 and 0.6 are a sensitivity band,
  # reported next to every headline count. Checked against the data in
  # 02_dge_gene_level.R (pooled-library estimates and the spread of a typical
  # expressed gene).
  bcv_primary    = 0.4,
  bcv_grid       = c(0.2, 0.4, 0.6),
  # NOISeq probability cutoff, chosen from the observed distribution
  # (11b_noiseq_threshold_curve.tsv): at this depth no-replicate NOISeq models
  # little technical noise and the conventional 0.9 keeps most genes.
  noiseq_prob    = 0.999,
  # A gene is FLAT on an axis only if its fold-change is small, not merely
  # non-significant.
  flat_lfc       = 0.5,
  fdr            = 0.05,
  lfc            = 1.0,      # |log2FC| >= 1: with n = 1, smaller changes are not defensible
  # 5'UTR / promoter usage
  utr_min_count  = 20,       # per-feature minimum, summed across samples
  utr_dprop      = 0.15,     # |change in relative usage| for a promoter shift
  utr_min_width  = 25,       # shorter 5'UTRs collect reads from the next exon under -O
  # splicing (rMATS and Shiba alike)
  psi_dpsi       = 0.10,
  psi_minreads   = 20,
  # enrichment
  gsea_fdr       = 0.05,
  ora_fdr        = 0.05,
  # The reference contrast is read at RAW p, deliberately lenient: its job is
  # to rule a compartment explanation out. It needs the complete GSEA tables
  # (gsea_full_<CAT>.csv); gsea_<CAT>.csv holds only sets with p <= 0.05.
  ruler_p        = 0.10
)

# ---- confound-control gene sets --------------------------------------------
# Removed from the count matrix before normalisation, and reported.
GLOBIN_GENES <- c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2",
                  "HBM", "HBQ1", "HBZ", "HBE1")
is_mt_gene <- function(gene_name) grepl("^MT-", gene_name)

# ---- lineage marker sets for composition scoring ---------------------------
# Marker-set scoring (GSVA/ssGSEA) rather than reference-based deconvolution:
# LM22 / quanTIseq / EPIC signatures come from PBMC and solid tumours and are
# out of domain for bone marrow with CD34+ progenitors and erythroid
# precursors. Every gene is listed here.
LINEAGE_MARKERS <- list(
  HSC_MPP        = c("CD34", "CRHBP", "AVP", "HLF", "MECOM", "PROM1", "MLLT3",
                     "HOPX", "NPR3"),
  Erythroid      = c("AHSP", "ALAS2", "SLC4A1", "GYPA", "KLF1", "EPB42", "CA1",
                     "GATA1", "RHAG", "SPTA1", "TFR2"),
  Granulo_precur = c("ELANE", "MPO", "PRTN3", "AZU1", "CTSG", "DEFA4", "RNASE2",
                     "MS4A3", "SERPINB10"),
  Neutrophil     = c("FCGR3B", "CSF3R", "CXCR2", "S100A8", "S100A9", "LCN2",
                     "CAMP", "LTF", "CEACAM8"),
  Monocyte       = c("LYZ", "CD14", "VCAN", "S100A12", "FCN1", "CSF1R", "CD68",
                     "MNDA"),
  T_cell         = c("CD3D", "CD3E", "CD3G", "IL7R", "TRAC", "CD2", "LCK",
                     "CD28", "THEMIS"),
  NK_cell        = c("GNLY", "NKG7", "KLRD1", "PRF1", "GZMB", "NCAM1", "KLRF1"),
  B_cell         = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "TCL1A",
                     "VPREB1", "DNTT"),
  Plasma_cell    = c("MZB1", "JCHAIN", "DERL3", "TNFRSF17", "PRDM1", "XBP1"),
  Megakaryocyte  = c("PF4", "PPBP", "ITGA2B", "GP9", "VWF", "TUBB1", "PLEK"),
  Dendritic      = c("FCER1A", "CD1C", "CLEC9A", "IRF8", "LILRA4", "CLEC4C")
)

# ---- genes of direct interest ----------------------------------------------
# Curated panels, not results: shown so that genes the project is about can be
# inspected whether or not they reach a threshold. Provenance for each below.
FOCUS_GENES <- list(
  # DNA methylation writers, erasers and readers (DNMT/TET/IDH families plus
  # canonical methyl-CpG readers)
  DNA_methylation = c("DNMT3A", "DNMT3B", "DNMT1", "TET1", "TET2", "TET3",
                      "IDH1", "IDH2", "UHRF1", "MBD1", "MBD2", "MECP2"),
  # Myeloid drivers on standard diagnostic myeloid NGS panels (ELN 2022 AML,
  # WHO/ICC myeloid classifications)
  Myeloid_drivers = c("ASXL1", "RUNX1", "SF3B1", "SRSF2", "U2AF1", "ZRSR2",
                      "EZH2", "STAG2", "BCOR", "TP53", "NPM1", "FLT3",
                      "CEBPA", "GATA2", "NRAS", "KRAS", "JAK2", "CBL",
                      "DDX41", "PPM1D", "CHEK2", "SETBP1", "ETV6", "PHF6"),
  # Core spliceosome components and splicing regulators
  Splicing_factor = c("SF3B1", "SRSF2", "U2AF1", "U2AF2", "ZRSR2", "SRRM1",
                      "SRRM2", "HNRNPU", "HNRNPA1", "RBFOX2", "PRPF8", "SNRNP200",
                      "LUC7L2", "DDX41"),
  # Canonical direct p53 transcriptional targets, intersected with the Hallmark
  # p53 pathway leading edge observed here
  p53_targets     = c("CDKN1A", "MDM2", "BAX", "BAK1", "ZMAT3", "AEN", "FDXR",
                      "RCHY1", "DDIT4", "TP53I3", "SESN1", "SESN2", "TNFRSF10B",
                      "TIGAR", "GADD45A", "BTG2", "TXNIP", "BTG1", "PPP1R15A")
)

FOCUS_PROVENANCE <- c(
  DNA_methylation = "DNMT/TET/IDH enzyme families plus canonical methyl-CpG readers.",
  Myeloid_drivers = "Genes on standard diagnostic myeloid NGS panels, per the 2022 European LeukemiaNet AML recommendations and current WHO/ICC myeloid classifications.",
  Splicing_factor = "Core spliceosome components and established splicing regulators.",
  p53_targets     = "Canonical direct p53 transcriptional targets, intersected with the MSigDB Hallmark p53 pathway leading edge observed in this dataset."
)

# ---- palettes --------------------------------------------------------------
PAL <- list(
  sample = setNames(.sheet$rna_colour, .sheet$sample),
  class  = c(primary = "#C44E52", context = "#8C8C8C", ruler = "#4C72B0"),
  up     = "#C44E52", down = "#4C72B0", ns = "#BDBDBD"
)
# clone legend colours, keyed by the short clone label of the sample sheet, untested first, then by VAF
CLONE_COLOURS <- unique(.sheet[nzchar(clone_ann), .(clone_ann, clone_colour, v = suppressWarnings(as.numeric(vaf_pct)))])[order(!is.na(v), v)]
CLONE_COLOURS <- setNames(CLONE_COLOURS$clone_colour, CLONE_COLOURS$clone_ann)

# ---- provenance stamp ------------------------------------------------------
run_stamp <- function() {
  paste0("generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), " | R ",
         getRversion(), " | edgeR ", as.character(packageVersion("edgeR")))
}

message("00_config.R loaded: ", nrow(SAMPLES), " samples, ",
        nrow(CONTRASTS), " contrasts (", length(PRIMARY_CONTRASTS), " primary).")
