################################################################################
# 04_enrichment.R  --  sections 8 and 9 of the standard workflow
#
#   Rscript 03_rnaseq/04_enrichment.R genes            # gene level (default)
#   Rscript 03_rnaseq/04_enrichment.R 5utr             # 5'UTR level
#   Rscript 03_rnaseq/04_enrichment.R genes HALLMARK   # one category only
#
# Section 8  ORA: GO-BP/CC/MF + KEGG + Hallmark, up and down separately
#            -> enrichment/<contrast>/<direction>/
#               enrich_<CAT>.csv, dotplot_<CAT>.pdf/png, cnetplot_<CAT>.pdf/png
#
# Section 9  GSEA: GO-BP/CC/MF + KEGG + Hallmark, on the whole ranked list
#            -> enrichment_gsea/<contrast>/
#               gsea_full_<CAT>.csv   every gene set tested, raw and adjusted p
#               gsea_<CAT>.csv        sets with raw p <= 0.05 (the figures use these)
#               dotplot_<CAT>, ridgeplot_<CAT>, curves/<CAT>/NN_<term>
#
# Then run 04b_enrichment_ma_plots.R: it draws a matched MA plot for every
# enriched term and the
# GO-BP figure across the two axes, from these CSVs, without recomputing.
#
# TWO IMPLEMENTATION DETAILS
#   1. clusterProfiler's pvalueCutoff drops every set above it from the output.
#      The reference-contrast filter in 05_pathway_ruler.R tests at raw
#      p < 0.10, so full tables (every set tested) are written as well.
#   2. GSEA samples at random, and seed = TRUE does not make repeated calls
#      identical. GSEA therefore runs serially with set.seed() before every
#      call (gsea_one() in dge_functions.R), which makes it identical run to run.
#
# ORA vs GSEA. ORA needs a hard DE cutoff, and with one library per timepoint
# those lists are short and uncertain. GSEA uses the whole ranked list and never
# has to decide whether a gene is "in". Where the two disagree, believe GSEA.
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(ggplot2); library(stringr); library(ggnewscale)
  library(DOSE); library(msigdbr); library(BiocParallel)
})

SCRIPTS <- "03_rnaseq"
source(file.path(SCRIPTS, "00_config.R"))
source(file.path(SCRIPTS, "dge_functions.R"))
register(SerialParam(), default = TRUE)   # reproducible GSEA -- see note 2 above

args  <- commandArgs(trailingOnly = TRUE)
LEVEL <- if (length(args)) args[1] else "genes"
stopifnot(LEVEL %in% c("genes", "5utr"))

# Optional second argument restricts which categories are written, so a single
# category can be redone without the whole (slow) run:
#   Rscript 03_rnaseq/04_enrichment.R genes HALLMARK
CATS <- if (length(args) > 1) strsplit(args[2], ",")[[1]] else c("BP", "CC", "MF", "KEGG", "HALLMARK")
message("categories: ", paste(CATS, collapse = ", "))

# ---- Hallmark ------------------------------------------------------------
# The standard workflow runs GO and KEGG only. Hallmark is added because GO splits
# coherent programmes (interferon response, p53 signalling) across many small,
# overlapping terms, none of which survives multiple testing on its own, while
# Hallmark carries each as a single set. A complement to GO, not a replacement.
.msig_H <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                    error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
entrez_col <- intersect(c("ncbi_gene", "entrez_gene"), names(.msig_H))[1]
HALLMARK_T2G <- data.frame(gs_name = .msig_H$gs_name,
                           entrez_gene = .msig_H[[entrez_col]],
                           stringsAsFactors = FALSE)
message("Hallmark: ", length(unique(HALLMARK_T2G$gs_name)), " gene sets")

OUT   <- if (LEVEL == "genes") PATHS$out_gene else PATHS$out_5utr
RDATA <- if (LEVEL == "genes") "dge_gene_level.RData" else "dge_5utr_level.RData"
if (!file.exists(file.path(OUT, RDATA)))
  stop("run the ", LEVEL, " DGE script first -- ", RDATA, " not found in ", OUT)

e <- new.env(); load(file.path(OUT, RDATA), envir = e)
de <- as.data.table(e$de)
message("enrichment for level '", LEVEL, "': ", uniqueN(de$gene), " features, ",
        uniqueN(de$contrast), " contrasts")

## ---- parameters (defaults as in the standard workflow) ---------------------
ora_fdr         <- THRESH$fdr    # DE selection: adjusted p
ora_lfc         <- THRESH$lfc    # DE selection: |log2FC|
term_p          <- 0.05          # ORA term significance
top_terms       <- 20            # terms per dotplot
split_dir       <- TRUE          # enrich up and down separately
gsea_term_p     <- 0.05          # raw-p cutoff for gsea_<CAT>.csv and its figures
gsea_top_dot    <- 15
gsea_top_curves <- 10
GSEA_SEED       <- 1234

prep <- function(cid) prep_contrast(de[contrast == cid][order(P.Value)])

################################################################################
# 8. Over-representation analysis
################################################################################
message("\n############ SECTION 8: ORA ############")
ora_summary <- list()

for (cid in CONTRASTS$id) {
  message("\n== Enrichment: ", cid, " (", CONTRASTS[cid, label], ") ==")
  tt       <- prep(cid)
  # as.character() strips the attributes na.omit() adds; with them attached,
  # enrichKEGG can silently ignore a custom universe
  universe <- as.character(unique(na.omit(tt$entrez)))
  # on/off calls rest on a single library's detection; they stay in the ranked
  # list for GSEA but are not fed to ORA
  sig      <- subset(tt, adj.P.Val < ora_fdr & abs(lfc_shrunk) >= ora_lfc &
                       !on_off & !is.na(entrez))

  gene_sets <- if (split_dir)
    list(up   = unique(sig$entrez[sig$lfc_shrunk > 0]),
         down = unique(sig$entrez[sig$lfc_shrunk < 0]))
  else list(all = unique(sig$entrez))

  fc <- tt$lfc_shrunk; names(fc) <- tt$symbol
  fc <- sort(fc, decreasing = TRUE)

  for (d in names(gene_sets)) {
    message("  direction: ", d)
    ora <- enrich_one(gene_sets[[d]], universe, term_p)
    if (is.null(ora)) { message("  skipped (too few genes)"); next }

    base <- file.path(OUT, "enrichment", cid, d)
    dir.create(base, showWarnings = FALSE, recursive = TRUE)
    dtitle <- paste0(CONTRASTS[cid, label], " [", d, "]")

    for (cat in CATS) {
      er <- ora[[cat]]
      if (is.null(er) || !nrow(as.data.frame(er))) next
      df <- as.data.frame(er)
      write.csv(df, file.path(base, paste0("enrich_", cat, ".csv")), row.names = FALSE)
      ora_summary[[length(ora_summary) + 1]] <-
        data.table(contrast = cid, direction = d, category = cat, n_terms = nrow(df))
      message("    ", cat, ": ", nrow(df), " terms")
      save_dotplot(er, paste0(dtitle, " - ", cat),
                   file.path(base, paste0("dotplot_", cat)), top_terms)
      save_cnetplot(er, fc, paste0(dtitle, " - ", cat),
                    file.path(base, paste0("cnetplot_", cat)))
    }
  }
}

if (length(ora_summary)) {
  ora_summary <- rbindlist(ora_summary)
  fwrite(ora_summary, file.path(OUT, "47_ora_summary.tsv"), sep = "\t")
  message("\n  ORA summary:"); print(dcast(ora_summary,
        contrast + direction ~ category, value.var = "n_terms", fill = 0))
}

################################################################################
# 9. GSEA (threshold-free)
################################################################################
message("\n############ SECTION 9: GSEA ############")
gsea_summary <- list()

for (cid in CONTRASTS$id) {
  message("\n== GSEA: ", cid, " (", CONTRASTS[cid, label], ") ==")
  tt    <- prep(cid)
  ranks <- rank_genes(tt)
  message("  ranked list: ", length(ranks), " genes")
  gsea  <- gsea_one(ranks, gsea_term_p = 1, seed = GSEA_SEED)   # every set: note 1

  base <- file.path(OUT, "enrichment_gsea", cid)
  dir.create(base, showWarnings = FALSE, recursive = TRUE)
  gtitle <- CONTRASTS[cid, label]

  for (cat in CATS) {
    gs <- gsea[[cat]]
    if (is.null(gs)) next
    n_tested <- nrow(gs@result)
    write.csv(as.data.frame(gs), file.path(base, paste0("gsea_full_", cat, ".csv")),
              row.names = FALSE)
    gs@result <- gs@result[!is.na(gs@result$pvalue) & gs@result$pvalue <= gsea_term_p, ,
                           drop = FALSE]
    df <- as.data.frame(gs)
    write.csv(df, file.path(base, paste0("gsea_", cat, ".csv")), row.names = FALSE)
    gsea_summary[[length(gsea_summary) + 1]] <-
      data.table(contrast = cid, category = cat, tested = n_tested,
                 raw_p05 = nrow(df), fdr05 = sum(df$p.adjust < THRESH$gsea_fdr),
                 up_fdr05 = sum(df$NES > 0 & df$p.adjust < THRESH$gsea_fdr),
                 down_fdr05 = sum(df$NES < 0 & df$p.adjust < THRESH$gsea_fdr))
    message("    ", cat, ": ", n_tested, " tested, ", nrow(df), " at raw p <= 0.05, ",
            sum(df$p.adjust < THRESH$gsea_fdr), " at FDR < 0.05")
    if (!nrow(df)) next

    save_gsea_dotplot(gs, paste0(gtitle, " - GSEA ", cat),
                      file.path(base, paste0("dotplot_", cat)), gsea_top_dot)
    save_ridgeplot(gs, paste0(gtitle, " - GSEA ", cat),
                   file.path(base, paste0("ridgeplot_", cat)), gsea_top_dot)
    save_gsea_curves(gs, cat, file.path(base, "curves", cat), gsea_top_curves)
  }
}

if (length(gsea_summary)) {
  gsea_summary <- rbindlist(gsea_summary)
  fwrite(gsea_summary, file.path(OUT, "48_gsea_summary.tsv"), sep = "\t")
  message("\n  GSEA summary:"); print(gsea_summary)
}

# all sets at raw p <= 0.05, every contrast, one table
gsea_files <- list.files(file.path(OUT, "enrichment_gsea"), pattern = "^gsea_[A-Z]+\\.csv$",
                         recursive = TRUE, full.names = TRUE)
if (length(gsea_files)) {
  allg <- rbindlist(lapply(gsea_files, function(f) {
    d <- as.data.table(read.csv(f))
    if (!nrow(d)) return(NULL)
    d[, `:=`(contrast = basename(dirname(f)),
             category = sub("^gsea_(.*)\\.csv$", "\\1", basename(f)))]
    d[, .(contrast, category, ID, Description, NES, pvalue, p.adjust, setSize)]
  }), fill = TRUE)
  fwrite(allg, file.path(OUT, "49_gsea_all_contrasts.tsv"), sep = "\t")
}

writeLines(c(run_stamp(),
             paste0("GSEA: serial (", class(bpparam())[1], "), set.seed(", GSEA_SEED,
                    ") before every call, eps = 0, BH"),
             paste0("KEGG gene sets downloaded by clusterProfiler at run time: ",
                    format(Sys.time(), "%Y-%m-%d %H:%M"))),
           file.path(OUT, "48b_enrichment_run_info.txt"))

message("\n== done ==  ", run_stamp())
message("ORA  -> ", file.path(OUT, "enrichment"))
message("GSEA -> ", file.path(OUT, "enrichment_gsea"))
message("next: Rscript 03_rnaseq/04b_enrichment_ma_plots.R ", LEVEL)
