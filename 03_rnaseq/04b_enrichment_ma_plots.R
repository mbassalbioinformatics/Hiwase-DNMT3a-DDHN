################################################################################
# 04b_enrichment_ma_plots.R  --  figures drawn from the saved enrichment CSVs
#
#   Rscript 03_rnaseq/04b_enrichment_ma_plots.R genes
#   Rscript 03_rnaseq/04b_enrichment_ma_plots.R 5utr
#
# Needs 04_enrichment.R to have been run for the level. Recomputes nothing.
#
# 1. A matched MA plot for every enriched term: the contrast's own MA plot with
#    the term's genes (ORA) or leading-edge genes (GSEA) highlighted.
#       enrichment/<contrast>/<dir>/ma_plots/<CAT>/NN_<term>.png/pdf
#       enrichment_gsea/<contrast>/ma_plots/<CAT>/<activated|suppressed>_NN_<term>.png/pdf
#    MA plots rather than volcanoes: with one library per group and a supplied
#    dispersion a volcano's two axes are almost the same quantity (see
#    plot_ma_contrast() in dge_functions.R). Earlier per-term figure folders
#    are cleared before redrawing.
#
# 2. 50_gsea_axes_GOBP.png: GO Biological Process terms across the within-patient
#    comparisons: the top terms by FDR from EACH axis, every tile coloured
#    from the full tables, and the subtitle states how many terms each axis has.
################################################################################

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(stringr) })
SCRIPTS <- "03_rnaseq"
source(file.path(SCRIPTS, "00_config.R"))
source(file.path(SCRIPTS, "dge_functions.R"))

args  <- commandArgs(trailingOnly = TRUE)
LEVEL <- if (length(args)) args[1] else "genes"
stopifnot(LEVEL %in% c("genes", "5utr"))
OUT   <- if (LEVEL == "genes") PATHS$out_gene else PATHS$out_5utr
RDATA <- if (LEVEL == "genes") "dge_gene_level.RData" else "dge_5utr_level.RData"
e <- new.env(); load(file.path(OUT, RDATA), envir = e)
de <- as.data.table(e$de)
CATS <- c("BP", "CC", "MF", "KEGG", "HALLMARK")
set.seed(1234)   # ggrepel label placement

################################################################################
# 1. Matched MA plots
################################################################################
old <- c(list.dirs(file.path(OUT, "enrichment"), recursive = TRUE),
         list.dirs(file.path(OUT, "enrichment_gsea"), recursive = TRUE))
old <- old[basename(old) %in% c("volcanoes", "ma_plots")]
if (length(old)) {
  message("removing ", length(old), " old per-term figure folders (volcanoes / ma_plots)")
  unlink(old, recursive = TRUE)
}

tally <- list()
for (cid in CONTRASTS$id) {
  tt <- de[contrast == cid]
  for (d in c("up", "down")) for (ct in CATS) {
    f <- file.path(OUT, "enrichment", cid, d, paste0("enrich_", ct, ".csv"))
    if (!file.exists(f)) next
    df <- read.csv(f)
    n  <- category_ma_plots(tt, df, ct, file.path(OUT, "enrichment", cid, d, "ma_plots", ct))
    tally[[length(tally) + 1]] <- data.table(contrast = cid, analysis = "ORA", direction = d,
                                             category = ct, terms_in_csv = nrow(df),
                                             plots = n)
  }
  for (ct in CATS) {
    f <- file.path(OUT, "enrichment_gsea", cid, paste0("gsea_", ct, ".csv"))
    if (!file.exists(f)) next
    df <- read.csv(f)
    n  <- gsea_term_ma_plots(tt, df, ct, file.path(OUT, "enrichment_gsea", cid, "ma_plots", ct))
    tally[[length(tally) + 1]] <- data.table(contrast = cid, analysis = "GSEA", direction = "both",
                                             category = ct, terms_in_csv = nrow(df),
                                             plots = n)
  }
  message("  ", cid, " done")
}
tally <- rbindlist(tally)
fwrite(tally, file.path(OUT, "47b_matched_ma_plot_counts.tsv"), sep = "\t")
message("\nmatched MA plots written:"); print(tally[, .(plots = sum(plots)), by = analysis])

################################################################################
# 2. GO-BP across the two axes
################################################################################
fig50 <- file.path(OUT, "50_gsea_axes_GOBP.png")
unlink(fig50)
read_full <- function(cid, ct) {
  f <- file.path(OUT, "enrichment_gsea", cid, paste0("gsea_full_", ct, ".csv"))
  if (!file.exists(f)) return(NULL)
  d <- as.data.table(read.csv(f)); d[, contrast := cid][]
}
ax <- rbindlist(lapply(PRIMARY_CONTRASTS, read_full, ct = "BP"), fill = TRUE)
if (nrow(ax)) {
  axes  <- PRIMARY_CONTRASTS[1:2]
  n_top <- 20
  n_sig <- sapply(axes, function(cid) ax[contrast == cid & p.adjust < THRESH$gsea_fdr, .N])
  pick  <- unique(unlist(lapply(axes, function(cid)
             ax[contrast == cid & p.adjust < THRESH$gsea_fdr][order(p.adjust)][, head(ID, n_top)])))
  if (length(pick)) {
    hm <- ax[ID %in% pick]
    ordr <- hm[contrast == axes[1]][order(NES), Description]
    hm[, term := factor(Description, levels = unique(c(ordr, setdiff(unique(Description), ordr))))]
    hm[, contrast_lab := factor(CONTRASTS[contrast, label],
                                levels = CONTRASTS[PRIMARY_CONTRASTS, label])]
    hm[, star := fifelse(p.adjust < THRESH$gsea_fdr, "*", "")]
    sub <- sprintf(paste0("Rows: the top %d GO-BP terms by FDR from each axis. Terms at FDR < 0.05: ",
                          "%s %d; %s %d%s. Every tile is coloured by NES whether significant ",
                          "or not; * = FDR < 0.05."),
                   n_top, CONTRASTS[axes[1], label], n_sig[1], CONTRASTS[axes[2], label], n_sig[2],
                   if (n_sig[2] == 0) ", so every row comes from the first axis" else "")
    p <- ggplot(hm, aes(contrast_lab, term, fill = NES)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = star), size = 4, vjust = .78) +
      scale_fill_distiller(palette = "RdBu", limits = c(-3, 3), oob = scales::squish) +
      scale_y_discrete(labels = function(x) str_wrap(x, 55)) +
      labs(x = NULL, y = NULL,
           title = "GO Biological Process GSEA across the within-patient comparisons",
           subtitle = str_wrap(sub, 105)) +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 18, hjust = 1),
            plot.subtitle = element_text(size = 8))
    ggsave(fig50, p, width = 10, height = max(5, length(pick) * 0.3 + 3), dpi = 300,
           limitsize = FALSE)
    message("wrote 50_gsea_axes_GOBP.png: ", length(pick), " terms; FDR<0.05 per axis: ",
            paste(n_sig, collapse = " / "))
  }
}
message("\n== done ==  ", run_stamp())
