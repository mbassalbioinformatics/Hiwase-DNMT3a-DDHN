#!/usr/bin/env Rscript
###############################################################################
# 05_locus_plot_gallery.R                                            2026-10-08
#
# Locus plots (locus_plot_functions.R) for:
#   focus/<panel>/   the RNA-seq report's curated panels (FOCUS_GENES in
#                    rnaseq/scripts/00_config.R): DNA methylation enzymes,
#                    myeloid drivers, splicing factors, p53 targets
#   top/<contrast>/  per comparison with RNA-seq, the genes where promoter
#                    methylation and expression move in opposite directions
#                    (hypo-up, hyper-down; meth_vs_expr/concordant_genes.tsv),
#                    ranked by |methylation change| x |log2FC|, excluding genes
#                    flagged by the reference comparison
#   clone/           genes concordant in BOTH large-clone-vs-donor comparisons
#                    (3.5 yr and DDHN vs donor), same direction, not in the
#                    reference: the methylation counterpart of the RNA-seq's
#                    clone-associated set
# Once the DMRs are in, the DMR-based plots (with the DMR boxed) come from the
# DMR script.
###############################################################################
suppressPackageStartupMessages(library(data.table))
source("config/paths.R")
source("04_integration/locus_plot_functions.R")
source("03_rnaseq/00_config.R")
OUT <- file.path(INT, "locus_plots")
N_TOP <- as.integer(Sys.getenv("N_TOP", "8"))
cache <- Sys.getenv("LOCUS_CACHE")                          # optional saved load_locus_data() result

L <- if (nzchar(cache) && file.exists(cache)) readRDS(cache) else load_locus_data(
  emseq = EMS,
  mane_gtf = MANE,
  cgi_file = file.path(RES, "cpgIslandExt.txt.gz"),
  rna_rdata = file.path(RNA, "05_processed_counts_genes", "dge_gene_level.RData"))
L$min_cov <- as.integer(Sys.getenv("MIN_COV", "1"))   # a CpG is drawn when it has at least this many reads
has_model <- unique(L$model[feature == "transcript", gene])

draw <- function(g, dir, title = NULL) {
  if (!g %in% has_model) return(invisible(FALSE))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  if (Sys.getenv("FORCE", "1") != "1" && file.exists(file.path(dir, paste0(g, ".png")))) return(invisible(TRUE))
  ok <- tryCatch({ ggsave(file.path(dir, paste0(g, ".png")), locus_plot(L, g, title = title),
                          width = 13, height = 5.4, dpi = 220); TRUE },
                 error = function(e) { message("  ", g, ": ", conditionMessage(e)); FALSE })
  invisible(ok)
}

## focus panels
for (pn in names(FOCUS_GENES)) {
  message("focus: ", pn)
  for (g in FOCUS_GENES[[pn]]) draw(g, file.path(OUT, "focus", pn))
}

## top concordant genes per comparison
conc <- fread(file.path(INT, "meth_vs_expr", "concordant_genes.tsv"))
conc <- conc[meth_compartment == FALSE]
for (cid in CONTRASTS$id) {
  top <- conc[contrast == cid][order(-abs(delta * lfc))][seq_len(min(.N, N_TOP))]
  message("top: ", cid, " (", nrow(top), ")")
  for (i in seq_len(nrow(top)))
    draw(top$gene[i], file.path(OUT, "top", cid),
         title = sprintf("%s   [%s: promoter %+.0f points, expression log2FC %+.2f]",
                         top$gene[i], CONTRASTS[cid, label], 100 * top$delta[i], top$lfc[i]))
}

## concordant in both large-clone-vs-donor comparisons, not in the reference
don <- CONTRASTS[class == "context", id]
both <- conc[contrast %in% don][, if (.N == 2 && uniqueN(class) == 1) .SD, by = gene]
cl <- unique(both[, .(gene, class, score = abs(delta * lfc))][order(-score)][, .(gene, class)])
fwrite(dcast(both, gene + class ~ contrast, value.var = c("delta", "lfc")), file.path(OUT, "clone_concordant_genes.tsv"), sep = "\t")
message("clone-concordant genes: ", nrow(cl))
for (i in seq_len(min(nrow(cl), 20))) draw(cl$gene[i], file.path(OUT, "clone"),
                                           title = sprintf("%s   [%s in both large-clone vs donor comparisons]", cl$gene[i], cl$class[i]))
message("written to ", OUT)
