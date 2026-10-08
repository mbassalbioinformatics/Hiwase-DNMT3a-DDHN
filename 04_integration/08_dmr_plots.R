#!/usr/bin/env Rscript
###############################################################################
# 08_dmr_plots.R                                                     2026-10-08
#
#   volcano/<comparison>.png  every region metilene tested: methylation change
#                             (g1 - g2) against -log10 q; q < 0.05 coloured and
#                             every one named (switch-region gene, else
#                             [nearest gene])
#   heatmap_all_dmrs.png      every metilene DMR from any comparison x every
#                             sample (timeline order), with which comparisons
#                             called it and its genomic context
#   heatmap_<comparison>.png  that comparison's DMRs x every sample
# Regenerated from the current tables on every run.
###############################################################################
source("config/paths.R")
source("04_integration/io_helpers.R")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(ComplexHeatmap); library(circlize)
})
setDTthreads(8)
DMRD <- file.path(EMS, "03_dmr_analysis")
COMP <- COMPARISONS
OUT  <- file.path(INT, "dmr"); dir.create(file.path(OUT, "volcano"), showWarnings = FALSE, recursive = TRUE)
source("04_integration/locus_plot_functions.R")           # LOCUS_SAMPLES
comp <- fread(COMP, colClasses = "character")
m <- fread(cmd = sprintf("gzip -dc '%s'", file.path(OUT, "dmr_master.tsv.gz")))[set == "metilene"]
m[, name := fifelse(!is.na(switch_gene), switch_gene, paste0("[", nearest_gene, "]"))]

## ---- volcano plots ------------------------------------------------------------------
set.seed(1)
for (k in seq_len(nrow(comp))) {
  id <- comp$id[k]
  x <- fread_metilene(DMRD, id, header = FALSE, select = c(1:5),
             col.names = c("chr", "start", "end", "q", "delta"))
  x[, y := -log10(pmax(q, 1e-300))]
  sig <- merge(x[q < 0.05], m[id == comp$id[k], .(chr, start, end, name, context)], by = c("chr", "start", "end"), all.x = TRUE)
  bg <- x[q >= 0.05][sample(.N, min(.N, 150000))]
  p <- ggplot() +
    geom_point(data = bg, aes(100 * delta, y), colour = "grey80", size = 0.25, alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
    geom_point(data = sig, aes(100 * delta, y, colour = delta > 0, shape = context %like% "switch"), size = 2) +
    geom_text_repel(data = sig, aes(100 * delta, y, label = name), size = 2.6, max.overlaps = Inf,
                    min.segment.length = 0, segment.size = 0.2, seed = 1) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"), labels = c(`TRUE` = "gained", `FALSE` = "lost"), name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), labels = c(`TRUE` = "switch region (promoter-exon1-intron1)", `FALSE` = "elsewhere"), name = NULL) +
    labs(x = sprintf("Methylation change, %s - %s (percentage points)", comp$g1[k], comp$g2[k]), y = "-log10 q (metilene)",
         title = sprintf("%s%s: %d DMRs at q < 0.05 (%d gained, %d lost)", comp$label[k],
                         ifelse(comp$status[k] == "provisional", " [provisional]", ""), nrow(sig), sum(sig$delta > 0), sum(sig$delta < 0)),
         subtitle = sprintf("%s regions tested (150,000 shown in grey); DMRs named by switch-region gene, else [nearest gene]",
                            format(nrow(x), big.mark = ","))) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(OUT, "volcano", paste0(id, ".png")), p, width = 11, height = 9, dpi = 250)
}

## ---- heatmaps -----------------------------------------------------------------------
smp <- LOCUS_SAMPLES$sample
lab <- setNames(LOCUS_SAMPLES$label, smp)
u <- unique(m[, c("dmr", "name", "context", "cgi_context", "chr", paste0("meth_", smp)), with = FALSE], by = "dmr")
called <- dcast(m[, .(dmr, id, v = 1L)], dmr ~ id, value.var = "v", fill = 0L)
u <- merge(u, called, by = "dmr")
mat <- as.matrix(u[, paste0("meth_", smp), with = FALSE]); colnames(mat) <- lab[smp]; rownames(mat) <- u$name
cmat <- as.matrix(u[, intersect(comp$id, names(u)), with = FALSE])
ctx_cols <- c(`switch region (promoter-exon1-intron1)` = "#7B3294", `gene body (beyond intron 1)` = "#C2A5CF", intergenic = "#E7D4E8")
draw_hm <- function(rows, file, title) {
  if (length(rows) < 2) return(invisible())
  ha <- rowAnnotation(context = u$context[rows], CGI = u$cgi_context[rows],
                      col = list(context = ctx_cols, CGI = c(CGI = "#1B7837", shore = "#7FBF7B", shelf = "#D9F0D3", `open sea` = "white")),
                      annotation_name_gp = gpar(fontsize = 8))
  hc <- Heatmap(cmat[rows, , drop = FALSE], name = "called in", col = c(`0` = "white", `1` = "black"),
                show_heatmap_legend = FALSE, column_names_gp = gpar(fontsize = 7), column_names_rot = 60,
                width = unit(ncol(cmat) * 3.2, "mm"), cluster_columns = FALSE)
  mm <- mat[rows, , drop = FALSE]
  imp <- mm; for (i in seq_len(nrow(imp))) imp[i, is.na(imp[i, ])] <- mean(imp[i, ], na.rm = TRUE)
  imp[is.na(imp)] <- 0.5                                   # rows with no value at all
  hc_rows <- hclust(dist(imp), method = "average")         # clustering on gaps filled with the row mean
  hm <- Heatmap(mm, name = "methylation", cluster_rows = hc_rows, col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), METH_COLS),
                cluster_columns = FALSE, column_split = factor(LOCUS_SAMPLES$case, levels = unique(LOCUS_SAMPLES$case)),
                show_row_names = length(rows) <= 120, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 8),
                column_names_rot = 40, right_annotation = ha, column_title_gp = gpar(fontsize = 9),
                na_col = "grey88")
  png(file, width = 11, height = max(6, min(26, 0.11 * length(rows) + 3)), units = "in", res = 250)
  draw(hm + hc, column_title = title, merge_legend = TRUE)
  invisible(dev.off())
}
draw_hm(seq_len(nrow(u)), file.path(OUT, "heatmap_all_dmrs.png"),
        sprintf("All %d metilene DMRs (q < 0.05) from any comparison; black = called in that comparison", nrow(u)))
for (i in comp$id) draw_hm(which(u[[i]] == 1L), file.path(OUT, paste0("heatmap_", i, ".png")),
                           sprintf("%s: %d DMRs", comp[id == i, label], sum(u[[i]] == 1L)))
message("written to ", OUT)
