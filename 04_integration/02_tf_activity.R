#!/usr/bin/env Rscript
###############################################################################
# 02_tf_activity.R                                                   2026-10-08
#
# Transcription-factor activity from the RNA-seq contrasts: decoupleR's
# univariate linear model (ulm) on each contrast's gene statistic, against the
# CollecTRI TF-target network. A TF's activity is read from the coordinated
# movement of its known targets, with sign (activation or repression), so it can
# change while the TF's own mRNA does not.
#
# Gene statistic: the edgeR t-like statistic in `t`, the same ranking the
# RNA-seq GSEA used.
# Network: CollecTRI, downloaded from the OmniPath web service on 8 Oct 2026
# (OmnipathR 4.0.0's own loaders failed) and processed as
# decoupleR::get_collectri(split_complexes = FALSE) does: JUN/FOS complexes ->
# AP1, REL/NFKB -> NFKB, mor = +1 for stimulation and -1 otherwise. A copy is
# kept in resources/ so the scores can be reproduced.
#
# The reference comparison applies as for pathways (RNA-seq 05_pathway_ruler.R):
# in a donor contrast (C4, C5) a TF is "compartment" if the reference (C6, day
# 30 vs donor) moves it the same way at p < 0.10, and "ruler-clean" otherwise.
#
# Output in integration/tf_activity/:
#   tf_activity_long.tsv    TF x contrast: score, p, FDR, n targets, ruler verdict
#   tf_activity_scores.tsv  TF x contrast scores
#   tf_activity_heatmap.png TFs at FDR < 0.05 in at least one contrast
#   tf_mrna_vs_activity_<contrast>.png  TF mRNA change against TF activity
###############################################################################
suppressPackageStartupMessages({
  library(data.table); library(decoupleR); library(ggplot2); library(ggrepel)
  library(ComplexHeatmap); library(circlize)
})
source("config/paths.R")
OUT  <- file.path(INT, "tf_activity"); dir.create(OUT, showWarnings = FALSE)
RAW  <- Sys.getenv("COLLECTRI_RAW")           # OmniPath interactions TSV (raw)
source("03_rnaseq/00_config.R")

## ---- network ------------------------------------------------------------------
net_file <- "resources/collectri_human_2026-10-08.tsv"         # pinned copy, in the repository
if (!file.exists(net_file)) {
  r <- fread(RAW)
  is_cx <- grepl("COMPLEX", r$source)
  r[is_cx, source_genesymbol := fifelse(grepl("JUN|FOS", source_genesymbol), "AP1",
                                fifelse(grepl("REL|NFKB", source_genesymbol), "NFKB", NA_character_))]
  r <- r[!is.na(source_genesymbol)]
  r <- r[!duplicated(r[, .(source_genesymbol, target_genesymbol)])]
  net <- r[, .(source = source_genesymbol, target = target_genesymbol,
               mor = fifelse(is_stimulation %in% c(TRUE, "True", 1), 1, -1))]
  fwrite(net, net_file, sep = "\t")
}
net <- fread(net_file)
message("CollecTRI: ", nrow(net), " edges, ", uniqueN(net$source), " TFs")

## ---- gene statistics ----------------------------------------------------------
e <- new.env(); load(file.path(RNA, "05_processed_counts_genes", "dge_gene_level.RData"), envir = e)
de <- as.data.table(e$de)
mat <- dcast(de[!is.na(symbol) & is.finite(t)], symbol ~ contrast, value.var = "t", fun.aggregate = function(v) v[1])
genes <- mat$symbol; mat <- as.matrix(mat[, -1]); rownames(mat) <- genes
mat <- mat[, CONTRASTS$id]

## ---- activity -----------------------------------------------------------------
res <- as.data.table(run_ulm(mat, net, .source = "source", .target = "target", .mor = "mor", minsize = 5))
setnames(res, c("source", "condition"), c("tf", "contrast"))
ntg <- net[target %in% genes, .N, by = .(tf = source)]
res <- merge(res[statistic == "ulm", .(tf, contrast, score, p = p_value)], ntg, by = "tf")
res[, fdr := p.adjust(p, "BH"), by = contrast]

rul <- res[contrast == RULER_CONTRAST, .(tf, ruler_score = score, ruler_p = p)]
res <- merge(res, rul, by = "tf")
don <- CONTRASTS[class == "context", id]
res[, verdict := fifelse(contrast %in% don,
         fifelse(ruler_p < THRESH$ruler_p & sign(ruler_score) == sign(score), "compartment", "ruler-clean"),
         NA_character_)]
res[, label := CONTRASTS[match(contrast, id), label]]
setorder(res, contrast, p)
fwrite(res, file.path(OUT, "tf_activity_long.tsv"), sep = "\t")
sc <- dcast(res, tf ~ contrast, value.var = "score")[, c("tf", CONTRASTS$id), with = FALSE]
fwrite(sc, file.path(OUT, "tf_activity_scores.tsv"), sep = "\t")

## ---- summary ------------------------------------------------------------------
message("\nTFs at FDR < 0.05 per contrast (up / down):")
print(res[fdr < 0.05, .(up = sum(score > 0), down = sum(score < 0)), by = .(contrast, label)])
message("\nRuler-clean TFs at FDR < 0.05 in BOTH donor contrasts, same direction:")
both <- res[contrast %in% don & fdr < 0.05 & verdict == "ruler-clean"][, if (.N == 2 && length(unique(sign(score))) == 1) .SD, by = tf]
print(dcast(both, tf ~ contrast, value.var = "score")[order(-abs(get(don[1])))])

## ---- heatmap ------------------------------------------------------------------
keep <- res[fdr < 0.05, unique(tf)]
top <- res[tf %in% keep, .(m = max(abs(score))), by = tf][order(-m)][seq_len(min(.N, 60)), tf]
hm <- as.matrix(sc[tf %in% top, -1]); rownames(hm) <- sc[tf %in% top, tf]
colnames(hm) <- CONTRASTS[match(colnames(hm), id), label]
pm <- as.matrix(dcast(res[tf %in% top], tf ~ contrast, value.var = "fdr")[, c("tf", CONTRASTS$id), with = FALSE][, -1])
rownames(pm) <- dcast(res[tf %in% top], tf ~ contrast, value.var = "fdr")$tf
pm <- pm[rownames(hm), , drop = FALSE]
lim <- quantile(abs(hm), 0.98, na.rm = TRUE)
png(file.path(OUT, "tf_activity_heatmap.png"), width = 9, height = max(6, 0.22 * nrow(hm) + 2.5), units = "in", res = 300)
draw(Heatmap(hm, name = "TF activity\n(ulm score)",
             col = colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B")),
             cluster_columns = FALSE, row_names_gp = gpar(fontsize = 8), column_names_rot = 30,
             column_names_gp = gpar(fontsize = 9),
             cell_fun = function(j, i, x, y, w, h, fill) if (!is.na(pm[i, j]) && pm[i, j] < 0.05)
               grid.text("*", x, y, gp = gpar(fontsize = 9)),
             column_title = "TF activity from target genes (CollecTRI, decoupleR ulm); * FDR < 0.05"))
invisible(dev.off())

## ---- TF mRNA against TF activity ----------------------------------------------
focus <- c("TP53", "MYC", "E2F1", "STAT1", "IRF1", "IRF7", "GATA1", "KLF1", "SPI1", "CEBPA",
           "RUNX1", "MEIS1", "HOXA9", "NFKB", "AP1", "TAL1", "GFI1", "MYB", "FOXO1", "HIF1A")
for (cid in c(don, "C2_transform")) {
  d <- merge(res[contrast == cid, .(tf, score, fdr)],
             de[contrast == cid, .(tf = symbol, lfc = lfc_shrunk)], by = "tf")
  d[, lab := fifelse(tf %in% focus | (fdr < 0.01 & abs(score) > quantile(abs(score), 0.97)), tf, "")]
  p <- ggplot(d, aes(lfc, score)) +
    geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
    geom_point(aes(colour = fdr < 0.05), size = 1.4, alpha = 0.7) +
    geom_text_repel(aes(label = lab), size = 3, max.overlaps = 40, min.segment.length = 0) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey60"), name = "activity FDR < 0.05") +
    labs(x = "TF mRNA change (shrunk log2FC)", y = "TF activity (ulm score)",
         title = paste0("TF mRNA vs activity: ", CONTRASTS[cid, label]),
         subtitle = "Activity is read from the TF's targets; it can move while the TF's own mRNA does not") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT, sprintf("tf_mrna_vs_activity_%s.png", cid)), p, width = 8, height = 6.5, dpi = 300)
}

## ---- TP53: directly activated vs repressed targets ---------------------------------
# CollecTRI's TP53 regulon includes cell-cycle genes p53 represses (via p21-DREAM).
# Where proliferation differs (marrow vs blood) they dominate the TP53 score, so the
# score is also computed on each half of the regulon separately.
tp <- net[source == "TP53"]
split_run <- function(n, lab) as.data.table(run_ulm(mat, n, .source = "source", .target = "target", .mor = "mor", minsize = 5))[
  statistic == "ulm" & source == "TP53", .(targets = lab, contrast = condition, score, p = p_value)]
tps <- rbind(split_run(net, sprintf("all (%d)", nrow(tp))),
             split_run(rbind(net[source != "TP53"], tp[mor > 0]), sprintf("activated only (%d)", sum(tp$mor > 0))),
             split_run(rbind(net[source != "TP53"], tp[mor < 0]), sprintf("repressed only (%d)", sum(tp$mor < 0))))
tps[, label := CONTRASTS[match(contrast, id), label]]
fwrite(tps, file.path(OUT, "tp53_activated_vs_repressed.tsv"), sep = "\t")
print(dcast(tps, targets ~ contrast, value.var = "score"))

message("written to ", OUT)
