#!/usr/bin/env Rscript
###############################################################################
# 23_starburst_programs.R                                            2026-10-19
#
# Starburst plots one program at a time. For each comparison with RNA-seq and
# each program of the Part A program figure (programs_shown.tsv) that changes in
# that comparison (FDR < 0.05 in methylation or in expression,
# 19_clinical_questions.R), every gene is drawn by its promoter methylation change
# (x) and expression change (y), as in the starbursts of 03_promoter_meth_vs_
# expression.R, all genes in grey and the program's genes coloured by quadrant and
# named. The programs whose methylation and expression move in opposite directions
# (methylation lost and expression up, or gained and down: the pattern expected
# when methylation regulates the genes) are marked in the index and the file name.
# Methylation in the program scores is the change of all tiles within 50 kb of the
# genes; the x axis here is the promoter alone (TSS -1,500 to +500).
#
# Output: integration/meth_vs_expr/starburst_programs/<comparison>/<program>.png
#         integration/meth_vs_expr/starburst_programs/index.tsv
###############################################################################
source("config/paths.R")
ID <- with(data.table::fread(SAMPLE_SHEET, colClasses = "character"), setNames(sample, key))   # role key -> sample ID (config/samples.tsv)

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(msigdbr) })
OUT <- file.path(INT, "meth_vs_expr", "starburst_programs")
CQD <- file.path(INT, "clinical")
D_MIN <- 0.10; LFC_MIN <- 0.5                       # the quadrant thresholds of 03_promoter_meth_vs_expression.R
PAIRS <- data.table(rna = c("C1_clone", "C2_transform", "C3_cumulative", "C4_DDHN_vs_donor", "C6_D30_vs_donor"),
                    meth = c(paste0(ID[["Y35"]], "_vs_", ID[["D30"]]), paste0(ID[["DX"]], "_vs_", ID[["Y35"]]), paste0(ID[["DX"]], "_vs_", ID[["D30"]]), paste0(ID[["DX"]], "_vs_", ID[["D23"]]), paste0(ID[["D30"]], "_vs_", ID[["D23"]])),
                    folder = c("3.5_years_vs_1_month", "DDHN_vs_3.5_years", "DDHN_vs_1_month", "DDHN_vs_donor_5_years",
                               "reference_1_month_vs_donor_5_years"))
smp <- fread(SAMPLE_SHEET, colClasses = "character")[methylation == "analysed"]
PAIRS[, title := vapply(meth, function(i) { p <- strsplit(i, "_vs_")[[1]]
  paste(smp[sample == p[1], plot_label], "vs", smp[sample == p[2], plot_label]) }, "")]

long <- fread(file.path(INT, "meth_vs_expr", "promoter_meth_vs_expr_long.tsv"))
prog <- fread(file.path(CQD, "programs.tsv"))
show <- fread(file.path(CQD, "programs_shown.tsv"))
H  <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))[, .(set = gs_name, gene = gene_symbol)]
C8 <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C8"))[grepl("^HAY_BONE_MARROW", gs_name), .(set = gs_name, gene = gene_symbol)]
PN <- fread(file.path(RES, "gene_panels.tsv"))[, .(set = paste0("PANEL_", panel), gene)]
SETS <- unique(rbind(H, C8, PN))
cols <- c(`hypo-up` = "#2166AC", `hyper-down` = "#B2182B", `hypo-down` = "#5E9BC9", `hyper-up` = "#E07B39", other = "grey40")
qname <- c(`hypo-up` = "methylation lost, expression up", `hyper-down` = "methylation gained, expression down",
           `hypo-down` = "methylation lost, expression down", `hyper-up` = "methylation gained, expression up", other = "smaller changes")
slug <- function(x) gsub("_+", "_", gsub("[^A-Za-z0-9+]+", "_", x))

idx <- list()
for (k in seq_len(nrow(PAIRS))) {
  pr <- PAIRS[k]
  d <- long[contrast == pr$rna & !is.na(delta) & is.finite(lfc)]
  pm <- prog[layer == "methylation" & id == pr$meth & set %in% show$set]
  pe <- prog[layer == "RNA-seq" & id == pr$rna & set %in% show$set]
  sets <- union(pm[padj < 0.05, set], pe[padj < 0.05, set])
  for (st in sets) {
    m <- pm[set == st][1]; r <- pe[set == st][1]
    opposite <- !is.na(m$NES) && !is.na(r$NES) && sign(m$NES) != sign(r$NES)
    g <- SETS[set == st, gene]
    dg <- d[mane %in% g | gene %in% g]
    if (!nrow(dg)) next
    dg[, lab := mane]
    lab_d <- dg[class != "other"]
    lbl <- show[set == st, label]
    subt <- sprintf(paste0("%d of the program's genes have promoter methylation and expression; %d lost methylation and rose, %d gained and fell, ",
                          "%d moved with methylation (hypo-down, hyper-up).\nProgram score: methylation NES %s (%s), expression NES %s (%s)%s. ",
                          "Named: the program's genes beyond both dashed lines (|change| >= %d points, |log2FC| >= %.1f)."),
                   nrow(dg), sum(dg$class == "hypo-up"), sum(dg$class == "hyper-down"), sum(dg$class %in% c("hypo-down", "hyper-up")),
                   sprintf("%+.2f", m$NES), if (!is.na(m$padj) && m$padj < 0.05) "FDR < 0.05" else "not significant",
                   sprintf("%+.2f", r$NES), if (!is.na(r$padj) && r$padj < 0.05) "FDR < 0.05" else "not significant",
                   if (opposite) "; methylation and expression move in opposite directions" else "", round(100 * D_MIN), LFC_MIN)
    subt <- paste(subt, "x axis: the promoter (TSS -1,500 to +500); the program score uses all tiles within 50 kb of each gene.")
    p <- ggplot(d, aes(100 * delta, lfc)) +
      geom_hline(yintercept = c(-LFC_MIN, LFC_MIN), colour = "grey75", linetype = 2) +
      geom_vline(xintercept = c(-D_MIN, D_MIN) * 100, colour = "grey75", linetype = 2) +
      geom_point(colour = "grey85", size = 0.5) +
      geom_point(data = dg, aes(colour = class), size = 2) +
      geom_text_repel(data = lab_d, aes(label = lab, colour = class), size = 3.1, max.overlaps = Inf, min.segment.length = 0,
                      segment.size = 0.25, segment.colour = "grey55", box.padding = 0.2, seed = 1, show.legend = FALSE) +
      scale_colour_manual(values = cols, labels = qname, name = NULL, drop = TRUE) +
      labs(x = "Promoter methylation change, first-named minus second-named sample (points)", y = "Expression change (shrunk log2FC)",
           title = sprintf("%s: %s", lbl, pr$title), subtitle = paste(strwrap(gsub("\n", " ", subt), 150), collapse = "\n")) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 9.5))
    f <- file.path(OUT, pr$folder, paste0(if (opposite) "opposite_" else "", slug(sub(" \\[.*", "", lbl)), ".png"))
    dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
    ggsave(f, p, width = 11, height = 9, dpi = 170)
    idx[[length(idx) + 1]] <- data.table(comparison = pr$title, rna_contrast = pr$rna, methylation_comparison = pr$meth, program = lbl,
      methylation_NES = round(m$NES, 2), methylation_FDR = signif(m$padj, 2), expression_NES = round(r$NES, 2), expression_FDR = signif(r$padj, 2),
      opposite_directions = opposite, genes = nrow(dg), hypo_up = sum(dg$class == "hypo-up"), hyper_down = sum(dg$class == "hyper-down"),
      hypo_up_genes = paste(dg[class == "hypo-up"][order(-abs(delta * lfc)), mane], collapse = ", "),
      hyper_down_genes = paste(dg[class == "hyper-down"][order(-abs(delta * lfc)), mane], collapse = ", "),
      file = substring(f, nchar(DATA_ROOT) + 2L))
  }
}
idx <- rbindlist(idx)
fwrite(idx, file.path(OUT, "index.tsv"), sep = "\t")
message(nrow(idx), " program starbursts (", sum(idx$opposite_directions), " with methylation and expression in opposite directions) in ", OUT)
