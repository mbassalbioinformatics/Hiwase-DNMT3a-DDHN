#!/usr/bin/env Rscript
###############################################################################
# 25_gene_of_note_clu.R                                              2026-10-19
#
# CLU (clusterin), the gene with the largest expression rise from the recipient
# at 3.5 years to DDHN diagnosis. Megakaryocytes make clusterin and pack it into
# platelet alpha-granules [Tschopp 1993], and bone-marrow RNA carries platelet and
# megakaryocyte transcripts in varying amounts, so CLU is set against six genes
# of that lineage (PF4, PPBP, ITGA2B, GP9, F13A1, TUBB1) in every RNA-seq library:
#   A  expression (CPM, log scale) of CLU and the six genes
#   B  CLU relative to the lineage: CLU CPM / geometric mean of the six (each + 1)
#   C  methylation at the CLU promoter and the tiles within 50 kb of its TSS
# plus a locus plot of the CLU promoter (locus_plot_functions.R).
#
# Output: integration/clinical/clu.tsv, clu_de.tsv, clu_relative.tsv, clu_methylation.tsv, clu_tiles.tsv, fig_clu.png;
#         integration/locus_plots/genes_of_note/CLU.png
###############################################################################
source("config/paths.R")
ID <- with(data.table::fread(SAMPLE_SHEET, colClasses = "character"), setNames(sample, key))   # role key -> sample ID (config/samples.tsv)

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
OUT <- file.path(INT, "clinical")
rd  <- function(f, ...) if (grepl("\\.gz$", f)) fread(cmd = sprintf("gzip -dc '%s'", f), ...) else fread(f, ...)
smp <- fread(SAMPLE_SHEET, colClasses = "character")[methylation == "analysed"]
GENE <- "CLU"; MK <- c("PF4", "PPBP", "ITGA2B", "GP9", "F13A1", "TUBB1")
RNA_S <- c(ID[["D30"]], ID[["Y35"]], ID[["DX"]], ID[["D23"]])
ALL_S <- c(ID[["D18"]], ID[["D30"]], ID[["Y35"]], ID[["DX"]], ID[["D23"]], ID[["T2A"]], ID[["T2DX"]])
short <- function(s) sub(" \\((DNMT3A|TET2)\\)$", "\n(\\1)", smp$plot_label[match(s, smp$sample)])

e <- new.env(); load(file.path(RNA, "05_processed_counts_genes", "dge_gene_level.RData"), envir = e)
cpm <- e$cpm_mat
stopifnot(all(c(GENE, MK) %in% rownames(cpm)), all(RNA_S %in% colnames(cpm)))
ex <- rbindlist(lapply(c(GENE, MK), function(g) data.table(gene = g, sample = RNA_S, cpm = as.numeric(cpm[g, RNA_S]))))
gm <- ex[gene %in% MK, .(lineage = exp(mean(log(cpm + 1)))), by = sample]
rel <- merge(ex[gene == GENE, .(sample, clu = cpm)], gm, by = "sample")[, ratio := (clu + 1) / lineage]
ed <- rd(file.path(RNA, "05_processed_counts_genes", "09_edgeR_all_bcv.tsv.gz"))[bcv == 0.4 & symbol %in% c(GENE, MK)]
fwrite(merge(dcast(ex, gene ~ sample, value.var = "cpm"),
             dcast(ed[, .(gene = symbol, contrast, lfc = round(lfc_shrunk, 2))], gene ~ contrast, value.var = "lfc"), by = "gene"),
       file.path(OUT, "clu.tsv"), sep = "\t")

# methylation: promoter and the tiles within 50 kb of the TSS
pm <- rd(file.path(INT, "region_methylation", "promoter_methylation.tsv.gz"))[id == GENE]
gn <- rd(file.path(INT, "region_methylation", "mane_select_genes.tsv.gz"))[gene == GENE][1]
ti <- rd(TILE_FILE)[grepl(paste0("^", gn$chr, ":"), id)]
ti[, start := as.integer(sub(".*:", "", id))][, mid := start + TILE_BP %/% 2L - 1L]
ti <- ti[abs(mid - gn$tss) <= 50e3]
ti[, kb := (mid - gn$tss) * fifelse(gn$strand == "+", 1, -1) / 1e3]
fwrite(ti[, c("id", "kb", ALL_S), with = FALSE], file.path(OUT, "clu_tiles.tsv"), sep = "\t")
ml <- rbind(data.table(region = "promoter (TSS -1,500 to +500)", sample = ALL_S, meth = as.numeric(pm[, ..ALL_S])),
            data.table(region = sprintf("mean of the %d tiles within 50 kb", nrow(ti)), sample = ALL_S,
                       meth = vapply(ALL_S, function(s) mean(ti[[s]], na.rm = TRUE), 0)))
ml[, lab := factor(short(sample), levels = short(ALL_S))]
fwrite(ml[, .(region, sample, methylation = round(meth, 4))], file.path(OUT, "clu_methylation.tsv"), sep = "\t")
fwrite(ed[symbol == GENE, .(gene = symbol, contrast, logFC = round(logFC, 3), log2FC_shrunk = round(lfc_shrunk, 3), FDR = signif(FDR, 3),
                            cpm_first = round(cpm_num, 2), cpm_second = round(cpm_den, 2))], file.path(OUT, "clu_de.tsv"), sep = "\t")
fwrite(rel[, .(sample, clu_cpm = round(clu, 2), lineage_geomean = round(lineage, 3), ratio = round(ratio, 3))], file.path(OUT, "clu_relative.tsv"), sep = "\t")

lv <- short(RNA_S)
ex[, `:=`(lab = factor(short(sample), levels = lv), clu = gene == GENE)]
pA <- ggplot(ex, aes(lab, cpm, group = gene)) +
  geom_line(data = ex[clu == FALSE], colour = "grey60", linewidth = 0.6) +
  geom_point(data = ex[clu == FALSE], colour = "grey60", size = 1.6) +
  geom_text(data = ex[clu == FALSE & sample == ID[["D23"]]], aes(label = gene), hjust = -0.25, size = 3, colour = "grey40") +
  geom_line(data = ex[clu == TRUE], colour = "#C0392B", linewidth = 1.4) +
  geom_point(data = ex[clu == TRUE], colour = "#C0392B", size = 3) +
  geom_text(data = ex[clu == TRUE & sample == ID[["D23"]]], aes(label = gene), hjust = -0.3, size = 3.8, colour = "#C0392B", fontface = "bold") +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10), breaks = c(0, 1, 3, 10, 30, 100, 300)) +
  scale_x_discrete(expand = expansion(add = c(0.3, 0.9))) +
  labs(x = NULL, y = "CPM (log scale near 0 made linear)", title = "A. CLU and six platelet / megakaryocyte genes",
       subtitle = "the 3.5-year library has almost none of the lineage's RNA, and almost no CLU") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 9.5))
rel[, lab := factor(short(sample), levels = lv)]
pB <- ggplot(rel, aes(lab, ratio)) + geom_col(fill = "#C0392B", width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", ratio)), vjust = -0.4, size = 3.4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "CLU / platelet-megakaryocyte genes", title = "B. CLU relative to the lineage",
       subtitle = "(CLU CPM + 1) / geometric mean of the six genes (CPM + 1)") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 9.5))
pC <- ggplot(ml, aes(lab, 100 * meth, fill = region)) + geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = c("#3F6E9A", "#9DB8D3"), name = NULL) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "methylation (%)", title = "C. CLU methylation in every library") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), axis.text.x = element_text(size = 8.5))
ggsave(file.path(OUT, "fig_clu.png"), (pA | pB) / pC + plot_layout(heights = c(1, 0.9)), width = 13, height = 10, dpi = 170)

# a locus plot of the CLU promoter, all libraries, with its expression
source("04_integration/locus_plot_functions.R")
L <- load_locus_data(emseq = EMS, mane_gtf = MANE, cgi_file = file.path(RES, "cpgIslandExt.txt.gz"),
                     rna_rdata = file.path(RNA, "05_processed_counts_genes", "dge_gene_level.RData"))
dir.create(file.path(INT, "locus_plots", "genes_of_note"), showWarnings = FALSE, recursive = TRUE)
p <- locus_plot(L, GENE, title = "CLU (clusterin)", subtitle = "Methylation per CpG (TAPS) at the CLU promoter")
ggsave(file.path(INT, "locus_plots", "genes_of_note", "CLU.png"), p, width = 14, height = attr(p, "height_in"), dpi = 170)
print(rel); print(ml)
message("written to ", OUT)
