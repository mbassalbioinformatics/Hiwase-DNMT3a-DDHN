#!/usr/bin/env Rscript
###############################################################################
# 03_promoter_meth_vs_expression.R                                   2026-10-08
#
# Promoter methylation change against expression change, gene by gene, for the
# six DNMT3A-case comparisons that have RNA-seq (C1-C6; same samples and same
# direction: first-named minus second). Every gene with a promoter value is
# used, not only DMRs or DE genes, so subtle expression changes count.
#
# Methylation acts as a switch, not a dial (the Mesa paper), so nothing here is
# a correlation. Per gene:
#   class   |delta| >= 0.10 and |log2FC| >= 0.5, by quadrant: hypo-up and
#           hyper-down (concordant with methylation silencing), hypo-down and
#           hyper-up (discordant)
#   switch  the promoter changes by >= 0.20 AND crosses 0.5, i.e. flips between
#           mostly unmethylated and mostly methylated
# Per comparison, two tests on the switches:
#   direction   among switches whose gene changes expression (|log2FC| >= 0.5),
#               is the change opposite to methylation more often than half the
#               time (binomial)?
#   enrichment  are switch genes more likely to change expression at all than
#               genes whose promoter does not switch (Fisher)?
# Reference: in the donor comparisons (C4, C5) a gene is flagged "compartment"
# when day 30 vs donor (C6) changes its promoter the same way by >= 0.10.
#
# Promoter = MANE Select TSS -1500/+500; methylation = mean of CpGs with >= 5
# reads (01_region_methylation.R). Expression = shrunk log2FC.
#
# Output in integration/meth_vs_expr/:
#   promoter_meth_vs_expr_long.tsv     every gene x comparison
#   concordant_genes.tsv               hypo-up and hyper-down genes
#   switch_tests.tsv                   per comparison
#   starburst_<contrast>_A_concordant.png   every gene in the two concordant quadrants named
#   starburst_<contrast>_B_discordant.png   every gene in the two discordant quadrants named
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
source("config/paths.R")
OUT <- file.path(INT, "meth_vs_expr"); dir.create(OUT, showWarnings = FALSE)
source("03_rnaseq/00_config.R")
D_MIN <- 0.10; LFC_MIN <- 0.5; SW_MIN <- 0.20

pm <- fread(cmd = sprintf("gzip -dc '%s'", file.path(INT, "region_methylation", "promoter_methylation.tsv.gz")))
mg <- fread(cmd = sprintf("gzip -dc '%s'", file.path(INT, "region_methylation", "mane_select_genes.tsv.gz")))
e <- new.env(); load(file.path(RNA, "05_processed_counts_genes", "dge_gene_level.RData"), envir = e)
de <- as.data.table(e$de)

# RNA genes to MANE genes: Ensembl ID first, then symbol
de[, ens := sub("\\..*", "", ensembl)]
map <- rbind(merge(unique(de[, .(gene, ens)]), mg[, .(ens = ensembl, mane = gene)], by = "ens")[, .(gene, mane)],
             merge(unique(de[, .(gene, symbol)]), mg[, .(symbol = gene, mane = gene)], by = "symbol")[, .(gene, mane)])
map <- map[!duplicated(gene)]
message(uniqueN(map$gene), " of ", uniqueN(de$gene), " RNA-seq genes have a MANE promoter")

long <- rbindlist(lapply(CONTRASTS$id, function(cid) {
  num <- CONTRASTS[cid, num]; den <- CONTRASTS[cid, den]
  m <- pm[, .(mane = id, meth_num = get(num), meth_den = get(den))]
  r <- de[contrast == cid, .(gene, lfc = lfc_shrunk, t, FDR, sig, on_off, rna_ruler_flagged = ruler_flagged,
                             cpm_num, cpm_den)]
  d <- merge(merge(r, map, by = "gene"), m, by = "mane")
  d[, `:=`(contrast = cid, num = num, den = den, delta = meth_num - meth_den)]
  d
}))
long[, class := fcase(delta <= -D_MIN & lfc >= LFC_MIN, "hypo-up",
                      delta >=  D_MIN & lfc <= -LFC_MIN, "hyper-down",
                      delta <= -D_MIN & lfc <= -LFC_MIN, "hypo-down",
                      delta >=  D_MIN & lfc >= LFC_MIN, "hyper-up",
                      default = "other")]
long[, switch := fcase(delta >= SW_MIN & meth_den < 0.5 & meth_num >= 0.5, "switch on (gain)",
                       delta <= -SW_MIN & meth_den >= 0.5 & meth_num < 0.5, "switch off (loss)",
                       default = "none")]
ref <- long[contrast == RULER_CONTRAST, .(gene, ref_delta = delta)]
long <- merge(long, ref, by = "gene", all.x = TRUE)
long[, meth_compartment := contrast %in% CONTRASTS[class == "context", id] &
                           !is.na(delta) & !is.na(ref_delta) &
                           abs(ref_delta) >= D_MIN & sign(ref_delta) == sign(delta)]
long[, label := CONTRASTS[match(contrast, id), label]]
long[, chrX := mane %in% mg[chr == "chrX", gene]]
setcolorder(long, c("gene", "mane", "contrast", "label", "num", "den", "meth_num", "meth_den", "delta",
                    "lfc", "t", "FDR", "sig", "class", "switch", "meth_compartment", "rna_ruler_flagged", "chrX"))
fwrite(long, file.path(OUT, "promoter_meth_vs_expr_long.tsv"), sep = "\t", na = "NA")
fwrite(long[class %in% c("hypo-up", "hyper-down")][order(contrast, -abs(delta * lfc))],
       file.path(OUT, "concordant_genes.tsv"), sep = "\t", na = "NA")

## ---- switch tests -------------------------------------------------------------------
tests <- rbindlist(lapply(CONTRASTS$id, function(cid) {
  d <- long[contrast == cid & !is.na(delta) & is.finite(lfc)]
  if (cid %in% CONTRASTS[class == "context", id]) d <- d[meth_compartment == FALSE]
  sw <- d[switch != "none"]
  chg <- sw[abs(lfc) >= LFC_MIN]
  conc <- chg[(switch == "switch on (gain)" & lfc < 0) | (switch == "switch off (loss)" & lfc > 0)]
  ft <- fisher.test(table(factor(d$switch != "none", c(TRUE, FALSE)), factor(abs(d$lfc) >= LFC_MIN, c(TRUE, FALSE))))
  data.table(contrast = cid, label = CONTRASTS[cid, label], genes = nrow(d),
             switches = nrow(sw), on = sum(sw$switch == "switch on (gain)"), off = sum(sw$switch == "switch off (loss)"),
             switches_with_expression_change = nrow(chg), concordant = nrow(conc),
             direction_p = if (nrow(chg)) binom.test(nrow(conc), nrow(chg), 0.5, alternative = "greater")$p.value else NA,
             pct_changed_switch = round(100 * mean(abs(sw$lfc) >= LFC_MIN), 1),
             pct_changed_other = round(100 * mean(abs(d[switch == "none"]$lfc) >= LFC_MIN), 1),
             enrichment_OR = unname(ft$estimate), enrichment_p = ft$p.value,
             concordant_genes = paste(conc[order(-abs(delta))]$gene, collapse = ", "))
}))
fwrite(tests, file.path(OUT, "switch_tests.tsv"), sep = "\t")
print(tests[, .(label, genes, switches, on, off, changed = switches_with_expression_change, concordant,
                direction_p = signif(direction_p, 2), pct_changed_switch, pct_changed_other,
                OR = round(enrichment_OR, 2), enrichment_p = signif(enrichment_p, 2))])

## ---- starburst plots: every gene in the coloured quadrants named ---------------------
cols <- c(`hypo-up` = "#2166AC", `hyper-down` = "#B2182B", `hypo-down` = "#7FB7D8", `hyper-up` = "#F09A6B", other = "grey82")
xg <- mg[chr == "chrX", gene]
for (cid in CONTRASTS$id) {
  d <- long[contrast == cid & !is.na(delta) & is.finite(lfc)]
  is_don <- cid %in% CONTRASTS[class == "context", id]
  s <- tests[contrast == cid]
  for (panel in c("A_concordant", "B_discordant")) {
    q <- if (panel == "A_concordant") c("hypo-up", "hyper-down") else c("hypo-down", "hyper-up")
    d[, lab := fifelse(class %in% q, paste0(gene, fifelse(mane %in% xg, " (X)", "")), "")]
    d[, face := fifelse(sig %in% TRUE, "bold", "plain")]
    p <- ggplot(d, aes(delta * 100, lfc)) +
      geom_hline(yintercept = c(-LFC_MIN, LFC_MIN), colour = "grey85", linetype = 2) +
      geom_vline(xintercept = c(-D_MIN, D_MIN) * 100, colour = "grey85", linetype = 2) +
      geom_point(data = d[class == "other"], colour = "grey82", size = 0.6) +
      { if (is_don) geom_point(data = d[class != "other"], aes(colour = class, shape = meth_compartment), size = 1.5)
        else geom_point(data = d[class != "other"], aes(colour = class), size = 1.5) } +
      geom_text_repel(data = d[lab != ""], aes(label = lab, fontface = face), size = 2.3, max.overlaps = Inf,
                      min.segment.length = 0, segment.size = 0.2, segment.colour = "grey45", box.padding = 0.18,
                      force = 2, seed = 1) +
      scale_colour_manual(values = cols, name = NULL) +
      { if (is_don) scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), name = "also in reference\n(day 30 vs donor)") } +
      labs(x = sprintf("Promoter methylation change, %s - %s (percentage points)", CONTRASTS[cid, num], CONTRASTS[cid, den]),
           y = "Expression change (shrunk log2FC)",
           title = sprintf("Promoter methylation vs expression: %s  (%s quadrants named)", CONTRASTS[cid, label],
                           ifelse(panel == "A_concordant", "concordant", "discordant")),
           subtitle = sprintf(paste0("%d genes. Promoter switches (change >= %d points crossing 50%%): %d; with an expression change %d, ",
                                     "of which %d opposite to methylation (binomial p %.2g). Bold = DE in the RNA-seq; (X) = X-linked"),
                              s$genes, round(100 * SW_MIN), s$switches, s$switches_with_expression_change, s$concordant, s$direction_p)) +
      theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 8.5))
    ggsave(file.path(OUT, sprintf("starburst_%s_%s.png", cid, panel)), p, width = 15, height = 11, dpi = 250)
  }
}
unlink(file.path(OUT, sprintf("starburst_%s.png", CONTRASTS$id)))          # the earlier single-plot version
message("written to ", OUT)
