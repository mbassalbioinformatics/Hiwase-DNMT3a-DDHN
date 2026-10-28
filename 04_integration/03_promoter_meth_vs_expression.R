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
#   switch  the promoter changes by >= 0.20 (20 points, the threshold supported
#           by the Mesa work) in either direction, whatever its starting level:
#           each gene turns on or off at its own methylation level, so no fixed
#           level (such as 50%) has to be crossed
# Per comparison, two tests on the switches:
#   direction   among switches whose gene changes expression (|log2FC| >= 0.5),
#               is the change opposite to methylation more often than half the
#               time (binomial)?
#   enrichment  are switch genes more likely to change expression at all than
#               genes whose promoter does not switch (Fisher)?
# Sensitivity: the same two tests with stricter definitions (>= 0.30; >= 0.20
# crossing 0.5; either of the two), to show how much each result depends on
# the definition.
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
#   switch_sensitivity.tsv             the two tests under each switch definition
#   starburst_<contrast>_A_concordant.png   all genes, then the two quadrants where expression moves
#                                           against methylation, zoomed, every gene named
#   starburst_<contrast>_B_discordant.png   the same for the two quadrants where it moves with methylation
###############################################################################
source("config/paths.R")

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })
OUT <- file.path(INT, "meth_vs_expr"); dir.create(OUT, showWarnings = FALSE)
source("03_rnaseq/00_config.R")
D_MIN <- 0.10; LFC_MIN <- 0.5; SW_MIN <- CUT   # CUT: config/settings.R (20 points)

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
long[, switch := fcase(delta >= SW_MIN, "switch on (gain)",
                       delta <= -SW_MIN, "switch off (loss)",
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

## ---- sensitivity to the switch definition -------------------------------------------
# The same two tests with stricter definitions: a larger change, a change that crosses
# 50% (mostly unmethylated <-> mostly methylated), and either of the two.
crosses <- function(d) pmin(d$meth_num, d$meth_den) < 0.5 & pmax(d$meth_num, d$meth_den) >= 0.5
DEFS <- list(`>= 20 points (used)`                  = function(d) abs(d$delta) >= SW_MIN,
             `>= 30 points`                         = function(d) abs(d$delta) >= 0.30,
             `>= 20 points, crossing 50%`           = function(d) abs(d$delta) >= SW_MIN & crosses(d),
             `>= 30 points, or >= 20 crossing 50%`  = function(d) abs(d$delta) >= 0.30 | (abs(d$delta) >= SW_MIN & crosses(d)))
sens <- rbindlist(lapply(CONTRASTS$id, function(cid) {
  d <- long[contrast == cid & !is.na(delta) & is.finite(lfc)]
  if (cid %in% CONTRASTS[class == "context", id]) d <- d[meth_compartment == FALSE]
  rbindlist(lapply(names(DEFS), function(nm) {
    s <- DEFS[[nm]](d)
    chg <- d[s & abs(lfc) >= LFC_MIN]
    conc <- chg[sign(lfc) == -sign(delta)]
    ft <- fisher.test(table(factor(s, c(TRUE, FALSE)), factor(abs(d$lfc) >= LFC_MIN, c(TRUE, FALSE))))
    data.table(contrast = cid, label = CONTRASTS[cid, label], definition = nm, switches = sum(s),
               switches_with_expression_change = nrow(chg), concordant = nrow(conc),
               direction_p = if (nrow(chg)) binom.test(nrow(conc), nrow(chg), 0.5, alternative = "greater")$p.value else NA_real_,
               enrichment_OR = unname(ft$estimate), enrichment_p = ft$p.value)
  }))
}))
stopifnot(isTRUE(all.equal(sens[definition == names(DEFS)[1], enrichment_p], tests$enrichment_p)))   # first definition = the tests above
fwrite(sens, file.path(OUT, "switch_sensitivity.tsv"), sep = "\t")
print(sens[, .(contrast, definition, switches, changed = switches_with_expression_change, concordant,
               direction_p = signif(direction_p, 2), OR = round(enrichment_OR, 2), enrichment_p = signif(enrichment_p, 2))])

## ---- starburst plots ---------------------------------------------------------------
# Per comparison two figures: A, the quadrants where expression moves against
# methylation (hypo-up, hyper-down); B, those where it moves with it (hypo-down,
# hyper-up). Each shows all genes small on the left, then each of its two
# quadrants zoomed, with every gene in the quadrant named at a readable size.
cols  <- c(`hypo-up` = "#2166AC", `hyper-down` = "#B2182B", `hypo-down` = "#5E9BC9", `hyper-up` = "#E07B39", other = "grey82")
qsign <- list(`hypo-up` = c(-1, 1), `hyper-down` = c(1, -1), `hypo-down` = c(-1, -1), `hyper-up` = c(1, 1))
qname <- c(`hypo-up` = "Methylation lost, expression up", `hyper-down` = "Methylation gained, expression down",
           `hypo-down` = "Methylation lost, expression down", `hyper-up` = "Methylation gained, expression up")
xg <- mg[chr == "chrX", gene]
for (cid in CONTRASTS$id) {
  d <- long[contrast == cid & !is.na(delta) & is.finite(lfc)]
  d[, lab := paste0(gene, fifelse(mane %in% xg, " (X)", ""))]
  d[, face := fifelse(sig %in% TRUE, "bold", "plain")]
  is_don <- cid %in% CONTRASTS[class == "context", id]
  s <- tests[contrast == cid]
  xlab <- "Promoter methylation change, first-named minus second-named sample (points)"
  for (panel in c("A_concordant", "B_discordant")) {
    q <- if (panel == "A_concordant") c("hypo-up", "hyper-down") else c("hypo-down", "hyper-up")
    ov <- ggplot(d, aes(delta * 100, lfc)) +
      geom_hline(yintercept = c(-LFC_MIN, LFC_MIN), colour = "grey75", linetype = 2) +
      geom_vline(xintercept = c(-D_MIN, D_MIN) * 100, colour = "grey75", linetype = 2) +
      geom_point(data = d[!class %in% q], colour = "grey82", size = 0.5) +
      geom_point(data = d[class %in% q], aes(colour = class), size = 1.2) +
      scale_colour_manual(values = cols, guide = "none") +
      labs(x = xlab, y = "Expression change (shrunk log2FC)", title = "All genes") +
      theme_bw(base_size = 12)
    zoom <- lapply(q, function(k) {
      dq <- d[class == k]
      # the panel reaches past the thresholds (dashed) so names near the busy corner have room
      xr <- if (qsign[[k]][1] < 0) c(min(dq$delta, -D_MIN) * 100 - 4, -D_MIN * 100 + 5) else c(D_MIN * 100 - 5, max(dq$delta, D_MIN) * 100 + 4)
      yr <- if (qsign[[k]][2] < 0) c(min(dq$lfc, -LFC_MIN) - 0.3, -LFC_MIN + 0.35) else c(LFC_MIN - 0.35, max(dq$lfc, LFC_MIN) + 0.3)
      ggplot(dq, aes(delta * 100, lfc)) +
        geom_hline(yintercept = qsign[[k]][2] * LFC_MIN, colour = "grey70", linetype = 2) +
        geom_vline(xintercept = qsign[[k]][1] * D_MIN * 100, colour = "grey70", linetype = 2) +
        { if (is_don) geom_point(aes(shape = meth_compartment), colour = cols[[k]], size = 2.2)
          else geom_point(colour = cols[[k]], size = 2.2) } +
        geom_text_repel(aes(label = lab, fontface = face), size = 3.4, max.overlaps = Inf, min.segment.length = 0,
                        segment.size = 0.25, segment.colour = "grey55", box.padding = 0.18, point.padding = 0.1,
                        force = 3, force_pull = 0.3, max.time = 5, max.iter = 1e5, seed = 1, xlim = xr, ylim = yr) +
        { if (is_don) scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), drop = FALSE,
                                         name = "Promoter also changes in the reference comparison (Recipient 1 month vs Donor 5 years post donation)") } +
        coord_cartesian(xlim = xr, ylim = yr, expand = FALSE) +
        labs(x = xlab, y = "Expression change (shrunk log2FC)", title = sprintf("%s: %d genes", qname[[k]], nrow(dq))) +
        theme_bw(base_size = 12) + theme(plot.title = element_text(colour = cols[[k]], face = "bold"), legend.position = "bottom")
    })
    fig <- (ov | zoom[[1]] | zoom[[2]]) + plot_layout(widths = c(0.55, 1, 1)) +
      plot_annotation(
        title = sprintf("Promoter methylation against expression: %s", CONTRASTS[cid, label]),
        subtitle = sprintf(paste0("%d genes. Coloured: methylation change >= %d points and |log2FC| >= %.1f. Every gene in the two quadrants where ",
                                  "expression moves %s methylation is named; bold = differentially expressed, (X) = X-linked.\n",
                                  "Promoter switches (change >= %d points): %d; %d change expression, %d against methylation (binomial p %.2g)."),
                           s$genes, round(100 * D_MIN), LFC_MIN, ifelse(panel == "A_concordant", "against", "with"),
                           round(100 * SW_MIN), s$switches, s$switches_with_expression_change, s$concordant, s$direction_p),
        theme = theme(plot.title = element_text(size = 18, face = "bold"), plot.subtitle = element_text(size = 12)))
    ggsave(file.path(OUT, sprintf("starburst_%s_%s.png", cid, panel)), fig, width = 26, height = 12.5, dpi = 200)
  }
}
unlink(file.path(OUT, sprintf("starburst_%s.png", CONTRASTS$id)))          # the earlier single-plot version
message("written to ", OUT)
