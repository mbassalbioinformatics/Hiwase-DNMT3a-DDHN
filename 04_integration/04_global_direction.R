#!/usr/bin/env Rscript
###############################################################################
# 04_global_direction.R                                              2026-10-08
#
# Which way does methylation move, and where, in each comparison? Uses 1-kb
# tiles (01_region_methylation.R) annotated by:
#   CpG-island context  CGI / shore (<= 2 kb) / shelf (2-4 kb) / open sea
#                       (UCSC cpgIslandExt, hg38)
#   ENCODE cCRE         PLS (promoter-like) > pELS (proximal enhancer) > dELS
#                       (distal enhancer) > CTCF-only > DNase-H3K4me3 > none
#                       (UCSC encodeCcreCombined, cell-type agnostic)
# A tile counts as lost or gained at |delta| >= 0.25. Single samples carry
# sampling noise, but noise is symmetric, so the balance of losses against
# gains (net) is the robust read-out: DNMT3A loss is expected to tip it towards
# losses, TET2 loss towards gains, especially at enhancers.
#
# Also: how alike the comparisons are, as Spearman correlations of their tile
# changes, over tiles that change by >= 0.10 in at least one comparison.
#
# Output in integration/global_direction/:
#   direction_by_context.tsv, direction_net.png, delta_correlation.png/.tsv
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ComplexHeatmap); library(circlize) })
setDTthreads(8)
source("config/paths.R")
COMP <- COMPARISONS
OUT  <- file.path(INT, "global_direction"); dir.create(OUT, showWarnings = FALSE)
D_CUT <- 0.25

tw <- fread(cmd = sprintf("gzip -dc '%s'", file.path(INT, "region_methylation", "tiles_1kb_methylation.tsv.gz")))
tw[, c("chr", "start") := tstrsplit(id, ":", fixed = TRUE)]
tw[, `:=`(start = as.integer(start))][, end := start + 999L]
comp <- fread(COMP, colClasses = "character")

## ---- annotation -----------------------------------------------------------------
cgi <- fread(cmd = sprintf("gzip -dc '%s'", file.path(RES, "cpgIslandExt.txt.gz")),
             select = 2:4, col.names = c("chr", "start", "end"))
cgi[, start := start + 1L]                                    # UCSC 0-based -> 1-based
ov_any <- function(a, b) {                                    # TRUE where tile a overlaps any interval of b
  setkey(b, chr, start, end)
  hit <- foverlaps(a[, .(id, chr, start, end)], b, nomatch = 0L, which = TRUE)
  a$id %in% a$id[unique(hit$xid)]
}
flank <- function(d, w0, w1) rbind(d[, .(chr, start = pmax(1L, start - w1), end = start - w0 - 1L)],
                                   d[, .(chr, start = end + w0 + 1L, end = end + w1)])[end >= start]
in_cgi   <- ov_any(tw, copy(cgi))
in_shore <- ov_any(tw, flank(cgi, 0L, 2000L))
in_shelf <- ov_any(tw, flank(cgi, 2000L, 4000L))
tw[, cgi_context := fifelse(in_cgi, "CGI", fifelse(in_shore, "shore", fifelse(in_shelf, "shelf", "open sea")))]

cc <- fread(file.path(RES, "encode_cCREs_hg38_ucsc.bed"), select = c(1, 2, 3, 11),
            col.names = c("chr", "start", "end", "label"))
cc[, start := start + 1L]
pri <- c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3")
tw[, ccre := "none"]
for (lab in rev(pri)) tw[ov_any(tw, cc[label == lab, .(chr, start, end)]), ccre := lab]
tw[, ccre := factor(ccre, levels = c(pri, "none"))]
message("tiles: ", nrow(tw), "; ", paste(names(table(tw$cgi_context)), table(tw$cgi_context), sep = " ", collapse = ", "))

## ---- direction per comparison and context ---------------------------------------
dir_tab <- rbindlist(lapply(seq_len(nrow(comp)), function(i) {
  g1 <- comp$g1[i]; g2 <- comp$g2[i]
  d <- tw[!is.na(get(g1)) & !is.na(get(g2)), .(cgi_context, ccre, delta = get(g1) - get(g2))]
  one <- function(x, ctx_type, ctx) data.table(
    id = comp$id[i], status = comp$status[i], label = comp$label[i], case = comp$case[i],
    context_type = ctx_type, context = ctx, tiles = nrow(x),
    lost = mean(x$delta <= -D_CUT), gained = mean(x$delta >= D_CUT), median_delta = median(x$delta))
  rbind(one(d, "all", "all"),
        rbindlist(lapply(c("CGI", "shore", "shelf", "open sea"), function(k) one(d[cgi_context == k], "CpG island", k))),
        rbindlist(lapply(c(pri, "none"), function(k) one(d[ccre == k], "ENCODE cCRE", k))))
}))
dir_tab[, net := gained - lost]
fwrite(dir_tab, file.path(OUT, "direction_by_context.tsv"), sep = "\t")
print(dir_tab[context == "all", .(id, status, tiles, lost = round(100 * lost, 2), gained = round(100 * gained, 2),
                                   net_pct = round(100 * net, 2), median_delta = round(median_delta, 4))])

d <- dir_tab[context_type != "all"]
d[, comparison := factor(paste0(id, ifelse(status == "provisional", " *", "")),
                         levels = rev(paste0(comp$id, ifelse(comp$status == "provisional", " *", ""))))]
d[, context := factor(context, levels = c("CGI", "shore", "shelf", "open sea", pri, "none"))]
p <- ggplot(d, aes(100 * net, comparison, fill = net > 0)) +
  geom_col() + geom_vline(xintercept = 0, colour = "grey40") +
  facet_wrap(~ context, nrow = 2) +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"), labels = c(`TRUE` = "more gains", `FALSE` = "more losses"), name = NULL) +
  labs(x = sprintf("Tiles gaining minus tiles losing methylation (%% of tiles; |change| >= %.2f)", D_CUT), y = NULL,
       title = "Direction of methylation change by genomic context",
       subtitle = "first-named sample minus second; * provisional comparison") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "direction_net.png"), p, width = 14, height = 7.5, dpi = 300)

## ---- how alike are the comparisons? --------------------------------------------
smp <- unique(c(comp$g1, comp$g2))
full <- tw[complete.cases(tw[, ..smp])]
dm <- sapply(seq_len(nrow(comp)), function(i) full[[comp$g1[i]]] - full[[comp$g2[i]]])
colnames(dm) <- comp$id
var_rows <- rowSums(abs(dm) >= 0.10) > 0
cm <- cor(dm[var_rows, ], method = "spearman")
fwrite(as.data.table(cm, keep.rownames = "id"), file.path(OUT, "delta_correlation.tsv"), sep = "\t")
png(file.path(OUT, "delta_correlation.png"), width = 8.5, height = 7.5, units = "in", res = 300)
draw(Heatmap(cm, name = "Spearman rho", col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
             cell_fun = function(j, i, x, y, w, h, fill) grid.text(sprintf("%.2f", cm[i, j]), x, y, gp = gpar(fontsize = 7)),
             row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
             column_title = sprintf("How alike the methylation changes are (%s variable 1-kb tiles)", format(sum(var_rows), big.mark = ","))))
invisible(dev.off())

## ---- do the TET2 case's changes fall where the DNMT3A clone's do? -------------------
# Tiles changing >= D_CUT in the TET2 comparison against each DNMT3A-case comparison
# (Fisher, all tiles with data in all 7 samples). Two DNMT3A-case comparisons that
# share no sample give the ceiling for a genuinely shared program.
full7 <- tw[complete.cases(tw[, ..smp])]
dd <- function(a, b) full7[[a]] - full7[[b]]
# query (overlap 0) and targets (overlap 1, 2, ...) from the comparisons table
has <- function(x) !is.na(x) & nzchar(x)
ov  <- comp[has(overlap)][order(as.integer(overlap))]
pairs <- setNames(lapply(seq_len(nrow(ov)), function(i) c(ov$g1[i], ov$g2[i])), ov$short)
ceil <- comp[has(ceiling)][order(as.integer(ceiling)), short]
D <- sapply(pairs, function(p) dd(p[1], p[2]))
ov_test <- function(a, b, dir) {
  x <- if (dir == "lost") D[, a] <= -D_CUT else D[, a] >= D_CUT
  y <- if (dir == "lost") D[, b] <= -D_CUT else D[, b] >= D_CUT
  ft <- fisher.test(table(factor(x, c(TRUE, FALSE)), factor(y, c(TRUE, FALSE))))
  data.table(a = a, b = b, direction = dir, tiles_a = sum(x), tiles_b = sum(y), shared = sum(x & y),
             expected = sum(x) * sum(y) / length(x), odds_ratio = unname(ft$estimate), p = ft$p.value)
}
ovt <- rbind(rbindlist(lapply(c("lost", "gained"), function(dr) rbindlist(lapply(names(pairs)[-1], function(b) ov_test(names(pairs)[1], b, dr))))),
             rbindlist(lapply(c("lost", "gained"), function(dr) ov_test(ceil[1], ceil[2], dr))))
fwrite(ovt, file.path(OUT, "overlap_tet2_vs_dnmt3a.tsv"), sep = "\t")
print(ovt[, .(a, b, direction, shared, expected = round(expected, 1), OR = round(odds_ratio, 1), p = signif(p, 2))])

## ---- inactive X: X-linked CpG-island promoters should sit near 0.5 in female cells --
pm <- fread(cmd = sprintf("gzip -dc '%s'", file.path(INT, "region_methylation", "promoter_methylation.tsv.gz")))
mg <- fread(cmd = sprintf("gzip -dc '%s'", file.path(INT, "region_methylation", "mane_select_genes.tsv.gz")))
mg[, `:=`(ps = fifelse(strand == "+", tss - 1500L, tss - 500L), pe = fifelse(strand == "+", tss + 500L, tss + 1500L))]
cg <- copy(cgi); setnames(cg, c("start", "end"), c("s", "e")); setkey(cg, chr, s, e)
cgp <- unique(foverlaps(mg[, .(gene, chr, s = ps, e = pe)], cg, by.x = c("chr", "s", "e"), nomatch = 0L)$gene)
x <- merge(pm, mg[, .(id = gene, chr)], by = "id")[id %in% cgp]
extra <- fread(SAMPLE_SHEET, colClasses = "character")[methylation == "extra", sample]
smp8 <- intersect(c(smp, extra), names(x))
xc <- rbindlist(lapply(smp8, function(s0) data.table(sample = s0,
  chrX_median = median(x[chr == "chrX"][[s0]], na.rm = TRUE),
  chrX_frac_0.3_0.7 = mean(between(x[chr == "chrX"][[s0]], 0.3, 0.7), na.rm = TRUE),
  autosome_median = median(x[chr != "chrX"][[s0]], na.rm = TRUE), n_chrX = sum(!is.na(x[chr == "chrX"][[s0]])))))
fwrite(xc, file.path(OUT, "chrX_cgi_promoters.tsv"), sep = "\t")
print(xc)

message("written to ", OUT)
