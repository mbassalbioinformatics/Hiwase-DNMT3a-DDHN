#!/usr/bin/env Rscript
###############################################################################
# 20_three_way_course.R                                              2026-10-17
#
# The DNMT3A patient's course as one three-way comparison, as the clinical team
# asked: recipient at 1 month (D30, day 30) -> 3.5 years (Y35) -> DDHN
# diagnosis (DX), all bone-marrow mononuclear cells. Every feature measured in
# all three samples gets a trajectory:
#   a = 3.5 years - 1 month, b = DDHN - 3.5 years, c = DDHN - 1 month (= a + b)
#   with the clone      |a| >= T and |b| < F      changed as the clone took over, then held
#   at diagnosis        |a| < F  and |b| >= T     changed only between 3.5 years and diagnosis
#   progressive         a, b the same sign, both >= F, |c| >= T
#   reversing           a, b opposite signs, both >= T (out and back)
#   stable              |a| < F and |b| < F
#   other               anything else (a partial change in one or both steps)
# T = CUT (config/settings.R: 20 points) for tiles (TILE_BP) and promoters alike, the
# report's threshold for a change; F = 10 points (flat). Direction: of a for "with
# the clone", of b for "at diagnosis", of c for "progressive", "lost then regained" or
# "gained then lost" for "reversing".
# Expression: the RNA-seq arm's classes for the same two steps (02_dge_gene_level.R:
# clone-associated, transformation-associated, progressive, reversing; robust DE with
# the other step flat, |log2FC| < 0.5).
# Check against the donors: a change "with the clone" could also be the 1-month marrow being unusual
# (regenerating, 33 days after the transplant). The two donor samples carry the clone in few cells
# (about 6% and 20%): where the change follows the clone they look like 1 month (within 10 points, and
# >= 15 points from the clone samples); where 1 month is the odd one they look like the clone samples.
# Genes: each tile goes to the nearest MANE TSS within 50 kb (as in 19); a gene's tiles
# in each class are counted; programs per class by over-representation (Fisher) of genes
# with >= 2 tiles of that class and direction, against genes with >= 5 tiles.
#
# Output: integration/clinical/three_way_*.tsv, fig_three_way.png
###############################################################################
source("config/paths.R")
ID <- with(data.table::fread(SAMPLE_SHEET, colClasses = "character"), setNames(sample, key))   # role key -> sample ID (config/samples.tsv)

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(msigdbr) })
OUT  <- file.path(INT, "clinical"); dir.create(OUT, showWarnings = FALSE)
rd   <- function(f, ...) if (grepl("\\.gz$", f)) fread(cmd = sprintf("gzip -dc '%s'", f), ...) else fread(f, ...)
S1 <- ID[["D30"]]; S2 <- ID[["Y35"]]; S3 <- ID[["DX"]]                 # 1 month, 3.5 years, DDHN diagnosis
FLAT <- 0.10; T_TILE <- CUT; T_PROM <- CUT; WINDOW <- 50e3
CLASSES <- c("with the clone", "at diagnosis", "progressive", "reversing", "other", "stable")

trajectory <- function(a, b, T) {
  c_ <- a + b
  cls <- fcase(abs(a) < FLAT & abs(b) < FLAT, "stable",
               abs(a) >= T & abs(b) < FLAT, "with the clone",
               abs(a) < FLAT & abs(b) >= T, "at diagnosis",
               sign(a) == sign(b) & abs(a) >= FLAT & abs(b) >= FLAT & abs(c_) >= T, "progressive",
               sign(a) != sign(b) & abs(a) >= T & abs(b) >= T, "reversing",
               default = "other")
  dir <- fcase(cls == "with the clone", fifelse(a < 0, "lost", "gained"),
               cls == "at diagnosis", fifelse(b < 0, "lost", "gained"),
               cls == "progressive", fifelse(c_ < 0, "lost", "gained"),
               cls == "reversing", fifelse(a < 0, "lost then regained", "gained then lost"),
               default = "")
  list(cls, dir)
}

## ---- genome tiles, with genomic context (as in 04_global_direction.R) ------------------------------
tw <- rd(TILE_FILE)[!is.na(get(S1)) & !is.na(get(S2)) & !is.na(get(S3))]
tw[, c("chr", "start") := tstrsplit(id, ":", fixed = TRUE)][, start := as.integer(start)][, end := start + TILE_BP - 1L]
cgi <- rd(file.path(RES, "cpgIslandExt.txt.gz"), select = 2:4, col.names = c("chr", "start", "end"))[, start := start + 1L]
ov_any <- function(a, b) { setkey(b, chr, start, end)
  hit <- foverlaps(a[, .(id, chr, start, end)], b, nomatch = 0L, which = TRUE); a$id %in% a$id[unique(hit$xid)] }
flank <- function(d, w0, w1) rbind(d[, .(chr, start = pmax(1L, start - w1), end = start - w0 - 1L)],
                                   d[, .(chr, start = end + w0 + 1L, end = end + w1)])[end >= start]
in_cgi <- ov_any(tw, copy(cgi)); in_shore <- ov_any(tw, flank(cgi, 0L, 2000L)); in_shelf <- ov_any(tw, flank(cgi, 2000L, 4000L))
tw[, cgi_context := fifelse(in_cgi, "CpG island", fifelse(in_shore, "shore", fifelse(in_shelf, "shelf", "open sea")))]
cc <- fread(file.path(RES, "encode_cCREs_hg38_ucsc.bed"), select = c(1, 2, 3, 11), col.names = c("chr", "start", "end", "label"))[, start := start + 1L]
pri <- c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3")
tw[, ccre := "none"]
for (lab in rev(pri)) tw[ov_any(tw, cc[label == lab, .(chr, start, end)]), ccre := lab]
tw[, `:=`(a = get(S2) - get(S1), b = get(S3) - get(S2))]
tw[, c("traj", "dir") := trajectory(a, b, T_TILE)]
tw[, region := fcase(ccre %in% c("pELS", "dELS"), "enhancer (pELS, dELS)", ccre == "PLS", "promoter-like (PLS)",
                     cgi_context %in% c("CpG island", "shore"), "other CpG island or shore", default = "other (open sea, shelf)")]
cat("tiles with data in all three samples:", format(nrow(tw), big.mark = ","), "\n")
anchor <- function(m1, clone_mean, donors) fcase(is.na(donors), "donors not measured",
  abs(donors - m1) < FLAT & abs(donors - clone_mean) >= 0.15, "follows the clone",
  abs(donors - clone_mean) < FLAT & abs(donors - m1) >= 0.15, "1 month stands out", default = "in between")
D1 <- ID[["D18"]]; D2 <- ID[["D23"]]                          # donor at donation, donor 5 years after donation
tw[, check := fifelse(traj %in% c("with the clone", "progressive"),
                      anchor(get(S1), (get(S2) + get(S3)) / 2, (get(D1) + get(D2)) / 2), NA_character_)]
ta <- tw[!is.na(check), .N, by = .(traj, dir, check)][order(traj, dir, check)]
print(ta)
fwrite(ta, file.path(OUT, "three_way_donor_check.tsv"), sep = "\t")
tt <- tw[, .N, by = .(traj, dir)][order(match(traj, CLASSES), dir)]
tt[, pct := 100 * N / nrow(tw)]
print(tt)
fwrite(tt, file.path(OUT, "three_way_tiles.tsv"), sep = "\t")
tc <- tw[traj %in% CLASSES[1:4], .N, by = .(traj, dir, region)]
tc <- merge(tc, tw[, .(all = .N), by = region], by = "region")[, pct_of_region := 100 * N / all]
fwrite(tc, file.path(OUT, "three_way_tiles_by_region.tsv"), sep = "\t")
# how many such patterns noise and the change of material alone give: the same rules on the reference
# comparison's samples are not possible (two samples), so the floor is the share of tiles that change by
# >= T_TILE between day 30 and the donor in 2023 (no disease, small clone at both ends)
ref_floor <- tw[!is.na(get(ID[["D23"]])), .(lost = mean(get(ID[["D30"]]) - get(ID[["D23"]]) <= -T_TILE),
                                            gained = mean(get(ID[["D30"]]) - get(ID[["D23"]]) >= T_TILE))]
fwrite(ref_floor, file.path(OUT, "three_way_reference_floor.tsv"), sep = "\t")

## ---- promoters ---------------------------------------------------------------------------------------
pm <- rd(file.path(INT, "region_methylation", "promoter_methylation.tsv.gz"))[!is.na(get(S1)) & !is.na(get(S2)) & !is.na(get(S3))]
pm[, `:=`(gene = id, m1 = get(S1), m2 = get(S2), m3 = get(S3))][, `:=`(a = m2 - m1, b = m3 - m2)]
pm[, c("traj", "dir") := trajectory(a, b, T_PROM)]
pm[, check := fifelse(traj %in% c("with the clone", "progressive"), anchor(m1, (m2 + m3) / 2, (get(D1) + get(D2)) / 2), NA_character_)]
fwrite(pm[!is.na(check), .N, by = .(traj, dir, check)][order(traj, dir, check)], file.path(OUT, "three_way_donor_check_promoters.tsv"), sep = "\t")
pt <- pm[, .N, by = .(traj, dir)][order(match(traj, CLASSES), dir)]
print(pt)
fwrite(pm[traj %in% CLASSES[1:4], .(gene, traj, dir, check, m1 = round(m1, 3), m2 = round(m2, 3), m3 = round(m3, 3),
                                      donor_donation = round(get(D1), 3), donor_5y = round(get(D2), 3))][order(match(traj, CLASSES), dir, -abs(m3 - m1))],
       file.path(OUT, "three_way_promoters.tsv"), sep = "\t")

## ---- genes: tiles per class, promoter class, expression class ------------------------------------------
gn <- rd(file.path(INT, "region_methylation", "mane_select_genes.tsv.gz"))[!duplicated(gene), .(gene, chr, tss)]
setkey(gn, chr, tss)
near <- gn[tw[, .(id, chr, tss = start + TILE_BP %/% 2L - 1L)], on = .(chr, tss), roll = "nearest", .(id, gene, dist = abs(x.tss - i.tss))]
tw <- merge(tw, near[!is.na(gene) & dist <= WINDOW, .(id, gene)], by = "id", all.x = TRUE)
# every tile in a trajectory class, with its gene: the locus plots of 22_part_a_locus_plots.R centre on these
fwrite(tw[traj %in% CLASSES[1:4], .(id, gene, traj, dir, check, a = round(a, 4), b = round(b, 4))],
       file.path(OUT, "three_way_tile_classes.tsv.gz"), sep = "\t")
gt <- tw[!is.na(gene), .(tiles = .N,
                         clone_lost = sum(traj == "with the clone" & dir == "lost"), clone_gained = sum(traj == "with the clone" & dir == "gained"),
                         dx_lost = sum(traj == "at diagnosis" & dir == "lost"), dx_gained = sum(traj == "at diagnosis" & dir == "gained"),
                         prog_lost = sum(traj == "progressive" & dir == "lost"), prog_gained = sum(traj == "progressive" & dir == "gained"),
                         rev = sum(traj == "reversing"),
                         follow_lost = sum(traj == "with the clone" & dir == "lost" & check %in% "follows the clone"),
                         follow_gained = sum(traj == "with the clone" & dir == "gained" & check %in% "follows the clone")), by = gene]
ax <- fread(file.path(RNA, "05_processed_counts_genes", "13_gene_axis_classification.tsv"))
e  <- new.env(); load(file.path(RNA, "05_processed_counts_genes", "dge_gene_level.RData"), envir = e)
cpm <- e$cpm_mat
gt <- merge(gt, pm[, .(gene, prom_traj = traj, prom_dir = dir, prom_m1 = m1, prom_m2 = m2, prom_m3 = m3)], by = "gene", all.x = TRUE)
gt <- merge(gt, ax[, .(gene, rna_class = class, lfc_clone, lfc_transform, ruler_suspect)], by = "gene", all.x = TRUE)
for (s in c(S1, S2, S3)) if (s %in% colnames(cpm)) gt[, (paste0("cpm_", s)) := cpm[match(gene, rownames(cpm)), s]]
fwrite(gt, file.path(OUT, "three_way_genes.tsv.gz"), sep = "\t")
rc <- ax[, .N, by = class][order(-N)]
print(rc)

## ---- programs per trajectory class (genes with >= 2 tiles of the class) ----------------------------------
H  <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))[, .(set = gs_name, gene = gene_symbol)]
C8 <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C8"))[grepl("^HAY_BONE_MARROW", gs_name), .(set = gs_name, gene = gene_symbol)]
PN <- fread(file.path(RES, "gene_panels.tsv"))[, .(set = paste0("PANEL_", panel), gene)]
SETS <- unique(rbind(H, C8, PN))
bgg <- gt[tiles >= 5, gene]
ora <- rbindlist(lapply(c("clone_lost", "clone_gained", "follow_lost", "follow_gained", "dx_lost", "dx_gained", "prog_lost", "prog_gained"), function(cl_) {
  hit <- gt[tiles >= 5 & get(cl_) >= 2, gene]
  if (length(hit) < 10) return(NULL)
  N <- length(bgg); K <- length(hit)
  r <- SETS[gene %in% bgg, .(n = .N, in_class = sum(gene %in% hit), genes = paste(head(intersect(gene, hit), 12), collapse = ", ")), by = set][n >= 10]
  r[, `:=`(class = cl_, hits = K, expected = n * K / N, p = phyper(in_class - 1, n, N - n, K, lower.tail = FALSE))]
  r[, `:=`(odds_ratio = (in_class / (K - in_class)) / ((n - in_class) / (N - n - K + in_class)), fdr = p.adjust(p, "BH"))]
  r
}))
fwrite(ora[order(class, p)], file.path(OUT, "three_way_programs.tsv"), sep = "\t")
print(ora[fdr < 0.05][order(class, p), .(class, set, in_class, n, OR = round(odds_ratio, 2), fdr = signif(fdr, 2))])

## ---- figure ---------------------------------------------------------------------------------------------------
LAB3 <- c("1 month", "3.5 years", "DDHN\ndiagnosis")
# named genes per class, for the report
nm_cls <- rbind(
  gt[clone_lost >= 2][order(-clone_lost, -tiles)][, .(class = "with the clone, lost", gene, tiles = clone_lost)],
  gt[clone_gained >= 2][order(-clone_gained, -tiles)][, .(class = "with the clone, gained", gene, tiles = clone_gained)],
  gt[dx_lost >= 1][order(-dx_lost, -tiles)][, .(class = "at diagnosis, lost", gene, tiles = dx_lost)],
  gt[dx_gained >= 1][order(-dx_gained, -tiles)][, .(class = "at diagnosis, gained", gene, tiles = dx_gained)],
  gt[prog_lost >= 2][order(-prog_lost, -tiles)][, .(class = "progressive, lost", gene, tiles = prog_lost)],
  gt[prog_gained >= 2][order(-prog_gained, -tiles)][, .(class = "progressive, gained", gene, tiles = prog_gained)],
  gt[follow_lost >= 2][order(-follow_lost, -tiles)][, .(class = "with the clone and following it, lost", gene, tiles = follow_lost)],
  gt[follow_gained >= 2][order(-follow_gained, -tiles)][, .(class = "with the clone and following it, gained", gene, tiles = follow_gained)])
fwrite(nm_cls, file.path(OUT, "three_way_named_genes.tsv"), sep = "\t")
fwrite(ax[class %in% c("clone-associated", "transformation-associated", "progressive", "reversing"),
          .(gene, class, lfc_clone = round(lfc_clone, 2), lfc_transform = round(lfc_transform, 2), ruler_suspect)][order(class, -abs(lfc_clone + lfc_transform))],
       file.path(OUT, "three_way_rna_genes.tsv"), sep = "\t")
# A: tiles per trajectory class and direction, split by region
fa <- tc[, .(N = sum(N)), by = .(traj, dir, region)]
fa[, traj := factor(traj, levels = CLASSES[1:4])]
fa[, dir2 := factor(dir, levels = c("lost", "gained", "lost then regained", "gained then lost"))]
fa[, region := factor(region, levels = c("enhancer (pELS, dELS)", "promoter-like (PLS)", "other CpG island or shore", "other (open sea, shelf)"))]
pA <- ggplot(fa, aes(dir2, N, fill = region)) + geom_col(width = 0.75) +
  geom_text(data = fa[, .(N = sum(N)), by = .(traj, dir2)], aes(dir2, N, label = format(N, big.mark = ",")), inherit.aes = FALSE, vjust = -0.3, size = 3.4) +
  facet_grid(~ traj, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("#D95F02", "#7570B3", "#1B9E77", "#BDBDBD"), name = NULL) +
  scale_y_continuous(labels = function(v) format(v, big.mark = ","), expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = paste(TILE_TXT, "tiles"), title = sprintf("A. How each %s tile moved across the three samples", TILE_TXT),
       subtitle = sprintf("tiles with data in all three: %s; change >= %s points, flat < 10 points", format(nrow(tw), big.mark = ","), CUT_TXT)) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.minor = element_blank(),
                                   plot.subtitle = element_text(size = 10))
# B2: the clone step checked against the donors
fb <- ta[traj == "with the clone"][, check := factor(check, levels = c("follows the clone", "in between", "1 month stands out", "donors not measured"))]
fb[, dir := factor(dir, levels = c("lost", "gained"))]
pD <- ggplot(fb, aes(dir, N, fill = check)) + geom_col(width = 0.7) +
  scale_fill_manual(values = c("follows the clone" = "#6A51A3", "in between" = "#BCBDDC", "1 month stands out" = "#FDAE6B", "donors not measured" = "#E0E0E0"), name = NULL) +
  scale_y_continuous(labels = function(v) format(v, big.mark = ","), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = paste(TILE_TXT, "tiles"), title = "B. The clone step against the donors", subtitle = "tiles changed with the clone: do the donors look like 1 month?") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 10))
# B: example loci: promoter methylation and expression at the three samples
# methylation lost with the clone and expression up (HOXA3, HOXA5, HOXB3); expression falling only at diagnosis with flat
# methylation (MEIS1); expression down with the clone and back up at diagnosis, methylation flat (ITGA2B)
EX <- c("HOXA3", "HOXA5", "HOXB3", "MEIS1", "ITGA2B")
ex <- rbind(
  pm[gene %in% EX, .(gene, layer = "promoter methylation (%)", s = LAB3, v = 100 * c(m1, m2, m3)), by = gene][, gene := NULL],
  rbindlist(lapply(EX, function(g) if (g %in% rownames(cpm)) data.table(gene = g, layer = "expression (CPM)", s = LAB3, v = as.numeric(cpm[g, c(S1, S2, S3)])))))
ex[, s := factor(s, levels = LAB3)][, gene := factor(gene, levels = EX)]
ex[, layer := factor(layer, levels = c("promoter methylation (%)", "expression (CPM)"))]
pB <- ggplot(ex, aes(s, v, group = gene, colour = gene)) + geom_line(linewidth = 1) + geom_point(size = 2.6) +
  facet_wrap(~ layer, scales = "free_y") +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(x = NULL, y = NULL, title = "C. Example genes across the three samples", subtitle = "left: promoter methylation; right: expression") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 10))
# C: expression classes
rc2 <- copy(rc)[class != "unchanged"]
rc2[, class := factor(class, levels = c("clone-associated", "transformation-associated", "progressive", "reversing", "ambiguous"))]
pC <- ggplot(rc2, aes(class, N)) + geom_col(fill = "#4C72B0", width = 0.7) + geom_text(aes(label = N), vjust = -0.3, size = 3.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "genes", title = "D. Genes by expression trajectory", subtitle = "RNA-seq: robust change in one step, the other flat") +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 10))
ggsave(file.path(OUT, "fig_three_way.png"), pA / (pD | pB | pC) + plot_layout(heights = c(1, 1.05)), width = 16, height = 12.5, dpi = 170)
message("written to ", OUT)
