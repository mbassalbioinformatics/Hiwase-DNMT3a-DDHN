#!/usr/bin/env Rscript
###############################################################################
# 21_tile_size_check.R                                               2026-10-19
#
# Why the genome tiles are 500 bp and a change is 20 points (config/settings.R).
# Until 19 Oct 2026 the tiles were 1 kb and a tile change 25 points. This script
# measures, from the same per-CpG reads and with the same rules as
# 01_region_methylation.R (CpGs with >= 5 reads; a tile needs >= 3 of them), how
# tiles of 250 bp, 500 bp, 1 kb and 2 kb behave at a 20- and a 25-point cut, in
# every comparison:
#   lost, gained      tiles crossing the cut (raw delta, as in 04_global_direction.R)
#   expected          tiles that sampling noise alone would carry across the cut:
#                     each tile's delta ~ N(median delta, (phi * s)^2), s from the
#                     binomial variance of its CpGs ((m + 0.5) / (t + 1) per CpG), phi
#                     the robust spread (1.4826 x MAD) of the centred z-scores, which
#                     scales the binomial s up to the observed spread (reads cover
#                     several CpGs, and libraries differ by more than counting)
#   real change       lost - expected (and gained - expected), as tiles and as Mb
#   DMRs picked up    metilene DMRs of the comparison overlapped by a tile called in
#                     the same direction: how well the tiles see focal changes
#   CpGs kept         share of the CpGs with >= 5 reads that lie in a tile with a value
# Small tiles see focal changes undiluted but hold fewer CpGs (noisier, and CpG-poor
# DNA drops out below 3 CpGs); large tiles average a focal change away.
#
# Output: integration/global_direction/tile_size_check.tsv, tile_size_info.tsv,
#         tile_size_check.png
###############################################################################
source("config/paths.R")
ID <- with(data.table::fread(SAMPLE_SHEET, colClasses = "character"), setNames(sample, key))   # role key -> sample ID (config/samples.tsv)

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setDTthreads(8)
OUT <- file.path(INT, "global_direction"); dir.create(OUT, showWarnings = FALSE)
SAMPLES <- c(ID[["D18"]], ID[["D30"]], ID[["Y35"]], ID[["DX"]], ID[["D23"]], ID[["T2A"]], ID[["T2DX"]])
KEEP <- paste0("chr", c(1:22, "X"))
MIN_COV <- 5L; MIN_CPG <- 3L
WIDTHS <- c(250L, 500L, 1000L, 2000L)
CUTS   <- c(0.20, 0.25)

read_cpgs <- function(s) {
  f <- file.path(EMS, paste0(s, "_msuite2"), "Msuite2.CpG.meth.call.xz")
  x <- fread(cmd = sprintf("xz -dc '%s'", f), select = c(1, 2, 4, 5, 8, 9),
             col.names = c("chr", "pos", "wC", "wT", "cC", "cT"))
  x <- x[chr %in% KEEP]
  x[, `:=`(m = wT + cT, t = wC + wT + cC + cT)]       # TAPS: T = methylated
  x[t >= MIN_COV, .(chr, pos, m, t)]
}

tiles <- setNames(vector("list", length(WIDTHS)), WIDTHS)
info <- list()
for (s in SAMPLES) {
  x <- read_cpgs(s)
  x[, `:=`(p = m / t, v = ((m + 0.5) / (t + 1)) * (1 - (m + 0.5) / (t + 1)) / t)]
  for (W in WIDTHS) {
    r <- x[, .(n5 = .N, f = sum(p), vs = sum(v)), by = .(chr, tile = (pos - 1L) %/% W)]
    info[[length(info) + 1]] <- data.table(sample = s, width = W, cpgs_5x = nrow(x), median_depth = as.numeric(median(x$t)),
                                           tiles_with_value = r[n5 >= MIN_CPG, .N], cpgs_kept = r[n5 >= MIN_CPG, sum(n5)] / nrow(x),
                                           median_cpgs = as.numeric(r[n5 >= MIN_CPG, median(n5)]),
                                           median_sd_points = 100 * r[n5 >= MIN_CPG, median(sqrt(vs) / n5)])
    r <- r[n5 >= MIN_CPG]
    tiles[[as.character(W)]][[s]] <- r[, .(chr, tile, sample = s, mean = f / n5, sd = sqrt(vs) / n5)]
  }
  message(sprintf("%-10s %.1f M CpGs at >= 5 reads", s, nrow(x) / 1e6))
  rm(x); gc(verbose = FALSE)
}
info <- rbindlist(info)
fwrite(info, file.path(OUT, "tile_size_info.tsv"), sep = "\t")

comp <- fread(COMPARISONS, colClasses = "character")
dmr  <- fread(cmd = sprintf("gzip -dc '%s'", file.path(INT, "dmr", "dmr_master.tsv.gz")))[set == "metilene", .(id, chr, start, end, direction)]
res <- list()
for (W in WIDTHS) {
  L <- rbindlist(tiles[[as.character(W)]])
  wm <- dcast(L, chr + tile ~ sample, value.var = "mean"); ws <- dcast(L, chr + tile ~ sample, value.var = "sd")
  for (k in seq_len(nrow(comp))) {
    a <- comp$g1[k]; b <- comp$g2[k]
    ok <- !is.na(wm[[a]]) & !is.na(wm[[b]])
    d <- wm[[a]][ok] - wm[[b]][ok]; s <- sqrt(ws[[a]][ok]^2 + ws[[b]][ok]^2)
    med <- median(d); phi <- mad((d - med) / s)
    pos <- wm[ok, .(chr, start = tile * W + 1L, end = (tile + 1L) * W)]
    dk <- dmr[id == comp$id[k]]
    for (cut in CUTS) {
      lost <- d <= -cut; gained <- d >= cut
      called <- rbind(pos[lost][, dir := "lost"], pos[gained][, dir := "gained"]); setkey(called, chr, start, end)
      ov <- foverlaps(dk[, .(chr, start, end, direction, row = .I)], called, nomatch = 0L)
      res[[length(res) + 1]] <- data.table(id = comp$id[k], label = comp$label[k], width = W, cut = cut, tiles = sum(ok),
        median_delta_points = 100 * med, phi = phi, lost = sum(lost), gained = sum(gained),
        exp_lost = sum(pnorm((-cut - med) / (phi * s))), exp_gained = sum(pnorm((med - cut) / (phi * s))),
        dmrs = nrow(dk), dmrs_picked_up = uniqueN(ov[dir == direction, row]))
    }
  }
  message("tiles of ", W, " bp done")
  rm(L, wm, ws); gc(verbose = FALSE)
}
res <- rbindlist(res)
res[, `:=`(noise_share = (exp_lost + exp_gained) / pmax(lost + gained, 1),
           real_lost = lost - exp_lost, real_gained = gained - exp_gained,
           real_lost_mb = (lost - exp_lost) * width / 1e6, real_gained_mb = (gained - exp_gained) * width / 1e6,
           dmr_share = dmrs_picked_up / pmax(dmrs, 1))]
fwrite(res, file.path(OUT, "tile_size_check.tsv"), sep = "\t")
print(res[id == paste0(ID[["Y35"]], "_vs_", ID[["D30"]]), .(width, cut, lost, gained, noise_pct = round(100 * noise_share, 1),
                                     real_lost_mb = round(real_lost_mb, 1), dmrs = paste0(dmrs_picked_up, "/", dmrs))])

## ---- figure: the four tile sizes, both cuts, for the main comparisons --------------------------------
KEY <- c(paste0(ID[["Y35"]], "_vs_", ID[["D30"]]), paste0(ID[["DX"]], "_vs_", ID[["D30"]]), paste0(ID[["DX"]], "_vs_", ID[["D23"]]), paste0(ID[["T2DX"]], "_vs_", ID[["T2A"]]), paste0(ID[["D30"]], "_vs_", ID[["D23"]]))
f <- res[id %in% KEY]
f[, comparison := factor(label[1], levels = comp[match(KEY, id), label]), by = id]
f[, cutf := factor(sprintf("%.0f points", 100 * cut), levels = c("20 points", "25 points"))]
fl <- rbind(f[, .(comparison, width, cutf, metric = "A. Methylation lost beyond noise (Mb)", v = real_lost_mb)],
            f[, .(comparison, width, cutf, metric = "B. Calls expected from noise alone (%)", v = 100 * noise_share)],
            f[, .(comparison, width, cutf, metric = "C. metilene DMRs picked up (%)", v = 100 * dmr_share)])
ck <- info[, .(v = 100 * mean(cpgs_kept)), by = width][, `:=`(metric = "D. CpGs kept in tiles with a value (%)", cutf = factor("20 points", levels = levels(f$cutf)),
                                                         comparison = factor("all libraries"))]
p <- ggplot(fl, aes(width, v, colour = comparison, linetype = cutf)) +
  annotate("rect", xmin = 420, xmax = 600, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_line(data = ck, colour = "grey30", linewidth = 0.8) + geom_point(data = ck, colour = "grey30", size = 2) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_x_log10(breaks = WIDTHS, labels = c("250 bp", "500 bp", "1 kb", "2 kb")) +
  scale_colour_manual(values = c("#6A51A3", "#9E9AC8", "#D95F02", "#E6AB02", "grey55"), name = NULL) +
  scale_linetype_manual(values = c(`20 points` = "solid", `25 points` = "22"), name = "a change of at least") +
  labs(x = "tile size (shaded: the size used)", y = NULL, title = "Which tile size and cut: the same reads cut four ways",
       subtitle = paste0("A: tiles lost minus those expected from sampling noise, times the tile size. B: share of the tiles called that noise alone would ",
                         "produce. C: share of the comparison's metilene DMRs overlapped by a tile called in the same direction. D: share of CpGs (>= 5 reads) ",
                         "in tiles with >= 3 of them")) +
  guides(colour = guide_legend(ncol = 3), linetype = guide_legend(ncol = 1)) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 9.5))
ggsave(file.path(OUT, "tile_size_check.png"), p, width = 17, height = 6.2, dpi = 170)
message("written to ", OUT)
