#!/usr/bin/env Rscript
###############################################################################
# 14_allelic_imbalance.R                                             2026-10-08
#
# Loss of heterozygosity (LOH) along the genome, from the same TAPS libraries.
# At CpG-destroying transversion SNPs the reads show a base that methylation can
# never produce, so other / (C + T + other) on each strand is the B-allele
# fraction (BAF) of the SNP (counts from 02_methylation_local/allele_balance_counts.py
# at the sites of the sample-identity check).
#
# Heterozygous sites are defined per person from the pooled reads of all that
# person's samples (pooled BAF 0.2-0.8, >= 12 reads), so a site stays in even
# when one sample has lost an allele; a person with a single library cannot be
# scanned this way and is left out. In each sample, a heterozygous site
# with >= 6 reads on the strand gives a BAF; in 10-Mb bins, the mean |BAF - 0.5|
# is compared with its expectation for a balanced heterozygote at the same depths
# (binomial, exact). Excess > 0.10 with z > 4 = allelic imbalance (LOH or a
# copy-number change in a large share of cells).
#
# Output in integration/vaf/: allelic_imbalance_bins.tsv, allelic_imbalance_segments.tsv,
# allelic_imbalance_genome.tsv, allelic_imbalance_flanks.tsv, baf_genome.png, baf_chr4_chr18_tet2_case.png
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source("config/paths.R")
IDD <- file.path(EMS, "00_sample_identity")
OUT <- file.path(INT, "vaf"); dir.create(OUT, showWarnings = FALSE)
MIN_STRAND <- 6L; MIN_POOLED <- 12L; BIN <- 10e6; EXCESS <- 0.10; ZMIN <- 4

.sheet <- fread(SAMPLE_SHEET, colClasses = "character")
SMP <- .sheet[methylation %in% c("analysed", "extra"), .(sample, person, label = plot_label)]
CHR <- paste0("chr", c(1:22, "X"))
# the TET2 case's two libraries: the flank check and the chr4/chr18 figure
TET2_CASE <- .sheet[methylation == "analysed" & case == unique(.sheet[methylation == "analysed", case])[2], sample]   # the second case

d <- rbindlist(lapply(SMP$sample, function(s) {
  x <- fread(file.path(IDD, paste0(s, ".strand_counts.tsv")))
  rbind(x[, .(sample = s, chr, pos, n = w_ct + w_other, alt = w_other)],
        x[, .(sample = s, chr, pos = pos + 1L, n = c_ct + c_other, alt = c_other)])
}))
d <- merge(d[chr %in% CHR], SMP[, .(sample, person)], by = "sample")
# leave out the short arms of the acrocentric chromosomes (repeats, poorly assembled; hg38 centromere ends)
acro <- data.table(chr = c("chr13", "chr14", "chr15", "chr21", "chr22"), p_end = c(18.1e6, 18.2e6, 20.5e6, 13.0e6, 15.1e6))
d <- d[!acro, on = .(chr, pos < p_end)]
d <- d[person %in% SMP[, .N, by = person][N >= 2, person]]
SMP <- SMP[sample %in% d$sample]

# heterozygous in the person (pooled over the person's samples)
het <- d[, .(pn = sum(n), pa = sum(alt)), by = .(person, chr, pos)][pn >= MIN_POOLED & pa / pn >= 0.2 & pa / pn <= 0.8]
b <- merge(d[n >= MIN_STRAND], het[, .(person, chr, pos)], by = c("person", "chr", "pos"))
b[, baf := alt / n]
message("heterozygous sites per sample (>= ", MIN_STRAND, " reads):")
print(b[, .N, by = sample])

# expected |BAF - 0.5| and its variance for a balanced heterozygote at depth n
nmax <- max(b$n)
ea <- vapply(seq_len(nmax), function(n) sum(dbinom(0:n, n, 0.5) * abs(0:n / n - 0.5)), 0)
va <- 0.25 / seq_len(nmax) - ea^2
b[, `:=`(dev = abs(baf - 0.5), e_dev = ea[n], v_dev = va[n])]

b[, bin := floor(pos / BIN)]
bins <- b[, .(sites = .N, obs = mean(dev), exp = mean(e_dev), se = sqrt(sum(v_dev)) / .N), by = .(sample, person, chr, bin)]
bins[, `:=`(excess = obs - exp, z = (obs - exp) / se, start = bin * BIN, end = (bin + 1) * BIN)]
bins[, imbalanced := sites >= 15 & excess > EXCESS & z > ZMIN]
bins <- merge(bins, SMP[, .(sample, label)], by = "sample")
bins[, chr := factor(chr, levels = CHR)]
setorder(bins, sample, chr, bin)
fwrite(bins[, .(sample, label, person, chr, start, end, sites, mean_abs_dev = round(obs, 3), expected = round(exp, 3),
                excess = round(excess, 3), z = round(z, 1), imbalanced)], file.path(OUT, "allelic_imbalance_bins.tsv"), sep = "\t")
gw <- b[, .(het_sites = .N, median_reads = as.double(median(n)), mean_abs_dev = round(mean(dev), 4),
            expected = round(mean(e_dev), 4), excess = round(mean(dev) - mean(e_dev), 4)), by = sample]
gw <- merge(SMP[, .(sample, label)], gw, by = "sample")[match(SMP$sample, sample)]
gw <- merge(gw, bins[, .(bins_tested = sum(sites >= 15), bins_imbalanced = sum(imbalanced)), by = sample], by = "sample")[match(SMP$sample, sample)]
fwrite(gw, file.path(OUT, "allelic_imbalance_genome.tsv"), sep = "\t")
cat("\ngenome-wide:\n"); print(gw)

# merge adjacent imbalanced bins into segments
seg <- bins[imbalanced == TRUE][, run := cumsum(c(1, diff(bin) != 1 | diff(as.integer(chr)) != 0)), by = sample]
seg <- seg[, .(start = min(start), end = max(end), bins = .N, sites = sum(sites),
               excess = round(weighted.mean(excess, sites), 3)), by = .(sample, label, chr, run)][, run := NULL]
fwrite(seg, file.path(OUT, "allelic_imbalance_segments.tsv"), sep = "\t")
cat("\nimbalanced segments (10-Mb bins merged):\n"); print(seg)

# the two common SNPs whose allele balance differs between the TET2-case samples (13_vaf_from_taps_reads.R):
# the nearest heterozygous CpG-SNPs on each side, to see how far any imbalance extends
FLANK <- data.table(snp = c("TET2 rs2454206", "SETBP1 rs3744825"), chr = c("chr4", "chr18"), pos = c(105275794L, 44952641L))
fk <- rbindlist(lapply(seq_len(nrow(FLANK)), function(i) {
  h <- het[person == SMP[sample == TET2_CASE[1], person] & chr == FLANK$chr[i]]
  up <- h[pos < FLANK$pos[i]][order(-pos)][1:3]; dn <- h[pos > FLANK$pos[i]][order(pos)][1:3]
  sites <- rbind(up[, side := "before"], dn[, side := "after"])[!is.na(pos)]
  w <- dcast(d[sample %in% TET2_CASE & chr == FLANK$chr[i] & pos %in% sites$pos], pos ~ sample,
             value.var = c("alt", "n"))
  merge(sites[, .(snp = FLANK$snp[i], chr, pos, side, distance_kb = round(abs(pos - FLANK$pos[i]) / 1e3, 1))], w, by = "pos")
}))
for (s_ in TET2_CASE) fk[, (paste0("baf_", s_)) := sprintf("%d/%d", get(paste0("alt_", s_)), get(paste0("n_", s_)))]
fk <- fk[order(snp, side == "after", distance_kb), c("snp", "chr", "pos", "side", "distance_kb", paste0("baf_", TET2_CASE)), with = FALSE]
fwrite(fk, file.path(OUT, "allelic_imbalance_flanks.tsv"), sep = "\t")
cat("\nnearest heterozygous CpG-SNPs around the two imbalanced SNPs (TET2 case):\n"); print(fk)

## ---- figures -----------------------------------------------------------------
# genome-wide BAF: chromosomes laid end to end
len <- b[, .(len = max(pos)), by = chr][match(CHR, chr)]
len[, off := c(0, cumsum(as.numeric(len))[-.N])]
b <- merge(b, len[, .(chr, off)], by = "chr"); b[, gx := off + pos]
bins <- merge(bins, len[, .(chr, off)], by = "chr")
b[, lab := factor(SMP$label[match(sample, SMP$sample)], levels = SMP$label)]
bins[, lab := factor(label, levels = SMP$label)]
mid <- len[, .(chr = sub("chr", "", chr), x = off + len / 2)]
shade <- len[seq(2, .N, 2), .(xmin = off, xmax = off + len)]
p <- ggplot(b, aes(gx, baf)) +
  geom_rect(data = shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey93") +
  geom_point(size = 0.15, alpha = 0.25, colour = "grey20") +
  geom_segment(data = bins[imbalanced == TRUE], aes(x = off + start, xend = off + end, y = 1.06, yend = 1.06),
               inherit.aes = FALSE, colour = "#B2182B", linewidth = 2.2) +
  facet_grid(lab ~ ., switch = "y") +
  scale_x_continuous(breaks = mid$x, labels = mid$chr, expand = c(0.005, 0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.1)) +
  labs(x = NULL, y = "B-allele fraction at heterozygous CpG-SNPs",
       title = "Allelic balance along the genome, from the TAPS reads",
       subtitle = paste0("each dot = one heterozygous SNP in one sample (>= ", MIN_STRAND,
                         " reads); red bar = 10-Mb bin with allelic imbalance (excess |BAF - 0.5| > ", EXCESS, ", z > ", ZMIN, ")")) +
  theme_bw(base_size = 12) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 10), strip.placement = "outside",
        panel.grid = element_blank(), axis.text.x = element_text(size = 8.5), panel.spacing.y = unit(2, "mm"))
ggsave(file.path(OUT, "baf_genome.png"), p, width = 16, height = 1.2 + 1.15 * nrow(SMP), dpi = 200)

z <- b[sample %in% TET2_CASE & chr %in% c("chr4", "chr18")][, chr := factor(chr, levels = c("chr4", "chr18"))]
zb <- bins[sample %in% TET2_CASE & chr %in% c("chr4", "chr18")][, chr := factor(as.character(chr), levels = c("chr4", "chr18"))]
mk <- data.table(chr = factor(c("chr4", "chr18"), levels = c("chr4", "chr18")), pos = c(105146e3, 44680e3), gene = c("TET2", "SETBP1"))
p2 <- ggplot(z, aes(pos / 1e6, baf)) +
  geom_vline(data = mk, aes(xintercept = pos / 1e6), colour = "#D95F02", linetype = 2) +
  geom_text(data = mk, aes(x = pos / 1e6, y = 1.08, label = gene), inherit.aes = FALSE, colour = "#D95F02",
            hjust = -0.1, size = 3.6, fontface = "italic") +
  geom_point(size = 0.6, alpha = 0.45) +
  geom_step(data = zb, aes(x = start / 1e6, y = 0.5 + obs), colour = "#B2182B", linewidth = 0.6, inherit.aes = FALSE) +
  geom_step(data = zb, aes(x = start / 1e6, y = 0.5 - obs), colour = "#B2182B", linewidth = 0.6, inherit.aes = FALSE) +
  facet_grid(lab ~ chr, scales = "free_x", space = "free_x") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1.12)) +
  labs(x = "position (Mb, hg38)", y = "B-allele fraction",
       title = "TET2 case: allelic balance on chromosomes 4 and 18",
       subtitle = "red lines = 0.5 +/- mean |BAF - 0.5| per 10-Mb bin; dashed = TET2 and SETBP1") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave(file.path(OUT, "baf_chr4_chr18_tet2_case.png"), p2, width = 13, height = 6.5, dpi = 200)
message("written to ", OUT)
