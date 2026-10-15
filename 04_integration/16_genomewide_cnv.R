#!/usr/bin/env Rscript
###############################################################################
# 16_genomewide_cnv.R                                                2026-10-15
#
# Copy number and allelic balance across the whole genome, from the GATK
# copy-number run on the TAPS reads (05_variants_server/04_cnv_loh.sh):
#   copy ratio      reads per 1-kb bin (mapping quality >= 10); bins kept when
#                   mappable (Umap k100 >= 0.9), outside segmental duplications
#                   (<= 0.5) and of ordinary GC (0.1-0.9); each library denoised
#                   against the panel of normal libraries (config/samples.tsv,
#                   column cnv_panel; a normal against the others)
#   allele balance  common SNPs (gnomAD AF 0.1-0.9; SNPs whose alternative allele
#                   creates a CpG left out); heterozygous sites from the person's
#                   first normal library (config/samples.tsv, column variant_role)
#   segments        GATK ModelSegments; gains and losses from CallCopyRatioSegments
# Here: a summary per library; every called gain or loss, classified by whether
# the person's own normal library shows it too (inherited) or not (acquired);
# the duplicated block inside KMT2A; and TET2 / SETBP1 in the TET2 case, where
# the 85-gene panel saw allelic imbalance at one SNP each (14_allelic_imbalance.R).
#
# Input   emseq/04_variants/cnv/ as copied from the server (environment variable CNV_DIR
#         overrides). The full allelic counts (allelic/) are optional: when
#         present, the reads around TET2 and SETBP1 are extracted once into
#         integration/cnv/tet2_case_allelic_counts.tsv.gz.
# Output  integration/cnv/: cnv_library_summary.tsv, cnv_called_segments.tsv,
#         cnv_low_maf_segments.tsv, kmt2a_block.tsv, kmt2a_bins.tsv, tet2_case_genes.tsv, cnv_run.tsv, cnv_inputs.tsv,
#         cnv_genome.png, cnv_kmt2a.png, cnv_tet2_case_genes.png
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("config/paths.R")                          # EMS, INT, MANE, SAMPLE_SHEET
CNV  <- Sys.getenv("CNV_DIR", file.path(EMS, "04_variants", "cnv"))
OUT  <- file.path(INT, "cnv"); dir.create(OUT, showWarnings = FALSE)
stopifnot(dir.exists(file.path(CNV, "segments")), dir.exists(file.path(CNV, "denoised")))

# libraries from the sample sheet; normal = in the copy-number panel of normals (column cnv_panel)
SMP <- fread(SAMPLE_SHEET, colClasses = "character")[methylation == "analysed",
         .(sample, person, case = paste(case, "case"), normal = cnv_panel == "yes", label = short_label)]
CASES <- unique(SMP$case)                        # e.g. "DNMT3A case", "TET2 case"
CASE_COL <- setNames(c("#6A51A3", "#D95F02")[seq_along(CASES)], CASES)
CHR  <- paste0("chr", c(1:22, "X"))
LV   <- SMP$label
TR   <- c("Copy ratio", "Allele fraction")
# short arms of the acrocentric chromosomes: repeats, few reliable bins (hg38 p-arm ends, as in 14_allelic_imbalance.R)
ACRO <- data.table(chr = c("chr13", "chr14", "chr15", "chr21", "chr22"), p_end = c(18.1e6, 18.2e6, 20.5e6, 13.0e6, 15.1e6))
in_acro <- function(ch, a) { p <- ACRO$p_end[match(ch, ACRO$chr)]; !is.na(p) & a < p }     # starts in the short arm or centromere
LOH_MAF <- 0.40      # a segment with minor-allele fraction below this would count as allelic imbalance
FLANK   <- 60e3      # KMT2A: bins this far either side of the block form the local baseline
ACQ     <- "clone library only: acquired"

gat <- function(f) {                      # GATK tables: SAM-style @ header lines, then a column header
  if (!file.exists(f) && file.exists(paste0(f, ".gz"))) f <- paste0(f, ".gz")
  if (!file.exists(f)) stop("missing input: ", f)
  inputs <<- c(inputs, f)
  fread(cmd = sprintf("%s '%s' | grep -v '^@'", if (grepl("\\.gz$", f)) "gzip -dc" else "cat", f))
}
inputs <- character(0)

## ---- read ------------------------------------------------------------------------
cr <- rbindlist(lapply(SMP$sample, function(s) gat(file.path(CNV, "denoised", paste0(s, ".denoisedCR.tsv")))[, sample := s]))
setnames(cr, c("CONTIG", "START", "END", "LOG2_COPY_RATIO"), c("chr", "start", "end", "l2"))
cr <- cr[chr %in% CHR]
cr[, `:=`(chr = factor(chr, levels = CHR), sample = factor(sample, levels = SMP$sample))]
setkey(cr, sample, chr, start)
mf  <- rbindlist(lapply(SMP$sample, function(s) gat(file.path(CNV, "segments", paste0(s, ".modelFinal.seg")))[, sample := s]))
cs  <- rbindlist(lapply(SMP$sample, function(s) gat(file.path(CNV, "segments", paste0(s, ".called.seg")))[, sample := s]))
het <- rbindlist(lapply(SMP$sample, function(s) gat(file.path(CNV, "segments", paste0(s, ".hets.tsv")))[, sample := s]))
setnames(mf, c("CONTIG", "START", "END"), c("chr", "start", "end"))
setnames(cs, c("CONTIG", "START", "END", "NUM_POINTS_COPY_RATIO", "MEAN_LOG2_COPY_RATIO", "CALL"), c("chr", "start", "end", "bins", "l2", "call"))
setnames(het, c("CONTIG", "POSITION", "REF_COUNT", "ALT_COUNT"), c("chr", "pos", "ref_n", "alt_n"))
mf <- mf[chr %in% CHR]; cs <- cs[chr %in% CHR]; het <- het[chr %in% CHR]
het[, `:=`(n = ref_n + alt_n, baf = alt_n / (ref_n + alt_n))]

# MANE Select exons and gene spans
gtf <- fread(cmd = sprintf("gzip -dc '%s'", MANE), sep = "\t", header = FALSE, quote = "",
             col.names = c("chr", "src", "type", "start", "end", "score", "strand", "frame", "attr"))
gtf <- gtf[type == "exon" & grepl('tag "MANE_Select"', attr, fixed = TRUE)]
gtf[, `:=`(gene = sub('.*gene_name "([^"]+)".*', "\\1", attr), tx = sub('.*transcript_id "([^"]+)".*', "\\1", attr),
           exon = as.integer(sub(".*exon_number ([0-9]+);.*", "\\1", attr)))]
gtf[, attr := NULL]
gene_spans <- gtf[, .(start = min(start), end = max(end), tx = tx[1]), by = .(gene, chr)]
inputs <- c(inputs, MANE)

## ---- per library -----------------------------------------------------------------
noise <- cr[, .(d = diff(l2)), by = .(sample, chr)][, .(bin_noise = round(mad(d) / sqrt(2), 3)), by = sample]
mf[, `:=`(mb = (end - start + 1) / 1e6, acro = in_acro(chr, start))]
low <- mf[NUM_POINTS_ALLELE_FRACTION >= 20 & !acro & MINOR_ALLELE_FRACTION_POSTERIOR_50 < LOH_MAF]
low <- merge(low, SMP[, .(sample, person, normal)], by = "sample")
# a low stretch in a clone library is shared when a normal library of the same person is also low there
low[, shared := mapply(function(s, p, ch, a, b) {
  nrm <- SMP[person == p & normal == TRUE & sample != s, sample]
  nrow(mf[sample %in% nrm & chr == ch & start <= b & end >= a & NUM_POINTS_ALLELE_FRACTION >= 20 &
          MINOR_ALLELE_FRACTION_POSTERIOR_50 < LOH_MAF + 0.02]) > 0 }, sample, person, chr, start, end)]
fwrite(low[, .(sample, chr, start, end, mb = round(mb, 3), het_points = NUM_POINTS_ALLELE_FRACTION,
               maf = MINOR_ALLELE_FRACTION_POSTERIOR_50, copy_ratio = round(2^LOG2_COPY_RATIO_POSTERIOR_50, 2),
               library = ifelse(normal, "normal", "clone"), shared_with_normal = ifelse(normal, NA, shared))],
       file.path(OUT, "cnv_low_maf_segments.tsv"), sep = "\t")
summ <- SMP[, .(sample, label, case)]
summ <- merge(summ, cr[, .(bins = .N), by = sample], by = "sample")
summ <- merge(summ, noise, by = "sample")
# expected |BAF - 0.5| of a balanced heterozygote at depth n (exact binomial)
ea <- vapply(1:500, function(n) sum(dbinom(0:n, n, 0.5) * abs(0:n / n - 0.5)), 0)
summ <- merge(summ, het[, .(het_sites = .N, het_median_reads = as.double(median(n))), by = sample], by = "sample")
summ <- merge(summ, het[n >= 6 & n <= 500, .(baf_dev = round(mean(abs(baf - 0.5)), 3), baf_dev_expected = round(mean(ea[n]), 3)), by = sample], by = "sample")
# how often a heterozygous SNP reads only one allele by chance (>= 12 reads): the base rate for a single imbalanced SNP
summ <- merge(summ, het[n >= 12, .(het_sites_ge12 = .N, one_allele_ge12 = sum(ref_n == 0 | alt_n == 0)), by = sample], by = "sample")
summ <- merge(summ, mf[NUM_POINTS_ALLELE_FRACTION >= 20, .(maf_typical = round(weighted.mean(MINOR_ALLELE_FRACTION_POSTERIOR_50, mb), 3)), by = sample], by = "sample")
summ <- merge(summ, mf[, .(segments = .N), by = sample], by = "sample")
summ <- merge(summ, low[, .(low_maf_segments = .N, low_maf_mb = round(sum(mb), 2), low_maf_largest_mb = round(max(mb), 2),
                            low_maf_lowest = round(min(MINOR_ALLELE_FRACTION_POSTERIOR_50), 3),
                            low_maf_acquired = sum(!normal & !shared)), by = sample], by = "sample", all.x = TRUE)
for (k in c("low_maf_segments", "low_maf_mb", "low_maf_largest_mb", "low_maf_acquired")) summ[is.na(get(k)), (k) := 0]

## ---- called gains and losses: inherited or acquired? -------------------------------
cs <- cs[call != "0"]
# compared in copy-ratio space: a bin without reads has log2 ~ -30 and would dominate a mean of log2 values
seg_cr <- function(s, ch, a, b) cr[.(s, ch)][start >= a & end <= b, mean(2^l2)]
for (t in SMP$sample) cs[, (paste0("cr_", t)) := mapply(function(ch, a, b) seg_cr(t, ch, a, b), chr, start, end)]
cs <- merge(cs, SMP[, .(sample, person, normal, label)], by = "sample")
cls <- function(i) {
  r <- cs[i]
  own <- r[[paste0("cr_", r$sample)]] - 1
  same_dir <- function(v) !is.na(v) && sign(v - 1) == sign(own) && abs(v - 1) >= 0.5 * abs(own)
  if (r$normal) return("normal library: inherited or technical")
  nrm <- SMP[person == r$person & normal == TRUE & sample != r$sample, sample]
  if (any(vapply(nrm, function(t) same_dir(r[[paste0("cr_", t)]]), TRUE))) return("also in the person's normal library: inherited")
  if (in_acro(r$chr, r$start)) return("short arm of an acrocentric chromosome: unreliable")
  ACQ
}
cs[, class := vapply(seq_len(.N), cls, "")]
cs[, genes := mapply(function(ch, a, b) { g <- gene_spans[chr == ch & start <= b & end >= a, gene]
                                          if (!length(g)) "" else paste(c(head(g, 4), if (length(g) > 4) "..."), collapse = ", ") },
                     chr, start, end)]
cs[, `:=`(size_kb = round((end - start + 1) / 1e3, 1), o1 = class != ACQ, o2 = match(chr, CHR))]
cs[, copy_ratio := round(vapply(seq_len(.N), function(i) cs[[paste0("cr_", sample[i])]][i], 0), 2)]
setorder(cs, o1, sample, o2, start)
lcols <- paste0("cr_", SMP$sample)
csw <- cs[, c("sample", "label", "chr", "start", "end", "size_kb", "bins", "l2", "copy_ratio", "call", "class", "genes", lcols), with = FALSE]
csw[, (lcols) := lapply(.SD, round, 2), .SDcols = lcols]
fwrite(csw, file.path(OUT, "cnv_called_segments.tsv"), sep = "\t")
summ <- merge(summ, cs[, .(called = .N, called_acquired = sum(class == ACQ)), by = sample], by = "sample", all.x = TRUE)
summ[is.na(called), `:=`(called = 0L, called_acquired = 0L)]
summ <- summ[match(SMP$sample, sample)]
fwrite(summ, file.path(OUT, "cnv_library_summary.tsv"), sep = "\t")
cat("per library:\n"); print(summ)
cat("\ncalled gains and losses by class:\n"); print(cs[, .N, by = class])
print(cs[class == ACQ, .(sample, chr, start, end, size_kb, copy_ratio, genes)])

## ---- the KMT2A block ---------------------------------------------------------------
km <- gene_spans[gene == "KMT2A"]
kb <- cs[class == ACQ & chr == km$chr & start <= km$end & end >= km$start & call == "+"]
stopifnot(nrow(kb) >= 1)
B0 <- max(kb$start); B1 <- min(kb$end)                      # the stretch every clone library calls
kex <- gtf[gene == "KMT2A"][order(exon)]
dup <- kex[start >= B0 & end <= B1, exon]
blk <- rbindlist(lapply(SMP$sample, function(s) {
  x <- cr[.(s, km$chr)]
  i <- x[start >= B0 & end <= B1, l2]
  f <- x[(end < B0 & start >= B0 - FLANK) | (start > B1 & end <= B1 + FLANK), l2]
  data.table(sample = s, bins = length(i), l2 = mean(i), se = sd(i) / sqrt(length(i)), flank_bins = length(f), flank_l2 = mean(f))
}))
blk[, `:=`(copy_ratio = 2^l2, ci_low = 2^(l2 - 1.96 * se), ci_high = 2^(l2 + 1.96 * se), flank_copy_ratio = 2^flank_l2)]
blk[, `:=`(chr = km$chr, block_start = B0, block_end = B1, block_kb = round((B1 - B0 + 1) / 1e3, 1), transcript = km$tx,
           first_exon = min(dup), last_exon = max(dup), exons_total = nrow(kex),
           upstream_exon_end = kex[exon == min(dup) - 1, end], downstream_exon_start = kex[exon == max(dup) + 1, start])]
vm <- fread(file.path(INT, "vaf", "clone_markers.tsv"), colClasses = list(character = "sample"))[
  protein == "p.Gly543Cys", .(sample, g543_alt = n_alt, g543_depth = depth)]
blk <- merge(SMP[, .(sample, label, case)], merge(blk, vm, by = "sample", all.x = TRUE), by = "sample")[match(SMP$sample, sample)]
num <- c("l2", "se", "flank_l2", "copy_ratio", "ci_low", "ci_high", "flank_copy_ratio")
blk[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(blk, file.path(OUT, "kmt2a_block.tsv"), sep = "\t")
fwrite(cr[chr == km$chr & start >= B0 - 80e3 & end <= B1 + 80e3, .(sample, start, end, l2 = round(l2, 3))], file.path(OUT, "kmt2a_bins.tsv"), sep = "\t")
cat(sprintf("\nKMT2A block %s:%d-%d (%.1f kb), exons %d-%d of %s\n", km$chr, B0, B1, (B1 - B0 + 1) / 1e3, min(dup), max(dup), km$tx))
print(blk[, .(label, bins, copy_ratio, ci_low, ci_high, flank_copy_ratio, g543 = paste0(g543_alt, "/", g543_depth))])

## ---- TET2 and SETBP1 in the TET2 case -------------------------------------------------
TC  <- SMP[case == CASES[2], sample]          # the TET2 case
SNP <- data.table(gene = c("TET2", "SETBP1"), snp = c("rs2454206", "rs3744825"), chr = c("chr4", "chr18"), pos = c(105275794L, 44952641L))
WIN <- 250e3
reg <- merge(SNP, gene_spans[, .(gene, g_start = start, g_end = end)], by = "gene")[match(SNP$gene, gene)]
reg[, `:=`(w0 = g_start - WIN, w1 = g_end + WIN)]
ac_file <- file.path(OUT, "tet2_case_allelic_counts.tsv.gz")
if (dir.exists(file.path(CNV, "allelic"))) {                 # extract once from the full allelic counts
  cond <- paste(sprintf('($1=="%s" && $2>=%.0f && $2<=%.0f)', reg$chr, reg$w0, reg$w1), collapse = " || ")
  ac <- rbindlist(lapply(TC, function(s) {
    f <- file.path(CNV, "allelic", paste0(s, ".allelicCounts.tsv")); inputs <<- c(inputs, f)
    fread(cmd = sprintf("awk -F'\\t' '%s' '%s'", cond, f), header = FALSE,
          col.names = c("chr", "pos", "ref_n", "alt_n", "ref_base", "alt_base"))[, sample := s]
  }))
  fwrite(ac, ac_file, sep = "\t")
}
ac <- if (file.exists(ac_file)) fread(cmd = sprintf("gzip -dc '%s'", ac_file), colClasses = list(character = "sample")) else NULL
if (!is.null(ac)) inputs <- c(inputs, ac_file)
tg <- rbindlist(lapply(seq_len(nrow(reg)), function(i) {
  r <- reg[i]
  rbindlist(lapply(TC, function(s) {
    x <- cr[.(s, r$chr)][start >= r$g_start & end <= r$g_end, l2]
    h <- het[sample == s & chr == r$chr & pos >= r$g_start - 50e3 & pos <= r$g_end + 50e3 & n >= 6 & n <= 500]
    a <- if (is.null(ac)) NULL else ac[sample == s & chr == r$chr & pos == r$pos]
    data.table(gene = r$gene, sample = s, gene_start = r$g_start, gene_end = r$g_end, gene_bins = length(x),
               gene_copy_ratio = round(2^mean(x), 3), lowest_10kb_copy_ratio = round(2^min(frollmean(x, 10), na.rm = TRUE), 3),
               het_sites = nrow(h), mean_abs_baf_dev = round(mean(abs(h$baf - 0.5)), 3), expected_dev = round(mean(ea[h$n]), 3),
               genome_baf_dev = summ[sample == s, baf_dev], genome_expected_dev = summ[sample == s, baf_dev_expected],
               snp = r$snp, snp_pos = r$pos,
               snp_alt_reads = if (is.null(a) || !nrow(a)) NA_character_ else sprintf("%d/%d", a$alt_n, a$ref_n + a$alt_n))
  }))
}))
fwrite(tg, file.path(OUT, "tet2_case_genes.tsv"), sep = "\t")
cat("\nTET2 case, TET2 and SETBP1:\n"); print(tg)

## ---- figure 1: the genome ------------------------------------------------------------
len <- cr[, .(len = max(end)), by = chr][order(chr)]
len[, off := c(0, cumsum(as.numeric(len))[-.N])]
gx <- function(ch, p) len$off[match(as.character(ch), as.character(len$chr))] + p
lab_it <- function(d) d[, `:=`(lab = factor(SMP$label[match(sample, SMP$sample)], levels = LV), trk = factor(trk, levels = TR))]
c1 <- cr[, .(l2 = median(l2)), by = .(sample, chr, w = start %/% 1e6)][, `:=`(x = gx(chr, (w + 0.5) * 1e6), y = pmin(2^l2, 2.35), trk = TR[1])]
set.seed(16)
h1 <- het[n >= 8][, .SD[base::sample(.N, min(.N, 25000L))], by = sample][, `:=`(x = gx(chr, pos), y = baf, trk = TR[2])]
s1 <- mf[, .(sample, x0 = gx(chr, start), x1 = gx(chr, end), cr = pmin(2^LOG2_COPY_RATIO_POSTERIOR_50, 2.35),
             maf = MINOR_ALLELE_FRACTION_POSTERIOR_50, naf = NUM_POINTS_ALLELE_FRACTION)]
sg <- rbind(s1[, .(sample, x0, x1, y = cr, trk = TR[1])],
            s1[naf >= 20, .(sample, x0, x1, y = maf, trk = TR[2])], s1[naf >= 20, .(sample, x0, x1, y = 1 - maf, trk = TR[2])])
m1 <- cs[, .(sample, x = gx(chr, (start + end) / 2), y = pmax(pmin(copy_ratio, 2.25), 0.08), call,
             kind = ifelse(class == ACQ, "acquired (clone libraries only)", "inherited or technical"),
             tag = ifelse(class == ACQ, sub(",.*", "", genes), ""), trk = TR[1])]
lab_it(c1); lab_it(h1); lab_it(sg); lab_it(m1)
sg <- merge(sg, SMP[, .(sample, case)], by = "sample")
lv_f <- factor(LV, levels = LV)
hl <- rbind(CJ(lab = lv_f, y = c(0.5, 1, 1.5))[, `:=`(trk = TR[1], lt = ifelse(y == 1, "solid", "dotted"))],
            CJ(lab = lv_f, y = 0.5)[, `:=`(trk = TR[2], lt = "solid")])
lims <- rbind(CJ(lab = lv_f, y = c(0, 2.4))[, trk := TR[1]], CJ(lab = lv_f, y = c(0, 1))[, trk := TR[2]])
hl[, trk := factor(trk, levels = TR)]; lims[, trk := factor(trk, levels = TR)]
mid <- len[, .(x = off + len / 2, chr = sub("chr", "", chr))]
shade <- len[seq(2, .N, 2), .(xmin = off, xmax = off + len)]
p1 <- ggplot() +
  geom_rect(data = shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey94") +
  geom_blank(data = lims, aes(x = 0, y = y)) +
  geom_hline(data = hl, aes(yintercept = y, linetype = lt), colour = "grey55", linewidth = 0.3) +
  geom_point(data = h1, aes(x, y), size = 0.05, alpha = 0.18, colour = "grey15") +
  geom_point(data = c1, aes(x, y), size = 0.35, alpha = 0.6, colour = "grey15") +
  geom_segment(data = sg, aes(x = x0, xend = x1, y = y, yend = y, colour = case), linewidth = 0.9) +
  geom_point(data = m1, aes(x, y, shape = call, fill = kind), size = 2.6, stroke = 0.3, colour = "grey10") +
  geom_text(data = m1[tag != ""], aes(x, y, label = tag), size = 3.3, hjust = -0.3, fontface = "italic", colour = "#B2182B") +
  facet_grid(lab + trk ~ ., scales = "free_y", switch = "y") +
  scale_linetype_identity() +
  scale_colour_manual(values = CASE_COL, guide = "none") +
  scale_shape_manual(values = c(`+` = 24, `-` = 25), labels = c(`+` = "called gain", `-` = "called loss"), name = NULL) +
  scale_fill_manual(values = c(`acquired (clone libraries only)` = "#B2182B", `inherited or technical` = "grey75"), name = NULL) +
  scale_x_continuous(breaks = mid$x, labels = mid$chr, expand = c(0.004, 0)) +
  scale_y_continuous(breaks = function(l) if (max(l) > 1.5) c(0.5, 1, 1.5, 2) else c(0, 0.5, 1)) +
  guides(shape = guide_legend(override.aes = list(fill = "grey60")), fill = guide_legend(override.aes = list(shape = 24))) +
  labs(x = NULL, y = NULL, title = "Copy number and allelic balance across the genome, every library",
       subtitle = paste0("Copy ratio: dots = median of the 1-kb bins in each 1-Mb window; lines = GATK segments (1 = two copies, ",
                         "1.5 = one extra copy in every cell, 0.5 = one copy lost in every cell).\nAllele fraction: dots = heterozygous ",
                         "common SNPs (25,000 per library); lines = the segments' minor-allele fraction and its mirror image.")) +
  theme_bw(base_size = 11) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9), strip.placement = "outside",
        strip.background = element_rect(fill = "grey97", colour = "grey80"), panel.grid = element_blank(),
        panel.spacing.y = unit(0.6, "mm"), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 8.5),
        legend.position = "bottom", plot.subtitle = element_text(size = 9.5))
ggsave(file.path(OUT, "cnv_genome.png"), p1, width = 16, height = 13, dpi = 170)

## ---- figure 2: KMT2A -------------------------------------------------------------------
X0 <- B0 - 65e3; X1 <- B1 + 75e3
kd <- merge(cr[chr == km$chr & start >= X0 & end <= X1, .(sample, start, end, l2)], SMP[, .(sample, label, case)], by = "sample")
kd[, `:=`(x = (start + end) / 2e3, y = pmin(2^l2, 2.8), lab = factor(label, levels = LV))]
setorder(kd, sample, start)
kd[, run := frollmean(2^l2, 5, align = "center"), by = sample]
pa <- ggplot(kd, aes(x, y)) +
  annotate("rect", xmin = B0 / 1e3, xmax = B1 / 1e3, ymin = -Inf, ymax = Inf, fill = "#B2182B", alpha = 0.10) +
  geom_hline(yintercept = 1, colour = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = c(1.5, 2), colour = "grey70", linewidth = 0.3, linetype = "dotted") +
  geom_point(aes(colour = case), size = 0.9, alpha = 0.75) +
  geom_line(aes(y = run), colour = "grey20", linewidth = 0.45, na.rm = TRUE) +
  facet_grid(lab ~ ., switch = "y") +
  scale_colour_manual(values = CASE_COL, guide = "none") +
  scale_x_continuous(limits = c(X0, X1) / 1e3, expand = c(0, 0)) +
  scale_y_continuous(breaks = c(1, 2), limits = c(0, 2.9)) +
  labs(x = NULL, y = "copy ratio (1-kb bins)", title = sprintf("KMT2A: a %.0f-kb block at about twice the normal level", (B1 - B0 + 1) / 1e3),
       subtitle = "dots = 1-kb bins; line = running mean of 5 bins; shaded = the block called in both clone libraries") +
  theme_bw(base_size = 11) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9), strip.placement = "outside",
        strip.background = element_rect(fill = "grey97", colour = "grey80"), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.subtitle = element_text(size = 9.5))
ex <- kex[end >= X0 & start <= X1]
ex[, `:=`(x0 = pmax(start, X0) / 1e3, x1 = pmin(end, X1) / 1e3)]
el <- ex[exon %in% unique(c(1L, min(dup), max(dup), max(dup) + 1L)), .(x = (x0 + x1) / 2, exon)]
pg <- ggplot() +
  annotate("segment", x = max(km$start, X0) / 1e3, xend = min(km$end, X1) / 1e3, y = 0, yend = 0, colour = "grey30", linewidth = 0.5) +
  geom_rect(data = ex, aes(xmin = x0, xmax = pmax(x1, x0 + 0.35), ymin = -0.32, ymax = 0.32), fill = "#2F5D8A") +
  geom_text(data = el, aes(x, 0.66, label = exon), size = 3.2) +
  annotate("text", x = X0 / 1e3 + 1, y = 0.66, hjust = 0, size = 3.2, label = "exon") +
  annotate("text", x = X0 / 1e3 + 1, y = -0.8, hjust = 0, size = 3.2, fontface = "italic",
           label = sprintf("KMT2A, MANE Select %s", sub("\\..*", "", km$tx))) +
  annotate("errorbar", xmin = B0 / 1e3, xmax = B1 / 1e3, y = -0.8, width = 0.3, orientation = "y", colour = "#B2182B", linewidth = 0.7) +
  annotate("text", x = (B0 + B1) / 2e3, y = -1.3, size = 3.3, colour = "#B2182B", label = sprintf("duplicated: exons %d-%d", min(dup), max(dup))) +
  scale_x_continuous(limits = c(X0, X1) / 1e3, expand = c(0, 0), labels = function(v) format(v, big.mark = ",")) +
  scale_y_continuous(limits = c(-1.55, 0.95)) +
  labs(x = sprintf("position on %s (kb, hg38)", km$chr), y = NULL) +
  theme_bw(base_size = 11) + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
bd <- copy(blk)[, g := ifelse(case == CASES[1] & !is.na(g543_depth), sprintf("\nDNMT3A G543C %d/%d", g543_alt, g543_depth), "")]
bd[, lab := factor(paste0(label, g), levels = rev(paste0(label, g)))]
pb <- ggplot(bd, aes(copy_ratio, lab, colour = case)) +
  geom_vline(xintercept = 1, colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = c(1.5, 2), colour = "grey70", linewidth = 0.3, linetype = "dotted") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.25, orientation = "y") +
  geom_point(size = 2.6) +
  scale_colour_manual(values = CASE_COL, name = NULL) +
  scale_x_continuous(limits = c(0.7, 2.3), breaks = c(1, 1.5, 2)) +
  labs(x = sprintf("copy ratio of the block\n(mean of %d bins, 95%% interval)", blk$bins[1]), y = NULL, title = "The block in each library",
       subtitle = "under each name: DNMT3A clone-mutation reads") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 9.5))
p2 <- ((pa / pg + plot_layout(heights = c(7, 1.5))) | pb) + plot_layout(widths = c(2.1, 1))
ggsave(file.path(OUT, "cnv_kmt2a.png"), p2, width = 15, height = 9.5, dpi = 200)

## ---- figure 3: TET2 and SETBP1 in the TET2 case --------------------------------------
TL <- LV[SMP$sample %in% TC]
td <- rbindlist(lapply(seq_len(nrow(reg)), function(i) { r <- reg[i]
  cr[sample %in% TC & chr == r$chr & start >= r$w0 & end <= r$w1, .(sample = as.character(sample), gene = r$gene, x = (start + end) / 2e6, l2)] }))
setorder(td, sample, gene, x)
td[, `:=`(y = pmin(2^l2, 2.4), run = frollapply(2^l2, 11, median, align = "center")), by = .(sample, gene)]
th <- rbindlist(lapply(seq_len(nrow(reg)), function(i) { r <- reg[i]
  het[sample %in% TC & chr == r$chr & pos >= r$w0 & pos <= r$w1 & n >= 6, .(sample, gene = r$gene, x = pos / 1e6, y = baf)] }))
lab3 <- function(d, tk) d[, `:=`(lab = factor(SMP$label[match(sample, SMP$sample)], levels = TL), trk = factor(tk, levels = TR),
                                 gene = factor(gene, levels = reg$gene))]
lab3(td, TR[1]); lab3(th, TR[2])
gfac <- function(g) factor(g, levels = reg$gene)
gs <- reg[, .(gene = gfac(gene), xmin = g_start / 1e6, xmax = g_end / 1e6)]
sn <- reg[, .(gene = gfac(gene), x = pos / 1e6, snp)]
sl <- sn[, .(gene, x, snp, lab = factor(TL[1], levels = TL), trk = factor(TR[1], levels = TR))]
h3 <- rbind(CJ(lab = factor(TL, levels = TL), y = 1)[, trk := TR[1]], CJ(lab = factor(TL, levels = TL), y = 0.5)[, trk := TR[2]])
h3[, trk := factor(trk, levels = TR)]
l3 <- rbind(CJ(lab = factor(TL, levels = TL), gene = gfac(reg$gene), y = c(0, 2.5))[, trk := TR[1]],
            CJ(lab = factor(TL, levels = TL), gene = gfac(reg$gene), y = c(0, 1))[, trk := TR[2]])
l3[, `:=`(trk = factor(trk, levels = TR), x = reg$pos[match(as.character(gene), reg$gene)] / 1e6)]
p3 <- ggplot() +
  geom_rect(data = gs, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "#D95F02", alpha = 0.10) +
  geom_vline(data = sn, aes(xintercept = x), colour = "#B2182B", linetype = "dashed", linewidth = 0.4) +
  geom_blank(data = l3, aes(x, y)) +
  geom_hline(data = h3, aes(yintercept = y), colour = "grey50", linewidth = 0.3) +
  geom_point(data = td, aes(x, y), size = 0.35, alpha = 0.35, colour = "grey20") +
  geom_line(data = td, aes(x, run), colour = "#D95F02", linewidth = 0.6, na.rm = TRUE) +
  geom_point(data = th, aes(x, y), size = 0.8, alpha = 0.6, colour = "grey15") +
  geom_text(data = sl, aes(x, 2.42, label = snp), size = 3.1, hjust = -0.08, colour = "#B2182B") +
  facet_grid(lab + trk ~ gene, scales = "free", switch = "y") +
  scale_x_continuous(labels = function(v) sprintf("%.2f", v)) +
  labs(x = "position (Mb, hg38)", y = NULL, title = "TET2 case: no copy loss at TET2 or SETBP1",
       subtitle = paste0("Copy ratio: dots = 1-kb bins, line = running median of 11 bins. Allele fraction: heterozygous common SNPs ",
                         "(>= 6 reads).\nShaded = the gene; dashed = the panel SNP that looked imbalanced in the DDHN.")) +
  theme_bw(base_size = 11) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9), strip.placement = "outside",
        strip.background = element_rect(fill = "grey97", colour = "grey80"), panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 9.5))
ggsave(file.path(OUT, "cnv_tet2_case_genes.png"), p3, width = 14, height = 8.2, dpi = 200)

## ---- the run: bins before and after filtering (from the GATK log, when copied) ------------
fl <- file.path(CNV, "logs", "filter_bins.log")
bl <- if (file.exists(fl)) regmatches(readLines(fl), regexpr("[0-9]+ / [0-9]+ intervals remain", readLines(fl))) else character(0)
bl <- if (length(bl)) as.numeric(strsplit(sub(" intervals remain", "", tail(bl, 1)), " / ")[[1]]) else c(NA, NA)
fwrite(data.table(bins_total = bl[2], bins_after_filter = bl[1], bins_used = summ[sample == SMP[normal == FALSE, sample][1], bins]), file.path(OUT, "cnv_run.tsv"), sep = "\t")
if (file.exists(fl)) inputs <- c(inputs, fl)

## ---- provenance -------------------------------------------------------------------------
fi <- file.info(unique(inputs))
fwrite(data.table(file = sub(paste0(DATA_ROOT, "/"), "", rownames(fi), fixed = TRUE), bytes = fi$size,
                  modified = format(fi$mtime, "%Y-%m-%d %H:%M")), file.path(OUT, "cnv_inputs.tsv"), sep = "\t")
message("written to ", OUT, " (input: ", CNV, ")")
