#!/usr/bin/env Rscript
###############################################################################
# 15_chimerism.R                                                    2026-10-13
#
# How much of each post-transplant sample could come from the patient's own
# (recipient) cells? The libraries hold only donor-derived genomes by genotype
# (Section 1); this measures how small any recipient share must be.
#
# Probe sites: CpG-destroying SNPs of the identity check, strand by strand
# (Watson reads give the C of the CpG, Crick reads the G; methylation never
# makes the "other" base, so other / (C + T + other) on the informative strand
# is the allele fraction). A site is informative for person X when X and the
# other person in the study (Y, unrelated) are opposite homozygotes. At such
# common SNPs an unrelated recipient carries the non-donor allele on about half
# of their chromosomes. Only sites where X has no "other" read at all (homozygous
# reference) are used: at homozygous-alternative sites about 10% of reads still
# show C/T, so the reverse direction has far too high a background. So
#     non-donor read fraction  ~  error + 0.5 x recipient share
# X's genotype is called without the library being tested (leave-one-out), so
# recipient reads in that library cannot hide a site. The error rate comes from
# the pure donor libraries (donor at donation and in 2023), which contain no
# recipient DNA at all.
#
# Inputs  emseq/00_sample_identity/<sample>.strand_counts.tsv (02_allele_balance_counts.py)
# Output  integration/vaf/chimerism.tsv
###############################################################################
suppressPackageStartupMessages(library(data.table))
source("config/paths.R")
OUT <- file.path(INT, "vaf")
# libraries from the sample sheet; the pure donor libraries (role "Donor, ...") give the error rate
SMP <- fread(SAMPLE_SHEET, colClasses = "character")[methylation == "analysed",
         .(sample, person, pure_donor = grepl("^Donor", role), label = plot_label)]
MIN_LOO <- 10L     # reads on the informative strand to call X homozygous reference without the tested library
MIN_Y   <- 12L     # reads in the other person's pooled libraries
Y_ALT   <- 0.8     # other-base fraction for Y to count as homozygous alternative (hom-alt sites read ~0.85-0.98)

## ---- per-strand loci ---------------------------------------------------------
x <- rbindlist(lapply(SMP$sample, function(s) {
  d <- fread(file.path(EMS, "00_sample_identity", paste0(s, ".strand_counts.tsv")))
  rbind(d[, .(sample = s, locus = paste0(chr, ":", pos, ":W"), ct = w_ct, oth = w_other)],
        d[, .(sample = s, locus = paste0(chr, ":", pos + 1L, ":C"), ct = c_ct, oth = c_other)])
}))
x <- merge(x, SMP[, .(sample, person)], by = "sample")
x[, n := ct + oth]

## ---- genotype of a person from a set of libraries -----------------------------
pool <- function(d) d[, .(n = sum(n), oth = sum(oth)), by = locus]

res <- rbindlist(lapply(seq_len(nrow(SMP)), function(i) {
  s <- SMP$sample[i]; X <- SMP$person[i]; Y <- setdiff(c("A", "B"), X)
  gy <- pool(x[person == Y])[n >= MIN_Y & oth / n >= Y_ALT, .(locus)]               # the other person: alternative allele
  gx <- pool(x[person == X & sample != s])[n >= MIN_LOO & oth == 0, .(locus)]       # X homozygous reference, leave-one-out
  inf <- merge(gx, gy, by = "locus")
  d <- merge(x[sample == s, .(locus, ct, oth, n)], inf, by = "locus")[n > 0]
  d[, nondonor := oth]                                                              # reads with the allele the donor lacks
  data.table(sample = s, person = X, sites = nrow(d), reads = sum(d$n), nondonor_reads = sum(d$nondonor),
             nondonor_fraction = sum(d$nondonor) / sum(d$n))
}))
res <- merge(SMP, res, by = c("sample", "person"), sort = FALSE)

## ---- error rate from the pure donor libraries; recipient share -----------------
err <- res[pure_donor == TRUE, sum(nondonor_reads) / sum(reads)]
res[, excess := nondonor_fraction - err]
# recipient share = 2 x excess; 95% upper bound from the exact binomial upper limit of the non-donor reads
res[, upper95 := pmax(0, 2 * (qbeta(0.975, nondonor_reads + 1, reads - nondonor_reads) - err))]
res[, recipient_share := pmax(0, 2 * excess)]
res[, binom_p := mapply(function(k, n) if (n > 0) binom.test(k, n, err, alternative = "greater")$p.value else NA_real_, nondonor_reads, reads)]
fwrite(res, file.path(OUT, "chimerism.tsv"), sep = "\t")
print(res[, .(label, sites, reads, nondonor_reads, nondonor_pct = round(100 * nondonor_fraction, 3),
              error_pct = round(100 * err, 3), recipient_share_pct = round(100 * recipient_share, 2),
              upper95_pct = round(100 * upper95, 2), p = signif(binom_p, 2))])
message("error rate (pure donor libraries): ", signif(100 * err, 3), "%")
