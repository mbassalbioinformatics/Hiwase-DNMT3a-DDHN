#!/usr/bin/env Rscript
###############################################################################
# 26_kmt2a_ptd_rnaseq.R                                              2026-10-19
#
# Can the RNA-seq see the KMT2A duplication? A tandem duplication of exons 2-10
# makes transcripts in which exon 10 is spliced back to exon 2. On the server, every
# RNA-seq read of the four libraries was searched for 30-base probes (15 + 15
# bases) across all 66 backward junctions (exon i -> exon j <= i, i = 2-12) and
# the 12 normal junctions of exons 1-13 (MANE Select ENST00000534358.8):
#   03_rnaseq/server/04_kmt2a_ptd_junctions.sh   the aligned reads (BAMs)
#   03_rnaseq/server/05_kmt2a_ptd_fastq.sh       all reads (trimmed FASTQs)
# This script summarises the FASTQ search. Reads across a junction can be PCR
# copies of one RNA molecule, so each junction's reads are also counted as
# molecules: distinct read starts relative to the junction (read file, strand,
# position of the probe in the read), one read per fragment. P(no PTD read): the
# chance of seeing no backward junction among n molecules across the exon-10
# splice donor if half of KMT2A transcripts carried the duplication (a PTD
# transcript has two exon-10 donors, one to exon 2, one to exon 11: 1/3 of donor
# reads would be PTD).
#
# Output: integration/cnv/kmt2a_ptd_rnaseq.tsv (per library),
#         integration/cnv/kmt2a_ptd_rnaseq_junctions.tsv (per library and junction)
###############################################################################
source("config/paths.R")
ID <- with(data.table::fread(SAMPLE_SHEET, colClasses = "character"), setNames(sample, key))   # role key -> sample ID (config/samples.tsv)

suppressPackageStartupMessages(library(data.table))
PTD <- file.path(RNA, "09_kmt2a_ptd"); FQ <- file.path(PTD, "fastq_search")
OUT <- file.path(INT, "cnv")
SAMPLES <- c(ID[["D30"]], ID[["Y35"]], ID[["DX"]], ID[["D23"]])
smp <- fread(SAMPLE_SHEET, colClasses = "character")[methylation == "analysed"]

summ   <- fread(file.path(FQ, "fastq_ptd_summary.tsv"), colClasses = list(character = "sample"))
counts <- fread(file.path(FQ, "fastq_junction_counts.tsv"), colClasses = list(character = "sample"))
probes <- fread(file.path(PTD, "probes.tsv"))
rc <- function(s) vapply(s, function(x) paste(rev(strsplit(chartr("ACGTN", "TGCAN", x), "")[[1]]), collapse = ""), "")

# molecules: distinct read starts relative to the junction, one read per fragment
mol <- rbindlist(lapply(SAMPLES, function(s) {
  r <- fread(file.path(FQ, s, "supporting_reads.tsv"), colClasses = "character")
  if (!nrow(r)) return(NULL)
  r <- merge(r, probes[, .(probe, p = toupper(sequence), q = toupper(reverse_complement))], by = "probe")
  r[, `:=`(seq = toupper(sequence), fwd = mapply(function(a, b) regexpr(b, a, fixed = TRUE)[1], toupper(sequence), p))]
  r[, rev_ := mapply(function(a, b) regexpr(b, a, fixed = TRUE)[1], seq, q)]
  r[, `:=`(strand = fifelse(fwd > 0, "+", "-"), start = fifelse(fwd > 0, fwd, nchar(seq) - rev_ - nchar(q) + 2L))]
  r <- r[order(file)][!duplicated(paste(probe, fragment))]          # one read per fragment
  r[, .(sample = s, molecules = uniqueN(paste(file, strand, start))), by = probe]
}))
jt <- merge(counts[, .(sample, probe, type, donor_exon, acceptor_exon, fragments, distinct_sequences)], mol,
            by = c("sample", "probe"), all.x = TRUE)[is.na(molecules), molecules := 0L]
fwrite(jt[order(sample, type, donor_exon, acceptor_exon)], file.path(OUT, "kmt2a_ptd_rnaseq_junctions.tsv"), sep = "\t")

e <- new.env(); load(file.path(RNA, "05_processed_counts_genes", "dge_gene_level.RData"), envir = e)
m10 <- jt[probe == "normal_e10_e11"]
res <- summ[, .(sample, reads_searched, ptd_fragments = PTD_fragments_all_probes)]
res <- merge(res, m10[, .(sample, normal_e10_e11_fragments = fragments, normal_e10_e11_molecules = molecules)], by = "sample")
res[, `:=`(p_no_ptd_read_if_half = (2 / 3)^normal_e10_e11_molecules,
           kmt2a_cpm = round(e$cpm_mat["KMT2A", sample], 1),
           name = smp$plot_label[match(sample, smp$sample)])]
res <- res[match(SAMPLES, sample)]
fwrite(res, file.path(OUT, "kmt2a_ptd_rnaseq.tsv"), sep = "\t")
print(res)
message("written to ", OUT)
