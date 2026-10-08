#!/usr/bin/env Rscript
###############################################################################
# 13_vaf_from_taps_reads.R                                           2026-10-08
#
# Variant allele fractions (VAF) of clonal-haematopoiesis driver mutations, read
# straight from the methylation (TAPS) libraries, so the clone can be followed in
# every methylation sample without separate DNA sequencing.
#
# Input: base counts per original DNA strand at every position of the coding
# exons (+/- 10 bp) of an 85-gene myeloid panel, from samtools mpileup of each
# Msuite2 alignment (server: 01_methylation_server/06_vaf_pileup.sh + vaf_pileup_counts.py;
# MAPQ >= 10, base quality >= 20).
#
# TAPS turns a methylated C into T on the strand that carries it: on the
# original top strand (OT) a T read may be a methylated C, and on the original
# bottom strand (OB) an A read (top-strand view) may be a methylated G. In blood
# this happens at CpGs (non-CpG conversion here ~0.25%, the error level), so an
# allele is counted on a strand only where TAPS can neither create nor hide it:
#   alt A   OT; OB too unless a G allele could sit in a CpG here
#   alt T   OB; OT too unless a C allele could sit in a CpG here
#   alt C   OB; OT too unless the next base can be G (alt C would be a CpG C)
#   alt G   OT; OB too unless the previous base can be C
# "Could": from the reference, or an allele seen on its clean strand (G on OT,
# C on OB) in any sample, so a SNP that creates a CpG next door is caught.
# VAF = alt reads / all reads on the strand(s) used.
#
# A sample carries a variant with >= 3 alt reads, VAF >= 10% and >= 8 reads.
# Germline or acquired is judged within a person (same genome): a variant whose
# VAF differs between samples of one person (Fisher exact test, p < 0.01) is an
# acquired (somatic) candidate, or allelic imbalance (LOH) when it is a germline
# SNP in another person; otherwise pooled VAF >= 85% = homozygous, >= 30% =
# heterozygous. A difference that reaches >= 85% in one sample, or sits on a
# known common SNP (KNOWN_SNPS, checked in dbSNP/ClinVar), is also allelic
# imbalance. Flags for artefacts: more than one alt allele at the position,
# another variant within 10 bp, a short tandem repeat or homopolymer around it,
# alt reads on one strand only (none on the other strand with >= 5 reads, or
# Fisher p < 0.05), low-level calls in more than one person, and support from a
# single sample with only 3 alt reads. Variants private to one sample are listed
# separately, with their read support per strand.
# Protein changes from the MANE Select v1.0 CDS.
# Indels are not counted (the pileup parser skips them), so frameshifts and
# other small insertions/deletions are not seen here.
#
# Output in integration/vaf/: strand_check.tsv, coverage_by_gene.tsv,
# variants_all.tsv, variants_acquired.tsv, clone_markers.tsv, hotspots.tsv,
# private_candidates.tsv, dnmt3a_tet2_positions_with_alt_reads.tsv, vaf_clones.png
###############################################################################
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38); library(ggplot2)
})
source("config/paths.R")
DMRD <- file.path(EMS, "03_dmr_analysis")
PANEL <- "01_methylation_server/vaf_panel_hg38.bed"
OUT  <- file.path(INT, "vaf"); dir.create(OUT, showWarnings = FALSE)
G    <- BSgenome.Hsapiens.UCSC.hg38

MIN_ALT <- 3L; MIN_VAF <- 0.10; MIN_DEPTH <- 8L; P_ACQ <- 0.01; P_STRAND <- 0.05

# common germline SNPs among the variants whose allele fraction differs within a person (looked up 2026-10-08)
KNOWN_SNPS <- data.table(chr = c("chr4", "chr18"), pos = c(105275794L, 44952641L), ref = c("A", "G"), alt = c("G", "A"),
                         known_snp = c("rs2454206 (TET2 p.Ile1762Val, common polymorphism)",
                                       "rs3744825 (SETBP1 p.Val1101Ile, gnomAD AF ~0.11-0.13, ClinVar benign)"))

# libraries in sample-sheet order (timeline order within each case); person = genome (genotype-confirmed)
SMP <- fread(SAMPLE_SHEET, colClasses = "character")[methylation %in% c("analysed", "extra"),
         .(sample, person, label = plot_label, case = fifelse(case == "none", "Not in design", paste(case, "case")))]
SMP[, order := .I]
CASE_LV  <- unique(SMP$case)
CASE_COL <- setNames(c("#6A51A3", "#D95F02", "#1B9E77")[seq_along(CASE_LV)], CASE_LV)
CASE_COL[grepl("^Not in design", CASE_LV)] <- "grey50"
# two-line labels for the figure x axis (horizontal text, no clipping): each label broken near its middle
half <- function(l) { w <- strsplit(l, " ")[[1]]; k <- ceiling(length(w) / 2)
  if (length(w) < 2) l else paste(paste(w[1:k], collapse = " "), paste(w[-(1:k)], collapse = " "), sep = "\n") }
SMP[, axis_label := vapply(label, half, "")]

## ---- counts, reference base and context -------------------------------------
x <- fread(cmd = sprintf("gzip -dc '%s'", file.path(DMRD, "vaf", "ALL.strand_counts.tsv.gz")), header = FALSE,
           col.names = c("sample", "chr", "pos", "OT_A", "OT_C", "OT_G", "OT_T", "OB_A", "OB_C", "OB_G", "OB_T"),
           colClasses = list(character = 1))
x <- x[sample %in% SMP$sample]
pos <- unique(x[, .(chr, pos)])
tri <- as.character(getSeq(G, GRanges(pos$chr, IRanges(pos$pos - 1L, pos$pos + 1L))))
pos[, `:=`(prev = substr(tri, 1, 1), ref = substr(tri, 2, 2), nxt = substr(tri, 3, 3))]
x <- merge(x, pos, by = c("chr", "pos"))
x <- x[ref %in% c("A", "C", "G", "T")]
x[, `:=`(OT = OT_A + OT_C + OT_G + OT_T, OB = OB_A + OB_C + OB_G + OB_T)]
message(format(nrow(x), big.mark = ","), " sample-positions over ", format(nrow(pos), big.mark = ","), " panel positions")

panel <- fread(PANEL, header = FALSE, col.names = c("chr", "start", "end", "gene"))
pg <- foverlaps(pos[, .(chr, s = pos - 1L, e = pos, pos)], setkey(copy(panel), chr, start, end),
                by.x = c("chr", "s", "e"), nomatch = 0L, mult = "first")[, .(chr, pos, gene)]

## ---- strand check: is the library directional, and is OT/OB assigned right? --
# At a CpG C, OT carries the methylation call (T = methylated) and OB the
# genotype (T only by sequencing error or a SNP); at a CpG G the reverse.
x[, `:=`(cpg_c = ref == "C" & nxt == "G", cpg_g = ref == "G" & prev == "C")]
sc <- x[, .(
  cpgC_OT_T_pct    = 100 * sum(OT_T[cpg_c]) / sum(OT_C[cpg_c] + OT_T[cpg_c]),
  cpgC_OB_T_pct    = 100 * sum(OB_T[cpg_c]) / sum(OB_C[cpg_c] + OB_T[cpg_c]),
  cpgG_OB_A_pct    = 100 * sum(OB_A[cpg_g]) / sum(OB_G[cpg_g] + OB_A[cpg_g]),
  cpgG_OT_A_pct    = 100 * sum(OT_A[cpg_g]) / sum(OT_G[cpg_g] + OT_A[cpg_g]),
  nonCpG_C_OT_T_pct = 100 * sum(OT_T[ref == "C" & !cpg_c]) / sum(OT_C[ref == "C" & !cpg_c] + OT_T[ref == "C" & !cpg_c]),
  median_reads_OT = as.double(median(OT)), median_reads_OB = as.double(median(OB))), by = sample]
ori <- fread(file.path(DMRD, "vaf", "read1_orientation_DNMT3A_region.tsv"), header = FALSE, colClasses = list(character = 1))
sc <- merge(sc, ori[, .(sample = V1, read1_forward = V3, read1_reverse = V5)], by = "sample")
sc <- sc[SMP[, .(sample, order)], on = "sample"][order(order)][, order := NULL]
num <- names(sc)[vapply(sc, is.double, TRUE)]
sc[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(sc, file.path(OUT, "strand_check.tsv"), sep = "\t")
print(sc)

## ---- allele counts on the strand(s) that report genotype ---------------------
# alleles that may be present: G read on OT and C read on OB are never TAPS products
al <- x[, .(c_alt = any(OB_C >= 2 & OB_C >= 0.1 * OB), g_alt = any(OT_G >= 2 & OT_G >= 0.1 * OT)), by = .(chr, pos)]
pos <- merge(pos, al, by = c("chr", "pos"))
pos[, `:=`(c_here = ref == "C" | c_alt, g_here = ref == "G" | g_alt)]
pos <- merge(pos, pos[, .(chr, pos = pos - 1L, g_next_seen = g_here)], by = c("chr", "pos"), all.x = TRUE)
pos <- merge(pos, pos[, .(chr, pos = pos + 1L, c_prev_seen = c_here)], by = c("chr", "pos"), all.x = TRUE)
pos[, `:=`(g_next = fifelse(is.na(g_next_seen), nxt == "G", g_next_seen),
           c_prev = fifelse(is.na(c_prev_seen), prev == "C", c_prev_seen))]
x <- merge(x, pos[, .(chr, pos, c_here, g_here, g_next, c_prev)], by = c("chr", "pos"))
long <- rbindlist(lapply(c("A", "C", "G", "T"), function(b) {
  d <- x[ref != b, .(sample, chr, pos, ref, alt = b, ot_alt = get(paste0("OT_", b)), ob_alt = get(paste0("OB_", b)),
                     OT, OB, c_here, g_here, g_next, c_prev)]
  d[, `:=`(use_ot = switch(b, A = TRUE, G = TRUE, T = !(c_here & g_next), C = !g_next),
           use_ob = switch(b, T = TRUE, C = TRUE, A = !(g_here & c_prev), G = !c_prev))]
  d[, `:=`(ot_alt = use_ot * ot_alt, ot_depth = use_ot * OT, ob_alt = use_ob * ob_alt, ob_depth = use_ob * OB)]
  d[, `:=`(n_alt = ot_alt + ob_alt, depth = ot_depth + ob_depth,
           strands = fifelse(use_ot & use_ob, "both", fifelse(use_ot, "top", "bottom")))]
  d[, .(sample, chr, pos, ref, alt, n_alt, depth, strands, ot_alt, ot_depth, ob_alt, ob_depth)]
}))
long[, vaf := fifelse(depth > 0, n_alt / depth, NA_real_)]
long[, called := n_alt >= MIN_ALT & vaf >= MIN_VAF & depth >= MIN_DEPTH]

# coverage on the informative strands (both strands for transversions, one for transitions)
cov <- merge(x[, .(sample, chr, pos, both = OT + OB, single = pmin(OT, OB))], pg, by = c("chr", "pos"))
covg <- cov[, .(median_both = as.double(median(both)), median_one_strand = as.double(median(single)),
                pct_ge10_both = round(100 * mean(both >= 10), 1)), by = .(gene, sample)]
covg <- dcast(covg, gene ~ sample, value.var = c("median_both", "median_one_strand", "pct_ge10_both"))
# panel genes without a single position passing the filters (no reads with MAPQ >= 10) are kept, as zeros
nog <- setdiff(unique(panel$gene), covg$gene)
if (length(nog)) {
  covg <- rbind(covg, data.table(gene = nog), fill = TRUE)
  for (c_ in setdiff(names(covg), "gene")) set(covg, which(is.na(covg[[c_]])), c_, 0)
}
setorder(covg, gene)
fwrite(covg, file.path(OUT, "coverage_by_gene.tsv"), sep = "\t")

## ---- variants called in at least one sample ---------------------------------
keys <- unique(long[called == TRUE, .(chr, pos, ref, alt)])
v <- long[keys, on = .(chr, pos, ref, alt)]
v <- merge(v, SMP[, .(sample, person, order)], by = "sample")
message(nrow(keys), " variants called (>= ", MIN_ALT, " alt reads, VAF >= ", 100 * MIN_VAF, "%, >= ", MIN_DEPTH,
        " reads) in at least one sample")

# artefact flags (per variant)
fl <- copy(keys)
fl[, multiallelic := .N > 1, by = .(chr, pos)]
cpos <- unique(keys[, .(chr, pos)])
fl[, clustered := mapply(function(ch, p) any(cpos$chr == ch & abs(cpos$pos - p) <= 10 & cpos$pos != p), chr, pos)]
win <- as.character(getSeq(G, GRanges(fl$chr, IRanges(fl$pos - 10L, fl$pos + 10L))))
fl[, tandem_repeat := grepl("A{6,}|C{6,}|G{6,}|T{6,}", win) | grepl("(..)\\1{3,}", win, perl = TRUE) |
                      grepl("(...)\\1{3,}", win, perl = TRUE)]
sb <- v[called == TRUE & strands == "both", .(ot_a = sum(ot_alt), ot_d = sum(ot_depth), ob_a = sum(ob_alt), ob_d = sum(ob_depth)),
        by = .(chr, pos, ref, alt)]
sb[, p_strand := mapply(function(a, b, c, d) fisher.test(matrix(c(a, b - a, c, d - c), 2))$p.value, ot_a, ot_d, ob_a, ob_d)]
fl <- merge(fl, sb[, .(chr, pos, ref, alt, p_strand)], by = c("chr", "pos", "ref", "alt"), all.x = TRUE)
# a real variant shows on both original strands; alt reads on one strand only (>= 5 reads on the other) mark
# strand-specific damage or conversion, the commonest artefact here
fl <- merge(fl, sb[, .(chr, pos, ref, alt, ot_a, ot_d, ob_a, ob_d)], by = c("chr", "pos", "ref", "alt"), all.x = TRUE)
fl[, one_strand := !is.na(p_strand) & (p_strand < P_STRAND | (ot_a == 0 & ot_d >= 5) | (ob_a == 0 & ob_d >= 5))]
fl[, c("ot_a", "ot_d", "ob_a", "ob_d") := NULL]

fisher_p <- function(a, d) {
  ok <- d >= MIN_DEPTH
  if (sum(ok) < 2) return(NA_real_)
  m <- rbind(a[ok], d[ok] - a[ok])
  if (sum(m[1, ]) == 0 || sum(m[2, ]) == 0) return(1)
  tryCatch(fisher.test(m, workspace = 2e7)$p.value,
           error = function(e) fisher.test(m, simulate.p.value = TRUE, B = 2e5)$p.value)
}
mn <- function(z) if (length(z)) min(z) else NA_real_
mx <- function(z) if (length(z)) max(z) else NA_real_
per <- v[, .(n_samples = .N, n_covered = sum(depth >= MIN_DEPTH), n_called = sum(called),
             alt_reads = sum(n_alt), reads = sum(depth),
             vaf_min = mn(vaf[depth >= MIN_DEPTH]), vaf_max = mx(vaf[depth >= MIN_DEPTH]),
             p_within = fisher_p(n_alt, depth)), by = .(chr, pos, ref, alt, person)]
per[, vaf_pooled := fifelse(reads > 0, alt_reads / reads, NA_real_)]
per[, status := fcase(
  n_called == 0, "not called",
  n_samples >= 2 & !is.na(p_within) & p_within < P_ACQ, "differs between samples",
  vaf_pooled >= 0.85, "germline, homozygous",
  vaf_pooled >= 0.30, "germline, heterozygous",
  default = "low level in every sample")]
# a variant that differs between one person's samples is an acquired candidate, unless it is a germline SNP in
# someone else: then the difference is allelic imbalance (loss of heterozygosity or a mixed sample)
germ <- per[grepl("^germline", status), .(germline_in = list(person)), by = .(chr, pos, ref, alt)]
low  <- per[n_called > 0 & vaf_pooled < 0.30, .(n_persons_low = .N), by = .(chr, pos, ref, alt)]

## ---- protein change from the MANE Select CDS --------------------------------
gtf <- fread(cmd = sprintf("gzip -dc '%s'", MANE), sep = "\t", header = FALSE, quote = "",
             col.names = c("chr", "src", "type", "start", "end", "score", "strand", "frame", "attr"))
gtf <- gtf[type %in% c("CDS", "stop_codon", "exon") & grepl('tag "MANE_Select"', attr, fixed = TRUE)]
gtf[, `:=`(gene = sub('.*gene_name "([^"]+)".*', "\\1", attr), tx = sub('.*transcript_id "([^"]+)".*', "\\1", attr))]
gtf <- gtf[gene %in% unique(panel$gene)]
tx_info <- rbindlist(lapply(split(gtf, gtf$tx), function(t) {
  s <- t$strand[1]
  cds <- reduce(GRanges(t$chr[1], IRanges(t[type %in% c("CDS", "stop_codon"), start], t[type %in% c("CDS", "stop_codon"), end])))
  cds <- cds[order(start(cds), decreasing = s == "-")]
  gpos <- unlist(lapply(seq_along(cds), function(i) if (s == "+") start(cds)[i]:end(cds)[i] else end(cds)[i]:start(cds)[i]))
  sq <- paste(as.character(getSeq(G, GRanges(t$chr[1], ranges(cds), strand = s))), collapse = "")
  data.table(gene = t$gene[1], tx = t$tx[1], chr = t$chr[1], strand = s, cds_seq = sq,
             gpos = list(gpos), ex_start = list(t[type == "exon", start]), ex_end = list(t[type == "exon", end]))
}))
ok <- tx_info[, .(ok = substr(cds_seq, 1, 3) == "ATG" & nchar(cds_seq) %% 3 == 0 &
                  as.character(translate(DNAString(substr(cds_seq, nchar(cds_seq) - 2, nchar(cds_seq))))) == "*"), by = .(gene, tx)]
message("MANE Select CDS for ", nrow(tx_info), " panel genes; ", sum(ok$ok), " start with ATG, end in a stop and are a multiple of 3",
        if (any(!ok$ok)) paste0(" (not: ", paste(ok[ok == FALSE, gene], collapse = ", "), ")") else "")

comp1 <- c(A = "T", C = "G", G = "C", T = "A")
aa3 <- function(a) ifelse(a == "*", "Ter", AMINO_ACID_CODE[a])
annotate_one <- function(chrom, pos, ref, alt, gname) {
  t <- tx_info[which(tx_info$gene == gname & tx_info$chr == chrom)]
  if (!nrow(t)) return(list(consequence = "no MANE Select transcript", hgvs_c = NA_character_, protein = NA_character_, aa_short = NA_character_))
  t <- t[1]; g <- t$gpos[[1]]; s <- t$strand
  i <- match(pos, g)
  r <- if (s == "+") ref else comp1[[ref]]; a <- if (s == "+") alt else comp1[[alt]]
  if (!is.na(i)) {
    stopifnot(substr(t$cds_seq, i, i) == r)
    k <- (i - 1L) %/% 3L; cod <- substr(t$cds_seq, 3L * k + 1L, 3L * k + 3L)
    alt_cod <- cod; substr(alt_cod, i - 3L * k, i - 3L * k) <- a
    ra <- as.character(translate(DNAString(cod))); aa <- as.character(translate(DNAString(alt_cod)))
    cons <- if (ra == aa) "synonymous" else if (aa == "*") "nonsense" else if (ra == "*") "stop lost" else if (k == 0) "start lost" else "missense"
    return(list(consequence = cons, hgvs_c = sprintf("c.%d%s>%s", i, r, a),
                protein = sprintf("p.%s%d%s", aa3(ra), k + 1L, if (ra == aa) "=" else aa3(aa)),
                aa_short = sprintf("%s%d%s", ra, k + 1L, if (ra == aa) "=" else aa)))
  }
  es <- t$ex_start[[1]]; ee <- t$ex_end[[1]]
  if (any(pos >= es & pos <= ee)) {
    up <- if (s == "+") pos < min(g) else pos > max(g)
    return(list(consequence = if (up) "5' UTR" else "3' UTR", hgvs_c = NA_character_, protein = NA_character_, aa_short = NA_character_))
  }
  d <- min(abs(c(es, ee) - pos))
  list(consequence = if (d <= 2) "splice site (intron +/-1-2)" else sprintf("intron, %d bp from the exon", d),
       hgvs_c = NA_character_, protein = NA_character_, aa_short = NA_character_)
}
# check: DNMT3A R882H (chr2:25,234,373 C>T in hg38) must come out as c.2645G>A p.Arg882His
chk <- annotate_one("chr2", 25234373L, "C", "T", "DNMT3A")
message("annotation check, DNMT3A chr2:25234373 C>T: ", chk$hgvs_c, " ", chk$protein)
stopifnot(identical(chk$hgvs_c, "c.2645G>A"), identical(chk$protein, "p.Arg882His"))

keys <- merge(keys, pg, by = c("chr", "pos"), all.x = TRUE)
ann <- keys[, annotate_one(chr, pos, ref, alt, gene), by = .(chr, pos, ref, alt, gene)]

## ---- hotspot codons: any non-reference read, every sample --------------------
# Codons recurrently mutated in clonal haematopoiesis and myeloid neoplasms. The genomic
# positions come from the MANE CDS, and the reference amino acid is checked.
HOT <- data.table(gene = c("DNMT3A", "IDH1", "IDH2", "IDH2", "JAK2", "SF3B1", "SF3B1", "SRSF2", "U2AF1", "U2AF1",
                           "NRAS", "NRAS", "NRAS", "KRAS", "KRAS", "KRAS", "TP53", "TP53", "TP53", "GNB1",
                           "FLT3", "KIT", "MPL", "SETBP1", "SETBP1"),
                  codon = c(882L, 132L, 140L, 172L, 617L, 666L, 700L, 95L, 34L, 157L, 12L, 13L, 61L, 12L, 13L, 61L,
                            175L, 248L, 273L, 57L, 835L, 816L, 515L, 868L, 870L),
                  aa = c("R", "R", "R", "R", "V", "K", "K", "P", "S", "Q", "G", "G", "Q", "G", "G", "Q",
                         "R", "R", "R", "K", "D", "D", "W", "D", "G"))
hot <- rbindlist(lapply(seq_len(nrow(HOT)), function(i) {
  t <- tx_info[gene == HOT$gene[i]][1]
  if (is.na(t$tx)) return(NULL)
  k <- HOT$codon[i]
  ref_aa <- as.character(translate(DNAString(substr(t$cds_seq, 3L * k - 2L, 3L * k))))
  if (ref_aa != HOT$aa[i]) { message("hotspot check: ", HOT$gene[i], " codon ", k, " is ", ref_aa, ", not ", HOT$aa[i]); return(NULL) }
  data.table(hotspot = paste0(HOT$gene[i], " ", HOT$aa[i], k), chr = t$chr, pos = t$gpos[[1]][(3L * k - 2L):(3L * k)])
}))
# per sample: the codon position and alt allele with the most alt reads (ties: the most reads)
hc <- merge(long, hot, by = c("chr", "pos"))[order(-n_alt, -depth)][, .SD[1], by = .(hotspot, sample)][, .(hotspot, sample, n_alt, depth)]
hc <- merge(hc, SMP[, .(sample, order)], by = "sample")[order(match(hotspot, unique(hot$hotspot)), order)]
hc[, cell := sprintf("%d/%d", n_alt, depth)]
hw <- dcast(hc, hotspot ~ sample, value.var = "cell")
# a codon without any read passing MAPQ >= 10 (U2AF1: falsely duplicated in GRCh38) is reported as such
miss <- setdiff(unique(hot$hotspot), hw$hotspot)
if (length(miss)) hw <- rbind(hw, data.table(hotspot = miss), fill = TRUE)
hw <- hw[match(unique(hot$hotspot), hotspot)]
for (s_ in setdiff(names(hw), "hotspot")) set(hw, which(is.na(hw[[s_]])), s_, "no reads")
setcolorder(hw, c("hotspot", intersect(SMP$sample, names(hw))))
fwrite(hw, file.path(OUT, "hotspots.tsv"), sep = "\t")
message(nrow(hw), " hotspot codons checked; codons with >= 2 non-reference reads in a sample: ",
        if (any(hc$n_alt >= 2)) paste(unique(hc[n_alt >= 2, hotspot]), collapse = ", ") else "none")

## ---- one row per variant and person -----------------------------------------
v[, cell := fifelse(depth > 0, sprintf("%d/%d (%.0f%%)", n_alt, depth, 100 * vaf), "0/0")]
wide <- dcast(v, chr + pos + ref + alt ~ sample, value.var = "cell")
setcolorder(wide, c("chr", "pos", "ref", "alt", intersect(SMP$sample, names(wide))))
strands <- unique(v[, .(chr, pos, ref, alt, strands)])
res <- per[status != "not called"]
res <- merge(res, germ, by = c("chr", "pos", "ref", "alt"), all.x = TRUE)
res[, germline_elsewhere := mapply(function(g, p) !is.null(g) && any(g != p), germline_in, person)]
res[, germline_in := NULL]
res <- merge(res, low, by = c("chr", "pos", "ref", "alt"), all.x = TRUE)
res <- merge(res, fl, by = c("chr", "pos", "ref", "alt"))
res[, recurrent_low := !is.na(n_persons_low) & n_persons_low >= 2]
mxa <- v[, .(max_alt = max(n_alt)), by = .(chr, pos, ref, alt, person)]
res <- merge(res, mxa, by = c("chr", "pos", "ref", "alt", "person"))
res[, thin := n_called == 1 & max_alt <= MIN_ALT]
res[, flags := apply(cbind(fifelse(multiallelic, "multi-allelic", ""), fifelse(clustered, "clustered", ""),
                           fifelse(tandem_repeat, "repeat", ""), fifelse(one_strand, "one strand", ""),
                           fifelse(recurrent_low, "low level in >1 person", ""),
                           fifelse(thin, "one sample, 3 alt reads", "")), 1,
                     function(z) paste(z[nzchar(z)], collapse = "; "))]
res <- merge(res, KNOWN_SNPS, by = c("chr", "pos", "ref", "alt"), all.x = TRUE)
res[, class := fcase(
  status == "differs between samples" & (germline_elsewhere | !is.na(known_snp) | vaf_max >= 0.85),
  "allelic imbalance at a germline SNP",
  status == "differs between samples", "acquired candidate",
  n_samples == 1 & grepl("^germline", status), paste0(status, " (single sample)"),
  n_samples == 1, "low level (single sample)",
  default = status)]
res[, confidence := fifelse(nzchar(flags), "low", "high")]
res <- merge(res, ann, by = c("chr", "pos", "ref", "alt"))
res <- merge(res, strands, by = c("chr", "pos", "ref", "alt"))
res <- merge(res, wide, by = c("chr", "pos", "ref", "alt"))
num <- c("vaf_min", "vaf_max", "vaf_pooled")
res[, (num) := lapply(.SD, function(z) round(100 * z, 1)), .SDcols = num]
res[, `:=`(p_within = signif(p_within, 3), p_strand = signif(p_strand, 3))]
setnames(res, num, paste0(num, "_pct"))
keep <- c("gene", "chr", "pos", "ref", "alt", "hgvs_c", "protein", "aa_short", "consequence", "known_snp", "person", "class",
          "confidence", "flags", "p_within", "p_strand", "strands", "n_samples", "n_covered", "n_called", "vaf_pooled_pct", "vaf_min_pct",
          "vaf_max_pct", intersect(SMP$sample, names(res)))
res <- res[, ..keep]
res <- res[order(person, -(class %in% c("acquired candidate", "allelic imbalance at a germline SNP")), confidence, gene, pos)]
fwrite(res, file.path(OUT, "variants_all.tsv"), sep = "\t")

acq <- res[class %in% c("acquired candidate", "allelic imbalance at a germline SNP")]
fwrite(acq, file.path(OUT, "variants_acquired.tsv"), sep = "\t")
cat("\nclass per person:\n"); print(dcast(res[, .N, by = .(person, class)], class ~ person, value.var = "N", fill = 0))
cat("\nVAF differs between a person's samples:\n")
print(acq[, c("gene", "chr", "pos", "ref", "alt", "aa_short", "consequence", "person", "class", "confidence", "flags", "p_within",
              intersect(SMP$sample, names(acq))), with = FALSE], nrows = 100)

## ---- variants private to one sample ------------------------------------------
# Called in one sample (>= 3 alt reads, VAF >= 10%, >= 8 reads), no alt read in any other sample of the same
# person (each with >= 6 reads), and fewer than 2 alt reads in every sample of the other people. Single-library
# calls at this depth sit near the noise floor; the read support is given so they can be weighed.
long_p <- merge(long, SMP[, .(sample, person)], by = "sample")
pc_ <- long_p[called == TRUE]
pv_ <- merge(pc_[, .(chr, pos, ref, alt, sample, person)], long_p[, .(chr, pos, ref, alt, o_sample = sample, o_person = person, o_alt = n_alt, o_depth = depth)],
             by = c("chr", "pos", "ref", "alt"), allow.cartesian = TRUE)[o_sample != sample]
ok_ <- pv_[, .(same_ok = all(o_alt[o_person == person] == 0) && all(o_depth[o_person == person] >= 6) && any(o_person == person),
               others_ok = all(o_alt[o_person != person] < 2)), by = .(chr, pos, ref, alt, sample)][same_ok & others_ok]
priv <- merge(pc_, ok_[, .(chr, pos, ref, alt, sample)], by = c("chr", "pos", "ref", "alt", "sample"))
priv <- merge(priv, pg, by = c("chr", "pos"), all.x = TRUE)
if (nrow(priv)) priv <- cbind(priv, rbindlist(lapply(seq_len(nrow(priv)), function(i)
  as.data.table(annotate_one(priv$chr[i], priv$pos[i], priv$ref[i], priv$alt[i], priv$gene[i])))))
priv <- merge(priv, fl[, .(chr, pos, ref, alt, tandem_repeat)], by = c("chr", "pos", "ref", "alt"), all.x = TRUE)
priv <- merge(priv, SMP[, .(sample, label, order)], by = "sample")
priv[, both_strands := fifelse(strands == "both", ot_alt > 0 & ob_alt > 0, NA)]
priv <- priv[order(-(both_strands %in% TRUE), -n_alt, order, chr, pos),
             .(sample, label, person, gene, chr, pos, ref, alt, hgvs_c, protein, aa_short, consequence,
               n_alt, depth, vaf_pct = round(100 * vaf, 1), strands, top_alt = ot_alt, top_reads = ot_depth,
               bottom_alt = ob_alt, bottom_reads = ob_depth, both_strands, repeat_context = tandem_repeat)]
fwrite(priv, file.path(OUT, "private_candidates.tsv"), sep = "\t")
cat("\nvariants private to one sample, by read support (alt reads on both strands / all):\n")
print(priv[, .(all = .N, both_strands = sum(both_strands %in% TRUE)), by = sample])
print(priv[both_strands %in% TRUE, .(sample, gene, protein, consequence, n_alt, depth, top_alt, top_reads, bottom_alt, bottom_reads, repeat_context)])

## ---- DNMT3A and TET2: every position with >= 2 alt reads in any sample -------
dt <- merge(long[n_alt >= 2 & vaf >= 0.05], pg[gene %in% c("DNMT3A", "TET2")], by = c("chr", "pos"))
dt <- merge(dt, SMP[, .(sample, person, order)], by = "sample")
dt <- dt[, .(gene, chr, pos, ref, alt, sample, person, n_alt, depth, vaf_pct = round(100 * vaf, 1), strands,
             in_variant_table = paste(chr, pos, ref, alt) %in% res[, paste(chr, pos, ref, alt)])][order(gene, pos, alt, sample)]
fwrite(dt, file.path(OUT, "dnmt3a_tet2_positions_with_alt_reads.tsv"), sep = "\t")

## ---- clone markers: high-confidence protein-changing candidates, every sample --
cp <- function(k, n, lo = TRUE) ifelse(n == 0, NA_real_, if (lo) ifelse(k == 0, 0, qbeta(0.025, k, n - k + 1)) else ifelse(k == n, 1, qbeta(0.975, k + 1, n - k)))
CODING <- c("missense", "nonsense", "start lost", "stop lost", "splice site (intron +/-1-2)")
mk <- unique(acq[confidence == "high" & consequence %in% CODING,
                 .(chr, pos, ref, alt, gene, protein, aa_short, consequence, known_snp, variant_person = person, class)])
pm <- priv[both_strands %in% TRUE & n_alt >= 4 & consequence %in% CODING]
pm <- pm[!mk, on = .(chr, pos, ref, alt)]
if (nrow(pm)) mk <- rbind(mk, pm[, .(chr, pos, ref, alt, gene, protein, aa_short,
                                      consequence, known_snp = NA_character_, variant_person = person,
                                      class = "candidate: one sample, both strands")])
cm <- merge(long[mk[, .(chr, pos, ref, alt)], on = .(chr, pos, ref, alt)], mk, by = c("chr", "pos", "ref", "alt"))
cm <- merge(cm, SMP, by = "sample")
cm[, `:=`(vaf_pct = round(100 * vaf, 1), ci_low_pct = round(100 * cp(n_alt, depth), 1), ci_high_pct = round(100 * cp(n_alt, depth, FALSE), 1))]
cm <- cm[order(gene, pos, order)]
fwrite(cm[, .(gene, chr, pos, ref, alt, protein, consequence, known_snp, class, variant_person, sample, label, sample_person = person,
              n_alt, depth, vaf_pct, ci_low_pct, ci_high_pct, strands, top_alt = ot_alt, top_reads = ot_depth,
              bottom_alt = ob_alt, bottom_reads = ob_depth)], file.path(OUT, "clone_markers.tsv"), sep = "\t")

## ---- figure: VAF of the clone markers in every sample ------------------------
if (nrow(cm)) {
  cm[, lab := factor(axis_label, levels = SMP$axis_label)]
  cm[, case := factor(case, levels = CASE_LV)]
  cm[, var := sprintf("%s %s: %s\n%s:%s %s>%s%s", gene, fifelse(is.na(aa_short), consequence, aa_short),
                      fcase(class == "acquired candidate", "acquired",
                            class == "allelic imbalance at a germline SNP", "germline SNP, allelic imbalance",
                            default = "candidate, one sample"),
                      chr, formatC(pos, big.mark = ",", format = "d"), ref, alt,
                      fifelse(is.na(known_snp), "", paste0(" (", sub(" .*", "", known_snp), ")")))]
  ord <- unique(cm[order(match(class, c("acquired candidate", "candidate: one sample, both strands", "allelic imbalance at a germline SNP")),
                         gene, pos), var])
  cm[, var := factor(var, levels = ord)]
  cm[, txt := sprintf("%d/%d", n_alt, depth)]
  nv <- uniqueN(cm$var)
  p <- ggplot(cm[depth > 0], aes(lab, 100 * vaf, colour = case)) +
    geom_errorbar(aes(ymin = 100 * cp(n_alt, depth), ymax = 100 * cp(n_alt, depth, FALSE)), width = 0.25) +
    geom_point(size = 2.6) +
    geom_text(aes(y = -10, label = txt), size = 3.2, colour = "grey25", show.legend = FALSE) +
    facet_wrap(~ var, ncol = 2) +
    scale_colour_manual(values = CASE_COL, name = NULL) +
    scale_y_continuous(limits = c(-15, 100), breaks = seq(0, 100, 25)) +
    labs(x = NULL, y = "VAF, % (95% interval)",
         title = "Acquired variants and allelic imbalance, read from the TAPS libraries",
         subtitle = paste("alt reads / reads on the genotype-informative strand(s) under each point; samples of other people",
                          "are negative controls; * provisional sample role", sep = "\n")) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(size = 9, lineheight = 0.9), legend.position = "top",
          strip.text = element_text(size = 10.5), panel.grid.minor = element_blank())
  ggsave(file.path(OUT, "vaf_clones.png"), p, width = 14, height = 1.6 + 3.7 * ceiling(nv / 2), dpi = 200, limitsize = FALSE)
}
message("written to ", OUT)
