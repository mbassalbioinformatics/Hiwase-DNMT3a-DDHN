#!/usr/bin/env Rscript
###############################################################################
# 06_dmr_master.R                                                    2026-10-08
#
# One table of every DMR, for three DMR sets per comparison:
#   metilene  q < 0.05 (primary)
#   both      metilene DMRs overlapping a Defiant DMR changing the same way
#   defiant   Defiant c0/d1 output (c0_CpN5_d1_G147_p0.05_P10)
# Direction is always g1 - g2 (first-named minus second). Defiant reports
# second minus first, so its difference is recomputed from the two means.
#
# For each DMR:
#   - methylation of every sample (mean of CpGs with >= 5 reads, TAPS reading)
#   - gene context from MANE Select v1.0. The switch region is promoter
#     (TSS -1500 to -1) + exon 1 + intron 1, where Mesa regulatory switches sit;
#     beyond it: other exons or introns of a gene, or intergenic. Nearest TSS and
#     signed distance (upstream negative)
#   - CpG-island context and ENCODE cCRE class
#   - the gene's expression change in the matching RNA-seq contrast (DNMT3A
#     comparisons with RNA-seq only), for the switch-region gene, else the
#     nearest gene within 50 kb
#
# Output in integration/dmr/:
#   dmr_master.tsv.gz   every DMR x set x comparison
#   dmr_counts.tsv      per comparison and set: DMRs gained / lost, by context
#   bed/<set>/<comparison>.bed  for HOMER annotatePeaks on the server (01_methylation_server/05_homer_annotation.sh)
###############################################################################
source("config/paths.R")
source("04_integration/io_helpers.R")
suppressPackageStartupMessages(library(data.table))
setDTthreads(8)
DMRD <- file.path(EMS, "03_dmr_analysis")
COMP <- COMPARISONS
OUT  <- file.path(INT, "dmr"); dir.create(OUT, showWarnings = FALSE)
source("04_integration/locus_plot_functions.R")       # LOCUS_SAMPLES
source("03_rnaseq/00_config.R")                  # CONTRASTS
comp <- fread(COMP, colClasses = "character")
cache <- Sys.getenv("LOCUS_CACHE")
L <- if (nzchar(cache) && file.exists(cache)) readRDS(cache) else load_locus_data(
  emseq = dirname(DMRD), mane_gtf = MANE, cgi_file = file.path(RES, "cpgIslandExt.txt.gz"),
  rna_rdata = file.path(RNA, "05_processed_counts_genes", "dge_gene_level.RData"))

## ---- DMR sets -------------------------------------------------------------------
read_metilene <- function(id) {
  x <- fread_metilene(DMRD, id, header = FALSE, col.names = c("chr", "start", "end", "q", "delta", "n_cpg", "p_mwu", "p_2dks", "mean_g1", "mean_g2"))
  x[q < 0.05]
}
read_defiant <- function(id, g1, g2) {
  f <- list.files(file.path(DMRD, "defiant", id), "_c0_CpN5_d1_G147_p0\\.05_P10\\.tsv$", full.names = TRUE)
  x <- fread(f, colClasses = list(character = 1))
  setnames(x, 1:6, c("chr", "start", "end", "n_cpg", "n_diff_cpg", "defiant_diff"))
  x[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
  num <- function(v) as.numeric(gsub("[^0-9.eE+-]", "", v))
  x[, `:=`(mean_g1 = num(get(g1)) / 100, mean_g2 = num(get(g2)) / 100)]
  x[, delta := mean_g1 - mean_g2]
  x[, .(chr, start, end, n_cpg, delta, mean_g1, mean_g2)]
}
ov <- function(a, b) {                                           # rows of a overlapping any row of b
  b <- copy(b)[, .(chr, start, end)]; setkey(b, chr, start, end)
  unique(foverlaps(a[, .(i = .I, chr, start, end)], b, by.x = c("chr", "start", "end"), nomatch = 0L)$i)
}

sets <- rbindlist(lapply(seq_len(nrow(comp)), function(k) {
  id <- comp$id[k]; g1 <- comp$g1[k]; g2 <- comp$g2[k]
  m <- read_metilene(id); d <- read_defiant(id, g1, g2)
  both <- m[unique(c(ov(m[delta > 0], d[delta > 0]) |> (\(i) which(m$delta > 0)[i])(),
                     ov(m[delta < 0], d[delta < 0]) |> (\(i) which(m$delta < 0)[i])()))]
  message(sprintf("%-22s metilene %4d  both %4d  defiant %6d", id, nrow(m), nrow(both), nrow(d)))
  rbind(m[, .(set = "metilene", chr, start, end, delta, q, n_cpg)],
        both[, .(set = "both", chr, start, end, delta, q, n_cpg)],
        d[, .(set = "defiant", chr, start, end, delta, q = NA_real_, n_cpg)])[, `:=`(id = id, g1 = g1, g2 = g2)]
}))
sets[, dmr := sprintf("%s:%d-%d", chr, start, end)]
sets[, direction := fifelse(delta > 0, "gained", "lost")]

## ---- methylation of every sample in every DMR (metilene and both sets) ---------------
regions <- unique(sets[set != "defiant", .(dmr, chr, start, end)])
for (s in LOCUS_SAMPLES$sample) {
  x <- L$cpg[[s]][t >= 5]
  r <- copy(regions)
  # mean of CpG fractions inside each region, via prefix sums per chromosome
  v <- rbindlist(lapply(unique(r$chr), function(ch) {
    xc <- x[chr == ch]; rc <- r[chr == ch]
    pre_f <- c(0, cumsum(xc$m / xc$t)); pre_n <- c(0, seq_len(nrow(xc)))
    i1 <- findInterval(rc$start - 1, xc$pos) + 1L; i2 <- findInterval(rc$end, xc$pos) + 1L
    n <- pre_n[i2] - pre_n[i1]
    data.table(dmr = rc$dmr, v = fifelse(n >= 1, (pre_f[i2] - pre_f[i1]) / n, NA_real_))
  }))
  regions <- merge(regions, setnames(v, "v", paste0("meth_", s)), by = "dmr", all.x = TRUE)
}

## ---- gene context ---------------------------------------------------------------
g <- fread(cmd = sprintf("gzip -dc '%s'", MANE), header = FALSE, sep = "\t", quote = "")
g <- g[grepl('tag "MANE_Select"', V9, fixed = TRUE) & V3 %in% c("transcript", "exon")]
g[, gene := sub('.*gene_name "([^"]+)".*', "\\1", V9)]
g[V3 == "exon", exon_number := as.integer(sub('.*exon_number ([0-9]+);.*', "\\1", V9))]
tx <- g[V3 == "transcript", .(gene, chr = V1, tx_start = V4, tx_end = V5, strand = V7)][!duplicated(gene)]
tx[, tss := fifelse(strand == "+", tx_start, tx_end)]
ex <- g[V3 == "exon", .(gene, chr = V1, start = V4, end = V5, strand = V7, exon_number)]
e1 <- ex[exon_number == 1][!duplicated(gene)]; e2 <- ex[exon_number == 2][!duplicated(gene)]
seg <- rbind(
  tx[, .(gene, chr, part = "promoter", start = fifelse(strand == "+", tss - 1500L, tss + 1L),
         end = fifelse(strand == "+", tss - 1L, tss + 1500L))],
  e1[, .(gene, chr, part = "exon 1", start, end)],
  merge(e1, e2, by = c("gene", "chr", "strand"), suffixes = c("1", "2"))[
    , .(gene, chr, part = "intron 1", start = fifelse(strand == "+", end1 + 1L, end2 + 1L),
        end = fifelse(strand == "+", start2 - 1L, start1 - 1L))][end >= start])
body <- tx[, .(gene, chr, start = tx_start, end = tx_end)]
setkey(seg, chr, start, end); setkey(body, chr, start, end)

r <- regions[, .(dmr, chr, start, end)]
hit <- foverlaps(r, seg, by.x = c("chr", "start", "end"), nomatch = 0L)
sw <- hit[, .(switch_gene = paste(unique(gene), collapse = ","),
              switch_part = paste(unique(part[order(match(part, c("promoter", "exon 1", "intron 1")))]), collapse = "+")), by = dmr]
bh <- foverlaps(r, body, by.x = c("chr", "start", "end"), nomatch = 0L)[, .(body_gene = paste(unique(gene), collapse = ",")), by = dmr]
r <- merge(merge(r, sw, by = "dmr", all.x = TRUE), bh, by = "dmr", all.x = TRUE)
# nearest TSS, signed distance (negative = upstream of the gene)
r[, mid := (start + end) %/% 2L]
nn <- rbindlist(lapply(unique(r$chr), function(ch) {
  t <- tx[chr == ch][order(tss)]; rr <- r[chr == ch]
  i <- findInterval(rr$mid, t$tss); i <- pmax(1L, i)
  j <- pmin(nrow(t), i + 1L)
  pick <- ifelse(abs(rr$mid - t$tss[i]) <= abs(rr$mid - t$tss[j]), i, j)
  d <- rr$mid - t$tss[pick]
  data.table(dmr = rr$dmr, nearest_gene = t$gene[pick], tss_distance = ifelse(t$strand[pick] == "+", d, -d))
}))
r <- merge(r, nn, by = "dmr")
r[, context := fifelse(!is.na(switch_gene), "switch region (promoter-exon1-intron1)",
               fifelse(!is.na(body_gene), "gene body (beyond intron 1)", "intergenic"))]

## CpG island context and ENCODE cCRE class
cgi <- copy(L$cgi); setkey(cgi, chr, start, end)
fl <- function(w0, w1) rbind(cgi[, .(chr, start = pmax(1L, start - w1), end = start - w0 - 1L)],
                             cgi[, .(chr, start = end + w0 + 1L, end = end + w1)])[end >= start]
on <- function(b) { b <- copy(b); setkey(b, chr, start, end); r$dmr %in% unique(foverlaps(r[, .(dmr, chr, start, end)], b, by.x = c("chr", "start", "end"), nomatch = 0L)$dmr) }
r[, cgi_context := fifelse(on(cgi), "CGI", fifelse(on(fl(0L, 2000L)), "shore", fifelse(on(fl(2000L, 4000L)), "shelf", "open sea")))]
cc <- fread(file.path(RES, "encode_cCREs_hg38_ucsc.bed"), select = c(1, 2, 3, 11), col.names = c("chr", "start", "end", "label"))
cc[, start := start + 1L]
r[, ccre := "none"]
for (lab in rev(c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"))) r[on(cc[label == lab, .(chr, start, end)]), ccre := lab]
r[, mid := NULL]

## ---- assemble -------------------------------------------------------------------
master <- merge(sets[set != "defiant"], merge(r, regions[, !c("chr", "start", "end")], by = "dmr"),
                by = c("dmr", "chr", "start", "end"))
master[, link_gene := fifelse(!is.na(switch_gene), sub(",.*", "", switch_gene),
                      fifelse(abs(tss_distance) <= 50000, nearest_gene, NA_character_))]
master[, link_type := fifelse(!is.na(switch_gene), "switch region", fifelse(!is.na(link_gene), "nearest TSS <= 50 kb", NA_character_))]

# expression of the linked gene in the RNA-seq contrast with the same two samples
e <- new.env(); load(file.path(RNA, "05_processed_counts_genes", "dge_gene_level.RData"), envir = e)
de <- as.data.table(e$de)
rmap <- CONTRASTS[, .(rna_contrast = id, g1 = num, g2 = den)]
master <- merge(master, rmap, by = c("g1", "g2"), all.x = TRUE)
master <- merge(master, de[, .(rna_contrast = contrast, link_gene = gene, rna_lfc = lfc_shrunk, rna_fdr = FDR,
                               rna_sig = sig, rna_on_off = on_off, rna_ruler_flagged = ruler_flagged)],
                by = c("rna_contrast", "link_gene"), all.x = TRUE)
master[, switch_concordant := fifelse(is.na(rna_lfc), NA, sign(rna_lfc) == -sign(delta))]
master <- merge(master, comp[, .(id, status, label, case)], by = "id")
setcolorder(master, c("set", "id", "label", "status", "case", "g1", "g2", "dmr", "chr", "start", "end", "direction",
                      "delta", "q", "n_cpg", "context", "switch_gene", "switch_part", "body_gene", "nearest_gene",
                      "tss_distance", "cgi_context", "ccre", "link_gene", "link_type", "rna_contrast", "rna_lfc",
                      "rna_fdr", "rna_sig", "rna_on_off", "rna_ruler_flagged", "switch_concordant"))
setorder(master, set, id, q)
fwrite(master, file.path(OUT, "dmr_master.tsv.gz"), sep = "\t", na = "NA")
fwrite(sets[set == "defiant", .(id, chr, start, end, delta, n_cpg)], file.path(OUT, "defiant_dmrs.tsv.gz"), sep = "\t")

## ---- counts ---------------------------------------------------------------------
cnt <- rbind(sets[, .(context = "all", n = .N, gained = sum(delta > 0), lost = sum(delta < 0)), by = .(set, id)],
             master[, .(n = .N, gained = sum(delta > 0), lost = sum(delta < 0)), by = .(set, id, context)])
cnt <- merge(cnt, comp[, .(id, label, status)], by = "id")
fwrite(cnt[order(set, id, context)], file.path(OUT, "dmr_counts.tsv"), sep = "\t")
print(dcast(cnt[context == "all"], id ~ set, value.var = "n")[match(comp$id, id)])
message("\nmetilene DMRs by context:")
print(dcast(master[set == "metilene"], id ~ context, fun.aggregate = length, value.var = "dmr")[match(comp$id, id)])

## ---- BEDs for HOMER annotatePeaks.pl on the server ----------------------------------
for (st in c("metilene", "both")) {
  dir.create(file.path(OUT, "bed", st), recursive = TRUE, showWarnings = FALSE)
  for (i in comp$id) {
    b <- master[set == st & id == i, .(chr, start, end, name = paste0(i, "_", direction, "_", seq_len(.N)), score = round(delta, 3), strand = "+")]
    fwrite(b, file.path(OUT, "bed", st, paste0(i, ".bed")), sep = "\t", col.names = FALSE)
  }
}
message("written to ", OUT)
