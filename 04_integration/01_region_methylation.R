#!/usr/bin/env Rscript
###############################################################################
# 01_region_methylation.R                                            2026-10-08
#
# Methylation per gene promoter, per gene body and per 1-kb genomic tile, for
# every sample, from Msuite2's per-CpG counts (TAPS: methylated = wT + cT,
# total = wC + wT + cC + cT). This does not depend on DMR
# calling: every gene gets a value, so subtle methylation changes can be set
# against expression changes.
#
# Regions, from MANE Select v1.5 (one transcript per gene, hg38):
#   promoter  TSS -1500 to +500, strand-aware
#   body      TSS +500 to the transcript end
#   tiles     1 kb, chr1-22 and chrX
# For each region and sample: CpGs covered, CpGs with >= 5 reads, methylated and
# total reads, and the mean methylation of the CpGs with >= 5 reads. That mean is
# NA when fewer than 3 such CpGs fall in the region.
#
# Output in integration/region_methylation/:
#   promoter_methylation.tsv.gz   gene x sample, fraction 0-1
#   promoter_counts.tsv.gz        the counts behind it, long format
#   body_methylation.tsv.gz       gene x sample
#   tiles_1kb_methylation.tsv.gz  tile x sample
###############################################################################
suppressPackageStartupMessages(library(data.table))
setDTthreads(8)

source("config/paths.R")
EMSEQ <- EMS
OUT   <- file.path(INT, "region_methylation")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# every library analysed, in sample-sheet order (timeline order within each case)
SAMPLES <- fread(SAMPLE_SHEET, colClasses = "character")[methylation %in% c("analysed", "extra"), sample]
KEEP    <- paste0("chr", c(1:22, "X"))
MIN_COV <- 5     # reads for a CpG to count towards the mean
MIN_CPG <- 3     # CpGs at MIN_COV for a region to get a value

## ---- regions ------------------------------------------------------------------
g <- fread(cmd = sprintf("gzip -dc '%s'", MANE), header = FALSE, sep = "\t", quote = "")
g <- g[V3 == "transcript" & grepl('tag "MANE_Select"', V9, fixed = TRUE) & V1 %in% KEEP]
attr_of <- function(x, key) sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x)
tx <- g[, .(chr = V1, tx_start = V4, tx_end = V5, strand = V7,
            gene = attr_of(V9, "gene_name"),
            ensembl = sub("\\..*", "", attr_of(V9, "gene_id")),
            transcript = attr_of(V9, "transcript_id"),
            gene_type = attr_of(V9, "gene_type"))]
tx[, tss := fifelse(strand == "+", tx_start, tx_end)]
tx <- tx[!duplicated(gene)]
prom <- tx[, .(id = gene, chr,
               start = fifelse(strand == "+", tss - 1500L, tss - 500L),
               end   = fifelse(strand == "+", tss + 500L,  tss + 1500L))]
body <- tx[, .(id = gene, chr,
               start = fifelse(strand == "+", tss + 501L, tx_start),
               end   = fifelse(strand == "+", tx_end,     tss - 501L))][end > start]
fwrite(tx, file.path(OUT, "mane_select_genes.tsv.gz"), sep = "\t")
message(nrow(tx), " MANE Select genes; ", nrow(body), " with a body beyond the promoter")

## ---- per-region sums from prefix sums ------------------------------------------
# For CpGs at positions p (sorted) with methylated m and total t reads, the sum
# over [start, end] is a difference of prefix sums at findInterval() positions.
region_sums <- function(x, reg) {
  rbindlist(lapply(intersect(KEEP, unique(reg$chr)), function(ch) {
    xc <- x[chr == ch]; r <- reg[chr == ch]
    if (!nrow(xc)) return(data.table(id = r$id, n = 0L, n5 = 0L, m = 0, t = 0, f5 = 0))
    ok5 <- xc$t >= MIN_COV
    pre <- function(v) c(0, cumsum(as.numeric(v)))
    cn <- pre(rep(1, nrow(xc))); cn5 <- pre(ok5); cm <- pre(xc$m); ct <- pre(xc$t)
    cf5 <- pre(ifelse(ok5, xc$m / xc$t, 0))
    i1 <- findInterval(r$start - 1, xc$pos) + 1L
    i2 <- findInterval(r$end, xc$pos) + 1L
    data.table(id = r$id, n = cn[i2] - cn[i1], n5 = cn5[i2] - cn5[i1],
               m = cm[i2] - cm[i1], t = ct[i2] - ct[i1], f5 = cf5[i2] - cf5[i1])
  }))
}

read_cpgs <- function(s) {
  f <- file.path(EMSEQ, paste0(s, "_msuite2"), "Msuite2.CpG.meth.call.xz")
  x <- fread(cmd = sprintf("xz -dc '%s'", f), select = c(1, 2, 4, 5, 8, 9),
             col.names = c("chr", "pos", "wC", "wT", "cC", "cT"))
  x <- x[chr %in% KEEP]
  x[, `:=`(m = wT + cT, t = wC + wT + cC + cT)]       # TAPS: T = methylated
  x <- x[t > 0, .(chr, pos, m, t)]
  setkey(x, chr, pos)
  x
}

prom_long <- list(); body_mat <- list(); tile_mat <- list(); prom_mat <- list()
for (s in SAMPLES) {
  t0 <- Sys.time()
  x <- read_cpgs(s)
  if (!length(tile_mat)) {                            # tiles from the first sample's extent
    ext <- x[, .(maxpos = max(pos)), by = chr]
    tiles <- ext[, .(start = seq(1L, maxpos, by = 1000L)), by = chr]
    tiles[, `:=`(end = start + 999L, id = paste0(chr, ":", start))]
  }
  p <- region_sums(x, prom); p[, sample := s]
  prom_long[[s]] <- p
  prom_mat[[s]] <- p[, .(id, v = fifelse(n5 >= MIN_CPG, f5 / n5, NA_real_))]
  b <- region_sums(x, body)
  body_mat[[s]] <- b[, .(id, v = fifelse(n5 >= MIN_CPG, f5 / n5, NA_real_))]
  tl <- region_sums(x, tiles)
  tile_mat[[s]] <- tl[, .(id, v = fifelse(n5 >= MIN_CPG, f5 / n5, NA_real_))]
  message(sprintf("%-10s %5.1f M CpGs  promoters with a value %d  mean promoter %.3f  mean tile %.3f  (%.0f s)",
                  s, nrow(x) / 1e6, sum(!is.na(prom_mat[[s]]$v)), mean(prom_mat[[s]]$v, na.rm = TRUE),
                  mean(tile_mat[[s]]$v, na.rm = TRUE), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  rm(x); gc(verbose = FALSE)
}

wide <- function(lst) {
  d <- Reduce(function(a, b) merge(a, b, by = "id", all = TRUE),
              lapply(names(lst), function(s) setnames(copy(lst[[s]]), "v", s)))
  d[, c("id", SAMPLES), with = FALSE]
}
fwrite(wide(prom_mat), file.path(OUT, "promoter_methylation.tsv.gz"), sep = "\t", na = "NA")
fwrite(rbindlist(prom_long), file.path(OUT, "promoter_counts.tsv.gz"), sep = "\t")
fwrite(wide(body_mat), file.path(OUT, "body_methylation.tsv.gz"), sep = "\t", na = "NA")
tw <- wide(tile_mat)
tw <- tw[rowSums(!is.na(tw[, ..SAMPLES])) > 0]
fwrite(tw, file.path(OUT, "tiles_1kb_methylation.tsv.gz"), sep = "\t", na = "NA")
message("written to ", OUT)
