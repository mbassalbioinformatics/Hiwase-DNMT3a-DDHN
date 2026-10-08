#!/usr/bin/env Rscript
###############################################################################
# 09_hic_loops.R                                                     2026-10-08
#
# Which DMRs sit in a chromatin loop anchor, what is at the other end, and is
# that partner region also differentially methylated (significantly or not)?
#
# Hi-C loops: John Pimanda's sorted human HSC, CMP, GMP and MEP (two replicates
# each; Projects/mss/encode_human/13_hic/*_Pimanda_*.bedpe.xz). The other files
# in that folder are not used.
# The call sets are permissive (HSC rep1: 990k loops, most > 1 Mb with < 10
# contacts) and the two replicates of a cell type share only 1-3% of loops, so
# requiring both replicates leaves almost nothing (23 of 641 DMRs). Instead:
#   - loops spanning 10 kb to 1 Mb (enhancer-promoter range), from all 8 files;
#   - for each DMR, the partner anchors (the other end of every loop with an
#     anchor on the DMR) are merged when within 10 kb of each other;
#   - a partner is kept when at least MIN_SUPPORT of the 8 files (cell type x
#     replicate) contain it; its support and cell types are reported.
#
# For each DMR-loop-partner (DMR sets: metilene, primary; all comparisons):
#   partner genes          MANE Select TSS within the partner anchor (+/- 2 kb)
#   partner DMR            a metilene DMR (q < 0.05) of the same comparison
#                          inside the partner anchor
#   partner best region    the most significant region metilene tested inside the
#                          partner anchor for that comparison, with its change and
#                          q, significant or not
#   partner gene RNA       log2FC in the matching RNA-seq contrast
#
# Output in integration/hic/:
#   dmr_loops.tsv          one row per DMR x cell type x partner
#   dmr_loops_summary.tsv  per comparison: DMRs looped, partners with genes / DMRs
###############################################################################
source("config/paths.R")
source("04_integration/io_helpers.R")
suppressPackageStartupMessages(library(data.table))
setDTthreads(8)
DMRD <- file.path(EMS, "03_dmr_analysis")
HIC  <- HIC_DIR
COMP <- COMPARISONS
OUT  <- file.path(INT, "hic"); dir.create(OUT, showWarnings = FALSE)
MIN_D <- 10000L; MAX_D <- 1000000L
KEEP <- paste0("chr", c(1:22, "X"))
comp <- fread(COMP, colClasses = "character")
source("03_rnaseq/00_config.R")

## ---- loops from all 8 files --------------------------------------------------------
MIN_SUPPORT <- as.integer(Sys.getenv("MIN_SUPPORT", "2"))
read_loops <- function(f) {
  x <- fread(cmd = sprintf("xz -dc '%s'", f), header = FALSE, skip = 2, select = c(1:6, 12),
             col.names = c("c1", "s1", "e1", "c2", "s2", "e2", "observed"))
  x[c1 %in% KEEP & c1 == c2 & (s2 - s1) >= MIN_D & (s2 - s1) <= MAX_D]
}
files <- list.files(HIC, "^(HSC|CMP|GMP|MEP)_rep[12]_Pimanda_.*\\.bedpe\\.xz$", full.names = TRUE)
loops <- rbindlist(lapply(files, function(f) {
  x <- read_loops(f); lab <- sub("_Pimanda.*", "", basename(f))
  message(sprintf("%-9s %7d loops in 10 kb-1 Mb", lab, nrow(x)))
  x[, `:=`(file = lab, cell_type = sub("_rep.*", "", lab))]
}))
loops[, loop := sprintf("%s:%d-%d|%d-%d|%s", c1, s1, e1, s2, e2, file)]

## ---- DMRs in anchors -----------------------------------------------------------------
m <- fread(cmd = sprintf("gzip -dc '%s'", file.path(INT, "dmr", "dmr_master.tsv.gz")))[set == "metilene"]
d <- m[, .(id, dmr, chr, start, end, delta, q, dmr_gene = fifelse(!is.na(switch_gene), switch_gene, nearest_gene), context)]
anch <- rbind(loops[, .(loop, file, cell_type, chr = c1, a_start = s1, a_end = e1, p_start = s2, p_end = e2, observed, side = 1L)],
              loops[, .(loop, file, cell_type, chr = c1, a_start = s2, a_end = e2, p_start = s1, p_end = e1, observed, side = 2L)])
setkey(anch, chr, a_start, a_end)
hit <- foverlaps(d, anch, by.x = c("chr", "start", "end"), nomatch = 0L)
message(nrow(hit), " DMR-loop hits (all files)")
# merge a DMR's partner anchors within 10 kb of each other; support = files containing the partner
setorder(hit, id, dmr, p_start)
hit[, grp := cumsum(c(TRUE, p_start[-1] - cummax(p_end)[-.N] > 10000L)), by = .(id, dmr)]
hit <- hit[, .(chr = chr[1], start = start[1], end = end[1], delta = delta[1], q = q[1], dmr_gene = dmr_gene[1],
               context = context[1], p_start = min(p_start), p_end = max(p_end), support = uniqueN(file),
               files = paste(sort(unique(file)), collapse = ","), cell_type = paste(sort(unique(cell_type)), collapse = ","),
               observed = max(observed), loop = paste(unique(loop), collapse = ";"), a_start = min(a_start)),
           by = .(id, dmr, grp)][support >= MIN_SUPPORT]
message("partners supported by >= ", MIN_SUPPORT, " files: ", nrow(hit), "; DMRs looped: ",
        uniqueN(hit[, paste(id, dmr)]), " of ", nrow(d))

## partner genes (TSS within the partner anchor +/- 2 kb)
tx <- fread(cmd = sprintf("gzip -dc '%s'", file.path(INT, "region_methylation", "mane_select_genes.tsv.gz")))
tss <- tx[, .(gene, chr, s = tss - 2000L, e = tss + 2000L)]; setkey(tss, chr, s, e)
pg <- foverlaps(unique(hit[, .(chr, p_start, p_end)]), tss, by.x = c("chr", "p_start", "p_end"), nomatch = 0L)[
  , .(partner_genes = paste(unique(gene), collapse = ",")), by = .(chr, p_start, p_end)]
hit <- merge(hit, pg, by = c("chr", "p_start", "p_end"), all.x = TRUE)

## partner methylation: significant DMR and best-tested region, per comparison
part <- rbindlist(lapply(unique(hit$id), function(i) {
  x <- fread_metilene(DMRD, i, header = FALSE, select = 1:6,
             col.names = c("chr", "start", "end", "q", "delta", "n_cpg"))
  p <- unique(hit[id == i, .(id, chr, p_start, p_end)])
  setkey(x, chr, start, end)
  o <- foverlaps(p, x, by.x = c("chr", "p_start", "p_end"), nomatch = 0L)
  o[, .(partner_n_tested = .N, partner_n_dmr = sum(q < 0.05),
        partner_best_delta = delta[which.min(q)], partner_best_q = min(q),
        partner_dmr_deltas = paste(sprintf("%+.2f", delta[q < 0.05]), collapse = ",")), by = .(id, chr, p_start, p_end)]
}))
hit <- merge(hit, part, by = c("id", "chr", "p_start", "p_end"), all.x = TRUE)

## partner gene expression in the matching RNA-seq contrast
e <- new.env(); load(file.path(RNA, "05_processed_counts_genes", "dge_gene_level.RData"), envir = e)
de <- as.data.table(e$de)
rmap <- merge(comp[, .(id, g1, g2)], CONTRASTS[, .(rna_contrast = id, g1 = num, g2 = den)], by = c("g1", "g2"))
hit <- merge(hit, rmap[, .(id, rna_contrast)], by = "id", all.x = TRUE)
expr_of <- function(genes, rc) {
  if (is.na(genes) || is.na(rc)) return(NA_character_)
  g <- strsplit(genes, ",")[[1]]
  x <- de[contrast == rc & gene %in% g]
  if (!nrow(x)) return(NA_character_)
  paste(sprintf("%s %+.2f%s", x$gene, x$lfc_shrunk, ifelse(x$sig, "*", "")), collapse = "; ")
}
hit[, partner_gene_rna := mapply(expr_of, partner_genes, rna_contrast)]
hit[, `:=`(loop_distance = abs(p_start - a_start))]
out <- hit[, .(id, dmr, dmr_gene, context, delta, q, support, cell_type, files, loop_distance, observed,
               partner = sprintf("%s:%d-%d", chr, p_start, p_end), partner_genes, partner_gene_rna,
               partner_n_tested, partner_n_dmr, partner_dmr_deltas, partner_best_delta, partner_best_q)]
setorder(out, id, q, -support)
fwrite(out, file.path(OUT, "dmr_loops.tsv"), sep = "\t", na = "NA")

summ <- out[, .(dmrs_looped = uniqueN(dmr), loops = .N, cell_types = paste(sort(unique(unlist(strsplit(cell_type, ",")))), collapse = ","),
                partners_with_genes = sum(!is.na(partner_genes)), partners_with_dmr = sum(partner_n_dmr > 0, na.rm = TRUE)), by = id]
summ <- merge(comp[, .(id, label, dmrs = sapply(id, function(i) nrow(d[id == i])))], summ, by = "id", all.x = TRUE)
fwrite(summ, file.path(OUT, "dmr_loops_summary.tsv"), sep = "\t", na = "NA")
print(summ[match(comp$id, id), .(label, dmrs, dmrs_looped, loops, partners_with_genes, partners_with_dmr)])
message("written to ", OUT)
