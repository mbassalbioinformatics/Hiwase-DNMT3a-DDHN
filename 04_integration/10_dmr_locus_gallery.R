#!/usr/bin/env Rscript
###############################################################################
# 10_dmr_locus_gallery.R                                             2026-10-08
#
# DMR-centred locus plots (Mesa-paper style): the DMR boxed where the Mesa
# figures box MM, its flanks, every gene in the window, CpG-by-CpG methylation
# in every sample, the linked gene's RNA-seq beside it, and, when the DMR sits in
# a Hi-C loop anchor, the loops to its partner regions with their methylation
# change below.
#
# One plot per DMR region (metilene, q < 0.05), titled with every comparison
# that calls it. Chosen, in this order:
#   switch_region/  DMRs in a gene's promoter-exon1-intron1 (where Mesa
#                   regulatory switches sit)
#   panel_genes/    DMRs linked to a panel gene (07_gene_panels.R) by switch
#                   region, nearest TSS within 50 kb, or a Hi-C loop partner
#   top/            the 10 most significant DMRs of each comparison not
#                   already plotted
# Window: DMR +/- 1.5 kb, widened to reach the gene's TSS when it is within
# 5 kb, then narrowed towards the DMR and that TSS (never cutting either) to at
# most 90 CpGs.
# The folder is cleared and every plot redrawn on each run; index.tsv lists them.
###############################################################################
suppressPackageStartupMessages(library(data.table))
source("config/paths.R")
source("04_integration/locus_plot_functions.R")
OUT <- file.path(INT, "locus_plots_dmr")
unlink(OUT, recursive = TRUE); dir.create(OUT)
cache <- Sys.getenv("LOCUS_CACHE")
L <- if (nzchar(cache) && file.exists(cache)) readRDS(cache) else load_locus_data(
  emseq = EMS,
  mane_gtf = MANE,
  cgi_file = file.path(RES, "cpgIslandExt.txt.gz"),
  rna_rdata = file.path(RNA, "05_processed_counts_genes", "dge_gene_level.RData"))
L$min_cov <- as.integer(Sys.getenv("MIN_COV", "2"))   # a CpG is drawn when it has at least this many reads

m <- fread(cmd = sprintf("gzip -dc '%s'", file.path(INT, "dmr", "dmr_master.tsv.gz")))[set == "metilene"]
lp <- fread(file.path(INT, "hic", "dmr_loops.tsv"))
pan <- fread(file.path(RES, "gene_panels.tsv"))
tx <- fread(cmd = sprintf("gzip -dc '%s'", file.path(INT, "region_methylation", "mane_select_genes.tsv.gz")))

## one row per DMR region
calls <- m[, .(calls = paste(sprintf("%s (%+.2f, q %.1e)", label, delta, q), collapse = "; "),
               best_q = min(q)), by = dmr]
u <- unique(m[order(q)], by = "dmr")[, .(dmr, chr, start, end, context, switch_gene, switch_part, body_gene, nearest_gene,
                                        tss_distance, link_gene)]
u <- merge(u, calls, by = "dmr")
loop_genes <- lp[!is.na(partner_genes), .(loop_genes = paste(unique(unlist(strsplit(partner_genes, ","))), collapse = ",")), by = dmr]
u <- merge(u, loop_genes, by = "dmr", all.x = TRUE)
genes_of <- function(i) unique(na.omit(c(strsplit(u$switch_gene[i], ",")[[1]],
                                         if (abs(u$tss_distance[i]) <= 50000) u$nearest_gene[i],
                                         strsplit(u$loop_genes[i], ",")[[1]])))
u[, panels := sapply(seq_len(.N), function(i) { g <- genes_of(i); p <- unique(pan[gene %in% g, panel]); if (length(p)) paste(p, collapse = ",") else NA_character_ })]
u[, category := fifelse(!is.na(switch_gene), "switch_region", fifelse(!is.na(panels), "panel_genes", NA_character_))]
top <- m[order(q), head(.SD, 10), by = id][, unique(dmr)]
u[is.na(category) & dmr %in% top, category := "top"]
u <- u[!is.na(category)]
message(nrow(u), " DMR regions to plot: ", paste(names(table(u$category)), table(u$category), collapse = ", "))

## window: DMR +/- 1.5 kb; the gene's TSS kept in view when it is within 5 kb.
# Narrowed towards the core (the DMR, and the TSS when near) until at most 90 CpGs
# remain; the core itself is never cut.
cpos <- L$cpg[[1]][t >= 1, .(chr, pos)]
window_for <- function(r) {
  core_a <- r$start - 150L; core_b <- r$end + 150L
  g <- sub(",.*", "", r$switch_gene)
  if (!is.na(g)) {
    t0 <- tx[gene == g, tss]
    if (length(t0) && abs(t0 - r$start) <= 5000) { core_a <- min(core_a, t0 - 300L); core_b <- max(core_b, t0 + 300L) }
  }
  a <- min(core_a, r$start - 1500L); b <- max(core_b, r$end + 1500L)
  repeat {
    n <- cpos[chr == r$chr & pos >= a & pos <= b, .N]
    if (n <= 90 || (core_a - a <= 50 && b - core_b <= 50)) break
    a <- a + as.integer((core_a - a) * 0.2); b <- b - as.integer((b - core_b) * 0.2)
  }
  list(r$chr, a, b)
}

idx <- list()
for (i in seq_len(nrow(u))) {
  r <- u[i]
  gene <- if (!is.na(r$switch_gene)) sub(",.*", "", r$switch_gene) else if (!is.na(r$link_gene)) r$link_gene else r$nearest_gene
  w <- window_for(r)
  part <- lp[dmr == r$dmr][!duplicated(partner)]
  where <- if (!is.na(r$switch_part)) paste0("switch region (", r$switch_part, ")") else
           if (grepl("^gene body", r$context)) paste0("gene body of ", sub(",.*", "", fifelse(is.na(r$body_gene), gene, r$body_gene)), ", beyond intron 1") else
           sprintf("intergenic, %.1f kb %s of the %s TSS", abs(r$tss_distance) / 1000,
                   ifelse(r$tss_distance < 0, "upstream", "downstream"), r$nearest_gene)
  title <- sprintf("%s  |  DMR %s  |  %s", gene, r$dmr, where)
  if (!is.na(r$panels)) title <- paste0(title, "\npanels: ", gsub(",", ", ", r$panels))
  sub_t <- paste0("Significant in: ", r$calls)
  f <- file.path(OUT, r$category, sprintf("%s_%s_%d-%d.png", gsub("[^A-Za-z0-9.-]", "_", gene), r$chr, r$start, r$end))
  dir.create(dirname(f), showWarnings = FALSE)
  ok <- tryCatch({
    p <- locus_plot(L, gene, region = w, highlight = c(r$start, r$end, "DMR"), title = title,
                    subtitle = paste(strwrap(sub_t, 190), collapse = "\n"), all_genes = TRUE)
    h <- attr(p, "height_in")
    lpn <- loop_panel(r$chr, r$start, r$end, part)
    if (!is.null(lpn)) {                       # stacked at exact heights, so the locus part keeps its layout
      p <- gridExtra::arrangeGrob(patchwork::patchworkGrob(p), ggplotGrob(lpn), ncol = 1,
                                  heights = grid::unit(c(h, attr(lpn, "height_in")), "in"))
      h <- h + attr(lpn, "height_in")
    }
    ggsave(f, p, width = 14, height = h, dpi = 160)
    TRUE
  }, error = function(e) { message("  ", r$dmr, ": ", conditionMessage(e)); FALSE })
  idx[[i]] <- data.table(category = r$category, gene, dmr = r$dmr, context = r$context, panels = r$panels,
                         loops = nrow(part), calls = r$calls, file = if (ok) sub(paste0(INT, "/"), "", f, fixed = TRUE) else NA)
}
fwrite(rbindlist(idx), file.path(OUT, "index.tsv"), sep = "\t", na = "NA")
message("written to ", OUT)
