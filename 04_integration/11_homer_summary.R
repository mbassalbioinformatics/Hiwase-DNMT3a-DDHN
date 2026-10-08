#!/usr/bin/env Rscript
###############################################################################
# 11_homer_summary.R                                                 2026-10-08
#
# Known-motif enrichment (HOMER findMotifsGenome.pl, run on the server by
# 01_methylation_server/04_homer_known_motifs.sh: known motifs, CpG-normalised
# background) in the metilene
# DMRs that gained or lost methylation, per comparison.
#
# With tens of DMRs per set, few motifs can pass multiple testing; the table
# keeps every motif at raw p <= 0.01 with its q-value, and the figure shows
# those motifs across all sets.
#
# Output in integration/homer/: known_motifs_p01.tsv, known_motifs_dotplot.png
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source("config/paths.R")
H    <- file.path(EMS, "03_dmr_analysis", "homer")
COMP <- COMPARISONS
OUT  <- file.path(INT, "homer"); dir.create(OUT, showWarnings = FALSE)
comp <- fread(COMP, colClasses = "character")

sets <- list.dirs(H, recursive = FALSE, full.names = TRUE)
sets <- sets[file.exists(file.path(sets, "knownResults.txt"))]
k <- rbindlist(lapply(sets, function(d) {
  x <- fread(file.path(d, "knownResults.txt"), sep = "\t", quote = "")
  setnames(x, c("motif", "consensus", "p", "logp", "q", "n_target", "pct_target", "n_bg", "pct_bg")[seq_len(ncol(x))])
  n <- as.integer(sub(".*of ([0-9]+)\\).*", "\\1", names(fread(file.path(d, "knownResults.txt"), nrows = 0))[6]))
  x[, `:=`(set = basename(d), n_dmrs = n,
           pct_target = as.numeric(sub("%", "", pct_target)), pct_bg = as.numeric(sub("%", "", pct_bg)))]
  x[, `:=`(id = sub("_(gained|lost)$", "", set), direction = sub(".*_", "", set), p = 10^logp)]
}))
k[, motif_short := sub("/.*", "", motif)]
k <- merge(k, comp[, .(id, label, status, case)], by = "id")
keep <- k[p <= 0.01]
fwrite(keep[order(id, direction, p), .(id, label, status, direction, n_dmrs, motif = motif_short, consensus,
                                        p = signif(p, 3), q, pct_target, pct_bg, enrichment = round(pct_target / pmax(pct_bg, 0.01), 2))],
       file.path(OUT, "known_motifs_p01.tsv"), sep = "\t")
message("sets: ", length(sets), "; motifs at p <= 0.01: ", nrow(keep), "; at q < 0.05: ", nrow(k[q < 0.05]))
print(keep[order(p)][1:min(.N, 40), .(label, direction, n = n_dmrs, motif = substr(motif_short, 1, 38), p = signif(p, 2), q,
                                       target = pct_target, bg = pct_bg)])

if (nrow(keep)) {
  mot <- unique(keep$motif_short)
  d <- k[motif_short %in% mot]
  d[, col := factor(paste0(label, ifelse(status == "provisional", " *", ""), "\n", direction, " (", n_dmrs, ")"))]
  p <- ggplot(d[p <= 0.05], aes(col, motif_short, size = -log10(p), colour = log2(pmax(pct_target, 0.1) / pmax(pct_bg, 0.1)))) +
    geom_point() + scale_colour_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", name = "log2 enrichment") +
    scale_size_continuous(range = c(1, 6), name = "-log10 p") +
    labs(x = NULL, y = NULL, title = "Known motifs in DMRs (HOMER; motifs at p <= 0.01 in at least one set; p <= 0.05 shown)",
         subtitle = "number of DMRs in brackets; * provisional comparison") +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 7))
  ggsave(file.path(OUT, "known_motifs_dotplot.png"), p, width = 16, height = max(5, 0.16 * length(mot) + 3), dpi = 220, limitsize = FALSE)
}
message("written to ", OUT)
