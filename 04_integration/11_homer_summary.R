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
# Second part: where the DMRs fall in the genome (HOMER annotatePeaks.pl hg38
# -annStats on the metilene and "both" DMR sets, run on the server by
# 01_methylation_server/05_homer_annotation.sh). HOMER assigns each DMR one category by priority
# (promoter-TSS = TSS -1 kb to +100 bp; TTS = -100 bp to +1 kb) and reports, per
# set, the observed / expected (genome share) ratio per category.
#
# Output in integration/homer/: known_motifs_p01.tsv, known_motifs_dotplot.png,
# dmr_annotation_counts.tsv, dmr_annotation_enrichment.tsv, dmr_annotation.png
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
         subtitle = paste0("number of DMRs in brackets", if (any(comp$status == "provisional")) "; * provisional comparison" else "")) +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 7))
  ggsave(file.path(OUT, "known_motifs_dotplot.png"), p, width = 16, height = max(5, 0.16 * length(mot) + 3), dpi = 220, limitsize = FALSE)
}

## ---- genomic annotation of the DMRs (HOMER annotatePeaks.pl) -------------------
A <- file.path(dirname(H), "homer_annotation")
if (dir.exists(A)) {
  CATS <- c(`promoter-TSS` = "Promoter-TSS", `5' UTR` = "5' UTR", exon = "Exon", intron = "Intron", `3' UTR` = "3' UTR",
            TTS = "TTS", `non-coding` = "Non-coding RNA", Intergenic = "Intergenic")
  an <- rbindlist(lapply(c("metilene", "both"), function(set) {
    fs <- list.files(file.path(A, set), pattern = "\\.annotated\\.txt$", full.names = TRUE)
    rbindlist(lapply(fs, function(f) {
      x <- fread(f, sep = "\t", quote = "", select = c(1, 2, 3, 4, 8, 9, 16), colClasses = "character")
      setnames(x, c("peak", "chr", "start", "end", "ann", "detail", "nearest_gene"))
      x[, `:=`(set = set, id = sub("\\.annotated\\.txt$", "", basename(f)))]
    }))
  }))
  an[, `:=`(direction = sub(".*_(gained|lost)_[0-9]+$", "\\1", peak), category = CATS[sub(" \\(.*", "", ann)],
            start = as.integer(start), end = as.integer(end))]
  # overlap with a CpG island (UCSC cpgIslandExt), whatever HOMER's category (HOMER labels a DMR "CpG" only when
  # no promoter, exon, UTR or TTS takes priority)
  cgi_tab <- fread(cmd = sprintf("gzip -dc '%s'", file.path(RES, "cpgIslandExt.txt.gz")), header = FALSE,
                   select = 2:4, col.names = c("chr", "start", "end"))
  cgi_tab[, start := start + 1L]                              # UCSC 0-based -> 1-based, as in 04
  setkey(cgi_tab, chr, start, end)
  an[, row := .I]
  hit <- foverlaps(an[, .(row, chr, start, end)], cgi_tab, by.x = c("chr", "start", "end"), nomatch = 0L)[, unique(row)]
  an[, cgi := row %in% hit][, row := NULL]
  stopifnot(!anyNA(an$category), all(an$direction %in% c("gained", "lost")))
  cnt <- an[, .N, by = .(id, set, direction, category)]
  cnt <- merge(cnt, an[, .(dmrs = .N, cgi = sum(cgi)), by = .(id, set, direction)], by = c("id", "set", "direction"))
  cnt[, pct := round(100 * N / dmrs, 1)]
  cnt <- merge(cnt, comp[, .(id, label, status)], by = "id")
  fwrite(cnt[order(set, id, direction, match(category, CATS))], file.path(OUT, "dmr_annotation_counts.tsv"), sep = "\t")

  # observed / expected per category from -annStats (block 1 = basic categories; CpG-Island from block 2)
  ann_stats <- function(f) {
    l <- readLines(f); h <- grep("^Annotation\t", l)
    rd <- function(rows) fread(text = paste(l[rows], collapse = "\n"), sep = "\t",
                               col.names = c("annotation", "dmrs", "size_bp", "log2_obs_exp", "logp"))
    rbind(rd(h[1]:(h[2] - 1))[, block := "basic"], rd(h[2]:length(l))[annotation == "CpG-Island"][, block := "detailed"])
  }
  en <- rbindlist(lapply(c("metilene", "both"), function(set) {
    fs <- list.files(file.path(A, set), pattern = "\\.annStats\\.txt$", full.names = TRUE)
    rbindlist(lapply(fs, function(f) ann_stats(f)[, `:=`(set = set, id = sub("\\.annStats\\.txt$", "", basename(f)))]))
  }))
  en <- merge(en, comp[, .(id, label, status)], by = "id")
  fwrite(en[order(set, id, block, -log2_obs_exp)], file.path(OUT, "dmr_annotation_enrichment.tsv"), sep = "\t")

  # figure: share of metilene DMRs per category, gained and lost, with the genome's share as reference
  short_lab <- function(l) {
    inner <- sub(".*\\((.* vs .*)\\).*", "\\1", l)
    s <- ifelse(inner != l, inner, sub(" \\(.*\\)$", "", l))
    s <- sub("^TET2 case: ", "TET2: ", sub("^Engraftment: ", "", s))
    sub("^(.)", "\\U\\1", s, perl = TRUE)
  }
  m <- cnt[set == "metilene"]
  m[, row := paste0(short_lab(label), ifelse(status == "provisional", " *", ""))]
  gen <- en[set == "metilene" & block == "basic", .(size_bp = size_bp[1]), by = annotation]
  gmap <- c(Promoter = "Promoter-TSS", `5UTR` = "5' UTR", Exon = "Exon", Intron = "Intron", `3UTR` = "3' UTR", TTS = "TTS",
            ncRNA = "Non-coding RNA", miRNA = "Non-coding RNA", snoRNA = "Non-coding RNA", rRNA = "Non-coding RNA",
            pseudo = "Intergenic", Intergenic = "Intergenic")
  gen <- gen[annotation %in% names(gmap), .(size_bp = sum(size_bp)), by = .(category = gmap[annotation])]
  gen[, pct := 100 * size_bp / sum(size_bp)]
  ref <- rbind(gen[, .(row = "Genome (share of bp)", direction = "gained", category, pct)],
               gen[, .(row = "Genome (share of bp)", direction = "lost", category, pct)])
  ord <- comp[match(unique(m$id), id)][order(match(id, comp$id)), id]
  rows <- c("Genome (share of bp)", rev(unique(m[match(ord, id), row])))
  dd <- rbind(m[, .(row, direction, category, pct)], ref)
  dd[, `:=`(row = factor(row, levels = rows), category = factor(category, levels = rev(CATS)),
            direction = factor(paste("Methylation", direction), levels = c("Methylation lost", "Methylation gained")))]
  cgi <- m[, .(txt = sprintf("%d DMRs, %d%% CGI", dmrs[1], round(100 * cgi[1] / dmrs[1]))), by = .(row, direction)]
  none <- CJ(row = unique(m$row), direction = c("gained", "lost"))[!cgi, on = .(row, direction)][, txt := "no DMRs"]
  cgi <- rbind(cgi, none)
  cgi[, `:=`(row = factor(row, levels = rows), direction = factor(paste("Methylation", direction), levels = levels(dd$direction)))]
  pal <- c(`Promoter-TSS` = "#B2182B", `5' UTR` = "#EF8A62", Exon = "#F4C27A", Intron = "#92C5DE", `3' UTR` = "#D1A9E0",
           TTS = "#9E9AC8", `Non-coding RNA` = "#A6D96A", Intergenic = "#BDBDBD")
  p <- ggplot(dd, aes(pct, row, fill = category)) +
    geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
    geom_text(data = cgi, aes(x = 102, y = row, label = txt), inherit.aes = FALSE, hjust = 0, size = 3.2) +
    facet_wrap(~ direction, ncol = 2) +
    scale_fill_manual(values = pal, breaks = unname(CATS), name = "HOMER category") +
    scale_x_continuous(limits = c(0, 135), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    labs(x = "% of DMRs", y = NULL,
         title = "Where the DMRs fall (metilene DMRs, HOMER annotatePeaks.pl, hg38)",
         subtitle = paste("right of each bar: number of DMRs and the share overlapping a CpG island (CGI);",
                          paste0("bottom bar: share of the genome in each category", if (any(comp$status == "provisional")) "; * provisional comparison" else ""), sep = "\n")) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
          axis.text.y = element_text(size = 10)) +
    guides(fill = guide_legend(nrow = 1))
  ggsave(file.path(OUT, "dmr_annotation.png"), p, width = 16, height = 2.2 + 0.42 * length(rows), dpi = 200)
  message("DMR annotation: ", nrow(an), " DMRs in ", uniqueN(an[, .(id, set)]), " sets")
}
message("written to ", OUT)

