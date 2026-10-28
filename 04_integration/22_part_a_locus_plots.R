#!/usr/bin/env Rscript
###############################################################################
# 22_part_a_locus_plots.R                                            2026-10-19
#
# A locus plot (locus_plot_functions.R) for every gene that the tables of Section 12
# of the reports name, so that each one can be looked at:
#   12.2, 12.3, 12.5   the leading genes (the first six) of every program in the
#                section's programs table: the displayed programs at FDR < 0.05 in
#                methylation, by |NES| (12.5 adds those at FDR < 0.05 in expression),
#                at most 12 -- the same rule as prog_tab() in 12_part_a.R
#   12.4         the genes with the most tiles of each trajectory class (the first
#                ten), as in the genes-and-programs table of 12.4
# Methylation lost and gained alike (lost: NES < 0 / a falling class).
# Each plot is centred on the gene's genome tile (TILE_BP) with the largest change
# in the table's direction, among the tiles within 50 kb of its TSS that the
# program analysis gave it (19_clinical_questions.R); for 12.4, the tile of that
# class with the largest change. The tile is boxed and labelled. The window is the
# tile +- 1 kb, shrunk until it holds at most 90 CpGs.
#
# Output: integration/locus_plots/part_a/<section>/<lost|gained>/<gene>.png
#         integration/locus_plots/part_a/index.tsv (gene, programs, tile, change)
###############################################################################
source("config/paths.R")
ID <- with(data.table::fread(SAMPLE_SHEET, colClasses = "character"), setNames(sample, key))   # role key -> sample ID (config/samples.tsv)

suppressPackageStartupMessages(library(data.table))
source("04_integration/locus_plot_functions.R")
CQD <- file.path(INT, "clinical")
OUT <- file.path(INT, "locus_plots", "part_a")
rd  <- function(f, ...) if (grepl("\\.gz$", f)) fread(cmd = sprintf("gzip -dc '%s'", f), ...) else fread(f, ...)
WINDOW <- 50e3; N_LEAD <- 6L; N_PROG <- 12L; N_A4 <- 10L

smp  <- fread(SAMPLE_SHEET, colClasses = "character")[methylation == "analysed"]
nm   <- function(s) smp[sample == s, plot_label]
prog <- fread(file.path(CQD, "programs.tsv"))
show <- fread(file.path(CQD, "programs_shown.tsv"))
SECT <- data.table(section = c("12.2", "12.3", "12.5"), id = c(paste0(ID[["D23"]], "_vs_", ID[["D18"]]), paste0(ID[["D30"]], "_vs_", ID[["D18"]]), paste0(ID[["DX"]], "_vs_", ID[["D23"]])),
                   rna = c(NA, NA, "C4_DDHN_vs_donor"),
                   folder = c("12.2_donor_5_years_vs_donation", "12.3_recipient_1_month_vs_donor_donation", "12.5_DDHN_vs_donor_5_years"))

## ---- the genes of the tables ----------------------------------------------------------------------
want <- rbindlist(lapply(seq_len(nrow(SECT)), function(k) {
  s <- SECT[k]
  m <- prog[layer == "methylation" & id == s$id & set %in% show$set & padj < 0.05][order(-abs(NES))]
  sets <- m$set
  if (!is.na(s$rna)) sets <- unique(c(sets, prog[layer == "RNA-seq" & id == s$rna & set %in% show$set & padj < 0.05][order(-abs(NES)), set]))
  sets <- head(sets, N_PROG)
  rbindlist(lapply(sets, function(st) {
    x <- prog[layer == "methylation" & id == s$id & set == st][1]
    data.table(section = s$section, folder = s$folder, id = s$id, what = show[set == st, label],
               nes = x$NES, direction = fifelse(x$NES < 0, "lost", "gained"),
               gene = head(strsplit(x$lead, ", ")[[1]], N_LEAD))
  }))
}))
twn <- fread(file.path(CQD, "three_way_named_genes.tsv"))
A4C <- data.table(class = c("with the clone and following it, lost", "with the clone and following it, gained", "at diagnosis, lost",
                            "at diagnosis, gained", "progressive, lost"),
                  traj = c("with the clone", "with the clone", "at diagnosis", "at diagnosis", "progressive"),
                  dir = c("lost", "gained", "lost", "gained", "lost"), follows = c(TRUE, TRUE, FALSE, FALSE, FALSE))
want <- rbind(want, rbindlist(lapply(seq_len(nrow(A4C)), function(k) {
  g <- head(twn[class == A4C$class[k], gene], N_A4)
  if (!length(g)) return(NULL)
  data.table(section = "12.4", folder = "12.4_recipient_1_month_3.5_years_DDHN", id = paste(ID[["D30"]], ID[["Y35"]], ID[["DX"]], sep = ","), what = A4C$class[k],
             nes = NA_real_, direction = A4C$dir[k], gene = g)
})))
message(nrow(want), " table entries, ", uniqueN(want[, .(folder, direction, gene)]), " plots")

## ---- tiles: each to the nearest MANE TSS within 50 kb (as in 19_clinical_questions.R) ------------------
ti <- rd(TILE_FILE)
ti[, `:=`(chr = sub(":.*", "", id), start = as.integer(sub(".*:", "", id)))][, mid := start + TILE_BP %/% 2L - 1L]
gn <- rd(file.path(INT, "region_methylation", "mane_select_genes.tsv.gz"))[!duplicated(gene), .(gene, chr, tss, strand)]
setkey(gn, chr, tss)
near <- gn[ti[, .(id, chr, tss = mid)], on = .(chr, tss), roll = "nearest", .(id, gene, dist = abs(x.tss - i.tss))]
ti <- merge(ti, near[!is.na(gene) & dist <= WINDOW, .(id, gene)], by = "id")
tcl <- rd(file.path(CQD, "three_way_tile_classes.tsv.gz"))

pick_tile <- function(w) {                    # the tile to centre on, with its change in points and a label
  g <- w$gene[1]; dirn <- w$direction[1]
  if (w$section[1] == "12.4") {
    k <- A4C[class == w$what[1]]
    x <- tcl[gene == g & traj == k$traj & dir == k$dir & (!k$follows | check %in% "follows the clone")]
    if (!nrow(x)) return(NULL)
    x[, ch := if (k$traj == "with the clone") a else if (k$traj == "at diagnosis") b else a + b]
    x <- x[order(if (dirn == "lost") ch else -ch)][1]
    step <- c(`with the clone` = "3.5 years vs 1 month", `at diagnosis` = "DDHN diagnosis vs 3.5 years", progressive = "DDHN diagnosis vs 1 month")[[k$traj]]
    return(list(id = x$id, change = 100 * x$ch, step = step))
  }
  p <- strsplit(w$id[1], "_vs_")[[1]]
  x <- ti[gene == g, .(id, ch = get(p[1]) - get(p[2]))][!is.na(ch)]
  if (!nrow(x)) return(NULL)
  x <- x[order(if (dirn == "lost") ch else -ch)][1]
  list(id = x$id, change = 100 * x$ch, step = paste(nm(p[1]), "vs", nm(p[2])))
}

## ---- plots ------------------------------------------------------------------------------------------
L <- load_locus_data(emseq = EMS, mane_gtf = MANE, cgi_file = file.path(RES, "cpgIslandExt.txt.gz"),
                     rna_rdata = file.path(RNA, "05_processed_counts_genes", "dge_gene_level.RData"))
has_model <- unique(L$model[feature == "transcript", gene])
# one plot per section, direction and gene; the first table row decides the tile (for A4: its class)
want[, label := fifelse(is.na(nes), what, sprintf("%s (NES %+.2f)", what, nes))]
jobs <- want[, .(what = what[1], nes = nes[1], in_tables = paste(unique(label), collapse = "; ")),
             by = .(section, folder, id, direction, gene)]
idx <- vector("list", nrow(jobs))
for (j in seq_len(nrow(jobs))) {
  w <- jobs[j]; g <- w$gene
  row <- data.table(section = w$section, direction = w$direction, gene = g, in_tables = w$in_tables, tile = NA_character_,
                    change_points = NA_real_, step = NA_character_, tss_distance_kb = NA_real_, file = NA_character_, note = "")
  tl <- if (g %in% has_model) pick_tile(w) else NULL
  if (is.null(tl)) { row$note <- if (g %in% has_model) "no tile of this gene in this direction" else "no MANE Select model"
                     idx[[j]] <- row; next }
  ch <- sub(":.*", "", tl$id); t0 <- as.integer(sub(".*:", "", tl$id)); t1 <- t0 + TILE_BP - 1L
  gx <- gn[gene == g][1]
  dkb <- ((t0 + t1) / 2 - gx$tss) * fifelse(gx$strand == "+", 1, -1) / 1e3
  half <- 1000L; cx <- (t0 + t1) %/% 2L
  cpg1 <- L$cpg[[1]][chr == ch]
  repeat { n <- cpg1[pos >= cx - TILE_BP %/% 2L - half & pos <= cx + TILE_BP %/% 2L + half, .N]
           if (n <= 90 || half <= 250) break; half <- as.integer(half * 0.85) }
  region <- list(ch, cx - TILE_BP %/% 2L - half, cx + TILE_BP %/% 2L + half)
  where <- if (abs(dkb) < 0.5) "at the TSS" else sprintf("%.1f kb %s of the TSS", abs(dkb), if (dkb < 0) "upstream" else "downstream")
  ttl <- sprintf("%s  -  Section %s, methylation %s", g, w$section, w$direction)
  subt <- paste0("In the table as: ", w$in_tables, "\n",
                sprintf("Boxed: the %s tile with the largest %s, %+.0f points (%s), %s", TILE_TXT,
                        if (w$direction == "lost") "loss" else "gain", tl$change, tl$step, where),
                "\nMethylation per CpG (TAPS)")
  f <- file.path(OUT, w$folder, w$direction, paste0(g, ".png"))
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
  ok <- tryCatch({
    p <- locus_plot(L, g, region = region, highlight = c(t0, t1, sprintf("%+.0f points", tl$change)), title = ttl, subtitle = subt, all_genes = TRUE)
    ggsave(f, p, width = 14, height = attr(p, "height_in"), dpi = 170); TRUE
  }, error = function(e) { message("  ", g, ": ", conditionMessage(e)); FALSE })
  row[, `:=`(tile = tl$id, change_points = round(tl$change, 1), step = tl$step, tss_distance_kb = round(dkb, 1),
             file = if (ok) substring(f, nchar(DATA_ROOT) + 2L) else NA_character_, note = if (ok) "" else "plot failed")]
  idx[[j]] <- row
}
idx <- rbindlist(idx)
fwrite(idx, file.path(OUT, "index.tsv"), sep = "\t")
message(sum(!is.na(idx$file)), " of ", nrow(idx), " plots written to ", OUT)
if (any(is.na(idx$file))) print(idx[is.na(file), .(section, direction, gene, note)])
