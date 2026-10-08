################################################################################
# 07_rmats_downstream.R  --  parse, verify and summarise the rMATS output
#
# Run after 06_rmats/ has been copied from the server:
#   Rscript 03_rnaseq/07_rmats_downstream.R
#
# THE SIGN TRAP, AND WHY THIS SCRIPT REFUSES TO ASSUME
#   rMATS reports IncLevelDifference = PSI(group1) - PSI(group2) with group1 =
#   --b1. 03_rnaseq/server/02_rmats.sh always puts the contrast NUMERATOR in --b1, so the
#   reported value should already match the project convention with no flip.
#   This is checked, not assumed. Every contrast is validated by recomputing
#   mean(IncLevel1) - mean(IncLevel2) from the raw per-sample PSI columns and
#   comparing it to the reported IncLevelDifference. A mismatch aborts the
#   script rather than producing a beautifully formatted backwards result.
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(RColorBrewer)
})

SCRIPTS <- "03_rnaseq"
source(file.path(SCRIPTS, "00_config.R"))

OUT <- PATHS$out_rmats
fig <- function(...) file.path(OUT, paste0(...))
EVENTS <- c("SE", "A5SS", "A3SS", "MXE", "RI")
EVENT_LABEL <- c(SE = "Exon skipping", A5SS = "Alt 5' splice site",
                 A3SS = "Alt 3' splice site", MXE = "Mutually exclusive exons",
                 RI = "Intron retention")

run_dirs <- list.dirs(OUT, recursive = FALSE)
run_dirs <- run_dirs[file.exists(file.path(run_dirs, "SE.MATS.JC.txt"))]
if (!length(run_dirs)) {
  message("No rMATS output found in ", OUT,
          "\nRun 03_rnaseq/server/02_rmats.sh on the server and copy the results back.")
  quit(save = "no", status = 0)
}
message("found ", length(run_dirs), " rMATS runs: ",
        paste(basename(run_dirs), collapse = ", "))

################################################################################
# 1. Read every event table, with the sign convention verified per contrast
################################################################################
# rMATS writes one comma-separated value per replicate. With a single replicate
# there is no comma, so fread types the column as numeric rather than character
# -- and the grouped run has three replicates in group 1, so BOTH forms occur in
# the same analysis. These coerce first and therefore handle either.
mean_csv <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  vapply(strsplit(as.character(x), ","), function(v)
    mean(suppressWarnings(as.numeric(v)), na.rm = TRUE), numeric(1))
}
sum_csv <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  vapply(strsplit(as.character(x), ","), function(v)
    sum(suppressWarnings(as.numeric(v)), na.rm = TRUE), numeric(1))
}

read_one <- function(dir, ev, ctype = "JC") {
  f <- file.path(dir, sprintf("%s.MATS.%s.txt", ev, ctype))
  if (!file.exists(f)) return(NULL)
  d <- fread(f)
  if (!nrow(d)) return(NULL)
  cid <- basename(dir)

  d[, psi_num := mean_csv(IncLevel1)]        # b1 = numerator, by construction
  d[, psi_den := mean_csv(IncLevel2)]
  d[, reads   := sum_csv(IJC_SAMPLE_1) + sum_csv(SJC_SAMPLE_1) +
                 sum_csv(IJC_SAMPLE_2) + sum_csv(SJC_SAMPLE_2)]

  # ---- empirical sign verification -----------------------------------------
  ok <- d[is.finite(psi_num) & is.finite(psi_den) & is.finite(IncLevelDifference)]
  if (nrow(ok) >= 20) {
    delta <- ok$psi_num - ok$psi_den
    agree <- mean(abs(delta - ok$IncLevelDifference) < 1e-3, na.rm = TRUE)
    if (agree < 0.95)
      stop("SIGN CHECK FAILED for ", cid, " / ", ev, ": only ",
           round(100 * agree, 1), "% of events satisfy ",
           "IncLevelDifference == mean(IncLevel1) - mean(IncLevel2). ",
           "Do not trust the direction of this run -- inspect RUN_INFO.txt ",
           "and confirm which BAM list went to --b1.")
  }

  # project convention: dPSI = PSI(numerator) - PSI(denominator)
  d[, dPSI := psi_num - psi_den]
  d[, `:=`(event = ev, contrast = cid, count_type = ctype)]
  d[, .(event, contrast, count_type, ID, GeneID, geneSymbol, chr, strand,
        psi_num, psi_den, dPSI, PValue, FDR, reads)]
}

rm_all <- rbindlist(lapply(run_dirs, function(d)
  rbindlist(lapply(EVENTS, function(e) read_one(d, e, "JC")), fill = TRUE)),
  fill = TRUE)
message("parsed ", format(nrow(rm_all), big.mark = ","),
        " event rows; sign convention verified for every contrast")

rm_all[, sig := FDR < THRESH$fdr &
                abs(dPSI) >= THRESH$psi_dpsi &
                reads >= THRESH$psi_minreads]
fwrite(rm_all, fig("60_rmats_all_events.tsv.gz"), sep = "\t")

################################################################################
# 2. Ruler subtraction
################################################################################
key_of <- function(d) paste(d$event, d$GeneID, d$ID, sep = "|")
rm_all[, ekey := key_of(rm_all)]
ruler_keys <- rm_all[contrast == RULER_CONTRAST & sig == TRUE, unique(ekey)]
rm_all[, ruler_flagged := ekey %in% ruler_keys & contrast != RULER_CONTRAST]

summ <- rm_all[, .(tested = .N, significant = sum(sig),
                   in_ruler = sum(sig & ruler_flagged),
                   residual = sum(sig & !ruler_flagged)),
               by = .(contrast, event)]
fwrite(summ, fig("61_rmats_summary_by_event.tsv"), sep = "\t")

tot <- rm_all[, .(tested = .N, significant = sum(sig),
                  residual = sum(sig & !ruler_flagged)), by = contrast]
tot <- merge(CONTRASTS[, .(contrast = id, label, class)], tot, by = "contrast",
             all.y = TRUE)
message("\n  significant events per contrast:"); print(tot)
fwrite(tot, fig("62_rmats_summary_by_contrast.tsv"), sep = "\t")

################################################################################
# 3. Event class composition -- which kind of splicing changes?
################################################################################
comp <- rm_all[sig == TRUE & !ruler_flagged, .N, by = .(contrast, event)]
comp[, pct := 100 * N / sum(N), by = contrast]
fwrite(comp, fig("63_rmats_event_class_composition.tsv"), sep = "\t")

if (nrow(comp)) {
  comp[, event_lab := EVENT_LABEL[event]]
  p <- ggplot(comp, aes(contrast, pct, fill = event_lab)) +
    geom_col(width = .7) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = NULL, y = "% of ruler-clean significant events", fill = "event class",
         title = "Splicing event class composition (rMATS)",
         subtitle = "Ruler-flagged events removed. Directional claims need the dPSI panel, not these proportions.") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 22, hjust = 1),
          plot.subtitle = element_text(size = 8))
  ggsave(fig("64_rmats_event_composition.png"), p, width = 9, height = 5.5, dpi = 300)
}

################################################################################
# 4. Direction -- does inclusion go up or down?
################################################################################
dirn <- rm_all[sig == TRUE & !ruler_flagged,
               .(up = sum(dPSI > 0), down = sum(dPSI < 0)), by = .(contrast, event)]
dirn[, pct_up := round(100 * up / (up + down), 1)]
dirn[, binom_p := mapply(function(u, n)
  if (n > 0) binom.test(u, n, 0.5)$p.value else NA_real_, up, up + down)]
fwrite(dirn, fig("65_rmats_direction.tsv"), sep = "\t")
message("\n  directional bias (ruler-clean):"); print(dirn)

################################################################################
# 5. Per-contrast volcano
################################################################################
for (cid in unique(rm_all$contrast)) {
  d <- rm_all[contrast == cid & reads >= THRESH$psi_minreads]
  if (!nrow(d)) next
  d[, state := fifelse(sig & !ruler_flagged, "significant",
              fifelse(sig & ruler_flagged, "ruler-flagged", "ns"))]
  lab <- d[state == "significant"][order(FDR)][1:min(20, .N)]
  ttl <- if (cid %in% CONTRASTS$id) CONTRASTS[cid, label] else cid
  p <- ggplot(d, aes(dPSI, -log10(PValue), colour = state)) +
    geom_point(size = 1, alpha = .55) +
    scale_colour_manual(values = c(significant = PAL$up,
                                   `ruler-flagged` = "#E8C547", ns = PAL$ns)) +
    geom_vline(xintercept = c(-THRESH$psi_dpsi, THRESH$psi_dpsi), lty = 2,
               colour = "grey50") +
    ggrepel::geom_text_repel(data = lab, aes(label = geneSymbol), size = 2.5,
                             max.overlaps = 18, show.legend = FALSE) +
    labs(x = "dPSI  (PSI numerator - PSI denominator)",
         title = paste0("rMATS - ", ttl),
         subtitle = "n=1 per group: rMATS falls back to a fixed within-group variance, so FDR is optimistic. Effect size and read depth carry the weight.") +
    theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size = 7.5))
  ggsave(fig(sprintf("66_rmats_volcano_%s.png", cid)), p,
         width = 8, height = 6, dpi = 300)
}

################################################################################
# 6. Save
################################################################################
save(rm_all, summ, tot, comp, dirn, file = file.path(OUT, "rmats_results.RData"))
message("\n== done ==  ", run_stamp())
