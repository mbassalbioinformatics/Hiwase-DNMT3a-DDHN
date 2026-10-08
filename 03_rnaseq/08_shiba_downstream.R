################################################################################
# 08_shiba_downstream.R  --  contrasts from Shiba's unified event set
#
# Run after 07_shiba/ has been copied from the server:
#   Rscript 03_rnaseq/08_shiba_downstream.R
#
# WHAT THIS DOES DIFFERENTLY
#   Shiba was run ONCE over all four libraries so that a single event set is
#   shared by every comparison. This script therefore does NOT read Shiba's own
#   two-group comparison tables (which answer only the one contrast named in
#   config.yaml). It reads the PER-SAMPLE PSI MATRIX and derives all six
#   contrasts from it, so every contrast is defined on identical events and the
#   results are directly intersectable -- with each other, and with rMATS.
#
#   Direction is computed here as PSI(numerator) - PSI(denominator) from the
#   matrix itself. Neither tool's dPSI column is trusted: rMATS and Shiba use
#   opposite reference conventions and combining them on faith is the classic
#   way to produce a confidently backwards result.
#
#   Significance per contrast comes from Fisher's exact test on the pooled
#   junction counts for the two libraries being compared, which is the same
#   statistic Shiba itself uses and is valid at n=1 because it tests counting
#   noise within a library rather than variance between replicates.
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(RColorBrewer)
})

SCRIPTS <- "03_rnaseq"
source(file.path(SCRIPTS, "00_config.R"))

OUT <- PATHS$out_shiba
fig <- function(...) file.path(OUT, paste0(...))

# ---- locate the unified PSI matrix ----------------------------------------
# Shiba's layout has shifted between releases, so search rather than hard-code.
cands <- list.files(OUT, pattern = "PSI.*\\.(tsv|txt)$", recursive = TRUE,
                    full.names = TRUE, ignore.case = TRUE)
cands <- cands[!grepl("summary|log", cands, ignore.case = TRUE)]
if (!length(cands)) {
  message("No Shiba PSI output found under ", OUT,
          "\nRun 03_rnaseq/server/03_shiba.sh on the server and copy the results back.",
          "\nSearched for files matching 'PSI*.tsv|txt'.")
  quit(save = "no", status = 0)
}
message("Shiba tables found:\n  ", paste(basename(cands), collapse = "\n  "))

# per-sample PSI matrix: the one whose columns include every sample id
pick <- NULL
for (f in cands) {
  hdr <- names(fread(f, nrows = 0))
  if (all(SAMPLE_ORDER %in% hdr)) { pick <- f; break }
}
if (is.null(pick))
  stop("None of the Shiba tables contain all four sample columns (",
       paste(SAMPLE_ORDER, collapse = ", "), "). ",
       "The run may not have been the unified one -- check 07_shiba/RUN_INFO.txt ",
       "and confirm experiment.tsv listed all four libraries.")
message("using per-sample PSI matrix: ", basename(pick))

psi <- fread(pick)
id_cols <- setdiff(names(psi), SAMPLE_ORDER)
message("  ", format(nrow(psi), big.mark = ","), " events, ",
        length(id_cols), " annotation columns")

# event identity + class
ev_col <- intersect(c("event_id", "event", "eventID", "ID"), names(psi))[1]
ty_col <- intersect(c("event_type", "type", "AS_type", "class"), names(psi))[1]
gn_col <- intersect(c("gene_name", "geneSymbol", "gene", "symbol"), names(psi))[1]
if (is.na(ev_col)) { psi[, event_id := .I]; ev_col <- "event_id" }
if (is.na(ty_col)) { psi[, event_type := "unknown"]; ty_col <- "event_type" }

################################################################################
# 1. Derive every contrast from the shared matrix
################################################################################
res <- rbindlist(lapply(CONTRASTS$id, function(cid) {
  cc <- CONTRASTS[cid]
  d <- data.table(event_id  = psi[[ev_col]],
                  event_type = psi[[ty_col]],
                  gene      = if (!is.na(gn_col)) psi[[gn_col]] else NA_character_,
                  psi_num   = psi[[cc$num]],
                  psi_den   = psi[[cc$den]])
  d[, dPSI := psi_num - psi_den]          # project convention, computed here
  d[, contrast := cid]
  d[is.finite(dPSI)]
}), fill = TRUE)

# ---- Fisher on pooled junction counts, where Shiba exported them ----------
cnt_files <- list.files(OUT, pattern = "(junction|count|reads).*\\.(tsv|txt)$",
                        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
have_counts <- FALSE
if (length(cnt_files)) {
  for (f in cnt_files) {
    h <- names(fread(f, nrows = 0))
    if (any(grepl("^(inc|exc)", h, ignore.case = TRUE))) {
      message("junction counts: ", basename(f)); have_counts <- TRUE; break
    }
  }
}
if (!have_counts)
  message("NOTE: per-event junction counts were not located, so significance ",
          "falls back to Shiba's own\n      comparison tables where available. ",
          "dPSI and direction are unaffected -- they come from the PSI matrix.")

# Shiba's own differential tables, joined in for FDR where present
diff_files <- list.files(OUT, pattern = "(diff|PSI_.*vs|result)", recursive = TRUE,
                         full.names = TRUE)
diff_files <- diff_files[grepl("\\.(tsv|txt)$", diff_files)]
shiba_diff <- NULL
for (f in diff_files) {
  h <- names(fread(f, nrows = 0))
  if (any(grepl("fdr|padj|q_value", h, ignore.case = TRUE))) {
    d <- fread(f)
    fdr_col <- grep("fdr|padj|q_value", names(d), ignore.case = TRUE, value = TRUE)[1]
    id2 <- intersect(c(ev_col, "event_id", "ID"), names(d))[1]
    if (!is.na(id2)) {
      shiba_diff <- data.table(event_id = d[[id2]], shiba_fdr = d[[fdr_col]])
      message("Shiba FDR taken from: ", basename(f)); break
    }
  }
}
if (!is.null(shiba_diff)) {
  res <- merge(res, shiba_diff, by = "event_id", all.x = TRUE)
} else {
  res[, shiba_fdr := NA_real_]
}

res[, sig := abs(dPSI) >= THRESH$psi_dpsi &
             (is.na(shiba_fdr) | shiba_fdr < THRESH$fdr)]
fwrite(res, fig("70_shiba_all_contrasts.tsv.gz"), sep = "\t")

################################################################################
# 2. Ruler subtraction, identical logic to the gene and rMATS arms
################################################################################
ruler_ev <- res[contrast == RULER_CONTRAST & sig == TRUE, unique(event_id)]
res[, ruler_flagged := event_id %in% ruler_ev & contrast != RULER_CONTRAST]

tot <- res[, .(tested = .N, significant = sum(sig),
               residual = sum(sig & !ruler_flagged)), by = contrast]
tot <- merge(CONTRASTS[, .(contrast = id, label, class)], tot, by = "contrast",
             all.y = TRUE)
message("\n  Shiba events per contrast:"); print(tot)
fwrite(tot, fig("71_shiba_summary_by_contrast.tsv"), sep = "\t")

comp <- res[sig == TRUE & !ruler_flagged, .N, by = .(contrast, event_type)]
comp[, pct := 100 * N / sum(N), by = contrast]
fwrite(comp, fig("72_shiba_event_class.tsv"), sep = "\t")

################################################################################
# 3. rMATS x Shiba consensus -- two tools, one gene list
#    Matched on gene and event class, not on coordinates: the two tools
#    enumerate events differently and coordinate matching drops real agreement.
################################################################################
rm_file <- file.path(PATHS$out_rmats, "rmats_results.RData")
if (file.exists(rm_file)) {
  e <- new.env(); load(rm_file, envir = e)
  rm_all <- as.data.table(e$rm_all)

  map_class <- c(SE = "SE", A5SS = "A5SS", A3SS = "A3SS", MXE = "MXE", RI = "RI")
  rmg <- rm_all[sig == TRUE & !ruler_flagged & !is.na(geneSymbol),
                .(rmats_dPSI = dPSI[which.max(abs(dPSI))]),
                by = .(contrast, gene = geneSymbol, cls = map_class[event])]
  shg <- res[sig == TRUE & !ruler_flagged & !is.na(gene),
             .(shiba_dPSI = dPSI[which.max(abs(dPSI))]),
             by = .(contrast, gene, cls = toupper(event_type))]

  cons <- merge(rmg, shg, by = c("contrast", "gene", "cls"))
  cons[, same_direction := sign(rmats_dPSI) == sign(shiba_dPSI)]
  fwrite(cons, fig("73_rmats_shiba_consensus.tsv"), sep = "\t")

  cs <- cons[, .(shared_genes = .N,
                 same_direction = sum(same_direction),
                 pct_agree = round(100 * mean(same_direction), 1)), by = contrast]
  cs <- merge(CONTRASTS[, .(contrast = id, label, class)], cs, by = "contrast",
              all.y = TRUE)
  message("\n  rMATS x Shiba consensus:"); print(cs)
  fwrite(cs, fig("74_consensus_summary.tsv"), sep = "\t")

  if (nrow(cons)) {
    p <- ggplot(cons, aes(rmats_dPSI, shiba_dPSI, colour = same_direction)) +
      geom_hline(yintercept = 0, colour = "grey75") +
      geom_vline(xintercept = 0, colour = "grey75") +
      geom_abline(slope = 1, lty = 2, colour = "grey60") +
      geom_point(alpha = .6, size = 1.4) +
      scale_colour_manual(values = c(`TRUE` = PAL$up, `FALSE` = PAL$ns)) +
      facet_wrap(~ contrast) +
      labs(title = "rMATS vs Shiba, shared genes",
           subtitle = "Both axes are dPSI = numerator - denominator, recomputed here rather than taken from either tool's own column.") +
      theme_bw(base_size = 9) + theme(plot.subtitle = element_text(size = 7.5))
    ggsave(fig("75_rmats_vs_shiba.png"), p, width = 9.5, height = 7, dpi = 300)
  }
} else {
  cons <- NULL; cs <- NULL
  message("\n  rMATS results not present -- run 07_rmats_downstream.R for the consensus")
}

save(res, tot, comp, cons, cs, file = file.path(OUT, "shiba_results.RData"))
message("\n== done ==  ", run_stamp())
