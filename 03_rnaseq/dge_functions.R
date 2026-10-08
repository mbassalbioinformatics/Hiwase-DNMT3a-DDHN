################################################################################
# dge_functions.R  --  shared pipeline, sourced by 02_ (genes) and 03_ (5'UTR)
#
# This follows the laboratory's standard limma-voom workflow ("the standard
# workflow" below) section by section, so the
# two levels produce the same family of outputs under the same numbering:
#
#   1  filtering                      -> 11_dge_library_sizes.png
#   2  DGEList + TMM
#   3  MDS over all samples           -> 12_dge_mds_*.png
#   4  per-comparison MDS             -> 13_mds_<contrast>.png
#   5  heatmaps, top 50 / top 100     -> 12_dge_heatmap_*.png
#   6  DGE, all contrasts             -> 17..22_dge_<contrast>.txt
#   7a MD plots per contrast          -> 21_MD_<contrast>.png
#   7b MA plots per contrast          -> 22_MAplot_<contrast>.png/pdf
#   8  ORA                            -> enrichment/<contrast>/<dir>/
#   9  GSEA                           -> enrichment_gsea/<contrast>/
#   8-9 matched MA plot per term      -> .../ma_plots/<CAT>/  (04b_enrichment_ma_plots.R)
#
#   7b and the matched per-term plots replace the standard workflow's volcano
#   plots; the reason is given above plot_ma_contrast().
#
# THREE PLACES THE STANDARD WORKFLOW COULD NOT BE FOLLOWED, AND WHY
#
# (a) limma-voom is replaced by edgeR with a fixed dispersion.
#     voom needs residual degrees of freedom to fit the mean-variance trend. The
#     design here is ~0 + sample with one library per sample, which is saturated
#     and leaves zero residual df, so voom has nothing to estimate from and its
#     weights would be meaningless. edgeR with a SUPPLIED dispersion is the
#     documented route for unreplicated designs. Everything downstream is
#     unchanged: topTags output is renamed to limma's column names (logFC,
#     AveExpr, t, P.Value, adj.P.Val) so the section 8 and 9 code from the
#     standard workflow works on it directly.
#
# (b) No patient blocking, and no batch term.
#     The standard workflow blocks on patient because every contrast there is
#     measured within patient. Here there is one patient and one donor, so there is nothing to
#     block on; the contrast structure itself carries the design.
#
# (c) Per-comparison MDS is over all samples, not the contrast pair.
#     An MDS of the two libraries in a contrast is two points and conveys
#     nothing. The per-comparison plots instead highlight the contrast pair
#     within the full four-sample space, which is the closest informative
#     analogue.
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma)
  library(RColorBrewer); library(pheatmap); library(ggplot2)
  library(clusterProfiler); library(org.Hs.eg.db)
  library(enrichplot); library(stringr); library(ggnewscale); library(viridis)
})

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

fig_height <- function(n_rows, per = 0.34, base = 2, min_h = 6)
  max(min_h, n_rows * per + base)

# sample annotation used by every plot below
ann_df <- function() {
  a <- SAMPLES[SAMPLE_ORDER]
  data.frame(
    timepoint = a$short,
    subject   = a$subject,
    disease   = a$disease,
    clone     = a$clone_ann,
    row.names = a$sample)
}

ann_colours <- function() list(
  timepoint = setNames(PAL$sample[SAMPLE_ORDER], SAMPLES[SAMPLE_ORDER, short]),
  subject   = c(patient = "#C44E52", donor = "#55A868"),
  disease   = c(no = "#D9D9D9", DDHN = "#C44E52"),
  # "not tested" is a mid grey, not a near-white: it also colours a point and its
  # label in the clone MDS plot, where a pale tint vanishes against the background
  clone     = CLONE_COLOURS)

################################################################################
# Section 1-2 : library sizes, normalisation diagnostics
################################################################################
plot_library_sizes <- function(dge, file) {
  png(file, width = 3000, height = 2000, units = "px", pointsize = 12, res = 300)
  par(mar = c(7, 5, 4, 2))
  barplot(dge$samples$lib.size / 1e6, names.arg = colnames(dge), las = 2,
          col = PAL$sample[colnames(dge)], ylab = "Library size (millions)")
  title("Barplot of library sizes (globin and MT removed)")
  dev.off()
}

plot_expression_boxplots <- function(log_raw, log_norm, file) {
  png(file, width = 6000, height = 3000, units = "px", pointsize = 12, res = 450)
  par(mfrow = c(1, 2), mar = c(7, 5, 4, 2))
  boxplot(log_raw, xlab = "", ylab = "Log2 CPM", las = 2,
          main = "Unnormalised logCPM", col = PAL$sample[colnames(log_raw)])
  abline(h = median(log_raw), col = "blue")
  boxplot(log_norm, xlab = "", ylab = "Log2 CPM", las = 2,
          main = "TMM-normalised logCPM", col = PAL$sample[colnames(log_norm)])
  abline(h = median(log_norm), col = "blue")
  dev.off()
}

# Replaces the standard workflow's voom mean-variance plot. voom cannot run on a
# saturated design (see header note a), so the analogous diagnostic is the
# mean-variance relationship of the counts with the assumed BCV overlaid.
plot_meanvar <- function(dge, file, bcv = THRESH$bcv_primary) {
  cpm_m <- cpm(dge, log = TRUE, prior.count = 3)
  mu <- rowMeans(cpm_m); sdv <- apply(cpm_m, 1, sd)
  png(file, width = 2500, height = 2000, units = "px", pointsize = 10, res = 300)
  par(mar = c(5, 5, 4, 2))
  smoothScatter(mu, sdv, xlab = "Average log2 CPM", ylab = "SD of log2 CPM",
                main = sprintf("Mean-variance relationship (assumed BCV = %.1f)", bcv))
  lo <- try(lowess(mu, sdv, f = .3), silent = TRUE)
  if (!inherits(lo, "try-error")) lines(lo, col = "red", lwd = 2)
  abline(h = bcv / log(2), col = "blue", lty = 2, lwd = 2)
  legend("topright", legend = c("lowess trend", sprintf("assumed BCV %.1f", bcv)),
         col = c("red", "blue"), lty = c(1, 2), lwd = 2, bty = "n")
  dev.off()
}

################################################################################
# Section 3-4 : MDS / PCA
################################################################################
mds_xy <- function(logcpm, dim = c(1, 2)) {
  m <- plotMDS(logcpm, plot = FALSE, dim.plot = dim)
  data.table(sample = colnames(logcpm), x = as.numeric(m$x), y = as.numeric(m$y),
             xlab = m$axislabel, dim = paste(dim, collapse = "-"))
}

.mds_plot <- function(d, colour_by, title, file, label = TRUE) {
  # "timepoint" is the display name; the column in SAMPLES is `short`. "clone"
  # uses the short legend form (clone_ann), which is what the palette is keyed by
  col <- c(timepoint = "short", subject = "subject", disease = "disease",
           clone = "clone_ann")[colour_by]
  if (is.na(col)) col <- colour_by
  a <- SAMPLES[d$sample]
  d[, grp := a[[col]]]
  d[, lab := a$short]
  pal <- switch(colour_by,
    subject   = ann_colours()$subject,
    disease   = ann_colours()$disease,
    clone     = ann_colours()$clone,
    timepoint = ann_colours()$timepoint,
    setNames(PAL$sample[d$sample], d$grp))
  p <- ggplot(d, aes(x, y, colour = grp)) +
    geom_hline(yintercept = 0, colour = "grey88") +
    geom_vline(xintercept = 0, colour = "grey88") +
    geom_point(size = 5) +
    { if (label) ggrepel::geom_text_repel(aes(label = lab), size = 3.6,
                                          show.legend = FALSE) } +
    scale_colour_manual(values = pal) +
    labs(x = "Leading logFC dim 1", y = "Leading logFC dim 2",
         colour = colour_by, title = title) +
    theme_bw(base_size = 11)
  ggsave(file, p, width = 6.5, height = 5.5, dpi = 450)
}

# Section 3: MDS over all samples, one plot per annotation
plot_mds_all <- function(logcpm, prefix) {
  d <- mds_xy(logcpm)
  .mds_plot(copy(d), "timepoint", "MDS, all samples (by timepoint)",
            paste0(prefix, "12_dge_mds_timepoint.png"))
  .mds_plot(copy(d), "subject",  "MDS, all samples (by subject)",
            paste0(prefix, "12_dge_mds_subject.png"))
  .mds_plot(copy(d), "disease",  "MDS, all samples (by disease state)",
            paste0(prefix, "12_dge_mds_disease.png"))
  .mds_plot(copy(d), "clone",    "MDS, all samples (by DNMT3A clone)",
            paste0(prefix, "12_dge_mds_clone.png"))

  # dimensions 2-3: with four libraries dim 1 is dominated by the donor, so the
  # structure among the patient timepoints only becomes visible further out
  d23 <- mds_xy(logcpm, dim = c(2, 3))
  .mds_plot(copy(d23), "timepoint", "MDS dimensions 2-3, all samples",
            paste0(prefix, "12_dge_mds_dim23.png"))

  # patient-only MDS: the donor's compartment difference otherwise swamps
  # everything and the timecourse cannot be seen at all
  pat <- logcpm[, PATIENT_SAMPLES]
  dp <- mds_xy(pat)
  .mds_plot(copy(dp), "timepoint", "MDS, patient timepoints only",
            paste0(prefix, "12_dge_mds_patientOnly.png"))
  invisible(d)
}

# Section 4: one plot per contrast, highlighting that contrast's two libraries
# within the full sample space
plot_mds_per_contrast <- function(logcpm, prefix) {
  d0 <- mds_xy(logcpm)
  for (i in seq_len(nrow(CONTRASTS))) {
    cc <- CONTRASTS[i]
    d <- copy(d0)
    d[, role := fifelse(sample == cc$num, "numerator",
                fifelse(sample == cc$den, "denominator", "other"))]
    d[, lab := SAMPLES[d$sample, short]]
    seg <- d[role != "other"]
    p <- ggplot(d, aes(x, y)) +
      geom_hline(yintercept = 0, colour = "grey88") +
      geom_vline(xintercept = 0, colour = "grey88") +
      { if (nrow(seg) == 2)
          annotate("segment", x = seg$x[1], y = seg$y[1],
                   xend = seg$x[2], yend = seg$y[2],
                   colour = "grey45", linetype = 2) } +
      geom_point(aes(colour = role, size = role)) +
      ggrepel::geom_text_repel(aes(label = lab), size = 3.6) +
      scale_colour_manual(values = c(numerator = PAL$up, denominator = PAL$down,
                                     other = "grey78")) +
      scale_size_manual(values = c(numerator = 5.5, denominator = 5.5, other = 3),
                        guide = "none") +
      labs(x = "Leading logFC dim 1", y = "Leading logFC dim 2",
           title = cc$label,
           subtitle = paste0(cc$class, " contrast")) +
      theme_bw(base_size = 11)
    ggsave(sprintf("%s13_mds_%s.png", prefix, cc$id), p,
           width = 6.5, height = 5.5, dpi = 450)
  }
}

################################################################################
# Section 5 : heatmaps of the most variable features
################################################################################
plot_heatmap <- function(logmat, n, file, w_in, h_in, fontsize_row = 6,
                         order_by_timepoint = FALSE, label = "features") {
  vg  <- apply(logmat, 1, var)
  sel <- names(sort(vg, decreasing = TRUE))[seq_len(min(n, nrow(logmat)))]
  hv  <- logmat[sel, , drop = FALSE]

  ac <- ann_df()
  cluster_cols <- TRUE
  if (order_by_timepoint) {
    ord <- match(SAMPLE_ORDER, colnames(hv))
    hv <- hv[, ord, drop = FALSE]
    cluster_cols <- FALSE
  }
  pheatmap(hv,
           scale             = "row",
           color             = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
           annotation_col    = ac,
           annotation_colors = ann_colours(),
           cluster_cols      = cluster_cols,
           show_rownames     = TRUE,
           show_colnames     = TRUE,
           fontsize_row      = fontsize_row,
           main              = paste("Top", n, "most variable", label),
           filename          = file, width = w_in, height = h_in)
}

plot_all_heatmaps <- function(logcpm, prefix, label = "genes") {
  plot_heatmap(logcpm,  50, paste0(prefix, "12_dge_heatmap_top50.png"),
               8, 12, fontsize_row = 7, label = label)
  plot_heatmap(logcpm, 100, paste0(prefix, "12_dge_heatmap_top100.png"),
               8, 18, fontsize_row = 5, label = label)
  plot_heatmap(logcpm,  50, paste0(prefix, "12_dge_heatmap_top50_byTimepoint.png"),
               8, 12, fontsize_row = 7, order_by_timepoint = TRUE, label = label)
  plot_heatmap(logcpm, 100, paste0(prefix, "12_dge_heatmap_top100_byTimepoint.png"),
               8, 18, fontsize_row = 5, order_by_timepoint = TRUE, label = label)
}

# heatmap of the DE features for one contrast
plot_contrast_heatmap <- function(logcpm, genes, cid, file, n = 60) {
  genes <- intersect(genes, rownames(logcpm))
  if (length(genes) < 3) return(invisible())
  genes <- head(genes, n)
  pheatmap(logcpm[genes, , drop = FALSE],
           scale = "row",
           color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
           annotation_col = ann_df(), annotation_colors = ann_colours(),
           cluster_cols = FALSE, show_rownames = TRUE,
           fontsize_row = max(4, min(8, 420 / length(genes))),
           main = paste0(CONTRASTS[cid, label], "\ntop ", length(genes),
                         " differential features"),
           filename = file, width = 7,
           height = max(4, length(genes) * 0.16 + 2))
}

################################################################################
# Section 6 : DGE. edgeR with fixed dispersion, output renamed to limma columns
#   so sections 8 and 9 read exactly as they did in the standard workflow.
################################################################################
fit_contrasts <- function(dge, bcv = THRESH$bcv_primary) {
  grp    <- factor(colnames(dge), levels = SAMPLE_ORDER)
  design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
  cm <- makeContrasts(contrasts = CONTRASTS[, sprintf("%s-%s", num, den)],
                      levels = design)
  colnames(cm) <- CONTRASTS$id
  list(fit = glmFit(dge, design, dispersion = bcv^2), design = design, cm = cm)
}

# topTable-equivalent. Column names deliberately match limma so that the ORA and
# GSEA code carried over from the standard workflow needs no modification.
top_table <- function(f, cid, lcpm1, cpm_mat) {
  lrt <- glmLRT(f$fit, contrast = f$cm[, cid])
  tt  <- as.data.table(topTags(lrt, n = Inf, sort.by = "none")$table,
                       keep.rownames = "gene")
  cc  <- CONTRASTS[cid]
  tt[, `:=`(
    AveExpr    = logCPM,
    t          = sign(lfc(lcpm1, gene, cc)) * sqrt(pmax(LR, 0)),
    P.Value    = PValue,
    adj.P.Val  = FDR,
    lfc_shrunk = lfc(lcpm1, gene, cc),
    cpm_num    = cpm_mat[cbind(gene, cc$num)],
    cpm_den    = cpm_mat[cbind(gene, cc$den)])]
  tt[, `:=`(max_cpm_pair = pmax(cpm_num, cpm_den),
            min_cpm_pair = pmin(cpm_num, cpm_den))]
  tt[, on_off := min_cpm_pair < 1 & max_cpm_pair >= 1]
  tt[, sig := adj.P.Val < THRESH$fdr & abs(lfc_shrunk) >= THRESH$lfc &
              max_cpm_pair >= 1]
  tt[, contrast := cid]
  tt[]
}
lfc <- function(lcpm1, gene, cc) lcpm1[cbind(gene, cc$num)] - lcpm1[cbind(gene, cc$den)]

################################################################################
# Section 7a : MD plots
################################################################################
plot_md_contrast <- function(tt, cid, file) {
  png(file, width = 3000, height = 3000, units = "px", pointsize = 14, res = 450)
  par(mar = c(5, 5, 4, 2))
  st <- ifelse(tt$sig & tt$lfc_shrunk > 0, 1, ifelse(tt$sig & tt$lfc_shrunk < 0, -1, 0))
  plot(tt$AveExpr, tt$lfc_shrunk, pch = 16, cex = .35, col = "grey72",
       xlab = "Average log2 CPM", ylab = "log2 fold-change",
       main = paste("MD:", CONTRASTS[cid, label]))
  points(tt$AveExpr[st ==  1], tt$lfc_shrunk[st ==  1], pch = 16, cex = .45, col = PAL$up)
  points(tt$AveExpr[st == -1], tt$lfc_shrunk[st == -1], pch = 16, cex = .45, col = PAL$down)
  abline(h = 0, col = "blue"); abline(h = c(-THRESH$lfc, THRESH$lfc), lty = 2, col = "grey40")
  legend("topright", bty = "n", pch = 16, col = c(PAL$up, PAL$down),
         legend = c(paste("up", sum(st == 1)), paste("down", sum(st == -1))))
  dev.off()
}

################################################################################
# Section 7b : per-contrast MA plot (replaces the standard workflow's volcano)
################################################################################
#
# WHY NOT A VOLCANO
#   A volcano plot works because effect size and statistical evidence are two
#   different things: a gene can move a lot unconvincingly, or a little very
#   convincingly, and the plot separates those. That separation requires an
#   estimate of variability.
#
#   Here there is none. With one library per group and a supplied dispersion the
#   p-value is very nearly a deterministic function of the fold-change and the
#   expression level. Measured on this data, among genes with more than 20 counts
#   in both libraries of a pair, the correlation between |log2FC| and -log10(P)
#   is 0.884 in the transformation contrast (14,849 genes) and 0.911-0.918 in the
#   other five. The two volcano axes are therefore close to the same axis, the
#   points collapse onto a manifold, and the characteristic arcs appear. That is
#   not a rendering fault to be recoloured -- the figure is showing a real
#   property of the statistic, and the property is that a volcano cannot say
#   anything here that the fold-change does not already say.
#
#   An MA plot avoids this because its axes are two genuinely independent
#   quantities: how strongly a gene is expressed, and how much it changed.
#
# GENES ABSENT FROM ONE LIBRARY ARE NOT PLOTTED
#   For a gene with zero counts in one library the fold-change is undefined, and
#   the value we would plot is set by the pseudo-count rather than by the data:
#   PF4 (2,540 vs 0 counts) sits at 4.46, 5.43 or 2.38 depending only on whether
#   the pseudo-count is 1, 0.5 or 5. Plotting that as a measured effect size is
#   not defensible, so these genes are excluded from the figure and reported in
#   their own table with the raw counts and no fold-change.
plot_ma_contrast <- function(tt, cid, prefix, n_label = 25) {
  d <- as.data.table(copy(tt))
  d[, ave := (log2(cpm_num + 1) + log2(cpm_den + 1)) / 2]
  d[, zero_in_pair := pmin(cpm_num, cpm_den) == 0]

  plotd <- d[zero_in_pair == FALSE]
  plotd[, state := fifelse(robust, "kept: passed all filters",
                  fifelse(sig | noiseq_sig, "called by one approach only",
                          "not called"))]
  plotd[, state := factor(state, levels = c("kept: passed all filters",
                                            "called by one approach only",
                                            "not called"))]
  pal <- c("kept: passed all filters"    = "#C0392B",
           "called by one approach only" = "#E8B647",
           "not called"                  = "#C8CDD2")
  lab <- plotd[robust == TRUE][order(adj.P.Val)][seq_len(min(n_label, .N))]

  p <- ggplot(plotd[order(state, decreasing = TRUE)], aes(ave, lfc_shrunk)) +
    geom_hline(yintercept = 0, colour = "grey35") +
    geom_hline(yintercept = c(-THRESH$lfc, THRESH$lfc), linetype = 2,
               colour = "grey55") +
    geom_point(aes(colour = state), size = 1.05, alpha = .6) +
    scale_colour_manual(values = pal, drop = FALSE, name = NULL) +
    { if (nrow(lab))
        ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3,
                                 max.overlaps = 22, min.segment.length = 0,
                                 colour = "grey15") } +
    labs(x = expression("Average expression, mean "*log[2]*"(CPM + 1) of the two libraries"),
         y = expression(Log[2]~"fold change  (numerator "-" denominator)"),
         title = CONTRASTS[cid, label],
         subtitle = sprintf(
           "%s genes kept; %s excluded (absent from one library, so no defined fold-change); dashed lines mark a %g-fold change",
           format(sum(d$robust), big.mark = ","),
           format(sum(d$zero_in_pair), big.mark = ","), 2^THRESH$lfc)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8.5),
          plot.subtitle = element_text(size = 8.4, colour = "grey30")) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

  ggsave(paste0(prefix, "22_MAplot_", cid, ".pdf"), p, width = 8.4, height = 6.6)
  ggsave(paste0(prefix, "22_MAplot_", cid, ".png"), p, width = 8.4, height = 6.6,
         dpi = 300)

  # the excluded genes, reported honestly with counts rather than a fold-change
  cc <- CONTRASTS[cid]
  z <- d[zero_in_pair == TRUE][order(-pmax(cpm_num, cpm_den))]
  if (nrow(z))
    fwrite(z[, .(gene,
                 present_in = fifelse(cpm_num > 0, cc$num, cc$den),
                 cpm_present = round(pmax(cpm_num, cpm_den), 2),
                 cpm_absent = 0,
                 P.Value, adj.P.Val)],
           paste0(prefix, "24_absent_from_one_library_", cid, ".tsv"), sep = "\t")
  invisible(nrow(z))
}

################################################################################
# Section 8 helpers : ORA
################################################################################
prep_contrast <- function(tt) {
  d <- copy(tt)
  d[, symbol := gene]
  suppressMessages(
    d[, entrez := mapIds(org.Hs.eg.db, symbol, "ENTREZID", "SYMBOL")])
  as.data.frame(d)
}

enrich_one <- function(genes, universe, term_p = 0.05) {
  message("    genes for ORA: ", length(genes))
  if (length(genes) < 5) return(NULL)
  go <- function(ont)
    tryCatch(enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
                      universe = universe, pAdjustMethod = "BH",
                      pvalueCutoff = term_p, qvalueCutoff = 0.2, readable = TRUE),
             error = function(e) { message("    GO ", ont, " failed: ", e$message); NULL })
  kegg <- tryCatch({
    k <- enrichKEGG(genes, organism = "hsa", universe = universe,
                    pAdjustMethod = "BH", pvalueCutoff = term_p)
    if (!is.null(k)) setReadable(k, org.Hs.eg.db, "ENTREZID") else NULL
  }, error = function(e) { message("    KEGG failed: ", e$message); NULL })
  hm <- if (exists("HALLMARK_T2G")) tryCatch({
    h <- enricher(genes, TERM2GENE = HALLMARK_T2G, universe = universe,
                  pAdjustMethod = "BH", pvalueCutoff = term_p, qvalueCutoff = 0.2)
    if (!is.null(h)) setReadable(h, org.Hs.eg.db, "ENTREZID") else NULL
  }, error = function(e) { message("    Hallmark failed: ", e$message); NULL }) else NULL
  list(BP = go("BP"), CC = go("CC"), MF = go("MF"), KEGG = kegg, HALLMARK = hm)
}

save_dotplot <- function(er, title, file, top_terms = 20) {
  if (is.null(er)) return(invisible())
  df <- as.data.frame(er); if (!nrow(df)) return(invisible())
  n_rows <- min(top_terms, nrow(df))
  h <- fig_height(n_rows, per = 0.40)
  p <- dotplot(er, showCategory = top_terms) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 55)) +
    ggtitle(title) + theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 8, lineheight = 0.85),
          plot.margin = margin(6, 6, 6, 10))
  ggsave(paste0(file, ".pdf"), p, width = 11, height = h, limitsize = FALSE)
  ggsave(paste0(file, ".png"), p, width = 11, height = h, dpi = 300, limitsize = FALSE)
}

save_cnetplot <- function(er, fc, title, file) {
  if (is.null(er) || !nrow(as.data.frame(er))) return(invisible())
  p <- tryCatch(cnetplot(er, showCategory = 7, foldChange = fc,
                         node_label = "category") + ggtitle(title),
                error = function(e)
                  tryCatch(cnetplot(er, showCategory = 7,
                                    color.params = list(foldChange = fc),
                                    node_label = "category") + ggtitle(title),
                           error = function(e2) NULL))
  if (is.null(p)) return(invisible())
  ggsave(paste0(file, ".pdf"), p, width = 11, height = 9)
  ggsave(paste0(file, ".png"), p, width = 11, height = 9, dpi = 400)
}

################################################################################
# Sections 8-9 : matched MA plot for each enriched term
#   The contrast's own MA plot, redrawn with the term's genes (ORA) or leading
#   edge (GSEA) highlighted, so the genes driving an enrichment and where they
#   sit can be seen directly. Replaces the standard workflow's matched volcano,
#   for the reason given above plot_ma_contrast(): with a supplied dispersion a
#   volcano's two axes are nearly the same axis. Drawn from the saved CSVs by
#   04b_enrichment_ma_plots.R. `tt` is one contrast's DGE table (needs gene,
#   lfc_shrunk, cpm_num, cpm_den, robust).
################################################################################
term_ma <- function(tt, term_genes, title, file, gene_label = "term genes",
                    n_label = 30) {
  d <- as.data.table(copy(tt))
  d[, in_term := gene %in% term_genes]
  n_absent <- d[in_term == TRUE & pmin(cpm_num, cpm_den) == 0, .N]
  d <- d[pmin(cpm_num, cpm_den) > 0]          # fold-change undefined otherwise
  if (!any(d$in_term)) return(invisible(FALSE))
  d[, ave := (log2(cpm_num + 1) + log2(cpm_den + 1)) / 2]
  up_lab <- "term gene, higher in first-named group"
  dn_lab <- "term gene, lower in first-named group"
  tg  <- d[in_term == TRUE]
  tg[, grp := factor(fifelse(lfc_shrunk > 0, up_lab, dn_lab), levels = c(up_lab, dn_lab))]
  lab <- tg[order(-abs(lfc_shrunk))][seq_len(min(n_label, .N))]
  sub <- paste0(nrow(tg), " ", gene_label, " plotted; black ring = passed all filters",
                if (n_absent) paste0("; ", n_absent, " absent from one library, not plotted") else "",
                if (nrow(tg) > n_label) paste0("; the ", n_label, " largest changes labelled") else "")
  p <- ggplot(d, aes(ave, lfc_shrunk)) +
    geom_hline(yintercept = 0, colour = "grey35") +
    geom_hline(yintercept = c(-THRESH$lfc, THRESH$lfc), linetype = 2, colour = "grey60") +
    geom_point(data = d[in_term == FALSE], colour = "#D5D9DD", size = .8, alpha = .5) +
    geom_point(data = tg, aes(colour = grp), size = 1.9, alpha = .9) +
    geom_point(data = tg[robust == TRUE], shape = 1, size = 3.3, stroke = .7, colour = "black") +
    scale_colour_manual(values = setNames(c("#C0392B", "#2E6DB4"), c(up_lab, dn_lab)),
                        drop = TRUE, name = NULL) +
    ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.8, max.overlaps = 30,
                             min.segment.length = 0, colour = "grey15") +
    labs(x = expression("Average expression, mean "*log[2]*"(CPM + 1) of the two libraries"),
         y = expression(Log[2]~"fold change  (first-named "-" second-named group)"),
         title = str_wrap(title, 85), subtitle = str_wrap(sub, 120)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10.5),
          plot.subtitle = element_text(size = 8.2, colour = "grey30"))
  ggsave(paste0(file, ".pdf"), p, width = 8.4, height = 7)
  ggsave(paste0(file, ".png"), p, width = 8.4, height = 7, dpi = 300)
  invisible(TRUE)
}

# one matched MA plot per ORA term (top `top_terms` in file order, i.e. by p)
category_ma_plots <- function(tt, df, category, outdir, top_terms = 20) {
  if (is.null(df) || !nrow(df)) return(0L)
  df <- head(as.data.frame(df), top_terms)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  n <- 0L
  for (j in seq_len(nrow(df))) {
    genes <- strsplit(df$geneID[j], "/")[[1]]
    ok <- term_ma(tt, genes, paste0(category, ": ", df$Description[j]),
                  file.path(outdir, sprintf("%02d_%s", j, safe_name(df$Description[j]))))
    n <- n + isTRUE(ok)
  }
  n
}

################################################################################
# Section 9 helpers : GSEA
################################################################################
rank_genes <- function(tt) {
  d <- tt[!is.na(tt$entrez) & is.finite(tt$t), ]
  d <- d[order(abs(d$t), decreasing = TRUE), ]
  d <- d[!duplicated(d$entrez), ]
  r <- d$t; names(r) <- d$entrez
  sort(r, decreasing = TRUE)
}

# GSEA for one ranked list, all five collections.
# fgsea samples, and clusterProfiler's seed = TRUE does NOT make repeated calls
# reproducible: two back-to-back runs on identical input differed by up to 0.05
# in NES. set.seed() immediately before every call, with BiocParallel registered
# as SerialParam by the calling script, makes results identical run to run.
# Use gsea_term_p = 1 to get every set tested (needed by the reference-contrast
# filter, which reads raw p for sets that did not reach 0.05).
gsea_one <- function(ranks, gsea_term_p = 0.05, seed = 1234) {
  go <- function(ont) {
    set.seed(seed)
    tryCatch(setReadable(
      gseGO(ranks, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
            pAdjustMethod = "BH", pvalueCutoff = gsea_term_p, eps = 0,
            seed = TRUE, verbose = FALSE), org.Hs.eg.db, "ENTREZID"),
      error = function(e) { message("    gseGO ", ont, " failed: ", e$message); NULL })
  }
  set.seed(seed)
  kegg <- tryCatch(setReadable(
    gseKEGG(ranks, organism = "hsa", pAdjustMethod = "BH",
            pvalueCutoff = gsea_term_p, eps = 0, seed = TRUE, verbose = FALSE),
    org.Hs.eg.db, "ENTREZID"),
    error = function(e) { message("    gseKEGG failed: ", e$message); NULL })
  hm <- NULL
  if (exists("HALLMARK_T2G")) {
    set.seed(seed)
    hm <- tryCatch(setReadable(
      GSEA(ranks, TERM2GENE = HALLMARK_T2G, pAdjustMethod = "BH",
           pvalueCutoff = gsea_term_p, eps = 0, seed = TRUE, verbose = FALSE),
      org.Hs.eg.db, "ENTREZID"),
      error = function(e) { message("    Hallmark GSEA failed: ", e$message); NULL })
  }
  list(BP = go("BP"), CC = go("CC"), MF = go("MF"), KEGG = kegg, HALLMARK = hm)
}

save_gsea_dotplot <- function(gs, title, file, gsea_top_dot = 15) {
  if (is.null(gs)) return(invisible())
  df <- as.data.frame(gs); if (!nrow(df)) return(invisible())
  n_rows <- min(gsea_top_dot, sum(df$NES > 0)) + min(gsea_top_dot, sum(df$NES < 0))
  h <- fig_height(n_rows, per = 0.40)
  p <- tryCatch(
    dotplot(gs, showCategory = gsea_top_dot, split = ".sign") +
      facet_grid(. ~ .sign) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 55)) +
      ggtitle(title) + theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"),
            axis.text.y = element_text(size = 8, lineheight = 0.85),
            panel.spacing = unit(1, "lines"),
            plot.margin = margin(6, 6, 6, 10)),
    error = function(e) NULL)
  if (is.null(p)) return(invisible())
  ggsave(paste0(file, ".pdf"), p, width = 13, height = h, limitsize = FALSE)
  ggsave(paste0(file, ".png"), p, width = 13, height = h, dpi = 300, limitsize = FALSE)
}

save_ridgeplot <- function(gs, title, file, gsea_top_dot = 15) {
  if (is.null(gs)) return(invisible())
  df <- as.data.frame(gs); if (!nrow(df)) return(invisible())
  h <- fig_height(min(gsea_top_dot, nrow(df)), per = 0.42)
  p <- tryCatch(
    ridgeplot(gs, showCategory = gsea_top_dot) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 55)) +
      ggtitle(title) + theme_bw(base_size = 10) +
      theme(axis.text.y = element_text(size = 8, lineheight = 0.85),
            plot.margin = margin(6, 6, 6, 10)),
    error = function(e) NULL)
  if (is.null(p)) return(invisible())
  ggsave(paste0(file, ".pdf"), p, width = 11, height = h, limitsize = FALSE)
  ggsave(paste0(file, ".png"), p, width = 11, height = h, dpi = 300, limitsize = FALSE)
}

save_gsea_curves <- function(gs, category, outdir, gsea_top_curves = 10) {
  if (is.null(gs)) return(invisible())
  df <- as.data.frame(gs); if (!nrow(df)) return(invisible())
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  o    <- order(df$p.adjust)
  up   <- head(df[o[df$NES[o] > 0], , drop = FALSE], gsea_top_curves)
  down <- head(df[o[df$NES[o] < 0], , drop = FALSE], gsea_top_curves)
  pick <- rbind(up, down)
  for (j in seq_len(nrow(pick))) {
    ttl <- sprintf("%s: %s (NES=%.2f, p.adj=%.1e)", category,
                   pick$Description[j], pick$NES[j], pick$p.adjust[j])
    p <- tryCatch(gseaplot2(gs, geneSetID = pick$ID[j], title = ttl,
                            pvalue_table = FALSE), error = function(e) NULL)
    if (is.null(p)) next
    fn <- file.path(outdir, sprintf("%02d_%s", j, safe_name(pick$Description[j])))
    ggsave(paste0(fn, ".pdf"), p, width = 8, height = 6)
    ggsave(paste0(fn, ".png"), p, width = 8, height = 6, dpi = 350)
  }
}

# one matched MA plot per GSEA set: top `top` per direction by adjusted p,
# highlighting the leading edge (the genes actually generating the score)
gsea_term_ma_plots <- function(tt, df, category, outdir, top = 15) {
  if (is.null(df) || !nrow(df)) return(0L)
  df <- as.data.frame(df); df <- df[!is.na(df$NES), ]
  if (!nrow(df)) return(0L)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  pick_dir <- function(pos) {
    sub <- df[if (pos) df$NES > 0 else df$NES < 0, , drop = FALSE]
    if (!nrow(sub)) return(sub)
    sub[order(sub$p.adjust), , drop = FALSE][seq_len(min(top, nrow(sub))), , drop = FALSE]
  }
  act <- pick_dir(TRUE); sup <- pick_dir(FALSE)
  pick <- rbind(act, sup); if (!nrow(pick)) return(0L)
  rank_in_dir <- c(seq_len(nrow(act)), seq_len(nrow(sup)))
  n <- 0L
  for (j in seq_len(nrow(pick))) {
    core <- strsplit(pick$core_enrichment[j], "/")[[1]]
    dirn <- if (pick$NES[j] > 0) "activated" else "suppressed"
    ttl  <- sprintf("%s: %s (NES = %.2f, %s, FDR = %.2g)", category, pick$Description[j],
                    pick$NES[j], dirn, pick$p.adjust[j])
    fn   <- file.path(outdir, sprintf("%s_%02d_%s", dirn, rank_in_dir[j],
                                      safe_name(pick$Description[j])))
    ok <- term_ma(tt, core, ttl, fn, gene_label = "leading-edge genes")
    n <- n + isTRUE(ok)
  }
  n
}
