################################################################################
# 02_dge_gene_level.R  --  gene-level differential expression
#
# Run:  Rscript 03_rnaseq/02_dge_gene_level.R
# Out:  05_processed_counts_genes/
#
# Enrichment (sections 8 and 9) lives in 04_enrichment.R so it can be re-run
# without repeating the DGE.
#
# Design note. There is one patient and one donor, so there is
# nothing to block on and no batch term is possible; the contrast structure
# carries the design instead. limma-voom is replaced by edgeR with a fixed
# dispersion because a ~0 + sample design with one library per sample is
# saturated and leaves voom no residual degrees of freedom to work with -- see
# the header of dge_functions.R.
#
# Contrasts (direction = numerator - denominator; positive = higher in first):
#   C1_clone          3.5yr vs D30   clonal establishment: the clone expands, no disease either end
#   C2_transform      DDHN vs 3.5yr  transformation: clone size flat while disease appears
#   C3_cumulative     DDHN vs D30    both axes summed
#   C4_DDHN_vs_donor  DDHN vs donor  compartment-confounded, context only
#   C5_2yr_vs_donor   3.5yr vs donor compartment-confounded, context only
#   C6_D30_vs_donor   D30  vs donor  RULER: no disease, small clone at both ends -> artefact floor
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma); library(Glimma)
  library(org.Hs.eg.db); library(RColorBrewer); library(pheatmap)
  library(ggplot2); library(viridis); library(NOISeq); library(GSVA)
})

SCRIPTS <- "03_rnaseq"
source(file.path(SCRIPTS, "00_config.R"))
source(file.path(SCRIPTS, "01_load_counts.R"))
source(file.path(SCRIPTS, "dge_functions.R"))

OUT <- PATHS$out_gene
fig <- function(...) file.path(OUT, paste0(...))
PFX <- paste0(OUT, "/")
set.seed(1234)

################################################################################
# 0. Load counts, collapse to gene level, build annotation
################################################################################
message("\n== 0. loading gene counts ==")
g <- load_gene_counts()

sym <- g$anno$gene_name
sym[is.na(sym) | sym == ""] <- g$anno$Geneid[is.na(sym) | sym == ""]

# Globin and mitochondrial transcripts removed BEFORE normalisation: globin runs
# 0.07% -> 44.4% and MT 2.1% -> 28.1% across these four libraries, so leaving them
# in would make TMM correct for erythrocyte content and RNA integrity rather
# than for sequencing depth. Both fractions are reported, not silently dropped.
strip <- strip_confounders(g$counts, sym)
counts_all <- strip$counts
anno_all   <- g$anno[strip$keep]
sym_all    <- sym[strip$keep]
fwrite(strip$fractions, fig("00_library_composition.tsv"), sep = "\t")

################################################################################
# 0b. Case-specific QC: composition, lineage, confound signatures
#     Kept ahead of the standard pipeline because in this dataset composition is
#     the leading alternative explanation for every gene-level result.
################################################################################
message("\n== 0b. composition / confound QC ==")
comp <- melt(strip$fractions[, .(sample, globin_pct, mt_pct, retained_pct)],
             id.vars = "sample", variable.name = "fraction", value.name = "pct")
comp[, fraction := factor(fraction, levels = c("retained_pct", "mt_pct", "globin_pct"),
                          labels = c("retained", "mitochondrial", "globin"))]
comp[, sample := factor(sample, levels = SAMPLE_ORDER)]
ggsave(fig("03_library_composition.png"),
  ggplot(comp, aes(sample, pct, fill = fraction)) + geom_col(width = .7) +
    scale_fill_manual(values = c(retained = "#9EC1CF", mitochondrial = "#DD8452",
                                 globin = "#C44E52")) +
    labs(y = "% of raw counts", x = NULL, title = "Library composition before filtering",
         subtitle = "Both gradients run along the D30->DDHN axis") +
    theme_bw(base_size = 11), width = 7, height = 4.5, dpi = 300)

################################################################################
# 1. Filtering (CPM > 1 in >= 2 samples)
################################################################################
message("\n== 1. filtering ==")
keep <- filter_expressed(counts_all)
counts <- counts_all[keep, , drop = FALSE]
anno   <- anno_all[keep]
gsym   <- sym_all[keep]
rownames(counts) <- make.unique(gsym)
check_sample_order(counts)

################################################################################
# 2. DGEList + TMM normalisation
################################################################################
message("\n== 2. DGEList + TMM ==")
dge <- DGEList(counts = counts,
               genes = data.frame(symbol = rownames(counts), ensembl = anno$Geneid))
dge <- normLibSizes(dge, method = "TMM")

log_raw_dgeLlist <- cpm(DGEList(counts), log = TRUE, prior.count = 3)
logcpm  <- cpm(dge, log = TRUE, prior.count = 3)
cpm_mat <- cpm(dge)
lcpm1   <- log2(cpm_mat + 1)   # bounded effect-size scale, see section 6

plot_library_sizes(dge, fig("11_dge_library_sizes.png"))
plot_expression_boxplots(log_raw_dgeLlist, logcpm, fig("11_dge_expression_boxplots.png"))
plot_meanvar(dge, fig("11_dge_mean_variance.png"))

norm_tab <- data.table(sample = colnames(dge), lib_size = dge$samples$lib.size,
                       norm_factor = round(dge$samples$norm.factors, 4))
print(norm_tab); fwrite(norm_tab, fig("01_normalisation_factors.tsv"), sep = "\t")
write.table(round(logcpm, 4), fig("rnaSeq_normalized_counts_collated.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# lineage + confound signature scoring (case-specific, see 0b)
lin_sets <- lapply(LINEAGE_MARKERS, function(x) intersect(x, rownames(logcpm)))
lin_sets <- lin_sets[sapply(lin_sets, length) >= 3]
lin_dt <- as.data.table(GSVA::gsva(GSVA::ssgseaParam(logcpm, lin_sets, normalize = FALSE),
                                   verbose = FALSE), keep.rownames = "lineage")
fwrite(lin_dt, fig("04_lineage_ssgsea_scores.tsv"), sep = "\t")
lm_long <- melt(lin_dt, id.vars = "lineage", variable.name = "sample", value.name = "score")
lm_long[, sample := factor(sample, levels = SAMPLE_ORDER, labels = SAMPLES[SAMPLE_ORDER, short])]
lm_long[, z := scale(score)[, 1], by = lineage]
ggsave(fig("05_lineage_scores.png"),
  ggplot(lm_long, aes(sample, lineage, fill = z)) + geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.0f", score)), size = 2.6, colour = "grey15") +
    scale_fill_distiller(palette = "RdBu", limits = c(-2, 2), oob = scales::squish) +
    labs(x = NULL, y = NULL, title = "Haematopoietic lineage content (ssGSEA)",
         subtitle = "Donor lacks the progenitor, granulocyte and erythroid content of all patient libraries") +
    theme_minimal(base_size = 11), width = 7.5, height = 5, dpi = 300)

conf_sets <- lapply(CONFOUND_SIGNATURES, function(x) intersect(x, rownames(logcpm)))
conf_sets <- conf_sets[sapply(conf_sets, length) >= 3]
conf_dt <- as.data.table(GSVA::gsva(GSVA::ssgseaParam(logcpm, conf_sets, normalize = FALSE),
                                    verbose = FALSE), keep.rownames = "signature")
fwrite(conf_dt, fig("04b_confound_signature_scores.tsv"), sep = "\t")
cf <- melt(conf_dt, id.vars = "signature", variable.name = "sample", value.name = "score")
cf[, sample := factor(sample, levels = SAMPLE_ORDER, labels = SAMPLES[SAMPLE_ORDER, short])]
cf[, z := scale(score)[, 1], by = signature]
ggsave(fig("04c_confound_signatures.png"),
  ggplot(cf, aes(sample, signature, fill = z)) + geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.0f", score)), size = 2.7, colour = "grey15") +
    scale_fill_distiller(palette = "RdBu", limits = c(-1.6, 1.6), oob = scales::squish) +
    labs(x = NULL, y = NULL, title = "Confound signatures: the competing explanations",
         subtitle = "D30 is a regenerating marrow 33 days post-transplant. Erythroid content rises across the trajectory.") +
    theme_minimal(base_size = 11), width = 8, height = 3.6, dpi = 300)

################################################################################
# 3. MDS over all samples
################################################################################
message("\n== 3. MDS, all samples ==")
plot_mds_all(logcpm, PFX)

cormat <- cor(logcpm, method = "spearman")
pheatmap(cormat, display_numbers = TRUE, number_format = "%.3f",
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         annotation_col = ann_df(), annotation_colors = ann_colours(),
         main = "Spearman correlation, logCPM",
         filename = fig("12_dge_sample_correlation.png"), width = 6.5, height = 5.5)

################################################################################
# 4. Per-comparison MDS plots
################################################################################
message("\n== 4. per-comparison MDS ==")
plot_mds_per_contrast(logcpm, PFX)

################################################################################
# 5. Heatmaps (top variable genes)
################################################################################
message("\n== 5. heatmaps ==")
plot_all_heatmaps(logcpm, PFX, label = "genes")

################################################################################
# 6. DGE -- one fit over all samples, six contrasts
################################################################################
message("\n== 6. differential expression ==")
f <- fit_contrasts(dge)

de <- rbindlist(lapply(CONTRASTS$id, function(cid) top_table(f, cid, lcpm1, cpm_mat)))
setkey(de, contrast, gene)

# --- BCV sensitivity band -------------------------------------------------
edger_all <- rbindlist(lapply(THRESH$bcv_grid, function(b) {
  fb <- fit_contrasts(dge, bcv = b)
  rbindlist(lapply(CONTRASTS$id, function(cid)
    top_table(fb, cid, lcpm1, cpm_mat)[, bcv := b]))
}))
fwrite(edger_all, fig("09_edgeR_all_bcv.tsv.gz"), sep = "\t")
sens <- dcast(edger_all[, .(n = sum(sig)), by = .(contrast, bcv)],
              contrast ~ bcv, value.var = "n")
setnames(sens, as.character(THRESH$bcv_grid), paste0("BCV_", THRESH$bcv_grid))
sens <- merge(CONTRASTS[, .(contrast = id, label, class)], sens, by = "contrast")
message("\n  features passing FDR<0.05 & |log2FC|>=1, by assumed dispersion:")
print(sens); fwrite(sens, fig("10_bcv_sensitivity.tsv"), sep = "\t")

# --- is the assumed dispersion defensible? --------------------------------
emp_bcv <- function(y, lab) {
  y <- normLibSizes(y, method = "TMM")
  d <- estimateDisp(y, design = matrix(1, ncol(y), 1), robust = TRUE)
  data.table(set = lab, n_libraries = ncol(y),
             common_bcv = round(sqrt(d$common.dispersion), 3),
             trended_bcv = round(median(sqrt(d$trended.dispersion)), 3),
             tagwise_bcv = round(median(sqrt(d$tagwise.dispersion)), 3))
}
bcv_cal <- rbind(emp_bcv(dge, "all four libraries"),
                 emp_bcv(dge[, PATIENT_SAMPLES], "patient timepoints only"))
bcv_cal[, assumed_primary := THRESH$bcv_primary]
print(bcv_cal); fwrite(bcv_cal, fig("10b_dispersion_calibration.tsv"), sep = "\t")

# A second summary of the same question. The pooled estimates above are
# dominated by the genes that differ most between libraries; this asks how much
# a TYPICAL expressed gene varies. Median SD of log2 CPM among genes averaging
# log2 CPM > 3, corrected for the downward bias of a median SD from so few
# libraries, and put on the BCV scale (x ln 2; Poisson noise is negligible at
# these counts). Also an upper bound, since the libraries genuinely differ.
typical_spread <- function(lc, lab, min_ave = 3) {
  ae <- rowMeans(lc); s <- apply(lc, 1, sd); k <- ae > min_ave
  df <- ncol(lc) - 1; corr <- sqrt(qchisq(0.5, df) / df)
  data.table(set = lab, n_libraries = ncol(lc), genes = sum(k),
             median_sd_log2cpm = round(median(s[k]), 3),
             implied_bcv = round(median(s[k]) / corr * log(2), 3),
             pct_genes_below_assumed = round(100 * mean(s[k] / corr * log(2) < THRESH$bcv_primary), 1))
}
spread_tab <- rbind(typical_spread(logcpm, "all four libraries"),
                    typical_spread(logcpm[, PATIENT_SAMPLES], "patient timepoints only"))
print(spread_tab); fwrite(spread_tab, fig("10c_typical_gene_spread.tsv"), sep = "\t")

# --- second approach: NOISeq ----------------------------------------------
message("\n  NOISeq (replicates = 'no')")
noiseq_one <- function(cid) {
  cc <- CONTRASTS[cid]
  sub <- counts[, c(cc$num, cc$den), drop = FALSE]
  fac <- data.frame(group = factor(c(cc$num, cc$den), levels = c(cc$den, cc$num)))
  rownames(fac) <- colnames(sub)
  r <- NOISeq::noiseq(NOISeq::readData(data = sub, factors = fac), factor = "group",
                      k = 0.5, norm = "tmm", replicates = "no", nss = 5,
                      lc = 1, pnr = 0.2, v = 0.02)
  out <- as.data.table(r@results[[1]], keep.rownames = "gene")
  out[, .(gene, noiseq_prob = prob, contrast = cid)]
}
noiseq_all <- rbindlist(lapply(CONTRASTS$id, noiseq_one), fill = TRUE)
fwrite(noiseq_all, fig("11_noiseq_all.tsv.gz"), sep = "\t")

prob_curve <- rbindlist(lapply(c(0.90, 0.95, 0.99, 0.999, 0.9999), function(q)
  noiseq_all[, .(cutoff = q, n = sum(noiseq_prob >= q, na.rm = TRUE)), by = contrast]))
prob_curve <- dcast(prob_curve, contrast ~ cutoff, value.var = "n")
print(prob_curve); fwrite(prob_curve, fig("11b_noiseq_threshold_curve.tsv"), sep = "\t")

de <- merge(de, noiseq_all, by = c("gene", "contrast"), all.x = TRUE)
de[, noiseq_sig := !is.na(noiseq_prob) & noiseq_prob >= THRESH$noiseq_prob]
de[, concordant := sig & noiseq_sig]
de[, robust := concordant & !on_off]

conc <- de[, .(edgeR = sum(sig), NOISeq = sum(noiseq_sig), both = sum(concordant),
               on_off = sum(concordant & on_off), robust = sum(robust)), by = contrast]
conc <- merge(CONTRASTS[, .(contrast = id, label, class)], conc, by = "contrast")
message("\n  agreement between the two approaches:"); print(conc)
fwrite(conc, fig("12_engine_concordance.tsv"), sep = "\t")

# --- the ruler ------------------------------------------------------------
ruler_genes <- de[contrast == RULER_CONTRAST & robust == TRUE, gene]
de[, ruler_flagged := gene %in% ruler_genes & contrast != RULER_CONTRAST]
ruler_tab <- de[robust == TRUE, .(robust = .N, in_ruler = sum(ruler_flagged),
                                  residual = sum(!ruler_flagged)), by = contrast]
ruler_tab[, pct_ruler := round(100 * in_ruler / robust, 1)]
ruler_tab <- merge(CONTRASTS[, .(contrast = id, label, class)], ruler_tab, by = "contrast")
message("\n  ruler subtraction:"); print(ruler_tab)
fwrite(ruler_tab, fig("12b_ruler_subtraction.tsv"), sep = "\t")

# --- one topTable per contrast -------------------------------------------
out_files <- setNames(sprintf("%02d_dge_%s.txt", 17:22, CONTRASTS$id), CONTRASTS$id)
for (cid in CONTRASTS$id) {
  tt <- de[contrast == cid][order(P.Value)]
  fwrite(tt[, .(gene, logFC, lfc_shrunk, AveExpr, t, P.Value, adj.P.Val,
                cpm_num, cpm_den, on_off, noiseq_prob, sig, concordant, robust,
                ruler_flagged)],
         fig(out_files[cid]), sep = "\t", quote = FALSE)
}

################################################################################
# 6b. Classification across the two axes (case-specific)
################################################################################
message("\n== 6b. axis classification ==")
wide <- dcast(de, gene ~ contrast, value.var = c("lfc_shrunk", "robust", "ruler_flagged"))
cls <- wide[, .(gene, lfc_clone = lfc_shrunk_C1_clone,
                lfc_transform = lfc_shrunk_C2_transform,
                lfc_cumul = lfc_shrunk_C3_cumulative,
                sig_clone = robust_C1_clone, sig_transform = robust_C2_transform,
                ruler_clone = ruler_flagged_C1_clone,
                ruler_transform = ruler_flagged_C2_transform)]
for (j in c("sig_clone", "sig_transform", "ruler_clone", "ruler_transform"))
  set(cls, which(is.na(cls[[j]])), j, FALSE)
cls[, flat_clone     := abs(lfc_clone)     < THRESH$flat_lfc]
cls[, flat_transform := abs(lfc_transform) < THRESH$flat_lfc]
cls[, class := fifelse(sig_clone & !sig_transform & flat_transform, "clone-associated",
              fifelse(!sig_clone & flat_clone & sig_transform, "transformation-associated",
              fifelse(sig_clone & sig_transform & sign(lfc_clone) == sign(lfc_transform), "progressive",
              fifelse(sig_clone & sig_transform, "reversing",
              fifelse(sig_clone | sig_transform, "ambiguous", "unchanged")))))]
lin_map <- rbindlist(lapply(names(LINEAGE_MARKERS), function(l)
  data.table(gene = LINEAGE_MARKERS[[l]], lineage = l)))
cls <- merge(cls, lin_map, by = "gene", all.x = TRUE)
cls[, ruler_suspect := (class == "clone-associated" & ruler_clone) |
                       (class == "transformation-associated" & ruler_transform)]
fwrite(cls, fig("13_gene_axis_classification.tsv"), sep = "\t")
cls_n <- cls[, .(N = .N, ruler_suspect = sum(ruler_suspect)), by = class][order(-N)]
cls_n[, clean := N - ruler_suspect]
print(cls_n); fwrite(cls_n, fig("14_axis_class_counts.tsv"), sep = "\t")

# The two axes share the 3.5-year library with opposite signs (C1 = 3.5yr - D30,
# C2 = DDHN - 3.5yr), so anything peculiar to that one library -- noise or real --
# moves a gene up on one axis and down on the other. A negative correlation is
# therefore partly built in and must not be read as biology.
r_all  <- cls[, cor(lfc_clone, lfc_transform, use = "complete.obs")]
r_move <- cls[class != "unchanged", cor(lfc_clone, lfc_transform, use = "complete.obs")]
ggsave(fig("16_axis_scatter.png"),
  ggplot(cls[class != "unchanged"], aes(lfc_clone, lfc_transform, colour = class)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(alpha = .6, size = 1.4) + scale_colour_brewer(palette = "Set1") +
    labs(x = "log2FC  clonal establishment  (3.5 years - 1 month)",
         y = "log2FC  transformation  (DDHN diagnosis - 3.5 years)",
         title = "The two axes are not independent",
         subtitle = sprintf(paste0("r = %.2f across all genes, %.2f among genes that moved.\n",
                                   "Both axes contain the 2-year library with opposite signs, so part of this is built in."),
                            r_all, r_move)) +
    theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 9)),
  width = 8.5, height = 6.8, dpi = 300)

################################################################################
# 7a. MD plots per contrast
################################################################################
message("\n== 7a. MD plots ==")
for (cid in CONTRASTS$id)
  plot_md_contrast(de[contrast == cid], cid, fig(sprintf("21_MD_%s.png", cid)))

################################################################################
# 7b. Per-contrast MA plots (replacing volcanoes -- see dge_functions.R)
################################################################################
message("\n== 7b. MA plots ==")
for (cid in CONTRASTS$id) plot_ma_contrast(de[contrast == cid], cid, PFX)

################################################################################
# 7c. Per-contrast heatmap of the differential features
################################################################################
message("\n== 7c. per-contrast heatmaps ==")
for (cid in CONTRASTS$id) {
  gs <- de[contrast == cid & robust == TRUE][order(adj.P.Val), gene]
  plot_contrast_heatmap(logcpm, gs, cid, fig(sprintf("23_heatmap_%s.png", cid)))
}

################################################################################
# 7d. Focus gene panels
################################################################################
message("\n== 7d. focus genes ==")
focus_dt <- rbindlist(lapply(names(FOCUS_GENES), function(k)
  data.table(panel = k, gene = FOCUS_GENES[[k]])))
focus_dt <- focus_dt[gene %in% rownames(logcpm)]
fwrite(merge(focus_dt, de[gene %in% focus_dt$gene], by = "gene")[order(panel, contrast, P.Value)],
       fig("31_focus_gene_dge.tsv"), sep = "\t")
for (pn in unique(focus_dt$panel)) {
  gs <- intersect(FOCUS_GENES[[pn]], rownames(logcpm))
  if (length(gs) < 2) next
  m <- logcpm[gs, , drop = FALSE]; colnames(m) <- SAMPLES[colnames(m), short]
  pheatmap(m, scale = "row", cluster_cols = FALSE,
           color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
           main = paste0(gsub("_", " ", pn), "  (row-scaled logCPM)"),
           filename = fig(sprintf("32_focus_%s.png", pn)),
           width = 5.5, height = max(3, length(gs) * 0.22 + 1.6))
}

################################################################################
# 8. Save state (enrichment runs from 04_enrichment.R)
################################################################################
save(dge, logcpm, lcpm1, cpm_mat, counts, anno, de, edger_all, noiseq_all,
     cls, cls_n, sens, conc, strip, lin_dt, conf_dt, cormat, norm_tab, bcv_cal,
     spread_tab, ruler_tab, prob_curve, f,
     file = file.path(OUT, "dge_gene_level.RData"))
message("\n== done ==  ", run_stamp(), "\noutputs in ", OUT)
