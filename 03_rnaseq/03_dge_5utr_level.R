################################################################################
# 03_dge_5utr_level.R  --  5'UTR counts, same pipeline as the gene level,
#                          plus differential promoter usage
#
# Run:  Rscript 03_rnaseq/03_dge_5utr_level.R
# Out:  05_processed_counts_5utr/
#
# STEP 0: PROMOTER CLUSTERING (cluster_utr_promoters() in 01_load_counts.R)
#   featureCounts was run with -O, so nested and overlapping 5'UTRs of one gene
#   count the same reads. Features shorter than THRESH$utr_min_width are
#   dropped, overlapping features within a gene are merged into clusters, and
#   each cluster takes its longest member's counts (members are never summed).
#   Everything below runs on the clusters.
#
# PART A (sections 0-7) mirrors 02_dge_gene_level.R exactly, on promoter
#   clusters summed to one value per GENE -- the gene's total 5'UTR signal.
#   Every figure and table produced for the gene level is produced here too,
#   under the same numbering. Each point in the Part A plots is a gene.
#
#   Summing the gene's disjoint promoter clusters gives its total 5' signal.
#
# PART B (section 10): differential PROMOTER USAGE between a gene's clusters.
#   Gene-level counting answers "how much does this gene fire?"; it cannot
#   answer "which start site does it fire from?", which is the question DNMT3A
#   and TET2 bear on. A call counts only if the switched promoter is detected in
#   both libraries of the pair -- the gene-level detection rule (section 10 note).
#
# DEPTH. A 5'UTR collects a few per cent of a gene's reads, so these libraries
#   are small (about 2.8-10.8 M reads after clustering) and 1 CPM is only 3-11
#   reads. Many 5'UTR-level changes rest on a handful of reads in one library;
#   read support is reported alongside the promoter counts.
#
# 5' COVERAGE. Degradation removes 5' coverage first, so 5'UTR signal is checked
#   per library against gene-level signal over the same genes, both from raw
#   counts.
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma)
  library(org.Hs.eg.db); library(RColorBrewer); library(pheatmap)
  library(ggplot2); library(NOISeq); library(GSVA)
})

SCRIPTS <- "03_rnaseq"
source(file.path(SCRIPTS, "00_config.R"))
source(file.path(SCRIPTS, "01_load_counts.R"))
source(file.path(SCRIPTS, "dge_functions.R"))

OUT <- PATHS$out_5utr
fig <- function(...) file.path(OUT, paste0(...))
PFX <- paste0(OUT, "/")
set.seed(1234)

################################################################################
# 0. Load 5'UTR counts and cluster them into promoters
################################################################################
message("\n== 0. loading 5'UTR counts ==")
u  <- load_utr_counts()
pc <- cluster_utr_promoters(u)
pc_qc <- pc$qc
fwrite(pc_qc, fig("00b_promoter_clustering.tsv"), sep = "\t")
fwrite(pc$members, fig("00c_promoter_cluster_membership.tsv.gz"), sep = "\t")

usym <- pc$anno$gene_name
usym[is.na(usym) | usym == ""] <- pc$anno$Geneid[is.na(usym) | usym == ""]
strip <- strip_confounders(pc$counts, usym)
feat_counts <- strip$counts
feat_anno   <- pc$anno[strip$keep]
fwrite(strip$fractions, fig("00_library_composition.tsv"), sep = "\t")

# ---- 5' coverage diagnostic ----------------------------------------------
# Like-for-like: raw clustered 5'UTR counts against raw gene counts, over the
# SAME genes (those carrying an annotated 5'UTR), both before any filtering.
utr_by_gene <- rowsum(pc$counts, group = pc$anno$Geneid, reorder = FALSE)
graw <- load_gene_counts(verbose = FALSE)$counts
common <- intersect(rownames(utr_by_gene), rownames(graw))
cov_bias <- data.table(
  sample      = SAMPLE_ORDER,
  utr_counts  = colSums(utr_by_gene[common, SAMPLE_ORDER, drop = FALSE]),
  gene_counts = colSums(graw[common, SAMPLE_ORDER, drop = FALSE]))
cov_bias[, utr_fraction_pct := round(100 * utr_counts / gene_counts, 2)]
cov_bias[, rel_to_median := round(utr_fraction_pct / median(utr_fraction_pct), 2)]
# flag only a substantial deficit
cov_bias[, five_prime_ok := rel_to_median > 0.6]
print(cov_bias); fwrite(cov_bias, fig("49_five_prime_coverage_bias.tsv"), sep = "\t")
DEGRADED <- cov_bias[five_prime_ok == FALSE, sample]
if (length(DEGRADED))
  message("  5'-depleted: ", paste(DEGRADED, collapse = ", "),
          " -- absolute 5'UTR activity unreliable in contrasts using them")

# ---- collapse to gene level for Part A ------------------------------------
gene_counts <- rowsum(feat_counts, group = feat_anno$Geneid, reorder = FALSE)
g2s <- feat_anno[!duplicated(Geneid), .(Geneid, gene_name)]
setkey(g2s, Geneid)
gsym0 <- g2s[rownames(gene_counts), gene_name]
gsym0[is.na(gsym0) | gsym0 == ""] <- rownames(gene_counts)[is.na(gsym0) | gsym0 == ""]
message(sprintf("  collapsed %s promoter clusters -> %s genes",
                format(nrow(feat_counts), big.mark = ","),
                format(nrow(gene_counts), big.mark = ",")))

################################################################################
# 1. Filtering
################################################################################
message("\n== 1. filtering ==")
keep <- filter_expressed(gene_counts)
counts <- gene_counts[keep, , drop = FALSE]
gsym   <- gsym0[keep]
ens    <- rownames(gene_counts)[keep]
rownames(counts) <- make.unique(gsym)
check_sample_order(counts)

################################################################################
# 2. DGEList + TMM
################################################################################
message("\n== 2. DGEList + TMM ==")
dge <- DGEList(counts = counts,
               genes = data.frame(symbol = rownames(counts), ensembl = ens))
dge <- normLibSizes(dge, method = "TMM")

log_raw_dgeLlist <- cpm(DGEList(counts), log = TRUE, prior.count = 3)
logcpm  <- cpm(dge, log = TRUE, prior.count = 3)
cpm_mat <- cpm(dge)
lcpm1   <- log2(cpm_mat + 1)

plot_library_sizes(dge, fig("11_dge_library_sizes.png"))
plot_expression_boxplots(log_raw_dgeLlist, logcpm, fig("11_dge_expression_boxplots.png"))
plot_meanvar(dge, fig("11_dge_mean_variance.png"))

norm_tab <- data.table(sample = colnames(dge), lib_size = dge$samples$lib.size,
                       norm_factor = round(dge$samples$norm.factors, 4))
print(norm_tab); fwrite(norm_tab, fig("01_normalisation_factors.tsv"), sep = "\t")
write.table(round(logcpm, 4), fig("rnaSeq_normalized_counts_collated.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

################################################################################
# 3-5. MDS, correlation, heatmaps
################################################################################
message("\n== 3. MDS, all samples ==")
plot_mds_all(logcpm, PFX)
cormat <- cor(logcpm, method = "spearman")
pheatmap(cormat, display_numbers = TRUE, number_format = "%.3f",
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         annotation_col = ann_df(), annotation_colors = ann_colours(),
         main = "Spearman correlation, 5'UTR logCPM",
         filename = fig("12_dge_sample_correlation.png"), width = 6.5, height = 5.5)

message("\n== 4. per-comparison MDS ==")
plot_mds_per_contrast(logcpm, PFX)

message("\n== 5. heatmaps ==")
plot_all_heatmaps(logcpm, PFX, label = "genes (5'UTR signal)")

################################################################################
# 6. DGE
################################################################################
message("\n== 6. differential 5'UTR activity ==")
f <- fit_contrasts(dge)
de <- rbindlist(lapply(CONTRASTS$id, function(cid) top_table(f, cid, lcpm1, cpm_mat)))
setkey(de, contrast, gene)

edger_all <- rbindlist(lapply(THRESH$bcv_grid, function(b) {
  fb <- fit_contrasts(dge, bcv = b)
  rbindlist(lapply(CONTRASTS$id, function(cid) top_table(fb, cid, lcpm1, cpm_mat)[, bcv := b]))
}))
fwrite(edger_all, fig("09_edgeR_all_bcv.tsv.gz"), sep = "\t")
sens <- dcast(edger_all[, .(n = sum(sig)), by = .(contrast, bcv)],
              contrast ~ bcv, value.var = "n")
setnames(sens, as.character(THRESH$bcv_grid), paste0("BCV_", THRESH$bcv_grid))
sens <- merge(CONTRASTS[, .(contrast = id, label, class)], sens, by = "contrast")
print(sens); fwrite(sens, fig("10_bcv_sensitivity.tsv"), sep = "\t")

message("\n  NOISeq")
noiseq_one <- function(cid) {
  cc <- CONTRASTS[cid]
  sub <- counts[, c(cc$num, cc$den), drop = FALSE]
  fac <- data.frame(group = factor(c(cc$num, cc$den), levels = c(cc$den, cc$num)))
  rownames(fac) <- colnames(sub)
  r <- NOISeq::noiseq(NOISeq::readData(data = sub, factors = fac), factor = "group",
                      k = 0.5, norm = "tmm", replicates = "no", nss = 5,
                      lc = 1, pnr = 0.2, v = 0.02)
  as.data.table(r@results[[1]], keep.rownames = "gene")[, .(gene, noiseq_prob = prob,
                                                            contrast = cid)]
}
noiseq_all <- rbindlist(lapply(CONTRASTS$id, noiseq_one), fill = TRUE)
fwrite(noiseq_all, fig("11_noiseq_all.tsv.gz"), sep = "\t")

de <- merge(de, noiseq_all, by = c("gene", "contrast"), all.x = TRUE)
de[, noiseq_sig := !is.na(noiseq_prob) & noiseq_prob >= THRESH$noiseq_prob]
de[, concordant := sig & noiseq_sig]
de[, robust := concordant & !on_off]
# flag every contrast that leans on a 5'-depleted library
de[, five_prime_confounded := CONTRASTS[contrast, num] %in% DEGRADED |
                              CONTRASTS[contrast, den] %in% DEGRADED]

conc <- de[, .(edgeR = sum(sig), NOISeq = sum(noiseq_sig), both = sum(concordant),
               robust = sum(robust), five_prime_confounded = five_prime_confounded[1]),
           by = contrast]
conc <- merge(CONTRASTS[, .(contrast = id, label, class)], conc, by = "contrast")
message("\n  agreement between the two approaches:"); print(conc)
fwrite(conc, fig("12_engine_concordance.tsv"), sep = "\t")

ruler_genes <- de[contrast == RULER_CONTRAST & robust == TRUE, gene]
de[, ruler_flagged := gene %in% ruler_genes & contrast != RULER_CONTRAST]
ruler_tab <- de[robust == TRUE, .(robust = .N, in_ruler = sum(ruler_flagged),
                                  residual = sum(!ruler_flagged)), by = contrast]
ruler_tab[, pct_ruler := round(100 * in_ruler / robust, 1)]
ruler_tab <- merge(CONTRASTS[, .(contrast = id, label, class)], ruler_tab, by = "contrast")
print(ruler_tab); fwrite(ruler_tab, fig("12b_ruler_subtraction.tsv"), sep = "\t")

out_files <- setNames(sprintf("%02d_dge_%s.txt", 17:22, CONTRASTS$id), CONTRASTS$id)
for (cid in CONTRASTS$id)
  fwrite(de[contrast == cid][order(P.Value),
            .(gene, logFC, lfc_shrunk, AveExpr, t, P.Value, adj.P.Val, cpm_num,
              cpm_den, on_off, noiseq_prob, sig, concordant, robust, ruler_flagged,
              five_prime_confounded)],
         fig(out_files[cid]), sep = "\t", quote = FALSE)

################################################################################
# 7. MD plots, MA plots, per-contrast heatmaps
################################################################################
message("\n== 7a. MD plots ==")
for (cid in CONTRASTS$id)
  plot_md_contrast(de[contrast == cid], cid, fig(sprintf("21_MD_%s.png", cid)))

message("\n== 7b. MA plots ==")
for (cid in CONTRASTS$id) plot_ma_contrast(de[contrast == cid], cid, PFX)

message("\n== 7c. per-contrast heatmaps ==")
for (cid in CONTRASTS$id)
  plot_contrast_heatmap(logcpm, de[contrast == cid & robust == TRUE][order(adj.P.Val), gene],
                        cid, fig(sprintf("23_heatmap_%s.png", cid)))

save(dge, logcpm, lcpm1, cpm_mat, counts, de, edger_all, noiseq_all, sens, conc,
     strip, cormat, norm_tab, ruler_tab, cov_bias, DEGRADED, f, pc_qc,
     file = file.path(OUT, "dge_5utr_level.RData"))

################################################################################
# PART B
# 10. Differential PROMOTER USAGE (between a gene's promoter clusters)
#
#   Relative usage is a within-gene ratio: invariant to library size and to
#   normalisation. It is NOT immune to cell composition -- different cell types
#   can use different promoters of the same gene, so a change in the cell mix
#   shifts usage without any cell changing -- which is why the reference
#   contrast is applied here too.
#
#   Tested with diffSpliceDGE on a fit at the SAME supplied dispersion as the
#   gene level (BCV 0.4): does a promoter's log fold-change differ from its
#   gene's? Its p-values are therefore conditional on the assumed dispersion
#   exactly as the gene-level ones are; they rank, they do not control error.
#   Calls are made per promoter (exon-level FDR); a gene counts as switching if
#   any of its promoters does. Gene-level Simes p-values are not used.
#
#   DETECTION RULE -- the same one the gene-level analysis uses. A call only
#   counts if the switched promoter is detected (>= THRESH$cpm CPM) in BOTH
#   libraries of the pair. Without it, 74-83% of switch calls were a promoter
#   with no reads at all in one library: an unreplicated on/off observation,
#   exactly what the gene-level 'robust' set excludes. Every call before this
#   rule is kept in 53_ (switched_any), and the read support of the calls is
#   tabulated in 55b_.
#
#   PURE SWITCH is judged on WHOLE-GENE expression from 02_dge_gene_level.R: the
#   gene shifts promoter usage while its whole-gene expression is not called
#   and changes by less than THRESH$flat_lfc.
################################################################################
message("\n############ PART B: promoter usage ############")

fk <- rowSums(feat_counts) >= THRESH$utr_min_count & filter_expressed(feat_counts, verbose = FALSE)
fcnt <- feat_counts[fk, , drop = FALSE]
fann <- feat_anno[fk]
fann[, n_kept := .N, by = Geneid]
message(sprintf("  %s promoter clusters across %s genes; %s genes retain >=2 promoters",
                format(nrow(fcnt), big.mark = ","),
                format(uniqueN(fann$Geneid), big.mark = ","),
                format(uniqueN(fann[n_kept >= 2, Geneid]), big.mark = ",")))

dgu <- DGEList(counts = fcnt, genes = as.data.frame(fann))
dgu <- normLibSizes(dgu, method = "TMM")
fitu <- glmFit(dgu, fit_contrasts(dgu)$design, dispersion = THRESH$bcv_primary^2)
cmu  <- fit_contrasts(dgu)$cm
cpmu <- cpm(dgu)

usage <- as.data.table(cpmu, keep.rownames = "feature_id")
usage <- merge(usage, fann[, .(feature_id, Geneid, gene_name)], by = "feature_id")
for (s in SAMPLE_ORDER) usage[, (paste0("u_", s)) := get(s) / sum(get(s)), by = Geneid]
fwrite(usage, fig("51_promoter_relative_usage.tsv.gz"), sep = "\t")

switchable <- fann[n_kept >= 2, which = TRUE]
fit_sw <- fitu[switchable, ]

ds_feat <- list()
for (cid in CONTRASTS$id) {
  ds <- diffSpliceDGE(fit_sw, contrast = cmu[, cid], geneid = "Geneid",
                      exonid = "feature_id", verbose = FALSE)
  ft <- as.data.table(topSpliceDGE(ds, test = "exon", n = Inf))
  cc <- CONTRASTS[cid]
  ft[, d_usage := usage[match(ft$feature_id, feature_id), get(paste0("u_", cc$num))] -
                  usage[match(ft$feature_id, feature_id), get(paste0("u_", cc$den))]]
  ft[, contrast := cid]
  ds_feat[[cid]] <- ft
}
ds_feat <- rbindlist(ds_feat, fill = TRUE)
ds_feat[, `:=`(num = CONTRASTS[contrast, num], den = CONTRASTS[contrast, den])]
ds_feat[, `:=`(cpm_num   = cpmu[cbind(feature_id, num)], cpm_den   = cpmu[cbind(feature_id, den)],
               reads_num = fcnt[cbind(feature_id, num)], reads_den = fcnt[cbind(feature_id, den)])]
ds_feat[, detected_both := pmin(cpm_num, cpm_den) >= THRESH$cpm]
# every call before the detection rule, kept for transparency
ds_feat[, switched_any := FDR < THRESH$fdr & abs(d_usage) >= THRESH$utr_dprop]
ds_feat[, switched := switched_any & detected_both]
fwrite(ds_feat, fig("53_promoter_switch_feature_level.tsv.gz"), sep = "\t")

read_support <- ds_feat[switched_any == TRUE, .(
  genes_all_calls               = uniqueN(Geneid),
  genes_detected_in_both        = uniqueN(Geneid[detected_both == TRUE]),
  genes_10_reads_in_both        = uniqueN(Geneid[pmin(reads_num, reads_den) >= 10]),
  calls_zero_in_one_library_pct = round(100 * mean(pmin(reads_num, reads_den) == 0), 1)),
  by = contrast]
read_support <- merge(CONTRASTS[, .(contrast = id, label, class)], read_support, by = "contrast")
message("\n  promoter switch calls by read support:"); print(read_support)
fwrite(read_support, fig("55b_promoter_switch_read_support.tsv"), sep = "\t")

switch_genes <- ds_feat[switched == TRUE,
                        .(n_features = .N,
                          max_d_usage = d_usage[which.max(abs(d_usage))],
                          top_feature = feature_id[which.max(abs(d_usage))]),
                        by = .(contrast, Geneid, gene_name)]
ruler_sw <- switch_genes[contrast == RULER_CONTRAST, unique(Geneid)]
switch_genes[, ruler_flagged := Geneid %in% ruler_sw & contrast != RULER_CONTRAST]
fwrite(switch_genes, fig("54_promoter_switching_genes.tsv"), sep = "\t")

sw_n <- switch_genes[, .(switching_genes = .N, in_ruler = sum(ruler_flagged),
                         residual = sum(!ruler_flagged)), by = contrast]
sw_n <- merge(CONTRASTS[, .(contrast = id, label, class)], sw_n, by = "contrast", all.x = TRUE)
for (j in c("switching_genes", "in_ruler", "residual"))
  set(sw_n, which(is.na(sw_n[[j]])), j, 0L)
message("\n  genes with a promoter usage shift:"); print(sw_n)
fwrite(sw_n, fig("55_promoter_switch_counts.tsv"), sep = "\t")

# switching without a change in whole-gene expression: the gene makes about the
# same amount of RNA from a different start site. Judged on the WHOLE-GENE
# analysis (02_dge_gene_level.R), not on 5'UTR signal; genes not tested at gene
# level are not counted as pure switches.
gl <- new.env(); load(file.path(PATHS$out_gene, "dge_gene_level.RData"), envir = gl)
gde <- as.data.table(gl$de)[, .(gene_name = gene, contrast, gene_lfc = lfc_shrunk, gene_sig = sig)]
pure <- merge(switch_genes, gde, by = c("gene_name", "contrast"), all.x = TRUE)
pure[, tested_at_gene_level := !is.na(gene_sig)]
pure[, pure_switch := tested_at_gene_level & !gene_sig & abs(gene_lfc) < THRESH$flat_lfc]
fwrite(pure, fig("56_promoter_switch_vs_gene_expression.tsv"), sep = "\t")
print(pure[, .(switching = .N, pure_switch = sum(pure_switch),
               not_tested_at_gene_level = sum(!tested_at_gene_level)), by = contrast])

ggsave(fig("59_promoter_switch_summary.png"),
  ggplot(sw_n, aes(factor(label, levels = CONTRASTS$label), switching_genes, fill = class)) +
    geom_col(width = .7) + geom_text(aes(label = switching_genes), vjust = -.3, size = 3.2) +
    scale_fill_manual(values = PAL$class) +
    scale_x_discrete(labels = function(x) sub(" \\(", "\n(", x)) +
    labs(x = NULL, y = "genes with a promoter usage shift",
         title = "Differential promoter usage across contrasts",
         subtitle = "Switched promoter detected (>= 1 CPM) in both libraries of the pair; the reference contrast (blue) returns switching too") +
    theme_bw(base_size = 10),
  width = 9, height = 5.5, dpi = 300)

# per-gene promoter usage plots for the top switchers
plot_gene_promoters <- function(gn, cid, file) {
  fs <- fann[gene_name == gn, feature_id]
  if (length(fs) < 2) return(invisible())
  d <- melt(usage[feature_id %in% fs, c("feature_id", paste0("u_", SAMPLE_ORDER)), with = FALSE],
            id.vars = "feature_id", variable.name = "sample", value.name = "usage")
  d[, sample := factor(sub("^u_", "", sample), levels = SAMPLE_ORDER,
                       labels = SAMPLES[SAMPLE_ORDER, short])]
  d[, promoter := sub("^[^:]+\\|[^:]+:", "", feature_id)]
  ggsave(file,
    ggplot(d, aes(sample, usage, fill = promoter)) + geom_col(width = .75) +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_viridis_d(option = "D", end = .92) +
      labs(y = "share of gene 5'UTR signal", x = NULL,
           title = paste0(gn, " - promoter usage"),
           subtitle = paste0("flagged in ", CONTRASTS[cid, label])) +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom", legend.text = element_text(size = 6)) +
      guides(fill = guide_legend(ncol = 2)),
    width = 7, height = 5.2, dpi = 300)
}
dir.create(fig("promoter_plots"), showWarnings = FALSE)
for (cid in CONTRASTS$id) {
  top <- switch_genes[contrast == cid & ruler_flagged == FALSE][order(-abs(max_d_usage))][1:min(8, .N)]
  for (i in seq_len(nrow(top)))
    plot_gene_promoters(top$gene_name[i], cid,
      fig(sprintf("promoter_plots/%s_%02d_%s.png", cid, i,
                  gsub("[^A-Za-z0-9]", "_", top$gene_name[i]))))
}

# BED export for DMR intersection
bed <- merge(ds_feat[switched == TRUE, .(feature_id, contrast, d_usage, FDR)],
             fann[, .(feature_id, chr, start, end, strand, gene_name)], by = "feature_id")
bed[, chr := sub("^chr", "", chr)]
setorder(bed, contrast, chr, start)
for (cid in unique(bed$contrast))
  fwrite(bed[contrast == cid, .(chr, start = start - 1L, end, name = feature_id,
                                score = round(pmin(1000, 1000 * abs(d_usage))), strand)],
         fig(sprintf("57_switching_promoters_%s.bed", cid)), sep = "\t", col.names = FALSE)
tss <- unique(bed[, .(chr, start, end, strand, gene_name, feature_id)])
tss[, tss_pos := fifelse(strand == "+", start, end)]
tss[, `:=`(w_start = pmax(0L, tss_pos - 1000L), w_end = tss_pos + 1000L)]
tss <- unique(tss, by = c("chr", "tss_pos", "strand")); setorder(tss, chr, w_start)
fwrite(tss[, .(chr, w_start, w_end, name = feature_id, score = 0, strand)],
       fig("58_switching_promoters_TSS_1kb.bed"), sep = "\t", col.names = FALSE)
message("  BED files written for DMR intersection")

save(dgu, fann, usage, ds_feat, switch_genes, sw_n, pure, read_support, cov_bias, DEGRADED, pc_qc,
     file = file.path(OUT, "promoter_usage.RData"))
message("\n== done ==  ", run_stamp(), "\noutputs in ", OUT)
