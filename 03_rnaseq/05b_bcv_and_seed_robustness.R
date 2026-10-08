################################################################################
# 05b_bcv_and_seed_robustness.R  --  do the headline results survive the two
#                                    things this analysis cannot pin down?
#
#   Rscript 03_rnaseq/05b_bcv_and_seed_robustness.R        (gene level; run after 05)
#
# 1. THE ASSUMED DISPERSION. Everything was fitted at BCV 0.2, 0.4 and 0.6
#    (02_dge_gene_level.R). Here: how many 'robust' genes survive at each, which
#    of the cleanest contrast's robust genes do, and whether the Hallmark results
#    hold when GSEA is re-ranked at 0.2 and 0.6.
# 2. GSEA SAMPLING. fgsea's p-values and NES normalisation come from random
#    sampling; the reported run uses one fixed seed. Here: 25 seeds at BCV 0.4
#    (that seed and 24 others) for the Hallmark collection, recording the range
#    of NES and adjusted p, and
#    how often p53 passes the full pathway-ruler criterion of 05_pathway_ruler.R
#    (FDR < 0.05 in both large-clone donor contrasts, same direction, ruler raw
#    p >= THRESH$ruler_p).
#
# Out (gene level):
#   84_hallmark_across_bcv.tsv          all sets x 6 contrasts x 3 BCVs, full tables
#   85_robust_genes_across_bcv.tsv      robust calls at each BCV per contrast
#   85b_cleanest_contrast_robust_genes_across_bcv.tsv
#   86_hallmark_seed_stability.tsv      NES / p / FDR range over 25 seeds
#   87_p53_verdict_stability.tsv        p53 criterion at each BCV and each seed
#   88_gsea_determinism_check.tsv       same seed twice -> identical results?
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
  library(msigdbr); library(BiocParallel)
})
SCRIPTS <- "03_rnaseq"
source(file.path(SCRIPTS, "00_config.R"))
source(file.path(SCRIPTS, "dge_functions.R"))
register(SerialParam(), default = TRUE)

OUT <- PATHS$out_gene
fig <- function(...) file.path(OUT, paste0(...))
e <- new.env(); load(file.path(OUT, "dge_gene_level.RData"), envir = e)
de <- as.data.table(e$de); ea <- as.data.table(e$edger_all)

.msig_H <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                    error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
entrez_col <- intersect(c("ncbi_gene", "entrez_gene"), names(.msig_H))[1]
HALLMARK_T2G <- data.frame(gs_name = .msig_H$gs_name, entrez_gene = .msig_H[[entrez_col]],
                           stringsAsFactors = FALSE)
SEED  <- 1234
# the seed behind the reported results plus 24 others, so the range shown
# always contains the values quoted in the report
SEEDS <- c(SEED, 1:24)
P53   <- "HALLMARK_P53_PATHWAY"
DONOR <- CONTRASTS[class == "context", id]

################################################################################
# 1. Genes: robust calls at each assumed dispersion
#    NOISeq and the on/off flag do not depend on the BCV; only the edgeR call does.
################################################################################
message("\n== 1. robust genes across BCV ==")
nq <- de[, .(gene, contrast, noiseq_sig, robust_04 = robust)]
rb <- rbindlist(lapply(THRESH$bcv_grid, function(b) {
  x <- merge(ea[bcv == b, .(gene, contrast, sig, on_off)], nq, by = c("gene", "contrast"))
  x[, robust_b := sig & noiseq_sig & !on_off]
  x[, .(bcv = b, robust = sum(robust_b), robust_at_primary = sum(robust_04),
        primary_set_still_robust = sum(robust_04 & robust_b)), by = contrast]
}))
rb <- merge(CONTRASTS[, .(contrast = id, label)], rb, by = "contrast")
print(rb); fwrite(rb, fig("85_robust_genes_across_bcv.tsv"), sep = "\t")

cc <- merge(ea[contrast == CLEANEST_CONTRAST, .(gene, bcv, sig, adj.P.Val)],
            de[contrast == CLEANEST_CONTRAST & robust == TRUE, .(gene, lfc_shrunk, cpm_num, cpm_den)],
            by = "gene")
cc <- dcast(cc, gene + lfc_shrunk + cpm_num + cpm_den ~ bcv, value.var = c("sig", "adj.P.Val"))
fwrite(cc[order(lfc_shrunk)], fig("85b_cleanest_contrast_robust_genes_across_bcv.tsv"), sep = "\t")

################################################################################
# 2. Hallmark GSEA re-ranked at each BCV (fixed seed), full tables
################################################################################
message("\n== 2. Hallmark across BCV ==")
# primary BCV ranks come from `de`, exactly as 04_enrichment.R builds them
rank_for <- function(b, cid) {
  tt <- if (isTRUE(all.equal(b, THRESH$bcv_primary))) de[contrast == cid][order(P.Value)]
        else ea[bcv == b & contrast == cid][order(P.Value)]
  rank_genes(prep_contrast(tt))
}
run_h <- function(ranks, seed) {
  set.seed(seed)
  as.data.table(as.data.frame(GSEA(ranks, TERM2GENE = HALLMARK_T2G, pAdjustMethod = "BH",
                                   pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE)))[
    , .(ID, NES, pvalue, p.adjust)]
}
RANKS <- list()
for (b in THRESH$bcv_grid) for (cid in CONTRASTS$id) RANKS[[paste(b, cid)]] <- rank_for(b, cid)
key <- function(b, cid) paste(b, cid)

across <- rbindlist(lapply(THRESH$bcv_grid, function(b) rbindlist(lapply(CONTRASTS$id, function(cid)
  run_h(RANKS[[key(b, cid)]], SEED)[, `:=`(bcv = b, contrast = cid)]))))
fwrite(across, fig("84_hallmark_across_bcv.tsv"), sep = "\t")

stored <- rbindlist(lapply(CONTRASTS$id, function(cid) {
  f <- file.path(OUT, "enrichment_gsea", cid, "gsea_full_HALLMARK.csv")
  if (!file.exists(f)) return(NULL)
  as.data.table(read.csv(f))[, .(ID, NES_04 = NES, contrast = cid)]
}))
if (nrow(stored)) {
  chk <- merge(across[bcv == THRESH$bcv_primary], stored, by = c("ID", "contrast"))
  dmax <- max(abs(chk$NES - chk$NES_04))
  message("  primary-BCV Hallmark vs 04_enrichment.R output: ", nrow(chk),
          " rows, max |dNES| = ", signif(dmax, 3))
  if (dmax > 1e-8) warning("Hallmark results here differ from 04_enrichment.R")
}

################################################################################
# 3. Determinism: the same seed twice must give identical results
################################################################################
message("\n== 3. determinism check ==")
r2 <- RANKS[[key(THRESH$bcv_primary, CLEANEST_CONTRAST)]]
h1 <- run_h(r2, SEED); h2 <- run_h(r2, SEED)
bp <- function() {
  set.seed(SEED)
  as.data.table(as.data.frame(gseGO(r2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                                    pAdjustMethod = "BH", pvalueCutoff = 1, eps = 0,
                                    seed = TRUE, verbose = FALSE)))[order(ID), .(ID, NES)]
}
b1 <- bp(); b2 <- bp()
det <- data.table(check = c("Hallmark, same seed twice", "GO-BP, same seed twice"),
                  contrast = CLEANEST_CONTRAST,
                  identical = c(identical(h1$NES, h2$NES), identical(b1$NES, b2$NES)),
                  max_abs_dNES = c(max(abs(h1$NES - h2$NES)), max(abs(b1$NES - b2$NES))))
print(det); fwrite(det, fig("88_gsea_determinism_check.tsv"), sep = "\t")

################################################################################
# 4. GSEA sampling: 25 seeds at the primary BCV
################################################################################
message("\n== 4. Hallmark over ", length(SEEDS), " seeds ==")
HEAD <- unique(c(CLEANEST_CONTRAST, DONOR, RULER_CONTRAST))
seeds <- rbindlist(lapply(SEEDS, function(s) rbindlist(lapply(HEAD, function(cid)
  run_h(RANKS[[key(THRESH$bcv_primary, cid)]], s)[, `:=`(seed = s, contrast = cid)]))))
stab <- seeds[, .(NES_min = round(min(NES), 3), NES_median = round(median(NES), 3),
                  NES_max = round(max(NES), 3), p_min = signif(min(pvalue), 3),
                  p_max = signif(max(pvalue), 3), fdr_min = signif(min(p.adjust), 3),
                  fdr_max = signif(max(p.adjust), 3),
                  pct_seeds_fdr05 = round(100 * mean(p.adjust < THRESH$gsea_fdr), 1)),
              by = .(contrast, ID)]
fwrite(stab[order(contrast, -abs(NES_median))], fig("86_hallmark_seed_stability.tsv"), sep = "\t")

################################################################################
# 5. The p53 criterion at each BCV and each seed
################################################################################
p53_row <- function(tab) {
  x <- tab[ID == P53]
  g <- function(cid, col) x[contrast == cid][[col]]
  data.table(DDHN_NES = g(DONOR[1], "NES"), DDHN_FDR = g(DONOR[1], "p.adjust"),
             yr2_NES = g(DONOR[2], "NES"), yr2_FDR = g(DONOR[2], "p.adjust"),
             ruler_NES = g(RULER_CONTRAST, "NES"), ruler_p = g(RULER_CONTRAST, "pvalue"),
             passes = g(DONOR[1], "p.adjust") < THRESH$gsea_fdr &
                      g(DONOR[2], "p.adjust") < THRESH$gsea_fdr &
                      sign(g(DONOR[1], "NES")) == sign(g(DONOR[2], "NES")) &
                      g(RULER_CONTRAST, "pvalue") >= THRESH$ruler_p)
}
v_bcv  <- rbindlist(lapply(THRESH$bcv_grid, function(b)
  p53_row(across[bcv == b])[, `:=`(run = "BCV", bcv = b, seed = SEED)]))
v_seed <- rbindlist(lapply(SEEDS, function(s)
  p53_row(seeds[seed == s])[, `:=`(run = "seed", bcv = THRESH$bcv_primary, seed = s)]))
v <- rbind(v_bcv, v_seed)
setcolorder(v, c("run", "bcv", "seed"))
fwrite(v, fig("87_p53_verdict_stability.tsv"), sep = "\t")
print(v_bcv)
message(sprintf("  p53 passes the pathway-ruler criterion in %d of %d seeds at BCV %.1f",
                sum(v_seed$passes), nrow(v_seed), THRESH$bcv_primary))
message(sprintf("  over seeds: DDHN NES %.2f to %.2f, FDR max %.3g; 3.5yr NES %.2f to %.2f, FDR max %.3g; ruler p min %.2g",
                min(v_seed$DDHN_NES), max(v_seed$DDHN_NES), max(v_seed$DDHN_FDR),
                min(v_seed$yr2_NES), max(v_seed$yr2_NES), max(v_seed$yr2_FDR), min(v_seed$ruler_p)))

message("\n== done ==  ", run_stamp())
