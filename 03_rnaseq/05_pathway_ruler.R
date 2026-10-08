################################################################################
# 05_pathway_ruler.R  --  apply the ruler at PATHWAY level
#
#   Rscript 03_rnaseq/05_pathway_ruler.R genes     (default)
#   Rscript 03_rnaseq/05_pathway_ruler.R 5utr
#
# WHY THIS EXISTS
#   The ruler contrast (D30 vs donor: no disease, a small clone at both ends)
#   is applied to pathways as it is to genes: a pathway that moves in
#   a donor contrast AND in the ruler is compartment; one that moves in the
#   donor contrasts but NOT in the ruler is not explained by compartment, and
#   is a candidate for something real.
#
# HOW THE RULER IS READ
#   The ruler is tested at RAW p < THRESH$ruler_p (0.10), deliberately lenient:
#   its job is to rule a compartment explanation OUT, so a generous bar makes
#   the surviving set smaller and safer.
#   That needs raw p for EVERY set in the ruler, including those nowhere near
#   significance. gsea_<CAT>.csv cannot supply it: clusterProfiler drops every
#   set above its pvalueCutoff (0.05), so the ruler is read from
#   gsea_full_<CAT>.csv, written by 04_enrichment.R.
#
# WHAT IT PRODUCES
#   80_pathway_ruler.tsv        every significant pathway x contrast, with the
#                               ruler's NES and raw p and the verdict
#   80b_pathway_ruler_summary.tsv   verdict counts per contrast
#   81_clone_associated_pathways.tsv   significant in both donor contrasts, same
#                                      direction, ruler-clean in both
#   82_pathway_ruler_heatmap.png       leukaemia-relevant Hallmark sets, all contrasts
#   83_leukaemia_themes.tsv     that panel as a table, significant or not
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(RColorBrewer); library(stringr)
})

SCRIPTS <- "03_rnaseq"
source(file.path(SCRIPTS, "00_config.R"))

args  <- commandArgs(trailingOnly = TRUE)
LEVEL <- if (length(args)) args[1] else "genes"
stopifnot(LEVEL %in% c("genes", "5utr"))
OUT <- if (LEVEL == "genes") PATHS$out_gene else PATHS$out_5utr
fig <- function(...) file.path(OUT, paste0(...))

CATS <- c("HALLMARK", "BP", "CC", "MF", "KEGG")

################################################################################
# 1. Read every GSEA table (sets at raw p <= 0.05)
################################################################################
message("\n== reading GSEA output (", LEVEL, ") ==")
gs <- rbindlist(lapply(CONTRASTS$id, function(cid)
  rbindlist(lapply(CATS, function(ct) {
    f <- file.path(OUT, "enrichment_gsea", cid, paste0("gsea_", ct, ".csv"))
    if (!file.exists(f)) return(NULL)
    d <- as.data.table(read.csv(f))
    if (!nrow(d)) return(NULL)
    d[, `:=`(contrast = cid, category = ct)]
    d[, .(contrast, category, ID, Description, setSize, NES, pvalue, p.adjust,
          core_enrichment)]
  }), fill = TRUE)), fill = TRUE)

if (!nrow(gs)) { message("no GSEA output found -- run 04_enrichment.R first"); quit(save = "no") }
message("  ", format(nrow(gs), big.mark = ","), " pathway x contrast rows")

################################################################################
# 2. Ruler verdict per pathway, from the FULL ruler tables
################################################################################
RUL     <- RULER_CONTRAST
RULER_P <- THRESH$ruler_p
sig <- gs[p.adjust < THRESH$gsea_fdr]

rfiles <- file.path(OUT, "enrichment_gsea", RUL, paste0("gsea_full_", CATS, ".csv"))
if (!all(file.exists(rfiles)))
  stop("full GSEA tables missing for the ruler (",
       paste(basename(rfiles[!file.exists(rfiles)]), collapse = ", "),
       "). gsea_<CAT>.csv only holds sets with p <= 0.05, so the ruler cannot be ",
       "read from it. Re-run 04_enrichment.R, which writes gsea_full_<CAT>.csv.")
rg <- rbindlist(lapply(rfiles, function(f)
  as.data.table(read.csv(f))[, .(ID, ruler_NES = NES, ruler_p = pvalue)]))
stopifnot(!anyDuplicated(rg$ID))
message("  ruler: ", nrow(rg), " sets tested; ", rg[ruler_p < RULER_P, .N],
        " at raw p < ", RULER_P)

res <- merge(sig[contrast != RUL], rg, by = "ID", all.x = TRUE)
res[, ruler_tested := !is.na(ruler_p)]
res[, in_ruler := ruler_tested & ruler_p < RULER_P]
res[, same_direction_as_ruler := in_ruler & sign(NES) == sign(ruler_NES)]

# verdict:
#   compartment          - in the ruler, same direction. The ruler explains it.
#   opposed-to-ruler     - in the ruler but the OTHER way. Not explained by compartment.
#   ruler-clean          - ruler raw p >= RULER_P.
#   not-tested-in-ruler  - the set was not tested in the ruler at all (should not occur)
res[, verdict := fifelse(!ruler_tested, "not-tested-in-ruler",
                 fifelse(same_direction_as_ruler, "compartment",
                 fifelse(in_ruler, "opposed-to-ruler", "ruler-clean")))]

fwrite(res[order(contrast, verdict, -abs(NES))], fig("80_pathway_ruler.tsv"), sep = "\t")

summ <- res[, .N, by = .(contrast, verdict)]
summ <- dcast(summ, contrast ~ verdict, value.var = "N", fill = 0)
summ <- merge(CONTRASTS[, .(contrast = id, label, class)], summ, by = "contrast")
message("\n  pathway-level ruler verdict:"); print(summ)
fwrite(summ, fig("80b_pathway_ruler_summary.tsv"), sep = "\t")

################################################################################
# 3. Pathways surviving the ruler in the DONOR contrasts
################################################################################
donor_con <- CONTRASTS[class == "context", id]
clean <- res[contrast %in% donor_con & verdict == "ruler-clean"]
clean_wide <- dcast(clean, ID + Description + category ~ contrast, value.var = "NES")
# require it in BOTH donor contrasts, same direction: present with the clone at
# 3.5yr AND at DDHN (large clone against the donor's small one), absent at day 30 (not tested; low clone burden)
both <- if (all(donor_con %in% names(clean_wide)))
  clean_wide[complete.cases(clean_wide[, ..donor_con])] else clean_wide[0]
if (nrow(both)) {
  both[, consistent := sign(get(donor_con[1])) == sign(get(donor_con[2]))]
  both <- both[consistent == TRUE]
  both <- merge(both, rg, by = "ID", all.x = TRUE)
}
fwrite(both, fig("81_clone_associated_pathways.tsv"), sep = "\t")
message("\n  pathways significant in BOTH donor contrasts, ruler-clean, ",
        "consistent direction: ", nrow(both))
if (nrow(both)) print(both[, .(Pathway = substr(gsub("^(HALLMARK_|GOBP_|KEGG_MEDICUS_)", "", Description), 1, 52),
                               category, DDHN = round(get(donor_con[1]), 2),
                               yr2 = round(get(donor_con[2]), 2),
                               ruler_NES = round(ruler_NES, 2), ruler_p = signif(ruler_p, 2))])

################################################################################
# 4. The leukaemia-relevant panel, reported whether significant or not
################################################################################
THEMES <- list(
  `p53 / DNA damage`     = c("HALLMARK_P53_PATHWAY", "HALLMARK_DNA_REPAIR",
                             "HALLMARK_APOPTOSIS", "HALLMARK_UV_RESPONSE_UP"),
  `Cell cycle`           = c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
                             "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
                             "HALLMARK_MITOTIC_SPINDLE"),
  `Cell killing / immune`= c("HALLMARK_INFLAMMATORY_RESPONSE",
                             "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                             "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                             "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                             "HALLMARK_ALLOGRAFT_REJECTION",
                             "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                             "HALLMARK_COMPLEMENT"),
  `Stemness / signalling`= c("HALLMARK_WNT_BETA_CATENIN_SIGNALING",
                             "HALLMARK_NOTCH_SIGNALING", "HALLMARK_HEDGEHOG_SIGNALING",
                             "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_MTORC1_SIGNALING",
                             "HALLMARK_KRAS_SIGNALING_UP"),
  `Erythroid content`    = c("HALLMARK_HEME_METABOLISM")
)
th <- rbindlist(lapply(names(THEMES), function(k) data.table(theme = k, ID = THEMES[[k]])))

panel <- merge(th, gs[category == "HALLMARK"], by = "ID", all.x = TRUE, allow.cartesian = TRUE)
panel[, sig := !is.na(p.adjust) & p.adjust < THRESH$gsea_fdr]
panel <- merge(panel, rg, by = "ID", all.x = TRUE)
fwrite(panel[order(theme, ID, contrast)], fig("83_leukaemia_themes.tsv"), sep = "\t")

pw <- dcast(panel[!is.na(contrast)], theme + ID ~ contrast, value.var = "NES")
message("\n  leukaemia-theme NES matrix (blank = raw p > 0.05 in that contrast):")
print(pw, nrows = 40)

################################################################################
# 5. Heatmap: NES across contrasts, with the ruler column marked
################################################################################
fig82 <- fig("82_pathway_ruler_heatmap.png"); unlink(fig82)
hm <- panel[!is.na(contrast) & !is.na(NES)]
if (nrow(hm)) {
  hm[, path := gsub("_", " ", gsub("^HALLMARK_", "", ID))]
  hm[, contrast_lab := factor(CONTRASTS[contrast, label], levels = CONTRASTS$label)]
  hm[, star := fifelse(sig, "*", "")]
  p <- ggplot(hm, aes(contrast_lab, path, fill = NES)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = star), size = 4.2, vjust = .78) +
    scale_fill_distiller(palette = "RdBu", limits = c(-3, 3), oob = scales::squish) +
    facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = NULL, y = NULL,
         title = "Leukaemia-relevant pathways across all contrasts",
         subtitle = str_wrap(paste0(
           "* FDR < 0.05. Blank = raw p > 0.05 in that comparison. The rightmost column is the ",
           "reference comparison (day 30 vs donor: no disease, small clone at both ends); for the verdict it is ",
           "tested at raw p < ", RULER_P, " from the full tables, so a blank there is not by ",
           "itself a pass."), 120)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 22, hjust = 1),
          strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.5),
          strip.placement = "outside",
          plot.subtitle = element_text(size = 7.8))
  ggsave(fig82, p, width = 10.5, height = max(6, uniqueN(hm$path) * 0.32 + 3.2), dpi = 300)
  message("  wrote 82_pathway_ruler_heatmap.png")
}

save(gs, rg, res, summ, both, panel, file = file.path(OUT, "pathway_ruler.RData"))
message("\n== done ==  ", run_stamp())
