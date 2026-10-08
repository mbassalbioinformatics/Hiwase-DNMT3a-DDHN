################################################################################
# 06_top_signatures.R  --  what actually comes top, with no gene selection by us
#
#   Rscript 03_rnaseq/06_top_signatures.R genes     (default)
#
# WHY THIS EXISTS
#   Curated themes and gene panels do not answer the obvious question:
#   if you take the enrichment output and simply rank it, what is at the top?
#   Presenting only pre-selected themes invites the reasonable objection that
#   the analysis found what it went looking for.
#
#   This script does no selection. It ranks every gene set in every contrast by
#   significance and reports the top of that ranking, alongside the reference
#   contrast verdict so a reader can see immediately which entries are tissue
#   compartment and which are not.
#
# WHAT IT PRODUCES
#   90_top_signatures_ranked.tsv      top 40 per contrast by FDR (all FDR < 0.05)
#   91_top_signatures_<contrast>.png  the 25 largest |NES| among those 40,
#                                     coloured by reference-contrast verdict
#   92_pathway_landscape.png          every Hallmark set, significant or not
#
# The reference verdict uses the same rule as 05_pathway_ruler.R: raw p <
# THRESH$ruler_p in the reference contrast, read from the FULL tables.
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(RColorBrewer); library(stringr)
})
SCRIPTS <- "03_rnaseq"
source(file.path(SCRIPTS, "00_config.R"))

args  <- commandArgs(trailingOnly = TRUE)
LEVEL <- if (length(args)) args[1] else "genes"
OUT <- if (LEVEL == "genes") PATHS$out_gene else PATHS$out_5utr
fig <- function(...) file.path(OUT, paste0(...))
CATS <- c("HALLMARK", "BP", "CC", "MF", "KEGG")

gs <- rbindlist(lapply(CONTRASTS$id, function(cid)
  rbindlist(lapply(CATS, function(ct) {
    f <- file.path(OUT, "enrichment_gsea", cid, paste0("gsea_", ct, ".csv"))
    if (!file.exists(f)) return(NULL)
    d <- as.data.table(read.csv(f)); if (!nrow(d)) return(NULL)
    d[, `:=`(contrast = cid, category = ct)]
    d[, .(contrast, category, ID, Description, setSize, NES, pvalue, p.adjust)]
  }), fill = TRUE)), fill = TRUE)
if (!nrow(gs)) { message("no GSEA output; run 04_enrichment.R first"); quit(save = "no") }

# ruler verdict, same rule as 05_pathway_ruler.R, from the full tables
RUL <- RULER_CONTRAST
rfiles <- file.path(OUT, "enrichment_gsea", RUL, paste0("gsea_full_", CATS, ".csv"))
if (!all(file.exists(rfiles)))
  stop("full GSEA tables missing for the reference contrast; re-run 04_enrichment.R")
ruler <- rbindlist(lapply(rfiles, function(f)
  as.data.table(read.csv(f))[, .(ID, ruler_NES = NES, ruler_p = pvalue)]))
ruler <- ruler[ruler_p < THRESH$ruler_p, .(ID, ruler_NES)]
gs <- merge(gs, ruler, by = "ID", all.x = TRUE)
gs[, verdict := fifelse(is.na(ruler_NES), "not seen in reference contrast",
                fifelse(sign(NES) == sign(ruler_NES),
                        "also in reference contrast (compartment)",
                        "opposite to reference contrast"))]
gs[contrast == RUL, verdict := "this IS the reference contrast"]

# For GO and KEGG the ID is an accession (GO:0035579, hsa03010) and the human
# readable name lives in Description; for MSigDB collections the ID *is* the
# name. Prefer Description whenever it is not just a repeat of the accession.
clean_name <- function(id, desc) {
  nm <- ifelse(!is.na(desc) & nzchar(desc) & desc != id, desc, id)
  gsub("_", " ", gsub("^(HALLMARK_|GOBP_|GOCC_|GOMF_|KEGG_MEDICUS_|REACTOME_)", "", nm))
}

top <- gs[p.adjust < 0.05][order(contrast, p.adjust)]
top <- top[, head(.SD, 40), by = contrast]
top[, Pathway := clean_name(ID, Description)]
fwrite(top[, .(contrast, category, Pathway, NES, pvalue, p.adjust, setSize, verdict)],
       fig("90_top_signatures_ranked.tsv"), sep = "\t")

message("\n== top signatures by significance, no pre-selection ==")
for (cid in CONTRASTS$id) {
  d <- top[contrast == cid]
  if (!nrow(d)) next
  cat("\n---", CONTRASTS[cid, label], "---\n")
  print(head(d[, .(Pathway = substr(Pathway, 1, 46), cat = category,
                   NES = round(NES, 2), FDR = signif(p.adjust, 2),
                   verdict = substr(verdict, 1, 34))], 12))
}

# ranked bar chart per contrast
for (cid in CONTRASTS$id) {
  f91 <- fig(sprintf("91_top_signatures_%s.png", cid)); unlink(f91)
  d <- top[contrast == cid][order(-abs(NES))][seq_len(min(25, .N))]
  if (!nrow(d)) next
  d[, Pathway := factor(Pathway, levels = rev(Pathway))]
  p <- ggplot(d, aes(NES, Pathway, fill = verdict)) +
    geom_col(width = .72) +
    geom_vline(xintercept = 0, colour = "grey40") +
    scale_fill_manual(values = c(
      "not seen in reference contrast"           = "#C0392B",
      "also in reference contrast (compartment)" = "#9AA5AE",
      "opposite to reference contrast"           = "#4C72B0",
      "this IS the reference contrast"           = "#6B7A88"), name = NULL) +
    labs(x = "Normalised enrichment score (positive = higher in the first group)",
         y = NULL, title = paste0("Top gene sets: ", CONTRASTS[cid, label]),
         subtitle = paste0("The 25 largest scores, in either direction, among the 40 lowest-FDR sets (all FDR < 0.05). No pre-selection.\n",
                           "Grey = the reference contrast returns it too (raw p < ", THRESH$ruler_p,
                           ", same direction), so tissue compartment explains it.")) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.subtitle = element_text(size = 7.8),
          axis.text.y = element_text(size = 8)) +
    guides(fill = guide_legend(nrow = 2))
  ggsave(f91, p, width = 10, height = max(5, nrow(d) * 0.30 + 2.8), dpi = 300)
}

################################################################################
# The full Hallmark landscape -- including sets that did NOT reach significance.
################################################################################
fig92 <- fig("92_pathway_landscape.png"); unlink(fig92)
hm <- gs[category == "HALLMARK"]
allsets <- unique(hm$ID)
grid <- CJ(ID = allsets, contrast = CONTRASTS$id, unique = TRUE)
hm <- merge(grid, hm, by = c("ID", "contrast"), all.x = TRUE)
hm[, Pathway := clean_name(ID, Description)]
hm[, Pathway := fifelse(is.na(Pathway), clean_name(ID, NA_character_), Pathway)]
hm[, sig := !is.na(p.adjust) & p.adjust < 0.05]
hm[, contrast_lab := factor(CONTRASTS[contrast, label], levels = CONTRASTS$label)]
# order by how strongly each set moves in the DDHN-vs-donor contrast
ord <- hm[contrast == "C4_DDHN_vs_donor", .(Pathway, o = fifelse(is.na(NES), 0, NES))]
ord <- unique(ord, by = "Pathway")
hm[, Pathway := factor(Pathway, levels = ord[order(o), Pathway])]

p <- ggplot(hm, aes(contrast_lab, Pathway, fill = NES)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = fifelse(sig, "*", "")), size = 3.6, vjust = .78) +
  scale_fill_distiller(palette = "RdBu", limits = c(-3, 3), oob = scales::squish,
                       na.value = "grey93") +
  labs(x = NULL, y = NULL,
       title = paste0("Hallmark gene sets reaching raw p <= 0.05 in at least one comparison (",
                      length(allsets), " of 50)"),
       subtitle = paste0("Grey = raw p > 0.05 in that comparison. * = FDR < 0.05.\n",
                         "Nothing that reached p <= 0.05 anywhere is filtered out, so sets that did NOT pass in a given comparison are visible alongside those that did.")) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 22, hjust = 1),
        axis.text.y = element_text(size = 7.4),
        plot.subtitle = element_text(size = 7.4))
ggsave(fig92, p, width = 10.5, height = max(8, uniqueN(hm$Pathway) * 0.21 + 3), dpi = 300)
message("\n  wrote 92_pathway_landscape.png (", uniqueN(hm$Pathway), " Hallmark sets)")

save(gs, top, file = file.path(OUT, "top_signatures.RData"))
message("\n== done ==  ", run_stamp())
