#!/usr/bin/env Rscript
###############################################################################
# 19_clinical_questions.R                                            2026-10-17
#
# The analyses the clinical team asked to see first (meeting and email of
# 16 October 2026), for Part A of the reports (12_build_reports.R):
#   the two patients' timelines: how big each clone was at each sample
#   donor at donation vs donor 2023      (the donor, five years apart)
#   donor at donation vs day 30          (engraftment)
#   day 30 -> 3.5 yr -> DDHN diagnosis   (the patient's course)
#   donor 2023 vs DDHN diagnosis         (donor and patient in the same months)
#   the programs that change, across all of them
#   TET2 (24.7 years from transplant to DDHN) against DNMT3A (4.5 years)
#
# Clone size: % of cells carrying a clone = 2 x the pooled VAF of its mutations
#   (17_genomewide_snvs.R, site counts; heterozygous mutations), or 2 x the clinical VAF.
# Growth: a clone making up a fraction f of cells has odds f / (1 - f). A clone with a
#   constant advantage changes its log-odds by the same amount s each year (logistic
#   growth, as most clonal-haematopoiesis clones do [Fabre 2022]); s = ln(odds2 / odds1) /
#   years, and while the clone is small it doubles every ln 2 / s years. The interval
#   carries the 95% binomial interval of the pooled VAF at the start (the end points are
#   large and tight). Between two samples s is the least average rate that gets the
#   clone from one size to the other: a clone that reached its final size earlier grew
#   faster.
# Programs from methylation: every genome tile (TILE_BP, config/settings.R) with data in
#   both samples is assigned to the nearest MANE Select TSS within 50 kb; a gene's score
#   is the mean change of its tiles (>= 5 tiles), so promoters and nearby enhancers count
#   together. Genes ranked by their change (largest gain first, largest loss last) go into
#   fgsea (multilevel, seed 1234) with Hallmark (50 sets), the bone-marrow cell types of
#   Hay et al. 2018 (MSigDB C8, 29 sets) and the 14 curated panels (resources/gene_panels.tsv).
#   Sign as everywhere in the reports: NES < 0 = methylation lost around the program's
#   genes, relative to other genes; NES > 0 = gained. The share of the program's tiles
#   that lost or gained >= CUT (20 points) gives the absolute direction.
#   Sensitivity: 20-kb window; promoter tiles only (within 2 kb of the TSS).
# Programs from RNA-seq: the same sets, fgsea on the edgeR statistic of each contrast
#   (signed square root of the likelihood ratio, BCV 0.4, as in 04_enrichment.R).
#
# Output: integration/clinical/
#   clone_timeline.tsv, clone_growth.tsv, comparisons_at_a_glance.tsv,
#   gene_scores_methylation.tsv.gz, programs.tsv, programs_sensitivity.tsv, programs_shown.tsv,
#   promoter_switches.tsv, programs_course.tsv, fig_timeline.png, fig_programs.png,
#   fig_programs_course.png, fig_tet2_vs_dnmt3a.png
###############################################################################
source("config/paths.R")
ID <- with(data.table::fread(SAMPLE_SHEET, colClasses = "character"), setNames(sample, key))   # role key -> sample ID (config/samples.tsv)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(msigdbr); library(fgsea)
})
OUT  <- file.path(INT, "clinical"); dir.create(OUT, showWarnings = FALSE)
rd   <- function(f, ...) if (grepl("\\.gz$", f)) fread(cmd = sprintf("gzip -dc '%s'", f), ...) else fread(f, ...)
SEED <- 1234L
WINDOW <- 50e3; MIN_TILES <- 5L; BIG <- CUT; SWITCH <- CUT

## ---- samples and comparisons ------------------------------------------------------------
smp <- fread(SAMPLE_SHEET, colClasses = "character")[methylation == "analysed"]
smp[, date := as.Date(sample_date)]
HSCT <- unique(smp[nzchar(hsct_date), .(case, hsct_date)])
HSCT <- setNames(as.Date(HSCT$hsct_date), HSCT$case)                        # transplant date per case (config/samples.tsv)
smp[, years := as.numeric(date - HSCT[case]) / 365.25]
stopifnot(!anyNA(smp$date))
# the clinical team's short names (email of 16 Oct 2026), on two lines for the figures
SHORT <- setNames(sub(" \\((DNMT3A|TET2)\\)$", "\n(\\1)", smp$plot_label), smp$sample)
# methylation comparisons in the order of the clinical team's questions; rna = the matching RNA-seq contrast
CMP <- data.table(
  id  = c(paste0(ID[["D23"]], "_vs_", ID[["D18"]]), paste0(ID[["D30"]], "_vs_", ID[["D18"]]), paste0(ID[["Y35"]], "_vs_", ID[["D30"]]), paste0(ID[["DX"]], "_vs_", ID[["Y35"]]), paste0(ID[["DX"]], "_vs_", ID[["D30"]]),
          paste0(ID[["DX"]], "_vs_", ID[["D23"]]), paste0(ID[["D30"]], "_vs_", ID[["D23"]]), paste0(ID[["T2DX"]], "_vs_", ID[["T2A"]])),
  q   = c("Q1", "Q2", "Q3", "Q3", "Q3", "Q4", "ref", "TET2"),
  lab = c("Donor: 5 years post donation vs donation", "Recipient 1 month vs donor at donation", "Recipient: 3.5 years vs 1 month",
          "Recipient: DDHN diagnosis vs 3.5 years", "Recipient: DDHN diagnosis vs 1 month", "Recipient DDHN diagnosis vs donor 5 years post donation",
          "Reference: recipient 1 month vs donor 5 years post donation", "TET2: recipient DDHN diagnosis vs 1 year"),
  rna = c(NA, NA, "C1_clone", "C2_transform", "C3_cumulative", "C4_DDHN_vs_donor", "C6_D30_vs_donor", NA))
CMP[, c("g1", "g2") := tstrsplit(id, "_vs_")]
CMP[, names := paste(smp[match(g1, sample), role], "vs", smp[match(g2, sample), role])]

## ---- 1. clone sizes over time -------------------------------------------------------------
pool <- fread(file.path(INT, "variants", "pooled.tsv"), colClasses = list(character = "sample"))
pool <- pool[grepl("^site counts", source)]
stopifnot(nrow(pool) > 0)
cells <- function(st, s) pool[set == st & sample == s, .(pct = 200 * vaf, lo = 200 * ci_low, hi = 200 * ci_high)]
clin  <- function(s) min(2 * as.numeric(smp[sample == s, vaf_pct]), 100)   # clinical VAF (%) x 2 = % of cells
mk <- function(what, st, ss) rbindlist(lapply(ss, function(s) cbind(data.table(clone = what, sample = s, years = smp$years[smp$sample == s]), cells(st, s))))
DN <- c(ID[["D18"]], ID[["D30"]], ID[["Y35"]], ID[["DX"]])
tl <- rbind(
  # the branch that became the DDHN: mutations of the 3.5-year and DDHN clone absent from the donor's blood in 2023
  mk("DNMT3A: the DDHN branch", "DNMT3A_clone_not_in_donor", DN),
  # the donor's own clone: the clone's oldest mutations, which the donor's 2023 blood also carries
  mk("DNMT3A: the donor's clone (oldest mutations)", "DNMT3A_clone_in_donor", DN),
  mk("TET2: the DDHN clone", "TET2_clone", c(ID[["T2A"]], ID[["T2DX"]])))
tl[, source := "genome-wide mutations, pooled"]
# clinical VAFs (x 2 = % of cells): the donor's own clone in 2023 cannot be measured from the pooled
# sites (they were chosen by its reads there), so it is taken from the clinical VAF
cv <- smp[nzchar(vaf_pct), .(sample, years, pct = sapply(sample, clin),
                              clone = fifelse(case == "TET2", "TET2: clinical VAF x 2", "DNMT3A: clinical VAF x 2"))]
tl <- rbind(tl, cv[, .(clone, sample, years, pct, lo = NA_real_, hi = NA_real_, source = "clinical VAF")], fill = TRUE)
ix <- match(tl$sample, smp$sample)
tl[, `:=`(name = smp$key[ix], short = SHORT[sample], date = format(smp$date[ix], "%d %b %Y"))]
fwrite(tl, file.path(OUT, "clone_timeline.tsv"), sep = "\t")

## ---- 2. growth rates ----------------------------------------------------------------------------
odds <- function(p) (p / 100) / (1 - p / 100)
grow <- function(what, a, b, f1, f1lo, f1hi, f2) {
  y <- smp[sample == b, years] - smp[sample == a, years]
  s <- log(odds(f2) / odds(f1)) / y
  data.table(clone = what, from = a, to = b, from_name = smp[sample == a, role], to_name = smp[sample == b, role],
             years = y, pct_from = f1, pct_to = f2, rate = s,
             rate_lo = log(odds(f2) / odds(f1hi)) / y, rate_hi = log(odds(f2) / odds(f1lo)) / y)
}
g_br <- tl[clone == "DNMT3A: the DDHN branch"]; g_tr <- tl[clone == "DNMT3A: the donor's clone (oldest mutations)"]
g_t2 <- tl[clone == "TET2: the DDHN clone"]
gr <- rbind(
  grow("DNMT3A DDHN branch, in the patient", ID[["D30"]], ID[["Y35"]], g_br[sample == ID[["D30"]], pct], g_br[sample == ID[["D30"]], lo],
       g_br[sample == ID[["D30"]], hi], g_br[sample == ID[["Y35"]], pct]),
  grow("DNMT3A clone, in the donor", ID[["D18"]], ID[["D23"]], g_tr[sample == ID[["D18"]], pct], g_tr[sample == ID[["D18"]], lo],
       g_tr[sample == ID[["D18"]], hi], clin(ID[["D23"]])),
  grow("TET2 DDHN clone, in the patient", ID[["T2A"]], ID[["T2DX"]], g_t2[sample == ID[["T2A"]], pct], g_t2[sample == ID[["T2A"]], lo],
       g_t2[sample == ID[["T2A"]], hi], g_t2[sample == ID[["T2DX"]], pct]))
gr[, `:=`(doubling_years = log(2) / rate, doubling_lo = log(2) / rate_hi, doubling_hi = log(2) / rate_lo,
          fold_per_year = exp(rate))]
fwrite(gr, file.path(OUT, "clone_growth.tsv"), sep = "\t")
print(gr[, .(clone, years = round(years, 2), pct_from = round(pct_from, 2), pct_to = round(pct_to, 1), rate = round(rate, 2),
             rate_lo = round(rate_lo, 2), rate_hi = round(rate_hi, 2), doubling_months = round(12 * doubling_years, 1))])

## ---- 3. programs: gene sets ----------------------------------------------------------------------
H  <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))[, .(set = gs_name, gene = gene_symbol)]
C8 <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C8"))[grepl("^HAY_BONE_MARROW", gs_name), .(set = gs_name, gene = gene_symbol)]
PN <- fread(file.path(RES, "gene_panels.tsv"))[, .(set = paste0("PANEL_", panel), gene)]
SETS <- unique(rbind(H, C8, PN)); SL <- split(SETS$gene, SETS$set)
msig_version <- tryCatch(as.character(packageVersion("msigdbdf")), error = function(e) as.character(packageVersion("msigdbr")))

# FIGS_ONLY=1 redraws the figures from the saved tables without re-running the enrichment
FIGS_ONLY <- Sys.getenv("FIGS_ONLY") == "1" && file.exists(file.path(OUT, "programs.tsv"))
ed <- rd(file.path(RNA, "05_processed_counts_genes", "09_edgeR_all_bcv.tsv.gz"))[bcv == 0.4 & is.finite(t) & nzchar(symbol)]
if (!FIGS_ONLY) {
  ## ---- 4. methylation around each gene ---------------------------------------------------------------
  ti <- rd(TILE_FILE)
  ti[, `:=`(chr = sub(":.*", "", id), mid = as.integer(sub(".*:", "", id)) + TILE_BP %/% 2L - 1L)]
  gn <- rd(file.path(INT, "region_methylation", "mane_select_genes.tsv.gz"))[!duplicated(gene), .(gene, chr, tss)]
  setkey(gn, chr, tss)
  near <- gn[ti[, .(id, chr, tss = mid)], on = .(chr, tss), roll = "nearest", .(id, gene, dist = abs(x.tss - i.tss))]
  ti <- merge(ti, near[!is.na(gene) & dist <= WINDOW], by = "id")
  cat(sprintf("tiles within %d kb of a MANE TSS: %s (genes %s)\n", WINDOW / 1e3, format(nrow(ti), big.mark = ","), uniqueN(ti$gene)))
  # a gene's score: the mean change of its tiles minus the median change of genes whose surroundings are
  # methylated to the same level (40 bins of the two samples' mean methylation). A library that reads every
  # CpG a little higher (the 3.5-year library: mean CpG methylation 80.4% against 78.9-79.5% for the others)
  # shifts methylated regions more than unmethylated ones; without this, that shift alone would rank programs
  # by how methylated their surroundings are. raw = the unadjusted mean change, for the absolute direction.
  score_genes <- function(g1, g2, max_dist = WINDOW, min_tiles = MIN_TILES) {
    d <- ti[dist <= max_dist, .(gene, delta = get(g1) - get(g2), level = (get(g1) + get(g2)) / 2)][!is.na(delta)]
    g <- d[, .(raw = mean(delta), level = mean(level), tiles = .N, lost = sum(delta <= -BIG), gained = sum(delta >= BIG)), by = gene][tiles >= min_tiles]
    g[, bin := cut(level, unique(quantile(level, seq(0, 1, length.out = 41))), include.lowest = TRUE)]
    g[, score := raw - median(raw), by = bin][, bin := NULL]
  }
  run_fgsea <- function(stat) {
    set.seed(SEED)
    f <- suppressWarnings(fgseaMultilevel(SL, stat, minSize = 10, maxSize = 2000, eps = 0, nPermSimple = 10000))
    f[, `:=`(lead = vapply(leadingEdge, function(x) paste(head(x, 15), collapse = ", "), ""),
             leading_edge = vapply(leadingEdge, paste, "", collapse = ", "))]          # lead: the first 15, for the reports
    f[, .(set = pathway, size, NES, pval, padj, lead, leading_edge)]
  }
  gscore <- list(); prog <- list(); sens <- list()
  for (i in seq_len(nrow(CMP))) {
    c1 <- CMP[i]
    g <- score_genes(c1$g1, c1$g2); g[, id := c1$id]; gscore[[c1$id]] <- g
    f <- run_fgsea(sort(setNames(g$score, g$gene), decreasing = TRUE))   # NES < 0: lost, NES > 0: gained
    # absolute direction: the share of the program's tiles (within the window) that lost or gained >= CUT
    # and the median raw change of the program's genes against the median gene, in points
    tt <- merge(SETS, g[, .(gene, tiles, lost, gained, raw)], by = "gene")[
      , .(tiles = sum(tiles), lost = sum(lost), gained = sum(gained), median_change = 100 * median(raw)), by = set]
    bg <- g[, .(lost = sum(lost) / sum(tiles), gained = sum(gained) / sum(tiles), median_change = 100 * median(raw))]
    f <- merge(f, tt, by = "set", all.x = TRUE)
    f[, `:=`(pct_lost = 100 * lost / tiles, pct_gained = 100 * gained / tiles, bg_lost = 100 * bg$lost, bg_gained = 100 * bg$gained,
             bg_median_change = bg$median_change)]
    prog[[paste0("M", i)]] <- f[, `:=`(layer = "methylation", id = c1$id)]
    for (v in list(list("window 20 kb", 20e3), list("promoter only (2 kb)", 2e3))) {
      g2 <- score_genes(c1$g1, c1$g2, v[[2]], if (v[[2]] <= 2e3) 1L else MIN_TILES)
      f2 <- run_fgsea(sort(setNames(g2$score, g2$gene), decreasing = TRUE))
      sens[[paste(i, v[[1]])]] <- merge(f[, .(set, NES_main = NES, padj_main = padj)], f2[, .(set, NES_alt = NES, padj_alt = padj)], by = "set")[
        , `:=`(id = c1$id, variant = v[[1]])]
    }
    cat(sprintf("%-22s genes %5d  sets FDR<0.05: %2d (lost %2d, gained %2d)\n", c1$id, nrow(g), sum(f$padj < 0.05, na.rm = TRUE),
                sum(f$padj < 0.05 & f$NES < 0, na.rm = TRUE), sum(f$padj < 0.05 & f$NES > 0, na.rm = TRUE)))
  }
  fwrite(rbindlist(gscore), file.path(OUT, "gene_scores_methylation.tsv.gz"), sep = "\t")

  ## ---- 5. programs from RNA-seq -------------------------------------------------------------------------
  for (cc in CMP[!is.na(rna), unique(rna)]) {
    d <- ed[contrast == cc][order(-abs(t))][!duplicated(symbol)]
    f <- run_fgsea(sort(setNames(d$t, d$symbol), decreasing = TRUE))
    prog[[cc]] <- f[, `:=`(layer = "RNA-seq", id = cc)]
    cat(sprintf("%-22s genes %5d  sets FDR<0.05: %2d (up %2d, down %2d)\n", cc, nrow(d), sum(f$padj < 0.05, na.rm = TRUE),
                sum(f$padj < 0.05 & f$NES > 0, na.rm = TRUE), sum(f$padj < 0.05 & f$NES < 0, na.rm = TRUE)))
  }
  prog <- rbindlist(prog, fill = TRUE)
  SENS <- rbindlist(sens)
  fwrite(prog, file.path(OUT, "programs.tsv"), sep = "\t")
  fwrite(SENS, file.path(OUT, "programs_sensitivity.tsv"), sep = "\t")
  print(SENS[, .(r = round(cor(NES_main, NES_alt, use = "complete.obs"), 2),
                 agree_sig = sum(padj_main < 0.05 & padj_alt < 0.05 & sign(NES_main) == sign(NES_alt), na.rm = TRUE),
                 sig_main = sum(padj_main < 0.05, na.rm = TRUE)), by = .(variant, id)])
} else {
  prog <- fread(file.path(OUT, "programs.tsv")); SENS <- fread(file.path(OUT, "programs_sensitivity.tsv"))
  message("FIGS_ONLY: programs read from ", OUT)
}

## ---- 6. promoter switches with gene names -------------------------------------------------------------
pm <- rd(file.path(INT, "region_methylation", "promoter_methylation.tsv.gz"))
sw <- rbindlist(lapply(seq_len(nrow(CMP)), function(i) {
  c1 <- CMP[i]; d <- data.table(gene = pm$id, m1 = pm[[c1$g1]], m2 = pm[[c1$g2]])[!is.na(m1) & !is.na(m2)]
  d[, delta := m1 - m2][abs(delta) >= SWITCH][, id := c1$id][order(-abs(delta))]
}))
fwrite(sw, file.path(OUT, "promoter_switches.tsv"), sep = "\t")

## ---- 7. each comparison at a glance -----------------------------------------------------------------
dirc <- fread(file.path(INT, "global_direction", "direction_by_context.tsv"))
cnt  <- fread(file.path(INT, "dmr", "dmr_counts.tsv"))
de   <- ed[, .(up = sum(sig & logFC > 0, na.rm = TRUE), down = sum(sig & logFC < 0, na.rm = TRUE)), by = contrast]
ctx  <- function(i, k, w) dirc[id == i & context == k, get(w)]
glance <- CMP[, .(id, q, lab, names, rna,
  tiles_lost = 100 * mapply(ctx, id, "all", "lost"), tiles_gained = 100 * mapply(ctx, id, "all", "gained"),
  dels_lost = 100 * mapply(ctx, id, "dELS", "lost"), dels_gained = 100 * mapply(ctx, id, "dELS", "gained"),
  cgi_lost = 100 * mapply(ctx, id, "CGI", "lost"), cgi_gained = 100 * mapply(ctx, id, "CGI", "gained"),
  sw_lost = vapply(id, function(x) sw[id == x & delta < 0, .N], 0L), sw_gained = vapply(id, function(x) sw[id == x & delta > 0, .N], 0L))]
cnt <- cnt[set == "metilene" & context == "all"]
stopifnot(all(CMP$id %in% cnt$id))
glance[, `:=`(dmr_lost = cnt$lost[match(glance$id, cnt$id)], dmr_gained = cnt$gained[match(glance$id, cnt$id)])]
glance[, `:=`(rna_up = de$up[match(glance$rna, de$contrast)], rna_down = de$down[match(glance$rna, de$contrast)])]
fwrite(glance, file.path(OUT, "comparisons_at_a_glance.tsv"), sep = "\t")
print(glance[, .(lab, tiles_lost = round(tiles_lost, 2), tiles_gained = round(tiles_gained, 2), dmr_lost, dmr_gained, sw_lost, sw_gained, rna_up, rna_down)])

## ---- figure: timeline -----------------------------------------------------------------------------------
COL <- c("DNMT3A: the DDHN branch" = "#6A51A3", "DNMT3A: the donor's clone (oldest mutations)" = "#B39DDB",
         "TET2: the DDHN clone" = "#D95F02")
# legend text: each clone with its copy-number status (Section 5: the KMT2A duplication is the only acquired change)
COL_LAB <- c("DNMT3A: the DDHN branch" = "DNMT3A: the DDHN branch, with the KMT2A duplication",
             "DNMT3A: the donor's clone (oldest mutations)" = "DNMT3A: the donor's clone (oldest mutations), without it",
             "TET2: the DDHN clone" = "TET2: the DDHN clone, no copy-number change")
pt <- copy(tl)[, `:=`(y = pct, ylo = lo, yhi = pmin(hi, 100))]
pt[source == "clinical VAF", clone := fifelse(grepl("TET2", clone), "TET2: the DDHN clone", "DNMT3A: the DDHN branch")]
pt[source == "clinical VAF" & sample == ID[["D23"]], clone := "DNMT3A: the donor's clone (oldest mutations)"]
pt[, meas := factor(source, levels = c("genome-wide mutations, pooled", "clinical VAF"),
                    labels = c("genome-wide clone mutations, pooled (2 x VAF, 95% interval)", "clinical VAF x 2"))]
SHP <- c(16, 4); names(SHP) <- levels(pt$meas)
# panel A: both patients on one time axis (years after their own transplant)
pa <- pt[clone %in% c("DNMT3A: the DDHN branch", "TET2: the DDHN clone") & !sample %in% c(ID[["D18"]], ID[["D23"]])]
dx <- smp[grepl("DDHN", role), .(years, case, txt = sprintf("%s DDHN diagnosis\n%.1f years", case, years))]
dx[, col := COL[fifelse(case == "TET2", "TET2: the DDHN clone", "DNMT3A: the DDHN branch")]]
lab_a <- gr[clone != "DNMT3A clone, in the donor",
  .(clone = fifelse(grepl("TET2", clone), "TET2: the DDHN clone", "DNMT3A: the DDHN branch"),
    txt = sprintf("%s\ndoubled about every %s", fifelse(grepl("TET2", clone), "TET2 DDHN clone", "DNMT3A DDHN branch"),
                  fifelse(doubling_years < 1, sprintf("%.0f months", 12 * doubling_years), sprintf("%.1f years", doubling_years))))]
lab_a[, `:=`(x = fifelse(grepl("TET2", clone), 11.5, 5.6), y = fifelse(grepl("TET2", clone), 1.6, 0.9))]
pl_a <- pa[source != "clinical VAF", .(years = mean(years), y = max(y), short = short[1], name = name[1], clone = clone[1]), by = sample]
# label positions placed by hand so that nothing overlaps: x, y of the text and its alignment
LP <- data.table(sample = c(ID[["D30"]], ID[["Y35"]], ID[["DX"]], ID[["T2A"]], ID[["T2DX"]]),
                 lx = c(0.6, 2.0, 6.4, 2.7, 21.6), ly = c(0.075, 165, 105, 0.165, 165), hj = c(0, 0.5, 0, 0, 0.5))
pl_a <- merge(pl_a, LP, by = "sample")
dx[, `:=`(lx = fifelse(case == "TET2", years - 0.25, years + 0.25), hj = fifelse(case == "TET2", 1, 0))]
p1 <- ggplot(pa[source != "clinical VAF"], aes(years, y, colour = clone)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = "dashed") +
  geom_vline(data = dx, aes(xintercept = years), colour = dx$col, linetype = "dotted") +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = ylo, ymax = yhi), width = 0.25, linewidth = 0.6, na.rm = TRUE) +
  geom_point(data = pa, aes(shape = meas), size = 3, stroke = 1.2) +
  geom_segment(data = pl_a, aes(x = years, y = y, xend = lx, yend = ly), colour = "grey55", linewidth = 0.4) +
  geom_label(data = pl_a, aes(lx, ly, label = short, hjust = hj), size = 3.1, lineheight = 0.9,
             linewidth = 0, fill = "white", show.legend = FALSE) +
  geom_label(data = lab_a, aes(x, y, label = txt, colour = clone), size = 3.5, hjust = 0, linewidth = 0.3, show.legend = FALSE, inherit.aes = FALSE) +
  annotate("text", x = 0.15, y = 330, label = "transplant", hjust = 0, size = 3.2, colour = "grey35") +
  geom_text(data = dx, aes(x = lx, y = 330, label = txt, hjust = hj), size = 3.1, lineheight = 0.9, colour = dx$col, inherit.aes = FALSE) +
  scale_colour_manual(values = COL, labels = COL_LAB, name = NULL) + scale_shape_manual(values = SHP, name = NULL) +
  scale_y_log10(breaks = c(0.1, 0.3, 1, 3, 10, 30, 100), labels = function(v) paste0(v, "%"), limits = c(0.05, 520)) +
  scale_x_continuous(breaks = seq(0, 25, 5), limits = c(-0.6, 26)) +
  labs(x = "years after the patient's own transplant", y = "share of cells carrying the clone (log scale)",
       title = "A. The two DDHN clones after transplant: a similar start, a very different speed") +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", legend.box = "vertical", panel.grid.minor = element_blank())
# panel B: the DNMT3A patient and the donor in detail
pb <- pt[grepl("^DNMT3A", clone)]
in_patient <- pb[source != "clinical VAF" & !sample %in% c(ID[["D23"]])]
donor_line <- rbind(pb[clone == "DNMT3A: the donor's clone (oldest mutations)" & sample == ID[["D18"]] & source != "clinical VAF"],
                    pb[source == "clinical VAF" & sample == ID[["D23"]]])
lab_b <- unique(pb[, .(years, short, name, sample)])[order(years)]
lab_b[, y := c(0.085, 0.045)[(seq_len(.N) - 1) %% 2 + 1]]
# each label sits under its sample; the donor 5 years post donation is moved right, clear of the DDHN leader
lab_b[, lx := years + fifelse(sample == ID[["D23"]], 0.22, 0)]
# dashed leaders from the labels up to the lowest point (or interval) of their sample; the 1-month label
# needs none, it sits on the transplant line next to its point
lab_b <- merge(lab_b, pb[, .(low = min(pmin(y, ylo, na.rm = TRUE))), by = sample], by = "sample")
p2 <- ggplot(pb, aes(years, y, colour = clone)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = "dashed") +
  geom_segment(data = lab_b[sample != ID[["D30"]]], aes(x = years, xend = years, y = y * 1.3, yend = low / 1.15), inherit.aes = FALSE,
               colour = "grey60", linetype = "dashed", linewidth = 0.4) +
  geom_line(data = in_patient, linewidth = 0.9) +
  geom_line(data = donor_line, linewidth = 0.9, linetype = "22") +
  geom_errorbar(aes(ymin = ylo, ymax = yhi), width = 0.08, linewidth = 0.6, na.rm = TRUE) +
  geom_point(aes(shape = meas), size = 3, stroke = 1.2) +
  # white-backed labels, so the transplant line does not run through the 1-month name
  geom_label(data = lab_b, aes(lx, y, label = short), inherit.aes = FALSE, size = 3, lineheight = 0.9, colour = "grey20",
             fill = "white", border.colour = NA, label.padding = unit(0.08, "lines")) +
  annotate("text", x = 1.9, y = 0.14, label = "dashed purple line: the donor's own clone, in the donor", size = 3.2, colour = "#8E7CC3") +
  scale_shape_manual(values = SHP, name = NULL) + scale_colour_manual(values = COL, labels = COL_LAB, name = NULL) +
  scale_y_log10(breaks = c(0.1, 0.3, 1, 3, 10, 30, 100), labels = function(v) paste0(v, "%"), limits = c(0.03, 160)) +
  scale_x_continuous(breaks = seq(-0.5, 5, 0.5), limits = c(-0.9, 5.8)) +
  labs(x = sprintf("years after the patient's transplant (%s); the donor's samples at their own dates",
                 format(HSCT[["DNMT3A"]], "%d %b %Y")), y = "share of cells carrying the clone (log scale)",
       title = "B. DNMT3A patient and donor: the donor's clone, and the branch of it that became the DDHN") +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", legend.box = "vertical", panel.grid.minor = element_blank())
ggsave(file.path(OUT, "fig_timeline.png"), p1 / p2 + plot_layout(heights = c(1, 1)), width = 13, height = 14.5, dpi = 170)

## ---- figure: programs ---------------------------------------------------------------------------------
SHOW <- data.table(set = c(
  "PANEL_HOX_MEIS_self_renewal", "PANEL_HSC_stemness", "PANEL_leukaemic_stem_cell", "HAY_BONE_MARROW_CD34_POS_HSC",
  "HAY_BONE_MARROW_CD34_POS_MULTILIN", "HAY_BONE_MARROW_CD34_POS_LMPP", "PANEL_progenitors",
  "HAY_BONE_MARROW_CD34_POS_GRAN", "HAY_BONE_MARROW_IMMATURE_NEUTROPHIL", "HAY_BONE_MARROW_NEUTROPHIL", "HAY_BONE_MARROW_MONOCYTE",
  "PANEL_granulocyte_monocyte", "HAY_BONE_MARROW_DENDRITIC_CELL",
  "HAY_BONE_MARROW_CD34_POS_ERP", "HAY_BONE_MARROW_EARLY_ERYTHROBLAST", "HAY_BONE_MARROW_ERYTHROBLAST", "PANEL_erythroid",
  "HAY_BONE_MARROW_CD34_POS_MKP", "PANEL_megakaryocyte",
  "HAY_BONE_MARROW_NAIVE_T_CELL", "HAY_BONE_MARROW_CD8_T_CELL", "HAY_BONE_MARROW_NK_CELLS", "HAY_BONE_MARROW_FOLLICULAR_B_CELL",
  "HAY_BONE_MARROW_PLASMA_CELL", "PANEL_lymphoid",
  "HALLMARK_MYC_TARGETS_V1", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_P53_PATHWAY", "PANEL_epigenetic_regulators", "PANEL_imprinted", "PANEL_splicing"),
  group = rep(c("Stem cells and\nself-renewal", "Myeloid", "Erythroid and\nmegakaryocyte", "Lymphoid", "Proliferation", "Inflammation\nand signalling", "Other"),
              c(7, 6, 6, 6, 3, 5, 5)))
nice <- function(s) {
  x <- sub("^HAY_BONE_MARROW_", "", s); x <- sub("^HALLMARK_", "", x); x <- sub("^PANEL_", "", x)
  x <- gsub("_", " ", x); x <- sub("CD34 POS ", "CD34+ ", x)
  x <- paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
  x <- sub("^Hox meis", "HOX/MEIS", x); x <- sub("^Hsc", "HSC", x); x <- sub("Cd34\\+ hsc", "CD34+ HSC", x); x <- sub("^Cd34\\+", "CD34+", x)
  x <- sub("lmpp", "LMPP", x); x <- sub("mkp", "MkP", x); x <- sub("erp", "ErP", x); x <- sub("gran", "granulocyte progenitor", x)
  x <- sub("multilin", "multilineage", x); x <- sub("Myc", "MYC", x); x <- sub("E2f", "E2F", x); x <- sub("G2m", "G2M", x)
  x <- sub("Tnfa", "TNFa", x); x <- sub("nfkb", "NF-kB", x); x <- sub("Wnt beta", "WNT/beta", x); x <- sub("P53", "p53", x)
  x <- sub("Nk cells", "NK cells", x); x <- sub("Cd8", "CD8", x); x <- sub(" t cell", " T cell", x); x <- sub(" b cell", " B cell", x)
  src <- fifelse(grepl("^HAY", s), " [Hay]", fifelse(grepl("^HALLMARK", s), " [H]", " [panel]"))
  paste0(x, src)
}
fwrite(SHOW[, .(set, group = gsub("\n", " ", group), label = nice(set))], file.path(OUT, "programs_shown.tsv"), sep = "\t")
pp <- merge(prog, SHOW, by = "set")
pp[, col := fifelse(layer == "methylation", CMP[match(pp$id, CMP$id), lab], CMP[match(pp$id, CMP$rna), lab])]
pp[, col := factor(col, levels = CMP$lab)]
pp[, `:=`(row = factor(nice(set), levels = rev(nice(SHOW$set))), group = factor(group, levels = unique(SHOW$group)),
          layer = factor(layer, levels = c("methylation", "RNA-seq"),
                         labels = c("Methylation around the genes\n(blue = lost, red = gained)", "Expression\n(green = down, purple = up)")))]
pp[, `:=`(sig = fifelse(!is.na(padj) & padj < 0.05, "*", ""), nes = pmax(pmin(NES, 2.8), -2.8))]
# two colour scales, so that methylation and expression cannot be mistaken for each other
p3 <- ggplot(pp, aes(col, row)) +
  geom_tile(data = pp[layer == levels(layer)[1]], aes(fill = nes), colour = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-2.8, 2.8),
                       name = "Methylation NES\n(< 0 lost,\n> 0 gained)", guide = guide_colourbar(order = 1)) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = pp[layer == levels(layer)[2]], aes(fill = nes), colour = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#1B7837", mid = "white", high = "#762A83", midpoint = 0, limits = c(-2.8, 2.8),
                       name = "Expression NES\n(< 0 down,\n> 0 up)", guide = guide_colourbar(order = 2)) +
  geom_text(aes(label = sig), size = 5, vjust = 0.78) +
  facet_grid(group ~ layer, scales = "free", space = "free") +
  labs(x = NULL, y = NULL, title = "Programs that change in each comparison: methylation and expression",
       subtitle = "NES from gene-set enrichment; * FDR < 0.05.\nRNA-seq exists only for the recipient at 1 month, 3.5 years and DDHN diagnosis and the donor 5 years post donation") +
  theme_bw(base_size = 11.5) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), strip.text.y = element_text(angle = 0, size = 10),
        panel.grid = element_blank(), plot.subtitle = element_text(size = 10), legend.position = "right")
ggsave(file.path(OUT, "fig_programs.png"), p3, width = 13.5, height = 14.5, dpi = 170)

## ---- figure: programs across the recipient's course (the three-way table of Part A as lines) -------------
# Each program as a line over the recipient at 1 month, 3.5 years and DDHN diagnosis, 1 month = 0.
# Methylation: the mean gene score of the program's genes (points; the level-adjusted change of the tiles
# within 50 kb that the enrichment ranks) for 3.5 years vs 1 month and DDHN diagnosis vs 1 month, minus the
# mean of all genes. Expression: the median log2 fold-change of the program's expressed genes (edgeR, BCV
# 0.4, logCPM > 0) for the same two contrasts, minus that of all expressed genes. Both relative to all genes,
# as the enrichment is. Filled points: that comparison's enrichment has FDR < 0.05. Programs: those
# of the figure above with FDR < 0.05 in at least one of the six comparisons of the course.
C_M <- c(`3.5 years` = paste0(ID[["Y35"]], "_vs_", ID[["D30"]]), `DDHN diagnosis` = paste0(ID[["DX"]], "_vs_", ID[["D30"]]), step = paste0(ID[["DX"]], "_vs_", ID[["Y35"]]))
C_R <- c(`3.5 years` = "C1_clone", `DDHN diagnosis` = "C3_cumulative", step = "C2_transform")
sig_any <- prog[(layer == "methylation" & id %in% C_M) | (layer == "RNA-seq" & id %in% C_R), .(any = any(padj < 0.05, na.rm = TRUE)), by = set][any == TRUE, set]
crs_sets <- intersect(SHOW$set, sig_any)
gsc <- rd(file.path(OUT, "gene_scores_methylation.tsv.gz"))[id %in% C_M[1:2]]
# relative to all genes, as the enrichment is: the program's mean gene score minus the mean of all genes (methylation), the
# program's median fold-change minus that of all expressed genes (expression)
mv <- merge(SETS[set %in% crs_sets], gsc[, .(gene, id, score)], by = "gene", allow.cartesian = TRUE)[, .(v = 100 * mean(score), genes = .N), by = .(set, id)]
mv <- merge(mv, gsc[, .(bg = 100 * mean(score)), by = id], by = "id")[, v := v - bg][, bg := NULL]
mv[, step := names(C_M)[match(id, C_M)]]
ev <- ed[contrast %in% C_R[1:2] & logCPM > 0][order(-abs(t))][!duplicated(paste(contrast, symbol))]
rv <- merge(SETS[set %in% crs_sets], ev[, .(gene = symbol, id = contrast, logFC)], by = "gene", allow.cartesian = TRUE)[, .(v = median(logFC), genes = .N), by = .(set, id)]
rv <- merge(rv, ev[, .(bg = median(logFC)), by = .(id = contrast)], by = "id")[, v := v - bg][, bg := NULL]
rv[, step := names(C_R)[match(id, C_R)]]
padj_of <- function(ly, s, i) prog[layer == ly][match(paste(s, i), paste(set, id)), padj]
mv[, padj := padj_of("methylation", set, id)]; rv[, padj := padj_of("RNA-seq", set, id)]
course <- rbind(cbind(rbind(data.table(set = crs_sets, step = "1 month", v = 0, padj = NA_real_), mv[, .(set, step, v, padj)]), layer = "methylation"),
                cbind(rbind(data.table(set = crs_sets, step = "1 month", v = 0, padj = NA_real_), rv[, .(set, step, v, padj)]), layer = "expression"))
course <- merge(course, SHOW, by = "set")
course[, `:=`(x = match(step, c("1 month", "3.5 years", "DDHN diagnosis")), lab = sub(" \\[.*", "", nice(set)),
              group = factor(group, levels = unique(SHOW$group)), sig = !is.na(padj) & padj < 0.05)]
course[, colr := factor(match(set, unique(set)) , levels = 1:8), by = group]
fwrite(course[, .(set, program = nice(set), group = gsub("\n", " ", group), layer, step, value = round(v, 4), padj)],
       file.path(OUT, "programs_course.tsv"), sep = "\t")
course_plot <- function(ly, ylab, ttl) {
  d <- course[layer == ly]
  ggplot(d, aes(x, v, group = set, colour = colr)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = sig), size = 2.4, stroke = 1, fill = "white") +
    ggrepel::geom_text_repel(data = d[x == 3], aes(label = lab), hjust = 0, direction = "y", nudge_x = 0.12, xlim = c(3.1, 5.2),
                             size = 2.9, segment.colour = "grey70", segment.size = 0.3, min.segment.length = 0, box.padding = 0.12, seed = 1) +
    facet_wrap(~ group, ncol = 1, scales = "free_y", strip.position = "left") +
    scale_x_continuous(breaks = 1:3, labels = c("1 month", "3.5 years", "DDHN\ndiagnosis"), limits = c(0.9, 5.2), expand = c(0, 0)) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "FDR < 0.05 against 1 month", `FALSE` = "not significant"), name = NULL) +
    scale_colour_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"), guide = "none") +
    labs(x = NULL, y = ylab, title = ttl) +
    theme_bw(base_size = 11) +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, size = 9.5), legend.position = "bottom",
          panel.grid.minor = element_blank(), plot.title = element_text(size = 12))
}
pcm <- course_plot("methylation", "change around the genes against 1 month, relative to all genes (points; below 0 = lost)",
                   "Methylation around the programs' genes")
pce <- course_plot("expression", "median log2 fold-change against 1 month, relative to all genes (below 0 = down)",
                   "Expression of the programs' genes")
ggsave(file.path(OUT, "fig_programs_course.png"),
       (pcm | pce) + plot_annotation(title = "Programs across the recipient's course: 1 month, 3.5 years and DDHN diagnosis",
                                     subtitle = "each line one program, 1 month = 0; filled points: the enrichment against 1 month has FDR < 0.05") +
         plot_layout(guides = "collect") & theme(legend.position = "bottom"),
       width = 15, height = 17, dpi = 170)

## ---- figure: TET2 vs DNMT3A -----------------------------------------------------------------------------
# A: growth rates as doubling times
gplot <- copy(gr)[, who := factor(clone, levels = rev(clone))]
pa4 <- ggplot(gplot, aes(12 * doubling_years, who, colour = who)) +
  geom_errorbar(aes(xmin = 12 * doubling_lo, xmax = 12 * doubling_hi), width = 0.2, linewidth = 0.8, orientation = "y") +
  geom_point(size = 4) +
  geom_text(aes(label = fifelse(doubling_years < 1, sprintf("%.1f months", 12 * doubling_years), sprintf("%.1f years", doubling_years))),
            vjust = -1.3, size = 3.6, show.legend = FALSE) +
  scale_x_log10(breaks = c(1, 3, 6, 12, 24, 60, 120, 240), labels = c("1 mo", "3 mo", "6 mo", "1 yr", "2 yr", "5 yr", "10 yr", "20 yr")) +
  scale_colour_manual(values = c("#D95F02", "#B39DDB", "#6A51A3"), guide = "none") +
  scale_y_discrete(labels = c("DNMT3A DDHN branch, in the patient" = "DNMT3A DDHN branch, in the patient\nwith the KMT2A duplication",
                              "DNMT3A clone, in the donor" = "DNMT3A clone, in the donor\nwithout the KMT2A duplication",
                              "TET2 DDHN clone, in the patient" = "TET2 DDHN clone, in the patient\nno copy-number change")) +
  labs(x = "time for the clone to double while small (log scale)", y = NULL,
       title = "A. How fast each clone grew", subtitle = "least average rate between the two samples; bar: 95% interval from the starting size") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 10))
# B: the programs of the two DDHNs against early post-transplant marrow, methylation
SEL <- c("PANEL_HOX_MEIS_self_renewal", "PANEL_HSC_stemness", "PANEL_leukaemic_stem_cell", "HAY_BONE_MARROW_CD34_POS_HSC", "PANEL_progenitors",
         "HAY_BONE_MARROW_CD34_POS_GRAN", "HAY_BONE_MARROW_IMMATURE_NEUTROPHIL", "HAY_BONE_MARROW_NEUTROPHIL", "HAY_BONE_MARROW_MONOCYTE",
         "PANEL_granulocyte_monocyte", "PANEL_erythroid", "PANEL_megakaryocyte", "PANEL_lymphoid", "HAY_BONE_MARROW_NK_CELLS")
tb <- prog[layer == "methylation" & id %in% c(paste0(ID[["DX"]], "_vs_", ID[["D30"]]), paste0(ID[["T2DX"]], "_vs_", ID[["T2A"]])) & set %in% SEL]
tb[, who := fifelse(id == paste0(ID[["DX"]], "_vs_", ID[["D30"]]), "DNMT3A: recipient DDHN diagnosis vs 1 month", "TET2: recipient DDHN diagnosis vs 1 year")]
tb[, row := factor(nice(set), levels = rev(nice(SEL)))]
pb4 <- ggplot(tb, aes(NES, row, fill = who)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = fifelse(padj < 0.05, "*", ""), x = NES + sign(NES) * 0.12), position = position_dodge(width = 0.75), vjust = 0.75, size = 5) +
  scale_fill_manual(values = c("DNMT3A: recipient DDHN diagnosis vs 1 month" = "#6A51A3", "TET2: recipient DDHN diagnosis vs 1 year" = "#D95F02"), name = NULL) +
  scale_x_continuous(limits = c(-2.9, 2.9)) +
  labs(x = "NES: methylation lost around the genes (left) or gained (right)", y = NULL,
       title = "B. Methylation programs of the two DDHNs against early post-transplant marrow", subtitle = "* FDR < 0.05") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 10))
# C: absolute direction for the two program families
fam <- data.table(set = SEL, family = fifelse(SEL %in% SEL[1:5], "stem / self-renewal", fifelse(SEL %in% SEL[6:10], "myeloid", "other")))
tc <- merge(prog[layer == "methylation" & id %in% c(paste0(ID[["DX"]], "_vs_", ID[["D30"]]), paste0(ID[["T2DX"]], "_vs_", ID[["T2A"]]))], fam[family != "other"], by = "set")
tc <- tc[, .(lost = 100 * sum(lost) / sum(tiles), gained = 100 * sum(gained) / sum(tiles), bg_lost = bg_lost[1], bg_gained = bg_gained[1]), by = .(id, family)]
tc <- rbind(tc[, .(id, family, lost, gained)], unique(tc[, .(id, family = "all genes", lost = bg_lost, gained = bg_gained)]))
tcl <- melt(tc, id.vars = c("id", "family"), variable.name = "direction", value.name = "pct")
tcl[, who := fifelse(id == paste0(ID[["DX"]], "_vs_", ID[["D30"]]), "DNMT3A: recipient DDHN diagnosis vs 1 month", "TET2: recipient DDHN diagnosis vs 1 year")]
tcl[, family := factor(family, levels = c("stem / self-renewal", "myeloid", "all genes"))]
tcl[, x := fifelse(direction == "lost", -pct, pct)]
pc4 <- ggplot(tcl, aes(x, family, fill = direction)) +
  geom_vline(xintercept = 0, colour = "grey40") + geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f%%", pct), hjust = fifelse(direction == "lost", 1.1, -0.1)), size = 3.3) +
  facet_wrap(~ who, ncol = 1) +
  scale_fill_manual(values = c(lost = "#2166AC", gained = "#B2182B"), name = sprintf("tiles that changed by >= %s points", CUT_TXT)) +
  scale_x_continuous(labels = function(v) paste0(abs(v), "%"), expand = expansion(mult = 0.35)) +
  labs(x = sprintf("share of the %s tiles around the genes (left: lost, right: gained)", TILE_TXT), y = NULL,
       title = "C. The same, as tiles that lost or gained methylation",
       subtitle = "stem: the first five programs of B; myeloid: the next five") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 9.5))
ggsave(file.path(OUT, "fig_tet2_vs_dnmt3a.png"), (pa4 / (pb4 | pc4)) + plot_layout(heights = c(0.55, 1.45)), width = 16, height = 13, dpi = 170)
writeLines(c(paste("msigdb", msig_version), paste("fgsea", as.character(packageVersion("fgsea"))), paste("built", format(Sys.time()))),
           file.path(OUT, "run_info.txt"))
message("written to ", OUT)
