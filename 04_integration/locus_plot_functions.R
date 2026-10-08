###############################################################################
# locus_plot_functions.R                                             2026-10-08
#
# Mesa-paper-style locus plots: gene model and CpG island on genomic
# coordinates, each CpG linked to a column of a per-sample methylation heatmap,
# and the RNA-seq expression of the gene beside it (as panels f and i in the
# Mesa figures).
#
#   load_locus_data()  CpG counts (TAPS reading) for the samples, MANE gene
#                      models, CpG islands and RNA-seq values; call once.
#                      A CpG is shown when it has >= min_cov reads (default 1).
#   locus_plot(L, gene, region = NULL, highlight = NULL, ...)
#                      region    c(chr, start, end); default: the promoter
#                                window around the MANE TSS, shrunk until it
#                                holds <= max_cpg CpGs
#                      highlight c(start, end, label): a box, e.g. a DMR, drawn
#                                where the Mesa figures draw MM
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

if (!exists("SAMPLE_SHEET")) source("config/paths.R")
.sheet <- fread(SAMPLE_SHEET, colClasses = "character")
LOCUS_SAMPLES <- .sheet[methylation == "analysed", .(sample, label = plot_label, case = paste(case, "case"))]
RNA_SAMPLES   <- .sheet[rna_seq == "yes", .(sample, label = rna_label, colour = rna_colour)]
METH_COLS <- c("#16325C", "#3F6E9A", "#F7D488", "#E8873A", "#C0392B")

load_locus_data <- function(emseq, mane_gtf, cgi_file, rna_rdata, min_cov = 1) {
  keep <- paste0("chr", c(1:22, "X"))
  cpg <- lapply(setNames(LOCUS_SAMPLES$sample, LOCUS_SAMPLES$sample), function(s) {
    f <- file.path(emseq, paste0(s, "_msuite2"), "Msuite2.CpG.meth.call.xz")
    x <- fread(cmd = sprintf("xz -dc '%s'", f), select = c(1, 2, 4, 5, 8, 9),
               col.names = c("chr", "pos", "wC", "wT", "cC", "cT"))
    x <- x[chr %in% keep, .(chr, pos, m = wT + cT, t = wC + wT + cC + cT)]   # TAPS: T = methylated
    setkey(x, chr, pos); x
  })
  g <- fread(cmd = sprintf("gzip -dc '%s'", mane_gtf), header = FALSE, sep = "\t", quote = "")
  g <- g[grepl('tag "MANE_Select"', V9, fixed = TRUE) & V3 %in% c("transcript", "exon", "CDS", "UTR")]
  g[, gene := sub('.*gene_name "([^"]+)".*', "\\1", V9)]
  model <- g[, .(gene, chr = V1, feature = V3, start = V4, end = V5, strand = V7)]
  cgi <- fread(cmd = sprintf("gzip -dc '%s'", cgi_file), select = 2:4, col.names = c("chr", "start", "end"))
  cgi[, start := start + 1L]
  e <- new.env(); load(rna_rdata, envir = e)
  list(cpg = cpg, model = model, cgi = cgi, cpm = e$cpm_mat, de = as.data.table(e$de), min_cov = min_cov)
}

.window_for <- function(L, gene, half = 1000, max_cpg = 90) {
  g0 <- gene
  tx <- L$model[gene == g0 & feature == "transcript"][1]
  if (is.na(tx$chr)) stop("no MANE transcript for ", gene)
  tss <- if (tx$strand == "+") tx$start else tx$end
  x <- L$cpg[[1]]
  repeat {                                    # shrink until the window is readable
    a <- tss - half; b <- tss + half
    n <- x[chr == tx$chr & pos >= a & pos <= b, .N]
    if (n <= max_cpg || half <= 200) break
    half <- round(half * 0.85)
  }
  list(chr = tx$chr, start = a, end = b, tss = tss, strand = tx$strand)
}

locus_plot <- function(L, gene, region = NULL, highlight = NULL, max_cpg = 90,
                       samples = LOCUS_SAMPLES$sample, rna_note = NULL, title = NULL, all_genes = FALSE,
                       subtitle = "Methylation per CpG (TAPS); * provisional sample role") {
  g0 <- gene
  w <- if (is.null(region)) .window_for(L, gene, max_cpg = max_cpg) else {
    tx <- L$model[gene == g0 & feature == "transcript"][1]
    list(chr = region[[1]], start = as.integer(region[[2]]), end = as.integer(region[[3]]),
         tss = if (!is.na(tx$chr)) ifelse(tx$strand == "+", tx$start, tx$end) else NA, strand = tx$strand)
  }
  span <- w$end - w$start
  nx <- function(p) (p - w$start) / span                      # genomic -> [0, 1]

  ## CpGs: union of positions with >= min_cov reads in any sample
  d <- rbindlist(lapply(samples, function(s)
    L$cpg[[s]][chr == w$chr & pos >= w$start & pos <= w$end][, sample := s]))
  pos_keep <- sort(unique(d[t >= L$min_cov, pos]))
  if (!length(pos_keep)) stop("no covered CpGs in the window for ", gene)
  d <- d[pos %in% pos_keep]
  d[, meth := fifelse(t >= L$min_cov, m / t, NA_real_)]
  cols <- data.table(pos = pos_keep, col = seq_along(pos_keep))
  ncol <- length(pos_keep)
  grid <- CJ(sample = samples, pos = pos_keep)
  grid <- merge(merge(grid, d[, .(sample, pos, meth)], by = c("sample", "pos"), all.x = TRUE), cols, by = "pos")
  lab <- LOCUS_SAMPLES[match(samples, sample)]
  grid[, row := factor(lab$label[match(sample, lab$sample)], levels = rev(lab$label))]
  grid[, xc := (col - 0.5) / ncol]

  ## gene model track
  m <- if (all_genes) L$model[chr == w$chr & end >= w$start & start <= w$end]
       else L$model[gene == g0 & chr == w$chr & end >= w$start & start <= w$end]
  glab <- m[feature == "transcript", .(x = nx((pmax(start, w$start) + pmin(end, w$end)) / 2), gene, strand)]
  ex <- m[feature %in% c("CDS", "UTR")]
  ex[, `:=`(h = fifelse(feature == "CDS", 0.18, 0.10), s = pmax(start, w$start), e = pmin(end, w$end))]
  txr <- m[feature == "transcript"]
  cg <- L$cgi[chr == w$chr & end >= w$start & start <= w$end][, `:=`(s = pmax(start, w$start), e = pmin(end, w$end))]
  brk <- pretty(c(w$start, w$end), n = 5); brk <- brk[brk >= w$start & brk <= w$end]
  gt <- ggplot() +
    { if (nrow(txr)) geom_segment(data = txr, aes(x = nx(pmax(start, w$start)), xend = nx(pmin(end, w$end)), y = 0.5, yend = 0.5), linewidth = 0.6) } +
    { if (all_genes && nrow(glab)) geom_text(data = glab, aes(x = x, y = 0.24, label = paste0(gene, ifelse(strand == "+", " >", " <"))),
                                             size = 2.8, fontface = "italic") } +
    { if (nrow(ex)) geom_rect(data = ex, aes(xmin = nx(s), xmax = nx(e), ymin = 0.5 - h, ymax = 0.5 + h,
                                              fill = feature), colour = "black", linewidth = 0.2) } +
    { if (nrow(cg)) geom_rect(data = cg, aes(xmin = nx(s), xmax = nx(e), ymin = 0.80, ymax = 0.88), fill = "#2E8B87") } +
    { if (nrow(cg)) annotate("text", x = nx(cg$s[1]), y = 0.995, label = "CpG island", hjust = 0, vjust = 1,
                              colour = "#1F6F6B", size = 2.8, fontface = "bold") } +
    { if (!is.null(highlight)) annotate("rect", xmin = nx(as.numeric(highlight[1])), xmax = nx(as.numeric(highlight[2])),
                                        ymin = 0, ymax = 1, fill = "#7B9BD1", alpha = 0.18, colour = "#5A7BC0") } +
    { if (!is.null(highlight) && length(highlight) > 2) annotate("text", x = mean(nx(as.numeric(highlight[1:2]))), y = 0.06,
                                                                  label = highlight[3], size = 3, fontface = "bold", colour = "#3A5BA0") } +
    { if (!is.na(w$tss) && w$tss >= w$start && w$tss <= w$end)
        list(annotate("segment", x = nx(w$tss), xend = nx(w$tss), y = 0.5, yend = 0.70),
             annotate("segment", x = nx(w$tss), xend = nx(w$tss) + ifelse(w$strand == "+", 0.04, -0.04),
                      y = 0.70, yend = 0.70, arrow = arrow(length = unit(0.08, "in"))),
             annotate("text", x = nx(w$tss) + ifelse(w$strand == "+", -0.005, 0.005), y = 0.70,
                      label = "TSS", size = 3, hjust = ifelse(w$strand == "+", 1, 0))) } +
    geom_segment(data = data.table(p = pos_keep), aes(x = nx(p), xend = nx(p), y = 0, yend = 0.08), linewidth = 0.25) +
    scale_fill_manual(values = c(CDS = "#B23A2F", UTR = "#D98C5F"), name = NULL) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = nx(brk),
                       labels = format(brk, big.mark = ","), position = "top") +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = sprintf("%s (hg38)%s", w$chr, ifelse(is.na(w$strand), "", paste0("   gene on ", w$strand, " strand"))), y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(), legend.position = "none")

  ## connector lines: CpG position -> heatmap column
  cn <- ggplot(cols[, .(x0 = nx(pos), x1 = (col - 0.5) / ncol)]) +
    geom_segment(aes(x = x0, xend = x1, y = 1, yend = 0), linewidth = 0.25, colour = "grey30") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void()

  ## heatmap
  num_brk <- unique(c(1, seq(5, ncol, by = 5)))
  hm <- ggplot(grid, aes(xc, row, fill = meth)) +
    geom_tile(width = 1 / ncol, colour = "white", linewidth = 0.15) +
    scale_fill_gradientn(colours = METH_COLS, limits = c(0, 1), na.value = "grey88", name = "Methylation") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = (num_brk - 0.5) / ncol, labels = num_brk) +
    labs(x = sprintf("CpG #  (%d CpGs; grey = %s)", ncol,
                     ifelse(L$min_cov <= 1, "no reads", sprintf("< %d reads", L$min_cov))), y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), legend.position = "bottom", legend.key.width = unit(0.5, "in"),
          legend.key.height = unit(0.12, "in"))

  ## RNA panel
  rp <- if (gene %in% rownames(L$cpm)) {
    r <- data.table(sample = RNA_SAMPLES$sample, cpm = as.numeric(L$cpm[gene, RNA_SAMPLES$sample]))
    r[, label := factor(RNA_SAMPLES$label, levels = RNA_SAMPLES$label)]
    st <- L$de[gene == g0 & contrast %in% c("C1_clone", "C2_transform", "C4_DDHN_vs_donor")]
    note <- if (is.null(rna_note)) paste(sprintf("%s: %+.2f (FDR %.2g)",
               c(C1_clone = "3.5yr vs D30", C2_transform = "DDHN vs 3.5yr", C4_DDHN_vs_donor = "DDHN vs donor")[st$contrast],
               st$lfc_shrunk, st$FDR), collapse = "\n") else rna_note
    ggplot(r, aes(label, cpm, fill = label)) + geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
      scale_fill_manual(values = setNames(RNA_SAMPLES$colour, RNA_SAMPLES$label), guide = "none") +
      labs(x = NULL, y = "RNA-seq (CPM)", title = "Expression", caption = paste0("log2FC (shrunk)\n", note)) +
      theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                                      plot.caption = element_text(hjust = 0, size = 7.5))
  } else ggplot() + annotate("text", x = 0, y = 0, label = "not in the\nRNA-seq data", size = 3.2) + theme_void()

  left <- gt / cn / hm + plot_layout(heights = c(1.2, 0.8, 0.22 * length(samples) + 0.35))
  (wrap_elements(left) | rp) + plot_layout(widths = c(4, 1)) +
    plot_annotation(title = if (is.null(title)) gene else title,
                    subtitle = subtitle,
                    theme = theme(plot.title = element_text(face = "bold.italic", size = 14)))
}


# Arc diagram: the DMR and the partner anchors it loops to (09_hic_loops.R), each partner
# labelled with its genes and the most significant change metilene tested inside it.
loop_panel <- function(dmr_chr, dmr_start, dmr_end, partners) {
  if (is.null(partners) || !nrow(partners)) return(NULL)
  p <- copy(partners)
  p[, ps := as.integer(sub(".*:(\\d+)-.*", "\\1", partner))][, pe := as.integer(sub(".*-(\\d+)$", "\\1", partner))]
  lo <- min(dmr_start, p$ps); hi <- max(dmr_end, p$pe); span <- hi - lo
  xm <- (dmr_start + dmr_end) / 2
  arcs <- rbindlist(lapply(seq_len(nrow(p)), function(i) {
    x2 <- (p$ps[i] + p$pe[i]) / 2; t <- seq(0, pi, length.out = 60)
    data.table(i = i, x = (xm + x2) / 2 + (x2 - xm) / 2 * cos(t), y = sin(t) * pmin(1, abs(x2 - xm) / span + 0.25))
  }))
  p[, lab := sprintf("%s%s\nbest tested change %+.2f (q %.2g)%s; %d/8 files: %s",
                     ifelse(is.na(partner_genes), "no TSS", partner_genes),
                     ifelse(is.na(partner_gene_rna), "", paste0(" [RNA ", partner_gene_rna, "]")),
                     partner_best_delta, partner_best_q, ifelse(partner_n_dmr > 0, " DMR", ""), support, cell_type)]
  ggplot() +
    geom_path(data = arcs, aes(x, y, group = i), colour = "#6A51A3", linewidth = 0.6) +
    annotate("rect", xmin = dmr_start, xmax = dmr_end + span * 0.004, ymin = -0.08, ymax = 0, fill = "#B2182B") +
    annotate("text", x = xm, y = -0.16, label = "DMR", size = 3, colour = "#B2182B", fontface = "bold") +
    geom_rect(data = p, aes(xmin = ps, xmax = pe, ymin = -0.08, ymax = 0,
                            fill = pmax(-0.5, pmin(0.5, partner_best_delta))), colour = "grey30", linewidth = 0.2) +
    geom_text(data = p, aes(x = (ps + pe) / 2, y = -0.2, label = lab), size = 2.4, vjust = 1, lineheight = 0.9) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-0.5, 0.5),
                         name = "partner\nchange") +
    scale_x_continuous(labels = function(v) format(v, big.mark = ","), expand = expansion(mult = 0.08)) +
    coord_cartesian(ylim = c(-0.55, 1.05), clip = "off") +
    labs(x = sprintf("%s: Hi-C loops (HSC/CMP/GMP/MEP, Pimanda; >= 2 of 8 files)", dmr_chr), y = NULL) +
    theme_minimal(base_size = 9) + theme(axis.text.y = element_blank(), panel.grid = element_blank(),
                                         legend.key.height = unit(0.15, "in"))
}
