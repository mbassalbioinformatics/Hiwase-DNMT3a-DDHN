###############################################################################
# locus_plot_functions.R                                             2026-10-09
#
# Mesa-paper-style locus plots: gene models and CpG islands on genomic
# coordinates, each CpG linked to a column of a per-sample methylation heatmap,
# and the RNA-seq expression of the gene beside it (as panels f and i in the
# Mesa figures).
#
#   load_locus_data()  CpG counts (TAPS reading) for the samples, MANE gene
#                      models, CpG islands and RNA-seq values; call once.
#                      A CpG is shown when it has >= min_cov reads (default 2).
#   locus_plot(L, gene, region = NULL, highlight = NULL, ...)
#                      region    c(chr, start, end); default: the promoter
#                                window around the MANE TSS, shrunk until it
#                                holds <= max_cpg CpGs
#                      highlight c(start, end, label): a box, e.g. a DMR, drawn
#                                where the Mesa figures draw MM
#                      The plot carries attr(, "height_in"), the height it needs
#                      at width 14 in, so text keeps one size in every plot.
#   loop_panel(...)    Hi-C loops from a DMR to its partner regions, with
#                      attr(, "height_in")
#
# The gene track is laid out in units of 0.1 inch so that no two texts meet:
# CpG ticks at the bottom; genes in lanes, one lane per overlapping gene (at
# most four; any further genes are named in the subtitle), each labelled above
# its own model and kept clear of the other genes and of the TSS arrow; the CpG
# island bar and its label above the lanes; the DMR label in a band of its own
# at the top.
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(ggrepel) })

if (!exists("SAMPLE_SHEET")) source("config/paths.R")
.sheet <- fread(SAMPLE_SHEET, colClasses = "character")
LOCUS_SAMPLES <- .sheet[methylation == "analysed", .(sample, label = plot_label, case = paste(case, "case"))]
RNA_SAMPLES   <- .sheet[rna_seq == "yes", .(sample, label = rna_label, colour = rna_colour)]
METH_COLS <- c("#16325C", "#3F6E9A", "#F7D488", "#E8873A", "#C0392B")
TRACK_UNIT <- 0.1        # inches per layout unit of the gene track
MAX_LANES  <- 4L
LABEL_SIZE <- 2.9        # mm; gene, CpG island and DMR labels

load_locus_data <- function(emseq, mane_gtf, cgi_file, rna_rdata, min_cov = 2) {
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

# Lanes for the gene track. The primary gene comes first, the rest by position.
# A gene occupies its visible extent plus its label (in panel fractions); it
# goes to the lowest lane where it meets nothing. The primary gene's label is
# moved off its TSS arrow.
.gene_lanes <- function(txr, primary, tss_x, char_w) {
  txr <- copy(txr)[order(gene != primary, s)]
  txr[, lw := nchar(lab) * char_w]
  txr[, lx := pmin(pmax((s + e) / 2, lw / 2 + 0.005), 1 - lw / 2 - 0.005)]
  i <- which(txr$gene == primary)
  if (!is.na(tss_x) && length(i) && abs(txr$lx[i] - tss_x) < txr$lw[i] / 2 + 0.05) {
    right <- tss_x + 0.05 + txr$lw[i] / 2; left <- tss_x - 0.05 - txr$lw[i] / 2
    if (right + txr$lw[i] / 2 <= 0.995) txr$lx[i] <- right else if (left - txr$lw[i] / 2 >= 0.005) txr$lx[i] <- left
  }
  txr[, `:=`(a = pmin(s, lx - lw / 2) - 0.012, b = pmax(e, lx + lw / 2) + 0.012, lane = NA_integer_)]
  occ <- vector("list", MAX_LANES)
  for (j in seq_len(nrow(txr))) for (k in seq_len(MAX_LANES)) {
    o <- occ[[k]]
    if (is.null(o) || all(txr$b[j] < o[, 1] | txr$a[j] > o[, 2])) {
      occ[[k]] <- rbind(o, c(txr$a[j], txr$b[j])); txr$lane[j] <- k; break
    }
  }
  txr
}

locus_plot <- function(L, gene, region = NULL, highlight = NULL, max_cpg = 90,
                       samples = LOCUS_SAMPLES$sample, rna_note = NULL, title = NULL, all_genes = FALSE,
                       subtitle = "Methylation per CpG (TAPS)", width_in = 14) {
  g0 <- gene
  w <- if (is.null(region)) .window_for(L, gene, max_cpg = max_cpg) else {
    tx <- L$model[gene == g0 & feature == "transcript"][1]
    list(chr = region[[1]], start = as.integer(region[[2]]), end = as.integer(region[[3]]),
         tss = if (!is.na(tx$chr) && tx$chr == region[[1]]) ifelse(tx$strand == "+", tx$start, tx$end) else NA,
         strand = tx$strand)
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

  ## gene track (y in layout units of TRACK_UNIT inch)
  m <- if (all_genes) L$model[chr == w$chr & end >= w$start & start <= w$end]
       else L$model[gene == g0 & chr == w$chr & end >= w$start & start <= w$end]
  txr <- unique(m[feature == "transcript", .(gene, strand, s = nx(pmax(start, w$start)), e = nx(pmin(end, w$end)))], by = "gene")
  txr[, lab := paste0(gene, ifelse(strand == "+", "  >", "  <"))]
  tss_x <- if (!is.na(w$tss) && w$tss >= w$start && w$tss <= w$end) nx(w$tss) else NA_real_
  char_w <- 0.0062 * 14 / width_in                            # label width per character, panel fraction
  hidden <- character(0)
  if (!nrow(txr)) txr <- data.table(gene = character(), strand = character(), s = numeric(), e = numeric(),
                                    lab = character(), lx = numeric(), lane = integer())
  if (nrow(txr)) {
    txr <- .gene_lanes(txr, g0, tss_x, char_w)
    hidden <- txr[is.na(lane), gene]; txr <- txr[!is.na(lane)]
  }
  nl <- if (nrow(txr)) max(txr$lane) else 1L
  txr[, y := 2.0 + (lane - 1) * 3.4]                          # lane centre
  ex <- merge(m[feature %in% c("CDS", "UTR"), .(gene, feature, s = nx(pmax(start, w$start)), e = nx(pmin(end, w$end)))],
              txr[, .(gene, y)], by = "gene")
  ex[, h := fifelse(feature == "CDS", 0.55, 0.32)]
  top_l <- 1.2 + nl * 3.4 - 0.4                                # top of the gene lanes
  H <- top_l + if (is.null(highlight)) 2.6 else 4.1            # track height, layout units
  cg <- L$cgi[chr == w$chr & end >= w$start & start <= w$end][, `:=`(s = nx(pmax(start, w$start)), e = nx(pmin(end, w$end)))]
  cg_lab_x <- if (nrow(cg)) cg$s[1] else NA_real_
  cg_hjust <- if (nrow(cg) && cg$s[1] + 10 * char_w > 0.995) 1 else 0
  if (cg_hjust == 1) cg_lab_x <- cg$e[1]
  py <- txr[gene == g0, y]
  brk <- pretty(c(w$start, w$end), n = 5); brk <- brk[brk >= w$start & brk <= w$end]
  hl <- if (!is.null(highlight)) as.numeric(highlight[1:2]) else NULL
  gt <- ggplot() +
    { if (!is.null(hl)) annotate("rect", xmin = nx(hl[1]), xmax = nx(hl[2]), ymin = 0, ymax = top_l + 0.9,
                                 fill = "#7B9BD1", alpha = 0.18, colour = "#5A7BC0") } +
    { if (!is.null(hl) && length(highlight) > 2)
        annotate("label", x = min(max(mean(nx(hl)), 0.03), 0.97), y = top_l + 3.15, label = highlight[3],
                 size = LABEL_SIZE, fontface = "bold", colour = "#3A5BA0", fill = "white") } +
    { if (nrow(txr)) geom_segment(data = txr, aes(x = s, xend = e, y = y, yend = y), linewidth = 0.6) } +
    { if (nrow(ex)) geom_rect(data = ex, aes(xmin = s, xmax = e, ymin = y - h, ymax = y + h, fill = feature),
                              colour = "black", linewidth = 0.2) } +
    { if (nrow(txr)) geom_text(data = txr, aes(x = lx, y = y + 0.72, label = lab), vjust = 0, size = LABEL_SIZE,
                               fontface = "italic") } +
    { if (!is.na(tss_x) && length(py))
        list(annotate("segment", x = tss_x, xend = tss_x, y = py, yend = py + 1.05, linewidth = 0.5),
             annotate("segment", x = tss_x, xend = tss_x + ifelse(w$strand == "+", 0.03, -0.03), y = py + 1.05,
                      yend = py + 1.05, linewidth = 0.5, arrow = arrow(length = unit(0.07, "in")))) } +
    { if (nrow(cg)) geom_rect(data = cg, aes(xmin = s, xmax = e, ymin = top_l + 0.3, ymax = top_l + 0.8), fill = "#2E8B87") } +
    { if (nrow(cg)) annotate("text", x = cg_lab_x, y = top_l + 1.25, label = "CpG island", hjust = cg_hjust, vjust = 0,
                             colour = "#1F6F6B", size = LABEL_SIZE, fontface = "bold") } +
    geom_segment(data = data.table(p = pos_keep), aes(x = nx(p), xend = nx(p), y = 0, yend = 0.6), linewidth = 0.25) +
    scale_fill_manual(values = c(CDS = "#B23A2F", UTR = "#D98C5F"), limits = c("CDS", "UTR"), name = NULL) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = nx(brk),
                       labels = format(brk, big.mark = ","), position = "top") +
    scale_y_continuous(limits = c(0, H), expand = c(0, 0)) +
    labs(x = sprintf("%s (hg38)", w$chr), y = NULL) +
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
      labs(x = NULL, y = "RNA-seq (CPM)", title = paste0(gene, " expression"), caption = paste0("log2FC (shrunk)\n", note)) +
      theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 35, hjust = 1),
                                      plot.title = element_text(face = "italic"),
                                      plot.caption = element_text(hjust = 0, size = 7.5))
  } else ggplot() + annotate("text", x = 0, y = 0, label = paste0(gene, ":\nnot in the\nRNA-seq data"), size = 3.2) + theme_void()

  track_in <- H * TRACK_UNIT
  heat_in  <- 0.24 * length(samples) + 0.2
  if (length(hidden)) subtitle <- paste0(subtitle, "\nAlso in this window, not drawn: ", paste(hidden, collapse = ", "))
  n_title <- 1 + lengths(regmatches(title %||% "", gregexpr("\n", title %||% "")))
  n_sub   <- 1 + lengths(regmatches(subtitle, gregexpr("\n", subtitle)))
  left <- gt / cn / hm + plot_layout(heights = c(track_in, 0.75, heat_in))
  # both columns wrapped whole: otherwise patchwork aligns the left column to the
  # expression panel's plotting area and squeezes the gene track
  p <- (wrap_elements(full = left) | wrap_elements(full = rp)) + plot_layout(widths = c(4, 1)) +
    plot_annotation(title = if (is.null(title)) gene else title, subtitle = subtitle,
                    theme = theme(plot.title = element_text(face = "bold.italic", size = 14)))
  attr(p, "height_in") <- track_in + 0.75 + heat_in + 1.45 + 0.28 * n_title + 0.2 * n_sub
  p
}
`%||%` <- function(a, b) if (is.null(a)) b else a


# Arc diagram: the DMR and the partner anchors it loops to (09_hic_loops.R), each partner
# labelled with its genes and the most significant change metilene tested inside it.
# Labels are repelled sideways so that neighbouring partners never overlap.
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
  short_genes <- function(g) {
    if (is.na(g) || !nzchar(g)) return("no gene TSS")
    v <- strsplit(g, ",")[[1]]
    if (length(v) > 4) paste0(paste(v[1:4], collapse = ", "), " +", length(v) - 4) else paste(v, collapse = ", ")
  }
  p[, tested := fifelse(is.na(partner_best_delta), "no region tested",
                        sprintf("best tested change %+.2f (q %s)%s", partner_best_delta, sprintf("%.2g", partner_best_q),
                                fifelse(!is.na(partner_n_dmr) & partner_n_dmr > 0, ", DMR", "")))]
  p[, lab := sprintf("%s%s\n%s\n%d of 8 files: %s",
                     vapply(partner_genes, short_genes, ""),
                     ifelse(is.na(partner_gene_rna), "", paste0("  [RNA ", partner_gene_rna, "]")),
                     tested, support, gsub(",", ", ", cell_type))]
  labs_dt <- rbind(p[, .(x = (ps + pe) / 2, lab, kind = "partner")], data.table(x = xm, lab = "DMR", kind = "dmr"))
  xl <- c(lo - 0.12 * span, hi + 0.12 * span)
  g <- ggplot() +
    geom_path(data = arcs, aes(x, y, group = i), colour = "#6A51A3", linewidth = 0.6) +
    annotate("rect", xmin = dmr_start, xmax = dmr_end + span * 0.004, ymin = -0.08, ymax = 0, fill = "#B2182B") +
    geom_rect(data = p, aes(xmin = ps, xmax = pe, ymin = -0.08, ymax = 0,
                            fill = pmax(-0.5, pmin(0.5, partner_best_delta))), colour = "grey30", linewidth = 0.2) +
    geom_text_repel(data = labs_dt, aes(x = x, y = -0.12, label = lab, colour = kind), size = 2.6, lineheight = 0.95,
                    direction = "x", nudge_y = -0.18, vjust = 1, hjust = 0.5, segment.size = 0.2, segment.colour = "grey60",
                    min.segment.length = 0.2, box.padding = 0.4, xlim = xl, seed = 1, max.overlaps = Inf) +
    scale_colour_manual(values = c(partner = "grey15", dmr = "#B2182B"), guide = "none") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-0.5, 0.5),
                         name = "partner\nchange") +
    scale_x_continuous(labels = function(v) format(v, big.mark = ","), limits = xl, expand = c(0, 0)) +
    coord_cartesian(ylim = c(-0.9, 1.05), clip = "off") +
    labs(x = sprintf("%s: Hi-C loops (HSC, CMP, GMP, MEP; GSE231486; loops in >= 2 of 8 files)", dmr_chr), y = NULL) +
    theme_minimal(base_size = 9) + theme(axis.text.y = element_blank(), panel.grid = element_blank(),
                                         legend.key.height = unit(0.15, "in"))
  attr(g, "height_in") <- 2.6
  g
}
