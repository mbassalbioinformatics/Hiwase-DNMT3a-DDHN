################################################################################
# 01_load_counts.R  --  featureCounts readers (functions only, no side effects)
#
# FEATURE ROWS
#   featureCounts was run with `-f`, so a row is a FEATURE (one GTF exon or one
#   5'UTR line), not a gene. A feature shared between N transcripts appears N
#   times with identical coordinates and identical counts. Summing the raw rows
#   inflates every gene by its isoform count.
#     GENES: 1,624,384 rows -> 674,142 unique exon coordinates -> 61,806 genes
#     5'UTR:   171,472 rows -> 111,594 unique 5'UTR coordinates -> 19,633 genes
#   (both loaders print these on every run -- quote the printed values)
#
#   The fix is de-duplicate on (Geneid, Chr, Start, End, Strand) THEN sum.
#
# WHAT THE GENE VALUE STILL IS NOT: a fragment count.
#   De-duplication removes exact duplicates only. With -f and -O a fragment
#   overlapping two exon features is counted to both -- every spliced fragment,
#   and every fragment where alternative exon boundaries overlap (47.5% of the
#   unique exons here overlap another exon of the same gene) -- and -M counts a
#   multi-mapper at each location. Summed counts divided by sequenced fragments
#   give a MINIMUM of 1.08-2.25 counts per fragment across the four libraries.
#   Consequences:
#     * Between-sample ratios cancel this only to first order. A gene's
#       multiplier depends on its exon structure AND on each library's fragment
#       length and coverage profile, so it is not guaranteed identical across
#       libraries. TMM absorbs the global part; the gene-specific remainder was
#       not measured.
#     * Comparisons BETWEEN genes within one library are invalid: a many-exon
#       gene is inflated more than a single-exon gene. No TPM-like use.
#     * The count model sees more counts than there were fragments, which
#       understates Poisson noise -- mainly relevant at low counts.
#   Counting genes without -f (featureCounts meta-feature mode) avoids all three.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
})

# ---- shared reader ---------------------------------------------------------
# Returns the raw featureCounts table with BAM paths renamed to sample IDs and
# column identity verified by name.
.read_fc <- function(path) {
  if (!file.exists(path)) stop("counts file not found: ", path)
  dt <- fread(path, skip = 1, sep = "\t", header = TRUE,
              showProgress = FALSE, data.table = TRUE)

  meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length",
                 "gene_name", "transcript_id")
  if (!all(meta_cols %in% names(dt)))
    stop("unexpected featureCounts header in ", basename(path), ": ",
         paste(setdiff(meta_cols, names(dt)), collapse = ", "), " missing")

  bam_cols <- setdiff(names(dt), meta_cols)
  # ".../<sample>_qual_trimmed.sorted.bam" -> "<sample>"
  sample_ids <- sub("_qual_trimmed\\.sorted\\.bam$", "", basename(bam_cols))
  setnames(dt, bam_cols, sample_ids)

  if (!setequal(sample_ids, SAMPLE_ORDER))
    stop("BAM columns in ", basename(path), " are (",
         paste(sample_ids, collapse = ", "),
         ") which does not match SAMPLE_ORDER (",
         paste(SAMPLE_ORDER, collapse = ", "), ")")

  # reorder to canonical order by NAME, never by position
  setcolorder(dt, c(meta_cols, SAMPLE_ORDER))
  dt[]
}

# ---- gene level ------------------------------------------------------------
# De-duplicate exons on coordinates within a gene, then sum to gene level.
# Returns list(counts = matrix [gene x sample], anno = data.table, qc = list).
load_gene_counts <- function(path = PATHS$counts_gene, verbose = TRUE) {
  dt <- .read_fc(path)
  n_raw <- nrow(dt)

  # one row per unique exon within a gene
  dt[, .dupkey := paste(Geneid, Chr, Start, End, Strand, sep = "|")]
  dt_u <- dt[!duplicated(.dupkey)]
  n_uniq <- nrow(dt_u)

  # sum unique exons -> gene
  gene_counts <- dt_u[, lapply(.SD, sum), by = .(Geneid), .SDcols = SAMPLE_ORDER]

  # gene-level annotation: gene_name (first, they are constant within Geneid),
  # exonic length = sum of unique exon lengths (a usable proxy for gene length)
  gene_anno <- dt_u[, .(gene_name = gene_name[1],
                        chr       = Chr[1],
                        strand    = Strand[1],
                        n_exons   = .N,
                        exonic_len = sum(Length)), by = .(Geneid)]

  stopifnot(identical(gene_counts$Geneid, gene_anno$Geneid))

  mat <- as.matrix(gene_counts[, ..SAMPLE_ORDER])
  rownames(mat) <- gene_counts$Geneid
  check_sample_order(mat)

  qc <- list(rows_raw = n_raw, rows_unique_exons = n_uniq,
             genes = nrow(mat),
             lib_size_raw = colSums(mat))

  if (verbose)
    message(sprintf("gene level: %s rows -> %s unique exons -> %s genes",
                    format(n_raw, big.mark = ","),
                    format(n_uniq, big.mark = ","),
                    format(nrow(mat), big.mark = ",")))

  list(counts = mat, anno = gene_anno, qc = qc)
}

# ---- 5'UTR / promoter level ------------------------------------------------
# Each DISTINCT 5'UTR coordinate within a gene is treated as one candidate
# promoter / TSS. Features are NOT collapsed to gene level -- collapsing is
# exactly what destroys the promoter-usage signal this table exists to measure.
#
# Why this table matters here more than in a normal RNA-seq project: promoter
# and first-exon CpG methylation is the substrate DNMT3A and TET2 act on, so
# 5'UTR usage is the direct transcriptional readout of the EM-seq hypothesis.
load_utr_counts <- function(path = PATHS$counts_5utr, verbose = TRUE) {
  dt <- .read_fc(path)
  n_raw <- nrow(dt)

  dt[, .dupkey := paste(Geneid, Chr, Start, End, Strand, sep = "|")]
  dt_u <- dt[!duplicated(.dupkey)]

  # Stable, human-readable AND unique feature id. The Geneid must be in the key:
  # distinct ENSG entries can share a gene_name and a 5'UTR coordinate (e.g.
  # ENSG00000284770 / ENSG00000285053 both named TBCE), and a symbol-only id
  # collides for those. A duplicated rowname makes matrix-by-name indexing
  # silently return the first match, which corrupts effect sizes rather than
  # raising an error -- so this is checked, not assumed.
  dt_u[, feature_id := sprintf("%s|%s:%s:%d-%d:%s",
                               fifelse(is.na(gene_name) | gene_name == "",
                                       Geneid, gene_name),
                               Geneid, Chr, Start, End, Strand)]
  if (anyDuplicated(dt_u$feature_id))
    stop("5'UTR feature_id is not unique after including Geneid -- ",
         sum(duplicated(dt_u$feature_id)), " collisions; fix the key before continuing")
  # transcripts sharing an identical 5'UTR are collapsed into one feature;
  # record how many, as a promoter supported by many isoforms is more credible
  n_tx <- dt[, .(n_tx = .N), by = .dupkey]
  dt_u <- merge(dt_u, n_tx, by = ".dupkey", sort = FALSE)

  mat <- as.matrix(dt_u[, ..SAMPLE_ORDER])
  rownames(mat) <- dt_u$feature_id
  check_sample_order(mat)

  anno <- dt_u[, .(feature_id, Geneid, gene_name, chr = Chr, start = Start,
                   end = End, strand = Strand, width = Length, n_tx)]

  # number of distinct promoters per gene -- only genes with >= 2 can switch
  anno[, n_utr_in_gene := .N, by = Geneid]

  qc <- list(rows_raw = n_raw, features = nrow(mat),
             genes = uniqueN(anno$Geneid),
             genes_multi_promoter = uniqueN(anno[n_utr_in_gene >= 2, Geneid]),
             lib_size_raw = colSums(mat))

  if (verbose)
    message(sprintf("5'UTR level: %s rows -> %s unique 5'UTR features across %s genes (%s with >=2 promoters)",
                    format(n_raw, big.mark = ","),
                    format(nrow(mat), big.mark = ","),
                    format(qc$genes, big.mark = ","),
                    format(qc$genes_multi_promoter, big.mark = ",")))

  list(counts = mat, anno = anno, qc = qc)
}

# ---- 5'UTR promoter clusters -----------------------------------------------
# Nested and overlapping 5'UTRs of one gene are NOT independent measurements.
# featureCounts was run with -O, so a fragment in the shared stretch is counted
# in every feature covering it. Before this step was added, on this dataset,
# 41% of within-gene feature pairs overlapped, 33% of tested features had a
# count vector identical to another feature of the same gene, and 42% of
# promoter-switch calls involved an overlapping feature. Treating those as
# separate promoters splits one promoter's reads between copies and distorts
# every within-gene share.
#
# Procedure (worked out on the HuH-6/SALL4 project):
#   1. drop features shorter than THRESH$utr_min_width. Under -O a feature of a
#      few bases collects reads from the neighbouring exon and marks no start site
#   2. within each gene and strand, merge overlapping features into clusters
#   3. represent each cluster by its LONGEST member and take that member's
#      counts. Members are never summed: they describe the same reads
#   4. the cluster's TSS is the representative's Start (+ strand) or End (-)
# The representative keeps its original feature_id, so IDs, the BED export and
# coordinate parsing downstream are unchanged.
#
# Merging gives up separating start sites whose 5'UTRs overlap. With short reads
# counted under -O that separation would rest only on the non-shared stretch.
cluster_utr_promoters <- function(u, min_width = THRESH$utr_min_width, verbose = TRUE) {
  a <- copy(u$anno)
  stopifnot(identical(a$feature_id, rownames(u$counts)))
  a[, dropped_short := width < min_width]
  k <- a[dropped_short == FALSE]

  # coordinates must be double: shift(fill = -Inf) on an integer column silently
  # fills NA, and cummax then collapses every gene into a single cluster.
  # data.table::shift is named explicitly: org.Hs.eg.db attaches IRanges, whose
  # shift() masks data.table's and has no `fill` argument.
  k[, `:=`(s = as.numeric(start), e = as.numeric(end))]
  setorder(k, Geneid, strand, s, e)
  k[, prev_end := data.table::shift(cummax(e), fill = -Inf), by = .(Geneid, strand)]
  k[, cluster := cumsum(s > prev_end), by = .(Geneid, strand)]
  stopifnot(!anyNA(k$cluster))
  k[, cluster_key := paste(Geneid, strand, cluster, sep = "#")]

  span <- k[, .(cluster_start = min(s), cluster_end = max(e), n_members = .N,
                n_tx_cluster = sum(n_tx)), by = cluster_key]
  # one representative per cluster: the longest member (first if tied)
  rep_rows <- k[, .I[which.max(width)], by = cluster_key]$V1
  rep <- merge(k[rep_rows], span, by = "cluster_key", sort = FALSE)
  rep[, tss := fifelse(strand == "+", s, e)]

  # ---- checks: fail loudly rather than return a subtly wrong table ----------
  stopifnot(!anyDuplicated(rep$feature_id), !anyDuplicated(rep$cluster_key))
  longest <- k[, .(max_w = max(width)), by = cluster_key]
  stopifnot(all(rep$width == longest$max_w[match(rep$cluster_key, longest$cluster_key)]))
  ordered <- rep[order(Geneid, strand, cluster_start)]
  ordered[, next_start := data.table::shift(cluster_start, type = "lead"), by = .(Geneid, strand)]
  stopifnot(ordered[!is.na(next_start), all(next_start > cluster_end)])   # clusters disjoint
  stopifnot(all(rep$width >= min_width))

  # keep the original file order of the representatives
  rep <- rep[order(match(feature_id, a$feature_id))]
  mat <- u$counts[rep$feature_id, , drop = FALSE]
  stopifnot(identical(rownames(mat), rep$feature_id))
  check_sample_order(mat)

  anno <- rep[, .(feature_id, Geneid, gene_name, chr, start, end, strand, width,
                  tss, cluster_start, cluster_end, n_members, n_tx = n_tx_cluster)]
  anno[, n_promoters_in_gene := .N, by = Geneid]

  # every input feature, and which cluster it went to (NA if dropped as short)
  k[, representative := rep$feature_id[match(cluster_key, rep$cluster_key)]]
  members <- merge(a[, .(feature_id, Geneid, gene_name, chr, start, end, strand,
                         width, dropped_short)],
                   k[, .(feature_id, representative)], by = "feature_id", all.x = TRUE)
  members[, is_representative := feature_id %in% rep$feature_id]

  qc <- data.table(
    step  = c("unique 5'UTR features", sprintf("shorter than %d bp (dropped)", min_width),
              "features clustered", "promoter clusters", "clusters with >1 member",
              "genes", "genes with >=2 clusters"),
    value = c(nrow(a), sum(a$dropped_short), nrow(k), nrow(rep), sum(rep$n_members > 1),
              uniqueN(rep$Geneid), uniqueN(anno[n_promoters_in_gene >= 2, Geneid])))

  if (verbose) {
    message("5'UTR promoter clustering:")
    print(qc)
  }
  list(counts = mat, anno = anno, members = members, qc = qc)
}

# ---- confound removal ------------------------------------------------------
# Strips globin and mitochondrial genes BEFORE any normalisation, and returns
# the fractions they occupied so they can be reported rather than hidden.
strip_confounders <- function(mat, gene_names, verbose = TRUE) {
  stopifnot(length(gene_names) == nrow(mat))
  tot <- colSums(mat)

  is_globin <- gene_names %in% GLOBIN_GENES
  is_mt     <- is_mt_gene(gene_names)

  frac <- data.table(
    sample      = colnames(mat),
    total_raw   = tot,
    globin_pct  = round(100 * colSums(mat[is_globin, , drop = FALSE]) / tot, 3),
    mt_pct      = round(100 * colSums(mat[is_mt,     , drop = FALSE]) / tot, 3)
  )
  frac[, retained_pct := round(100 - globin_pct - mt_pct, 2)]

  keep <- !(is_globin | is_mt)
  if (verbose) {
    message("removing ", sum(is_globin), " globin and ", sum(is_mt),
            " mitochondrial features before normalisation")
    print(frac)
  }

  list(counts = mat[keep, , drop = FALSE], keep = keep, fractions = frac)
}

# ---- expression filter -----------------------------------------------------
# CPM computed on the post-strip matrix so globin does not steal library size.
filter_expressed <- function(mat, cpm_cut = THRESH$cpm,
                             min_samp = THRESH$min_samples, verbose = TRUE) {
  cpm_mat <- edgeR::cpm(mat)
  keep <- rowSums(cpm_mat > cpm_cut) >= min_samp
  if (verbose)
    message(sprintf("expression filter (CPM > %.1f in >= %d samples): %s of %s features retained",
                    cpm_cut, min_samp,
                    format(sum(keep), big.mark = ","),
                    format(length(keep), big.mark = ",")))
  keep
}
