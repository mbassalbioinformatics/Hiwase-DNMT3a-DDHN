#!/usr/bin/env Rscript
###############################################################################
# 17_genomewide_snvs.R                                               2026-10-15
#
# Point mutations and small insertions/deletions across the whole genome, from
# the Mutect2 run on the TAPS reads (05_variants_server/03_somatic_calls.sh, one run
# per person, roles from config/samples.tsv, column variant_role; DNMT3A case:
# normals = donor at donation and day 30, tumours = 3.5 yr, DDHN and donor 2023;
# TET2 case: normal = the earlier sample, tumour = DDHN).
#
# FilterMutectCalls assumes ordinary DNA. TAPS adds one-strand C>T / G>A changes
# (unmethylated C converted at about 0.2%), so the PASS calls are filtered further:
#   1 strand rule  a C>T needs >= 2 alt reads in F2R1 orientation (the original
#                  bottom strand), a G>A >= 2 in F1R2, in a library with >= 2 alt
#                  reads; other changes need >= 2 alt reads
#   2 a dinucleotide change with >= 90% of its alt reads in one orientation is removed
#   3 gnomAD allele frequency > 0.1% (POPAF < 3) removed
#   4 outside uniquely mappable sequence (Umap k100) or inside a segmental duplication removed
# Then: carrier patterns and VAFs; protein-changing calls (Funcotator, GENCODE 43;
# the 85-gene myeloid panel flagged); and a pooled test of whether each clone's
# mutations are already present in the earlier samples:
#   sites    SNVs with >= 3 alt reads and VAF >= 0.25 in every clone library, the
#            strand rule met, gnomAD AF <= 0.1%, median alt mapping quality >= 40,
#            VAF < 0.2 in the normal libraries; FILTER PASS, normal_artifact or
#            contamination (FilterMutectCalls marks a call normal_artifact exactly
#            when an earlier sample has alt reads, so PASS alone would hide them)
#   counts   orientation counts (base quality >= 20), informative strand for C>T / G>A
#   subsets  DNMT3A case: seen in the donor's 2023 blood (>= 2 alt reads: shared with
#            the donor's own clone) or not seen there (0 alt reads of >= 15)
#   background  the other person's libraries at the same sites, from
#            05_variants_server/06_site_counts.sh (emseq/04_variants/site_counts/) when it has run;
#            until then the TET2 case's earlier sample at its own clone sites
#
# Input   integration/variants/calls_<case>.tsv.gz (17_extract_calls.py, run from
#         here when missing or older than the VCFs in emseq/04_variants/vcf/)
# Output  integration/variants/: funnel.tsv, calls_final.tsv.gz, carriers.tsv,
#         vaf_summary.tsv, clonal_counts.tsv, coding.tsv, pooled.tsv, pooled_sites_other_person.tsv, site_counts_sites.tsv (for the server),
#         snv_vaf.png, snv_pooled.png
###############################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("config/paths.R")                          # EMS, INT, SAMPLE_SHEET
VCFD <- file.path(EMS, "04_variants", "vcf")
SCD  <- file.path(EMS, "04_variants", "site_counts")
OUT  <- file.path(INT, "variants"); dir.create(OUT, showWarnings = FALSE)
# ---- libraries, cases, labels and persons from the sample sheet (config/samples.tsv)
.sheet <- fread(SAMPLE_SHEET, colClasses = "character")[methylation == "analysed"]
ID <- setNames(.sheet$sample, .sheet$key)
CASE1 <- unique(.sheet$case)[1]; CASE2 <- unique(.sheet$case)[2]
LAB <- setNames(.sheet$short_label, .sheet$sample)
PERSON <- setNames(.sheet$person, .sheet$sample)
.cl <- .sheet[clone_library == "yes"]; .two <- .cl[, .N, by = case][N == 2, case]
.c2 <- .cl[case == .two, short_label]; .c1 <- .cl[case != .two, short_label][1]
CL <- list(both = paste0("clone: ", .c2[1], " and ", .c2[2]), a = paste(.c2[1], "only"), b = paste(.c2[2], "only"),
           donor = paste(.sheet[donor_clone == "yes", short_label], "only"),
           one = paste0("clone: ", .c1, ", VAF >= 0.25"), one_low = paste0(.c1, ", VAF < 0.25"))    # carrier classes
# ----
CASES <- c(CASE1, CASE2)
CASE_COL <- setNames(c("#6A51A3", "#D95F02"), paste(CASES, "case"))
DONOR_CLONE <- .sheet[donor_clone == "yes", sample]   # the library that holds the donor's own clone: splits the clone's mutations into shared and branch-only
OKF <- c("PASS", "normal_artifact", "contamination")      # FILTER states allowed in the pooled test
MIN_POPAF <- 3; MIN_MMQ <- 40; CLONE_ALT <- 3L; CLONE_VAF <- 0.25; NORMAL_MAX_VAF <- 0.2
CODING <- c("MISSENSE", "NONSENSE", "NONSTOP", "FRAME_SHIFT_INS", "FRAME_SHIFT_DEL", "IN_FRAME_INS", "IN_FRAME_DEL",
            "SPLICE_SITE", "START_CODON_SNP", "START_CODON_INS", "START_CODON_DEL")
PANEL <- unique(fread("01_methylation_server/vaf_panel_hg38.bed", header = FALSE)$V4)

## ---- the call tables (extracted again when the VCFs are newer) ---------------------------
vcfs <- list.files(VCFD, "\\.vcf\\.gz$", full.names = TRUE)
tabs <- file.path(OUT, paste0("calls_", CASES, ".tsv.gz"))
if (!all(file.exists(tabs)) || min(file.mtime(tabs)) < max(file.mtime(vcfs))) {
  st <- system2("python3", c(shQuote(file.path("04_integration", "17_extract_calls.py")), shQuote(VCFD), shQuote(OUT)))
  if (st != 0) stop("17_extract_calls.py failed")
}
read_case <- function(case) {
  f <- file.path(OUT, paste0("calls_", case, ".tsv.gz"))
  meta <- system(sprintf("gzip -dc '%s' | head -1", f), intern = TRUE)
  roles <- lapply(strsplit(sub("^#", "", meta), ";")[[1]], function(kv) strsplit(sub(".*=", "", kv), ",")[[1]])
  names(roles) <- sub("=.*", "", strsplit(sub("^#", "", meta), ";")[[1]])
  d <- fread(cmd = sprintf("gzip -dc '%s' | tail -n +2", f), colClasses = list(character = c("chr", "ref", "alt", "filter", "gene", "vclass", "protein")))
  d[, `:=`(case = case, id = .I)]
  setnames(d, "alt", "alt_allele")
  libs <- c(roles$normals, roles$tumours)
  long <- rbindlist(lapply(libs, function(s) d[, .(id, sample = s, a = get(paste0(s, "_alt")), n = get(paste0(s, "_dp")),
                                                  a1 = get(paste0(s, "_f1")), n1 = get(paste0(s, "_f1_dp")),
                                                  a2 = get(paste0(s, "_f2")), n2 = get(paste0(s, "_f2_dp")))]))
  long <- merge(long, d[, .(id, ref, alt_allele)], by = "id")
  long[, `:=`(ct = ref == "C" & alt_allele == "T", ga = ref == "G" & alt_allele == "A")]
  # informative reads: the strand TAPS cannot convert for C>T / G>A, both orientations otherwise
  long[, `:=`(ia = fifelse(ct, a2, fifelse(ga, a1, a1 + a2)), idp = fifelse(ct, n2, fifelse(ga, n1, n1 + n2)),
              sr = fifelse(ct, a2, fifelse(ga, a1, a)))]                  # strand-rule count (all reads for other changes)
  long[, c("ref", "alt_allele", "ct", "ga") := NULL]
  long[, `:=`(role = fifelse(sample %in% roles$normals, "normal", "tumour"), vaf = fifelse(n > 0, a / n, NA_real_))]
  d[, (grep(paste0("^(", paste(libs, collapse = "|"), ")_"), names(d), value = TRUE)) := NULL]   # per-library columns now in long
  list(d = d, long = long, roles = roles)
}
X <- setNames(lapply(CASES, read_case), CASES)
CASE_PERSON <- vapply(CASES, function(cs) PERSON[[X[[cs]]$roles$tumours[1]]], "")
DC <- CASES[vapply(CASES, function(cs) DONOR_CLONE %in% X[[cs]]$roles$tumours, TRUE)]   # the case of the donor-clone library

## ---- the filter funnel and the final calls ----------------------------------------------
fin <- list(); funnel <- list()
for (case in CASES) {
  d <- X[[case]]$d; L <- X[[case]]$long
  p <- d[pass == 1]
  tl <- L[role == "tumour" & id %in% p$id]
  k1 <- tl[a >= 2 & sr >= 2, unique(id)]
  ori <- tl[, .(F1 = sum(a1), F2 = sum(a2)), by = id]
  bad_mnv <- merge(p[type == "MNV", .(id)], ori, by = "id")[pmax(F1, F2) >= 0.9 * (F1 + F2), id]
  p[, popaf := fifelse(is.na(popaf), 6, popaf)]
  s1 <- p[id %in% k1]; s2 <- s1[!id %in% bad_mnv]; s3 <- s2[popaf >= MIN_POPAF]; s4 <- s3[umap100 == 1 & segdup == 0]
  funnel[[case]] <- data.table(case = case, step = c("PASS (FilterMutectCalls)", "strand rule", "dinucleotides seen in one orientation removed",
                                                     "gnomAD allele frequency > 0.1% removed", "outside unique mappability or in a segmental duplication removed"),
                               calls = c(nrow(p), nrow(s1), nrow(s2), nrow(s3), nrow(s4)))
  f <- copy(s4)
  w <- dcast(L[id %in% f$id, .(id, sample, cell = paste0(a, "/", n))], id ~ sample, value.var = "cell")
  car <- L[id %in% f$id & role == "tumour" & a >= 2, .(carriers = paste(sample, collapse = "+")), by = id]
  f <- merge(merge(f, w, by = "id"), car, by = "id", all.x = TRUE)
  f[is.na(carriers), carriers := ""]
  cl <- X[[case]]$roles$clone
  vv <- dcast(L[id %in% f$id & sample %in% c(cl, X[[case]]$roles$tumours), .(id, sample, vaf)], id ~ sample, value.var = "vaf")
  setnames(vv, setdiff(names(vv), "id"), paste0("vaf_", setdiff(names(vv), "id")))
  f <- merge(f, vv, by = "id")
  if (length(cl) == 2) {
    f[, class := fcase(grepl(cl[1], carriers) & grepl(cl[2], carriers), CL$both,
                       grepl(cl[1], carriers), CL$a, grepl(cl[2], carriers), CL$b,
                       grepl(DONOR_CLONE, carriers), CL$donor, default = "none with >= 2 reads")]
  } else {
    f[, class := fcase(grepl(cl[1], carriers) & get(paste0("vaf_", cl[1])) >= CLONE_VAF, CL$one,
                       grepl(cl[1], carriers), CL$one_low, default = "none with >= 2 reads")]
  }
  fin[[case]] <- f
}
funnel <- rbindlist(funnel)
fwrite(funnel, file.path(OUT, "funnel.tsv"), sep = "\t")
print(funnel)

## ---- carriers, VAFs, protein-changing calls -------------------------------------------------
carriers <- rbindlist(lapply(CASES, function(case) fin[[case]][, .(calls = .N, snv = sum(type == "SNV"), indel = sum(type == "indel"),
                                                                     mnv = sum(type == "MNV")), by = class][, case := case]))
setcolorder(carriers, "case")
fwrite(carriers, file.path(OUT, "carriers.tsv"), sep = "\t")
print(carriers)
vs <- rbindlist(lapply(CASES, function(case) {
  L <- X[[case]]$long
  L[id %in% fin[[case]]$id & role == "tumour" & a >= 2, .(calls = .N, median_vaf = round(median(vaf), 3), vaf_ge_025 = sum(vaf >= 0.25),
                                                          median_depth = as.double(median(n))), by = sample][, case := case]
}))
fwrite(vs, file.path(OUT, "vaf_summary.tsv"), sep = "\t")
coding <- rbindlist(lapply(CASES, function(case) {
  f <- fin[[case]][vclass %in% CODING]
  libs <- c(X[[case]]$roles$normals, X[[case]]$roles$tumours)
  f[, .(case, gene, vclass, protein, change = paste0(chr, ":", pos, " ", ref, ">", alt_allele), class, panel = gene %in% PANEL,
        reads = do.call(paste, c(lapply(libs, function(s) paste0(LAB[s], " ", get(s))), sep = "; ")))]
}))
setorder(coding, case, -panel, class, gene)
fwrite(coding, file.path(OUT, "coding.tsv"), sep = "\t")
cat("\nprotein-changing calls:", coding[, .N, by = case][, paste(case, N, collapse = ", "), ], "; in the 85-gene panel:",
    coding[panel == TRUE, paste(case, gene, protein, collapse = "; ")], "\n")
fwrite(rbindlist(lapply(CASES, function(case) fin[[case]][, .(case, chr, pos, ref, alt = alt_allele, type, class, carriers, gene, vclass, protein)])),
       file.path(OUT, "calls_final.tsv.gz"), sep = "\t")

## ---- pooled test: the clone's mutations in the earlier samples ------------------------------
pool_sites <- rbindlist(lapply(CASES, function(case) {
  d <- X[[case]]$d; L <- X[[case]]$long; r <- X[[case]]$roles
  ok <- d[type == "SNV" & fifelse(is.na(popaf), 6, popaf) >= MIN_POPAF & !is.na(mmq_alt) & mmq_alt >= MIN_MMQ &
          vapply(strsplit(filter, ";"), function(x) all(x %in% OKF), TRUE), id]
  germ <- L[id %in% ok & role == "normal" & n > 0 & vaf >= NORMAL_MAX_VAF, unique(id)]
  ok <- setdiff(ok, germ)
  big <- L[id %in% ok & a >= CLONE_ALT & n > 0 & vaf >= CLONE_VAF & ia >= 1]
  clone <- big[sample %in% r$clone, .N, by = id][N == length(r$clone), id]
  out <- d[id %in% clone, .(case, id, chr, pos, ref, alt = alt_allele, set = paste0(case, "_clone"))]
  if (case == DC) {
    dn <- L[id %in% clone & sample == DONOR_CLONE]
    out[id %in% dn[a >= 2, id], set := paste0(case, "_clone_in_donor")]
    out[id %in% dn[a == 0 & n >= 15, id], set := paste0(case, "_clone_not_in_donor")]
    nob <- L[id %in% ok & sample %in% r$clone & a == 0 & n >= 10]
    for (s in r$clone) {
      o <- setdiff(r$clone, s)
      only <- intersect(big[sample == s, id], nob[sample == o, id])
      out <- rbind(out, d[id %in% only, .(case, id, chr, pos, ref, alt = alt_allele, set = paste0(case, "_", s, "_only"))])
    }
  }
  out
}))
fwrite(pool_sites[, .(set, chr, pos, ref, alt)], file.path(OUT, "site_counts_sites.tsv"), sep = "\t")
cat("\npooled-test sites:\n"); print(pool_sites[, .N, by = set])
binci <- function(x, n) if (n > 0) binom.test(x, n)$conf.int else c(NA_real_, NA_real_)
# server read counts at the same sites in every library (05_variants_server/06_site_counts.sh). A change that shows in the
# other person's libraries (>= 2 alt reads, any orientation) cannot be the clone's mutation: it is a systematic
# artefact or an inherited polymorphism, and the site leaves the test. What the other person's libraries show at
# the remaining sites is the error background.
sc_file <- file.path(SCD, "site_counts.tsv")
pooled_sc <- NULL; seen <- NULL
if (file.exists(sc_file)) {
  sc <- fread(sc_file, colClasses = list(character = c("set", "sample", "chr", "ref", "alt")))
  sc[, `:=`(person = PERSON[sample], set_person = CASE_PERSON[sub("_.*", "", set)],
            depth_all = f1_ref + f1_alt + f1_other + f2_ref + f2_alt + f2_other)]
  seen <- sc[person != set_person, .(other_alt = sum(f1_alt + f2_alt), other_depth = sum(depth_all)), by = .(set, chr, pos, ref, alt)]
  seen[, seen_in_other := other_alt >= 2]
  fwrite(seen, file.path(OUT, "pooled_sites_other_person.tsv"), sep = "\t")
  sc <- merge(sc, seen[, .(set, chr, pos, ref, alt, seen_in_other)], by = c("set", "chr", "pos", "ref", "alt"))[seen_in_other == FALSE]
  sc[, `:=`(ct = ref == "C" & alt == "T", ga = ref == "G" & alt == "A")]
  sc[, `:=`(ia = fifelse(ct, f2_alt, fifelse(ga, f1_alt, f1_alt + f2_alt)),
            idp = fifelse(ct, f2_ref + f2_alt + f2_other, fifelse(ga, f1_ref + f1_alt + f1_other, depth_all)))]
  sc_all <- rbind(sc, copy(sc[grepl(paste0("^", DC, "_clone"), set)])[, set := paste0(DC, "_clone_all")])
  pooled_sc <- sc_all[, .(sites = .N, alt = sum(ia), depth = sum(idp)), by = .(set, sample)][, source := "site counts (server)"]
  dropped <- merge(pool_sites, seen[seen_in_other == TRUE, .(set, chr, pos, ref, alt)], by = c("set", "chr", "pos", "ref", "alt"))
  pool_sites[, kept := !id %in% dropped$id | case != dropped$case[match(id, dropped$id)]]
  cat("\nsites seen in the other person's libraries (left out of the pooled test):\n"); print(seen[, .(sites = .N, seen = sum(seen_in_other)), by = set])
} else pool_sites[, kept := TRUE]
pooled_vcf <- rbindlist(lapply(CASES, function(cc) {
  case <- cc; L <- X[[cc]]$long; ps <- pool_sites[pool_sites$case == cc & kept == TRUE]
  m <- merge(L[id %in% ps$id], ps[, .(id, set)], by = "id")
  sets <- unique(c(ps$set, if (case == DC) paste0(DC, "_clone_all")))
  rbindlist(lapply(sets, function(st) {
    mm <- if (st == paste0(DC, "_clone_all")) m[grepl(paste0("^", DC, "_clone"), set)] else m[set == st]
    mm[, .(sites = uniqueN(id), alt = sum(ia), depth = sum(idp)), by = sample][, `:=`(set = st, source = "Mutect2 orientation counts")]
  }))
}))
pooled <- rbind(pooled_vcf, pooled_sc, fill = TRUE)
pooled[, c("ci_low", "ci_high") := as.list(binci(alt, depth)), by = seq_len(nrow(pooled))]
pooled[, `:=`(vaf = alt / depth, label = LAB[sample], person = PERSON[sample])]
pooled[, set_person := CASE_PERSON[sub("_.*", "", set)]]
pooled[, background := person != set_person]                      # the other person's libraries: error background
setcolorder(pooled, c("set", "source", "sample", "label", "sites", "alt", "depth", "vaf", "ci_low", "ci_high", "background"))
fwrite(pooled[, -c("person", "set_person")], file.path(OUT, "pooled.tsv"), sep = "\t")
cat("\npooled informative-strand VAF (%):\n")
print(dcast(pooled, set + source ~ label, value.var = "vaf", fun.aggregate = function(x) if (length(x)) round(100 * x[1], 2) else NA_real_))

## ---- clonal calls: VAF >= 0.25 in every clone library, not seen in the other person -----------
clonal <- rbindlist(lapply(CASES, function(case) {
  f <- copy(fin[[case]]); cl <- X[[case]]$roles$clone
  f[, big := Reduce(`&`, lapply(cl, function(s) get(paste0("vaf_", s)) >= CLONE_VAF & grepl(s, carriers)))]
  if (!is.null(seen)) {
    sx <- seen[grepl(paste0("^", case), set) & seen_in_other == TRUE, .(chr, pos, ref, alt_allele = alt)]
    f[sx, on = .(chr, pos, ref, alt_allele), seen := TRUE]
  }
  if (!"seen" %in% names(f)) f[, seen := NA]
  cs <- case
  f[big == TRUE, .(case = cs, clonal = sum(is.na(seen)), snv = sum(is.na(seen) & type == "SNV"), indel = sum(is.na(seen) & type == "indel"),
                   mnv = sum(is.na(seen) & type == "MNV"), seen_in_other_person = sum(!is.na(seen)))]
}))
fwrite(clonal, file.path(OUT, "clonal_counts.tsv"), sep = "\t")
cat("\nclonal calls (VAF >= 0.25 in every clone library):\n"); print(clonal)

## ---- figure: VAFs --------------------------------------------------------------------------
fa <- fin[[DC]][class %in% c(CL$both, CL$a, CL$b)]
fa[, class := factor(class, levels = c(CL$both, CL$a, CL$b))]
cl2 <- X[[DC]]$roles$clone
pa <- ggplot(fa, aes(.data[[paste0("vaf_", cl2[1])]], .data[[paste0("vaf_", cl2[2])]], colour = class)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(size = 0.9, alpha = 0.5) +
  scale_colour_manual(values = setNames(c("#6A51A3", "#1B9E77", "#E7298A"), c(CL$both, CL$a, CL$b)), name = NULL) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "VAF at 3.5 years", y = "VAF at DDHN", title = "DNMT3A case: each mutation at 3.5 years and at DDHN",
       subtitle = "final calls with >= 2 alt reads in either sample") +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 9.5))
hv <- rbindlist(lapply(CASES, function(case) X[[case]]$long[id %in% fin[[case]]$id & role == "tumour" & a >= 2, .(sample, vaf)][, case := paste(case, "case")]))
hv[, lab := factor(LAB[sample], levels = LAB)]
pb <- ggplot(hv, aes(vaf, fill = case)) +
  geom_histogram(binwidth = 0.025, boundary = 0, colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0.5, colour = "grey50", linetype = "dotted") +
  facet_wrap(~ lab, ncol = 1, scales = "free_y", strip.position = "right") +
  scale_fill_manual(values = CASE_COL, guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "VAF (share of reads with the mutation)", y = "mutations", title = "VAF of every final call, per library",
       subtitle = "calls with >= 2 alt reads in that library; dotted line: 50%") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0), plot.subtitle = element_text(size = 9.5))
ggsave(file.path(OUT, "snv_vaf.png"), pa + pb + plot_layout(widths = c(1.15, 1)), width = 14, height = 7.5, dpi = 200)

## ---- figure: the pooled test -----------------------------------------------------------------
# A set defined by the donor's 2023 reads says nothing about the donor's 2023 library itself, so it is not shown there.
OC <- setdiff(CASES, DC)            # the other case(s): one panel each
EARLY <- c(X[[DC]]$roles$normals, DONOR_CLONE, unlist(lapply(OC, function(cs) X[[cs]]$roles$normals)))
PANELS <- data.table(set = c(paste0(DC, c("_clone_all", "_clone_in_donor", "_clone_not_in_donor")), paste0(OC, "_clone")),
                     title = c(paste(DC, c("clone: all its mutations", "clone: also in the donor's 2023 blood",
                                           "clone: not seen in the donor's 2023 blood")), paste(OC, "clone: all its mutations")),
                     libs = c(list(c(X[[DC]]$roles$normals, DONOR_CLONE), X[[DC]]$roles$normals, X[[DC]]$roles$normals),
                              lapply(OC, function(cs) X[[cs]]$roles$normals)))
src <- if (!is.null(pooled_sc)) "site counts (server)" else "Mutect2 orientation counts"
pp <- rbindlist(lapply(seq_len(nrow(PANELS)), function(i) {
  x <- pooled[source == src & set == PANELS$set[i]]
  own <- x[sample %in% PANELS$libs[[i]], .(what = LAB[sample], alt, depth)]
  bg <- x[background == TRUE, .(what = "other person's libraries (background)", alt = sum(alt), depth = sum(depth))]
  n <- pool_sites[kept == TRUE & (set == PANELS$set[i] | (PANELS$set[i] == paste0(DC, "_clone_all") & grepl(paste0("^", DC, "_clone"), set))), .N]
  rbind(own, if (nrow(bg) && bg$depth > 0) bg)[, panel := sprintf("%s\n(%s sites)", PANELS$title[i], format(n, big.mark = ","))]
}))
pp[, c("lo", "hi") := as.list(binci(alt, depth)), by = seq_len(nrow(pp))]
pp[, `:=`(vaf = alt / depth, panel = factor(panel, levels = unique(panel)),
          what = factor(what, levels = c(LAB[EARLY], "other person's libraries (background)")))]
FLOOR <- 1e-4                                                     # 0.01%: a zero is drawn on the floor of the log scale
pc <- ggplot(pp, aes(what, 100 * pmax(vaf, FLOOR), colour = what)) +
  geom_hline(yintercept = 100 * FLOOR, colour = "grey85") +
  geom_errorbar(aes(ymin = 100 * pmax(lo, FLOOR), ymax = 100 * pmax(hi, FLOOR)), width = 0.2) +
  geom_point(size = 2.6) +
  geom_text(aes(y = 100 * pmax(hi, FLOOR), label = paste0(alt, "/", prettyNum(depth, big.mark = ","))), vjust = -0.6, size = 3, colour = "grey25") +
  facet_wrap(~ panel, nrow = 1, scales = "free_x") +
  scale_y_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1, 3, 10), labels = function(v) paste0(v, "%"), expand = expansion(mult = c(0.05, 0.12))) +
  scale_colour_manual(values = setNames(c(rep_len(c("#6A51A3", "#9E9AC8", "#3F007D"), length(EARLY) - length(OC)), rep("#D95F02", length(OC)), "grey45"),
                                        c(LAB[EARLY], "other person's libraries (background)")), guide = "none") +
  labs(x = NULL, y = "share of reads carrying the clone's mutations (pooled, log scale)",
       title = "Are the clones' mutations already present in the earlier samples?",
       subtitle = paste0("Each point pools every site of the set (alt reads / reads on the informative strand); bars: 95% interval. ",
                         if (src == "site counts (server)") paste0("Grey: the same sites in the other person's libraries, where the mutations cannot be.\n",
                                                               "Sites where the other person's libraries show the change (artefacts) are left out.") else
                           "Background (the other person's libraries) not yet measured.")) +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid.minor = element_blank(),
                                   plot.subtitle = element_text(size = 9.5))
ggsave(file.path(OUT, "snv_pooled.png"), pc, width = 14, height = 6.5, dpi = 200)
message("written to ", OUT, if (is.null(pooled_sc)) " (site counts from the server not yet available: background from the TET2 case only)" else "")
