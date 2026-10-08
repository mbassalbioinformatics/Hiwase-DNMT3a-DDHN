# Methods

## 1. Design

Two cases of donor-cell-derived haematological neoplasm (DDHN) after allogeneic
haematopoietic stem cell transplantation.

- **Case 1 (DNMT3A).** Recipient bone marrow at day 30 after transplant, at 3.5 years
  and at DDHN diagnosis, and donor samples at donation and later. TAPS methylation
  for all five; RNA-seq for day 30, 3.5 years, DDHN and the later donor sample.
- **Case 2 (TET2).** Two samples sharing one genome, the DDHN and an earlier sample.
  TAPS methylation only.

Between day 30 and 3.5 years the clone expands with no disease ("clonal
establishment"); between 3.5 years and diagnosis clone size is flat while disease
appears ("transformation"). The two axes are analysed separately. Comparisons are
pairwise and named `first_vs_second`; every difference is first minus second.

With one library per timepoint there are no replicates. p-values and FDRs therefore
rank features; confidence comes from agreement between independent comparisons, from
a **reference comparison** (day 30 against the donor: no disease and a small clone at
both ends, so it measures tissue and engraftment effects) and from recurrence in the
second case.

## 2. TAPS methylation

**Libraries** were prepared with Watchmaker TAPS+ chemistry and sequenced paired-end;
the reads of all sequencing runs of a library were combined.

**Alignment and calls.** Msuite2 v2.2.0 with
`-x hg38 -m TAPS -3 --minsize 15 --maxins 1500 --cut-r1-head 8 --cut-r2-head 8`
(three-letter alignment; 8 bases trimmed from the 5' end of each read). In TAPS a
methylated cytosine is read as T and an unmethylated one as C. Methylation at a CpG is
the share of T among C and T calls, both strands combined, from Msuite2's per-CpG base
counts (`02_taps_inputs.py`). chr1-22 and chrX are analysed. TAPS reads
5-methylcytosine and 5-hydroxymethylcytosine together.

**Sample identity** (`sample_identity_check.py`). At a CpG where a person carries a
SNP that destroys the C or the G, reads show a base neither methylation state can
produce; the share of such reads is a genotype (0, 0.5, 1) that methylation does not
affect. Candidate sites (>= 10 reads, other-base share >= 0.25 in any library) were
genotyped in every library (0: < 0.10; 1: 0.25-0.75; 2: > 0.90) and libraries compared
pairwise as the share of shared sites with the same genotype.

**Driver mutations** (`06_vaf_pileup.sh`). Base counts over the coding exons (+/- 10 bp)
of 85 myeloid and clonal-haematopoiesis driver genes, MAPQ >= 10, base quality >= 20,
read 1 and read 2 counted separately so that each base is assigned to the original top
(read 1 forward, read 2 reverse) or bottom strand. At a CpG, C>T is read on the strand
that does not carry the C.

## 3. Differentially methylated regions

**metilene** (primary): per comparison, the two samples' per-CpG methylation joined
with `bedtools unionbedg`, then `metilene -m 3 -d 0.1 -M 147 -f 1` (>= 3 CpGs, mean
difference >= 0.1, <= 147 bp between CpGs, de novo segmentation). DMRs at q < 0.05.

**Defiant** (second caller): `defiant -b -c 0 -d 1 -G 147 -CpN 5 -p 0.05 -R`. The
"both" set holds metilene DMRs overlapping a Defiant DMR changing the same way.

**Motifs.** HOMER `findMotifsGenome.pl <bed> hg38 <out> -size given -nomotif -cpg`
(known motifs, CpG-normalised background) in DMRs that lost and gained methylation,
separately, for comparisons with >= 5 DMRs in that direction.

**Annotation.** HOMER `annotatePeaks.pl <bed> hg38 -annStats`, and in
`04_integration/06_dmr_master.R`: a DMR lies in a gene's regulatory switch region when
it overlaps the promoter (TSS -1,500 to -1), exon 1 or intron 1 of the MANE Select v1.0
transcript; otherwise in the gene body or intergenic. Each DMR is linked to its
switch-region gene, else to the nearest TSS within 50 kb, and annotated with CpG-island
context and ENCODE cCRE class.

## 4. RNA-seq

**Libraries.** Stranded poly(A) libraries (Tecan/NuGEN Universal Plus mRNA-Seq),
MGI DNBSEQ-G400, paired-end 100 bp. Read 1 is in the sense orientation.

**Quantification.** featureCounts (Subread 2.1.1) against Ensembl 107:
`-p --countReadPairs -M -O -J -f -d 20 -s 1 -g gene_id`, at gene level and at
5'UTR level (`-t five_prime_utr`). With `-f` a row is an exon feature; rows were
de-duplicated on (gene, chromosome, start, end, strand) before summing per gene.
Haemoglobin and mitochondrial transcripts were removed before normalisation and
reported; features above 1 CPM in at least two libraries were tested.

**Differential expression.** TMM normalisation; edgeR with a fixed dispersion
(biological coefficient of variation 0.4, with 0.2 and 0.6 as a sensitivity band),
a `~ 0 + sample` design and likelihood-ratio tests for six contrasts. NOISeq
(no-replicate mode, q >= 0.999) as a second engine. A gene is changed at FDR < 0.05,
|log2FC| >= 1 and >= 1 CPM in one library of the pair; the headline set also requires
NOISeq agreement and excludes on/off calls (detected in one library only). Effect sizes
are log2(CPM + 1) differences. Genes are classified across the two axes (clone,
transformation, progressive, reversing), with "flat" requiring |log2FC| < 0.5.

**Composition and confounding.** Lineage content scored per library by ssGSEA (GSVA)
over curated marker sets; proliferation, oxidative phosphorylation, erythroid and
interferon signatures reported beside the results.

**Promoter usage.** Overlapping 5'UTRs of a gene merged into promoter clusters (each
represented by its longest member; 5'UTRs < 25 bp dropped); relative usage tested with
`diffSpliceDGE` at the same fixed dispersion; a shift needs FDR < 0.05, a change in
usage >= 0.15 and the promoter detected in both libraries.

**Splicing.** rMATS-turbo 4.3.0 (`--libType fr-secondstrand`, `--variable-read-length`,
`--novelSS`; numerator as `--b1`) for the six contrasts and one grouped run, and Shiba
(one run over all libraries, Fisher's exact test on junction reads). Events at
FDR < 0.05, |dPSI| >= 0.10 and >= 20 junction reads; direction recomputed from the
per-sample PSI for both tools.

**Enrichment.** clusterProfiler over-representation (GO, KEGG, MSigDB Hallmark) and
GSEA on the signed square root of the likelihood-ratio statistic, run serially with the
seed set before every call. Donor contrasts are filtered with the reference contrast:
a gene set significant in a donor contrast is "compartment" when the reference contrast
moves it the same way at raw p < 0.10.

## 5. Integration

- **Region methylation.** For each MANE Select gene: promoter (TSS -1,500 to +500),
  gene body, and the genome in 1-kb tiles; mean of CpGs with >= 5 reads, if >= 3 such
  CpGs.
- **Direction by context.** 1-kb tiles changing by >= 0.25, summarised by CpG-island
  context (island; shore <= 2 kb; shelf 2-4 kb; open sea) and ENCODE cCRE class (PLS,
  pELS, dELS, CTCF-only, DNase-H3K4me3, none). Overlap between comparisons by Fisher's
  exact test over tiles with data in all libraries.
- **Methylation and expression: switches, not correlation.** A promoter switches when it
  changes by >= 0.20 and crosses 0.5. Enrichment of expression change among switched
  genes (Fisher) and its direction (binomial); starburst plots name every gene in each
  quadrant (|change| >= 0.10, |log2FC| >= 0.5).
- **TF activity.** decoupleR univariate linear model on each contrast's edgeR
  t-statistic with the CollecTRI network; FDR per contrast; donor contrasts read against
  the reference contrast.
- **Gene panels.** Fourteen panels from MSigDB (C2, C5 GO:BP, C8, Hallmark) and curated
  lists, every gene with its source (`07_gene_panels.R`).
- **Chromatin loops.** Hi-C loops of sorted HSC, CMP, GMP and MEP (two replicates each),
  10 kb-1 Mb, pooled; partner anchors merged within 10 kb and kept when present in >= 2
  of the 8 files. For each DMR in an anchor: the partner's genes, any DMR there, and the
  most significant change metilene tested inside it.
- **Locus plots.** Gene model and CpG islands, each CpG linked to a column of a
  per-sample heatmap (every CpG with >= 1 read in any library), the DMR boxed, the
  gene's RNA-seq beside it and Hi-C partners below.

## 6. References

- Liu Y et al. Bisulfite-free direct detection of 5-methylcytosine and 5-hydroxymethylcytosine at base resolution. *Nat Biotechnol* 2019;37:424-9.
- Li L et al. Msuite2: all-in-one DNA methylation data analysis toolkit. *Comput Struct Biotechnol J* 2022;20:1271-6.
- Jühling F et al. metilene: fast and sensitive calling of differentially methylated regions. *Genome Res* 2016;26:256-62.
- Condon DE et al. Defiant: (DMRs: easy, fast, identification and ANnoTation). *BMC Bioinformatics* 2018;19:31.
- Heinz S et al. Simple combinations of lineage-determining transcription factors prime cis-regulatory elements. *Mol Cell* 2010;38:576-89.
- Chen Y et al. edgeR v4. *Nucleic Acids Res* 2025;53:gkaf018.
- Badia-i-Mompel P et al. decoupleR. *Bioinform Adv* 2022;2:vbac016.
- Müller-Dott S et al. Expanding the coverage of regulons from high-confidence prior knowledge (CollecTRI). *Nucleic Acids Res* 2023;51:10934-49.
- Morales J et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research (MANE). *Nature* 2022;604:310-5.
- ENCODE Project Consortium et al. Expanded encyclopaedias of DNA elements in the human and mouse genomes. *Nature* 2020;583:699-710.
- Liberzon A et al. The MSigDB hallmark gene set collection. *Cell Syst* 2015;1:417-25.
- Brenet F et al. DNA methylation of the first exon is tightly linked to transcriptional silencing. *PLoS One* 2011;6:e14524.
- Anastasiadi D et al. Consistent inverse correlation between DNA methylation of the first intron and gene expression. *Epigenetics Chromatin* 2018;11:37.
