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

**DNA** was extracted from bone-marrow or peripheral-blood mononuclear cells (or, for one
donor sample, haematopoietic progenitor cells collected from peripheral blood) with the
Qiagen QIAamp DNA Blood Mini Kit; where RNA-seq was also done, DNA and RNA came from the
same cells.

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

**Chimerism** (`04_integration/15_chimerism.R`). At identity SNP sites where the two people
of the study are opposite homozygotes, reads carrying the allele the donor lacks were counted
in every library (donor genotype called leave-one-out); the pure donor libraries give the
error rate, and the recipient share is about twice the excess over it.

**Library depth** (`cpg_depth_summary.py`). Per library, the CpGs on chr1-22 and chrX with
>= 5 informative reads (the depth every methylation analysis requires of a CpG) and
their mean methylation.

**Driver mutations** (`06_vaf_pileup.sh`, `04_integration/13_vaf_from_taps_reads.R`).
Base counts over the coding exons (+/- 10 bp) of 85 myeloid and clonal-haematopoiesis
driver genes (samtools mpileup, MAPQ >= 10, base quality >= 20), read 1 and read 2
counted separately so that each base is assigned to the original top (read 1 forward,
read 2 reverse) or bottom strand. TAPS turns a methylated C into T on the strand that
carries it, so on the top strand a T read can be a methylated C and on the bottom strand
an A read (top-strand view) a methylated G. Each allele is counted only where conversion
can neither create nor hide it: alt A on the top strand, and on the bottom strand unless
a G allele there could be a CpG G; alt T on the bottom strand, and on the top strand
unless a C allele could be a CpG C; alt C on the bottom strand, and on the top strand
unless the next base can be G; alt G on the top strand, and on the bottom strand unless
the previous base can be C ("can": from the reference, or an allele seen on its clean
strand in any library, which catches a SNP that creates a CpG next door). VAF = alt
reads / reads on the strands used. A library carries a variant with >= 3 alt reads,
VAF >= 10% and >= 8 reads. Within a person (same genome), a variant whose VAF differs
between libraries (Fisher's exact test, p < 0.01) is acquired, or allelic imbalance when
it is a germline SNP elsewhere, a known common SNP, or reaches >= 85%. Flags: alt reads
on one strand only (none on a strand with >= 5 reads, or Fisher p < 0.05), a
homopolymer >= 6 or di-/trinucleotide repeat x 4 within 10 bp, another variant within
10 bp, more than one alt allele, one library with only 3 alt reads. Changes seen in a
single library are listed with their per-strand support. 25 hotspot codons are checked
in every library. Protein changes come from the MANE Select v1.5 CDS (checked: DNMT3A
chr2:25,234,373 C>T = c.2645G>A p.Arg882His). Insertions and deletions are not counted.

**Allelic balance** (`allele_balance_counts.py`, `04_integration/14_allelic_imbalance.R`).
At the transversion CpG-SNPs of the identity check, the share of "other" bases on each
strand is the B-allele fraction (BAF). Heterozygous sites are defined per person from
the pooled libraries (BAF 0.2-0.8, >= 12 reads); per library, sites with >= 6 reads on
the strand. In 10-Mb bins (>= 15 sites; acrocentric short arms excluded) the mean
|BAF - 0.5| is compared with its exact binomial expectation for a balanced heterozygote
at the same depths; excess > 0.10 with z > 4 marks allelic imbalance (LOH or a
copy-number change in most cells).

## 3. Differentially methylated regions

**metilene** (primary): per comparison, the two samples' per-CpG methylation joined
with `bedtools unionbedg`, then `metilene -m 3 -d 0.1 -M 147 -f 1` (>= 3 CpGs, mean
difference >= 0.1, <= 147 bp between CpGs, de novo segmentation). DMRs at q < 0.05.

**Defiant** (second caller): `defiant -b -c 0 -d 1 -G 147 -CpN 5 -p 0.05 -R`. The
"both" set holds metilene DMRs overlapping a Defiant DMR changing the same way.

**Motifs.** HOMER v4.11.1 `findMotifsGenome.pl <bed> hg38 <out> -size given -nomotif -cpg`
(known motifs, CpG-normalised background) in DMRs that lost and gained methylation,
separately, for comparisons with >= 5 DMRs in that direction.

**Annotation.** HOMER v4.11.1 `annotatePeaks.pl <bed> hg38 -annStats` (one category per DMR by
priority, promoter-TSS = TSS -1 kb to +100 bp; observed / expected per category;
summarised with UCSC CpG-island overlap in `04_integration/11_homer_summary.R`), and in
`04_integration/06_dmr_master.R`: a DMR lies in a gene's regulatory switch region when
it overlaps the promoter (TSS -1,500 to -1), exon 1 or intron 1 of the MANE Select v1.5
transcript; otherwise in the gene body or intergenic. Each DMR is linked to its
switch-region gene, else to the nearest TSS within 50 kb, and annotated with CpG-island
context and ENCODE cCRE class.

## 4. RNA-seq

**Libraries.** Stranded poly(A) libraries (Tecan/NuGEN Universal Plus mRNA-Seq),
MGI DNBSEQ-G400, paired-end 100 bp. Read 1 is in the sense orientation.

**Alignment.** STAR 2.7.10b, to GRCh38 as distributed in the 10x Genomics Cell Ranger
reference 3.0.0 (alignment run outside this repository; the analysis starts from the BAMs).

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
  changes by >= 0.20 (20 points), up or down, at any level: each gene turns on or off at
  its own methylation level, so no fixed level has to be crossed. Enrichment of
  expression change among switched genes (Fisher) and its direction (binomial), repeated
  with stricter definitions (>= 0.30; >= 0.20 crossing 0.5; either) as a sensitivity
  check; starburst plots name every gene in each quadrant (|change| >= 0.10,
  |log2FC| >= 0.5).
- **TF activity.** decoupleR univariate linear model on each contrast's edgeR
  t-statistic with the CollecTRI network; FDR per contrast; donor contrasts read against
  the reference contrast.
- **Gene panels.** Fourteen panels from MSigDB (C2, C5 GO:BP, C8, Hallmark) and curated
  lists, every gene with its source (`07_gene_panels.R`).
- **Chromatin loops.** Hi-C of sorted HSC/MPP, CMP, GMP and MEP (two replicates each; GEO
  GSE231486, Subramanian 2023), reprocessed into loop calls; loops of 10 kb-1 Mb, pooled; partner anchors merged within 10 kb and kept when present in >= 2
  of the 8 files. For each DMR in an anchor: the partner's genes, any DMR there, and the
  most significant change metilene tested inside it.
- **Locus plots.** Gene models (overlapping genes on separate lanes) and CpG islands, each
  CpG linked to a column of a per-sample heatmap (every CpG with >= 2 reads in any library),
  the DMR boxed, the gene's RNA-seq beside it and Hi-C partners below. Window: the DMR
  +/- 1.5 kb, extended to the gene's TSS when it lies within 5 kb.

## 6. Genome-wide variants from the TAPS reads

Run on the server with `05_variants_server/` (GATK 4.6.2.0 in Docker; `run_variants.sh`
runs steps 02-05) and analysed in `04_integration/16_genomewide_cnv.R` and
`17_genomewide_snvs.R`. Library roles come from `config/samples.tsv` (`variant_role`,
`cnv_panel`, `clone_library`, `donor_clone`).

**Masked alignments** (`02_prepare_bams.sh`, `taps_mask_methylation.py`). At a reference
CpG the C is read by the original top strand and the G by the original bottom strand; in
TAPS these bases report methylation, not genotype, and would look like C>T / G>A variants
everywhere. In the Msuite2 alignments they were set to N with base quality 0, so the
callers take the genotype at a CpG from the other strand; secondary and supplementary
alignments were dropped and duplicates marked (GATK MarkDuplicates). Msuite2's aligner
gives a mapping quality of at most 42, below several GATK defaults (20-50), so every
step used mapping quality >= 10.

**Point mutations** (`03_somatic_calls.sh`, `05_annotate_collect.sh`;
`04_integration/17_extract_calls.py`, `17_genomewide_snvs.R`). Mutect2, one multi-sample
run per case in 25-Mb chunks: the libraries with little or none of the clone as normals
(in the DNMT3A case the donation sample and day 30, so that differences of material
between the progenitor-cell graft and mononuclear cells cancel), the person's other
libraries as tumours; 1000 Genomes panel of normals and gnomAD germline resource;
read-orientation model (LearnReadOrientationModel); contamination of each library
against the case's first normal (GetPileupSummaries, CalculateContamination);
FilterMutectCalls. TAPS also converts about 0.2% of unmethylated C outside CpGs, on one
strand only, so the PASS calls were filtered further: (1) strand rule: a C>T needs >= 2
alt reads in F2R1 orientation (original bottom strand) and a G>A >= 2 in F1R2, in a
library with >= 2 alt reads; other changes need >= 2 alt reads; (2) a dinucleotide change
with >= 90% of its alt reads in one orientation is removed; (3) gnomAD allele frequency
> 0.1% (POPAF < 3) is removed; (4) calls outside uniquely mappable sequence (Umap k100)
or inside segmental duplications are removed. Calls were annotated with Funcotator
(GENCODE 43; data sources v1.8). A call is clonal when it has >= 2 alt reads and
VAF >= 0.25 in every library in which the clone is large (`clone_library`) and is not
seen in the other person.

**Pooled clone test** (`06_site_counts.sh`, `taps_site_counts.py`, `17_genomewide_snvs.R`).
Whether a clone's mutations were already present in earlier samples, at a level no
single site can show, was tested by pooling sites: SNVs with >= 3 alt reads, VAF >= 0.25
and alt reads on the informative strand in every clone library, gnomAD AF <= 0.1%,
median alt mapping quality >= 40 and VAF < 0.2 in every normal library; FILTER PASS,
normal_artifact or contamination (FilterMutectCalls marks a call normal_artifact when a
normal has alt reads, so PASS alone would hide the mutations sought). Reads at these
sites were counted in every library by orientation (base quality >= 20, mapping quality
>= 10, duplicates skipped, overlapping mates counted once), on the informative strand
for C>T and G>A. The other person's libraries, which cannot carry these mutations, give
the error background; sites with >= 2 alt reads there were removed as systematic
artefacts. In the DNMT3A case the clone's mutations were also split by the donor's later
blood (`donor_clone`): seen there (>= 2 alt reads; shared with the donor's own clone) or
not (0 alt reads of >= 15; branch-only), and by presence in one clone library but not
the other (0 alt reads of >= 10). Per library: alt reads / depth summed over the sites,
with an exact binomial 95% interval; for heterozygous mutations the share of cells is
about twice the share of reads.

**Copy number and allelic balance** (`04_cnv_loh.sh`, `taps_snp_sites.py`,
`04_integration/16_genomewide_cnv.R`). GATK somatic copy-number workflow: read counts in
1-kb bins, kept when mappable (Umap k100 >= 0.9), with segmental-duplication content
<= 0.5 and GC 0.1-0.9 (FilterIntervals); denoised against a panel of normals
(`cnv_panel`; each panel library against the others); allelic counts (base quality
>= 20, so the masked bases are ignored) at common gnomAD SNPs (AF 0.1-0.9) without SNPs
whose alternative allele creates a CpG (TAPS reads that allele's methylated C as T on one
strand); ModelSegments with the heterozygous sites of the person's first normal library;
CallCopyRatioSegments. A called gain or loss is inherited (or technical) when a normal
library of the same person departs from copy ratio 1 in the same direction by at least
half as much (mean copy ratio of the segment's bins, not of log2 values, since a bin
without reads has log2 near -30), and acquired otherwise; calls starting on the short arm
of an acrocentric chromosome are not interpreted. Allelic imbalance per library: mean
|BAF - 0.5| at heterozygous sites against its exact binomial expectation at the same
depths, and segments with >= 20 heterozygous sites and minor-allele fraction < 0.40,
compared with the person's normal libraries. For the gain inside KMT2A, the stretch called
in every clone library was measured in every library (mean log2 copy ratio of its bins,
95% interval) against 60 kb of flanking bins, and its exons taken from the MANE Select
transcript.

## 7. References

- Liu Y et al. Bisulfite-free direct detection of 5-methylcytosine and 5-hydroxymethylcytosine at base resolution. *Nat Biotechnol* 2019;37:424-9.
- Li L et al. Msuite2: all-in-one DNA methylation data analysis toolkit. *Comput Struct Biotechnol J* 2022;20:1271-6.
- Jühling F et al. metilene: fast and sensitive calling of differentially methylated regions. *Genome Res* 2016;26:256-62.
- Condon DE et al. Defiant: (DMRs: easy, fast, identification and ANnoTation). *BMC Bioinformatics* 2018;19:31.
- Heinz S et al. Simple combinations of lineage-determining transcription factors prime cis-regulatory elements. *Mol Cell* 2010;38:576-89.
- Danecek P et al. Twelve years of SAMtools and BCFtools. *GigaScience* 2021;10:giab008.
- Chen Y et al. edgeR v4. *Nucleic Acids Res* 2025;53:gkaf018.
- Badia-i-Mompel P et al. decoupleR. *Bioinform Adv* 2022;2:vbac016.
- Müller-Dott S et al. Expanding the coverage of regulons from high-confidence prior knowledge (CollecTRI). *Nucleic Acids Res* 2023;51:10934-49.
- Morales J et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research (MANE). *Nature* 2022;604:310-5.
- ENCODE Project Consortium et al. Expanded encyclopaedias of DNA elements in the human and mouse genomes. *Nature* 2020;583:699-710.
- Liberzon A et al. The MSigDB hallmark gene set collection. *Cell Syst* 2015;1:417-25.
- Brenet F et al. DNA methylation of the first exon is tightly linked to transcriptional silencing. *PLoS One* 2011;6:e14524.
- Anastasiadi D et al. Consistent inverse correlation between DNA methylation of the first intron and gene expression. *Epigenetics Chromatin* 2018;11:37.
- Subramanian S et al. Genome-wide transcription factor-binding maps reveal cell-specific changes in the regulatory architecture of human HSPCs. *Blood* 2023;142:1448-62. (GEO GSE231486)
- McKenna A et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res* 2010;20:1297-303.
- Cibulskis K et al. Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples. *Nat Biotechnol* 2013;31:213-9.
- Benjamin D et al. Calling somatic SNVs and indels with Mutect2. *bioRxiv* 2019; doi:10.1101/861054.
- Karczewski KJ et al. The mutational constraint spectrum quantified from variation in 141,456 humans (gnomAD). *Nature* 2020;581:434-43.
- Karimzadeh M et al. Umap and Bismap: quantifying genome and methylome mappability. *Nucleic Acids Res* 2018;46:e120.
