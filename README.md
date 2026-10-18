# Hiwase-DNMT3a-DDHN

DNA methylation (TAPS) and RNA-seq of donor-cell-derived haematological neoplasms.

Analysis code for a study of donor-cell-derived haematological neoplasms (DDHN)
after allogeneic haematopoietic stem cell transplantation. Serial samples from two
cases, a DNMT3A case and a TET2 case, were profiled by TAPS methylation sequencing;
the DNMT3A case was also profiled by RNA-seq. The two are integrated to ask which
programs change as the clone expands and at transformation, and which changes are
attributable to DNMT3A and which to TET2.

No data are included. Sample details are kept in local configuration files that are
never committed (see [config/README.md](config/README.md)).

## Pipeline

```
TAPS FASTQ ── 01_msuite2.sh ──► alignment + per-CpG C/T counts (Msuite2, -m TAPS)
                │
                ├─ 02_taps_inputs.py ──► methylation per CpG (T / (C + T))
                │       └─ 03_dmr_calling.sh ──► DMRs (metilene, primary; Defiant)
                │               ├─ 04_homer_known_motifs.sh ──► motifs in DMRs lost / gained
                │               └─ 05_homer_annotation.sh   ──► HOMER genomic annotation
                ├─ sample_identity_check.py ──► genotype concordance between libraries
                │       └─ allele_balance_counts.py ──► allele balance at the same SNPs (LOH scan)
                ├─ cpg_depth_summary.py ──► CpGs with >= 5 reads and their mean methylation, per library
                ├─ 06_vaf_pileup.sh ──► strand-aware base counts at driver genes (VAFs)
                └─ 05_variants_server/02 ──► methylation-masked BAMs (CpG bases that report methylation → N)
                        ├─ 03_somatic_calls.sh ──► genome-wide SNVs and indels (Mutect2, one run per case)
                        ├─ 04_cnv_loh.sh ──► copy number and allelic balance (GATK somatic CNV)
                        ├─ 05_annotate_collect.sh ──► mappability flags, Funcotator, transfer folder
                        └─ 06_site_counts.sh ──► orientation-split read counts at the clone mutations

RNA BAM ── 01_featurecounts.sh ──► gene and 5'UTR counts ──► 03_rnaseq/02-06 (edgeR, NOISeq,
       │                                                    promoter usage, ORA, GSEA)
       ├─ 02_rmats.sh ──► 03_rnaseq/07_rmats_downstream.R
       └─ 03_shiba.sh ──► 03_rnaseq/08_shiba_downstream.R

methylation + RNA-seq ── 04_integration/01-17 ──► region methylation, direction of change by
       genomic context, promoter switches against expression, TF activity, DMR annotation,
       gene panels, Hi-C loops, locus plots, driver mutations and VAFs from the TAPS reads,
       allelic balance along the genome, recipient DNA (chimerism), genome-wide copy
       number and point mutations (and whether the clone was already in earlier samples)
```

## Layout

| Folder | Runs on | Contents |
|---|---|---|
| `config/` | both | path and sample templates; copy each `*.example.*` to the real name and fill it in |
| `01_methylation_server/` | server | Msuite2, TAPS inputs, metilene / Defiant, HOMER, VAF pileups, tool versions |
| `02_methylation_local/` | analysis computer | sample identity check; allele balance at the identity SNPs; CpG depth summary |
| `03_rnaseq/server/` | server | featureCounts, rMATS, Shiba |
| `03_rnaseq/` | analysis computer | RNA-seq analysis (`00_config.R`, `01`-`08`, `dge_functions.R`) |
| `04_integration/` | analysis computer | methylation and RNA-seq integration (`01`-`11`); mutations and VAFs (`13`); allelic balance (`14`); chimerism (`15`); genome-wide copy number (`16`) and point mutations (`17`) |
| `05_variants_server/` | server | genome-wide variants from the TAPS reads: masked BAMs, Mutect2, GATK copy number, annotation, site counts (`run_variants.sh` runs 02-05) |
| `resources/` | analysis computer | pinned CollecTRI network; download notes for the other annotation |
| `docs/METHODS.md` | | methods |
| `environment/versions.md` | | software versions |

## Running

Every script is run from the repository root.

1. **Configure.** Copy `config/paths.example.sh` to `config/paths.sh` (server) and
   `config/paths.example.R` to `config/paths.R` (analysis computer); copy
   `config/samples.example.tsv`, `config/comparisons.example.tsv` and
   `config/rnaseq_contrasts.example.tsv` to `samples.tsv`, `comparisons.tsv` and
   `rnaseq_contrasts.tsv` and fill them in.
2. **Methylation, on the server.**
   ```
   bash 01_methylation_server/01_msuite2.sh                 # writes the Msuite2 batch script; run it
   bash 01_methylation_server/03_dmr_calling.sh             # TAPS inputs, metilene, Defiant
   bash 01_methylation_server/04_homer_known_motifs.sh
   bash 01_methylation_server/06_vaf_pileup.sh
   ```
   Copy `<sample>_msuite2/Msuite2.CpG.meth.call.xz` and the DMR output folder (with
   `vaf/` from `06_vaf_pileup.sh`) to the analysis computer (`emseq/` under `DATA_ROOT`).
   `bash 01_methylation_server/00_tool_versions.sh` records the server's tool versions.
3. **RNA-seq, on the server:** `03_rnaseq/server/01_featurecounts.sh`, `02_rmats.sh`,
   `03_shiba.sh`; copy the counts and the rMATS / Shiba folders to `rnaseq/`.
4. **Genome-wide variants, on the server** (after Msuite2; bash >= 4.4, samtools, Docker):
   ```
   bash 05_variants_server/01_inspect_bams.sh                # read-only check of the alignments
   bash 05_variants_server/run_variants.sh --test chr21      # rehearsal on one contig
   nohup bash 05_variants_server/run_variants.sh > variants_run.log 2>&1 &
   ```
   Copy `<VAR_OUT>/transfer/` to `emseq/04_variants/`. After the first run of
   `04_integration/17_genomewide_snvs.R`, copy its `site_counts_sites.tsv` to the
   server, run `bash 05_variants_server/06_site_counts.sh <that file>`, copy
   `<VAR_OUT>/site_counts/` to `emseq/04_variants/site_counts/` and run 17 again.
5. **Analysis computer.**
   ```
   python3 02_methylation_local/sample_identity_check.py <DATA_ROOT>/emseq <DATA_ROOT>/emseq/00_sample_identity
   python3 02_methylation_local/allele_balance_counts.py <DATA_ROOT>/emseq <DATA_ROOT>/emseq/00_sample_identity
   python3 02_methylation_local/cpg_depth_summary.py <DATA_ROOT>/emseq <DATA_ROOT>/emseq/00_sample_identity
   bash run_rnaseq.sh
   bash run_integration.sh
   ```
6. **HOMER annotation, on the server:** copy `integration/dmr/bed/` (written by
   `04_integration/06_dmr_master.R`) to the server and run
   `01_methylation_server/05_homer_annotation.sh`.

## Conventions

- Comparisons are named `first_vs_second`; every difference is first-named minus
  second-named, in methylation and RNA-seq alike.
- Methylation is TAPS: a methylated C is read as T. Methylation is the share of T among
  C and T calls at a CpG, both strands combined.
- With one library per timepoint, p-values and FDRs rank features; confidence comes from
  agreement between independent comparisons, a reference comparison (the recipient at
  1 month against the donor 5 years after donation) and recurrence in the second case. See [docs/METHODS.md](docs/METHODS.md).
