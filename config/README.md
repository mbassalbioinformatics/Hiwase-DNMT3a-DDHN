# Configuration

Each file here has a committed template (`*.example.*`) and a local copy with the
real values. The local copies hold sample details and are listed in `.gitignore`;
never commit them.

| Local file | Template | Used by |
|---|---|---|
| `paths.sh` | `paths.example.sh` | server scripts (`01_methylation_server/`, `03_rnaseq/server/`) |
| `paths.R` | `paths.example.R` | every R script |
| `samples.tsv` | `samples.example.tsv` | everything that needs the sample list |
| `comparisons.tsv` | `comparisons.example.tsv` | DMR calling and the integration |
| `rnaseq_contrasts.tsv` | `rnaseq_contrasts.example.tsv` | RNA-seq analysis, rMATS |

## samples.tsv

One row per sequenced library, tab-separated, every column read as text.

| Column | Meaning |
|---|---|
| `sample` | library ID, as used in file names (`<sample>_msuite2/`, `<sample>_qual_trimmed.sorted.bam`) |
| `key` | short role key used in code: `D18` donor at donation, `D30` day 30, `Y35` 3.5 years, `DX` DDHN diagnosis, `D23` donor later sample, `T2A` / `T2DX` second case, `X` extra, `NEG` negative control |
| `case`, `person` | case (gene) and genotype group from the identity check |
| `uln` | sequencing library number, found in the FASTQ file names |
| `methylation` | `analysed` (in every analysis), `extra` (identity check and region methylation only), `control` (Msuite2 only) |
| `rna_seq` | `yes` if the library has RNA-seq |
| `plot_label`, `short_label` | labels for locus plots and tables |
| `role`, `status`, `note` | description; `final` or `provisional` |
| `sample_date`, `hsct_date` | ISO dates |
| `material` | e.g. BMNC, PBMNC |
| `vaf_gene`, `vaf_pct`, `vaf_note` | clone VAF from clinical testing, or why there is none |
| `subject`, `tp_label`, `short`, `disease`, `clone`, `clone_ann`, `clone_colour`, `vaf_source`, `rna_label`, `rna_colour` | RNA-seq sample annotation and colours |

## comparisons.tsv

`id` (`g1_vs_g2`), `g1`, `g2`, `case`, `status`, `label`, `note`, then:

| Column | Meaning |
|---|---|
| `role` | what the comparison measures (e.g. `transformation`, `clonal_establishment`, `reference`) |
| `large_vs_small` | `yes` when a large clone is compared with a small one |
| `rna_contrast` | the matching RNA-seq contrast id, if any |
| `short` | short name used in the overlap table of `04_integration/04_global_direction.R` |
| `overlap` | `0` for the query comparison of that overlap test, `1`, `2`, ... for the comparisons it is tested against, in order |
| `ceiling` | `1` and `2` mark two comparisons sharing no sample, tested against each other as the ceiling for a shared program |

## rnaseq_contrasts.tsv

`id`, `num`, `den`, `label`, `class` (`primary`, `context`, `ruler`), `meaning`.
Direction is always numerator minus denominator.
