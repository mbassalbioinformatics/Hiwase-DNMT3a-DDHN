# Software

Versions used for the analysis.

## Analysis computer

- R 4.6.0, Bioconductor 3.23
- Python 3.10.14

| R package | version |
|---|---|
| data.table | 1.18.6.1 |
| ggplot2 | 4.0.3 |
| ggrepel | 0.9.8 |
| patchwork | 1.3.2 |
| ComplexHeatmap | 2.28.0 |
| circlize | 0.4.18 |
| decoupleR | 2.17.0 |
| msigdbr | 26.1.1 |
| jsonlite | 2.0.0 |
| edgeR | 4.10.5 |
| limma | 3.68.5 |
| NOISeq | 2.56.0 |
| GSVA | 2.6.6 |
| clusterProfiler | 4.20.0 |
| enrichplot | 1.32.1 |
| DOSE | 4.6.0 |
| enrichit | 0.2.6 |
| org.Hs.eg.db | 3.23.1 |
| BiocParallel | 1.46.0 |
| pheatmap | 1.0.13 |
| RColorBrewer | 1.1.3 |
| viridis | 0.6.5 |
| Glimma | 2.22.1 |
| stringr | 1.6.0 |

## Server

Recorded with `bash 01_methylation_server/00_tool_versions.sh` (PATH and every conda environment).

| Tool | Version | Note |
|---|---|---|
| Msuite2 | 2.2.0 | |
| bowtie2 | 2.5.4 | used by Msuite2 |
| samtools | 1.21 | VAF pileups |
| bedtools | 2.31.1 | metilene inputs |
| metilene | not printed by the binary | |
| Defiant | not printed by the binary | |
| HOMER | 4.11.1 or 5.1 | two installs; the one used is to be confirmed |
| featureCounts (Subread) | 2.1.1 | |
| rMATS-turbo | 4.3.0 | |
| Shiba | 0.8.2 | |
| STAR | 2.7.10b | installed; RNA-seq aligner to be confirmed from the BAM header |
| Python | 3.9.18 | |
| R | 4.6.1 | |
