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

| Tool | Version |
|---|---|
| Msuite2 | 2.2.0 |
| featureCounts (Subread) | 2.1.1 |
| rMATS-turbo | 4.3.0 |
| Shiba | to be confirmed (configuration written for 0.8.2) |
| metilene | to be recorded |
| Defiant | to be recorded |
| HOMER | to be recorded |
| bedtools | to be recorded |
| samtools | to be recorded |

Record the server versions with `bash 01_methylation_server/00_tool_versions.sh`, which looks for
each tool on PATH and in every conda environment.
