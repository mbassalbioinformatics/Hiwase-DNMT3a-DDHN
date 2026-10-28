## config/settings.R -- analysis settings that several scripts share; change them here only.
## Sourced at the end of config/paths.R, so every script that sources paths.R has them.

## Genome tiles: consecutive, non-overlapping windows of TILE_BP bp on chr1-22 and chrX, starting at
## position 1. A tile's methylation is the mean of its CpGs with >= 5 reads, if it has >= 3 of them
## (01_region_methylation.R). 500 bp rather than 1 kb: methylation changes are focal (metilene DMRs
## have a median length of 250-490 bp per comparison), so 1-kb tiles dilute them; 500-bp tiles find
## about three times more change beyond noise and pick up more of the DMRs, and still hold 90% of the
## CpGs (250-bp tiles: 73%). The evidence: 04_integration/21_tile_size_check.R.
TILE_BP   <- 500L
TILE_FILE <- file.path(INT, "region_methylation", sprintf("tiles_%dbp_methylation.tsv.gz", TILE_BP))
TILE_TXT  <- sprintf("%d-bp", TILE_BP)            # for figure text, e.g. "500-bp tiles"

## A change: a tile or a promoter gains or loses methylation when it moves by at least CUT (0.20 = 20
## points), up or down, at any level. Stricter definitions appear only as sensitivity checks.
CUT     <- 0.20
CUT_TXT <- sprintf("%.0f", 100 * CUT)             # "20", for "20 points"
