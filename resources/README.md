# Annotation resources

Kept in the repository:

- `collectri_human_2026-10-08.tsv`: the CollecTRI TF-target network used for TF
  activity (`04_integration/02_tf_activity.R`), downloaded from OmniPath on
  8 October 2026 and processed as `decoupleR::get_collectri(split_complexes = FALSE)`
  does (JUN/FOS complexes to AP1, REL/NFKB to NFKB; mode of regulation +1 for
  stimulation, -1 otherwise). Source query:
  `https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&organisms=9606&loops=yes&fields=type,curation_effort&license=academic`.
  To rebuild it from a fresh download, delete the file and run
  `COLLECTRI_RAW=<downloaded.tsv> Rscript 04_integration/02_tf_activity.R`.

Downloaded into `RES` (see `config/paths.R`):

- `cpgIslandExt.txt.gz`: UCSC hg38 CpG islands,
  `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz`
- `encode_cCREs_hg38_ucsc.bed`: ENCODE candidate cis-regulatory elements, from the UCSC
  track `encodeCcreCombined.bb`
  (`https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb`),
  converted to BED, e.g. with pyBigWig:
  ```python
  import pyBigWig
  bb = pyBigWig.open("encodeCcreCombined.bb")
  with open("encode_cCREs_hg38_ucsc.bed", "w") as out:
      for chrom in bb.chroms():
          for s, e, rest in bb.entries(chrom, 0, bb.chroms(chrom)) or []:
              out.write(f"{chrom}\t{s}\t{e}\t{rest}\n")
  ```
- `gene_panels.tsv`, `gene_panels_provenance.tsv`: written by `04_integration/07_gene_panels.R`.

Elsewhere (paths in `config/paths.R`):

- MANE Select v1.0 GTF (`MANE.GRCh38.v1.0.ensembl_genomic.gtf.gz`, NCBI/EMBL-EBI).
- Hi-C loop calls of sorted human HSC, CMP, GMP and MEP (two replicates each), from the
  Pimanda laboratory (`*_Pimanda_*.bedpe.xz`).
- HOMER with the hg38 genome installed (server).
- Ensembl 107 GTF for RNA-seq quantification (server).
