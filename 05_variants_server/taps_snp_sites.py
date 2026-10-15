#!/usr/bin/env python3
"""
taps_snp_sites.py <in.vcf.gz> <ref.fa> <out.vcf.gz>                2026-10-13

SNP sites for allelic counts on TAPS libraries (used by 04_cnv_loh.sh).
Drops biallelic SNPs whose alternative allele creates a CpG that the reference
lacks: an alternative C before a reference G, or an alternative G after a
reference C. On that allele the new CpG is methylated like any other, and TAPS
reads its methylated C as T on one strand, so the alternative allele would be
under-counted and the allele fraction biased. Reference CpG sites are kept:
02_prepare_bams.sh masks their methylation-reporting bases (base quality 0)
and the other strand reads the genotype.
"""
import sys

import pysam


def main():
    vin_path, ref_path, out_path = sys.argv[1:4]
    fa = pysam.FastaFile(ref_path)
    vin = pysam.VariantFile(vin_path)
    vout = pysam.VariantFile(out_path, "wz", header=vin.header)
    kept = dropped = 0
    chrom = seq = None
    for rec in vin:
        if rec.chrom != chrom:
            chrom = rec.chrom
            seq = fa.fetch(chrom).upper()
        p = rec.pos - 1
        alt = rec.alts[0].upper() if rec.alts else ""
        prev = seq[p - 1] if p > 0 else "N"
        nxt = seq[p + 1] if p + 1 < len(seq) else "N"
        if (alt == "C" and nxt == "G") or (alt == "G" and prev == "C"):
            dropped += 1
            continue
        vout.write(rec)
        kept += 1
    vout.close()
    print(f"{kept} sites kept, {dropped} dropped (alternative allele creates a CpG)")


if __name__ == "__main__":
    main()
