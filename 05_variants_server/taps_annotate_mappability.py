#!/usr/bin/env python3
"""
taps_annotate_mappability.py <in.vcf.gz> <umap.bed.gz> <segdup.bed.gz> <out.vcf.gz>     2026-10-14

Marks each call (used by 05_annotate_collect.sh) with two INFO flags from
the tracks the copy-number step already uses:
  UMAP100  the position lies in a region where every 100-base read is unique (Umap k100)
  SEGDUP   the position lies in a segmental duplication
Calls outside UMAP100 or inside SEGDUP sit where reads from similar copies can map
wrongly; in the rehearsals these regions held most of the calls whose alternative
reads mapped poorly. The flags are only added; nothing is removed here.
The BED files are read once, streamed (they need not be indexed or sorted the
same way as the VCF); coordinates are 0-based half-open, VCF positions 1-based.
"""
import bisect
import gzip
import sys

import pysam


def marked(bed_path, pos_by_contig):
    hit = set()
    with gzip.open(bed_path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t", 3)
            ps = pos_by_contig.get(f[0])
            if not ps:
                continue
            s, e = int(f[1]), int(f[2])
            i = bisect.bisect_left(ps, s + 1)          # first 1-based position > s
            while i < len(ps) and ps[i] <= e:          # 1-based position <= e
                hit.add((f[0], ps[i]))
                i += 1
    return hit


def main():
    vin_path, umap, segdup, out_path = sys.argv[1:5]
    vin = pysam.VariantFile(vin_path)
    pos = {}
    for rec in vin:
        pos.setdefault(rec.chrom, set()).add(rec.pos)
    pos = {c: sorted(p) for c, p in pos.items()}
    in_umap, in_segdup = marked(umap, pos), marked(segdup, pos)

    hdr = vin.header.copy()
    hdr.add_line('##INFO=<ID=UMAP100,Number=0,Type=Flag,Description="Position in a region where every 100-base read is unique (Umap k100)">')
    hdr.add_line('##INFO=<ID=SEGDUP,Number=0,Type=Flag,Description="Position in a segmental duplication">')
    vin = pysam.VariantFile(vin_path)
    n = u = d = 0
    with pysam.VariantFile(out_path, "wz", header=hdr) as vout:
        for rec in vin:
            rec.translate(hdr)
            key = (rec.chrom, rec.pos)
            if key in in_umap:
                rec.info["UMAP100"] = True
                u += 1
            if key in in_segdup:
                rec.info["SEGDUP"] = True
                d += 1
            vout.write(rec)
            n += 1
    pysam.tabix_index(out_path, preset="vcf", force=True)
    print(f"{n} calls: {u} in uniquely mappable regions (UMAP100), {d} in segmental duplications (SEGDUP)")


if __name__ == "__main__":
    main()
