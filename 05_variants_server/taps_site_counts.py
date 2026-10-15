#!/usr/bin/env python3
"""
taps_site_counts.py <sites.tsv> <out.tsv> <name=bam> [<name=bam> ...]          2026-10-15

Reads at given SNV sites in every library, split by read orientation, as
Mutect2 counts them (used by 06_site_counts.sh):
  F1R2 = read 1 forward / read 2 reverse = original top strand (OT)
  F2R1 = read 1 reverse / read 2 forward = original bottom strand (OB)
TAPS converts a C to T only on the strand that carries it, so a C>T is judged
on F2R1 reads and a G>A on F1R2 reads; other changes on both.
Bases need base quality >= 20 and the read mapping quality >= 10; duplicates,
secondary and QC-failed reads are skipped; overlapping mates are counted once.

sites.tsv: tab-separated with a header: set, chr, pos, ref, alt (1-based position)
out.tsv:   set, chr, pos, ref, alt, sample, f1_ref, f1_alt, f1_other, f2_ref, f2_alt, f2_other
"""
import sys

import pysam

MIN_BQ, MIN_MQ = 20, 10


def main():
    sites_path, out_path = sys.argv[1:3]
    bams = dict(a.split("=", 1) for a in sys.argv[3:])
    with open(sites_path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        sites = [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]
    with open(out_path + ".tmp", "w") as out:
        out.write("\t".join(["set", "chr", "pos", "ref", "alt", "sample", "f1_ref", "f1_alt", "f1_other",
                             "f2_ref", "f2_alt", "f2_other"]) + "\n")
        for name, path in bams.items():
            bam = pysam.AlignmentFile(path)
            for s in sites:
                p0 = int(s["pos"]) - 1
                c = {k: 0 for k in ("f1_ref", "f1_alt", "f1_other", "f2_ref", "f2_alt", "f2_other")}
                for col in bam.pileup(s["chr"], p0, p0 + 1, truncate=True, min_base_quality=MIN_BQ,
                                      min_mapping_quality=MIN_MQ, ignore_overlaps=True, stepper="all",
                                      flag_filter=0x4 | 0x100 | 0x200 | 0x400 | 0x800):
                    for pr in col.pileups:
                        if pr.is_del or pr.is_refskip or pr.query_position is None:
                            continue
                        r = pr.alignment
                        # filters applied here: the pileup's own base-quality option does not act with stepper="all"
                        if r.mapping_quality < MIN_MQ or r.flag & (0x4 | 0x100 | 0x200 | 0x400 | 0x800):
                            continue
                        if r.query_qualities[pr.query_position] < MIN_BQ:
                            continue
                        b = r.query_sequence[pr.query_position]
                        if b == "N":
                            continue
                        f1 = (r.is_read1 and not r.is_reverse) or (r.is_read2 and r.is_reverse)
                        k = ("f1_" if f1 else "f2_") + ("ref" if b == s["ref"] else "alt" if b == s["alt"] else "other")
                        c[k] += 1
                out.write("\t".join([s["set"], s["chr"], s["pos"], s["ref"], s["alt"], name] + [str(c[k]) for k in
                                     ("f1_ref", "f1_alt", "f1_other", "f2_ref", "f2_alt", "f2_other")]) + "\n")
            print(f"{name}: {len(sites)} sites", flush=True)
    import os
    os.replace(out_path + ".tmp", out_path)


if __name__ == "__main__":
    main()
