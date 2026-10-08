#!/usr/bin/env python3
###############################################################################
# vaf_pileup_counts.py                                              2026-10-08
#
# Reads `samtools mpileup` output (no reference given, so read bases are printed
# as letters: upper case = read aligned forward, lower case = reverse) for read 1
# and read 2 separately, and writes base counts per original DNA strand.
#
# Directional TAPS libraries:
#   OT (original top)    = read 1 forward + read 2 reverse
#   OB (original bottom) = read 1 reverse + read 2 forward
# TAPS turns a methylated C into T on the strand that carries it, so at a CpG:
#   top-strand C:    the OT reads show C/T by methylation; OB reads show G/A by genotype only
#   bottom-strand C: the OB reads show G/A (top-strand view) by methylation; OT reads show
#                    the genotype only
# The caller (analysis done on the Mac) picks the strand that reports genotype.
#
# USAGE  vaf_pileup_counts.py <sample> <read1.mpileup> <read2.mpileup>  > counts.tsv
# Output: sample chr pos OT_A OT_C OT_G OT_T OB_A OB_C OB_G OB_T
###############################################################################
import sys
from collections import defaultdict

BASES = "ACGT"


def parse_bases(s):
    """Counts of forward (upper) and reverse (lower) A/C/G/T in an mpileup base string."""
    fw = dict.fromkeys(BASES, 0); rv = dict.fromkeys(BASES, 0)
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == "^":                       # read start: next char is mapping quality
            i += 2; continue
        if c in "+-":                      # indel: +3ACG or -2TT
            j = i + 1
            while j < n and s[j].isdigit():
                j += 1
            i = j + int(s[i + 1:j]); continue
        if c in "ACGT":
            fw[c] += 1
        elif c in "acgt":
            rv[c.upper()] += 1
        i += 1                             # '$', '*', '>', '<', 'N', 'n' are skipped
    return fw, rv


def read_pileup(path, table, mate):
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            fw, rv = parse_bases(p[4])
            key = (p[0], int(p[1]))
            ot, ob = table[key]
            if mate == 1:                  # read 1: forward = OT, reverse = OB
                for b in BASES: ot[b] += fw[b]; ob[b] += rv[b]
            else:                          # read 2: forward = OB, reverse = OT
                for b in BASES: ob[b] += fw[b]; ot[b] += rv[b]


if __name__ == "__main__":
    sample, r1, r2 = sys.argv[1:4]
    table = defaultdict(lambda: (dict.fromkeys(BASES, 0), dict.fromkeys(BASES, 0)))
    read_pileup(r1, table, 1)
    read_pileup(r2, table, 2)
    out = sys.stdout
    for (chrom, pos) in sorted(table, key=lambda k: (k[0], k[1])):
        ot, ob = table[(chrom, pos)]
        if sum(ot.values()) + sum(ob.values()) == 0:
            continue
        out.write("\t".join([sample, chrom, str(pos)] + [str(ot[b]) for b in BASES] + [str(ob[b]) for b in BASES]) + "\n")
