#!/usr/bin/env python3
###############################################################################
# allele_balance_counts.py
#
# Strand-resolved allele counts at the CpG-destroying SNP sites found by
# sample_identity_check.py (union_sites.txt), for an allelic-imbalance (loss
# of heterozygosity) scan along the genome.
#
# At a reference CpG, the Watson reads report the C (C/T = methylation, "other"
# = a SNP at the C) and the Crick reads report the G of the CpG (C/T on the
# Crick strand = methylation, "other" = a SNP at the G). Methylation never makes
# an "other" base, so other / (C + T + other) on each strand is the allele
# fraction of a transversion SNP at that base (B-allele fraction, BAF).
#
# USAGE  python3 02_methylation_local/allele_balance_counts.py <emseq_dir> <identity_dir> [sample_sheet]
#        (<identity_dir> is the output folder of sample_identity_check.py; the libraries
#        are those marked analysed or extra in the sample sheet, default config/samples.tsv)
# Output in <identity_dir>: <sample>.strand_counts.tsv
#        chr  pos  w_ct  w_other  c_ct  c_other   (pos = CpG C, 1-based; the
#        Crick counts belong to pos + 1)
###############################################################################
import csv, os, subprocess, sys
from multiprocessing import Pool


def run(job):
    s, emseq, out, union = job
    p = subprocess.Popen(["xz", "-dc", os.path.join(emseq, f"{s}_msuite2", "Msuite2.CpG.meth.call.xz")],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    n = 0
    with open(os.path.join(out, f"{s}.strand_counts.tsv"), "w") as o:
        o.write("chr\tpos\tw_ct\tw_other\tc_ct\tc_other\n")
        for line in p.stdout:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if f"{c[0]}:{c[1]}" in union:
                o.write(f"{c[0]}\t{c[1]}\t{int(c[3]) + int(c[4])}\t{c[5]}\t{int(c[7]) + int(c[8])}\t{c[9]}\n")
                n += 1
    return s, n


if __name__ == "__main__":
    emseq, out = sys.argv[1], sys.argv[2]
    sheet = sys.argv[3] if len(sys.argv) > 3 else "config/samples.tsv"
    SAMPLES = [r["sample"] for r in csv.DictReader(open(sheet), delimiter="\t") if r["methylation"] in ("analysed", "extra")]
    union = set(l.strip() for l in open(os.path.join(out, "union_sites.txt")) if l.strip())
    with Pool(len(SAMPLES)) as p:
        for s, n in p.map(run, [(s, emseq, out, union) for s in SAMPLES]):
            print(s, n, "sites")
