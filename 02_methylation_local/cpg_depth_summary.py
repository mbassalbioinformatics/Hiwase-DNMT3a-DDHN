#!/usr/bin/env python3
###############################################################################
# cpg_depth_summary.py
#
# Per library: CpGs on chr1-22 and chrX with >= 5 informative reads (C + T on
# both strands, the depth every methylation analysis requires of a CpG), and
# the mean methylation of those CpGs (TAPS reading: T = methylated). Feeds the
# library table of the reports.
#
# USAGE  python3 02_methylation_local/cpg_depth_summary.py <emseq_dir> <identity_dir> [sample_sheet]
#        (the libraries are those marked analysed in the sample sheet,
#        default config/samples.tsv)
# Reads  <emseq_dir>/<sample>_msuite2/Msuite2.CpG.meth.call.xz
# Output <identity_dir>/cpg_depth_summary.tsv   sample  cpgs_min5  mean_methylation_min5
###############################################################################
import csv, os, subprocess, sys
from multiprocessing import Pool

MIN_READS = 5
KEEP = {f"chr{c}" for c in list(range(1, 23)) + ["X"]}


def summarise(job):
    s, emseq = job
    p = subprocess.Popen(["xz", "-dc", os.path.join(emseq, f"{s}_msuite2", "Msuite2.CpG.meth.call.xz")],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    n = 0; tot = 0.0
    for line in p.stdout:
        if line[0] == "#":
            continue
        c = line.split("\t", 10)
        if c[0] not in KEEP:
            continue
        wC, wT, cC, cT = int(c[3]), int(c[4]), int(c[7]), int(c[8])
        ct = wC + wT + cC + cT
        if ct >= MIN_READS:
            n += 1; tot += (wT + cT) / ct
    if p.wait() != 0:
        sys.exit(f"{s}: xz failed")
    return s, n, tot / n if n else float("nan")


if __name__ == "__main__":
    emseq, out = sys.argv[1], sys.argv[2]
    sheet = sys.argv[3] if len(sys.argv) > 3 else "config/samples.tsv"
    SAMPLES = [r["sample"] for r in csv.DictReader(open(sheet), delimiter="\t") if r["methylation"] == "analysed"]
    with Pool(len(SAMPLES)) as pool:
        res = pool.map(summarise, [(s, emseq) for s in SAMPLES])
    with open(os.path.join(out, "cpg_depth_summary.tsv"), "w") as f:
        f.write(f"sample\tcpgs_min{MIN_READS}\tmean_methylation_min{MIN_READS}\n")
        for s, n, m in res:
            f.write(f"{s}\t{n}\t{m:.4f}\n")
            print(s, n, f"{m:.4f}")
