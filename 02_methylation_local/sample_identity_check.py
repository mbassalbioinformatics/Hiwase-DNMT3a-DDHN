#!/usr/bin/env python3
###############################################################################
# sample_identity_check.py
#
# Which libraries come from the same person? Genotype fingerprint from
# CpG-destroying SNPs in Msuite2's per-CpG calls (Msuite2.CpG.meth.call.xz).
#
# At a reference CpG the reads show C or T (methylation, under any conversion
# chemistry) unless the person carries a SNP that destroys the CpG; then the
# reads show another base (column wOther / cOther). The fraction of "other" bases
# is therefore a genotype readout that methylation does not touch:
#   ~0 = reference, ~0.5 = heterozygous, ~1 = homozygous.
#
#   pass 1   per sample: candidate sites (>= 10 reads, other fraction >= 0.25),
#            plus coverage, mean CpG methylation (T / (C + T), TAPS) and the share
#            of chrX CpGs at intermediate methylation (X inactivation in females)
#   union    all candidate sites from every sample
#   pass 2   per sample: reads and other-base reads at every union site
#   compare  genotype per site (0: < 0.10, 1: 0.25-0.75, 2: > 0.90, else none);
#            for each pair, the share of sites with the same genotype among sites
#            with >= 10 reads and a genotype in both. Same person ~0.99; unrelated
#            people ~0.6.
#
# USAGE  python3 02_methylation_local/sample_identity_check.py <emseq_dir> <out_dir> [sample_sheet]
#        (<emseq_dir> holds <sample>_msuite2/Msuite2.CpG.meth.call.xz; the libraries
#        are those marked analysed or extra in the sample sheet, default
#        config/samples.tsv)
# Output in <out_dir>: <sample>.json (pass 1 summary), union_sites.txt,
#        <sample>.union_counts.tsv, identity_concordance.tsv (pairwise)
###############################################################################
import collections, csv, itertools, json, os, subprocess, sys
from multiprocessing import Pool

MIN_READS, CAND_FRAC = 10, 0.25


def stream(path):
    p = subprocess.Popen(["xz", "-dc", path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        if line[0] != "#":
            yield line.rstrip("\n").split("\t")


def pass1(job):
    s, emseq, out = job
    cov = collections.Counter(); n = 0; tot = 0.0; xm = [0, 0]; cand = []
    for c in stream(os.path.join(emseq, f"{s}_msuite2", "Msuite2.CpG.meth.call.xz")):
        wC, wT, wO, cC, cT, cO = (int(c[i]) for i in (3, 4, 5, 7, 8, 9))
        allr = wC + wT + wO + cC + cT + cO
        cov[c[0]] += allr
        ct = wC + wT + cC + cT
        if ct >= MIN_READS:
            m = (wT + cT) / ct; n += 1; tot += m          # TAPS: T = methylated
            if c[0] == "chrX":
                xm[1] += 1; xm[0] += 0.3 <= m <= 0.7
        if allr >= MIN_READS and (wO + cO) / allr >= CAND_FRAC:
            cand.append(f"{c[0]}:{c[1]}")
    json.dump({"sample": s, "coverage_by_chr": cov, "cpgs_10x": n, "mean_methylation": tot / n,
               "chrX_intermediate": xm, "candidates": len(cand)}, open(os.path.join(out, f"{s}.json"), "w"))
    return cand


def pass2(job):
    s, emseq, out, union = job
    with open(os.path.join(out, f"{s}.union_counts.tsv"), "w") as o:
        for c in stream(os.path.join(emseq, f"{s}_msuite2", "Msuite2.CpG.meth.call.xz")):
            k = f"{c[0]}:{c[1]}"
            if k in union:
                allr = sum(int(c[i]) for i in (3, 4, 5, 7, 8, 9)); oth = int(c[5]) + int(c[9])
                o.write(f"{k}\t{allr}\t{oth}\n")
    return s


def genotype(allr, oth):
    if allr < MIN_READS:
        return None
    f = oth / allr
    return 0 if f < 0.10 else 1 if 0.25 <= f <= 0.75 else 2 if f > 0.90 else None


if __name__ == "__main__":
    emseq, out = sys.argv[1], sys.argv[2]
    sheet = sys.argv[3] if len(sys.argv) > 3 else "config/samples.tsv"
    SAMPLES = [r["sample"] for r in csv.DictReader(open(sheet), delimiter="\t") if r["methylation"] in ("analysed", "extra")]
    os.makedirs(out, exist_ok=True)
    with Pool(len(SAMPLES)) as p:
        cands = p.map(pass1, [(s, emseq, out) for s in SAMPLES])
    union = set().union(*cands)
    open(os.path.join(out, "union_sites.txt"), "w").write("\n".join(sorted(union)) + "\n")
    with Pool(len(SAMPLES)) as p:
        p.map(pass2, [(s, emseq, out, union) for s in SAMPLES])
    gt = {}
    for s in SAMPLES:
        g = {}
        for line in open(os.path.join(out, f"{s}.union_counts.tsv")):
            k, a, o = line.rstrip("\n").split("\t")
            v = genotype(int(a), int(o))
            if v is not None:
                g[k] = v
        gt[s] = g
    with open(os.path.join(out, "identity_concordance.tsv"), "w") as o:
        o.write("sample_a\tsample_b\tsites_compared\tconcordance\n")
        for a, b in itertools.combinations(SAMPLES, 2):
            shared = gt[a].keys() & gt[b].keys()
            same = sum(gt[a][k] == gt[b][k] for k in shared)
            o.write(f"{a}\t{b}\t{len(shared)}\t{same / len(shared) if shared else float('nan'):.4f}\n")
    print(open(os.path.join(out, "identity_concordance.tsv")).read())
