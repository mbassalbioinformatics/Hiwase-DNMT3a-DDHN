#!/usr/bin/env python3
###############################################################################
# 02_taps_inputs.py                                                    (server)
#
# Inputs for DMR calling from Msuite2's per-CpG base counts
# (Msuite2.CpG.meth.call.xz: #chr Locus Total wC wT wOther Context cC cT cOther).
# In TAPS a methylated C is read as T, so per CpG, both strands combined:
#     methylated   = wT + cT
#     unmethylated = wC + cC
#     methylation  = methylated / (methylated + unmethylated)
# CpGs without a C or T call are skipped. Only chr1-22 and chrX are written
# (all libraries are female; this drops chrY, chrM, the spike-ins and the
# alt/random/Un contigs). Start = Locus, end = Locus + 1.
#
# Writes, per sample, xz-compressed and without a header:
#   <sample>.CpG.meth.TAPS.DefiantFormat3.bedgraph.xz  chr start end fraction   (metilene, Defiant)
#   <sample>.CpG.meth.TAPS.MDFormat.bedgraph.xz        chr start end % methylated unmethylated
#
# --check also compares every CpG with Msuite2's own bedgraph
# (Msuite2.CpG.meth.bedgraph.xz, written by Msuite2 run with -m TAPS): same CpG,
# same value. A sample that fails is not written.
#
# USAGE
#   python3 02_taps_inputs.py <msuite2_dir> <sample> [...] [--out DIR] [--check]
###############################################################################
import os, subprocess, sys
from multiprocessing import Pool

KEEP = {f"chr{c}" for c in list(range(1, 23)) + ["X"]}


def num(x, places):
    """Fixed notation without trailing zeros."""
    return f"{x:.{places}f}".rstrip("0").rstrip(".")


def reformat(job):
    ms_dir, sample, out_dir, check = job
    d = os.path.join(ms_dir, f"{sample}_msuite2")
    call = os.path.join(d, "Msuite2.CpG.meth.call.xz")
    own = os.path.join(d, "Msuite2.CpG.meth.bedgraph.xz")
    if not os.path.exists(call) or (check and not os.path.exists(own)):
        return f"{sample}\tFAILED\tmissing {'call file' if not os.path.exists(call) else 'Msuite2 bedgraph'} in {d}"
    od = out_dir or d
    outs = [os.path.join(od, f"{sample}.CpG.meth.TAPS.{k}.bedgraph.xz") for k in ("MDFormat", "DefiantFormat3")]
    tmps = [o + ".tmp" for o in outs]

    pc = subprocess.Popen(["xz", "-dc", call], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    pb = subprocess.Popen(["xz", "-dc", own], stdout=subprocess.PIPE, text=True, bufsize=1 << 20) if check else None
    fmd, fdf = (open(t, "wb") for t in tmps)
    xmd = subprocess.Popen(["xz", "-6", "-T2", "-c"], stdin=subprocess.PIPE, stdout=fmd, text=True)
    xdf = subprocess.Popen(["xz", "-6", "-T2", "-c"], stdin=subprocess.PIPE, stdout=fdf, text=True)
    n = written = skipped = bad = cov10 = 0
    s10 = 0.0
    for line in pc.stdout:
        if line[0] == "#":
            continue
        c = line.split("\t")
        unmeth = int(c[3]) + int(c[7])            # wC + cC
        meth = int(c[4]) + int(c[8])              # wT + cT
        if meth + unmeth == 0:
            skipped += 1
            continue
        n += 1
        m = 100.0 * meth / (meth + unmeth)
        if check:
            b = pb.stdout.readline().split("\t")
            if len(b) < 4 or b[0] != c[0] or int(b[2]) != int(c[1]) or abs(float(b[3]) - m) > 0.001:
                bad += 1
                if bad <= 3:
                    print(f"  {sample}: mismatch at {c[0]}:{c[1]}: Msuite2 line {b}", file=sys.stderr)
        if c[0] not in KEEP:
            continue
        pct = round(m, 2)
        start = int(c[1])
        xmd.stdin.write(f"{c[0]}\t{start}\t{start + 1}\t{num(pct, 2)}\t{meth}\t{unmeth}\n")
        xdf.stdin.write(f"{c[0]}\t{start}\t{start + 1}\t{num(pct / 100, 4)}\n")
        written += 1
        if meth + unmeth >= 10:
            cov10 += 1
            s10 += m
    extra = sum(1 for _ in pb.stdout) if check else 0
    xmd.stdin.close(); xdf.stdin.close()
    ok = (xmd.wait() == 0) & (xdf.wait() == 0) & (pc.wait() == 0) & ((pb.wait() == 0) if check else True) \
         & (bad == 0) & (extra == 0) & (written > 0)
    fmd.close(); fdf.close()
    if not ok:
        for t in tmps:
            os.remove(t)
        return f"{sample}\tFAILED\tCpGs {n}\tmismatches {bad}\tunmatched Msuite2 lines {extra}\tnothing written"
    for t, o in zip(tmps, outs):
        os.rename(t, o)
    checked = f"checked against Msuite2 ({n} CpGs, 0 mismatches)" if check else f"{n} CpGs"
    return (f"{sample}\tok\t{checked}\twritten chr1-22,X {written}\tno C/T {skipped}\t"
            f"mean methylation, CpGs >=10 reads: {s10 / cov10 if cov10 else float('nan'):.2f}%")


if __name__ == "__main__":
    a = sys.argv[1:]
    out, check = None, "--check" in a
    a = [x for x in a if x != "--check"]
    if "--out" in a:
        i = a.index("--out"); out = a[i + 1]; del a[i:i + 2]
        os.makedirs(out, exist_ok=True)
    if len(a) < 2:
        sys.exit("usage: 02_taps_inputs.py <msuite2_dir> <sample> [...] [--out DIR] [--check]")
    ms_dir, samples = a[0], a[1:]
    with Pool(min(len(samples), 8)) as p:
        res = p.map(reformat, [(ms_dir, s, out, check) for s in samples])
    print("\n".join(res))
    sys.exit(0 if all("\tok\t" in r for r in res) else 1)
