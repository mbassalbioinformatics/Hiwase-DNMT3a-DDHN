#!/usr/bin/env python3
###############################################################################
# taps_mask_methylation.py                                          2026-10-13
#
# Makes Msuite2's TAPS alignments usable by ordinary variant callers (Mutect2,
# GATK CNV / allelic counts).
#
# TAPS turns a methylated C into T on the DNA strand that carries it. At a
# reference CpG the C is read by the original top strand (OT) and the G by the
# original bottom strand (OB); those bases report methylation, not genotype, and
# would look like C>T / G>A variants everywhere. This script turns exactly those
# bases into N with base quality 0, so callers take the genotype at a CpG from
# the other strand. Both are needed: Mutect2 raises any base quality below 18
# to 6 before genotyping, so a quality-0 T would still count as weak evidence
# for T, while an N matches every haplotype equally and carries no evidence;
# pileup tools (GetPileupSummaries, CollectAllelicCounts) skip N. The MD tag of
# a masked read is dropped (it would disagree with the N bases). Every other
# base is untouched; TAPS background conversion outside CpGs (about 0.2% of C)
# is handled after calling, by the strand check of the C>T / G>A calls.
#
# Directional library (checked on the DNMT3A region, 8 Oct 2026):
#   OT = read 1 forward or read 2 reverse;  OB = read 1 reverse or read 2 forward
#
# One contig per call, so a sample runs as many processes in parallel. The
# output header lists the contigs of the GATK reference dictionary in its order
# (some hg38 dictionaries are in lexical order) and one read group named after
# the sample. Secondary and supplementary alignments are dropped; a read whose
# mate lies on a contig outside the dictionary keeps its alignment and its mate
# is marked unmapped.
#
# USAGE  taps_mask_methylation.py <in.bam> <ref.fa> <ref.dict> <contig> <sample> <out.bam>
###############################################################################
import sys
import numpy as np
import pysam

MASK_VERSION = "2"   # in the @PG header line; 02_prepare_bams.sh rebuilds BAMs masked by an older version (1 = quality 0 only)


def dict_contigs(path):
    sq = []
    for line in open(path):
        if line.startswith("@SQ"):
            f = dict(x.split(":", 1) for x in line.rstrip("\n").split("\t")[1:])
            sq.append({"SN": f["SN"], "LN": int(f["LN"])})
    return sq


def cpg_c_positions(fa, contig):
    """0-based position of the C of every CG on the reference contig."""
    seq = fa.fetch(contig).upper()
    pos, i = [], seq.find("CG")
    while i != -1:
        pos.append(i)
        i = seq.find("CG", i + 1)
    return np.array(pos, dtype=np.int64)


def main():
    inbam, ref, dictf, contig, sample, out = sys.argv[1:7]
    sq = dict_contigs(dictf)
    new_tid = {d["SN"]: i for i, d in enumerate(sq)}
    if contig not in new_tid:
        sys.exit(f"{contig} is not in {dictf}")
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq,
           "RG": [{"ID": sample, "SM": sample, "LB": sample, "PL": "ILLUMINA"}],
           "PG": [{"ID": "taps_mask_methylation", "PN": "taps_mask_methylation.py", "VN": MASK_VERSION, "CL": " ".join(sys.argv)}]}
    src = pysam.AlignmentFile(inbam, "rb")
    if contig not in src.references:
        sys.exit(f"{contig} is not in {inbam}")
    src_names = src.references
    dst = pysam.AlignmentFile(out, "wb", header=hdr)
    fa = pysam.FastaFile(ref)
    cpg = cpg_c_positions(fa, contig)
    ctid = new_tid[contig]
    n_reads = n_masked = n_mate_off = 0

    for r in src.fetch(contig):
        if r.is_secondary or r.is_supplementary:
            continue
        n_reads += 1
        if not r.is_unmapped and r.query_qualities is not None:
            ot = ((r.is_read1 and not r.is_reverse) or (r.is_read2 and r.is_reverse)) if r.is_paired else not r.is_reverse
            s, e = r.reference_start, r.reference_end
            lo, hi = np.searchsorted(cpg, s - 1), np.searchsorted(cpg, e)
            if hi > lo:
                tgt = cpg[lo:hi] if ot else cpg[lo:hi] + 1           # C on the top strand, G on the bottom strand
                tgt = tgt[(tgt >= s) & (tgt < e)]
                if len(tgt):
                    cig = r.cigartuples
                    if all(op in (0, 4, 5, 7, 8) for op, _ in cig):  # no insertion, deletion or skip
                        lead = cig[0][1] if cig[0][0] == 4 else (cig[1][1] if cig[0][0] == 5 and len(cig) > 1 and cig[1][0] == 4 else 0)
                        qps = [int(rp) - s + lead for rp in tgt]
                    else:
                        rq = {rp: qp for qp, rp in r.get_aligned_pairs(matches_only=True)}
                        qps = [rq[int(rp)] for rp in tgt if int(rp) in rq]
                    if qps:
                        q = r.query_qualities                       # read before the sequence is replaced (that clears it)
                        seq = bytearray(r.query_sequence, "ascii")
                        for qp in qps:
                            seq[qp] = 78                            # N
                            q[qp] = 0
                        r.query_sequence = seq.decode("ascii")
                        r.query_qualities = q
                        if r.has_tag("MD"):
                            r.set_tag("MD", None)
                        n_masked += len(qps)
        # reference ids in the new header
        r.reference_id = ctid
        if r.is_paired and r.next_reference_id >= 0:
            mname = src_names[r.next_reference_id]
            if mname in new_tid:
                r.next_reference_id = new_tid[mname]
            else:
                r.next_reference_id = -1
                r.next_reference_start = -1
                r.template_length = 0
                r.mate_is_unmapped = True
                r.mate_is_reverse = False
                r.is_proper_pair = False
                n_mate_off += 1
        r.set_tag("RG", sample)
        dst.write(r)
    dst.close()
    print(f"{sample}\t{contig}\treads {n_reads}\tmasked bases {n_masked}\tmates off-dictionary {n_mate_off}", file=sys.stderr)


if __name__ == "__main__":
    main()
