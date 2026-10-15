#!/usr/bin/env python3
"""
17_extract_calls.py                                                  2026-10-15

Turns the genome-wide Mutect2 VCFs from the server (05_variants_server/03 and 05,
copied to emseq/04_variants/vcf/) into one table per case for 17_genomewide_snvs.R.

One row per alternative allele (multi-allelic records split), with per library:
  alt, dp          reads for the allele / all reads (FORMAT AD)
  f1, f1_dp        alt reads / all reads in F1R2 orientation (original top strand, OT)
  f2, f2_dp        the same in F2R1 orientation (original bottom strand, OB)
TAPS turns a C into T only on the strand that carries it, so a C>T change is
checked on F2R1 reads and a G>A change on F1R2 reads (17_genomewide_snvs.R).

Rows kept:
  - every PASS record (the call set and its filter funnel)
  - non-PASS records where at least one clone library has >= 3 alt reads at
    VAF >= 0.15: candidates for the pooled test of whether the clone's mutations
    are already present in earlier samples. FilterMutectCalls marks a call
    'normal_artifact' exactly when an earlier (normal) sample has alt reads, so
    the test cannot be run on PASS calls alone.
Flags from 05_variants_server/05_annotate_collect.sh (PASS calls only): UMAP100, SEGDUP and
the Funcotator gene, class and protein change (GENCODE 43).

usage: 17_extract_calls.py <vcf dir> <out dir>
"""
import csv
import gzip
import os
import sys

import pysam

def clone_libraries(sheet):
    """libraries in which the DDHN clone is large (config/samples.tsv, clone_library = yes), by case"""
    out = {}
    with open(sheet) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("methylation") == "analysed" and r.get("clone_library") == "yes":
                out.setdefault(r["case"], []).append(r["sample"])
    return out


CLONE = clone_libraries(os.environ.get("SAMPLE_SHEET", "config/samples.tsv"))
MIN_ALT, MIN_VAF = 3, 0.15


def roles(vf):
    nrm = [r.value for r in vf.header.records if r.key == "normal_sample"]
    tum = [r.value for r in vf.header.records if r.key == "tumor_sample"]
    return nrm, tum


def flags(path):
    """(chr, pos, ref, alt) -> UMAP100, SEGDUP from the PASS + mappability VCF"""
    out = {}
    for rec in pysam.VariantFile(path):
        for a in rec.alts:
            out[(rec.chrom, rec.pos, rec.ref, a)] = ("UMAP100" in rec.info, "SEGDUP" in rec.info)
    return out


def funco(path):
    """(chr, pos, ref, alt) -> gene, class, protein change (first transcript Funcotator reports)"""
    vf = pysam.VariantFile(path)
    fields = vf.header.info["FUNCOTATION"].description.split("Funcotation fields are: ")[-1].split("|")
    gi, ci, pi = (fields.index(f"Gencode_43_{k}") for k in ("hugoSymbol", "variantClassification", "proteinChange"))
    out = {}
    for rec in vf:
        for a, f in zip(rec.alts, rec.info["FUNCOTATION"]):
            x = f.strip("[]").split("|")
            out[(rec.chrom, rec.pos, rec.ref, a)] = (x[gi], x[ci], x[pi])
    return out


def kind(ref, alt):
    if len(ref) != len(alt):
        return "indel"
    return "SNV" if len(ref) == 1 else "MNV"


def main():
    vdir, odir = sys.argv[1:3]
    os.makedirs(odir, exist_ok=True)
    for case, clone in CLONE.items():
        vf = pysam.VariantFile(os.path.join(vdir, f"{case}.filtered.vcf.gz"))
        nrm, tum = roles(vf)
        libs = nrm + tum
        assert set(clone) <= set(tum), (case, clone, tum)
        fl = flags(os.path.join(vdir, f"{case}.pass.mappability.vcf.gz"))
        fu = funco(os.path.join(vdir, f"{case}.funcotated.vcf.gz"))
        cols = ["chr", "pos", "ref", "alt", "type", "filter", "pass", "popaf", "mmq_alt", "tlod", "umap100", "segdup",
                "gene", "vclass", "protein"] + [f"{s}_{k}" for s in libs for k in ("alt", "dp", "f1", "f1_dp", "f2", "f2_dp")]
        n_rec = n_pass = n_kept = 0
        out = os.path.join(odir, f"calls_{case}.tsv.gz")
        with gzip.open(out + ".tmp", "wt") as fh:
            fh.write("#normals=" + ",".join(nrm) + ";tumours=" + ",".join(tum) + ";clone=" + ",".join(clone) + "\n")
            fh.write("\t".join(cols) + "\n")
            for rec in vf:
                n_rec += 1
                filt = ";".join(rec.filter.keys()) or "."
                ok = filt == "PASS"
                n_pass += ok
                popaf = rec.info.get("POPAF")
                mmq = rec.info.get("MMQ")
                tlod = rec.info.get("TLOD")
                for i, a in enumerate(rec.alts, start=1):
                    per = []
                    for s in libs:
                        smp = rec.samples[s]
                        ad, f1, f2 = smp["AD"], smp["F1R2"], smp["F2R1"]
                        alt = ad[i] or 0
                        per.append((alt, sum(x or 0 for x in ad), f1[i] or 0, (f1[0] or 0) + (f1[i] or 0),
                                    f2[i] or 0, (f2[0] or 0) + (f2[i] or 0)))
                    if not ok:
                        big = any(p[0] >= MIN_ALT and p[1] and p[0] / p[1] >= MIN_VAF for s, p in zip(libs, per) if s in clone)
                        if not big:
                            continue
                    key = (rec.chrom, rec.pos, rec.ref, a)
                    u, d = fl.get(key, ("", ""))
                    g, c, pc = fu.get(key, ("", "", ""))
                    row = [rec.chrom, rec.pos, rec.ref, a, kind(rec.ref, a), filt, int(ok),
                           popaf[i - 1] if popaf else "", mmq[i] if mmq and len(mmq) > i else "",
                           tlod[i - 1] if tlod else "", "" if u == "" else int(u), "" if d == "" else int(d), g, c, pc]
                    row += [x for p in per for x in p]
                    fh.write("\t".join(str(x) for x in row) + "\n")
                    n_kept += 1
        os.replace(out + ".tmp", out)
        print(f"{case}: {n_rec} records, {n_pass} PASS; {n_kept} rows written to {out}")


if __name__ == "__main__":
    main()
