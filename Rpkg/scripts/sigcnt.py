#!/usr/bin/env python3
"""Bipartition counts from read signatures, each read counted once per set.

Reads come from inst/python/dexseq_count_sig.py, which writes for every sample a
<out>.sig file: gene, the exonic parts a read touches, the introns it crosses,
and the number of reads with that signature. A set of parts (or of introns) is
counted once per read by summing the signatures that touch it. Summing per-part
counts instead counts a read once for every part of the set it overlaps.

Two steps:

  split  stream the per-sample .sig files into gene buckets
         sigcnt.py split --samples samples.tsv --outdir buckets [--nbuckets 256]
         samples.tsv: path<TAB>sample<TAB>group, one line per sample

  count  for each bucket, count every bipartition side of its genes
         sigcnt.py count --buckets buckets --samples samples.tsv --bipartdir <filtered>
                         --junctiondir <junction sets> --type internal|terminal
                         --outdir <counts> [--cores N]

`count` writes, per gene, <gene>.bipartition.exsj.<type>.cnt.txt in the merged
exon+junction layout exontest.R reads (merge_exon_sj_counts.R writes the same),
<gene>.bipartition.txt with the event number and the mean counts of each set (the
annotation exoncnt.R writes, which exontest.R joins on gene and event), and
measured_junctions.txt. As exoncnt.R does, a gene is written only if some event
has reads on a distinct set. A side is
counted on its exonic distinct parts when it has any (source "exon"), otherwise
on its distinct introns from bipartition_junction_sets.R (source "sj"); S is
counted on ref_ex_part. A read touching both a side and S counts in both, as it
does with per-part counts.
"""
import argparse, collections, gzip, os, sys, zlib
from multiprocessing import Pool
import numpy as np


def read_samples(path):
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and f[0]:
            out.append((f[0], f[1], f[2]))
    return out


def bucket_of(gene, n):
    return zlib.crc32(gene.encode()) % n


## ---- split ------------------------------------------------------------------
def split(a):
    samples = read_samples(a.samples)
    os.makedirs(a.outdir, exist_ok=True)
    fh = [gzip.open(os.path.join(a.outdir, "bucket%03d.tsv.gz" % i), "wt", compresslevel=1)
          for i in range(a.nbuckets)]
    for k, (path, sample, group) in enumerate(samples):
        with open(path) as f:
            next(f)
            for line in f:
                gene, parts, junc, cnt = line.rstrip("\n").split("\t")
                gene = gene.strip('"')
                fh[bucket_of(gene, a.nbuckets)].write(
                    "%s\t%d\t%s\t%s\t%s\n" % (gene, k, parts.replace('"', ""), junc, cnt))
        if (k + 1) % 50 == 0:
            sys.stderr.write("%d of %d samples\n" % (k + 1, len(samples)))
    for f in fh:
        f.close()


## ---- count ------------------------------------------------------------------
def parts_of(cell):
    if cell in ("", "NA"):
        return None
    return frozenset(int(x[1:]) for x in cell.split(","))


def introns_of(cell):
    """junction ids chr:start:end, "|" or "," separated; the end may be the intron's
    last base or the next exon's first, so both are accepted (as sum_sj_counts does)"""
    if cell in ("", "NA"):
        return None
    out = set()
    for j in cell.replace("|", ",").split(","):
        f = j.strip().split(":")
        if len(f) >= 3:
            s, e = int(f[1]), int(f[2])
            out.add("%d-%d" % (s, e))
            out.add("%d-%d" % (s, e - 1))
    return frozenset(out)


def read_table(path):
    with open(path) as f:
        h = f.readline().rstrip("\n").split("\t")
        return [dict(zip(h, l.rstrip("\n").split("\t"))) for l in f]


def count_bucket(job):
    b, a, nsamp, samples, genes = job
    ## gene -> signature -> per-sample counts
    sig = collections.defaultdict(dict)
    with gzip.open(os.path.join(a.buckets, "bucket%03d.tsv.gz" % b), "rt") as f:
        for line in f:
            gene, k, parts, junc, cnt = line.rstrip("\n").split("\t")
            key = (parts, junc)
            v = sig[gene].get(key)
            if v is None:
                v = sig[gene][key] = np.zeros(nsamp)
            v[int(k)] += int(cnt)
    n_genes = n_sj = 0
    measured = []
    ## every gene with a bipartition file in this bucket; written only if it has reads
    for gene in genes:
        sg = sig.get(gene, {})
        bp = os.path.join(a.bipartdir, "%s.bipartition.%s.txt" % (gene, a.type))
        rows = read_table(bp)
        if not rows:
            continue
        jfile = os.path.join(a.junctiondir, "%s.bipartition.%s.junctions.txt" % (gene, a.type))
        jrows = {int(r["event"]): r for r in read_table(jfile)} if os.path.exists(jfile) else {}
        keys = list(sg)
        kparts = [frozenset(int(x) for x in p.split("+")) for p, _ in keys]
        kjunc = [frozenset(j.split(",")) if j else frozenset() for _, j in keys]
        mat = np.vstack([sg[k] for k in keys]) if keys else np.zeros((0, nsamp))
        cache = {}

        def by_parts(ps):
            if ps not in cache:
                hit = [i for i, kp in enumerate(kparts) if kp & ps]
                cache[ps] = mat[hit].sum(axis=0) if hit else np.zeros(nsamp)
            return cache[ps]

        def by_introns(js):
            if ("j", js) not in cache:
                hit = [i for i, kj in enumerate(kjunc) if kj & js]
                cache[("j", js)] = mat[hit].sum(axis=0) if hit else np.zeros(nsamp)
            return cache[("j", js)]

        lines, ann, gene_measured, any_reads = [], [], [], False
        hdr = ("gene\tevent\tsource\tsink\tref_ex_part\tsetdiff1\tsetdiff2\tsample\tref\t"
               "diff1\tdiff2\tgroups\tdiff1_source\tdiff2_source\tintron_distinct1\t"
               "intron_distinct2\tintron_shared\n")
        for ev, r in enumerate(rows, start=1):
            jr = jrows.get(ev, {})
            cells = {}
            for k in (1, 2):
                sd = r.get("setdiff%d" % k, "")
                ps = parts_of(sd)
                ij = jr.get("intron_distinct%d" % k, "NA")
                if ps:
                    cells[k] = (by_parts(ps), "exon")
                elif introns_of(ij):
                    cells[k] = (by_introns(introns_of(ij)), "sj")
                    gene_measured.append("%s\t%d\tdiff%d_vs_ref\t%s" % (gene, ev, k, ij))
                else:
                    cells[k] = (None, "exon")
            sref = parts_of(r.get("ref_ex_part", ""))
            ref = by_parts(sref) if sref else None
            sd1 = r.get("setdiff1", "") or "NA"
            sd2 = r.get("setdiff2", "") or "NA"
            mean = lambda v: "NaN" if v is None else repr(float(v.mean()))
            if any(c[0] is not None and c[0].sum() > 0 for c in cells.values()):
                any_reads = True
            ann.append([str(ev), gene, r["source"], r["sink"], r.get("ref_ex_part", "") or "NA", sd1, sd2,
                        r.get("transcripts1", ""), r.get("transcripts2", ""), r.get("path1", ""),
                        r.get("path2", ""), str(ev), mean(ref), mean(cells[1][0]), mean(cells[2][0])])
            for i, (path, sample, group) in enumerate(samples):
                v = [("NA" if c is None else "%d" % c[i]) for c in
                     (ref, cells[1][0], cells[2][0])]
                lines.append("\t".join([gene, str(ev), r["source"], r["sink"],
                                    r.get("ref_ex_part", "") or "NA", sd1, sd2, sample,
                                    v[0], v[1], v[2], group, cells[1][1], cells[2][1],
                                    jr.get("intron_distinct1", "NA"),
                                    jr.get("intron_distinct2", "NA"),
                                    jr.get("intron_shared", "NA")]) + "\n")
        if not any_reads:
            continue
        with open(os.path.join(a.outdir, "%s.bipartition.exsj.%s.cnt.txt" % (gene, a.type)), "w") as fo:
            fo.write(hdr)
            fo.writelines(lines)
        with open(os.path.join(a.outdir, "%s.bipartition.txt" % gene), "w") as fo:
            fo.write("gene\tsource\tsink\tref_ex_part\tsetdiff1\tsetdiff2\ttranscripts1\ttranscripts2\t"
                     "path1\tpath2\tevent\tref_mean\tdiff1_mean\tdiff2_mean\n")
            for x in ann:
                fo.write("\t".join(x) + "\n")
        measured.extend(gene_measured)
        n_sj += len(gene_measured)
        n_genes += 1
    return n_genes, n_sj, measured


def count(a):
    samples = read_samples(a.samples)
    os.makedirs(a.outdir, exist_ok=True)
    nb = len([f for f in os.listdir(a.buckets) if f.startswith("bucket")])
    suffix = ".bipartition.%s.txt" % a.type
    by_bucket = collections.defaultdict(list)
    for f in os.listdir(a.bipartdir):
        if f.endswith(suffix):
            g = f[:-len(suffix)]
            by_bucket[bucket_of(g, nb)].append(g)
    jobs = [(b, a, len(samples), samples, sorted(by_bucket[b])) for b in range(nb)]
    with Pool(a.cores) as p:
        res = p.map(count_bucket, jobs, chunksize=1)
    with open(os.path.join(a.outdir, "measured_junctions.txt"), "w") as fo:
        fo.write("gene\tevent\tcomparison\tjunctions\n")
        for _, _, m in res:
            for line in m:
                fo.write(line + "\n")
    sys.stderr.write("genes written: %d; junction-measured sides: %d\n" %
                     (sum(r[0] for r in res), sum(r[1] for r in res)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("split")
    s.add_argument("--samples", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--nbuckets", type=int, default=256)
    c = sub.add_parser("count")
    for k in ("buckets", "samples", "bipartdir", "junctiondir", "type", "outdir"):
        c.add_argument("--" + k, required=True)
    c.add_argument("--cores", type=int, default=8)
    a = ap.parse_args()
    split(a) if a.cmd == "split" else count(a)
