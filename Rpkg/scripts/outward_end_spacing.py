#!/usr/bin/env python3
"""How many distinct long-read ends does a gene have, and how far apart?

WHY THIS EXISTS. `exclusive` is null at 3' and 1.88 at 5' in the walk
comparison, and six separate checks agree the 3' signal is absent without any
of them saying why. This measures the thing that explains it.

The count of distinct ends per gene is the SAME at both ends. What differs is
spacing: 5' ends cluster at a median 34 bp apart, 3' ends sit a median 2.2 kb
apart. The matching window is +/-W = 100 bp, so at 5' the window spans a
cluster and the terminus is judged as a promoter REGION -- `exclusive` fires if
any model in the cluster is lineage-restricted. At 3' the window holds one model
end, so `exclusive` must be satisfied by that single model, and it fails
whenever the transcript is present in both lineages at differing abundance,
which is what alternative polyadenylation mostly is.

So `exclusive` is a region-level test at 5' and a site-level test at 3'. The two
ends' figures are therefore not strictly comparable and the 3' value is a floor,
not an estimate. Judging a 5' cluster as one region is deliberate: it is the
right granularity for an alternative-TSS claim.

The model -> gene map is rebuilt here from the GTF. outward_walk_longread.py
dropped its own copy when `share` was removed, since nothing else needed it.

Usage: IG=1 python3 -I outward_end_spacing.py
"""
import os, re, sys, collections, contextlib, io

MAIN = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "outward_walk_longread.py")
sys.argv = [MAIN]
G = {"__name__": "_walk_side", "__file__": MAIN}
sys.stderr.write("reusing %s (its stdout discarded)\n" % MAIN)
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(open(MAIN).read(), MAIN, "exec"), G)
mcoord = G["mcoord"]; cnt = G["cnt"]; O = G["O"]

mg = {}
for ln in open("%s/OUT.transcript_models.gtf" % O):
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "transcript": continue
    t = re.search(r'transcript_id "([^"]+)"', f[8])
    g = re.search(r'gene_id "([^"]+)"', f[8])
    if t and g: mg[t.group(1)] = g.group(1)
print("models with a gene id: %d | models with counts: %d" % (len(mg), len(cnt)))

five = collections.defaultdict(set); three = collections.defaultdict(set)
for m in cnt:
    g = mg.get(m); co = mcoord.get(m)
    if g and co:
        five[g].add(co[2]); three[g].add(co[3])

def summarise(d):
    k = sorted(len(v) for v in d.values() if v)
    gaps = []
    for v in d.values():
        if len(v) < 2: continue
        s = sorted(v); gaps += [s[i+1] - s[i] for i in range(len(s) - 1)]
    return k, sorted(gaps)

print("\n### distinct IsoQuant model ends per gene")
print("  %-6s %8s %9s %9s %9s %12s" % ("end", "genes", "mean", "median", "p90", ">=3 ends"))
for lab, d in (("5'", five), ("3'", three)):
    k, _ = summarise(d)
    print("  %-6s %8d %9.2f %9d %9d %11.1f%%" %
          (lab, len(k), sum(k)/len(k), k[len(k)//2], k[int(.9*len(k))],
           100*sum(1 for x in k if x >= 3)/len(k)))

print("\n### spacing between consecutive ends of the same gene (bp)")
print("  %-6s %10s %9s %9s %12s" % ("end", "n gaps", "median", "p90", "<100bp"))
for lab, d in (("5'", five), ("3'", three)):
    _, gaps = summarise(d)
    print("  %-6s %10d %9d %9d %11.1f%%" %
          (lab, len(gaps), gaps[len(gaps)//2], gaps[int(.9*len(gaps))],
           100*sum(1 for x in gaps if x < 100)/len(gaps)))
