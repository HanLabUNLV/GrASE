#!/usr/bin/env python3
"""List outward flank bins that cannot be measured.

dexseq_count discards a read that overlaps features of two genes as ambiguous,
so a flank bin covered completely by other genes' bins in the same annotation
always counts 0. score_outward_walk.R takes this list as --unmeasurable so such
a flank is not read as an empty flank.

Usage:
  python3 ig_flank_overlap.py <intronparts_manifest.tsv> <out.tsv> [--stranded yes|no|reverse]
    --stranded  as given to dexseq_count; when not "no", only bins on the same
                strand compete for a read (default: reverse)
    --flanks_only
                consider only other genes' flank bins as covering, for flanks
                counted with their own annotation
Output: TSV with columns gene, bin_id
"""
import argparse, bisect, collections

ap = argparse.ArgumentParser()
ap.add_argument("manifest")
ap.add_argument("out")
ap.add_argument("--stranded", default="reverse", choices=["yes", "no", "reverse"])
ap.add_argument("--flanks_only", action="store_true")
a = ap.parse_args()

def is_flank(b):
    return b.startswith("IG_")

rows = []
with open(a.manifest) as fh:
    ix = {k: i for i, k in enumerate(fh.readline().rstrip("\n").split("\t"))}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        rows.append((f[ix["chrom"]], f[ix["strand"]], int(f[ix["start"]]), int(f[ix["end"]]),
                     f[ix["gene"]], f[ix["bin_id"]]))

groups = collections.defaultdict(list)
for r in rows:
    if a.flanks_only and not is_flank(r[5]):
        continue
    groups[(r[0], r[1] if a.stranded != "no" else "")].append(r)

out = []
for v in groups.values():
    v.sort(key=lambda x: x[2])
    starts = [x[2] for x in v]
    maxlen = max(x[3] - x[2] for x in v) + 1
    for a_ in v:
        if not is_flank(a_[5]):
            continue
        lo = bisect.bisect_left(starts, a_[2] - maxlen)
        hi = bisect.bisect_right(starts, a_[3])
        iv = sorted((max(b[2], a_[2]), min(b[3], a_[3])) for b in v[lo:hi]
                    if b[4] != a_[4] and b[3] >= a_[2] and b[2] <= a_[3])
        pos = a_[2]
        for s, e in iv:
            if s > pos:
                break
            pos = max(pos, e + 1)
        if pos > a_[3]:
            out.append((a_[4], a_[5]))

with open(a.out, "w") as fh:
    fh.write("gene\tbin_id\n")
    for g, b in sorted(set(out)):
        fh.write("%s\t%s\n" % (g, b))
print("flank bins: %d   unmeasurable: %d" % (sum(1 for r in rows if is_flank(r[5])), len(set(out))))
