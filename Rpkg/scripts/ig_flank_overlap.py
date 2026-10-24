import collections, os
D = os.environ.get("DICE", "/mnt/data1/home/mirahan/DICE")
bins = []
with open(f"{D}/intronparts_ig/intronparts_manifest.tsv") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {k: i for i, k in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        b = f[ix["bin_id"]]
        if "IG_below" not in b and "IG_above" not in b:
            continue
        bins.append((f[ix["chrom"]], int(f[ix["start"]]), int(f[ix["end"]]),
                     f[ix["gene"]], b))
print("IG bins: %d" % len(bins))

byc = collections.defaultdict(list)
for b in bins:
    byc[b[0]].append(b)

hit = set()           # (gene, bin) overlapped by ANOTHER gene's IG bin
ov_bp = collections.defaultdict(int)
for c, v in byc.items():
    v.sort(key=lambda x: x[1])
    for i, a in enumerate(v):
        for b in v[i+1:]:
            if b[1] > a[2]:
                break
            if a[3] == b[3]:
                continue          # same gene: not ambiguous
            n = min(a[2], b[2]) - max(a[1], b[1]) + 1
            if n <= 0:
                continue
            hit.add((a[3], a[4])); hit.add((b[3], b[4]))
            ov_bp[(a[3], a[4])] += n
            ov_bp[(b[3], b[4])] += n

print("IG bins overlapped by another gene's IG bin: %d (%.1f%%)"
      % (len(hit), 100.0*len(hit)/len(bins)))

# how much of the overlapped bin's 300 bp is lost
full = sum(1 for k, n in ov_bp.items() if n >= 300)
half = sum(1 for k, n in ov_bp.items() if 150 <= n < 300)
print("  of those, fully covered (>=300bp): %d (%.1f%% of overlapped)" % (full, 100.0*full/len(hit)))
print("            >=half covered:          %d (%.1f%%)" % (full+half, 100.0*(full+half)/len(hit)))

# per gene: does the gene have at least one clean flank
g_all = collections.defaultdict(int)
g_hit = collections.defaultdict(int)
for c, s, e, g, b in bins:
    g_all[g] += 1
    if (g, b) in hit:
        g_hit[g] += 1
genes = len(g_all)
both_clean = sum(1 for g in g_all if g_hit[g] == 0)
one_clean  = sum(1 for g in g_all if g_hit[g] == 1)
none_clean = sum(1 for g in g_all if g_hit[g] == g_all[g])
print("genes: %d   both flanks clean %d (%.1f%%)   one flank hit %d (%.1f%%)   both hit %d (%.1f%%)"
      % (genes, both_clean, 100.0*both_clean/genes, one_clean, 100.0*one_clean/genes,
         none_clean, 100.0*none_clean/genes))
