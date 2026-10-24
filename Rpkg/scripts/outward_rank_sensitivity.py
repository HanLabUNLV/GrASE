#!/usr/bin/env python3
"""SIDE ANALYSIS: is the 3' stops/continues inversion an artefact of the rank
BINNING rather than a property of the data?

outward_walk_longread.py standardises on bins 0,1,2,3,4,5-9,10+. Long-read end
support falls very steeply with rank at the 3' end (75.8% at rank 0 down to
5.3% at 10+), and `10+` pools rank 10 with rank 80. Called termini sit at the
SHALLOW end of that bin while the background parts filling it run much deeper,
so the bin's mean rate overestimates the expected for any called group, and
overestimates it most for the group with the most mass in the bin. At 3' that
group is walk_continues (median rank 11 against 8 for stops_measured), which is
exactly the direction needed to manufacture the observed inversion.

This recomputes obs/exp with NO binning -- expected summed over EXACT rank --
and prints the three schemes side by side. Everything else is reused from
outward_walk_longread.py.

Usage: IG=1 python3 -I outward_rank_sensitivity.py
"""
import os, sys, csv, collections, contextlib, io

HERE  = os.path.dirname(os.path.abspath(__file__))
MAIN  = os.path.join(HERE, "outward_walk_longread.py")
SPLIT = os.path.join(HERE, "outward_stops_split.py")
sys.argv = [MAIN] + sys.argv[1:]
G = {"__name__": "_walk_side", "__file__": MAIN}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(open(MAIN).read(), MAIN, "exec"), G)
for nm in ("build", "near", "verd", "meta", "plen", "npart", "chrom_of",
           "rank_anchor", "strand_of", "DEFS", "PAIRS", "IG", "have"):
    globals()[nm] = G[nm]

MEASURED  = ("drop: outward tested, uncalled", "drop: outward empty")
UNPOWERED = ("drop: outward near-empty", "drop: outward untestable")
def stratum(r):
    pc = r["part_called"] == "TRUE"; tc = r["tile_called"] == "TRUE"
    if not pc:  return "no_claim"
    if tc:      return "walk_continues"
    tv = r["term_verdict"]
    if tv in MEASURED:  return "stops_measured"
    if tv in UNPOWERED: return "stops_unpowered"
    return None
GRPS = ["stops_measured", "stops_unpowered", "walk_stops_pooled",
        "walk_continues", "no_claim"]

# three binnings of the same exact ranks
def b_coarse(r): return str(r) if r < 5 else ("5-9" if r < 10 else "10+")
def b_fine(r):
    if r < 5:  return str(r)
    if r < 10: return "5-9"
    if r < 15: return "10-14"
    if r < 20: return "15-19"
    if r < 30: return "20-29"
    if r < 50: return "30-49"
    return "50+"
def b_exact(r): return str(r)
SCHEMES = (("coarse (published)", b_coarse), ("fine", b_fine), ("exact rank", b_exact))

# BACKGROUND: every exonic part of every gene in the panel, at its exact rank,
# duplicated across schemes. Same loop as the main script's, which bins as it
# goes and so cannot be reused here.
bgn = collections.Counter(); bgh = collections.Counter()
for trt, ref in PAIRS:
    ds = build(trt, ref)
    for g in {r["gene"] for r in verd.get("%s_vs_%s" % (trt, ref), [])} & set(plen):
        P = sorted(plen[g], key=int); N = len(P); strand = strand_of[g]
        for kind, e in (("5prime", "5"), ("3prime", "3")):
            low = (kind == "5prime") == (strand == "+")
            for idx, pn in enumerate(P):
                s0, e0 = plen[g][pn]
                pos = s0 if low else e0
                rk = idx if low else N-1-idx
                bgn[(kind, rk)] += 1
                for d in DEFS:
                    if near(ds[(d, e)], chrom_of[g], strand, pos): bgh[(kind, d, rk)] += 1

obs = collections.Counter(); nrk = collections.Counter(); ncell = collections.Counter()
obs_rk = collections.Counter()
for trt, ref in PAIRS:
    ds = build(trt, ref)
    for r in verd.get("%s_vs_%s" % (trt, ref), []):
        key = (r["gene"], r["event"], r["side"])
        anc = rank_anchor.get((key, int(r["trank"])))
        if anc is None: continue
        k = (r["gene"], r["event"], r["side"], anc)
        if k not in have or k not in meta: continue
        g = r["gene"]
        if g not in plen: continue
        co = plen[g].get(anc)   # rank_anchor is canonical %03d
        if not co: continue
        grp = stratum(r)
        if grp is None: continue
        kind, strand = meta[k]
        low = (kind == "5prime") == (strand == "+")
        pos = co[0] if low else co[1]
        rk = max((int(anc)-1) if low else (npart[g]-int(anc)), 0)
        e = "5" if kind == "5prime" else "3"
        groups = [grp] + (["walk_stops_pooled"] if grp.startswith("stops_") else [])
        for gr in groups:
            ncell[(kind, gr)] += 1; nrk[(kind, gr, rk)] += 1
            for d in DEFS:
                if near(ds[(d, e)], chrom_of[g], strand, pos):
                    obs[(kind, gr, d)] += 1; obs_rk[(kind, gr, d, rk)] += 1

def oe(kind, gr, d, binf):
    # collapse the exact-rank background to this scheme's bins, then spend the
    # group's own exact-rank counts against the bin rate it would have been given
    bn = collections.Counter(); bh = collections.Counter()
    for (kk, rk), n in bgn.items():
        if kk != kind: continue
        bn[binf(rk)] += n; bh[binf(rk)] += bgh[(kind, d, rk)]
    exp = 0.0
    for (kk, gg, rk), n in nrk.items():
        if kk != kind or gg != gr: continue
        b = binf(rk)
        exp += n * (bh[b]/bn[b] if bn[b] else 0)
    o = obs[(kind, gr, d)]
    return o, exp, (o/exp) if exp else float("nan")

print("\n### rank-binning sensitivity of the stops/continues ratio   IG=%s" % IG)
for name, binf in SCHEMES:
    print("\n  -- %s" % name)
    print("     %-8s %-18s %8s | %-24s" % ("end", "group", "n", "o/e  exists / share / excl"))
    R = {}
    for kind in ("5prime", "3prime"):
        for gr in GRPS:
            n = ncell[(kind, gr)]
            if not n: continue
            vals = [oe(kind, gr, d, binf)[2] for d in DEFS]
            R[(kind, gr)] = vals
            print("     %-8s %-18s %8d | %6.2f %6.2f %6.2f" % (kind, gr, n, *vals))
    print("     %-8s %-18s %8s | %-24s" % ("", "o/e RATIO (invalid)", "", ""))
    for kind in ("5prime", "3prime"):
        b = R.get((kind, "walk_continues"))
        for gr in ("stops_measured", "stops_unpowered", "walk_stops_pooled"):
            a = R.get((kind, gr))
            if not a or not b: continue
            print("     %-8s %-18s %8d | %6.2fx %5.2fx %5.2fx" %
                  (kind, gr, ncell[(kind, gr)], *[x/y for x, y in zip(a, b)]))

## THE CEILING PROBLEM, and a Mantel-Haenszel estimator that does not have it.
##
## obs/exp is a ratio of RATES, so its headroom depends on the baseline. At 3'
## rank 0 the long-read background is 75.8%, which caps any o/e at 1/0.758 =
## 1.32x; at rank 10+ the background is 5.3% and o/e can reach 19x. Pooling obs
## and exp as counts across ranks therefore weights the result toward whichever
## group sits in the high-headroom ranks. At 3' that is walk_continues (median
## rank 11), while the walk_stops set piles up near the gene end (27.2% of
## stops_unpowered at rank 0) where NO enrichment is expressible. That is a
## ceiling effect, not a binning artefact, and finer bins cannot fix it.
##
## A Mantel-Haenszel odds ratio stratified by EXACT rank compares the two groups
## within each rank and combines on the odds scale, which is not bounded by the
## baseline rate. This is the same estimator the locality analysis already uses
## (stratified by gene there, by rank here).
print("\n### Mantel-Haenszel OR, stops vs continues, stratified by EXACT rank")
print("    OR > 1 means the walk_stops group has the long-read end more often AT THE")
print("    SAME POSITION IN THE GENE. Unlike o/e it is not capped by the baseline.")
print("  %-8s %-18s %8s | %-26s" % ("end", "group", "n", "MH OR  exists / share / excl"))
for kind in ("5prime", "3prime"):
    for gr in ("stops_measured", "stops_unpowered", "walk_stops_pooled"):
        if not ncell[(kind, gr)]: continue
        cols = []
        for d in DEFS:
            num = den = 0.0
            ranks_here = {rk for (kk, gg, rk) in nrk if kk == kind and gg in (gr, "walk_continues")}
            for rk in ranks_here:
                a = obs_rk[(kind, gr, d, rk)]; n1 = nrk[(kind, gr, rk)]
                c = obs_rk[(kind, "walk_continues", d, rk)]; n2 = nrk[(kind, "walk_continues", rk)]
                b = n1 - a; dd = n2 - c; n = n1 + n2
                if n == 0: continue
                num += a * dd / n; den += b * c / n
            cols.append(num/den if den else float("nan"))
        print("  %-8s %-18s %8d | %8.2f %8.2f %8.2f" % (kind, gr, ncell[(kind, gr)], *cols))

print("\n### why: headroom by rank at each end (background rate for `exists`)")
print("  %-8s %6s %10s %10s %12s" % ("end", "rank", "bg rate", "max o/e", "n bg parts"))
for kind in ("5prime", "3prime"):
    for rk in (0, 1, 2, 5, 10, 20):
        n = bgn[(kind, rk)]
        if not n: continue
        rate = bgh[(kind, "exists", rk)] / n
        print("  %-8s %6d %9.1f%% %10s %12d" %
              (kind, rk, 100*rate, ("%.2fx" % (1/rate)) if rate else "inf", n))
print("\n### where each group sits: share of its termini at ranks 0-1 (the saturated zone)")
for kind in ("5prime", "3prime"):
    for gr in GRPS:
        n = ncell[(kind, gr)]
        if not n: continue
        lo = sum(v for (kk, gg, rk), v in nrk.items() if kk == kind and gg == gr and rk <= 1)
        print("  %-8s %-18s %8d  ranks 0-1: %5.1f%%" % (kind, gr, n, 100*lo/n))
