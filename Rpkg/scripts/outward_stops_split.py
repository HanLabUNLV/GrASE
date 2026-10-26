#!/usr/bin/env python3
"""SIDE ANALYSIS, not part of outward_walk_longread.py.

TWO things, and the second one supersedes a test that was wrong.

(1) Split `walk_stops` by whether its verdict rests on a MEASURED outward flank
    or an UNPOWERED one. `walk_stops` is `part_called AND NOT flank_called`, and
    "not called" has two very different causes: the flank was quantified and the
    test said no (a measurement), or it never cleared the read floor / was never
    tested (an absence of measurement scored as a boundary). At 3' the second
    kind dominates, so the pooled figure may be reporting the unpowered half.

(2) Compare `walk_stops` against `walk_continues` with a MANTEL-HAENSZEL odds
    ratio stratified by EXACT rank, which is the correct estimator for that
    comparison. It replaces the ratio-of-two-obs/exp-values this script used
    before 2026-10-22; why that is invalid is in outward_walk_longread.py's
    header and is not restated here. obs/exp against the background is still
    printed for the two strata, because "is this group enriched at all" is a
    fair question; what must not be done is divide two of those numbers.

Strata of `walk_stops` come from the scorer's own `term_verdict`, verified to
partition it exactly:
    MEASURED    drop: outward tested, uncalled   flank tested, not called
                drop: outward empty              flank measured at exactly zero
    UNPOWERED   drop: outward near-empty         flank below the read floor
                drop: outward untestable         flank never tested
    (excluded)  outward flank unmeasurable        flank ambiguity, positive FALSE

References, the matching window and the background are reused from
outward_walk_longread.py by executing it with its stdout discarded, so there is
one definition of the reference, not two. Standardisation here is on EXACT
rank rather than that script's bins.

NON-REDUNDANCY RULE. Only the two walk_stops strata are printed in the o/e and
MH tables: outward_walk_longread.py already prints walk_stops, walk_continues
and no_claim panel-wide, and repeating them here produced two copies of the
same figures. The whole-group rows are still ACCUMULATED, because the
headroom block needs their rank distribution.

Usage: IG=1 python3 -I outward_stops_split.py [ratio] [window] [dshare] [minpres] [mingene] [mode]
"""
import os, sys, csv, math, collections, contextlib, io, bisect

MAIN = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "outward_walk_longread.py")
sys.argv = [MAIN] + sys.argv[1:]
G = {"__name__": "_walk_side", "__file__": MAIN}
sys.stderr.write("reusing %s (its stdout discarded)\n" % MAIN)
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(open(MAIN).read(), MAIN, "exec"), G)
for nm in ("build", "near", "hit", "verd", "meta", "plen", "npart", "have",
           "chrom_of", "pair", "strand_of", "DEFS", "PAIRS", "IG",
           "mh_or"):
    globals()[nm] = G[nm]

MEASURED  = ("drop: outward tested, uncalled", "drop: outward empty")
UNPOWERED = ("drop: outward near-empty", "drop: outward untestable")
def stratum(r):
    pc = r["part_called"] == "TRUE"; tc = r["flank_called"] == "TRUE"
    if not pc:  return "no_claim"
    if tc:      return "walk_continues"
    tv = r["term_verdict"]
    if tv in MEASURED:  return "stops_measured"
    if tv in UNPOWERED: return "stops_unpowered"
    return None
GRPS  = ["stops_measured", "stops_unpowered", "walk_stops", "walk_continues", "no_claim"]
## PRINTED strata. GRPS is what gets accumulated; these two lists are what
## gets printed, and they exclude whatever outward_walk_longread.py prints.
OE_GRPS = ["stops_measured", "stops_unpowered"]
STOPS = ["stops_measured", "stops_unpowered"]
ENDS  = ("5prime", "3prime")

## BACKGROUND at EXACT rank: every exonic part of every gene in the analysis.
## This is the "expected" in obs/exp -- a random exonic part at the same
## position in the gene, NOT any verdict group.
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
                    ## dpi None: a background part carries no claim to match
                    if hit(ds, d, e, chrom_of[g], strand, pos, None):
                        bgh[(kind, d, rk)] += 1

def nearest_signed(dmap, c, st, pos, kind):
    """Signed distance to the nearest long-read end, POSITIVE = outward. `near()`
    is binary within W bp, so absence and displacement score alike; sign and
    magnitude separate them."""
    v = dmap.get((c, st))
    if not v: return None
    i = bisect.bisect_left(v, pos)
    cand = [v[j] for j in (i-1, i) if 0 <= j < len(v)]
    if not cand: return None
    p = min(cand, key=lambda q: abs(q - pos))
    return (p - pos) * (1 if ((kind == "3prime") == (st == "+")) else -1)

n_g = collections.Counter(); h_g = collections.Counter()
n_rk = collections.Counter(); h_rk = collections.Counter()
vmix = collections.Counter(); sig = collections.defaultdict(list)
rk_of = collections.defaultdict(list)

for trt, ref in PAIRS:
    ds = build(trt, ref)
    for r in verd.get("%s_vs_%s" % (trt, ref), []):
        key = (r["gene"], r["event"], r["side"])
        anc = "%03d" % int(r["anchor"])
        k = (r["gene"], r["event"], r["side"], anc)
        if k not in have or k not in meta: continue
        g = r["gene"]
        if g not in plen: continue
        co = plen[g].get(anc)   # canonical %03d
        if not co: continue
        grp = stratum(r)
        if grp is None: continue
        kind, strand = meta[k]
        low = (kind == "5prime") == (strand == "+")
        pos = co[0] if low else co[1]
        rk  = max((int(anc)-1) if low else (npart[g]-int(anc)), 0)
        e   = "5" if kind == "5prime" else "3"
        ## direction of the call; `hd` not `hit`, which is the resolver itself
        try: dpi = float(r["part_dpi"])
        except (KeyError, ValueError): dpi = None
        hd = {d: hit(ds, d, e, chrom_of[g], strand, pos, dpi) for d in DEFS}
        # pool the two stops strata into `walk_stops` as well, so the published
        # quantity and its two halves are computed the same way
        for gr in ([grp] + (["walk_stops"] if grp.startswith("stops_") else [])):
            n_g[(kind, gr)] += 1; n_rk[(kind, gr, rk)] += 1
            for d in DEFS:
                if hd[d]: h_g[(kind, gr, d)] += 1; h_rk[(kind, gr, d, rk)] += 1
        vmix[(kind, grp, r["term_verdict"])] += 1
        rk_of[(kind, grp)].append(rk)
        ## the exists index is split by direction now; take the claimed half
        xi = ("exists_trt" if (dpi or 0) > 0 else "exists_ref") if dpi else "exists_trt"
        sd = nearest_signed(ds[(xi, e)], chrom_of[g], strand, pos, kind)
        if sd is not None: sig[(kind, grp)].append(sd)

def oe(kind, gr, d):
    obs = h_g[(kind, gr, d)]
    exp = sum(n * (bgh[(kind, d, rk)] / bgn[(kind, rk)] if bgn[(kind, rk)] else 0)
              for (kk, gg, rk), n in n_rk.items() if kk == kind and gg == gr)
    return obs, exp, (obs / exp) if exp else float("nan")

## mh_or() comes from outward_walk_longread.py, where the same test compares
## walk_stops with walk_continues panel-wide. Here it is applied to the two
## halves of walk_stops separately, which is the point of this script.

print("\n### walk_stops split by whether the outward flank was MEASURED   IG=%s" % IG)
print("### enrichment against the rank-matched BACKGROUND (a random exonic part at")
print("### the same rank). DO NOT divide two of these -- see the MH block below.")
print("### walk_stops / walk_continues / no_claim are in outward_walk_longread.py.")
print(("  %-8s %-16s %8s" + " | %-22s" * len(DEFS)) %
      (("end", "stratum", "n") + tuple("`%s`" % d for d in DEFS)))
print("  %-8s %-16s %8s | %6s %6s %7s | %6s %6s %7s | %6s %6s %7s" %
      ("", "", "", "obs", "exp", "o/e", "obs", "exp", "o/e", "obs", "exp", "o/e"))
for kind in ENDS:
    for gr in OE_GRPS:
        n = n_g[(kind, gr)]
        if not n: continue
        out = []
        for d in DEFS:
            obs, exp, r = oe(kind, gr, d)
            out += ["%5.1f%%" % (100*obs/n), "%5.1f%%" % (100*exp/n), "%6.2f" % r]
        print(("  %-8s %-16s %8d" + " | %s %s %s" * len(DEFS))
              % (kind, gr, n, *out))

print("\n### THE TEST: Mantel-Haenszel OR, stops vs walk_continues, stratified by")
print("### EXACT rank. OR > 1 means the walk_stops group carries the long-read end")
print("### more often AT THE SAME POSITION IN THE GENE. Not capped by the")
print("### baseline rate, which is why it replaces the ratio of two o/e values.")
for d in DEFS:
    print("  -- reference: %s" % d)
    print("     %-8s %-16s %8s %9s %-16s %-16s" %
          ("end", "stratum", "n", "MH OR", "95% CI", "Breslow-Day"))
    for kind in ENDS:
        for gr in STOPS:
            if not n_g[(kind, gr)]: continue
            OR, lo, hi, df, p = mh_or(n_rk, h_rk, kind, gr, "walk_continues", d)
            print("     %-8s %-16s %8d %9.2f  [%5.2f, %5.2f]   df=%-4d p=%.1e%s" %
                  (kind, gr, n_g[(kind, gr)], OR, lo, hi, df, p,
                   "  heterogeneous" if (p == p and p < 0.05) else ""))

print("\n### why the ratio of o/e values failed: headroom by rank (`exists`)")
print("  %-8s %6s %10s %10s %12s" % ("end", "rank", "bg rate", "max o/e", "n bg parts"))
for kind in ENDS:
    for rk in (0, 1, 2, 5, 10, 20):
        n = bgn[(kind, rk)]
        if not n: continue
        rate = bgh[(kind, "exists", rk)] / n
        print("  %-8s %6d %9.1f%% %10s %12d" %
              (kind, rk, 100*rate, ("%.2fx" % (1/rate)) if rate else "inf", n))
print("  the groups do not share a rank distribution, so the cap bites unevenly:")
for kind in ENDS:
    for gr in GRPS:
        n = n_g[(kind, gr)]
        if not n: continue
        lo = sum(v for (kk, gg, rk), v in n_rk.items() if kk == kind and gg == gr and rk <= 1)
        print("  %-8s %-16s %8d  ranks 0-1: %5.1f%%" % (kind, gr, n, 100*lo/n))

## The per-stratum DENSITY table (median unit_mean / unit_length for terminus
## and flank) was removed with the depth criterion on 2026-10-23. What it showed
## still holds and is the reason the two strata must stay split: walk_continues
## needs a flank ABOVE the read floor and unpowered stops needs one BELOW it, so
## the strata are selected in opposite directions on coverage. That is a
## statement about selection, not a criterion, and the stratum sizes below carry
## it without computing a density.

def q(v, p):
    v = sorted(v); return v[min(len(v)-1, int(p*len(v)))]
print("\n### signed distance to the NEAREST long-read `exists` end (+ = outward):")
print("### separates ABSENCE of a boundary from DISPLACEMENT of one")
print("  %-8s %-16s %8s %8s %8s %8s %8s %8s %8s" %
      ("end", "stratum", "n", "p10", "median", "p90", "|d|<100", "|d|<1kb", "outward"))
for kind in ENDS:
    for gr in GRPS:
        v = sig[(kind, gr)]
        if gr == "walk_stops" or not v: continue
        print("  %-8s %-16s %8d %8d %8d %8d %7.1f%% %7.1f%% %7.1f%%" %
              (kind, gr, len(v), q(v, .10), q(v, .50), q(v, .90),
               100*sum(abs(x) < 100 for x in v)/len(v),
               100*sum(abs(x) < 1000 for x in v)/len(v),
               100*sum(x > 0 for x in v)/len(v)))

print("\n### rank composition (0 = the gene's outermost part at that end)")
print("  %-8s %-16s %8s %8s %8s %9s" % ("end", "stratum", "n", "median", "p90", "rank0 %"))
for kind in ENDS:
    for gr in GRPS:
        v = rk_of[(kind, gr)]
        if gr == "walk_stops" or not v: continue
        print("  %-8s %-16s %8d %8d %8d %8.1f%%" %
              (kind, gr, len(v), q(v, .50), q(v, .90), 100*sum(x == 0 for x in v)/len(v)))

print("\n### which term_verdict feeds each stops stratum")
for kind in ENDS:
    for gr in ("stops_measured", "stops_unpowered"):
        tot = n_g[(kind, gr)]
        if not tot: continue
        for tv in sorted(v for (kk, gg, v) in vmix if kk == kind and gg == gr):
            n = vmix[(kind, gr, tv)]
            print("  %-8s %-16s %-32s %8d  %5.1f%%" % (kind, gr, tv, n, 100*n/tot))
