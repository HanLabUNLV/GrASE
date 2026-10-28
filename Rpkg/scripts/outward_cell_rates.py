#!/usr/bin/env python3
"""Long-read support of walk_stops against walk_continues INSIDE each matched
cell, rather than pooled.

WHY THIS EXISTS. The Mantel-Haenszel OR in outward_walk_longread.py is a single
pooled summary, and Breslow-Day rejects for every one of them, so each OR is an
average over strata that genuinely differ. A pooled average cannot distinguish
an effect present in most strata from a thin majority of strata, and those are
very different claims about the walk. This prints the rate in each cell so the
distinction is visible, and tests it with a sign test over cells, which assumes
nothing about a common OR.

The cell is (rank bin, d_near bin) -- the same cell the expectation in
outward_walk_longread.py is matched on, so this table and that o/e table are on
the same footing. rbin, dbin and dnear are imported from that script, not
redefined.

WHY THE o/e TABLE MUST NOT BE READ AS A GROUP COMPARISON, which this makes
concrete. o/e weights each cell by its expected count, n times background.
walk_stops sits disproportionately in near-boundary cells where background is
saturated -- at 5' rbin 0 / d_near 0-99 the background is 71.4% and walk_stops
scores 67.9%, a ratio below 1 -- and with several hundred termini that cell
carries a large share of the weight. walk_continues sits further out where
background is a few percent and ratios are large. So the two groups are scored
against different baselines under a weighting that penalises whichever group is
nearer the boundary, which is why their o/e values can be equal while the
within-cell rates differ by a factor of two.

MIN_N (default 40) is the smallest cell reported and counted. Smaller cells are
dominated by single termini; the tally is printed at two thresholds so the
choice is visible rather than load-bearing.

Usage: IG=1 python3 -I outward_cell_rates.py [min_n]
"""
import os, sys, csv, math, collections, contextlib, io

MIN_N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
MAIN = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "outward_walk_longread.py")
sys.argv = [MAIN]
G = {"__name__": "_walk_side", "__file__": MAIN}
sys.stderr.write("reusing %s (its stdout discarded)\n" % MAIN)
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(open(MAIN).read(), MAIN, "exec"), G)
for nm in ("build", "near", "verd", "meta", "plen", "npart", "have", "chrom_of",
           "strand_of", "DEFS", "PAIRS", "IG", "rbin", "dbin", "dnear",
           "bg_n", "bg_h", "mh_or", "hit", "PAN"):
    globals()[nm] = G[nm]

ENDS = ("5prime", "3prime")
DBL  = ["0-99", "100-299", "300-999", "1k-3k", "3k-10k",
        "10k-30k", "30k-100k", "100k+"]

def binom_two_sided(k, n):
    """Exact two-sided sign test against p=0.5. No scipy in this tree."""
    if n == 0: return float("nan")
    c = [math.comb(n, i) for i in range(n + 1)]
    tot = float(sum(c))
    tgt = c[k]
    return min(1.0, sum(x for x in c if x <= tgt * (1 + 1e-12)) / tot)

## GENE-LEVEL LONG-READ COVERAGE IS A THIRD CONFOUND, and it is not optional.
## In outward_called_vs_uncalled.py, controlling it moved the called-vs-uncalled
## `exists` odds ratio from 4.30 to 1.58 and erased `share` and `exclusive`
## entirely. walk_stops and walk_continues are both inside called genes so they
## are more comparable than that comparison was, but "more comparable" is not
## "controlled". These counters add an EXPRESSION band to the cell and restrict
## the references to expressed genes.
##
## EXPRESSION IS MEASURED ON SHORT READS, PER TERMINUS. IsoQuant gene totals
## were used here until 2026-10-26; they come from the same models that form the
## outcome, sit downstream of whether the terminus is really used, and IsoQuant
## is not an established quantifier.
## PER TERMINUS, not per gene. A terminus is an exonic part and DEXSeq counts
## every part, so the band is the coverage of the exact feature under test. The
## gene-level mean of baseMean used here until 2026-10-26 is dominated by the
## shared ref_ex_part and varies with a gene's bipartition count; in the
## called-vs-uncalled comparison it selected controls BETTER covered at the
## terminus than the cases and masked the effect entirely.
_LIN = {"%s_vs_%s" % (t, r): (t, r) for t, r in PAIRS}
_sx = {}
for _r in csv.DictReader(open(os.path.join(PAN, "part_shortread_expr.tsv")),
                         delimiter="\t"):
    _sx[(_r["gene"], _r["part"])] = _r
GBL = ["0", "<0.5", "0.5-2", "2-10", "10-50", "50-200", ">=200"]
_EB = [0.5, 2.0, 10.0, 50.0, 200.0]
def _cov(g, anc, c):
    r = _sx.get((g, anc))
    if r is None: return None
    t, rf = _LIN[c]
    return max(float(r[t]), float(r[rf]))
def gband(g, anc, c):
    v = _cov(g, anc, c)
    if v is None: return "NA"
    if v == 0: return GBL[0]
    for i, b in enumerate(_EB):
        if v < b: return GBL[i + 1]
    return GBL[-1]
def evaluable(g, anc, c):
    v = _cov(g, anc, c)
    return v is not None and v > 0
nk3 = collections.Counter(); hk3 = collections.Counter()   # cell x expression band
nkE = collections.Counter(); hkE = collections.Counter()   # expressed genes only
n_c = collections.Counter(); h_c = collections.Counter()
for trt, ref in PAIRS:
    ds = build(trt, ref)
    for r in verd.get("%s_vs_%s" % (trt, ref), []):
        g = r["gene"]
        if g not in plen: continue
        anc = "%03d" % int(r["anchor"])
        k = (g, r["event"], r["side"], anc)
        if k not in have or k not in meta: continue
        co = plen[g].get(anc)
        if not co: continue
        kind, strand = meta[k]
        e = "5" if kind == "5prime" else "3"
        low = (kind == "5prime") == (strand == "+")
        pos = co[0] if low else co[1]
        rk = (int(anc) - 1) if low else (npart[g] - int(anc))
        grp = r["walk"]
        if grp == "no_claim": continue
        cc = "%s_vs_%s" % (trt, ref)
        c = (rbin(max(rk, 0)), dbin(dnear(g, pos, low)))
        c3 = c + (gband(g, anc, cc),)
        ev = evaluable(g, anc, cc)
        n_c[(kind, grp, c)] += 1
        nk3[(kind, grp, c3)] += 1
        if ev: nkE[(kind, grp, c)] += 1
        try: dpi = float(r["part_dpi"])
        except (KeyError, ValueError): dpi = None
        for d in DEFS:
            if hit(ds, d, e, chrom_of[g], strand, pos, dpi):
                h_c[(kind, grp, d, c)] += 1
                hk3[(kind, grp, d, c3)] += 1
                if ev: hkE[(kind, grp, d, c)] += 1

def tally(kind, d, min_n):
    hi = lo = 0; ns = nc = 0
    for c in sorted({cc for (kk, gg, cc) in n_c if kk == kind}):
        a = n_c[(kind, "walk_stops", c)]; b = n_c[(kind, "walk_continues", c)]
        if a < min_n or b < min_n: continue
        ra = h_c[(kind, "walk_stops", d, c)] / a
        rb = h_c[(kind, "walk_continues", d, c)] / b
        if ra > rb: hi += 1
        else: lo += 1
        ns += a; nc += b
    return hi, lo, ns, nc

for kind in ENDS:
    print("\n### %s -- `exists` rate inside each (rank bin, d_near bin) cell, "
          "cells with >=%d in both" % (kind, MIN_N))
    print("  %-6s %-10s %8s %8s %9s %9s %9s  %s" %
          ("rbin", "d_near", "n stops", "n cont", "stops", "cont", "bg", "stops higher"))
    for c in sorted({cc for (kk, gg, cc) in n_c if kk == kind}):
        a = n_c[(kind, "walk_stops", c)]; b = n_c[(kind, "walk_continues", c)]
        if a < MIN_N or b < MIN_N: continue
        ra = h_c[(kind, "walk_stops", "exists", c)] / a
        rb = h_c[(kind, "walk_continues", "exists", c)] / b
        bg = bg_h[(kind, "exists", c)] / bg_n[(kind, c)] if bg_n[(kind, c)] else float("nan")
        print("  %-6s %-10s %8d %8d %8.1f%% %8.1f%% %8.1f%%  %s" %
              (c[0], DBL[c[1]], a, b, 100*ra, 100*rb, 100*bg,
               "yes" if ra > rb else "NO"))

print("\n### SIGN TEST over cells. Counts cells, not termini, so no common odds")
print("### ratio is assumed -- which matters because Breslow-Day rejects.")
print("  %-8s %-10s %6s %6s %6s %10s %9s %9s" %
      ("end", "reference", "min n", "stops", "cont", "two-sided p", "n stops", "n cont"))
for kind in ENDS:
    for d in DEFS:
        for mn in (MIN_N, 100):
            hi, lo, ns, nc = tally(kind, d, mn)
            print("  %-8s %-10s %6d %6d %6d %10.4f %9d %9d" %
                  (kind, d, mn, hi, lo, binom_two_sided(hi, hi + lo), ns, nc))


## ---- the same comparison with GENE-LEVEL LONG-READ COVERAGE controlled ----
print("\n### walk_stops vs walk_continues, MH OR, coverage controlled")
print("### `as shipped` stratifies on exact rank only; the others add d_near and")
print("### then the short-read expression band, or restrict to expressed genes.")
print("  %-8s %-10s %-26s %-26s %-26s" %
      ("end", "reference", "rank x d_near", "+ expr band", "expressed genes only"))
for kind in ("5prime", "3prime"):
    for d in DEFS:
        out = []
        for nk, hk in ((n_c, h_c), (nk3, hk3), (nkE, hkE)):
            OR, lo, hi, df, p = mh_or(nk, hk, kind, "walk_stops", "walk_continues", d)
            out.append("%5.2f [%5.2f, %5.2f]" % (OR, lo, hi)
                       if OR == OR else "            n/a")
        print("  %-8s %-10s %-26s %-26s %-26s" % (kind, d, out[0], out[1], out[2]))

print("\n### and by short-read expression band alone")
print("  %-8s %-10s %s" % ("end", "reference", "".join("%16s" % b for b in GBL)))
for kind in ("5prime", "3prime"):
    for d in DEFS:
        row = []
        for gb in GBL:
            sub = collections.Counter(); subh = collections.Counter()
            for (kk, gg, c3), v in nk3.items():
                if kk == kind and c3[2] == gb: sub[(kk, gg, c3[:2])] += v
            for (kk, gg, dd, c3), v in hk3.items():
                if kk == kind and dd == d and c3[2] == gb: subh[(kk, gg, dd, c3[:2])] += v
            n = sum(v for (kk, gg, _), v in sub.items() if kk == kind)
            OR, lo, hi, df, p = mh_or(sub, subh, kind, "walk_stops", "walk_continues", d)
            row.append("  n/a" if n < 200 or OR != OR else "%5.2f (%d)" % (OR, n))
        print("  %-8s %-10s %s" % (kind, d, "".join("%16s" % x for x in row)))
