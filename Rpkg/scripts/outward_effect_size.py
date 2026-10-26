#!/usr/bin/env python3
"""Does long-read support of a GrASE terminal call rise with its effect size,
and does the DIFFERENTIAL reference beat a matched control?

TWO QUESTIONS left open by outward_called_vs_uncalled.py, which pooled all calls
and used `exists` only.

(1) If the calls are real, the largest-|delta_pi| calls should separate from
    their matched controls more than the smallest. A flat profile would mean the
    enrichment is not tracking what the test measures.

(2) `share` and `exclusive` are properties of the LONG-READ data at a position
    for a lineage pair, not of the GrASE call, so a control CAN be scored on
    them -- under the contrast of the case it is matched to. Case and control
    then differ only in being called, and `exclusive` is the reference that
    speaks to what a terminal call actually asserts.

GEOMETRY IS REUSED, NOT REBUILT. It execs outward_called_vs_uncalled.py with
stdout discarded and takes gplen/gchrom/gnpart, the per-terminus call status and
place() from it. outward_walk_longread.py's own `plen` covers called genes only
and must not be used here; that script's header explains why.

UNIT is (gene, terminus, end, contrast): the same terminus called in two
contrasts is two tests, because the reference differs between them.

Restricted to the 6 contrasts in PAIRS; the long-read panel has no NK.
`significant` is exontest's side-level flag, NOT the scorer's call condition --
it admits |delta_pi| below 0.1 -- so the bins start at 0.

Usage: IG=1 python3 -I outward_effect_size.py [seed]
"""
import os, sys, csv, math, random, collections, contextlib, io

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 1
CVU = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "outward_called_vs_uncalled.py")
sys.argv = [CVU]
G = {"__name__": "_cvu_side", "__file__": CVU}
sys.stderr.write("reusing %s (its stdout discarded)\n" % CVU)
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(open(CVU).read(), CVU, "exec"), G)
for nm in ("near", "hit", "build", "PAIRS", "DEFS", "gplen", "gchrom",
           "gnpart", "place", "tstrand", "ctl", "PAN", "UNITS", "sx", "xbin",
           "EXL"):
    globals()[nm] = G[nm]

PANEL = {"%s_vs_%s" % (t, r): (t, r) for t, r in PAIRS}

## EXPRESSION IS CONTROLLED ON SHORT READS, and effect size is read at fixed
## expression. IsoQuant gene totals were used here until 2026-10-26: they come
## from the same models that form the outcome, sit downstream of whether the
## terminus is really used, and IsoQuant is not an established quantifier.
## The confound they were covering is real and worth recording: as |delta_pi|
## rises the panel's coverage of the gene falls while SHORT-READ d_support stays
## flat (mean 124 to 150 reads across the bins). So GrASE's own evidence does
## not weaken with effect size, only the panel's ability to check it -- which
## depresses `exists` and inflates `exclusive`, since exclusive needs b == 0 and
## sparse data supplies that for free.
## bands are the TERMINUS'S OWN DEXSeq coverage, from outward_called_vs_uncalled
GBL = EXL
def gbin(g, anc, c):
    i = xbin(g, anc, c)
    return "NA" if i is None else GBL[i]
EB = [0.0, 0.1, 0.2, 0.3, 0.5, 1.01]
EBL = ["<0.1", "0.1-0.2", "0.2-0.3", "0.3-0.5", ">=0.5"]
def ebin(d):
    for i in range(len(EB) - 1):
        if EB[i] <= d < EB[i + 1]: return i
    return len(EBL) - 1

## (gene,event,side,anchor) -> end, to resolve verdict rows
kind_of = {}
for r in csv.DictReader(open(UNITS), delimiter="\t"):
    if r["anchor_kind"] != "terminus_self": continue
    kind_of[(r["gene"], r["event"], r["side"], "%03d" % int(r["anchor"]))] = r["kind"]

## control pool, per contrast: termini no called side reaches, matched on end,
## rank, d_near AND short-read expression
def build_pool(c):
    pl = collections.defaultdict(list)
    for k in ctl:
        g, anc, kind = k
        if g not in gplen or anc not in gplen[g]: continue
        xb = xbin(g, anc, c)
        if xb is None: continue
        pos, cell = place(g, anc, kind, tstrand[k])
        pl[(kind, cell, xb)].append((g, anc, tstrand[k], pos))
    return pl

## THE POSITIVE IS part_called AND THE EFFECT SIZE IS part_dpi, both from the
## verdict table. Until 2026-10-26 this read side-level significance from
## called_effect.tsv and expanded each called side to ALL of its termini, which
## is the inflated set outward_called_vs_uncalled.py abandoned: it includes the
## candidate starts the walk exists to reject, and they sit at gene-outermost
## positions where long-read ends are common. It also binned on the SIDE's
## delta_pi, which is not the terminus's effect size. With both fixed this is a
## stratification of that script's table, not a separate analysis.
cases = collections.defaultdict(list)
seen = set()
for r in csv.DictReader(open(os.path.join(PAN, "outward_walk_verdict.tsv")),
                        delimiter="\t"):
    if r["part_called"] != "TRUE": continue
    c = r["contrast"]
    if c not in PANEL: continue
    anc = "%03d" % int(r["anchor"])
    kd = kind_of.get((r["gene"], r["event"], r["side"], anc))
    if kd is None: continue
    k = (c, r["gene"], anc, kd)
    if k in seen: continue
    seen.add(k)
    try: dpi = float(r["part_dpi"])
    except ValueError: continue
    st = tstrand.get((r["gene"], anc, kd))
    if st is None: continue
    cases[c].append((r["gene"], anc, kd, st, ebin(abs(dpi)), dpi))

def mcnemar(b, c):
    if b + c == 0: return float("nan"), float("nan")
    chi = (abs(b - c) - 1) ** 2 / (b + c)
    return chi, math.erfc(math.sqrt(chi / 2.0))

rng = random.Random(SEED)
acc = collections.Counter(); nomatch = collections.Counter()
for c, (trt, ref) in PANEL.items():
    ds = build(trt, ref)
    pool = build_pool(c)
    for (g, anc, kind, strand, eb, dpi) in cases.get(c, ()):
        if g not in gplen or anc not in gplen[g]: continue
        pos, cell = place(g, anc, kind, strand)
        xb = xbin(g, anc, c)
        if xb is None:
            nomatch[(kind, "no coverage row")] += 1; continue
        cand = pool.get((kind, cell, xb), ())
        if not cand:
            nomatch[(kind, "no control in cell")] += 1; continue
        cg, canc, cst, cpos = cand[rng.randrange(len(cand))]
        e = "5" if kind == "5prime" else "3"
        acc[(kind, eb, "n")] += 1
        gb = gbin(g, anc, c); acc[(kind, "G", gb, eb, "n")] += 1
        for d in DEFS:
            ## both scored in the CASE's claimed direction
            h1 = hit(ds, d, e, gchrom[g], strand, pos, dpi)
            h0 = hit(ds, d, e, gchrom[cg], cst, cpos, dpi)
            if h1: acc[(kind, eb, d, "case")] += 1
            if h0: acc[(kind, eb, d, "ctl")] += 1
            if h1 and not h0:
                acc[(kind, eb, d, "b")] += 1; acc[(kind, "G", gb, eb, d, "b")] += 1
            if h0 and not h1:
                acc[(kind, eb, d, "c")] += 1; acc[(kind, "G", gb, eb, d, "c")] += 1

print("dropped cases: " + (", ".join("%s/%s %d" % (k[0], k[1], v)
      for k, v in sorted(nomatch.items())) if nomatch else "0"))
print("Matched on end, rank, d_near and the terminus's own short-read\n"
      "coverage; references resolved in the case's claimed direction.")
for kind in ("5prime", "3prime"):
    for d in DEFS:
        print("\n### %s -- `%s`, called against matched uncalled (clean) control" % (kind, d))
        print("  %-9s %8s %9s %9s %7s %7s %8s %10s" %
              ("|delta_pi|", "pairs", "called", "control", "b", "c", "OR", "chi2"))
        for eb, lab in enumerate(EBL):
            n = acc[(kind, eb, "n")]
            if not n: continue
            b = acc[(kind, eb, d, "b")]; cc = acc[(kind, eb, d, "c")]
            chi, p = mcnemar(b, cc)
            print("  %-9s %8d %8.1f%% %8.1f%% %7d %7d %8.2f %10.0f" %
                  (lab, n, 100*acc[(kind,eb,d,"case")]/n,
                   100*acc[(kind,eb,d,"ctl")]/n, b, cc,
                   (b/cc) if cc else float("inf"), chi))


## the effect-size trend at FIXED long-read coverage
for kind in ("5prime", "3prime"):
    for d in DEFS:
        print("\n### %s -- `%s`, OR by |delta_pi| WITHIN long-read coverage band" % (kind, d))
        print("  %-10s %s" % ("expr band", "".join("%14s" % e for e in EBL)))
        for gb in GBL:
            row = []
            for eb in range(len(EBL)):
                n = acc[(kind, "G", gb, eb, "n")]
                b = acc[(kind, "G", gb, eb, d, "b")]; cc = acc[(kind, "G", gb, eb, d, "c")]
                row.append("  n/a" if n < 100 else
                           ("%5.2f (%d)" % ((b/cc) if cc else float("inf"), n)))
            print("  %-10s %s" % (gb, "".join("%14s" % x for x in row)))
