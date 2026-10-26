#!/usr/bin/env python3
"""Do GrASE terminal CALLS coincide with long-read transcript ends more often
than terminal sides it did not call, at matched position?

WHY THIS EXISTS. Every other long-read comparison here is between two walk
verdicts INSIDE already-called sides -- walk_stops against walk_continues -- so
it tests the walk's verdict, not the call, which both groups share. This is the
first comparison against sides GrASE did not call.

IT BUILDS ITS OWN GEOMETRY, and must. outward_walk_longread.py's `plen` is built
over `strand_of`, which it fills from the side table -- called sides only, 7,585
genes. Reusing it silently drops every uncalled terminus in the other ~12,400
genes and leaves a control pool drawn entirely from genes that have calls. The
per-gene DEXSeq GFFs cover 60,669 genes, so the coordinates are available; this
reads them for every gene in the units table.

THE POSITIVE IS part_called, THE TERMINUS UNIT ITSELF. Not "a terminus of a
called side": a called side contributes its whole ordered run of candidate
starts t_1 .. t_n, and the walk exists precisely to rule most of them out, so
treating all of them as called inflates the positive set with the termini the
method rejects. part_called comes from outward_walk_verdict.tsv, where the
terminus was tested as its own unit under the call condition. A terminus counts
as positive if part_called is TRUE in at least one contrast, matching the
contrast-free reference below.

The control pool is unchanged: termini that no called side reaches, from the
all-sides units table, deduplicated per (gene, anchor, end).

EXPRESSION IS CONTROLLED PER TERMINUS, ON SHORT READS. A terminal call can only
be made where there is coverage, so called termini are better covered than a
random uncalled terminus; matching on position alone leaves that in, and it is
large -- this comparison gave `exists` 4.30 with no expression control at all.
The control variable is the TERMINUS'S OWN mean DEXSeq count, max over the
contrast's two lineages, from part_shortread_expr.tsv. A terminus is an exonic
part and DEXSeq counts every part whether GrASE called it or not, so a control
terminus has a count exactly as a called one does.
Two measures were used before this and are both worse. IsoQuant gene totals come
from the same models that form the outcome, sit downstream of whether the
terminus is really used, and IsoQuant is not an established quantifier. A
gene-level mean of baseMean over the gene's tested sides is dominated by the
shared ref_ex_part, which appears in many rows, and varies with how many
bipartitions a gene has rather than with its expression.

THE CONTROL IS SAMPLED, NOT POOLED. For each called terminus a control is drawn
from the uncalled set in the SAME (rank bin, d_near bin) cell, so the two match
exactly on both positional variables. The pooled Mantel-Haenszel alternative
assumes a common odds ratio and Breslow-Day rejects that throughout this data; a
matched pair needs no such assumption. The test is McNemar on discordant pairs.

TWO CONTROL DEFINITIONS, because "uncalled" is ambiguous:
  any     any terminus no called side reaches
  clean   the same, restricted to genes with NO called terminal side at all, so
          a real alternative TSS elsewhere in the gene cannot put a model end
          near the control by coincidence

REFERENCES. `exists` is unioned over all 6 PAIRS, covering all four lineages, so
it is contrast-free: is there a transcript end here at all. `share` and
`exclusive` are properties of the long-read data for a lineage pair, not of the
GrASE call, so they are reported by outward_effect_size.py, which scores each
control under its case's own contrast.

Geometry from outward_units.allsides.tsv (enumerate_outward_units.R --all_sides).
Never point build_outward_counts.R at that file.

Usage: IG=1 python3 -I outward_called_vs_uncalled.py [seed]
"""
import os, re, sys, csv, math, random, collections, contextlib, io

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 1
MAIN = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "outward_walk_longread.py")
sys.argv = [MAIN]
G = {"__name__": "_walk_side", "__file__": MAIN}
sys.stderr.write("reusing %s (its stdout discarded)\n" % MAIN)
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(open(MAIN).read(), MAIN, "exec"), G)
for nm in ("build", "near", "hit", "PAIRS", "DEFS", "rbin", "dbin", "D", "PAN"):
    globals()[nm] = G[nm]

UNITS = os.path.join(PAN, "outward_units.allsides.tsv")

## ---- pass 1: geometry of every terminus, and which SIDES reach it ---------
tstrand = {}; kind_of = {}
for r in csv.DictReader(open(UNITS), delimiter="\t"):
    if r["anchor_kind"] != "terminus_self": continue
    anc = "%03d" % int(r["anchor"])
    k = (r["gene"], anc, r["kind"])
    tstrand[k] = r["strand"]
    ## (gene,event,side,anchor) -> end, to resolve the verdict table's rows
    kind_of[(r["gene"], r["event"], r["side"], anc)] = r["kind"]

## ---- pass 2: the positive set is part_called, per terminus ---------------
pos = set()
VERD = os.path.join(PAN, "outward_walk_verdict.tsv")
for r in csv.DictReader(open(VERD), delimiter="\t"):
    if r["part_called"] != "TRUE": continue
    anc = "%03d" % int(r["anchor"])
    kd = kind_of.get((r["gene"], r["event"], r["side"], anc))
    if kd is None: continue
    pos.add((r["gene"], anc, kd))

## a terminus is a control only if NO called side reaches it, and it is not
## itself a positive under any contrast
reached = set()
for r in csv.DictReader(open(UNITS), delimiter="\t"):
    if r["anchor_kind"] != "terminus_self" or r["called"] != "TRUE": continue
    reached.add((r["gene"], "%03d" % int(r["anchor"]), r["kind"]))
tcall = {}
for k in tstrand:
    if k in pos: tcall[k] = True
    elif k in reached: continue          # candidate start of a called side: neither
    else: tcall[k] = False
sys.stderr.write("positives (part_called) %d | controls %d | excluded candidate "
                 "starts of called sides %d\n" %
                 (sum(1 for v in tcall.values() if v),
                  sum(1 for v in tcall.values() if not v),
                  len(reached - pos)))
genes = {k[0] for k in tcall}
sys.stderr.write("termini %d over %d genes\n" % (len(tcall), len(genes)))

## ---- full geometry, every gene, from the per-gene DEXSeq GFFs -------------
gplen = {}; gchrom = {}; gnpart = {}
RX = re.compile(r'exonic_part_number "(\d+)"')
for g in genes:
    gf = "%s/dexseq.gff/%s.dexseq.gff" % (D, g)
    if not os.path.exists(gf): continue
    d = {}
    for ln in open(gf):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exonic_part": continue
        m = RX.search(f[8])
        if not m: continue
        d[m.group(1)] = (int(f[3]), int(f[4]))
        gchrom[g] = f[0]
    if d: gplen[g] = d; gnpart[g] = len(d)
sys.stderr.write("geometry for %d of %d genes\n" % (len(gplen), len(genes)))

def place(g, anc, kind, strand):
    """(pos, cell) -- the same construction outward_walk_longread.py uses."""
    low = (kind == "5prime") == (strand == "+")
    s0, e0 = gplen[g][anc]
    pos = s0 if low else e0
    rk = (int(anc) - 1) if low else (gnpart[g] - int(anc))
    lo = min(a for a, b in gplen[g].values()); hi = max(b for a, b in gplen[g].values())
    dn = max((pos - lo) if low else (hi - pos), 0)
    return pos, (rbin(max(rk, 0)), dbin(dn))

## ---- per-terminus short-read coverage, the control variable --------------
PANEL = {"%s_vs_%s" % (t, r): (t, r) for t, r in PAIRS}
sx = {}
for r in csv.DictReader(open(os.path.join(PAN, "part_shortread_expr.tsv")),
                        delimiter="\t"):
    sx[(r["gene"], r["part"])] = r
## exact zero is its own bin: no reads at all is not the same as a few
EXL = ["0", "<0.5", "0.5-2", "2-10", "10-50", "50-200", ">=200"]
EXB = [0.5, 2.0, 10.0, 50.0, 200.0]
def xbin(g, anc, c):
    r = sx.get((g, anc))
    if r is None: return None
    t, rf = PANEL[c]
    ## max over the two lineages, matching add_support()'s d_support convention:
    ## a pooled mean dilutes a terminus on in one condition and off in the other
    v = max(float(r[t]), float(r[rf]))
    if v == 0: return 0
    for i, b in enumerate(EXB):
        if v < b: return i + 1
    return len(EXL) - 1
sys.stderr.write("per-part coverage for %d (gene, part) pairs\n" % len(sx))

## ---- positives: part_called, PER CONTRAST, carrying the claimed direction --
## The reference is directional, so the unit must be (terminus, contrast): the
## same terminus called in two contrasts makes two claims, with their own signs
## and their own reference. A union over contrasts has no direction and cannot
## be scored against a directional reference, which is what this script did
## until 2026-10-26.
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
    cases[c].append((r["gene"], anc, kd, tstrand.get((r["gene"], anc, kd)), dpi))

## ---- controls: termini no called side reaches ----------------------------
reached = set()
for r in csv.DictReader(open(UNITS), delimiter="\t"):
    if r["anchor_kind"] != "terminus_self" or r["called"] != "TRUE": continue
    reached.add((r["gene"], "%03d" % int(r["anchor"]), r["kind"]))
ctl = [k for k in tstrand if k not in reached]
sys.stderr.write("cases %d over %d contrasts | control termini %d\n" %
                 (sum(len(v) for v in cases.values()), len(cases), len(ctl)))

def mcnemar(b, c):
    """chi2 with continuity correction. erfc underflows above chi2 ~1500, so
    chi2 is printed and a p of 0 means below double precision."""
    if b + c == 0: return float("nan"), float("nan")
    chi = (abs(b - c) - 1) ** 2 / (b + c)
    return chi, math.erfc(math.sqrt(chi / 2.0))

rng = random.Random(SEED)
acc = collections.Counter(); nomatch = collections.Counter()
for c, (trt, ref) in PANEL.items():
    ds = build(trt, ref)
    ## control pool for THIS contrast: matched on end, rank, distance, expression
    pool = collections.defaultdict(list)
    for k in ctl:
        g, anc, kd = k
        if g not in gplen or anc not in gplen[g]: continue
        xb = xbin(g, anc, c)
        if xb is None: continue
        pos, cell = place(g, anc, kd, tstrand[k])
        pool[(kd, cell, xb)].append((g, anc, tstrand[k], pos))
    for (g, anc, kd, st, dpi) in cases.get(c, ()):
        if st is None or g not in gplen or anc not in gplen[g]: continue
        xb = xbin(g, anc, c)
        if xb is None: nomatch[(kd, "no coverage row")] += 1; continue
        pos, cell = place(g, anc, kd, st)
        cand = pool.get((kd, cell, xb), ())
        if not cand: nomatch[(kd, "no control in cell")] += 1; continue
        cg, canc, cst, cpos = cand[rng.randrange(len(cand))]
        e = "5" if kd == "5prime" else "3"
        acc[(kd, "n")] += 1
        for d in DEFS:
            ## the CONTROL is scored under the case's contrast AND the case's
            ## claimed direction: the pair must differ only in being called.
            h1 = hit(ds, d, e, gchrom[g], st, pos, dpi)
            h0 = hit(ds, d, e, gchrom[cg], cst, cpos, dpi)
            if h1: acc[(kd, d, "case")] += 1
            if h0: acc[(kd, d, "ctl")] += 1
            if h1 and not h0: acc[(kd, d, "b")] += 1
            if h0 and not h1: acc[(kd, d, "c")] += 1

print("dropped cases: " + (", ".join("%s/%s %d" % (k[0], k[1], v)
      for k, v in sorted(nomatch.items())) if nomatch else "0"))
print("\n### called termini vs uncalled, matched on end, rank, d_near and the")
print("### TERMINUS'S OWN short-read coverage; references resolved in the")
print("### CLAIMED direction")
print("  %-8s %-10s %9s %9s %9s %7s %7s %8s %10s" %
      ("end", "reference", "pairs", "called", "control", "b", "c", "OR", "chi2"))
for kind in ("5prime", "3prime"):
    n = acc[(kind, "n")]
    if not n: continue
    for d in DEFS:
        b = acc[(kind, d, "b")]; cc = acc[(kind, d, "c")]
        chi, pv = mcnemar(b, cc)
        print("  %-8s %-10s %9d %8.1f%% %8.1f%% %7d %7d %8.2f %10.0f" %
              (kind, d, n, 100*acc[(kind,d,"case")]/n, 100*acc[(kind,d,"ctl")]/n,
               b, cc, (b/cc) if cc else float("inf"), chi))
