#!/usr/bin/env python3
"""The outward walk's groups, scored against the IsoQuant long reads.

UNIT is (terminus, contrast). The walk partitions called termini three ways:

  walk_stops      part_called AND NOT flank_called
  walk_continues  part_called AND flank_called
  no_claim        NOT part_called

and each is scored against three long-read references, all evaluated as "an
IsoQuant transcript model whose matching end lies within +/-W bp of the
terminus", differing in what the model must do:

  exists      a + b > 0                       a transcript ENDS here
  exclusive   (a >= MINPRES and b == 0) or (b >= MINPRES and a == 0)

`exists` answers an existence question and `exclusive` a differential one.

THE `share` REFERENCE WAS REMOVED. It was
  |a/gtot_trt - b/gtot_ref| >= DSHARE, gene >= MINGENE in both
which uses IsoQuant counts as a QUANTIFIER -- it normalises by a gene total and
thresholds the difference of two ratios. IsoQuant is not established as an
accurate quantifier, so nothing here may depend on its magnitudes. `exists` and
`exclusive` use it only to detect a model end, with MINPRES as a presence floor,
which is the use it can carry. gtot is gone from this tree entirely; anything
needing an expression control uses the terminus's own DEXSeq count from
part_shortread_expr.tsv instead.

Rates are reported against a rank-matched within-gene background (a
random exonic part of the same gene at the same rank), and group-vs-group
comparison is a Mantel-Haenszel OR stratified by EXACT rank -- never a ratio
of two obs/exp values, which is invalid because obs/exp is a rate ratio whose
headroom is 1/background and the groups do not share a rank distribution.

THE DEPTH CRITERION (C1) WAS REMOVED 2026-10-23. It was
  density(outward flank)/density(terminus) < RATIO
a boundary-EXISTENCE claim, and it is no longer computed anywhere: it fired on
only 37.0% of IsoQuant-confirmed 5' boundaries and 15.3% at 3' exonic flanks,
it was a threshold rather than a test, and GrASE claims differential usage of
an annotated terminus, not that the terminus is real. The retired analyses are
under dead_scripts/ with the figures they produced. What remains here is the
walk's own discrimination, which is a differential criterion throughout.

Gene edge: with the intergenic flanks built (IG=1) a gene-edge terminus HAS an
outward flank, so those termini are INCLUDED. Without them the flank is absent
and the walk cannot judge the side.

Usage: IG=1 python3 -I outward_walk_longread.py [window] [minpres] [mode]
"""
import os, re, csv, math, bisect, collections, sys
from scipy.stats import chi2

W       = int(sys.argv[1])   if len(sys.argv) > 1 else 100
MINPRES = int(sys.argv[2])   if len(sys.argv) > 2 else 2
MODE    = sys.argv[3]        if len(sys.argv) > 3 else "pooled"
D = os.environ.get("DICE", "/mnt/data1/home/mirahan/DICE")
IG  = os.environ.get("IG", "0") == "1"
MAN = f"{D}/intronparts_ig/intronparts_manifest.tsv" if IG else \
      f"{D}/intronparts/intronparts_manifest.tsv"
PAN = f"{D}/intronparts/lineage/ig" if IG else f"{D}/intronparts/lineage"
sys.stderr.write("IG=%s  manifest=%s  panel=%s\n" % (IG, MAN, PAN))
TR = {"DRR465268":"CD4","DRR465269":"CD4","DRR465270":"CD4","DRR465271":"CD4",
 "DRR465272":"CD4","DRR465273":"CD4","DRR465274":"CD4","DRR465275":"CD4",
 "DRR465276":"CD4","DRR465277":"CD4","DRR465278":"CD4","DRR465279":"CD8",
 "DRR465280":"CD8","DRR465281":"CD8","DRR465295":"CD8","DRR465282":"B",
 "DRR465283":"B","DRR465284":"B","DRR465285":"B","DRR465288":"MONO.CLASSIC",
 "DRR465289":"MONO.CLASSIC","DRR465290":"MONO.CLASSIC","DRR465291":"MONO.CLASSIC"}
# MATCHED mode uses the ONE TRAILS subset that corresponds to each DICE group,
# so the long-read contrast is the same contrast DICE ran. Pooled mode uses all
# 11 CD4 subsets (Tregs, Th1/2/17, Tfh, memory ...), which asks a broader
# question than naive-CD4-vs-B and may be why C2+ anti-associates. The cost is
# depth: one library per lineage instead of up to eleven.
MATCHED = {"DRR465268":"CD4", "DRR465295":"CD8",
           "DRR465282":"B",   "DRR465288":"MONO.CLASSIC"}
if MODE == "matched": TR = MATCHED
LIN = ["CD4","CD8","B","MONO.CLASSIC"]
PAIRS = [("CD4","B"),("CD8","B"),("MONO.CLASSIC","B"),
         ("CD8","CD4"),("MONO.CLASSIC","CD4"),("MONO.CLASSIC","CD8")]
DEFS = ("exists","exclusive")          # the REPORTED references
## BOTH REFERENCES ARE DIRECTIONAL, so each index is built in two halves. A call
## says which lineage the terminus is used in -- the sign of part_dpi -- and only
## long-read evidence in THAT lineage corroborates it. Until 2026-10-26 both
## references ignored direction, which credits wrong-direction coincidences and
## dilutes a real effect toward the null; with roughly half the hits available in
## the wrong lineage that is a large dilution, not a rounding.
##
##   exists     the transcript is seen in the claimed lineage       a > 0
##   exclusive  and is wholly absent from the other       a >= MINPRES, b == 0
##
## They are a hierarchy, not alternatives: exclusive implies exists. The _trt
## and _ref halves are the same test with the lineages swapped, selected by the
## sign of dpi. The BACKGROUND has no call and so no direction; there the union
## of the two halves is used, which for exists reduces to the old a + b > 0.
IDX  = ("exists_trt","exists_ref","excl_trt","excl_ref")

O = f"{D}/longread/isoquant_TBM/OUT"
mcoord = {}
for ln in open(f"{O}/OUT.transcript_models.gtf"):
    f = ln.split("\t")
    if len(f) < 9 or f[2] != "transcript": continue
    t = re.search(r'transcript_id "([^"]+)"', f[8])
    if not t: continue
    c, s, e, st = f[0], int(f[3]), int(f[4]), f[6]
    mcoord[t.group(1)] = (c, st, s if st == "+" else e, e if st == "+" else s)
cnt = {}
with open(f"{O}/OUT.transcript_model_grouped_counts.tsv") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")[1:]
    lof = [TR.get(h) for h in hdr]
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        v = dict.fromkeys(LIN, 0.0)
        for val, l in zip(p[1:], lof):
            if l: v[l] += float(val)
        if sum(v.values()) > 0: cnt[p[0]] = v
tot = {l: sum(v[l] for v in cnt.values()) for l in LIN}
print(f"  mode {MODE} | libraries {len(TR)} | models with counts {len(cnt)} | window+/-{W}")
print("  lineage totals: " + "  ".join("%s %.0f" % (l, tot[l]) for l in LIN))

def build(trt, ref):
    out = {(d, e): collections.defaultdict(list) for d in IDX for e in ("5","3")}
    for m, v in cnt.items():
        co = mcoord.get(m)
        if not co: continue
        c, st, five, three = co
        a, b = v[trt], v[ref]
        ok = {"exists_trt": a > 0, "exists_ref": b > 0}
        ok["excl_trt"] = a >= MINPRES and b == 0
        ok["excl_ref"] = b >= MINPRES and a == 0
        for d in IDX:
            if ok[d]:
                out[(d,"5")][(c,st)].append(five); out[(d,"3")][(c,st)].append(three)
    for k in out:
        for kk in out[k]: out[k][kk] = sorted(set(out[k][kk]))
    return out

def hit(ds, d, e, c, st, pos, dpi):
    """One reference, resolved for direction. Both are directional: `exists`
    asks whether a transcript ends here IN THE LINEAGE THE CALL CLAIMS, and
    `exclusive` additionally that it is absent from the other. dpi is the
    terminus's part_dpi and its sign picks the lineage. With dpi None there is
    no claim to match -- an uncalled background part -- and the union of the two
    halves is used, which for exists is the old a + b > 0."""
    pre = "exists" if d == "exists" else "excl"
    if dpi is None:
        return (near(ds[(pre + "_trt", e)], c, st, pos) or
                near(ds[(pre + "_ref", e)], c, st, pos))
    if dpi == 0: return False
    return near(ds[(pre + ("_trt" if dpi > 0 else "_ref"), e)], c, st, pos)

def near(d, c, st, pos):
    v = d.get((c, st))
    if not v: return False
    i = bisect.bisect_left(v, pos - W)
    return i < len(v) and v[i] <= pos + W

rows = list(csv.DictReader(open(f"{PAN}/outward_side_table.tsv"), delimiter="\t"))
pair = collections.defaultdict(dict); meta = {}; strand_of = {}
for r in rows:
    ak = r["anchor_kind"]
    if ak not in ("terminus","terminus_self"): continue
    ## ONE CANONICAL ANCHOR FORM, `%03d`, matching the GFF's exonic_part_number
    ## and enumerate_outward_units.R's pid(). Normalise on read: `046` and `46`
    ## have both been on disk, and a mismatch makes every terminus fail the
    ## `k in have` guard silently rather than loudly.
    anc0 = "%03d" % int(r["anchor"])
    k = (r["gene"], r["event"], r["side"], anc0)
    pair[k][ak] = (r["unit_type"], r["unit_id"], float(r["unit_mean"]))
    meta[k] = (r["kind"], r["strand"]); strand_of.setdefault(r["gene"], r["strand"])
plen = {}; npart = {}; chrom_of = {}
for g in sorted(strand_of):
    gf = f"{D}/dexseq.gff/{g}.dexseq.gff"
    if not os.path.exists(gf): continue
    d = {}
    for ln in open(gf):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exonic_part": continue
        d[re.search(r'exonic_part_number "(\d+)"', f[8]).group(1)] = (int(f[3]), int(f[4]))
        chrom_of[g] = f[0]
    if d: plen[g] = d; npart[g] = len(d)
## THE UNIVERSE. A terminus is judgeable here when BOTH of its anchor rows
## resolved -- the terminus part itself and the flank immediately outward --
## because the walk's verdict needs both. This replaces the old `c1`
## membership test, which required those same two rows AND a positive
## terminus density, so the universe is now slightly larger: a terminus with
## zero coverage is no longer excluded as a side effect of computing a depth
## ratio. Inside part_called that changes nothing, since a called part has
## coverage by definition; for no_claim it does, so no_claim counts here are
## NOT comparable to the pre-2026-10-23 tables.
have = {k for k, d in pair.items()
        if "terminus" in d and "terminus_self" in d and k[0] in npart}
print("  judgeable termini (both anchor rows resolved): %d of %d" %
      (len(have), len(pair)))
def rbin(r): return str(r) if r < 5 else ("5-9" if r < 10 else "10+")

## THE EXPECTATION IS MATCHED ON RANK AND ON DISTANCE. rbin alone is not
## enough: inside an exact rank, parts >=1kb from the gene boundary carry the
## long-read end at 0.11 the odds of parts <1kb (5' `exists`), because DEXSeq
## flattening splits a first exon into several parts that all sit within 99 bp
## of the boundary and rank spreads them over ranks 0-3. So a background part
## is comparable only when it matches BOTH. d_near is measured to the boundary
## the terminus faces; the opposite boundary controls position only indirectly
## and is not used.
DB = [0, 100, 300, 1000, 3000, 10000, 30000, 100000]
def dbin(x): return bisect.bisect_right(DB, x) - 1
span_of = {g: (min(a for a, b in d.values()), max(b for a, b in d.values()))
           for g, d in plen.items()}
def dnear(g, pos, low):
    lo, hi = span_of[g]
    return max((pos - lo) if low else (hi - pos), 0)
verd = collections.defaultdict(list)
for r in csv.DictReader(open(f"{PAN}/outward_walk_verdict.tsv"), delimiter="\t"):
    # `gene_edge` stays TRUE geometrically even once a flank exists, so it must
    # NOT be used to skip rows in IG mode -- that would discard exactly the
    # termini the flanks were built to make testable.
    if IG or r["gene_edge"] != "TRUE": verd[r["contrast"]].append(r)

cells = collections.Counter(); hits = collections.Counter(); rank_n = collections.Counter()
bg_n = collections.Counter(); bg_h = collections.Counter()
# EXACT-rank counters, for the MH test at the end. The binned ones above are
# for the o/e tables; MH must not be binned, see the comment on that block.
n_rk = collections.Counter(); h_rk = collections.Counter()
# Generic label counters for the CORRESPONDENCE block at the end. Each criterion
# is a separate binary split, so one terminus contributes to several labels.
lab_n = collections.Counter(); lab_h = collections.Counter()
for trt, ref in PAIRS:
    ds = build(trt, ref)
    for r in verd.get(f"{trt}_vs_{ref}", []):
        key = (r["gene"], r["event"], r["side"])
        ## score_outward_walk.R carries `anchor` itself, so the part number
        ## needs no reconstruction from `rank`.
        anc = "%03d" % int(r["anchor"])
        k = (r["gene"], r["event"], r["side"], anc)
        if k not in have or k not in meta: continue
        g = r["gene"]
        if g not in plen: continue
        co = plen[g].get(anc)   # already canonical
        if not co: continue
        kind, strand = meta[k]
        low = (kind == "5prime") == (strand == "+")
        pos = co[0] if low else co[1]
        rk = (int(anc)-1) if low else (npart[g]-int(anc))
        pc = r["part_called"] == "TRUE"; tc = r["flank_called"] == "TRUE"
        ## the call's claimed lineage. NA for an uncalled terminus, which has no
        ## claim, so its references fall back to the non-directional union.
        try: dpi = float(r["part_dpi"])
        except (KeyError, ValueError): dpi = None
        # THREE-WAY, not binary. C2+ is `part_called AND NOT flank_called`, so it
        # conditions on the NEIGHBOURHOOD being differentially quiet -- while the
        # Q2 reference asks whether a differentially-used boundary is nearby,
        # which is likelier where differential activity is dense. Splitting C2-
        # into its two halves shows whether the inversion is that conditioning:
        #   walk_stops      part_called & NOT flank_called   = C2+
        #   walk_continues  part_called & flank_called       = signal continues outward
        #   no_claim        NOT part_called                 = terminus not differentially used
        grp = "walk_stops" if (pc and not tc) else ("walk_continues" if (pc and tc) else "no_claim")
        e = "5" if kind == "5prime" else "3"
        cell = (kind, grp)
        cells[cell] += 1
        rank_n[(cell, (rbin(max(rk,0)), dbin(dnear(g, pos, low))))] += 1
        n_rk[(kind, grp, max(rk,0))] += 1
        ## LABELS. walk_stops is split by WHICH branch of u_dead fired,
        ## because `not flank_called` has two very different causes: the flank
        ## was quantified and the test said no (differential), or it never
        ## cleared the floor (depth). Only `drop: outward tested, uncalled` is
        ## differential, so the shipped walk_stops is a mixture and has to be
        ## reported beside its halves rather than instead of them.
        rkx = max(rk, 0); tv = r["term_verdict"]
        labs = [grp]
        if pc and not tc:
            labs.append("walkS:diff" if tv == "drop: outward tested, uncalled"
                        else "walkS:depth")
        for L in labs: lab_n[(kind, L, rkx)] += 1
        for d in DEFS:
            if hit(ds, d, e, chrom_of[g], strand, pos, dpi):
                hits[(cell, d)] += 1; h_rk[(kind, grp, d, rkx)] += 1
                for L in labs: lab_h[(kind, L, d, rkx)] += 1
    for g in {r["gene"] for r in verd.get(f"{trt}_vs_{ref}", [])} & set(plen):
        P = sorted(plen[g], key=int); N = len(P); strand = strand_of[g]
        for kind, e in (("5prime","5"), ("3prime","3")):
            low = (kind == "5prime") == (strand == "+")
            for idx, pn in enumerate(P):
                s0, e0 = plen[g][pn]
                pos = s0 if low else e0
                b = (rbin(idx if low else N-1-idx), dbin(dnear(g, pos, low)))
                bg_n[(kind, b)] += 1
                for d in DEFS:
                    ## dpi None: a background part carries no claim to match
                    if hit(ds, d, e, chrom_of[g], strand, pos, None):
                        bg_h[(kind, d, b)] += 1

## Enrichment of each group against the MATCHED BACKGROUND -- a random exonic
## part at the same rank AND the same distance from the gene boundary it
## faces, NOT another group. These are read one at a
## time and must NOT be divided by one another: the group-vs-group question is
## the Mantel-Haenszel test below. obs and exp are pooled as COUNTS and
## divided once, never by averaging two obs/exp values.
def pooled(kind, grp):
    n = cells[(kind, grp)]
    res = {}
    cs = [c for (gk, c) in rank_n if gk == (kind, grp)]
    for d in DEFS:
        obs = hits[((kind, grp), d)]
        exp = sum(rank_n[((kind, grp), c)] *
                  (bg_h[(kind, d, c)] / bg_n[(kind, c)] if bg_n[(kind, c)] else 0)
                  for c in cs)
        res[d] = (obs, exp, (obs / exp) if exp else float('nan'))
    return n, res

print("\n### each walk group against the rank- and distance-matched background")
print("  %-8s %-14s %9s | %-22s | %-22s" %
      ("end", "grp", "n", "exists (Q1)", "exclusive (Q2)"))
POOL = {}
for kind in ("5prime", "3prime"):
    for grp in ("walk_stops", "walk_continues", "no_claim"):
        n, res = pooled(kind, grp)
        POOL[(kind, grp)] = (n, res)
        if not n: continue
        out = []
        for d in DEFS:
            obs, exp, oe = res[d]
            out += ["%5.1f%%" % (100*obs/n), "%5.1f%%" % (100*exp/n), "%6.2f" % oe]
        ## the row follows DEFS; it assumed three references until `share` went
        print(("  %-8s %-14s %9d" + " | %s %s %s" * len(DEFS))
              % (kind, grp, n, *out))

## THE WALK'S OWN DISCRIMINATION, as a Mantel-Haenszel odds ratio stratified by
## EXACT rank.
##
## Until 2026-10-22 this block divided walk_stops' o/e by walk_continues' o/e.
## That is not a valid test. o/e is a ratio of RATES against a common background,
## so its headroom depends on the baseline: at 3' rank 0 the background is ~76%,
## which caps any o/e at 1.32x, while at rank 20 it is ~6% and o/e can reach 17x.
## The two groups do not share a rank distribution -- 31% of 3' walk_stops sits
## at ranks 0-1 against 8% of walk_continues -- so dividing their o/e values
## hands the comparison to whichever group holds the headroom. At 3' that is
## walk_continues, and the division alone inverted the result: the figure this
## block used to print was 0.84x, against 1.61x for the same data compared
## within rank on the odds scale. The 5' direction was unaffected (1.20x -> 2.30x).
##
## MH compares the two groups inside each exact rank and combines on the odds
## scale, which has no ceiling, and it needs no background group at all --
## walk_continues IS the comparison. CI is Robins-Breslow-Greenland. Breslow-Day
## tests the common-OR assumption MH makes; it rejects here, so the single OR is
## a summary of rank-specific ratios that genuinely differ, not a constant
## effect. Quote it as such, and read the per-rank detail from
## outward_stops_split.py before leaning on any one value.
def mh_or(nk, hk, kind, g1, g2, d):
    R = S = 0.0; vp = vm = vq = 0.0; tabs = []
    for rk in sorted({r for (kk, gg, r) in nk if kk == kind and gg in (g1, g2)}):
        n1 = nk[(kind, g1, rk)]; a = hk[(kind, g1, d, rk)]; b = n1 - a
        n2 = nk[(kind, g2, rk)]; c = hk[(kind, g2, d, rk)]; e2 = n2 - c
        n = n1 + n2
        if n1 == 0 or n2 == 0: continue
        R += a*e2/n; S += b*c/n
        P = (a+e2)/n; Q = (b+c)/n
        vp += P*(a*e2/n); vm += P*(b*c/n) + Q*(a*e2/n); vq += Q*(b*c/n)
        tabs.append((a, b, c, e2, n1, n2, a+c, n))
    if R <= 0 or S <= 0: return (float("nan"),)*3 + (0, float("nan"))
    OR  = R/S
    var = vp/(2*R*R) + vm/(2*R*S) + vq/(2*S*S)
    se  = 1.96*math.sqrt(var)
    bd = 0.0; df = 0
    for (a, b, c, e2, n1, n2, m1, n) in tabs:
        # expected a under the common OR: (OR-1)A^2 - A[OR(n1+m1)+(n2-m1)] + OR*n1*m1 = 0
        if m1 == 0 or m1 == n: continue
        Aq = OR - 1.0; Bq = -(OR*(n1+m1) + (n2-m1)); Cq = OR*n1*m1
        if abs(Aq) < 1e-12:
            A = n1*m1/n
        else:
            disc = Bq*Bq - 4*Aq*Cq
            if disc < 0: continue
            cand = [x for x in ((-Bq + sg*math.sqrt(disc))/(2*Aq) for sg in (1, -1))
                    if max(0.0, m1-n2) - 1e-9 <= x <= min(n1, m1) + 1e-9]
            if not cand: continue
            A = cand[0]
        A = min(max(A, 1e-9), min(n1, m1) - 1e-9)
        B, C, E = n1 - A, m1 - A, n2 - m1 + A
        if min(B, C, E) <= 0: continue
        bd += (a - A)**2 * (1.0/A + 1.0/B + 1.0/C + 1.0/E); df += 1
    df = max(df - 1, 0)
    return OR, math.exp(math.log(OR)-se), math.exp(math.log(OR)+se), df, \
           (chi2.sf(bd, df) if df > 0 else float("nan"))

print("\n### THE TEST: walk_stops vs walk_continues, Mantel-Haenszel OR stratified")
print("### by EXACT rank. OR > 1 means the walk_stops group carries the long-read")
print("### end more often AT THE SAME POSITION IN THE GENE.")
print("  %-8s %-10s %9s %-16s %-12s" % ("end", "reference", "MH OR", "95% CI", "Breslow-Day"))
for kind in ("5prime", "3prime"):
    for d in DEFS:
        OR, lo, hi, df, p = mh_or(n_rk, h_rk, kind, "walk_stops", "walk_continues", d)
        print("  %-8s %-10s %9.2f  [%5.2f, %5.2f]   df=%-4d p=%.1e%s" %
              (kind, d, OR, lo, hi, df, p,
               "  heterogeneous" if (p == p and p < 0.05) else ""))
print("  (n stops=%d, continues=%d at 5prime; %d and %d at 3prime)" % (
      POOL[("5prime", "walk_stops")][0], POOL[("5prime", "walk_continues")][0],
      POOL[("3prime", "walk_stops")][0], POOL[("3prime", "walk_continues")][0]))


## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## WHICH BRANCH OF THE WALK CARRIES THE CONCORDANCE. walk_stops is a mixture
## of a differential branch (the flank was tested and not called) and a depth
## branch (the flank never cleared the floor). They are separate claims, so the
## shipped figure is reported beside both halves and never instead of them.
## All three rows use the same comparison group, walk_continues, and the same
## estimator, so they are directly comparable.
print("\n### THE CORRESPONDENCE: each criterion against the reference that matches it")
print("### MH OR stratified by exact rank. Diagonal cells are the intended tests;")
print("### off-diagonal cells are the discriminant check.")
print("  %-8s %-34s %-10s %9s %-16s %-12s" %
      ("end", "criterion (contrast)", "reference", "MH OR", "95% CI", "Breslow-Day"))
## Four criteria, not two. The walk as shipped is the last row and is a MIXTURE
## of the two above it, which is why it has to be reported beside them rather
## than instead of them.
CORR = [
    ("walk depth branches",          "walkS:depth", "walk_continues"),
    ("walk differential branch",     "walkS:diff",  "walk_continues"),
    ("walk as shipped (combined)",   "walk_stops",  "walk_continues"),
]
for kind in ("5prime", "3prime"):
    for nm, g1, g2 in CORR:
        for d in DEFS:
            OR, lo, hi, df, p = mh_or(lab_n, lab_h, kind, g1, g2, d)
            if OR != OR: continue
            print("  %-8s %-34s %-10s %9.2f  [%5.2f, %5.2f]   df=%-4d p=%.1e%s" %
                  (kind, nm, d, OR, lo, hi, df, p,
                   "  het" if (p == p and p < 0.05) else ""))
    print()
