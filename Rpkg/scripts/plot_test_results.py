#!/usr/bin/env python3
"""Three-panel figure for one GrASE bipartition: gene model, the nested tests at
that locus, and the path proportion across conditions.

Nested paths are represented as exonic-part SETS rather than as routes through
the graph. A path through an outer bubble fans out when inner bubbles sit inside
it, so it has no single line; but the tested quantity pi = y_D / (y_D + y_S) is
defined on the part sets D and S, which are unambiguous at any nesting depth.
Panel B therefore draws one row per bipartition, ordered by span containment, so
nesting reads as intervals inside intervals.

Example
-------
  python3 plot_bubble_example.py --gene ENSG00000108349 --event 14 \
      --analysis ~/DICE/grase_lineage --kind TSSTTS \
      --conditions B,CD4,CD8,NK,MONO.CLASSIC \
      --symbol CASC3 --note "EJC core component" \
      --out ~/GrASE/media/fig_casc3_bubble.png

Requires: matplotlib, numpy.  Reads only files GrASE already produces.
"""
import argparse, csv, glob, os, sys
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

CB = {"S": "#4477AA", "D1": "#CCBB44", "D2": "#EE6677", "other": "#DDDDDD"}

def read_graph(graphml_path):
    """(node name -> genomic position, (from_name,to_name) -> ex_or_in) from
    the v_name/v_position/e_ex_or_in attributes GrASE writes into its
    GraphML files."""
    ns = "{http://graphml.graphdrawing.org/xmlns}"
    root = ET.parse(graphml_path).getroot()
    key_name = key_pos = key_exin = None
    for k in root.findall(f"{ns}key"):
        an = k.get("attr.name")
        if an == "name": key_name = k.get("id")
        elif an == "position": key_pos = k.get("id")
        elif an == "ex_or_in": key_exin = k.get("id")
    id2name, pos = {}, {}
    for node in root.iter(f"{ns}node"):
        nid = node.get("id")
        nm = pv = None
        for d in node.findall(f"{ns}data"):
            if d.get("key") == key_name: nm = d.text
            if d.get("key") == key_pos: pv = d.text
        if nm is not None:
            id2name[nid] = nm; pos[nm] = pv
    edge_type = {}
    for edge in root.iter(f"{ns}edge"):
        et = None
        for d in edge.findall(f"{ns}data"):
            if d.get("key") == key_exin: et = d.text
        sn, tn = id2name.get(edge.get("source")), id2name.get(edge.get("target"))
        # a node pair can carry both a coarse edge (ex/in/R/L) and a finer
        # ex_part edge covering the same span; vertex-path hops are always
        # coarse, so never let an ex_part edge overwrite a coarse one
        if sn is not None and tn is not None and et != "ex_part":
            edge_type[(sn, tn)] = et
    return pos, edge_type

def parts(cell):
    if cell in ("", "NA", "NaN", None): return []
    return [x.strip() for x in cell.split(",") if x.strip() and x.strip() != "NA"]

def set_length(cell, coords):
    """Summed length in bases of a comma-separated exonic part set.

    None when the set is empty (a junction-substituted side is a point feature
    with no length) or when no part resolves against the GFF.
    """
    ls = [coords[p][1] - coords[p][0] + 1 for p in parts(cell) if p in coords]
    return sum(ls) if ls else None


def pi_perbase(pi, len_d, len_s):
    """Raw-count pi -> per-base pi.

    pi as tested is a COUNT ratio y_D / (y_D + y_S), which the beta-binomial
    likelihood requires but which is not a molecular proportion when D and S
    differ in length: a longer set collects proportionally more reads at equal
    molar concentration. Dividing each count by its feature length first gives

        (D/len_d) / (D/len_d + S/len_s)

    which reduces to the closed form below, so no recounting is needed. Under
    uniform coverage this reads as the fraction of MOLECULES following the
    distinct path, and is therefore comparable across partitions and genes.
    Same form as grase::pi_perbase().
    """
    den = pi * len_s + (1.0 - pi) * len_d
    return (pi * len_s / den) if den > 0 else None


def resolve_gene(gene, gffdir):
    hits = sorted(glob.glob(os.path.join(gffdir, f"{gene}*.dexseq.gff")))
    if not hits:
        sys.exit(f"no dexseq gff for {gene} in {gffdir}")
    return os.path.basename(hits[0]).split(".dexseq")[0], hits[0]

def read_coords(gff):
    coords, strand = {}, "+"
    for ln in open(gff):
        c = ln.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "exonic_part": continue
        try: n = "E" + c[8].split('exonic_part_number "')[1].split('"')[0]
        except IndexError: continue
        coords[n] = (int(c[3]), int(c[4])); strand = c[6]
    return coords, strand

def read_bipartitions(bipdir, gv, kinds):
    out = []
    for k in kinds:
        p = os.path.join(bipdir, f"{gv}.bipartition.{k}.txt")
        if not os.path.exists(p): continue
        for r in csv.DictReader(open(p), delimiter="\t"):
            r["_kind"] = k; out.append(r)
    return out

def read_sj_meta(sjdir, gv, gene_base):
    """event -> {intron_distinct1, intron_distinct2, intron_shared}.

    Only the merged exon+SJ count files carry these; they are constant per
    (gene,event) but written on every sample row, so the first row per event
    wins and the rest are skipped.
    """
    if not sjdir: return {}
    d = os.path.expanduser(sjdir)
    cand = [os.path.join(d, f"{gv}.bipartition.exoncnt.txt")]
    cand += sorted(glob.glob(os.path.join(d, f"{gene_base}.*.bipartition.exoncnt.txt")))
    for f in cand:
        if not os.path.exists(f): continue
        out = {}
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                ev = r.get("event")
                if ev is None or ev in out: continue
                if "intron_distinct1" not in r: return {}   # exonic-only file
                out[ev] = {k: r.get(k, "") for k in
                           ("intron_distinct1", "intron_distinct2", "intron_shared")}
        return out
    return {}


def split_junctions(cell):
    """'chr1:100:200,chr1:300:400' -> [(100,200),(300,400)]; NA/empty -> []."""
    if cell is None: return []
    cell = cell.strip()
    if cell in ("", "NA"): return []
    out = []
    for j in cell.split(","):
        f = j.strip().split(":")
        if len(f) < 3: continue
        try: out.append((int(f[-2]), int(f[-1])))
        except ValueError: continue
    return out


def read_gene_rows(andir, kind, gene_base):
    """every annotated test row for this gene, any event/comparison/contrast"""
    rows = []
    for cand in (os.path.join(andir, f"test_bipartition.{kind}_betabinom_EBmap.annotated.txt"),
                 os.path.join(andir, "exontest_results",
                              f"test_bipartition.{kind}_betabinom_EBmap.annotated.txt")):
        if not os.path.exists(cand): continue
        with open(cand) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["gene"].split(".")[0] == gene_base: rows.append(r)
        break
    return rows

def is_sig(r):
    """Read the call exontest.R made. This script does not decide significance.

    The `significant` column comes from add_significant() / is_significant() in
    the grase package, which is the single authoritative rule: padj, the signed
    lfc_diff_net > delta test, the source-aware dpi floor (min_dpi for exonic
    sides, min_dpi_sj for junction-sourced ones) and the read-support floor.
    A visualisation must not re-derive any of that -- thresholds passed here
    could and did disagree with exontest.R's --padj_threshold, marking calls
    SIG that the pipeline does not make.
    """
    v = r.get("significant")
    if v is None or v == "":
        sys.exit("no `significant` column in the annotated table: re-run "
                 "exontest.R to produce the calls before plotting them")
    return v == "TRUE"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True, help="Ensembl gene id, version optional")
    ap.add_argument("--event", required=True)
    ap.add_argument("--analysis", required=True, help="e.g. ~/DICE/grase_lineage")
    ap.add_argument("--kind", default="TSSTTS", choices=["TSSTTS", "internal"])
    ap.add_argument("--conditions", required=True, help="comma-separated, in plot order")
    ap.add_argument("--comparison", default=None,
                    help="diff1_vs_ref | diff2_vs_ref; default = most significant")
    ap.add_argument("--gffdir", default="~/DICE/dexseq.gff")
    ap.add_argument("--raw", action="store_true",
                    help="plot the raw count ratio y_D/(y_D+y_S) instead of the "
                         "length-normalized per-base proportion (the default). "
                         "Raw is what the test is fitted on; per-base is what is "
                         "comparable across partitions and genes.")
    ap.add_argument("--bipdir", default="~/DICE/bipartition.filtered")
    ap.add_argument("--feature_contrast", default=None,
                    help="contrast to quote in the panel-D caption, e.g. "
                         "MONO.CLASSIC_vs_CD8. Default: the CALLED contrast "
                         "with the largest |delta_pi|.")
    ap.add_argument("--mark_parts", default="",
                    help="annotate named exonic parts in panel A with a bracket and "
                         "label, e.g. 'E019:A,E020:A,E021:B,E022:B,E024:C' to mark the "
                         "CD45 variable exons. Parts sharing a label are bracketed "
                         "together. Use to tie exonic parts to named regions the part "
                         "numbering does not carry -- canonical exon names (CD45 A/B/C) "
                         "or a named intron whose span the parts happen to tile "
                         "(XBP1's 26-nt IRE1 intron is exactly E008+E009).")
    ap.add_argument("--sjcounts", default="",
                    help="merged exon+SJ count dir (per-gene *.exoncnt.txt with "
                         "intron_distinct1/2 columns). When given, sides whose "
                         "exonic distinct set is EMPTY -- the ones that exist only "
                         "because split reads were substituted -- get their splice "
                         "junctions drawn as arcs in panel A. Without it the figure "
                         "is unchanged, so exonic runs are unaffected.")
    ap.add_argument("--graphmldir", default="~/DICE/graphml.v34",
                    help="node name -> genomic position, for panel E's routes")
    ap.add_argument("--stack_kinds", default="internal,TSSTTS",
                    help="which bipartition files to draw in panel B")
    ap.add_argument("--max_stack", type=int, default=0,
                    help="if >0, panel B shows only the N bipartitions whose spans "
                         "are nearest the focal test (dense loci can exceed 90 rows)")
    ap.add_argument("--symbol", default=None); ap.add_argument("--note", default="")
    ap.add_argument("--desc", default=None,
                    help="panel D description; default is derived, e.g. "
                         "'retained-intron path E006 against shared E005'")
    ap.add_argument("--ylab", default=None,
                    help="default: 'path proportion $\\pi$ of distinct D1|D2'")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    gffdir = os.path.expanduser(a.gffdir); bipdir = os.path.expanduser(a.bipdir)
    graphmldir = os.path.expanduser(a.graphmldir)
    andir  = os.path.expanduser(a.analysis)
    base   = a.gene.split(".")[0]
    gv, gff = resolve_gene(base, gffdir)
    coords, strand = read_coords(gff)
    names = sorted(coords, key=lambda e: int(e[1:])); idx = {e: i for i, e in enumerate(names)}
    conds = [c.strip() for c in a.conditions.split(",")]

    # genomic position -> panel x-coordinate, on the same exonic-part scale
    # used by panels A/B, so panel E's routes line up with "where the exons are"
    node_pos, edge_type = {}, {}
    graphml_path = os.path.join(graphmldir, f"{gv}.graphml")
    if os.path.exists(graphml_path):
        node_pos, edge_type = read_graph(graphml_path)
    boundary_x = {}
    for e in names:
        s, en = coords[e]
        boundary_x[s] = idx[e]
        boundary_x[en] = idx[e] + 1
    # the x-axis is always ascending GENOMIC coordinate (per the "genomic
    # order ->" axis label), regardless of strand; R (transcription start)
    # and L (transcription end) sit on opposite sides of that axis depending
    # on strand -- R is leftmost/lowest-coordinate on + strand genes, but
    # rightmost/highest-coordinate on - strand genes (transcription runs
    # toward decreasing genomic position), so their fixed x needs to flip.
    def node_x(name):
        if name == "R": return len(names) + 0.4 if strand == "-" else -0.4
        if name == "L": return -0.4 if strand == "-" else len(names) + 0.4
        try: v = int(node_pos.get(name))
        except (TypeError, ValueError): return None
        return boundary_x.get(v)
    # boundary_x only knows exact exonic-part boundaries. A splice junction's
    # donor/acceptor usually sit ON such boundaries, but not always (and an
    # intronic position never does), so interpolate on the same scale instead of
    # requiring an exact hit.
    _bx = sorted(boundary_x.items())
    def gx(pos):
        if not _bx: return None
        if pos <= _bx[0][0]:  return _bx[0][1]
        if pos >= _bx[-1][0]: return _bx[-1][1]
        for i in range(1, len(_bx)):
            g1, x1 = _bx[i - 1]; g2, x2 = _bx[i]
            if g1 <= pos <= g2:
                if g2 == g1: return x1
                return x1 + (x2 - x1) * (pos - g1) / (g2 - g1)
        return None

    def route_junctions(nodes):
        """genomic (lo,hi) for the INTRON hops on this one route.

        intron_distinct is a UNION over the side's routes (intron_functions.R:
        `in_tx1 <- ekey %in% pkey(pairs1)`), not an intersection, so a junction
        listed for a side need not be on every route of it. Membership must be
        tested per route or the arcs claim more than the data says.
        """
        out = set()
        for nf, nt in zip(nodes[:-1], nodes[1:]):
            if edge_type.get((nf, nt)) != "in": continue
            try: pf, pt = int(node_pos.get(nf)), int(node_pos.get(nt))
            except (TypeError, ValueError): continue
            out.add((min(pf, pt), max(pf, pt)))
        return out

    def hop_parts(n_from, n_to):
        """exonic parts spanned by one coarse hop, empty unless it is an
        'ex' edge (an 'in'/intron hop skips the exonic parts in between).
        Node genomic positions increase along the vertex path on + strand
        genes but decrease on - strand genes, so compare by min/max, not
        by assuming n_from's position is smaller."""
        if edge_type.get((n_from, n_to)) != "ex": return []
        try: p_from, p_to = int(node_pos.get(n_from)), int(node_pos.get(n_to))
        except (TypeError, ValueError): return []
        lo, hi = min(p_from, p_to), max(p_from, p_to)
        return [e for e in names if coords[e][0] >= lo and coords[e][1] <= hi]

    sj_meta = read_sj_meta(a.sjcounts, gv, base)

    stack_kinds = [k.strip() for k in a.stack_kinds.split(",")]
    all_rows = {k: read_gene_rows(andir, k, base) for k in stack_kinds}
    if a.kind not in all_rows: all_rows[a.kind] = read_gene_rows(andir, a.kind, base)
    rows = [r for r in all_rows[a.kind] if r["event"] == str(a.event)]
    if not rows: sys.exit(f"no tests for {base} event {a.event} ({a.kind}) in {andir}")
    # default comparison: the one whose effect a figure should show -- largest
    # |delta_pi| among rows that reach padj<0.05, falling back to smallest padj.
    # Selecting on padj alone can pick a comparison that is significant but has a
    # negligible shift, which makes a misleading panel C.
    if a.comparison:
        comp = a.comparison
    else:
        def eff(r):
            try: return abs(float(r["delta_pi"])) if float(r["padj"]) < 0.05 else -1.0
            except (ValueError, TypeError, KeyError): return -1.0
        # padj is the STRING "NA" for any side exontest.R did not test -- the
        # read-support filter (--min_reads) sets p.value <- NA before adjustment,
        # so this is common, not rare. "NA" is truthy, so `or 1` does not catch
        # it and float() raises; sort those to the back instead.
        def padj_or_1(r):
            try: return float(r.get("padj"))
            except (ValueError, TypeError): return 1.0
        cand = max(rows, key=eff)
        comp = cand["comparison"] if eff(cand) > 0 else \
               min(rows, key=padj_or_1)["comparison"]
    rows = [r for r in rows if r["comparison"] == comp]
    if not rows: sys.exit(f"no rows for comparison {comp}")

    focal = rows[0]
    S, D1, D2 = focal["ref_ex_part"], focal["setdiff1"], focal["setdiff2"]
    which_D  = "D1" if "diff1" in comp else "D2"
    tested_D = D1 if which_D == "D1" else D2
    pi, best = {}, {}
    for r in rows:
        c = r["contrast"]
        if "_vs_" not in c: continue
        t, _, rf = c.partition("_vs_")
        try: pit, pir, pa = float(r["pi_trt"]), float(r["pi_ref"]), float(r["padj"])
        except (ValueError, TypeError): continue
        pi[t] = pit; pi[rf] = pir
        best[c] = pa
    missing = [c for c in conds if c not in pi]
    if missing: sys.exit(f"no pi for conditions: {missing}")

    # Length-normalize unless --raw. The tested quantity is the count ratio, but
    # D and S routinely differ in length, which inflates it by the length ratio
    # and makes it incomparable between partitions. Per-base is the default so
    # the axis agrees with what the text of a case study quotes.
    perbase = not a.raw
    len_D = set_length(tested_D, coords)
    len_S = set_length(S, coords)
    if perbase and (len_D is None or len_S is None):
        # a junction-substituted side has no length, so per-base is undefined
        why = "distinct set is a junction (no length)" if len_D is None \
              else "reference set did not resolve against the GFF"
        print(f"  per-base pi undefined: {why}; plotting the raw ratio")
        perbase = False
    if perbase:
        conv = {c: pi_perbase(v, len_D, len_S) for c, v in pi.items()}
        if any(v is None for v in conv.values()):
            print("  per-base pi undefined for some conditions; plotting raw")
        else:
            pi = conv
            print(f"  per-base pi: len_D {len_D} bp, len_S {len_S} bp")
    # the full call rule, matching panel B's stars and add_significant(): padj
    # AND lfc_diff_net > 0 AND |delta_pi| >= dpi. Counting padj alone overstated
    # support (XBP1 event 5: 7 of 10 by padj, 4 under the rule).
    sig_contrasts = [r["contrast"] for r in rows
                     if r["contrast"] in best and r["comparison"] == comp
                     and is_sig(r)]
    nsig = len(sig_contrasts)

    bips = read_bipartitions(bipdir, gv, [k.strip() for k in a.stack_kinds.split(",")])
    stack = []
    for b in bips:
        pp = parts(b["ref_ex_part"]) + parts(b["setdiff1"]) + parts(b["setdiff2"])
        ii = [idx[p] for p in pp if p in idx]
        if ii: stack.append((min(ii), max(ii), b))
    def _key(a_, b_, c_):
        return (tuple(sorted(parts(a_))), tuple(sorted(parts(b_))), tuple(sorted(parts(c_))))
    focal_key = _key(S, D1, D2)
    matched = next((b for b in bips
                     if _key(b["ref_ex_part"], b["setdiff1"], b["setdiff2"]) == focal_key), None)
    def route_list(cell):
        return [p.strip().split("-") for p in cell.split(",") if p.strip()]
    def n_tx(cell):
        return len([t for t in cell.split(",") if t.strip()])
    e_groups = ([("D1", nodes) for nodes in route_list(matched["path1"])] +
                [("D2", nodes) for nodes in route_list(matched["path2"])]) if matched else []
    e_n_rows = len(e_groups)
    # S can sit on an edge feeding INTO the bubble's source rather than on any
    # hop within the listed routes (e.g. when the source is itself a
    # convergence point with several alternate upstream starts) -- find any
    # such adjacent exonic parts so panel E can still show them.
    # S can sit just OUTSIDE the listed routes: on an edge feeding into the
    # bubble's source (ZCRB1's case -- source is a convergence of alternate
    # starts) or on an edge leaving the bubble's sink (CD226's case -- sink
    # is where routes converge before a shared downstream exon). Check both.
    adjacent_S = set()
    if matched:
        src, snk = matched["source"], matched["sink"]
        for (u, v), et in edge_type.items():
            if et == "ex" and (v == src or u == snk):
                for e in hop_parts(u, v):
                    if e in parts(S):
                        adjacent_S.add(e)
    # per bipartition, which distinct set(s) drove a significant call (D1, D2,
    # or both -- a bipartition can be significant on one comparison, the other,
    # or both across the ten contrasts)
    sig_D = {}
    for k, rs in all_rows.items():
        for r in rs:
            if is_sig(r):
                bkey = _key(r["ref_ex_part"], r["setdiff1"], r["setdiff2"])
                wd = "D1" if "diff1" in r["comparison"] else "D2"
                sig_D.setdefault(bkey, set()).add(wd)
    sig_keys = set(sig_D)
    stack.sort(key=lambda r: (r[0], -(r[1] - r[0])))
    if a.max_stack and len(stack) > a.max_stack:
        fi = [i for i, (lo, hi, b) in enumerate(stack)
              if _key(b["ref_ex_part"], b["setdiff1"], b["setdiff2"]) == focal_key]
        c = fi[0] if fi else len(stack) // 2
        half = a.max_stack // 2
        lo_i, hi_i = max(0, c - half), min(len(stack), max(0, c - half) + a.max_stack)
        n_drop = len(stack) - (hi_i - lo_i)
        stack = stack[lo_i:hi_i]
        print(f"note: panel B truncated to {len(stack)} of {len(stack)+n_drop} bipartitions")

    # ---------------- draw ----------------
    # size panel E's row so its rows sit at the same physical pitch as panel
    # B's (same box height, same row-to-row spacing): both use a 1.0-data-unit
    # row pitch with a 0.9-unit top/bottom margin, so equal inches-per-data-unit
    # across the two panels requires ratio_E / range_E == ratio_B / range_B.
    ratio_A, ratio_B, ratio_CD = 0.42, 2.05, 1.15
    range_B, range_E = len(stack) + 0.9, max(e_n_rows, 1) + 0.9
    ratio_E = ratio_B * range_E / range_B
    # keep A/B/C/D at the same absolute size as before (inches per ratio-unit
    # unchanged); only the figure's total height grows or shrinks to fit E
    old_sum, old_fig_h = ratio_A + ratio_B + 1.3 + ratio_CD, 8.6
    per_unit_h = old_fig_h / old_sum
    new_sum = ratio_A + ratio_B + ratio_E + ratio_CD
    fig = plt.figure(figsize=(11.5, per_unit_h * new_sum))
    gs = fig.add_gridspec(4, 2, height_ratios=[ratio_A, ratio_B, ratio_E, ratio_CD],
                          width_ratios=[3.1, 1], hspace=.5, wspace=.28)
    sym = a.symbol or gv

    axA = fig.add_subplot(gs[0, :])
    EH = 1.0 / 6.0          # exon height: a gene model reads better as a thin track
    for e in names:
        col = CB["other"]
        if e in parts(S): col = CB["S"]
        elif e in parts(D1): col = CB["D1"]
        elif e in parts(D2): col = CB["D2"]
        axA.add_patch(Rectangle((idx[e], 0), .86, EH, facecolor=col,
                                edgecolor="#555555", lw=.4))
    axA.plot([0, len(names)], [EH / 2, EH / 2], color="#888888", lw=.7, zorder=0)
    lbl = list(dict.fromkeys(parts(S)[:1] + parts(D2)[:1] + parts(D1)[:1]))
    for j, e in enumerate(dict.fromkeys(lbl)):
        if e not in idx: continue
        dy = EH + .30 if j % 2 == 0 else EH + .13
        axA.annotate(e, (idx[e] + .43, dy), ha="center", fontsize=6.0)
        axA.plot([idx[e] + .43] * 2, [EH + .02, dy - .03], color="#777777", lw=.5)
    # Canonical-exon brackets. The exonic-part numbering carries no biological
    # names, and mapping a part to e.g. "CD45 variable exon A" by eye is exactly
    # how PTPRC_analysis_notes.md came to misassign two of three (E015 is exon 3,
    # not variable A). Marks are therefore explicit and caller-supplied.
    top = EH + .52
    if a.mark_parts:
        lab_parts = {}
        for tok in a.mark_parts.split(","):
            if ":" not in tok: continue
            pn, lb = tok.split(":", 1)
            pn, lb = pn.strip(), lb.strip()
            if pn in idx: lab_parts.setdefault(lb, []).append(pn)
        y0 = EH + .66
        for lb, pp in lab_parts.items():
            lo = min(idx[q] for q in pp)
            hi = max(idx[q] for q in pp) + .86
            axA.plot([lo, hi], [y0, y0], color="#333333", lw=1.0,
                     solid_capstyle="butt", zorder=6)
            for xe in (lo, hi):
                axA.plot([xe, xe], [y0 - .05, y0], color="#333333", lw=1.0, zorder=6)
            axA.annotate(lb, ((lo + hi) / 2, y0 + .03), ha="center", va="bottom",
                         fontsize=7.5, fontweight="bold", color="#222222", zorder=6)
            top = max(top, y0 + .30)
        if lab_parts:
            axA.annotate("brackets: named annotation regions",
                         (.012, -.055), xycoords="axes fraction", va="top",
                         fontsize=6.6, color="#555555")

    axA.set_xlim(-.6, len(names)); axA.set_ylim(-.34, top); axA.axis("off")
    axA.set_title(f"A   {sym} exonic parts; the tested bipartition ({strand} strand)",
                  loc="left", fontsize=10.5)
    h = [Rectangle((0, 0), 1, 1, facecolor=CB[k]) for k in ("S", "D2", "D1", "other")]
    def _setlab(x, cell=None):
        pp = parts(x)
        if pp: return ",".join(pp)
        # empty exonic set: with junction data loaded, say so rather than "none",
        # which reads as if the side were untested
        if cell and sj_meta.get(str(a.event)):
            nj = len(split_junctions(sj_meta[str(a.event)].get(cell)))
            if nj: return f"{nj} junction{'s' if nj > 1 else ''}, no exonic part"
        return "none"
    axA.legend(h, [f"shared S ({_setlab(S, 'intron_shared')})",
                   f"distinct D2 ({_setlab(D2, 'intron_distinct2')})",
                   f"distinct D1 ({_setlab(D1, 'intron_distinct1')})",
                   "not in this test"],
               loc="lower left", bbox_to_anchor=(0, -.30), ncol=4, fontsize=7.2,
               frameon=False, handlelength=1.1)

    axB = fig.add_subplot(gs[1, :])
    n_marked = 0
    for y, (lo, hi, b) in enumerate(stack):
        yy = len(stack) - y
        bkey = _key(b["ref_ex_part"], b["setdiff1"], b["setdiff2"])
        if bkey == focal_key:
            axB.add_patch(Rectangle((-.6, yy - .46), len(names) + .6, .92,
                                    facecolor="#F2F2F2", edgecolor="none", zorder=0))
        axB.plot([lo, hi + .86], [yy, yy], color="#BBBBBB", lw=.8, zorder=1)
        if bkey in sig_keys:
            star_lbl = "*" + "/".join(sorted(sig_D[bkey]))  # e.g. "*D1", "*D2", "*D1/D2"
            axB.annotate(star_lbl, (-0.38, yy), ha="right", va="center",
                         fontsize=7.5, color="#222222", zorder=4,
                         annotation_clip=False)
            n_marked += 1
        for key, col in (("ref_ex_part", "S"), ("setdiff1", "D1"), ("setdiff2", "D2")):
            for p in parts(b[key]):
                if p in idx:
                    axB.add_patch(Rectangle((idx[p], yy - .32), .86, .64,
                                            facecolor=CB[col], lw=0, zorder=2))
        if _key(b["ref_ex_part"], b["setdiff1"], b["setdiff2"]) == focal_key:
            axB.annotate("tested here", xy=(hi + .95, yy), xycoords="data",
                         xytext=(26, 0), textcoords="offset points",
                         fontsize=8, va="center", ha="left",
                         arrowprops=dict(arrowstyle="->", lw=.9, color="#333333",
                                         shrinkA=1, shrinkB=3), zorder=4)
    axB.set_xlim(-.6, len(names)); axB.set_ylim(.2, len(stack) + 1.1)
    axB.set_yticks([]); axB.set_xticks([])
    for s_ in ("top", "right", "left", "bottom"): axB.spines[s_].set_visible(False)
    axB.set_title(f"B   All {len(stack)} bipartitions tested at this locus, "
                  f"ordered by span containment", loc="left", fontsize=10.5)
    axB.set_xlabel("exonic parts, genomic order  ->", fontsize=8.5)
    axB.annotate(f"each row = one tested bipartition; rows nest where spans are contained\n"
                 f"*D1 / *D2 / *D1/D2 = significant (exontest.R call rule) "
                 f"in at least one contrast, labelled by "
                 f"which distinct set drove it [{n_marked} of {len(stack)} shown]",
                 (.012, .01), xycoords="axes fraction", fontsize=7.8, color="#555555")

    axC = fig.add_subplot(gs[3, 0])
    vals = [pi[c] for c in conds]
    axC.plot(range(len(conds)), vals, "-o", color=CB["D2"], lw=2, ms=7, mec="white", mew=1.2)
    span = max(vals) - min(vals)
    for i, v in enumerate(vals):
        axC.annotate(f"{v:.3f}", (i, v + max(span, .05) * .11), ha="center", fontsize=8)
    axC.set_xticks(range(len(conds)))
    axC.set_xticklabels([c.replace(".CLASSIC", "").replace("MONO", "MONO") for c in conds],
                        fontsize=9)
    axC.set_ylim(min(vals) - max(span, .05) * .30, max(vals) + max(span, .05) * .30)
    # wrap long y-labels so they are not clipped at the figure edge
    yl = a.ylab or ((f"per-base path proportion $\\pi$\nof distinct {which_D}")
                    if perbase else
                    (f"path proportion $\\pi$ (raw)\nof distinct {which_D}"))
    if len(yl) > 26 and "\n" not in yl:
        w = yl.split(" "); half = len(yl) // 2; cur = 0; out = []
        for t in w:
            out.append(t); cur += len(t) + 1
            if cur >= half and t is not w[-1]: out.append("\n"); cur = 0
        yl = " ".join(out).replace(" \n ", "\n")
    axC.set_ylabel(yl, fontsize=9, labelpad=6); axC.grid(axis="y", alpha=.25)
    for s_ in ("top", "right"): axC.spines[s_].set_visible(False)
    axC.set_title("D   Path proportion across conditions", loc="left", fontsize=10.5)

    axD = fig.add_subplot(gs[3, 1]); axD.axis("off")
    axD.text(0, .95, sym + (f"\n({a.note})" if a.note else ""), fontsize=9.5, va="top")
    desc = a.desc or f"distinct {which_D} {tested_D}\nagainst shared S {S}"
    axD.text(0, .64, desc, fontsize=8.5, va="top")
    # Report pi for a CALLED contrast and name it. Previously this printed the
    # first and last entries of --conditions, which is an arbitrary pair: for
    # XBP1 it showed B -> MONO.CLASSIC while the only call was NK vs B, so the
    # figure quoted numbers from a comparison that was not significant.
    if sig_contrasts:
        # Quote the called contrast with the LARGEST |delta_pi|, not the first
        # in file order, which is arbitrary. --feature_contrast overrides.
        if a.feature_contrast and a.feature_contrast in sig_contrasts:
            shown = a.feature_contrast
        else:
            if a.feature_contrast:
                sys.stderr.write(
                    f"note: --feature_contrast {a.feature_contrast} is not among the "
                    f"called contrasts {sig_contrasts}; using the largest effect\n")
            dpi_of = {}
            for r in rows:
                if r["contrast"] in sig_contrasts and r["comparison"] == comp:
                    try: dpi_of[r["contrast"]] = abs(float(r["delta_pi"]))
                    except (ValueError, TypeError): pass
            shown = max(sig_contrasts, key=lambda c: dpi_of.get(c, 0.0))
        t_c, _, r_c = shown.partition("_vs_")
        head = (f"$\\pi$ {pi[r_c]:.3f} ({r_c}) $\\rightarrow$ {pi[t_c]:.3f} ({t_c})"
                f"   [{shown.replace('_vs_', ' vs ')}]")
        if nsig > 1:
            head += f"\nand {nsig - 1} other called contrast(s)"
    else:
        lo_c, hi_c = conds[0], conds[-1]
        head = (f"$\\pi$ {pi[lo_c]:.3f} ({lo_c}) $\\rightarrow$ {pi[hi_c]:.3f} ({hi_c})"
                f"   [no contrast called; range shown]")
    axD.text(0, .34, head + f"\nsignificant (exontest.R call) in {nsig} of "
                            f"{len(best)} contrasts",
             fontsize=8.5, va="top")

    # A side whose exonic distinct set is EMPTY contributes no exonic parts to
    # the test: what is counted for it is the substituted split reads. Drawing
    # its routes as solid exonic boxes therefore shows the wrong quantity -- the
    # boxes are parts the route TRAVERSES, not parts the test COUNTS. Below,
    # such a side is drawn hollow and its counted junctions are drawn on the row.
    sj_side = {}
    if sj_meta.get(str(a.event)):
        for grp, side_parts, cell in (("D1", D1, "intron_distinct1"),
                                      ("D2", D2, "intron_distinct2")):
            if parts(side_parts): continue
            jj = split_junctions(sj_meta[str(a.event)].get(cell))
            if jj: sj_side[grp] = jj

    axE = fig.add_subplot(gs[2, :])
    if matched is None:
        axE.axis("off")
        axE.annotate("no matching bipartition row found for panel E",
                     (.02, .5), xycoords="axes fraction", fontsize=8, color="#888888")
    else:
        groups, n_rows = e_groups, e_n_rows
        n_tx1, n_tx2 = n_tx(matched.get("transcripts1", "")), n_tx(matched.get("transcripts2", ""))
        axE.set_xlim(-.6, len(names))  # same x-range as panels A/B, so exons line up
        axE.set_ylim(.2, n_rows + .9)
        # ratio_E above was sized so this axis has the same inches-per-data-unit
        # as panel B, so reusing B's own row pitch (1.0) and box height (.64)
        # here reproduces the same physical row spacing and box size as panel B.
        box_h = .64
        for y, (grp, nodes) in enumerate(groups):
            yy = n_rows - y
            # R/L have no real genomic position on this axis; a route that
            # jumps straight to L (or starts at R) does not actually traverse
            # the intervening exonic parts, so it must not be drawn as if it
            # covered every exon out to the panel edge.
            real = [(n, node_x(n)) for n in nodes if n not in ("R", "L")]
            real = [(n, x) for n, x in real if x is not None]
            if not real:
                continue
            xs_ok = [x for _, x in real]
            # S can sit just outside the listed route (on the edge feeding
            # into the source) -- fold its box position into the extent so
            # the connecting line and R/L stubs still reach past it
            adj_xs = [idx[e] for e in adjacent_S]
            all_xs = xs_ok + adj_xs
            outer_lo, outer_hi = min(all_xs), max(all_xs) + .86
            axE.plot([outer_lo, outer_hi], [yy, yy], color="#BBBBBB",
                     lw=.8, zorder=1)
            for e in adjacent_S:
                axE.add_patch(Rectangle((idx[e], yy - box_h / 2), .86, box_h,
                                        facecolor=CB["S"], lw=0, zorder=2))
            # exonic-part boxes, same size/position as panel B; blue wherever
            # the route carries the shared S part, else the route's D1/D2 color
            hollow = grp in sj_side
            for n_from, n_to in zip(nodes[:-1], nodes[1:]):
                for e in hop_parts(n_from, n_to):
                    if e in parts(S):
                        axE.add_patch(Rectangle((idx[e], yy - box_h / 2), .86, box_h,
                                                facecolor=CB["S"], lw=0, zorder=2))
                    elif hollow:
                        # traversed but NOT counted: this side's counted signal is
                        # the junction drawn below, not these parts
                        axE.add_patch(Rectangle((idx[e], yy - box_h / 2), .86, box_h,
                                                facecolor="none", edgecolor=CB[grp],
                                                lw=.6, zorder=2))
                    else:
                        axE.add_patch(Rectangle((idx[e], yy - box_h / 2), .86, box_h,
                                                facecolor=CB[grp], lw=0, zorder=2))
            # Draw a junction only on the routes that actually traverse it.
            # Height is capped just above the box top (.32) so an arc cannot
            # reach into the row above.
            if hollow:
                rjs = route_junctions(nodes)
                for (j0, j1) in sj_side[grp]:
                    if (min(j0, j1), max(j0, j1)) not in rjs: continue
                    x0, x1 = gx(j0), gx(j1)
                    if x0 is None or x1 is None: continue
                    if x1 < x0: x0, x1 = x1, x0
                    t = np.linspace(0, np.pi, 40)
                    axE.plot(x0 + (x1 - x0) * (1 - np.cos(t)) / 2,
                             yy + .34 * np.sin(t),
                             color=CB[grp], lw=1.4, solid_capstyle="round", zorder=3)
            # node numbering runs with genomic x on + strand genes but against
            # it on - strand genes, so which physical edge (outer_lo/outer_hi)
            # is the route's start vs end flips with strand
            start_edge, end_edge = (outer_hi, outer_lo) if strand == "-" else (outer_lo, outer_hi)
            if nodes[0] == "R":
                r_x = node_x("R")
                # finer dots: dot diameter scales with lw for ":" style, and
                # an explicit dash pattern keeps them from merging at this width
                axE.plot(sorted([r_x, start_edge]), [yy, yy], ":", color=CB[grp],
                         lw=0.9, dashes=(0.8, 1.6), dash_capstyle="round", zorder=2)
                if r_x <= start_edge:
                    axE.annotate("R ->", (r_x, yy), xytext=(-3, 0), ha="right",
                                 va="center", textcoords="offset points",
                                 fontsize=5.5, color="#555555")
                else:
                    axE.annotate("<- R", (r_x, yy), xytext=(3, 0), ha="left",
                                 va="center", textcoords="offset points",
                                 fontsize=5.5, color="#555555")
            end_label = nodes[-1]
            if end_label == "L":
                l_x = node_x("L")
                axE.plot(sorted([end_edge, l_x]), [yy, yy], ":", color=CB[grp],
                         lw=0.9, dashes=(0.8, 1.6), dash_capstyle="round", zorder=2)
                if l_x >= end_edge:
                    axE.annotate("-> L", (l_x, yy), xytext=(3, 0), ha="left",
                                 va="center", textcoords="offset points",
                                 fontsize=5.5, color="#555555")
                else:
                    axE.annotate("L <-", (l_x, yy), xytext=(-3, 0), ha="right",
                                 va="center", textcoords="offset points",
                                 fontsize=5.5, color="#555555")
            # Interior node numbers (the bubble's sink when it is not L) are
            # NOT labelled: they are internal graph bookkeeping, carry no
            # meaning for a reader, and were the most distracting mark on the
            # panel. Only R and L are annotated, since those say something --
            # the route reaches the transcript start or end.
            axE.annotate(grp, (-.9, yy), ha="right", va="center",
                         fontsize=5.5, color=CB[grp], annotation_clip=False)
        axE.set_yticks([]); axE.set_xticks([])
        for s_ in ("top", "right", "left", "bottom"): axE.spines[s_].set_visible(False)
        n_d1 = sum(1 for grp, _ in groups if grp == "D1")
        n_d2 = sum(1 for grp, _ in groups if grp == "D2")
        axE.set_title(f"C   Tested bipartition {matched['source']}-{matched['sink']} as graph "
                      f"routes along the exons: source (blue) diverges into D1 vs D2 routes",
                      loc="left", fontsize=10.5)
        xlab = (f"exonic parts, genomic order  ->   (D1: {n_d1} route(s), "
                f"{n_tx1} transcript(s) total; D2: {n_d2} route(s), "
                f"{n_tx2} transcript(s) total)")
        if sj_side:
            xlab += ("\n" + " / ".join(sorted(sj_side)) +
                     " has no exonic distinct set: hollow boxes are parts the route "
                     "TRAVERSES, arcs are the counted junctions, drawn on the routes "
                     "that carry them (the side's junction set is a union over its "
                     "routes, so not every route carries every junction)")
        axE.set_xlabel(xlab, fontsize=8.5)

    os.makedirs(os.path.dirname(os.path.expanduser(a.out)) or ".", exist_ok=True)
    fig.savefig(os.path.expanduser(a.out), dpi=200, bbox_inches="tight")
    print(f"wrote {a.out}   ({len(stack)} bipartitions, comparison {comp}, {nsig} sig contrasts)")

if __name__ == "__main__":
    main()
