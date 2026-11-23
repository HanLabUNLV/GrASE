#!/usr/bin/env python3
"""
Per-sample counts of the called events, row for row with the called results.

Reads an exontest.R annotated result and the per-gene count files it was tested
on (<countdir>/<gene>.<prefix>.cnt.txt), and writes, one file per contrast
with a call, <out_prefix>.<contrast>.txt:

      Two rows per significant row of the result, in the result's order: the
      tested side D (setdiff1 or setdiff2, from `comparison`) and the reference
      S (ref_ex_part). Columns: contrast, gene, event, comparison, count (D or
      S), parts, source (the count file's diff1/2_source: exon or sj; S is
      always an exonic part), then one column per sample of the contrast's two
      groups (treatment first, samples sorted). Significant row i of that
      contrast is rows 2i-1, 2i.

Stops if a called event has no count rows, a sample of its contrast is missing,
a (gene, event, sample) occurs twice, an event's part or source columns differ
by sample, or a comparison is not diff1_vs_ref / diff2_vs_ref.

Usage: python3 called_counts.py --results FILE --countdir DIR --prefix P
         --samples FILE --out_prefix PATH
  --samples  TSV <signature file> <sample> <group> (stage 30's samples.tsv)
"""
import argparse, csv, os, subprocess, sys, tempfile

COLS = ['event', 'sample', 'groups', 'ref', 'diff1', 'diff2']
EVENT_COLS = ['ref_ex_part', 'setdiff1', 'setdiff2', 'diff1_source', 'diff2_source']
SIDE = {'diff1_vs_ref': ('diff1', 'setdiff1', 'diff1_source'),
        'diff2_vs_ref': ('diff2', 'setdiff2', 'diff2_source')}

ap = argparse.ArgumentParser()
for a in ('results', 'countdir', 'prefix', 'samples', 'out_prefix'):
    ap.add_argument('--' + a, required=True)
o = ap.parse_args()
csv.field_size_limit(sys.maxsize)

by_group = {}
for row in csv.reader(open(o.samples), delimiter='\t'):
    by_group.setdefault(row[2], []).append(row[1])

# the called rows, in result order
calls = []
with open(o.results) as fh:
    r = csv.reader(fh, delimiter='\t'); h = next(r)
    ic, ig, ie, icmp, isg = (h.index(c) for c in ('contrast', 'gene', 'event', 'comparison', 'significant'))
    for row in r:
        if row[isg] != 'TRUE': continue
        if row[icmp] not in SIDE: sys.exit('unexpected comparison %r' % row[icmp])
        calls.append((row[ic], row[ig], row[ie], row[icmp]))
cols_of, tmp = {}, {}
outdir = os.path.dirname(os.path.abspath(o.out_prefix))
for c in dict.fromkeys(x[0] for x in calls):
    t, _, ref = c.partition('_vs_')
    if t not in by_group or ref not in by_group: sys.exit('contrast %s: a group has no samples in --samples' % c)
    cols_of[c] = sorted(by_group[t]) + sorted(by_group[ref])
    tmp[c] = tempfile.NamedTemporaryFile('w', delete=False, dir=outdir)
need = {}
for i, (c, g, e, _) in enumerate(calls):
    need.setdefault(g, {}).setdefault(e, []).append(i)
used = {g for c in cols_of for g in c.split('_vs_')}

for g in sorted(need):
    f = os.path.join(o.countdir, '%s.%s.cnt.txt' % (g, o.prefix))
    ev = need[g]; cnt = {e: {} for e in ev}; info = {}
    with open(f) as fh:
        r = csv.reader(fh, delimiter='\t'); h = next(r); ix = {c: h.index(c) for c in COLS + EVENT_COLS}
        for row in r:
            e = row[ix['event']]
            if e not in ev or row[ix['groups']] not in used: continue
            s = row[ix['sample']]
            if s in cnt[e]: sys.exit('%s event %s: sample %s twice' % (g, e, s))
            cnt[e][s] = {k: row[ix[k]] for k in ('ref', 'diff1', 'diff2')}
            ei = {k: row[ix[k]] for k in EVENT_COLS}
            if info.setdefault(e, ei) != ei: sys.exit('%s event %s: event columns differ by sample' % (g, e))
    for e, idx in ev.items():
        if not cnt[e]: sys.exit('%s event %s: called but no count rows in %s' % (g, e, f))
        for i in idx:
            c, _, _, cmp_ = calls[i]; k, parts, src = SIDE[cmp_]
            miss = [s for s in cols_of[c] if s not in cnt[e]]
            if miss: sys.exit('%s event %s: %d samples of %s missing, e.g. %s' % (g, e, len(miss), c, miss[0]))
            key = [c, g, e, cmp_]
            d = key + ['D', info[e][parts], info[e][src]] + [cnt[e][s][k] for s in cols_of[c]]
            sr = key + ['S', info[e]['ref_ex_part'], 'exon'] + [cnt[e][s]['ref'] for s in cols_of[c]]
            tmp[c].write('%d\t%s\n%d\t%s\n' % (2 * i, '\t'.join(d), 2 * i + 1, '\t'.join(sr)))
# back into result order, one file per contrast
for c, t in tmp.items():
    t.close()
    with open('%s.%s.txt' % (o.out_prefix, c), 'w') as w:
        w.write('\t'.join(['contrast', 'gene', 'event', 'comparison', 'count', 'parts', 'source'] + cols_of[c]) + '\n')
        w.flush()
        p = subprocess.run('LC_ALL=C sort -t "\t" -k1,1n -S 2G -T %s %s | cut -f2-' % (outdir, t.name),
                           shell=True, stdout=w)
    os.unlink(t.name)
    if p.returncode: sys.exit('sort failed')
    print('%s: %d calls, %d samples' % (c, sum(x[0] == c for x in calls), len(cols_of[c])))
print('called rows %d (%d events, %d genes)' % (
    len(calls), sum(len(v) for v in need.values()), len(need)))
