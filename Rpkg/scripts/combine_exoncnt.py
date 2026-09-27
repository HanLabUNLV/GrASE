#!/usr/bin/env python3
"""
Concatenate per-gene exoncnt files into the combined master exontest reads.

ALIGNS BY COLUMN NAME. The per-gene files do NOT all share one column order --
on the DICE TSS/TTS arm 2,311 files put diff1_source/diff2_source before the
intron_* columns and 689 put them after. A plain `tail -n +2` concatenation
(as in run_exontest_merged_rerun.sh) therefore interleaves junction strings
into the diff2_source column and produces a silently corrupt master.

Emits only the columns exontest needs, which also drops the intron_* columns:
they are 22.5% of the bytes and exontest never reads them.

Streams, so memory is O(one line).

Usage: python3 combine_exoncnt.py <countdir> <out.txt> [--keep-intron]
"""
import sys, os, glob

COLS = ['gene', 'event', 'source', 'sink', 'ref_ex_part', 'setdiff1',
        'setdiff2', 'sample', 'ref', 'diff1', 'diff2', 'groups',
        'diff1_source', 'diff2_source']
INTRON = ['intron_distinct1', 'intron_distinct2', 'intron_shared']

def main():
    countdir, out = sys.argv[1], sys.argv[2]
    cols = COLS + (INTRON if '--keep-intron' in sys.argv else [])
    files = sorted(glob.glob(os.path.join(countdir, '*.bipartition.exoncnt.txt')))
    print('per-gene files: %d' % len(files))
    orders, nrow, nskip = {}, 0, 0
    with open(out, 'w') as o:
        o.write('\t'.join(cols) + '\n')
        for fn in files:
            with open(fn) as fh:
                hdr = fh.readline().rstrip('\n').split('\t')
                orders['\t'.join(hdr)] = orders.get('\t'.join(hdr), 0) + 1
                idx = []
                for c in cols:
                    idx.append(hdr.index(c) if c in hdr else None)
                for line in fh:
                    p = line.rstrip('\n').split('\t')
                    if len(p) != len(hdr):
                        nskip += 1
                        continue
                    o.write('\t'.join('NA' if i is None else p[i] for i in idx) + '\n')
                    nrow += 1
    print('distinct column orders encountered: %d' % len(orders))
    for k, v in sorted(orders.items(), key=lambda x: -x[1]):
        print('  %5d files: %s' % (v, k[:110]))
    print('rows written: %d   malformed rows skipped: %d' % (nrow, nskip))
    print('columns: %d' % len(cols))

if __name__ == '__main__':
    main()
