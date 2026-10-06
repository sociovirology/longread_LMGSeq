#!/usr/bin/env python3
"""
Profile alignment lengths in usearch .b6 output, and (optionally) re-derive
per-segment strain calls under a minimum alignment-length filter.

Two modes:

  profile   Report the distribution of alignment lengths per segment, so you
            can see how much of the local-mode signal is short junk (e.g. the
            12 nt matches) versus full-length matches.

  sweep     Re-compute per-segment majority strain calls at a series of
            minimum alignment-length thresholds, and report how reassortment
            rate / genotype completeness respond. This is the discriminating
            test for whether the ~99% reassortment estimate survives once
            short spurious matches are removed.

.b6 is BLAST tabular format (usearch -blast6out):
  1 query  2 target  3 %id  4 alnlen  5 mism  6 gapopen
  7 qlo  8 qhi  9 tlo  10 thi  11 evalue  12 bits

Target IDs are expected as STRAIN_SEGMENT (e.g. PAN99_NS, CA09_PB2).

Usage:
  # One file, length distribution per segment
  python3 b6_length_profile.py profile --b6 plate02_well17_90_merged.b6

  # Whole directory, length distribution pooled
  python3 b6_length_profile.py profile --b6dir genotyping_outputs_local/b6/

  # Threshold sweep across a plate
  python3 b6_length_profile.py sweep --b6dir genotyping_outputs_local/b6/ \\
      --thresholds 0 50 100 200 400 600 800 1000

  # Emit a filtered summary file in the same format as the pipeline's
  # *_strain_assignment_output_all_samples.txt, at one threshold
  python3 b6_length_profile.py sweep --b6dir genotyping_outputs_local/b6/ \\
      --thresholds 200 --emit-summary filtered_min200.txt
"""

import argparse
import glob
import os
import statistics
from collections import defaultdict


def parse_b6(path, min_len=0, min_id=0.0, best_hit_only=True):
    """Yield (read_id, strain, segment, pct_id, aln_len, bits) for hits passing filters.

    best_hit_only: usearch can emit several hits per read. Keep only the
    highest-bitscore hit per read, applied AFTER the length/identity filter,
    so that filtering out a short spurious top hit lets the next-best
    legitimate hit be used rather than discarding the read entirely.
    """
    best = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 12:
                continue
            read_id = parts[0]
            target = parts[1].split(',')[0]
            try:
                pct_id = float(parts[2])
                aln_len = int(parts[3])
                bits = float(parts[11])
            except ValueError:
                continue
            if aln_len < min_len or pct_id < min_id:
                continue
            if '_' not in target:
                continue
            strain, segment = target.split('_', 1)
            rec = (read_id, strain, segment, pct_id, aln_len, bits)
            if not best_hit_only:
                yield rec
            else:
                prev = best.get(read_id)
                if prev is None or bits > prev[5]:
                    best[read_id] = rec
    if best_hit_only:
        for rec in best.values():
            yield rec


def well_name(path):
    return os.path.basename(path).replace('_90_merged.b6', '').replace('.b6', '')


def cmd_profile(args):
    files = [args.b6] if args.b6 else sorted(glob.glob(os.path.join(args.b6dir, '*.b6')))
    if not files:
        print('No .b6 files found.')
        return

    per_seg = defaultdict(list)
    total_hits = 0
    for path in files:
        for _, strain, segment, pct_id, aln_len, _ in parse_b6(path, best_hit_only=True):
            per_seg[segment].append(aln_len)
            total_hits += 1

    print(f'Files: {len(files)}   Best-hit alignments: {total_hits:,}\n')
    SEG_ORDER = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']
    bins = [(0, 50), (50, 100), (100, 200), (200, 400), (400, 800), (800, 1500), (1500, 10 ** 9)]
    header = f'{"seg":<5}{"n":>9}{"median":>8}{"mean":>8}' + ''.join(
        f'{f"{lo}-{hi if hi < 10**9 else ""}":>10}' for lo, hi in bins)
    print(header)
    print('-' * len(header))
    for seg in SEG_ORDER + [s for s in sorted(per_seg) if s not in SEG_ORDER]:
        lens = per_seg.get(seg)
        if not lens:
            continue
        lens_sorted = sorted(lens)
        row = f'{seg:<5}{len(lens):>9,}{lens_sorted[len(lens)//2]:>8}{statistics.mean(lens):>8.0f}'
        for lo, hi in bins:
            c = sum(1 for L in lens if lo <= L < hi)
            row += f'{c/len(lens)*100:>9.1f}%'
        print(row)

    allv = [L for v in per_seg.values() for L in v]
    allv.sort()
    print(f'\nOverall: median={allv[len(allv)//2]} '
          f'p05={allv[int(len(allv)*0.05)]} p95={allv[int(len(allv)*0.95)]}')
    for t in (12, 25, 50, 100, 200):
        print(f'  alignments < {t:>4} nt: {sum(1 for L in allv if L < t):>9,} '
              f'({sum(1 for L in allv if L < t)/len(allv)*100:.1f}%)')


def calls_at_threshold(files, min_len, min_id):
    """Return {well: {segment: (strain, max_for_seg, total_seg)}}"""
    out = {}
    for path in files:
        w = well_name(path)
        counts = defaultdict(lambda: defaultdict(int))  # segment -> strain -> n
        for _, strain, segment, _, _, _ in parse_b6(path, min_len=min_len, min_id=min_id):
            counts[segment][strain] += 1
        wd = {}
        for segment, strains in counts.items():
            top = max(strains.items(), key=lambda kv: kv[1])
            wd[segment] = (top[0], top[1], sum(strains.values()))
        out[w] = wd
    return out


def cmd_sweep(args):
    files = sorted(glob.glob(os.path.join(args.b6dir, '*.b6'))) if args.b6dir else [args.b6]
    if not files:
        print('No .b6 files found.')
        return
    print(f'Files: {len(files)}\n')
    print(f'{"min_len":>8}{"wells":>7}{"usable":>8}{"complete8":>11}{"%compl":>8}'
          f'{"reass_all":>11}{"reass_c8":>10}{"med_prop":>10}{"calls<0.75":>12}')
    print('-' * 85)

    for t in args.thresholds:
        calls = calls_at_threshold(files, t, args.min_id)
        usable = {w: d for w, d in calls.items() if d}
        complete = {w: d for w, d in usable.items() if len(d) == 8}

        def reass(dd):
            return sum(1 for d in dd.values() if len({v[0] for v in d.values()}) > 1)

        props = [v[1] / v[2] for d in usable.values() for v in d.values() if v[2]]
        props.sort()
        med = props[len(props) // 2] if props else float('nan')
        lowconf = sum(1 for p in props if p < 0.75)

        ra = reass(usable) / len(usable) if usable else float('nan')
        rc = reass(complete) / len(complete) if complete else float('nan')
        print(f'{t:>8}{len(calls):>7}{len(usable):>8}{len(complete):>11}'
              f'{len(complete)/len(usable)*100 if usable else 0:>7.0f}%'
              f'{ra:>11.3f}{rc:>10.3f}{med:>10.3f}{lowconf:>12}')

        if args.emit_summary and len(args.thresholds) == 1:
            with open(args.emit_summary, 'w') as f:
                for w in sorted(calls):
                    d = calls[w]
                    if not d:
                        f.write(f'{w}_90_merged.b6\tNONE\tNONE\t0\t0\t0\t0\n')
                        continue
                    for segment, (strain, mx, tot) in d.items():
                        f.write(f'{w}_90_merged.b6\t{segment}\t{strain}_{segment}\t'
                                f'{mx}\t{tot}\t0\t0\n')
            print(f'\nWrote filtered summary to {args.emit_summary}')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    p1 = sub.add_parser('profile', help='alignment-length distribution per segment')
    g1 = p1.add_mutually_exclusive_group(required=True)
    g1.add_argument('--b6')
    g1.add_argument('--b6dir')
    p1.set_defaults(func=cmd_profile)

    p2 = sub.add_parser('sweep', help='re-derive calls across min-length thresholds')
    g2 = p2.add_mutually_exclusive_group(required=True)
    g2.add_argument('--b6')
    g2.add_argument('--b6dir')
    p2.add_argument('--thresholds', type=int, nargs='+',
                    default=[0, 50, 100, 200, 400, 600, 800, 1000])
    p2.add_argument('--min-id', type=float, default=0.0,
                    help='minimum percent identity (0-100), default 0')
    p2.add_argument('--emit-summary', default=None,
                    help='with a single threshold, write a pipeline-format summary file')
    p2.set_defaults(func=cmd_sweep)

    args = ap.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
