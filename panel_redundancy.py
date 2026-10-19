#!/usr/bin/env python3
"""Within-reference redundancy check for the HAB panel.

danaSeq maps with `minimap2 -ax map-ont --secondary=no` and then keeps only
alignments with MAPQ >= 1 (bin/filter_minimap2.awk). When two or more
sequences inside one reference are (near-)identical, a read that falls in
the shared region has equally good hits and gets MAPQ 0, so it is dropped.
This script finds those sequence pairs by aligning each built <ref>.fna
against itself:

    minimap2 -t 4 -cx asm20 -X <ref>.fna <ref>.fna

and keeps pairs with alignment identity >= --min-id over >= --min-cov of
the shorter sequence. Output: data/redundancy_pairs.tsv (consumed by
hab_inventory.py).

Usage (after build_hab_refs.py):
    python3 panel_redundancy.py --build-dir build
"""
import argparse
import os
import subprocess
import sys

import hab_panel as hp


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--build-dir", default=os.path.join(here, "build"))
    ap.add_argument("--minimap2", default=os.environ.get("MINIMAP2", hp.MINIMAP2_DEFAULT))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-id", type=float, default=0.995)
    ap.add_argument("--min-cov", type=float, default=0.9)
    ap.add_argument("--out", default=os.path.join(here, "data", "redundancy_pairs.tsv"))
    args = ap.parse_args()

    out_rows = []
    for ref in hp.HAB_REFS:
        fna = os.path.join(args.build_dir, ref, f"{ref}.fna")
        if not os.path.exists(fna):
            sys.exit(f"missing {fna}; run build_hab_refs.py first")
        paf = subprocess.run(
            [args.minimap2, "-t", str(args.threads), "-cx", "asm20", "-X", fna, fna],
            capture_output=True, text=True, check=True).stdout
        best = {}
        for line in paf.splitlines():
            p = line.split("\t")
            q, qlen, t, tlen = p[0], int(p[1]), p[5], int(p[6])
            if q == t:
                continue
            nmatch, alen = int(p[9]), int(p[10])
            ident = nmatch / alen
            cov = alen / min(qlen, tlen)
            if ident < args.min_id or cov < args.min_cov:
                continue
            a, b = sorted([(q, qlen), (t, tlen)])
            key = (a[0], b[0])
            if key not in best or alen > best[key][2]:
                best[key] = (a[1], b[1], alen, ident, min(cov, 1.0))
        for (a, b), (alen_a, alen_b, alen, ident, cov) in sorted(best.items()):
            out_rows.append((ref, a, alen_a, b, alen_b, alen, f"{ident:.4f}", f"{cov:.2f}"))
        print(f"[redundancy] {ref}: {len(best)} pairs >= {args.min_id:.1%} id "
              f"over >= {args.min_cov:.0%} of the shorter sequence", file=sys.stderr)

    with open(args.out, "w") as f:
        f.write(f"# minimap2 -cx asm20 -X self-alignment; identity >= {args.min_id} "
                f"over >= {args.min_cov} of the shorter sequence\n")
        f.write("ref\tseqid_a\tlength_a\tseqid_b\tlength_b\taln_len\tidentity\tcov_shorter\n")
        for r in out_rows:
            f.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
