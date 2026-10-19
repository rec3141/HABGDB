#!/usr/bin/env python3
"""Rebuild the danaSeq HAB mapping references from the committed accession
lists, and optionally verify them against an existing (live) panel.

For each refs/<ref>/accessions.tsv this:

  1. fetches every sequence from NCBI nuccore via EFetch (FASTA), using
     seq_start/seq_stop for :start-stop subranges, with retries and
     exponential backoff (NCBI eutils is flaky);
  2. hard-masks the 1-based ranges in the `mask_n` column to N, then
     checks each sequence against the committed seq_sha1 / length;
  3. writes <outdir>/<ref>/<ref>.fna in accessions.tsv order (duplicates
     kept — the live panel has some, see README);
  4. builds the minimap2 index exactly as mapping_refs/PROTOCOL.md says:
         minimap2 -t 4 -d <ref>.idx <ref>.fna
  5. writes <ref>.contigs.json (same format as project_zebra/
     build_contig_index.py) and copies meta.json + source.txt.

The output dir defaults to ./build (gitignored). Writing into the live
/data/scratch/refdbs/mapping_refs is refused unless --allow-live is given.

With --compare-to DIR, each built reference is compared with
DIR/<ref>/: sequence-level identity of the .fna (same seqids, same
sequences, ignoring order and line wrapping), header descriptions,
contigs.json content, meta.json / source.txt bytes, and .idx bytes
(informational). Exit status is non-zero if any sequence or contigs.json
differs.

Usage:
    python3 build_hab_refs.py                         # all refs → ./build
    python3 build_hab_refs.py --outdir /tmp/hab --refs microcystin nodularin
    python3 build_hab_refs.py --outdir /tmp/hab \\
        --compare-to /data/scratch/refdbs/mapping_refs \\
        --report data/reproducibility_check.tsv
"""
import argparse
import collections
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

import hab_panel as hp

LIVE_DIR = "/data/scratch/refdbs/mapping_refs"


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_one(ref, refs_dir, outdir, minimap2, threads, no_index):
    src = os.path.join(refs_dir, ref)
    rows = hp.read_accessions(os.path.join(src, "accessions.tsv"))
    dst = os.path.join(outdir, ref)
    os.makedirs(dst, exist_ok=True)
    print(f"[build] {ref}: fetching {len(rows)} sequences "
          f"({len({r['seqid'] for r in rows})} unique)", file=sys.stderr)
    got = hp.fetch_rows(rows)

    records, problems = [], []
    for r in rows:
        sid, desc, seq = got[r["seqid"]]
        mask = hp.parse_mask(r.get("mask_n"))
        if mask:
            # Reproduce the live panel's hard-masking (annotated transposase
            # CDSs replaced by N in three cyr records — see README).
            seq = hp.apply_mask(seq, mask)
        if r["length"] and len(seq) != r["length"]:
            problems.append(f"{sid}: length {len(seq)} != committed {r['length']}")
        if r["seq_sha1"] and hp.seq_sha1(seq) != r["seq_sha1"]:
            problems.append(f"{sid}: sequence sha1 differs from committed value "
                            f"(NCBI record changed?)")
        records.append((sid, desc, seq))
    for p in problems:
        print(f"  [build] WARNING {ref}: {p}", file=sys.stderr)

    fna = os.path.join(dst, f"{ref}.fna")
    hp.write_fasta(fna, records)
    hp.write_contigs_json(os.path.join(dst, f"{ref}.contigs.json"), records)
    for fn in ("meta.json", "source.txt"):
        shutil.copy2(os.path.join(src, fn), os.path.join(dst, fn))
    if not no_index:
        idx = os.path.join(dst, f"{ref}.idx")
        # PROTOCOL.md Step 2, verbatim.
        subprocess.run([minimap2, "-t", str(threads), "-d", idx, fna],
                       check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    total = sum(len(s) for _, _, s in records)
    print(f"[build] {ref}: wrote {len(records)} sequences, {total:,} bp"
          f"{'' if no_index else ' + index'}", file=sys.stderr)
    return problems


def compare_one(ref, built_dir, live_dir):
    """Return an ordered dict of check → (ok: bool|None, detail)."""
    b = os.path.join(built_dir, ref)
    l = os.path.join(live_dir, ref)
    res = collections.OrderedDict()
    if not os.path.isdir(l):
        res["live_dir"] = (False, f"{l} missing")
        return res

    bf = hp.read_fasta(os.path.join(b, f"{ref}.fna"))
    lf = hp.read_fasta(os.path.join(l, f"{ref}.fna"))
    bset = collections.Counter((sid, hp.seq_sha1(seq)) for sid, _, seq in bf)
    lset = collections.Counter((sid, hp.seq_sha1(seq)) for sid, _, seq in lf)
    only_b = sorted((bset - lset).elements())
    only_l = sorted((lset - bset).elements())
    res["fna_sequences"] = (
        not only_b and not only_l,
        f"{sum(bset.values())} built vs {sum(lset.values())} live records"
        + (f"; only in build: {[s for s, _ in only_b]}" if only_b else "")
        + (f"; only in live: {[s for s, _ in only_l]}" if only_l else ""))
    res["fna_order"] = (
        [s for s, _, _ in bf] == [s for s, _, _ in lf],
        "same record order" if [s for s, _, _ in bf] == [s for s, _, _ in lf]
        else "record order differs")
    bd = {(sid, d) for sid, d, _ in bf}
    ld = {(sid, d) for sid, d, _ in lf}
    diff_desc = sorted({sid for sid, _ in bd ^ ld})
    res["fna_headers"] = (not diff_desc,
                          "identical descriptions" if not diff_desc
                          else f"description differs for {diff_desc}")

    with open(os.path.join(b, f"{ref}.contigs.json")) as f:
        bc_raw = f.read()
    with open(os.path.join(l, f"{ref}.contigs.json")) as f:
        lc_raw = f.read()
    same_content = json.loads(bc_raw) == json.loads(lc_raw)
    res["contigs_json"] = (same_content,
                           ("identical content" if same_content else "content differs")
                           + ("; byte-identical" if bc_raw == lc_raw else
                              "; not byte-identical"))

    for fn in ("meta.json", "source.txt"):
        same = sha256_file(os.path.join(b, fn)) == sha256_file(os.path.join(l, fn))
        res[fn] = (same, "byte-identical" if same else "differs")

    bi, li = os.path.join(b, f"{ref}.idx"), os.path.join(l, f"{ref}.idx")
    if os.path.exists(bi) and os.path.exists(li):
        same = sha256_file(bi) == sha256_file(li)
        res["idx"] = (same, "byte-identical" if same else
                      f"differs ({os.path.getsize(bi)} vs {os.path.getsize(li)} bytes)")
    else:
        res["idx"] = (None, "not compared (index missing)")
    return res


def minimap2_version(minimap2):
    try:
        return subprocess.run([minimap2, "--version"], capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception as e:
        return f"unknown ({e})"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--refs-dir", default=os.path.join(here, "refs"),
                    help="Directory of refs/<ref>/accessions.tsv (default: ./refs)")
    ap.add_argument("--outdir", default=os.path.join(here, "build"),
                    help="Output directory (default: ./build; never the live panel)")
    ap.add_argument("--refs", nargs="*", default=None,
                    help=f"Subset of references (default: {' '.join(hp.HAB_REFS)})")
    ap.add_argument("--minimap2", default=os.environ.get("MINIMAP2", hp.MINIMAP2_DEFAULT))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--no-index", action="store_true",
                    help="Skip the minimap2 index build")
    ap.add_argument("--compare-to", metavar="DIR",
                    help=f"Verify the build against DIR/<ref>/ (e.g. {LIVE_DIR}); read-only")
    ap.add_argument("--report", metavar="TSV",
                    help="Write the --compare-to results as a TSV")
    ap.add_argument("--allow-live", action="store_true",
                    help=f"Permit --outdir {LIVE_DIR} (overwrites the live panel!)")
    args = ap.parse_args()

    out_real = os.path.realpath(args.outdir)
    live_real = os.path.realpath(LIVE_DIR)
    if (out_real == live_real or out_real.startswith(live_real + os.sep)) \
            and not args.allow_live:
        sys.exit(f"refusing to write into the live panel {LIVE_DIR}; "
                 f"build elsewhere and copy, or pass --allow-live")

    refs = args.refs or hp.HAB_REFS
    if not args.no_index and not os.access(args.minimap2, os.X_OK):
        sys.exit(f"minimap2 not executable: {args.minimap2} (set --minimap2 / $MINIMAP2)")
    mm2 = minimap2_version(args.minimap2) if not args.no_index else None
    os.makedirs(args.outdir, exist_ok=True)

    all_problems = {}
    for ref in refs:
        all_problems[ref] = build_one(ref, args.refs_dir, args.outdir,
                                      args.minimap2, args.threads, args.no_index)

    manifest = {
        "built": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "minimap2": args.minimap2, "minimap2_version": mm2,
        "index_cmd": f"minimap2 -t {args.threads} -d <ref>.idx <ref>.fna",
        "refs": {ref: {fn: sha256_file(os.path.join(args.outdir, ref, fn))
                       for fn in sorted(os.listdir(os.path.join(args.outdir, ref)))}
                 for ref in refs},
        "warnings": {r: p for r, p in all_problems.items() if p},
    }
    with open(os.path.join(args.outdir, "build_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
        f.write("\n")

    if not args.compare_to:
        return
    print(f"\n[verify] {args.outdir} vs {args.compare_to}", file=sys.stderr)
    failed = False
    lines = ["ref\tcheck\tresult\tdetail"]
    for ref in refs:
        res = compare_one(ref, args.outdir, args.compare_to)
        for check, (ok, detail) in res.items():
            word = "n/a" if ok is None else ("PASS" if ok else "DIFF")
            lines.append(f"{ref}\t{check}\t{word}\t{detail}")
            print(f"  {ref:20s} {check:14s} {word:4s}  {detail}", file=sys.stderr)
            if check in ("fna_sequences", "contigs_json") and ok is False:
                failed = True
    if args.report:
        with open(args.report, "w") as f:
            f.write(f"# build_hab_refs.py --compare-to {args.compare_to} "
                    f"({time.strftime('%Y-%m-%d')}, {mm2})\n")
            f.write("\n".join(lines) + "\n")
    print(f"[verify] {'FAILED' if failed else 'all sequence + contigs.json checks passed'}",
          file=sys.stderr)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
