#!/usr/bin/env python3
"""Snapshot the live danaSeq HAB mapping references into version control.

Reads (never writes) the hand-built panel at
/data/scratch/refdbs/mapping_refs/<ref>/ for every subdir whose meta.json
has "category": "HAB", parses the accession.version (plus any
:start-stop subrange) of every sequence from the .fna headers, joins NCBI
esummary metadata (current organism, taxid, full record length), and
writes:

    refs/<ref>/accessions.tsv       ordered sequence list (source of truth
                                    for build_hab_refs.py)
    refs/<ref>/meta.json            verbatim copy of the live file
    refs/<ref>/source.txt           verbatim copy of the live file
    refs/<ref>/<ref>.contigs.json   verbatim copy (reference for verify)
    data/ncbi_esummary.json         trimmed esummary docsums, keyed by
                                    accession.version
    data/live_panel_checksums.tsv   size / mtime / sha256 of every live file
    data/subrange_features.tsv      annotated features of every subranged
                                    record, marked in / OUT / CLIP relative
                                    to the extracted :start-stop window
    data/transposase_features.tsv   transposase CDSs in every panel record
                                    and whether the panel keeps, trims or
                                    hard-masks them

Live sequences are compared with the current NCBI records; differences
that are pure N-masking are recorded in the `mask_n` column (1-based,
inclusive, panel coordinates) so build_hab_refs.py can reproduce them.

This is a one-off snapshot tool: re-running it overwrites the curated
accession lists with whatever is currently live. Day-to-day edits to the
panel should go into refs/<ref>/accessions.tsv directly, followed by
build_hab_refs.py.
"""
import argparse
import glob
import hashlib
import json
import os
import re
import shutil
import sys
import time

import hab_panel as hp

LIVE_DEFAULT = "/data/scratch/refdbs/mapping_refs"

# Curated tier overrides where the NCBI title alone misclassifies the
# sequence. Each carries the evidence it rests on.
TIER_OVERRIDES = {
    "AY210783.2": ("A", "override: source.txt lists AY210783.2 as the "
                        "'complete nda cluster (canonical)'; the NCBI title "
                        "enumerates ndaA-ndaI CDSs without the word 'cluster'"),
    "CP080598.1:3534366-3547937": ("A", "override: source.txt — gnt BGC "
                        "region (gntJ..gntT, incl. 500 bp flanks) extracted "
                        "from the ITEP-024 chromosome with EFetch "
                        "seq_start/seq_stop"),
}

ESUMMARY_KEEP = ["accessionversion", "title", "organism", "taxid", "slen",
                 "moltype", "topology", "sourcedb", "createdate",
                 "updatedate", "strain", "biomol", "extra"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def discover_hab_refs(live):
    refs = []
    for mj in sorted(glob.glob(os.path.join(live, "*", "meta.json"))):
        with open(mj) as f:
            md = json.load(f)
        if (md.get("category") or "").split("-")[0].upper() == "HAB":
            refs.append(os.path.basename(os.path.dirname(mj)))
    return refs


def build_rows(ref, records, summ, ncbi):
    rows = []
    first_by_sha = {}
    first_by_id = {}
    for i, (sid, desc, seq) in enumerate(records, 1):
        acc, start, stop = hp.split_seqid(sid)
        doc = summ[acc]
        tier = hp.tier_for_title(desc)
        tier_note = ""
        if sid in TIER_OVERRIDES:
            tier, tier_note = TIER_OVERRIDES[sid]
        elif acc in TIER_OVERRIDES:
            tier, tier_note = TIER_OVERRIDES[acc]
        notes = []
        # Any difference from the NCBI record must be explainable as
        # hard-masking (bases → N); anything else is a hard error.
        nseq = ncbi[sid][2]
        mask = []
        if nseq != seq:
            mask = hp.derive_mask(nseq, seq)
            if mask is None:
                raise SystemExit(f"{ref} {sid}: live sequence differs from NCBI "
                                 f"by more than N-masking; cannot record it")
            notes.append(f"{sum(b - a + 1 for a, b in mask)} bp hard-masked to N "
                         f"vs NCBI (undocumented in meta.json/source.txt)")
        if start:
            notes.append(f"subrange {start}-{stop} of a {doc['slen']} bp record")
        elif len(seq) != doc["slen"]:
            notes.append(f"length {len(seq)} != NCBI slen {doc['slen']}")
        sha = hp.seq_sha1(seq)
        if sid in first_by_id:
            notes.append(f"exact duplicate of row {first_by_id[sid]} "
                         f"(same seqid listed twice)")
        elif sha in first_by_sha:
            j, other = first_by_sha[sha]
            notes.append(f"sequence identical to row {j} ({other})")
        first_by_id.setdefault(sid, i)
        first_by_sha.setdefault(sha, (i, sid))
        n_amb = sum(1 for c in nseq if c not in "ACGT")
        if n_amb:
            notes.append(f"{n_amb} non-ACGT bases in the NCBI record")
        if len(seq) < 200:
            notes.append(f"only {len(seq)} bp")
        org = doc["organism"]
        t_genus, o_genus = desc.split()[0], org.split()[0]
        if t_genus != o_genus:
            notes.append(f"title genus {t_genus}, NCBI taxonomy genus {o_genus}")
        elif not desc.startswith(org):
            notes.append(f"title organism differs from NCBI organism '{org}'")
        rows.append({
            "order": i, "seqid": sid, "accession": acc,
            "seq_start": start, "seq_stop": stop, "length": len(seq),
            "record_length": doc["slen"], "organism": org,
            "taxid": doc["taxid"], "tier": tier, "tier_note": tier_note,
            "mask_n": hp.format_mask(mask),
            "seq_sha1": sha, "title": desc, "note": "; ".join(notes),
        })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--live", default=LIVE_DEFAULT,
                    help=f"Live mapping_refs dir, read-only (default {LIVE_DEFAULT})")
    ap.add_argument("--repo", default=here)
    args = ap.parse_args()

    refs = discover_hab_refs(args.live)
    print(f"[record] HAB refs in {args.live}: {', '.join(refs)}", file=sys.stderr)
    if sorted(refs) != sorted(hp.HAB_REFS):
        print(f"[record] WARNING: HAB refs differ from hab_panel.HAB_REFS "
              f"({hp.HAB_REFS})", file=sys.stderr)

    per_ref = {}
    all_accs = []
    for ref in refs:
        recs = hp.read_fasta(os.path.join(args.live, ref, f"{ref}.fna"))
        per_ref[ref] = recs
        all_accs += [hp.split_seqid(sid)[0] for sid, _, _ in recs]

    summ = hp.esummary(all_accs)

    # Fetch the NCBI version of every sequence so live/NCBI differences
    # (hard-masking) are captured in the mask_n column.
    ncbi = {}
    for ref in refs:
        rows = []
        for sid, _, _ in per_ref[ref]:
            acc, start, stop = hp.split_seqid(sid)
            rows.append({"seqid": sid, "accession": acc,
                         "seq_start": start, "seq_stop": stop})
        ncbi.update(hp.fetch_rows(rows))

    checks = []
    panel_rows = {}
    for ref in refs:
        src = os.path.join(args.live, ref)
        dst = os.path.join(args.repo, "refs", ref)
        os.makedirs(dst, exist_ok=True)
        rows = build_rows(ref, per_ref[ref], summ, ncbi)
        panel_rows[ref] = rows
        hp.write_accessions(os.path.join(dst, "accessions.tsv"), rows)
        for fn in ("meta.json", "source.txt", f"{ref}.contigs.json"):
            shutil.copy2(os.path.join(src, fn), os.path.join(dst, fn))
        for fn in sorted(os.listdir(src)):
            p = os.path.join(src, fn)
            st = os.stat(p)
            checks.append((ref, fn, st.st_size,
                           time.strftime("%Y-%m-%dT%H:%M:%S",
                                         time.localtime(st.st_mtime)),
                           sha256(p)))
        print(f"[record] {ref}: {len(rows)} sequences, "
              f"{sum(r['length'] for r in rows):,} bp", file=sys.stderr)

    with open(os.path.join(args.repo, "data", "ncbi_esummary.json"), "w") as f:
        json.dump({"retrieved": time.strftime("%Y-%m-%d"),
                   "records": {a: {k: summ[a].get(k) for k in ESUMMARY_KEEP}
                               for a in sorted(summ)}},
                  f, indent=1, sort_keys=False)
        f.write("\n")
    with open(os.path.join(args.repo, "data", "live_panel_checksums.tsv"), "w") as f:
        f.write("ref\tfile\tbytes\tmtime\tsha256\n")
        for row in checks:
            f.write("\t".join(map(str, row)) + "\n")
    fts = {}
    for ref in refs:
        for r in panel_rows[ref]:
            if r["accession"] not in fts:
                fts[r["accession"]] = hp.efetch_feature_table(r["accession"])

    # Which annotated features does each :start-stop extraction keep?
    # gene / CDS / rRNA / tRNA features carrying a gene, product or
    # locus_tag; for records > 1 Mb only features within 5 kb of the window.
    with open(os.path.join(args.repo, "data", "subrange_features.tsv"), "w") as f:
        f.write("ref\tseqid\taccession\tstatus\tfeature\tlo\thi\tgene\t"
                "product\tlocus_tag\n")
        for ref in refs:
            for r in panel_rows[ref]:
                acc, start, stop = r["accession"], r["seq_start"], r["seq_stop"]
                if not start:
                    continue
                for ft in fts[acc]:
                    if ft["key"] not in ("CDS", "rRNA", "tRNA", "gene"):
                        continue
                    if not (ft["gene"] or ft["product"] or ft["locus_tag"]):
                        continue
                    if summ[acc]["slen"] > 1_000_000 and (
                            ft["hi"] < start - 5000 or ft["lo"] > stop + 5000):
                        continue
                    if ft["lo"] >= start and ft["hi"] <= stop:
                        status = "in"
                    elif ft["hi"] < start or ft["lo"] > stop:
                        status = "OUT"
                    else:
                        status = "CLIP"
                    f.write("\t".join(map(str, [
                        ref, r["seqid"], acc, status, ft["key"], ft["lo"],
                        ft["hi"], ft["gene"], ft["product"], ft["locus_tag"]])) + "\n")

    # Every CDS annotated as a transposase / IS element in any panel
    # sequence, and what the panel does with it (kept, outside the
    # extracted subrange, or hard-masked).
    with open(os.path.join(args.repo, "data", "transposase_features.tsv"), "w") as f:
        f.write("ref\torder\tseqid\tlo\thi\tgene\tproduct\tpanel_status\n")
        for ref in refs:
            for r in panel_rows[ref]:
                acc, start, stop = r["accession"], r["seq_start"], r["seq_stop"]
                off = (start - 1) if start else 0
                lim = stop if start else summ[acc]["slen"]
                mask = hp.parse_mask(r["mask_n"])
                for ft in fts[acc]:
                    if ft["key"] != "CDS" or not re.search(
                            r"transpos|insertion element|\bIS\d", ft["product"], re.I):
                        continue
                    if summ[acc]["slen"] > 1_000_000 and (
                            ft["hi"] < (start or 0) or ft["lo"] > lim):
                        continue
                    if ft["hi"] <= off or ft["lo"] > lim:
                        status = "outside subrange"
                    elif ft["lo"] <= off or ft["hi"] > lim:
                        status = "clipped by subrange"
                    else:
                        lo, hi = ft["lo"] - off, ft["hi"] - off
                        if any(a <= lo and hi <= b for a, b in mask):
                            status = "hard-masked to N"
                        elif any(a <= hi and lo <= b for a, b in mask):
                            status = "partially masked"
                        else:
                            status = "in panel (unmasked)"
                    f.write("\t".join(map(str, [ref, r["order"], r["seqid"],
                            ft["lo"], ft["hi"], ft["gene"], ft["product"],
                            status])) + "\n")
    print("[record] done.", file=sys.stderr)


if __name__ == "__main__":
    main()
