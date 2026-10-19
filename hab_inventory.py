#!/usr/bin/env python3
"""HAB toxin reference panel → per-toxin inventory.

Reads the committed panel records (refs/<ref>/accessions.tsv, meta.json,
source.txt) plus the vendored QA tables in data/ and renders:

    hab_inventory.tsv   one row per toxin reference
    hab_inventory.json  per-toxin summary + every panel sequence
    hab_inventory.html  per-toxin cards + filterable / sortable sequence
                        table for github.io

Everything on the page is derived from those files: producer genera are
the NCBI-taxonomy genera of the panel sequences (data/ncbi_esummary.json),
gaps and exclusions are quoted from the reference's own meta.json /
source.txt, and QA flags are computed (stated vs actual counts, accessions
named in meta.json but absent from the .fna, subranges that drop annotated
genes, hard-masked bases, near-identical sequences).

No network access: re-run record_live_panel.py / panel_redundancy.py to
refresh the inputs.
"""
import argparse
import collections
import html
import json
import os
import re
import sys
import time

import hab_panel as hp

# Compound label + chemistry class exactly as the danaSeq SPA shows them
# (nanopore_live/viz/src/views/HABsView.svelte, HAB_COMPOUNDS).
COMPOUND = {
    "microcystin":        ("Microcystin",        "cyclic heptapeptide"),
    "cylindrospermopsin": ("Cylindrospermopsin", "guanidine alkaloid"),
    "saxitoxin":          ("Saxitoxin",          "PSP alkaloid"),
    "anatoxin_a":         ("Anatoxin-a",         "tropane alkaloid"),
    "nodularin":          ("Nodularin",          "cyclic pentapeptide"),
    "guanitoxin":         ("Guanitoxin",         "phosphoguanidine ester"),
    "aetokthonotoxin":    ("Aetokthonotoxin",    "bis-chlorinated indole alkaloid"),
}

# Known gaps and inclusion / exclusion decisions, quoted or closely
# paraphrased from each reference's meta.json / source.txt (the `src`
# field says which). Only additions from outside those files are the
# NCBI lookups explicitly marked "NCBI, retrieved <date>".
GAPS = {
    "microcystin": [
        ("Non-cyanobacterial mcy orthologs not yet included; add from MAG "
         "assembly when available.", "meta.json notes"),
    ],
    "cylindrospermopsin": [],
    "saxitoxin": [
        ("Gemmata massiliana (planctomycete with sxt orthologs) deferred — "
         "add as separate ref when investigating non-cyano PSP producers.",
         "meta.json notes"),
    ],
    "anatoxin_a": [
        ("Dolichospermum spp. are primary Lake Winnipeg producers but no "
         "standalone ana cluster deposited yet — add from MAG when assembled. "
         "(The panel contains no Dolichospermum sequence.)", "meta.json notes"),
    ],
    "nodularin": [
        ("Baltic/brackish specialist; not yet documented in Great Lakes or "
         "Lake Winnipeg. Present for completeness and ballast water / "
         "migrating bird monitoring.", "meta.json notes"),
    ],
    "guanitoxin": [
        ("Non-cyanobacterial gnt-like reads already detected in Lake "
         "Winnipeg/Great Lakes samples but not yet assembled. When MAG "
         "contigs are assembled, add them here and rebuild.", "meta.json notes"),
        ("Single producer only (Sphaerospermopsis torques-reginae ITEP-024).",
         "accessions.tsv"),
    ],
    "aetokthonotoxin": [
        ("No complete genome deposited in NCBI as of 2026-06-11 — rebuild "
         "when available.", "meta.json notes"),
        ("Scaffold-level A. hydrillicola assemblies do exist (e.g. "
         "GCF_017591595.2, strain Thurmond2011, released 2023-09-18); none is "
         "at Complete Genome level.", "NCBI Datasets, retrieved 2026-10-19"),
    ],
}

# Inconsistencies a regex can't see, each with its evidence.
MANUAL_FLAGS = {
    "guanitoxin": [
        "meta.json still describes the reference as the chromosome "
        "(CP080598.1, ~5 Mb) and its notes say 'The chromosome ... contains the "
        "gnt cluster embedded'; the panel actually holds the 13,572 bp gnt "
        "region (source.txt: 'Previous FASTA was the whole 5.3 Mb chromosome; "
        "trimmed to BGC only'). meta.json was not updated.",
    ],
    "nodularin": [
        "meta.json source/notes ('AY210783 (NSOR10) + AF227755 (CCY9414)', "
        "'2 complete nda cluster sequences') predate the source.txt revision: "
        "AF227755 is listed there as excluded and the panel has 6 sequences "
        "(AY210783.2 + 5 Nostoc).",
    ],
}

EXCLUSIONS = {
    "microcystin": [],
    "cylindrospermopsin": [
        ("FJ418586.4 (Oscillatoria sp. PCC 6506 cyr cluster) included here as "
         "an additional diverse producer; it was taken out of the anatoxin-a "
         "search results.", "meta.json source/notes"),
    ],
    "saxitoxin": [
        ("Gemmata massiliana chromosome excluded (non-cyano ortholog, too large "
         "for this panel).", "meta.json source; source.txt"),
    ],
    "anatoxin_a": [
        ("FJ418586.4 excluded (cyr, added to cylindrospermopsin instead).",
         "meta.json source"),
    ],
    "nodularin": [
        ("AF227755.1 excluded — Sinorhizobium sp. C5 16S rRNA, wrong sequence "
         "(not nda).", "source.txt"),
        ("KJ462515 excluded — Pseudomonas putida J5, UNVERIFIED record; would "
         "cause false positives.", "source.txt"),
        ("DQ842505 excluded — Nodularia spumigena AV1 ndaA, 112 bp only (too "
         "short for reliable mapping).", "source.txt"),
    ],
    "guanitoxin": [
        ("16S sequences from ITEP-026 removed to avoid cross-species false "
         "positives.", "meta.json notes"),
        ("Whole 5.3 Mb chromosome replaced by the gnt BGC region "
         "(CP080598.1:3534366-3547937, gntJ..gntT + 500 bp flanks).",
         "source.txt"),
    ],
    "aetokthonotoxin": [
        ("16S/ITS sequences excluded.", "source.txt"),
    ],
}

_ACC_MENTION = re.compile(r"\b([A-Z]{4,6}\d{8,10}|[A-Z]{1,2}_?\d{5,6})(?:\.\d+)?\b")
_COUNT_MENTION = re.compile(r"\b(\d+)[- ](?:[\w-]+ ){0,3}?sequences?\b", re.I)


def load_tsv(path, skip_comments=True):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        lines = [l for l in f if not (skip_comments and l.startswith("#"))]
    hdr = lines[0].rstrip("\n").split("\t")
    return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in lines[1:] if l.strip()]


def redundancy_groups(pairs, ref):
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            x = parent[x]
        return x

    for p in pairs:
        if p["ref"] == ref:
            parent[find(p["seqid_a"])] = find(p["seqid_b"])
    groups = collections.defaultdict(list)
    for x in list(parent):
        groups[find(x)].append(x)
    return sorted((sorted(g) for g in groups.values() if len(g) > 1),
                  key=lambda g: -len(g))


def summarise(ref, refs_dir, pairs, subfeat):
    d = os.path.join(refs_dir, ref)
    rows = hp.read_accessions(os.path.join(d, "accessions.tsv"))
    with open(os.path.join(d, "meta.json")) as f:
        meta = json.load(f)
    with open(os.path.join(d, "source.txt")) as f:
        source_txt = f.read()
    with open(os.path.join(d, f"{ref}.contigs.json")) as f:
        contigs = json.load(f)

    uniq_ids = list(dict.fromkeys(r["seqid"] for r in rows))
    uniq_seqs = {r["seq_sha1"] for r in rows}
    first = {}
    for r in rows:
        first.setdefault(r["seqid"], r)

    genera = collections.Counter()
    genus_alias = {}
    for r in first.values():
        g = r["organism"].split()[0]
        genera[g] += 1
        tg = r["title"].split()[0]
        if tg != g:
            genus_alias.setdefault(g, set()).add(tg)
    producers = [{"genus": g, "n": n, "title_genus": sorted(genus_alias.get(g, []))}
                 for g, n in genera.most_common()]
    species = sorted({" ".join(r["organism"].split()[:2]) for r in first.values()})

    tiers = collections.Counter(r["tier"] for r in first.values())

    meta_text = " ".join(str(meta.get(k, "")) for k in ("source", "notes"))
    panel_bases = {r["accession"].split(".")[0] for r in rows}
    named = []
    for m in _ACC_MENTION.finditer(meta_text):
        base = m.group(1)
        if base not in named:
            named.append(base)
    named_in = [a for a in named if a in panel_bases]
    named_absent = [a for a in named if a not in panel_bases]

    flags = []
    for label, text in (("meta.json", meta_text), ("source.txt", source_txt)):
        for m in _COUNT_MENTION.finditer(text):
            n = int(m.group(1))
            if n not in (len(rows), len(uniq_ids)):
                flags.append(f"{label} states {m.group(0)!r} but the panel has "
                             f"{len(rows)} records ({len(uniq_ids)} unique seqids)")
                break
    if named_absent:
        flags.append("named in meta.json but not in the panel: "
                     + ", ".join(named_absent))
    dups = len(rows) - len(uniq_ids)
    if dups:
        dup_ids = [s for s, n in collections.Counter(r["seqid"] for r in rows).items() if n > 1]
        flags.append(f"{dups} seqid(s) listed twice: {', '.join(dup_ids)} "
                     f"(duplicate contig names in contigs.json; reads on them get MAPQ 0)")
    ident = [r for r in rows if "identical to row" in (r["note"] or "")]
    for r in ident:
        flags.append(f"{r['seqid']}: {r['note'].split(';')[0]}")
    masked = [r for r in rows if r["mask_n"]]
    if masked:
        nbp = sum(b - a + 1 for r in masked for a, b in hp.parse_mask(r["mask_n"]))
        flags.append(f"{nbp:,} bp in {len(masked)} sequence(s) hard-masked to N "
                     f"relative to NCBI ({', '.join(r['seqid'] for r in masked)}) — "
                     f"exactly the annotated transposase CDSs; not documented in "
                     f"meta.json/source.txt")
    for r in first.values():
        if not r["seq_start"]:
            continue
        dropped = [f for f in subfeat if f["seqid"] == r["seqid"]
                   and f["status"] != "in" and f["feature"] == "gene" and f["gene"]]
        if dropped and r["record_length"] < 1_000_000:
            flags.append(f"{r['seqid']} omits annotated gene(s) "
                         + ", ".join(f"{f['gene']} ({f['lo']}-{f['hi']})" for f in dropped)
                         + f" of the {r['record_length']:,} bp record")
    flags += MANUAL_FLAGS.get(ref, [])
    groups = redundancy_groups(pairs, ref)
    in_groups = sum(len(g) for g in groups)
    if groups:
        flags.append(f"{in_groups} of {len(uniq_ids)} unique seqids fall in "
                     f"{len(groups)} near-identical group(s) (largest: "
                     f"{len(groups[0])}); reads in the shared regions can get "
                     f"MAPQ 0 and are then dropped by danaSeq's mapq>=1 filter")

    return {
        "id": ref,
        "compound": COMPOUND[ref][0],
        "class": COMPOUND[ref][1],
        "common_name": meta.get("common_name", ""),
        "category": meta.get("category", ""),
        "added": meta.get("added", ""),
        "meta_species": meta.get("species", ""),
        "n_records": len(rows),
        "n_unique_seqids": len(uniq_ids),
        "n_unique_sequences": len(uniq_seqs),
        "total_bp": contigs["total_length"],
        "tiers": {t: tiers.get(t, 0) for t in "ABCD"},
        "producer_genera": producers,
        "producer_species": species,
        "named_in_meta": named_in,
        "complete_bgcs": [r["seqid"] for r in first.values() if r["tier"] == "A"],
        "gaps": [{"text": t, "source": s} for t, s in GAPS.get(ref, [])],
        "exclusions": [{"text": t, "source": s} for t, s in EXCLUSIONS.get(ref, [])],
        "qa_flags": flags,
        "redundancy_groups": groups,
        "meta_notes": meta.get("notes", ""),
        "meta_source": meta.get("source", ""),
        "sequences": rows,
    }


def acc_link(seqid):
    acc = seqid.split(":")[0]
    return (f'<a href="{hp.NCBI_NUCCORE_URL.format(acc=acc)}">'
            f'{html.escape(seqid)}</a>')


def write_html(path, toxins, generated):
    esc = html.escape
    n_seq = sum(t["n_records"] for t in toxins)
    tier_tot = collections.Counter()
    for t in toxins:
        for r in {r["seqid"]: r for r in t["sequences"]}.values():
            tier_tot[r["tier"]] += 1
    n_uniq = sum(tier_tot.values())
    tier_pills = " ".join(
        f'<span class="tier-pill" style="background:{hp.TIER_META[k][1]}1a;'
        f'color:{hp.TIER_META[k][1]}"><b>{k}</b> {hp.TIER_META[k][0]} &middot; '
        f'<b>{tier_tot[k]}</b> ({tier_tot[k]/n_uniq*100:.0f}%)</span>'
        for k in "ABCD")

    def tier_badge(k):
        c = hp.TIER_META[k][1]
        return (f'<span class="badge" style="background:{c}1a;color:{c}" '
                f'title="{hp.TIER_META[k][0]}">{k}</span>')

    summary_rows = []
    for t in toxins:
        gen = ", ".join(
            f'<em>{esc(p["genus"])}</em> ({p["n"]})'
            + (f' <span class="muted small">[title: {esc("/".join(p["title_genus"]))}]</span>'
               if p["title_genus"] else "")
            for p in t["producer_genera"])
        tiers = " ".join(f'{tier_badge(k)}&nbsp;{t["tiers"][k]}'
                         for k in "ABCD" if t["tiers"][k])
        summary_rows.append(f"""
        <tr>
          <td data-sort="{esc(t['compound'].lower())}"><a href="#{t['id']}"><b>{esc(t['compound'])}</b></a><br><span class="muted small">{esc(t['class'])}</span></td>
          <td class="small">{gen}</td>
          <td class="mono right" data-sort="{t['n_records']}">{t['n_records']}<br><span class="muted small" title="unique seqids / unique sequences">{t['n_unique_seqids']} / {t['n_unique_sequences']}</span></td>
          <td class="mono right" data-sort="{t['total_bp']}">{hp.fmt_bp(t['total_bp'])}</td>
          <td class="small" data-sort="{t['tiers']['A']}">{tiers}</td>
          <td class="mono right" data-sort="{len(t['qa_flags'])}">{len(t['qa_flags']) or ''}</td>
        </tr>""")

    cards = []
    for t in toxins:
        def lst(items, fmt):
            return ("<ul>" + "".join(f"<li>{fmt(i)}</li>" for i in items) + "</ul>"
                    if items else '<p class="muted small">None recorded.</p>')
        src_fmt = lambda g: f'{esc(g["text"])} <span class="muted small">({esc(g["source"])})</span>'
        red = lst(t["redundancy_groups"],
                  lambda g: f'{len(g)} sequences: <span class="mono small">'
                            + ", ".join(esc(s) for s in g) + "</span>")
        cards.append(f"""
<section class="card" id="{t['id']}">
  <h2>{esc(t['compound'])} <span class="muted small">&middot; <code>{t['id']}</code> &middot; {esc(t['class'])}</span></h2>
  <p class="small"><b>meta.json:</b> {esc(t['common_name'])} &middot; category <code>{esc(t['category'])}</code> &middot; added {esc(t['added'])}<br>
  <b>species field:</b> <em>{esc(t['meta_species'])}</em></p>
  <p class="small"><b>Panel:</b> {t['n_records']} records ({t['n_unique_seqids']} unique seqids, {t['n_unique_sequences']} unique sequences), {t['total_bp']:,} bp &middot;
  <b>NCBI organisms:</b> {', '.join('<em>' + esc(s) + '</em>' for s in t['producer_species'])}</p>
  <p class="small"><b>Named in meta.json:</b> {', '.join(acc_link(a) for a in t['named_in_meta']) or '—'}<br>
  <b>Complete BGCs (tier A):</b> {', '.join(acc_link(a) for a in t['complete_bgcs']) or '—'}</p>
  <h3>Known gaps</h3>{lst(t['gaps'], src_fmt)}
  <h3>Inclusion / exclusion decisions</h3>{lst(t['exclusions'], src_fmt)}
  <h3>QA flags (computed)</h3>{lst(t['qa_flags'], esc)}
  <h3>Near-identical sequence groups (&ge;99.5% id over &ge;90% of the shorter)</h3>{red}
  <details><summary class="small">meta.json source / notes (verbatim)</summary>
    <p class="small"><b>source:</b> {esc(t['meta_source'])}</p>
    <p class="small"><b>notes:</b> {esc(t['meta_notes'])}</p></details>
</section>""")

    seq_rows = []
    for t in toxins:
        for r in t["sequences"]:
            note = r["note"] or ""
            if r["mask_n"]:
                note = (f"masked N: {r['mask_n']}; " + note) if note else f"masked N: {r['mask_n']}"
            seq_rows.append(f"""
        <tr data-ref="{t['id']}" data-tier="{r['tier']}">
          <td data-sort="{t['id']}">{esc(t['compound'])}</td>
          <td class="mono right" data-sort="{r['order']}">{r['order']}</td>
          <td data-sort="{'ABCD'.index(r['tier'])}">{tier_badge(r['tier'])}</td>
          <td class="mono small" data-sort="{esc(r['seqid'])}">{acc_link(r['seqid'])}</td>
          <td class="small" data-sort="{esc(r['organism'].lower())}"><a href="{hp.NCBI_TAXONOMY_BROWSER_URL.format(tax_id=r['taxid'])}"><em>{esc(r['organism'])}</em></a></td>
          <td class="mono right" data-sort="{r['length']}">{hp.fmt_bp(r['length'])}</td>
          <td class="small">{esc(r['title'])}</td>
          <td class="small muted">{esc(note)}{('<br>' + esc(r['tier_note'])) if r['tier_note'] else ''}</td>
        </tr>""")

    ref_options = "".join(f'<option value="{t["id"]}">{esc(t["compound"])} ({t["n_records"]})</option>'
                          for t in toxins)
    tier_options = "".join(f'<option value="{k}">{k} — {hp.TIER_META[k][0]}</option>' for k in "ABCD")

    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>HAB toxin reference panel &rarr; inventory</title>
<style>
  body {{
    font: 15px/1.5 -apple-system, BlinkMacSystemFont, "Inter", "Segoe UI", Roboto, sans-serif;
    color: #1f2937; max-width: 64rem; margin: 2rem auto; padding: 0 1rem;
  }}
  h1 {{ font-size: 1.5rem; margin-bottom: 0.25rem; }}
  h2 {{ font-size: 1.15rem; margin: 0 0 0.4rem; color: #b45309; }}
  h3 {{ font-size: 0.95rem; margin: 1rem 0 0.25rem; color: #1e40af; }}
  .meta {{ color: #475569; font-size: 0.85rem; margin-bottom: 1rem; }}
  .tier-summary {{ margin: 1rem 0 1.5rem; display: flex; flex-wrap: wrap; gap: 0.5rem; }}
  .tier-pill {{ display: inline-block; padding: 0.4rem 0.8rem; border-radius: 6px; font-size: 0.9rem; }}
  .filter-bar {{ margin: 1rem 0; display: flex; flex-wrap: wrap; gap: 0.5rem; align-items: center; }}
  .filter-bar input, .filter-bar select {{
    padding: 0.4rem 0.7rem; font-size: 0.9rem; border: 1px solid #cbd5e1; border-radius: 4px;
  }}
  .filter-bar input {{ flex: 1; min-width: 18rem; max-width: 28rem; }}
  .table-wrap {{ overflow-x: auto; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.9rem; }}
  th, td {{ text-align: left; padding: 0.4rem 0.6rem; border-bottom: 1px solid #e2e8f0; vertical-align: top; }}
  th {{
    font-weight: 600; color: #1e40af; background: #f8fafc; position: sticky; top: 0;
    cursor: pointer; user-select: none;
  }}
  th:hover {{ background: #eef2ff; }}
  th.sorted-asc::after  {{ content: " \\25B2"; opacity: 0.6; font-size: 0.8em; }}
  th.sorted-desc::after {{ content: " \\25BC"; opacity: 0.6; font-size: 0.8em; }}
  .mono {{ font-family: "JetBrains Mono", ui-monospace, monospace; font-size: 0.85em; }}
  .right {{ text-align: right; }}
  .small {{ font-size: 0.85em; }}
  .muted {{ color: #64748b; }}
  .badge {{ display: inline-block; padding: 0.05em 0.55em; border-radius: 9999px; font-weight: 600; font-size: 0.85em; }}
  .card {{ border: 1px solid #e2e8f0; border-radius: 8px; padding: 1rem 1.25rem; background: #fcfcfd; margin: 1rem 0; }}
  .card ul {{ margin: 0.25rem 0; padding-left: 1.25rem; font-size: 0.9rem; }}
  code {{ font-family: "JetBrains Mono", ui-monospace, monospace; background: #f1f5f9; padding: 0.1em 0.35em; border-radius: 4px; font-size: 0.85em; }}
  a {{ color: #1e40af; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
  footer {{ margin-top: 2rem; color: #475569; font-size: 0.8rem; }}
</style>
</head><body>
<h1>HAB toxin reference panel &rarr; inventory</h1>
<div class="meta">Generated {generated} &middot; {len(toxins)} toxin references, {n_seq} panel records ({n_uniq} unique seqids) &middot;
the <code>category: "HAB"</code> references consumed by danaSeq <code>nanopore_live --mapping_refs</code>, hand-built 2026-06-11 from NCBI nuccore searches and
reproduced from the accession lists in <a href="https://github.com/rec3141/HABGDB/tree/main/refs">refs/</a>. Organisms are current NCBI Taxonomy names
(<a href="https://eutils.ncbi.nlm.nih.gov/">Entrez esummary</a>); where the record title still uses an older genus it is shown in brackets.
Outputs: <a href="hab_inventory.tsv">TSV</a> &middot; <a href="hab_inventory.json">JSON</a>.</div>

<div class="tier-summary">{tier_pills}</div>

<div class="table-wrap"><table id="summary">
  <thead><tr>
    <th onclick="sortBy('summary',0)">Toxin</th>
    <th onclick="sortBy('summary',1)">Producer genera in panel (NCBI taxonomy, n seqids)</th>
    <th class="right" onclick="sortBy('summary',2)" title="records; unique seqids / unique sequences">Records</th>
    <th class="right" onclick="sortBy('summary',3)">Total bp</th>
    <th onclick="sortBy('summary',4)">Tiers</th>
    <th class="right" onclick="sortBy('summary',5)">QA flags</th>
  </tr></thead>
  <tbody>{''.join(summary_rows)}</tbody>
</table></div>

<p class="small muted">Tiers (per sequence, from the NCBI record title; overrides noted per row):
{' &middot; '.join(f'<b>{k}</b> {hp.TIER_META[k][0]}' for k in 'ABCD')}.</p>

{''.join(cards)}

<h2 style="margin-top:2rem">All panel sequences</h2>
<div class="filter-bar">
  <input id="filter" placeholder="Filter (accession, organism, title, note)" oninput="applyFilters()">
  <select id="refFilter" onchange="applyFilters()"><option value="">All toxins</option>{ref_options}</select>
  <select id="tierFilter" onchange="applyFilters()"><option value="">All tiers</option>{tier_options}</select>
</div>
<div class="table-wrap"><table id="seqs">
  <thead><tr>
    <th onclick="sortBy('seqs',0)">Toxin</th>
    <th class="right" onclick="sortBy('seqs',1)">#</th>
    <th onclick="sortBy('seqs',2)">Tier</th>
    <th onclick="sortBy('seqs',3)">Seqid</th>
    <th onclick="sortBy('seqs',4)">Organism (NCBI)</th>
    <th class="right" onclick="sortBy('seqs',5)">Length</th>
    <th onclick="sortBy('seqs',6)">Title</th>
    <th onclick="sortBy('seqs',7)">Notes</th>
  </tr></thead>
  <tbody id="rows">{''.join(seq_rows)}</tbody>
</table></div>

<footer>
  Built by <a href="https://github.com/rec3141/HABGDB">HABGDB</a> &middot;
  sequence metadata from NCBI Entrez esummary / feature tables &middot;
  rebuild the references with <code>build_hab_refs.py</code>.
</footer>

<script>
  function applyFilters() {{
    const q   = document.getElementById('filter').value.toLowerCase();
    const ref = document.getElementById('refFilter').value;
    const tr  = document.getElementById('tierFilter').value;
    for (const r of document.querySelectorAll('#rows tr')) {{
      let show = true;
      if (q   && !r.textContent.toLowerCase().includes(q)) show = false;
      if (ref && r.dataset.ref  !== ref) show = false;
      if (tr  && r.dataset.tier !== tr)  show = false;
      r.style.display = show ? '' : 'none';
    }}
  }}
  const sortState = {{}};
  function sortBy(tid, col) {{
    const st = sortState[tid] || (sortState[tid] = {{ col: null, dir: 1 }});
    if (st.col === col) st.dir = -st.dir; else {{ st.col = col; st.dir = 1; }}
    const table = document.getElementById(tid);
    const tbody = table.querySelector('tbody');
    const rows = Array.from(tbody.querySelectorAll('tr'));
    rows.sort((a, b) => {{
      const av = a.children[col].dataset.sort ?? a.children[col].textContent;
      const bv = b.children[col].dataset.sort ?? b.children[col].textContent;
      if (av === '' && bv !== '') return  1;
      if (bv === '' && av !== '') return -1;
      const an = parseFloat(av), bn = parseFloat(bv);
      if (!isNaN(an) && !isNaN(bn)) return (an - bn) * st.dir;
      return av.localeCompare(bv) * st.dir;
    }});
    rows.forEach(r => tbody.appendChild(r));
    table.querySelectorAll('th').forEach((th, i) => {{
      th.classList.toggle('sorted-asc',  i === col && st.dir > 0);
      th.classList.toggle('sorted-desc', i === col && st.dir < 0);
    }});
  }}
</script>
</body></html>
"""
    with open(path, "w") as f:
        f.write(page)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--refs-dir", default=os.path.join(here, "refs"))
    ap.add_argument("--data-dir", default=os.path.join(here, "data"))
    ap.add_argument("--out-tsv",  default=os.path.join(here, "hab_inventory.tsv"))
    ap.add_argument("--out-html", default=os.path.join(here, "hab_inventory.html"))
    ap.add_argument("--out-json", default=os.path.join(here, "hab_inventory.json"))
    args = ap.parse_args()

    pairs = load_tsv(os.path.join(args.data_dir, "redundancy_pairs.tsv"))
    subfeat = load_tsv(os.path.join(args.data_dir, "subrange_features.tsv"))
    toxins = [summarise(ref, args.refs_dir, pairs, subfeat) for ref in hp.HAB_REFS]
    generated = time.strftime("%Y-%m-%d")

    with open(args.out_tsv, "w") as f:
        cols = ["id", "compound", "class", "n_records", "n_unique_seqids",
                "n_unique_sequences", "total_bp", "tier_A", "tier_B", "tier_C",
                "tier_D", "producer_genera", "named_in_meta", "complete_bgcs",
                "known_gaps", "exclusions", "qa_flags"]
        f.write("\t".join(cols) + "\n")
        for t in toxins:
            row = [t["id"], t["compound"], t["class"], t["n_records"],
                   t["n_unique_seqids"], t["n_unique_sequences"], t["total_bp"],
                   *(t["tiers"][k] for k in "ABCD"),
                   "; ".join(f'{p["genus"]} ({p["n"]})' for p in t["producer_genera"]),
                   ", ".join(t["named_in_meta"]), ", ".join(t["complete_bgcs"]),
                   " | ".join(g["text"] for g in t["gaps"]),
                   " | ".join(g["text"] for g in t["exclusions"]),
                   " | ".join(t["qa_flags"])]
            f.write("\t".join(str(x).replace("\t", " ") for x in row) + "\n")

    with open(args.out_json, "w") as f:
        json.dump({"generated": generated,
                   "source": "HABGDB refs/<ref>/accessions.tsv + meta.json + source.txt",
                   "tiers": {k: v[0] for k, v in hp.TIER_META.items()},
                   "toxins": toxins}, f, indent=1)
        f.write("\n")

    write_html(args.out_html, toxins, generated)
    for t in toxins:
        print(f"[inventory] {t['id']:20s} {t['n_records']:3d} records "
              f"{t['total_bp']:>9,} bp  tiers {t['tiers']}  flags {len(t['qa_flags'])}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
