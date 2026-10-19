#!/usr/bin/env python3
"""Shared helpers for the HABGDB reference panel.

Used by record_live_panel.py (snapshot the live panel into accession
lists), build_hab_refs.py (rebuild the panel from those lists) and
hab_inventory.py (render the inventory pages).

    Tier A  Complete biosynthetic gene cluster (BGC)
    Tier B  Partial BGC / multi-gene fragment
    Tier C  Single gene, complete CDS
    Tier D  Single-gene fragment (partial CDS / partial sequence)

Tiers are assigned per sequence from the NCBI record title (see
`tier_for_title`), with a small number of curated overrides recorded in
the `tier_note` column of each refs/<ref>/accessions.tsv.

NCBI E-utilities are flaky (eutils timed out repeatedly on 2026-06-11 and
again during HABGDB's first build), so every request goes through
`eutils_get`, which retries with exponential backoff on network errors,
HTTP 429/5xx, and the JSON/text error bodies eutils returns with HTTP 200.
"""

import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "HABGDB-build/1.0 (https://github.com/rec3141/HABGDB)"
NCBI_NUCCORE_URL = "https://www.ncbi.nlm.nih.gov/nuccore/{acc}"
NCBI_TAXONOMY_BROWSER_URL = (
    "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={tax_id}"
)

# The seven HAB references in the danaSeq panel, in HABsView.svelte order.
HAB_REFS = [
    "microcystin",
    "cylindrospermopsin",
    "saxitoxin",
    "anatoxin_a",
    "nodularin",
    "guanitoxin",
    "aetokthonotoxin",
]

TIER_META = {
    "A": ("Complete BGC",               "#047857"),
    "B": ("Partial BGC / multi-gene",   "#0e7490"),
    "C": ("Single gene, complete CDS",  "#ca8a04"),
    "D": ("Single-gene fragment",       "#be123c"),
}

ACCESSION_COLUMNS = [
    "order", "seqid", "accession", "seq_start", "seq_stop", "length",
    "record_length", "organism", "taxid", "tier", "tier_note", "mask_n",
    "seq_sha1",
    "title", "note",
]

MINIMAP2_DEFAULT = "/matika/apps/minimap2/minimap2"

# Courtesy delay between eutils calls (NCBI: <=3 req/s without an API key).
_MIN_INTERVAL = 0.4 if not os.environ.get("NCBI_API_KEY") else 0.12
_last_call = [0.0]


# ---------------------------------------------------------------------------
# E-utilities
# ---------------------------------------------------------------------------

def eutils_get(endpoint, params, *, attempts=6, base_backoff=2.0, timeout=60,
               validate=None):
    """POST to an eutils endpoint and return the decoded body.

    Retries up to `attempts` times with exponential backoff (2, 4, 8, 16,
    32 s) on: URLError / socket timeouts, HTTP 429 and 5xx, empty bodies,
    eutils error payloads, and anything `validate(body)` rejects (it should
    raise ValueError). Raises RuntimeError after the last attempt."""
    params = dict(params)
    if os.environ.get("NCBI_API_KEY"):
        params["api_key"] = os.environ["NCBI_API_KEY"]
    params.setdefault("tool", "HABGDB")
    data = urllib.parse.urlencode(params).encode()
    url = f"{EUTILS}/{endpoint}"
    last_err = None
    for attempt in range(1, attempts + 1):
        wait = _MIN_INTERVAL - (time.time() - _last_call[0])
        if wait > 0:
            time.sleep(wait)
        _last_call[0] = time.time()
        try:
            req = urllib.request.Request(url, data=data,
                                         headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                body = r.read().decode("utf-8", "replace")
            if not body.strip():
                raise ValueError("empty response")
            head = body.lstrip()[:400]
            if (head.startswith("{") and '"error"' in head.lower()) or \
               "<ERROR>" in head or "Error occurred" in head:
                raise ValueError(f"eutils error body: {head[:200]!r}")
            if validate is not None:
                validate(body)
            return body
        except urllib.error.HTTPError as e:
            last_err = f"HTTP {e.code}"
            if e.code not in (429, 500, 502, 503, 504):
                raise
        except (urllib.error.URLError, TimeoutError, ConnectionError,
                ValueError, OSError) as e:
            last_err = f"{type(e).__name__}: {e}"
        if attempt < attempts:
            backoff = base_backoff * 2 ** (attempt - 1)
            print(f"  [eutils] {endpoint} attempt {attempt}/{attempts} failed "
                  f"({last_err}); retrying in {backoff:.0f}s", file=sys.stderr)
            time.sleep(backoff)
    raise RuntimeError(f"eutils {endpoint} failed after {attempts} attempts: "
                       f"{last_err}")


def esummary(accessions, batch=100):
    """accession.version → esummary docsum dict (nuccore, JSON)."""
    out = {}
    accs = list(dict.fromkeys(accessions))
    for i in range(0, len(accs), batch):
        chunk = accs[i:i + batch]
        body = eutils_get("esummary.fcgi",
                          {"db": "nuccore", "id": ",".join(chunk),
                           "retmode": "json"})
        res = json.loads(body).get("result", {})
        for uid in res.get("uids", []):
            r = res[uid]
            out[r["accessionversion"]] = r
        missing = [a for a in chunk if a not in out]
        if missing:
            raise RuntimeError(f"esummary returned no record for {missing}")
    return out


def efetch_fasta(accession, seq_start=None, seq_stop=None, *, attempts=6):
    """Fetch one nuccore record (optionally a subrange) as FASTA text.

    Validated: exactly one record whose seqid is `accession` (or
    `accession:start-stop` for subranges) and a non-empty sequence."""
    expect = accession if not seq_start else f"{accession}:{seq_start}-{seq_stop}"
    params = {"db": "nuccore", "id": accession, "rettype": "fasta",
              "retmode": "text"}
    if seq_start:
        params.update(seq_start=str(seq_start), seq_stop=str(seq_stop),
                      strand="1")

    def _validate(body):
        recs = parse_fasta_text(body)
        if len(recs) != 1:
            raise ValueError(f"expected 1 record, got {len(recs)}")
        if recs[0][0] != expect:
            raise ValueError(f"expected seqid {expect}, got {recs[0][0]}")
        if not recs[0][2]:
            raise ValueError("empty sequence")

    body = eutils_get("efetch.fcgi", params, attempts=attempts,
                      validate=_validate)
    return parse_fasta_text(body)[0]


def efetch_feature_table(accession):
    """Fetch a nuccore 5-column feature table and return
    [{start, stop, lo, hi, partial, key, gene, product, locus_tag}] for
    gene/CDS/rRNA-type features (one entry per feature)."""
    def _validate(body):
        if not body.lstrip().startswith(">Feature"):
            raise ValueError("not a feature table")
    body = eutils_get("efetch.fcgi", {"db": "nuccore", "id": accession,
                                      "rettype": "ft", "retmode": "text"},
                      validate=_validate, timeout=180)
    feats = []
    cur = None
    for line in body.splitlines():
        if line.startswith(">"):
            continue
        p = line.split("\t")
        if p[0] and len(p) >= 3 and p[2]:
            a, b = p[0], p[1]
            start, stop = int(a.strip("<>")), int(b.strip("<>"))
            cur = {"start": start, "stop": stop, "lo": min(start, stop),
                   "hi": max(start, stop), "key": p[2],
                   "partial": "<" in a or ">" in a or "<" in b or ">" in b,
                   "gene": "", "product": "", "locus_tag": ""}
            feats.append(cur)
        elif p[0] and cur is not None and len(p) >= 2:
            # additional span of a multi-interval feature
            try:
                v = [int(x.strip("<>")) for x in p[:2]]
                cur["lo"] = min(cur["lo"], *v)
                cur["hi"] = max(cur["hi"], *v)
            except ValueError:
                pass
        elif len(p) >= 5 and cur is not None and p[3] in ("gene", "product", "locus_tag"):
            if not cur[p[3]]:
                cur[p[3]] = p[4]
    return feats


def fetch_rows(rows):
    """seqid → (seqid, description, seq). Full records are fetched in
    batches; subranges one at a time (seq_start applies per request)."""
    got = {}
    full = [r["seqid"] for r in rows if not r["seq_start"]]
    full = list(dict.fromkeys(full))
    for i in range(0, len(full), 50):
        chunk = full[i:i + 50]
        want = set(chunk)

        def _validate(body, want=want):
            ids = {sid for sid, _, seq in parse_fasta_text(body) if seq}
            if ids != want:
                raise ValueError(f"batch mismatch: missing {sorted(want - ids)} "
                                 f"extra {sorted(ids - want)}")

        body = eutils_get("efetch.fcgi",
                             {"db": "nuccore", "id": ",".join(chunk),
                              "rettype": "fasta", "retmode": "text"},
                             validate=_validate, timeout=120)
        for rec in parse_fasta_text(body):
            got[rec[0]] = rec
    for r in rows:
        if r["seq_start"] and r["seqid"] not in got:
            got[r["seqid"]] = efetch_fasta(r["accession"], r["seq_start"],
                                              r["seq_stop"])
    return got



# ---------------------------------------------------------------------------
# Hard-masking (mask_n column)
# ---------------------------------------------------------------------------

def parse_mask(spec):
    """'35345-36265,100-200' → [(35345, 36265), (100, 200)] (1-based,
    inclusive, in panel-sequence coordinates)."""
    out = []
    for part in (spec or "").split(","):
        part = part.strip()
        if part:
            a, b = part.split("-")
            out.append((int(a), int(b)))
    return out


def format_mask(ranges):
    return ",".join(f"{a}-{b}" for a, b in ranges)


def apply_mask(seq, ranges):
    s = list(seq)
    for a, b in ranges:
        s[a - 1:b] = "N" * (b - a + 1)
    return "".join(s)


def derive_mask(ncbi_seq, panel_seq):
    """Return the N-runs that turn ncbi_seq into panel_seq, or None if the
    two differ in any other way (length, or a non-N substitution)."""
    if len(ncbi_seq) != len(panel_seq):
        return None
    runs = []
    for i, (x, y) in enumerate(zip(ncbi_seq, panel_seq), 1):
        if x == y:
            continue
        if y != "N":
            return None
        if runs and runs[-1][1] == i - 1:
            runs[-1][1] = i
        else:
            runs.append([i, i])
    return [tuple(r) for r in runs]


# ---------------------------------------------------------------------------
# FASTA / contigs.json
# ---------------------------------------------------------------------------

def parse_fasta_text(text):
    """[(seqid, description, SEQUENCE)] — sequence uppercased, unwrapped."""
    recs = []
    sid = desc = None
    buf = []
    for line in text.splitlines():
        line = line.rstrip("\r")
        if line.startswith(">"):
            if sid is not None:
                recs.append((sid, desc, "".join(buf).upper()))
            parts = line[1:].split(maxsplit=1)
            sid = parts[0] if parts else ""
            desc = parts[1] if len(parts) > 1 else ""
            buf = []
        elif line.strip():
            buf.append(line.strip())
    if sid is not None:
        recs.append((sid, desc, "".join(buf).upper()))
    return recs


def read_fasta(path):
    with open(path) as f:
        return parse_fasta_text(f.read())


def write_fasta(path, records, width=70):
    """Write [(seqid, description, seq)] wrapped at `width` (efetch's own
    line width)."""
    with open(path, "w") as f:
        for sid, desc, seq in records:
            f.write(f">{sid} {desc}\n" if desc else f">{sid}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + "\n")
            f.write("\n")  # efetch separates records with a blank line


def contigs_payload(records):
    """Same structure as project_zebra/build_contig_index.py: first
    whitespace-delimited token of the header, per-record length, running
    offset in file order (duplicates kept)."""
    contigs = []
    total = 0
    for sid, _desc, seq in records:
        contigs.append({"name": sid, "length": len(seq), "offset": total})
        total += len(seq)
    return {"total_length": total, "contigs": contigs}


def write_contigs_json(path, records):
    payload = contigs_payload(records)
    # json.dump defaults (", " / ": " separators, no trailing newline) —
    # byte-identical to build_contig_index.py's output.
    with open(path, "w") as f:
        json.dump(payload, f)
    return payload


def seq_sha1(seq):
    return hashlib.sha1(seq.upper().encode()).hexdigest()


# ---------------------------------------------------------------------------
# seqid parsing / tiers
# ---------------------------------------------------------------------------

_SUBRANGE_RE = re.compile(r"^(?P<acc>[A-Z]{1,6}_?\d+\.\d+):(?P<start>\d+)-(?P<stop>\d+)$")
_ACC_RE = re.compile(r"^[A-Z]{1,6}_?\d+\.\d+$")


def split_seqid(seqid):
    """'DQ787200.1:1-30300' → ('DQ787200.1', 1, 30300); 'AY212249.1' →
    ('AY212249.1', None, None)."""
    m = _SUBRANGE_RE.match(seqid)
    if m:
        return m["acc"], int(m["start"]), int(m["stop"])
    if _ACC_RE.match(seqid):
        return seqid, None, None
    raise ValueError(f"unrecognised seqid {seqid!r} (expected accession.version"
                     f" or accession.version:start-stop)")


def tier_for_title(title):
    """Classify a sequence by its NCBI title.

    A  title mentions a gene cluster and not 'cluster, partial'
    B  'cluster, partial', or more than one gene named
    C  single gene, 'complete cds'
    D  anything else (partial cds / partial sequence)"""
    t = title.lower()
    if "cluster" in t and "cluster, partial" not in t:
        return "A"
    if "cluster, partial" in t:
        return "B"
    n_genes = len(re.findall(r"\bgenes?\b", t))
    if " genes" in t or n_genes >= 2:
        return "B"
    if "complete cds" in t:
        return "C"
    return "D"


# ---------------------------------------------------------------------------
# accessions.tsv I/O
# ---------------------------------------------------------------------------

def read_accessions(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    for r in rows:
        r["order"] = int(r["order"])
        for k in ("seq_start", "seq_stop", "length", "record_length", "taxid"):
            r[k] = int(r[k]) if r.get(k) else None
    rows.sort(key=lambda r: r["order"])
    return rows


def write_accessions(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ACCESSION_COLUMNS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k))
                        for k in ACCESSION_COLUMNS})


def fmt_bp(n):
    if not n:
        return ""
    if n >= 1_000_000:
        return f"{n/1e6:.2f} Mb"
    if n >= 1_000:
        return f"{n/1e3:.1f} kb"
    return f"{n} bp"
