# HABGDB

**Harmful Algal Bloom Genomic Database**: a version-controlled record
of the cyanotoxin reference panel that danaSeq's `nanopore_live` pipeline
uses to map reads for shotgun-metagenomic HAB surveillance.

This is a sibling of [AISGDB](https://github.com/rec3141/AISGDB) (aquatic
invasive species) and [SARGDB](https://github.com/rec3141/SARGDB) (species
at risk). All three describe reference directories with the same
`mapping_refs/<name>/` layout. AISGDB and SARGDB tier *species* by the
NCBI genome resources available for them. HABGDB tiers the *sequences* of
a toxin biosynthetic-gene-cluster (BGC) panel by how complete each BGC is.

Browse the rendered docs at **<https://rec3141.github.io/HABGDB/>**.

## Contents

| File | What it is |
|---|---|
| [`index.html`](index.html) | Landing page for github.io |
| [`hab_inventory.html`](hab_inventory.html) | For each toxin: producer genera, panel contents, key accessions, known gaps, exclusion decisions and computed QA flags. Also a filterable, sortable table of all 102 panel sequences |
| [`hab_inventory.{tsv,json}`](hab_inventory.tsv) | Machine-readable versions of the inventory (the JSON includes every sequence row) |
| [`hab_inventory.py`](hab_inventory.py) | Renders the inventory from `refs/` + `data/`, with no network access |
| [`refs/<ref>/accessions.tsv`](refs/) | **Source of truth.** The ordered list of sequences in each reference: seqid (accession.version, plus `:start-stop` for subranges), NCBI organism / taxid / record length, BGC tier, N-mask ranges, SHA-1 of the panel sequence, NCBI title, and notes |
| `refs/<ref>/meta.json`, `source.txt` | Verbatim copies of the live files, holding the original search queries, notes, and inclusion / exclusion decisions |
| `refs/<ref>/<ref>.contigs.json` | Verbatim copy of the live contig sidecar, used as the expected value when verifying a rebuild |
| [`build_hab_refs.py`](build_hab_refs.py) | Rebuilds `<ref>.fna`, `<ref>.idx`, `<ref>.contigs.json` and `meta.json` from NCBI. With `--compare-to` it also verifies the rebuild against an existing panel |
| [`record_live_panel.py`](record_live_panel.py) | One-off snapshot tool that produced `refs/` and `data/` from the live panel (read-only on the live side) |
| [`panel_redundancy.py`](panel_redundancy.py) | Self-alignment of each built reference to find near-identical sequences (a MAPQ-0 risk, see below) |
| [`hab_panel.py`](hab_panel.py) | Shared helpers: E-utilities with retries and backoff, FASTA / contigs.json I/O, the tier rule, N-masking |
| [`data/`](data/) | Snapshots and QA tables (see below) |

## The panel

The panel was built on **2026-06-11** from NCBI nuccore searches, by hand, directly in
`/data/scratch/refdbs/mapping_refs/` and with no version control. It covers
seven toxin references, which are the subdirectories whose `meta.json` has
`"category": "HAB"`:

| Reference | Compound class (HABsView) | Records | bp | Producer genera in panel (NCBI taxonomy) |
|---|---|---:|---:|---|
| `microcystin` | cyclic heptapeptide | 25 (23 unique) | 1,107,436 | *Microcystis*, *Anabaena*, *Phormidium*, *Kaarinaea* (title: *Nostoc*), *Planktothrix* |
| `cylindrospermopsin` | guanidine alkaloid | 17 | 317,922 | *Cylindrospermopsis*, *Raphidiopsis*, *Aphanizomenon*, *Kamptonema* (title: *Oscillatoria*) |
| `saxitoxin` | PSP alkaloid | 34 | 65,608 | *Dolichospermum* (title: *Anabaena*), *Cylindrospermopsis*, *Aphanizomenon*, *Microseira* (title: *Lyngbya*) |
| `anatoxin_a` | tropane alkaloid | 17 | 44,388 | *Cuspidothrix*, *Kamptonema* (title: *Oscillatoria*) |
| `nodularin` | cyclic pentapeptide | 6 | 120,982 | *Nostoc*, *Nodularia* |
| `guanitoxin` | phosphoguanidine ester | 1 | 13,572 | *Sphaerospermopsis* |
| `aetokthonotoxin` | bis-chlorinated indole alkaloid | 2 | 46,169 | *Aetokthonos* |

Each reference directory holds `<ref>.fna`, `<ref>.idx` (minimap2 index),
`<ref>.contigs.json`, `meta.json` and `source.txt`, in the formats defined
by `mapping_refs/PROTOCOL.md`. Five sequences are NCBI subranges fetched
with EFetch `seq_start`/`seq_stop`, which is why their FASTA headers look
like `DQ787200.1:1-30300`. Three cylindrospermopsin records were
hard-masked to N after download (see QA findings).

### BGC tier definitions

Every sequence is tiered from its NCBI record title (`hab_panel.tier_for_title`).
Two curated overrides are recorded in the `tier_note` column.

| Tier | Definition | Panel count (unique seqids) |
|---|---|---:|
| A | Complete BGC: the title names a gene cluster and not "cluster, partial" | 35 |
| B | Partial BGC or multi-gene fragment | 4 |
| C | Single gene, complete CDS | 16 |
| D | Single-gene fragment (partial cds / partial sequence) | 45 |

The overrides are AY210783.2, which source.txt calls the "complete nda
cluster (canonical)" although its title only lists the ndaA–ndaI CDSs, and
the guanitoxin CP080598.1 subrange, whose title is just "chromosome".

## Rebuilding the references

```bash
# Fetch everything from NCBI and build into ./build (gitignored).
python3 build_hab_refs.py

# Any output dir; a subset of references.
python3 build_hab_refs.py --outdir /path/to/out --refs microcystin nodularin

# Rebuild and verify against the live panel (read-only).
python3 build_hab_refs.py --outdir /path/to/out \
    --compare-to /data/scratch/refdbs/mapping_refs \
    --report data/reproducibility_check.tsv
```

For each reference, `build_hab_refs.py` does the following:

1. Fetches every sequence with EFetch (FASTA). Full records are fetched in
   batches and subranges one at a time. Each request is POSTed through
   `hab_panel.eutils_get`, which retries 6 times with 2/4/8/16/32 s
   backoff on network errors, HTTP 429/5xx, empty bodies, eutils error
   payloads, and responses with the wrong record set. NCBI eutils was
   timing out on the day the panel was built.
2. Applies the `mask_n` ranges and checks each sequence's length and SHA-1
   against `accessions.tsv`. If NCBI has changed a record, it prints a
   warning instead of failing silently.
3. Writes `<ref>.fna` in `accessions.tsv` order, keeping duplicates, with
   70-column wrapping like efetch.
4. Writes `<ref>.contigs.json`, byte-compatible with
   `project_zebra/build_contig_index.py`, and copies `meta.json` and
   `source.txt`.
5. Builds the index exactly as PROTOCOL.md Step 2 specifies:
   `minimap2 -t 4 -d <ref>.idx <ref>.fna`. The default binary is
   `/matika/apps/minimap2/minimap2`; override it with `--minimap2` or
   `$MINIMAP2`.
6. Writes `build_manifest.json`, which records the minimap2 version and
   the sha256 of every output.

The default `--outdir` is `./build`. The script **refuses** to write into
`/data/scratch/refdbs/mapping_refs` unless you pass `--allow-live`. To
update the live panel, build somewhere else, verify, and then copy the
reference directory across.

### Reproducibility check (2026-10-19)

All seven references were rebuilt from NCBI into a temporary directory
and compared against the live panel. Results are in
[`data/reproducibility_check.tsv`](data/reproducibility_check.tsv).

| Check | Result (all 7 references) |
|---|---|
| `.fna` sequences: same seqids, same sequences, order and wrapping ignored | PASS |
| `.fna` record order and header descriptions | identical |
| `contigs.json` | identical content, and byte-identical |
| `meta.json`, `source.txt` | byte-identical |
| `.idx` (minimap2 2.28-r1209) | byte-identical |

On the first pass, three cylindrospermopsin sequences did not match. The
live copies had been hard-masked to N, which `record_live_panel.py` now
captures in the `mask_n` column. After that change every reference
rebuilds byte-for-byte.

## Plugging into danaSeq

The pipeline does not care where a reference came from. It needs a
directory of `<name>/<name>.idx` subdirectories:

```bash
nextflow run /data/dev/danaSeq/nanopore_live/main.nf \
    --input        /path/to/run \
    --outdir       /data/scratch/nanopore_live/<run-id> \
    --mapping_refs /data/scratch/refdbs/mapping_refs \
    --run_db_integration -resume
```

- `main.nf` fans out across every subdirectory that contains
  `<name>/<name>.idx`. `MAP_REFERENCE` runs
  `minimap2 -ax map-ont --secondary=no`, and `bin/filter_minimap2.awk`
  then keeps hits with `mapq >= 1 && aligned_len >= 10`.
- `viz/preprocess/preprocess.py` `build_mapping_payloads` reads each
  reference's `meta.json` `category`. It normalises the value with
  `.split('-')[0].upper()` and accepts only `AIS`, `HAB` or `SAR`.
  **Use `"HAB"`, not `"HABs"`**: "HABs" becomes "HABS", is not recognised,
  and the reference falls back to the AIS view. With `HAB`, the payload is
  written as `data/hab_<ref>.json`, using `<ref>.contigs.json` for the
  position histogram.
- `viz/src/views/HABsView.svelte` `HAB_COMPOUNDS` lists the seven
  reference ids. Each `id` must match the reference directory, the `.idx`
  basename and the `hab_<id>.json` suffix. A new HAB reference needs a new
  entry there before it appears in the HABs tab.
- The standalone backfill route (`project_zebra/run-nanopore-mapping.sh`,
  then `load_mapping_to_duckdb.py`, then `preprocess_habs.py`) consumes the
  same files. `log_hab_backfill.txt` shows that it was started for
  all seven references on 2026-06-15.

## QA findings on the live panel

These findings are computed by the scripts here and listed per toxin on
the inventory page. None of them is fixed in the panel itself:
HABGDB records the panel as it is.

- **Duplicates (microcystin).** `MK870090.1` and `AY212249.1` each appear
  twice, and `AJ536156.1` has the same sequence as `AY212249.1` (both are
  *Anabaena* sp. 90). As a result, `contigs.json` contains duplicate
  contig names. In cylindrospermopsin, `HQ682085.1` is identical to
  `HQ682086.1`, and both are 71 bp. In saxitoxin, `EU629185.1` is identical
  to `EU629189.1`.
- **Near-identical sequences and MAPQ 0.** `panel_redundancy.py` finds
  groups at ≥99.5 % identity over ≥90 % of the shorter sequence. These
  include 11 *Microcystis* mcy clusters (AF183408 plus 10 OQ2910xx), 14
  of the 17 cyr sequences, 12 of the 14 *Cuspidothrix* anaC fragments, and
  single-gene sxt/ana fragments contained in the T3 / PCC 6506 clusters.
  minimap2 can give a read that falls in such a shared region MAPQ 0, and
  danaSeq's `mapq >= 1` filter then drops it, so per-producer resolution
  inside these groups is limited.
- **Undocumented masking (cylindrospermopsin).** In `KJ139745.1`,
  `EU140798.1` and `FJ418586.4`, 3,786 bp are N in the panel but bases at
  NCBI. The masked ranges are exactly the annotated transposase CDSs
  ([`data/transposase_features.tsv`](data/transposase_features.tsv)).
  Neither `meta.json` nor `source.txt` mentions this.
- **Subranges that drop annotated genes**
  ([`data/subrange_features.tsv`](data/subrange_features.tsv)):
  - `DQ787200.1:1-30300` (saxitoxin, *C. raciborskii* T3) stops before the
    transposase at 30,612. It therefore leaves out **sxtY and sxtZ**
    (32,612–34,677) and an OmpR-family regulator of the 37,606 bp record.
  - `FJ477836.2:5100-33100` (anatoxin-a) leaves out **anaH** (33,188–33,658,
    annotated as a transposase), along with flanking transposases.
  - `GQ385961.1` and `AF183408.1` trim only flanking genes.
  - The guanitoxin window matches source.txt exactly: gntJ starts at
    3,534,866 and gntT ends at 3,547,437, each with 500 bp flanks.

  Every transposase CDS in the panel is either outside a subrange or
  masked, except the atn1 transposase in AY212249.1 / AJ536156.1 and one
  in AY210783.2.
- **Stale `meta.json` text.** Each of these disagrees with the panel or
  with its own `source.txt`:
  - cylindrospermopsin says "33-sequence" (the panel has 17).
  - saxitoxin says "35 sequences" (34) and mentions soil metagenome contig
    `CBFDRI010000008`, which is not in the `.fna`.
  - nodularin gives its source as "AY210783 + AF227755" and "2 complete
    nda cluster sequences". The panel has 6 sequences, and source.txt
    excludes AF227755 as a *Sinorhizobium* 16S record.
  - guanitoxin still describes the ~5 Mb chromosome, although the panel
    holds the 13,572 bp gnt region.
- **Renamed taxa.** NCBI Taxonomy has since moved several producers:
  *Anabaena circinalis* → *Dolichospermum circinale*, *Lyngbya wollei* →
  *Microseira wollei*, *Oscillatoria* sp. PCC 6506 → *Kamptonema* sp.
  PCC 6506, and *Nostoc* sp. 152 → *Kaarinaea lacus* PCC 9237. The FASTA
  titles keep the old names. The `OQ291093.1` title says FACHB-1388, while
  the NCBI organism is FACHB-1338.
- **Aetokthonotoxin.** `meta.json` notes that no *complete* genome
  exists. That is still true at Complete Genome level, but scaffold-level
  *A. hydrillicola* Thurmond2011 assemblies are in NCBI
  (e.g. GCF_017591595.2, released 2023-09-18).

### Data files

| File | Source | What it is |
|---|---|---|
| [`data/ncbi_esummary.json`](data/ncbi_esummary.json) | NCBI Entrez esummary (nuccore), retrieved 2026-10-19 | Trimmed docsums for all 100 panel accessions: title, organism, taxid, slen, create/update dates |
| [`data/live_panel_checksums.tsv`](data/live_panel_checksums.tsv) | `/data/scratch/refdbs/mapping_refs/<ref>/` | Size, mtime and sha256 of every live file at snapshot time, as provenance for the 2026-06-11 build |
| [`data/subrange_features.tsv`](data/subrange_features.tsv) | NCBI feature tables (`rettype=ft`) | Annotated genes and CDSs of every subranged record, marked in / OUT / CLIP relative to the extracted window |
| [`data/transposase_features.tsv`](data/transposase_features.tsv) | NCBI feature tables | Every transposase CDS in a panel record, and whether the panel keeps, trims or masks it |
| [`data/redundancy_pairs.tsv`](data/redundancy_pairs.tsv) | `panel_redundancy.py` (minimap2 asm20 self-alignment) | Sequence pairs at ≥99.5 % identity over ≥90 % of the shorter sequence, within each reference |
| [`data/reproducibility_check.tsv`](data/reproducibility_check.tsv) | `build_hab_refs.py --compare-to` | Per-reference results of rebuild vs live |

## Re-running

```bash
# Rebuild and verify (network: NCBI eutils).
python3 build_hab_refs.py --outdir build --compare-to /data/scratch/refdbs/mapping_refs \
    --report data/reproducibility_check.tsv
# Redundancy table (needs the build).
python3 panel_redundancy.py --build-dir build
# Inventory pages (offline).
python3 hab_inventory.py
```

`record_live_panel.py` is the tool that created `refs/` from the live
panel. Re-running it **overwrites** the committed accession lists with
whatever is currently live. To change the panel, edit
`refs/<ref>/accessions.tsv` (and `meta.json` / `source.txt`), rebuild,
and commit.

## Related projects

- **[AISGDB](https://github.com/rec3141/AISGDB)**: aquatic invasive
  species inventory (GLANSIS + federal SOR/2015-121).
- **[SARGDB](https://github.com/rec3141/SARGDB)**: species-at-risk
  inventory (federal SARA + MB ESEA).
- **[danaSeq](https://github.com/rec3141/danaSeq)**: the nanopore
  pipeline that consumes these reference directories. Its HABs view
  (`/<run>/#habs`) shows per-sample hit counts, identity histograms and
  per-producer position histograms for each toxin.
- **[microscape.app](https://microscape.app)**: the SPA frontend hosting
  the AIS, SAR and HAB views.

## License

CC BY 4.0 for documents; MIT for code. See individual files for headers.
