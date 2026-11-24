# Caching Cheatsheet

A single-page reference for every cache the fl.AI gene-classification pipeline
reads from or writes to. Use it to figure out **where data lives**, **when it
gets reused**, and **how to invalidate it**.

---

## TL;DR

| Layer | What it stores | Where it lives | Survives across runs? | Bust it with |
| --- | --- | --- | --- | --- |
| Run store | Per-gene final results (hits, summaries, OpenAI usage) | `<input_dir>/.batch_state/run_store_<species>.json` | Yes | `--force-all` or `CACHE_SCHEMA_VERSION` bump |
| Batch checkpoints | Per-CSV progress + per-batch artifacts | `<input_dir>/.batch_state/<species>_<csv_stem>_<hash>/` | Yes | `--force-all`, input fingerprint change, or schema bump |
| PubMed metadata + full-text method | Title/abstract/year/DOI per PMID, fetched-fulltext method per PMID | SQLite `pubmed_cache.sqlite` (Turbo or `.flai_cache/`) | Yes, lab-shared | `cache validate --apply`, delete DB |
| DIOPT responses | Raw JSON per `(FBgn, taxon, filter)` | `--diopt-cache-dir` ↦ `$FLAI_CACHE_DIR/diopt` ↦ `~/.cache/flai-gene-classification/diopt` | Yes | Delete files; errors auto-expire (TTL) |
| Species reference data | HGNC/MGI tables, `gene2pubmed.gz`, `generifs_basic.gz`, UniProt TSV | `--{human,mouse}-data-dir` ↦ `$FLAI_CACHE_DIR/{human,mouse}-data` ↦ `~/.cache/flai-gene-classification/{human,mouse}-data` | Yes | `--refresh-human-data` / `--refresh-mouse-data` |
| FlyBase data | `fb_synonym*`, `entity_publication*`, `fbrf_pmid_pmcid_doi*` TSVs | `--flybase-data-dir` ↦ `$FLYBASE_DATA_DIR` ↦ Turbo `Data/FlyBase` ↦ repo `Data/FlyBase` | Yes, lab-shared | Re-download (handled by `flybase_data.py` bootstrap) |
| FlyBase Chado XML + citation index | `chado_FBgn.xml.gz` (about 6 GB) and `synonym_evidence_<release>.tsv.gz` (publications per gene per symbol) | repo `Data/FlyBase_Chado/` | No (gitignored) | `python HelperScripts/build_synonym_evidence.py` (resumes the download, rebuilds the index) |
| OpenAI pricing | Model price table (incl. `cached_input_per_1m`) | `Data/openai_pricing.csv` | Yes (committed) | Soft-run auto-refreshes from LiteLLM JSON |
| Soft-run profiles | Empirical token averages + per-species keyword pass rates | `.flai_system/soft_run_reference_profiles.json` | Yes (committed) | Species pass rates auto-update when a real run sees more references than the stored calibration; token-average profiles are rebuilt with `HelperScripts/calibrate_soft_run_profile.py` |
| In-memory caches | Per-process `lru_cache`s + dict caches around the SQLite DB | RAM only | No | Process exit |
| OpenAI prompt cache | Server-side prefix cache (counted in `cached_input_tokens`) | OpenAI side, not local | Yes (OpenAI-managed) | Not user-controlled |

---

## 1. Run store — per-gene result cache

- **File**: `<input_directory>/.batch_state/run_store_<species>.json` (legacy
  `run_store.json` is read once for migration).
- **Schema version**: `CACHE_SCHEMA_VERSION = 3`. Stores load empty if the
  on-disk version doesn't match.
- **Key per gene**:
  `v<schema>|<species>|<gene_id>|<reference_limit>|<sha256(keywords)[:16]>`
  (`_make_gene_store_key`). Changing keywords, reference limit, or species
  produces a new key and forces reprocessing.
- **What's stored per gene**: `gene_record`, `hit_info`, summaries (with chunk
  counts/word counts), and `openai_usage` (`requests` + `totals`, including
  `cached_input_tokens` and `reasoning_output_tokens`).
- **Concurrency**: writes are guarded by `flock` on
  `run_store_<species>.json.lock`; the writer reloads the disk store and merges
  before saving so parallel jobs on the same input directory don't clobber each
  other.
- **Reuse rule**: a batch reuses a gene's record iff its store key already
  exists AND `--force-all` is not set.

## 2. Batch checkpoints — per-CSV resume state

- **Directory**:
  `<input_directory>/.batch_state/<species>_<csv_stem>_<sha1(species|abs_csv_path)[:10]>/`
- **Files**:
  - `state.json` — `version`, `csv_path`, `input_fingerprint`,
    `total_batches`, `completed_batches`, `batch_files`, processed/reused
    counters.
  - `batch_<index>.json` — full per-batch payload (records used to rebuild the
    workbook).
- **Input fingerprint** (`_make_input_fingerprint`): sha256 over
  `{schema_version, species, sorted gene_ids, normalized keywords,
  reference_limit, diopt_filter}`. Any change resets that CSV's batch state.
- **Resume rule**: completed batch indices are skipped on the next run unless
  `--force-all` is passed.

## 3. PubMed metadata + full-text method (SQLite)

- **DB resolution** (`HelperScripts/pubmed_cache_db.py`,
  `resolve_pubmed_cache_db_path`):
  1. Explicit path argument
  2. `FLAI_PUBMED_CACHE_DB`
  3. Lab Turbo: `<resolve_pubmed_cache_dir()>/pubmed_cache.sqlite` (i.e.
     `${FLAI_PUBMED_CACHE_DIR or PUBMED_CACHE_DIR}` ↦ `<FLAI_DATA_DIR>/PubMed Cache`
     ↦ Turbo `umms-rallada/UM Lab Users/Aadish/Data/PubMed Cache`)
  4. Local fallback: `.flai_cache/pubmed_cache.sqlite`
- **Tables**: `pubmed_metadata`, `fulltext_method`, `cache_meta`
  (schema version 1).
- **Writes**: serialized via an `fcntl`-locked sentinel file
  (`*.sqlite.lock`); SQLite uses `busy_timeout=60s`, `synchronous=NORMAL`.
- **Upserts**: only overwrite a column when the incoming value is non-empty
  (so partial fetches don't blank-out completer rows).
- **In-process wrappers**: `_pmid_title_abstract_cache` and
  `_fulltext_method_cache` mirror lookups for the lifetime of the process.
- **Validate / migrate**:
  ```bash
  python flai_gene_classification.py cache validate --apply
  ```
  Runs `HelperScripts/validate_pubmed_cache.py` to ingest legacy CSV caches and
  PRAGMA-integrity-check the SQLite DB.

## 4. DIOPT response cache

- **Location**: `--diopt-cache-dir` ↦ `$FLAI_DIOPT_CACHE_DIR` ↦
  `$FLAI_CACHE_DIR/diopt` ↦ `~/.cache/flai-gene-classification/diopt`.
- **Filename pattern**: `<FBgn>_<output_taxon>_<filter>.json`.
- **What's cached**: the first JSON response from any of DIOPT's candidate
  endpoints (the client also remembers the winning endpoint per
  `(taxon, filter)` for the rest of the process).
- **Error caching**: failed lookups are written as
  `{"error": ..., "cached_at": <epoch>}` and re-tried after
  `FLAI_DIOPT_ERROR_CACHE_TTL_SECONDS` (default 6 hours; set to `0` to disable).
- **Atomic writes**: temp file + `replace()` to avoid partial JSON.

## 5. Human / mouse reference data caches

- **Location**:
  `--{human,mouse}-data-dir` ↦ `$FLAI_CACHE_DIR/{human,mouse}-data` ↦
  `~/.cache/flai-gene-classification/{human,mouse}-data`.
- **Files downloaded on first use** (via
  `species_data_utils.download_if_needed`):
  - Human: `hgnc_complete_set.txt`, `gene2pubmed.gz`, `generifs_basic.gz`,
    `uniprot_sprot_human.tsv`.
  - Mouse: MGI marker report + the equivalent NCBI/UniProt files.
- **Reuse rule**: any existing file is kept unless `refresh=True`. Use
  `--refresh-human-data` / `--refresh-mouse-data` (or
  `human_gene_data.configure(refresh=True)` /
  `mouse_gene_data.configure(refresh=True)`) to re-download.
- **In-process**: `load_gene_catalog` is `@lru_cache(maxsize=1)`; the configure
  hooks call `cache_clear()` when the data dir changes.

## 6. FlyBase data store

- **Location** (`HelperScripts/flybase_data.py`):
  `--flybase-data-dir` ↦ `$FLYBASE_DATA_DIR` ↦ `<FLAI_DATA_DIR>/FlyBase` ↦
  Turbo `umms-rallada/UM Lab Users/Aadish/Data/FlyBase` ↦ repo `Data/FlyBase`.
- **Required files**: `fb_synonym*`, `entity_publication*`,
  `fbrf_pmid_pmcid_doi*` (TSV/TSV.GZ).
- **Bootstrap**: missing files are downloaded from the pinned FlyBase release
  URLs by `_download_file()`.

### FlyBase Chado XML and citation index

- **Location**: `Data/FlyBase_Chado/` in the repo (gitignored): `chado_FBgn.xml.gz`,
  `chado_FBgn.dtd`, and `synonym_evidence_<release>.tsv.gz`.
- **Purpose**: `GetFBgnIDs.py` uses the index only to resolve a symbol shared by
  several genes (counts distinct FlyBase publications per gene per symbol).
  No network or LLM call at mapping time.
- **Build**: `python HelperScripts/build_synonym_evidence.py` downloads the XML
  (resumable) and streams it into the index; run it again to refresh for a new
  release (change `FLYBASE_RELEASE` in that file). The older
  `.flai_system/fbgn_collision_api_cache.json` is no longer used.

## 7. OpenAI pricing CSV

- **Path**: `Data/openai_pricing.csv` (committed seed file).
- **Columns**: `model, input_per_1m, cached_input_per_1m, output_per_1m,
  source, fetched_at_utc`.
- **Refresh**: `_refresh_openai_pricing_csv()` pulls
  `LITELLM_PRICING_URL` and rewrites the CSV; soft-run uses this when invoked
  with the refresh path. Cached-input pricing comes from LiteLLM's
  `cache_read_input_token_cost` (or `input_cost_per_token_batches` as
  fallback).

## 8. Soft-run reference profiles

- **Path**: `.flai_system/soft_run_reference_profiles.json` (committed).
- **Contents**:
  - `profiles[]`: empirical token averages (input/output per reference, per
    gene, requests per reference, etc.). A profile keyed on a
    `(summary_model, classification_model)` pair is used when it matches the
    configured models; otherwise the model-independent `"default": true`
    profile is used and repriced at the current models' rates.
  - `species_keyword_pass_rates`: per-species `keyword_pass_rate` plus the
    sample size (`calibration_limited_references`,
    `calibration_metadata_keyword_matches`) and source CSV.
- **Rebuilding profiles**: `HelperScripts/calibrate_soft_run_profile.py` pools
  per-request usage from saved run stores (`--dry-run` to preview, `--default`
  to force the model-independent profile). It leaves species pass rates alone.
- **Update rule**: `_maybe_update_species_keyword_pass_rate` only overwrites a
  species entry if a real run observed at least as many limited references as
  the stored calibration.

## 9. In-memory caches (per process)

- `_pmid_title_abstract_cache` and `_fulltext_method_cache` — read-through
  caches in front of the PubMed SQLite DB.
- `_meta_cache` — PMC metadata results in `metadata_resolver` flows.
- `@lru_cache(maxsize=1)` decorators in `flai_gene_classification.py`
  (FlyBase loaders, soft-run system data) and in
  `HelperScripts/{human,mouse}_gene_data.py` (`load_gene_catalog`).
- All of these die with the process — they don't persist to disk.

## 10. OpenAI prompt cache (server-side)

- The OpenAI Responses API auto-caches prompt prefixes; cached input tokens
  are reported in `usage.input_tokens_details.cached_tokens` and aggregated
  into the run-store totals as `cached_input_tokens`.
- The `cached_input_per_1m` column in `Data/openai_pricing.csv` is what
  cost reporting uses for those tokens.
- This cache is OpenAI-managed; we don't control TTL or invalidation.

---

## Environment variables (override defaults)

| Variable | Effect |
| --- | --- |
| `FLAI_DATA_DIR` | Project data root (used by FlyBase + PubMed cache fallbacks). |
| `FLAI_CACHE_DIR` | Root for human/mouse data and DIOPT cache when no explicit dir is given. |
| `FLAI_PUBMED_CACHE_DIR` / `PUBMED_CACHE_DIR` | Directory containing the PubMed SQLite DB. |
| `FLAI_PUBMED_CACHE_DB` | Explicit path to the SQLite file (highest priority after CLI). |
| `FLYBASE_DATA_DIR` | Explicit FlyBase TSV directory. |
| `FLAI_DIOPT_CACHE_DIR` | DIOPT JSON cache directory (also honored via `--diopt-cache-dir`). |
| `FLAI_DIOPT_ERROR_CACHE_TTL_SECONDS` | TTL for cached DIOPT error payloads (default `21600`). |

## Invalidation cookbook

- **Reprocess one input directory from scratch**: `--force-all`.
- **Reprocess one CSV only**: change keywords, `--reference-limit`, species, or
  `--diopt-filter` (any of these flips its `input_fingerprint`); or delete its
  `<input_dir>/.batch_state/<species>_<csv_stem>_<hash>/` folder.
- **Discard all gene-level cached results**: bump `CACHE_SCHEMA_VERSION` (in
  `flai_gene_classification.py`); the next load returns the empty default.
- **Force fresh DIOPT lookups**: delete the `.json` files in the DIOPT cache
  dir (errors auto-expire after the TTL).
- **Refresh species data**: `--refresh-human-data` / `--refresh-mouse-data`.
- **Audit / repair PubMed cache**:
  `python flai_gene_classification.py cache validate --apply`.
- **Refresh OpenAI pricing**: re-run any soft-run, or rerun
  `_refresh_openai_pricing_csv()` directly.
