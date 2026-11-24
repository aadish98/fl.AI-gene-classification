# fl.AI Gene Classification

Standalone gene classification pipeline that classifies fly, human, or mouse gene sets, optionally maps fly genes to human or mouse orthologs, finds species-scoped supporting literature, summarizes evidence, and classifies each gene into user-defined categories.

## Documentation

All non-code documentation lives under [`Docs/`](Docs/):

- [`Docs/Decks/fl.AI_Data_Flow_Deck.pptx`](Docs/Decks/fl.AI_Data_Flow_Deck.pptx) — non-technical walkthrough of how data moves through the pipeline, including per-organism (fly / human / mouse) flowchart slides.
- [`Docs/Figures/reference_selection_flowchart.png`](Docs/Figures/reference_selection_flowchart.png) — visual of how candidate references are gathered, ranked, and filtered.
- [`Docs/Scripts/`](Docs/Scripts/) — scripts that regenerate the documentation assets above.

## High-Level Flow

- Read every input CSV in a directory and extract fly gene symbols from a configurable column.
- Native human and mouse runs can read symbols directly and resolve them through HGNC or MGI without FlyBase or DIOPT.
- Convert symbols to `FBgn` identifiers using the project helper script.
- Optional: map fly genes to human or mouse orthologs with DIOPT.
- Build a species-specific gene set (`fly`, `human`, or `mouse`) with canonical symbols, IDs, and synonyms.
- Collect candidate references from species-specific sources.
- Merge and deduplicate references, then rank by source consensus and recency.
- Fetch metadata and keep references whose title/abstract matches your keywords.
- Generate organism-strict LLM summaries and classifications: fly, human, and mouse runs each restrict evidence to the selected organism's gene and skip references that only describe another organism's ortholog/paralog/homolog.
- Pull full text when possible, summarize gene-specific evidence, classify into keyword categories, and export one Excel per input CSV.

## Example End-To-End Flow

Assume one gene (`gene X`, `FBgn0009999`) and keywords `sleep,circadian,synapse`.

1. Reference collection gathers candidate papers from FlyBase, PubMed, and Europe PMC.
2. Merge and ranking collapse duplicates and prioritize papers seen in multiple sources.
3. Keyword filtering removes papers with no relevant keyword evidence in title or abstract.
4. Full-text retrieval and summarization produce gene-specific function and phenotype summaries.
5. Aggregation produces a final classification plus detailed paper-level exports.

## Requirements

- Python 3.10+ recommended
- Dependencies in `requirements.txt`
- Environment variables:
  - `OPENAI_API_KEY` (required for summarization/classification)
  - `OPENAI_SUMMARY_MODEL` (optional, defaults to `gpt-6-luna`)
  - `OPENAI_CLASSIFICATION_MODEL` (optional, defaults to `gpt-6.1-sol`)
  - `OPENAI_SUMMARY_REASONING_EFFORT` (optional, defaults to `low`)
  - `OPENAI_CLASSIFICATION_REASONING_EFFORT` (optional, defaults to `high`)
  - `NCBI_API_KEY` (optional, improves NCBI throughput)
  - `UNPAYWALL_TOKEN` (optional, defaults to configured email)

## Data Dependencies

This script expects access to shared project data on the lab storage mount:

- macOS: `/Volumes/umms-rallada`
- ARC/HPC Linux: `/nfs/turbo/umms-rallada`

It uses shared FlyBase reference files and shared PubMed/full-text caches under that mount. The FlyBase Chado XML download and its citation index live in the repo-local, gitignored `Data/FlyBase_Chado/` folder. Human and mouse modes also use cached HGNC/MGI, NCBI, GeneRIF, UniProt, and DIOPT data.

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run: Prepared, Parallel-Safe Workflow

```bash
python flai_gene_classification.py prepare /path/to/input_dir \
  --keywords "sleep,circadian,synapse" \
  --reference-limit 500 \
  --input-gene-col ext_gene
```

This single-writer preparation step converts input fly symbols to `FBgn` IDs and validates the source CSVs. After preparation, classification jobs can skip source CSV mutation.

### How Fly Symbols Become FBgn IDs

`GetFBgnIDs.py` (also run by the preparation step) maps each symbol **exactly as written**; it never edits your gene column and never calls OpenAI or any web API. Matching is case-sensitive, in this order:

1. The current symbol in the FlyBase synonym table.
2. A symbol synonym in that table (`mapping_note`: `matched via synonym of <symbol>`).
3. If several genes share the synonym (e.g. `Trh`), the gene whose symbol FlyBase publications use most, **only when the lead is clear** (at least 3 publications and 2x the runner-up). The note lists the counts, e.g. `synonym shared by 3 genes; FlyBase citations: Trhn 41, trh 6, Hn 2 → Trhn`. Otherwise `flybase_gene_id` is `-` and the note lists every candidate: `ambiguous synonym: ...`.
4. Only if nothing matched, Greek letters are spelled out in a lookup key (`14-3-3ζ` -> `14-3-3zeta`) because FlyBase's table spells them out; the note says so.

Step 3 reads a local citation index built once from the FlyBase Chado XML gene file (about 6 GB download, stored under `Data/FlyBase_Chado/`, which is gitignored):

```bash
python HelperScripts/build_synonym_evidence.py
```

Without the index, shared synonyms stay unmapped and the note says to build it. The thresholds are `MIN_WINNER_PUBS` and `WIN_RATIO` in `HelperScripts/flybase_synonym_evidence.py`.

### Fly Literature

```bash
python flai_gene_classification.py classify /path/to/input_dir \
  --species fly \
  --keywords "sleep,circadian,synapse" \
  --reference-limit 500
```

This classifies prepared fly inputs and writes `<input_name>_classification.xlsx`.

### Fly to Human Orthologs

```bash
python flai_gene_classification.py classify /path/to/input_dir \
  --species human \
  --keywords "sleep,circadian,synapse" \
  --reference-limit 500
```

The pipeline writes `Human/<input_name>_human.csv` files containing the fly-to-human DIOPT crosswalk, then classifies each unique human ortholog and exports `Human/<input_name>_human_classification.xlsx`.

### Native Human Genes

For a CSV containing human symbols, use matching input and target species. The native path does not run FlyBase conversion, DIOPT, or emit fly-trace columns:

```bash
python flai_gene_classification.py prepare /path/to/human_input \
  --species human \
  --input-species human \
  --input-gene-col human_gene_symbol

python flai_gene_classification.py classify /path/to/human_input \
  --species human \
  --input-species human \
  --keywords "sleep,circadian,Alzheimer" \
  --reference-limit 300
```

The native human input must contain `human_entrez_gene_id` after `prepare` (or provide a resolvable symbol column to `prepare`). Native output is written as `<input_name>_human_classification.xlsx`.

### Fly to Mouse Orthologs

```bash
python flai_gene_classification.py classify /path/to/input_dir \
  --species mouse \
  --keywords "sleep,circadian,synapse" \
  --reference-limit 500
```

The pipeline writes `Mouse/<input_name>_mouse.csv` files containing the fly-to-mouse DIOPT crosswalk, then classifies each unique mouse ortholog and exports `Mouse/<input_name>_mouse_classification.xlsx`.

### Estimate and Cache Maintenance

```bash
python flai_gene_classification.py estimate /path/to/input_dir --species human
python flai_gene_classification.py cache validate --apply
```

The PubMed metadata cache now uses SQLite. It prefers `FLAI_PUBMED_CACHE_DB`, then the shared Turbo PubMed cache directory, then a local gitignored `.flai_cache/pubmed_cache.sqlite` fallback.

## Ortholog Defaults

- DIOPT output species:
  - human: `9606`
  - mouse: `10090`
- Default DIOPT filter: `exclude_low_score_2`, corresponding to DIOPT's "score > 2, unless only match score is 1 or 2" filter.
- All surviving orthologs are retained. One fly gene can produce multiple human or mouse rows.
- Fly traceability columns are carried into target-species outputs: input fly symbol, primary FlyBase symbol, FlyBase ID, DIOPT score, best-match flag, and supporting algorithms.

## Reference Sources

Fly mode uses:

- FlyBase local reference tables
- PubMed search scoped to Drosophila / fruit fly
- Europe PMC search scoped to `Drosophila melanogaster` via structured organism/MeSH filters

Human mode uses:

- HGNC complete set for symbols/synonyms and Entrez/HGNC IDs
- NCBI `gene2pubmed` filtered to `tax_id == 9606`
- NCBI GeneRIF filtered to `tax_id == 9606`
- UniProt Swiss-Prot reviewed human entries
- PubMed and Europe PMC searches scoped to `Homo sapiens`

Mouse mode uses:

- MGI marker reports for symbols/synonyms and Entrez/MGI IDs
- NCBI `gene2pubmed` filtered to `tax_id == 10090`
- NCBI GeneRIF filtered to `tax_id == 10090`
- UniProt Swiss-Prot reviewed mouse entries
- PubMed and Europe PMC searches scoped to `Mus musculus`

## Output Workbook

Each output workbook includes:

- `Gene Set`: the full input gene set and, for ortholog modes, the fly-to-human or fly-to-mouse crosswalk. Unmapped fly genes are included with `Status = unmapped`.
- `Classification`: one row per classified species-specific gene.
- `Reference Summaries`: high-quality reference summaries with source labels; fly traceability columns are included only for fly-to-human or fly-to-mouse ortholog modes.

Optional:

- Add `--force-all` to ignore batch checkpoints and recompute all genes.
- Add `--refresh-human-data` or `--refresh-mouse-data` to re-download species data caches.

## Input Expectations

- The input path should be a directory containing one or more CSV files.
- Each CSV is treated as an independent gene set.
- Fly input symbols default to `ext_gene`; native human and mouse runs can use a configured symbol column with `--input-gene-col`.
- Output workbooks are written next to each input file as `<input_name>_classification.xlsx`, `<input_name>_human_classification.xlsx`, or `<input_name>_mouse_classification.xlsx`.

