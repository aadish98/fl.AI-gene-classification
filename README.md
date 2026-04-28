# fl.AI Gene Classification

Standalone gene classification pipeline that converts input fly gene symbols to FlyBase IDs, optionally maps them to human or mouse orthologs, finds species-scoped supporting literature, summarizes evidence, and classifies each gene into user-defined categories.

## Documentation

All non-code documentation lives under [`Docs/`](Docs/):

- [`Docs/Decks/fl.AI_Data_Flow_Deck.pptx`](Docs/Decks/fl.AI_Data_Flow_Deck.pptx) — non-technical walkthrough of how data moves through the pipeline, including per-organism (fly / human / mouse) flowchart slides.
- [`Docs/Figures/reference_selection_flowchart.png`](Docs/Figures/reference_selection_flowchart.png) — visual of how candidate references are gathered, ranked, and filtered.
- [`Docs/Scripts/`](Docs/Scripts/) — scripts that regenerate the documentation assets above.

## High-Level Flow

- Read every input CSV in a directory and extract fly gene symbols from a configurable column.
- Convert symbols to `FBgn` identifiers using the project helper script.
- Optional: map fly genes to human or mouse orthologs with DIOPT.
- Build a species-specific gene set (`fly`, `human`, or `mouse`) with canonical symbols, IDs, and synonyms.
- Collect candidate references from species-specific sources.
- Merge and deduplicate references, then rank by source consensus and recency.
- Fetch metadata and keep references whose title/abstract matches your keywords.
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
  - `OPENAI_SUMMARY_MODEL` (optional, defaults to `gpt-5.4-nano`)
  - `OPENAI_CLASSIFICATION_MODEL` (optional, defaults to `gpt-5.4`)
  - `OPENAI_SUMMARY_REASONING_EFFORT` (optional, defaults to `low`)
  - `OPENAI_CLASSIFICATION_REASONING_EFFORT` (optional, defaults to `high`)
  - `NCBI_API_KEY` (optional, improves NCBI throughput)
  - `UNPAYWALL_TOKEN` (optional, defaults to configured email)

## Data Dependencies

This script expects access to shared project data on the lab storage mount:

- macOS: `/Volumes/umms-rallada`
- ARC/HPC Linux: `/nfs/turbo/umms-rallada`

It uses shared FlyBase reference files and shared PubMed/full-text caches under that mount. Human and mouse modes also use cached HGNC/MGI, NCBI, GeneRIF, UniProt, and DIOPT data.

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run: Fly Literature

```bash
python flai_gene_classification.py /path/to/input_dir \
  --keywords "sleep,circadian,synapse" \
  --reference-limit 500 \
  --input-gene-col ext_gene
```

This preserves the original fly workflow: input fly symbols -> `FBgn` IDs -> fly literature -> `<input_name>_classification.xlsx`.

## Run: Fly to Human Orthologs

```bash
python flai_gene_classification.py /path/to/input_dir \
  --orthologs human \
  --keywords "sleep,circadian,synapse" \
  --reference-limit 500 \
  --input-gene-col ext_gene
```

The pipeline writes sibling `<input_name>_human.csv` files containing the fly-to-human DIOPT crosswalk, then classifies each unique human ortholog and exports `<input_name>_human_classification.xlsx`.

## Run: Fly to Mouse Orthologs

```bash
python flai_gene_classification.py /path/to/input_dir \
  --orthologs mouse \
  --keywords "sleep,circadian,synapse" \
  --reference-limit 500 \
  --input-gene-col ext_gene
```

The pipeline writes sibling `<input_name>_mouse.csv` files containing the fly-to-mouse DIOPT crosswalk, then classifies each unique mouse ortholog and exports `<input_name>_mouse_classification.xlsx`.

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
- Europe PMC search scoped to Drosophila / fruit fly

Human mode uses:

- HGNC complete set for symbols/synonyms and Entrez/HGNC IDs
- NCBI `gene2pubmed` filtered to `tax_id == 9606`
- NCBI GeneRIF filtered to `tax_id == 9606`
- UniProt Swiss-Prot reviewed human entries
- PubMed and Europe PMC searches scoped to Homo sapiens

Mouse mode uses:

- MGI marker reports for symbols/synonyms and Entrez/MGI IDs
- NCBI `gene2pubmed` filtered to `tax_id == 10090`
- NCBI GeneRIF filtered to `tax_id == 10090`
- UniProt Swiss-Prot reviewed mouse entries
- PubMed and Europe PMC searches scoped to Mus musculus

## Output Workbook

Each output workbook includes:

- `Gene Set`: the full input gene set and, for ortholog modes, the fly-to-human or fly-to-mouse crosswalk. Unmapped fly genes are included with `Status = unmapped`.
- `Classification`: one row per classified species-specific gene.
- `Reference Summaries`: high-quality reference summaries with source labels and fly traceability columns for ortholog modes.

Optional:

- Add `--force-all` to ignore batch checkpoints and recompute all genes.
- Add `--refresh-human-data` or `--refresh-mouse-data` to re-download species data caches.

## Input Expectations

- The input path should be a directory containing one or more CSV files.
- Each CSV is treated as an independent gene set.
- The gene symbol column defaults to `ext_gene`, but can be changed with `--input-gene-col`.
- Output workbooks are written next to each input file as `<input_name>_classification.xlsx`, `<input_name>_human_classification.xlsx`, or `<input_name>_mouse_classification.xlsx`.

