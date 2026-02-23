# fl.AI-CLI

Standalone Fly Gene Classification pipeline that converts input gene symbols to FlyBase IDs, finds supporting literature, summarizes evidence, and classifies each gene into user-defined categories.

## Pipeline Flowchart

![Fly Gene Classification Pipeline](fly_gene_pipeline_high_level_flowchart.png)

## High-Level Flow (8 bullets)

- Read every input CSV in a directory and extract gene symbols from a configurable column.
- Convert symbols to `FBgn` identifiers using the project helper script.
- Validate genes and load FlyBase synonym/reference tables.
- Collect candidate references from FlyBase, PubMed, and Europe PMC.
- Merge and deduplicate references, then rank by source consensus and recency.
- Fetch metadata and keep only references whose title/abstract matches your keywords.
- Pull full text when possible (with DOI/PMC fallbacks), then summarize gene-specific evidence.
- Aggregate per-gene summaries, classify into keyword categories, and export one Excel per input CSV.

## Example Flowchart Walkthrough (dummy numbers for gene X)

Assume one gene (`gene X`, `FBgn0009999`) and keywords:
`["sleep", "circadian", "synapse"]`

1. **Reference collection**
   - FlyBase returns 9 PMCIDs.
   - PubMed search returns 14 PMCIDs.
   - Europe PMC search returns 11 PMCIDs.
   - Combined raw total = 34 references.

2. **Merge + deduplicate + rank**
   - 34 raw references collapse to 21 unique PMCIDs.
   - Priority is higher when a PMCID is found by multiple sources (e.g., FlyBase + PubMed), then by newer publication year.
   - With `--reference-limit 15`, only top 15 proceed.

3. **Keyword filtering (title/abstract)**
   - Of 15 references, 10 contain at least one keyword match.
   - 5 references are dropped for no relevant keyword evidence.

4. **Full-text retrieval + summarization**
   - Full text is retrieved for 7 of 10 references.
   - 3 references fall back to title+abstract summaries.
   - Quality control keeps 6 high-quality summaries (non-"No evidence found").

5. **Aggregation + classification**
   - The 6 summaries are concatenated into one evidence block for `gene X`.
   - Model output example:
     - Category: `sleep; circadian`
     - Confidence: `84`
     - Rationale: short evidence-backed explanation from aggregated summaries.

6. **Export**
   - Output file: `<input_name>_classification.xlsx`
   - Sheets:
     - `Gene Set`
     - `Classification`
     - `Reference Summaries`

## Requirements

- Python 3.10+ recommended
- Dependencies in `requirements.txt`
- Environment variables:
  - `OPENAI_API_KEY` (required for summarization/classification)
  - `NCBI_API_KEY` (optional, improves NCBI throughput)
  - `UNPAYWALL_TOKEN` (optional, defaults to configured email)

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run

```bash
python flai-gene-classification.py /path/to/input_dir \
  --keywords "sleep,circadian,synapse" \
  --reference-limit 500 \
  --input-gene-col ext_gene
```

Optional:

- Add `--force-all` to ignore batch checkpoints and recompute all genes.

