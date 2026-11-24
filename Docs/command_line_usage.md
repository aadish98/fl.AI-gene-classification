# Command Line Usage

## Cheat Sheet

Run the pipeline from the repository root:

```bash
python flai_gene_classification.py <input_directory> [options]
```

The current subcommand form is:

```bash
python flai_gene_classification.py classify <input_directory> --species fly|human|mouse [options]
python flai_gene_classification.py estimate <input_directory> --species fly|human|mouse [options]
```

| Argument | Default | Purpose |
| --- | --- | --- |
| `input_directory` | required | Directory containing input CSV files. |
| `--keywords`, `-k` | `""` | Comma-separated classification keywords, such as `circadian,sleep,rhythm`. |
| `--reference-limit`, `-r` | `500` | Maximum candidate references retained per gene before metadata filtering. |
| `--input-gene-col` | `ext_gene` | Source gene-symbol column used during FBgn conversion. |
| `--input-species` | `fly` | Species represented by the input CSVs. Match `--species` for native human or mouse input. |
| `--flybase-data-dir` | auto-detected | Optional FlyBase data directory. |
| `--force-all` | off | Ignore checkpoint/run-store reuse and estimate or process everything again. |
| `--soft-run` | off | Estimate OpenAI API usage and cost from metadata only. Does not call OpenAI or fetch full paper text. |
| `--orthologs` | `none` | Optionally map fly genes to `human` or `mouse` orthologs before classification. LLM summaries and classifications are organism-strict for the selected species. |
| `--species` | `fly` | Target species for classification. With matching `--input-species`, runs natively without ortholog mapping. |
| `--diopt-filter` | `exclude_low_score_2` | DIOPT filter used for ortholog lookup. |
| `--human-data-dir` | `$FLAI_CACHE_DIR/human-data` or `~/.cache/flai-gene-classification/human-data` | Cache directory for human HGNC/gene2pubmed/GeneRIF/UniProt data. |
| `--mouse-data-dir` | `$FLAI_CACHE_DIR/mouse-data` or `~/.cache/flai-gene-classification/mouse-data` | Cache directory for mouse MGI/gene2pubmed/GeneRIF/UniProt data. |
| `--refresh-human-data` | off | Re-download cached human reference data before processing. |
| `--refresh-mouse-data` | off | Re-download cached mouse reference data before processing. |
| `--diopt-cache-dir` | auto/cache default | Optional cache directory for DIOPT API responses. |
| `--diopt-workers` | `8` | Concurrent DIOPT lookup workers for ortholog mapping. |

## Basic Fly Classification

Classify fly genes using the default input gene column and reference limit:

```bash
python flai_gene_classification.py ./InputCSVs
```

## Native Human Classification

Prepare a directory containing human symbols directly. Preparation resolves symbols through the HGNC catalog and writes `human_entrez_gene_id`, `human_gene_symbol`, and `hgnc_id` into the source CSV:

```bash
python flai_gene_classification.py prepare ./HumanInput \
  --species human \
  --input-species human \
  --input-gene-col human_gene_symbol
```

Then estimate or classify without FlyBase or DIOPT:

```bash
python flai_gene_classification.py estimate ./HumanInput \
  --species human \
  --input-species human \
  --keywords "sleep,circadian,Alzheimer" \
  --reference-limit 20

python flai_gene_classification.py classify ./HumanInput \
  --species human \
  --input-species human \
  --keywords "sleep,circadian,Alzheimer" \
  --reference-limit 300
```

## Keyword-Scoped Classification

Limit classification categories and literature search/filtering to a keyword set:

```bash
python flai_gene_classification.py ./InputCSVs \
  --keywords "circadian,sleep,rhythm"
```

## Lower Reference Limit

Reduce the number of retained candidate references per gene before metadata filtering:

```bash
python flai_gene_classification.py ./InputCSVs \
  --keywords "metabolism,feeding" \
  --reference-limit 100
```

## Soft-Run Cost Estimate

Estimate OpenAI usage and cost without calling OpenAI or fetching paper bodies:

```bash
python flai_gene_classification.py ./InputCSVs \
  --keywords "circadian,sleep,rhythm" \
  --reference-limit 100 \
  --soft-run
```

The soft-run uses a fast candidate-count proxy with an empirical species keyword-pass rate. It does not call OpenAI, fetch full paper text, or resolve every candidate title/abstract.

Species-aware prompt and reference-filter changes are cache-versioned. The
first real run after a cache schema bump regenerates older cached summaries;
use `--force-all` when you want to explicitly discard all checkpoint reuse.

## Force a Fresh Estimate or Run

Ignore existing checkpoint/run-store reuse:

```bash
python flai_gene_classification.py ./InputCSVs \
  --keywords "circadian,sleep,rhythm" \
  --reference-limit 100 \
  --soft-run \
  --force-all
```

Use the same flag without `--soft-run` to reprocess a real classification run from scratch.

## Human Ortholog Mode

Map fly genes to human orthologs before estimating or classifying:

```bash
python flai_gene_classification.py ./InputCSVs \
  --orthologs human \
  --keywords "neurodegeneration,sleep" \
  --soft-run
```

Fly-to-human runs summarize and classify only human-gene evidence. Europe PMC,
PubMed, Gene2PubMed, GeneRIF, and UniProt discovery are scoped to `Homo
sapiens`, and LLM summaries skip references that only describe non-human
orthologs.

Native human runs use the same organism-strict evidence filters but skip the fly-to-human mapping step entirely.

## Mouse Ortholog Mode

Map fly genes to mouse orthologs:

```bash
python flai_gene_classification.py ./InputCSVs \
  --orthologs mouse \
  --diopt-filter "exclude_low_score_2" \
  --keywords "development,behavior"
```

Mouse runs summarize and classify only mouse-gene evidence. Europe PMC,
PubMed, Gene2PubMed, GeneRIF, and UniProt discovery are scoped to `Mus
musculus`, and LLM summaries skip references that only describe non-mouse
orthologs.

## Fly Symbol to FBgn Mapping

The preparation step (and `python GetFBgnIDs.py <input_dir> [gene_column]` on its own) adds `flybase_gene_id` and `mapping_note` columns. The gene column is never changed. Symbols are matched exactly (case-sensitive): current symbol, then synonym, then, for a synonym shared by several genes, the gene FlyBase publications cite most when the lead is clear. No OpenAI or web-API call is made.

| `mapping_note` | `flybase_gene_id` | Meaning |
| --- | --- | --- |
| (empty) | FBgn | Exact current symbol. |
| `matched via synonym of <symbol>` | FBgn | Exact synonym of exactly one gene. |
| `synonym shared by N genes; FlyBase citations: ... → <symbol>` | FBgn | Shared synonym with a clear citation leader. |
| `ambiguous synonym: <symbol> <FBgn> (<pubs>), ...` | `-` | Shared synonym, no clear leader; candidates listed. |
| `ambiguous synonym (evidence index missing; ...)` | `-` | Citation index not built; run `python HelperScripts/build_synonym_evidence.py`. |
| `Greek letters spelled out for lookup` | FBgn | Matched only after spelling out Greek letters in a lookup key. |
| (empty) | `-` | No match. |

The citation index is built once from the FlyBase Chado XML (about 6 GB, `Data/FlyBase_Chado/`, gitignored).

## Custom Data Directories

Use explicit cache/data locations for reproducible runs:

```bash
python flai_gene_classification.py ./InputCSVs \
  --flybase-data-dir ./Data/FlyBase \
  --human-data-dir ./Data/Human \
  --diopt-cache-dir ./Data/DIOPT \
  --orthologs human
```

## Faster Ortholog Lookup

Increase DIOPT lookup concurrency for larger ortholog runs:

```bash
python flai_gene_classification.py ./InputCSVs \
  --orthologs human \
  --diopt-workers 16 \
  --diopt-cache-dir ./Data/DIOPT \
  --soft-run
```

Use higher worker counts cautiously; very high values can hit remote API limits.
