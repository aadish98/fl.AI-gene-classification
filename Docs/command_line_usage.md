# Command Line Usage

## Cheat Sheet

Run the pipeline from the repository root:

```bash
python flai_gene_classification.py <input_directory> [options]
```

| Argument | Default | Purpose |
| --- | --- | --- |
| `input_directory` | required | Directory containing input CSV files. |
| `--keywords`, `-k` | `""` | Comma-separated classification keywords, such as `circadian,sleep,rhythm`. |
| `--reference-limit`, `-r` | `500` | Maximum candidate references retained per gene before metadata filtering. |
| `--input-gene-col` | `ext_gene` | Source gene-symbol column used during FBgn conversion. |
| `--flybase-data-dir` | auto-detected | Optional FlyBase data directory. |
| `--force-all` | off | Ignore checkpoint/run-store reuse and estimate or process everything again. |
| `--soft-run` | off | Estimate OpenAI API usage and cost from metadata only. Does not call OpenAI or fetch full paper text. |
| `--orthologs` | `none` | Optionally map fly genes to `human` or `mouse` orthologs before classification. |
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

## Mouse Ortholog Mode

Map fly genes to mouse orthologs:

```bash
python flai_gene_classification.py ./InputCSVs \
  --orthologs mouse \
  --diopt-filter "exclude_low_score_2" \
  --keywords "development,behavior"
```

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
