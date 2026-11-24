# OpenAI Model History

Archive of models this pipeline used as defaults and has since retired. Current
docs describe only the current defaults (`gpt-6-luna` for summaries,
`gpt-6.1-sol` for classification); anything that needs the old numbers lives here.

## Retired defaults

| Model | Role | Retired | Why | Replaced by |
| --- | --- | --- | --- | --- |
| `gpt-5.4-nano` | Per-chunk summaries and multi-chunk rewrite | 2026-10-07 (this pipeline) | OpenAI deprecated it on 2026-10-01; API removal 2027-04-01 | `gpt-6-luna` |
| `gpt-5.4` | Per-gene final classification | 2026-10-07 (this pipeline) | Not deprecated; replaced for lower price | `gpt-6.1-sol` |

Source for the deprecation: <https://developers.openai.com/api/docs/deprecations>.

## Prices (USD per 1M tokens)

From `Data/openai_pricing.csv` (LiteLLM, fetched 2026-10-07).

| Stage | Retired model | Input / Cached / Output | Current model | Input / Cached / Output | Change |
| --- | --- | --- | --- | --- | --- |
| Summary | `gpt-5.4-nano` | $0.20 / $0.02 / $1.25 | `gpt-6-luna` | $0.10 / $0.01 / $0.50 | -50% / -50% / -60% |
| Classification | `gpt-5.4` | $2.50 / $0.25 / $15.00 | `gpt-6.1-sol` | $2.00 / $0.10 / $10.00 | -20% / -60% / -33% |

## Expected cost at the switch

Uses the pooled token averages below, no cache discount.

| | Retired defaults | Current defaults | Change |
| --- | --- | --- | --- |
| Summary, per 1,000 references | $4.71 | $2.30 | -51% |
| Classification, per 1,000 genes | $5.98 | $4.20 | -30% |
| Example: CSW-Large-sized run (15,399 refs, 557 genes) | ~$75.93 | ~$37.80 | -50% |

Token averages were measured on the retired models. Output includes reasoning
tokens, which the current models may produce in different amounts. Re-run
`HelperScripts/calibrate_soft_run_profile.py` after the first runs on the
current models to get a measured profile.

## Calibration profiles measured on the retired models

### Original profile (Tx-Omics_Conserved-Small, `--reference-limit 50`)

Based on the saved fly-to-human and fly-to-mouse runs, `gpt-5.4-nano` for
summaries and `gpt-5.4` for classification. Output split was inferred from total
cost and model prices.

- 147 references queried; 26 genes reached GPT classification; 389 model requests.
- `gpt-5.4-nano` used 2,366,000 input tokens; `gpt-5.4` used 20,558.8 input tokens.
- Combined output: 98,010 tokens; cost about $0.75.

| Metric | Value |
| --- | --- |
| Summary requests per reference | 2.47 |
| Summary input tokens per reference | 16,095 |
| Summary output tokens per reference | 616 |
| Classification input tokens per gene | 791 |
| Classification output tokens per gene | 288 |

### Pooled per-request telemetry (replaced the original profile)

Pooled from 10 run stores under `Data/**/.batch_state/` on the retired models
(stores with fewer than 50 references excluded):
100,010 references and 3,245 classified genes.

| Metric | Original profile | Pooled |
| --- | --- | --- |
| Summary requests per reference | 2.47 | 3.05 |
| Summary input tokens per reference | 16,095 | 20,848 |
| Summary output tokens per reference | 616 | 436 |
| Classification input tokens per gene | 791 | 639 |
| Classification output tokens per gene | 288 | 292 |

Per run, summary input per reference ranged from about 16.5K to 22.5K.
Cached-input share was about 1.8% for summaries and 0% for classification.

## Worked-example comparison

The methodology worked example (100 refs, 50 genes, fixed token caps) totalled
**$3.53** at retired-model prices (summary $0.20 / $1.25, classification
$2.50 / $15.00). At current prices it totals **$2.33**.
