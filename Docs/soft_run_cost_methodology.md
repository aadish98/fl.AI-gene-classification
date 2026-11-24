# Soft-Run Cost Methodology

This document explains how the soft-run cost estimate is calculated. It is a budgeting estimate for a planned run, not a guarantee of final spend.

## What Soft-Run Does

Soft-run estimates OpenAI API usage without sending any OpenAI requests. It also does not fetch full paper text. By default it uses candidate-reference counts and an empirical keyword-pass rate, so it does not resolve every candidate paper title/abstract.

The estimate is intentionally conservative because it cannot observe future model behavior, generated reasoning tokens, full-text chunk counts, or whether a model would reject a reference as not useful.

## Reference Count

For each gene, soft-run estimates references that would pass metadata-level filtering. Let:

- $c^{(g)}$ be the number of candidate references for gene $g$.
- $\rho_s$ be the empirical keyword-pass rate for species $s$.
- $m^{(g)}$ be the estimated number of keyword-passing references for gene $g$.
- $n^{(g)}$ be the capped number of references estimated to be OpenAI-bound for gene $g$.

Soft-run first applies the same per-gene reference limit as the real pipeline:

$$
\ell^{(g)} = \min(c^{(g)},\ L)
$$

Then it estimates keyword-passing references:

$$
m^{(g)} = \ell^{(g)} \rho_s
$$

The final per-gene OpenAI-bound estimate is:

$$
n^{(g)} = \min(m^{(g)},\ L)
$$

The total OpenAI-bound reference count is:

$$
N_{\text{refs}} = \sum_{g} n^{(g)}
$$

This cap prevents a single gene with many metadata matches from dominating the estimate.

Reference collection uses the same organism-scoped source adapters as a real
run. In particular, Europe PMC is filtered with structured `ORGANISM` / `MESH`
clauses for fly, human, and mouse before the candidate count is estimated, so
the soft-run does not budget for off-organism EPMC-only references that the
real pipeline would no longer summarize.

## Keyword-Pass Calibration

Species keyword-pass rates live in `.flai_system/soft_run_reference_profiles.json` under `species_keyword_pass_rates`. If no species rate exists, soft-run falls back to `1.0`, which makes the reference estimate a candidate-count upper bound.

The initial human and mouse seed values come from `Data/Tx-Omics_Conserved-Small`:

- Human: `84 / (16 * 50) = 0.105`.
- Mouse: `63 / (16 * 50) = 0.07875`.

These seeds are conservative because the denominator assumes every mapped ortholog gene had the full `--reference-limit 50` candidate references.

Completed real runs update a species pass rate only when the new run has at least as many limited references as the stored calibration:

$$
\rho_s = \frac{\text{metadata keyword matches}}{\text{limited references}}
$$

## Gene Classification Count

A gene is counted as classification-bound if it has at least one estimated OpenAI-bound reference:

$$
N_{\text{genes}} = \left|\{g : n^{(g)} \ge 1\}\right|
$$

## Token Constants

Soft-run uses four constants:

- $\tau^{\text{in}}_{\text{ref}}$: estimated input tokens per OpenAI-bound reference.
- $\tau^{\text{out}}_{\text{ref}}$: capped output tokens per OpenAI-bound reference.
- $\tau^{\text{in}}_{\text{cls}}$: estimated input tokens per gene classification.
- $\tau^{\text{out}}_{\text{cls}}$: capped output tokens per gene classification.

Output uses the configured generation cap rather than an average because reasoning tokens and visible output tokens are only known after a real response. The capped output estimate is therefore an upper-bound component, not expected spend.

The live summary/classification prompts include target-organism context and a
`species_specific` summary field. This adds a small fixed prompt overhead, but
`SHARED_PROMPT_OVERHEAD_WORDS=800` remains the runtime chunk-budget reserve and
the fixed soft-run constants above remain conservative for budgeting.

## Empirical Reference Averages

When a calibration profile is available, soft-run also prints a reference-average planning estimate alongside the cap-based upper bound.

Calibration profiles live in `.flai_system/soft_run_reference_profiles.json`. They are project-level system calibration inputs, not run outputs, so they are kept outside `Data/`.

Profile selection:

1. A profile keyed on the exact configured (summary model, classification model) pair, if one exists.
2. Otherwise the model-independent default profile (`"default": true`). Token counts depend on the prompts and papers rather than the model, so the default profile is repriced at the current models' rates, and the report says so.
3. If neither exists, only the upper-bound estimate is printed.

The default profile is pooled from per-request telemetry in 10 saved run stores (100,010 references, 3,245 classified genes):

- Summary requests per queried reference: 3.05.
- Summary input tokens per queried reference: 20,848.
- Summary output tokens per queried reference: 436.
- Classification input tokens per GPT-classified gene: 639.
- Classification output tokens per GPT-classified gene: 292.

Rebuild it with `HelperScripts/calibrate_soft_run_profile.py` (`--dry-run` to preview). When run stores exist for the currently configured model pair, the script writes a profile keyed on that pair instead. Output tokens include reasoning tokens, which differ by model, so a pair-specific profile is more accurate once available. The earlier profile and the models it was measured on are in [`model_history.md`](model_history.md).

New real runs record successful per-request OpenAI usage in `.batch_state/run_store.json` under each cached gene record:

- `openai_usage.requests`: one entry per successful OpenAI request, with request type, model, reasoning effort, reference/gene context, and the returned usage payload.
- `openai_usage.totals`: per-gene aggregate request count and token totals.
- `summaries[].chunk_count` and `summaries[].chunk_word_counts`: per-reference full-text chunking metrics used to explain why actual reference input can exceed the fixed soft-run constant.

These fields are what `HelperScripts/calibrate_soft_run_profile.py` reads to build calibration profiles.

## Token Totals

The estimated input token total is:

$$
T_{\text{in}} =
N_{\text{refs}} \cdot \tau^{\text{in}}_{\text{ref}}
+ N_{\text{genes}} \cdot \tau^{\text{in}}_{\text{cls}}
$$

The capped output token total is:

$$
T_{\text{out}} =
N_{\text{refs}} \cdot \tau^{\text{out}}_{\text{ref}}
+ N_{\text{genes}} \cdot \tau^{\text{out}}_{\text{cls}}
$$

## Cost Formula

Let:

- $p^{\text{in}}_{s}$ be the summary model input price per 1M tokens.
- $p^{\text{out}}_{s}$ be the summary model output price per 1M tokens.
- $p^{\text{in}}_{c}$ be the classification model input price per 1M tokens.
- $p^{\text{out}}_{c}$ be the classification model output price per 1M tokens.

The estimated total upper-bound cost is:

$$
C =
\frac{N_{\text{refs}} \tau^{\text{in}}_{\text{ref}}}{10^{6}} p^{\text{in}}_{s}
+ \frac{N_{\text{refs}} \tau^{\text{out}}_{\text{ref}}}{10^{6}} p^{\text{out}}_{s}
+ \frac{N_{\text{genes}} \tau^{\text{in}}_{\text{cls}}}{10^{6}} p^{\text{in}}_{c}
+ \frac{N_{\text{genes}} \tau^{\text{out}}_{\text{cls}}}{10^{6}} p^{\text{out}}_{c}
$$

Soft-run reports this as three separate values:

- Estimated input cost
- Capped output cost
- Estimated total upper-bound cost

## Pricing Source

Pricing is refreshed from an external maintained pricing table before each soft-run. The refreshed table is stored locally as a CSV so model pricing can be reused if the network is unavailable later.

The refresh is a plain public HTTPS download. It does not call OpenAI and does not require an OpenAI API key.

## Known Limitations

Soft-run does not fetch full paper text, so it does not model how many full-text chunks a paper would produce.

Soft-run uses an empirical keyword-pass proxy instead of resolving every candidate paper's metadata, so the reference count is an estimate. It also does not know which references the model would later judge as high quality, so it does not model early stopping based on high-quality evidence.

Soft-run estimates usage from fixed token constants. It is designed for budgeting and run planning, not penny-level accounting.

## Worked Example

Assume:

- $N_{\text{refs}} = 100$
- $N_{\text{genes}} = 50$
- $\tau^{\text{in}}_{\text{ref}} = 6500$
- $\tau^{\text{out}}_{\text{ref}} = 1200$
- $\tau^{\text{in}}_{\text{cls}} = 2000$
- $\tau^{\text{out}}_{\text{cls}} = 4000$
- $p^{\text{in}}_{s} = 0.10$
- $p^{\text{out}}_{s} = 0.50$
- $p^{\text{in}}_{c} = 2.00$
- $p^{\text{out}}_{c} = 10.00$

Then:

$$
T_{\text{in}} = 100 \cdot 6500 + 50 \cdot 2000 = 750000
$$

$$
T_{\text{out}} = 100 \cdot 1200 + 50 \cdot 4000 = 320000
$$

The cost components are:

$$
C_{\text{summary,in}} = \frac{100 \cdot 6500}{10^6} \cdot 0.10 = 0.065
$$

$$
C_{\text{summary,out}} = \frac{100 \cdot 1200}{10^6} \cdot 0.50 = 0.06
$$

$$
C_{\text{classification,in}} = \frac{50 \cdot 2000}{10^6} \cdot 2.00 = 0.20
$$

$$
C_{\text{classification,out}} = \frac{50 \cdot 4000}{10^6} \cdot 10.00 = 2.00
$$

So:

$$
C = 0.065 + 0.06 + 0.20 + 2.00 = 2.325
$$

In this example, the estimated total upper-bound cost is about $2.33. The same example at the previous defaults' prices is compared in [`model_history.md`](model_history.md).
