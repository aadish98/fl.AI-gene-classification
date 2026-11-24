# OpenAI API Usage in the fl.AI Gene Classification Pipeline

This document is an exhaustive reference of every place the pipeline calls the
OpenAI API, the exact prompts (system + user) sent, the structured-output schema
expected back, the default model, and how each call can be overridden via
environment variables.

For every call site we show **both** the f-string template with `{placeholder}`
fields **and** a fully-rendered "live example" that mirrors what is actually
sent on the wire — for fly, human, and mouse where the species adapter changes
the input, and for the species-agnostic rewrite call once. The live examples
use real symbols / IDs / synonyms from the gene catalogs the pipeline loads
(FlyBase, HGNC, MGI / NCBI Gene); the title/abstract/chunk bodies are abridged
illustrations so they remain readable.

Everything in this document lives in `flai_gene_classification.py`. No other
module in the repo (helper scripts for FlyBase, DIOPT orthologs, HGNC/MGI
mappers, PubMed/PMC fetchers, Unpaywall, `summarize_gene_set.py`, the soft-run
estimator, etc.) makes any OpenAI request. In particular:

- **Ortholog mapping** is done by `HelperScripts/get_orthologs.py` /
  `HelperScripts/diopt_client.py` against the DIOPT REST API, **not** GPT.
- **Gene-set summary** (`summarize_gene_set.py`) is pure pandas / openpyxl.
- **Soft-run cost estimation** is metadata-only and never calls OpenAI.

---

## 1. Global configuration

### 1.1 Default models

Defined at the top of `flai_gene_classification.py`:

```71:74:flai_gene_classification.py
DEFAULT_SUMMARY_MODEL = "gpt-6-luna"
DEFAULT_CLASSIFICATION_MODEL = "gpt-6.2-sol"
DEFAULT_SUMMARY_REASONING_EFFORT = "medium"
DEFAULT_CLASSIFICATION_REASONING_EFFORT = "high"
```

| Purpose | Default model | Default reasoning effort |
| --- | --- | --- |
| Per-chunk reference summarization | `gpt-6-luna` | `medium` |
| Multi-chunk paper-level rewrite | `gpt-6-luna` | `medium` |
| Per-gene final classification | `gpt-6.2-sol` | `high` |

Reasoning effort is only attached to the request if the resolved model name
is in the GPT-5-or-later family (matches `gpt-<N>` with `N >= 5`,
case-insensitive, e.g. `gpt-6-luna`, `gpt-6.2-sol`) and the effort string is
non-empty. If
the OpenAI server rejects the parameter (`Unsupported parameter:
'reasoning.effort'`), the call is automatically retried once without the
`reasoning` field. See `_get_reasoning_config` and the retry block in
`_parse_structured_completion`.

### Model pricing

Prices are read from `Data/openai_pricing.csv` (refreshed from LiteLLM on every
`--soft-run`).

| Model (USD / 1M tokens) | Input | Cached input | Output |
| --- | --- | --- | --- |
| `gpt-6-luna` (summary default) | $0.10 | $0.01 | $0.50 |
| `gpt-6.2-sol` (classification default) | $2.00 | $0.10 | $10.00 |

Expected cost from the pooled run telemetry (~20,848 input / ~436 output tokens
per reference; ~639 input / ~292 output tokens per classified gene; no cache
discount): about **$2.30 per 1,000 references** and **$4.20 per 1,000
classified genes**. Output tokens include reasoning tokens, so actual cost
depends on how many the models emit at the configured effort.

`gpt-6.2-sol` accepts `low`, `medium`, `high`, `xhigh` and `max` reasoning
effort, but not `none` or `minimal`.

Retired defaults, their prices, and the old cost comparison are in
[`model_history.md`](model_history.md).

### 1.2 Environment-variable overrides

Resolution order is **explicit argument → task-specific env var → default**.

| Env var | Effect | Used by |
| --- | --- | --- |
| `OPENAI_API_KEY` | Required to authenticate `OpenAI(api_key=...)` (lazy `_get_openai_client`). | All calls |
| `OPENAI_SUMMARY_MODEL` | Override summary / rewrite model. | `_get_summary_model` |
| `OPENAI_CLASSIFICATION_MODEL` | Override classification model. | `_get_classification_model` |
| `OPENAI_SUMMARY_REASONING_EFFORT` | Override reasoning effort for summary + rewrite calls. | `_get_reasoning_effort` |
| `OPENAI_CLASSIFICATION_REASONING_EFFORT` | Override reasoning effort for classification calls. | `_get_reasoning_effort` |

There are **no CLI flags** for model selection; everything goes through env
vars (or by editing the constants).

### 1.3 Token / chunk budgets

```81:84:flai_gene_classification.py
INPUT_TOKENS_PER_REF = 6500
OUTPUT_CAP_PER_REF = 1200
INPUT_TOKENS_PER_CLASSIFICATION = 2000
OUTPUT_CAP_PER_CLASSIFICATION = 4000
```

```215:217:flai_gene_classification.py
MAX_PROMPT_WORDS = 5000
SHARED_PROMPT_OVERHEAD_WORDS = 800
MIN_CHUNK_WORDS = 750
```

`OUTPUT_CAP_PER_REF` is sent as `max_output_tokens` for both summary and
rewrite calls; `OUTPUT_CAP_PER_CLASSIFICATION` is used for the primary
classification call (the short-context retry uses `max_output_tokens=2000`).
`INPUT_TOKENS_PER_REF` / `INPUT_TOKENS_PER_CLASSIFICATION` are only used by
the soft-run cost estimator, not by live requests.

Reference text is split with `_build_shared_prompt_chunks`: the shared
prompt header (system + metadata) is reserved as `SHARED_PROMPT_OVERHEAD_WORDS
+ words(metadata_blob)`, leaving up to `MAX_PROMPT_WORDS - overhead` words per
chunk (with a floor of `MIN_CHUNK_WORDS`).

### 1.4 The single API entry point

Every live call goes through `_parse_structured_completion`, which uses the
**Responses API with structured parsing** (`client.responses.parse`):

```2173:2199:flai_gene_classification.py
            request_kwargs: dict[str, Any] = {
                "model": model,
                "input": input_items,
                "text_format": response_format,
                "max_output_tokens": max_output_tokens,
            }
            if instructions:
                request_kwargs["instructions"] = instructions
            reasoning = _get_reasoning_config(model, reasoning_effort=reasoning_effort)
            if reasoning is not None:
                request_kwargs["reasoning"] = reasoning

            response = _get_openai_client().responses.parse(**request_kwargs)
            parsed = getattr(response, "output_parsed", None)
```

Behaviour:

- Chat-style messages are split by `_build_responses_input`: every `system`
  message is concatenated into the `instructions` field; everything else
  becomes `input` items with their original `role`.
- Up to **3 attempts** with `1.5 * (attempt + 1)` second backoff between
  attempts.
- A 401 / `Unauthorized` / `Incorrect API key` failure raises immediately as
  `OpenAI authentication failed: ...`.
- Token usage from each successful attempt is recorded via
  `_record_openai_usage_event` for downstream cost reporting.

---

## 2. Pipeline feature 1 — Per-chunk reference summarization

**Function:** `summarize_reference_function`
**Caller:** `process_gene_set` loop, once per `(gene, paper, chunk)` triple.
**Default model:** `gpt-6-luna` (override `OPENAI_SUMMARY_MODEL`).
**Reasoning effort:** `medium` (override `OPENAI_SUMMARY_REASONING_EFFORT`).
**`max_output_tokens`:** `OUTPUT_CAP_PER_REF` = `1200`.
**Structured response model:** `FunctionPhenotypeSummary`

```1909:1914:flai_gene_classification.py
class FunctionPhenotypeSummary(BaseModel):
    function: str = ""
    phenotypes: str = ""
    skip_reference: bool = False
    skip_reason: str = ""
    species_specific: bool = True
```

### 2.1 System prompt (verbatim)

```text
You are an expert biomedical assistant focused on the {scientific_name} gene `{gene_symbol}` (NCBI taxon ID {taxon_id}). Use ONLY the provided title, abstract, and text chunk. Do not describe orthologs, paralogs, or homologs in any other species. If the chunk only discusses a non-{display_name} ortholog/paralog, or contains no direct evidence about the {display_name} gene `{gene_symbol}`, set skip_reference=true and species_specific=false.
```

### 2.2 User prompt template

`{display_name}`, `{scientific_name}`, `{taxon_id}`, `{gene_symbol}`,
`{gene_id}`, `{synonyms_str}` (sorted unique, comma-separated), `{title}`,
`{abstract}`, and `{text_chunk}` are substituted. The exact body sent (see
`_build_summary_messages`):

```text
Target species: {scientific_name} (NCBI taxon ID {taxon_id})
Target gene: {gene_symbol} (Gene ID: {gene_id}; Synonyms: {synonyms_str})

Return JSON with fields:
- function: concise summary of {display_name} gene {gene_symbol} function evidence from this chunk only
- phenotypes: concise summary of phenotype evidence of {display_name} gene {gene_symbol} from this chunk only
- skip_reference: true if this chunk does not contain direct {display_name}-specific evidence for {gene_symbol}
- skip_reason: brief reason when skip_reference is true
- species_specific: true only if the chunk provides evidence explicitly about the {display_name} gene {gene_symbol}; false for evidence only about another organism's ortholog/paralog/homolog

Rules:
- Always name {gene_symbol} explicitly
- Be factual and use only the provided text
- Describe only {display_name} biology; ignore orthologs, paralogs, or homologs in other species
- If evidence is generic, cross-species without {display_name}-specific evidence, or only about a non-{display_name} ortholog, set skip_reference=true and species_specific=false
- Set species_specific=true only when the chunk's evidence is explicitly about the {display_name} gene `{gene_symbol}`
- If one field is absent but the chunk still contains evidence, return an empty string for that field

Title: {title or ''}
Abstract: {abstract or ''}

Text chunk:
{text_chunk or ''}
```

If `species_specific=false` and the model did **not** already mark
`skip_reference=true`, the post-processor forces `skip_reference=true` with a
species-specific skip reason. If both `function` and `phenotypes` come back
empty and the model did **not** mark `skip_reference=true`, the post-processor
in `summarize_reference_function` forces `skip_reference=true` with reason
`"No function or phenotype evidence returned"`. On any non-auth exception the
function returns `{function: "", phenotypes: "", skip_reference: true,
skip_reason: "API error: <e>"}` so the pipeline keeps going.

### 2.3 What text is sent

`text_chunk` comes from `_build_shared_prompt_chunks` and corresponds to one
piece of the paper's full text (PMC OA / Europe PMC / Unpaywall / similar)
sized to fit `MAX_PROMPT_WORDS - overhead`. If the gene only has a GeneRIF /
title+abstract evidence (no full text), the chunk loop is still executed once
but with an empty body or with the GeneRIF-augmented text appended (see
`process_gene_set` around the `GENERIF:` concatenation).

### 2.4 Live examples — fully-rendered per-chunk summary prompts

The system prompt is species-aware. It changes based on `display_name`,
`scientific_name`, `taxon_id`, `gene_symbol`, `gene_id` (the species'
`gene.gene_id`: FlyBase ID for fly, NCBI Entrez Gene ID for human and mouse),
and the synonym list (sorted, comma-joined by `_build_summary_messages`).

#### 2.4.1 Fly — `per` (period)

`gene_symbol="per"`, `gene_id="FBgn0003068"`,
`synonyms=["CG2647", "dPER", "PER", "Per", "dPer", "period", "per"]`.

System message:

```text
You are an expert biomedical assistant focused on the Drosophila melanogaster gene `per` (NCBI taxon ID 7227). Use ONLY the provided title, abstract, and text chunk. Do not describe orthologs, paralogs, or homologs in any other species. If the chunk only discusses a non-Fly ortholog/paralog, or contains no direct evidence about the Fly gene `per`, set skip_reference=true and species_specific=false.
```

User message:

```text
Target species: Drosophila melanogaster (NCBI taxon ID 7227)
Target gene: per (Gene ID: FBgn0003068; Synonyms: CG2647, PER, Per, dPER, dPer, per, period)

Return JSON with fields:
- function: concise summary of Fly gene per function evidence from this chunk only
- phenotypes: concise summary of phenotype evidence of Fly gene per from this chunk only
- skip_reference: true if this chunk does not contain direct Fly-specific evidence for per
- skip_reason: brief reason when skip_reference is true
- species_specific: true only if the chunk provides evidence explicitly about the Fly gene per; false for evidence only about another organism's ortholog/paralog/homolog

Rules:
- Always name per explicitly
- Be factual and use only the provided text
- Describe only Fly biology; ignore orthologs, paralogs, or homologs in other species
- If evidence is generic, cross-species without Fly-specific evidence, or only about a non-Fly ortholog, set skip_reference=true and species_specific=false
- Set species_specific=true only when the chunk's evidence is explicitly about the Fly gene `per`
- If one field is absent but the chunk still contains evidence, return an empty string for that field

Title: Molecular analysis of the period locus in Drosophila melanogaster
Abstract: [abstract abridged for example] We characterise per mutant alleles per^S, per^L, and per^0 and show that per encodes a circadian clock component required for normal locomotor rhythmicity in adult flies.

Text chunk:
[text chunk abridged for example] In wild-type Canton-S flies, per mRNA cycles with a ~24 h period in the heads of adults entrained to LD 12:12. Loss-of-function per^0 animals are arrhythmic in DD; per^S shortens the free-running period to ~19 h; per^L lengthens it to ~28 h. Rescue with a wild-type per transgene restores rhythmicity. ...
```

#### 2.4.2 Human — `PER1`

`gene_symbol="PER1"`, `gene_id="5187"` (Entrez Gene ID),
`synonyms=["KIAA0482", "PER", "PER1", "RIGUI", "hPER", "hPER1"]`.

```text
Target species: Homo sapiens (NCBI taxon ID 9606)
Target gene: PER1 (Gene ID: 5187; Synonyms: KIAA0482, PER, PER1, RIGUI, hPER, hPER1)

Return JSON with fields:
- function: concise summary of Human gene PER1 function evidence from this chunk only
- phenotypes: concise summary of phenotype evidence of Human gene PER1 from this chunk only
- skip_reference: true if this chunk does not contain direct Human-specific evidence for PER1
- skip_reason: brief reason when skip_reference is true
- species_specific: true only if the chunk provides evidence explicitly about the Human gene PER1; false for evidence only about another organism's ortholog/paralog/homolog

Rules:
- Always name PER1 explicitly
- Be factual and use only the provided text
- Describe only Human biology; ignore orthologs, paralogs, or homologs in other species
- If evidence is generic, cross-species without Human-specific evidence, or only about a non-Human ortholog, set skip_reference=true and species_specific=false
- Set species_specific=true only when the chunk's evidence is explicitly about the Human gene `PER1`
- If one field is absent but the chunk still contains evidence, return an empty string for that field

Title: A PER1 polymorphism associated with shorter sleep duration in humans
Abstract: [abstract abridged for example] We genotyped PER1 SNPs in N=512 healthy adults phenotyped for habitual sleep duration by actigraphy and self-report.

Text chunk:
[text chunk abridged for example] PER1 rs2735611 minor-allele carriers showed ~22 min shorter total sleep time and reduced slow-wave activity in NREM relative to non-carriers. PER1 mRNA in peripheral blood was rhythmic with peak expression at ZT2-4. ...
```

#### 2.4.3 Mouse — `Per1`

`gene_symbol="Per1"`, `gene_id="18626"` (Entrez Gene ID),
`synonyms=["Per", "Per1", "Rigui", "m-Rigui", "mPer1"]`.

```text
Target species: Mus musculus (NCBI taxon ID 10090)
Target gene: Per1 (Gene ID: 18626; Synonyms: Per, Per1, Rigui, m-Rigui, mPer1)

Return JSON with fields:
- function: concise summary of Mouse gene Per1 function evidence from this chunk only
- phenotypes: concise summary of phenotype evidence of Mouse gene Per1 from this chunk only
- skip_reference: true if this chunk does not contain direct Mouse-specific evidence for Per1
- skip_reason: brief reason when skip_reference is true
- species_specific: true only if the chunk provides evidence explicitly about the Mouse gene Per1; false for evidence only about another organism's ortholog/paralog/homolog

Rules:
- Always name Per1 explicitly
- Be factual and use only the provided text
- Describe only Mouse biology; ignore orthologs, paralogs, or homologs in other species
- If evidence is generic, cross-species without Mouse-specific evidence, or only about a non-Mouse ortholog, set skip_reference=true and species_specific=false
- Set species_specific=true only when the chunk's evidence is explicitly about the Mouse gene `Per1`
- If one field is absent but the chunk still contains evidence, return an empty string for that field

Title: Per1 deletion in mice perturbs circadian sleep architecture
Abstract: [abstract abridged for example] Per1^-/- mice were phenotyped by EEG/EMG over 3 baseline days plus a 6 h sleep deprivation challenge.

Text chunk:
[text chunk abridged for example] Compared to wild-type littermates, Per1^-/- mice showed reduced NREM sleep consolidation, shorter free-running period (~23.4 h) in DD, and attenuated rebound delta power after sleep deprivation. Per1 transcript was undetectable in Per1^-/- SCN by qPCR. ...
```

A representative successful response from the model (independent of species)
matches the `FunctionPhenotypeSummary` schema, e.g.:

```json
{
  "function": "Per1 acts as a core component of the circadian transcriptional feedback loop required for normal NREM consolidation.",
  "phenotypes": "Per1^-/- mice show reduced NREM consolidation, shortened free-running period (~23.4 h), and attenuated delta-power rebound after 6 h sleep deprivation.",
  "skip_reference": false,
  "skip_reason": "",
  "species_specific": true
}
```

If the chunk had no Mouse-specific `Per1` evidence, or only described a
non-mouse ortholog, the model would return `skip_reference=true`,
`species_specific=false`, and a one-line `skip_reason`.

---

## 3. Pipeline feature 2 — Multi-chunk paper-level rewrite

**Function:** `rewrite_function_summary`
**Caller:** `process_gene_set`, only when **>=2** non-skipped chunk summaries
exist for the same paper. Single-chunk papers skip this call and reuse the
chunk summary directly.
**Default model:** `gpt-6-luna` (override `OPENAI_SUMMARY_MODEL`).
**Reasoning effort:** `medium` (override `OPENAI_SUMMARY_REASONING_EFFORT`).
**`max_output_tokens`:** `OUTPUT_CAP_PER_REF` = `1200`.
**Structured response model:** `FinalFunctionPhenotypeSummary`

```1916:1918:flai_gene_classification.py
class FinalFunctionPhenotypeSummary(BaseModel):
    function: str = ""
    phenotypes: str = ""
```

### 3.1 System prompt (verbatim)

```text
You are an expert biomedical assistant focused on the {scientific_name} gene `{gene_symbol}` (NCBI taxon ID {taxon_id}). Merge chunk-level evidence summaries into a single cohesive paper-level summary. Use ONLY the provided chunk summaries and do not introduce evidence about non-{display_name} orthologs, paralogs, or homologs.
```

### 3.2 User prompt template

`merged_input` is the per-chunk concatenation:

```text
Chunk 1
Function: {function or 'None'}
Phenotypes: {phenotypes or 'None'}

Chunk 2
Function: ...
Phenotypes: ...
```

The user message itself:

```text
Rewrite these chunk-level summaries for the {display_name} gene {gene_symbol} (Gene ID: {gene_id}) into one cohesive paper summary.

Return JSON with fields:
- function
- phenotypes

Rules:
- Always name {gene_symbol} explicitly
- Summarize only {display_name} biology for {gene_symbol}
- Preserve factual content only
- Remove redundancy across chunks
- Do not introduce new claims
- Do not include evidence about non-{display_name} orthologs/paralogs/homologs, even if input chunk summaries imply it

Chunk summaries:
{merged_input}
```

Failure handling: if the call raises (other than auth errors) the function
falls back to a deterministic concatenation of the per-chunk `function` /
`phenotypes` strings via `_join_unique_texts`. Empty model output for either
field also falls back to that deterministic concatenation.

### 3.3 Live example — paper-level rewrite for mouse `Per1`

Suppose the paper produced two chunk-level summaries because the full text
was split across two prompt chunks. The caller passes `gene.symbol="Per1"`,
`gene.gene_id="18626"`, and the two `chunk_summary` dicts that came back
from §2 (this prompt does **not** receive synonyms, title, abstract, or
the original chunk text — only the chunk-level structured outputs).

System message:

```text
You are an expert biomedical assistant focused on the Mus musculus gene `Per1` (NCBI taxon ID 10090). Merge chunk-level evidence summaries into a single cohesive paper-level summary. Use ONLY the provided chunk summaries and do not introduce evidence about non-Mouse orthologs, paralogs, or homologs.
```

User message:

```text
Rewrite these chunk-level summaries for the Mouse gene Per1 (Gene ID: 18626) into one cohesive paper summary.

Return JSON with fields:
- function
- phenotypes

Rules:
- Always name Per1 explicitly
- Summarize only Mouse biology for Per1
- Preserve factual content only
- Remove redundancy across chunks
- Do not introduce new claims
- Do not include evidence about non-Mouse orthologs/paralogs/homologs, even if input chunk summaries imply it

Chunk summaries:
Chunk 1
Function: Per1 is described as a core component of the mammalian circadian transcriptional feedback loop, expressed rhythmically in the SCN.
Phenotypes: Per1^-/- mice exhibit a shortened free-running period (~23.4 h) in constant darkness.

Chunk 2
Function: Per1 transcripts are induced by light pulses in the SCN and required for phase resetting of the circadian oscillator.
Phenotypes: Per1^-/- mice show reduced NREM sleep consolidation and attenuated delta-power rebound after 6 h of sleep deprivation.
```

A successful `FinalFunctionPhenotypeSummary` response merges both chunks
without introducing new claims, e.g.:

```json
{
  "function": "Per1 is a core component of the mammalian circadian transcriptional feedback loop, expressed rhythmically in the SCN and induced by light pulses, where it is required for phase resetting of the oscillator.",
  "phenotypes": "Per1^-/- mice exhibit a shortened free-running period (~23.4 h) in DD, reduced NREM sleep consolidation, and attenuated delta-power rebound after 6 h sleep deprivation."
}
```

If the model errors or returns empty fields, the pipeline falls back to a
deterministic merge of the two `Function:` / `Phenotypes:` strings via
`_join_unique_texts` and never blocks the gene.

---

## 4. Pipeline feature 3 — Per-gene final classification

**Function:** `classify_gene_from_text`
**Caller:** `process_gene_set`, once per gene that has at least one
high-quality reference summary. The aggregated input is
`"\n\n".join(high_quality_summary_blocks)`.
**Default model:** `gpt-6.2-sol` (override `OPENAI_CLASSIFICATION_MODEL`).
**Reasoning effort:** `high` (override `OPENAI_CLASSIFICATION_REASONING_EFFORT`).
**`max_output_tokens`:** `OUTPUT_CAP_PER_CLASSIFICATION` = `4000`
(short-context retry uses `2000`).
**Structured response model:** `GeneClassification`

```2470:2474:flai_gene_classification.py
class GeneClassification(BaseModel):
    gene: str
    category: Optional[List[str]]
    confidence: int
    rationale: str
```

### 4.1 System prompt (verbatim)

```text
You are an expert biomedical assistant focused on the {scientific_name} gene `{gene_symbol}` (NCBI taxon ID {taxon_id}). Use ONLY the provided text. Classify only evidence pertaining to the {display_name} gene. If evidence is weak, ambiguous, or predominantly about a non-target organism, return category as empty list [] or None.
```

### 4.2 User prompt template

`cats_line` is the comma-joined list of allowed keyword categories (the
literal string `None` if no keywords were supplied). The full text is
hard-truncated to **200 000 characters** before being embedded.

```text
Classify the {display_name} gene '{gene_symbol}' ({scientific_name}, NCBI taxon ID {taxon_id}). Allowed categories: [{cats_line}]. You may classify into 1 category, N categories, or None (return empty list [] or None). If the text predominantly describes a non-{display_name} ortholog/paralog/homolog, return an empty category list with rationale noting cross-species text. Return JSON with fields: gene, category (list of strings, can be empty or None), confidence (0-100), rationale (1-2 lines).

TEXT START
{full_text[:200000]}
TEXT END
```

### 4.3 Post-processing

- `category` is filtered against the user-supplied keyword set
  (`_normalize_category_output`); anything not in `valid_cats` is dropped, so
  the model cannot invent labels.
- `confidence` is clamped to `0..100` (`_safe_confidence`).
- `rationale` is coerced to string.

### 4.4 Length-error retry path

If the first attempt raises an error containing
`"length limit was reached"` or `"maximum context length"`, the function
**retries once** with the same prompt template but using
`full_text[:50000]`, `max_output_tokens=2000`, and the same
`OPENAI_CLASSIFICATION_REASONING_EFFORT`. The retry has its own usage event
tagged `request_type="classification_retry_short_context"`. Any other error
returns `{category: [], confidence: 0, rationale: "Error: <e>"}` so the
pipeline can continue.

### 4.5 No-evidence short-circuit (no API call)

If a gene ends up with **zero** high-quality summary chunks, the pipeline
skips the classification call entirely and writes:

```text
category: "None"
confidence: 0
rationale: "No relevant references found."
classified_by: "Self"
```

This is a pure fallback; no OpenAI request is made for that gene.

### 4.6 Live examples — fully-rendered classification prompts

The classification call receives `agg_text = "\n\n".join(chunks)`, where
each `chunk` is the per-paper rendering produced by
`_format_reference_summary` — i.e. `"- Function: ...\n- Phenotypes: ..."`
(the same text written into the `Reference Summary` column of the output
workbook). For the SLURM batch jobs in `BATCH_SCRIPTS/process_txomics_*`
the `keywords_list` resolved by the CLI is `["sleep", "circadian"]`, so
the rendered `cats_line` is `sleep, circadian`. The system prompt is
species-aware:

```text
You are an expert biomedical assistant focused on the {scientific_name} gene `{gene_symbol}` (NCBI taxon ID {taxon_id}). Use ONLY the provided text. Classify only evidence pertaining to the {display_name} gene. If evidence is weak, ambiguous, or predominantly about a non-target organism, return category as empty list [] or None.
```

#### 4.6.1 Fly — `per`

`gene_symbol="per"`, `keywords_list=["sleep", "circadian"]`, two
high-quality reference summaries aggregated.

User message:

```text
Classify the Fly gene 'per' (Drosophila melanogaster, NCBI taxon ID 7227). Allowed categories: [sleep, circadian]. You may classify into 1 category, N categories, or None (return empty list [] or None). If the text predominantly describes a non-Fly ortholog/paralog/homolog, return an empty category list with rationale noting cross-species text. Return JSON with fields: gene, category (list of strings, can be empty or None), confidence (0-100), rationale (1-2 lines).

TEXT START
- Function: per encodes a core circadian clock protein in Drosophila melanogaster, dimerising with TIM in a transcriptional/translational negative feedback loop driving rhythmic expression of clock genes.
- Phenotypes: per^0 null flies are arrhythmic in DD; per^S shortens the free-running locomotor period to ~19 h and per^L lengthens it to ~28 h. Sleep duration in adult flies is reduced in per^0 mutants relative to Canton-S controls.

- Function: per cycles in the lateral neurons (LNvs) of the adult Drosophila brain and is required for the molecular circadian oscillator in central pacemaker neurons.
- Phenotypes: Targeted RNAi knockdown of per in LNvs phenocopies per^0 arrhythmicity in locomotor activity assays.
TEXT END
```

Expected `GeneClassification` payload:

```json
{
  "gene": "per",
  "category": ["sleep", "circadian"],
  "confidence": 95,
  "rationale": "Multiple references describe per as a core circadian clock gene in Drosophila with arrhythmic locomotor phenotype in null mutants and reduced sleep duration in per^0 flies."
}
```

#### 4.6.2 Human — `PER1`

`gene_symbol="PER1"`, `keywords_list=["sleep", "circadian"]`.

```text
Classify the Human gene 'PER1' (Homo sapiens, NCBI taxon ID 9606). Allowed categories: [sleep, circadian]. You may classify into 1 category, N categories, or None (return empty list [] or None). If the text predominantly describes a non-Human ortholog/paralog/homolog, return an empty category list with rationale noting cross-species text. Return JSON with fields: gene, category (list of strings, can be empty or None), confidence (0-100), rationale (1-2 lines).

TEXT START
- Function: PER1 is a core mammalian circadian transcriptional repressor that dimerises with CRY proteins and inhibits CLOCK/BMAL1-driven transcription of clock-controlled genes.
- Phenotypes: A PER1 missense variant (S662G-equivalent in PER2 family literature) and PER1 promoter SNPs are associated with familial advanced sleep-phase syndrome and altered habitual sleep timing in genome-wide studies.

- Function: PER1 mRNA in human peripheral blood mononuclear cells exhibits a 24 h rhythm with peak expression in the early subjective day, mirroring SCN-driven oscillations.
- Phenotypes: Reduced PER1 amplitude in PBMCs is observed in shift workers and correlates with reported sleep disturbance scores.
TEXT END
```

Expected payload:

```json
{
  "gene": "PER1",
  "category": ["sleep", "circadian"],
  "confidence": 90,
  "rationale": "Human PER1 shows rhythmic PBMC expression and PER1 variants associate with familial advanced sleep-phase syndrome, supporting roles in both sleep and circadian regulation."
}
```

#### 4.6.3 Mouse — `Per1`

`gene_symbol="Per1"`, `keywords_list=["sleep", "circadian"]`.

```text
Classify the Mouse gene 'Per1' (Mus musculus, NCBI taxon ID 10090). Allowed categories: [sleep, circadian]. You may classify into 1 category, N categories, or None (return empty list [] or None). If the text predominantly describes a non-Mouse ortholog/paralog/homolog, return an empty category list with rationale noting cross-species text. Return JSON with fields: gene, category (list of strings, can be empty or None), confidence (0-100), rationale (1-2 lines).

TEXT START
- Function: Per1 is a core component of the mammalian circadian transcriptional feedback loop, expressed rhythmically in the SCN and induced by light pulses, where it is required for phase resetting of the oscillator.
- Phenotypes: Per1^-/- mice exhibit a shortened free-running period (~23.4 h) in DD, reduced NREM sleep consolidation, and attenuated delta-power rebound after 6 h sleep deprivation.

- Function: Per1 expression in cortical neurons is sleep-wake-state dependent and induced by sleep deprivation in a homeostatic manner.
- Phenotypes: Per1/Per2 double-mutant mice show severe disruption of sleep homeostasis with abolished delta-power response to 6 h sleep deprivation.
TEXT END
```

Expected payload:

```json
{
  "gene": "Per1",
  "category": ["sleep", "circadian"],
  "confidence": 92,
  "rationale": "Per1 knockouts show both circadian period changes and disrupted sleep homeostasis with attenuated delta-power rebound after sleep deprivation, supporting both labels."
}
```

After post-processing (`_normalize_category_output`), any element of
`category` outside the user-supplied keyword set would be silently dropped
— so the model cannot, for example, append `arousal` and have it survive
into the output workbook.

---

## 5. Usage accounting

Every successful Responses-API attempt is recorded into a per-gene list of
events via `_record_openai_usage_event`, with the `request_type` tag used
throughout the codebase:

| `request_type` | Source | Notes |
| --- | --- | --- |
| `summary_chunk` | `summarize_reference_function` | One per chunk attempted. Includes `chunk_index`, `total_chunks`, `chunk_words`, `paper_id`. |
| `summary_rewrite` | `rewrite_function_summary` | One per multi-chunk paper. Includes `chunk_count`. |
| `classification` | `classify_gene_from_text` | First attempt. Includes `high_quality_reference_count`, `classification_text_words`. |
| `classification_retry_short_context` | `classify_gene_from_text` | Truncated retry. |

`_summarize_openai_usage_events` aggregates these into per-gene totals
(`request_count`, `input_tokens`, `output_tokens`, `total_tokens`,
`cached_input_tokens`, `reasoning_output_tokens`) which are then written into
the run's classification output and used for cost reporting against
`Data/openai_pricing.csv`.

---

## 6. What is **not** an OpenAI call

For completeness, these adjacent components in the same pipeline never invoke
the OpenAI API:

- **Fly → human / mouse ortholog mapping**:
  `HelperScripts/get_orthologs.py` calls the DIOPT JSON API directly; the
  resulting CSVs are then re-fed into the same classification pipeline using
  `human` / `mouse` species adapters with the same prompts as above.
- **Gene symbol → ID resolution**: FlyBase TSVs for fly (plus a local
  FlyBase citation index for shared synonyms), HGNC for human, MGI for mouse —
  all local lookup tables; no OpenAI or web-API calls.
- **Reference discovery**: PubMed (metapub), Europe PMC, NCBI ELink/EFetch,
  GeneRIF, and PMC OA via `metadata_resolver.py` and the full-text fetchers.
  Europe PMC queries are organism-scoped for all three species using structured
  `ORGANISM` / `MESH` clauses (`Drosophila melanogaster`, `Homo sapiens`, or
  `Mus musculus`) before any OpenAI summarization is attempted.
- **Soft-run cost estimation** (`flai_gene_classification.py --soft-run` and
  `_plan_soft_run_gene_set`): planning only, never opens a network connection
  to OpenAI; pricing is pulled from LiteLLM's public JSON
  (`LITELLM_PRICING_URL`) into `Data/openai_pricing.csv`.
- **Excel rollup** (`summarize_gene_set.py`): pure pandas / openpyxl.

---

## 7. Quick reference — all live OpenAI requests

| # | Feature | Function | Model env var (default) | Reasoning effort env var (default) | Output cap | Schema |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Per-chunk reference summary | `summarize_reference_function` | `OPENAI_SUMMARY_MODEL` (`gpt-6-luna`) | `OPENAI_SUMMARY_REASONING_EFFORT` (`medium`) | 1200 | `FunctionPhenotypeSummary` |
| 2 | Multi-chunk paper rewrite | `rewrite_function_summary` | `OPENAI_SUMMARY_MODEL` (`gpt-6-luna`) | `OPENAI_SUMMARY_REASONING_EFFORT` (`medium`) | 1200 | `FinalFunctionPhenotypeSummary` |
| 3 | Final classification | `classify_gene_from_text` | `OPENAI_CLASSIFICATION_MODEL` (`gpt-6.2-sol`) | `OPENAI_CLASSIFICATION_REASONING_EFFORT` (`high`) | 4000 (retry: 2000) | `GeneClassification` |

All three go through `_parse_structured_completion → client.responses.parse`
with up to 3 retries and per-attempt token-usage logging.
