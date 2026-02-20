# Phenotype Matcher v2

Maps free-text clinical phenotype descriptions to standardised ontology identifiers: **HPO**, **MONDO**, **OMIM**, and **Orphanet**.

The tool handles noisy clinical language: misspellings, inflected forms, abbreviations, negation (`"no seizures"`, `"ruled out"`), severity modifiers (`"severe"`, `"mild"`), and coordination structures (`"short and malformed fingers"`).

---

## Requirements

### Ontology files

| File | Description | Default path |
|---|---|---|
| `hp.obo` | HPO ontology | `ontology/hp.obo` |
| `mondo.obo` | MONDO ontology | `ontology/mondo.obo` |
| `phenotype.hpoa` | HPO disease associations (OMIM + Orphanet) | `data/phenotype.hpoa` |

These files are not bundled. Download them from:
- HPO: https://hpo.jax.org/app/data/ontology
- MONDO: https://mondo.monarchinitiative.org/
- phenotype.hpoa: https://hpo.jax.org/app/data/annotations

### Environment

```bash
export OPENROUTER_API_KEY="your_key_here"
```

An OpenRouter API key is required for LLM-based disambiguation and validation. A free-tier key is sufficient with the default `deepseek` model. Get one at https://openrouter.ai/.

The SapBERT embedding model (~400 MB) is downloaded from HuggingFace on first use and cached by sentence-transformers. Subsequent runs use the disk cache. A GPU is strongly recommended for building the embedding index; matching at inference time runs on CPU.

### Installation

```bash
cd tools/phenotype_matcher_v2
pip install -e .
# or
uv pip install -e .
```

---

## Quick start

```python
from phenotype_matcher_v2 import PhenotypeMatcher, MatcherConfig

cfg = MatcherConfig(
    hpo_path="ontology/hp.obo",
    mondo_path="ontology/mondo.obo",
    hpoa_path="data/phenotype.hpoa",
)
matcher = PhenotypeMatcher(cfg)

out = matcher.match("severe intellectual disability, no seizures")

print(out.get_hpo_ids())
# ['HP:0001249']   (Intellectual disability — present)

print([(p.hpo_id, p.label) for p in out.get_excluded_phenotypes()])
# [('HP:0001250', 'Seizure')]   (negated by "no")

print([(p.hpo_id, p.severity_id) for p in out.get_present_phenotypes()])
# [('HP:0001249', 'HP:0012828')]  (severe)
```

---

## Use cases

### 1. Normalise a comma-separated phenotype list

```python
out = matcher.match("hypotonia, ataxia, feeding difficulties, absent speech")

for p in out.phenotypes:
    print(p.hpo_id, p.label, "excluded=" + str(p.excluded))
```

The tokeniser hard-segments on commas, semicolons, and full stops, then matches each term independently through the full pipeline.

### 2. Detect negated phenotypes

```python
out = matcher.match("seizures ruled out, hypotonia present")

present  = out.get_hpo_ids()            # non-excluded
excluded = out.get_hpo_ids(excluded=True)
```

Negation is detected by searching a five-word window around each matched span for triggers including: `no`, `not`, `without`, `absence of`, `ruled out`, `negative for`, `denied`, and others. The phenotype is returned with `excluded=True` rather than omitted, so downstream tools can distinguish "absent" from "not mentioned".

### 3. Capture severity modifiers

```python
out = matcher.match("profound intellectual disability, mild hypotonia")

for p in out.phenotypes:
    print(p.hpo_id, p.label, "severity:", p.severity_label, p.severity_id)
# HP:0001249  Intellectual disability  severity: Profound  HP:0012829
# HP:0001252  Hypotonia               severity: Mild      HP:0012825
```

Severity modifiers are descendants of `HP:0012824` (Severity) in the HPO. The modifier closest to the phenotype span in the surrounding word window is selected.

### 4. Match disease names to MONDO / OMIM / Orphanet

```python
out = matcher.match("Dravet syndrome")

print(out.get_mondo_ids())     # ['MONDO:0100079', 'MONDO:0100135']
print(out.get_orphanet_ids())  # ['ORPHA:33069']

out2 = matcher.match("phenylketonuria")

print(out2.get_omim_ids())     # ['OMIM:261600']
print(out2.get_mondo_ids())    # ['MONDO:0009861']
```

Disease names are loaded from `phenotype.hpoa` (OMIM/Orphanet) and `mondo.obo`. Orphanet IDs are also extracted from MONDO xrefs, so a single disease name can return IDs across all three namespaces simultaneously.

### 5. Match informal or misspelled terms

```python
# Informal term — matched via SapBERT semantic similarity
out = matcher.match("fits and falls")
print(out.get_hpo_ids())   # seizure-related HPO terms

# Misspelling — matched via edit-distance-1 fuzzy lookup
out = matcher.match("Seizuree")
print(out.get_hpo_ids())   # ['HP:0001250']
```

Phase 2e uses SapBERT embeddings to find semantically similar ontology labels when no lexical match exists. Phase 1c / 2c use Levenshtein distance 1 to catch single-character typos, with LLM validation to filter false positives.

### 6. Process clinical notes in bulk

```python
import json

notes = [
    "Patient presents with seizures and hypotonia.",
    "No intellectual disability. Mild ataxia.",
    "Dravet syndrome diagnosed at age 6 months.",
]

results = []
for note in notes:
    out = matcher.match(note)
    results.append(out.to_dict())

with open("results.json", "w") as f:
    json.dump(results, f, indent=2)
```

The matcher is initialised once and is thread-safe for concurrent `match()` calls. For high-throughput use, instantiate a single `PhenotypeMatcher` and reuse it.

### 7. Extract HPO terms for phenopacket generation

```python
from phenopackets import Phenopacket, PhenotypicFeature, OntologyClass

out = matcher.match("severe intellectual disability, absent speech, no seizures")

features = []
for p in out.phenotypes:
    feature = PhenotypicFeature(
        type=OntologyClass(id=p.hpo_id, label=p.label),
        excluded=p.excluded,
    )
    if p.severity_id:
        feature.severity.CopyFrom(OntologyClass(id=p.severity_id, label=p.severity_label))
    features.append(feature)
```

### 8. Adjust matching aggressiveness

```python
from phenotype_matcher_v2 import MatcherConfig, PhenotypeMatcher

# High precision: raise the ANN similarity floor, shrink the relative window
cfg = MatcherConfig(
    hpo_path="ontology/hp.obo",
    theta_abs=0.80,    # stricter absolute floor (default 0.70)
    delta=0.02,        # tighter relative window (default 0.05)
    ann_top_k=3,       # fewer ANN candidates (default 5)
)

# High recall: lower the floor, widen the window
cfg = MatcherConfig(
    hpo_path="ontology/hp.obo",
    theta_abs=0.60,
    delta=0.10,
    ann_top_k=10,
)
```

### 9. Use a different LLM or embedding model

```python
cfg = MatcherConfig(
    hpo_path="ontology/hp.obo",
    embedding_model="fast",      # all-MiniLM-L6-v2 (fast, lower accuracy)
    llm_model="accurate",        # anthropic/claude-sonnet-4-5
)

# Or pass a fully-qualified model string
cfg = MatcherConfig(
    hpo_path="ontology/hp.obo",
    embedding_model="pritamdeka/S-PubMedBert-MS-MARCO",
    llm_model="openai/gpt-4o",
)
```

Available embedding model presets: `sapbert` (default), `fast`, `balanced`, `accurate`, `medical`.
Available LLM model presets: `deepseek` (default, free), `gpt4oss`, `accurate`, `gemini`.

### 10. Command-line usage

See the [Command-line interface](#command-line-interface) section below for full documentation.

---

## Configuration reference

All parameters are on `MatcherConfig`:

| Parameter | Default | Description |
|---|---|---|
| `hpo_path` | `"ontology/hp.obo"` | Path to HPO OBO file |
| `mondo_path` | `"ontology/mondo.obo"` | Path to MONDO OBO file |
| `hpoa_path` | `"data/phenotype.hpoa"` | Path to HPO annotations file |
| `embedding_model` | `"sapbert"` | Embedding model preset or HuggingFace name |
| `llm_model` | `"deepseek"` | LLM model preset or OpenRouter model string |
| `theta_abs` | `0.70` | ANN absolute similarity floor |
| `delta` | `0.05` | ANN relative margin from top score |
| `ann_top_k` | `5` | Maximum ANN candidates per token |
| `det_neg_win` | `5` | Negation detection window in words |
| `det_mod_win` | `4` | Modifier detection window in words |
| `cache_dir` | `"data/phenotype_matcher_v2_cache"` | Directory for embedding cache |
| `device` | `"cpu"` | Torch device for embeddings (`"cpu"` or `"cuda"`) |
| `api_key` | `None` | OpenRouter API key (falls back to `OPENROUTER_API_KEY`) |
| `debug` | `False` | Print progress messages |

---

## Output schema

`PhenotypeOutput` has two lists and helper methods:

```python
out.phenotypes          # List[PhenotypeMatch]
out.diseases            # List[DiseaseMatch]
out.get_hpo_ids()       # List[str]  — present HPO IDs
out.get_hpo_ids(excluded=True)  # List[str] — excluded HPO IDs
out.get_present_phenotypes()    # List[PhenotypeMatch]
out.get_excluded_phenotypes()   # List[PhenotypeMatch]
out.get_mondo_ids()     # List[str]
out.get_omim_ids()      # List[str]
out.get_orphanet_ids()  # List[str]
out.to_dict()           # JSON-serialisable dict
```

`PhenotypeMatch` fields:

| Field | Type | Description |
|---|---|---|
| `hpo_id` | `str` | HPO identifier, e.g. `"HP:0001250"` |
| `label` | `str` | Canonical label, e.g. `"Seizure"` |
| `excluded` | `bool` | `True` if negated in the input text |
| `severity_id` | `Optional[str]` | HPO severity term ID, e.g. `"HP:0012828"` |
| `severity_label` | `Optional[str]` | Severity label, e.g. `"Severe"` |
| `confidence` | `float` | Always `1.0` in v2 |
| `matched_by` | `str` | Always `"v2"` |

`DiseaseMatch` fields:

| Field | Type | Description |
|---|---|---|
| `mondo_id` | `Optional[str]` | MONDO identifier |
| `mondo_label` | `Optional[str]` | MONDO canonical label |
| `omim_phenotype_ids` | `List[str]` | OMIM phenotype IDs |
| `omim_labels` | `List[str]` | OMIM disease names |
| `orphanet_ids` | `List[str]` | Orphanet IDs |
| `confidence` | `float` | Always `1.0` in v2 |
| `matched_by` | `str` | Always `"v2"` |

---

## Performance notes

- **First run**: loading HPO + MONDO via pronto and building SapBERT embeddings takes 10–30 minutes depending on hardware. Embeddings are cached to `cache_dir` and reused on subsequent runs.
- **Subsequent runs**: ontology loading takes ~30 seconds; matching is fast (milliseconds per term for lexical paths, ~50 ms per token if the ANN path is reached).
- **LLM calls**: each untrusted match (stemmed, fuzzy, ANN) makes one `llm_val_batch` call. The default DeepSeek model on the free OpenRouter tier has rate limits; for bulk processing, consider a paid tier or an alternative model.
- **GPU**: recommended for embedding index construction. Set `device="cuda"` in `MatcherConfig`. Inference (`match()` calls) runs on CPU.

---

## Command-line interface

After installation (`pip install -e .`), the `phenotype-matcher` command is available.

```
phenotype-matcher {match,batch,info,check} [options]
```

All subcommands share the same ontology path and model arguments (see [Common options](#common-options) below).

---

### `phenotype-matcher match`

Match one or more text strings to ontology IDs.

```
phenotype-matcher match [TEXT ...] [--input FILE] [--format json|tsv|pretty]
                        [--present-only] [--hpo-only] [--diseases-only]
                        [common options]
```

**Input sources** (can be combined):
- Positional arguments: one or more quoted strings
- `--input FILE`: one text per line from a file
- `--input -`: one text per line from stdin

**Output formats:**

| Format | Description |
|---|---|
| `json` | JSON object (or array for multiple inputs). Default. |
| `tsv` | Tab-separated: `input`, `hpo_present`, `hpo_excluded`, `severity`, `mondo_ids`, `omim_ids`, `orphanet_ids`. IDs are semicolon-separated. Severity as `HPO_ID=SEVERITY_ID` pairs. |
| `pretty` | Human-readable text with labelled sections. |

**Output filters:**

| Flag | Effect |
|---|---|
| `--present-only` | Omit excluded (negated) phenotypes |
| `--hpo-only` | Omit disease matches |
| `--diseases-only` | Omit phenotype matches |

**Examples:**

```bash
# Single term, default JSON output
phenotype-matcher match "seizures, hypotonia"

# Multiple terms, human-readable
phenotype-matcher match "no fever" "severe ataxia" --format pretty

# File input, TSV output, redirect to file
phenotype-matcher match --input terms.txt --format tsv > results.tsv

# Stdin pipeline
echo "intellectual disability" | phenotype-matcher match --input -

# HPO IDs only, no diseases, no excluded terms
phenotype-matcher match "seizures, no fever" --hpo-only --present-only --format tsv

# Fast model, cheaper LLM, GPU
phenotype-matcher match "ataxia" --model fast --llm gpt4oss --device cuda

# Non-default ontology paths
phenotype-matcher match "hypotonia" \
  --hpo /data/ontology/hp.obo \
  --mondo /data/ontology/mondo.obo \
  --hpoa /data/phenotype.hpoa \
  --cache-dir /data/cache
```

**TSV output columns:**

| Column | Content |
|---|---|
| `input` | Original input text |
| `hpo_present` | Semicolon-separated HPO IDs of present phenotypes |
| `hpo_excluded` | Semicolon-separated HPO IDs of excluded (negated) phenotypes |
| `severity` | Semicolon-separated `HPO_ID=SEVERITY_ID` pairs |
| `mondo_ids` | Semicolon-separated MONDO IDs |
| `omim_ids` | Semicolon-separated OMIM IDs |
| `orphanet_ids` | Semicolon-separated Orphanet IDs |

---

### `phenotype-matcher batch`

Process a TSV or CSV file, matching the text in a named column. All original columns are preserved and six new columns are appended.

```
phenotype-matcher batch INPUT [--column COL] [--output FILE]
                        [--sep tab|comma] [common options]
```

| Argument | Description |
|---|---|
| `INPUT` | Input file (TSV or CSV) |
| `--column COL` | Column name containing phenotype text (required) |
| `--output FILE` | Output file (default: stdout) |
| `--sep tab\|comma` | Input delimiter (default: tab) |

The output is always TSV regardless of the input delimiter.

Progress is reported to stderr every 10 rows, or on every row with `--debug`.

**Examples:**

```bash
# Basic batch run
phenotype-matcher batch patients.tsv --column phenotypes --output results.tsv

# CSV input
phenotype-matcher batch cohort.csv --sep comma --column clinical_description \
  --output cohort_normalised.tsv

# With faster models for large batches
phenotype-matcher batch large_dataset.tsv --column findings \
  --model fast --llm gpt4oss --output large_results.tsv

# Stdout (pipe to downstream tool)
phenotype-matcher batch patients.tsv --column phenotypes | \
  python filter_hpo.py > filtered.tsv
```

**Added output columns** (same as `match --format tsv`):

`hpo_present`, `hpo_excluded`, `severity`, `mondo_ids`, `omim_ids`, `orphanet_ids`

Rows with an empty phenotype text column pass through unchanged with empty values in the new columns.

---

### `phenotype-matcher info`

Look up a term by ontology ID or label string. Shows the canonical label, all synonyms, namespace, and cross-references.

```
phenotype-matcher info QUERY [common options]
```

`QUERY` is either:
- An ontology ID: `HP:0001250`, `MONDO:0100079`, `OMIM:261600`, `ORPHA:33069`
- A label or partial label: `"seizure"`, `"dravet"`, `"intellectual"`

For ID queries, synonyms and namespace information are printed. For label queries, all matching term IDs are listed. For partial label queries (no exact match), up to 20 labels containing the query string are shown.

**Examples:**

```bash
# Look up by HPO ID
phenotype-matcher info HP:0001250

# Look up by MONDO ID (shows Orphanet cross-references)
phenotype-matcher info MONDO:0100079

# Look up by OMIM ID
phenotype-matcher info OMIM:261600

# Exact label lookup (lists all matching IDs)
phenotype-matcher info "seizure"

# Partial label search
phenotype-matcher info "dravet"
phenotype-matcher info "intellectual"
```

**Example output (ID query):**

```
ID     : HP:0001250
Label  : Seizure
Synonyms (3):
  convulsion
  fit
  seizures
Namespace: HPO phenotype (descendant of HP:0000118)
```

**Example output (label query):**

```
Exact match for 'seizure':
  HP:0001250              Seizure
```

---

### `phenotype-matcher check`

Verify that all required files, environment variables, and Python packages are present.

```
phenotype-matcher check [common options]
```

Reports:
- Presence and size of each ontology file
- Whether `OPENROUTER_API_KEY` is set
- Whether the embedding cache for the selected model exists
- Version of each installed Python package

Exit code is 0 if everything is in order, 1 if any required item is missing.

```bash
# Check default configuration
phenotype-matcher check

# Check with non-default paths and a specific model
phenotype-matcher check \
  --hpo /data/hp.obo \
  --mondo /data/mondo.obo \
  --hpoa /data/phenotype.hpoa \
  --model fast
```

**Example output:**

```
Ontology files:
  [OK]  hp.obo   (HPO): ontology/hp.obo  (52,847 KB)
  [OK]  mondo.obo (MONDO): ontology/mondo.obo  (98,231 KB)
  [OK]  phenotype.hpoa (OMIM/Orphanet): data/phenotype.hpoa  (35,012 KB)

API key:
  [OK]  OPENROUTER_API_KEY set  (sk-or-v1...)

Embedding cache:
  [OK]  sapbert (cambridgeltl/SapBERT-from-PubMedBERT-fulltext): cached  (1,204 MB)

Python packages:
  [OK]  pronto  (2.5.8)
  [OK]  pyahocorasick  (2.1.0)
  [OK]  rapidfuzz  (3.9.3)
  [OK]  spacy  (3.7.4)
  [OK]  sentence-transformers  (3.0.1)
  [OK]  numpy  (1.26.4)
  [OK]  requests  (2.31.0)
  [OK]  torch  (2.3.0)

Configuration OK.
```

---

### Common options

All subcommands accept these options:

**Ontology paths:**

| Option | Default | Description |
|---|---|---|
| `--hpo PATH` | `ontology/hp.obo` | Path to HPO OBO file |
| `--mondo PATH` | `ontology/mondo.obo` | Path to MONDO OBO file |
| `--hpoa PATH` | `data/phenotype.hpoa` | Path to HPO disease annotations |

**Model selection:**

| Option | Default | Description |
|---|---|---|
| `--model` | `sapbert` | Embedding model: `sapbert`, `fast`, `balanced`, `accurate`, `medical` |
| `--llm` | `deepseek` | LLM model: `deepseek`, `gpt4oss`, `accurate`, `gemini` |
| `--device` | `cpu` | Torch device: `cpu` or `cuda` |

**Runtime:**

| Option | Default | Description |
|---|---|---|
| `--cache-dir DIR` | `data/phenotype_matcher_v2_cache` | Embedding cache directory |
| `--api-key KEY` | — | OpenRouter API key (overrides `OPENROUTER_API_KEY`) |
| `--debug` | off | Print progress messages to stderr |
