# PAVS Normalization - Quick Start Guide

## Setup (One-time)

```bash
# 1. Set API key
export OPENROUTER_API_KEY="your_key_here"

# 2. Install dependencies
cd /home/leechuck/Public/software/pavs
uv sync

# 3. Create logs directory
mkdir -p logs
```

---

## Standard Pipeline (Fast, Production)

### Process All Files (Parallel)
```bash
uv run python scripts/process_all_phenotypes.py --parallel --overwrite
```

### Process Specific File
```bash
uv run python scripts/process_all_phenotypes.py --files ahmed-variants.tsv --overwrite
```

### Background Processing (All Files in Parallel)
```bash
# Ahmed files
nohup uv run python scripts/process_all_phenotypes.py --files ahmed-variants.tsv,ahmed-pmid28454995.tsv --overwrite > logs/ahmed.log 2>&1 &

# Other files
nohup uv run python scripts/process_all_phenotypes.py --files fawzan-variants.tsv --overwrite > logs/fawzan-variants.log 2>&1 &
nohup uv run python scripts/process_all_phenotypes.py --files marwa-variants.tsv --overwrite > logs/marwa-variants.log 2>&1 &
nohup uv run python scripts/process_all_phenotypes.py --files PMC6562004.tsv --overwrite > logs/PMC6562004.log 2>&1 &
nohup uv run python scripts/process_all_phenotypes.py --files PMC7082194.tsv --overwrite > logs/PMC7082194.log 2>&1 &
nohup uv run python scripts/process_all_phenotypes.py --files ddd-diagnoses.tsv --overwrite > logs/ddd-diagnoses.log 2>&1 &

# Monitor
tail -f logs/*.log
```

### Test (10 rows)
```bash
uv run python scripts/process_all_phenotypes.py 10
```

**Output**: `data/processed/*_normalized.tsv`

---

## GraphRAG Pipeline (High Accuracy, Research)

**Key Differences from Standard:**
- ✅ Dual-branch retrieval (Phenotypic abnormality vs Severity separate)
- ✅ Multi-phenotype support (one term → multiple HPO codes)
- ✅ File-specific splitting (handles `/`, `|`, embedded HPO IDs)
- ✅ Severity labels from ontology (not LLM - prevents hallucination)
- ✅ Pre-negation detection + LLM confirmation



### Default (Balanced Settings)
```bash
uv run python scripts/process_all_phenotypes_graphrag.py --parallel --overwrite
```

### Maximum Accuracy
```bash
uv run python scripts/process_all_phenotypes_graphrag.py \
  --embedding-model accurate \
  --llm accurate \
  --top-k 10 \
  --parallel \
  --overwrite
```

### Medical Domain (PubMed Text)
```bash
uv run python scripts/process_all_phenotypes_graphrag.py \
  --embedding-model medical \
  --top-k 12 \
  --llm fast \
  --parallel
```

### Test (5 rows)
```bash
uv run python scripts/process_all_phenotypes_graphrag.py 5
```

**Output**: `data/processed_graphrag/*_normalized.tsv`

---

## Testing

```bash
# Run test suite
uv run python scripts/test_phenotypes_processing.py
```

Tests include:
- ✅ Disambiguation logic
- ✅ Variant classification normalization
- ✅ NAN gene handling (SCN11A hallucination prevention)
- ✅ File structure (36 columns)
- ✅ Content validation (IDs, OMIM, ACMG codes)

---

## Common Options

| Option | Description | Example |
|--------|-------------|---------|
| `--files` | Process specific files | `--files ahmed-variants.tsv` |
| `--overwrite` | Delete existing and reprocess | `--overwrite` |
| `--parallel` | Process files concurrently | `--parallel` |
| `10` | Limit to N rows per file | `10` |

**GraphRAG only:**
| Option | Description | Values |
|--------|-------------|--------|
| `--embedding-model` | Sentence transformer | `fast`, `balanced`, `accurate`, `medical`, `biobert` |
| `--top-k` | Number of candidates | `5` (default), `10`, `15` |
| `--llm` | LLM for selection | `fast`, `accurate`, `balanced`, `cheap` |

---

## Monitoring Progress

```bash
# Check running processes
ps aux | grep process_all_phenotypes

# Check logs
tail -f logs/*.log

# Check output files
ls -lh data/processed/
ls -lh data/processed_graphrag/

# Count processed rows
wc -l data/processed/*_normalized.tsv
```

---

## Output Schema

36 columns including:
- `id`, `source_file`
- `sex_id`, `sex_label`, `age`
- `hpo_ids`, `hpo_labels` (phenotypes)
- `hpo_excluded_ids`, `hpo_excluded_labels` (negated)
- `hpo_severity_ids`, `hpo_severity_labels` (severity modifiers)
- `disease_omim_gene_ids`, `disease_omim_phenotype_ids`, `disease_omim_labels`
- `disease_mondo_ids`, `disease_mondo_labels`
- `gene_hgnc_ids`, `gene_symbols`
- `variant_hgvs`, `variant_validated`, `variant_classification`, `variant_evidence`
- `publication_id`
- `genotype_id`, `genotype_label`
- **`solved`** (true/false - case solved status)

---

## Which Method to Use?

| Your Need | Use This |
|-----------|----------|
| **Fast, production** | Standard pipeline with `--parallel` |
| **Highest accuracy** | GraphRAG with `--embedding-model accurate --top-k 10` |
| **Medical literature** | GraphRAG with `--embedding-model medical` |
| **Testing** | Either with limit: `10` |
| **Resume interrupted** | Just re-run same command (auto-resumes) |
