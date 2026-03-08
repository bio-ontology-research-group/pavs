# GraphRAG vs Standard Normalization: Comparison Guide

## Overview

PAVS offers two normalization pipelines optimized for different use cases:

| Feature | Standard Pipeline | GraphRAG Pipeline |
|---------|------------------|-------------------|
| **Script** | `process_all_phenotypes.py` | `process_all_phenotypes_graphrag.py` |
| **Method** | Fuzzy matching + LLM validation | Semantic embeddings + graph + LLM |
| **Speed** | ~16s per row | Slower (embedding overhead) |
| **Accuracy** | High (validated) | Highest (semantic + structural) |
| **Best For** | Production, large datasets | Research, complex phenotypes |
| **Output Dir** | `data/processed/` | `data/processed_graphrag/` |

---

## Standard Pipeline

### How It Works

1. **Candidate Generation**: RapidFuzz fuzzy string matching against HPO/MONDO names
2. **LLM Extraction**: LLM maps terms to ontologies with fuzzy candidates as hints
3. **LLM Validation**: Second LLM call strictly validates each extracted term against source text
4. **Filtering**: Remove hallucinations, over-specific terms, invalid codes

### Key Features

✅ **Two-step validation** prevents hallucinations  
✅ **Fast** - Provider routing (Clarifai, Chutes, Google Vertex AI)  
✅ **Production-ready** - Tested on 7 diverse datasets  
✅ **Incremental** - Line-by-line writing, auto-resume  
✅ **Parallel files** - Process up to 7 files simultaneously  

### Configuration

- **Model**: `openai/gpt-oss-120b`
- **Providers**: Clarifai → Chutes → Google Vertex AI
- **Retry**: Up to 3 attempts for failed API calls

### When to Use

- ✅ Production pipelines requiring speed
- ✅ Large datasets (thousands of cases)
- ✅ Well-formatted clinical text
- ✅ Standard phenotype terminology

### Example Commands

```bash
# Process all files in parallel
uv run python scripts/process_all_phenotypes.py --parallel

# Process specific files in background
nohup uv run python scripts/process_all_phenotypes.py --files ahmed-pmid28454995.tsv --overwrite > logs/ahmed.log 2>&1 &

# Test with 10 rows
uv run python scripts/process_all_phenotypes.py 10
```

---

## GraphRAG Pipeline

### How It Works

1. **Semantic Retrieval**: 
   - Encode clinical term using sentence transformer
   - Compare against 35,000+ pre-computed embeddings
   - Retrieve top-K semantically similar candidates
   
2. **Graph Expansion**:
   - Add parent terms from ontology hierarchy
   - Include definitions and synonyms
   - Build rich structural context
   
3. **LLM Selection**:
   - LLM sees original term + candidates + graph structure
   - Selects best match considering semantics AND structure
   - Handles negations and severity modifiers

### Key Features

✅ **Semantic understanding** - Captures meaning, not just strings  
✅ **Graph-aware** - Considers ontology relationships  
✅ **Highly customizable** - Multiple embedding and LLM models  
✅ **Tunable** - Adjust top-K for precision/recall trade-off  
✅ **Medical domain models** - Specialized PubMed/BioNLP embeddings  

### Configuration

**Embedding Models** (5 options):
- `fast` - all-MiniLM-L6-v2 (80MB, fastest)
- `balanced` - all-mpnet-base-v2 (420MB, **recommended**)
- `accurate` - BAAI/bge-large-en-v1.5 (1.3GB, best)
- `medical` - S-PubMedBert-MS-MARCO (PubMed-trained)
- `biobert` - biobert-base-cased-v1.2 (biomedical)

**LLM Models** (4 options):
- `fast` - openai/gpt-oss-120b (good balance)
- `accurate` - anthropic/claude-3.5-sonnet (highest quality)
- `balanced` - google/gemini-2.0-flash-exp:free
- `cheap` - deepseek/deepseek-chat

**Top-K**: Adjustable (default: 5, recommended: 5-15)

### When to Use

- ✅ Research datasets requiring highest accuracy
- ✅ Complex or ambiguous clinical descriptions
- ✅ Rare phenotypes with many synonyms
- ✅ Medical literature text (use `--embedding-model medical`)
- ✅ Quality over speed

### Example Commands

```bash
# Balanced setup (recommended for research)
uv run python scripts/process_all_phenotypes_graphrag.py \
  --embedding-model balanced \
  --top-k 10

# Maximum accuracy setup
uv run python scripts/process_all_phenotypes_graphrag.py \
  --embedding-model accurate \
  --llm accurate \
  --top-k 15

# Medical domain specialized
uv run python scripts/process_all_phenotypes_graphrag.py \
  --embedding-model medical \
  --top-k 12 \
  --llm fast

# Fast testing
uv run python scripts/process_all_phenotypes_graphrag.py 10 --embedding-model fast

# Parallel processing
uv run python scripts/process_all_phenotypes_graphrag.py --parallel --overwrite
```

---

## Technical Comparison

### Candidate Generation

| Method | Standard | GraphRAG |
|--------|----------|----------|
| **Algorithm** | RapidFuzz token_sort_ratio | Cosine similarity on embeddings |
| **Candidates** | 10 per term | Configurable (5-15) |
| **Advantages** | Fast, no GPU needed | Semantic understanding |
| **Limitations** | String-based only | Requires embedding computation |

### Accuracy Mechanisms

| Mechanism | Standard | GraphRAG |
|-----------|----------|----------|
| **Hallucination prevention** | Two-step LLM validation | Graph structure constrains choices |
| **Ambiguity resolution** | Context-based rules + gene info | Semantic similarity + definitions |
| **Negation detection** | LLM in both steps | LLM with graph context |
| **Severity extraction** | LLM validation against HP:0012824 | LLM with severity candidates |

### Performance

**Standard Pipeline** (Case 6 test):
- Input: "feeding difficulties, noisy breathing, ventricular arrhythmia, cardiac arrest"
- Output: 4/4 phenotypes correct
- Hallucinations: 0
- Time: ~16s per row

**GraphRAG Pipeline** (estimated):
- Input: Same
- Output: 4/4 phenotypes correct (with definitions)
- Hallucinations: 0 (graph-constrained)
- Time: ~25-40s per row (depends on embedding model)

---

## Decision Matrix

| Your Situation | Recommended Pipeline | Configuration |
|----------------|---------------------|---------------|
| Processing 10,000+ cases | **Standard** | `--parallel` |
| Publishing research paper | **GraphRAG** | `--embedding-model accurate --top-k 10` |
| Clinical text from PubMed | **GraphRAG** | `--embedding-model medical --top-k 12` |
| Need results in 1 hour | **Standard** | Default settings |
| Maximum quality, time flexible | **GraphRAG** | `--embedding-model accurate --llm accurate --top-k 15` |
| Testing/development | **Standard** | `10` (limit to 10 rows) |
| Mixed rare/common phenotypes | **GraphRAG** | `--embedding-model balanced --top-k 10` |

---

## Common Workflows

### Production Deployment
```bash
# Standard method, parallel processing, all files
uv run python scripts/process_all_phenotypes.py --parallel --overwrite
```

### Research Publication
```bash
# GraphRAG with high accuracy
uv run python scripts/process_all_phenotypes_graphrag.py \
  --embedding-model accurate \
  --llm accurate \
  --top-k 10 \
  --parallel
```

### Quick Testing
```bash
# Standard method, 5 rows per file
uv run python scripts/process_all_phenotypes.py 5
```

### Resuming Interrupted Job
```bash
# Both methods auto-resume - just re-run the same command
uv run python scripts/process_all_phenotypes.py --parallel
# OR
uv run python scripts/process_all_phenotypes_graphrag.py --parallel
```

---

## Combining Both Methods

You can use both pipelines and compare results:

```bash
# Run standard method
uv run python scripts/process_all_phenotypes.py --parallel

# Run GraphRAG method  
uv run python scripts/process_all_phenotypes_graphrag.py --parallel --embedding-model balanced

# Compare outputs
diff data/processed/ahmed_normalized.tsv data/processed_graphrag/ahmed_normalized.tsv
```

This allows you to:
- Validate results across methods
- Use standard for speed, GraphRAG for quality checks
- Identify cases where methods disagree (may need manual review)
