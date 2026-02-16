## Phenotype Matcher

A standalone tool for matching clinical phenotype descriptions to standardized ontology identifiers using Graph RAG (Retrieval-Augmented Generation).

### Quick Start

```bash
# 1. Navigate to the package directory
cd /path/to/pavs/tools/phenotype_matcher

# 2. Install with UV (recommended)
uv pip install -e .

# 3. Set your OpenRouter API key
export OPENROUTER_API_KEY="your-api-key-here"

# 4. Download MONDO ontology
mkdir -p ../../ontology
wget http://purl.obolibrary.org/obo/mondo.obo -O ../../ontology/mondo.obo

# 5. Build indices (optional but recommended for performance)
cd ../..  # Back to project root
.venv/bin/phenotype-index build --ontology hpo --model balanced

# 6. Run a simple match
.venv/bin/phenotype-matcher "severe intellectual disability" --output-format summary

# Or use uv run (from project root)
uv run phenotype-matcher "severe intellectual disability"
```

**Note:** After installation, run commands from the **project root** (`/path/to/pavs/`), not from `tools/phenotype_matcher/`.

### Overview

Phenotype Matcher maps free-text clinical descriptions to:
- **HPO** (Human Phenotype Ontology) - for phenotypic features
- **MONDO** (Disease Ontology) - for diseases
- **OMIM** (Online Mendelian Inheritance in Man) - for diseases and genes
- **OrphaNet** - for rare diseases (via MONDO cross-references)

The tool uses a sophisticated Graph RAG architecture that combines:
1. **Semantic embeddings** - for finding similar terms
2. **Ontology graph structure** - for understanding hierarchical relationships
3. **LLM reasoning** - for contextual selection and multi-phenotype extraction

### Key Features

- **Multi-phenotype extraction**: Handles descriptions with multiple phenotypes (e.g., "seizures, hypotonia, and feeding difficulties")
- **Negation detection**: Identifies excluded phenotypes (e.g., "no cardiac abnormalities")
- **Severity modifiers**: Extracts severity (e.g., "severe", "mild") as separate HPO terms
- **Dual-branch retrieval**: Separates phenotypes from severity modifiers for better accuracy
- **Configurable models**: Choose from fast, balanced, or accurate embedding and LLM models

### Installation

#### Option 1: Using UV (Recommended)

[UV](https://github.com/astral-sh/uv) is a fast Python package installer and resolver.

```bash
# Install the tool with UV
cd tools/phenotype_matcher
uv pip install -e .

# Or run directly without installing (UV will handle dependencies)
uv run phenotype-matcher "severe intellectual disability"

# Run Python scripts with UV
uv run python -c "from phenotype_matcher import PhenotypeMatcher; print('Imported successfully')"
```

#### Option 2: Using pip (Standard Python)

```bash
cd tools/phenotype_matcher
pip install -e .

# Optional: Install with progress bar support
pip install -e ".[progress]"
```

#### Dependencies

- Python >= 3.8
- pyhpo >= 3.2.0
- pronto >= 2.5.0
- sentence-transformers >= 2.2.0
- scikit-learn >= 1.0.0
- requests >= 2.28.0
- torch >= 1.13.0

#### Required Data

The tool requires HPO and MONDO ontology files:
- `ontology/hp.obo` (HPO ontology) - automatically downloaded by pyhpo on first use
- `ontology/mondo.obo` (MONDO ontology) - download from http://purl.obolibrary.org/obo/mondo.obo

```bash
# Download MONDO ontology
mkdir -p ontology
wget http://purl.obolibrary.org/obo/mondo.obo -O ontology/mondo.obo
```

### Index Management

The tool uses embedding indices for fast semantic search. You can pre-build indices for different models and ontologies.

#### Building Indices

```bash
# Build indices for all ontologies with balanced model (recommended for first-time setup)
phenotype-index build --ontology all --model balanced

# Build all combinations (HPO, MONDO, OMIM × all models)
phenotype-index build --ontology all --model all

# Build specific ontology/model combination
phenotype-index build --ontology hpo --model medical
phenotype-index build --ontology mondo --model accurate

# Rebuild (force overwrite) existing indices
phenotype-index rebuild --ontology hpo --model balanced
phenotype-index rebuild --ontology all --model all  # Rebuild everything

# Use CPU if CUDA is slow or hanging
phenotype-index build --ontology hpo --model balanced --device cpu
```

**Performance Tips:**
- **HPO** (~14,500 terms): 2-3 minutes on CPU, <1 min on GPU
- **MONDO** (~25,000 terms): 4-5 minutes on CPU, 1-2 min on GPU  
- **OMIM** (~8,000 terms): 1-2 minutes on CPU, <1 min on GPU
- Use `--device cpu` if CUDA hangs or causes issues
- Start with `fast` model, then build others as needed
- Building is one-time cost - subsequent loads are instant

#### Managing Indices

```bash
# List all existing indices
phenotype-index list

# Get detailed information about an index
phenotype-index info --ontology hpo --model balanced

# Delete an index
phenotype-index delete --ontology hpo --model fast
```

**Example output of `phenotype-index list`:**
```
Existing Indices:

HPO:
  ✓ all-MiniLM-L6-v2 (14523 terms, 6.8 MB)
  ✓ all-mpnet-base-v2 (14523 terms, 42.3 MB)
  ✓ pritamdeka_S-PubMedBert-MS-MARCO (14523 terms, 42.3 MB)

MONDO:
  ✓ all-mpnet-base-v2 (25341 terms, 73.5 MB)

OMIM:
  (none)
```

**Index Storage:**

Indices are organized by ontology and model:
```
data/graph_rag_cache/
├── hpo/
│   ├── embeddings_all-MiniLM-L6-v2.pkl
│   ├── embeddings_all-mpnet-base-v2.pkl
│   └── embeddings_pritamdeka_S-PubMedBert-MS-MARCO.pkl
├── mondo/
│   └── embeddings_all-mpnet-base-v2.pkl
└── omim/
    └── embeddings_all-mpnet-base-v2.pkl
```

**What gets indexed?**

For each ontology term, the embedding includes:
- **Primary label** (e.g., "Seizure")
- **Definition** (full clinical description)
- **ALL synonyms** (no filtering - every synonym is included)

Example embedding text for HP:0001250 (Seizure):
```
Name: Seizure. Definition: A seizure is an intermittent abnormality of the 
central nervous system due to a sudden, excessive, disorderly discharge of 
brain neurons. Synonyms: Seizures, Epileptic seizure, Convulsions.
```

Example for HP:0000303 (Mandibular prognathia) with 28 synonyms:
```
Name: Mandibular prognathia. Definition: [...]. Synonyms: Anterior 
positioned lower jaw, Increased anterior-posterior length of mandible, 
Increased forward projection of chin, Increased forward projection of jaw, 
[...all 28 synonyms included...].
```

**Synonym statistics (HPO):**
- Mean: 1.21 synonyms per term
- Median: 1 synonym per term  
- Maximum: 28 synonyms (HP:0000303)
- **100% of synonyms are indexed** for every term

This comprehensive text embedding allows the tool to match:
- Official term names ("Seizure")
- Clinical descriptions ("intermittent abnormality")
- Layperson terms ("Convulsions")
- **ALL synonyms and alternative phrasings** (no limits)

**Why pre-build indices?**
- First run without indices: 5-10 minutes (computes embeddings)
- Subsequent runs with indices: ~1 second (loads from cache)
- Pre-building allows quick switching between models
- Parallel processing benefits from pre-built indices
- **All synonyms are indexed** for comprehensive matching

### Usage

See `example_usage.py` for comprehensive examples of all features.

#### Python API

**With UV:**

```bash
# Run Python script with UV
uv run python my_script.py

# Or use UV's inline script execution
uv run python -c "
from phenotype_matcher import PhenotypeMatcher, PhenotypeInput
matcher = PhenotypeMatcher()
output = matcher.match(PhenotypeInput(text='seizures'))
print(output.get_hpo_ids())
"
```

**Python code (works with both UV and standard Python):**

```python
from phenotype_matcher import PhenotypeMatcher, PhenotypeInput

# Initialize matcher
matcher = PhenotypeMatcher()

# Match phenotypes
input_data = PhenotypeInput(text="severe intellectual disability and seizures")
output = matcher.match(input_data)

# Access results
for pheno in output.phenotypes:
    print(f"{pheno.hpo_id}: {pheno.label}")
    if pheno.severity_label:
        print(f"  Severity: {pheno.severity_label}")

# Get specific identifiers
hpo_ids = output.get_hpo_ids()  # Present phenotypes
excluded_ids = output.get_hpo_ids(excluded=True)  # Negated phenotypes
omim_ids = output.get_omim_ids()
mondo_ids = output.get_mondo_ids()
```

#### Command Line

**With UV (recommended):**

```bash
# Basic usage
uv run phenotype-matcher "severe intellectual disability"

# Multiple phenotypes (comma-separated)
uv run phenotype-matcher "seizures, feeding difficulties, hypotonia"

# Negation detection
uv run phenotype-matcher "no cardiac abnormalities, normal vision"

# Provide context for disambiguation
uv run phenotype-matcher "ASD" --context "cardiac defect"

# Different output formats
uv run phenotype-matcher "seizures" --output-format hpo_ids
uv run phenotype-matcher "seizures" --output-format tsv
uv run phenotype-matcher "seizures" --output-format summary

# Batch processing from file
uv run phenotype-matcher --input phenotypes.txt --output results.jsonl

# Multiple separate phenotypes (batch mode via command line)
uv run phenotype-matcher -b "seizures" -b "hypotonia" -b "feeding difficulties"

# Parallel processing with 4 threads (faster for large batches)
uv run phenotype-matcher --input phenotypes.txt --output results.jsonl -j 4 --progress

# Combine batch and parallelization with progress bar
uv run phenotype-matcher -b "seizures" -b "hypotonia" -b "cough" -j 2 --progress

# Use different models
uv run phenotype-matcher "global developmental delay" \
  --embedding-model medical \
  --llm-model accurate
```

**With standard Python:**

```bash
# Basic usage (after pip install -e .)
phenotype-matcher "severe intellectual disability"

# Multiple separate phenotypes (batch mode)
phenotype-matcher -b "seizures" -b "hypotonia" -b "feeding difficulties"

# Batch processing from file
phenotype-matcher --input phenotypes.txt --output results.jsonl

# Parallel processing (4 threads)
phenotype-matcher --input phenotypes.txt --output results.jsonl -j 4 --progress

# Different output formats
phenotype-matcher --input phenotypes.txt --output-format tsv --output results.tsv -j 4
phenotype-matcher -b "seizures" -b "hypotonia" --output-format hpo_ids
```

**Batch Processing Features:**

- **Multiple inputs**: Use `-b` multiple times or `--input file.txt`
- **Parallelization**: Use `-j N` to process with N threads (default: 1)
- **Progress bar**: Add `--progress` flag (requires `pip install tqdm`)
- **Thread-safe writing**: Results written immediately as completed
- **Line-by-line output**: Each result on a new line (JSONL format for JSON)

**Performance Example:**
```bash
# Sequential (1 thread) - processes one at a time
phenotype-matcher --input 100_phenotypes.txt -j 1  # ~400 seconds

# Parallel (4 threads) - 4x faster
phenotype-matcher --input 100_phenotypes.txt -j 4 --progress  # ~100 seconds
```

#### Testing

The tool includes a comprehensive test suite with difficult cases to validate the matcher:

```bash
# Run all test cases (with UV)
uv run phenotype-matcher --test

# Run all test cases (standard Python)
phenotype-matcher --test

# Run a specific test case
phenotype-matcher --test-case case_1

# List available test cases
uv run python -c "from phenotype_matcher import test_cases; print('\n'.join([f'{c[\"id\"]}: {c[\"name\"]}' for c in test_cases.get_all_test_cases()]))"
```

**Test cases include:**
- `case_1`: Complex cardiac case with typo ("Caught" → "Cough")
- `case_2`: Skeletal abnormalities with typo ("numerus" → "humerus")
- `case_3`: Multiple phenotypes with conjunctions
- `case_4`: Mixed negation and assertion
- `case_5`: Multiple phenotypes with different severity levels

All test cases are verified against `hp.obo` for ground truth accuracy.

### Configuration

#### Embedding Models

| Model | Size | Speed | Best For |
|-------|------|-------|----------|
| `fast` | 80MB | Fastest | Testing, large datasets |
| `balanced` | 420MB | Medium | General use |
| `accurate` | 1.3GB | Slower | Research, complex cases |
| `medical` | 768 dim | Medium | Medical text, PubMed-trained |
| `biobert` | 768 dim | Medium | Biomedical literature |
| `sapbert` | 768 dim | Medium | **Biomedical entity linking (RECOMMENDED for medical ontology matching)** |

**SapBERT** (Self-Alignment Pre-training for BERT) is specifically trained for biomedical entity linking and ontology matching. It excels at:
- Mapping clinical terms to HPO/MONDO/OMIM
- Handling synonyms and lexical variations
- Medical domain understanding
- **Recommended for production use with medical ontologies**

#### LLM Models

| Model | Speed | Cost | Best For |
|-------|-------|------|----------|
| `fast` | Fast | Low | **Production (recommended)** |
| `accurate` | Slow | Higher | Research, critical accuracy |
| `balanced` | Medium | Free | Testing |
| `cheap` | Fast | Lowest | Large-scale processing |

#### Python Configuration

```python
from phenotype_matcher import PhenotypeMatcher, MatcherConfig

config = MatcherConfig(
    embedding_model="sapbert",   # or "fast", "balanced", "accurate", "medical", "biobert"
    llm_model="accurate",        # or "fast", "balanced", "cheap"
    top_k_phenotype=10,          # Number of phenotype candidates
    top_k_severity=5,            # Number of severity candidates
    cache_dir="data/embeddings", # Cache directory
    device="cpu",                # or "cuda"
    debug=True                   # Enable debug logging
)

matcher = PhenotypeMatcher(config)
```

### Output Schema

#### PhenotypeOutput

```python
{
    "phenotypes": [
        {
            "hpo_id": "HP:0001250",
            "label": "Seizure",
            "excluded": false,
            "severity_id": "HP:0012828",
            "severity_label": "Severe",
            "confidence": 0.95
        },
        {
            "hpo_id": "HP:0001249",
            "label": "Intellectual disability",
            "excluded": false,
            "severity_id": "HP:0012828",
            "severity_label": "Severe",
            "confidence": 0.93
        }
    ],
    "diseases": [
        {
            "mondo_id": "MONDO:0001234",
            "mondo_label": "Epileptic encephalopathy",
            "omim_gene_ids": ["OMIM:123456"],
            "omim_phenotype_ids": ["OMIM:654321"],
            "omim_labels": ["Epileptic encephalopathy 1"],
            "orphanet_ids": ["ORPHA:1234"],
            "confidence": 0.88
        }
    ],
    "raw_input": "severe seizures and intellectual disability",
    "processing_metadata": {
        "terms_processed": 1,
        "llm_calls": 1,
        "processing_time_seconds": 2.5
    }
}
```

### Examples

#### Example 1: Single phenotype with severity

```python
input_data = PhenotypeInput(text="severe intellectual disability")
output = matcher.match(input_data)

# Output:
# phenotypes: [
#   {
#     "hpo_id": "HP:0001249",
#     "label": "Intellectual disability",
#     "excluded": false,
#     "severity_id": "HP:0012828",
#     "severity_label": "Severe"
#   }
# ]
```

#### Example 2: Multiple phenotypes

```python
input_data = PhenotypeInput(text="seizures, hypotonia, and feeding difficulties")
output = matcher.match(input_data)

# Output:
# phenotypes: [
#   {"hpo_id": "HP:0001250", "label": "Seizure", "excluded": false},
#   {"hpo_id": "HP:0001252", "label": "Hypotonia", "excluded": false},
#   {"hpo_id": "HP:0011968", "label": "Feeding difficulties", "excluded": false}
# ]
```

#### Example 3: Negation detection

```python
input_data = PhenotypeInput(text="no cardiac abnormalities, normal vision")
output = matcher.match(input_data)

# Output:
# phenotypes: [
#   {"hpo_id": "HP:0001627", "label": "Abnormal heart morphology", "excluded": true},
#   {"hpo_id": "HP:0000504", "label": "Abnormality of vision", "excluded": true}
# ]
```

#### Example 4: Context-based disambiguation

```python
# "ASD" can mean "Autism Spectrum Disorder" or "Atrial Septal Defect"
input_data = PhenotypeInput(text="ASD", context="cardiac defect, heart murmur")
output = matcher.match(input_data)

# With cardiac context, matches to:
# {"hpo_id": "HP:0001631", "label": "Atrial septal defect"}
```

### Architecture

#### Dual-Branch Retrieval

The tool searches two separate HPO branches:

1. **Phenotypic abnormality branch** (HP:0000118)
   - Contains all clinical phenotypes (~14,000+ terms)
   - Used for main phenotype mapping

2. **Severity modifier branch** (HP:0012824)
   - Contains severity terms (Mild, Moderate, Severe, Profound)
   - Searched separately to ensure severity candidates are always available

#### Processing Workflow (3-Step Architecture)

```
Input: "short femur and numerus with absent radius, ASD"
           ↓
    [1. NER - Extract Individual Terms]
           ↓
    LLM extracts: ["short femur", "short numerus", "absent radius", "ASD"]
    With metadata: negation, modifiers, etc.
           ↓
    [2. Acronym Expansion + RAG Retrieval]
           ↓
    "ASD" → Expand to ["Atrial septal defect", "Autism spectrum disorder"]
    For each term/expansion:
      ┌──────────────┬─────────────┐
      ↓              ↓              
  Phenotype      Severity
  Branch         Branch
  (HP:0000118)   (HP:0012824)
      ↓              ↓
  Top-K          Top-K
  candidates     severities
      └──────────────┬─────────────┘
                     ↓
    [3. LLM Validation + Disambiguation]
           ↓
    Context: other phenotypes, gene hint
    Specificity rule: prefer general over too-specific
    Disambiguation: cardiac context → Atrial septal defect
           ↓
    Output: HP:0003097 (Short femur), HP:0003974 (Absent radius), 
            HP:0001631 (Atrial septal defect)
```

#### Key Design Decisions

1. **NER-first extraction (NEW)**: LLM extracts individual phenotypes before RAG search. Solves complex splitting problems like "short femur and humerus" → 2 phenotypes.

2. **Acronym expansion (NEW)**: Hard-coded dictionary expands ambiguous acronyms (ASD, DD, CHD, etc.) so RAG can find all possible matches, then disambiguation selects the correct one.

3. **Specificity control (NEW)**: LLM instructed to NOT select overly specific terms (e.g., "hypotonia" should NOT match "Episodic generalized hypotonia"). General terms are acceptable.

4. **Labels only from LLM, IDs from ontology**: LLM returns labels, which are looked up in the ontology to get IDs. This prevents hallucination of non-existent identifiers.

5. **Context-based disambiguation**: Uses co-occurring phenotypes and gene hints to resolve ambiguous terms (e.g., "ASD" + cardiac phenotypes → Atrial septal defect).

6. **Multi-phenotype support**: A single clinical description can map to multiple HPO terms.

4. **File-specific splitting**: Different data sources use different separators (`,`, `|`, `/`, etc.) - configurable via `split_by` parameter.

### API Reference

See inline documentation in:
- `phenotype_matcher/schemas.py` - Data schemas
- `phenotype_matcher/matcher.py` - Main matcher class
- `phenotype_matcher/cli.py` - Command-line interface

### Performance

- **First run**: 5-10 minutes (computing embeddings for ~35,000 ontology terms)
- **Subsequent runs**: ~1 second initialization (loading from cache)
- **Per-term matching**: ~2-5 seconds (depending on LLM model)
- **Cache size**: ~100-500MB (depending on embedding model)

### Environment Variables

- `OPENROUTER_API_KEY`: Required for LLM API access (unless provided via `--api-key` or config)

**Setting the API key:**

```bash
# Linux/macOS
export OPENROUTER_API_KEY="your-api-key-here"

# Or add to your shell profile (~/.bashrc, ~/.zshrc)
echo 'export OPENROUTER_API_KEY="your-api-key-here"' >> ~/.bashrc

# Windows (PowerShell)
$env:OPENROUTER_API_KEY="your-api-key-here"

# Windows (Command Prompt)
set OPENROUTER_API_KEY=your-api-key-here
```

### Running Without Installation

You can also run the tool directly with UV without installing:

```bash
# From the project root
cd tools/phenotype_matcher

# Run directly with UV (it will handle dependencies)
uv run --with pyhpo --with pronto --with sentence-transformers \
  --with scikit-learn --with requests --with torch \
  python -m phenotype_matcher.cli "severe intellectual disability"

# Or create a simple wrapper script
cat > run_matcher.sh << 'EOF'
#!/bin/bash
cd tools/phenotype_matcher
uv run --with pyhpo --with pronto --with sentence-transformers \
  --with scikit-learn --with requests --with torch \
  python -m phenotype_matcher.cli "$@"
EOF
chmod +x run_matcher.sh

# Then use it
./run_matcher.sh "seizures" --output-format hpo_ids
```

### Troubleshooting

#### ModuleNotFoundError: No module named 'phenotype_matcher'

This error occurs when the package isn't installed correctly. **Solution:**

```bash
# 1. Navigate to the package directory
cd /path/to/pavs/tools/phenotype_matcher

# 2. Install the package
uv pip install -e .

# 3. Run commands from the PROJECT ROOT (not tools/phenotype_matcher/)
cd /path/to/pavs  # Back to project root

# 4. Use the venv binaries or uv run
.venv/bin/phenotype-matcher "seizures"
# OR
uv run phenotype-matcher "seizures"
```

**Key points:**
- Install from `tools/phenotype_matcher/` directory
- Run commands from the **project root** directory
- Use `.venv/bin/phenotype-matcher` or `uv run phenotype-matcher`

#### CUDA errors
If you get CUDA errors, the tool will automatically fall back to CPU. To force CPU:
```bash
.venv/bin/phenotype-matcher "..." --device cpu
```

#### Missing ontology files
Ensure the following files exist:
- `ontology/mondo.obo` - Download from http://purl.obolibrary.org/obo/mondo.obo
- HPO is automatically downloaded by pyhpo

#### API rate limits
If you hit OpenRouter rate limits, try:
- Using the `cheap` LLM model
- Adding delays between requests
- Caching results

### Citation

If you use this tool in your research, please cite:

```
[Citation information to be added]
```

### License

MIT License - see LICENSE file for details

### Contributing

Contributions are welcome! Please see CONTRIBUTING.md for guidelines.

### Support

For issues and questions:
- GitHub Issues: [link]
- Documentation: [link]
- Email: [contact email]
