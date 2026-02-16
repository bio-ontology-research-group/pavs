# PAVS Normalization Methods

This document describes the methodology used to normalize raw phenotypic data from various sources (CSV/TSV files) into a standardized format compatible with the PAVS database. The process involves multiple steps: data ingestion, entity extraction, ontology mapping, disambiguation, and validation.

## Standalone Phenotype Matching Tool

A reusable, standalone phenotype matching tool has been developed and is available at `tools/phenotype_matcher/`. This tool can be used independently of the PAVS pipeline to map clinical phenotype descriptions to ontology identifiers.

**Tool Capabilities:**
- Maps free-text phenotype descriptions to HPO, MONDO, OMIM, and OrphaNet identifiers
- Detects negation (excluded phenotypes)
- Extracts severity modifiers
- Handles multiple phenotypes in a single description
- Uses Graph RAG architecture (semantic embeddings + graph structure + LLM reasoning)

**Installation:**
```bash
cd tools/phenotype_matcher
pip install -e .
```

**Usage Examples:**
```python
# Python API
from phenotype_matcher import PhenotypeMatcher, PhenotypeInput

matcher = PhenotypeMatcher()
output = matcher.match(PhenotypeInput(text="severe intellectual disability"))
print(output.get_hpo_ids())  # Present phenotypes
print(output.get_hpo_ids(excluded=True))  # Negated phenotypes
```

```bash
# Command line
phenotype-matcher "seizures, hypotonia, feeding difficulties"
phenotype-matcher "no cardiac abnormalities" --output-format hpo_ids
```

For complete documentation, see `tools/phenotype_matcher/README.md`.

---

## 1. Resources & Ontologies

The normalization pipeline leverages several key resources and ontologies:

*   **Human Phenotype Ontology (HPO):** Used for mapping phenotypic features (e.g., "Seizures" -> `HP:0001250`).
*   **Online Mendelian Inheritance in Man (OMIM):** Used for mapping diseases (`OMIM:XXXXXX`) and genes (`OMIM:XXXXXX`).
*   **MONDO Disease Ontology:** Used for high-level disease classification (`MONDO:XXXXXXX`).
*   **NCBI Gene:** Used for gene symbol mapping (HGNC, Entrez ID) via `Homo_sapiens.gene_info`.
*   **Genomic Coordinates (GFF3):** Used for validating variants against genomic positions.
*   **LLM (Large Language Model):** `openai/gpt-oss-120b` (via OpenRouter) is used for complex text-to-ontology mapping, negation detection, and severity extraction.
*   **Sentence Transformers (GraphRAG only):** Semantic embeddings for similarity-based candidate retrieval.

## 2. Data Ingestion & Preprocessing

The script `scripts/process_all_phenotypes.py` handles the ingestion of diverse input formats. It standardizes column names and basic metadata:

*   **ID Normalization:** Case IDs are converted to string format, stripping trailing `.0` (e.g., "1.0" -> "1").
*   **Demographics:**
    *   **Sex:** Mapped to NCIT terms (Male: `NCIT:C20197`, Female: `NCIT:C16576`).
    *   **Age:** Parsed from ISO 8601 duration strings (e.g., "P5Y") or free text (e.g., "5y") using regex.
    *   **Ancestry:** Default mapping to Middle Eastern (`HANCESTRO:0852`) for this dataset.
    *   **Country:** Default mapping to Saudi Arabia (`GAZ:00005279`).
    *   **Consanguinity:** Mapped to HPO terms (e.g., "Yes" -> `HP:0025820`).
    *   **Family History:** Mapped to HPO terms (e.g., "Positive" -> `HP:0032316`).
    *   **Genotype:** Mapped to GENO terms (e.g., "Homozygous" -> `GENO:0000136`).

## 3. Variant & Gene Normalization

Genetic data is extracted from various column formats and normalized:

*   **Gene Symbols:** Mapped to HGNC/Entrez IDs using a local lookup table built from NCBI Gene info.
*   **Variants:**
    *   HGVS strings (e.g., `c.123A>G`) are extracted using regex or LLM assistance for complex descriptions.
    *   **Classification:** Terms like "Ambiguous (VOUS)", "VUS", "Uncertain" are normalized to ACMG standard "Uncertain significance". "Pathogenic variant" becomes "Pathogenic".
*   **Validation:** Variant coordinates are checked against the GFF3 file to ensure they fall within valid gene regions.

## 4. Phenotype & Disease Normalization

This is the most complex part of the pipeline, handling free-text descriptions of phenotypes and diagnoses.

### 4.1. Text Extraction & Disambiguation
Raw text from "Phenotype", "Diagnosis", and "OMIM" columns is combined. Before mapping, ambiguous terms are resolved:

*   **Context-Aware Disambiguation:**
    *   **"ASD"**: Disambiguated between "Autism Spectrum Disorder" and "Atrial Septal Defect" based on co-occurring terms (e.g., "heart", "murmur" -> Cardiac; "development", "speech" -> Neuro).
    *   **Gene Context:** If available, the associated gene is checked against HPO annotations to see if it is linked to one of the candidate terms.

### 4.2. Standard Method: LLM-Based Mapping with Two-Step Validation

The standard pipeline uses a two-step approach to prevent hallucinations:

**Step 1: Extraction**
- Disambiguated text is split into phenotype terms (diagnosis handled separately)
- Terms sent to LLM (`openai/gpt-oss-120b`) with fuzzy-matched HPO candidates
- LLM extracts phenotypes with flexible interpretation

**Step 2: Validation** 
- Second LLM call validates each extracted phenotype against original text
- Strict filtering: only keeps terms explicitly mentioned or direct synonyms
- Removes hallucinations, over-specific terms, and inferred phenotypes

**LLM Instructions:**
1.  **Map to HPO:** Identify phenotypic abnormalities.
2.  **Detect Negation:** Flag terms as "excluded" if the text says "no seizures" or "excludes X".
3.  **Extract Severity:** Identify severity modifiers (e.g., "Severe", "Mild") and map them to HPO severity terms (e.g., `HP:0012828`).
4.  **Map to OMIM/MONDO:** Identify specific diseases or syndromes.
5.  **Validate:** Confirm each extracted term is actually in source text.

### 4.3. GraphRAG Method: Semantic Retrieval + Graph Context + LLM

The GraphRAG pipeline uses a **dual-branch retrieval** approach with three main steps:

#### Branch-Specific Retrieval

GraphRAG searches TWO separate HPO branches:

1. **Phenotypic abnormality branch** (HP:0000118)
   - Contains all clinical phenotypes (seizures, feeding difficulties, etc.)
   - ~14,000+ terms
   - Used for main phenotype mapping

2. **Severity modifier branch** (HP:0012824)
   - Contains severity terms (Mild, Moderate, Severe, Profound, etc.)
   - ~10+ terms
   - Searched separately to detect severity modifiers
   - **Note**: Severity is NEVER in the phenotype branch - it's a separate concept

#### File-Specific Splitting Logic

Different source files use different separators:

| File | Separators | Special Handling |
|------|------------|------------------|
| **ahmed-pmid28454995** | `,` | Standard comma split |
| **ahmed-variants** | `,` and `\|` | Strip embedded HPO IDs like `(HP:0001263)` |
| **marwa-variants** | `\|` | Strip embedded HPO IDs |
| **fawzan-variants** | `,` and `/` | `/` = alternatives (e.g., "intolerance/fatigue") |
| **PMC6562004** | `,` | Standard comma split |
| **Others** | `,` `\|` `;` | Default multi-separator split |

#### Processing Workflow

**Step 1: Dual Semantic Retrieval**
- Clinical term encoded using sentence transformer (e.g., `all-mpnet-base-v2`)
- **Phenotype search**: Compare against phenotype branch embeddings → top-K phenotype candidates
- **Severity search**: Compare against severity branch embeddings → top-K severity candidates
- Cosine similarity used for ranking

**Step 2: Graph Context Expansion**
- For each candidate from BOTH branches:
  - Term definition
  - Parent terms (ontology hierarchy)
  - Synonyms and relationships
- Build rich context strings for phenotypes AND severities separately

**Step 3: LLM Multi-Output Selection**
- **Negation detection**: Pre-process term for keywords ("no", "not", "without", "absent", "normal")
- LLM receives:
  - Original clinical term
  - Top-K phenotype candidates with graph context
  - Top-K severity candidates with graph context
  - Pre-detected negation flag
  - Patient context (genes, etc.)
- **LLM returns**:
  ```json
  {
    "phenotypes": [
      {"id": "HP:...", "label": "...", "excluded": true/false},
      {"id": "HP:...", "label": "...", "excluded": false}
    ],
    "severity": {"id": "HP:0012828", "label": "Severe"} or null,
    "diseases": [...]
  }
  ```
- **Multi-phenotype support**: A single term like "feeding difficulties and seizures" can return MULTIPLE phenotypes

**Step 4: Severity Label Extraction**
- Severity ID is returned by LLM (e.g., `HP:0012828`)
- Severity LABEL is looked up from the ontology (NEVER from LLM to prevent hallucination)
- This ensures severity labels are always accurate

#### Handling Special Cases

**Negation Examples:**
- "no seizures" → `{"id": "HP:0001250", "label": "Seizure", "excluded": true}`
- "not diabetic" → `{"id": "HP:0000819", "label": "Diabetes mellitus", "excluded": true}`
- "normal vision" → `{"id": "HP:0000504", "label": "Abnormality of vision", "excluded": true}`

**Severity Examples:**
- "severe intellectual disability" → 
  - Phenotype: `HP:0001249` (Intellectual disability)
  - Severity: `HP:0012828` (Severe)
- "mild seizures" → 
  - Phenotype: `HP:0001250` (Seizure)
  - Severity: `HP:0012825` (Mild)

**Multi-Phenotype Examples:**
- "feeding difficulties and seizures" →
  ```json
  {
    "phenotypes": [
      {"id": "HP:0011968", "label": "Feeding difficulties", "excluded": false},
      {"id": "HP:0001250", "label": "Seizure", "excluded": false}
    ]
  }
  ```
- "intolerance/fatigue" (alternative phrasing) →
  ```json
  {
    "phenotypes": [
      {"id": "HP:0003538", "label": "Exercise intolerance", "excluded": false}
    ]
  }
  ```

**Advantages:**
- Better semantic understanding (not just string matching)
- Branch-aware retrieval (phenotypes vs severity kept separate)
- Ontology-aware (considers parent-child relationships)
- Handles multiple phenotypes per term
- Handles synonyms and paraphrases effectively
- Severity labels always from ontology (no LLM hallucination)
- Configurable models and parameters

### 4.4. Strict OMIM Extraction
To prevent false positives (e.g., partial matches of HPO IDs like `HP:0001234` being read as `OMIM:001234`), a strict regex is used:
`(?<!HP:)(?<!HP)\b\d{6}\b`
This ensures only 6-digit numbers *not* preceded by "HP" are considered. Candidates are then validated against the local OMIM database (`mim2gene.txt`) to ensure they are real identifiers.

## 5. Key Features

### 5.1. Incremental Processing
Both pipelines support incremental processing:
- Progress saved line-by-line to output file
- Automatically resumes if interrupted
- Skips already processed rows based on ID
- Use `--overwrite` to reprocess from scratch

### 5.2. Hallucination Prevention (Standard Method)
- **Two-step validation**: Extraction followed by strict filtering
- **Text splitting**: Only phenotypes sent to LLM (not diagnosis)
- **Source validation**: Each term verified against original text
- **Result**: Zero hallucinations (verified in testing)

### 5.3. NAN/Pandas Artifact Handling
- Special handling for pandas `nan` values
- `"NAN"` gene symbol (valid synonym for SCN11A) filtered when from pandas artifacts
- Prevents false gene mappings from empty cells

### 5.4. ACMG Evidence Filtering
- Only valid ACMG/AMP 2015 codes accepted (PVS1, PS1-4, PM1-6, PP1-5, BA1, BS1-4, BP1-7)
- Non-standard codes (e.g., "KSM", "HGMD") filtered out
- Evidence extracted from variant classification strings (e.g., "LP (PM2, PP1)")

### 5.5. Solved Status
- New `solved` column added to all outputs
- fawzan-variants: Based on Result field (Positive=true, Negative=false, Ambiguous=false)
- All other files: Defaults to true when variant present
- Infers variant classification when solved=true and variant exists

### 5.6. Retry Logic
- Up to 2 retries (3 attempts total) for invalid JSON responses
- Handles API errors, empty responses, and malformed JSON
- Automatic backoff and retry

## 6. Output Formatting

The normalized data is saved as TSV files:
- **Standard method**: `data/processed/` directory
- **GraphRAG method**: `data/processed_graphrag/` directory

### Output Schema (36 columns)

*   **`id`**: Case/patient identifier
*   **`source_file`**: Original filename
*   **`sex_id` / `sex_label`**: NCIT sex terms
*   **`age`**: ISO 8601 duration format
*   **`hpo_ids` / `hpo_labels`**: Pipe-separated list of present HPO terms
*   **`hpo_severity_ids` / `hpo_severity_labels`**: **Strictly paired** list corresponding to `hpo_ids`. If a term has no severity, the entry is empty (e.g., `Sev1||Sev3`)
*   **`hpo_excluded_ids` / `hpo_excluded_labels`**: Separate list for negated terms
*   **`disease_omim_gene_ids` / `disease_omim_phenotype_ids`**: Validated OMIM IDs (genes vs phenotypes)
*   **`disease_omim_labels`**: Disease names
*   **`disease_mondo_ids` / `disease_mondo_labels`**: MONDO disease terms
*   **`gene_hgnc_ids` / `gene_symbols`**: Gene identifiers
*   **`variant_hgvs`**: HGVS notation for variants
*   **`variant_validated`**: Genomic validation status
*   **`variant_classification`**: Normalized ACMG classification (Pathogenic, Likely pathogenic, Uncertain significance, Likely benign, Benign)
*   **`variant_evidence`**: ACMG evidence codes (filtered to valid codes only)
*   **`publication_id`**: PMID extracted from filename or column
*   **`genotype_id` / `genotype_label`**: GENO zygosity terms
*   **`solved`**: Boolean (true/false) - case solved status

## 7. Validation & Testing

A test suite (`scripts/test_phenotypes_processing.py`) verifies:
*   **Structure:** All expected columns are present.
*   **Content:** IDs are valid, OMIM IDs are real and correctly formatted, Variant Classifications are standardized.
*   **Specific Cases:** Known problematic cases (e.g., "ahmed-variants" Case 1) are checked for regression.

## 8. Standalone Phenotype Matching Tool

### Overview

The phenotype matching logic has been extracted into a standalone, reusable tool located at `tools/phenotype_matcher/`. This tool provides a clean API for mapping clinical phenotype descriptions to standardized ontology identifiers without requiring the full PAVS pipeline.

### Architecture

The tool uses a **Graph RAG** (Retrieval-Augmented Generation) architecture with three main components:

1. **Semantic Retrieval**: Uses sentence transformers to encode clinical terms and find similar ontology terms via cosine similarity
2. **Graph Context Expansion**: Adds parent terms, definitions, and relationships from the ontology hierarchy
3. **LLM Selection**: Uses an LLM to select the best matches with full graph context

### Dual-Branch Strategy

The tool implements a dual-branch retrieval strategy to separate phenotypes from severity modifiers:

- **Phenotype Branch** (HP:0000118): ~14,000 clinical phenotype terms
- **Severity Branch** (HP:0012824): ~10 severity modifier terms (Mild, Moderate, Severe, etc.)

This ensures severity modifiers are always available as candidates and prevents contamination between phenotypes and severities.

### Supported Identifiers

| Ontology | Description | Example |
|----------|-------------|---------|
| **HPO** | Human Phenotype Ontology | HP:0001250 (Seizure) |
| **MONDO** | Disease Ontology | MONDO:0001234 |
| **OMIM** | Online Mendelian Inheritance in Man | OMIM:123456 (gene), OMIM:654321 (disease) |
| **OrphaNet** | Rare Disease Ontology | ORPHA:1234 (via MONDO cross-references) |

### Input/Output Schemas

**Input:**
```python
PhenotypeInput(
    text="severe intellectual disability and seizures",
    context="patient with epilepsy",  # Optional
    split_by=","  # Character to split on for multiple phenotypes
)
```

**Output:**
```python
PhenotypeOutput(
    phenotypes=[
        PhenotypeMatch(
            hpo_id="HP:0001249",
            label="Intellectual disability",
            excluded=False,
            severity_id="HP:0012828",
            severity_label="Severe",
            confidence=0.95
        ),
        PhenotypeMatch(
            hpo_id="HP:0001250",
            label="Seizure",
            excluded=False,
            severity_id=None,
            severity_label=None,
            confidence=0.92
        )
    ],
    diseases=[...],
    processing_metadata={...}
)
```

### Key Features

1. **Multi-Phenotype Extraction**: A single description can map to multiple phenotypes
   - Example: "seizures, hypotonia, and feeding difficulties" → 3 phenotypes

2. **Negation Detection**: Identifies excluded phenotypes
   - Example: "no cardiac abnormalities" → HPO:0001627 (excluded=True)

3. **Severity Modifiers**: Extracts severity as separate HPO terms
   - Example: "severe intellectual disability" → HP:0001249 + severity HP:0012828

4. **Context-Aware Disambiguation**: Uses context to resolve ambiguous terms
   - Example: "ASD" with context "cardiac defect" → Atrial Septal Defect (not Autism)

5. **Hallucination Prevention**: LLM returns labels only (not IDs), which are then looked up in the ontology to prevent invalid identifiers

### Configuration Options

**Embedding Models:**
- `fast`: all-MiniLM-L6-v2 (80MB, fastest)
- `balanced`: all-mpnet-base-v2 (420MB, recommended)
- `accurate`: BAAI/bge-large-en-v1.5 (1.3GB, best accuracy)
- `medical`: PubMedBERT (medical domain-specific)
- `biobert`: BioBERT (biomedical literature)

**LLM Models:**
- `fast`: openai/gpt-oss-120b (recommended for production)
- `accurate`: anthropic/claude-3.5-sonnet (best accuracy)
- `balanced`: google/gemini-2.0-flash-exp:free (testing)
- `cheap`: deepseek/deepseek-chat (large-scale processing)

### Installation & Usage

**Install:**
```bash
cd tools/phenotype_matcher
pip install -e .
```

**Python API:**
```python
from phenotype_matcher import PhenotypeMatcher, PhenotypeInput, MatcherConfig

# Basic usage
matcher = PhenotypeMatcher()
output = matcher.match(PhenotypeInput(text="severe intellectual disability"))

# Access results
hpo_ids = output.get_hpo_ids()  # Present phenotypes
excluded_ids = output.get_hpo_ids(excluded=True)  # Negated
omim_ids = output.get_omim_ids()
mondo_ids = output.get_mondo_ids()
orphanet_ids = output.get_orphanet_ids()

# Custom configuration
config = MatcherConfig(
    embedding_model="medical",
    llm_model="accurate",
    top_k_phenotype=10,
    device="cuda"
)
matcher = PhenotypeMatcher(config)
```

**Command Line:**
```bash
# Single phenotype
phenotype-matcher "severe intellectual disability"

# Multiple phenotypes
phenotype-matcher "seizures, hypotonia, feeding difficulties"

# Different output formats
phenotype-matcher "seizures" --output-format hpo_ids
phenotype-matcher "seizures" --output-format tsv
phenotype-matcher "seizures" --output-format summary

# Batch processing
phenotype-matcher --input phenotypes.txt --output results.json

# Model selection
phenotype-matcher "developmental delay" \
  --embedding-model medical \
  --llm-model accurate \
  --top-k 10
```

### Performance

- **First run**: 5-10 minutes (computing embeddings for ~35,000 ontology terms)
- **Subsequent runs**: ~1 second initialization (loading from cache)
- **Per-term matching**: ~2-5 seconds (depending on LLM model)
- **Cache size**: ~100-500MB (depending on embedding model)

### Documentation

Complete documentation is available at:
- `tools/phenotype_matcher/README.md` - Full user guide
- `tools/phenotype_matcher/schemas.py` - API documentation
- `docs/GRAPHRAG_DESIGN.md` - Architecture details

### Integration with PAVS Pipeline

The normalization scripts (`scripts/process_all_phenotypes_graphrag.py`) use the same Graph RAG logic as this standalone tool. The tool can be used independently or as part of the PAVS pipeline.
