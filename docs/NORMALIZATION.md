# PAVS: Normalization Pipeline

This document describes the pipeline for standardizing clinical phenotype and variant datasets from the 7 source TSV files in `data/phenotypes/` into normalized formats.

## Pipeline Overview

All normalization scripts are in the `normalization/` directory. They use `phenotype_matcher_v2` (the current production tool in `tools/phenotype_matcher_v2/`) for HPO/MONDO/OMIM/Orphanet mapping.

The pipeline is coordinated by `normalization/combine_normalize_phenotypes.py`, which calls source-specific normalizers and writes the unified output.

### Core Technologies

- **`phenotype_matcher_v2`**: HPO/MONDO/OMIM/Orphanet mapping via Aho-Corasick exact match, stemmed fuzzy match, SapBERT ANN semantic search, and LLM disambiguation. Negation and severity-modifier detection are integrated at the span level. See `tools/phenotype_matcher_v2/ALGORITHM.md` for the full algorithm. See `docs/NORMALIZATION_METHODS.md` for source-specific parsing details including HGMD variant handling for fawzan-variants.
- **Rapidfuzz**: Fuzzy string matching for candidate retrieval.
- **OpenRouter (Gemini/DeepSeek)**: Resolves complex phenotype splits and performs semantic selection between ontology candidates.
- **Ontologies**:
    - **HPO**: `ontology/hp.obo` (Phenotypes & Mode of Inheritance)
    - **MONDO**: `ontology/mondo.obo` (Disease Ontology)
    - **GENO**: Genotype Ontology (Zygosity)
    - **NCIT**: NCI Thesaurus (Sex)
    - **OMIM/ORPHA**: from `data/phenotype.hpoa` (Diseases)
    - **Entrez**: `data/Homo_sapiens.gene_info` (Genes)

## How to Run the Normalization

Requires `OPENROUTER_API_KEY` to be set in the environment.

```bash
uv sync
export OPENROUTER_API_KEY="your_key_here"
uv pip install -e tools/phenotype_matcher_v2

# Normalize all 7 source files → combined_normalized.tsv
uv run python normalization/combine_normalize_phenotypes.py \
    --workers 4 \
    --output data/combined_normalized.tsv \
    --acmg-obo ontology/acmg_criteria.obo

# Smoke test (5 rows only):
uv run python normalization/combine_normalize_phenotypes.py \
    --limit 5 --workers 1 \
    --output /tmp/test_combined.tsv \
    --acmg-obo /tmp/test_acmg.obo
```

To run individual source normalizers directly (for testing or partial re-processing):

```bash
# Ahmed dataset
uv run python normalization/normalize_ahmed.py

# Fawzan dataset
uv run python normalization/normalize_fawzan.py

# Marwa dataset
uv run python normalization/normalize_marwa.py

# Ahmed PMID dataset
uv run python normalization/normalize_ahmed_pmid28454995.py

# PMC7082194 dataset
uv run python normalization/normalize_pmc7082194.py

# PMC6562004 dataset
uv run python normalization/normalize_pmc6562004.py

# DDD dataset
uv run python normalization/normalize_ddd.py
```

## Source-Specific Normalizers

### 1. Ahmed Normalization (`normalize_ahmed*.py`)

**Source**: `data/phenotypes/ahmed*` (A: 291 rows, B: 234 rows)
**Focus**: Complex phenotypic sentences, mode of inheritance mapping, and merging multiple clinical descriptions.

Standardized fields:
- `normalized_hpos`: Features from sentences, `Diagnosis` and `Additional clinical phenotype` columns.
- `normalized_variants_hgvs`: Constructed from transcript/cDNA columns.
- `normalized_gene_entrez`: NCBI Entrez Gene ID.
- `normalized_zygosity_geno`: GENO ID.
- `normalized_moi_hpo`: HPO ID for Mode of Inheritance.
- `normalized_omim_data`: JSON string with disease ID, label, and linked gene.

### 2. Fawzan Normalization (`normalize_fawzan.py`)

**Source**: `data/phenotypes/fawzan-variants.tsv` (1,024 rows, source letter F)
**Focus**: Variant nomenclature (HGVS) and HGMD field parsing.

Standardized fields:
- `normalized_variants_hgvs`: Normalized HGVS strings (e.g., `NM_001194998(CEP152):c.2021G>T`).
- `normalized_genes_entrez`: Mapped Entrez IDs.
- `normalized_zygosity_geno`: GENO IDs (cleaned of variant debris).
- `normalized_hgmd_diseases`: Disease IDs mapped from the HGMD info string.
- `normalized_dbsnp`: Extracted rsIDs.

HGMD `Info=` strings are parsed for `ID=` (CM accession), `Ref=`/`Alt=` alleles, `dbSNP=` rsID, and `Disease=`. GRCh38 coordinates are obtained via the Ensembl variation API using the rsID. Annotation uses the VEP region endpoint (`/vep/human/region/{chrom}:{pos}-{pos}:1/{alt}`). See `docs/NORMALIZATION_METHODS.md` for details.

### 3. Marwa Normalization (`normalize_marwa.py`)

**Source**: `data/phenotypes/marwa-variants.tsv` (1,421 rows, source letter M)
**Focus**: Disentangling mixed phenotype and disease descriptions.

Standardized fields:
- `normalized_hpos`: Phenotypic features.
- `normalized_diseases`: Syndrome/Disease identifiers (OMIM/ORPHA).
- `normalized_disease_labels`: Official disease names from the annotation store.

### 3. PMC7082194 Normalization (`normalize_pmc7082194.py`)

**Source**: `data/phenotypes/PMC7082194.tsv` (522 rows, source letter Q)
**Focus**: Extracting phenotypes from non-standard trailing columns and standardizing inheritance to GENO.

Standardized fields:
- `normalized_hpos`: Phenotypic features extracted from the `Unnamed:13` column.
- `normalized_variants_hgvs`: Standardized variant nomenclature (converting 1-letter to 3-letter HGVS p.).
- `normalized_zygosity_geno`: Inheritance terms (het/homo) mapped to GENO.

### 4. PMC6562004 Normalization (`normalize_pmc6562004.py`)

**Source**: `data/phenotypes/PMC6562004.tsv` (2,218 rows, source letter P)
**Focus**: Multi-gene/variant columns, skips first comment row.

### 5. DDD Normalization (`normalize_ddd.py`)

**Source**: `data/phenotypes/ddd-diagnoses.tsv` (1,856 rows, source letter D)
**Focus**: Gene-disease associations with semicolon-delimited HPO IDs and allelic mode.

## Shared Utilities

`normalization/normalization_utils.py` — singleton `NormalizationUtils` class:

```python
from normalization.normalization_utils import NormalizationUtils
nu = NormalizationUtils()

# Normalize sex to NCIT
sex_id, sex_label = nu.normalize_sex("Male")  # → NCIT:C20197, Male

# Normalize zygosity to GENO
geno_id, geno_label = nu.normalize_genotype("homozygous")  # → GENO:0000136, homozygous

# Normalize consanguinity
cons_id, cons_label = nu.normalize_consanguinity("yes")  # → HP:0001006, yes
```

The class loads HPO, MONDO, GENO, NCIT, OMIM, and Entrez gene info at initialization.

## Output

**`data/combined_normalized.tsv`** — 38 columns, PAVS:XNNNNNNN IDs, pipe-delimited multi-value fields.

See `docs/PIPELINE.md` for the full 38-column schema.

## Verification

Normalization quality is checked via `analysis/verify_standardization.py`, which verifies:

1. **ID Validity**: Ensures all HPO and GENO IDs exist in their respective ontologies.
2. **HGVS Syntax**: Validates nomenclature prefixes.
3. **Semantic Count**: Compares extracted feature counts with original text structure.
