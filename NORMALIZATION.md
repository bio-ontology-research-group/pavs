# PAVS: Data Normalization Documentation

This document describes the pipeline for standardizing clinical phenotype and variant datasets (Ahmed, Fawzan, and Marwa) into normalized formats.

## Pipeline Overview

All normalization scripts are located in the `scripts/` directory and utilize a combination of exact ontology matching, Graph-RAG retrieval, and LLM-based selection (Gemini 2.0 Flash) to ensure high-quality clinical mappings.

### Core Technologies
- **Graph-RAG**: Uses the HPO hierarchy (`networkx`) to provide context (parents/children) to the LLM during mapping.
- **Rapidfuzz**: Performs fuzzy string matching for initial candidate retrieval.
- **OpenRouter (Gemini 2.0 Flash)**: Resolves complex phenotype splits and performs semantic selection between ontology candidates.
- **Ontologies**:
    - **HPO**: `ontology/hp.obo` (Phenotypes & Mode of Inheritance)
    - **GENO**: Genotype Ontology (Zygosity)
    - **OMIM/ORPHA**: `data/phenotype.hpoa` (Diseases)
    - **Entrez**: `data/Homo_sapiens.gene_info` (Genes)

## Normalized Datasets

The outputs are stored in `data/cleaned_and_normalized/`. All files preserve their original columns for traceability while adding a set of `normalized_` or `norm_` fields.

### 1. Ahmed Normalization (`ahmed_normalized.tsv`)
**Focus**: Complex phenotypic sentences and mode of inheritance mapping.
- **Standardized Fields**:
    - `normalized_hpos`: List of validated HPO IDs.
    - `normalized_hpo_labels`: Official HPO names.
    - `normalized_gene_entrez`: NCBI Entrez Gene ID.
    - `normalized_zygosity_geno`: GENO ID.
    - `normalized_moi_hpo`: HPO ID for Mode of Inheritance.
    - `normalized_omim_data`: JSON string with disease ID, label, and linked gene.

### 2. Fawzan Normalization (`fawzan_normalized.tsv`)
**Focus**: Variant nomenclature (HGVS) and HGMD field parsing.
- **Standardized Fields**:
    - `normalized_variants_hgvs`: Normalized HGVS strings (e.g., `NM_001194998(CEP152):c.2021G>T`).
    - `normalized_genes_entrez`: Mapped Entrez IDs.
    - `normalized_zygosity_geno`: GENO IDs (cleaned of variant debris).
    - `normalized_hgmd_diseases`: Disease IDs mapped from the HGMD info string via RAG.
    - `normalized_dbsnp`: Extracted rsIDs.

### 3. Marwa Normalization (`marwa_normalized.tsv`)
**Focus**: Disentangling mixed phenotype and disease descriptions.
- **Standardized Fields**:
    - `normalized_hpos`: Phenotypic features.
    - `normalized_diseases`: Syndrome/Disease identifiers (OMIM/ORPHA).
    - `normalized_disease_labels`: Official disease names from the annotation store.

### 4. Ahmed PMID Normalization (`ahmed_pmid28454995_normalized.tsv`)
**Focus**: Merging multiple clinical descriptions and mapping variants from separate columns.
- **Standardized Fields**:
    - `normalized_hpos`: Merged features from `Diagnosis` and `Additional clinical phenotype`.
    - `normalized_variants_hgvs`: Constructed HGVS strings from transcript and cDNA columns.
    - `normalized_genes_entrez`: Mapped variant genes to Entrez IDs.
    - `normalized_diseases`: Parsed OMIM identifiers.

### 5. PMC7082194 Normalization (`pmc7082194_normalized.tsv`)
**Focus**: Extracting phenotypes from non-standard trailing columns and standardizing inheritance to GENO.
- **Standardized Fields**:
    - `normalized_hpos`: Phenotypic features extracted from the final data column.
    - `normalized_variants_hgvs`: Standardized variant nomenclature.
    - `normalized_zygosity_geno`: Inheritance terms (het/homo) mapped to GENO.

## How to Run the Normalization

Requires `OPENROUTER_API_KEY` to be set in the environment.

```bash
# Normalize Ahmed dataset
uv run --with pandas --with rapidfuzz --with networkx python scripts/normalize_ahmed.py

# Normalize Fawzan dataset
uv run --with pandas --with rapidfuzz --with networkx python scripts/normalize_fawzan.py

# Normalize Marwa dataset
uv run --with pandas --with rapidfuzz --with networkx python scripts/normalize_marwa.py

# Normalize Ahmed PMID dataset
uv run --with pandas --with rapidfuzz --with networkx python scripts/normalize_ahmed_pmid28454995.py

# Normalize PMC7082194 dataset
uv run --with pandas --with rapidfuzz --with networkx python scripts/normalize_pmc7082194.py
```

## Standalone Phenotype Matching Tool

A standalone, reusable phenotype matching tool is now available at `tools/phenotype_matcher/`. This tool can be used independently to map clinical phenotype descriptions to standardized ontology identifiers.

### Features
- **Input**: Free-text phenotype descriptions (single or multiple phenotypes)
- **Output**: HPO, MONDO, OMIM, and OrphaNet identifiers with negation and severity modifiers
- **Method**: Graph RAG (semantic embeddings + graph structure + LLM reasoning)

### Installation

```bash
cd tools/phenotype_matcher
pip install -e .
```

### Basic Usage

**Python API:**
```python
from phenotype_matcher import PhenotypeMatcher, PhenotypeInput

matcher = PhenotypeMatcher()
input_data = PhenotypeInput(text="severe intellectual disability and seizures")
output = matcher.match(input_data)

for pheno in output.phenotypes:
    print(f"{pheno.hpo_id}: {pheno.label}")
    if pheno.severity_label:
        print(f"  Severity: {pheno.severity_label}")
```

**Command Line:**
```bash
phenotype-matcher "severe intellectual disability"
phenotype-matcher "seizures, hypotonia" --output-format tsv
```

See `tools/phenotype_matcher/README.md` for comprehensive documentation.

## Verification

The normalization quality is monitored via `scripts/verify_standardization.py`, which checks for:
1.  **Semantic Count Mismatches**: Compares extracted feature counts with original text structure.
2.  **ID Validity**: Ensures all HPO and GENO IDs exist in their respective ontologies.
3.  **HGVS Syntax**: Validates nomenclature prefixes.
