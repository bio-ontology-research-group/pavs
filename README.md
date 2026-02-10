# PAVS: Phenotype-Associated Variants in Saudi Populations

A pipeline to integrate clinical genetic data from disparate sources into a unified, standardized dataset and generate [GA4GH Phenopackets v2.0.2](https://phenopacket-schema.readthedocs.io/en/latest/).

## Overview

This project unifies data from multiple Excel/TSV sources, normalizes clinical terms to Human Phenotype Ontology (HPO), standardizes variants to HGVS-like syntax, and exports them as Phenopackets for use in Middle Eastern and Saudi Arabian population studies.

## Data Sources

The master dataset is integrated from the following files in `data/`:
- `439_2017_1821_MOESM1_ESM.tsv`: Literature-mined data.
- `Manually curated Data.tsv`: Clinical data from collaborators.
- `Variant list.tsv`: Additional curated variant records.

## Pipeline Workflow

The pipeline consists of three main stages, all executed using `uv` for dependency management:

### 1. Data Processing & Integration
Cleans raw source files, handles inconsistent column naming (e.g., "Zyogsity"), extracts basic demographics (Sex, Age), and maps zygosity to [Genotype Ontology (GENO)](http://purl.obolibrary.org/obo/geno.owl) terms.
```bash
uv run scripts/process_data.py
```
**Output:** `master_data_v1.tsv`

### 2. HPO Normalization (LLM-Assisted)
Extracts unique clinical terms and maps them to HPO IDs using a local OBO lookup followed by an LLM-based disambiguation (Gemini 2.0 Flash) for ambiguous or complex terms.
```bash
# Extract and map terms
uv run scripts/map_text_to_hpo.py
# (Optional) Run LLM batch processor if new terms are added
uv run scripts/hpo_llm_map.py
```
**Final Step:** Reconcile all mappings and clean variant strings:
```bash
uv run scripts/finalize_master.py
```
**Output:** `master_data_final.tsv` and `hpo_mappings_final.json`

### 3. Phenopacket Generation
Generates Phenopacket v2.0.2 JSON files.
- **Normalization:** Variants are normalized to `Transcript:c.Change` or `Gene:c.Change`.
- **Metadata:** Includes HANCESTRO:0014 (Middle Eastern) and GAZ:00000570 (Saudi Arabia) tags.
- **Taxonomy:** Standardized to NCBITaxon:9606 (Homo sapiens).

```bash
uv run scripts/generate_phenopackets.py
```
**Output:** 2,692 JSON files in `phenopackets/generated/` and a `phenopacket_store.summary.tsv` manifest.

## How to Regenerate Everything

To perform a clean run of the entire pipeline:

1. **Ensure Resources are present:**
   - `ontology/hp.obo` (Human Phenotype Ontology)

2. **Execute Pipeline:**
   ```bash
   # Unify and clean raw data
   uv run scripts/process_data.py
   
   # Reconcile HPO mappings and finalize master TSV
   uv run scripts/finalize_master.py
   
   # Generate Phenopackets
   uv run scripts/generate_phenopackets.py
   ```

## Standardized Formats

- **Phenotypes:** Human Phenotype Ontology (HP:XXXXXXX)
- **Zygosity:** Genotype Ontology (GENO:XXXXXXX)
- **Variants:** HGVS-like (`Transcript:c.XXX`)
- **Population:** Middle Eastern (HANCESTRO:0014), Saudi Arabia (GAZ:00000570)
- **Reference Build:** GRCh37 (default for this dataset)

## Metadata
- **Creator:** Robert Hoehndorf
- **ORCID:** https://orcid.org/0000-0001-8149-5890
- **Schema Version:** 2.0.2