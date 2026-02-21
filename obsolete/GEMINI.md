# PAVS: Phenotypic and Variant Standardization Pipeline

## Project Overview
This project is a clinical data processing pipeline designed to ingest heterogeneous genomic and phenotypic datasets (originally in Excel/CSV), standardize them to international ontologies, and export them as **Phenopackets v2.0**. It specifically focuses on cases from the Middle East/Saudi Arabia.

### Key Technologies
- **Python**: Core logic for data parsing and transformation.
- **Pandas**: Data manipulation and integrated master sheet generation.
- **Ontologies**: HPO (Phenotypes), OMIM (Diseases), HANCESTRO (Ancestry), GAZ (Geography), GENO (Zygosity).
- **LLM Integration**: Uses Gemini 1.5/2.0 Flash via OpenRouter for resolving ambiguous clinical text to HPO terms.
- **uv**: Fast Python package manager for executing scripts with specific dependencies.

### Architecture
1.  **Ingestion**: Converts Excel sheets in `data/` to TSV.
2.  **Normalization**:
    *   `process_data.py`: Initial extraction and ID generation.
    *   `hpo_llm_map.py`: RAG-based LLM mapping for ambiguous phenotypes.
    *   `finalize_master.py`: Integration of all mappings into `master_data_final.tsv`.
3.  **Export**: `generate_phenopackets.py` creates Phenopacket v2.0 JSON files.

## Pipeline Execution
All scripts should be run using `uv` to ensure consistent dependency management.

### 1. Initial Processing
Extracts obvious IDs and prepares an ambiguity queue for LLM mapping.
```bash
uv run --with pandas python scripts/process_data.py
```

### 2. HPO Normalization (LLM)
Maps free-text clinical descriptions to HPO IDs. Requires `OPENROUTER_API_KEY`.
```bash
export OPENROUTER_API_KEY="your_key_here"
uv run --with rapidfuzz --with requests --with pandas python scripts/hpo_llm_map.py
```

### 3. Finalization
Integrates LLM results and cleans variant nomenclature.
```bash
uv run --with pandas python scripts/finalize_master.py
```

### 4. Phenopacket Generation
Generates Phenopacket v2.0 JSON files and a store summary.
```bash
uv run --with pandas python scripts/generate_phenopackets.py
```

## Development Conventions
- **Provenance**: Always maintain `Source` and `Original_ID` in derived datasets.
- **IDs**: Master IDs follow the pattern `PAVS_MASTER_NNNN`.
- **Genomics**: Variants are standardized to HGVS strings; default genome build is **GRCh37**.
- **Population Context**: All phenopackets are tagged with `HANCESTRO:0014` (Middle Eastern) and `GAZ:00000570` (Saudi Arabia).
- **Metadata**: Generated artifacts must credit "Robert Hoehndorf" (ORCID: 0000-0001-8149-5890).

## Directory Structure
- `data/`: Source TSV files and original Excel backups (`data/excel_files`).
- `ontology/`: HPO OBO/OWL files.
- `phenopackets/`: Publicly available store (`0.1.26`) and generated packets (`generated/`).
- `scripts/`: Python pipeline scripts.
