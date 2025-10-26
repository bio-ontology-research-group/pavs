# PAVS Database Pipeline

Phenotype-Associated Variants in Saudi Populations (PAVS) - A standardized pipeline for parsing and structuring clinical data from a CSV file into standardized [Phenopackets](https://phenopacket-schema.readthedocs.io/en/latest/) format with automated HPO annotation.


## Quick Start

### Installation
```bash
git clone https://github.com/your-username/pavs.git
cd pavs
pip install -r requirements.txt
```


### Resources Setup

The `resources/` directory contains:
- **`hp.obo`** - HPO ontology file  
  Download: http://purl.obolibrary.org/obo/hp.obo
- **`en_core_sci_sm-0.5.4/`** - spaCy biomedical NLP model  
  Download: https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_core_sci_sm-0.5.4.tar.gz

---

## Pipeline Scripts

### 1. `process_and_merge_pavs.py` - Data Integration

Merges and standardizes data from multiple sources.

**Usage:**
```bash
python scripts/process_and_merge_pavs.py
```

**What it does:**
- Integrates 3 data sources (literature, clinical data, hospital collaborators)
- Standardizes column names and formats
- Extracts HPO IDs from text using regex
- Removes duplicates
- Generates unified PAVS IDs (PAVS0001-PAVS2579)

**Input:**
- Multiple source files in `data/` directory

**Output:**
- `data/PAVS_final_data.tsv` - Merged and standardized dataset

---

### 2. `map_text_to_hpo.py` - HPO Annotation

Annotates free-text phenotypes with HPO terms using NLP.

**Usage:**
```bash
python scripts/map_text_to_hpo.py \
    -i data/phenotype_text.tsv \
    -o data/hpo_annotations.json \
    --hpo resources/hp.obo
```

**Parameters:**
- `-i, --input` - Input TSV file (required)
- `-o, --output` - Output JSON file (required)
- `--hpo` - Path to HPO .obo file (default: `resources/hp.obo`)


**What it does:**
- Uses spaCy biomedical model (`en_core_sci_sm`) for entity recognition
- Detects and excludes negated findings using NegSpaCy
- Fuzzy matches clinical terms to HPO 
- Combines with regex-extracted HPO IDs
- Deduplicates HPO terms

**Input Format:**

TSV file with columns:
- `case_id` - Patient/case identifier
- `phenotype_text` - Free-text phenotype description

**Output Format:**
```json
{
  "CASE001": ["HP:0001250", "HP:0001263"],
  "CASE002": ["HP:0002104", "HP:0000252"]
}
```
**Requirements:**
- `resources/hp.obo` - HPO ontology
- `resources/en_core_sci_sm-0.5.4/` - spaCy model

---

### 3. `convert_to_phenopackets.py` - Format Conversion

Converts processed data to Phenopackets v2 format.

**Usage:**
```bash
python scripts/convert_to_phenopackets.py \
    -i data/PAVS_final_data.tsv \
    -o data/PAVS_phenopackets.json \
    --hpo resources/hp.obo
```

**Parameters:**
- `-i, --input` - Input TSV file (default: `PAVS_final_data.tsv`)
- `-o, --output` - Output file (default: `PAVS_phenopackets.json`)
- `-f, --format` - Output format: `json`, `json-lines`, or `individual` (default: `json`)
- `--hpo` - Path to HPO .obo file (default: `resources/hp.obo`)

**Output Formats:**
```bash
# Single JSON file (default)
python scripts/convert_to_phenopackets.py -o PAVS_phenopackets.json

# JSON Lines (one phenopacket per line)
python scripts/convert_to_phenopackets.py -f json-lines -o PAVS_phenopackets.json

# Individual files (one per patient)
python scripts/convert_to_phenopackets.py -f individual -o PAVS_phenopackets/
```

**What it does:**
- Converts to Phenopackets v2 schema
- Populates HPO term labels from ontology (99.9% coverage)
- Extracts gene symbols from variant descriptions (69.8% coverage)
- Extracts OMIM IDs from diagnosis, phenotypes, and comments
- Parses age to ISO8601 duration format (P22Y, P6M, etc.)
- Maps ACMG pathogenicity classifications
- Adds comprehensive metadata

**Output Statistics:**
```
Total phenopackets: 2,579
Phenotypic features: 11,130 HPO terms (99.9% with labels)
Genomic variants: 2,527 (69.8% with gene symbols)
Clinical diagnoses: 1,985 (77.0% coverage)
Age information: 543 (21.1%)
Sex specified: 1,850 (71.7%)
```

---

### 4. `validate_phenopackets.py` - Validation

Validates Phenopackets schema compliance and data quality.

**Usage:**
```bash
python scripts/validate_phenopackets.py data/PAVS_phenopackets.json
```

**Parameters:**
- `file` - Phenopackets JSON file to validate (required)
- `-v, --verbose` - Show detailed warnings (optional)


---

## Data Files

### Input: `data/PAVS_final_data.tsv`

Processed dataset with standardized columns:

| Column | Description | Example |
|--------|-------------|---------|
| `ID` | Unique identifier | PAVS0001 |
| `sex` | Patient sex | male, female |
| `age` | Age at diagnosis | 22Y, 6M, infant |
| `phenotypicFeatureIds` | Comma-separated HPO IDs | HP:0001250,HP:0001263 |
| `genomicVariants` | Semicolon-separated variants | SCN1A:c.2398G>A;BRCA1:p.Arg1751Ter |
| `diagnosis` | Clinical diagnosis | Dravet syndrome |
| `consanguinityStatus` | Consanguinity status | yes, no, unknown |
| `familyId` | Family identifier | FAM001 |
| `procedure` | Test type | WES, WGS, Panel |
| `procedureStrategy` | Testing strategy | Trio, Singleton |
| `diagnosticComment` | Additional clinical notes | Free text |
| `zygosityStatus` | Variant zygosity | homozygous, heterozygous |
| `variantInterpretation` | ACMG classification | Pathogenic, Likely pathogenic, VUS |
| `dataSourceType` | Data source | literature, clinical_collaborator |
| `externalReference` | Publication reference | PMID:12345678, DOI:10.1234/... |



## Citation



