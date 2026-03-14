# PAVS Knowledge Graph

Data processing pipeline and RDF generation for the PAVS (Phenotypic and Variant Standardization) knowledge graph of Saudi Arabian rare disease patients.

## Overview

This repository contains:
- **Source data**: 7 TSV files with clinical phenotypes and genomic variants from Saudi rare disease cohorts
- **Normalization pipeline**: Scripts to standardize phenotypes to HPO, diseases to MONDO/OMIM, variants to HGVS
- **Annotation pipeline**: Enrichment with VEP, ClinVar, gnomAD, ClinGen, mouse phenotypes, GO, GTEx
- **Phenopacket generation**: GA4GH Phenopackets v2.0.2 output
- **RDF generation**: Turtle (.ttl) files for Virtuoso triple store

## Pipeline Overview

```
data/phenotypes/          (7 source TSV files)
      │
      ▼  combine_normalize_phenotypes.py
      │
combined_normalized.tsv   (38 columns, ~7,500 rows)
      │
      ▼  annotate_variants.py
      │
combined_annotated.tsv    (77+ columns)
      │
      ▼  generate_phenopackets_v2.py
      │
phenopackets/generated_v2/  +  PAVS_phenopackets.json
      │
      ▼  compute_hpo_ic.py  +  generate_rdf.py
      │
rdf_output/*.ttl  →  Virtuoso triple store
```

## Quick Start

### Prerequisites

```bash
# Install uv (recommended) or use pip
curl -LsSf https://astral.sh/uv/install.sh | sh

# Install phenotype-matcher library
git clone git@github.com:bio-ontology-research-group/phenotype-matcher.git
cd phenotype-matcher
pip install -e .
cd ..

# Clone this repository
git clone git@github.com:bio-ontology-research-group/pavs-knowledge-graph.git
cd pavs-knowledge-graph

# Install dependencies
uv sync
# or: pip install -r requirements.txt

# Set OpenRouter API key for normalization
export OPENROUTER_API_KEY="your_key_here"
```

### Download Reference Data

```bash
# Download all external reference databases (~5 GB)
bash scripts/download_reference_data.sh
```

This downloads to `data/reference/`:
- `phenotype.hpoa` - HPO disease annotations (~35 MB)
- `genes_to_disease.txt` - Gene→disease associations
- `gnomad.v4.1.constraint_metrics.tsv` - pLI/LOEUF scores
- `clinvar.vcf.gz` - ClinVar variants (GRCh38)
- `erepo.tabbed.txt` - ClinGen expert-curated variants
- `goa_human.gaf.gz` - GO annotations
- `Homo_sapiens.gene_info` - NCBI gene info
- `HMD_HumanPhenotype.rpt` + `MGI_GenePheno.rpt` - Mouse phenotype homologs
- `E-GTEX-8/` - GTEx v8 tissue expression
- `allele-frequencies/` - Saudi population allele frequencies

### Download Ontologies

```bash
# HPO
wget https://purl.obolibrary.org/obo/hp.obo -P ontology/

# MONDO
wget http://purl.obolibrary.org/obo/mondo.obo -P ontology/

# Additional ontologies (optional, used for normalization)
wget http://purl.obolibrary.org/obo/geno.obo -P ontology/
wget http://purl.obolibrary.org/obo/ncit.obo -P ontology/
wget http://purl.obolibrary.org/obo/mp.obo -P ontology/
wget http://purl.obolibrary.org/obo/go-basic.obo -P ontology/
```

## Running the Pipeline

### Step 1: Normalize Source Data

Combine all 7 source TSVs and normalize phenotypes to HPO, diseases to MONDO/OMIM, genotypes to GENO, sex to NCIT:

```bash
uv run python normalization/combine_normalize_phenotypes.py \
    --workers 4 \
    --output data/combined_normalized.tsv \
    --acmg-obo ontology/acmg_criteria.obo
```

**Output**: `data/combined_normalized.tsv` (38 columns, PAVS:XNNNNNNN IDs)

**Smoke test** (5 rows only):
```bash
uv run python normalization/combine_normalize_phenotypes.py \
    --limit 5 --workers 1 \
    --output /tmp/test_combined.tsv \
    --acmg-obo /tmp/test_acmg.obo
```

### Step 2: Annotate Variants and Genes

Enrich each row with variant-level and gene-level annotations:

```bash
uv run python normalization/annotate_variants.py \
    --input  data/combined_normalized.tsv \
    --output data/combined_annotated.tsv
```

**Output**: `data/combined_annotated.tsv` (~77 columns)

**Smoke test** (10 rows, skip VEP REST calls):
```bash
uv run python normalization/annotate_variants.py \
    --input  data/combined_normalized.tsv \
    --output /tmp/test_annotated.tsv \
    --limit 10 --no-vep
```

Variant annotations added:
- Ensembl VEP: consequence, impact, SIFT, PolyPhen, gnomAD AFs, rsID
- ClinVar: clinical significance, disease, allele ID
- ClinGen: assertion, ACMG codes
- Saudi population: AC/AN from 302 WGS individuals

Gene annotations added:
- Gene-disease associations from OMIM
- pLI/LOEUF constraint metrics from gnomAD
- Mouse phenotypes from MGI
- GO terms (biological process, molecular function, cellular component)
- GTEx tissue expression (percentile rank)

### Step 3: Generate Phenopackets v2

Convert to GA4GH Phenopackets v2.0.2 JSON:

```bash
uv run python intake/generate_phenopackets_v2.py \
    --input data/combined_annotated.tsv \
    --output-dir phenopackets/generated_v2 \
    --combined data/PAVS_phenopackets.json \
    --zip data/PAVS_phenopackets.zip
```

**Outputs**:
- `phenopackets/generated_v2/` - one JSON per case
- `data/PAVS_phenopackets.json` - combined array (for upload to pavs-phenopackets repo)
- `data/PAVS_phenopackets.zip` - download bundle

### Step 4: Compute HPO Information Content

Calculate IC values for semantic similarity:

```bash
uv run python intake/compute_hpo_ic.py \
    --hpo ontology/hp.obo \
    --hpoa data/reference/phenotype.hpoa \
    --output rdf_output/hpo_ic.ttl
```

**Output**: `rdf_output/hpo_ic.ttl` (IC values + HPO hierarchy + labels)

### Step 5: Generate All RDF

Convert to Turtle format for Virtuoso:

```bash
uv run python intake/generate_rdf.py \
    --phenopackets data/PAVS_phenopackets.json \
    --annotated data/combined_annotated.tsv \
    --literature-dir phenopackets/0.1.26 \
    --output-dir rdf_output/
```

**Outputs** (in `rdf_output/`):
- `cases.ttl` - Saudi case records (~192K lines)
- `genes.ttl` - Gene annotations (~103K lines)
- `hpoa.ttl` - Disease→HPO (~2.2M lines)
- `hpo_ic.ttl` - IC values + hierarchy (~56K lines)
- `literature.ttl` - Literature phenopackets (~455K lines)

### Step 6: Load into Virtuoso

See [pavs-website](https://github.com/bio-ontology-research-group/pavs-website) for Virtuoso deployment. The RDF files can be mounted as a Docker volume or loaded via `isql`:

```bash
# Manual load via isql (if Virtuoso is running)
isql virtuoso:1111 dba pavs_dba /path/to/load_ttl.sql
```

## Source Data

| File | Letter | Rows | Description |
|---|---|---|---|
| `ahmed-pmid28454995.tsv` | A, B | 525 | HPOs, Gene, Variant, Protein, Zygosity, ACMG |
| `fawzan-variants.tsv` | F | 1,024 | Phenotype, Variant(s), Zyogsity [sic] |
| `marwa-variants.tsv` | M | 1,421 | phenotypes (pipe+HPO), variants (pipe) |
| `PMC6562004.tsv` | P | 2,218 | Gene(s), Variant(s) |
| `PMC7082194.tsv` | Q | 522 | Unnamed:13=phenotype, protein 1-letter AAs |
| `ddd-diagnoses.tsv` | D | 1,856 | gene-disease assoc, hpo_ids (semicolons), allelic_mode |

**Total**: 7,566 cases (5,710 Saudi, 1,856 non-Saudi DDD controls)

### Case IDs

Format: `PAVS:XNNNNNNN` where X is the source letter (A/B/F/M/P/Q/D) and NNNNNNN is a 7-digit zero-padded sequence.

## Key Conventions

- **Genome build**: GRCh38 for VCF coordinates; GRCh37 for legacy HGVS strings
- **Population tags**: Saudi cases tagged with `HANCESTRO:0852` (Middle Eastern) and `GAZ:00005279` (Saudi Arabia)
- **Multi-value fields**: pipe-delimited (`|`) in TSV outputs
- **Provenance**: `source_file` and `source_id` columns preserved throughout

## Development

### Adding a New Source Dataset

1. Add TSV file to `data/phenotypes/`
2. Create normalizer in `normalization/normalize_NEWSOURCE.py`
3. Update `normalization/combine_normalize_phenotypes.py` to import and call the new normalizer
4. Choose a unique source letter and update case ID generation
5. Re-run pipeline from Step 1

### Updating Annotations

Reference databases are downloaded by `scripts/download_reference_data.sh`. To update:

1. Edit download URLs in the script
2. Re-run `bash scripts/download_reference_data.sh`
3. Re-run annotation pipeline from Step 2

## Docker Support

For reproducible pipeline execution:

```bash
# Build image
docker build -t pavs-pipeline -f Dockerfile .

# Run normalization
docker run --rm -v $(pwd):/app pavs-pipeline \
    python normalization/combine_normalize_phenotypes.py \
    --output /app/data/combined_normalized.tsv

# Run annotation
docker run --rm -v $(pwd):/app pavs-pipeline \
    python normalization/annotate_variants.py \
    --input /app/data/combined_normalized.tsv \
    --output /app/data/combined_annotated.tsv
```

## Testing

```bash
# Run test suite
pytest tests/

# Validate phenopackets
uv run python intake/validate_shex.py \
    --phenopackets data/PAVS_phenopackets.json \
    --schema ontology/pavs_shapes.shex
```

## Documentation

- **[PIPELINE.md](docs/PIPELINE.md)** - Complete pipeline documentation with field schemas
- **[NORMALIZATION_METHODS.md](docs/NORMALIZATION_METHODS.md)** - Algorithm comparison and validation results
- **[RDF_SCHEMA.md](docs/RDF_SCHEMA.md)** - RDF schema and SPARQL query examples
- **[DATA_MODEL.md](docs/DATA_MODEL.md)** - Data model and ontology mappings

## Citation

If you use this data or pipeline, please cite:

```
Abdelhakim M, Althagafi A, Schofield PN, Hoehndorf R. PAVS: Phenotypic and 
Variant Standardization of Saudi Arabian rare disease patients. [In preparation]
```

## License

- **Code** (normalization/intake/scripts): GNU General Public License v3.0
- **Data** (data/phenotypes/, outputs): Creative Commons Attribution 4.0 International (CC-BY-4.0)

See [LICENSE](LICENSE) and [LICENSE-DATA](LICENSE-DATA) for details.

## Related Projects

- [PAVS](https://github.com/bio-ontology-research-group/pavs) - Main PAVS repository
- [Phenotype Matcher](https://github.com/bio-ontology-research-group/phenotype-matcher) - Normalization library
- [PAVS Website](https://github.com/bio-ontology-research-group/pavs-website) - Web interface
- [PAVS Phenopackets](https://github.com/bio-ontology-research-group/pavs-phenopackets) - Generated phenopackets
