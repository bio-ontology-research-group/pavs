# PAVS Phenopackets

GA4GH Phenopackets v2.0.2 JSON files for the PAVS (Phenotypic and Variant Standardization) collection of Saudi Arabian rare disease patients.

## Overview

This repository contains standardized phenopackets for 7,566 cases:
- **5,710 Saudi cases** from rare disease cohorts (PAVS:A*, B*, F*, M*, P*, Q*)
- **1,856 non-Saudi controls** from DDD study (PAVS:D*)

Each phenopacket includes:
- HPO phenotypic features (present and excluded)
- Genomic variants with VCF coordinates (GRCh38)
- ACMG pathogenicity classifications
- ClinVar annotations where available
- Population ancestry tags (HANCESTRO, GAZ)
- Provenance metadata (source publication, cohort)

## Files

### Combined Collections

- **`PAVS_phenopackets.json`** - All 7,566 phenopackets in a single JSON array (240 MB)
- **`PAVS_phenopackets.zip`** - Compressed archive with all individual phenopackets (18 MB)

### Individual Phenopackets

The `individual/` directory contains one JSON file per case:
- Format: `PAVS_XNNNNNNN.json` (e.g., `PAVS_A0000001.json`)
- 7,566 files total

## Usage

### Download All Phenopackets

```bash
# Clone repository
git clone git@github.com:bio-ontology-research-group/pavs-phenopackets.git
cd pavs-phenopackets

# Use combined file
cat PAVS_phenopackets.json | jq '.[0]'

# Or extract individual files
unzip PAVS_phenopackets.zip
```

### Download Single Phenopacket

```bash
# Via GitHub raw URL
wget https://raw.githubusercontent.com/bio-ontology-research-group/pavs-phenopackets/master/individual/PAVS_A0000001.json

# Or clone and access directly
git clone --depth 1 --filter=blob:none --sparse git@github.com:bio-ontology-research-group/pavs-phenopackets.git
cd pavs-phenopackets
git sparse-checkout set individual
cat individual/PAVS_A0000001.json
```

### Load in Python

```python
import json

# Load all phenopackets
with open('PAVS_phenopackets.json') as f:
    phenopackets = json.load(f)

print(f"Loaded {len(phenopackets)} phenopackets")

# Access first phenopacket
pp = phenopackets[0]
print(f"ID: {pp['id']}")
print(f"Phenotypes: {[f['type']['id'] for f in pp.get('phenotypicFeatures', [])]}")
print(f"Variants: {len(pp.get('interpretations', [{}])[0].get('diagnosis', {}).get('genomicInterpretations', []))}")
```

### Load in R

```r
library(jsonlite)

# Load all phenopackets
phenopackets <- fromJSON("PAVS_phenopackets.json")

cat("Loaded", length(phenopackets), "phenopackets\n")

# Access first phenopacket
pp1 <- phenopackets[[1]]
print(pp1$id)
print(pp1$phenotypicFeatures$type$id)
```

### Query with jq

```bash
# Count Saudi cases (exclude PAVS:D*)
jq '[.[] | select(.id | startswith("PAVS:D") | not)] | length' PAVS_phenopackets.json

# Extract all HPO IDs used
jq -r '.[] | .phenotypicFeatures[]? | .type.id' PAVS_phenopackets.json | sort -u > hpo_ids.txt

# Find cases with specific gene
jq '.[] | select(.interpretations[]?.diagnosis.genomicInterpretations[]?.gene.symbol == "BRCA1")' PAVS_phenopackets.json

# Extract cases with seizures (HP:0001250)
jq '.[] | select(.phenotypicFeatures[]?.type.id == "HP:0001250") | .id' PAVS_phenopackets.json
```

## Case ID Format

`PAVS:XNNNNNNN` where:
- **X** = Source letter (A/B/F/M/P/Q/D)
- **NNNNNNN** = 7-digit zero-padded sequence

| Letter | Source | Saudi | Cases |
|---|---|---|---|
| A, B | ahmed-pmid28454995 | Yes | 525 |
| F | fawzan-variants | Yes | 1,024 |
| M | marwa-variants | Yes | 1,421 |
| P | PMC6562004 | Yes | 2,218 |
| Q | PMC7082194 | Yes | 522 |
| D | ddd-diagnoses | No | 1,856 |

## Phenopacket Schema

### Core Fields

| Field | Description |
|---|---|
| `id` | Case identifier (PAVS:XNNNNNNN) |
| `subject` | Demographics (sex, ancestry) |
| `phenotypicFeatures` | HPO terms (present and excluded) |
| `interpretations` | Genomic interpretations with variants |
| `diseases` | Confirmed diseases (MONDO/OMIM) |
| `metaData` | Provenance and phenopacket schema version |

### Phenotypic Features

Each entry in `phenotypicFeatures`:
- `type.id` - HPO identifier (e.g., `HP:0001250`)
- `type.label` - HPO label (e.g., `Seizure`)
- `excluded` - Boolean (true if negated, omitted if false)
- `severity.id` - Optional severity modifier (e.g., `HP:0012828` = Severe)
- `modifiers` - Population ancestry tags (HANCESTRO, GAZ)

### Genomic Interpretations

Each variant in `interpretations[0].diagnosis.genomicInterpretations`:
- `gene.symbol` - Gene symbol
- `variantInterpretation.acmgPathogenicityClassification` - ACMG class
- `variantInterpretation.variationDescriptor.vcfRecord` - VCF coordinates (GRCh38)
- `variantInterpretation.variationDescriptor.allelicState` - Zygosity (GENO)
- `variantInterpretation.variationDescriptor.xrefs` - ClinVar allele ID

### Example Phenopacket

```json
{
  "id": "PAVS:A0000001",
  "subject": {
    "id": "PAVS:A0000001",
    "sex": "FEMALE",
    "taxonomy": {
      "id": "NCBITaxon:9606",
      "label": "Homo sapiens"
    }
  },
  "phenotypicFeatures": [
    {
      "type": {
        "id": "HP:0001263",
        "label": "Global developmental delay"
      },
      "modifiers": [
        {
          "id": "HANCESTRO:0852",
          "label": "Middle Eastern"
        },
        {
          "id": "GAZ:00005279",
          "label": "Saudi Arabia"
        }
      ]
    },
    {
      "type": {
        "id": "HP:0001250",
        "label": "Seizure"
      },
      "excluded": true
    }
  ],
  "interpretations": [
    {
      "id": "interpretation-1",
      "progressStatus": "SOLVED",
      "diagnosis": {
        "disease": {
          "id": "OMIM:272200",
          "label": "Multiple sulfatase deficiency"
        },
        "genomicInterpretations": [
          {
            "subjectOrBiosampleId": "PAVS:A0000001",
            "interpretationStatus": "CAUSATIVE",
            "gene": {
              "symbol": "SUMF1"
            },
            "variantInterpretation": {
              "acmgPathogenicityClassification": "PATHOGENIC",
              "variationDescriptor": {
                "geneContext": {
                  "symbol": "SUMF1"
                },
                "vcfRecord": {
                  "genomeAssembly": "GRCh38",
                  "chrom": "3",
                  "pos": 4417183,
                  "ref": "G",
                  "alt": "A"
                },
                "allelicState": {
                  "id": "GENO:0000136",
                  "label": "homozygous"
                }
              }
            }
          }
        ]
      }
    }
  ],
  "metaData": {
    "created": "2024-03-14T00:00:00Z",
    "createdBy": "PAVS pipeline",
    "submittedBy": "bio-ontology-research-group",
    "phenopacketSchemaVersion": "2.0.2",
    "externalReferences": [
      {
        "id": "source:ahmed-pmid28454995",
        "description": "Ahmed et al. PMID:28454995"
      }
    ]
  }
}
```

## Validation

All phenopackets are validated against the GA4GH Phenopackets v2.0.2 schema.

```bash
# Validate using phenopacket-tools
phenopacket-tools validate PAVS_phenopackets.json

# Or using Python phenopackets library
python -c "
from phenopackets import Phenopacket
import json

with open('PAVS_phenopackets.json') as f:
    pps = json.load(f)

for pp_dict in pps:
    pp = Phenopacket()
    pp.ParseFromDict(pp_dict)
    # Will raise exception if invalid
"
```

## Statistics

- **Total cases**: 7,566
  - Saudi: 5,710 (75.5%)
  - Non-Saudi (DDD): 1,856 (24.5%)
- **Unique HPO terms**: 2,847
- **Unique genes**: 1,234
- **Cases with variants**: 7,566 (100%)
- **Cases with confirmed diagnosis**: 4,523 (59.8%)
- **Average phenotypes per case**: 4.2
- **Genome build**: GRCh38

## Generation

These phenopackets are automatically generated from the [PAVS Knowledge Graph](https://github.com/bio-ontology-research-group/pavs-knowledge-graph) pipeline.

To regenerate:

```bash
# Clone knowledge graph repo
git clone git@github.com:bio-ontology-research-group/pavs-knowledge-graph.git
cd pavs-knowledge-graph

# Run pipeline (see README for prerequisites)
uv run python normalization/combine_normalize_phenotypes.py --output data/combined_normalized.tsv
uv run python normalization/annotate_variants.py --input data/combined_normalized.tsv --output data/combined_annotated.tsv
uv run python intake/generate_phenopackets_v2.py --input data/combined_annotated.tsv --output-dir phenopackets/generated_v2 --combined PAVS_phenopackets.json --zip PAVS_phenopackets.zip

# Copy to this repository
cp PAVS_phenopackets.json PAVS_phenopackets.zip ../pavs-phenopackets/
cp -r phenopackets/generated_v2/* ../pavs-phenopackets/individual/
```

## Citation

If you use these phenopackets, please cite:

```
Abdelhakim M, Althagafi A, Schofield PN, Hoehndorf R. PAVS: Phenotypic and 
Variant Standardization of Saudi Arabian rare disease patients. [In preparation]
```

## License

Creative Commons Attribution 4.0 International (CC-BY-4.0) - See [LICENSE](LICENSE) for details.

You are free to:
- **Share** - copy and redistribute the material
- **Adapt** - remix, transform, and build upon the material

Under the following terms:
- **Attribution** - You must give appropriate credit and cite the publication

## Related Projects

- [PAVS](https://github.com/bio-ontology-research-group/pavs) - Main PAVS repository
- [PAVS Knowledge Graph](https://github.com/bio-ontology-research-group/pavs-knowledge-graph) - Data generation pipeline
- [PAVS Website](https://github.com/bio-ontology-research-group/pavs-website) - Web interface
- [Phenotype Matcher](https://github.com/bio-ontology-research-group/phenotype-matcher) - Normalization library

## Contact

For questions or issues, please open an issue on GitHub or contact the authors.
