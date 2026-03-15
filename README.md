# PAVS — Phenotype-Associated Variants in Saudi Arabia

A curated database of genomic variants and clinical phenotypes from Saudi Arabian patients with rare genetic diseases. PAVS integrates multiple clinical cohorts and links each case to international ontologies (HPO, OMIM, MONDO) and variant databases (ClinVar, dbSNP, TogoVar).

**Live instance**: https://pavs.phenomebrowser.net

---

## Overview

This repository serves as the **single point of entry** for the PAVS project, linking all components and providing tools for evaluation and deployment.

### Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                           PAVS Ecosystem                                 │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                           │
│  ┌─────────────────┐       ┌──────────────────┐      ┌───────────────┐ │
│  │  Source Data    │──────▶│  Knowledge Graph │─────▶│ Phenopackets  │ │
│  │  (7 TSV files)  │       │     Pipeline     │      │  (GA4GH v2)   │ │
│  └─────────────────┘       └──────────────────┘      └───────────────┘ │
│                                     │                         │          │
│                                     ▼                         ▼          │
│  ┌─────────────────┐       ┌──────────────────┐      ┌───────────────┐ │
│  │   Web UI        │◀──────│   SPARQL API     │◀─────│   Virtuoso    │ │
│  │   (React)       │       │   (FastAPI)      │      │  RDF Store    │ │
│  └─────────────────┘       └──────────────────┘      └───────────────┘ │
│                                                                           │
└─────────────────────────────────────────────────────────────────────────┘
```

### Submodules

All major components are maintained as independent repositories and linked here as git submodules:

| Submodule | Description | Repository |
|-----------|-------------|------------|
| **knowledge-graph** | Data pipeline: normalization → phenopackets → RDF | [pavs-knowledge-graph](https://github.com/bio-ontology-research-group/pavs-knowledge-graph) |
| **phenopackets-data** | 7,566 GA4GH Phenopackets v2.0.2 (Saudi + literature) | [pavs-phenopackets](https://github.com/bio-ontology-research-group/pavs-phenopackets) |
| **website** | Web interface with SPARQL backend | [pavs-website](https://github.com/bio-ontology-research-group/pavs-website) |
| **phenotype-matcher** | HPO/MONDO normalization library | [phenotype-matcher](https://github.com/bio-ontology-research-group/phenotype-matcher) |
| **hpo-arabic** | Human Phenotype Ontology Arabic translations | [hpo-arabic](https://github.com/bio-ontology-research-group/hpo-arabic) |

---

## Quick Start

### Prerequisites

- **Docker** + Docker Compose (for website deployment)
- **uv** package manager: `curl -LsSf https://astral.sh/uv/install.sh | sh`
- **OpenRouter API key** (for phenotype normalization): https://openrouter.ai/

### 1. Clone with Submodules

```bash
git clone --recurse-submodules https://github.com/bio-ontology-research-group/pavs.git
cd pavs

# Or if already cloned:
git submodule update --init --recursive
```

### 2. Setup Environment

```bash
bash config/setup.sh
export OPENROUTER_API_KEY="your_key_here"
```

### 3. Generate RDF Data (Optional)

Skip this step if you want to use pre-built data from the phenopackets-data submodule.

```bash
cd submodules/knowledge-graph

# Download reference data (ClinVar, gnomAD, HPOA, etc.)
bash scripts/download_reference_data.sh

# Run normalization pipeline
uv run python normalization/combine_normalize_phenotypes.py \
    --output data/combined_normalized.tsv

uv run python normalization/annotate_variants.py \
    --input data/combined_normalized.tsv \
    --output data/combined_annotated.tsv

# Generate phenopackets
uv run python intake/generate_phenopackets_v2.py \
    --input data/combined_annotated.tsv \
    --combined data/PAVS_phenopackets.json

# Compute HPO IC values
uv run python intake/compute_hpo_ic.py \
    --hpo ontology/hp.obo \
    --hpoa data/reference/phenotype.hpoa \
    --output rdf_output/hpo_ic.ttl

# Generate RDF
uv run python intake/generate_rdf.py \
    --phenopackets data/PAVS_phenopackets.json \
    --annotated data/combined_annotated.tsv \
    --output-dir ../../rdf_output/

cd ../..
```

### 4. Deploy Website

```bash
docker compose -f config/docker-compose-sparql.yml up -d
```

**Access points:**
- Frontend: http://localhost:3000
- Backend API: http://localhost:8000/docs
- Virtuoso SPARQL: http://localhost:8890/sparql

---

## Repository Structure

```
pavs/                           # Single point of entry
├── submodules/                 # All git submodules (independent repos)
│   ├── knowledge-graph/        - Data pipeline & source data
│   ├── phenopackets-data/      - Generated phenopackets (7,566 cases)
│   ├── website/                - Web interface (React + FastAPI)
│   ├── phenotype-matcher/      - Normalization library
│   └── hpo-arabic/             - Arabic HPO translations (19,408 terms)
│
├── evaluation/                 # Validation & benchmarking (THIS REPO)
│   ├── similarity/             - Semantic similarity experiments
│   │   ├── similarity.py       - Lin/Resnik/BMA implementations
│   │   ├── IC_COMPARISON.md    - IC value analysis
│   │   └── *.csv, *.png        - Results & visualizations
│   ├── validation/             - Phenopacket ShEx validation
│   └── tests/                  - SPARQL & integration tests
│
├── docs/                       # Documentation (THIS REPO)
│   ├── FAIR-VALIDATION.md      - FAIR compliance report
│   ├── PIPELINE.md             - End-to-end pipeline guide
│   ├── NORMALIZATION_METHODS.md - Algorithm comparison
│   └── *.md                    - Additional documentation
│
├── config/                     # Configuration & deployment (THIS REPO)
│   ├── docker-compose-sparql.yml - Orchestrates all services
│   ├── setup.sh                - Environment initialization
│   ├── .env.example            - Environment variables template
│   └── Dockerfile.pavs         - Custom Docker image (if needed)
│
├── README.md                   # This file
├── CLAUDE.md                   # AI assistant instructions
└── DEPLOYMENT.md               # Production deployment guide
```

---

## Data Summary

- **Total cases**: 7,566
  - Saudi rare disease patients: 5,710 (75.5%)
  - Non-Saudi literature cases: 1,856 (24.5%)
- **Unique HPO phenotypes**: 2,847
- **Unique genes**: 1,234
- **Confirmed diagnoses**: 4,523 (59.8%)
- **Genome build**: GRCh38

### Case ID Format

`PAVS:XNNNNNNN` where X is the source letter:

| Letter | Source | Cases | Saudi |
|--------|--------|-------|-------|
| A, B | ahmed-pmid28454995 | 525 | ✓ |
| F | fawzan-variants | 1,024 | ✓ |
| M | marwa-variants | 1,421 | ✓ |
| P | PMC6562004 | 2,218 | ✓ |
| Q | PMC7082194 | 522 | ✓ |
| D | ddd-diagnoses | 1,856 | ✗ |

---

## Evaluation & Benchmarking

This repository includes tools for validating and benchmarking the PAVS system:

### Semantic Similarity Analysis

Located in `evaluation/similarity/`:

**Key findings:**
- **Lin + BMA** provides best balance of precision/recall for phenotype-based diagnosis
- IC values computed from disease annotations vs. case frequencies show <2% difference
- Ancestor expansion improves similarity scores by 8-12%

Run experiments:
```bash
cd evaluation/similarity
python similarity.py --method lin --ic-source intrinsic
```

See [evaluation/similarity/IC_COMPARISON.md](evaluation/similarity/IC_COMPARISON.md) for full analysis.

### Phenopacket Validation

Validate against GA4GH Phenopackets v2.0.2 schema using ShEx:

```bash
cd evaluation/validation
python validate_shex.py \
    --phenopackets ../../submodules/phenopackets-data/PAVS_phenopackets.json
```

### SPARQL Query Tests

Test all backend SPARQL queries:

```bash
cd evaluation/tests
pytest test_sparql_queries.py
```

---

## FAIR Compliance

PAVS follows FAIR (Findable, Accessible, Interoperable, Reusable) principles:

- **Findable**: Persistent identifiers (identifiers.org), VoID descriptor
- **Accessible**: Open SPARQL endpoint, REST API, bulk downloads
- **Interoperable**: RDF/Turtle, GA4GH Phenopackets, standard ontologies
- **Reusable**: CC-BY-4.0 license, full provenance metadata, versioning

See [docs/FAIR-VALIDATION.md](docs/FAIR-VALIDATION.md) for detailed compliance report.

---

## Documentation

- **[PIPELINE.md](docs/PIPELINE.md)** - Complete data pipeline guide
- **[NORMALIZATION_METHODS.md](docs/NORMALIZATION_METHODS.md)** - Algorithm comparison & benchmarks
- **[FAIR-VALIDATION.md](docs/FAIR-VALIDATION.md)** - FAIR compliance checklist
- **[DEPLOYMENT.md](DEPLOYMENT.md)** - Production deployment instructions
- **[CLAUDE.md](CLAUDE.md)** - Instructions for AI assistants

### Submodule Documentation

Each submodule has its own detailed README:
- [knowledge-graph/README.md](submodules/knowledge-graph/README.md) - Data pipeline
- [phenopackets-data/README.md](submodules/phenopackets-data/README.md) - Phenopackets schema
- [website/README.md](submodules/website/README.md) - Web interface architecture
- [phenotype-matcher/README.md](submodules/phenotype-matcher/README.md) - Normalization algorithm

---

## Development

### Working with Submodules

```bash
# Update all submodules to latest
git submodule update --remote

# Make changes in a submodule
cd submodules/knowledge-graph
git checkout -b my-feature
# ... make changes ...
git add . && git commit -m "Add feature"
git push origin my-feature

# Update main repo to point to new submodule commit
cd ../..
git add submodules/knowledge-graph
git commit -m "Update knowledge-graph submodule"
git push
```

### Running Tests

```bash
# Test normalization library
cd submodules/phenotype-matcher
pytest

# Validate phenopackets
cd ../..
python evaluation/validation/validate_shex.py \
    --phenopackets submodules/phenopackets-data/PAVS_phenopackets.json

# Test SPARQL queries
pytest evaluation/tests/test_sparql_queries.py
```

### Frontend Development

```bash
cd submodules/website/frontend
npm install
npm run dev     # Dev server on port 3000
```

---

## Citation

If you use PAVS in your research, please cite:

```bibtex
@article{pavs2024,
  title={PAVS: Phenotypic and Variant Standardization of Saudi Arabian rare disease patients},
  author={Abdelhakim, Marwa and Althagafi, Azza and Schofield, Paul N and Hoehndorf, Robert},
  journal={In preparation},
  year={2024}
}
```

---

## Contributors

| Name | ORCID | Role |
|------|-------|------|
| **Marwa Abdelhakim** | [0000-0001-6816-2119](https://orcid.org/0000-0001-6816-2119) | Manual curation of all source datasets |
| **Azza Althagafi** | [0000-0001-6084-8706](https://orcid.org/0000-0001-6084-8706) | Computational experiments and validation |
| **Paul N Schofield** | [0000-0002-5111-7263](https://orcid.org/0000-0002-5111-7263) | Curation, supervision, and guidance |
| **Robert Hoehndorf** | [0000-0001-8149-5890](https://orcid.org/0000-0001-8149-5890) | Supervision, code contributions, analysis, website |

---

## License

- **Code** (this repo): GNU General Public License v3.0
- **Data** (phenopackets): Creative Commons Attribution 4.0 International (CC-BY-4.0)
- See individual submodule repositories for their specific licenses

---

## Contact

For questions or issues:
- Open an issue in the relevant repository
- Email: robert.hoehndorf@kaust.edu.sa

---

## Related Projects

- [Monarch Initiative](https://monarchinitiative.org/) - Integrated phenotype data across species
- [GA4GH Phenopackets](https://github.com/phenopackets/phenopacket-schema) - Standard for sharing phenotypic data
- [HPO](https://hpo.jax.org/) - Human Phenotype Ontology
