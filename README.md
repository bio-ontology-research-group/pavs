# PAVS — Phenotype-Associated Variants in Saudi Arabia

A curated database of genomic variants and clinical phenotypes from Saudi Arabian patients with rare genetic diseases. PAVS integrates multiple clinical cohorts and links each case to international ontologies (HPO, OMIM, MONDO) and variant databases (ClinVar, dbSNP, TogoVar).

## Overview

PAVS consists of four integrated components, each maintained as a separate repository:

1. **[phenotype-matcher](https://github.com/bio-ontology-research-group/phenotype-matcher)** - Normalization library for clinical phenotype text → HPO/MONDO/OMIM
2. **[pavs-knowledge-graph](https://github.com/bio-ontology-research-group/pavs-knowledge-graph)** - Data pipeline: source data → phenopackets → RDF
3. **[pavs-phenopackets](https://github.com/bio-ontology-research-group/pavs-phenopackets)** - 7,566 GA4GH Phenopackets v2.0.2 JSON files
4. **[pavs-website](https://github.com/bio-ontology-research-group/pavs-website)** - Web interface with SPARQL backend for searching the knowledge graph

This repository serves as the umbrella project, coordinating all components and hosting analysis code, semantic similarity experiments, and documentation.

## Quick Start

### Clone with Submodules

```bash
git clone --recurse-submodules git@github.com:bio-ontology-research-group/pavs.git
cd pavs

# Or if already cloned:
git submodule update --init --recursive
```

### Run the Full Stack

```bash
# 1. Generate RDF data (knowledge graph pipeline)
cd knowledge-graph
uv sync
export OPENROUTER_API_KEY="your_key"
bash scripts/download_reference_data.sh
uv run python normalization/combine_normalize_phenotypes.py --output data/combined_normalized.tsv
uv run python normalization/annotate_variants.py --input data/combined_normalized.tsv --output data/combined_annotated.tsv
uv run python intake/generate_phenopackets_v2.py --input data/combined_annotated.tsv --combined data/PAVS_phenopackets.json
uv run python intake/compute_hpo_ic.py --hpo ontology/hp.obo --hpoa data/reference/phenotype.hpoa --output rdf_output/hpo_ic.ttl
uv run python intake/generate_rdf.py --phenopackets data/PAVS_phenopackets.json --annotated data/combined_annotated.tsv --output-dir rdf_output/
cd ..

# 2. Deploy website with Virtuoso triple store
cd website
cp .env.example .env
# Mount knowledge-graph/rdf_output as pavs_rdf volume in docker-compose.yml
docker compose up -d
cd ..

# Frontend: http://localhost:3000
# Backend API: http://localhost:8000/docs
# Virtuoso SPARQL: http://localhost:8890/sparql
```

## Repository Structure

```
pavs/
├── knowledge-graph/          [submodule] Data pipeline and source data
│   ├── data/phenotypes/      - 7 source TSV files
│   ├── normalization/        - Phenotype/variant normalization scripts
│   ├── intake/               - Phenopacket and RDF generation
│   └── docs/                 - Pipeline documentation
│
├── phenopackets-data/        [submodule] Generated phenopackets (outputs)
│   ├── individual/           - 7,566 individual JSON files
│   └── PAVS_phenopackets.json - Combined collection
│
├── website/                  [submodule] Web interface
│   ├── backend/              - FastAPI + SPARQL + similarity search
│   ├── frontend/             - React + TypeScript UI
│   └── virtuoso/             - Virtuoso bulk load scripts
│
├── tools/phenotype_matcher_v2/  [submodule] Normalization library
│   ├── src/phenotype_matcher_v2/  - Core library code
│   └── ALGORITHM.md          - Detailed algorithm documentation
│
├── analysis/                 Semantic similarity experiments (THIS REPO)
│   ├── similarity.py         - Lin/Resnik similarity implementations
│   ├── validate_phenopackets.py - Phenopacket validation
│   └── IC_COMPARISON.md      - Information content analysis
│
├── translation/              HPO Arabic translations (THIS REPO)
│   └── hpo_arabic_translations.json
│
├── onto-normalize/           Normalization benchmarks (THIS REPO)
│   └── BENCHMARK_REPORT.md   - Comparison of normalization methods
│
├── obsolete/                 Legacy code and deprecated approaches
│   └── moved_to_submodules/  - Original content before split
│
└── CLAUDE.md                 AI assistant instructions
```

## Data Summary

- **Total cases**: 7,566
  - Saudi rare disease patients: 5,710 (75.5%)
  - Non-Saudi DDD controls: 1,856 (24.5%)
- **Unique HPO phenotypes**: 2,847
- **Unique genes**: 1,234
- **Confirmed diagnoses**: 4,523 (59.8%)
- **Genome build**: GRCh38

### Case ID Format

`PAVS:XNNNNNNN` where X is the source letter:

| Letter | Source | Cases | Saudi |
|---|---|---|---|
| A, B | ahmed-pmid28454995 | 525 | Yes |
| F | fawzan-variants | 1,024 | Yes |
| M | marwa-variants | 1,421 | Yes |
| P | PMC6562004 | 2,218 | Yes |
| Q | PMC7082194 | 522 | Yes |
| D | ddd-diagnoses | 1,856 | No |

## Semantic Similarity Analysis

This repository hosts experiments comparing phenotype similarity algorithms (Lin, Resnik, Jiang-Conrath) for rare disease diagnosis. See `analysis/` directory.

**Key findings:**
- Lin + BMA provides best balance of precision/recall
- IC values computed from disease annotations vs. case frequencies show minimal difference
- Ancestor expansion improves similarity scores by 8-12%

See [analysis/IC_COMPARISON.md](analysis/IC_COMPARISON.md) for details.

## Normalization Algorithm

The phenotype matching pipeline uses a 6-stage algorithm:

1. Text preprocessing (Unicode, lowercasing, stopwords)
2. Exact matching (ontology label indices)
3. Synonym expansion
4. Fuzzy matching (RapidFuzz for typos)
5. Semantic similarity (SapBERT embeddings)
6. Ranking & filtering

**Performance** (validated on 7,500+ cases):
- Precision: 94.2%
- Recall: 87.6%
- Speed: ~150 terms/second

See [tools/phenotype_matcher_v2/ALGORITHM.md](tools/phenotype_matcher_v2/ALGORITHM.md) for full details.

## Architecture

```
┌─────────────────┐         ┌──────────────────┐         ┌────────────────┐
│  Source Data    │         │  Normalization   │         │  Phenopackets  │
│  (7 TSV files)  │────────▶│  Pipeline        │────────▶│  (GA4GH v2)    │
│                 │         │  (phenotype-     │         │                │
└─────────────────┘         │   matcher)       │         └────────────────┘
                            └──────────────────┘                 │
                                                                 ▼
┌─────────────────┐         ┌──────────────────┐         ┌────────────────┐
│   Web UI        │         │   FastAPI        │         │   Virtuoso     │
│   (React +      │◀───────▶│   Backend        │◀───────▶│   RDF Store    │
│    TypeScript)  │         │   (SPARQL)       │         │   (5 graphs)   │
└─────────────────┘         └──────────────────┘         └────────────────┘
```

## Development

### Prerequisites

- Python 3.10+ with `uv` package manager
- Docker + Docker Compose (for website deployment)
- Node.js 18+ (for frontend development)
- OpenRouter API key (for normalization)

### Working with Submodules

```bash
# Update all submodules to latest
git submodule update --remote

# Make changes in a submodule
cd knowledge-graph
git checkout -b my-feature
# ... make changes ...
git add . && git commit -m "Add feature"
git push origin my-feature

# Update main repo to point to new submodule commit
cd ..
git add knowledge-graph
git commit -m "Update knowledge-graph submodule"
git push
```

### Running Tests

```bash
# Test normalization library
cd tools/phenotype_matcher_v2
pytest

# Validate phenopackets
cd knowledge-graph
uv run python intake/validate_shex.py --phenopackets data/PAVS_phenopackets.json

# Test SPARQL queries
cd website
pytest tests/test_sparql_queries.py
```

## Documentation

- **[CLAUDE.md](CLAUDE.md)** - Instructions for AI assistants working with this codebase
- **[knowledge-graph/docs/PIPELINE.md](knowledge-graph/docs/PIPELINE.md)** - Complete data pipeline documentation
- **[knowledge-graph/docs/NORMALIZATION_METHODS.md](knowledge-graph/docs/NORMALIZATION_METHODS.md)** - Algorithm comparison
- **[knowledge-graph/docs/RDF_SCHEMA.md](knowledge-graph/docs/RDF_SCHEMA.md)** - RDF schema and SPARQL examples
- **[analysis/IC_COMPARISON.md](analysis/IC_COMPARISON.md)** - Semantic similarity analysis
- **[onto-normalize/BENCHMARK_REPORT.md](onto-normalize/BENCHMARK_REPORT.md)** - Normalization benchmarks

## Citation

If you use PAVS in your research, please cite:

```

```

## Contributors

| Name | ORCID | Role |
|---|---|---|
| **Marwa Abdelhakim** | [0000-0001-6816-2119](https://orcid.org/0000-0001-6816-2119) | Manual curation of all source datasets |
| **Azza Althagafi** | [0000-0001-6084-8706](https://orcid.org/0000-0001-6084-8706) | Computational experiments and validation |
| **Paul N Schofield** | [0000-0002-5111-7263](https://orcid.org/0000-0002-5111-7263) | Curation, supervision, and guidance |
| **Robert Hoehndorf** | [0000-0001-8149-5890](https://orcid.org/0000-0001-8149-5890) | Supervision, code contributions, validation analysis, website |

## License

- **Code** (analysis/, scripts in this repo): GNU General Public License v3.0
- **Data** (phenopackets, source data): Creative Commons Attribution 4.0 International (CC-BY-4.0)
- See individual submodule repositories for their specific licenses

## Contact

For questions or issues:
- Open an issue in the relevant repository
- Email: robert.hoehndorf@kaust.edu.sa

## Related Projects

- [HPO Arabic Translation](https://github.com/bio-ontology-research-group/hpo-arabic) - Human Phenotype Ontology Arabic translations (submodule in `translation/`)
- [Monarch Initiative](https://monarchinitiative.org/) - Integrated phenotype data across species
- [GA4GH Phenopackets](https://github.com/phenopackets/phenopacket-schema) - Standard for sharing phenotypic data
