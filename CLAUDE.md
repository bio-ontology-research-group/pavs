# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PAVS (Phenotypic and Variant Standardization) is a curated database of genomic variants and clinical phenotypes from Saudi Arabian patients with rare genetic diseases. The SPARQL-based web stack is the production system; all data lives in a Virtuoso triple store.

Always use `uv` to run Python scripts.

## Repository Structure

```
normalization/     # Data normalization pipeline (normalize_*.py, normalization_utils.py)
intake/            # Format conversion (phenopackets, RDF, Virtuoso loading)
website/
  backend/         # FastAPI SPARQL proxy (was backend_sparql/)
  frontend/        # React + TypeScript frontend (was frontend_sparql/)
  virtuoso/        # Virtuoso isql bulk-load script
analysis/          # Phenotype similarity analysis, validation, plots
tools/
  phenotype_matcher_v2/  # Standalone phenotype matching library
docs/              # All project documentation
scripts/           # Legacy utility scripts
data/              # Source data and reference files
  phenotypes/      # 7 source TSV files (intake only from here)
ontology/          # External reference ontologies
rdf_output/        # Generated RDF/Turtle files
phenopackets/      # Generated phenopacket JSON files
obsolete/          # Deprecated code (SQLite stack, phenotype_matcher v1)
```

## Key Conventions

- **Case IDs**: `PAVS:XNNNNNNN` (X = source letter A/B/F/M/P/Q/D, 7-digit seq)
- **Genome build**: GRCh38 for VCF coordinates; GRCh37 for legacy HGVS strings
- **Population tags**: Saudi cases use `HANCESTRO:0852` (Middle Eastern) and `GAZ:00005279` (Saudi Arabia)
- **Provenance**: Maintain `source_file` and `source_id` in all derived datasets
- **Multi-value fields**: pipe-delimited in TSV outputs
- **Data intake**: Only from files in `data/phenotypes/` (7 source TSVs)
- **Reference data**: External databases in `data/reference/` (phenotype.hpoa, clinvar, gnomad, etc.) — download with `bash scripts/download_reference_data.sh`

## Running the SPARQL Stack

```bash
# Full stack (Virtuoso + loader + backend + frontend):
docker compose -f docker-compose-sparql.yml up -d

# After code changes to website/backend/ or website/frontend/ only:
docker compose -f docker-compose-sparql.yml build backend frontend
docker compose -f docker-compose-sparql.yml up -d backend frontend

# Health check:
curl http://localhost:8000/api/health

# Frontend: http://localhost:3000
# Backend API docs: http://localhost:8000/docs
# Virtuoso SPARQL UI: http://localhost:8890/sparql
```

Backend startup takes ~15 seconds (loads IC values, ancestor cache, case HPO sets, disease labels from HPOA).

## SPARQL Stack Architecture

### Named graphs

| Graph URI | Contents |
|---|---|
| `graph/cases` | Saudi PAVS cases: phenotypes, variants, demographics |
| `graph/genes` | Gene→disease, pLI/LOEUF, GO, GTEx |
| `graph/hpoa` | Disease→HPO from `data/phenotype.hpoa` |
| `graph/hpo-ic` | IC values + HPO hierarchy + `rdfs:label` (needed for autocomplete) |
| `graph/literature` | 9,588 non-Saudi literature phenopackets |

### In-memory caches (loaded at startup)

- `ic_cache` — `{HP:NNNNNNN → float}` from graph/hpo-ic
- `ancestor_cache` — `{HP:NNNNNNN → set of ancestor HP IDs}` (transitive closure)
- `case_hpo_cache` — list of `{id, hpo_ids, is_saudi, gene, disease, source}` for similarity scoring
- `disease_label_cache` — `{"OMIM:272200" → "Multiple sulfatase deficiency"}` loaded from `data/phenotype.hpoa`

### Key files

| File | Role |
|---|---|
| `website/backend/main.py` | FastAPI app, startup caches, all endpoints |
| `website/backend/sparql_queries.py` | All SPARQL query templates |
| `website/backend/similarity.py` | Lin/Resnik + BMA (pure Python, no pyhpo) |
| `intake/generate_rdf.py` | Converts source data → Turtle (.ttl) |
| `intake/compute_hpo_ic.py` | Computes HPO IC values → `rdf_output/hpo_ic.ttl` |
| `website/virtuoso/load_ttl.sql` | isql bulk-load script (`ld_dir + rdf_loader_run`) |
| `docker-compose-sparql.yml` | Orchestrates virtuoso + init + loader + backend + frontend |
| `website/frontend/src/components/` | React components for each tab |
| `website/frontend/nginx.conf` | Proxies `/api/` → backend:8000 |

### Docker service order

`virtuoso` → `init` (generate RDF, --skip-load) → `loader` (isql bulk load) → `backend` → `frontend`

## SPARQL Stack — API Endpoints

| Endpoint | Notes |
|---|---|
| `GET /api/health` | Status check |
| `GET /api/search/hpo?q=` | HPO autocomplete via rdfs:label in graph/hpo-ic |
| `POST /api/search/phenotype` | Lin+BMA similarity, in-memory; body: `{hpo_ids, method, limit, include_disease_phenotypes, include_non_saudi}` |
| `GET /api/search/gene?q=` | Gene symbol substring search |
| `GET /api/search/variant` | By gene / rsid / hgvs / acmg |
| `GET /api/search/disease?q=` | Disease label substring |
| `GET /api/case/{id}` | Full case: properties, phenotypes (with labels), variants, suggested_diseases |
| `GET /api/togovar-search?chrom=&pos=` | Proxies to TogoVar POST API, 302 redirect to variant page |
| `GET /api/phenopacket/{id}/download` | Individual phenopacket JSON |
| `GET /api/phenopackets/download-all` | PAVS_phenopackets.zip |

### Case detail response — key fields

- `phenotypes`: `[{id, label}]` — HPO terms with labels from graph/hpo-ic (all terms, not just the first)
- `excluded_phenotypes`: same format for negated HPO terms
- `suggested_diseases`: `[{id, label, url, gene}]` — shown only when `diseaseLabel` is absent AND gene has 1–3 OMIM associations; labels from `disease_label_cache`
- `properties.progressStatus`: mapped from Phenopackets v2 enum (SOLVED→"Solved", UNSOLVED→"Unsolved", IN_PROGRESS/UNKNOWN_PROGRESS→"Unknown")
- `properties.isSaudi`: normalised to `"true"/"false"` string (Virtuoso returns `"1"` for `true^^xsd:boolean`)
- `variants[].togovar_url`: path `/api/togovar-search?chrom=X&pos=Y` (proxied by nginx → backend)

## SPARQL Stack — Known Issues / Gotchas

- **Virtuoso boolean encoding**: SPARQL JSON returns `"1"` not `"true"` for `true^^xsd:boolean`. Code normalises with `in ("1", "true")`.
- **Turtle syntax**: bare booleans (`true`/`false`) and numeric literals cannot have `^^xsd:type` annotations in Turtle. Generate RDF without them.
- **Virtuoso DirsAllowed**: The `/rdf_import` mount must be in `VIRT_Parameters_DirsAllowed` env var for bulk loading. On existing volumes, set it with Virtuoso's `inifile` tool.
- **10 MB HTTP limit**: Virtuoso's SPARQL HTTP endpoint rejects payloads > 10 MB. Use `isql` + `ld_dir + rdf_loader_run` for large TTL files.
- **Lin similarity > 1.0**: Can occur if ancestor IC values are inconsistent with child IC values due to imperfect disease annotation propagation in phenotype.hpoa. Scores are displayed as raw decimals (not %) to avoid confusion.
- **TogoVar**: Uses `POST https://grch38.togovar.org/api/search/variant` with `{"query":{"location":{"chromosome":"16","position":48258198}}}`. Chromosome must be bare number string ("16", not "chr16"). Saudi-specific variants are often absent from the Japanese cohort database.

## RDF Generation Pipeline

```bash
# Step 1: HPO IC values + hierarchy + labels
uv run python intake/compute_hpo_ic.py \
    --hpo ontology/hp.obo --hpoa data/phenotype.hpoa \
    --output rdf_output/hpo_ic.ttl

# Step 2: All case/gene/hpoa/literature RDF
uv run python intake/generate_rdf.py \
    --phenopackets data/PAVS_phenopackets.json \
    --annotated data/combined_annotated.tsv \
    --literature-dir phenopackets/0.1.26 \
    --output-dir rdf_output/

# Step 3: Load (handled by Docker loader service, or manually):
# isql virtuoso:1111 dba pavs_dba /load/load_ttl.sql
```

Output TTL file sizes: cases ~192K lines, genes ~103K, hpoa ~2.2M, hpo_ic ~56K, literature ~455K.

## Data Processing Pipeline (upstream)

```bash
uv sync
export OPENROUTER_API_KEY="your_key_here"

# Combine and normalise all 7 source TSVs → combined_normalized.tsv (38 columns)
uv run python normalization/combine_normalize_phenotypes.py \
    --workers 4 --output data/combined_normalized.tsv \
    --acmg-obo ontology/acmg_criteria.obo

# Smoke test (5 rows):
uv run python normalization/combine_normalize_phenotypes.py \
    --limit 5 --workers 1 --output /tmp/test_combined.tsv \
    --acmg-obo /tmp/test_acmg.obo
```

Source files in `data/phenotypes/`:

| File | Letter | Description |
|---|---|---|
| `ahmed-pmid28454995.tsv` | A/B | HPOs, Gene, Variant, Protein, Zygosity, ACMG |
| `fawzan-variants.tsv` | F | Phenotype, Variant(s), Zyogsity [sic] |
| `marwa-variants.tsv` | M | phenotypes (pipe+HPO), variants (pipe) |
| `PMC6562004.tsv` | P | skiprows=1, Gene(s), Variant(s) |
| `PMC7082194.tsv` | Q | Unnamed:13=phenotype, protein 1-letter AAs |
| `ddd-diagnoses.tsv` | D | gene-disease assoc, hpo_ids (semicolons), allelic_mode |

## Frontend Development

```bash
cd website/frontend
npm install
npm run dev     # Dev server on port 3000, proxies /api/ → localhost:8000
npm run build   # Production build to build/
```

i18n files: `src/i18n/en.json` and `src/i18n/ar.json`. Adding a new UI string requires updating both files.
