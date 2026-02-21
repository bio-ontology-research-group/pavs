# PAVS — Phenotypic and Variant Standardization

A curated database of genomic variants and phenotypes from Saudi Arabian patients with rare genetic diseases.  Data integrates multiple clinical cohorts and links each case to international ontologies (HPO, OMIM, MONDO) and variant databases (ClinVar, dbSNP, TogoVar).

---

## Two stacks

| | **SPARQL stack** (current) | **SQLite stack** (original) |
|---|---|---|
| Storage | Virtuoso triple store | SQLite + phenopacket JSON |
| Search | SPARQL + in-memory HPO similarity | pyhpo + SQLite |
| Compose file | `docker-compose-sparql.yml` | `docker-compose.yml` |
| Backend | `backend_sparql/` | `backend/` |
| Frontend | `frontend_sparql/` | `frontend/` |

---

## SPARQL stack — quick start

```bash
# 1. Generate RDF from source data (only needed when data changes)
uv run python scripts/compute_hpo_ic.py \
    --hpo ontology/hp.obo --hpoa data/phenotype.hpoa \
    --output rdf_output/hpo_ic.ttl

uv run python scripts/generate_rdf.py \
    --phenopackets data/PAVS_phenopackets.json \
    --annotated data/combined_annotated.tsv \
    --literature-dir phenopackets/0.1.26 \
    --output-dir rdf_output/

# 2. Start Virtuoso, load data, start API + frontend
docker compose -f docker-compose-sparql.yml up -d

# Rebuild after code changes to backend_sparql/ or frontend_sparql/
docker compose -f docker-compose-sparql.yml build backend frontend
docker compose -f docker-compose-sparql.yml up -d backend frontend
```

Frontend: http://localhost:3000
Backend API + docs: http://localhost:8000/docs
Virtuoso SPARQL UI: http://localhost:8890/sparql

---

## SPARQL stack — architecture

```
Virtuoso (port 8890)          backend_sparql/ (port 8000)     frontend_sparql/ (port 3000)
  graph/cases     ──SPARQL──▶  FastAPI                  ──▶  React + TypeScript
  graph/genes                    similarity.py (Lin/BMA)        5 tabs:
  graph/hpoa                     sparql_queries.py              • Phenotype Search
  graph/hpo-ic                   /api/togovar-search proxy      • Variant Lookup
  graph/literature                                              • Disease Browser
                                                               • Gene Browser
                                                               • About (Markdown)
```

### Named graphs in Virtuoso

| Graph | Contents | Triples |
|---|---|---|
| `graph/cases` | Saudi PAVS cases (phenopackets, variants, demographics) | 163 K |
| `graph/genes` | Gene→disease, constraint metrics (pLI, LOEUF), GO, GTEx | 82 K |
| `graph/hpoa` | Disease→HPO annotations from `phenotype.hpoa` | 1.6 M |
| `graph/hpo-ic` | Pre-computed HPO IC values + ontology hierarchy + labels | 56 K |
| `graph/literature` | 9,588 non-Saudi literature phenopackets | 407 K |

### Case IDs

Cases use the format `PAVS:XNNNNNNN` where X is the source letter (A=ahmed, F=fawzan, etc.).

### Backend startup sequence

1. Fetch all HPO IC values from Virtuoso → `ic_cache`
2. Fetch HPO parent closure → `ancestor_cache`
3. Fetch all case HPO sets → `case_hpo_cache` (used for in-memory similarity scoring)
4. Load `data/phenotype.hpoa` → `disease_label_cache` (12,958 OMIM disease names)

### Phenotype similarity

Lin + BMA (Best Match Average) by default; Resnik + BMA available.  Scores are **not** percentages — Lin values are in [0, 1]; Resnik values are unbounded IC units.  Both query and target HPO sets are expanded with ancestor terms before scoring.

### TogoVar integration

The `/api/togovar-search?chrom=X&pos=Y` endpoint proxies to `POST https://grch38.togovar.org/api/search/variant`.  If found it redirects to `https://grch38.togovar.org/variant/{tgv_id}`.  If not found (e.g. Saudi-specific variants absent from Japanese cohorts) it redirects to the TogoVar homepage.

### Suggested diseases

For Saudi cases without a confirmed disease, if the causative gene has 1–3 associated OMIM diseases, a clearly-marked suggestion is shown in the case detail view (dashed amber box, italic text). Not shown when the gene has ≥4 diseases (too ambiguous).

---

## SPARQL stack — data pipeline

### Prerequisites

```bash
uv sync
# Required data files (not in git):
#   ontology/hp.obo          — HPO ontology
#   data/phenotype.hpoa      — HPO disease annotations (~35 MB)
#   data/combined_annotated.tsv  — processed PAVS cases
#   data/genes_to_disease.txt    — gene→disease associations
#   phenopackets/0.1.26/     — literature phenopackets (non-Saudi)
```

### Step 1 — Compute HPO information content

```bash
uv run python scripts/compute_hpo_ic.py \
    --hpo ontology/hp.obo \
    --hpoa data/phenotype.hpoa \
    --output rdf_output/hpo_ic.ttl
```

Parses `hp.obo` for hierarchy, propagates disease annotations upward through ancestor terms, computes `IC(t) = -log2(freq(t) / total_diseases)`, emits `rdfs:label` triples for HPO autocomplete.

### Step 2 — Generate all RDF

```bash
uv run python scripts/generate_rdf.py \
    --phenopackets data/PAVS_phenopackets.json \
    --annotated data/combined_annotated.tsv \
    --literature-dir phenopackets/0.1.26 \
    --output-dir rdf_output/
```

Outputs: `rdf_output/{cases,genes,hpoa,literature}.ttl`

### Step 3 — Load into Virtuoso (handled by Docker)

The `loader` service in `docker-compose-sparql.yml` runs `isql` with `virtuoso/load_ttl.sql`, which uses Virtuoso's `ld_dir + rdf_loader_run` bulk loader (no file-size limit, unlike the HTTP SPARQL endpoint).

### Regenerating after data changes

```bash
# Regenerate TTL files, then reload Virtuoso:
docker compose -f docker-compose-sparql.yml run --rm init --force
docker compose -f docker-compose-sparql.yml restart loader
```

---

## SPARQL stack — API endpoints

| Method | Endpoint | Description |
|---|---|---|
| GET | `/api/health` | IC terms, ancestor terms, case cache size |
| GET | `/api/search/hpo?q=` | HPO autocomplete (label substring match) |
| POST | `/api/search/phenotype` | Rank cases by HPO similarity (Lin+BMA default) |
| GET | `/api/search/gene?q=` | Search cases by gene symbol |
| GET | `/api/search/variant` | Lookup by gene / rsID / HGVS / ACMG class |
| GET | `/api/search/disease?q=` | Search by disease label substring |
| GET | `/api/case/{id}` | Full case detail: demographics, phenotypes+labels, variants, suggested diseases |
| GET | `/api/gene/{symbol}` | Gene detail from genes graph |
| GET | `/api/togovar-search?chrom=&pos=` | Proxy to TogoVar API → redirect to variant page |
| GET | `/api/phenopacket/{id}/download` | Download individual phenopacket JSON |
| GET | `/api/phenopackets/download-all` | Download `PAVS_phenopackets.zip` |

### POST `/api/search/phenotype` request body

```json
{
  "hpo_ids": ["HP:0001263", "HP:0001250"],
  "method": "lin",
  "limit": 20,
  "include_disease_phenotypes": false,
  "include_non_saudi": false
}
```

### GET `/api/case/{id}` response shape

```json
{
  "id": "PAVS:A0000001",
  "properties": { "sex": "...", "consanguinity": "yes", "progressStatus": "Solved", "isSaudi": "true", ... },
  "phenotypes": [{ "id": "HP:0001263", "label": "Intellectual disability" }],
  "excluded_phenotypes": [],
  "suggested_diseases": [],
  "variants": [{ "gene": "SUMF1", "hgvsC": "...", "togovar_url": "/api/togovar-search?chrom=3&pos=4417183", ... }],
  "phenopacket_download_url": "/api/phenopacket/PAVS:A0000001/download"
}
```

`suggested_diseases` is populated only when `diseaseLabel` is absent and the gene has 1–3 OMIM associations:

```json
"suggested_diseases": [{
  "id": "OMIM:618792",
  "label": "Epileptic encephalopathy, early infantile, 84",
  "url": "https://omim.org/entry/618792",
  "gene": "UGDH"
}]
```

---

## Original SQLite stack

```bash
# Start (port 8000 + 3000)
uv run manage.py start
uv run manage.py status / stop / logs / restart

# Docker (port 8000)
docker compose up --build
```

See `backend/` and `frontend/` for source. API docs at `http://localhost:8000/docs`.

---

## Data processing pipeline (upstream of both stacks)

```bash
uv sync
export OPENROUTER_API_KEY="your_key_here"

# Normalize source TSVs → HPO/OMIM/gene-normalized TSVs
uv run python scripts/combine_normalize_phenotypes.py \
    --workers 4 --output data/combined_normalized.tsv \
    --acmg-obo ontology/acmg_criteria.obo

# Smoke test (5 rows):
uv run python scripts/combine_normalize_phenotypes.py \
    --limit 5 --workers 1 --output /tmp/test_combined.tsv \
    --acmg-obo /tmp/test_acmg.obo
```

Source files in `data/phenotypes/`: ahmed-variants.tsv (A), ahmed-pmid28454995.tsv (B), fawzan-variants.tsv (F), marwa-variants.tsv (M), PMC6562004.tsv (P), PMC7082194.tsv (Q), ddd-diagnoses.tsv (D).

Output: `data/combined_normalized.tsv` (38 columns, PAVS:XNNNNNNN IDs, pipe-delimited multi-value fields).

Full pipeline documentation (source schemas, output columns, normalization rules, phenopacket structure, `phenotype_matcher_v2` algorithms): **`scripts/PIPELINE.md`**

Normalization toolkit documentation: **`tools/phenotype_matcher_v2/README.md`**
