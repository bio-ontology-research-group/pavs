# PAVS — Phenotypic and Variant Standardization

A curated database of genomic variants and phenotypes from Saudi Arabian patients with rare genetic diseases. Data integrates multiple clinical cohorts and links each case to international ontologies (HPO, OMIM, MONDO) and variant databases (ClinVar, dbSNP, TogoVar).

---

## Contributors

| Name | ORCID | Role |
|---|---|---|
| **Marwa Abdelhakim** | [0000-0001-6816-2119](https://orcid.org/0000-0001-6816-2119) | Manual curation of all source datasets |
| **Azza Althagafi** | [0000-0001-6084-8706](https://orcid.org/0000-0001-6084-8706) | Computational experiments and validation |
| **Paul N Schofield** | [0000-0002-5111-7263](https://orcid.org/0000-0002-5111-7263) | Curation, supervision, and guidance |
| **Robert Hoehndorf** | [0000-0001-8149-5890](https://orcid.org/0000-0001-8149-5890) | Supervision, code contributions, validation analysis, website |

---

## Quick start

```bash
# Start Virtuoso + loader + API + frontend
docker compose -f docker-compose-sparql.yml up -d

# Rebuild after code changes to backend_sparql/ or frontend_sparql/ only
docker compose -f docker-compose-sparql.yml build backend frontend
docker compose -f docker-compose-sparql.yml up -d backend frontend
```

| Service | URL |
|---|---|
| Frontend | http://localhost:3000 |
| Backend API + docs | http://localhost:8000/docs |
| Virtuoso SPARQL UI | http://localhost:8890/sparql |

---

## Architecture

```
Virtuoso (port 8890)          backend_sparql/ (port 8000)     frontend_sparql/ (port 3000)
  graph/cases     ──SPARQL──▶  FastAPI                  ──▶  React + TypeScript
  graph/genes                    similarity.py (Lin/BMA)        5 tabs:
  graph/hpoa                     sparql_queries.py              • Phenotype Search
  graph/hpo-ic                   /api/togovar-search proxy      • Variant Lookup
  graph/literature                                              • Disease Browser
                                                               • Gene Browser
                                                               • About
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

Cases use the format `PAVS:XNNNNNNN` where X is the source letter:

| Letter | Source dataset |
|---|---|
| A | ahmed-variants |
| B | ahmed-pmid28454995 |
| F | fawzan-variants |
| M | marwa-variants |
| P | PMC6562004 |
| Q | PMC7082194 |
| D | ddd-diagnoses (non-Saudi) |

### Backend startup sequence

At startup the backend loads four in-memory caches from Virtuoso and from disk:

1. `ic_cache` — HPO term → IC float (12,725 terms) from `graph/hpo-ic`
2. `ancestor_cache` — HPO term → set of ancestor HP IDs, transitive closure (19,408 terms)
3. `case_hpo_cache` — all cases with HPO sets (17,149 cases) for in-memory similarity scoring
4. `disease_label_cache` — `"OMIM:272200" → "Multiple sulfatase deficiency"` from `data/phenotype.hpoa` (12,958 entries)

Startup takes ~15 seconds.

### Phenotype similarity

Lin + BMA (Best Match Average) by default; Resnik + BMA available. Scores are **not** percentages — Lin values are in [0, 1]; Resnik values are unbounded IC units. Both query and target HPO sets are expanded with ancestor terms before scoring.

### TogoVar integration

The `/api/togovar-search?chrom=X&pos=Y` endpoint proxies to `POST https://grch38.togovar.org/api/search/variant`. If found it redirects to the TogoVar variant page. If not found (e.g. Saudi-specific variants absent from the Japanese cohort database) it redirects to the TogoVar homepage.

### Suggested diseases

For Saudi cases without a confirmed disease, if the causative gene has 1–3 associated OMIM diseases, a clearly-marked suggestion is shown in the case detail view (dashed amber box). Not shown when the gene has ≥4 diseases (too ambiguous).

---

## API endpoints

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
  "properties": { "sex": "...", "consanguinity": "yes", "progressStatus": "Solved", "isSaudi": "true" },
  "phenotypes": [{ "id": "HP:0001263", "label": "Intellectual disability" }],
  "excluded_phenotypes": [],
  "suggested_diseases": [],
  "variants": [{ "gene": "SUMF1", "hgvsC": "...", "togovar_url": "/api/togovar-search?chrom=3&pos=4417183" }],
  "phenopacket_download_url": "/api/phenopacket/PAVS:A0000001/download"
}
```

---

## Data processing pipeline

The full pipeline converts raw source TSVs to Phenopackets v2.0.2 JSON and then to RDF for Virtuoso. All scripts are run with `uv`.

```
data/phenotypes/  (7 source TSVs)
      │
      ▼  normalization/combine_normalize_phenotypes.py
      │
data/combined_normalized.tsv  (38 columns, ~7,500 rows)
      │
      ▼  normalization/annotate_variants.py
      │
data/combined_annotated.tsv  (77+ columns)
      │
      ▼  intake/generate_phenopackets_v2.py
      │
      ├─▶ phenopackets/generated_v2/  (individual JSON per case)
      ├─▶ data/PAVS_phenopackets.json (combined array)
      └─▶ data/PAVS_phenopackets.zip  (download bundle)
            │
            ▼  intake/compute_hpo_ic.py + intake/generate_rdf.py
            │
      rdf_output/*.ttl  →  Virtuoso triple store  →  website
```

### Prerequisites

```bash
uv sync
export OPENROUTER_API_KEY="your_key_here"

# Install the phenotype normalization toolkit
uv pip install -e tools/phenotype_matcher_v2

# Download all external reference data (to data/reference/):
bash scripts/download_reference_data.sh
```

Required ontology files (not in git):

| Path | Description | Source |
|---|---|---|
| `ontology/hp.obo` | Human Phenotype Ontology | https://hpo.jax.org/app/data/ontology |
| `ontology/mondo.obo` | MONDO Disease Ontology | https://mondo.monarchinitiative.org/ |

Required reference data (in `data/reference/`, downloaded by script):

| File | Description |
|---|---|
| `phenotype.hpoa` | HPO disease annotations (~35 MB) |
| `genes_to_disease.txt` | Gene→disease associations |
| `gnomad.v4.1.constraint_metrics.tsv` | pLI / LOEUF constraint scores |
| `clinvar.vcf.gz` | ClinVar variant database (GRCh38) |
| `erepo.tabbed.txt` | ClinGen expert-curated variants |
| `goa_human.gaf.gz` | GO annotations for human genes |
| `Homo_sapiens.gene_info` | NCBI gene info |
| `HMD_HumanPhenotype.rpt` | Mouse phenotype homologs (MGI) |
| `MGI_GenePheno.rpt` | MGI genotype-phenotype associations |
| `E-GTEX-8/` | GTEx v8 tissue expression (EBI) |
| `allele-frequencies/` | Saudi population AF — **not public**, contact authors |

Other required data:

| Path | Description |
|---|---|
| `phenopackets/0.1.26/` | Literature phenopackets (9,588 non-Saudi cases) — contact authors |

### Step 1 — Normalize source TSVs

```bash
uv run python normalization/combine_normalize_phenotypes.py \
    --workers 4 \
    --output data/combined_normalized.tsv \
    --acmg-obo ontology/acmg_criteria.obo

# Smoke test (5 rows only):
uv run python normalization/combine_normalize_phenotypes.py \
    --limit 5 --workers 1 \
    --output /tmp/test_combined.tsv \
    --acmg-obo /tmp/test_acmg.obo
```

Reads all 7 source files in `data/phenotypes/`, normalizes each field to ontology identifiers using `phenotype_matcher_v2`, and writes `data/combined_normalized.tsv` (38 columns, PAVS:XNNNNNNN IDs, pipe-delimited multi-value fields). Also generates `ontology/acmg_criteria.obo`.

Source files (all in `data/phenotypes/`):

| File | Source letter | Rows |
|---|---|---|
| `ahmed-variants.tsv` | A | 291 |
| `ahmed-pmid28454995.tsv` | B | 234 |
| `fawzan-variants.tsv` | F | 1,024 |
| `marwa-variants.tsv` | M | 1,421 |
| `PMC6562004.tsv` | P | 2,218 |
| `PMC7082194.tsv` | Q | 522 |
| `ddd-diagnoses.tsv` | D | 1,856 |

Each source has a dedicated normalizer in `normalization/normalize_*.py`; see **`docs/PIPELINE.md`** for the full field-by-field schema.

### Step 2 — Annotate variants and genes

```bash
uv run python normalization/annotate_variants.py \
    --input  data/combined_normalized.tsv \
    --output data/combined_annotated.tsv

# Quick smoke test (10 rows, skip VEP REST calls):
uv run python normalization/annotate_variants.py \
    --input  data/combined_normalized.tsv \
    --output /tmp/test_annotated.tsv \
    --limit 10 --no-vep
```

Enriches each row with variant-level and gene-level annotations from multiple sources:

**Variant-level** (adds ~17 columns):
- Ensembl VEP REST — consequence, impact, HGVSc/p, SIFT, PolyPhen, gnomAD AFs, rsID, domains
- `data/clinvar.vcf.gz` — clinical significance, disease, ClinVar allele ID
- `data/erepo.tabbed.txt` — ClinGen expert curation (assertion, ACMG codes)
- `data/allele-frequencies/` — Saudi population AF (AC/AN from 302 WGS individuals)

**Gene-level** (adds ~10 columns):
- `data/genes_to_disease.txt` — OMIM/disease associations
- `data/gnomad.v4.1.constraint_metrics.tsv` — pLI, LOEUF
- `data/HMD_HumanPhenotype.rpt` + `data/MGI_GenePheno.rpt` + `ontology/mp.obo` — mouse phenotypes
- `data/goa_human.gaf.gz` + `ontology/go-basic.obo` — GO biological process / MF / CC
- `data/E-GTEX-8/` — tissue expression (GTEx v8, percentile rank)

Output: `data/combined_annotated.tsv` (~77 columns) — the authoritative input for all downstream steps.

### Step 3 — Generate Phenopackets v2

```bash
uv run python intake/generate_phenopackets_v2.py \
    --input data/combined_annotated.tsv \
    --output-dir phenopackets/generated_v2 \
    --combined data/PAVS_phenopackets.json \
    --zip data/PAVS_phenopackets.zip
```

Outputs `phenopackets/generated_v2/` (one JSON per case), `data/PAVS_phenopackets.json` (combined array), and `data/PAVS_phenopackets.zip` (download bundle).

### Step 4 — Compute HPO information content

```bash
uv run python intake/compute_hpo_ic.py \
    --hpo ontology/hp.obo \
    --hpoa data/phenotype.hpoa \
    --output rdf_output/hpo_ic.ttl
```

Parses `hp.obo` for hierarchy, propagates disease annotations through ancestors, computes `IC(t) = -log2(freq(t) / total_diseases)`, emits `rdfs:label` triples needed for HPO autocomplete. Output: `rdf_output/hpo_ic.ttl`.

### Step 5 — Generate all RDF

```bash
uv run python intake/generate_rdf.py \
    --phenopackets data/PAVS_phenopackets.json \
    --annotated data/combined_annotated.tsv \
    --literature-dir phenopackets/0.1.26 \
    --output-dir rdf_output/
```

Outputs: `rdf_output/{cases,genes,hpoa,literature}.ttl`.

### Step 6 — Load into Virtuoso (handled by Docker)

The `loader` service in `docker-compose-sparql.yml` runs `isql` with `virtuoso/load_ttl.sql`, which uses Virtuoso's `ld_dir + rdf_loader_run` bulk loader (no 10 MB file-size limit). Named graphs created: `graph/cases`, `graph/genes`, `graph/hpoa`, `graph/hpo-ic`, `graph/literature`.

### Regenerating after data changes

```bash
docker compose -f docker-compose-sparql.yml run --rm init --force
docker compose -f docker-compose-sparql.yml restart loader
```

---

## Normalization toolkit

`tools/phenotype_matcher_v2/` is the standalone, reusable phenotype matching library used throughout the pipeline.

```bash
# Install
uv pip install -e tools/phenotype_matcher_v2

# Python API
from phenotype_matcher_v2 import PhenotypeMatcher
matcher = PhenotypeMatcher()
results = matcher.match("severe intellectual disability and seizures")

# CLI
phenotype-matcher "severe intellectual disability and seizures"
phenotype-matcher "seizures, hypotonia" --output-format tsv
```

The library maps free-text clinical descriptions to HPO, MONDO, OMIM, and Orphanet identifiers, detecting negation ("no seizures") and severity modifiers ("profound"). Algorithm details: `tools/phenotype_matcher_v2/ALGORITHM.md`.

---

## Source-specific normalizers

Each source dataset has a dedicated normalizer in `normalization/`:

| Script | Source | Key handling |
|---|---|---|
| `normalization/normalize_ahmed.py` | ahmed-variants.tsv (A) | Complex phenotypic sentences, ACMG classification |
| `normalization/normalize_ahmed_pmid28454995.py` | ahmed-pmid28454995.tsv (B) | Merges Diagnosis + Additional clinical phenotype columns |
| `normalization/normalize_fawzan.py` | fawzan-variants.tsv (F) | Variant nomenclature (HGVS), HGMD field parsing |
| `normalization/normalize_marwa.py` | marwa-variants.tsv (M) | Disentangles mixed phenotype and disease descriptions |
| `normalization/normalize_pmc6562004.py` | PMC6562004.tsv (P) | skiprows=1 (comment line), multi-gene/variant columns |
| `normalization/normalize_pmc7082194.py` | PMC7082194.tsv (Q) | Non-standard phenotype column (Unnamed:13), 1-letter AAs |
| `normalization/normalize_ddd.py` | ddd-diagnoses.tsv (D) | Gene-disease associations, semicolon-delimited HPO IDs |

All normalizers share `normalization/normalization_utils.py` — a singleton `NormalizationUtils` class that loads HPO, MONDO, GENO, NCIT, and OMIM at initialisation and provides `normalize_sex()`, `normalize_genotype()`, `normalize_consanguinity()`, and related helpers.

---

## Key conventions

- **Case IDs**: `PAVS:XNNNNNNN` (X = source letter, 7-digit zero-padded sequence)
- **Genome build**: GRCh38 for VCF coordinates; GRCh37 for legacy HGVS strings
- **Population tags**: Saudi cases tagged with `HANCESTRO:0852` (Middle Eastern) and `GAZ:00005279` (Saudi Arabia)
- **Multi-value fields**: pipe-delimited (`|`) in all TSV outputs
- **Provenance**: `source_file` and `source_id` columns preserved in all derived datasets

---

## Virtuoso gotchas

- **Boolean encoding**: SPARQL JSON returns `"1"` not `"true"` for `true^^xsd:boolean`. Backend normalises with `in ("1", "true")`.
- **Turtle syntax**: Bare `true`/`false` and numeric literals cannot have `^^xsd:type` annotations in Turtle.
- **10 MB HTTP limit**: Virtuoso's SPARQL HTTP endpoint rejects payloads > 10 MB. Use `isql` + `ld_dir + rdf_loader_run` for large TTL files.
- **DirsAllowed**: The `/rdf_import` mount must appear in `VIRT_Parameters_DirsAllowed` env var.

---

## Full pipeline documentation

**`docs/PIPELINE.md`** — complete source schemas, output column reference (38 + 77 columns), normalization rules, phenopacket structure, ACMG OBO generation.

**`tools/phenotype_matcher_v2/ALGORITHM.md`** — detailed description of the 6-stage phenotype matching algorithm.

**`docs/NORMALIZATION_METHODS.md`** — algorithm comparison and validation results.

**`docs/NORMALIZATION.md`** — normalization pipeline execution guide.
