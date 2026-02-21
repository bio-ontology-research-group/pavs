# backend_sparql — Change Log

## 2026-02 — SPARQL stack implementation & bug fixes

### Architecture

Replaced the original SQLite/pyhpo backend (`backend/`) with a SPARQL-native backend
(`backend_sparql/`) backed by a Virtuoso triple store.  All persistent data lives in five
named graphs; the Python layer only holds pre-loaded caches for phenotype similarity
computation.

Files:
- `main.py` — FastAPI app; startup caches; all API endpoints
- `sparql_queries.py` — All SPARQL query templates (no inline SPARQL elsewhere)
- `similarity.py` — Lin/Resnik + BMA, pure Python, no pyhpo dependency

### Bug fixes applied

#### Multi-phenotype data loss (Issue #6)
**Symptom**: Only one phenotype shown per case in the UI.
**Root cause**: `case_data["properties"][prop_short] = val` in `get_case()` overwrote
duplicate keys — `hasPhenotype` rows from the SPARQL result were all collapsed to the
last one.
**Fix**: Accumulate `hasPhenotype` and `hasExcludedPhenotype` into separate lists
(`phenotype_uris`, `excluded_phenotype_uris`) before building the properties dict.
Return them as `phenotypes: [{id, label}]` and `excluded_phenotypes: [{id, label}]`.

#### HPO labels missing (Issue #1)
**Symptom**: Phenotypes displayed as bare IDs ("HP:0001263") with no human-readable label.
**Fix**: Added `get_hpo_labels(hpo_uris)` SPARQL query that fetches `rdfs:label` from
`graph/hpo-ic`.  Called after collecting phenotype URIs; result stored in `label_map`.
`compute_hpo_ic.py` was extended to emit `rdfs:label` triples during TTL generation.

#### isSaudi always false (Issue from earlier session)
**Symptom**: Saudi cases not recognised; ancestry badge not shown.
**Root cause**: Virtuoso encodes `true^^xsd:boolean` as the string `"1"` in SPARQL JSON
responses, not `"true"`.
**Fix**: Normalise with `if is_saudi in ("1", "true", "True")` at the end of startup
loading and in `get_case()`.

#### progressStatus not human-readable (Issue #4)
**Symptom**: Status shown as "IN_PROGRESS" (Phenopackets v2 enum value).
**Fix**: Added `_STATUS_MAP` dict; maps SOLVED→"Solved", UNSOLVED→"Unsolved",
IN_PROGRESS→"Unknown", UNKNOWN_PROGRESS→"Unknown".

#### Score displayed as percentage (Issue #3)
**Symptom**: Score "118.4%" for Lin similarity (should be ≤ 1.0; Resnik is unbounded).
**Fix**: Frontend `PhenotypeSearch.tsx` was showing `(score * 100).toFixed(1) + "%"`.
Changed to `score.toFixed(4)` (raw decimal).  Lin BMA is nominally [0,1] but can
fractionally exceed 1.0 due to IC value inconsistencies between ic_cache and
ancestor_cache when loaded from the same hpo_ic.ttl.

#### TogoVar links broken (Issue #5 — two-stage fix)
**Stage 1 symptom**: Links using wrong domain (`togovar.org` instead of
`grch38.togovar.org`) and wrong URL format (`?query=` GET parameter).
**Stage 1 fix**: Switched to `grch38.togovar.org/search?term=chrom:pos`.  Still broken
because TogoVar's `/search` endpoint returns 501 for all GET requests.

**Stage 2 diagnosis**: TogoVar's API requires
`POST https://grch38.togovar.org/api/search/variant` with JSON body
`{"query":{"location":{"chromosome":"16","position":48258198}}}`.
Chromosome must be a bare string number ("16"), not "chr16".

**Stage 2 fix**: Added `/api/togovar-search?chrom=X&pos=Y` backend proxy endpoint.
- Strips "chr" prefix from chromosome
- POSTs to `https://grch38.togovar.org/api/search/variant`
- Extracts `data[0].id` (e.g. `tgv47264307`)
- Returns `302` redirect to `https://grch38.togovar.org/variant/{tgv_id}`
- Returns `302` to `https://grch38.togovar.org/` (homepage) if variant not found
  (expected for Saudi-specific variants absent from Japanese cohorts)

`_togovar_url()` now returns `/api/togovar-search?chrom=X&pos=Y`, which is served
through nginx (`/api/` → backend:8000) in both production Docker and Vite dev proxy.

### New features

#### Suggested OMIM diseases
For Saudi cases without a confirmed `diseaseLabel`:
- Looks up related OMIM diseases for the causative gene via `graph/genes`
- Only shown if the gene has 1–3 associated diseases (≥4 → too ambiguous, suppressed)
- Labels fetched from `disease_label_cache` (loaded from `data/phenotype.hpoa` at
  startup, 12,958 disease names)
- Response field: `suggested_diseases: [{id, label, url, gene}]`
- Frontend: dashed amber border, italic text, "⚠ Suggested (not a diagnosis)" warning

### Turtle / Virtuoso ingestion fixes (pipeline)

- `generate_rdf.py`: removed all `^^xsd:boolean` and `^^xsd:type` annotations from
  Turtle literals (bare `true`/`false` and numeric literals only)
- `compute_hpo_ic.py`: added `rdfs:label` triples from `name:` fields in `hp.obo`
- `docker-compose-sparql.yml`: added `loader` service using `isql` + `ld_dir +
  rdf_loader_run` to bypass the 10 MB HTTP SPARQL endpoint limit
- `virtuoso/load_ttl.sql`: isql script clears old graphs, calls `ld_dir` for each TTL,
  then `rdf_loader_run()` + `checkpoint`
- `VIRT_Parameters_DirsAllowed`: must include `/rdf_import` to allow bulk loading from
  the mounted volume
