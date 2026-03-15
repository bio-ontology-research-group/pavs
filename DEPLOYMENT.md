# PAVS Deployment History

## 2026-03-14: Complete KG Rebuild and IRI Migration

### Overview
Complete rebuild of the PAVS knowledge graph with migration of all IRIs from `http://pavs.kaust.edu.sa/` to `https://pavs.phenomebrowser.net/`.

### Changes Made

#### 1. IRI Migration
Migrated all resource identifiers from HTTP to HTTPS with actual domain:
- **Old**: `http://pavs.kaust.edu.sa/`
- **New**: `https://pavs.phenomebrowser.net/`

**Affected Components:**
- Named graph URIs (cases, genes, hpoa, hpo-ic, literature, metadata)
- Resource IRIs (cases, variants, genes)
- Ontology namespace (pavs:, pav:)
- Backend SPARQL query prefixes
- Metadata dataset identifier

#### 2. RDF Data Regeneration
Generated fresh RDF files with new IRIs:
- `cases.ttl` - 175,969 triples (7,510 phenopackets: 5,654 Saudi)
- `genes.ttl` - 147,331 triples (2,523 genes)
- `hpoa.ttl` - 1,613,546 triples (disease-phenotype associations)
- `hpo_ic.ttl` - 56,424 triples (IC values + hierarchy + 19,934 labels)
- `literature.ttl` - 407,246 triples (9,588 literature phenopackets)
- `metadata.ttl` - 163 triples (FAIR-compliant metadata)

#### 3. Metadata Enhancement
Added comprehensive FAIR-compliant metadata including:
- VoID (Vocabulary of Interlinked Datasets)
- DCAT (Data Catalog Vocabulary)
- Dublin Core terms
- PROV-O provenance
- PAV versioning
- Schema.org markup
- CC-BY 4.0 license
- Version: 1.0.0

#### 4. Code Updates
**knowledge-graph/intake/**
- `compute_hpo_ic.py` - Updated namespace prefixes
- `generate_rdf.py` - Updated namespace prefixes
- `generate_metadata.py` - Updated base URI and website URLs
- `load_virtuoso.py` - Updated graph URIs

**website/backend/**
- `sparql_queries.py` - Updated graph constants and SPARQL prefixes

**website/virtuoso/**
- `load_ttl.sql` - Updated graph URIs and added metadata.ttl

#### 5. Deployment
**Server**: onto.phenomecentral.org at `/data/pavs/`

**Process:**
1. Stopped all containers
2. Removed Virtuoso volume for clean state
3. Rebuilt backend container with updated code
4. Started fresh stack with new RDF data

**Services:**
- Frontend: `http://onto:20000` → `https://pavs.phenomebrowser.net/`
- Backend API: `http://onto:20001`
- Virtuoso SPARQL: `http://onto:20002` → `https://pavs.phenomebrowser.net/sparql`
- Virtuoso SQL: `http://onto:20003`

### Verification Results

#### Backend Initialization
- ✅ 12,725 HPO IC values loaded
- ✅ 19,408 HPO ancestor sets built
- ✅ 17,201 cases cached for similarity search
- ✅ 12,958 disease labels loaded
- ✅ 19,408 Arabic translations loaded

#### API Endpoints Tested
- ✅ `/api/health` - All services operational
- ✅ `/api/case/{id}` - Returns correct IRIs
- ✅ `/api/search/hpo` - HPO autocomplete working
- ✅ `/api/search/phenotype` - Similarity search functioning
- ✅ `/sparql` - Public SPARQL endpoint accessible

#### IRI Validation
**Example case (PAVS:A0000001):**
- Case IRI: `https://pavs.phenomebrowser.net/data/PAVS_A0000001`
- Type: `https://pavs.phenomebrowser.net/ontology/Case`
- Variant: `https://pavs.phenomebrowser.net/data/var_PAVS_A0000001_SUMF1`

**Named Graphs:**
```sparql
https://pavs.phenomebrowser.net/graph/cases       (175,969 triples)
https://pavs.phenomebrowser.net/graph/genes       (147,331 triples)
https://pavs.phenomebrowser.net/graph/hpo-ic      (56,424 triples)
https://pavs.phenomebrowser.net/graph/hpoa        (1,613,546 triples)
https://pavs.phenomebrowser.net/graph/literature  (407,246 triples)
https://pavs.phenomebrowser.net/graph/metadata    (163 triples)
```

### Statistics

**Data Loaded:**
- Total triples: ~2.4 million
- Cases: 17,201 total
  - Saudi cases: 5,654
  - Literature cases: 9,588
  - Other cases: 1,959
- Genes: 2,523
- HPO terms with IC: 12,725
- Disease-phenotype associations: 280,003

**Coverage:**
- HPO ancestor terms: 19,408
- Saudi-specific term counts: 2,855
- Diseases in HPOA: 12,958

### Access Points

**Public URLs:**
- Website: https://pavs.phenomebrowser.net/
- API: https://pavs.phenomebrowser.net/api/
- SPARQL: https://pavs.phenomebrowser.net/sparql
- API Docs: https://pavs.phenomebrowser.net/docs

**SPARQL Example Queries:**

List all graphs:
```sparql
SELECT ?g (COUNT(*) AS ?triples) 
WHERE { GRAPH ?g { ?s ?p ?o } } 
GROUP BY ?g ORDER BY ?g
```

Get dataset metadata:
```sparql
PREFIX dct: <http://purl.org/dc/terms/>
SELECT ?title ?version ?license
WHERE {
  GRAPH <https://pavs.phenomebrowser.net/graph/metadata> {
    <https://pavs.phenomebrowser.net/dataset>
      dct:title ?title ;
      dct:hasVersion ?version ;
      dct:license ?license
  }
}
```

Query case by identifier:
```sparql
PREFIX dc: <http://purl.org/dc/terms/>
PREFIX pavs: <https://pavs.phenomebrowser.net/ontology/>
SELECT ?case ?type ?variant
WHERE {
  GRAPH <https://pavs.phenomebrowser.net/graph/cases> {
    ?case dc:identifier "PAVS:A0000001" ;
          a ?type ;
          pavs:hasVariant ?variant
  }
}
```

### Notes

- Previous deployment backed up at `/data/pavs_backup_20260314_095915/` on onto server
- Old IRI graphs (`http://pavs.kaust.edu.sa/graph/*`) were removed from Virtuoso
- No changes made to nginx configuration (routing remains intact)
- All port mappings unchanged from previous deployment
- Local development environment also updated with new IRIs

## 2026-03-14: Security Dependency Updates

### Overview
Patched high-severity vulnerabilities in `pillow` and `protobuf` identified in March 2026.

### Vulnerabilities Fixed
- **Pillow (CVE-2026-25990)**: Out-of-bounds write when loading PSD images.
- **Protobuf (CVE-2026-0994)**: JSON recursion depth bypass and potential Denial of Service.

### Changes Made
- **pyproject.toml**: Added `tool.uv.override-dependencies` for `protobuf>=6.33.5` to bypass restrictive version pins in transitive dependencies (`pyphetools`).
- **uv.lock**: 
    - Updated `pillow` from `12.1.0` to `12.1.1`.
    - Updated `protobuf` from `3.20.3` to `6.33.5`.
- **Environment**: Synchronized local venv using `uv sync`.

### Related Issues

- Fixed metadata.ttl Turtle syntax error (missing angle brackets around dataset URI)
- Updated docker-compose volume mounts to use main directory paths
- Rebuilt backend container to pick up new graph URIs

## 2026-03-15: Website Improvements and Bug Fixes

### Overview
Improved phenotype search functionality, cohort filtering, and diagnosis suggestions.

### Changes Made

#### 1. Phenotype Search Enhancements
- **"Show Only Diagnosed Cases" Filter**: Added a new filter to exclude cases without a confirmed causative gene.
- **Removed Redundant Options**:
    - Removed "ClinVar" cohort from search (no cases currently loaded in production VCF).
    - Removed "Expand with OMIM/ClinVar disease HPO annotations" checkbox as requested.
- **ACMG Visibility & Translation Fix**:
    - Fixed an issue where Arabic labels in the ACMG class dropdown were invisible due to CSS rendering quirks.
    - Ensured consistent translation of ACMG classes in search results using the `localizeAcmg` helper.
- **Back Button Support**: Refactored Variant Lookup tabs and HPO Browser selection to be URL-driven. Enhanced HPO Browser to also track tree expansions in the URL, ensuring the browser "Back" button correctly restores the tree state (expanded/collapsed nodes).

#### 2. Diagnosis Logic Improvements
- **Suggested Diagnoses**: Refined logic to only show gene-based suggested diseases when no explicit clinical diagnosis exists.
- **Global Availability**: Suggestions are now shown for all cohorts (not just Saudi) when the diagnosis is missing.

#### 3. Cohort Filtering Fixes
- **Robust Saudi Identification**: Updated backend to check both `is_saudi` flags and source name strings (e.g., "Saudi cohort").
- **Case Fixes**: Fixed `PAVS:P0000919` being misclassified as Literature/DDD; it now correctly appears in Saudi results.

#### 4. Internationalization
- Added English and Arabic translations for "Show only diagnosed cases".
- Removed translations for decommissioned search options.

### Deployment
- **Files Updated on onto**: 
    - `website/backend/main.py`
    - `website/frontend/src/components/PhenotypeSearch.tsx`
    - `website/frontend/src/i18n/en.json`
    - `website/frontend/src/i18n/ar.json`
- **Containers Rebuilt**: `pavs-backend-1`, `pavs-frontend-1`

### Future Maintenance

To regenerate RDF with updates:
```bash
cd knowledge-graph
uv run python intake/prepare_all.py --force --skip-load \
  --input ../data/combined_annotated-gpt4oss.tsv \
  --hpo ../ontology/hp.obo \
  --hpoa ../data/reference/phenotype.hpoa \
  --literature-dir ../phenopackets/0.1.26

# Generate metadata
uv run python intake/generate_metadata.py \
  --version 1.0.0 \
  --output rdf_output/metadata.ttl

# Fix metadata URI (temporary until script is fixed)
sed -i 's|^https://pavs.phenomebrowser.net/dataset$|<https://pavs.phenomebrowser.net/dataset>|' \
  rdf_output/metadata.ttl

# Copy to main directory
cp rdf_output/*.ttl ../rdf_output/
```

To reload Virtuoso on server:
```bash
ssh onto
cd /data/pavs
docker compose -f docker-compose-sparql.yml down -v  # Clean slate
docker compose -f docker-compose-sparql.yml up -d    # Reload all
```
