# PAVS Data Generation Pipeline

This document describes the full pipeline from raw source TSVs to Phenopackets v2.0 JSON and the Virtuoso RDF triple store.

```
data/phenotypes/         (7 raw TSVs)
       │
       ▼
combine_normalize_phenotypes.py
       │
       ▼ data/combined_normalized.tsv  (38 columns, ~7 500 rows)
       │
       │  (manual annotation / VEP annotation)
       │
       ▼ data/combined_annotated.tsv   (77+ columns)
       │
       ▼
generate_phenopackets_v2.py
       │
       ├─▶ phenopackets/generated_v2/PAVS_MASTER_NNNN.json  (individual)
       ├─▶ data/PAVS_phenopackets.json                       (combined array)
       └─▶ data/PAVS_phenopackets.zip                        (download bundle)
               │
               ▼
         generate_rdf.py + compute_hpo_ic.py
               │
               ▼ rdf_output/*.ttl  → Virtuoso triple store
```

---

## Step 0 — Install and configure

```bash
# From repo root
uv sync
export OPENROUTER_API_KEY="your_key_here"

# Install the normalization toolkit
uv pip install -e tools/phenotype_matcher_v2

# Verify environment
uv run phenotype-matcher check
```

Required ontology files (not in git):

| Path | Source |
|------|--------|
| `ontology/hp.obo` | HPO: https://hpo.jax.org/app/data/ontology |
| `ontology/mondo.obo` | MONDO: https://mondo.monarchinitiative.org/ |
| `data/phenotype.hpoa` | HPO annotations: https://hpo.jax.org/app/data/annotations |

---

## Step 1 — Normalize source TSVs (`combine_normalize_phenotypes.py`)

Reads all 7 source files in `data/phenotypes/`, normalizes each field to standard
ontology identifiers, and writes a single combined TSV.

```bash
uv run python normalization/combine_normalize_phenotypes.py \
    --workers 4 \
    --output data/combined_normalized.tsv \
    --acmg-obo ontology/acmg_criteria.obo

# Smoke test (5 rows only)
uv run python normalization/combine_normalize_phenotypes.py \
    --limit 5 --workers 1 \
    --output /tmp/test_combined.tsv \
    --acmg-obo /tmp/test_acmg.obo
```

### Source files

| File | Source letter | Rows | Notes |
|------|--------------|------|-------|
| `ahmed*` | A/B | 525 | HPOs, Gene, Variant, Protein, Zygosity, ACMG Classification |
| `fawzan-variants.tsv` | F | 1,024 | Phenotype, Variant(s), Zyogsity (typo in source) |
| `marwa-variants.tsv` | M | 1,421 | phenotypes (pipe+HPO), variants (pipe) |
| `PMC6562004.tsv` | P | 2,218 | skiprows=1 (comment line), Gene(s), Variant(s) |
| `PMC7082194.tsv` | Q | 522 | `Unnamed:13` = phenotype, protein 1-letter AAs |
| `ddd-diagnoses.tsv` | D | 1,856 | gene-disease association, hpo_ids (semicolons), allelic_mode |

All Saudi cohorts (A, B, F, M, P, Q) share the same pipeline. D is the DDD non-Saudi cohort.

### PAVS IDs

Each output row gets a globally unique `PAVS:XNNNNNNN` identifier where X is the source letter and NNNNNNN is a zero-padded 7-digit sequence number within that source.

### Normalization performed

| Field | Method |
|-------|--------|
| Sex | Regex → NCIT:C20197 (Male) / NCIT:C16576 (Female) |
| Consanguinity | Keyword regex → yes/no/unknown |
| Inheritance | Regex mapping to HPO inheritance terms (HP:0000006 AD, HP:0000007 AR, etc.) |
| Zygosity | Regex → GENO:0000136 homozygous / GENO:0000135 heterozygous / GENO:0000134 hemizygous / GENO:0000402 compound het |
| Phenotypes | `phenotype_matcher_v2` (see below) → HPO IDs; negated terms → excluded |
| HGVS c. / g. | `preprocess_hgvs()` normalisation (Unicode minus, operator spacing) + regex validation |
| 1-letter amino acids | Converted to 3-letter HGVS p. notation |
| ACMG classification | Keyword matching → Pathogenic / Likely pathogenic / VUS / Likely benign / Benign |
| Disease (OMIM/MONDO) | Via matcher or direct column mapping |

### Output schema (38 columns)

All multi-value fields are **pipe-delimited** (`|`).

| Column | Content |
|--------|---------|
| `pavs_id` | `PAVS:XNNNNNNN` |
| `source_file` | Source dataset name (e.g. `ahmed-variants`) |
| `source_id` | Original row identifier from source |
| `sex_id` | NCIT URI |
| `sex_label` | `Male` / `Female` / `Unknown` |
| `age_iso` | ISO 8601 duration (e.g. `P3Y`) |
| `consanguinity_id` | `HP:0001006` or empty |
| `consanguinity_label` | `yes` / `no` / `unknown` |
| `family_history_id` | |
| `family_history_label` | |
| `test_type` | Sequencing type (WES, WGS, panel…) |
| `result_status` | `SOLVED` / `IN_PROGRESS` |
| `inheritance_id` | HPO inheritance term ID |
| `inheritance_label` | HPO inheritance label |
| `phenotypes_present_ids` | Pipe-delimited HPO IDs (present) |
| `phenotypes_present_labels` | Pipe-delimited HPO labels |
| `phenotypes_excluded_ids` | Pipe-delimited HPO IDs (negated/excluded) |
| `phenotypes_excluded_labels` | |
| `phenotypes_modifier_ids` | HPO severity/modifier IDs |
| `phenotypes_modifier_labels` | |
| `phenotypes_raw` | Original free-text phenotype string |
| `gene_symbol` | HGNC symbol |
| `gene_entrez_id` | NCBI Gene ID |
| `variant_cdna` | HGVS c. notation |
| `variant_protein` | HGVS p. notation |
| `variant_genomic_grch38` | HGVS g. notation (GRCh38) |
| `variant_hgvs_valid` | `true` / `false` (regex validation result) |
| `zygosity_id` | GENO ontology ID |
| `zygosity_label` | |
| `acmg_classification` | Normalized ACMG class string |
| `acmg_evidence` | ACMG evidence codes (PM2, PP3…) |
| `disease_mondo_id` | MONDO ID if available |
| `disease_mondo_label` | |
| `disease_omim_ids` | Pipe-delimited OMIM IDs |
| `disease_omim_labels` | |
| `disease_orphanet_ids` | Pipe-delimited Orphanet IDs |
| `pmid` | PubMed ID |
| `notes` | Free-text notes |

### Side effect: `ontology/acmg_criteria.obo`

The script also generates a minimal OBO file for the 28 standard ACMG criteria
(PVS1, PS1–PS4, PM1–PM6, PP1–PP5, BA1, BS1–BS4, BP1–BP7), with strength and
direction annotations. This is used by downstream RDF generation to provide
ACMG evidence-code URIs.

---

## Step 2 — Annotation (manual / VEP)

Between steps 1 and 3, `data/combined_normalized.tsv` is augmented with:

- VEP annotation columns (`vep_hgvsc`, `vep_hgvsp`, `vep_hgvsg`, `vep_consequence`, `vep_sift`, `vep_polyphen`)
- Population allele frequency (`saudi_af`, `gnomad_af`, `saudi_ac`, `saudi_an`)
- ClinVar columns (`clinvar_allele_id`, `clinvar_sig`)
- Gene constraint metrics (`pli`, `loeuf`)
- rsID (`vep_rsid`)

The result is `data/combined_annotated.tsv` (~77 columns), which is the authoritative
input for phenopacket generation.

---

## Step 3 — Generate Phenopackets v2 (`generate_phenopackets_v2.py`)

Converts `data/combined_annotated.tsv` to Phenopackets v2.0.2 JSON.

```bash
uv run python intake/generate_phenopackets_v2.py \
    --input data/combined_annotated.tsv \
    --output-dir phenopackets/generated_v2 \
    --combined data/PAVS_phenopackets.json \
    --zip data/PAVS_phenopackets.zip
```

### Outputs

| Path | Content |
|------|---------|
| `phenopackets/generated_v2/PAVS_MASTER_NNNN.json` | One file per case |
| `data/PAVS_phenopackets.json` | All Saudi cases as a JSON array |
| `data/PAVS_phenopackets.zip` | ZIP bundle for website download |

### Phenopacket structure

Each phenopacket includes:

```json
{
  "id": "PAVS:A0000001",
  "subject": {
    "id": "PAVS:A0000001",
    "sex": "MALE",
    "timeAtLastEncounter": { "age": { "iso8601duration": "P5Y" } }
  },
  "phenotypicFeatures": [
    { "type": { "id": "HP:0001249", "label": "Intellectual disability" }, "excluded": false },
    { "type": { "id": "HP:0001250", "label": "Seizure" }, "excluded": true },
    { "type": { "id": "HANCESTRO:0852", "label": "Middle Eastern" }, "excluded": false }
  ],
  "diseases": [
    { "term": { "id": "OMIM:272200", "label": "Multiple sulfatase deficiency" } }
  ],
  "interpretations": [{
    "id": "PAVS:A0000001",
    "progressStatus": "SOLVED",
    "diagnosis": {
      "disease": { "id": "OMIM:272200", "label": "Multiple sulfatase deficiency" },
      "genomicInterpretations": [{
        "subjectOrBiosampleId": "PAVS:A0000001",
        "interpretationStatus": "CAUSATIVE",
        "variantInterpretation": {
          "acmgPathogenicityClassification": "PATHOGENIC",
          "variationDescriptor": {
            "id": "PAVS:A0000001_SUMF1",
            "geneContext": { "symbol": "SUMF1", "valueId": "NCBIGene:285362" },
            "expressions": [
              { "syntax": "hgvs.c", "value": "NM_182760.3:c.785A>G" },
              { "syntax": "hgvs.g", "value": "NC_000003.12:g.4417183A>G" },
              { "syntax": "hgvs.p", "value": "p.Gln262Arg" }
            ],
            "vcfRecord": {
              "genomeAssembly": "hg38",
              "chrom": "chr3",
              "pos": 4417183,
              "ref": "A",
              "alt": "G"
            },
            "allelicState": { "id": "GENO:0000136", "label": "homozygous" }
          }
        }
      }]
    }
  }],
  "metaData": {
    "created": "2025-01-01",
    "createdBy": "PAVS pipeline",
    "submittedBy": "Robert Hoehndorf",
    "phenopacketSchemaVersion": "2.0",
    "resources": [ /* HP, GENO, OMIM, MONDO, HANCESTRO, GAZ */ ],
    "externalReferences": [
      { "id": "pavs:source:ahmed-variants", "description": "Source dataset: ahmed-variants" }
    ]
  }
}
```

### Key decisions

**Saudi ancestry annotation** — All Saudi sources (A, B, F, M, P, Q) receive a
`HANCESTRO:0852` ("Middle Eastern") phenotypicFeature with `excluded: false`.
DDD cases (D) do not receive this annotation.

**vcfRecord** — Extracted from the `vep_hgvsg` column by parsing
`NC_XXXXXX.YY:g.POSREF>ALT` notation. Chromosome number is mapped to `chrN` format
(NC_000023 → chrX, NC_000024 → chrY).

**ACMG normalization** — Free-text values (`"Likely pathogenic"`, `"VUS"`) are
normalized to Phenopackets v2 enum strings (`LIKELY_PATHOGENIC`,
`UNCERTAIN_SIGNIFICANCE`).

**progressStatus** — Set to `SOLVED` if `result_status == "SOLVED"`, otherwise
`IN_PROGRESS`.

**Disease preference** — MONDO ID is used if available; otherwise OMIM IDs are
used. Multiple OMIM IDs create multiple disease entries.

---

## Normalization toolkit — `phenotype_matcher_v2`

The phenotype normalization toolkit used in step 1. It maps free-text clinical
descriptions to HPO, MONDO, OMIM, and Orphanet identifiers.

Full documentation: `tools/tools/phenotype_matcher_v2/README.md`

### Install

```bash
uv pip install -e tools/phenotype_matcher_v2
```

### Quick example

```python
from phenotype_matcher_v2 import PhenotypeMatcher, MatcherConfig

cfg = MatcherConfig(hpo_path="ontology/hp.obo")
m = PhenotypeMatcher(cfg)

out = m.match("severe intellectual disability, no seizures")
print(out.get_hpo_ids())            # ['HP:0001249']
print(out.get_hpo_ids(excluded=True))   # ['HP:0001250']
print(out.phenotypes[0].severity_id)    # 'HP:0012828'  (Severe)
```

### 6-algorithm matching pipeline

`match(s)` runs the following phases in order, merging results into a single output set:

| Phase | Algorithm | Method |
|-------|-----------|--------|
| 1a | Whole-string exact | Direct lookup in `exact_map` (trusted) |
| 1b | Whole-string stemmed | Lemmatized lookup in `stemmed_map` (untrusted → LLM validation) |
| 1c | Whole-string fuzzy | Levenshtein-1 over all ontology labels (untrusted → LLM validation) |
| 2a–d | Per-token exact/stemmed/fuzzy/ANN | Same as above applied to each segment after hard segmentation |
| 2e | AC automaton span matching | Aho-Corasick over `exact_map`; finds matches at any position within token |

Untrusted matches (1b, 1c, 2b–e) are always validated by `llm_val_batch` before being accepted. Trusted matches (1a, 2a) bypass LLM validation.

### Negation and modifier detection

For each matched span, a word-window search checks for:
- **Negation** (±5 words): `no`, `not`, `without`, `absence of`, `ruled out`, `negative for`, `resolved`, `denied`, and others. Negated phenotypes are returned with `excluded=True`.
- **Severity modifiers** (±4 words): descendants of `HP:0012824` in HPO (Mild, Moderate, Severe, Profound…).

### Tokenizer (Algorithm 1)

Hard segmentation on `,  ;  .  ' but '`, followed by coordination expansion:
- `"Adj1 and Adj2 Noun"` → `["Adj1 Noun", "Adj2 Noun"]`
- `"Noun with A and B"` → `["Noun with A", "Noun with B"]`

### OntologyIndex

Loaded once per `PhenotypeMatcher` instance; cached to disk for subsequent runs.

| Structure | Contents |
|-----------|---------|
| `exact_map` | label.lower() → [term_ids] (all HPO/MONDO/OMIM/Orphanet labels + synonyms) |
| `stemmed_map` | frozenset(lemmas) → [term_ids] |
| `modifier_map` | modifier_label.lower() → modifier_term_id |
| `ac_automaton` | Aho-Corasick automaton over all keys of `exact_map` |
| `embeddings` | Per-label SapBERT embeddings (np.ndarray, cached to pkl) |
| `phenotype_ids` | Descendants of HP:0000118 (Phenotypic abnormality) |
| `modifier_ids` | Descendants of HP:0012824 (Severity) |

### Embedding models

| Preset | HuggingFace name | Speed | Accuracy |
|--------|----------------|-------|---------|
| `sapbert` (default) | `cambridgeltl/SapBERT-from-PubMedBERT-fulltext` | medium | best for clinical |
| `fast` | `all-MiniLM-L6-v2` | fast | lower |
| `balanced` | `all-mpnet-base-v2` | medium | good |
| `accurate` | `BAAI/bge-large-en-v1.5` | slow | high general |
| `medical` | `pritamdeka/S-PubMedBert-MS-MARCO` | medium | high medical |

### LLM models (via OpenRouter)

| Preset | Model | Notes |
|--------|-------|-------|
| `deepseek` (default) | `deepseek/deepseek-chat` | Free tier available |
| `gpt4oss` | `openai/gpt-4o-mini` | Low cost |
| `accurate` | `anthropic/claude-sonnet-4-5` | Highest accuracy |
| `gemini` | `google/gemini-2.0-flash-exp:free` | Free |

### CLI

```bash
# Single match, pretty output
phenotype-matcher match "severe intellectual disability, no seizures" --format pretty

# Batch process a TSV column
phenotype-matcher batch patients.tsv --column phenotypes --output results.tsv \
    --workers 4

# Ontology term lookup
phenotype-matcher info HP:0001250
phenotype-matcher info "dravet"

# Verify setup
phenotype-matcher check
```

---

## Step 4 — Generate RDF and load Virtuoso

After phenopackets are generated, the SPARQL stack pipeline continues:

```bash
# Compute HPO information content values
uv run python intake/compute_hpo_ic.py \
    --hpo ontology/hp.obo \
    --hpoa data/phenotype.hpoa \
    --output rdf_output/hpo_ic.ttl

# Generate all RDF from phenopackets + annotated TSV + literature
uv run python intake/generate_rdf.py \
    --phenopackets data/PAVS_phenopackets.json \
    --annotated data/combined_annotated.tsv \
    --literature-dir phenopackets/0.1.26 \
    --output-dir rdf_output/

# Load into Virtuoso (handled by Docker loader service)
docker compose -f docker-compose-sparql.yml up -d
```

See `README.md` for full SPARQL stack documentation.
