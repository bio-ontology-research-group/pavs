# PAVS Normalization Methods

This document describes the methodology used to normalize raw phenotypic and variant data from 7 source TSV files into a standardized format compatible with the PAVS database. For a step-by-step run guide, see `docs/NORMALIZATION.md`.

---

## 1. Resources & Ontologies

| Resource | Path | Used for |
|---|---|---|
| HPO (`hp.obo`) | `ontology/hp.obo` | Phenotype term mapping |
| MONDO (`mondo.obo`) | `ontology/mondo.obo` | Disease ontology |
| HPOA (`phenotype.hpoa`) | `data/phenotype.hpoa` | OMIM / Orphanet disease-phenotype associations |
| NCBI Gene info | `data/Homo_sapiens.gene_info` | Gene symbol → Entrez ID |
| SapBERT embeddings | `data/phenotype_matcher_v2_cache/` | Semantic ANN matching |
| ACMG criteria OBO | `ontology/acmg_criteria.obo` | ACMG evidence code validation |

**LLM**: OpenRouter API (DeepSeek by default; configurable). Used for candidate disambiguation and batch validation during HPO matching. Requires `OPENROUTER_API_KEY`.

---

## 2. HPO / Disease Matching — `phenotype_matcher_v2`

All phenotype-to-HPO/MONDO/OMIM/Orphanet mapping is performed by `tools/phenotype_matcher_v2/`. This is the single, current production tool (the obsolete GraphRAG v1 tool in `tools/phenotype_matcher/` is no longer used).

### Installation

```bash
uv pip install -e tools/phenotype_matcher_v2
```

### Matching pipeline overview

A free-text phenotype string is processed through six phases in sequence. Full algorithm details are in `tools/phenotype_matcher_v2/ALGORITHM.md`; the key points are:

| Phase | Method | Notes |
|---|---|---|
| 1a | Exact whole-string lookup | Trusted; no LLM |
| 1b | Stemmed whole-string (spaCy lemmas) | Untrusted; LLM validates |
| 1c | Fuzzy whole-string (edit distance 1) | Untrusted; LLM validates |
| 2a | Exact token lookup | Per token after hard segmentation |
| 2b | Stemmed token | Per token |
| 2c | Fuzzy token | Per token, skipped if 2a/2b succeeded |
| 2d | Aho-Corasick substring | Finds embedded terms (e.g. "seizure" inside longer text) |
| 2e | ANN semantic (SapBERT cosine) | Fallback for paraphrases and abbreviations |

**Tokenizer**: Hard-segments on commas, semicolons, and full stops. Also handles coordination patterns like `"short and malformed fingers"` (Pattern A), `"seizures with fever and rash"` (Pattern B), and multi-term `"and"` splits like `"severe microcephaly and epilepsy"` (Pattern C — only when left side is multi-word, to preserve compound terms like `"rod and cone dystrophy"`).

### Negation detection

Negation is detected per matched span in Phase 2d (word-window around span) and per full token in Phase 2e (whole-token search).

**Trigger vocabulary** (`V_NEG`): `no`, `not`, `never`, `none`, `nor`, `neither`, `without`, `w/o`, `denies`, `denied`, `denying`, `declines`, `absence of`, `absent`, `lack of`, `missing`, `free of`, `excluded`, `excludes`, `ruled out`, `rules out`, `negative for`, `negative`, `unremarkable for`, `resolution of`, `resolved`, `disappeared`.

**Critical rule — absence-encoding labels**: HPO terms whose canonical label starts with `"absent "`, `"absence of "`, `"missing "`, `"loss of "`, or `"lack of "` already represent the negative finding. Applying an additional negation trigger to such a term would create a double-negation. **Rule**: if the matched label starts with any of these prefixes, negation detection always returns `False` regardless of context. This applies in both Phase 2d and Phase 2e.

Examples:

| Input | Matched HPO label | Result |
|---|---|---|
| `"no speech"` | `"Absent speech"` (HP:0002371) | **present** (not excluded) |
| `"lack of eye brows"` | `"Absent eyebrows"` (HP:0000561) | **present** (not excluded) |
| `"no seizures"` | `"Seizure"` (HP:0001250) | **excluded** |

**Parenthetical qualifiers**: Context like `"microcephaly (not congenital)"` — the `(not congenital)` portion is a subtype qualifier, not a direct negation of microcephaly. Context windows strip parenthetical content before searching for triggers.

**Cross-clause scope**: Post-context is truncated at the first conjunction (`and`, `or`, `with`, `but`) to prevent negation from one clause affecting terms in a later clause. Example: `"short femur and absent radius"` — `"absent"` does not negate `"short femur"`.

### Modifier detection

Severity modifiers (HP:0012824 descendants: Mild, Moderate, Severe, Profound, etc.) are detected in the **pre-context only** (not combined pre+post) of a matched span. This prevents modifiers from an earlier term leaking onto a subsequent term in a list.

**Severity-in-label suppression**: if the matched HPO label itself already encodes a severity (e.g., `"Severe intellectual disability"` contains `"severe"`), the modifier is suppressed to avoid double-encoding.

### Parent-term subsumption

When both a parent and a descendant HPO term are matched (e.g., `"Intellectual disability"` HP:0001249 and `"Severe intellectual disability"` HP:0010864), the parent is suppressed from the output. This is computed using the `hpo_ancestors` map built at index load time.

---

## 3. Data Ingestion & Preprocessing

### 3.1 Common preprocessing (all sources)

- **ID normalization**: strip trailing `.0` from numeric IDs.
- **Sex** → NCIT terms (Male: `NCIT:C20197`, Female: `NCIT:C16576`).
- **Zygosity** → GENO terms (Homozygous: `GENO:0000136`, Heterozygous: `GENO:0000135`, Compound heterozygous: `GENO:0000402`).
- **Consanguinity** → HPO terms (Yes: `HP:0001006`).
- **Age** → ISO 8601 duration (e.g., `P5Y`) via regex parsing of free text.
- **Ancestry**: `HANCESTRO:0852` (Middle Eastern) for all Saudi cases.
- **Country**: `GAZ:00005279` (Saudi Arabia).

### 3.2 Ahmed (A/B)

- **Sources**: `ahmed*` (A/B, 525 rows total).
- **Phenotypes**: Merged from free-text sentences, `Diagnosis`, and `Additional clinical phenotype` columns; embedded `(HP:XXXXXXX)` IDs stripped before matching.
- **Variants**: HGVS-style `Gene:NM_XXXXXX:exon:c.X:p.X` or constructed from `c.DNA` and transcript columns.
- **ACMG codes**: extracted from classification strings like `"LP (PM2, PP1)"`.

### 3.3 Fawzan (`fawzan-variants.tsv`, letter F)

**Source**: 1,024 rows. Most complex source due to HGMD variant encoding.

#### Phenotype splitting

Comma- and slash-delimited. Slashes denote alternatives (`"intolerance/fatigue"`) and are split into separate tokens.

#### HGMD variant parsing

Variants may appear in a compound `Info=CLASS=DM;MUT=ALT;ID=CM981338;Ref=C;Alt=T;Start=76885871;...;dbSNP=rs111033201;Disease=Usher syndrome 1b|` string.

Parsed fields:
- `ID=` — HGMD CM accession (e.g. `CM981338`)
- `Ref=` / `Alt=` — reference and alternate alleles
- `dbSNP=` — rsID (e.g. `rs111033201`)
- `Disease=` — HGMD disease label
- `Start=` — GRCh37 start position (may be stale; GRCh38 looked up from Ensembl)

**GRCh38 coordinate lookup**: If the HGMD Info block contains a `dbSNP=` rsID, the Ensembl variation REST API is queried:
```
GET https://rest.ensembl.org/variation/human/{rs_id}?content-type=application/json
```
The GRCh38 `seq_region_name` (chromosome) and `start` position are extracted from the `mappings` array. Alleles from the HGMD `Ref`/`Alt` fields are used (not from Ensembl, since Ensembl may return multi-allelic strings).

**Duplicate variant prevention**: If a case already has a variant entry (from the `Variant(s)` column) and HGMD Info provides coordinates, the genomic field of the existing entry is updated in-place rather than appending a second variant.

#### VEP annotation for HGMD variants

When genomic coordinates `Chr{chrom}(GRCh38):g.{pos}{ref}>{alt}` are known, VEP is called via the region endpoint:
```
GET https://rest.ensembl.org/vep/human/region/{chrom}:{pos}-{pos}:1/{alt}?hgvs=1
```
This returns VEP annotations including HGVS `c.` and `p.` strings, consequence, SIFT, and PolyPhen-2 scores.

#### Disease name matching for HGMD

HGMD disease labels (e.g., `"Usher syndrome 1b"`) often use Arabic numerals where OMIM/MONDO use Roman numerals or vice versa. The standard `token_sort_ratio ≥ 80` threshold may fail for such cases. A fallback `partial_ratio ≥ 88` pass is applied when the primary match fails, allowing partial-string matching for disease family names.

### 3.5 Marwa (`marwa-variants.tsv`, letter M)

- Phenotypes: pipe-delimited; embedded HPO IDs like `HP:0001263|phenotype name` — HPO IDs extracted directly when present, free text sent to matcher otherwise.
- Variants: pipe-delimited.

### 3.6 PMC6562004 (`PMC6562004.tsv`, letter P)

- Skip first comment row (`skiprows=1`).
- Multi-gene / multi-variant columns.

### 3.7 PMC7082194 (`PMC7082194.tsv`, letter Q)

- Phenotypes in `Unnamed:13` column.
- 1-letter amino acid codes in variant protein strings converted to 3-letter HGVS `p.` notation.

### 3.8 DDD (`ddd-diagnoses.tsv`, letter D)

- HPO IDs: semicolon-delimited in `hpo_ids` column — used directly without text matching.
- Allelic mode → GENO zygosity.

---

## 4. Variant Normalization

### HGVS extraction

HGVS strings (e.g., `c.2005C>T`, `p.R669X`) are extracted via regex or HGMD Info field parsing. Format: `Transcript(Gene):c.X:p.X` when possible.

### VEP annotation

After normalization, variants with known genomic coordinates are annotated via the Ensembl VEP REST API. Two endpoint modes:

| Key type | Endpoint | Example |
|---|---|---|
| HGVS notation | `/vep/human/hgvs/{hgvs}` | `NM_000260.4:c.2005C>T` |
| Genomic region (from HGMD) | `/vep/human/region/{chrom}:{pos}-{pos}:1/{alt}` | `11:77174825-77174825:1/T` |

VEP populates: `variant_hgvs` (canonical HGVS), `variant_consequence`, `variant_sift`, `variant_polyphen`, `variant_gnomad_af`.

### ACMG classification

Classification strings normalized to:
- `"Pathogenic"`, `"Likely pathogenic"`, `"Uncertain significance"`, `"Likely benign"`, `"Benign"`

ACMG evidence codes: only valid 2015 AMP/ACMG codes accepted (PVS1, PS1–4, PM1–6, PP1–5, BA1, BS1–4, BP1–7). Non-standard codes (e.g., `"KSM"`, `"HGMD"`) are filtered out.

---

## 5. Disease Normalization

Disease strings from clinical diagnoses are mapped to OMIM / MONDO identifiers using:

1. **Exact match** against MONDO term names and synonyms.
2. **Fuzzy match** via `rapidfuzz.token_sort_ratio ≥ 80` against MONDO names.
3. **Partial-ratio fallback** via `rapidfuzz.partial_ratio ≥ 88` for cases where the disease family name partially matches (e.g., HGMD `"Usher syndrome 1b"` → MONDO `"Usher syndrome, type I"`).
4. **OMIM ID extraction**: strict regex `(?<!HP:)(?<!HP)\b\d{6}\b` to avoid false positives from HPO IDs.

---

## 6. Shared Utilities (`normalization_utils.py`)

`NormalizationUtils` is a singleton class loaded once per process:

```python
from normalization.normalization_utils import NormalizationUtils
nu = NormalizationUtils()

sex_id, sex_label = nu.normalize_sex("Male")
geno_id, geno_label = nu.normalize_genotype("homozygous")
cons_id, cons_label = nu.normalize_consanguinity("yes")
```

Loads HPO, MONDO, GENO, NCIT, OMIM, and Entrez gene data at initialization.

---

## 7. Output

**`data/combined_normalized.tsv`** — 38 columns, `PAVS:XNNNNNNN` IDs, pipe-delimited multi-value fields.

See `docs/PIPELINE.md` for the full column schema.

Key multi-value columns and pairing rules:

| Column pair | Notes |
|---|---|
| `hpo_ids` / `hpo_labels` | Present phenotype terms; parallel lists |
| `hpo_severity_ids` / `hpo_severity_labels` | Severity per phenotype; empty entry if no severity (strictly paired) |
| `hpo_excluded_ids` / `hpo_excluded_labels` | Negated/absent phenotypes; separate from present |
| `variant_hgvs` / `variant_genomic` / `variant_classification` | Pipe-delimited variant data |
| `genotype_id` / `genotype_label` | GENO zygosity terms |

---

## 8. Incremental Processing

The main script (`combine_normalize_phenotypes.py`) supports:
- Resume from interruption: skips already-processed IDs based on output file
- `--limit N`: process only the first N rows (smoke test)
- `--workers N`: parallel processing across source files
- `--overwrite`: reprocess from scratch

---

## 9. Verification

Normalization quality is checked via `analysis/verify_standardization.py`:
1. **ID validity**: all HPO and GENO IDs exist in their ontologies.
2. **HGVS syntax**: nomenclature prefix validation.
3. **Semantic count**: extracted feature count compared with original text structure.
