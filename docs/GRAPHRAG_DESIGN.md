# GraphRAG Design: Key Decisions and Implementation

## Core Design Principles

### 1. Dual-Branch Retrieval Strategy

**Problem**: Severity modifiers (Mild, Moderate, Severe) are in a completely different HPO branch (HP:0012824) than phenotypes (HP:0000118).

**Solution**: Search **two separate branches** independently:

```
Clinical term: "severe intellectual disability"
                    ↓
    ┌───────────────┴───────────────┐
    ↓                               ↓
Search HP:0000118              Search HP:0012824
(Phenotypic abnormality)       (Severity)
    ↓                               ↓
Top-K phenotypes               Top-K severities
- Intellectual disability      - Severe
- Developmental delay          - Moderate  
- Cognitive impairment         - Mild
    ↓                               ↓
        LLM selects from BOTH lists
                    ↓
    Phenotype: HP:0001249 + Severity: HP:0012828
```

**Why this matters**:
- If you search ALL HPO at once, severity terms rarely appear in top-K
- Semantic similarity between "severe ID" and "Severe" is lower than between "severe ID" and other phenotypes
- Dual-branch ensures severity candidates are always available

### 2. Multi-Phenotype Output

**Problem**: A single term can describe multiple phenotypes.

**Examples**:
- "feeding difficulties and seizures" = 2 phenotypes
- "hypotonia, muscle weakness, scoliosis" = 3 phenotypes (but split by commas first)
- "developmental delay with ataxia" = 2 phenotypes

**Solution**: LLM returns an **array** of phenotypes:

```json
{
  "phenotypes": [
    {"id": "HP:0011968", "label": "Feeding difficulties", "excluded": false},
    {"id": "HP:0001250", "label": "Seizure", "excluded": false}
  ],
  "severity": null
}
```

**Processing logic**:
```python
for pheno in result.get("phenotypes", []):
    if pheno.get("excluded"):
        hpo_excluded_ids.append(pheno["id"])
    else:
        hpo_ids.append(pheno["id"])
```

### 3. Severity Label from Ontology, Never LLM

**Problem**: LLMs can hallucinate severity labels.

**Wrong approach** ❌:
```json
{"severity": {"id": "HP:0012828", "label": "Severe"}}  // LLM provides label
```

**Correct approach** ✅:
```python
# LLM provides only the ID
severity_id = result.get("severity", {}).get("id")

# Look up the label from ontology
severity_term = Ontology.get_hpo_object(severity_id)
severity_label = severity_term.name  # From ontology, not LLM
```

**Why**: The ontology is the source of truth. LLM might say "Very Severe" when HPO uses "Profound".

### 4. File-Specific Splitting Logic

**Problem**: Each source file uses different separators and formats.

**Investigation findings**:

| File | Format | Separators | Special Handling |
|------|--------|------------|------------------|
| ahmed-pmid28454995 | Free text | `,` | Standard |
| ahmed-variants | HPO IDs embedded | `,` `\|` | Strip `(HP:XXXXXXX)` before processing |
| marwa-variants | HPO IDs embedded | `\|` | Strip `(HP:XXXXXXX)` before processing |
| fawzan-variants | Free text with alternatives | `,` `/` | `/` = "or" (split both) |
| PMC6562004 | Free text, complex | `,` | May need phrase extraction |
| PMC7082194 | Free text | `,` | Standard |
| ddd-diagnoses | Already HPO IDs | N/A | Skip (already normalized) |

**Implementation**:

```python
if filename == "marwa-variants.tsv":
    # "Retinitis pigmentosa(HP:0000510)| Rod cone dystrophy(HP:0000510)"
    # → Extract text, strip IDs, split by |
    terms_with_ids = re.split(r"\|", p_val)
    raw_terms = [re.sub(r'\(HP:\d+\)', '', t).strip() for t in terms_with_ids]

elif filename == "fawzan-variants.tsv":
    # "exercise intolerance/easy fatigue, muscle weakness"
    # → Split by comma, then split by / for alternatives
    initial = re.split(r"[,;]\s*", p_val)
    expanded = []
    for t in initial:
        if "/" in t:
            expanded.extend(t.split("/"))
        else:
            expanded.append(t)
    raw_terms = expanded
```

### 5. Negation Detection: Pre-Processing + LLM Confirmation

**Two-stage approach**:

**Stage 1: Pre-detection** (fast, rule-based):
```python
excluded = any(neg in term.lower() for neg in ["no ", "not ", "without ", "absent ", "normal "])
```

**Stage 2: LLM Confirmation** (context-aware):
```
Prompt: "The term {"CONTAINS" if excluded else "DOES NOT contain"} negation words..."
LLM: Returns excluded: true/false in JSON
```

**Why both?**:
- Pre-detection catches obvious cases quickly
- LLM confirmation handles context (e.g., "not severe" vs "not present")
- Pre-detection hint helps LLM focus on the right interpretation

**Edge cases handled**:
- "no severe seizures" → excluded=true (negated), NOT severity="none"
- "not diabetic" → excluded=true
- "normal vision" → excluded=true (implies abnormality is absent)

### 6. Handling Multiple Outputs from Single Term

**Scenario**: "feeding difficulties and seizures" should map to 2 HPO terms.

**GraphRAG handles this properly**:

```python
# LLM returns multiple phenotypes
result = {
  "phenotypes": [
    {"id": "HP:0011968", "label": "Feeding difficulties", "excluded": false},
    {"id": "HP:0001250", "label": "Seizure", "excluded": false}
  ]
}

# Processing extracts ALL phenotypes
for pheno in result.get("phenotypes", []):
    hpo_ids.append(pheno["id"])
    hpo_labels.append(pheno["label"])

# Output: hpo_ids = "HP:0011968|HP:0001250"
```

**Standard pipeline comparison**:
- Standard: Splits first ("feeding difficulties", "seizures"), then normalizes each separately
- GraphRAG: Can detect multiple phenotypes in a single term PLUS handle the splits

**Both work, but GraphRAG handles complex cases better**:
- "developmental delay with ataxia" (Standard: 1 term → may miss ataxia, GraphRAG: 1 term → 2 phenotypes)

---

## Implementation Details

### Descendant Set Building

On initialization, GraphRAG builds complete descendant sets:

```python
# Get ALL terms under HP:0000118 (Phenotypic abnormality)
phenotype_root = Ontology.get_hpo_object("HP:0000118")
phenotype_descendants = {child.id for child in phenotype_root.all_children}
# Result: ~14,000 HPO IDs

# Get ALL terms under HP:0012824 (Severity)  
severity_root = Ontology.get_hpo_object("HP:0012824")
severity_descendants = {child.id for child in severity_root.all_children}
# Result: ~10 severity IDs
```

**Filtering during retrieval**:
```python
phenotype_indices = [
    i for i, t in enumerate(self.terms)
    if t["source"] == "HPO" and t["id"] in self.phenotype_descendants
]
```

This ensures:
- Only phenotypes searched when looking for phenotypes
- Only severities searched when looking for severities
- No cross-contamination

### LLM Output Schema

**GraphRAG output** (allows multiple phenotypes):
```json
{
  "phenotypes": [
    {"id": "HP:...", "label": "...", "excluded": false},
    ...
  ],
  "severity": {"id": "HP:0012828", "label": "Severe"} or null,
  "diseases": [
    {"mondo_id": "MONDO:...", "omim_id": "OMIM:...", ...}
  ]
}
```

**Standard pipeline output** (one term → one phenotype):
```json
{
  "term1": {
    "hpo": {"id": "HP:...", "label": "...", "excluded": false},
    "omim": {"id": "None", ...}
  }
}
```

**Why the difference**:
- Standard: Splits BEFORE normalization, each split gets one output
- GraphRAG: Can detect multiple phenotypes WITHIN a term

---

## Validation & Quality Control

### Severity Validation

**Rule**: Severity IDs MUST be descendants of HP:0012824.

```python
# After LLM returns severity ID
severity_id = result.get("severity", {}).get("id")
if severity_id and severity_id in self.severity_descendants:
    # Valid - use it
    # Look up label from ontology
    severity_label = Ontology.get_hpo_object(severity_id).name
else:
    # Invalid or not a severity term - ignore
    severity_id = None
```

### Negation Validation

**Two checks**:
1. Pre-detection: `any(neg in term for neg in ["no ", "not ", ...])`
2. LLM confirmation: Returns `excluded: true/false`

**Conflict resolution**: If pre-detection says yes but LLM says no → trust LLM (context-aware)

---

## Performance Optimizations

### Embedding Cache

- First run: Computes embeddings for all 35,000+ terms (~5-10 minutes)
- Subsequent runs: Loads from cache (`data/graph_rag_cache/`) (~1 second)
- Cache per model (different embeddings for different sentence transformers)

### Branch Filtering

- Pre-filtering reduces search space:
  - Phenotype search: ~14,000 terms instead of 35,000
  - Severity search: ~10 terms instead of 35,000
- Faster cosine similarity computation
- Better semantic precision (no mixing apples and oranges)

### Incremental Writing

- Each row written to disk immediately after processing
- If interrupted, resume skips already-processed rows
- No data loss on crashes

---

## Comparison Example: "severe feeding difficulties"

### Standard Pipeline
```
1. Split: ["severe feeding difficulties"]
2. Fuzzy match: Get 10 candidates for "severe feeding difficulties"
   → HP:0011968 (Feeding difficulties), HP:0001508 (Failure to thrive), ...
3. LLM: Select best match → HP:0011968
4. Validate: Confirm in source text → ✓
5. Severity: LLM extracts "severe" separately → HP:0012828
   Result: hpo_ids=HP:0011968, hpo_severity_ids=HP:0012828
```

### GraphRAG Pipeline
```
1. Pre-detect negation: NO
2. Retrieve phenotypes (HP:0000118 branch):
   → HP:0011968 (Feeding difficulties) [0.92]
   → HP:0001508 (Failure to thrive) [0.75]
   → ...
3. Retrieve severities (HP:0012824 branch):
   → HP:0012828 (Severe) [0.89]
   → HP:0012829 (Profound) [0.72]
   → ...
4. LLM selects from BOTH:
   {
     "phenotypes": [{"id": "HP:0011968", "label": "Feeding difficulties", "excluded": false}],
     "severity": {"id": "HP:0012828", "label": "Severe"}
   }
5. Look up severity label from ontology: "Severe" ✓
   Result: hpo_ids=HP:0011968, hpo_severity_ids=HP:0012828, hpo_severity_labels=Severe
```

**Both produce same output, but GraphRAG**:
- More principled (branch-aware)
- Better handles ambiguous phrasing
- Can detect multiple phenotypes in one term
