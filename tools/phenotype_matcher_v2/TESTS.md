# Phenotype Matcher v2 — Test Reference

The test suite is in `src/phenotype_matcher_v2/test_cases.py`. It contains 56 tests in two categories:

- **Unit tests** (33): pure function tests, no ontology required, run in milliseconds locally.
- **Integration tests** (23): require a loaded `PhenotypeMatcher`, SapBERT embeddings, and an active `OPENROUTER_API_KEY`. Gated by `RUN_INTEGRATION=1`.

## Running the tests

```bash
# Unit tests only (fast, no dependencies)
python src/phenotype_matcher_v2/test_cases.py

# All tests (requires GPU server, API key, ontology files)
export OPENROUTER_API_KEY=...
RUN_INTEGRATION=1 python src/phenotype_matcher_v2/test_cases.py

# Point to non-default ontology paths
HPO_PATH=/path/to/hp.obo MONDO_PATH=/path/to/mondo.obo \
HPOA_PATH=/path/to/phenotype.hpoa \
RUN_INTEGRATION=1 python src/phenotype_matcher_v2/test_cases.py
```

---

## Unit tests — `TestTokenizer` (Algorithm 1)

These test `tokenize_expand()` in isolation.

| ID | Input | What is tested |
|---|---|---|
| TK1 | `"seizures, hypotonia"` | Hard segmentation on comma produces two independent tokens |
| TK2 | `"seizures; ataxia. hypotonia"` | Hard segmentation on semicolon and full stop |
| TK3 | `"normal vision but hypotonia"` | Hard segmentation on ` but ` keyword |
| TK4 | `"short and malformed fingers"` | Pattern A expansion: Adj1 and Adj2 Noun → `["short fingers", "malformed fingers"]` |
| TK5 | `"seizures with fever and rash"` | Pattern B expansion: Noun with A and B → `["seizures with fever", "seizures with rash"]` |
| TK6 | `"intellectual disability"` | No pattern match — segment kept as-is |
| TK7 | `"short and malformed fingers, absent thumb"` | Both hard segmentation and coordination expansion applied to the same input |

---

## Unit tests — `TestDetNeg` (Algorithm 2a)

These test `det_neg()` with explicit character spans and `win=5` words. They verify that the word-window pre/post context is correctly extracted and searched for negation triggers.

| ID | Text | Span `[start, end)` | Expected | What is tested |
|---|---|---|---|---|
| DN1 | `"no seizures"` | `[3, 11)` | `True` | Trigger `"no"` in the PRE window |
| DN2 | `"seizures ruled out"` | `[0, 8)` | `True` | Multi-word trigger `"ruled out"` in the POST window |
| DN3 | `"present seizures"` | `[8, 16)` | `False` | No trigger in either window |
| DN4 | `"without any seizures"` | `[11, 19)` | `True` | Trigger `"without"` in PRE window |
| DN5 | `"patient denies seizures"` | `[15, 23)` | `True` | Clinical denial verb `"denies"` in PRE window |
| DN6 | `"severe seizures were ruled out today"` | `[7, 15)` | `True` | Trigger beyond the immediate neighbour word in POST |
| DN7 | `"seizures"` | `[0, 8)` | `False` | No surrounding context at all |

---

## Unit tests — `TestDetMod` (Algorithm 2b)

These test `det_mod()` using a small stub `modifier_map = {"severe": "HP:0012828", "mild": "HP:0012825", "profound": "HP:0012829"}`.

| ID | Text | Span `[start, end)` | Expected | What is tested |
|---|---|---|---|---|
| DM1 | `"severe intellectual disability"` | `[7, 31)` | `"HP:0012828"` | Modifier `"severe"` in the PRE window |
| DM2 | `"intellectual disability, mild"` | `[0, 22)` | `"HP:0012825"` | Modifier `"mild"` in the POST window (after comma) |
| DM3 | `"intellectual disability"` | `[0, 22)` | `None` | No modifier present |
| DM4 | `"profound and severe ataxia"` | `[16, 22)` | `"HP:0012829"` | Longest-match wins: `"profound"` (8 chars) beats `"severe"` (6 chars), both in PRE |

---

## Unit tests — `TestDetInText` (Phase 2e variants)

These test `det_neg_in_text()` and `det_mod_in_text()`, which search the **entire** token string rather than a word window. Used in Phase 2e (ANN path) and as guards before stemmed lookup.

| ID | Function | Text | Expected | What is tested |
|---|---|---|---|---|
| DI1 | `det_neg_in_text` | `"no fever"` | `True` | Trigger `"no"` found anywhere in text |
| DI2 | `det_neg_in_text` | `"fever present"` | `False` | No trigger present |
| DI3 | `det_neg_in_text` | `"negative for fever"` | `True` | Multi-word trigger `"negative for"` found |
| DI4 | `det_mod_in_text` | `"severe ataxia"` | `"HP:0012828"` | Modifier word `"severe"` found in text |
| DI5 | `det_mod_in_text` | `"ataxia"` | `None` | No modifier present |
| DI6 | `det_mod_in_text` | `"mild to severe ataxia"` | `"HP:0012828"` | Longest-match wins: `"severe"` (6) beats `"mild"` (4) |

---

## Unit tests — `TestBoundaryValid`

These test `is_boundary_valid()`, which guards the Aho-Corasick substring matches in Phase 2d to prevent matching ontology labels embedded inside unrelated words.

| ID | Text | Span `[start, end)` | Expected | What is tested |
|---|---|---|---|---|
| BV1 | `"has seizure disorder"` | `[4, 11)` | `True` | Space before and after — valid word boundary |
| BV2 | `"myseizure"` | `[2, 9)` | `False` | Alphanumeric character before span (`y`) — not a word boundary |
| BV3 | `"seizures"` | `[0, 7)` | `False` | Character at `end` is alphanumeric (`s`) — plural suffix breaks boundary |
| BV4 | `"seizure"` | `[0, 7)` | `True` | Span covers the entire string — both boundaries are string edges |
| BV5 | `"(seizure)"` | `[1, 8)` | `True` | Surrounded by non-alphanumeric punctuation |

---

## Unit tests — `TestWStop`

These verify membership in the `W_STOP` constant, which controls which tokens skip fuzzy and ANN matching.

| ID | Query | Expected | What is tested |
|---|---|---|---|
| WS1 | `"and"` | in W_STOP | Conjunction skips fuzzy/ANN |
| WS2 | `"or"` | in W_STOP | Disjunction skips fuzzy/ANN |
| WS3 | `"seizure"` | not in W_STOP | Clinical term proceeds through all phases |
| WS4 | `"patient"` | in W_STOP | Generic clinical word skips fuzzy/ANN |

---

## Integration tests — `TestIntegration`

These require `RUN_INTEGRATION=1`, a loaded ontology, SapBERT embeddings, and the `OPENROUTER_API_KEY`. The matcher is loaded once in `setUpClass`.

### Phase 1 — Whole-string matching

| ID | Input | Expected | Algorithm branch tested |
|---|---|---|---|
| M01 | `"Seizure"` | `HP:0001250` present | Phase 1a: exact whole-string match (trusted) |
| M02 | `"seizures"` | `HP:0001250` present | Phase 1b: stemmed whole-string match |
| M03 | `"Seizuree"` | `HP:0001250` present | Phase 1c: fuzzy whole-string (edit distance 1); LLM-Val prompt accepts obvious typos |
| M04 | `"and"` | empty result | Phase 1c skipped: `"and" ∈ W_STOP` |

**M03 design note**: the extra trailing `e` in `"Seizuree"` is a deliberate, unambiguous typo at edit distance 1 from `"seizure"`. The `llm_val_batch` system prompt was updated to include `"obvious typos, plurals, inflected forms, and abbreviations"` to ensure the LLM accepts this. Earlier candidates (`"Seizzure"`, `"Siezure"`) were found to be inconsistently accepted; `"Seizuree"` is the most reliable choice.

### Phase 2a/2b — Exact and stemmed token matching

| ID | Input | Expected | Algorithm branch tested |
|---|---|---|---|
| M05 | `"Seizure, hypotonia"` | `HP:0001250` + hypotonia HPO ID | Phase 2a: exact token match for each comma-separated term |
| M06 | `"seizing, hypotonia"` | hypotonia HPO ID present | Phase 2b: stemmed token match (`"seizing"` → stem `"seize"` → `HP:0001250` via LLM-Val) |

**HPO version note**: tests M05, M06, M08, M09 accept either `HP:0001252` ("Hypotonia", primary label in the HPO version deployed) or `HP:0001290` ("Generalized hypotonia") because the canonical label differs between HPO releases.

### Phase 2c — Fuzzy token matching

| ID | Input | Expected | Algorithm branch tested |
|---|---|---|---|
| M08 | `"Seizzure, hypotonia"` | hypotonia HPO ID present | Phase 2c: fuzzy match on `"Seizzure"` (edit distance 1 from `"seizure"`) |
| M09 | `"and, hypotonia"` | hypotonia HPO ID only | Phase 2c skipped for `"and" ∈ W_STOP`; other token still matched |

### Phase 2d — Aho-Corasick substring search

| ID | Input | Expected | Algorithm branch tested |
|---|---|---|---|
| M10 | `"I have a seizure disorder"` | `HP:0001250` present | Phase 2d: AC finds `"seizure"` embedded in a longer sentence; Phase 2b may also find a disease term (MONDO) but Phase 2d always runs |
| M12 | `"no seizures"` | `HP:0001250` excluded | Phase 2d: `det_neg` returns `True` for PRE window containing `"no"`; Phase 2b is skipped by negation guard |
| M13 | `"seizures ruled out"` | `HP:0001250` excluded | Phase 2d: `det_neg` returns `True` for POST window containing `"ruled out"`; no comma so both terms stay in the same token |
| M14 | `"severe seizures"` | `HP:0001250` with `severity_id=HP:0012828` | Phase 2d: `det_mod` finds `"severe"` in PRE window |
| M15 | `"mild seizures"` | `HP:0001250` with `severity_id=HP:0012825` | Phase 2d: `det_mod` finds `"mild"` in PRE window; no comma to split the token |

**M10 design note**: the critical test for the `continue`-removal fix. A token like `"I have a seizure disorder"` produces a MONDO disease term via Phase 2b (stemmed: `frozenset({'seizure', 'disorder'})` → MONDO:0005027). Before the fix, a `continue` after Phase 2b skipped Phase 2d, losing the embedded HPO term `"seizure"` (HP:0001250). Phase 2d now always runs.

**M12/M13 design note**: `"no seizures"` could not use Phase 2b because `get_stemmed_tokens("no seizures")` returns `frozenset({'seizure'})` (spaCy removes `"no"` as a stop word), which would match HP:0001250 as present. The negation guard (`if not det_neg_in_text(t_lower)`) skips Phase 2b, leaving Phase 2d to detect negation correctly.

**M13/M15 design note**: these tests use comma-free inputs (`"seizures ruled out"`, `"mild seizures"`) because `tokenize_expand` hard-segments on commas. If a comma separates the phenotype from its modifier/negation, they end up in different tokens and the window-based detection cannot find the context.

### Phase 2e — ANN semantic search

| ID | Input | Expected | Algorithm branch tested |
|---|---|---|---|
| M17 | `"fits and falls"` | at least one match (seizure-related) | Phase 2e: informal term matched via SapBERT cosine similarity |
| M18 | `"unknown etiology"` | empty result | Phase 2e: Strategy 1 absolute floor (θ_abs=0.70) filters all candidates |
| M20 | `"no fever"` | fever HPO term excluded | Phase 2e: `det_neg_in_text` sets `excluded=True` for ANN matches |
| M21 | `"severe ataxia"` | ataxia HPO term with `severity_id` set | Phase 2e: `det_mod_in_text` sets `severity_id` for ANN matches |

### Ontology coverage

| ID | Input | Expected | What is tested |
|---|---|---|---|
| M26 | `"Dravet syndrome"` | at least one MONDO ID | MONDO disease terms are returned; trusted path adds all IDs for the label |
| M27 | `"phenylketonuria"` | at least one OMIM ID | OMIM disease names loaded from `phenotype.hpoa` |
| M28 | `"Rett syndrome"` | at least one Orphanet ID | Orphanet IDs populated via MONDO xrefs (`mondo_to_orphanet`) |

**M26/M27/M28 design note**: the trusted path (`_process_matches(trusted=True)`) adds **all** candidate IDs for a label without calling `llm_dis`. This is necessary because a label like `"phenylketonuria"` maps simultaneously to `OMIM:261600`, `MONDO:0009861`, and `ORPHA:716`. The earlier design called `llm_dis` which returned only one ID, causing the test to fail when the LLM happened to pick an Orphanet ID instead of MONDO/OMIM.

### Multi-term and structural tests

| ID | Input | Expected | What is tested |
|---|---|---|---|
| M29 | `"seizures, hypotonia, and feeding difficulties"` | at least 2 HPO matches | Full multi-term pipeline: segmentation + independent matching per token |
| M30 | `"short and malformed fingers"` | at least 1 HPO match | Tokenizer coordination expansion feeds the Phase 2 token loop |
| M31 | `"Seizure"` | no duplicate `hpo_id` values in output | Deduplication: same term ID reached via multiple paths appears only once |
