# Phenotype Matcher v2 — Algorithm Reference

This document describes the complete matching pipeline, every data structure, and every decision point.

---

## Overview

The matcher maps a free-text clinical string to standardised ontology identifiers across four namespaces:

- **HPO** (`HP:`) — Human Phenotype Ontology phenotype terms
- **MONDO** (`MONDO:`) — Mondo Disease Ontology disease terms
- **OMIM** (`OMIM:`) — OMIM disease phenotype entries
- **Orphanet** (`ORPHA:`) — Orphanet rare disease entries

The pipeline consists of six algorithms executed in sequence:

| # | Algorithm | Module | Function |
|---|---|---|---|
| 1 | Tokenise & Expand | `tokenizer.py` | `tokenize_expand()` |
| 2a | Negation detection (window) | `detection.py` | `det_neg()` |
| 2b | Modifier detection (window) | `detection.py` | `det_mod()` |
| 2e | Negation detection (whole text) | `detection.py` | `det_neg_in_text()` |
| 2e | Modifier detection (whole text) | `detection.py` | `det_mod_in_text()` |
| 3 | LLM disambiguation | `llm.py` | `llm_dis()` |
| 4 | LLM batched validation | `llm.py` | `llm_val_batch()` |
| 5 | Process matches | `matcher.py` | `_process_matches()` |
| 6 | Full match pipeline | `matcher.py` | `match()` |

---

## Ontology Index

All matching lookup structures are built once at startup by `OntologyIndex` (`ontology_index.py`).

### Loaded ontologies

| Source | Format | Content |
|---|---|---|
| `hp.obo` | OBO via pronto | HPO terms: names, synonyms, hierarchy |
| `mondo.obo` | OBO via pronto | MONDO terms: names, synonyms, Orphanet xrefs |
| `phenotype.hpoa` | TSV | OMIM and Orphanet disease names (columns 0–1) |

### Built structures

```
primary_label    : Dict[str, str]
    term_id → canonical name (first name registered wins)

exact_map        : Dict[str, List[str]]
    label.lower() → [term_ids]
    Contains every name and synonym for every loaded term.

stemmed_map      : Dict[FrozenSet[str], List[str]]
    frozenset(spaCy lemmas) → [term_ids]
    Keys built from every key in exact_map.

modifier_map     : Dict[str, str]
    modifier_label.lower() → modifier_term_id
    Restricted to descendants of HP:0012824 (Severity).

phenotype_ids    : Set[str]
    All descendants of HP:0000118 (Phenotypic abnormality).

modifier_ids     : Set[str]
    All descendants of HP:0012824 (Severity).

ac_automaton     : ahocorasick.Automaton
    Aho-Corasick automaton over all keys of exact_map.
    Enables O(n) multi-pattern substring search.

embeddings       : np.ndarray  (shape: N × D)
    L2-normalised SapBERT vectors for every label in exact_map.
    Cached to disk as a pickle file after first build.

embedding_labels : List[str]
    Label string for row i of embeddings.

mondo_to_orphanet: Dict[str, List[str]]
    MONDO term_id → [ORPHA: ids] from MONDO xrefs.
```

### OMIM / Orphanet loading

Lines in `phenotype.hpoa` beginning with `#` are skipped. For each remaining line, `parts[0]` is the disease ID and `parts[1]` is the disease name. Only IDs with prefix `OMIM:` or `ORPHA:` are loaded. Each disease ID is registered once (first occurrence wins).

### Orphanet via MONDO xrefs

For each MONDO term, pronto xrefs are scanned for entries beginning with `Orphanet:` or `ORPHA:`. These are normalised to the `ORPHA:` prefix and stored in `mondo_to_orphanet`. This populates the `orphanet_ids` field of `DiseaseMatch` entries.

---

## Algorithm 1 — Tokenise & Expand (`tokenize_expand`)

**Input**: raw string `text`
**Output**: `List[str]` of tokens

### Steps

1. **Hard segmentation**: split `text` on the regex `[,;.]|\s+but\s+` (comma, semicolon, full stop, or ` but `). Each resulting segment is stripped of leading/trailing whitespace.

2. **Coordination expansion**: for each non-empty segment, attempt two patterns in order:

   - **Pattern A** — `(\w+) and (\w+) (.+)`:
     "Adj1 and Adj2 Noun" → emit `"Adj1 Noun"` and `"Adj2 Noun"`
     Example: `"short and malformed fingers"` → `["short fingers", "malformed fingers"]`

   - **Pattern B** — `(.+?) with (\w+) and (.+)`:
     "Noun with A and B" → emit `"Noun with A"` and `"Noun with B"`
     Example: `"seizures with fever and rash"` → `["seizures with fever", "seizures with rash"]`

   - **No match**: keep the segment as-is.

3. **Deduplication**: tokens are added to the result only if not already seen (order preserved).

### Design note

Hard segmentation on commas is intentional: negation and modifier context must be within the same token as the phenotype term. A comma signals a new independent clause. Negation across commas (e.g. `"seizures, ruled out"`) is therefore not detected; the test cases use comma-free forms (e.g. `"seizures ruled out"`) to keep negation within a single token.

---

## Algorithm 2 — Negation and Modifier Detection

### Negation trigger vocabulary (`V_NEG`)

The following strings are negation triggers (matched as substrings after lower-casing):

```
no, not, never, none, nor, neither,
without, w/o,
denies, denied, denying, declines,
absence of, absent, lack of, missing, free of,
excluded, excludes, ruled out, rules out,
negative for, negative,
unremarkable for,
resolution of, resolved, disappeared
```

Multi-word triggers (e.g. `"ruled out"`, `"absence of"`) are tested before their sub-strings to prevent partial matches from shadowing them.

### Word-window helpers

`_split_words(text)` splits `text` on whitespace and returns `(start, end, word)` triples.

`_word_windows(text, char_start, char_end, win)`:
1. Finds which word indices overlap the character span `[char_start, char_end)`.
2. Returns `(pre, post)` where `pre` is the `win` words before the span and `post` is the `win` words after, each joined by spaces.

### Algorithm 2a — `det_neg(text, start, end, win=5)`

Returns `True` if any `V_NEG` trigger appears as a substring of the concatenated pre/post word windows (case-insensitive). The `win` parameter is words, defaulting to 5.

### Algorithm 2b — `det_mod(text, start, end, modifier_map, win=4)`

Searches the concatenated pre/post word windows (case-insensitive) for the **longest** key present in `modifier_map`. Returns the corresponding modifier ID, or `None`. Longest-match ensures `"profound"` beats `"mild"` when both appear.

### `is_boundary_valid(text, start, end)`

Returns `True` iff the character immediately before `start` (if any) and the character at `end` (if any) are both **not** alphanumeric. This prevents matching `"seizure"` inside `"myseizuredisorder"`. The `end` parameter is exclusive (Python slice convention).

### Phase 2e variants — `det_neg_in_text` / `det_mod_in_text`

Used when the token is passed directly to the ANN (Phase 2e), where there is no specific character span. These functions search the entire token string rather than a word window. `det_mod_in_text` still applies longest-match.

These variants are also used as **guards** before Phase 1b and Phase 2b (stemmed lookup): because spaCy's lemmatiser removes stop words such as `"no"`, `get_stemmed_tokens("no seizures")` returns `frozenset({'seizure'})`, which would incorrectly match HP:0001250 as present. The guard `if not det_neg_in_text(...)` skips the stemmed path when negation is present; Phase 2d (Aho-Corasick with `det_neg`) handles such tokens correctly instead.

---

## Algorithm 3 — LLM Disambiguation (`llm_dis`)

**Invoked when**: a label maps to more than one candidate ID (ambiguous label) and the match is **untrusted** (stemmed, fuzzy, or ANN path).

**Input**: `p_cand` (list of candidate IDs), `context` (token string), `primary_label` dict, optional `cue` hint
**Output**: singleton `Set[str]` with the chosen ID, or empty set on failure

**Prompt strategy**:
- System: "Select the single best concept label that matches the clinical context."
- User: clinical context + bullet list of `Label (ID)` pairs + optional cue

The LLM returns the exact label string. The function maps it back to the ID. Falls back to partial substring matching if the exact string is not found. Falls back to empty set on any failure.

---

## Algorithm 4 — Batched LLM Validation (`llm_val_batch`)

**Invoked for**: all untrusted matches (stemmed, fuzzy) — one call per token, validating all candidate labels at once.

**Input**: `p_list` (candidate IDs), `text` (token string), `primary_label` dict
**Output**: `Set[str]` of IDs whose labels are present or implied in `text`

**Prompt strategy**:
- System: "Determine which of the listed medical concepts are present or clearly implied in the given clinical text, including obvious typos, plurals, inflected forms, and abbreviations."
- User: clinical text + bulleted candidate labels
- Expected response: `{"valid": ["exact label 1", ...]}`

The "including obvious typos" clause is critical: fuzzy matches reach `llm_val_batch` with text that may contain misspellings. Without this clause the LLM rejects valid matches.

On parse failure the function falls back to substring matching of each label in the raw response string.

---

## Algorithm 5 — Process Matches (`_process_matches`)

**Input**: `labels` (list of label strings from `exact_map`), `ctx` (surrounding text), `trusted` flag
**Output**: `Set[Tuple[str, Optional[str], bool]]` — set of `(term_id, modifier_id, negated)` tuples

### Trusted path (exact matches)

All IDs for each label are added directly, without any LLM call. This preserves all namespaces: a disease name like `"phenylketonuria"` maps to `OMIM:261600`, `MONDO:0009861`, and `ORPHA:716` simultaneously — all three are included.

### Untrusted path (stemmed, fuzzy)

1. If `|p_cand| > 1`: call `llm_dis` to pick the best ID. Fall back to all candidates if `llm_dis` returns empty.
2. Call `llm_val_batch` to confirm which IDs are actually present/implied in `ctx`.
3. Add only the validated IDs.

In both paths, `modifier_id` and `negated` are set to `None`/`False` here — the caller (Phase 2d or 2e) overlays the correct values.

---

## Algorithm 6 — Full Match Pipeline (`match`)

**Input**: raw string `s`
**Output**: `PhenotypeOutput`

The accumulator `O` is a set of `(term_id, modifier_id, negated)` tuples. Deduplication is by tuple identity; the same ID with different negation flags produces separate entries.

### Phase 1 — Whole string

Applied to the full input string `s` before tokenisation.

**1a — Exact whole-string (trusted)**
If `s.lower()` is a key in `exact_map`, call `_process_matches([s_lower], s, trusted=True)` and add results to `O`.

**1b — Stemmed whole-string (untrusted)**
Guard: skip entirely if `det_neg_in_text(s_lower)` is `True`.
Compute `stems = get_stemmed_tokens(s)`. If `stems` is a key in `stemmed_map`, collect term IDs not already found by 1a. For each such ID, look up its primary label and call `_process_matches` with `trusted=False`. Only the primary label (not all synonyms) is used per ID to minimise LLM calls.

**1c — Fuzzy whole-string (untrusted)**
Skip if `s_lower ∈ W_STOP`.
Use `rapidfuzz` with `Levenshtein.distance` and `score_cutoff=1` to find labels at edit distance exactly 1. Call `_process_matches` with `trusted=False` for any hits not already covered.

### Phase 2 — Token loop

For each token `t` produced by `tokenize_expand(s)`:

**2a — Exact token (trusted)**
If `t_lower ∈ exact_map`, call `_process_matches([t_lower], t, trusted=True)`. Results go into `O_new`.

**2b — Stemmed token (untrusted)**
Guard: skip if `det_neg_in_text(t_lower)` is `True`.
Same logic as 1b but applied to token `t`. Results go into `O_new`.

After 2a/2b: merge `O_new` into `O` if non-empty. **Phase 2d always runs regardless**. This is critical: a token like `"I have a seizure disorder"` may produce a MONDO disease term via 2b while the AC automaton in 2d finds the embedded HPO phenotype term `"seizure"` (HP:0001250). If 2d were skipped after 2b, the phenotype term would be lost.

**2c — Fuzzy token (untrusted)**
Skip if `O_new` is non-empty (2a/2b already found something) or `t_lower ∈ W_STOP`.
Find labels at edit distance 1 via `_fuzzy_lookup`. Call `_process_matches`. If results found, merge into `O` and `continue` to the next token.

**2d — Aho-Corasick substring search**
Iterate all matches of the AC automaton over `t_lower`. For each match at `[start_idx, end_idx]` (inclusive):
1. Check `is_boundary_valid(t_lower, start_idx, end_idx + 1)`. Skip if invalid.
2. Compute `mod_id = det_mod(t_lower, start_idx, end_idx + 1, modifier_map, win=det_mod_win)`.
3. Compute `negated = det_neg(t_lower, start_idx, end_idx + 1, win=det_neg_win)`.
4. If `|p_cand| > 1`: call `llm_dis` to disambiguate; fall back to all candidates.
5. Add `(pid, mod_id, negated)` to `O` for each resolved ID.
6. Set `found_phrase = True`.

If `found_phrase`, `continue` to next token (skip Phase 2e).

**2e — ANN semantic search**
Encode `t` with the sentence transformer. Compute cosine similarity against all pre-built label embeddings. Retrieve top-k candidates and apply **Strategy 1** thresholding:

1. Drop all candidates with score below `θ_abs = 0.70`.
2. Let `S_max` be the highest remaining score. Keep only candidates within `[S_max − δ, S_max]` where `δ = 0.05`.
3. Cap at `k = 5` candidates.

For each surviving label, apply `det_neg_in_text(t_lower)` and `det_mod_in_text(t_lower, modifier_map)` to the whole token. Add `(pid, mod_id, negated)` to `O` for each ID in `exact_map[label]`.

### Output construction (`_build_output`)

Route each `(term_id, mod_id, negated)` tuple to the appropriate output structure:

| Condition | Output |
|---|---|
| `term_id.startswith("HP:")` and `term_id ∈ phenotype_ids` | `PhenotypeMatch(hpo_id, label, excluded=negated, severity_id=mod_id, ...)` |
| `term_id.startswith("MONDO:")` | `DiseaseMatch(mondo_id=term_id, orphanet_ids=mondo_to_orphanet.get(term_id, []))` |
| `term_id.startswith("OMIM:")` | `DiseaseMatch(omim_phenotype_ids=[term_id], omim_labels=[label])` |
| `term_id.startswith("ORPHA:")` | `DiseaseMatch(orphanet_ids=[term_id])` |

Deduplication is by `term_id` within each output list.

---

## Stop-word Set (`W_STOP`)

Tokens in `W_STOP` skip Phase 1c (fuzzy whole-string), Phase 2c (fuzzy token), and Phase 2e (ANN) to prevent meaningless semantic matches on generic words. Members include common English function words and generic clinical vocabulary:

```
and, or,
age, sex, pain, type, the, a, an, is, are, was, were, be, been,
have, has, had, do, does, did,
to, of, in, on, at, by, for, with, about, as, into, from, out,
case, report, patient, diagnosis, clinical, finding, feature,
sign, symptom, present, status
```

---

## Strategy 1 — Hybrid Dynamic Thresholding

Applied to ANN scores in Phase 2e.

```
floor = max(θ_abs, S_max − δ)
keep  = {label : score[label] ≥ floor}[:cap]
```

With defaults `θ_abs = 0.70`, `δ = 0.05`, `cap = 5`:
- Absolute floor prevents garbage matches when the best score is low.
- Relative margin keeps near-synonyms (e.g. `"kidney failure"` and `"renal failure"`) while discarding distant candidates.
- Cap bounds the number of results per token.

---

## LLM Integration

All LLM calls go to [OpenRouter](https://openrouter.ai/) via a single `query_llm(system, user, api_key, model)` function. The API key is read from the `OPENROUTER_API_KEY` environment variable or passed via `MatcherConfig.api_key`.

Available model presets:

| Preset | Model | Notes |
|---|---|---|
| `deepseek` | `deepseek/deepseek-chat` | Default; free tier via OpenRouter |
| `gpt4oss` | `openai/gpt-4o-mini` | Cheap, good quality |
| `accurate` | `anthropic/claude-sonnet-4-5` | Highest accuracy |
| `gemini` | `google/gemini-2.0-flash-exp:free` | Free alternative |

Responses are parsed as JSON. Markdown code fences are stripped. On any network or parse failure the LLM functions return empty results and the pipeline continues without LLM input.

---

## Embedding Models

SapBERT (`cambridgeltl/SapBERT-from-PubMedBERT-fulltext`) is the default. Embeddings are built once and cached to disk as a pickle file in `cache_dir`. The cache file name encodes the model name.

Available presets:

| Preset | Model | Use case |
|---|---|---|
| `sapbert` | `cambridgeltl/SapBERT-from-PubMedBERT-fulltext` | Default; trained on biomedical synonyms |
| `fast` | `all-MiniLM-L6-v2` | Fast, lower accuracy |
| `balanced` | `all-mpnet-base-v2` | Good general purpose |
| `accurate` | `BAAI/bge-large-en-v1.5` | High accuracy, slower |
| `medical` | `pritamdeka/S-PubMedBert-MS-MARCO` | Medical domain |
