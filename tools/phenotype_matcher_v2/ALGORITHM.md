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

hpo_ancestors    : Dict[str, Set[str]]
    term_id → set of all ancestor HP IDs (excludes term itself).
    Used for parent-term subsumption in _build_output().

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

2. **Coordination expansion**: for each non-empty segment, attempt three patterns in order:

   - **Pattern A** — `(\w+) and (\w+) (.+)`:
     "Adj1 and Adj2 Noun" → emit `"Adj1 Noun"` and `"Adj2 Noun"`
     Example: `"short and malformed fingers"` → `["short fingers", "malformed fingers"]`

   - **Pattern B** — `(.+?) with (\w+) and (.+)`:
     "Noun with A and B" → emit `"Noun with A"` and `"Noun with B"`
     Example: `"seizures with fever and rash"` → `["seizures with fever", "seizures with rash"]`

   - **Pattern C** — `(.+?) and (.+)` where the left side is multi-word:
     Splits "PhenotypeA and PhenotypeB" when both sides are independent terms.
     Guard: left side must contain a space (prevents splitting compound terms like `"rod and cone dystrophy"`).
     Example: `"severe microcephaly and epilepsy"` → `["severe microcephaly", "epilepsy"]`

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

### Absence-encoding label prefixes (`_NEG_ENCODED_LABEL_STARTS`)

HPO labels that already encode a negative/absent finding:

```
"absent ", "absence of ", "missing ", "loss of ", "lack of "
```

When the matched term's canonical label starts with one of these prefixes, the term **is** the negative finding — applying an external negation trigger would create a double-negation, which is exceedingly rare in clinical text. Example: `"no speech"` semantically matches `"Absent speech"` (HP:0002371). The external `"no"` caused the ANN match; it must not additionally mark the term as excluded.

**Rule**: if `matched_label.lower().startswith(any prefix in _NEG_ENCODED_LABEL_STARTS)`, `det_neg()` returns `False` immediately without examining context windows.

### Word-window helpers

`_split_words(text)` splits `text` on whitespace and returns `(start, end, word)` triples.

`_word_windows(text, char_start, char_end, win)`:
1. Strips parenthetical content from the pre- and post-context portions only (not the matched span itself), using `re.sub(r'\([^)]*\)', ' ', ...)`. This prevents qualifiers like `"(not congenital)"` from triggering negation for the surrounding term.
2. Returns `(pre, post)` where `pre` is the `win` words before the span and `post` is the `win` words after, each joined by spaces.

### Algorithm 2a — `det_neg(text, start, end, win=5, matched_label="")`

Checks whether the matched span at `[start, end)` is negated:

1. **Absence-encoding early exit**: if `matched_label` starts with any prefix in `_NEG_ENCODED_LABEL_STARTS`, return `False` immediately (prevents double-negation of terms like "Absent speech").

2. Compute `(pre, post)` via `_word_windows`.

3. **Cross-clause truncation**: apply `_truncate_post_at_conjunction(post)` — truncate the post-context at the first occurrence of `and`, `or`, `with`, or `but`. This prevents negation from leaking across independent clauses (e.g., `"short femur and absent radius"` — for `"short femur"`, the post-context is truncated before `"and"`, so `"absent"` in the next clause is not seen).

4. Form `combined_lower = (pre + " " + post).lower()`.

5. For each trigger in `V_NEG` (longest first):
   - If trigger is found in `combined_lower`:
     - **Label-bypass**: if `matched_label` is provided and the same trigger appears in `matched_label.lower()`, skip (the trigger is part of the term's own name, not an external negation). Example: `"absent"` in `"Absent speech"` must not trigger negation for `"no speech"`.
     - Otherwise return `True`.

6. Return `False`.

### Algorithm 2b — `det_mod(text, start, end, modifier_map, win=4)`

Searches the **pre-context only** (not combined pre+post) of the `win`-word window for the **longest** key present in `modifier_map`. Returns the corresponding modifier ID, or `None`.

**Pre-context only**: modifiers semantically precede the clinical term (e.g., `"severe microcephaly"`). Searching post-context caused modifiers from an earlier term to leak onto subsequent terms in a list (e.g., `"moderate developmental delay, autistic behavior"` — `"moderate"` is in the post-context window for `"autistic behavior"` when using combined pre+post).

### `is_boundary_valid(text, start, end)`

Returns `True` iff the character immediately before `start` (if any) and the character at `end` (if any) are both **not** alphanumeric. This prevents matching `"seizure"` inside `"myseizuredisorder"`. The `end` parameter is exclusive (Python slice convention).

### Phase 2e variants — `det_neg_in_text` / `det_mod_in_text`

Used when the token is passed directly to the ANN (Phase 2e), where there is no specific character span. These functions search the entire token string rather than a word window. `det_mod_in_text` still applies longest-match.

These variants are also used as **guards** before Phase 1b and Phase 2b (stemmed lookup): because spaCy's lemmatiser removes stop words such as `"no"`, `get_stemmed_tokens("no seizures")` returns `frozenset({'seizure'})`, which would incorrectly match HP:0001250 as present. The guard `if not det_neg_in_text(...)` skips the stemmed path when negation is present; Phase 2d (Aho-Corasick with `det_neg`) handles such tokens correctly instead.

**Note**: In Phase 2e, `det_neg_in_text` returns a single `negated` flag for the whole token. However, each surviving ANN label is then individually checked against `_NEG_ENCODED_LABEL_STARTS`: if the label already encodes absence, `lbl_negated` is forced to `False` for that specific match, overriding the token-level flag. This is essential for tokens like `"no speech"` which semantically match `"Absent speech"`.

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
2. Retrieve the primary label `label` for the matched term.
3. Compute `mod_id = det_mod(t_lower, start_idx, end_idx + 1, modifier_map, win=det_mod_win)`.
4. Compute `negated = det_neg(t_lower, start_idx, end_idx + 1, win=det_neg_win, matched_label=label)`. The `matched_label` parameter enables the absence-encoding bypass and label-bypass rules (see Algorithm 2a).
5. If `|p_cand| > 1`: call `llm_dis` to disambiguate; fall back to all candidates.
6. Add `(pid, mod_id, negated)` to `O` for each resolved ID.
7. Set `found_phrase = True`.

If `found_phrase`, `continue` to next token (skip Phase 2e).

**2e — ANN semantic search**
Encode `t` with the sentence transformer. Compute cosine similarity against all pre-built label embeddings. Retrieve top-k candidates and apply **Strategy 1** thresholding:

1. Drop all candidates with score below `θ_abs = 0.70`.
2. Let `S_max` be the highest remaining score. Keep only candidates within `[S_max − δ, S_max]` where `δ = 0.05`.
3. Cap at `k = 5` candidates.

Compute a single token-level `negated = det_neg_in_text(t_lower)` and `mod_id = det_mod_in_text(t_lower, modifier_map)`.

For each surviving label:
- Compute `lbl_negated = negated`.
- **Absence-encoding per-label override**: if `negated` is `True` and `label.lower()` starts with any prefix in `_NEG_ENCODED_LABEL_STARTS`, set `lbl_negated = False`. This ensures that tokens like `"no speech"` which semantically match `"Absent speech"` are not doubly-negated.
- Add `(pid, mod_id, lbl_negated)` to `O` for each ID in `exact_map[label]`.

### Output construction (`_build_output`)

Converts the `(term_id, mod_id, negated)` accumulator `O` into a `PhenotypeOutput`.

**Step 1 — Severity-in-label suppression**

Before routing, check whether a matched HPO phenotype's primary label already contains a severity word (e.g., `"Severe intellectual disability"` contains `"severe"`). If so, clear `mod_id` to avoid double-encoding severity. Recognised severity words: `severe`, `mild`, `moderate`, `profound`, `borderline`, `marked`.

**Step 2 — Parent-term subsumption**

After collecting all present (non-negated) HPO term IDs, suppress any term that is an ancestor of another matched term. For example, if both `HP:0001249` (Intellectual disability) and `HP:0010864` (Severe intellectual disability) are matched:
- `HP:0010864` is a descendant of `HP:0001249`.
- `HP:0001249` is in `hpo_ancestors[HP:0010864]`.
- Therefore `HP:0001249` is added to the suppressed set and omitted from the output.

This is computed using `idx.hpo_ancestors`: for each present term `tid`, if any other present term `other` has `tid` in its ancestor set, `tid` is suppressed.

**Step 3 — Routing**

Route each surviving tuple to the appropriate output structure:

| Condition | Output |
|---|---|
| `term_id.startswith("HP:")` and `term_id ∈ phenotype_ids` (not suppressed) | `PhenotypeMatch(hpo_id, label, excluded=negated, severity_id=mod_id, ...)` |
| `term_id.startswith("MONDO:")` | `DiseaseMatch(mondo_id=term_id, orphanet_ids=mondo_to_orphanet.get(term_id, []))` |
| `term_id.startswith("OMIM:")` | `DiseaseMatch(omim_phenotype_ids=[term_id], omim_labels=[label])` |
| `term_id.startswith("ORPHA:")` | `DiseaseMatch(orphanet_ids=[term_id])` |

Deduplication is by `term_id` within each output list.

### LLM post-hoc validation (`_llm_validate`)

If `cfg.api_key` is set, an optional post-hoc pass is applied after `_build_output`:

- **Negation validation** (`cfg.llm_validate_negation=True`): for each excluded phenotype, the LLM confirms whether the term is truly absent in context. If the LLM disagrees, the phenotype is flipped to present.
- **Modifier validation** (`cfg.llm_validate_modifiers=True`): for each phenotype with a severity modifier, the LLM confirms whether the modifier correctly describes that phenotype's severity in context. If not, the modifier is cleared.

Disabled (no-op) when `api_key` is absent so the pipeline stays fully offline by default. Results are cached per `(term_label, context)` pair to minimise API calls.

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

---

## Known Edge Cases and Rules

### Double-negation prevention

Clinical terms whose HPO label starts with an absence-encoding prefix (`"absent "`, `"absence of "`, `"missing "`, `"loss of "`, `"lack of "`) represent the negative finding as the canonical phenotype. Examples:

| Input text | ANN match | Correct behaviour |
|---|---|---|
| `"no speech"` | `"Absent speech"` (HP:0002371) | present (not excluded) |
| `"lack of eye brows"` | `"Absent eyebrows"` (HP:0000561) | present (not excluded) |
| `"loss of vision"` | `"Loss of vision"` (HP:0000572) | present (not excluded) |

The `_NEG_ENCODED_LABEL_STARTS` guard in `det_neg()` and the per-label override in Phase 2e both implement this rule.

### Parenthetical qualifiers

Expressions like `"microcephaly (not congenital)"` contain a subtype qualifier in parentheses. This qualifier must not trigger negation for the surrounding term. `_word_windows()` strips parenthetical content before computing context windows.

### Cross-clause negation

`"short femur and numerus with absent radius and tibia"` — for the span `"short femur"`, the post-context is truncated before `"and"`, preventing `"absent"` from the next clause from marking the femur term as negated.

### Modifier scope

`"moderate developmental delay, autistic behavior"` — `"moderate"` is in the pre-context of `"autistic behavior"` when the window is 4 words. `det_mod()` uses **pre-context only**, so it does not see `"moderate"` in post-context of earlier terms. Comma segmentation also limits cross-term leakage.

### Severity double-encoding

`"Severe intellectual disability"` matches HP:0010864 directly (the HPO term for severe intellectual disability). Applying `"severe"` as a modifier on top would double-encode severity. `_build_output()` suppresses `mod_id` when the matched label already contains a severity word.

### Parent-term redundancy

When both `"Intellectual disability"` (HP:0001249) and `"Severe intellectual disability"` (HP:0010864) are matched, the parent HP:0001249 is suppressed because HP:0010864 is in its descendant set.
