"""
Algorithm 2: negation and modifier detection.

det_neg(text, start, end, win)           — word-window negation (Algorithm 2a)
det_mod(text, start, end, win, mod_map)  — word-window modifier (Algorithm 2b)
det_neg_in_text(text)                    — whole-token negation  (Phase 2e)
det_mod_in_text(text, mod_map)           — whole-token modifier  (Phase 2e)
is_boundary_valid(text, start, end)      — word-boundary check
"""

import re
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Negation trigger vocabulary
# ---------------------------------------------------------------------------
V_NEG = {
    "no",
    "not",
    "never",
    "none",
    "nor",
    "neither",
    "without",
    "w/o",
    "denies",
    "denied",
    "denying",
    "declines",
    "absence of",
    "absent",
    "lack of",
    "missing",
    "free of",
    "excluded",
    "excludes",
    "ruled out",
    "rules out",
    "negative for",
    "negative",
    "unremarkable for",
    "resolution of",
    "resolved",
    "disappeared",
}

# Sort longest first so multi-word triggers are tried before their sub-strings.
_V_NEG_SORTED: List[str] = sorted(V_NEG, key=len, reverse=True)

# Clause boundary words that prevent cross-clause negation propagation.
_CLAUSE_BOUNDARIES = re.compile(r'\b(and|or|with|but)\b', re.IGNORECASE)

# HPO label prefixes that already encode absence.  When the matched term's label
# starts with one of these, the term IS the negative finding — applying an
# external negation trigger would create a double-negation, which is exceedingly
# rare in clinical text.  E.g. "no speech" → "Absent speech" should be PRESENT.
_NEG_ENCODED_LABEL_STARTS = (
    "absent ",
    "absence of ",
    "missing ",
    "loss of ",
    "lack of ",
)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _trigger_in_context(trigger: str, context: str) -> bool:
    """
    Check whether *trigger* appears in *context* with word-boundary semantics
    for single-word triggers to avoid false matches inside longer words
    (e.g. 'no' inside 'membranous').
    Multi-word triggers use plain substring matching (they are specific enough).
    """
    if ' ' in trigger:
        return trigger in context
    return bool(re.search(r'\b' + re.escape(trigger) + r'\b', context))


def _split_words(text: str) -> List[Tuple[int, int, str]]:
    """
    Return list of (word_start, word_end, word) tuples for every non-space
    run in *text*.  word_end is the index one past the last character.
    """
    return [(m.start(), m.end(), m.group()) for m in re.finditer(r"\S+", text)]


def _word_windows(
    text: str, char_start: int, char_end: int, win: int
) -> Tuple[str, str]:
    """
    Return (pre_context, post_context) each as a space-joined string of
    up to *win* words before / after the matched span [char_start, char_end).

    Fix 3b: Parenthetical content like "(not congenital)" is stripped from the
    context windows so it doesn't trigger false negation.
    """
    # Strip parenthetical qualifiers from context portions only (not the span itself)
    pre_text = re.sub(r'\([^)]*\)', ' ', text[:char_start])
    post_text = re.sub(r'\([^)]*\)', ' ', text[char_end:])
    pre_words = _split_words(pre_text)
    post_words = _split_words(post_text)
    pre_context = " ".join(w for _, _, w in pre_words[-win:])
    post_context = " ".join(w for _, _, w in post_words[:win])
    return pre_context, post_context


def _truncate_post_at_conjunction(post_context: str) -> str:
    """
    Fix 4: Truncate post-context at the first cross-clause boundary word
    (and/or/with/but) to prevent negation from leaking across clauses.
    e.g. "short femur and numerus with absent radius" — the post-context
    for "short femur" is truncated before "and", so "absent" is not seen.
    """
    m = _CLAUSE_BOUNDARIES.search(post_context)
    if m:
        return post_context[:m.start()].strip()
    return post_context


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def is_boundary_valid(text: str, start: int, end: int) -> bool:
    """
    Return True iff the characters immediately before *start* and at *end*
    are NOT alphanumeric (or are at the string boundary).

    *start* is inclusive, *end* is exclusive (i.e., text[start:end] is the
    matched span).
    """
    before_ok = (start == 0) or (not text[start - 1].isalnum())
    after_ok = (end >= len(text)) or (not text[end].isalnum())
    return before_ok and after_ok


def det_neg(
    text: str,
    start: int,
    end: int,
    win: int = 5,
    matched_label: str = "",
) -> bool:
    """
    Algorithm 2a — Negation detection via word window.

    Searches the *win*-word pre- and post-context of the span [start, end)
    in *text* for any trigger in V_NEG.

    Fix 3a: *matched_label* — if provided, triggers that appear in the
    matched term's own label are skipped (e.g. "absent" in "Absent speech"
    must not trigger negation for "no speech").

    Fix 4: Post-context is truncated at the first conjunction (and/or/with/but)
    to prevent cross-clause negation leakage.
    """
    # If the matched label already encodes absence, never mark as negated.
    # "no speech" → "Absent speech": the label IS the negative finding;
    # the external "no" is what caused the ANN match, not a separate negation.
    if matched_label:
        ll = matched_label.lower()
        if any(ll.startswith(pfx) for pfx in _NEG_ENCODED_LABEL_STARTS):
            return False

    pre, post = _word_windows(text, start, end, win)
    # Truncate post at conjunction to prevent cross-clause negation (Fix 4)
    post = _truncate_post_at_conjunction(post)
    combined_lower = (pre + " " + post).lower()
    label_lower = matched_label.lower()
    for trigger in _V_NEG_SORTED:
        if _trigger_in_context(trigger, combined_lower):
            # Fix 3a: Skip trigger if it's part of the matched term's own name
            if label_lower and _trigger_in_context(trigger, label_lower):
                continue
            return True
    return False


def det_mod(
    text: str,
    start: int,
    end: int,
    modifier_map: Dict[str, str],
    win: int = 4,
) -> Optional[str]:
    """
    Algorithm 2b — Modifier detection via word window.

    Fix 1a: Searches *win*-word PRE-context only (not pre+post) for the
    *longest* matching key in *modifier_map* (case-insensitive). Modifiers
    semantically precede the clinical term, so searching post-context caused
    modifiers from earlier terms to leak onto subsequent terms.

    Returns the corresponding modifier ID or None.
    """
    pre, _post = _word_windows(text, start, end, win)
    span = pre.lower()  # Fix 1a: pre-context only

    best_label: Optional[str] = None
    best_len = -1
    for label in modifier_map:
        if _trigger_in_context(label, span) and len(label) > best_len:
            best_label = label
            best_len = len(label)

    return modifier_map[best_label] if best_label is not None else None


def det_neg_in_text(text: str) -> bool:
    """
    Phase 2e variant: search the entire *text* for any negation trigger.
    """
    t = text.lower()
    for trigger in _V_NEG_SORTED:
        if _trigger_in_context(trigger, t):
            return True
    return False


def det_mod_in_text(
    text: str,
    modifier_map: Dict[str, str],
) -> Optional[str]:
    """
    Phase 2e variant: return the longest modifier label found anywhere in *text*.
    """
    t = text.lower()
    best_label: Optional[str] = None
    best_len = -1
    for label in modifier_map:
        if _trigger_in_context(label, t) and len(label) > best_len:
            best_label = label
            best_len = len(label)
    return modifier_map[best_label] if best_label is not None else None
