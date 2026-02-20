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


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

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
    Find which word indices overlap [char_start, char_end) in *text*.
    Return (pre_context, post_context) each as a space-joined string of
    up to *win* words before / after the matched span.
    """
    words = _split_words(text)
    # Identify span words (overlap: ws < char_end and we > char_start)
    span_indices = [
        i for i, (ws, we, _) in enumerate(words)
        if ws < char_end and we > char_start
    ]

    if not span_indices:
        return "", ""

    first_span = span_indices[0]
    last_span = span_indices[-1]

    pre_words = [w for _, _, w in words[max(0, first_span - win): first_span]]
    post_words = [w for _, _, w in words[last_span + 1: last_span + 1 + win]]

    return " ".join(pre_words), " ".join(post_words)


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


def det_neg(text: str, start: int, end: int, win: int = 5) -> bool:
    """
    Algorithm 2a — Negation detection via word window.

    Searches the *win*-word pre- and post-context of the span [start, end)
    in *text* for any trigger in V_NEG.
    """
    pre, post = _word_windows(text, start, end, win)
    combined_lower = (pre + " " + post).lower()
    for trigger in _V_NEG_SORTED:
        if trigger in combined_lower:
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

    Searches *win*-word pre+post context for the *longest* matching key in
    *modifier_map* (case-insensitive).  Returns the corresponding modifier ID
    or None.
    """
    pre, post = _word_windows(text, start, end, win)
    span = (pre + " " + post).lower()

    best_label: Optional[str] = None
    best_len = -1
    for label in modifier_map:
        if label in span and len(label) > best_len:
            best_label = label
            best_len = len(label)

    return modifier_map[best_label] if best_label is not None else None


def det_neg_in_text(text: str) -> bool:
    """
    Phase 2e variant: search the entire *text* for any negation trigger.
    """
    t = text.lower()
    for trigger in _V_NEG_SORTED:
        if trigger in t:
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
        if label in t and len(label) > best_len:
            best_label = label
            best_len = len(label)
    return modifier_map[best_label] if best_label is not None else None
