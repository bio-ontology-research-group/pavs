"""
Normalization and cleaning utilities for phenotype matching v2.
Copied from v1 and kept as-is.
"""

import re
import pickle
import os
from typing import Any, Dict, FrozenSet

try:
    import spacy
    nlp = spacy.load("en_core_web_sm", disable=["parser", "ner"])
except Exception:
    nlp = None

_persistent_cache: Dict[str, FrozenSet[str]] = {}


def load_persistent_cache(cache_path: str) -> None:
    """Load persistent lemmatization cache from disk."""
    global _persistent_cache
    if os.path.exists(cache_path):
        try:
            with open(cache_path, "rb") as f:
                _persistent_cache = pickle.load(f)
        except Exception:
            _persistent_cache = {}


def save_persistent_cache(cache_path: str) -> None:
    """Save persistent lemmatization cache to disk."""
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    try:
        with open(cache_path, "wb") as f:
            pickle.dump(_persistent_cache, f)
    except Exception:
        pass


def get_stemmed_tokens(text: Any) -> FrozenSet[str]:
    """
    Return a frozenset of lemmatized, non-stop tokens from *text*.
    Uses an in-memory cache keyed on the lowercased input.
    """
    if not text or not isinstance(text, str):
        return frozenset()

    text_clean = text.lower().strip()
    if text_clean in _persistent_cache:
        return _persistent_cache[text_clean]

    processed = re.sub(r"[^\w\s]", " ", text_clean)

    if nlp:
        try:
            doc = nlp(processed)
            result = frozenset(
                token.lemma_ for token in doc if not token.is_stop and not token.is_space
            )
        except Exception:
            result = frozenset(t for t in processed.split() if t)
    else:
        result = frozenset(t for t in processed.split() if t)

    _persistent_cache[text_clean] = result
    return result


def compare_stemmed(text1: str, text2: str) -> bool:
    """Return True iff the two strings share the same stemmed token set."""
    return get_stemmed_tokens(text1) == get_stemmed_tokens(text2)
