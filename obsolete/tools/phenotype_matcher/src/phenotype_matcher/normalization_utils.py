"""
Normalization and cleaning utilities for phenotype matching.
"""

import re
import spacy
import pickle
import os
from functools import lru_cache
from typing import List, Set, Optional, FrozenSet, Dict, Any

# Load spacy for lemmatization
try:
    nlp = spacy.load("en_core_web_sm", disable=["parser", "ner"])
except:
    # Fallback if model is not found
    nlp = None

_persistent_cache: Dict[str, FrozenSet[str]] = {}

def load_persistent_cache(cache_path: str):
    """Load persistent cache from disk."""
    global _persistent_cache
    if os.path.exists(cache_path):
        try:
            with open(cache_path, "rb") as f:
                _persistent_cache = pickle.load(f)
        except:
            _persistent_cache = {}

def save_persistent_cache(cache_path: str):
    """Save persistent cache to disk."""
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    try:
        with open(cache_path, "wb") as f:
            pickle.dump(_persistent_cache, f)
    except:
        pass

def get_stemmed_tokens(text: Any) -> FrozenSet[str]:
    """
    Get a set of stemmed (lemmatized) tokens from text.
    Uses an in-memory persistent cache for performance.
    """
    if not text or not isinstance(text, str):
        return frozenset()
    
    text_clean = text.lower().strip()
    if text_clean in _persistent_cache:
        return _persistent_cache[text_clean]
        
    # Remove punctuation
    processed_text = re.sub(r"[^\w\s]", " ", text_clean)
    
    if nlp:
        try:
            doc = nlp(processed_text)
            result = frozenset({token.lemma_ for token in doc if not token.is_stop and not token.is_space})
        except:
            result = frozenset({t for t in processed_text.split() if t})
    else:
        # Fallback to simple split
        result = frozenset({t for t in processed_text.split() if t})
        
    _persistent_cache[text_clean] = result
    return result

def compare_stemmed(text1: str, text2: str) -> bool:
    """Compare two strings using their stemmed tokens."""
    return get_stemmed_tokens(text1) == get_stemmed_tokens(text2)
