"""
Algorithm 1: tokenize_expand(text) → List[str]

Hard segmentation + coordination expansion for phenotype strings.
"""

import re
from typing import List


# Pattern 1: "Adj1 and Adj2 Noun" → ["Adj1 Noun", "Adj2 Noun"]
_ADJ_AND_ADJ_NOUN = re.compile(r"^(\w+)\s+and\s+(\w+)\s+(.+)$", re.IGNORECASE)

# Pattern 2: "Noun with A and B" → ["Noun with A", "Noun with B"]
_NOUN_WITH_A_AND_B = re.compile(r"^(.+?)\s+with\s+(\w+)\s+and\s+(.+)$", re.IGNORECASE)

# Pattern 3: "Multi-word phrase and another phrase" → split into two phrases.
# Guard: left side must be multi-word (has a space) so we don't split
# compound terms like "rod and cone dystrophy".
_NOUN_AND_NOUN = re.compile(r'^(.+?)\s+and\s+(.+)$', re.IGNORECASE)

# Hard-segmentation delimiter: comma, semicolon, period, or " but "
_HARD_SEG = re.compile(r"[,;.]|\s+but\s+", re.IGNORECASE)


def tokenize_expand(text: str) -> List[str]:
    """
    Algorithm 1: segment *text* and expand coordinated structures.

    Steps:
    1. Split on hard delimiters (,  ;  .  ' but ').
    2. For each segment, try coordination patterns:
       - Pattern 1: Adj1 and Adj2 Noun  → emit two tokens
       - Pattern 2: Noun with A and B   → emit two tokens
       - Pattern 3: Multi-word and B    → emit two tokens (last resort;
         left side must be multi-word to avoid splitting compounds)
       - Otherwise                      → keep as-is.
    3. Return non-empty, deduplicated list (order preserved).
    """
    segments = [s.strip() for s in _HARD_SEG.split(text)]
    seen: set = set()
    result: List[str] = []

    for seg in segments:
        if not seg:
            continue

        expansions: List[str] = []

        m1 = _ADJ_AND_ADJ_NOUN.match(seg)
        if m1:
            a1, a2, noun = m1.group(1), m1.group(2), m1.group(3).strip()
            expansions = [f"{a1} {noun}", f"{a2} {noun}"]
        else:
            m2 = _NOUN_WITH_A_AND_B.match(seg)
            if m2:
                noun_part = m2.group(1).strip()
                a = m2.group(2).strip()
                b = m2.group(3).strip()
                expansions = [f"{noun_part} with {a}", f"{noun_part} with {b}"]
            else:
                # Fix 5 — Pattern 3: split "multi-word phrase and another phrase"
                # Guard: left side must contain a space (multi-word) to avoid
                # splitting compound terms like "rod and cone dystrophy".
                m3 = _NOUN_AND_NOUN.match(seg)
                if m3 and ' ' in m3.group(1):
                    left = m3.group(1).strip()
                    right = m3.group(2).strip()
                    # Recursively expand each side (handles nested "and")
                    sub_left = tokenize_expand(left)
                    sub_right = tokenize_expand(right)
                    # Emit sub-results directly and skip normal append
                    for token in sub_left + sub_right:
                        token = token.strip()
                        if token and token not in seen:
                            seen.add(token)
                            result.append(token)
                    continue
                else:
                    expansions = [seg]

        for token in expansions:
            token = token.strip()
            if token and token not in seen:
                seen.add(token)
                result.append(token)

    return result
