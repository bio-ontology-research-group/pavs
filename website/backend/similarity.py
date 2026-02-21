"""
similarity.py — HPO phenotype similarity (Lin + BMA, no pyhpo dependency).

IC values and ancestor sets are pre-loaded from Virtuoso at startup.
"""

from __future__ import annotations
import math
from typing import Dict, List, Set


def _mica_ic(a: str, b: str, ic: Dict[str, float], ancestors: Dict[str, Set[str]]) -> float:
    """Information content of the Most Informative Common Ancestor."""
    anc_a = ancestors.get(a, {a})
    anc_b = ancestors.get(b, {b})
    common = anc_a & anc_b
    if not common:
        return 0.0
    return max(ic.get(t, 0.0) for t in common)


def resnik_sim(a: str, b: str, ic: Dict[str, float], ancestors: Dict[str, Set[str]]) -> float:
    return _mica_ic(a, b, ic, ancestors)


def lin_sim(a: str, b: str, ic: Dict[str, float], ancestors: Dict[str, Set[str]]) -> float:
    mica = _mica_ic(a, b, ic, ancestors)
    denom = ic.get(a, 0.0) + ic.get(b, 0.0)
    if denom == 0.0:
        return 0.0
    return (2.0 * mica) / denom


def _bma(
    q_terms: List[str],
    t_terms: List[str],
    ic: Dict[str, float],
    ancestors: Dict[str, Set[str]],
    method: str,
) -> float:
    """Best-Match Average in one direction (q → t)."""
    sim_fn = lin_sim if method == "lin" else resnik_sim
    scores = []
    for q in q_terms:
        best = max((sim_fn(q, t, ic, ancestors) for t in t_terms), default=0.0)
        scores.append(best)
    return sum(scores) / len(scores) if scores else 0.0


def bma_similarity(
    query_hpos: List[str],
    target_hpos: List[str],
    ic: Dict[str, float],
    ancestors: Dict[str, Set[str]],
    method: str = "lin",
) -> float:
    """
    Symmetric Best-Match Average similarity.
    method: "lin" (default) or "resnik"
    """
    if not query_hpos or not target_hpos:
        return 0.0

    # Expand with ancestors for richer matching
    def expand(terms):
        expanded = set()
        for t in terms:
            expanded.add(t)
            expanded |= ancestors.get(t, set())
        return list(expanded)

    q_exp = expand(query_hpos)
    t_exp = expand(target_hpos)

    q_to_t = _bma(q_exp, t_exp, ic, ancestors, method)
    t_to_q = _bma(t_exp, q_exp, ic, ancestors, method)
    return (q_to_t + t_to_q) / 2.0
