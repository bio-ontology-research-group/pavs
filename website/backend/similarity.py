"""
similarity.py — HPO phenotype similarity (Lin + BMA, no pyhpo dependency).

IC values and ancestor sets are pre-loaded from Virtuoso at startup.
"""

from __future__ import annotations
from typing import Dict, List, Set, Tuple

# Module-level MICA cache: canonical (a, b) tuple → max IC of common ancestors.
# Populated lazily; shared across all calls in the process.
_mica_cache: Dict[Tuple[str, str], float] = {}


def _mica_ic_cached(
    a: str, b: str,
    ic: Dict[str, float],
    anc_a: Set[str],
    anc_b: Set[str],
) -> float:
    """Return MICA IC, using the module-level cache (canonical key ordering)."""
    key = (a, b) if a <= b else (b, a)
    cached = _mica_cache.get(key)
    if cached is not None:
        return cached
    common = anc_a & anc_b
    val = max((ic.get(t, 0.0) for t in common), default=0.0) if common else 0.0
    _mica_cache[key] = val
    return val


def _mica_ic(a: str, b: str, ic: Dict[str, float], ancestors: Dict[str, Set[str]]) -> float:
    """Information content of the Most Informative Common Ancestor."""
    anc_a = ancestors.get(a, {a})
    anc_b = ancestors.get(b, {b})
    return _mica_ic_cached(a, b, ic, anc_a, anc_b)


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


def _bma_fast(
    q_exp: List[str],
    t_exp: List[str],
    ic: Dict[str, float],
    ancestors: Dict[str, Set[str]],
    method: str,
) -> float:
    """Fast BMA using pre-expanded term lists and cached MICA IC values."""
    use_lin = method == "lin"
    scores = []
    for q in q_exp:
        anc_q = ancestors.get(q, {q})
        ic_q = ic.get(q, 0.0)
        best = 0.0
        for t in t_exp:
            anc_t = ancestors.get(t, {t})
            mica = _mica_ic_cached(q, t, ic, anc_q, anc_t)
            if use_lin:
                denom = ic_q + ic.get(t, 0.0)
                val = (2.0 * mica) / denom if denom > 0.0 else 0.0
            else:
                val = mica
            if val > best:
                best = val
        scores.append(best)
    return sum(scores) / len(scores) if scores else 0.0


def expand_hpos(terms: List[str], ancestors: Dict[str, Set[str]]) -> List[str]:
    """Expand HPO terms with their ancestors."""
    expanded: Set[str] = set()
    for t in terms:
        expanded.add(t)
        expanded |= ancestors.get(t, set())
    return list(expanded)


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

    q_exp = expand_hpos(query_hpos, ancestors)
    t_exp = expand_hpos(target_hpos, ancestors)

    q_to_t = _bma_fast(q_exp, t_exp, ic, ancestors, method)
    t_to_q = _bma_fast(t_exp, q_exp, ic, ancestors, method)
    return (q_to_t + t_to_q) / 2.0


def bma_similarity_fast(
    q_exp: List[str],
    t_exp: List[str],
    ic: Dict[str, float],
    ancestors: Dict[str, Set[str]],
    method: str = "lin",
) -> float:
    """
    Symmetric BMA using pre-expanded term lists.
    Call expand_hpos() on both sides before using this.
    """
    if not q_exp or not t_exp:
        return 0.0
    q_to_t = _bma_fast(q_exp, t_exp, ic, ancestors, method)
    t_to_q = _bma_fast(t_exp, q_exp, ic, ancestors, method)
    return (q_to_t + t_to_q) / 2.0
