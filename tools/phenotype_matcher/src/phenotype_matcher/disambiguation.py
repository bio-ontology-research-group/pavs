"""
Disambiguation logic for ambiguous phenotype terms.

Handles cases where a term could match multiple distinct phenotypes
using contextual clues from co-occurring phenotypes and optional gene hints.
"""

from typing import List, Dict, Any, Optional, Set
import logging
from pyhpo import Ontology


logger = logging.getLogger(__name__)


def get_top_level_category(hpo_id: str) -> Optional[str]:
    """
    Get the top-level phenotypic abnormality category for an HPO term.

    Returns the first-level child of HP:0000118 (Phenotypic abnormality)
    that this term belongs to.

    Args:
        hpo_id: HPO identifier (e.g., "HP:0001250")

    Returns:
        Top-level category ID or None if not found
    """
    try:
        term = Ontology.get_hpo_object(hpo_id)
        if not term:
            return None

        # Get all ancestors
        ancestors = set()
        current_terms = [term]
        visited = set()

        while current_terms:
            current = current_terms.pop()
            if current.id in visited:
                continue
            visited.add(current.id)
            ancestors.add(current.id)
            current_terms.extend(current.parents)

        # Find direct children of HP:0000118 (Phenotypic abnormality)
        root = Ontology.get_hpo_object("HP:0000118")
        if root:
            for child in root.children:
                if child.id in ancestors:
                    return child.id

        return None
    except Exception as e:
        logger.warning(f"Could not get top-level category for {hpo_id}: {e}")
        return None


def get_category_name(hpo_id: str) -> str:
    """Get human-readable name for a category."""
    try:
        term = Ontology.get_hpo_object(hpo_id)
        return term.name if term else hpo_id
    except:
        return hpo_id


def detect_ambiguous_candidates(
    candidates: List[Dict[str, Any]],
    similarity_threshold: float = 0.3,
    category_diversity_threshold: int = 2,
) -> bool:
    """
    Detect if candidate matches are ambiguous.

    Ambiguity indicators:
    1. Candidates from multiple top-level HPO categories
    2. Large spread in semantic similarity scores

    Args:
        candidates: List of candidate dictionaries with 'term' and 'score'
        similarity_threshold: Min difference in scores to indicate ambiguity
        category_diversity_threshold: Min number of different categories

    Returns:
        True if candidates appear ambiguous
    """
    if len(candidates) < 2:
        return False

    # Check category diversity
    categories = set()
    for cand in candidates:
        term_id = cand["term"]["id"]
        cat = get_top_level_category(term_id)
        if cat:
            categories.add(cat)

    if len(categories) >= category_diversity_threshold:
        logger.info(
            f"Ambiguity detected: {len(categories)} different top-level categories"
        )
        return True

    # Check similarity score spread
    scores = [c["score"] for c in candidates]
    if len(scores) >= 2:
        score_spread = max(scores) - min(scores)
        if score_spread < similarity_threshold:
            logger.info(f"Ambiguity detected: Low score spread ({score_spread:.3f})")
            return True

    return False


def build_context_from_phenotypes(phenotype_labels: List[str]) -> str:
    """
    Build contextual description from co-occurring phenotypes.

    Args:
        phenotype_labels: List of other phenotypes in the description

    Returns:
        Contextual description string
    """
    if not phenotype_labels:
        return ""

    # Remove duplicates and limit length
    unique_labels = list(set(phenotype_labels))[:10]
    return f"co-occurring phenotypes: {', '.join(unique_labels)}"


def get_gene_associated_phenotypes(gene_symbol: str) -> Set[str]:
    """
    Get HPO terms associated with a gene.

    This is a placeholder - in a real implementation, this would query
    a gene-phenotype database (e.g., from HPO annotations).

    Args:
        gene_symbol: Gene symbol (e.g., "SCN1A")

    Returns:
        Set of HPO IDs associated with the gene
    """
    # Placeholder - would query actual gene-phenotype associations
    # For now, return empty set
    # TODO: Integrate with HPO gene annotations or phenotype.hpoa
    logger.debug(f"Gene-phenotype lookup for {gene_symbol} not yet implemented")
    return set()


def score_candidate_with_context(
    candidate: Dict[str, Any],
    context_phenotypes: List[str],
    gene_hint: Optional[str] = None,
) -> float:
    """
    Score a candidate match using contextual information.

    Args:
        candidate: Candidate dictionary with 'term' and 'score'
        context_phenotypes: Other phenotypes in the description
        gene_hint: Optional gene symbol for disambiguation

    Returns:
        Adjusted score (0.0 to 1.0)
    """
    base_score = candidate["score"]
    term_id = candidate["term"]["id"]

    # Get top-level category for this candidate
    candidate_category = get_top_level_category(term_id)

    # Score boost based on context phenotypes
    context_boost = 0.0
    if context_phenotypes and candidate_category:
        # Check if context phenotypes share the same category
        matching_categories = 0
        for pheno_label in context_phenotypes:
            # This is simplified - in practice, we'd need to look up each phenotype's category
            # For now, we use category name matching
            cat_name = get_category_name(candidate_category).lower()
            if any(word in pheno_label.lower() for word in cat_name.split()):
                matching_categories += 1

        if matching_categories > 0:
            context_boost = min(0.2, matching_categories * 0.05)
            logger.debug(f"Context boost for {term_id}: +{context_boost:.3f}")

    # Score boost based on gene hint
    gene_boost = 0.0
    if gene_hint:
        gene_phenotypes = get_gene_associated_phenotypes(gene_hint)
        if term_id in gene_phenotypes:
            gene_boost = 0.15
            logger.debug(
                f"Gene boost for {term_id} (gene={gene_hint}): +{gene_boost:.3f}"
            )

    # Combine scores (never exceed 1.0)
    adjusted_score = min(1.0, base_score + context_boost + gene_boost)

    return adjusted_score


def disambiguate_candidates(
    candidates: List[Dict[str, Any]],
    query_text: str,
    context_phenotypes: Optional[List[str]] = None,
    gene_hint: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """
    Disambiguate ambiguous candidate matches using context.

    Args:
        candidates: List of candidate dictionaries with 'term' and 'score'
        query_text: Original query text
        context_phenotypes: Other phenotypes in the description (for context)
        gene_hint: Optional gene symbol for disambiguation

    Returns:
        Re-ranked candidates with adjusted scores
    """
    if not candidates:
        return candidates

    # Check if disambiguation is needed
    if not detect_ambiguous_candidates(candidates):
        logger.debug("No ambiguity detected, using original ranking")
        return candidates

    logger.info(f"Disambiguating {len(candidates)} candidates for '{query_text}'")

    # Log categories for debugging
    categories = {}
    for cand in candidates:
        term_id = cand["term"]["id"]
        cat = get_top_level_category(term_id)
        if cat:
            cat_name = get_category_name(cat)
            categories[term_id] = cat_name
            logger.debug(
                f"  {term_id} ({cand['term']['name']}): {cat_name} [score={cand['score']:.3f}]"
            )

    # Rescore with context
    context_phenos = context_phenotypes or []
    reranked = []

    for cand in candidates:
        adjusted_score = score_candidate_with_context(
            cand,
            context_phenos,
            gene_hint,
        )

        # Create new candidate dict with adjusted score
        new_cand = cand.copy()
        new_cand["original_score"] = cand["score"]
        new_cand["adjusted_score"] = adjusted_score
        new_cand["score"] = adjusted_score  # Update main score

        if gene_hint:
            new_cand["gene_hint_used"] = gene_hint

        reranked.append(new_cand)

    # Re-sort by adjusted score
    reranked.sort(key=lambda x: x["adjusted_score"], reverse=True)

    # Log final ranking
    logger.info("After disambiguation:")
    for i, cand in enumerate(reranked[:5], 1):
        term = cand["term"]
        orig = cand.get("original_score", 0)
        adj = cand.get("adjusted_score", 0)
        cat = categories.get(term["id"], "Unknown")
        logger.info(
            f"  {i}. {term['id']} ({term['name']}): {orig:.3f} → {adj:.3f} [{cat}]"
        )

    return reranked


def build_disambiguation_prompt_context(
    candidates: List[Dict[str, Any]],
    query_text: str,
    context_phenotypes: Optional[List[str]] = None,
    gene_hint: Optional[str] = None,
) -> str:
    """
    Build additional context for LLM prompt when disambiguation is needed.

    Args:
        candidates: List of candidate dictionaries
        query_text: Original query text
        context_phenotypes: Other phenotypes in description
        gene_hint: Optional gene symbol

    Returns:
        Additional prompt text for disambiguation
    """
    if not detect_ambiguous_candidates(candidates):
        return ""

    # Build disambiguation guidance
    parts = ["\n**DISAMBIGUATION REQUIRED**"]
    parts.append(
        f"The term '{query_text}' is ambiguous and could match multiple phenotypes:"
    )

    # List categories
    categories = {}
    for cand in candidates[:5]:
        term_id = cand["term"]["id"]
        term_name = cand["term"]["name"]
        cat = get_top_level_category(term_id)
        if cat:
            cat_name = get_category_name(cat)
            categories.setdefault(cat_name, []).append(f"{term_name} ({term_id})")

    for cat_name, terms in categories.items():
        parts.append(f"  - {cat_name}: {', '.join(terms)}")

    # Add context clues
    if context_phenotypes:
        parts.append(
            f"\n**Context from other phenotypes:** {', '.join(context_phenotypes[:5])}"
        )
        parts.append(
            "Use these co-occurring phenotypes to determine which category is most appropriate."
        )

    if gene_hint:
        parts.append(f"\n**Gene hint:** {gene_hint}")
        parts.append(
            f"Consider phenotypes typically associated with {gene_hint} mutations."
        )
        parts.append(
            "NOTE: The gene hint is only for disambiguation - do NOT override good matches."
        )

    parts.append("\nSelect the phenotype that best fits the overall clinical context.")

    return "\n".join(parts)
