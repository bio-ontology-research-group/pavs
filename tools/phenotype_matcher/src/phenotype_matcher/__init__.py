"""
Phenotype Matcher - A standalone tool for matching clinical phenotype descriptions to ontology identifiers.

This tool uses Graph RAG (Retrieval-Augmented Generation) with semantic embeddings and
ontology graph structure to map free-text phenotype descriptions to standardized identifiers:
- HPO (Human Phenotype Ontology)
- MONDO (Disease Ontology)
- OMIM (Online Mendelian Inheritance in Man)
- OrphaNet (via MONDO cross-references)

The tool supports:
- Multiple phenotype mentions in a single description
- Negation detection (excluded phenotypes)
- Severity modifier extraction
- High-accuracy semantic matching with contextual understanding
- Disambiguation of ambiguous terms using context and gene hints
"""

from .matcher import PhenotypeMatcher
from .schemas import (
    PhenotypeInput,
    PhenotypeOutput,
    PhenotypeMatch,
    DiseaseMatch,
    MatcherConfig,
)
from . import test_cases
from . import disambiguation

__version__ = "1.0.0"
__all__ = [
    "PhenotypeMatcher",
    "PhenotypeInput",
    "PhenotypeOutput",
    "PhenotypeMatch",
    "DiseaseMatch",
    "MatcherConfig",
    "test_cases",
    "disambiguation",
]
