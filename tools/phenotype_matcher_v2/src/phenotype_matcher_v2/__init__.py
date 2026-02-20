"""
Phenotype Matcher v2

Formal algorithm for mapping clinical phenotype descriptions to standardized
ontology identifiers: HPO, MONDO, OMIM, and Orphanet.

Quick start::

    from phenotype_matcher_v2 import PhenotypeMatcher, MatcherConfig

    cfg = MatcherConfig(hpo_path="ontology/hp.obo")
    m = PhenotypeMatcher(cfg)
    out = m.match("severe intellectual disability, no seizures")
    print(out.get_hpo_ids())
    print(out.get_excluded_phenotypes())
"""

from .matcher import PhenotypeMatcher
from .schemas import (
    DiseaseMatch,
    MatcherConfig,
    PhenotypeInput,
    PhenotypeMatch,
    PhenotypeOutput,
)

__version__ = "2.0.0"
__all__ = [
    "PhenotypeMatcher",
    "MatcherConfig",
    "PhenotypeInput",
    "PhenotypeOutput",
    "PhenotypeMatch",
    "DiseaseMatch",
]
