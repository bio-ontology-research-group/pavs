"""
Data schemas for Phenotype Matcher v2.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any


@dataclass
class MatcherConfig:
    """
    Configuration for PhenotypeMatcher v2.

    Embedding model presets:
        "sapbert"  → cambridgeltl/SapBERT-from-PubMedBERT-fulltext  [default]
        "fast"     → all-MiniLM-L6-v2
        "balanced" → all-mpnet-base-v2
        "accurate" → BAAI/bge-large-en-v1.5
        "medical"  → pritamdeka/S-PubMedBert-MS-MARCO

    LLM model presets:
        "deepseek" → deepseek/deepseek-chat          [default, free tier]
        "gpt4oss"  → openai/gpt-4o-mini              [cheap]
        "accurate" → anthropic/claude-sonnet-4-5
        "gemini"   → google/gemini-2.0-flash-exp:free
    """

    # Ontology paths
    hpo_path: str = "ontology/hp.obo"
    mondo_path: str = "ontology/mondo.obo"
    hpoa_path: str = "data/phenotype.hpoa"

    # Models
    embedding_model: str = "sapbert"
    llm_model: str = "deepseek"

    # Strategy 1 thresholds
    theta_abs: float = 0.70       # absolute floor for ANN scores
    delta: float = 0.05           # relative margin: keep within [S_max - delta, S_max]
    ann_top_k: int = 5            # cardinality cap for ANN candidates

    # Detection windows (in words)
    det_neg_win: int = 5
    det_mod_win: int = 4

    # Cache directory for embeddings
    cache_dir: str = "data/phenotype_matcher_v2_cache"

    # Device
    device: str = "cpu"

    # API
    api_key: Optional[str] = None

    # Debug
    debug: bool = False

    # Match-level cache: maximum number of distinct input strings to cache
    # (0 = disabled). Useful when the input TSV contains many identical strings.
    match_cache_size: int = 10_000


@dataclass
class PhenotypeMatch:
    """A matched HPO phenotype term."""

    hpo_id: str
    label: str
    excluded: bool = False
    severity_id: Optional[str] = None
    severity_label: Optional[str] = None
    confidence: float = 0.0
    matched_by: str = "unknown"


@dataclass
class DiseaseMatch:
    """A matched disease term (MONDO, OMIM, or Orphanet)."""

    mondo_id: Optional[str] = None
    mondo_label: Optional[str] = None
    omim_gene_ids: List[str] = field(default_factory=list)
    omim_phenotype_ids: List[str] = field(default_factory=list)
    omim_labels: List[str] = field(default_factory=list)
    orphanet_ids: List[str] = field(default_factory=list)
    confidence: float = 0.0
    matched_by: str = "unknown"


@dataclass
class PhenotypeInput:
    """Input for phenotype matching."""

    text: str
    context: Optional[str] = None
    gene_hint: Optional[str] = None
    split_by: Optional[str] = ","


@dataclass
class PhenotypeOutput:
    """Output from phenotype matching."""

    phenotypes: List[PhenotypeMatch] = field(default_factory=list)
    diseases: List[DiseaseMatch] = field(default_factory=list)
    raw_input: str = ""
    processing_metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "phenotypes": [
                {
                    "hpo_id": p.hpo_id,
                    "label": p.label,
                    "excluded": p.excluded,
                    "severity_id": p.severity_id,
                    "severity_label": p.severity_label,
                    "confidence": p.confidence,
                    "matched_by": p.matched_by,
                }
                for p in self.phenotypes
            ],
            "diseases": [
                {
                    "mondo_id": d.mondo_id,
                    "mondo_label": d.mondo_label,
                    "omim_gene_ids": d.omim_gene_ids,
                    "omim_phenotype_ids": d.omim_phenotype_ids,
                    "omim_labels": d.omim_labels,
                    "orphanet_ids": d.orphanet_ids,
                    "confidence": d.confidence,
                    "matched_by": d.matched_by,
                }
                for d in self.diseases
            ],
            "raw_input": self.raw_input,
            "processing_metadata": self.processing_metadata,
        }

    def get_present_phenotypes(self) -> List[PhenotypeMatch]:
        return [p for p in self.phenotypes if not p.excluded]

    def get_excluded_phenotypes(self) -> List[PhenotypeMatch]:
        return [p for p in self.phenotypes if p.excluded]

    def get_hpo_ids(self, excluded: bool = False) -> List[str]:
        if excluded:
            return [p.hpo_id for p in self.phenotypes if p.excluded]
        return [p.hpo_id for p in self.phenotypes if not p.excluded]

    def get_omim_ids(self) -> List[str]:
        omim_ids = []
        for d in self.diseases:
            omim_ids.extend(d.omim_gene_ids)
            omim_ids.extend(d.omim_phenotype_ids)
        return list(set(omim_ids))

    def get_mondo_ids(self) -> List[str]:
        return [d.mondo_id for d in self.diseases if d.mondo_id]

    def get_orphanet_ids(self) -> List[str]:
        ids = []
        for d in self.diseases:
            ids.extend(d.orphanet_ids)
        return list(set(ids))
