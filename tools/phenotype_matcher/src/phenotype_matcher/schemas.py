"""
Data schemas for the Phenotype Matcher tool.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any


@dataclass
class MatcherConfig:
    """
    Configuration for the PhenotypeMatcher.

    Attributes:
        embedding_model: Sentence transformer model for semantic embeddings
            - "fast": all-MiniLM-L6-v2 (384 dim, 80MB)
            - "balanced": all-mpnet-base-v2 (768 dim, 420MB) [default]
            - "accurate": BAAI/bge-large-en-v1.5 (1024 dim, 1.3GB)
            - "medical": pritamdeka/S-PubMedBert-MS-MARCO
            - "biobert": dmis-lab/biobert-base-cased-v1.2
            - "sapbert": cambridgeltl/SapBERT-from-PubMedBERT-fulltext (RECOMMENDED for medical ontology matching)
            - Custom model name from sentence-transformers

        llm_model: LLM model for semantic selection
            - "fast": openai/gpt-oss-120b [default]
            - "accurate": anthropic/claude-3.5-sonnet
            - "balanced": google/gemini-2.0-flash-exp:free
            - Custom model name (OpenRouter format)

        top_k_phenotype: Number of phenotype candidates to retrieve (default: 5)
        top_k_severity: Number of severity candidates to retrieve (default: 5)
        top_k_disease: Number of disease candidates to retrieve (default: 5)

        cache_dir: Directory for caching embeddings (default: data/graph_rag_cache)
        device: Compute device ("cpu" or "cuda", default: "cpu")

        api_key: OpenRouter API key (if not set, uses OPENROUTER_API_KEY env var)
        use_ner: Use LLM-based NER for term extraction (default: True, recommended)
        expand_acronyms: Expand medical acronyms during search (default: True)
        debug: Enable debug logging (default: False)
    """

    embedding_model: str = "balanced"
    llm_model: str = "fast"
    top_k_phenotype: int = 5
    top_k_severity: int = 5
    top_k_disease: int = 5
    cache_dir: str = "data/graph_rag_cache"
    device: str = "cpu"
    api_key: Optional[str] = None
    use_ner: bool = True
    expand_acronyms: bool = True
    debug: bool = False


@dataclass
class PhenotypeMatch:
    """
    A matched phenotype term.

    Attributes:
        hpo_id: HPO identifier (e.g., "HP:0001250")
        label: Human-readable phenotype name (e.g., "Seizure")
        excluded: Whether this phenotype is negated/excluded (default: False)
        severity_id: HPO severity modifier ID (e.g., "HP:0012828" for "Severe")
        severity_label: Human-readable severity (e.g., "Severe", "Mild")
        confidence: Semantic similarity score (0.0-1.0)
    """

    hpo_id: str
    label: str
    excluded: bool = False
    severity_id: Optional[str] = None
    severity_label: Optional[str] = None
    confidence: float = 0.0


@dataclass
class DiseaseMatch:
    """
    A matched disease term.

    Attributes:
        mondo_id: MONDO identifier (e.g., "MONDO:0001234")
        mondo_label: MONDO disease name
        omim_gene_ids: List of associated OMIM gene IDs (e.g., ["OMIM:123456"])
        omim_phenotype_ids: List of associated OMIM disease IDs (e.g., ["OMIM:654321"])
        omim_labels: List of OMIM disease names
        orphanet_ids: List of OrphaNet identifiers (if available via MONDO xrefs)
        confidence: Semantic similarity score (0.0-1.0)
    """

    mondo_id: Optional[str] = None
    mondo_label: Optional[str] = None
    omim_gene_ids: List[str] = field(default_factory=list)
    omim_phenotype_ids: List[str] = field(default_factory=list)
    omim_labels: List[str] = field(default_factory=list)
    orphanet_ids: List[str] = field(default_factory=list)
    confidence: float = 0.0


@dataclass
class PhenotypeInput:
    """
    Input for phenotype matching.

    Attributes:
        text: Clinical phenotype description (may contain multiple phenotypes)
              Examples:
              - "severe intellectual disability"
              - "seizures, feeding difficulties, and hypotonia"
              - "no cardiac abnormalities"

        context: Optional additional context (e.g., patient info, clinical notes)
                 Can help with disambiguation (e.g., "ASD" in cardiac vs neuro context)

        gene_hint: Optional gene symbol to aid disambiguation (e.g., "SCN1A", "MECP2")
                   Used as a SOFT cue only - helps choose between ambiguous matches
                   but NEVER overrides good semantic matches.
                   Example: "ASD" + gene_hint="NKX2-5" → favors Atrial Septal Defect
                           "ASD" + gene_hint="SHANK3" → favors Autism Spectrum Disorder

        split_by: Character(s) to split the text on for multiple phenotypes
                  Default: "," (comma)
                  Set to None to treat entire text as single phenotype
    """

    text: str
    context: Optional[str] = None
    gene_hint: Optional[str] = None
    split_by: Optional[str] = ","


@dataclass
class PhenotypeOutput:
    """
    Output from phenotype matching.

    Attributes:
        phenotypes: List of matched phenotype terms (both present and excluded)
        diseases: List of matched disease terms

        raw_input: Original input text
        processing_metadata: Additional processing information
            - terms_processed: Number of individual terms processed
            - llm_calls: Number of LLM API calls made
            - processing_time_seconds: Total processing time
    """

    phenotypes: List[PhenotypeMatch] = field(default_factory=list)
    diseases: List[DiseaseMatch] = field(default_factory=list)

    raw_input: str = ""
    processing_metadata: Dict[str, Any] = field(default_factory=dict)
    ner_extracted_terms: List[Dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert output to dictionary format."""
        return {
            "phenotypes": [
                {
                    "hpo_id": p.hpo_id,
                    "label": p.label,
                    "excluded": p.excluded,
                    "severity_id": p.severity_id,
                    "severity_label": p.severity_label,
                    "confidence": p.confidence,
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
                }
                for d in self.diseases
            ],
            "raw_input": self.raw_input,
            "processing_metadata": self.processing_metadata,
        }

    def get_present_phenotypes(self) -> List[PhenotypeMatch]:
        """Get only non-excluded phenotypes."""
        return [p for p in self.phenotypes if not p.excluded]

    def get_excluded_phenotypes(self) -> List[PhenotypeMatch]:
        """Get only excluded/negated phenotypes."""
        return [p for p in self.phenotypes if p.excluded]

    def get_hpo_ids(self, excluded: bool = False) -> List[str]:
        """
        Get list of HPO IDs.

        Args:
            excluded: If True, return excluded phenotypes. If False, return present phenotypes.
        """
        if excluded:
            return [p.hpo_id for p in self.phenotypes if p.excluded]
        else:
            return [p.hpo_id for p in self.phenotypes if not p.excluded]

    def get_omim_ids(self) -> List[str]:
        """Get all OMIM IDs (both gene and phenotype)."""
        omim_ids = []
        for disease in self.diseases:
            omim_ids.extend(disease.omim_gene_ids)
            omim_ids.extend(disease.omim_phenotype_ids)
        return list(set(omim_ids))  # Remove duplicates

    def get_mondo_ids(self) -> List[str]:
        """Get all MONDO IDs."""
        return [d.mondo_id for d in self.diseases if d.mondo_id]

    def get_orphanet_ids(self) -> List[str]:
        """Get all OrphaNet IDs."""
        orphanet_ids = []
        for disease in self.diseases:
            orphanet_ids.extend(disease.orphanet_ids)
        return list(set(orphanet_ids))  # Remove duplicates
