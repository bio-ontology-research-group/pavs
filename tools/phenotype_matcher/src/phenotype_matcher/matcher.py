"""
Main PhenotypeMatcher class implementing Graph RAG for phenotype normalization.
"""

import os
import json
import pickle
import logging
import time
import re
from typing import List, Dict, Any, Optional
import numpy as np
from pyhpo import Ontology
import pronto
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity
import requests

from .schemas import (
    PhenotypeInput,
    PhenotypeOutput,
    PhenotypeMatch,
    DiseaseMatch,
    MatcherConfig,
)
from . import acronyms
from . import ner


# Embedding model presets
EMBEDDING_MODELS = {
    "fast": "all-MiniLM-L6-v2",  # 384 dim, 80MB
    "balanced": "all-mpnet-base-v2",  # 768 dim, 420MB (default)
    "accurate": "BAAI/bge-large-en-v1.5",  # 1024 dim, 1.3GB
    "medical": "pritamdeka/S-PubMedBert-MS-MARCO",  # Medical domain
    "biobert": "dmis-lab/biobert-base-cased-v1.2",  # Biomedical
    "sapbert": "cambridgeltl/SapBERT-from-PubMedBERT-fulltext",  # Biomedical entity linking (RECOMMENDED for medical)
}

# LLM model presets
LLM_MODELS = {
    "fast": "openai/gpt-oss-120b",
    "accurate": "anthropic/claude-3.5-sonnet",
    "balanced": "google/gemini-2.0-flash-exp:free",
    "cheap": "deepseek/deepseek-chat",
}


class PhenotypeMatcher:
    """
    Phenotype matching tool using Graph RAG.

    Maps clinical phenotype descriptions to standardized ontology identifiers
    using semantic embeddings, graph structure, and LLM reasoning.

    Example:
        >>> from phenotype_matcher import PhenotypeMatcher, PhenotypeInput
        >>>
        >>> # Initialize matcher
        >>> matcher = PhenotypeMatcher()
        >>>
        >>> # Match phenotypes
        >>> input_data = PhenotypeInput(text="severe intellectual disability and seizures")
        >>> output = matcher.match(input_data)
        >>>
        >>> # Access results
        >>> for pheno in output.phenotypes:
        >>>     print(f"{pheno.hpo_id}: {pheno.label}")
        >>>     if pheno.severity_label:
        >>>         print(f"  Severity: {pheno.severity_label}")
    """

    def __init__(self, config: Optional[MatcherConfig] = None):
        """
        Initialize the phenotype matcher.

        Args:
            config: Configuration object. If None, uses defaults.
        """
        self.config = config or MatcherConfig()

        # Resolve model names from presets
        self.embedding_model_name = EMBEDDING_MODELS.get(
            self.config.embedding_model, self.config.embedding_model
        )
        self.llm_model_name = LLM_MODELS.get(
            self.config.llm_model, self.config.llm_model
        )

        # Setup logging
        log_level = logging.DEBUG if self.config.debug else logging.INFO
        logging.basicConfig(
            level=log_level, format="%(asctime)s - %(levelname)s - %(message)s"
        )
        self.logger = logging.getLogger(__name__)

        # Create cache directory
        os.makedirs(self.config.cache_dir, exist_ok=True)

        # Initialize components
        self.sentence_model = None
        self.embeddings = None
        self.terms = []  # List of ontology terms
        self.term_id_map = {}  # ID -> index mapping
        self.phenotype_descendants = set()  # HP:0000118 descendants
        self.severity_descendants = set()  # HP:0012824 descendants

        # API key
        self.api_key = self.config.api_key or os.getenv("OPENROUTER_API_KEY")

        # Load ontologies and embeddings
        self._initialize()

    def _initialize(self):
        """Load ontology data and embeddings."""
        self.logger.info("Initializing PhenotypeMatcher...")
        self._load_ontologies()
        self._build_descendant_sets()
        self._load_or_compute_embeddings()
        self.logger.info("Initialization complete.")

    def _load_ontologies(self):
        """Load HPO and MONDO ontologies."""
        self.logger.info("Loading ontologies...")

        # Load HPO
        try:
            _ = Ontology()  # Initialize pyhpo
            for term in Ontology:
                self.terms.append(
                    {
                        "id": term.id,
                        "name": term.name,
                        "definition": term.definition or "",
                        "synonyms": list(term.synonym)
                        if hasattr(term, "synonym")
                        else [],
                        "parents": [p.id for p in term.parents],
                        "children": [c.id for c in term.children],
                        "source": "HPO",
                    }
                )
            self.logger.info(f"Loaded {len(self.terms)} HPO terms")
        except Exception as e:
            self.logger.error(f"Failed to load HPO: {e}")

        # Load MONDO
        mondo_path = os.path.join("ontology", "mondo.obo")
        if os.path.exists(mondo_path):
            try:
                ms = pronto.Ontology(mondo_path)
                mondo_count = 0
                for term in ms.terms():
                    if term.id.startswith("MONDO:"):
                        self.terms.append(
                            {
                                "id": term.id,
                                "name": term.name,
                                "definition": str(term.definition)
                                if term.definition
                                else "",
                                "synonyms": [s.description for s in term.synonyms],
                                "parents": [
                                    p.id for p in term.superclasses(distance=1).to_set()
                                ],
                                "children": [
                                    c.id for c in term.subclasses(distance=1).to_set()
                                ],
                                "source": "MONDO",
                            }
                        )
                        mondo_count += 1
                self.logger.info(f"Loaded {mondo_count} MONDO terms")
            except Exception as e:
                self.logger.error(f"Failed to load MONDO: {e}")
        else:
            self.logger.warning(f"MONDO ontology not found at {mondo_path}")

        # Build ID map
        self.term_id_map = {t["id"]: i for i, t in enumerate(self.terms)}
        self.logger.info(f"Total terms loaded: {len(self.terms)}")

    def _build_descendant_sets(self):
        """Build descendant sets for phenotype and severity branches."""
        self.logger.info("Building descendant sets...")

        try:
            # Phenotypic abnormality branch (HP:0000118)
            phenotype_root = Ontology.get_hpo_object("HP:0000118")
            if phenotype_root:
                self.phenotype_descendants = {phenotype_root.id}

                def add_pheno_descendants(term):
                    for child in term.children:
                        if child.id not in self.phenotype_descendants:
                            self.phenotype_descendants.add(child.id)
                            add_pheno_descendants(child)

                add_pheno_descendants(phenotype_root)

            # Severity branch (HP:0012824)
            severity_root = Ontology.get_hpo_object("HP:0012824")
            if severity_root:
                self.severity_descendants = {severity_root.id}

                def add_sev_descendants(term):
                    for child in term.children:
                        if child.id not in self.severity_descendants:
                            self.severity_descendants.add(child.id)
                            add_sev_descendants(child)

                add_sev_descendants(severity_root)

            self.logger.info(
                f"Phenotype branch: {len(self.phenotype_descendants)} terms"
            )
            self.logger.info(f"Severity branch: {len(self.severity_descendants)} terms")

        except Exception as e:
            self.logger.warning(f"Could not build descendant sets: {e}")

    def _load_or_compute_embeddings(self):
        """Load cached embeddings or compute them."""
        cache_file = os.path.join(
            self.config.cache_dir,
            f"embeddings_{self.embedding_model_name.replace('/', '_')}.pkl",
        )

        if os.path.exists(cache_file):
            self.logger.info(f"Loading cached embeddings from {cache_file}")
            try:
                with open(cache_file, "rb") as f:
                    data = pickle.load(f)
                    self.embeddings = data["embeddings"]

                    # Validate cache
                    if len(self.embeddings) != len(self.terms):
                        self.logger.warning("Cache mismatch. Recomputing embeddings...")
                        self._compute_embeddings(cache_file)
            except Exception as e:
                self.logger.error(f"Failed to load cache: {e}. Recomputing...")
                self._compute_embeddings(cache_file)
        else:
            self._compute_embeddings(cache_file)

    def _compute_embeddings(self, cache_path: str):
        """Compute embeddings for all ontology terms."""
        self.logger.info(f"Computing embeddings using {self.embedding_model_name}...")

        try:
            self.sentence_model = SentenceTransformer(
                self.embedding_model_name, device=self.config.device
            )
        except Exception as e:
            if self.config.device == "cuda":
                self.logger.warning(f"CUDA failed: {e}. Falling back to CPU.")
                self.config.device = "cpu"
                self.sentence_model = SentenceTransformer(
                    self.embedding_model_name, device="cpu"
                )
            else:
                raise

        # Prepare text: "Name: [name]. Definition: [def]. Synonyms: [syns]"
        # CRITICAL: Include ALL synonyms for comprehensive matching
        texts = []
        for t in self.terms:
            txt = f"Name: {t['name']}."
            if t["definition"]:
                txt += f" Definition: {t['definition']}."
            if t["synonyms"]:
                # Include ALL synonyms (no filtering - 100% coverage)
                # Max synonyms in HPO: 28, so this is not a memory issue
                synonyms_text = ", ".join(t["synonyms"])
                txt += f" Synonyms: {synonyms_text}."
            texts.append(txt)

        # Compute embeddings
        self.embeddings = self.sentence_model.encode(
            texts,
            show_progress_bar=True,
            convert_to_numpy=True,
            normalize_embeddings=True,
        )

        # Save to cache
        with open(cache_path, "wb") as f:
            pickle.dump({"embeddings": self.embeddings}, f)

        self.logger.info("Embeddings computed and cached.")

    def match(self, input_data: PhenotypeInput) -> PhenotypeOutput:
        """
        Match phenotype description to ontology identifiers.

        Uses 3-step process:
        1. NER: Extract individual phenotype terms using LLM
        2. RAG: Retrieve candidates (with acronym expansion)
        3. Validation: LLM selects best matches with disambiguation

        Args:
            input_data: Input containing phenotype description

        Returns:
            PhenotypeOutput with matched phenotypes and diseases
        """
        start_time = time.time()

        # STEP 1: NER - Extract individual phenotype terms using LLM (if enabled)
        if self.config.use_ner and self.api_key:
            self.logger.info("Step 1: Extracting phenotype terms using NER...")
            extracted_terms = ner.extract_phenotype_terms_llm(
                input_data.text,
                self.api_key,
                self.llm_model_name,
            )
            llm_calls = 1  # NER call
        else:
            extracted_terms = []
            llm_calls = 0

        # Fallback to simple splitting if NER disabled or fails
        if not extracted_terms:
            self.logger.warning("NER failed, falling back to comma splitting")
            if input_data.split_by:
                terms = [
                    {
                        "term": t.strip(),
                        "modifiers": [],
                        "excluded": False,
                        "original_span": t.strip(),
                    }
                    for t in input_data.text.split(input_data.split_by)
                    if t.strip()
                ]
            else:
                terms = [
                    {
                        "term": input_data.text.strip(),
                        "modifiers": [],
                        "excluded": False,
                        "original_span": input_data.text.strip(),
                    }
                ]
            extracted_terms = terms

        self.logger.info(f"Extracted {len(extracted_terms)} phenotype terms")

        all_phenotypes = []
        all_diseases = []

        # Collect all term strings for context
        all_term_strings = [t.get("term", "") for t in extracted_terms]

        # STEP 2 & 3: For each extracted term, do RAG + validation
        for term_data in extracted_terms:
            term = term_data.get("term", "")
            modifiers = term_data.get("modifiers", []) or []  # Ensure it's a list
            excluded_by_ner = term_data.get("excluded", False)

            if not term:
                continue

            self.logger.debug(
                f"Processing term: '{term}', modifiers: {modifiers}, excluded: {excluded_by_ner}"
            )

            # Process the term
            result = self._normalize_term(
                term,
                input_data.context or "",
                gene_hint=input_data.gene_hint,
                other_terms=all_term_strings,
                ner_excluded=excluded_by_ner,
                ner_modifiers=modifiers,
            )
            llm_calls += 1

            # Convert LLM response to PhenotypeMatch objects
            phenotype_labels = result.get("phenotype_labels", [])
            severity_label = result.get("severity_label")
            excluded = result.get(
                "excluded", excluded_by_ner
            )  # Use NER exclusion if LLM didn't provide

            # Use NER modifiers if LLM didn't extract severity
            if not severity_label and modifiers:
                # Check if any NER modifier is a severity term
                severity_keywords = ["mild", "moderate", "severe", "profound"]
                for mod in modifiers:
                    if mod.lower() in severity_keywords:
                        severity_label = mod.capitalize()
                        break

            for pheno_label in phenotype_labels:
                # Look up HPO ID from label
                hpo_id = self._lookup_hpo_id(pheno_label)
                if hpo_id:
                    severity_id = None
                    if severity_label:
                        severity_id = self._lookup_severity_id(severity_label)

                    all_phenotypes.append(
                        PhenotypeMatch(
                            hpo_id=hpo_id,
                            label=pheno_label,
                            excluded=excluded,
                            severity_id=severity_id,
                            severity_label=severity_label if severity_id else None,
                            confidence=0.9,  # Placeholder - could use actual scores
                        )
                    )

            # Process disease labels
            disease_labels = result.get("disease_labels", [])
            for disease in disease_labels:
                mondo_label = disease.get("mondo")
                omim_label = disease.get("omim")

                disease_match = DiseaseMatch()

                if mondo_label:
                    mondo_id = self._lookup_mondo_id(mondo_label)
                    if mondo_id:
                        disease_match.mondo_id = mondo_id
                        disease_match.mondo_label = mondo_label

                if omim_label:
                    # OMIM lookup would go here
                    disease_match.omim_labels.append(omim_label)

                if disease_match.mondo_id or disease_match.omim_labels:
                    all_diseases.append(disease_match)

        # Create output
        output = PhenotypeOutput(
            phenotypes=all_phenotypes,
            diseases=all_diseases,
            raw_input=input_data.text,
            ner_extracted_terms=extracted_terms,  # Store NER output for debugging
            processing_metadata={
                "terms_processed": len(extracted_terms),
                "llm_calls": llm_calls,
                "processing_time_seconds": time.time() - start_time,
                "embedding_model": self.embedding_model_name,
                "llm_model": self.llm_model_name,
                "ner_used": self.config.use_ner,
                "acronym_expansion_used": self.config.expand_acronyms,
            },
        )

        return output

    def _normalize_term(
        self,
        term: str,
        context: str = "",
        gene_hint: Optional[str] = None,
        other_terms: Optional[List[str]] = None,
        ner_excluded: bool = False,
        ner_modifiers: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Normalize a single term using Graph RAG with acronym expansion.

        Args:
            term: The phenotype term to normalize
            context: Additional context
            gene_hint: Optional gene symbol for disambiguation (SOFT cue only)
            other_terms: Other phenotype terms in the description (for context)
            ner_excluded: Exclusion status from NER (if available)
            ner_modifiers: Modifiers extracted by NER (if available)

        Returns dict with phenotype_labels, severity_label, excluded, disease_labels.
        """
        # Use NER exclusion if provided, otherwise detect from text
        if ner_excluded:
            excluded = True
        else:
            excluded = any(
                neg in term.lower()
                for neg in ["no ", "not ", "without ", "absent ", "normal "]
            )

        # STEP 2A: Check for acronym expansion (if enabled)
        search_terms = [term]
        acronym_context = ""
        if self.config.expand_acronyms and acronyms.should_expand_for_search(term):
            expansions = acronyms.expand_acronym(term)
            if expansions:
                self.logger.info(f"Expanding acronym '{term}' to: {expansions}")
                search_terms = expansions
                acronym_context = f"\n**ACRONYM DISAMBIGUATION**: '{term}' could mean: {', '.join(expansions)}. Use context to choose the correct interpretation."

        # STEP 2B: Retrieve phenotype candidates (possibly from multiple expansions)
        all_pheno_candidates = []
        seen_ids = set()

        for search_term in search_terms:
            candidates = self._retrieve_phenotype_candidates(
                search_term, k=self.config.top_k_phenotype
            )
            # Merge candidates, avoiding duplicates
            for cand in candidates:
                term_id = cand["term"]["id"]
                if term_id not in seen_ids:
                    seen_ids.add(term_id)
                    all_pheno_candidates.append(cand)

        # Sort by score and take top-k
        all_pheno_candidates.sort(key=lambda x: x["score"], reverse=True)
        pheno_candidates = all_pheno_candidates[
            : self.config.top_k_phenotype * 2
        ]  # Allow more for disambiguation
        pheno_context = self._build_graph_context(pheno_candidates)

        # Retrieve severity candidates
        sev_candidates = self._retrieve_severity_candidates(
            term, k=self.config.top_k_severity
        )
        sev_context = (
            self._build_graph_context(sev_candidates)
            if sev_candidates
            else "No severity modifiers found."
        )

        # Build context information for disambiguation
        context_info = []
        if context:
            context_info.append(f"Patient Context: {context}")

        if other_terms and len(other_terms) > 1:
            # Show other phenotypes in the description for disambiguation
            other_phenos = [t for t in other_terms if t != term]
            if other_phenos:
                context_info.append(
                    f"Other phenotypes in description: {', '.join(other_phenos[:5])}"
                )

        if gene_hint:
            context_info.append(
                f"Gene hint: {gene_hint} (use ONLY as SOFT cue for disambiguation, NEVER override good matches)"
            )

        context_block = (
            "\n".join(context_info)
            if context_info
            else "No additional context provided."
        )

        # Query LLM
        has_multiple = any(
            word in term.lower() for word in [" and ", " with ", " or ", ","]
        )

        prompt = f"""You are an expert clinical geneticist. Map the clinical term "{term}" to standard ontology labels.

{context_block}
{acronym_context}

CRITICAL RULES:
1. You MUST select ONLY from the provided candidate NAMES/LABELS below
2. DO NOT output any IDs (HP:, MONDO:, OMIM:) - ONLY output the exact label/name from candidates
3. This term {"MAY CONTAIN MULTIPLE PHENOTYPES" if has_multiple else "likely describes a single phenotype"}
4. Extract ALL phenotypes mentioned - a single term can map to 0, 1, or MULTIPLE phenotype labels
5. If gene hint provided, use it ONLY for disambiguation between similar candidates, NEVER to override clear semantic matches

SPECIFICITY RULES (CRITICAL):
- DO NOT select terms that are MORE SPECIFIC than what is described
  Example: "hypotonia" should NOT match "Episodic generalized hypotonia" (too specific)
  Example: "seizures" should NOT match "Absence seizures" unless "absence" is mentioned
- It is ACCEPTABLE to select terms that are MORE GENERAL (less precise)
  Example: "absence seizures" CAN match "Seizure" (more general, still correct)
- Only select specific terms if the specificity is explicitly mentioned in the clinical text

**Phenotype Candidates** (select from these names only):
{pheno_context}

**Severity Candidates** (select from these names only):
{sev_context}

**Negation Detection**: The term {"CONTAINS" if excluded else "DOES NOT contain"} negation words.

Instructions:
1. Carefully read the term and identify EACH distinct phenotype mentioned
2. For EACH phenotype, select the best matching label from Phenotype Candidates
3. Prefer general terms over overly specific ones (unless specificity is in the text)
4. If ambiguous (e.g., "ASD"), use the other phenotypes and gene hint to disambiguate
5. If severity detected (mild, moderate, severe, profound), select from Severity Candidates
6. Mark excluded: {str(excluded).lower()}

Output ONLY labels, NEVER IDs:
{{
  "phenotype_labels": ["exact label from candidates", ...],
  "severity_label": "exact label from candidates" or null,
  "excluded": {str(excluded).lower()},
  "disease_labels": [{{"mondo": "label", "omim": "label"}}, ...] or []
}}

FORBIDDEN: DO NOT output HP:XXXXXXX, MONDO:XXXXXXX, or OMIM:XXXXXX. Only output labels.
"""

        return self._query_llm(prompt)

    def _retrieve_phenotype_candidates(self, query: str, k: int = 5) -> List[Dict]:
        """Retrieve top-k phenotype candidates from HP:0000118 branch."""
        if self.sentence_model is None:
            self.sentence_model = SentenceTransformer(
                self.embedding_model_name, device=self.config.device
            )

        # Filter to phenotype branch
        pheno_indices = [
            i
            for i, t in enumerate(self.terms)
            if t["source"] == "HPO" and t["id"] in self.phenotype_descendants
        ]

        if not pheno_indices:
            return []

        # Encode query
        query_emb = self.sentence_model.encode(
            [query], convert_to_numpy=True, normalize_embeddings=True
        )

        # Compute similarities
        pheno_embs = self.embeddings[pheno_indices]
        scores = cosine_similarity(query_emb, pheno_embs)[0]
        top_k_local = np.argsort(scores)[::-1][:k]

        # Build results
        results = []
        for local_idx in top_k_local:
            global_idx = pheno_indices[local_idx]
            term = self.terms[global_idx]
            results.append({"term": term, "score": float(scores[local_idx])})

        return results

    def _retrieve_severity_candidates(self, query: str, k: int = 5) -> List[Dict]:
        """Retrieve top-k severity candidates from HP:0012824 branch."""
        if self.sentence_model is None:
            self.sentence_model = SentenceTransformer(
                self.embedding_model_name, device=self.config.device
            )

        # Filter to severity branch
        sev_indices = [
            i
            for i, t in enumerate(self.terms)
            if t["source"] == "HPO" and t["id"] in self.severity_descendants
        ]

        if not sev_indices:
            return []

        # Encode query
        query_emb = self.sentence_model.encode(
            [query], convert_to_numpy=True, normalize_embeddings=True
        )

        # Compute similarities
        sev_embs = self.embeddings[sev_indices]
        scores = cosine_similarity(query_emb, sev_embs)[0]
        top_k_local = np.argsort(scores)[::-1][:k]

        # Build results
        results = []
        for local_idx in top_k_local:
            global_idx = sev_indices[local_idx]
            term = self.terms[global_idx]
            results.append({"term": term, "score": float(scores[local_idx])})

        return results

    def _build_graph_context(self, candidates: List[Dict]) -> str:
        """Build graph context string from candidates."""
        context_str = ""
        seen_ids = set()

        for cand in candidates:
            t = cand["term"]
            if t["id"] in seen_ids:
                continue
            seen_ids.add(t["id"])

            context_str += (
                f"- **{t['name']}** ({t['id']}) [Score: {cand['score']:.4f}]\n"
            )
            if t["definition"]:
                context_str += f"  Definition: {t['definition']}\n"

            # Add parents for context
            parents = []
            for pid in t["parents"]:
                if pid in self.term_id_map:
                    p = self.terms[self.term_id_map[pid]]
                    parents.append(f"{p['name']} ({p['id']})")

            if parents:
                context_str += f"  Parents: {', '.join(parents[:3])}\n"
            context_str += "\n"

        return context_str

    def _query_llm(self, prompt: str) -> Dict[str, Any]:
        """Query LLM via OpenRouter API."""
        if not self.api_key:
            self.logger.error("OpenRouter API key not set")
            return {}

        try:
            resp = requests.post(
                "https://openrouter.ai/api/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {self.api_key}",
                    "Content-Type": "application/json",
                },
                json={
                    "model": self.llm_model_name,
                    "messages": [{"role": "user", "content": prompt}],
                    "response_format": {"type": "json_object"},
                },
                timeout=60,
            )

            if resp.status_code != 200:
                self.logger.error(f"LLM API Error: {resp.text}")
                return {}

            content = resp.json()["choices"][0]["message"]["content"]
            # Clean markdown code blocks
            content = re.sub(r"```json\n?|\n?```", "", content).strip()
            return json.loads(content)

        except Exception as e:
            self.logger.error(f"LLM query failed: {e}")
            return {}

    def _lookup_hpo_id(self, label: str) -> Optional[str]:
        """Look up HPO ID from label."""
        label_lower = label.lower().strip()

        # Exact match
        for term in self.terms:
            if term["source"] == "HPO" and term["name"].lower() == label_lower:
                return term["id"]

        # Synonym match
        for term in self.terms:
            if term["source"] == "HPO":
                for syn in term.get("synonyms", []):
                    if syn.lower() == label_lower:
                        return term["id"]

        return None

    def _lookup_severity_id(self, label: str) -> Optional[str]:
        """Look up severity HPO ID from label."""
        label_lower = label.lower().strip()

        for term in self.terms:
            if term["source"] == "HPO" and term["id"] in self.severity_descendants:
                if term["name"].lower() == label_lower:
                    return term["id"]

        return None

    def _lookup_mondo_id(self, label: str) -> Optional[str]:
        """Look up MONDO ID from label."""
        label_lower = label.lower().strip()

        for term in self.terms:
            if term["source"] == "MONDO" and term["name"].lower() == label_lower:
                return term["id"]

        return None
