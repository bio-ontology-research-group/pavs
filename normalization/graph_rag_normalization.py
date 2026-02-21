import os
import json
import pickle
import logging
import numpy as np
import pandas as pd
from tqdm import tqdm
from pyhpo import Ontology
import pronto
import torch
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity
import requests
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class GraphRAG:
    def __init__(
        self,
        model_name="all-MiniLM-L6-v2",
        llm_model="google/gemini-2.5-flash",
        device=None,
        cache_dir="data/graph_rag_cache",
        debug=False,
    ):
        self.model_name = model_name
        self.llm_model = llm_model
        # Force CPU if device not specified - avoid old GPU compatibility issues
        self.device = device if device is not None else "cpu"
        self.debug = debug
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)

        self.model = None
        self.embeddings = None
        self.terms = []  # List of dicts: {'id': '...', 'name': '...', 'definition': '...'}
        self.term_id_map = {}  # id -> index in self.terms

        # Initialize
        self._load_ontology_data()
        self._load_or_compute_embeddings()

    def _load_ontology_data(self):
        logging.info("Loading HPO and MONDO ontologies...")

        # HPO
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
        except Exception as e:
            logging.error(f"Failed to load HPO: {e}")

        # MONDO
        mondo_path = "ontology/mondo.obo"
        if os.path.exists(mondo_path):
            try:
                ms = pronto.Ontology(mondo_path)
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
            except Exception as e:
                logging.error(f"Failed to load MONDO: {e}")

        self.term_id_map = {t["id"]: i for i, t in enumerate(self.terms)}

        # Build descendant sets for filtering by branch
        self._build_descendant_sets()

        logging.info(f"Loaded {len(self.terms)} terms.")

    def _build_descendant_sets(self):
        """Build sets of all descendants for key HPO branches."""
        # Use pyhpo for proper hierarchy traversal
        try:
            _ = Ontology()  # Ensure initialized

            # Get all descendants of Phenotypic abnormality (HP:0000118)
            phenotype_root = Ontology.get_hpo_object("HP:0000118")
            self.phenotype_descendants = set()
            if phenotype_root:
                self.phenotype_descendants = {phenotype_root.id}

                # Use children property (recursive)
                def add_descendants(term):
                    for child in term.children:
                        if child.id not in self.phenotype_descendants:
                            self.phenotype_descendants.add(child.id)
                            add_descendants(child)

                add_descendants(phenotype_root)

            # Get all descendants of Severity (HP:0012824)
            severity_root = Ontology.get_hpo_object("HP:0012824")
            self.severity_descendants = set()
            if severity_root:
                self.severity_descendants = {severity_root.id}

                def add_sev_descendants(term):
                    for child in term.children:
                        if child.id not in self.severity_descendants:
                            self.severity_descendants.add(child.id)
                            add_sev_descendants(child)

                add_sev_descendants(severity_root)

            logging.info(f"Phenotype branch: {len(self.phenotype_descendants)} terms")
            logging.info(f"Severity branch: {len(self.severity_descendants)} terms")
        except Exception as e:
            logging.warning(
                f"Could not build descendant sets: {e}. Using all HPO terms."
            )
            self.phenotype_descendants = set()
            self.severity_descendants = set()
            if severity_root:
                self.severity_descendants = {severity_root.id}
                for child in severity_root.all_children:
                    self.severity_descendants.add(child.id)

            logging.info(f"Phenotype branch: {len(self.phenotype_descendants)} terms")
            logging.info(f"Severity branch: {len(self.severity_descendants)} terms")
        except Exception as e:
            logging.warning(
                f"Could not build descendant sets: {e}. Using all HPO terms."
            )
            self.phenotype_descendants = set()
            self.severity_descendants = set()

    def _load_or_compute_embeddings(self):
        emb_path = os.path.join(
            self.cache_dir, f"embeddings_{self.model_name.replace('/', '_')}.pkl"
        )

        if os.path.exists(emb_path):
            logging.info(f"Loading cached embeddings from {emb_path}...")
            with open(emb_path, "rb") as f:
                data = pickle.load(f)
                self.embeddings = data["embeddings"]
                # Basic validation to ensure cache matches loaded terms
                if len(self.embeddings) != len(self.terms):
                    logging.warning("Cached embeddings count mismatch. Recomputing...")
                    self._compute_embeddings(emb_path)
        else:
            self._compute_embeddings(emb_path)

    def _compute_embeddings(self, save_path):
        try:
            self._run_embedding_computation(save_path, self.device)
        except Exception as e:
            if self.device == "cuda":
                logging.warning(
                    f"Embedding computation failed on CUDA: {e}. Falling back to CPU."
                )
                self.device = "cpu"
                self._run_embedding_computation(save_path, self.device)
            else:
                raise e

    def _run_embedding_computation(self, save_path, device):
        logging.info(f"Computing embeddings using {self.model_name} on {device}...")
        self.model = SentenceTransformer(self.model_name, device=device)

        # Prepare text for embedding: "Name: [NAME]. Definition: [DEF]. Synonyms: [SYNS]"
        texts = []
        for t in self.terms:
            txt = f"Name: {t['name']}."
            if t["definition"]:
                txt += f" Definition: {t['definition']}."
            if t["synonyms"]:
                txt += f" Synonyms: {', '.join(t['synonyms'][:5])}."
            texts.append(txt)

        self.embeddings = self.model.encode(
            texts,
            show_progress_bar=True,
            convert_to_numpy=True,
            normalize_embeddings=True,
        )

        with open(save_path, "wb") as f:
            pickle.dump({"embeddings": self.embeddings}, f)
        logging.info("Embeddings computed and saved.")

    def retrieve_candidates(self, query, k=5):
        if self.model is None:
            self.model = SentenceTransformer(self.model_name, device=self.device)

        query_emb = self.model.encode(
            [query], convert_to_numpy=True, normalize_embeddings=True
        )

        # Cosine similarity
        scores = cosine_similarity(query_emb, self.embeddings)[0]

        # Get top k indices
        top_k_indices = np.argsort(scores)[::-1][:k]

        results = []
        for idx in top_k_indices:
            term = self.terms[idx]
            results.append({"term": term, "score": float(scores[idx])})
        return results

    def get_graph_context(self, candidates):
        """
        Expands candidates with immediate graph neighborhood (parents).
        Returns a formatted string for the LLM.
        """
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
                context_str += (
                    f"  Parents: {', '.join(parents[:3])}\n"  # Limit to 3 parents
                )
            context_str += "\n"

        return context_str

    def retrieve_phenotype_candidates(self, query, k=5):
        """Retrieve candidates from Phenotypic abnormality branch only (HP:0000118)."""
        if self.model is None:
            self.model = SentenceTransformer(self.model_name, device=self.device)

        # Filter to only phenotype terms using descendant set
        if hasattr(self, "phenotype_descendants") and self.phenotype_descendants:
            phenotype_indices = [
                i
                for i, t in enumerate(self.terms)
                if t["source"] == "HPO" and t["id"] in self.phenotype_descendants
            ]
        else:
            # Fallback: use all HPO terms
            phenotype_indices = [
                i for i, t in enumerate(self.terms) if t["source"] == "HPO"
            ]

        if not phenotype_indices:
            return []

        query_emb = self.model.encode(
            [query], convert_to_numpy=True, normalize_embeddings=True
        )
        phenotype_embs = self.embeddings[phenotype_indices]

        scores = cosine_similarity(query_emb, phenotype_embs)[0]
        top_k_local = np.argsort(scores)[::-1][:k]

        results = []
        for local_idx in top_k_local:
            global_idx = phenotype_indices[local_idx]
            term = self.terms[global_idx]
            results.append({"term": term, "score": float(scores[local_idx])})
        return results

    def retrieve_severity_candidates(self, query, k=5):
        """Retrieve candidates from Severity modifier branch only (HP:0012824)."""
        if self.model is None:
            self.model = SentenceTransformer(self.model_name, device=self.device)

        # Filter to only severity terms using descendant set
        if hasattr(self, "severity_descendants") and self.severity_descendants:
            severity_indices = [
                i
                for i, t in enumerate(self.terms)
                if t["source"] == "HPO" and t["id"] in self.severity_descendants
            ]
        else:
            return []

        if not severity_indices:
            return []

        query_emb = self.model.encode(
            [query], convert_to_numpy=True, normalize_embeddings=True
        )
        severity_embs = self.embeddings[severity_indices]

        scores = cosine_similarity(query_emb, severity_embs)[0]
        top_k_local = np.argsort(scores)[::-1][:k]

        results = []
        for local_idx in top_k_local:
            global_idx = severity_indices[local_idx]
            term = self.terms[global_idx]
            results.append({"term": term, "score": float(scores[local_idx])})
        return results

    def retrieve_severity_candidates(self, query, k=5):
        """Retrieve candidates from Severity modifier branch only (HP:0012824)."""
        if self.model is None:
            self.model = SentenceTransformer(self.model_name, device=self.device)

        # Filter to only severity terms (descendants of HP:0012824)
        severity_indices = [
            i
            for i, t in enumerate(self.terms)
            if t["source"] == "HPO"
            and (t["id"] == "HP:0012824" or "HP:0012824" in str(t.get("parents", [])))
        ]

        if not severity_indices:
            return []

        query_emb = self.model.encode(
            [query], convert_to_numpy=True, normalize_embeddings=True
        )
        severity_embs = self.embeddings[severity_indices]

        scores = cosine_similarity(query_emb, severity_embs)[0]
        top_k_local = np.argsort(scores)[::-1][:k]

        results = []
        for local_idx in top_k_local:
            global_idx = severity_indices[local_idx]
            term = self.terms[global_idx]
            results.append({"term": term, "score": float(scores[local_idx])})
        return results

    def normalize_term(self, term, context_text="", k=5):
        """
        Full GraphRAG workflow for a single term.
        Returns potentially multiple phenotypes if term contains multiple features.
        """
        # Detect negation
        excluded = any(
            neg in term.lower()
            for neg in ["no ", "not ", "without ", "absent ", "normal "]
        )

        if self.debug:
            print(f"\n{'=' * 80}")
            print(f'DEBUG: Query term: "{term}"')
            print(f"DEBUG: Negation detected: {excluded}")
            print(f"{'=' * 80}")

        # 1. Retrieve phenotype candidates
        # If term contains "and"/"with", retrieve candidates BOTH ways
        has_conjunctions = any(
            conj in term.lower() for conj in [" and ", " with ", " or "]
        )

        if has_conjunctions:
            if self.debug:
                print(f"\nDEBUG: Detected conjunctions in term")

            # (a) Retrieve for entire phrase
            whole_phrase_candidates = self.retrieve_phenotype_candidates(term, k=k)

            # (b) Split into sub-phrases and retrieve for each
            sub_phrases = re.split(
                r"\s+and\s+|\s+with\s+|\s+or\s+", term, flags=re.IGNORECASE
            )
            if self.debug:
                print(f"DEBUG: Sub-phrases:")
                for i, sp in enumerate(sub_phrases, 1):
                    print(f'  {i}. "{sp.strip()}"')

            # Retrieve candidates for each sub-phrase
            sub_phrase_candidates = []
            seen_ids = set()
            for sub_phrase in sub_phrases:
                sp_clean = sub_phrase.strip().rstrip(".")
                if len(sp_clean) < 3:
                    continue
                sp_candidates = self.retrieve_phenotype_candidates(sp_clean, k=k)

                if self.debug:
                    print(f'\nDEBUG: Candidates for "{sp_clean}":')
                    for i, cand in enumerate(sp_candidates[:5], 1):
                        print(
                            f"    {i}. [{cand['score']:.3f}] {cand['term']['id']}: {cand['term']['name']}"
                        )

                # Add to combined list, avoiding duplicates
                for cand in sp_candidates:
                    if cand["term"]["id"] not in seen_ids:
                        seen_ids.add(cand["term"]["id"])
                        sub_phrase_candidates.append(cand)

            # Combine both sets: whole phrase + sub-phrases
            combined_candidates = []
            combined_seen = set()

            # Add whole phrase candidates first (keep context)
            for cand in whole_phrase_candidates:
                if cand["term"]["id"] not in combined_seen:
                    combined_seen.add(cand["term"]["id"])
                    combined_candidates.append(cand)

            # Add sub-phrase candidates (fill in missing terms)
            for cand in sub_phrase_candidates:
                if cand["term"]["id"] not in combined_seen:
                    combined_seen.add(cand["term"]["id"])
                    combined_candidates.append(cand)

            # Sort by score and take top k*3 (need more since we have multiple phenotypes)
            combined_candidates.sort(key=lambda x: x["score"], reverse=True)
            phenotype_candidates = combined_candidates[: k * 3]
        else:
            # Single phenotype - retrieve normally
            phenotype_candidates = self.retrieve_phenotype_candidates(term, k=k)

        if self.debug:
            print(
                f"\nDEBUG: Top-{len(phenotype_candidates)} Phenotype Candidates (combined):"
            )
            for i, cand in enumerate(phenotype_candidates, 1):
                t = cand["term"]
                print(f"  {i}. [{cand['score']:.3f}] {t['id']}: {t['name']}")
                if t.get("synonyms"):
                    print(f"     Synonyms: {', '.join(t['synonyms'][:3])}")

        phenotype_context = self.get_graph_context(phenotype_candidates)

        # 2. Retrieve severity candidates
        severity_candidates = self.retrieve_severity_candidates(term, k=min(k, 5))

        if self.debug and severity_candidates:
            print(f"\nDEBUG: Top-{min(k, 5)} Severity Candidates:")
            for i, cand in enumerate(severity_candidates, 1):
                t = cand["term"]
                print(f"  {i}. [{cand['score']:.3f}] {t['id']}: {t['name']}")

        severity_context = (
            self.get_graph_context(severity_candidates)
            if severity_candidates
            else "No severity modifiers found."
        )

        # 3. LLM Selection - LABELS ONLY, NO IDs
        # Detect if term contains multiple phenotypes
        has_multiple_indicator = any(
            word in term.lower() for word in [" and ", " with ", " or ", ","]
        )

        prompt = f"""You are an expert clinical geneticist. Map the clinical term "{term}" to standard ontology labels.

Patient Context: {context_text}

CRITICAL RULES:
1. You MUST select ONLY from the provided candidate NAMES/LABELS below
2. DO NOT output any IDs (HP:, MONDO:, OMIM:) - ONLY output the exact label/name from candidates
3. This term {"MAY CONTAIN MULTIPLE PHENOTYPES" if has_multiple_indicator else "likely describes a single phenotype"} (look for words like "and", "with")
4. Extract ALL phenotypes mentioned - a single term can map to 0, 1, or MULTIPLE phenotype labels

**Phenotype Candidates** (select from these names only):
{phenotype_context}

**Severity Candidates** (select from these names only):
{severity_context}

**Negation Detection**: The term {"CONTAINS" if excluded else "DOES NOT contain"} negation words.

Instructions:
1. Carefully read the term and identify EACH distinct phenotype mentioned
2. For EACH phenotype, select the best matching label from Phenotype Candidates
3. If the term contains "and", "with", or similar conjunctions, likely indicates MULTIPLE phenotypes
4. If severity detected (mild, moderate, severe, profound), select from Severity Candidates
5. Mark excluded: {str(excluded).lower()}

Examples:
- "scaly skin and hydronephrosis with reflux" → ["Scaling skin", "Hydronephrosis", "Vesicoureteral reflux"] (3 phenotypes)
- "severe intellectual disability" → ["Intellectual disability"] + severity: "Severe"
- "global developmental delay" → ["Global developmental delay"] (1 phenotype)

Output ONLY labels, NEVER IDs:
{{
  "phenotype_labels": ["exact label from candidates", ...],
  "severity_label": "exact label from candidates" or null,
  "excluded": {str(excluded).lower()},
  "disease_labels": [{{"mondo": "label", "omim": "label"}}, ...] or []
}}

FORBIDDEN: DO NOT output HP:XXXXXXX, MONDO:XXXXXXX, or OMIM:XXXXXX. Only output labels from candidates.
"""
        result = self._query_llm(prompt)

        if self.debug:
            import json

            print(f"\nDEBUG: LLM Response:")
            print(json.dumps(result, indent=2))
            print(f"{'=' * 80}\n")

        return result

    def _query_llm(self, prompt):
        api_key = os.getenv("OPENROUTER_API_KEY")
        if not api_key:
            logging.error("OPENROUTER_API_KEY not set")
            return {}

        try:
            resp = requests.post(
                "https://openrouter.ai/api/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {api_key}",
                    "Content-Type": "application/json",
                },
                json={
                    "model": self.llm_model,
                    "messages": [{"role": "user", "content": prompt}],
                    "response_format": {"type": "json_object"},
                },
                timeout=60,
            )
            if resp.status_code != 200:
                logging.error(f"LLM API Error: {resp.text}")
                return {}

            content = resp.json()["choices"][0]["message"]["content"]
            # Clean markdown code blocks if present
            content = re.sub(r"```json\n?|\n?```", "", content).strip()
            return json.loads(content)

        except Exception as e:
            logging.error(f"LLM Query Failed: {e}")
            return {}


if __name__ == "__main__":
    # Test run
    rag = GraphRAG()
    term = "Atrial septal defect"
    print(f"Testing normalization for: {term}")
    res = rag.normalize_term(term)
    print(json.dumps(res, indent=2))
