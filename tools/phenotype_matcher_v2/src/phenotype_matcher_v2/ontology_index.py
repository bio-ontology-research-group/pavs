"""
OntologyIndex: loads HPO / MONDO / OMIM / Orphanet, builds all look-up structures.

Loaded structures
-----------------
primary_label   : Dict[str, str]               term_id  → canonical name
exact_map       : Dict[str, List[str]]         label.lower() → [term_ids]
stemmed_map     : Dict[FrozenSet, List[str]]   frozenset(lemmas) → [term_ids]
modifier_map    : Dict[str, str]               modifier_label.lower() → modifier_term_id
phenotype_ids   : Set[str]                     descendants of HP:0000118
modifier_ids    : Set[str]                     descendants of HP:0012824
ac_automaton    : ahocorasick.Automaton        over all keys of exact_map
embeddings      : np.ndarray                   per-label SapBERT embeddings
embedding_labels: List[str]                    label string per row of embeddings
mondo_to_orphanet: Dict[str, List[str]]        MONDO id → [Orphanet ids]
"""

import os
import pickle
import re
from collections import defaultdict
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

import numpy as np

try:
    import pronto
except ImportError:
    pronto = None  # type: ignore

try:
    import ahocorasick
except ImportError:
    ahocorasick = None  # type: ignore

from .normalization_utils import get_stemmed_tokens


# HPO root IDs
_HPO_PHENOTYPE_ROOT = "HP:0000118"   # Phenotypic abnormality
_HPO_SEVERITY_ROOT  = "HP:0012824"  # Severity

# Prefixes for disease entries in phenotype.hpoa
_DISEASE_PREFIXES = ("OMIM:", "ORPHA:")


class OntologyIndex:
    """Load and index HPO, MONDO, OMIM, and Orphanet ontologies."""

    def __init__(
        self,
        hpo_path: str = "ontology/hp.obo",
        mondo_path: str = "ontology/mondo.obo",
        hpoa_path: str = "data/phenotype.hpoa",
        embedding_model_name: str = "cambridgeltl/SapBERT-from-PubMedBERT-fulltext",
        cache_dir: str = "data/phenotype_matcher_v2_cache",
        device: str = "cpu",
        debug: bool = False,
    ):
        self.hpo_path = hpo_path
        self.mondo_path = mondo_path
        self.hpoa_path = hpoa_path
        self.embedding_model_name = embedding_model_name
        self.cache_dir = cache_dir
        self.device = device
        self.debug = debug

        # Structures populated by _build()
        self.primary_label: Dict[str, str] = {}
        self.exact_map: Dict[str, List[str]] = defaultdict(list)
        self.stemmed_map: Dict[FrozenSet, List[str]] = defaultdict(list)
        self.modifier_map: Dict[str, str] = {}
        self.phenotype_ids: Set[str] = set()
        self.modifier_ids: Set[str] = set()
        self.mondo_to_orphanet: Dict[str, List[str]] = defaultdict(list)
        self.ac_automaton = None
        self.embeddings: Optional[np.ndarray] = None
        self.embedding_labels: List[str] = []
        self.hpo_ancestors: Dict[str, Set[str]] = {}
        self.phenotype_word_set: Set[str] = set()

        os.makedirs(cache_dir, exist_ok=True)
        self._build()

    # ------------------------------------------------------------------
    # Build pipeline
    # ------------------------------------------------------------------

    def _build(self) -> None:
        if self.debug:
            print("[OntologyIndex] Loading HPO …")
        self._load_hpo()

        if self.debug:
            print("[OntologyIndex] Generating abnormality synonyms …")
        self._generate_abnormality_synonyms()

        if self.debug:
            print("[OntologyIndex] Loading MONDO …")
        self._load_mondo()

        if self.debug:
            print("[OntologyIndex] Loading OMIM/Orphanet from phenotype.hpoa …")
        self._load_hpoa()

        if self.debug:
            print("[OntologyIndex] Building stemmed map …")
        self._build_stemmed_map()

        if self.debug:
            print("[OntologyIndex] Building Aho-Corasick automaton …")
        self._build_ac_automaton()

        if self.debug:
            print("[OntologyIndex] Building embeddings …")
        self._build_embeddings()

        if self.debug:
            print("[OntologyIndex] Building phenotype word set …")
        self._build_phenotype_word_set()

        if self.debug:
            print("[OntologyIndex] Ready.")

    # ------------------------------------------------------------------
    # HPO loading
    # ------------------------------------------------------------------

    def _load_hpo(self) -> None:
        if pronto is None:
            raise ImportError("pronto is required: pip install pronto")

        ont = pronto.Ontology(self.hpo_path)

        # Collect descendant sets for phenotype root and severity root
        try:
            pheno_root = ont[_HPO_PHENOTYPE_ROOT]
            self.phenotype_ids = {t.id for t in pheno_root.subclasses()}
        except (KeyError, AttributeError):
            self.phenotype_ids = set()

        try:
            sev_root = ont[_HPO_SEVERITY_ROOT]
            self.modifier_ids = {t.id for t in sev_root.subclasses()}
        except (KeyError, AttributeError):
            self.modifier_ids = set()

        for term in ont.terms():
            if not term.id.startswith("HP:"):
                continue
            name = term.name or ""
            if not name:
                continue

            self._add_label(term.id, name)

            for syn in term.synonyms:
                if syn.description:
                    self._add_label(term.id, syn.description)

            # Build modifier_map for severity descendants
            if term.id in self.modifier_ids:
                self.modifier_map[name.lower()] = term.id
                for syn in term.synonyms:
                    if syn.description:
                        self.modifier_map[syn.description.lower()] = term.id

        # Fix 2a: Build ancestor map for parent subsumption in _build_output()
        for term in ont.terms():
            if not term.id.startswith("HP:"):
                continue
            # superclasses() returns the term itself + all ancestors; exclude self
            self.hpo_ancestors[term.id] = {
                t.id for t in term.superclasses() if t.id != term.id
            }

    # ------------------------------------------------------------------
    # Synthetic synonyms
    # ------------------------------------------------------------------

    def _generate_abnormality_synonyms(self) -> None:
        """
        1. For every HPO label matching "Abnormality of [the] X" — including
        synonyms, not just primary labels — add "X abnormality" and
        "X abnormalities" as synthetic synonyms.

        2. For every label containing "ganglia" or "ganglion", generate the
        opposite form (e.g. "basal ganglion" from "basal ganglia") as a
        synthetic synonym.

        NOTE: delete the embeddings cache after this change so the new
        synonym labels are included in the ANN index.
        """
        # --- Part 1: Abnormality of X ---
        pattern = re.compile(r"^abnormality of (?:the )?(.+)$", re.I)
        seen_pairs: set = set()

        for label_lower, term_ids in list(self.exact_map.items()):
            m = pattern.match(label_lower)
            if not m:
                continue
            noun = m.group(1).strip()

            noun_forms = {noun}
            words = noun.split()
            last = words[-1] if words else ""
            if last.endswith("s") and not last.endswith("ss") and len(last) > 3:
                singular = last[:-1]
                noun_forms.add(
                    " ".join(words[:-1] + [singular]) if len(words) > 1 else singular
                )

            for term_id in term_ids:
                if not term_id.startswith("HP:"):
                    continue
                for n in noun_forms:
                    for suffix in ("abnormality", "abnormalities"):
                        syn = f"{n} {suffix}"
                        pair = (syn, term_id)
                        if pair in seen_pairs:
                            continue
                        seen_pairs.add(pair)
                        self._add_label(term_id, syn)

        # --- Part 2: Ganglion/Ganglia expansion ---
        for label_lower, term_ids in list(self.exact_map.items()):
            if "ganglia" in label_lower or "ganglion" in label_lower:
                for term_id in term_ids:
                    # Switch ganglia <-> ganglion
                    syn = label_lower.replace("ganglia", "@@@").replace("ganglion", "ganglia").replace("@@@", "ganglion")
                    if syn != label_lower:
                        self._add_label(term_id, syn)

    # ------------------------------------------------------------------
    # MONDO loading
    # ------------------------------------------------------------------

    def _load_mondo(self) -> None:
        if not os.path.exists(self.mondo_path):
            if self.debug:
                print(f"[OntologyIndex] MONDO path not found: {self.mondo_path}")
            return

        if pronto is None:
            raise ImportError("pronto is required: pip install pronto")

        ont = pronto.Ontology(self.mondo_path)

        for term in ont.terms():
            if not term.id.startswith("MONDO:"):
                continue
            name = term.name or ""
            if not name:
                continue

            self._add_label(term.id, name)

            for syn in term.synonyms:
                if syn.description:
                    self._add_label(term.id, syn.description)

            # Collect Orphanet xrefs
            for xref in term.xrefs:
                xid = xref.id if hasattr(xref, "id") else str(xref)
                if xid.startswith("Orphanet:") or xid.startswith("ORPHA:"):
                    orphanet_id = re.sub(r"^Orphanet:", "ORPHA:", xid)
                    self.mondo_to_orphanet[term.id].append(orphanet_id)

    # ------------------------------------------------------------------
    # OMIM / Orphanet from phenotype.hpoa
    # ------------------------------------------------------------------

    def _load_hpoa(self) -> None:
        if not os.path.exists(self.hpoa_path):
            if self.debug:
                print(f"[OntologyIndex] phenotype.hpoa not found: {self.hpoa_path}")
            return

        seen: Set[str] = set()
        with open(self.hpoa_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                disease_id = parts[0].strip()
                disease_name = parts[1].strip()

                if not disease_id.startswith(_DISEASE_PREFIXES):
                    continue
                if not disease_name:
                    continue
                if disease_id in seen:
                    continue
                seen.add(disease_id)

                self._add_label(disease_id, disease_name)

    # ------------------------------------------------------------------
    # Helper: add one (term_id, label) pair to exact_map + primary_label
    # ------------------------------------------------------------------

    def _add_label(self, term_id: str, label: str) -> None:
        key = label.lower().strip()
        if not key:
            return
        if term_id not in self.exact_map.get(key, []):
            self.exact_map[key].append(term_id)
        # Register canonical name (first one wins)
        if term_id not in self.primary_label:
            self.primary_label[term_id] = label

    # ------------------------------------------------------------------
    # Stemmed map
    # ------------------------------------------------------------------

    def _build_stemmed_map(self) -> None:
        for label_lower, ids in self.exact_map.items():
            stems = get_stemmed_tokens(label_lower)
            if stems:
                for tid in ids:
                    if tid not in self.stemmed_map[stems]:
                        self.stemmed_map[stems].append(tid)

    # ------------------------------------------------------------------
    # Aho-Corasick automaton
    # ------------------------------------------------------------------

    def _build_ac_automaton(self) -> None:
        if ahocorasick is None:
            if self.debug:
                print("[OntologyIndex] pyahocorasick not available; AC disabled.")
            return

        A = ahocorasick.Automaton()
        for label_lower in self.exact_map:
            A.add_word(label_lower, label_lower)
        A.make_automaton()
        self.ac_automaton = A

    # ------------------------------------------------------------------
    # Phenotype word set (for NER n-gram pre-filtering)
    # ------------------------------------------------------------------

    _PHENOTYPE_WORD_STOPWORDS = frozenset({
        "of", "the", "and", "with", "abnormality", "abnormalities",
        "increased", "decreased", "abnormal", "type", "that", "this",
        "for", "from", "are", "was", "were", "been", "has", "have",
    })

    def _build_phenotype_word_set(self) -> None:
        """
        Extract single words (>=4 chars, excluding stopwords) from all
        HPO-related exact_map keys. Used as a cheap pre-filter: only n-grams
        containing at least one word in this set proceed to embedding.
        """
        words: Set[str] = set()
        for label_lower, term_ids in self.exact_map.items():
            # Only include labels that map to at least one HPO phenotype term
            if not any(tid.startswith("HP:") for tid in term_ids):
                continue
            for word in label_lower.split():
                if len(word) >= 4 and word not in self._PHENOTYPE_WORD_STOPWORDS:
                    words.add(word)
        self.phenotype_word_set = words
        if self.debug:
            print(f"[OntologyIndex] phenotype_word_set: {len(words)} words")

    # ------------------------------------------------------------------
    # Embeddings
    # ------------------------------------------------------------------

    def _build_embeddings(self) -> None:
        cache_file = os.path.join(
            self.cache_dir,
            f"embeddings_{self.embedding_model_name.replace('/', '_')}.pkl",
        )

        if os.path.exists(cache_file):
            if self.debug:
                print(f"[OntologyIndex] Loading embeddings from cache: {cache_file}")
            with open(cache_file, "rb") as f:
                data = pickle.load(f)
            self.embeddings = data["embeddings"]
            self.embedding_labels = data["labels"]
            return

        try:
            from sentence_transformers import SentenceTransformer
        except ImportError:
            if self.debug:
                print("[OntologyIndex] sentence-transformers not available; ANN disabled.")
            return

        labels = list(self.exact_map.keys())
        if not labels:
            return

        if self.debug:
            print(
                f"[OntologyIndex] Encoding {len(labels)} labels with "
                f"{self.embedding_model_name} …"
            )

        model = SentenceTransformer(self.embedding_model_name, device=self.device)
        vecs = model.encode(
            labels,
            batch_size=256,
            show_progress_bar=self.debug,
            convert_to_numpy=True,
            normalize_embeddings=True,
        )
        self.embeddings = vecs
        self.embedding_labels = labels

        with open(cache_file, "wb") as f:
            pickle.dump({"embeddings": vecs, "labels": labels}, f)

        if self.debug:
            print(f"[OntologyIndex] Embeddings cached to {cache_file}")
