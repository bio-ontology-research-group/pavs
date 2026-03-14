"""
Main PhenotypeMatcher class implementing Graph RAG for phenotype normalization.
"""

import os
import json
import pickle
import logging
import time
import re
import tempfile
import shutil
import traceback
from typing import List, Dict, Any, Optional, Set
import numpy as np
from pyhpo import Ontology
import pronto
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity
import requests
from tqdm import tqdm

from .schemas import (
    PhenotypeInput,
    PhenotypeOutput,
    PhenotypeMatch,
    DiseaseMatch,
    MatcherConfig,
)
from . import acronyms
from . import ner
from . import normalization_utils


# Embedding model presets
EMBEDDING_MODELS = {
    "fast": "all-MiniLM-L6-v2",
    "balanced": "all-mpnet-base-v2",
    "accurate": "BAAI/bge-large-en-v1.5",
    "sapbert": "cambridgeltl/SapBERT-from-PubMedBERT-fulltext",
}

# LLM model presets
LLM_MODELS = {
    "fast": "openai/gpt-oss-120b",
    "accurate": "anthropic/claude-3.5-sonnet",
    "balanced": "google/gemini-2.0-flash-exp:free",
}


class PhenotypeMatcher:
    def __init__(self, config: Optional[MatcherConfig] = None):
        self.config = config or MatcherConfig()
        self.embedding_model_name = EMBEDDING_MODELS.get(self.config.embedding_model, self.config.embedding_model)
        self.llm_model_name = LLM_MODELS.get(self.config.llm_model, self.config.llm_model)
        log_level = logging.DEBUG if self.config.debug else logging.INFO
        logging.basicConfig(level=log_level, format="%(asctime)s - %(levelname)s - %(message)s")
        self.logger = logging.getLogger(__name__)
        os.makedirs(self.config.cache_dir, exist_ok=True)
        self.sentence_model, self.embeddings, self.terms, self.term_id_map = None, None, [], {}
        self.phenotype_descendants, self.severity_descendants = set(), set()
        self.api_key = self.config.api_key or os.getenv("OPENROUTER_API_KEY")
        self._initialize()

    def _initialize(self):
        self.logger.info("Initializing PhenotypeMatcher...")
        self._load_ontologies()
        self._build_descendant_sets()
        self._load_or_compute_embeddings()
        self.logger.info("Initialization complete.")

    def _load_ontologies(self):
        self.logger.info("Loading ontologies...")
        stemmed_cache_path = os.path.join(self.config.cache_dir, "stemmed_cache.pkl")
        normalization_utils.load_persistent_cache(stemmed_cache_path)
        try:
            hpo_dir = self.config.hpo_data_dir or (os.path.dirname(self.config.hpo_path) if os.path.isfile(self.config.hpo_path) else self.config.hpo_path)
            req = ["hp.obo", "phenotype.hpoa", "genes_to_phenotype.txt"]
            missing = [f for f in req if not os.path.exists(os.path.join(hpo_dir, f))]
            if missing and "hp.obo" not in missing:
                temp_hpo_dir = tempfile.mkdtemp(prefix="hpo_data_")
                os.symlink(os.path.abspath(self.config.hpo_path), os.path.join(temp_hpo_dir, "hp.obo"))
                for f in ["phenotype.hpoa", "genes_to_phenotype.txt"]:
                    found = False
                    for d in [hpo_dir, "data", "../data", "../../data", "/mnt/data1/pavs/data"]:
                        src = os.path.join(d, f)
                        if os.path.exists(src):
                            os.symlink(os.path.abspath(src), os.path.join(temp_hpo_dir, f))
                            found = True; break
                hpo_dir = temp_hpo_dir
            _ = Ontology(hpo_dir)
            for term in Ontology:
                syns = list(term.synonym) if hasattr(term, "synonym") else []
                self.terms.append({"id": str(term.id), "name": str(term.name), "synonyms": [str(s) for s in syns], "source": "HPO"})
            self.logger.info(f"Loaded {len(self.terms)} HPO terms")
        except Exception as e: self.logger.error(f"Failed to load HPO: {e}")

        mondo_path = self.config.mondo_path
        if os.path.exists(mondo_path):
            try:
                ms = pronto.Ontology(mondo_path)
                for term in ms.terms():
                    if term.id.startswith("MONDO:"):
                        self.terms.append({"id": str(term.id), "name": str(term.name), "synonyms": [str(s.description) for s in term.synonyms], "source": "MONDO"})
            except Exception as e: self.logger.error(f"Failed to load MONDO: {e}")

        self.term_id_map = {t["id"]: i for i, t in enumerate(self.terms)}
        self.exact_map, self.stemmed_map = {}, {}
        for i, t in enumerate(tqdm(self.terms, desc="Building lookup maps")):
            for name in [t["name"]] + t["synonyms"]:
                nl = str(name).lower().strip()
                if not nl: continue
                if nl not in self.exact_map: self.exact_map[nl] = []
                if i not in self.exact_map[nl]: self.exact_map[nl].append(i)
                st = normalization_utils.get_stemmed_tokens(nl)
                if st:
                    if st not in self.stemmed_map: self.stemmed_map[st] = []
                    if i not in self.stemmed_map[st]: self.stemmed_map[st].append(i)
        normalization_utils.save_persistent_cache(stemmed_cache_path)

    def _build_descendant_sets(self):
        try:
            pheno_root = Ontology.get_hpo_object("HP:0000118")
            if pheno_root:
                self.phenotype_descendants = {str(pheno_root.id)}
                def add_d(t):
                    for c in t.children:
                        if str(c.id) not in self.phenotype_descendants: self.phenotype_descendants.add(str(c.id)); add_d(c)
                add_d(pheno_root)
            sev_root = Ontology.get_hpo_object("HP:0012824")
            if sev_root:
                self.severity_descendants = {str(sev_root.id)}
                def add_s(t):
                    for c in t.children:
                        if str(c.id) not in self.severity_descendants: self.severity_descendants.add(str(c.id)); add_s(c)
                add_s(sev_root)
        except: pass

    def _load_or_compute_embeddings(self):
        cache_file = os.path.join(self.config.cache_dir, f"embeddings_{self.embedding_model_name.replace('/', '_')}.pkl")
        if os.path.exists(cache_file):
            try:
                with open(cache_file, "rb") as f:
                    data = pickle.load(f); self.embeddings = data["embeddings"]
                    if len(self.embeddings) != len(self.terms): self._compute_embeddings(cache_file)
            except: self._compute_embeddings(cache_file)
        else: self._compute_embeddings(cache_file)

    def _compute_embeddings(self, cache_path: str):
        self.logger.info(f"Computing embeddings...")
        self.sentence_model = SentenceTransformer(self.embedding_model_name, device=self.config.device)
        texts = [f"Name: {t['name']}. Synonyms: {', '.join(t['synonyms'])}." for t in self.terms]
        self.embeddings = self.sentence_model.encode(texts, show_progress_bar=True, convert_to_numpy=True, normalize_embeddings=True)
        with open(cache_path, "wb") as f: pickle.dump({"embeddings": self.embeddings}, f)

    def _add_phenotype(self, phenos: List[PhenotypeMatch], t_obj: Dict, excluded: bool = False, conf: float = 1.0, sev_id: str = None, sev_lab: str = None, method: str = "literal"):
        if not any(p.hpo_id == t_obj["id"] for p in phenos):
            phenos.append(PhenotypeMatch(hpo_id=t_obj["id"], label=t_obj["name"], excluded=excluded, confidence=conf, severity_id=sev_id, severity_label=sev_lab, matched_by=method))

    def _add_disease(self, diseases: List[DiseaseMatch], t_obj: Dict, confidence: float = 1.0, method: str = "literal"):
        if not any(d.mondo_id == t_obj["id"] for d in diseases):
            diseases.append(DiseaseMatch(mondo_id=t_obj["id"], mondo_label=t_obj["name"], confidence=confidence, matched_by=method))

    def _is_negated(self, text: str, term_name: str) -> bool:
        tl, tn = text.lower(), term_name.lower()
        if tn not in tl: return False
        pre = tl.split(tn)[0].strip().split()
        if not pre: return False
        return any(w.strip(",.;") in {"no", "not", "without", "normal", "negative", "denies"} for w in pre[-3:])

    def _find_literal_matches(self, text: str) -> List[Dict]:
        if not isinstance(text, str): return []
        cl = re.sub(r"[^\w\s]", " ", text.lower()); words = cl.split(); matches, seen = [], set()
        for length in range(min(len(words), 10), 0, -1):
            for i in range(len(words) - length + 1):
                ngram = " ".join(words[i : i + length])
                if not ngram or len(ngram) < 3: continue
                # REDUCED FALSE POSITIVES
                if len(ngram) < 5 and ngram.upper() not in acronyms.MEDICAL_ACRONYMS:
                    if ngram.lower() in {"caught", "noisy", "short", "absent", "long", "small", "large", "mild", "severe", "present"}: continue
                for t_idx in self.exact_map.get(ngram, []):
                    t = self.terms[t_idx]
                    if t["id"] not in seen: matches.append(t); seen.add(t["id"])
                st = normalization_utils.get_stemmed_tokens(ngram)
                if st and st in self.stemmed_map:
                    for t_idx in self.stemmed_map[st]:
                        t = self.terms[t_idx]
                        if t["id"] not in seen: matches.append(t); seen.add(t["id"])
        return matches

    def match(self, input_data: PhenotypeInput) -> PhenotypeOutput:
        start_time = time.time(); llm_calls = 0
        self.logger.info("Layer 1: Relaxed NER...")
        phrases = ner.relaxed_ner(input_data.text)
        all_phenos, all_diseases, todo_llm = [], [], []
        try:
            for p_data in phrases:
                term = p_data["term"]; tl = str(term).lower().strip(); orig_text = p_data["context"]
                if not tl: continue
                
                # Perfect Match
                perfect_indices = self.exact_map.get(tl, [])
                if not perfect_indices:
                    st = normalization_utils.get_stemmed_tokens(tl)
                    if st: perfect_indices = self.stemmed_map.get(st, [])
                
                h_perf = [self.terms[i] for i in perfect_indices if self.terms[i]["id"].startswith("HP:")]
                m_perf = [self.terms[i] for i in perfect_indices if self.terms[i]["id"].startswith("MONDO:")]
                
                h_done, m_done = False, False
                if len(h_perf) == 1 and tl.upper() not in acronyms.MEDICAL_ACRONYMS:
                    self._add_phenotype(all_phenos, h_perf[0], excluded=self._is_negated(tl, h_perf[0]["name"]), method="perfect_match")
                    h_done = True
                if len(m_perf) == 1:
                    self._add_disease(all_diseases, m_perf[0], method="perfect_match")
                    m_done = True
                if h_done and m_done: continue

                # Literal Inclusion
                literal = self._find_literal_matches(tl)
                h_lit = [t for t in literal if t["id"].startswith("HP:") and t["id"] in self.phenotype_descendants]
                m_lit = [t for t in literal if t["id"].startswith("MONDO:")]
                
                # Filter redundant literal matches within the phrase
                h_final = []
                for t1 in h_lit:
                    if not any(t2["id"] != t1["id"] and t1["name"].lower() in t2["name"].lower() for t2 in h_lit):
                        h_final.append(t1)
                
                added = False
                for t in h_final:
                    if t["name"].upper() not in acronyms.MEDICAL_ACRONYMS:
                        self._add_phenotype(all_phenos, t, excluded=self._is_negated(tl, t["name"]), method="literal_inclusion")
                        added = True
                for t in m_lit:
                    self._add_disease(all_diseases, t, method="literal_inclusion"); added = True
                
                if any(t["name"].upper() in acronyms.MEDICAL_ACRONYMS for t in h_lit) or (not added and not literal):
                    todo_llm.append({"phrase_data": p_data, "candidates": h_lit + m_lit, "type": "disambiguation" if any(t["name"].upper() in acronyms.MEDICAL_ACRONYMS for t in h_lit) else "rag"})
            
            if todo_llm:
                pbar = tqdm(todo_llm, desc="LLM Processing", leave=False, disable=len(todo_llm) <= 1)
                for item in pbar:
                    t_str = str(item["phrase_data"]["term"])
                    if item["type"] == "rag" and not item["candidates"]:
                        item["candidates"] = [c["term"] for c in self._retrieve_phenotype_candidates(t_str, k=self.config.top_k_phenotype)]
                    res = self._call_llm_layered(t_str, input_data.text, item["candidates"], item["type"], input_data.gene_hint)
                    llm_calls += 1
                    if res and isinstance(res, dict):
                        severity_label = res.get("severity_label")
                        severity_id = self._lookup_severity_id(severity_label) if severity_label else None
                        for pl in res.get("phenotype_labels", []):
                            hid = self._lookup_hpo_id(pl)
                            if hid:
                                csid, cslab = severity_id, severity_label
                                if not csid:
                                    for sw, sh in [("mild", "HP:0012825"), ("moderate", "HP:0012826"), ("severe", "HP:0012828"), ("profound", "HP:0012829")]:
                                        if sw in pl.lower() or sw in t_str.lower(): csid, cslab = sh, sw.capitalize(); break
                                self._add_phenotype(all_phenos, {"id": hid, "name": pl}, res.get("excluded", False), 0.8, csid, cslab, method="LLM_"+item["type"])
                        for d in res.get("disease_labels", []):
                            if isinstance(d, dict) and (ml := d.get("mondo")):
                                if (mid := self._lookup_mondo_id(ml)): self._add_disease(all_diseases, {"id": mid, "name": ml}, 0.8, method="LLM_"+item["type"])
            
            # FINAL CROSS-PHRASE FILTER: remove broader or incorrect terms
            filtered_phenos = []
            for p1 in all_phenos:
                # If we have a more specific version of the same term (e.g. "Severe intellectual disability" vs "Intellectual disability"), keep only specific
                # OR if we have mutually exclusive terms due to context errors
                is_redundant = False
                for p2 in all_phenos:
                    if p1.hpo_id != p2.hpo_id:
                        # Simple name-based redundancy check
                        if p1.label.lower() in p2.label.lower():
                            is_redundant = True; break
                if not is_redundant: filtered_phenos.append(p1)
            all_phenos = filtered_phenos

        except Exception as e: self.logger.error(f"Matcher failed: {e}\n{traceback.format_exc()}")
        return PhenotypeOutput(phenotypes=all_phenos, diseases=all_diseases, raw_input=input_data.text, ner_extracted_terms=phrases,
                               processing_metadata={"terms_processed": len(phrases), "llm_calls": llm_calls, "processing_time_seconds": time.time() - start_time, "embedding_model": self.embedding_model_name, "llm_model": self.llm_model_name})

    def _call_llm_layered(self, term: str, text: str, cands: List[Dict], mtype: str, gene: str = None) -> Dict[str, Any]:
        cb = "".join([f"- {c['name']} ({c['id']})\n" + (f"  Definition: {c['definition']}\n" if c.get('definition') else "") for c in cands])
        if mtype == "rag":
            prompt = f"Identify terms equivalent to '{term}' in '{text}'. Rules: match ALL words, ignore 'PICU' etc. JSON: {{phenotype_labels: [], severity_label: str|null, excluded: bool}}\nCandidates:\n{cb}"
        else:
            prompt = f"Disambiguate '{term}' in '{text}' (Gene: {gene}). ASD heart -> Atrial septal defect; ASD behavior -> Autistic behavior. JSON: {{phenotype_labels: [], severity_label: str|null, excluded: bool}}\nCandidates:\n{cb}"
        return self._query_llm(prompt)

    def _retrieve_phenotype_candidates(self, query: str, k: int = 5) -> List[Dict]:
        if self.sentence_model is None: self.sentence_model = SentenceTransformer(self.embedding_model_name, device=self.config.device)
        idx = [i for i, t in enumerate(self.terms) if t["source"] == "HPO" and t["id"] in self.phenotype_descendants]
        if not idx: return []
        scores = cosine_similarity(self.sentence_model.encode([str(query)], convert_to_numpy=True, normalize_embeddings=True), self.embeddings[idx])[0]
        return [{"term": self.terms[idx[i]], "score": float(scores[i])} for i in np.argsort(scores)[::-1][:k]]

    def _query_llm(self, prompt: str) -> Dict[str, Any]:
        if not self.api_key: return {}
        try:
            resp = requests.post("https://openrouter.ai/api/v1/chat/completions", headers={"Authorization": f"Bearer {self.api_key}", "Content-Type": "application/json"},
                                json={"model": self.llm_model_name, "messages": [{"role": "user", "content": prompt}], "response_format": {"type": "json_object"}}, timeout=60)
            c = resp.json()["choices"][0]["message"]["content"].strip()
            if c.startswith("```"): c = re.sub(r"^```(?:json)?\n", "", c); c = re.sub(r"\n```$", "", c)
            s, e = c.find("{"), c.rfind("}")
            return json.loads(c[s:e+1]) if s != -1 else {}
        except: return {}

    def _lookup_hpo_id(self, label: Any) -> Optional[str]:
        if isinstance(label, dict): label = label.get("label") or label.get("name")
        if not isinstance(label, str): return None
        l = label.lower().strip()
        for t in self.terms:
            if t["source"] == "HPO" and (t["name"].lower() == l or any(s.lower() == l for s in t.get("synonyms", []))): return t["id"]
        return None

    def _lookup_severity_id(self, label: Any) -> Optional[str]:
        if isinstance(label, dict): label = label.get("label") or label.get("name")
        if not isinstance(label, str): return None
        l = label.lower().strip()
        for t in self.terms:
            if t["source"] == "HPO" and t["id"] in self.severity_descendants and t["name"].lower() == l: return t["id"]
        return None

    def _lookup_mondo_id(self, label: Any) -> Optional[str]:
        if isinstance(label, dict): label = label.get("label") or label.get("name")
        if not isinstance(label, str): return None
        l = label.lower().strip()
        for t in self.terms:
            if t["source"] == "MONDO" and t["name"].lower() == l: return t["id"]
        return None
