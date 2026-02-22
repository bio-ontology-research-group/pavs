"""
PhenotypeMatcher v2 — Algorithms 5 & 6.

Algorithm 5: process_matches(l_found, ctx, trusted) → Set[Tuple]
Algorithm 6: match(s: str) → PhenotypeOutput
"""

import os
import threading
from collections import OrderedDict
from dataclasses import replace as dc_replace
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

from .detection import (
    _NEG_ENCODED_LABEL_STARTS,
    det_mod,
    det_mod_in_text,
    det_neg,
    det_neg_in_text,
    is_boundary_valid,
)
from .llm import llm_check_modifier, llm_check_negation, llm_dis, llm_val_batch
from .normalization_utils import get_stemmed_tokens
from .ontology_index import OntologyIndex
from .schemas import (
    DiseaseMatch,
    MatcherConfig,
    PhenotypeMatch,
    PhenotypeOutput,
)
from .tokenizer import tokenize_expand

# ---------------------------------------------------------------------------
# Model preset maps
# ---------------------------------------------------------------------------

EMBEDDING_MODELS: Dict[str, str] = {
    "sapbert":  "cambridgeltl/SapBERT-from-PubMedBERT-fulltext",
    "fast":     "all-MiniLM-L6-v2",
    "balanced": "all-mpnet-base-v2",
    "accurate": "BAAI/bge-large-en-v1.5",
    "medical":  "pritamdeka/S-PubMedBert-MS-MARCO",
}

LLM_MODELS: Dict[str, str] = {
    "deepseek": "deepseek/deepseek-chat",
    "gpt4oss":  "openai/gpt-oss-120b",
    "accurate": "anthropic/claude-sonnet-3-5",
    "gemini":   "google/gemini-2.0-flash-exp:free",
}

# ---------------------------------------------------------------------------
# Stop-word set (skip fuzzy matching for these)
# ---------------------------------------------------------------------------

W_STOP: Set[str] = {
    "and", "or",
    "age", "sex", "pain", "type", "the", "a", "an", "is", "are", "was",
    "were", "be", "been", "have", "has", "had", "do", "does", "did",
    "to", "of", "in", "on", "at", "by", "for", "with", "about", "as",
    "into", "from", "out", "case", "report", "patient",
    "diagnosis", "clinical", "finding", "feature", "sign", "symptom",
    "present", "status",
}

# Fix 1b: Severity words that, when present in a term label, suppress the
# rule-based severity modifier (to avoid double-encoding e.g. "Severe
# intellectual disability" + modifier "severe").
_SEVERITY_WORDS = frozenset([
    "severe", "mild", "moderate", "profound", "borderline", "marked",
])

# Tuple fields: (term_id: str, modifier_id: Optional[str], negated: bool)
_MatchTuple = Tuple[str, Optional[str], bool]


class PhenotypeMatcher:
    """
    Maps free-text clinical descriptions to HPO / MONDO / OMIM / Orphanet IDs.

    Usage::

        cfg = MatcherConfig(hpo_path="ontology/hp.obo")
        m = PhenotypeMatcher(cfg)
        out = m.match("severe intellectual disability, no seizures")
        print(out.get_hpo_ids())
    """

    def __init__(self, config: Optional[MatcherConfig] = None):
        self.cfg = config or MatcherConfig()

        # Resolve model presets
        emb_name = EMBEDDING_MODELS.get(self.cfg.embedding_model, self.cfg.embedding_model)
        self._llm_model = LLM_MODELS.get(self.cfg.llm_model, self.cfg.llm_model)
        self._api_key = self.cfg.api_key or os.environ.get("OPENROUTER_API_KEY")

        # Load ontology index
        self._idx = OntologyIndex(
            hpo_path=self.cfg.hpo_path,
            mondo_path=self.cfg.mondo_path,
            hpoa_path=self.cfg.hpoa_path,
            embedding_model_name=emb_name,
            cache_dir=self.cfg.cache_dir,
            device=self.cfg.device,
            debug=self.cfg.debug,
        )

        # Lazy-load sentence transformer (only if needed and not already loaded
        # inside OntologyIndex)
        self._encoder = None
        self._encoder_lock = threading.Lock()

        # Shared LLM response cache: keyed by (function_tag, frozenset(ids), text, model)
        # Thread-safe under the GIL for individual get/set operations.
        self._llm_cache: dict = {}

        # Match-level cache: keyed by stripped input string → PhenotypeOutput
        # Uses OrderedDict for O(1) LRU eviction.
        self._match_cache: OrderedDict = OrderedDict()

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def match(self, s: str) -> PhenotypeOutput:
        """
        Algorithm 6 — match a raw clinical string.

        Returns a PhenotypeOutput with deduplicated HPO / disease hits.
        """
        # Match-level cache: return immediately for repeated identical inputs
        if self.cfg.match_cache_size > 0:
            s_key = s.strip()
            if s_key in self._match_cache:
                return self._match_cache[s_key]

        O: Set[_MatchTuple] = set()
        idx = self._idx

        # ----------------------------------------------------------------
        # Phase 1 — whole string
        # ----------------------------------------------------------------
        s_lower = s.lower().strip()

        # 1a — exact whole-string match (trusted)
        if s_lower in idx.exact_map:
            O |= self._process_matches([s_lower], s, trusted=True)

        # 1b — stemmed whole-string match (untrusted)
        # Skip entirely if the string contains a negation trigger — the
        # stemmer removes stop words like "no", so a stemmed lookup of
        # "no seizures" would incorrectly match HP:0001250 as present.
        # Phase 2d handles negation correctly via the AC automaton.
        if not det_neg_in_text(s_lower):
            stems = get_stemmed_tokens(s)
            if stems in idx.stemmed_map:
                exact_ids = set(idx.exact_map.get(s_lower, []))
                stem_ids = [tid for tid in idx.stemmed_map[stems]
                            if tid not in exact_ids]
                if stem_ids:
                    # Use ONE representative label per term ID (the primary label).
                    # This avoids calling llm_val_batch dozens of times for all synonyms.
                    seen_labels: set = set()
                    stem_labels_primary: List[str] = []
                    for tid in stem_ids:
                        plab = idx.primary_label.get(tid, "").lower()
                        if plab and plab in idx.exact_map and plab not in seen_labels:
                            seen_labels.add(plab)
                            stem_labels_primary.append(plab)
                    if stem_labels_primary:
                        O |= self._process_matches(stem_labels_primary, s, trusted=False)

        # 1c — fuzzy whole-string match (untrusted), skip if in W_STOP
        if s_lower not in W_STOP:
            fuzzy_labels = self._fuzzy_lookup(s_lower)
            # Exclude labels already covered
            covered = set(idx.exact_map.get(s_lower, []))
            fuzzy_labels = [lbl for lbl in fuzzy_labels if lbl not in [
                l2 for l2 in idx.exact_map if any(t in covered for t in idx.exact_map[l2])
            ]]
            if fuzzy_labels:
                O |= self._process_matches(fuzzy_labels, s, trusted=False)

        # ----------------------------------------------------------------
        # Phase 2 — token loop
        # ----------------------------------------------------------------
        for t in tokenize_expand(s):
            t_lower = t.lower().strip()
            if not t_lower:
                continue

            # 2a — exact token match (trusted)
            O_new: Set[_MatchTuple] = set()
            if t_lower in idx.exact_map:
                O_new |= self._process_matches([t_lower], t, trusted=True)

            # 2b — stemmed token match (untrusted)
            # Skip if the token contains negation (e.g. "no seizures"):
            # the stemmer strips "no", so we would incorrectly match as present.
            # Phase 2d handles such tokens via the AC automaton with det_neg.
            if not det_neg_in_text(t_lower):
                t_stems = get_stemmed_tokens(t)
                if t_stems in idx.stemmed_map:
                    t_exact_ids = set(idx.exact_map.get(t_lower, []))
                    t_stem_ids = [tid for tid in idx.stemmed_map[t_stems]
                                  if tid not in t_exact_ids]
                    if t_stem_ids:
                        seen_t: set = set()
                        t_stem_labels_primary: List[str] = []
                        for tid in t_stem_ids:
                            plab = idx.primary_label.get(tid, "").lower()
                            if plab and plab in idx.exact_map and plab not in seen_t:
                                seen_t.add(plab)
                                t_stem_labels_primary.append(plab)
                        if t_stem_labels_primary:
                            O_new |= self._process_matches(
                                t_stem_labels_primary, t, trusted=False
                            )

            # Merge 2a/2b results into O regardless, but don't skip Phase 2d.
            # Rationale: 2b may find a DISEASE term (e.g. MONDO) for a token
            # like "I have a seizure disorder" but the AC automaton in 2d is still
            # needed to find the embedded HPO PHENOTYPE term "seizure".
            # Phase 2c (fuzzy) is skipped when 2a/2b already found something
            # because it adds lower-quality edit-distance matches.
            if O_new:
                O |= O_new

            # 2c — fuzzy token match (untrusted), skip W_STOP and skip when
            # 2a/2b already yielded results.
            if not O_new and t_lower not in W_STOP:
                fuzzy_tok = self._fuzzy_lookup(t_lower)
                if fuzzy_tok:
                    O_fuz = self._process_matches(fuzzy_tok, t, trusted=False)
                    if O_fuz:
                        O |= O_fuz
                        continue

            # 2d — Aho-Corasick substring search
            found_phrase = False
            if idx.ac_automaton is not None:
                for end_idx, label in idx.ac_automaton.iter(t_lower):
                    start_idx = end_idx - len(label) + 1
                    # end_idx here is the last char index (inclusive)
                    # is_boundary_valid uses exclusive end
                    if not is_boundary_valid(t_lower, start_idx, end_idx + 1):
                        continue

                    mod_id = det_mod(
                        t_lower, start_idx, end_idx + 1,
                        idx.modifier_map, win=self.cfg.det_mod_win,
                    )
                    # Fix 3a: pass matched_label so triggers embedded in the
                    # term name (e.g. "absent" in "Absent speech") are skipped.
                    negated = det_neg(
                        t_lower, start_idx, end_idx + 1,
                        win=self.cfg.det_neg_win,
                        matched_label=label,
                    )

                    p_cand = idx.exact_map[label]
                    if len(p_cand) > 1:
                        chosen = llm_dis(
                            p_cand, t, idx.primary_label,
                            api_key=self._api_key, model=self._llm_model,
                            cache=self._llm_cache,
                        )
                        p_res = chosen if chosen else set(p_cand)
                    else:
                        p_res = set(p_cand)

                    for pid in p_res:
                        O.add((pid, mod_id, negated))
                    found_phrase = True

            if found_phrase:
                continue

            # 2e — ANN / semantic match
            ann_labels, ann_scores = self._ann_search(t, k=self.cfg.ann_top_k)
            if ann_labels:
                # Strategy 1: floor + relative margin + cap
                filtered = self._strategy1_filter(
                    ann_labels, ann_scores,
                    theta_abs=self.cfg.theta_abs,
                    delta=self.cfg.delta,
                    cap=self.cfg.ann_top_k,
                )
                if filtered:
                    mod_id = det_mod_in_text(t_lower, idx.modifier_map)
                    negated = det_neg_in_text(t_lower)
                    for lbl in filtered:
                        # If the ANN-matched label already encodes absence
                        # (e.g. "Absent speech"), the negation trigger in the
                        # token (e.g. "no" in "no speech") is what caused the
                        # semantic match — not an external negation of the term.
                        # Override negated=False for such labels.
                        lbl_negated = negated
                        if negated and any(
                            lbl.lower().startswith(pfx)
                            for pfx in _NEG_ENCODED_LABEL_STARTS
                        ):
                            lbl_negated = False
                        for pid in idx.exact_map.get(lbl, []):
                            O.add((pid, mod_id, lbl_negated))

        output = self._build_output(O, s)

        # Fix 6c: LLM-based post-hoc validation of negation and modifier
        # decisions (only when api_key is configured).
        output.phenotypes = self._llm_validate(output.phenotypes, s)

        # Store in match-level cache with bounded eviction
        if self.cfg.match_cache_size > 0:
            s_key = s.strip()
            if s_key not in self._match_cache:
                if len(self._match_cache) >= self.cfg.match_cache_size:
                    self._match_cache.popitem(last=False)
                self._match_cache[s_key] = output

        return output

    # ------------------------------------------------------------------
    # Algorithm 5: process_matches
    # ------------------------------------------------------------------

    def _process_matches(
        self,
        labels: List[str],
        ctx: str,
        trusted: bool,
    ) -> Set[_MatchTuple]:
        """
        Algorithm 5 — resolve labels to (term_id, mod_id, negated) tuples.

        *trusted* paths skip LLM-Val; untrusted paths call llm_val_batch.
        """
        idx = self._idx
        result: Set[_MatchTuple] = set()

        for label in labels:
            p_cand = list(idx.exact_map.get(label, []))
            if not p_cand:
                continue

            if trusted:
                # Trusted path (exact match): add ALL IDs for this label without
                # LLM disambiguation.  This preserves every namespace (MONDO,
                # OMIM, ORPHA) for multi-ID labels such as disease names.
                for pid in p_cand:
                    result.add((pid, None, False))
            else:
                # Untrusted path (stemmed / fuzzy): disambiguate if ambiguous,
                # then validate with LLM.
                if len(p_cand) > 1:
                    p_res = llm_dis(
                        p_cand, ctx, idx.primary_label,
                        api_key=self._api_key, model=self._llm_model,
                        cache=self._llm_cache,
                    )
                    if not p_res:
                        p_res = set(p_cand)
                else:
                    p_res = set(p_cand)

                valid = llm_val_batch(
                    list(p_res), ctx, idx.primary_label,
                    api_key=self._api_key, model=self._llm_model,
                    cache=self._llm_cache,
                )
                for pid in valid:
                    result.add((pid, None, False))

        return result

    # ------------------------------------------------------------------
    # ANN search
    # ------------------------------------------------------------------

    def _ann_search(self, text: str, k: int) -> Tuple[List[str], List[float]]:
        """
        Encode *text* with the sentence transformer and return top-k
        (label, cosine_score) pairs from the pre-built embedding index.
        """
        idx = self._idx
        if idx.embeddings is None or not idx.embedding_labels:
            return [], []

        vec = self._encode(text)
        if vec is None:
            return [], []

        # Cosine similarity (embeddings are L2-normalised)
        sims = idx.embeddings @ vec  # shape (N,)
        top_indices = np.argsort(sims)[::-1][:k]

        labels = [idx.embedding_labels[i] for i in top_indices]
        scores = [float(sims[i]) for i in top_indices]
        return labels, scores

    def _encode(self, text: str):
        """Lazy-load encoder and return a normalised embedding vector.

        Uses double-checked locking so that concurrent threads do not load
        the SentenceTransformer model simultaneously when --workers > 1.
        """
        if self._encoder is None:
            with self._encoder_lock:
                if self._encoder is None:
                    try:
                        from sentence_transformers import SentenceTransformer
                        emb_name = EMBEDDING_MODELS.get(
                            self.cfg.embedding_model, self.cfg.embedding_model
                        )
                        self._encoder = SentenceTransformer(emb_name, device=self.cfg.device)
                    except Exception:
                        return None

        vec = self._encoder.encode(
            [text],
            convert_to_numpy=True,
            normalize_embeddings=True,
        )[0]
        return vec

    # ------------------------------------------------------------------
    # Strategy 1 filter
    # ------------------------------------------------------------------

    @staticmethod
    def _strategy1_filter(
        labels: List[str],
        scores: List[float],
        theta_abs: float = 0.70,
        delta: float = 0.05,
        cap: int = 5,
    ) -> List[str]:
        """
        Hybrid Dynamic Thresholding (Strategy 1):
        1. Drop candidates below absolute floor θ_abs.
        2. Keep only those within [S_max − δ, S_max].
        3. Cap at *cap* results.
        """
        if not scores:
            return []

        s_max = max(scores)
        floor = max(theta_abs, s_max - delta)

        filtered = [
            lbl for lbl, sc in zip(labels, scores)
            if sc >= floor
        ]
        return filtered[:cap]

    # ------------------------------------------------------------------
    # Fuzzy lookup (edit distance 1)
    # ------------------------------------------------------------------

    def _fuzzy_lookup(self, s: str) -> List[str]:
        """
        Return labels in exact_map at Levenshtein distance exactly 1 from *s*.
        """
        try:
            from rapidfuzz.process import extract
            from rapidfuzz.distance import Levenshtein

            hits = extract(
                s,
                list(self._idx.exact_map.keys()),
                scorer=Levenshtein.distance,
                score_cutoff=1,
                limit=20,
            )
            # extract returns (match, score, idx); keep score == 1
            return [h[0] for h in hits if h[1] == 1]
        except Exception:
            return []

    # ------------------------------------------------------------------
    # Fix 6c: LLM post-hoc validation
    # ------------------------------------------------------------------

    def _llm_validate(
        self, phenotypes: List[PhenotypeMatch], original_text: str
    ) -> List[PhenotypeMatch]:
        """
        Fix 6c — Post-hoc LLM validation of negation and modifier decisions.

        Only runs when an API key is configured (self._api_key is not None).
        For each excluded phenotype, confirms it is truly excluded; if the LLM
        disagrees, flips it to present.  For each modifier, confirms it applies;
        if the LLM disagrees, removes the modifier.
        """
        if not self._api_key:
            return phenotypes

        validated: List[PhenotypeMatch] = []
        for pm in phenotypes:
            if self.cfg.llm_validate_negation and pm.excluded:
                confirmed = llm_check_negation(
                    pm.label,
                    original_text,
                    api_key=self._api_key,
                    model=self._llm_model,
                    cache=self._llm_cache,
                )
                if not confirmed:
                    # LLM says term is NOT excluded — flip to present
                    pm = dc_replace(pm, excluded=False)

            if (
                self.cfg.llm_validate_modifiers
                and not pm.excluded
                and pm.severity_id
                and pm.severity_label
            ):
                confirmed = llm_check_modifier(
                    pm.label,
                    pm.severity_label,
                    original_text,
                    api_key=self._api_key,
                    model=self._llm_model,
                    cache=self._llm_cache,
                )
                if not confirmed:
                    pm = dc_replace(pm, severity_id=None, severity_label=None)

            validated.append(pm)
        return validated

    # ------------------------------------------------------------------
    # Build PhenotypeOutput from tuple set
    # ------------------------------------------------------------------

    def _build_output(self, O: Set[_MatchTuple], raw_input: str) -> PhenotypeOutput:
        idx = self._idx
        phenotypes: List[PhenotypeMatch] = []
        diseases: List[DiseaseMatch] = []
        seen_pheno: Set[str] = set()
        seen_disease: Set[str] = set()

        # Fix 2b: Compute parent subsumption — suppress a term if a more
        # specific (descendant) matched term is also present and not negated.
        all_hp_present = {
            term_id
            for (term_id, _, negated) in O
            if term_id.startswith("HP:") and term_id in idx.phenotype_ids and not negated
        }
        suppressed: Set[str] = set()
        for tid in all_hp_present:
            for other in all_hp_present:
                if other != tid and tid in idx.hpo_ancestors.get(other, set()):
                    # tid is an ancestor of `other` → suppress the parent
                    suppressed.add(tid)
                    break

        for (term_id, mod_id, negated) in O:
            if term_id.startswith("HP:") and term_id in idx.phenotype_ids:
                # Fix 2b: Skip terms subsumed by a more specific matched term
                if not negated and term_id in suppressed:
                    continue
                if term_id in seen_pheno:
                    continue
                seen_pheno.add(term_id)
                label = idx.primary_label.get(term_id, term_id)

                # Fix 1b: Suppress modifier when severity is already encoded
                # in the term label (e.g. "Severe intellectual disability").
                if mod_id and any(w in label.lower() for w in _SEVERITY_WORDS):
                    mod_id = None

                sev_label = idx.primary_label.get(mod_id) if mod_id else None
                phenotypes.append(PhenotypeMatch(
                    hpo_id=term_id,
                    label=label,
                    excluded=negated,
                    severity_id=mod_id,
                    severity_label=sev_label,
                    confidence=1.0,
                    matched_by="v2",
                ))

            elif term_id.startswith("MONDO:"):
                key = term_id
                if key in seen_disease:
                    continue
                seen_disease.add(key)
                label = idx.primary_label.get(term_id, term_id)
                orphanet_ids = idx.mondo_to_orphanet.get(term_id, [])
                diseases.append(DiseaseMatch(
                    mondo_id=term_id,
                    mondo_label=label,
                    orphanet_ids=orphanet_ids,
                    confidence=1.0,
                    matched_by="v2",
                ))

            elif term_id.startswith("OMIM:"):
                key = term_id
                if key in seen_disease:
                    continue
                seen_disease.add(key)
                label = idx.primary_label.get(term_id, term_id)
                diseases.append(DiseaseMatch(
                    omim_phenotype_ids=[term_id],
                    omim_labels=[label],
                    confidence=1.0,
                    matched_by="v2",
                ))

            elif term_id.startswith("ORPHA:"):
                key = term_id
                if key in seen_disease:
                    continue
                seen_disease.add(key)
                diseases.append(DiseaseMatch(
                    orphanet_ids=[term_id],
                    confidence=1.0,
                    matched_by="v2",
                ))

        return PhenotypeOutput(
            phenotypes=phenotypes,
            diseases=diseases,
            raw_input=raw_input,
        )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> None:
    import argparse, json as _json

    parser = argparse.ArgumentParser(description="Phenotype Matcher v2")
    parser.add_argument("text", nargs="+", help="Clinical phenotype text")
    parser.add_argument("--hpo", default="ontology/hp.obo")
    parser.add_argument("--mondo", default="ontology/mondo.obo")
    parser.add_argument("--hpoa", default="data/phenotype.hpoa")
    parser.add_argument("--model", default="sapbert")
    parser.add_argument("--llm", default="deepseek")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    cfg = MatcherConfig(
        hpo_path=args.hpo,
        mondo_path=args.mondo,
        hpoa_path=args.hpoa,
        embedding_model=args.model,
        llm_model=args.llm,
        device=args.device,
        debug=args.debug,
    )
    matcher = PhenotypeMatcher(cfg)
    text = " ".join(args.text)
    result = matcher.match(text)
    print(_json.dumps(result.to_dict(), indent=2))


if __name__ == "__main__":
    main()
