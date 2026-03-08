"""
NER extractor for clinical narratives.

Identifies candidate phenotype-bearing text spans from long clinical text
using a combination of:
  1. Aho-Corasick exact substring matching (reuses OntologyIndex.ac_automaton)
  2. Sliding n-gram window + embedding similarity landscape + 2D peak detection

Produces *text spans only* — all normalization, negation, modifier detection,
and LLM validation happen downstream in PhenotypeMatcher.match().
"""

import math
import re
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Set, Tuple

import numpy as np

from .detection import is_boundary_valid
from .ontology_index import OntologyIndex
from .schemas import MatcherConfig


@dataclass
class SentenceSpan:
    """A sentence with its character offsets in the original text."""
    text: str
    start: int  # char offset in original text
    end: int


@dataclass
class CandidateSpan:
    """A candidate phenotype-bearing text span."""
    text: str
    source: str  # "ac" or "peak"


# Simple sentence boundary pattern: split on period followed by space+uppercase,
# or on newlines. This avoids requiring spaCy just for sentence splitting.
_SENT_SPLIT = re.compile(r'(?<=[.!?])\s+(?=[A-Z])|[\n\r]+')


def _sentence_split(text: str) -> List[SentenceSpan]:
    """Split text into sentence spans using regex heuristics."""
    spans: List[SentenceSpan] = []
    parts = _SENT_SPLIT.split(text)
    offset = 0
    for part in parts:
        part_stripped = part.strip()
        if not part_stripped:
            offset += len(part)
            continue
        # Find actual start in original text
        idx = text.find(part_stripped, offset)
        if idx == -1:
            idx = offset
        spans.append(SentenceSpan(
            text=part_stripped,
            start=idx,
            end=idx + len(part_stripped),
        ))
        offset = idx + len(part_stripped)
    return spans


class NerExtractor:
    """
    Extract candidate phenotype spans from clinical narrative text.

    Uses two complementary strategies:
    - AC automaton scan: catches exact HPO label substrings in running prose
    - Sliding n-gram window: computes embedding similarity landscape and
      detects peaks that signal phenotype mentions even below the
      normalization threshold
    """

    def __init__(
        self,
        ontology_index: OntologyIndex,
        config: MatcherConfig,
        encoder_fn: Callable[[List[str]], np.ndarray],
    ):
        """
        Parameters
        ----------
        ontology_index : OntologyIndex
            Pre-built ontology index with ac_automaton, embeddings, etc.
        config : MatcherConfig
            Configuration with NER parameters.
        encoder_fn : callable
            Batch-encodes list of strings → normalized numpy matrix (N × D).
        """
        self._idx = ontology_index
        self._cfg = config
        self._encode = encoder_fn

    def extract(self, text: str) -> List[CandidateSpan]:
        """
        Extract candidate phenotype spans from clinical narrative text.

        Returns a deduplicated list of CandidateSpan objects.
        """
        sentences = _sentence_split(text)
        candidates: List[CandidateSpan] = []

        for sent in sentences:
            # Strategy 1: AC automaton scan
            candidates.extend(self._ac_scan(sent.text))

            # Strategy 2: Sliding n-gram similarity landscape + peak detection
            candidates.extend(self._ngram_peaks(sent.text))

        # Deduplicate
        return self._deduplicate(candidates)

    # ------------------------------------------------------------------
    # Strategy 1: AC automaton scan
    # ------------------------------------------------------------------

    def _ac_scan(self, sentence: str) -> List[CandidateSpan]:
        """Run Aho-Corasick automaton on a sentence to find exact label matches."""
        idx = self._idx
        if idx.ac_automaton is None:
            return []

        results: List[CandidateSpan] = []
        sent_lower = sentence.lower()

        for end_idx, label in idx.ac_automaton.iter(sent_lower):
            start_idx = end_idx - len(label) + 1
            if not is_boundary_valid(sent_lower, start_idx, end_idx + 1):
                continue
            # Extract the original-case text from the sentence
            span_text = sentence[start_idx:end_idx + 1]
            results.append(CandidateSpan(text=span_text, source="ac"))

        return results

    # ------------------------------------------------------------------
    # Strategy 2: Sliding n-gram + similarity landscape + peak detection
    # ------------------------------------------------------------------

    def _ngram_peaks(self, sentence: str) -> List[CandidateSpan]:
        """
        Generate n-grams, compute similarity landscape, detect peaks.

        The similarity landscape is 2D: position (x) × window size (y).
        Peak detection finds local maxima that stand out above the sentence
        median similarity, favoring longer (more specific) matches.
        """
        idx = self._idx
        cfg = self._cfg

        if idx.embeddings is None:
            return []

        words = sentence.split()
        n_words = len(words)
        if n_words == 0:
            return []

        min_w = cfg.ner_min_ngram_words
        max_w = min(cfg.ner_max_ngram_words, n_words)
        if min_w > max_w:
            return []

        pheno_words = idx.phenotype_word_set

        # Generate n-grams and track their grid coordinates
        ngram_texts: List[str] = []
        ngram_coords: List[Tuple[int, int]] = []  # (window_size_idx, start_pos)

        for w in range(min_w, max_w + 1):
            for i in range(n_words - w + 1):
                ngram = " ".join(words[i:i + w])
                # Pre-filter: at least one word must be in phenotype_word_set
                ngram_lower_words = ngram.lower().split()
                if not any(wd in pheno_words for wd in ngram_lower_words):
                    continue
                ngram_texts.append(ngram)
                ngram_coords.append((w - min_w, i))

        if not ngram_texts:
            return []

        # Batch-encode all surviving n-grams
        ngram_embeddings = self._encode(ngram_texts)  # (M × D)
        if ngram_embeddings is None or len(ngram_embeddings) == 0:
            return []

        # Compute max similarity to any HPO label for each n-gram
        # sims shape: (M × N) where N = number of HPO label embeddings
        sims = ngram_embeddings @ idx.embeddings.T
        max_sims = np.max(sims, axis=1)         # (M,)
        best_hpo_idx = np.argmax(sims, axis=1)  # (M,)

        # Build 2D landscape grid
        n_window_sizes = max_w - min_w + 1
        n_positions = n_words  # max possible start position (some cols may be empty)

        # Initialize with -inf (no data)
        landscape = np.full((n_window_sizes, n_positions), -np.inf)
        hpo_idx_grid = np.full((n_window_sizes, n_positions), -1, dtype=int)
        text_grid: Dict[Tuple[int, int], str] = {}

        for m, (wi, pi) in enumerate(ngram_coords):
            landscape[wi, pi] = float(max_sims[m])
            hpo_idx_grid[wi, pi] = int(best_hpo_idx[m])
            text_grid[(wi, pi)] = ngram_texts[m]

        # Length-biased scoring
        scored = np.copy(landscape)
        for wi in range(n_window_sizes):
            w = wi + min_w
            bonus = cfg.ner_length_bonus * math.log(max(w, 1))
            mask = landscape[wi] > -np.inf
            scored[wi, mask] += bonus

        # Adaptive threshold: median of all valid max_sim values
        valid_sims = max_sims[np.isfinite(max_sims)]
        if len(valid_sims) == 0:
            return []
        noise_floor = float(np.median(valid_sims)) + cfg.ner_peak_min_delta

        # 2D peak detection: find local maxima in 8-connected neighborhood
        peaks: List[Tuple[int, int, float]] = []  # (wi, pi, scored_value)

        for wi in range(n_window_sizes):
            for pi in range(n_positions):
                val = scored[wi, pi]
                if val <= -np.inf:
                    continue
                # Must exceed noise floor (based on raw similarity, not scored)
                if landscape[wi, pi] < noise_floor:
                    continue

                # Check 8-connected neighbors
                is_peak = True
                for dw in (-1, 0, 1):
                    for dp in (-1, 0, 1):
                        if dw == 0 and dp == 0:
                            continue
                        nw, np_ = wi + dw, pi + dp
                        if 0 <= nw < n_window_sizes and 0 <= np_ < n_positions:
                            if scored[nw, np_] > val:
                                is_peak = False
                                break
                    if not is_peak:
                        break

                if is_peak:
                    peaks.append((wi, pi, val))

        if not peaks:
            return []

        # Specificity preference: when peaks overlap in text position,
        # prefer the one whose best HPO match is more specific (descendant)
        peaks_by_pos = self._resolve_overlapping_peaks(
            peaks, hpo_idx_grid, text_grid, min_w, idx
        )

        # Convert surviving peaks to CandidateSpans
        results: List[CandidateSpan] = []
        for wi, pi in peaks_by_pos:
            ngram_text = text_grid.get((wi, pi))
            if ngram_text:
                results.append(CandidateSpan(text=ngram_text, source="peak"))

        return results

    def _resolve_overlapping_peaks(
        self,
        peaks: List[Tuple[int, int, float]],
        hpo_idx_grid: np.ndarray,
        text_grid: Dict[Tuple[int, int], str],
        min_w: int,
        idx: OntologyIndex,
    ) -> List[Tuple[int, int]]:
        """
        Resolve overlapping peaks using ontology specificity.

        When two peaks overlap in word positions, prefer the one whose
        best-matching HPO term is more specific (descendant). If neither
        subsumes the other, keep both.
        """
        # Sort peaks by scored value (descending) for greedy selection
        peaks_sorted = sorted(peaks, key=lambda x: x[2], reverse=True)
        surviving: List[Tuple[int, int]] = []
        covered_positions: List[Tuple[int, int, int, int]] = []  # (wi, pi, start_word, end_word)

        for wi, pi, score in peaks_sorted:
            w = wi + min_w
            start_word = pi
            end_word = pi + w  # exclusive

            # Check overlap with already-accepted peaks
            overlaps_with = None
            for j, (owi, opi, ostart, oend) in enumerate(covered_positions):
                if start_word < oend and end_word > ostart:
                    overlaps_with = j
                    break

            if overlaps_with is None:
                # No overlap — accept
                surviving.append((wi, pi))
                covered_positions.append((wi, pi, start_word, end_word))
            else:
                # Overlap — check ontology specificity
                owi, opi, ostart, oend = covered_positions[overlaps_with]
                existing_hpo_i = int(hpo_idx_grid[owi, opi])
                new_hpo_i = int(hpo_idx_grid[wi, pi])

                if existing_hpo_i < 0 or new_hpo_i < 0:
                    continue  # no HPO match info

                existing_label = idx.embedding_labels[existing_hpo_i]
                new_label = idx.embedding_labels[new_hpo_i]

                existing_ids = idx.exact_map.get(existing_label, [])
                new_ids = idx.exact_map.get(new_label, [])

                # Check if new peak's HPO term is more specific
                new_is_more_specific = False
                existing_is_more_specific = False

                for nid in new_ids:
                    if not nid.startswith("HP:"):
                        continue
                    ancestors = idx.hpo_ancestors.get(nid, set())
                    for eid in existing_ids:
                        if eid in ancestors:
                            new_is_more_specific = True
                            break
                    if new_is_more_specific:
                        break

                if not new_is_more_specific:
                    for eid in existing_ids:
                        if not eid.startswith("HP:"):
                            continue
                        ancestors = idx.hpo_ancestors.get(eid, set())
                        for nid in new_ids:
                            if nid in ancestors:
                                existing_is_more_specific = True
                                break
                        if existing_is_more_specific:
                            break

                if new_is_more_specific and not existing_is_more_specific:
                    # Replace existing with new (more specific)
                    surviving[overlaps_with] = (wi, pi)
                    covered_positions[overlaps_with] = (wi, pi, start_word, end_word)
                elif not new_is_more_specific and not existing_is_more_specific:
                    # Different branches — keep both
                    surviving.append((wi, pi))
                    covered_positions.append((wi, pi, start_word, end_word))
                # else: existing is more specific or equal — keep existing

        return surviving

    # ------------------------------------------------------------------
    # Deduplication
    # ------------------------------------------------------------------

    @staticmethod
    def _deduplicate(candidates: List[CandidateSpan]) -> List[CandidateSpan]:
        """
        Deduplicate candidates:
        1. Exact text dedup (case-insensitive)
        2. Containment: if span A's text is a substring of span B's text, keep only B
        """
        if not candidates:
            return []

        # 1. Exact text dedup
        seen: Dict[str, CandidateSpan] = {}
        for c in candidates:
            key = c.text.lower().strip()
            if key not in seen:
                seen[key] = c

        unique = list(seen.values())

        # 2. Containment: remove spans that are substrings of other spans
        # Sort by text length (longest first)
        unique.sort(key=lambda c: len(c.text), reverse=True)
        result: List[CandidateSpan] = []
        result_texts_lower: List[str] = []

        for c in unique:
            c_lower = c.text.lower().strip()
            # Check if this text is contained in any already-accepted longer span
            is_contained = False
            for longer in result_texts_lower:
                if c_lower in longer:
                    is_contained = True
                    break
            if not is_contained:
                result.append(c)
                result_texts_lower.append(c_lower)

        return result
