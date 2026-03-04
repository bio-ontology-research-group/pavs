#!/usr/bin/env python3
"""
Tune NER extractor parameters (ner_peak_min_delta, ner_length_bonus) on
RAG-HPO CSC/GSC datasets.

Grid-searches parameter combinations, caching the expensive similarity
landscape computation so only peak detection + normalization need to rerun.

Usage:
    # Quick grid search on 20 cases:
    uv run python onto-normalize/tune_ner.py --limit 20

    # Full grid search:
    uv run python onto-normalize/tune_ner.py --dataset csc

    # Manual exploration with specific parameters:
    uv run python onto-normalize/tune_ner.py --peak-min-delta 0.03 --length-bonus 0.01 --limit 10

    # Visualize similarity landscape for a single case:
    uv run python onto-normalize/tune_ner.py --visualize --case-idx 0 --dataset csc
"""

import argparse
import json
import sys
import os
import time
import math
from pathlib import Path
from itertools import product

import numpy as np

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "tools" / "phenotype_matcher_v2" / "src"))

from phenotype_matcher_v2 import PhenotypeMatcher, MatcherConfig
from phenotype_matcher_v2.ner_extractor import NerExtractor, _sentence_split


def load_benchmark(cases_path: str, annotations_path: str):
    """Load benchmark cases and gold standard annotations."""
    with open(cases_path) as f:
        cases = json.load(f)
    with open(annotations_path) as f:
        annotations = json.load(f)
    return cases, annotations


def evaluate_case(predicted_hpo_ids: set, gold_hpo_ids: set):
    """Compute precision, recall, F1 for a single case."""
    tp = len(predicted_hpo_ids & gold_hpo_ids)
    fp = len(predicted_hpo_ids - gold_hpo_ids)
    fn = len(gold_hpo_ids - predicted_hpo_ids)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {"tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1}


def run_grid_search(
    cases, annotations, matcher, peak_deltas, length_bonuses, no_llm=True
):
    """
    Grid search over (peak_min_delta, length_bonus) combinations.

    Returns a list of dicts with parameters and macro-averaged P/R/F1.
    """
    results = []

    for delta, bonus in product(peak_deltas, length_bonuses):
        # Update config
        matcher.cfg.ner_peak_min_delta = delta
        matcher.cfg.ner_length_bonus = bonus

        # Reset NER extractor so it picks up new config
        if hasattr(matcher, "_ner_extractor"):
            matcher._ner_extractor = None

        total_tp, total_fp, total_fn = 0, 0, 0
        case_metrics = []

        for case in cases:
            case_id = case["case_id"]
            clinical_note = case["clinical_note"]
            gold_hpo_ids = set(annotations.get(case_id, []))
            if not gold_hpo_ids:
                continue

            try:
                output = matcher.match_narrative(clinical_note)
                predicted = set(output.get_hpo_ids())
            except Exception as e:
                print(f"  ERROR {case_id}: {e}", file=sys.stderr)
                predicted = set()

            metrics = evaluate_case(predicted, gold_hpo_ids)
            case_metrics.append(metrics)
            total_tp += metrics["tp"]
            total_fp += metrics["fp"]
            total_fn += metrics["fn"]

        n = len(case_metrics)
        if n == 0:
            continue

        micro_p = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0.0
        micro_r = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0.0
        micro_f1 = 2 * micro_p * micro_r / (micro_p + micro_r) if (micro_p + micro_r) > 0 else 0.0

        macro_p = sum(m["precision"] for m in case_metrics) / n
        macro_r = sum(m["recall"] for m in case_metrics) / n
        macro_f1 = sum(m["f1"] for m in case_metrics) / n

        result = {
            "peak_min_delta": delta,
            "length_bonus": bonus,
            "n_cases": n,
            "micro_precision": round(micro_p, 4),
            "micro_recall": round(micro_r, 4),
            "micro_f1": round(micro_f1, 4),
            "macro_precision": round(macro_p, 4),
            "macro_recall": round(macro_r, 4),
            "macro_f1": round(macro_f1, 4),
        }
        results.append(result)

        print(
            f"  delta={delta:.3f} bonus={bonus:.3f}  "
            f"micro P={micro_p:.3f} R={micro_r:.3f} F1={micro_f1:.3f}  "
            f"macro P={macro_p:.3f} R={macro_r:.3f} F1={macro_f1:.3f}",
            flush=True,
        )

    return results


def visualize_landscape(case, matcher):
    """
    Dump the similarity landscape for a single case to stdout.

    Shows the 2D grid of max cosine similarities for each (window_size, position)
    and highlights detected peaks.
    """
    clinical_note = case["clinical_note"]
    cfg = matcher.cfg
    idx = matcher._idx

    # Ensure encoder is loaded
    if matcher._encoder is None:
        matcher._batch_encode(["test"])

    sentences = _sentence_split(clinical_note)

    print(f"\nCase: {case['case_id']}")
    print(f"Text length: {len(clinical_note)} chars, {len(clinical_note.split())} words")
    print(f"Sentences: {len(sentences)}")
    print(f"Config: peak_min_delta={cfg.ner_peak_min_delta}, length_bonus={cfg.ner_length_bonus}")
    print()

    for si, sent in enumerate(sentences):
        words = sent.text.split()
        n_words = len(words)
        if n_words == 0:
            continue

        print(f"--- Sentence {si+1}: {sent.text[:100]}{'...' if len(sent.text) > 100 else ''}")
        print(f"    Words: {n_words}")

        min_w = cfg.ner_min_ngram_words
        max_w = min(cfg.ner_max_ngram_words, n_words)
        if min_w > max_w:
            print("    (too short for n-gram analysis)")
            continue

        pheno_words = idx.phenotype_word_set

        # Generate n-grams
        ngram_texts = []
        ngram_coords = []

        for w in range(min_w, max_w + 1):
            for i in range(n_words - w + 1):
                ngram = " ".join(words[i:i + w])
                ngram_lower_words = ngram.lower().split()
                if not any(wd in pheno_words for wd in ngram_lower_words):
                    continue
                ngram_texts.append(ngram)
                ngram_coords.append((w - min_w, i))

        if not ngram_texts:
            print("    (no phenotype-relevant n-grams)")
            continue

        # Batch encode
        ngram_embeddings = matcher._batch_encode(ngram_texts)
        if ngram_embeddings is None:
            continue

        # Compute similarities
        sims = ngram_embeddings @ idx.embeddings.T
        max_sims = np.max(sims, axis=1)
        best_hpo_idx = np.argmax(sims, axis=1)

        # Build landscape
        n_window_sizes = max_w - min_w + 1
        landscape = np.full((n_window_sizes, n_words), -np.inf)

        for m, (wi, pi) in enumerate(ngram_coords):
            landscape[wi, pi] = float(max_sims[m])

        valid_sims = max_sims[np.isfinite(max_sims)]
        if len(valid_sims) == 0:
            continue

        median_sim = float(np.median(valid_sims))
        noise_floor = median_sim + cfg.ner_peak_min_delta

        print(f"    Median sim: {median_sim:.3f}, Noise floor: {noise_floor:.3f}")
        print(f"    N-grams after filter: {len(ngram_texts)}")

        # Print landscape
        print(f"    {'w\\pos':<6}", end="")
        for pi in range(min(n_words, 20)):
            print(f"{pi:>6}", end="")
        print()

        for wi in range(n_window_sizes):
            w = wi + min_w
            print(f"    w={w:<3d} ", end="")
            for pi in range(min(n_words, 20)):
                val = landscape[wi, pi]
                if val <= -np.inf:
                    print(f"{'  -  ':>6}", end="")
                elif val >= noise_floor:
                    print(f"{val:>6.2f}", end="")
                else:
                    print(f"{val:>6.2f}", end="")
            print()

        # Show top peaks
        scored = np.copy(landscape)
        for wi in range(n_window_sizes):
            w = wi + min_w
            bonus = cfg.ner_length_bonus * math.log(max(w, 1))
            mask = landscape[wi] > -np.inf
            scored[wi, mask] += bonus

        peaks = []
        for wi in range(n_window_sizes):
            for pi in range(n_words):
                val = scored[wi, pi]
                if val <= -np.inf or landscape[wi, pi] < noise_floor:
                    continue
                is_peak = True
                for dw in (-1, 0, 1):
                    for dp in (-1, 0, 1):
                        if dw == 0 and dp == 0:
                            continue
                        nw, np_ = wi + dw, pi + dp
                        if 0 <= nw < n_window_sizes and 0 <= np_ < n_words:
                            if scored[nw, np_] > val:
                                is_peak = False
                                break
                    if not is_peak:
                        break
                if is_peak:
                    w = wi + min_w
                    ngram = " ".join(words[pi:pi + w])
                    # Find best HPO match
                    for m, (cwi, cpi) in enumerate(ngram_coords):
                        if cwi == wi and cpi == pi:
                            hpo_i = int(best_hpo_idx[m])
                            hpo_label = idx.embedding_labels[hpo_i]
                            hpo_ids = idx.exact_map.get(hpo_label, [])
                            peaks.append((val, landscape[wi, pi], ngram, hpo_label, hpo_ids))
                            break

        if peaks:
            peaks.sort(key=lambda x: x[0], reverse=True)
            print(f"\n    Peaks ({len(peaks)}):")
            for score, raw_sim, ngram, hpo_label, hpo_ids in peaks[:10]:
                ids_str = ",".join(hpo_ids[:2])
                print(f"      score={score:.3f} sim={raw_sim:.3f} "
                      f"'{ngram}' → '{hpo_label}' [{ids_str}]")
        print()


def main():
    parser = argparse.ArgumentParser(description="Tune NER extractor parameters")
    parser.add_argument("--no-llm", action="store_true", default=True,
                        help="Disable LLM validation (default: True)")
    parser.add_argument("--limit", type=int, default=0,
                        help="Limit number of cases (0=all)")
    parser.add_argument("--dataset", choices=["csc", "gsc", "both"], default="csc",
                        help="Which dataset to run")

    # Manual parameter exploration
    parser.add_argument("--peak-min-delta", type=float, nargs="+", default=None,
                        help="Peak min delta values to try (space-separated)")
    parser.add_argument("--length-bonus", type=float, nargs="+", default=None,
                        help="Length bonus values to try (space-separated)")

    # Visualization mode
    parser.add_argument("--visualize", action="store_true",
                        help="Visualize similarity landscape for a single case")
    parser.add_argument("--case-idx", type=int, default=0,
                        help="Case index for --visualize mode")

    args = parser.parse_args()

    base_dir = Path(__file__).parent

    print("Initializing phenotype_matcher_v2 [offline, NER tuning]...", flush=True)
    cfg = MatcherConfig(
        hpo_path=str(project_root / "ontology" / "hp.obo"),
        mondo_path=str(project_root / "ontology" / "mondo.obo"),
        hpoa_path=str(project_root / "data" / "phenotype.hpoa"),
        embedding_model="sapbert",
        device="cpu",
        api_key=None,
        llm_model="deepseek",
        debug=False,
        match_cache_size=10_000,
    )
    matcher = PhenotypeMatcher(cfg)
    print("Matcher initialized.", flush=True)

    # Default grid
    if args.peak_min_delta is None:
        peak_deltas = [0.01, 0.02, 0.03, 0.05, 0.07, 0.10]
    else:
        peak_deltas = args.peak_min_delta

    if args.length_bonus is None:
        length_bonuses = [0.00, 0.01, 0.02, 0.03, 0.05]
    else:
        length_bonuses = args.length_bonus

    for ds_key in (["csc", "gsc"] if args.dataset == "both" else [args.dataset]):
        cases, annotations = load_benchmark(
            str(base_dir / f"{ds_key}_cases.json"),
            str(base_dir / f"{ds_key}_annotations.json"),
        )
        if args.limit > 0:
            cases = cases[:args.limit]

        if args.visualize:
            idx = min(args.case_idx, len(cases) - 1)
            visualize_landscape(cases[idx], matcher)
            continue

        print(f"\n{'='*60}")
        print(f"Grid search on {ds_key.upper()}: {len(cases)} cases")
        print(f"Parameters: delta={peak_deltas}, bonus={length_bonuses}")
        print(f"Total combinations: {len(peak_deltas) * len(length_bonuses)}")
        print(f"{'='*60}\n")

        start = time.time()
        grid_results = run_grid_search(
            cases, annotations, matcher, peak_deltas, length_bonuses
        )
        elapsed = time.time() - start

        print(f"\nGrid search completed in {elapsed:.1f}s")

        # Find best by macro F1
        if grid_results:
            best = max(grid_results, key=lambda r: r["macro_f1"])
            print(f"\nBest macro F1: {best['macro_f1']:.4f}")
            print(f"  delta={best['peak_min_delta']:.3f}, bonus={best['length_bonus']:.3f}")
            print(f"  P={best['macro_precision']:.4f} R={best['macro_recall']:.4f}")

            # Find best precision at recall >= 0.40
            high_recall = [r for r in grid_results if r["macro_recall"] >= 0.40]
            if high_recall:
                best_prec = max(high_recall, key=lambda r: r["macro_precision"])
                print(f"\nBest precision (recall >= 0.40): P={best_prec['macro_precision']:.4f}")
                print(f"  delta={best_prec['peak_min_delta']:.3f}, bonus={best_prec['length_bonus']:.3f}")
                print(f"  R={best_prec['macro_recall']:.4f} F1={best_prec['macro_f1']:.4f}")

        # Save results
        output_file = base_dir / f"tune_ner_{ds_key}_results.json"
        with open(output_file, "w") as f:
            json.dump(grid_results, f, indent=2)
        print(f"\nResults saved to: {output_file}")


if __name__ == "__main__":
    main()
