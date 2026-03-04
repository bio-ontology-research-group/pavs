#!/usr/bin/env python3
"""
Benchmark phenotype_matcher_v2 on RAG-HPO CSC and GSC datasets.

Evaluates HPO term extraction from clinical narratives against gold standard
annotations, computing precision, recall, and F1 per case and overall.

Matches are evaluated as exact HPO ID matches (same as RAG-HPO paper methodology).
"""

import argparse
import json
import sys
import os
import time
from pathlib import Path
from collections import defaultdict

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "tools" / "phenotype_matcher_v2" / "src"))

from phenotype_matcher_v2 import PhenotypeMatcher, MatcherConfig


def load_benchmark(cases_path: str, annotations_path: str):
    """Load benchmark cases and gold standard annotations."""
    with open(cases_path) as f:
        cases = json.load(f)
    with open(annotations_path) as f:
        annotations = json.load(f)
    return cases, annotations


def evaluate_case(predicted_hpo_ids: set, gold_hpo_ids: set):
    """Compute TP, FP, FN, precision, recall, F1 for a single case."""
    tp = len(predicted_hpo_ids & gold_hpo_ids)
    fp = len(predicted_hpo_ids - gold_hpo_ids)
    fn = len(gold_hpo_ids - predicted_hpo_ids)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return {
        "tp": tp, "fp": fp, "fn": fn,
        "precision": precision, "recall": recall, "f1": f1
    }


def run_benchmark(dataset_name: str, cases: list, annotations: dict, matcher: PhenotypeMatcher, narrative: bool = False):
    """Run benchmark on a dataset and return per-case and aggregate results."""
    results = []
    total_tp, total_fp, total_fn = 0, 0, 0
    skipped = 0

    print(f"\n{'='*60}")
    print(f"Benchmarking on {dataset_name}: {len(cases)} cases")
    print(f"{'='*60}", flush=True)

    _bench_start = time.time()
    for i, case in enumerate(cases):
        case_id = case["case_id"]
        clinical_note = case["clinical_note"]

        # Get gold standard
        gold_hpo_ids = set(annotations.get(case_id, []))
        if not gold_hpo_ids:
            skipped += 1
            continue

        # Run matcher
        try:
            match_fn = matcher.match_narrative if narrative else matcher.match
            output = match_fn(clinical_note)
            predicted_hpo_ids = set(output.get_hpo_ids())  # present (non-excluded) HPO IDs
        except Exception as e:
            print(f"  ERROR on case {case_id}: {e}")
            predicted_hpo_ids = set()

        # Evaluate
        metrics = evaluate_case(predicted_hpo_ids, gold_hpo_ids)
        metrics["case_id"] = case_id
        metrics["n_predicted"] = len(predicted_hpo_ids)
        metrics["n_gold"] = len(gold_hpo_ids)
        results.append(metrics)

        total_tp += metrics["tp"]
        total_fp += metrics["fp"]
        total_fn += metrics["fn"]

        # Per-case progress
        elapsed = time.time() - _bench_start
        avg_per_case = elapsed / (len(results)) if results else 0
        remaining = avg_per_case * (len(cases) - i - 1)
        print(f"  [{i+1:3d}/{len(cases)}] case={case_id:>4s}  "
              f"P={metrics['precision']:.2f} R={metrics['recall']:.2f} F1={metrics['f1']:.2f}  "
              f"pred={len(predicted_hpo_ids):2d} gold={len(gold_hpo_ids):2d}  "
              f"[{elapsed:.0f}s elapsed, ~{remaining:.0f}s left]",
              flush=True)

    # Micro-averaged metrics
    micro_precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0.0
    micro_recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0.0
    micro_f1 = 2 * micro_precision * micro_recall / (micro_precision + micro_recall) if (micro_precision + micro_recall) > 0 else 0.0

    # Macro-averaged metrics
    n = len(results)
    macro_precision = sum(r["precision"] for r in results) / n if n > 0 else 0.0
    macro_recall = sum(r["recall"] for r in results) / n if n > 0 else 0.0
    macro_f1 = sum(r["f1"] for r in results) / n if n > 0 else 0.0

    summary = {
        "dataset": dataset_name,
        "n_cases": len(results),
        "skipped": skipped,
        "total_tp": total_tp,
        "total_fp": total_fp,
        "total_fn": total_fn,
        "micro_precision": round(micro_precision, 4),
        "micro_recall": round(micro_recall, 4),
        "micro_f1": round(micro_f1, 4),
        "macro_precision": round(macro_precision, 4),
        "macro_recall": round(macro_recall, 4),
        "macro_f1": round(macro_f1, 4),
    }

    print(f"\n--- {dataset_name} Summary ---")
    print(f"Cases evaluated: {len(results)} (skipped: {skipped})")
    print(f"Total: TP={total_tp}, FP={total_fp}, FN={total_fn}")
    print(f"Micro-avg: P={micro_precision:.4f} R={micro_recall:.4f} F1={micro_f1:.4f}")
    print(f"Macro-avg: P={macro_precision:.4f} R={macro_recall:.4f} F1={macro_f1:.4f}")

    return results, summary


def print_comparison(our_results: dict, published_path: str):
    """Print comparison table with published results."""
    with open(published_path) as f:
        published = json.load(f)

    for dataset_key, dataset_label in [("csc", "CSC (Case Study Cohort)"), ("gsc", "GSC (Gold Standard Corpora)")]:
        print(f"\n{'='*70}")
        print(f"Comparison: {dataset_label}")
        print(f"{'='*70}")
        print(f"{'Tool':<35s} {'Precision':>10s} {'Recall':>10s} {'F1':>10s}")
        print("-" * 70)

        # Our results (macro-averaged to match RAG-HPO paper methodology)
        if dataset_key in our_results:
            r = our_results[dataset_key]
            print(f"{'** phenotype_matcher_v2 **':<35s} {r['macro_precision']:>10.4f} {r['macro_recall']:>10.4f} {r['macro_f1']:>10.4f}")

        # Published results
        if dataset_key in published:
            for tool, r in published[dataset_key].items():
                print(f"{tool:<35s} {r['precision']:>10.4f} {r['recall']:>10.4f} {r['f1']:>10.4f}")


def main():
    parser = argparse.ArgumentParser(description="Benchmark phenotype_matcher_v2")
    parser.add_argument("--no-llm", action="store_true",
                        help="Disable LLM validation (faster, offline-only)")
    parser.add_argument("--limit", type=int, default=0,
                        help="Limit number of cases per dataset (0=all)")
    parser.add_argument("--dataset", choices=["csc", "gsc", "both"], default="both",
                        help="Which dataset to run")
    parser.add_argument("--narrative", action="store_true",
                        help="Use NER-based narrative mode (match_narrative)")
    args = parser.parse_args()

    base_dir = Path(__file__).parent

    api_key = None if args.no_llm else os.environ.get("OPENROUTER_API_KEY")
    mode = "offline (no LLM)" if args.no_llm or not api_key else "with LLM validation"
    narrative_str = " + narrative NER" if args.narrative else ""
    print(f"Initializing phenotype_matcher_v2 [{mode}{narrative_str}]...", flush=True)
    cfg = MatcherConfig(
        hpo_path=str(project_root / "ontology" / "hp.obo"),
        mondo_path=str(project_root / "ontology" / "mondo.obo"),
        hpoa_path=str(project_root / "data" / "phenotype.hpoa"),
        embedding_model="sapbert",
        device="cpu",
        api_key=api_key,
        llm_model="deepseek",
        debug=False,
        match_cache_size=10_000,
    )
    matcher = PhenotypeMatcher(cfg)
    print("Matcher initialized.", flush=True)

    our_results = {}
    all_per_case = {}

    if args.dataset in ("csc", "both"):
        csc_cases, csc_annotations = load_benchmark(
            str(base_dir / "csc_cases.json"),
            str(base_dir / "csc_annotations.json"),
        )
        if args.limit > 0:
            csc_cases = csc_cases[:args.limit]
        start = time.time()
        csc_per_case, csc_summary = run_benchmark("CSC", csc_cases, csc_annotations, matcher, narrative=args.narrative)
        csc_time = time.time() - start
        print(f"CSC time: {csc_time:.1f}s ({csc_time/len(csc_per_case):.1f}s/case)")
        our_results["csc"] = csc_summary
        all_per_case["csc_per_case"] = csc_per_case

    if args.dataset in ("gsc", "both"):
        gsc_cases, gsc_annotations = load_benchmark(
            str(base_dir / "gsc_cases.json"),
            str(base_dir / "gsc_annotations.json"),
        )
        if args.limit > 0:
            gsc_cases = gsc_cases[:args.limit]
        start = time.time()
        gsc_per_case, gsc_summary = run_benchmark("GSC", gsc_cases, gsc_annotations, matcher, narrative=args.narrative)
        gsc_time = time.time() - start
        print(f"GSC time: {gsc_time:.1f}s ({gsc_time/len(gsc_per_case):.1f}s/case)")
        our_results["gsc"] = gsc_summary
        all_per_case["gsc_per_case"] = gsc_per_case

    # Save detailed results
    with open(base_dir / "benchmark_results.json", "w") as f:
        json.dump({
            "summary": our_results,
            **all_per_case,
        }, f, indent=2)

    # Print comparison
    print_comparison(our_results, str(base_dir / "published_results.json"))

    print(f"\nDetailed results saved to: {base_dir / 'benchmark_results.json'}")


if __name__ == "__main__":
    main()
