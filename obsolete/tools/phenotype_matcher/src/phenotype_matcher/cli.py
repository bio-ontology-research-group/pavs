"""
Command-line interface for the Phenotype Matcher tool.
"""

import argparse
import json
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional, List
from pathlib import Path

from .matcher import PhenotypeMatcher
from .schemas import PhenotypeInput, MatcherConfig
from . import test_cases


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Phenotype Matcher - Map clinical phenotype descriptions to ontology identifiers",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Match a single phenotype
  phenotype-matcher "severe intellectual disability"
  
  # Match multiple separate phenotypes (batch mode)
  phenotype-matcher -b "seizures" -b "hypotonia" -b "feeding difficulties"
  
  # Batch processing from file
  phenotype-matcher --input phenotypes.txt --output results.jsonl
  
  # Parallel processing with 4 threads
  phenotype-matcher --input phenotypes.txt --output results.jsonl -j 4 --progress
  
  # Batch with different models
  phenotype-matcher -b "seizures" -b "hypotonia" --embedding-model medical --llm-model accurate
  
  # Output in different formats
  phenotype-matcher --input phenotypes.txt --output-format tsv --output results.tsv
  phenotype-matcher -b "seizures" -b "hypotonia" --output-format hpo_ids
  
  # Provide context for disambiguation
  phenotype-matcher "ASD" --context "cardiac defect"
  
  # Run all difficult test cases
  phenotype-matcher --test
  
  # Run a specific test case
  phenotype-matcher --test-case case_1
        """,
    )

    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("text", nargs="?", help="Phenotype description to match")
    input_group.add_argument(
        "-i",
        "--input",
        type=str,
        help="Input file with one phenotype description per line",
    )
    input_group.add_argument(
        "--test",
        action="store_true",
        help="Run difficult test cases to validate the matcher",
    )
    input_group.add_argument(
        "--test-case",
        type=str,
        help="Run a specific test case by ID (e.g., case_1, case_2)",
    )

    # Batch input option (can be used with text argument)
    parser.add_argument(
        "-b",
        "--batch",
        action="append",
        help="Add phenotype description to batch (can be used multiple times)",
    )

    # Output options
    parser.add_argument(
        "-o", "--output", type=str, help="Output file (default: stdout)"
    )
    parser.add_argument(
        "--output-format",
        choices=["json", "hpo_ids", "tsv", "summary"],
        default="json",
        help="Output format (default: json)",
    )

    # Context
    parser.add_argument(
        "-c",
        "--context",
        type=str,
        help="Additional context for disambiguation (e.g., patient info)",
    )
    parser.add_argument(
        "-g",
        "--gene-hint",
        type=str,
        help="Gene symbol for disambiguation (SOFT cue only, e.g., SCN1A, NKX2-5)",
    )
    parser.add_argument(
        "--split-by",
        type=str,
        default=",",
        help="Character to split input text on (default: comma)",
    )

    # Model configuration
    parser.add_argument(
        "--embedding-model",
        choices=["fast", "balanced", "accurate", "medical", "biobert", "sapbert"],
        default="balanced",
        help="Embedding model preset (default: balanced)",
    )
    parser.add_argument(
        "--llm-model",
        choices=["fast", "accurate", "balanced", "cheap"],
        default="fast",
        help="LLM model preset (default: fast)",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=5,
        help="Number of candidates to retrieve (default: 5)",
    )

    # Other options
    parser.add_argument(
        "--cache-dir",
        type=str,
        default="data/graph_rag_cache",
        help="Cache directory for embeddings (default: data/graph_rag_cache)",
    )
    parser.add_argument(
        "--device",
        choices=["cpu", "cuda"],
        default="cpu",
        help="Compute device (default: cpu)",
    )
    parser.add_argument(
        "--api-key",
        type=str,
        help="OpenRouter API key (default: use OPENROUTER_API_KEY env var)",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")

    # Advanced options
    parser.add_argument(
        "--no-ner",
        action="store_true",
        help="Disable LLM-based NER for term extraction (use simple comma splitting)",
    )
    parser.add_argument(
        "--no-acronym-expansion",
        action="store_true",
        help="Disable acronym expansion (e.g., ASD, DD, CHD)",
    )

    # Parallelization options
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel threads for batch processing (default: 1, no parallelization)",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="Show progress bar during batch processing",
    )

    args = parser.parse_args()

    # Create config
    config = MatcherConfig(
        embedding_model=args.embedding_model,
        llm_model=args.llm_model,
        top_k_phenotype=args.top_k,
        top_k_severity=args.top_k,
        top_k_disease=args.top_k,
        cache_dir=args.cache_dir,
        device=args.device,
        api_key=args.api_key,
        use_ner=not args.no_ner,
        expand_acronyms=not args.no_acronym_expansion,
        debug=args.debug,
    )

    # Handle test mode
    if args.test or args.test_case:
        run_test_cases(args.test_case, config)
        return

    # Collect all inputs to process
    inputs_to_process = []

    if args.input:
        # Load from file
        inputs_to_process = load_inputs_from_file(args.input)
    elif args.batch:
        # Multiple batch inputs from command line
        inputs_to_process = args.batch
        if args.text:
            # Include the positional argument too
            inputs_to_process.insert(0, args.text)
    elif args.text:
        # Single input
        inputs_to_process = [args.text]

    if not inputs_to_process:
        print("Error: No input provided", file=sys.stderr)
        sys.exit(1)

    # Process inputs (with optional parallelization)
    if args.jobs > 1 and len(inputs_to_process) > 1:
        # Parallel processing
        results = process_batch_parallel(
            inputs_to_process,
            config,
            args.context,
            args.gene_hint,
            args.split_by,
            args.jobs,
            args.output,
            args.output_format,
            args.progress,
        )
    else:
        # Sequential processing
        results = process_batch_sequential(
            inputs_to_process,
            config,
            args.context,
            args.gene_hint,
            args.split_by,
            args.output,
            args.output_format,
            args.progress,
        )

    # If no output file specified, print results
    if not args.output:
        output_text = format_output(results, args.output_format)
        print(output_text)


def load_inputs_from_file(input_file: str) -> List[str]:
    """Load phenotype descriptions from a file (one per line)."""
    inputs = []
    with open(input_file, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                inputs.append(line)
    return inputs


def write_result_to_file(
    output_file: str, result, format_type: str, file_lock: threading.Lock
):
    """Thread-safe write of a single result to output file."""
    # Format single result
    if format_type == "json":
        line = json.dumps(result.to_dict())
    elif format_type == "hpo_ids":
        line = "|".join(result.get_hpo_ids(excluded=False))
    elif format_type == "tsv":
        hpo_ids = "|".join(result.get_hpo_ids(excluded=False))
        hpo_labels = "|".join([p.label for p in result.get_present_phenotypes()])
        excluded_ids = "|".join(result.get_hpo_ids(excluded=True))
        severity_ids = "|".join(
            [p.severity_id for p in result.phenotypes if p.severity_id]
        )
        omim_ids = "|".join(result.get_omim_ids())
        mondo_ids = "|".join(result.get_mondo_ids())
        line = f"{result.raw_input}\t{hpo_ids}\t{hpo_labels}\t{excluded_ids}\t{severity_ids}\t{omim_ids}\t{mondo_ids}"
    else:
        line = json.dumps(result.to_dict())

    # Thread-safe write
    with file_lock:
        with open(output_file, "a") as f:
            f.write(line + "\n")


def process_single_input(
    text: str,
    matcher: PhenotypeMatcher,
    context: Optional[str],
    gene_hint: Optional[str],
    split_by: str,
) -> tuple:
    """Process a single input and return (index, result, error)."""
    try:
        input_data = PhenotypeInput(
            text=text,
            context=context,
            gene_hint=gene_hint,
            split_by=split_by if split_by != "none" else None,
        )
        result = matcher.match(input_data)
        return (text, result, None)
    except Exception as e:
        return (text, None, str(e))


def process_batch_sequential(
    inputs: List[str],
    config: MatcherConfig,
    context: Optional[str],
    gene_hint: Optional[str],
    split_by: str,
    output_file: Optional[str],
    output_format: str,
    show_progress: bool,
) -> list:
    """Process inputs sequentially with optional progress display."""
    # Initialize matcher
    matcher = PhenotypeMatcher(config)

    results = []
    file_lock = threading.Lock()

    # Initialize output file if specified
    if output_file:
        # Write header for TSV format
        if output_format == "tsv":
            with open(output_file, "w") as f:
                f.write(
                    "input\thpo_ids\thpo_labels\texcluded_ids\tseverity_ids\tomim_ids\tmondo_ids\n"
                )
        else:
            # Clear file
            open(output_file, "w").close()

    # Process with optional progress
    if show_progress:
        try:
            from tqdm import tqdm

            iterator = tqdm(inputs, desc="Processing", unit="phenotype")
        except ImportError:
            print("Warning: tqdm not installed, progress bar disabled", file=sys.stderr)
            iterator = inputs
    else:
        iterator = inputs

    for text in iterator:
        _, result, error = process_single_input(
            text, matcher, context, gene_hint, split_by
        )

        if error:
            print(f"Error processing '{text}': {error}", file=sys.stderr)
            continue

        results.append(result)

        # Write to file immediately if specified
        if output_file and result:
            write_result_to_file(output_file, result, output_format, file_lock)

    return results


def process_batch_parallel(
    inputs: List[str],
    config: MatcherConfig,
    context: Optional[str],
    gene_hint: Optional[str],
    split_by: str,
    num_jobs: int,
    output_file: Optional[str],
    output_format: str,
    show_progress: bool,
) -> list:
    """Process inputs in parallel with thread-safe file writing."""
    # Initialize matcher (shared by all threads)
    matcher = PhenotypeMatcher(config)

    results = []
    file_lock = threading.Lock()

    # Initialize output file if specified
    if output_file:
        # Write header for TSV format
        if output_format == "tsv":
            with open(output_file, "w") as f:
                f.write(
                    "input\thpo_ids\thpo_labels\texcluded_ids\tseverity_ids\tomim_ids\tmondo_ids\n"
                )
        else:
            # Clear file
            open(output_file, "w").close()

    print(f"Processing {len(inputs)} inputs with {num_jobs} parallel threads...")

    # Process in parallel
    with ThreadPoolExecutor(max_workers=num_jobs) as executor:
        # Submit all tasks
        future_to_input = {
            executor.submit(
                process_single_input, text, matcher, context, gene_hint, split_by
            ): text
            for text in inputs
        }

        # Process completed tasks with optional progress
        if show_progress:
            try:
                from tqdm import tqdm

                iterator = tqdm(
                    as_completed(future_to_input),
                    total=len(inputs),
                    desc="Processing",
                    unit="phenotype",
                )
            except ImportError:
                print(
                    "Warning: tqdm not installed, progress bar disabled",
                    file=sys.stderr,
                )
                iterator = as_completed(future_to_input)
        else:
            iterator = as_completed(future_to_input)

        for future in iterator:
            text, result, error = future.result()

            if error:
                print(f"Error processing '{text}': {error}", file=sys.stderr)
                continue

            results.append(result)

            # Write to file immediately if specified (thread-safe)
            if output_file and result:
                write_result_to_file(output_file, result, output_format, file_lock)

    print(f"Completed processing {len(results)} phenotypes")
    return results


def process_batch_file(
    matcher: PhenotypeMatcher, input_file: str, context: Optional[str], split_by: str
) -> list:
    """Process batch file with one phenotype per line (legacy function)."""
    results = []

    with open(input_file, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            try:
                input_data = PhenotypeInput(
                    text=line,
                    context=context,
                    split_by=split_by if split_by != "none" else None,
                )
                output = matcher.match(input_data)
                results.append(output)
            except Exception as e:
                print(f"Error processing line {line_num}: {e}", file=sys.stderr)

    return results


def format_output(results: list, format_type: str) -> str:
    """Format results according to specified format."""
    if format_type == "json":
        # Full JSON output
        output_dicts = [r.to_dict() for r in results]
        return json.dumps(output_dicts, indent=2)

    elif format_type == "hpo_ids":
        # Just HPO IDs (present only)
        all_ids = []
        for result in results:
            all_ids.extend(result.get_hpo_ids(excluded=False))
        return "\n".join(all_ids)

    elif format_type == "tsv":
        # TSV format
        lines = [
            "input\thpo_ids\thpo_labels\texcluded_ids\tseverity_ids\tomim_ids\tmondo_ids"
        ]
        for result in results:
            hpo_ids = "|".join(result.get_hpo_ids(excluded=False))
            hpo_labels = "|".join([p.label for p in result.get_present_phenotypes()])
            excluded_ids = "|".join(result.get_hpo_ids(excluded=True))
            severity_ids = "|".join(
                [p.severity_id for p in result.phenotypes if p.severity_id]
            )
            omim_ids = "|".join(result.get_omim_ids())
            mondo_ids = "|".join(result.get_mondo_ids())

            lines.append(
                f"{result.raw_input}\t{hpo_ids}\t{hpo_labels}\t{excluded_ids}\t{severity_ids}\t{omim_ids}\t{mondo_ids}"
            )
        return "\n".join(lines)

    elif format_type == "summary":
        # Human-readable summary
        lines = []
        for i, result in enumerate(results, 1):
            if len(results) > 1:
                lines.append(f"\n=== Result {i} ===")
            lines.append(f"Input: {result.raw_input}")

            # Show NER extracted terms if available
            if hasattr(result, "ner_extracted_terms") and result.ner_extracted_terms:
                lines.append(
                    f"\nNER Extracted ({len(result.ner_extracted_terms)} terms):"
                )
                for ner_term in result.ner_extracted_terms:
                    term_str = ner_term.get("term", "")
                    modifiers = ner_term.get("modifiers", [])
                    excluded = ner_term.get("excluded", False)
                    mod_str = f" [{', '.join(modifiers)}]" if modifiers else ""
                    excl_str = " (excluded)" if excluded else ""
                    lines.append(f"  • '{term_str}'{mod_str}{excl_str}")

            lines.append(
                f"\nPresent Phenotypes ({len(result.get_present_phenotypes())}):"
            )
            for pheno in result.get_present_phenotypes():
                severity_str = (
                    f" [{pheno.severity_label}]" if pheno.severity_label else ""
                )
                lines.append(f"  - {pheno.hpo_id}: {pheno.label}{severity_str}")

            if result.get_excluded_phenotypes():
                lines.append(
                    f"\nExcluded Phenotypes ({len(result.get_excluded_phenotypes())}):"
                )
                for pheno in result.get_excluded_phenotypes():
                    lines.append(f"  - {pheno.hpo_id}: {pheno.label}")

            if result.diseases:
                lines.append(f"\nDiseases ({len(result.diseases)}):")
                for disease in result.diseases:
                    if disease.mondo_id:
                        lines.append(f"  - {disease.mondo_id}: {disease.mondo_label}")
                    for omim_label in disease.omim_labels:
                        lines.append(f"  - OMIM: {omim_label}")

            lines.append(
                f"\nProcessing: {result.processing_metadata.get('processing_time_seconds', 0):.2f}s"
            )

        return "\n".join(lines)

    else:
        return json.dumps([r.to_dict() for r in results], indent=2)


def run_test_cases(specific_case_id: Optional[str], config: MatcherConfig):
    """Run difficult test cases to validate the matcher."""
    print("Phenotype Matcher - Test Suite")
    print("=" * 80)

    # Initialize matcher
    print("\nInitializing matcher...")
    try:
        matcher = PhenotypeMatcher(config)
    except Exception as e:
        print(f"Error initializing matcher: {e}", file=sys.stderr)
        sys.exit(1)

    # Get test cases
    if specific_case_id:
        case = test_cases.get_test_case(specific_case_id)
        if not case:
            print(f"Error: Test case '{specific_case_id}' not found", file=sys.stderr)
            print(f"\nAvailable test cases:")
            for c in test_cases.get_all_test_cases():
                print(f"  - {c['id']}: {c['name']}")
            sys.exit(1)
        cases_to_run = [case]
    else:
        cases_to_run = test_cases.get_all_test_cases()

    # Run tests
    total_cases = len(cases_to_run)
    passed = 0
    failed = 0

    for i, case in enumerate(cases_to_run, 1):
        print(f"\n[{i}/{total_cases}] Running test case: {case['name']}")
        print(f"Description: {case['description']}")

        # Run the matcher
        input_data = PhenotypeInput(
            text=case["input"],
            context=case.get("context"),
            gene_hint=case.get("gene_hint"),
        )
        try:
            output = matcher.match(input_data)
        except Exception as e:
            print(f"ERROR: Matcher failed with exception: {e}")
            failed += 1
            continue

        # Validate results
        results = test_cases.validate_results(case["id"], output)
        test_cases.print_validation_results(results, output)

        if results["pass"]:
            passed += 1
        else:
            failed += 1

    # Summary
    print(f"\n{'=' * 80}")
    print(f"Test Summary")
    print(f"{'=' * 80}")
    print(f"Total cases: {total_cases}")
    print(
        f"Passed: {passed} ({100 * passed // total_cases if total_cases > 0 else 0}%)"
    )
    print(f"Failed: {failed}")
    print(f"{'=' * 80}")

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
