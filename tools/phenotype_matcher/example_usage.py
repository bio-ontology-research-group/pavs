#!/usr/bin/env python3
"""
Example usage of the Phenotype Matcher tool.

This script demonstrates various ways to use the phenotype matcher.
Run with: uv run python example_usage.py
"""

import os
from phenotype_matcher import PhenotypeMatcher, PhenotypeInput, MatcherConfig


def example_basic_usage():
    """Example 1: Basic usage with default configuration."""
    print("=" * 80)
    print("Example 1: Basic Usage")
    print("=" * 80)

    # Initialize matcher with default config
    matcher = PhenotypeMatcher()

    # Match a simple phenotype
    input_data = PhenotypeInput(text="seizure")
    output = matcher.match(input_data)

    print(f"\nInput: {output.raw_input}")
    print(f"\nMatched {len(output.phenotypes)} phenotype(s):")
    for pheno in output.phenotypes:
        print(f"  - {pheno.hpo_id}: {pheno.label}")
        if pheno.excluded:
            print(f"    (Excluded/Negated)")
        if pheno.severity_label:
            print(f"    Severity: {pheno.severity_label} ({pheno.severity_id})")

    print(
        f"\nProcessing time: {output.processing_metadata.get('processing_time_seconds', 0):.2f}s"
    )


def example_severity_detection():
    """Example 2: Phenotype with severity modifier."""
    print("\n" + "=" * 80)
    print("Example 2: Severity Detection")
    print("=" * 80)

    matcher = PhenotypeMatcher()

    input_data = PhenotypeInput(text="severe intellectual disability")
    output = matcher.match(input_data)

    print(f"\nInput: {output.raw_input}")
    print(f"\nMatched {len(output.phenotypes)} phenotype(s):")
    for pheno in output.phenotypes:
        print(f"  - {pheno.hpo_id}: {pheno.label}")
        if pheno.severity_label:
            print(f"    Severity: {pheno.severity_label} ({pheno.severity_id})")


def example_multiple_phenotypes():
    """Example 3: Multiple phenotypes in one description."""
    print("\n" + "=" * 80)
    print("Example 3: Multiple Phenotypes")
    print("=" * 80)

    matcher = PhenotypeMatcher()

    input_data = PhenotypeInput(text="seizures, hypotonia, and feeding difficulties")
    output = matcher.match(input_data)

    print(f"\nInput: {output.raw_input}")
    print(f"\nMatched {len(output.phenotypes)} phenotype(s):")
    for i, pheno in enumerate(output.phenotypes, 1):
        print(f"  {i}. {pheno.hpo_id}: {pheno.label}")


def example_negation():
    """Example 4: Negation detection."""
    print("\n" + "=" * 80)
    print("Example 4: Negation Detection")
    print("=" * 80)

    matcher = PhenotypeMatcher()

    input_data = PhenotypeInput(text="no cardiac abnormalities, normal vision")
    output = matcher.match(input_data)

    print(f"\nInput: {output.raw_input}")

    present = output.get_present_phenotypes()
    excluded = output.get_excluded_phenotypes()

    print(f"\nPresent phenotypes: {len(present)}")
    for pheno in present:
        print(f"  - {pheno.hpo_id}: {pheno.label}")

    print(f"\nExcluded/Negated phenotypes: {len(excluded)}")
    for pheno in excluded:
        print(f"  - {pheno.hpo_id}: {pheno.label}")


def example_custom_config():
    """Example 5: Custom configuration."""
    print("\n" + "=" * 80)
    print("Example 5: Custom Configuration")
    print("=" * 80)

    # Custom configuration
    config = MatcherConfig(
        embedding_model="fast",  # Use faster model
        llm_model="fast",  # Use faster LLM
        top_k_phenotype=3,  # Only retrieve 3 candidates
        device="cpu",  # Force CPU
        debug=False,  # Disable debug logging
    )

    matcher = PhenotypeMatcher(config)

    input_data = PhenotypeInput(text="global developmental delay")
    output = matcher.match(input_data)

    print(f"\nConfiguration:")
    print(f"  Embedding model: {config.embedding_model}")
    print(f"  LLM model: {config.llm_model}")
    print(f"  Top-K: {config.top_k_phenotype}")

    print(f"\nInput: {output.raw_input}")
    print(f"\nMatched {len(output.phenotypes)} phenotype(s):")
    for pheno in output.phenotypes:
        print(f"  - {pheno.hpo_id}: {pheno.label}")


def example_output_formats():
    """Example 6: Different output formats."""
    print("\n" + "=" * 80)
    print("Example 6: Output Formats")
    print("=" * 80)

    matcher = PhenotypeMatcher()

    input_data = PhenotypeInput(text="seizures and hypotonia")
    output = matcher.match(input_data)

    print(f"\nInput: {output.raw_input}")

    # Get HPO IDs only
    hpo_ids = output.get_hpo_ids()
    print(f"\nHPO IDs (present): {hpo_ids}")

    # Get excluded IDs
    excluded_ids = output.get_hpo_ids(excluded=True)
    print(f"HPO IDs (excluded): {excluded_ids}")

    # Get OMIM IDs
    omim_ids = output.get_omim_ids()
    print(f"OMIM IDs: {omim_ids}")

    # Get MONDO IDs
    mondo_ids = output.get_mondo_ids()
    print(f"MONDO IDs: {mondo_ids}")

    # Get as dictionary
    output_dict = output.to_dict()
    print(f"\nDictionary format has keys: {list(output_dict.keys())}")


def example_context_disambiguation():
    """Example 7: Using context for disambiguation."""
    print("\n" + "=" * 80)
    print("Example 7: Context-Based Disambiguation")
    print("=" * 80)

    matcher = PhenotypeMatcher()

    # ASD can mean "Autism Spectrum Disorder" or "Atrial Septal Defect"
    # Context helps disambiguate

    print("\nCase 1: ASD with cardiac context")
    input_data = PhenotypeInput(
        text="ASD", context="cardiac defect, heart murmur, congenital heart disease"
    )
    output = matcher.match(input_data)

    print(f"  Input: {output.raw_input}")
    print(f"  Context: {input_data.context}")
    for pheno in output.phenotypes:
        print(f"  Matched: {pheno.hpo_id}: {pheno.label}")

    print("\nCase 2: ASD with neurodevelopmental context")
    input_data = PhenotypeInput(
        text="ASD", context="developmental delay, speech delay, social difficulties"
    )
    output = matcher.match(input_data)

    print(f"  Input: {output.raw_input}")
    print(f"  Context: {input_data.context}")
    for pheno in output.phenotypes:
        print(f"  Matched: {pheno.hpo_id}: {pheno.label}")


def main():
    """Run all examples."""
    # Check for API key
    if not os.getenv("OPENROUTER_API_KEY"):
        print("WARNING: OPENROUTER_API_KEY not set!")
        print("Some examples may fail without an API key.")
        print("Set it with: export OPENROUTER_API_KEY='your-key-here'")
        print()

    try:
        # Run examples
        example_basic_usage()
        example_severity_detection()
        example_multiple_phenotypes()
        example_negation()
        example_custom_config()
        example_output_formats()
        example_context_disambiguation()

        print("\n" + "=" * 80)
        print("All examples completed successfully!")
        print("=" * 80)

    except Exception as e:
        print(f"\n\nError running examples: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
