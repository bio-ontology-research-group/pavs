#!/usr/bin/env python3
"""
GraphRAG-based phenotype normalization pipeline.

This script uses semantic embeddings + graph structure + LLM reasoning
for high-accuracy phenotype and disease normalization.

Usage:
    python scripts/process_all_phenotypes_graphrag.py --help
"""

import pandas as pd
import re
import os
import json
import argparse
import sys
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from normalization.graph_rag_normalization import GraphRAG

OUTPUT_DIR = "data/processed_graphrag"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Recommended sentence transformer models
EMBEDDING_MODELS = {
    "fast": "all-MiniLM-L6-v2",  # Fast, 384 dim, 80MB
    "balanced": "all-mpnet-base-v2",  # Balanced, 768 dim, 420MB (recommended)
    "accurate": "BAAI/bge-large-en-v1.5",  # High accuracy, 1024 dim, 1.3GB
    "medical": "pritamdeka/S-PubMedBert-MS-MARCO",  # Medical domain, 768 dim
    "biobert": "dmis-lab/biobert-base-cased-v1.2",  # Biomedical, 768 dim
    "sapbert": "cambridgeltl/SapBERT-from-PubMedBERT-fulltext",  # Biomedical ontology alignment, 768 dim
}

# Recommended LLM models
LLM_MODELS = {
    "fast": "openai/gpt-oss-120b",
    "accurate": "anthropic/claude-3.5-sonnet",
    "balanced": "google/gemini-2.0-flash-exp:free",
    "cheap": "deepseek/deepseek-chat",
}

COLUMNS = [
    "id",
    "source_file",
    "sex_id",
    "sex_label",
    "age",
    "ancestry_id",
    "ancestry_label",
    "country_id",
    "country_label",
    "consanguinity_id",
    "consanguinity_label",
    "family_history_id",
    "family_history_label",
    "hpo_ids",
    "hpo_labels",
    "hpo_excluded_ids",
    "hpo_excluded_labels",
    "hpo_severity_ids",
    "hpo_severity_labels",
    "disease_omim_gene_ids",
    "disease_omim_phenotype_ids",
    "disease_omim_labels",
    "disease_omim_excluded_ids",
    "disease_mondo_ids",
    "disease_mondo_labels",
    "disease_mondo_excluded_ids",
    "gene_hgnc_ids",
    "gene_symbols",
    "variant_hgvs",
    "variant_validated",
    "variant_classification",
    "variant_evidence",
    "publication_id",
    "genotype_id",
    "genotype_label",
    "solved",
]


def process_file_graphrag(filename, graph_rag, limit=None, top_k=5):
    """Process a single file using GraphRAG normalization."""
    path = f"data/phenotypes/{filename}"
    if not os.path.exists(path):
        print(f"File not found: {path}")
        return

    print(f"\n🔬 Processing {filename} with GraphRAG...")
    print(f"  Settings: top_k={top_k}, model={graph_rag.model_name}")

    # Load file (simplified - handle only common cases for now)
    df = pd.read_csv(path, sep="\t", comment="#").rename(
        columns=lambda x: str(x).strip()
    )

    if limit:
        df = df.head(limit)

    # Check for resume
    output_path = os.path.join(OUTPUT_DIR, filename.replace(".tsv", "_normalized.tsv"))
    processed_ids = set()

    if os.path.exists(output_path):
        try:
            existing_df = pd.read_csv(output_path, sep="\t")
            processed_ids = set(existing_df["id"].astype(str))
            print(f"  📂 Resuming: {len(processed_ids)} rows already processed")
        except Exception as e:
            print(f"  ⚠️  Could not read existing file: {e}")

    file_mode = "a" if os.path.exists(output_path) else "w"
    write_header = not os.path.exists(output_path)
    rows_processed = 0

    for idx, row in tqdm(df.iterrows(), total=len(df), desc=filename):
        rec = {k: "" for k in COLUMNS}
        rec["source_file"] = filename

        # Get ID
        raw_id = row.get("Case", row.get("ID", row.get("ID_p", idx)))
        try:
            f_id = float(raw_id)
            if f_id.is_integer():
                rec["id"] = str(int(f_id))
            else:
                rec["id"] = str(raw_id)
        except:
            rec["id"] = str(raw_id)

        # Skip if already processed
        if rec["id"] in processed_ids:
            continue

        # Extract phenotype text
        p_val = str(
            row.get(
                "Phenotype",
                row.get(
                    "HPOs",
                    row.get(
                        "Additional clinical phenotype",
                        row.get("phenotypes", row.get("hpo_ids", "")),
                    ),
                ),
            )
        )

        # File-specific splitting logic
        raw_terms = []
        if p_val and p_val != "nan":
            if filename == "marwa-variants.tsv":
                # Uses | separator, already has HPO IDs - extract just the text before (HP:...)
                terms_with_ids = re.split(r"\|", p_val)
                raw_terms = [
                    re.sub(r"\(HP:\d+\)", "", t).strip()
                    for t in terms_with_ids
                    if t.strip()
                ]
            elif filename == "ahmed-variants.tsv":
                # Has embedded HPO IDs - keep complex terms intact
                # Don't split by "and"/"with" - let LLM extract multiple phenotypes from complex terms
                terms_with_ids = re.split(r"[,|]\s*", p_val)
                raw_terms = [
                    re.sub(r"\(HP:\d+\)", "", t).strip()
                    for t in terms_with_ids
                    if t.strip()
                ]
            elif filename == "fawzan-variants.tsv":
                # Uses commas and / for alternatives
                raw_terms = [
                    t.strip()
                    for t in re.split(r"[,;]\s*", p_val)
                    if t.strip() and t.lower() != "nan"
                ]
                # Further split by / for alternatives like "intolerance/fatigue"
                expanded = []
                for t in raw_terms:
                    if "/" in t and not any(x in t.lower() for x in ["w/", "c/", "s/"]):
                        expanded.extend([x.strip() for x in t.split("/")])
                    else:
                        expanded.append(t)
                raw_terms = expanded
            else:
                # Default: split by comma, semicolon, pipe
                raw_terms = [
                    t.strip()
                    for t in re.split(r"[,|;]\s*", p_val)
                    if t.strip() and t.lower() != "nan"
                ]

            # Process each term with GraphRAG
            hpo_ids = []
            hpo_labels = []
            hpo_excluded_ids = []
            hpo_excluded_labels = []
            hpo_severity_ids = []
            hpo_severity_labels = []
            omim_gene_ids = set()
            omim_pheno_ids = set()
            omim_labels = set()
            mondo_ids = set()
            mondo_labels = set()

            for term in raw_terms:
                if not term or len(term) < 3:
                    continue

                # Use GraphRAG to normalize
                result = graph_rag.normalize_term(term, context_text="", k=top_k)

                if result:
                    # Handle case where LLM returns list or dict
                    if isinstance(result, list):
                        # LLM returned array - skip for now
                        print(
                            f"Warning: LLM returned list for term '{term}', expected dict"
                        )
                        continue

                    if not isinstance(result, dict):
                        print(
                            f"Warning: Unexpected result type {type(result)} for term '{term}'"
                        )
                        continue

                    # Get labels from LLM (NO IDs)
                    phenotype_labels = result.get("phenotype_labels", [])
                    is_excluded = result.get("excluded", False)
                    severity_label = result.get("severity_label")

                    # Handle single label as list
                    if isinstance(phenotype_labels, str):
                        phenotype_labels = (
                            [phenotype_labels] if phenotype_labels != "null" else []
                        )

                    # Look up HPO IDs from labels using ontology
                    for label in phenotype_labels:
                        if not label or label in ["None", "null"]:
                            continue

                        # Look up ID from ontology
                        from pyhpo import Ontology

                        Ontology()

                        # Find matching term
                        matches = [
                            t for t in Ontology if t.name.lower() == label.lower()
                        ]
                        if matches:
                            pheno_id = matches[0].id
                            pheno_label = matches[0].name  # Use ontology label, not LLM

                            if is_excluded:
                                hpo_excluded_ids.append(pheno_id)
                                hpo_excluded_labels.append(pheno_label)
                            else:
                                hpo_ids.append(pheno_id)
                                hpo_labels.append(pheno_label)
                        else:
                            print(
                                f"Warning: Could not find HPO term for label '{label}'"
                            )

                    # Handle severity - look up ID from label
                    if severity_label and severity_label not in ["None", "null", None]:
                        matches = [
                            t
                            for t in Ontology
                            if t.name.lower() == severity_label.lower()
                        ]
                        if matches:
                            sev_id = matches[0].id
                            sev_label = matches[0].name  # From ontology
                            # Validate it's actually a severity term
                            if sev_id.startswith("HP:"):
                                # TODO: Add validation that it's under HP:0012824
                                if hpo_ids:  # Pair with last phenotype
                                    hpo_severity_ids.append(sev_id)
                                    hpo_severity_labels.append(sev_label)

                    # Handle diseases (OMIM/MONDO) from LLM response
                    disease_labels = result.get("disease_labels", [])
                    if isinstance(disease_labels, list):
                        for disease in disease_labels:
                            if isinstance(disease, dict):
                                omim_label = disease.get("omim")
                                mondo_label = disease.get("mondo")

                                # Look up IDs from labels
                                # Note: Would need utils.lookup_omim_id_by_label() and lookup_mondo_id_by_label()
                                # For now, just store labels
                                if omim_label and omim_label not in ["None", "null"]:
                                    omim_labels.add(omim_label)
                                if mondo_label and mondo_label not in ["None", "null"]:
                                    mondo_labels.add(mondo_label)

            rec["hpo_ids"] = "|".join(hpo_ids)
            rec["hpo_labels"] = "|".join(hpo_labels)
            rec["hpo_excluded_ids"] = "|".join(hpo_excluded_ids)
            rec["hpo_excluded_labels"] = "|".join(hpo_excluded_labels)
            rec["hpo_severity_ids"] = "|".join(hpo_severity_ids)
            rec["hpo_severity_labels"] = "|".join(hpo_severity_labels)
            rec["disease_omim_labels"] = "|".join(omim_labels)
            rec["disease_mondo_labels"] = "|".join(mondo_labels)

        # Write incrementally
        row_df = pd.DataFrame([rec])
        row_df.to_csv(
            output_path,
            sep="\t",
            mode=file_mode,
            header=write_header,
            index=False,
        )

        if write_header:
            write_header = False
            file_mode = "a"

        rows_processed += 1

    print(f"  ✅ Processed {rows_processed} new rows")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="GraphRAG-based phenotype normalization with semantic embeddings + graph structure.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all files with default settings
  python scripts/process_all_phenotypes_graphrag.py
  
  # Process specific file with custom settings
  python scripts/process_all_phenotypes_graphrag.py --files ahmed-variants.tsv --top-k 10 --embedding-model accurate
  
  # Use medical-domain embeddings with high k
  python scripts/process_all_phenotypes_graphrag.py --embedding-model medical --top-k 15 --llm accurate
  
  # Process in parallel with overwrite
  python scripts/process_all_phenotypes_graphrag.py --parallel --overwrite --limit 100
        """,
    )

    # File selection
    parser.add_argument(
        "--files",
        nargs="+",
        help="Specific files to process (default: all)",
        default=None,
    )
    parser.add_argument(
        "limit",
        type=int,
        nargs="?",
        default=None,
        help="Limit number of records per file (for testing)",
    )

    # GraphRAG parameters
    parser.add_argument(
        "--embedding-model",
        choices=list(EMBEDDING_MODELS.keys()) + list(EMBEDDING_MODELS.values()),
        default="balanced",
        help="Sentence transformer model for embeddings (default: balanced = all-mpnet-base-v2)",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=5,
        help="Number of top candidates to retrieve (default: 5)",
    )
    parser.add_argument(
        "--llm",
        choices=list(LLM_MODELS.keys()) + list(LLM_MODELS.values()),
        default="fast",
        help="LLM model for final selection (default: fast = openai/gpt-oss-120b)",
    )

    # Processing options
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete existing normalized files and reprocess from scratch",
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Process files in parallel (up to 7 files simultaneously)",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Show debug output: query terms, retrieved candidates, and LLM responses",
    )

    args = parser.parse_args()

    # Resolve model names
    embedding_model = EMBEDDING_MODELS.get(args.embedding_model, args.embedding_model)
    llm_model = LLM_MODELS.get(args.llm, args.llm)

    print("=" * 80)
    print("GraphRAG Phenotype Normalization Pipeline")
    print("=" * 80)
    print(f"Embedding Model: {embedding_model}")
    print(f"LLM Model: {llm_model}")
    print(f"Top-K Candidates: {args.top_k}")
    print("=" * 80)

    # Initialize GraphRAG
    print("\n🔧 Initializing GraphRAG...")
    graph_rag = GraphRAG(
        model_name=embedding_model, llm_model=llm_model, debug=args.debug
    )

    # File list
    all_files = [
        "ahmed-variants.tsv",
        "ahmed-pmid28454995.tsv",
        "fawzan-variants.tsv",
        "marwa-variants.tsv",
        "PMC6562004.tsv",
        "PMC7082194.tsv",
        "ddd-diagnoses.tsv",
    ]

    files_to_process = args.files if args.files else all_files
    files_to_process = [f for f in files_to_process if f in all_files]

    if not files_to_process:
        print("❌ No valid files specified. Available files:")
        for f in all_files:
            print(f"  - {f}")
        sys.exit(1)

    # Handle overwrite
    if args.overwrite:
        for f in files_to_process:
            output_file = os.path.join(OUTPUT_DIR, f.replace(".tsv", "_normalized.tsv"))
            if os.path.exists(output_file):
                print(f"🗑️  Deleting existing: {output_file}")
                os.remove(output_file)

    # Check which files need processing
    files_needing_work = []
    for f in files_to_process:
        output_file = os.path.join(OUTPUT_DIR, f.replace(".tsv", "_normalized.tsv"))
        if not os.path.exists(output_file):
            files_needing_work.append(f)
        else:
            print(f"✓ Skipping {f} (already processed, use --overwrite to reprocess)")

    if not files_needing_work:
        print("\n✅ All files already processed! Use --overwrite to reprocess.")
        sys.exit(0)

    print(f"\n📋 Processing {len(files_needing_work)} file(s)...")

    # Process files
    if args.parallel and len(files_needing_work) > 1:
        print("🚀 Parallel mode: Processing files concurrently")
        with ThreadPoolExecutor(
            max_workers=min(len(files_needing_work), 7)
        ) as executor:
            futures = {
                executor.submit(
                    process_file_graphrag, f, graph_rag, args.limit, args.top_k
                ): f
                for f in files_needing_work
            }
            for future in as_completed(futures):
                filename = futures[future]
                try:
                    future.result()
                    print(f"✅ Completed: {filename}")
                except Exception as e:
                    print(f"❌ Error processing {filename}: {e}")
    else:
        for f in files_needing_work:
            process_file_graphrag(f, graph_rag, args.limit, args.top_k)
            print(f"✅ Completed: {f}")

    print("\n" + "=" * 80)
    print("✅ GraphRAG Processing Complete!")
    print(f"📁 Output directory: {OUTPUT_DIR}")
    print("=" * 80)
