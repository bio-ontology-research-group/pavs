#!/usr/bin/env python3
"""
validate_shex.py — Validate generated PAVS RDF against ShEx shapes.

Loads all Turtle files from rdf_output/ and validates instances of each
PAVS type against the corresponding ShEx shape.

Usage:
    uv run python intake/validate_shex.py \
        --shapes ontology/pavs_shapes.shex \
        --rdf-dir rdf_output/
"""

import argparse
import sys
from pathlib import Path

from pyshex import ShExEvaluator
from rdflib import RDF, Graph, Namespace

PAVS = Namespace("http://pavs.kaust.edu.sa/ontology/")

# Map each PAVS class to its ShEx shape name
SHAPE_MAP = {
    PAVS.Case: "http://pavs.kaust.edu.sa/ontology/CaseShape",
    PAVS.GenomicVariant: "http://pavs.kaust.edu.sa/ontology/GenomicVariantShape",
    PAVS.GeneRecord: "http://pavs.kaust.edu.sa/ontology/GeneRecordShape",
    PAVS.DiseaseHPOAssociation: "http://pavs.kaust.edu.sa/ontology/DiseaseHPOAssociationShape",
}


def load_rdf(rdf_dir: str) -> Graph:
    """Load all .ttl files from rdf_dir into a single graph."""
    g = Graph()
    rdf_path = Path(rdf_dir)
    ttl_files = sorted(rdf_path.glob("*.ttl"))
    if not ttl_files:
        print(f"ERROR: No .ttl files found in {rdf_dir}")
        sys.exit(1)

    for ttl in ttl_files:
        # Skip hpo_ic.ttl — it contains HPO hierarchy, not PAVS instances
        if ttl.name == "hpo_ic.ttl":
            continue
        print(f"  Loading {ttl.name} …")
        g.parse(str(ttl), format="turtle")

    print(f"  Total triples: {len(g)}")
    return g


def validate(shapes_file: str, rdf_dir: str, max_errors: int = 50) -> bool:
    """Validate all PAVS instances against ShEx shapes."""
    print(f"Loading ShEx shapes from {shapes_file}")
    shapes_text = Path(shapes_file).read_text(encoding="utf-8")

    print(f"Loading RDF from {rdf_dir}")
    g = load_rdf(rdf_dir)

    all_passed = True
    total_validated = 0
    total_failed = 0

    for pavs_class, shape_uri in SHAPE_MAP.items():
        class_label = str(pavs_class).split("/")[-1]
        instances = list(g.subjects(RDF.type, pavs_class))

        if not instances:
            print(f"\n  {class_label}: no instances found — skipping")
            continue

        print(f"\n  Validating {len(instances)} {class_label} instances against {shape_uri.split('/')[-1]} …")

        errors = 0
        for i, instance in enumerate(instances):
            evaluator = ShExEvaluator(
                rdf=g,
                schema=shapes_text,
                focus=str(instance),
                start=shape_uri,
            )
            results = evaluator.evaluate()

            for result in results:
                if not result.result:
                    errors += 1
                    total_failed += 1
                    if errors <= max_errors:
                        print(f"    FAIL: {instance}")
                        # Truncate long reason strings
                        reason = result.reason or "unknown"
                        if len(reason) > 200:
                            reason = reason[:200] + "…"
                        print(f"          {reason}")

            total_validated += 1
            if (i + 1) % 500 == 0:
                print(f"    … {i + 1}/{len(instances)} validated ({errors} failures so far)")

        if errors == 0:
            print(f"    PASS: all {len(instances)} instances valid")
        else:
            all_passed = False
            if errors > max_errors:
                print(f"    … and {errors - max_errors} more failures (truncated)")
            print(f"    FAIL: {errors}/{len(instances)} instances invalid")

    print(f"\n{'=' * 60}")
    print(f"Total validated: {total_validated}")
    print(f"Total failures:  {total_failed}")
    print(f"Result:          {'PASS' if all_passed else 'FAIL'}")
    print(f"{'=' * 60}")

    return all_passed


def main():
    parser = argparse.ArgumentParser(description="Validate PAVS RDF against ShEx shapes")
    parser.add_argument("--shapes", default="ontology/pavs_shapes.shex",
                        help="Path to ShEx shapes file")
    parser.add_argument("--rdf-dir", default="rdf_output/",
                        help="Directory containing .ttl files")
    parser.add_argument("--max-errors", type=int, default=50,
                        help="Max errors to display per shape (default: 50)")
    args = parser.parse_args()

    ok = validate(args.shapes, args.rdf_dir, args.max_errors)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
