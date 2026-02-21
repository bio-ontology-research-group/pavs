#!/usr/bin/env python3
"""
prepare_all.py — Full PAVS data preparation pipeline.

Steps (each skipped if output already exists unless --force):
  1. generate_phenopackets_v2  → phenopackets/generated_v2/ + data/PAVS_phenopackets.json
  2. compute_hpo_ic            → rdf_output/hpo_ic.ttl
  3. generate_rdf              → rdf_output/*.ttl
  4. load_virtuoso             → Virtuoso triple store (waits for endpoint)

Usage:
    uv run python intake/prepare_all.py [--force] [--endpoint http://localhost:8890]
"""

import argparse
import subprocess
import sys
import time
import urllib.request
from pathlib import Path


def log(msg: str):
    print(f"\n{'='*60}\n  {msg}\n{'='*60}", flush=True)


def run(cmd: list[str]):
    print(f"  $ {' '.join(cmd)}", flush=True)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"  ERROR: command failed (exit {result.returncode})", flush=True)
        sys.exit(result.returncode)


def wait_for_virtuoso(endpoint: str, max_seconds: int = 120):
    print(f"  Waiting for Virtuoso at {endpoint} …", flush=True)
    deadline = time.time() + max_seconds
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(f"{endpoint}/sparql", timeout=5):
                print("  Virtuoso is up!", flush=True)
                return
        except Exception:
            time.sleep(2)
    print(f"  ERROR: Virtuoso not available after {max_seconds}s", flush=True)
    sys.exit(1)


def needs_run(paths: list[str]) -> bool:
    """Return True if any of the expected output paths is missing."""
    return any(not Path(p).exists() for p in paths)


def main():
    parser = argparse.ArgumentParser(description="Full PAVS data prep + ingestion pipeline")
    parser.add_argument("--force", action="store_true", help="Re-run all steps even if outputs exist")
    parser.add_argument("--endpoint", default="http://localhost:8890", help="Virtuoso base URL")
    parser.add_argument("--skip-load", action="store_true", help="Skip Virtuoso loading step")
    # Pass-through args for individual scripts
    parser.add_argument("--input", default="data/combined_annotated.tsv")
    parser.add_argument("--hpo", default="ontology/hp.obo")
    parser.add_argument("--hpoa", default="data/reference/phenotype.hpoa")
    parser.add_argument("--literature-dir", default="phenopackets/0.1.26")
    args = parser.parse_args()

    python = [sys.executable]

    # ------------------------------------------------------------------
    # Step 1: Generate phenopackets v2
    # ------------------------------------------------------------------
    phenopackets_json = "data/PAVS_phenopackets.json"
    phenopackets_zip  = "data/PAVS_phenopackets.zip"

    if args.force or needs_run([phenopackets_json, phenopackets_zip]):
        log("Step 1/4 — Generating Phenopackets v2 from combined_annotated.tsv")
        run(python + [
            "intake/generate_phenopackets_v2.py",
            "--input",       args.input,
            "--output-dir",  "phenopackets/generated_v2",
            "--combined",    phenopackets_json,
            "--zip",         phenopackets_zip,
        ])
    else:
        log("Step 1/4 — Phenopackets already exist, skipping (use --force to regenerate)")

    # ------------------------------------------------------------------
    # Step 2: Compute HPO IC values
    # ------------------------------------------------------------------
    hpo_ic_ttl = "rdf_output/hpo_ic.ttl"

    if args.force or needs_run([hpo_ic_ttl]):
        log("Step 2/4 — Computing HPO information content")
        Path("rdf_output").mkdir(exist_ok=True)
        run(python + [
            "intake/compute_hpo_ic.py",
            "--hpo",    args.hpo,
            "--hpoa",   args.hpoa,
            "--output", hpo_ic_ttl,
        ])
    else:
        log("Step 2/4 — hpo_ic.ttl already exists, skipping")

    # ------------------------------------------------------------------
    # Step 3: Generate RDF (cases, genes, HPOA, literature)
    # ------------------------------------------------------------------
    rdf_outputs = [
        "rdf_output/cases.ttl",
        "rdf_output/genes.ttl",
        "rdf_output/hpoa.ttl",
        "rdf_output/literature.ttl",
    ]

    if args.force or needs_run(rdf_outputs):
        log("Step 3/4 — Generating RDF Turtle files")
        run(python + [
            "intake/generate_rdf.py",
            "--phenopackets",   phenopackets_json,
            "--annotated",      args.input,
            "--hpoa",           args.hpoa,
            "--literature-dir", args.literature_dir,
            "--output-dir",     "rdf_output/",
        ])
    else:
        log("Step 3/4 — RDF files already exist, skipping")

    # ------------------------------------------------------------------
    # Step 4: Load into Virtuoso
    # ------------------------------------------------------------------
    if args.skip_load:
        log("Step 4/4 — Virtuoso load skipped (--skip-load)")
    else:
        log("Step 4/4 — Loading RDF into Virtuoso")
        wait_for_virtuoso(args.endpoint)
        run(python + [
            "intake/load_virtuoso.py",
            "--endpoint",              args.endpoint,
            "--input-dir",             "rdf_output/",
            "--virtuoso-import-dir",   "/rdf_import",
        ])

    log("All steps complete.")


if __name__ == "__main__":
    main()
