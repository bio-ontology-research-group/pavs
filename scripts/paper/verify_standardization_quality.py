#!/usr/bin/env python3
"""
verify_standardization_quality.py — Verify HPO and HGVS standardization quality
for the PAVS paper.

Checks:
1. HPO ID validity: all assigned HPO IDs exist in hp.obo
2. HGVS format validity: variant strings match expected HGVS patterns
3. Gene symbol validity: gene symbols exist in NCBI gene info
4. Disease ID validity: OMIM/MONDO IDs match expected patterns

Usage:
    uv run python scripts/paper/verify_standardization_quality.py
"""

import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent


def load_hpo_ids(obo_path):
    """Load all valid HP IDs from hp.obo."""
    ids = set()
    with open(obo_path) as f:
        for line in f:
            if line.startswith("id: HP:"):
                ids.add(line.strip().replace("id: ", ""))
    return ids


def load_gene_symbols(gene_info_path):
    """Load valid gene symbols from NCBI gene info."""
    symbols = set()
    with open(gene_info_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) > 2:
                symbols.add(parts[2])  # Symbol column
    return symbols


def main():
    tsv_path = REPO / "data" / "combined_normalized.tsv"
    obo_path = REPO / "ontology" / "hp.obo"
    gene_info_path = REPO / "data" / "Homo_sapiens.gene_info"

    if not tsv_path.exists():
        print(f"ERROR: {tsv_path} not found")
        sys.exit(1)

    df = pd.read_csv(tsv_path, sep="\t", dtype=str, low_memory=False)
    print(f"Loaded {len(df)} rows from {tsv_path.name}")
    print(f"Columns: {len(df.columns)}")
    print()

    results = {}

    # ---------------------------------------------------------------
    # 1. HPO ID validity
    # ---------------------------------------------------------------
    print("=" * 60)
    print("1. HPO ID VALIDITY")
    print("=" * 60)

    if obo_path.exists():
        valid_hpo = load_hpo_ids(obo_path)
        print(f"Loaded {len(valid_hpo)} valid HPO terms from hp.obo")

        total_hpo = 0
        valid_count = 0
        invalid_ids = []

        for col in ["phenotypes_present_ids", "phenotypes_excluded_ids"]:
            for val in df[col].dropna():
                for hpo_id in str(val).split("|"):
                    hpo_id = hpo_id.strip()
                    if not hpo_id or hpo_id == "nan":
                        continue
                    if hpo_id.startswith("HP:"):
                        total_hpo += 1
                        if hpo_id in valid_hpo:
                            valid_count += 1
                        else:
                            if len(invalid_ids) < 20:
                                invalid_ids.append(hpo_id)

        pct = (valid_count / total_hpo * 100) if total_hpo > 0 else 0
        print(f"Total HPO IDs assigned: {total_hpo}")
        print(f"Valid HPO IDs: {valid_count} ({pct:.2f}%)")
        print(f"Invalid HPO IDs: {total_hpo - valid_count}")
        if invalid_ids:
            print(f"First invalid IDs: {invalid_ids[:10]}")
        results["hpo_total"] = total_hpo
        results["hpo_valid"] = valid_count
        results["hpo_valid_pct"] = round(pct, 2)
    else:
        print(f"WARNING: {obo_path} not found, skipping HPO validation")

    print()

    # ---------------------------------------------------------------
    # 2. HGVS format validity
    # ---------------------------------------------------------------
    print("=" * 60)
    print("2. HGVS FORMAT VALIDITY")
    print("=" * 60)

    # Patterns for valid HGVS
    hgvs_c_pattern = re.compile(r"^(NM_|ENST)\S+:c\.\S+")
    hgvs_g_pattern = re.compile(r"^(NC_|chr)\S+:g\.\S+")
    hgvs_p_pattern = re.compile(r"^(NP_|p\.)\S+")
    # Also accept bare c./g./p. notation
    hgvs_bare_pattern = re.compile(r"^[cgpnm]\.\S+")

    for col_name, label in [("variant_cdna", "HGVS c."),
                             ("variant_genomic_grch38", "HGVS g."),
                             ("variant_protein", "HGVS p.")]:
        total = 0
        valid = 0
        invalid_examples = []

        for val in df[col_name].dropna():
            for v in str(val).split("|"):
                v = v.strip()
                if not v or v == "nan":
                    continue
                total += 1
                if (hgvs_c_pattern.match(v) or hgvs_g_pattern.match(v) or
                    hgvs_p_pattern.match(v) or hgvs_bare_pattern.match(v) or
                    v.startswith("NC_") or v.startswith("NM_") or
                    v.startswith("NP_") or v.startswith("NG_")):
                    valid += 1
                else:
                    if len(invalid_examples) < 5:
                        invalid_examples.append(v)

        pct = (valid / total * 100) if total > 0 else 0
        print(f"{label}: {valid}/{total} valid ({pct:.1f}%)")
        if invalid_examples:
            print(f"  Examples of non-standard: {invalid_examples[:5]}")
        results[f"{col_name}_total"] = total
        results[f"{col_name}_valid"] = valid
        results[f"{col_name}_valid_pct"] = round(pct, 2)

    # Also check the variant_hgvs_valid column (regex validation from pipeline)
    if "variant_hgvs_valid" in df.columns:
        valid_flags = df["variant_hgvs_valid"].dropna()
        true_count = (valid_flags == "true").sum()
        false_count = (valid_flags == "false").sum()
        total_flags = true_count + false_count
        pct = (true_count / total_flags * 100) if total_flags > 0 else 0
        print(f"\nPipeline HGVS validation flag: {true_count}/{total_flags} passed ({pct:.1f}%)")
        results["pipeline_hgvs_valid"] = int(true_count)
        results["pipeline_hgvs_total"] = int(total_flags)
        results["pipeline_hgvs_valid_pct"] = round(pct, 2)

    print()

    # ---------------------------------------------------------------
    # 3. Gene symbol validity
    # ---------------------------------------------------------------
    print("=" * 60)
    print("3. GENE SYMBOL VALIDITY")
    print("=" * 60)

    if gene_info_path.exists():
        valid_genes = load_gene_symbols(gene_info_path)
        print(f"Loaded {len(valid_genes)} valid gene symbols from Homo_sapiens.gene_info")

        total_genes = 0
        valid_gene_count = 0
        invalid_genes = []

        for val in df["gene_symbol"].dropna():
            for g in str(val).split("|"):
                g = g.strip()
                if not g or g == "nan":
                    continue
                total_genes += 1
                if g in valid_genes:
                    valid_gene_count += 1
                else:
                    if g not in invalid_genes and len(invalid_genes) < 20:
                        invalid_genes.append(g)

        pct = (valid_gene_count / total_genes * 100) if total_genes > 0 else 0
        print(f"Total gene symbols: {total_genes}")
        print(f"Valid (in NCBI): {valid_gene_count} ({pct:.1f}%)")
        if invalid_genes:
            print(f"Unrecognized symbols ({len(invalid_genes)}): {invalid_genes[:10]}")
        results["gene_total"] = total_genes
        results["gene_valid"] = valid_gene_count
        results["gene_valid_pct"] = round(pct, 2)
    else:
        print(f"WARNING: {gene_info_path} not found, skipping gene validation")

    print()

    # ---------------------------------------------------------------
    # 4. Disease ID validity
    # ---------------------------------------------------------------
    print("=" * 60)
    print("4. DISEASE ID VALIDITY")
    print("=" * 60)

    omim_pattern = re.compile(r"^OMIM:\d+$")
    mondo_pattern = re.compile(r"^MONDO:\d+$")
    orpha_pattern = re.compile(r"^(ORPHA|Orphanet):\d+$")

    total_disease = 0
    valid_disease = 0

    for col in ["disease_omim_ids", "disease_mondo_id", "disease_orphanet_ids"]:
        for val in df[col].dropna():
            for d in str(val).split("|"):
                d = d.strip()
                if not d or d == "nan":
                    continue
                total_disease += 1
                if omim_pattern.match(d) or mondo_pattern.match(d) or orpha_pattern.match(d):
                    valid_disease += 1

    pct = (valid_disease / total_disease * 100) if total_disease > 0 else 0
    print(f"Total disease IDs: {total_disease}")
    print(f"Valid format: {valid_disease} ({pct:.1f}%)")
    results["disease_total"] = total_disease
    results["disease_valid"] = valid_disease
    results["disease_valid_pct"] = round(pct, 2)

    print()

    # ---------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------
    print("=" * 60)
    print("SUMMARY FOR PAPER")
    print("=" * 60)
    print()
    if "hpo_total" in results:
        print(f"HPO terms: {results['hpo_valid']}/{results['hpo_total']} valid ({results['hpo_valid_pct']}%)")
    if "pipeline_hgvs_total" in results:
        print(f"HGVS (pipeline flag): {results['pipeline_hgvs_valid']}/{results['pipeline_hgvs_total']} passed ({results['pipeline_hgvs_valid_pct']}%)")
    if "gene_total" in results:
        print(f"Gene symbols: {results['gene_valid']}/{results['gene_total']} valid ({results['gene_valid_pct']}%)")
    if "disease_total" in results:
        print(f"Disease IDs: {results['disease_valid']}/{results['disease_total']} valid ({results['disease_valid_pct']}%)")


if __name__ == "__main__":
    main()
