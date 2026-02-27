#!/usr/bin/env python3
"""
compute_paper_statistics.py — Compute all statistics needed for the PAVS paper.

Outputs a summary to stdout and writes detailed stats to scripts/paper/paper_stats.json.

Usage:
    uv run python scripts/paper/compute_paper_statistics.py
"""

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent

def safe(val):
    return val is not None and str(val).strip() not in ("", "nan", "None", "NaN")


def main():
    stats = {}

    # -----------------------------------------------------------------------
    # 1. Source file row counts
    # -----------------------------------------------------------------------
    phenotypes_dir = REPO / "data" / "phenotypes"
    source_stats = {}
    for f in sorted(phenotypes_dir.glob("*.tsv")):
        skiprows = 1 if f.name == "PMC6562004.tsv" else 0
        df = pd.read_csv(f, sep="\t", dtype=str, low_memory=False, skiprows=skiprows)
        source_stats[f.stem] = len(df)
    stats["source_files"] = source_stats
    total_source_rows = sum(source_stats.values())
    stats["total_source_rows"] = total_source_rows

    # -----------------------------------------------------------------------
    # 2. Combined normalized TSV stats
    # -----------------------------------------------------------------------
    norm_path = REPO / "data" / "combined_normalized.tsv"
    if norm_path.exists():
        df_norm = pd.read_csv(norm_path, sep="\t", dtype=str, low_memory=False)
        stats["normalized_rows"] = len(df_norm)
        stats["normalized_columns"] = len(df_norm.columns)

        # Count cases by source
        source_counts = df_norm["source_file"].value_counts().to_dict()
        stats["cases_by_source"] = source_counts

        # Saudi vs non-Saudi
        saudi_sources = {"ahmed-variants", "ahmed-pmid28454995", "fawzan-variants",
                         "marwa-variants", "PMC6562004", "PMC7082194"}
        saudi_mask = df_norm["source_file"].isin(saudi_sources)
        stats["saudi_cases"] = int(saudi_mask.sum())
        stats["non_saudi_cases"] = int((~saudi_mask).sum())

        # Unique genes
        all_genes = set()
        for g in df_norm["gene_symbol"].dropna():
            for gs in str(g).split("|"):
                gs = gs.strip()
                if gs and gs != "nan":
                    all_genes.add(gs)
        stats["unique_genes"] = len(all_genes)

        # Unique diseases (OMIM + MONDO)
        all_diseases = set()
        for col in ["disease_omim_ids", "disease_mondo_id"]:
            for val in df_norm[col].dropna():
                for d in str(val).split("|"):
                    d = d.strip()
                    if d and d != "nan":
                        all_diseases.add(d)
        stats["unique_diseases"] = len(all_diseases)

        # HPO stats
        hpo_counts = []
        for val in df_norm["phenotypes_present_ids"].dropna():
            ids = [x.strip() for x in str(val).split("|") if x.strip().startswith("HP:")]
            hpo_counts.append(len(ids))
        if hpo_counts:
            stats["mean_hpo_per_case"] = round(sum(hpo_counts) / len(hpo_counts), 2)
            stats["median_hpo_per_case"] = sorted(hpo_counts)[len(hpo_counts) // 2]
            stats["max_hpo_per_case"] = max(hpo_counts)
            stats["cases_with_hpo"] = sum(1 for c in hpo_counts if c > 0)
            stats["cases_without_hpo"] = sum(1 for c in hpo_counts if c == 0)

        # Sex distribution
        sex_counts = df_norm["sex_label"].value_counts().to_dict()
        stats["sex_distribution"] = sex_counts

        # Zygosity distribution
        zyg_counts = df_norm["zygosity_label"].value_counts().to_dict()
        stats["zygosity_distribution"] = zyg_counts

        # ACMG distribution
        acmg_counts = df_norm["acmg_classification"].value_counts().to_dict()
        stats["acmg_distribution"] = acmg_counts

        # Consanguinity
        consang_counts = df_norm["consanguinity_label"].value_counts().to_dict()
        stats["consanguinity_distribution"] = consang_counts

        # Excluded phenotypes
        excluded_count = 0
        for val in df_norm["phenotypes_excluded_ids"].dropna():
            ids = [x.strip() for x in str(val).split("|") if x.strip().startswith("HP:")]
            excluded_count += len(ids)
        stats["total_excluded_phenotypes"] = excluded_count

        # Variants with HGVS
        hgvs_c_count = df_norm["variant_cdna"].dropna().apply(lambda x: safe(x)).sum()
        hgvs_g_count = df_norm["variant_genomic_grch38"].dropna().apply(lambda x: safe(x)).sum()
        stats["variants_with_hgvs_c"] = int(hgvs_c_count)
        stats["variants_with_hgvs_g"] = int(hgvs_g_count)

        # Top genes by frequency
        gene_counter = Counter()
        for g in df_norm["gene_symbol"].dropna():
            for gs in str(g).split("|"):
                gs = gs.strip()
                if gs and gs != "nan":
                    gene_counter[gs] += 1
        stats["top_20_genes"] = gene_counter.most_common(20)

        # Result status
        result_counts = df_norm["result_status"].value_counts().to_dict()
        stats["result_status"] = result_counts

    # -----------------------------------------------------------------------
    # 3. Combined annotated TSV stats (the VEP-augmented version)
    # -----------------------------------------------------------------------
    ann_path = REPO / "data" / "combined_annotated.tsv"
    if ann_path.exists():
        df_ann = pd.read_csv(ann_path, sep="\t", dtype=str, low_memory=False)
        stats["annotated_rows"] = len(df_ann)
        stats["annotated_columns"] = len(df_ann.columns)

        # VEP annotation coverage
        vep_cols = ["vep_consequence", "vep_rsid", "vep_sift", "vep_polyphen",
                    "vep_cadd_phred", "vep_revel", "vep_alphamissense"]
        vep_coverage = {}
        for col in vep_cols:
            if col in df_ann.columns:
                vep_coverage[col] = int(df_ann[col].dropna().apply(lambda x: safe(x)).sum())
        stats["vep_annotation_coverage"] = vep_coverage

        # ClinVar coverage
        if "clinvar_allele_id" in df_ann.columns:
            stats["clinvar_annotated"] = int(df_ann["clinvar_allele_id"].dropna().apply(lambda x: safe(x)).sum())

    # -----------------------------------------------------------------------
    # 4. Phenopackets stats
    # -----------------------------------------------------------------------
    pp_path = REPO / "data" / "PAVS_phenopackets.json"
    if pp_path.exists():
        with open(pp_path) as f:
            pps = json.load(f)
        stats["total_phenopackets"] = len(pps)

        # Count with diseases, variants, phenotypes
        with_disease = sum(1 for p in pps if p.get("diseases"))
        with_variant = sum(1 for p in pps if p.get("interpretations"))
        with_pheno = sum(1 for p in pps
                         if any(f.get("type", {}).get("id", "").startswith("HP:")
                                for f in p.get("phenotypicFeatures", [])))
        stats["phenopackets_with_disease"] = with_disease
        stats["phenopackets_with_variant"] = with_variant
        stats["phenopackets_with_phenotype"] = with_pheno

    # -----------------------------------------------------------------------
    # 5. RDF / Knowledge graph stats
    # -----------------------------------------------------------------------
    rdf_dir = REPO / "rdf_output"
    if rdf_dir.exists():
        rdf_sizes = {}
        for f in sorted(rdf_dir.glob("*.ttl")):
            line_count = sum(1 for _ in open(f))
            rdf_sizes[f.name] = line_count
        stats["rdf_line_counts"] = rdf_sizes

    # -----------------------------------------------------------------------
    # 6. Literature phenopackets
    # -----------------------------------------------------------------------
    lit_dir = REPO / "phenopackets" / "0.1.26"
    if lit_dir.exists():
        lit_count = len(list(lit_dir.rglob("*.json")))
        stats["literature_phenopackets"] = lit_count

    # -----------------------------------------------------------------------
    # 7. Analysis results
    # -----------------------------------------------------------------------
    analysis_dir = REPO / "analysis"

    # Comprehensive results
    comp_path = analysis_dir / "phenotype_similarity_comprehensive.csv"
    if comp_path.exists():
        df_comp = pd.read_csv(comp_path)
        stats["analysis_total_evaluations"] = len(df_comp)

        # Best method: extrinsic resnik
        for ic_type in ["intrinsic", "extrinsic"]:
            for method in ["lin", "resnik"]:
                subset = df_comp[(df_comp["ic_type"] == ic_type) & (df_comp["sim_method"] == method)]
                for group in subset["group"].unique():
                    grp = subset[subset["group"] == group]
                    n = len(grp)
                    hits1 = (grp["rank"] == 1).sum() / n * 100 if n > 0 else 0
                    hits10 = (grp["rank"] <= 10).sum() / n * 100 if n > 0 else 0
                    mrr = (1.0 / grp["rank"]).mean() if n > 0 else 0
                    key = f"perf_{ic_type}_{method}_{group}"
                    stats[key] = {
                        "n": n,
                        "hits_at_1": round(hits1, 2),
                        "hits_at_10": round(hits10, 2),
                        "mrr": round(mrr, 4),
                    }

    # Density stats
    density_path = analysis_dir / "phenotypic_density_ic_stats.csv"
    if density_path.exists():
        df_density = pd.read_csv(density_path)
        stats["density_stats"] = df_density.to_dict(orient="records")

    # -----------------------------------------------------------------------
    # Output
    # -----------------------------------------------------------------------
    out_path = REPO / "scripts" / "paper" / "paper_stats.json"
    with open(out_path, "w") as f:
        json.dump(stats, f, indent=2, default=str)
    print(f"Stats written to {out_path}")

    # Print summary
    print("\n=== PAVS Paper Statistics ===\n")
    print(f"Source files: {len(source_stats)}")
    for name, count in source_stats.items():
        print(f"  {name}: {count} rows")
    print(f"Total source rows: {total_source_rows}")

    if "normalized_rows" in stats:
        print(f"\nNormalized TSV: {stats['normalized_rows']} rows × {stats['normalized_columns']} columns")
        print(f"Saudi cases: {stats.get('saudi_cases', '?')}")
        print(f"Non-Saudi (DDD): {stats.get('non_saudi_cases', '?')}")
        print(f"Unique genes: {stats.get('unique_genes', '?')}")
        print(f"Unique diseases: {stats.get('unique_diseases', '?')}")
        print(f"Mean HPO/case: {stats.get('mean_hpo_per_case', '?')}")
        print(f"Max HPO/case: {stats.get('max_hpo_per_case', '?')}")
        print(f"Cases with HPO: {stats.get('cases_with_hpo', '?')}")
        print(f"\nSex: {stats.get('sex_distribution', {})}")
        print(f"Zygosity: {stats.get('zygosity_distribution', {})}")
        print(f"ACMG: {stats.get('acmg_distribution', {})}")
        print(f"Consanguinity: {stats.get('consanguinity_distribution', {})}")
        print(f"Result status: {stats.get('result_status', {})}")

        if "top_20_genes" in stats:
            print("\nTop 20 genes:")
            for gene, count in stats["top_20_genes"]:
                print(f"  {gene}: {count}")

    if "annotated_rows" in stats:
        print(f"\nAnnotated TSV: {stats['annotated_rows']} rows × {stats['annotated_columns']} columns")
        print(f"VEP coverage: {stats.get('vep_annotation_coverage', {})}")
        print(f"ClinVar annotated: {stats.get('clinvar_annotated', '?')}")

    if "total_phenopackets" in stats:
        print(f"\nPhenopackets: {stats['total_phenopackets']}")
        print(f"  With disease: {stats.get('phenopackets_with_disease', '?')}")
        print(f"  With variant: {stats.get('phenopackets_with_variant', '?')}")
        print(f"  With phenotype: {stats.get('phenopackets_with_phenotype', '?')}")

    if "rdf_line_counts" in stats:
        print(f"\nRDF files:")
        for name, lines in stats["rdf_line_counts"].items():
            print(f"  {name}: {lines:,} lines")

    if "literature_phenopackets" in stats:
        print(f"\nLiterature phenopackets: {stats['literature_phenopackets']}")

    # Analysis performance
    for key in sorted(stats.keys()):
        if key.startswith("perf_extrinsic_resnik_"):
            group = key.replace("perf_extrinsic_resnik_", "")
            v = stats[key]
            print(f"\nExtrinsic Resnik — {group}: n={v['n']}, Hits@1={v['hits_at_1']}%, Hits@10={v['hits_at_10']}%, MRR={v['mrr']}")


if __name__ == "__main__":
    main()
