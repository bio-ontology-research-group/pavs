#!/usr/bin/env python3
"""Regenerate Tables 5 and 6 with the population/source cohort split.

Reviewer 2 asked us to separate the population axis (Saudi vs. non-Saudi) from
the phenotype-source axis (literature-curated case reports vs. clinically
derived records), because the original "PAVS clinical" cohort merged Saudi
clinical sources, Saudi literature case reports, and the mixed-population Ziats
cohort.

This script splits the Saudi/mixed records by their PAVS source-letter prefix:

    Saudi-clinical    : A, B, F, P   (Alfares, Monies 2017, Monies 2019)
    Saudi-literature  : M            (Abdelhakim curated case reports)
    Mixed             : Q            (Ziats et al., mixed population)
    DDD               : D            (non-Saudi research cohort)
    Literature        : Phenopacket Store (worldwide case reports)

Table 5 (prioritization) is recomputed from the per-case ranking results, using
exactly the same metric definitions as cohort_prioritization_analysis.py:
    Hits@k   = mean(rank <= k)
    AUC      = mean((num_candidates - rank) / (num_candidates - 1))
    MRR      = mean(1 / rank)

Run:
    python evaluation/similarity/scripts/cohort_split_analysis.py \
        --results evaluation/similarity/phenotype_similarity_comprehensive.csv
"""
import argparse
import re

import numpy as np
import pandas as pd

PREFIX_TO_COHORT = {
    "A": "Saudi-clinical", "B": "Saudi-clinical",
    "F": "Saudi-clinical", "P": "Saudi-clinical",
    "M": "Saudi-literature",
    "Q": "Mixed",
    "D": "DDD",
}
COHORT_ORDER = ["Literature", "Saudi-literature", "DDD", "Mixed", "Saudi-clinical"]


def prefix(pavs_id):
    m = re.match(r"PAVS:([A-Za-z]+)", str(pavs_id))
    return m.group(1) if m else None


def assign_cohort(row):
    if row["group"] == "Literature":
        return "Literature"
    p = prefix(row["pavs_id"])
    return PREFIX_TO_COHORT.get(p, f"other:{p}")


def metrics(g):
    r = g["rank"].values
    nc = g["num_candidates"].values
    return {
        "n": len(g),
        "AUC": np.mean((nc - r) / (nc - 1)),
        "MRR": np.mean(1.0 / r),
        "Hits@1": np.mean(r <= 1) * 100,
        "Hits@10": np.mean(r <= 10) * 100,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results",
                    default="evaluation/similarity/phenotype_similarity_comprehensive.csv")
    args = ap.parse_args()

    df = pd.read_csv(args.results)
    df["cohort"] = df.apply(assign_cohort, axis=1)

    print("Table 5 (split) | cohort | method | n | AUC | MRR | Hits@1 | Hits@10")
    rows = []
    for c in COHORT_ORDER:
        for ic in ["extrinsic", "intrinsic"]:
            for sm in ["resnik", "lin"]:
                g = df[(df.cohort == c) & (df.ic_type == ic) & (df.sim_method == sm)]
                if not len(g):
                    continue
                m = metrics(g)
                rows.append({"Cohort": c, "Method": f"{ic.capitalize()} {sm.capitalize()}", **m})
                print(f"  {c:18s} {ic[:3]}-{sm:6s} "
                      f"{m['n']:5d} {m['AUC']:.4f} {m['MRR']:.4f} "
                      f"{m['Hits@1']:.2f} {m['Hits@10']:.2f}")
        print()
    pd.DataFrame(rows).to_csv(
        "evaluation/similarity/cohort_split_prioritization.csv", index=False)
    print("Wrote evaluation/similarity/cohort_split_prioritization.csv")


if __name__ == "__main__":
    main()
