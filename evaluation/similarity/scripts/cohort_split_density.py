#!/usr/bin/env python3
"""Regenerate Table 6 (phenotypic density + IC) with the population/source split.

Mirrors calculate_density_stats.py (extrinsic IC from the HPO gene-to-phenotype
reference; cases with >= 1 HPO term only), but splits the Saudi/mixed cohorts by
PAVS source-letter prefix into Saudi-clinical (A,B,F,P), Saudi-literature (M),
and Mixed (Q), alongside DDD. Literature (Phenopacket Store) is unchanged and
reported from the existing pipeline.

For final manuscript values, point --obo / --g2p at the canonical HPO release
used elsewhere in the paper (2024-08-13); the bundled pyhpo data gives values
that agree with the published aggregates to two decimals.

Run:
    python evaluation/similarity/scripts/cohort_split_density.py \
        --phenopackets data/PAVS_phenopackets.json \
        --obo ontology/hp.obo \
        --g2p data/reference/genes_to_phenotype.txt
"""
import argparse
import json
import math
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

sys.setrecursionlimit(1_000_000)

PREFIX_TO_COHORT = {
    "A": "Saudi-clinical", "B": "Saudi-clinical",
    "F": "Saudi-clinical", "P": "Saudi-clinical",
    "M": "Saudi-literature", "Q": "Mixed", "D": "DDD",
}
ORDER = ["Saudi-clinical", "Saudi-literature", "Mixed", "DDD", "Literature"]


def parse_obo(path):
    terms = {}
    t = None
    for line in open(path):
        line = line.strip()
        if line == "[Term]":
            t = {"id": None, "is_a": []}
        elif line.startswith("id: ") and t is not None:
            t["id"] = line[4:]
        elif line.startswith("is_a: ") and t is not None:
            t["is_a"].append(line[6:].split("!")[0].strip())
        elif line == "" and t is not None:
            if t["id"]:
                terms[t["id"]] = t
            t = None
    return terms


def build_ancestors(terms):
    anc = {}

    def get(t):
        if t in anc:
            return anc[t]
        r = {t}
        for p in terms.get(t, {}).get("is_a", []):
            r |= get(p)
        anc[t] = r
        return r

    for t in list(terms):
        get(t)
    return anc


def extrinsic_ic(terms, anc, g2p):
    df = pd.read_csv(g2p, sep="\t", comment="#")
    gcol = "gene_symbol" if "gene_symbol" in df.columns else df.columns[1]
    hcol = "hpo_id" if "hpo_id" in df.columns else \
        [c for c in df.columns if "hpo" in c.lower()][0]
    g2h = defaultdict(set)
    for _, r in df.iterrows():
        h = str(r[hcol])
        if h in terms:
            g2h[str(r[gcol])].add(h)
    tot = len(g2h)
    cnt = {t: 0 for t in terms}
    for hs in g2h.values():
        ga = set()
        for h in hs:
            ga |= anc.get(h, set())
        for a in ga:
            if a in cnt:
                cnt[a] += 1
    return {t: (-math.log(c / tot) if c > 0 else 0.0) for t, c in cnt.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phenopackets", default="data/PAVS_phenopackets.json")
    ap.add_argument("--obo", default="ontology/hp.obo")
    ap.add_argument("--g2p", default="data/reference/genes_to_phenotype.txt")
    ap.add_argument("--store", default=None,
                    help="Phenopacket Store dir; if given, adds the Literature cohort")
    ap.add_argument("--present-only", action="store_true",
                    help="Count only non-excluded (present) HPO features. By default "
                         "ALL HPO features are counted, matching similarity.py, which "
                         "builds the patient HPO set from every phenotypicFeature "
                         "regardless of its excluded flag.")
    args = ap.parse_args()

    terms = parse_obo(args.obo)
    anc = build_ancestors(terms)
    ic = extrinsic_ic(terms, anc, args.g2p)

    packets = json.load(open(args.phenopackets))
    nh = defaultdict(list)
    ai = defaultdict(list)
    for pp in packets:
        m = re.match(r"PAVS:([A-Za-z]+)", pp["id"])
        cohort = PREFIX_TO_COHORT.get(m.group(1)) if m else None
        if not cohort:
            continue
        hpos = [f["type"]["id"] for f in pp.get("phenotypicFeatures", [])
                if f.get("type", {}).get("id", "").startswith("HP:")
                and (not args.present_only or not f.get("excluded", False))
                and f["type"]["id"] in terms]
        if hpos:  # match calculate_density_stats.py `if hpos`
            nh[cohort].append(len(hpos))
            ai[cohort].append(float(np.mean([ic.get(h, 0.0) for h in hpos])))

    if args.store:
        import os
        for root, _, files in os.walk(args.store):
            for fn in files:
                if not fn.endswith(".json"):
                    continue
                try:
                    sp = json.load(open(os.path.join(root, fn)))
                except Exception:
                    continue
                hpos = [f["type"]["id"] for f in sp.get("phenotypicFeatures", [])
                        if f.get("type", {}).get("id", "").startswith("HP:")
                        and (not args.present_only or not f.get("excluded", False))
                        and f["type"]["id"] in terms]
                if hpos:
                    nh["Literature"].append(len(hpos))
                    ai["Literature"].append(float(np.mean([ic.get(h, 0.0) for h in hpos])))

    print(f"{'Cohort':18s} {'n':>5s} {'MeanHPO':>8s} {'SD':>6s} {'Max':>4s} "
          f"{'MeanIC':>7s} {'IC SD':>6s}")
    for c in ORDER:
        v, i = nh[c], ai[c]
        print(f"{c:18s} {len(v):5d} {np.mean(v):8.2f} {np.std(v, ddof=1):6.2f} "
              f"{max(v):4d} {np.mean(i):7.3f} {np.std(i, ddof=1):6.3f}")


if __name__ == "__main__":
    main()
