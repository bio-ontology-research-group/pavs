#!/usr/bin/env python3
"""Build an expert-grading set for negation and severity-modifier accuracy.

Reviewer 2 noted that the manual review evaluated HPO-ID correctness but not
the accuracy of the LLM-validated negation and modifier assignments, which can
flip the meaning of an annotation even when the HPO ID is right.

This script extracts, from the final phenopackets, every negated phenotype
assignment (``excluded == true``) and a reproducible sample of modifier
assignments, each paired with the raw clinical text that produced it, into a
CSV that two domain experts grade as correct / incorrect.  Because negation is
rare in the corpus (a few dozen assignments), every negation can be graded
exhaustively; modifiers are sampled.

Outputs:
    negation_eval_set.csv   -- all negated assignments + raw context (grade all)
    modifier_eval_set.csv   -- sampled modifier assignments + raw context

Run:
    python evaluation/similarity/scripts/sample_negation_modifier_eval.py \
        --phenopackets data/PAVS_phenopackets.json \
        --raw data/master_data_final.tsv \
        --out evaluation/similarity
"""
import argparse
import csv
import json
import os
import re

SAMPLE_SEED = 20260609   # fixed for reproducibility
MODIFIER_SAMPLE = 100


def load_raw_context(raw_tsv):
    """Map pavs_id -> raw clinical phenotype text, if a raw table is available."""
    ctx = {}
    if not raw_tsv or not os.path.exists(raw_tsv):
        return ctx
    with open(raw_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        id_col = next((c for c in reader.fieldnames if c.lower() in
                       ("pavs_id", "id", "case_id")), None)
        raw_col = next((c for c in reader.fieldnames
                        if "raw" in c.lower() and "phenot" in c.lower()), None)
        if not id_col:
            return ctx
        for row in reader:
            ctx[row[id_col]] = row.get(raw_col, "") if raw_col else ""
    return ctx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phenopackets", default="data/PAVS_phenopackets.json")
    ap.add_argument("--raw", default="data/master_data_final.tsv")
    ap.add_argument("--out", default="evaluation/similarity")
    args = ap.parse_args()

    with open(args.phenopackets) as fh:
        packets = json.load(fh)
    ctx = load_raw_context(args.raw)

    negations, modifiers = [], []
    for pp in packets:
        pid = pp["id"]
        raw = ctx.get(pid, "")
        for feat in pp.get("phenotypicFeatures", []):
            term = feat.get("type", {})
            tid = term.get("id", "")
            if not tid.startswith("HP:"):
                continue
            base = {
                "pavs_id": pid,
                "hpo_id": tid,
                "hpo_label": term.get("label", ""),
                "raw_phenotype_text": raw,
            }
            if feat.get("excluded", False):
                negations.append(dict(base, assignment="NEGATED",
                                      expert_correct="", expert_comment=""))
            for mod in feat.get("modifiers", []) or []:
                modifiers.append(dict(base,
                                      modifier_id=mod.get("id", ""),
                                      modifier_label=mod.get("label", ""),
                                      expert_correct="", expert_comment=""))

    # reproducible modifier sample
    import random
    rng = random.Random(SAMPLE_SEED)
    mod_sample = (modifiers if len(modifiers) <= MODIFIER_SAMPLE
                  else rng.sample(modifiers, MODIFIER_SAMPLE))

    os.makedirs(args.out, exist_ok=True)
    neg_path = os.path.join(args.out, "negation_eval_set.csv")
    mod_path = os.path.join(args.out, "modifier_eval_set.csv")

    def dump(path, rows, fields):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

    dump(neg_path, negations,
         ["pavs_id", "hpo_id", "hpo_label", "assignment",
          "raw_phenotype_text", "expert_correct", "expert_comment"])
    dump(mod_path, mod_sample,
         ["pavs_id", "hpo_id", "hpo_label", "modifier_id", "modifier_label",
          "raw_phenotype_text", "expert_correct", "expert_comment"])

    print(f"Negated assignments (grade ALL): {len(negations)} -> {neg_path}")
    print(f"Modifier assignments total: {len(modifiers)}; "
          f"sampled {len(mod_sample)} -> {mod_path}")
    print("\nExperts: fill 'expert_correct' with 1 (correct) or 0 (incorrect). "
          "Accuracy = mean of expert_correct.")


if __name__ == "__main__":
    main()
