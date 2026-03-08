#!/usr/bin/env python3
"""
compute_hpo_ic.py — Compute HPO information content from phenotype.hpoa
and write Turtle triples for Virtuoso ingestion.

Usage:
    uv run python scripts/compute_hpo_ic.py \
        --hpo ontology/hp.obo \
        --hpoa data/phenotype.hpoa \
        --output rdf_output/hpo_ic.ttl
"""

import argparse
import json
import math
import re
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Parse hp.obo to build the parent→children and term→ancestors maps
# ---------------------------------------------------------------------------

def parse_hpo_obo(obo_path: str) -> tuple[dict, dict, dict]:
    """Return (parents, children, names) dicts mapping HP:NNNNNNN → ..."""
    parents: dict[str, set] = defaultdict(set)
    children: dict[str, set] = defaultdict(set)
    names: dict[str, str] = {}

    current_id = None
    in_term = False

    with open(obo_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if line == "[Term]":
                in_term = True
                current_id = None
            elif line.startswith("[") and line != "[Term]":
                in_term = False
            elif in_term:
                if line.startswith("id: "):
                    raw = line[4:].strip()
                    current_id = raw.replace("HP:", "HP:")
                elif line.startswith("name: ") and current_id:
                    names[current_id] = line[6:].strip()
                elif line.startswith("is_a: ") and current_id:
                    parent = line[6:].split("!")[0].strip()
                    parents[current_id].add(parent)
                    children[parent].add(current_id)
                elif line.startswith("is_obsolete: true"):
                    current_id = None

    return parents, children, names


def build_ancestor_sets(parents: dict) -> dict[str, set]:
    """For each term return the set of all ancestors (including self)."""
    ancestors: dict[str, set] = {}

    def get_ancestors(term: str) -> set:
        if term in ancestors:
            return ancestors[term]
        anc = {term}
        for p in parents.get(term, []):
            anc |= get_ancestors(p)
        ancestors[term] = anc
        return anc

    all_terms = set(parents.keys())
    for p_set in parents.values():
        all_terms |= p_set

    for t in all_terms:
        get_ancestors(t)

    return ancestors


# ---------------------------------------------------------------------------
# Parse phenotype.hpoa
# ---------------------------------------------------------------------------

def parse_hpoa(hpoa_path: str) -> dict[str, set]:
    """Return disease_id → set of HP term IDs."""
    disease_hpos: dict[str, set] = defaultdict(set)
    with open(hpoa_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("database_id"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            disease_id = parts[0]   # e.g. OMIM:272200
            qualifier = parts[2]    # "NOT" or ""
            hpo_id = parts[3]       # HP:0001263
            if qualifier != "NOT" and hpo_id.startswith("HP:"):
                disease_hpos[disease_id].add(hpo_id)
    return disease_hpos


# ---------------------------------------------------------------------------
# Compute IC
# ---------------------------------------------------------------------------

def compute_ic(hpoa_path: str, ancestors: dict[str, set]) -> dict[str, float]:
    """Compute IC(t) = -log2(freq(t) / total_diseases)."""
    disease_hpos = parse_hpoa(hpoa_path)
    total_diseases = len(disease_hpos)
    if total_diseases == 0:
        raise ValueError("No disease associations found in HPOA file")

    # Propagate each term up to ancestors
    term_disease_count: dict[str, set] = defaultdict(set)
    for disease_id, hpo_set in disease_hpos.items():
        for hpo in hpo_set:
            for anc in ancestors.get(hpo, {hpo}):
                term_disease_count[anc].add(disease_id)

    ic: dict[str, float] = {}
    for term, diseases in term_disease_count.items():
        freq = len(diseases) / total_diseases
        ic[term] = -math.log2(freq) if freq > 0 else 0.0

    return ic


# ---------------------------------------------------------------------------
# Write Turtle
# ---------------------------------------------------------------------------

PREFIXES = """\
@prefix pavs:  <http://pavs.kaust.edu.sa/ontology/> .
@prefix hp:    <http://purl.obolibrary.org/obo/HP_> .
@prefix rdfs:  <http://www.w3.org/2000/01/rdf-schema#> .
@prefix xsd:   <http://www.w3.org/2001/XMLSchema#> .

"""


def hp_uri(hpo_id: str) -> str:
    """Convert HP:0001263 → hp:0001263 (using prefix hp:)."""
    return "hp:" + hpo_id.replace("HP:", "")


def load_arabic_labels(json_path: str) -> dict[str, str]:
    """Load hpo_arabic_translations.json → {HP:ID → arabic_technical_name}."""
    result: dict[str, str] = {}
    try:
        with open(json_path, encoding="utf-8") as fh:
            terms = json.load(fh)
        for term in terms:
            hp_id = term.get("id", "")
            label = term.get("arabic_technical_name", "").strip()
            if hp_id and label:
                result[hp_id] = label
        print(f"  Loaded {len(result)} Arabic HPO labels from {json_path}")
    except Exception as e:
        print(f"  Warning: could not load Arabic labels: {e}")
    return result


def write_ttl(ic: dict[str, float], parents: dict, children: dict, names: dict,
              output_path: str, arabic: dict[str, str] | None = None):
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(PREFIXES)

        # IC values
        for term, value in sorted(ic.items()):
            uri = hp_uri(term)
            fh.write(f"{uri}  pavs:informationContent  {value:.6f} .\n")

        # English labels
        fh.write("\n# HPO term labels\n")
        for term, name in sorted(names.items()):
            uri = hp_uri(term)
            escaped = name.replace('\\', '\\\\').replace('"', '\\"')
            fh.write(f'{uri}  rdfs:label  "{escaped}" .\n')

        # Arabic labels (language-tagged so SPARQL can match Arabic queries)
        if arabic:
            fh.write("\n# Arabic HPO term labels\n")
            for term in sorted(arabic.keys()):
                uri = hp_uri(term)
                escaped = arabic[term].replace('\\', '\\\\').replace('"', '\\"')
                fh.write(f'{uri}  rdfs:label  "{escaped}"@ar .\n')

        # Parent/child hierarchy (rdfs:subClassOf)
        fh.write("\n# HPO hierarchy\n")
        for term, parent_set in sorted(parents.items()):
            uri = hp_uri(term)
            for parent in sorted(parent_set):
                pu = hp_uri(parent)
                fh.write(f"{uri}  rdfs:subClassOf  {pu} .\n")

    ar_count = len(arabic) if arabic else 0
    print(f"Wrote {len(ic)} IC values + {len(names)} EN labels + {ar_count} AR labels + hierarchy to {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compute HPO IC values → Turtle")
    parser.add_argument("--hpo", default="ontology/hp.obo", help="hp.obo path")
    parser.add_argument("--hpoa", default="data/reference/phenotype.hpoa", help="phenotype.hpoa path")
    parser.add_argument("--output", default="rdf_output/hpo_ic.ttl", help="Output TTL path")
    parser.add_argument("--arabic", default="translation/hpo_arabic_translations.json",
                        help="Arabic HPO translations JSON (optional)")
    args = parser.parse_args()

    print("Parsing hp.obo …")
    parents, children, names = parse_hpo_obo(args.hpo)
    print(f"  {len(parents)} HPO terms with parents, {len(names)} term names")

    print("Building ancestor sets …")
    ancestors = build_ancestor_sets(parents)

    print("Computing IC from phenotype.hpoa …")
    ic = compute_ic(args.hpoa, ancestors)
    print(f"  {len(ic)} terms with IC values")

    arabic: dict[str, str] = {}
    if args.arabic and Path(args.arabic).exists():
        print("Loading Arabic HPO labels …")
        arabic = load_arabic_labels(args.arabic)

    write_ttl(ic, parents, children, names, args.output, arabic or None)


if __name__ == "__main__":
    main()
