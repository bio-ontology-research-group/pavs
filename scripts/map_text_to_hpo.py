#!/usr/bin/env python3

import argparse
import spacy
from negspacy.negation import Negex
from negspacy.termsets import termset
import pronto
from rapidfuzz import process, fuzz
import unicodedata
import json
import csv

def normalize(text):
    """Normalize text to lowercased, stripped, and Unicode-normalized string."""
    return unicodedata.normalize("NFKC", text).lower().strip()

def load_hpo_obo(hpo_path):
    """
    Load HPO ontology from .obo file and build a dict of normalized names/synonyms → (HPO ID, Label).
    """
    print(f"Loading HPO from {hpo_path} …")
    hpo = pronto.Ontology(hpo_path)
    name_to_id = {}
    for term in hpo.terms():
        if not term.id.startswith("HP:"):
            continue
        all_names = [term.name] + [syn.description for syn in term.synonyms]
        for name in all_names:
            name_to_id[normalize(name)] = (term.id, term.name)
    print(f"Loaded {len(name_to_id):,} names/synonyms → HPO terms.")
    return name_to_id

def extract_phrases(text, nlp):
    """
    Extract candidate phrases from text: hybrid of NER entities (non-negated) + noun chunks,
    filtering out irrelevant terms.
    """
    doc = nlp(text)
    phrases = set()

    # Named entities
    for ent in doc.ents:
        if not ent._.negex:
            phrases.add(normalize(ent.text))
        else:
            print(f"Skipping negated: {ent.text}")

    # Noun chunks
    for chunk in doc.noun_chunks:
        phrases.add(normalize(chunk.text))

    # Filter out unwanted phrases
    blacklist = {
        "patient", "patients", "man", "woman", "male", "female",
        "the patient", "evidence", "person", "subject"
    }
    phrases = [p for p in phrases if len(p) > 3 and p not in blacklist]

    return sorted(phrases)

def match_to_hpo(phrase, name_to_id, threshold=60): #was 70 
    """
    Match a phrase to HPO:
    - exact match if available
    - otherwise fallback to fuzzy matching
    Returns HPO ID or None.
    """
    # Exact match
    if phrase in name_to_id:
        hpo_id, _ = name_to_id[phrase]
        return hpo_id

    # Fuzzy match
    candidates = process.extract(
        phrase, name_to_id.keys(), scorer=fuzz.token_sort_ratio, limit=3
    )
    if candidates and candidates[0][1] >= threshold:
        best, _, _ = candidates[0]
        hpo_id, _ = name_to_id[best]
        return hpo_id
    return None

def process_patient(text, nlp, name_to_id):
    """
    Process a single patient's text → list of HPO IDs.
    """
    phrases = extract_phrases(text, nlp)
    hpo_ids = set()
    for phrase in phrases:
        hpo_id = match_to_hpo(phrase, name_to_id)
        if hpo_id:
            hpo_ids.add(hpo_id)
    return sorted(hpo_ids)

def main():
    parser = argparse.ArgumentParser(
        description="Map TSV of patients + clinical text to JSON of patient → HPO IDs."
    )
    parser.add_argument("-i", "--input", required=True, help="Input TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output JSON file")
    parser.add_argument("--hpo", default="resources/hp.obo", help="Path to HPO .obo file")
    args = parser.parse_args()

    # Load HPO
    name_to_id = load_hpo_obo(args.hpo)

    # Load NLP
    print("Loading NLP model …")
    nlp = spacy.load("resources/en_core_sci_sm-0.5.4/")
    
    ts = termset("en_clinical")
    nlp.add_pipe("negex", config={"neg_termset": ts.get_patterns()})

    # Read TSV & process each patient
    result = {}
    missing_samples = []
    with open(args.input, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) < 2:
                continue
            patient_id, text = row[0], row[1]
            print(f"\nProcessing {patient_id}")
            hpo_list = process_patient(text, nlp, name_to_id)
            result[patient_id] = hpo_list
            print(f"{patient_id} → {hpo_list}")

            if len(hpo_list)<1:
                missing_samples.append(patient_id)

    # Write JSON
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"\n✅ All done — results saved to {args.output}")
    print(f'Missing HPO: {len(missing_samples)}')
    print(missing_samples)

if __name__ == "__main__":
    main()
