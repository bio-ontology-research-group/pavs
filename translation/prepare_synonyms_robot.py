import json
import csv
import argparse
import os
import pronto

def prepare_synonyms_robot(json_path, obo_path, output_path):
    """Convert the translation JSON and HPO OBO to ROBOT template format for synonyms."""
    if not os.path.exists(json_path):
        print(f"Error: {json_path} not found.")
        return
    if not os.path.exists(obo_path):
        print(f"Error: {obo_path} not found.")
        return

    print(f"Loading HPO from {obo_path}...")
    ontology = pronto.Ontology(obo_path)
    
    with open(json_path, 'r', encoding='utf-8') as f:
        translations = json.load(f)
    
    print(f"Processing {len(translations)} synonyms...")
    
    # ROBOT template headers
    # Row 1: Human-readable headers
    # Row 2: ROBOT template strings
    fieldnames = ['ID', 'LABEL', 'Related Synonym']
    
    count = 0
    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='	')
        writer.writerow({'ID': 'ID', 'LABEL': 'LABEL', 'Related Synonym': 'Related Synonym'})
        writer.writerow({'ID': 'ID', 'LABEL': 'LABEL', 'Related Synonym': 'A oboInOwl:hasRelatedSynonym'})
        
        for trans in translations:
            hp_id = trans['id']
            if hp_id not in ontology:
                continue
            
            term = ontology[hp_id]
            ar_lay = trans.get('arabic_layperson_synonym', '').strip()
            en_label = term.name
            
            if ar_lay and ar_lay != en_label:
                writer.writerow({
                    'ID': hp_id,
                    'LABEL': en_label,
                    'Related Synonym': ar_lay
                })
                count += 1
                
    print(f"Successfully exported {count} synonyms to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Prepare HPO Arabic synonyms in ROBOT template format.")
    parser.add_argument("--json", default="translation/hpo_arabic_translations.json", help="Path to input JSON")
    parser.add_argument("--obo", default="ontology/hp.obo", help="Path to original hp.obo")
    parser.add_argument("--output", default="translation/hp-ar-synonyms.robot.tsv", help="Path to output ROBOT TSV")
    
    args = parser.parse_args()
    
    prepare_synonyms_robot(args.json, args.obo, args.output)

if __name__ == "__main__":
    main()
