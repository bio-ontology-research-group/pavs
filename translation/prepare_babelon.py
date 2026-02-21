import json
import csv
import argparse
import os
import datetime
import pronto

def prepare_babelon(json_path, obo_path, output_path):
    """Convert the translation JSON and HPO OBO to Babelon TSV format."""
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
    
    print(f"Processing {len(translations)} translations...")
    
    # Babelon headers
    fieldnames = [
        'source_language', 'source_value', 'subject_id', 'predicate_id', 
        'translation_language', 'translation_value', 'translation_status', 
        'translator', 'translator_expertise', 'translation_date', 'comment'
    ]
    
    today = datetime.date.today().isoformat()
    
    count = 0
    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='	')
        writer.writeheader()
        
        for trans in translations:
            hp_id = trans['id']
            if hp_id not in ontology:
                continue
            
            term = ontology[hp_id]
            
            # 1. Label
            ar_label = trans.get('arabic_technical_name', '').strip()
            en_label = term.name
            
            if ar_label and ar_label != en_label:
                writer.writerow({
                    'source_language': 'en',
                    'source_value': en_label,
                    'subject_id': hp_id,
                    'predicate_id': 'rdfs:label',
                    'translation_language': 'ar',
                    'translation_value': ar_label,
                    'translation_status': 'CANDIDATE',
                    'translator': 'PAVS-Pipeline (GPT-4o)',
                    'translator_expertise': 'ALGORITHM',
                    'translation_date': today,
                    'comment': 'LLM-generated via PAVS pipeline'
                })
                count += 1
            
            # 2. Definition
            ar_def = trans.get('arabic_definition', '').strip()
            en_def = ""
            if term.definition:
                en_def = str(term.definition).strip()
            
            if ar_def and en_def:
                writer.writerow({
                    'source_language': 'en',
                    'source_value': en_def,
                    'subject_id': hp_id,
                    'predicate_id': 'IAO:0000115',
                    'translation_language': 'ar',
                    'translation_value': ar_def,
                    'translation_status': 'CANDIDATE',
                    'translator': 'PAVS-Pipeline (GPT-4o)',
                    'translator_expertise': 'ALGORITHM',
                    'translation_date': today,
                    'comment': 'LLM-generated via PAVS pipeline'
                })
                count += 1
                
    print(f"Successfully exported {count} Babelon rows to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Prepare HPO Arabic translations in Babelon format.")
    parser.add_argument("--json", default="translation/hpo_arabic_translations.json", help="Path to input JSON")
    parser.add_argument("--obo", default="ontology/hp.obo", help="Path to original hp.obo")
    parser.add_argument("--output", default="translation/hp-ar.babelon.tsv", help="Path to output Babelon TSV")
    
    args = parser.parse_args()
    
    prepare_babelon(args.json, args.obo, args.output)

if __name__ == "__main__":
    main()
