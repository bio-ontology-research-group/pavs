import json
import csv
import argparse
import os

def json_to_tsv(json_path, tsv_path):
    """Convert the translation JSON to a TSV file."""
    if not os.path.exists(json_path):
        print(f"Error: {json_path} not found.")
        return

    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    if not data:
        print("No data to export.")
        return

    fieldnames = ['id', 'english_technical_name', 'arabic_technical_name', 'arabic_layperson_synonym', 'arabic_definition']
    
    with open(tsv_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in data:
            filtered_row = {k: row.get(k, '') for k in fieldnames}
            writer.writerow(filtered_row)
    print(f"Successfully exported TSV to {tsv_path}")

def update_obo(json_path, input_obo, output_obo):
    """
    Update the OBO file by inserting Arabic synonyms.
    Uses line-by-line processing to preserve original formatting and handle large files.
    """
    if not os.path.exists(json_path):
        print(f"Error: {json_path} not found.")
        return
    if not os.path.exists(input_obo):
        print(f"Error: {input_obo} not found.")
        return

    with open(json_path, 'r', encoding='utf-8') as f:
        translations = json.load(f)
    
    trans_map = {t['id']: t for t in translations}
    
    print(f"Processing OBO from {input_obo}...")
    
    with open(input_obo, 'r', encoding='iso-8859-1') as f_in:
        with open(output_obo, 'w', encoding='utf-8') as f_out:
            current_id = None
            for line in f_in:
                f_out.write(line)
                line_s = line.strip()
                
                if line_s.startswith('id: '):
                    current_id = line_s[4:].split(' ! ')[0]
                
                if line_s.startswith('name: ') and current_id in trans_map:
                    trans = trans_map[current_id]
                    
                    ar_tech = trans.get('arabic_technical_name', '').strip()
                    en_tech = trans.get('english_technical_name', '').strip()
                    if ar_tech and ar_tech != en_tech:
                        f_out.write(f'synonym: "{ar_tech}" EXACT [PAVS:AR]\n')
                    
                    ar_lay = trans.get('arabic_layperson_synonym', '').strip()
                    if ar_lay and ar_lay not in [en_tech, ar_tech]:
                        f_out.write(f'synonym: "{ar_lay}" RELATED [PAVS:AR]\n')
                    
                    ar_def = trans.get('arabic_definition', '').strip()
                    if ar_def:
                        f_out.write(f'comment: Arabic Definition: {ar_def}\n')

    print(f"Successfully generated updated OBO: {output_obo}")

def main():
    parser = argparse.ArgumentParser(description="Export HPO Arabic translations to TSV and updated OBO.")
    parser.add_argument("--json", default="translation/hpo_arabic_translations.json", help="Path to input JSON")
    parser.add_argument("--obo-in", default="ontology/hp.obo", help="Path to original hp.obo")
    parser.add_argument("--obo-out", default="translation/hp-ar.obo", help="Path to output hp-ar.obo")
    parser.add_argument("--tsv", default="translation/hpo_arabic_translations.tsv", help="Path to output TSV")
    
    args = parser.parse_args()
    
    json_to_tsv(args.json, args.tsv)
    update_obo(args.json, args.obo_in, args.obo_out)

if __name__ == "__main__":
    main()
