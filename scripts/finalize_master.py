import pandas as pd
import json
import re
import os

def load_mappings():
    if os.path.exists('hpo_mappings_final.json'):
        with open('hpo_mappings_final.json', 'r') as f:
            return json.load(f)
    elif os.path.exists('hpo_mappings_partial.json'):
        with open('hpo_mappings_partial.json', 'r') as f:
            return json.load(f)
    return {}

def clean_variant(v):
    if not isinstance(v, str) or v == 'nan' or not v:
        return ""
    parts = v.split('|')
    cleaned_parts = []
    for p in parts:
        nm = re.search(r'(NM_\d+\.\d+)', p)
        c_dna = re.search(r'(c\.[^:\s,;]+)', p)
        if nm and c_dna:
            cleaned_parts.append(f"{nm.group(1)}:{c_dna.group(1)}")
        else:
            cleaned_parts.append(p.strip())
    return "|".join(cleaned_parts)

def map_zygosity(z):
    if not isinstance(z, str): return ""
    z = z.lower()
    if 'hom' in z: return "GENO:0000136" # homozygous
    if 'het' in z: return "GENO:0000135" # heterozygous
    if 'hemi' in z: return "GENO:0000606" # hemizygous
    return ""

def main():
    df = pd.read_csv('master_data_v1.tsv', sep='\t')
    mappings = load_mappings()
    
    final_records = []
    
    for _, row in df.iterrows():
        hpos = set(str(row['Phenotypes_HPO']).split('|')) if pd.notna(row['Phenotypes_HPO']) and row['Phenotypes_HPO'] != "" else set()
        unmatched = str(row['Unmatched_Phenotypes']).split('|') if pd.notna(row['Unmatched_Phenotypes']) and row['Unmatched_Phenotypes'] != "" else []
        
        for term in unmatched:
            mapped_id = mappings.get(term)
            if mapped_id and mapped_id != "None" and mapped_id.startswith("HP:"):
                hpos.add(mapped_id)
        
        hpos = {h for h in hpos if h.startswith('HP:')}
        
        final_records.append({
            'Master_ID': row['Master_ID'],
            'Source': row['Source'],
            'Original_ID': row['Original_ID'],
            'Sex': row['Sex'],
            'Age': row['Age'],
            'Zygosity_GENO': map_zygosity(row['Zygosity']),
            'Phenotypes_HPO': "|".join(sorted(list(hpos))),
            'Variants_HGVS': clean_variant(row['Variants']),
            'Genome_Build': 'GRCh37',
            'Diseases': row['Diseases'],
            'Original_Unmatched': row['Unmatched_Phenotypes']
        })
    
    final_df = pd.DataFrame(final_records)
    final_df.to_csv('master_data_final.tsv', sep='\t', index=False)
    print(f"Final master data saved with {len(final_df)} records (v2.0 fields included).")

if __name__ == "__main__":
    main()