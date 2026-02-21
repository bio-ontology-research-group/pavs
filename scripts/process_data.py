import pandas as pd
import pathlib
import json
import re

def parse_hpo_obo(obo_path):
    hpo_map = {} # name/synonym -> ID
    hpo_id_to_name = {} # ID -> name
    
    with open(obo_path, 'r') as f:
        current_id = None
        current_name = None
        for line in f:
            line = line.strip()
            if line.startswith('id: HP:'):
                current_id = line.split('id: ')[1].split(' ! ')[0]
            elif line.startswith('name: '):
                current_name = line.split('name: ')[1]
                hpo_id_to_name[current_id] = current_name
                hpo_map[current_name.lower()] = current_id
            elif line.startswith('synonym: '):
                synonym = re.findall(r'\"(.+?)\"', line)
                if synonym:
                    hpo_map[synonym[0].lower()] = current_id
    return hpo_map, hpo_id_to_name

def extract_hpo_from_text(text, hpo_map):
    if not isinstance(text, str):
        return [], []
    
    # Pre-processing: Expand abbreviations
    text = re.sub(r'\bDD\b', 'Developmental delay', text)
    text = re.sub(r'\bGDD\b', 'Global developmental delay', text)
    
    found_ids = re.findall(r'HP:\d{7}', text)
    terms = re.split(r'[,|;]\s*', text)
    matched_ids = []
    unmatched_terms = []
    for term in terms:
        clean_term = term.strip().lower()
        # Remove parenthetical HPO IDs
        clean_term = re.sub(r'\(hp:\d{7}\)', '', clean_term).strip()
        if clean_term in hpo_map:
            matched_ids.append(hpo_map[clean_term])
        elif not any(hp_id in term for hp_id in found_ids):
            if clean_term and len(clean_term) > 2:
                unmatched_terms.append(term.strip())
    return list(set(found_ids + matched_ids)), list(set(unmatched_terms))

def normalize_variants(row, source):
    variants = []
    if source == '439_2017':
        v1 = row.get('Variant(s)')
        if pd.notna(v1) and str(v1) != '0': variants.append(str(v1))
        v2 = row.get('Variant(s).1')
        if pd.notna(v2) and str(v2) != '0': variants.append(str(v2))
    elif source == 'Manually':
        v = row.get('variants')
        if pd.notna(v): variants.append(str(v))
    elif source == 'Variant':
        for i in ['', ' 2', ' 3', ' 4']:
            g = row.get(f'Gene{i}'.strip())
            v = row.get(f'Variant{i}'.replace(' ', ''))
            if pd.notna(g) and pd.notna(v):
                variants.append(f"{g}:{v}")
    return "|".join(variants)

def map_sex(sex_str):
    if not isinstance(sex_str, str): return "UNKNOWN_SEX"
    sex_str = sex_str.upper()
    if 'F' in sex_str: return "FEMALE"
    if 'M' in sex_str: return "MALE"
    return "UNKNOWN_SEX"

def main():
    hpo_map, hpo_id_to_name = parse_hpo_obo('ontology/hp.obo')
    master_records = []
    
    # Process 439_2017...
    df1 = pd.read_csv('data/439_2017_1821_MOESM1_ESM.tsv', sep='\t')
    df1.columns = df1.columns.str.strip()
    for idx, row in df1.iterrows():
        hpo_ids, unmatched = extract_hpo_from_text(row['Phenotype'], hpo_map)
        master_records.append({
            'Source': '439_2017_1821_MOESM1_ESM.tsv',
            'Original_ID': row['ID'],
            'Sex': map_sex(row.get('Gender')),
            'Age': row.get('Age'),
            'Zygosity': row.get('Zyogsity'), # Note the typo in source
            'Phenotypes_HPO': "|".join(hpo_ids),
            'Variants': normalize_variants(row, '439_2017'),
            'Diseases': '',
            'Unmatched_Phenotypes': "|".join(unmatched)
        })

    # Process Manually Curated
    df2 = pd.read_csv('data/Manually curated Data.tsv', sep='\t')
    df2.columns = df2.columns.str.strip()
    for idx, row in df2.iterrows():
        hpo_ids, unmatched = extract_hpo_from_text(row['phenotypes'], hpo_map)
        omim = "|".join(re.findall(r'OMIM:\d+', str(row['phenotypes'])))
        master_records.append({
            'Source': 'Manually curated Data.tsv',
            'Original_ID': row['ID'],
            'Sex': 'UNKNOWN_SEX',
            'Age': '',
            'Zygosity': '',
            'Phenotypes_HPO': "|".join(hpo_ids),
            'Variants': normalize_variants(row, 'Manually'),
            'Diseases': omim,
            'Unmatched_Phenotypes': "|".join(unmatched)
        })

    # Process Variant list
    df3 = pd.read_csv('data/Variant list.tsv', sep='\t')
    df3.columns = df3.columns.str.strip()
    for idx, row in df3.iterrows():
        hpo_ids, unmatched = extract_hpo_from_text(row['HPOs'], hpo_map)
        omim = str(row.get('Omim', ''))
        master_records.append({
            'Source': 'Variant list.tsv',
            'Original_ID': row['Case'],
            'Sex': map_sex(row.get('Gender')),
            'Age': row.get('DOB'),
            'Zygosity': row.get('Zygosity'),
            'Phenotypes_HPO': "|".join(hpo_ids),
            'Variants': normalize_variants(row, 'Variant'),
            'Diseases': omim,
            'Unmatched_Phenotypes': "|".join(unmatched)
        })

    master_df = pd.DataFrame(master_records)
    master_df.insert(0, 'Master_ID', [f'PAVS_MASTER_{i+1:04d}' for i in range(len(master_df))])
    master_df.to_csv('master_data_v1.tsv', sep='\t', index=False)
    print(f"Updated master_data_v1.tsv with sex and zygosity. Total: {len(master_df)}")

if __name__ == "__main__":
    main()