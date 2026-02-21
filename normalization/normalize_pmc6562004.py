import pandas as pd
import re
import os
import json
from normalization.normalization_utils import get_utils
from tqdm import tqdm

utils = get_utils()
OUTPUT_DIR = 'data/processed'
os.makedirs(OUTPUT_DIR, exist_ok=True)

COLUMNS = [
    'id', 'source_file', 
    'sex_id', 'sex_label', 
    'age', 
    'ancestry_id', 'ancestry_label', 
    'country_id', 'country_label', 
    'consanguinity_id', 'consanguinity_label', 
    'family_history_id', 'family_history_label', 
    'hpo_ids', 'hpo_labels', 
    'hpo_excluded_ids', 'hpo_excluded_labels', 
    'disease_omim_ids', 'disease_omim_labels', 
    'disease_omim_excluded_ids', 'disease_omim_excluded_labels',
    'disease_mondo_ids', 'disease_mondo_labels', 
    'gene_ids', 'gene_symbols', 
    'variant_hgvs', 
    'variant_classification', 
    'variant_evidence',
    'genotype_id', 'genotype_label'
]

def parse_age(age_str):
    if pd.isna(age_str): return ""
    s = str(age_str).lower().strip()
    y_match = re.search(r'(\d+)\s*years?', s)
    m_match = re.search(r'(\d+)\s*month?', s)
    
    y = int(y_match.group(1)) if y_match else 0
    m = int(m_match.group(1)) if m_match else 0
    
    if y == 0 and m == 0:
        if s.isdigit(): return f"P{s}Y"
        return ""
    
    res = "P"
    if y > 0: res += f"{y}Y"
    if m > 0: res += f"{m}M"
    return res

def process_pmc6562004():
    input_path = 'data/phenotypes/PMC6562004.tsv'
    if not os.path.exists(input_path): return
    print(f"Processing {input_path}...")
    
    try:
        # Skip comment lines if any, check head
        # The file starts with #Primary indication... but that's a comment
        # Pandas read_csv handles comment='#'
        df = pd.read_csv(input_path, sep='	', comment='#')
        df.columns = df.columns.str.strip()
        
        records = []
        for _, row in tqdm(df.iterrows(), total=len(df)):
            rec = {k: '' for k in COLUMNS}
            rec['id'] = str(row['ID'])
            rec['source_file'] = 'PMC6562004.tsv'
            
            sid, slabel = utils.normalize_sex(row.get('Gender'))
            rec['sex_id'] = sid
            rec['sex_label'] = slabel
            
            rec['age'] = parse_age(row.get('Age'))
            
            aid, alabel = utils.normalize_ancestry()
            rec['ancestry_id'] = aid
            rec['ancestry_label'] = alabel
            cid, clabel = utils.normalize_country()
            rec['country_id'] = cid
            rec['country_label'] = clabel
            
            con_id, con_label = utils.normalize_consanguinity(row.get('Consanguinity'))
            rec['consanguinity_id'] = con_id
            rec['consanguinity_label'] = con_label
            
            fh_id, fh_label = utils.normalize_family_history(row.get('Family hx'))
            rec['family_history_id'] = fh_id
            rec['family_history_label'] = fh_label
            
            genes_found = set()
            variants_hgvs = []
            
            g_raw = row.get('Gene(s)') # Sometimes 'Gene' or 'Gene(s)'
            if pd.isna(g_raw): g_raw = row.get('Gene')
            
            if pd.notna(g_raw):
                g_str = str(g_raw).strip()
                g_norm = utils.map_gene(g_str)
                if g_norm:
                    genes_found.add(g_norm)
                    g_sym = g_norm[1]
                    
                    v_raw = row.get('Variant(s)')
                    if pd.notna(v_raw):
                        # Clean up "NM_...:c...."
                        v_str = str(v_raw).strip()
                        # Often contains "Hom" or other notes?
                        # Using LLM for robust extraction since format varies
                        extracted = utils.extract_variant_hgvs(f"{g_sym} {v_str}", gene_context=g_sym)
                        if extracted:
                            variants_hgvs.append(extracted)
            
            rec['gene_ids'] = "|".join([x[0] for x in genes_found])
            rec['gene_symbols'] = "|".join([x[1] for x in genes_found])
            rec['variant_hgvs'] = "|".join(variants_hgvs)
            
            zyg = row.get('Zygosity')
            if pd.isna(zyg): zyg = row.get('Inheritance AR/AD/XL/NA') # Fallback if zygosity missing?
            gid, glabel = utils.normalize_genotype(zyg)
            rec['genotype_id'] = gid
            rec['genotype_label'] = glabel
            
            pheno_text = str(row.get('Phenotype', ''))
            
            present_hpos = set()
            excluded_hpos = set()
            
            if pheno_text and pheno_text.lower() != 'nan':
                 llm_res_str = utils.llm_map_terms([pheno_text], context_str=f"Gene: {rec['gene_symbols']}")
                 try:
                    mapping = json.loads(llm_res_str)
                    for term_key, data in mapping.items():
                        if data.get('excluded'):
                            if data.get('hpo_id') and data['hpo_id'] != 'None':
                                excluded_hpos.add(data['hpo_id'])
                        else:
                            if data.get('hpo_id') and data['hpo_id'] != 'None':
                                present_hpos.add(data['hpo_id'])
                 except: pass

            rec['hpo_ids'] = "|".join(sorted(list(present_hpos)))
            labels = []
            for h in sorted(list(present_hpos)):
                match = next((t for t in utils.hpo_terms if t['id'] == h), None)
                labels.append(match['name'] if match else h)
            rec['hpo_labels'] = "|".join(labels)
            
            rec['hpo_excluded_ids'] = "|".join(sorted(list(excluded_hpos)))
            
            records.append(rec)

        pd.DataFrame(records).to_csv(os.path.join(OUTPUT_DIR, 'PMC6562004_normalized.tsv'), sep='	', index=False)
        print("Finished PMC6562004.")
        
    except Exception as e:
        print(f"Error processing PMC6562004: {e}")

if __name__ == "__main__":
    process_pmc6562004()
