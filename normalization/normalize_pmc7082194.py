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
    y_match = re.search(r'(\d+)\s*y', s)
    m_match = re.search(r'(\d+)\s*m', s)
    
    y = int(y_match.group(1)) if y_match else 0
    m = int(m_match.group(1)) if m_match else 0
    
    if y == 0 and m == 0:
        if s.isdigit(): return f"P{s}Y"
        return ""
    
    res = "P"
    if y > 0: res += f"{y}Y"
    if m > 0: res += f"{m}M"
    return res

def process_pmc7082194():
    input_path = 'data/phenotypes/PMC7082194.tsv'
    if not os.path.exists(input_path): return
    print(f"Processing {input_path}...")
    
    try:
        df = pd.read_csv(input_path, sep='\t')
        
        # Split into two DFs (Left: 0-9, Right: 10-13)
        df_gen = df.iloc[:, :10].copy()
        df_phen = df.iloc[:, 10:14].copy()
        df_phen.columns = ['ID', 'Sex', 'Age', 'Phenotype']
        
        df_gen['ID_clean'] = pd.to_numeric(df_gen['ID'], errors='coerce').fillna(-1).astype(int).astype(str)
        df_phen['ID_clean'] = pd.to_numeric(df_phen['ID'], errors='coerce').fillna(-1).astype(int).astype(str)
        
        merged = df_phen.merge(df_gen, on='ID_clean', how='left', suffixes=('_p', '_g'))
        
        records = []
        for _, row in tqdm(merged.iterrows(), total=len(merged)):
            if row['ID_clean'] == '-1': continue
            
            rec = {k: '' for k in COLUMNS}
            rec['id'] = str(row['ID_clean'])
            rec['source_file'] = 'PMC7082194.tsv'
            
            sid, slabel = utils.normalize_sex(row.get('Sex_p'))
            rec['sex_id'] = sid
            rec['sex_label'] = slabel
            
            rec['age'] = parse_age(row.get('Age_p'))
            
            aid, alabel = utils.normalize_ancestry()
            rec['ancestry_id'] = aid
            rec['ancestry_label'] = alabel
            cid, clabel = utils.normalize_country()
            rec['country_id'] = cid
            rec['country_label'] = clabel
            
            genes_found = set()
            variants_hgvs = []
            
            g_raw = row.get('Gene(s)')
            g_sym = ""
            if pd.notna(g_raw):
                g_str = str(g_raw).strip()
                g_norm = utils.map_gene(g_str)
                if g_norm:
                    genes_found.add(g_norm)
                    g_sym = g_norm[1]
                    
                    hgvs_raw = row.get('HGVS Format')
                    if pd.notna(hgvs_raw):
                        v_str = str(hgvs_raw).strip()
                        if v_str.startswith('NM_') or v_str.startswith('c.'):
                             if not v_str.startswith(g_sym) and ':' in v_str:
                                 variants_hgvs.append(f"{g_sym}:{v_str.split(':')[1]}")
                             elif ':' in v_str:
                                 variants_hgvs.append(f"{g_sym}:{v_str.split(':', 1)[1]}")
                             else:
                                 variants_hgvs.append(f"{g_sym}:{v_str}")
                    else:
                        nuc = row.get('Nucleotide Change')
                        if pd.notna(nuc):
                             variants_hgvs.append(f"{g_sym}:{str(nuc).strip()}")
            
            rec['gene_ids'] = "|".join([x[0] for x in genes_found])
            rec['gene_symbols'] = "|".join([x[1] for x in genes_found])
            rec['variant_hgvs'] = "|".join(variants_hgvs)
            
            inh = str(row.get('Inheritance (if known)', '')).lower()
            if 'homo' in inh:
                rec['genotype_id'] = 'GENO:0000136'
                rec['genotype_label'] = 'homozygous'
            elif 'het' in inh:
                rec['genotype_id'] = 'GENO:0000135'
                rec['genotype_label'] = 'heterozygous'
            else:
                rec['genotype_id'] = 'GENO:0000137'
                rec['genotype_label'] = 'unspecified zygosity'

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

        pd.DataFrame(records).to_csv(os.path.join(OUTPUT_DIR, 'PMC7082194_normalized.tsv'), sep='\t', index=False)
        print("Finished PMC7082194.")
        
    except Exception as e:
        print(f"Error processing PMC7082194: {e}")

if __name__ == "__main__":
    process_pmc7082194()
