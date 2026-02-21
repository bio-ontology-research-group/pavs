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

def clean_text(text):
    if pd.isna(text) or str(text).strip() == '?': return ""
    return str(text).strip()

def process_fawzan():
    input_path = 'data/phenotypes/fawzan-variants.tsv'
    if not os.path.exists(input_path): return
    print(f"Processing {input_path}...")
    
    df = pd.read_csv(input_path, sep='\t')
    df.columns = df.columns.str.strip()
    records = []
    
    for _, row in tqdm(df.iterrows(), total=len(df)):
        rec = {k: '' for k in COLUMNS}
        rec['id'] = clean_text(row['ID'])
        rec['source_file'] = 'fawzan-variants.tsv'
        
        # Sex
        sid, slabel = utils.normalize_sex(row.get('Gender'))
        rec['sex_id'] = sid
        rec['sex_label'] = slabel
        
        # Age
        age_val = clean_text(row.get('Age'))
        if age_val and age_val.replace('.', '', 1).isdigit():
            rec['age'] = f"P{age_val}Y"
        elif age_val:
            rec['age'] = age_val # Keep original if not simple number
            
        # Ancestry/Country
        aid, alabel = utils.normalize_ancestry()
        rec['ancestry_id'] = aid
        rec['ancestry_label'] = alabel
        cid, clabel = utils.normalize_country()
        rec['country_id'] = cid
        rec['country_label'] = clabel
        
        # Consanguinity
        con_id, con_label = utils.normalize_consanguinity(row.get('Consanguinity'))
        rec['consanguinity_id'] = con_id
        rec['consanguinity_label'] = con_label
        
        # Family History
        fh_id, fh_label = utils.normalize_family_history(row.get('Family hx'))
        rec['family_history_id'] = fh_id
        rec['family_history_label'] = fh_label
        
        # Genotype
        gid, glabel = utils.normalize_genotype(row.get('Zyogsity')) # Note typo in source column name
        rec['genotype_id'] = gid
        rec['genotype_label'] = glabel
        
        # Variants
        # Format: Gene:Transcript:Exon:c.:p.
        # e.g. CEP152:NM_001194998:exon16:c.2021G>T:p.C674F
        genes_found = set()
        variants_hgvs = []
        
        def parse_variant(v_str):
            if not v_str: return
            parts = v_str.split(':')
            if len(parts) >= 4:
                gene = parts[0]
                # Transcript might be parts[1]
                c_dna = next((p for p in parts if p.startswith('c.')), None)
                
                g_norm = utils.map_gene(gene)
                g_sym = g_norm[1] if g_norm else gene.upper()
                g_id = g_norm[0] if g_norm else ''
                if g_id: genes_found.add((g_id, g_sym))
                
                if c_dna and g_sym:
                    variants_hgvs.append(f"{g_sym}:{c_dna}")
                else:
                    variants_hgvs.append(f"{g_sym}:{v_str}") # Fallback
            else:
                 # Try simple parse if just Gene:c.
                 m = re.match(r'^([^:]+):(c\..+)', v_str)
                 if m:
                     gene, c_dna = m.groups()
                     g_norm = utils.map_gene(gene)
                     g_sym = g_norm[1] if g_norm else gene.upper()
                     g_id = g_norm[0] if g_norm else ''
                     if g_id: genes_found.add((g_id, g_sym))
                     variants_hgvs.append(f"{g_sym}:{c_dna}")

        parse_variant(clean_text(row.get('Variant(s)')))
        parse_variant(clean_text(row.get('Variant(s).1')))
        
        rec['gene_ids'] = "|".join([x[0] for x in genes_found])
        rec['gene_symbols'] = "|".join([x[1] for x in genes_found])
        rec['variant_hgvs'] = "|".join(variants_hgvs)
        
        # Phenotypes
        pheno_text = clean_text(row.get('Phenotype'))
        
        present_hpos = set()
        excluded_hpos = set()
        
        if pheno_text:
            llm_res_str = utils.llm_map_terms([pheno_text], context_str=f"Gene: {rec['gene_symbols']}")
            try:
                mapping = json.loads(llm_res_str)
                # Key is likely full pheno_text if passed as list of 1
                for term_key, data in mapping.items():
                    # Check if 'pheno_text' was split by LLM or returned as one block
                    # If returned as object with keys, process values
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
        
    pd.DataFrame(records).to_csv(os.path.join(OUTPUT_DIR, 'fawzan_normalized.tsv'), sep='\t', index=False)
    print("Finished fawzan.")

if __name__ == "__main__":
    process_fawzan()
