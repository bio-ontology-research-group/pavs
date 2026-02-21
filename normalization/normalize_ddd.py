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

def process_ddd():
    input_path = 'data/phenotypes/ddd-diagnoses.tsv'
    if not os.path.exists(input_path): return
    print(f"Processing {input_path}...")
    
    try:
        df = pd.read_csv(input_path, sep='	')
        
        records = []
        for idx, row in tqdm(df.iterrows(), total=len(df)):
            rec = {k: '' for k in COLUMNS}
            rec['id'] = f"DDD_ENTRY_{idx+1}"
            rec['source_file'] = 'ddd-diagnoses.tsv'
            
            # Sex/Age/Ancestry/Country/Consanguinity/FH: Empty
            
            # Genes
            g_str = str(row['gene']).strip()
            g_norm = utils.map_gene(g_str)
            if g_norm:
                rec['gene_ids'] = g_norm[0]
                rec['gene_symbols'] = g_norm[1]
            else:
                rec['gene_symbols'] = g_str
            
            # Variant Info (as classification/evidence)
            rec['variant_classification'] = str(row.get('category', ''))
            rec['variant_evidence'] = str(row.get('mutation_consequence', ''))
            
            # Genotype
            am = str(row.get('allelic_mode', '')).lower()
            if 'biallelic' in am:
                rec['genotype_id'] = 'GENO:0000136'
                rec['genotype_label'] = 'homozygous'
            elif 'monoallelic' in am:
                rec['genotype_id'] = 'GENO:0000135'
                rec['genotype_label'] = 'heterozygous'
            elif 'hemizygous' in am:
                rec['genotype_id'] = 'GENO:0000134'
                rec['genotype_label'] = 'hemizygous'
            elif 'x-linked' in am:
                 rec['genotype_id'] = 'GENO:0000134'
                 rec['genotype_label'] = 'hemizygous' # Assumption for X-linked if sex unknown? Or just leave it?
            
            # Phenotypes
            raw_hpo = str(row.get('hpo_ids', ''))
            hpo_ids = set()
            for h in re.findall(r'HP:\d{7}', raw_hpo):
                hpo_ids.add(h)
            
            rec['hpo_ids'] = "|".join(sorted(list(hpo_ids)))
            labels = []
            for h in sorted(list(hpo_ids)):
                match = next((t for t in utils.hpo_terms if t['id'] == h), None)
                labels.append(match['name'] if match else h)
            rec['hpo_labels'] = "|".join(labels)
            
            # Disease
            oid = str(row.get('omim_id', ''))
            if oid and oid != 'nan':
                rec['disease_omim_ids'] = f"OMIM:{oid}"
                rec['disease_omim_labels'] = str(row.get('syndrome', ''))
            
            records.append(rec)

        pd.DataFrame(records).to_csv(os.path.join(OUTPUT_DIR, 'ddd_diagnoses_normalized.tsv'), sep='	', index=False)
        print("Finished DDD.")
        
    except Exception as e:
        print(f"Error processing DDD: {e}")

if __name__ == "__main__":
    process_ddd()
