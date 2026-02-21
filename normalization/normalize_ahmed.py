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
    'id', 'source_file', 'sex_id', 'sex_label', 'age', 
    'ancestry_id', 'ancestry_label', 'country_id', 'country_label', 
    'consanguinity_id', 'consanguinity_label', 'family_history_id', 'family_history_label', 
    'hpo_ids', 'hpo_labels', 'hpo_excluded_ids', 'hpo_excluded_labels', 
    'disease_omim_ids', 'disease_omim_labels', 'disease_omim_excluded_ids', 'disease_omim_excluded_labels',
    'disease_mondo_ids', 'disease_mondo_labels', 'disease_mondo_excluded_ids', 'disease_mondo_excluded_labels',
    'gene_ids', 'gene_symbols', 'variant_hgvs', 'variant_classification', 'variant_evidence', 'genotype_id', 'genotype_label'
]

def process_ahmed_variants():
    input_path = 'data/phenotypes/ahmed-variants.tsv'
    if not os.path.exists(input_path): return
    print(f"Processing {input_path}...")
    df = pd.read_csv(input_path, sep='\t').rename(columns=lambda x: x.strip())
    records = []
    
    for _, row in tqdm(df.iterrows(), total=len(df)):
        rec = {k: '' for k in COLUMNS}
        rec['id'] = str(row.get('Case', ''))
        rec['source_file'] = 'ahmed-variants.tsv'
        rec['sex_id'], rec['sex_label'] = utils.normalize_sex(row.get('Gender'))
        rec['ancestry_id'], rec['ancestry_label'] = utils.normalize_ancestry()
        rec['country_id'], rec['country_label'] = utils.normalize_country()
        rec['consanguinity_id'], rec['consanguinity_label'] = utils.normalize_consanguinity(row.get('Consanguinity'))
        rec['family_history_id'], rec['family_history_label'] = utils.normalize_family_history(row.get('Family History'))
        rec['genotype_id'], rec['genotype_label'] = utils.normalize_genotype(row.get('Zygosity'))
        
        genes = set()
        vars_hgvs = []
        def add_var(g, v, p):
            if pd.isna(g): return
            gn = utils.map_gene(g)
            if gn:
                genes.add(gn)
                v_str = str(v).strip() if pd.notna(v) else ""
                if v_str.startswith('c.'): vars_hgvs.append(f"{gn[1]}:{v_str}")
                else: vars_hgvs.append(utils.extract_variant_hgvs(f"{g} {v} {p}", gn[1]))
        
        add_var(row.get('Gene'), row.get('Variant'), row.get('Protein'))
        add_var(row.get('Gene 2'), row.get('Variant 2'), row.get('Protein 2'))
        
        rec['gene_ids'] = "|".join([x[0] for x in genes])
        rec['gene_symbols'] = "|".join([x[1] for x in genes])
        rec['variant_hgvs'] = "|".join(vars_hgvs)
        rec['variant_classification'] = str(row.get('ACMG Classification', ''))
        
        # Phenotypes/Diseases via LLM
        raw_text = str(row.get('HPOs', ''))
        mapping = utils.llm_normalize([raw_text], context=rec['gene_symbols'])
        
        hpo_p, hpo_l, hpo_e, hpo_el = [], [], [], []
        omim_p, omim_l, omim_e, omim_el = [], [], [], []
        mondo_p, mondo_l, mondo_e, mondo_el = [], [], [], []
        
        for term, data in mapping.items():
            h = data.get('hpo', {})
            if h.get('id') and h['id'] != 'None':
                if h.get('excluded'): hpo_e.append(h['id']); hpo_el.append(h['label'])
                else: hpo_p.append(h['id']); hpo_l.append(h['label'])
            
            o = data.get('omim', {})
            if o.get('id') and o['id'] != 'None':
                if o.get('excluded'): omim_e.append(o['id']); omim_el.append(o['label'])
                else: omim_p.append(o['id']); omim_l.append(o['label'])

            m = data.get('mondo', {})
            if m.get('id') and m['id'] != 'None':
                if m.get('excluded'): mondo_e.append(m['id']); mondo_el.append(m['label'])
                else: mondo_p.append(m['id']); mondo_l.append(m['label'])

        rec['hpo_ids'], rec['hpo_labels'] = "|".join(hpo_p), "|".join(hpo_l)
        rec['hpo_excluded_ids'], rec['hpo_excluded_labels'] = "|".join(hpo_e), "|".join(hpo_el)
        rec['disease_omim_ids'], rec['disease_omim_labels'] = "|".join(omim_p), "|".join(omim_l)
        rec['disease_mondo_ids'], rec['disease_mondo_labels'] = "|".join(mondo_p), "|".join(mondo_l)
        
        records.append(rec)
    pd.DataFrame(records).to_csv(os.path.join(OUTPUT_DIR, 'ahmed_variants_normalized.tsv'), sep='\t', index=False)

if __name__ == "__main__":
    process_ahmed_variants()
