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

def process_marwa():
    input_path = 'data/phenotypes/marwa-variants.tsv'
    if not os.path.exists(input_path): return
    print(f"Processing {input_path}...")
    
    df = pd.read_csv(input_path, sep='\t')
    df.columns = df.columns.str.strip()
    records = []
    
    for _, row in tqdm(df.iterrows(), total=len(df)):
        rec = {k: '' for k in COLUMNS}
        rec['id'] = clean_text(row['ID'])
        rec['source_file'] = 'marwa-variants.tsv'
        
        # Sex: Not provided in Marwa file head
        rec['sex_id'], rec['sex_label'] = utils.normalize_sex('Unknown')
        
        # Age: Not provided
        
        # Ancestry/Country
        aid, alabel = utils.normalize_ancestry()
        rec['ancestry_id'] = aid
        rec['ancestry_label'] = alabel
        cid, clabel = utils.normalize_country()
        rec['country_id'] = cid
        rec['country_label'] = clabel
        
        # Consanguinity: Not explicitly provided
        
        # Variants: "variants" column
        # e.g. "CDHR1:NM_033100.3:c.476C>A; p.Ala159Glu"
        # e.g. "chr2:112,309,575–112,815,205 Hom Del of ex8" -> Needs LLM
        v_str = clean_text(row['variants'])
        
        genes_found = set()
        variants_hgvs = []
        genotype_hint = ""
        
        if v_str:
            # Try splitting by colon
            parts = v_str.split(':')
            if len(parts) >= 3 and (parts[2].startswith('c.') or 'NM_' in v_str):
                gene = parts[0]
                g_norm = utils.map_gene(gene)
                g_sym = g_norm[1] if g_norm else gene.upper()
                g_id = g_norm[0] if g_norm else ''
                if g_id: genes_found.add((g_id, g_sym))
                
                # Extract c.
                m_c = re.search(r'(c\.[^;\s]+)', v_str)
                if m_c:
                    variants_hgvs.append(f"{g_sym}:{m_c.group(1)}")
                else:
                    variants_hgvs.append(f"{g_sym}:{v_str}")
            else:
                # LLM
                extracted = utils.extract_variant_hgvs(v_str)
                if extracted:
                    # Try to parse gene from extracted
                    m_g = re.match(r'^([^:]+):', extracted)
                    if m_g:
                        gene = m_g.group(1)
                        g_norm = utils.map_gene(gene)
                        if g_norm: genes_found.add(g_norm)
                    variants_hgvs.append(extracted)
            
            # Infer genotype from text if possible ("Hom", "Het")
            if 'Hom' in v_str: genotype_hint = 'homozygous'
            elif 'Het' in v_str: genotype_hint = 'heterozygous'
            
        rec['gene_ids'] = "|".join([x[0] for x in genes_found])
        rec['gene_symbols'] = "|".join([x[1] for x in genes_found])
        rec['variant_hgvs'] = "|".join(variants_hgvs)
        
        if genotype_hint:
            gid, glabel = utils.normalize_genotype(genotype_hint)
            rec['genotype_id'] = gid
            rec['genotype_label'] = glabel
        
        # Phenotypes: "phenotypes" column
        # e.g. "macular focal atrophy (HP:0007401)| mild intraretinal pigment (HP:0007703)..."
        raw_hpo = clean_text(row['phenotypes'])
        
        # Extract explicit IDs
        explicit_ids = re.findall(r'HP:\d{7}', raw_hpo)
        explicit_omim = re.findall(r'OMIM:\d+', raw_hpo)
        
        clean_pheno = re.sub(r'\(?HP:\d{7}\)?', '', raw_hpo)
        clean_pheno = re.sub(r'\(?OMIM:\d+\)?', '', clean_pheno)
        
        present_hpos = set(explicit_ids)
        present_diseases = set(explicit_omim)
        excluded_hpos = set()
        
        if clean_pheno:
            llm_res_str = utils.llm_map_terms([clean_pheno], context_str=f"Gene: {rec['gene_symbols']}")
            try:
                mapping = json.loads(llm_res_str)
                for term_key, data in mapping.items():
                    if data.get('excluded'):
                         if data.get('hpo_id') and data['hpo_id'] != 'None':
                            excluded_hpos.add(data['hpo_id'])
                    else:
                        if data.get('hpo_id') and data['hpo_id'] != 'None':
                            present_hpos.add(data['hpo_id'])
                        if data.get('disease_id'):
                            present_diseases.add(data['disease_id'])
            except: pass
            
        rec['hpo_ids'] = "|".join(sorted(list(present_hpos)))
        labels = []
        for h in sorted(list(present_hpos)):
            match = next((t for t in utils.hpo_terms if t['id'] == h), None)
            labels.append(match['name'] if match else h)
        rec['hpo_labels'] = "|".join(labels)
        
        rec['hpo_excluded_ids'] = "|".join(sorted(list(excluded_hpos)))
        rec['disease_omim_ids'] = "|".join(sorted(list(present_diseases)))
        
        records.append(rec)
        
    pd.DataFrame(records).to_csv(os.path.join(OUTPUT_DIR, 'marwa_normalized.tsv'), sep='\t', index=False)
    print("Finished marwa.")

if __name__ == "__main__":
    process_marwa()
