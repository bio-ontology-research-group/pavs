import pandas as pd
import json
import re
import os

def load_mappings():
    if os.path.exists('hpo_mappings_final.json'):
        with open('hpo_mappings_final.json', 'r') as f:
            return json.load(f)
    return {}

def load_patient_mappings():
    if os.path.exists('patient_hpo_mappings.json'):
        with open('patient_hpo_mappings.json', 'r') as f:
            return json.load(f)
    return {}

def load_valid_genes():
    genes = set()
    path = 'data/genes_to_phenotype.txt'
    if os.path.exists(path):
        try:
            df = pd.read_csv(path, sep='\t')
            if 'gene_symbol' in df.columns:
                genes = set(df['gene_symbol'].unique())
            else:
                df = pd.read_csv(path, sep='\t', header=None)
                genes = set(df[1].unique())
        except: pass
    return genes

def load_valid_omims():
    omims = set()
    path = 'data/phenotype.hpoa'
    if os.path.exists(path):
        try:
            with open(path, 'r') as f:
                for line in f:
                    if line.startswith('OMIM:'):
                        omims.add(line.split('\t')[0])
        except: pass
    return omims

def clean_variant(v, valid_genes):
    if not isinstance(v, str) or v == 'nan' or not v: return ""
    parts = v.split('|')
    cleaned_parts = []
    for p in parts:
        p = p.strip()
        if not p: continue
        gene_match = re.search(r'^([^:]+)', p)
        gene_symbol = "not known"
        rest = p
        if gene_match:
            potential_gene = gene_match.group(1).strip()
            if potential_gene in valid_genes:
                gene_symbol = potential_gene
                rest = p[len(potential_gene):].lstrip(':')
        nm = re.search(r'(NM_\d+(\.\d+)?)', rest)
        c_dna = re.search(r'(c\.[^:\s,;]+)', rest)
        if nm and c_dna: variant_part = f"{nm.group(1)}:{c_dna.group(1)}"
        elif c_dna: variant_part = c_dna.group(1)
        else: variant_part = rest.split(' ')[0].split('(')[0].rstrip(':')
        cleaned_parts.append(f"{gene_symbol}:{variant_part}")
    return "|".join(cleaned_parts)

def main():
    df = pd.read_csv('master_data_v1.tsv', sep='\t')
    global_mappings = load_mappings()
    patient_mappings = load_patient_mappings()
    valid_genes = load_valid_genes()
    valid_omims = load_valid_omims()
    
    final_records = []
    
    for _, row in df.iterrows():
        p_id = str(row['Master_ID'])
        
        # High confidence HPOs from initial parsing
        hpos = []
        if pd.notna(row['Phenotypes_HPO']) and row['Phenotypes_HPO'] != "":
            for h in str(row['Phenotypes_HPO']).split('|'):
                if h.startswith('HP:'):
                    hpos.append({'id': h, 'severity': 'None', 'excluded': False})
        
        # Unmatched terms resolution
        unmatched = str(row['Unmatched_Phenotypes']).split('|') if pd.notna(row['Unmatched_Phenotypes']) and row['Unmatched_Phenotypes'] != "" else []
        p_map = patient_mappings.get(p_id, {})
        
        for term in unmatched:
            res = p_map.get(term)
            if res:
                if isinstance(res, dict):
                    if res.get('hpo_id') and res['hpo_id'] != "None":
                        hpos.append({
                            'id': res['hpo_id'],
                            'severity': res.get('severity_id', 'None'),
                            'excluded': res.get('excluded', False)
                        })
                elif isinstance(res, str) and res.startswith('HP:'):
                    hpos.append({'id': res, 'severity': 'None', 'excluded': False})
            else:
                # Fallback to global
                g_id = global_mappings.get(term)
                if g_id and g_id != "None" and g_id.startswith('HP:'):
                    hpos.append({'id': g_id, 'severity': 'None', 'excluded': False})
        
        # Deduplicate and separate present vs excluded
        present = []
        excluded = []
        for h in hpos:
            entry = h['id']
            if h['severity'] and h['severity'] != "None":
                entry += f"({h['severity']})"
            
            if h['excluded']:
                excluded.append(entry)
            else:
                present.append(entry)
        
        final_records.append({
            'Master_ID': p_id,
            'Source': row['Source'],
            'Original_ID': row['Original_ID'],
            'Sex': row['Sex'],
            'Age': row['Age'],
            'Zygosity_GENO': row['Zygosity'], # Keep for now
            'Phenotypes_Present_HPO': "|".join(sorted(list(set(present)))),
            'Phenotypes_Excluded_HPO': "|".join(sorted(list(set(excluded)))),
            'Variants_HGVS': clean_variant(row['Variants'], valid_genes),
            'Genome_Build': 'GRCh37',
            'Diseases': row['Diseases']
        })
    
    final_df = pd.DataFrame(final_records)
    final_df.to_csv('master_data_final.tsv', sep='\t', index=False)
    print(f"Final master data saved with {len(final_df)} records.")

if __name__ == "__main__":
    main()
