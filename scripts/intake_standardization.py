import pandas as pd
import re
import json
import os
import ast

# Ontology Mappings
GENO_MAPPING = {
    'het': 'GENO:0000135',
    'heterozygous': 'GENO:0000135',
    'homo': 'GENO:0000136',
    'homozygous': 'GENO:0000136',
    'hemi': 'GENO:0000134',
    'hemizygous': 'GENO:0000134',
    'unknown': 'GENO:0000137'
}

GENO_LABELS = {
    'GENO:0000135': 'heterozygous',
    'GENO:0000136': 'homozygous',
    'GENO:0000134': 'hemizygous',
    'GENO:0000137': 'unknown zygosity'
}

SEX_MAPPING = {
    'm': 'MALE',
    'male': 'MALE',
    'f': 'FEMALE',
    'female': 'FEMALE',
    '?': 'UNKNOWN_SEX',
    'unknown': 'UNKNOWN_SEX'
}

DEFAULT_ANCESTRY = "HANCESTRO:0852" # Middle Eastern
DEFAULT_ANCESTRY_LABEL = "Middle Eastern population"
DEFAULT_COUNTRY = "GAZ:00005279"    # Saudi Arabia
DEFAULT_COUNTRY_LABEL = "Saudi Arabia"

def load_hpo_labels(obo_path):
    labels = {}
    if not os.path.exists(obo_path): return labels
    with open(obo_path, 'r') as f:
        current_id = None
        for line in f:
            if line.startswith('id: HP:'):
                current_id = line.split('id: ')[1].strip()
            elif line.startswith('name: ') and current_id:
                labels[current_id] = line.split('name: ')[1].strip()
                current_id = None
    return labels

def load_disease_labels(hpoa_path):
    labels = {}
    if not os.path.exists(hpoa_path): return labels
    # database_id     disease_name
    df = pd.read_csv(hpoa_path, sep='\t', comment='#')
    for _, row in df.iterrows():
        labels[row['database_id']] = row['disease_name']
    return labels

HPO_LABELS = load_hpo_labels('ontology/hp.obo')
DISEASE_LABELS = load_disease_labels('data/phenotype.hpoa')

def get_hpo_labels(ids):
    return [HPO_LABELS.get(i, "Unknown Phenotype") for i in ids]

def get_disease_labels(ids):
    return [DISEASE_LABELS.get(i, "Unknown Disease") for i in ids]

def clean_hpo_ids(text):
    if pd.isna(text): return []
    return sorted(list(set(re.findall(r'HP:(\d{7})', str(text)))))

def format_hpo(ids):
    return [f"HP:{i}" for i in ids]

def clean_omim_ids(text):
    if pd.isna(text): return []
    explicit = re.findall(r'OMIM:(\d{6})', str(text), re.I)
    return sorted(list(set([f"OMIM:{i}" for i in explicit])))

def clean_orpha_ids(text):
    if pd.isna(text): return []
    return sorted(list(set([f"ORPHA:{i}" for i in re.findall(r'ORPHA:?\s*(\d+)', str(text), re.I)])))

def normalize_zygosity(text):
    if pd.isna(text): return GENO_MAPPING['unknown']
    val = str(text).lower().strip()
    return GENO_MAPPING.get(val, GENO_MAPPING['unknown'])

def normalize_sex(text):
    if pd.isna(text): return SEX_MAPPING['unknown']
    val = str(text).lower().strip()
    return SEX_MAPPING.get(val, SEX_MAPPING['unknown'])

def extract_gene_hgvs(variant_str):
    if pd.isna(variant_str) or variant_str == '0' or variant_str == '':
        return None, None
    variant_str = str(variant_str).strip()
    if re.match(r'^(NC_|NM_|NP_|NG_|NR_|LRG_)', variant_str):
        return None, variant_str
    m1 = re.match(r'([^:]+):([^:]+):(?:exon\d+:)?(c\.[^:]+)(?::(p\.[^:]+))?', variant_str)
    if m1: return m1.group(1), m1.group(3)
    m2 = re.match(r'([^:]+):([^:]+):(c\.[^; ]+)', variant_str)
    if m2: return m2.group(1), m2.group(3)
    m3 = re.match(r'([^:]+):.*(c\.[^: (]+)', variant_str)
    if m3: return m3.group(1), m3.group(2)
    if variant_str.startswith('c.') or variant_str.startswith('g.'):
        return None, variant_str
    return None, variant_str

def split_phenotypes(text, delimiter=','):
    if pd.isna(text): return [], ""
    clean_text = re.sub(r'\(?HP:\d{7}\)?', '', str(text))
    original_terms = [t.strip() for t in clean_text.split(delimiter) if t.strip()]
    raw_terms = [t.strip() for t in str(text).split(delimiter) if t.strip()]
    justification = f"Original semantic terms (approx): {len(original_terms)}."
    return raw_terms, justification

def process_fawzan_file(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df.columns = [c.strip() for c in df.columns]
    results = []
    for _, row in df.iterrows():
        gene1, hgvs1 = extract_gene_hgvs(row['Variant(s)'])
        gene2, hgvs2 = extract_gene_hgvs(row['Variant(s).1'])
        genes = list(set([g for g in [gene1, gene2] if g]))
        variants = list(set([v for v in [hgvs1, hgvs2] if v]))
        hpo_ids = format_hpo(clean_hpo_ids(row['Phenotype']))
        raw_terms, just = split_phenotypes(row['Phenotype'], delimiter=',')
        low_conf = [t for t in raw_terms if not re.search(r'HP:\d{7}', t)]
        zyg = normalize_zygosity(row['Zyogsity'])
        disease_ids = sorted(list(set(clean_omim_ids(row['Phenotype']) + clean_orpha_ids(row['Phenotype']))))
        
        entry = row.to_dict() # Keep original columns
        entry.update({
            'norm_original_id': str(row['ID']),
            'norm_sex': normalize_sex(row['Gender']),
            'norm_phenotypes_high_conf': hpo_ids,
            'norm_phenotypes_labels': get_hpo_labels(hpo_ids),
            'norm_phenotypes_low_conf': low_conf,
            'norm_genes': genes,
            'norm_variants_hgvs': variants,
            'norm_zygosity': zyg,
            'norm_zygosity_label': GENO_LABELS.get(zyg, 'unknown'),
            'norm_diseases_high_conf': disease_ids,
            'norm_diseases_labels': get_disease_labels(disease_ids),
            'norm_ancestry': DEFAULT_ANCESTRY,
            'norm_ancestry_label': DEFAULT_ANCESTRY_LABEL,
            'norm_country': DEFAULT_COUNTRY,
            'norm_country_label': DEFAULT_COUNTRY_LABEL
        })
        results.append(entry)
    return results

def process_ahmed_file(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df.columns = [c.strip() for c in df.columns]
    results = []
    for _, row in df.iterrows():
        gene = str(row['Gene']).strip() if not pd.isna(row['Gene']) else None
        hgvs = str(row['Variant']).strip() if not pd.isna(row['Variant']) else None
        if hgvs and not hgvs.startswith('c.') and not hgvs.startswith('g.') and not hgvs.startswith('NC_'):
            hgvs = f"c.{hgvs}"
        hpo_ids = format_hpo(clean_hpo_ids(row['HPOs']))
        raw_terms = [t.strip() for t in re.split(r',', str(row['HPOs'])) if t.strip()]
        low_conf = [t for t in raw_terms if not re.search(r'HP:\d{7}', t)]
        omim_ids = []
        if not pd.isna(row['Omim']):
            o_val = str(row['Omim']).strip()
            digits = re.findall(r'\b\d{6}\b', o_val)
            omim_ids = [f"OMIM:{d}" for d in digits]
        orpha_ids = clean_orpha_ids(row['HPOs'])
        disease_ids = sorted(list(set(omim_ids + orpha_ids)))
        zyg = normalize_zygosity(row['Zygosity'])
        
        entry = row.to_dict()
        entry.update({
            'norm_original_id': str(row['Case']),
            'norm_sex': normalize_sex(row['Gender']),
            'norm_phenotypes_high_conf': hpo_ids,
            'norm_phenotypes_labels': get_hpo_labels(hpo_ids),
            'norm_phenotypes_low_conf': low_conf,
            'norm_genes': [gene] if gene else [],
            'norm_variants_hgvs': [hgvs] if hgvs else [],
            'norm_zygosity': zyg,
            'norm_zygosity_label': GENO_LABELS.get(zyg, 'unknown'),
            'norm_diseases_high_conf': disease_ids,
            'norm_diseases_labels': get_disease_labels(disease_ids),
            'norm_ancestry': DEFAULT_ANCESTRY,
            'norm_ancestry_label': DEFAULT_ANCESTRY_LABEL,
            'norm_country': DEFAULT_COUNTRY,
            'norm_country_label': DEFAULT_COUNTRY_LABEL
        })
        results.append(entry)
    return results

def process_marwa_file(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df.columns = [c.strip() for c in df.columns]
    results = []
    for _, row in df.iterrows():
        hpo_ids = format_hpo(clean_hpo_ids(row['phenotypes']))
        raw_terms = [t.strip() for t in str(row['phenotypes']).split('|') if t.strip()]
        low_conf = [t for t in raw_terms if not re.search(r'HP:\d{7}', t)]
        gene, hgvs = extract_gene_hgvs(row['variants'])
        zygosity = GENO_MAPPING['unknown']
        if not pd.isna(row['variants']):
            v_lower = str(row['variants']).lower()
            if 'hom' in v_lower: zygosity = GENO_MAPPING['homozygous']
            elif 'het' in v_lower: zygosity = GENO_MAPPING['heterozygous']
            elif 'hemi' in v_lower: zygosity = GENO_MAPPING['hemizygous']
        omim_ids = clean_omim_ids(row['phenotypes'])
        orpha_ids = clean_orpha_ids(row['phenotypes'])
        disease_ids = sorted(list(set(omim_ids + orpha_ids)))
        
        entry = row.to_dict()
        entry.update({
            'norm_original_id': str(row['ID']),
            'norm_phenotypes_high_conf': hpo_ids,
            'norm_phenotypes_labels': get_hpo_labels(hpo_ids),
            'norm_phenotypes_low_conf': low_conf,
            'norm_genes': [gene] if gene else [],
            'norm_variants_hgvs': [hgvs] if hgvs else [],
            'norm_zygosity': zygosity,
            'norm_zygosity_label': GENO_LABELS.get(zygosity, 'unknown'),
            'norm_diseases_high_conf': disease_ids,
            'norm_diseases_labels': get_disease_labels(disease_ids),
            'norm_ancestry': DEFAULT_ANCESTRY,
            'norm_ancestry_label': DEFAULT_ANCESTRY_LABEL,
            'norm_country': DEFAULT_COUNTRY,
            'norm_country_label': DEFAULT_COUNTRY_LABEL
        })
        results.append(entry)
    return results

def main():
    paths = {'fawzan': 'data/phenotypes/fawzan-variants.tsv', 'ahmed': 'data/phenotypes/ahmed-variants.tsv', 'marwa': 'data/phenotypes/marwa-variants.tsv'}
    os.makedirs('data/cleaned', exist_ok=True)
    
    for key, path in paths.items():
        if os.path.exists(path):
            if key == 'fawzan': res = process_fawzan_file(path)
            elif key == 'ahmed': res = process_ahmed_file(path)
            elif key == 'marwa': res = process_marwa_file(path)
            pd.DataFrame(res).to_csv(f'data/cleaned/{key}_cleaned.tsv', sep='\t', index=False)
            print(f"Processed {key} file: {len(res)} entries")

if __name__ == "__main__":
    main()
