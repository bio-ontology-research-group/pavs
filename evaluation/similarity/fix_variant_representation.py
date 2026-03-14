import pandas as pd
import re
import os
import ast
import json
import requests

# --- Configuration ---
OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY")
MODEL_NAME = "google/gemini-2.0-flash-001"

# --- Zygosity Mappings ---
GENO_MAPPING = {
    'het': 'GENO:0000135', 'hom': 'GENO:0000136', 'hemi': 'GENO:0000134',
    'heterozygous': 'GENO:0000135', 'homozygous': 'GENO:0000136',
    'hemizygous': 'GENO:0000134', 'unknown': 'GENO:0000137',
    'comphet': 'GENO:0000402', 'compound heterozygous': 'GENO:0000402'
}

def load_genes(gene_info_path):
    gene_to_id = {}
    if not os.path.exists(gene_info_path): return gene_to_id
    df = pd.read_csv(gene_info_path, sep='\t')
    for _, row in df.iterrows():
        symbol = str(row['Symbol']).lower()
        gene_id = str(row['GeneID'])
        gene_to_id[symbol] = gene_id
        synonyms = str(row['Synonyms']).split('|')
        for s in synonyms:
            if s != '-': gene_to_id[s.lower()] = gene_id
    return gene_to_id

GENE_TO_ID = load_genes('data/Homo_sapiens.gene_info')

def get_entrez(symbol):
    if not symbol: return "Unknown"
    return GENE_TO_ID.get(str(symbol).strip().lower(), "Unknown")

def llm_parse_variants(text):
    if not OPENROUTER_API_KEY: return None
    prompt = f"""Analyze this clinical variant description. Extract all genes, variants (HGVS if possible), and the overall genotype (zygosity).
    Identify if it's a compound heterozygote (2 variants in 1 gene) or digenic (variants in 2+ genes).
    Text: "{text}"
    Return ONLY a JSON object:
    {{
        "variants": ["HGVS1", "HGVS2"],
        "genes": ["GENE1", "GENE2"],
        "overall_zygosity": "GENO:ID",
        "is_compound_het": true,
        "is_digenic": false
    }}
    Use GENO:0000402 for compound het, GENO:0000135 for simple het, GENO:0000136 for homo.
    """
    try:
        response = requests.post("https://openrouter.ai/api/v1/chat/completions",
            headers={"Authorization": f"Bearer {OPENROUTER_API_KEY}"},
            json={"model": MODEL_NAME, "messages": [{"role": "user", "content": prompt}], "response_format": {"type": "json_object"}}, timeout=15)
        return response.json()['choices'][0]['message']['content']
    except: return None

def fix_ahmed():
    file_path = 'data/cleaned_and_normalized/ahmed_normalized.tsv'
    if not os.path.exists(file_path): return
    df = pd.read_csv(file_path, sep='\t')
    for i, row in df.iterrows():
        genes, variants = [], []
        for suffix in ['', ' 2', ' 3', ' 4']:
            g_col, v_col = f'Gene{suffix}'.strip(), f'Variant{suffix}'.strip()
            if g_col in df.columns and v_col in df.columns:
                g, v = str(row[g_col]).strip(), str(row[v_col]).strip()
                if g and g != 'nan' and v and v != 'nan':
                    genes.append(g); variants.append(f"c.{v}" if not v.startswith('c.') else v)
        if len(variants) > 1:
            unique_genes = list(set(genes))
            zyg = 'GENO:0000402' if len(unique_genes) == 1 else 'GENO:0000135'
            df.at[i, 'normalized_variants_hgvs'] = str(variants)
            df.at[i, 'normalized_genes_entrez'] = str([get_entrez(g) for g in unique_genes])
            df.at[i, 'normalized_zygosity_geno'] = zyg
    df.to_csv(file_path, sep='\t', index=False)
    print("Fixed ahmed_normalized.tsv")

def fix_fawzan():
    file_path = 'data/cleaned_and_normalized/fawzan_normalized.tsv'
    if not os.path.exists(file_path): return
    df = pd.read_csv(file_path, sep='\t')
    for i, row in df.iterrows():
        try: vars_list = ast.literal_eval(str(row['normalized_variants_hgvs']))
        except: vars_list = []
        if len(vars_list) > 1:
            genes = []
            for v in vars_list:
                m = re.search(r'\((.*?)\)', v)
                if m: genes.append(m.group(1))
            if len(set(genes)) == 1: df.at[i, 'normalized_zygosity_geno'] = str(['GENO:0000402'])
    df.to_csv(file_path, sep='\t', index=False)
    print("Fixed fawzan_normalized.tsv")

def fix_marwa():
    file_path = 'data/cleaned_and_normalized/marwa_normalized.tsv'
    if not os.path.exists(file_path): return
    df = pd.read_csv(file_path, sep='\t')
    for i, row in df.iterrows():
        v_text = str(row['variants'])
        if ';' in v_text:
            res_str = llm_parse_variants(v_text)
            if res_str:
                try:
                    data = json.loads(res_str)
                    df.at[i, 'normalized_variants_hgvs'] = str(data.get('variants', []))
                    g_list = data.get('genes', [])
                    df.at[i, 'norm_genes'] = str(g_list)
                    df.at[i, 'normalized_genes_entrez'] = str([get_entrez(g) for g in g_list])
                    df.at[i, 'norm_zygosity'] = data.get('overall_zygosity', 'GENO:0000137')
                except: continue
    df.to_csv(file_path, sep='\t', index=False)
    print("Fixed marwa_normalized.tsv")

if __name__ == "__main__":
    fix_ahmed()
    fix_fawzan()
    fix_marwa()