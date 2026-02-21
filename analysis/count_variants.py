import pandas as pd
import ast
import os

def count_stats(file_path):
    if not os.path.exists(file_path):
        return None
    
    df = pd.read_csv(file_path, sep='	')
    
    comp_het_count = 0
    digenic_count = 0
    
    # Identify relevant columns
    gene_col = None
    zyg_col = None
    
    for col in df.columns:
        if 'genes_entrez' in col: gene_col = col
        if 'zygosity_geno' in col: zyg_col = col
    
    # Fallback for marwa which might have different names
    if not gene_col and 'norm_genes' in df.columns: gene_col = 'normalized_genes_entrez'
    if not zyg_col and 'norm_zygosity' in df.columns: zyg_col = 'norm_zygosity'

    for _, row in df.iterrows():
        # Digenic check: > 1 unique gene
        genes = []
        g_val = str(row.get(gene_col, '[]'))
        try:
            if g_val.startswith('['):
                genes = ast.literal_eval(g_val)
            else:
                genes = [g_val] if g_val != 'nan' else []
        except:
            genes = []
        
        # Filter Unknowns and duplicates
        genes = [g for g in genes if g and str(g).lower() != 'unknown']
        if len(set(genes)) > 1:
            digenic_count += 1
            
        # Comp Het check: GENO:0000402
        z_val = str(row.get(zyg_col, ''))
        if 'GENO:0000402' in z_val:
            comp_het_count += 1
            
    return {'comp_het': comp_het_count, 'digenic': digenic_count, 'total': len(df)}

files = [
    'data/cleaned_and_normalized/ahmed_normalized.tsv',
    'data/cleaned_and_normalized/fawzan_normalized.tsv',
    'data/cleaned_and_normalized/marwa_normalized.tsv',
    'data/cleaned_and_normalized/ahmed_pmid28454995_normalized.tsv',
    'data/cleaned_and_normalized/pmc7082194_normalized.tsv'
]

print(f"{'File':<40} | {'Comp Het':<10} | {'Digenic':<10} | {'Total'}")
print("-" * 80)
for f in files:
    stats = count_stats(f)
    if stats:
        name = os.path.basename(f)
        print(f"{name:<40} | {stats['comp_het']:<10} | {stats['digenic']:<10} | {stats['total']}")
