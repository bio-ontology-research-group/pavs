import pandas as pd
import numpy as np
import os
import json
import math

def parse_hpo_obo(obo_path):
    terms = {}
    with open(obo_path, 'r') as f:
        term = None
        for line in f:
            line = line.strip()
            if line == "[Term]":
                term = {'id': None, 'is_a': []}
            elif line.startswith("id: ") and term is not None:
                term['id'] = line[4:]
            elif line.startswith("is_a: ") and term is not None:
                parent = line[6:].split('!')[0].strip()
                term['is_a'].append(parent)
            elif line == "" and term is not None:
                if term['id']:
                    terms[term['id']] = term
                term = None
    return terms

def build_ancestors(terms):
    ancestors = {}
    def get_ancestors(term_id):
        if term_id in ancestors: return ancestors[term_id]
        res = {term_id}
        if term_id in terms:
            for parent in terms[term_id]['is_a']:
                res.update(get_ancestors(parent))
        ancestors[term_id] = res
        return res
    for t_id in terms: get_ancestors(t_id)
    return {k: frozenset(v) for k, v in ancestors.items()}

def compute_ic_extrinsic(terms, ancestors, gene_to_hpos):
    term_counts = {t_id: 0 for t_id in terms}
    total_genes = len(gene_to_hpos)
    for gene, hpos in gene_to_hpos.items():
        all_gene_anc = set()
        for h_id in hpos:
            if h_id in ancestors: all_gene_anc.update(ancestors[h_id])
        for anc in all_gene_anc:
            if anc in term_counts: term_counts[anc] += 1
    return {t_id: (-math.log(count / total_genes) if count > 0 else 0.0) for t_id, count in term_counts.items()}

def main():
    terms = parse_hpo_obo("ontology/hp.obo")
    ancestors = build_ancestors(terms)
    
    # Build Gene-to-HPO reference for Extrinsic IC
    gene_to_hpos = {}
    df_gp = pd.read_csv("data/reference/genes_to_phenotype.txt", sep='\t', comment='#')
    for _, row in df_gp.iterrows():
        gene, hpo_id = str(row['gene_symbol']), str(row['hpo_id'])
        if hpo_id in terms:
            if gene not in gene_to_hpos: gene_to_hpos[gene] = []
            gene_to_hpos[gene].append(hpo_id)
    gene_to_hpos = {k: set(v) for k, v in gene_to_hpos.items()}
    ic_extrinsic = compute_ic_extrinsic(terms, ancestors, gene_to_hpos)

    # Prepare cohorts
    cohort_data = []

    # 1. DDD and Saudi from TSV
    df_cases = pd.read_csv('data/combined_annotated-gpt4oss.tsv', sep='\t')
    for _, row in df_cases.iterrows():
        hpos = [h.strip() for h in str(row['phenotypes_present_ids']).replace('|', ',').split(',') if h.strip().startswith('HP:')]
        hpos = [h for h in hpos if h in terms]
        if hpos:
            group = 'DDD' if 'ddd' in str(row['source_file']).lower() else 'Saudi'
            avg_ic = np.mean([ic_extrinsic.get(h, 0.0) for h in hpos])
            cohort_data.append({'group': group, 'num_hpos': len(hpos), 'avg_ic': avg_ic})

    # 2. Literature from Phenopackets
    lit_dir = "phenopackets/0.1.26/"
    for root, _, files in os.walk(lit_dir):
        for file in files:
            if file.endswith(".json"):
                try:
                    with open(os.path.join(root, file)) as f:
                        data = json.load(f)
                        hpos = [feat['type']['id'] for feat in data.get('phenotypicFeatures', []) if feat.get('type', {}).get('id') in terms]
                        if hpos:
                            avg_ic = np.mean([ic_extrinsic.get(h, 0.0) for h in hpos])
                            cohort_data.append({'group': 'Literature', 'num_hpos': len(hpos), 'avg_ic': avg_ic})
                except: continue

    df_stats = pd.DataFrame(cohort_data)
    summary = df_stats.groupby('group').agg({
        'num_hpos': ['mean', 'std', 'max'],
        'avg_ic': ['mean', 'std']
    }).round(3)
    
    print(summary)
    summary.to_csv('analysis/phenotypic_density_ic_stats.csv')
    
    markdown_table = summary.to_markdown()
    
    with open('analysis/PHENOTYPIC_STATS.md', 'w') as f:
        f.write('# Phenotypic Density and Information Content Stats\n\n')
        f.write('This table compares the number of HPO terms per case and the average Extrinsic Information Content (IC) per annotation across cohorts.\n\n')
        f.write(markdown_table)
        f.write('\n\n**Note**: Average IC represents the mean specificity of the terms used. Higher IC indicates more specific/rare phenotypic descriptions.\n')

if __name__ == "__main__":
    main()
