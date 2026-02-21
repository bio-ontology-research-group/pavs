
import pandas as pd
import numpy as np
import os
import sys
import time
import math
import json
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from functools import lru_cache
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, auc

# --- HPO Parser ---

def parse_hpo_obo(obo_path):
    print(f"Parsing HPO from {obo_path}...")
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
    print("Building ancestor sets...")
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

def build_descendants(terms):
    print("Counting descendants...")
    full_desc = {}
    children = {t_id: [] for t_id in terms}
    for t_id, t_info in terms.items():
        for parent in t_info['is_a']:
            if parent in children: children[parent].append(t_id)
    def get_desc(t_id):
        if t_id in full_desc: return full_desc[t_id]
        res = {t_id}
        for child in children[t_id]: res.update(get_desc(child))
        full_desc[t_id] = res
        return res
    for t_id in terms: get_desc(t_id)
    return {t_id: len(res) for t_id, res in full_desc.items()}

def compute_ic(desc_counts, total):
    print("Computing IC...")
    return {t_id: -math.log(count / total) for t_id, count in desc_counts.items()}

# Global variables for workers
_ic = None
_ancestors = None
_gene_to_hpos = None
_all_genes = None

def init_worker(ic, ancestors, gene_to_hpos, all_genes):
    global _ic, _ancestors, _gene_to_hpos, _all_genes
    _ic = ic
    _ancestors = ancestors
    _gene_to_hpos = gene_to_hpos
    _all_genes = all_genes

@lru_cache(maxsize=100000)
def lin_sim(t1, t2):
    if t1 == t2: return 1.0
    if t1 not in _ancestors or t2 not in _ancestors: return 0.0
    common = _ancestors[t1] & _ancestors[t2]
    if not common: return 0.0
    max_ic = max(_ic.get(c, 0.0) for c in common)
    denom = _ic.get(t1, 0) + _ic.get(t2, 0)
    if denom == 0: return 0.0
    return (2.0 * max_ic) / denom

def bma_sim(set1, set2):
    if not set1 or not set2: return 0.0
    sum1 = 0
    for t1 in set1:
        max_s = 0
        for t2 in set2:
            s = lin_sim(t1, t2)
            if s > max_s: max_s = s
        sum1 += max_s
    sum2 = 0
    for t2 in set2:
        max_s = 0
        for t1 in set1:
            s = lin_sim(t1, t2)
            if s > max_s: max_s = s
        sum2 += max_s
    return (sum1 / len(set1) + sum2 / len(set2)) / 2.0

def process_case(args):
    pavs_id, case_hpos, correct_gene, group = args
    similarities = []
    for gene in _all_genes:
        sim = bma_sim(case_hpos, _gene_to_hpos[gene])
        similarities.append((gene, sim))
    similarities.sort(key=lambda x: (-x[1], x[0]))
    rank = -1
    score = 0
    for i, (gene, sim) in enumerate(similarities):
        if gene == correct_gene:
            rank = i + 1
            score = sim
            break
    return {'pavs_id': pavs_id, 'gene': correct_gene, 'rank': rank, 'score': score, 'group': group, 'num_candidates': len(similarities)}

def parse_literature_phenopackets(phenopackets_dir, terms):
    print(f"Parsing literature phenopackets from {phenopackets_dir}...")
    lit_cases = []
    for root, dirs, files in os.walk(phenopackets_dir):
        for file in files:
            if file.endswith(".json"):
                path = os.path.join(root, file)
                try:
                    with open(path, 'r') as f:
                        data = json.load(f)
                        pavs_id = data.get('id', file)
                        hpos = []
                        for feature in data.get('phenotypicFeatures', []):
                            h_id = feature.get('type', {}).get('id')
                            if h_id and h_id in terms:
                                hpos.append(h_id)
                        
                        # Find gene symbol from interpretations
                        gene_symbol = None
                        for interpretation in data.get('interpretations', []):
                            diagnosis = interpretation.get('diagnosis', {})
                            for gi in diagnosis.get('genomicInterpretations', []):
                                gene_symbol = gi.get('variantInterpretation', {}).get('variationDescriptor', {}).get('geneContext', {}).get('symbol')
                                if gene_symbol: break
                            if gene_symbol: break
                        
                        if hpos and gene_symbol:
                            lit_cases.append((pavs_id, tuple(sorted(list(set(hpos)))), gene_symbol, 'Literature'))
                except:
                    continue
    return lit_cases

def main():
    start_all = time.time()
    terms = parse_hpo_obo("ontology/hp.obo")
    ancestors = build_ancestors(terms)
    desc_counts = build_descendants(terms)
    ic = compute_ic(desc_counts, len(terms))
    
    print("Loading gene mappings...")
    gene_to_hpos = {}
    df_gp = pd.read_csv("data/genes_to_phenotype.txt", sep='\t', comment='#')
    for _, row in df_gp.iterrows():
        gene, hpo_id = str(row['gene_symbol']), str(row['hpo_id'])
        if hpo_id in terms:
            if gene not in gene_to_hpos: gene_to_hpos[gene] = []
            gene_to_hpos[gene].append(hpo_id)
    gene_to_hpos = {k: tuple(sorted(list(set(v)))) for k, v in gene_to_hpos.items()}
    all_genes = sorted(list(gene_to_hpos.keys()))
    
    tasks = []
    # 1. Saudi vs DDD cases
    print("Loading clinical cases from data/combined_annotated.tsv...")
    df_cases = pd.read_csv('data/combined_annotated.tsv', sep='\t')
    for _, row in df_cases.iterrows():
        gene = str(row['gene_symbol'])
        if gene in gene_to_hpos:
            hpos = [h.strip() for h in str(row['phenotypes_present_ids']).replace('|', ',').split(',') if h.strip().startswith('HP:')]
            hpos = tuple(sorted(list(set([h for h in hpos if h in terms]))))
            if hpos:
                group = 'DDD' if 'ddd' in str(row['source_file']).lower() else 'Saudi'
                tasks.append((row['pavs_id'], hpos, gene, group))
                
    # 2. Literature Phenopackets
    lit_cases = parse_literature_phenopackets("phenopackets/0.1.26/", terms)
    for l_id, l_hpos, l_gene, l_group in lit_cases:
        if l_gene in gene_to_hpos:
            tasks.append((l_id, l_hpos, l_gene, l_group))
            
    print(f"Total tasks to analyze: {len(tasks)}")
    
    num_procs = cpu_count()
    print(f"Analyzing with {num_procs} cores...")
    results = []
    with Pool(processes=num_procs, initializer=init_worker, initargs=(ic, ancestors, gene_to_hpos, all_genes)) as pool:
        for res in tqdm(pool.imap_unordered(process_case, tasks), total=len(tasks)):
            if res:
                results.append(res)
    
    results_df = pd.DataFrame(results)
    os.makedirs('analysis', exist_ok=True)
    results_df.to_csv('analysis/phenotype_similarity_full.csv', index=False)
    
    print(f"Total time: {(time.time() - start_all)/60:.2f} minutes.")
    generate_summary_and_plots(results_df)
    create_readme(results_df)

def generate_summary_and_plots(df):
    groups = df['group'].unique()
    lines = ["Phenotype Semantic Similarity Analysis Summary\n", "============================================\n\n"]
    
    for label in sorted(groups):
        g_df = df[df['group'] == label]
        if g_df.empty: continue
        ranks = g_df['rank'].values
        num_cand = g_df['num_candidates'].values
        aucs = (num_cand - ranks) / (num_cand - 1)
        mrr = np.mean(1.0 / ranks)
        
        lines.append(f"--- {label} (n={len(g_df)}) ---\n")
        lines.append(f"Mean Rank: {np.mean(ranks):.2f}\n")
        lines.append(f"Median Rank: {np.median(ranks):.2f}\n")
        lines.append(f"Hits@1: {np.mean(ranks <= 1):.4f}\n")
        lines.append(f"Hits@10: {np.mean(ranks <= 10):.4f}\n")
        lines.append(f"Hits@50: {np.mean(ranks <= 50):.4f}\n")
        lines.append(f"Mean AUC: {np.mean(aucs):.4f}\n")
        lines.append(f"Mean AUPR proxy (MRR): {mrr:.4f}\n\n")
        
        # ROC Plot
        plt.figure(figsize=(8, 6))
        sorted_ranks = np.sort(ranks)
        unique_ranks, counts = np.unique(sorted_ranks, return_counts=True)
        tpr = np.cumsum(counts) / len(ranks)
        fpr = (unique_ranks - 1) / (num_cand[0] - 1)
        fpr = np.insert(fpr, 0, 0); tpr = np.insert(tpr, 0, 0)
        plt.plot(fpr, tpr, label=f'AUC = {np.mean(aucs):.4f}')
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(f'ROC Curve - {label}')
        plt.legend(); plt.grid(True, linestyle=':', alpha=0.6)
        plt.savefig(f'analysis/roc_curve_{label.lower()}.png'); plt.close()

    with open('analysis/interpretation_full.txt', 'w') as f: f.writelines(lines)

def create_readme(df):
    with open('analysis/README.md', 'w') as f:
        f.write("# Phenotype Semantic Similarity Analysis\n\n")
        f.write("## Overview\n")
        f.write("This analysis compares Saudi clinical cases against the Deciphering Developmental Disorders (DDD) ")
        f.write("cohort and literature-curated Phenopackets (public store 0.1.26). The goal is to evaluate ")
        f.write("how well clinical phenotypes prioritize the causative gene using semantic similarity.\n\n")
        f.write("## Methodology\n")
        f.write("- **Similarity Metric**: Lin's semantic similarity with Best Match Average (BMA) aggregation.\n")
        f.write("- **Ontology**: Human Phenotype Ontology (HPO).\n")
        f.write("- **Reference**: HPO associations from `genes_to_phenotype.txt`.\n")
        f.write("- **Ranking**: Each case is ranked against all genes in the reference dataset (~5,200 genes).\n\n")
        f.write("## Cohorts\n")
        f.write("- **Saudi**: Cases from Saudi Arabian clinical studies (Alkuraya, Marwa, Fawzan, etc.).\n")
        f.write("- **DDD**: Cases from the DDD project as a baseline comparison.\n")
        f.write("- **Literature**: High-quality manually curated Phenopackets from the public store.\n\n")
        f.write("## Performance Metrics\n")
        f.write("- **AUC**: Area Under the ROC Curve, representing the probability that a true gene is ranked higher than a random one.\n")
        f.write("- **Hits@k**: Proportion of cases where the true gene is in the top k ranks.\n")
        f.write("- **MRR**: Mean Reciprocal Rank, a proxy for AUPR in ranking tasks.\n\n")
        f.write("## Results Summary\n")
        f.write("See `interpretation_full.txt` for detailed metrics and `roc_curve_*.png` for performance curves.\n")

if __name__ == "__main__":
    main()
