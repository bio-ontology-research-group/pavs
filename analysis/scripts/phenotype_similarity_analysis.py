
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
    return {t_id: res for t_id, res in full_desc.items()}

def compute_ic_intrinsic(terms, desc_map):
    print("Computing Intrinsic IC...")
    total = len(terms)
    return {t_id: -math.log(len(desc_set) / total) for t_id, desc_set in desc_map.items()}

def compute_ic_extrinsic(terms, ancestors, gene_to_hpos):
    print("Computing Extrinsic IC (over genes)...")
    term_counts = {t_id: 0 for t_id in terms}
    total_genes = len(gene_to_hpos)
    for gene, hpos in gene_to_hpos.items():
        all_gene_anc = set()
        for h_id in hpos:
            if h_id in ancestors: all_gene_anc.update(ancestors[h_id])
        for anc in all_gene_anc:
            if anc in term_counts: term_counts[anc] += 1
    return {t_id: (-math.log(count / total_genes) if count > 0 else 0.0) for t_id, count in term_counts.items()}

# Global variables for workers
_ic = None
_ancestors = None
_gene_to_hpos = None
_all_genes = None
_method = 'lin'

def init_worker(ic, ancestors, gene_to_hpos, all_genes, method):
    global _ic, _ancestors, _gene_to_hpos, _all_genes, _method
    _ic = ic
    _ancestors = ancestors
    _gene_to_hpos = gene_to_hpos
    _all_genes = all_genes
    _method = method

@lru_cache(maxsize=200000)
def pair_sim(t1, t2):
    if t1 not in _ancestors or t2 not in _ancestors: return 0.0
    common = _ancestors[t1] & _ancestors[t2]
    if not common: return 0.0
    max_ic = max(_ic.get(c, 0.0) for c in common)
    if _method == 'resnik':
        return max_ic
    # Lin
    denom = _ic.get(t1, 0) + _ic.get(t2, 0)
    return (2.0 * max_ic) / denom if denom > 0 else 0.0

def bma_sim(set1, set2):
    if not set1 or not set2: return 0.0
    sum1 = 0
    for t1 in set1:
        ms = 0
        for t2 in set2:
            s = pair_sim(t1, t2)
            if s > ms: ms = s
        sum1 += ms
    sum2 = 0
    for t2 in set2:
        ms = 0
        for t1 in set1:
            s = pair_sim(t1, t2)
            if s > ms: ms = s
        sum2 += ms
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

def run_analysis(tasks, ic, ancestors, gene_to_hpos, all_genes, method, ic_type, num_cores):
    print(f"Starting {ic_type} {method} analysis...")
    pair_sim.cache_clear()
    results = []
    with Pool(processes=num_cores, initializer=init_worker, initargs=(ic, ancestors, gene_to_hpos, all_genes, method)) as pool:
        for res in tqdm(pool.imap_unordered(process_case, tasks), total=len(tasks), desc=f"{ic_type}-{method}"):
            if res: results.append(res)
    df = pd.DataFrame(results)
    df['ic_type'] = ic_type
    df['sim_method'] = method
    return df

def main():
    num_cores = cpu_count()
    terms = parse_hpo_obo("ontology/hp.obo")
    ancestors = build_ancestors(terms)
    desc_map = build_descendants(terms)
    
    print("Loading gene mappings...")
    gene_to_hpos = {}
    df_gp = pd.read_csv("data/reference/genes_to_phenotype.txt", sep='\t', comment='#')
    for _, row in df_gp.iterrows():
        gene, hpo_id = str(row['gene_symbol']), str(row['hpo_id'])
        if hpo_id in terms:
            if gene not in gene_to_hpos: gene_to_hpos[gene] = []
            gene_to_hpos[gene].append(hpo_id)
    gene_to_hpos = {k: tuple(sorted(list(set(v)))) for k, v in gene_to_hpos.items()}
    all_genes = sorted(list(gene_to_hpos.keys()))
    
    ic_intrinsic = compute_ic_intrinsic(terms, desc_map)
    ic_extrinsic = compute_ic_extrinsic(terms, ancestors, gene_to_hpos)
    
    tasks = []
    df_cases = pd.read_csv('data/combined_annotated.tsv', sep='\t')
    for _, row in df_cases.iterrows():
        gene = str(row['gene_symbol'])
        if gene in gene_to_hpos:
            hpos = [h.strip() for h in str(row['phenotypes_present_ids']).replace('|', ',').split(',') if h.strip().startswith('HP:')]
            hpos = tuple(sorted(list(set([h for h in hpos if h in terms]))))
            if hpos:
                group = 'DDD' if 'ddd' in str(row['source_file']).lower() else 'Saudi'
                tasks.append((row['pavs_id'], hpos, gene, group))
                
    lit_dir = "phenopackets/0.1.26/"
    for root, _, files in os.walk(lit_dir):
        for file in files:
            if file.endswith(".json"):
                try:
                    with open(os.path.join(root, file)) as f:
                        data = json.load(f)
                        hpos = [f['type']['id'] for f in data.get('phenotypicFeatures', []) if f.get('type', {}).get('id') in terms]
                        g_sym = None
                        for interp in data.get('interpretations', []):
                            for gi in interp.get('diagnosis', {}).get('genomicInterpretations', []):
                                g_sym = gi.get('variantInterpretation', {}).get('variationDescriptor', {}).get('geneContext', {}).get('symbol')
                                if g_sym: break
                            if g_sym: break
                        if hpos and g_sym and g_sym in gene_to_hpos:
                            tasks.append((data.get('id', file), tuple(sorted(list(set(hpos)))), g_sym, 'Literature'))
                except: continue

    all_results = []
    # 4 combinations
    for ic_type, ic_vals in [('intrinsic', ic_intrinsic), ('extrinsic', ic_extrinsic)]:
        for sim_method in ['lin', 'resnik']:
            res = run_analysis(tasks, ic_vals, ancestors, gene_to_hpos, all_genes, sim_method, ic_type, num_cores)
            res.to_csv(f'analysis/results_{ic_type}_{sim_method}.csv', index=False)
            all_results.append(res)
    
    final_df = pd.concat(all_results)
    final_df.to_csv('analysis/phenotype_similarity_comprehensive.csv', index=False)
    print("Comprehensive results saved to analysis/phenotype_similarity_comprehensive.csv")

if __name__ == "__main__":
    main()
