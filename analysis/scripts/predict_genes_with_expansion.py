#!/usr/bin/env python3
"""
Gene prediction with HPO expansion (Resnik + BMA)
"""

import pandas as pd
import pronto
from collections import defaultdict
import numpy as np
from tqdm import tqdm
import pickle
import os
import json
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve
import matplotlib.pyplot as plt
import re
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')


# ==================== HPO Expansion ====================

def expand_hpo_with_ancestors(hpo_terms, ont, max_levels=6):
    """
    Expand HPO terms by adding ancestors up to max_levels
    
    Args:
        hpo_terms: List of patient HPO IDs
        ont: Pronto ontology object
        max_levels: How many ancestor levels to include (default: 3)
    
    Returns:
        Expanded list of HPO terms (original + ancestors)
    """
    expanded = set(hpo_terms)
    
    for term_id in hpo_terms:
        try:
            term = ont[term_id]
            # Get ancestors up to max_levels
            level = 0
            for ancestor in term.superclasses():
                if level >= max_levels:
                    break
                if ancestor.id.startswith('HP:'):
                    expanded.add(ancestor.id)
                    level += 1
        except:
            # If term not found, just keep the original
            pass
    
    return list(expanded)


# ==================== Data Loading ====================

def load_hpo_ontology(hpo_file='resources/hp.obo'):
    """Load HPO ontology"""
    print(f"Loading HPO ontology from {hpo_file}...")
    ont = pronto.Ontology(hpo_file)
    print(f"Loaded HPO ontology")
    return ont


def load_hpo_gene_annotations(annotation_file='resources/phenotype_to_genes.txt'):
    """Load HPO-gene annotations"""
    print(f"Loading HPO-gene annotations from {annotation_file}...")
    
    gene_to_hpo = defaultdict(set)
    hpo_to_genes = defaultdict(set)
    
    with open(annotation_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            hpo_id = parts[0]
            gene_symbol = parts[3]
            
            if hpo_id.startswith('HP:') and gene_symbol:
                gene_to_hpo[gene_symbol].add(hpo_id)
                hpo_to_genes[hpo_id].add(gene_symbol)
    
    # Convert to regular dict for pickling
    gene_to_hpo = {k: list(v) for k, v in gene_to_hpo.items()}
    hpo_to_genes = {k: list(v) for k, v in hpo_to_genes.items()}
    
    print(f"Loaded {len(gene_to_hpo)} genes")
    return gene_to_hpo, hpo_to_genes


# ==================== IC Calculation ====================

def calculate_ic(ont, hpo_to_genes):
    """Calculate Information Content for all HPO terms"""
    print("Calculating Information Content...")
    
    all_genes = set()
    for genes in hpo_to_genes.values():
        all_genes.update(genes)
    
    total_genes = len(all_genes)
    ic = {}
    
    for term in ont.terms():
        if not term.id.startswith('HP:'):
            continue
        
        genes_for_term = set()
        if term.id in hpo_to_genes:
            genes_for_term.update(hpo_to_genes[term.id])
        
        for child in term.subclasses():
            if child.id in hpo_to_genes:
                genes_for_term.update(hpo_to_genes[child.id])
        
        if len(genes_for_term) > 0:
            p = len(genes_for_term) / total_genes
            ic[term.id] = -np.log(p)
        else:
            ic[term.id] = 0.0
    
    print(f"Calculated IC for {len(ic)} terms")
    return ic


# ==================== Ancestor Precomputation ====================

def precompute_ancestors(ont, gene_to_hpo, cache_file='ancestor_cache.pkl'):
    """Precompute ancestors for all HPO terms"""
    if os.path.exists(cache_file):
        print(f"Loading ancestor cache from {cache_file}...")
        with open(cache_file, 'rb') as f:
            ancestor_cache = pickle.load(f)
        print(f"Loaded ancestor cache")
        return ancestor_cache
    
    print("Precomputing ancestors for all HPO terms...")
    
    all_hpo_terms = set()
    for hpo_set in gene_to_hpo.values():
        all_hpo_terms.update(hpo_set)
    
    ancestor_cache = {}
    
    for term_id in tqdm(all_hpo_terms, desc="Computing ancestors"):
        try:
            term = ont[term_id]
            ancestors = frozenset([term_id] + [t.id for t in term.superclasses()])
            ancestor_cache[term_id] = ancestors
        except:
            ancestor_cache[term_id] = frozenset([term_id])
    
    with open(cache_file, 'wb') as f:
        pickle.dump(ancestor_cache, f)
    
    print(f"Cached ancestors for {len(ancestor_cache)} terms")
    return ancestor_cache


# ==================== Resnik Similarity ====================

def resnik_similarity(term1_id, term2_id, ancestor_cache, ic):
    """
    Calculate Resnik similarity between two HPO terms
    (UNCHANGED - original implementation)
    """
    ancestors1 = ancestor_cache.get(term1_id, frozenset([term1_id]))
    ancestors2 = ancestor_cache.get(term2_id, frozenset([term2_id]))
    
    common = ancestors1 & ancestors2
    
    if not common:
        return 0.0
    
    max_ic = max(ic.get(ancestor, 0.0) for ancestor in common)
    return max_ic


def bma_similarity(patient_hpo_ids, gene_hpo_ids, ancestor_cache, ic):
    """
    Calculate Best Match Average similarity
    (UNCHANGED - original implementation)
    """
    if not patient_hpo_ids or not gene_hpo_ids:
        return 0.0
    
    total_sim = sum(
        max(
            resnik_similarity(p_term, g_term, ancestor_cache, ic)
            for g_term in gene_hpo_ids
        )
        for p_term in patient_hpo_ids
    )
    
    return total_sim / len(patient_hpo_ids)


# ==================== Gene Prediction with HPO Expansion ====================

def predict_genes_with_expansion(patient_hpo_ids, all_genes, gene_to_hpo, ont, ancestor_cache, ic):
    """
    Predict genes with HPO expansion
    
    CHANGE: Expands patient HPO terms before similarity calculation
    """
    if not patient_hpo_ids:
        return [], np.array([])
    
    # *** KEY CHANGE: Expand patient HPO terms ***
    patient_hpo_expanded = expand_hpo_with_ancestors(patient_hpo_ids, ont, max_levels=3)
    
    # Calculate scores using original Resnik + BMA (but with expanded terms)
    scores = np.array([
        bma_similarity(patient_hpo_expanded, gene_to_hpo[gene], ancestor_cache, ic)
        for gene in all_genes
    ])
    
    # Get top 20
    if len(scores) > 20:
        top_indices = np.argpartition(scores, -20)[-20:]
        top_indices = top_indices[np.argsort(scores[top_indices])[::-1]]
    else:
        top_indices = np.argsort(scores)[::-1]
    
    ranked_genes = [(all_genes[i], scores[i]) for i in top_indices[:20]]
    
    return ranked_genes, scores


# ==================== Parallel Worker ====================

def predict_single_case_parallel(args):
    """Worker function for parallel processing"""
    case_data, all_genes, gene_to_hpo, ont, ancestor_cache, ic = args
    
    case_id = case_data['case_id']
    true_gene = case_data['true_gene']
    patient_hpo_ids = case_data['patient_hpo_ids']
    
    if not patient_hpo_ids:
        return None
    
    # Predict with HPO expansion
    top_20_ranked, all_scores = predict_genes_with_expansion(
        patient_hpo_ids, all_genes, gene_to_hpo, ont, ancestor_cache, ic
    )
    
    # Calculate metrics
    true_gene_idx = all_genes.index(true_gene)
    true_gene_score = all_scores[true_gene_idx]
    true_gene_rank = (all_scores > true_gene_score).sum() + 1
    
    top_genes = [g for g, s in top_20_ranked]
    
    result = {
        'case_id': case_id,
        'true_gene': true_gene,
        'true_gene_rank': int(true_gene_rank),
        'true_gene_score': float(true_gene_score),
        'in_top_1': true_gene in top_genes[:1],
        'in_top_5': true_gene in top_genes[:5],
        'in_top_10': true_gene in top_genes[:10],
        'in_top_20': true_gene in top_genes[:20],
        'top_5_predictions': ','.join([g for g, s in top_20_ranked[:5]]),
        'num_patient_hpo': len(patient_hpo_ids),
        'num_patient_hpo_expanded': len(expand_hpo_with_ancestors(patient_hpo_ids, ont, max_levels=3)),
        'num_gene_hpo': len(gene_to_hpo[true_gene])
    }
    
    # For AUC
    y_true = np.zeros(len(all_genes))
    y_true[true_gene_idx] = 1
    
    return result, y_true, all_scores


# ==================== Helper Functions ====================

def extract_gene_from_variant(variant_str):
    """Extract gene symbol from variant description"""
    if pd.isna(variant_str) or variant_str == 'Not reported':
        return None
    
    variant_str = str(variant_str).strip()
    
    match = re.match(r'^([A-Z][A-Z0-9]+):', variant_str)
    if match:
        gene = match.group(1)
        if not any(gene.startswith(p) for p in ['NM_', 'NR_', 'NC_', 'CHR']):
            return gene
    
    match = re.match(r'^([A-Z][A-Z0-9]+)\s*[\(\[]', variant_str)
    if match:
        return match.group(1)
    
    match = re.match(r'^([A-Z][A-Z0-9]+)\s+[cpgn]\.', variant_str)
    if match:
        return match.group(1)
    
    return None


def parse_hpo_ids(hpo_str):
    """Parse comma-separated HPO IDs"""
    if pd.isna(hpo_str) or hpo_str == 'Not reported':
        return []
    
    return [x.strip() for x in str(hpo_str).split(',') if x.strip().startswith('HP:')]


# ==================== Plotting ====================

def plot_roc_pr_curves(y_true, y_scores, roc_auc, pr_auc):
    """Plot ROC and Precision-Recall curves"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # ROC Curve
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    ax1.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.4f})')
    ax1.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random')
    ax1.set_xlim([0.0, 1.0])
    ax1.set_ylim([0.0, 1.05])
    ax1.set_xlabel('False Positive Rate')
    ax1.set_ylabel('True Positive Rate')
    ax1.set_title('ROC Curve')
    ax1.legend(loc="lower right")
    ax1.grid(alpha=0.3)
    
    # Precision-Recall Curve
    precision, recall, _ = precision_recall_curve(y_true, y_scores)
    ax2.plot(recall, precision, color='darkorange', lw=2, label=f'PR curve (AUC = {pr_auc:.4f})')
    ax2.set_xlim([0.0, 1.0])
    ax2.set_ylim([0.0, 1.05])
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Precision')
    ax2.set_title('Precision-Recall Curve')
    ax2.legend(loc="lower left")
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('roc_pr_curves.png', dpi=300, bbox_inches='tight')
    print("Saved ROC and PR curves to roc_pr_curves.png")
    plt.close()


def plot_rank_distribution(results_df):
    """Plot distribution of true gene ranks"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    ranks = results_df['true_gene_rank']
    
    ax1.hist(ranks, bins=50, edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Rank of True Gene')
    ax1.set_ylabel('Number of Cases')
    ax1.set_title('Distribution of True Gene Ranks')
    ax1.axvline(ranks.median(), color='red', linestyle='--', linewidth=2, label=f'Median = {ranks.median():.0f}')
    ax1.axvline(ranks.mean(), color='orange', linestyle='--', linewidth=2, label=f'Mean = {ranks.mean():.1f}')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    sorted_ranks = np.sort(ranks)
    cumulative = np.arange(1, len(sorted_ranks) + 1) / len(sorted_ranks) * 100
    ax2.plot(sorted_ranks, cumulative, linewidth=2)
    ax2.set_xlabel('Rank of True Gene')
    ax2.set_ylabel('Cumulative Percentage (%)')
    ax2.set_title('Cumulative Distribution of Ranks')
    ax2.grid(alpha=0.3)
    
    for k in [1, 5, 10, 20]:
        pct = (ranks <= k).sum() / len(ranks) * 100
        ax2.axvline(k, color='red', linestyle='--', alpha=0.5)
        ax2.text(k, pct + 2, f'Top-{k}\n{pct:.1f}%', ha='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('rank_distribution.png', dpi=300, bbox_inches='tight')
    print("Saved rank distribution to rank_distribution.png")
    plt.close()


# ==================== Main Pipeline ====================

def run_gene_prediction_with_expansion(data_file='data/PAVS_final_data.tsv',
                                       hpo_file='resources/hp.obo',
                                       annotation_file='resources/phenotype_to_genes.txt',
                                       output_file='gene_predictions_expanded.tsv',
                                       n_cores=None):
    
    if n_cores is None:
        n_cores = max(1, cpu_count() - 1)
    
    print(f"Using {n_cores} CPU cores")
    
    ont = load_hpo_ontology(hpo_file)
    gene_to_hpo, hpo_to_genes = load_hpo_gene_annotations(annotation_file)
    all_genes = sorted(gene_to_hpo.keys())
    
    ic = calculate_ic(ont, hpo_to_genes)
    ancestor_cache = precompute_ancestors(ont, gene_to_hpo)
    
    data = pd.read_csv(data_file, sep='\t')
    print(f"Loaded {len(data)} cases")
    
    data_with_both = data[
        (data['phenotypicFeatureIds'] != 'Not reported') & 
        (data['genomicVariants'] != 'Not reported')
    ].copy()
    
    data_with_both['true_gene'] = data_with_both['genomicVariants'].apply(extract_gene_from_variant)
    data_with_gene = data_with_both[data_with_both['true_gene'].notna()].copy()
    data_final = data_with_gene[data_with_gene['true_gene'].isin(all_genes)].copy()
    
    print(f"Cases to evaluate: {len(data_final)}")
    
    print("\n" + "="*60)
    print("Running Predictions with HPO Expansion")
    print("="*60)
    
    case_args = []
    for idx, row in data_final.iterrows():
        case_data = {
            'case_id': row['ID'],
            'true_gene': row['true_gene'],
            'patient_hpo_ids': parse_hpo_ids(row['phenotypicFeatureIds'])
        }
        case_args.append((case_data, all_genes, gene_to_hpo, ont, ancestor_cache, ic))
    
    results = []
    all_y_true = []
    all_y_scores = []
    
    with Pool(n_cores) as pool:
        for result in tqdm(
            pool.imap(predict_single_case_parallel, case_args),
            total=len(case_args),
            desc="Predicting"
        ):
            if result is not None:
                case_result, y_true, y_scores = result
                results.append(case_result)
                all_y_true.extend(y_true)
                all_y_scores.extend(y_scores)
    
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nResults saved to {output_file}")
    
    # Calculate metrics
    print("\n" + "="*60)
    print("Performance Metrics")
    print("="*60)
    
    roc_auc = roc_auc_score(all_y_true, all_y_scores)
    pr_auc = average_precision_score(all_y_true, all_y_scores)
    mrr = (1.0 / results_df['true_gene_rank']).mean()
    
    print(f"\nROC-AUC: {roc_auc:.4f}")
    print(f"PR-AUC:  {pr_auc:.4f}")
    print(f"MRR:     {mrr:.4f}")
    
    print("\nTop-k Accuracy:")
    total = len(results_df)
    for k in [1, 5, 10, 20, 50]:
        count = results_df[f'in_top_{k}'].sum() if k <= 20 else (results_df['true_gene_rank'] <= k).sum()
        print(f"  Top {k:2d}: {count:4d} / {total} ({count/total*100:5.1f}%)")
    
    print(f"\nMedian rank: {results_df['true_gene_rank'].median():.0f}")
    print(f"Mean rank:   {results_df['true_gene_rank'].mean():.1f}")
    
    # HPO expansion statistics
    print("\n" + "="*60)
    print("HPO Expansion Statistics")
    print("="*60)
    avg_original = results_df['num_patient_hpo'].mean()
    avg_expanded = results_df['num_patient_hpo_expanded'].mean()
    print(f"Average patient HPO terms (original): {avg_original:.1f}")
    print(f"Average patient HPO terms (expanded): {avg_expanded:.1f}")
    print(f"Expansion factor: {avg_expanded/avg_original:.2f}x")
    
    plot_roc_pr_curves(all_y_true, all_y_scores, roc_auc, pr_auc)
    plot_rank_distribution(results_df)
    
    summary = {
        'roc_auc': float(roc_auc),
        'pr_auc': float(pr_auc),
        'mrr': float(mrr),
        'top_1': float(results_df['in_top_1'].sum() / total),
        'top_5': float(results_df['in_top_5'].sum() / total),
        'top_10': float(results_df['in_top_10'].sum() / total),
        'median_rank': float(results_df['true_gene_rank'].median()),
        'mean_rank': float(results_df['true_gene_rank'].mean()),
        'avg_hpo_original': float(avg_original),
        'avg_hpo_expanded': float(avg_expanded)
    }
    
    with open('prediction_metrics.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\nSummary metrics saved to prediction_metrics.json")
    
    return results_df


# ==================== Main ====================

def main():
    """Main entry point"""
    results = run_gene_prediction_with_expansion(
        data_file='data/PAVS_final_data.tsv',
        hpo_file='resources/hp.obo',
        annotation_file='resources/phenotype_to_genes.txt',
        output_file='gene_predictions_expanded_6_level.tsv',
        n_cores=None
    )


if __name__ == "__main__":
    main()