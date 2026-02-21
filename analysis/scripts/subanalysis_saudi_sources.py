import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def subanalysis_saudi_sources(results_path, annotated_path, output_dir):
    print(f"Loading results from {results_path}...")
    df = pd.read_csv(results_path)
    
    # Map PAVS ID to Source Name
    print(f"Loading source mapping from {annotated_path}...")
    df_map = pd.read_csv(annotated_path, sep='\t', usecols=['pavs_id', 'source_file'])
    id_to_source = dict(zip(df_map['pavs_id'], df_map['source_file']))
    
    # Extract source name
    df['source_file'] = df['pavs_id'].map(id_to_source)
    
    # Filter only Saudi sources and exclude Literature and DDD
    saudi_df = df[df['group'] == 'Saudi'].copy()
    saudi_df = saudi_df.dropna(subset=['source_file'])
    
    os.makedirs(output_dir, exist_ok=True)
    
    stats = []
    
    # Setup for plotting
    plt.figure(1, figsize=(10, 8))
    plt.title('ROC Curves by Saudi Source', fontsize=14)
    
    plt.figure(2, figsize=(10, 8))
    plt.title('Precision-Recall Curves by Saudi Source', fontsize=14)
    
    for source in sorted(saudi_df['source_file'].unique()):
        g_df = saudi_df[saudi_df['source_file'] == source]
        if g_df.empty: continue
        
        ranks = g_df['rank'].values
        num_cand = g_df['num_candidates'].values
        
        # Metrics
        hits_at_1 = np.mean(ranks <= 1)
        hits_at_10 = np.mean(ranks <= 10)
        hits_at_50 = np.mean(ranks <= 50)
        mean_rank = np.mean(ranks)
        median_rank = np.median(ranks)
        aucs = (num_cand - ranks) / (num_cand - 1)
        mean_auc = np.mean(aucs)
        mrr = np.mean(1.0 / ranks)
        
        stats.append({
            'Source': source,
            'n': len(g_df),
            'Mean Rank': mean_rank,
            'Median Rank': median_rank,
            'Hits@1': hits_at_1,
            'Hits@10': hits_at_10,
            'Hits@50': hits_at_50,
            'Mean AUC': mean_auc,
            'MRR': mrr
        })
        
        # ROC Plotting
        sorted_ranks = np.sort(ranks)
        unique_ranks, counts = np.unique(sorted_ranks, return_counts=True)
        tpr = np.cumsum(counts) / len(ranks)
        fpr = (unique_ranks - 1) / (num_cand[0] - 1)
        fpr = np.insert(fpr, 0, 0); tpr = np.insert(tpr, 0, 0)
        
        plt.figure(1)
        plt.plot(fpr, tpr, label=f'{source} (AUC={mean_auc:.3f})', linewidth=2)
        
        # PR Plotting
        max_k = 1000
        ks = np.arange(1, max_k + 1)
        recalls = []
        precisions = []
        for k in ks:
            hits = np.sum(ranks <= k)
            recalls.append(hits / len(ranks))
            precisions.append(hits / (k * len(ranks)))
        
        plt.figure(2)
        plt.plot(recalls, precisions, label=f'{source}', linewidth=2)

    # Save ROC
    plt.figure(1)
    plt.plot([0, 1], [0, 1], 'k--', alpha=0.5)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right', fontsize=9)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.savefig(os.path.join(output_dir, 'saudi_sources_roc.png'), dpi=300)
    plt.close()
    
    # Save PR
    plt.figure(2)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc='upper right', fontsize=9)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.ylim([0, 1.05])
    plt.xlim([0, 1.05])
    plt.savefig(os.path.join(output_dir, 'saudi_sources_pr.png'), dpi=300)
    plt.close()
    
    # Save Stats Table
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(os.path.join(output_dir, 'saudi_sources_stats.csv'), index=False)
    
    # Print Markdown table
    print("\n| Source | n | Mean Rank | Median Rank | Hits@1 | Hits@10 | Mean AUC | MRR |")
    print("| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |")
    for s in stats:
        print(f"| {s['Source']} | {s['n']} | {s[ 'Mean Rank']:.1f} | {s['Median Rank']:.1f} | {s['Hits@1']:.2%} | {s['Hits@10']:.2%} | {s['Mean AUC']:.4f} | {s['MRR']:.4f} |")

if __name__ == "__main__":
    subanalysis_saudi_sources('analysis/phenotype_similarity_full.csv', 'data/combined_annotated.tsv', 'analysis/')
