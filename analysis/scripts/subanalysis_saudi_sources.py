
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def subanalysis_comprehensive(results_path, annotated_path, output_dir):
    print(f"Loading comprehensive results from {results_path}...")
    df_all = pd.read_csv(results_path)
    
    print(f"Loading source mapping from {annotated_path}...")
    df_map = pd.read_csv(annotated_path, sep='\t', usecols=['pavs_id', 'source_file'])
    id_to_source = dict(zip(df_map['pavs_id'], df_map['source_file']))
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Define the 4 combinations
    combinations = [
        ('intrinsic', 'lin'),
        ('intrinsic', 'resnik'),
        ('extrinsic', 'lin'),
        ('extrinsic', 'resnik')
    ]
    
    full_stats = []
    
    for ic_type, sim_method in combinations:
        print(f"Analyzing {ic_type}-{sim_method}...")
        df = df_all[(df_all['ic_type'] == ic_type) & (df_all['sim_method'] == sim_method)].copy()
        df['source_file'] = df['pavs_id'].map(id_to_source)
        
        # Filter only Saudi sources for this sub-analysis
        saudi_df = df[df['group'] == 'Saudi'].copy()
        saudi_df = saudi_df.dropna(subset=['source_file'])
        
        # Plotting Setup
        plt.figure(figsize=(12, 8))
        plt.title(f'ROC Curves by Saudi Source ({ic_type.capitalize()} {sim_method.capitalize()})', fontsize=14)
        
        for source in sorted(saudi_df['source_file'].unique()):
            g_df = saudi_df[saudi_df['source_file'] == source]
            if g_df.empty: continue
            
            ranks = g_df['rank'].values
            num_cand = g_df['num_candidates'].values
            
            # Metrics
            hits_at_1 = np.mean(ranks <= 1)
            hits_at_10 = np.mean(ranks <= 10)
            mean_auc = np.mean((num_cand - ranks) / (num_cand - 1))
            mrr = np.mean(1.0 / ranks)
            
            full_stats.append({
                'IC Type': ic_type,
                'Similarity': sim_method,
                'Source': source,
                'n': len(g_df),
                'Hits@1': hits_at_1,
                'Hits@10': hits_at_10,
                'Mean AUC': mean_auc,
                'MRR': mrr
            })
            
            # ROC Plotting
            sorted_ranks = np.sort(ranks)
            unique_ranks, counts = np.unique(sorted_ranks, return_counts=True)
            tpr = np.cumsum(counts) / len(ranks)
            fpr = (unique_ranks - 1) / (num_cand[0] - 1)
            fpr = np.insert(fpr, 0, 0); tpr = np.insert(tpr, 0, 0)
            plt.plot(fpr, tpr, label=f'{source} (AUC={mean_auc:.3f})', linewidth=1.5)

        plt.plot([0, 1], [0, 1], 'k--', alpha=0.5)
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.legend(loc='lower right', fontsize=8)
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.savefig(os.path.join(output_dir, f'saudi_sources_roc_{ic_type}_{sim_method}.png'), dpi=300)
        plt.close()

    # Save comprehensive stats
    stats_df = pd.DataFrame(full_stats)
    stats_df.to_csv(os.path.join(output_dir, 'saudi_sources_comprehensive_stats.csv'), index=False)
    
    # Generate a Markdown summary table for the top combination (Extrinsic Resnik usually best)
    print("\n### Comparative Summary for Saudi Sources (Extrinsic Resnik)")
    best_df = stats_df[(stats_df['IC Type'] == 'extrinsic') & (stats_df['Similarity'] == 'resnik')]
    print("| Source | n | Hits@1 | Hits@10 | Mean AUC | MRR |")
    print("| :--- | :--- | :--- | :--- | :--- | :--- |")
    for _, s in best_df.sort_values('MRR', ascending=False).iterrows():
        print(f"| {s['Source']} | {s['n']} | {s['Hits@1']:.2%} | {s['Hits@10']:.2%} | {s['Mean AUC']:.4f} | {s['MRR']:.4f} |")

    # Overall comparison of methods for the Saudi cohort
    print("\n### Method Comparison (Averaged over all Saudi cases)")
    method_cmp = stats_df.groupby(['IC Type', 'Similarity']).agg({
        'Hits@1': 'mean',
        'Hits@10': 'mean',
        'Mean AUC': 'mean',
        'MRR': 'mean'
    }).reset_index()
    print("| IC Type | Similarity | Hits@1 | Hits@10 | Mean AUC | MRR |")
    print("| :--- | :--- | :--- | :--- | :--- | :--- |")
    for _, s in method_cmp.iterrows():
        print(f"| {s['IC Type']} | {s['Similarity']} | {s['Hits@1']:.2%} | {s['Hits@10']:.2%} | {s['Mean AUC']:.4f} | {s['MRR']:.4f} |")

if __name__ == "__main__":
    subanalysis_comprehensive('analysis/phenotype_similarity_comprehensive.csv', 'data/combined_annotated.tsv', 'analysis/')
