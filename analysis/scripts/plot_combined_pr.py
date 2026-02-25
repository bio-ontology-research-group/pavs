
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def plot_combined_pr(results_path, output_path):
    print(f"Loading results from {results_path}...")
    df = pd.read_csv(results_path)
    
    plt.figure(figsize=(10, 8))
    
    for label in sorted(df['group'].unique()):
        g_df = df[df['group'] == label]
        ranks = g_df['rank'].values
        num_cand = g_df['num_candidates'].values
        
        # In this task, there is 1 positive gene and (num_cand - 1) negative genes.
        # We can simulate the Precision-Recall curve based on the ranks.
        # For a single case with rank r:
        # Precision(k) = 1/k if r <= k, else 0
        # Recall(k) = 1 if r <= k, else 0
        
        # Aggregating over cases:
        # Recall(k) = Proportion of cases with rank <= k
        # Precision(k) = Mean(1/k for cases with rank <= k, 0 otherwise) 
        #              = (Count of cases with rank <= k) / (Sum of k over all cases) ? No.
        # Precision(k) = (Count of cases with rank <= k) / (k * Total cases)
        
        max_k = 1000 # Let's plot up to rank 1000
        ks = np.arange(1, max_k + 1)
        recalls = []
        precisions = []
        
        for k in ks:
            hits = np.sum(ranks <= k)
            recall = hits / len(ranks)
            precision = hits / (k * len(ranks))
            recalls.append(recall)
            precisions.append(precision)
            
        plt.plot(recalls, precisions, label=f'{label}', linewidth=2)
    
    plt.xlabel('Recall', fontsize=12)
    plt.ylabel('Precision', fontsize=12)
    plt.title('Precision-Recall Curves (Ranking up to Rank 1000)', fontsize=14)
    plt.legend(loc='upper right', fontsize=10)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.ylim([0, 1.05])
    plt.xlim([0, 1.05])
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"PR curves saved to {output_path}")

if __name__ == "__main__":
    plot_combined_pr('analysis/phenotype_similarity_comprehensive.csv', 'analysis/combined_pr_curves.png')
