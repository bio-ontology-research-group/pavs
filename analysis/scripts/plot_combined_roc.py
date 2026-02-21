
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def plot_combined_roc(results_path, output_path):
    print(f"Loading results from {results_path}...")
    df = pd.read_csv(results_path)
    
    plt.figure(figsize=(10, 8))
    
    for label in sorted(df['group'].unique()):
        g_df = df[df['group'] == label]
        ranks = g_df['rank'].values
        num_cand = g_df['num_candidates'].values
        
        # Calculate AUC for label
        aucs = (num_cand - ranks) / (num_cand - 1)
        mean_auc = np.mean(aucs)
        
        # ROC Curve
        sorted_ranks = np.sort(ranks)
        unique_ranks, counts = np.unique(sorted_ranks, return_counts=True)
        tpr = np.cumsum(counts) / len(ranks)
        fpr = (unique_ranks - 1) / (num_cand[0] - 1)
        
        # Add point (0,0) and (1,1)
        fpr = np.insert(fpr, 0, 0)
        tpr = np.insert(tpr, 0, 0)
        if fpr[-1] < 1.0:
            fpr = np.append(fpr, 1.0)
            tpr = np.append(tpr, 1.0)
            
        plt.plot(fpr, tpr, label=f'{label} (AUC = {mean_auc:.4f})', linewidth=2)
    
    plt.plot([0, 1], [0, 1], 'k--', alpha=0.5)
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title('Combined Phenotype Similarity ROC Curves', fontsize=14)
    plt.legend(loc='lower right', fontsize=10)
    plt.grid(True, linestyle=':', alpha=0.6)
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Combined ROC curve saved to {output_path}")

if __name__ == "__main__":
    plot_combined_roc('analysis/phenotype_similarity_full.csv', 'analysis/combined_roc_curves.png')
