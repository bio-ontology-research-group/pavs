import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def full_cohort_analysis(results_path, output_dir):
    print(f"Loading results from {results_path}...")
    df_all = pd.read_csv(results_path)
    
    os.makedirs(output_dir, exist_ok=True)
    
    combinations = [
        ('intrinsic', 'lin'),
        ('intrinsic', 'resnik'),
        ('extrinsic', 'lin'),
        ('extrinsic', 'resnik')
    ]
    
    cohorts = ['Saudi', 'DDD', 'Literature']
    
    full_stats = []
    
    for ic_type, sim_method in combinations:
        print(f"Analyzing {ic_type}-{sim_method}...")
        df = df_all[(df_all['ic_type'] == ic_type) & (df_all['sim_method'] == sim_method)].copy()
        
        # ROC plot per combination — one curve per cohort
        plt.figure(figsize=(10, 7))
        plt.title(f'ROC Curves by Cohort ({ic_type.capitalize()} {sim_method.capitalize()})', fontsize=14)
        
        for cohort in cohorts:
            g_df = df[df['group'] == cohort].copy()
            if g_df.empty:
                print(f"  No data for cohort: {cohort}")
                continue
            
            ranks = g_df['rank'].values
            num_cand = g_df['num_candidates'].values
            
            # Metrics
            hits_at_1  = np.mean(ranks <= 1)
            hits_at_10 = np.mean(ranks <= 10)
            mean_auc   = np.mean((num_cand - ranks) / (num_cand - 1))
            mrr        = np.mean(1.0 / ranks)
            
            full_stats.append({
                'IC Type':    ic_type,
                'Similarity': sim_method,
                'Cohort':     cohort,
                'n':          len(g_df),
                'Hits@1':     hits_at_1,
                'Hits@10':    hits_at_10,
                'Mean AUC':   mean_auc,
                'MRR':        mrr
            })
            
            # ROC curve
            sorted_ranks = np.sort(ranks)
            unique_ranks, counts = np.unique(sorted_ranks, return_counts=True)
            tpr = np.cumsum(counts) / len(ranks)
            fpr = (unique_ranks - 1) / (num_cand[0] - 1)
            fpr = np.insert(fpr, 0, 0)
            tpr = np.insert(tpr, 0, 0)
            plt.plot(fpr, tpr, label=f'{cohort} (AUC={mean_auc:.4f}, n={len(g_df)})', linewidth=2)
        
        plt.plot([0, 1], [0, 1], 'k--', alpha=0.5)
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.legend(loc='lower right', fontsize=10)
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'cohorts_roc_{ic_type}_{sim_method}.png'), dpi=300)
        plt.close()
    
    # Save full stats
    stats_df = pd.DataFrame(full_stats)
    stats_df.to_csv(os.path.join(output_dir, 'cohorts_comprehensive_stats.csv'), index=False)
    
    # Print summary table — all 4 methods x 3 cohorts
    print("\n### Full Results: All Combinations x All Cohorts")
    print("| Cohort | Method | n | AUC | MRR | Hits@1 (%) | Hits@10 (%) |")
    print("| :--- | :--- | :--- | :--- | :--- | :--- | :--- |")
    for _, s in stats_df.sort_values(['Cohort', 'IC Type', 'Similarity']).iterrows():
        method = f"{s['IC Type'].capitalize()} {s['Similarity'].capitalize()}"
        print(f"| {s['Cohort']} | {method} | "
              f"{s['n']} | {s['Mean AUC']:.4f} | {s['MRR']:.4f} | "
              f"{s['Hits@1']*100:.2f} | {s['Hits@10']*100:.2f} |")

# Save markdown table to file
    lines = []
    
    lines.append("### Full Results: All Combinations x All Cohorts\n")
    lines.append("| Cohort | Method | n | AUC | MRR | Hits@1 (%) | Hits@10 (%) |\n")
    lines.append("| :--- | :--- | :--- | :--- | :--- | :--- | :--- |\n")
    for _, s in stats_df.sort_values(['Cohort', 'IC Type', 'Similarity']).iterrows():
        method = f"{s['IC Type'].capitalize()} {s['Similarity'].capitalize()}"
        lines.append(f"| {s['Cohort']} | {method} | "
                     f"{s['n']} | {s['Mean AUC']:.4f} | {s['MRR']:.4f} | "
                     f"{s['Hits@1']*100:.2f} | {s['Hits@10']*100:.2f} |\n")
    
    lines.append("\n### Best Method per Cohort (by MRR)\n")
    lines.append("| Cohort | Method | n | AUC | MRR | Hits@1 (%) | Hits@10 (%) |\n")
    lines.append("| :--- | :--- | :--- | :--- | :--- | :--- | :--- |\n")
    best = stats_df.loc[stats_df.groupby('Cohort')['MRR'].idxmax()]
    for _, s in best.iterrows():
        method = f"{s['IC Type'].capitalize()} {s['Similarity'].capitalize()}"
        lines.append(f"| {s['Cohort']} | {method} | "
                     f"{s['n']} | {s['Mean AUC']:.4f} | {s['MRR']:.4f} | "
                     f"{s['Hits@1']*100:.2f} | {s['Hits@10']*100:.2f} |\n")
    
    md_path = os.path.join(output_dir, 'cohorts_summary.md')
    with open(md_path, 'w') as f:
        f.writelines(lines)
    print(f"\nMarkdown summary saved to {md_path}")

if __name__ == "__main__":
    full_cohort_analysis(
        results_path='analysis/phenotype_similarity_comprehensive.csv',
        output_dir='analysis/'
    )

