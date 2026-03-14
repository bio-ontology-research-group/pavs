# Phenotype Semantic Similarity Analysis

## Overview
This analysis compares Saudi clinical cases against the Deciphering Developmental Disorders (DDD) cohort and literature-curated Phenopackets (public store 0.1.26). The goal is to evaluate how well clinical phenotypes prioritize the causative gene using semantic similarity.

## Methodology
- **Similarity Metrics**: Lin and Resnik semantic similarity.
- **Aggregation**: Best Match Average (BMA).
- **Ontology**: Human Phenotype Ontology (HPO).
- **IC Calculation**: Both Intrinsic (topology-based) and Extrinsic (gene-based) Information Content.
- **Ranking**: Each case is ranked against all genes in the reference dataset (~5,200 genes).

## Performance Breakdown

We compare four combinations of Information Content (IC) and Similarity Measures. **Extrinsic Resnik** consistently provides the best prioritization performance.

| Cohort | Method | n | Mean AUC | MRR | Hits@1 | Hits@10 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Literature** | Extrinsic Resnik | 8,887 | **0.9828** | **0.6542** | **57.78%** | **77.87%** |
| **Literature** | Extrinsic Lin | 8,887 | 0.9811 | 0.6465 | 56.82% | 77.91% |
| **DDD** | Extrinsic Resnik | 1,443 | **0.9951** | **0.7063** | **62.02%** | **85.45%** |
| **DDD** | Extrinsic Lin | 1,443 | 0.9941 | 0.6841 | 60.29% | 83.99% |
| **Saudi** | Extrinsic Resnik | 3,088 | **0.8915** | **0.0676** | **3.69%** | **12.47%** |
| **Saudi** | Intrinsic Resnik | 3,088 | 0.8896 | 0.0585 | 2.91% | 11.24% |

### Saudi Source Sub-analysis (Extrinsic Resnik)

The Saudi cohort consists of six distinct sources. Performance using the top-performing **Extrinsic Resnik** method:

| Source | n | Hits@1 | Hits@10 | Mean AUC | MRR |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **marwa-variants** | 1,018 | 8.35% | 21.91% | 0.9315 | 0.1287 |
| **ahmed-variants** | 198 | 2.02% | 8.59% | 0.8599 | 0.0423 |
| **PMC6562004** | 1,213 | 1.57% | 8.74% | 0.8715 | 0.0402 |
| **PMC7082194** | 122 | 1.64% | 5.74% | 0.8933 | 0.0344 |
| **fawzan-variants** | 453 | 0.66% | 5.96% | 0.8685 | 0.0306 |
| **ahmed-pmid28454995** | 84 | 1.19% | 5.95% | 0.8909 | 0.0296 |

## Key Files
- `phenotype_similarity_comprehensive.csv`: Raw results for all 4 methods.
- `cohorts_comprehensive_stats.csv`: Full metrics (AUC, MRR, Hits@1, Hits@10) for all 4 IC/similarity combinations across all three cohorts
- `saudi_sources_comprehensive_stats.csv`: Statistical breakdown for all Saudi sub-cohorts and methods.
- `saudi_sources_roc_extrinsic_resnik.png`: ROC curves for Saudi sources (top method).
- `combined_roc_curves.png`: Comparative ROC curves across major cohorts.
- `METHOD.md`: Detailed documentation of the computational methodology.

## Summary Observations
1. **Method Optimization**: Using Extrinsic (gene-based) IC with Resnik similarity significantly improves prioritization accuracy, especially Hits@1 and MRR, by focusing on clinically informative terms.
2. **Data Density Impact**: The Literature (Avg 21.5 HPO terms) and DDD (Avg 15.3 terms) cohorts significantly outperform the Saudi cohort (Avg 3.3 terms).
3. **Internal Variability**: Within the Saudi cohort, the Alkuraya dataset (marwa-variants) shows the highest quality phenotypic descriptions, achieving a 21.9% Hits@10.
