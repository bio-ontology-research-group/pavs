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

| Cohort | Method | Mean AUC | Mean AUPR (MRR) | Hits@1 | Hits@10 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **DDD** | Extrinsic Resnik | **0.9941** | **0.7134** | **60.50%** | **85.17%** |
| **DDD** | Extrinsic Lin | 0.9941 | 0.6841 | 55.86% | 83.99% |
| **Literature** | Extrinsic Resnik | **0.9813** | **0.6712** | **56.36%** | **79.60%** |
| **Literature** | Extrinsic Lin | 0.9811 | 0.6465 | 51.60% | 77.91% |
| **Saudi** | Extrinsic Resnik | 0.8861 | **0.0625** | **3.35%** | **11.29%** |
| **Saudi** | Intrinsic Resnik | **0.8880** | 0.0574 | 2.81% | 10.97% |

### Saudi Source Sub-analysis (Extrinsic Resnik)

The Saudi cohort consists of six distinct sources. Performance using the top-performing **Extrinsic Resnik** method:

| Source | n | Hits@1 | Hits@10 | Mean AUC | MRR |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **marwa-variants** | 1,064 | 7.89% | 20.68% | 0.9307 | 0.1238 |
| **ahmed-variants** | 200 | 2.50% | 9.00% | 0.8628 | 0.0482 |
| **PMC7082194** | 122 | 1.64% | 7.38% | 0.8956 | 0.0455 |
| **PMC6562004** | 1,228 | 1.63% | 8.71% | 0.8728 | 0.0423 |
| **fawzan-variants** | 460 | 1.09% | 6.74% | 0.8654 | 0.0346 |
| **ahmed-pmid28454995** | 88 | 1.14% | 3.41% | 0.8892 | 0.0251 |

## Key Files
- `phenotype_similarity_comprehensive.csv`: Raw results for all 4 methods.
- `cohorts_comprehensive_stats.csv`: Full metrics (AUC, MRR, Hits@1, Hits@10) for all 4 IC/similarity combinations across all three cohorts
- `saudi_sources_comprehensive_stats.csv`: Statistical breakdown for all Saudi sub-cohorts and methods.
- `saudi_sources_roc_extrinsic_resnik.png`: ROC curves for Saudi sources (top method).
- `combined_roc_curves.png`: Comparative ROC curves across major cohorts.
- `METHOD.md`: Detailed documentation of the computational methodology.

## Summary Observations
1. **Method Optimization**: Using Extrinsic (gene-based) IC with Resnik similarity significantly improves prioritization accuracy, especially Hits@1 and MRR, by focusing on clinically informative terms.
2. **DDD/Literature vs. Saudi**: Curated datasets show exceptional alignment. The Saudi cohort presents greater real-world challenge, with lower Hits@1 indicating phenotypic heterogeneity or sparse clinical descriptions.
3. **Internal Variability**: Within the Saudi cohort, the Alkuraya dataset (marwa-variants) shows the highest quality phenotypic descriptions, achieving a 20.7% Hits@10.
