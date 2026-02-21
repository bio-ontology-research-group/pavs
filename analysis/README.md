# Phenotype Semantic Similarity Analysis

## Overview
This analysis compares Saudi clinical cases against the Deciphering Developmental Disorders (DDD) cohort and literature-curated Phenopackets (public store 0.1.26). The goal is to evaluate how well clinical phenotypes prioritize the causative gene using semantic similarity.

## Methodology
- **Similarity Metric**: Lin's semantic similarity with Best Match Average (BMA) aggregation.
- **Ontology**: Human Phenotype Ontology (HPO).
- **Reference**: HPO associations from `genes_to_phenotype.txt`.
- **Ranking**: Each case is ranked against all genes in the reference dataset (~5,200 genes).

## Cohorts
- **Saudi** (n=3,162): Cases from Saudi Arabian clinical studies (Alkuraya, Marwa, Fawzan, etc.).
- **DDD** (n=1,443): Cases from the DDD project as a baseline comparison.
- **Literature** (n=8,887): High-quality manually curated Phenopackets from the public store.

## Performance Breakdown

| Cohort | n | Mean Rank | Median Rank | Hits@1 | Hits@10 | Hits@50 | Mean AUC | Mean AUPR (MRR) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **DDD** | 1,443 | 39.76 | 1.00 | 55.86% | 80.32% | 90.23% | 0.9925 | 0.6436 |
| **Literature** | 8,887 | 106.21 | 1.00 | 51.60% | 72.65% | 83.31% | 0.9797 | 0.5920 |
| **Saudi** | 3,162 | 567.64 | 198.00 | 2.81% | 10.97% | 26.47% | 0.8909 | 0.0574 |

### Saudi Source Sub-analysis

The Saudi cohort consists of six distinct sources with varying phenotypic data quality and specificity:

| Source | n | Mean Rank | Median Rank | Hits@1 | Hits@10 | Mean AUC | Mean AUPR (MRR) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **marwa-variants** | 1,064 | 364.6 | 93.0 | 5.36% | 18.42% | 0.9300 | 0.0982 |
| **PMC7082194** | 122 | 492.9 | 175.0 | 1.64% | 8.20% | 0.9053 | 0.0438 |
| **ahmed-pmid28454995** | 88 | 582.5 | 239.5 | 1.14% | 4.55% | 0.8880 | 0.0278 |
| **PMC6562004** | 1,228 | 680.6 | 329.5 | 1.55% | 7.17% | 0.8691 | 0.0367 |
| **fawzan-variants** | 460 | 696.1 | 297.0 | 0.87% | 7.39% | 0.8661 | 0.0308 |
| **ahmed-variants** | 200 | 698.1 | 298.5 | 3.00% | 7.50% | 0.8657 | 0.0501 |

## Key Files
- `phenotype_similarity_full.csv`: Raw results for every case.
- `interpretation_full.txt`: Detailed performance metrics for each major cohort.
- `saudi_sources_stats.csv`: Performance metrics for each Saudi source.
- `combined_roc_curves.png`: Comparative ROC curves (Saudi vs DDD vs Literature).
- `combined_pr_curves.png`: Comparative Precision-Recall curves.
- `saudi_sources_roc.png`: ROC curves for the 6 Saudi sources.
- `saudi_sources_pr.png`: PR curves for the 6 Saudi sources.
- `METHOD.md`: Detailed documentation of the computational methodology.
- `roc_curve_*.png`: Individual cohort ROC curves.

## Performance Metrics
- **AUC**: Area Under the ROC Curve, representing the probability that a true gene is ranked higher than a random one.
- **Hits@k**: Proportion of cases where the true gene is in the top k ranks.
- **MRR**: Mean Reciprocal Rank, a proxy for AUPR in ranking tasks.

## Results Summary
See `interpretation_full.txt` for detailed metrics and `roc_curve_*.png` for performance curves.

### Summary Observations
1. **DDD Performance**: The DDD cohort shows very high phenotypic alignment (Mean AUC: 0.9925, Hits@1: 55.86%), indicating high specificity in clinical descriptions or strong curation.
2. **Literature Performance**: Manually curated phenopackets from the public store perform exceptionally well (Mean AUC: 0.9797, Hits@1: 51.60%).
3. **Saudi Cohort Performance**: The Saudi cases show lower phenotypic alignment (Mean AUC: 0.8909, Hits@1: 2.81%). This suggests greater phenotypic heterogeneity, potential under-reporting in clinical text, or variants in genes that exhibit atypical presentations in the Saudi population.

The performance gap between curated literature/DDD and the Saudi cohort highlights the need for standardized phenotypic data collection in clinical settings.
