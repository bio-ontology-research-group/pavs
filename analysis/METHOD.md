# Methodology: Phenotype Semantic Similarity Analysis

This document describes the computational approach used to evaluate the phenotypic alignment between clinical cases and their causative genes.

## 1. Data Sources
- **Ontology**: Human Phenotype Ontology (HPO), `hp.obo` (2026-01-08).
- **Gene-Phenotype Associations**: `genes_to_phenotype.txt`, providing the reference set of HPO terms for each gene (~5,200 genes).
- **Cohorts**:
  - **Saudi**: Clinical data from Middle Eastern datasets.
  - **DDD**: Deciphering Developmental Disorders cohort data.
  - **Literature**: Manually curated Phenopackets from the public store (v0.1.26).

## 2. Information Content (IC) Calculation
We implement two distinct methods for computing Information Content (IC):

### 2.1 Intrinsic Resnik IC
This method relies solely on the topology of the HPO directed acyclic graph (DAG).
The IC of a term $t$ is defined as:
$IC_{intrinsic}(t) = -\ln\left(\frac{descendants(t)}{total\_terms}\right)$
where $descendants(t)$ is the number of terms in the sub-ontology rooted at $t$, and $total\_terms$ is the total number of terms in the HPO.

### 2.2 Extrinsic (Gene-based) IC
Since the task is to rank **genes**, this method computes IC over the distribution of phenotypes in the gene-to-phenotype reference set. This reflects the informativeness of a term relative to the known disease-gene landscape.
For each term $t$, the frequency $f(t)$ is the number of genes associated with $t$ or any of its descendants.
The IC is defined as:
$IC_{extrinsic}(t) = -\ln\left(\frac{f(t)}{Total\_Genes}\right)$
where $Total\_Genes$ is the total number of unique genes in the reference database (5,193).

## 3. Pairwise Similarity: Lin's Measure

The similarity between two individual HPO terms $t_1$ and $t_2$ is calculated using **Lin's similarity**, which normalizes the IC of the Most Informative Common Ancestor (MICA) by the IC of the individual terms:
$Sim_{Lin}(t_1, t_2) = \frac{2 \times IC(MICA(t_1, t_2))}{IC(t_1) + IC(t_2)}$
where $MICA(t_1, t_2)$ is the common ancestor of $t_1$ and $t_2$ with the highest IC value.

## 4. Group Similarity: Best Match Average (BMA)
To compare a set of clinical phenotypes $P$ (from a patient) with a set of gene-associated phenotypes $G$ (from the reference), we use the **Best Match Average (BMA)** aggregation:
$Sim_{BMA}(P, G) = \frac{1}{2} \left( \frac{\sum_{p \in P} \max_{g \in G} Sim(p, g)}{|P|} + \frac{\sum_{g \in G} \max_{p \in P} Sim(g, p)}{|G|} \right)$

This symmetric measure captures how well each clinical term is "explained" by the gene's known phenotypes and vice versa.

## 5. Evaluation Framework
For each case:
1. The patient's HPO set is extracted.
2. The BMA similarity is calculated against the HPO sets of **all 5,193 genes** in the reference database.
3. All genes are ranked by their similarity scores (descending).
4. The rank of the **actual causative gene** is recorded.

### Performance Metrics
- **Mean AUC (Area Under ROC)**: Probability that the true gene is ranked higher than a random gene.
- **Hits@k**: Percentage of cases where the true gene is in the top $k$ positions (e.g., $k=1, 10, 50$).
- **Mean Reciprocal Rank (MRR)**: Average of $1/rank$ across all cases, providing a proxy for Area Under the Precision-Recall (AUPR) curve.

## 6. Sub-analysis by Source
The Saudi cohort is further subdivided into its six original source datasets (e.g., Alkuraya, Marwa). The same evaluation methodology described above is applied independently to each sub-cohort to assess variability in phenotypic data quality and specificity across different studies.
