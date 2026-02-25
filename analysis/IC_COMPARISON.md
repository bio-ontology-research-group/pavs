# Comparative Analysis: Intrinsic vs Extrinsic IC

## Intrinsic Resnik
- **DDD**: AUC=0.9929, MRR=0.6240, Hits@10=78.59%
- **Literature**: AUC=0.9814, MRR=0.6003, Hits@10=73.23%
- **Saudi**: AUC=0.8896, MRR=0.0585, Hits@10=11.24%

## Extrinsic (Gene-based) Resnik
- **DDD**: AUC=0.9951, MRR=0.7063, Hits@10=85.45%
- **Literature**: AUC=0.9828, MRR=0.6542, Hits@10=77.87%
- **Saudi**: AUC=0.8915, MRR=0.0676, Hits@10=12.47%

## Summary
Extrinsic Information Content (IC) calculated over the gene-to-phenotype distribution consistently outperforms intrinsic (topology-based) IC across all cohorts. The improvement is most pronounced in the **MRR (Mean Reciprocal Rank)**, which captures the precision of the top-ranked causative gene.
