# About PAVS

## What is PAVS?

The **Phenotypic and Variant Standardization (PAVS)** database is a curated resource
of genomic variants and associated clinical phenotypes from Saudi Arabian patients with
rare genetic diseases. It integrates data from multiple clinical cohorts and links each
case to international ontologies (HPO, OMIM, MONDO) and variant databases (ClinVar, dbSNP).

## Data Sources

- **PAVS Saudi cohort** — variants curated from clinical studies at KAUST (shown in green)
- **DDD study** — Deciphering Developmental Disorders, non-Saudi cases (shown in gray)
- **Literature phenopackets** — published case reports from the GA4GH phenopacket corpus (shown in indigo)

## Phenotype Similarity

Cases are searchable by phenotype using information-content-based similarity (**Lin + BMA**
by default, or Resnik + BMA). Similarity scores are computed using the Human Phenotype
Ontology (HPO) and disease–phenotype associations from the HPO Annotation file (phenotype.hpoa).

The Lin similarity between two HPO terms *a* and *b* is:

> Lin(a, b) = 2 × IC(MICA) / (IC(a) + IC(b))

where MICA is the Most Informative Common Ancestor. Best-Match Average (BMA) symmetrises
this over two phenotype sets.

## Variant Annotation

Variants are annotated using the Ensembl Variant Effect Predictor (VEP) with:
- gnomAD allele frequencies
- SIFT and PolyPhen-2 pathogenicity scores
- ClinVar classifications
- Saudi cohort allele frequencies (novel to this resource)

Every variant with a known rsID or HGVS g. notation has a **TogoVar** link for instant
cross-referencing against Japanese and global variant databases.

## Citation

If you use PAVS in your research, please cite:

> [Authors]. PAVS: Phenotypic and Variant Standardization for Saudi Arabian Rare
> Disease Cases. [Journal], [Year]. DOI: [TODO]

## Contact

**Robert Hoehndorf** (ORCID: [0000-0001-8149-5890](https://orcid.org/0000-0001-8149-5890))
King Abdullah University of Science and Technology (KAUST)

## Acknowledgements

[TODO: add funding sources, collaborators, and data contributors]

## License

Data and code are released under [TODO: specify license].
