# About PAVS

## What is PAVS?

The **Phenotype-Associated Variants in Saudi Arabia (PAVS)** database is a curated resource
of genomic variants and associated clinical phenotypes from Saudi Arabian patients with
rare genetic diseases. It integrates data from multiple clinical cohorts and links each
case to international ontologies (HPO, OMIM, MONDO) and variant databases (ClinVar, dbSNP).

## Data Sources

- **PAVS Saudi cohort** — variants curated from clinical studies at KAUST (shown in green)
- **DDD study** — Deciphering Developmental Disorders, non-Saudi cases (shown in gray)
- **Literature phenopackets** — published case reports from the GA4GH phenopacket corpus (shown in indigo)

## Phenotype Similarity

Cases are searchable by phenotype using information-content-based similarity (**Lin + BMA**
by default, or **Resnik + BMA**). Similarity scores are computed using the Human Phenotype
Ontology (HPO) and disease–phenotype associations from the HPO Annotation file (phenotype.hpoa).

The **information content (IC)** of an HPO term *t* quantifies how specific it is,
based on its frequency across disease annotations:

> IC(*t*) = −log₂ P(*t*)

where P(*t*) is the proportion of cases annotated with *t* or any of its descendants.
Rare, specific terms have high IC; broad, common terms have low IC.

The **Lin similarity** between two HPO terms *a* and *b* is:

> Lin(*a*, *b*) = 2 × IC(MICA) / (IC(*a*) + IC(*b*))

where MICA is the Most Informative Common Ancestor of *a* and *b*.

The **Resnik similarity** uses only the MICA information content:

> Resnik(*a*, *b*) = IC(MICA)

**Best-Match Average (BMA)** symmetrises term-level scores across two phenotype sets *A* and *B*:

> BMA(*A*, *B*) = [ Σ_{a∈A} max_{b∈B} sim(*a*,*b*) + Σ_{b∈B} max_{a∈A} sim(*a*,*b*) ] / (|*A*| + |*B*|)

## Variant Annotation

Variants are annotated using the Ensembl Variant Effect Predictor (VEP) with:
- gnomAD allele frequencies
- SIFT and PolyPhen-2 pathogenicity scores
- ClinVar classifications
- Saudi cohort allele frequencies (novel to this resource)

Every variant with a known rsID or HGVS g. notation has a **TogoVar** link for instant
cross-referencing against Japanese and global variant databases.

## Citation

## Contact

**Robert Hoehndorf** (ORCID: [0000-0001-8149-5890](https://orcid.org/0000-0001-8149-5890))
King Abdullah University of Science and Technology (KAUST)

## Acknowledgements

This work has been supported by funding from King Abdullah University of Science and
Technology (KAUST) Office of Sponsored Research (OSR) under Award No. URF/1/5041-01-01,
REI/1/5235-01-01, REI/1/4938-01-01, and REI/1/5659-01-01.

This work was supported by funding from King Abdullah University of Science and Technology
(KAUST) — KAUST Center of Excellence for Smart Health (KCSH), award no. 5932, and by
funding from King Abdullah University of Science and Technology (KAUST) — Center of
Excellence for Generative AI, award no. 5940.

We acknowledge support from the KAUST Supercomputing Laboratory.

## License

Data are released under [Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/).

Code is released under the [GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0.html)
at [github.com/bio-ontology-research-group/pavs](https://github.com/bio-ontology-research-group/pavs).
