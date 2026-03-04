# Bridging clinical narratives and structured phenotypes with large language models and sentence transformers

**Authors:** Jihao Cai, Guozhuang Li, Yongxin Yang, Kexin Xu, Sen Zhao, Timothy Hospedales, Lina Zhao, Jianle Yang, Zhihong Wu, Terry Jianguo Zhang, Zefu Chen, Nan Wu

**Journal:** Journal of Genetics and Genomics (J Genet Genomics)
**Published:** February 14, 2026 (Online ahead of print)
**PMID:** 41698530
**DOI:** 10.1016/j.jgg.2026.02.009
**URL:** https://www.sciencedirect.com/science/article/pii/S1673852726000536

**Keywords:** Electronic health records; Human phenotype ontology; Large language model; Mendelian disorders; Sentence transformers

## Abstract

The researchers address automated phenotyping from electronic health records by proposing LEAP (LLM-Enhanced Automated Phenotyping). Their two-stage framework combines language models for text extraction with a sentence-transformer fine-tuned on 5,330,557 instances for mapping to Human Phenotype Ontology identifiers. The method handles lengthy clinical narratives while producing valid HPO outputs. Testing showed relative improvements of 19.68%-412.68% in precision and 44.14%-298.77% in F1 score compared with existing tools.

## Tool Availability

- LEAP web interface: https://phenogemini.org/extract
- Part of the PhenoGemini platform for molecular diagnosis of Mendelian disorders

## Related Work

LEAP is the phenotype extraction component of PhenoGemini, which also includes:
- Gene prioritization via "phenotypic twins" matching
- Database of 382,474 individuals with 1,252,565 phenotype-genotype associations
- See: Research Square preprint rs-8387316/v1

## Benchmark Context

LEAP was compared against existing HPO extraction tools including Doc2HPO, ClinPhen, and FastHPOCR.
Standard benchmarks in this field include:
- Case Study Cohort (CSC): 112 published case reports, 1,794 HPO terms
- Gold Standard Corpora (GSC): 114 curated OMIM entries, 1,013 HPO terms
(As used in RAG-HPO evaluation: doi:10.1186/s13073-025-01521-w)
