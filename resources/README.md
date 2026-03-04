# Resources: Clinical and Ontological Assets

This directory serves as a central repository for foundational clinical guidelines, research papers, and models used by the PAVS project.

## Foundational Guidelines & Papers

- **`nihms697486.pdf` / `nihms697486.txt`**: ACMG (American College of Medical Genetics and Genomics) 2015 variant interpretation guidelines. This is the global gold standard for standardizing clinical variant reporting.
- **`LEAP_paper.md`**: A summary of the LEAP (LLM-Enhanced Automated Phenotyping) paper published in the Journal of Genetics and Genomics (2026), which serves as a benchmark for our phenotype normalization tools.

## Models & Ontology Files

- **`en_core_sci_sm-0.5.4/`**: A spaCy model specialized for scientific and clinical text, used for entity recognition.
- **`hp.obo`**: The Human Phenotype Ontology (HPO) file used for phenotype standardization.
- **`phenotype_to_genes.txt`**: Mapping file connecting HPO phenotype terms to causative genes.

## External Tools

- **`phenopacket-tools/`**: Scripts and utilities for generating and validating GA4GH Phenopackets v2.0.
