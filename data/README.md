# Data Directory

## Structure

```
data/
├── phenotypes/          # Source TSV files — input to the normalization pipeline
│   ├── ahmed-variants.tsv           (A, 291 rows)
│   ├── ahmed-pmid28454995.tsv       (B, 234 rows)
│   ├── fawzan-variants.tsv          (F, 1,024 rows)
│   ├── marwa-variants.tsv           (M, 1,421 rows)
│   ├── PMC6562004.tsv               (P, 2,218 rows)
│   ├── PMC7082194.tsv               (Q, 522 rows)
│   └── ddd-diagnoses.tsv            (D, 1,856 rows — non-Saudi DDD cohort)
│
├── reference/           # External reference databases — download separately
│   ├── phenotype.hpoa               # HPO disease annotations
│   ├── genes_to_disease.txt         # Gene→disease associations (HPO)
│   ├── genes_to_phenotype.txt       # Gene→phenotype associations (HPO)
│   ├── Homo_sapiens.gene_info       # NCBI gene info
│   ├── mim2gene.txt                 # OMIM→Entrez gene mapping
│   ├── clinvar.vcf.gz               # ClinVar variant database (GRCh38)
│   ├── erepo.tabbed.txt             # ClinGen expert-curated variants
│   ├── goa_human.gaf.gz             # GO annotations for human genes
│   ├── gnomad.v4.1.constraint_metrics.tsv  # pLI / LOEUF constraint scores
│   ├── HMD_HumanPhenotype.rpt       # MGI human-mouse homologs + phenotypes
│   ├── MGI_GenePheno.rpt            # MGI genotype-phenotype associations
│   ├── GCF_000001405.40_GRCh38.p14_genomic.gff  # GRCh38 genome annotation
│   ├── E-GTEX-8/                    # GTEx v8 tissue expression (EBI Expression Atlas)
│   └── allele-frequencies/          # Saudi population allele frequencies (302 WGS)
│                                    # Publicly available: https://doi.org/10.6084/m9.figshare.28059686.v1
│                                    # Paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC11761371/
│
├── combined_normalized.tsv   # Generated: Step 1 output (38 columns, ~7,500 rows)
├── combined_annotated.tsv    # Generated: Step 2 output (77+ columns)
├── PAVS_phenopackets.json    # Generated: Step 3 output (all cases as JSON array)
├── PAVS_phenopackets.zip     # Generated: Step 3 output (download bundle)
│
├── annotation_cache/    # Local cache for VEP REST results (generated)
├── excel_files/         # Original Excel source backups
└── README.md            # This file
```

## Downloading Reference Data

Run the download script from the repository root:

```bash
bash scripts/download_reference_data.sh
```

This downloads all publicly available reference files to `data/reference/`.

The Saudi allele frequencies (`data/reference/allele-frequencies/`) are publicly available from Figshare (https://doi.org/10.6084/m9.figshare.28059686.v1). They were published as part of the paper PMC11761371 (https://pmc.ncbi.nlm.nih.gov/articles/PMC11761371/).

## Generated Files

The generated output files (`combined_normalized.tsv`, `combined_annotated.tsv`, `PAVS_phenopackets.json`, `PAVS_phenopackets.zip`) are produced by the pipeline:

```
normalization/combine_normalize_phenotypes.py → data/combined_normalized.tsv
normalization/annotate_variants.py            → data/combined_annotated.tsv
intake/generate_phenopackets_v2.py            → data/PAVS_phenopackets.json/.zip
```

See the main `README.md` and `docs/PIPELINE.md` for full pipeline instructions.
