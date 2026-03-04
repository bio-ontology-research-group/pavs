# Onto-Normalize: Phenotype Extraction Benchmarks

This directory contains the benchmarking infrastructure and results for evaluating the PAVS project's phenotype normalization tool, **phenotype_matcher_v2**, against industry-standard datasets and state-of-the-art tools.

## Contents

- **`benchmark.py`**: The main execution script for running `phenotype_matcher_v2` on the CSC and GSC datasets.
- **`BENCHMARK_REPORT.md`**: A detailed report comparing our results with published benchmarks from RAG-HPO and LEAP (Cai et al., 2026).
- **`benchmark_results.json`**: Raw output from the benchmark runs.
- **`csc_cases.json` / `csc_annotations.json`**: Input narratives and gold-standard HPO annotations for the Case Study Cohort (112 cases).
- **`gsc_cases.json` / `gsc_annotations.json`**: Input narratives and gold-standard HPO annotations for the Gold Standard Corpora (114 OMIM cases).
- **`published_results.json`**: Aggregated performance metrics (Precision, Recall, F1) for competing tools (ClinPhen, Doc2HPO, RAG-HPO, etc.).
- **`HPO_addons.csv`**: Supplementary HPO mappings used to tune the normalization process.
- **`tune_ner.py`**: Utility for optimizing Named Entity Recognition (NER) thresholds.

## Key Performance (Macro-F1)

| Dataset | phenotype_matcher_v2 (PAVS) | Next Best Non-LLM Tool | State-of-the-art (LLM) |
|---------|-----------------------------|------------------------|-----------------------|
| **CSC** | **0.51**                    | 0.48 (FastHPOCR)       | 0.78 (RAG-HPO 70B)    |
| **GSC** | **0.56**                    | 0.55 (Doc2HPO)         | 0.71 (RAG-HPO 70B)    |

*Note: `phenotype_matcher_v2` achieves the highest precision (0.81) on the GSC dataset among all tested tools.*

## Usage

To reproduce the benchmarks, ensure `phenotype_matcher_v2` is installed in your environment, then run:

```bash
uv run python benchmark.py
```
