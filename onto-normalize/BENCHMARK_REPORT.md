# Benchmark Report: phenotype_matcher_v2 on HPO Extraction Benchmarks

## Context

This benchmark evaluates **phenotype_matcher_v2** (the PAVS project's phenotype normalization tool) against published results from the field of clinical phenotype extraction. The benchmarks come from the RAG-HPO paper (doi:10.1186/s13073-025-01521-w), which are the same standard datasets used to evaluate tools including LEAP (Cai et al., JGG 2026, doi:10.1016/j.jgg.2026.02.009).

## Datasets

| Dataset | Cases | Gold HPO Terms | Source |
|---------|-------|---------------|--------|
| **CSC** (Case Study Cohort) | 112 | 1,789 | Published case reports from BMJ/Oxford Medical Case Reports |
| **GSC** (Gold Standard Corpora) | 114 | 1,012 | Curated OMIM disease descriptions |

## Configuration

- **phenotype_matcher_v2** settings: SapBERT embeddings, offline mode (no LLM validation), default thresholds (theta_abs=0.70, delta=0.05)
- **Evaluation**: Exact HPO ID matching (same methodology as RAG-HPO paper). Macro-averaged precision, recall, F1 across cases.

## Results: CSC (Case Study Cohort)

| Tool | Precision | Recall | F1 | Notes |
|------|-----------|--------|----|-------|
| RAG-HPO + LLaMa 3.1-70B | 0.8076 | 0.7595 | **0.7772** | Best overall |
| RAG-HPO + LLaMa 4-Scout | 0.6535 | 0.6890 | 0.6599 | |
| RAG-HPO + LLaMa 3-70B | 0.7096 | 0.6146 | 0.6435 | |
| **phenotype_matcher_v2** | **0.6678** | **0.4229** | **0.5071** | **This tool (no LLM)** |
| RAG-HPO + LLaMa 3.1-8B | 0.5338 | 0.4595 | 0.4782 | |
| FastHPOCR | 0.5296 | 0.4648 | 0.4781 | |
| Doc2HPO | 0.6936 | 0.3564 | 0.4561 | |
| ClinPhen | 0.6268 | 0.3532 | 0.4265 | |
| RAG-HPO + Mistral 24B | 0.7340 | 0.3128 | 0.4233 | |
| RAG-HPO + DeepSeek R1 | 0.3695 | 0.1086 | 0.1615 | |

## Results: GSC (Gold Standard Corpora)

| Tool | Precision | Recall | F1 | Notes |
|------|-----------|--------|----|-------|
| RAG-HPO + LLaMa 3-70B | 0.6904 | 0.7657 | **0.7082** | Best overall |
| FastHPOCR | 0.7357 | 0.6983 | 0.6998 | |
| RAG-HPO + LLaMa 4-Scout | 0.6238 | 0.7856 | 0.6793 | |
| **phenotype_matcher_v2** | **0.8053** | **0.4568** | **0.5649** | **This tool (no LLM)** |
| Doc2HPO | 0.7998 | 0.4351 | 0.5464 | |
| ClinPhen | 0.6097 | 0.4161 | 0.4613 | |

## Analysis

### Strengths of phenotype_matcher_v2

1. **High precision**: On GSC, achieves the **highest precision (0.8053)** of all tools tested, tied with Doc2HPO. On CSC, precision (0.6678) is competitive with RAG-HPO LLaMa 3-70B and above FastHPOCR, ClinPhen, and several RAG-HPO variants.

2. **No LLM dependency**: Unlike RAG-HPO (which requires a 70B parameter LLM) and LEAP (which requires OpenRouter API calls), phenotype_matcher_v2 runs fully offline using SapBERT embeddings + rule-based matching. This makes it fast, reproducible, and free.

3. **Beats established baselines**: Outperforms ClinPhen and Doc2HPO on both CSC and GSC by F1 score. Also beats FastHPOCR on CSC and is competitive on GSC.

4. **Beats small LLM variants**: Outperforms RAG-HPO + LLaMa 3.1-8B, RAG-HPO + Mistral 24B, and RAG-HPO + DeepSeek R1 on CSC.

### Weaknesses

1. **Lower recall**: The main gap vs. top-performing LLM tools is recall (0.42 CSC, 0.46 GSC). The tool matches fewer HPO terms per case than the gold standard. This is expected for a precision-oriented, offline tool.

2. **Gap to large LLMs**: RAG-HPO + LLaMa 3.1-70B achieves 0.78 F1 on CSC vs. our 0.51. The large LLM's ability to interpret clinical context gives it a significant advantage in recall.

### Key Insight

phenotype_matcher_v2 is a **high-precision, moderate-recall** tool. It predicts fewer HPO terms but those it predicts are usually correct. This is the right profile for:
- Automated pipeline use (where false positives are costly)
- Pre-populating annotations for human review
- Combining with other tools in an ensemble

### Potential Improvements

- **Enable LLM validation**: Running with `--no-llm` was used for reproducibility. Enabling LLM validation (OpenRouter) should improve recall by catching fuzzy/ambiguous matches.
- **Tune thresholds**: Lowering `theta_abs` from 0.70 to 0.60 would increase recall at some precision cost.
- **Ancestor-aware matching**: The current evaluation uses exact HPO ID matching. An ancestor-aware evaluation (where predicting a parent/child of the gold term counts as partial credit) would likely improve scores, as phenotype_matcher_v2 sometimes predicts a slightly more general/specific term.

## Runtime

| Dataset | Cases | Total Time | Per Case |
|---------|-------|-----------|----------|
| CSC | 112 | 4,842s (~81 min) | 43.2s |
| GSC | 114 | 2,457s (~41 min) | 21.5s |

Run on CPU (no GPU). GSC cases are shorter (OMIM descriptions) than CSC (full clinical reports), explaining the faster per-case time.

## Ranking Summary (by F1)

### CSC
1. RAG-HPO + LLaMa 3.1-70B (0.78)
2. RAG-HPO + LLaMa 4-Scout (0.66)
3. RAG-HPO + LLaMa 3-70B (0.64)
4. **phenotype_matcher_v2 (0.51)** -- best non-LLM tool
5. FastHPOCR (0.48)
6. RAG-HPO + LLaMa 3.1-8B (0.48)
7. Doc2HPO (0.46)
8. ClinPhen (0.43)

### GSC
1. RAG-HPO + LLaMa 3-70B (0.71)
2. FastHPOCR (0.70)
3. RAG-HPO + LLaMa 4-Scout (0.68)
4. **phenotype_matcher_v2 (0.56)** -- highest precision overall
5. Doc2HPO (0.55)
6. ClinPhen (0.46)

## Files

- `benchmark.py` — Benchmark script
- `csc_cases.json`, `gsc_cases.json` — Input clinical narratives
- `csc_annotations.json`, `gsc_annotations.json` — Gold standard HPO annotations
- `published_results.json` — Published tool results for comparison
- `benchmark_results.json` — Detailed per-case results
- `RAG-HPO_Tests_and_Data_Analysis.xlsx` — Original supplementary data
