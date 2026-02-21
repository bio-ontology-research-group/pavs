# HPO Arabic Translation Pipeline

This directory contains tools for translating the Human Phenotype Ontology (HPO) from English to Arabic using the OpenRouter API (GPT-4o).

## Tools

### 1. `translate_hpo_ar.py`
The core translation script.
- **Batching**: Processes multiple terms per API call to reduce costs.
- **Resume Support**: Automatically detects existing translations in the output JSON and skips them.
- **Strict Schema**: Uses JSON Schema enforcement to ensure consistent output.
- **Cost Tracking**: Displays real-time estimated cost in the progress bar.
- **Contextual Translation**: Provides the LLM with full definitions, parent terms, and all synonyms for high-accuracy medical translation.

**Usage:**
```bash
export OPENROUTER_API_KEY="your_key"
uv run --with aiohttp --with pronto --with tqdm python translation/translate_hpo_ar.py 
    --input ../ontology/hp.obo 
    --output hpo_arabic_translations.json 
    --batch-size 10
```

### 2. `export_results.py`
Converts the JSON output into distribution formats.
- **TSV Export**: Generates a tab-separated file for easy spreadsheet viewing.
- **OBO Generation**: Creates an updated `hp-ar.obo` file.
    - Arabic names are added as `synonym: "..." EXACT [PAVS:AR]`
    - Layperson terms are added as `synonym: "..." RELATED [PAVS:AR]`
    - Arabic definitions are added as `comment: Arabic Definition: ...`

**Usage:**
```bash
python3 translation/export_results.py --json hpo_arabic_translations.json --obo-in ../ontology/hp.obo
```

## Output Files
- `hpo_arabic_translations.json`: Raw LLM output.
- `hpo_arabic_translations.tsv`: Tabular version.
- `hp-ar.obo`: Standard OBO file enriched with Arabic language synonyms.
