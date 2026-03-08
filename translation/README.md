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

## Translation Conventions

### 1. Arabic Sentence Structure
For HPO terms, always state the abnormality/problem first, followed by the anatomical location.
- **Rule**: `[Abnormality] + في + [Location]`
- **Examples**:
    - *Fingernail hypoplasia* → نقص تنسج في أظافر اليد
    - *Renal tubular dysfunction* → خلل في الأنابيب الكلوية
    - *Cardiac rhythm abnormality* → اضطراب في نظم القلب

**Exception**: When using an adjective like "Metaphyseal" (كردوسي):
- *Metaphyseal dysplasia* → خلل تنسجي كردوسي

### 2. Medical Glossary
| English | Arabic |
| :--- | :--- |
| Metaphysis | كردوس |
| Metaphyseal | كردوسي |
| Epiphysis | مشاشة |
| Dysplasia | خلل تنسجي |
| Hyperplasia | فرط تنسج |
| Hypertrophy | تضخم |
| Atresia | رتق |
| Malacia | تلين |
| -pathy | اعتلال |
| Dystrophy | حثل |
| Atrophy | ضمور |

### 3. Abnormality Translation Rules
- **Physical Structure**: Use **شذوذ** (e.g., bone shape).
- **Functional Deficit**: Use **خلل** (e.g., enzyme activity).
- **Complex System/Pattern**: Use **اضطراب** (e.g., heart rhythm).
- **Morphogenesis**:
    - *Hypoplasia* → نقص تنسج في ...
    - *Aplasia* → عدم التنسج في ...
    - *Agenesis* → عدم التكوّن في ...
