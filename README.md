# Data Parsing Pipeline

This repository contains scripts for parsing and structuring clinical data from a CSV file.

## Description

The main script, `scripts/data_parser.py`, reads a manually curated data file (e.g., tab-separated CSV) and processes the `phenotypes` and `variants` columns. These columns contain complex, semi-structured text with multiple entries.

The script performs the following actions:
- Parses HPO terms, OMIM IDs, and free-text phenotype descriptions from the `phenotypes` column.
- Parses HGVS strings and other variant notations from the `variants` column.
- Creates new columns in the dataset (`parsed_hpo_ids`, `parsed_omim_ids`, `parsed_pheno_text`, `parsed_variants`) to store the structured data.
- Saves the processed data to a new CSV file.

## Installation

To install the required Python libraries, run the following command:

```bash
pip install -r requirements.txt
```

## Usage

To run the data parser script, you must provide the path to the input file.

```bash
python scripts/data_parser.py <path_to_input_file.csv>
```

For example:
```bash
python scripts/data_parser.py "data/Manually curated Data.csv"
```

You can also specify an output file path using the `--output_csv_path` argument. If not provided, it will default to `parsed_data.csv` in the current directory.

```bash
python scripts/data_parser.py "data/Manually curated Data.csv" --output_csv_path "output/structured_data.csv"
```

For more information on command-line arguments, use the `--help` flag:
```bash
python scripts/data_parser.py --help
```
