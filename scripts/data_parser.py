# SCRIPT 1: data_parser.py
# Description: This script reads the manually curated CSV file, parses the complex
# 'phenotypes' and 'variants' columns, and outputs a structured pandas DataFrame.

import pandas as pd
import re
import hgvs.parser
import logging
import argparse
import os
from tqdm import tqdm

# --- Setup Logging ---
# Suppress verbose logging from the hgvs library to keep output clean
logging.basicConfig()
logging.getLogger("hgvs").setLevel(logging.WARNING)

# --- Regular Expressions for Parsing ---
# Regex to find HPO terms, e.g., (HP:0000510)
hpo_pattern = re.compile(r'\((HP:\d{7})\)')
# Regex to find OMIM disease terms, e.g., (OMIM:204000)
omim_pattern = re.compile(r'\((OMIM:\d{6})\)')
# Regex to find free-text terms that are followed by an HPO ID
# This captures the label associated with an ID.
hpo_label_pattern = re.compile(r'([\w\s;-]+)\s*\((HP:\d{7})\)')
omim_label_pattern = re.compile(r'([\w\s;-]+)\s*\((OMIM:\d{6})\)')

# --- HGVS Parser ---
# Initialize the parser once to be reused, as initialization is expensive.
hgvs_parser = hgvs.parser.Parser()


def build_hpo_mapper(obo_path):
    """
    Parses the hp.obo file to build a mapping from term labels and synonyms
    to HPO IDs. Matching is case-insensitive.
    """
    if not os.path.exists(obo_path):
        logging.error(f"HPO OBO file not found at: {obo_path}")
        raise FileNotFoundError(f"HPO OBO file not found at: {obo_path}")

    mapper = {}
    with open(obo_path, 'r', encoding='utf-8') as f:
        current_id = None
        for line in f:
            line = line.strip()
            if line == "[Term]":
                current_id = None
            elif line.startswith('id:'):
                current_id = line.split('id: ')[1]
            elif line.startswith('name:'):
                if current_id:
                    name = line.split('name: ')[1]
                    mapper[name.lower()] = current_id
            elif line.startswith('synonym:'):
                if current_id:
                    match = re.search(r'"([^"]+)"', line)
                    if match:
                        synonym = match.group(1)
                        mapper[synonym.lower()] = current_id
    
    logging.info(f"Built HPO mapper with {len(mapper)} labels/synonyms from {obo_path}.")
    return mapper


def parse_phenotypes(phenotype_string, hpo_mapper=None):
    """
    Parses a string from the 'phenotypes' column to extract HPO terms, OMIM IDs,
    and any remaining free-text descriptions. If an hpo_mapper is provided,
    it attempts to map free text to HPO IDs.
    """
    if not isinstance(phenotype_string, str):
        return [], [], []

    # Normalize delimiters, treating '|' and ',' as separators.
    phenotype_string = phenotype_string.replace(',', '|')
    terms = [term.strip() for term in phenotype_string.split('|') if term.strip()]

    hpo_ids = set()
    omim_ids = set()
    unmapped_text = set()

    for term in terms:
        hpo_label_match = hpo_label_pattern.search(term)
        omim_label_match = omim_label_pattern.search(term)

        # Case 1: Term has an explicit HPO or OMIM ID, e.g., "Seizures (HP:0001250)"
        if hpo_label_match:
            _, hpo_id = hpo_label_match.groups()
            hpo_ids.add(hpo_id)
        elif omim_label_match:
            _, omim_id = omim_label_match.groups()
            omim_ids.add(omim_id)
        else:
            # Case 2: Term might be just an ID, e.g., "(HP:0001250)"
            hpo_matches = hpo_pattern.findall(term)
            omim_matches = omim_pattern.findall(term)
            
            hpo_ids.update(hpo_matches)
            omim_ids.update(omim_matches)
            
            # Case 3: Term is free text, e.g., "Seizures"
            # Clean out any parenthesized content that might have been missed
            cleaned_term = re.sub(r'\([^)]*\)', '', term).strip()
            if cleaned_term and not hpo_matches and not omim_matches:
                # Attempt to map the free text to an HPO ID
                if hpo_mapper:
                    mapped_id = hpo_mapper.get(cleaned_term.lower())
                    if mapped_id:
                        hpo_ids.add(mapped_id)
                    else:
                        # If mapping fails, add to unmapped text
                        unmapped_text.add(cleaned_term)
                else:
                    # If no mapper, it's just unmapped text
                    unmapped_text.add(cleaned_term)
                
    return sorted(list(hpo_ids)), sorted(list(omim_ids)), sorted(list(unmapped_text))


def split_variant_expressions(variant_string):
    """
    Splits a variant string that may contain multiple HGVS expressions.
    Common patterns include:
    - Gene:transcript:genomic_change:protein_change
    - Multiple expressions separated by commas or semicolons
    - Mixed genomic (g.) and protein (p.) coordinates
    """
    expressions = []
    
    # First, check if it's a compound expression with colons
    if ':' in variant_string:
        parts = variant_string.split(':')
        
        # Try to identify gene symbol (usually first part if it doesn't start with NM_ or c. or g. or p.)
        gene_symbol = None
        transcript = None
        
        for i, part in enumerate(parts):
            part = part.strip()
            
            # Gene symbol is typically uppercase letters/numbers, not starting with known prefixes
            if i == 0 and not any(part.startswith(prefix) for prefix in ['NM_', 'NR_', 'c.', 'g.', 'p.', 'chr']):
                gene_symbol = part
            # Transcript usually starts with NM_ or NR_
            elif part.startswith(('NM_', 'NR_')):
                transcript = part
            # HGVS expressions start with c., g., p., etc.
            elif any(part.startswith(prefix) for prefix in ['c.', 'g.', 'p.', 'n.', 'r.']):
                # Construct full HGVS expression
                if transcript and part.startswith(('c.', 'n.', 'r.')):
                    expressions.append(f"{transcript}:{part}")
                elif gene_symbol and part.startswith('p.'):
                    expressions.append(f"{gene_symbol}:{part}")
                else:
                    expressions.append(part)
    
    # Also check for expressions separated by commas or semicolons
    if not expressions:
        # Split by common delimiters
        potential_variants = re.split(r'[,;]\s*', variant_string)
        for var in potential_variants:
            var = var.strip()
            if var:
                expressions.append(var)
    
    # If still no expressions, return the original
    if not expressions:
        expressions = [variant_string]
    
    return expressions


def parse_variants(variant_string):
    """
    Parses a string from the 'variants' column to extract HGVS strings and
    other variant descriptions. Handles multiple HGVS expressions for the same variant.
    """
    if not isinstance(variant_string, str):
        return []

    # Normalize delimiters and split into individual variant strings
    variant_string = variant_string.replace('|', ';')
    variants = [v.strip() for v in variant_string.split(';') if v.strip()]
    
    all_parsed_variants = []
    
    for var in variants:
        # Split potential multiple expressions
        expressions = split_variant_expressions(var)
        
        parsed_expressions = []
        for expr in expressions:
            try:
                # Attempt to parse as a standard HGVS string
                parsed_var = hgvs_parser.parse_hgvs_variant(expr)
                # Successfully parsed, store as a standardized string
                parsed_expressions.append(str(parsed_var))
            except hgvs.exceptions.HGVSParseError:
                # If parsing fails, store the original string as a description
                # This handles non-standard formats like chromosomal coordinates
                parsed_expressions.append(expr)
        
        # Join related expressions with comma (they represent the same variant)
        if parsed_expressions:
            all_parsed_variants.append(','.join(parsed_expressions))
            
    return all_parsed_variants


def main(input_csv_path, output_csv_path, hpo_obo_path=None):
    """
    Main function to execute the parsing pipeline.
    """
    hpo_mapper = None
    if hpo_obo_path:
        print(f"Building HPO mapper from {hpo_obo_path}...")
        hpo_mapper = build_hpo_mapper(hpo_obo_path)

    print(f"Reading data from {input_csv_path}...")
    df = pd.read_csv(input_csv_path, sep='\t')

    print(f"Found {len(df)} rows to process.")
    
    # Create new columns for the parsed data
    df['parsed_hpo_ids'] = ''
    df['parsed_omim_ids'] = ''
    df['parsed_pheno_text'] = ''
    df['parsed_variants'] = ''

    print("Parsing phenotypes and variants for each row...")
    for index, row in tqdm(df.iterrows(), total=df.shape[0], desc="Parsing rows"):
        # Parse phenotypes
        hpo_ids, omim_ids, unmapped_text = parse_phenotypes(row['phenotypes'], hpo_mapper=hpo_mapper)
        df.at[index, 'parsed_hpo_ids'] = ';'.join(hpo_ids)
        df.at[index, 'parsed_omim_ids'] = ';'.join(omim_ids)
        df.at[index, 'parsed_pheno_text'] = ';'.join(unmapped_text)
        
        # Parse variants
        variants = parse_variants(row['variants'])
        df.at[index, 'parsed_variants'] = ';'.join(variants)

    # Select and reorder columns for the output file
    output_df = df
    
    print(f"Saving structured data to {output_csv_path}...")
    output_df.to_csv(output_csv_path, index=False)
    print("Parsing complete.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Reads a CSV file, parses 'phenotypes' and 'variants' columns, and outputs a structured CSV."
    )
    parser.add_argument(
        "input_csv_path",
        help="Path to the input CSV file."
    )
    parser.add_argument(
        "--output_csv_path",
        default="parsed_data.csv",
        help="Path for the output structured CSV file (default: parsed_data.csv)."
    )
    parser.add_argument(
        "--hpo_obo_path",
        default=None,
        help="Path to the hp.obo file to enable mapping of free-text phenotypes to HPO IDs."
    )
    args = parser.parse_args()

    # Run the pipeline
    main(args.input_csv_path, args.output_csv_path, args.hpo_obo_path)
