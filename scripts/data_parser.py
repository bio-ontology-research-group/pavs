# SCRIPT 1: data_parser.py
# Description: This script reads the manually curated CSV file, parses the complex
# 'phenotypes' and 'variants' columns, and outputs a structured pandas DataFrame.

import pandas as pd
import re
import hgvs.parser
import logging

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
hpo_label_pattern = re.compile(r'([\w\s-]+)\s*\((HP:\d{7})\)')

def parse_phenotypes(phenotype_string):
    """
    Parses a string from the 'phenotypes' column to extract HPO terms, OMIM IDs,
    and any remaining free-text descriptions.
    """
    if not isinstance(phenotype_string, str):
        return [], [], []

    # Normalize delimiters
    phenotype_string = phenotype_string.replace('|', ';').replace(',', ';')
    terms = [term.strip() for term in phenotype_string.split(';') if term.strip()]

    hpo_ids = set()
    omim_ids = set()
    free_text = set()

    for term in terms:
        # Find all HPO IDs in the term
        hpo_matches = hpo_pattern.findall(term)
        if hpo_matches:
            hpo_ids.update(hpo_matches)

        # Find all OMIM IDs in the term
        omim_matches = omim_pattern.findall(term)
        if omim_matches:
            omim_ids.update(omim_matches)

        # If no standard identifiers are found, treat it as free text
        if not hpo_matches and not omim_matches:
            # Clean up the term by removing any parenthetical content without IDs
            cleaned_term = re.sub(r'\([^)]*\)', '', term).strip()
            if cleaned_term:
                free_text.add(cleaned_term)
                
    return sorted(list(hpo_ids)), sorted(list(omim_ids)), sorted(list(free_text))


def parse_variants(variant_string):
    """
    Parses a string from the 'variants' column to extract HGVS strings and
    other variant descriptions.
    """
    if not isinstance(variant_string, str):
        return []

    # Initialize the HGVS parser
    hp = hgvs.parser.Parser()
    
    # Normalize delimiters and split into individual variant strings
    variant_string = variant_string.replace('|', ';').replace(',', ';')
    variants = [v.strip() for v in variant_string.split(';') if v.strip()]
    
    parsed_variants = []
    for var in variants:
        try:
            # Attempt to parse as a standard HGVS string
            parsed_var = hp.parse_hgvs_variant(var)
            # Successfully parsed, store as a standardized string
            parsed_variants.append(str(parsed_var))
        except hgvs.exceptions.HGVSParseError:
            # If parsing fails, store the original string as a description
            # This handles non-standard formats like chromosomal coordinates
            parsed_variants.append(var)
            
    return parsed_variants


def main(input_csv_path, output_csv_path):
    """
    Main function to execute the parsing pipeline.
    """
    print(f"Reading data from {input_csv_path}...")
    df = pd.read_csv(input_csv_path)
    
    # Create new columns for the parsed data
    df['parsed_hpo_ids'] = ''
    df['parsed_omim_ids'] = ''
    df['parsed_pheno_text'] = ''
    df['parsed_variants'] = ''

    print("Parsing phenotypes and variants for each row...")
    for index, row in df.iterrows():
        # Parse phenotypes
        hpo_ids, omim_ids, free_text = parse_phenotypes(row['phenotypes'])
        df.at[index, 'parsed_hpo_ids'] = ';'.join(hpo_ids)
        df.at[index, 'parsed_omim_ids'] = ';'.join(omim_ids)
        df.at[index, 'parsed_pheno_text'] = ';'.join(free_text)
        
        # Parse variants
        variants = parse_variants(row['variants'])
        df.at[index, 'parsed_variants'] = ';'.join(variants)

    # Select and reorder columns for the output file
    output_df = df
    
    print(f"Saving structured data to {output_csv_path}...")
    output_df.to_csv(output_csv_path, index=False)
    print("Parsing complete.")

if __name__ == '__main__':
    # Define file paths
    INPUT_FILE = 'a.csv'
    OUTPUT_FILE = 'parsed_saudi_variants.csv'
    
    # Run the pipeline
    main(INPUT_FILE, OUTPUT_FILE)
