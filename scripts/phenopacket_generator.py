# Description: This script reads the parsed data and generates a collection of
# GA4GH PhenoPacket JSON files using the pyphetools library.

import pandas as pd
from pyphetools.creation import Individual, MetaData, HpTerm, Disease, Citation, HgvsVariant
from phenopackets.schema.v2.core.interpretation_pb2 import VariantInterpretation
from pyphetools.validation import ContentValidator
from google.protobuf.json_format import MessageToJson
import os
import re
import argparse

def extract_gene_symbol(variant_string):
    """
    Extracts the gene symbol from a variant string.
    Common patterns:
    - GENE:NM_xxx:c.xxx
    - GENE:p.xxx
    - Just returns None if no clear gene symbol
    """
    if ':' in variant_string:
        parts = variant_string.split(':')
        # First part is often the gene symbol if it doesn't start with known prefixes
        first_part = parts[0].strip()
        if not any(first_part.startswith(prefix) for prefix in ['NM_', 'NR_', 'c.', 'g.', 'p.', 'chr', 'NC_']):
            # Likely a gene symbol
            return first_part
    return None

def parse_variant_components(variant_string):
    """
    Parses a variant string to extract gene symbol, transcript, and HGVS expressions.
    Returns a dict with 'symbol', 'transcript', 'genomic', 'protein', and 'other' components.
    """
    components = {
        'symbol': None,
        'transcript': None,
        'genomic': None,
        'protein': None,
        'cdna': None,
        'other': []
    }
    
    # Handle comma-separated expressions (same variant, different representations)
    expressions = [e.strip() for e in variant_string.split(',')]
    
    for expr in expressions:
        if ':' in expr:
            parts = expr.split(':', 2)
            
            # Check first part
            if len(parts) >= 1:
                first = parts[0].strip()
                # Gene symbol check
                if not any(first.startswith(prefix) for prefix in ['NM_', 'NR_', 'c.', 'g.', 'p.', 'chr', 'NC_']):
                    components['symbol'] = first
                elif first.startswith(('NM_', 'NR_')):
                    components['transcript'] = first
                    
            # Check second part if exists
            if len(parts) >= 2:
                second = parts[1].strip()
                if second.startswith(('NM_', 'NR_')):
                    components['transcript'] = second
                elif second.startswith('c.'):
                    components['cdna'] = second
                elif second.startswith('g.'):
                    components['genomic'] = second
                elif second.startswith('p.'):
                    components['protein'] = second
                    
            # Check third part if exists
            if len(parts) >= 3:
                third = parts[2].strip()
                if third.startswith('c.'):
                    components['cdna'] = third
                elif third.startswith('g.'):
                    components['genomic'] = third
                elif third.startswith('p.'):
                    components['protein'] = third
        else:
            # Single expression without colons
            if expr.startswith('p.'):
                components['protein'] = expr
            elif expr.startswith('g.'):
                components['genomic'] = expr
            elif expr.startswith('c.'):
                components['cdna'] = expr
            elif expr.startswith(('NM_', 'NR_')):
                components['transcript'] = expr
            else:
                components['other'].append(expr)
    
    return components

def create_phenopackets(parsed_data_path, output_dir):
    """
    Generates PhenoPacket JSON files from a structured DataFrame.
    """
    df = pd.read_csv(parsed_data_path).astype(str)
    df.columns = df.columns.str.strip()

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # --- Create Metadata ---
    # This is a critical step for FAIR compliance. It documents the resources used.
    # Replace with the actual versions you are using.
    created_by = "Marwa Abdelhakim"
    meta = MetaData(created_by=created_by)
    # Add HPO version used for annotation
    meta.hpo('2024-04-19')
    # Add MONDO for diseases (even if using OMIM, MONDO is often the bridge)
    meta.mondo('2024-04-01')

    # We will validate that each phenopacket has at least one HPO term.
    # Set min_hpo=0 to allow generation of phenopackets for individuals with no HPO terms.
    # Set min_disease=0 to allow generation for individuals with no disease annotation.
    validator = ContentValidator(min_hpo=0)
    
    print(f"Generating {len(df)} PhenoPackets...")

    for _, row in df.iterrows():
        individual_id = row['ID']

        # --- 1. Create HPO Terms ---
        hpo_terms = []
        if pd.notna(row['parsed_hpo_ids']) and row['parsed_hpo_ids']!= 'nan':
            hpo_ids = row['parsed_hpo_ids'].split(';')
            hpo_labels = []
            if pd.notna(row['parsed_pheno_text']) and row['parsed_pheno_text'] != 'nan':
                hpo_labels = row['parsed_pheno_text'].split(';')
            
            # If lengths match, we can pair them. Otherwise, we fall back to using IDs as labels.
            if len(hpo_ids) == len(hpo_labels):
                for hpo_id, label in zip(hpo_ids, hpo_labels):
                    term = HpTerm(hpo_id=hpo_id, label=label)
                    hpo_terms.append(term)
            else:
                for hpo_id in hpo_ids:
                    # Fallback for mismatched labels and IDs from the parser
                    term = HpTerm(hpo_id=hpo_id, label=hpo_id)
                    hpo_terms.append(term)

        # --- 2. Add Disease Diagnosis ---
        disease_obj = None
        if pd.notna(row['parsed_omim_ids']) and row['parsed_omim_ids']!= 'nan':
            omim_ids = row['parsed_omim_ids'].split(';')
            if omim_ids:
                # pyphetools Individual takes a single disease, so we use the first one.
                disease_obj = Disease(disease_id=omim_ids[0], disease_label=omim_ids[0])

        # --- 3. Add Genetic Interpretations ---
        variant_interpretations = []
        if pd.notna(row['parsed_variants']) and row['parsed_variants']!= 'nan':
            # Split by semicolon to get individual variants
            variant_strings = row['parsed_variants'].split(';')
            
            for var_string in variant_strings:
                var_string = var_string.strip()
                if not var_string:
                    continue
                    
                # Parse the variant components
                components = parse_variant_components(var_string)
                
                # Extract or infer the gene symbol
                symbol = components['symbol']
                transcript = components['transcript']
                
                # Build the HGVS expression
                # Prefer cDNA if we have transcript, otherwise use genomic
                hgvs_expr = None
                if transcript and components['cdna']:
                    hgvs_expr = f"{transcript}:{components['cdna']}"
                elif components['genomic']:
                    hgvs_expr = components['genomic']
                elif components['cdna']:
                    hgvs_expr = components['cdna']
                elif components['protein']:
                    # Protein-only variant
                    hgvs_expr = components['protein']
                else:
                    # Use the original string as fallback
                    hgvs_expr = var_string
                
                # We need to provide assembly and a placeholder vcf_d for HgvsVariant.
                # We will assume hg38, but this should be confirmed.
                assembly = 'hg38'
                vcf_d = {'chr': 'N/A', 'pos': 0, 'ref': 'N/A', 'alt': 'N/A'} # Placeholder
                
                # Create the variant with gene symbol
                variant = HgvsVariant(
                    assembly=assembly,
                    vcf_d=vcf_d,
                    symbol=symbol,  # This should add the gene symbol to the variant
                    transcript=transcript,
                    g_hgvs=hgvs_expr if components['genomic'] else None,
                    c_hgvs=components['cdna'] if components['cdna'] else None,
                    p_hgvs=components['protein'] if components['protein'] else None
                )
                
                # Get the GA4GH VariantInterpretation message
                variant_interpretation = variant.to_ga4gh_variant_interpretation()
                variant_interpretations.append(variant_interpretation)

        # --- 4. Add Provenance (Citation) ---
        citation = None
        if pd.notna(row['reference']) and row['reference'].startswith('http'):
            pubmed_id_match = re.search(r'(\d+)$', row['reference'].strip('/'))
            if pubmed_id_match:
                pmid = f"PMID:{pubmed_id_match.group(1)}"
                # The Citation object requires a title. We'll use the reference URL as a placeholder.
                citation = Citation(pmid=pmid, title=row['reference'])
        
        # --- 5. Assemble the Individual ---
        ind = Individual(
            individual_id=individual_id,
            sex="UNKNOWN",
            hpo_terms=hpo_terms,
            disease=disease_obj,
            interpretation_list=variant_interpretations,
            citation=citation
        )

        # --- 6. Create the Phenopacket from the Individual ---
        phenopacket = ind.to_ga4gh_phenopacket(
            metadata=meta,
            phenopacket_id=f"PAVS_{individual_id}"
        )

        # If no disease was found in the input, remove the default one added by pyphetools
        if disease_obj is None:
            del phenopacket.diseases[:]

        # --- 7. Validate and Save ---
        # It's good practice to validate each phenopacket.
        # The validator is created once outside the loop.
        errors = validator.validate_phenopacket(phenopacket)

        # The ContentValidator has an undocumented requirement for a disease.
        # We will filter out this specific error for individuals where we do not expect a disease annotation.
        if disease_obj is None:
            errors = [e for e in errors if "disease annotation" not in str(e)]

        if errors:
            print(f"Validation errors for {individual_id}: {errors}")
            continue

        # Save to JSON file
        json_string = MessageToJson(phenopacket)
        output_path = os.path.join(output_dir, f"PAVS_{individual_id}.json")
        with open(output_path, 'w') as f:
            f.write(json_string)

    print(f"Successfully generated {len(os.listdir(output_dir))} phenopackets in '{output_dir}'.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate GA4GH Phenopackets from a parsed data file.")
    parser.add_argument('parsed_data_path', help="Path to the input CSV file from data_parser.py.")
    parser.add_argument('--output_dir', default='phenopackets', help="Path to the output directory for PhenoPacket JSON files (default: 'phenopackets').")
    
    args = parser.parse_args()
    
    create_phenopackets(args.parsed_data_path, args.output_dir)
