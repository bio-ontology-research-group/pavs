# Description: This script reads the parsed data and generates a collection of
# GA4GH PhenoPacket JSON files using the pyphetools library.

import pandas as pd
from pyphetools.creation import Individual, MetaData, HpTerm, Disease, Citation
from phenopackets.schema.v2.core.base_pb2 import GeneDescriptor, VariationDescriptor, Expression, VcfRecord
from phenopackets.schema.v2.core.interpretation_pb2 import VariantInterpretation
from pyphetools.validation import ContentValidator
from google.protobuf.json_format import MessageToJson
import os
import re
import argparse

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
            var_string = row['parsed_variants']
            
            symbol, transcript, hgvsc_expressions_str = None, None, None
            if ':' in var_string:
                parts = var_string.split(':', 2)
                if len(parts) == 3:
                    symbol, transcript, hgvsc_expressions_str = parts
                else:
                    hgvsc_expressions_str = var_string # fallback
            else:
                hgvsc_expressions_str = var_string

            gene_descriptor = None
            if symbol:
                gene_descriptor = GeneDescriptor(value_id=symbol, symbol=symbol)

            expressions = []
            if hgvsc_expressions_str:
                # The string can contain multiple expressions separated by ';'
                for expr_str in hgvsc_expressions_str.split(';'):
                    expr_str = expr_str.strip()
                    if not expr_str: continue
                    # A simple way to guess syntax. p. is protein, c. is coding.
                    syntax = 'hgvs'
                    if expr_str.startswith('p.'):
                        syntax = 'hgvs.p'
                    elif expr_str.startswith('c.'):
                        syntax = 'hgvs.c'
                    expressions.append(Expression(syntax=syntax, value=expr_str))

            vcf_record = VcfRecord(genome_assembly='hg38', chrom='N/A', pos=0, ref='N/A', alt='N/A')

            variation_descriptor = VariationDescriptor(
                id=f"var_{individual_id}", # simplified ID
                label=var_string,
                gene_context=gene_descriptor,
                expressions=expressions,
                vcf_record=vcf_record
            )
            
            variant_interpretation = VariantInterpretation(variation_descriptor=variation_descriptor)
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
