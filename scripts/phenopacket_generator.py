# Description: This script reads the parsed data and generates a collection of
# GA4GH PhenoPacket JSON files using the pyphetools library.

import pandas as pd
from pyphetools.creation import Individual, MetaData, HpTerm, Disease, Citation, HgvsVariant
from phenopackets.schema.v2.phenopackets_pb2 import (
    ExternalReference,
)
from phenopackets.schema.v2.core.interpretation_pb2 import VariantInterpretation
from pyphetools.validation import PhenopacketValidator
from google.protobuf.json_format import MessageToJson
import os
import re
import argparse

def create_phenopackets(parsed_data_path, output_dir):
    """
    Generates PhenoPacket JSON files from a structured DataFrame.
    """
    df = pd.read_csv(parsed_data_path).astype(str)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # --- Create Metadata ---
    # This is a critical step for FAIR compliance. It documents the resources used.
    # Replace with the actual versions you are using.
    created_by = "PAVS_Curation_Team"
    meta = MetaData(created_by=created_by)
    meta.add_default_metadata() # Adds basic creation timestamp
    # Add HPO version used for annotation
    meta.add_ontology('hp', 'Human Phenotype Ontology', '2024-04-19', 'http://purl.obolibrary.org/obo/hp.owl')
    # Add MONDO for diseases (even if using OMIM, MONDO is often the bridge)
    meta.add_ontology('mondo', 'Monarch Disease Ontology', '2024-04-01', 'http://purl.obolibrary.org/obo/mondo.obo')
    
    print(f"Generating {len(df)} PhenoPackets...")

    for _, row in df.iterrows():
        individual_id = row['ID']

        # --- 1. Create HPO Terms ---
        hpo_terms = []
        if pd.notna(row['parsed_hpo_ids']) and row['parsed_hpo_ids']!= 'nan':
            hpo_ids = row['parsed_hpo_ids'].split(';')
            for hpo_id in hpo_ids:
                # pyphetools requires a label. We will use the ID as the label for now.
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
            variants = row['parsed_variants'].split(';')
            for var_string in variants:
                # We need to provide assembly and a placeholder vcf_d for HgvsVariant.
                # We will assume hg38, but this should be confirmed.
                assembly = 'hg38'
                vcf_d = {'chr': 'N/A', 'pos': 0, 'ref': 'N/A', 'alt': 'N/A'} # Placeholder

                # HgvsVariant needs a transcript and hgvsc notation.
                # We will create a simple HgvsVariant object.
                # A more robust solution might use a library like 'hgvs' to parse the string.
                if ':' in var_string:
                    parts = var_string.split(':', 2)
                    if len(parts) == 3:
                        # Assumes format like GENE:TRANSCRIPT:c.CHANGE
                        symbol, transcript, hgvsc = parts
                        variant = HgvsVariant(assembly=assembly, vcf_d=vcf_d, symbol=symbol, transcript=transcript, hgvsc=hgvsc)
                    else:
                        # Fallback for other formats
                        variant = HgvsVariant(assembly=assembly, vcf_d=vcf_d, g_hgvs=var_string)
                else:
                    # Handle non-HGVS variants if necessary, or just label them
                    variant = HgvsVariant(assembly=assembly, vcf_d=vcf_d, g_hgvs=var_string)

                # Get the GA4GH VariantInterpretation message
                variant_interpretation = variant.to_ga4gh_variant_interpretation()
                variant_interpretations.append(variant_interpretation)

        # --- 4. Add Provenance (Citation and External Reference) ---
        citation = None
        external_references = []
        if pd.notna(row['reference']) and row['reference'].startswith('http'):
            pubmed_id_match = re.search(r'(\d+)$', row['reference'].strip('/'))
            if pubmed_id_match:
                pmid = pubmed_id_match.group(1)
                # For pyphetools Individual object
                citation = Citation(pmid=pmid)
                # For GA4GH Phenopacket object
                ext_ref = ExternalReference(id=f"PMID:{pmid}", reference=row['reference'])
                external_references.append(ext_ref)
        
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

        # Add external references to the phenopacket
        if external_references:
            phenopacket.external_references.extend(external_references)

        # --- 7. Validate and Save ---
        # It's good practice to validate each phenopacket
        validator = PhenopacketValidator(phenopacket_string=MessageToJson(phenopacket))
        errors = validator.validate()
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
