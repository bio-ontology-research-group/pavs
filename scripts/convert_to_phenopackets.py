#!/usr/bin/env python3
"""
Convert PAVS final data to Phenopackets v2 format (JSON).
Enhanced version with:
- HPO term labels from ontology
- Gene symbol extraction from variants
- OMIM ID extraction from diagnosis
- Improved age parsing
- Better variant parsing
"""

import pandas as pd
import json
from datetime import datetime
import argparse
import pronto
import re


def load_hpo_ontology(hpo_path='resources/hp.obo'):
    """
    Load HPO ontology and create a mapping of HPO ID -> label.
    """
    print(f"Loading HPO ontology from {hpo_path}...")
    
    try:
        hpo = pronto.Ontology(hpo_path)
        hpo_labels = {}
        
        for term in hpo.terms():
            if term.id.startswith('HP:'):
                hpo_labels[term.id] = term.name
        
        print(f"✅ Loaded {len(hpo_labels):,} HPO terms with labels")
        return hpo_labels
        
    except Exception as e:
        print(f"⚠️  Could not load HPO ontology: {e}")
        print("   Continuing without HPO labels...")
        return {}


def parse_age(age_str):
    """
    Enhanced age parsing with support for multiple formats.
    
    Examples:
        '22Y' -> 'P22Y'
        '6 months' -> 'P6M'
        '3 years 6 months' -> 'P3Y6M'
        'infant' -> 'P6M'
        'newborn' -> 'P0D'
    
    Returns:
        ISO8601 duration string or None
    """
    if pd.isna(age_str) or age_str == 'Not reported' or age_str == '':
        return None
    
    age_str = str(age_str).strip().upper()
    
    # Already in ISO format
    if age_str.startswith('P'):
        return age_str
    
    # Handle text descriptions
    text_ages = {
        'NEWBORN': 'P0D',
        'NEONATE': 'P0D',
        'NEONATAL': 'P0D',
        'INFANT': 'P6M',
        'CHILD': 'P5Y',
        'CHILDHOOD': 'P5Y',
        'ADOLESCENT': 'P15Y',
        'ADOLESCENCE': 'P15Y',
        'ADULT': 'P25Y',
        'ADULTHOOD': 'P25Y',
    }
    
    if age_str in text_ages:
        return text_ages[age_str]
    
    # Parse numeric values with regex
    iso_age = 'P'
    
    # Years - handles "22Y", "22 years", "22 y", "22 year"
    years_match = re.search(r'(\d+)\s*(?:YEARS?|Y)', age_str)
    if years_match:
        iso_age += f"{years_match.group(1)}Y"
    
    # Months - handles "6M", "6 months", "6 month"
    months_match = re.search(r'(\d+)\s*(?:MONTHS?|M)(?!Y)', age_str)  # Negative lookahead to avoid matching MY
    if months_match:
        iso_age += f"{months_match.group(1)}M"
    
    # Days - handles "15D", "15 days", "15 day"
    days_match = re.search(r'(\d+)\s*(?:DAYS?|D)', age_str)
    if days_match:
        iso_age += f"{days_match.group(1)}D"
    
    # Weeks - convert to days (7 days per week)
    weeks_match = re.search(r'(\d+)\s*(?:WEEKS?|W)', age_str)
    if weeks_match:
        days = int(weeks_match.group(1)) * 7
        iso_age += f"{days}D"
    
    return iso_age if len(iso_age) > 1 else None


def parse_sex(sex_str):
    """
    Convert sex string to Phenopackets Sex enum.
    
    Returns: UNKNOWN_SEX, FEMALE, MALE, OTHER_SEX
    """
    if pd.isna(sex_str) or sex_str == 'Not reported':
        return "UNKNOWN_SEX"
    
    sex_str = str(sex_str).strip().lower()
    
    # Male variants
    if sex_str in ['male', 'm', 'man', 'boy']:
        return "MALE"
    
    # Female variants
    elif sex_str in ['female', 'f', 'woman', 'girl']:
        return "FEMALE"
    
    # Other
    elif sex_str in ['other', 'intersex']:
        return "OTHER_SEX"
    
    # Unknown
    else:
        return "UNKNOWN_SEX"


def parse_hpo_terms(hpo_str, hpo_labels):
    """
    Parse comma-separated HPO IDs into list of PhenotypicFeature objects.
    Includes labels from HPO ontology.
    
    Args:
        hpo_str: Comma-separated HPO IDs (e.g., "HP:0001250,HP:0001263")
        hpo_labels: Dictionary mapping HPO ID to label
    
    Returns:
        List of phenotypic feature dictionaries
    """
    if pd.isna(hpo_str) or hpo_str == 'Not reported' or hpo_str == '':
        return []
    
    # Split by comma and filter valid HPO IDs
    hpo_ids = [x.strip() for x in str(hpo_str).split(',') if x.strip().startswith('HP:')]
    
    features = []
    for hpo_id in hpo_ids:
        # Get label from ontology
        label = hpo_labels.get(hpo_id, f"Unknown term: {hpo_id}")
        
        features.append({
            "type": {
                "id": hpo_id,
                "label": label
            }
        })
    
    return features


def extract_gene_from_variant(variant_text):
    """
    Extract gene symbol from variant text.
    
    Handles various formats:
        - GENE:c.123A>G
        - GENE:NM_123456.1:c.123A>G
        - GENE p.Arg123Ter
        - GENE (c.123A>G)
    
    Args:
        variant_text: Variant description string
    
    Returns:
        Gene symbol string or None
    """
    if not variant_text or pd.isna(variant_text):
        return None
    
    variant_text = str(variant_text).strip()
    
    # Pattern 1: GENE:variant or GENE:transcript:variant
    match = re.match(r'^([A-Z][A-Z0-9]+):', variant_text)
    if match:
        gene = match.group(1)
        # Exclude common prefixes that aren't gene symbols
        if not any(gene.startswith(prefix) for prefix in ['NM_', 'NR_', 'NC_', 'CHR']):
            return gene
    
    # Pattern 2: GENE (variant) or GENE [variant]
    match = re.match(r'^([A-Z][A-Z0-9]+)\s*[\(\[]', variant_text)
    if match:
        return match.group(1)
    
    # Pattern 3: GENE followed by space and variant notation
    match = re.match(r'^([A-Z][A-Z0-9]+)\s+[cpgn]\.', variant_text)
    if match:
        return match.group(1)
    
    return None


def extract_omim_id(diagnosis_str):
    """
    Extract OMIM ID from diagnosis text.
    
    Handles formats:
        - Disease (OMIM:123456)
        - OMIM:123456
        - Disease [OMIM 123456]
        - Disease, OMIM 123456
    
    Args:
        diagnosis_str: Diagnosis/disease description
    
    Returns:
        OMIM ID string (e.g., "OMIM:123456") or None
    """
    if pd.isna(diagnosis_str) or diagnosis_str == 'Not reported':
        return None
    
    diagnosis_str = str(diagnosis_str)
    
    # Pattern 1: OMIM:123456
    match = re.search(r'OMIM:\s*(\d{6})', diagnosis_str, re.IGNORECASE)
    if match:
        return f"OMIM:{match.group(1)}"
    
    # Pattern 2: OMIM 123456 (without colon)
    match = re.search(r'OMIM\s+(\d{6})', diagnosis_str, re.IGNORECASE)
    if match:
        return f"OMIM:{match.group(1)}"
    
    # Pattern 3: #123456 (OMIM-style numbering)
    match = re.search(r'#(\d{6})', diagnosis_str)
    if match:
        return f"OMIM:{match.group(1)}"
    
    return None


def parse_acmg_classification(classification_str):
    """
    Map pathogenicity classification to Phenopackets ACMG enum.
    
    Args:
        classification_str: Classification string
    
    Returns:
        ACMG classification enum value
    """
    if pd.isna(classification_str) or classification_str == 'Not reported':
        return "NOT_PROVIDED"
    
    classification_map = {
        'pathogenic': 'PATHOGENIC',
        'likely pathogenic': 'LIKELY_PATHOGENIC',
        'uncertain significance': 'UNCERTAIN_SIGNIFICANCE',
        'vus': 'UNCERTAIN_SIGNIFICANCE',
        'likely benign': 'LIKELY_BENIGN',
        'benign': 'BENIGN'
    }
    
    classification_lower = str(classification_str).lower().strip()
    
    for key, value in classification_map.items():
        if key in classification_lower:
            return value
    
    return "NOT_PROVIDED"


def parse_variants(variant_str, subject_id, pathogenicity_str='Not reported'):
    """
    Parse variant string into correct Phenopackets v2 Interpretation structure.
    Enhanced with gene symbol extraction and ACMG classification.
    
    Args:
        variant_str: Semicolon-separated variant descriptions
        subject_id: Subject/patient identifier
        pathogenicity_str: ACMG pathogenicity classification
    
    Returns:
        List containing single interpretation dictionary, or empty list
    """
    if pd.isna(variant_str) or str(variant_str).strip() in ['Not reported', '']:
        return []
    
    # Split by semicolon and clean
    variant_parts = [v.strip() for v in str(variant_str).split(';') if v.strip()]
    
    if not variant_parts:
        return []
    
    # Parse ACMG classification once (applies to all variants)
    acmg_class = parse_acmg_classification(pathogenicity_str)
    
    # Build genomic interpretations
    genomic_interpretations = []
    
    for idx, variant_text in enumerate(variant_parts, 1):
        # Extract gene symbol
        gene_symbol = extract_gene_from_variant(variant_text)
        
        # Build variation descriptor
        variation_descriptor = {
            "id": f"{subject_id}_variant{idx}",
            "variation": {
                "text": variant_text
            }
        }
        
        # Add gene context if gene symbol found
        if gene_symbol:
            variation_descriptor["geneContext"] = {
                "valueId": f"HGNC:{gene_symbol}",
                "symbol": gene_symbol
            }
        
        # Build genomic interpretation
        genomic_interp = {
            "subjectOrBiosampleId": str(subject_id),
            "interpretationStatus": "UNKNOWN_STATUS",
            "variantInterpretation": {
                "variationDescriptor": variation_descriptor,
                "acmgPathogenicityClassification": acmg_class
            }
        }
        
        genomic_interpretations.append(genomic_interp)
    
    if not genomic_interpretations:
        return []
    
    # Create interpretation with required id field
    return [{
        "id": f"{subject_id}_interpretation",
        "progressStatus": "COMPLETED",
        "diagnosis": {
            "genomicInterpretations": genomic_interpretations
        }
    }]


def create_phenopacket(row, index, hpo_labels):
    """
    Create a Phenopacket v2 object from a data row.
    Enhanced version with all improvements.
    
    Args:
        row: Pandas Series containing patient data
        index: Row index (for tracking)
        hpo_labels: Dictionary mapping HPO ID to label
    
    Returns:
        Phenopacket dictionary
    """
    phenopacket_id = str(row['ID'])
    
    # Base phenopacket structure
    phenopacket = {
        "id": phenopacket_id,
        "subject": {
            "id": phenopacket_id,
        },
        "phenotypicFeatures": parse_hpo_terms(row['phenotypicFeatureIds'], hpo_labels),
        "metaData": {
            "created": datetime.now().isoformat() + "Z",
            "createdBy": "PAVS Data Pipeline",
            "resources": [
                {
                    "id": "hp",
                    "name": "Human Phenotype Ontology",
                    "url": "http://purl.obolibrary.org/obo/hp.owl",
                    "version": "2024-04-26",
                    "namespacePrefix": "HP",
                    "iriPrefix": "http://purl.obolibrary.org/obo/HP_"
                },
                {
                    "id": "geno",
                    "name": "Genotype Ontology",
                    "url": "http://purl.obolibrary.org/obo/geno.owl",
                    "version": "2020-03-08",
                    "namespacePrefix": "GENO",
                    "iriPrefix": "http://purl.obolibrary.org/obo/GENO_"
                },
                {
                    "id": "omim",
                    "name": "Online Mendelian Inheritance in Man",
                    "url": "https://www.omim.org",
                    "version": "2024-01",
                    "namespacePrefix": "OMIM",
                    "iriPrefix": "https://omim.org/entry/"
                }
            ],
            "phenopacketSchemaVersion": "2.0"
        }
    }
    
    # Add sex
    sex = parse_sex(row['sex'])
    if sex != "UNKNOWN_SEX":
        phenopacket["subject"]["sex"] = sex
    
    # Add age (enhanced parsing)
    age = parse_age(row['age'])
    if age:
        phenopacket["subject"]["timeAtLastEncounter"] = {
            "age": {
                "iso8601duration": age
            }
        }
    
    # Add genomic interpretations (enhanced with gene symbols and ACMG)
    variant_value = row.get('genomicVariants', 'Not reported')
    pathogenicity_value = row.get('variantInterpretation', 'Not reported')
    
    if pd.notna(variant_value) and str(variant_value).strip() not in ['Not reported', '', 'nan']:
        interpretations = parse_variants(variant_value, phenopacket_id, pathogenicity_value)
        if interpretations:
            phenopacket["interpretations"] = interpretations
    
    # Add diagnosis (enhanced with OMIM ID extraction)
    diagnosis_value = row.get('diagnosticComment', 'Not reported')
    if diagnosis_value != 'Not reported' and pd.notna(diagnosis_value):
        omim_id = extract_omim_id(diagnosis_value)
            
        phenopacket["diseases"] = [{
            "term": {
                "id": omim_id if omim_id else "",
                "label": str(diagnosis_value)
            }
        }]
            
    # Add external references
    phenopacket["metaData"]["externalReferences"] = []
    
    if row['externalReference'] != 'Not reported':
        refs = str(row['externalReference']).split(',')
        for ref in refs:
            ref = ref.strip()
            if ref:
                phenopacket["metaData"]["externalReferences"].append({
                    "id": ref,
                    "reference": ref,
                    "description": f"Source: {row['dataSourceType']}"
                })
    
    # Store PAVS-specific metadata
    pavs_metadata = {}
    
    if row.get('consanguinityStatus', 'Not reported') != 'Not reported':
        pavs_metadata['consanguinity'] = row['consanguinityStatus']
    
    if row.get('familyId', 'Not reported') != 'Not reported':
        pavs_metadata['familyId'] = row['familyId']
    
    if row.get('totalFamilyMembers', 'Not reported') != 'Not reported':
        pavs_metadata['totalFamilyMembers'] = row['totalFamilyMembers']
    
    if row.get('totalCohortMembers', 'Not reported') != 'Not reported':
        pavs_metadata['totalCohortMembers'] = row['totalCohortMembers']
    
    if row.get('procedure', 'Not reported') != 'Not reported':
        pavs_metadata['testingProcedure'] = row['procedure']
    
    if row.get('procedureStrategy', 'Not reported') != 'Not reported':
        pavs_metadata['testingStrategy'] = row['procedureStrategy']
    
    if row.get('diagnosticComment', 'Not reported') != 'Not reported':
        pavs_metadata['diagnosticComment'] = row['diagnosticComment']
    
    if row.get('zygosityStatus', 'Not reported') != 'Not reported':
        pavs_metadata['zygosityStatus'] = row['zygosityStatus']
    
    if row.get('dataSourceType', 'Not reported') != 'Not reported':
        pavs_metadata['dataSourceType'] = row['dataSourceType']
    
    # Add PAVS metadata to metaData section
    if pavs_metadata:
        phenopacket["metaData"]["pavsData"] = pavs_metadata
    
    return phenopacket


def convert_to_phenopackets(input_file, output_file, output_format='json', hpo_file='resources/hp.obo'):
    """
    Convert PAVS TSV data to Phenopackets v2 format.
    Enhanced version with all improvements.
    
    Args:
        input_file: Path to PAVS_final_data.tsv
        output_file: Path to output file
        output_format: 'json' for single JSON file, 'json-lines' for JSONL, 'individual' for separate files
        hpo_file: Path to HPO .obo file
    """
    print("=" * 60)
    print("Converting PAVS Data to Phenopackets v2 Format (Enhanced)")
    print("=" * 60)
    
    # Load HPO ontology
    hpo_labels = load_hpo_ontology(hpo_file)
    
    # Read data
    print(f"\nReading data from {input_file}...")
    data = pd.read_csv(input_file, sep='\t')
    print(f"✅ Loaded {len(data)} records")
    
    # Convert each row to a phenopacket
    print("\nConverting to Phenopackets...")
    phenopackets = []
    errors = []
    
    for idx, row in data.iterrows():
        try:
            phenopacket = create_phenopacket(row, idx, hpo_labels)
            phenopackets.append(phenopacket)
            
            if (idx + 1) % 100 == 0:
                print(f"   Processed {idx + 1}/{len(data)} records...")
                
        except Exception as e:
            error_msg = f"Error processing record {row['ID']}: {str(e)}"
            print(f"⚠️  {error_msg}")
            errors.append({'ID': row['ID'], 'error': str(e)})
            continue
    
    print(f"✅ Successfully converted {len(phenopackets)} records")
    
    if errors:
        print(f"⚠️  {len(errors)} records had errors")
        error_df = pd.DataFrame(errors)
        error_df.to_csv('phenopacket_conversion_errors.tsv', sep='\t', index=False)
        print("   Error details saved to: phenopacket_conversion_errors.tsv")
    
    # Save output
    print(f"\nSaving to {output_file}...")
    
    if output_format == 'json':
        # Single JSON file with array of phenopackets
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(phenopackets, f, indent=2, ensure_ascii=False)
        print(f"✅ Saved as single JSON file")
        
    elif output_format == 'json-lines':
        # JSON Lines format (one phenopacket per line)
        with open(output_file, 'w', encoding='utf-8') as f:
            for pp in phenopackets:
                f.write(json.dumps(pp, ensure_ascii=False) + '\n')
        print(f"✅ Saved as JSON Lines file")
        
    elif output_format == 'individual':
        # Individual JSON files
        import os
        output_dir = output_file.replace('.json', '_individual')
        os.makedirs(output_dir, exist_ok=True)
        
        for pp in phenopackets:
            individual_file = os.path.join(output_dir, f"{pp['id']}.json")
            with open(individual_file, 'w', encoding='utf-8') as f:
                json.dump(pp, f, indent=2, ensure_ascii=False)
        
        print(f"✅ Saved {len(phenopackets)} individual JSON files to {output_dir}/")
    
    # Print summary statistics
    print("\n" + "=" * 60)
    print("Conversion Summary")
    print("=" * 60)
    
    # Count phenopackets with various features
    with_phenotypes = sum(1 for pp in phenopackets if pp.get('phenotypicFeatures'))
    with_variants = sum(1 for pp in phenopackets if pp.get('interpretations'))
    with_diagnosis = sum(1 for pp in phenopackets if pp.get('diseases'))
    with_age = sum(1 for pp in phenopackets if 'timeAtLastEncounter' in pp.get('subject', {}))
    with_sex = sum(1 for pp in phenopackets if pp.get('subject', {}).get('sex') not in ['UNKNOWN_SEX', None])
    
    # Count features with labels
    total_features = sum(len(pp.get('phenotypicFeatures', [])) for pp in phenopackets)
    features_with_labels = sum(
        1 for pp in phenopackets 
        for feature in pp.get('phenotypicFeatures', [])
        if feature.get('type', {}).get('label') and 
           not feature.get('type', {}).get('label', '').startswith('Unknown term:')
    )
    
    # Count variants with gene symbols
    total_variants = sum(
        len(interp.get('diagnosis', {}).get('genomicInterpretations', []))
        for pp in phenopackets
        for interp in pp.get('interpretations', [])
    )
    
    variants_with_genes = sum(
        1 for pp in phenopackets
        for interp in pp.get('interpretations', [])
        for gi in interp.get('diagnosis', {}).get('genomicInterpretations', [])
        if 'geneContext' in gi.get('variantInterpretation', {}).get('variationDescriptor', {})
    )
    
    # Count diseases with OMIM IDs
    diseases_with_omim = sum(
        1 for pp in phenopackets
        for disease in pp.get('diseases', [])
        if disease.get('term', {}).get('id', '').startswith('OMIM:')
    )
    
    print(f"Total phenopackets: {len(phenopackets)}")
    print(f"\nSubject information:")
    print(f"  With sex specified: {with_sex} ({with_sex/len(phenopackets)*100:.1f}%)")
    print(f"  With age information: {with_age} ({with_age/len(phenopackets)*100:.1f}%)")
    
    print(f"\nPhenotypic features:")
    print(f"  Cases with phenotypes: {with_phenotypes} ({with_phenotypes/len(phenopackets)*100:.1f}%)")
    if total_features > 0:
        print(f"  Total HPO features: {total_features}")
        print(f"  Features with labels: {features_with_labels} ({features_with_labels/total_features*100:.1f}%)")
    
    print(f"\nGenomic variants:")
    print(f"  Cases with variants: {with_variants} ({with_variants/len(phenopackets)*100:.1f}%)")
    if total_variants > 0:
        print(f"  Total variants: {total_variants}")
        print(f"  Variants with gene symbols: {variants_with_genes} ({variants_with_genes/total_variants*100:.1f}%)")
    
    print(f"\nDiseases:")
    print(f"  Cases with diagnosis: {with_diagnosis} ({with_diagnosis/len(phenopackets)*100:.1f}%)")
    if with_diagnosis > 0:
        print(f"  Diagnoses with OMIM IDs: {diseases_with_omim} ({diseases_with_omim/with_diagnosis*100:.1f}%)")
    
    print("\n" + "=" * 60)
    print("Conversion Complete!")
    print("=" * 60)
    
    # Show sample outputs
    if phenopackets:
        print(f"\nSample phenopacket ID: {phenopackets[0]['id']}")
        
        if phenopackets[0].get('phenotypicFeatures'):
            sample_feature = phenopackets[0]['phenotypicFeatures'][0]
            print(f"\nSample HPO term:")
            print(f"  ID: {sample_feature['type']['id']}")
            print(f"  Label: {sample_feature['type']['label']}")
        
        if phenopackets[0].get('interpretations'):
            gi = phenopackets[0]['interpretations'][0].get('diagnosis', {}).get('genomicInterpretations', [])
            if gi:
                var_desc = gi[0].get('variantInterpretation', {}).get('variationDescriptor', {})
                print(f"\nSample variant:")
                print(f"  Text: {var_desc.get('variation', {}).get('text')}")
                if 'geneContext' in var_desc:
                    print(f"  Gene: {var_desc['geneContext'].get('symbol')}")
    
    return phenopackets


def main():
    parser = argparse.ArgumentParser(
        description="Convert PAVS data to Phenopackets v2 format (Enhanced version)"
    )
    parser.add_argument(
        '-i', '--input',
        default='PAVS_final_data.tsv',
        help='Input TSV file (default: PAVS_final_data.tsv)'
    )
    parser.add_argument(
        '-o', '--output',
        default='PAVS_phenopackets.json',
        help='Output file (default: PAVS_phenopackets.json)'
    )
    parser.add_argument(
        '-f', '--format',
        choices=['json', 'json-lines', 'individual'],
        default='json',
        help='Output format: json (single file), json-lines (JSONL), or individual (separate files)'
    )
    parser.add_argument(
        '--hpo',
        default='resources/hp.obo',
        help='Path to HPO .obo file (default: resources/hp.obo)'
    )
    
    args = parser.parse_args()
    
    convert_to_phenopackets(args.input, args.output, args.format, args.hpo)


if __name__ == "__main__":
    main()