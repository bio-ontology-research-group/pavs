#!/usr/bin/env python3
"""
Convert PAVS final data to Phenopackets v2 format (JSON).
Includes HPO term labels from ontology.
"""

import pandas as pd
import json
from datetime import datetime
import argparse
import pronto


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
    Parse age string into ISO8601 duration format.
    Examples: '22Y' -> 'P22Y', '6M' -> 'P6M', '3Y6M' -> 'P3Y6M'
    """
    if pd.isna(age_str) or age_str == 'Not reported' or age_str == '':
        return None
    
    age_str = str(age_str).strip().upper()
    
    # Already in ISO format
    if age_str.startswith('P'):
        return age_str
    
    # Parse common formats
    iso_age = 'P'
    
    # Years
    if 'Y' in age_str:
        years = age_str.split('Y')[0].strip()
        if years.isdigit():
            iso_age += f"{years}Y"
        age_str = age_str.split('Y')[1] if 'Y' in age_str else ''
    
    # Months
    if 'M' in age_str:
        months = age_str.split('M')[0].strip()
        if months.isdigit():
            iso_age += f"{months}M"
    
    # Days
    if 'D' in age_str:
        days = age_str.split('D')[0].strip()
        if days.isdigit():
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
    
    if sex_str in ['male', 'm']:
        return "MALE"
    elif sex_str in ['female', 'f']:
        return "FEMALE"
    elif sex_str in ['other']:
        return "OTHER_SEX"
    else:
        return "UNKNOWN_SEX"


def parse_hpo_terms(hpo_str, hpo_labels):
    """
    Parse comma-separated HPO IDs into list of PhenotypicFeature objects.
    Now includes labels from HPO ontology.
    """
    if pd.isna(hpo_str) or hpo_str == 'Not reported' or hpo_str == '':
        return []
    
    hpo_ids = [x.strip() for x in str(hpo_str).split(',') if x.strip().startswith('HP:')]
    
    features = []
    for hpo_id in hpo_ids:
        # Get label from ontology, or use descriptive text if not found
        label = hpo_labels.get(hpo_id, f"Unknown term: {hpo_id}")
        
        features.append({
            "type": {
                "id": hpo_id,
                "label": label
            }
        })
    
    return features

def parse_variants(variant_str, subject_id):
    """
    Parse variant string into correct Phenopackets v2 Interpretation structure.
    """
    if pd.isna(variant_str) or str(variant_str).strip() in ['Not reported', '']:
        return []
    
    # Split and clean variants
    variant_parts = [v.strip() for v in str(variant_str).split(';') if v.strip()]
    
    if not variant_parts:
        return []
    
    # Build genomic interpretations
    genomic_interpretations = []
    for idx, variant_text in enumerate(variant_parts, 1):
        genomic_interpretations.append({
            "subjectOrBiosampleId": str(subject_id),
            "interpretationStatus": "UNKNOWN_STATUS",
            "variantInterpretation": {
                "variationDescriptor": {
                    "id": f"{subject_id}_variant{idx}",
                    "variation": {
                        "text": variant_text
                    }
                },
                "acmgPathogenicityClassification": "NOT_PROVIDED"
            }
        })
    
    if not genomic_interpretations:
        return []
    
    # Create interpretation with REQUIRED id field
    return [{
        "id": f"{subject_id}_interpretation1",
        "progressStatus": "COMPLETED",
        "diagnosis": {
            "genomicInterpretations": genomic_interpretations
        }
    }]

def parse_variants1(variant_str):
    """
    Parse variant string into list of genomic interpretations.
    This is a simplified parser - you may need to enhance based on your variant format.
    """
    if pd.isna(variant_str) or variant_str == 'Not reported' or variant_str == '':
        return []
    
    # Split by semicolon (based on your data format)
    variant_parts = [x.strip() for x in str(variant_str).split(';') if x.strip()]
    
    interpretations = []
    for variant in variant_parts:
        if not variant:
            continue
            
        # Basic structure -
        interpretation = {
            "subjectOrBiosampleId": "",
            "interpretationStatus": "UNKNOWN_STATUS",
            "variantInterpretation": {
                "variationDescriptor": {
                    "id": "",
                    "variation": {
                        "text": variant  
                    }
                }
            }
        }
        interpretations.append(interpretation)
    
    return interpretations


def parse_zygosity(zygosity_str):
    """
    Parse zygosity string.
    Common values: homozygous, heterozygous, hemizygous, compound heterozygous
    """
    if pd.isna(zygosity_str) or zygosity_str == 'Not reported':
        return None
    
    zygosity_map = {
        'homozygous': 'GENO:0000136',
        'heterozygous': 'GENO:0000135',
        'hemizygous': 'GENO:0000134',
        'compound heterozygous': 'GENO:0000402'
    }
    
    zygosity_lower = str(zygosity_str).lower().strip()
    
    for key, value in zygosity_map.items():
        if key in zygosity_lower:
            return {
                "id": value,
                "label": key
            }
    
    return None


def create_phenopacket(row, index, hpo_labels):
    """
    Create a Phenopacket v2 object from a data row.
    Now accepts hpo_labels dictionary for term labels.
    """
    #phenopacket_id = row['ID']
    phenopacket_id = str(row['ID'])
    
    # Base phenopacket structure
    phenopacket = {
        "id": phenopacket_id,
        "subject": {
            "id": phenopacket_id,
        },
        "phenotypicFeatures": parse_hpo_terms(row['phenotypicFeatureIds'], hpo_labels),
        "interpretations": [],
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
                }
            ],
            "phenopacketSchemaVersion": "2.0"
        }
    }
    
    # Add sex
    sex = parse_sex(row['sex'])
    if sex != "UNKNOWN_SEX":
        phenopacket["subject"]["sex"] = sex
    
    # Add age
    age = parse_age(row['age'])
    if age:
        phenopacket["subject"]["timeAtLastEncounter"] = {
            "age": {
                "iso8601duration": age
            }
        }

    # Add genomic interpretations (variants)
    variant_value = row.get('genomicVariants', 'Not reported')
    if pd.notna(variant_value) and str(variant_value).strip() not in ['Not reported', '', 'nan']:
        interpretations = parse_variants(variant_value, phenopacket_id)
        if interpretations:
            phenopacket["interpretations"] = interpretations
    '''
    # Add genomic interpretations (variants)
    if row['genomicVariants'] != 'Not reported':
        interpretations = parse_variants(row['genomicVariants'])
        if interpretations:
            phenopacket["interpretations"] = interpretations
    '''

    # Add diagnosis
    if row['diagnosis'] != 'Not reported' and pd.notna(row['diagnosis']):
        phenopacket["diseases"] = [{
            "term": {
                "id": "",  # Add OMIM/disease ontology ID if available
                "label": str(row['diagnosis'])
            }
        }]
    
    # Add additional metadata as extensions
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
    
    # Add custom fields that don't fit standard schema
    if 'extensions' not in phenopacket:
        phenopacket['extensions'] = []
    
    # Store additional PAVS-specific data
    pavs_extension = {
        "type": "PAVS_clinical_data",
        "value": {}
    }
    
    if row['consanguinityStatus'] != 'Not reported':
        pavs_extension["value"]["consanguinity"] = row['consanguinityStatus']
    
    if row['familyId'] != 'Not reported':
        pavs_extension["value"]["familyId"] = row['familyId']
    
    if row['totalFamilyMembers'] != 'Not reported':
        pavs_extension["value"]["totalFamilyMembers"] = row['totalFamilyMembers']
    
    if row['totalCohortMembers'] != 'Not reported':
        pavs_extension["value"]["totalCohortMembers"] = row['totalCohortMembers']
    
    if row['procedure'] != 'Not reported':
        pavs_extension["value"]["testingProcedure"] = row['procedure']
    
    if row['procedureStrategy'] != 'Not reported':
        pavs_extension["value"]["testingStrategy"] = row['procedureStrategy']
    
    if row['diagnosticComment'] != 'Not reported':
        pavs_extension["value"]["diagnosticComment"] = row['diagnosticComment']
    
    if row['zygosityStatus'] != 'Not reported':
        pavs_extension["value"]["zygosityStatus"] = row['zygosityStatus']
    
    if row['variantInterpretation'] != 'Not reported':
        pavs_extension["value"]["variantInterpretation"] = row['variantInterpretation']
    
    if row['dataSourceType'] != 'Not reported':
        pavs_extension["value"]["dataSourceType"] = row['dataSourceType']
    
    if pavs_extension["value"]:
        phenopacket['extensions'].append(pavs_extension)
    
    return phenopacket


def convert_to_phenopackets(input_file, output_file, output_format='json', hpo_file='resources/hp.obo'):
    """
    Convert PAVS TSV data to Phenopackets format.
    
    Args:
        input_file: Path to PAVS_final_data.tsv
        output_file: Path to output file
        output_format: 'json' for single JSON file, 'json-lines' for JSONL, 'individual' for separate files
        hpo_file: Path to HPO .obo file
    """
    print("=" * 60)
    print("Converting PAVS Data to Phenopackets v2 Format")
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
    
    for idx, row in data.iterrows():
        try:
            phenopacket = create_phenopacket(row, idx, hpo_labels)
            phenopackets.append(phenopacket)
            
            if (idx + 1) % 100 == 0:
                print(f"   Processed {idx + 1}/{len(data)} records...")
                
        except Exception as e:
            print(f"⚠️  Error processing record {row['ID']}: {e}")
            continue
    
    print(f"✅ Successfully converted {len(phenopackets)} records")
    
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
    
    # Count features with labels
    total_features = sum(len(pp.get('phenotypicFeatures', [])) for pp in phenopackets)
    features_with_labels = sum(
        1 for pp in phenopackets 
        for feature in pp.get('phenotypicFeatures', [])
        if feature.get('type', {}).get('label') and 
           not feature.get('type', {}).get('label', '').startswith('Unknown term:')
    )
    
    print(f"Total phenopackets: {len(phenopackets)}")
    print(f"With phenotypic features: {with_phenotypes} ({with_phenotypes/len(phenopackets)*100:.1f}%)")
    if total_features > 0:
        print(f"Total HPO features: {total_features}")
        print(f"HPO features with labels: {features_with_labels} ({features_with_labels/total_features*100:.1f}%)")
    print(f"With genomic variants: {with_variants} ({with_variants/len(phenopackets)*100:.1f}%)")
    print(f"With diagnosis: {with_diagnosis} ({with_diagnosis/len(phenopackets)*100:.1f}%)")
    print(f"With age information: {with_age} ({with_age/len(phenopackets)*100:.1f}%)")
    
    print("\n" + "=" * 60)
    print("Conversion Complete!")
    print("=" * 60)
    
    # Show sample phenotypic feature with label
    if phenopackets and phenopackets[0].get('phenotypicFeatures'):
        sample_feature = phenopackets[0]['phenotypicFeatures'][0]
        print(f"\nSample HPO term with label:")
        print(f"  ID: {sample_feature['type']['id']}")
        print(f"  Label: {sample_feature['type']['label']}")
    
    return phenopackets


def main():
    parser = argparse.ArgumentParser(
        description="Convert PAVS data to Phenopackets v2 format"
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