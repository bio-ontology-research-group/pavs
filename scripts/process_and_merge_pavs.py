import pandas as pd
import subprocess
import json
import re
import warnings

# Hide that specific warning
warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")

def process_literature_data(dataset):
    """Process literature dataset (PAVS & DDIEM)"""
    print("Processing literature data...")
    
    # Load the dataset
    data = pd.read_csv(dataset, encoding='latin-1', sep='\t')
    data = data.dropna(how='all')
    
    # Rename columns
    data.rename(columns={
        'id': 'ID',
        'patient_gender': 'gender',
        'patient_age': 'age'
    }, inplace=True)
    
    # Filter out specific reference
    data = data[data['reference'] != 'https://pubmed.ncbi.nlm.nih.gov/28600779/']
    
    # Set IDs
    data['old_id'] = data['ID']
    data['ID'] = range(1, len(data) + 1)
    data = data.drop_duplicates()
    
    # Add data source type
    data['dataSourceType'] = 'literature'
    # Extract HPO IDs (if any)
    data['HPO_ID'] = data['phenotypes'].apply(lambda x: ','.join(list(dict.fromkeys(re.findall(r'HP:\d+', x)))))

    
    # Select final columns
    data = data[['ID', 'test', 'test_strategy', 'gender', 'age', 'consanguinity', 
                 'family_id', 'number_of_family_members', 'cohort_size', 'phenotypes', 
                 'result', 'result_comment', 'variants', 'zygosity', 'pathogenicity', 
                 'reference', 'dataSourceType', 'HPO_ID']]
    
    print(f"Literature data processed: {data.shape}")
    return data


def process_clinical_literature_data(dataset):
    """Process clinical literature dataset (Dr Fawzan)"""
    print("Processing clinical literature data...")
    
    # Load the dataset
    data = pd.read_csv(dataset, encoding='latin-1', sep='\t')
    data = data.dropna(how='all')
    
    
    # Set IDs
    data['old_id'] = data['ID']
    data['ID'] = range(1, len(data) + 1)
    data = data.drop_duplicates()
    
    # Rename columns
    data.rename(columns={
        'id': 'ID',
        'Testing Strategy': 'test_strategy',
        'Gender': 'gender',
        'Age': 'age',
        'Test': 'test',
        'Consanguinity': 'consanguinity',
        'Family hx': 'result_comment',
        'Result': 'result',
        'Zyogsity': 'zygosity',
        'Phenotype': 'phenotypes',    
    }, inplace=True)
    
    # Process variants and other fields
    data['variants'] = data['Variant(s)']
    data['pathogenicity'] = data['HGMD']
    data['reference'] = 'https://pubmed.ncbi.nlm.nih.gov/28600779/'
    data['family_id'] = data['ID']
    data['dataSourceType'] = 'clinical_literature'
    data['number_of_family_members'] = 'Not reported'
    data['cohort_size'] = 'Not reported'

    # Extract HPO IDs (if any)
    data['HPO_ID'] = data['phenotypes'].apply(lambda x: ','.join(list(dict.fromkeys(re.findall(r'HP:\d+', x)))))
    
    # Select final columns
    data = data[['ID', 'test', 'test_strategy', 'gender', 'age', 'consanguinity', 
                 'family_id', 'number_of_family_members', 'cohort_size', 'phenotypes', 
                 'result', 'result_comment', 'variants', 'zygosity', 'pathogenicity', 
                 'reference', 'dataSourceType','HPO_ID']]
    
    print(f"Clinical literature data processed: {data.shape}")
    return data


def process_hospital_collaborator_data(dataset):
    """Process hospital collaborator dataset (National Guards)"""
    print("Processing hospital collaborator data...")
    
    # Load the dataset
    data = pd.read_excel(dataset, engine="openpyxl")    

    # Rename columns
    data = data.rename(columns={
        'Case': 'ID',
        'Test': 'test',
        ' Test Type': 'test_strategy',
        'Gender': 'gender',
        'DOB': 'age',
        'Consanguinity ': 'consanguinity',
        'Results Internal': 'result',
        'Zygosity': 'zygosity',
        'HPOs': 'phenotypes',    
        'Family ID': 'family_id'})
        
    # Combine variant columns
    variant_columns = ['Gene', 'Variant ', 'Gene 2', 'Variant 2  ', 'Gene 3', 'Variant 3  ', 'Gene 4', 'Variant4']
    data['variants'] = data[variant_columns].fillna('').agg(';'.join, axis=1)

    # Combine pathogenicity columns
    pathogenicity_columns = ['pathogenicity ', 'pathogenicity 2 ', 'pathogenicity 3', 'pathogenicity 4']
    data['pathogenicity'] = data[pathogenicity_columns].fillna('').agg(';'.join, axis=1)
    
    # Combine comment fields
    data['result_comment'] = data.apply(
        lambda row: ';'.join([
            f"Comments: {row['Comments']}" if pd.notna(row['Comments']) else '',
            f"Inheritance: {row['Inheritance']}" if pd.notna(row['Inheritance']) else '',
            f"OMIM: {row['Omim']}" if pd.notna(row['Omim']) else ''
        ]).strip(';'), axis=1
    )
    
    # Set metadata fields
    data['reference'] = 'hospital_collaborator, https://link.springer.com/article/10.1186/s12920-020-00743-8 , https://onlinelibrary.wiley.com/doi/full/10.1111/cge.13842'
    data['dataSourceType'] = 'clinical_collaborator'
    data['number_of_family_members'] = data['test_strategy']
    data['cohort_size'] = data['test_strategy']
    
    # Extract HPO IDs
    data['HPO_ID'] = data['phenotypes'].apply(lambda x: ','.join(list(dict.fromkeys(re.findall(r'HP:\d+', x)))))

    
    # Select final columns
    data = data[['ID', 'test', 'test_strategy', 'gender', 'age', 'consanguinity', 
                 'family_id', 'number_of_family_members', 'cohort_size', 'phenotypes', 
                 'result', 'result_comment', 'variants', 'zygosity', 'pathogenicity', 
                 'reference', 'dataSourceType', 'HPO_ID']]
    
    print(f"Hospital collaborator data processed: {data.shape}")
    
    return data

def create_combined_phenotype_file(data1, data2, output_file='data/id_pheno_combined.tsv'):
    """
    Combine phenotype text from data1 and data2 into one TSV file for HPO annotation.
    Add prefixes to IDs to keep them unique.
    """
    print("Creating combined phenotype file for annotation...")
    
    # Extract ID and phenotypes from data1 with prefix
    pheno1 = data1[['ID', 'phenotypes']].copy()
    pheno1['ID'] = 'LIT_' + pheno1['ID'].astype(str)
    pheno1['source'] = 'literature'
    
    # Extract ID and phenotypes from data2 with prefix
    pheno2 = data2[['ID', 'phenotypes']].copy()
    pheno2['ID'] = 'CLIT_' + pheno2['ID'].astype(str)
    pheno2['source'] = 'clinical_literature'
    
    # Combine both
    combined = pd.concat([pheno1, pheno2], ignore_index=True)
    
    # Save as TSV (only ID and phenotypes columns for the annotation script)
    combined[['ID', 'phenotypes']].to_csv(output_file, sep='\t', index=False)
    
    print(f"✅ Combined phenotype file saved: {output_file}")
    print(f"   Total records: {len(combined)}")
    print(f"   Literature records: {len(pheno1)} (prefixed with LIT_)")
    print(f"   Clinical literature records: {len(pheno2)} (prefixed with CLIT_)")
    
    return output_file


def merge_hpo_ids(data, annotation_file, id_prefix):
    """
    Merge HPO annotations with existing HPO IDs from regex extraction.
    Keep only unique HPO IDs.
    
    Args:
        data: DataFrame to merge annotations into
        annotation_file: JSON file with annotations
        id_prefix: Prefix used for IDs in annotation file (e.g., 'LIT_' or 'CLIT_')
    """
    print(f"\nMerging HPO annotations for {id_prefix} records...")
    
    # Load annotations
    with open(annotation_file, 'r') as f:
        annotations = json.load(f)
    
    # Convert annotations to DataFrame, filtering by prefix
    hpo_list = []
    for pid, hpo_terms in annotations.items():
        if pid.startswith(id_prefix):
            # Remove prefix to match original IDs
            original_id = pid.replace(id_prefix, '', 1)
            hpo_list.append({
                'ID': original_id,
                'HPO_ID_annotated': ','.join(hpo_terms)
            })
    
    hpo_df = pd.DataFrame(hpo_list)
    
    # Merge with main data
    data['ID'] = data['ID'].astype(str)
    data = data.merge(hpo_df, on='ID', how='left')
    
    # Combine and deduplicate HPO IDs
    def combine_unique_hpo(row):
        ids = set()
        
        # Add regex-extracted HPO IDs (if exist)
        if pd.notna(row.get('HPO_ID')) and row['HPO_ID']:
            ids.update([x.strip() for x in str(row['HPO_ID']).split(',') if x.strip()])
        
        # Add annotated HPO IDs
        if pd.notna(row.get('HPO_ID_annotated')) and row['HPO_ID_annotated']:
            ids.update([x.strip() for x in row['HPO_ID_annotated'].split(',') if x.strip()])
        
        # Remove empty strings and invalid IDs
        ids = {x for x in ids if x and x.startswith('HP:')}
        
        return ','.join(sorted(ids)) if ids else ''
    
    data['HPO_ID'] = data.apply(combine_unique_hpo, axis=1)
    
    # Drop temporary column
    if 'HPO_ID_annotated' in data.columns:
        data = data.drop(columns=['HPO_ID_annotated'])
    
    print(f"✅ HPO IDs merged and deduplicated for {id_prefix} records")
    
    return data

def run_hpo_annotation(input_tsv, output_json, hpo_script="map_text_to_hpo.py"):
    """
    Run the HPO annotation script.
    """
    print(f"\nRunning HPO annotation on {input_tsv}...")
    
    cmd = [
        "python3", hpo_script,
        "-i", input_tsv,
        "-o", output_json,
        "--hpo", "resources/hp.obo"
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"✅ HPO annotations saved to {output_json}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running annotation script: {e}")
        return False


def annotate_and_merge_hpo(data1, data2, data3):
    """
    Complete workflow: create combined file, annotate, and merge back.
    """
    print("\n" + "=" * 60)
    print("HPO Annotation Workflow")
    print("=" * 60)
    
    # Step 1: Extract HPO IDs from existing text (regex)
    print("\nStep 1: Extracting HPO IDs from text using regex...")
    for data, name in [(data1, 'data1'), (data2, 'data2'), (data3, 'data3')]:
        data['HPO_ID'] = data['phenotypes'].apply(
            lambda x: ','.join(list(dict.fromkeys(re.findall(r'HP:\d+', x)))) if pd.notna(x) else ''
        )
        existing_count = (data['HPO_ID'] != '').sum()
        print(f"   {name}: {existing_count} records with existing HPO IDs")
    
    # Step 2: Create combined phenotype file for data1 and data2 with unique IDs
    combined_file = create_combined_phenotype_file(data1, data2)
    
    # Step 3: Run HPO annotation
    annotation_output = 'data/hpo_annotations_combined.json'
    success = run_hpo_annotation(combined_file, annotation_output)
    
    if not success:
        print("⚠️  HPO annotation failed. Continuing with regex-extracted IDs only.")
        return data1, data2, data3
    
    # Step 4: Merge annotations back into data1 and data2 using their respective prefixes
    print("\nStep 4: Merging annotations back into datasets...")
    data1 = merge_hpo_ids(data1, annotation_output, id_prefix='LIT_')
    data2 = merge_hpo_ids(data2, annotation_output, id_prefix='CLIT_')
    
    # Step 5: Report statistics
    print("\n" + "=" * 60)
    print("HPO Annotation Statistics")
    print("=" * 60)
    
    for data, name in [(data1, 'Literature'), (data2, 'Clinical Literature'), (data3, 'Hospital Collaborator')]:
        total = len(data)
        with_hpo = (data['HPO_ID'] != '').sum()
        without_hpo = total - with_hpo
        print(f"{name}:")
        print(f"   Total records: {total}")
        print(f"   With HPO IDs: {with_hpo} ({with_hpo/total*100:.1f}%)")
        print(f"   Without HPO IDs: {without_hpo} ({without_hpo/total*100:.1f}%)")
        
        # Sample of combined HPO IDs
        if with_hpo > 0:
            sample = data[data['HPO_ID'] != '']['HPO_ID'].iloc[0]
            num_hpo = len(sample.split(',')) if sample else 0
            print(f"   Sample (first record): {num_hpo} HPO terms")
    
    return data1, data2, data3


def main():
    """Main execution function"""
    print("=" * 60)
    print("Starting data processing pipeline")
    print("=" * 60)
    
    # Process each dataset
    data1 = process_literature_data("../data/PAVS_DDIEM_lit.tsv")
    data2 = process_clinical_literature_data("../data/439_2017_1821_MOESM1_ESM.txt")
    data3 = process_hospital_collaborator_data("../data/Variant_list_National_guards_transcripts.xlsx")
    
    # Annotate and merge HPO terms (combines regex + NLP annotation)
    data1, data2, data3 = annotate_and_merge_hpo(data1, data2, data3)
    
    print("\n" + "=" * 60)
    print("Combining all datasets")
    print("=" * 60)
    
    # Print dataset sizes
    print(f'Literature size: {data1.shape}')
    print(f'Clinical literature size: {data2.shape}')
    print(f'Clinical collaborator size: {data3.shape}')
    
    # Combine all datasets
    final_data = pd.concat([data1, data2, data3], ignore_index=True)
    print(f'Total combined: {final_data.shape}')
    
    # Generate new IDs
    final_data['ID'] = ['PAVS' + str(i+1) for i in range(len(final_data))]
    final_data = final_data.drop_duplicates()
    
    # Rename columns to align with Phenopackets schema
    new_column_names = {
        'ID': 'ID',
        'test': 'procedure',
        'test_strategy': 'procedureStrategy',
        'gender': 'sex',
        'age': 'age',  
        'consanguinity': 'consanguinityStatus',
        'family_id': 'familyId',
        'number_of_family_members': 'totalFamilyMembers',
        'cohort_size': 'totalCohortMembers',
        'phenotypes': 'phenotypicFeatures',
        'result': 'diagnosis',
        'result_comment': 'diagnosticComment',
        'variants': 'genomicVariants',
        'zygosity': 'zygosityStatus',
        'pathogenicity': 'variantInterpretation',
        'HPO_ID': 'phenotypicFeatureIds',
        'dataSourceType': 'dataSourceType',
        'reference': 'externalReference',
    }
    
    final_data.rename(columns=new_column_names, inplace=True)
    
    # Organize columns in logical order
    final_data = final_data[[
        'ID',
        'sex',
        'age',
        'consanguinityStatus',
        'familyId',
        'totalFamilyMembers',
        'totalCohortMembers',
        'phenotypicFeatures',
        'phenotypicFeatureIds',
        'procedure',
        'procedureStrategy',
        'diagnosis',
        'diagnosticComment',
        'genomicVariants',
        'zygosityStatus',
        'variantInterpretation',
        'dataSourceType',
        'externalReference'
    ]]
    
    # Replace NaN values
    final_data = final_data.fillna('Not reported')
    
    print("\n" + "=" * 60)
    print("Saving final dataset")
    print("=" * 60)
    
    # Filter out records without HPO terms
    records_without_hpo = (final_data['phenotypicFeatureIds'] == 'Not reported').sum()
    print(f"Records without HPO IDs: {records_without_hpo}")
    
    #final_data = final_data[final_data['phenotypicFeatureIds'] != 'Not reported']
    
    # Save final data
    final_data.to_csv('PAVS_final_data.tsv', sep='\t', index=False)
    
    print(f'Final data saved: {final_data.shape}')
    print("\n" + "=" * 60)
    print("Processing complete!")
    print("=" * 60)
    print(f"\nFirst few rows:\n{final_data.head()}")


if __name__ == "__main__":
    main()