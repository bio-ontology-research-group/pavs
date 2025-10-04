# Import necessary libraries
import os
import json
import pandas as pd
from tqdm.notebook import tqdm
from pyhpo import Ontology, HPOSet

# --- Download necessary HPO files ---
# pyhpo automatically downloads the required files on first instantiation
# if they are not found in its data directory.
print("Initializing HPO Ontology... (This may download files on first run)")
_ = Ontology()
print("Ontology initialized successfully.")

# --- Configuration ---
PHENOPACKET_DIR = 'phenopackets'

# The Ontology object, once initialized, contains all disease annotations.
# We can access them directly.
omim_diseases = Ontology.omim_diseases

print(f"Loaded {len(omim_diseases)} OMIM disease profiles from the HPO annotations.")

# For faster lookup, create a dictionary mapping OMIM ID to the HPOSet
omim_profiles = {disease.id: disease.hpo_set for disease in omim_diseases}

def calculate_similarity_ranking(patient_hpo_set, omim_disease_profiles):
    """
    Calculates the similarity between a patient's HPOSet and all OMIM diseases,
    returning a ranked list of diseases.
    
    Args:
        patient_hpo_set (HPOSet): The set of HPO terms for the patient.
        omim_disease_profiles (dict): A dictionary of OMIM IDs to HPOSets.
        
    Returns:
        pandas.DataFrame: A DataFrame with OMIM IDs and their similarity scores,
                          sorted in descending order.
    """
    if not patient_hpo_set:
        return pd.DataFrame(columns=['omim_id', 'similarity_score'])
        
    results =
    for omim_id, disease_hpo_set in omim_disease_profiles.items():
        # The core similarity calculation using pyhpo's default method
        score = patient_hpo_set.similarity(disease_hpo_set)
        results.append({'omim_id': omim_id, 'similarity_score': score})
        
    # Create a DataFrame and sort by score
    ranked_df = pd.DataFrame(results)
    ranked_df = ranked_df.sort_values(by='similarity_score', ascending=False).reset_index(drop=True)
    
    return ranked_df


def run_validation(phenopacket_dir, omim_profiles):
    """
    Iterates through all PhenoPackets, performs similarity analysis,
    and evaluates the ranking of the true diagnosis.
    """
    phenopacket_files = [f for f in os.listdir(phenopacket_dir) if f.endswith('.json')]
    
    results =
    
    for filename in tqdm(phenopacket_files, desc="Validating PhenoPackets"):
        filepath = os.path.join(phenopacket_dir, filename)
        with open(filepath, 'r') as f:
            data = json.load(f)
            
        # Extract patient HPO terms
        patient_hpo_ids = [pf['type']['id'] for pf in data.get('phenotypicFeatures',)]
        if not patient_hpo_ids:
            continue # Skip if no phenotypes are recorded
            
        patient_hpo_set = HPOSet.from_queries(patient_hpo_ids)
        
        # Extract ground truth diagnosis (taking the first one if multiple)
        ground_truth_diseases = [d['term']['id'] for d in data.get('diseases',)]
        if not ground_truth_diseases:
            continue # Skip if no diagnosis is recorded
            
        ground_truth_omim_id = ground_truth_diseases
        
        # Get the ranked list of candidate diseases
        ranked_diseases = calculate_similarity_ranking(patient_hpo_set, omim_profiles)
        
        # Find the rank of the true diagnosis
        rank_info = ranked_diseases[ranked_diseases['omim_id'] == ground_truth_omim_id]
        
        if not rank_info.empty:
            rank = rank_info.index + 1 # index is 0-based, rank is 1-based
            score = rank_info['similarity_score'].iloc
        else:
            rank = float('inf') # Not found
            score = 0.0
            
        results.append({
            'phenopacket_id': data['id'],
            'ground_truth_omim': ground_truth_omim_id,
            'rank': rank,
            'score': score,
            'num_hpo_terms': len(patient_hpo_ids)
        })
        
    return pd.DataFrame(results)

# --- Execute the validation ---
validation_results_df = run_validation(PHENOPACKET_DIR, omim_profiles)

print("\nValidation Results Summary:")
print(validation_results_df.head())


# --- Calculate Top-N Accuracies ---
total_cases = len(validation_results_df)
if total_cases > 0:
    top1_accuracy = (validation_results_df['rank'] == 1).sum() / total_cases * 100
    top5_accuracy = (validation_results_df['rank'] <= 5).sum() / total_cases * 100
    top10_accuracy = (validation_results_df['rank'] <= 10).sum() / total_cases * 100
    median_rank = validation_results_df['rank'].median()

    print("\n--- Overall Performance Metrics ---")
    print(f"Total cases evaluated: {total_cases}")
    print(f"Top-1 Accuracy: {top1_accuracy:.2f}%")
    print(f"Top-5 Accuracy: {top5_accuracy:.2f}%")
    print(f"Top-10 Accuracy: {top10_accuracy:.2f}%")
    print(f"Median Rank of Correct Diagnosis: {median_rank}")
else:
    print("No valid cases with both phenotypes and diagnoses were found to evaluate.")

# Save detailed results for further analysis or inclusion in supplementary materials
validation_results_df.to_csv('phenotypic_validation_results.csv', index=False)
