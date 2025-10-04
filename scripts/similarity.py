# Import necessary libraries
import os
import json
import argparse
import logging
import pandas as pd
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from pyhpo import Ontology, HPOSet


def initialize_hpo_ontology():
    """
    Initializes the HPO Ontology and loads OMIM disease profiles.
    pyhpo automatically downloads required files if they are not found.
    """
    logging.info("Initializing HPO Ontology... (This may download files on first run)")
    _ = Ontology()
    logging.info("Ontology initialized successfully.")
    
    omim_diseases = Ontology.omim_diseases
    logging.info(f"Loaded {len(omim_diseases)} OMIM disease profiles from HPO annotations.")
    
    # For faster lookup, create a dictionary mapping OMIM ID to the HPOSet
    omim_profiles = {disease.id: disease.hpo_set for disease in omim_diseases}
    return omim_profiles

def calculate_similarity_ranking(patient_hpo_set, omim_disease_profiles, similarity_method='resnik'):
    """
    Calculates the similarity between a patient's HPOSet and all OMIM diseases,
    returning a ranked list of diseases.
    
    Args:
        patient_hpo_set (HPOSet): The set of HPO terms for the patient.
        omim_disease_profiles (dict): A dictionary of OMIM IDs to HPOSets.
        similarity_method (str): The similarity algorithm to use ('resnik', 'lin', 'jc', 'rel').
        
    Returns:
        pandas.DataFrame: A DataFrame with OMIM IDs and their similarity scores,
                          sorted in descending order.
    """
    if not patient_hpo_set:
        return pd.DataFrame(columns=['omim_id', 'similarity_score'])
        
    results = []
    for omim_id, disease_hpo_set in omim_disease_profiles.items():
        score = patient_hpo_set.similarity(disease_hpo_set, kind=similarity_method)
        results.append({'omim_id': omim_id, 'similarity_score': score})
        
    ranked_df = pd.DataFrame(results)
    ranked_df = ranked_df.sort_values(by='similarity_score', ascending=False).reset_index(drop=True)
    
    return ranked_df

def run_validation(phenopacket_dir, omim_profiles, similarity_method):
    """
    Iterates through all PhenoPackets, performs similarity analysis,
    and evaluates the ranking of the true diagnosis.
    
    Args:
        phenopacket_dir (str): Path to the directory containing PhenoPacket JSON files.
        omim_profiles (dict): A dictionary of OMIM IDs to HPOSets.
        similarity_method (str): The similarity algorithm to use.
        
    Returns:
        pandas.DataFrame: A DataFrame with detailed validation results for each phenopacket.
    """
    phenopacket_files = [f for f in os.listdir(phenopacket_dir) if f.endswith('.json')]
    
    results = []
    
    with logging_redirect_tqdm():
        for filename in tqdm(phenopacket_files, desc="Validating PhenoPackets"):
            filepath = os.path.join(phenopacket_dir, filename)
            with open(filepath, 'r') as f:
                data = json.load(f)
                
            patient_hpo_ids = [pf['type']['id'] for pf in data.get('phenotypicFeatures', [])]
            if not patient_hpo_ids:
                logging.warning(f"Skipping {data['id']}: No phenotypic features found.")
                continue
                
            patient_hpo_set = HPOSet.from_queries(patient_hpo_ids)
            
            ground_truth_diseases = [d['term']['id'] for d in data.get('diseases', [])]
            if not ground_truth_diseases:
                logging.warning(f"Skipping {data['id']}: No diagnosis found.")
                continue
                
            ground_truth_omim_id = ground_truth_diseases[0]
            
            ranked_diseases = calculate_similarity_ranking(patient_hpo_set, omim_profiles, similarity_method=similarity_method)
            
            rank_info = ranked_diseases[ranked_diseases['omim_id'] == ground_truth_omim_id]
            
            if not rank_info.empty:
                rank = rank_info.index[0] + 1
                score = rank_info['similarity_score'].iloc[0]
            else:
                rank = float('inf')
                score = 0.0
                
            results.append({
                'phenopacket_id': data['id'],
                'ground_truth_omim': ground_truth_omim_id,
                'rank': rank,
                'score': score,
                'num_hpo_terms': len(patient_hpo_ids)
            })
        
    return pd.DataFrame(results)

def calculate_performance_metrics(results_df):
    """
    Calculates overall performance metrics from the validation results.
    
    Args:
        results_df (pandas.DataFrame): The detailed validation results.
        
    Returns:
        dict: A dictionary containing performance metrics.
    """
    total_cases = len(results_df)
    if total_cases == 0:
        return {
            "total_cases_evaluated": 0,
            "message": "No valid cases with both phenotypes and diagnoses were found."
        }

    top1_accuracy = (results_df['rank'] == 1).sum() / total_cases * 100
    top5_accuracy = (results_df['rank'] <= 5).sum() / total_cases * 100
    top10_accuracy = (results_df['rank'] <= 10).sum() / total_cases * 100
    median_rank = results_df['rank'].median()

    metrics = {
        "total_cases_evaluated": total_cases,
        "top1_accuracy_percent": round(top1_accuracy, 2),
        "top5_accuracy_percent": round(top5_accuracy, 2),
        "top10_accuracy_percent": round(top10_accuracy, 2),
        "median_rank": median_rank
    }
    return metrics

def main(phenopacket_dir, output_json_path, similarity_method):
    """
    Main function to run the phenotypic similarity validation pipeline.
    """
    omim_profiles = initialize_hpo_ontology()
    
    validation_results_df = run_validation(phenopacket_dir, omim_profiles, similarity_method)
    
    if validation_results_df.empty:
        logging.warning("Validation produced no results. Exiting.")
        return

    # Save detailed results to JSON
    logging.info(f"Saving detailed validation results to {output_json_path}")
    validation_results_df.to_json(output_json_path, orient='records', indent=4)
    
    # Calculate and display performance metrics
    performance_metrics = calculate_performance_metrics(validation_results_df)
    
    # Print metrics to stdout
    print("\n--- Overall Performance Metrics ---")
    print(json.dumps(performance_metrics, indent=4))
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phenotypic Similarity Validation Script")
    parser.add_argument(
        "--phenopacket_dir",
        type=str,
        default="phenopackets",
        help="Directory containing the PhenoPacket JSON files."
    )
    parser.add_argument(
        "--output_json",
        type=str,
        default="phenotypic_validation_results.json",
        help="Path to save the detailed validation results JSON file."
    )
    parser.add_argument(
        "--similarity_method",
        type=str,
        default="resnik",
        choices=['resnik', 'lin', 'jc', 'rel'],
        help="The similarity algorithm to use. Defaults to 'resnik'."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging (DEBUG level)."
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress info-level logging, showing only warnings and errors."
    )
    parser.add_argument(
        "--log-file",
        type=str,
        default=None,
        help="Path to a file to write logs to. By default, logs are written to stderr."
    )
    
    args = parser.parse_args()

    # --- Configure Logging ---
    log_level = logging.INFO
    if args.verbose:
        log_level = logging.DEBUG
    elif args.quiet:
        log_level = logging.WARNING

    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    
    if args.log_file:
        logging.basicConfig(filename=args.log_file, level=log_level, format=log_format, filemode='w')
    else:
        logging.basicConfig(level=log_level, format=log_format)
    
    main(
        phenopacket_dir=args.phenopacket_dir,
        output_json_path=args.output_json,
        similarity_method=args.similarity_method,
    )
