# Import necessary libraries
import os
import json
import argparse
import logging
import pandas as pd
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from pyhpo import Ontology, HPOSet
from sklearn.metrics import roc_auc_score


def load_gene_hpo_profiles():
    """
    Initializes the HPO Ontology and creates a mapping from gene symbols to their
    associated HPO term sets.
    """
    logging.info("Initializing HPO Ontology... (This may download files on first run)")
    _ = Ontology()
    logging.info("Ontology initialized successfully.")
    
    gene_profiles = {}
    logging.info("Building gene-to-phenotype profiles...")
    for gene in tqdm(Ontology.genes, desc="Processing Genes"):
        if not gene.hpo:
            continue
        
        # Convert integer HPO IDs from the gene object to HP:XXXXXXX string format
        hpo_ids_as_str = [Ontology[hpo_id].id for hpo_id in gene.hpo]
        
        # Create an HPOSet for similarity calculations
        gene_profiles[gene.name] = HPOSet.from_queries(hpo_ids_as_str)
        
    logging.info(f"Built profiles for {len(gene_profiles)} genes with phenotype associations.")
    return gene_profiles

def calculate_similarity_ranking(patient_hpo_set, entity_profiles, similarity_method='resnik'):
    """
    Calculates the similarity between a patient's HPOSet and all entity profiles (e.g., genes),
    returning a ranked list of entities.
    
    Args:
        patient_hpo_set (HPOSet): The set of HPO terms for the patient.
        entity_profiles (dict): A dictionary of entity IDs (e.g., gene symbols) to HPOSets.
        similarity_method (str): The similarity algorithm to use ('resnik', 'lin', 'jc', 'rel').
        
    Returns:
        pandas.DataFrame: A DataFrame with entity IDs and their similarity scores,
                          sorted in descending order.
    """
    if not patient_hpo_set:
        return pd.DataFrame(columns=['entity_id', 'similarity_score'])
        
    results = []
    for entity_id, entity_hpo_set in entity_profiles.items():
        score = patient_hpo_set.similarity(entity_hpo_set, method=similarity_method)
        results.append({'entity_id': entity_id, 'similarity_score': score})
        
    ranked_df = pd.DataFrame(results)
    ranked_df = ranked_df.sort_values(by='similarity_score', ascending=False).reset_index(drop=True)
    
    return ranked_df

def run_validation(phenopacket_dir, gene_profiles, similarity_method):
    """
    Iterates through all PhenoPackets, performs similarity analysis against gene profiles,
    and evaluates the ranking of the true causal gene.
    
    Args:
        phenopacket_dir (str): Path to the directory containing PhenoPacket JSON files.
        gene_profiles (dict): A dictionary of gene symbols to HPOSets.
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
                
            # --- Extract Patient Phenotypes ---
            patient_hpo_ids = [pf['type']['id'] for pf in data.get('phenotypicFeatures', [])]
            if not patient_hpo_ids:
                logging.warning(f"Skipping {data['id']}: No phenotypic features found.")
                continue
            patient_hpo_set = HPOSet.from_queries(patient_hpo_ids)
            
            # --- Extract Ground Truth Gene ---
            ground_truth_gene = None
            interpretations = data.get('interpretations', [])
            if interpretations:
                # Path: Interpretation -> diagnosis -> genomicInterpretations -> variantInterpretation -> variationDescriptor -> geneContext
                try:
                    gene_context = interpretations[0]['diagnosis']['genomicInterpretations'][0]['variantInterpretation']['variationDescriptor']['geneContext']
                    ground_truth_gene = gene_context.get('symbol')
                except (KeyError, IndexError):
                    pass  # Will be caught by the check below

            if ground_truth_gene is None:
                logging.warning(f"Skipping phenopacket '{data['id']}': Could not extract ground truth gene symbol from interpretations.")
                continue

            if not ground_truth_gene:
                logging.warning(f"Skipping phenopacket '{data['id']}': Ground truth gene symbol is present but empty.")
                continue

            if ground_truth_gene not in gene_profiles:
                logging.warning(f"Skipping {data['id']}: Ground truth gene '{ground_truth_gene}' not in HPO gene profiles.")
                continue

            # --- Rank all genes by phenotype similarity ---
            ranked_genes_df = calculate_similarity_ranking(patient_hpo_set, gene_profiles, similarity_method)
            
            # --- Find Rank of the correct gene ---
            rank_info = ranked_genes_df[ranked_genes_df['entity_id'] == ground_truth_gene]
            rank = rank_info.index[0] + 1 if not rank_info.empty else float('inf')
            score = rank_info['similarity_score'].iloc[0] if not rank_info.empty else 0.0

            # --- Calculate ROC AUC for this ranking ---
            y_true = (ranked_genes_df['entity_id'] == ground_truth_gene).astype(int).tolist()
            y_score = ranked_genes_df['similarity_score'].tolist()
            
            roc_auc = float('nan')
            if len(set(y_true)) > 1: # Requires both positive and negative samples
                roc_auc = roc_auc_score(y_true, y_score)
            else:
                logging.warning(f"Cannot calculate ROC AUC for {data['id']}: only one class present.")

            results.append({
                'phenopacket_id': data['id'],
                'ground_truth_gene': ground_truth_gene,
                'rank': rank,
                'score': score,
                'num_hpo_terms': len(patient_hpo_ids),
                'roc_auc': roc_auc
            })
        
    return pd.DataFrame(results)

def calculate_performance_metrics(results_df):
    """
    Calculates overall performance metrics (Hits@k, ROC AUC) from the validation results.
    
    Args:
        results_df (pandas.DataFrame): The detailed validation results.
        
    Returns:
        dict: A dictionary containing performance metrics.
    """
    total_cases = len(results_df)
    if total_cases == 0:
        return {
            "total_cases_evaluated": 0,
            "message": "No valid cases were found to evaluate."
        }

    # Calculate Hits@k
    hits_at_1 = (results_df['rank'] == 1).sum() / total_cases
    hits_at_10 = (results_df['rank'] <= 10).sum() / total_cases
    hits_at_100 = (results_df['rank'] <= 100).sum() / total_cases
    
    # Calculate Mean ROC AUC, ignoring NaNs
    mean_roc_auc = results_df['roc_auc'].mean()
    median_rank = results_df['rank'].median()

    metrics = {
        "total_cases_evaluated": total_cases,
        "hits_at_1": round(hits_at_1, 4),
        "hits_at_10": round(hits_at_10, 4),
        "hits_at_100": round(hits_at_100, 4),
        "mean_roc_auc": round(mean_roc_auc, 4),
        "median_rank": median_rank
    }
    return metrics

def main(phenopacket_dir, output_json_path, similarity_method):
    """
    Main function to run the phenotypic similarity validation pipeline.
    """
    gene_profiles = load_gene_hpo_profiles()
    
    if not gene_profiles:
        logging.error("Gene profiles could not be loaded. Exiting.")
        return

    validation_results_df = run_validation(phenopacket_dir, gene_profiles, similarity_method)
    
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
