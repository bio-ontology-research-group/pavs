import pandas as pd
import re
import ast
import os

def load_hpo_ids(obo_path):
    if not os.path.exists(obo_path):
        print(f"Warning: {obo_path} not found.")
        return set()
    ids = set()
    with open(obo_path, 'r') as f:
        for line in f:
            if line.startswith('id: HP:'):
                ids.add(line.split('id: ')[1].strip())
    return ids

def verify_hgvs(hgvs_list):
    if not hgvs_list: return True, []
    invalid = []
    prefixes = (r'c\.', r'g\.', r'm\.', r'n\.', r'p\.', r'NC_', r'NM_', r'NP_', r'NG_', r'NR_', r'LRG_')
    pattern = '^(' + '|'.join(prefixes) + ')'
    for h in hgvs_list:
        if not re.match(pattern, str(h)):
            invalid.append(h)
    return len(invalid) == 0, invalid

def main():
    hpo_ids = load_hpo_ids('ontology/hp.obo')
    data_path = 'data/cleaned/all_cleaned_normalized.tsv'
    if not os.path.exists(data_path):
        print(f"Error: {data_path} not found.")
        return

    df = pd.read_csv(data_path, sep='\t')
    
    # Converters
    for col in ['phenotypes_high_conf', 'phenotypes_low_conf', 'genes', 'variants_hgvs', 'diseases_high_conf']:
        df[col] = df[col].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) and x.startswith('[') else x)

    semantic_mismatches = 0
    issues = []
    
    for i, row in df.iterrows():
        # 1. HPO Validity
        p_high = row['phenotypes_high_conf']
        for h in p_high:
            if h not in hpo_ids:
                issues.append(f"Row {i}: Invalid HPO ID {h}")

        # 2. Semantic Count Check
        # Extract the estimated count from justification
        m = re.search(r'(\d+)', str(row['phenotype_count_justification']))
        if m:
            est_count = int(m.group(1))
            actual_extracted = len(p_high) + len(row['phenotypes_low_conf'])
            if actual_extracted > est_count + 2: # Allow some buffer for splitting
                semantic_mismatches += 1
                if semantic_mismatches < 10:
                    issues.append(f"Row {i}: Semantic Mismatch! Est: {est_count}, Actual: {actual_extracted}. Justification: {row['phenotype_count_justification']}")

        # 3. HGVS
        ok, invalid = verify_hgvs(row['variants_hgvs'])
        if not ok:
            issues.append(f"Row {i}: Possible invalid HGVS: {invalid}")

        # 4. Disease IDs
        for d in row['diseases_high_conf']:
            if not (d.startswith('OMIM:') or d.startswith('ORPHA:')):
                issues.append(f"Row {i}: Invalid Disease ID {d}")

    print(f"Verification completed for {len(df)} rows.")
    print(f"Found {len(issues)} specific issues and {semantic_mismatches} total semantic count mismatches.")
    for issue in issues[:20]:
        print(issue)

if __name__ == "__main__":
    main()
