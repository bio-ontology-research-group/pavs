
import pandas as pd
import os

PROCESSED_DIR = 'data/processed'

def check_file(filename, expected_ids):
    path = os.path.join(PROCESSED_DIR, filename)
    if not os.path.exists(path):
        print(f"FAIL: {filename} missing.")
        return
    
    df = pd.read_csv(path, sep='	')
    print(f"
--- Validating {filename} ---")
    print(f"Rows: {len(df)}, Columns: {len(df.columns)}")
    
    # Check structure
    expected_cols = ['id', 'source_file', 'sex_id', 'hpo_ids', 'gene_symbols']
    missing_cols = [c for c in expected_cols if c not in df.columns]
    if missing_cols:
        print(f"FAIL: Missing columns {missing_cols}")
    else:
        print("PASS: Structure OK")

    # Check content for specific IDs
    for eid in expected_ids:
        row = df[df['id'].astype(str) == str(eid)]
        if row.empty:
            print(f"FAIL: ID {eid} not found in {filename}")
            continue
        
        row = row.iloc[0]
        print(f"Checking ID {eid}:")
        print(f"  Sex: {row['sex_id']} ({row['sex_label']})")
        print(f"  HPO IDs: {row['hpo_ids'][:50]}...")
        print(f"  Genes: {row['gene_symbols']}")
        
        # Validation logic
        success = True
        if pd.isna(row['sex_id']) or 'NCIT' not in str(row['sex_id']):
            if filename != 'ddd_diagnoses_normalized.tsv':
                print(f"  WARN: Sex ID '{row['sex_id']}' not normalized to NCIT")
        
        if filename != 'ddd_diagnoses_normalized.tsv':
            if 'HANCESTRO' not in str(row['ancestry_id']):
                print(f"  FAIL: Ancestry not normalized")
                success = False
        
        if success:
            print(f"  PASS: Content looks valid")

if __name__ == "__main__":
    check_file('ahmed-variants_normalized.tsv', ['1', '2'])
    check_file('ahmed-pmid28454995_normalized.tsv', ['1', '2'])
    check_file('fawzan-variants_normalized.tsv', ['16-2457', '16-2459'])
    check_file('marwa-variants_normalized.tsv', ['PAVS1', 'PAVS2'])
    check_file('PMC6562004_normalized.tsv', ['17-1399', '17-2526'])
    check_file('PMC7082194_normalized.tsv', ['1', '2'])
    check_file('ddd-diagnoses_normalized.tsv', ['DDD_ENTRY_1', 'DDD_ENTRY_2'])
