import pandas as pd
import os

NCIT_MAPPING = {
    'MALE': 'NCIT:C20197',
    'FEMALE': 'NCIT:C16576',
    'UNKNOWN_SEX': 'NCIT:C17998',
    'M': 'NCIT:C20197',
    'F': 'NCIT:C16576',
    'U': 'NCIT:C17998',
    'nan': 'NCIT:C17998'
}

def normalize_sex_value(val):
    val = str(val).upper().strip()
    return NCIT_MAPPING.get(val, 'NCIT:C17998')

files = [
    'data/cleaned_and_normalized/ahmed_normalized.tsv',
    'data/cleaned_and_normalized/fawzan_normalized.tsv',
    'data/cleaned_and_normalized/marwa_normalized.tsv',
    'data/cleaned_and_normalized/ahmed_pmid28454995_normalized.tsv',
    'data/cleaned_and_normalized/pmc7082194_normalized.tsv'
]

for f in files:
    if os.path.exists(f):
        df = pd.read_csv(f, sep='	')
        
        # Identify gender column
        target_col = None
        for col in ['norm_sex', 'Gender', 'Sex', 'normalized_sex']:
            if col in df.columns:
                target_col = col
                break
        
        if target_col:
            df['norm_sex'] = df[target_col].apply(normalize_sex_value)
        else:
            # If no sex column exists, add it as unknown
            df['norm_sex'] = 'NCIT:C17998'
            
        df.to_csv(f, sep='	', index=False)
        print(f"Updated sex to NCIT in {f}")
