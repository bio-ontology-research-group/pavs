import pandas as pd
df = pd.read_csv('data/phenotypes/ahmed-pmid28454995.tsv', sep='	')
for i, col in enumerate(df.columns):
    print(f"{i}: {repr(col)}")
