import pandas as pd
import json
import re

# ── 1. Load HPOA ─────────────────────────────────────────────────────────────
print("Loading HPOA...")
hpoa = pd.read_csv('phenotype.hpoa', sep='\t', comment='#', low_memory=False,
                   names=['DatabaseID', 'DiseaseName', 'Qualifier',
                          'HPO_ID', 'Reference', 'Evidence',
                          'Onset', 'Frequency', 'Sex', 'Modifier',
                          'Aspect', 'Biocuration'])

def extract_pmids_from_series(series):
    pmids = set()
    for val in series.dropna():
        found = re.findall(r'PMID:(\d+)', str(val), re.IGNORECASE)
        pmids.update(found)
    return pmids

hpoa_pmids = extract_pmids_from_series(hpoa['Reference'])
print(f"Unique PMIDs in HPOA: {len(hpoa_pmids)}")

# ── 2. Load PAVS cases TSV for source labels ──────────────────────────────────
df = pd.read_csv('PAVS_cases.tsv', sep='\t')
saudi_ids = set(df[df['source'] == 'PAVS-Saudi']['case_id'])
ddd_ids   = set(df[df['source'] == 'DDD']['case_id'])
mixed_ids = set(df[df['source'] == 'PAVS-mixed']['case_id'])

print(f"\nSaudi cases: {len(saudi_ids)}")
print(f"DDD cases: {len(ddd_ids)}")
print(f"Mixed cases: {len(mixed_ids)}")

# ── 3. Extract PMIDs from phenopackets JSON ───────────────────────────────────
print("\nLoading PAVS_phenopackets.json...")
with open('PAVS_phenopackets.json') as f:
    data = json.load(f)

saudi_pmids  = set()
ddd_pmids    = set()
mixed_pmids  = set()

saudi_case_pmid_list = []
ddd_case_pmid_list   = []

for case in data:
    case_id = case.get('id', '')
    refs = case.get('metaData', {}).get('externalReferences', [])
    
    # Extract PMIDs from externalReferences
    case_pmids = set()
    for ref in refs:
        ref_id = ref.get('id', '')
        match = re.search(r'PMID:(\d+)', ref_id, re.IGNORECASE)
        if match:
            case_pmids.add(match.group(1))
    
    # Classify by source
    if case_id in saudi_ids:
        saudi_pmids.update(case_pmids)
        saudi_case_pmid_list.extend(list(case_pmids))
    elif case_id in ddd_ids:
        ddd_pmids.update(case_pmids)
        ddd_case_pmid_list.extend(list(case_pmids))
    elif case_id in mixed_ids:
        mixed_pmids.update(case_pmids)

print(f"Unique PMIDs in Saudi cases: {len(saudi_pmids)}")
print(f"Unique PMIDs in DDD cases: {len(ddd_pmids)}")
print(f"Unique PMIDs in Mixed cases: {len(mixed_pmids)}")

# ── 4. Extract PMIDs from literature.ttl ─────────────────────────────────────
print("\nExtracting PMIDs from literature.ttl...")
lit_pmids = set()
lit_case_pmid_list = []

with open('literature.ttl', 'r') as f:
    for line in f:
        if 'pav:lit_PMID_' in line:
            found = re.findall(r'pav:lit_PMID_(\d+)', line)
            for pmid in found:
                lit_pmids.add(pmid)
                lit_case_pmid_list.append(pmid)

print(f"Unique PMIDs in Literature cohort: {len(lit_pmids)}")
print(f"Total Literature case-PMID entries: {len(lit_case_pmid_list)}")

# ── 5. Calculate PMID-level overlap ──────────────────────────────────────────
lit_overlap   = lit_pmids   & hpoa_pmids
ddd_overlap   = ddd_pmids   & hpoa_pmids
saudi_overlap = saudi_pmids & hpoa_pmids

print(f"\n=== PMID-LEVEL OVERLAP ===")
print(f"Literature : {len(lit_overlap):4d} / {len(lit_pmids):4d} PMIDs in HPOA "
      f"({len(lit_overlap)/max(len(lit_pmids),1)*100:.1f}%)")
print(f"DDD        : {len(ddd_overlap):4d} / {len(ddd_pmids):4d} PMIDs in HPOA "
      f"({len(ddd_overlap)/max(len(ddd_pmids),1)*100:.1f}%)")
print(f"Saudi      : {len(saudi_overlap):4d} / {len(saudi_pmids):4d} PMIDs in HPOA "
      f"({len(saudi_overlap)/max(len(saudi_pmids),1)*100:.1f}%)")

# ── 6. Case-level overlap ─────────────────────────────────────────────────────
lit_cases_overlap   = sum(1 for p in lit_case_pmid_list   if p in hpoa_pmids)
ddd_cases_overlap   = sum(1 for p in ddd_case_pmid_list   if p in hpoa_pmids)
saudi_cases_overlap = sum(1 for p in saudi_case_pmid_list if p in hpoa_pmids)

print(f"\n=== CASE-LEVEL OVERLAP (cases whose PMID appears in HPOA) ===")
print(f"Literature : {lit_cases_overlap:4d} / {len(lit_case_pmid_list):4d} cases "
      f"({lit_cases_overlap/max(len(lit_case_pmid_list),1)*100:.1f}%)")
print(f"DDD        : {ddd_cases_overlap:4d} / {len(ddd_case_pmid_list):4d} cases "
      f"({ddd_cases_overlap/max(len(ddd_case_pmid_list),1)*100:.1f}%)")
print(f"Saudi      : {saudi_cases_overlap:4d} / {len(saudi_case_pmid_list):4d} cases "
      f"({saudi_cases_overlap/max(len(saudi_case_pmid_list),1)*100:.1f}%)")

results = pd.DataFrame([
    {
        'cohort': 'Literature',
        'total_pmids': len(lit_pmids),
        'pmids_in_hpoa': len(lit_overlap),
        'pmid_overlap_pct': f"{len(lit_overlap)/max(len(lit_pmids),1)*100:.1f}%",
        'total_cases': len(lit_case_pmid_list),
        'cases_with_pmid_in_hpoa': lit_cases_overlap,
        'case_overlap_pct': f"{lit_cases_overlap/max(len(lit_case_pmid_list),1)*100:.1f}%"
    },
    {
        'cohort': 'DDD',
        'total_pmids': len(ddd_pmids),
        'pmids_in_hpoa': len(ddd_overlap),
        'pmid_overlap_pct': f"{len(ddd_overlap)/max(len(ddd_pmids),1)*100:.1f}%",
        'total_cases': len(ddd_case_pmid_list),
        'cases_with_pmid_in_hpoa': ddd_cases_overlap,
        'case_overlap_pct': f"{ddd_cases_overlap/max(len(ddd_case_pmid_list),1)*100:.1f}%"
    },
    {
        'cohort': 'Saudi',
        'total_pmids': len(saudi_pmids),
        'pmids_in_hpoa': len(saudi_overlap),
        'pmid_overlap_pct': f"{len(saudi_overlap)/max(len(saudi_pmids),1)*100:.1f}%",
        'total_cases': len(saudi_case_pmid_list),
        'cases_with_pmid_in_hpoa': saudi_cases_overlap,
        'case_overlap_pct': f"{saudi_cases_overlap/max(len(saudi_case_pmid_list),1)*100:.1f}%"
    }
])

results.to_csv('hpoa_overlap_results.csv', index=False)
print("\nResults saved to hpoa_overlap_results.csv")
print(results.to_string(index=False))

with open('lit_overlapping_pmids.txt', 'w') as f:
    for pmid in sorted(lit_overlap):
        f.write(f"PMID:{pmid}\n")

with open('ddd_overlapping_pmids.txt', 'w') as f:
    for pmid in sorted(ddd_overlap):
        f.write(f"PMID:{pmid}\n")

print("Overlapping PMIDs saved to lit_overlapping_pmids.txt and ddd_overlapping_pmids.txt")