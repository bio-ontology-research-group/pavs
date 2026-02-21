import re
import pandas as pd
import os

def load_valid_genes():
    genes = set()
    path = 'data/genes_to_phenotype.txt'
    if os.path.exists(path):
        # The file might have headers or not, let's be safe
        try:
            df = pd.read_csv(path, sep='	')
            if 'gene_symbol' in df.columns:
                genes = set(df['gene_symbol'].unique())
            else:
                # If no header, maybe it's the 2nd column
                df = pd.read_csv(path, sep='	', header=None)
                genes = set(df[1].unique())
        except Exception as e:
            print(f"Error loading genes: {e}")
    return genes

def load_valid_omims():
    omims = set()
    path = 'data/phenotype.hpoa'
    if os.path.exists(path):
        try:
            with open(path, 'r') as f:
                for line in f:
                    if line.startswith('OMIM:'):
                        omims.add(line.split('	')[0])
        except Exception as e:
            print(f"Error loading OMIMs: {e}")
    return omims

def validate_and_clean_variant(v, valid_genes):
    if not isinstance(v, str) or v == 'nan' or not v:
        return ""
    parts = v.split('|')
    cleaned_parts = []
    for p in parts:
        p = p.strip()
        if not p: continue
        
        # Try to extract gene
        gene_match = re.search(r'^([^:]+)', p)
        gene_symbol = "not known"
        rest = p
        
        if gene_match:
            potential_gene = gene_match.group(1).strip()
            if potential_gene in valid_genes:
                gene_symbol = potential_gene
                rest = p[len(potential_gene):].lstrip(':')
            elif potential_gene.startswith('NM_'):
                # Gene is missing, first part is transcript
                gene_symbol = "not known"
                rest = p
            else:
                # Invalid gene symbol
                gene_symbol = "not known"
                rest = p[len(potential_gene):].lstrip(':')
        
        # Now clean the rest to just Transcript:c.xxx if possible, or just c.xxx
        nm = re.search(r'(NM_\d+(\.\d+)?)', rest)
        c_dna = re.search(r'(c\.[^:\s,;]+)', rest)
        
        if nm and c_dna:
            variant_part = f"{nm.group(1)}:{c_dna.group(1)}"
        elif c_dna:
            variant_part = c_dna.group(1)
        else:
            # If no c.dna found, just keep the rest but maybe it's junk
            variant_part = rest.split(' ')[0].split('(')[0].rstrip(':')
            
        cleaned_parts.append(f"{gene_symbol}:{variant_part}")
    
    return "|".join(cleaned_parts)

def validate_diseases(d, valid_omims):
    if not isinstance(d, str) or not d or d == 'nan':
        return "not known"
    
    d_ids = re.findall(r'(OMIM:\d+)', d)
    valid_found = [oid for oid in d_ids if oid in valid_omims]
    
    if valid_found:
        return "|".join(valid_found)
    else:
        return "not known"

# Mock data for testing
valid_genes = {'CEP152', 'ASXL3', 'COL6A2'}
valid_omims = {'OMIM:243400', 'OMIM:613287'}

test_variants = [
    "CEP152:NM_001194998:exon16:c.2021G>T:p.C674F",
    "INVALIDGENE:NM_001194998:c.2021G>T",
    "NM_030632:c.4117_4118insAGAT",
    "ASXL3:c.4117_4118insAGAT (het)",
    "COL6A2:NM_058175:exon26:c.2422+1G>A | INVALID:c.123A>G"
]

test_diseases = [
    "OMIM:243400",
    "OMIM:999999", # Invalid
    "ORPHA:123",   # Not OMIM
    "OMIM:243400 | ORPHA:456",
    ""
]

print("Testing Variants:")
for v in test_variants:
    print(f"Original: {v}")
    print(f"Cleaned:  {validate_and_clean_variant(v, valid_genes)}")

print("\nTesting Diseases:")
for d in test_diseases:
    print(f"Original: {d}")
    print(f"Cleaned:  {validate_diseases(d, valid_omims)}")
