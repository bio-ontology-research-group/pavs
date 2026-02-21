import pandas as pd
import re
import os
import sys
import json
import requests
from rapidfuzz import process, utils
import networkx as nx

# --- Configuration ---
OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY")
MODEL_NAME = "google/gemini-2.0-flash-001"

# --- Load Lookups ---

def load_hpo(obo_path):
    hpo_graph = nx.DiGraph()
    id_to_name, name_to_ids = {}, {}
    if not os.path.exists(obo_path): return hpo_graph, id_to_name, name_to_ids
    with open(obo_path, 'r') as f:
        current_id = None
        for line in f:
            line = line.strip()
            if line.startswith('id: HP:'):
                current_id = line.split('id: ')[1].split(' ! ')[0]
                hpo_graph.add_node(current_id)
            elif line.startswith('name: ') and current_id:
                name = line.split('name: ')[1]
                id_to_name[current_id] = name
                name_to_ids.setdefault(name.lower(), []).append(current_id)
            elif line.startswith('is_a: HP:') and current_id:
                parent_id = line.split('is_a: ')[1].split(' ! ')[0]
                hpo_graph.add_edge(parent_id, current_id) 
            elif line.startswith('synonym: ') and current_id:
                synonym = re.findall(r'"(.+?)"', line)
                if synonym: name_to_ids.setdefault(synonym[0].lower(), []).append(current_id)
    return hpo_graph, id_to_name, name_to_ids

def load_genes(gene_info_path):
    gene_to_id = {}
    if not os.path.exists(gene_info_path): return gene_to_id
    df = pd.read_csv(gene_info_path, sep='\t')
    for _, row in df.iterrows():
        symbol = str(row['Symbol'])
        gene_id = str(row['GeneID'])
        gene_to_id[symbol.lower()] = gene_id
        synonyms = str(row['Synonyms']).split('|')
        for s in synonyms:
            if s != '-': gene_to_id[s.lower()] = gene_id
    return gene_to_id

def load_disease_labels(hpoa_path):
    id_to_name, name_to_id = {}, {}
    if not os.path.exists(hpoa_path): return id_to_name, name_to_id
    df = pd.read_csv(hpoa_path, sep='\t', comment='#')
    for _, row in df.iterrows():
        d_id = str(row['database_id'])
        d_name = str(row['disease_name'])
        id_to_name[d_id] = d_name
        name_to_id[d_name.lower()] = d_id
    return id_to_name, name_to_id

# --- Initialization ---
print("Loading ontologies...")
HPO_GRAPH, HPO_ID_TO_NAME, HPO_NAME_TO_IDS = load_hpo('ontology/hp.obo')
GENE_TO_ID = load_genes('data/Homo_sapiens.gene_info')
DISEASE_ID_TO_NAME, DISEASE_NAME_TO_ID = load_disease_labels('data/phenotype.hpoa')

GENO_MAPPING = {
    'heterozygous': 'GENO:0000135', 'homozygous': 'GENO:0000136', 'hemizygous': 'GENO:0000134',
    'het': 'GENO:0000135', 'hom': 'GENO:0000136', 'hemi': 'GENO:0000134',
    'homo': 'GENO:0000136', 'unknown': 'GENO:0000137'
}

# --- Helpers ---

def get_hpo_candidates(query_term, limit=5):
    all_names = list(HPO_NAME_TO_IDS.keys())
    matches = process.extract(query_term.lower(), all_names, limit=limit * 2)
    candidates = []
    seen_ids = set()
    for match_name, score, index in matches:
        for hpo_id in HPO_NAME_TO_IDS[match_name]:
            if hpo_id not in seen_ids:
                name = HPO_ID_TO_NAME.get(hpo_id, "Unknown")
                candidates.append({'id': hpo_id, 'name': name, 'score': score})
                seen_ids.add(hpo_id)
            if len(candidates) >= limit: break
        if len(candidates) >= limit: break
    return candidates

def llm_choose_hpo(query_term, candidates):
    if not OPENROUTER_API_KEY or not candidates: return candidates[0]['id'] if candidates else None
    cand_str = "\n".join([f"{i+1}. {c['id']} ({c['name']})" for i, c in enumerate(candidates)])
    prompt = f"Map clinical term: '{query_term}'.\nCandidates:\n{cand_str}\nRespond with ONLY the ID or 'None'."
    try:
        response = requests.post("https://openrouter.ai/api/v1/chat/completions",
            headers={"Authorization": f"Bearer {OPENROUTER_API_KEY}"},
            json={"model": MODEL_NAME, "messages": [{"role": "user", "content": prompt}]}, timeout=10)
        res = response.json()['choices'][0]['message']['content'].strip()
        match = re.search(r'HP:\d{7}', res)
        return match.group(0) if match else None
    except: return candidates[0]['id']

def llm_split_phenotypes(text):
    if not OPENROUTER_API_KEY: return [t.strip() for t in text.split(',') if t.strip()]
    prompt = f"Split into distinct phenotypic features. Text: \"{text}\". Return ONLY a JSON list of strings."
    try:
        response = requests.post("https://openrouter.ai/api/v1/chat/completions",
            headers={"Authorization": f"Bearer {OPENROUTER_API_KEY}"},
            json={"model": MODEL_NAME, "messages": [{"role": "user", "content": prompt}], "response_format": {"type": "json_object"}}, timeout=10)
        content = response.json()['choices'][0]['message']['content']
        match = re.search(r'\[.*\]', content, re.DOTALL)
        return json.loads(match.group(0)) if match else [text]
    except: return [t.strip() for t in text.split(',') if t.strip()]

def normalize_ahmed_pmid():
    input_file = 'data/phenotypes/ahmed-pmid28454995.tsv'
    output_file = 'data/cleaned_and_normalized/ahmed_pmid28454995_normalized.tsv'
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    df = pd.read_csv(input_file, sep='\t')
    
    results = []
    for i, row in df.iterrows():
        # 1. Phenotypes
        diag = str(row['Diagnosis']) if not pd.isna(row['Diagnosis']) else ""
        add_p = str(row['Additional clinical phenotype']) if not pd.isna(row['Additional clinical phenotype']) else ""
        full_text = f"{diag}, {add_p}".strip(", ")
        
        p_terms = llm_split_phenotypes(full_text)
        f_hpos, f_labels = [], []
        for term in p_terms:
            if term.lower() in HPO_NAME_TO_IDS:
                h_id = HPO_NAME_TO_IDS[term.lower()][0]
                f_hpos.append(h_id); f_labels.append(HPO_ID_TO_NAME[h_id])
            else:
                cands = get_hpo_candidates(term)
                h_id = (cands[0]['id'] if cands and cands[0]['score'] > 98 else llm_choose_hpo(term, cands))
                if h_id: f_hpos.append(h_id); f_labels.append(HPO_ID_TO_NAME.get(h_id, "Unknown"))

        # 2. Variants
        norm_variants = []
        norm_genes = []
        
        # Variant 1
        g1 = str(row['Variant 1 Gene']).strip()
        t1 = str(row['Transcript ID']).strip()
        c1 = str(row['c.DNA']).strip()
        z1 = str(row['Variant 1 Allele State']).lower().strip()
        
        # Variant 2 (note double space and leading space from inspection)
        g2 = str(row['Variant 2  Gene']).strip()
        t2 = str(row['Transcript ID.1']).strip()
        c2 = str(row[' c.DNA']).strip()
        z2 = str(row['Variant 2 Allele State']).lower().strip()

        if g1 and g1 != 'nan' and c1 and c1 != 'nan':
            if not c1.startswith('c.'): c1 = f"c.{c1}"
            norm_variants.append(f"{t1}({g1}):{c1}" if t1 and t1 != 'nan' else f"{g1}:{c1}")
            norm_genes.append(g1)
            
        if g2 and g2 != 'nan' and c2 and c2 != 'nan':
            if not c2.startswith('c.'): c2 = f"c.{c2}"
            norm_variants.append(f"{t2}({g2}):{c2}" if t2 and t2 != 'nan' else f"{g2}:{c2}")
            norm_genes.append(g2)

        # 3. Zygosity & Digenic logic
        overall_zyg = "GENO:0000137"
        if len(norm_genes) == 1:
            overall_zyg = GENO_MAPPING.get(z1, "GENO:0000137")
            if len(norm_variants) > 1: overall_zyg = "GENO:0000402" # Compound Het
        elif len(norm_genes) > 1:
            overall_zyg = "GENO:0000135" # Digenic represented as multiple hets or complex

        # 4. OMIM
        omim_raw = str(row['OMIM'])
        norm_diseases = []
        if omim_raw and omim_raw != 'nan':
            for o_id in re.findall(r'\d{6}', omim_raw):
                norm_diseases.append(f"OMIM:{o_id}")

        entry = row.to_dict()
        entry.update({
            'normalized_hpos': f_hpos,
            'normalized_hpo_labels': f_labels,
            'normalized_variants_hgvs': norm_variants,
            'normalized_genes_entrez': [GENE_TO_ID.get(g.lower(), "Unknown") for g in norm_genes],
            'normalized_zygosity_geno': overall_zyg,
            'normalized_diseases': norm_diseases,
            'normalized_disease_labels': [DISEASE_ID_TO_NAME.get(d, "Unknown") for d in norm_diseases],
            'reference': 'PMID:28454995',
            'pathogenicity': 'PATHOGENIC',
            'norm_ancestry': 'HANCESTRO:0852',
            'norm_country': 'GAZ:00005279'
        })
        results.append(entry)
        if i % 20 == 0: print(f"Processed {i} rows...")

    pd.DataFrame(results).to_csv(output_file, sep='\t', index=False)
    print(f"Ahmed PMID Normalization complete. Output: {output_file}")

if __name__ == "__main__":
    normalize_ahmed_pmid()
