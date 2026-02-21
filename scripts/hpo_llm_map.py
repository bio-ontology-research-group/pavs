import json
import os
import requests
import pandas as pd
from rapidfuzz import process, utils
import time
import re
import sys
from pyhpo import Ontology, HPOSet

# Initialize HPO Ontology
print("Initializing HPO Ontology...")
Ontology()
hpo_choices = [{'name': term.name, 'id': term.id} for term in Ontology]
hpo_names = [x['name'] for x in hpo_choices]

SEVERITY_MAP = {
    "HP:0012825": "Mild",
    "HP:0012826": "Moderate",
    "HP:0012828": "Severe",
    "HP:0012829": "Profound",
    "HP:0012827": "Borderline"
}

def get_candidates(term, context_hpo_set=None, limit=10):
    results = process.extract(term, hpo_names, processor=utils.default_process, limit=limit * 2)
    candidates = []
    for match_name, score, index in results:
        candidates.append({'name': match_name, 'id': hpo_choices[index]['id'], 'score': score})
    
    if context_hpo_set and candidates:
        for cand in candidates:
            try:
                cand_int = int(cand['id'].replace('HP:', ''))
                cand_set = HPOSet.from_queries([cand_int])
                sim = context_hpo_set.similarity(cand_set)
                cand['combined_score'] = cand['score'] + (sim * 50)
            except:
                cand['combined_score'] = cand['score']
        candidates.sort(key=lambda x: x['combined_score'], reverse=True)
    
    return candidates[:limit]

def call_gemini(batch_data):
    api_key = os.getenv("OPENROUTER_API_KEY")
    if not api_key:
        print("Error: OPENROUTER_API_KEY not set.")
        return None
        
    url = "https://openrouter.ai/api/v1/chat/completions"
    
    prompt = """You are a clinical geneticist expert in the Human Phenotype Ontology (HPO). 
Map clinical descriptions to HPO terms, identifying severity and negation.

Output Format:
Return a JSON object where keys are the original terms and values are objects with:
{
  "hpo_id": "HP:XXXXXXX" or "None",
  "severity_id": "HP:XXXXXXX" (one of: HP:0012825 [Mild], HP:0012826 [Moderate], HP:0012828 [Severe], HP:0012829 [Profound], HP:0012827 [Borderline]) or "None",
  "excluded": true/false (true if the text indicates the phenotype is NOT present, e.g., 'no seizures', 'not diabetic')
}

Guidelines:
1. Disambiguation: Use context (Matched HPOs and Genes) to resolve terms like 'ASD'.
2. Severity: Identify words like 'mild', 'severe', 'moderate' and map to the specific severity HPO ID.
3. Negation: If the term is negated (e.g. 'no...', 'not...', 'negative for...'), set "excluded": true and find the HPO ID for the phenotype itself.
4. If a term is not a phenotype, set hpo_id to "None".

"""
    
    for item in batch_data:
        term = item['term']
        candidates = item['candidates']
        context_hpos = item['context_hpos']
        genes = item['genes']
        gene_phenotypes = item['gene_phenotypes']
        
        cand_str = "\n".join([f"- {c['id']}: {c['name']}" for c in candidates])
        context_str = f"Context (Matched HPOs): {', '.join(context_hpos)}\n"
        context_str += f"Context (Genes): {', '.join(genes)}\n"
        if gene_phenotypes:
            context_str += f"Context (Gene phenotypes): {', '.join(gene_phenotypes[:15])}\n"
            
        prompt += f"Term: {term}\n{context_str}Candidates:\n{cand_str}\n\n"

    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json"
    }
    
    data = {
        "model": "google/gemini-2.0-flash-001",
        "messages": [{"role": "user", "content": prompt}],
        "response_format": {"type": "json_object"}
    }
    
    try:
        response = requests.post(url, headers=headers, json=data)
        response_json = response.json()
        if 'choices' in response_json:
            return response_json['choices'][0]['message']['content']
        else:
            print(f"API Error: {response_json}")
            return None
    except Exception as e:
        print(f"Error: {e}")
        return None

def main():
    if not os.path.exists('master_data_v1.tsv'):
        print("Error: master_data_v1.tsv not found.")
        return

    df = pd.read_csv('master_data_v1.tsv', sep='\t')
    
    patient_specific_mappings = {}
    if os.path.exists('patient_hpo_mappings.json'):
        with open('patient_hpo_mappings.json', 'r') as f:
            patient_specific_mappings = json.load(f)

    items_to_map = []
    
    for _, row in df.iterrows():
        p_id = str(row['Master_ID'])
        unmatched_str = str(row['Unmatched_Phenotypes'])
        if not unmatched_str or unmatched_str == 'nan':
            continue
            
        unmatched_terms = [t.strip() for t in unmatched_str.split('|') if t.strip()]
        
        context_hpo_ids = str(row['Phenotypes_HPO']).split('|') if pd.notna(row['Phenotypes_HPO']) else []
        context_hpo_ints = []
        context_hpo_names = []
        for h_id in context_hpo_ids:
            try:
                h_int = int(h_id.replace('HP:', ''))
                context_hpo_ints.append(h_int)
                context_hpo_names.append(Ontology[h_int].name)
            except: pass
        
        context_hpo_set = HPOSet.from_queries(context_hpo_ints) if context_hpo_ints else None
        
        # Get gene info from Variants
        genes = []
        if pd.notna(row['Variants']):
            for v in str(row['Variants']).split('|'):
                m = re.match(r'^([^:]+)', v.strip())
                if m: genes.append(m.group(1).strip())
        
        gene_phenotypes = []
        for symbol in genes:
            try:
                for g in Ontology.genes:
                    if g.name == symbol:
                        gene_phenotypes.extend([Ontology[h].name for h in g.hpo])
                        break
            except: pass

        for term in unmatched_terms:
            if p_id in patient_specific_mappings and term in patient_specific_mappings[p_id]:
                # If existing mapping is old style (just a string), we should re-map it
                if isinstance(patient_specific_mappings[p_id][term], str):
                    pass 
                else:
                    continue
            
            candidates = get_candidates(term, context_hpo_set)
            if not candidates: continue

            items_to_map.append({
                'patient_id': p_id,
                'term': term,
                'candidates': candidates,
                'context_hpos': context_hpo_names,
                'genes': list(set(genes)),
                'gene_phenotypes': list(set(gene_phenotypes))
            })

    print(f"Total patient-term pairs to map: {len(items_to_map)}")
    
    batch_size = 5
    for i in range(0, len(items_to_map), batch_size):
        batch = items_to_map[i:i+batch_size]
        print(f"Processing batch {i//batch_size + 1}/{(len(items_to_map)//batch_size)+1}...")
        
        res = call_gemini(batch)
        if res:
            try:
                res_clean = re.sub(r'```json\n|\n```', '', res).strip()
                batch_map = json.loads(res_clean)
                for item in batch:
                    term = item['term']
                    p_id = item['patient_id']
                    if term in batch_map:
                        if p_id not in patient_specific_mappings:
                            patient_specific_mappings[p_id] = {}
                        patient_specific_mappings[p_id][term] = batch_map[term]
            except Exception as e:
                print(f"Failed to parse JSON: {e}")
        
        with open('patient_hpo_mappings.json', 'w') as f:
            json.dump(patient_specific_mappings, f, indent=2)
        time.sleep(1)

    print("Mapping complete!")

if __name__ == "__main__":
    main()
