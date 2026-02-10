
import json
import os
import requests
import pandas as pd
from rapidfuzz import process, utils
import time
import re
import sys

# Load HPO
def parse_hpo_obo(obo_path):
    hpo_list = [] # List of (name, id)
    with open(obo_path, 'r') as f:
        current_id = None
        current_name = None
        for line in f:
            line = line.strip()
            if line.startswith('id: HP:'):
                current_id = line.split('id: ')[1].split(' ! ')[0]
            elif line.startswith('name: '):
                current_name = line.split('name: ')[1]
                hpo_list.append({'name': current_name, 'id': current_id})
            elif line.startswith('synonym: '):
                synonym = re.findall(r'\"(.+?)\"', line)
                if synonym:
                    hpo_list.append({'name': synonym[0], 'id': current_id})
    return hpo_list

def get_candidates(term, hpo_list, limit=10):
    choices = [x['name'] for x in hpo_list]
    results = process.extract(term, choices, processor=utils.default_process, limit=limit)
    candidates = []
    for match_name, score, index in results:
        candidates.append({'name': match_name, 'id': hpo_list[index]['id'], 'score': score})
    return candidates

def call_gemini(batch, hpo_list):
    api_key = os.getenv("OPENROUTER_API_KEY")
    url = "https://openrouter.ai/api/v1/chat/completions"
    
    prompt = """You are a clinical geneticist. Map the following clinical descriptions to the most appropriate Human Phenotype Ontology (HPO) ID from the provided candidates. 

Guidelines:
1. Select the ID that most specifically represents the phenotypic abnormality described.
2. If the description contains multiple distinct phenotypes, select the most prominent one or the one that encompasses the others.
3. IMPORTANT: If the description is about family history (e.g., 'no similar FH', 'consanguineous') or is not a phenotypic abnormality, return 'None'.
4. If none of the candidates are a good match, return 'None'.

"""
    
    for term in batch:
        candidates = get_candidates(term, hpo_list)
        cand_str = "\n".join([f"- {c['id']}: {c['name']}" for c in candidates])
        prompt += f"Term: {term}\nCandidates:\n{cand_str}\nSelected ID: \n\n"

    prompt += "Return ONLY a JSON object where keys are terms and values are the selected HPO ID (e.g. {'Term A': 'HP:0001234'})."

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
    hpo_list = parse_hpo_obo('ontology/hp.obo')
    with open('unique_unmatched_terms.json', 'r') as f:
        terms = json.load(f)
    
    mappings = {}
    if os.path.exists('hpo_mappings_partial.json'):
        with open('hpo_mappings_partial.json', 'r') as f:
            mappings = json.load(f)
    
    # Filter terms already mapped
    terms_to_map = [t for t in terms if t not in mappings]
    print(f"Resuming mapping. {len(terms_to_map)} terms left out of {len(terms)}.")
    
    batch_size = 10
    for i in range(0, len(terms_to_map), batch_size):
        batch = terms_to_map[i:i+batch_size]
        print(f"Processing batch {i//batch_size + 1}/{(len(terms_to_map)//batch_size)+1}...")
        sys.stdout.flush()
        
        res = call_gemini(batch, hpo_list)
        if res:
            try:
                # Clean potential markdown from response
                res_clean = re.sub(r'```json\n|\n```', '', res).strip()
                batch_map = json.loads(res_clean)
                mappings.update(batch_map)
            except Exception as e:
                print(f"Failed to parse JSON for batch {i}: {e}. Response was: {res}")
        
        # Save progress every batch for safety
        with open('hpo_mappings_partial.json', 'w') as f:
            json.dump(mappings, f, indent=2)
            
        time.sleep(0.5)

    with open('hpo_mappings_final.json', 'w') as f:
        json.dump(mappings, f, indent=2)
    print("Mapping complete!")

if __name__ == "__main__":
    main()
