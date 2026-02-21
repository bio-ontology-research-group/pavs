from pyhpo import Ontology, HPOSet
import os

# Initialize Ontology
print("Loading HPO Ontology (this takes ~15-20s)...")
try:
    Ontology()
    print(f"Ontology loaded. {len(Ontology)} terms found.")
except Exception as e:
    print(f"Failed to load ontology: {e}")

def get_hpo_set(hpo_ids):
    if not hpo_ids:
        return None
    try:
        if isinstance(hpo_ids, str):
            if not hpo_ids: return None
            hpo_ids = [x.strip() for x in hpo_ids.split(",") if x.strip()]
        if not hpo_ids: return None
        return HPOSet.from_queries(hpo_ids)
    except Exception as e:
        return None

def calculate_similarity(query_hpo_ids, target_hpo_set, method="resnik"):
    try:
        if not target_hpo_set or not query_hpo_ids:
            return 0.0
        query_set = HPOSet.from_queries(query_hpo_ids)
        return query_set.similarity(target_hpo_set, method=method)
    except Exception as e:
        return 0.0

def search_hpo_terms(query, limit=10):
    results = []
    query = query.lower()
    count = 0
    for term in Ontology:
        if query in term.name.lower() or query in term.id.lower():
            results.append({"id": term.id, "name": term.name})
            count += 1
            if count >= limit:
                break
    return results