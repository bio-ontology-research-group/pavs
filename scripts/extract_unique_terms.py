
import json

with open('ambiguity_queue.json', 'r') as f:
    queue = json.load(f)

unique_terms = set()
for item in queue:
    for term in item['terms']:
        unique_terms.add(term.strip())

print(f"Total unique unmatched terms: {len(unique_terms)}")
with open('unique_unmatched_terms.json', 'w') as f:
    json.dump(list(unique_terms), f, indent=2)
