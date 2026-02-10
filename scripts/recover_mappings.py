
import json
import re
import os

def recover():
    mappings = {}
    
    # Load what we already have
    if os.path.exists('hpo_mappings_partial.json'):
        with open('hpo_mappings_partial.json', 'r') as f:
            mappings.update(json.load(f))
            print(f"Loaded {len(mappings)} from partial file.")

    # Scrape the log for JSON blocks
    if os.path.exists('hpo_mapping_v3.log'):
        with open('hpo_mapping_v3.log', 'r') as f:
            content = f.read()
            # Find everything between [ ] or { } that looks like JSON
            json_blocks = re.findall(r'(\[[\s\S]*?\]|\{[\s\S]*?\})', content)
            
            for block in json_blocks:
                try:
                    data = json.loads(block)
                    if isinstance(data, dict):
                        mappings.update(data)
                    elif isinstance(data, list):
                        for item in data:
                            if isinstance(item, dict):
                                mappings.update(item)
                except:
                    continue
    
    print(f"Total mappings after recovery: {len(mappings)}")
    with open('hpo_mappings_final.json', 'w') as f:
        json.dump(mappings, f, indent=2)

if __name__ == "__main__":
    recover()
