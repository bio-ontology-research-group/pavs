import os

replacements = {
    "HANCESTRO:0852": "HANCESTRO:0852",
    "GAZ:00005279": "GAZ:00005279",
    "GENO:0000134": "GENO:0000134",
    "GENO:0000137": "GENO:0000137",
    "GENO:0000511": "GENO:0000511"
}

target_dir = 'data/cleaned_and_normalized'
for filename in os.listdir(target_dir):
    if filename.endswith('.tsv'):
        path = os.path.join(target_dir, filename)
        with open(path, 'r') as f:
            content = f.read()
        
        for old, new in replacements.items():
            content = content.replace(old, new)
            
        with open(path, 'w') as f:
            f.write(content)
        print(f"Patched IDs in {filename}")
