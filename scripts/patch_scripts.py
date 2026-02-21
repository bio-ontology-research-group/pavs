import os

replacements = {
    "HANCESTRO:0852": "HANCESTRO:0852",
    "GAZ:00005279": "GAZ:00005279",
    "GENO:0000134": "GENO:0000134",
    "GENO:0000137": "GENO:0000137",
    "GENO:0000511": "GENO:0000511",
    "HP:0025820": "HP:0025820"
}

script_dir = 'scripts'
for filename in os.listdir(script_dir):
    if filename.endswith('.py'):
        path = os.path.join(script_dir, filename)
        with open(path, 'r') as f:
            content = f.read()
        
        changed = False
        for old, new in replacements.items():
            if old in content:
                content = content.replace(old, new)
                changed = True
        
        if changed:
            with open(path, 'w') as f:
                f.write(content)
            print(f"Updated IDs in script: {filename}")
