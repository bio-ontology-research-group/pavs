import os
import subprocess

files = [
    ('data/cleaned_and_normalized/ahmed_normalized.tsv', 'phenopackets/generated/ahmed'),
    ('data/cleaned_and_normalized/fawzan_normalized.tsv', 'phenopackets/generated/fawzan'),
    ('data/cleaned_and_normalized/marwa_normalized.tsv', 'phenopackets/generated/marwa'),
    ('data/cleaned_and_normalized/ahmed_pmid28454995_normalized.tsv', 'phenopackets/generated/pmid28454995'),
    ('data/cleaned_and_normalized/pmc7082194_normalized.tsv', 'phenopackets/generated/pmid31618753')
]

for input_file, output_dir in files:
    if os.path.exists(input_file):
        print(f"Generating phenopackets for {input_file} into {output_dir}...")
        # Ensure output_dir is fresh
        if os.path.exists(output_dir):
            import shutil
            shutil.rmtree(output_dir)
        
        subprocess.run([
            "uv", "run", "python", "scripts/generate_phenopackets_v2.py",
            "--input", input_file,
            "--output_dir", output_dir
        ])
    else:
        print(f"Skipping {input_file} (not found)")