#!/bin/bash
# Example batch and parallel processing with the Phenotype Matcher

echo "======================================================================"
echo "Phenotype Matcher - Batch and Parallel Processing Examples"
echo "======================================================================"
echo ""

# Create a sample input file
cat > /tmp/sample_phenotypes.txt << 'EOF'
# Sample phenotype descriptions (one per line)
severe intellectual disability
seizures and hypotonia
feeding difficulties
no cardiac abnormalities
short femur and absent radius
polyhydramnios
cough and noisy breathing
ventricular arrhythmia
global developmental delay
profound hearing loss
EOF

echo "Created sample input file: /tmp/sample_phenotypes.txt"
echo ""
echo "Contents:"
cat /tmp/sample_phenotypes.txt
echo ""
echo "======================================================================"

# Example 1: Sequential processing
echo ""
echo "Example 1: Sequential Processing (1 thread)"
echo "----------------------------------------------------------------------"
echo "Command: phenotype-matcher --input /tmp/sample_phenotypes.txt --output-format tsv -j 1"
echo ""
# phenotype-matcher --input /tmp/sample_phenotypes.txt --output-format tsv -j 1

# Example 2: Parallel processing with 4 threads
echo ""
echo "Example 2: Parallel Processing (4 threads)"
echo "----------------------------------------------------------------------"
echo "Command: phenotype-matcher --input /tmp/sample_phenotypes.txt --output /tmp/results.jsonl -j 4 --progress"
echo ""
# phenotype-matcher --input /tmp/sample_phenotypes.txt --output /tmp/results.jsonl -j 4 --progress

# Example 3: Batch mode from command line
echo ""
echo "Example 3: Batch Mode (command line arguments)"
echo "----------------------------------------------------------------------"
echo "Command: phenotype-matcher -b 'seizures' -b 'hypotonia' -b 'cough' --output-format hpo_ids"
echo ""
# phenotype-matcher -b "seizures" -b "hypotonia" -b "cough" --output-format hpo_ids

# Example 4: Parallel TSV output
echo ""
echo "Example 4: Parallel TSV Output"
echo "----------------------------------------------------------------------"
echo "Command: phenotype-matcher --input /tmp/sample_phenotypes.txt --output /tmp/results.tsv --output-format tsv -j 4"
echo ""
# phenotype-matcher --input /tmp/sample_phenotypes.txt --output /tmp/results.tsv --output-format tsv -j 4

# Example 5: With UV
echo ""
echo "Example 5: Using UV (recommended)"
echo "----------------------------------------------------------------------"
echo "Command: uv run phenotype-matcher --input /tmp/sample_phenotypes.txt -j 4 --progress"
echo ""
# uv run phenotype-matcher --input /tmp/sample_phenotypes.txt -j 4 --progress

echo ""
echo "======================================================================"
echo "To run these examples, uncomment the command lines in this script"
echo "Make sure you have set OPENROUTER_API_KEY environment variable"
echo "======================================================================"
