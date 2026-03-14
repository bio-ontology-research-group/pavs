#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# PAVS Setup Script
# Initializes all submodules and prepares the development environment
# ─────────────────────────────────────────────────────────────────────────────

set -e

echo "=== PAVS Setup Script ==="
echo

# Check if git submodules are initialized
echo "1. Checking submodules..."
if [ ! -f "knowledge-graph/.git" ]; then
    echo "   Initializing submodules..."
    git submodule update --init --recursive
else
    echo "   ✓ Submodules already initialized"
fi
echo

# Check for required tools
echo "2. Checking prerequisites..."
command -v docker >/dev/null 2>&1 || { echo "   ✗ Docker not found. Please install Docker."; exit 1; }
echo "   ✓ Docker found"

command -v uv >/dev/null 2>&1 || { echo "   ✗ uv not found. Install with: curl -LsSf https://astral.sh/uv/install.sh | sh"; exit 1; }
echo "   ✓ uv found"

if [ ! -d "knowledge-graph/ontology" ]; then
    mkdir -p knowledge-graph/ontology
fi
echo

# Check for ontology files
echo "3. Checking ontology files..."
if [ ! -f "knowledge-graph/ontology/hp.obo" ]; then
    echo "   Downloading HPO ontology..."
    wget -q https://purl.obolibrary.org/obo/hp.obo -P knowledge-graph/ontology/
    echo "   ✓ Downloaded hp.obo"
else
    echo "   ✓ hp.obo exists"
fi

if [ ! -f "knowledge-graph/ontology/mondo.obo" ]; then
    echo "   Downloading MONDO ontology..."
    wget -q http://purl.obolibrary.org/obo/mondo.obo -P knowledge-graph/ontology/
    echo "   ✓ Downloaded mondo.obo"
else
    echo "   ✓ mondo.obo exists"
fi
echo

# Check for API key
echo "4. Checking environment..."
if [ -z "$OPENROUTER_API_KEY" ]; then
    echo "   ⚠ OPENROUTER_API_KEY not set"
    echo "     Get a free key at https://openrouter.ai/"
    echo "     Then: export OPENROUTER_API_KEY='your_key'"
else
    echo "   ✓ OPENROUTER_API_KEY is set"
fi
echo

# Install phenotype-matcher
echo "5. Installing phenotype-matcher library..."
if [ ! -d "tools/phenotype_matcher_v2/.git" ]; then
    echo "   ✗ Submodule not initialized. Run: git submodule update --init"
    exit 1
fi

cd tools/phenotype_matcher_v2
if pip show phenotype-matcher-v2 >/dev/null 2>&1; then
    echo "   ✓ phenotype-matcher-v2 already installed"
else
    echo "   Installing in development mode..."
    pip install -e . >/dev/null 2>&1
    echo "   ✓ Installed phenotype-matcher-v2"
fi
cd ../..
echo

# Install knowledge-graph dependencies
echo "6. Installing knowledge-graph dependencies..."
cd knowledge-graph
uv sync --quiet
echo "   ✓ Dependencies installed"
cd ..
echo

# Create directories
echo "7. Creating output directories..."
mkdir -p rdf_output
mkdir -p knowledge-graph/rdf_output
echo "   ✓ Directories created"
echo

echo "=== Setup Complete ==="
echo
echo "Next steps:"
echo "  1. Download reference data:"
echo "       cd knowledge-graph && bash scripts/download_reference_data.sh"
echo
echo "  2. Generate RDF data:"
echo "       cd knowledge-graph"
echo "       uv run python normalization/combine_normalize_phenotypes.py --output data/combined_normalized.tsv"
echo "       uv run python normalization/annotate_variants.py --input data/combined_normalized.tsv --output data/combined_annotated.tsv"
echo "       uv run python intake/generate_phenopackets_v2.py --input data/combined_annotated.tsv --combined data/PAVS_phenopackets.json"
echo "       uv run python intake/compute_hpo_ic.py --hpo ontology/hp.obo --hpoa data/reference/phenotype.hpoa --output rdf_output/hpo_ic.ttl"
echo "       uv run python intake/generate_rdf.py --phenopackets data/PAVS_phenopackets.json --annotated data/combined_annotated.tsv --output-dir ../rdf_output/"
echo "       cd .."
echo
echo "  3. Start the website:"
echo "       docker compose -f docker-compose-sparql.yml up -d"
echo "       # Frontend: http://localhost:3000"
echo "       # Backend API: http://localhost:8000/docs"
echo "       # SPARQL: http://localhost:8890/sparql"
echo
