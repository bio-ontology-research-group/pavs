# ─────────────────────────────────────────────────────────────────────────────
# Dockerfile for PAVS Knowledge Graph Pipeline
# Unified image for normalization, annotation, phenopacket generation, and RDF
# ─────────────────────────────────────────────────────────────────────────────
FROM python:3.12-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install uv for faster dependency management
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.cargo/bin:${PATH}"

# Copy dependency files
COPY pyproject.toml uv.lock requirements.txt ./

# Install Python dependencies
RUN uv sync --frozen || pip install -r requirements.txt

# Install phenotype-matcher from GitHub
RUN pip install git+https://github.com/bio-ontology-research-group/phenotype-matcher.git

# Copy application code
COPY normalization/ ./normalization/
COPY intake/ ./intake/
COPY scripts/ ./scripts/

# Copy source data (small TSV files only)
COPY data/phenotypes/ ./data/phenotypes/

# Create directories for outputs
RUN mkdir -p data/reference ontology rdf_output phenopackets/generated_v2

# Default command shows available scripts
CMD ["echo", "Available commands: combine_normalize_phenotypes.py, annotate_variants.py, generate_phenopackets_v2.py, generate_rdf.py"]
