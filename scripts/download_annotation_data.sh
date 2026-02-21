#!/usr/bin/env bash
# download_annotation_data.sh — one-time downloads for variant annotation enrichment
# Run from the repository root: bash scripts/download_annotation_data.sh
# All files that are already present are skipped.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${REPO_ROOT}/data"
ONT_DIR="${REPO_ROOT}/ontology"

mkdir -p "${DATA_DIR}" "${ONT_DIR}" "${DATA_DIR}/annotation_cache"

# Helper: skip if file already exists
need() {
    local path="$1"
    [[ ! -f "${path}" ]]
}

echo "=== PAVS annotation data downloader ==="
echo "Repo root: ${REPO_ROOT}"
echo ""

# ----------------------------------------------------------------
# 1. MGI human–mouse ortholog file (~1 MB)
#    Maps human HGNC symbol → mouse MGI gene ID.
#    Required to bridge genes_to_disease → MGI_GenePheno phenotypes.
# ----------------------------------------------------------------
HMD="${DATA_DIR}/HMD_HumanPhenotype.rpt"
if need "${HMD}"; then
    echo "[1/4] Downloading MGI HMD_HumanPhenotype.rpt ..."
    wget -q -O "${HMD}" \
        "https://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt"
    echo "      → ${HMD}"
else
    echo "[1/4] HMD_HumanPhenotype.rpt already present — skipping."
fi

# ----------------------------------------------------------------
# 2. GO ontology (go-basic.obo, ~40 MB)
#    Required for GO ID → term label mapping.
# ----------------------------------------------------------------
GO_OBO="${ONT_DIR}/go-basic.obo"
if need "${GO_OBO}"; then
    echo "[2/4] Downloading go-basic.obo ..."
    wget -q -O "${GO_OBO}" \
        "http://purl.obolibrary.org/obo/go/go-basic.obo"
    echo "      → ${GO_OBO}"
else
    echo "[2/4] go-basic.obo already present — skipping."
fi

# ----------------------------------------------------------------
# 3. gnomAD v4.1 constraint scores (~50 MB)
#    Provides pLI and LOEUF per gene.
# ----------------------------------------------------------------
GNOMAD="${DATA_DIR}/gnomad.v4.1.constraint_metrics.tsv"
if need "${GNOMAD}"; then
    echo "[3/4] Downloading gnomAD v4.1 constraint metrics ..."
    wget -q -P "${DATA_DIR}" \
        "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
    echo "      → ${GNOMAD}"
else
    echo "[3/4] gnomAD constraint already present — skipping."
fi

# ----------------------------------------------------------------
# 4. Tabix-index the ClinVar VCF for fast positional lookup (optional)
#    Requires tabix to be installed (from htslib).
#    The annotation script works without this index (slower scan mode),
#    so this step is advisory only.
# ----------------------------------------------------------------
CLINVAR_TBI="${DATA_DIR}/clinvar.vcf.gz.tbi"
if need "${CLINVAR_TBI}"; then
    if command -v tabix &>/dev/null; then
        echo "[4/4] Creating tabix index for clinvar.vcf.gz ..."
        tabix -p vcf "${DATA_DIR}/clinvar.vcf.gz"
        echo "      → ${CLINVAR_TBI}"
    else
        echo "[4/4] tabix not found — ClinVar will be scanned without index (slower)."
        echo "      Install htslib and re-run to enable fast positional lookup."
    fi
else
    echo "[4/4] clinvar.vcf.gz.tbi already present — skipping."
fi

# ----------------------------------------------------------------
# Optional VEP plugin files (commented out — very large)
# ----------------------------------------------------------------
cat << 'EOF'

=== Optional large VEP plugin files (NOT downloaded automatically) ===

# REVEL missense score (~2 GB):
#   wget -P data/vep_plugins/ https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip

# AlphaMissense GRCh38 (~8 GB):
#   wget -P data/vep_plugins/ https://zenodo.org/record/8208688/files/AlphaMissense_hg38.tsv.gz

# SpliceAI GRCh38 SNVs (~30 GB):
#   wget -P data/vep_plugins/ https://basespace.illumina.com/...

# CADD GRCh38 (~80 GB):
#   wget -P data/vep_plugins/ https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz

# VEP GRCh38 offline cache (~15 GB):
#   docker pull ensemblorg/ensembl-vep
#   docker run --rm -v $(pwd)/data/vep_cache:/opt/vep/.vep ensemblorg/ensembl-vep \
#     INSTALL.pl -a cf -s homo_sapiens -y GRCh38 --CACHEDIR /opt/vep/.vep

EOF

echo ""
echo "=== Done. Run 'uv run python scripts/annotate_variants.py --help' to start annotation. ==="
