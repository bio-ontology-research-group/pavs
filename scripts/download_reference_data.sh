#!/usr/bin/env bash
# download_reference_data.sh — Download all external reference files required by PAVS.
#
# Run from the repository root:
#   bash scripts/download_reference_data.sh
#
# Files are downloaded to data/reference/ (created if it doesn't exist).
# Large files (gnomAD, GTEx, GRCh38 GFF) are downloaded only if not already present.

set -euo pipefail

REFDIR="data/reference"
mkdir -p "$REFDIR"

log() { echo "==> $*"; }
skip_if_exists() {
    if [ -f "$1" ] || [ -d "$1" ]; then
        echo "    Already exists: $1 (skipping)"
        return 0
    fi
    return 1
}

# ── HPO Disease Annotations ────────────────────────────────────────────────────
log "HPO disease annotations (phenotype.hpoa)"
skip_if_exists "$REFDIR/phenotype.hpoa" || \
    curl -fSL "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype.hpoa" \
         -o "$REFDIR/phenotype.hpoa"

# ── HPO Gene Annotations ───────────────────────────────────────────────────────
log "HPO gene-to-disease (genes_to_disease.txt)"
skip_if_exists "$REFDIR/genes_to_disease.txt" || \
    curl -fSL "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_disease.txt" \
         -o "$REFDIR/genes_to_disease.txt"

log "HPO genes-to-phenotype (genes_to_phenotype.txt)"
skip_if_exists "$REFDIR/genes_to_phenotype.txt" || \
    curl -fSL "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt" \
         -o "$REFDIR/genes_to_phenotype.txt"

# ── OMIM Gene mapping ──────────────────────────────────────────────────────────
log "OMIM mim2gene.txt"
skip_if_exists "$REFDIR/mim2gene.txt" || \
    curl -fSL "https://omim.org/static/omim/data/mim2gene.txt" \
         -o "$REFDIR/mim2gene.txt"

# ── NCBI Gene Info ─────────────────────────────────────────────────────────────
log "NCBI gene info (Homo_sapiens.gene_info)"
skip_if_exists "$REFDIR/Homo_sapiens.gene_info" || {
    curl -fSL "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz" \
         -o "$REFDIR/Homo_sapiens.gene_info.gz"
    gunzip "$REFDIR/Homo_sapiens.gene_info.gz"
}

# ── ClinVar ────────────────────────────────────────────────────────────────────
log "ClinVar VCF (GRCh38) — large file, may take a while"
skip_if_exists "$REFDIR/clinvar.vcf.gz" || \
    curl -fSL "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" \
         -o "$REFDIR/clinvar.vcf.gz"

# ── ClinGen (erepo) ────────────────────────────────────────────────────────────
log "ClinGen expert-curated variant repository (erepo.tabbed.txt)"
skip_if_exists "$REFDIR/erepo.tabbed.txt" || \
    curl -fSL "https://erepo.clinicalgenome.org/evrepo/api/classifications/export?type=tabbed" \
         -o "$REFDIR/erepo.tabbed.txt"

# ── gnomAD Constraint Metrics ─────────────────────────────────────────────────
log "gnomAD v4.1 constraint metrics (pLI, LOEUF) — large file"
skip_if_exists "$REFDIR/gnomad.v4.1.constraint_metrics.tsv" || \
    curl -fSL "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv" \
         -o "$REFDIR/gnomad.v4.1.constraint_metrics.tsv"

# ── GO Annotations ─────────────────────────────────────────────────────────────
log "Gene Ontology human annotations (goa_human.gaf.gz)"
skip_if_exists "$REFDIR/goa_human.gaf.gz" || \
    curl -fSL "https://current.geneontology.org/annotations/goa_human.gaf.gz" \
         -o "$REFDIR/goa_human.gaf.gz"

# ── Mouse Phenotype Data (MGI/IMPC) ───────────────────────────────────────────
log "MGI HMD_HumanPhenotype.rpt (human-mouse homologs with phenotype)"
skip_if_exists "$REFDIR/HMD_HumanPhenotype.rpt" || \
    curl -fSL "https://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt" \
         -o "$REFDIR/HMD_HumanPhenotype.rpt"

log "MGI GenePheno.rpt (genotype-phenotype associations)"
skip_if_exists "$REFDIR/MGI_GenePheno.rpt" || \
    curl -fSL "https://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt" \
         -o "$REFDIR/MGI_GenePheno.rpt"

# ── GRCh38 Genomic GFF ────────────────────────────────────────────────────────
log "GRCh38 genomic GFF (GCF_000001405.40_GRCh38.p14) — very large file"
GFF="$REFDIR/GCF_000001405.40_GRCh38.p14_genomic.gff"
skip_if_exists "$GFF" || {
    curl -fSL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz" \
         -o "${GFF}.gz"
    gunzip "${GFF}.gz"
}

# ── GTEx v8 Expression ─────────────────────────────────────────────────────────
log "GTEx v8 expression data (E-GTEX-8) — download from EBI Expression Atlas"
GTEX_DIR="$REFDIR/E-GTEX-8"
skip_if_exists "$GTEX_DIR" || {
    mkdir -p "$GTEX_DIR"
    BASE="https://www.ebi.ac.uk/gxa/experiments-content/E-GTEX-8/resources/ExperimentDownloadSupplier.RnaSeqBaseline"
    curl -fSL "$BASE/tpms.tsv.gz" -o "$GTEX_DIR/E-GTEX-8-tpms.tsv.gz"
    gunzip "$GTEX_DIR/E-GTEX-8-tpms.tsv.gz"
    curl -fSL "https://www.ebi.ac.uk/gxa/experiments/E-GTEX-8/download?fileType=experiment-design" \
         -o "$GTEX_DIR/E-GTEX-8-condensed-sdrf.tsv"
    curl -fSL "https://www.ebi.ac.uk/gxa/json/experiments/E-GTEX-8" \
         -o "$GTEX_DIR/E-GTEX-8-configuration.xml" 2>/dev/null || true
}

# ── Saudi Allele Frequencies ───────────────────────────────────────────────────
log ""
log "NOTE: Saudi allele frequencies (data/reference/allele-frequencies/) are"
log "      not publicly available. Contact the authors for access."
log "      These are derived from 302 WGS individuals."

log ""
log "All public reference files downloaded to $REFDIR/"
log "Next steps:"
log "  1. Obtain data/reference/allele-frequencies/ from the authors"
log "  2. Run: uv run python normalization/combine_normalize_phenotypes.py ..."
log "  3. Run: uv run python normalization/annotate_variants.py ..."
