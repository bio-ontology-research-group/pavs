#!/usr/bin/env python3
"""
generate_phenopackets_v2.py — Generate Phenopackets v2.0.2 from combined_annotated.tsv.

Key features:
- Reads data/combined_annotated.tsv (7,510 rows × 77+ columns)
- Emits full vcfRecord (chrom/pos/ref/alt from vep_hgvsg)
- Sets acmgPathogenicityClassification from acmg_classification column
- Adds HANCESTRO:0852 ("Middle Eastern") to all Saudi cases in phenotypicFeatures
- Adds ClinVar allele ID to variantInterpretation where available
- Stores source badge in metaData.externalReferences

Usage:
    uv run python scripts/generate_phenopackets_v2.py \
        --input data/combined_annotated.tsv \
        --output-dir phenopackets/generated_v2 \
        --combined data/PAVS_phenopackets.json \
        --zip data/PAVS_phenopackets.zip
"""

import argparse
import ast
import datetime
import json
import os
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SAUDI_SOURCES = {
    "ahmed-variants", "ahmed-pmid28454995", "fawzan-variants",
    "marwa-variants", "PMC6562004", "PMC7082194",
}

ACMG_NORM = {
    "Pathogenic": "PATHOGENIC",
    "Likely pathogenic": "LIKELY_PATHOGENIC",
    "Uncertain significance": "UNCERTAIN_SIGNIFICANCE",
    "Uncertain Significance": "UNCERTAIN_SIGNIFICANCE",
    "VUS": "UNCERTAIN_SIGNIFICANCE",
    "Likely benign": "LIKELY_BENIGN",
    "Benign": "BENIGN",
    "pathogenic": "PATHOGENIC",
    "likely pathogenic": "LIKELY_PATHOGENIC",
}

ZYGOSITY_LABEL = {
    "GENO:0000136": "homozygous",
    "GENO:0000135": "heterozygous",
    "GENO:0000134": "hemizygous",
}

META_DATE = datetime.date.today().isoformat()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def safe(val) -> bool:
    return val is not None and str(val).strip() not in ("", "nan", "None", "NaN")


def parse_pipe(val) -> list:
    if not safe(val):
        return []
    return [x.strip() for x in str(val).split("|") if x.strip()]


def parse_hgvsg(hgvsg: str):
    """Extract chrom, pos, ref, alt from HGVS g. notation NC_000003.12:g.4417183A>G"""
    if not safe(hgvsg):
        return None
    m = re.match(r"NC_(\d+)\.\d+:g\.(\d+)([A-Z]+)>([A-Z]+)", str(hgvsg))
    if not m:
        return None
    chrom_num = int(m.group(1))
    pos = int(m.group(2))
    ref = m.group(3)
    alt = m.group(4)
    if chrom_num == 23:
        chrom = "chrX"
    elif chrom_num == 24:
        chrom = "chrY"
    else:
        chrom = f"chr{chrom_num}"
    return {"genomeAssembly": "hg38", "chrom": chrom, "pos": pos, "ref": ref, "alt": alt}


def acmg_code(val: str) -> str:
    return ACMG_NORM.get(str(val).strip(), "UNCERTAIN_SIGNIFICANCE")


def build_metadata(source: str) -> dict:
    resources = [
        {
            "id": "hp",
            "name": "Human Phenotype Ontology",
            "url": "https://hpo.jax.org",
            "version": "2024-08-13",
            "namespacePrefix": "HP",
            "iriPrefix": "http://purl.obolibrary.org/obo/HP_",
        },
        {
            "id": "geno",
            "name": "Genotype Ontology",
            "url": "https://github.com/oborel/obo-relations",
            "namespacePrefix": "GENO",
            "iriPrefix": "http://purl.obolibrary.org/obo/GENO_",
        },
        {
            "id": "omim",
            "name": "Online Mendelian Inheritance in Man",
            "url": "https://www.omim.org",
            "namespacePrefix": "OMIM",
            "iriPrefix": "https://www.omim.org/entry/",
        },
        {
            "id": "mondo",
            "name": "Monarch Disease Ontology",
            "url": "http://purl.obolibrary.org/obo/mondo.obo",
            "namespacePrefix": "MONDO",
            "iriPrefix": "http://purl.obolibrary.org/obo/MONDO_",
        },
        {
            "id": "hancestro",
            "name": "Human Ancestry Ontology",
            "url": "https://github.com/EBISPOT/ancestro",
            "namespacePrefix": "HANCESTRO",
            "iriPrefix": "http://purl.obolibrary.org/obo/HANCESTRO_",
        },
        {
            "id": "gaz",
            "name": "Gazetteer",
            "url": "http://purl.obolibrary.org/obo/gaz.obo",
            "namespacePrefix": "GAZ",
            "iriPrefix": "http://purl.obolibrary.org/obo/GAZ_",
        },
    ]
    return {
        "created": META_DATE,
        "createdBy": "PAVS pipeline",
        "submittedBy": "Robert Hoehndorf",
        "resources": resources,
        "phenopacketSchemaVersion": "2.0",
        "externalReferences": [
            {"id": f"pavs:source:{source}", "description": f"Source dataset: {source}"}
        ],
    }


def build_phenopacket(row: pd.Series) -> dict:
    pavs_id = str(row["pavs_id"]).strip()
    source = str(row.get("source_file", "unknown")).strip()
    is_saudi = source in SAUDI_SOURCES

    # Subject
    sex = str(row.get("sex_id", "")).strip()
    sex_map = {
        "NCIT:C20197": "MALE",
        "NCIT:C16576": "FEMALE",
    }
    sex_label = sex_map.get(sex, "UNKNOWN_SEX")

    subject: dict = {
        "id": pavs_id,
        "sex": sex_label,
    }
    if safe(row.get("age_iso")):
        subject["timeAtLastEncounter"] = {
            "age": {"iso8601duration": str(row["age_iso"]).strip()}
        }

    # Phenotypic features
    phenotypic_features = []

    # Present HPO terms
    present_ids = parse_pipe(row.get("phenotypes_present_ids", ""))
    present_labels = parse_pipe(row.get("phenotypes_present_labels", ""))
    for i, hpo_id in enumerate(present_ids):
        if not hpo_id.startswith("HP:"):
            continue
        label = present_labels[i] if i < len(present_labels) else hpo_id
        feat: dict = {
            "type": {"id": hpo_id, "label": label},
            "excluded": False,
        }
        # Modifier
        modifier_ids = parse_pipe(row.get("phenotypes_modifier_ids", ""))
        modifier_labels = parse_pipe(row.get("phenotypes_modifier_labels", ""))
        if modifier_ids:
            feat["modifiers"] = [
                {"id": mid, "label": (modifier_labels[j] if j < len(modifier_labels) else mid)}
                for j, mid in enumerate(modifier_ids)
            ]
        phenotypic_features.append(feat)

    # Excluded HPO terms
    excluded_ids = parse_pipe(row.get("phenotypes_excluded_ids", ""))
    excluded_labels = parse_pipe(row.get("phenotypes_excluded_labels", ""))
    for i, hpo_id in enumerate(excluded_ids):
        if not hpo_id.startswith("HP:"):
            continue
        label = excluded_labels[i] if i < len(excluded_labels) else hpo_id
        phenotypic_features.append({
            "type": {"id": hpo_id, "label": label},
            "excluded": True,
        })

    # Ancestry: HANCESTRO:0852 for Saudi cases
    if is_saudi:
        phenotypic_features.append({
            "type": {"id": "HANCESTRO:0852", "label": "Middle Eastern"},
            "excluded": False,
        })

    # Diseases
    diseases = []
    omim_ids = parse_pipe(row.get("disease_omim_ids", ""))
    omim_labels = parse_pipe(row.get("disease_omim_labels", ""))
    mondo_id = str(row.get("disease_mondo_id", "")).strip()
    mondo_label = str(row.get("disease_mondo_label", "")).strip()

    if safe(mondo_id):
        diseases.append({
            "term": {"id": mondo_id, "label": mondo_label if safe(mondo_label) else mondo_id}
        })
    elif omim_ids:
        for i, oid in enumerate(omim_ids):
            if safe(oid):
                ol = omim_labels[i] if i < len(omim_labels) else oid
                diseases.append({"term": {"id": oid, "label": ol if safe(ol) else oid}})

    # Variant interpretation
    gene = str(row.get("gene_symbol", "")).strip()
    hgvsc = str(row.get("vep_hgvsc", row.get("variant_cdna", ""))).strip()
    hgvsp = str(row.get("vep_hgvsp", row.get("variant_protein", ""))).strip()
    hgvsg = str(row.get("vep_hgvsg", row.get("variant_genomic_grch38", ""))).strip()
    acmg_raw = str(row.get("acmg_classification", "")).strip()
    acmg_label = acmg_code(acmg_raw) if safe(acmg_raw) else "UNCERTAIN_SIGNIFICANCE"
    result_status = str(row.get("result_status", "")).strip().upper()
    progress = "SOLVED" if result_status == "SOLVED" else "IN_PROGRESS"

    zyg_id = str(row.get("zygosity_id", "")).strip()
    zyg_label = ZYGOSITY_LABEL.get(zyg_id, "heterozygous")

    gene_entrez = str(row.get("gene_entrez_id", "")).strip()

    # Build variation descriptor
    expressions = []
    if safe(hgvsc):
        expressions.append({"syntax": "hgvs.c", "value": hgvsc})
    if safe(hgvsg):
        expressions.append({"syntax": "hgvs.g", "value": hgvsg})
    if safe(hgvsp):
        expressions.append({"syntax": "hgvs.p", "value": hgvsp})

    vcf_rec = parse_hgvsg(hgvsg)

    var_desc: dict = {
        "id": f"{pavs_id}_{gene}",
        "moleculeContext": "genomic",
    }
    if safe(gene):
        gene_ctx: dict = {"symbol": gene}
        if safe(gene_entrez):
            gene_ctx["valueId"] = f"NCBIGene:{gene_entrez}"
        var_desc["geneContext"] = gene_ctx
    if expressions:
        var_desc["expressions"] = expressions
    if vcf_rec:
        var_desc["vcfRecord"] = vcf_rec

    # Zygosity as allelic state
    if zyg_id:
        var_desc["allelicState"] = {
            "id": zyg_id,
            "label": zyg_label,
        }

    # ClinVar
    clinvar_id = str(row.get("clinvar_allele_id", "")).strip()
    clinvar_sig = str(row.get("clinvar_sig", "")).strip()
    if safe(clinvar_id):
        var_desc["vcfRecord"] = var_desc.get("vcfRecord") or {}
        # Add ClinVar reference to the descriptor
        var_desc["xrefs"] = [f"ClinVar:{clinvar_id}"]
        if safe(clinvar_sig):
            var_desc["clinvarId"] = clinvar_id

    interp_status = "CAUSATIVE" if acmg_label in ("PATHOGENIC", "LIKELY_PATHOGENIC") else "CONTRIBUTORY"

    variant_interpretation: dict = {
        "acmgPathogenicityClassification": acmg_label,
        "therapeuticActionability": "UNKNOWN_ACTIONABILITY",
        "variationDescriptor": var_desc,
    }

    genomic_interpretation: dict = {
        "subjectOrBiosampleId": pavs_id,
        "interpretationStatus": interp_status,
        "variantInterpretation": variant_interpretation,
    }

    # Inheritance
    inheritance_id = str(row.get("inheritance_id", "")).strip()
    inheritance_label = str(row.get("inheritance_label", "")).strip()
    disease_with_inheritance = diseases.copy()
    if disease_with_inheritance and safe(inheritance_id):
        disease_with_inheritance[0]["modeOfInheritance"] = {
            "id": inheritance_id,
            "label": inheritance_label,
        }

    # Build interpretation
    interp_obj = {
        "id": f"interp-{pavs_id}",
        "progressStatus": progress,
    }
    if disease_with_inheritance or safe(gene):
        diag: dict = {}
        if disease_with_inheritance:
            diag["disease"] = disease_with_inheritance[0]["term"]
        diag["genomicInterpretations"] = [genomic_interpretation]
        interp_obj["diagnosis"] = diag

    # Individual extensions for consanguinity
    extensions = []
    consanguinity_id = str(row.get("consanguinity_id", "")).strip()
    consanguinity_label = str(row.get("consanguinity_label", "")).strip()
    family_history_id = str(row.get("family_history_id", "")).strip()
    family_history_label = str(row.get("family_history_label", "")).strip()

    if safe(consanguinity_id):
        extensions.append({
            "name": "consanguinity",
            "value": {"ontologyClass": {"id": consanguinity_id, "label": consanguinity_label}},
        })
    if safe(family_history_id):
        extensions.append({
            "name": "familyHistory",
            "value": {"ontologyClass": {"id": family_history_id, "label": family_history_label}},
        })
    if is_saudi:
        extensions.append({"name": "geographicOrigin", "value": {"id": "GAZ:00005279", "label": "Saudi Arabia"}})

    if extensions:
        subject["extensions"] = extensions

    # Assemble phenopacket
    pp: dict = {
        "id": pavs_id,
        "subject": subject,
        "phenotypicFeatures": phenotypic_features,
        "metaData": build_metadata(source),
    }

    if disease_with_inheritance:
        pp["diseases"] = disease_with_inheritance

    if safe(gene):
        pp["interpretations"] = [interp_obj]

    # PMID
    pmid = str(row.get("pmid", "")).strip()
    if safe(pmid):
        pp["metaData"]["externalReferences"].append({
            "id": f"PMID:{pmid}",
            "reference": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        })

    return pp


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate Phenopackets v2.0.2 from combined_annotated.tsv")
    parser.add_argument("--input", default="data/combined_annotated.tsv")
    parser.add_argument("--output-dir", default="phenopackets/generated_v2")
    parser.add_argument("--combined", default="data/PAVS_phenopackets.json")
    parser.add_argument("--zip", default="data/PAVS_phenopackets.zip")
    parser.add_argument("--limit", type=int, default=None, help="Process only first N rows (for testing)")
    parser.add_argument("--saudi-only", action="store_true", help="Only write Saudi phenopackets to combined JSON")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.input} …")
    df = pd.read_csv(args.input, sep="\t", dtype=str, low_memory=False)
    print(f"  {len(df)} rows")

    if args.limit:
        df = df.head(args.limit)
        print(f"  Limited to {len(df)} rows")

    all_packets = []
    saudi_packets = []
    errors = 0

    for idx, row in df.iterrows():
        try:
            pp = build_phenopacket(row)
        except Exception as e:
            print(f"  ERROR row {idx} ({row.get('pavs_id', '?')}): {e}")
            errors += 1
            continue

        pid = pp["id"]
        # Write individual file
        out_path = out_dir / f"{pid}.json"
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump(pp, fh, indent=2, ensure_ascii=False)

        all_packets.append(pp)
        source = str(row.get("source_file", "")).strip()
        if source in SAUDI_SOURCES:
            saudi_packets.append(pp)

        if (idx + 1) % 500 == 0:
            print(f"  … {idx + 1} rows processed")

    print(f"\nGenerated {len(all_packets)} phenopackets ({len(saudi_packets)} Saudi), {errors} errors")
    print(f"Individual files → {out_dir}/")

    # Write combined JSON
    combined_path = Path(args.combined)
    combined_path.parent.mkdir(parents=True, exist_ok=True)
    target = saudi_packets if args.saudi_only else all_packets
    with open(combined_path, "w", encoding="utf-8") as fh:
        json.dump(target, fh, indent=2, ensure_ascii=False)
    print(f"Combined JSON ({len(target)} phenopackets) → {combined_path}")

    # Write ZIP
    zip_path = Path(args.zip)
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for pp in target:
            pid = pp["id"]
            data = json.dumps(pp, indent=2, ensure_ascii=False)
            zf.writestr(f"{pid}.json", data)
    print(f"ZIP ({len(target)} files) → {zip_path}")


if __name__ == "__main__":
    main()
