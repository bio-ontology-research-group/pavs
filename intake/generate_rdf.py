#!/usr/bin/env python3
"""
generate_rdf.py — Convert all PAVS data sources to Turtle (.ttl) files
ready for bulk-loading into Virtuoso.

Outputs:
  rdf_output/cases.ttl       — Saudi case records (from PAVS_phenopackets.json)
  rdf_output/genes.ttl       — Gene annotations (from combined_annotated.tsv)
  rdf_output/hpoa.ttl        — Disease→HPO from phenotype.hpoa
  rdf_output/ddd.ttl         — DDD non-Saudi variants
  rdf_output/literature.ttl  — 9,588 literature phenopackets

Usage:
    uv run python scripts/generate_rdf.py \
        --phenopackets data/PAVS_phenopackets.json \
        --annotated data/combined_annotated.tsv \
        --literature-dir phenopackets/0.1.26 \
        --output-dir rdf_output/
"""

import argparse
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Namespace prefixes (shared across all output files)
# ---------------------------------------------------------------------------

PREFIXES = """\
@prefix pavs:    <http://pavs.kaust.edu.sa/ontology/> .
@prefix pav:     <http://pavs.kaust.edu.sa/data/> .
@prefix hp:      <http://purl.obolibrary.org/obo/HP_> .
@prefix mondo:   <http://purl.obolibrary.org/obo/MONDO_> .
@prefix geno:    <http://purl.obolibrary.org/obo/GENO_> .
@prefix hancestro: <http://purl.obolibrary.org/obo/HANCESTRO_> .
@prefix gaz:     <http://purl.obolibrary.org/obo/GAZ_> .
@prefix omim:    <https://omim.org/entry/> .
@prefix ncbigene: <http://identifiers.org/ncbigene/> .
@prefix dbsnp:   <http://identifiers.org/dbsnp/> .
@prefix clinvar: <http://identifiers.org/clinvar/> .
@prefix hgnc:    <http://identifiers.org/hgnc.symbol/> .
@prefix dc:      <http://purl.org/dc/terms/> .
@prefix rdfs:    <http://www.w3.org/2000/01/rdf-schema#> .
@prefix xsd:     <http://www.w3.org/2001/XMLSchema#> .
@prefix owl:     <http://www.w3.org/2002/07/owl#> .

"""

# Saudi source files (non-DDD, non-Literature)
SAUDI_SOURCES = {
    "ahmed-variants", "ahmed-pmid28454995", "fawzan-variants",
    "marwa-variants", "PMC6562004", "PMC7082194",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def ttl_str(s: str) -> str:
    """Escape a string for Turtle."""
    s = str(s).replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n").replace("\r", "")
    return f'"{s}"'


def ttl_uri(uri: str) -> str:
    """Wrap a full URI in angle brackets."""
    return f"<{uri}>"


def safe(val) -> bool:
    """Return True if value is non-null and non-empty."""
    return val is not None and str(val).strip() not in ("", "nan", "None")


def hp_uri(hpo_id: str) -> str:
    """HP:0001263 → hp:0001263"""
    return "hp:" + hpo_id.replace("HP:", "")


def omim_uri(omim_id: str) -> str:
    """OMIM:272200 → omim:272200"""
    return "omim:" + omim_id.replace("OMIM:", "")


def mondo_uri(mondo_id: str) -> str:
    """MONDO:0015286 → mondo:0015286"""
    return "mondo:" + mondo_id.replace("MONDO:", "")


def geno_uri(geno_id: str) -> str:
    """GENO:0000136 → geno:0000136"""
    return "geno:" + geno_id.replace("GENO:", "")


def case_uri(pavs_id: str) -> str:
    """PAVS_MASTER_0001 → pav:PAVS_MASTER_0001"""
    safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", pavs_id)
    return f"pav:{safe_id}"


def var_uri(pavs_id: str, gene: str) -> str:
    safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", pavs_id)
    safe_gene = re.sub(r"[^A-Za-z0-9_.-]", "_", str(gene))
    return f"pav:var_{safe_id}_{safe_gene}"


# ---------------------------------------------------------------------------
# Generate cases.ttl from PAVS_phenopackets.json
# ---------------------------------------------------------------------------

def generate_cases_ttl(phenopackets_json: str, output_path: str):
    print(f"Loading phenopackets from {phenopackets_json} …")
    with open(phenopackets_json, encoding="utf-8") as fh:
        packets = json.load(fh)

    if isinstance(packets, dict):
        # Might be wrapped
        packets = packets.get("phenopackets", list(packets.values()))

    print(f"  {len(packets)} phenopackets found")

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(PREFIXES)

        for pp in packets:
            pid = pp.get("id", "")
            if not pid:
                continue

            cu = case_uri(pid)

            # Source
            source = "unknown"
            meta = pp.get("metaData", {})
            for ext in meta.get("externalReferences", []):
                if ext.get("id", "").startswith("pavs:source:"):
                    source = ext["id"].replace("pavs:source:", "")
            if source == "unknown":
                # Try to infer from id prefix
                if pid.startswith("PAVS:A"):
                    source = "ahmed-variants"
                elif pid.startswith("PAVS:B"):
                    source = "ahmed-pmid28454995"
                elif pid.startswith("PAVS:F"):
                    source = "fawzan-variants"
                elif pid.startswith("PAVS:M"):
                    source = "marwa-variants"
                elif pid.startswith("PAVS:P"):
                    source = "PMC6562004"
                elif pid.startswith("PAVS:Q"):
                    source = "PMC7082194"
                elif pid.startswith("PAVS:D"):
                    source = "ddd-diagnoses"

            is_saudi = source in SAUDI_SOURCES or (
                "PAVS_MASTER" in pid and source != "ddd-diagnoses"
            )
            is_saudi_str = "true" if is_saudi else "false"

            fh.write(f"\n{cu}\n")
            fh.write(f"  a pavs:Case ;\n")
            fh.write(f"  dc:identifier {ttl_str(pid)} ;\n")
            fh.write(f"  pavs:source {ttl_str(source)} ;\n")
            fh.write(f"  pavs:isSaudi {is_saudi_str} ;\n")

            # Ancestry for Saudi cases
            if is_saudi:
                fh.write(f"  pavs:ancestry hancestro:0852 ;\n")
                fh.write(f"  pavs:geographicOrigin gaz:00005279 ;\n")

            # Subject sex
            subject = pp.get("subject", {})
            sex = subject.get("sex", "UNKNOWN_SEX")
            sex_map = {
                "MALE": "http://purl.obolibrary.org/obo/NCIT_C20197",
                "FEMALE": "http://purl.obolibrary.org/obo/NCIT_C16576",
            }
            if sex in sex_map:
                fh.write(f"  pavs:sex {ttl_uri(sex_map[sex])} ;\n")

            # Phenotypic features (HPO)
            for feat in pp.get("phenotypicFeatures", []):
                ftype = feat.get("type", {})
                fid = ftype.get("id", "")
                excluded = feat.get("excluded", False)
                if fid.startswith("HP:"):
                    pred = "pavs:hasExcludedPhenotype" if excluded else "pavs:hasPhenotype"
                    fh.write(f"  {pred} {hp_uri(fid)} ;\n")

            # Diseases
            for disease in pp.get("diseases", []):
                did = disease.get("term", {}).get("id", "")
                dlabel = disease.get("term", {}).get("label", "")
                if did.startswith("OMIM:"):
                    fh.write(f"  pavs:hasDisease {omim_uri(did)} ;\n")
                elif did.startswith("MONDO:"):
                    fh.write(f"  pavs:hasDisease {mondo_uri(did)} ;\n")
                if dlabel:
                    fh.write(f"  pavs:diseaseLabel {ttl_str(dlabel)} ;\n")

            # Interpretations → variants
            for interp in pp.get("interpretations", []):
                progress = interp.get("progressStatus", "")
                if progress:
                    fh.write(f"  pavs:progressStatus {ttl_str(progress)} ;\n")
                for gi in interp.get("diagnosis", {}).get("genomicInterpretations", []):
                    vi = gi.get("variantInterpretation", {})
                    vd = vi.get("variationDescriptor", {})
                    gene = vd.get("geneContext", {}).get("symbol", "unknown")
                    vu = var_uri(pid, gene)
                    fh.write(f"  pavs:hasVariant {vu} ;\n")

            fh.write(f"  rdfs:label {ttl_str(pid)} .\n")

            # Now write the variant triples
            for interp in pp.get("interpretations", []):
                for gi in interp.get("diagnosis", {}).get("genomicInterpretations", []):
                    vi = gi.get("variantInterpretation", {})
                    vd = vi.get("variationDescriptor", {})
                    gene_ctx = vd.get("geneContext", {})
                    gene_sym = gene_ctx.get("symbol", "unknown")
                    gene_val_id = gene_ctx.get("valueId", "")

                    acmg = vi.get("acmgPathogenicityClassification", "")
                    vu = var_uri(pid, gene_sym)

                    fh.write(f"\n{vu}\n")
                    fh.write(f"  a pavs:Variant ;\n")
                    fh.write(f"  pavs:affectsGene hgnc:{gene_sym} ;\n")
                    if acmg:
                        fh.write(f"  pavs:acmgClass {ttl_str(acmg)} ;\n")

                    # HGVS expressions
                    for expr in vd.get("expressions", []):
                        syntax = expr.get("syntax", "")
                        value = expr.get("value", "")
                        if syntax == "hgvs.c" and value:
                            fh.write(f"  pavs:hgvsC {ttl_str(value)} ;\n")
                        elif syntax == "hgvs.g" and value:
                            fh.write(f"  pavs:hgvsG {ttl_str(value)} ;\n")

                    # VCF record
                    vcf = vd.get("vcfRecord", {})
                    if vcf.get("chrom"):
                        fh.write(f"  pavs:vcfChrom {ttl_str(vcf['chrom'])} ;\n")
                        fh.write(f"  pavs:vcfPos {vcf.get('pos', 0)} ;\n")
                        fh.write(f"  pavs:vcfRef {ttl_str(vcf.get('ref', ''))} ;\n")
                        fh.write(f"  pavs:vcfAlt {ttl_str(vcf.get('alt', ''))} ;\n")

                    fh.write(f"  rdfs:label {ttl_str(gene_sym)} .\n")

    print(f"Wrote cases.ttl → {output_path}")


# ---------------------------------------------------------------------------
# Generate genes.ttl from combined_annotated.tsv
# ---------------------------------------------------------------------------

def generate_genes_ttl(annotated_tsv: str, output_path: str):
    print(f"Loading annotated data from {annotated_tsv} …")
    df = pd.read_csv(annotated_tsv, sep="\t", dtype=str, low_memory=False)
    print(f"  {len(df)} rows, {len(df.columns)} columns")

    genes_seen: set = set()

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(PREFIXES)

        for _, row in df.iterrows():
            gene = row.get("gene_symbol", "")
            if not safe(gene) or gene in genes_seen:
                continue
            genes_seen.add(gene)

            fh.write(f"\nhgnc:{gene}\n")
            fh.write(f"  a pavs:Gene ;\n")
            fh.write(f"  rdfs:label {ttl_str(gene)} ;\n")

            entrez = row.get("gene_entrez_id", "")
            if safe(entrez):
                fh.write(f"  pavs:ncbiGeneId ncbigene:{entrez} ;\n")

            pli = row.get("gnomad_pli", "")
            if safe(pli):
                fh.write(f"  pavs:pli {pli} ;\n")

            loeuf = row.get("gnomad_loeuf", "")
            if safe(loeuf):
                fh.write(f"  pavs:loeuf {loeuf} ;\n")

            go_bp = row.get("go_biological_process", "")
            if safe(go_bp):
                fh.write(f"  pavs:goProcess {ttl_str(go_bp)} ;\n")

            expressed = row.get("expressed_in", "")
            if safe(expressed):
                fh.write(f"  pavs:expressedIn {ttl_str(expressed[:500])} ;\n")

            disease_omim = row.get("gene_disease_omim", "")
            if safe(disease_omim):
                for d in str(disease_omim).split("|"):
                    d = d.strip()
                    if d.startswith("OMIM:"):
                        fh.write(f"  pavs:relatedDisease {omim_uri(d)} ;\n")

            fh.write(f"  rdfs:seeAlso {ttl_uri(f'http://identifiers.org/hgnc.symbol/{gene}')} .\n")

        # Also write variant-level data linking back to cases
        fh.write("\n# Variant annotations from combined_annotated.tsv\n")
        for _, row in df.iterrows():
            pavs_id = row.get("pavs_id", "")
            gene = row.get("gene_symbol", "")
            if not safe(pavs_id) or not safe(gene):
                continue

            cu = case_uri(pavs_id)
            vu = var_uri(pavs_id, gene)

            rs_id = row.get("vep_rsid", "")
            hgvsc = row.get("vep_hgvsc", row.get("variant_cdna", ""))
            hgvsp = row.get("vep_hgvsp", row.get("variant_protein", ""))
            hgvsg = row.get("vep_hgvsg", row.get("variant_genomic_grch38", ""))
            zyg_id = row.get("zygosity_id", "")
            zyg_label = row.get("zygosity_label", "")
            acmg = row.get("acmg_classification", "")
            consequence = row.get("vep_consequence", "")
            sift = row.get("vep_sift", "")
            polyphen = row.get("vep_polyphen", "")
            gnomad_af = row.get("vep_af_gnomad", "")
            af_saudi = row.get("af_saudi", "")
            ac_saudi = row.get("ac_saudi", "")
            clinvar_id = row.get("clinvar_allele_id", "")
            clinvar_sig = row.get("clinvar_sig", "")
            source = row.get("source_file", "")
            pmid = row.get("pmid", "")

            is_saudi = source in SAUDI_SOURCES

            fh.write(f"\n{vu}\n")
            fh.write(f"  a pavs:Variant ;\n")
            fh.write(f"  pavs:affectsGene hgnc:{gene} ;\n")
            fh.write(f"  pavs:isSaudi {str(is_saudi).lower()} ;\n")

            if safe(hgvsc):
                fh.write(f"  pavs:hgvsC {ttl_str(hgvsc)} ;\n")
            if safe(hgvsp):
                fh.write(f"  pavs:hgvsP {ttl_str(hgvsp)} ;\n")
            if safe(hgvsg):
                fh.write(f"  pavs:hgvsG {ttl_str(hgvsg)} ;\n")
            if safe(rs_id):
                rs_num = str(rs_id).replace("rs", "").strip()
                if rs_num.isdigit():
                    fh.write(f"  pavs:rsId dbsnp:rs{rs_num} ;\n")
                    fh.write(f"  rdfs:seeAlso dbsnp:rs{rs_num} ;\n")
            if safe(zyg_id):
                fh.write(f"  pavs:zygosity {geno_uri(zyg_id)} ;\n")
            if safe(acmg):
                fh.write(f"  pavs:acmgClass {ttl_str(acmg)} ;\n")
            if safe(consequence):
                fh.write(f"  pavs:vepConsequence {ttl_str(consequence)} ;\n")
            if safe(sift):
                fh.write(f"  pavs:sift {ttl_str(sift)} ;\n")
            if safe(polyphen):
                fh.write(f"  pavs:polyphen {ttl_str(polyphen)} ;\n")
            if safe(gnomad_af):
                fh.write(f"  pavs:gnomadAF {gnomad_af} ;\n")
            if safe(af_saudi):
                fh.write(f"  pavs:saudiAF {af_saudi} ;\n")
            if safe(ac_saudi):
                fh.write(f"  pavs:saudiAC {ac_saudi} ;\n")
            if safe(clinvar_id):
                fh.write(f"  pavs:clinvarId {ttl_str(clinvar_id)} ;\n")
                fh.write(f"  rdfs:seeAlso clinvar:{clinvar_id} ;\n")
            if safe(clinvar_sig):
                fh.write(f"  pavs:clinvarSig {ttl_str(clinvar_sig)} ;\n")
            if safe(pmid):
                fh.write(f"  pavs:pmid {ttl_str(str(pmid).strip())} ;\n")

            fh.write(f"  rdfs:label {ttl_str(gene)} .\n")

    print(f"Wrote genes.ttl → {output_path} ({len(genes_seen)} genes)")


# ---------------------------------------------------------------------------
# Generate hpoa.ttl from phenotype.hpoa
# ---------------------------------------------------------------------------

def generate_hpoa_ttl(hpoa_path: str, output_path: str):
    print(f"Loading HPOA from {hpoa_path} …")
    count = 0
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(PREFIXES)

        with open(hpoa_path, encoding="utf-8") as hfh:
            for line in hfh:
                if line.startswith("#") or line.startswith("database_id"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                disease_id = parts[0]
                qualifier = parts[2]
                hpo_id = parts[3]
                evidence = parts[5] if len(parts) > 5 else ""
                freq = parts[7] if len(parts) > 7 else ""

                if qualifier == "NOT" or not hpo_id.startswith("HP:"):
                    continue

                assoc_uri = f"pav:hpoa_{disease_id.replace(':', '_')}_{hpo_id.replace(':', '_')}"
                disease_uri = ""
                if disease_id.startswith("OMIM:"):
                    disease_uri = omim_uri(disease_id)
                elif disease_id.startswith("MONDO:"):
                    disease_uri = mondo_uri(disease_id)
                else:
                    disease_uri = ttl_uri(f"http://identifiers.org/{disease_id}")

                fh.write(f"\n{assoc_uri}\n")
                fh.write(f"  a pavs:DiseaseHPOAssociation ;\n")
                fh.write(f"  pavs:disease {disease_uri} ;\n")
                fh.write(f"  pavs:hpoTerm {hp_uri(hpo_id)} ;\n")
                if evidence:
                    fh.write(f"  pavs:evidence {ttl_str(evidence)} ;\n")
                if freq:
                    fh.write(f"  pavs:frequency {ttl_str(freq)} ;\n")
                fh.write(f"  rdfs:label {ttl_str(f'{disease_id}→{hpo_id}')} .\n")
                count += 1

    print(f"Wrote hpoa.ttl → {output_path} ({count} associations)")


# ---------------------------------------------------------------------------
# Generate literature.ttl from phenopackets/0.1.26/
# ---------------------------------------------------------------------------

def generate_literature_ttl(lit_dir: str, output_path: str):
    lit_path = Path(lit_dir)
    json_files = list(lit_path.rglob("*.json"))
    print(f"Found {len(json_files)} literature phenopackets in {lit_dir}")

    count = 0
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(PREFIXES)

        for jf in json_files:
            try:
                with open(jf, encoding="utf-8") as jfh:
                    pp = json.load(jfh)
            except Exception as e:
                print(f"  WARN: could not parse {jf}: {e}")
                continue

            pid = pp.get("id", jf.stem)
            if not pid:
                continue

            safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", pid)
            cu = f"pav:lit_{safe_id}"

            fh.write(f"\n{cu}\n")
            fh.write(f"  a pavs:Case ;\n")
            fh.write(f"  dc:identifier {ttl_str(pid)} ;\n")
            fh.write(f"  pavs:source {ttl_str('Literature')} ;\n")
            fh.write(f"  pavs:isSaudi false ;\n")

            # Sex
            subject = pp.get("subject", {})
            sex = subject.get("sex", "")
            sex_map = {
                "MALE": "http://purl.obolibrary.org/obo/NCIT_C20197",
                "FEMALE": "http://purl.obolibrary.org/obo/NCIT_C16576",
            }
            if sex in sex_map:
                fh.write(f"  pavs:sex {ttl_uri(sex_map[sex])} ;\n")

            # Phenotypes
            for feat in pp.get("phenotypicFeatures", []):
                ftype = feat.get("type", {})
                fid = ftype.get("id", "")
                excluded = feat.get("excluded", False)
                if fid.startswith("HP:"):
                    pred = "pavs:hasExcludedPhenotype" if excluded else "pavs:hasPhenotype"
                    fh.write(f"  {pred} {hp_uri(fid)} ;\n")

            # Diseases
            for disease in pp.get("diseases", []):
                did = disease.get("term", {}).get("id", "")
                dlabel = disease.get("term", {}).get("label", "")
                if did.startswith("OMIM:"):
                    fh.write(f"  pavs:hasDisease {omim_uri(did)} ;\n")
                elif did.startswith("MONDO:"):
                    fh.write(f"  pavs:hasDisease {mondo_uri(did)} ;\n")
                if dlabel:
                    fh.write(f"  pavs:diseaseLabel {ttl_str(dlabel)} ;\n")

            # Interpretations → variants
            for interp in pp.get("interpretations", []):
                progress = interp.get("progressStatus", "")
                if progress:
                    fh.write(f"  pavs:progressStatus {ttl_str(progress)} ;\n")

                for gi in interp.get("diagnosis", {}).get("genomicInterpretations", []):
                    vi = gi.get("variantInterpretation", {})
                    vd = vi.get("variationDescriptor", {})
                    gene_ctx = vd.get("geneContext", {})
                    gene_sym = gene_ctx.get("symbol", "unknown")
                    acmg = vi.get("acmgPathogenicityClassification", "")

                    v_safe = re.sub(r"[^A-Za-z0-9_.-]", "_", pid)
                    g_safe = re.sub(r"[^A-Za-z0-9_.-]", "_", gene_sym)
                    vu = f"pav:litvar_{v_safe}_{g_safe}"

                    fh.write(f"  pavs:hasVariant {vu} ;\n")

                    # Write variant block after
                    # (collect and write below to avoid interrupting subject block)

            fh.write(f"  rdfs:label {ttl_str(pid)} .\n")

            # Now write variant triples
            for interp in pp.get("interpretations", []):
                for gi in interp.get("diagnosis", {}).get("genomicInterpretations", []):
                    vi = gi.get("variantInterpretation", {})
                    vd = vi.get("variationDescriptor", {})
                    gene_ctx = vd.get("geneContext", {})
                    gene_sym = gene_ctx.get("symbol", "unknown")
                    acmg = vi.get("acmgPathogenicityClassification", "")

                    v_safe = re.sub(r"[^A-Za-z0-9_.-]", "_", pid)
                    g_safe = re.sub(r"[^A-Za-z0-9_.-]", "_", gene_sym)
                    vu = f"pav:litvar_{v_safe}_{g_safe}"

                    fh.write(f"\n{vu}\n")
                    fh.write(f"  a pavs:Variant ;\n")
                    fh.write(f"  pavs:affectsGene hgnc:{gene_sym} ;\n")
                    fh.write(f"  pavs:isSaudi false ;\n")
                    if acmg:
                        fh.write(f"  pavs:acmgClass {ttl_str(acmg)} ;\n")

                    for expr in vd.get("expressions", []):
                        syntax = expr.get("syntax", "")
                        value = expr.get("value", "")
                        if syntax == "hgvs.c" and value:
                            fh.write(f"  pavs:hgvsC {ttl_str(value)} ;\n")
                        elif syntax == "hgvs.g" and value:
                            fh.write(f"  pavs:hgvsG {ttl_str(value)} ;\n")

                    vcf = vd.get("vcfRecord", {})
                    if vcf.get("chrom"):
                        fh.write(f"  pavs:vcfChrom {ttl_str(vcf['chrom'])} ;\n")
                        pos = vcf.get("pos", 0)
                        fh.write(f"  pavs:vcfPos {pos} ;\n")
                        fh.write(f"  pavs:vcfRef {ttl_str(vcf.get('ref', ''))} ;\n")
                        fh.write(f"  pavs:vcfAlt {ttl_str(vcf.get('alt', ''))} ;\n")

                    fh.write(f"  rdfs:label {ttl_str(gene_sym)} .\n")

            count += 1
            if count % 1000 == 0:
                print(f"  … {count} literature phenopackets written")

    print(f"Wrote literature.ttl → {output_path} ({count} phenopackets)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate RDF from PAVS data sources")
    parser.add_argument("--phenopackets", default="data/PAVS_phenopackets.json")
    parser.add_argument("--annotated", default="data/combined_annotated.tsv")
    parser.add_argument("--hpoa", default="data/reference/phenotype.hpoa")
    parser.add_argument("--literature-dir", default="phenopackets/0.1.26")
    parser.add_argument("--output-dir", default="rdf_output/")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    generate_cases_ttl(args.phenopackets, str(out / "cases.ttl"))
    generate_genes_ttl(args.annotated, str(out / "genes.ttl"))
    generate_hpoa_ttl(args.hpoa, str(out / "hpoa.ttl"))
    generate_literature_ttl(args.literature_dir, str(out / "literature.ttl"))

    print("\nAll TTL files generated successfully.")


if __name__ == "__main__":
    main()
