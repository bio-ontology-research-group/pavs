#!/usr/bin/env python3
"""
combine_normalize_phenotypes.py

Reads all 7 source files in data/phenotypes/, normalizes every field to
standard ontologies (HPO, GENO, NCIT, MONDO, OMIM, HGVS), and writes a
single combined TSV: data/combined_normalized.tsv

Also generates: ontology/acmg_criteria.obo

Usage:
    uv run python scripts/combine_normalize_phenotypes.py \\
        --output data/combined_normalized.tsv \\
        --acmg-obo ontology/acmg_criteria.obo \\
        --workers 4 \\
        [--limit N]
"""

import argparse
import csv
import os
import re
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple

import pandas as pd
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHENOTYPES_DIR = "data/phenotypes"

OUTPUT_COLUMNS = [
    "pavs_id",
    "source_file",
    "source_id",
    "sex_id",
    "sex_label",
    "age_iso",
    "consanguinity_id",
    "consanguinity_label",
    "family_history_id",
    "family_history_label",
    "test_type",
    "result_status",
    "inheritance_id",
    "inheritance_label",
    "phenotypes_present_ids",
    "phenotypes_present_labels",
    "phenotypes_excluded_ids",
    "phenotypes_excluded_labels",
    "phenotypes_modifier_ids",
    "phenotypes_modifier_labels",
    "phenotypes_raw",
    "gene_symbol",
    "gene_entrez_id",
    "variant_cdna",
    "variant_protein",
    "variant_genomic_grch38",
    "variant_hgvs_valid",
    "zygosity_id",
    "zygosity_label",
    "acmg_classification",
    "acmg_evidence",
    "disease_mondo_id",
    "disease_mondo_label",
    "disease_omim_ids",
    "disease_omim_labels",
    "disease_orphanet_ids",
    "pmid",
    "notes",
]

# Source file letter codes for PAVS IDs
SOURCE_LETTERS = {
    "ahmed-variants": "A",
    "ahmed-pmid28454995": "B",
    "fawzan-variants": "F",
    "marwa-variants": "M",
    "PMC6562004": "P",
    "PMC7082194": "Q",
    "ddd-diagnoses": "D",
}

# Amino acid 1-letter to 3-letter conversion
AMINO_ACID_MAP = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter", "X": "Xaa",
}

# HPO inheritance mapping
INHERITANCE_MAP = {
    # AR / autosomal recessive
    r"\bAR\b": ("HP:0000007", "Autosomal recessive inheritance"),
    r"autosomal[\s-]*recessive": ("HP:0000007", "Autosomal recessive inheritance"),
    r"\bbiallelic\b": ("HP:0000007", "Autosomal recessive inheritance"),
    # AD / autosomal dominant
    r"\bAD\b": ("HP:0000006", "Autosomal dominant inheritance"),
    r"autosomal[\s-]*dominant": ("HP:0000006", "Autosomal dominant inheritance"),
    r"\bmonoallelic\b": ("HP:0000006", "Autosomal dominant inheritance"),
    # XLD / X-linked dominant (before generic XL)
    r"\bXLD\b": ("HP:0001423", "X-linked dominant inheritance"),
    r"x[\s-]*linked[\s-]*dominant": ("HP:0001423", "X-linked dominant inheritance"),
    # XLR / X-linked recessive (before generic XL)
    r"\bXLR\b": ("HP:0001419", "X-linked recessive inheritance"),
    r"x[\s-]*linked[\s-]*recessive": ("HP:0001419", "X-linked recessive inheritance"),
    # Generic X-linked
    r"\bXL\b": ("HP:0001417", "X-linked inheritance"),
    r"x[\s-]*linked": ("HP:0001417", "X-linked inheritance"),
    # Mitochondrial
    r"\bMT\b": ("HP:0001427", "Mitochondrial inheritance"),
    r"mitochondrial": ("HP:0001427", "Mitochondrial inheritance"),
    # Sporadic
    r"\bsporadic\b": ("HP:0003745", "Sporadic"),
    # Y-linked
    r"y[\s-]*linked": ("HP:0001450", "Y-linked inheritance"),
    # De novo
    r"de[\s-]*novo": ("HP:0025352", "Typically de novo"),
    # Semidominant
    r"semidominant": ("HP:0032113", "Semidominant inheritance"),
}

# PMC6562004 primary indication map (numeric code → text)
PRIMARY_INDICATION_MAP = {
    "1": "Developmental Delay/Intellectual Disability",
    "2": "Dysmorphism/Congenital Malformations",
    "3": "Neuromuscular",
    "4": "Prenatal",
    "5": "Immunodeficiency",
    "6": "Liver Disease",
    "7": "Metabolic",
    "8": "Pulmonary",
    "9": "Retinal dystrophy",
    "10": "Deafness",
    "11": "Others",
}

# Zygosity keywords → GENO IDs
ZYGOSITY_KEYWORD_MAP = [
    (r"\bhemizygous\b", "GENO:0000134", "hemizygous"),
    (r"\bhemi\b", "GENO:0000134", "hemizygous"),
    (r"\bhomozygous\b", "GENO:0000136", "homozygous"),
    (r"\bhomo\b", "GENO:0000136", "homozygous"),
    (r"\bhom\b", "GENO:0000136", "homozygous"),
    (r"\bheterozygous\b", "GENO:0000135", "heterozygous"),
    (r"\bhet\b", "GENO:0000135", "heterozygous"),
    (r"\bcompound[\s-]*het", "GENO:0000402", "compound heterozygous"),
]

# ACMG criteria definitions for OBO generation
ACMG_CRITERIA = [
    ("PVS1", "Null variant in a gene where loss of function is a known mechanism of disease.",
     "strength=very_strong; direction=pathogenic"),
    ("PS1", "Same amino acid change as a previously established pathogenic variant.",
     "strength=strong; direction=pathogenic"),
    ("PS2", "De novo variant (both maternity and paternity confirmed) in a patient with the disease and no family history.",
     "strength=strong; direction=pathogenic"),
    ("PS3", "Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product.",
     "strength=strong; direction=pathogenic"),
    ("PS4", "The prevalence of the variant in affected individuals is significantly increased compared with the prevalence in controls.",
     "strength=strong; direction=pathogenic"),
    ("PM1", "Located in a mutational hot spot and/or critical and well-established functional domain without benign variation.",
     "strength=moderate; direction=pathogenic"),
    ("PM2", "Absent from controls (or at extremely low frequency if recessive) in population databases.",
     "strength=moderate; direction=pathogenic"),
    ("PM3", "For recessive disorders, detected in trans with a pathogenic variant.",
     "strength=moderate; direction=pathogenic"),
    ("PM4", "Protein length changes as a result of in-frame deletions/insertions in a nonrepeat region or stop-loss variants.",
     "strength=moderate; direction=pathogenic"),
    ("PM5", "Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before.",
     "strength=moderate; direction=pathogenic"),
    ("PM6", "Assumed de novo, but without confirmation of paternity and maternity.",
     "strength=moderate; direction=pathogenic"),
    ("PP1", "Cosegregation with disease in multiple affected family members in a gene definitively known to cause the disease.",
     "strength=supporting; direction=pathogenic"),
    ("PP2", "Missense variant in a gene that has a low rate of benign missense variation and in which missense variants are a common mechanism of disease.",
     "strength=supporting; direction=pathogenic"),
    ("PP3", "Multiple lines of computational evidence support a deleterious effect on the gene or gene product.",
     "strength=supporting; direction=pathogenic"),
    ("PP4", "Patient's phenotype or family history is highly specific for a disease with a single genetic etiology.",
     "strength=supporting; direction=pathogenic"),
    ("PP5", "Reputable source recently reports variant as pathogenic, but the evidence is not available to the laboratory to perform an independent evaluation.",
     "strength=supporting; direction=pathogenic"),
    ("BA1", "Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium.",
     "strength=stand_alone; direction=benign"),
    ("BS1", "Allele frequency is greater than expected for disorder.",
     "strength=strong; direction=benign"),
    ("BS2", "Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked (hemizygous) disorder with full penetrance expected at an early age.",
     "strength=strong; direction=benign"),
    ("BS3", "Well-established in vitro or in vivo functional studies show no damaging effect on protein function or splicing.",
     "strength=strong; direction=benign"),
    ("BS4", "Lack of segregation in affected members of a family.",
     "strength=strong; direction=benign"),
    ("BP1", "Missense variant in a gene for which primarily truncating variants are known to cause disease.",
     "strength=supporting; direction=benign"),
    ("BP2", "Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder or observed in cis with a pathogenic variant in any inheritance pattern.",
     "strength=supporting; direction=benign"),
    ("BP3", "In-frame deletions/insertions in a repetitive region without a known function.",
     "strength=supporting; direction=benign"),
    ("BP4", "Multiple lines of computational evidence suggest no impact on gene or gene product.",
     "strength=supporting; direction=benign"),
    ("BP5", "Variant found in a case with an alternate molecular basis for disease.",
     "strength=supporting; direction=benign"),
    ("BP6", "Reputable source recently reports variant as benign, but the evidence is not available to the laboratory to perform an independent evaluation.",
     "strength=supporting; direction=benign"),
    ("BP7", "A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the splice consensus sequence nor the creation of a new splice site.",
     "strength=supporting; direction=benign"),
]


# ---------------------------------------------------------------------------
# HGVS Validator
# ---------------------------------------------------------------------------

def preprocess_hgvs(s: str) -> str:
    """
    Normalize common formatting issues in HGVS strings before validation.

    Fixes:
    - Unicode/typographic minus characters → ASCII hyphen
    - Leading/trailing parentheses and spaces
    - Non-standard operators: "->" → ">", "/" → ">"
    - Spaces around operators within HGVS notation (e.g. "c.80 G > A" → "c.80G>A")
    """
    if not s:
        return s
    # Unicode/typographic minus variants → ASCII hyphen
    s = s.translate(str.maketrans("\u2010\u2011\u2012\u2013\u2014", "-----"))
    # Strip leading/trailing whitespace and parentheses
    s = s.lstrip("( ").rstrip(") ")
    # Replace non-standard operators
    s = re.sub(r"\s*->\s*", ">", s)   # -> → >
    s = re.sub(r"([A-Za-z0-9])\s*/\s*([A-Za-z])", r"\1>\2", s)  # allele slash: A/G → A>G
    # Remove spaces around operators within HGVS notation
    s = re.sub(r"(c|g|m|n|r)\.\s+", r"\1.", s)     # space after prefix dot: "c. 80" → "c.80"
    s = re.sub(r"([A-Za-z0-9])\s+>", r"\1>", s)    # space before >
    s = re.sub(r">\s+([A-Za-z])", r">\1", s)        # space after >
    s = re.sub(r"\s+([A-Za-z])>", r"\1>", s)        # leading space before nucleotide+>
    return s.strip()


class HgvsValidator:
    """
    HGVS format validation. Tries biocommons.hgvs.parser first;
    falls back to regex if unavailable.
    """

    def __init__(self):
        self._parser = None
        try:
            from biocommons.hgvs.parser import Parser as HgvsParser
            self._parser = HgvsParser()
        except Exception:
            pass

    def validate(self, variant_str: str) -> bool:
        """Validate a variant string, preprocessing it first and handling composite strings."""
        if not variant_str or not variant_str.strip():
            return False
        s = preprocess_hgvs(variant_str.strip())
        # Try each comma-separated component (handles "c.X,p.Y,NC_..." composites)
        for part in [p.strip() for p in s.split(",")]:
            if not part:
                continue
            if self._validate_single(part):
                return True
        return False

    def _validate_single(self, s: str) -> bool:
        """Validate a single (non-composite) HGVS string."""
        if self._parser is not None:
            try:
                self._parser.parse(s)
                return True
            except Exception:
                pass
        # Regex fallback
        # Strip transcript prefix: NM_000291.3:c.758T>C → c.758T>C
        if ":" in s:
            s = s.rsplit(":", 1)[-1].strip()
        cdna_re = re.compile(r"^c\.\d+[\+\-]?\d*[A-Za-z>delinsdup_\[\];*]+", re.I)
        prot_re = re.compile(r"^p\.[A-Z][a-z]{2}\d+", re.I)
        gen_re  = re.compile(r"^g\.\d+[A-Za-z>delinsdup_]+", re.I)
        mito_re = re.compile(r"^m\.\d+[A-Za-z>delinsdup_]+", re.I)
        return bool(cdna_re.match(s) or prot_re.match(s) or gen_re.match(s) or mito_re.match(s))


# ---------------------------------------------------------------------------
# Normalization helpers
# ---------------------------------------------------------------------------

def normalize_age(text) -> str:
    """Convert various age text formats to ISO 8601 duration string."""
    if not text or str(text).strip() in ("", "nan", "?", "N/A", "NA", "Unknown"):
        return ""
    t = str(text).strip()

    # "{N}Y/{N}M/{N}D" format: e.g. "22Y" or "2Y6M" or "2Y6M15D"
    m = re.match(r"^(\d+)[Yy](?:(\d+)[Mm])?(?:(\d+)[Dd])?$", t)
    if m:
        y, mo, d = m.group(1), m.group(2), m.group(3)
        result = "P" + y + "Y"
        if mo:
            result += mo + "M"
        if d:
            result += d + "D"
        return result

    # "{N}y {N}m" format: "2y 11m"
    m = re.match(r"^(\d+)[Yy]\s+(\d+)[Mm]$", t)
    if m:
        return f"P{m.group(1)}Y{m.group(2)}M"

    # "{N}y" format: "5y"
    m = re.match(r"^(\d+)[Yy]$", t)
    if m:
        return f"P{m.group(1)}Y"

    # "N years N months" format
    m = re.match(r"^(\d+)\s+years?\s+(\d+)\s+months?", t, re.I)
    if m:
        return f"P{m.group(1)}Y{m.group(2)}M"

    # "N years" format
    m = re.match(r"^(\d+)\s+years?", t, re.I)
    if m:
        return f"P{m.group(1)}Y"

    # "N months" format
    m = re.match(r"^(\d+)\s+months?", t, re.I)
    if m:
        return f"P{m.group(1)}M"

    # "N Year" or "N Years" (from PMC6562004)
    m = re.match(r"^(\d+)\s+[Yy]ear", t, re.I)
    if m:
        return f"P{m.group(1)}Y"

    # "N Month" (from PMC6562004)
    m = re.match(r"^(\d+)\s+[Mm]onth", t, re.I)
    if m:
        return f"P{m.group(1)}M"

    # Plain integer → assume years
    m = re.match(r"^(\d+)$", t)
    if m:
        return f"P{m.group(1)}Y"

    return ""


def normalize_inheritance(text) -> Tuple[List[str], List[str]]:
    """
    Map inheritance text to HPO IDs/labels using INHERITANCE_MAP.
    Returns (hpo_ids, hpo_labels). Multiple values are possible.
    Zygosity keywords (hom, het, hemi) are NOT mapped to HPO inheritance.
    """
    if not text or str(text).strip() in ("", "nan", "N/A", "NA"):
        return [], []

    t = str(text).strip()
    # Split on common delimiters
    tokens = re.split(r"[|;,/\n]", t)

    hpo_ids: List[str] = []
    hpo_labels: List[str] = []
    seen: set = set()

    for tok in tokens:
        tok = tok.strip()
        if not tok:
            continue
        tok_lower = tok.lower()
        # Skip pure zygosity keywords
        if re.match(r"^(hom|het|hemi|homo|homozygous|heterozygous|hemizygous)$", tok_lower):
            continue

        # Try INHERITANCE_MAP patterns
        for pattern, (hpo_id, label) in INHERITANCE_MAP.items():
            if re.search(pattern, tok, re.I):
                if hpo_id not in seen:
                    seen.add(hpo_id)
                    hpo_ids.append(hpo_id)
                    hpo_labels.append(label)
                break  # First matching pattern wins for this token

    return hpo_ids, hpo_labels


def normalize_result_status(text) -> str:
    """Map raw result/status text to one of: Solved / Unsolved / Undetermined / Unknown."""
    if not text or str(text).strip() in ("", "nan", "N/A"):
        return "Unknown"
    t = str(text).strip().lower()
    if any(x in t for x in ["positive", "solved", " y", "=y"]) and "negative" not in t and "unsolved" not in t:
        return "Solved"
    if any(x in t for x in ["negative", "unsolved", " n", "=n"]) and "positive" not in t and "solved" not in t:
        return "Unsolved"
    if any(x in t for x in ["ambiguous", "inconclusive", "unclear"]):
        return "Undetermined"
    # Check single character values
    if t.strip() == "y":
        return "Solved"
    if t.strip() == "n":
        return "Unsolved"
    if t.strip() == "solved":
        return "Solved"
    if t.strip() == "unsolved":
        return "Unsolved"
    return "Unknown"


def normalize_test_type(text) -> str:
    """Map test type text to: WES / WGS / gene_panel / microarray / other."""
    if not text or str(text).strip() in ("", "nan"):
        return ""
    t = str(text).strip().lower()
    if "wes" in t or "exome" in t:
        return "WES"
    if "wgs" in t or "genome" in t or "whole genome" in t:
        return "WGS"
    if "panel" in t:
        return "gene_panel"
    if "microarray" in t or "array" in t or "snp" in t or "chip" in t:
        return "microarray"
    return "other"


def convert_protein_1letter_to_3letter(p_str: str) -> str:
    """
    Convert single-letter amino acid HGVS protein notation to 3-letter.
    e.g. p.R200W → p.Arg200Trp; p.R58X → p.Arg58Ter
    Also handles compound: p.R200W;p.N197S
    """
    if not p_str or not p_str.strip():
        return p_str

    def convert_one(s: str) -> str:
        s = s.strip()
        # Already 3-letter: p.Arg200Trp etc.
        if re.match(r"^p\.[A-Z][a-z]{2}", s):
            return s
        # 1-letter: p.R200W
        m = re.match(r"^(p\.)([A-Z\*])(\d+)([A-Z\*?=]?)$", s)
        if m:
            prefix, aa1, pos, aa2 = m.groups()
            aa1_3 = AMINO_ACID_MAP.get(aa1, aa1)
            aa2_3 = AMINO_ACID_MAP.get(aa2, aa2) if aa2 else ""
            return f"{prefix}{aa1_3}{pos}{aa2_3}"
        # frameshift: p.K1373fs
        m = re.match(r"^(p\.)([A-Z\*])(\d+)(fs.*)$", s)
        if m:
            prefix, aa1, pos, suffix = m.groups()
            aa1_3 = AMINO_ACID_MAP.get(aa1, aa1)
            return f"{prefix}{aa1_3}{pos}{suffix}"
        return s

    # Handle semicolon-separated
    if ";" in p_str:
        parts = p_str.split(";")
        return ";".join(convert_one(p) for p in parts)
    return convert_one(p_str)


def extract_pmid(text) -> List[str]:
    """Extract PubMed IDs from text (URLs or plain integers)."""
    if not text or str(text).strip() in ("", "nan"):
        return []
    t = str(text)
    # Extract all numeric sequences that look like PMIDs (5-9 digits)
    pmids = re.findall(r"(?:pubmed/|PMID[:\s]*|^)(\d{5,9})", t, re.I)
    if not pmids:
        # Try plain integers
        pmids = re.findall(r"\b(\d{5,9})\b", t)
    return list(dict.fromkeys(pmids))  # deduplicate preserving order


def clean_gene_symbol(raw: str) -> List[str]:
    """
    Pre-normalize a raw gene field to a list of candidate HGNC symbols.

    Handles:
    - Multi-gene comma/semicolon-separated fields: "PKLR, WRAP53" → ["PKLR", "WRAP53"]
    - Parenthetical chromosome annotations: "DEAF1 (CHR11" → ["DEAF1"]
    - Trailing RefSeq transcript IDs: "DVL2 NM_004422" → ["DVL2"]
    - Fields that are clearly not gene symbols (e.g. ".", "-", genome builds)
    """
    raw = str(raw).strip()
    if not raw or raw in (".", "-", "N/A", "NA", ""):
        return []
    # Split multi-gene fields (comma or semicolon)
    parts = re.split(r"[,;]", raw)
    results = []
    for p in parts:
        p = p.strip()
        if not p:
            continue
        # Remove parenthetical annotations: "DEAF1 (CHR11" → "DEAF1"
        p = re.sub(r"\s*\(.*", "", p).strip()
        # Remove trailing RefSeq/transcript IDs: "DVL2 NM_004422" → "DVL2"
        p = re.split(r"\s+(?:NM_|NR_|NG_|NC_|XM_|XR_)", p)[0].strip()
        # Must look like a gene symbol: letters, digits, hyphens only, starts with letter
        if re.match(r"^[A-Za-z][A-Za-z0-9\-]+$", p):
            results.append(p.upper())
    return results


def normalize_gene(raw_symbol: str, genes_dict: dict) -> Tuple[str, str]:
    """
    Look up gene symbol in NormalizationUtils.genes dict.
    Returns (symbol, entrez_id) as pipe-delimited strings for multi-gene inputs.
    Falls back to cleaned candidate without Entrez ID on lookup failure.
    """
    if not raw_symbol or str(raw_symbol).strip() in ("", "nan"):
        return ("", "")
    candidates = clean_gene_symbol(str(raw_symbol))
    if not candidates:
        return ("", "")
    symbols: List[str] = []
    entrezs: List[str] = []
    for sym in candidates:
        info = genes_dict.get(sym)
        if info:
            symbols.append(info.get("symbol", sym))
            entrezs.append(str(info.get("entrez", "")))
        else:
            symbols.append(sym)
            entrezs.append("")
    return ("|".join(symbols), "|".join(entrezs))


def normalize_zygosity(text) -> Tuple[str, str]:
    """Map zygosity text to GENO ID and label."""
    if not text or str(text).strip() in ("", "nan"):
        return ("GENO:0000137", "unspecified zygosity")
    t = str(text).strip().lower()
    for pattern, geno_id, label in ZYGOSITY_KEYWORD_MAP:
        if re.search(pattern, t, re.I):
            return (geno_id, label)
    return ("GENO:0000137", "unspecified zygosity")


def extract_zygosity_from_inheritance(text) -> Optional[Tuple[str, str]]:
    """
    Extract zygosity keywords from inheritance text.
    Returns (geno_id, label) if found, else None.
    """
    if not text:
        return None
    t = str(text).lower()
    for pattern, geno_id, label in ZYGOSITY_KEYWORD_MAP:
        if re.search(pattern, t, re.I):
            return (geno_id, label)
    return None


def infer_compound_het(variant_entries: List[Tuple[str, str]]) -> List[str]:
    """
    Given list of (gene, zygosity_id) for one patient,
    returns updated zygosity_ids: replaces GENO:0000135 with GENO:0000402
    for variants in the same gene where ≥2 are heterozygous.
    """
    from collections import Counter
    gene_counts: Counter = Counter()
    for gene, zyg in variant_entries:
        if gene and zyg == "GENO:0000135":
            gene_counts[gene.upper()] += 1

    result = []
    for gene, zyg in variant_entries:
        if gene and gene.upper() in gene_counts and gene_counts[gene.upper()] >= 2 and zyg == "GENO:0000135":
            result.append("GENO:0000402")
        else:
            result.append(zyg)
    return result


def parse_variant_string(raw: str) -> Tuple[str, str, str, str]:
    """
    Parse various variant string formats.
    Returns (cdna, protein, genomic, notes_fragment).
    Handles:
    - GENE:NM_xxxxx:ExonN:c.xxx:p.xxx
    - NM_xxxxx:c.xxx:p.xxx
    - NC_xxxxx.N:g.xxxxA>G
    - ChrX(GRCh38):g.xxxxA>G
    - c.xxx:p.xxx
    """
    if not raw or str(raw).strip() in ("", "nan"):
        return ("", "", "", "")

    raw = str(raw).strip()
    cdna = protein = genomic = notes_frag = ""

    # NC_ format (genomic reference): NC_000012.11:g.49425536T>C
    if re.match(r"NC_\d+\.\d+:", raw):
        m = re.search(r"(g\.\d+[A-Za-z>delinsdup_]+)", raw)
        if m:
            genomic = raw  # keep full string as genomic
        else:
            notes_frag = raw
        return (cdna, protein, genomic, notes_frag)

    # ChrX(GRCh38):g.xxxxxxA>G format
    if re.match(r"Chr[0-9XY]+\([A-Za-z0-9]+\):g\.", raw, re.I):
        genomic = raw
        return (cdna, protein, genomic, notes_frag)

    # Parse colon-separated fields
    parts = [p.strip() for p in raw.split(":")]

    # Extract c. and p. notations from parts
    for part in parts:
        if part.startswith("c.") and not cdna:
            cdna = part
        elif part.startswith("p.") and not protein:
            protein = part
        elif part.startswith("g.") and not genomic:
            genomic = part

    # If nothing found, save raw to notes
    if not cdna and not protein and not genomic:
        notes_frag = raw

    return (cdna, protein, genomic, notes_frag)


def safe_str(val) -> str:
    """Convert value to string, returning '' for NaN/None."""
    if val is None:
        return ""
    s = str(val).strip()
    if s.lower() in ("nan", "none", "null", "na", "n/a"):
        return ""
    return s


def pipe_join(items: List[str]) -> str:
    """Join non-empty items with pipe separator."""
    return "|".join(x for x in items if x and x.strip())


# ---------------------------------------------------------------------------
# File-specific parsers → list of intermediate dicts
# ---------------------------------------------------------------------------

def _empty_row() -> dict:
    """Return empty intermediate dict."""
    return {
        "_source_id": "",
        "_raw_sex": "",
        "_raw_age": "",
        "_raw_consanguinity": "",
        "_raw_family_history": "",
        "_raw_test": "",
        "_raw_result": "",
        "_raw_inheritance": [],
        "_raw_phenotype": "",
        "_raw_variants": [],  # list of dicts: {gene, cdna, protein, genomic, zygosity, acmg_class, acmg_evidence, omim}
        "_raw_disease": "",
        "_raw_omim": "",
        "_pmid": [],
        "_notes": "",
        "_source_file": "",
    }


def parse_ahmed_variants(path: str, limit: Optional[int] = None) -> List[dict]:
    """Parse ahmed-variants.tsv."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    if limit:
        df = df.head(limit)

    rows = []
    # Variant slot column groups (with actual column names from the file)
    variant_slots = [
        {
            "gene": "Gene",
            "cdna": "Variant ",
            "protein": "Protein",
            "zygosity": "Zygosity",
            "pathogenicity": "pathogenicity ",
            "inheritance": "Inheritance",
            "omim": "Omim",
            "acmg": "ACMG Classification",
        },
        {
            "gene": "Gene 2",
            "cdna": "Variant 2  ",
            "protein": "Protein 2",
            "zygosity": "Zygosity 2",
            "pathogenicity": "pathogenicity 2 ",
            "inheritance": "Inheritance.1",
            "omim": "Omim.1",
            "acmg": "ACMG Classification.1",
        },
        {
            "gene": "Gene 3",
            "cdna": "Variant 3  ",
            "protein": "Protein 3",
            "zygosity": "Zygosity 3",
            "pathogenicity": "pathogenicity 3",
            "inheritance": "Inheritance.2",
            "omim": "OMIM",
            "acmg": "ACMG Classification.2",
        },
        {
            "gene": "Gene 4",
            "cdna": "Variant4",
            "protein": "Protein 4",
            "zygosity": "Zygosity 4",
            "pathogenicity": "pathogenicity 4",
            "inheritance": "Inheritance.3",
            "omim": "OMIM.1",
            "acmg": "ACMG Classification.3",
        },
    ]

    cols = set(df.columns)

    for _, row in df.iterrows():
        rec = _empty_row()
        rec["_source_file"] = "ahmed-variants"
        rec["_source_id"] = safe_str(row.get("Case", ""))
        rec["_raw_sex"] = safe_str(row.get("Gender", ""))
        # DOB is omitted per plan
        rec["_raw_consanguinity"] = safe_str(row.get("Consanguinity ", ""))
        rec["_raw_family_history"] = safe_str(row.get("Family History", ""))
        rec["_raw_test"] = safe_str(row.get(" Test Type", row.get("Test", "")))
        rec["_raw_result"] = safe_str(row.get("Results Internal", ""))
        rec["_raw_phenotype"] = safe_str(row.get("HPOs", ""))

        variants = []
        inheritance_vals = []
        for slot in variant_slots:
            gene_col = slot["gene"]
            if gene_col not in cols:
                continue
            gene = safe_str(row.get(gene_col, ""))
            if not gene:
                continue
            cdna = safe_str(row.get(slot["cdna"], ""))
            protein = safe_str(row.get(slot["protein"], ""))
            zyg = safe_str(row.get(slot["zygosity"], ""))
            path_val = safe_str(row.get(slot["pathogenicity"], ""))
            inh = safe_str(row.get(slot["inheritance"], ""))
            omim = safe_str(row.get(slot["omim"], ""))
            acmg = safe_str(row.get(slot["acmg"], ""))
            if inh:
                inheritance_vals.append(inh)
            variants.append({
                "gene": gene,
                "cdna": cdna,
                "protein": protein,
                "genomic": "",
                "zygosity": zyg,
                "acmg_class": path_val if path_val else acmg,
                "acmg_evidence": acmg,
                "omim": omim,
            })

        rec["_raw_variants"] = variants
        rec["_raw_inheritance"] = inheritance_vals
        rec["_raw_omim"] = safe_str(row.get("Omim", ""))
        rows.append(rec)

    return rows


def parse_ahmed_pmid(path: str, limit: Optional[int] = None) -> List[dict]:
    """Parse ahmed-pmid28454995.tsv."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    if limit:
        df = df.head(limit)

    # Skip PolyPhen-2, SIFT, Align GVGD, ESP, ExAC, dbSNP, frequency columns
    skip_patterns = ["polyphen", "sift", "align gvgd", "esp", "exac", "dbsnp",
                     "allele frequency", "maf", "unnamed"]

    rows = []
    variant_slots = [
        {
            "gene": "Variant 1 Gene",
            "transcript": "Transcript ID",
            "gdna": "g.DNA",
            "cdna": "c.DNA",
            "protein": "Amino Acid",
            "zygosity": "Variant 1 Allele State",
            "pmid_col": "References PMID",
        },
        {
            "gene": "Variant 2  Gene",
            "transcript": "Transcript ID.1",
            "gdna": "g.DNA.1",
            "cdna": " c.DNA",
            "protein": "Amino Acid.1",
            "zygosity": "Variant 2 Allele State",
            "pmid_col": "References",
        },
    ]

    cols = set(df.columns)

    for _, row in df.iterrows():
        rec = _empty_row()
        rec["_source_file"] = "ahmed-pmid28454995"
        rec["_source_id"] = safe_str(row.get("Case", ""))
        rec["_raw_sex"] = ""  # No sex column
        rec["_raw_consanguinity"] = ""
        rec["_raw_family_history"] = ""
        rec["_raw_test"] = "WES"  # All WES per paper
        rec["_raw_result"] = ""
        rec["_raw_phenotype"] = safe_str(row.get("Additional clinical phenotype", ""))
        rec["_raw_disease"] = safe_str(row.get("Diagnosis", ""))
        rec["_raw_omim"] = safe_str(row.get("OMIM", ""))

        variants = []
        pmids = []
        for slot in variant_slots:
            gene_col = slot["gene"]
            if gene_col not in cols:
                continue
            gene = safe_str(row.get(gene_col, ""))
            if not gene:
                continue
            gdna = safe_str(row.get(slot["gdna"], ""))
            cdna = safe_str(row.get(slot["cdna"], ""))
            protein = safe_str(row.get(slot["protein"], ""))
            zyg = safe_str(row.get(slot["zygosity"], ""))
            pmid_raw = safe_str(row.get(slot["pmid_col"], ""))
            pmids.extend(extract_pmid(pmid_raw))
            # gdna is GRCh38 format e.g. ChrX(GRCh38):g.78123223G>C
            variants.append({
                "gene": gene,
                "cdna": cdna,
                "protein": protein,
                "genomic": gdna,
                "zygosity": zyg,
                "acmg_class": "",
                "acmg_evidence": "",
                "omim": "",
            })

        rec["_raw_variants"] = variants
        rec["_pmid"] = list(dict.fromkeys(pmids))
        rows.append(rec)

    return rows


def parse_fawzan(path: str, limit: Optional[int] = None) -> List[dict]:
    """Parse fawzan-variants.tsv."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    if limit:
        df = df.head(limit)

    rows = []
    for _, row in df.iterrows():
        rec = _empty_row()
        rec["_source_file"] = "fawzan-variants"
        rec["_source_id"] = safe_str(row.get("ID", ""))
        rec["_raw_sex"] = safe_str(row.get("Gender", ""))
        rec["_raw_age"] = safe_str(row.get("Age", ""))
        rec["_raw_consanguinity"] = safe_str(row.get("Consanguinity", ""))
        rec["_raw_family_history"] = safe_str(row.get("Family hx", ""))
        rec["_raw_test"] = safe_str(row.get("Test", ""))
        raw_result = safe_str(row.get("Result", ""))
        rec["_raw_result"] = raw_result
        rec["_raw_phenotype"] = safe_str(row.get("Phenotype", ""))

        # Note HGMD accession in notes
        hgmd = safe_str(row.get("HGMD", ""))
        if hgmd:
            rec["_notes"] = f"HGMD:{hgmd}"

        variants = []
        zyg_raw = safe_str(row.get("Zyogsity", ""))  # note typo in header

        for var_col in ["Variant(s)", "Variant(s).1"]:
            var_raw = safe_str(row.get(var_col, ""))
            if not var_raw:
                continue

            # Parse GENE:NM_xxxxx:ExonN:c.xxx:p.xxx format
            cdna, protein, genomic, notes_frag = parse_variant_string(var_raw)
            gene = ""
            parts = [p.strip() for p in var_raw.split(":")]
            if parts and not parts[0].startswith(("c.", "p.", "g.", "NC_", "NM_")):
                gene = parts[0]

            if notes_frag and not cdna and not genomic:
                if rec["_notes"]:
                    rec["_notes"] += "; " + notes_frag
                else:
                    rec["_notes"] = notes_frag

            variants.append({
                "gene": gene,
                "cdna": cdna,
                "protein": protein,
                "genomic": genomic,
                "zygosity": zyg_raw,
                "acmg_class": "",
                "acmg_evidence": "",
                "omim": "",
            })

        rec["_raw_variants"] = variants
        rows.append(rec)

    return rows


def parse_marwa(path: str, limit: Optional[int] = None) -> List[dict]:
    """Parse marwa-variants.tsv."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    if limit:
        df = df.head(limit)

    rows = []
    for _, row in df.iterrows():
        rec = _empty_row()
        rec["_source_file"] = "marwa-variants"
        rec["_source_id"] = safe_str(row.get("id", ""))
        rec["_raw_sex"] = safe_str(row.get("patient_gender", ""))
        rec["_raw_age"] = safe_str(row.get("patient_age", ""))
        rec["_raw_consanguinity"] = safe_str(row.get("consanguinity", ""))
        rec["_raw_test"] = safe_str(row.get("test", ""))
        rec["_raw_result"] = safe_str(row.get("result", ""))
        rec["_raw_phenotype"] = safe_str(row.get("phenotypes", ""))

        # Extract PMID from reference column (PubMed URLs)
        ref = safe_str(row.get("reference", ""))
        rec["_pmid"] = extract_pmid(ref)

        # Parse variants (pipe-delimited GENE:NM_:c. format)
        variants_raw = safe_str(row.get("variants", ""))
        zygosity_raw = safe_str(row.get("zygosity", ""))
        zyg_parts = [z.strip() for z in zygosity_raw.split("|")] if zygosity_raw else []

        variants = []
        if variants_raw:
            var_parts = [v.strip() for v in variants_raw.split("|")]
            for i, vp in enumerate(var_parts):
                if not vp:
                    continue
                cdna, protein, genomic, notes_frag = parse_variant_string(vp)
                gene = ""
                if not notes_frag:
                    # Variant is parseable — safe to extract gene from first colon-token
                    parts = [p.strip() for p in vp.split(":")]
                    if parts and not parts[0].startswith(("c.", "p.", "g.", "NC_", "NM_")):
                        candidate = parts[0]
                        # Handle "GENE NM_xxx" format — take only first word
                        if " " in candidate:
                            first_tok = candidate.split()[0]
                            candidate = first_tok if re.match(r"^[A-Z][A-Z0-9-]+$", first_tok) else ""
                        gene = candidate
                zyg = zyg_parts[i] if i < len(zyg_parts) else ""
                if notes_frag and not cdna and not genomic:
                    if rec["_notes"]:
                        rec["_notes"] += "; " + notes_frag
                    else:
                        rec["_notes"] = notes_frag
                variants.append({
                    "gene": gene,
                    "cdna": cdna,
                    "protein": protein,
                    "genomic": genomic,
                    "zygosity": zyg,
                    "acmg_class": safe_str(row.get("pathogenicity", "")),
                    "acmg_evidence": "",
                    "omim": "",
                })

        rec["_raw_variants"] = variants
        rows.append(rec)

    return rows


def parse_pmc6562004(path: str, limit: Optional[int] = None) -> List[dict]:
    """Parse PMC6562004.tsv (has comment row as first line)."""
    # The first row starts with '#' — use skiprows=1
    df = pd.read_csv(path, sep="\t", skiprows=1, dtype=str)
    if limit:
        df = df.head(limit)

    rows = []
    for _, row in df.iterrows():
        rec = _empty_row()
        rec["_source_file"] = "PMC6562004"
        rec["_source_id"] = safe_str(row.get("ID", ""))
        rec["_raw_sex"] = safe_str(row.get("Gender", ""))
        rec["_raw_age"] = safe_str(row.get("Age", ""))
        rec["_raw_consanguinity"] = safe_str(row.get("Consanguinity", ""))
        rec["_raw_family_history"] = safe_str(row.get("Family hx", ""))
        rec["_raw_test"] = safe_str(row.get("Test", ""))
        solved = safe_str(row.get("Solved Y/N", ""))
        rec["_raw_result"] = "Y" if solved.upper() == "Y" else ("N" if solved.upper() == "N" else "")
        rec["_raw_phenotype"] = safe_str(row.get("Phenotype", ""))

        # Primary indication
        pi = safe_str(row.get("Primary indication", ""))
        if pi:
            pi_text = PRIMARY_INDICATION_MAP.get(pi.strip(), pi.strip())
            rec["_notes"] = f"Primary indication: {pi_text}"

        # Parse variant: Gene(s) and Variant(s) columns
        # Variant(s) format: NM_003242.5:c.1301T>A:p.Met434Lys
        gene_raw = safe_str(row.get("Gene(s)", ""))
        var_raw = safe_str(row.get("Variant(s)", ""))
        zyg_raw = safe_str(row.get("Zygosity", ""))
        inh_raw = safe_str(row.get("Inheritance AR/AD/XL/NA", ""))

        variants = []
        if gene_raw and var_raw:
            # May have pipe/semicolon separated for dual/triple locus
            genes = re.split(r"[;|]", gene_raw)
            vars_parts = re.split(r"[;|]", var_raw)
            zyg_parts = re.split(r"[;|]", zyg_raw) if zyg_raw else []
            for i, (g, v) in enumerate(zip(genes, vars_parts)):
                g = g.strip()
                v = v.strip()
                if not g and not v:
                    continue
                cdna, protein, genomic, notes_frag = parse_variant_string(v)
                zyg = zyg_parts[i].strip() if i < len(zyg_parts) else zyg_raw
                if notes_frag and not cdna and not genomic:
                    if rec["_notes"]:
                        rec["_notes"] += "; " + notes_frag
                    else:
                        rec["_notes"] = notes_frag
                variants.append({
                    "gene": g,
                    "cdna": cdna,
                    "protein": protein,
                    "genomic": genomic,
                    "zygosity": zyg,
                    "acmg_class": "",
                    "acmg_evidence": "",
                    "omim": "",
                })

        rec["_raw_variants"] = variants
        rec["_raw_inheritance"] = [inh_raw] if inh_raw else []
        rows.append(rec)

    return rows


def parse_pmc7082194(path: str, limit: Optional[int] = None) -> List[dict]:
    """Parse PMC7082194.tsv (14 cols, col13 = clinical phenotype free text)."""
    df = pd.read_csv(path, sep="\t", dtype=str, header=0)
    if limit:
        df = df.head(limit)

    rows = []
    for _, row in df.iterrows():
        rec = _empty_row()
        rec["_source_file"] = "PMC7082194"
        rec["_source_id"] = safe_str(row.get("ID", ""))
        rec["_raw_sex"] = safe_str(row.get("Sex", ""))
        rec["_raw_age"] = safe_str(row.get("Age (years, months)", ""))
        rec["_raw_phenotype"] = safe_str(row.get("Unnamed: 13", ""))  # col 13 = clinical phenotype

        # PMID from Ref column
        ref = safe_str(row.get("Ref ", row.get("Ref", "")))
        rec["_pmid"] = extract_pmid(ref)

        # Note that source data is hg19
        rec["_notes"] = "[source: hg19]"

        # Gene symbol from Gene(s) column
        gene_raw = safe_str(row.get("Gene(s)", ""))

        # Try to extract valid gene symbol
        gene = ""
        # If it looks like a chromosomal description, put in notes
        if re.search(r"[Kk][Bb]\s+duplication|deletion|chr\d", gene_raw, re.I):
            if rec["_notes"]:
                rec["_notes"] += f"; gene_raw: {gene_raw}"
            else:
                rec["_notes"] = f"gene_raw: {gene_raw}"
        elif re.match(r"^[A-Z][A-Z0-9]+\d*$", gene_raw):
            gene = gene_raw

        # Nucleotide Change (c.), Protein Change (p.), HGVS Format
        cdna_raw = safe_str(row.get("Nucleotide Change", ""))
        protein_raw = safe_str(row.get("Protein Change", ""))
        hgvs_raw = safe_str(row.get("HGVS Format", ""))

        # Parse c. notation from Nucleotide Change
        cdna = ""
        if cdna_raw:
            # May be multiple: "maternal: c.172C>T; paternal: c.590A>G"
            cdna_matches = re.findall(r"c\.[^\s;,]+", cdna_raw)
            cdna = ";".join(cdna_matches) if cdna_matches else cdna_raw

        # Parse p. notation from Protein Change
        protein = ""
        if protein_raw:
            prot_matches = re.findall(r"p\.[^\s;,]+", protein_raw)
            if prot_matches:
                # Convert 1-letter to 3-letter
                protein = ";".join(convert_protein_1letter_to_3letter(p) for p in prot_matches)
            else:
                protein = protein_raw

        # Inheritance: extract zygosity keywords
        inh_raw = safe_str(row.get("Inheritance (if known)", ""))
        inh_list = [inh_raw] if inh_raw else []

        variants = []
        if gene or cdna or protein:
            # Handle compound variants (multiple c. in cdna)
            cdna_parts = cdna.split(";") if cdna else [""]
            prot_parts = protein.split(";") if protein else [""]
            max_vars = max(len(cdna_parts), len(prot_parts))
            for i in range(max_vars):
                c = cdna_parts[i] if i < len(cdna_parts) else ""
                p = prot_parts[i] if i < len(prot_parts) else ""
                variants.append({
                    "gene": gene,
                    "cdna": c.strip(),
                    "protein": p.strip(),
                    "genomic": "",  # hg19, noted above
                    "zygosity": inh_raw,  # zygosity extracted in normalize_one_row
                    "acmg_class": "",
                    "acmg_evidence": "",
                    "omim": "",
                })

        rec["_raw_variants"] = variants
        rec["_raw_inheritance"] = inh_list
        rows.append(rec)

    return rows


def parse_ddd_diagnoses(path: str, limit: Optional[int] = None) -> List[dict]:
    """
    Parse ddd-diagnoses.tsv.
    Each row is a gene-disease-phenotype association, not a patient record.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    if limit:
        df = df.head(limit)

    rows = []
    for idx, row in df.iterrows():
        rec = _empty_row()
        rec["_source_file"] = "ddd-diagnoses"
        rec["_source_id"] = str(idx)
        # No patient fields
        rec["_raw_sex"] = ""
        rec["_raw_age"] = ""
        rec["_raw_consanguinity"] = ""
        rec["_raw_family_history"] = ""

        # Gene
        gene = safe_str(row.get("gene", ""))
        # HPO IDs: semicolon+space delimited, direct passthrough
        hpo_raw = safe_str(row.get("hpo_ids", ""))
        # Disease
        omim_raw = safe_str(row.get("omim_id", ""))
        syndrome = safe_str(row.get("syndrome", ""))
        # Inheritance: allelic_mode
        allelic_mode = safe_str(row.get("allelic_mode", ""))
        # PMID
        pubmed_raw = safe_str(row.get("pubmed_ids", ""))

        rec["_raw_phenotype"] = hpo_raw  # Direct HPO IDs
        rec["_raw_inheritance"] = [allelic_mode] if allelic_mode else []
        rec["_raw_omim"] = omim_raw
        rec["_raw_disease"] = syndrome
        rec["_pmid"] = extract_pmid(pubmed_raw)

        # No variants for DDD
        # Gene goes in variants list with no variant info for gene_symbol output
        rec["_raw_variants"] = [{
            "gene": gene,
            "cdna": "",
            "protein": "",
            "genomic": "",
            "zygosity": "",
            "acmg_class": "",
            "acmg_evidence": "",
            "omim": omim_raw,
        }]

        rows.append(rec)

    return rows


# ---------------------------------------------------------------------------
# Disease normalization
# ---------------------------------------------------------------------------

def normalize_disease(
    raw_disease: str,
    raw_omim: str,
    nu,
) -> Tuple[List[str], List[str], List[str], List[str], List[str]]:
    """
    Returns (mondo_ids, mondo_labels, omim_ids, omim_labels, orphanet_ids).
    """
    mondo_ids: List[str] = []
    mondo_labels: List[str] = []
    omim_ids: List[str] = []
    omim_labels: List[str] = []
    orphanet_ids: List[str] = []

    # Parse OMIM IDs from raw_omim
    if raw_omim:
        # May be "300653 & 609340" or "300653" etc.
        nums = re.findall(r"\d{5,6}", raw_omim)
        for n in nums:
            omim_ids.append(f"OMIM:{n}")
            # Get label from disease_map
            label = nu.disease_map.get(f"OMIM:{n}", nu.disease_map.get(n, ""))
            omim_labels.append(label)

    # Try MONDO lookup for disease name
    if raw_disease and nu.mondo_names:
        from rapidfuzz import process as rfprocess, fuzz
        results = rfprocess.extract(raw_disease, nu.mondo_names, limit=1, scorer=fuzz.token_sort_ratio)
        if results and results[0][1] >= 80:
            idx = results[0][2]
            mondo_ids.append(nu.mondo_list[idx]["id"])
            mondo_labels.append(nu.mondo_list[idx]["name"])

    return mondo_ids, mondo_labels, omim_ids, omim_labels, orphanet_ids


# ---------------------------------------------------------------------------
# Phenotype normalization for DDD (direct HPO passthrough)
# ---------------------------------------------------------------------------

def normalize_ddd_phenotypes(hpo_raw: str, hpo_label_map=None) -> Tuple[List[str], List[str], List[str], List[str], List[str], List[str]]:
    """
    Parse semicolon-delimited HPO IDs from ddd-diagnoses hpo_ids field.
    Returns (present_ids, present_labels, excluded_ids, excluded_labels, modifier_ids, modifier_labels).
    """
    present_ids: List[str] = []
    present_labels: List[str] = []
    excluded_ids: List[str] = []
    excluded_labels: List[str] = []
    modifier_ids: List[str] = []
    modifier_labels: List[str] = []

    if not hpo_raw:
        return present_ids, present_labels, excluded_ids, excluded_labels, modifier_ids, modifier_labels

    parts = re.split(r"[;|]", hpo_raw)
    for part in parts:
        part = part.strip()
        if re.match(r"^HP:\d{7}$", part):
            present_ids.append(part)
            present_labels.append(hpo_label_map.get(part, "") if hpo_label_map else "")
            modifier_ids.append("")
            modifier_labels.append("")

    return present_ids, present_labels, excluded_ids, excluded_labels, modifier_ids, modifier_labels


def normalize_marwa_phenotypes(
    phenotype_raw: str,
    matcher,
    hpo_label_map=None,
) -> Tuple[List[str], List[str], List[str], List[str], List[str], List[str]]:
    """
    Parse marwa phenotypes: split on |, passthrough HP:XXXXXXX, match free text.
    Returns (present_ids, present_labels, excluded_ids, excluded_labels, modifier_ids, modifier_labels).
    """
    present_ids: List[str] = []
    present_labels: List[str] = []
    excluded_ids: List[str] = []
    excluded_labels: List[str] = []
    modifier_ids: List[str] = []
    modifier_labels: List[str] = []

    if not phenotype_raw:
        return present_ids, present_labels, excluded_ids, excluded_labels, modifier_ids, modifier_labels

    # Extract HPO IDs from inline text like "macular focal atrophy (HP:0007401)"
    # Also split on |
    parts = [p.strip() for p in phenotype_raw.split("|")]
    direct_hpo = set()
    free_text_parts = []

    for part in parts:
        # Find all HPO IDs in the part
        hpo_matches = re.findall(r"HP:\d{7}", part)
        if hpo_matches:
            direct_hpo.update(hpo_matches)
        else:
            # Remove HPO ID annotations if any embedded
            clean = re.sub(r"\(HP:\d{7}\)", "", part).strip()
            if clean:
                free_text_parts.append(clean)

    # Add direct HPO IDs
    for hpo_id in sorted(direct_hpo):
        present_ids.append(hpo_id)
        present_labels.append(hpo_label_map.get(hpo_id, "") if hpo_label_map else "")
        modifier_ids.append("")
        modifier_labels.append("")

    # Match free text via phenotype matcher
    if free_text_parts and matcher:
        combined_text = ", ".join(free_text_parts)
        try:
            result = matcher.match(combined_text)
            for pm in result.phenotypes:
                if pm.excluded:
                    excluded_ids.append(pm.hpo_id)
                    excluded_labels.append(pm.label)
                else:
                    present_ids.append(pm.hpo_id)
                    present_labels.append(pm.label)
                    modifier_ids.append(pm.severity_id or "")
                    modifier_labels.append(pm.severity_label or "")
        except Exception:
            pass

    return present_ids, present_labels, excluded_ids, excluded_labels, modifier_ids, modifier_labels


# ---------------------------------------------------------------------------
# normalize_one_row: maps intermediate dict → final 37-column output dict
# ---------------------------------------------------------------------------

def normalize_one_row(
    rec: dict,
    matcher,
    nu,
    hgvs_validator: HgvsValidator,
    hpo_label_map=None,
) -> dict:
    """Normalize one intermediate record to the 37-column output schema."""
    out = {col: "" for col in OUTPUT_COLUMNS}

    out["source_file"] = rec.get("_source_file", "")
    out["source_id"] = rec.get("_source_id", "")
    out["notes"] = rec.get("_notes", "")

    # Sex
    sex_id, sex_label = nu.normalize_sex(rec.get("_raw_sex", ""))
    out["sex_id"] = sex_id
    out["sex_label"] = sex_label

    # Age
    out["age_iso"] = normalize_age(rec.get("_raw_age", ""))

    # Consanguinity
    cons_id, cons_label = nu.normalize_consanguinity(rec.get("_raw_consanguinity", ""))
    out["consanguinity_id"] = cons_id
    out["consanguinity_label"] = cons_label

    # Family history
    fh_id, fh_label = nu.normalize_family_history(rec.get("_raw_family_history", ""))
    out["family_history_id"] = fh_id
    out["family_history_label"] = fh_label

    # Test type
    out["test_type"] = normalize_test_type(rec.get("_raw_test", ""))

    # Result status
    out["result_status"] = normalize_result_status(rec.get("_raw_result", ""))

    # Inheritance
    inh_vals = rec.get("_raw_inheritance", [])
    if isinstance(inh_vals, str):
        inh_vals = [inh_vals]
    all_inh_ids: List[str] = []
    all_inh_labels: List[str] = []
    seen_inh: set = set()
    for inh_text in inh_vals:
        ids, labels = normalize_inheritance(inh_text)
        for i, l in zip(ids, labels):
            if i not in seen_inh:
                seen_inh.add(i)
                all_inh_ids.append(i)
                all_inh_labels.append(l)
    out["inheritance_id"] = pipe_join(all_inh_ids)
    out["inheritance_label"] = pipe_join(all_inh_labels)

    # Phenotypes
    raw_phenotype = rec.get("_raw_phenotype", "")
    out["phenotypes_raw"] = raw_phenotype

    source_file = rec.get("_source_file", "")

    present_ids: List[str] = []
    present_labels: List[str] = []
    excluded_ids: List[str] = []
    excluded_labels: List[str] = []
    modifier_ids: List[str] = []
    modifier_labels: List[str] = []

    if source_file == "ddd-diagnoses":
        # Direct HPO passthrough
        pi, pl, ei, el, mi, ml = normalize_ddd_phenotypes(raw_phenotype, hpo_label_map)
        present_ids, present_labels = pi, pl
        excluded_ids, excluded_labels = ei, el
        modifier_ids, modifier_labels = mi, ml
    elif source_file == "marwa-variants":
        # Mixed HPO IDs + free text
        pi, pl, ei, el, mi, ml = normalize_marwa_phenotypes(raw_phenotype, matcher, hpo_label_map)
        present_ids, present_labels = pi, pl
        excluded_ids, excluded_labels = ei, el
        modifier_ids, modifier_labels = mi, ml
    elif raw_phenotype and matcher:
        # General free text via phenotype matcher
        # For ahmed-variants, extract any HPO IDs inline first
        inline_hpo = re.findall(r"HP:\d{7}", raw_phenotype)
        # Remove inline HPO annotations for cleaner text
        clean_text = re.sub(r"\(HP:\d{7}\)", "", raw_phenotype).strip()
        clean_text = re.sub(r"HP:\d{7}", "", clean_text).strip()

        # Add direct HPO IDs
        for hpo_id in inline_hpo:
            if hpo_id not in present_ids:
                present_ids.append(hpo_id)
                present_labels.append(hpo_label_map.get(hpo_id, "") if hpo_label_map else "")
                modifier_ids.append("")
                modifier_labels.append("")

        # Match remaining text
        if clean_text:
            try:
                result = matcher.match(clean_text)
                for pm in result.phenotypes:
                    if pm.excluded:
                        excluded_ids.append(pm.hpo_id)
                        excluded_labels.append(pm.label)
                    else:
                        if pm.hpo_id not in present_ids:
                            present_ids.append(pm.hpo_id)
                            present_labels.append(pm.label)
                            modifier_ids.append(pm.severity_id or "")
                            modifier_labels.append(pm.severity_label or "")
            except Exception:
                pass

    out["phenotypes_present_ids"] = pipe_join(present_ids)
    out["phenotypes_present_labels"] = pipe_join(present_labels)
    out["phenotypes_excluded_ids"] = pipe_join(excluded_ids)
    out["phenotypes_excluded_labels"] = pipe_join(excluded_labels)
    out["phenotypes_modifier_ids"] = pipe_join(modifier_ids)
    out["phenotypes_modifier_labels"] = pipe_join(modifier_labels)

    # Variants
    raw_variants = rec.get("_raw_variants", [])
    gene_symbols: List[str] = []
    gene_entrez: List[str] = []
    variant_cdna: List[str] = []
    variant_protein: List[str] = []
    variant_genomic: List[str] = []
    variant_hgvs_valid: List[str] = []
    zygosity_ids: List[str] = []
    zygosity_labels: List[str] = []
    acmg_classifications: List[str] = []
    acmg_evidences: List[str] = []
    notes_parts: List[str] = [out["notes"]] if out["notes"] else []

    for vd in raw_variants:
        gene_raw = vd.get("gene", "")
        sym, entrez = normalize_gene(gene_raw, nu.genes)
        gene_symbols.append(sym)
        gene_entrez.append(entrez)

        # cDNA
        cdna = vd.get("cdna", "")
        if cdna:
            valid = hgvs_validator.validate(cdna)
            variant_cdna.append(cdna)
            variant_hgvs_valid.append(str(valid))
            if not valid:
                notes_parts.append(f"invalid_hgvs: {cdna}")
        else:
            variant_cdna.append("")
            variant_hgvs_valid.append("")

        # Protein (convert 1-letter if needed)
        protein = vd.get("protein", "")
        if protein:
            protein = convert_protein_1letter_to_3letter(protein)
        variant_protein.append(protein)

        # Genomic
        genomic = vd.get("genomic", "")
        variant_genomic.append(genomic)

        # Zygosity - extract from zygosity field or inheritance text
        zyg_raw = vd.get("zygosity", "")
        # Also check if inheritance text has zygosity keyword
        inh_zyg = extract_zygosity_from_inheritance(zyg_raw)
        if inh_zyg:
            zyg_id, zyg_label = inh_zyg
        else:
            zyg_id, zyg_label = normalize_zygosity(zyg_raw)
        zygosity_ids.append(zyg_id)
        zygosity_labels.append(zyg_label)

        # ACMG classification + evidence
        # acmg_ev_raw contains the full ACMG Classification string (class + evidence)
        # acmg_class_raw may be a shorter pathogenicity value (just class)
        acmg_class_raw = vd.get("acmg_class", "")
        acmg_ev_raw = vd.get("acmg_evidence", "")

        norm_class = ""
        ev_codes = ""
        if acmg_ev_raw:
            # Full ACMG string has both class and evidence codes
            norm_class, ev_codes = nu.normalize_variant_classification(acmg_ev_raw)
        if not norm_class and acmg_class_raw:
            # Fall back to pathogenicity-only field for class
            norm_class, _ = nu.normalize_variant_classification(acmg_class_raw)
        acmg_classifications.append(norm_class)
        acmg_evidences.append(ev_codes)

    # Compound heterozygous detection
    if zygosity_ids:
        variant_gene_zyg = list(zip(gene_symbols, zygosity_ids))
        zygosity_ids = infer_compound_het(variant_gene_zyg)
        # Update labels accordingly
        zygosity_labels = []
        for zyg_id in zygosity_ids:
            zyg_label_map = {
                "GENO:0000134": "hemizygous",
                "GENO:0000135": "heterozygous",
                "GENO:0000136": "homozygous",
                "GENO:0000137": "unspecified zygosity",
                "GENO:0000402": "compound heterozygous",
            }
            zygosity_labels.append(zyg_label_map.get(zyg_id, zyg_id))

    out["gene_symbol"] = pipe_join(gene_symbols)
    out["gene_entrez_id"] = pipe_join(gene_entrez)
    out["variant_cdna"] = pipe_join(variant_cdna)
    out["variant_protein"] = pipe_join(variant_protein)
    out["variant_genomic_grch38"] = pipe_join(variant_genomic)
    out["variant_hgvs_valid"] = pipe_join(variant_hgvs_valid)
    out["zygosity_id"] = pipe_join(zygosity_ids)
    out["zygosity_label"] = pipe_join(zygosity_labels)
    out["acmg_classification"] = pipe_join(acmg_classifications)
    out["acmg_evidence"] = pipe_join(acmg_evidences)

    # Disease normalization
    raw_disease = rec.get("_raw_disease", "")
    raw_omim = rec.get("_raw_omim", "")
    if not raw_omim:
        # Try OMIM from variants
        omim_from_variants = [vd.get("omim", "") for vd in raw_variants if vd.get("omim", "")]
        raw_omim = "|".join(omim_from_variants)

    mondo_ids, mondo_labels, omim_ids, omim_labels, orphanet_ids = normalize_disease(
        raw_disease, raw_omim, nu
    )
    out["disease_mondo_id"] = pipe_join(mondo_ids)
    out["disease_mondo_label"] = pipe_join(mondo_labels)
    out["disease_omim_ids"] = pipe_join(omim_ids)
    out["disease_omim_labels"] = pipe_join(omim_labels)
    out["disease_orphanet_ids"] = pipe_join(orphanet_ids)

    # PMID
    pmids = rec.get("_pmid", [])
    out["pmid"] = pipe_join(pmids)

    # Notes
    out["notes"] = "; ".join(p for p in notes_parts if p)

    return out


# ---------------------------------------------------------------------------
# ACMG OBO generation
# ---------------------------------------------------------------------------

def build_acmg_obo(out_path: str):
    """Generate ontology/acmg_criteria.obo with 27 ACMG/AMP 2015 criteria."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True) if os.path.dirname(out_path) else None
    with open(out_path, "w") as f:
        f.write("format-version: 1.2\n")
        f.write("data-version: acmg_amp_2015\n")
        f.write("ontology: acmg\n\n")
        f.write("[Term]\n")
        f.write("id: ACMG:0000000\n")
        f.write("name: ACMG evidence criterion\n")
        f.write('def: "An evidence criterion from the ACMG/AMP variant classification guidelines (Richards et al. 2015, PMID:25741868)." []\n\n')

        for code, definition, comment in ACMG_CRITERIA:
            f.write("[Term]\n")
            f.write(f"id: ACMG:{code}\n")
            f.write(f"name: {code}\n")
            f.write(f'def: "{definition}" []\n')
            f.write(f"comment: {comment}\n")
            f.write("is_a: ACMG:0000000\n\n")

    print(f"ACMG OBO written: {out_path} ({len(ACMG_CRITERIA)} criteria + 1 root = {len(ACMG_CRITERIA)+1} terms)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Combine and normalize 7 phenotype source files into a single TSV."
    )
    parser.add_argument(
        "--output", default="data/combined_normalized.tsv",
        help="Output TSV path (default: data/combined_normalized.tsv)"
    )
    parser.add_argument(
        "--acmg-obo", default="ontology/acmg_criteria.obo",
        help="Path for ACMG OBO output (default: ontology/acmg_criteria.obo)"
    )
    parser.add_argument(
        "--workers", type=int, default=4,
        help="Number of parallel workers (default: 4)"
    )
    parser.add_argument(
        "--limit", type=int, default=None,
        help="Process only first N rows per file (for testing)"
    )
    parser.add_argument(
        "--llm-model", default="deepseek",
        help="LLM model preset for PhenotypeMatcher (default: deepseek). "
             "Options: deepseek, gpt4oss, accurate, gemini"
    )
    args = parser.parse_args()

    # Initialize NormalizationUtils
    print("Initializing NormalizationUtils...")
    sys.path.insert(0, os.path.dirname(__file__))
    from normalization.normalization_utils import NormalizationUtils
    nu = NormalizationUtils()
    nu.initialize()

    # Initialize PhenotypeMatcher (shared, thread-safe)
    print("Initializing PhenotypeMatcher v2...")
    from phenotype_matcher_v2 import PhenotypeMatcher, MatcherConfig
    matcher = PhenotypeMatcher(MatcherConfig(
        hpo_path="ontology/hp.obo",
        mondo_path="ontology/mondo.obo",
        hpoa_path="data/phenotype.hpoa",
        llm_model=args.llm_model,
        match_cache_size=10_000,
    ))
    print(f"Using LLM model preset: {args.llm_model}")

    # Build HPO label map from loaded terms
    hpo_label_map = {t["id"]: t["name"] for t in nu.hpo_terms} if hasattr(nu, "hpo_terms") and nu.hpo_terms else {}

    # Initialize HGVS validator
    hgvs_validator = HgvsValidator()

    # Parse all 7 files
    print("Parsing source files...")
    phenotypes_dir = PHENOTYPES_DIR
    all_rows: List[dict] = []
    source_counts: dict = {}  # source → count for ID generation

    parsers = [
        ("ahmed-variants", "A", parse_ahmed_variants),
        ("ahmed-pmid28454995", "B", parse_ahmed_pmid),
        ("fawzan-variants", "F", parse_fawzan),
        ("marwa-variants", "M", parse_marwa),
        ("PMC6562004", "P", parse_pmc6562004),
        ("PMC7082194", "Q", parse_pmc7082194),
        ("ddd-diagnoses", "D", parse_ddd_diagnoses),
    ]

    for stem, letter, parse_fn in parsers:
        path = os.path.join(phenotypes_dir, f"{stem}.tsv")
        if not os.path.exists(path):
            print(f"  WARNING: {path} not found, skipping")
            continue
        print(f"  Parsing {stem}.tsv...")
        try:
            rows = parse_fn(path, limit=args.limit)
        except Exception as e:
            print(f"  ERROR parsing {stem}: {e}")
            traceback.print_exc()
            rows = []
        # Assign PAVS IDs
        for i, row in enumerate(rows, start=1):
            row["_pavs_id"] = f"PAVS:{letter}{i:07d}"
        source_counts[stem] = len(rows)
        all_rows.extend(rows)
        print(f"    → {len(rows)} rows")

    print(f"\nTotal intermediate rows: {len(all_rows)}")
    print("Normalizing rows (parallel)...")

    # Normalize all rows in parallel
    results: List[dict] = [None] * len(all_rows)

    def normalize_with_id(i_row):
        i, row = i_row
        try:
            out = normalize_one_row(row, matcher, nu, hgvs_validator, hpo_label_map)
            out["pavs_id"] = row.get("_pavs_id", "")
            return i, out
        except Exception as e:
            print(f"  ERROR in row {i}: {e}")
            traceback.print_exc()
            out = {col: "" for col in OUTPUT_COLUMNS}
            out["pavs_id"] = row.get("_pavs_id", "")
            out["source_file"] = row.get("_source_file", "")
            out["source_id"] = row.get("_source_id", "")
            out["notes"] = f"normalization_error: {e}"
            return i, out

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(normalize_with_id, (i, row)): i
                   for i, row in enumerate(all_rows)}
        for future in tqdm(as_completed(futures), total=len(all_rows), desc="Normalizing"):
            i, out = future.result()
            results[i] = out

    # Write output TSV
    print(f"\nWriting output to {args.output}...")
    os.makedirs(os.path.dirname(args.output), exist_ok=True) if os.path.dirname(args.output) else None
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in results:
            if row is not None:
                writer.writerow(row)

    print(f"Written {len(results)} rows ({len(results)+1} lines including header)")

    # Generate ACMG OBO
    print(f"\nGenerating ACMG OBO at {args.acmg_obo}...")
    build_acmg_obo(args.acmg_obo)

    # Print summary
    print("\nSummary:")
    for stem, letter, _ in parsers:
        count = source_counts.get(stem, 0)
        print(f"  {letter} ({stem}): {count} rows")
    print(f"  Total: {len(all_rows)} rows")
    print(f"\nDone. Output: {args.output}")


if __name__ == "__main__":
    main()
