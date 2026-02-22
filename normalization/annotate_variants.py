#!/usr/bin/env python3
"""Variant & Gene Annotation Enrichment for PAVS combined_normalized.tsv

Reads  data/combined_normalized.tsv  (7,510 rows, 38 columns)
Writes data/combined_annotated.tsv   (same rows, +17 variant + 10 gene columns)

Sources used (all local unless marked with [REST]):
  Variant-level
    [REST] Ensembl VEP — consequence, impact, HGVSc/p, SIFT, PolyPhen, gnomAD AFs, rsID, domains
    data/clinvar.vcf.gz        — clinical significance, disease, ClinVar allele ID
    data/erepo.tabbed.txt      — ClinGen expert curation (assertion, ACMG codes)
    data/allele-frequencies/   — Saudi population AF (AC/AN from 302 WGS individuals)
  Gene-level
    data/genes_to_disease.txt  — OMIM/disease associations
    data/gnomad.v4.1.constraint_metrics.tsv  — pLI, LOEUF (download with download_annotation_data.sh)
    data/HMD_HumanPhenotype.rpt + data/MGI_GenePheno.rpt + ontology/mp.obo  — mouse phenotypes
    data/goa_human.gaf.gz + ontology/go-basic.obo  — GO biological process / MF / CC
    data/E-GTEX-8/             — tissue expression (GTEx v8, percentile rank)

Usage:
    uv run python scripts/annotate_variants.py \\
        --input  data/combined_normalized.tsv \\
        --output data/combined_annotated.tsv

    # Quick smoke test (10 rows, no VEP REST):
    uv run python scripts/annotate_variants.py \\
        --input  data/combined_normalized.tsv \\
        --output /tmp/test_annotated.tsv \\
        --limit 10 --no-vep
"""

import argparse
import csv
import gzip
import json
import logging
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests
from urllib.parse import quote

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
REF = DATA / "reference"   # external reference databases
ONT = REPO / "ontology"
CACHE_DIR = DATA / "annotation_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# GRCh38 chromosome number → RefSeq NC accession
GRCH38_NC = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12",
    "4": "NC_000004.12", "5": "NC_000005.10", "6": "NC_000006.12",
    "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12",
    "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9",  "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
    "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
    "22": "NC_000022.11", "X": "NC_000023.11",  "Y": "NC_000024.10",
    "MT": "NC_012920.1",
}


def _build_hgvsg(chrom: str, pos: str, ref: str, alt: str) -> str:
    """Construct GRCh38 HGVS g. notation from VEP-resolved coordinates."""
    nc = GRCH38_NC.get(str(chrom).upper(), GRCH38_NC.get(str(chrom), ""))
    if not nc or not pos:
        return ""
    try:
        pos_i = int(pos)
    except (ValueError, TypeError):
        return ""
    ref = ref or ""
    alt = alt or ""
    if not ref and not alt:
        return ""
    # SNV
    if len(ref) == 1 and len(alt) == 1 and ref != "-" and alt != "-":
        return f"{nc}:g.{pos_i}{ref}>{alt}"
    # Deletion (alt represented as "-" or empty)
    if alt == "-" or (ref and not alt):
        if len(ref) == 1:
            return f"{nc}:g.{pos_i}del"
        return f"{nc}:g.{pos_i}_{pos_i + len(ref) - 1}del"
    # Insertion (ref represented as "-" or empty)
    if ref == "-" or (alt and not ref):
        return f"{nc}:g.{pos_i - 1}_{pos_i}ins{alt}"
    # Multi-base substitution / complex indel
    if len(ref) > 1 or len(alt) > 1:
        return f"{nc}:g.{pos_i}_{pos_i + len(ref) - 1}delins{alt}"
    return ""


# ---------------------------------------------------------------------------
# Chromosome normalisation helpers
# ---------------------------------------------------------------------------

_CHR_STRIP = re.compile(r"^[Cc]hr", re.I)


def norm_chrom(c: str) -> str:
    """Strip 'chr'/'Chr' prefix; handle NC_* accession → number."""
    c = _CHR_STRIP.sub("", c.strip())
    # NC_000001.11 → 1, NC_000023.11 → X, NC_000024.10 → Y
    if c.startswith("NC_"):
        nc_map = {f"{i:06d}": str(i) for i in range(1, 23)}
        nc_map["000023"] = "X"
        nc_map["000024"] = "Y"
        m = re.match(r"NC_(\d{6})\.\d+", c)
        if m:
            return nc_map.get(m.group(1), c)
    return c


# ---------------------------------------------------------------------------
# Stage 1: Extract unique variants and resolve coordinates
# ---------------------------------------------------------------------------

_GENOMIC_RE = re.compile(
    r"[Cc]hr([0-9XYMTxy]+)(?:\(GRCh38\))?:g\.(\d+)([A-Z]+)>([A-Z]+)",
    re.I,
)
_GENOMIC_DEL = re.compile(
    r"[Cc]hr([0-9XYMTxy]+)(?:\(GRCh38\))?:g\.(\d+)(?:_(\d+))?del([A-Z]*)",
    re.I,
)
_GENOMIC_INS = re.compile(
    r"[Cc]hr([0-9XYMTxy]+)(?:\(GRCh38\))?:g\.(\d+)_(\d+)ins([A-Z]+)",
    re.I,
)


def parse_genomic_hgvs(s: str):
    """Return (chrom_no_chr, pos_1based, ref, alt) or None."""
    m = _GENOMIC_RE.search(s)
    if m:
        return m.group(1).upper(), int(m.group(2)), m.group(3), m.group(4)
    m = _GENOMIC_DEL.search(s)
    if m:
        # Simple deletion — return as HGVS-style; VEP handles it
        return m.group(1).upper(), int(m.group(2)), None, None
    return None


def extract_variant_keys(df: pd.DataFrame):
    """Return a mapping {pavs_id: list_of_variant_dicts} and a set of unique HGVS keys
    suitable for Ensembl VEP REST batching.

    Each variant_dict has keys:
        'hgvs'  — query string sent to VEP (e.g. "SUMF1:c.785A>G")
        'coords' — (chrom, pos, ref, alt) if parsed from genomic column, else None
    """
    variant_map = {}  # pavs_id → list[dict]
    unique_hgvs = {}  # hgvs_string → (chrom, pos, ref, alt) or None

    for _, row in df.iterrows():
        pid = row["pavs_id"]
        gene = str(row.get("gene_symbol", "") or "").strip()
        cdna_raw = str(row.get("variant_cdna", "") or "").strip()
        genomic_raw = str(row.get("variant_genomic_grch38", "") or "").strip()

        variants = []

        # Split multi-variant rows (pipe-delimited)
        cdna_parts = [v.strip() for v in cdna_raw.split("|") if v.strip()] or [""]
        genomic_parts = [v.strip() for v in genomic_raw.split("|") if v.strip()] or [""] * len(cdna_parts)
        gene_parts = [g.strip() for g in gene.split("|")] if "|" in gene else [gene] * len(cdna_parts)

        for i, cdna in enumerate(cdna_parts):
            g = gene_parts[i] if i < len(gene_parts) else gene
            gcoord = genomic_parts[i] if i < len(genomic_parts) else ""

            coords = None
            hgvs = None

            # Try genomic coord first — use region endpoint key for SNVs with known coords
            if gcoord and gcoord not in ("nan", ""):
                parsed = parse_genomic_hgvs(gcoord)
                if parsed and parsed[2]:  # has chrom, pos, ref, alt
                    coords = parsed
                    chrom, pos, ref, alt = parsed
                    # Use region-style key so run_vep_rest picks the region endpoint.
                    # Format: __region__:11:77174825-77174825:1/T
                    hgvs = f"{_VEP_REGION_PREFIX}{chrom}:{pos}-{pos}:1/{alt}"

            # Fall back to GENE:c.xxx (HGVS endpoint)
            if hgvs is None and cdna and cdna not in ("nan", "") and g and g not in ("nan", ""):
                # Take first gene if compound symbol
                g0 = g.split("/")[0].strip()
                hgvs = f"{g0}:{cdna}"

            if hgvs:
                if hgvs not in unique_hgvs:
                    unique_hgvs[hgvs] = coords
                variants.append({"hgvs": hgvs, "coords": coords, "gene": g, "cdna": cdna})

        variant_map[pid] = variants

    return variant_map, unique_hgvs


# ---------------------------------------------------------------------------
# Stage 2: Ensembl VEP REST
# ---------------------------------------------------------------------------

VEP_GET_BASE = "https://rest.ensembl.org/vep/human/hgvs"
VEP_REGION_BASE = "https://rest.ensembl.org/vep/human/region"
VEP_GET_PARAMS = {
    "canonical": 1,
    "hgvs": 1,
    "domains": 1,
    "sift": "b",
    "polyphen": "b",
    "af": 1,
}
VEP_SLEEP = 0.1  # 10 req/second (Ensembl allows up to 15/s)

# Region-style key prefix used in unique_hgvs dict to distinguish from HGVS keys
_VEP_REGION_PREFIX = "__region__:"

MUTALYZER_BASE = "https://mutalyzer.nl/api/normalize"
MUTALYZER_SLEEP = 0.4  # ~2.5 req/sec


def run_vep_rest(unique_hgvs: dict, cache_file: Path) -> dict:
    """Query Ensembl VEP REST (GET per variant) for all unique HGVS strings.
    Returns {hgvs_string: annotation_dict}.
    Results are cached incrementally to cache_file — safe to interrupt and resume.
    """
    # Load existing cache
    cache = {}
    if cache_file.exists():
        log.info("Loading existing VEP cache from %s", cache_file)
        with open(cache_file) as f:
            cache = json.load(f)

    # Determine which variants still need querying
    to_query = [h for h in unique_hgvs if h and h not in cache]
    if not to_query:
        log.info("All %d variants already in VEP cache — skipping REST queries.", len(cache))
        return cache

    log.info("VEP GET: %d variants to query (%d already cached) ...",
             len(to_query), len(cache))

    session = requests.Session()
    session.headers.update({"Accept": "application/json"})

    new_since_save = 0
    for i, hgvs in enumerate(to_query, 1):
        # Region-style keys use the VEP region endpoint; all others use HGVS endpoint.
        if hgvs.startswith(_VEP_REGION_PREFIX):
            region_path = hgvs[len(_VEP_REGION_PREFIX):]  # e.g. "11:77174825-77174825:1/T"
            url = f"{VEP_REGION_BASE}/{quote(region_path, safe=':/-')}"
        else:
            url = f"{VEP_GET_BASE}/{quote(hgvs, safe='')}"
        for attempt in range(3):
            try:
                r = session.get(url, params=VEP_GET_PARAMS, timeout=30)
                if r.status_code == 200:
                    data = r.json()
                    entry = data[0] if isinstance(data, list) and data else {}
                    cache[hgvs] = _parse_vep_entry(entry) if entry else {}
                    break
                elif r.status_code == 400:
                    # HGVS not recognised — cache as empty so we skip on future runs
                    log.debug("VEP 400 (unrecognised HGVS): %s", hgvs)
                    cache[hgvs] = {}
                    break
                elif r.status_code == 429:
                    wait = int(r.headers.get("Retry-After", 30))
                    log.warning("VEP rate-limited — sleeping %ds", wait)
                    time.sleep(wait)
                    # retry
                else:
                    log.warning("VEP HTTP %s for %s: %s", r.status_code, hgvs, r.text[:200])
                    time.sleep(2 ** attempt)
            except requests.RequestException as exc:
                log.warning("VEP request error (attempt %d) for %s: %s", attempt + 1, hgvs, exc)
                time.sleep(2 ** attempt)

        new_since_save += 1

        # Incremental save every 100 new variants (safe to Ctrl-C and resume)
        if new_since_save >= 100:
            with open(cache_file, "w") as f:
                json.dump(cache, f)
            log.info("  [%d/%d] Saved VEP cache (%d total entries)", i, len(to_query), len(cache))
            new_since_save = 0

        time.sleep(VEP_SLEEP)

    # Final save
    with open(cache_file, "w") as f:
        json.dump(cache, f)
    filled = sum(1 for v in cache.values() if v)
    log.info("VEP REST complete: %d / %d variants annotated (cache: %d entries)",
             filled, len(to_query), len(cache))
    return cache


def _parse_vep_entry(e: dict) -> dict:
    """Extract the fields we care about from a single VEP response entry."""
    out = {
        "vep_consequence": "",
        "vep_impact": "",
        "vep_hgvsc": "",
        "vep_hgvsp": "",
        "vep_hgvsg": "",
        "vep_rsid": "",
        "vep_sift": "",
        "vep_polyphen": "",
        "vep_af_gnomad": "",
        "vep_af_gnomad_afr": "",
        "vep_af_gnomad_eas": "",
        "vep_af_gnomad_mid": "",
        "vep_af_gnomad_sas": "",
        "vep_domains": "",
        "vep_cadd_phred": "",
        "vep_revel": "",
        "vep_alphamissense": "",
        "vep_spliceai_max": "",
        # resolved coordinates
        "_chrom": "",
        "_pos": "",
        "_end": "",
        "_ref": "",
        "_alt": "",
    }

    # Genomic coordinates (GRCh38 from Ensembl REST)
    out["_chrom"] = str(e.get("seq_region_name", ""))
    out["_pos"] = str(e.get("start", ""))
    out["_end"] = str(e.get("end", ""))
    allele = e.get("allele_string", "")
    if "/" in allele:
        parts = allele.split("/", 1)
        out["_ref"] = parts[0]
        out["_alt"] = parts[1]

    # GRCh38 HGVS g. notation constructed from VEP coordinates
    out["vep_hgvsg"] = _build_hgvsg(out["_chrom"], out["_pos"], out["_ref"], out["_alt"])

    # Most severe consequence
    out["vep_consequence"] = e.get("most_severe_consequence", "")

    # Canonical transcript
    for tc in e.get("transcript_consequences", []):
        if tc.get("canonical") == 1:
            out["vep_impact"] = tc.get("impact", "")
            out["vep_hgvsc"] = tc.get("hgvsc", "")
            out["vep_hgvsp"] = tc.get("hgvsp", "")
            # SIFT: prediction(score)
            sift_pred = tc.get("sift_prediction", "")
            sift_score = tc.get("sift_score", "")
            if sift_pred:
                out["vep_sift"] = f"{sift_pred}({sift_score})" if sift_score != "" else sift_pred
            # PolyPhen
            pp_pred = tc.get("polyphen_prediction", "")
            pp_score = tc.get("polyphen_score", "")
            if pp_pred:
                out["vep_polyphen"] = f"{pp_pred}({pp_score})" if pp_score != "" else pp_pred
            # Domains: list of [db, domain] pairs
            domains_raw = tc.get("domains", [])
            dom_strs = [f"{d.get('db','')}/{d.get('name','')}" for d in domains_raw if d.get("db") in ("Pfam", "PROSITE", "PRINTS")]
            out["vep_domains"] = "|".join(dict.fromkeys(dom_strs))
            break  # canonical only

    # Colocated variants: rsID + gnomAD AFs
    for cv in e.get("colocated_variants", []):
        cv_id = cv.get("id", "")
        if cv_id.startswith("rs") and not out["vep_rsid"]:
            out["vep_rsid"] = cv_id
        # gnomAD AFs from colocated_variants frequencies
        freqs = cv.get("frequencies", {})
        for alt_key, freq_dict in freqs.items():
            if not out["vep_af_gnomad"] and "gnomade" in freq_dict:
                out["vep_af_gnomad"] = str(freq_dict.get("gnomade", ""))
            if not out["vep_af_gnomad_afr"] and "gnomade_afr" in freq_dict:
                out["vep_af_gnomad_afr"] = str(freq_dict["gnomade_afr"])
            if not out["vep_af_gnomad_eas"] and "gnomade_eas" in freq_dict:
                out["vep_af_gnomad_eas"] = str(freq_dict["gnomade_eas"])
            if not out["vep_af_gnomad_mid"] and "gnomade_mid" in freq_dict:
                out["vep_af_gnomad_mid"] = str(freq_dict["gnomade_mid"])
            if not out["vep_af_gnomad_sas"] and "gnomade_sas" in freq_dict:
                out["vep_af_gnomad_sas"] = str(freq_dict["gnomade_sas"])

    return out


# ---------------------------------------------------------------------------
# Stage 2b: Mutalyzer normalization for unrecognized HGVS strings
# ---------------------------------------------------------------------------

def run_mutalyzer(unrecognized: list, cache_file: Path) -> dict:
    """For each unrecognized HGVS string, query Mutalyzer to get NM_ transcript candidates.

    Returns {original_hgvs: [list_of_NM_candidates]}.
    Candidates are ordered: RefSeq Select first, then other NM_/NR_ transcripts.
    Empty list means Mutalyzer could not resolve the variant.
    Results are cached to cache_file.
    """
    cache = {}
    if cache_file.exists():
        log.info("Loading Mutalyzer cache from %s", cache_file)
        with open(cache_file) as f:
            cache = json.load(f)

    to_query = [h for h in unrecognized if h not in cache]
    if not to_query:
        log.info("All %d variants already in Mutalyzer cache", len(cache))
        return cache

    log.info("Mutalyzer: normalizing %d unrecognized variants ...", len(to_query))
    session = requests.Session()
    session.headers.update({"Accept": "application/json"})

    for i, hgvs in enumerate(to_query, 1):
        url = f"{MUTALYZER_BASE}/{quote(hgvs, safe='')}"
        try:
            r = session.get(url, timeout=30)
            data = r.json()
        except Exception as exc:
            log.warning("Mutalyzer error for %s: %s", hgvs, exc)
            cache[hgvs] = []
            time.sleep(MUTALYZER_SLEEP)
            continue

        # Successful normalization (no errors)
        normalized = data.get("normalized_description")
        if normalized and not data.get("custom", {}).get("errors"):
            cache[hgvs] = [normalized]
            time.sleep(MUTALYZER_SLEEP)
            continue

        # EGENEASREFERENCEID: gene name used instead of transcript — extract candidates
        custom = data.get("custom", {})
        errors = custom.get("errors", [])
        candidates = []
        for err in errors:
            if err.get("code") != "EGENEASREFERENCEID":
                continue
            options = err.get("options", {})
            # Extract variant suffix after first ":"
            parts = hgvs.split(":", 1)
            suffix = parts[1].lstrip(":") if len(parts) > 1 else ""
            seen = set()
            nm_select, nm_other = [], []
            for build_key, opts in options.items():
                for opt in opts:
                    tid = opt.get("transcript_id", "")
                    if tid in seen:
                        continue
                    seen.add(tid)
                    if not (tid.startswith("NM_") or tid.startswith("NR_")):
                        continue
                    candidate = f"{tid}:{suffix}"
                    if opt.get("tag") == "RefSeq Select":
                        nm_select.append(candidate)
                    else:
                        nm_other.append(candidate)
            candidates = nm_select + nm_other
            break

        cache[hgvs] = candidates
        time.sleep(MUTALYZER_SLEEP)

    with open(cache_file, "w") as f:
        json.dump(cache, f)
    resolved = sum(1 for v in cache.values() if v)
    log.info("Mutalyzer: %d / %d variants have transcript candidates", resolved, len(cache))
    return cache


# ---------------------------------------------------------------------------
# Stage 2c: VEP second pass using Mutalyzer transcript candidates
# ---------------------------------------------------------------------------

def run_vep_second_pass(mutalyzer_map: dict, vep_cache: dict,
                        vep_cache_file: Path, session: requests.Session) -> dict:
    """Try VEP GET with each Mutalyzer transcript candidate for previously-empty entries.

    Results are stored in vep_cache under the *original* HGVS key.
    Saves vep_cache_file every 50 new annotations.
    Returns updated vep_cache.
    """
    to_query = {
        orig: candidates
        for orig, candidates in mutalyzer_map.items()
        if candidates and not vep_cache.get(orig)
    }
    if not to_query:
        log.info("VEP second pass: nothing new to query")
        return vep_cache

    log.info("VEP second pass: %d variants with Mutalyzer candidates ...", len(to_query))
    new_since_save = 0

    for i, (orig_hgvs, candidates) in enumerate(to_query.items(), 1):
        found = False
        for candidate in candidates:
            url = f"{VEP_GET_BASE}/{quote(candidate, safe='')}"
            for attempt in range(3):
                try:
                    r = session.get(url, params=VEP_GET_PARAMS, timeout=30)
                    if r.status_code == 200:
                        data = r.json()
                        entry = data[0] if isinstance(data, list) and data else {}
                        if entry:
                            vep_cache[orig_hgvs] = _parse_vep_entry(entry)
                            log.debug("VEP 2nd pass: %s → %s", orig_hgvs, candidate)
                            found = True
                        break
                    elif r.status_code == 400:
                        break  # try next candidate
                    elif r.status_code == 429:
                        wait = int(r.headers.get("Retry-After", 30))
                        log.warning("VEP rate-limited — sleeping %ds", wait)
                        time.sleep(wait)
                    else:
                        log.warning("VEP HTTP %s for %s", r.status_code, candidate)
                        time.sleep(2 ** attempt)
                except requests.RequestException as exc:
                    log.warning("VEP error (attempt %d): %s", attempt + 1, exc)
                    time.sleep(2 ** attempt)
            if found:
                break
            time.sleep(VEP_SLEEP)

        if not found and orig_hgvs not in vep_cache:
            vep_cache[orig_hgvs] = {}  # exhausted all candidates

        new_since_save += 1
        if new_since_save >= 50:
            with open(vep_cache_file, "w") as f:
                json.dump(vep_cache, f)
            annotated = sum(1 for v in vep_cache.values() if v)
            log.info("  [%d/%d] VEP 2nd pass saved (total annotated: %d)",
                     i, len(to_query), annotated)
            new_since_save = 0
        time.sleep(VEP_SLEEP)

    with open(vep_cache_file, "w") as f:
        json.dump(vep_cache, f)
    annotated = sum(1 for v in vep_cache.values() if v)
    log.info("VEP second pass complete: %d / %d total variants annotated",
             annotated, len(vep_cache))
    return vep_cache


# ---------------------------------------------------------------------------
# Stage 3: ClinVar VCF
# ---------------------------------------------------------------------------

def load_clinvar(vcf_path: Path, target_positions: set) -> tuple:
    """Scan clinvar.vcf.gz; return (clinvar_index, clinvar_varid_index).

    clinvar_index:   {(chrom_no_chr, pos, ref, alt) → annotation_dict}
    clinvar_varid_index: {variation_id_str → (chrom, pos, ref, alt)}
    """
    log.info("Scanning ClinVar VCF for %d target positions ...", len(target_positions))
    clinvar_index = {}
    clinvar_varid_index = {}
    if not target_positions:
        return clinvar_index, clinvar_varid_index

    def _parse_info(info_str: str) -> dict:
        d = {}
        for item in info_str.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                d[k] = v
        return d

    opener = gzip.open if str(vcf_path).endswith(".gz") else open

    # If we have a tabix index, use pysam for targeted fetch
    try:
        import pysam
        tbi = Path(str(vcf_path) + ".tbi")
        if tbi.exists():
            log.info("  Using tabix index for fast ClinVar lookup ...")
            vcf = pysam.VariantFile(str(vcf_path))
            for chrom, pos in target_positions:
                try:
                    for rec in vcf.fetch(contig=chrom, start=pos - 1, stop=pos):
                        for alt in (rec.alts or []):
                            key = (rec.chrom.lstrip("chr"), rec.pos, rec.ref, alt)
                            varid = str(rec.info.get("ALLELEID", ""))
                            clinvar_index[key] = _build_clinvar_entry(rec, alt)
                            if varid:
                                clinvar_varid_index[varid] = key
                except ValueError:
                    pass
            return clinvar_index, clinvar_varid_index
    except Exception:
        pass

    # Fallback: scan whole file (filtered by target chrom:pos)
    target_pos_set = {(c, p) for c, p in target_positions}
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if len(parts) < 8:
                continue
            chrom = parts[0].lstrip("chr")
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if (chrom, pos) not in target_pos_set:
                continue
            ref = parts[3]
            alts = parts[4].split(",")
            info = _parse_info(parts[7])
            for alt in alts:
                key = (chrom, pos, ref, alt)
                varid = info.get("ALLELEID", "")
                entry = {
                    "clinvar_allele_id": varid,
                    "clinvar_sig": info.get("CLNSIG", "").replace("_", " "),
                    "clinvar_review_status": info.get("CLNREVSTAT", "").replace("_", " "),
                    "clinvar_disease": info.get("CLNDN", "").replace("_", " "),
                    "clinvar_disease_db": info.get("CLNDISDB", ""),
                }
                clinvar_index[key] = entry
                if varid:
                    clinvar_varid_index[varid] = key

    log.info("  ClinVar: %d entries loaded", len(clinvar_index))
    return clinvar_index, clinvar_varid_index


def _build_clinvar_entry(rec, alt: str) -> dict:
    """Build ClinVar annotation dict from pysam VariantRecord."""
    def _get(key):
        v = rec.info.get(key)
        if isinstance(v, (list, tuple)):
            return "|".join(str(x) for x in v)
        return str(v) if v is not None else ""

    return {
        "clinvar_allele_id": _get("ALLELEID"),
        "clinvar_sig": _get("CLNSIG").replace("_", " "),
        "clinvar_review_status": _get("CLNREVSTAT").replace("_", " "),
        "clinvar_disease": _get("CLNDN").replace("_", " "),
        "clinvar_disease_db": _get("CLNDISDB"),
    }


# ---------------------------------------------------------------------------
# Stage 4: ClinGen Evidence Repository
# ---------------------------------------------------------------------------

def load_clingen(erepo_path: Path) -> dict:
    """Return {clinvar_variation_id: annotation_dict}."""
    log.info("Loading ClinGen evidence repository ...")
    index = {}
    with open(erepo_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            varid = row.get("ClinVar Variation Id", "").strip()
            if varid:
                index[varid] = {
                    "clingen_assertion": row.get("Assertion", ""),
                    "clingen_acmg_codes": row.get("Applied Evidence Codes (Met)", ""),
                    "clingen_disease": row.get("Disease", ""),
                    "clingen_mondo_id": row.get("Mondo Id", ""),
                }
    log.info("  ClinGen: %d entries", len(index))
    return index


# ---------------------------------------------------------------------------
# Stage 5: Saudi allele frequencies
# ---------------------------------------------------------------------------

def load_saudi_af(af_dir: Path, target_positions: set) -> dict:
    """Scan Saudi AF files; return {(chrom_no_chr, pos, ref, alt) → (AC, AN)}.

    Saudi files have 6 columns: chrom(chr-prefixed) pos ref alt AN AC
    target_positions: set of (chrom_no_chr, pos) tuples
    """
    index = {}
    files = sorted(af_dir.glob("*Saudi*Allele*Freq*.txt"))
    if not files:
        log.warning("No Saudi allele frequency files found in %s", af_dir)
        return index

    if not target_positions:
        return index

    # Build fast lookup set with and without chr prefix
    pos_set = {(f"chr{c}", p) for c, p in target_positions} | {(c, p) for c, p in target_positions}

    log.info("Scanning %d Saudi AF file(s) for %d target positions ...", len(files), len(target_positions))
    for fpath in files:
        matches = 0
        with open(fpath) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                raw_chrom, pos_s, ref, alt = parts[0], parts[1], parts[2], parts[3]
                try:
                    pos = int(pos_s)
                except ValueError:
                    continue
                if (raw_chrom, pos) not in pos_set:
                    continue
                try:
                    an, ac = int(parts[4]), int(parts[5])
                except ValueError:
                    continue
                chrom = norm_chrom(raw_chrom)
                key = (chrom, pos, ref, alt)
                index[key] = (ac, an)
                matches += 1
        log.info("  %s: %d matches", fpath.name, matches)

    log.info("  Saudi AF total: %d entries", len(index))
    return index


# ---------------------------------------------------------------------------
# Stage 6A: Gene–disease (OMIM)
# ---------------------------------------------------------------------------

def load_gene_disease(path: Path) -> dict:
    """Return {gene_symbol: (omim_ids_pipe, types_pipe)}."""
    index = defaultdict(lambda: ([], []))
    if not path.exists():
        log.warning("genes_to_disease.txt not found: %s", path)
        return {}
    with open(path) as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            gene = parts[1].strip()
            atype = parts[2].strip()
            disease_id = parts[3].strip()  # e.g. OMIM:123456
            if atype == "MENDELIAN" and disease_id.startswith("OMIM:"):
                omim = disease_id
                lst_ids, lst_types = index[gene]
                if omim not in lst_ids:
                    lst_ids.append(omim)
                    lst_types.append(atype)
    return {g: ("|".join(ids), "|".join(types)) for g, (ids, types) in index.items()}


# ---------------------------------------------------------------------------
# Stage 6B: gnomAD constraint (pLI, LOEUF)
# ---------------------------------------------------------------------------

def load_gnomad_constraint(path: Path) -> dict:
    """Return {gene_symbol: {gnomad_pli, gnomad_loeuf}}."""
    if not path.exists():
        log.warning("gnomAD constraint file not found: %s — run download_annotation_data.sh", path)
        return {}
    log.info("Loading gnomAD constraint scores ...")
    index = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row.get("gene", "").strip()
            if not gene:
                continue
            # Only keep canonical gene entries (not transcript-level)
            if row.get("mane_select", "").strip() not in ("", "true", "1") and \
               row.get("canonical", "").strip() not in ("", "true", "1"):
                continue
            pli = row.get("lof.pLI", row.get("pLI", "")).strip()
            loeuf = row.get("lof.oe_ci.upper", row.get("oe_lof_upper", "")).strip()
            if gene not in index:
                index[gene] = {"gnomad_pli": pli, "gnomad_loeuf": loeuf}
    log.info("  gnomAD constraint: %d genes", len(index))
    return index


# ---------------------------------------------------------------------------
# Stage 6C: Mouse phenotypes (MGI)
# ---------------------------------------------------------------------------

def load_mouse_phenotypes(hmd_path: Path, gp_path: Path, mp_obo_path: Path) -> dict:
    """Return {human_gene_symbol: (mp_ids_pipe, mp_labels_pipe)}."""
    # Load MP labels
    mp_labels = {}
    if mp_obo_path.exists():
        log.info("Loading MP ontology ...")
        try:
            import pronto
            ont = pronto.Ontology(str(mp_obo_path))
            for term in ont.terms():
                mp_labels[term.id] = term.name
        except Exception as exc:
            log.warning("Could not load mp.obo: %s", exc)
    else:
        log.warning("mp.obo not found: %s", mp_obo_path)

    if not hmd_path.exists():
        log.warning("HMD_HumanPhenotype.rpt not found — run download_annotation_data.sh")
        return {}
    if not gp_path.exists():
        log.warning("MGI_GenePheno.rpt not found: %s", gp_path)
        return {}

    # Build human gene → set of mouse MGI gene IDs from HMD_HumanPhenotype.rpt
    # Columns (tab-sep): human_symbol | human_entrez | mouse_symbol | MGI_marker_id | human_chr | mouse_chr
    log.info("Loading HMD human–mouse ortholog map ...")
    human_to_mgi = defaultdict(set)
    with open(hmd_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            human_sym = parts[0].strip()
            mgi_id = parts[3].strip()
            if human_sym and mgi_id.startswith("MGI:"):
                human_to_mgi[human_sym].add(mgi_id)

    # Build MGI gene ID → set of MP terms from MGI_GenePheno.rpt
    # Columns: allele_composition | allele_symbol | genotype_MGI_id | background | MP_id | pubmed_id | gene_MGI_id | genotype_MGI_id2
    log.info("Loading MGI_GenePheno.rpt ...")
    mgi_to_mp = defaultdict(set)
    with open(gp_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            mp_id = parts[4].strip()
            gene_mgi = parts[6].strip()
            if mp_id.startswith("MP:") and gene_mgi.startswith("MGI:"):
                mgi_to_mp[gene_mgi].add(mp_id)

    # Combine
    result = {}
    for human_sym, mgi_ids in human_to_mgi.items():
        mp_ids = set()
        for mgi_id in mgi_ids:
            mp_ids.update(mgi_to_mp.get(mgi_id, set()))
        if mp_ids:
            labels = [mp_labels.get(mid, mid) for mid in sorted(mp_ids)]
            result[human_sym] = ("|".join(sorted(mp_ids)), "|".join(labels))

    log.info("  Mouse phenotype map: %d human genes with MP annotations", len(result))
    return result


# ---------------------------------------------------------------------------
# Stage 6D: GO annotations
# ---------------------------------------------------------------------------

def load_go_annotations(gaf_path: Path, go_obo_path: Path) -> dict:
    """Return {gene_symbol: {'P': [...], 'F': [...], 'C': [...]}} with term labels."""
    go_labels = {}
    if go_obo_path.exists():
        log.info("Loading GO ontology from %s ...", go_obo_path)
        try:
            import pronto
            ont = pronto.Ontology(str(go_obo_path))
            for term in ont.terms():
                go_labels[term.id] = term.name
        except Exception as exc:
            log.warning("Could not load go-basic.obo: %s", exc)
    else:
        log.warning("go-basic.obo not found at %s — run download_annotation_data.sh", go_obo_path)

    if not gaf_path.exists():
        log.warning("goa_human.gaf.gz not found: %s", gaf_path)
        return {}

    log.info("Loading GO annotations from %s ...", gaf_path)
    gene_go = defaultdict(lambda: {"P": set(), "F": set(), "C": set()})

    with gzip.open(gaf_path, "rt") as f:
        for line in f:
            if line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            gene_sym = parts[2].strip()
            go_id = parts[4].strip()
            aspect = parts[8].strip()  # P/F/C
            qualifier = parts[3].strip()
            if "NOT" in qualifier:
                continue
            if aspect in ("P", "F", "C"):
                gene_go[gene_sym][aspect].add(go_id)

    result = {}
    for gene, aspects in gene_go.items():
        entry = {}
        for asp, ids in aspects.items():
            values = []
            for goid in ids:
                label = go_labels.get(goid)
                values.append(label if label else goid)  # fall back to raw ID
            entry[asp] = values
        result[gene] = entry

    log.info("  GO annotations: %d genes", len(result))
    return result


# ---------------------------------------------------------------------------
# Stage 6E: GTEx expression
# ---------------------------------------------------------------------------

def load_gtex_expression(data_dir: Path) -> dict:
    """Return {gene_name: 'tissue1:pct|tissue2:pct|...'} for tissues at >=75th percentile.

    Uses percentile rank within each gene across 53 tissues.
    Expression value = median of the 5-value distribution per cell.
    """
    tpms_path = data_dir / "E-GTEX-8" / "E-GTEX-8-tpms.tsv"
    xml_path = data_dir / "E-GTEX-8" / "E-GTEX-8-configuration.xml"

    if not tpms_path.exists() or not xml_path.exists():
        log.warning("GTEx files not found in %s", data_dir / "E-GTEX-8")
        return {}

    # Parse configuration XML: g1..g53 → tissue label
    tree = ET.parse(str(xml_path))
    col_to_tissue = {}
    for ag in tree.getroot().findall(".//assay_group"):
        gid = ag.get("id", "")
        label = ag.get("label", "").replace(" ", "_")
        if gid:
            col_to_tissue[gid] = label

    log.info("Loading GTEx expression (%d tissues) ...", len(col_to_tissue))

    result = {}
    with open(tpms_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        # Map column index → tissue label
        col_idx_to_tissue = {}
        for idx, col_name in enumerate(header):
            if col_name in col_to_tissue:
                col_idx_to_tissue[idx] = col_to_tissue[col_name]

        tissue_cols = sorted(col_idx_to_tissue.keys())
        n_tissues = len(tissue_cols)
        if n_tissues == 0:
            log.warning("No tissue columns matched in GTEx tpms header")
            return {}

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            gene_name = parts[1].strip()
            if not gene_name:
                continue

            # Extract median (index 2 of 5) for each tissue
            medians = []
            for idx in tissue_cols:
                if idx >= len(parts):
                    medians.append(0.0)
                    continue
                vals_str = parts[idx].strip()
                if not vals_str:
                    medians.append(0.0)
                    continue
                vals = vals_str.split(",")
                try:
                    medians.append(float(vals[2]) if len(vals) >= 3 else float(vals[0]))
                except (ValueError, IndexError):
                    medians.append(0.0)

            # Compute percentile rank within this gene across all 53 tissues
            sorted_vals = sorted(medians)
            n = len(sorted_vals)
            if n == 0 or max(medians) == 0.0:
                continue

            pcts = []
            for i, idx in enumerate(tissue_cols):
                m = medians[i]
                rank = sum(1 for v in sorted_vals if v <= m)
                pct = round(rank / n * 100)
                if pct >= 75:
                    tissue = col_idx_to_tissue[idx]
                    pcts.append((pct, tissue))

            if pcts:
                pcts.sort(key=lambda x: -x[0])
                result[gene_name] = "|".join(f"{t}:{p}" for p, t in pcts)

    log.info("  GTEx: expression data for %d genes", len(result))
    return result


# ---------------------------------------------------------------------------
# Stage 7: Join and write
# ---------------------------------------------------------------------------

NEW_VARIANT_COLS = [
    "vep_consequence", "vep_impact", "vep_hgvsc", "vep_hgvsp", "vep_hgvsg",
    "vep_rsid", "vep_sift", "vep_polyphen",
    "vep_af_gnomad", "vep_af_gnomad_afr", "vep_af_gnomad_eas",
    "vep_af_gnomad_mid", "vep_af_gnomad_sas", "vep_domains",
    "vep_cadd_phred", "vep_revel", "vep_alphamissense", "vep_spliceai_max",
    "clinvar_allele_id", "clinvar_sig", "clinvar_review_status",
    "clinvar_disease", "clinvar_disease_db",
    "clingen_assertion", "clingen_acmg_codes",
    "clingen_disease", "clingen_mondo_id",
    "af_saudi", "ac_saudi",
]

NEW_GENE_COLS = [
    "gene_disease_omim", "gene_disease_type",
    "gnomad_pli", "gnomad_loeuf",
    "mouse_mp_ids", "mouse_mp_labels",
    "go_biological_process", "go_molecular_function", "go_cellular_component",
    "expressed_in",
]

EMPTY_VARIANT = {c: "" for c in NEW_VARIANT_COLS}
EMPTY_GENE = {c: "" for c in NEW_GENE_COLS}


def annotate_row(row, variant_map, vep_results, clinvar_index, clinvar_varid_index,
                 clingen_index, saudi_index, gene_disease, gnomad_constraint,
                 mouse_mp, go_anns, gtex_expr) -> dict:
    """Return new annotation columns for a single dataframe row."""
    out = {}
    out.update(EMPTY_VARIANT)
    out.update(EMPTY_GENE)

    pid = row["pavs_id"]
    gene = str(row.get("gene_symbol", "") or "").split("|")[0].strip()

    variants = variant_map.get(pid, [])

    # --- Variant-level --- annotate the first variant with known coords
    for v in variants:
        hgvs = v.get("hgvs", "")
        vep = vep_results.get(hgvs, {})
        chrom = vep.get("_chrom") or (v["coords"][0] if v.get("coords") else "")
        pos_s = vep.get("_pos") or (str(v["coords"][1]) if v.get("coords") else "")
        ref = vep.get("_ref") or (v["coords"][2] if v.get("coords") and v["coords"][2] else "")
        alt = vep.get("_alt") or (v["coords"][3] if v.get("coords") and v["coords"][3] else "")

        # Copy VEP fields
        for col in NEW_VARIANT_COLS:
            if col in vep and not out[col]:
                out[col] = str(vep.get(col, ""))

        # Fallback: construct vep_hgvsg from coordinates if not already set
        # (needed for existing cache entries that pre-date this field)
        if not out["vep_hgvsg"] and chrom and pos_s:
            out["vep_hgvsg"] = _build_hgvsg(chrom, pos_s, ref, alt)

        # ClinVar lookup by (chrom, pos, ref, alt)
        if chrom and pos_s and ref and alt:
            try:
                pos_i = int(pos_s)
            except ValueError:
                pos_i = 0
            cvar_key = (norm_chrom(chrom), pos_i, ref, alt)
            cv_entry = clinvar_index.get(cvar_key, {})
            for col in ["clinvar_allele_id", "clinvar_sig", "clinvar_review_status",
                        "clinvar_disease", "clinvar_disease_db"]:
                if not out[col]:
                    out[col] = cv_entry.get(col, "")

            # Saudi AF
            saudi_key = (norm_chrom(chrom), pos_i, ref, alt)
            saudi = saudi_index.get(saudi_key)
            if saudi and not out["af_saudi"]:
                ac, an = saudi
                out["ac_saudi"] = str(ac)
                out["af_saudi"] = f"{ac/an:.6f}" if an else ""

        # ClinGen via ClinVar variation ID
        cv_varid = out.get("clinvar_allele_id", "")
        if cv_varid:
            # clinvar_varid_index is keyed by ALLELEID (same field)
            cg = clingen_index.get(cv_varid, {})
            for col in ["clingen_assertion", "clingen_acmg_codes",
                        "clingen_disease", "clingen_mondo_id"]:
                if not out[col]:
                    out[col] = cg.get(col, "")

        break  # use annotations from first variant only (multi-gene rows handled above)

    # --- Gene-level ---
    first_gene = gene.split("/")[0].strip()

    gd = gene_disease.get(first_gene, ("", ""))
    out["gene_disease_omim"] = gd[0]
    out["gene_disease_type"] = gd[1]

    gc = gnomad_constraint.get(first_gene, {})
    out["gnomad_pli"] = gc.get("gnomad_pli", "")
    out["gnomad_loeuf"] = gc.get("gnomad_loeuf", "")

    mp = mouse_mp.get(first_gene)
    if mp:
        out["mouse_mp_ids"] = mp[0]
        out["mouse_mp_labels"] = mp[1]

    go = go_anns.get(first_gene, {})
    out["go_biological_process"] = "|".join(go.get("P", [])[:20])  # cap at 20 terms
    out["go_molecular_function"] = "|".join(go.get("F", [])[:20])
    out["go_cellular_component"] = "|".join(go.get("C", [])[:20])

    out["expressed_in"] = gtex_expr.get(first_gene, "")

    return out


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", default=str(DATA / "combined_normalized.tsv"),
                   help="Input TSV (default: data/combined_normalized.tsv)")
    p.add_argument("--output", default=str(DATA / "combined_annotated.tsv"),
                   help="Output TSV (default: data/combined_annotated.tsv)")
    p.add_argument("--limit", type=int, default=0,
                   help="Process only first N rows (0 = all)")
    p.add_argument("--no-vep", action="store_true",
                   help="Skip Ensembl VEP REST queries (use cached results only)")
    p.add_argument("--no-cache", action="store_true",
                   help="Ignore existing VEP cache and re-query")
    return p.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    vep_cache = CACHE_DIR / "ensembl_vep_rest.json"

    log.info("Reading %s ...", input_path)
    df = pd.read_csv(input_path, sep="\t", dtype=str, keep_default_na=False)
    if args.limit:
        df = df.head(args.limit)
    log.info("  %d rows, %d columns", len(df), len(df.columns))

    # ------ Stage 1: extract variant identifiers ------
    variant_map, unique_hgvs = extract_variant_keys(df)
    log.info("Unique variant HGVS strings: %d", len(unique_hgvs))

    # ------ Stage 2: VEP REST ------
    if args.no_cache and vep_cache.exists():
        vep_cache.unlink()

    mutalyzer_cache = CACHE_DIR / "mutalyzer_normalize.json"

    if args.no_vep:
        vep_results = {}
        if vep_cache.exists():
            log.info("Loading existing VEP cache (--no-vep skips new queries).")
            with open(vep_cache) as f:
                vep_results = json.load(f)
    else:
        vep_results = run_vep_rest(unique_hgvs, vep_cache)

        # ------ Stage 2b: Mutalyzer normalization for unrecognized variants ------
        unrecognized = [h for h, v in vep_results.items() if not v]
        log.info("VEP unrecognized: %d variants → Mutalyzer normalization ...", len(unrecognized))
        mutalyzer_map = run_mutalyzer(unrecognized, mutalyzer_cache)

        # ------ Stage 2c: VEP second pass with Mutalyzer candidates ------
        vep_session = requests.Session()
        vep_session.headers.update({"Accept": "application/json"})
        vep_results = run_vep_second_pass(mutalyzer_map, vep_results, vep_cache, vep_session)

    # Collect resolved coordinates from VEP results
    coord_set = set()
    for hgvs_str, anns in vep_results.items():
        chrom = anns.get("_chrom", "")
        pos_s = anns.get("_pos", "")
        if chrom and pos_s:
            try:
                coord_set.add((norm_chrom(chrom), int(pos_s)))
            except ValueError:
                pass
    # Also add coords parsed directly from genomic column
    for pid, variants in variant_map.items():
        for v in variants:
            if v.get("coords") and v["coords"][2]:
                coord_set.add((v["coords"][0], v["coords"][1]))

    log.info("Resolved genomic positions: %d", len(coord_set))

    # ------ Stage 3: ClinVar ------
    clinvar_path = REF / "clinvar.vcf.gz"
    if clinvar_path.exists():
        clinvar_index, clinvar_varid_index = load_clinvar(clinvar_path, coord_set)
    else:
        log.warning("clinvar.vcf.gz not found — skipping ClinVar annotations")
        clinvar_index, clinvar_varid_index = {}, {}

    # ------ Stage 4: ClinGen ------
    erepo_path = REF / "erepo.tabbed.txt"
    if erepo_path.exists():
        clingen_index = load_clingen(erepo_path)
    else:
        log.warning("erepo.tabbed.txt not found — skipping ClinGen annotations")
        clingen_index = {}

    # ------ Stage 5: Saudi AF ------
    af_dir = REF / "allele-frequencies"
    if af_dir.exists():
        saudi_index = load_saudi_af(af_dir, coord_set)
    else:
        log.warning("data/allele-frequencies/ not found — skipping Saudi AF")
        saudi_index = {}

    # ------ Stage 6: Gene-level ------
    gene_disease = load_gene_disease(REF / "genes_to_disease.txt")
    gnomad_constraint = load_gnomad_constraint(REF / "gnomad.v4.1.constraint_metrics.tsv")
    mouse_mp = load_mouse_phenotypes(
        REF / "HMD_HumanPhenotype.rpt",
        REF / "MGI_GenePheno.rpt",
        ONT / "mp.obo",
    )
    go_anns = load_go_annotations(REF / "goa_human.gaf.gz", ONT / "go-basic.obo")
    gtex_expr = load_gtex_expression(REF)

    # ------ Stage 7: Annotate and write ------
    log.info("Annotating %d rows ...", len(df))
    new_cols = {c: [] for c in NEW_VARIANT_COLS + NEW_GENE_COLS}

    for _, row in df.iterrows():
        ann = annotate_row(row, variant_map, vep_results, clinvar_index,
                           clinvar_varid_index, clingen_index, saudi_index,
                           gene_disease, gnomad_constraint, mouse_mp,
                           go_anns, gtex_expr)
        for col in new_cols:
            new_cols[col].append(ann.get(col, ""))

    for col, vals in new_cols.items():
        df[col] = vals

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    log.info("Wrote %d rows × %d columns → %s", len(df), len(df.columns), output_path)

    # Quick summary
    total = len(df)
    for col in ["vep_consequence", "vep_hgvsg", "clinvar_sig", "af_saudi", "expressed_in", "go_biological_process"]:
        filled = (df[col] != "").sum()
        log.info("  %-30s %d / %d (%.0f%%)", col, filled, total, 100 * filled / total if total else 0)


if __name__ == "__main__":
    main()
