"""
backend_sparql/main.py — FastAPI SPARQL proxy with HPO similarity engine.

Startup sequence:
1. Fetch IC values from Virtuoso → ic_cache
2. Fetch HPO parent closure → ancestor_cache
3. Fetch all case HPO sets → case_hpo_cache

All displayable data still fetched from SPARQL at query time.

Environment variables:
  SPARQL_ENDPOINT          (default: http://localhost:8890/sparql)
  SPARQL_UPDATE_ENDPOINT   (default: http://localhost:8890/sparql-auth)
"""

from __future__ import annotations

import json
import logging
import os
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, RedirectResponse, StreamingResponse
from pydantic import BaseModel

from .similarity import bma_similarity
from . import sparql_queries as Q

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("pavs")

SPARQL_ENDPOINT = os.environ.get("SPARQL_ENDPOINT", "http://localhost:8890/sparql")
DATA_DIR = Path(os.environ.get("DATA_DIR", "data"))

app = FastAPI(title="PAVS SPARQL Search API", version="2.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ---------------------------------------------------------------------------
# In-memory caches (loaded at startup)
# ---------------------------------------------------------------------------
ic_cache: Dict[str, float] = {}
ancestor_cache: Dict[str, Set[str]] = {}
case_hpo_cache: List[Dict[str, Any]] = []
disease_label_cache: Dict[str, str] = {}   # "OMIM:272200" → "Multiple sulfatase deficiency"


# ---------------------------------------------------------------------------
# SPARQL helpers
# ---------------------------------------------------------------------------

def sparql_select(query: str, retries: int = 3) -> List[Dict[str, Any]]:
    """Execute a SELECT query and return list of binding dicts."""
    encoded = urllib.parse.urlencode({"query": query, "format": "application/sparql-results+json"})
    url = f"{SPARQL_ENDPOINT}?{encoded}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                data = json.loads(resp.read())
                return [
                    {k: v.get("value") for k, v in row.items()}
                    for row in data.get("results", {}).get("bindings", [])
                ]
        except Exception as e:
            if attempt == retries - 1:
                log.error(f"SPARQL query failed: {e}\nQuery: {query[:200]}")
                raise
            time.sleep(1 + attempt)
    return []


def _hp_id(uri_or_id: str) -> str:
    """Normalize HP URI or curie to HP:NNNNNNN format."""
    if uri_or_id.startswith("HP:"):
        return uri_or_id
    # http://purl.obolibrary.org/obo/HP_0001263 → HP:0001263
    if "HP_" in uri_or_id:
        return "HP:" + uri_or_id.split("HP_")[-1]
    # hp:0001263 → HP:0001263
    if uri_or_id.startswith("hp:"):
        return "HP:" + uri_or_id[3:]
    return uri_or_id


# ---------------------------------------------------------------------------
# Startup: load caches
# ---------------------------------------------------------------------------

def _build_ancestor_sets(direct_parents: Dict[str, Set[str]]) -> Dict[str, Set[str]]:
    """Compute transitive ancestor closure via iterative BFS."""
    ancestors: Dict[str, Set[str]] = {}
    all_terms = set(direct_parents.keys())
    for ps in direct_parents.values():
        all_terms |= ps

    def get_anc(term: str) -> Set[str]:
        if term in ancestors:
            return ancestors[term]
        anc: Set[str] = {term}
        for p in direct_parents.get(term, set()):
            anc |= get_anc(p)
        ancestors[term] = anc
        return anc

    for t in all_terms:
        get_anc(t)
    return ancestors


def _load_disease_labels(hpoa_path: Path) -> Dict[str, str]:
    """Parse phenotype.hpoa and return {OMIM:NNNNNN → disease_name} dict."""
    labels: Dict[str, str] = {}
    try:
        with open(hpoa_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("database_id"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    disease_id = parts[0].strip()   # e.g. OMIM:272200
                    name = parts[1].strip()
                    if disease_id and name and disease_id not in labels:
                        labels[disease_id] = name
    except Exception as e:
        log.warning(f"Could not load disease labels from {hpoa_path}: {e}")
    return labels


@app.on_event("startup")
async def startup():
    global ic_cache, ancestor_cache, case_hpo_cache, disease_label_cache

    log.info("Loading HPO IC values from Virtuoso …")
    try:
        rows = sparql_select(Q.load_all_ic())
        for row in rows:
            term_uri = row.get("term", "")
            ic_val = row.get("ic", "")
            if term_uri and ic_val:
                hp_id = _hp_id(term_uri)
                try:
                    ic_cache[hp_id] = float(ic_val)
                except ValueError:
                    pass
        log.info(f"  Loaded {len(ic_cache)} IC values")
    except Exception as e:
        log.warning(f"IC load failed (Virtuoso may not be ready): {e}")

    log.info("Loading HPO ancestor sets …")
    try:
        rows = sparql_select(Q.load_all_ancestors())
        direct_parents: Dict[str, Set[str]] = defaultdict(set)
        for row in rows:
            child = _hp_id(row.get("child", ""))
            parent = _hp_id(row.get("parent", ""))
            if child and parent:
                direct_parents[child].add(parent)
        ancestor_cache = _build_ancestor_sets(direct_parents)
        log.info(f"  Loaded ancestor sets for {len(ancestor_cache)} terms")
    except Exception as e:
        log.warning(f"Ancestor load failed: {e}")

    log.info("Loading case HPO sets …")
    try:
        rows = sparql_select(Q.load_case_hpo_sets(include_non_saudi=True))
        case_hpo_cache.clear()
        for row in rows:
            hpos_raw = row.get("hpos", "")
            hpo_ids = [_hp_id(h) for h in hpos_raw.split("|") if h] if hpos_raw else []
            case_hpo_cache.append({
                "id": row.get("id", ""),
                "case_uri": row.get("case", ""),
                "gene": row.get("gene", ""),
                "disease": row.get("disease", ""),
                "source": row.get("source", "unknown"),
                "is_saudi": row.get("isSaudi", "0") in ("true", "True", "1"),
                "hpo_ids": hpo_ids,
            })
        log.info(f"  Loaded {len(case_hpo_cache)} cases into similarity cache")
    except Exception as e:
        log.warning(f"Case HPO load failed: {e}")

    log.info("Loading disease labels from phenotype.hpoa …")
    hpoa_path = DATA_DIR / "phenotype.hpoa"
    disease_label_cache.update(_load_disease_labels(hpoa_path))
    log.info(f"  Loaded {len(disease_label_cache)} disease labels")

    log.info("PAVS backend ready.")


# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------

class PhenotypeSearchRequest(BaseModel):
    hpo_ids: List[str]
    method: str = "lin"            # "lin" (default) or "resnik"
    limit: int = 20
    include_disease_phenotypes: bool = False
    include_non_saudi: bool = False


class VariantSearchRequest(BaseModel):
    gene: Optional[str] = None
    rsid: Optional[str] = None
    hgvs: Optional[str] = None
    acmg_class: Optional[str] = None
    limit: int = 100


# ---------------------------------------------------------------------------
# API endpoints
# ---------------------------------------------------------------------------

@app.get("/api/health")
def health():
    return {
        "status": "ok",
        "ic_terms": len(ic_cache),
        "ancestor_terms": len(ancestor_cache),
        "cases_cached": len(case_hpo_cache),
    }


@app.get("/api/search/hpo")
def hpo_autocomplete(q: str = Query(..., min_length=2)):
    """HPO term autocomplete via SPARQL."""
    try:
        rows = sparql_select(Q.hpo_autocomplete_simple(q, limit=20))
        results = []
        seen = set()
        for row in rows:
            hid = row.get("id", "")
            label = row.get("label", "")
            if hid and hid not in seen:
                seen.add(hid)
                results.append({"id": hid, "label": label or hid})
        return results
    except Exception as e:
        log.error(f"HPO autocomplete error: {e}")
        return []


@app.post("/api/search/phenotype")
def search_by_phenotype(req: PhenotypeSearchRequest):
    """Rank all cases by HPO similarity (Lin+BMA by default)."""
    if not req.hpo_ids:
        raise HTTPException(400, "hpo_ids must be non-empty")

    method = req.method if req.method in ("lin", "resnik") else "lin"

    results = []
    for case in case_hpo_cache:
        if not req.include_non_saudi and not case["is_saudi"]:
            continue

        target_hpos = case["hpo_ids"]
        if not target_hpos:
            continue

        # Optionally expand with disease HPO terms
        if req.include_disease_phenotypes and case.get("disease"):
            try:
                disease_uri = f"https://omim.org/entry/{case['disease'].replace('OMIM:', '')}"
                hpoa_rows = sparql_select(Q.hpoa_for_disease(disease_uri))
                disease_hpos = [_hp_id(r.get("hpo", "")) for r in hpoa_rows if r.get("hpo")]
                target_hpos = list(set(target_hpos) | set(disease_hpos))
            except Exception:
                pass

        score = bma_similarity(req.hpo_ids, target_hpos, ic_cache, ancestor_cache, method)
        if score > 0:
            results.append({**case, "score": round(score, 6)})

    results.sort(key=lambda x: x["score"], reverse=True)
    return results[: req.limit]


@app.get("/api/search/gene")
def search_by_gene(q: str = Query(..., min_length=1), limit: int = 100):
    try:
        rows = sparql_select(Q.search_by_gene(q, limit))
        return _format_variant_results(rows)
    except Exception as e:
        log.error(f"Gene search error: {e}")
        raise HTTPException(500, str(e))


@app.get("/api/search/variant")
def search_variant(
    gene: Optional[str] = None,
    rsid: Optional[str] = None,
    hgvs: Optional[str] = None,
    acmg: Optional[str] = None,
    limit: int = 100,
):
    try:
        if rsid:
            rows = sparql_select(Q.search_by_rsid(rsid, limit))
        elif hgvs:
            rows = sparql_select(Q.search_by_hgvs(hgvs, limit))
        elif gene:
            rows = sparql_select(Q.search_by_gene(gene, limit))
        elif acmg:
            rows = sparql_select(Q.search_by_acmg(acmg, limit))
        else:
            raise HTTPException(400, "Provide gene, rsid, hgvs, or acmg parameter")
        results = _format_variant_results(rows)
        return results
    except HTTPException:
        raise
    except Exception as e:
        log.error(f"Variant search error: {e}")
        raise HTTPException(500, str(e))


_STATUS_MAP = {
    "SOLVED": "Solved",
    "UNSOLVED": "Unsolved",
    "IN_PROGRESS": "Unknown",
    "UNKNOWN_PROGRESS": "Unknown",
}


def _togovar_url(row: Dict) -> Optional[str]:
    """Return a /api/togovar-search proxy path when we have chrom+pos."""
    chrom = row.get("vcfChrom") or row.get("vcfchrom")
    pos = row.get("vcfPos") or row.get("vcfpos")
    if chrom and pos:
        # Strip "chr" prefix — TogoVar expects bare chromosome numbers/letters
        chrom_clean = chrom.lstrip("chr") if chrom.lower().startswith("chr") else chrom
        return f"/api/togovar-search?chrom={chrom_clean}&pos={pos}"
    return None


def _format_variant_results(rows: List[Dict]) -> List[Dict]:
    results = []
    for row in rows:
        r = {k: v for k, v in row.items()}
        tv = _togovar_url(row)
        if tv:
            r["togovar_url"] = tv
        results.append(r)
    return results


@app.get("/api/search/disease")
def search_disease(q: str = Query(..., min_length=2), limit: int = 50):
    try:
        rows = sparql_select(Q.search_by_disease(q, limit))
        return rows
    except Exception as e:
        log.error(f"Disease search error: {e}")
        raise HTTPException(500, str(e))


@app.get("/api/case/{case_id:path}")
def get_case(case_id: str):
    """Return full case details: demographics, phenotypes, variants."""
    try:
        # Properties — accumulate multi-valued props properly
        props_rows = sparql_select(Q.get_case_detail_select(case_id))
        case_data: Dict[str, Any] = {"id": case_id, "properties": {}}
        phenotype_uris: List[str] = []
        excluded_phenotype_uris: List[str] = []

        for row in props_rows:
            prop = row.get("prop", "")
            val = row.get("val", "")
            if not prop or val is None:
                continue
            prop_short = prop.split("/")[-1].split("#")[-1]
            if prop_short == "hasPhenotype":
                phenotype_uris.append(val)
            elif prop_short == "hasExcludedPhenotype":
                excluded_phenotype_uris.append(val)
            else:
                case_data["properties"][prop_short] = val

        # Normalise isSaudi boolean ("1" from Virtuoso xsd:boolean)
        is_saudi = case_data["properties"].get("isSaudi")
        if is_saudi in ("1", "true", "True"):
            case_data["properties"]["isSaudi"] = "true"

        # Map progressStatus to human-readable string
        raw_status = case_data["properties"].get("progressStatus", "")
        if raw_status:
            case_data["properties"]["progressStatus"] = _STATUS_MAP.get(raw_status, raw_status)

        # Look up HPO labels for phenotypes
        all_uris = phenotype_uris + excluded_phenotype_uris
        label_map: Dict[str, str] = {}
        if all_uris:
            try:
                label_rows = sparql_select(Q.get_hpo_labels(all_uris))
                for lr in label_rows:
                    term_uri = lr.get("term", "")
                    lbl = lr.get("label", "")
                    if term_uri and lbl:
                        label_map[term_uri] = lbl
            except Exception as le:
                log.warning(f"HPO label lookup failed: {le}")

        def uri_to_hpo(uri: str) -> Dict[str, str]:
            hid = "HP:" + uri.split("HP_")[-1] if "HP_" in uri else uri
            return {"id": hid, "label": label_map.get(uri, "")}

        case_data["phenotypes"] = [uri_to_hpo(u) for u in phenotype_uris]
        case_data["excluded_phenotypes"] = [uri_to_hpo(u) for u in excluded_phenotype_uris]

        # Variants
        var_rows = sparql_select(Q.get_case_variants(case_id))
        variants = []
        for row in var_rows:
            v = {k: val for k, val in row.items() if val}
            tv = _togovar_url(row)
            if tv:
                v["togovar_url"] = tv
            variants.append(v)
        case_data["variants"] = variants

        # Suggested OMIM diseases (when case has no diagnosed disease)
        # Only shown for Saudi cases with 1-3 diseases for the causative gene.
        case_data["suggested_diseases"] = []
        has_disease = bool(case_data["properties"].get("diseaseLabel"))
        is_saudi = case_data["properties"].get("isSaudi") == "true"
        if is_saudi and not has_disease and variants:
            seen_genes: set = set()
            for var in variants:
                gene = var.get("gene", "")
                if not gene or gene in seen_genes:
                    continue
                seen_genes.add(gene)
                try:
                    disease_rows = sparql_select(Q.get_gene_diseases(gene))
                    disease_uris = [r.get("disease", "") for r in disease_rows if r.get("disease")]
                    if 1 <= len(disease_uris) <= 3:
                        for uri in disease_uris:
                            omim_num = uri.split("/")[-1]
                            omim_id = f"OMIM:{omim_num}"
                            label = disease_label_cache.get(omim_id, omim_id)
                            case_data["suggested_diseases"].append({
                                "id": omim_id,
                                "label": label,
                                "url": uri,
                                "gene": gene,
                            })
                except Exception as se:
                    log.warning(f"Gene disease lookup failed for {gene}: {se}")

        # Phenopacket download link
        case_data["phenopacket_download_url"] = f"/api/phenopacket/{case_id}/download"

        return case_data
    except Exception as e:
        log.error(f"Case detail error: {e}")
        raise HTTPException(500, str(e))


@app.get("/api/gene/{symbol}")
def get_gene(symbol: str):
    """Return gene details from the genes graph."""
    try:
        rows = sparql_select(Q.get_gene_detail(symbol))
        props: Dict[str, Any] = {"symbol": symbol}
        for row in rows:
            prop = row.get("prop", "")
            val = row.get("val", "")
            prop_short = prop.split("/")[-1].split("#")[-1]
            props[prop_short] = val
        props["hgnc_url"] = f"http://identifiers.org/hgnc.symbol/{symbol}"
        return props
    except Exception as e:
        log.error(f"Gene detail error: {e}")
        raise HTTPException(500, str(e))


@app.get("/api/phenopacket/{case_id:path}/download")
def download_phenopacket(case_id: str):
    """Download the phenopacket JSON for a specific case."""
    # Try generated_v2 directory first, then generated
    for subdir in ["generated_v2", "generated"]:
        path = Path("phenopackets") / subdir / f"{case_id}.json"
        if path.exists():
            return FileResponse(
                path,
                media_type="application/json",
                filename=f"{case_id}.json",
            )
    raise HTTPException(404, f"Phenopacket not found: {case_id}")


@app.get("/api/togovar-search")
def togovar_search(chrom: str = Query(...), pos: int = Query(...)):
    """Proxy variant lookup to TogoVar API and redirect to the variant page.

    TogoVar expects chromosome without 'chr' prefix (e.g. '3', '16', 'X').
    Returns 302 redirect to https://grch38.togovar.org/variant/{tgv_id} if found,
    or 404 if the variant is not in the TogoVar database.
    """
    chrom_clean = chrom.lstrip("chr") if chrom.lower().startswith("chr") else chrom
    payload = json.dumps({
        "query": {"location": {"chromosome": chrom_clean, "position": pos}},
        "limit": 1,
    }).encode()
    req_obj = urllib.request.Request(
        "https://grch38.togovar.org/api/search/variant",
        data=payload,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req_obj, timeout=10) as resp:
            body = json.loads(resp.read())
            results = body.get("data", [])
            if results:
                tgv_id = results[0].get("id", "")
                if tgv_id:
                    return RedirectResponse(
                        f"https://grch38.togovar.org/variant/{tgv_id}",
                        status_code=302,
                    )
    except Exception as e:
        log.warning(f"TogoVar API call failed for {chrom_clean}:{pos} — {e}")
    # Variant not in TogoVar (e.g. population-specific); send to homepage
    return RedirectResponse("https://grch38.togovar.org/", status_code=302)


@app.get("/api/phenopackets/download-all")
def download_all_phenopackets():
    """Stream the PAVS_phenopackets.zip bundle."""
    zip_path = DATA_DIR / "PAVS_phenopackets.zip"
    if not zip_path.exists():
        raise HTTPException(404, "Combined phenopackets zip not available")
    return FileResponse(
        zip_path,
        media_type="application/zip",
        filename="PAVS_phenopackets.zip",
    )
