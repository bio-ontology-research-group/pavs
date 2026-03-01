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

import gzip
import json
import logging
import os
import re
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

from .similarity import bma_similarity, bma_similarity_fast, expand_hpos
from . import sparql_queries as Q

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("pavs")

SPARQL_ENDPOINT = os.environ.get("SPARQL_ENDPOINT", "http://localhost:8890/sparql")
DATA_DIR = Path(os.environ.get("DATA_DIR", "data"))
HPO_TRANSLATION_PATH = Path(os.environ.get(
    "HPO_TRANSLATION_PATH",
    str(DATA_DIR.parent / "translation" / "hpo_arabic_translations.json"),
))
HPO_OBO_PATH = Path(os.environ.get(
    "HPO_OBO_PATH",
    str(DATA_DIR.parent / "ontology" / "hp.obo"),
))

app = FastAPI(title="PAVS SPARQL Search API", version="2.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ---------------------------------------------------------------------------
# Constants & Mappings
# ---------------------------------------------------------------------------

ZYGOSITY_MAP = {
    "http://purl.obolibrary.org/obo/GENO_0000135": {
        "en": "Heterozygous",
        "ar": "متباين الزيجوت",
        "desc_en": "Individual has two different alleles at a particular locus.",
        "desc_ar": "الفرد لديه أليلين مختلفين في موضع معين."
    },
    "http://purl.obolibrary.org/obo/GENO_0000136": {
        "en": "Homozygous",
        "ar": "متماثل الزيجوت",
        "desc_en": "Individual has two identical alleles at a particular locus.",
        "desc_ar": "الفرد لديه أليلين متطابقين في موضع معين."
    },
    "http://purl.obolibrary.org/obo/GENO_0000134": {
        "en": "Hemizygous",
        "ar": "أحادي الزيجوت",
        "desc_en": "Individual has only one allele at a particular locus (e.g., X-linked in males).",
        "desc_ar": "الفرد لديه أليل واحد فقط في موضع معين (على سبيل المثال، المرتبط بـ X في الذكور)."
    }
}

CONSEQUENCE_MAP = {
    "missense_variant": {"en": "Missense Variant", "ar": "تغير في تسلسل البروتين"},
    "stop_gained": {"en": "Stop Gained", "ar": "زيادة توقف"},
    "frameshift_variant": {"en": "Frameshift Variant", "ar": "تغيير إطار القراءة"},
    "splice_acceptor_variant": {"en": "Splice Acceptor Variant", "ar": "متغير مستقبل الوصل"},
    "splice_donor_variant": {"en": "Splice Donor Variant", "ar": "متغير مانح الوصل"},
    "synonymous_variant": {"en": "Synonymous Variant", "ar": "تغير صامت"},
    "inframe_insertion": {"en": "Inframe Insertion", "ar": "إدخال في الإطار"},
    "inframe_deletion": {"en": "Inframe Deletion", "ar": "حذف في الإطار"},
    "stop_lost": {"en": "Stop Lost", "ar": "فقدان توقف"},
    "start_lost": {"en": "Start Lost", "ar": "فقدان بدء"},
}

SIFT_MAP = {
    "deleterious": {"en": "Deleterious", "ar": "ضار"},
    "tolerated": {"en": "Tolerated", "ar": "محتمل"},
}

POLYPHEN_MAP = {
    "probably_damaging": {"en": "Probably Damaging", "ar": "من المرجح أن يكون ضاراً"},
    "possibly_damaging": {"en": "Possibly Damaging", "ar": "من المحتمل أن يكون ضاراً"},
    "benign": {"en": "Benign", "ar": "حميد"},
    "unknown": {"en": "Unknown", "ar": "غير معروف"},
}

CONSANGUINITY_MAP = {
    "yes": {"en": "Yes", "ar": "نعم"},
    "no": {"en": "No", "ar": "لا"},
    "Offspring of consanguineous parents": {"en": "Offspring of consanguineous parents", "ar": "من أبناء الأقارب"},
    "Consanguineous": {"en": "Consanguineous", "ar": "زواج أقارب"},
}

FIELD_DESCRIPTIONS = {
    "zygosity": {
        "en": "The genomic state of a variant (e.g., whether it is present on one or both copies of a chromosome).",
        "ar": "الحالة الجينية للمتغير (على سبيل المثال، ما إذا كان موجودًا في نسخة واحدة أو كلتا نسختي الكروموسوم)."
    },
    "consequence": {
        "en": "The predicted effect of the variant on the protein sequence (e.g., missense, stop-gain).",
        "ar": "التأثير المتوقع للمتغير على تسلسل البروتين (مثل التغير في تسلسل البروتين، زيادة توقف)."
    },
    "sift": {
        "en": "SIFT (Sorting Intolerant From Tolerant) predicts whether an amino acid substitution affects protein function.",
        "ar": "يتنبأ SIFT (فرز غير المتسامح من المتسامح) ما إذا كان استبدال الحمض الأميني يؤثر على وظيفة البروتين."
    },
    "polyphen": {
        "en": "PolyPhen-2 (Polymorphism Phenotyping v2) predicts the possible impact of an amino acid substitution on the structure and function of a protein.",
        "ar": "يتنبأ PolyPhen-2 (تنميط تعدد الأشكال الإصدار الثاني) بالتأثير المحتمل لاستبدال الحمض الأميني على بنية ووظيفة البروتين."
    }
}

# ---------------------------------------------------------------------------
# In-memory caches (loaded at startup)
# ---------------------------------------------------------------------------
ic_cache: Dict[str, float] = {}
ancestor_cache: Dict[str, Set[str]] = {}
descendant_cache: Dict[str, Set[str]] = {}    # HP:ID → set of all descendant HP IDs (incl. self)
case_hpo_cache: List[Dict[str, Any]] = []
disease_label_cache: Dict[str, str] = {}   # "OMIM:272200" → "Multiple sulfatase deficiency"
gene_disease_cache: Dict[str, List[str]] = {}  # "SUMF1" → ["Multiple sulfatase deficiency"]
hpo_ar_cache: Dict[str, Dict[str, str]] = {}  # "HP:0001263" → {arabic_label, arabic_layperson, arabic_definition}
hpo_obo_cache: Dict[str, Dict] = {}           # "HP:0001263" → {definition, layperson_synonyms, synonyms}
children_cache: Dict[str, Set[str]] = {}       # HP:ID → set of direct child HP IDs
term_case_count: Dict[str, int] = {}           # HP:ID → propagated case count (all cohorts)
term_saudi_count: Dict[str, int] = {}          # HP:ID → propagated case count (Saudi only)
disease_hpo_cache: Dict[str, List[str]] = {}   # "OMIM:272200" → [HP:NNNNNNN, ...]


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


def _load_disease_hpos(hpoa_path: Path) -> Dict[str, List[str]]:
    """Parse phenotype.hpoa and return {OMIM:NNNNNN → [HP:NNNNNNN, ...]} dict (excludes NOT qualifiers)."""
    result: Dict[str, List[str]] = defaultdict(list)
    try:
        with open(hpoa_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("database_id"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                disease_id = parts[0].strip()   # e.g. OMIM:272200
                qualifier = parts[2].strip()    # e.g. "NOT" or ""
                hpo_id = parts[3].strip()       # e.g. HP:0001263
                if disease_id and hpo_id.startswith("HP:") and qualifier != "NOT":
                    if hpo_id not in result[disease_id]:
                        result[disease_id].append(hpo_id)
        log.info(f"  Loaded HPO annotations for {len(result)} diseases from phenotype.hpoa")
    except Exception as e:
        log.warning(f"Could not load disease HPO annotations from {hpoa_path}: {e}")
    return dict(result)


def _vcf_info_val(info_str: str, key: str) -> Optional[str]:
    """Extract a value from a semicolon-delimited VCF INFO string."""
    prefix = f"{key}="
    for field in info_str.split(";"):
        if field.startswith(prefix):
            return field[len(prefix):]
    return None


def _load_clinvar_cases(
    vcf_path: Path,
    disease_hpo_cache: Dict[str, List[str]],
    disease_label_cache: Dict[str, str],
    ancestor_cache: Dict[str, Set[str]],
) -> List[Dict[str, Any]]:
    """Parse clinvar.vcf.gz and build one synthetic case per unique (gene, OMIM_disease) pair
    for pathogenic / likely-pathogenic variants. HPO phenotypes come from disease_hpo_cache."""
    PATHOGENIC_TERMS = {"Pathogenic", "Likely_pathogenic"}
    seen: Dict[tuple, bool] = {}
    cases: List[Dict[str, Any]] = []
    skipped_no_hpo = 0

    try:
        open_fn = gzip.open if str(vcf_path).endswith(".gz") else open
        with open_fn(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t", 8)
                if len(parts) < 8:
                    continue
                info_str = parts[7]

                clnsig = _vcf_info_val(info_str, "CLNSIG")
                if not clnsig:
                    continue
                sig_parts = set(re.split(r"[/|,]", clnsig))
                if not (sig_parts & PATHOGENIC_TERMS):
                    continue

                geneinfo = _vcf_info_val(info_str, "GENEINFO")
                if not geneinfo or geneinfo == ".":
                    continue
                genes = [g.split(":")[0] for g in geneinfo.split("|") if ":" in g]
                if not genes:
                    continue

                clndisdb = _vcf_info_val(info_str, "CLNDISDB")
                if not clndisdb or clndisdb == ".":
                    continue
                omim_ids = []
                for entry in re.split(r"[|]", clndisdb):
                    for db_ref in entry.split(","):
                        db_ref = db_ref.strip()
                        if db_ref.startswith("OMIM:") and not db_ref.startswith("OMIM_"):
                            num = db_ref[5:]
                            if num.isdigit():
                                omim_ids.append(f"OMIM:{num}")

                for gene in genes:
                    for omim_id in omim_ids:
                        key = (gene, omim_id)
                        if key in seen:
                            continue
                        seen[key] = True
                        hpo_ids = disease_hpo_cache.get(omim_id, [])
                        if not hpo_ids:
                            skipped_no_hpo += 1
                            continue
                        hpo_expanded = expand_hpos(hpo_ids, ancestor_cache)
                        disease_label = disease_label_cache.get(omim_id, omim_id)
                        cases.append({
                            "id": f"ClinVar:{gene}:{omim_id}",
                            "case_uri": "",
                            "gene": gene,
                            "disease": omim_id,          # OMIM ID for HPO lookup
                            "disease_label": disease_label,  # human-readable name for display
                            "source": "ClinVar",
                            "is_saudi": False,
                            "hpo_ids": hpo_ids,
                            "hpo_expanded": hpo_expanded,
                        })
        log.info(
            f"  Loaded {len(cases)} ClinVar (gene, disease) cases "
            f"({skipped_no_hpo} pairs skipped — no HPOA phenotypes)"
        )
    except Exception as e:
        log.warning(f"Could not load ClinVar cases from {vcf_path}: {e}")
    return cases


def _load_gene_disease(tsv_path: Path, disease_labels: Dict[str, str]) -> Dict[str, List[str]]:
    """Parse genes_to_disease.txt → {gene_symbol → [disease_label, ...]}."""
    result: Dict[str, List[str]] = defaultdict(list)
    try:
        with open(tsv_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("ncbi_gene_id") or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                gene = parts[1].strip()
                disease_id = parts[3].strip()   # e.g. OMIM:272200
                if gene and disease_id:
                    label = disease_labels.get(disease_id, disease_id)
                    if label not in result[gene]:
                        result[gene].append(label)
        log.info(f"  Loaded gene→disease for {len(result)} genes")
    except Exception as e:
        log.warning(f"Could not load gene→disease from {tsv_path}: {e}")
    return dict(result)


def _load_hpo_arabic(json_path: Path) -> Dict[str, Dict[str, str]]:
    """Load hpo_arabic_translations.json → {HP:ID → {arabic_label, arabic_layperson, arabic_definition}}."""
    result: Dict[str, Dict[str, str]] = {}
    try:
        with open(json_path, encoding="utf-8") as fh:
            terms = json.load(fh)
        for term in terms:
            hp_id = term.get("id", "")
            if hp_id:
                result[hp_id] = {
                    "arabic_label": term.get("arabic_technical_name", ""),
                    "arabic_layperson": term.get("arabic_layperson_synonym", ""),
                    "arabic_definition": term.get("arabic_definition", ""),
                }
        log.info(f"  Loaded {len(result)} Arabic HPO translations")
    except Exception as e:
        log.warning(f"Could not load Arabic HPO translations from {json_path}: {e}")
    return result


def _parse_hpo_obo(obo_path: Path) -> Dict[str, Dict]:
    """Parse hp.obo and return {HP:ID → {definition, layperson_synonyms, synonyms}}."""
    data: Dict[str, Dict] = {}
    current: Dict[str, Any] = {}
    in_term = False
    try:
        with open(obo_path, encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip()
                if line == "[Term]":
                    if current.get("id", "").startswith("HP:"):
                        data[current["id"]] = current
                    current = {}
                    in_term = True
                elif line == "" and in_term:
                    pass  # blank line within stanza — keep going
                elif not in_term:
                    continue
                elif line.startswith("id: "):
                    current["id"] = line[4:].strip()
                elif line.startswith("name: "):
                    current["name"] = line[6:].strip()
                elif line.startswith("def: "):
                    raw = line[5:].strip()
                    if raw.startswith('"'):
                        end_q = raw.rfind('"', 1)
                        if end_q > 0:
                            current["definition"] = raw[1:end_q]
                elif line.startswith("synonym: "):
                    raw = line[9:].strip()
                    if raw.startswith('"'):
                        end_q = raw.rfind('"', 1)
                        if end_q > 0:
                            syn_val = raw[1:end_q]
                            scope_rest = raw[end_q + 1:].strip()
                            if "layperson" in scope_rest:
                                current.setdefault("layperson_synonyms", []).append(syn_val)
                            else:
                                # Only keep EXACT synonyms that differ from name
                                if "EXACT" in scope_rest:
                                    current.setdefault("synonyms", []).append(syn_val)
        # Flush last term
        if current.get("id", "").startswith("HP:"):
            data[current["id"]] = current
        log.info(f"  Parsed {len(data)} terms from hp.obo")
    except Exception as e:
        log.warning(f"Could not parse hp.obo from {obo_path}: {e}")
    return data


@app.on_event("startup")
async def startup():
    global ic_cache, ancestor_cache, descendant_cache, case_hpo_cache, disease_label_cache, gene_disease_cache, children_cache, term_case_count, term_saudi_count, disease_hpo_cache

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

        # Build children_cache (parent → direct children) by inverting direct_parents
        ch: Dict[str, Set[str]] = defaultdict(set)
        for child, parents in direct_parents.items():
            for p in parents:
                ch[p].add(child)
        children_cache.update(ch)
        log.info(f"  Built children cache for {len(children_cache)} parent terms")
    except Exception as e:
        log.warning(f"Ancestor load failed: {e}")

    log.info("Building descendant cache from ancestor sets …")
    try:
        dd: Dict[str, Set[str]] = defaultdict(set)
        for child, ancs in ancestor_cache.items():
            for anc in ancs:
                dd[anc].add(child)
        descendant_cache.update(dd)
        log.info(f"  Built descendant cache for {len(descendant_cache)} terms")
    except Exception as e:
        log.warning(f"Descendant cache build failed: {e}")

    log.info("Loading case HPO sets …")
    try:
        rows = sparql_select(Q.load_case_hpo_sets(include_non_saudi=True))
        case_hpo_cache.clear()
        for row in rows:
            hpos_raw = row.get("hpos", "")
            hpo_ids = [_hp_id(h) for h in hpos_raw.split("|") if h] if hpos_raw else []
            hpo_expanded = expand_hpos(hpo_ids, ancestor_cache)
            case_hpo_cache.append({
                "id": row.get("id", ""),
                "case_uri": row.get("case", ""),
                "gene": row.get("gene", ""),
                "disease": row.get("disease", ""),
                "source": row.get("source", "unknown"),
                "is_saudi": row.get("isSaudi", "0") in ("true", "True", "1"),
                "hpo_ids": hpo_ids,
                "hpo_expanded": hpo_expanded,
            })
        log.info(f"  Loaded {len(case_hpo_cache)} cases into similarity cache (with pre-expanded HPO sets)")

        # Build term_case_count (all cohorts) and term_saudi_count (Saudi only)
        from collections import Counter
        cnt: Counter = Counter()
        saudi_cnt: Counter = Counter()
        for case in case_hpo_cache:
            for term in case["hpo_expanded"]:
                cnt[term] += 1
                if case.get("is_saudi"):
                    saudi_cnt[term] += 1
        term_case_count.update(cnt)
        term_saudi_count.update(saudi_cnt)
        log.info(f"  Built case counts for {len(term_case_count)} terms ({len(term_saudi_count)} Saudi)")
    except Exception as e:
        log.warning(f"Case HPO load failed: {e}")

    log.info("Loading disease labels and HPO annotations from phenotype.hpoa …")
    hpoa_path = DATA_DIR / "reference" / "phenotype.hpoa"
    disease_label_cache.update(_load_disease_labels(hpoa_path))
    log.info(f"  Loaded {len(disease_label_cache)} disease labels")
    disease_hpo_cache.update(_load_disease_hpos(hpoa_path))
    log.info(f"  Loaded HPO annotations for {len(disease_hpo_cache)} diseases")

    log.info("Loading ClinVar pathogenic cases …")
    clinvar_vcf = DATA_DIR / "reference" / "clinvar.vcf.gz"
    if clinvar_vcf.exists():
        clinvar_cases = _load_clinvar_cases(clinvar_vcf, disease_hpo_cache, disease_label_cache, ancestor_cache)
        case_hpo_cache.extend(clinvar_cases)
        log.info(f"  case_hpo_cache now has {len(case_hpo_cache)} entries (incl. ClinVar)")
    else:
        log.warning(f"  ClinVar VCF not found at {clinvar_vcf}, skipping")

    log.info("Loading gene→disease map …")
    gene_disease_tsv = DATA_DIR / "reference" / "genes_to_disease.txt"
    gene_disease_cache.update(_load_gene_disease(gene_disease_tsv, disease_label_cache))

    log.info("Loading Arabic HPO translations …")
    hpo_ar_cache.update(_load_hpo_arabic(HPO_TRANSLATION_PATH))

    log.info("Parsing hp.obo for definitions and synonyms …")
    hpo_obo_cache.update(_parse_hpo_obo(HPO_OBO_PATH))

    log.info("PAVS backend ready.")


# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------

class PhenotypeSearchRequest(BaseModel):
    hpo_ids: List[str]
    method: str = "lin"            # "lin" (default) or "resnik"
    limit: int = 20
    include_disease_phenotypes: bool = False
    include_saudi: bool = True
    include_ddd: bool = False
    include_literature: bool = False
    include_clinvar: bool = False


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


@app.get("/api/hpo/{hp_id:path}")
def get_hpo_term(hp_id: str):
    """Return enriched HPO term info: English definition, layperson synonyms, Arabic translations."""
    # Normalise HP:0001263 / HP_0001263 / 0001263
    hp_id = hp_id.strip()
    if hp_id.startswith("HP_"):
        hp_id = "HP:" + hp_id[3:]
    elif not hp_id.startswith("HP:"):
        hp_id = "HP:" + hp_id

    obo = hpo_obo_cache.get(hp_id, {})
    ar = hpo_ar_cache.get(hp_id, {})

    # Deduplicate synonyms: exclude the term name itself
    term_name = obo.get("name", "")
    synonyms = [s for s in obo.get("synonyms", []) if s != term_name]
    layperson_synonyms = obo.get("layperson_synonyms", [])

    return {
        "id": hp_id,
        "definition": obo.get("definition", ""),
        "layperson_synonyms": layperson_synonyms,
        "synonyms": synonyms,
        "arabic_label": ar.get("arabic_label", ""),
        "arabic_layperson": ar.get("arabic_layperson", ""),
        "arabic_definition": ar.get("arabic_definition", ""),
    }


@app.get("/api/hpo-children/{hp_id}")
def get_hpo_children(hp_id: str):
    """Return direct children of an HPO term with propagated case counts.
    Used by the HPO tree browser for lazy-loading."""
    hp_id = hp_id.strip()
    if hp_id.startswith("HP_"):
        hp_id = "HP:" + hp_id[3:]
    elif not hp_id.startswith("HP:"):
        hp_id = "HP:" + hp_id

    child_ids = children_cache.get(hp_id, set())
    result = []
    for cid in child_ids:
        obo = hpo_obo_cache.get(cid, {})
        ar = hpo_ar_cache.get(cid, {})
        result.append({
            "id": cid,
            "label": obo.get("name", cid),
            "definition": obo.get("definition", ""),
            "arabic_label": ar.get("arabic_label", ""),
            "case_count": term_case_count.get(cid, 0),
            "saudi_case_count": term_saudi_count.get(cid, 0),
            "child_count": len(children_cache.get(cid, set())),
        })
    result.sort(key=lambda x: x["saudi_case_count"], reverse=True)
    return result


@app.post("/api/search/phenotype")
def search_by_phenotype(req: PhenotypeSearchRequest):
    """Rank all cases by HPO similarity (Lin+BMA by default)."""
    if not req.hpo_ids:
        raise HTTPException(400, "hpo_ids must be non-empty")

    method = req.method if req.method in ("lin", "resnik") else "lin"

    # Pre-expand query terms once
    q_exp = expand_hpos(req.hpo_ids, ancestor_cache)

    results = []
    for case in case_hpo_cache:
        src = case.get("source", "")
        is_saudi = case.get("is_saudi", False)
        is_ddd = "ddd" in src.lower()
        is_clinvar = src == "ClinVar"
        is_lit = not is_saudi and not is_ddd and not is_clinvar

        if is_saudi and not req.include_saudi:
            continue
        if is_ddd and not req.include_ddd:
            continue
        if is_lit and not req.include_literature:
            continue
        if is_clinvar and not req.include_clinvar:
            continue

        # Optionally expand with disease HPO terms (modifies target only, uses in-memory cache)
        if req.include_disease_phenotypes and case.get("disease"):
            disease_hpos = disease_hpo_cache.get(case["disease"], [])
            if disease_hpos:
                extra = list(set(case["hpo_ids"]) | set(disease_hpos))
                t_exp = expand_hpos(extra, ancestor_cache)
            else:
                t_exp = case.get("hpo_expanded") or expand_hpos(case["hpo_ids"], ancestor_cache)
        else:
            t_exp = case.get("hpo_expanded") or expand_hpos(case["hpo_ids"], ancestor_cache)

        if not t_exp:
            continue

        score = bma_similarity_fast(q_exp, t_exp, ic_cache, ancestor_cache, method)
        if score > 0:
            gene = case.get("gene", "")
            disease_raw = case.get("disease", "")
            # Prefer an explicit label; for ClinVar cases use disease_label; otherwise raw value
            disease_display = case.get("disease_label") or disease_raw

            # Suggested disease for Saudi cases that have no explicit diagnosis
            suggested_disease = ""
            if is_saudi and not disease_raw and gene:
                assoc = gene_disease_cache.get(gene, [])
                if 1 <= len(assoc) <= 3:
                    suggested_disease = "; ".join(assoc)

            results.append({
                "id": case["id"],
                "gene": gene,
                "disease": disease_display,
                "suggested_disease": suggested_disease,
                "source": case["source"],
                "is_saudi": case["is_saudi"],
                "hpo_ids": case["hpo_ids"],
                "score": round(score, 6),
            })

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
        disease_uris: List[str] = []

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
            elif prop_short == "hasDisease":
                disease_uris.append(val)
            else:
                case_data["properties"][prop_short] = val

        # Normalise isSaudi boolean ("1" from Virtuoso xsd:boolean)
        is_saudi = case_data["properties"].get("isSaudi")
        if is_saudi in ("1", "true", "True"):
            case_data["properties"]["isSaudi"] = "true"

        # Map progressStatus to human-readable string (keep raw key for frontend to translate)
        # but normalize for consistency if needed.
        raw_status = case_data["properties"].get("progressStatus", "")
        # We keep it as is, or we could map to a standard set of keys.
        # But for now, let's just make sure Solved/Unsolved are available.

        # Localize consanguinity
        raw_consang = case_data["properties"].get("consanguinity", "")
        if raw_consang in CONSANGUINITY_MAP:
            case_data["properties"]["consanguinity_label"] = CONSANGUINITY_MAP[raw_consang]["en"]
            case_data["properties"]["consanguinity_ar"] = CONSANGUINITY_MAP[raw_consang]["ar"]

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

        def uri_to_disease(uri: str) -> Dict[str, str]:
            # omim: https://omim.org/entry/123456 -> OMIM:123456
            if "omim.org/entry/" in uri:
                did = "OMIM:" + uri.split("/")[-1]
            # mondo: http://purl.obolibrary.org/obo/MONDO_0001234 -> MONDO:0001234
            elif "MONDO_" in uri:
                did = "MONDO:" + uri.split("MONDO_")[-1]
            # orphanet: http://www.orpha.net/ORDO/Orphanet_123 -> Orphanet:123
            elif "ORDO/Orphanet_" in uri:
                did = "Orphanet:" + uri.split("Orphanet_")[-1]
            else:
                did = uri.split("/")[-1]
            
            return {"id": did, "label": disease_label_cache.get(did, ""), "url": uri}

        case_data["diseases"] = [uri_to_disease(u) for u in disease_uris]

        # Variants — merge multiple SPARQL rows for the same variant URI
        var_rows = sparql_select(Q.get_case_variants(case_id))
        variants_dict: Dict[str, Dict[str, Any]] = {}
        for row in var_rows:
            v_uri = row.get("v", "")
            if not v_uri:
                continue
            if v_uri not in variants_dict:
                variants_dict[v_uri] = {"v": v_uri}
            
            # Merge all available fields into the existing dict
            for k, val in row.items():
                if val and k not in variants_dict[v_uri]:
                    variants_dict[v_uri][k] = val
            
            # Derive TogoVar link if possible
            if "togovar_url" not in variants_dict[v_uri]:
                tv = _togovar_url(row)
                if tv:
                    variants_dict[v_uri]["togovar_url"] = tv
        
        case_data["variants"] = list(variants_dict.values())

        # Enrich variant labels and translations
        for var in case_data["variants"]:
            # Zygosity
            z_uri = var.get("zygosity")
            if z_uri and z_uri in ZYGOSITY_MAP:
                var["zygosity_label"] = ZYGOSITY_MAP[z_uri]["en"]
                var["zygosity_ar"] = ZYGOSITY_MAP[z_uri]["ar"]
            
            # Consequence (raw: "missense_variant")
            cons_raw = var.get("consequence")
            if cons_raw:
                c_map = CONSEQUENCE_MAP.get(cons_raw)
                if c_map:
                    var["consequence_label"] = c_map["en"]
                    var["consequence_ar"] = c_map["ar"]
                else:
                    # Fallback: simple title case and same for AR if no mapping
                    var["consequence_label"] = cons_raw.replace("_", " ").title()
                    var["consequence_ar"] = cons_raw

            # SIFT (raw: "deleterious(0)")
            sift_raw = var.get("sift")
            if sift_raw:
                import re
                m = re.match(r"^([a-z_]+)\(([\d\.]+)\)$", sift_raw.lower())
                if m:
                    label, score = m.groups()
                    s_map = SIFT_MAP.get(label)
                    if s_map:
                        var["sift_label"] = s_map["en"]
                        var["sift_ar"] = s_map["ar"]
                    var["sift_score"] = score

            # PolyPhen (raw: "probably_damaging(0.99)")
            pp_raw = var.get("polyphen")
            if pp_raw:
                import re
                m = re.match(r"^([a-z_]+)\(([\d\.]+)\)$", pp_raw.lower())
                if m:
                    label, score = m.groups()
                    p_map = POLYPHEN_MAP.get(label)
                    if p_map:
                        var["polyphen_label"] = p_map["en"]
                        var["polyphen_ar"] = p_map["ar"]
                    var["polyphen_score"] = score

        # Field descriptions for the "i" buttons
        case_data["field_metadata"] = FIELD_DESCRIPTIONS

        # Suggested OMIM diseases (when case has no diagnosed disease)
        # Only shown for Saudi cases with 1-3 diseases for the causative gene.
        case_data["suggested_diseases"] = []
        has_disease = bool(case_data["properties"].get("diseaseLabel"))
        is_saudi = case_data["properties"].get("isSaudi") == "true"
        if is_saudi and not has_disease and case_data["variants"]:
            seen_genes: set = set()
            for var in case_data["variants"]:
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


@app.get("/api/gene/{symbol}/diseases")
def gene_diseases(symbol: str):
    """Return diseases associated with a gene and their HPO phenotypes from HPOA."""
    try:
        rows = sparql_select(Q.get_gene_diseases_hpo(symbol))
        # Group by disease URI
        disease_map: Dict[str, Dict] = {}
        for row in rows:
            uri = row.get("diseaseUri", "")
            if not uri:
                continue
            if uri not in disease_map:
                omim_num = uri.split("/")[-1]
                omim_id = f"OMIM:{omim_num}"
                disease_map[uri] = {
                    "disease_id": omim_id,
                    "label": disease_label_cache.get(omim_id, omim_id),
                    "omim_url": uri,
                    "clinvar_url": f"https://www.ncbi.nlm.nih.gov/clinvar/?term={symbol}[gene]",
                    "hpo_terms": [],
                }
            hp_id = row.get("hpo", "")
            label = row.get("hpoLabel", "")
            if hp_id:
                disease_map[uri]["hpo_terms"].append({"id": hp_id, "label": label or hp_id})
        return list(disease_map.values())
    except Exception as e:
        log.error(f"Gene diseases error for {symbol}: {e}")
        raise HTTPException(500, str(e))


@app.get("/api/genes")
def list_genes():
    """List all genes from Saudi cases with case counts."""
    try:
        rows = sparql_select(Q.list_all_genes())
        result = []
        for r in rows:
            sym = r.get("gene", "")
            cnt = r.get("saudiCount", "0")
            if sym:
                try:
                    result.append({"symbol": sym, "saudi_count": int(cnt)})
                except (ValueError, TypeError):
                    result.append({"symbol": sym, "saudi_count": 0})
        return result
    except Exception as e:
        log.error(f"Gene list error: {e}")
        raise HTTPException(500, str(e))


@app.get("/api/gene/{symbol}/cases")
def gene_cases(
    symbol: str,
    include_ddd: bool = False,
    include_literature: bool = False,
):
    """Return cases with variants in a specific gene, with optional DDD/Literature sources."""
    try:
        graphs_list = [f"<{Q.GRAPH_CASES}>"]
        if include_ddd:
            graphs_list.append(f"<{Q.GRAPH_DDD}>")
        if include_literature:
            graphs_list.append(f"<{Q.GRAPH_LIT}>")
        graphs_str = " ".join(graphs_list)
        rows = sparql_select(Q.get_gene_cases(symbol, graphs_str))
        # Fallback disease labels from genes_to_disease.txt
        fallback_diseases = gene_disease_cache.get(symbol, [])
        fallback_label = "; ".join(fallback_diseases) if fallback_diseases else ""
        results = []
        for row in rows:
            r = {k: v for k, v in row.items() if v is not None}
            r["gene"] = symbol
            if not r.get("disease") and fallback_label:
                r["disease"] = fallback_label
            tv = _togovar_url(row)
            if tv:
                r["togovar_url"] = tv
            results.append(r)
        return results
    except Exception as e:
        log.error(f"Gene cases error for {symbol}: {e}")
        raise HTTPException(500, str(e))


@app.get("/api/phenotype/{hp_id:path}/cases")
def phenotype_cases(
    hp_id: str,
    include_children: bool = True,
    include_ddd: bool = False,
    include_literature: bool = False,
    limit: int = 500,
):
    """Return cases that have (or whose descendants have) the given HPO term.
    Uses the in-memory case_hpo_cache + descendant_cache — no SPARQL at query time.
    """
    # Normalise the HP ID
    hp_id = hp_id.strip()
    if hp_id.startswith("HP_"):
        hp_id = "HP:" + hp_id[3:]
    elif not hp_id.startswith("HP:"):
        hp_id = "HP:" + hp_id

    # Build the set of HP IDs to match (term itself + optionally its descendants)
    match_ids: Set[str] = {hp_id}
    if include_children:
        match_ids |= descendant_cache.get(hp_id, set())

    results = []
    for case in case_hpo_cache:
        source = case.get("source", "")
        is_saudi = case.get("is_saudi", False)
        is_ddd = "ddd" in source.lower()
        is_lit = not is_saudi and not is_ddd

        if not is_saudi and not (include_ddd and is_ddd) and not (include_literature and is_lit):
            continue

        case_hpos = set(case["hpo_ids"])
        if case_hpos & match_ids:
            results.append({
                "id": case["id"],
                "gene": case["gene"],
                "disease": case["disease"],
                "source": source,
                "is_saudi": is_saudi,
                "hpo_ids": case["hpo_ids"],
            })
        if len(results) >= limit:
            break

    return results


def _sparql_check_and_run(query: str):
    """Validate and execute a SELECT query; return (vars_list, rows)."""
    if not query:
        raise HTTPException(400, "query required")
    upper = query.lstrip().upper()
    if "SELECT" not in upper:
        raise HTTPException(400, "Only SELECT queries are supported")
    try:
        rows = sparql_select(query)
        rows = rows[:2000]
        vars_list = list(rows[0].keys()) if rows else []
        return vars_list, rows
    except Exception as e:
        raise HTTPException(500, str(e))


@app.get("/api/sparql")
def sparql_get(
    query: str = Query(..., description="SPARQL SELECT query"),
    format: str = Query("json", description="Response format: json | tsv | csv"),
):
    """Execute a SPARQL SELECT query and return results in the requested format.

    Programmatic access example:
      curl 'http://localhost:8000/api/sparql?format=tsv&query=SELECT+...'
    """
    vars_list, rows = _sparql_check_and_run(query)

    if format == "tsv":
        lines = ["\t".join(vars_list)]
        for row in rows:
            lines.append("\t".join(str(row.get(v, "")) for v in vars_list))
        content = "\n".join(lines)
        return StreamingResponse(
            iter([content]),
            media_type="text/tab-separated-values",
            headers={"Content-Disposition": "attachment; filename=sparql_results.tsv"},
        )
    elif format == "csv":
        import csv
        import io
        out = io.StringIO()
        writer = csv.writer(out)
        writer.writerow(vars_list)
        for row in rows:
            writer.writerow([str(row.get(v, "")) for v in vars_list])
        return StreamingResponse(
            iter([out.getvalue()]),
            media_type="text/csv",
            headers={"Content-Disposition": "attachment; filename=sparql_results.csv"},
        )
    else:
        return {"vars": vars_list, "rows": rows, "count": len(rows)}


@app.post("/api/sparql")
def sparql_proxy(body: dict):
    """Execute a user-supplied SPARQL SELECT query against Virtuoso and return rows.
    Restricted to SELECT queries; max 2000 rows returned."""
    query = (body.get("query") or "").strip()
    vars_list, rows = _sparql_check_and_run(query)
    return {"vars": vars_list, "rows": rows, "count": len(rows)}


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
