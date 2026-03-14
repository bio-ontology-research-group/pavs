"""
Tests for the SPARQL Explorer example queries.

Runs each example query via the backend /api/sparql endpoint and asserts:
  - The query returns at least one row (non-empty results)
  - Known-true values appear in the results (spot-checks)

Run with:
  uv run pytest tests/test_sparql_queries.py -v

Requires the backend to be running at PAVS_API (default: http://localhost:8000).
"""

import os
import pytest
import requests

API = os.environ.get("PAVS_API", "http://localhost:8000")
ENDPOINT = f"{API}/api/sparql"

PREFIXES = """PREFIX pavs:  <http://pavs.kaust.edu.sa/ontology/>
PREFIX pav:   <http://pavs.kaust.edu.sa/data/>
PREFIX hgnc:  <http://identifiers.org/hgnc.symbol/>
PREFIX hp:    <http://purl.obolibrary.org/obo/HP_>
PREFIX dc:    <http://purl.org/dc/terms/>
PREFIX rdfs:  <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd:   <http://www.w3.org/2001/XMLSchema#>

"""


def run_query(query: str) -> dict:
    """Execute a SPARQL query via POST and return the parsed response dict."""
    resp = requests.post(ENDPOINT, json={"query": query}, timeout=30)
    resp.raise_for_status()
    return resp.json()


def run_query_get(query: str, fmt: str = "json") -> requests.Response:
    """Execute a SPARQL query via GET."""
    resp = requests.get(ENDPOINT, params={"query": query, "format": fmt}, timeout=30)
    resp.raise_for_status()
    return resp


# ──────────────────────────────────────────────────────────────────────────────
# Connectivity
# ──────────────────────────────────────────────────────────────────────────────

def test_backend_health():
    """Backend must be reachable."""
    resp = requests.get(f"{API}/api/health", timeout=10)
    assert resp.status_code == 200


# ──────────────────────────────────────────────────────────────────────────────
# Query 1: Most common HPO phenotypes (Saudi cases)
# ──────────────────────────────────────────────────────────────────────────────

QUERY_COMMON_PHENOTYPES = PREFIXES + """SELECT ?label (COUNT(DISTINCT ?case) AS ?n) WHERE {
  GRAPH <http://pavs.kaust.edu.sa/graph/cases> {
    ?case a pavs:Case ;
          pavs:hasPhenotype ?hpo .
  }
  GRAPH <http://pavs.kaust.edu.sa/graph/hpo-ic> {
    ?hpo rdfs:label ?label .
  }
}
GROUP BY ?label
ORDER BY DESC(?n)
LIMIT 20"""


def test_common_phenotypes_returns_results():
    data = run_query(QUERY_COMMON_PHENOTYPES)
    assert data["count"] > 0, "Expected at least one HPO phenotype row"


def test_common_phenotypes_vars():
    data = run_query(QUERY_COMMON_PHENOTYPES)
    assert "label" in data["vars"]
    assert "n" in data["vars"]


def test_common_phenotypes_counts_are_positive():
    data = run_query(QUERY_COMMON_PHENOTYPES)
    for row in data["rows"]:
        assert int(row["n"]) > 0, f"Expected positive count, got {row['n']}"


# ──────────────────────────────────────────────────────────────────────────────
# Query 2: Saudi cases for gene SUMF1
# ──────────────────────────────────────────────────────────────────────────────

QUERY_SUMF1 = PREFIXES + """SELECT ?id ?hgvsC ?acmg ?disease WHERE {
  GRAPH <http://pavs.kaust.edu.sa/graph/cases> {
    ?case a pavs:Case ;
          dc:identifier ?id ;
          pavs:hasVariant ?v .
    ?v pavs:affectsGene hgnc:SUMF1 .
    OPTIONAL { ?v pavs:hgvsC ?hgvsC }
    OPTIONAL { ?v pavs:acmgClass ?acmg }
    OPTIONAL { ?case pavs:diseaseLabel ?disease }
  }
}
ORDER BY ?id"""


def test_sumf1_returns_results():
    data = run_query(QUERY_SUMF1)
    assert data["count"] > 0, "Expected SUMF1 cases"


def test_sumf1_contains_known_case():
    """PAVS:A0000001 is a known SUMF1 case."""
    data = run_query(QUERY_SUMF1)
    ids = [row["id"] for row in data["rows"]]
    assert "PAVS:A0000001" in ids, f"Expected PAVS:A0000001 in SUMF1 results, got: {ids[:5]}"


def test_sumf1_known_case_has_disease():
    """PAVS:A0000001 should have a disease label containing 'sulfatid' or 'sulfat'."""
    data = run_query(QUERY_SUMF1)
    for row in data["rows"]:
        if row["id"] == "PAVS:A0000001":
            disease = (row.get("disease") or "").lower()
            assert "sulfat" in disease, f"Expected sulfatidosis-related disease, got: {disease}"
            return
    pytest.fail("PAVS:A0000001 not found in SUMF1 results")


def test_sumf1_known_case_has_hgvsc():
    """PAVS:A0000001 should have hgvsC containing 'c.785'."""
    data = run_query(QUERY_SUMF1)
    for row in data["rows"]:
        if row["id"] == "PAVS:A0000001":
            assert "c.785" in (row.get("hgvsC") or ""), \
                f"Expected c.785 in hgvsC, got: {row.get('hgvsC')}"
            return
    pytest.fail("PAVS:A0000001 not found in SUMF1 results")


# ──────────────────────────────────────────────────────────────────────────────
# Query 3: Pathogenic variants with disease labels
# ──────────────────────────────────────────────────────────────────────────────

QUERY_PATHOGENIC = PREFIXES + """SELECT ?id ?gene ?hgvsC ?acmg ?disease WHERE {
  GRAPH <http://pavs.kaust.edu.sa/graph/cases> {
    ?case a pavs:Case ;
          dc:identifier ?id ;
          pavs:diseaseLabel ?disease ;
          pavs:hasVariant ?v .
    ?v pavs:acmgClass ?acmg ;
       pavs:affectsGene ?gUri .
    FILTER(CONTAINS(LCASE(?acmg), "pathogenic"))
    OPTIONAL { ?v pavs:hgvsC ?hgvsC }
    BIND(STRAFTER(STR(?gUri), "hgnc.symbol/") AS ?gene)
  }
}
ORDER BY ?id
LIMIT 50"""


def test_pathogenic_returns_results():
    data = run_query(QUERY_PATHOGENIC)
    assert data["count"] > 0, "Expected pathogenic variant rows"


def test_pathogenic_acmg_values_contain_pathogenic():
    data = run_query(QUERY_PATHOGENIC)
    for row in data["rows"]:
        acmg = (row.get("acmg") or "").lower()
        assert "pathogenic" in acmg, f"Non-pathogenic ACMG value slipped through: {acmg}"


def test_pathogenic_has_gene_and_disease():
    data = run_query(QUERY_PATHOGENIC)
    for row in data["rows"]:
        assert row.get("gene"), f"Row missing gene: {row}"
        assert row.get("disease"), f"Row missing disease: {row}"


# ──────────────────────────────────────────────────────────────────────────────
# Query 4: Most constrained genes (LOEUF < 0.15)
# ──────────────────────────────────────────────────────────────────────────────

QUERY_CONSTRAINED = PREFIXES + """SELECT ?gene ?loeuf WHERE {
  GRAPH <http://pavs.kaust.edu.sa/graph/genes> {
    ?geneUri pavs:loeuf ?loeuf .
    BIND(STRAFTER(STR(?geneUri), "hgnc.symbol/") AS ?gene)
  }
  FILTER(xsd:decimal(?loeuf) < 0.15)
}
ORDER BY xsd:decimal(?loeuf)
LIMIT 30"""


def test_constrained_genes_returns_results():
    data = run_query(QUERY_CONSTRAINED)
    assert data["count"] > 0, "Expected constrained genes"


def test_constrained_genes_loeuf_below_threshold():
    data = run_query(QUERY_CONSTRAINED)
    for row in data["rows"]:
        loeuf = float(row["loeuf"])
        assert loeuf < 0.15, f"LOEUF {loeuf} exceeds threshold for gene {row.get('gene')}"


def test_constrained_genes_known_entry():
    """MED13L is a well-known constrained gene with very low LOEUF."""
    data = run_query(QUERY_CONSTRAINED)
    genes = [row["gene"] for row in data["rows"]]
    assert "MED13L" in genes, f"Expected MED13L in constrained genes, got top: {genes[:5]}"


# ──────────────────────────────────────────────────────────────────────────────
# Query 5: HPO terms for OMIM:272200 (Multiple sulfatase deficiency)
# ──────────────────────────────────────────────────────────────────────────────

QUERY_DISEASE_HPO = PREFIXES + """SELECT ?hpo ?label WHERE {
  GRAPH <http://pavs.kaust.edu.sa/graph/hpoa> {
    ?assoc pavs:disease <https://omim.org/entry/272200> ;
           pavs:hpoTerm ?hpoTerm .
    BIND(REPLACE(STR(?hpoTerm), ".*HP_", "HP:") AS ?hpo)
  }
  OPTIONAL {
    GRAPH <http://pavs.kaust.edu.sa/graph/hpo-ic> {
      ?hpoTerm rdfs:label ?label .
    }
  }
}"""


def test_disease_hpo_returns_results():
    data = run_query(QUERY_DISEASE_HPO)
    assert data["count"] > 0, "Expected HPO annotations for OMIM:272200"


def test_disease_hpo_contains_intellectual_disability():
    """OMIM:272200 (Multiple sulfatase deficiency) is annotated with intellectual disability."""
    data = run_query(QUERY_DISEASE_HPO)
    labels = [str(row.get("label") or "").lower() for row in data["rows"]]
    found = any("intellectual" in l or "developmental" in l or "disability" in l for l in labels)
    assert found, f"Expected intellectual disability phenotype for OMIM:272200, labels: {labels[:5]}"


def test_disease_hpo_ids_look_like_hp_terms():
    data = run_query(QUERY_DISEASE_HPO)
    for row in data["rows"]:
        hp = row.get("hpo", "")
        assert hp.startswith("HP:"), f"Expected HP:NNNNNNN format, got: {hp}"


# ──────────────────────────────────────────────────────────────────────────────
# Query 6: Cases with consanguineous parents
# ──────────────────────────────────────────────────────────────────────────────

QUERY_CONSANGUINEOUS = PREFIXES + """SELECT ?id ?gene ?disease WHERE {
  GRAPH <http://pavs.kaust.edu.sa/graph/cases> {
    ?case a pavs:Case ;
          dc:identifier ?id ;
          pavs:consanguinity "Offspring of consanguineous parents" .
    OPTIONAL { ?case pavs:diseaseLabel ?disease }
    OPTIONAL {
      ?case pavs:hasVariant ?v .
      ?v pavs:affectsGene ?gUri .
      BIND(STRAFTER(STR(?gUri), "hgnc.symbol/") AS ?gene)
    }
  }
}
LIMIT 50"""


def test_consanguineous_returns_results():
    data = run_query(QUERY_CONSANGUINEOUS)
    assert data["count"] > 0, "Expected consanguineous case rows"


def test_consanguineous_known_case():
    """PAVS:A0000001 is a known consanguineous SUMF1 case."""
    data = run_query(QUERY_CONSANGUINEOUS)
    ids = [row["id"] for row in data["rows"]]
    assert "PAVS:A0000001" in ids, f"Expected PAVS:A0000001 in consanguineous cases, got: {ids[:5]}"


def test_consanguineous_ids_start_with_pavs():
    data = run_query(QUERY_CONSANGUINEOUS)
    for row in data["rows"]:
        assert row["id"].startswith("PAVS:"), f"Unexpected case ID: {row['id']}"


# ──────────────────────────────────────────────────────────────────────────────
# GET endpoint + format tests
# ──────────────────────────────────────────────────────────────────────────────

SIMPLE_QUERY = PREFIXES + """SELECT ?id WHERE {
  GRAPH <http://pavs.kaust.edu.sa/graph/cases> {
    ?case a pavs:Case ;
          dc:identifier ?id .
  }
}
LIMIT 5"""


def test_get_endpoint_json():
    resp = run_query_get(SIMPLE_QUERY, fmt="json")
    data = resp.json()
    assert data["count"] == 5
    assert "id" in data["vars"]


def test_get_endpoint_tsv():
    resp = run_query_get(SIMPLE_QUERY, fmt="tsv")
    assert "text/tab-separated-values" in resp.headers.get("content-type", "")
    lines = resp.text.strip().split("\n")
    assert lines[0] == "id", f"Expected 'id' header, got: {lines[0]}"
    assert len(lines) == 6  # header + 5 rows
    for line in lines[1:]:
        assert line.startswith("PAVS:"), f"Expected PAVS: ID, got: {line}"


def test_get_endpoint_csv():
    resp = run_query_get(SIMPLE_QUERY, fmt="csv")
    assert "text/csv" in resp.headers.get("content-type", "")
    # csv module writes \r\n; strip both to be safe
    lines = [l.strip() for l in resp.text.strip().splitlines()]
    assert lines[0] == "id"
    assert len(lines) == 6


def test_get_endpoint_rejects_non_select():
    resp = requests.get(ENDPOINT, params={"query": "INSERT DATA { <x> <y> <z> }", "format": "json"}, timeout=10)
    assert resp.status_code == 400


def test_post_endpoint_rejects_non_select():
    resp = requests.post(ENDPOINT, json={"query": "DROP GRAPH <http://example.org/>"}, timeout=10)
    assert resp.status_code == 400
