"""
sparql_queries.py — SPARQL query templates for the PAVS backend.
"""

GRAPH_CASES    = "http://pavs.kaust.edu.sa/graph/cases"
GRAPH_GENES    = "http://pavs.kaust.edu.sa/graph/genes"
GRAPH_HPOA     = "http://pavs.kaust.edu.sa/graph/hpoa"
GRAPH_HPO_IC   = "http://pavs.kaust.edu.sa/graph/hpo-ic"
GRAPH_LIT      = "http://pavs.kaust.edu.sa/graph/literature"
GRAPH_DDD      = "http://pavs.kaust.edu.sa/graph/ddd"
GRAPH_HP_OWL   = "http://purl.obolibrary.org/obo/hp.owl"

PREFIXES = """
PREFIX pavs:    <http://pavs.kaust.edu.sa/ontology/>
PREFIX pav:     <http://pavs.kaust.edu.sa/data/>
PREFIX hp:      <http://purl.obolibrary.org/obo/HP_>
PREFIX mondo:   <http://purl.obolibrary.org/obo/MONDO_>
PREFIX omim:    <https://omim.org/entry/>
PREFIX hgnc:    <http://identifiers.org/hgnc.symbol/>
PREFIX dbsnp:   <http://identifiers.org/dbsnp/>
PREFIX clinvar: <http://identifiers.org/clinvar/>
PREFIX dc:      <http://purl.org/dc/terms/>
PREFIX rdfs:    <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd:     <http://www.w3.org/2001/XMLSchema#>
PREFIX owl:     <http://www.w3.org/2002/07/owl#>
"""


def hpo_autocomplete(q: str, limit: int = 20) -> str:
    """Autocomplete HPO terms by label substring."""
    q_lower = q.lower().replace('"', '\\"')
    return f"""
{PREFIXES}
SELECT DISTINCT ?id ?label WHERE {{
  GRAPH <{GRAPH_HPO_IC}> {{
    ?term rdfs:subClassOf* ?ancestor .
  }}
  GRAPH <{GRAPH_HP_OWL}> {{
    ?term rdfs:label ?label .
    FILTER(CONTAINS(LCASE(STR(?label)), "{q_lower}"))
  }}
  BIND(REPLACE(STR(?term), ".*HP_", "HP:") AS ?id)
  FILTER(STRSTARTS(?id, "HP:"))
}} LIMIT {limit}
"""


def hpo_autocomplete_simple(q: str, limit: int = 20) -> str:
    """Simpler HPO autocomplete using only the IC graph labels."""
    q_lower = q.lower().replace('"', '\\"')
    return f"""
{PREFIXES}
SELECT DISTINCT ?id ?label WHERE {{
  GRAPH <{GRAPH_HPO_IC}> {{
    ?term rdfs:subClassOf ?parent .
    OPTIONAL {{ ?term rdfs:label ?label }}
  }}
  FILTER(CONTAINS(LCASE(COALESCE(STR(?label), "")), "{q_lower}"))
  BIND(REPLACE(STR(?term), ".*HP_", "HP:") AS ?id)
  FILTER(STRSTARTS(?id, "HP:"))
}} LIMIT {limit}
"""


def load_all_ic() -> str:
    """Fetch all IC values from Virtuoso."""
    return f"""
{PREFIXES}
SELECT ?term ?ic WHERE {{
  GRAPH <{GRAPH_HPO_IC}> {{
    ?term pavs:informationContent ?ic .
  }}
}}
"""


def load_all_ancestors() -> str:
    """Fetch the HPO parent closure (subClassOf* is computed in Python)."""
    return f"""
{PREFIXES}
SELECT ?child ?parent WHERE {{
  GRAPH <{GRAPH_HPO_IC}> {{
    ?child rdfs:subClassOf ?parent .
  }}
}}
"""


def load_case_hpo_sets(include_non_saudi: bool = True) -> str:
    """Load all cases with their HPO sets for similarity scoring."""
    graphs = f"<{GRAPH_CASES}>"
    if include_non_saudi:
        graphs += f" <{GRAPH_LIT}> <{GRAPH_DDD}>"

    filter_clause = "" if include_non_saudi else 'FILTER(?isSaudi = true^^xsd:boolean)'

    return f"""
{PREFIXES}
SELECT ?case ?id ?gene ?disease ?source ?isSaudi
       (GROUP_CONCAT(DISTINCT ?hpo ; separator="|") AS ?hpos)
WHERE {{
  GRAPH ?g {{
    ?case a pavs:Case ;
          dc:identifier ?id ;
          pavs:source ?source ;
          pavs:isSaudi ?isSaudi .
    {filter_clause}
    OPTIONAL {{ ?case pavs:hasPhenotype ?hpoTerm .
               BIND(REPLACE(STR(?hpoTerm), ".*HP_", "HP:") AS ?hpo) }}
    OPTIONAL {{ ?case pavs:hasVariant ?v .
               ?v pavs:affectsGene ?geneUri .
               BIND(STRAFTER(STR(?geneUri), "hgnc.symbol/") AS ?gene) }}
    OPTIONAL {{ ?case pavs:diseaseLabel ?disease }}
  }}
  VALUES ?g {{ {graphs} }}
}}
GROUP BY ?case ?id ?gene ?disease ?source ?isSaudi
"""


def search_by_gene(gene: str, limit: int = 100) -> str:
    g_lower = gene.lower().replace('"', '\\"')
    return f"""
{PREFIXES}
SELECT DISTINCT ?case ?id ?disease ?source ?isSaudi ?hgvsC ?hgvsG ?acmg ?saudiAF WHERE {{
  GRAPH ?g {{
    ?case a pavs:Case ;
          dc:identifier ?id ;
          pavs:source ?source ;
          pavs:isSaudi ?isSaudi ;
          pavs:hasVariant ?v .
    ?v pavs:affectsGene ?geneUri .
    FILTER(CONTAINS(LCASE(STRAFTER(STR(?geneUri), "hgnc.symbol/")), "{g_lower}"))
    OPTIONAL {{ ?case pavs:diseaseLabel ?disease }}
    OPTIONAL {{ ?v pavs:hgvsC ?hgvsC }}
    OPTIONAL {{ ?v pavs:hgvsG ?hgvsG }}
    OPTIONAL {{ ?v pavs:acmgClass ?acmg }}
    OPTIONAL {{ ?v pavs:saudiAF ?saudiAF }}
  }}
}} LIMIT {limit}
"""


def search_by_rsid(rsid: str, limit: int = 100) -> str:
    rs_clean = rsid.strip().replace('"', '\\"')
    # Normalize: rs12345 → 12345 for URI
    rs_num = rs_clean.replace("rs", "").strip()
    return f"""
{PREFIXES}
SELECT DISTINCT ?case ?id ?disease ?source ?isSaudi ?hgvsC ?hgvsG ?acmg ?saudiAF ?gene WHERE {{
  GRAPH ?g {{
    ?case a pavs:Case ;
          dc:identifier ?id ;
          pavs:source ?source ;
          pavs:isSaudi ?isSaudi ;
          pavs:hasVariant ?v .
    ?v pavs:rsId dbsnp:rs{rs_num} .
    OPTIONAL {{ ?case pavs:diseaseLabel ?disease }}
    OPTIONAL {{ ?v pavs:hgvsC ?hgvsC }}
    OPTIONAL {{ ?v pavs:hgvsG ?hgvsG }}
    OPTIONAL {{ ?v pavs:acmgClass ?acmg }}
    OPTIONAL {{ ?v pavs:saudiAF ?saudiAF }}
    OPTIONAL {{ ?v pavs:affectsGene ?gUri .
               BIND(STRAFTER(STR(?gUri), "hgnc.symbol/") AS ?gene) }}
  }}
}} LIMIT {limit}
"""


def search_by_hgvs(hgvs: str, limit: int = 100) -> str:
    h = hgvs.strip().replace('"', '\\"')
    return f"""
{PREFIXES}
SELECT DISTINCT ?case ?id ?disease ?source ?isSaudi ?hgvsC ?hgvsG ?acmg ?saudiAF ?gene WHERE {{
  GRAPH ?g {{
    ?case a pavs:Case ;
          dc:identifier ?id ;
          pavs:source ?source ;
          pavs:isSaudi ?isSaudi ;
          pavs:hasVariant ?v .
    {{ ?v pavs:hgvsC ?hgvsC . FILTER(CONTAINS(?hgvsC, "{h}")) }}
    UNION
    {{ ?v pavs:hgvsG ?hgvsG . FILTER(CONTAINS(?hgvsG, "{h}")) }}
    OPTIONAL {{ ?case pavs:diseaseLabel ?disease }}
    OPTIONAL {{ ?v pavs:hgvsG ?hgvsG }}
    OPTIONAL {{ ?v pavs:acmgClass ?acmg }}
    OPTIONAL {{ ?v pavs:saudiAF ?saudiAF }}
    OPTIONAL {{ ?v pavs:affectsGene ?gUri .
               BIND(STRAFTER(STR(?gUri), "hgnc.symbol/") AS ?gene) }}
  }}
}} LIMIT {limit}
"""


def search_by_acmg(acmg_class: str, limit: int = 200) -> str:
    a = acmg_class.strip().replace('"', '\\"')
    return f"""
{PREFIXES}
SELECT DISTINCT ?case ?id ?disease ?source ?isSaudi ?hgvsC ?gene WHERE {{
  GRAPH ?g {{
    ?case a pavs:Case ;
          dc:identifier ?id ;
          pavs:source ?source ;
          pavs:isSaudi ?isSaudi ;
          pavs:hasVariant ?v .
    ?v pavs:acmgClass ?acmg .
    FILTER(CONTAINS(LCASE(?acmg), LCASE("{a}")))
    OPTIONAL {{ ?case pavs:diseaseLabel ?disease }}
    OPTIONAL {{ ?v pavs:hgvsC ?hgvsC }}
    OPTIONAL {{ ?v pavs:affectsGene ?gUri .
               BIND(STRAFTER(STR(?gUri), "hgnc.symbol/") AS ?gene) }}
  }}
}} LIMIT {limit}
"""


def search_by_disease(q: str, limit: int = 50) -> str:
    q_lower = q.lower().replace('"', '\\"')
    return f"""
{PREFIXES}
SELECT DISTINCT ?case ?id ?disease ?source ?isSaudi ?gene WHERE {{
  GRAPH ?g {{
    ?case a pavs:Case ;
          dc:identifier ?id ;
          pavs:source ?source ;
          pavs:isSaudi ?isSaudi ;
          pavs:diseaseLabel ?disease .
    FILTER(CONTAINS(LCASE(?disease), "{q_lower}"))
    OPTIONAL {{ ?case pavs:hasVariant ?v .
               ?v pavs:affectsGene ?gUri .
               BIND(STRAFTER(STR(?gUri), "hgnc.symbol/") AS ?gene) }}
  }}
}} LIMIT {limit}
"""


def get_case_detail(case_id: str) -> str:
    safe_id = case_id.replace('"', '\\"')
    return f"""
{PREFIXES}
CONSTRUCT {{
  ?case ?p ?o .
  ?v ?vp ?vo .
}}
WHERE {{
  GRAPH ?g {{
    ?case dc:identifier "{safe_id}" ;
          ?p ?o .
    OPTIONAL {{
      ?case pavs:hasVariant ?v .
      ?v ?vp ?vo .
    }}
  }}
}}
"""


def get_case_detail_select(case_id: str) -> str:
    safe_id = case_id.replace('"', '\\"')
    return f"""
{PREFIXES}
SELECT ?prop ?val WHERE {{
  GRAPH ?g {{
    ?case dc:identifier "{safe_id}" ;
          ?prop ?val .
  }}
}}
"""


def get_case_variants(case_id: str) -> str:
    safe_id = case_id.replace('"', '\\"')
    return f"""
{PREFIXES}
SELECT ?v ?gene ?hgvsC ?hgvsG ?hgvsP ?rsId ?zygosity ?acmg ?consequence
       ?sift ?polyphen ?gnomadAF ?saudiAF ?clinvarId ?clinvarSig
       ?vcfChrom ?vcfPos WHERE {{
  GRAPH ?g {{
    ?case dc:identifier "{safe_id}" ;
          pavs:hasVariant ?v .
    OPTIONAL {{ ?v pavs:affectsGene ?gUri .
               BIND(STRAFTER(STR(?gUri), "hgnc.symbol/") AS ?gene) }}
    OPTIONAL {{ ?v pavs:hgvsC ?hgvsC }}
    OPTIONAL {{ ?v pavs:hgvsG ?hgvsG }}
    OPTIONAL {{ ?v pavs:hgvsP ?hgvsP }}
    OPTIONAL {{ ?v pavs:rsId ?rsIdUri .
               BIND(STRAFTER(STR(?rsIdUri), "dbsnp/") AS ?rsId) }}
    OPTIONAL {{ ?v pavs:zygosity ?zygosity }}
    OPTIONAL {{ ?v pavs:acmgClass ?acmg }}
    OPTIONAL {{ ?v pavs:vepConsequence ?consequence }}
    OPTIONAL {{ ?v pavs:sift ?sift }}
    OPTIONAL {{ ?v pavs:polyphen ?polyphen }}
    OPTIONAL {{ ?v pavs:gnomadAF ?gnomadAF }}
    OPTIONAL {{ ?v pavs:saudiAF ?saudiAF }}
    OPTIONAL {{ ?v pavs:clinvarId ?clinvarId }}
    OPTIONAL {{ ?v pavs:clinvarSig ?clinvarSig }}
    OPTIONAL {{ ?v pavs:vcfChrom ?vcfChrom }}
    OPTIONAL {{ ?v pavs:vcfPos ?vcfPos }}
  }}
}}
"""


def get_hpo_labels(hpo_uris: list) -> str:
    """Fetch rdfs:label for a list of HP term URIs from the hpo-ic graph."""
    if not hpo_uris:
        return ""
    values = " ".join(f"<{uri}>" for uri in hpo_uris if uri)
    if not values:
        return ""
    return f"""
{PREFIXES}
SELECT ?term ?label WHERE {{
  GRAPH <{GRAPH_HPO_IC}> {{
    VALUES ?term {{ {values} }}
    OPTIONAL {{ ?term rdfs:label ?label }}
  }}
}}
"""


def get_gene_detail(gene_symbol: str) -> str:
    safe_sym = gene_symbol.replace('"', '\\"')
    return f"""
{PREFIXES}
SELECT ?prop ?val WHERE {{
  GRAPH <{GRAPH_GENES}> {{
    hgnc:{safe_sym} ?prop ?val .
  }}
}}
"""


def get_gene_diseases(gene_symbol: str) -> str:
    """Return all OMIM disease URIs linked to a gene from the genes graph."""
    safe = gene_symbol.replace('"', '\\"').replace("'", "\\'")
    return f"""
{PREFIXES}
SELECT ?disease WHERE {{
  GRAPH <{GRAPH_GENES}> {{
    hgnc:{safe} pavs:relatedDisease ?disease .
  }}
}}
"""


def hpoa_for_disease(disease_uri: str) -> str:
    return f"""
{PREFIXES}
SELECT ?hpo WHERE {{
  GRAPH <{GRAPH_HPOA}> {{
    ?assoc pavs:disease <{disease_uri}> ;
           pavs:hpoTerm ?hpoTerm .
    BIND(REPLACE(STR(?hpoTerm), ".*HP_", "HP:") AS ?hpo)
  }}
}}
"""
