# PAVS reviewer-revision: linear walkthrough and proof of work

*2026-06-10T03:44:58Z by Showboat 0.6.1*
<!-- showboat-id: eee4f304-8e04-4786-b653-42d851b5f98b -->

This document re-derives every number changed in the revised manuscript directly from the data, verifies the live SPARQL endpoint, and records one data-deployment issue found during the walkthrough (the deployed knowledge graph is missing variant annotations that the canonical generator produces). Run `uvx showboat verify REVISION_WALKTHROUGH.md` to re-execute everything.

## 0. Setup — ensure the Phenopacket Store (v0.1.26) is available.

```bash
test -d /tmp/ppktstore/0.1.26 || { mkdir -p /tmp/ppktstore && curl -sL -o /tmp/ppktstore/all.zip https://github.com/monarch-initiative/phenopacket-store/releases/download/0.1.26/all_phenopackets.zip && (cd /tmp/ppktstore && unzip -q -o all.zip); }; echo "PAVS phenopackets: $(.venv/bin/python -c "import json;print(len(json.load(open(\"data/PAVS_phenopackets.json\"))))")"; echo "Store phenopackets: $(find /tmp/ppktstore/0.1.26 -name "*.json" | wc -l)"
```

```output
PAVS phenopackets: 7510
Store phenopackets: 9588
```

## 1. Cohort split (Reviewer 2.1b) — Table 5 prioritization, reproduced from per-case ranking results; matches published merged-Saudi numbers exactly.

```bash
.venv/bin/python evaluation/similarity/scripts/cohort_split_analysis.py 2>&1 | grep -E 'ext-resnik|Wrote'
```

```output
  Literature         ext-resnik  8887 0.9828 0.6542 57.78 77.87
  Saudi-literature   ext-resnik  1018 0.9315 0.1287 8.35 21.91
  DDD                ext-resnik  1443 0.9951 0.7063 62.02 85.45
  Mixed              ext-resnik   122 0.8933 0.0344 1.64 5.74
  Saudi-clinical     ext-resnik  1948 0.8705 0.0377 1.39 7.96
Wrote evaluation/similarity/cohort_split_prioritization.csv
```

Source, not population, drives top-rank accuracy: Saudi-clinical (1.39%) approx Mixed non-Saudi clinical (1.64%), while Saudi case reports (8.35%) are six-fold higher. All AUCs high (0.87-0.99).

## 2. Density split (Table 6) — counts all features (matching similarity.py); reproduces published DDD/Literature means; identical term count (~3.8) across clinical cohorts, differing IC (case reports 2.95 vs clinical 2.46).

```bash
.venv/bin/python evaluation/similarity/scripts/cohort_split_density.py --obo .venv/lib/python3.10/site-packages/pyhpo/data/hp.obo --g2p .venv/lib/python3.10/site-packages/pyhpo/data/genes_to_phenotype.txt --store /tmp/ppktstore/0.1.26 2>&1 | tail -7
```

```output
Cohort                 n  MeanHPO     SD  Max  MeanIC  IC SD
Saudi-clinical      3253     3.75   2.62   26   2.457  1.118
Saudi-literature    1244     3.76   3.43   32   2.947  1.301
Mixed                505     3.81   2.27   15   2.615  1.225
DDD                 1456    19.33  15.82  117   3.234  0.821
Literature          9587    21.47  21.87  192   3.470  0.932
```

## 3. Phenotype-less cases (Reviewer 2.2) — 1,052 (14%) have no HPO term; most remain genotype-informative.

```bash
.venv/bin/python - <<'PY'
import json,re
from collections import Counter
d=json.load(open('data/PAVS_phenopackets.json'))
cm={'A':'Saudi-clinical','B':'Saudi-clinical','F':'Saudi-clinical','P':'Saudi-clinical','M':'Saudi-literature','Q':'Mixed','D':'DDD'}
noph=[];hasvar=0;hasdis=0
for pp in d:
    hp=[f for f in pp.get('phenotypicFeatures',[]) if f.get('type',{}).get('id','').startswith('HP:') and not f.get('excluded',False)]
    if not hp:
        p=re.match(r'PAVS:([A-Za-z]+)',pp['id']).group(1)
        noph.append(cm.get(p,p)); hasvar+=bool(pp.get('interpretations')); hasdis+=bool(pp.get('diseases'))
print('phenotype-less cases:',len(noph),'of 7510 (%.1f%%)'%(100*len(noph)/7510))
print('by cohort:',dict(Counter(noph)))
print('with variant interpretation:',hasvar,'(%.1f%%)'%(100*hasvar/len(noph)),'| with disease:',hasdis)
PY
```

```output
phenotype-less cases: 1053 of 7510 (14.0%)
by cohort: {'Saudi-clinical': 458, 'Saudi-literature': 178, 'Mixed': 17, 'DDD': 400}
with variant interpretation: 912 (86.6%) | with disease: 523
```

## 4. Negation & modifier prevalence (Reviewer 2.3) — rare in the corpus; a curation set for expert grading is generated (44 negations graded exhaustively + 100 sampled modifiers).

```bash
.venv/bin/python - <<'PY'
import json
d=json.load(open('data/PAVS_phenopackets.json'))
tot=neg=mod=0
for pp in d:
    for f in pp.get('phenotypicFeatures',[]):
        if not f.get('type',{}).get('id','').startswith('HP:'): continue
        tot+=1; neg+=bool(f.get('excluded',False)); mod+=bool(f.get('modifiers'))
print(f'total HP assignments: {tot}')
print(f'negated (excluded): {neg} ({100*neg/tot:.1f}%)')
print(f'with severity modifier: {mod} ({100*mod/tot:.1f}%)')
PY
```

```output
total HP assignments: 47294
negated (excluded): 44 (0.1%)
with severity modifier: 898 (1.9%)
```

```bash
.venv/bin/python evaluation/similarity/scripts/sample_negation_modifier_eval.py --raw combined_normalized-gpt4oss.tsv --out /tmp/nmeval 2>&1 | grep -E 'Negated|Modifier'; echo "curation rows: negations=$(grep -c . /tmp/nmeval/negation_eval_set.csv | tr -d '\r'), counting header"
```

```output
Negated assignments (grade ALL): 44 -> /tmp/nmeval/negation_eval_set.csv
Modifier assignments total: 958; sampled 100 -> /tmp/nmeval/modifier_eval_set.csv
curation rows: negations=45, counting header
```

## 5. Cross-source overlap with the Phenopacket Store (Reviewer 1.1) — publication-level upper bound on duplication.

```bash
.venv/bin/python evaluation/similarity/scripts/crosssource_pmid_overlap.py --pavs data/PAVS_phenopackets.json --store /tmp/ppktstore/0.1.26 2>&1 | tail -8
```

```output
Saudi case reports (prefix M): 1422
  with >=1 source PMID:        1361
  distinct source PMIDs:       103
Store phenopackets scanned:    9588  (distinct PMIDs: 1584)

Overlapping PMIDs:             5
Saudi case reports sharing a Store PMID: 129 (9.1% of M cases)
PMIDs: ['24369382', '25539947', '28600779', '29297947', '32431071']
```

## 6. Live SPARQL — the website's corrected example query (originally queried pavs:zygosity, which the deployed graph lacks; see section 8). Reproduces the manuscript's ELAC2/ATP7B statistic.

```bash
Q="PREFIX pavs: <https://pavs.phenomebrowser.net/ontology/>
SELECT ?gene (COUNT(DISTINCT ?case) AS ?nSaudiCases) WHERE {
  GRAPH <https://pavs.phenomebrowser.net/graph/cases> {
    ?case a pavs:Case ; pavs:isSaudi true ; pavs:hasVariant ?v .
    ?v pavs:affectsGene ?gUri . BIND(STRAFTER(STR(?gUri),\"hgnc.symbol/\") AS ?gene) }
} GROUP BY ?gene ORDER BY DESC(?nSaudiCases) LIMIT 5"
curl -sG --max-time 40 https://pavs.phenomebrowser.net/sparql --data-urlencode "query=$Q" --data-urlencode "format=csv" | tr -d "\r" | head -8
```

```output
"gene","nSaudiCases"
"ELAC2",50
"ATP7B",42
"TULP1",39
"ACADVL",35
"ISCA2",33
```

## 7. SourceBadge ↔ data consistency — the badge keys are the pavs:source strings in the graph; all are covered, so every case renders a correct population x provenance badge. (Note: PMC7082194/Ziats is tagged isSaudi=true in the data but treated as Mixed-population per the manuscript; the badge follows the manuscript.)

```bash
Q="PREFIX pavs: <https://pavs.phenomebrowser.net/ontology/>
SELECT ?source (COUNT(?c) AS ?n) WHERE {
  GRAPH <https://pavs.phenomebrowser.net/graph/cases> { ?c a pavs:Case ; pavs:source ?source }
} GROUP BY ?source ORDER BY DESC(?n)"
curl -sG --max-time 40 https://pavs.phenomebrowser.net/sparql --data-urlencode "query=$Q" --data-urlencode "format=csv" | tr -d "\r" | head -10
```

```output
"source","n"
"PMC6562004",2217
"ddd-diagnoses",1856
"marwa-variants",1422
"fawzan-variants",1019
"PMC7082194",522
"ahmed-variants",251
"ahmed-pmid28454995",223
```

## 8. ISSUE FOUND — stale deployed knowledge graph. The generator (pavs-knowledge-graph/intake/generate_rdf.py) emits pavs:zygosity, vepConsequence, sift, polyphen, gnomadAF, rsId, clinvarSig; the backend queries them (OPTIONAL); but the DEPLOYED graph carries only acmgClass. The manuscript's variant-annotation counts and the case-detail variant panel are not reflected live. Data is correct; deployment is behind. Remedy: regenerate RDF with generate_rdf.py and reload Virtuoso. (A redeploy will flip the zeros below to non-zero, which verify will flag.)

```bash
Q="PREFIX pavs: <https://pavs.phenomebrowser.net/ontology/>
SELECT (COUNT(?z) AS ?zygosity) (COUNT(?vep) AS ?vep) (COUNT(?g) AS ?gnomadAF) (COUNT(?rs) AS ?rsId) (COUNT(?a) AS ?acmgClass) WHERE {
  GRAPH <https://pavs.phenomebrowser.net/graph/cases> { ?v a pavs:GenomicVariant .
    OPTIONAL{?v pavs:zygosity ?z} OPTIONAL{?v pavs:vepConsequence ?vep} OPTIONAL{?v pavs:gnomadAF ?g}
    OPTIONAL{?v pavs:rsId ?rs} OPTIONAL{?v pavs:acmgClass ?a} } }"
curl -sG --max-time 40 https://pavs.phenomebrowser.net/sparql --data-urlencode "query=$Q" --data-urlencode "format=csv" | tr -d "\r" | head -4
```

```output
"zygosity","vep","gnomadAF","rsId","acmgClass"
0,0,0,0,6093
```
