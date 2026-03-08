# PAVS RDF Schema Documentation

The PAVS project uses a custom RDF schema to represent clinical genomics data, standardized to international ontologies.

## Namespace Prefixes

| Prefix | URI |
|--------|-----|
| `pavs:` | `http://pavs.kaust.edu.sa/ontology/` |
| `pav:` | `http://pavs.kaust.edu.sa/data/` |
| `hp:` | `http://purl.obolibrary.org/obo/HP_` |
| `mondo:` | `http://purl.obolibrary.org/obo/MONDO_` |
| `hancestro:` | `http://purl.obolibrary.org/obo/HANCESTRO_` |
| `gaz:` | `http://purl.obolibrary.org/obo/GAZ_` |
| `omim:` | `https://omim.org/entry/` |
| `hgnc:` | `http://identifiers.org/hgnc.symbol/` |

---

## Classes

### `pavs:Case`
Represents an individual clinical case (patient).
- `dc:identifier`: Unique identifier (e.g., `PAVS:A0000001`).
- `pavs:source`: Source dataset name.
- `pavs:isSaudi`: Boolean flag for Saudi cohort membership.
- `pavs:sex`: Patient sex (NCIT URI).
- `pavs:hasPhenotype`: Link to HPO term (present).
- `pavs:hasExcludedPhenotype`: Link to HPO term (negated).
- `pavs:hasDisease`: Link to OMIM or MONDO disease.
- `pavs:hasVariant`: Link to a `pavs:Variant` instance.
- `pavs:ancestry`: Link to HANCESTRO term (e.g., `hancestro:0852` for Middle Eastern).
- `pavs:geographicOrigin`: Link to GAZ term (e.g., `gaz:00005279` for Saudi Arabia).

### `pavs:Variant`
Represents a specific genomic variant identified in a case.
- `pavs:affectsGene`: Link to `pavs:Gene` (HGNC symbol URI).
- `pavs:hgvsC`, `pavs:hgvsP`, `pavs:hgvsG`: HGVS notations.
- `pavs:acmgClass`: ACMG pathogenicity classification string.
- `pavs:vepConsequence`: VEP most severe consequence.
- `pavs:vepImpact`: VEP impact (HIGH, MODERATE, etc.).
- `pavs:sift`, `pavs:polyphen`: Functional predictions.
- `pavs:gnomadAF`: Global gnomAD Allele Frequency.
- `pavs:gnomadAF_MID`: Middle Eastern gnomAD AF.
- `pavs:gnomadAF_SAS`: South Asian gnomAD AF.
- `pavs:saudiAF`, `pavs:saudiAC`: Local Saudi cohort Allele Frequency and Count.
- `pavs:caddPhred`, `pavs:revelScore`, `pavs:alphaMissense`, `pavs:spliceAI`: Pathogenicity scores.
- `pavs:clinvarId`, `pavs:clinvarSig`, `pavs:clinvarReview`: ClinVar details.
- `pavs:clingenAssertion`: ClinGen expert curation.

### `pavs:Gene`
Represents a human gene.
- `pavs:ncbiGeneId`: Link to NCBI Gene.
- `pavs:pli`, `pavs:loeuf`: Constraint metrics.
- `pavs:goProcess`, `pavs:goMolecularFunction`, `pavs:goCellularComponent`: Gene Ontology annotations.
- `pavs:mousePhenotypes`: Concatenated MGI mouse phenotype labels.
- `pavs:expressedIn`: Tissues with high expression (GTEx).
- `pavs:relatedDisease`: Link to associated OMIM diseases.

---

## Graphs

| Graph URI | Description |
|-----------|-------------|
| `http://pavs.kaust.edu.sa/graph/cases` | Saudi case records and variants. |
| `http://pavs.kaust.edu.sa/graph/genes` | Gene-level annotations and cross-references. |
| `http://pavs.kaust.edu.sa/graph/literature` | Literature-derived phenopackets (non-Saudi). |
| `http://pavs.kaust.edu.sa/graph/hpoa` | Global Disease-HPO associations (from HPOA). |
| `http://pavs.kaust.edu.sa/graph/hpo-ic` | HPO Information Content values and hierarchy labels. |
