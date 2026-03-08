# PAVS Data Model

The PAVS RDF data model is grounded in the [Semanticscience Integrated Ontology (SIO)](https://github.com/MaastrichtU-IDS/semanticscience). Every PAVS class is a subclass of a SIO class, and every PAVS property is a sub-property of a SIO property. The ontology is defined in `ontology/pavs.ttl`.

## Ontology File

- **`ontology/pavs.ttl`** — OWL ontology declaring all PAVS classes and properties with SIO superclasses/superproperties
- **`ontology/pavs_shapes.shex`** — ShEx shapes for validating generated RDF
- **`ontology/sio.owl`** — SIO ontology (imported by pavs.ttl)

## Class Hierarchy

| PAVS Class | SIO Superclass | SIO ID | Description |
|---|---|---|---|
| `pavs:Case` | `sio:patient` | SIO_000393 | Patient case with phenotypes, variants, demographics |
| `pavs:GenomicVariant` | `sio:genomic sequence variant` | SIO_001381 | Genomic variant observed in a patient |
| `pavs:GeneRecord` | `sio:gene` | SIO_010035 | Gene with constraint scores, GO, expression, diseases |
| `pavs:DiseaseHPOAssociation` | `sio:association` | SIO_000897 | Disease-to-HPO annotation from HPOA |
| `pavs:HPOTerm` | `sio:quality` | SIO_000005 | HPO phenotype term |

## Property Hierarchy

### Object Properties

All object properties are sub-properties of SIO object properties.

| PAVS Property | SIO Super-property | Used on | Target |
|---|---|---|---|
| `pavs:hasPhenotype` | `sio:has phenotype` (SIO_001279) | Case | HP term IRI |
| `pavs:hasExcludedPhenotype` | `sio:has attribute` (SIO_000008) | Case | HP term IRI |
| `pavs:hasVariant` | `sio:has attribute` (SIO_000008) | Case | GenomicVariant |
| `pavs:hasDisease` | `sio:has attribute` (SIO_000008) | Case | OMIM/MONDO IRI |
| `pavs:affectsGene` | `sio:refers to` (SIO_000628) | GenomicVariant | Gene IRI |
| `pavs:sex` | `sio:has attribute` (SIO_000008) | Case | NCIT sex IRI |
| `pavs:ancestry` | `sio:has attribute` (SIO_000008) | Case | HANCESTRO IRI |
| `pavs:geographicOrigin` | `sio:has attribute` (SIO_000008) | Case | GAZ IRI |
| `pavs:zygosity` | `sio:has attribute` (SIO_000008) | GenomicVariant | GENO IRI |
| `pavs:rsId` | `sio:refers to` (SIO_000628) | GenomicVariant | dbSNP IRI |
| `pavs:relatedDisease` | `sio:refers to` (SIO_000628) | GeneRecord | OMIM IRI |
| `pavs:disease` | `sio:refers to` (SIO_000628) | DiseaseHPOAssociation | Disease IRI |
| `pavs:hpoTerm` | `sio:refers to` (SIO_000628) | DiseaseHPOAssociation | HP IRI |
| `pavs:ncbiGeneId` | `sio:refers to` (SIO_000628) | GeneRecord | NCBI Gene IRI |

### Datatype Properties

All datatype properties are sub-properties of `sio:has value` (SIO_000300). Standard properties (`dc:identifier`, `rdfs:label`, `rdfs:seeAlso`, `rdfs:subClassOf`) are used alongside SIO-grounded properties.

| Property | Type | Used on | Description |
|---|---|---|---|
| `pavs:source` | string | Case, GenomicVariant | Source file name |
| `pavs:isSaudi` | boolean | Case, GenomicVariant | Saudi origin flag |
| `pavs:age` | string | Case | ISO 8601 duration |
| `pavs:diseaseLabel` | string | Case | Disease display name |
| `pavs:progressStatus` | string | Case | SOLVED/UNSOLVED/UNKNOWN |
| `pavs:pmid` | string | Case, GenomicVariant | PubMed ID |
| `pavs:consanguinity` | string | Case | Consanguinity label |
| `pavs:familyHistory` | string | Case | Family history label |
| `pavs:hgvsC` | string | GenomicVariant | HGVS c. notation |
| `pavs:hgvsP` | string | GenomicVariant | HGVS p. notation |
| `pavs:hgvsG` | string | GenomicVariant | HGVS g. notation |
| `pavs:acmgClass` | string | GenomicVariant | ACMG classification |
| `pavs:vcfChrom` | string | GenomicVariant | VCF chromosome |
| `pavs:vcfPos` | integer | GenomicVariant | VCF position |
| `pavs:vcfRef` | string | GenomicVariant | VCF ref allele |
| `pavs:vcfAlt` | string | GenomicVariant | VCF alt allele |
| `pavs:vepConsequence` | string | GenomicVariant | VEP consequence |
| `pavs:vepImpact` | string | GenomicVariant | VEP impact |
| `pavs:sift` | string | GenomicVariant | SIFT prediction |
| `pavs:polyphen` | string | GenomicVariant | PolyPhen-2 prediction |
| `pavs:gnomadAF` | double | GenomicVariant | gnomAD global AF |
| `pavs:gnomadAF_MID` | double | GenomicVariant | gnomAD Middle Eastern AF |
| `pavs:gnomadAF_SAS` | double | GenomicVariant | gnomAD South Asian AF |
| `pavs:caddPhred` | double | GenomicVariant | CADD phred score |
| `pavs:revelScore` | double | GenomicVariant | REVEL score |
| `pavs:alphaMissense` | string | GenomicVariant | AlphaMissense prediction |
| `pavs:spliceAI` | double | GenomicVariant | SpliceAI max delta |
| `pavs:saudiAF` | double | GenomicVariant | Saudi cohort AF |
| `pavs:saudiAC` | integer | GenomicVariant | Saudi cohort allele count |
| `pavs:clinvarId` | string | GenomicVariant | ClinVar allele ID |
| `pavs:clinvarSig` | string | GenomicVariant | ClinVar significance |
| `pavs:clinvarReview` | string | GenomicVariant | ClinVar review status |
| `pavs:clingenAssertion` | string | GenomicVariant | ClinGen assertion |
| `pavs:pli` | double | GeneRecord | pLI score |
| `pavs:loeuf` | double | GeneRecord | LOEUF score |
| `pavs:goProcess` | string | GeneRecord | GO biological process |
| `pavs:goMolecularFunction` | string | GeneRecord | GO molecular function |
| `pavs:goCellularComponent` | string | GeneRecord | GO cellular component |
| `pavs:mousePhenotypes` | string | GeneRecord | Mouse MP labels |
| `pavs:expressedIn` | string | GeneRecord | GTEx expression |
| `pavs:informationContent` | double | HPOTerm | HPO IC value |
| `pavs:evidence` | string | DiseaseHPOAssociation | Evidence code |
| `pavs:frequency` | string | DiseaseHPOAssociation | Frequency annotation |

## Named Graphs

| Graph URI | Contents | Source |
|---|---|---|
| `graph/cases` | Saudi PAVS cases | `rdf_output/cases.ttl` |
| `graph/genes` | Gene records + variant data | `rdf_output/genes.ttl` |
| `graph/hpoa` | Disease-HPO associations | `rdf_output/hpoa.ttl` |
| `graph/hpo-ic` | HPO IC + hierarchy + labels | `rdf_output/hpo_ic.ttl` |
| `graph/literature` | Literature phenopackets | `rdf_output/literature.ttl` |

## Validation

ShEx shapes in `ontology/pavs_shapes.shex` validate all generated RDF:

```bash
uv run python intake/validate_shex.py \
    --shapes ontology/pavs_shapes.shex \
    --rdf-dir rdf_output/
```

## ELK Reasoner Compatibility

The ontology uses only OWL 2 EL constructs (SubClassOf, SubObjectPropertyOf) and is compatible with the ELK reasoner. Load `ontology/pavs.ttl` with `ontology/sio.owl` in Protege and classify with ELK.

## Design Decisions

1. **Minimal SPARQL disruption**: All `pavs:` property URIs are unchanged in the instance data. The SIO grounding is declared in the ontology file only, so SPARQL queries and backend code need no changes.

2. **Class renames**: `pavs:Variant` was renamed to `pavs:GenomicVariant` and `pavs:Gene` to `pavs:GeneRecord` to better align with SIO superclass semantics.

3. **Standard properties preserved**: `dc:identifier`, `rdfs:label`, `rdfs:seeAlso`, and `rdfs:subClassOf` are used as-is alongside SIO-grounded properties.
