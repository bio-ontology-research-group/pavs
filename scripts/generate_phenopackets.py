import pandas as pd
import json
import os
import datetime
import re

def parse_hpo_obo(obo_path):
    hpo_id_to_name = {}
    if not os.path.exists(obo_path): return {}
    with open(obo_path, 'r') as f:
        current_id = None
        for line in f:
            line = line.strip()
            if line.startswith('id: HP:'):
                current_id = line.split('id: ')[1].split(' ! ')[0]
            elif line.startswith('name: '):
                hpo_id_to_name[current_id] = line.split('name: ')[1]
    return hpo_id_to_name

def create_phenopacket(row, hpo_labels):
    pp_id = row['Master_ID']
    subject_id = str(row['Original_ID'])
    
    # Subject (Individual) v2.0
    subject = {
        "id": subject_id,
        "sex": row['Sex'] if row['Sex'] in ["MALE", "FEMALE", "UNKNOWN_SEX"] else "UNKNOWN_SEX",
        "taxonomy": {
            "id": "NCBITaxon:9606",
            "label": "Homo sapiens"
        }
    }
    
    # Phenotypic Features
    phenotypes = []
    if pd.notna(row['Phenotypes_HPO']) and row['Phenotypes_HPO']:
        ids = str(row['Phenotypes_HPO']).split('|')
        for h_id in ids:
            if h_id.startswith('HP:'):
                phenotypes.append({
                    "type": {
                        "id": h_id,
                        "label": hpo_labels.get(h_id, "Unknown Phenotype")
                    }
                })
    
    # Ancestry and Location as phenotypic features
    phenotypes.append({
        "type": {"id": "HANCESTRO:0014", "label": "Middle Eastern, North African or Persian"}
    })
    phenotypes.append({
        "type": {"id": "GAZ:00000570", "label": "Saudi Arabia"}
    })
            
    # Diseases
    diseases = []
    if pd.notna(row['Diseases']) and row['Diseases']:
        d_ids = re.findall(r'(OMIM:\d+|ORPHA:\d+)', str(row['Diseases']))
        for d_id in d_ids:
            diseases.append({
                "term": {"id": d_id, "label": "Unknown Disease"}
            })

    # Interpretations (Variants) v2.0
    interpretations = []
    variants_summary = []
    if pd.notna(row['Variants_HGVS']) and row['Variants_HGVS']:
        v_strings = str(row['Variants_HGVS']).split('|')
        genomic_interpretations = []
        
        zygosity_id = row['Zygosity_GENO'] if pd.notna(row['Zygosity_GENO']) and row['Zygosity_GENO'] else "GENO:0000133"
        zygosity_label = {"GENO:0000135": "heterozygous", "GENO:0000136": "homozygous", "GENO:0000606": "hemizygous"}.get(zygosity_id, "unknown zygosity")

        for v_str in v_strings:
            v_str = v_str.strip()
            # Try to extract Gene, Transcript, and HGVS.c
            # Format 1: Gene:Transcript:exon:c.xxx:p.xxx
            # Format 2: Gene:c.xxx
            # Format 3: Transcript:c.xxx
            
            gene_symbol = "Unknown"
            hgvs_value = v_str
            
            # Extract Gene Symbol (usually first part before colon)
            gene_match = re.search(r'^([^:]+)', v_str)
            if gene_match:
                gene_symbol = gene_match.group(1).strip()

            # Try to find a transcript (NM_...) and the c. change
            transcript_match = re.search(r'(NM_\d+(\.\d+)?)', v_str)
            c_change_match = re.search(r'(c\.[^: ]+)', v_str)
            
            if transcript_match and c_change_match:
                hgvs_value = f"{transcript_match.group(1)}:{c_change_match.group(1)}"
            elif c_change_match:
                # If no transcript, use Gene:c.xxx
                hgvs_value = f"{gene_symbol}:{c_change_match.group(1)}"
            
            # Clean up hgvs_value (remove trailing protein or notes)
            hgvs_value = hgvs_value.split(' ')[0].split('(')[0].rstrip(':')

            genomic_interpretations.append({
                "subjectOrBiosampleId": subject_id,
                "interpretationStatus": "CAUSATIVE",
                "variantInterpretation": {
                    "variationDescriptor": {
                        "id": f"var_{pp_id}_{gene_symbol}_{len(genomic_interpretations)}",
                        "geneContext": {"symbol": gene_symbol},
                        "expressions": [{"syntax": "hgvs.c", "value": hgvs_value}],
                        "moleculeContext": "genomic",
                        "allelicState": {"id": zygosity_id, "label": zygosity_label}
                    }
                }
            })
            variants_summary.append({'gene': gene_symbol, 'hgvs': hgvs_value})
        
        if genomic_interpretations:
            interpretations.append({
                "id": f"int_{pp_id}",
                "progressStatus": "SOLVED",
                "diagnosis": {
                    "genomicInterpretations": genomic_interpretations
                }
            })

    # MetaData v2.0
    metadata = {
        "created": datetime.datetime.utcnow().isoformat() + "Z",
        "createdBy": "Robert Hoehndorf",
        "resources": [
            {"id": "hp", "name": "human phenotype ontology", "url": "http://purl.obolibrary.org/obo/hp.owl", "namespacePrefix": "HP"},
            {"id": "omim", "name": "An Online Catalog of Human Genes and Genetic Disorders", "url": "https://www.omim.org", "namespacePrefix": "OMIM"},
            {"id": "hancestro", "name": "Human Ancestry Ontology", "url": "http://purl.obolibrary.org/obo/hancestro.owl", "namespacePrefix": "HANCESTRO"},
            {"id": "gaz", "name": "Gazetteer", "url": "http://purl.obolibrary.org/obo/gaz.owl", "namespacePrefix": "GAZ"},
            {"id": "geno", "name": "Genotype Ontology", "url": "http://purl.obolibrary.org/obo/geno.owl", "namespacePrefix": "GENO"},
            {"id": "hgnc", "name": "HUGO Gene Nomenclature Committee", "url": "https://www.genenames.org", "namespacePrefix": "HGNC"}
        ],
        "phenopacketSchemaVersion": "2.0.2",
        "externalReferences": [
            {"id": "ORCID:0000-0001-8149-5890", "description": "Robert Hoehndorf"}
        ]
    }

    phenopacket = {
        "id": pp_id,
        "subject": subject,
        "phenotypicFeatures": phenotypes,
        "diseases": diseases,
        "interpretations": interpretations,
        "metaData": metadata
    }
    return phenopacket, variants_summary

def main():
    hpo_labels = parse_hpo_obo('ontology/hp.obo')
    df = pd.read_csv('master_data_final.tsv', sep='\t')
    output_dir = 'phenopackets/generated'
    os.makedirs(output_dir, exist_ok=True)
    summary_data = []
    
    for _, row in df.iterrows():
        pp, variants = create_phenopacket(row, hpo_labels)
        file_path = os.path.join(output_dir, f"{row['Master_ID']}.json")
        with open(file_path, 'w') as f:
            json.dump(pp, f, indent=2)
            
        disease_id = pp['diseases'][0]['term']['id'] if pp['diseases'] else ""
        gene = variants[0]['gene'] if variants else ""
        v1 = variants[0]['hgvs'] if variants else ""
        v2 = variants[1]['hgvs'] if len(variants) > 1 else ""
        
        summary_data.append({
            'patient_id': row['Master_ID'],
            'cohort': row['Source'],
            'disease_id': disease_id,
            'disease': "",
            'gene': gene,
            'allele_1': v1,
            'allele_2': v2,
            'PMID': "",
            'filename': f"generated/{row['Master_ID']}.json"
        })
        
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv('phenopackets/generated/phenopacket_store.summary.tsv', sep='\t', index=False)
    print(f"Updated {len(df)} phenopackets to v2.0 in {output_dir}")

if __name__ == "__main__":
    main()