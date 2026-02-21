import json
import os
import glob
from sqlalchemy.ext.asyncio import AsyncSession
from .models import Phenopacket
from sqlalchemy.future import select

def extract_gene(data):
    for interpretation in data.get('interpretations', []):
        diagnosis = interpretation.get('diagnosis')
        if not diagnosis: continue
        for genomic_interp in diagnosis.get('genomicInterpretations', []):
            variant_interp = genomic_interp.get('variantInterpretation')
            if not variant_interp: continue
            variation_desc = variant_interp.get('variationDescriptor')
            if not variation_desc: continue
            gene_context = variation_desc.get('geneContext')
            if gene_context and gene_context.get('symbol'):
                return gene_context.get('symbol')
    return None

def extract_disease(data):
    diseases = data.get('diseases', [])
    if diseases:
        d = diseases[0].get('term')
        if d:
            return d.get('id'), d.get('label')
    for interpretation in data.get('interpretations', []):
        diagnosis = interpretation.get('diagnosis')
        if diagnosis and diagnosis.get('disease'):
            d = diagnosis.get('disease')
            return d.get('id'), d.get('label')
    return None, None

def extract_phenotypes(data):
    return [pf['type']['id'] for pf in data.get('phenotypicFeatures', []) if 'type' in pf]

# Global caches for faster resolution
phenopacket_cache = []
disease_to_phenotypes = {}  # Map disease ID -> list of HPO IDs

async def load_data(db: AsyncSession):
    global disease_to_phenotypes
    # Clear existing records to allow "reloading"
    from sqlalchemy import delete
    await db.execute(delete(Phenopacket))
    await db.commit()
    print("Database cleared for reload.")

    # Load disease name mapping and phenotypes from phenotype.hpoa
    disease_names = {}
    disease_to_phenotypes = {}
    hpoa_path = "data/phenotype.hpoa"
    if os.path.exists(hpoa_path):
        print(f"Loading disease metadata from {hpoa_path}...")
        with open(hpoa_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                parts = line.split('\t')
                if len(parts) >= 4:
                    did = parts[0].strip()
                    dname = parts[1].strip()
                    hpo_id = parts[3].strip()
                    
                    disease_names[did] = dname
                    if did not in disease_to_phenotypes:
                        disease_to_phenotypes[did] = []
                    if hpo_id.startswith("HP:"):
                        disease_to_phenotypes[did].append(hpo_id)
                    
                    # Also store without prefix for easier lookup
                    if ':' in did:
                        pure_id = did.split(':')[1]
                        disease_names[pure_id] = dname
                        if pure_id not in disease_to_phenotypes:
                            disease_to_phenotypes[pure_id] = []
                        disease_to_phenotypes[pure_id].append(hpo_id)

    # Supplemental names from Mondo if possible
    mondo_path = "ontology/mondo.obo"
    if os.path.exists(mondo_path):
        print(f"Loading supplemental names from {mondo_path}...")
        current_id = None
        current_name = None
        with open(mondo_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("id: "):
                    current_id = line[4:]
                elif line.startswith("name: "):
                    current_name = line[6:]
                    if current_id: disease_names[current_id] = current_name
                elif line.startswith("xref: OMIM:"):
                    # xref: OMIM:100100 {source="MONDO:equivalentTo"}
                    parts = line.split(' ')
                    if len(parts) >= 2:
                        omim_id = parts[1]
                        if current_name:
                            disease_names[omim_id] = current_name
                            if ':' in omim_id:
                                disease_names[omim_id.split(':')[1]] = current_name
    
    # Load gene-to-disease associations
    gene_to_dis = {}
    g2d_path = "data/genes_to_disease.txt"
    if os.path.exists(g2d_path):
        print(f"Loading gene-to-disease associations from {g2d_path}...")
        with open(g2d_path, 'r') as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    gsym = parts[1]
                    disid = parts[3]
                    if gsym not in gene_to_dis:
                        gene_to_dis[gsym] = []
                    gene_to_dis[gsym].append(disid)

    # Load gene ID mapping from Homo_sapiens.gene_info
    gene_id_to_symbol = {}
    gene_info_path = "data/Homo_sapiens.gene_info"
    if os.path.exists(gene_info_path):
        print(f"Loading gene info from {gene_info_path}...")
        with open(gene_info_path, 'r') as f:
            header = f.readline()
            for line in f:
                parts = line.split('\t')
                if len(parts) >= 3:
                    gene_id_to_symbol[parts[1]] = parts[2]

    def clean_gene_symbol(symbol):
        if not symbol: return None
        symbol_str = str(symbol).strip()
        
        # If it's a pure number, try to map from Entrez ID
        if symbol_str.isdigit():
            return gene_id_to_symbol.get(symbol_str, symbol_str)
            
        # Handle genomic ranges and coordinates
        # 1. Underscore range: 103645506_106261806
        if '_' in symbol_str and any(c.isdigit() for c in symbol_str):
            return "Genomic Range"
        # 2. Hyphen range: 86218116-86872832
        if '-' in symbol_str and any(c.isdigit() for c in symbol_str):
            # Check if it looks like a range (starts/ends with digit or is mostly digits)
            parts = symbol_str.split('-')
            if all(p.isdigit() for p in parts if p):
                 return "Genomic Range"
        
        # Handle weird symbols
        if symbol_str in ["=", "?", "."] or symbol_str.startswith("?"):
            return "Unknown"
            
        return symbol_str

    def extract_gene_cleaned(data):
        for interpretation in data.get('interpretations', []):
            diagnosis = interpretation.get('diagnosis')
            if not diagnosis: continue
            for genomic_interp in diagnosis.get('genomicInterpretations', []):
                variant_interp = genomic_interp.get('variantInterpretation')
                if not variant_interp: continue
                variation_desc = variant_interp.get('variationDescriptor')
                if not variation_desc: continue
                gene_context = variation_desc.get('geneContext')
                if gene_context and gene_context.get('symbol'):
                    return clean_gene_symbol(gene_context.get('symbol'))
        return None

    def resolve_disease(gene, dis_id, dis_name):
        is_suggested = 0
        
        # Normalize dis_id: if it's just a number, assume OMIM
        if dis_id and str(dis_id).isdigit():
            dis_id = f"OMIM:{dis_id}"

        # If we have a disease name from packet, use it
        if dis_name and str(dis_name).lower() not in ["unknown", "unknown disease"]:
            return dis_id, dis_name, 0
        
        # If we have id but no name (or name is "unknown"), look up name
        if dis_id:
            if dis_id in disease_names:
                return dis_id, disease_names[dis_id], 0
            # Try without prefix
            if ':' in str(dis_id):
                pure_id = str(dis_id).split(':')[1]
                if pure_id in disease_names:
                    return dis_id, disease_names[pure_id], 0
        
        # If both missing or unknown, try suggest ALL from gene
        if gene and gene in gene_to_dis:
            suggested_ids = gene_to_dis[gene]
            # Get all names, fallback to ID if name not found
            names = [disease_names.get(sid, sid) for sid in suggested_ids]
            
            return "|".join(suggested_ids), "|".join(names), 1
            
        return dis_id or "UNKNOWN", dis_name or "Unknown Disease", 0

    # Load Generated Phenopackets (ONLY source)
    generated_dir = "phenopackets/generated"
    if os.path.exists(generated_dir):
        print(f"Indexing Generated data from {generated_dir}...")
        json_files = glob.glob(os.path.join(generated_dir, "**/*.json"), recursive=True)
        print(f"Found {len(json_files)} generated files.")
        for file_path in json_files:
            with open(file_path, 'r') as f:
                try:
                    item = json.load(f)
                except:
                    continue
                gene = extract_gene_cleaned(item)
                dis_id, dis_name = extract_disease(item)
                
                dis_id, dis_name, is_suggested = resolve_disease(gene, dis_id, dis_name)
                phenos = extract_phenotypes(item)
                
                # Get clinical phenotypes from hpoa for all resolved disease IDs
                all_disease_phenos = set()
                for did in str(dis_id).split('|'):
                    d_phenos = disease_to_phenotypes.get(did, [])
                    if not d_phenos and ':' in str(did):
                        d_phenos = disease_to_phenotypes.get(str(did).split(':')[1], [])
                    for dp in d_phenos:
                        all_disease_phenos.add(dp)
                
                db_item = Phenopacket(
                    external_id=item.get('id'),
                    gene_symbol=gene,
                    disease_id=dis_id,
                    disease_name=dis_name,
                    disease_is_suggested=is_suggested,
                    phenotypes=",".join(phenos),
                    disease_phenotypes=",".join(list(all_disease_phenos)),
                    source="Generated",
                    file_path=file_path,
                    raw_json=item
                )
                db.add(db_item)
        await db.commit()
    print("All data indexed successfully.")