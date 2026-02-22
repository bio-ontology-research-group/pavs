import pandas as pd
import json
import os
import re
import requests
import time
import pickle
from rapidfuzz import process, utils, fuzz
from pyhpo import Ontology
import pronto
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configuration
GENE_INFO_PATH = "data/reference/Homo_sapiens.gene_info"
HPO_OBO_PATH = "ontology/hp.obo"
MONDO_OBO_PATH = "ontology/mondo.obo"
GENO_OBO_PATH = "ontology/geno.obo"
OMIM_HPOA_PATH = "data/reference/phenotype.hpoa"
MIM2GENE_PATH = "data/reference/mim2gene.txt"
GFF_PATH = "data/reference/GCF_000001405.40_GRCh38.p14_genomic.gff"
GEMINI_MODEL = "openai/gpt-oss-120b"
# Provider preferences for faster routing (in order of preference)
PROVIDER_PREFERENCES = ["Clarifai", "Chutes", "Google Vertex AI"]

# Valid ACMG evidence codes from ACMG/AMP 2015 guidelines
VALID_ACMG_CODES = {
    "PVS1",
    "PS1",
    "PS2",
    "PS3",
    "PS4",
    "PM1",
    "PM2",
    "PM3",
    "PM4",
    "PM5",
    "PM6",
    "PP1",
    "PP2",
    "PP3",
    "PP4",
    "PP5",
    "BA1",
    "BS1",
    "BS2",
    "BS3",
    "BS4",
    "BP1",
    "BP2",
    "BP3",
    "BP4",
    "BP5",
    "BP6",
    "BP7",
}


class NormalizationUtils:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(NormalizationUtils, cls).__new__(cls)
            cls._instance.initialized = False
        return cls._instance

    def initialize(self):
        if self.initialized:
            return
        print("Initializing Normalization Utils...")
        try:
            Ontology()
            self.hpo_terms = [{"id": t.id, "name": t.name} for t in Ontology]
            self.hpo_names = [t["name"] for t in self.hpo_terms]
            self.hpo_gene_map = {g.name: g for g in Ontology.genes}
        except:
            self.hpo_terms, self.hpo_names = [], []
            self.hpo_gene_map = {}
        self.mondo = {}
        self.omim_to_mondo: dict = {}
        if os.path.exists(MONDO_OBO_PATH):
            try:
                ms = pronto.Ontology(MONDO_OBO_PATH)
                for term in ms.terms():
                    if term.id.startswith("MONDO:"):
                        self.mondo[term.id] = term.name
                        for xref in term.xrefs:
                            if xref.id.startswith("OMIM:"):
                                self.omim_to_mondo[xref.id] = term.id
            except:
                pass
        self.mondo_list = [{"id": k, "name": v} for k, v in self.mondo.items()]
        self.mondo_names = [x["name"] for x in self.mondo_list]
        self.disease_map = {}
        if os.path.exists(OMIM_HPOA_PATH):
            try:
                df = pd.read_csv(
                    OMIM_HPOA_PATH,
                    sep="\t",
                    comment="#",
                    usecols=[0, 1],
                    names=["database_id", "disease_name"],
                )
                self.disease_map = dict(zip(df["database_id"], df["disease_name"]))
            except:
                pass
        self.disease_list = [{"id": k, "name": v} for k, v in self.disease_map.items()]
        self.disease_names = [x["name"] for x in self.disease_list]
        self.mim_info = {}
        if os.path.exists(MIM2GENE_PATH):
            try:
                df_mim = pd.read_csv(
                    MIM2GENE_PATH,
                    sep="\t",
                    comment="#",
                    header=None,
                    usecols=[0, 1],
                    names=["mim", "type"],
                )
                self.mim_info = dict(zip(df_mim["mim"].astype(str), df_mim["type"]))
            except:
                pass
        self.genes = {}
        cache_path = "data/gene_map_v2.pkl"
        if os.path.exists(cache_path):
            with open(cache_path, "rb") as f:
                self.genes = pickle.load(f)
        elif os.path.exists(GENE_INFO_PATH):
            try:
                for chunk in pd.read_csv(
                    GENE_INFO_PATH,
                    sep="\t",
                    usecols=["GeneID", "Symbol", "dbXrefs", "Synonyms"],
                    chunksize=100000,
                ):
                    for _, row in chunk.iterrows():
                        gid, sym, dbx = (
                            str(row["GeneID"]),
                            str(row["Symbol"]).upper(),
                            str(row["dbXrefs"]),
                        )
                        hgnc = re.search(r"HGNC:HGNC:(\d+)", dbx)
                        hgnc_id = f"HGNC:{hgnc.group(1)}" if hgnc else ""
                        data = {"hgnc": hgnc_id, "entrez": gid, "symbol": sym}
                        self.genes[sym] = data
                        syns = str(row["Synonyms"])
                        if syns != "-":
                            for s in syns.split("|"):
                                if s.upper() not in self.genes:
                                    self.genes[s.upper()] = data
                with open(cache_path, "wb") as f:
                    pickle.dump(self.genes, f)
            except:
                pass
        self.gene_coords = {}
        if os.path.exists(GFF_PATH):
            print("Indexing GFF...")
            try:
                with open(GFF_PATH, "r") as f:
                    for line in f:
                        if line.startswith("#"):
                            continue
                        p = line.split("\t")
                        if len(p) > 8 and p[2] == "gene":
                            m = re.search(r"Name=([^;]+)", p[8])
                            if m:
                                self.gene_coords[m.group(1).upper()] = True
            except:
                pass
        self.initialized = True
        self.graph_rag = None

    def enable_graph_rag(self):
        try:
            from normalization.graph_rag_normalization import GraphRAG

            print("Enabling GraphRAG...")
            self.graph_rag = GraphRAG()
        except ImportError as e:
            print(f"Could not import GraphRAG: {e}")

    def classify_omim(self, omim_id):
        oid = str(omim_id).replace("OMIM:", "")
        mtype = self.mim_info.get(oid, "unknown")
        if "phenotype" in mtype or "disorder" in mtype:
            return "disease"
        if "gene" in mtype:
            return "gene"
        return "unknown"

    def normalize_sex(self, text):
        m = {
            "male": ("NCIT:C20197", "Male"),
            "m": ("NCIT:C20197", "Male"),
            "female": ("NCIT:C16576", "Female"),
            "f": ("NCIT:C16576", "Female"),
        }
        return m.get(str(text).lower().strip(), ("NCIT:C17998", "Unknown"))

    def normalize_genotype(self, text):
        m = {
            "homozygous": ("GENO:0000136", "homozygous"),
            "heterozygous": ("GENO:0000135", "heterozygous"),
            "hemizygous": ("GENO:0000134", "hemizygous"),
        }
        t = str(text).lower()
        for k, v in m.items():
            if k[:3] in t:
                return v
        return ("GENO:0000137", "unspecified zygosity")

    def normalize_ancestry(self):
        return ("HANCESTRO:0852", "Middle Eastern")

    def normalize_country(self):
        return ("GAZ:00005279", "Saudi Arabia")

    def normalize_consanguinity(self, text):
        if not text or pd.isna(text):
            return ("", "")
        if any(x in str(text).lower() for x in ["yes", "positive", "consanguineous"]):
            return ("HP:0025820", "Offspring of consanguineous parents")
        return ("", "")

    def normalize_family_history(self, text):
        if not text or pd.isna(text):
            return ("", "")
        if any(x in str(text).lower() for x in ["yes", "positive"]):
            return ("HP:0032316", "Positive family history")
        return ("", "")

    def filter_acmg_evidence(self, evidence_str):
        """Filter evidence string to only include valid ACMG codes."""
        if not evidence_str or pd.isna(evidence_str):
            return ""

        # Split by common delimiters
        codes = re.split(r"[,;\s]+", str(evidence_str).strip())
        valid_codes = []

        for code in codes:
            code_clean = code.strip().upper()
            if code_clean in VALID_ACMG_CODES:
                valid_codes.append(code_clean)

        return ", ".join(valid_codes)

    def normalize_variant_classification(self, text):
        if not text or pd.isna(text):
            return ("", "")
        t = str(text).lower().strip()
        evidence = ""

        # Robust parsing for evidence
        # Case 5: "P  PVS1, PM2, PP3) " (missing opening paren)
        # Standard: "LP (PM2, PP1)"

        # Try finding the start of criteria (usually uppercase letters followed by numbers, e.g. PVS1, PM2)
        # But classification is also letters.
        # Split by '(' if present
        if "(" in text:
            parts = text.split("(", 1)
            classification = parts[0].strip()
            evidence = parts[1].replace(")", "").strip()
        elif ")" in text:
            # Case 5 scenario: "P PVS1...)"
            parts = text.split(")", 1)
            # This splits at the end. Not helpful.
            # Use regex to find the break between Class and Evidence (criteria starting with P, B, etc?)
            # Let's try splitting by the first occurrence of a criteria pattern?
            m = re.search(r"\s+(PVS1|PS\d|PM\d|PP\d|BA1|BS\d|BP\d)", text)
            if m:
                classification = text[: m.start()].strip()
                evidence = text[m.start() :].replace(")", "").strip()
            else:
                classification = text.replace(")", "").strip()
        else:
            classification = text.strip()

        t_class = classification.lower()
        normalized_class = classification.title()  # Default

        # Standardize Classification
        # Handle "Positive" and "Negative" from fawzan-variants
        if t_class in ["positive", "negative"]:
            # These are not variant classifications, they're test results
            # We'll handle this in the main script - return empty for now
            normalized_class = ""
        elif (
            any(x in t_class for x in ["pathogenic", "pathogenicity", "p "])
            or t_class == "p"
        ):
            if "likely" in t_class or "lp" in t_class.split():
                normalized_class = "Likely pathogenic"
            else:
                normalized_class = "Pathogenic"
        elif any(x in t_class for x in ["benign", "lb"]):
            if "likely" in t_class or "lb" in t_class.split():
                normalized_class = "Likely benign"
            else:
                normalized_class = "Benign"
        elif any(
            x in t_class for x in ["vus", "uncertain", "unknown", "ambiguous", "vous"]
        ):
            normalized_class = "Uncertain significance"
        elif t_class == "lp":
            normalized_class = "Likely pathogenic"
        elif t_class == "lb":
            normalized_class = "Likely benign"

        # Filter evidence to only valid ACMG codes
        evidence = self.filter_acmg_evidence(evidence)

        return (normalized_class, evidence)

    def validate_severity(self, hpo_id):
        if not hpo_id or not hpo_id.startswith("HP:"):
            return False
        try:
            # Check if subclass of Severity (HP:0012824)
            term = Ontology.get_hpo_object(hpo_id)
            severity_root = Ontology.get_hpo_object("HP:0012824")
            if term and severity_root and severity_root in term.superclasses:
                return True
        except:
            pass
        return False

    def validate_omim_against_text(self, omim_id, original_text):
        # Strict validation: Check if the numeric ID exists in the original text (unless explicit column)
        # This prevents LLM from hallucinating IDs not present in the input
        oid_num = omim_id.replace("OMIM:", "")
        if oid_num in original_text:
            return True
        return False

    def disambiguate_term(self, term, other_terms_text, gene_symbol):
        # Specific handler for ASD
        if term.upper() == "ASD":
            candidates = [
                {"id": "HP:0000729", "name": "Autism spectrum disorder"},
                {"id": "HP:0001631", "name": "Atrial septal defect"},
            ]

            # 1. Semantic Similarity with other terms (if available)
            # This is hard because we only have text here usually.
            # Simple keyword matching in context:
            context = str(other_terms_text).lower()
            if any(
                x in context for x in ["heart", "cardiac", "septal", "murmur", "defect"]
            ):
                return "Atrial septal defect"
            if any(
                x in context
                for x in ["development", "behavior", "speech", "social", "delay"]
            ):
                return "Autism spectrum disorder"

            # 2. Gene Context
            if gene_symbol and self.initialized:
                try:
                    # Check if gene is associated with either term in HPO
                    gene = self.hpo_gene_map.get(
                        gene_symbol.upper()
                    )  # Need to build/access this
                    if gene:
                        # This part requires a gene-to-hpo map which Ontology.genes provides
                        # but accessing it by symbol efficiently:
                        g_obj = next(
                            (
                                g
                                for g in Ontology.genes
                                if g.name == gene_symbol.upper()
                            ),
                            None,
                        )
                        if g_obj:
                            for term_id in [c["id"] for c in candidates]:
                                t_obj = Ontology.get_hpo_object(term_id)
                                if t_obj and t_obj in g_obj.hpo:
                                    return next(
                                        c["name"]
                                        for c in candidates
                                        if c["id"] == term_id
                                    )
                except:
                    pass

            # Default fallback (most common in this dataset?)
            # PMC6562004 has lots of developmental delay, so ASD -> Autism is likely default
            return "Autism spectrum disorder"

        return term

    def map_gene(self, text):
        return self.genes.get(str(text).strip().upper())

    def validate_variant(self, gene_sym):
        return self.gene_coords.get(str(gene_sym).upper(), False)

    def get_hpo_candidates(self, query):
        if not self.hpo_names:
            return []

        # Use multiple scorers and combine results
        # token_set_ratio is better for partial matches
        res1 = process.extract(
            query, self.hpo_names, limit=15, scorer=fuzz.token_set_ratio
        )
        # WRatio for balanced scoring
        res2 = process.extract(query, self.hpo_names, limit=15, scorer=fuzz.WRatio)

        # Combine and deduplicate
        seen = set()
        combined = []
        for name, score, idx in list(res1) + list(res2):
            if idx not in seen:
                seen.add(idx)
                combined.append((name, score, idx))

        # Sort by score and take top 10
        combined.sort(key=lambda x: x[1], reverse=True)

        return [
            {"id": self.hpo_terms[idx]["id"], "name": name, "score": score}
            for name, score, idx in combined[:10]
        ]

    def lookup_hpo_id_by_label(self, label):
        """Look up HPO ID from label/name. Returns None if not found."""
        if not label or label in ["None", "null", ""]:
            return None

        label_lower = label.lower().strip()

        # First try exact match
        for term in self.hpo_terms:
            if term["name"].lower() == label_lower:
                return term["id"]

        # Then try fuzzy match (high threshold)
        matches = process.extract(label, self.hpo_names, limit=1, scorer=fuzz.ratio)
        if matches and matches[0][1] >= 95:  # 95% similarity required
            idx = matches[0][2]
            return self.hpo_terms[idx]["id"]

        return None

    def lookup_omim_id_by_label(self, label):
        """Look up OMIM ID from label/name. Returns None if not found."""
        if not label or label in ["None", "null", ""]:
            return None

        # Search in disease_map
        label_lower = label.lower().strip()
        for omim_id, info in self.disease_map.items():
            if info.get("name", "").lower() == label_lower:
                return f"OMIM:{omim_id}"

        return None

    def lookup_mondo_id_by_label(self, label):
        """Look up MONDO ID from label/name. Returns None if not found."""
        if not label or label in ["None", "null", ""]:
            return None

        # Search in mondo_diseases
        label_lower = label.lower().strip()
        for disease in self.mondo_diseases:
            if disease["name"].lower() == label_lower:
                return disease["id"]

        return None

    def get_disease_candidates(self, query):
        all_d, all_n = (
            self.disease_list + self.mondo_list,
            self.disease_names + self.mondo_names,
        )
        res = process.extract(query, all_n, limit=5, scorer=fuzz.token_sort_ratio)
        return [{"id": all_d[idx]["id"], "name": name} for name, score, idx in res]

    def llm_normalize(self, terms, context="", max_workers=10):
        # GraphRAG path
        if getattr(self, "graph_rag", None):
            results = {}
            for term in terms:
                if not term or term == "nan":
                    continue
                try:
                    res = self.graph_rag.normalize_term(term, context_text=context)
                    if res:
                        # GraphRAG returns single object: {"hpo": ..., "omim": ...}
                        # llm_normalize expects {term: {"hpo": ..., ...}}
                        results[term] = res
                except Exception as e:
                    print(f"GraphRAG error for {term}: {e}")
            return results

        api_key = os.getenv("OPENROUTER_API_KEY")
        if not api_key or not terms:
            return {}
        data_in = [
            {
                "term": t,
                "hpo_c": self.get_hpo_candidates(t),
                "disease_c": self.get_disease_candidates(t),
            }
            for t in terms
            if t and t != "nan"
        ]
        if not data_in:
            return {}

        prompt = f"""You are a clinical geneticist. Map clinical terms to standard ontology labels.

Context (patient genes): {context}

CRITICAL RULES:
1. You MUST select ONLY from the provided candidate labels in hpo_c, disease_c
2. Output ONLY the LABEL/NAME, NEVER invent or output IDs (HP:, OMIM:, MONDO:)
3. A single term can map to MULTIPLE phenotypes if it contains words like "and", "with", or describes multiple features
4. Extract ALL phenotypes mentioned in each term - do not just pick one
5. Detect negation (no, not, without, absent) and mark excluded: true
6. If no good match in candidates, return "None"

Output JSON structure where each KEY is the original term:
{{
  "term": {{
    "hpo_labels": ["exact label from hpo_c", ...],
    "hpo_excluded": false,
    "omim_labels": ["label from disease_c", ...] or "None",
    "mondo_labels": ["label from disease_c", ...] or "None"
  }}
}}

Examples:
- Single phenotype: "feeding difficulties" → {{"hpo_labels": ["Feeding difficulties"]}}
- Multiple phenotypes: "scaly skin and hydronephrosis with reflux" → {{"hpo_labels": ["Scaling skin", "Hydronephrosis", "Vesicoureteral reflux"]}}
- Negation: "no seizures" → {{"hpo_labels": ["Seizure"], "hpo_excluded": true}}

IMPORTANT: Look for conjunctions (and, with, or) as indicators of MULTIPLE phenotypes in one term. Extract ALL of them.

FORBIDDEN: DO NOT output any IDs (HP:XXXXXXX, OMIM:XXXXXX, MONDO:XXXXXXX). Only output labels that exist in the candidates.

Terms with candidates:
{json.dumps(data_in)}
"""
        # Retry up to 2 times for invalid JSON or API errors
        result = None
        for attempt in range(3):
            try:
                request_data = {
                    "model": GEMINI_MODEL,
                    "messages": [{"role": "user", "content": prompt}],
                }
                if PROVIDER_PREFERENCES:
                    request_data["provider"] = {"order": PROVIDER_PREFERENCES}

                resp = requests.post(
                    "https://openrouter.ai/api/v1/chat/completions",
                    headers={
                        "Authorization": f"Bearer {api_key}",
                        "Content-Type": "application/json",
                    },
                    json=request_data,
                    timeout=60,
                )
                if resp.status_code != 200:
                    if attempt < 2:
                        print(
                            f"Variant API Error {resp.status_code}, retrying... (attempt {attempt + 1}/3)"
                        )
                        time.sleep(1)
                        continue
                    else:
                        print(f"Variant API Error {resp.status_code}: {resp.text}")
                        return ""

                raw_content = resp.json()["choices"][0]["message"]["content"]
                content = re.sub(r"```json\n?|\n?```", "", raw_content).strip()

                # Try to parse JSON
                parsed = json.loads(content)
                result = parsed
                break  # Success - exit retry loop

            except json.JSONDecodeError as e:
                if attempt < 2:
                    print(
                        f"Variant JSON decode error, retrying... (attempt {attempt + 1}/3)"
                    )
                    time.sleep(1)
                    continue
                else:
                    print(f"Variant JSON decode error after 3 attempts: {e}")
                    return ""
            except Exception as e:
                if attempt < 2:
                    print(
                        f"Variant extraction error, retrying... (attempt {attempt + 1}/3): {str(e)[:50]}"
                    )
                    time.sleep(1)
                    continue
                else:
                    print(f"Variant extraction error: {e}")
                    return ""

        if not result:
            return {}

        # VALIDATION STEP: Filter out hallucinated phenotypes
        # Ask LLM to verify each HPO term is actually mentioned in the source
        validated_result = self._validate_llm_results(result, terms, context, api_key)

        final_result = {}
        for k, v in validated_result.items():
            if isinstance(v.get("hpo"), list):
                final_result[k] = v
            else:
                final_result[k] = v
        return final_result

    def _validate_llm_results(self, llm_results, original_terms, context, api_key):
        """
        Validation step: Filter out hallucinated phenotypes.
        For each extracted HPO term, verify it's actually mentioned in the source text.
        """
        if not llm_results:
            return {}

        # Build validation prompt - validate LABELS not IDs
        validation_data = []
        for term, data in llm_results.items():
            # New format: hpo_labels
            hpo_labels = data.get("hpo_labels", data.get("hpo"))
            if hpo_labels and hpo_labels != "None":
                if isinstance(hpo_labels, list):
                    for label in hpo_labels:
                        if isinstance(label, dict):
                            label = label.get("label")
                        if label and label != "None":
                            validation_data.append(
                                {"original_term": term, "hpo_label": label}
                            )
                elif isinstance(hpo_labels, str):
                    validation_data.append(
                        {"original_term": term, "hpo_label": hpo_labels}
                    )
                elif isinstance(hpo_labels, dict):
                    # Old format support
                    label = hpo_labels.get("label")
                    if label and label != "None":
                        validation_data.append(
                            {"original_term": term, "hpo_label": label}
                        )

        if not validation_data:
            return llm_results

        validation_prompt = f"""You are validating phenotype extraction. Check if each HPO label is ACTUALLY mentioned in the original clinical text.

Original clinical terms: {original_terms}

For each HPO label below, return "true" ONLY if the term or a DIRECT synonym appears in the original text.
Return "false" if the term is:
- A more specific subtype not explicitly mentioned
- A related but different phenotype
- An interpretation or inference not directly stated
- About a different body system entirely

STRICT matching examples:
- Text: "feeding difficulties" → "Feeding difficulties" = TRUE, "Difficulty standing" = FALSE
- Text: "noisy breathing" → "Noisy breathing" = TRUE, "Habitual mouth breathing" = FALSE
- Text: "ventricular arrhythmia" → "Ventricular arrhythmia" = TRUE, "Paroxysmal ventricular tachycardia" = FALSE
- Text: "cardiac arrest" → "Cardiac arrest" = TRUE, "Sudden cardiac death" = FALSE

Return JSON array:
[
  {{"original_term": "text", "hpo_label": "label", "is_valid": true}},
  {{"original_term": "text", "hpo_label": "label", "is_valid": false}}
]

Validate these:
{json.dumps(validation_data)}
"""

        # Retry validation up to 2 times
        validation_results = None
        for attempt in range(3):
            try:
                request_data = {
                    "model": GEMINI_MODEL,
                    "messages": [{"role": "user", "content": validation_prompt}],
                }
                if PROVIDER_PREFERENCES:
                    request_data["provider"] = {"order": PROVIDER_PREFERENCES}

                resp = requests.post(
                    "https://openrouter.ai/api/v1/chat/completions",
                    headers={
                        "Authorization": f"Bearer {api_key}",
                        "Content-Type": "application/json",
                    },
                    json=request_data,
                    timeout=60,
                )

                if resp.status_code != 200:
                    if attempt < 2:
                        print(
                            f"Validation API Error {resp.status_code}, retrying... (attempt {attempt + 1}/3)"
                        )
                        time.sleep(1)
                        continue
                    else:
                        print(f"Validation API Error: {resp.text}")
                        return llm_results  # Return original if validation fails

                raw_validation = resp.json()["choices"][0]["message"]["content"]
                if not raw_validation or not raw_validation.strip():
                    if attempt < 2:
                        print(
                            f"Empty validation response, retrying... (attempt {attempt + 1}/3)"
                        )
                        time.sleep(1)
                        continue
                    else:
                        print(
                            "Empty validation response after retries, skipping validation"
                        )
                        return llm_results

                validation_content = re.sub(
                    r"```json\n?|\n?```",
                    "",
                    raw_validation,
                ).strip()

                if not validation_content:
                    if attempt < 2:
                        print(
                            f"Empty validation content, retrying... (attempt {attempt + 1}/3)"
                        )
                        time.sleep(1)
                        continue
                    else:
                        print(
                            "Empty validation content after retries, skipping validation"
                        )
                        return llm_results

                validation_results = json.loads(validation_content)

                # Handle different response formats
                if isinstance(validation_results, dict):
                    # If it's a dict, try to get the array from common keys
                    validation_results = validation_results.get(
                        "results", validation_results.get("validations", [])
                    )

                if not isinstance(validation_results, list):
                    if attempt < 2:
                        print(
                            f"Unexpected validation format, retrying... (attempt {attempt + 1}/3)"
                        )
                        time.sleep(1)
                        continue
                    else:
                        print(
                            f"Unexpected validation response format: {type(validation_results)}"
                        )
                        return llm_results

                # Filter out invalid results - use LABELS not IDs
                valid_hpo_labels = set()
                for v in validation_results:
                    if isinstance(v, dict) and v.get("is_valid"):
                        valid_hpo_labels.add(v.get("hpo_label"))

                # Rebuild llm_results with only valid HPO labels
                filtered_results = {}
                for term, data in llm_results.items():
                    filtered_data = data.copy()

                    # New format: hpo_labels
                    hpo_labels = data.get("hpo_labels", data.get("hpo"))

                    if isinstance(hpo_labels, list):
                        # Filter to only valid labels
                        filtered_labels = []
                        for label in hpo_labels:
                            if isinstance(label, dict):
                                label_text = label.get("label")
                            else:
                                label_text = label

                            if label_text in valid_hpo_labels:
                                filtered_labels.append(label_text)

                        filtered_data["hpo_labels"] = (
                            filtered_labels if filtered_labels else "None"
                        )
                    elif isinstance(hpo_labels, str):
                        filtered_data["hpo_labels"] = (
                            hpo_labels if hpo_labels in valid_hpo_labels else "None"
                        )
                    elif isinstance(hpo_labels, dict):
                        # Old format support
                        label = hpo_labels.get("label")
                        if label in valid_hpo_labels:
                            filtered_data["hpo_labels"] = [label]
                        else:
                            filtered_data["hpo_labels"] = "None"
                    else:
                        filtered_data["hpo_labels"] = "None"

                    filtered_results[term] = filtered_data

                return filtered_results

            except (json.JSONDecodeError, KeyError, Exception) as e:
                if attempt < 2:
                    print(
                        f"Validation error, retrying... (attempt {attempt + 1}/3): {str(e)[:50]}"
                    )
                    time.sleep(1)
                    continue
                else:
                    print(f"Validation Error after retries: {e}")
                    return llm_results

        return llm_results  # Fallback if all retries failed

    def extract_variant_hgvs(self, text, gene=""):
        api_key = os.getenv("OPENROUTER_API_KEY")
        if not api_key or not text or text == "nan":
            return ""
        prompt = f"Extract standard HGVS (e.g. Gene:c.123A>G or Gene:p.Ala123Gly) from text: '{text}'. Gene context: '{gene}'. Return JSON: {{\"variants\": [\"...\"]}}. Only include variants EXPLICITLY in the text. Do not invent genes."
        try:
            request_data = {
                "model": GEMINI_MODEL,
                "messages": [{"role": "user", "content": prompt}],
                "response_format": {"type": "json_object"},
            }
            if PROVIDER_PREFERENCES:
                request_data["provider"] = {"order": PROVIDER_PREFERENCES}

            resp = requests.post(
                "https://openrouter.ai/api/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {api_key}",
                    "Content-Type": "application/json",
                },
                json=request_data,
                timeout=60,
            )
            content = re.sub(
                r"```json\n?|\n?```",
                "",
                resp.json()["choices"][0]["message"]["content"],
            ).strip()

            raw_variants = json.loads(content).get("variants", [])
            valid_variants = []

            for v in raw_variants:
                # Validation: Gene part must match context or be in text
                if ":" in v:
                    v_gene = v.split(":")[0].strip().upper()
                    # Check 1: Matches provided context
                    if gene and v_gene == str(gene).strip().upper():
                        valid_variants.append(v)
                        continue

                    # Check 2: Exists in original text (strict check)
                    if v_gene in str(text).upper():
                        valid_variants.append(v)
                        continue

                    # Check 3: Is it an alias of the context?
                    if gene and self.initialized:
                        mapped_context = self.map_gene(gene)
                        mapped_v = self.map_gene(v_gene)
                        if (
                            mapped_context
                            and mapped_v
                            and mapped_context["hgnc"] == mapped_v["hgnc"]
                        ):
                            valid_variants.append(v)
                            continue
                else:
                    # No gene part? Only valid if we can prepend context
                    if gene:
                        valid_variants.append(f"{gene}:{v}")

            return "|".join(valid_variants)
        except Exception as e:
            print(f"Extract Variant Error: {e}")
            return text


def get_utils():
    u = NormalizationUtils()
    u.initialize()
    return u
