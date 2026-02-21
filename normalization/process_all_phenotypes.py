import pandas as pd
import re
import os
import json
from normalization.normalization_utils import get_utils
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

utils = get_utils()
OUTPUT_DIR = "data/processed"
os.makedirs(OUTPUT_DIR, exist_ok=True)

COLUMNS = [
    "id",
    "source_file",
    "sex_id",
    "sex_label",
    "age",
    "ancestry_id",
    "ancestry_label",
    "country_id",
    "country_label",
    "consanguinity_id",
    "consanguinity_label",
    "family_history_id",
    "family_history_label",
    "hpo_ids",
    "hpo_labels",
    "hpo_excluded_ids",
    "hpo_excluded_labels",
    "hpo_severity_ids",
    "hpo_severity_labels",
    "disease_omim_gene_ids",
    "disease_omim_phenotype_ids",
    "disease_omim_labels",
    "disease_omim_excluded_ids",
    "disease_mondo_ids",
    "disease_mondo_labels",
    "disease_mondo_excluded_ids",
    "gene_hgnc_ids",
    "gene_symbols",
    "variant_hgvs",
    "variant_validated",
    "variant_classification",
    "variant_evidence",
    "publication_id",
    "genotype_id",
    "genotype_label",
    "solved",
]


def parse_iso_age(age_str):
    if pd.isna(age_str):
        return ""
    s = str(age_str).lower().strip()
    if s == "?" or not s or s == "nan":
        return ""
    y = re.search(r"(\d+)\s*(y|year)", s)
    m = re.search(r"(\d+)\s*(m|month)", s)
    if not y and not m and s.isdigit():
        return f"P{s}Y"
    res = "P"
    if y:
        res += f"{y.group(1)}Y"
    if m:
        res += f"{m.group(1)}M"
    return res if res != "P" else ""


def map_results(
    mapping, rec, explicit_hpos=None, explicit_omims=None, full_text_context=""
):
    # Store phenotypes as dictionary to maintain pairing
    # key: hpo_id, value: {label, severity_id, severity_label, excluded}
    phenotypes = {}

    # Store OMIMs and MONDOs (sets are fine for these as they don't have modifiers in this context)
    omim_g, omim_p_set, omim_l, omim_e = set(), set(), set(), set()
    mondo_p, mondo_l, mondo_e = set(), set(), set()

    # Process explicit OMIMs first
    for oid in explicit_omims or []:
        otype = utils.classify_omim(oid)
        if otype == "gene":
            omim_g.add(oid)
        else:
            omim_p_set.add(oid)
        if name := utils.disease_map.get(oid):
            omim_l.add(name)

    # Process explicit HPOs (default: present, no severity)
    for hid in explicit_hpos or []:
        phenotypes[hid] = {
            "label": next((t["name"] for t in utils.hpo_terms if t["id"] == hid), hid),
            "excluded": False,
            "severity_id": "",
            "severity_label": "",
        }

    if mapping and isinstance(mapping, dict):
        for term, data in mapping.items():
            if not isinstance(data, dict):
                continue

            # HPO handling - LLM returns LABELS ONLY, we look up IDs
            hpo_labels_data = data.get(
                "hpo_labels", data.get("hpo")
            )  # Support both formats
            is_excluded = data.get("hpo_excluded", False)

            # Normalize exclusion flag
            if isinstance(is_excluded, str):
                is_excluded = is_excluded.lower() in ["true", "yes"]

            # Normalize to list
            if isinstance(hpo_labels_data, str):
                hpo_labels_data = [hpo_labels_data] if hpo_labels_data != "None" else []
            elif isinstance(hpo_labels_data, dict):
                # Old format support (will be removed after full migration)
                if hpo_labels_data.get("id"):
                    hpo_labels_data = [hpo_labels_data]
                else:
                    hpo_labels_data = []
            elif not isinstance(hpo_labels_data, list):
                hpo_labels_data = []

            for h in hpo_labels_data:
                # If old format (dict with id), extract label
                if isinstance(h, dict):
                    label = h.get("label")
                    # Check if LLM also provided ID (old behavior) - we'll verify it
                    llm_id = (
                        h.get("id") if str(h.get("id", "")).startswith("HP:") else None
                    )
                    is_excluded = h.get("excluded", is_excluded)
                    if isinstance(is_excluded, str):
                        is_excluded = is_excluded.lower() in ["true", "yes"]
                else:
                    # New format: just label string
                    label = h
                    llm_id = None

                if not label or label in ["None", "null", ""]:
                    continue

                # Look up correct ID from ontology by label
                hid = utils.lookup_hpo_id_by_label(label)

                if not hid:
                    print(
                        f"Warning: Could not find HPO ID for label '{label}' - skipping"
                    )
                    continue

                # Verify if LLM also provided an ID (old behavior)
                if llm_id and llm_id != hid:
                    print(
                        f"Warning: LLM ID mismatch for '{label}': LLM={llm_id}, Actual={hid}"
                    )
                    # Use ontology ID, not LLM ID

                # Handle Severity - STRICT VALIDATION
                sev_id = ""
                sev_label = ""

                # Extract severity from dict format only (not string)
                if isinstance(h, dict):
                    sev = h.get("severity")
                    if sev and str(sev).startswith("HP:"):
                        if utils.validate_severity(sev):
                            sev_id = sev
                            sev_label = next(
                                (t["name"] for t in utils.hpo_terms if t["id"] == sev),
                                sev,
                            )
                        else:
                            # Invalid severity ID - ignore
                            pass
                    elif sev and str(sev) != "None":
                        # Sometimes LLM gives text like "Severe" - look up ID
                        sev_id_lookup = utils.lookup_hpo_id_by_label(str(sev))
                        if sev_id_lookup and utils.validate_severity(sev_id_lookup):
                            sev_id = sev_id_lookup
                            sev_label = str(sev)

                # Update/Overwrite existing entry
                phenotypes[hid] = {
                    "label": label,
                    "excluded": is_excluded,
                    "severity_id": sev_id,
                    "severity_label": sev_label,
                }

            # OMIM handling - STRICT VALIDATION against text if not explicit
            o = data.get("omim")
            if isinstance(o, dict) and o.get("id") and str(o["id"]).startswith("OMIM:"):
                oid = o["id"]
                # Only accept LLM-derived OMIM if it appears in the text OR is already explicit
                if (
                    oid in omim_g or oid in omim_p_set
                ) or utils.validate_omim_against_text(oid, full_text_context):
                    label = o.get("label")
                    if not label or label == "None":
                        label = utils.disease_map.get(oid, oid)

                    if o.get("excluded"):
                        omim_e.add(oid)
                    else:
                        otype = utils.classify_omim(oid)
                        (omim_g if otype == "gene" else omim_p_set).add(oid)
                        omim_l.add(label)

            # MONDO handling
            m = data.get("mondo")
            if (
                isinstance(m, dict)
                and m.get("id")
                and str(m["id"]).startswith("MONDO:")
            ):
                mid = m["id"]
                label = m.get("label")
                if not label or label == "None":
                    label = utils.mondo.get(mid, mid)

                if m.get("excluded"):
                    mondo_e.add(mid)
                else:
                    mondo_p.add(mid)
                    mondo_l.add(label)

    # Separate and Sort Phenotypes
    present_list = sorted([k for k, v in phenotypes.items() if not v["excluded"]])
    excluded_list = sorted([k for k, v in phenotypes.items() if v["excluded"]])

    # Construct Paired Lists
    p_ids, p_labels, p_sev_ids, p_sev_labels = [], [], [], []
    for hid in present_list:
        p_ids.append(hid)
        p_labels.append(phenotypes[hid]["label"])
        p_sev_ids.append(phenotypes[hid]["severity_id"])
        p_sev_labels.append(phenotypes[hid]["severity_label"])

    e_ids, e_labels = [], []
    for hid in excluded_list:
        e_ids.append(hid)
        e_labels.append(phenotypes[hid]["label"])

    rec.update(
        {
            "hpo_ids": "|".join(p_ids),
            "hpo_labels": "|".join(p_labels),
            "hpo_excluded_ids": "|".join(e_ids),
            "hpo_excluded_labels": "|".join(e_labels),
            "hpo_severity_ids": "|".join(p_sev_ids),
            "hpo_severity_labels": "|".join(p_sev_labels),
            "disease_omim_gene_ids": "|".join(sorted(list(omim_g))),
            "disease_omim_phenotype_ids": "|".join(sorted(list(omim_p_set))),
            "disease_omim_labels": "|".join(sorted(list(set(filter(None, omim_l))))),
            "disease_omim_excluded_ids": "|".join(sorted(list(omim_e))),
            "disease_mondo_ids": "|".join(sorted(list(mondo_p))),
            "disease_mondo_labels": "|".join(sorted(list(set(filter(None, mondo_l))))),
            "disease_mondo_excluded_ids": "|".join(sorted(list(mondo_e))),
        }
    )


def get_publication_id(row, filename):
    if filename == "ahmed-pmid28454995.tsv":
        if (
            pd.notna(row.get("References PMID"))
            and str(row.get("References PMID")).strip()
        ):
            return f"PMID:{str(row.get('References PMID')).strip()}"
        return "PMID:28454995"
    elif filename == "marwa-variants.tsv":
        ref = str(row.get("reference", ""))
        m = re.search(r"pubmed/(\d+)", ref) or re.search(
            r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", ref
        )
        if m:
            return f"PMID:{m.group(1)}"
    elif filename == "PMC6562004.tsv":
        return "PMID:31182893"
    elif filename == "PMC7082194.tsv":
        return "PMID:32256298"
    elif filename == "ddd-diagnoses.tsv":
        pmids = str(row.get("pubmed_ids", ""))
        if pmids and pmids != "nan":
            return "|".join(
                [f"PMID:{p.strip()}" for p in pmids.split("|") if p.strip()]
            )
    return ""


def process_file(filename, limit=None, indices=None):
    path = f"data/phenotypes/{filename}"
    if not os.path.exists(path):
        return
    print(f"\nProcessing {filename}...")

    # Pre-read PMC6562004 mapping
    pmc6562004_map = {}
    if filename == "PMC6562004.tsv":
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                header_line = f.readline().strip()
            if header_line.startswith("#Primary indication:"):
                parts = header_line.replace("#Primary indication:", "").split("/")
                for p in parts:
                    m = re.search(r"(.+?)\s*\((\d+)\)", p)
                    if m:
                        pmc6562004_map[m.group(2)] = m.group(1).strip()
        except Exception as e:
            print(f"Warning: Could not parse PMC6562004 header: {e}")

    if filename == "ahmed-pmid28454995.tsv":
        df = pd.read_csv(path, sep="\t")
        df.columns = [
            "Case",
            "Diagnosis",
            "OMIM",
            "Phenotype",
            "Variant 1 Gene",
            "Transcript ID",
            "g.DNA",
            "c.DNA",
            "Amino Acid",
            "Variant 1 Allele State",
            "dbSNP",
            "ESP",
            "ExAC",
            "ClinVar",
            "AF Controls",
            "AF Internal",
            "PolyPhen-2",
            "SIFT",
            "Align GVGD",
            "PMID",
            "SEP",
            "Variant 2 Gene",
            "Transcript ID 2",
            "g.DNA 2",
            "c.DNA 2",
            "Amino Acid 2",
            "Variant 2 Allele State",
            "dbSNP 2",
            "ESP 2",
            "MAF ExAC 2",
            "AF Controls 2",
            "AF Internal 2",
            "ClinVar 2",
            "PolyPhen-2 2",
            "SIFT 2",
            "Align GVGD 2",
            "References 2",
        ]
    else:
        df = pd.read_csv(path, sep="\t", comment="#").rename(
            columns=lambda x: str(x).strip()
        )

    if indices:
        df = df.iloc[indices]
    elif limit:
        df = df.head(limit)
    if filename == "PMC7082194.tsv":
        df_gen = df.iloc[:, :10].copy()
        df_phen = df.iloc[:, 10:14].copy()
        df_phen.columns = ["ID_p", "Sex_p", "Age_p", "Phenotype_p"]
        df_gen["ID_c"] = (
            pd.to_numeric(df_gen["ID"], errors="coerce")
            .fillna(-1)
            .astype(int)
            .astype(str)
        )
        df_phen["ID_c"] = (
            pd.to_numeric(df_phen["ID_p"], errors="coerce")
            .fillna(-1)
            .astype(int)
            .astype(str)
        )
        df = df_phen.merge(df_gen, on="ID_c", how="left")

    # Incremental processing: load existing output and skip already processed rows
    output_path = os.path.join(OUTPUT_DIR, filename.replace(".tsv", "_normalized.tsv"))
    processed_ids = set()

    if os.path.exists(output_path):
        try:
            existing_df = pd.read_csv(output_path, sep="\t")
            processed_ids = set(existing_df["id"].astype(str))
            print(f"  Resuming: {len(processed_ids)} rows already processed")
        except Exception as e:
            print(f"  Warning: Could not read existing file: {e}")

    # Open file in append mode if it exists, otherwise write header
    file_mode = "a" if os.path.exists(output_path) else "w"
    write_header = not os.path.exists(output_path)

    records = []
    rows_processed = 0

    for idx, row in tqdm(df.iterrows(), total=len(df)):
        rec = {k: "" for k in COLUMNS}
        rec["source_file"] = filename

        raw_id = row.get("Case", row.get("ID", row.get("ID_p", idx)))
        try:
            # Try to convert to float then int to remove .0
            f_id = float(raw_id)
            if f_id.is_integer():
                rec["id"] = str(int(f_id))
            else:
                rec["id"] = str(raw_id)
        except:
            rec["id"] = str(raw_id)

        # Skip if already processed (incremental resume)
        if rec["id"] in processed_ids:
            continue

        # Metadata
        rec["sex_id"], rec["sex_label"] = utils.normalize_sex(
            row.get("Gender", row.get("Sex", row.get("Sex_p", "")))
        )
        if filename == "ahmed-pmid28454995.tsv":
            if any(
                "hemi" in str(row.get(f"Variant {i} Allele State", "")).lower()
                for i in [1, 2]
            ):
                rec["sex_id"], rec["sex_label"] = utils.normalize_sex("Male")

        rec["age"] = parse_iso_age(
            row.get("Age", row.get("Age_p", row.get("Age (years, months)", "")))
        )
        if filename != "ddd-diagnoses.tsv":
            rec["ancestry_id"], rec["ancestry_label"] = utils.normalize_ancestry()
            rec["country_id"], rec["country_label"] = utils.normalize_country()
        rec["consanguinity_id"], rec["consanguinity_label"] = (
            utils.normalize_consanguinity(
                row.get("Consanguinity", row.get("Consanguinity ", ""))
            )
        )
        rec["family_history_id"], rec["family_history_label"] = (
            utils.normalize_family_history(
                row.get("Family hx", row.get("Family History", ""))
            )
        )
        rec["genotype_id"], rec["genotype_label"] = utils.normalize_genotype(
            row.get(
                "Zygosity",
                row.get(
                    "Zyogsity",
                    row.get("allelic_mode", row.get("Variant 1 Allele State", "")),
                ),
            )
        )

        # Variants & Genes
        genes, v_hgvs, v_valid = {}, [], []
        v_texts = []
        if filename == "fawzan-variants.tsv":
            v_texts = [str(row.get("Variant(s)", "")), str(row.get("Variant(s).1", ""))]
        elif filename == "ahmed-variants.tsv":
            for s in ["", " 2", " 3"]:
                if pd.notna(g := row.get(f"Gene{s}".strip())):
                    v_texts.append(f"{g} {row.get(f'Variant{s}'.strip(), '')}")
        elif filename == "marwa-variants.tsv":
            v_texts = [str(row.get("variants", ""))]
        elif filename == "PMC6562004.tsv":
            v_texts = [f"{row.get('Gene(s)', '')} {row.get('Variant(s)', '')}"]
        elif filename == "PMC7082194.tsv":
            v_texts = [
                f"{row.get('Gene(s)', '')} {row.get('HGVS Format', '')}",
                f"{row.get('Gene(s)', '')} {row.get('Nucleotide Change', '')}",
            ]
        elif filename == "ahmed-pmid28454995.tsv":
            v_texts = [
                f"{row.get('Variant 1 Gene', '')} {row.get('c.DNA', '')}",
                f"{row.get('Variant 2 Gene', '')} {row.get('c.DNA 2', '')}",
            ]
        elif filename == "ddd-diagnoses.tsv":
            v_texts = [str(row.get("gene", ""))]

        for vt in v_texts:
            vt = str(vt).strip()
            if not vt or vt.lower() in ["nan", "0", "none"]:
                continue
            m_g = re.match(r"^([A-Z0-9-]+)", vt.upper())
            g_sym = m_g.group(1) if m_g else vt.split(":")[0].strip().upper()
            # Skip if extracted gene symbol is 'NAN' (pandas nan artifact)
            # NAN is a valid gene synonym for SCN11A, but we don't want to map it from nan values
            if g_sym == "NAN":
                continue
            if gn := utils.map_gene(g_sym):
                genes[gn["entrez"]] = gn
                v_hgvs.append(utils.extract_variant_hgvs(vt, gn["symbol"]))
                v_valid.append(
                    "True" if utils.validate_variant(gn["symbol"]) else "False"
                )
            else:
                if v_res := utils.extract_variant_hgvs(vt):
                    v_hgvs.append(v_res)
                    for ev in v_res.split("|"):
                        if egn := utils.map_gene(ev.split(":")[0]):
                            genes[egn["entrez"]] = egn
                            v_valid.append(
                                "True"
                                if utils.validate_variant(egn["symbol"])
                                else "False"
                            )

        # Get Result field for solved status (fawzan-variants.tsv has "Positive"/"Negative"/"Ambiguous")
        result_field = str(row.get("Result", row.get("category", "")))

        # Determine "solved" status (lowercase true/false for boolean consistency)
        solved = ""
        result_lower = result_field.lower().strip()
        if "positive" in result_lower:
            solved = "true"
        elif "negative" in result_lower:
            solved = "false"
        elif any(x in result_lower for x in ["ambiguous", "vous", "vus", "uncertain"]):
            solved = "false"

        # Get variant classification
        # For fawzan-variants, use Result field; for others, use ACMG Classification
        if filename == "fawzan-variants.tsv":
            v_class, v_evidence = utils.normalize_variant_classification(result_field)
        else:
            v_class, v_evidence = utils.normalize_variant_classification(
                str(row.get("ACMG Classification", row.get("category", "")))
            )

        # Special handling for fawzan-variants: infer classification from solved status if needed
        if not v_class and solved == "true" and (v_hgvs and any(v_hgvs)):
            v_class = "Pathogenic"

        # Default solved=true for all files if variant is present (unless Result field says otherwise)
        if not solved and (v_hgvs and any(v_hgvs)):
            solved = "true"

        rec.update(
            {
                "publication_id": get_publication_id(row, filename),
                "gene_hgnc_ids": "|".join([x["hgnc"] for x in genes.values()]),
                "gene_symbols": "|".join([x["symbol"] for x in genes.values()]),
                "variant_hgvs": "|".join(v_hgvs),
                "variant_validated": "|".join(v_valid),
                "variant_classification": v_class,
                "variant_evidence": v_evidence,
                "solved": solved,
            }
        )

        # Phenotypes & Diseases
        p_val = str(
            row.get(
                "Phenotype",
                row.get(
                    "HPOs",
                    row.get(
                        "Additional clinical phenotype",
                        row.get(
                            "phenotypes", row.get("Phenotype_p", row.get("hpo_ids", ""))
                        ),
                    ),
                ),
            )
        )
        d_val = str(
            row.get("Diagnosis", row.get("Clinical Snydrome", row.get("syndrome", "")))
        )

        if filename == "PMC6562004.tsv":
            ind_code = str(row.get("Primary indication", ""))
            if ind_code in pmc6562004_map:
                d_val = f"{d_val} {pmc6562004_map[ind_code]}".strip()

        # Robust OMIM extraction
        omim_cols = [c for c in df.columns if "mim" in c.lower()]
        o_vals = [
            str(row[c])
            for c in omim_cols
            if pd.notna(row[c]) and str(row[c]).strip().lower() not in ["", "nan"]
        ]
        o_val = o_vals[0] if o_vals else ""

        full_text = f"Phenotypes: {p_val}. Diagnosis: {d_val}. OMIM: {o_val}."
        explicit_hpos = re.findall(r"HP:\d{7}", full_text)

        explicit_omims = []
        # 1. From OMIM column (o_val)
        potential_omims = re.findall(r"\d{6}", o_val)
        # 2. From Diagnosis/Phenotype text (strict: 6 digits NOT preceded by HP:)
        potential_omims.extend(
            re.findall(r"(?<!HP:)(?<!HP)\b\d{6}\b", f"{d_val} {p_val}")
        )

        for oid in set(potential_omims):
            if oid in utils.mim_info or oid in utils.disease_map:
                explicit_omims.append(f"OMIM:{oid}")

        # Split ONLY phenotype text into terms for LLM
        # Diagnosis (d_val) is handled separately via OMIM/MONDO extraction
        # DO NOT send diagnosis to LLM for HPO normalization - it causes hallucinations
        raw_terms = []
        if p_val and p_val != "nan":
            initial_split = [
                t.strip()
                for t in re.split(r"[,|;]\s*", p_val)
                if t.strip() and t.lower() != "nan"
            ]

            # Extract phenotypes from complex descriptive phrases
            for term in initial_split:
                # If term starts with non-phenotype context words, extract phenotypes from it
                if re.match(
                    r"^(picu|icu|admission|admitted|caught|diagnosed|presented|suspected)\b",
                    term.lower(),
                ):
                    # Extract phenotypes mentioned after "due to", "with", "and", etc.
                    # e.g., "PICU addition due to ventricular arrhythmia and post cardiac arrest"
                    sub_terms = re.split(
                        r"\s+(?:and|due to|with|following|after)\s+",
                        term,
                        flags=re.IGNORECASE,
                    )
                    for st in sub_terms:
                        st_clean = re.sub(
                            r"^(due to|with|and|following|after|post)\s+",
                            "",
                            st,
                            flags=re.IGNORECASE,
                        ).strip()
                        # Only keep if it doesn't start with context words
                        if st_clean and not re.match(
                            r"^(picu|icu|admission|addition|caught|admitted)\b",
                            st_clean.lower(),
                        ):
                            raw_terms.append(st_clean)
                else:
                    # Keep the term as-is
                    raw_terms.append(term)

        terms = []
        for t in raw_terms:
            t_disambiguated = utils.disambiguate_term(
                t, full_text, rec.get("gene_symbols")
            )
            terms.append(t_disambiguated)

        if terms:
            map_results(
                utils.llm_normalize(terms, context=rec["gene_symbols"]),
                rec,
                explicit_hpos,
                explicit_omims,
                full_text_context=full_text,
            )

        records.append(rec)
        rows_processed += 1

        # Incremental write: append row immediately to file
        row_df = pd.DataFrame([rec])
        row_df.to_csv(
            output_path,
            sep="\t",
            mode=file_mode,
            header=write_header,
            index=False,
        )
        # After first write, switch to append mode without header
        if write_header:
            write_header = False
            file_mode = "a"

    print(
        f"  ✅ Processed {rows_processed} new rows (total in file: {rows_processed + len(processed_ids)})"
    )


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(
        description="Process phenotype files with incremental progress."
    )
    parser.add_argument(
        "limit",
        type=int,
        nargs="?",
        default=None,
        help="Limit number of records per file",
    )
    parser.add_argument(
        "--files",
        nargs="+",
        help="Specific files to process (default: all)",
        default=None,
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete existing normalized files and reprocess from scratch",
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Process files in parallel (up to 7 files simultaneously)",
    )
    parser.add_argument(
        "--graph-rag", action="store_true", help="Enable GraphRAG normalization"
    )

    args = parser.parse_args()

    if args.graph_rag:
        utils.enable_graph_rag()

    all_files = [
        "ahmed-variants.tsv",
        "ahmed-pmid28454995.tsv",
        "fawzan-variants.tsv",
        "marwa-variants.tsv",
        "PMC6562004.tsv",
        "PMC7082194.tsv",
        "ddd-diagnoses.tsv",
    ]

    files_to_process = args.files if args.files else all_files

    # Filter to only valid files
    files_to_process = [f for f in files_to_process if f in all_files]

    if not files_to_process:
        print("No valid files specified. Available files:")
        for f in all_files:
            print(f"  - {f}")
        sys.exit(1)

    # Handle overwrite
    if args.overwrite:
        for f in files_to_process:
            output_file = os.path.join(OUTPUT_DIR, f.replace(".tsv", "_normalized.tsv"))
            if os.path.exists(output_file):
                print(f"🗑️  Deleting existing: {output_file}")
                os.remove(output_file)

    # Check which files need processing (incremental)
    files_needing_work = []
    for f in files_to_process:
        output_file = os.path.join(OUTPUT_DIR, f.replace(".tsv", "_normalized.tsv"))
        if not os.path.exists(output_file):
            files_needing_work.append(f)
        else:
            print(f"✓ Skipping {f} (already processed, use --overwrite to reprocess)")

    if not files_needing_work:
        print("\n✅ All files already processed! Use --overwrite to reprocess.")
        sys.exit(0)

    print(f"\n📋 Processing {len(files_needing_work)} file(s)...")

    # Process files
    if args.parallel and len(files_needing_work) > 1:
        print("🚀 Parallel mode: Processing files concurrently")
        with ThreadPoolExecutor(
            max_workers=min(len(files_needing_work), 7)
        ) as executor:
            futures = {
                executor.submit(process_file, f, args.limit): f
                for f in files_needing_work
            }
            for future in as_completed(futures):
                filename = futures[future]
                try:
                    future.result()
                    print(f"✅ Completed: {filename}")
                except Exception as e:
                    print(f"❌ Error processing {filename}: {e}")
    else:
        for f in files_needing_work:
            process_file(f, limit=args.limit)
            print(f"✅ Completed: {f}")
