"""
Prepare HPO Arabic babelon TSV for submission to hpo-translations repo.

- All terms from hp-ar.babelon.tsv are included (CANDIDATE status).
- Terms in reviewed-translations.tsv with "Reviewed" or "Reviewed can be ..."
  are promoted to OFFICIAL status and use the reviewed translation value.
- Output uses the 7-column upstream format matching the hpo-translations repo.
"""

import csv
import sys

REVIEWED_STATUSES = {"Reviewed", "Reviewed can be improved", "Reviewed can be modified"}

BABELON_IN  = "translation/hp-ar.babelon.tsv"
REVIEWED_IN = "translation/reviewed-translations.tsv"
OUTPUT      = "translation/hp-ar-submission.babelon.tsv"

# Upstream column order (7 columns)
FIELDNAMES = [
    "source_language",
    "translation_language",
    "subject_id",
    "predicate_id",
    "source_value",
    "translation_value",
    "translation_status",
]

def load_reviewed(path):
    """Return dict: hp_id -> {label, label_status, layman, layman_status, defn, defn_status}"""
    reviewed = {}
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            hp_id = row["id"].strip()
            if not hp_id:
                continue
            reviewed[hp_id] = {
                "label":        row.get("arabic_technical_name", "").strip(),
                "label_status": row.get("Reviewed", "").strip(),
                "defn":         row.get("arabic_definition", "").strip(),
                "defn_status":  row.get("Reviewed def", "").strip(),
            }
    return reviewed

def main():
    reviewed = load_reviewed(REVIEWED_IN)

    # Track stats
    total = 0
    promoted_label = 0
    promoted_defn  = 0

    with open(BABELON_IN, encoding="utf-8") as fin, \
         open(OUTPUT, "w", encoding="utf-8", newline="") as fout:

        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()

        for row in reader:
            hp_id     = row["subject_id"].strip()
            predicate = row["predicate_id"].strip()
            src_lang  = row.get("source_language", "en").strip()
            src_val   = row.get("source_value", "").strip()
            tr_val    = row.get("translation_value", "").strip()
            status    = "CANDIDATE"

            rev = reviewed.get(hp_id)
            if rev:
                if predicate == "rdfs:label" and rev["label_status"] in REVIEWED_STATUSES:
                    reviewed_val = rev["label"]
                    if reviewed_val:
                        tr_val = reviewed_val
                        status = "OFFICIAL"
                        promoted_label += 1
                elif predicate == "IAO:0000115" and rev["defn_status"] in REVIEWED_STATUSES:
                    reviewed_val = rev["defn"]
                    if reviewed_val:
                        tr_val = reviewed_val
                        status = "OFFICIAL"
                        promoted_defn += 1

            writer.writerow({
                "source_language":    src_lang,
                "translation_language": "ar",
                "subject_id":         hp_id,
                "predicate_id":       predicate,
                "source_value":       src_val,
                "translation_value":  tr_val,
                "translation_status": status,
            })
            total += 1

    print(f"Written {total} rows to {OUTPUT}")
    print(f"  OFFICIAL (label):      {promoted_label}")
    print(f"  OFFICIAL (definition): {promoted_defn}")
    print(f"  CANDIDATE:             {total - promoted_label - promoted_defn}")

if __name__ == "__main__":
    main()
