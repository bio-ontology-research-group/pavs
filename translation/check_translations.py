import json
import re
import collections
from typing import List, Dict, Any

def contains_latin(text: str) -> bool:
    """Check if a string contains Latin (English) characters."""
    if not text:
        return False
    return bool(re.search(r'[a-zA-Z]', text))

def check_translations(json_path: str):
    with open(json_path, 'r', encoding='utf-8') as f:
        translations = json.load(f)

    print(f"Analyzing {len(translations)} terms...\n")

    english_in_arabic = []
    suspicious_lengths = []
    
    # For consistency checking
    # Common prefixes/words to track
    keywords = ["Abnormal", "Delayed", "Increased", "Decreased", "Congenital", "Abnormality"]
    consistency_map = collections.defaultdict(lambda: collections.defaultdict(int))

    for trans in translations:
        hp_id = trans['id']
        en_name = trans.get('english_technical_name', '')
        ar_name = trans.get('arabic_technical_name', '')
        ar_def = trans.get('arabic_definition', '')

        # 1. English Detection
        if contains_latin(ar_name) or contains_latin(ar_def):
            # Only flag if it's more than just a stray character or if it looks like a word
            if re.search(r'[a-zA-Z]{3,}', ar_name) or re.search(r'[a-zA-Z]{3,}', ar_def):
                english_in_arabic.append({
                    'id': hp_id,
                    'en': en_name,
                    'ar_name': ar_name,
                    'ar_def': ar_def
                })

        # 2. Length Check (Ratio)
        if en_name and ar_name:
            ratio = len(ar_name) / len(en_name)
            if ratio < 0.4 or ratio > 2.5:
                suspicious_lengths.append({
                    'id': hp_id,
                    'en': en_name,
                    'ar': ar_name,
                    'ratio': round(ratio, 2)
                })

        # 3. Consistency Prep
        for kw in keywords:
            if en_name.lower().startswith(kw.lower()):
                # Get the first word of the Arabic translation
                parts = ar_name.split()
                ar_first_word = parts[0] if parts else ""
                if ar_first_word:
                    consistency_map[kw][ar_first_word] += 1

    # --- REPORTING ---

    print("--- 1. ENGLISH DETECTED IN ARABIC FIELDS ---")
    if not english_in_arabic:
        print("None found.")
    else:
        for item in english_in_arabic[:15]:
            print(f"[{item['id']}] {item['en']} -> {item['ar_name']}")
        if len(english_in_arabic) > 15:
            print(f"... and {len(english_in_arabic) - 15} more.")
    print("\n")

    print("--- 2. SUSPICIOUS LENGTH RATIOS (AR/EN) ---")
    if not suspicious_lengths:
        print("None found.")
    else:
        sorted_lengths = sorted(suspicious_lengths, key=lambda x: abs(1-x['ratio']), reverse=True)
        for item in sorted_lengths[:15]:
            print(f"[{item['id']}] Ratio {item['ratio']}: '{item['en']}' -> '{item['ar']}'")
        if len(suspicious_lengths) > 15:
            print(f"... and {len(suspicious_lengths) - 15} more.")
    print("\n")

    print("--- 3. CONSISTENCY CHECK (Common Prefixes) ---")
    for kw, variations in consistency_map.items():
        print(f"English: '{kw}'")
        sorted_vars = sorted(variations.items(), key=lambda x: x[1], reverse=True)
        for ar_word, count in sorted_vars[:5]:
            print(f"  - {ar_word}: {count} times")
        if len(sorted_vars) > 5:
            print(f"  - ... and {len(sorted_vars) - 5} other variations.")
        print()

if __name__ == "__main__":
    check_translations("translation/hpo_arabic_translations.json")
