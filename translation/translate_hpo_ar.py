import argparse
import asyncio
import json
import logging
import os
import random
from typing import List, Dict, Any, Optional

import aiohttp
import pronto

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

OPENROUTER_URL = "https://openrouter.ai/api/v1/chat/completions"

SYSTEM_PROMPT = """You are a bilingual medical geneticist and expert translator. Translate the following JSON array of Human Phenotype Ontology (HPO) terms to Arabic.

STRICT RULES:
1. For 'arabic_technical_name': Provide the standard medical Arabic translation. If (and only if) no standard Arabic term exists in medical literature, keep the English name.
2. For 'arabic_definition': ALWAYS translate the 'english_definition' into clear, descriptive medical Arabic. Never return English for the definition if an English definition was provided.
3. For 'arabic_layperson_synonym': Provide a simplified Arabic version that a non-medical person would understand. If uncertain, you may use a simplified version of the technical name.
4. Output strictly in the specified JSON format. Ensure the 'translations' array has the same length as the input.

ARABIC SENTENCE STRUCTURE FOR HPO TERMS:
In Arabic HPO translations, always state the problem first, then the anatomical location. This is the reverse of the English structure.
Rule: Translate as: [Arabic term for abnormality] + في + [anatomical location]. Do not mirror English word order.
Examples:
- Fingernail hypoplasia → نقص تنسج في أظافر اليد
- Renal tubular dysfunction → خلل في الأنابيب الكلوية
- Cardiac rhythm abnormality → اضطراب في نظم القلب

Exception for adjectives:
- Metaphyseal dysplasia → خلل تنسجي كردوسي (Abnormality + Adjective)

GLOSSARY (Use these EXACT terms):
- Atresia: رتق
- Malacia / -malacia: تلين
- Metaphysis: كردوس
- Metaphyseal (adjective): كردوسي
- interface: الواجهة / السطح البيني
- Serrated: مسنّن
- -cele: قيلة
- -pathy: اعتلال
- Dystrophy: حثل
- Atrophy: ضمور
- Mesomelic: متوسط الأطراف
- Intraepithelial: داخل الظهارة
- Adenoid: اللحمية
- Anteversion: انكفاء أمامي
- Retroversion: انكفاء خلفي   
- Polyps: سلائل
- Capital femoral: رأس عظم الفخذ
- Epiphysis: مشاشة
- Notch: ثلمة
- Cleft / Fissure: شق
- Hyperplasia: فرط تنسج
- Hypertrophy: تضخم
- Neoplasia: تنشّؤ
- Dysplasia: خلل تنسجي

SPECIAL RULES FOR TERMS CONTAINING "Abnormal" OR "Abnormality":
Apply the following decision rules IN ORDER, stopping at the first matching rule:

Rule 1 — Physical Structure:
  If the term describes the physical shape, size, or morphology of an anatomical structure (organ, bone, tissue, cell):
  → Use شذوذ (Shudhūdh) as the translation for "Abnormal/Abnormality".
  → If the structural change is at the histological or microscopic level, use شذوذ نسيجي.
  → If Rule 1 applies, do NOT proceed to further rules.

Rule 2 — Measurable Functional Deficit:
  If the abnormality is a specific, measurable failure quantified as an output, level, activity, or pressure:
  → Use خلل (Khalal).
  → Append وظيفي only when the context explicitly denotes organ or physiological function.
  → Omit وظيفي for enzyme activity and biochemical levels.
  → If Rule 2 applies, do NOT proceed to further rules.

Rule 3 — Complex System, Pattern, or Behavior:
  If the term describes a complex pattern, rhythm, cycle, behavioral pattern, or a disturbance involving multiple interacting factors:
  → Use اضطراب (Iḍitārāb).
  → If Rule 3 applies, do NOT proceed to further rules.

Rule 4 — Default Descriptive Quality (use sparingly):
  For basic sensory or simple descriptive traits (color, odor, appearance, basic lab finding):
  → Use غير طبيعي (Ghayr Ṭabīʼi).
  → Re-examine Rules 1–3 before using this fallback.

SPECIAL RULES FOR MORPHOGENESIS TERMS:
These rules take STRICT PRIORITY over all other rules.

Hypoplasia (underdevelopment / reduced size of a normally-formed structure):
  → arabic_technical_name MUST begin with: نقص تنسج
  → Pattern: نقص تنسج + في + [anatomical structure]
  → NEVER use نقص التطور, نقص التخلق, ضمور, صغر, or any other variant.

Aplasia (complete congenital absence of tissue/structure formation):
  → arabic_technical_name MUST begin with: عدم التنسج
  → Pattern: عدم التنسج + في + [anatomical structure]
  → NEVER use عدم التكون, غياب, انعدام, غَيْبَة, فقد تطوري, or any other variant.

Agenesis (embryological failure of an organ to develop):
  → arabic_technical_name MUST begin with: عدم التكوّن
  → Pattern: عدم التكوّن + في + [anatomical structure]
  → NEVER use غياب, عدم تشكّل, عدم نمو, or any other variant.

Combined Aplasia/Hypoplasia terms:
  → Use: عدم التنسج أو نقص التنسج + في + [anatomical structure]
"""

RESPONSE_SCHEMA = {
    "type": "json_schema",
    "json_schema": {
        "name": "hpo_translation_batch",
        "strict": True,
        "schema": {
            "type": "object",
            "properties": {
                "translations": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "id": {"type": "string"},
                            "arabic_technical_name": {"type": "string"},
                            "arabic_layperson_synonym": {"type": "string"},
                            "arabic_definition": {"type": "string"}
                        },
                        "required": ["id", "arabic_technical_name", "arabic_layperson_synonym", "arabic_definition"],
                        "additionalProperties": False
                    }
                }
            },
            "required": ["translations"],
            "additionalProperties": False
        }
    }
}

class UsageTracker:
    def __init__(self, price_per_1m_prompt=0.15, price_per_1m_completion=0.60):
        self.total_prompt_tokens = 0
        self.total_completion_tokens = 0
        self.price_per_1m_prompt = price_per_1m_prompt
        self.price_per_1m_completion = price_per_1m_completion

    def add_usage(self, usage: Dict[str, Any]):
        self.total_prompt_tokens += usage.get("prompt_tokens", 0)
        self.total_completion_tokens += usage.get("completion_tokens", 0)

    @property
    def total_cost(self):
        prompt_cost = (self.total_prompt_tokens / 1_000_000) * self.price_per_1m_prompt
        comp_cost = (self.total_completion_tokens / 1_000_000) * self.price_per_1m_completion
        return prompt_cost + comp_cost

    def __str__(self):
        return (f"Tokens: [Prompt: {self.total_prompt_tokens}, Completion: {self.total_completion_tokens}] "
                f"Est. Cost: ${self.total_cost:.4f}")

usage_tracker = UsageTracker()

async def translate_batch(
    session: aiohttp.ClientSession,
    batch_payload: List[Dict[str, Any]],
    api_key: str,
    model: str,
    semaphore: asyncio.Semaphore,
    max_retries: int = 5
) -> Optional[List[Dict[str, Any]]]:
    """Translate a batch of HPO terms using OpenRouter."""
    payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": json.dumps(batch_payload)}
        ],
        "response_format": RESPONSE_SCHEMA
    }

    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "X-Title": "PAVS HPO Batch Translation Tool"
    }

    async with semaphore:
        for attempt in range(max_retries):
            try:
                async with session.post(OPENROUTER_URL, headers=headers, json=payload, timeout=300) as resp:
                    if resp.status == 200:
                        data = await resp.json()
                        usage_tracker.add_usage(data.get("usage", {}))
                        content_str = data['choices'][0]['message']['content']
                        content = json.loads(content_str)
                        return content.get('translations', [])
                    elif resp.status == 429:
                        wait_time = (2 ** attempt) + random.uniform(0.5, 1.5)
                        logger.warning(f"Rate limited (429). Waiting {wait_time:.2f}s...")
                        await asyncio.sleep(wait_time)
                    else:
                        resp_text = await resp.text()
                        logger.error(f"API error {resp.status}: {resp_text}")
                        if resp.status >= 500:
                            wait_time = (2 ** attempt) + random.uniform(0.5, 1.5)
                            await asyncio.sleep(wait_time)
                        else:
                            break
            except Exception as e:
                logger.error(f"Request exception: {str(e)}")
                wait_time = (2 ** (attempt + 1)) + random.uniform(0.5, 1.5)
                await asyncio.sleep(wait_time)
        
        return None

def chunk_list(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

async def process_ontology(args):
    api_key = args.api_key or os.environ.get("OPENROUTER_API_KEY")
    if not api_key:
        logger.error("OpenRouter API key not provided.")
        return
    
    # Update tracker prices if provided
    usage_tracker.price_per_1m_prompt = args.price_prompt
    usage_tracker.price_per_1m_completion = args.price_completion

    existing_translations = {}
    if os.path.exists(args.output):
        try:
            with open(args.output, 'r', encoding='utf-8') as f:
                content = json.load(f)
                if isinstance(content, list):
                    existing_translations = {t['id']: t for t in content}
                    logger.info(f"Resuming: Loaded {len(existing_translations)} existing translations.")
        except Exception as e:
            logger.warning(f"Could not load existing output file: {e}")

    # Force re-translation of terms whose English name contains the filter substring
    if args.retranslate_filter:
        removed = {
            tid for tid, t in existing_translations.items()
            if args.retranslate_filter.lower() in t.get('english_technical_name', '').lower()
        }
        for tid in removed:
            del existing_translations[tid]
        logger.info(f"--retranslate-filter '{args.retranslate_filter}': forcing re-translation of {len(removed)} terms.")

    logger.info(f"Loading ontology from {args.input}...")
    try:
        onto = pronto.Ontology(args.input)
    except Exception as e:
        logger.error(f"Failed to load ontology: {e}")
        return
    
    terms_to_translate = []
    all_terms_map = {}
    
    # Load target IDs if provided
    target_ids = None
    if args.target_ids and os.path.exists(args.target_ids):
        with open(args.target_ids, 'r') as f:
            target_ids = {line.strip() for line in f if line.strip().startswith('HP:')}
        logger.info(f"Filtering for {len(target_ids)} target IDs from {args.target_ids}")

    for term in onto.terms():
        if term.obsolete:
            continue
            
        # If target_ids provided, only process those
        if target_ids is not None and term.id not in target_ids:
            continue

        all_synonyms = [s.description for s in term.synonyms]
        layperson_synonyms = [s.description for s in term.synonyms if s.type and s.type.id == 'layperson']
        parent_names = [p.name for p in term.superclasses(distance=1) if p.name and p.id != term.id]
        
        # Clean definition text (remove metadata brackets like [PMID:123])
        raw_def = str(term.definition) if term.definition else ""
        # Remove trailing metadata like " ... [PMID:123, ...]"
        if ' [' in raw_def:
            # We want to keep everything before the last metadata block that starts with ' ['
            # Usually, the metadata is at the end.
            parts = raw_def.split(' [')
            # If there's multiple blocks, we try to take all but the last one if it looks like metadata
            if len(parts) > 1 and (':' in parts[-1] or 'http' in parts[-1]):
                 definition = ' ['.join(parts[:-1]).strip().strip('"').strip("'")
            else:
                 definition = raw_def.strip().strip('"').strip("'")
        else:
            definition = raw_def.strip().strip('"').strip("'")
        
        term_payload = {
            "id": term.id,
            "english_technical_name": term.name,
            "all_english_synonyms": all_synonyms,
            "english_layperson_synonyms": layperson_synonyms,
            "english_definition": definition,
            "context_parent_terms": parent_names
        }
        all_terms_map[term.id] = term_payload
        
        # If target_ids provided, we FORCE re-translation of those terms
        if target_ids is not None:
            if term.id in target_ids:
                terms_to_translate.append(term_payload)
        elif term.id not in existing_translations:
            terms_to_translate.append(term_payload)

    logger.info(f"Terms remaining: {len(terms_to_translate)}")
    
    if len(terms_to_translate) == 0:
        logger.info("All terms already translated.")
        return

    if args.limit:
        terms_to_translate = terms_to_translate[:args.limit]

    batches = list(chunk_list(terms_to_translate, args.batch_size))
    logger.info(f"Processing {len(batches)} batches (size {args.batch_size}).")

    semaphore = asyncio.Semaphore(args.concurrency)
    async with aiohttp.ClientSession() as session:
        tasks = [
            translate_batch(session, batch, api_key, args.model, semaphore)
            for batch in batches
        ]
        
        try:
            from tqdm.asyncio import tqdm
            pbar = tqdm(total=len(batches), desc="Translating batches")
            for f in asyncio.as_completed(tasks):
                batch_res = await f
                if batch_res:
                    for item in batch_res:
                        tid = item['id']
                        if tid in all_terms_map:
                            existing_translations[tid] = {
                                "id": tid,
                                "english_technical_name": all_terms_map[tid]['english_technical_name'],
                                "arabic_technical_name": item['arabic_technical_name'],
                                "arabic_layperson_synonym": item['arabic_layperson_synonym'],
                                "arabic_definition": item['arabic_definition']
                            }
                    with open(args.output, 'w', encoding='utf-8') as f_out:
                        json.dump(list(existing_translations.values()), f_out, ensure_ascii=False, indent=2)
                pbar.set_postfix({"cost": f"${usage_tracker.total_cost:.4f}"})
                pbar.update(1)
            pbar.close()
        except ImportError:
            for batch_res in await asyncio.gather(*tasks):
                if batch_res:
                    for item in batch_res:
                        tid = item['id']
                        if tid in all_terms_map:
                            existing_translations[tid] = {
                                "id": tid,
                                "english_technical_name": all_terms_map[tid]['english_technical_name'],
                                "arabic_technical_name": item['arabic_technical_name'],
                                "arabic_layperson_synonym": item['arabic_layperson_synonym'],
                                "arabic_definition": item['arabic_definition']
                            }
            with open(args.output, 'w', encoding='utf-8') as f_out:
                json.dump(list(existing_translations.values()), f_out, ensure_ascii=False, indent=2)

    logger.info(f"Complete. Total: {len(existing_translations)}")
    logger.info(f"API Usage: {usage_tracker}")

def main():
    parser = argparse.ArgumentParser(description="Batch translate HPO terms with improved definition handling.")
    parser.add_argument("--input", default="ontology/hp.obo", help="Path to hp.obo")
    parser.add_argument("--output", default="translation/hpo_arabic_translations.json", help="Path to output JSON file")
    parser.add_argument("--model", default="openai/gpt-4o", help="OpenRouter model string")
    parser.add_argument("--api-key", help="OpenRouter API key")
    parser.add_argument("--concurrency", type=int, default=5, help="Number of concurrent API requests")
    parser.add_argument("--batch-size", type=int, default=10, help="Number of terms per batch (lowered for long definitions)")
    parser.add_argument("--limit", type=int, help="Limit terms for this specific run")
    parser.add_argument("--price-prompt", type=float, default=0.15, help="Price per 1M prompt tokens")
    parser.add_argument("--price-completion", type=float, default=0.60, help="Price per 1M completion tokens")
    parser.add_argument("--target-ids", help="Optional: File containing HPO IDs to focus on (one per line)")
    parser.add_argument("--retranslate-filter", help="Force re-translation of already-translated terms whose English name contains this substring (e.g. 'Abnormal')")
    
    args = parser.parse_args()
    asyncio.run(process_ontology(args))

if __name__ == "__main__":
    main()
