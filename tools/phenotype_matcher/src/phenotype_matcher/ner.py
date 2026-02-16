"""
Named Entity Recognition for phenotype extraction.

Uses LLM to extract individual phenotype terms from complex descriptions
before performing RAG search. This solves the splitting problem.
"""

import json
import re
import logging
from typing import List, Dict, Any, Optional
import requests

logger = logging.getLogger(__name__)


def extract_phenotype_terms_llm(
    text: str,
    api_key: str,
    llm_model: str = "openai/gpt-oss-120b",
) -> List[Dict[str, Any]]:
    """
    Use LLM to extract individual phenotype terms from text.

    This solves the complex splitting problem by using an LLM to identify
    distinct phenotypes, their modifiers, and negation status.

    Args:
        text: Clinical description (may contain multiple phenotypes)
        api_key: OpenRouter API key
        llm_model: LLM model to use

    Returns:
        List of extracted terms with metadata:
        [
            {
                "term": "intellectual disability",
                "modifiers": ["severe"],
                "excluded": false,
                "original_span": "severe intellectual disability"
            },
            {
                "term": "seizures",
                "modifiers": [],
                "excluded": true,
                "original_span": "no seizures"
            },
            ...
        ]
    """
    prompt = f"""MEDICAL PHENOTYPE SPAN DETECTOR

Your task: Split the input text into individual phenotype mentions. Follow the examples EXACTLY.

===== STEP-BY-STEP PROCESS =====
1. Find ALL delimiters: commas (,) and conjunctions (and, with, or)
2. Split text at each delimiter
3. For anatomy patterns "ADJ NOUN1 and NOUN2", create "ADJ NOUN1" + "ADJ NOUN2"
4. Extract modifiers (severe, mild, profound, post, pre)
5. Detect negation (no, not, without, normal)

===== FEW-SHOT EXAMPLES =====

Example 1: Basic splitting (split on , and "and" and "with")
Input: "seizures, hypotonia, and feeding difficulties with developmental delay"
Step 1: Find delimiters: [,] [and] [with]
Step 2: Split: ["seizures", "hypotonia", "feeding difficulties", "developmental delay"]
Output:
[
  {{"term": "seizures", "modifiers": [], "excluded": false}},
  {{"term": "hypotonia", "modifiers": [], "excluded": false}},
  {{"term": "feeding difficulties", "modifiers": [], "excluded": false}},
  {{"term": "developmental delay", "modifiers": [], "excluded": false}}
]

Example 2: Anatomical pattern (CARRY FORWARD adjectives)
Input: "short femur and numerus with absent radius and tibia"
Step 1: Split on "and" and "with": ["short femur", "numerus", "absent radius", "tibia"]
Step 2: Carry forward "short" to "numerus" (pattern: "short NOUN1 and NOUN2")
Step 3: Carry forward "absent" to "tibia" (pattern: "absent NOUN1 and NOUN2")
Output:
[
  {{"term": "short femur", "modifiers": [], "excluded": false}},
  {{"term": "short numerus", "modifiers": [], "excluded": false}},
  {{"term": "absent radius", "modifiers": [], "excluded": false}},
  {{"term": "absent tibia", "modifiers": [], "excluded": false}}
]

Example 3: Negation detection (remove negation words, set excluded=true)
Input: "no seizures, normal vision, but hypotonia and intellectual disability present"
Step 1: Split on [,] [and]: ["no seizures", "normal vision", "but hypotonia", "intellectual disability present"]
Step 2: Detect negation: "no" in "no seizures", "normal" in "normal vision"
Step 3: Remove negation words and "but", "present"
Output:
[
  {{"term": "seizures", "modifiers": [], "excluded": true}},
  {{"term": "vision", "modifiers": [], "excluded": true}},
  {{"term": "hypotonia", "modifiers": [], "excluded": false}},
  {{"term": "intellectual disability", "modifiers": [], "excluded": false}}
]

Example 4: Severity modifiers (extract and remove from term)
Input: "severe intellectual disability, mild hypotonia, profound hearing loss"
Step 1: Split on [,]: ["severe intellectual disability", "mild hypotonia", "profound hearing loss"]
Step 2: Extract modifiers: "severe", "mild", "profound"
Step 3: Remove modifiers from terms
Output:
[
  {{"term": "intellectual disability", "modifiers": ["severe"], "excluded": false}},
  {{"term": "hypotonia", "modifiers": ["mild"], "excluded": false}},
  {{"term": "hearing loss", "modifiers": ["profound"], "excluded": false}}
]

Example 5: Acronyms (keep unchanged)
Input: "ASD, heart murmur, cyanosis"
Step 1: Split on [,]: ["ASD", "heart murmur", "cyanosis"]
Step 2: Keep acronyms as-is (do NOT expand ASD to "Atrial septal defect")
Output:
[
  {{"term": "ASD", "modifiers": [], "excluded": false}},
  {{"term": "heart murmur", "modifiers": [], "excluded": false}},
  {{"term": "cyanosis", "modifiers": [], "excluded": false}}
]

Example 6: Complex case (ignore non-phenotypes, extract modifiers)
Input: "Caught, feeding difficulties, noisy breathing, PICU addition due to ventricular arrhythmia and post cardiac arrest"
Step 1: Split on [,] [and]: ["Caught", "feeding difficulties", "noisy breathing", "PICU addition due to ventricular arrhythmia", "post cardiac arrest"]
Step 2: Remove "PICU addition due to" (location phrase)
Step 3: Extract "post" as modifier
Output:
[
  {{"term": "Caught", "modifiers": [], "excluded": false}},
  {{"term": "feeding difficulties", "modifiers": [], "excluded": false}},
  {{"term": "noisy breathing", "modifiers": [], "excluded": false}},
  {{"term": "ventricular arrhythmia", "modifiers": [], "excluded": false}},
  {{"term": "cardiac arrest", "modifiers": ["post"], "excluded": false}}
]

===== NOW YOUR TURN =====

Input text to process: "{text}"

ALGORITHM (apply these steps in order):

STEP 1: SPLIT text at delimiters
- Delimiters: "," (comma), " and ", " with ", " or ", " but "
- Example: "A, B and C with D" → ["A", "B", "C", "D"]
- **DO THIS MECHANICALLY - split at EVERY delimiter**

STEP 2: EXPAND anatomical patterns (if ADJ precedes "NOUN1 and NOUN2")
- If you see pattern "ADJ NOUN1 and NOUN2" in step 1, the split gave you ["ADJ NOUN1", "NOUN2"]
- Carry forward the ADJ: "NOUN2" becomes "ADJ NOUN2"
- Example: after step 1 you have ["short femur", "numerus"] → change to ["short femur", "short numerus"]
- Example: after step 1 you have ["absent radius", "tibia"] → change to ["absent radius", "absent tibia"]

STEP 3: CLEAN each term
- Remove filler words: "but", "present", "noted", "due to", "addition"
- Remove location phrases: "PICU", "ICU", etc.

STEP 4: EXTRACT modifiers from each term
- Severity modifiers: "severe", "mild", "moderate", "profound" → put in modifiers list, remove from term
- Temporal modifiers: "post", "pre", "chronic", "episodic" → put in modifiers list, remove from term
- Example: "severe intellectual disability" → term="intellectual disability", modifiers=["severe"]

STEP 5: DETECT negation for each term  
- Negation words: "no", "not", "without", "normal" → set excluded=true, remove word from term
- Example: "no seizures" → term="seizures", excluded=true
- Exception: "absent" in anatomy is part of the phenotype, NOT negation

Output JSON with "result" key containing the array:
{{"result": [...]}}

CRITICAL: You MUST split at EVERY "and", "with", "or", and "," in step 1. Do not group them!
"""

    try:
        resp = requests.post(
            "https://openrouter.ai/api/v1/chat/completions",
            headers={
                "Authorization": f"Bearer {api_key}",
                "Content-Type": "application/json",
            },
            json={
                "model": llm_model,
                "messages": [{"role": "user", "content": prompt}],
                "response_format": {"type": "json_object"},
            },
            timeout=60,
        )

        if resp.status_code != 200:
            logger.error(f"LLM NER API Error: {resp.text}")
            return []

        content = resp.json()["choices"][0]["message"]["content"]
        # Clean markdown code blocks
        content = re.sub(r"```json\n?|\n?```", "", content).strip()
        data = json.loads(content)

        # Handle both array and object with various keys
        if isinstance(data, list):
            return data
        elif isinstance(data, dict):
            # Try common key names
            for key in ["result", "phenotypes", "terms", "extractions", "entities"]:
                if key in data and isinstance(data[key], list):
                    logger.info(f"NER returned data in '{key}' field")
                    return data[key]

            logger.warning(f"Unexpected NER response format: {data}")
            return []
        else:
            logger.warning(f"Unexpected NER response type: {type(data)}")
            return []

    except Exception as e:
        logger.error(f"LLM NER failed: {e}")
        return []


def simple_split_fallback(text: str) -> List[str]:
    """
    Fallback splitting if LLM NER fails.

    Uses simple comma splitting as before.
    """
    return [t.strip() for t in text.split(",") if t.strip()]
