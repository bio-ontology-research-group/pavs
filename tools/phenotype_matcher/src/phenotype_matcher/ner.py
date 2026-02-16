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
    prompt = f"""You are a medical span detection system. Your ONLY job is to identify phenotype spans in clinical text.

TASK: Detect phenotype spans (exact text boundaries) and their properties.

===== FEW-SHOT EXAMPLES =====

Example 1: Basic conjunction splitting
Input: "seizures, hypotonia, and feeding difficulties with developmental delay"
Output:
[
  {{"span": "seizures", "term": "seizures", "modifiers": [], "excluded": false}},
  {{"span": "hypotonia", "term": "hypotonia", "modifiers": [], "excluded": false}},
  {{"span": "feeding difficulties", "term": "feeding difficulties", "modifiers": [], "excluded": false}},
  {{"span": "developmental delay", "term": "developmental delay", "modifiers": [], "excluded": false}}
]

Example 2: Anatomical expansion ("and" creates multiple instances)
Input: "short femur and numerus with absent radius and tibia"
Output:
[
  {{"span": "short femur", "term": "short femur", "modifiers": [], "excluded": false}},
  {{"span": "numerus", "term": "short numerus", "modifiers": [], "excluded": false}},
  {{"span": "absent radius", "term": "absent radius", "modifiers": [], "excluded": false}},
  {{"span": "tibia", "term": "absent tibia", "modifiers": [], "excluded": false}}
]

Example 3: Negation detection
Input: "no seizures, normal vision, but hypotonia and intellectual disability present"
Output:
[
  {{"span": "no seizures", "term": "seizures", "modifiers": [], "excluded": true}},
  {{"span": "normal vision", "term": "vision", "modifiers": [], "excluded": true}},
  {{"span": "hypotonia", "term": "hypotonia", "modifiers": [], "excluded": false}},
  {{"span": "intellectual disability", "term": "intellectual disability", "modifiers": [], "excluded": false}}
]

Example 4: Severity modifiers
Input: "severe intellectual disability, mild hypotonia, profound hearing loss"
Output:
[
  {{"span": "severe intellectual disability", "term": "intellectual disability", "modifiers": ["severe"], "excluded": false}},
  {{"span": "mild hypotonia", "term": "hypotonia", "modifiers": ["mild"], "excluded": false}},
  {{"span": "profound hearing loss", "term": "hearing loss", "modifiers": ["profound"], "excluded": false}}
]

Example 5: Acronyms (keep as-is)
Input: "ASD, heart murmur, cyanosis"
Output:
[
  {{"span": "ASD", "term": "ASD", "modifiers": [], "excluded": false}},
  {{"span": "heart murmur", "term": "heart murmur", "modifiers": [], "excluded": false}},
  {{"span": "cyanosis", "term": "cyanosis", "modifiers": [], "excluded": false}}
]

Example 6: Complex case with modifiers and negation
Input: "Caught, feeding difficulties, noisy breathing, PICU addition due to ventricular arrhythmia and post cardiac arrest"
Output:
[
  {{"span": "Caught", "term": "Caught", "modifiers": [], "excluded": false}},
  {{"span": "feeding difficulties", "term": "feeding difficulties", "modifiers": [], "excluded": false}},
  {{"span": "noisy breathing", "term": "noisy breathing", "modifiers": [], "excluded": false}},
  {{"span": "ventricular arrhythmia", "term": "ventricular arrhythmia", "modifiers": [], "excluded": false}},
  {{"span": "post cardiac arrest", "term": "cardiac arrest", "modifiers": ["post"], "excluded": false}}
]
(Note: "PICU addition" is ignored as it's a location/procedure, not a phenotype)

===== YOUR TASK =====

Input: "{text}"

RULES:
1. Detect ALL phenotype spans (exact text boundaries from input)
2. For each span, extract:
   - span: exact text from input
   - term: core phenotype (without negation/severity words)
   - modifiers: list of modifiers (severe, mild, profound, post, etc.)
   - excluded: true if negated (no, not, without, normal)
3. Split on conjunctions: "and", "with", "or", ","
4. For anatomy patterns like "short A and B", expand to "short A" + "short B"
5. Keep acronyms unchanged (ASD, DD, VSD, etc.)
6. Ignore locations (PICU, ICU, ER), measurements (IQ), procedures

Output as JSON array (use "result" key):
{{"result": [...]}}
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
