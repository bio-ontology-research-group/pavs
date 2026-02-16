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
    prompt = f"""You are a medical NER (Named Entity Recognition) system. Extract individual phenotype terms from this clinical description.

Clinical text: "{text}"

Your task (BE VERY CAREFUL with splitting):
1. Identify EACH distinct phenotype mentioned
2. For EACH phenotype, extract:
   - The core phenotype term (without modifiers)
   - Any severity modifiers (mild, moderate, severe, profound)
   - Negation status (is it excluded/absent/normal?)
   - The original span of text

3. Handle complex anatomy correctly:
   - "short femur and humerus" → TWO phenotypes: "short femur", "short humerus"
   - "absent radius and tibia" → TWO phenotypes: "absent radius", "absent tibia"
   - "feeding difficulties with developmental delay" → TWO phenotypes: "feeding difficulties", "developmental delay"
   - "hypotonia and intellectual disability" → TWO phenotypes: "hypotonia", "intellectual disability"

4. Handle conjunctions correctly:
   - Split on: "and", "with", "or", ","
   - "A, B, and C with D" → FOUR phenotypes: "A", "B", "C", "D"

5. Detect negation keywords: "no", "not", "without", "normal" (but NOT "absent" in anatomical context like "absent radius")

6. Keep acronyms as-is (e.g., "ASD" stays "ASD", not expanded)

7. Ignore non-phenotypes: locations (PICU, ICU), measurements (IQ), procedures

Output JSON array:
[
  {{
    "term": "core phenotype without modifiers",
    "modifiers": ["severity", "temporal", etc.],
    "excluded": true/false,
    "original_span": "exact text from input"
  }},
  ...
]

CRITICAL RULES:
- Extract ALL individual phenotypes separately
- "A and B" → TWO phenotypes, not one
- "feeding difficulties with developmental delay" → TWO phenotypes
- Keep acronyms as-is (ASD, DD, etc.) - DO NOT expand them
- For anatomy, expand patterns: "short femur and humerus" → "short femur" + "short humerus"

EXAMPLES:

Example 1:
Input: "seizures, hypotonia, and feeding difficulties with developmental delay"
Output: [
  {{"term": "seizures", "modifiers": [], "excluded": false}},
  {{"term": "hypotonia", "modifiers": [], "excluded": false}},
  {{"term": "feeding difficulties", "modifiers": [], "excluded": false}},
  {{"term": "developmental delay", "modifiers": [], "excluded": false}}
]

Example 2:
Input: "short femur and numerus with absent radius and tibia"
Output: [
  {{"term": "short femur", "modifiers": [], "excluded": false}},
  {{"term": "short numerus", "modifiers": [], "excluded": false}},
  {{"term": "absent radius", "modifiers": [], "excluded": false}},
  {{"term": "absent tibia", "modifiers": [], "excluded": false}}
]

Example 3:
Input: "no seizures, normal vision, but hypotonia"
Output: [
  {{"term": "seizures", "modifiers": [], "excluded": true}},
  {{"term": "vision abnormality", "modifiers": ["normal"], "excluded": true}},
  {{"term": "hypotonia", "modifiers": [], "excluded": false}}
]

Example 4:
Input: "ASD, heart murmur"
Output: [
  {{"term": "ASD", "modifiers": [], "excluded": false}},
  {{"term": "heart murmur", "modifiers": [], "excluded": false}}
]
(Note: ASD kept as-is, NOT expanded to "Atrial septal defect")
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

        # Handle both array and object with "phenotypes" key
        if isinstance(data, list):
            return data
        elif isinstance(data, dict) and "phenotypes" in data:
            return data["phenotypes"]
        elif isinstance(data, dict) and "terms" in data:
            return data["terms"]
        else:
            logger.warning(f"Unexpected NER response format: {data}")
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
