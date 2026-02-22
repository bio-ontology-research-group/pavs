"""
LLM helpers for Phenotype Matcher v2.

Algorithm 3: llm_dis  — disambiguate among candidate concept labels
Algorithm 4: llm_val_batch — validate which of a set of labels are present in text
"""

import json
import os
from typing import Any, Dict, List, Optional, Set

import requests


# ---------------------------------------------------------------------------
# Low-level HTTP call
# ---------------------------------------------------------------------------

def query_llm(
    system: str,
    user: str,
    api_key: Optional[str] = None,
    model: str = "deepseek/deepseek-chat",
) -> Any:
    """
    POST a chat request to OpenRouter and return the parsed content.

    Returns a dict if the response looks like JSON, otherwise a str.
    Raises RuntimeError on HTTP / parsing failure.
    """
    key = api_key or os.environ.get("OPENROUTER_API_KEY", "")
    if not key:
        raise RuntimeError(
            "No API key: set OPENROUTER_API_KEY or pass api_key= to MatcherConfig."
        )

    payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
        "temperature": 0.0,
    }
    headers = {
        "Authorization": f"Bearer {key}",
        "Content-Type": "application/json",
    }

    resp = requests.post(
        "https://openrouter.ai/api/v1/chat/completions",
        json=payload,
        headers=headers,
        timeout=30,
    )
    resp.raise_for_status()
    content = resp.json()["choices"][0]["message"]["content"].strip()

    # Try to parse as JSON
    try:
        # Strip markdown code fences if present
        stripped = content.strip("`").strip()
        if stripped.startswith("json"):
            stripped = stripped[4:].strip()
        return json.loads(stripped)
    except (json.JSONDecodeError, ValueError):
        return content


# ---------------------------------------------------------------------------
# Algorithm 3: Disambiguation
# ---------------------------------------------------------------------------

def llm_dis(
    p_cand: List[str],
    context: str,
    primary_label: Dict[str, str],
    cue: Optional[str] = None,
    api_key: Optional[str] = None,
    model: str = "deepseek/deepseek-chat",
    cache: Optional[dict] = None,
) -> Set[str]:
    """
    Algorithm 3 — LLM disambiguation.

    Given a list of candidate term IDs (*p_cand*), ask the LLM to pick the
    one best matching the clinical *context*.  Returns a singleton set with
    the chosen ID, or empty set on failure.

    *cue* is an optional extra hint (e.g. a gene symbol).
    *cache* is an optional shared dict for memoisation across calls.
    """
    if not p_cand:
        return set()

    cache_key = ("dis", frozenset(p_cand), context.lower().strip(), model)
    if cache is not None and cache_key in cache:
        return cache[cache_key]

    label_to_id: Dict[str, str] = {}
    lines: List[str] = []
    for pid in p_cand:
        lbl = primary_label.get(pid, pid)
        label_to_id[lbl] = pid
        lines.append(f"  - {lbl}  ({pid})")

    candidates_text = "\n".join(lines)
    cue_text = f"\nAdditional hint: {cue}" if cue else ""

    system = (
        "You are a clinical ontology expert. "
        "Select the single best concept label that matches the clinical context. "
        "Return ONLY the exact label string, nothing else."
    )
    user = (
        f"Clinical context: {context}{cue_text}\n\n"
        f"Candidate concepts:\n{candidates_text}\n\n"
        "Return the exact Label of the best match."
    )

    try:
        result = query_llm(system, user, api_key=api_key, model=model)
    except Exception:
        return set()

    ret: Set[str] = set()
    if isinstance(result, str):
        sel = result.strip().strip('"').strip("'")
        if sel in label_to_id:
            ret = {label_to_id[sel]}
        else:
            # Partial match fallback
            for lbl, pid in label_to_id.items():
                if lbl.lower() in sel.lower() or sel.lower() in lbl.lower():
                    ret = {pid}
                    break

    if cache is not None:
        cache[cache_key] = ret
    return ret


# ---------------------------------------------------------------------------
# Algorithm 4: Batched validation
# ---------------------------------------------------------------------------

def llm_val_batch(
    p_list: List[str],
    text: str,
    primary_label: Dict[str, str],
    api_key: Optional[str] = None,
    model: str = "deepseek/deepseek-chat",
    cache: Optional[dict] = None,
) -> Set[str]:
    """
    Algorithm 4 — Batched LLM validation.

    Given a list of candidate term IDs, ask the LLM which labels are actually
    present (or implied) in *text*.  Returns a set of valid IDs.

    The LLM is expected to return JSON: {"valid": ["exact label", ...]}
    *cache* is an optional shared dict for memoisation across calls.
    """
    if not p_list:
        return set()

    cache_key = ("val", frozenset(p_list), text.lower().strip(), model)
    if cache is not None and cache_key in cache:
        return cache[cache_key]

    label_to_id: Dict[str, str] = {}
    bullet_lines: List[str] = []
    for pid in p_list:
        lbl = primary_label.get(pid, pid)
        label_to_id[lbl] = pid
        bullet_lines.append(f"  • {lbl}")

    bullets = "\n".join(bullet_lines)

    system = (
        "You are a clinical NLP expert. "
        "Determine which of the listed medical concepts are present or clearly "
        "implied in the given clinical text, including obvious typos, plurals, "
        "inflected forms, and abbreviations. "
        'Return JSON exactly as: {"valid": ["exact label 1", "exact label 2", ...]}'
    )
    user = (
        f"Clinical text: {text}\n\n"
        f"Candidate concepts:\n{bullets}\n\n"
        'Return JSON {"valid": [labels present or clearly implied in the text]}.'
    )

    try:
        result = query_llm(system, user, api_key=api_key, model=model)
    except Exception:
        return set()

    valid_ids: Set[str] = set()
    if isinstance(result, dict):
        valid_labels = result.get("valid", [])
        if isinstance(valid_labels, list):
            for lbl in valid_labels:
                lbl_str = str(lbl).strip()
                if lbl_str in label_to_id:
                    valid_ids.add(label_to_id[lbl_str])
    elif isinstance(result, str):
        # Try to recover from malformed response
        for lbl, pid in label_to_id.items():
            if lbl.lower() in result.lower():
                valid_ids.add(pid)

    if cache is not None:
        cache[cache_key] = valid_ids
    return valid_ids


# ---------------------------------------------------------------------------
# Fix 6a: LLM-based post-hoc negation and modifier validation
# ---------------------------------------------------------------------------

def llm_check_negation(
    term_label: str,
    context: str,
    api_key: Optional[str] = None,
    model: str = "deepseek/deepseek-chat",
    cache: Optional[dict] = None,
) -> bool:
    """
    Fix 6a — Return True if *term_label* is expressed as absent/excluded in
    *context*.  Used to validate (and potentially flip) rule-based negation
    decisions.
    """
    cache_key = ("llm_neg", term_label.lower(), context.lower().strip(), model)
    if cache is not None and cache_key in cache:
        return cache[cache_key]

    system = (
        "You are a clinical NLP expert. "
        "Answer only with the single word: true or false"
    )
    user = (
        f"In the clinical text: '{context}'\n"
        f"Is the phenotype '{term_label}' expressed as absent, excluded, "
        f"or not present?\n"
        "Reply only with: true or false"
    )

    try:
        raw = query_llm(system, user, api_key=api_key, model=model)
        result = str(raw).strip().lower().startswith("true")
    except Exception:
        result = True  # conservative: keep negation on error

    if cache is not None:
        cache[cache_key] = result
    return result


def llm_check_modifier(
    term_label: str,
    modifier_label: str,
    context: str,
    api_key: Optional[str] = None,
    model: str = "deepseek/deepseek-chat",
    cache: Optional[dict] = None,
) -> bool:
    """
    Fix 6a — Return True if *modifier_label* correctly describes the severity
    or degree of *term_label* in *context*.  Used to validate (and potentially
    drop) rule-based modifier assignments.
    """
    cache_key = (
        "llm_mod",
        term_label.lower(),
        modifier_label.lower(),
        context.lower().strip(),
        model,
    )
    if cache is not None and cache_key in cache:
        return cache[cache_key]

    system = (
        "You are a clinical NLP expert. "
        "Answer only with the single word: true or false"
    )
    user = (
        f"In the clinical text: '{context}'\n"
        f"Is '{modifier_label}' the severity or degree of '{term_label}'?\n"
        "Reply only with: true or false"
    )

    try:
        raw = query_llm(system, user, api_key=api_key, model=model)
        result = str(raw).strip().lower().startswith("true")
    except Exception:
        result = True  # conservative: keep modifier on error

    if cache is not None:
        cache[cache_key] = result
    return result
