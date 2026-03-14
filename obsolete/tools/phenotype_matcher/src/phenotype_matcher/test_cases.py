"""
Difficult test cases for the Phenotype Matcher tool.

These cases test the tool's ability to handle:
- Typos and misspellings
- Complex multi-phenotype descriptions
- Medical terminology variations
- Long, unstructured clinical descriptions
"""

from typing import List, Dict, Any, Optional
from pyhpo import Ontology


# Initialize ontology for similarity checks
try:
    if not Ontology.hpo_count():
        _ = Ontology()
except:
    pass


def are_terms_related(id1: str, id2: str, max_distance: int = 1) -> bool:
    """
    Check if two HPO terms are related (parent/child within max_distance).
    """
    if id1 == id2:
        return True
        
    try:
        t1 = Ontology.get_hpo_object(id1)
        t2 = Ontology.get_hpo_object(id2)
        
        if not t1 or not t2:
            return False
            
        # Check parents of t1
        parents = t1.parents
        if t2 in parents:
            return True
            
        # Check children of t1
        children = t1.children
        if t2 in children:
            return True
            
        # For max_distance > 1, we could do more, but let's start with 1
        return False
    except:
        return False


# Ground truth verified against hp.obo
DIFFICULT_TEST_CASES = [
    {
        "id": "case_1",
        "name": "Complex cardiac case with typo",
        "input": "Caught, feeding difficulties, noisy breathing, PICU addition due to ventricular arrhythmia and post cardiac arrest",
        "description": "Tests typo correction (Caught→Cough), multiple phenotypes, and clinical context parsing",
        "expected_phenotypes": [
            {
                "hpo_id": "HP:0012735",
                "label": "Cough",
                "comment": "Should detect 'Caught' as typo for 'Cough'",
                "optional": True,
            },
            {
                "hpo_id": "HP:0011968",
                "label": "Feeding difficulties",
                "comment": "Direct match",
            },
            {
                "hpo_id": "HP:0010307",
                "label": "Stridor",
                "comment": "'Noisy breathing' is a synonym/layperson term for Stridor",
                "alternatives": [
                    "HP:0030829"
                ],  # Abnormal breath sound is also acceptable
            },
            {
                "hpo_id": "HP:0004308",
                "label": "Ventricular arrhythmia",
                "comment": "Direct match from clinical context",
                "alternatives": [
                    "HP:0006682"
                ],  # Premature ventricular contraction is also acceptable
            },
            {
                "hpo_id": "HP:0001695",
                "label": "Cardiac arrest",
                "comment": "Direct match from 'post cardiac arrest'",
            },
        ],
        "minimum_expected": 3,  # At least 3 of 5 should be detected (typo is hard)
        "notes": "Typo detection is LLM-dependent. 'PICU addition' should be ignored. Abnormal breath sound is acceptable alternative to Stridor.",
    },
    {
        "id": "case_2",
        "name": "Skeletal abnormalities with typo",
        "input": "short femur and numerus with absent radius and tibia, polyhydramnios",
        "description": "Tests typo correction (numerus→humerus), multiple skeletal phenotypes, and prenatal finding",
        "expected_phenotypes": [
            {"hpo_id": "HP:0003097", "label": "Short femur", "comment": "Direct match"},
            {
                "hpo_id": "HP:0005792",
                "label": "Short humerus",
                "comment": "Should detect 'numerus' as typo for 'humerus'",
            },
            {
                "hpo_id": "HP:0003974",
                "label": "Absent radius",
                "comment": "Direct match",
            },
            {
                "hpo_id": "HP:0009556",
                "label": "Absent tibia",
                "comment": "Direct match",
            },
            {
                "hpo_id": "HP:0001561",
                "label": "Polyhydramnios",
                "comment": "Prenatal finding - direct match",
            },
        ],
        "minimum_expected": 5,  # All 5 should be detected
        "notes": "Tests handling of anatomical terms and typo correction in medical terminology.",
    },
    {
        "id": "case_3",
        "name": "Multiple phenotypes with conjunctions",
        "input": "seizures, hypotonia, and feeding difficulties with developmental delay",
        "description": "Tests handling of multiple conjunctions and compound descriptions",
        "expected_phenotypes": [
            {
                "hpo_id": "HP:0001250",
                "label": "Seizure",
                "comment": "Direct match (seizures→Seizure)",
            },
            {"hpo_id": "HP:0001252", "label": "Hypotonia", "comment": "Direct match"},
            {
                "hpo_id": "HP:0011968",
                "label": "Feeding difficulties",
                "comment": "Direct match",
            },
            {
                "hpo_id": "HP:0001263",
                "label": "Global developmental delay",
                "comment": "'developmental delay' should map to this",
            },
        ],
        "minimum_expected": 4,
        "notes": "Tests proper splitting on commas and 'and' conjunctions.",
    },
    {
        "id": "case_4",
        "name": "Negation with multiple phenotypes",
        "input": "no seizures, normal vision, but hypotonia and intellectual disability present",
        "description": "Tests mixed negation and assertion in one description",
        "expected_phenotypes": [
            {
                "hpo_id": "HP:0001250",
                "label": "Seizure",
                "excluded": True,
                "comment": "'no seizures' should be marked as excluded",
            },
            {
                "hpo_id": "HP:0000504",
                "label": "Abnormality of vision",
                "excluded": True,
                "comment": "'normal vision' implies excluded vision abnormality",
            },
            {
                "hpo_id": "HP:0001252",
                "label": "Hypotonia",
                "excluded": False,
                "comment": "Present phenotype",
            },
            {
                "hpo_id": "HP:0001249",
                "label": "Intellectual disability",
                "excluded": False,
                "comment": "Present phenotype",
            },
        ],
        "minimum_expected": 3,
        "notes": "Tests ability to handle mixed negation and assertion. LLM must correctly parse 'but' conjunction.",
    },
    {
        "id": "case_5",
        "name": "Severity with multiple phenotypes",
        "input": "severe intellectual disability, mild hypotonia, profound hearing loss",
        "description": "Tests severity extraction for multiple phenotypes with different severities",
        "expected_phenotypes": [
            {
                "hpo_id": "HP:0001249",
                "label": "Intellectual disability",
                "severity_id": "HP:0012828",
                "severity_label": "Severe",
                "comment": "NER extracts severity as modifier",
                "alternatives": ["HP:0010864"],  # Pre-coordinated term also acceptable
            },
            {
                "hpo_id": "HP:0001252",
                "label": "Hypotonia",
                "severity_id": "HP:0012825",
                "severity_label": "Mild",
                "comment": "NER extracts mild as modifier",
            },
            {
                "hpo_id": "HP:0000365",
                "label": "Hearing impairment",
                "severity_id": "HP:0012829",
                "severity_label": "Profound",
                "comment": "NER extracts profound as modifier, 'hearing loss' → 'Hearing impairment'",
                "alternatives": ["HP:0012715"],  # Pre-coordinated term also acceptable
            },
        ],
        "minimum_expected": 2,
        "notes": "With NER, severity modifiers are extracted separately. Both base term + modifier and pre-coordinated terms are acceptable.",
    },
    {
        "id": "case_6",
        "name": "Ambiguous term - ASD with cardiac context",
        "input": "ASD, heart murmur, cyanosis",
        "description": "Tests disambiguation of 'ASD' (Atrial Septal Defect vs Autism) using cardiac context",
        "context": "heart murmur, cyanosis",
        "expected_phenotypes": [
            {
                "hpo_id": "HP:0001631",
                "label": "Atrial septal defect",
                "comment": "ASD should resolve to cardiac defect (not autism) given cardiac context",
            },
            {
                "hpo_id": "HP:0030148",
                "label": "Heart murmur",
                "comment": "Direct match",
            },
            {
                "hpo_id": "HP:0000961",
                "label": "Cyanosis",
                "comment": "Direct match",
            },
        ],
        "minimum_expected": 2,  # Heart murmur and cyanosis should definitely match
        "notes": "CRITICAL: 'ASD' with cardiac phenotypes should match Atrial Septal Defect (HP:0001631), NOT Autism Spectrum Disorder (HP:0000729). Tests category-based disambiguation.",
        "disambiguation_test": True,
        "correct_hpo_id": "HP:0001631",  # What ASD SHOULD match
        "incorrect_hpo_id": "HP:0000729",  # What ASD should NOT match
    },
    {
        "id": "case_7",
        "name": "Ambiguous term - ASD with neurodevelopmental context",
        "input": "ASD, developmental delay, speech delay",
        "description": "Tests disambiguation of 'ASD' (Autism vs Atrial Septal Defect) using neuro context",
        "context": "developmental delay, speech delay",
        "expected_phenotypes": [
            {
                "hpo_id": "HP:0000729",
                "label": "Autistic behavior",  # or "Autism spectrum disorder"
                "comment": "ASD should resolve to autism (not cardiac) given neurodevelopmental context",
            },
            {
                "hpo_id": "HP:0001263",
                "label": "Global developmental delay",
                "comment": "Developmental delay",
                "optional": True,
            },
            {
                "hpo_id": "HP:0000750",
                "label": "Delayed speech and language development",
                "comment": "Speech delay",
                "optional": True,
            },
        ],
        "minimum_expected": 1,  # At minimum, should get ASD correct
        "notes": "CRITICAL: 'ASD' with neurodevelopmental phenotypes should match Autism (HP:0000729), NOT Atrial Septal Defect (HP:0001631). Tests category-based disambiguation.",
        "disambiguation_test": True,
        "correct_category": "Abnormality of the nervous system",
        "incorrect_category": "Abnormality of the cardiovascular system",
    },
    {
        "id": "case_8",
        "name": "Ambiguous term - ASD with gene hint (NKX2-5, cardiac)",
        "input": "ASD",
        "description": "Tests disambiguation of 'ASD' using gene hint for cardiac gene",
        "gene_hint": "NKX2-5",  # Cardiac transcription factor
        "expected_phenotypes": [
            {
                "hpo_id": "HP:0001631",
                "label": "Atrial septal defect",
                "comment": "NKX2-5 is associated with cardiac defects, should favor ASD=Atrial Septal Defect",
            },
        ],
        "minimum_expected": 1,
        "notes": "CRITICAL: Gene hint 'NKX2-5' (cardiac transcription factor) should favor Atrial Septal Defect. Gene hint should be used for disambiguation ONLY, not override.",
        "disambiguation_test": True,
        "correct_category": "Abnormality of the cardiovascular system",
        "gene_hint_test": True,
    },
    {
        "id": "case_9",
        "name": "Ambiguous term - ASD with gene hint (SHANK3, autism)",
        "input": "ASD",
        "description": "Tests disambiguation of 'ASD' using gene hint for autism gene",
        "gene_hint": "SHANK3",  # Autism-associated gene
        "expected_phenotypes": [
            {
                "hpo_id": "HP:0000729",
                "label": "Autistic behavior",
                "comment": "SHANK3 is associated with autism, should favor ASD=Autism Spectrum Disorder",
            },
        ],
        "minimum_expected": 1,
        "notes": "CRITICAL: Gene hint 'SHANK3' (autism gene) should favor Autism Spectrum Disorder. Gene hint should be soft cue only.",
        "disambiguation_test": True,
        "correct_category": "Abnormality of the nervous system",
        "gene_hint_test": True,
    },
    {
        "id": "case_10",
        "name": "Gene hint should not override good matches",
        "input": "clear autistic behavior and repetitive movements",
        "description": "Tests that gene hint does NOT override clear semantic matches",
        "gene_hint": "NKX2-5",  # Cardiac gene - should NOT override autism match
        "expected_phenotypes": [
            {
                "hpo_id": "HP:0000729",
                "label": "Autistic behavior",
                "comment": "Clear autism description should match autism despite cardiac gene hint",
            },
            {
                "hpo_id": "HP:0000733",
                "label": "Stereotypy",  # or related to repetitive movements
                "comment": "Repetitive movements",
                "optional": True,
            },
        ],
        "minimum_expected": 1,
        "notes": "CRITICAL: Gene hint (NKX2-5, cardiac) should NOT override the clear semantic match to autism. This tests that gene hints are SOFT cues only.",
        "gene_hint_test": True,
        "soft_cue_test": True,
    },
]


def get_test_case(case_id: str) -> Dict[str, Any]:
    """Get a specific test case by ID."""
    for case in DIFFICULT_TEST_CASES:
        if case["id"] == case_id:
            return case
    return None


def get_all_test_cases() -> List[Dict[str, Any]]:
    """Get all test cases."""
    return DIFFICULT_TEST_CASES


def validate_results(case_id: str, output) -> Dict[str, Any]:
    """
    Validate matcher output against expected results for a test case.

    Args:
        case_id: Test case identifier
        output: PhenotypeOutput from matcher

    Returns:
        Dictionary with validation results
    """
    case = get_test_case(case_id)
    if not case:
        return {"error": f"Test case {case_id} not found"}

    results = {
        "case_id": case_id,
        "case_name": case["name"],
        "input": case["input"],
        "total_expected": len(case["expected_phenotypes"]),
        "minimum_expected": case["minimum_expected"],
        "total_matched": len(output.phenotypes),
        "matches": [],
        "missing": [],
        "unexpected": [],
        "pass": False,
    }

    # Check each expected phenotype
    matched_ids = {p.hpo_id for p in output.phenotypes}

    for expected in case["expected_phenotypes"]:
        expected_id = expected["hpo_id"]
        expected_label = expected["label"]
        is_optional = expected.get("optional", False)
        alternatives = expected.get("alternatives", [])
        found = False

        # Check if expected ID or any alternative was matched
        for matched in output.phenotypes:
            # Check for exact or alternative match
            is_exact = matched.hpo_id == expected_id or matched.hpo_id in alternatives
            
            # Check for close match in hierarchy
            is_close = are_terms_related(matched.hpo_id, expected_id, max_distance=1)
            
            if is_exact or is_close:
                found = True
                match_info = {
                    "hpo_id": matched.hpo_id,
                    "expected_id": expected_id,
                    "expected_label": expected_label,
                    "matched_label": matched.label,
                    "comment": expected.get("comment", ""),
                    "correct": True,
                    "close_match": is_close and not is_exact
                }

                if matched.hpo_id in alternatives:
                    match_info["alternative_match"] = True

                # Check exclusion if specified
                if "excluded" in expected:
                    if matched.excluded == expected["excluded"]:
                        match_info["exclusion_correct"] = True
                        match_info["expected_excluded"] = expected["excluded"]
                        match_info["actual_excluded"] = matched.excluded
                    else:
                        match_info["exclusion_correct"] = False
                        match_info["expected_excluded"] = expected["excluded"]
                        match_info["actual_excluded"] = matched.excluded

                # Check severity if specified
                if "severity_id" in expected:
                    # Case 1: Exact severity match
                    if matched.severity_id == expected["severity_id"]:
                        match_info["severity_correct"] = True
                    # Case 2: Pre-coordinated term (ID is in alternatives)
                    elif matched.hpo_id in alternatives:
                        match_info["severity_correct"] = True
                    # Case 3: Matched label contains severity word
                    elif expected.get("severity_label") and expected["severity_label"].lower() in matched.label.lower():
                        match_info["severity_correct"] = True
                    else:
                        match_info["severity_correct"] = False
                        match_info["expected_severity"] = expected.get("severity_label")
                        match_info["actual_severity"] = matched.severity_label

                results["matches"].append(match_info)
                break

        if not found and not is_optional:
            # Only report as missing if not optional
            results["missing"].append(
                {
                    "hpo_id": expected_id,
                    "label": expected_label,
                    "comment": expected.get("comment", ""),
                    "optional": is_optional,
                }
            )

    # Check for unexpected matches
    expected_ids = {e["hpo_id"] for e in case["expected_phenotypes"]}
    # Also include alternatives
    for expected in case["expected_phenotypes"]:
        expected_ids.update(expected.get("alternatives", []))

    for matched in output.phenotypes:
        if matched.hpo_id not in expected_ids:
            results["unexpected"].append(
                {
                    "hpo_id": matched.hpo_id,
                    "label": matched.label,
                    "excluded": matched.excluded,
                }
            )

    # Determine pass/fail
    results["pass"] = len(results["matches"]) >= case["minimum_expected"]
    results["notes"] = case.get("notes", "")

    return results


def print_validation_results(results: Dict[str, Any], output=None):
    """Pretty print validation results."""
    print(f"\n{'=' * 80}")
    print(f"Test Case: {results['case_name']} ({results['case_id']})")
    print(f"{'=' * 80}")
    print(f"Input: {results['input']}")

    # Show NER extracted terms if available
    if output and hasattr(output, "ner_extracted_terms") and output.ner_extracted_terms:
        print(f"\nNER Extracted Terms ({len(output.ner_extracted_terms)}):")
        for i, ner_term in enumerate(output.ner_extracted_terms, 1):
            term_str = ner_term.get("term", "")
            modifiers = ner_term.get("modifiers", [])
            excluded = ner_term.get("excluded", False)
            mod_str = f" [{', '.join(modifiers)}]" if modifiers else ""
            excl_str = " (EXCLUDED)" if excluded else ""
            print(f"  {i}. '{term_str}'{mod_str}{excl_str}")

    print(
        f"\nExpected: {results['total_expected']} phenotypes (minimum: {results['minimum_expected']})"
    )
    print(f"Matched: {results['total_matched']} phenotypes")
    print(f"Correct matches: {len(results['matches'])}")

    if results["matches"]:
        print(f"\n✓ Correct Matches ({len(results['matches'])}):")
        for match in results["matches"]:
            status = "✓"
            details = []

            if "exclusion_correct" in match:
                if match["exclusion_correct"]:
                    # Show exclusion status when it matches expectation
                    if match.get("expected_excluded"):
                        details.append("excluded=True ✓")
                else:
                    status = "⚠"
                    details.append(
                        f"exclusion: expected={match['expected_excluded']}, got={match['actual_excluded']}"
                    )

            if "severity_correct" in match:
                if match["severity_correct"]:
                    # Show severity when it matches
                    details.append(f"severity={match.get('expected_severity')} ✓")
                else:
                    status = "⚠"
                    details.append(
                        f"severity: expected={match['expected_severity']}, got={match['actual_severity']}"
                    )

            detail_str = f" ({'; '.join(details)})" if details else ""
            match_type_str = " (CLOSE MATCH)" if match.get("close_match") else ""
            print(f"  {status} {match['hpo_id']}: {match['matched_label']}{match_type_str}{detail_str}")
            if match.get("close_match"):
                print(f"     └─ Expected: {match['expected_id']}: {match['expected_label']}")
            if match.get("comment"):
                print(f"     └─ {match['comment']}")

    if results["missing"]:
        print(f"\n✗ Missing Expected Phenotypes ({len(results['missing'])}):")
        for missing in results["missing"]:
            print(f"  ✗ {missing['hpo_id']}: {missing['label']}")
            if missing.get("comment"):
                print(f"     └─ {missing['comment']}")

    if results["unexpected"]:
        print(f"\n⚠ Unexpected Matches ({len(results['unexpected'])}):")
        for unexpected in results["unexpected"]:
            excl_str = " (excluded)" if unexpected.get("excluded") else ""
            print(f"  ⚠ {unexpected['hpo_id']}: {unexpected['label']}{excl_str}")

    if results.get("notes"):
        print(f"\nNotes: {results['notes']}")

    status = "PASS ✓" if results["pass"] else "FAIL ✗"
    print(f"\nResult: {status}")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    """Print all test cases for inspection."""
    print("Phenotype Matcher - Difficult Test Cases")
    print("=" * 80)
    print(f"\nTotal test cases: {len(DIFFICULT_TEST_CASES)}\n")

    for case in DIFFICULT_TEST_CASES:
        print(f"Case {case['id']}: {case['name']}")
        print(f"  Input: {case['input']}")
        print(f"  Expected phenotypes: {len(case['expected_phenotypes'])}")
        print(f"  Minimum to pass: {case['minimum_expected']}")
        print(f"  Description: {case['description']}")
        print()
